//! Tape evaluators: batched points, intervals, batched spatial gradients, and
//! full forward-mode derivatives (x, y, z, t, every parameter).

use crate::V3;
use crate::graph::{BinOp, UnOp, Var, bin_f64, dbin, defined, div_f64, dun, un_f64};
use crate::interval::{self, I, add_o, mul_o, sub_o};
use crate::tape::{Op, Tape};

/// Values of the non-spatial inputs.
#[derive(Clone, Debug, PartialEq)]
pub struct Ctx {
    pub t: f64,
    pub params: Vec<f64>,
    /// When set, interval evaluation treats time as this whole range (point
    /// evaluation still uses `t`): rates over a step of a motion
    /// (ink_kernel::motion) bound how far things move in it.
    pub t_span: Option<(f64, f64)>,
    /// Params interval evaluation treats as ranges (tolerances), so a proof
    /// covers every value in each span (point evaluation still uses
    /// `params`, which proofs keep inside the spans).
    pub p_spans: ParamSpans,
}

impl Ctx {
    /// The params at these values, at moment t (no spans).
    pub const fn new(t: f64, params: Vec<f64>) -> Ctx {
        Ctx { t, params, t_span: None, p_spans: ParamSpans::NONE }
    }

    /// The context as evaluating `tapes` reads it: t, the spans, and the
    /// values of the params they read or that are ranged; every other param
    /// NaN, and none past the last kept. A cached answer is computed from
    /// this context and keyed by all of it (`fingerprint`), so it cannot
    /// read what its key leaves out. None: a scope with derived params,
    /// whose derivations are code (`Derive`), not data.
    pub fn project(&self, tapes: &[&Tape]) -> Option<Ctx> {
        if !self.p_spans.derive.0.is_empty() {
            return None;
        }
        let keep: Vec<u32> = tapes.iter().flat_map(|t| t.params_read()).chain(self.p_spans.list.iter().map(|s| s.k)).collect();
        let mut params = vec![f64::NAN; keep.iter().max().map_or(0, |&k| k as usize + 1)];
        for k in keep {
            params[k as usize] = self.params.get(k as usize).copied().unwrap_or(f64::NAN);
        }
        Some(Ctx { params, ..self.clone() })
    }

    /// Hash all of the context a question reads into `fp` (of a projected
    /// one: `project`): t and its span where a tape reads t (an answer that
    /// names a moment gives the asker's: query::clearance), the
    /// params by their numbers in it
    /// (`names`); a ranged param by the model's number too, which a
    /// witness's `params` give. False (nothing hashed): a scope with
    /// derived params.
    pub fn fingerprint(&self, fp: &mut crate::graph::Fingerprint, names: &crate::tape::Names) -> bool {
        // Every field, named: one added must be hashed too.
        let Ctx { t, params, t_span, p_spans } = self;
        // `mask` is the ranged params' numbers' bits, hashed with them.
        let ParamSpans { mask: _, list, derive } = p_spans;
        if !derive.0.is_empty() {
            return false;
        }
        fp.word(names.reads_t as u64);
        if names.reads_t {
            fp.f64(*t).word(t_span.is_some() as u64);
            if let Some((a, b)) = t_span {
                fp.f64(*a).f64(*b);
            }
        }
        fp.word(list.len() as u64);
        for &ParamSpan { k, span, band, derived } in list {
            fp.word(k as u64).word(names.param(k) as u64).f64(span.0).f64(span.1).f64(band.0).f64(band.1).word(derived as u64);
        }
        fp.word(names.params().len() as u64);
        names.params().iter().for_each(|&k| _ = fp.f64(params.get(k as usize).copied().unwrap_or(f64::NAN)));
        true
    }

    /// Param k's range: its span, else its value.
    pub fn param(&self, k: u32) -> I {
        match self.p_spans.get(k) {
            Some((lo, hi)) => I::new(lo, hi),
            None => I::point(self.params[k as usize]),
        }
    }

    /// Whether t and every ranged param lie in their spans (where a witness
    /// must be).
    pub fn in_scope(&self) -> bool {
        self.t_span.is_none_or(|(a, b)| a <= self.t && self.t <= b)
            && self.p_spans.list.iter().all(|s| s.span.0 <= self.params[s.k as usize] && self.params[s.k as usize] <= s.span.1)
    }

    /// The same, with t and the ranged params moved to the middles of their
    /// spans where they lie outside them (a nominal value may lie outside
    /// its tolerance band).
    pub fn into_scope(&self) -> Ctx {
        let mut c = self.clone();
        if let Some(s) = c.t_span
            && !(s.0 <= c.t && c.t <= s.1)
        {
            c.t = span_mid(s);
        }
        for s in &self.p_spans.list {
            let v = &mut c.params[s.k as usize];
            if !(s.span.0 <= *v && *v <= s.span.1) {
                *v = span_mid(s.span);
            }
        }
        c.derive();
        c
    }

    /// One point of the scope: t, and every ranged param narrowed to its
    /// value, still ranged (so an expression reading it keeps its meaning
    /// in real numbers: Tape::fixed). Confirming a witness happens here. A
    /// derived param's span is then its enclosure at that point (it is
    /// confirmed for every value there, so for its own).
    pub fn point(&self) -> Ctx {
        let mut c = Ctx { t_span: None, ..self.clone() };
        // A derived span encloses its param over the others' spans; at a
        // point outside them (a witness must lie inside: in_scope), start
        // from its band.
        let outside = self.p_spans.list.iter().any(|s| !s.derived && !(s.span.0..=s.span.1).contains(&self.params[s.k as usize]));
        for s in c.p_spans.list.iter_mut() {
            if !s.derived {
                let v = self.params[s.k as usize];
                s.span = (v, v);
            } else if outside {
                s.span = s.band;
            }
        }
        c.derive();
        c
    }

    /// Param k's span narrowed to `span` (inside it), its value at the middle.
    pub fn with_param(&self, k: u32, span: (f64, f64)) -> Ctx {
        let mut c = self.clone();
        let mut widened = false;
        if let Some(s) = c.p_spans.list.iter_mut().find(|s| s.k == k) {
            debug_assert!(!s.derived, "a derived param follows the others");
            widened = !(s.span.0 <= span.0 && span.1 <= s.span.1);
            s.span = span;
            c.params[k as usize] = span_mid(span);
        }
        // Derived spans enclose their params over the spans they were
        // derived over, so over a wider one they start from their bands.
        if widened {
            for s in c.p_spans.list.iter_mut().filter(|s| s.derived) {
                s.span = s.band;
            }
        }
        c.derive();
        c
    }

    /// Param k (independent) ranged over `band` alone, its span the whole
    /// band and its value the middle: a piece of an input's range asked as
    /// a question of its own, which a search may halve PARAM_SPLITS times
    /// (ink_kernel::motion::over_input). Derived params are derived again.
    pub fn with_band(&self, k: u32, band: (f64, f64)) -> Ctx {
        let mut c = self.clone();
        let s = ParamSpan { k, span: band, band, derived: false };
        match c.p_spans.list.iter_mut().find(|s| s.k == k) {
            Some(x) => {
                debug_assert!(!x.derived, "a derived param follows the others");
                *x = s;
            }
            None => {
                c.p_spans.list.push(s);
                c.p_spans.list.sort_by_key(|s| s.k);
                c.p_spans.mask |= 1 << (k % 64);
            }
        }
        if c.params.len() <= k as usize {
            c.params.resize(k as usize + 1, f64::NAN);
        }
        c.params[k as usize] = span_mid(band);
        for s in c.p_spans.list.iter_mut().filter(|s| s.derived) {
            s.span = s.band;
        }
        c.derive();
        c
    }

    /// Every ranged param at the middle of its span (still ranged); t as it
    /// is. A point of the scope, so derived params are derived there.
    pub fn param_mids(&self) -> Ctx {
        let mut c = self.clone();
        for s in c.p_spans.list.iter_mut().filter(|s| !s.derived) {
            let m = span_mid(s.span);
            c.params[s.k as usize] = m;
            s.span = (m, m);
        }
        c.derive();
        c
    }

    /// Derive the derived params again over the spans of the others (see
    /// ParamSpans): after a change to those.
    pub fn derive(&mut self) {
        if self.p_spans.derive.0.is_empty() {
            return;
        }
        let ds = self.p_spans.derive.clone();
        for d in &ds.0 {
            d.derive(self);
        }
    }
}

/// Params given as ranges (docs/IR.md, Tolerances): each one's index, its
/// span, and the band the span is part of (searches halve spans, and stop
/// at a fraction of the band). An expression that reads a ranged param
/// means its value in real numbers, even where its span has narrowed to a
/// single value: the f64 fast path for constants and params (`Tape::fixed`)
/// skips it, so a proof and its witness are about the same expression.
///
/// Some ranged params are derived from the others (a solve's unknowns,
/// the root of its equations: `Derive`). Their spans enclose their values
/// over the others' spans; searches never halve them, but narrow the
/// params they are derived from, and every context a search makes of a
/// piece of the scope, or of one point of it (`Ctx::with_param`,
/// `param_mids`, `point`, `into_scope`), derives them again over it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParamSpans {
    /// Bit k % 64 of each ranged param k (see Tape::pmask): a prefilter;
    /// `list` says which are.
    mask: u64,
    /// Sorted by index.
    list: Vec<ParamSpan>,
    /// What derives the derived ones, in order (each may read those
    /// before it).
    derive: Derivations,
}

/// A ranged param: its index, its span now, and its band; whether it is
/// derived from the others (its band then the enclosure over theirs).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParamSpan {
    pub k: u32,
    pub span: (f64, f64),
    pub band: (f64, f64),
    pub derived: bool,
}

/// Ranged params derived from the others over a scope: design's `solve`
/// over ranged params, whose unknowns are the root of its equations,
/// certified to be one and enclosed for every value of the params it reads
/// (docs/IR.md, Solving constraints). What it sets must hold for every
/// value of those params in their spans (outward rounding throughout):
/// proofs cover its spans, and a witness is confirmed over the span it
/// gives at one point of the others.
pub trait Derive: Send + Sync {
    /// The params it derives.
    fn owns(&self) -> &[u32];
    /// The ranged params it reads (independent, or derived before it).
    fn reads(&self) -> &[u32];
    /// Set the spans of its params in `ctx.p_spans` to enclosures of their
    /// values for every value of the params it reads in their spans there,
    /// and their values in `ctx.params` to ones inside those spans (their
    /// values at the params' values, as far as they can be found). Their
    /// spans in `ctx` already enclose those values (the enclosure may start
    /// from them).
    fn derive(&self, ctx: &mut Ctx);
    /// Each of its params' derivatives along each param it reads (`reads`
    /// order), at ctx's values, in floats (0 where unknown): for searches
    /// that push params toward a witness, never for a bound.
    fn rates(&self, ctx: &Ctx) -> Vec<Vec<f64>>;
}

/// A scope's derivations, compared by identity.
#[derive(Clone, Default)]
pub struct Derivations(pub Vec<std::sync::Arc<dyn Derive>>);

impl std::fmt::Debug for Derivations {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Derivations({})", self.0.len())
    }
}

impl PartialEq for Derivations {
    fn eq(&self, other: &Derivations) -> bool {
        self.0.len() == other.0.len() && self.0.iter().zip(&other.0).all(|(a, b)| std::sync::Arc::ptr_eq(a, b))
    }
}

impl ParamSpans {
    pub const NONE: ParamSpans = ParamSpans { mask: 0, list: Vec::new(), derive: Derivations(Vec::new()) };

    /// These params ranged over these bands (lo <= hi, finite).
    pub fn new(bands: impl IntoIterator<Item = (u32, (f64, f64))>) -> ParamSpans {
        let mut list: Vec<ParamSpan> = bands.into_iter().map(|(k, b)| ParamSpan { k, span: b, band: b, derived: false }).collect();
        list.sort_by_key(|s| s.k);
        list.dedup_by_key(|s| s.k);
        ParamSpans { mask: list.iter().fold(0, |m, s| m | 1 << (s.k % 64)), list, derive: Derivations::default() }
    }

    /// With `d`'s params derived, each over `bands` (by index: their
    /// enclosures over the spans of the params it reads), after those
    /// derived already. Its params must not be ranged yet.
    pub fn derived(self, d: std::sync::Arc<dyn Derive>, bands: &[(u32, (f64, f64))]) -> ParamSpans {
        let mut s = self.with_derived(bands);
        s.derive.0.push(d);
        s
    }

    /// With these params ranged over `bands` as derived ones, with nothing
    /// deriving them: what a derivation evaluates its own equations in.
    pub fn with_derived(mut self, bands: &[(u32, (f64, f64))]) -> ParamSpans {
        for &(k, b) in bands {
            assert!(self.get(k).is_none(), "param {k} is ranged already");
            self.list.push(ParamSpan { k, span: b, band: b, derived: true });
            self.mask |= 1 << (k % 64);
        }
        self.list.sort_by_key(|s| s.k);
        self
    }

    /// With param k ranged over `band` too, independent (t ranged over a
    /// motion: ink-lang's motion::ranged); the derivations kept, as none
    /// reads it. It must not be ranged yet.
    pub fn with_band(mut self, k: u32, band: (f64, f64)) -> ParamSpans {
        assert!(self.get(k).is_none(), "param {k} is ranged already");
        self.list.push(ParamSpan { k, span: band, band, derived: false });
        self.mask |= 1 << (k % 64);
        self.list.sort_by_key(|s| s.k);
        self
    }

    /// Param k's span set to `span`, as it is (a derivation's result: no
    /// other span follows it).
    pub fn set_span(&mut self, k: u32, span: (f64, f64)) {
        if let Some(s) = self.list.iter_mut().find(|s| s.k == k) {
            s.span = span;
        }
    }

    pub fn is_empty(&self) -> bool {
        self.list.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &ParamSpan> {
        self.list.iter()
    }

    /// Param k's span, if it is ranged.
    pub fn get(&self, k: u32) -> Option<(f64, f64)> {
        if self.mask & 1 << (k % 64) == 0 {
            return None;
        }
        self.list.iter().find(|s| s.k == k).map(|s| s.span)
    }

    /// Whether param k is ranged and derived from the others.
    pub fn is_derived(&self, k: u32) -> bool {
        self.mask & 1 << (k % 64) != 0 && self.list.iter().any(|s| s.k == k && s.derived)
    }

    /// The independent ranged params that ranged param `k` stands for: k
    /// itself, or for a derived one, those it is derived from
    /// (transitively), sorted.
    pub fn sources(&self, k: u32) -> Vec<u32> {
        let mut out = Vec::new();
        let mut todo = vec![k];
        while let Some(k) = todo.pop() {
            if !self.is_derived(k) {
                out.push(k);
                continue;
            }
            if let Some(d) = self.derive.0.iter().find(|d| d.owns().contains(&k)) {
                todo.extend(d.reads().iter().copied());
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    /// `grad`, derivatives along every param (by index) at ctx's values,
    /// with each derived param's share passed on to the params it is
    /// derived from (the chain rule, by Derive::rates: floats, for witness
    /// searches only).
    pub fn chain(&self, ctx: &Ctx, grad: &mut [f64]) {
        for d in self.derive.0.iter().rev() {
            for (row, &u) in d.rates(ctx).iter().zip(d.owns()) {
                let g = grad[u as usize];
                if g == 0.0 {
                    continue;
                }
                for (&r, &du) in d.reads().iter().zip(row) {
                    grad[r as usize] += g * du;
                }
            }
        }
    }

    /// The independent ranged params' values in `params` (a witness's: the
    /// derived ones follow from them).
    pub fn values(&self, params: &[f64]) -> Vec<(u32, f64)> {
        self.list.iter().filter(|s| !s.derived).map(|s| (s.k, params[s.k as usize])).collect()
    }
}

/// Lane width for batched evaluation (8 AVX-512 vectors of f64).
pub const L: usize = 64;

pub fn eval_points(tape: &Tape, ctx: &Ctx, pts: &[V3], out: &mut Vec<f64>) {
    // A few points (a cell's samples) in narrow lanes, as eval_grads.
    if pts.len() <= 8 { points_in::<8>(tape, ctx, pts, out) } else { points_in::<L>(tape, ctx, pts, out) }
}

/// `eval_points` in chunks of N points.
fn points_in<const N: usize>(tape: &Tape, ctx: &Ctx, pts: &[V3], out: &mut Vec<f64>) {
    out.clear();
    out.reserve(pts.len());
    let mut s = vec![[0.0f64; N]; tape.nslots];
    let (mut x, mut y, mut z) = ([0.0; N], [0.0; N], [0.0; N]);
    for chunk in pts.chunks(N) {
        for k in 0..N {
            let p = chunk.get(k).copied().unwrap_or([0.0; 3]);
            x[k] = p[0];
            y[k] = p[1];
            z[k] = p[2];
        }
        run_lanes(tape, ctx, &x, &y, &z, &mut s);
        out.extend_from_slice(&s[tape.out_slot()][..chunk.len()]);
    }
}

/// `eval_points` split across all cores (for large point sets).
pub fn eval_points_par(tape: &Tape, ctx: &Ctx, pts: &[V3]) -> Vec<f64> {
    use rayon::prelude::*;
    pts.par_chunks(4096)
        .flat_map_iter(|c| {
            let mut out = Vec::with_capacity(c.len());
            eval_points(tape, ctx, c, &mut out);
            out
        })
        .collect()
}

pub fn eval_point(tape: &Tape, ctx: &Ctx, p: V3) -> f64 {
    let mut s = vec![[0.0f64; 1]; tape.nslots];
    run_lanes(tape, ctx, &[p[0]], &[p[1]], &[p[2]], &mut s);
    s[tape.out_slot()][0]
}

/// The tape over N points at once, one lane each (L for batches, 1 for a
/// lone point: the same code, so the same values).
fn run_lanes<const N: usize>(tape: &Tape, ctx: &Ctx, x: &[f64; N], y: &[f64; N], z: &[f64; N], s: &mut [[f64; N]]) {
    for ins in &tape.code {
        let o = ins.out as usize;
        match ins.op {
            Op::Const(v) => s[o] = [v; N],
            Op::Var(Var::X) => s[o] = *x,
            Op::Var(Var::Y) => s[o] = *y,
            Op::Var(Var::Z) => s[o] = *z,
            Op::Var(Var::T) => s[o] = [ctx.t; N],
            Op::Var(Var::Param(p)) => s[o] = [ctx.params[p as usize]; N],
            Op::Un(op) => {
                let a = s[ins.a as usize];
                let r = &mut s[o];
                match op {
                    UnOp::Neg => (0..N).for_each(|k| r[k] = -a[k]),
                    UnOp::Abs => (0..N).for_each(|k| r[k] = a[k].abs()),
                    UnOp::Sqrt => (0..N).for_each(|k| r[k] = a[k].max(0.0).sqrt()),
                    UnOp::Square => (0..N).for_each(|k| r[k] = a[k] * a[k]),
                    UnOp::Floor => (0..N).for_each(|k| r[k] = a[k].floor()),
                    _ => (0..N).for_each(|k| r[k] = un_f64(op, a[k])),
                }
            }
            Op::Bin(op) => {
                let a = s[ins.a as usize];
                let b = s[ins.b as usize];
                let r = &mut s[o];
                match op {
                    BinOp::Add => (0..N).for_each(|k| r[k] = defined(a[k] + b[k])),
                    BinOp::Sub => (0..N).for_each(|k| r[k] = defined(a[k] - b[k])),
                    BinOp::Mul => (0..N).for_each(|k| r[k] = defined(a[k] * b[k])),
                    BinOp::Div => (0..N).for_each(|k| r[k] = div_f64(a[k], b[k])),
                    BinOp::Min => (0..N).for_each(|k| r[k] = a[k].min(b[k])),
                    BinOp::Max => (0..N).for_each(|k| r[k] = a[k].max(b[k])),
                    _ => (0..N).for_each(|k| r[k] = bin_f64(op, a[k], b[k])),
                }
            }
            Op::Select => {
                let c = s[ins.a as usize];
                let a = s[ins.b as usize];
                let b = s[ins.c as usize];
                let r = &mut s[o];
                (0..N).for_each(|k| r[k] = if c[k] < 0.0 { a[k] } else { b[k] });
            }
            Op::Tag(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real => s[o] = s[ins.a as usize],
            Op::Frame => s[o] = [0.0; N],
            Op::Grid(gi, m) => {
                let (x, y, z) = (s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]);
                let g = &tape.grids[gi as usize];
                let r = &mut s[o];
                (0..N).for_each(|k| r[k] = g.component([x[k], y[k], z[k]], m));
            }
        }
    }
}

/// Whether instruction i is a constants-and-params value taken in floats
/// (Tape::fixed): one that reads no ranged param. A ranged one means its
/// value in real numbers, bounded like any other. The params' bits
/// (Tape::pmask) rule most out at once; where they meet the ranged ones',
/// the params it reads (Tape::params_of) decide, so a param sharing only a
/// ranged one's bit (64 apart) leaves it in floats.
#[inline]
fn fixed(tape: &Tape, ctx: &Ctx, i: usize) -> bool {
    tape.fixed[i] && (tape.pmask[i] & ctx.p_spans.mask == 0 || tape.params_of(i).iter().all(|&k| ctx.p_spans.get(k).is_none()))
}

/// A constants-and-params instruction in floats, from its argument values.
/// Params are fixed for a query, so this is the value a def would have: the
/// interval bound of such an expression is this point, not an enclosure that
/// widens through long chains (and through atan2 or floor near a cut).
fn fixed_value(ins: &crate::tape::Ins, arg: impl Fn(u32) -> f64, ctx: &Ctx) -> f64 {
    match ins.op {
        Op::Const(v) => v,
        Op::Var(Var::Param(p)) => ctx.params[p as usize],
        Op::Un(op) => un_f64(op, arg(ins.a)),
        Op::Bin(op) => bin_f64(op, arg(ins.a), arg(ins.b)),
        Op::Select => {
            if arg(ins.a) < 0.0 {
                arg(ins.b)
            } else {
                arg(ins.c)
            }
        }
        _ => arg(ins.a),
    }
}

/// Bound the tape over a box. If `choices` is given, record which branch of each
/// min/max/select always wins (indexed by instruction).
pub fn eval_interval(tape: &Tape, ctx: &Ctx, b: [I; 3], choices: Option<&mut [u8]>) -> I {
    eval_interval_in(tape, ctx, b, choices, &mut Vec::new())
}

/// `eval_interval` with the caller's storage for the slots (reused across
/// calls that evaluate many small boxes or points).
pub fn eval_interval_in(tape: &Tape, ctx: &Ctx, b: [I; 3], mut choices: Option<&mut [u8]>, s: &mut Vec<I>) -> I {
    s.clear();
    s.resize(tape.nslots, I::point(0.0));
    let mut slopes = Slopes::of(tape, tape.nslots);
    let (n, mut i) = (tape.code.len(), 0);
    for &j in &tape.jumps {
        let j = j as usize;
        // One inside the window of one before it was bounded there.
        if j >= i {
            plain_run(tape, ctx, b, i, j, s, &mut slopes, choices.as_deref_mut());
            i = fork(tape, &mut Plain { tape, ctx, b, s, slopes: &mut slopes }, j, n, choices.as_deref_mut());
        }
    }
    plain_run(tape, ctx, b, i, n, s, &mut slopes, choices);
    s[tape.out_slot()]
}

/// `forked` over instructions `from..to`, which hold no jump, with Plain's
/// step written out: the slots an argument the loop keeps in registers,
/// not reloaded through the run's fields.
#[inline(always)]
fn plain_run(tape: &Tape, ctx: &Ctx, b: [I; 3], from: usize, to: usize, s: &mut [I], slopes: &mut Option<Slopes>, mut choices: Option<&mut [u8]>) {
    for (i, ins) in tape.code[from..to].iter().enumerate() {
        let i = from + i;
        let mut v = interval_ins(tape, ctx, i, ins, s, b, choices.as_deref_mut());
        if let Some(sl) = slopes {
            v = sl.tighten(tape, ctx, i, ins, b, v);
        }
        s[ins.out as usize] = v;
    }
}

/// An interval evaluator's state over a region, which a jump forks
/// (`forked`): its slots, and whatever it keeps beside them.
trait Fork {
    type Out: Copy;
    /// Instruction `i`'s arguments as they are before it is stepped (its
    /// output's slot may be one's): what a side of a jump there is taken
    /// from.
    type Args;
    /// Instruction `i` bounded and stored.
    fn step(&mut self, i: usize, choices: Option<&mut [u8]>) -> I;
    fn args(&self, i: usize) -> Self::Args;
    /// The bounds of the first two arguments.
    fn bounds(args: &Self::Args) -> (I, I);
    fn ctx(&self) -> &Ctx;
    /// Instruction `i`'s value taken as side `k` (0 below, 1 above) of its
    /// jump, bounded by `side` where it takes that side (jump_sides).
    fn side(&mut self, i: usize, args: &Self::Args, k: usize, side: I);
    /// What instruction `i`'s slot holds.
    fn out(&self, i: usize) -> Self::Out;
    /// Instruction `i`'s slot set to hold `o`.
    fn put(&mut self, i: usize, o: Self::Out);
    /// Instruction `i`'s slot bounding a value either of two bound.
    fn join(&mut self, i: usize, a: Self::Out, b: Self::Out);
    /// Whether every value is read after the run (an evaluator over the
    /// tape's SSA form), not only the output: then every value a fork sets
    /// is joined (the jump's own, taken as each side, and each its window
    /// computes), not only the one it joins at.
    fn every_value(&self) -> bool {
        false
    }
}

/// Instructions `from..to` of the tape's code by `run`. A jump the region
/// straddles once (`jump_sides`) is a choice between two sides, each
/// continuous over the region, as a min is between its operands: what
/// follows it is bounded once with each side, up to where the two bounds
/// join without loss (Tape::rejoin), and their hull taken there. So a polar
/// repeat's copy is bounded as one copy turned near its seam, not as every
/// angle, and an angle near atan2's cut as near pi or near -pi, not as
/// anything between. A jump whose join lies past `to` (one also read
/// outside the side it is in) is bounded across it, as is one straddled
/// more than once (the cell is halved). A choice of a min, max or select
/// holds where both sides make it.
fn forked<F: Fork>(tape: &Tape, run: &mut F, from: usize, to: usize, mut choices: Option<&mut [u8]>) {
    let mut i = from;
    let first = tape.jumps.partition_point(|&j| (j as usize) < from);
    for &j in &tape.jumps[first..] {
        let j = j as usize;
        if j >= to {
            break;
        }
        if j >= i {
            (i..j).for_each(|k| _ = run.step(k, choices.as_deref_mut()));
            i = fork(tape, run, j, to, choices.as_deref_mut());
        }
    }
    (i..to).for_each(|k| _ = run.step(k, choices.as_deref_mut()));
}

/// `forked` from a jump instruction `i` (Tape::jumps): bounded, and where
/// the region straddles it once (`jump_sides`), its window bounded once
/// with each side and joined; unless its join lies past `to`. Where the run
/// goes on. Out of line, as jumps are rare: the loop over the rest stays as
/// tight as one with none.
#[inline(never)]
fn fork<F: Fork>(tape: &Tape, run: &mut F, i: usize, to: usize, mut choices: Option<&mut [u8]>) -> usize {
    // Before the output's slot (which may be an argument's) is written.
    let args = run.args(i);
    let (y, x) = F::bounds(&args);
    let v = run.step(i, choices.as_deref_mut());
    let Some(sides) = jump_sides(tape, run.ctx(), i, &tape.code[i], y, x, v) else { return i + 1 };
    let r = tape.rejoin(i);
    if r <= i || r >= to {
        return i + 1;
    }
    // What the window's slots held before it: the values it reads from
    // before the jump, which the first side may overwrite.
    let keep: Vec<F::Out> = (i + 1..=r).map(|k| run.out(k)).collect();
    run.side(i, &args, 0, sides[0]);
    forked(tape, run, i + 1, r + 1, choices.as_deref_mut());
    // The value it joins at, and where every value is read, every value
    // the sides set: the jump's own (one side each) and its window's.
    let first = run.out(r);
    let window: Vec<F::Out> = if run.every_value() { (i..r).map(|k| run.out(k)).collect() } else { Vec::new() };
    let picked = choices.as_deref().map(|ch| ch[i + 1..=r].to_vec());
    (i + 1..=r).zip(&keep).for_each(|(k, &o)| run.put(k, o));
    run.side(i, &args, 1, sides[1]);
    forked(tape, run, i + 1, r + 1, choices.as_deref_mut());
    for (k, a) in (i..r).zip(window) {
        let b = run.out(k);
        run.join(k, a, b);
    }
    let second = run.out(r);
    run.join(r, first, second);
    if let (Some(ch), Some(p)) = (choices, picked) {
        ch[i + 1..=r].iter_mut().zip(p).for_each(|(c, p)| *c = if *c == p { p } else { crate::tape::BOTH });
    }
    r + 1
}

/// eval_interval's state.
struct Plain<'a> {
    tape: &'a Tape,
    ctx: &'a Ctx,
    b: [I; 3],
    s: &'a mut Vec<I>,
    slopes: &'a mut Option<Slopes>,
}

impl Fork for Plain<'_> {
    type Out = (I, Option<interval::G>);
    type Args = (I, I);
    #[inline(always)]
    fn step(&mut self, i: usize, choices: Option<&mut [u8]>) -> I {
        let ins = &self.tape.code[i];
        let mut v = interval_ins(self.tape, self.ctx, i, ins, self.s, self.b, choices);
        if let Some(sl) = self.slopes {
            v = sl.tighten(self.tape, self.ctx, i, ins, self.b, v);
        }
        self.s[ins.out as usize] = v;
        v
    }
    fn args(&self, i: usize) -> (I, I) {
        let ins = &self.tape.code[i];
        (self.s[ins.a as usize], self.s[ins.b as usize])
    }
    fn bounds(args: &(I, I)) -> (I, I) {
        *args
    }
    fn ctx(&self) -> &Ctx {
        self.ctx
    }
    fn side(&mut self, i: usize, _: &(I, I), _: usize, side: I) {
        let o = self.tape.code[i].out as usize;
        self.s[o] = side;
        if let Some(sl) = self.slopes {
            sl.0[o].v = interval::meet(sl.0[o].v, side);
        }
    }
    fn out(&self, i: usize) -> Self::Out {
        let o = self.tape.code[i].out as usize;
        (self.s[o], self.slopes.as_ref().map(|sl| sl.0[o]))
    }
    fn put(&mut self, i: usize, (v, g): Self::Out) {
        let o = self.tape.code[i].out as usize;
        self.s[o] = v;
        if let (Some(sl), Some(g)) = (&mut *self.slopes, g) {
            sl.0[o] = g;
        }
    }
    fn join(&mut self, i: usize, a: Self::Out, b: Self::Out) {
        let o = self.tape.code[i].out as usize;
        self.s[o] = interval::hull(a.0, b.0);
        if let (Some(sl), Some(ga), Some(gb)) = (&mut *self.slopes, a.1, b.1) {
            sl.0[o] = ga.hull(&gb);
        }
    }
}

/// The two sides of a jump instruction `i` straddles once over the box
/// (bound `v`; an atan2's arguments' bounds `y`, `x`), each a bound of its value
/// where it takes that side, continuous there: a floor whose value is one
/// of two integers, each; atan2 across its cut away from the origin (x < 0
/// throughout), near -pi where y < 0 and near pi where y >= 0 (atan2 is odd
/// in y). None for any other instruction or box.
#[inline(never)]
fn jump_sides(tape: &Tape, ctx: &Ctx, i: usize, ins: &crate::tape::Ins, y: I, x: I, v: I) -> Option<[I; 2]> {
    if fixed(tape, ctx, i) {
        return None;
    }
    match ins.op {
        Op::Un(UnOp::Floor) if v.hi - v.lo == 1.0 => Some([I::point(v.lo), I::point(v.hi)]),
        Op::Bin(BinOp::Atan2) => {
            if !(x.hi < 0.0 && y.lo < 0.0 && y.hi >= 0.0) {
                return None;
            }
            let atan2 = |y: I| interval::bin(BinOp::Atan2, y, x).0;
            let below = atan2(I::new(0.0, -y.lo));
            Some([I::new(-below.hi, -below.lo), atan2(I::new(0.0, y.hi))])
        }
        _ => None,
    }
}

/// The middle of a span (of time, or of a ranged param), whose bounds
/// scope_run builds on: in the span and finite for any finite span (halved
/// before the sum, which cannot overflow: 0.5 (t0 + t1) is infinite for
/// t0 = 1e308, t1 = 1.1e308).
pub fn span_mid((t0, t1): (f64, f64)) -> f64 {
    I { lo: t0, hi: t1 }.mid()
}

/// eval_interval, and also first order in the inputs of the scope the tape
/// reads (time over its span, the ranged params: Inputs::of): a bound over
/// the box for every value of them that blurs each value only by how fast
/// it really changes with them (see scope_run). Where a param is read twice
/// (a block `(box a ...)` placed at a / 2, whose fixed face reads a twice),
/// plain bounds lose the dependency and blur the value by the param's whole
/// band; here it moves by its rate, exactly 0 for that face. Over a step of
/// a motion, a face sliding along itself moves by its rate, 0, not by how
/// far the part slides.
pub fn eval_interval_scope(tape: &Tape, ctx: &Ctx, b: [I; 3], choices: Option<&mut [u8]>) -> I {
    let inputs = Inputs::of([tape], ctx);
    if inputs.is_empty() {
        return eval_interval(tape, ctx, b, choices);
    }
    scope_run(tape, false, ctx, &inputs, b, choices, false).v[tape.out_slot()]
}

/// Conditions eval_interval_cases splits on at once (at most 3^3 cases).
const MAX_CASE_CONDITIONS: usize = 3;

/// eval_interval_scope, met with the hull of its bounds in each case of
/// the selects it leaves undecided: those whose conditions are one value
/// or its negation take their sides together at every point (c < 0, c =
/// 0 or c > 0), so each case is bounded with the branches it takes alone;
/// so are floors whose argument crosses one integer, one argument's alike.
/// A switch read several times (a min's derivative, `Graph::diff`'s
/// select, in each term of a product's; a trapezoid written as clamped
/// pieces, one ramp's speed ending as the next one's starts) is then not
/// bounded by pairing a branch of one read with the other branch of
/// another. Sound: at any point the tape is the case its conditions'
/// signs pick. Over more than MAX_CASE_CONDITIONS conditions,
/// eval_interval_scope's bound alone.
pub fn eval_interval_cases(tape: &Tape, ctx: &Ctx, b: [I; 3]) -> I {
    use crate::tape::{BOTH, FIRST, SECOND};
    let mut ch = vec![BOTH; tape.code.len()];
    let v = eval_interval_scope(tape, ctx, b, Some(&mut ch));
    // Undecided selects by condition, through tags and negations (x - y
    // is -(y - x), in floats too: a min's derivative reads t1 - t where a
    // clamp of t - t1 reads that): (its instruction, whether it reads the
    // condition negated).
    let ssa = &tape.ssa;
    let condition = |mut c: usize| {
        let mut neg = false;
        loop {
            match ssa[c].op {
                Op::Tag(_) => c = ssa[c].a as usize,
                Op::Un(UnOp::Neg) => (c, neg) = (ssa[c].a as usize, !neg),
                _ => break,
            }
        }
        let (x, y) = (ssa[c].a as usize, ssa[c].b as usize);
        match ssa[c].op {
            Op::Bin(BinOp::Sub) => ((x.min(y), x.max(y)), neg ^ (x > y)),
            _ => ((c, usize::MAX), neg),
        }
    };
    // What the value reads: not a derivative's mark (a drive's selects).
    let mut read = vec![false; ssa.len()];
    read[ssa.len() - 1] = true;
    for i in (0..ssa.len()).rev() {
        if read[i] {
            let (args, k) = if let Op::Deriv(_) = ssa[i].op { ([ssa[i].a, 0, 0], 1) } else { ssa[i].args() };
            args[..k].iter().for_each(|&a| read[a as usize] = true);
        }
    }
    let mut groups: Vec<((usize, usize), Vec<(usize, bool)>)> = Vec::new();
    for (i, ins) in ssa.iter().enumerate() {
        if !read[i] || !matches!(ins.op, Op::Select) || ch[i] != BOTH || fixed(tape, ctx, i) || same_branches(tape, i) {
            continue;
        }
        let (c, neg) = condition(ins.a as usize);
        match groups.iter_mut().find(|g| g.0 == c) {
            Some(g) => g.1.push((i, neg)),
            None => groups.push((c, vec![(i, neg)])),
        }
    }
    // Each group's cases, as what they force: the side each of its selects
    // takes where the condition is negative, zero or positive (a select
    // takes its first branch where what it reads is negative).
    type Case = (Vec<(usize, u8)>, Vec<(usize, f64)>);
    let mut cases: Vec<Vec<Case>> = groups
        .iter()
        .map(|(_, sel)| {
            let mut cs: Vec<Case> = Vec::new();
            for sign in [-1.0, 0.0, 1.0] {
                let c: Vec<(usize, u8)> = sel.iter().map(|&(i, neg)| (i, if (if neg { -sign } else { sign }) < 0.0 { FIRST } else { SECOND })).collect();
                if !cs.iter().any(|x| x.0 == c) {
                    cs.push((c, Vec::new()));
                }
            }
            cs
        })
        .collect();
    // Floors whose argument's bounds cross one integer n are n - 1 below
    // it and n from it, every floor of one argument alike (a wrap by
    // `mod`, a - m floor(a / m), at its seam).
    let mut floors: Vec<(usize, f64, Vec<usize>)> = Vec::new();
    for (i, ins) in ssa.iter().enumerate() {
        if !read[i] || !matches!(ins.op, Op::Un(UnOp::Floor)) || fixed(tape, ctx, i) {
            continue;
        }
        if let Some(f) = floors.iter_mut().find(|f| f.0 == ins.a as usize) {
            f.2.push(i);
            continue;
        }
        let a = eval_interval_scope(&tape.cone(ins.a as usize), ctx, b, None);
        if a.hi.floor() == a.lo.floor() + 1.0 {
            floors.push((ins.a as usize, a.hi.floor(), vec![i]));
        }
    }
    for (_, n, is) in &floors {
        cases.push([n - 1.0, *n].map(|k| (Vec::new(), is.iter().map(|&i| (i, k)).collect())).to_vec());
    }
    if cases.is_empty() || cases.len() > MAX_CASE_CONDITIONS {
        return v;
    }
    let mut hull: Option<I> = None;
    let mut pick = vec![0usize; cases.len()];
    loop {
        let (mut forced, mut values) = (vec![BOTH; ssa.len()], Vec::new());
        for (g, &k) in pick.iter().enumerate() {
            for &(i, side) in &cases[g][k].0 {
                forced[i] = side;
            }
            values.extend_from_slice(&cases[g][k].1);
        }
        let w = eval_interval_scope(&tape.with_consts(&values).simplify(&forced), ctx, b, None);
        hull = Some(hull.map_or(w, |h| interval::hull(h, w)));
        // The next case, as an odometer over the groups.
        let Some(g) = (0..pick.len()).find(|&g| pick[g] + 1 < cases[g].len()) else { break };
        pick[g] += 1;
        pick[..g].iter_mut().for_each(|k| *k = 0);
    }
    hull.map_or(v, |h| interval::meet(v, h))
}

/// eval_interval_scope, and how much each input of the scope the tape
/// reads (time over its span, each ranged param) blurs it: the bound's rate
/// along the input over the box and scope times how far the input reaches
/// from its middle (rounded up; infinite where the rate is unbounded), the
/// input's share of what the bound adds to its value at their centre.
/// Halving an input at least halves its share, so a search halves the one
/// with the largest. Every input read is measured: MAX_INPUTS at a time,
/// time in the first, each pass first order in its own (the rest ranged in
/// its centre, so its bound holds for every value of them), the bound the
/// meet of the passes' and the branch choices the first's. A pass per
/// MAX_INPUTS inputs the tape reads; none for one it does not.
pub fn eval_interval_scope_blur(tape: &Tape, ctx: &Ctx, b: [I; 3], mut choices: Option<&mut [u8]>) -> (I, Vec<(Var, f64)>) {
    let read = Inputs::read([tape], ctx);
    if read.is_empty() {
        return (eval_interval(tape, ctx, b, choices), Vec::new());
    }
    let (mut v, mut blur) = (I::ENTIRE, Vec::with_capacity(read.len()));
    for inputs in read.passes() {
        let sc = scope_run(tape, false, ctx, &inputs, b, choices.take(), true);
        let (o, n, hw) = (tape.out_slot(), inputs.len(), inputs.half_widths(ctx));
        v = interval::meet(v, sc.v[o]);
        let axes = inputs.t.then_some(Var::T).into_iter().chain(inputs.params.iter().map(|&k| Var::Param(k)));
        blur.extend(axes.enumerate().map(|(j, x)| {
            let r = interval::mag(sc.d[o * n + j]);
            (x, if r == 0.0 { 0.0 } else { interval::mul_up(r, hw[j]) })
        }));
    }
    (v, blur)
}

/// The tape's bound over the box and the scope, and its rates along
/// `inputs` (time, ranged params; at most MAX_INPUTS), each bounded over the
/// box and every value of the scope by the derivative rules (interval::d_un
/// and the rest: ENTIRE where the value may jump). A solve's certified
/// roots are built on these (the interval Jacobian of its equations).
pub fn eval_interval_rates(tape: &Tape, ctx: &Ctx, b: [I; 3], inputs: &Inputs) -> (I, Vec<I>) {
    assert!(inputs.len() <= MAX_INPUTS, "at most MAX_INPUTS rates at once");
    if inputs.is_empty() {
        return (eval_interval(tape, ctx, b, None), Vec::new());
    }
    let sc = scope_run(tape, false, ctx, inputs, b, None, true);
    let (o, n) = (tape.out_slot(), inputs.len());
    (sc.v[o], sc.d[o * n..(o + 1) * n].to_vec())
}

/// Inputs beyond space a first-order bound takes at once (time and ranged
/// params). Any more ranged params stay ranged in its centre (see Inputs).
pub const MAX_INPUTS: usize = 8;

/// The inputs beyond space a first-order bound is taken in (scope_run): its
/// axes, time first. A ranged param that is not one (beyond MAX_INPUTS)
/// keeps its span in the bound's centre, so the bound still holds for every
/// value of it: the centre is then a box along it, and the bound is first
/// order only in the axes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Inputs {
    /// Time (ctx.t_span).
    pub t: bool,
    /// Ranged params, by index.
    pub params: Vec<u32>,
}

impl Inputs {
    /// Those of ctx's scope that `tapes` read, at most MAX_INPUTS: time,
    /// where ctx spans it (a step of a motion), and each ranged param.
    pub fn of<'a>(tapes: impl IntoIterator<Item = &'a Tape>, ctx: &Ctx) -> Inputs {
        let mut all = Inputs::read(tapes, ctx);
        all.params.truncate(MAX_INPUTS - all.t as usize);
        all
    }

    /// Every input of ctx's scope that `tapes` read: time where ctx spans
    /// more than one moment, and the ranged params, in order.
    fn read<'a>(tapes: impl IntoIterator<Item = &'a Tape>, ctx: &Ctx) -> Inputs {
        let spans_t = ctx.t_span.is_some_and(|(a, b)| a < b);
        let (mut t, mut params) = (false, Vec::new());
        if spans_t || !ctx.p_spans.is_empty() {
            for tape in tapes {
                for ins in &tape.ssa {
                    match ins.op {
                        Op::Var(Var::T) => t |= spans_t,
                        Op::Var(Var::Param(k)) if ctx.p_spans.get(k).is_some() => params.push(k),
                        _ => {}
                    }
                }
            }
            params.sort_unstable();
            params.dedup();
        }
        Inputs { t, params }
    }

    /// These in passes of at most MAX_INPUTS, time in the first.
    fn passes(&self) -> Vec<Inputs> {
        let first = MAX_INPUTS - self.t as usize;
        let mut out = vec![Inputs { t: self.t, params: self.params.iter().copied().take(first).collect() }];
        out.extend(self.params.get(first..).unwrap_or(&[]).chunks(MAX_INPUTS).map(|c| Inputs { t: false, params: c.to_vec() }));
        out.retain(|x| !x.is_empty());
        out
    }

    pub fn len(&self) -> usize {
        self.t as usize + self.params.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The axis of this input: t, or param k.
    fn axis(&self, v: Var) -> Option<usize> {
        match v {
            Var::T if self.t => Some(0),
            Var::Param(k) => self.params.iter().position(|&p| p == k).map(|j| j + self.t as usize),
            _ => None,
        }
    }

    /// ctx at their centre: t at the middle of its span (no span), each
    /// param at the middle of its span, still ranged (so it keeps its
    /// meaning in real numbers: Ctx::point).
    fn centre(&self, ctx: &Ctx) -> Ctx {
        let mut c = ctx.clone();
        if self.t
            && let Some(s) = c.t_span.take()
        {
            c.t = span_mid(s);
        }
        for s in c.p_spans.list.iter_mut().filter(|s| self.params.contains(&s.k)) {
            let m = span_mid(s.span);
            s.span = (m, m);
            c.params[s.k as usize] = m;
        }
        c
    }

    /// The distance from each one's centre to the far end of its span,
    /// rounded up (the centre may round to one end of a narrow span); 0
    /// past the last.
    pub fn half_widths(&self, ctx: &Ctx) -> [f64; MAX_INPUTS] {
        let hw = |s: (f64, f64)| {
            let m = span_mid(s);
            (s.1 - m).max(m - s.0).next_up()
        };
        let t = self.t.then(|| hw(ctx.t_span.expect("a time span")));
        let mut out = [0.0; MAX_INPUTS];
        for (o, h) in out.iter_mut().zip(t.into_iter().chain(self.params.iter().map(|&k| hw(ctx.p_spans.get(k).expect("a ranged param"))))) {
            *o = h;
        }
        out
    }
}

/// Instruction `i` (`ins`, in the tape's slot code or its SSA form)'s bound
/// from its arguments' bounds in `s`.
fn interval_ins(tape: &Tape, ctx: &Ctx, i: usize, ins: &crate::tape::Ins, s: &[I], b: [I; 3], choices: Option<&mut [u8]>) -> I {
    if fixed(tape, ctx, i) {
        return I::point(fixed_value(ins, |k| s[k as usize].lo, ctx));
    }
    match ins.op {
        Op::Const(v) => I::point(v),
        Op::Var(Var::X) => b[0],
        Op::Var(Var::Y) => b[1],
        Op::Var(Var::Z) => b[2],
        Op::Var(Var::T) => ctx.t_span.map_or(I::point(ctx.t), |(a, b)| I::new(a, b)),
        Op::Var(Var::Param(p)) => ctx.param(p),
        Op::Un(op) => interval::un(op, s[ins.a as usize]),
        // A value times itself is its square (never negative).
        Op::Bin(BinOp::Mul) if ins.a == ins.b => interval::un(UnOp::Square, s[ins.a as usize]),
        Op::Bin(op) => {
            let (v, c) = interval::bin(op, s[ins.a as usize], s[ins.b as usize]);
            if let Some(ch) = choices {
                ch[i] = c;
            }
            v
        }
        Op::Select => {
            let (v, c) = interval::select(s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]);
            if let Some(ch) = choices {
                ch[i] = c;
            }
            v
        }
        Op::Tag(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real => s[ins.a as usize],
        Op::Frame => I::point(0.0),
        Op::Grid(gi, m) => {
            let (g, b) = (&tape.grids[gi as usize], [s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]]);
            g.component_bound(b, m)
        }
    }
}

/// A value linear in x, y, z: coefficients and constant term, as intervals
/// that hold the real-number ones (a coefficient may be the range of
/// another value, as in a rotation by an angle that varies).
#[derive(Clone, Copy)]
struct Linear {
    a: [I; 3],
    d: I,
}

/// The region of `eval_interval_linear`: a box, or a far sector: the points
/// c + r (s q) (componentwise) for s in `s` and q in `q`, which lie on the
/// faces of the cube [-1, 1]^3. A sector reaches infinity when s does.
#[derive(Clone, Copy, Debug)]
pub enum Region {
    Box([I; 3]),
    Sector { c: V3, r: V3, q: [I; 3], s: I },
}

impl Region {
    /// The coordinates' ranges.
    pub fn bounds(&self) -> [I; 3] {
        match *self {
            Region::Box(b) => b,
            Region::Sector { c, r, q, s } => std::array::from_fn(|k| add_o(I::point(c[k]), mul_o(I::point(r[k]), mul_o(s, q[k])))),
        }
    }
}

impl Linear {
    fn constant(d: I) -> Linear {
        Linear { a: [I::point(0.0); 3], d }
    }
    fn is_constant(&self) -> bool {
        self.a.iter().all(|c| c.lo == 0.0 && c.hi == 0.0)
    }
    fn map(&self, f: impl Fn(I) -> I) -> Linear {
        Linear { a: self.a.map(&f), d: f(self.d) }
    }
    fn zip(&self, o: &Linear, f: impl Fn(I, I) -> I) -> Linear {
        Linear { a: std::array::from_fn(|k| f(self.a[k], o.a[k])), d: f(self.d, o.d) }
    }
    /// Its range over the region. Over a sector, d + a.c + s (a.(r q)): the
    /// distance s multiplies the direction's sum once, so a rotated
    /// coordinate is one-signed on a narrow sector however far it reaches.
    fn range(&self, region: &Region) -> I {
        match *region {
            Region::Box(b) => (0..3).fold(self.d, |acc, k| add_o(acc, mul_o(self.a[k], b[k]))),
            Region::Sector { c, r, q, s } => {
                let near = (0..3).fold(self.d, |acc, k| add_o(acc, mul_o(self.a[k], I::point(c[k]))));
                let dir = (0..3).fold(I::point(0.0), |acc, k| add_o(acc, mul_o(mul_o(self.a[k], I::point(r[k])), q[k])));
                add_o(near, mul_o(s, dir))
            }
        }
    }
}

/// `eval_interval` over a region, where values linear in x, y, z (sums and
/// multiples of the coordinates: rotations of rotations, projections) are
/// also bounded from their collected coefficients, and each value takes the
/// tighter bound. Plain interval arithmetic loses the dependency between
/// such terms, which over an unbounded region loses everything (the far side
/// of a rotated part); collected, and over a far sector, a rotation is exact
/// up to the sector's width. Sound: both bounds hold for the real-number
/// value.
pub fn eval_interval_linear(tape: &Tape, ctx: &Ctx, region: &Region, choices: Option<&mut [u8]>) -> I {
    let mut run = Lin { tape, ctx, region, b: region.bounds(), s: vec![I::point(0.0); tape.nslots], lin: vec![None; tape.nslots], slopes: Slopes::of(tape, tape.nslots) };
    forked(tape, &mut run, 0, tape.code.len(), choices);
    run.s[tape.out_slot()]
}

/// eval_interval_linear's state.
struct Lin<'a> {
    tape: &'a Tape,
    ctx: &'a Ctx,
    region: &'a Region,
    b: [I; 3],
    s: Vec<I>,
    lin: Vec<Option<Linear>>,
    slopes: Option<Slopes>,
}

impl Fork for Lin<'_> {
    type Out = (I, Option<Linear>, Option<interval::G>);
    type Args = (I, I);
    fn args(&self, i: usize) -> (I, I) {
        let ins = &self.tape.code[i];
        (self.s[ins.a as usize], self.s[ins.b as usize])
    }
    fn bounds(args: &(I, I)) -> (I, I) {
        *args
    }
    fn ctx(&self) -> &Ctx {
        self.ctx
    }
    #[inline(always)]
    fn step(&mut self, i: usize, choices: Option<&mut [u8]>) -> I {
        let (tape, ctx, region, b) = (self.tape, self.ctx, self.region, self.b);
        let (s, lin) = (&self.s, &self.lin);
        let ins = &tape.code[i];
        let mut v = interval_ins(tape, ctx, i, ins, s, b, choices);
        if let Some(sl) = &mut self.slopes {
            v = sl.tighten(tape, ctx, i, ins, b, v);
        }
        let arg = |k: u32| lin[k as usize];
        // A value's linear form, or its range as a constant term (to scale
        // another's coefficients by).
        let scale = |k: u32| lin[k as usize].filter(|l| l.is_constant()).map_or(s[k as usize], |l| l.d);
        // Time is a constant term too: it does not vary over the region.
        let l = if tape.fixed[i] || matches!(ins.op, Op::Const(_) | Op::Var(Var::Param(_) | Var::T)) {
            Some(Linear::constant(v))
        } else {
            match ins.op {
                Op::Var(Var::X) | Op::Var(Var::Y) | Op::Var(Var::Z) => {
                    let k = match ins.op {
                        Op::Var(Var::X) => 0,
                        Op::Var(Var::Y) => 1,
                        _ => 2,
                    };
                    let mut l = Linear::constant(I::point(0.0));
                    l.a[k] = I::point(1.0);
                    Some(l)
                }
                Op::Un(UnOp::Neg) => arg(ins.a).map(|l| l.map(|c| I::new(-c.hi, -c.lo))),
                Op::Tag(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real => arg(ins.a),
                Op::Bin(BinOp::Add) => arg(ins.a).zip(arg(ins.b)).map(|(p, q)| p.zip(&q, add_o)),
                Op::Bin(BinOp::Sub) => arg(ins.a).zip(arg(ins.b)).map(|(p, q)| p.zip(&q, sub_o)),
                // A linear value times another value: the other's range
                // scales its coefficients (exact for a constant).
                Op::Bin(BinOp::Mul) => match (arg(ins.a), arg(ins.b)) {
                    (Some(p), q) if q.is_none_or(|q| q.is_constant()) => Some(p.map(|c| mul_o(c, scale(ins.b)))),
                    (p, Some(q)) if p.is_none_or(|p| p.is_constant()) => Some(q.map(|c| mul_o(c, scale(ins.a)))),
                    _ => None,
                },
                Op::Bin(BinOp::Div) => match arg(ins.a) {
                    Some(p) if arg(ins.b).is_none_or(|q| q.is_constant()) && (scale(ins.b).lo > 0.0 || scale(ins.b).hi < 0.0) => {
                        Some(p.map(|c| interval::bin(BinOp::Div, c, scale(ins.b)).0))
                    }
                    _ => None,
                },
                // One-signed over the region: the form itself, or negated.
                Op::Un(UnOp::Abs) => arg(ins.a).and_then(|l| match s[ins.a as usize] {
                    x if x.lo >= 0.0 => Some(l),
                    x if x.hi <= 0.0 => Some(l.map(|c| I::new(-c.hi, -c.lo))),
                    _ => None,
                }),
                // A member that wins everywhere over the region: its form.
                Op::Bin(BinOp::Min | BinOp::Max) => {
                    let (p, q) = (s[ins.a as usize], s[ins.b as usize]);
                    let first = match ins.op {
                        Op::Bin(BinOp::Min) => (p.hi <= q.lo, q.hi <= p.lo),
                        _ => (p.lo >= q.hi, q.lo >= p.hi),
                    };
                    match first {
                        (true, _) => arg(ins.a),
                        (_, true) => arg(ins.b),
                        _ => None,
                    }
                }
                _ => None,
            }
        };
        // A value that does not depend on x, y, z (of time, say) is a
        // constant term.
        let l = l.or_else(|| {
            let (args, k) = ins.args();
            (k > 0 && args[..k].iter().all(|&a| lin[a as usize].is_some_and(|l| l.is_constant()))).then(|| Linear::constant(v))
        });
        if let Some(l) = l {
            v = interval::meet(v, l.range(region));
        }
        self.s[ins.out as usize] = v;
        self.lin[ins.out as usize] = l;
        v
    }
    fn side(&mut self, i: usize, _: &(I, I), _: usize, side: I) {
        let o = self.tape.code[i].out as usize;
        self.s[o] = side;
        // A floor's side is one number over the region.
        self.lin[o] = (side.lo == side.hi).then(|| Linear::constant(side));
        if let Some(sl) = &mut self.slopes {
            sl.0[o].v = interval::meet(sl.0[o].v, side);
        }
    }
    fn out(&self, i: usize) -> Self::Out {
        let o = self.tape.code[i].out as usize;
        (self.s[o], self.lin[o], self.slopes.as_ref().map(|sl| sl.0[o]))
    }
    fn put(&mut self, i: usize, (v, l, g): Self::Out) {
        let o = self.tape.code[i].out as usize;
        self.s[o] = v;
        self.lin[o] = l;
        if let (Some(sl), Some(g)) = (&mut self.slopes, g) {
            sl.0[o] = g;
        }
    }
    fn join(&mut self, i: usize, a: Self::Out, b: Self::Out) {
        let o = self.tape.code[i].out as usize;
        self.s[o] = interval::hull(a.0, b.0);
        // Coefficients that hold either form's (a value is one or the other).
        self.lin[o] = a.1.zip(b.1).map(|(p, q)| p.zip(&q, interval::hull));
        if let (Some(sl), Some(ga), Some(gb)) = (&mut self.slopes, a.2, b.2) {
            sl.0[o] = ga.hull(&gb);
        }
    }
}

/// eval_interval over a scope beyond space, first order in its `inputs`
/// (time, ranged params), run over the tape's slot code, or with `every`
/// its SSA form (every value kept): every slot's bound over the box for
/// every value of the inputs, its bound at their centre, and, with `rates`,
/// its rates along them.
///
/// Plain interval bounds over a time span blur a moving shape along its
/// motion: a hole turning about its own axis bounds as an annulus, its
/// rotated coordinates as a box around the arc they sweep. Over a ranged
/// param they lose the dependency between two reads of it: a block `(box a
/// ...)` placed at a / 2 has its fixed face blurred by a's band. So each
/// value's bound here is also the mean value form in the inputs u: its
/// bound at their centre c (Inputs::centre) plus, for each input j, its
/// rate d/du_j over the box and the whole scope (a convex box, so it holds
/// on the segment from c to any point of it) times (u_j - c_j). Each value
/// keeps the intersection of that and its plain bound, and passes it on
/// (both contain every value, so their intersection does). The rates come
/// from the derivative rules shared with the gradient bounds
/// (interval::d_un, d_bin, ...), which are unbounded wherever the value may
/// jump or its slope is infinite (floor across an integer, a division by
/// something that may be 0, atan2 across its cut, a select between
/// different fields, ln and sqrt at 0, a sampled field's slope across a
/// cell face or about its source), as the mean value theorem needs; there
/// the plain bound stands, and no bound at the centre is computed. A
/// Lipschitz promise bounds slopes in space, not rates in the inputs: its
/// value's rates are its formula's. Values that do not depend on the inputs
/// are bounded once. A rotation of the hole about its axis then blurs it by
/// the spread of its rate over the box, which shrinks with the box: d/dt of
/// its distance is 0 at every point, and its bound over the box is about
/// the box's width times the turning speed (not the radius times it).
///
/// A jump the box and scope straddle once is a choice of two sides
/// (`forked`), as in eval_interval, each side the jump's branch continued
/// over the whole region (a floor's one integer; atan2's branch near pi or
/// near -pi, interval::atan2_branch), so continuous there, with bounded
/// rates (atan2's formula; a floor's 0). The tape with the jump taken as
/// one side is the tape wherever the jump takes that side, so each side is
/// bounded first order over the whole region, and the two joined where the
/// tape joins them (Tape::rejoin). The joined value's rates are the hull
/// of the sides' where no jump it went through moves with the inputs (each
/// point then takes one side for every value of them); elsewhere it may
/// jump along them, and its rates are unbounded. A turning gear's teeth
/// over a step are bounded as one tooth near its seam, turning.
///
/// `choices`, as in eval_interval, records the branches decided over the
/// whole box and scope.
fn scope_run(tape: &Tape, every: bool, ctx: &Ctx, inputs: &Inputs, b: [I; 3], choices: Option<&mut [u8]>, rates: bool) -> Scoped {
    match inputs.len() {
        0 | 1 => scope_run_n::<1>(tape, every, ctx, inputs, b, choices, rates),
        2 => scope_run_n::<2>(tape, every, ctx, inputs, b, choices, rates),
        3 | 4 => scope_run_n::<4>(tape, every, ctx, inputs, b, choices, rates),
        _ => scope_run_n::<MAX_INPUTS>(tape, every, ctx, inputs, b, choices, rates),
    }
}

/// What scope_run bounds, per slot: the value over the box and scope (`v`),
/// at the inputs' centre (`m`), and its rates along them (`d`, `inputs.len()`
/// per slot, when asked for: 0 where the value does not depend on them,
/// ENTIRE where unbounded).
pub struct Scoped {
    pub v: Vec<I>,
    pub m: Vec<I>,
    pub d: Vec<I>,
}

/// scope_run over every value of the tape's SSA form (a max's terms: the
/// clearance proofs' first-order forms, ink-kernel's `centred`), with the
/// rates along `inputs` where `rates`.
pub fn eval_interval_scope_every(tape: &Tape, ctx: &Ctx, inputs: &Inputs, b: [I; 3], rates: bool) -> Scoped {
    scope_run(tape, true, ctx, inputs, b, None, rates)
}

/// scope_run with N rates per value (N at least inputs.len(); the rest are
/// 0).
fn scope_run_n<const N: usize>(tape: &Tape, every: bool, ctx: &Ctx, inputs: &Inputs, b: [I; 3], choices: Option<&mut [u8]>, rates: bool) -> Scoped {
    let (code, n) = if every { (&tape.ssa[..], tape.ssa.len()) } else { (&tape.code[..], tape.nslots) };
    // The branches decided, the caller's or our own.
    let mut own;
    let ch: &mut [u8] = match choices {
        Some(c) => c,
        None => {
            own = vec![crate::tape::BOTH; code.len()];
            &mut own
        }
    };
    let zero = I::point(0.0);
    let mut run = Scope::<N> {
        tape,
        code,
        ctx,
        mid: inputs.centre(ctx),
        inputs,
        hw: inputs.half_widths(ctx),
        b,
        v: vec![zero; n],
        m: vec![zero; n],
        d: vec![None; n],
        crossed: vec![false; n],
        slopes: Slopes::of(tape, n),
        every,
    };
    forked(tape, &mut run, 0, code.len(), Some(ch));
    let k = inputs.len();
    let mut flat = Vec::with_capacity(if rates { n * k } else { 0 });
    if rates {
        for x in &run.d {
            flat.extend_from_slice(&x.unwrap_or([zero; N])[..k]);
        }
    }
    Scoped { v: run.v, m: run.m, d: flat }
}

/// scope_run's state.
struct Scope<'a, const N: usize> {
    tape: &'a Tape,
    code: &'a [crate::tape::Ins],
    ctx: &'a Ctx,
    /// ctx at the inputs' centre. Ranged params that are not inputs stay
    /// ranged: the value at the centre must hold for every value of them.
    mid: Ctx,
    inputs: &'a Inputs,
    /// |u_j - c_j| over the scope, rounded up.
    hw: [f64; MAX_INPUTS],
    b: [I; 3],
    /// Per slot: the bound over the scope; a bound at the centre (where a
    /// rate is unbounded, the bound over the scope, which holds there too);
    /// and the rates over the scope (None where the value does not depend
    /// on the inputs, ENTIRE where unbounded).
    v: Vec<I>,
    m: Vec<I>,
    d: Vec<Option<[I; N]>>,
    /// Per slot: whether the inputs may cross a jump it is bounded on one
    /// side of. Each side's rates are bounded; the value itself may jump
    /// along them.
    crossed: Vec<bool>,
    /// Gradient bounds hold over the whole scope (each partial at every
    /// moment and param value), so they tighten the bound over it.
    slopes: Option<Slopes>,
    every: bool,
}

/// How far the mean value form reaches from a value at the inputs' centre
/// by rates `du`: the sum of |d/du_j| hw_j (u_j - c_j lies in [-hw_j,
/// hw_j]), rounded up; an input whose rate is exactly 0 adds nothing.
/// Infinite where a rate is unbounded.
fn reach<const N: usize>(du: &[I; N], hw: &[f64; MAX_INPUTS], k: usize) -> f64 {
    (0..k).fold(0.0, |e, j| {
        let r = interval::mag(du[j]);
        if r == 0.0 { e } else { interval::add_up(e, interval::mul_up(r, hw[j])) }
    })
}

/// A bound `over` met with a value at the centre `m` within `e` (both hold
/// every value, so they meet).
fn around(over: I, m: I, e: f64) -> I {
    interval::meet(over, I::new((m.lo - e).next_down(), (m.hi + e).next_up()))
}

/// What a value's slot holds in scope_run: its bounds over the scope and at
/// the centre, rates, whether the inputs may cross a jump it is on one side
/// of, and its gradient bounds.
type ScopeSlot<const N: usize> = (I, I, Option<[I; N]>, bool, Option<interval::G>);

impl<const N: usize> Fork for Scope<'_, N> {
    type Out = ScopeSlot<N>;
    /// Each argument's bounds over the scope and at the centre, rates, and
    /// whether the inputs may cross a jump it is on one side of.
    type Args = [(I, I, Option<[I; N]>, bool); 2];

    #[inline(always)]
    fn step(&mut self, i: usize, choices: Option<&mut [u8]>) -> I {
        let (tape, ctx, ins) = (self.tape, self.ctx, &self.code[i]);
        let ch = choices.expect("scope_run passes its choices");
        let o = ins.out as usize;
        let mut over = interval_ins(tape, ctx, i, ins, &self.v, self.b, Some(&mut *ch));
        if let Some(sl) = &mut self.slopes {
            over = sl.tighten(tape, ctx, i, ins, self.b, over);
        }
        let (args, na) = ins.args();
        let args = &args[..na];
        let axis = match ins.op {
            Op::Var(x) => self.inputs.axis(x),
            _ => None,
        };
        let (v, d) = (&self.v, &self.d);
        if fixed(tape, ctx, i) || !(axis.is_some() || args.iter().any(|&k| d[k as usize].is_some())) {
            (self.v[o], self.m[o], self.d[o], self.crossed[o]) = (over, over, None, false);
            return over;
        }
        let zero = I::point(0.0);
        let (a, bb, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
        let da = |k: usize| d[k].unwrap_or([zero; N]);
        let crossed = args.iter().any(|&k| self.crossed[k as usize]);
        let du = if let Some(j) = axis {
            Some(std::array::from_fn(|k| I::point(if k == j { 1.0 } else { 0.0 })))
        } else if args.iter().all(|&k| da(k as usize) == [zero; N]) {
            // A function of values constant in the inputs is constant in them.
            Some([zero; N])
        } else {
            match ins.op {
                Op::Un(op) => interval::d_un(op, v[a], over, da(a)),
                Op::Bin(op) => interval::d_bin(op, v[a], v[bb], over, ch[i], da(a), da(bb)),
                Op::Select => interval::d_select(ch[i], same_branches(tape, i), da(bb), da(c)),
                // A mark of `a` as the derivative along t of `b`: where b
                // may jump (its rate along t has no bound), so may a, by
                // the jump's impulse, whatever a's formula. (One side of a
                // jump t may cross has bounded rates, but every value read
                // after it is joined with unbounded ones: `join`.)
                Op::Deriv(3) => {
                    let mut r = da(a);
                    if let Some(j) = self.inputs.axis(Var::T)
                        && !(da(bb)[j].lo.is_finite() && da(bb)[j].hi.is_finite())
                    {
                        r[j] = I::ENTIRE;
                    }
                    Some(r)
                }
                // The value of `a` (a promise bounds slopes in space, not
                // rates in time or params).
                Op::Tag(_) | Op::Lipschitz | Op::Deriv(_) | Op::Read | Op::Real => Some(da(a)),
                Op::Grid(gi, mk) => Some(interval::d_grid(&tape.grids[gi as usize], mk, [v[a], v[bb], v[c]], [da(a), da(bb), da(c)])),
                Op::Frame => Some([zero; N]),
                Op::Const(_) | Op::Var(_) => unreachable!(),
            }
        };
        // Along an input no argument moves along, neither does the value,
        // whatever the rule gives along the others: a bore's sqrt, its slope
        // unbounded on the axis, has no rate along the bore's diameter. (An
        // unbounded rate along another input leaves the value bounded over
        // the scope alone, but what reads it keeps its rate along this one,
        // 0 or the diameter's: a search halves only the inputs that blur
        // it, not the diameter's band to its finest beside a turn.)
        let mut du = du.unwrap_or([I::ENTIRE; N]);
        for (j, r) in du.iter_mut().enumerate() {
            if axis.is_none() && args.iter().all(|&k| da(k as usize)[j] == zero) {
                *r = zero;
            }
        }
        let e = reach(&du, &self.hw, self.inputs.len());
        if !e.is_finite() {
            (self.v[o], self.m[o], self.d[o], self.crossed[o]) = (over, over, Some(du), crossed);
            return over;
        }
        let at_mid = interval_ins(tape, &self.mid, i, ins, &self.m, self.b, None);
        let r = around(over, at_mid, e);
        (self.v[o], self.m[o], self.d[o], self.crossed[o]) = (r, at_mid, Some(du), crossed);
        r
    }
    fn args(&self, i: usize) -> Self::Args {
        let ins = &self.code[i];
        [ins.a, ins.b].map(|k| {
            let k = k as usize;
            (self.v[k], self.m[k], self.d[k], self.crossed[k])
        })
    }
    fn bounds(args: &Self::Args) -> (I, I) {
        (args[0].0, args[1].0)
    }
    fn ctx(&self) -> &Ctx {
        self.ctx
    }
    fn side(&mut self, i: usize, args: &Self::Args, k: usize, side: I) {
        let o = self.code[i].out as usize;
        let zero = [I::point(0.0); N];
        // Whether an argument moves with the inputs (or the inputs may
        // cross a jump it is on one side of).
        let moves = |&(_, _, d, c): &(I, I, Option<[I; N]>, bool)| c || d.is_some_and(|d| d != zero);
        let [y, x] = args;
        let (v, m, d, crossed) = match self.code[i].op {
            // One integer over the region: constant in the inputs.
            Op::Un(UnOp::Floor) => (side, side, y.2.map(|_| zero), moves(y)),
            // The branch continued over the region, at the centre from its
            // arguments' bounds there (within those over the scope: both
            // hold there), its rates atan2's.
            Op::Bin(BinOp::Atan2) if y.2.is_none() && x.2.is_none() => {
                let v = interval::atan2_branch(y.0, x.0, k == 1);
                (v, v, None, false)
            }
            Op::Bin(BinOp::Atan2) => {
                let up = k == 1;
                let v = interval::atan2_branch(y.0, x.0, up);
                let (yc, xc) = (interval::meet(y.1, y.0), interval::meet(x.1, x.0));
                let m = if xc.hi < 0.0 { interval::meet(v, interval::atan2_branch(yc, xc, up)) } else { v };
                let du = interval::d_atan2(y.0, x.0, y.2.unwrap_or(zero), x.2.unwrap_or(zero)).unwrap_or([I::ENTIRE; N]);
                let e = reach(&du, &self.hw, self.inputs.len());
                let crossed = moves(y) | moves(x);
                if e.is_finite() { (around(v, m, e), m, Some(du), crossed) } else { (v, v, Some([I::ENTIRE; N]), crossed) }
            }
            _ => unreachable!("jump_sides takes floors and atan2"),
        };
        (self.v[o], self.m[o], self.d[o], self.crossed[o]) = (v, m, d, crossed);
        // This side's range, not a meet: the bound across the cut is not
        // of this branch.
        if let Some(sl) = &mut self.slopes {
            sl.0[o].v = v;
        }
    }
    fn out(&self, i: usize) -> Self::Out {
        let o = self.code[i].out as usize;
        (self.v[o], self.m[o], self.d[o], self.crossed[o], self.slopes.as_ref().map(|sl| sl.0[o]))
    }
    fn put(&mut self, i: usize, (v, m, d, c, g): Self::Out) {
        let o = self.code[i].out as usize;
        (self.v[o], self.m[o], self.d[o], self.crossed[o]) = (v, m, d, c);
        if let (Some(sl), Some(g)) = (&mut self.slopes, g) {
            sl.0[o] = g;
        }
    }
    fn join(&mut self, i: usize, a: Self::Out, b: Self::Out) {
        let o = self.code[i].out as usize;
        let (v, crossed, zero) = (interval::hull(a.0, b.0), a.3 || b.3, [I::point(0.0); N]);
        let (m, d) = match (a.2, b.2) {
            (None, None) => (interval::hull(a.1, b.1), None),
            // Each point takes one side for every value of the inputs: its
            // rates are that side's.
            (da, db) if !crossed => {
                let (da, db) = (da.unwrap_or(zero), db.unwrap_or(zero));
                (interval::hull(a.1, b.1), Some(std::array::from_fn(|j| interval::hull(da[j], db[j]))))
            }
            // It may jump along the inputs.
            _ => (v, Some([I::ENTIRE; N])),
        };
        (self.v[o], self.m[o], self.d[o], self.crossed[o]) = (v, m, d, crossed);
        if let (Some(sl), Some(ga), Some(gb)) = (&mut self.slopes, a.4, b.4) {
            sl.0[o] = ga.hull(&gb);
        }
    }
    fn every_value(&self) -> bool {
        self.every
    }
}

/// Value and spatial gradient `[f, df/dx, df/dy, df/dz]` at each point.
pub fn eval_grads(tape: &Tape, ctx: &Ctx, pts: &[V3]) -> Vec<[f64; 4]> {
    // A few points (a cell's samples) in narrow lanes: the slots are
    // allocated, and each op run, per lane.
    if pts.len() <= 8 { grads_in::<8>(tape, ctx, pts) } else { grads_in::<L>(tape, ctx, pts) }
}

/// `eval_grads` in chunks of N points.
fn grads_in<const N: usize>(tape: &Tape, ctx: &Ctx, pts: &[V3]) -> Vec<[f64; 4]> {
    let mut out = Vec::with_capacity(pts.len());
    let mut s = vec![[[0.0f64; N]; 4]; tape.nslots];
    for chunk in pts.chunks(N) {
        for (i, ins) in tape.code.iter().enumerate() {
            let o = ins.out as usize;
            match ins.op {
                Op::Const(v) => s[o] = [[v; N], [0.0; N], [0.0; N], [0.0; N]],
                Op::Var(v @ (Var::X | Var::Y | Var::Z)) => {
                    let axis = match v {
                        Var::X => 0,
                        Var::Y => 1,
                        _ => 2,
                    };
                    let r = &mut s[o];
                    *r = [[0.0; N]; 4];
                    for (k, p) in chunk.iter().enumerate() {
                        r[0][k] = p[axis];
                    }
                    r[axis + 1] = [1.0; N];
                }
                Op::Var(Var::T) => s[o] = [[ctx.t; N], [0.0; N], [0.0; N], [0.0; N]],
                Op::Var(Var::Param(p)) => s[o] = [[ctx.params[p as usize]; N], [0.0; N], [0.0; N], [0.0; N]],
                Op::Un(op) => {
                    let a = s[ins.a as usize];
                    let r = &mut s[o];
                    for k in 0..N {
                        let (v, fa) = dun(op, a[0][k]);
                        r[0][k] = v;
                        for j in 1..4 {
                            r[j][k] = mul0(fa, a[j][k]);
                        }
                    }
                }
                Op::Bin(op) => {
                    let a = s[ins.a as usize];
                    let b = s[ins.b as usize];
                    let r = &mut s[o];
                    for k in 0..N {
                        let (v, fa, fb) = dbin(op, a[0][k], b[0][k]);
                        r[0][k] = v;
                        for j in 1..4 {
                            r[j][k] = defined(mul0(fa, a[j][k]) + mul0(fb, b[j][k]));
                        }
                    }
                }
                Op::Select => {
                    let c = s[ins.a as usize];
                    let a = s[ins.b as usize];
                    let b = s[ins.c as usize];
                    let r = &mut s[o];
                    for k in 0..N {
                        let src = if c[0][k] < 0.0 { &a } else { &b };
                        for j in 0..4 {
                            r[j][k] = src[j][k];
                        }
                    }
                }
                Op::Read => {
                    // An op of its coordinates: its partials (its form's
                    // gradient there) chained with theirs.
                    let f = &tape.code[tape.ssa[i].b as usize];
                    let (cx, cy, cz) = (s[f.a as usize], s[f.b as usize], s[f.c as usize]);
                    let at: Vec<V3> = (0..chunk.len()).map(|k| [cx[0][k], cy[0][k], cz[0][k]]).collect();
                    let gr = eval_grads(&tape.forms[ins.c as usize], ctx, &at);
                    let v = s[ins.a as usize][0];
                    let r = &mut s[o];
                    r[0] = v;
                    for (k, gr) in gr.iter().enumerate() {
                        for j in 1..4 {
                            r[j][k] = chain([gr[1], gr[2], gr[3]], [cx[j][k], cy[j][k], cz[j][k]]);
                        }
                    }
                }
                Op::Tag(_) | Op::Lipschitz | Op::Deriv(_) | Op::Real => s[o] = s[ins.a as usize],
                Op::Frame => s[o] = [[0.0; N]; 4],
                Op::Grid(gi, m) => {
                    // Chain rule through the coordinate expressions.
                    let (cx, cy, cz) = (s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]);
                    let g = &tape.grids[gi as usize];
                    let r = &mut s[o];
                    for k in 0..N {
                        let (v, gr) = g.component_grad([cx[0][k], cy[0][k], cz[0][k]], m);
                        r[0][k] = v;
                        for j in 1..4 {
                            r[j][k] = chain(gr, [cx[j][k], cy[j][k], cz[j][k]]);
                        }
                    }
                }
            }
        }
        let r = &s[tape.out_slot()];
        for k in 0..chunk.len() {
            out.push([r[0][k], r[1][k], r[2][k], r[3][k]]);
        }
    }
    out
}

/// The chain rule's product: `0 * inf = 0` (a zero partial means "no
/// dependence"), as the value rules have it.
#[inline]
fn mul0(f: f64, d: f64) -> f64 {
    if f == 0.0 || d == 0.0 { 0.0 } else { f * d }
}

/// The chain rule through three coordinates: partials `f` times their
/// derivatives `d`, summed in order under the value rules (each sum
/// defined), as graph::Graph::diff writes it.
#[inline]
fn chain(f: V3, d: V3) -> f64 {
    defined(defined(mul0(f[0], d[0]) + mul0(f[1], d[1])) + mul0(f[2], d[2]))
}

/// Value and every partial derivative at one point:
/// `[f, d/dx, d/dy, d/dz, d/dt, d/dparam0, ...]`.
pub fn eval_dual(tape: &Tape, ctx: &Ctx, p: V3) -> Vec<f64> {
    let nd = 4 + ctx.params.len();
    let w = nd + 1;
    let (mut s, mut r) = (vec![0.0f64; tape.nslots * w], vec![0.0; w]);
    for (i, ins) in tape.code.iter().enumerate() {
        let o = ins.out as usize * w;
        r.fill(0.0);
        match ins.op {
            Op::Const(v) => r[0] = v,
            Op::Var(v) => {
                let (val, dir) = match v {
                    Var::X => (p[0], 0),
                    Var::Y => (p[1], 1),
                    Var::Z => (p[2], 2),
                    Var::T => (ctx.t, 3),
                    Var::Param(i) => (ctx.params[i as usize], 4 + i as usize),
                };
                r[0] = val;
                r[1 + dir] = 1.0;
            }
            Op::Un(op) => {
                let a = &s[ins.a as usize * w..][..w];
                let (v, fa) = dun(op, a[0]);
                r[0] = v;
                for j in 1..w {
                    r[j] = mul0(fa, a[j]);
                }
            }
            Op::Bin(op) => {
                let a = &s[ins.a as usize * w..][..w];
                let b = &s[ins.b as usize * w..][..w];
                let (v, fa, fb) = dbin(op, a[0], b[0]);
                r[0] = v;
                for j in 1..w {
                    r[j] = defined(mul0(fa, a[j]) + mul0(fb, b[j]));
                }
            }
            Op::Select => {
                let c = s[ins.a as usize * w];
                let src = if c < 0.0 { ins.b } else { ins.c } as usize * w;
                r.copy_from_slice(&s[src..src + w]);
            }
            Op::Read => {
                // As eval_grads.
                let f = &tape.code[tape.ssa[i].b as usize];
                let (a, b, c) = (&s[f.a as usize * w..][..w], &s[f.b as usize * w..][..w], &s[f.c as usize * w..][..w]);
                let gr = eval_grads(&tape.forms[ins.c as usize], ctx, &[[a[0], b[0], c[0]]])[0];
                r[0] = s[ins.a as usize * w];
                for j in 1..w {
                    r[j] = chain([gr[1], gr[2], gr[3]], [a[j], b[j], c[j]]);
                }
            }
            Op::Tag(_) | Op::Lipschitz | Op::Deriv(_) | Op::Real => r.copy_from_slice(&s[ins.a as usize * w..][..w]),
            Op::Frame => {}
            Op::Grid(gi, m) => {
                let (a, b, c) = (&s[ins.a as usize * w..][..w], &s[ins.b as usize * w..][..w], &s[ins.c as usize * w..][..w]);
                let (v, gr) = tape.grids[gi as usize].component_grad([a[0], b[0], c[0]], m);
                r[0] = v;
                for j in 1..w {
                    r[j] = chain(gr, [a[j], b[j], c[j]]);
                }
            }
        }
        s[o..o + w].copy_from_slice(&r);
    }
    let o = tape.out_slot() * w;
    s[o..o + w].to_vec()
}

/// Construction names (graph tag ids) on the active path at `p`, outermost
/// first. From the root, follow the branch of each min/max/select that decides
/// the value (ties take the first operand), and every operand of other ops,
/// collecting the names passed. On a surface this says which named
/// construction the point belongs to: implicit faces without B-rep.
pub fn provenance(tape: &Tape, ctx: &Ctx, p: V3) -> Vec<u32> {
    let v = values(tape, ctx, p);
    names_on(tape, |ins| {
        let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
        match ins.op {
            Op::Bin(BinOp::Min) => Follow::One(if v[a] <= v[b] || v[b].is_nan() { a } else { b }),
            Op::Bin(BinOp::Max | BinOp::Emax) => Follow::One(if v[a] >= v[b] || v[b].is_nan() { a } else { b }),
            Op::Select => Follow::One(if v[a] < 0.0 { b } else { c }),
            _ => Follow::All,
        }
    })
}

/// The names the value at `p` depends on: provenance's walk, but a min or
/// max follows its losing side too where that side comes within `reach` (a
/// distance) of the winner's value by its own slope (both are on the
/// surface there: at an edge between them, or on a face they share), and
/// a product does not follow a factor whose co-factor is 0 (a blend's side
/// weighted out). Slopes are measured `reach` along each corner of a
/// tetrahedron about p, so scaling either side, or the field, by any
/// positive factor reads the same.
pub fn decided_names(tape: &Tape, ctx: &Ctx, p: V3, reach: f64) -> Vec<u32> {
    let v = values(tape, ctx, p);
    let around: Vec<Vec<f64>> =
        [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]].iter().map(|d: &V3| values(tape, ctx, std::array::from_fn(|k| p[k] + reach * d[k]))).collect();
    // At least the slope of instruction i's value (each corner is sqrt 3
    // reach away; a gradient's component along one of them is at least
    // 1 / sqrt 3 of it).
    let slope = |i: usize| around.iter().map(|w| (w[i] - v[i]).abs()).fold(0.0, f64::max) / reach;
    names_on(tape, |ins| {
        let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
        let side = |first: bool| {
            let lost = if first { b } else { a };
            if !((v[a] - v[b]).abs() > reach * slope(lost)) { Follow::All } else { Follow::One(if first { a } else { b }) }
        };
        match ins.op {
            Op::Bin(BinOp::Min) => side(v[a] < v[b]),
            Op::Bin(BinOp::Max | BinOp::Emax) => side(v[a] > v[b]),
            Op::Bin(BinOp::Mul) => match (v[a] != 0.0, v[b] != 0.0) {
                (true, true) => Follow::All,
                (false, false) => Follow::None,
                (a_on, _) => Follow::One(if a_on { b } else { a }),
            },
            Op::Select => Follow::One(if v[a] < 0.0 { b } else { c }),
            _ => Follow::All,
        }
    })
}

/// Every instruction's value at `p`.
fn values(tape: &Tape, ctx: &Ctx, p: V3) -> Vec<f64> {
    let n = tape.ssa.len();
    let mut v = vec![0.0f64; n];
    for (i, ins) in tape.ssa.iter().enumerate() {
        let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
        v[i] = match ins.op {
            Op::Const(x) => x,
            Op::Var(Var::X) => p[0],
            Op::Var(Var::Y) => p[1],
            Op::Var(Var::Z) => p[2],
            Op::Var(Var::T) => ctx.t,
            Op::Var(Var::Param(k)) => ctx.params[k as usize],
            Op::Un(op) => un_f64(op, v[a]),
            Op::Bin(op) => bin_f64(op, v[a], v[b]),
            Op::Select => {
                if v[a] < 0.0 {
                    v[b]
                } else {
                    v[c]
                }
            }
            Op::Tag(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real => v[a],
            Op::Frame => 0.0,
            Op::Grid(gi, m) => tape.grids[gi as usize].component([v[a], v[b], v[c]], m),
        };
    }
    v
}

/// The names on every path the value can take: `provenance` with both
/// branches of each min, max and select the tape holds. For a tape
/// simplified over a region (Tape::simplify drops a branch interval
/// arithmetic proved cannot decide the value anywhere in it), every name
/// whose surface lies in the region, and maybe others.
pub fn names_on_paths(tape: &Tape) -> Vec<u32> {
    names_on(tape, |_| Follow::All)
}

/// Which operands of a min, max, select or product the walk follows.
enum Follow {
    All,
    One(usize),
    None,
}

/// The walk from the root `provenance`, `decided_names` and
/// `names_on_paths` share: `pick` gives the operands an op follows (All:
/// every operand, but a select's condition, which names nothing), and a
/// name provenance stops at (Tape::opaque) hides those inside it.
fn names_on(tape: &Tape, pick: impl Fn(&crate::tape::Ins) -> Follow) -> Vec<u32> {
    let n = tape.ssa.len();
    let mut out = Vec::new();
    let mut seen = vec![false; n];
    let mut stack = vec![n - 1];
    while let Some(i) = stack.pop() {
        if seen[i] {
            continue;
        }
        seen[i] = true;
        let ins = &tape.ssa[i];
        let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
        match pick(ins) {
            Follow::One(k) => {
                stack.push(k);
                continue;
            }
            Follow::None => continue,
            Follow::All => {}
        }
        // Push the later operands first, so the first one's names come out first.
        match ins.op {
            Op::Tag(t) => {
                if !out.contains(&t) {
                    out.push(t);
                }
                if !tape.opaque.contains(&t) {
                    stack.push(a);
                }
            }
            Op::Select => stack.extend([c, b]),
            Op::Un(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real => stack.push(a),
            Op::Grid(..) => stack.extend([c, b, a]),
            Op::Bin(_) => stack.extend([b, a]),
            Op::Const(_) | Op::Var(_) | Op::Frame => {}
        }
    }
    out
}

/// An undefined case an evaluation passed through (docs/IR.md, Fields, Value rules: its
/// value is 0, and it is reported, never silent): what it was, where in the
/// tape (an ssa index), the nearest named construction around it, how many
/// of the points hit it, and the first such point.
#[derive(Clone, Debug)]
pub struct Undefined {
    pub what: &'static str,
    pub ins: usize,
    pub name: Option<u32>,
    pub count: usize,
    pub at: V3,
}

/// The undefined cases evaluating the tape at `pts` passes through, one per
/// kind and instruction, most frequent first. Only cases on the active path
/// count: those the value at the point depends on (not the branch a select
/// did not take, nor the losing side of a min or max), so a guarded
/// division is not reported.
pub fn undefined_cases(tape: &Tape, ctx: &Ctx, pts: &[V3]) -> Vec<Undefined> {
    let n = tape.ssa.len();
    // The nearest tag above each instruction (consumers come after operands).
    let mut owner: Vec<Option<u32>> = vec![None; n];
    for i in (0..n).rev() {
        let ins = &tape.ssa[i];
        let mine = if let Op::Tag(t) = ins.op { Some(t) } else { owner[i] };
        let (args, k) = ins.args();
        for &a in &args[..k] {
            if owner[a as usize].is_none() {
                owner[a as usize] = mine;
            }
        }
    }
    let mut found: Vec<Undefined> = Vec::new();
    let mut v = vec![0.0f64; n];
    let mut bad: Vec<Option<&'static str>> = vec![None; n];
    let mut active = vec![false; n];
    for &p in pts {
        for (i, ins) in tape.ssa.iter().enumerate() {
            let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
            let what = match ins.op {
                Op::Bin(op) => crate::graph::undefined_bin(op, v[a], v[b]),
                Op::Un(op) => crate::graph::undefined_un(op, v[a]),
                _ => None,
            };
            bad[i] = what;
            v[i] = match ins.op {
                Op::Const(x) => x,
                Op::Var(Var::X) => p[0],
                Op::Var(Var::Y) => p[1],
                Op::Var(Var::Z) => p[2],
                Op::Var(Var::T) => ctx.t,
                Op::Var(Var::Param(k)) => ctx.params[k as usize],
                Op::Un(op) => un_f64(op, v[a]),
                Op::Bin(op) => bin_f64(op, v[a], v[b]),
                Op::Select => {
                    if v[a] < 0.0 {
                        v[b]
                    } else {
                        v[c]
                    }
                }
                Op::Tag(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real => v[a],
                Op::Frame => 0.0,
                Op::Grid(gi, m) => tape.grids[gi as usize].component([v[a], v[b], v[c]], m),
            };
        }
        // The active path, from the output (consumers come after operands).
        active.iter_mut().for_each(|x| *x = false);
        active[n - 1] = true;
        for i in (0..n).rev() {
            if !active[i] {
                continue;
            }
            let ins = &tape.ssa[i];
            let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
            match ins.op {
                Op::Select => {
                    active[a] = true;
                    active[if v[a] < 0.0 { b } else { c }] = true;
                }
                Op::Bin(BinOp::Min) => {
                    active[a] |= v[a] <= v[b];
                    active[b] |= v[b] <= v[a];
                }
                Op::Bin(BinOp::Max) => {
                    active[a] |= v[a] >= v[b];
                    active[b] |= v[b] >= v[a];
                }
                _ => {
                    let (args, k) = ins.args();
                    args[..k].iter().for_each(|&x| active[x as usize] = true);
                }
            }
        }
        for i in 0..n {
            if let (true, Some(what)) = (active[i], bad[i]) {
                match found.iter_mut().find(|u| u.ins == i && u.what == what) {
                    Some(u) => u.count += 1,
                    None => found.push(Undefined { what, ins: i, name: owner[i], count: 1, at: p }),
                }
            }
        }
    }
    found.sort_by_key(|u| std::cmp::Reverse(u.count));
    found
}

/// Value and gradient bounds over a box (see interval::G). The gradient bound
/// is a Lipschitz constant for the field there: |f(p) - f(q)| <= L |p - q|.
pub fn eval_interval_grad(tape: &Tape, ctx: &Ctx, b: [I; 3]) -> interval::G {
    grad_run(tape, &tape.code, tape.nslots, ctx, b)[tape.out_slot()]
}

/// `eval_interval_grad`'s bounds of every value of the tape's SSA form (a
/// max's terms: the clearance proofs' first-order forms, ink-kernel's
/// `centred`).
pub fn eval_interval_grad_every(tape: &Tape, ctx: &Ctx, b: [I; 3]) -> Vec<interval::G> {
    grad_run(tape, &tape.ssa, tape.ssa.len(), ctx, b)
}

/// Whether select instruction `i`'s two branches are the same field, up to
/// names (tags): then it is continuous where its condition changes sign.
fn same_branches(tape: &Tape, i: usize) -> bool {
    let ssa = &tape.ssa;
    let strip = |mut j: usize| {
        while let Op::Tag(_) = ssa[j].op {
            j = ssa[j].a as usize;
        }
        j
    };
    strip(ssa[i].b as usize) == strip(ssa[i].c as usize)
}

/// Interval forward-mode AD over `code` (the tape's slot code, or its SSA
/// form to keep every value): the bounds in every slot.
fn grad_run(tape: &Tape, code: &[crate::tape::Ins], nslots: usize, ctx: &Ctx, b: [I; 3]) -> Vec<interval::G> {
    let mut s = vec![interval::G::constant(0.0); nslots];
    let mut f = vec![0u8; nslots];
    for (i, ins) in code.iter().enumerate() {
        let g = grad_ins(tape, ctx, i, ins, &s, &f, b);
        f[ins.out as usize] = real_flags(ins, &s, &f, &g);
        s[ins.out as usize] = g;
    }
    s
}

// What `real_flags` says of a value: it reads x, y or z; or, in space alone,
// it is not the real number it stands for.
const SPATIAL: u8 = 1;
const UNREAL: u8 = 2;

/// The flags of instruction `ins`'s value `g`, from its arguments' in `f`.
/// A promise is made for its constants as computed, in real numbers: a
/// value that does not read x, y or z and is infinite, or comes of an
/// undefined case (inf / inf, a quotient by a range through 0), is not one,
/// and nor is what is computed from it. Its params' and t's values are this
/// query's, a warp's substitutions included, so the check is of the
/// numbers the promise is used with (a param at 1e200 made a k of
/// 1 + inf / inf, 1.001 at every real value, 0.001).
fn real_flags(ins: &crate::tape::Ins, s: &[interval::G], f: &[u8], g: &interval::G) -> u8 {
    let (args, n) = ins.args();
    let of = |k: usize| f[args[k] as usize];
    let mut out = match ins.op {
        Op::Var(Var::X | Var::Y | Var::Z) => return SPATIAL,
        // Of the taken branch, where the condition decides one.
        Op::Select => {
            let c = s[ins.a as usize].v;
            of(0)
                | if c.hi < 0.0 {
                    of(1)
                } else if c.lo >= 0.0 {
                    of(2)
                } else {
                    of(1) | of(2)
                }
        }
        _ => (0..n).fold(0, |acc, k| acc | of(k)),
    };
    if out & SPATIAL == 0 {
        let arg = |k: usize| s[args[k] as usize].v;
        let undefined = match ins.op {
            Op::Bin(BinOp::Div) => arg(1).lo <= 0.0 && arg(1).hi >= 0.0,
            Op::Bin(BinOp::Pow) => arg(0).lo < 0.0 && !(arg(1).lo == arg(1).hi && arg(1).lo.fract() == 0.0),
            Op::Un(UnOp::Ln) => arg(0).lo < 0.0,
            _ => false,
        };
        if undefined || !(g.v.lo.is_finite() && g.v.hi.is_finite()) {
            out |= UNREAL;
        }
    }
    out
}

/// Instruction `i`'s value and gradient bounds (interval::G) from its
/// arguments' in `s` (and their flags, `real_flags`, in `f`).
fn grad_ins(tape: &Tape, ctx: &Ctx, i: usize, ins: &crate::tape::Ins, s: &[interval::G], f: &[u8], b: [I; 3]) -> interval::G {
    use interval::G;
    if fixed(tape, ctx, i) {
        return G::constant(fixed_value(ins, |k| s[k as usize].v.lo, ctx));
    }
    match ins.op {
        Op::Const(v) => G::constant(v),
        Op::Var(Var::X) => G::coordinate(b[0], 0),
        Op::Var(Var::Y) => G::coordinate(b[1], 1),
        Op::Var(Var::Z) => G::coordinate(b[2], 2),
        Op::Var(Var::T) => G { v: ctx.t_span.map_or(I::point(ctx.t), |(a, b)| I::new(a, b)), d: [I::point(0.0); 3], l: 0.0, s2: None, m: interval::NO_GRAM },
        // Constant in space at every value of its span (a nonnegative
        // constant is a trivial sum of squares, as in G::constant).
        Op::Var(Var::Param(p)) => {
            let v = ctx.param(p);
            G { v, d: [I::point(0.0); 3], l: 0.0, s2: (v.lo >= 0.0).then_some(0.0), m: interval::NO_GRAM }
        }
        Op::Un(op) => interval::g_un(op, s[ins.a as usize]),
        // A value times itself is its square (as interval_ins): so a sum
        // of products a a (std's norm, sqrt (dot a a)) is a sum of
        // squares, whose square root's slope is bounded (G::s2, G::m).
        Op::Bin(BinOp::Mul) if ins.a == ins.b => interval::g_un(UnOp::Square, s[ins.a as usize]),
        Op::Bin(op) => interval::g_bin(op, s[ins.a as usize], s[ins.b as usize]),
        Op::Select => interval::g_select(s[ins.a as usize], s[ins.b as usize], s[ins.c as usize], same_branches(tape, i)),
        Op::Tag(_) | Op::Read | Op::Real => s[ins.a as usize],
        // Trusted only where its field's, k's and frame's constants are real.
        Op::Lipschitz if (f[ins.a as usize] | f[ins.b as usize] | f[ins.c as usize]) & UNREAL != 0 => s[ins.a as usize],
        Op::Lipschitz => interval::g_lipschitz(s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]),
        Op::Frame => interval::g_frame([s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]]),
        // Sampled fields: value bounds only (no derivative bounds yet).
        Op::Grid(gi, m) => interval::g_grid(&tape.grids[gi as usize], m, [s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]]),
        // A derivative: its formula's bounds, its value also within the
        // field's partial along the axis (along t, the formula's alone).
        Op::Deriv(k) => {
            let a = s[ins.a as usize];
            if k < 3 { G { v: interval::meet(a.v, s[ins.b as usize].d[k as usize]), ..a } } else { a }
        }
    }
}

/// Gradient bounds beside an evaluator's value bounds, for the instructions
/// whose bounds a derivative's range reads (Tape::slopes): a derivative of
/// a field along x, y or z lies within the field's partial along it
/// (interval::G), which keeps relations the formula's own bounds lose.
/// The partial of sqrt(u1^2 + u2^2 + ...) is at most sqrt(sum (du_i/dx)^2)
/// (Cauchy-Schwarz), so d/dx of a distance to a point is within [-1, 1]
/// over any box, where its formula, x / sqrt(x^2 + y^2) near the point, is
/// a quotient of two ranges through 0.
///
/// Like every range, both bounds are of the derivative's value in real
/// numbers, so they meet. The partials bound the derivative rules' real
/// value at every point, their conventions too (a tie's mean lies in the
/// hull of both sides, and a quotient by 0 is 0, inside every bound that
/// holds near it), and never rest on Lipschitz promises. Not the f64
/// evaluator's value: at (2.7e-162, 0) the squares in the formula above
/// underflow and floats give 1.21, where the slope is 1. The formula's
/// parts of constants and params alone mean their f64 values, and the
/// partials hold those too: each is one float operation on values its
/// operands' bounds hold, or (the product rule's sum of products) the
/// operations the partials take, rounded outward; either way inside, since
/// one rounding of a value inside a bound with float ends stays inside it.
/// The exception is a tie's mean, da/2 + db/2, whose halvings round below
/// 2^-1021: interval::tie_slack makes room for it.
struct Slopes(Vec<interval::G>, Vec<u8>);

impl Slopes {
    /// For an evaluator over `n` slots (the tape's code, or its SSA form).
    fn of(tape: &Tape, n: usize) -> Option<Slopes> {
        (!tape.slopes.is_empty()).then(|| Slopes(vec![interval::G::constant(0.0); n], vec![0; n]))
    }

    /// Instruction `i`'s bound `v` (from interval_ins over the box `b`),
    /// within these: a derivative's, within its field's partial; one whose
    /// gradient bounds a derivative reads, within their value, which it
    /// then passes on with them.
    fn tighten(&mut self, tape: &Tape, ctx: &Ctx, i: usize, ins: &crate::tape::Ins, b: [I; 3], v: I) -> I {
        let mut v = v;
        if let Op::Deriv(k) = ins.op
            && k < 3
        {
            v = interval::meet(v, self.0[ins.b as usize].d[k as usize]);
        }
        if tape.slopes[i] {
            let g = grad_ins(tape, ctx, i, ins, &self.0, &self.1, b);
            self.1[ins.out as usize] = real_flags(ins, &self.0, &self.1, &g);
            v = interval::meet(v, g.v);
            self.0[ins.out as usize] = interval::G { v, ..g };
        }
        v
    }
}
