//! Evaluation inputs and the spans and derivations of ranged parameters.

use crate::graph::{Var, bin_f64, un_f64};
use crate::interval::I;
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

    /// What evaluating `tapes` reads of the context: the params they
    /// read (Tape::params_read), and the derivations of those that are
    /// derived, in order, with what each reads in turn (what it derives,
    /// reads and its equations read), sorted.
    pub fn read_by(&self, tapes: &[&Tape]) -> (Vec<u32>, Vec<std::sync::Arc<dyn Derive>>) {
        let mut read: Vec<u32> = tapes.iter().flat_map(|t| t.params_read()).collect();
        // Each derivation reads only those before it: one pass from the last.
        let mut derive = Vec::new();
        for d in self.p_spans.derive.0.iter().rev() {
            if d.owns().iter().any(|k| read.contains(k)) {
                read.extend(d.owns().iter().chain(d.reads()).copied().chain(d.tapes().iter().flat_map(|t| t.params_read())));
                derive.insert(0, d.clone());
            }
        }
        read.sort_unstable();
        read.dedup();
        (read, derive)
    }

    /// The context as evaluating `tapes` reads it (`read_by`): t and its
    /// span, the ranged params read and what derives them, and the values
    /// of the params read; every other param NaN and unranged, and none
    /// past the last kept. A cached answer is computed from this context
    /// and keyed by all of it (`fingerprint`), so it cannot read what its
    /// key leaves out; a witness found in it names only the ranged params
    /// read (the answer is the same at any value of the others).
    pub fn project(&self, tapes: &[&Tape]) -> Ctx {
        let (keep, derive) = self.read_by(tapes);
        let mut params = vec![f64::NAN; keep.last().map_or(0, |&k| k as usize + 1)];
        for &k in &keep {
            params[k as usize] = self.params.get(k as usize).copied().unwrap_or(f64::NAN);
        }
        let list: Vec<ParamSpan> = self.p_spans.list.iter().filter(|s| keep.binary_search(&s.k).is_ok()).copied().collect();
        let mask = list.iter().fold(0, |m, s| m | 1 << (s.k % 64));
        Ctx { params, p_spans: ParamSpans { mask, list, derive: Derivations(derive) }, ..self.clone() }
    }

    /// Hash all of the context a question reads into `fp` (of a projected
    /// one: `project`): t and its span where a tape reads t (an answer that
    /// names a moment gives the asker's: query::clearance), the
    /// params by their numbers in it
    /// (`names`); a ranged param by the model's number too, which a
    /// witness's `params` give; and what derives the derived ones, as data
    /// (Derive::fingerprint). False: a derivation is not data.
    pub fn fingerprint(&self, fp: &mut crate::graph::Fingerprint, names: &crate::tape::Names) -> bool {
        // Every field, named: one added must be hashed too.
        let Ctx { t, params, t_span, p_spans } = self;
        // `mask` is the ranged params' numbers' bits, hashed with them.
        let ParamSpans { mask: _, list, derive } = p_spans;
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
        fp.word(derive.0.len() as u64);
        derive.0.iter().all(|d| {
            fp.word(d.tapes().len() as u64);
            d.tapes().iter().for_each(|t| t.fingerprint(fp, names));
            d.fingerprint(fp, names)
        })
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

/// Ranged params derived from the others over a scope: a solve over
/// ranged params (no form reads one now: the solver left, to return as a
/// plug), whose unknowns are the root of its equations, certified to be one and enclosed for
/// every value of the params it reads. What it sets must hold for every
/// value of those params in their spans (outward rounding throughout):
/// proofs cover its spans, and a witness is confirmed over the span it
/// gives at one point of the others.
pub trait Derive: Send + Sync {
    /// The params it derives.
    fn owns(&self) -> &[u32];
    /// The equations it evaluates in a context: what it reads of one is
    /// their params, `reads`, its own, and t where they read it.
    fn tapes(&self) -> &[Tape] {
        &[]
    }
    /// Hash everything else it derives from into `fp` (params numbered by
    /// `names`), so that two derivations hashed alike derive alike from
    /// what `tapes` says they read: a cached answer over its params is
    /// keyed so (Ctx::fingerprint). False: it is not data (an answer
    /// reading its params is not cached).
    fn fingerprint(&self, _fp: &mut crate::graph::Fingerprint, _names: &crate::tape::Names) -> bool {
        false
    }
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
    /// Whether it bounds its params' rates (`rate_bounds`): a first-order
    /// bound then reads them through the params they are derived from
    /// (scope_run), not as inputs of their own.
    fn chains(&self) -> bool {
        false
    }
    /// For each of its params, an enclosure of its value at the middle of
    /// the spans of the params it reads in ctx, and bounds on its rate
    /// along each of them (`reads` order) over all of their spans: it lies
    /// in the first plus the rates times how far each read param is from
    /// its middle. None: not bounded there.
    fn rate_bounds(&self, _ctx: &Ctx) -> Option<Vec<(I, Vec<I>)>> {
        None
    }
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

    /// Independent param k's span and band narrowed to `span` (a quotient
    /// by a period: ink_core::period::quotient).
    pub fn narrow(&mut self, k: u32, span: (f64, f64)) {
        if let Some(s) = self.list.iter_mut().find(|s| s.k == k && !s.derived) {
            (s.span, s.band) = (span, span);
        }
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

    /// What derives param k, where it is derived, bounds its rates
    /// (Derive::chains) and reads only independent ranged params: a
    /// first-order bound reads k through them.
    pub fn chained(&self, k: u32) -> Option<&std::sync::Arc<dyn Derive>> {
        if !self.is_derived(k) {
            return None;
        }
        let d = self.derive.0.iter().find(|d| d.owns().contains(&k))?;
        (d.chains() && d.reads().iter().all(|&r| self.get(r).is_some() && !self.is_derived(r))).then_some(d)
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

    /// Every param ranged here or that a derivation owns or reads,
    /// sorted: those `chain` passes derivatives between.
    pub fn along(&self) -> Vec<u32> {
        let mut ks: Vec<u32> = self.list.iter().map(|s| s.k).collect();
        for d in &self.derive.0 {
            ks.extend(d.owns().iter().chain(d.reads()).copied());
        }
        ks.sort_unstable();
        ks.dedup();
        ks
    }

    /// The independent ranged params' values in `params` (a witness's: the
    /// derived ones follow from them).
    pub fn values(&self, params: &[f64]) -> Vec<(u32, f64)> {
        self.list.iter().filter(|s| !s.derived).map(|s| (s.k, params[s.k as usize])).collect()
    }
}

/// The middle of a span (of time, or of a ranged param), whose bounds
/// scope_run builds on: in the span and finite for any finite span (halved
/// before the sum, which cannot overflow: 0.5 (t0 + t1) is infinite for
/// t0 = 1e308, t1 = 1.1e308).
pub fn span_mid((t0, t1): (f64, f64)) -> f64 {
    I { lo: t0, hi: t1 }.mid()
}

/// ctx with input `var` at `at` alone: time no span, a ranged param a
/// span of one value (still ranged, so it keeps its meaning in real
/// numbers: Inputs::centre).
pub(super) fn input_at(ctx: &Ctx, var: Var, at: f64) -> Ctx {
    let mut c = ctx.clone();
    match var {
        Var::T => (c.t, c.t_span) = (at, None),
        Var::Param(k) => {
            c.params[k as usize] = at;
            c.p_spans.list.iter_mut().filter(|s| s.k == k).for_each(|s| s.span = (at, at));
        }
        _ => unreachable!("an input of the scope"),
    }
    c
}

/// Inputs beyond space a first-order bound takes at once (time and ranged
/// params). Any more ranged params stay ranged in its centre (see Inputs).
pub const MAX_INPUTS: usize = 8;

/// The inputs beyond space a first-order bound is taken in (scope_run): its
/// axes, time first. A ranged param that is not one (beyond MAX_INPUTS)
/// keeps its span in the bound's centre, so the bound still holds for every
/// value of it: the centre is then a box along it, and the bound is first
/// order only in the axes. A derived param whose rates are bounded
/// (ParamSpans::chained: a solve on t's unknowns over a piece of t) is no
/// input: it is read through those it is derived from, at its rates along
/// them, so values it closes against them (a loop the solve closes) keep
/// their dependency.
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
    pub(super) fn read<'a>(tapes: impl IntoIterator<Item = &'a Tape>, ctx: &Ctx) -> Inputs {
        let spans_t = ctx.t_span.is_some_and(|(a, b)| a < b);
        let (mut t, mut params) = (false, Vec::new());
        if spans_t || !ctx.p_spans.is_empty() {
            let mut stack: Vec<&Tape> = tapes.into_iter().collect();
            let mut seen = Vec::new();
            while let Some(tape) = stack.pop() {
                let ptr = tape as *const Tape;
                if seen.contains(&ptr) {
                    continue;
                }
                seen.push(ptr);
                for ins in &tape.ssa {
                    match ins.op {
                        Op::Var(Var::T) => t |= spans_t,
                        Op::Var(Var::Param(k)) if ctx.p_spans.get(k).is_some() => match ctx.p_spans.chained(k).filter(|d| rated(d.as_ref(), ctx).is_some()) {
                            Some(d) => params.extend(d.reads()),
                            None => params.push(k),
                        },
                        Op::CornerReach(fa, fb) => {
                            stack.push(&tape.forms[fa as usize]);
                            stack.push(&tape.forms[fb as usize]);
                        }
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
    pub(super) fn passes(&self) -> Vec<Inputs> {
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
    pub(super) fn axis(&self, v: Var) -> Option<usize> {
        match v {
            Var::T if self.t => Some(0),
            Var::Param(k) => self.params.iter().position(|&p| p == k).map(|j| j + self.t as usize),
            _ => None,
        }
    }

    /// ctx at their centre: t at the middle of its span (no span), each
    /// param at the middle of its span, still ranged (so it keeps its
    /// meaning in real numbers: Ctx::point).
    pub(super) fn centre(&self, ctx: &Ctx) -> Ctx {
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

/// `d`'s rate bounds over ctx's spans (Derive::rate_bounds), where every
/// rate is bounded: then its params are read through those it reads
/// (Inputs::read, `chained`); else each is an input of its own.
pub(super) fn rated(d: &dyn Derive, ctx: &Ctx) -> Option<Vec<(I, Vec<I>)>> {
    d.rate_bounds(ctx).filter(|b| b.iter().flat_map(|x| &x.1).all(|r| r.lo.is_finite() && r.hi.is_finite()))
}

/// Whether instruction i is a constants-and-params value taken in floats
/// (Tape::fixed): one that reads no ranged param. A ranged one means its
/// value in real numbers, bounded like any other. The params' bits
/// (Tape::pmask) rule most out at once; where they meet the ranged ones',
/// the params it reads (Tape::params_of) decide, so a param sharing only a
/// ranged one's bit (64 apart) leaves it in floats.
#[inline]
pub(crate) fn fixed(tape: &Tape, ctx: &Ctx, i: usize) -> bool {
    tape.fixed[i] && (tape.pmask[i] & ctx.p_spans.mask == 0 || tape.params_of(i).iter().all(|&k| ctx.p_spans.get(k).is_none()))
}

/// A constants-and-params instruction in floats, from its argument values.
/// Params are fixed for a query, so this is the value a def would have: the
/// interval bound of such an expression is this point, not an enclosure that
/// widens through long chains (and through atan2 or floor near a cut).
pub(crate) fn fixed_value(ins: &crate::tape::Ins, arg: impl Fn(u32) -> f64, ctx: &Ctx) -> f64 {
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
