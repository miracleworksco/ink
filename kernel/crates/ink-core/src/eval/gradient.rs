//! Spatial gradient bounds and tightening of derivative ranges.

use super::context::{Ctx, fixed, fixed_value};
use super::corner::corner_reach_interval;
use super::fork::{Fork, Gap, PAIRED, forked_live, same_branches};
use super::promise::promise_valid;
use super::range::{interval_ins, remainder_within};
use crate::graph::{BinOp, UnOp, Var};
use crate::interval::{self, I};
use crate::tape::{Op, Tape};

/// Value and gradient bounds over a box (see interval::G). The gradient bound
/// is a Lipschitz constant for the field there: |f(p) - f(q)| <= L |p - q|.
pub fn eval_interval_grad(tape: &Tape, ctx: &Ctx, b: [I; 3]) -> interval::G {
    grad_run(tape, &tape.code, tape.nslots, false, ctx, b)[tape.out_slot()]
}

/// `eval_interval_grad`'s bounds of every value of the tape's SSA form (a
/// max's terms: the clearance proofs' first-order forms, ink-kernel's
/// `centred`).
pub fn eval_interval_grad_every(tape: &Tape, ctx: &Ctx, b: [I; 3]) -> Vec<interval::G> {
    grad_run(tape, &tape.ssa, tape.ssa.len(), true, ctx, b)
}

/// Interval forward-mode AD over `code` (the tape's slot code, or its SSA
/// form to keep every value, `every`): the bounds in every slot.
///
/// A jump the box straddles once is a choice of two sides (`forked`), as in
/// eval_interval, each its branch continued over the box (a floor's one
/// integer, slope 0; atan2's branch near pi or near -pi, its slope atan2's
/// formula); a floor straddling several integers too, each side a run of
/// them (slope 0 between its crossings), its seam checked at all its
/// crossings at once (a belt's teeth, 2 mm apart, in a box of 3). Each side's bounds hold for the tape with the jump taken so,
/// over the whole box, so for the field wherever it takes that side. Where
/// the field is continuous across the seam, its slope on any segment is
/// one side's or the other's, so the sides' joined bounds hold: where the
/// two sides are one function over the box (their difference exact and 0:
/// atan2's branches a whole turn apart, read by a sine or cosine, or in
/// degrees as turns less their floor; Grad::gap), or give one exact value
/// where they meet (Grad::seam: a wrap read by its distance from 1/2, a
/// gear's angle from its nearest tooth's centre line).
/// Elsewhere the field may jump where the sides meet, so its partials are
/// unbounded, but along an axis the jump's arguments do not move along
/// (each line along it takes one side throughout).
fn grad_run(tape: &Tape, code: &[crate::tape::Ins], nslots: usize, every: bool, ctx: &Ctx, b: [I; 3]) -> Vec<interval::G> {
    let mut run = Grad { tape, code, ctx, b, s: vec![interval::G::constant(0.0); nslots], every, g: Vec::new(), pairs: Vec::new() };
    forked_live(tape, &mut run, None);
    run.s
}

/// A fork grad_run is inside: its jump `i`, where it joins `r`, whether its
/// second side is running, the axes its jump's arguments are constant
/// along over the box, and which values of its window i..=r its sides give
/// alike where they meet (Grad::seam).
struct GradPair {
    i: usize,
    r: usize,
    second: bool,
    flat: [bool; 3],
    seam: Vec<bool>,
}

/// grad_run's state.
struct Grad<'a> {
    tape: &'a Tape,
    code: &'a [crate::tape::Ins],
    ctx: &'a Ctx,
    b: [I; 3],
    s: Vec<interval::G>,
    every: bool,
    /// Per slot, for each fork in `pairs` (the innermost PAIRED) whose
    /// second side is running: the value there less the first side's at
    /// each point, where exact (Gap::exact; None: the same; a Gap of
    /// ENTIRE: unknown). Empty until a second side runs.
    g: Vec<[Option<Gap>; PAIRED]>,
    pairs: Vec<GradPair>,
}

impl Grad<'_> {
    /// Instruction `ins`'s difference between the sides of fork `k`, from
    /// its arguments': exact (Gap::exact), through sums, and through
    /// products and quotients by a value the same on both sides; else
    /// unknown. None where it is 0.
    fn gap(&self, ins: &crate::tape::Ins, k: usize) -> Option<Gap> {
        let (args, na) = ins.args();
        let g = |s: u32| self.g[s as usize][k];
        if args[..na].iter().all(|&s| g(s).is_none()) {
            return None;
        }
        let zero = I::point(0.0);
        let gap = |s: u32| g(s).unwrap_or(Gap::of(zero));
        let v = |s: u32| self.s[s as usize].v;
        let unknown = Gap::of(I::ENTIRE);
        let x = Gap::exact(ins, g).unwrap_or(match ins.op {
            Op::Bin(op @ (BinOp::Add | BinOp::Sub)) => Gap::sum(op, gap(ins.a), gap(ins.b)),
            Op::Un(UnOp::Neg) => Gap::sum(BinOp::Sub, Gap::of(zero), gap(ins.a)),
            Op::Bin(BinOp::Mul) if g(ins.b).is_none() && gap(ins.a).turns == 0 => Gap::of(interval::mul_o(gap(ins.a).v, v(ins.b))),
            Op::Bin(BinOp::Mul) if g(ins.a).is_none() && gap(ins.b).turns == 0 => Gap::of(interval::mul_o(v(ins.a), gap(ins.b).v)),
            Op::Bin(BinOp::Div) if g(ins.b).is_none() && gap(ins.a).turns == 0 => Gap::of(interval::bin(BinOp::Div, gap(ins.a).v, v(ins.b)).0),
            // floor(u + n) = floor(u) + n for a whole n.
            Op::Un(UnOp::Floor) => Grad::whole(gap(ins.a)).unwrap_or(unknown),
            // The value of `a` (a promise's constants bound slopes, not it).
            Op::Tag(_) | Op::Read | Op::Real | Op::Lipschitz if args[1..na].iter().all(|&s| g(s).is_none()) => gap(ins.a),
            _ => unknown,
        });
        (x != Gap::of(zero)).then_some(x)
    }

    /// `x` where it is a whole number exactly.
    fn whole(x: Gap) -> Option<Gap> {
        (x.turns == 0 && x.v.lo == x.v.hi && x.v.lo == x.v.lo.trunc()).then_some(x)
    }

    /// Which values of a floor `i`'s window i..=r its two sides give alike
    /// where they meet: the window bounded where the floor's argument is a
    /// whole number m it crosses, the floor m - 1 on one side and m on
    /// the other, every other value as over the box; m any of `m`, the
    /// whole numbers it crosses (one, or a run where it crosses several:
    /// a value is alike at all of them by one bound of them all). A value is alike where
    /// its arguments are, or where both sides bound it by one number, by an
    /// operation whose bound is a point only where it is exact (a sum,
    /// product, quotient, min, max, negation, magnitude or floor): there
    /// it is that number on both sides at every point where they meet.
    /// Where the floor's argument is a / p, p one nonzero number over the
    /// box, a is m p there, so a less p times the floor (std's mod) is p on
    /// the side below and 0 on the side above, exactly. A value is
    /// opposite where the sides bound it so by numbers v and -v, and odd
    /// and even functions carry it (Seam::of): a polar repeat's angle in
    /// its wedge is p/2 one side and -p/2 the other, its turned y
    /// opposite, x alike, so a copy built symmetric about its wedge's
    /// middle (of |y|, y y, ...) is alike.
    fn seam(&self, i: usize, r: usize, m: I) -> Vec<bool> {
        let ins = &self.code[i];
        let ssa = &self.tape.ssa;
        let quotient = &ssa[ssa[i].a as usize];
        let (a, p) = if matches!(quotient.op, Op::Bin(BinOp::Div)) { (Some(quotient.a), Some(quotient.b)) } else { (None, None) };
        // Per window value: p, where it is p times the floor.
        let mut times = vec![None; r - i + 1];
        let mut s0: Vec<I> = self.s.iter().map(|g| g.v).collect();
        s0[ins.a as usize] = m;
        let mut s1 = s0.clone();
        let o = ins.out as usize;
        // Exact: whole numbers under 2^52 (jump_sides).
        (s0[o], s1[o]) = (I::new(m.lo - 1.0, m.hi - 1.0), m);
        let mut sides = vec![Seam::Alike; s0.len()];
        sides[o] = Seam::Differ;
        let mut alike = vec![false; r - i + 1];
        for j in i + 1..=r {
            let x = &self.code[j];
            let (mut v0, mut v1) =
                (interval_ins(self.tape, self.ctx, j, self.code, &s0, self.b, None, None, false), interval_ins(self.tape, self.ctx, j, self.code, &s1, self.b, None, None, false));
            let y = &ssa[j];
            match y.op {
                Op::Bin(BinOp::Mul) if Some(y.a) == p && y.b as usize == i || Some(y.b) == p && y.a as usize == i => {
                    let v = s0[if Some(y.a) == p { x.a } else { x.b } as usize];
                    times[j - i] = (v.lo == v.hi && v.lo != 0.0).then_some(v);
                }
                Op::Bin(BinOp::Sub) if Some(y.a) == a && (y.b as usize) > i && (y.b as usize) < j => {
                    if let Some(v) = times[y.b as usize - i] {
                        (v0, v1) = (v, I::point(0.0));
                    }
                }
                _ => {}
            }
            let exact = matches!(x.op, Op::Un(UnOp::Neg | UnOp::Abs | UnOp::Floor) | Op::Bin(BinOp::Add | BinOp::Sub | BinOp::Mul | BinOp::Div | BinOp::Min | BinOp::Max));
            let points = exact && v0.lo == v0.hi && v1.lo == v1.hi;
            let side = match Seam::of(self.tape, self.code, j, &sides) {
                Seam::Alike => Seam::Alike,
                _ if points && v0 == v1 => Seam::Alike,
                _ if points && v0.lo == -v1.lo => Seam::Opposite,
                side => side,
            };
            let o = x.out as usize;
            (s0[o], s1[o], sides[o]) = (v0, v1, side);
            alike[j - i] = side == Seam::Alike;
        }
        alike
    }
}

/// How a value of a floor's window compares on the floor's two sides
/// where they meet (Grad::seam), at every such point: the same, one the
/// other's negation, or neither shown.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Seam {
    Alike,
    Opposite,
    Differ,
}

impl Seam {
    /// Instruction `x`'s, from its arguments' in `sides`: alike where they
    /// are all alike; else through odd functions (negation, sine, a
    /// product or quotient by an alike value: x / 0 is 0 either way) and
    /// even ones (magnitude, square, cosine, a product of two opposites),
    /// each an identity in real numbers, and a promise (its field's value).
    pub(super) fn of(tape: &Tape, code: &[crate::tape::Ins], j: usize, sides: &[Seam]) -> Seam {
        use Seam::*;
        let x = &code[j];
        if matches!(x.op, Op::CornerReach(..)) {
            let f = &code[tape.ssa[j].b as usize];
            return if [x.a, x.c, f.a, f.b, f.c].into_iter().all(|k| sides[k as usize] == Alike) { Alike } else { Differ };
        }
        let (args, na) = x.args();
        let side = |k: usize| sides[args[k] as usize];
        if (0..na).all(|k| side(k) == Alike) {
            return Alike;
        }
        match x.op {
            Op::Un(UnOp::Neg | UnOp::Sin | UnOp::Rad | UnOp::Deg) | Op::Tag(_) | Op::Real => side(0),
            Op::Un(UnOp::Abs | UnOp::Square | UnOp::Cos) if side(0) == Opposite => Alike,
            Op::Bin(BinOp::Mul | BinOp::Div) if side(0) != Differ && side(1) != Differ => {
                if side(0) == side(1) {
                    Alike
                } else {
                    Opposite
                }
            }
            Op::Bin(BinOp::Add | BinOp::Sub) if side(0) == Opposite && side(1) == Opposite => Opposite,
            // A promise's value is its field's; its frame's is 0.
            Op::Lipschitz => side(0),
            Op::Frame => Alike,
            _ => Differ,
        }
    }
}

impl Fork for Grad<'_> {
    type Out = (interval::G, [Option<Gap>; PAIRED]);
    type Args = (interval::G, interval::G);

    fn step(&mut self, i: usize, _: Option<&mut [u8]>) -> I {
        let ins = &self.code[i];
        let o = ins.out as usize;
        let x = grad_ins(self.tape, self.ctx, i, self.code, &self.s, self.b, !self.every);
        // Its differences from the forks' first sides, from its arguments'
        // before the output's slot (maybe one's) is written.
        if !self.g.is_empty() {
            let mut d = self.g[o];
            for (k, p) in self.pairs.iter().enumerate().take(PAIRED) {
                if p.second {
                    d[k] = self.gap(ins, k);
                }
            }
            self.g[o] = d;
        }
        self.s[o] = x;
        x.v
    }
    fn args(&self, i: usize) -> Self::Args {
        let ins = &self.code[i];
        (self.s[ins.a as usize], self.s[ins.b as usize])
    }
    fn bounds(args: &Self::Args) -> (I, I) {
        (args.0.v, args.1.v)
    }
    fn ctx(&self) -> &Ctx {
        self.ctx
    }
    fn code(&self) -> &[crate::tape::Ins] {
        self.code
    }
    fn value(&self, k: u32) -> I {
        self.s[k as usize].v
    }
    fn side(&mut self, i: usize, (y, x): &Self::Args, k: usize, side: I) {
        use interval::G;
        let ins = &self.code[i];
        let o = ins.out as usize;
        let zero = I::point(0.0);
        let (s, flat) = match ins.op {
            Op::Un(UnOp::Floor) => (G { v: side, d: [zero; 3], l: 0.0, s2: None, m: interval::NO_GRAM }, y.d.map(|d| d == zero)),
            Op::Bin(BinOp::Atan2) => {
                let d = interval::d_atan2(y.v, x.v, y.d, x.d).unwrap_or([I::ENTIRE; 3]);
                (G { v: side, d, l: f64::INFINITY, s2: None, m: interval::NO_GRAM }, std::array::from_fn(|j| y.d[j] == zero && x.d[j] == zero))
            }
            _ => unreachable!("jump_sides takes floors and atan2"),
        };
        self.s[o] = s;
        // The outer forks' differences: none where its arguments have none.
        // The outer forks' differences: a floor's argument's where whole
        // (each side's integer moves with it), else none where its
        // arguments have none.
        let outer = if k == 0 { self.pairs.len() } else { self.pairs.len() - 1 };
        if !self.g.is_empty() {
            for q in 0..outer.min(PAIRED) {
                if self.pairs[q].second {
                    let (ga, gb) = (self.g[ins.a as usize][q], self.g[ins.b as usize][q]);
                    self.g[o][q] = match ins.op {
                        Op::Un(UnOp::Floor) => ga.map(|x| Grad::whole(x).unwrap_or(Gap::of(I::ENTIRE))),
                        _ => (ga.is_some() || gb.is_some()).then_some(Gap::of(I::ENTIRE)),
                    };
                }
            }
        }
        if k == 0 {
            let r = self.tape.rejoin(i);
            let seam = match ins.op {
                Op::Un(UnOp::Floor) => self.seam(i, r, I::new(side.lo + 1.0, side.hi + 1.0)),
                _ => Vec::new(),
            };
            self.pairs.push(GradPair { i, r, second: false, flat, seam });
        } else {
            self.pairs[outer].second = true;
            if outer < PAIRED {
                if self.g.is_empty() {
                    self.g = vec![[None; PAIRED]; self.s.len()];
                }
                self.g.iter_mut().for_each(|g| g[outer] = None);
                self.g[o][outer] = Some(match ins.op {
                    // A run of integers less another: not one number.
                    Op::Un(UnOp::Floor) if side.lo != side.hi => Gap::of(I::ENTIRE),
                    Op::Un(UnOp::Floor) => Gap::of(I::point(1.0)),
                    _ => Gap { v: zero, turns: 1 },
                });
            }
        }
    }
    fn out(&self, i: usize) -> Self::Out {
        let o = self.code[i].out as usize;
        (self.s[o], self.g.get(o).copied().unwrap_or([None; PAIRED]))
    }
    fn put(&mut self, i: usize, (s, g): Self::Out) {
        let o = self.code[i].out as usize;
        self.s[o] = s;
        if let Some(x) = self.g.get_mut(o) {
            *x = g;
        }
    }
    /// The fork joined is the innermost.
    fn join(&mut self, i: usize, a: Self::Out, b: Self::Out) {
        let o = self.code[i].out as usize;
        let depth = self.pairs.len() - 1;
        let p = &self.pairs[depth];
        let same = (depth < PAIRED && b.1[depth].is_none()) || p.seam.get(i - p.i) == Some(&true);
        self.s[o] = if same {
            a.0.hull(&b.0)
        } else {
            let d = std::array::from_fn(|k| if p.flat[k] { interval::hull(a.0.d[k], b.0.d[k]) } else { I::ENTIRE });
            interval::G { v: interval::hull(a.0.v, b.0.v), d, l: f64::INFINITY, s2: None, m: interval::NO_GRAM }
        };
        // The outer forks' differences: either side's.
        if let Some(g) = self.g.get_mut(o) {
            *g = std::array::from_fn(|k| match (a.1[k], b.1[k]) {
                _ if k >= depth => None,
                (x, y) if x == y => x,
                _ => Some(Gap::of(I::ENTIRE)),
            });
        }
        if i == p.r {
            self.pairs.pop();
        }
    }
    fn every_value(&self) -> bool {
        self.every
    }
    fn many(&self) -> bool {
        true
    }
}

/// Instruction `i`'s value and gradient bounds (interval::G) from its
/// arguments' in `s` (and their flags, `real_flags`, in `f`).
fn grad_ins(tape: &Tape, ctx: &Ctx, i: usize, code: &[crate::tape::Ins], s: &[interval::G], b: [I; 3], slots: bool) -> interval::G {
    use interval::G;
    let ins = &code[i];
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
        Op::Bin(op) => {
            let g = interval::g_bin(op, s[ins.a as usize], s[ins.b as usize]);
            if slots && op == BinOp::Sub { G { v: remainder_within(tape, i, |k| s[k].v, g.v), ..g } } else { g }
        }
        Op::Select => interval::g_select(s[ins.a as usize], s[ins.b as usize], s[ins.c as usize], same_branches(tape, i)),
        Op::Tag(_) | Op::Read | Op::Real => s[ins.a as usize],
        // Constant in space, somewhere in its range (as a ranged param).
        Op::Within(_) => {
            let v = I::new(s[ins.a as usize].v.lo, s[ins.b as usize].v.hi);
            G { v, d: [I::point(0.0); 3], l: 0.0, s2: (v.lo >= 0.0).then_some(0.0), m: interval::NO_GRAM }
        }
        // Trusted only where its field's, k's and frame's constants are real.
        Op::Lipschitz if !promise_valid(tape, ctx, i, b) => s[ins.a as usize],
        Op::Lipschitz => interval::g_lipschitz(s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]),
        Op::Frame => interval::g_frame([s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]]),
        // Sampled fields: value bounds only (no derivative bounds yet).
        Op::Grid(gi, m) => interval::g_grid(&tape.grids[gi as usize], m, [s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]]),
        Op::CornerReach(fa, fb) => {
            let f = &code[tape.ssa[i].b as usize];
            let q = [s[f.a as usize].v, s[f.b as usize].v, s[f.c as usize].v];
            let v = corner_reach_interval(&tape.forms[fa as usize], &tape.forms[fb as usize], ctx, s[ins.a as usize].v, q, s[ins.c as usize].v);
            G { v, d: [I::ENTIRE; 3], l: f64::INFINITY, s2: None, m: interval::NO_GRAM }
        }
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
pub(super) struct Slopes(pub(super) Vec<interval::G>);

impl Slopes {
    /// For an evaluator over `n` slots (the tape's code, or its SSA form).
    pub(super) fn of(tape: &Tape, n: usize) -> Option<Slopes> {
        (!tape.slopes.is_empty()).then(|| Slopes(vec![interval::G::constant(0.0); n]))
    }

    /// Instruction `i`'s bound `v` (from interval_ins over the box `b`),
    /// within these: a derivative's, within its field's partial; one whose
    /// gradient bounds a derivative reads, within their value, which it
    /// then passes on with them. `slots` as interval_ins's.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn tighten(&mut self, tape: &Tape, ctx: &Ctx, i: usize, code: &[crate::tape::Ins], b: [I; 3], v: I, slots: bool) -> I {
        let ins = &code[i];
        let mut v = v;
        if let Op::Deriv(k) = ins.op
            && k < 3
        {
            v = interval::meet(v, self.0[ins.b as usize].d[k as usize]);
        }
        if tape.slopes[i] {
            let g = grad_ins(tape, ctx, i, code, &self.0, b, slots);
            v = interval::meet(v, g.v);
            self.0[ins.out as usize] = interval::G { v, ..g };
        }
        v
    }
}
