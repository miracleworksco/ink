//! First-order bounds over time and ranged parameters, with rates and jumps.

use super::cases::moment;
use super::context::{Ctx, Inputs, MAX_INPUTS, fixed, input_at, rated};
use super::fork::{Fork, Gap, PAIRED, forked, forked_live, same_branches};
use super::gradient::Slopes;
use super::join;
use super::promise::Frames;
use super::range::{eval_interval, interval_ins};
use crate::graph::{BinOp, UnOp, Var};
use crate::interval::{self, I, sub_o};
use crate::tape::{Op, Tape};

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
    scope_run(tape, false, ctx, &inputs, b, choices, false, true).v[tape.out_slot()]
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
        let sc = scope_run(tape, false, ctx, &inputs, b, choices.take(), true, true);
        let (o, n, hw) = (tape.out_slot(), inputs.len(), inputs.half_widths(ctx));
        v = interval::meet(v, sc.v[o]);
        let axes = inputs.t.then_some(Var::T).into_iter().chain(inputs.params.iter().map(|&k| Var::Param(k)));
        // A jump is every input's share (any may move its seam), as were
        // unbounded rates.
        blur.extend(axes.enumerate().map(|(j, x)| {
            let r = interval::mag(sc.d[o * n + j]);
            (x, interval::add_up(if r == 0.0 { 0.0 } else { interval::mul_up(r, hw[j]) }, sc.j[o]))
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
    let (v, d, j) = eval_interval_rates_jump(tape, ctx, b, inputs);
    (v, if j > 0.0 { vec![I::ENTIRE; d.len()] } else { d })
}

/// Prove a scalar expression continuous on a closed time span by comparing
/// exact values at its switches. This can settle expressions that the
/// interval rate bound leaves undecided, such as a switching factor times a
/// factor that is zero right at the switch. False includes unproved cases.
pub fn continuous_in_time(tape: &Tape, ctx: &Ctx, b: [I; 3], span: (f64, f64)) -> bool {
    join::continuous(tape, ctx, b, span).is_some()
}

/// eval_interval_rates, with where the value may jump along the inputs
/// (a seam they move, its sides joined: Scoped) each side's rates and a
/// bound on the jump's height: it moves from its value at the inputs'
/// centre by at most the rates times how far each input moves plus that.
pub fn eval_interval_rates_jump(tape: &Tape, ctx: &Ctx, b: [I; 3], inputs: &Inputs) -> (I, Vec<I>, f64) {
    assert!(inputs.len() <= MAX_INPUTS, "at most MAX_INPUTS rates at once");
    if inputs.is_empty() {
        return (eval_interval(tape, ctx, b, None), Vec::new(), 0.0);
    }
    let sc = scope_run(tape, false, ctx, inputs, b, None, true, true);
    let (o, n) = (tape.out_slot(), inputs.len());
    (sc.v[o], sc.d[o * n..(o + 1) * n].to_vec(), sc.j[o])
}

/// eval_interval_rates in one run of the tape, every jump bounded across
/// (not forked: `forked`), so a value the box and scope straddle a jump
/// of has unbounded rates. Forked, a jump whose seam the inputs move
/// jumps there by its sides' difference, so its rates are unbounded too
/// unless that is exactly 0; one whose seam they do not move has its
/// branches' rates, which this does not bound. So the rates are
/// eval_interval_rates', or unbounded where only a fork bounds them. A
/// fork bounds its window twice, so forks nested in each other cost a
/// doubling apiece; one run costs what the tape is long: for a question
/// one run must answer (whether a value moves at all: ink-kernel's
/// motion::still).
pub fn eval_interval_rates_across(tape: &Tape, ctx: &Ctx, b: [I; 3], inputs: &Inputs) -> (I, Vec<I>) {
    assert!(inputs.len() <= MAX_INPUTS, "at most MAX_INPUTS rates at once");
    let sc = scope_run(tape, false, ctx, inputs, b, None, true, false);
    let (o, n) = (tape.out_slot(), inputs.len());
    (sc.v[o], if sc.j[o] > 0.0 { vec![I::ENTIRE; n] } else { sc.d[o * n..(o + 1) * n].to_vec() })
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
/// of the sides'. Where the inputs move the seam it takes one side at some
/// values of them and the other at others, so it may also jump along them,
/// by at most the second side less the first at a point: that difference
/// is carried through the window like a rate along one more input (each
/// value's from its arguments' by its slope over the way between the
/// sides' values; sine and cosine first take whole turns off it), met with
/// the two sides' bounds' difference, and joined into the value's jump
/// (Scoped::j), which goes on through what reads it by its slope as rates
/// do. A value moves from the centre by at most its rates times how far
/// each input moves plus its jump. A turning gear's teeth over a step are
/// bounded as one tooth near its seam, turning; a polar repeat turned
/// across atan2's cut is one copy on both sides (a whole turn apart), so
/// it jumps by rounding only. Callers that need a value continuous (a
/// speed, a followed root, a Jacobian) take a jump as unbounded rates
/// (eval_interval_rates, Scoped::rate).
///
/// `choices`, as in eval_interval, records the branches decided over the
/// whole box and scope.
#[allow(clippy::too_many_arguments)]
fn scope_run(tape: &Tape, every: bool, ctx: &Ctx, inputs: &Inputs, b: [I; 3], choices: Option<&mut [u8]>, rates: bool, fork: bool) -> Scoped {
    match inputs.len() {
        0 | 1 => scope_run_n::<1>(tape, every, ctx, inputs, b, choices, rates, fork),
        2 => scope_run_n::<2>(tape, every, ctx, inputs, b, choices, rates, fork),
        3 | 4 => scope_run_n::<4>(tape, every, ctx, inputs, b, choices, rates, fork),
        _ => scope_run_n::<MAX_INPUTS>(tape, every, ctx, inputs, b, choices, rates, fork),
    }
}

/// What scope_run bounds, per slot: the value over the box and scope (`v`),
/// at the inputs' centre (`m`), its rates along them (`d`, `inputs.len()`
/// per slot, when asked for: 0 where the value does not depend on them,
/// ENTIRE where unbounded), and how far it may jump along them (`j`: 0
/// where it is continuous in them). The value moves from the centre by at
/// most its rates times how far each input moves, plus its jump.
pub struct Scoped {
    pub v: Vec<I>,
    pub m: Vec<I>,
    pub d: Vec<I>,
    pub j: Vec<f64>,
}

impl Scoped {
    /// Slot `o`'s rate along input `k` of `n`, as a bound on how fast the
    /// value moves: ENTIRE where it may jump (its rates bound each side's
    /// slope, not a jump between them).
    pub fn rate(&self, o: usize, n: usize, k: usize) -> I {
        if self.j[o] > 0.0 { I::ENTIRE } else { self.d[o * n + k] }
    }
}

/// scope_run over every value of the tape's SSA form (a max's terms: the
/// clearance proofs' first-order forms, ink-kernel's `centred`), with the
/// rates along `inputs` where `rates`.
pub fn eval_interval_scope_every(tape: &Tape, ctx: &Ctx, inputs: &Inputs, b: [I; 3], rates: bool) -> Scoped {
    scope_run(tape, true, ctx, inputs, b, None, rates, true)
}

/// The derived params `tape` reads through the inputs they are derived
/// from (ParamSpans::chained, Derive::rate_bounds): each one's enclosure at
/// the inputs' centre and its rates along them. One whose rates are not
/// all bounded is an input of its own (Inputs::read); one that reads a
/// param not among the inputs is read as a value ranged over its span, as
/// a ranged param that is no input is; one the inputs name stays one.
fn chained<const N: usize>(tape: &Tape, ctx: &Ctx, inputs: &Inputs) -> Vec<(u32, I, [I; N])> {
    let mut out = Vec::new();
    let mut done: Vec<*const ()> = Vec::new();
    for ins in &tape.ssa {
        let Op::Var(Var::Param(k)) = ins.op else { continue };
        // One asked for as an input (rates along it wanted) is one.
        let Some(d) = ctx.p_spans.chained(k).filter(|d| d.owns().iter().all(|&u| inputs.axis(Var::Param(u)).is_none())) else { continue };
        let id = std::sync::Arc::as_ptr(d) as *const ();
        if done.contains(&id) {
            continue;
        }
        done.push(id);
        let Some(axes) = d.reads().iter().map(|&r| inputs.axis(Var::Param(r))).collect::<Option<Vec<usize>>>() else { continue };
        let Some(bounds) = rated(d.as_ref(), ctx) else { continue };
        for (&u, (m, rates)) in d.owns().iter().zip(bounds) {
            let mut du = [I::point(0.0); N];
            for (&a, r) in axes.iter().zip(rates) {
                du[a] = r;
            }
            out.push((u, m, du));
        }
    }
    out
}

/// scope_run with N rates per value (N at least inputs.len(); the rest are
/// 0), its jumps forked where `fork` (forked_live), else every one bounded
/// across.
#[allow(clippy::too_many_arguments)]
fn scope_run_n<const N: usize>(tape: &Tape, every: bool, ctx: &Ctx, inputs: &Inputs, b: [I; 3], choices: Option<&mut [u8]>, rates: bool, fork: bool) -> Scoped {
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
        j: vec![0.0; n],
        g: Vec::new(),
        pairs: Vec::new(),
        slopes: Slopes::of(tape, n),
        every,
        fv: Frames::default(),
        next_gap: None,
        fm: Frames::default(),
        chain: chained(tape, ctx, inputs),
    };
    if fork {
        forked_live(tape, &mut run, Some(ch));
    } else {
        forked(tape, &mut run, 0, code.len(), Some(ch), &vec![false; code.len()]);
    }
    let k = inputs.len();
    let mut flat = Vec::with_capacity(if rates { n * k } else { 0 });
    if rates {
        for x in &run.d {
            flat.extend_from_slice(&x.unwrap_or([zero; N])[..k]);
        }
    }
    Scoped { v: run.v, m: run.m, d: flat, j: run.j }
}

/// A fork scope_run is inside: its jump `i`, where it joins `r`, whether
/// the inputs move its seam, whether its second side is running, and its
/// first side's bound of each value of its window i..=r (over each side of
/// any fork inside it).
struct Pair {
    i: usize,
    r: usize,
    moving: bool,
    second: bool,
    first: Vec<Option<I>>,
    /// Whether a derivative along t is marked of a value that differs
    /// between the sides: it has the jump's impulse.
    impulse: bool,
    /// The jump's second side less its first where it switches: a floor's
    /// next integer, atan2's whole turn, a select's branches' difference at
    /// the Moment it switches.
    gap: Gap,
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
    /// Per slot: how far the value may jump along the inputs, where they
    /// move a seam it was bounded on each side of (0: continuous in them).
    /// Its rates bound each side's slope; it moves from its value at the
    /// centre by at most them times how far each input moves, plus this.
    j: Vec<f64>,
    /// Per slot, for each fork in `pairs` whose second side is running:
    /// the value there less its value on the first side at each point of
    /// the region (None: the same). Empty until a second side runs.
    g: Vec<[Option<Gap>; PAIRED]>,
    pairs: Vec<Pair>,
    /// Gradient bounds hold over the whole scope (each partial at every
    /// moment and param value), so they tighten the bound over it.
    slopes: Option<Slopes>,
    every: bool,
    /// The frames' extents over the scope and at the centre (interval_ins).
    fv: Frames,
    /// The select forked next: its branches' difference at its Moment
    /// (select_sides), for its Pair.
    next_gap: Option<Gap>,
    fm: Frames,
    /// Derived params read through the inputs they are derived from
    /// (ParamSpans::chained): each one's enclosure at the inputs' centre
    /// and its rates along them.
    chain: Vec<(u32, I, [I; N])>,
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

/// `x` less the whole turns (2 pi each) nearest its middle: what moves the
/// sine and cosine of a value moved by `x` as far.
fn off_turns(x: I) -> I {
    use std::f64::consts::TAU;
    let k = (x.mid() / TAU).round();
    if k == 0.0 || !k.is_finite() {
        return x;
    }
    interval::sub_o(x, interval::mul_o(I::point(k), I::new(TAU, TAU.next_up())))
}

/// What a value's slot holds in scope_run: its bounds over the scope and at
/// the centre, rates, how far it may jump, its gradient bounds, and its
/// differences from the first sides of the forks paired.
type ScopeSlot<const N: usize> = (I, I, Option<[I; N]>, f64, Option<interval::G>, [Option<Gap>; PAIRED]);

impl<const N: usize> Scope<'_, N> {
    /// The derivative rules for instruction `i` (its bound `over`, its
    /// branch `choice`), its arguments' partials `d` by slot along M
    /// inputs, time the `t`-th: its partials, None where unbounded.
    fn rule<const M: usize>(&self, i: usize, over: I, choice: u8, d: impl Fn(usize) -> [I; M], t: Option<usize>) -> Option<[I; M]> {
        let (tape, ins, v) = (self.tape, &self.code[i], &self.v);
        let (a, bb, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
        match ins.op {
            Op::Un(op) => interval::d_un(op, v[a], over, d(a)),
            Op::Bin(op) => interval::d_bin(op, v[a], v[bb], over, choice, d(a), d(bb)),
            Op::Select => interval::d_select(choice, same_branches(tape, i), d(bb), d(c)),
            // A mark of `a` as the derivative along t of `b`: where b may
            // jump (its rate along t has no bound, or a seam t moves joins
            // its sides), so may a, by the jump's impulse, whatever a's
            // formula.
            Op::Deriv(3) => {
                let mut r = d(a);
                if let Some(j) = t
                    && !(d(bb)[j].lo.is_finite() && d(bb)[j].hi.is_finite() && self.j[bb] == 0.0)
                {
                    r[j] = I::ENTIRE;
                }
                Some(r)
            }
            // The value of `a` (a promise bounds slopes in space, not
            // rates in time or params).
            Op::Tag(_) | Op::Lipschitz | Op::Deriv(_) | Op::Read | Op::Real => Some(d(a)),
            Op::Grid(gi, mk) => Some(interval::d_grid(&tape.grids[gi as usize], mk, [v[a], v[bb], v[c]], [d(a), d(bb), d(c)])),
            Op::Frame | Op::Within(_) => Some([I::point(0.0); M]),
            Op::CornerReach(..) => None,
            Op::Const(_) | Op::Var(_) => unreachable!(),
        }
    }

    /// Instruction `i`'s value `v` (bound over the region on the sides
    /// running) less its value at each point on the first side of fork
    /// `k` (whose second is running): None where no argument differs.
    /// Else by its slope over each argument's values on the way from the
    /// first side's to the second's, as a rate along an input that takes
    /// one to the other (sine and cosine taking whole turns off the
    /// difference first: equal there); met with v less the first side's
    /// bound of it.
    fn diff(&self, i: usize, k: usize, v: I, vals: [I; 3]) -> Option<Gap> {
        let ins = &self.code[i];
        let (args, na) = ins.args();
        let g = |s: u32| self.g[s as usize][k];
        if args[..na].iter().all(|&s| g(s).is_none()) {
            return None;
        }
        let zero = I::point(0.0);
        let gs = |s: u32| g(s).map_or(zero, Gap::bound);
        if let Some(x) = Gap::exact(ins, g) {
            return Some(x);
        }
        // An argument's values on the way: the second side's less [0, 1]
        // times the difference.
        let (a, b, c) = (ins.a, ins.b, ins.c);
        let way = |s: u32, x: I| {
            let w = if s == a {
                vals[0]
            } else if s == b {
                vals[1]
            } else {
                vals[2]
            };
            interval::sub_o(w, interval::hull(zero, x))
        };
        let slope = match ins.op {
            // floor(u + x) - floor(u) lies in [floor(x), ceil(x)].
            Op::Un(UnOp::Floor) => Some([I::new(gs(a).lo.floor(), gs(a).hi.ceil())]),
            // Whole turns move neither (tan's period is half a turn).
            Op::Un(op @ (UnOp::Sin | UnOp::Cos | UnOp::Tan)) => {
                let x = off_turns(g(a).map_or(zero, |x| x.v));
                let w = way(a, x);
                interval::d_un(op, w, interval::un(op, w), [x])
            }
            Op::Un(op) => {
                let w = way(a, gs(a));
                interval::d_un(op, w, interval::un(op, w), [gs(a)])
            }
            // Over a divisor the same at each point: exact where the
            // quotient is.
            Op::Bin(BinOp::Div) if g(b).is_none() => Some([interval::bin(BinOp::Div, gs(a), vals[1]).0]),
            Op::Bin(op) => {
                let (wa, wb) = (way(a, gs(a)), way(b, gs(b)));
                interval::d_bin(op, wa, wb, interval::bin(op, wa, wb).0, crate::tape::BOTH, [gs(a)], [gs(b)])
            }
            // The condition the same at each point: its branch's difference.
            Op::Select if g(a).is_none() => Some([interval::hull(gs(b), gs(c))]),
            Op::Tag(_) | Op::Lipschitz | Op::Deriv(_) | Op::Read | Op::Real => Some([gs(a)]),
            Op::Grid(gi, mk) => {
                let w = [way(a, gs(a)), way(b, gs(b)), way(c, gs(c))];
                Some(interval::d_grid(&self.tape.grids[gi as usize], mk, w, [[gs(a)], [gs(b)], [gs(c)]]))
            }
            _ => None,
        };
        let p = &self.pairs[k];
        let sides = p.first[i - p.i].map(|f| interval::sub_o(v, f));
        Some(Gap::of(match (slope, sides) {
            (Some([x]), Some(y)) => interval::meet(x, y),
            (Some([x]), None) => x,
            (None, Some(y)) => y,
            (None, None) => I::ENTIRE,
        }))
    }

    /// Instruction `i`'s value `v`, just bounded, as the forks paired see
    /// it (those of `pairs[..upto]` whose window holds it): recorded where
    /// a first side runs, and given its arguments' bounds `vals` (as they
    /// were before it), its difference from it bounded where a second
    /// does.
    fn pair(&mut self, i: usize, v: I, upto: usize, vals: Option<[I; 3]>) {
        let o = self.code[i].out as usize;
        for k in 0..upto {
            let p = &self.pairs[k];
            if i < p.i || i > p.r {
                continue;
            }
            if !p.second {
                let at = i - p.i;
                let f = &mut self.pairs[k].first[at];
                *f = Some(f.map_or(v, |f| interval::hull(f, v)));
            } else if let Some(vals) = vals
                && k < PAIRED
            {
                let ins = &self.code[i];
                if matches!(ins.op, Op::Deriv(3)) && self.g[ins.b as usize][k].is_some_and(|x| x.bound() != I::point(0.0)) {
                    self.pairs[k].impulse = true;
                }
                self.g[o][k] = self.diff(i, k, v, vals);
            }
        }
    }

    /// Instruction `i` (not constant in the inputs) bounded first order in
    /// them: its value over the region `over`, branch `choice`.
    fn first_order(&mut self, i: usize, over: I, choice: u8, axis: Option<usize>) -> I {
        let (tape, ins) = (self.tape, &self.code[i]);
        let o = ins.out as usize;
        let (args, na) = ins.args();
        let args = &args[..na];
        let zero = I::point(0.0);
        let da = |k: usize| self.d[k].unwrap_or([zero; N]);
        let du = if let Some(j) = axis {
            Some(std::array::from_fn(|k| I::point(if k == j { 1.0 } else { 0.0 })))
        } else if args.iter().all(|&k| da(k as usize) == [zero; N]) {
            // A function of values constant in the inputs is constant in them.
            Some([zero; N])
        } else {
            self.rule(i, over, choice, da, self.inputs.axis(Var::T))
        };
        // How far it may jump: its arguments' jumps through its slope, as
        // rates along one more input reaching 1 (each side's value at a
        // point is that of a branch of it, each continuous).
        let jump = if args.iter().all(|&k| self.j[k as usize] == 0.0) {
            0.0
        } else {
            let j = |k: usize| [I::new(-self.j[k], self.j[k])];
            self.rule(i, over, choice, j, None).map_or(f64::INFINITY, |[x]| interval::mag(x))
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
        let e = interval::add_up(reach(&du, &self.hw, self.inputs.len()), jump);
        if !e.is_finite() {
            (self.v[o], self.m[o], self.d[o], self.j[o]) = (over, over, Some(du), jump);
            return over;
        }
        let at_mid = interval_ins(tape, &self.mid, i, self.code, &self.m, self.b, None, Some(&mut self.fm), false);
        let r = around(over, at_mid, e);
        (self.v[o], self.m[o], self.d[o], self.j[o]) = (r, at_mid, Some(du), jump);
        r
    }
}

impl<const N: usize> Fork for Scope<'_, N> {
    type Out = ScopeSlot<N>;
    /// Each argument's bounds over the scope and at the centre, rates, how
    /// far it may jump, and its differences from the paired forks' first
    /// sides.
    type Args = [(I, I, Option<[I; N]>, f64, [Option<Gap>; PAIRED]); 3];

    #[inline(always)]
    fn step(&mut self, i: usize, choices: Option<&mut [u8]>) -> I {
        let (tape, ctx, ins) = (self.tape, self.ctx, &self.code[i]);
        let ch = choices.expect("scope_run passes its choices");
        let o = ins.out as usize;
        let mut over = interval_ins(tape, ctx, i, self.code, &self.v, self.b, Some(&mut *ch), Some(&mut self.fv), false);
        if let Some(sl) = &mut self.slopes {
            over = sl.tighten(tape, ctx, i, self.code, self.b, over, false);
        }
        let (args, na) = ins.args();
        let axis = match ins.op {
            Op::Var(x) => self.inputs.axis(x),
            _ => None,
        };
        // The arguments' bounds before the output's slot (maybe one's) is
        // written, for the forks paired.
        let vals = (!self.pairs.is_empty()).then(|| [ins.a, ins.b, ins.c].map(|k| self.v[k as usize]));
        let link = match ins.op {
            Op::Var(Var::Param(k)) => self.chain.iter().find(|c| c.0 == k).map(|c| (c.1, c.2)),
            _ => None,
        };
        let (d, j) = (&self.d, &self.j);
        let r = if let Op::CornerReach(fa, fb) = ins.op {
            // Sub-tapes `a` and `b` may read t or ranged params directly, and
            // frame coordinates `fx, fy, fz` may move even though `Op::Frame`
            // itself carries zero rate. Where `guard >= 0` over the whole scope,
            // the op is identically 0; otherwise any input varied by `r`,
            // `guard`, `fx/fy/fz`, or either sub-tape gets an unbounded rate.
            if self.v[ins.c as usize].lo >= 0.0 {
                (self.v[o], self.m[o], self.d[o], self.j[o]) = (over, over, None, 0.0);
            } else {
                let zero = I::point(0.0);
                let f = &self.code[tape.ssa[i].b as usize];
                let deps = [ins.a, ins.c, f.a, f.b, f.c];
                let sub = Inputs::read([&*tape.forms[fa as usize], &*tape.forms[fb as usize]], ctx);
                let any_jump = deps.iter().any(|&s| j[s as usize] > 0.0);
                let mut du = [zero; N];
                let mut any_varies = false;
                for (k, slot) in du.iter_mut().enumerate().take(self.inputs.len()) {
                    let sub_reads = if self.inputs.t && k == 0 { sub.t } else { sub.params.contains(&self.inputs.params[k - self.inputs.t as usize]) };
                    let dep_moves = deps.iter().any(|&s| j[s as usize] > 0.0 || d[s as usize].is_some_and(|dr| dr[k] != zero));
                    if sub_reads || dep_moves {
                        *slot = I::ENTIRE;
                        any_varies = true;
                    }
                }
                if any_varies || any_jump {
                    (self.v[o], self.m[o], self.d[o], self.j[o]) = (over, over, Some(du), if any_jump { f64::INFINITY } else { 0.0 });
                } else {
                    (self.v[o], self.m[o], self.d[o], self.j[o]) = (over, over, None, 0.0);
                }
            }
            over
        } else if let Some((m, du)) = link {
            // Derived: from its value at the centre by its rates along the
            // inputs it is derived from.
            let e = reach(&du, &self.hw, self.inputs.len());
            let v = if e.is_finite() { around(over, m, e) } else { over };
            let m = if e.is_finite() { interval::meet(m, over) } else { over };
            (self.v[o], self.m[o], self.d[o], self.j[o]) = (v, m, Some(du), 0.0);
            v
        } else if fixed(tape, ctx, i) || !(axis.is_some() || args[..na].iter().any(|&k| d[k as usize].is_some() || j[k as usize] > 0.0)) {
            (self.v[o], self.m[o], self.d[o], self.j[o]) = (over, over, None, 0.0);
            over
        } else {
            self.first_order(i, over, ch[i], axis)
        };
        if !self.pairs.is_empty() {
            self.pair(i, r, self.pairs.len(), vals);
        }
        r
    }
    fn args(&self, i: usize) -> Self::Args {
        let ins = &self.code[i];
        [ins.a, ins.b, ins.c].map(|k| {
            let k = k as usize;
            (self.v[k], self.m[k], self.d[k], self.j[k], self.g.get(k).copied().unwrap_or([None; PAIRED]))
        })
    }
    /// A select whose condition switches at a Moment, or whose branches
    /// agree exactly at its affine clock's zero, forks: each side its branch over the region, their
    /// difference where it switches its branches' there, bounded at the
    /// moment alone (exact where they are: whole turns apart, read by a
    /// turn, it does not jump).
    fn select_sides(&mut self, i: usize, args: &Self::Args) -> Option<[I; 2]> {
        let (tape, c) = (self.tape, args[0].0);
        if !(c.lo < 0.0 && 0.0 <= c.hi) || same_branches(tape, i) || fixed(tape, self.ctx, i) {
            return None;
        }
        let ins = &tape.ssa[i];
        let gap = moment(tape, self.ctx, ins.a as usize, self.b).map(|m| {
            let at = input_at(self.ctx, m.var, m.at);
            let branch = |k: u32| eval_interval(&tape.cone(k as usize), &at, self.b, None);
            sub_o(branch(ins.c), branch(ins.b))
        });
        let zero = I::point(0.0);
        self.next_gap = Some(Gap::of(if gap == Some(zero) || join::equal(tape, self.ctx, self.b, i).is_some() { zero } else { gap? }));
        Some([args[1].0, args[2].0])
    }
    /// A select on the condition of a select forked around it, joining
    /// where that one does, takes the side running: each side of that
    /// fork is then the tape with both taken as that branch over the whole
    /// region (continuous there, bounded first order as any side is), and
    /// the tape at each point is one of the two, as both read the
    /// condition's sign there. Joining where it does, no value reading
    /// this one but the join is read past the window, so the second side
    /// overwrites nothing the rest of the run reads (the fork restores
    /// the window's slots between its sides and joins at its end). Its
    /// difference between the sides, where the seam is, is its branches'
    /// there (select_sides: each at the moment, as the second side takes
    /// it), as the fork's own is, plus its first branch's difference
    /// between the sides there, where it reads what the fork sets
    /// (`(select c (select c y 2) 2)`: y before the moment and 2 from it,
    /// its branches both 2 there).
    fn ride(&mut self, i: usize, args: &Self::Args) -> bool {
        let tape = self.tape;
        if !matches!(tape.ssa[i].op, Op::Select) {
            return false;
        }
        let (c, r) = (tape.ssa[i].a, tape.rejoin(i));
        let Some(q) = self.pairs.iter().rposition(|p| p.r == r && p.i < i && matches!(tape.ssa[p.i].op, Op::Select) && tape.ssa[p.i].a == c) else { return false };
        // Not straddled, or one branch: the step decided it.
        if self.select_sides(i, args).is_none() {
            return false;
        }
        let gap = self.next_gap.take().expect("select_sides set it");
        let k = self.pairs[q].second as usize;
        let o = self.code[i].out as usize;
        let w = &args[1 + k];
        let v = w.0;
        (self.v[o], self.m[o], self.d[o], self.j[o]) = (w.0, w.1, w.2, w.3);
        if let Some(sl) = &mut self.slopes {
            sl.0[o].v = v;
        }
        // The forks it is inside, as Fork::side records them; for its own,
        // its branches' difference.
        let y = &args[0];
        for p in 0..self.pairs.len() {
            let pair = &self.pairs[p];
            let at = i - pair.i;
            if i > pair.r {
                continue;
            }
            if !pair.second {
                let f = &mut self.pairs[p].first[at];
                *f = Some(f.map_or(v, |f| interval::hull(f, v)));
            } else if p < PAIRED {
                self.g[o][p] = if p == q {
                    Some(args[1].4[q].map_or(gap, |b| Gap::sum(BinOp::Add, gap, b)))
                } else if y.4[p].is_none() {
                    w.4[p]
                } else {
                    Some(Gap::of(pair.first[at].map_or(I::ENTIRE, |f| interval::sub_o(v, f))))
                };
            }
        }
        true
    }
    fn bounds(args: &Self::Args) -> (I, I) {
        (args[0].0, args[1].0)
    }
    fn ctx(&self) -> &Ctx {
        self.ctx
    }
    fn code(&self) -> &[crate::tape::Ins] {
        self.code
    }
    fn value(&self, k: u32) -> I {
        self.v[k as usize]
    }
    fn side(&mut self, i: usize, args: &Self::Args, k: usize, side: I) {
        let (ins, zero) = (&self.code[i], [I::point(0.0); N]);
        let o = ins.out as usize;
        // Whether an argument moves with the inputs (or may jump along them).
        let moves = |&(_, _, d, j, _): &(I, I, Option<[I; N]>, f64, [Option<Gap>; PAIRED])| j > 0.0 || d.is_some_and(|d| d != zero);
        let [y, x, _] = args;
        let (v, m, d, j, moving) = match ins.op {
            // Its branch, as it is over the region; its seam moving where
            // its condition does.
            Op::Select => {
                let w = &args[1 + k];
                (w.0, w.1, w.2, w.3, moves(y))
            }
            // One integer over the region: constant in the inputs.
            Op::Un(UnOp::Floor) => (side, side, y.2.map(|_| zero), 0.0, moves(y)),
            // The branch continued over the region, at the centre from its
            // arguments' bounds there (within those over the scope: both
            // hold there), its rates atan2's, and its arguments' jumps
            // through its slope.
            Op::Bin(BinOp::Atan2) if !moves(y) && !moves(x) && y.2.is_none() && x.2.is_none() => {
                let v = interval::atan2_branch(y.0, x.0, k == 1);
                (v, v, None, 0.0, false)
            }
            Op::Bin(BinOp::Atan2) => {
                let up = k == 1;
                let v = interval::atan2_branch(y.0, x.0, up);
                let (yc, xc) = (interval::meet(y.1, y.0), interval::meet(x.1, x.0));
                let m = if xc.hi < 0.0 { interval::meet(v, interval::atan2_branch(yc, xc, up)) } else { v };
                let du = interval::d_atan2(y.0, x.0, y.2.unwrap_or(zero), x.2.unwrap_or(zero)).unwrap_or([I::ENTIRE; N]);
                let j = if y.3 == 0.0 && x.3 == 0.0 {
                    0.0
                } else {
                    interval::d_atan2(y.0, x.0, [I::new(-y.3, y.3)], [I::new(-x.3, x.3)]).map_or(f64::INFINITY, |[x]| interval::mag(x))
                };
                let e = interval::add_up(reach(&du, &self.hw, self.inputs.len()), j);
                let moving = moves(y) | moves(x);
                if e.is_finite() { (around(v, m, e), m, Some(du), j, moving) } else { (v, v, Some(du), j, moving) }
            }
            _ => unreachable!("jump_sides takes floors and atan2, select_sides selects"),
        };
        (self.v[o], self.m[o], self.d[o], self.j[o]) = (v, m, d, j);
        // This side's range, not a meet: the bound across the cut is not
        // of this branch.
        if let Some(sl) = &mut self.slopes {
            sl.0[o].v = v;
        }
        // The forks this one is inside: its value recorded, or (by its
        // arguments' differences and its value less the first side's) its
        // difference bounded.
        let outer = if k == 0 { self.pairs.len() } else { self.pairs.len() - 1 };
        for q in 0..outer {
            let p = &self.pairs[q];
            if !p.second {
                let at = i - p.i;
                let f = &mut self.pairs[q].first[at];
                *f = Some(f.map_or(v, |f| interval::hull(f, v)));
            } else if q < PAIRED {
                self.g[o][q] = match ins.op {
                    // Its condition the same at each point: its branch's.
                    Op::Select if y.4[q].is_none() => args[1 + k].4[q],
                    Op::Un(UnOp::Floor) if y.4[q].is_none() => None,
                    Op::Bin(BinOp::Atan2) if y.4[q].is_none() && x.4[q].is_none() => None,
                    _ => {
                        let sides = p.first[i - p.i].map_or(I::ENTIRE, |f| interval::sub_o(v, f));
                        Some(Gap::of(match (ins.op, y.4[q].map(Gap::bound)) {
                            (Op::Un(UnOp::Floor), Some(x)) => interval::meet(sides, I::new(x.lo.floor(), x.hi.ceil())),
                            _ => sides,
                        }))
                    }
                };
            }
        }
        // This fork: its first side recorded; as its second starts, every
        // difference from the first cleared but the jump's own (the next
        // integer; atan2's branch a whole turn up).
        let depth = outer;
        if k == 0 {
            let r = self.tape.rejoin(i);
            let mut first = vec![None; r - i + 1];
            first[0] = Some(v);
            let gap = match ins.op {
                Op::Un(UnOp::Floor) => Gap::of(I::point(1.0)),
                Op::Select => self.next_gap.take().expect("select_sides set it"),
                _ => Gap { v: I::point(0.0), turns: 1 },
            };
            self.pairs.push(Pair { i, r, moving, second: false, first, impulse: false, gap });
        } else {
            self.pairs[depth].second = true;
            if depth < PAIRED {
                if self.g.is_empty() {
                    self.g = vec![[None; PAIRED]; self.v.len()];
                }
                self.g.iter_mut().for_each(|g| g[depth] = None);
                self.g[o][depth] = Some(self.pairs[depth].gap);
            }
        }
    }
    fn out(&self, i: usize) -> Self::Out {
        let o = self.code[i].out as usize;
        let g = self.g.get(o).copied().unwrap_or([None; PAIRED]);
        (self.v[o], self.m[o], self.d[o], self.j[o], self.slopes.as_ref().map(|sl| sl.0[o]), g)
    }
    fn put(&mut self, i: usize, (v, m, d, j, s, g): Self::Out) {
        let o = self.code[i].out as usize;
        (self.v[o], self.m[o], self.d[o], self.j[o]) = (v, m, d, j);
        if let (Some(sl), Some(s)) = (&mut self.slopes, s) {
            sl.0[o] = s;
        }
        if let Some(x) = self.g.get_mut(o) {
            *x = g;
        }
    }
    /// The fork joined is the innermost paired. Where the inputs move its
    /// seam, a value is one side's at some of their values and the other's
    /// at others: it may jump between them, by at most the second side
    /// less the first at a point (the pair's difference), on top of either
    /// side's own jumps. Each side's rates bound its slope.
    fn join(&mut self, i: usize, a: Self::Out, b: Self::Out) {
        let o = self.code[i].out as usize;
        let zero = I::point(0.0);
        let depth = self.pairs.len() - 1;
        let p = &self.pairs[depth];
        let gap = match if depth < PAIRED { b.5[depth] } else { Some(Gap::of(I::ENTIRE)) } {
            Some(x) if p.moving => interval::mag(interval::meet(x.bound(), interval::sub_o(b.0, a.0))),
            _ => 0.0,
        };
        let j = interval::add_up(a.3.max(b.3), gap);
        let (v, m) = (interval::hull(a.0, b.0), interval::hull(a.1, b.1));
        let mut d = match (a.2, b.2) {
            (None, None) if j == 0.0 => None,
            (da, db) => {
                let (da, db) = (da.unwrap_or([zero; N]), db.unwrap_or([zero; N]));
                Some(std::array::from_fn(|k| interval::hull(da[k], db[k])))
            }
        };
        // A speed of what jumps where t moves the seam: unbounded along t
        // (each side's constant too).
        if let (Some(t), true) = (self.inputs.axis(Var::T), p.moving && p.impulse) {
            d.get_or_insert([zero; N])[t] = I::ENTIRE;
        }
        (self.v[o], self.m[o], self.d[o], self.j[o]) = (v, m, d, j);
        if let (Some(sl), Some(ga), Some(gb)) = (&mut self.slopes, a.4, b.4) {
            sl.0[o] = ga.hull(&gb);
        }
        // The outer forks' differences: either side's.
        if let Some(g) = self.g.get_mut(o) {
            *g = std::array::from_fn(|k| match (a.5[k], b.5[k]) {
                _ if k >= depth => None,
                (None, None) => None,
                (Some(x), Some(y)) if x.turns == y.turns => Some(Gap { v: interval::hull(x.v, y.v), turns: x.turns }),
                (x, y) => Some(Gap::of(interval::hull(x.map_or(zero, Gap::bound), y.map_or(zero, Gap::bound)))),
            });
        }
        let r = p.r;
        self.pair(i, v, depth, None);
        if i == r {
            self.pairs.pop();
        }
    }
    fn every_value(&self) -> bool {
        self.every
    }
}
