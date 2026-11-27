//! Forking interval evaluation across jumps and tracking differences between sides.

use super::context::{Ctx, fixed};
use crate::graph::{BinOp, UnOp};
use crate::interval::{self, I};
use crate::tape::{Op, Tape};

/// An interval evaluator's state over a region, which a jump forks
/// (`forked`): its slots, and whatever it keeps beside them.
pub(super) trait Fork {
    type Out: Copy;
    /// Instruction `i`'s arguments as they are before it is stepped (its
    /// output's slot may be one's): what a side of a jump there is taken
    /// from.
    type Args;
    /// Instruction `i` bounded and stored.
    fn step(&mut self, i: usize, choices: Option<&mut [u8]>) -> I;
    /// The instructions it runs (the tape's code, or its SSA form).
    fn code(&self) -> &[crate::tape::Ins];
    /// The bound slot `k` holds.
    fn value(&self, k: u32) -> I;
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
    /// The sides of select `i` (its arguments `args`), as a jump, where
    /// this evaluator forks one: None, but for scope_run's.
    fn select_sides(&mut self, _i: usize, _args: &Self::Args) -> Option<[I; 2]> {
        None
    }
    /// Whether select `i` (its arguments `args`) is decided by a fork it is
    /// inside, and if so it is set to that fork's side: one on the same
    /// condition joining where it does (Tape::rejoin). At each point both
    /// read one sign of it, so each side of the fork is the tape with
    /// every select on it taking that side: forking each again would bound
    /// combinations of sides no point takes, twice as many per select.
    /// None but scope_run's forks selects.
    fn ride(&mut self, _i: usize, _args: &Self::Args) -> bool {
        false
    }
    /// Whether every value is read after the run (an evaluator over the
    /// tape's SSA form), not only the output: then every value a fork sets
    /// is joined (the jump's own, taken as each side, and each its window
    /// computes), not only the one it joins at.
    fn every_value(&self) -> bool {
        false
    }
    /// Whether a floor the box straddles more than once is forked too, its
    /// sides runs of whole numbers (jump_sides): for an evaluator that
    /// gains by it (Grad: a slope bound needs only that the field is
    /// continuous at each crossing), not one bounding values, where the
    /// runs' hull is the floor's own bound.
    fn many(&self) -> bool {
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
/// holds where both sides make it. Only jumps `live` are forked (all,
/// where it is empty; see forked_live).
pub(super) fn forked<F: Fork>(tape: &Tape, run: &mut F, from: usize, to: usize, mut choices: Option<&mut [u8]>, live: &[bool]) {
    let mut i = from;
    let first = tape.jumps.partition_point(|&j| (j as usize) < from);
    for &j in &tape.jumps[first..] {
        let j = j as usize;
        if j >= to {
            break;
        }
        if j >= i {
            (i..j).for_each(|k| _ = run.step(k, choices.as_deref_mut()));
            i = fork(tape, run, j, to, choices.as_deref_mut(), live);
        }
    }
    (i..to).for_each(|k| _ = run.step(k, choices.as_deref_mut()));
}

/// `forked` from a jump instruction `i` (Tape::jumps): bounded, and where
/// the region straddles it once (`jump_sides`; a select, `select_sides`), its window bounded once
/// with each side and joined; unless its join lies past `to`. Where the run
/// goes on. Out of line, as jumps are rare: the loop over the rest stays as
/// tight as one with none.
#[inline(never)]
pub(super) fn fork<F: Fork>(tape: &Tape, run: &mut F, i: usize, to: usize, mut choices: Option<&mut [u8]>, live: &[bool]) -> usize {
    // Before the output's slot (which may be an argument's) is written.
    let args = run.args(i);
    let (y, x) = F::bounds(&args);
    let v = run.step(i, choices.as_deref_mut());
    if live.get(i) == Some(&false) || run.ride(i, &args) {
        return i + 1;
    }
    // Forks inside a fork's window bound it twice each: nested past FORK_DEPTH they would double the
    // work at every level (a polar repeat of rollers, each read by a paint that reads them all, nested
    // one per roller), and are bounded across instead, as a jump straddled twice is: a wider bound.
    if FORK_DEPTH_NOW.with(|d| d.get()) >= FORK_DEPTH {
        return i + 1;
    }
    let Some((sides, r)) = sides_of(tape, run, i, to, &args, y, x, v) else { return i + 1 };
    FORKED.with(|c| c.set(c.get() + 1));
    FORK_DEPTH_NOW.with(|d| d.set(d.get() + 1));
    // What the window's slots held before it: the values it reads from
    // before the jump, which the first side may overwrite. One per slot
    // (many of the window's instructions share one).
    let mut slots: Vec<(u32, usize)> = (i + 1..=r).map(|k| (run.code()[k].out, k)).collect();
    slots.sort_unstable();
    slots.dedup_by_key(|x| x.0);
    let keep: Vec<(usize, F::Out)> = slots.iter().map(|&(_, k)| (k, run.out(k))).collect();
    run.side(i, &args, 0, sides[0]);
    forked(tape, run, i + 1, r + 1, choices.as_deref_mut(), live);
    // The value it joins at, and where every value is read, every value
    // the sides set: the jump's own (one side each) and its window's.
    let first = run.out(r);
    let window: Vec<F::Out> = if run.every_value() { (i..r).map(|k| run.out(k)).collect() } else { Vec::new() };
    let picked = choices.as_deref().map(|ch| ch[i + 1..=r].to_vec());
    keep.iter().for_each(|&(k, o)| run.put(k, o));
    run.side(i, &args, 1, sides[1]);
    forked(tape, run, i + 1, r + 1, choices.as_deref_mut(), live);
    for (k, a) in (i..r).zip(window) {
        let b = run.out(k);
        run.join(k, a, b);
    }
    FORK_DEPTH_NOW.with(|d| d.set(d.get() - 1));
    let second = run.out(r);
    run.join(r, first, second);
    if let (Some(ch), Some(p)) = (choices, picked) {
        ch[i + 1..=r].iter_mut().zip(p).for_each(|(c, p)| *c = if *c == p { p } else { crate::tape::BOTH });
    }
    r + 1
}

thread_local! {
    static FORKED: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
    /// How many forks the one being bounded on this thread sits inside.
    static FORK_DEPTH_NOW: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
}

/// The most forks one may sit inside: past it a jump is bounded across (`fork`). Not tuned, and
/// nothing measured separates it from its neighbours: on the catalog's sample and the drafted outlines
/// (581 files, 10,621 claims) 1, 2, 4 and no limit answer every claim alike, and their times differ
/// by less than the machine's noise; what grows with it is the worst case, a mecanum wheel's paint
/// (89, 108, 150 and 234 s at 0, 2, 4 and 6 deep; no limit, hours).
const FORK_DEPTH: u32 = 2;

/// How many jumps interval bounds have forked on this thread (each bounds
/// its window twice): the work a fork costs, for tests to hold.
pub fn forked_count() -> u64 {
    FORKED.with(|c| c.get())
}

/// Where jump `i` (bounds of its arguments `y`, `x`, of its value `v`) is
/// forked in a run up to `to`: its sides (jump_sides) and the instruction
/// they join at, before `to`.
#[allow(clippy::too_many_arguments)]
pub(super) fn forks(tape: &Tape, ctx: &Ctx, i: usize, to: usize, y: I, x: I, v: I, many: bool) -> Option<([I; 2], usize)> {
    let r = tape.rejoin(i);
    if r <= i || r >= to {
        return None;
    }
    jump_sides(tape, ctx, i, &tape.code[i], y, x, v, many).map(|s| (s, r))
}

/// `forks` for jump `i` of `run`, a select's sides from the run
/// (Fork::select_sides): it joins at itself where it is the output, its
/// sides' gap its jump.
#[allow(clippy::too_many_arguments)]
fn sides_of<F: Fork>(tape: &Tape, run: &mut F, i: usize, to: usize, args: &F::Args, y: I, x: I, v: I) -> Option<([I; 2], usize)> {
    if !matches!(tape.code[i].op, Op::Select) {
        return forks(tape, run.ctx(), i, to, y, x, v, run.many());
    }
    let r = tape.rejoin(i);
    if r < i || r >= to {
        return None;
    }
    run.select_sides(i, args).map(|s| (s, r))
}

/// `forked` over all of `run`'s code, its jumps forked only where the
/// output reads them. A fork bounds its window twice, so forks inside it
/// four times, and a window reaching the output (a jump's value read twice,
/// one read through a negation: a screw's thread, and the same thread cut
/// from its nut) holds every jump after it: a cell beside a row of tapped
/// holes, on the half-plane each hole's atan2 cuts at, forks once per hole,
/// nested. So the first pass bounds every jump across, noting the min, max
/// and select branches it decides. A jump the output reads only through
/// branches decided against it (Tape::live) is bounded across again: no
/// value the output reads depends on it, and the bounds a fork narrows
/// keep the same branches decided, so the output's bound is the forked
/// one. The rest are forked in a second pass; where none straddled, the
/// first is the run. Where forking every jump would take at most two runs
/// (Tape::forked_steps), as many as the first pass and one more, all are
/// forked at once.
pub(super) fn forked_live<F: Fork>(tape: &Tape, run: &mut F, mut choices: Option<&mut [u8]>) {
    use crate::tape::BOTH;
    let n = run.code().len();
    if tape.forked_steps() <= 2 * n as u64 {
        forked(tape, run, 0, n, choices, &[]);
        return;
    }
    let (mut picks, mut straddled, mut next) = (vec![BOTH; n], Vec::new(), 0);
    for i in 0..n {
        let ins = run.code()[i];
        let v = |k: u32| run.value(k);
        match ins.op {
            Op::Bin(op @ (BinOp::Min | BinOp::Max | BinOp::Emax)) => picks[i] = interval::bin(op, v(ins.a), v(ins.b)).1,
            Op::Select => picks[i] = interval::select(v(ins.a), v(ins.b), v(ins.c)).1,
            _ => {}
        }
        if tape.jumps.get(next) == Some(&(i as u32)) {
            next += 1;
            let args = run.args(i);
            let (y, x) = F::bounds(&args);
            let v = run.step(i, choices.as_deref_mut());
            if sides_of(tape, run, i, n, &args, y, x, v).is_some() {
                straddled.push(i);
            }
        } else {
            run.step(i, choices.as_deref_mut());
        }
    }
    if straddled.is_empty() {
        return;
    }
    let live = tape.live(&picks);
    if straddled.iter().any(|&j| live[j]) {
        forked(tape, run, 0, n, choices, &live);
    }
}

/// The two sides of a jump instruction `i` straddles once over the box
/// (bound `v`; an atan2's arguments' bounds `y`, `x`), each a bound of its value
/// where it takes that side, continuous there: a floor whose value is one
/// of two integers, each; atan2 across its cut away from the origin (x < 0
/// throughout), near -pi where y < 0 and near pi where y >= 0 (atan2 is odd
/// in y). None for any other instruction or box.
///
/// A floor's integers are its bound's ends rounded inward: a bound met with
/// another (a first-order one, rounded outward) may end between integers,
/// [-5e-324, 1] for a floor of 0 or 1, and a side at -5e-324 bounds what
/// reads it at a value the floor never takes (atan2 of it near -pi, not pi).
///
/// With `many`, a floor straddling several integers lo..=hi (whole numbers
/// a float holds exactly, under 2^52) has the sides lo..=hi-1 and
/// lo+1..=hi: each a bound of its value wherever it takes one of them,
/// constant between its crossings.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn jump_sides(tape: &Tape, ctx: &Ctx, i: usize, ins: &crate::tape::Ins, y: I, x: I, v: I, many: bool) -> Option<[I; 2]> {
    if fixed(tape, ctx, i) {
        return None;
    }
    match ins.op {
        Op::Un(UnOp::Floor) => {
            // + 0.0: one zero (the ceiling of -5e-324 is -0).
            let (lo, hi) = (v.lo.ceil() + 0.0, v.hi.floor() + 0.0);
            if hi - lo == 1.0 {
                return Some([I::point(lo), I::point(hi)]);
            }
            let whole = |k: f64| k.abs() < (1u64 << 52) as f64;
            (many && hi - lo > 1.0 && whole(lo) && whole(hi)).then(|| [I::new(lo, hi - 1.0), I::new(lo + 1.0, hi)])
        }
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

/// Forks whose sides are paired (the innermost PAIRED a run is inside).
pub(super) const PAIRED: usize = 4;

/// A difference between a fork's sides: `v` plus `turns` whole turns (2 pi
/// each, pi exact; a whole number, counted in integers so no sum of them
/// rounds), so a whole turn of degrees read in radians (UnOp::Rad) or
/// atan2's branches stay exact to the sine and cosine that read them.
/// What is not shown a whole number of turns goes into `v`, rounded
/// outward.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Gap {
    pub(super) v: I,
    pub(super) turns: i64,
}

impl Gap {
    pub(super) fn of(v: I) -> Gap {
        Gap { v, turns: 0 }
    }
    /// As an interval (the count and 2 pi rounded outward).
    pub(super) fn bound(self) -> I {
        use std::f64::consts::TAU;
        if self.turns == 0 {
            return self.v;
        }
        let k = self.turns as f64;
        let k = if k as i128 == self.turns as i128 { I::point(k) } else { I::new(k.next_down(), k.next_up()) };
        interval::add_o(self.v, interval::mul_o(k, I::new(TAU, TAU.next_up())))
    }
    /// The whole turns `d` degrees are, where d is a float that is 360 k
    /// for a whole k (shown in integers, not by a rounded quotient: 1e20
    /// is 280 degrees past a whole turn).
    fn of_degrees(d: f64) -> Option<i64> {
        if !(d.is_finite() && d == d.trunc() && d.abs() < 1e38) {
            return None;
        }
        let d = d as i128;
        if d % 360 != 0 {
            return None;
        }
        i64::try_from(d / 360).ok()
    }
    /// `turns` turns in degrees, where that is a float.
    fn degrees(turns: i64) -> Option<f64> {
        let d = 360 * turns as i128;
        let f = d as f64;
        (f as i128 == d).then_some(f)
    }
    /// a + b or a - b: the counts added where their sum is a whole number
    /// an i64 holds, else both as bounds.
    pub(super) fn sum(op: BinOp, a: Gap, b: Gap) -> Gap {
        let turns = if op == BinOp::Add { a.turns.checked_add(b.turns) } else { a.turns.checked_sub(b.turns) };
        match turns {
            Some(turns) => Gap { v: interval::bin(op, a.v, b.v).0, turns },
            None => Gap::of(interval::bin(op, a.bound(), b.bound()).0),
        }
    }

    /// Instruction `ins`'s difference between the sides, from its
    /// arguments' (`g`, None where the same), where whole turns carry it
    /// exactly: through sums; from degrees that are a whole number of 360s
    /// into radians and back; and none left by a sine, cosine or tangent of
    /// whole turns alone. None where no rule here gives it.
    pub(super) fn exact(ins: &crate::tape::Ins, g: impl Fn(u32) -> Option<Gap>) -> Option<Gap> {
        let zero = I::point(0.0);
        let turns = |s: u32| g(s).map_or(0, |x| x.turns);
        let gap = |s: u32| g(s).unwrap_or(Gap::of(zero));
        match ins.op {
            Op::Un(UnOp::Rad) => {
                let x = gap(ins.a);
                if x.turns == 0
                    && x.v.lo == x.v.hi
                    && let Some(turns) = Gap::of_degrees(x.v.lo)
                {
                    return Some(Gap { v: zero, turns });
                }
                None
            }
            Op::Un(UnOp::Deg) if gap(ins.a).v == zero && turns(ins.a) != 0 => Gap::degrees(turns(ins.a)).map(|d| Gap::of(I::point(d))),
            Op::Un(UnOp::Sin | UnOp::Cos | UnOp::Tan) if gap(ins.a).v == zero => Some(Gap::of(zero)),
            Op::Bin(op @ (BinOp::Add | BinOp::Sub)) if turns(ins.a) != 0 || turns(ins.b) != 0 => Some(Gap::sum(op, gap(ins.a), gap(ins.b))),
            Op::Un(UnOp::Neg) if turns(ins.a) != 0 => Some(Gap::sum(BinOp::Sub, Gap::of(zero), gap(ins.a))),
            Op::Tag(_) | Op::Lipschitz | Op::Deriv(_) | Op::Read | Op::Real if turns(ins.a) != 0 => g(ins.a),
            _ => None,
        }
    }
}

/// Whether select instruction `i`'s two branches are the same field, up to
/// names (tags): then it is continuous where its condition changes sign.
pub(super) fn same_branches(tape: &Tape, i: usize) -> bool {
    let ssa = &tape.ssa;
    let strip = |mut j: usize| {
        while let Op::Tag(_) = ssa[j].op {
            j = ssa[j].a as usize;
        }
        j
    };
    strip(ssa[i].b as usize) == strip(ssa[i].c as usize)
}
