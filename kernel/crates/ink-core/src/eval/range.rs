//! Plain interval evaluation, with branch choices and optional clamps.

use super::context::{Ctx, fixed, fixed_value};
use super::corner::corner_reach_interval;
use super::fork::{Fork, fork, forks};
use super::gradient::Slopes;
use super::promise::{Frames, promised};
use crate::graph::{BinOp, UnOp, Var};
use crate::interval::{self, I};
use crate::tape::{Op, Tape};

/// Bound the tape over a box. If `choices` is given, record which branch of each
/// min/max/select always wins (indexed by instruction).
pub fn eval_interval(tape: &Tape, ctx: &Ctx, b: [I; 3], choices: Option<&mut [u8]>) -> I {
    thread_local! {
        static SLOTS: std::cell::Cell<Vec<I>> = const { std::cell::Cell::new(Vec::new()) };
    }
    let mut slots = SLOTS.with(std::cell::Cell::take);
    let out = eval_interval_in(tape, ctx, b, choices, &mut slots);
    if slots.capacity() * std::mem::size_of::<I>() <= 1 << 20 {
        SLOTS.with(|scratch| scratch.set(slots));
    }
    out
}

/// Bound the tape over the points of a box where each clamped instruction's
/// value lies in its interval (`clamps`: instruction, interval): each step
/// as `eval_interval`'s (across a jump too, not forked), each clamped value
/// intersected with its interval. So over a box split by the sign of a
/// value the tape reads in several places, each part is bounded with that
/// value's sign known everywhere it is read (bounds of the box alone lose
/// that the reads agree). None where a clamp leaves a value none: no point
/// of the box is in that part. `vals` gets each instruction's bound.
pub fn eval_interval_clamped(tape: &Tape, ctx: &Ctx, b: [I; 3], clamps: &[(usize, I)], vals: &mut Vec<I>) -> Option<I> {
    let mut s = vec![I::point(0.0); tape.nslots];
    vals.clear();
    for (i, ins) in tape.code.iter().enumerate() {
        let mut v = interval_ins(tape, ctx, i, &tape.code, &s, b, None, None, true);
        for &(k, c) in clamps {
            if k == i {
                v = I { lo: v.lo.max(c.lo), hi: v.hi.min(c.hi) };
                if !(v.lo <= v.hi) {
                    return None;
                }
            }
        }
        vals.push(v);
        s[ins.out as usize] = v;
    }
    Some(s[tape.out_slot()])
}

/// `eval_interval` with the caller's storage for the slots (reused across
/// calls that evaluate many small boxes or points).
pub fn eval_interval_in(tape: &Tape, ctx: &Ctx, b: [I; 3], choices: Option<&mut [u8]>, s: &mut Vec<I>) -> I {
    s.clear();
    s.resize(tape.nslots, I::point(0.0));
    let (mut slopes, mut frames) = (Slopes::of(tape, tape.nslots), Frames::default());
    let n = tape.code.len();
    if tape.forked_steps() <= 2 * n as u64 {
        return forked_now(tape, ctx, b, choices, s, &mut slopes, &mut frames, &[]);
    }
    // First every jump bounded across, the branches it decides kept (the
    // caller's choices or our own): see forked_live.
    let mut own = Vec::new();
    let ch: &mut [u8] = match choices {
        Some(c) => c,
        None => {
            own.resize(n, crate::tape::BOTH);
            &mut own
        }
    };
    let (mut i, mut straddled) = (0, Vec::new());
    for &j in &tape.jumps {
        let j = j as usize;
        plain_run(tape, ctx, b, i, j, s, &mut slopes, &mut frames, Some(&mut *ch));
        let mut run = Plain { tape, ctx, b, s, slopes: &mut slopes, frames: &mut frames };
        let (y, x) = run.args(j);
        let v = run.step(j, Some(&mut *ch));
        if forks(tape, ctx, j, n, y, x, v, false).is_some() {
            straddled.push(j);
        }
        i = j + 1;
    }
    plain_run(tape, ctx, b, i, n, s, &mut slopes, &mut frames, Some(&mut *ch));
    if straddled.is_empty() {
        return s[tape.out_slot()];
    }
    let live = tape.live(ch);
    if !straddled.iter().any(|&j| live[j]) {
        return s[tape.out_slot()];
    }
    forked_now(tape, ctx, b, Some(ch), s, &mut slopes, &mut frames, &live)
}

/// eval_interval_in's run forking the jumps `live` (all, where empty).
#[allow(clippy::too_many_arguments)]
fn forked_now(tape: &Tape, ctx: &Ctx, b: [I; 3], mut choices: Option<&mut [u8]>, s: &mut Vec<I>, slopes: &mut Option<Slopes>, frames: &mut Frames, live: &[bool]) -> I {
    let (n, mut i) = (tape.code.len(), 0);
    for &j in &tape.jumps {
        let j = j as usize;
        // One inside the window of one before it was bounded there.
        if j >= i {
            plain_run(tape, ctx, b, i, j, s, slopes, frames, choices.as_deref_mut());
            i = fork(tape, &mut Plain { tape, ctx, b, s, slopes, frames }, j, n, choices.as_deref_mut(), live);
        }
    }
    plain_run(tape, ctx, b, i, n, s, slopes, frames, choices);
    s[tape.out_slot()]
}

/// `forked` over instructions `from..to`, which hold no jump, with Plain's
/// step written out: the slots an argument the loop keeps in registers,
/// not reloaded through the run's fields. One copy, its steps inline
/// (`interval_ins`).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn plain_run(tape: &Tape, ctx: &Ctx, b: [I; 3], from: usize, to: usize, s: &mut [I], slopes: &mut Option<Slopes>, frames: &mut Frames, mut choices: Option<&mut [u8]>) {
    for (i, ins) in tape.code[from..to].iter().enumerate() {
        let i = from + i;
        let mut v = interval_ins(tape, ctx, i, &tape.code, s, b, choices.as_deref_mut(), Some(frames), true);
        if let Some(sl) = slopes {
            v = sl.tighten(tape, ctx, i, &tape.code, b, v, true);
        }
        s[ins.out as usize] = v;
    }
}

/// eval_interval's state.
struct Plain<'a> {
    tape: &'a Tape,
    ctx: &'a Ctx,
    b: [I; 3],
    s: &'a mut Vec<I>,
    slopes: &'a mut Option<Slopes>,
    frames: &'a mut Frames,
}

impl Fork for Plain<'_> {
    type Out = (I, Option<interval::G>);
    type Args = (I, I);
    #[inline(always)]
    fn step(&mut self, i: usize, choices: Option<&mut [u8]>) -> I {
        let ins = &self.tape.code[i];
        let mut v = interval_ins(self.tape, self.ctx, i, &self.tape.code, self.s, self.b, choices, Some(self.frames), true);
        if let Some(sl) = self.slopes {
            v = sl.tighten(self.tape, self.ctx, i, &self.tape.code, self.b, v, true);
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
    fn code(&self) -> &[crate::tape::Ins] {
        &self.tape.code
    }
    fn value(&self, k: u32) -> I {
        self.s[k as usize]
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

/// Bound `v` of instruction `i`, met with what it is where it is a
/// remainder a - p floor(a / p) (Tape::remainder, `slot` the bound in a
/// slot) of p one nonzero number: between 0 and p, as in real numbers
/// floor(a / p) = k puts a between k p and (k + 1) p. Bounds of the
/// formula read a twice and lose that where the floor takes several
/// values.
pub(super) fn remainder_within(tape: &Tape, i: usize, slot: impl Fn(usize) -> I, v: I) -> I {
    match tape.remainder(i).map(|k| slot(k as usize)) {
        Some(p) if p.lo == p.hi && p.lo != 0.0 && p.lo.is_finite() => interval::meet(v, I::new(p.lo.min(0.0), p.lo.max(0.0))),
        _ => v,
    }
}

/// Instruction `i` (`ins`, in the tape's slot code or its SSA form)'s bound
/// from its arguments' bounds in `s`. `frames`, by slot beside `s`, keeps
/// each frame's extent over the box, the diagonal of its coordinates'
/// bounds (set where a frame is bounded, read where a promise is): apart
/// from the values, so no rule of theirs touches it. None: a promise is
/// read as its value. `slots`: `s` holds the tape's code's slots, so a
/// min or max of a and a + c reads c's (Tape::apart) and is ordered by
/// its sign. Inline in the loops that run it an instruction at a time:
/// called, each instruction saved and restored 18 registers on an M4, a
/// third of its own time on the Falcon swerve check.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) fn interval_ins(
    tape: &Tape,
    ctx: &Ctx,
    i: usize,
    code: &[crate::tape::Ins],
    s: &[I],
    b: [I; 3],
    choices: Option<&mut [u8]>,
    frames: Option<&mut Frames>,
    slots: bool,
) -> I {
    let ins = &code[i];
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
            let (a, b) = (s[ins.a as usize], s[ins.b as usize]);
            // a and a + c: ordered by c's sign, where bounds of each
            // (both reading a) may overlap.
            let apart = if slots && matches!(op, BinOp::Min | BinOp::Max) { tape.apart(i) } else { None };
            let order = apart.map_or(0.0, |(c, neg)| {
                let c = s[c as usize];
                let d = if neg { I::new(-c.hi, -c.lo) } else { c };
                if d.lo > 0.0 {
                    1.0
                } else if d.hi < 0.0 {
                    -1.0
                } else {
                    0.0
                }
            });
            let (v, c) = match (op, order) {
                (BinOp::Min, 1.0) | (BinOp::Max, -1.0) => (a, crate::tape::FIRST),
                (BinOp::Min, -1.0) | (BinOp::Max, 1.0) => (b, crate::tape::SECOND),
                _ => interval::bin(op, a, b),
            };
            if let Some(ch) = choices {
                ch[i] = c;
            }
            if slots && op == BinOp::Sub { remainder_within(tape, i, |k| s[k], v) } else { v }
        }
        Op::Select => {
            let (v, c) = interval::select(s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]);
            if let Some(ch) = choices {
                ch[i] = c;
            }
            v
        }
        Op::Lipschitz => {
            let (v, k) = (s[ins.a as usize], s[ins.b as usize]);
            // Where its promise could bound it narrower: that bound is at
            // least k times its frame's extent wide (`promised`).
            if let Some(f) = frames
                && k.lo > 0.0
                && f.widths.get(ins.c as usize).is_some_and(|e| v.hi - v.lo > k.lo * e)
            {
                promised(tape, ctx, i, b, f).map_or(v, |p| interval::meet(v, p))
            } else {
                v
            }
        }
        Op::Tag(_) | Op::Read | Op::Deriv(_) | Op::Real => s[ins.a as usize],
        Op::Within(_) => I::new(s[ins.a as usize].lo, s[ins.b as usize].hi),
        Op::Frame => {
            if let Some(f) = frames {
                if f.widths.len() < s.len() {
                    f.widths.resize(s.len(), 0.0);
                }
                f.widths[ins.out as usize] = interval::norm([ins.a, ins.b, ins.c].map(|k| s[k as usize].hi - s[k as usize].lo));
            }
            I::point(0.0)
        }
        Op::Grid(gi, m) => {
            let (g, b) = (&tape.grids[gi as usize], [s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]]);
            g.component_bound(b, m)
        }
        Op::CornerReach(fa, fb) => {
            let f = &code[tape.ssa[i].b as usize];
            let q = [s[f.a as usize], s[f.b as usize], s[f.c as usize]];
            corner_reach_interval(&tape.forms[fa as usize], &tape.forms[fb as usize], ctx, s[ins.a as usize], q, s[ins.c as usize])
        }
    }
}
