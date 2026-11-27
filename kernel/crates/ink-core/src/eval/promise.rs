//! Bounds from Lipschitz promises and checks that a promise applies.

use super::context::Ctx;
use super::range::interval_ins;
use crate::graph::Var;
use crate::interval::{self, I, add_o, sub_o};
use crate::tape::{Op, Tape};

/// A run's frame widths by live slot and its promises' values by SSA
/// instruction. The latter hold throughout forks: a promise reads its own
/// plain bounds for this region and its centre, regardless of a fork's side.
#[derive(Default)]
pub(super) struct Frames {
    pub(super) widths: Vec<f64>,
    values: Option<PromiseMemo>,
}

impl Drop for Frames {
    fn drop(&mut self) {
        if let Some(mut values) = self.values.take() {
            // Keep one reusable allocation per worker, within a small
            // fixed budget. Each run still starts with no values read.
            let bytes = (values.at.bounds.capacity() + values.over.bounds.capacity()) * std::mem::size_of::<I>();
            if bytes <= 8 << 20 {
                values.at.read.fill(false);
                values.over.read.fill(false);
                values.walk.clear(); // a run a panic cut short leaves its walk
                let _ = PROMISE_MEMO.try_with(|scratch| scratch.set(values));
            }
        }
    }
}

#[derive(Default)]
struct PromiseMemo {
    over: PromiseValues,
    at: PromiseValues,
    walk: Vec<(u32, bool)>,
}

thread_local! {
    static PROMISE_MEMO: std::cell::Cell<PromiseMemo> = Default::default();
}

/// Plain bounds computed on demand, once per instruction and region even
/// when several promises read it. Values remain in SSA order so a shared
/// operand can be reused wherever it is asked for.
#[derive(Default)]
struct PromiseValues {
    bounds: Vec<I>,
    read: Vec<bool>,
}

impl PromiseValues {
    fn read(&mut self, tape: &Tape, ctx: &Ctx, b: [I; 3], root: u32, walk: &mut Vec<(u32, bool)>) -> I {
        let n = root as usize + 1;
        if self.bounds.len() < n {
            self.bounds.resize(n, I::point(0.0));
            self.read.resize(n, false);
        }
        if !self.read[root as usize] {
            walk.clear();
            walk.push((root, false));
            while let Some((j, operands_ready)) = walk.pop() {
                if self.read[j as usize] {
                    continue;
                }
                let ins = &tape.ssa[j as usize];
                if operands_ready {
                    self.bounds[j as usize] = cone_ins(tape, ctx, j, ins, &self.bounds, b);
                    self.read[j as usize] = true;
                } else {
                    let (args, k) = ins.args();
                    walk.push((j, true));
                    walk.extend(args[..k].iter().filter(|&&a| !self.read[a as usize]).map(|&a| (a, false)));
                }
            }
        }
        self.bounds[root as usize]
    }
}

/// Lipschitz instruction `i`'s bound over the box `b` from its promise:
/// `(lipschitz k f)` moves f by at most k times as far as its frame's
/// coordinates move, so over the box f lies within k d of its value at the
/// box's centre, d the farthest the coordinates over the box lie from
/// theirs at the centre. Its own formula's bound (`v`) can stay wide however
/// small the box (a normal divided by the length of a short chord whose
/// two ends each move: a drafted polygon's line through two points of a
/// small arc), where this falls with the box. It is at least k times the
/// frame's extent wide (d is at least half the diagonal of the coordinates'
/// bounds over the box), so it is tried where the formula's is wider. Read again
/// from the tape's SSA form, in plain intervals (no jump forked, a
/// promise within f read as its value): k and the coordinates
/// over the box, f and the coordinates at the centre. None where it promises
/// nothing: k not positive, its constants not real over the box
/// (`promise_valid`), or a bound not finite.
#[cold]
#[inline(never)]
pub(super) fn promised(tape: &Tape, ctx: &Ctx, i: usize, b: [I; 3], frames: &mut Frames) -> Option<I> {
    let centre = b.map(|x| I::point(x.mid()));
    if !centre.iter().all(|c| c.lo.is_finite()) {
        return None;
    }
    let values = frames.values.get_or_insert_with(|| PROMISE_MEMO.with(std::cell::Cell::take));
    let ins = tape.ssa[i];
    debug_assert!(matches!(ins.op, Op::Lipschitz), "a promise");
    let f = tape.ssa[ins.c as usize];
    let coords = [f.a, f.b, f.c];
    let k = values.over.read(tape, ctx, b, ins.b, &mut values.walk);
    if !(k.lo > 0.0 && promise_valid(tape, ctx, i, b)) {
        return None;
    }
    let over = coords.map(|j| values.over.read(tape, ctx, b, j, &mut values.walk));
    let at = coords.map(|j| values.at.read(tape, ctx, centre, j, &mut values.walk));
    let d = interval::norm(std::array::from_fn(|j| interval::mag(sub_o(over[j], at[j]))));
    let r = interval::mul_up(k.hi, d);
    let f = add_o(values.at.read(tape, ctx, centre, ins.a, &mut values.walk), I::new(-r, r));
    (f.lo.is_finite() && f.hi.is_finite()).then_some(f)
}

/// Instruction `j` of a promise's cone (`ins`, its arguments renumbered
/// into `s`) bounded over `b`: a promise within it read as its value.
fn cone_ins(tape: &Tape, ctx: &Ctx, j: u32, ins: &crate::tape::Ins, s: &[I], b: [I; 3]) -> I {
    match ins.op {
        Op::Lipschitz => s[ins.a as usize],
        _ => interval_ins(tape, ctx, j as usize, &tape.ssa, s, b, None, None, false),
    }
}

/// Whether promise `i` (a Lipschitz instruction of the tape's SSA form)
/// holds over the box `b`: its field's, k's and frame's constants are real
/// numbers there (`real_flags`, `promise_holds`). Every bound that trusts a
/// promise asks this: values' (`promised`) and gradients' (`grad_ins`).
/// Each value the flags read is bounded over the box, in plain intervals:
/// a select's condition, and each value that reads no x, y or z through
/// what it takes (the flags of its arguments say; a select reads the branch
/// its condition decides).
pub(super) fn promise_valid(tape: &Tape, ctx: &Ctx, i: usize, b: [I; 3]) -> bool {
    // k read at the frame's image bounds no slope along a way that leaves it (Promise::leaps).
    if tape.promise(i).leaps {
        return false;
    }
    if tape.promise_always_real(i) {
        return true;
    }
    let p = tape.promise(i);
    // Each cone instruction's bound and flags at its index in the SSA form: the cone's own
    // positions are written before they are read (its arguments come before it), the rest never.
    // Taken out of the thread's scratch while in use: a nested call takes an empty one.
    let (mut s, mut f) = PROMISE_SLOTS.with(std::cell::Cell::take);
    if s.len() < tape.ssa.len() {
        s.resize(tape.ssa.len(), I::point(0.0));
        f.resize(tape.ssa.len(), 0);
    }
    for (n, &j) in p.ins.iter().enumerate() {
        let ins = tape.ssa[j as usize];
        let j = j as usize;
        s[j] = I::point(0.0);
        if p.cond[n] || real_flags(&ins, |k| s[k as usize], &f, I::point(0.0)) & SPATIAL == 0 {
            s[j] = cone_ins(tape, ctx, j as u32, &ins, &s, b);
        }
        f[j] = real_flags(&ins, |k| s[k as usize], &f, s[j]);
    }
    let holds = promise_holds(&tape.ssa[i], &f);
    PROMISE_SLOTS.with(|c| c.set((s, f)));
    holds
}

thread_local! {
    /// promise_valid's slots, by SSA index (its cone's positions only are read).
    static PROMISE_SLOTS: std::cell::Cell<(Vec<I>, Vec<u8>)> = const { std::cell::Cell::new((Vec::new(), Vec::new())) };
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
/// `s` holds the arguments' bounds, `v` the value's.
fn real_flags(ins: &crate::tape::Ins, s: impl Fn(u32) -> I, f: &[u8], v: I) -> u8 {
    let (args, n) = ins.args();
    let of = |k: usize| f[args[k] as usize];
    let mut out = match ins.op {
        Op::Var(Var::X | Var::Y | Var::Z) => return SPATIAL,
        // Of the taken branch, where the condition decides one.
        Op::Select => {
            let c = s(ins.a);
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
        let arg = |k: usize| {
            let r = s(args[k]);
            [r.lo, r.hi]
        };
        // The one table of operator domains (graph::undefined_un_over).
        let undefined = match ins.op {
            Op::Un(op) => crate::graph::undefined_un_over(op, arg(0)).is_some(),
            Op::Bin(op) => crate::graph::undefined_bin_over(op, arg(0), arg(1)).is_some(),
            _ => false,
        };
        if undefined || !(v.lo.is_finite() && v.hi.is_finite()) {
            out |= UNREAL;
        }
    }
    out
}

/// Whether a promise, Lipschitz instruction `ins`, holds over a region: its
/// field's, k's and frame's constants are real numbers there (their
/// `real_flags` in `f`). Gradient bounds and value bounds trust it alike.
fn promise_holds(ins: &crate::tape::Ins, f: &[u8]) -> bool {
    (f[ins.a as usize] | f[ins.b as usize] | f[ins.c as usize]) & UNREAL == 0
}
