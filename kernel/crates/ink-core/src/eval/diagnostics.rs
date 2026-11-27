//! Construction names and undefined operations on active tape paths.

use super::context::Ctx;
use super::corner::corner_reach_point;
use crate::V3;
use crate::graph::{BinOp, UnOp, Var, bin_f64, un_f64};
use crate::tape::{Op, Tape};

/// Construction names (graph tag ids) on the active path at `p`, outermost
/// first. From the root, follow the branch of each min/max/select that decides
/// the value (on a tie, the operand that is the other under a name, or the
/// one the other computes from and leaves as it is, else the first), and
/// every operand of other ops,
/// collecting the names passed. On a surface this says which named
/// construction the point belongs to: implicit faces without B-rep.
pub fn provenance(tape: &Tape, ctx: &Ctx, p: V3) -> Vec<u32> {
    let v = values(tape, ctx, p);
    names_on(tape, |ins| {
        let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
        // A side that is the other under a name names it; one that computes
        // from it and leaves it as it is does not.
        let tie = |first: usize, second: usize| match (holds(tape, &v, first, second), holds(tape, &v, second, first)) {
            (Holds::Computed, _) | (_, Holds::Plain) => second,
            _ => first,
        };
        match ins.op {
            Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax) if v[a] == v[b] => Follow::One(tie(a, b)),
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
            // Equal to what it is built from, a side leaves it as it is:
            // only that names the value.
            Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax) if v[a] == v[b] && holds(tape, &v, a, b) == Holds::Computed => Follow::One(b),
            Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax) if v[a] == v[b] && holds(tape, &v, b, a) == Holds::Computed => Follow::One(a),
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

/// How instruction `from` holds `to`'s value at the point whose values are
/// `v`: along its active path (a min's, max's or select's deciding operand,
/// through names and identities), `to` itself (Plain: `from` is it under a
/// name, which forms that face), or through an operation that computes
/// from it (Computed: a construction on it, which a tie says leaves it as
/// it is there), or not at all.
#[derive(Clone, Copy, PartialEq)]
enum Holds {
    No,
    Plain,
    Computed,
}

fn holds(tape: &Tape, v: &[f64], from: usize, to: usize) -> Holds {
    holds_in(tape, v, from, to, &mut std::collections::HashMap::new())
}

/// `holds`, each instruction once (`memo`).
fn holds_in(tape: &Tape, v: &[f64], from: usize, to: usize, memo: &mut std::collections::HashMap<usize, Holds>) -> Holds {
    if from == to {
        return Holds::Plain;
    }
    if from < to {
        return Holds::No;
    }
    if let Some(&h) = memo.get(&from) {
        return h;
    }
    let ins = &tape.ssa[from];
    let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
    let h = match ins.op {
        Op::Tag(_) | Op::Real | Op::Lipschitz | Op::Deriv(_) => holds_in(tape, v, a, to, memo),
        Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax) if v[a] == v[b] => match (holds_in(tape, v, a, to, memo), holds_in(tape, v, b, to, memo)) {
            (Holds::Plain, _) | (_, Holds::Plain) => Holds::Plain,
            (Holds::Computed, _) | (_, Holds::Computed) => Holds::Computed,
            _ => Holds::No,
        },
        Op::Bin(BinOp::Min) => holds_in(tape, v, if v[a] < v[b] { a } else { b }, to, memo),
        Op::Bin(BinOp::Max | BinOp::Emax) => holds_in(tape, v, if v[a] > v[b] { a } else { b }, to, memo),
        Op::Select => holds_in(tape, v, if v[a] < 0.0 { b } else { c }, to, memo),
        _ if reads(tape, from, to) => Holds::Computed,
        _ => Holds::No,
    };
    memo.insert(from, h);
    h
}

/// Whether instruction `from` reads `to`'s value (through any of its
/// arguments).
fn reads(tape: &Tape, from: usize, to: usize) -> bool {
    let mut stack = vec![from];
    let mut seen = std::collections::HashSet::new();
    while let Some(i) = stack.pop() {
        let (args, k) = tape.ssa[i].args();
        for &j in &args[..k] {
            let j = j as usize;
            if j == to {
                return true;
            }
            if j > to && seen.insert(j) {
                stack.push(j);
            }
        }
    }
    false
}

/// Every instruction's value at `p`.
pub(super) fn values(tape: &Tape, ctx: &Ctx, p: V3) -> Vec<f64> {
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
            Op::Within(_) => v[c],
            Op::Frame => 0.0,
            Op::Grid(gi, m) => tape.grids[gi as usize].component([v[a], v[b], v[c]], m),
            Op::CornerReach(fa, fb) => {
                if v[c] < 0.0 {
                    let f = &tape.ssa[b];
                    let q = [v[f.a as usize], v[f.b as usize], v[f.c as usize]];
                    corner_reach_point(&tape.forms[fa as usize], &tape.forms[fb as usize], ctx, v[a], q)
                } else {
                    0.0
                }
            }
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
            Op::Const(_) | Op::Var(_) | Op::Frame | Op::Within(_) | Op::CornerReach(..) => {}
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
    use rayon::prelude::*;
    let mut found = if n >= 1024 && pts.len() >= 128 {
        // Independent evaluations, combined in sample order to retain the
        // first witness and the order of cases whose counts tie.
        let chunks: Vec<_> = pts.par_chunks(64).map(|pts| undefined_chunk(tape, ctx, pts, &owner)).collect();
        let mut found: Vec<Undefined> = Vec::new();
        for u in chunks.into_iter().flatten() {
            if let Some(first) = found.iter_mut().find(|v| v.ins == u.ins && v.what == u.what) {
                first.count += u.count;
            } else {
                found.push(u);
            }
        }
        found
    } else {
        undefined_chunk(tape, ctx, pts, &owner)
    };
    found.sort_by_key(|u| std::cmp::Reverse(u.count));
    found
}

/// Points evaluated together by `undefined_chunk`: each instruction's op
/// is dispatched once for the lot, and their active paths are walked
/// together as one bit per point.
const W: usize = 16;
const _: () = assert!(W <= 64);

/// `undefined_cases` over `pts`, uncombined: W points at a time, each
/// point's values and active path its own (its lane), so each case and its
/// count and first point are as one point at a time would find them.
fn undefined_chunk(tape: &Tape, ctx: &Ctx, pts: &[V3], owner: &[Option<u32>]) -> Vec<Undefined> {
    let n = tape.ssa.len();
    let mut found: Vec<Undefined> = Vec::new();
    let mut v = vec![[0.0f64; W]; n];
    let mut bad: Vec<(usize, usize, &'static str)> = Vec::new();
    let mut active = vec![0u64; n];
    for chunk in pts.chunks(W) {
        let lanes = chunk.len();
        values_lanes(tape, ctx, chunk, &mut v);
        // The active path, from the output (consumers come after operands).
        active.iter_mut().for_each(|x| *x = 0);
        active[n - 1] = u64::MAX >> (64 - lanes);
        for i in (0..n).rev() {
            let m = active[i];
            if m == 0 {
                continue;
            }
            let ins = &tape.ssa[i];
            let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
            match ins.op {
                Op::Bin(op) => each(m, |k| {
                    if let Some(what) = crate::graph::undefined_bin(op, v[a][k], v[b][k]) {
                        bad.push((k, i, what));
                    }
                }),
                Op::Un(op) => each(m, |k| {
                    if let Some(what) = crate::graph::undefined_un(op, v[a][k]) {
                        bad.push((k, i, what));
                    }
                }),
                _ => {}
            }
            match ins.op {
                Op::Select => {
                    let below = mask(lanes, |k| v[a][k] < 0.0);
                    active[a] |= m;
                    active[b] |= m & below;
                    active[c] |= m & !below;
                }
                Op::Bin(BinOp::Min) => {
                    active[a] |= m & mask(lanes, |k| v[a][k] <= v[b][k]);
                    active[b] |= m & mask(lanes, |k| v[b][k] <= v[a][k]);
                }
                Op::Bin(BinOp::Max) => {
                    active[a] |= m & mask(lanes, |k| v[a][k] >= v[b][k]);
                    active[b] |= m & mask(lanes, |k| v[b][k] >= v[a][k]);
                }
                _ => {
                    let (args, k) = ins.args();
                    args[..k].iter().for_each(|&x| active[x as usize] |= m);
                }
            }
        }
        // Point by point, and at each point in the order the forward tape
        // reaches its cases, though activity is found walking backwards.
        bad.sort_unstable_by_key(|&(k, i, _)| (k, i));
        for (k, i, what) in bad.drain(..) {
            match found.iter_mut().find(|u| u.ins == i && u.what == what) {
                Some(u) => u.count += 1,
                None => found.push(Undefined { what, ins: i, name: owner[i], count: 1, at: chunk[k] }),
            }
        }
    }
    found
}

/// Each lane of `r` as `f` of its index.
#[inline(always)]
fn fill(r: &mut [f64; W], f: impl Fn(usize) -> f64) {
    (0..W).for_each(|k| r[k] = f(k));
}

/// `f` of each lane whose bit is set in `m`.
#[inline(always)]
fn each(mut m: u64, mut f: impl FnMut(usize)) {
    while m != 0 {
        f(m.trailing_zeros() as usize);
        m &= m - 1;
    }
}

/// The lanes below `lanes` where `f` holds, one bit each.
#[inline(always)]
fn mask(lanes: usize, f: impl Fn(usize) -> bool) -> u64 {
    (0..lanes).fold(0, |m, k| m | ((f(k) as u64) << k))
}

/// Every instruction's value at each of `pts` (up to W, one lane each; the
/// lanes past them repeat the last): `values`, the same scalar rules lane by
/// lane, so the same values to the bit.
fn values_lanes(tape: &Tape, ctx: &Ctx, pts: &[V3], v: &mut [[f64; W]]) {
    let at = |k: usize| pts[k.min(pts.len() - 1)];
    for (i, ins) in tape.ssa.iter().enumerate() {
        let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
        // Operands come before their instruction.
        let (before, rest) = v.split_at_mut(i);
        let r = &mut rest[0];
        match ins.op {
            Op::Const(x) => *r = [x; W],
            Op::Var(Var::X) => fill(r, |k| at(k)[0]),
            Op::Var(Var::Y) => fill(r, |k| at(k)[1]),
            Op::Var(Var::Z) => fill(r, |k| at(k)[2]),
            Op::Var(Var::T) => *r = [ctx.t; W],
            Op::Var(Var::Param(k)) => *r = [ctx.params[k as usize]; W],
            Op::Un(op) => {
                let x = &before[a];
                match op {
                    UnOp::Neg => (0..W).for_each(|k| r[k] = un_f64(UnOp::Neg, x[k])),
                    UnOp::Abs => (0..W).for_each(|k| r[k] = un_f64(UnOp::Abs, x[k])),
                    UnOp::Sqrt => (0..W).for_each(|k| r[k] = un_f64(UnOp::Sqrt, x[k])),
                    UnOp::Square => (0..W).for_each(|k| r[k] = un_f64(UnOp::Square, x[k])),
                    _ => (0..W).for_each(|k| r[k] = un_f64(op, x[k])),
                }
            }
            Op::Bin(op) => {
                let (x, y) = (&before[a], &before[b]);
                match op {
                    BinOp::Add => (0..W).for_each(|k| r[k] = bin_f64(BinOp::Add, x[k], y[k])),
                    BinOp::Sub => (0..W).for_each(|k| r[k] = bin_f64(BinOp::Sub, x[k], y[k])),
                    BinOp::Mul => (0..W).for_each(|k| r[k] = bin_f64(BinOp::Mul, x[k], y[k])),
                    BinOp::Div => (0..W).for_each(|k| r[k] = bin_f64(BinOp::Div, x[k], y[k])),
                    BinOp::Min => (0..W).for_each(|k| r[k] = bin_f64(BinOp::Min, x[k], y[k])),
                    BinOp::Max => (0..W).for_each(|k| r[k] = bin_f64(BinOp::Max, x[k], y[k])),
                    _ => (0..W).for_each(|k| r[k] = bin_f64(op, x[k], y[k])),
                }
            }
            Op::Select => {
                let (x, y, z) = (&before[a], &before[b], &before[c]);
                (0..W).for_each(|k| r[k] = if x[k] < 0.0 { y[k] } else { z[k] });
            }
            Op::Tag(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real => *r = before[a],
            Op::Within(_) => *r = before[c],
            Op::Frame => *r = [0.0; W],
            Op::Grid(gi, m) => {
                let g = &tape.grids[gi as usize];
                fill(r, |k| g.component([before[a][k], before[b][k], before[c][k]], m));
            }
            Op::CornerReach(fa, fb) => {
                let f = &tape.ssa[b];
                let (ta, tb) = (&tape.forms[fa as usize], &tape.forms[fb as usize]);
                fill(r, |k| {
                    if before[c][k] < 0.0 {
                        let q = [before[f.a as usize][k], before[f.b as usize][k], before[f.c as usize][k]];
                        corner_reach_point(ta, tb, ctx, before[a][k], q)
                    } else {
                        0.0
                    }
                });
            }
        }
    }
}
