//! The batched point evaluators against one point at a time: values,
//! gradients and undefined cases the same to the bit, however many points
//! share a batch and wherever an output reuses an operand's slot.

use super::diagnostics::{undefined_cases, values};
use super::{Ctx, L, eval_dual, eval_grads, eval_point, eval_points};
use crate::V3;
use crate::graph::{BinOp, Graph, UnOp, Var};
use crate::tape::{Op, Tape};

/// A field with each undefined case the walk reports, some guarded (a
/// select's untaken branch, a min's or max's losing side, ties between
/// them), named, and long enough that outputs reuse operand slots.
fn field() -> (Graph, Tape) {
    let mut g = Graph::new();
    let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
    let one = g.constant(1.0);
    let half = g.constant(0.5);
    let ratio = g.bin(BinOp::Div, y, x);
    let ratio = g.tag("ratio", ratio);
    let zz = g.un(UnOp::Square, z);
    let w = g.sub(one, zz);
    let log = g.un(UnOp::Ln, w);
    let log = g.tag("log", log);
    // ln x where it wins the min (x <= y), or ties it.
    let lx = g.un(UnOp::Ln, x);
    let low = g.min(lx, y);
    // 1 / (y - 1) only where x < 0 selects it, ln z only where not.
    let ym = g.sub(y, one);
    let inv = g.bin(BinOp::Div, one, ym);
    let lz = g.un(UnOp::Ln, z);
    let lz = g.tag("else", lz);
    let guarded = g.select(x, inv, lz);
    let guarded = g.tag("guarded", guarded);
    let root = g.bin(BinOp::Pow, z, half);
    let high = g.max(root, x);
    let mut f = g.add(ratio, log);
    for t in [low, guarded, high] {
        let s = g.mul(f, t);
        f = g.add(s, t);
    }
    let sq = g.un(UnOp::Sqrt, f);
    let f = g.sub(sq, x);
    let tape = Tape::new(&g, f);
    (g, tape)
}

/// Points on a grid of zeros, signs and ties, more than a batch, not a
/// whole number of them.
fn points() -> Vec<V3> {
    let v = [-1.0, -0.0, 0.0, 0.5, 1.0, 2.0];
    let mut pts = Vec::new();
    for a in v {
        for b in v {
            for c in v {
                pts.push([a, b, c]);
            }
        }
    }
    pts.extend([[1.0, 0.0, 1.0]; 3]);
    pts
}

#[test]
fn the_test_field_reuses_operand_slots() {
    let (_, tape) = field();
    assert!(tape.code.iter().any(|ins| matches!(ins.op, Op::Bin(_) | Op::Un(_)) && (ins.out == ins.a || ins.out == ins.b)));
}

#[test]
fn batched_gradients_are_each_points_own() {
    let (_, tape) = field();
    let ctx = Ctx::new(0.0, vec![]);
    let pts = points();
    assert!(pts.len() > L && !pts.len().is_multiple_of(L));
    // Every lane width eval_grads picks, and points whose batch is full or not.
    for n in [1, 2, 8, 9, 16, 17, 32, 33, L, pts.len()] {
        let got = eval_grads(&tape, &ctx, &pts[..n]);
        for (p, g) in pts[..n].iter().zip(got) {
            let d = eval_dual(&tape, &ctx, *p);
            assert_eq!(g.map(f64::to_bits), [d[0], d[1], d[2], d[3]].map(f64::to_bits), "at {p:?} in {n}");
        }
    }
}

#[test]
fn batched_values_are_each_points_own() {
    let (_, tape) = field();
    let ctx = Ctx::new(0.0, vec![]);
    let pts = points();
    for n in [1, 5, 8, 9, L, pts.len()] {
        let mut got = Vec::new();
        eval_points(&tape, &ctx, &pts[..n], &mut got);
        for (p, v) in pts[..n].iter().zip(got) {
            assert_eq!(v.to_bits(), eval_point(&tape, &ctx, *p).to_bits(), "at {p:?}");
            assert_eq!(v.to_bits(), values(&tape, &ctx, *p).last().unwrap().to_bits(), "at {p:?}");
        }
    }
}

/// `undefined_cases` one point at a time: the walk it batches, as it was
/// written before it was.
fn undefined_one_by_one(tape: &Tape, ctx: &Ctx, pts: &[V3]) -> Vec<(&'static str, usize, Option<u32>, usize, [u64; 3])> {
    let n = tape.ssa.len();
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
    let mut found: Vec<(&'static str, usize, Option<u32>, usize, V3)> = Vec::new();
    for &p in pts {
        let v = values(tape, ctx, p);
        let mut active = vec![false; n];
        let mut bad = Vec::new();
        active[n - 1] = true;
        for i in (0..n).rev() {
            if !active[i] {
                continue;
            }
            let ins = &tape.ssa[i];
            let (a, b, c) = (ins.a as usize, ins.b as usize, ins.c as usize);
            let what = match ins.op {
                Op::Bin(op) => crate::graph::undefined_bin(op, v[a], v[b]),
                Op::Un(op) => crate::graph::undefined_un(op, v[a]),
                _ => None,
            };
            if let Some(what) = what {
                bad.push((i, what));
            }
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
        for (i, what) in bad.into_iter().rev() {
            match found.iter_mut().find(|u| u.1 == i && u.0 == what) {
                Some(u) => u.3 += 1,
                None => found.push((what, i, owner[i], 1, p)),
            }
        }
    }
    found.sort_by_key(|u| std::cmp::Reverse(u.3));
    found.into_iter().map(|(w, i, o, c, p)| (w, i, o, c, p.map(f64::to_bits))).collect()
}

#[test]
fn batched_undefined_cases_are_each_points_own() {
    let (g, tape) = field();
    let ctx = Ctx::new(0.0, vec![]);
    let pts = points();
    let mut kinds = std::collections::HashSet::new();
    // In grid order, reversed and shuffled: which case a point meets first,
    // and so the order of cases whose counts tie, differs between them.
    let orders = [pts.clone(), pts.iter().rev().copied().collect(), (0..pts.len()).map(|k| pts[k * 7 % pts.len()]).collect::<Vec<_>>()];
    for pts in &orders {
        for n in [1, 2, 15, 16, 17, 64, 65, pts.len()] {
            let got: Vec<_> = undefined_cases(&tape, &ctx, &pts[..n]).into_iter().map(|u| (u.what, u.ins, u.name, u.count, u.at.map(f64::to_bits))).collect();
            assert_eq!(got, undefined_one_by_one(&tape, &ctx, &pts[..n]), "{n} points");
            kinds.extend(got.iter().map(|u| (u.0, u.2.map(|k| g.tag_names[k as usize].clone()))));
        }
    }
    // Each kind is met, and each of the select's guarded branches where it
    // is taken.
    for what in ["division by zero", "ln of a negative", "a negative base to a power that is not a whole number"] {
        assert!(kinds.iter().any(|k| k.0 == what), "{what}: {kinds:?}");
    }
    assert!(kinds.contains(&("division by zero", Some("guarded".to_string()))), "{kinds:?}");
    assert!(kinds.contains(&("ln of a negative", Some("else".to_string()))), "{kinds:?}");
}
