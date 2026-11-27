//! Switch moments and interval bounds for correlated branch choices.

use super::context::{Ctx, Inputs, fixed, input_at, span_mid};
use super::fork::same_branches;
use super::point::eval_point;
use super::range::eval_interval;
use super::scope::{eval_interval_rates, eval_interval_scope};
use crate::graph::{BinOp, UnOp, Var};
use crate::interval::{self, I};
use crate::tape::{Op, Tape};

/// Where a select's condition changes sign over the box and the scope,
/// found by its value: at `at`, a float, of the one input `var` of the
/// scope it reads (time, or a ranged param that is not derived; anything
/// spatial it reads a point of the box), along which it is strictly
/// monotone throughout (its rate's bound keeps one sign: rising or
/// falling) and exactly 0 there (bounded at that value alone: [0, 0], so
/// 0 in real numbers). It is then negative on one side of `at` and
/// positive on the other, at every point of the box and the scope: two
/// conditions with one moment switch together, however they are written.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Moment {
    pub var: Var,
    pub at: f64,
    pub rising: bool,
}

/// The Moment condition `c` (an ssa index) switches at, over the box `b`
/// and ctx's scope; None where it is not one (it reads no input or more
/// than one, or space along a side of the box, is not shown monotone, or
/// has no float its bound at is [0, 0]).
pub fn moment(tape: &Tape, ctx: &Ctx, c: usize, b: [I; 3]) -> Option<Moment> {
    let sides = (0..3).filter(|&k| b[k].lo != b[k].hi).fold(0, |m, k| m | 1 << k);
    if fixed(tape, ctx, c) || tape.axes_read(c) & sides != 0 {
        return None;
    }
    let cone = tape.cone(c);
    let mut var = None;
    for ins in &cone.ssa {
        let Op::Var(v) = ins.op else { continue };
        let input = match v {
            Var::X | Var::Y | Var::Z => {
                let x = b[match v {
                    Var::X => 0,
                    Var::Y => 1,
                    _ => 2,
                }];
                if x.lo != x.hi {
                    return None;
                }
                continue;
            }
            Var::T => ctx.t_span.is_some_and(|(a, b)| a < b),
            Var::Param(k) => ctx.p_spans.get(k).is_some_and(|(a, b)| a < b && !ctx.p_spans.is_derived(k)),
        };
        if input && var.replace(v).is_some_and(|w| w != v) {
            return None;
        }
    }
    let var = var?;
    let span = match var {
        Var::T => ctx.t_span?,
        Var::Param(k) => ctx.p_spans.get(k)?,
        _ => return None,
    };
    let inputs = match var {
        Var::T => Inputs { t: true, params: Vec::new() },
        Var::Param(k) => Inputs { t: false, params: vec![k] },
        _ => return None,
    };
    let rate = eval_interval_rates(&cone, ctx, b, &inputs).1[0];
    // A slope with one sign on either side of a pole does not establish
    // monotonicity across it. The bound must cover the whole interval.
    if !rate.lo.is_finite() || !rate.hi.is_finite() {
        return None;
    }
    let rising = if rate.lo > 0.0 {
        true
    } else if rate.hi < 0.0 {
        false
    } else {
        return None;
    };
    // Halved by its value in floats to a pair of neighbouring floats (or
    // one where it is 0), then each float about them checked in bounds.
    let p = [b[0].lo, b[1].lo, b[2].lo];
    let f = |u: f64| eval_point(&cone, &input_at(ctx, var, u), p);
    let (mut lo, mut hi) = span;
    for _ in 0..2100 {
        let m = span_mid((lo, hi));
        if !(lo < m && m < hi) {
            break;
        }
        let v = f(m);
        if v.is_nan() {
            return None;
        }
        if v == 0.0 {
            (lo, hi) = (m, m);
            break;
        }
        if (v < 0.0) == rising {
            lo = m;
        } else {
            hi = m;
        }
    }
    let exact = |u: f64| {
        let v = eval_interval(&cone, &input_at(ctx, var, u), b, None);
        v.lo == 0.0 && v.hi == 0.0
    };
    let near = [lo.next_down(), lo, hi, hi.next_up()];
    let at = near.into_iter().find(|&u| span.0 <= u && u <= span.1 && exact(u))?;
    Some(Moment { var, at, rising })
}

/// What eval_interval_cases groups selects' conditions by: their form
/// (one value or its negation: the pair of a difference's operands, or
/// the value and usize::MAX), or the Moment they switch at.
#[derive(Clone, Copy, PartialEq)]
enum Switch {
    Form((usize, usize)),
    At(Var, u64),
}

/// Conditions eval_interval_cases splits on at once (at most 3^3 cases).
const MAX_CASE_CONDITIONS: usize = 3;

/// eval_interval_scope, met with the hull of its bounds in each case of
/// the selects it leaves undecided: those whose conditions are one value
/// or its negation, or switch at one Moment, take their sides together at
/// every point (c < 0, c = 0 or c > 0), so each case is bounded with the
/// branches it takes alone;
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
    let mut groups: Vec<(Switch, Vec<(usize, bool)>)> = Vec::new();
    for (i, ins) in ssa.iter().enumerate() {
        if !read[i] || !matches!(ins.op, Op::Select) || ch[i] != BOTH || fixed(tape, ctx, i) || same_branches(tape, i) {
            continue;
        }
        let (c, neg) = condition(ins.a as usize);
        match groups.iter_mut().find(|g| g.0 == Switch::Form(c)) {
            Some(g) => g.1.push((i, neg)),
            None => groups.push((Switch::Form(c), vec![(i, neg)])),
        }
    }
    // Conditions written apart that switch at one Moment are one group,
    // its sign that of the input less the moment (a condition falling
    // through it reads it negated).
    if groups.len() > 1 {
        let mut by_value: Vec<(Switch, Vec<(usize, bool)>)> = Vec::new();
        for (key, sel) in groups {
            let (s, flip) = match moment(tape, ctx, ssa[sel[0].0].a as usize, b) {
                Some(m) => (Switch::At(m.var, m.at.to_bits()), m.rising == sel[0].1),
                None => (key, false),
            };
            let sel = sel.into_iter().map(|(i, neg)| (i, neg ^ flip));
            match by_value.iter_mut().find(|g| g.0 == s) {
                Some(g) => g.1.extend(sel),
                None => by_value.push((s, sel.collect())),
            }
        }
        groups = by_value;
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
    let mut floors: Vec<(usize, [f64; 2], Vec<usize>)> = Vec::new();
    for (i, ins) in ssa.iter().enumerate() {
        if !read[i] || !matches!(ins.op, Op::Un(UnOp::Floor)) || fixed(tape, ctx, i) {
            continue;
        }
        if let Some(f) = floors.iter_mut().find(|f| f.0 == ins.a as usize) {
            f.2.push(i);
            continue;
        }
        let a = eval_interval_scope(&tape.cone(ins.a as usize), ctx, b, None);
        let (lo, hi) = (a.lo.floor(), a.hi.floor());
        if hi - lo == 1.0 {
            floors.push((ins.a as usize, [lo, hi], vec![i]));
        }
    }
    for (_, ends, is) in &floors {
        cases.push(ends.map(|k| (Vec::new(), is.iter().map(|&i| (i, k)).collect())).to_vec());
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
