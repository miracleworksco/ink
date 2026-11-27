//! A step's loads in time: where they jump, and their integrals over a step.

use crate::eval::Ctx;
use crate::tape::Tape;
use ink_core::{eval, interval};

/// The moments in (0, duration] where one of `tapes`, functions of t
/// (anywhere in the box `over`), may jump (its rate over a span unbounded:
/// a select or floor on t switching), ascending: each isolated by
/// bisection to a step between adjacent floats, such steps that touch
/// taken as one switch (rounding widens a bound by a float or two around
/// it), and given as that span, `lo` before the jump and `hi` after it
/// (one at the duration, too: the step ending there is the one before
/// it). None past `most` of them, or past the evaluations isolating that
/// many takes (two a halving, a halving for each float between the span's
/// ends at most): a stretch whose rate has no bound anywhere is no jump.
pub(crate) fn jumps_in_t(tapes: &[&Tape], ctx: &Ctx, over: [interval::I; 3], duration: f64, most: usize) -> Option<Vec<(f64, f64)>> {
    let c = Ctx { p_spans: eval::ParamSpans::NONE, ..ctx.clone() };
    let inputs = eval::Inputs { t: true, params: Vec::new() };
    const HALVINGS: usize = (f64::MANTISSA_DIGITS as i32 + f64::MAX_EXP - f64::MIN_EXP) as usize;
    let mut found: Vec<(f64, f64)> = Vec::new();
    for tape in tapes {
        let mut budget = 2 * HALVINGS * (most + 1);
        // Its steps come in order (the left half searched first), so each
        // joins the switch before it where they touch.
        let mut mine: Vec<(f64, f64)> = Vec::new();
        let mut open = vec![(0.0, duration)];
        while let Some(sp) = open.pop() {
            budget = budget.checked_sub(1)?;
            let span = Ctx { t: eval::span_mid(sp), t_span: Some(sp), ..c.clone() };
            let (_, rate, jump) = eval::eval_interval_rates_jump(tape, &span, over, &inputs);
            if jump <= 0.0 && rate.first().is_none_or(|r| r.lo.is_finite() && r.hi.is_finite()) {
                continue;
            }
            let m = eval::span_mid(sp);
            if sp.0 < m && m < sp.1 {
                open.extend([(m, sp.1), (sp.0, m)]);
                continue;
            }
            match mine.last_mut() {
                Some(last) if sp.0 <= last.1 => last.1 = last.1.max(sp.1),
                _ => mine.push(sp),
            }
            if mine.len() > most {
                return None;
            }
        }
        found.extend(mine);
    }
    found.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut merged: Vec<(f64, f64)> = Vec::new();
    for sp in found {
        match merged.last_mut() {
            Some(last) if sp.0 <= last.1 => last.1 = last.1.max(sp.1),
            _ => merged.push(sp),
        }
    }
    let out: Vec<(f64, f64)> = merged.into_iter().filter(|&sp| eval::span_mid(sp) > 0.0).collect();
    (out.len() <= most).then_some(out)
}

/// How closely `integral_in_t` bounds a step's integral: the spread of
/// its bounds under this share of their magnitude.
pub(crate) const INTEGRAL_TOL: f64 = 1e-3;
/// The pieces `integral_in_t` splits a step into at most.
pub(crate) const MAX_PIECES: usize = 1 << 14;

/// The integral of `tape`, a function of t alone (a transient load's
/// total), from a to b, and the spread of its bounds as a share of their
/// magnitude. Each piece of [a, b] is bounded by interval evaluation over
/// its times (`Ctx::t_span`), so no pulse between samples is missed: its
/// bound covers all of it. The piece whose bound is widest is halved until
/// the bounds' spread is under INTEGRAL_TOL of their magnitude (or there
/// are MAX_PIECES), and each is read at its middle, within its bound. Err
/// where a bound is not finite.
pub(crate) fn integral_in_t(tape: &Tape, ctx: &Ctx, a: f64, b: f64) -> Result<(f64, f64), String> {
    use std::collections::BinaryHeap;
    // A piece: (its bound's width times its length, a, b, bound).
    struct Piece(f64, f64, f64, interval::I);
    impl PartialEq for Piece {
        fn eq(&self, o: &Self) -> bool {
            self.0.total_cmp(&o.0).is_eq()
        }
    }
    impl Eq for Piece {}
    impl PartialOrd for Piece {
        fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(o))
        }
    }
    impl Ord for Piece {
        fn cmp(&self, o: &Self) -> std::cmp::Ordering {
            self.0.total_cmp(&o.0)
        }
    }
    // At the params' values: a toleranced param's span would keep every
    // piece's bound wide.
    let (mut scratch, mut c) = (Vec::new(), Ctx { p_spans: eval::ParamSpans::NONE, ..ctx.clone() });
    let mut piece = |a: f64, b: f64| -> Result<Piece, String> {
        (c.t, c.t_span) = (0.5 * (a + b), Some((a, b)));
        let i = eval::eval_interval_in(tape, &c, [interval::I::point(0.0); 3], None, &mut scratch);
        if !(i.lo.is_finite() && i.hi.is_finite() && i.lo <= i.hi) {
            return Err(format!("its value is not bounded over t from {a} to {b}"));
        }
        Ok(Piece((i.hi - i.lo) * (b - a), a, b, i))
    };
    let size = |p: &Piece| p.3.lo.abs().max(p.3.hi.abs()) * (p.2 - p.1);
    let first = piece(a, b)?;
    let (mut width, mut mag) = (first.0, size(&first));
    let mut heap = BinaryHeap::from([first]);
    while width > INTEGRAL_TOL * mag && heap.len() < MAX_PIECES {
        let p = heap.pop().expect("a piece");
        let m = 0.5 * (p.1 + p.2);
        if !(p.1 < m && m < p.2) {
            heap.push(p);
            break;
        }
        let (l, r) = (piece(p.1, m)?, piece(m, p.2)?);
        width += l.0 + r.0 - p.0;
        mag += size(&l) + size(&r) - size(&p);
        heap.push(l);
        heap.push(r);
    }
    // Summed afresh: the running sums drift.
    let (width, mag): (f64, f64) = (heap.iter().map(|p| p.0).sum(), heap.iter().map(size).sum());
    // A value that passes through an undefined case at a moment read (the language takes it as 0: a
    // logarithm of a negative, a division by zero), or is infinite there, has no answer: refused,
    // saying what and when. (A square root of a negative is defined: the root of 0.)
    let mut c = ctx.clone();
    let mut value = 0.0;
    for p in heap.iter() {
        for t in [p.1, 0.5 * (p.1 + p.2), p.2] {
            c.t = t;
            let v = eval::eval_point(tape, &c, [0.0; 3]);
            if let Some(u) = eval::undefined_cases(tape, &c, &[[0.0; 3]]).first() {
                return Err(format!("its value passes through {} at t = {t}, which has no value", u.what));
            }
            if !v.is_finite() {
                return Err(format!("its value is infinite at t = {t}"));
            }
            if t == 0.5 * (p.1 + p.2) {
                value += v.clamp(p.3.lo, p.3.hi) * (p.2 - p.1);
            }
        }
    }
    Ok((value, if mag > 0.0 { width / mag } else { 0.0 }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink_core::graph;

    #[test]
    fn a_steps_heat_is_its_sources_integral_however_short_the_pulse() {
        // A 10 J pulse 0.02 s wide at t = 2.5, inside a 6 s step: read at
        // the step's end (or its middle) it put no heat in.
        let mut g = graph::Graph::new();
        let (t, c) = (g.var(graph::Var::T), g.constant(2.5));
        let (a, b) = (g.sub(t, c), g.sub(c, t));
        let d = g.max(a, b);
        let (one, k) = (g.constant(1.0), g.constant(100.0));
        let kd = g.mul(k, d);
        let tri = g.sub(one, kd);
        let zero = g.constant(0.0);
        let pos = g.max(zero, tri);
        let w = g.constant(1000.0);
        let p = g.mul(w, pos);
        let tape = Tape::new(&g, p);
        let ctx = Ctx::new(0.0, Vec::new());
        let (e, spread) = integral_in_t(&tape, &ctx, 0.0, 6.0).unwrap();
        assert!((e - 10.0).abs() <= 10.0 * INTEGRAL_TOL && spread <= INTEGRAL_TOL, "{e} J, spread {spread}");
        // A step past it takes none, exactly.
        assert_eq!(integral_in_t(&tape, &ctx, 6.0, 12.0).unwrap(), (0.0, 0.0));
        // Unbounded power is refused, not sampled.
        let inv = g.bin(graph::BinOp::Div, one, a);
        assert!(integral_in_t(&Tape::new(&g, inv), &ctx, 0.0, 6.0).is_err());
    }
}
