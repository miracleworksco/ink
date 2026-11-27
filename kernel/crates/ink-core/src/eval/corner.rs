//! Point values and bounds for corner reach.

use super::context::Ctx;
use super::gradient::eval_interval_grad;
use super::point::{eval_grads, eval_point};
use super::range::eval_interval;
use crate::V3;
use crate::graph::{BinOp, UnOp, defined, div_f64};
use crate::interval::{self, I};
use crate::tape::Tape;

pub(super) fn corner_reach_point(ta: &Tape, tb: &Tape, ctx: &Ctx, r: f64, q: V3) -> f64 {
    let ga = eval_grads(ta, ctx, &[q])[0];
    let gb = eval_grads(tb, ctx, &[q])[0];
    let na = (ga[1] * ga[1] + ga[2] * ga[2] + ga[3] * ga[3]).sqrt();
    let nb = (gb[1] * gb[1] + gb[2] * gb[2] + gb[3] * gb[3]).sqrt();
    let ra = r - ga[0];
    let rb = r - gb[0];
    let qa = [q[0] + rb * div_f64(gb[1], nb), q[1] + rb * div_f64(gb[2], nb), q[2] + rb * div_f64(gb[3], nb)];
    let qb = [q[0] + ra * div_f64(ga[1], na), q[1] + ra * div_f64(ga[2], na), q[2] + ra * div_f64(ga[3], na)];
    let am = eval_point(ta, ctx, qa);
    let bm = eval_point(tb, ctx, qb);
    defined(r - am.max(bm))
}

fn prune_corner(t: &std::sync::Arc<Tape>, ctx: &Ctx, b: [I; 3]) -> std::sync::Arc<Tape> {
    if t.len() <= 32 {
        return t.clone();
    }
    let mut ch = vec![crate::tape::BOTH; t.len()];
    eval_interval(t, ctx, b, Some(&mut ch));
    t.narrowed(&ch)
}

/// Sound box enclosing `q + (r - g.v) * (g.d / |g.d|)` over `q` (where `g.d / 0 = 0`),
/// even when `|g.d|` contains `0` (at a crease/corner of a sub-shape).
/// Uses the first-order centered form around `mid(q)` (`v_mid = eval_interval(sub, mid)`)
/// so along a flat face `q[k]` and `-v * u[k]` cancel instead of doubling box width
/// (e.g. a flush face where `reach == 0` stays tight on large boxes).
fn corner_shifted_box(sub: &Tape, ctx: &Ctx, q: [I; 3], r: I, g: interval::G) -> Option<[I; 3]> {
    if !(g.v.lo.is_finite() && g.v.hi.is_finite() && g.d.iter().all(|c| c.lo.is_finite() && c.hi.is_finite())) {
        return None;
    }
    let s2 = interval::add_o(interval::add_o(interval::un(UnOp::Square, g.d[0]), interval::un(UnOp::Square, g.d[1])), interval::un(UnOp::Square, g.d[2]));
    let n = interval::un(UnOp::Sqrt, s2);
    let dr = interval::sub_o(r, g.v);
    let mid = [0, 1, 2].map(|k| (0.5 * (q[k].lo + q[k].hi)).clamp(q[k].lo, q[k].hi));
    let delta = [0, 1, 2].map(|k| interval::sub_o(q[k], I::point(mid[k])));
    let v_mid = eval_interval(sub, ctx, mid.map(I::point), None);
    let dr_mid = (v_mid.lo.is_finite() && v_mid.hi.is_finite()).then(|| interval::sub_o(r, v_mid));
    let mut out = [I::ENTIRE; 3];
    for k in 0..3 {
        // |div_f64(g.d[k], n)| <= 1 + 2e-15 for every finite gradient vector
        // (and 0 when n == 0 or n == inf), with the same sign as g.d[k].
        let mut u_lo: f64 = if g.d[k].lo >= 0.0 { 0.0 } else { -1.000000000001 };
        let mut u_hi: f64 = if g.d[k].hi <= 0.0 { 0.0 } else { 1.000000000001 };
        if n.lo > 0.0 && n.hi.is_finite() {
            let qn = interval::div_o(g.d[k], n);
            if qn.lo.is_finite() && qn.hi.is_finite() {
                u_lo = u_lo.max(qn.lo);
                u_hi = u_hi.min(qn.hi);
            }
        }
        if !(u_lo <= u_hi) {
            (u_lo, u_hi) = (-1.000000000001, 1.000000000001);
        }
        let u = I::new(u_lo, u_hi);
        let plain = interval::add_o(q[k], interval::mul_o(dr, u));
        let c = if let Some(drm) = dr_mid {
            let mut cen = interval::add_o(I::point(mid[k]), interval::mul_o(drm, u));
            for j in 0..3 {
                let du = interval::mul_o(g.d[j], u);
                let coeff = if j == k { interval::sub_o(I::point(1.0), du) } else { interval::un(UnOp::Neg, du) };
                cen = interval::add_o(cen, interval::mul_o(coeff, delta[j]));
            }
            let met = interval::meet(plain, cen);
            if met.lo <= met.hi { met } else { plain }
        } else {
            plain
        };
        if !(c.lo.is_finite() && c.hi.is_finite()) {
            return None;
        }
        out[k] = c;
    }
    Some(out)
}

pub(super) fn corner_reach_interval(ta: &std::sync::Arc<Tape>, tb: &std::sync::Arc<Tape>, ctx: &Ctx, r: I, q: [I; 3], guard: I) -> I {
    if guard.lo >= 0.0 {
        return I::point(0.0);
    }
    let reach = if q.iter().all(|c| c.lo.is_finite() && c.hi.is_finite()) && r.lo.is_finite() && r.hi.is_finite() {
        let ta_q = prune_corner(ta, ctx, q);
        let tb_q = prune_corner(tb, ctx, q);
        let ga = eval_interval_grad(&ta_q, ctx, q);
        let gb = eval_interval_grad(&tb_q, ctx, q);
        let qa = corner_shifted_box(&tb_q, ctx, q, r, gb);
        let qb = corner_shifted_box(&ta_q, ctx, q, r, ga);
        let am = qa.map_or(I::ENTIRE, |b| eval_interval(ta, ctx, b, None));
        let bm = qb.map_or(I::ENTIRE, |b| eval_interval(tb, ctx, b, None));
        interval::sub_o(r, interval::bin(BinOp::Max, am, bm).0)
    } else {
        I::ENTIRE
    };
    if guard.hi < 0.0 { reach } else { interval::hull(reach, I::point(0.0)) }
}
