//! The first crossing of a field's zero set along a ray.

use crate::octree::Local;
use ink_core::V3;
use ink_core::eval::{eval_interval, eval_interval_in, eval_point, span_mid};
use ink_core::interval::{I, add_o, mul_o, sub_o};
use ink_core::sign::inside;

/// A side of a field's zero set: inside (the field < 0) or outside (>= 0).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Inside,
    Outside,
}

impl Side {
    pub fn holds(self, v: f64) -> bool {
        match self {
            Side::Inside => inside(v),
            Side::Outside => !inside(v),
        }
    }
}

/// Where a ray reached a side of a field (ray_to): at t, and, where
/// bounds rounded outward do not show a point within its eps after t on
/// that side, how far from 0 the field's real value may be there
/// (`rounding`: the ray runs in the surface, so only rounding decides the
/// side, and it may be on either).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayHit {
    pub t: f64,
    pub rounding: Option<f64>,
}

/// The first t in [t0, t1] where o + t dir is on `side` of the field.
/// Interval arithmetic over segments of the ray, nearest first, each with
/// the tape `local` keeps for the smallest block holding it, rules out the
/// stretches wholly on the other side, as does the field's value at a
/// segment's middle (bounded outward) and its rate along the ray (a mean
/// value bound, finite only where the field is continuous: no undecided
/// select, floor or atan2 across a jump); a segment neither rules out is
/// halved down to `eps`, where the values at its ends decide. So a stretch
/// of `side` longer than eps is never stepped over, and a touch that does
/// not cross is not a crossing. A segment whose middle's bound holds 0 and
/// whose rate times its half-length is under that bound's width is a ray
/// lying in the surface: halving it decides nothing, so it ends the search,
/// a hit where it crosses if bounds show that one sure, else from its
/// start, within rounding of 0 there (RayHit::rounding). The crossing found
/// is bisected to the precision of floats: t is on `side`, in floats, and
/// a hit is sure where bounds show a point of it within eps after t.
/// None: no point of the ray is on that side, but for stretches shorter
/// than eps.
pub fn ray_to(local: &Local, o: V3, dir: V3, range: (f64, f64), side: Side, eps: f64) -> Option<RayHit> {
    ray_to_within(local, o, dir, range, side, eps, f64::INFINITY)
}

/// A ray's evaluations, counted on its Local when it ends
/// (Local::ray_evals).
struct Evals<'l, 'a>(&'l Local<'a>, std::cell::Cell<usize>);

impl Evals<'_, '_> {
    fn add(&self, n: usize) {
        self.1.set(self.1.get() + n);
    }
}

impl Drop for Evals<'_, '_> {
    fn drop(&mut self) {
        self.0.count_ray_evals(self.1.get());
    }
}

/// ray_to's answer where it is at most `stop`; else None, or a hit past
/// `stop`, searched no further than the segment holding `stop`. Segments
/// are searched in order along the ray, so those before `stop` are the
/// ones ray_to searches: a hit at most `stop` is the one it finds.
pub fn ray_to_within(local: &Local, o: V3, dir: V3, (t0, t1): (f64, f64), side: Side, eps: f64, stop: f64) -> Option<RayHit> {
    let ctx = local.ctx();
    let evals = Evals(local, std::cell::Cell::new(0));
    let at = |t: f64| [o[0] + t * dir[0], o[1] + t * dir[1], o[2] + t * dir[2]];
    let f = |t: f64| {
        evals.add(1);
        let p = at(t);
        eval_point(&local.tape_over(&p.map(I::point)), ctx, p)
    };
    let eps = eps.max(1e-9 * t1.abs());
    // Bounds `r` show: no point on `side`, every point on it.
    let decide = |r: I| match side {
        Side::Inside => (r.none_inside(), r.all_inside()),
        Side::Outside => (r.all_inside(), r.none_inside()),
    };
    // The bounds at t: of the ray's real point there, and its float.
    let point = |t: f64| {
        evals.add(1);
        let (p, b) = (at(t), segment_box(o, dir, t, t));
        let b = std::array::from_fn(|k| I::new(b[k].lo.min(p[k]), b[k].hi.max(p[k])));
        eval_interval(&local.tape_over(&b), ctx, b, None)
    };
    // t, sure where bounds show it, or a point up to eps after it and
    // within the ray, on `side` (2^-20 eps first).
    let hit = |t: f64| {
        let r = point(t);
        let sure = decide(r).1 || (0..=20).rev().map(|i| t + eps * 0.5f64.powi(i)).filter(|&s| s <= t1).any(|s| decide(point(s)).1);
        RayHit { t, rounding: (!sure).then(|| r.lo.abs().max(r.hi.abs())) }
    };
    let start = f(t0);
    if side.holds(start) {
        return Some(hit(t0));
    }
    // Segments from t0, each twice as long as the one before, so a ray that
    // crosses soon is searched in the small blocks near where it starts;
    // each depth-first, nearest half first. The first is as long as the
    // field's value at t0 (a distance field's gap to the surface; a field
    // scaled by k costs log2 k more segments), at least eps.
    let (mut from, mut len) = (t0, if start.is_finite() { start.abs().max(eps) } else { eps });
    let (mut stack, mut slots) = (Vec::new(), Vec::new());
    loop {
        let (a, b) = match stack.pop() {
            Some(s) => s,
            None if from < t1 => {
                let s = (from, (from + len).min(t1));
                (from, len) = (s.1, 2.0 * len);
                s
            }
            None => return None,
        };
        if a > stop {
            return None;
        }
        evals.add(1);
        let bx = segment_box(o, dir, a, b);
        let tape = local.tape_over(&bx);
        let r = eval_interval_in(&tape, ctx, bx, None, &mut slots);
        let (mut none, mut all) = decide(r);
        if !none && !all {
            // Its middle's bound and its rate along the ray: the real field
            // over the segment lies within their mean value bound.
            evals.add(2);
            let g = ink_core::eval::eval_interval_grad(&tape, ctx, bx);
            let rate = (0..3).fold(I::point(0.0), |s, k| add_o(s, mul_o(g.d[k], I::point(dir[k]))));
            let m = span_mid((a, b));
            let mid = eval_interval_in(&tape, ctx, segment_box(o, dir, m, m), None, &mut slots);
            if rate.lo.is_finite() && rate.hi.is_finite() && mid.lo.is_finite() && mid.hi.is_finite() {
                let e = add_o(mid, mul_o(rate, sub_o(I::new(a, b), I::point(m))));
                (none, all) = decide(e);
                let undecided = !(decide(mid).0 || decide(mid).1);
                if !none && !all && undecided && ink_core::interval::mag(rate) * (b - a) * 0.5 <= mid.hi - mid.lo {
                    // Where it crosses in floats, if it does, and sure there
                    // where bounds show it (`hit`); else from a, to rounding.
                    let (fa, fb) = (f(a), f(b));
                    let t = match (side.holds(fa), side.holds(fb)) {
                        (false, true) => crossing(&f, side, (a, fa), (b, fb)),
                        _ => a,
                    };
                    return Some(match hit(t) {
                        h @ RayHit { rounding: None, .. } => h,
                        _ => RayHit { t: a, rounding: Some(e.lo.abs().max(e.hi.abs())) },
                    });
                }
            }
        }
        if none {
            continue;
        }
        if b - a <= eps.max(b.abs() * 1e-12) || all {
            // On that side somewhere in [a, b]: find the first crossing.
            let (fa, fb) = (f(a), f(b));
            if side.holds(fa) {
                return Some(hit(a));
            }
            if !side.holds(fb) {
                // No crossing at the ends; refine this stretch further.
                if b - a > eps {
                    let m = span_mid((a, b));
                    stack.push((m, b));
                    stack.push((a, m));
                }
                continue;
            }
            return Some(hit(crossing(&f, side, (a, fa), (b, fb))));
        }
        let m = span_mid((a, b));
        stack.push((m, b));
        stack.push((a, m));
    }
}

/// The crossing between a (off `side`) and b (on it), to the precision of
/// floats: the end on `side`. False position, Illinois's way (an end kept
/// twice has its value halved), and a bisection after any step that leaves
/// more than half the bracket, so it takes at most about twice the steps of
/// bisection and, where the field is smooth, a handful.
pub(super) fn crossing(f: &dyn Fn(f64) -> f64, side: Side, (mut a, mut fa): (f64, f64), (mut b, mut fb): (f64, f64)) -> f64 {
    let (mut kept, mut halve) = (0i8, false);
    loop {
        let (w, mid) = (b - a, span_mid((a, b)));
        if mid <= a || mid >= b {
            break;
        }
        let m = if halve { mid } else { a + w * (fa / (fa - fb)) };
        let m = if m > a && m < b { m } else { mid };
        let fm = f(m);
        if side.holds(fm) {
            (b, fb) = (m, fm);
            if kept == -1 {
                fa *= 0.5;
            }
            kept = -1;
        } else {
            (a, fa) = (m, fm);
            if kept == 1 {
                fb *= 0.5;
            }
            kept = 1;
        }
        halve = b - a > 0.5 * w;
    }
    b
}

/// A box holding the points o + t dir for t in [t0, t1], rounded outward
/// (nearest ends can miss a sliver of the segment, and the crossing in it).
pub(super) fn segment_box(o: V3, dir: V3, t0: f64, t1: f64) -> [I; 3] {
    std::array::from_fn(|k| add_o(I::point(o[k]), mul_o(I::new(t0, t1), I::point(dir[k]))))
}
