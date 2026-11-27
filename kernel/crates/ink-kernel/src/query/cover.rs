//! Proving a box stays within the tolerance of the outside of a set.

use ink_core::V3;
use ink_core::eval::{Ctx, eval_point, eval_points};
use ink_core::interval::I;
use ink_core::sign::inside;
use ink_core::tape::Tape;

/// The tolerance's rule for a set (sign 1: an emptiness claim's solid, f
/// below 0; -1: a void count's air, f not below 0, `ink_core::sign`): a
/// box of it settles to the tolerance where every point of the box lies
/// within the tolerance of a point shown outside the set (`covered`), so
/// none of it is deeper (a set no thicker than twice the tolerance passes
/// as absent). `region`: the set's box, past which there is none of it;
/// `asked`, the region asked (all of space for both, for a void count,
/// whose air goes on past its shape's box).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Cover {
    pub tol: f64,
    pub region: [I; 3],
    pub asked: [I; 3],
    pub sign: f64,
    /// The space axes the field reads (Tape::space_read): along one it
    /// reads not, no region ends.
    pub read: [bool; 3],
    /// The space axes the box's own field reads not (`covered`): along
    /// them it is the same all through the box, so a point there is as far
    /// outside as one beside it, and a point is shown outside along the
    /// box's whole span of them.
    pub still: [bool; 3],
}

impl Cover {
    /// Box b grown by the tolerance, within the region asked.
    pub(crate) fn grown(&self, b: [I; 3]) -> [I; 3] {
        std::array::from_fn(|k| I::new((b[k].lo - self.tol).next_down().max(self.asked[k].lo), (b[k].hi + self.tol).next_up().min(self.asked[k].hi)))
    }

    /// Whether every point of box b (sides `side`) lies within the
    /// tolerance of a point shown outside the set (`out`: by bounds there
    /// rounded outward, over the whole scope, and the point), so none of the
    /// box's part of the set is deeper than it. Needs no slope: it holds for any
    /// field. The points may lie beyond the box, within the tolerance of it
    /// (`wide` gives the field there): a box inside a contact thinner than twice the
    /// tolerance is settled from outside it. `pts` are the box's samples,
    /// the centres of its 4 x 4 x 4 sub-boxes (`box_samples`), `v` the field
    /// there (`sub` holds only in the box; `out` is asked of a point, as a
    /// box of its still axes' spans, and the bounds there).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn covered<'w>(
        &self,
        sub: &Tape,
        wide: &dyn Fn() -> &'w Tape,
        out: &dyn Fn(I, &Ctx, [I; 3]) -> bool,
        ctx: &Ctx,
        b: [I; 3],
        side: V3,
        pts: &[V3],
        v: &[f64],
    ) -> bool {
        if self.tol == 0.0 {
            return false;
        }
        let (tol, sign) = (self.tol, self.sign);
        // Outside the set at a point (sign times the field there, `x`), by
        // sign::inside: the solid's outside is not inside, the air's is.
        let outside = |x: f64| if sign > 0.0 { !inside(x) } else { inside(-x) };
        let grown = self.grown(b);
        let (read, still) = (self.read, self.still);
        let within = |p: &V3, r: &[I; 3]| (0..3).all(|k| r[k].contains(p[k]));
        let far = |q: V3, c: V3, half: V3, b: [I; 3]| far(q, c, half, b, still);
        let side: V3 = std::array::from_fn(|k| if still[k] { 0.0 } else { side[k] });
        // Shown outside: in the box by its own tape; past the region asked
        // or the material's box (a proven bound), by them; else in the
        // grown box by `wide`; each along the box's span of its still axes
        // (no step moves a point along them).
        let past = |q: &V3| (0..3).any(|k| read[k] && !(self.asked[k].contains(q[k]) && self.region[k].contains(q[k])));
        let shown = |q: &V3, slots: &mut Vec<I>| -> bool {
            let t = if within(q, &b) {
                sub
            } else if past(q) {
                return true;
            } else if within(q, &grown) {
                wide()
            } else {
                return false;
            };
            let at: [I; 3] = std::array::from_fn(|k| if still[k] { b[k] } else { I::point(q[k]) });
            out(ink_core::eval::eval_interval_in(t, ctx, at, None, slots), ctx, at)
        };
        // A point shown outside within the tolerance of every point of the
        // box about c (field x there): c itself, else a step from c up the
        // field's slope to where it would reach 0 were it a distance, else
        // the farthest the box allows up and down the slope and along each
        // axis (the nearest outside may not be up the slope: the side of an
        // overlap whose field is the lower, scaled). The steps are guesses;
        // the check is not.
        let near = |c: V3, x: f64, half: V3, slots: &mut Vec<I>| -> Option<V3> {
            if outside(x) && far(c, c, half, b) <= tol && shown(&c, slots) {
                return Some(c);
            }
            // A step this long keeps every point of the box about c within
            // the tolerance of where it lands, with room for far's rounding
            // (along one axis read alone, the farthest is reach + its half).
            let reach = 0.99 * (tol - ink_core::interval::norm(std::array::from_fn(|k| if still[k] { 0.0 } else { half[k] })));
            if !(reach > 0.0) {
                return None;
            }
            let [_, gx, gy, gz] = ink_core::eval::eval_grads(sub, ctx, &[c])[0].map(|x| sign * x);
            // (A slope that is 0 or not finite gives no step up it: its
            // points are NaN, which no test below may take as near or
            // outside, so they are dropped here.)
            let n = (gx * gx + gy * gy + gz * gz).sqrt();
            let step = ((-x / n).max(0.0) * 1.01 + 0.01 * reach).min(reach);
            let up = [gx / n, gy / n, gz / n];
            let axis = |k: usize, s: f64| -> V3 { std::array::from_fn(|a| if a == k { s } else { 0.0 }) };
            let ways = [
                (up, step),
                (up, reach),
                (up.map(|u| -u), reach),
                (axis(0, 1.0), reach),
                (axis(0, -1.0), reach),
                (axis(1, 1.0), reach),
                (axis(1, -1.0), reach),
                (axis(2, 1.0), reach),
                (axis(2, -1.0), reach),
            ];
            ways.into_iter()
                .filter(|(u, _)| (0..3).all(|k| !still[k] || u[k] == 0.0))
                .map(|(u, s)| [c[0] + s * u[0], c[1] + s * u[1], c[2] + s * u[2]])
                .filter(|q| q.iter().all(|x| x.is_finite()))
                .find(|q| {
                    far(*q, c, half, b) <= tol && (past(q) || within(q, &grown) && outside(sign * eval_point(if within(q, &b) { sub } else { wide() }, ctx, *q))) && shown(q, slots)
                })
        };
        let mut slots = Vec::new();
        // First the fewest boxes: the n x n x n whose corners are no more
        // than the tolerance from their centres (one for a box under about
        // the tolerance across, 8 for one twice that).
        let n = (side.iter().fold(0.0f64, |m, &x| m.max(x)) * 3f64.sqrt() / (2.0 * tol)).ceil().max(1.0);
        if n < 4.0 {
            let n = n as usize;
            let part: V3 = side.map(|x| x / n as f64);
            let mut grid = Vec::with_capacity(n * n * n);
            for a in 0..n {
                for c in 0..n {
                    for d in 0..n {
                        grid.push([b[0].lo + part[0] * (d as f64 + 0.5), b[1].lo + part[1] * (c as f64 + 0.5), b[2].lo + part[2] * (a as f64 + 0.5)]);
                    }
                }
            }
            let mut g = Vec::new();
            eval_points(sub, ctx, &grid, &mut g);
            g.iter_mut().for_each(|x| *x *= sign);
            let half = part.map(|x| 0.5 * x);
            if grid.iter().zip(&g).all(|(q, &x)| near(*q, x, half, &mut slots).is_some()) {
                return true;
            }
        }
        // Else each sub-box by a point already shown outside where one is
        // near enough (one covers its neighbours), else a sample outside
        // (its own centre first, then the nearest the box's centre), else
        // a step from its centre. The deepest first: a box that cannot be
        // covered fails at once.
        let v: Vec<f64> = v.iter().map(|x| sign * x).collect();
        let centre: V3 = std::array::from_fn(|k| 0.5 * b[k].lo + 0.5 * b[k].hi);
        let mut order: Vec<usize> = (0..pts.len()).filter(|&k| outside(v[k])).collect();
        let d2 = |k: usize| (0..3).map(|a| (pts[k][a] - centre[a]).powi(2)).sum::<f64>();
        order.sort_by(|&i, &j| d2(i).total_cmp(&d2(j)));
        let half: V3 = side.map(|x| 0.125 * x);
        // Samples outside, confirmed (lazily: None not yet asked).
        let mut out: Vec<Option<bool>> = vec![None; pts.len()];
        let mut shown_pts: Vec<V3> = Vec::new();
        let mut deepest: Vec<usize> = (0..pts.len()).collect();
        deepest.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
        deepest.into_iter().all(|j| {
            if shown_pts.iter().any(|&q| far(q, pts[j], half, b) <= tol) {
                return true;
            }
            let sample = std::iter::once(j)
                .chain(order.iter().copied().filter(|&k| k != j))
                .find(|&k| outside(v[k]) && far(pts[k], pts[j], half, b) <= tol && *out[k].get_or_insert_with(|| shown(&pts[k], &mut slots)));
            match sample.map(|k| pts[k]).or_else(|| near(pts[j], v[j], half, &mut slots)) {
                Some(q) => {
                    shown_pts.push(q);
                    true
                }
                None => false,
            }
        })
    }
}

/// The farthest a point of the box about c (half-sides `half`) is from q,
/// rounded up, with room for the coordinates' rounding (sub-boxes tile box
/// b to within a few ulps of its coordinates), along the axes not `still`
/// (along which q stands for every point of the box's span).
pub(super) fn far(q: V3, c: V3, half: V3, b: [I; 3], still: [bool; 3]) -> f64 {
    ink_core::interval::norm(std::array::from_fn(|k| {
        if still[k] {
            return 0.0;
        }
        let ulps = 16.0 * f64::EPSILON * c[k].abs().max(q[k].abs()).max(b[k].lo.abs()).max(b[k].hi.abs());
        ((q[k] - c[k]).abs().next_up() + half[k] + ulps).next_up().next_up()
    }))
}
