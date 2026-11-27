//! Clearance and minimum distance between solids, with witness pairs.

use super::bounds::{focus, max_extent};
use super::empty::prove_empty;
use super::proof::{Confirm, Ended, Verdict, at_moment, box_samples, confirm_inside, contract_to, in_span, scope_shares, share, signed};
use super::ray::{Side, crossing};
use crate::scope::{self, At};
use crate::search;
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_interval, eval_interval_scope, eval_point, eval_points};
use ink_core::interval::I;
use ink_core::sign::inside;
use ink_core::tape::{BOTH, Tape};
use std::sync::Arc;

/// An upper bound on the distance between two points (rounded up). A
/// difference of 0 is exact: rounded up, the distance between a point and
/// itself would be 1.5e-323, reported as an overlapping clearance's value.
pub fn distance_up(p: V3, q: V3) -> f64 {
    ink_core::interval::norm(std::array::from_fn(|k| match (p[k] - q[k]).abs() {
        0.0 => 0.0,
        d => d.next_up(),
    }))
}

/// Newton steps onto a field's zero set from p, until a step is shorter
/// than `tol` (a length: the field's own scale does not matter). None: the
/// steps did not settle (no gradient, or a field that is not a distance
/// sends them back and forth across a crease).
pub fn to_surface(tape: &Tape, ctx: &Ctx, mut p: V3, tol: f64) -> Option<V3> {
    for _ in 0..30 {
        let g = eval_grads(tape, ctx, &[p])[0];
        let n2 = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
        if g[0] == 0.0 {
            return Some(p);
        }
        let k = g[0] / n2;
        if !(n2 > 0.0 && k.is_finite()) {
            return None;
        }
        p = [p[0] - k * g[1], p[1] - k * g[2], p[2] - k * g[3]];
        if !(k.abs() * n2.sqrt() > tol) {
            return p.iter().all(|x| x.is_finite()).then_some(p);
        }
    }
    None
}

/// A pair of points confirmed inside A and B, for a distance search that
/// found none (`min_distance` walks from B along A's slope, which a field
/// flat away from its solid, a union with a constant, does not have): a
/// point proven inside each (the emptiness search's witness, `finest` its
/// cells, over the boxes `ab` and `bb`), the segment between them cut where
/// it enters each solid, then moved closer (`closer_pair`). Its distance,
/// rounded up, is an upper bound on the solids' minimum distance: the
/// search goes on from it. None: a solid has no point found inside.
pub fn inside_pair(ta: &Arc<Tape>, tb: &Arc<Tape>, ctx: &Ctx, ab: (V3, V3), bb: (V3, V3), finest: f64) -> Option<(f64, V3, V3)> {
    let scale = max_extent(ab.0, ab.1).max(max_extent(bb.0, bb.1));
    let point_in = |t: &Arc<Tape>, (lo, hi): (V3, V3)| match prove_empty(t, ctx, lo, hi, None, finest, 0.0).verdict {
        Verdict::Disproven { witness, .. } => confirm_inside(t, ctx, witness, scale),
        _ => None,
    };
    let (p, q) = (point_in(ta, ab)?, point_in(tb, bb)?);
    // Where the segment from x (outside t's solid) to y (inside) first
    // enters it, confirmed; else y.
    let enter = |t: &Tape, x: V3, y: V3| {
        let at = |s: f64| std::array::from_fn(|k| x[k] + s * (y[k] - x[k]));
        let f = |s: f64| eval_point(t, ctx, at(s));
        match f(0.0) {
            fx if !inside(fx) => confirm_inside(t, ctx, at(crossing(&f, Side::Inside, (0.0, fx), (1.0, f(1.0)))), scale).unwrap_or(y),
            _ => y,
        }
    };
    let p = enter(ta, q, p);
    let q = enter(tb, p, q);
    // Each in turn slid along its surface toward the other (`slide`),
    // while that brings them closer.
    let (mut p, mut q) = (p, q);
    let mut dist = distance_up(p, q);
    for _ in 0..64 {
        let p2 = slide(ta, ctx, p, q, scale);
        let q2 = slide(tb, ctx, q, p2, scale);
        let d2 = distance_up(p2, q2);
        if !(d2 < dist) {
            break;
        }
        (dist, p, q) = (d2, p2, q2);
    }
    let (d2, p2, q2) = closer_pair(ta, tb, ctx, p, q, scale);
    Some(if d2 < dist { (d2, p2, q2) } else { (dist, p, q) })
}

/// p, confirmed inside t's solid at its surface, moved along the surface
/// toward x: the part of x - p across t's slope at p, a step back onto the
/// surface (to_surface) and confirmed inside, halved until it comes closer
/// to x; p where none does. Only t's slope at its own surface is read, so
/// a field flat away from its solid moves it.
fn slide(t: &Tape, ctx: &Ctx, p: V3, x: V3, scale: f64) -> V3 {
    let g = eval_grads(t, ctx, &[p])[0];
    let n2 = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
    if !(n2 > 0.0 && n2.is_finite()) {
        return p;
    }
    let to: V3 = std::array::from_fn(|k| x[k] - p[k]);
    let along = (to[0] * g[1] + to[1] * g[2] + to[2] * g[3]) / n2;
    let across: V3 = std::array::from_fn(|k| to[k] - along * g[k + 1]);
    let gap = distance_up(p, x);
    let mut h = 1.0;
    for _ in 0..20 {
        let y = std::array::from_fn(|k| p[k] + h * across[k]);
        if let Some(y) = to_surface(t, ctx, y, 1e-12 * scale).and_then(|y| confirm_inside(t, ctx, y, scale))
            && distance_up(y, x) < gap
        {
            return y;
        }
        h *= 0.5;
    }
    p
}

/// A pair of points confirmed inside A and B (p, q) moved closer by
/// alternating projections (q onto A, then that onto B), a move kept only
/// while both points stay confirmed inside (`confirm_inside`, `scale` a
/// length) and the pair comes closer: the pair and its distance rounded
/// up, an upper bound on the solids' minimum distance.
pub fn closer_pair(ta: &Tape, tb: &Tape, ctx: &Ctx, mut p: V3, mut q: V3, scale: f64) -> (f64, V3, V3) {
    let tol = 1e-12 * scale;
    let mut dist = distance_up(p, q);
    for _ in 0..20 {
        let Some(p2) = to_surface(ta, ctx, q, tol).and_then(|x| confirm_inside(ta, ctx, x, scale)) else { break };
        let Some(q2) = to_surface(tb, ctx, p2, tol).and_then(|x| confirm_inside(tb, ctx, x, scale)) else { break };
        let d2 = distance_up(p2, q2);
        if !(d2 < dist) {
            break;
        }
        (dist, p, q) = (d2, p2, q2);
    }
    (dist, p, q)
}

/// Samples a box seeks a witness pair from (those nearest both shapes first).
const WITNESS_STARTS: usize = 8;

/// Every point within d of the box, rounded outward (rounded to nearest,
/// lo - d may land above the true value and drop points within d: a bound
/// over it then certifies nothing).
fn grow(bx: [I; 3], d: f64) -> [I; 3] {
    bx.map(|x| I::new((x.lo - d).next_down(), (x.hi + d).next_up()))
}

/// A box's centre (halved first, so it cannot overflow; on an unbounded
/// side a finite end, else 0) and the distance from it to the box's far
/// side on each axis, rounded up: the centre may round to one end of a
/// narrow side.
fn centre_reach(bx: [I; 3]) -> (V3, V3) {
    let c: V3 = bx.map(|x| match x.mid() {
        m if m.is_finite() => m,
        _ if x.lo.is_finite() => x.lo,
        _ if x.hi.is_finite() => x.hi,
        _ => 0.0,
    });
    (c, std::array::from_fn(|k| (c[k] - bx[k].lo).max(bx[k].hi - c[k]).next_up()))
}

/// The share of a field's slope bound along axis k (0 to 1; 1 where the
/// slope is unbounded): times the box's half-width there, how far along k
/// it spreads the field's bounds, as a length.
fn slope_share(g: &ink_core::interval::G, l: f64, k: usize) -> f64 {
    let m = ink_core::interval::mag(g.d[k]);
    if !(m > 0.0) {
        0.0
    } else if !(m.is_finite() && l.is_finite() && l > 0.0) {
        1.0
    } else {
        m / l
    }
}

/// A clearance d as a question of the one search: no point of B (the
/// covered shape, `b`) in a box within d of A. `a` is A's field near the
/// region, where witness pairs are sought with `b`; each box keeps A
/// simplified over it grown by d and B over it.
///
/// A box is settled where any test passes, each sound for any fields:
/// - b >= 0 over it (no point of B), by its bounds;
/// - a > 0 over it grown by d, by its bounds;
/// - a >= d L over it, L bounding |grad a| over it grown by d, or the
///   mean-value form a(centre) >= L (r + d), r its half-diagonal;
/// - max(b, a - L d) >= 0 by its first-order bound, whose pairs see B's
///   surface and A's facing each other across the gap, so the boxes need
///   only be about sqrt(margin x radius), where the tests of a alone need
///   boxes as small as the margin;
/// - over the scope's axes (ranged params, time over a step of a motion),
///   a's or b's bound first order in them;
/// - where the box is flat (`Flat`: the shapes there are sets in fewer
///   axes times spans along the rest, an extrusion's), the gap between
///   their spans at least d, else the tests above of the sets, asked for
///   the clearance the gap leaves.
///
/// Open, it proves every point of B in it at least a / L from A (its
/// bound, below d; flat, that and the gap in quadrature), and its smear on each axis is a's share of its slope
/// bound along it times the half-width: the spread a's bounds take along
/// it, as a length (B's slope would cost a bound per box for no fewer
/// boxes). A ranged param's, and time's, is its blur of a's bound over L,
/// or of b's over b's slope bound where more: a step's moments are halved
/// where the motion blurs a box more than its size does (a face sliding
/// along itself, little). A witness pair is at the box's moment.
struct Clear<'a> {
    a: &'a Arc<Tape>,
    b: &'a Arc<Tape>,
    d: f64,
    /// The region the boxes cover, and its size (the length witness points
    /// are confirmed at).
    region: [I; 3],
    scale: f64,
}

/// A point of B, a point of A within d of it, and their moment and ranged
/// params (the same for both).
type Near = (V3, V3, At);

/// What a clearance box keeps for its halves: A and B simplified over it
/// (A over it grown by d), where witness pairs are sought; and where it
/// is flat, what the tests ask instead.
type Pair = (Arc<Tape>, Arc<Tape>, Option<Flat>);

/// How finely a flat box's spans are found: their ends to 2^-SPAN_DEPTH of
/// its side along them (outward), which no halving refines.
const SPAN_DEPTH: u32 = 40;

/// A clearance box as a question in fewer axes. Where each shape's field
/// where it is read (B over the box, A over it grown by d), for its sign
/// alone (`signed`), is a max of terms each reading axis k alone or not
/// reading it, for each such k (an extrusion's span; a cuboid's sides),
/// the shape there is the set of the rest of its terms times its span
/// along those axes. Two such points are then at least sqrt(e^2 + g^2)
/// apart, e the rests' distance in the other axes and g the gap between
/// the spans (each found by `contract_to`, outward). So the box is asked
/// of the rests, for a clearance of sqrt(d^2 - g^2) (settled where g >= d),
/// and never halved along those axes, which the rests read not; its halves
/// keep it, their rests the box's simplified, their spans the box's.
/// Not taken where the spans' terms read the scope (a span over ranged
/// params narrows as they are halved), nor along an axis whose taking
/// leaves a shape no rest (a cuboid is taken along all but one).
#[derive(Clone)]
struct Flat {
    a: Arc<Tape>,
    b: Arc<Tape>,
    still: [bool; 3],
    /// The gaps between the spans, squared and summed, rounded down
    /// (infinite: a shape has no point there).
    gap2: f64,
}

impl Flat {
    /// Box bx flat, A (`sa`) read over `grown`, B over bx; None where it
    /// is not.
    fn of(sa: &Arc<Tape>, sb: &Arc<Tape>, ctx: &Ctx, bx: [I; 3], grown: [I; 3]) -> Option<Flat> {
        let (ta, tb) = (signed(sa, ctx, grown), signed(sb, ctx, bx));
        let (ca, cb) = (ta.conjuncts(), tb.conjuncts());
        let alone = |cs: &[(usize, bool, u8)], k: usize| cs.iter().all(|c| c.2 & 1 << k == 0 || c.2 == 1 << k);
        let mut still: [bool; 3] = std::array::from_fn(|k| alone(&ca, k) && alone(&cb, k));
        // Each shape's terms: (its rest's, its spans').
        let split = |cs: &[(usize, bool, u8)], still: [bool; 3]| -> (Vec<(usize, bool)>, Vec<(usize, bool)>) {
            let (span, rest): (Vec<&(usize, bool, u8)>, Vec<_>) = cs.iter().partition(|c| c.2.count_ones() == 1 && still[c.2.trailing_zeros() as usize]);
            (rest.iter().map(|c| (c.0, c.1)).collect(), span.iter().map(|c| (c.0, c.1)).collect())
        };
        let ((ra, pa), (rb, pb)) = loop {
            let (a, b) = (split(&ca, still), split(&cb, still));
            if !a.0.is_empty() && !b.0.is_empty() {
                break (a, b);
            }
            still[still.iter().position(|&s| s)?] = false;
        };
        if !still.contains(&true) {
            return None;
        }
        let scoped = |t: &Tape, terms: &[(usize, bool)]| {
            !terms.is_empty() && {
                let t = t.cone_terms(terms);
                !scope::read_params(&t, ctx).is_empty() || scope::spans_t(&t, ctx)
            }
        };
        if scoped(&ta, &pa) || scoped(&tb, &pb) {
            return None;
        }
        let depths = still.map(|s| if s { SPAN_DEPTH } else { 0 });
        let gap2 = match (contract_to(&ta, ctx, grown, depths), contract_to(&tb, ctx, bx, depths)) {
            (Some(na), Some(nb)) => (0..3).filter(|&k| still[k]).fold(0.0, |sum: f64, k| {
                let g = (nb[k].lo - na[k].hi).max(na[k].lo - nb[k].hi).next_down();
                if g > 0.0 { (sum + (g * g).next_down()).next_down() } else { sum }
            }),
            _ => f64::INFINITY,
        };
        Some(Flat { a: Arc::new(ta.cone_terms(&ra)), b: Arc::new(tb.cone_terms(&rb)), still, gap2 })
    }
}

impl search::Question for Clear<'_> {
    type State = Pair;
    type Witness = Near;

    fn test(&self, n: &search::Node<Self::State>) -> search::Test<Self::State, Near> {
        let (bx, ctx) = (n.b, &*n.ctx);
        match &n.state {
            (fa, fb, Some(f)) => self.ask(bx, ctx, (fa, fb), (&f.a, &f.b), Some(f)),
            // The shapes' own tests first; a box they leave open is asked
            // again of its rests where it is flat.
            (fa, fb, None) => match self.near(bx, ctx, self.d, fa, fb) {
                None => search::Test::Settled,
                Some((sa, sb)) => match self.ask(bx, ctx, (&sa, &sb), (&sa, &sb), None) {
                    open @ search::Test::Open { .. } => match Flat::of(&sa, &sb, ctx, bx, grow(bx, self.d)) {
                        Some(f) => self.ask(bx, ctx, (&sa, &sb), (&f.a, &f.b), Some(&f)),
                        None => open,
                    },
                    t => t,
                },
            },
        }
    }

    /// From samples of B in the box (its 4 x 4 x 4), nearest A first,
    /// Newton steps along grad a toward A's surface, keeping within d;
    /// else, from the samples nearest both shapes, a point of B at or
    /// under its surface and one of A beside it, moved closer by
    /// alternating projections (closer_pair): a gap under d that no sample
    /// lies within d of. Over ranged params, the same with each param
    /// pushed to the end of its span that lowers both fields at the sample
    /// (scope::pushed).
    fn seek(&self, n: &search::Node<Self::State>) -> Option<Near> {
        let (bx, ctx) = (n.b, &*n.ctx);
        let (sa, sb, _) = &n.state;
        if !bx.iter().all(|x| x.lo.is_finite() && x.hi.is_finite()) {
            return None;
        }
        let pts = box_samples(bx, 4);
        let (mut vb, mut va) = (Vec::new(), Vec::new());
        eval_points(sb, ctx, &pts, &mut vb);
        eval_points(sa, ctx, &pts, &mut va);
        let mut order: Vec<usize> = (0..pts.len()).filter(|&k| pts[k].iter().all(|x| x.is_finite()) && (0..3).all(|j| bx[j].contains(pts[k][j]))).collect();
        // A sample where a >= L d has no point of A within d (L bounds
        // |grad a| on the box grown by d): no walk starts there.
        let reach = ink_core::interval::mul_up(ink_core::interval::lipschitz(&ink_core::eval::eval_interval_grad(sa, ctx, grow(bx, self.d))), self.d);
        order.sort_by(|&i, &j| va[i].total_cmp(&va[j]));
        let walks = order.iter().filter(|&&k| inside(vb[k]) && !(va[k] >= reach)).take(WITNESS_STARTS);
        if let Some(w) = walks.copied().find_map(|k| self.walk(sa, ctx, pts[k])) {
            return Some(w);
        }
        order.sort_by(|&i, &j| va[i].max(vb[i]).total_cmp(&va[j].max(vb[j])));
        if let Some(w) = order.iter().take(WITNESS_STARTS).find_map(|&k| self.pair((sa, sb), ctx, pts[k])) {
            return Some(w);
        }
        if scope::read_params(sa, ctx).is_empty() && scope::read_params(sb, ctx).is_empty() {
            return None;
        }
        order.into_iter().take(4).find_map(|k| {
            let pushed = scope::pushed(&[sa, sb], ctx, pts[k])?.0;
            self.walk(sa, &pushed, pts[k]).or_else(|| self.pair((sa, sb), &pushed, pts[k]))
        })
    }
}

impl Clear<'_> {
    /// A and B over box bx (A over it grown by d), each simplified there;
    /// None where either has no point there (a field never negative within
    /// d of the box has no material there: its zero set has no interior, a
    /// polygon's hollow turned about an edge on the axis).
    fn near(&self, bx: [I; 3], ctx: &Ctx, d: f64, sa: &Arc<Tape>, sb: &Arc<Tape>) -> Option<(Arc<Tape>, Arc<Tape>)> {
        let mut chb = vec![BOTH; sb.len()];
        if eval_interval(sb, ctx, bx, Some(&mut chb)).none_inside() {
            return None;
        }
        let mut cha = vec![BOTH; sa.len()];
        if eval_interval(sa, ctx, grow(bx, d), Some(&mut cha)).none_inside() {
            return None;
        }
        Some((sa.narrowed(&cha), sb.narrowed(&chb)))
    }

    /// Box bx's tests (`Clear`), of `asked`: the shapes simplified over it
    /// (`near`), or where it is flat its rests, for the clearance its gap
    /// leaves. `full` is the shapes as simplified, kept for seeking pairs.
    fn ask(&self, bx: [I; 3], ctx: &Ctx, full: (&Arc<Tape>, &Arc<Tape>), asked: (&Arc<Tape>, &Arc<Tape>), flat: Option<&Flat>) -> search::Test<Pair, Near> {
        use ink_core::interval::{add_up, lipschitz, mul_up, norm, sqrt_up};
        use search::Test::{Open, Settled};
        let (still, gap2) = flat.map_or(([false; 3], 0.0), |f| (f.still, f.gap2));
        let d2 = mul_up(self.d, self.d);
        if gap2 >= d2 {
            return Settled;
        }
        let d = if gap2 > 0.0 { sqrt_up(add_up(d2, -gap2)).min(self.d) } else { self.d };
        let (sa, sb) = match flat {
            Some(_) => match self.near(bx, ctx, d, asked.0, asked.1) {
                Some(x) => x,
                None => return Settled,
            },
            None => (asked.0.clone(), asked.1.clone()),
        };
        let grown = grow(bx, d);
        let (c, hw) = centre_reach(bx);
        let hw: V3 = std::array::from_fn(|k| if still[k] { 0.0 } else { hw[k] });
        let r = norm(hw);

        // The centre values as point intervals (outward-rounded, so sound),
        // first order in the ranged params: the mean-value forms stay tight
        // where a formula has dependencies that inflate plain bounds
        // (helices, rotations).
        let ga = ink_core::eval::eval_interval_grad(&sa, ctx, grown);
        let la = lipschitz(&ga);
        let ac = eval_interval_scope(&sa, ctx, c.map(I::point), None).lo;
        let alo = eval_interval(&sa, ctx, bx, None).lo;
        if la.is_finite() && ac > 0.0 && (d > 0.0 && alo >= mul_up(d, la) || ac >= mul_up(la, add_up(r, d))) {
            return Settled;
        }
        if la.is_finite() && crate::centred::centred_at_least(&[(&sb, 0.0), (&sa, mul_up(la, d))], ctx, bx, 0.0) {
            return Settled;
        }
        // Over the scope's axes (ranged params; time over a step of a
        // motion): bounds first order in them, and how much each blurs
        // them.
        let (mut params, mut t) = (Vec::new(), 0.0);
        let mut ps = scope::splittable(&sa, ctx, scope::PARAM_SPLITS);
        ps.extend(scope::splittable(&sb, ctx, scope::PARAM_SPLITS));
        let ts = scope::spans_t(&sa, ctx) || scope::spans_t(&sb, ctx);
        if !ps.is_empty() || ts {
            let (ra, blur_a) = ink_core::eval::eval_interval_scope_blur(&sa, ctx, grown, None);
            let (rb, blur_b) = ink_core::eval::eval_interval_scope_blur(&sb, ctx, bx, None);
            if ra.none_inside() || rb.none_inside() {
                return Settled;
            }
            // Each axis's blur of A's and of B's bounds as a length (over
            // each one's slope bound), the larger.
            let lb = lipschitz(&ink_core::eval::eval_interval_grad(&sb, ctx, bx));
            let (shares, tb) = scope_shares(ctx, &[(&blur_a, la), (&blur_b, lb)]);
            ps.sort_unstable();
            ps.dedup();
            params = ps.iter().map(|&k| (k, share(&shares, k))).collect();
            t = if ts { tb } else { 0.0 };
        }
        // Every point of B here is at least a / L from A, where a > 0 there:
        // a(p) >= alo, or a(centre) - L r, and A's field rises no faster
        // than L within d of the box.
        let bound = if la.is_finite() && la > 0.0 { (alo / la).next_down().max(((ac / la).next_down() - r.next_up()).next_down()) } else { 0.0 };
        let space: V3 = std::array::from_fn(|k| match slope_share(&ga, la, k) {
            0.0 => 0.0,
            s => s * hw[k],
        });
        // Flat, the rests' bound and the gap in quadrature.
        let bound = bound.clamp(0.0, d);
        let bound = if gap2 > 0.0 { ((bound * bound).next_down() + gap2).next_down().max(0.0).sqrt().next_down().clamp(0.0, self.d) } else { bound };
        let space = std::array::from_fn(|k: usize| if still[k] { 0.0 } else { space[k] });
        let state = match flat {
            Some(_) => (full.0.clone(), full.1.clone(), Some(Flat { a: sa, b: sb, still, gap2 })),
            None => (sa, sb, None),
        };
        Open { bound, smear: search::Smear { space, params, t, still }, state }
    }

    /// From q in the box, Newton steps toward A's surface along the box's
    /// tape `sa` (`t`: A over the box grown by d, B over the box),
    /// overshooting by a sliver of d so a step lands inside: a point of A
    /// within d of q (so where `sa` is A), confirmed with B's and A's
    /// whole fields by bounds rounded outward?
    fn walk(&self, sa: &Arc<Tape>, ctx: &Ctx, q: V3) -> Option<Near> {
        let mut p = q;
        for _ in 0..12 {
            let g = ink_core::eval::eval_grads(sa, ctx, &[p])[0];
            if inside(g[0]) {
                return self.confirmed(ctx, q, p);
            }
            let n2 = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
            if !(n2 > 0.0) {
                return None;
            }
            let k = (g[0] + 1e-3 * self.d) / n2;
            p = [p[0] - k * g[1], p[1] - k * g[2], p[2] - k * g[3]];
        }
        None
    }

    /// A pair of points closer than d near x, along the box's tapes: x
    /// onto B, that onto A, moved closer in turn (closer_pair); then each
    /// a quarter of the room left (d less their gap) deeper into its shape
    /// along its slope, so the pair stays within d and inside where the
    /// scope moves the shapes a little (a touches claim over ranged
    /// params). Confirmed with the whole fields (`confirmed`), the deeper
    /// pair first.
    fn pair(&self, (sa, sb): (&Arc<Tape>, &Arc<Tape>), ctx: &Ctx, x: V3) -> Option<Near> {
        let onto = |t: &Tape, p: V3| match Confirm::new(ctx).below(t, p, 0.0) {
            true => Some(p),
            false => to_surface(t, ctx, p, 1e-12 * self.scale).and_then(|y| confirm_inside(t, ctx, y, self.scale)),
        };
        let q = onto(sb, x)?;
        let p = onto(sa, q)?;
        let (gap, p, q) = closer_pair(sa, sb, ctx, p, q, self.scale);
        if !(gap < self.d) {
            return None;
        }
        let room = 0.25 * (self.d - gap);
        let deeper = |t: &Tape, x: V3| {
            let g = eval_grads(t, ctx, &[x])[0];
            let n = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
            [x[0] - room * g[1] / n, x[1] - room * g[2] / n, x[2] - room * g[3] / n]
        };
        self.confirmed(ctx, deeper(sb, q), deeper(sa, p)).or_else(|| self.confirmed(ctx, q, p))
    }

    /// q in the region and B, p in A and within d of q, by B's and A's
    /// whole fields (A's holds within d of the region) and bounds rounded
    /// outward: a witness pair.
    fn confirmed(&self, ctx: &Ctx, q: V3, p: V3) -> Option<Near> {
        let mut confirm = Confirm::new(ctx);
        let ok = q.iter().chain(&p).all(|v| v.is_finite()) && (0..3).all(|k| self.region[k].contains(q[k])) && distance_up(p, q) < self.d;
        (ok && confirm.below(self.b, q, 0.0) && confirm.below(self.a, p, 0.0)).then(|| (q, p, At::of(ctx)))
    }
}

/// One way round for `clearance`: boxes cover `covered` in lo..hi (all of
/// it within d of the other shape; infinite ends: all of space, tested
/// whole), and `grown` is grown by d around each. `swapped`: the covered
/// shape is the claim's first (a pair found is reported the claim's way
/// round).
pub struct Way {
    pub grown: Arc<Tape>,
    pub covered: Arc<Tape>,
    pub lo: V3,
    pub hi: V3,
    pub swapped: bool,
}

/// No point of the second shape within d of the first, exactly: the one
/// search (crate::search) of each way round, its boxes halved down to
/// `finest`, testing at most `budget` boxes between them.
///
/// Either way decides it, and their costs can differ by orders of
/// magnitude: where one shape's field is far from a distance within d of
/// the other (a part cut by a half-space reads 0.1 at points 0.3 away),
/// the boxes covering the other stay undecided at any size. So each round
/// goes to the way with the least volume undecided (the first on a tie),
/// and a way that ends undecided leaves the rest to the other. A way's box
/// stuck at `finest` (the gap there within about it of d, or a field that
/// reads well under the distance near the other shape: a box of slabs
/// beside an edge) ends the way's proof, as best first it is the box
/// nearest a violation. Undecided both ways, both hunt a pair closer than d
/// in the rest (search::Search::hunt), in turn as before, each to its cap
/// or the budget (`hunt`); none found, it gives the way that proved more:
/// the accuracy it holds to (`attained`: every point of the covered
/// shape at least d less that from the other), its boxes left and the box
/// that ended it. The answer is a function of the question; it is cached
/// by it. Returns the verdict, a pair closer than d (the claim's way
/// round), and how finely it looked (search::Search::res).
pub fn clearance(ways: &[Way], ctx: &Ctx, d: f64, finest: f64, budget: usize) -> (Verdict, Option<(V3, V3)>, f64) {
    let tapes: Vec<&Tape> = ways.iter().flat_map(|w| [&*w.grown, &*w.covered]).collect();
    let region = ways.first().map_or(([0.0; 3], [0.0; 3]), |w| (w.lo, w.hi));
    let inputs = |fp: &mut ink_core::graph::Fingerprint| {
        for w in ways {
            w.lo.iter().chain(&w.hi).for_each(|v| _ = fp.f64(*v));
            fp.word(w.swapped as u64);
        }
        fp.f64(d).f64(finest).word(budget as u64);
    };
    let (v, pair, h) = crate::cache::memo_region("clearance", &tapes, ctx, region, inputs, |c| {
        let (v, pair, h, _) = clearance_now(ways, c, d, finest, budget);
        (v, pair, h)
    });
    (at_moment(v, ctx, &tapes), pair, h)
}

/// `clearance`, uncached, and the boxes it tested.
pub(super) fn clearance_now(ways: &[Way], ctx: &Ctx, d: f64, finest: f64, budget: usize) -> (Verdict, Option<(V3, V3)>, f64, usize) {
    let ctx = Arc::new(in_span(ctx));
    // The grown shape within d of the region (its cells clipped to it,
    // each grown by d): all of it that matters.
    let near: Vec<Arc<Tape>> = ways.iter().map(|w| focus(&w.grown, &ctx, w.lo.map(|x| (x - d).next_down()), w.hi.map(|x| (x + d).next_up()))).collect();
    let qs: Vec<Clear> = ways
        .iter()
        .zip(&near)
        .map(|(w, a)| Clear { a, b: &w.covered, d, region: std::array::from_fn(|k| I::new(w.lo[k], w.hi[k])), scale: max_extent(w.lo, w.hi).min(f64::MAX) })
        .collect();
    let mut searches: Vec<search::Search<Pair, Near>> = ways
        .iter()
        .zip(&qs)
        .map(|(w, q)| {
            let b = std::array::from_fn(|k| I::new(w.lo[k], w.hi[k]));
            search::Search::new(q, vec![search::Node { b, ctx: ctx.clone(), state: (q.a.clone(), w.covered.clone(), None) }], search::Precision::holds(finest))
        })
        .collect();
    let mut ended: Vec<Option<search::Undecided>> = vec![None; ways.len()];
    let swap = |w: &Way, (q, p): (V3, V3)| if w.swapped { (p, q) } else { (q, p) };
    let tested = |s: &[search::Search<Pair, Near>]| s.iter().map(|s| s.tested()).sum::<usize>();
    for hunt in [false, true] {
        // Neither way proves it: unless stopped, both hunt a pair closer
        // than d (search::Search::hunt), in turn as before.
        if hunt {
            if ended.iter().flatten().any(|u| u.stopped()) {
                break;
            }
            searches.iter_mut().for_each(|s| s.hunt());
            ended.iter_mut().for_each(|u| *u = None);
        }
        while let Some(k) = (0..ways.len()).filter(|&k| ended[k].is_none()).min_by(|&i, &j| searches[i].volume().total_cmp(&searches[j].volume())) {
            let others: usize = searches.iter().enumerate().filter(|x| x.0 != k).map(|x| x.1.tested()).sum();
            match searches[k].round(&qs[k], budget.saturating_sub(others)) {
                None => {}
                Some(search::Outcome::Proven) => return (Verdict::Proven, None, searches[k].res(), tested(&searches)),
                Some(search::Outcome::Found((q, p, at))) => {
                    let v = Verdict::Disproven { witness: q, t: at.t, value: distance_up(p, q), params: at.params };
                    return (v, Some(swap(&ways[k], (q, p))), searches[k].res(), tested(&searches));
                }
                Some(search::Outcome::Unknown(u)) => ended[k] = Some(u),
            }
        }
    }
    // Both ways' hunts: their cells, stopped by the time limit if either
    // was, by the refinement limit if both were, else by a cap or the budget.
    let hunts: Vec<search::Hunt> = ended.iter().flatten().filter_map(|u| u.hunt).collect();
    let hunt = (!hunts.is_empty()).then(|| {
        use search::HuntEnd::{Budget, Cap, Stuck, Time};
        let by = |e| hunts.iter().any(|h| h.stopped == e);
        let stopped = if by(Time) {
            Time
        } else if hunts.iter().all(|h| h.stopped == Stuck) {
            Stuck
        } else if by(Cap) {
            Cap
        } else {
            Budget
        };
        search::Hunt { cells: hunts.iter().map(|h| h.cells).sum(), stopped }
    });
    // The way that proved the larger clearance (a stop ends both).
    let Some((k, u)) =
        ended.into_iter().enumerate().filter_map(|(k, u)| Some((k, u?))).max_by(|x, y| x.1.stopped().cmp(&y.1.stopped()).then(x.1.bound.total_cmp(&y.1.bound)).then(y.0.cmp(&x.0)))
    else {
        return (Verdict::Unknown { cells: 0, volume: 0.0, ended: Ended::Finest, cell: None, attained: None, hunt: None }, None, finest, tested(&searches));
    };
    let attained = (u.bound > 0.0).then(|| (d - u.bound).next_up());
    (Verdict::Unknown { cells: u.cells, volume: u.volume, ended: u.ended, cell: u.cell, attained, hunt }, None, searches[k].res(), tested(&searches))
}

/// The minimum distance between two shapes, bracketed (`min_distance`):
/// every point of each at least `lower` from the other, proven; a pair
/// confirmed inside them (on the first, on the second, the claim's way
/// round) and its distance rounded up, `upper`; why the search ended short
/// of the accuracy (None: it reached it); and how finely it looked.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct MinDistance {
    pub lower: f64,
    pub upper: Option<(f64, V3, V3)>,
    pub ended: Option<Ended>,
    pub res: f64,
    /// Boxes tested.
    pub tested: usize,
}

/// The minimum distance between the shapes of `ways` (`clearance`'s, each
/// way round), to `tol`: one search, clearance's question asked of a d
/// that shrinks. It starts at `far` (no pair is farther apart); each pair
/// a box's seek finds closer than d, moved closer by alternating
/// projections (closer_pair), is the upper bound, and d becomes it less
/// `tol`; it starts from `seed`'s where given (a pair confirmed inside
/// both, the claim's way round: (distance, on the first, on the second)). The boxes already tested keep what they proved (a box settled at
/// a d is settled at any smaller one, and one whose bound reaches the new d
/// is settled by it: Search::settle_from), so the search goes on where it
/// was. Proven at d, the distance is within [d, upper]; ended undecided,
/// `lower` is the bound the way that proved more attained. Boxes halved to
/// `finest`, at most `budget` tested. Cached by the question.
pub fn min_distance(ways: &[Way], ctx: &Ctx, far: f64, seed: Option<(f64, V3, V3)>, tol: f64, finest: f64, budget: usize) -> MinDistance {
    let tapes: Vec<&Tape> = ways.iter().flat_map(|w| [&*w.grown, &*w.covered]).collect();
    let region = ways.first().map_or(([0.0; 3], [0.0; 3]), |w| (w.lo, w.hi));
    let inputs = |fp: &mut ink_core::graph::Fingerprint| {
        for w in ways {
            w.lo.iter().chain(&w.hi).for_each(|v| _ = fp.f64(*v));
            fp.word(w.swapped as u64);
        }
        fp.f64(far).f64(tol).f64(finest).word(budget as u64);
        if let Some((g, p, q)) = seed {
            [g].iter().chain(&p).chain(&q).for_each(|v| _ = fp.f64(*v));
        }
    };
    crate::cache::memo_region("min distance", &tapes, ctx, region, inputs, |c| min_distance_now(ways, c, far, seed, tol, finest, budget))
}

fn min_distance_now(ways: &[Way], ctx: &Ctx, far: f64, seed: Option<(f64, V3, V3)>, tol: f64, finest: f64, budget: usize) -> MinDistance {
    let ctx = Arc::new(in_span(ctx));
    let mut upper = seed;
    let mut d = seed.map_or(far, |u| u.0 - tol);
    if !(d > 0.0) {
        return MinDistance { lower: 0.0, upper, ended: None, res: f64::INFINITY, tested: 0 };
    }
    let near: Vec<Arc<Tape>> = ways.iter().map(|w| focus(&w.grown, &ctx, w.lo.map(|x| (x - d).next_down()), w.hi.map(|x| (x + d).next_up()))).collect();
    let region = |w: &Way| std::array::from_fn(|k| I::new(w.lo[k], w.hi[k]));
    let scale = |w: &Way| max_extent(w.lo, w.hi).min(f64::MAX);
    let ask = |k: usize, d: f64| Clear { a: &near[k], b: &ways[k].covered, d, region: region(&ways[k]), scale: scale(&ways[k]) };
    let mut searches: Vec<search::Search<Pair, Near>> = ways
        .iter()
        .enumerate()
        .map(|(k, w)| {
            search::Search::new(
                &ask(k, d),
                vec![search::Node { b: region(w), ctx: ctx.clone(), state: (near[k].clone(), w.covered.clone(), None) }],
                search::Precision::holds(finest),
            )
        })
        .collect();
    let mut ended: Vec<Option<search::Undecided>> = vec![None; ways.len()];
    let tested = |s: &[search::Search<Pair, Near>]| s.iter().map(|s| s.tested()).sum::<usize>();
    let res = |s: &[search::Search<Pair, Near>]| s.iter().map(|s| s.res()).fold(f64::INFINITY, f64::min);
    while let Some(k) = (0..ways.len()).filter(|&k| ended[k].is_none()).min_by(|&i, &j| searches[i].volume().total_cmp(&searches[j].volume())) {
        let others: usize = searches.iter().enumerate().filter(|x| x.0 != k).map(|x| x.1.tested()).sum();
        match searches[k].round(&ask(k, d), budget.saturating_sub(others)) {
            None => {}
            Some(search::Outcome::Proven) => return MinDistance { lower: d.max(0.0), upper, ended: None, res: res(&searches), tested: tested(&searches) },
            Some(search::Outcome::Found((q, p, _))) => {
                let w = &ways[k];
                let (gap, p, q) = closer_pair(&w.grown, &w.covered, &ctx, p, q, scale(w));
                // The claim's way round: on the first, on the second.
                let pair = if w.swapped { (gap, q, p) } else { (gap, p, q) };
                if upper.is_none_or(|u| gap < u.0) {
                    upper = Some(pair);
                }
                d = gap - tol;
                if !(d > 0.0) {
                    return MinDistance { lower: 0.0, upper, ended: None, res: res(&searches), tested: tested(&searches) };
                }
                searches.iter_mut().for_each(|s| s.settle_from(d));
            }
            Some(search::Outcome::Unknown(u)) => ended[k] = Some(u),
        }
    }
    // The way that proved more (a stop ends both).
    let u = ended.into_iter().flatten().max_by(|x, y| x.stopped().cmp(&y.stopped()).then(x.bound.total_cmp(&y.bound)));
    let lower = u.as_ref().map_or(0.0, |u| u.bound.clamp(0.0, d.max(0.0)));
    MinDistance { lower, upper, ended: Some(u.map_or(Ended::Finest, |u| u.ended)), res: res(&searches), tested: tested(&searches) }
}
