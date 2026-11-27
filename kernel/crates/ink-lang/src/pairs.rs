//! Pairs of a target's items: those whose boxes can meet (the broad phase
//! of `ink interfere`'s proofs), and a pair asked in its own frame
//! (`framed`), so every pair placed alike asks the same question.

use crate::view::ModelView;
use crate::{Doc, Item, Link, Xf};
use ink_core::V3;
use ink_core::eval::Ctx;
use ink_core::graph::{BinOp, UnOp};
use ink_core::interval::{I, add_o, add_up, bin, mag, mul_o, mul_up, sqrt_up, sub_o};
use ink_core::tape::Tape;
use ink_kernel::query;
use ink_kernel::rigid::Rigid;

/// The pairs of items whose boxes overlap (each item's box and whether it
/// is proven bounded given), and how many pairs their boxes set apart:
/// all, and for each of `kinds`, those both of whose items are of it.
/// Over a motion, `over` gives an item's boxes over pieces of it where it
/// moves (crate::motion::item_boxes; empty: its box holds it throughout):
/// a pair whose boxes meet at no moment is apart too (motion::apart_over).
pub fn candidate_pairs(bounds: &[Option<(V3, V3, f64)>], whole: &[bool], over: &[Vec<((f64, f64), (V3, V3))>], kinds: &[&[bool]]) -> (Vec<(usize, usize)>, usize, Vec<usize>) {
    // Broad phase: sort by min x and sweep; only boxes that overlap on every
    // axis become candidates. Everything else is proven apart by its box.
    // Items not proven bounded have no box that holds them: they are
    // candidates against every other item.
    let throughout: Vec<[((f64, f64), (V3, V3)); 1]> = bounds.iter().map(|b| [((f64::NEG_INFINITY, f64::INFINITY), b.map_or(([0.0; 3], [0.0; 3]), |b| (b.0, b.1)))]).collect();
    let pieces = |i: usize| -> &[((f64, f64), (V3, V3))] { over.get(i).filter(|p| !p.is_empty()).map_or(&throughout[i], |p| p) };
    let moves = |i: usize| over.get(i).is_some_and(|p| !p.is_empty());
    let mut order: Vec<usize> = (0..bounds.len()).filter(|&i| whole[i] && bounds[i].is_some()).collect();
    order.sort_by(|&i, &j| bounds[i].unwrap().0[0].total_cmp(&bounds[j].unwrap().0[0]));
    let mut candidates = Vec::new();
    for (k, &i) in order.iter().enumerate() {
        let bi = bounds[i].unwrap();
        for &j in &order[k + 1..] {
            let bj = bounds[j].unwrap();
            if bj.0[0] > bi.1[0] {
                break;
            }
            if (1..3).all(|d| bi.0[d] <= bj.1[d] && bj.0[d] <= bi.1[d]) && !((moves(i) || moves(j)) && crate::motion::apart_over(pieces(i), pieces(j))) {
                candidates.push((i.min(j), i.max(j)));
            }
        }
    }
    let pairs = |n: usize| n * (n.max(1) - 1) / 2;
    let boxed_apart = pairs(order.len()) - candidates.len();
    // Of those, the pairs apart whose items are both of each kind.
    let apart_among = |kind: &[bool]| pairs(order.iter().filter(|&&i| kind[i]).count()) - candidates.iter().filter(|&&(i, j)| kind[i] && kind[j]).count();
    let among = kinds.iter().map(|k| apart_among(k)).collect();
    for i in (0..bounds.len()).filter(|&i| !whole[i]) {
        for j in (0..bounds.len()).filter(|&j| j != i && !(whole[j] && bounds[j].is_none())) {
            candidates.push((i.min(j), i.max(j)));
        }
    }
    candidates.sort();
    candidates.dedup();
    (candidates, boxed_apart, among)
}

/// The share of the accuracy a pair asked in a shared frame (`framed`)
/// leaves the frame's rounding: a fixed share, so the question is the same
/// wherever the frame lies.
pub const FRAME_SHARE: f64 = 1.0 / (1u64 << 20) as f64;

/// A pair of items asked in its own frame (`framed`): that of one of its
/// parts, the other and its zones placed relative to it, so every pair
/// placed alike asks the same question, wherever it lies.
pub struct Framed {
    /// Each item's field in the frame, and its box there (its body's,
    /// moved: ModelView::local_at).
    pub items: [(u32, (V3, V3, f64)); 2],
    /// Its zones in the frame, their union.
    pub zone: Option<u32>,
    /// Where they can meet: the intersection of their boxes, padded (as
    /// interfere::pair_region); apart, the box between them.
    pub region: (V3, V3),
    /// The frame's placement in the world, its numbers, and its entries'
    /// bounds over the view's scope (R row by row, then t).
    pub at: Rigid,
    entries: [I; 12],
    /// How far each field's solid read back in the frame from the world
    /// is from the question's: a point inside one has a point inside the
    /// other within this.
    pub error: f64,
    /// Least and greatest stretch of lengths by the frame's placement, and
    /// |det R|.
    pub stretch: (f64, f64),
    pub det: I,
    /// R^T R = I over the scope: R's inverse is its transpose.
    exact_turn: bool,
}

impl Framed {
    /// The accuracy a pair asked to `tol` is asked to in its frame.
    pub fn asked(tol: f64) -> f64 {
        tol * (1.0 - FRAME_SHARE)
    }

    /// A point of the frame in the world (by the numbers), and how far it
    /// may be from where the frame's placement puts it.
    pub fn world(&self, c: V3) -> (V3, f64) {
        let p = self.at.apply(c);
        let mut e = 0.0;
        for (i, &pi) in p.iter().enumerate() {
            let mut x = self.entries[9 + i];
            for (j, &cj) in c.iter().enumerate() {
                x = add_o(x, mul_o(self.entries[3 * i + j], I::point(cj)));
            }
            let d = mag(sub_o(x, I::point(pi)));
            e = add_up(e, mul_up(d, d));
        }
        (p, sqrt_up(e))
    }

    /// The half-side of a cube about a witness in the frame that, inside
    /// the pair there, puts the cube of half-side `tol` about its world
    /// point (`off` from where the frame puts it: `world`) inside the
    /// world's pair: a ball of radius h in the question holds one of h less
    /// the error in the world's read back. That cube lies in the ball of
    /// radius √3 tol, which the frame stretches by at least its least; an
    /// exact turn (R^T R = I) carries it back to a cube whose box has
    /// half-side tol times R's greatest column sum of |entries| (1 for a
    /// quarter turn), and the offset to a ball of radius `off`.
    pub fn cube(&self, tol: f64, off: f64) -> f64 {
        let ball = mul_up(tol, sqrt_up(3.0));
        let general = add_up(bin(BinOp::Div, I::point(add_up(ball, off)), I::point(self.stretch.0)).0.hi, self.error);
        if !self.exact_turn {
            return general;
        }
        let column = (0..3).map(|j| (0..3).fold(0.0, |s, i| add_up(s, mag(self.entries[3 * i + j])))).fold(0.0, f64::max);
        general.min(add_up(add_up(mul_up(tol, column), off), self.error))
    }
}

/// The pair `a`, `b` of a target's items (with `zones`, each a field in
/// the coordinates of the instance path it names) asked in a frame of its
/// own. That of the deepest instance holding both and each zone's path:
/// each field placed there by the links below it (Q), the same nodes in
/// every copy of that instance; and where those are numbers, the frame of
/// one of the parts, the other's and the zones' Q relative to it, rounded
/// (`anchored`), so pairs placed alike anywhere are one question. The
/// world places each by all its links (W), and the frame lies at S, the
/// links above it (and the anchor's below). Read back in the frame, the world
/// reads a field at ρ_W(S c) and the question at ρ_Q(c) (ρ_X(p) = R^T (p -
/// t), as Doc::warped reads a placement): they differ by Δ(c) =
/// (R_W^T R_S - R_Q^T) c + R_W^T (t_S - t_W) + R_Q^T t_Q, bounded by interval
/// arithmetic on every entry's bounds over the region (grown by twice the
/// accuracy, where either's material near it lies). Reading c + (R_Q^T +
/// M)^-1 Δ(c) instead of c turns one's value into the other's, a move of
/// at most |Δ| / (σ_min(R_Q) - |M|): that is each field's `error`, so an
/// intersection of them that holds a ball of radius r in one holds the
/// ball of radius r less it about the same centre in the other. A union of
/// fields moved apart does not (its seams open): a zone with a part in each
/// end's coordinates is asked in the world. S stretches lengths by its
/// singular values, within its skew of 1 squared (rigid::turn_bounds). So
/// the pair holds in the world where the question holds to `asked` and
/// σ_max (asked + error) <= `tol`; a witness's cube is moved by `cube`.
/// Over a motion that moves the frame (the pair rides one rigid group:
/// crate::motion::rigid_group), nothing below it moves, so the question
/// holds at every moment; then M and Δ's other term, and S's stretch, are
/// bounded over the whole span, each entry read as a polynomial
/// (rigid::frame_delta, turn_bounds): S and W read the motion's terms alike,
/// which cancel, where each entry's bounds over the span would differ by
/// the whole travel. None: asked in the world, where it is placed exactly
/// as the question places it (S the identity), an item is no rigid body
/// bounded in its frame, an entry is not finite or not read, or the error
/// does not fit the share.
pub fn framed(v: &mut ModelView, a: &Item, b: &Item, zones: &[(u32, &str)], tol: f64) -> Option<Framed> {
    if !(tol > 0.0) {
        return None;
    }
    let mut k = a.links.iter().zip(&b.links).take_while(|(x, y)| x.path == y.path).count();
    let mut ends = Vec::with_capacity(zones.len());
    for &(f, end) in zones {
        let (it, m) = [a, b].into_iter().find_map(|it| it.links.iter().position(|l| l.path == end).map(|m| (it, m + 1)))?;
        k = k.min(m);
        ends.push((f, it, m));
    }
    let fold = |doc: &mut Doc, from: Xf, links: &[Link]| links.iter().fold(from, |x, l| doc.compose(&l.at, &x));
    // An item's links run from the world down (the target's own prefix
    // among them: Doc::items).
    let id = v.doc.identity();
    let mut s = fold(&mut v.doc, id, &a.links[..k]);
    // Each field read: its own, where the question places it, where the
    // world does.
    let mut fields: Vec<(u32, Xf, Xf)> = Vec::with_capacity(2 + ends.len());
    for it in [a, b] {
        let (_, w) = it.instance.as_ref()?;
        let q = fold(&mut v.doc, id, &it.links[k..]);
        fields.push((it.shape, q, *w));
    }
    for &(f, it, m) in &ends {
        let q = fold(&mut v.doc, id, &it.links[k..m]);
        let w = fold(&mut v.doc, id, &it.links[..m]);
        fields.push((f, q, w));
    }
    if let Some((n, q)) = anchored(v, [a, b], &fields, tol) {
        s = v.doc.compose(&fields[n].1, &s);
        for (f, q) in fields.iter_mut().zip(q) {
            f.1 = q;
        }
    }
    if s == id {
        return None;
    }
    // The span of the motion moving the frame, if one does.
    let entries_of = |x: &Xf| -> [u32; 12] { std::array::from_fn(|j| if j < 9 { x.r[j / 3][j % 3] } else { x.t[j - 9] }) };
    let moving = v.motion.filter(|&span| crate::claims::moves(v, &entries_of(&s), span));
    let mut items = [(0, ([0.0; 3], [0.0; 3], 0.0)); 2];
    for (n, it) in [a, b].into_iter().enumerate() {
        let (body, _) = it.instance.as_ref()?;
        let l = v.local_at(body, &fields[n].1)?;
        let (lo, hi, h) = v.bounds(l.shape)?;
        let (lo, hi) = l.world_box(lo, hi);
        items[n] = (v.doc.warped(it.shape, &fields[n].1), (lo, hi, h));
    }
    // Boxes apart in the frame (their world boxes, turned, met) leave the
    // pair nothing in common: the box between them is its region.
    let (ba, bb) = (items[0].1, items[1].1);
    let (lo, hi): (V3, V3) = (std::array::from_fn(|d| ba.0[d].max(bb.0[d])), std::array::from_fn(|d| ba.1[d].min(bb.1[d])));
    let (lo, hi): (V3, V3) = (std::array::from_fn(|d| lo[d].min(hi[d])), std::array::from_fn(|d| lo[d].max(hi[d])));
    let pad = 1e-6 * query::max_extent(lo, hi).max(1e-9);
    let region = (lo.map(|x| x - pad), hi.map(|x| x + pad));
    let mut zone = None;
    for (f, q, _) in &fields[2..] {
        let z = v.doc.warped(*f, q);
        zone = Some(zone.map_or(z, |u| v.doc.g.min(u, z)));
    }
    // Each entry's bounds over the view's scope.
    let bound = |v: &ModelView, n: u32| -> Option<I> {
        let e = match v.doc.g.as_const(n) {
            Some(c) => I::point(c),
            None => ink_core::eval::eval_interval(&Tape::new(&v.doc.g, n), &v.ctx, [I::ENTIRE; 3], None),
        };
        (e.lo.is_finite() && e.hi.is_finite()).then_some(e)
    };
    let read = |v: &ModelView, x: &Xf| -> Option<([[I; 3]; 3], [I; 3])> {
        let mut r = [[I::point(0.0); 3]; 3];
        let mut t = [I::point(0.0); 3];
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = bound(v, x.r[i][j])?;
            }
            t[i] = bound(v, x.t[i])?;
        }
        Some((r, t))
    };
    let (rs, ts) = read(v, &s)?;
    // S's turn: at the view's moment, or over the motion moving it.
    let scope = match moving {
        Some(span) => Ctx { t_span: Some(span), ..v.ctx.clone() },
        None => v.ctx.clone(),
    };
    let turn = ink_kernel::rigid::turn_bounds(&v.doc.g, &s.r, &scope)?;
    let grow = 2.0 * tol;
    let far: V3 = std::array::from_fn(|d| sub_o(I::point(region.0[d]), I::point(grow)).lo.abs().max(add_o(I::point(region.1[d]), I::point(grow)).hi.abs()));
    let least = |skew: f64| sub_o(I::point(1.0), I::point(skew)).lo.max(0.0).sqrt().next_down();
    let mut error = 0.0f64;
    for (_, q, w) in &fields {
        let skew = ink_kernel::rigid::turn_bounds(&v.doc.g, &q.r, &v.ctx)?.skew;
        // Δ = M c + d's terms: at the view's moment, from each entry's
        // bounds; over a motion moving the frame, read as polynomials.
        let (m, d): ([I; 9], [I; 3]) = match moving {
            Some(_) => ink_kernel::rigid::frame_delta(&v.doc.g, &scope, &entries_of(&s), &entries_of(q), &entries_of(w))?,
            None => {
                let ((rq, tq), (rw, tw)) = (read(v, q)?, read(v, w)?);
                let d = std::array::from_fn(|i| (0..3).fold(I::point(0.0), |d, k| add_o(add_o(d, mul_o(rw[k][i], sub_o(ts[k], tw[k]))), mul_o(rq[k][i], tq[k]))));
                let m = std::array::from_fn(|ij| {
                    let (i, j) = (ij / 3, ij % 3);
                    sub_o((0..3).fold(I::point(0.0), |m, k| add_o(m, mul_o(rw[k][i], rs[k][j]))), rq[j][i])
                });
                (m, d)
            }
        };
        let (mut delta, mut norm) = (0.0, 0.0);
        for i in 0..3 {
            let mut e = mag(d[i]);
            for j in 0..3 {
                let m = mag(m[3 * i + j]);
                e = add_up(e, mul_up(m, far[j]));
                norm = add_up(norm, mul_up(m, m));
            }
            delta = add_up(delta, mul_up(e, e));
        }
        let room = sub_o(I::point(least(skew)), I::point(sqrt_up(norm))).lo;
        if !(room > 0.0) {
            return None;
        }
        error = error.max(bin(BinOp::Div, I::point(sqrt_up(delta)), I::point(room)).0.hi);
    }
    let det = ink_core::interval::un(UnOp::Abs, turn.det);
    let stretch = (least(turn.skew), sqrt_up(add_up(1.0, turn.skew)));
    if !(det.lo > 0.0 && stretch.0 > 0.0 && mul_up(stretch.1, add_up(Framed::asked(tol), error)) <= tol) {
        return None;
    }
    let mut entries = [I::point(0.0); 12];
    for i in 0..3 {
        entries[3 * i..3 * i + 3].copy_from_slice(&rs[i]);
        entries[9 + i] = ts[i];
    }
    let num = |n: u32| v.doc.scalar(n, &v.ctx);
    let at = Rigid::from_entries(std::array::from_fn(|j| if j < 9 { num(s.r[j / 3][j % 3]) } else { num(s.t[j - 9]) }));
    Some(Framed { items, zone, region, at, entries, error, stretch, det, exact_turn: turn.skew == 0.0 })
}

/// The pair's own frame (`framed`): that of its anchor, the smaller of its
/// parts (by its body's box; the first of a tie), where every field's
/// placement below the shared instance is numbers (`ModelView::numbers`)
/// and the anchor's a rigid motion. Each field is then placed relative to
/// the anchor, Q_n^-1 Q = (R_n^T R, R_n^T (t - t_n)), every entry rounded
/// to a grid of powers of two: t's to the largest at most 2^-24 of the
/// accuracy, R's to the largest at most that over the anchor's reach (the
/// largest |coordinate| of its box, which holds the pair's region). So pairs
/// placed alike, whose floats differ by the rounding of where each lies (a
/// screw seated on each of 19 nuts), ask one question; the grid's rounding
/// and the floats' are `framed`'s error, charged to the accuracy. The
/// anchor's own placement is the identity. Returns the anchor and each
/// field's placement; None: the frame stays the shared instance's.
fn anchored(v: &mut ModelView, ends: [&Item; 2], fields: &[(u32, Xf, Xf)], tol: f64) -> Option<(usize, Vec<Xf>)> {
    let mut boxes = [([0.0; 3], [0.0; 3]); 2];
    let mut bodies = [""; 2];
    for (n, it) in ends.into_iter().enumerate() {
        bodies[n] = &it.instance.as_ref()?.0;
        let (lo, hi, _) = v.bounds(v.doc.body(bodies[n])?.shape)?;
        boxes[n] = (lo, hi);
    }
    let extent = |(lo, hi): (V3, V3)| query::max_extent(lo, hi);
    let n = usize::from(extent(boxes[1]) < extent(boxes[0]));
    v.local_at(bodies[n], &fields[n].1)?;
    let nums: Vec<Rigid> = fields.iter().map(|(_, q, _)| v.numbers(q)).collect::<Option<_>>()?;
    let anchor = nums[n];
    let (lo, hi) = boxes[n];
    let reach = lo.iter().chain(&hi).fold(tol, |m, x| m.max(x.abs()));
    let fine = tol / (1u64 << 24) as f64;
    // The largest power of two at most x (positive, normal).
    let below = |x: f64| (x >= f64::MIN_POSITIVE).then(|| f64::from_bits(x.to_bits() & 0xfff0_0000_0000_0000));
    let (gt, gr) = (below(fine)?, below(fine / reach)?);
    let id = v.doc.identity();
    let g = &mut v.doc.g;
    let mut on = |x: f64, step: f64| g.constant((x / step).round() * step + 0.0);
    let q = (0..nums.len())
        .map(|f| match f == n {
            true => id,
            false => {
                let x = nums[f];
                let r = std::array::from_fn(|i| std::array::from_fn(|j| on((0..3).map(|k| anchor.r[k][i] * x.r[k][j]).sum(), gr)));
                let t = std::array::from_fn(|i| on((0..3).map(|k| anchor.r[k][i] * (x.t[k] - anchor.t[k])).sum(), gt));
                Xf { r, t }
            }
        })
        .collect();
    Some((n, q))
}
