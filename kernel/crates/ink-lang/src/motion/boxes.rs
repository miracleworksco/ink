use super::*;
use std::collections::HashMap;

/// An item's box over the view's motion and whether it is proven bounded
/// (as ModelView::item_bounds and item_whole at a moment): a box holding
/// the item at every moment, for an interfere's broad phase, and where it
/// moves rigidly, its boxes over pieces of the span (`swept`: in order of
/// t, covering it), so that pairs whose boxes meet at no moment are apart
/// (pairs::candidate_pairs).
pub type ItemBox = (Option<(V3, V3, f64)>, bool, Vec<((f64, f64), (V3, V3))>);

/// Pieces of the span an item that moves rigidly is boxed over (`swept`).
const SWEPT_PIECES: usize = 64;

/// Classify a placed field or contact zone before a box or a single pose
/// can settle its pairs. A derivative can look still on both sides of a
/// jump and give no placement at the jump itself.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum FieldMotion {
    Still,
    Moving,
    UnboundedSpeed,
}

impl FieldMotion {
    pub(crate) fn of(v: &mut ModelView, f: u32) -> Self {
        if jumps(v, &Claim::Empty(f), &[f]).is_some() {
            Self::UnboundedSpeed
        } else if v.motion.is_some_and(|span| claims::moves(v, &[f], span)) {
            Self::Moving
        } else {
            Self::Still
        }
    }

    pub(crate) fn of_item(v: &mut ModelView, it: &crate::Item) -> Self {
        let f = v.doc.placed_field(it);
        Self::of(v, f)
    }
}

/// Each of `items`' ItemBox. A placed body proven bounded in its frame
/// (over the span, where its own field reads t) has its boxes over
/// SWEPT_PIECES uniform pieces of the span (`swept`), and their hull; any
/// other that moves has its proven bounds with t ranged over the whole span
/// (the one search halves t where it blurs them). One that reads a sampled
/// or solved value has no box: a candidate against every other item.
pub fn item_boxes(v: &mut ModelView, items: &[crate::Item]) -> Vec<ItemBox> {
    let motions: Vec<FieldMotion> = items.iter().map(|it| FieldMotion::of_item(v, it)).collect();
    boxes_of(v, items, &motions)
}

/// `item_boxes`, given how each item moves (`motions`). The moving items'
/// searches are one scope's (t over the whole span), so they are found
/// side by side, each once: every body's own box (copies of one body share
/// it), every isometry's sweep, then every other's bounds. Each is what it
/// is found alone, whatever the threads.
pub(crate) fn boxes_of(v: &mut ModelView, items: &[crate::Item], motions: &[FieldMotion]) -> Vec<ItemBox> {
    use rayon::prelude::*;
    let mut out: Vec<Option<ItemBox>> = (0..items.len()).map(|_| None).collect();
    let over = v.motion.map(|span| Ctx { t: ink_core::eval::span_mid(span), t_span: Some(span), ..v.ctx.clone() });
    // Each moving item's placed field, and its body and the isometry
    // placing it where it is one over the span.
    let mut moving: Vec<(usize, u32, Option<(u32, Isometry)>)> = Vec::new();
    for (k, (it, &motion)) in items.iter().zip(motions).enumerate() {
        out[k] = match motion {
            FieldMotion::Still => Some((v.item_bounds(it), v.item_whole(it), Vec::new())),
            FieldMotion::UnboundedSpeed => Some((None, false, Vec::new())),
            FieldMotion::Moving => {
                let over = over.as_ref().expect("a moving item spans time");
                let f = v.doc.placed_field(it);
                if v.doc.inexact(&[f], &[]).any() {
                    Some((None, false, Vec::new()))
                } else {
                    // Moved forward, R q + T, only by an isometry: the field
                    // reads R^T (p - T), whose solid R^-T moves, not R.
                    let iso = v.doc.placed.get(&f).map(|p| (p.shape, p.at)).and_then(|(s, at)| Isometry::over(&v.doc.g, &at.r, &at.t, over).map(|iso| (s, iso)));
                    moving.push((k, f, iso));
                    None
                }
            }
        };
    }
    if let Some(over) = over.filter(|_| !moving.is_empty()) {
        let span = over.t_span.expect("over a span");
        // Each body's own box: at rest, or over the span where its own
        // field reads t (a belt's teeth running round its loop).
        let mut bodies: Vec<u32> = moving.iter().filter_map(|m| m.2.map(|(s, _)| s)).collect();
        bodies.sort_unstable();
        bodies.dedup();
        let own: HashMap<u32, Option<(V3, V3)>> =
            v.at_point(over.clone(), |v| bodies.iter().copied().zip(found(v, &bodies)).map(|(s, (w, b))| (s, b.filter(|_| w).map(|b| (b.0, b.1)))).collect());
        let (g, ctx) = (&v.doc.g, &v.ctx);
        let pieces: Vec<Option<Vec<((f64, f64), (V3, V3))>>> =
            moving.par_iter().map(|(_, _, iso)| iso.as_ref().and_then(|(s, iso)| swept(g, ctx, iso, own[s]?, span, true))).collect();
        let mut world = Vec::new();
        for ((k, f, _), pieces) in moving.iter().zip(pieces) {
            match pieces {
                Some(pieces) => {
                    let (lo, hi) = pieces.iter().map(|p| p.1).reduce(ink_core::v3::box_union).expect("pieces");
                    out[*k] = Some((Some((lo, hi, 0.0)), true, pieces));
                }
                None => world.push((*k, *f)),
            }
        }
        // The rest: their placed fields' bounds over the span.
        let fields: Vec<u32> = world.iter().map(|w| w.1).collect();
        let boxes = v.at_point(over, |v| found(v, &fields));
        for ((k, _), (w, b)) in world.into_iter().zip(boxes) {
            out[k] = Some((b, w, Vec::new()));
        }
    }
    out.into_iter().map(|b| b.expect("every item boxed")).collect()
}

/// Whether each of `fs` is proven bounded, and its bounds where it is
/// (ModelView::whole, bounds), found at once, in parallel.
fn found(v: &mut ModelView, fs: &[u32]) -> Vec<(bool, Option<(V3, V3, f64)>)> {
    v.reach_all(fs);
    let whole: Vec<u32> = fs.iter().copied().filter(|&f| v.whole(f)).collect();
    v.bounds_all(&whole);
    fs.iter().map(|&f| if v.whole(f) { (true, v.bounds(f)) } else { (false, None) }).collect()
}

/// Boxes holding a body (its box (lo, hi) in its own frame) placed by
/// `at` (an isometry over the span: so its solid is where R p + T puts it)
/// over SWEPT_PIECES uniform pieces of `span`, each holding it at every
/// moment of its piece: R p + T with R's and T's entries bounded over the
/// piece. None: an entry not bounded.
pub(super) fn swept(g: &ink_core::graph::Graph, ctx: &Ctx, at: &Isometry, (lo, hi): (V3, V3), span: (f64, f64), refine: bool) -> Option<Vec<((f64, f64), (V3, V3))>> {
    use interval::I;
    let entry = |n: u32| Tape::new(g, n);
    let r: Vec<Tape> = at.r().iter().flatten().map(|&n| entry(n)).collect();
    let t: Vec<Tape> = at.t().iter().map(|&n| entry(n)).collect();
    let p: [I; 3] = std::array::from_fn(|k| I::new(lo[k], hi[k]));
    let o = [I::point(0.0); 3];
    let size = |b: &(V3, V3)| (0..3).map(|k| b.1[k] - b.0[k]).fold(0.0, f64::max);
    let own = size(&(lo, hi));
    let boxed = |s: (f64, f64)| -> Option<(V3, V3)> {
        let c = Ctx { t: ink_core::eval::span_mid(s), t_span: Some(s), ..ctx.clone() };
        let mut b = ([0.0; 3], [0.0; 3]);
        for j in 0..3 {
            let mut x = ink_core::eval::eval_interval(&t[j], &c, o, None);
            for k in 0..3 {
                x = interval::bin(BinOp::Add, x, interval::bin(BinOp::Mul, ink_core::eval::eval_interval(&r[3 * j + k], &c, o, None), p[k]).0).0;
            }
            if !(x.lo.is_finite() && x.hi.is_finite()) {
                return None;
            }
            (b.0[j], b.1[j]) = (x.lo, x.hi);
        }
        Some(b)
    };
    // SWEPT_PIECES uniform pieces, each (`refine`) halved while its box
    // reaches past the body's own size by more than BLUR_SHARE of it, down
    // to a SWEPT_FINEST-th of the span.
    let n = SWEPT_PIECES;
    let finest = (span.1 - span.0) / SWEPT_FINEST as f64;
    let mut todo: Vec<(f64, f64)> = (0..n)
        .rev()
        .map(|i| (span.0 + (span.1 - span.0) * i as f64 / n as f64, if i + 1 == n { span.1 } else { span.0 + (span.1 - span.0) * (i + 1) as f64 / n as f64 }))
        .collect();
    let mut out = Vec::new();
    while let Some(s) = todo.pop() {
        let b = boxed(s)?;
        let m = ink_core::eval::span_mid(s);
        if refine && size(&b) > own * (1.0 + BLUR_SHARE) && s.1 - s.0 > finest && s.0 < m && m < s.1 {
            todo.extend([(m, s.1), (s.0, m)]);
            continue;
        }
        out.push((s, b));
    }
    Some(out)
}

/// The share of a body's size its box over a piece of t may reach past
/// it before the piece is halved (`swept`).
const BLUR_SHARE: f64 = 1.0 / 8.0;

/// The finest piece of the span `swept` halves to: this share of it.
const SWEPT_FINEST: usize = 1 << 10;

/// Whether two items' boxes over pieces of a motion (`swept`, in order of
/// t, covering the span; empty: its box holds it throughout) meet at no
/// moment: each pair of pieces that share a moment, in one pass.
pub fn apart_over(a: &[((f64, f64), (V3, V3))], b: &[((f64, f64), (V3, V3))]) -> bool {
    let meet = |x: &(V3, V3), y: &(V3, V3)| (0..3).all(|k| x.0[k] <= y.1[k] && y.0[k] <= x.1[k]);
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        let ((sa, ba), (sb, bb)) = (&a[i], &b[j]);
        if sb.0 <= sa.1 && sa.0 <= sb.1 && meet(ba, bb) {
            return false;
        }
        // The piece ending first is done with (closed: at a tie, the other
        // is still met against the next, which shares its end).
        if sa.1 <= sb.1 { i += 1 } else { j += 1 }
    }
    true
}
