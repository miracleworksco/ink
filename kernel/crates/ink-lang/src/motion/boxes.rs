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
/// (over the span, where its own field reads t), placed rigidly (`Sweep`),
/// has its boxes over SWEPT_PIECES uniform pieces of the span (`swept`),
/// and their hull; any other that moves has its proven bounds with t
/// ranged over the whole span (the one search halves t where it blurs
/// them). One that reads a sampled or solved value has no box: a candidate
/// against every other item.
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
    // Each moving item's placed field, and its body and the placement
    // swept where it is rigid over the span; each link's turn read once.
    let mut moving: Vec<(usize, u32, Option<(u32, Sweep)>)> = Vec::new();
    let mut turns: HashMap<Xf, Option<TurnBounds>> = HashMap::new();
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
                    let sweep = v.doc.placed.get(&f).map(|p| (p.shape, p.at)).and_then(|(s, at)| Some((s, sweep(v, it, &at, over, &mut turns)?)));
                    moving.push((k, f, sweep));
                    None
                }
            }
        };
    }
    if let Some(over) = over.filter(|_| !moving.is_empty()) {
        let span = over.t_span.expect("over a span");
        // Each body's own box: at rest, or over the span where its own
        // field reads t (a belt's teeth running round its loop).
        let mut bodies: Vec<u32> = moving.iter().filter_map(|m| m.2.as_ref().map(|(s, _)| *s)).collect();
        bodies.sort_unstable();
        bodies.dedup();
        let own: HashMap<u32, Option<(V3, V3)>> =
            v.at_point(over.clone(), |v| bodies.iter().copied().zip(found(v, &bodies)).map(|(s, (w, b))| (s, b.filter(|_| w).map(|b| (b.0, b.1)))).collect());
        let (g, ctx) = (&v.doc.g, &v.ctx);
        let pieces: Vec<Option<Vec<((f64, f64), (V3, V3))>>> =
            moving.par_iter().map(|(_, _, sweep)| sweep.as_ref().and_then(|(s, sweep)| swept(g, ctx, sweep, own[s]?, span, true))).collect();
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

/// A placement p = R q + T a body's box is swept by (`swept`): R's and
/// T's entries, and how far its solid may lie from where they move it.
/// The field it places reads R^T (p - T), whose solid is T + R^-T q =
/// T + R (q + d), |d| at most `slack` max |q| (TurnBounds::slack): the
/// body's box widened by that is swept. 0 for an isometry.
pub(super) struct Sweep {
    pub r: [[u32; 3]; 3],
    pub t: [u32; 3],
    pub slack: f64,
}

impl From<&Isometry> for Sweep {
    fn from(iso: &Isometry) -> Self {
        Sweep { r: *iso.r(), t: *iso.t(), slack: 0.0 }
    }
}

/// Item `it`'s placement `at` swept over the span (`over`): proven an
/// isometry (rigid::Isometry), else proven near one a factor at a time.
/// Where `at` is its links composed (Item::links, as Doc::items composes
/// them: the same nodes), and a joint's link its anchor, motion and the
/// placement below it composed (Doc::instance_placement: the same nodes),
/// R is their turns' product: with every factor's entries finite over the
/// span and its turn read alone (turn_bounds; each once, `turns`), its
/// slack is theirs composed (TurnBounds::composed). So a turn about an
/// axis tilted by floats, under other turns, whose composed entries are
/// too many terms to read or whose angles do not cancel, is as near a
/// rotation as its factors are. None: neither is proven.
fn sweep(v: &mut ModelView, it: &crate::Item, at: &Xf, over: &Ctx, turns: &mut HashMap<Xf, Option<TurnBounds>>) -> Option<Sweep> {
    if let Some(iso) = Isometry::over(&v.doc.g, &at.r, &at.t, over) {
        return Some(Sweep::from(&iso));
    }
    let id = v.doc.identity();
    if it.links.is_empty() || it.links.iter().fold(id, |x, l| v.doc.compose(&l.at, &x)) != *at {
        return None;
    }
    let mut factors = Vec::with_capacity(3 * it.links.len());
    for l in &it.links {
        match &l.joint {
            Some(j) => {
                let moved = v.doc.compose(&j.below, &j.motion);
                if v.doc.compose(&moved, &j.anchor) != l.at {
                    return None;
                }
                factors.extend([j.anchor, j.motion, j.below]);
            }
            None => factors.push(l.at),
        }
    }
    let g = &v.doc.g;
    let read = |f: &Xf| ink_kernel::rigid::finite_over(g, f.r.iter().flatten().chain(&f.t), over).then(|| ink_kernel::rigid::turn_bounds(g, &f.r, over)).flatten();
    let bounds: Vec<TurnBounds> = factors.iter().map(|f| *turns.entry(*f).or_insert_with(|| read(f))).collect::<Option<_>>()?;
    Some(Sweep { r: at.r, t: at.t, slack: TurnBounds::composed(&bounds).slack()? })
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
/// `at` (rigid over the span: its solid is where R q + T puts its box,
/// widened by `at.slack`) over SWEPT_PIECES uniform pieces of `span`, each
/// holding it at every moment of its piece: R q + T with R's and T's
/// entries bounded over the piece. None: an entry not bounded.
pub(super) fn swept(g: &ink_core::graph::Graph, ctx: &Ctx, at: &Sweep, (lo, hi): (V3, V3), span: (f64, f64), refine: bool) -> Option<Vec<((f64, f64), (V3, V3))>> {
    use interval::I;
    let entry = |n: u32| Tape::new(g, n);
    let r: Vec<Tape> = at.r.iter().flatten().map(|&n| entry(n)).collect();
    let t: Vec<Tape> = at.t.iter().map(|&n| entry(n)).collect();
    let d = interval::mul_up(at.slack, lo.into_iter().chain(hi).map(f64::abs).fold(0.0, f64::max));
    let (lo, hi) = match d == 0.0 {
        true => (lo, hi),
        false => (lo.map(|x| interval::sub_o(I::point(x), I::point(d)).lo), hi.map(|x| interval::add_o(I::point(x), I::point(d)).hi)),
    };
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
