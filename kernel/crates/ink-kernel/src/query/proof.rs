//! Verdicts and the shared tests and sampling for proof searches.

use crate::search;
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_interval};
use ink_core::interval::I;
use ink_core::tape::{BOTH, Tape};
use std::sync::Arc;

// The instructions evaluated on this thread by `extent`'s tests and covers
// and by an emptiness search's reads of a and b (`Except`), and the
// longest tape among them: the tests' measure of their work.
#[cfg(test)]
thread_local!(pub(super) static WORK: std::cell::Cell<(usize, usize)> = const { std::cell::Cell::new((0, 0)) });

/// Count an evaluation of `tape` toward `WORK` (tests only).
pub(super) fn worked(tape: &Tape) {
    worked_len(tape.len());
}

/// `worked` for a tape of `_n` instructions.
pub(super) fn worked_len(_n: usize) {
    #[cfg(test)]
    WORK.with(|w| w.set((w.get().0 + _n, w.get().1.max(_n))));
}

#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Verdict {
    Proven,
    /// A point inside, at time t and the ranged params' values (none when
    /// no param is ranged), with the field's value there.
    Disproven {
        witness: V3,
        t: f64,
        value: f64,
        params: Vec<(u32, f64)>,
    },
    /// Boxes left undecided and their volume, and one of them, (lo, hi):
    /// the one that ended the proof (stuck at the finest cells), else the
    /// first left in its order (none where it had none); why the proof
    /// ended (`ended`); the hunt for a violation after it, none found, and
    /// why it stopped (`hunt`: search::Search::hunt); and the accuracy it holds to
    /// (`attained`: nothing violates it deeper than that, from the boxes
    /// left: no point of an emptiness's solid deeper; a clearance d, no
    /// pair closer than d less it; None: none proven, or the time limit
    /// stopped it).
    Unknown {
        cells: usize,
        volume: f64,
        ended: Ended,
        cell: Option<[V3; 2]>,
        attained: Option<f64>,
        hunt: Option<search::Hunt>,
    },
}

/// Why a search ended undecided.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Ended {
    /// A box undecided with every side at its finest cells (`cell`), the
    /// first the search reached: the claim cannot be proven at them.
    Finest,
    /// It tested the boxes it may (BOXES, or the budget asked), or would
    /// have tested over four times as many, or past an eighth of them over
    /// twice as many (search::Search::round).
    Limit,
    /// The call's time limit passed (crate::stop): the boxes are those it
    /// left.
    Stopped,
}

/// Confirms candidate witnesses: whether a tape is below a level at a point
/// (at time ctx.t) by bounds rounded outward. Searches find witnesses with
/// the f64 evaluator, whose value can be a rounding below the level where
/// the real value is above it (x x - (2 x - 1) rounds to 0 at x = 1 + 2^-27,
/// really 2^-54); a witness counts only once this confirms it. ctx.t must
/// lie in the span a simplified tape was simplified over. Made once per leaf
/// or search: the point context and the evaluator's storage are reused.
pub(super) struct Confirm {
    at: Ctx,
    in_span: bool,
    slots: Vec<I>,
}

impl Confirm {
    /// At ctx's point (ctx.t, ctx.params), which must lie in its spans; a
    /// ranged param keeps the meaning it has in the proof (Ctx::point).
    pub(super) fn new(ctx: &Ctx) -> Confirm {
        Confirm { at: ctx.point(), in_span: ctx.in_scope(), slots: Vec::new() }
    }
    pub(super) fn below(&mut self, tape: &Tape, p: V3, level: f64) -> bool {
        self.in_span && ink_core::eval::eval_interval_in(tape, &self.at, p.map(I::point), None, &mut self.slots).hi < level
    }
    /// Every point of the box is inside (below 0).
    pub(super) fn inside(&mut self, tape: &Tape, b: [I; 3]) -> bool {
        self.in_span && ink_core::eval::eval_interval_in(tape, &self.at, b, None, &mut self.slots).all_inside()
    }
}

/// `ctx` with its moment and ranged params inside their spans: a witness
/// carries ctx.t and ctx.params, which must be in what the proof covered.
pub(super) fn in_span(ctx: &Ctx) -> Ctx {
    ctx.into_scope()
}

/// The tape over box b for its sign alone (Tape::sign_choices): below, at
/// and above 0 where `t` is, there; the same where no choice is made.
pub(crate) fn signed(t: &Arc<Tape>, ctx: &Ctx, b: [I; 3]) -> Arc<Tape> {
    let mut vals = Vec::new();
    let mut ch = vec![BOTH; t.len()];
    match ink_core::eval::eval_interval_clamped(t, ctx, b, &[], &mut vals) {
        Some(_) if t.sign_choices(&vals, &mut ch) => t.narrowed(&ch),
        _ => t.clone(),
    }
}

/// Each input's blur of the bounds `blurs` gives (eval_interval_scope_blur:
/// every input of the scope a bound reads), each over its scale (a slope
/// bound, to make it a length; 1, as it is), the largest: each ranged
/// param's, sorted by param, and time's. A derived param's (a solve's
/// unknown) counts for those it is derived from, which are what halve.
pub(super) fn scope_shares(ctx: &Ctx, blurs: &[(&Vec<(ink_core::graph::Var, f64)>, f64)]) -> (Vec<(u32, f64)>, f64) {
    use ink_core::graph::Var;
    let (mut shares, mut t): (Vec<(u32, f64)>, f64) = (Vec::new(), 0.0);
    for &(blur, l) in blurs {
        for &(v, x) in blur.iter() {
            let x = if x == 0.0 {
                0.0
            } else if l.is_finite() && l > 0.0 {
                x / l
            } else {
                f64::INFINITY
            };
            match v {
                // NaN: as large as any.
                Var::T => t = if x.is_nan() || t.is_nan() { f64::NAN } else { t.max(x) },
                Var::Param(j) if ctx.p_spans.is_derived(j) => shares.extend(ctx.p_spans.sources(j).into_iter().map(|k| (k, x))),
                Var::Param(j) => shares.push((j, x)),
                _ => {}
            }
        }
    }
    shares.sort_unstable_by(|p, q| p.0.cmp(&q.0).then(q.1.total_cmp(&p.1)));
    shares.dedup_by_key(|x| x.0);
    (shares, t)
}

/// Param k's share (`scope_shares`): none where no bound reads it.
pub(super) fn share(shares: &[(u32, f64)], k: u32) -> f64 {
    shares.binary_search_by_key(&k, |x| x.0).map_or(0.0, |i| shares[i].1)
}

/// `k^3` sample points at the centres of a box's sub-boxes, `k` along each
/// side.
pub(crate) fn box_samples(bx: [I; 3], k: usize) -> Vec<V3> {
    let s: V3 = std::array::from_fn(|j| (bx[j].hi - bx[j].lo) / k as f64);
    let mut out = Vec::with_capacity(k * k * k);
    for a in 0..k {
        for b in 0..k {
            for c in 0..k {
                out.push([bx[0].lo + s[0] * (c as f64 + 0.5), bx[1].lo + s[1] * (b as f64 + 0.5), bx[2].lo + s[2] * (a as f64 + 0.5)]);
            }
        }
    }
    out
}

/// The moment a cached answer's witness is at: the one asked (in_span)
/// where no tape of the question reads t, its cache key leaving the moment
/// out (the answer is the same at any); None where one does: the key holds
/// the moment and its span, and the witness keeps its own (a search over a
/// span finds it at a moment of it).
pub fn asked_moment(ctx: &Ctx, tapes: &[&Tape]) -> Option<f64> {
    (!ink_core::tape::Names::of(tapes, ctx).reads_t).then(|| in_span(ctx).t)
}

/// A verdict's witness at the moment `asked_moment` gives.
pub(crate) fn at_moment(mut v: Verdict, ctx: &Ctx, tapes: &[&Tape]) -> Verdict {
    if let Verdict::Disproven { t, .. } = &mut v
        && let Some(m) = asked_moment(ctx, tapes)
    {
        *t = m;
    }
    v
}

/// A point of the solid {tape < 0} at p or just under its surface there,
/// confirmed by bounds rounded outward: p itself, else p stepped along the
/// field's descent to a depth of 1e-12 `scale` (a length), growing tenfold
/// to 1e-6 `scale`. None: no point is confirmed (no gradient at p, or the
/// field is not below 0 just past it).
pub fn confirm_inside(tape: &Tape, ctx: &Ctx, p: V3, scale: f64) -> Option<V3> {
    let mut confirm = Confirm::new(ctx);
    let ok = |c: &mut Confirm, x: V3| x.iter().all(|v| v.is_finite()) && c.below(tape, x, 0.0);
    if ok(&mut confirm, p) {
        return Some(p);
    }
    let g = eval_grads(tape, ctx, &[p])[0];
    let n = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
    if !(n > 0.0 && n.is_finite()) {
        return None;
    }
    let mut depth = 1e-12 * scale;
    for _ in 0..7 {
        // Newton's step onto the surface, then depth further along -grad.
        let k = (g[0] / n + depth) / n;
        let x = [p[0] - k * g[1], p[1] - k * g[2], p[2] - k * g[3]];
        if ok(&mut confirm, x) {
            return Some(x);
        }
        depth *= 10.0;
    }
    None
}

/// Box b narrowed along each axis to where the tape can be at most 0 by
/// its terms reading that axis alone (Tape::conjuncts: a slab's span), each
/// end found to 1/256 of the side by halving; None where it is nowhere.
/// The tape is above 0 at every point of b past the box returned.
pub(super) fn contract(t: &Tape, ctx: &Ctx, b: [I; 3]) -> Option<[I; 3]> {
    contract_to(t, ctx, b, [8; 3])
}

/// `contract`, each end along axis k found to 2^-depths[k] of the side (0:
/// that axis kept as it is).
pub(super) fn contract_to(t: &Tape, ctx: &Ctx, b: [I; 3], depths: [u32; 3]) -> Option<[I; 3]> {
    // The first end where `may` holds in x, to `depth` halvings, from the
    // low end (`up`) or the high.
    fn end(may: &dyn Fn(I) -> bool, x: I, depth: u32, up: bool) -> Option<f64> {
        if !may(x) {
            return None;
        }
        if depth == 0 || !(x.lo < x.mid() && x.mid() < x.hi) {
            return Some(if up { x.lo } else { x.hi });
        }
        let (l, h) = (I::new(x.lo, x.mid()), I::new(x.mid(), x.hi));
        let (first, then) = if up { (l, h) } else { (h, l) };
        end(may, first, depth - 1, up).or_else(|| end(may, then, depth - 1, up))
    }
    let mut b = b;
    for (i, below, axes) in t.conjuncts() {
        let Some(k) = (0..3).find(|&k| axes == 1 << k && depths[k] > 0) else { continue };
        let term = t.cone(i);
        // Whether the term can be at most 0 (at least 0, not `below`) with
        // axis k in x (NaN: it can).
        let may = |x: I| {
            let mut c = b;
            c[k] = x;
            let v = eval_interval(&term, ctx, c, None);
            if below { !(v.lo > 0.0) } else { !(v.hi < 0.0) }
        };
        let lo = end(&may, b[k], depths[k], true)?;
        let hi = end(&may, b[k], depths[k], false)?;
        b[k] = I::new(lo, hi);
    }
    Some(b)
}
