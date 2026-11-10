//! Integrals over material (volume, mass, mass moments, a claim's
//! integrals): proven bounds, and an estimate within them, from the one
//! search (search::Search) over boxes of the body's box.
//!
//! Each box is bounded on its own. One the material fills: each
//! integrand's range times its volume, or tighter, its centred form
//! (exact for an integrand linear across it). One the surface crosses: by
//! the mean value theorem f(x) = f(c) + sum_k a_k (x_k - c_k), with c its
//! centre and each a_k in the range of the field's slope along x_k over the
//! box (which covers both sides of a kink), so the material holds the
//! points where the largest such sum is below 0 and lies within those where
//! the least is: in each octant about c those are a box cut by a plane,
//! whose volumes have a closed form (`cut`). A flat face is bounded
//! exactly in a box of any size, so a box need not fit inside a wall:
//! once halved until it holds one face of it, its share is known. Boxes
//! wait widest bounds first and are halved on the axis the field changes
//! most along (a plate's across its thickness, not along it), until the
//! first integral's bounds are within the accuracy asked of it.

mod cells;
pub use cells::PASS_CELLS;

use crate::search::{self, Node, Precision, Smear, Test};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_interval, eval_interval_grad, eval_points};
use ink_core::graph::{BinOp, UnOp};
use ink_core::interval::{I, add_o, mag, mul_o, sub_o, un};
use ink_core::tape::{BOTH, FIRST, Op, SECOND, Tape};
use rayon::prelude::*;
use std::collections::HashSet;
use std::sync::Arc;

/// How far, relative to it, a mass, a volume or a published mass's
/// density is measured to: its proven bounds within this of its estimate.
pub const ACCURACY: f64 = 0.02;

/// Boxes a measurement may test: one that would test more stops there,
/// unsettled, with the bounds it reached (a few seconds' work).
pub const BOXES: usize = 1 << 20;

/// How far a measurement goes: until the first integral's bounds are
/// within this share of its estimate either way, or with every box halved to sides of
/// at most this (mm) where its bounds still differ.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Until {
    Accuracy(f64),
    Res(f64),
}

/// Integrals as measured: each one's estimate (`sums`) and proven bounds
/// (`lo`, `hi`; the estimate lies within them), the smallest side of a box
/// tested (`res`), the boxes tested, whether it reached what it was asked
/// (`settled`; else why not, `why`), where its bounds were out of reach,
/// its estimate's last move (`sampled`: twice that move, relative to it;
/// `estimate`), which its spread is then, and how far each estimate is
/// known to (`moved`: twice what its last refinement moved it; 0 where
/// none was made).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Measured {
    pub sums: Vec<f64>,
    pub lo: Vec<f64>,
    pub hi: Vec<f64>,
    pub res: f64,
    pub boxes: usize,
    pub settled: bool,
    pub why: Option<String>,
    pub sampled: Option<f64>,
    pub moved: Vec<f64>,
}

impl Measured {
    /// How far integral k can be from its estimate: to the farther of its
    /// bounds.
    pub fn error(&self, k: usize) -> f64 {
        (self.sums[k] - self.lo[k]).max(self.hi[k] - self.sums[k])
    }

    /// Of `measured` with moments: rho * [1, x, y, z, xx, yy, zz, xy, yz,
    /// zx], then (the second integrand's, where it is 1) the volume.
    pub fn moments(&self) -> [f64; 10] {
        std::array::from_fn(|k| self.sums[k])
    }
    pub fn volume(&self) -> f64 {
        self.sums[10]
    }

    /// How far the first integral can be from its estimate, relative to
    /// it: proven, or where its bounds were out of reach, its sampled
    /// spread (`sampled`; None: it is 0).
    pub fn relative_error(&self) -> Option<f64> {
        self.sums.first().filter(|&&v| v != 0.0).map(|&v| self.sampled.unwrap_or(self.error(0) / v.abs()))
    }

    /// Why it did not settle (None: it did).
    pub fn unsettled(&self) -> Option<String> {
        (!self.settled).then(|| self.why.clone().unwrap_or_default())
    }
}

/// `integrate`, cached as the kernel's other answers are.
pub fn measured(shape: &Arc<Tape>, gs: &[Option<&Tape>], moments: bool, ctx: &Ctx, lo: V3, hi: V3, until: Until) -> Measured {
    let tapes: Vec<&Tape> = std::iter::once(&**shape).chain(gs.iter().flatten().copied()).collect();
    crate::cache::memo_region(
        "measured",
        &tapes,
        ctx,
        (lo, hi),
        |fp| {
            let (kind, v) = match until {
                Until::Accuracy(a) => (0, a),
                Until::Res(r) => (1, r),
            };
            fp.word(kind).word(v.to_bits()).word(gs.len() as u64).word(moments as u64);
            gs.iter().for_each(|g| _ = fp.word(g.is_some() as u64));
        },
        |ctx| integrate(shape, gs, moments, ctx, lo, hi, until),
    )
}

/// The terms of mass moments: [1, x, y, z, xx, yy, zz, xy, yz, zx].
pub(crate) const MOMENTS: [fn(V3) -> f64; 10] =
    [|_| 1.0, |p| p[0], |p| p[1], |p| p[2], |p| p[0] * p[0], |p| p[1] * p[1], |p| p[2] * p[2], |p| p[0] * p[1], |p| p[1] * p[2], |p| p[2] * p[0]];

/// The moments' terms' ranges over a box, rounded outward.
fn moment_ranges(b: [I; 3]) -> [I; 10] {
    let sq = |x| un(UnOp::Square, x);
    [I::point(1.0), b[0], b[1], b[2], sq(b[0]), sq(b[1]), sq(b[2]), mul_o(b[0], b[1]), mul_o(b[1], b[2]), mul_o(b[2], b[0])]
}

/// A box waiting: its tape, simplified over it, and its share of each
/// integral (estimates, then lower bounds, then upper), which its halves
/// replace.
#[derive(Clone)]
struct Held {
    tape: Arc<Tape>,
    share: Arc<[f64]>,
}

/// The question: integrands `gs` (None: 1) over shape < 0, the first's
/// times the moments' terms where `moments`.
struct Integrals<'a> {
    shape: &'a Arc<Tape>,
    gs: &'a [Option<&'a Tape>],
    moments: bool,
    ctx: &'a Ctx,
    /// Integrals in all.
    n: usize,
}

/// A box tested: its share (empty: none), whether its bounds are as
/// tight as rounding leaves them, each axis's smear, and its tape.
struct Tested {
    share: Vec<f64>,
    settled: bool,
    smear: [f64; 3],
    tape: Arc<Tape>,
}

impl Integrals<'_> {
    fn terms(&self, i: usize) -> &'static [fn(V3) -> f64] {
        if i == 0 && self.moments { &MOMENTS } else { &MOMENTS[..1] }
    }

    /// Each integral from its integrand's values `v` (each at `pts`, the
    /// first's times the moments' terms) weighted `w`.
    fn at(&self, i: usize, pts: &[V3], v: &[f64], w: &[f64], out: &mut Vec<f64>) {
        self.spread_at(i, pts, v, w, &[], out)
    }

    /// `at`, each point's material spread about it along each axis by
    /// `var` (its second moment there per unit, none where empty): the
    /// squared terms (xx, yy, zz) take it.
    fn spread_at(&self, i: usize, pts: &[V3], v: &[f64], w: &[f64], var: &[V3], out: &mut Vec<f64>) {
        for (j, term) in self.terms(i).iter().enumerate() {
            let own = |p: usize| if (4..7).contains(&j) && !var.is_empty() { var[p][j - 4] } else { 0.0 };
            out.push(v.iter().zip(pts).zip(w).enumerate().map(|(p, ((a, x), b))| a * (term(*x) + own(p)) * b).sum());
        }
    }

    /// A box's share of each integral, with bounds, `vol` the bounds of
    /// the material's volume in it: `full`, the material fills it. In a
    /// full box an integrand that jumps across one select's condition (a
    /// print's walls and infill) is each side's integral over the box
    /// where the condition puts it, bounded as the material's surface is
    /// (`slope_volume`) and estimated as it is.
    /// A box the surface crosses estimated again from `k` samples along
    /// each axis (None: one the material fills, whose estimate stands): a
    /// plane's cut misses a curved surface by its curvature times a
    /// sample's side squared, from below where the surface is convex (a
    /// hole's), so at twice the samples a quarter as much.
    fn finer(&self, tape: &Arc<Tape>, b: [I; 3], k: [usize; 3]) -> Option<Vec<f64>> {
        if eval_interval(tape, self.ctx, b, None).all_inside() {
            return None;
        }
        let sides: V3 = std::array::from_fn(|a| (b[a].hi - b[a].lo) / k[a] as f64);
        let mut pts = Vec::with_capacity(k[0] * k[1] * k[2]);
        for i in 0..k[2] {
            for j in 0..k[1] {
                for l in 0..k[0] {
                    pts.push([b[0].lo + sides[0] * (l as f64 + 0.5), b[1].lo + sides[1] * (j as f64 + 0.5), b[2].lo + sides[2] * (i as f64 + 0.5)]);
                }
            }
        }
        let wv = sides.iter().product::<f64>();
        let (mut w, mut at, mut var) = (Vec::new(), Vec::new(), Vec::new());
        for ((f, d, v), p) in shares(tape, self.ctx, &pts, sides, 3).into_iter().zip(&pts) {
            w.push(f * wv);
            at.push(std::array::from_fn(|k| p[k] + d[k]));
            var.push(v);
        }
        let mut est = Vec::with_capacity(self.n);
        let mut v = Vec::new();
        for (i, g) in self.gs.iter().enumerate() {
            v.clear();
            match g {
                Some(t) => eval_points(t, self.ctx, &pts, &mut v),
                None => v.resize(pts.len(), 1.0),
            }
            self.spread_at(i, &at, &v, &w, &var, &mut est);
        }
        Some(est)
    }

    fn share(&self, tape: &Arc<Tape>, b: [I; 3], full: bool, vol: I) -> (Vec<f64>, bool, [f64; 3]) {
        let n = self.n;
        let mut s = vec![0.0; 3 * n];
        let pts = crate::query::box_samples(b, 2);
        let wv = b.iter().map(|x| x.hi - x.lo).product::<f64>() / 8.0;
        let sides = b.map(|x| (x.hi - x.lo) / 2.0);
        // Each sample weighted by the material in its eighth of the box,
        // the field taken as linear there.
        // Each sample's material, and where in its eighth it lies.
        let fill = |t: &Tape| -> (Vec<f64>, Vec<V3>, Vec<V3>) {
            let (mut w, mut at, mut var) = (Vec::new(), Vec::new(), Vec::new());
            for ((f, d, v), p) in shares(t, self.ctx, &pts, sides, 3).into_iter().zip(&pts) {
                w.push(f * wv);
                at.push(std::array::from_fn(|k| p[k] + d[k]));
                var.push(v);
            }
            (w, at, var)
        };
        let values = |t: Option<&Tape>, at: &[V3]| -> Vec<f64> {
            let mut v = Vec::new();
            match t {
                Some(t) => eval_points(t, self.ctx, at, &mut v),
                None => v.resize(at.len(), 1.0),
            }
            v
        };
        // Gauss's 2x2x2 points: exact for an integrand of degree 3.
        let c = b.map(|x| x.mid());
        let g = 0.5 / 3f64.sqrt();
        let gp: Vec<V3> = (0..8).map(|k| std::array::from_fn(|a| c[a] + if (k >> a) & 1 == 1 { g } else { -g } * (b[a].hi - b[a].lo))).collect();
        let surface = (!full).then(|| fill(tape));
        let terms = moment_ranges(b);
        let mut tight = true;
        let mut smear = [0.0; 3];
        let mut est = Vec::with_capacity(n);
        let mut k = 0;
        for (i, g) in self.gs.iter().enumerate() {
            let range = g.map_or(I::point(1.0), |t| eval_interval(t, self.ctx, b, None));
            // The integrand's own integral's bounds, and the estimates.
            let first = match (g, &surface) {
                (_, Some((w, at, var))) => {
                    self.spread_at(i, at, &values(*g, &pts), w, var, &mut est);
                    mul_o(range, vol)
                }
                (None, None) => {
                    self.at(i, &gp, &[1.0; 8], &[wv; 8], &mut est);
                    vol
                }
                (Some(t), None) => match jump(t, self.ctx, b) {
                    None => {
                        self.at(i, &gp, &values(Some(t), &gp), &[wv; 8], &mut est);
                        let (centred, d) = centred(t, self.ctx, b, vol);
                        let rel = mag(range).max(f64::MIN_POSITIVE);
                        for a in 0..3 {
                            smear[a] += mag(d[a]) * (b[a].hi - b[a].lo) * 0.5 / rel;
                        }
                        let r = mul_o(range, vol);
                        I::new(r.lo.max(centred.lo), r.hi.min(centred.hi))
                    }
                    Some(j) => {
                        // On the second side throughout, and the first's
                        // difference where the condition is below 0.
                        let (va, vb, frac) = (values(Some(&j.a), &pts), values(Some(&j.b), &pts), fill(&j.cond).0);
                        let v: Vec<f64> = (0..8).map(|p| vb[p] + (va[p] - vb[p]) * frac[p] / wv).collect();
                        self.at(i, &pts, &v, &[wv; 8], &mut est);
                        let (rb, (centred, _)) = (eval_interval(&j.b, self.ctx, b, None), centred(&j.b, self.ctx, b, vol));
                        let rb = mul_o(rb, vol);
                        let on_b = I::new(rb.lo.max(centred.lo), rb.hi.min(centred.hi));
                        let diff = sub_o(eval_interval(&j.a, self.ctx, b, None), eval_interval(&j.b, self.ctx, b, None));
                        let dc = eval_interval_grad(&j.cond, self.ctx, b).d;
                        let fc = eval_interval(&j.cond, self.ctx, c.map(I::point), None);
                        let (vl, vh) = slope_volume(fc, dc, b);
                        let vh = vh.min(vol.hi);
                        let below = I::new(vl.max(0.0).min(vh), vh);
                        let rel = mag(eval_interval(&j.cond, self.ctx, b, None)).max(f64::MIN_POSITIVE);
                        for a in 0..3 {
                            let m = mag(dc[a]) * (b[a].hi - b[a].lo) * 0.5 / rel;
                            smear[a] += if m.is_finite() { m } else { f64::NAN };
                        }
                        let r = add_o(on_b, mul_o(diff, below));
                        let plain = mul_o(range, vol);
                        I::new(r.lo.max(plain.lo), r.hi.min(plain.hi))
                    }
                },
            };
            // Settled where rounding is all that is left of its bounds.
            tight &= first.hi - first.lo <= 1e-12 * mag(first);
            for j in 0..self.terms(i).len() {
                let r = if j == 0 { first } else { mul_o(mul_o(range, terms[j]), vol) };
                s[k] = est[k].clamp(r.lo, r.hi);
                (s[n + k], s[2 * n + k]) = (r.lo, r.hi);
                k += 1;
            }
        }
        (s, tight, smear)
    }

    /// A box tested on `tape`: no material (None), or its share.
    fn test(&self, tape: &Arc<Tape>, b: [I; 3]) -> Option<Tested> {
        let mut ch = vec![BOTH; tape.len()];
        let r = eval_interval(tape, self.ctx, b, Some(&mut ch));
        if r.none_inside() {
            return None;
        }
        let side = b.map(|x| sub_o(I::point(x.hi), I::point(x.lo)));
        let vol = mul_o(mul_o(side[0], side[1]), side[2]);
        if r.all_inside() {
            let (share, settled, smear) = self.share(tape, b, true, vol);
            return Some(Tested { share, settled, smear, tape: tape.clone() });
        }
        let sub = if ch.iter().any(|&c| c != BOTH) { Arc::new(tape.simplify(&ch)) } else { tape.clone() };
        let d = eval_interval_grad(&sub, self.ctx, b).d;
        let c = b.map(|x| x.mid());
        let fc = eval_interval(&sub, self.ctx, c.map(I::point), None);
        let (vlo, vhi) = slope_volume(fc, d, b);
        let vh = vhi.min(vol.hi);
        let v = I::new(vlo.max(0.0).min(vh), vh);
        let (share, settled, _) = self.share(&sub, b, false, v);
        // The field's change along each axis over the box: NaN, where its
        // slope is unbounded (a jump), as large as any.
        let smear = std::array::from_fn(|k| {
            let m = mag(d[k]) * (b[k].hi - b[k].lo) * 0.5;
            if m.is_finite() { m } else { f64::NAN }
        });
        Some(Tested { share, settled, smear, tape: sub })
    }
}

/// An integrand that jumps in a box: across the one select there whose
/// condition its bounds leave undecided between branches that are not one
/// field. Its condition, and the integrand with the select on each branch.
struct Jump {
    cond: Tape,
    a: Tape,
    b: Tape,
}

/// `g`'s jump in box `b` (None: none, or more than one select undecided).
fn jump(g: &Tape, ctx: &Ctx, b: [I; 3]) -> Option<Jump> {
    let mut ch = vec![BOTH; g.len()];
    eval_interval(g, ctx, b, Some(&mut ch));
    let g = g.simplify(&ch);
    let strip = |mut j: usize| {
        while let Op::Tag(_) = g.ssa[j].op {
            j = g.ssa[j].a as usize;
        }
        j
    };
    let mut selects = (0..g.len()).filter(|&i| matches!(g.ssa[i].op, Op::Select) && strip(g.ssa[i].b as usize) != strip(g.ssa[i].c as usize));
    let (Some(i), None) = (selects.next(), selects.next()) else { return None };
    let mut pick = vec![BOTH; g.len()];
    pick[i] = FIRST;
    let a = g.simplify(&pick);
    pick[i] = SECOND;
    Some(Jump { cond: g.cone(g.ssa[i].a as usize), a, b: g.simplify(&pick) })
}

/// Bounds of the integral of `g` over a box full of material (volume
/// `vol`): its centred form, g(c) vol plus, per axis, the slope's range
/// times the integral of x_k - c_k where positive, less where negative
/// (exact for g linear in x, y, z: a centre of mass's moment), and the
/// slopes' ranges. A g
/// that jumps in the box (a select, a floor) has unbounded slopes there:
/// its plain bound is the tighter.
fn centred(g: &Tape, ctx: &Ctx, b: [I; 3], vol: I) -> (I, [I; 3]) {
    let side: [I; 3] = b.map(|x| sub_o(I::point(x.hi), I::point(x.lo)));
    let gg = eval_interval_grad(g, ctx, b);
    let centre: V3 = b.map(|x| x.mid());
    let mut centred = mul_o(vol, eval_interval(g, ctx, centre.map(I::point), None));
    for k in 0..3 {
        // The integral of x_k - c_k where positive and where negative: the
        // cross-section times half the squared distance to each face.
        let area = mul_o(side[(k + 1) % 3], side[(k + 2) % 3]);
        let half = |d: I| mul_o(mul_o(area, mul_o(d, d)), I::point(0.5));
        let up = half(sub_o(I::point(b[k].hi), I::point(centre[k])));
        let down = half(sub_o(I::point(centre[k]), I::point(b[k].lo)));
        centred = add_o(centred, sub_o(mul_o(gg.d[k], up), mul_o(gg.d[k], down)));
    }
    (centred, gg.d)
}

/// Bounds of the volume where f < 0 in box `b`, from f at its centre c
/// (`fc`, rounded outward) and the ranges of its slopes over the box
/// (`d`): f(x) lies between f(c) plus the least and the most of sum_k a_k
/// (x_k - c_k) over a_k in d_k, each linear in an octant about c, so the
/// material holds where the most is below 0, and lies where the least is.
/// Exact for a flat face, to rounding.
fn slope_volume(fc: I, d: [I; 3], b: [I; 3]) -> (f64, f64) {
    let c: V3 = b.map(|x| x.mid());
    let (mut lo, mut hi) = (I::point(0.0), I::point(0.0));
    for o in 0..8 {
        let up: [bool; 3] = std::array::from_fn(|k| (o >> k) & 1 == 1);
        // Along u_k = +-(x_k - c_k), from 0 to the box's face.
        let h: [I; 3] = std::array::from_fn(|k| if up[k] { sub_o(I::point(b[k].hi), I::point(c[k])) } else { sub_o(I::point(c[k]), I::point(b[k].lo)) });
        let least: [f64; 3] = std::array::from_fn(|k| if up[k] { d[k].lo } else { -d[k].hi });
        let most: [f64; 3] = std::array::from_fn(|k| if up[k] { d[k].hi } else { -d[k].lo });
        hi = add_o(hi, I::point(cut(fc.lo, least, h, true)));
        lo = add_o(lo, I::point(cut(fc.hi, most, h, false)));
    }
    (lo.lo, hi.hi)
}

/// Terms of a sum of uniforms (`uniform_sum_below`) under this share of
/// the largest are taken by their range: the formula's rounding grows as
/// the inverse square of the smallest term's share.
const TAU: f64 = 1e-3;

/// The formula's rounding (`uniform_sum_below`), at most 200 eps / TAU^2
/// of a probability, rounded up.
const ROUNDING: f64 = 1e-7;

/// The volume of the points u of [0, h_0] x [0, h_1] x [0, h_2] where
/// r + sum_k a_k u_k < 0: above it where `upper`, else below it (the
/// sides h, rounded outward, at their larger or smaller ends).
fn cut(r: f64, a: [f64; 3], h: [I; 3], upper: bool) -> f64 {
    let h: [f64; 3] = h.map(|x| if upper { x.hi } else { x.lo.max(0.0) });
    let (pick, nudge) = if upper { ((|x: I| x.hi) as fn(I) -> f64, 1.0 + 4.0 * f64::EPSILON) } else { ((|x: I| x.lo) as fn(I) -> f64, 1.0 - 4.0 * f64::EPSILON) };
    let full = h[0] * h[1] * h[2] * nudge;
    if !(full > 0.0) {
        return 0.0;
    }
    if !r.is_finite() || a.iter().any(|x| !x.is_finite()) {
        return if upper { full } else { 0.0 };
    }
    // The region grows as r and each b_k = a_k h_k fall: an upper bound
    // takes them low. With u_k = h_k Y_k, it is sum_k b_k Y_k < -r, Y in
    // [0, 1]^3; Y_k -> 1 - Y_k where b_k < 0 makes every term
    // nonnegative: sum_k |b_k| Y_k < -r - sum of the negative b_k.
    let low = |x: I| if upper { x.lo } else { x.hi };
    let b: [f64; 3] = std::array::from_fn(|k| low(mul_o(I::point(a[k]), I::point(h[k]))));
    let mut s = I::point(-r);
    for &bk in b.iter().filter(|&&x| x < 0.0) {
        s = sub_o(s, I::point(bk));
    }
    let p = uniform_sum_below(b.map(f64::abs), pick(s), upper);
    pick(mul_o(I::point(p), I::point(full)))
}

/// The probability that sum_k c_k Y_k < s, the Y_k uniform on [0, 1] and
/// independent, c_k >= 0: above it where `upper`, else below it. A term
/// under TAU of the largest is taken by its range (an upper bound leaves
/// it out; a lower takes it at its largest). The rest by inclusion and
/// exclusion over the corners of the box of Y, in units of the largest
/// term, from the nearer end (the
/// sum is symmetric about its middle), within ROUNDING.
fn uniform_sum_below(c: [f64; 3], s: f64, upper: bool) -> f64 {
    let top = c.iter().copied().fold(0.0, f64::max);
    let mut s = I::point(s);
    let mut kept: Vec<f64> = Vec::with_capacity(3);
    for &x in &c {
        if x > TAU * top {
            kept.push(x);
        } else if !upper {
            s = sub_o(s, I::point(x));
        }
    }
    let s = if upper { s.hi } else { s.lo };
    if !(s > 0.0) {
        return 0.0;
    }
    let total = kept.iter().sum::<f64>();
    if s >= total * (1.0 + 8.0 * f64::EPSILON) {
        return 1.0;
    }
    let m = kept.len();
    if m == 0 {
        return 1.0;
    }
    // In units of the largest term, so the powers below neither underflow
    // nor overflow (terms of 1e-108 cube to subnormals): each quotient
    // rounded the way that widens the bound (the region grows as the
    // terms fall and s rises).
    let (up, down): (fn(f64) -> f64, fn(f64) -> f64) = (f64::next_up, f64::next_down);
    let (toward_terms, toward_s) = if upper { (down, up) } else { (up, down) };
    let unit = kept.iter().copied().fold(0.0, f64::max);
    let kept: Vec<f64> = kept.iter().map(|&x| toward_terms(x / unit)).collect();
    let (s, total) = (toward_s(s / unit), kept.iter().sum::<f64>());
    let cdf = |s: f64| {
        let mut acc = 0.0;
        for set in 0..1usize << m {
            let shift: f64 = (0..m).filter(|k| set >> k & 1 == 1).map(|k| kept[k]).sum();
            let x = s - shift;
            if x > 0.0 {
                let t = x.powi(m as i32);
                acc += if set.count_ones() % 2 == 0 { t } else { -t };
            }
        }
        let fact = [1.0, 1.0, 2.0, 6.0][m];
        acc / (fact * kept.iter().product::<f64>())
    };
    let p = if s > 0.5 * total { 1.0 - cdf(total - s) } else { cdf(s) };
    if upper { (p + ROUNDING).min(1.0) } else { (p - ROUNDING).max(0.0) }
}

/// The share of each of the boxes of sides `s` centred on `pts` that is
/// material, the field taken as linear there (`fraction`), and how far
/// from the centre that material's centroid lies (`centroid`); but where
/// the field is an intersection or a union (a max or min of branches: a
/// box's edge, a hole's rim), its branches' combined (`combined`), to
/// `depth` levels of choices: one plane's cut misreads a box's edge at
/// every size.
fn shares(t: &Tape, ctx: &Ctx, pts: &[V3], s: V3, depth: u32) -> Vec<(f64, V3, V3)> {
    // Through names and slope promises, which keep the value.
    let strip = |mut j: usize| {
        while let Op::Tag(_) | Op::Lipschitz = t.ssa[j].op {
            j = t.ssa[j].a as usize;
        }
        j
    };
    let ins = t.ssa[strip(t.len() - 1)];
    let leaf = || {
        eval_grads(t, ctx, pts)
            .iter()
            .map(|f| {
                let share = fraction(f, s);
                let (c, v) = centroid(f, s, share);
                (share, c, v)
            })
            .collect()
    };
    let intersection = match ins.op {
        Op::Bin(BinOp::Max | BinOp::Emax) if depth > 0 => true,
        Op::Bin(BinOp::Min) if depth > 0 => false,
        _ => return leaf(),
    };
    let (a, b) = (t.cone(ins.a as usize), t.cone(ins.b as usize));
    let (ga, gb) = (eval_grads(&a, ctx, pts), eval_grads(&b, ctx, pts));
    let (fa, fb) = (shares(&a, ctx, pts, s, depth - 1), shares(&b, ctx, pts, s, depth - 1));
    (0..pts.len()).map(|i| combined(fa[i], fb[i], &ga[i], &gb[i], intersection)).collect()
}

/// Where in a box of sides `s` the material a plane (the field's value and
/// gradient at its centre, `f`) leaves it, `share` of the box, lies, from
/// its centre, and its second moment about there per unit along each
/// axis: along each axis, the slab that share of the box makes across the
/// gradient's part along it (exact for a face across an axis).
fn centroid(f: &[f64; 4], s: V3, share: f64) -> (V3, V3) {
    let range: f64 = (0..3).map(|k| (f[k + 1] * s[k]).abs()).sum();
    if !(range > 0.0 && range.is_finite()) {
        return ([0.0; 3], s.map(|x| x * x / 12.0));
    }
    let part = |k: usize| (f[k + 1] * s[k]).abs() / range;
    let long = |k: usize| s[k] * (1.0 - (1.0 - share) * part(k));
    (std::array::from_fn(|k| -(1.0 - share) * 0.5 * s[k] * f[k + 1] * s[k] / range), std::array::from_fn(|k| long(k) * long(k) / 12.0))
}

/// The share of a box in both of two regions (`intersection`), or in
/// either, from each one's share and centroid (`fa`, `fb`) and its field's
/// gradient there: their product, as for faces across each other (a box's
/// edge, as cells see it exactly; their centroids' offsets add), but of
/// faces nearly facing each other (a slab), what they overlap by, and of
/// faces nearly alike, the lesser. A union is the complement of its
/// complements' intersection, its centroid their shares' mean; its spread
/// along an axis the wider, an intersection's the narrower.
fn combined((fa, ca, va): (f64, V3, V3), (fb, cb, vb): (f64, V3, V3), ga: &[f64; 4], gb: &[f64; 4], intersection: bool) -> (f64, V3, V3) {
    let (na, nb) = ((ga[1] * ga[1] + ga[2] * ga[2] + ga[3] * ga[3]).sqrt(), (gb[1] * gb[1] + gb[2] * gb[2] + gb[3] * gb[3]).sqrt());
    let cos = (ga[1] * gb[1] + ga[2] * gb[2] + ga[3] * gb[3]) / (na * nb);
    let both = |a: f64, b: f64| match cos {
        c if c < -0.5 => (a + b - 1.0).max(0.0),
        c if c > 0.5 => a.min(b),
        _ => a * b,
    };
    let f = if intersection { both(fa, fb) } else { 1.0 - both(1.0 - fa, 1.0 - fb) };
    let f = if f.is_finite() {
        f.clamp(0.0, 1.0)
    } else if intersection {
        fa.min(fb)
    } else {
        fa.max(fb)
    };
    let c: V3 = match intersection {
        true => std::array::from_fn(|k| ca[k] + cb[k]),
        false if fa + fb > 0.0 => std::array::from_fn(|k| (fa * ca[k] + fb * cb[k]) / (fa + fb)),
        false => [0.0; 3],
    };
    // Along each axis the material spans what both leave it, or either.
    let v: V3 = std::array::from_fn(|k| if intersection { va[k].min(vb[k]) } else { va[k].max(vb[k]) });
    (f, c, v)
}

/// The share of a box of sides `s` centred on a point that is material,
/// the field taken as linear there (`f`: its value and gradient): the
/// box cut by a plane.
fn fraction(f: &[f64; 4], s: V3) -> f64 {
    if !f.iter().all(|x| x.is_finite()) {
        return if f[0] < 0.0 { 1.0 } else { 0.0 };
    }
    // f + sum_k g_k s_k (Y_k - 1/2) < 0, Y uniform on [0, 1]^3.
    let b: V3 = std::array::from_fn(|k| f[k + 1] * s[k]);
    let start = -f[0] + 0.5 * b.iter().sum::<f64>() - b.iter().filter(|&&x| x < 0.0).sum::<f64>();
    let (up, down) = (uniform_sum_below(b.map(f64::abs), start, true), uniform_sum_below(b.map(f64::abs), start, false));
    0.5 * (up + down)
}

/// Integrals of each of `gs` (None: 1) over the region shape < 0 within
/// `lo..hi`, with `moments`, the first's times [1, x, y, z, xx, yy, zz,
/// xy, yz, zx] in its place (mass moments, rho the first): estimates and
/// proven bounds (the module's doc), measured `until` the first's bounds
/// are within an accuracy of it, or every box whose bounds differ is
/// halved to a size. Boxes wait widest bounds first (the integrals' widths,
/// each relative to its bounds over the whole box), BOXES at most, or
/// until the call's time limit; either leaves it unsettled.
pub fn integrate(shape: &Arc<Tape>, gs: &[Option<&Tape>], moments: bool, ctx: &Ctx, lo: V3, hi: V3, until: Until) -> Measured {
    integrate_within(shape, gs, moments, ctx, lo, hi, until, BOXES)
}

/// `integrate`, testing at most `boxes` boxes.
#[allow(clippy::too_many_arguments)]
fn integrate_within(shape: &Arc<Tape>, gs: &[Option<&Tape>], moments: bool, ctx: &Ctx, lo: V3, hi: V3, until: Until, boxes: usize) -> Measured {
    let n = gs.len() + if moments && !gs.is_empty() { 9 } else { 0 };
    let q = Integrals { shape, gs, moments, ctx, n };
    let b: [I; 3] = std::array::from_fn(|k| I::new(lo[k], hi[k]));
    let done = |sums: Vec<f64>, lo: Vec<f64>, hi: Vec<f64>, res: f64, boxes: usize, why: Option<String>, moved: Vec<f64>| Measured {
        moved,
        sampled: None,
        sums,
        lo,
        hi,
        res,
        boxes,
        settled: why.is_none(),
        why,
    };
    let side = crate::query::max_extent(lo, hi);
    let Some(root) = q.test(q.shape, b).filter(|_| n > 0) else { return done(vec![0.0; n], vec![0.0; n], vec![0.0; n], side, 1, None, vec![0.0; n]) };
    // Each integral's scale: its bounds over the whole box. The integrals
    // of the integrands themselves (not the moments' x, y, z terms, near 0
    // about the origin) order the boxes.
    let firsts: Vec<usize> = (0..gs.len()).map(|i| if moments && i > 0 { i + 9 } else { i }).collect();
    let scale: Vec<f64> = (0..n).map(|k| root.share[n + k].abs().max(root.share[2 * n + k].abs())).collect();
    let width = |s: &[f64]| firsts.iter().filter(|&&k| scale[k] > 0.0).map(|&k| (s[2 * n + k] - s[n + k]) / scale[k]).sum::<f64>();
    let (mut settled, mut open) = (vec![0.0; 3 * n], vec![0.0; 3 * n]);
    let add = |to: &mut [f64], s: &[f64], sign: f64| to.iter_mut().zip(s).for_each(|(t, x)| *t += sign * x);
    let ctx_arc = Arc::new(ctx.clone());
    let place = |t: Tested, open: &mut Vec<f64>, settled: &mut Vec<f64>| -> Test<Held, ()> {
        if t.settled {
            add(settled, &t.share, 1.0);
            return Test::Settled;
        }
        add(open, &t.share, 1.0);
        let bound = -width(&t.share);
        Test::Open { bound, smear: Smear { space: t.smear, ..Default::default() }, state: Held { tape: t.tape, share: t.share.into() } }
    };
    let first = place(root, &mut open, &mut settled);
    let finest = match until {
        Until::Res(r) => r,
        Until::Accuracy(_) => 0.0,
    };
    let mut s =
        search::Search::start(vec![(Node { b, ctx: ctx_arc, state: Held { tape: shape.clone(), share: Arc::from(vec![]) } }, first)], Precision { finest, splits: 0, maps: true });
    let reached = |open: &[f64], settled: &[f64]| match until {
        Until::Accuracy(a) => {
            let (est, l, h) = (open[0] + settled[0], open[n] + settled[n], open[2 * n] + settled[2 * n]);
            (est - l).max(h - est) <= a * est.abs()
        }
        Until::Res(_) => false,
    };
    let mut why = None;
    // At each doubling of the boxes tested: the first integral's error
    // and estimate.
    let (mut marks, mut mark): (Vec<(f64, f64)>, usize) = (Vec::new(), MARK);
    loop {
        if reached(&open, &settled) {
            break;
        }
        if let Until::Accuracy(a) = until
            && s.tested() >= mark
        {
            mark *= 2;
            let (est, l, h) = (open[0] + settled[0], open[n] + settled[n], open[2 * n] + settled[2 * n]);
            marks.push(((est - l).max(h - est), est));
            if out_of_reach(&marks, a, s.tested()) {
                // The bounds so far, from the settled boxes' and those
                // left's shares, the halved ones' not taken away.
                let mut total = settled.clone();
                s.left().for_each(|nd| add(&mut total, &nd.state.share, 1.0));
                let (l, h): (Vec<f64>, Vec<f64>) = (0..n).map(|k| widened((total[n + k], total[2 * n + k]), s.tested())).unzip();
                return cells::estimate(shape, gs, moments, ctx, lo, hi, a, (l, h));
            }
        }
        if s.tested() >= boxes {
            why = Some(format!("after {boxes} boxes"));
            break;
        }
        let most = (s.tested() / 16).clamp(search::BATCH, 1 << 16);
        let outcome = s.round_by(usize::MAX >> 4, most, |nodes| {
            // Each box halved gives its share up to its halves.
            let mut seen = HashSet::new();
            for nd in nodes {
                if seen.insert(Arc::as_ptr(&nd.state.share) as *const f64 as usize) {
                    add(&mut open, &nd.state.share, -1.0);
                }
            }
            let tested: Vec<Option<Tested>> = nodes.par_iter().map(|nd| q.test(&nd.state.tape, nd.b)).collect();
            tested.into_iter().map(|t| t.map_or(Test::Settled, |t| place(t, &mut open, &mut settled))).collect()
        });
        if let Some(o) = outcome {
            if let search::Outcome::Unknown(u) = o
                && u.stopped()
            {
                why = Some("the call's time limit passed".into());
            }
            break;
        }
    }
    // The totals again, from the settled boxes' and those left's shares,
    // added without the halved ones' taken away, and widened by the
    // rounding of a sum of as many terms.
    let mut total = settled;
    let mut terms = 1.0;
    // The boxes left the surface crosses, estimated again finer (their
    // bounds stand): where an accuracy was asked, at 4 samples a side;
    // where a res was, those still wider than it (the search stopped
    // short of it) at samples of half of it, as boxes of it take theirs,
    // or as near as the accuracy path's samples, 64 a box left, allow.
    let left: Vec<&Node<Held>> = s.left().collect();
    // Samples per axis of a box at `step` apart (at least the box's own 2),
    // counted in floats, which a fine res cannot overflow.
    let per_axis = |b: [I; 3], step: f64| b.map(|x| (2.0 * (x.hi - x.lo) / step).ceil().max(2.0));
    let ks: Vec<[usize; 3]> = match until {
        Until::Accuracy(_) => vec![[4; 3]; left.len()],
        Until::Res(r) => {
            let budget = 64.0 * left.len().max(1) as f64;
            let total = |step: f64| left.iter().map(|nd| per_axis(nd.b, step)).filter(|k| *k != [2.0; 3]).map(|k| k[0] * k[1] * k[2]).sum::<f64>();
            let mut step = r;
            while !(total(step) <= budget) {
                step *= 1.25;
            }
            left.iter().map(|nd| per_axis(nd.b, step).map(|k| k as usize)).collect()
        }
    };
    let finer: Vec<Option<Vec<f64>>> = left.par_iter().zip(&ks).map(|(nd, &k)| (k != [2; 3]).then(|| q.finer(&nd.state.tape, nd.b, k)).flatten()).collect();
    // How far that moved each estimate, twice over: what it is known to.
    let mut moved = vec![0.0; n];
    for (nd, f) in left.iter().zip(finer) {
        add(&mut total, &nd.state.share, 1.0);
        terms += 1.0;
        if let Some(e) = f {
            for k in 0..n {
                let change = e[k].clamp(nd.state.share[n + k], nd.state.share[2 * n + k]) - nd.state.share[k];
                total[k] += change;
                moved[k] += 2.0 * change.abs();
            }
        }
    }
    let slack = |k: usize| (terms + s.tested() as f64) * f64::EPSILON * (total[n + k].abs() + total[2 * n + k].abs());
    let lo: Vec<f64> = (0..n).map(|k| total[n + k] - slack(k)).collect();
    let hi: Vec<f64> = (0..n).map(|k| total[2 * n + k] + slack(k)).collect();
    let sums: Vec<f64> = (0..n).map(|k| total[k].clamp(lo[k], hi[k])).collect();
    let res = s.res().min(side);
    let boxes = s.tested() + 1;
    let error = (sums[0] - lo[0]).max(hi[0] - sums[0]);
    let why = match until {
        Until::Accuracy(a) if error > a * sums[0].abs() => {
            Some(format!("it may be {:.2}% from its estimate, {}", 100.0 * error / sums[0].abs(), why.as_deref().unwrap_or("its bounds no tighter")))
        }
        _ => why,
    };
    done(sums, lo, hi, res, boxes, why, moved)
}

/// Boxes tested before the first of `integrate`'s marks, and between a
/// mark and the next, as many as before it: a proof given up costs about
/// four times that (`out_of_reach` reads two rates).
const MARK: usize = 1 << 13;

/// Boxes a proof may be projected to need (`out_of_reach`): one that
/// would take more costs more than the estimate on cells it replaces,
/// but where walls are thin (goBILDA's 1200 mm goTUBE, projected to need
/// 192 thousand at 33 thousand, 3.5 s of CPU against the cells' 21 s).
const PROOF: usize = BOXES >> 2;

/// Whether a measurement's bounds are out of its reach, from its marks
/// (the first integral's error and estimate at each doubling of the boxes
/// tested, `tested` so far), accuracy `a`: at the rates the last two
/// doublings narrowed them, each would take more than PROOF boxes to come
/// within `a` (a field whose slopes do not narrow with its boxes: a jump,
/// an edge on every box). One rate alone understates a proof that has yet
/// to reach its walls (a tilted sheet: 853 thousand, then 213 thousand).
/// Then it is estimated on cells (`cells::estimate`) at the cost that
/// takes.
fn out_of_reach(marks: &[(f64, f64)], a: f64, tested: usize) -> bool {
    // The error falls as the boxes to the power `rate`.
    let needed = |(e0, _): (f64, f64), (e1, v1): (f64, f64), tested: f64| {
        let rate = (e0 / e1).ln() / 2f64.ln();
        if rate > 0.0 { tested * (e1 / (a * v1.abs())).powf(1.0 / rate) } else { f64::INFINITY }
    };
    let [.., m0, m1, m2] = *marks else { return false };
    needed(m0, m1, tested as f64 / 2.0) > PROOF as f64 && needed(m1, m2, tested as f64) > PROOF as f64
}

/// Bounds summed over `tested` boxes, widened by the rounding of a sum of
/// as many terms.
fn widened((lo, hi): (f64, f64), tested: usize) -> (f64, f64) {
    let slack = (tested as f64 + 1.0) * f64::EPSILON * (lo.abs() + hi.abs());
    (lo - slack, hi + slack)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink_core::graph::Var;

    /// A box's bounds hold the material a fine grid of points counts in
    /// it, and its estimate is near the count: a sphere, a tilted plane and
    /// an egg crate, over boxes all about their surfaces.
    #[test]
    fn a_boxs_bounds_hold_its_material() {
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let sq = |b: &mut crate::test_fields::B, v| b.un(UnOp::Square, v);
        let (xx, yy, zz) = (sq(&mut b, x), sq(&mut b, y), sq(&mut b, z));
        let r2 = b.bin(BinOp::Add, xx, yy);
        let r2 = b.bin(BinOp::Add, r2, zz);
        let r = b.un(UnOp::Sqrt, r2);
        let three = b.c(3.0);
        let sphere = b.bin(BinOp::Sub, r, three);
        let (k, h) = (b.c(0.3), b.c(0.7));
        let kx = b.bin(BinOp::Mul, k, x);
        let hy = b.bin(BinOp::Mul, h, y);
        let plane = b.bin(BinOp::Add, kx, hy);
        let plane = b.bin(BinOp::Sub, plane, z);
        let f = b.c(5.0);
        let (fx, fy) = (b.bin(BinOp::Mul, f, x), b.bin(BinOp::Mul, f, y));
        let (sx, sy) = (b.un(UnOp::Sin, fx), b.un(UnOp::Sin, fy));
        let crate_ = b.bin(BinOp::Mul, sx, sy);
        let ctx = Ctx::new(0.0, vec![]);
        let mut seed = 12345u64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for f in [sphere, plane, crate_] {
            let t = Arc::new(b.tape(f));
            let q = Integrals { shape: &t, gs: &[None], moments: false, ctx: &ctx, n: 1 };
            for _ in 0..300 {
                let c: V3 = std::array::from_fn(|_| 6.0 * rnd() - 3.0);
                let w: V3 = std::array::from_fn(|_| 0.05 + 1.5 * rnd());
                let bx: [I; 3] = std::array::from_fn(|k| I::new(c[k] - w[k] / 2.0, c[k] + w[k] / 2.0));
                let Some(got) = q.test(&t, bx) else { continue };
                let n: usize = 40;
                let mut v = Vec::new();
                let pts: Vec<V3> = (0..n * n * n).map(|i| std::array::from_fn(|k| bx[k].lo + w[k] * (((i / n.pow(k as u32)) % n) as f64 + 0.5) / n as f64)).collect();
                eval_points(&t, &ctx, &pts, &mut v);
                let vol = w[0] * w[1] * w[2];
                let counted = v.iter().filter(|&&x| x < 0.0).count() as f64 / pts.len() as f64 * vol;
                // A grid of n a side counts a face to within about 3 / n of the box.
                let slack = 3.0 / n as f64 * vol;
                let (est, lo, hi) = (got.share[0], got.share[1], got.share[2]);
                assert!(lo - slack <= counted && counted <= hi + slack, "{bx:?}: {lo} {hi}, counted {counted}");
                assert!(lo <= est && est <= hi && (est - counted).abs() <= (hi - lo).max(slack), "{bx:?}: {est} in {lo} {hi}, counted {counted}");
            }
        }
    }

    /// A measurement at a res its boxes stop short of estimates the boxes
    /// left at that res: a ball's, whose plane cuts in boxes far wider
    /// read its curved surface short of where it is, all one way.
    #[test]
    fn boxes_left_wider_than_the_res_are_estimated_at_it() {
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let (xx, yy, zz) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y), b.un(UnOp::Square, z));
        let r2 = b.bin(BinOp::Add, xx, yy);
        let r2 = b.bin(BinOp::Add, r2, zz);
        let r = b.un(UnOp::Sqrt, r2);
        let one = b.c(1.0);
        let ball = b.bin(BinOp::Sub, r, one);
        let t = Arc::new(b.tape(ball));
        let ctx = Ctx::new(0.0, vec![]);
        let exact = 4.0 / 3.0 * std::f64::consts::PI;
        let m = integrate_within(&t, &[None], false, &ctx, [-1.5; 3], [1.5; 3], Until::Res(0.01), 4096);
        assert!(!m.settled && m.lo[0] <= exact && exact <= m.hi[0], "{m:?}");
        // Its boxes' own 2 x 2 x 2 samples read it 9e-4 over; at the res, 1.5e-4.
        let off = (m.sums[0] - exact) / exact;
        assert!(off.abs() < 3e-4, "{off}: {m:?}");
    }

    /// A res far finer than the boxes left could reach asks for more
    /// samples than any count holds: they are coarsened to the budget, not
    /// counted past it (a ball of radius 1000 at 1e-5 mm).
    #[test]
    fn a_tiny_res_on_a_capped_search_is_sampled_within_budget() {
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let (xx, yy, zz) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y), b.un(UnOp::Square, z));
        let r2 = b.bin(BinOp::Add, xx, yy);
        let r2 = b.bin(BinOp::Add, r2, zz);
        let r = b.un(UnOp::Sqrt, r2);
        let radius = b.c(1000.0);
        let ball = b.bin(BinOp::Sub, r, radius);
        let t = Arc::new(b.tape(ball));
        let m = integrate_within(&t, &[None], false, &Ctx::new(0.0, vec![]), [-1500.0; 3], [1500.0; 3], Until::Res(1e-5), 4096);
        let exact = 4.0 / 3.0 * std::f64::consts::PI * 1e9;
        assert!(!m.settled && m.lo[0] <= exact && exact <= m.hi[0], "{m:?}");
    }

    /// A print's density, walls 1.6 mm thick at 1 and infill at 0.2 (a
    /// select on the shape grown 1.6 inward), is proven, not sampled: a
    /// box the infill's boundary crosses is bounded by its condition's
    /// slope, as a box the surface crosses is (a ball of radius 20).
    #[test]
    fn walls_and_infill_are_proven() {
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let (xx, yy, zz) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y), b.un(UnOp::Square, z));
        let r2 = b.bin(BinOp::Add, xx, yy);
        let r2 = b.bin(BinOp::Add, r2, zz);
        let r = b.un(UnOp::Sqrt, r2);
        let radius = b.c(20.0);
        let ball = b.bin(BinOp::Sub, r, radius);
        let wall = b.c(1.6);
        let inner = b.bin(BinOp::Add, ball, wall);
        let (sparse, solid) = (b.c(0.2), b.c(1.0));
        let rho = b.0.select(inner, sparse, solid);
        let (shape, rho) = (Arc::new(b.tape(ball)), b.tape(rho));
        let ctx = Ctx::new(0.0, vec![]);
        let m = integrate(&shape, &[Some(&rho), None], false, &ctx, [-20.0; 3], [20.0; 3], Until::Accuracy(ACCURACY));
        let ball = |r: f64| 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
        let exact = ball(20.0) - 0.8 * ball(18.4);
        assert!(m.settled && m.sampled.is_none() && m.lo[0] <= exact && exact <= m.hi[0] && m.error(0) <= ACCURACY * m.sums[0], "{m:?}");
        assert!(m.boxes < 60_000, "{} boxes", m.boxes);
    }

    /// A field scaled far below 1 (1e-100 (x + y + z - 4e-9) over a cube
    /// 4e-8 across) is bounded as one of ordinary size: the plane's terms,
    /// 1e-108 each, cubed would round to subnormals.
    #[test]
    fn a_tiny_field_is_bounded_as_any() {
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let sum = b.bin(BinOp::Add, x, y);
        let sum = b.bin(BinOp::Add, sum, z);
        let shift = b.c(4e-9);
        let plane = b.bin(BinOp::Sub, sum, shift);
        let scale = b.c(1e-100);
        let f = b.bin(BinOp::Mul, scale, plane);
        let shape = Arc::new(b.tape(f));
        let m = integrate(&shape, &[None], false, &Ctx::new(0.0, vec![]), [-2e-8; 3], [2e-8; 3], Until::Accuracy(ACCURACY));
        // The cube's volume below the plane: half of it, and the slab
        // between x + y + z = 0 and 4e-9, less its corners' cut-offs.
        let (side, d) = (4e-8f64, 4e-9f64);
        let (u, c) = (d + 1.5 * side, side);
        let simplex = |t: f64| if t > 0.0 { t.powi(3) / 6.0 } else { 0.0 };
        let below = simplex(u) - 3.0 * simplex(u - c) + 3.0 * simplex(u - 2.0 * c) - simplex(u - 3.0 * c);
        assert!((below / 3.67786667e-23 - 1.0).abs() < 1e-6, "{below}");
        assert!(m.lo[0] <= below && below <= m.hi[0], "{m:?}");
    }

    /// A box, off the cells' grid, weighs what it is to a millionth: a box
    /// its edges cross takes each face's cut, combined (one plane's cut, of
    /// a distance curved round the edge, read a 100 mm cube 0.1% heavy).
    #[test]
    fn a_box_weighs_what_it_is() {
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let (sx, sy, sz) = (b.slab(x, 0.3, 50.0), b.slab(y, 0.1, 50.0), b.slab(z, 0.7, 50.0));
        // std's box: the exact distance, curved round its edges outside.
        let xy = b.bin(BinOp::Emax, sx, sy);
        let cube = b.bin(BinOp::Emax, xy, sz);
        let shape = Arc::new(b.tape(cube));
        let m = integrate(&shape, &[None], false, &Ctx::new(0.0, vec![]), [-51.0; 3], [52.0; 3], Until::Accuracy(ACCURACY));
        assert!((m.sums[0] / 1e6 - 1.0).abs() < 1e-6, "{m:?}");
    }

    /// The volume of a box cut by a plane, from its closed form, against
    /// counting a fine grid of points: within ROUNDING of the box where
    /// the plane is tilted, exact where it is not.
    #[test]
    fn a_box_cut_by_a_plane_is_bounded_both_ways() {
        let h = [I::point(1.0), I::point(2.0), I::point(0.5)];
        for (r, a) in [(-0.3, [1.0, 0.5, -2.0]), (0.2, [-1.0, -0.25, 0.0]), (-1.0, [0.0, 0.0, 1.0]), (-0.7, [0.3, 0.2, 1e-6])] {
            let (up, down) = (cut(r, a, h, true), cut(r, a, h, false));
            let k = 200;
            let mut count = 0usize;
            for i in 0..k {
                for j in 0..k {
                    for l in 0..k {
                        let u = [(i as f64 + 0.5) / k as f64, 2.0 * (j as f64 + 0.5) / k as f64, 0.5 * (l as f64 + 0.5) / k as f64];
                        count += (r + a[0] * u[0] + a[1] * u[1] + a[2] * u[2] < 0.0) as usize;
                    }
                }
            }
            let grid = count as f64 / (k * k * k) as f64;
            assert!(down <= up && up - down <= 3e-3 && (down - 1e-2..=up + 1e-2).contains(&grid), "{r} {a:?}: {down} {up} counted {grid}");
        }
    }
}
