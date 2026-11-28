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
//! once halved until it holds one face of it, its share is known. Where
//! the field joins terms (a max or min, a magnitude: an edge, a thin
//! wall), each term is bounded so on its own slopes and the volumes are
//! joined (`branch_volume`), so a box need not part the faces either; where
//! it picks one (a select), the one picked there is. Boxes
//! wait widest bounds first and are halved on the axis the field changes
//! most along (a plate's across its thickness, not along it), until the
//! first integral's bounds are within the accuracy asked of it.

mod cells;
pub use cells::PASS_CELLS;

use crate::search::{self, Node, Precision, Smear, Test};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_interval_grad, eval_interval_grad_every, eval_point, eval_points};
use ink_core::graph::Fingerprint;
use ink_core::graph::{BinOp, UnOp};
use ink_core::interval::{G, I, add_o, div_o, mag, mul_o, sub_o, un};
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
/// within this share of its estimate either way (`Moments`: and, with
/// moments, each moment's within it of its own scale, for what reads a
/// centre or an inertia; without moments, either way, each other
/// integrand's within it of its bounds over the region, for what reads a
/// mean), or with every box halved to sides of at most this (mm) where its
/// bounds still differ.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Until {
    Accuracy(f64),
    Moments(f64),
    Res(f64),
}

/// Integrals as measured: each one's estimate (`sums`) and proven bounds
/// (`lo`, `hi`; the estimate lies within them), the smallest side of a box
/// tested (`res`), the boxes tested, whether it reached what it was asked
/// (`settled`; else why not, `why`), where its bounds were out of reach,
/// its estimate's last move (`sampled`: twice that move, relative to it;
/// `estimate`), which its spread is then, and how far each estimate is
/// known to (`moved`: twice what its last refinement moved it; 0 where
/// none was made). Moments are about `about`, the region's centre, so a
/// body keeps its own however far it is from the origin.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Measured {
    pub sums: Vec<f64>,
    pub about: V3,
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
    /// zx] about `about`, then (the second integrand's, where it is 1) the
    /// volume.
    pub fn moments(&self) -> crate::rigid::Moments {
        crate::rigid::Moments { about: self.about, m: std::array::from_fn(|k| self.sums[k]) }
    }

    pub fn volume(&self) -> f64 {
        self.sums[10]
    }

    /// Why it did not settle (None: it did).
    pub fn unsettled(&self) -> Option<String> {
        (!self.settled).then(|| self.why.clone().unwrap_or_default())
    }
}

/// `integrate`, cached as the kernel's other answers are. An integrand the
/// same at every point (a density) is a factor: the shape's own integrals
/// are measured, and scaled (`uniform`), so a body's volume and its moments
/// at the density that volume gives are one question. Its moments go on
/// from where its volume's search stopped (`Kept`).
pub fn measured(shape: &Arc<Tape>, gs: &[Option<&Tape>], moments: bool, ctx: &Ctx, lo: V3, hi: V3, until: Until) -> Measured {
    if let Some(m) = uniform(shape, gs, moments, ctx, lo, hi, until) {
        return m;
    }
    let tapes: Vec<&Tape> = std::iter::once(&**shape).chain(gs.iter().flatten().copied()).collect();
    let inputs = |fp: &mut Fingerprint, kind: u64, v: f64| {
        fp.word(kind).word(v.to_bits()).word(gs.len() as u64).word(moments as u64);
        gs.iter().for_each(|g| _ = fp.word(g.is_some() as u64));
    };
    let (kind, v) = match until {
        Until::Accuracy(a) => (0, a),
        Until::Res(r) => (1, r),
        Until::Moments(a) => (2, a),
    };
    // An accuracy's search and the moments' that goes on from it are one search of one question:
    // kept by the question, less where it stops (the kind of until). A paused search holds per-box
    // tapes indexing the caller's raw param slots, so the key also pins each tape's own param
    // numbers (which `cache::question` canonicalizes by first appearance).
    let kept = match until {
        Until::Accuracy(a) | Until::Moments(a) if moments && !gs.is_empty() => {
            let ctx = ctx.project(&tapes);
            let mut fp = Fingerprint::default();
            fp.str("measured, kept");
            lo.iter().chain(&hi).for_each(|x| _ = fp.f64(*x));
            inputs(&mut fp, 0, a);
            for t in &tapes {
                let ps = t.params_read();
                fp.word(ps.len() as u64);
                ps.iter().for_each(|&k| _ = fp.word(k as u64));
                for inst in &t.ssa {
                    if let Op::Var(ink_core::graph::Var::Param(k)) = inst.op {
                        fp.word(k as u64);
                    }
                }
                fp.word(u64::MAX);
            }
            crate::cache::question(&mut fp, &tapes, &ctx).map(|_| fp.value())
        }
        _ => None,
    };
    crate::cache::memo_region("measured", &tapes, ctx, (lo, hi), |fp| inputs(fp, kind, v), |ctx| integrate_kept(shape, gs, moments, ctx, lo, hi, until, kept))
}

/// `measured` where integrands are the same at every point of the region
/// (read no x, y or z): the shape's integrals with those integrands 1,
/// each of theirs times its value (estimates) and its range over the region
/// (bounds, rounded outward), as each box's would be. None: there is none,
/// or one's value or range is not finite.
fn uniform(shape: &Arc<Tape>, gs: &[Option<&Tape>], moments: bool, ctx: &Ctx, lo: V3, hi: V3, until: Until) -> Option<Measured> {
    let b: [I; 3] = std::array::from_fn(|k| I::new(lo[k], hi[k]));
    let c = about(lo, hi);
    let factors: Vec<Option<(f64, I)>> = gs
        .iter()
        .map(|g| {
            let g = (*g)?;
            let plain = g.ssa.iter().all(|x| match x.op {
                Op::Var(v) => !matches!(v, ink_core::graph::Var::X | ink_core::graph::Var::Y | ink_core::graph::Var::Z),
                Op::Const(_) | Op::Un(_) | Op::Bin(_) | Op::Select | Op::Tag(_) => true,
                _ => false,
            });
            let (v, r) = (eval_point(g, ctx, c), eval_interval(g, ctx, b, None));
            (plain && v.is_finite() && r.lo.is_finite() && r.hi.is_finite()).then_some((v, r))
        })
        .collect();
    if factors.iter().all(Option::is_none) {
        return None;
    }
    let ones: Vec<Option<&Tape>> = gs.iter().zip(&factors).map(|(g, f)| if f.is_some() { None } else { *g }).collect();
    let mut m = measured(shape, &ones, moments, ctx, lo, hi, until);
    let mut k = 0;
    for (i, f) in factors.iter().enumerate() {
        let terms = if i == 0 && moments { MOMENTS.len() } else { 1 };
        if let Some((v, r)) = *f {
            for j in k..k + terms {
                let w = mul_o(r, I::new(m.lo[j], m.hi[j]));
                (m.lo[j], m.hi[j]) = (w.lo, w.hi);
                m.sums[j] = (v * m.sums[j]).max(w.lo).min(w.hi);
                m.moved[j] *= v.abs();
            }
        }
        k += terms;
    }
    Some(m)
}

/// The point moments are taken about: the region's centre, so that a body
/// far from the origin keeps its own moments (about the origin, its
/// second moments are its distance squared, which rounds them away).
fn about(lo: V3, hi: V3) -> V3 {
    std::array::from_fn(|k| 0.5 * lo[k] + 0.5 * hi[k])
}

/// The terms of mass moments: [1, x, y, z, xx, yy, zz, xy, yz, zx] (of a
/// point less `about`).
pub(crate) const MOMENTS: [fn(V3) -> f64; 10] =
    [|_| 1.0, |p| p[0], |p| p[1], |p| p[2], |p| p[0] * p[0], |p| p[1] * p[1], |p| p[2] * p[2], |p| p[0] * p[1], |p| p[1] * p[2], |p| p[2] * p[0]];

/// The moments' terms' integrals over a box full of material, in closed
/// form (each a product of one per axis: the box's side, (h^2 - l^2) / 2,
/// (h^3 - l^3) / 3 of its ends l, h less `about`), rounded outward. A
/// box's moments bounded by the terms' ranges times its volume are a
/// box's own size wide, and a full box is never split: a cylinder's
/// inertia was proven to 47% where its mass was to 1.4%.
fn moment_integrals(b: [I; 3], about: V3) -> [I; 10] {
    let per: [[I; 3]; 3] = std::array::from_fn(|k| {
        // Each end less `about` as an interval holding it: the box's own ends rounded outward would hold
        // more than the box, and a full box's integral would be bounded too low as well as too high.
        let (l, h) = (sub_o(I::point(b[k].lo), I::point(about[k])), sub_o(I::point(b[k].hi), I::point(about[k])));
        let side = sub_o(I::point(b[k].hi), I::point(b[k].lo));
        let sq = |v| un(UnOp::Square, v);
        [side, mul_o(mul_o(side, add_o(h, l)), I::point(0.5)), div_o(mul_o(side, add_o(add_o(sq(h), mul_o(h, l)), sq(l))), I::point(3.0))]
    });
    const POW: [[usize; 3]; 10] = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [2, 0, 0], [0, 2, 0], [0, 0, 2], [1, 1, 0], [0, 1, 1], [1, 0, 1]];
    POW.map(|p| mul_o(mul_o(per[0][p[0]], per[1][p[1]]), per[2][p[2]]))
}

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
    /// What moments are about (`about`).
    about: V3,
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
            out.push(v.iter().zip(pts).zip(w).enumerate().map(|(p, ((a, x), b))| a * (term(std::array::from_fn(|k| x[k] - self.about[k])) + own(p)) * b).sum());
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
        let fill = |t: &Arc<Tape>| -> (Vec<f64>, Vec<V3>, Vec<V3>) {
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
        let terms = moment_ranges(std::array::from_fn(|k| sub_o(b[k], I::point(self.about[k]))));
        let whole = (full && self.moments).then(|| moment_integrals(b, self.about));
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
                let plain = mul_o(mul_o(range, terms[j]), vol);
                // A full box's term integrated exactly, times the integrand's range where that holds: it
                // is one value (a density), or the term keeps one sign over the box.
                let exact = whole.filter(|_| range.lo == range.hi || terms[j].lo >= 0.0 || terms[j].hi <= 0.0).map(|w| mul_o(range, w[j]));
                let r = match (j, exact) {
                    (0, _) => first,
                    (_, Some(e)) => I::new(e.lo.max(plain.lo), e.hi.min(plain.hi)),
                    (_, None) => plain,
                };
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
        let sub = tape.narrowed_much(&ch);
        let d = eval_interval_grad(&sub, self.ctx, b).d;
        let c = b.map(|x| x.mid());
        let fc = eval_interval(&sub, self.ctx, c.map(I::point), None);
        let (vlo, vhi) = slope_volume(fc, d, b);
        let vh = vhi.min(vol.hi);
        let v = I::new(vlo.max(0.0).min(vh), vh);
        let v = match branch_volume(&sub, self.ctx, b, vol) {
            Some(w) if v.lo.max(w.lo) <= v.hi.min(w.hi) => I::new(v.lo.max(w.lo), v.hi.min(w.hi)),
            _ => v,
        };
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
    cond: Arc<Tape>,
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
    Some(Jump { cond: Arc::new(g.cone(g.ssa[i].a as usize)), a, b: g.simplify(&pick) })
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

/// Bounds of the volume where `t` is below 0 in box `b` (volume `vol`),
/// from its terms (None: it has none, no max, min, magnitude or select). An edge
/// or a thin wall is where two faces meet: the field's slopes over a box
/// holding it span both faces' (a kink), so its own bounds
/// (`slope_volume`) leave out where the faces are, but each face's slopes
/// alone are tight. So through maxes (an emax is below 0 where its max
/// is), mins, negations and names, each term's volume is bounded on its
/// own slopes, and they are joined: the points where a max is below 0 are
/// those where both terms are (A and B), so at least |A| + |B| - |box| and
/// at most the lesser; a min's, where either is, at least the greater and
/// at most |A| + |B|; a negation's, where its term is above 0, bounded on
/// the negated term's slopes. A value a |g| + q, a a constant (through
/// sums, negations and constant scales: a tooth's flank, a slab about its
/// middle), is the max of a g + q and -a g + q where a is above 0 (the
/// min, below), each of whose slopes is its parts'. A select is the branch
/// its condition picks over the box (a tooth's side, on the sign of its
/// angle's sine), whose terms are followed in turn; one its condition
/// leaves undecided, below 0 at least where both branches are and at most
/// where either is. Each value's own bounds (below 0 throughout or
/// nowhere, or on its own slopes) hold too.
fn branch_volume(t: &Tape, ctx: &Ctx, b: [I; 3], vol: I) -> Option<I> {
    if !t.ssa.iter().any(|x| matches!(x.op, Op::Bin(BinOp::Max | BinOp::Min | BinOp::Emax) | Op::Un(UnOp::Abs) | Op::Select)) {
        return None;
    }
    let over = eval_interval_grad_every(t, ctx, b);
    let at = eval_interval_grad_every(t, ctx, b.map(|x| I::point(x.mid())));
    let terms = Terms { t, b, vol, over, at };
    let mut memo = vec![[None; 2]; t.ssa.len()];
    Some(terms.below(t.ssa.len() - 1, false, 64, &mut memo))
}

/// `branch_volume`'s tape and box, and every value's bounds over the box
/// (`over`) and at its centre (`at`).
struct Terms<'a> {
    t: &'a Tape,
    b: [I; 3],
    vol: I,
    over: Vec<G>,
    at: Vec<G>,
}

/// Values a kink's search (`Terms::kink`) may visit.
const KINK_SEARCH: usize = 64;

impl Terms<'_> {
    /// The volume where a value is below 0, from its value at the centre,
    /// its slopes and its range over the box: never NaN (an unbounded
    /// slope leaves it the whole box), which the joins' bounds rely on.
    fn own(&self, fc: I, d: [I; 3], range: I) -> I {
        if range.hi < 0.0 {
            return self.vol;
        }
        if range.lo >= 0.0 {
            return I::point(0.0);
        }
        let (lo, hi) = slope_volume(fc, d, self.b);
        let hi = hi.min(self.vol.hi);
        I::new(lo.max(0.0).min(hi), hi)
    }

    /// Where both of two sets of these volumes are.
    fn both(&self, a: I, c: I) -> I {
        I::new(sub_o(add_o(I::point(a.lo), I::point(c.lo)), I::point(self.vol.hi)).lo.max(0.0), a.hi.min(c.hi))
    }

    /// Where either is.
    fn either(&self, a: I, c: I) -> I {
        I::new(a.lo.max(c.lo), add_o(I::point(a.hi), I::point(c.hi)).hi.min(self.vol.hi))
    }

    /// The volume where value j is below 0 (above 0, where `flip`), its
    /// terms followed `depth` deep.
    fn below(&self, j: usize, flip: bool, depth: u32, memo: &mut Vec<[Option<I>; 2]>) -> I {
        if let Some(x) = memo[j][flip as usize] {
            return x;
        }
        let sign = |x: I| if flip { un(UnOp::Neg, x) } else { x };
        let (g, c) = (&self.over[j], &self.at[j]);
        let mine = self.own(sign(c.v), g.d.map(sign), sign(g.v));
        let ins = self.t.ssa[j];
        let (a, e) = (ins.a as usize, ins.b as usize);
        let joined = if depth == 0 || mine.lo == mine.hi {
            None
        } else {
            match ins.op {
                Op::Tag(_) | Op::Lipschitz => Some(self.below(a, flip, depth - 1, memo)),
                Op::Un(UnOp::Neg) => Some(self.below(a, !flip, depth - 1, memo)),
                // Below 0 where both are; above 0, where either is (a min the other way about).
                Op::Bin(BinOp::Max | BinOp::Emax | BinOp::Min) => {
                    let (x, y) = (self.below(a, flip, depth - 1, memo), self.below(e, flip, depth - 1, memo));
                    Some(if matches!(ins.op, Op::Bin(BinOp::Min)) == flip { self.both(x, y) } else { self.either(x, y) })
                }
                // The branch the condition's range picks (b below 0, c not), or between where both are and
                // where either is.
                Op::Select => {
                    let (cond, c) = (self.over[a].v, ins.c as usize);
                    Some(if cond.hi < 0.0 {
                        self.below(e, flip, depth - 1, memo)
                    } else if cond.lo >= 0.0 {
                        self.below(c, flip, depth - 1, memo)
                    } else {
                        let (x, y) = (self.below(e, flip, depth - 1, memo), self.below(c, flip, depth - 1, memo));
                        I::new(self.both(x, y).lo, self.either(x, y).hi)
                    })
                }
                _ => self.kink(j, flip),
            }
        };
        let x = match joined {
            Some(w) if mine.lo.max(w.lo) <= mine.hi.min(w.hi) => I::new(mine.lo.max(w.lo), mine.hi.min(w.hi)),
            _ => mine,
        };
        memo[j][flip as usize] = Some(x);
        x
    }

    /// Value j (negated where `flip`) as a |g| + q, g's range over the box
    /// on both sides of 0 (a kink) and a one nonzero number: the volume
    /// where it is below 0, from those of its sides a g + q and -a g + q.
    fn kink(&self, j: usize, flip: bool) -> Option<I> {
        let one = I::point(if flip { -1.0 } else { 1.0 });
        let mut rest = Vec::new();
        let mut budget = KINK_SEARCH;
        let (m, g, a) = self.path(j, one, &mut rest, &mut budget)?;
        // q must not read the magnitude: it would then not be q.
        if rest.iter().any(|&(k, _)| self.reads(k, m)) {
            return None;
        }
        // s a g + q (q alone where s is 0).
        let side = |s: f64| {
            let a = mul_o(a, I::point(s));
            let sum = |f: &dyn Fn(&G) -> I, of: &[G]| rest.iter().fold(if s == 0.0 { I::point(0.0) } else { mul_o(a, f(&of[g])) }, |acc, &(k, c)| add_o(acc, mul_o(c, f(&of[k]))));
            let d = std::array::from_fn(|k| sum(&|x| x.d[k], &self.over));
            self.own(sum(&|x| x.v, &self.at), d, sum(&|x| x.v, &self.over))
        };
        let (x, y, q) = (side(1.0), side(-1.0), side(0.0));
        let up = |u: f64, w: f64| add_o(I::point(u), I::point(w));
        // Where q is below 0 one side is (the lesser of a g and -a g is at most 0), so with a above 0 the
        // sides together hold at least q's, and where both are, they hold at most |A| + |B| less that (a
        // slab |g| - h: exactly its sides' less the box). With a below 0, where both are q is, so at most
        // q's, and either holds at least |A| + |B| less that.
        Some(if a.lo > 0.0 {
            let w = self.both(x, y);
            I::new(w.lo, w.hi.min(sub_o(up(x.hi, y.hi), I::point(q.lo)).hi.max(w.lo)))
        } else {
            let w = self.either(x, y);
            I::new(w.lo.max(sub_o(up(x.lo, y.lo), I::point(q.hi)).lo.min(w.hi)), w.hi)
        })
    }

    /// The magnitude value j is `coef` times plus a sum of other values
    /// (pushed to `rest`, each with its coefficient), through sums,
    /// negations, names and products and quotients by one nonzero number:
    /// the magnitude, its argument and coefficient. None: none whose
    /// argument's range is on both sides of 0, with a coefficient of one
    /// sign, within `budget` values.
    fn path(&self, j: usize, coef: I, rest: &mut Vec<(usize, I)>, budget: &mut usize) -> Option<(usize, usize, I)> {
        if *budget == 0 {
            return None;
        }
        *budget -= 1;
        let ins = self.t.ssa[j];
        let (a, b) = (ins.a as usize, ins.b as usize);
        // One nonzero number over the box.
        let number = |k: usize| {
            let v = self.over[k].v;
            (v.lo == v.hi && v.lo != 0.0 && v.lo.is_finite()).then_some(v)
        };
        let either = |rest: &mut Vec<(usize, I)>, budget: &mut usize, (p, cp): (usize, I), (o, co): (usize, I)| {
            let n = rest.len();
            rest.push((o, co));
            if let Some(x) = self.path(p, cp, rest, budget) {
                return Some(x);
            }
            rest.truncate(n);
            None
        };
        match ins.op {
            Op::Tag(_) | Op::Lipschitz => self.path(a, coef, rest, budget),
            Op::Un(UnOp::Neg) => self.path(a, un(UnOp::Neg, coef), rest, budget),
            Op::Un(UnOp::Abs) => {
                let r = self.over[a].v;
                (r.lo < 0.0 && r.hi > 0.0 && (coef.lo > 0.0 || coef.hi < 0.0)).then_some((j, a, coef))
            }
            Op::Bin(BinOp::Mul) if ins.a != ins.b => match (number(a), number(b)) {
                (_, Some(c)) => self.path(a, mul_o(coef, c), rest, budget),
                (Some(c), _) => self.path(b, mul_o(coef, c), rest, budget),
                _ => None,
            },
            Op::Bin(BinOp::Div) => number(b).and_then(|c| self.path(a, div_o(coef, c), rest, budget)),
            Op::Bin(op @ (BinOp::Add | BinOp::Sub)) => {
                let other = if op == BinOp::Sub { un(UnOp::Neg, coef) } else { coef };
                either(rest, budget, (a, coef), (b, other)).or_else(|| either(rest, budget, (b, other), (a, coef)))
            }
            _ => None,
        }
    }

    /// Whether value k reads value m (m among the values it is made of).
    fn reads(&self, k: usize, m: usize) -> bool {
        let mut seen = HashSet::new();
        let mut stack = vec![k];
        while let Some(i) = stack.pop() {
            if i == m {
                return true;
            }
            // A value reads only those before it.
            if i < m || !seen.insert(i) {
                continue;
            }
            let (args, n) = self.t.ssa[i].args();
            stack.extend(args[..n].iter().map(|&x| x as usize));
        }
        false
    }
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
fn shares(t: &Arc<Tape>, ctx: &Ctx, pts: &[V3], s: V3, depth: u32) -> Vec<(f64, V3, V3)> {
    shares_with(t, pts, s, depth, |roots| ink_core::eval::eval_grads_at(t, ctx, pts, roots))
}

/// Which branches contribute to the estimate, evaluated on the same tape:
/// the gradient at a branch is the one its separate cone would give. A
/// child shared by two branches needs only one reading.
fn shares_with(t: &Tape, pts: &[V3], s: V3, depth: u32, grads: impl FnOnce(&[usize]) -> Vec<Vec<[f64; 4]>>) -> Vec<(f64, V3, V3)> {
    struct Branch {
        gradient: usize,
        children: Option<(usize, usize, bool)>,
    }
    fn visit(t: &Tape, mut i: usize, depth: u32, branches: &mut Vec<Branch>) -> usize {
        let gradient = i;
        while let Op::Tag(_) | Op::Lipschitz = t.ssa[i].op {
            i = t.ssa[i].a as usize;
        }
        let ins = t.ssa[i];
        let children = match ins.op {
            Op::Bin(op @ (BinOp::Min | BinOp::Max | BinOp::Emax)) if depth > 0 => {
                Some((visit(t, ins.a as usize, depth - 1, branches), visit(t, ins.b as usize, depth - 1, branches), op != BinOp::Min))
            }
            _ => None,
        };
        branches.push(Branch { gradient, children });
        branches.len() - 1
    }
    let mut branches = Vec::new();
    let root = visit(t, t.len() - 1, depth, &mut branches);
    // A combined root's own gradient is never read.
    let read = if branches[root].children.is_some() { root } else { root + 1 };
    let mut roots: Vec<_> = branches[..read].iter().map(|b| b.gradient).collect();
    roots.sort_unstable();
    roots.dedup();
    let gradients = grads(&roots);
    for b in &mut branches[..read] {
        b.gradient = roots.binary_search(&b.gradient).unwrap();
    }
    let mut values = Vec::with_capacity(branches.len());
    (0..pts.len())
        .map(|p| {
            values.clear();
            for b in &branches {
                values.push(match b.children {
                    Some((a, b, intersection)) => combined(values[a], values[b], &gradients[branches[a].gradient][p], &gradients[branches[b].gradient][p], intersection),
                    None => {
                        let g = &gradients[b.gradient][p];
                        let f = fraction(g, s);
                        let (c, v) = centroid(g, s, f);
                        (f, c, v)
                    }
                });
            }
            values[root]
        })
        .collect()
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
    integrate_within(shape, gs, moments, ctx, lo, hi, until, BOXES, None, false).0
}

/// `integrate`, its search kept under `key` (`Kept`): where an accuracy's
/// stops, for the moments' to go on from.
#[allow(clippy::too_many_arguments)]
fn integrate_kept(shape: &Arc<Tape>, gs: &[Option<&Tape>], moments: bool, ctx: &Ctx, lo: V3, hi: V3, until: Until, key: Option<u128>) -> Measured {
    let resumed = match (until, key) {
        (Until::Moments(_), Some(k)) => match kept(k) {
            Some(Kept::Done(m)) => return m,
            Some(Kept::Paused(p)) => Some(p),
            None => None,
        },
        _ => None,
    };
    let keeps = key.is_some() && matches!(until, Until::Accuracy(_));
    let (m, k) = integrate_within(shape, gs, moments, ctx, lo, hi, until, BOXES, resumed, keeps);
    if let (Some(key), Some(k)) = (key, k) {
        keep(key, k);
    }
    m
}

/// A measure with moments to an accuracy (`Until::Accuracy`), as its search
/// left it: the moments' measure of the same question to the same accuracy
/// (`Until::Moments`) is that search gone on until its moments are within
/// it too, as the two test the same boxes in the same order until the first
/// stops, so it goes on from there (`Paused`), or, where the first stopped
/// for any other reason than its accuracy, it is the same answer (`Done`).
/// A body's volume is measured for the density its published mass gives,
/// then its moments at that density: one search, not two.
enum Kept {
    Paused(Paused),
    Done(Measured),
}

/// A search where its first integral reached its accuracy: its boxes left,
/// the open and settled boxes' shares, and its marks (`integrate_within`).
struct Paused {
    s: search::Search<Held, ()>,
    open: Vec<f64>,
    settled: Vec<f64>,
    marks: Vec<(usize, f64, f64)>,
    moment_marks: Vec<(usize, f64, f64)>,
    mark: usize,
}

/// Searches kept, the latest last, by question: a body's moments are asked
/// just after its volume, so a few are enough, and each holds its boxes.
static KEPT: std::sync::Mutex<Vec<(u128, Kept)>> = std::sync::Mutex::new(Vec::new());
const KEEP: usize = 2;

fn keep(key: u128, k: Kept) {
    let mut dropped = Vec::new();
    {
        let mut kept = KEPT.lock().unwrap_or_else(|e| e.into_inner());
        let mut i = 0;
        while i < kept.len() {
            if kept[i].0 == key {
                dropped.push(kept.remove(i));
            } else {
                i += 1;
            }
        }
        if kept.len() >= KEEP {
            dropped.push(kept.remove(0));
        }
        kept.push((key, k));
    }
    drop(dropped);
}

fn kept(key: u128) -> Option<Kept> {
    let mut kept = KEPT.lock().unwrap_or_else(|e| e.into_inner());
    let i = kept.iter().position(|(q, _)| *q == key)?;
    Some(kept.remove(i).1)
}

/// `integrate`, testing at most `boxes` boxes, going on from `resumed`
/// (an accuracy's search of the same question, for its moments), and,
/// where `keeps` (an accuracy's, with moments), with what its moments'
/// measure is to go on from (`Kept`; None: it was stopped by the call's
/// time limit, which another call need not reach).
#[allow(clippy::too_many_arguments)]
fn integrate_within(
    shape: &Arc<Tape>,
    gs: &[Option<&Tape>],
    moments: bool,
    ctx: &Ctx,
    lo: V3,
    hi: V3,
    until: Until,
    boxes: usize,
    resumed: Option<Paused>,
    keeps: bool,
) -> (Measured, Option<Kept>) {
    let keeps = keeps && moments && !gs.is_empty() && matches!(until, Until::Accuracy(_));
    let n = gs.len() + if moments && !gs.is_empty() { 9 } else { 0 };
    let q = Integrals { shape, gs, moments, about: about(lo, hi), ctx, n };
    let b: [I; 3] = std::array::from_fn(|k| I::new(lo[k], hi[k]));
    let done = |sums: Vec<f64>, lo: Vec<f64>, hi: Vec<f64>, res: f64, boxes: usize, why: Option<String>, moved: Vec<f64>| Measured {
        about: q.about,
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
    let Some(root) = q.test(q.shape, b).filter(|_| n > 0) else { return (done(vec![0.0; n], vec![0.0; n], vec![0.0; n], side, 1, None, vec![0.0; n]), None) };
    // Each integral's scale: its bounds over the whole box. The integrals
    // of the integrands themselves (not the moments' x, y, z terms, near 0
    // about the origin) order the boxes.
    let firsts: Vec<usize> = (0..gs.len()).map(|i| if moments && i > 0 { i + 9 } else { i }).collect();
    let scale: Vec<f64> = (0..n).map(|k| root.share[n + k].abs().max(root.share[2 * n + k].abs())).collect();
    // With moments, the boxes holding most of the second moments' doubt go first too, whatever is asked:
    // ordered by the mass's doubt alone, six catalog parts' mass proofs gave up at 13 to 27% (a goBILDA
    // bracket 27.6 to 47.5 g), and ordered so they reach 2% (40.0 to 41.1 g), no slower.
    let ordered: Vec<usize> = firsts.iter().copied().chain(if moments && !gs.is_empty() { 4..7 } else { 0..0 }).collect();
    let width = |s: &[f64]| ordered.iter().filter(|&&k| scale[k] > 0.0).map(|&k| (s[2 * n + k] - s[n + k]) / scale[k]).sum::<f64>();
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
    let finest = match until {
        Until::Res(r) => r,
        Until::Accuracy(_) | Until::Moments(_) => 0.0,
    };
    // At each doubling of the boxes tested: the first integral's error and estimate, and the moments'
    // largest share of the accuracy (as an error of that share of the accuracy, on 1).
    let Paused { mut s, mut open, mut settled, mut marks, mut moment_marks, mut mark } = resumed.unwrap_or_else(|| {
        let (mut settled, mut open) = (vec![0.0; 3 * n], vec![0.0; 3 * n]);
        let first = place(root, &mut open, &mut settled);
        let s = search::Search::start(
            vec![(Node { b, ctx: ctx_arc, state: Held { tape: shape.clone(), share: Arc::from(vec![]) } }, first)],
            Precision { finest, splits: 0, maps: true, closing: false },
        );
        Paused { s, open, settled, marks: Vec::new(), moment_marks: Vec::new(), mark: MARK }
    });
    // The first integral's error against the accuracy (reached at 1 or under), and, with moments, the
    // largest of each moment's against its own scale: a second moment's, itself; a first moment's, the
    // mass's and that axis's second moment's geometric mean (what moves the centre by that share of its
    // radius of gyration); a product's, its two second moments'. Without moments, the largest of each
    // other integrand's against its scale (its bounds over the whole region, which order the boxes): a
    // mean's moments, which its claim decides on with its weight. A weight proven in two boxes (a cube's
    // faces are exact in any box) left them where those two boxes put them: a centre of mass 0 to 10.
    let reached = |open: &[f64], settled: &[f64], a: f64| -> (f64, f64) {
        let at = |k: usize| (open[k] + settled[k], open[n + k] + settled[n + k], open[2 * n + k] + settled[2 * n + k]);
        let share = |k: usize, scale: f64| {
            let (est, l, h) = at(k);
            (est - l).max(h - est) / (a * scale)
        };
        let size = |k: usize| at(k).0.abs();
        let moment = match moments && !gs.is_empty() {
            false => firsts[1.min(firsts.len())..].iter().filter(|&&k| scale[k] > 0.0).map(|&k| share(k, scale[k])).fold(0.0, f64::max),
            true => (4..7)
                .map(|k| share(k, size(k)))
                .chain((1..4).map(|k| share(k, (size(0) * size(k + 3)).sqrt())))
                .chain([(7, 4, 5), (8, 5, 6), (9, 6, 4)].map(|(k, i, j)| share(k, (size(i) * size(j)).sqrt())))
                .fold(0.0, f64::max),
        };
        (share(0, size(0)), moment)
    };
    let mut why = None;
    // Whether it stopped where its first integral reached its accuracy.
    let mut reached_it = false;
    loop {
        // An accuracy's measure with moments does not wait for them, but marks them as the moments'
        // would, which is then that measure gone on.
        let (first, others) = match until {
            Until::Accuracy(a) | Until::Moments(a) => reached(&open, &settled, a),
            Until::Res(_) => (f64::INFINITY, 0.0),
        };
        let moment = if matches!(until, Until::Accuracy(_)) && moments { 0.0 } else { others };
        if first <= 1.0 && moment <= 1.0 {
            reached_it = true;
            break;
        }
        if let Until::Accuracy(a) | Until::Moments(a) = until
            && s.tested() >= mark
        {
            mark *= 2;
            let (est, l, h) = (open[0] + settled[0], open[n] + settled[n], open[2 * n] + settled[2 * n]);
            marks.push((s.tested(), (est - l).max(h - est), est));
            moment_marks.push((s.tested(), others * a, 1.0));
            if first > 1.0 && out_of_reach(&marks, a) {
                // The bounds so far, from the settled boxes' and those
                // left's shares, the halved ones' not taken away.
                let mut total = settled.clone();
                s.left().for_each(|nd| add(&mut total, &nd.state.share, 1.0));
                let (l, h): (Vec<f64>, Vec<f64>) = (0..n).map(|k| widened((total[n + k], total[2 * n + k]), s.tested())).unzip();
                drop(s);
                let m = cells::estimate(shape, gs, moments, ctx, lo, hi, a, (l, h));
                let k = (keeps && !crate::stop::reached()).then(|| Kept::Done(m.clone()));
                return (m, k);
            }
            // The mass reached and the moments out of reach as the mass would be: their bounds as far as
            // they came, each proven.
            if first <= 1.0 && out_of_reach(&moment_marks, a) {
                break;
            }
        }
        if s.tested() >= boxes {
            // (The moments short of their accuracy leave the measure settled: their bounds hold.)
            if first > 1.0 {
                why = Some(format!("after {boxes} boxes"));
            }
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
    let mut total = settled.clone();
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
        Until::Accuracy(_) | Until::Moments(_) => vec![[4; 3]; left.len()],
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
        Until::Accuracy(a) | Until::Moments(a) if error > a * sums[0].abs() => {
            Some(format!("it may be {:.2}% from its estimate, {}", 100.0 * error / sums[0].abs(), why.as_deref().unwrap_or("its bounds no tighter")))
        }
        _ => why,
    };
    let m = done(sums, lo, hi, res, boxes, why, moved);
    let k = (keeps && !crate::stop::reached()).then(|| if reached_it { Kept::Paused(Paused { s, open, settled, marks, moment_marks, mark }) } else { Kept::Done(m.clone()) });
    (m, k)
}

/// Boxes tested before the first of `integrate`'s marks, and between a
/// mark and the next, as many as before it: a proof given up costs about
/// four times that (`out_of_reach` reads two rates), or eight (three).
const MARK: usize = 1 << 13;

/// Boxes a proof may be projected to need (`out_of_reach`): one that
/// would take more costs more than the estimate on cells it replaces,
/// but where walls are thin (goBILDA's 1200 mm goTUBE, projected to need
/// 192 thousand at 33 thousand, 3.5 s of CPU against the cells' 21 s).
const PROOF: usize = BOXES >> 2;

/// Whether a measurement's bounds are out of its reach, from its marks
/// (the boxes tested, the first integral's error and its estimate, at
/// each doubling of the boxes tested or past it), accuracy `a`: at the
/// rates the last two of them narrowed the bounds (each against the boxes
/// it took), each would take more than PROOF boxes to come within `a` (a
/// field whose slopes do not narrow with its boxes: a jump, an edge on
/// every box). One rate alone understates a proof that has yet to reach
/// its walls (a tilted sheet: 853 thousand, then 213 thousand), and two
/// do where the later is faster, its boxes still coming down to its
/// smallest features (a toothed pulley's 0.15 mm fillets: 0.49, 0.53,
/// then 0.84, given up at 35 thousand boxes and proven at 150 thousand):
/// then the rate before them must project it past PROOF too. Then it is
/// estimated on cells (`cells::estimate`) at the cost that takes.
fn out_of_reach(marks: &[(usize, f64, f64)], a: f64) -> bool {
    // The error falls as the boxes to the power `rate`.
    let rate = |(t0, e0, _): (usize, f64, f64), (t1, e1, _): (usize, f64, f64)| (e0 / e1).ln() / (t1 as f64 / t0 as f64).ln();
    let needed = |m0, m1: (usize, f64, f64)| {
        let r = rate(m0, m1);
        if r > 0.0 { m1.0 as f64 * (m1.1 / (a * m1.2.abs())).powf(1.0 / r) } else { f64::INFINITY }
    };
    let [.., m0, m1, m2] = *marks else { return false };
    let over = |m0, m1| needed(m0, m1) > PROOF as f64;
    over(m0, m1) && over(m1, m2) && (!(rate(m1, m2) > rate(m0, m1)) || matches!(*marks, [.., m, _, _, _] if over(m, m0)))
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

    #[test]
    fn a_shared_edge_reads_its_gradients_in_one_pass() {
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let xy = b.bin(BinOp::Max, x, y);
        let xz = b.bin(BinOp::Max, x, z);
        let field = b.bin(BinOp::Min, xy, xz);
        let tape = b.tape(field);
        let pts: Vec<_> = (0..73).map(|i| [i as f64 / 32.0 - 1.0, 0.25, -0.25]).collect();
        let ctx = Ctx::new(0.0, vec![]);
        let sides = [1.0; 3];
        let got = shares_with(&tape, &pts, sides, 3, |roots| {
            // x, y, z and the two intersections. The common x occurs twice.
            assert_eq!(roots.len(), 5);
            let g = ink_core::eval::eval_grads_at(&tape, &ctx, &pts, roots);
            for (&i, row) in roots.iter().zip(&g) {
                assert_eq!(*row, ink_core::eval::eval_grads(&tape.cone(i), &ctx, &pts));
            }
            g
        });
        let gradients = |f| ink_core::eval::eval_grads(&b.tape(f), &ctx, &pts);
        let [x, y, z, xy, xz] = [x, y, z, xy, xz].map(gradients);
        let leaf = |g| {
            let f = fraction(g, sides);
            let (c, v) = centroid(g, sides, f);
            (f, c, v)
        };
        for (i, got) in got.iter().enumerate() {
            let a = combined(leaf(&x[i]), leaf(&y[i]), &x[i], &y[i], true);
            let b = combined(leaf(&x[i]), leaf(&z[i]), &x[i], &z[i], true);
            assert_eq!(*got, combined(a, b, &xy[i], &xz[i], false));
        }
    }

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
            let q = Integrals { shape: &t, gs: &[None], moments: false, about: [0.0; 3], ctx: &ctx, n: 1 };
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

    /// A box's bounds from its field's terms (`branch_volume`) hold its
    /// material: random maxes, mins, emaxes, negations, selects and magnitudes
    /// (a |g| + q, a slab about its middle, a wall about its mid-surface)
    /// of planes, balls, cylinders, an egg crate and std's thread, in random
    /// boxes, against a count on a grid. And a thin wall whose two faces
    /// share a box is bounded to the grid's count, where its own slopes
    /// (spanning both faces) left out where it is.
    #[test]
    fn a_boxs_bounds_from_its_terms_hold_its_material() {
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let mut seed = 987654321u64;
        let mut rnd = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let ctx = Ctx::new(0.0, vec![]);
        fn leaf(b: &mut crate::test_fields::B, [x, y, z]: [ink_core::graph::NodeId; 3], kind: usize, r: &mut dyn FnMut() -> f64) -> ink_core::graph::NodeId {
            let mut c = |b: &mut crate::test_fields::B, lo: f64, hi: f64| b.c(lo + (hi - lo) * r());
            match kind {
                0 => {
                    let (a, e, f, d) = (c(b, -1.0, 1.0), c(b, -1.0, 1.0), c(b, -1.0, 1.0), c(b, -1.0, 1.0));
                    let (ax, ey, fz) = (b.bin(BinOp::Mul, a, x), b.bin(BinOp::Mul, e, y), b.bin(BinOp::Mul, f, z));
                    let s = b.bin(BinOp::Add, ax, ey);
                    let s = b.bin(BinOp::Add, s, fz);
                    b.bin(BinOp::Sub, s, d)
                }
                1 | 2 => {
                    let (cx, cy, cz, rad) = (c(b, -2.0, 2.0), c(b, -2.0, 2.0), c(b, -2.0, 2.0), c(b, 0.5, 2.5));
                    let (dx, dy, dz) = (b.bin(BinOp::Sub, x, cx), b.bin(BinOp::Sub, y, cy), b.bin(BinOp::Sub, z, cz));
                    let (dx, dy, dz) = (b.un(UnOp::Square, dx), b.un(UnOp::Square, dy), b.un(UnOp::Square, dz));
                    let s = b.bin(BinOp::Add, dx, dy);
                    // A ball, or a cylinder about z.
                    let s = if kind == 1 { b.bin(BinOp::Add, s, dz) } else { s };
                    let s = b.un(UnOp::Sqrt, s);
                    b.bin(BinOp::Sub, s, rad)
                }
                3 => {
                    let (f, k) = (b.c(3.0), c(b, -0.5, 0.5));
                    let (fx, fy) = (b.bin(BinOp::Mul, f, x), b.bin(BinOp::Mul, f, y));
                    let (sx, sy) = (b.un(UnOp::Sin, fx), b.un(UnOp::Sin, fy));
                    let s = b.bin(BinOp::Mul, sx, sy);
                    b.bin(BinOp::Sub, s, k)
                }
                _ => b.thread_iso([x, y, z], 6.0, 1.0),
            }
        }
        fn tree(b: &mut crate::test_fields::B, at: [ink_core::graph::NodeId; 3], depth: u32, r: &mut dyn FnMut() -> f64) -> ink_core::graph::NodeId {
            let pick = (r() * 10.0) as usize;
            if depth == 0 || pick >= 8 {
                let kind = (r() * 5.0) as usize;
                return leaf(b, at, kind, r);
            }
            let u = tree(b, at, depth - 1, r);
            match pick {
                0 | 1 => {
                    let v = tree(b, at, depth - 1, r);
                    b.bin(if pick == 0 { BinOp::Max } else { BinOp::Min }, u, v)
                }
                2 => {
                    let v = tree(b, at, depth - 1, r);
                    b.bin(BinOp::Emax, u, v)
                }
                3 => b.un(UnOp::Neg, u),
                // A select between two terms on a third's sign.
                7 => {
                    let (v, cond) = (tree(b, at, depth - 1, r), tree(b, at, depth - 1, r));
                    b.0.select(cond, u, v)
                }
                // a |u| + q: a wall about u's zero, q a term or a constant.
                _ => {
                    let a = b.c(if pick == 6 { -(0.2 + r()) } else { 0.2 + r() });
                    let m = b.un(UnOp::Abs, u);
                    let m = b.bin(BinOp::Mul, a, m);
                    let q = if pick == 5 { tree(b, at, depth - 1, r) } else { b.c(-0.3 * r()) };
                    let q2 = b.c(0.5);
                    let q = b.bin(BinOp::Mul, q2, q);
                    if r() < 0.5 { b.bin(BinOp::Add, m, q) } else { b.bin(BinOp::Sub, q, m) }
                }
            }
        }
        let (mut tested, mut tighter) = (0, 0);
        for _ in 0..60 {
            let f = tree(&mut b, [x, y, z], 3, &mut rnd);
            let t = Arc::new(b.tape(f));
            let q = Integrals { shape: &t, gs: &[None], moments: false, about: [0.0; 3], ctx: &ctx, n: 1 };
            for _ in 0..40 {
                let c: V3 = std::array::from_fn(|_| 5.0 * rnd() - 2.5);
                let w: V3 = std::array::from_fn(|_| 0.05 + 1.2 * rnd());
                let bx: [I; 3] = std::array::from_fn(|k| I::new(c[k] - w[k] / 2.0, c[k] + w[k] / 2.0));
                let Some(got) = q.test(&t, bx) else { continue };
                let n: usize = 40;
                let mut v = Vec::new();
                let pts: Vec<V3> = (0..n * n * n).map(|i| std::array::from_fn(|k| bx[k].lo + w[k] * (((i / n.pow(k as u32)) % n) as f64 + 0.5) / n as f64)).collect();
                eval_points(&t, &ctx, &pts, &mut v);
                let vol = w[0] * w[1] * w[2];
                let counted = v.iter().filter(|&&x| x < 0.0).count() as f64 / pts.len() as f64 * vol;
                let slack = 3.0 / n as f64 * vol;
                let (lo, hi) = (got.share[1], got.share[2]);
                assert!(lo - slack <= counted && counted <= hi + slack, "{bx:?}: {lo} {hi}, counted {counted}");
                tested += 1;
                // Tighter than the field's own slopes leave it.
                let side = bx.map(|x| sub_o(I::point(x.hi), I::point(x.lo)));
                let whole = mul_o(mul_o(side[0], side[1]), side[2]);
                let mut ch = vec![BOTH; t.len()];
                eval_interval(&t, &ctx, bx, Some(&mut ch));
                let sub = t.narrowed_much(&ch);
                let (l, h) = slope_volume(eval_interval(&sub, &ctx, bx.map(|x| I::point(x.mid())), None), eval_interval_grad(&sub, &ctx, bx).d, bx);
                tighter += (hi - lo < 0.9 * (h.min(whole.hi) - l.max(0.0))) as usize;
            }
        }
        assert!(tested > 500 && tighter > tested / 20, "{tested} boxes, {tighter} tighter");
        // A wall 0.2 thick about x = 0.31 (as |x - 0.31| - 0.1), in a box 1 wide: its slopes along x span
        // -1..1, its own bounds 0..1, its sides' the wall (to the cut's rounding).
        let wall = b.slab(x, 0.31, 0.1);
        let t = Arc::new(b.tape(wall));
        let q = Integrals { shape: &t, gs: &[None], moments: false, about: [0.0; 3], ctx: &ctx, n: 1 };
        let got = q.test(&t, [I::new(0.0, 1.0); 3]).unwrap();
        assert!((got.share[1] - 0.2).abs() < 1e-6 && (got.share[2] - 0.2).abs() < 1e-6, "{:?}", got.share);
        // The same wall as the branch a select picks throughout the box (y below 2), the other a ball
        // crossing it: bounded as the wall alone, where the select's own slopes left it 0..1.
        let (cond, ball) = (b.c(2.0), leaf(&mut b, [x, y, z], 1, &mut || 0.5));
        let cond = b.bin(BinOp::Sub, y, cond);
        let picked = b.0.select(cond, wall, ball);
        let t = Arc::new(b.tape(picked));
        let q = Integrals { shape: &t, gs: &[None], moments: false, about: [0.0; 3], ctx: &ctx, n: 1 };
        let got = q.test(&t, [I::new(0.0, 1.0); 3]).unwrap();
        assert!((got.share[1] - 0.2).abs() < 1e-6 && (got.share[2] - 0.2).abs() < 1e-6, "{:?}", got.share);
    }

    /// A proof is given up where its marks' last two rates each project it
    /// past PROOF boxes and the later is no faster, or, where it is faster
    /// (its boxes still coming down to its smallest features), the rate
    /// before those projects it past PROOF too. Rates are each against the
    /// boxes its marks were taken at. A toothed rim on a bore at 2% (its
    /// error, relative, at 8, 17 and 33 thousand boxes: 9.2, 7.2 and 5.3%,
    /// then 3.2% at 67 thousand, proven at 134 thousand) was given up at 33
    /// thousand; a rate that stays about 0.5 (a field that jumps across a
    /// face) is given up at the third mark, and one that rises once at the
    /// fourth.
    #[test]
    fn a_proof_is_given_up_only_where_its_rates_stop_rising() {
        let rim = [(8221, 0.092156, 1.0), (16613, 0.071998, 1.0), (33235, 0.052626, 1.0), (67069, 0.03205, 1.0)];
        assert!(!out_of_reach(&rim[..3], ACCURACY) && !out_of_reach(&rim, ACCURACY));
        let at = |rates: &[f64]| {
            let mut m = vec![(1usize << 13, 0.2, 1.0)];
            for r in rates {
                let &(t, e, _) = m.last().unwrap();
                m.push((2 * t, e * 2f64.powf(-r), 1.0));
            }
            m
        };
        assert!(out_of_reach(&at(&[0.5, 0.48]), ACCURACY));
        assert!(!out_of_reach(&at(&[0.5, 0.55]), ACCURACY) && out_of_reach(&at(&[0.5, 0.48, 0.55]), ACCURACY));
        // A round that took the boxes to four times: its rate, 0.48 against them, read as one doubling's
        // would be 0.95, which projects 257 thousand, under PROOF.
        let past = [(8192, 0.2, 1.0), (16384, 0.2 * 2f64.powf(-0.5), 1.0), (65536, 0.2 * 2f64.powf(-1.45), 1.0)];
        assert!(out_of_reach(&past, ACCURACY));
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
        let m = integrate_within(&t, &[None], false, &ctx, [-1.5; 3], [1.5; 3], Until::Res(0.01), 4096, None, false).0;
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
        let m = integrate_within(&t, &[None], false, &Ctx::new(0.0, vec![]), [-1500.0; 3], [1500.0; 3], Until::Res(1e-5), 4096, None, false).0;
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

    /// A cube 2 mm across centred at (10, 20, 30), and one at 1e8: its
    /// moments are measured about its box's centre, within their bounds,
    /// so far out it keeps its own (about the origin its distance squared
    /// would round them away).
    #[test]
    fn moments_are_measured_about_the_body() {
        for c in [[10.0, 20.0, 30.0], [1e8; 3]] {
            let mut b = crate::test_fields::B::new();
            let axes = [Var::X, Var::Y, Var::Z].map(|v| b.v(v));
            let faces: Vec<_> = (0..3)
                .map(|k| {
                    let at = b.c(c[k]);
                    let d = b.bin(BinOp::Sub, axes[k], at);
                    let a = b.un(UnOp::Abs, d);
                    let one = b.c(1.0);
                    b.bin(BinOp::Sub, a, one)
                })
                .collect();
            let m01 = b.bin(BinOp::Max, faces[0], faces[1]);
            let cube = b.bin(BinOp::Max, m01, faces[2]);
            let shape = Arc::new(b.tape(cube));
            let (lo, hi) = (c.map(|v| v - 1.5), c.map(|v| v + 1.5));
            let m = integrate(&shape, &[None], true, &Ctx::new(0.0, vec![]), lo, hi, Until::Accuracy(ACCURACY));
            // About its centre: mass 8, no first moment, 8/3 squared (to
            // what its samples' coordinates keep, 1.5e-8 at 1e8), each within
            // its bounds.
            let want = [8.0, 0.0, 0.0, 0.0, 8.0 / 3.0, 8.0 / 3.0, 8.0 / 3.0, 0.0, 0.0, 0.0];
            assert_eq!(m.about, c);
            for k in 0..10 {
                assert!((m.sums[k] - want[k]).abs() <= 1e-7 && m.lo[k] <= want[k] + 1e-7 && want[k] - 1e-7 <= m.hi[k], "{k}: {m:?}");
            }
        }
    }

    /// A cylinder r 7.3, h 41.9 and a ball r 10, asked for their moments
    /// (`Until::Moments`): their mass, first and
    /// second moments proven to the accuracy asked, each against its own
    /// scale, and holding the exact ones. A full box's moments were bounded
    /// by their range times its volume, and the search stopped at the mass:
    /// the cylinder's second moments were proven to 47%.
    #[test]
    fn a_cylinders_and_a_balls_moments_are_proven_to_the_accuracy() {
        use std::f64::consts::PI;
        let (r, h, z0) = (7.3, 41.9, 5.1);
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let (xx, yy) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y));
        let rr = b.bin(BinOp::Add, xx, yy);
        let rho = b.un(UnOp::Sqrt, rr);
        let rc = b.c(r);
        let side = b.bin(BinOp::Sub, rho, rc);
        let ends = b.slab(z, z0, h / 2.0);
        let cyl = b.bin(BinOp::Max, side, ends);
        let cyl = Arc::new(b.tape(cyl));
        let m = PI * r * r * h;
        let cylinder = (cyl, [-8.0, -8.0, z0 - 21.5], [8.0, 8.0, z0 + 21.5], [m, 0.0, 0.0, 0.0, m * r * r / 4.0, m * r * r / 4.0, m * h * h / 12.0, 0.0, 0.0, 0.0]);
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let sq = [x, y, z].map(|a| b.un(UnOp::Square, a));
        let s = b.bin(BinOp::Add, sq[0], sq[1]);
        let s = b.bin(BinOp::Add, s, sq[2]);
        let len = b.un(UnOp::Sqrt, s);
        let ten = b.c(10.0);
        let ball = b.bin(BinOp::Sub, len, ten);
        let ball = Arc::new(b.tape(ball));
        let v = 4.0 / 3.0 * PI * 1000.0;
        let ball = (ball, [-11.0; 3], [11.0; 3], [v, 0.0, 0.0, 0.0, v * 20.0, v * 20.0, v * 20.0, 0.0, 0.0, 0.0]);
        for (shape, lo, hi, want) in [cylinder, ball] {
            // Asked for the mass alone it stops at the mass, as cheaply as before; asked for its moments, on.
            let mass = integrate(&shape, &[None], true, &Ctx::new(0.0, vec![]), lo, hi, Until::Accuracy(ACCURACY));
            let m = integrate(&shape, &[None], true, &Ctx::new(0.0, vec![]), lo, hi, Until::Moments(ACCURACY));
            assert!(mass.error(0) <= ACCURACY * want[0] && mass.boxes < m.boxes, "{} boxes for the mass, {} for its moments", mass.boxes, m.boxes);
            let scale = |k: usize| match k {
                0 | 4..=6 => want[k],
                1..=3 => (want[0] * want[k + 3]).sqrt(),
                _ => [(4, 5), (5, 6), (6, 4)].map(|(i, j)| (want[i] * want[j]).sqrt())[k - 7],
            };
            assert!(m.settled, "{m:?}");
            for k in 0..10 {
                assert!(m.lo[k] <= want[k] && want[k] <= m.hi[k], "{k}: {} not in [{}, {}]", want[k], m.lo[k], m.hi[k]);
                assert!(m.error(k) <= 1.01 * ACCURACY * scale(k), "{k}: {} against {}", m.error(k), ACCURACY * scale(k));
            }
        }
    }

    /// A cylinder r 7.3, h 41.9 about z = 5.1.
    fn cylinder() -> (Arc<Tape>, V3, V3) {
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let (xx, yy) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y));
        let rr = b.bin(BinOp::Add, xx, yy);
        let r = b.un(UnOp::Sqrt, rr);
        let rc = b.c(7.3);
        let side = b.bin(BinOp::Sub, r, rc);
        let ends = b.slab(z, 5.1, 41.9 / 2.0);
        let cyl = b.bin(BinOp::Max, side, ends);
        (Arc::new(b.tape(cyl)), [-8.0, -8.0, 5.1 - 21.5], [8.0, 8.0, 5.1 + 21.5])
    }

    /// A body's volume to an accuracy, measured with its moments, and then
    /// its moments to that accuracy, going on from where the volume's search
    /// stopped (`Kept`): each the same, field for field, as a search of its
    /// own (a cylinder, whose moments take more boxes than its volume).
    #[test]
    fn moments_go_on_from_the_volumes_search() {
        let (shape, lo, hi) = cylinder();
        let (ctx, gs) = (Ctx::new(0.0, vec![]), [None, None]);
        for a in [ACCURACY, 2e-3, 5e-4, 1e-4] {
            let own = |until| integrate(&shape, &gs, true, &ctx, lo, hi, until);
            let (volume, kept) = integrate_within(&shape, &gs, true, &ctx, lo, hi, Until::Accuracy(a), BOXES, None, true);
            assert_eq!(format!("{volume:?}"), format!("{:?}", own(Until::Accuracy(a))));
            let resumed = match kept {
                Some(Kept::Paused(p)) => Some(p),
                Some(Kept::Done(m)) => {
                    assert_eq!(format!("{m:?}"), format!("{:?}", own(Until::Moments(a))));
                    continue;
                }
                None => panic!("its search was not kept"),
            };
            let (moments, again) = integrate_within(&shape, &gs, true, &ctx, lo, hi, Until::Moments(a), BOXES, resumed, true);
            assert!(again.is_none() && moments.boxes > volume.boxes, "{} boxes, then {}", volume.boxes, moments.boxes);
            assert_eq!(format!("{moments:?}"), format!("{:?}", own(Until::Moments(a))));
        }
    }

    /// Two bodies whose shapes agree up to param numbering (`Param(0)` vs
    /// `Param(1)`, or swapped `(Param(0), Param(1))`) must not resume each
    /// other's kept volume search: `Paused` holds tapes indexing the original
    /// param slots.
    #[test]
    fn kept_search_does_not_collide_across_param_slots() {
        let cyl_param = |p_r: u32, p_h: Option<u32>| -> Arc<Tape> {
            let mut b = crate::test_fields::B::new();
            let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
            let (xx, yy) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y));
            let rr = b.bin(BinOp::Add, xx, yy);
            let r = b.un(UnOp::Sqrt, rr);
            let rp = b.v(Var::Param(p_r));
            let side = b.bin(BinOp::Sub, r, rp);
            let ends = match p_h {
                Some(k) => {
                    let hp = b.v(Var::Param(k));
                    let az = b.un(UnOp::Abs, z);
                    b.bin(BinOp::Sub, az, hp)
                }
                None => b.slab(z, 5.1, 41.9 / 2.0),
            };
            let cyl = b.bin(BinOp::Max, side, ends);
            Arc::new(b.tape(cyl))
        };
        let (lo, hi) = ([-8.0, -8.0, 5.1 - 21.5], [8.0, 8.0, 5.1 + 21.5]);
        let gs = [None, None];
        let a = 2e-3;
        // Order 1: A reads p0, B reads p1; volume(A), volume(B), inertia(A).
        let (ta, tb) = (cyl_param(0, None), cyl_param(1, None));
        let (ctx_a, ctx_b) = (Ctx::new(0.0, vec![7.3]), Ctx::new(0.0, vec![0.0, 7.3]));
        let _ = measured(&ta, &gs, true, &ctx_a, lo, hi, Until::Accuracy(a));
        let _ = measured(&tb, &gs, true, &ctx_b, lo, hi, Until::Accuracy(a));
        let ma = measured(&ta, &gs, true, &ctx_a, lo, hi, Until::Moments(a));
        assert_eq!(format!("{ma:?}"), format!("{:?}", integrate(&ta, &gs, true, &ctx_a, lo, hi, Until::Moments(a))));
        // Order 2: A reads p1, B reads p0; volume(A), volume(B), inertia(A).
        let _ = measured(&tb, &gs, true, &ctx_b, lo, hi, Until::Accuracy(a));
        let _ = measured(&ta, &gs, true, &ctx_a, lo, hi, Until::Accuracy(a));
        let mb = measured(&tb, &gs, true, &ctx_b, lo, hi, Until::Moments(a));
        assert_eq!(format!("{mb:?}"), format!("{:?}", integrate(&tb, &gs, true, &ctx_b, lo, hi, Until::Moments(a))));
        // Order 3: same params_read {0, 1}, swapped order in SSA.
        let (lo2, hi2) = ([-8.0, -8.0, -21.5], [8.0, 8.0, 21.5]);
        let (tc, td) = (cyl_param(0, Some(1)), cyl_param(1, Some(0)));
        let (ctx_c, ctx_d) = (Ctx::new(0.0, vec![7.3, 20.95]), Ctx::new(0.0, vec![20.95, 7.3]));
        let _ = measured(&tc, &gs, true, &ctx_c, lo2, hi2, Until::Accuracy(a));
        let _ = measured(&td, &gs, true, &ctx_d, lo2, hi2, Until::Accuracy(a));
        let mc = measured(&tc, &gs, true, &ctx_c, lo2, hi2, Until::Moments(a));
        assert_eq!(format!("{mc:?}"), format!("{:?}", integrate(&tc, &gs, true, &ctx_c, lo2, hi2, Until::Moments(a))));
    }

    /// A density the same at every point (a param, 2.7) is a factor of the
    /// body's own integrals: its moments' bounds those of the cylinder's
    /// shape times it, rounded outward, holding the exact ones; one that
    /// reads x is not.
    #[test]
    fn a_uniform_density_is_a_factor() {
        use std::f64::consts::PI;
        let (shape, lo, hi) = cylinder();
        let mut b = crate::test_fields::B::new();
        let p = b.v(Var::Param(0));
        let rho = b.tape(p);
        let x = b.v(Var::X);
        let graded = b.bin(BinOp::Add, p, x);
        let graded = b.tape(graded);
        let ctx = Ctx::new(0.0, vec![2.7]);
        let until = Until::Moments(ACCURACY);
        let m = uniform(&shape, &[Some(&rho), None], true, &ctx, lo, hi, until).expect("a param is the same everywhere");
        let own = measured(&shape, &[None, None], true, &ctx, lo, hi, until);
        let (r, h) = (7.3, 41.9);
        let v = PI * r * r * h;
        let mass = 2.7 * v;
        let want = [mass, 0.0, 0.0, 0.0, mass * r * r / 4.0, mass * r * r / 4.0, mass * h * h / 12.0, 0.0, 0.0, 0.0, v];
        for k in 0..11 {
            let w = if k < 10 { mul_o(I::point(2.7), I::new(own.lo[k], own.hi[k])) } else { I::new(own.lo[k], own.hi[k]) };
            assert!(m.lo[k] == w.lo && m.hi[k] == w.hi && m.lo[k] <= want[k] && want[k] <= m.hi[k], "{k}: {} not in [{}, {}]", want[k], m.lo[k], m.hi[k]);
            assert!(m.lo[k] <= m.sums[k] && m.sums[k] <= m.hi[k], "{k}: {m:?}");
        }
        assert!(m.error(0) <= ACCURACY * mass, "{m:?}");
        assert!(uniform(&shape, &[Some(&graded), None], true, &ctx, lo, hi, until).is_none());
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
