//! A solve over ranged params: its unknowns as certified roots (docs/IR.md,
//! Solving constraints, "Over ranged params").
//!
//! Over a box P of the ranged params its equations F(u, p) = 0 read (a
//! tolerance's band, or what an earlier solve derives), the Krawczyk
//! operator
//!
//! ```text
//! K(X, P) = c - Y F(c, P) + (I - Y J(X, P)) (X - c)
//! ```
//!
//! (c a point of the box X, Y any matrix: here the inverse of the Jacobian
//! at c; F(c, P) the equations' bounds at u = c over P, J(X, P) their
//! derivative ranges in u over X x P by the kernel's rules, which are
//! unbounded where a value may jump; every operation rounded outward) holds
//! every root in X of F(., p), for each p in P: for a root u, u = u - Y
//! F(u, p), and F(u, p) - F(c, p) is a matrix of J(X, P) times u - c (the
//! mean value theorem, row by row, along the segment from c to u, which
//! lies in X). Where K(X, P) lies in the interior of X, F(., p) has
//! exactly one root in X for each p in P: u - Y F(u, p) maps X into itself
//! (Brouwer's fixed-point theorem gives a fixed point, a root, once Y is
//! nonsingular), and the strict inclusion makes Y and every matrix in
//! J(X, P) nonsingular (Rump), so two roots cannot share X. That is the
//! certificate, made once over the whole box. Pieces of it and single
//! points then narrow the enclosure by K again, intersected with the box it
//! started from (the one root lies in both), which needs no new
//! certificate.
//!
//! At a point (no ranges) the same test shows a solution regular
//! (`Roots::regular`): K maps a box about it into its interior, so the
//! box holds one root and the Jacobian has full rank on all of it. With
//! more equations than unknowns Y is a left inverse of J at c, and the
//! argument holds for the roots of Y F alone: the box holds one root of
//! Y F and J has full rank on it, but no root of F is shown (u = 0 and
//! u = 1/8, weighted alike, have the one root 1/16 of Y F and none of F).
//! So a certificate that a root exists (`certify`, `enclose`) needs as
//! many equations as unknowns; `regular` gives the rank alone there, at
//! a root the solver measured.
//!
//! With t ranged over a step of a motion (`Roots::enclose`), the box holds
//! one root at every moment of it: where it also holds the root followed
//! to the step's start, the root at each moment of the step is that one,
//! moved along t (Follow). The roots in the box at its moments
//! form one curve: there is one root at each, and it moves continuously,
//! because the roots in a closed box of a continuous function form a
//! closed set.
//!
//! Over a step, what the equations read from solves on t moves along t
//! (`Moving`): F(c, w(t), t) is bounded at the step's middle moment, w in
//! its enclosure there, plus its rate along t over the box (dF/dt + dF/dw
//! w', w' bounded by `Roots::rate`) times t's distance from the middle (the
//! mean value theorem along t: the path stays in the box). That is first
//! order along the path, where the box alone (t and w ranged apart) is
//! first order in each: a rate read over a step is bounded by how far it
//! moves, not by how far its root does.

use crate::eval::{self, Ctx, Derive, Inputs, MAX_INPUTS, span_mid};
use crate::graph::{BinOp, Var};
use crate::tape::{Op, Tape};
use ink_core::interval::{self, I, add_o, mul_o, sub_o};
use std::collections::HashMap;
use std::sync::Mutex;

/// The equations do not read x, y, z: any box will do.
const B0: [I; 3] = [I { lo: 0.0, hi: 0.0 }; 3];

/// Inflations of the first guess tried before giving up (epsilon-inflation).
const INFLATIONS: usize = 16;

/// Narrowings of an enclosure per piece of the scope (each one K again).
const NARROWINGS: usize = 6;

fn width(a: I) -> f64 {
    a.hi - a.lo
}
fn meet(a: I, b: I) -> Option<I> {
    let m = I { lo: a.lo.max(b.lo), hi: a.hi.min(b.hi) };
    (m.lo <= m.hi).then_some(m)
}
fn clamp(v: f64, x: I) -> f64 {
    if v.is_finite() { v.clamp(x.lo, x.hi) } else { span_mid((x.lo, x.hi)) }
}

/// A param a solve on t sets, read over a span of t: its index, its
/// enclosure at the span's middle moment, and a bound on its rate along t
/// over the span (Roots::rate).
pub type Moving = (u32, I, I);

/// A solve's unknowns derived over ranged params: the one root of its
/// equations in `band` for every value of the params it reads in `over`.
pub struct Roots {
    unknowns: Vec<u32>,
    /// The ranged params its equations read, besides its unknowns.
    reads: Vec<u32>,
    tapes: Vec<Tape>,
    /// The moment it is certified at, where the equations read t.
    t: Option<f64>,
    /// The box it is certified over: each read param's span.
    over: Vec<(f64, f64)>,
    /// The root's enclosure over all of it.
    pub band: Vec<I>,
    /// Y of the certificate, for points where the Jacobian is singular
    /// (any Y keeps K an enclosure): n x m, m the equations.
    y0: Vec<Vec<f64>>,
    /// Enclosures and values found, by the spans, starting box and values
    /// they were derived from (searches derive one piece several times).
    cache: Mutex<HashMap<Vec<u64>, (Vec<I>, Vec<f64>)>>,
    /// Params read that move along t over ctx's span, at bounded rates.
    moving: Vec<Moving>,
}

impl Roots {
    /// Certify the root of `tapes` = 0 in `unknowns` over the spans of the
    /// ranged params of `ctx` they read, near the solution `solve` gives at
    /// the spans' middles (warm-started from ctx's values). Err: why not.
    pub fn certify(unknowns: &[usize], tapes: Vec<Tape>, uses_t: bool, ctx: &Ctx, solve: impl Fn(&mut [f64], f64) -> Result<(), String>) -> Result<Roots, String> {
        let unknowns: Vec<u32> = unknowns.iter().map(|&k| k as u32).collect();
        square(tapes.len(), unknowns.len())?;
        let mut reads: Vec<u32> = tapes
            .iter()
            .flat_map(|t| t.ssa.iter())
            .filter_map(|i| match i.op {
                Op::Var(Var::Param(k)) if ctx.p_spans.get(k).is_some() && !unknowns.contains(&k) => Some(k),
                _ => None,
            })
            .collect();
        reads.sort_unstable();
        reads.dedup();
        let over: Vec<(f64, f64)> = reads.iter().map(|&k| ctx.p_spans.get(k).expect("ranged")).collect();
        // Solved at the middle of the box, from the values.
        let mut mid = ctx.clone();
        for (&k, &s) in reads.iter().zip(&over) {
            mid.params[k as usize] = span_mid(s);
        }
        let mut params = mid.params.clone();
        solve(&mut params, ctx.t).map_err(|m| format!("at the middle of the ranges: {m}"))?;
        mid.params = params;
        let c0: Vec<f64> = unknowns.iter().map(|&k| mid.params[k as usize]).collect();
        mid.p_spans = mid.p_spans.with_derived(&unknowns.iter().zip(&c0).map(|(&k, &c)| (k, (c, c))).collect::<Vec<_>>());
        let mut r = Roots { unknowns, reads, tapes, t: uses_t.then_some(ctx.t), over, band: Vec::new(), y0: Vec::new(), cache: Mutex::new(HashMap::new()), moving: Vec::new() };
        let (c, y) = r.newton(&mid, c0);
        let Some(y) = y else {
            return Err("its Jacobian is singular at the solution at the middle of the ranges".into());
        };
        let Some(band) = r.contained(&mid, &c, &y, None, None) else {
            return Err(format!(
                "no single root is certified over the ranges (Krawczyk's test failed after {INFLATIONS} inflations). Two roots may lie close, the Jacobian may be singular in the range, or the ranges may be too wide"
            ));
        };
        (r.band, r.y0) = (band, y);
        // Narrowed over the whole box, as over any piece of it.
        let (band, _) = r.narrow(&mid, r.band.clone(), c);
        r.band = band;
        // Where the values lie in the ranges, the solution there (a numerical
        // method's, to its accuracy) must be the root certified, narrowed to
        // those values: a solve at the values on another branch would make
        // the ranges' claims about another configuration. (A check of the
        // model's consistency, not part of any proof.)
        let at_values = r.reads.iter().zip(&r.over).all(|(&k, s)| s.0 <= ctx.params[k as usize] && ctx.params[k as usize] <= s.1);
        if at_values {
            let mut point = mid.clone();
            for &k in &r.reads {
                let v = ctx.params[k as usize];
                point.p_spans.set_span(k, (v, v));
                point.params[k as usize] = v;
            }
            let solved: Vec<f64> = r.unknowns.iter().map(|&k| ctx.params[k as usize]).collect();
            let (e, _) = r.narrow(&point, r.band.clone(), solved.clone());
            if !e.iter().zip(&solved).all(|(e, &v)| (v - span_mid((e.lo, e.hi))).abs() <= 1e-6 * v.abs().max(1.0) + (e.hi - e.lo)) {
                return Err("its solution at the params' values is not the root certified over their ranges. The guesses or equations pick another root at those values".into());
            }
        }
        Ok(r)
    }

    /// ctx with the unknowns over `x` (their values `c`).
    fn at(&self, ctx: &Ctx, x: &[I], c: &[f64]) -> Ctx {
        let mut s = ctx.clone();
        for (i, &u) in self.unknowns.iter().enumerate() {
            s.p_spans.set_span(u, (x[i].lo, x[i].hi));
            s.params[u as usize] = c[i];
        }
        s
    }

    /// The equations' bounds over ctx's spans (first order in them), each
    /// case of a switch they read bounded apart (eval_interval_cases): an
    /// angle wrapped by `floor` at its seam is bounded as continuous. Over
    /// a span of t with params moving along it, met with their bounds along
    /// that path (the module's).
    fn residuals(&self, ctx: &Ctx) -> Vec<I> {
        let plain = self.tapes.iter().map(|t| eval::eval_interval_cases(t, ctx, B0));
        let Some(span) = ctx.t_span.filter(|_| !self.moving.is_empty() && self.moving.len() < MAX_INPUTS) else { return plain.collect() };
        let tm = span_mid(span);
        let mut mid = Ctx { t: tm, t_span: None, ..ctx.clone() };
        for &(k, x, _) in &self.moving {
            mid.p_spans.set_span(k, (x.lo, x.hi));
            mid.params[k as usize] = clamp(mid.params[k as usize], x);
        }
        let from_mid = sub_o(I::new(span.0, span.1), I::point(tm));
        let inputs = Inputs { t: true, params: self.moving.iter().map(|m| m.0).collect() };
        plain
            .zip(&self.tapes)
            .map(|(plain, tape)| {
                let (_, d) = eval::eval_interval_rates(tape, ctx, B0, &inputs);
                let rate = self.moving.iter().zip(&d[1..]).fold(d[0], |r, (m, dw)| add_o(r, mul_o(*dw, m.2)));
                interval::meet(plain, add_o(eval::eval_interval_cases(tape, &mid, B0), mul_o(rate, from_mid)))
            })
            .collect()
    }

    /// The equations' derivative ranges in the unknowns over ctx's spans
    /// (row per equation).
    fn jacobian(&self, ctx: &Ctx) -> Vec<Vec<I>> {
        let n = self.unknowns.len();
        let mut jac = vec![vec![I::ENTIRE; n]; self.tapes.len()];
        for (row, t) in jac.iter_mut().zip(&self.tapes) {
            for (ci, chunk) in self.unknowns.chunks(MAX_INPUTS).enumerate() {
                let (_, d) = eval::eval_interval_rates(t, ctx, B0, &Inputs { t: false, params: chunk.to_vec() });
                for (j, dj) in d.into_iter().enumerate() {
                    row[ci * MAX_INPUTS + j] = dj;
                }
            }
        }
        jac
    }

    /// The equations' values and Jacobian in the unknowns at `u` (ctx's
    /// values otherwise), in floats.
    fn point(&self, ctx: &Ctx, u: &[f64]) -> (Vec<f64>, Vec<Vec<f64>>) {
        let mut c = ctx.clone();
        for (i, &k) in self.unknowns.iter().enumerate() {
            c.params[k as usize] = u[i];
        }
        let mut f = Vec::new();
        let mut jac = Vec::new();
        for t in &self.tapes {
            let d = eval::eval_dual(t, &c, [0.0; 3]);
            f.push(d[0]);
            jac.push(self.unknowns.iter().map(|&k| d[5 + k as usize]).collect());
        }
        (f, jac)
    }

    /// Newton's method from `u` at ctx's values (Gauss-Newton with more
    /// equations): the root it reaches (a point: only a guess, which K
    /// checks) and the Jacobian's left inverse there.
    fn newton(&self, ctx: &Ctx, mut u: Vec<f64>) -> (Vec<f64>, Option<Vec<Vec<f64>>>) {
        for _ in 0..12 {
            let (f, j) = self.point(ctx, &u);
            let Some(y) = left_inverse(&j) else { return (u, None) };
            let step: Vec<f64> = y.iter().map(|row| row.iter().zip(&f).map(|(a, b)| a * b).sum()).collect();
            if !step.iter().all(|s| s.is_finite()) {
                return (u, Some(y));
            }
            let small = step.iter().zip(&u).all(|(s, v)| s.abs() <= 4.0 * f64::EPSILON * v.abs());
            for (v, s) in u.iter_mut().zip(&step) {
                *v -= s;
            }
            if small {
                break;
            }
        }
        let (_, j) = self.point(ctx, &u);
        let y = left_inverse(&j);
        (u, y)
    }

    /// K(X, P) with centre c (in X) and matrix y, ctx holding P; None where
    /// unbounded (a derivative range that is).
    fn krawczyk(&self, ctx: &Ctx, x: &[I], c: &[f64], y: &[Vec<f64>]) -> Option<Vec<I>> {
        let (n, rows) = (c.len(), self.tapes.len());
        let fc = self.residuals(&self.at(ctx, &c.iter().map(|&v| I::point(v)).collect::<Vec<_>>(), c));
        let jx = self.jacobian(&self.at(ctx, x, c));
        let d: Vec<I> = (0..n).map(|j| sub_o(x[j], I::point(c[j]))).collect();
        let mut k = Vec::with_capacity(n);
        for i in 0..n {
            let yf = (0..rows).fold(I::point(0.0), |s, j| add_o(s, mul_o(I::point(y[i][j]), fc[j])));
            let mut acc = sub_o(I::point(c[i]), yf);
            for j in 0..n {
                let m = (0..rows).fold(I::point(if i == j { 1.0 } else { 0.0 }), |m, l| sub_o(m, mul_o(I::point(y[i][l]), jx[l][j])));
                acc = add_o(acc, mul_o(m, d[j]));
            }
            if !(acc.lo.is_finite() && acc.hi.is_finite()) {
                return None;
            }
            k.push(acc);
        }
        Some(k)
    }

    /// A box X holding c (and `holding`, if given; inside `within`, if
    /// given) that K(X, P) maps into its interior, P ctx's spans, and K(X,
    /// P): the first guess c - Y F(c, P), inflated (a tenth of its width)
    /// until it is one, at most INFLATIONS times. Then for every p in P the
    /// equations have a root in X, the one root there, and K(X, P) holds it
    /// (the module's argument). None: no such box found.
    fn contained(&self, ctx: &Ctx, c: &[f64], y: &[Vec<f64>], within: Option<&[I]>, holding: Option<&[I]>) -> Option<Vec<I>> {
        let mut x = self.newton_box(ctx, c, y);
        for _ in 0..INFLATIONS {
            let mut xe = inflate(&x, c, holding);
            if let Some(w) = within {
                // c lies in `within` (the caller's), so each meet holds it.
                debug_assert!(c.iter().zip(w).all(|(&ci, wi)| wi.lo <= ci && ci <= wi.hi));
                xe = xe.iter().zip(w).map(|(&a, &b)| meet(a, b)).collect::<Option<Vec<I>>>()?;
            }
            let k = self.krawczyk(ctx, &xe, c, y)?;
            if k.iter().zip(&xe).all(|(k, x)| x.lo < k.lo && k.hi < x.hi) {
                return Some(k);
            }
            x = k;
        }
        None
    }

    /// c - Y F(c, P), P ctx's spans: where Newton's step from c goes over
    /// them.
    fn newton_box(&self, ctx: &Ctx, c: &[f64], y: &[Vec<f64>]) -> Vec<I> {
        let fc = self.residuals(&self.at(ctx, &c.iter().map(|&v| I::point(v)).collect::<Vec<_>>(), c));
        (0..c.len()).map(|i| sub_o(I::point(c[i]), (0..fc.len()).fold(I::point(0.0), |s, j| add_o(s, mul_o(I::point(y[i][j]), fc[j]))))).collect()
    }

    /// Whether the equations' rates along the unknowns over x (at c) and
    /// along t over ctx's span are bounded: nothing they read switches
    /// there (a select or floor, whose rate the kernel leaves unbounded
    /// across its jump).
    fn smooth(&self, ctx: &Ctx, x: &[I], c: &[f64]) -> bool {
        let at = self.at(ctx, x, c);
        let bounded = |d: &I| d.lo.is_finite() && d.hi.is_finite();
        let along_t = Inputs { t: true, params: Vec::new() };
        self.jacobian(&at).iter().flatten().all(bounded)
            && (at.t_span.is_none() || self.tapes.iter().all(|t| eval::eval_interval_rates(t, &at, B0, &along_t).1.iter().all(bounded)))
    }

    /// The root's enclosure over ctx's spans of the params it reads,
    /// narrowed from `x` (which holds it for all of them) by K, and the
    /// root at ctx's values (by Newton's method from `c`, inside the
    /// enclosure).
    ///
    /// Where K over x does not narrow it (x straddles a kink of the
    /// equations, where their slopes over x make I - Y J too wide), a box
    /// of its own inside x is tried (`contained`): it holds a root for every
    /// value of the spans, and x holds only one (it lies in the band
    /// certified), so it holds that one.
    fn narrow(&self, ctx: &Ctx, mut x: Vec<I>, mut c: Vec<f64>) -> (Vec<I>, Vec<f64>) {
        for i in 0..NARROWINGS {
            let (u, y) = self.newton(ctx, c.iter().zip(&x).map(|(&v, &xi)| clamp(v, xi)).collect());
            c = u.iter().zip(&x).map(|(&v, &xi)| clamp(v, xi)).collect();
            let y = y.unwrap_or_else(|| self.y0.clone());
            let Some(k) = self.krawczyk(ctx, &x, &c, &y) else { break };
            // Both hold the root; never trust rounding to keep them meeting.
            let Some(m) = k.iter().zip(&x).map(|(&k, &x)| meet(k, x)).collect::<Option<Vec<I>>>() else { break };
            let gain = m.iter().zip(&x).any(|(a, b)| width(*a) < 0.5 * width(*b));
            x = m;
            if !gain {
                // Its centre in x (K's mean value form needs it in the box).
                c = c.iter().zip(&x).map(|(&v, &xi)| clamp(v, xi)).collect();
                if i == 0
                    && let Some(k) = self.contained(ctx, &c, &y, Some(&x), None)
                    && let Some(m) = k.iter().zip(&x).map(|(&k, &x)| meet(k, x)).collect::<Option<Vec<I>>>()
                {
                    x = m;
                    continue;
                }
                break;
            }
        }
        let c = c.iter().zip(&x).map(|(&v, &xi)| clamp(v, xi)).collect();
        (x, c)
    }
}

impl Roots {
    /// Equations to find one root of, at moment t, with no ranges of their
    /// own (no reads, no band), reading params `moving` along t.
    fn of(unknowns: &[usize], tapes: Vec<Tape>, t: f64, moving: &[Moving]) -> Roots {
        Roots {
            unknowns: unknowns.iter().map(|&k| k as u32).collect(),
            reads: Vec::new(),
            tapes,
            t: Some(t),
            over: Vec::new(),
            band: Vec::new(),
            y0: Vec::new(),
            cache: Mutex::new(HashMap::new()),
            moving: moving.to_vec(),
        }
    }

    /// ctx with the unknowns as spans at their values, which K sets
    /// (Roots::at).
    fn spanned(&self, ctx: &Ctx) -> Ctx {
        let at: Vec<(u32, (f64, f64))> = self.unknowns.iter().map(|&k| (k, (ctx.params[k as usize], ctx.params[k as usize]))).collect();
        Ctx { p_spans: ctx.p_spans.clone().with_derived(&at), ..ctx.clone() }
    }

    /// The root about ctx's values of the unknowns, over ctx's span of t
    /// and spans of params (at one moment and value where it has none):
    /// Newton's method from them at ctx's point, c, and K(X, P) of a box X
    /// holding c (and `from`, if given) that K maps into its interior. So
    /// at every moment and value of the spans the equations have one root
    /// in X, which K(X, P) holds, and the Jacobian has full rank on X (the
    /// module's argument); `moving`, params among ctx's spans that move
    /// along t at bounded rates. Err: why not, as where the equations are
    /// not as many as the unknowns.
    pub fn enclose(unknowns: &[usize], tapes: Vec<Tape>, ctx: &Ctx, from: Option<&[I]>, moving: &[Moving]) -> Result<(Vec<f64>, Vec<I>), String> {
        square(tapes.len(), unknowns.len())?;
        Roots::contain(unknowns, tapes, ctx, from, moving)
    }

    /// `enclose` for any number of equations: with more than the unknowns,
    /// the box about the root of Y F (the module's argument).
    fn contain(unknowns: &[usize], tapes: Vec<Tape>, ctx: &Ctx, from: Option<&[I]>, moving: &[Moving]) -> Result<(Vec<f64>, Vec<I>), String> {
        let r = Roots::of(unknowns, tapes, ctx.t, moving);
        let ctx = &r.spanned(ctx);
        let (c, y) = r.newton(ctx, unknowns.iter().map(|&k| ctx.params[k]).collect());
        let Some(y) = y else { return Err("its Jacobian in its unknowns is singular at the root".into()) };
        if let Some(k) = r.contained(ctx, &c, &y, None, from) {
            return Ok((c, k));
        }
        // Why: where nothing switches over the first box tried, two roots
        // meet (or the one there ends); else the equations jump there, or
        // their bounds cannot show them continuous.
        let why = match r.smooth(ctx, &inflate(&r.newton_box(ctx, &c, &y), &c, from), &c) {
            true => "a dead point, where two roots meet and the Jacobian in its unknowns is singular",
            false => {
                "its equations switch there (a select or floor reading t or its unknowns), and their bounds do not show them continuous \
                      (a jump; more than three switches at one moment; a switch in its unknowns)"
            }
        };
        Err(format!("no single root with a Jacobian of full rank is certified about the solution (Krawczyk's test failed after {INFLATIONS} inflations): {why}"))
    }

    /// Whether `tapes` (a solve's equations in `unknowns`) are continuous
    /// along t over ctx's span for every value of the unknowns in `tube`:
    /// each one's rate along t, plus its rate along each param `moving`
    /// along t times that param's rate, bounded, so it moves by at most
    /// that times the time elapsed. A box holding one root at every moment
    /// of a span (`enclose`) is a path of that root only where they are:
    /// where they jump, roots either side of the jump may share the box
    /// (the box about a crank's root over a step where four wrapped angles
    /// jump together held the roots 1.329 before and 1.183 after). A jump,
    /// a seam t moves (eval::scope_run joins its sides with unbounded
    /// rates, even where they agree there) and an
    /// infinite slope alike leave the rate unbounded: not shown continuous.
    pub fn continuous(unknowns: &[usize], tapes: &[Tape], ctx: &Ctx, tube: &[I], moving: &[Moving]) -> bool {
        if ctx.t_span.is_none() {
            return true;
        }
        let r = Roots::of(unknowns, tapes.to_vec(), ctx.t, moving);
        let c: Vec<f64> = unknowns.iter().zip(tube).map(|(&k, &x)| clamp(ctx.params[k], x)).collect();
        let at = r.at(&r.spanned(ctx), tube, &c);
        let ks: Vec<u32> = moving.iter().map(|m| m.0).collect();
        let first = ks.len().min(MAX_INPUTS - 1);
        tapes.iter().all(|tape| {
            let (_, d) = eval::eval_interval_rates(tape, &at, B0, &Inputs { t: true, params: ks[..first].to_vec() });
            let mut along = d[1..].to_vec();
            for chunk in ks[first..].chunks(MAX_INPUTS) {
                along.extend(eval::eval_interval_rates(tape, &at, B0, &Inputs { t: false, params: chunk.to_vec() }).1);
            }
            let total = moving.iter().zip(&along).fold(d[0], |s, (m, dw)| add_o(s, mul_o(*dw, m.2)));
            total.lo.is_finite() && total.hi.is_finite()
        })
    }

    /// The one root in `tube` (a box `enclose` certified over spans holding
    /// ctx's) at ctx: Newton's method from ctx's values, clamped into it,
    /// narrowed by K (Roots::narrow); its value and enclosure.
    pub fn inside(unknowns: &[usize], tapes: Vec<Tape>, ctx: &Ctx, tube: &[I]) -> Result<(Vec<f64>, Vec<I>), String> {
        let mut r = Roots::of(unknowns, tapes, ctx.t, &[]);
        let ctx = &r.spanned(ctx);
        let c: Vec<f64> = unknowns.iter().zip(tube).map(|(&k, &x)| clamp(ctx.params[k], x)).collect();
        r.y0 = left_inverse(&r.point(ctx, &c).1).ok_or("its Jacobian in its unknowns is singular in the box certified to hold one root")?;
        let (x, c) = r.narrow(ctx, tube.to_vec(), c);
        Ok((c, x))
    }

    /// The regular root of `tapes` = 0 the unknowns' values at ctx (a
    /// point: no ranges, no span of t) stand for: Newton's method from
    /// them, the root it reaches, where a box about it is one K maps into
    /// its interior (`enclose`). With more equations than unknowns, Gauss-
    /// Newton's point, the Jacobian of full rank about it, but F's being 0
    /// there not shown (the module's argument): a root only as far as the
    /// solver measured one there. Err: why not.
    pub fn regular(unknowns: &[usize], tapes: Vec<Tape>, ctx: &Ctx) -> Result<Vec<f64>, String> {
        Roots::contain(unknowns, tapes, ctx, None, &[]).map(|r| r.0)
    }

    /// A bound on the rate along t, over ctx's span of t, of the root of
    /// `tapes` in `unknowns` that ctx's spans of them hold (one at each
    /// moment, its Jacobian regular): by the implicit function theorem J
    /// x' = -(dF/dt + dF/dw w'), J = dF/dx, each over the box, w the params
    /// `moving` along t at their rates; every solution of that interval
    /// system enclosed (`linear`). None: not bounded (a switch, or J not
    /// shown regular over the box).
    pub fn rate(unknowns: &[usize], tapes: &[Tape], ctx: &Ctx, moving: &[Moving]) -> Option<Vec<I>> {
        let n = unknowns.len();
        if tapes.len() != n {
            return None;
        }
        let along: Vec<u32> = unknowns.iter().map(|&k| k as u32).chain(moving.iter().map(|m| m.0)).collect();
        let (mut jac, mut rhs) = (Vec::with_capacity(n), Vec::with_capacity(n));
        for tape in tapes {
            let mut d = Vec::with_capacity(along.len());
            let mut dt = I::ENTIRE;
            for chunk in along.chunks(MAX_INPUTS - 1) {
                let (_, r) = eval::eval_interval_rates(tape, ctx, B0, &Inputs { t: true, params: chunk.to_vec() });
                dt = r[0];
                d.extend_from_slice(&r[1..]);
            }
            let total = moving.iter().zip(&d[n..]).fold(dt, |s, (m, dw)| add_o(s, mul_o(*dw, m.2)));
            rhs.push(I::new(-total.hi, -total.lo));
            jac.push(d[..n].to_vec());
        }
        linear(&jac, &rhs)
    }
}

/// Every solution of a z = b for each matrix in `a` and vector in `b` (n x
/// n), enclosed: with y the inverse of a's middle, z = y b + (I - y a) z,
/// so where |I - y a| < rho < 1 (the infinity norm) every a in it is
/// regular and |z| <= |y b| / (1 - rho); K, from that box, encloses z
/// again. None: a not shown regular so.
fn linear(a: &[Vec<I>], b: &[I]) -> Option<Vec<I>> {
    let n = b.len();
    let y = invert(&a.iter().map(|r| r.iter().map(I::mid).collect()).collect::<Vec<_>>())?;
    let yb: Vec<I> = (0..n).map(|i| (0..n).fold(I::point(0.0), |s, l| add_o(s, mul_o(I::point(y[i][l]), b[l])))).collect();
    let m: Vec<Vec<I>> =
        (0..n).map(|i| (0..n).map(|j| (0..n).fold(I::point(if i == j { 1.0 } else { 0.0 }), |s, l| sub_o(s, mul_o(I::point(y[i][l]), a[l][j])))).collect()).collect();
    let rho = m.iter().map(|r| r.iter().fold(0.0, |s, &e| interval::add_up(s, interval::mag(e)))).fold(0.0, f64::max);
    let beta = yb.iter().map(|&e| interval::mag(e)).fold(0.0, f64::max);
    let r = interval::bin(BinOp::Div, I::point(beta), I::point(sub_o(I::point(1.0), I::point(rho)).lo)).0.hi;
    if !(rho < 1.0 && r.is_finite()) {
        return None;
    }
    let mut z = vec![I::new(-r, r); n];
    for _ in 0..2 {
        z = (0..n).map(|i| interval::meet(z[i], (0..n).fold(yb[i], |s, j| add_o(s, mul_o(m[i][j], z[j]))))).collect();
    }
    Some(z)
}

impl Derive for Roots {
    fn owns(&self) -> &[u32] {
        &self.unknowns
    }

    /// -J_u^-1 J_p at ctx's values (the implicit function theorem).
    fn rates(&self, ctx: &Ctx) -> Vec<Vec<f64>> {
        let n = self.unknowns.len();
        let mut ju = Vec::new();
        let mut jp = Vec::new();
        for t in &self.tapes {
            let d = eval::eval_dual(t, ctx, [0.0; 3]);
            ju.push(self.unknowns.iter().map(|&k| d[5 + k as usize]).collect::<Vec<f64>>());
            jp.push(self.reads.iter().map(|&k| d[5 + k as usize]).collect::<Vec<f64>>());
        }
        let Some(y) = invert(&ju) else { return vec![vec![0.0; self.reads.len()]; n] };
        (0..n).map(|i| (0..self.reads.len()).map(|j| -(0..n).map(|l| y[i][l] * jp[l][j]).sum::<f64>()).collect()).collect()
    }

    fn reads(&self) -> &[u32] {
        &self.reads
    }

    fn derive(&self, ctx: &mut Ctx) {
        let p: Vec<(f64, f64)> = self.reads.iter().map(|&k| ctx.p_spans.get(k).unwrap_or((ctx.params[k as usize], ctx.params[k as usize]))).collect();
        // Outside the box certified (or at another moment), no root is
        // known: every value.
        let certified = p.iter().zip(&self.over).all(|(s, o)| o.0 <= s.0 && s.1 <= o.1) && self.t.is_none_or(|t| t == ctx.t && ctx.t_span.is_none());
        if !certified {
            for &u in &self.unknowns {
                ctx.p_spans.set_span(u, (f64::NEG_INFINITY, f64::INFINITY));
            }
            return;
        }
        // The spans in ctx hold the root over p (ParamSpans), as the band
        // does over all of it.
        let x: Vec<I> = self.unknowns.iter().zip(&self.band).map(|(&u, &b)| ctx.p_spans.get(u).and_then(|(lo, hi)| meet(I { lo, hi }, b)).unwrap_or(b)).collect();
        let values: Vec<f64> = self.reads.iter().map(|&k| ctx.params[k as usize]).collect();
        let key: Vec<u64> = p.iter().flat_map(|s| [s.0, s.1]).chain(x.iter().flat_map(|i| [i.lo, i.hi])).chain(values.iter().copied()).map(f64::to_bits).collect();
        let hit = self.cache.lock().unwrap().get(&key).cloned();
        let (x, c) = hit.unwrap_or_else(|| {
            let c: Vec<f64> = self.unknowns.iter().map(|&u| ctx.params[u as usize]).collect();
            let r = self.narrow(ctx, x, c);
            let mut cache = self.cache.lock().unwrap();
            if cache.len() > 1 << 16 {
                cache.clear();
            }
            cache.insert(key, r.clone());
            r
        });
        for ((&u, xi), ci) in self.unknowns.iter().zip(x).zip(c) {
            ctx.p_spans.set_span(u, (xi.lo, xi.hi));
            ctx.params[u as usize] = ci;
        }
    }
}

/// x holding c (which K's mean value form needs) and `holding`, if given,
/// widened by a tenth of its width each way, rounded outward.
fn inflate(x: &[I], c: &[f64], holding: Option<&[I]>) -> Vec<I> {
    x.iter()
        .zip(c)
        .enumerate()
        .map(|(i, (xi, &ci))| {
            let h = holding.map_or(I::point(ci), |h| h[i]);
            let (lo, hi) = (xi.lo.min(ci).min(h.lo), xi.hi.max(ci).max(h.hi));
            let e = 0.1 * (hi - lo) + 1e-12 * lo.abs().max(hi.abs()) + 1e-300;
            I::new((lo - e).next_down(), (hi + e).next_up())
        })
        .collect()
}

/// Err unless there are as many equations as unknowns: a root of a
/// weighted sum of more equations is not a root of each (the module's
/// argument).
fn square(equations: usize, unknowns: usize) -> Result<(), String> {
    match equations == unknowns {
        true => Ok(()),
        false => Err(format!("it has {equations} equations for {unknowns} unknowns. A certified root needs as many equations as unknowns (least squares has no proven range)")),
    }
}

/// A left inverse of a m x n matrix of full column rank, m >= n: its
/// inverse where square, else (a^T a)^-1 a^T; None where singular.
fn left_inverse(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let (m, n) = (a.len(), a.first().map_or(0, Vec::len));
    if m == n {
        return invert(a);
    }
    let ata: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| (0..m).map(|r| a[r][i] * a[r][j]).sum()).collect()).collect();
    let inv = invert(&ata)?;
    Some((0..n).map(|i| (0..m).map(|r| (0..n).map(|j| inv[i][j] * a[r][j]).sum()).collect()).collect())
}

/// The inverse of a small square matrix (Gauss-Jordan, partial pivoting);
/// None where singular.
fn invert(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut m: Vec<Vec<f64>> = a.iter().enumerate().map(|(i, r)| r.iter().copied().chain((0..n).map(|j| if i == j { 1.0 } else { 0.0 })).collect()).collect();
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| m[i][c].abs().total_cmp(&m[j][c].abs()))?;
        if !(m[p][c].abs() > 1e-300) {
            return None;
        }
        m.swap(c, p);
        let d = m[c][c];
        for v in m[c].iter_mut() {
            *v /= d;
        }
        for r in 0..n {
            if r != c {
                let f = m[r][c];
                if f != 0.0 {
                    for k in 0..2 * n {
                        m[r][k] -= f * m[c][k];
                    }
                }
            }
        }
    }
    let inv: Vec<Vec<f64>> = m.into_iter().map(|r| r[n..].to_vec()).collect();
    inv.iter().flatten().all(|v| v.is_finite()).then_some(inv)
}
