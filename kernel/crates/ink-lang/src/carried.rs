//! Measures of a body that only moves rigidly over a motion (docs/IR.md,
//! Motion): a thickness asked once, where it is the same at every moment;
//! an integral or a mean taken once in the body's own frame and carried by
//! its placement, an expression of t ranged over the span.
//!
//! A body moves rigidly where its world field is a shape that reads no t
//! placed (an instance's placement) or warped (a `move`, a `rotate-z`) by
//! steps each proven an isometry over the span (rigid::Isometry: a
//! rotation exactly, its translation finite). Its shape there is one
//! question, answered once. A shape that changes with t (it deforms, or is
//! made of parts that move apart) has no such form: it is asked as before,
//! and is not settled between the moments read.

use crate::claims::{self, Check, ClaimOptions, Evidence, Status};
use crate::motion::{Motion, Stepped};
use crate::view::ModelView;
use crate::{Claim, InkError};
use ink_core::eval::{Ctx, eval_interval_scope, span_mid};
use ink_core::graph::{BinOp, Graph, Node, UnOp, Var};
use ink_core::interval::{self, I, add_o};
use ink_core::tape::Tape;
use ink_kernel::query;
use ink_kernel::ranged::{Comparison, Status as S};
use ink_kernel::rigid::Isometry;

type M3 = [[u32; 3]; 3];

/// A world field as a shape that reads no t, `shape`, moved rigidly at
/// every moment of a scope: p = r q + t takes its points to the world's.
struct Rigidly {
    shape: u32,
    r: M3,
    t: [u32; 3],
    /// The field reads no t: it is `shape`, where it is.
    still: bool,
}

/// The most placements and warps read through to a shape that reads no t.
const DEPTH: usize = 16;

/// The most pieces of the span an integral's or a mean's bounds are asked
/// over (`halved`): the same work on every machine.
const PIECES: usize = 4096;

fn xyz(x: Var) -> bool {
    matches!(x, Var::X | Var::Y | Var::Z)
}

fn zero(g: &Graph, n: u32) -> bool {
    g.as_const(n) == Some(0.0)
}

// What is built of a placement's and an integrand's numbers is real
// arithmetic, as the field it stands for is where it reads x, y, z or t
// (Graph::frame_bin): a sum of fixed numbers is that sum, bounded, never
// its rounded float, so (1e16 + 1) - 1e16 is 1. Exact numbers fold
// exactly.
fn mul(g: &mut Graph, a: u32, b: u32) -> u32 {
    g.frame_bin(BinOp::Mul, a, b)
}

fn add(g: &mut Graph, a: u32, b: u32) -> u32 {
    g.frame_bin(BinOp::Add, a, b)
}

fn sub(g: &mut Graph, a: u32, b: u32) -> u32 {
    g.frame_bin(BinOp::Sub, a, b)
}

fn neg(g: &mut Graph, a: u32) -> u32 {
    g.un_real(UnOp::Neg, a)
}

/// An expression affine in x, y and z: its coefficients of them and its
/// constant term (None: 0), each reading none of x, y and z.
type Affine = [Option<u32>; 4];

/// `n` read as affine in x, y and z, exactly as it is written (a sum, a
/// difference, a product or a quotient by what reads none of them), in
/// real numbers (`mul`), each part of it that reads none of them a number
/// of its own (`leaves`) and each divisor one too (`divisors`); None where
/// it is not so written. The reading is n's value only where every one of
/// them is finite and no divisor is 0 (`finite`): a term times a pole, less
/// the same term, is no 0 where the pole is infinite.
fn affine(g: &mut Graph, n: u32, leaves: &mut Vec<u32>, divisors: &mut Vec<u32>) -> Option<Affine> {
    if !g.depends_on(n, xyz) {
        leaves.push(n);
        return Some([None, None, None, Some(n)]);
    }
    let both = |g: &mut Graph, a: Affine, b: Affine, op: fn(&mut Graph, u32, u32) -> u32, lone: fn(&mut Graph, u32) -> u32| -> Affine {
        std::array::from_fn(|i| match (a[i], b[i]) {
            (Some(x), Some(y)) => Some(op(g, x, y)),
            (Some(x), None) => Some(x),
            (None, Some(y)) => Some(lone(g, y)),
            (None, None) => None,
        })
    };
    match g.nodes[n as usize] {
        Node::Var(Var::X) => Some([Some(g.constant(1.0)), None, None, None]),
        Node::Var(Var::Y) => Some([None, Some(g.constant(1.0)), None, None]),
        Node::Var(Var::Z) => Some([None, None, Some(g.constant(1.0)), None]),
        Node::Tag(_, a) => affine(g, a, leaves, divisors),
        Node::Un(UnOp::Neg, a) => Some(affine(g, a, leaves, divisors)?.map(|e| e.map(|e| neg(g, e)))),
        Node::Bin(BinOp::Add, a, b) => {
            let (a, b) = (affine(g, a, leaves, divisors)?, affine(g, b, leaves, divisors)?);
            Some(both(g, a, b, add, |_, y| y))
        }
        Node::Bin(BinOp::Sub, a, b) => {
            let (a, b) = (affine(g, a, leaves, divisors)?, affine(g, b, leaves, divisors)?);
            Some(both(g, a, b, sub, neg))
        }
        Node::Bin(BinOp::Mul, a, b) if !g.depends_on(a, xyz) => {
            leaves.push(a);
            Some(affine(g, b, leaves, divisors)?.map(|e| e.map(|e| mul(g, a, e))))
        }
        Node::Bin(BinOp::Mul, a, b) if !g.depends_on(b, xyz) => {
            leaves.push(b);
            Some(affine(g, a, leaves, divisors)?.map(|e| e.map(|e| mul(g, e, b))))
        }
        Node::Bin(BinOp::Div, a, b) if !g.depends_on(b, xyz) => {
            divisors.push(b);
            Some(affine(g, a, leaves, divisors)?.map(|e| e.map(|e| g.frame_bin(BinOp::Div, e, b))))
        }
        _ => None,
    }
}

/// `affine` of `n`, where its reading is its value at every point of the
/// scope `over`: every number it is made of finite there, and no divisor
/// 0 (bounds rounded outward).
fn affine_over(g: &mut Graph, n: u32, over: &Ctx) -> Option<Affine> {
    let (mut leaves, mut divisors) = (Vec::new(), Vec::new());
    let a = affine(g, n, &mut leaves, &mut divisors)?;
    let apart = |d: &u32| {
        let x = eval_interval_scope(&Tape::new(g, *d), over, [I::point(0.0); 3], None);
        x.lo > 0.0 || x.hi < 0.0
    };
    (ink_kernel::rigid::finite_over(g, leaves.iter().chain(&divisors), over) && divisors.iter().all(apart)).then_some(a)
}

/// `f` as a shape that reads no t moved rigidly over the scope `over`
/// (`Rigidly`): through its placement (Doc::placement_of), a warp of it
/// (Doc::warps) whose coordinates are A p + b with A a rotation (the shape
/// moved by p = A^T (q - b)), or a name on it, each step proven an
/// isometry over the scope. None: it is not shown so.
fn rigidly(v: &mut ModelView, f: u32, over: &Ctx, depth: usize) -> Option<Rigidly> {
    if !v.doc.g.depends_on(f, |x| x == Var::T) {
        let g = &mut v.doc.g;
        let (z, o) = (g.constant(0.0), g.constant(1.0));
        return Some(Rigidly { shape: f, r: [[o, z, z], [z, o, z], [z, z, o]], t: [z; 3], still: true });
    }
    let depth = depth.checked_sub(1)?;
    let (inner, ra, ta) = if let Some(p) = v.doc.placement_of(f) {
        (p.shape, p.at.r, p.at.t)
    } else if let Some(&(inner, map)) = v.doc.warps.get(&f) {
        let g = &mut v.doc.g;
        let mut rows = [[None; 4]; 3];
        for (row, &n) in rows.iter_mut().zip(&map) {
            *row = affine_over(g, n, over)?;
        }
        let z = g.constant(0.0);
        // p = A^T q - A^T b.
        let r: M3 = std::array::from_fn(|i| std::array::from_fn(|j| rows[j][i].unwrap_or(z)));
        let t: [u32; 3] = std::array::from_fn(|i| {
            let mut s = z;
            for (j, row) in rows.iter().enumerate() {
                if let Some(b) = row[3] {
                    let m = mul(g, r[i][j], b);
                    s = add(g, s, m);
                }
            }
            neg(g, s)
        });
        (inner, r, t)
    } else if let Node::Tag(_, a) = v.doc.g.nodes[f as usize] {
        return rigidly(v, a, over, depth);
    } else {
        return None;
    };
    Isometry::over(&v.doc.g, &ra, &ta, over)?;
    let c = rigidly(v, inner, over, depth)?;
    // p = ra (r q + t) + ta.
    let g = &mut v.doc.g;
    let dot = |g: &mut Graph, a: [u32; 3], b: [u32; 3]| {
        let mut s = g.constant(0.0);
        for k in 0..3 {
            let m = mul(g, a[k], b[k]);
            s = add(g, s, m);
        }
        s
    };
    let col = |m: &M3, j: usize| [m[0][j], m[1][j], m[2][j]];
    let r: M3 = std::array::from_fn(|i| std::array::from_fn(|j| dot(g, ra[i], col(&c.r, j))));
    let t: [u32; 3] = std::array::from_fn(|i| {
        let d = dot(g, ra[i], c.t);
        add(g, d, ta[i])
    });
    Some(Rigidly { shape: c.shape, r, t, still: false })
}

/// A claim over a motion `span` of a body that only moves rigidly, where
/// what it measures is the same of the body moved so or is carried by its
/// placement: a thickness every way (Package::rigid_metric), asked once at
/// a moment of the span; an integral or a mean of what is affine in x, y
/// and z (a centre of mass, a mass), its own frame's integrals measured
/// once and carried by its placement with t ranged over the span
/// (`moments`). None: not so (a body that deforms, an integrand of other
/// fields of the world or of t, a bound that moves, what a package
/// solves), or undecided there:
/// it is asked as any claim over a motion is.
pub(crate) fn over_span(v: &mut ModelView, claim: &Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions) -> Result<Option<Check>, InkError> {
    if bound.is_some() || !(span.0 < span.1) {
        return Ok(None);
    }
    // What a package solves is asked as before: one on t is set moment by moment, and another's error
    // band moves the claim (claims::decide's estimated).
    let roots = v.doc.claim_roots(claim);
    if !v.doc.solves_read(&roots).is_empty() {
        return Ok(None);
    }
    let over = Ctx { t: span_mid(span), t_span: Some(span), ..v.ctx.clone() };
    // A moment of the span: the view's where it is in it.
    let moment = if span.0 <= v.ctx.t && v.ctx.t <= span.1 { v.ctx.t } else { span.0 };
    let c = match claim {
        Claim::Metric { metric, .. } if let Some(f) = metric.shape.filter(|_| v.doc.registry.rigid_metric(metric.metric)) => {
            if rigidly(v, f, &over, DEPTH).is_none() {
                return Ok(None);
            }
            v.set_moment(moment)?;
            let c = claims::decide(v, claim.clone(), opts)?;
            let status = match c.status {
                Status::Holds => S::Proven,
                Status::Fails => S::Disproven,
                Status::Unknown => S::Unknown,
            };
            let t = (c.status == Status::Fails).then_some(moment);
            let m = Motion { t, steps: 1, ..Motion::none(status, span) };
            Check { motion: Some(Stepped::rigid(m)), note: c.note.or(Some(ONCE)), ..c }
        }
        Claim::Integral { parts, cmp, bound, .. } => {
            let parts: Vec<(u32, Option<u32>, Vec<u32>)> = parts.iter().map(|(s, gs)| (*s, None, gs.clone())).collect();
            let c = Comparison { cmp: *cmp, bound: *bound, tol: claims::claim_accuracy(v, claim, opts.accuracy) };
            let Some(c) = moments(v, &parts, false, c, span, moment, &over, opts.res) else { return Ok(None) };
            c
        }
        Claim::Mean { parts, cmp, bound, .. } => {
            let parts: Vec<(u32, Option<u32>, Vec<u32>)> = parts.iter().map(|(s, w, hs)| (*s, Some(*w), hs.clone())).collect();
            let c = Comparison { cmp: *cmp, bound: *bound, tol: claims::claim_accuracy(v, claim, opts.accuracy) };
            let Some(c) = moments(v, &parts, true, c, span, moment, &over, opts.res) else { return Ok(None) };
            c
        }
        _ => return Ok(None),
    };
    let k = v.t_param();
    let mut ranged = c.ranged.clone();
    ranged.push((k, span));
    Ok(Some(Check { ranged, ..c }))
}

const ONCE: &str = "the body only moves rigidly over the motion (each placement an isometry at every moment), and what is measured is the same of it moved so: asked once, at a moment of the span";

/// What a part's integral is made of: a function of t, times one of the
/// part's integrals in its own frame (its shape's, and the integrand's
/// index there).
type Term = (u32, (usize, usize));

/// An integral claim's (`mean`: a mean's) parts, each its world shape, its
/// weight (None: 1) and integrands, over a motion `span` (`over`): in each
/// part's own frame (`rigidly`), its weight there reads no t and each
/// integrand h is affine in x, y and z (h = c0 + c.p), so the integral of
/// w h over the part at a moment is (c0 + c.t) W + (r^T c).M, W and M its
/// own frame's integrals of w and w q, the same at every moment. A part
/// that does not move takes any integrand that reads no t as it is. Those
/// integrals are measured once, to proven bounds (refined as an integral
/// claim is: claims::measure_steps), and the claim's value bounded over
/// pieces of the span by interval bounds of the terms of t (`halved`).
/// Holds where it holds over every piece; fails over a piece where it
/// fails at every moment of it (its middle the moment given). None: a part
/// not so, a shape not proven bounded, what it reads ranged (its integrals
/// are measured at one value), or undecided.
#[allow(clippy::too_many_arguments)]
fn moments(v: &mut ModelView, parts: &[(u32, Option<u32>, Vec<u32>)], mean: bool, c: Comparison, span: (f64, f64), moment: f64, over: &Ctx, res: Option<f64>) -> Option<Check> {
    // Each shape measured, with its integrands (None: 1).
    let mut measured: Vec<(u32, Vec<Option<u32>>)> = Vec::new();
    let mut at = |shape: u32, g: Option<u32>| -> (usize, usize) {
        let si = measured.iter().position(|m| m.0 == shape).unwrap_or_else(|| {
            measured.push((shape, Vec::new()));
            measured.len() - 1
        });
        let gs = &mut measured[si].1;
        let ii = gs.iter().position(|&h| h == g).unwrap_or_else(|| {
            gs.push(g);
            gs.len() - 1
        });
        (si, ii)
    };
    let n = parts.first().map_or(0, |p| p.2.len());
    // For each integrand (and, of a mean, its weight first), its terms.
    let mut terms: Vec<Vec<Term>> = vec![Vec::new(); n + usize::from(mean)];
    // Every part's own shape, whether or not a term reads it (an integrand of 0 over a half-space is no
    // 0 the bounds show): proven bounded, and reading no ranged param.
    let mut own: Vec<u32> = Vec::new();
    for (shape, w, hs) in parts {
        let p = rigidly(v, *shape, over, DEPTH)?;
        own.push(p.shape);
        let g = &mut v.doc.g;
        let one = g.constant(1.0);
        // Its weight in its own frame: a moving one reads no x, y or z.
        if w.is_some_and(|w| g.depends_on(w, |x| x == Var::T) || (!p.still && g.depends_on(w, xyz))) {
            return None;
        }
        let weighted = |g: &mut Graph, h: u32| w.map_or(h, |w| g.mul(w, h));
        if mean {
            terms[0].push((one, at(p.shape, *w)));
        }
        for (i, &h) in hs.iter().enumerate() {
            let list = &mut terms[i + usize::from(mean)];
            if p.still && !g.depends_on(h, |x| x == Var::T) {
                let wh = weighted(g, h);
                list.push((one, at(p.shape, Some(wh))));
                continue;
            }
            let [cx, cy, cz, c0] = affine_over(g, h, over)?;
            let cs = [cx, cy, cz];
            // c0 + c.t, times W.
            let mut a0 = c0.unwrap_or_else(|| g.constant(0.0));
            for (i, ci) in cs.iter().enumerate() {
                if let Some(ci) = *ci {
                    let m = mul(g, ci, p.t[i]);
                    a0 = add(g, a0, m);
                }
            }
            if !zero(g, a0) {
                list.push((a0, at(p.shape, *w)));
            }
            // (r^T c)_k, times the k-th of M.
            for k in 0..3 {
                let mut a = g.constant(0.0);
                for (i, ci) in cs.iter().enumerate() {
                    if let Some(ci) = *ci {
                        let m = mul(g, ci, p.r[i][k]);
                        a = add(g, a, m);
                    }
                }
                if !zero(g, a) {
                    let q = g.var([Var::X, Var::Y, Var::Z][k]);
                    let wq = weighted(g, q);
                    list.push((a, at(p.shape, Some(wq))));
                }
            }
        }
    }
    // What is measured once reads no ranged param: it would be measured at one value of it.
    let read: Vec<u32> = measured.iter().flat_map(|(s, gs)| std::iter::once(*s).chain(gs.iter().flatten().copied())).chain(own.iter().copied()).collect();
    if !v.ranged_reads(&read).is_empty() || read.iter().any(|&f| v.doc.g.depends_on(f, |x| x == Var::T)) {
        return None;
    }
    if !own.iter().all(|&s| v.whole(s)) {
        return None;
    }
    let mut regions = Vec::new();
    for (shape, gs) in &measured {
        if !v.whole(*shape) {
            return None;
        }
        let tapes: Vec<Option<Tape>> = gs.iter().map(|g| g.map(|g| Tape::new(&v.doc.g, g))).collect();
        regions.push((v.bounds(*shape).map(|(lo, hi, _)| (v.clipped(*shape), lo, hi)), tapes));
    }
    let alphas: Vec<Vec<(Tape, (usize, usize))>> = terms.iter().map(|ts| ts.iter().map(|&(a, ix)| (Tape::new(&v.doc.g, a), ix)).collect()).collect();
    let nodes: Vec<Vec<u32>> = terms.iter().map(|ts| ts.iter().map(|t| t.0).collect()).collect();
    for until in claims::measure_steps(res) {
        // Each shape's integrals: their estimates and proven bounds.
        let mut sums: Vec<Vec<(f64, I)>> = Vec::new();
        let mut coarsest: f64 = 0.0;
        for (region, tapes) in &regions {
            match region {
                None => sums.push(vec![(0.0, I::point(0.0)); tapes.len()]),
                Some((tape, lo, hi)) => {
                    let gs: Vec<Option<&Tape>> = tapes.iter().map(Option::as_ref).collect();
                    let m = query::measured(tape, &gs, false, &v.ctx, *lo, *hi, until);
                    coarsest = coarsest.max(m.res);
                    sums.push((0..tapes.len()).map(|k| (m.sums[k], I::new(m.lo[k], m.hi[k]))).collect());
                }
            }
        }
        if let Some((status, t, pieces)) = halved(&alphas, &sums, mean, c, over, span) {
            // The values at the moment given: a failure's, else the span's.
            let at_t = t.unwrap_or(moment);
            let ctx = Ctx { t: at_t, t_span: None, ..v.ctx.clone() };
            let value = |v: &ModelView, k: usize| -> f64 { nodes[k].iter().zip(&terms[k]).map(|(&a, &(_, (si, ii)))| v.doc.scalar(a, &ctx) * sums[si][ii].0).sum() };
            let values: Vec<f64> = (0..terms.len()).map(|k| value(v, k)).collect();
            let at_piece = |x: I| {
                if c.holds_on(x) {
                    Status::Holds
                } else if c.fails_on(x) {
                    Status::Fails
                } else {
                    Status::Unknown
                }
            };
            let each = |k: usize| match t {
                None => Status::Holds,
                Some(t) => at_piece(piece_value(&alphas, &sums, mean, k, &Ctx { t, t_span: Some((t, t)), ..over.clone() })),
            };
            let evidence = match mean {
                false => Evidence::Integral { integrals: values.iter().enumerate().map(|(k, &x)| (x, each(k))).collect(), res: coarsest, sampled: None },
                true => {
                    Evidence::Mean { weight: values[0], means: values[1..].iter().enumerate().map(|(k, &x)| (x / values[0], each(k + 1))).collect(), res: coarsest, sampled: None }
                }
            };
            let m = Motion { t, steps: pieces, ..Motion::none(if status == Status::Holds { S::Proven } else { S::Disproven }, span) };
            return Some(Check { motion: Some(Stepped { ranged: true, ..Stepped::rigid(m) }), ..claims::checked(status, evidence) });
        }
    }
    None
}

/// Integrand k's value over the scope `ctx` (a piece of the span): its
/// terms' functions of t bounded there, times their integrals' bounds,
/// summed; of a mean, over the weight's (integrand 0).
fn piece_value(alphas: &[Vec<(Tape, (usize, usize))>], sums: &[Vec<(f64, I)>], mean: bool, k: usize, ctx: &Ctx) -> I {
    let total = |k: usize| {
        alphas[k].iter().fold(I::point(0.0), |s, (a, (si, ii))| {
            let x = eval_interval_scope(a, ctx, [I::point(0.0); 3], None);
            if !(x.lo.is_finite() && x.hi.is_finite()) {
                return I::ENTIRE;
            }
            add_o(s, interval::bin(BinOp::Mul, x, sums[*si][*ii].1).0)
        })
    };
    let x = total(k);
    match mean && k > 0 {
        false => x,
        true => {
            let w = total(0);
            if !(w.lo > 0.0) {
                return I::ENTIRE;
            }
            interval::bin(BinOp::Div, x, w).0
        }
    }
}

/// The claim over `span` from its integrals' bounds (`piece_value`),
/// asked over pieces of it, each undecided one halved, at most PIECES in
/// all: Holds where every piece holds, Fails over the first piece every
/// moment of which fails (its middle the moment), with how many were
/// asked; None where pieces are left undecided.
fn halved(alphas: &[Vec<(Tape, (usize, usize))>], sums: &[Vec<(f64, I)>], mean: bool, c: Comparison, over: &Ctx, span: (f64, f64)) -> Option<(Status, Option<f64>, usize)> {
    let first = usize::from(mean);
    let mut queue = std::collections::VecDeque::from([span]);
    let mut asked = 0;
    while let Some(p) = queue.pop_front() {
        asked += 1;
        let ctx = Ctx { t: span_mid(p), t_span: Some(p), ..over.clone() };
        let values: Vec<I> = (first..alphas.len()).map(|k| piece_value(alphas, sums, mean, k, &ctx)).collect();
        if values.iter().any(|&x| c.fails_on(x)) {
            return Some((Status::Fails, Some(span_mid(p)), asked));
        }
        if values.iter().all(|&x| c.holds_on(x)) {
            continue;
        }
        let mid = span_mid(p);
        if asked + queue.len() + 2 > PIECES || !(p.0 < mid && mid < p.1) {
            return None;
        }
        queue.extend([(p.0, mid), (mid, p.1)]);
    }
    Some((Status::Holds, None, asked))
}
