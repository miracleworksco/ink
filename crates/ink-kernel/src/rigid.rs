//! Rigid placements with numbers, and mass properties: a body's bounds and
//! moments computed once in its own frame, then moved (instancing). And
//! placements made of fields (of t and params): how far one is from rigid
//! over a scope (`turn_bounds`), and how a frame's reading of one differs
//! from the world's (`frame_delta`).

use crate::query::sym_eigen;
use ink_core::V3;
use ink_core::eval::Ctx;
use ink_core::graph::{BinOp, Graph, Node, NodeId, UnOp, Var, bin_f64, un_f64};
use ink_core::interval::{self, I, add_o, mul_o, sub_o};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

/// A rigid placement: p_world = r p_local + t.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rigid {
    pub r: [[f64; 3]; 3],
    pub t: V3,
}

impl Rigid {
    pub const IDENTITY: Rigid = Rigid { r: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], t: [0.0; 3] };

    /// A placement from its twelve numbers: r row by row, then t.
    pub fn from_entries(v: [f64; 12]) -> Rigid {
        Rigid { r: [[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]], t: [v[9], v[10], v[11]] }
    }

    /// A point moved (rounded to nearest: for meshes and estimates, not bounds).
    pub fn apply(&self, p: V3) -> V3 {
        std::array::from_fn(|i| self.r[i][0] * p[0] + self.r[i][1] * p[1] + self.r[i][2] * p[2] + self.t[i])
    }

    /// World box of a local box, rounded outward, so a proven outer bound
    /// stays one under the forward map: each coordinate's least and greatest
    /// t + r x over the box, term by term (in real numbers exactly the box of its corners),
    /// every product and sum rounded down or up. Exact steps stay exact, so
    /// an identity placement keeps the box. An instance's world field uses
    /// R^T as its inverse map: callers must first widen for R^-T - R.
    pub fn aabb(&self, lo: V3, hi: V3) -> (V3, V3) {
        let (mut a, mut b) = (self.t, self.t);
        for i in 0..3 {
            for k in 0..3 {
                let (p, q) = (mul_dir(self.r[i][k], lo[k]), mul_dir(self.r[i][k], hi[k]));
                a[i] = add_dir(a[i], p.0.min(q.0)).0;
                b[i] = add_dir(b[i], p.1.max(q.1)).1;
            }
        }
        (a, b)
    }

    /// Moments [m, mx, my, mz, xx, yy, zz, xy, yz, zx] (see `query::moments`)
    /// of a body moved by this placement: exact in real numbers (the centre
    /// of mass moves, the second moments about it rotate). Estimates, rounded
    /// to nearest.
    pub fn moments(&self, m: &[f64; 10]) -> [f64; 10] {
        let mass = m[0];
        if mass == 0.0 {
            return *m;
        }
        let c = [m[1] / mass, m[2] / mass, m[3] / mass];
        let s = [
            [m[4] - mass * c[0] * c[0], m[7] - mass * c[0] * c[1], m[9] - mass * c[2] * c[0]],
            [m[7] - mass * c[0] * c[1], m[5] - mass * c[1] * c[1], m[8] - mass * c[1] * c[2]],
            [m[9] - mass * c[2] * c[0], m[8] - mass * c[1] * c[2], m[6] - mass * c[2] * c[2]],
        ];
        let r = self.r;
        let sw: [[f64; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| (0..3).map(|l| r[i][k] * s[k][l] * r[j][l]).sum::<f64>()).sum()));
        let cw = self.apply(c);
        let raw = |i: usize, j: usize| sw[i][j] + mass * cw[i] * cw[j];
        [mass, mass * cw[0], mass * cw[1], mass * cw[2], raw(0, 0), raw(1, 1), raw(2, 2), raw(0, 1), raw(1, 2), raw(2, 0)]
    }
}

/// a b rounded down and up. Exact by the fma residual away from underflow
/// and overflow (elsewhere a step each way); 0 times anything is 0, as in
/// interval products.
fn mul_dir(a: f64, b: f64) -> (f64, f64) {
    if a == 0.0 || b == 0.0 {
        return (0.0, 0.0);
    }
    let p = a * b;
    let e = if p.is_finite() && p.abs() >= f64::MIN_POSITIVE * 2f64.powi(54) { a.mul_add(b, -p) } else { f64::NAN };
    (if e >= 0.0 { p } else { p.next_down() }, if e <= 0.0 { p } else { p.next_up() })
}

/// a + b rounded down and up: the two-sum error is exact unless the sum
/// overflows. An infinite term gives an infinite sum (and inf - inf no
/// bound).
fn add_dir(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    if s.is_nan() {
        return (f64::NEG_INFINITY, f64::INFINITY);
    }
    if a.is_infinite() || b.is_infinite() {
        return (s, s);
    }
    let bb = s - a;
    let e = if s.is_finite() { (a - (s - bb)) + (b - bb) } else { f64::NAN };
    (if e >= 0.0 { s } else { s.next_down() }, if e <= 0.0 { s } else { s.next_up() })
}

/// Mass properties from moments: the inertia tensor is about the centre of
/// mass; principal moments ascend, each with its axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProps {
    pub mass: f64,
    pub center: V3,
    pub inertia: [[f64; 3]; 3],
    pub principal_moments: [f64; 3],
    pub principal_axes: [V3; 3],
}

/// Mass properties of moments [m, mx, my, mz, xx, yy, zz, xy, yz, zx]; None
/// when there is no mass.
pub fn mass_props(m: &[f64; 10]) -> Option<MassProps> {
    let mass = m[0];
    if mass <= 0.0 {
        return None;
    }
    let c = [m[1] / mass, m[2] / mass, m[3] / mass];
    // Second moments about the centre of mass.
    let (sxx, syy, szz) = (m[4] - mass * c[0] * c[0], m[5] - mass * c[1] * c[1], m[6] - mass * c[2] * c[2]);
    let (sxy, syz, szx) = (m[7] - mass * c[0] * c[1], m[8] - mass * c[1] * c[2], m[9] - mass * c[2] * c[0]);
    let inertia = [[syy + szz, -sxy, -szx], [-sxy, sxx + szz, -syz], [-szx, -syz, sxx + syy]];
    let (vals, vecs) = sym_eigen(inertia);
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| vals[i].total_cmp(&vals[j]));
    Some(MassProps { mass, center: c, inertia, principal_moments: order.map(|i| vals[i]), principal_axes: order.map(|i| [vecs[0][i], vecs[1][i], vecs[2][i]]) })
}

/// How far a placement's turn is from a rotation over a scope: see
/// `turn_bounds`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TurnBounds {
    /// Holds det r.
    pub det: I,
    /// At least the largest row sum of |r^T r - I| (its infinity norm).
    pub skew: f64,
}

impl TurnBounds {
    /// Bound on |(R^T R)^-1 - I| from the Neumann series. None when
    /// invertibility is not established. Exact zero stays zero.
    pub fn slack(&self) -> Option<f64> {
        if !(self.skew < 1.0) {
            return None;
        }
        Some(if self.skew == 0.0 { 0.0 } else { interval::bin(BinOp::Div, I::point(self.skew), sub_o(I::point(1.0), I::point(self.skew))).0.hi })
    }
}

/// A placement p -> r p + t proven an isometry at every point of a scope
/// (every moment and param value of the Ctx it was proven over): r a proper
/// rotation exactly (`turn_bounds`: det exactly 1, r^T r exactly I) and
/// every entry of r and t finite there (outward bounds). So the field it
/// places, f(r^T (p - t)), holds f's solid moved by p = r q + t, a rigid
/// bijection of space, and a box of f's moves as that map moves it. A
/// matrix proven a rotation alone is not one: a translation infinite at a
/// moment (a move by ln((t - 0.5)^2)) places nothing there. Made only by
/// `over`, so a weaker check cannot stand in for it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Isometry {
    r: [[NodeId; 3]; 3],
    t: [NodeId; 3],
}

impl Isometry {
    /// `(r, t)` proven an isometry over `scope`, else None.
    pub fn over(g: &Graph, r: &[[NodeId; 3]; 3], t: &[NodeId; 3], scope: &Ctx) -> Option<Isometry> {
        if !finite_over(g, r.iter().flatten().chain(t), scope) {
            return None;
        }
        let b = turn_bounds(g, r, scope)?;
        (b.skew == 0.0 && b.det.lo == 1.0 && b.det.hi == 1.0).then_some(Isometry { r: *r, t: *t })
    }

    pub fn r(&self) -> &[[NodeId; 3]; 3] {
        &self.r
    }

    pub fn t(&self) -> &[NodeId; 3] {
        &self.t
    }
}

/// Whether each of `nodes` is finite at every point of `scope`: its
/// bounds there (over t's span and the ranged params), rounded outward.
pub fn finite_over<'a>(g: &Graph, nodes: impl IntoIterator<Item = &'a NodeId>, scope: &Ctx) -> bool {
    nodes.into_iter().all(|&n| {
        let x = ink_core::eval::eval_interval_scope(&ink_core::tape::Tape::new(g, n), scope, [I::point(0.0); 3], None);
        x.lo.is_finite() && x.hi.is_finite()
    })
}

/// Bounds, over every moment and param value of `ctx`'s scope, on how far
/// a placement's turn `r` (p_world = r p + t), whose entries are fields of
/// t and params, is from a rotation: its determinant, and r^T r - I. None
/// where an entry has too many terms to read or too high a power, a fixed
/// value is not finite, an atom other than a sine or cosine may be
/// infinite over the scope, or a bound is not finite.
///
/// Interval arithmetic cannot see that a turn through a ranged angle is a
/// rotation: over a quarter turn, cos^2 + sin^2 bounds [0, 2]. So each
/// entry is read as a polynomial in the cosines and sines of its angles
/// (any other node is a variable of its own: an atom), with coefficients
/// bounded by intervals rounded outward and exact where the floats are, and
/// sin^2 is written 1 - cos^2 until no sine is squared. Only what reads t
/// or a ranged param is read so: an expression of constants and fixed
/// params means its f64 value, so it is that float, a coefficient, before
/// any algebra ((p + 1) - p is 0 for p = 1e16, as the evaluator has it, not
/// 1). What is read means its value in real numbers; with every atom finite
/// the polynomials' algebra is that arithmetic (inf - inf is 0, not what
/// cancelling says), and sin^2 + cos^2 = 1 wherever the angle is finite
/// (sin and cos of an infinite angle are 0), so it is used for an angle
/// only where its bounds over the scope are finite: the polynomial's value
/// is the entry's at every point of the scope. Its bounds are each term's,
/// summed, rounded outward. A turn whose entries are exact (about an axis
/// of the frame, by a ranged angle: cos, sin, 0 and 1) then has det exactly
/// 1 and r^T r exactly I; a turn through a fixed angle, or composed with
/// one (a turn of 30 degrees is a matrix of floats, which no rotation is
/// exactly), is within a few ulps of them, as far as its floats are.
pub fn turn_bounds(g: &Graph, r: &[[NodeId; 3]; 3], ctx: &Ctx) -> Option<TurnBounds> {
    let mut read = Reading::new(g, ctx);
    let mut e = Vec::with_capacity(9);
    for row in r {
        for &n in row {
            e.push(read.poly(n)?);
        }
    }
    let at = |i: usize, j: usize| &e[3 * i + j];
    let minor = |a: (usize, usize), b: (usize, usize), c: (usize, usize), d: (usize, usize)| at(a.0, a.1).mul(at(b.0, b.1))?.add(&at(c.0, c.1).mul(at(d.0, d.1))?, -1.0);
    let det = at(0, 0)
        .mul(&minor((1, 1), (2, 2), (1, 2), (2, 1))?)?
        .add(&at(0, 1).mul(&minor((1, 0), (2, 2), (1, 2), (2, 0))?)?, -1.0)?
        .add(&at(0, 2).mul(&minor((1, 0), (2, 1), (1, 1), (2, 0))?)?, 1.0)?;
    // r^T r - I, on and above the diagonal.
    let mut gram = Vec::with_capacity(6);
    for i in 0..3 {
        for j in i..3 {
            let mut s = Poly::constant(I::point(if i == j { -1.0 } else { 0.0 }));
            for k in 0..3 {
                s = s.add(&at(k, i).mul(at(k, j))?, 1.0)?;
            }
            gram.push(s);
        }
    }
    let finite = read.finite()?;
    let mut atoms = HashMap::new();
    let mut value = |p: Poly| bound(g, ctx, &finite, &mut atoms, p);
    let det = value(det)?;
    let gram: Vec<I> = gram.into_iter().map(&mut value).collect::<Option<_>>()?;
    // Row i of the symmetric matrix from its upper triangle (0 1 2 / 3 4 / 5).
    let index = |i: usize, j: usize| {
        let (i, j) = (i.min(j), i.max(j));
        [0, 3, 5][i] + j - i
    };
    let skew = (0..3).map(|i| (0..3).fold(0.0, |s, j| interval::add_up(s, interval::mag(gram[index(i, j)])))).fold(0.0, f64::max);
    skew.is_finite().then_some(TurnBounds { det, skew })
}

/// A polynomial's bounds over `ctx`'s scope, reduced (sin^2 = 1 - cos^2 for
/// the `finite` angles), each atom bounded once (`atoms`), rounded outward.
/// None where they are not finite.
fn bound(g: &Graph, ctx: &Ctx, finite: &HashSet<NodeId>, atoms: &mut HashMap<Atom, I>, p: Poly) -> Option<I> {
    let p = p.reduced(&|n| finite.contains(&n))?;
    let mut of = |a: Atom| -> I { *atoms.entry(a).or_insert_with(|| atom_bounds(a, &ink_core::tape::Tape::new(g, a.node()), ctx)) };
    p.value(&mut of)
}

/// An atom's bounds over ctx's scope, its node's tape given.
fn atom_bounds(a: Atom, node: &ink_core::tape::Tape, ctx: &Ctx) -> I {
    let x = ink_core::eval::eval_interval(node, ctx, [I::ENTIRE; 3], None);
    match a {
        Atom::Cos(_) => interval::un(UnOp::Cos, x),
        Atom::Sin(_) => interval::un(UnOp::Sin, x),
        Atom::Other(_) => x,
    }
}

/// Fields of t and params (placement entries) rewritten as the polynomials
/// `turn_bounds` reads them as, in the sines and cosines of their angles
/// (sin^2 written 1 - cos^2), so a turn composed with its own inverse is
/// read as what is left: R_a^T R_a, a turn's sines and cosines summed in
/// products that interval bounds over a piece of a motion blur by the turn,
/// is exactly I, and R_a^T (t_a + R_a p - t_a) is p. Equal to the entries in
/// real numbers wherever ctx's scope reaches (as turn_bounds' reading is).
/// None where an entry is not so read, or a coefficient is not one float.
pub fn simplified(g: &mut Graph, nodes: &[NodeId], ctx: &Ctx) -> Option<Vec<NodeId>> {
    let polys: Vec<Poly> = {
        let mut read = Reading::new(g, ctx);
        let ps: Vec<Poly> = nodes.iter().map(|&n| read.poly(n)).collect::<Option<_>>()?;
        let finite = read.finite()?;
        ps.into_iter().map(|p| p.reduced(&|n| finite.contains(&n))).collect::<Option<_>>()?
    };
    polys
        .iter()
        .map(|p| {
            let mut sum: Option<NodeId> = None;
            for (m, c) in &p.0 {
                if c.lo != c.hi {
                    return None;
                }
                let mut term = g.constant(c.lo);
                for &(a, k) in m {
                    let base = match a {
                        Atom::Cos(n) => g.un(UnOp::Cos, n),
                        Atom::Sin(n) => g.un(UnOp::Sin, n),
                        Atom::Other(n) => n,
                    };
                    for _ in 0..k {
                        term = g.mul(term, base);
                    }
                }
                sum = Some(match sum {
                    None => term,
                    Some(s) => g.add(s, term),
                });
            }
            Some(sum.unwrap_or_else(|| g.constant(0.0)))
        })
        .collect()
}

/// How a field read in a frame differs from the world's read back in it,
/// over every value of `ctx`'s scope: a frame placed at S, the field placed
/// in it by Q and in the world by W (each p_world = R p + T: R row major,
/// then T), read at ρ_W(S c) and ρ_Q(c) (ρ_X(p) = R^T (p - T)), differ by
/// Δ(c) = M c + d, M = R_W^T R_S - R_Q^T and d = R_W^T (T_S - T_W) + R_Q^T
/// T_Q. Each entry is read as a polynomial (see `turn_bounds`) and bounded
/// over the scope: a frame and a world that move together (parts riding a
/// carriage, S and W reading its slide's terms alike) cancel them, where
/// interval bounds of each entry over a motion would differ by its whole
/// travel. M (row major) and d. None: not read, or not bounded.
pub fn frame_delta(g: &Graph, ctx: &Ctx, s: &[NodeId; 12], q: &[NodeId; 12], w: &[NodeId; 12]) -> Option<([I; 9], [I; 3])> {
    let mut read = Reading::new(g, ctx);
    let mut polys = |ns: &[NodeId; 12]| -> Option<Vec<Poly>> { ns.iter().map(|&n| read.poly(n)).collect() };
    let (s, q, w) = (polys(s)?, polys(q)?, polys(w)?);
    let r = |p: &[Poly], i: usize, j: usize| p[3 * i + j].clone();
    let t = |p: &[Poly], i: usize| p[9 + i].clone();
    let mut m = Vec::with_capacity(9);
    for i in 0..3 {
        for j in 0..3 {
            let mut x = Poly::default().add(&r(&q, j, i), -1.0)?;
            for k in 0..3 {
                x = x.add(&r(&w, k, i).mul(&r(&s, k, j))?, 1.0)?;
            }
            m.push(x);
        }
    }
    let mut d = Vec::with_capacity(3);
    for i in 0..3 {
        let mut x = Poly::default();
        for k in 0..3 {
            x = x.add(&r(&w, k, i).mul(&t(&s, k).add(&t(&w, k), -1.0)?)?, 1.0)?;
            x = x.add(&r(&q, k, i).mul(&t(&q, k))?, 1.0)?;
        }
        d.push(x);
    }
    let finite = read.finite()?;
    let mut atoms = HashMap::new();
    let mut value = |p: Poly| bound(g, ctx, &finite, &mut atoms, p);
    let m: Vec<I> = m.into_iter().map(&mut value).collect::<Option<_>>()?;
    let d: Vec<I> = d.into_iter().map(&mut value).collect::<Option<_>>()?;
    Some((m.try_into().ok()?, d.try_into().ok()?))
}

/// The value of a field of t and params (an entry of a placement) where it
/// is one number over all of `ctx`'s scope, exactly: where its polynomial
/// (see `turn_bounds`) is a constant with no rounding. A turn about z
/// through a ranged angle, as Rodrigues' formula builds it, is (1 - cos) +
/// cos on the axis, which interval arithmetic over a whole turn bounds by
/// [-1, 3]: this is 1. An expression of constants and fixed params is its
/// f64 value.
pub fn constant_over(g: &Graph, n: NodeId, ctx: &Ctx) -> Option<f64> {
    let mut read = Reading::new(g, ctx);
    let p = read.poly(n)?;
    let finite = read.finite()?;
    let p = p.reduced(&|n| finite.contains(&n))?;
    match p.0.len() {
        0 => Some(0.0),
        1 => p.0.get(&Vec::new()).filter(|c| c.lo == c.hi).map(|c| c.lo),
        _ => None,
    }
}

/// A field of t and params that is one number wherever its angles are
/// finite, to within a few of its floats' roundings: (c, d), its
/// polynomial's (see `turn_bounds`) constant term, and a bound on how far
/// it is from c: the constant's own width and each other term's
/// coefficient, where each such term is a product of sines and cosines
/// (each at most 1). A body placed on a turning frame, seen from that
/// frame (a pulley on a motor's shaft frame), reads its turn's cos^2 +
/// sin^2 and the floats of the frames between: c, and d a few ulps. None
/// where another term reads any other atom (t, a param ranged in ctx), or
/// it cannot be read.
pub fn nearly_constant(g: &Graph, n: NodeId, ctx: &Ctx) -> Option<(f64, f64)> {
    let mut read = Reading::new(g, ctx);
    let p = read.poly(n)?;
    let finite = read.finite()?;
    let p = p.reduced(&|n| finite.contains(&n))?;
    let (mut c, mut d) = (0.0, 0.0);
    for (m, k) in &p.0 {
        if m.is_empty() {
            c = k.mid();
            d = interval::add_up(d, interval::mag(sub_o(*k, I::point(c))));
        } else if m.iter().all(|&(a, _)| matches!(a, Atom::Cos(x) | Atom::Sin(x) if finite.contains(&x))) {
            d = interval::add_up(d, interval::mag(*k));
        } else {
            return None;
        }
    }
    d.is_finite().then_some((c, d))
}

/// How a turn about z winds with its angles, exactly: `(x, y)` its first
/// column's entries (cos and sin of the turn), read as polynomials (see
/// `turn_bounds`) whose only atoms are the cosines and sines of angles
/// a_i, and integers m_i with x + i y = w prod_i (cos a_i + i sin a_i)^m_i
/// for a constant w, identically: the turn is sum_i m_i a_i plus a
/// constant. A placement composing a turn by a twice (`(rot-z q) (rot-z
/// q)`) winds 2 with it. None where that is not shown: another atom, more
/// than three angles or a winding past 64, or a remainder that is not
/// exactly constant (a float coefficient left where terms should cancel).
pub fn winding(g: &Graph, x: NodeId, y: NodeId, ctx: &Ctx) -> Option<Vec<(NodeId, i32)>> {
    winding_near(g, x, y, ctx).filter(|w| w.1 == 0.0).map(|w| w.0)
}

/// `winding` where the remainder need not be exactly constant: the winding
/// whose remainder's other terms are least, with the sum of their
/// coefficients' magnitudes (each atom a sine or cosine, at most 1: how far
/// the remainder is from its constant), 0 where it is exact. A turn about
/// an axis tilted e from the frame's winds as the frame's would, its
/// remainder of order e^2. None as `winding`, or where two windings are
/// equally near.
pub fn winding_near(g: &Graph, x: NodeId, y: NodeId, ctx: &Ctx) -> Option<(Vec<(NodeId, i32)>, f64)> {
    let mut read = Reading::new(g, ctx);
    let (px, py) = (read.poly(x)?, read.poly(y)?);
    let finite = read.finite()?;
    let (px, py) = (px.reduced(&|n| finite.contains(&n))?, py.reduced(&|n| finite.contains(&n))?);
    // Each angle's degree, the most its winding can be.
    let mut degree: BTreeMap<NodeId, u32> = BTreeMap::new();
    for m in px.0.keys().chain(py.0.keys()) {
        let mut d: BTreeMap<NodeId, u32> = BTreeMap::new();
        for &(a, k) in m {
            match a {
                Atom::Cos(n) | Atom::Sin(n) if finite.contains(&n) => *d.entry(n).or_default() += k as u32,
                _ => return None,
            }
        }
        for (n, k) in d {
            let e = degree.entry(n).or_default();
            *e = (*e).max(k);
        }
    }
    let angles: Vec<(NodeId, u32)> = degree.into_iter().collect();
    if angles.len() > 3 || angles.iter().any(|&(_, d)| d == 0 || d > 64) {
        return None;
    }
    let one = Poly::constant(I::point(1.0));
    // (re, im) products, exactly.
    let cmul = |(a, b): &(Poly, Poly), (c, d): &(Poly, Poly)| -> Option<(Poly, Poly)> { Some((a.mul(c)?.add(&b.mul(d)?, -1.0)?, a.mul(d)?.add(&b.mul(c)?, 1.0)?)) };
    // How far a remainder is from its constant: its other terms'
    // coefficients, rounded up.
    let off = |p: &Poly| p.0.iter().filter(|(m, _)| !m.is_empty()).fold(0.0, |s, (_, c)| interval::add_up(s, interval::mag(*c)));
    let mut found: Option<(Vec<(NodeId, i32)>, f64)> = None;
    let mut tied = false;
    // Each angle wound +d or -d: the sign choice that leaves the first
    // column nearest a constant times the product.
    for signs in 0..1u32 << angles.len() {
        let mut w = (one.clone(), Poly::default());
        for (i, &(a, d)) in angles.iter().enumerate() {
            let sin = if signs >> i & 1 == 1 { Poly::default().add(&Poly::atom(Atom::Sin(a)), -1.0)? } else { Poly::atom(Atom::Sin(a)) };
            let base = (Poly::atom(Atom::Cos(a)), sin);
            for _ in 0..d {
                w = cmul(&w, &base)?;
            }
        }
        // (x + i y) times w's conjugate.
        let (wr, wi) = w;
        let re = px.mul(&wr)?.add(&py.mul(&wi)?, 1.0)?.reduced(&|n| finite.contains(&n))?;
        let im = py.mul(&wr)?.add(&px.mul(&wi)?, -1.0)?.reduced(&|n| finite.contains(&n))?;
        let w = re.0.contains_key(&Vec::new()) || im.0.contains_key(&Vec::new());
        let e = interval::add_up(off(&re), off(&im));
        if !(w && e.is_finite()) || found.as_ref().is_some_and(|f| f.1 < e) {
            continue;
        }
        tied = found.as_ref().is_some_and(|f| f.1 == e);
        found = Some((angles.iter().enumerate().map(|(i, &(a, d))| (a, if signs >> i & 1 == 1 { -(d as i32) } else { d as i32 })).collect(), e));
    }
    found.filter(|_| !tied)
}

/// A variable of a placement's polynomials: the cosine or sine of an angle
/// (the node they are of), or any other node.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Atom {
    Cos(NodeId),
    Sin(NodeId),
    Other(NodeId),
}

impl Atom {
    fn node(self) -> NodeId {
        match self {
            Atom::Cos(n) | Atom::Sin(n) | Atom::Other(n) => n,
        }
    }
}

/// A monomial: its atoms with their powers, sorted by atom.
type Mono = Vec<(Atom, u8)>;

/// A polynomial: monomials with their coefficients' bounds (none exactly 0).
#[derive(Clone, Debug, Default)]
struct Poly(BTreeMap<Mono, I>);

/// Terms a polynomial may have: past them a placement is not read. A
/// power past u8 is not read either (`mono_mul`).
const POLY_TERMS: usize = 4096;

impl Poly {
    fn constant(c: I) -> Poly {
        let mut p = Poly::default();
        p.term(Vec::new(), c);
        p
    }

    fn atom(a: Atom) -> Poly {
        Poly(BTreeMap::from([(vec![(a, 1)], I::point(1.0))]))
    }

    /// Adds c m.
    fn term(&mut self, m: Mono, c: I) {
        let s = match self.0.get(&m) {
            Some(&d) => add_i(d, c),
            None => c,
        };
        if s.lo == 0.0 && s.hi == 0.0 {
            self.0.remove(&m);
        } else {
            self.0.insert(m, s);
        }
    }

    /// self + sign o (sign 1 or -1); None past POLY_TERMS.
    fn add(&self, o: &Poly, sign: f64) -> Option<Poly> {
        let mut p = self.clone();
        for (m, &c) in &o.0 {
            p.term(m.clone(), if sign < 0.0 { neg(c) } else { c });
        }
        (p.0.len() <= POLY_TERMS).then_some(p)
    }

    /// None past POLY_TERMS, or where a power passes u8.
    fn mul(&self, o: &Poly) -> Option<Poly> {
        let mut p = Poly::default();
        for (m, &c) in &self.0 {
            for (n, &d) in &o.0 {
                p.term(mono_mul(m, n)?, mul_i(c, d));
                if p.0.len() > POLY_TERMS {
                    return None;
                }
            }
        }
        Some(p)
    }

    /// Its bounds, each atom's given by `of`: each term's, summed, rounded
    /// outward. None where they are not finite.
    fn value(&self, of: &mut dyn FnMut(Atom) -> I) -> Option<I> {
        let v = self.0.iter().fold(I::point(0.0), |s, (m, c)| {
            let term = m.iter().fold(*c, |t, &(a, k)| mul_o(t, ipow(of(a), k)));
            add_o(s, term)
        });
        (v.lo.is_finite() && v.hi.is_finite()).then_some(v)
    }

    /// sin^2 written 1 - cos^2, for each angle `finite` holds of, until no
    /// such sine is squared. None past POLY_TERMS, or where a power passes
    /// u8.
    fn reduced(self, finite: &dyn Fn(NodeId) -> bool) -> Option<Poly> {
        let mut out = Poly::default();
        let mut work: Vec<(Mono, I)> = self.0.into_iter().collect();
        while let Some((m, c)) = work.pop() {
            match m.iter().position(|&(a, k)| matches!(a, Atom::Sin(n) if k >= 2 && finite(n))) {
                None => out.term(m, c),
                Some(i) => {
                    let Atom::Sin(n) = m[i].0 else { unreachable!() };
                    let mut less = m;
                    less[i].1 -= 2;
                    if less[i].1 == 0 {
                        less.remove(i);
                    }
                    let cos2 = mono_mul(&less, &vec![(Atom::Cos(n), 2)])?;
                    work.push((less, c));
                    work.push((cos2, neg(c)));
                }
            }
            if out.0.len() + work.len() > POLY_TERMS {
                return None;
            }
        }
        Some(out)
    }
}

fn neg(c: I) -> I {
    I::new(-c.hi, -c.lo)
}

/// A sum and a product of coefficients' bounds, rounded outward, and exact
/// where the floats are (1 - cos + cos is 1, not 1 and a cosine times
/// 1e-324, which a cancelling turn's polynomial must not keep).
fn add_i(a: I, b: I) -> I {
    I::new(add_dir(a.lo, b.lo).0, add_dir(a.hi, b.hi).1)
}

fn mul_i(a: I, b: I) -> I {
    let p = [mul_dir(a.lo, b.lo), mul_dir(a.lo, b.hi), mul_dir(a.hi, b.lo), mul_dir(a.hi, b.hi)];
    I::new(p.iter().map(|x| x.0).fold(f64::INFINITY, f64::min), p.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max))
}

/// m n; None where a power passes u8 (it is not read: a wrapped power
/// would be another polynomial).
fn mono_mul(m: &Mono, n: &Mono) -> Option<Mono> {
    let mut out: Mono = Vec::with_capacity(m.len() + n.len());
    let (mut i, mut j) = (0, 0);
    while i < m.len() || j < n.len() {
        match (m.get(i), n.get(j)) {
            (Some(&(a, k)), Some(&(b, l))) if a == b => {
                out.push((a, k.checked_add(l)?));
                (i, j) = (i + 1, j + 1);
            }
            (Some(&x), Some(&(b, _))) if x.0 < b => {
                out.push(x);
                i += 1;
            }
            (Some(_), Some(&y)) | (None, Some(&y)) => {
                out.push(y);
                j += 1;
            }
            (Some(&x), None) => {
                out.push(x);
                i += 1;
            }
            (None, None) => unreachable!(),
        }
    }
    Some(out)
}

/// Bounds of x^k, rounded outward (squares, so an even power is not negative).
fn ipow(x: I, k: u8) -> I {
    match k {
        0 => I::point(1.0),
        1 => x,
        _ if k.is_multiple_of(2) => ipow(interval::un(UnOp::Square, x), k / 2),
        _ => mul_o(x, ipow(x, k - 1)),
    }
}

/// Fields of t and params read as polynomials in their atoms, with every
/// atom the reading met (some may cancel out of the polynomials), over
/// `ctx`'s scope.
struct Reading<'a> {
    g: &'a Graph,
    ctx: &'a Ctx,
    memo: HashMap<NodeId, Poly>,
    /// Each node's f64 value where it is fixed (`Reading::fixed`).
    fixed: HashMap<NodeId, Option<f64>>,
    atoms: BTreeSet<Atom>,
}

impl<'a> Reading<'a> {
    fn new(g: &'a Graph, ctx: &'a Ctx) -> Self {
        Reading { g, ctx, memo: HashMap::new(), fixed: HashMap::new(), atoms: BTreeSet::new() }
    }

    /// The f64 value of an expression of constants and params none of which
    /// is ranged in `ctx` (it means that value), as the evaluators compute
    /// it (graph::un_f64 and bin_f64, and eval's fixed values); None for
    /// any other.
    fn fixed(&mut self, n: NodeId) -> Option<f64> {
        if let Some(&v) = self.fixed.get(&n) {
            return v;
        }
        let v = match self.g.nodes[n as usize] {
            Node::Const(c) => Some(c),
            Node::Var(Var::Param(k)) if self.ctx.p_spans.get(k).is_none() => self.ctx.params.get(k as usize).copied(),
            Node::Un(op, a) => self.fixed(a).map(|a| un_f64(op, a)),
            Node::Bin(op, a, b) => match (self.fixed(a), self.fixed(b)) {
                (Some(a), Some(b)) => Some(bin_f64(op, a, b)),
                _ => None,
            },
            Node::Select(c, a, b) => match (self.fixed(c), self.fixed(a), self.fixed(b)) {
                (Some(c), Some(a), Some(b)) => Some(if c < 0.0 { a } else { b }),
                _ => None,
            },
            Node::Tag(_, a) => self.fixed(a),
            _ => None,
        };
        self.fixed.insert(n, v);
        v
    }

    /// A field as a polynomial: an expression of constants and fixed
    /// params is its f64 value; else sums, differences, products,
    /// negations and squares of those and of atoms (a name is its
    /// expression), and quotients by a fixed value other than 0 (products
    /// by its reciprocal, bounded: a belt's shift written (/ (* -5 x) 24)),
    /// and anything else is an atom. None where a fixed value is not
    /// finite, past POLY_TERMS, or where a power passes u8.
    fn poly(&mut self, n: NodeId) -> Option<Poly> {
        if let Some(p) = self.memo.get(&n) {
            return Some(p.clone());
        }
        let atom = |r: &mut Self, a: Atom| {
            r.atoms.insert(a);
            Poly::atom(a)
        };
        if let Some(v) = self.fixed(n) {
            return v.is_finite().then(|| Poly::constant(I::point(v)));
        }
        let p = match self.g.nodes[n as usize] {
            Node::Un(UnOp::Neg, a) => Poly::default().add(&self.poly(a)?, -1.0)?,
            Node::Un(UnOp::Square, a) => {
                let a = self.poly(a)?;
                a.mul(&a)?
            }
            Node::Un(UnOp::Cos, a) => atom(self, Atom::Cos(a)),
            Node::Un(UnOp::Sin, a) => atom(self, Atom::Sin(a)),
            Node::Bin(BinOp::Add, a, b) => self.poly(a)?.add(&self.poly(b)?, 1.0)?,
            Node::Bin(BinOp::Sub, a, b) => self.poly(a)?.add(&self.poly(b)?, -1.0)?,
            Node::Bin(BinOp::Mul, a, b) => self.poly(a)?.mul(&self.poly(b)?)?,
            Node::Bin(BinOp::Div, a, b) if self.fixed(b).is_some_and(|c| c != 0.0) => {
                let c = self.fixed(b).expect("fixed");
                self.poly(a)?.mul(&Poly::constant(interval::bin(BinOp::Div, I::point(1.0), I::point(c)).0))?
            }
            Node::Tag(_, a) => self.poly(a)?,
            _ => atom(self, Atom::Other(n)),
        };
        self.memo.insert(n, p.clone());
        Some(p)
    }

    /// The angles whose bounds over the scope are finite, so sin^2 + cos^2
    /// = 1 at each of its points. None where another atom may be infinite:
    /// the polynomials' algebra ((x + 1) - x = 1) is the fields' value only
    /// where every atom is finite (inf - inf is 0).
    fn finite(&self) -> Option<HashSet<NodeId>> {
        let bounds = |n: NodeId| ink_core::eval::eval_interval(&ink_core::tape::Tape::new(self.g, n), self.ctx, [I::ENTIRE; 3], None);
        let finite = |b: I| b.lo.is_finite() && b.hi.is_finite();
        let mut angles = HashSet::new();
        for &a in &self.atoms {
            match a {
                Atom::Other(n) if !finite(bounds(n)) => return None,
                Atom::Cos(n) | Atom::Sin(n) if finite(bounds(n)) => {
                    angles.insert(n);
                }
                _ => {}
            }
        }
        Some(angles)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Raw moments of a uniform box of mass `m`, centre `c`, half sizes `h`.
    fn box_moments(m: f64, c: V3, h: V3) -> [f64; 10] {
        let s = |i: usize| m * (h[i] * h[i] / 3.0 + c[i] * c[i]);
        [m, m * c[0], m * c[1], m * c[2], s(0), s(1), s(2), m * c[0] * c[1], m * c[1] * c[2], m * c[2] * c[0]]
    }

    fn close(a: &[f64], b: &[f64], tol: f64) -> bool {
        a.iter().zip(b).all(|(x, y)| (x - y).abs() <= tol * (1.0 + x.abs().max(y.abs())))
    }

    /// A turn seen from itself: R(a)^T R(a) and R(a)^T (t + R(a) p - t),
    /// summed products of a turning angle's sine and cosine, which interval
    /// bounds over a piece of the turn blur by its angle (the frame of a
    /// moving part: ink-lang's motion::moving_frame), are rewritten as I and
    /// p exactly; R(a)^T R(a + b) as R(b). Each rewritten entry equals its
    /// entry at sampled moments.
    #[test]
    fn a_turn_composed_with_its_inverse_is_rewritten_as_what_is_left() {
        let mut g = Graph::new();
        let t = g.var(Var::T);
        let k = g.constant(6.0);
        let a = g.mul(k, t);
        let k2 = g.constant(0.3);
        let s2 = g.un(UnOp::Sin, t);
        let b = g.mul(k2, s2);
        let ab = g.add(a, b);
        let turn = |g: &mut Graph, x: NodeId| -> [[NodeId; 2]; 2] {
            let (c, s) = (g.un(UnOp::Cos, x), g.un(UnOp::Sin, x));
            let ns = g.neg(s);
            [[c, ns], [s, c]]
        };
        let (ra, rab) = (turn(&mut g, a), turn(&mut g, ab));
        let (tx, ty) = (g.constant(15.0), g.constant(-4.0));
        let dot = |g: &mut Graph, u: [NodeId; 2], w: [NodeId; 2]| {
            let (p, q) = (g.mul(u[0], w[0]), g.mul(u[1], w[1]));
            g.add(p, q)
        };
        let col = |m: &[[NodeId; 2]; 2], j: usize| [m[0][j], m[1][j]];
        // R(a)^T R(a), R(a)^T R(a + b), and R(a)^T (R(a) (15, -4)).
        let mut entries = Vec::new();
        for m in [&ra, &rab] {
            for i in 0..2 {
                for j in 0..2 {
                    entries.push(dot(&mut g, col(&ra, i), col(m, j)));
                }
            }
        }
        let moved = [dot(&mut g, ra[0], [tx, ty]), dot(&mut g, ra[1], [tx, ty])];
        for i in 0..2 {
            entries.push(dot(&mut g, col(&ra, i), moved));
        }
        let ctx = Ctx { t_span: Some((0.0, 1.0)), ..Ctx::new(0.5, vec![]) };
        let out = simplified(&mut g, &entries, &ctx).expect("read as polynomials");
        let konst = |n: NodeId| g.as_const(n);
        assert_eq!([0, 1, 2, 3].map(|i| konst(out[i])), [Some(1.0), Some(0.0), Some(0.0), Some(1.0)]);
        assert_eq!([konst(out[8]), konst(out[9])], [Some(15.0), Some(-4.0)]);
        // R(b): reads no turn of a.
        let (cb, sb) = (g.un(UnOp::Cos, b), g.un(UnOp::Sin, b));
        let nsb = g.neg(sb);
        for (i, want) in [cb, nsb, sb, cb].into_iter().enumerate() {
            let tape = ink_core::tape::Tape::new(&g, out[4 + i]);
            let w = ink_core::tape::Tape::new(&g, want);
            for q in 0..=20 {
                let c = Ctx::new(q as f64 / 20.0, vec![]);
                let (x, y) = (ink_core::eval::eval_point(&tape, &c, [0.0; 3]), ink_core::eval::eval_point(&w, &c, [0.0; 3]));
                assert!((x - y).abs() < 1e-12, "entry {i} at t {}: {x} vs {y}", c.t);
            }
        }
        for (e, o) in entries.iter().zip(&out) {
            let (te, to) = (ink_core::tape::Tape::new(&g, *e), ink_core::tape::Tape::new(&g, *o));
            for q in 0..=20 {
                let c = Ctx::new(q as f64 / 20.0, vec![]);
                let (x, y) = (ink_core::eval::eval_point(&te, &c, [0.0; 3]), ink_core::eval::eval_point(&to, &c, [0.0; 3]));
                assert!((x - y).abs() < 1e-12, "at t {}: {x} vs {y}", c.t);
            }
        }
    }

    #[test]
    fn moved_moments_are_the_moved_body_s() {
        // A 6 x 4 x 2 box off its frame's origin, turned 90 degrees about z
        // and moved: the same box, its half sizes swapped in x and y.
        let local = box_moments(2.5, [3.0, 1.0, 0.5], [3.0, 2.0, 1.0]);
        let rg = Rigid { r: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], t: [1.0, 2.0, 3.0] };
        let world = box_moments(2.5, [0.0, 5.0, 3.5], [2.0, 3.0, 1.0]);
        let moved = rg.moments(&local);
        assert!(close(&moved, &world, 1e-14), "{moved:?} vs {world:?}");
        // Inertia about the centre of mass: a box's, diagonal; principal
        // moments ascending, each with its axis.
        let p = mass_props(&moved).unwrap();
        let (a, b, c) = (2.0, 3.0, 1.0);
        let want = [2.5 * (b * b + c * c) / 3.0, 2.5 * (a * a + c * c) / 3.0, 2.5 * (a * a + b * b) / 3.0];
        assert!(close(&p.center, &[0.0, 5.0, 3.5], 1e-14));
        assert!(close(&[p.inertia[0][0], p.inertia[1][1], p.inertia[2][2]], &want, 1e-12), "{:?}", p.inertia);
        let mut sorted = want;
        sorted.sort_by(f64::total_cmp);
        assert!(close(&p.principal_moments, &sorted, 1e-12), "{:?}", p.principal_moments);
        // The smallest moment is about the longest side (y here).
        assert!(p.principal_axes[0][1].abs() > 1.0 - 1e-12, "{:?}", p.principal_axes);
        // No mass: nothing to move, no properties.
        assert_eq!(rg.moments(&[0.0; 10]), [0.0; 10]);
        assert!(mass_props(&[0.0; 10]).is_none());
    }

    #[test]
    fn moved_box_rounds_outward() {
        // 0.1 + 0.7 rounds to nearest below the real sum of those doubles, so
        // moving the corner by t = 0.1 to nearest would cut the box short.
        let rg = Rigid { t: [0.1, 0.0, 0.0], ..Rigid::IDENTITY };
        let (lo, hi) = rg.aabb([0.7; 3], [0.7; 3]);
        assert!(lo[0] <= 0.1 + 0.7 && hi[0] > 0.1 + 0.7, "{lo:?} {hi:?}");
        // A turned box: the box of its turned corners, to the rounding.
        let (c, s) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
        let rg = Rigid { r: [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], t: [5.0, -2.0, 7.0] };
        let (blo, bhi) = ([-1.0, 0.0, -2.0], [4.0, 3.0, 2.0]);
        let (lo, hi) = rg.aabb(blo, bhi);
        let corners: Vec<V3> = (0..8).map(|k| rg.apply([[blo[0], bhi[0]][k & 1], [blo[1], bhi[1]][(k >> 1) & 1], [blo[2], bhi[2]][(k >> 2) & 1]])).collect();
        for i in 0..3 {
            let (a, b) = corners.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), w| (a.min(w[i]), b.max(w[i])));
            assert!(lo[i] <= a && b <= hi[i] && a - lo[i] < 1e-14 && hi[i] - b < 1e-14, "{i}: {lo:?} {hi:?} vs {a} {b}");
        }
        // Exact steps stay exact: the identity keeps a box, a quarter turn
        // and a whole move swap and shift it; an unbounded side stays
        // unbounded, and a zero entry ignores it.
        let b = ([f64::NEG_INFINITY, 0.25, -3.0], [1.0, 1.0, 1.5]);
        assert_eq!(Rigid::IDENTITY.aabb(b.0, b.1), b);
        let quarter = Rigid { r: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], t: [1.0, 2.0, 3.0] };
        assert_eq!(quarter.aabb(b.0, b.1), ([0.0, f64::NEG_INFINITY, 0.0], [0.75, 3.0, 4.5]));
    }

    /// A turn's bounds over a scope. Interval arithmetic
    /// alone bounds cos^2 + sin^2 over a whole turn by [0, 2]: read as
    /// polynomials with sin^2 = 1 - cos^2, a turn about z through a ranged
    /// angle has det exactly 1 and r^T r exactly I. A turn of floats is no
    /// rotation: one whose determinant is known exactly (cos 30 rounded to
    /// 26 bits, so a^2 + b^2 is a float, not 1), composed with a ranged turn
    /// about x, has its det enclosed (a bound of 1 would miss it), and
    /// r^T r - I at least as far from 0 as it is at the angle 0. A stretch keeps
    /// its own det.
    #[test]
    fn turn_bounds_hold_rotations_and_their_floats() {
        use ink_core::eval::ParamSpans;
        let mut g = Graph::new();
        let (zero, one) = (g.constant(0.0), g.constant(1.0));
        // A turn about axis i by the param, as the language's `rot-x`..
        // `rot-z` make it: cos and sin, and (1 - cos) + cos on the axis.
        let turn = |g: &mut Graph, axis: usize| {
            let a = g.var(Var::Param(0));
            let k = g.constant(std::f64::consts::PI / 180.0);
            let r = g.mul(a, k);
            let (c, s) = (g.un(UnOp::Cos, r), g.un(UnOp::Sin, r));
            let (omc, ns) = (g.sub(one, c), g.neg(s));
            let on = g.add(omc, c);
            let (i, j) = ((axis + 1) % 3, (axis + 2) % 3);
            let mut m = [[zero; 3]; 3];
            m[axis][axis] = on;
            (m[i][i], m[i][j], m[j][i], m[j][j]) = (c, ns, s, c);
            m
        };
        let product = |g: &mut Graph, a: [[NodeId; 3]; 3], b: [[NodeId; 3]; 3]| -> [[NodeId; 3]; 3] {
            std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    let mut e = g.constant(0.0);
                    for k in 0..3 {
                        let p = g.mul(a[i][k], b[k][j]);
                        e = g.add(e, p);
                    }
                    e
                })
            })
        };
        let whole = Ctx { p_spans: ParamSpans::new([(0, (0.0, 360.0))]), ..Ctx::new(0.0, vec![0.0]) };
        let rz = turn(&mut g, 2);
        assert_eq!(turn_bounds(&g, &rz, &whole), Some(TurnBounds { det: I::point(1.0), skew: 0.0 }));
        let cs = g.mul(rz[0][0], rz[0][0]);
        let ss = g.mul(rz[1][0], rz[1][0]);
        let plain = g.add(cs, ss);
        let plain = ink_core::eval::eval_interval(&ink_core::tape::Tape::new(&g, plain), &whole, [I::ENTIRE; 3], None);
        assert!(plain.lo <= 0.0 && plain.hi >= 2.0, "{plain:?}");
        // a has 26 significant bits, so a^2 + 1/4 is a float: exact.
        let a = f64::from_bits(30f64.to_radians().cos().to_bits() & !((1u64 << 27) - 1));
        let exact = a * a + 0.25;
        assert!(exact != 1.0 && a.mul_add(a, 0.25) == exact, "{exact}");
        let (ca, cb, nb) = (g.constant(a), g.constant(0.5), g.constant(-0.5));
        let floats = [[ca, nb, zero], [cb, ca, zero], [zero, zero, one]];
        let rx = turn(&mut g, 0);
        let r = product(&mut g, floats, rx);
        let b = turn_bounds(&g, &r, &whole).unwrap();
        assert!(b.det.contains(exact) && b.det.hi - b.det.lo < 1e-14, "{b:?}");
        assert!(b.skew >= (exact - 1.0).abs() && b.skew < 1e-7, "{b:?}");
        // The turn's axis entry, (1 - cos) + cos, is 1 over the scope. A
        // sum that cancels is its constant only where what cancels is
        // finite: e^p overflows for p near 1000, and (inf + 1) - inf is 0.
        assert_eq!(constant_over(&g, rz[2][2], &whole), Some(1.0));
        let p = g.var(Var::Param(0));
        let e = g.un(UnOp::Exp, p);
        let e1 = g.add(e, one);
        let cancels = g.sub(e1, e);
        let up_to = |hi: f64| Ctx { p_spans: ParamSpans::new([(0, (0.0, hi))]), ..Ctx::new(0.0, vec![0.0]) };
        assert_eq!((constant_over(&g, cancels, &up_to(1.0)), constant_over(&g, cancels, &up_to(1000.0))), (Some(1.0), None));
        // Stretched 2 along x: det 2, and r^T r - I is 3 on its first row.
        let (two, z) = (g.constant(2.0), zero);
        let stretch = [[two, z, z], [z, one, z], [z, z, one]];
        let r = product(&mut g, stretch, rz);
        let b = turn_bounds(&g, &r, &whole).unwrap();
        assert!(b.det.contains(2.0) && b.det.hi - b.det.lo < 1e-14 && b.skew >= 3.0, "{b:?}");
        // An expression of constants and fixed params means its f64 value, so
        // a placement reads it as that float before any algebra: (p + 1) - p
        // is 0 for p = 1e16, as the evaluator has it, and 1 only where p is
        // ranged (its value in real numbers). A turn through a fixed angle of
        // 1 radian is the floats cos 1 and sin 1, whose squares sum to 1 +
        // 4.8e-17, not a rotation; through a ranged angle (at the same value)
        // it is one.
        {
            use ink_core::eval::eval_interval;
            use ink_core::tape::Tape;
            let mut g = Graph::new();
            let (zero, one) = (g.constant(0.0), g.constant(1.0));
            let p = g.var(Var::Param(0));
            let p1 = g.add(p, one);
            let d = g.sub(p1, p);
            let at = |v: f64, ranged: bool| Ctx { p_spans: ParamSpans::new(ranged.then_some((0, (v, v)))), ..Ctx::new(0.0, vec![v]) };
            let bound = |n: NodeId, ctx: &Ctx| eval_interval(&Tape::new(&g, n), ctx, [I::ENTIRE; 3], None);
            assert_eq!((bound(d, &at(1e16, false)), constant_over(&g, d, &at(1e16, false))), (I::point(0.0), Some(0.0)));
            assert!(bound(d, &at(1e16, true)).contains(1.0) && constant_over(&g, d, &at(1e16, true)) == Some(1.0));
            let (c, s) = (g.un(UnOp::Cos, p), g.un(UnOp::Sin, p));
            let ns = g.neg(s);
            let r = [[c, ns, zero], [s, c, zero], [zero, zero, one]];
            // cos^2 1 + sin^2 1 - 1 of the floats: each square exactly as two
            // floats (fma), the leading halves' sum exactly as two (two-sum), and
            // its leading half minus 1 exact (Sterbenz).
            let (cf, sf) = (1f64.cos(), 1f64.sin());
            let (a1, b1) = (cf * cf, sf * sf);
            let (a2, b2) = (cf.mul_add(cf, -a1), sf.mul_add(sf, -b1));
            let s1 = a1 + b1;
            let e1 = (a1 - (s1 - (s1 - a1))) + (b1 - (s1 - a1));
            let defect = (s1 - 1.0) + e1 + a2 + b2;
            assert!(defect > 4.8e-17 && defect < 4.9e-17, "{defect}");
            let b = turn_bounds(&g, &r, &at(1.0, false)).unwrap();
            assert!(b.det.lo <= 1.0 && b.det.hi > 1.0 && b.det.hi - b.det.lo < 1e-15 && b.skew >= defect && b.skew < 1e-15, "{b:?}");
            assert_eq!(turn_bounds(&g, &r, &at(1.0, true)), Some(TurnBounds { det: I::point(1.0), skew: 0.0 }));
        }
    }

    #[test]
    fn quotients_by_constants_are_read() {
        use ink_core::eval::ParamSpans;
        let mut g = Graph::new();
        let a = g.var(Var::Param(0));
        let (m5, c24, k) = (g.constant(-5.0), g.constant(24.0), g.constant(5.0 / 24.0));
        let n = g.mul(m5, a);
        let q = g.bin(BinOp::Div, n, c24);
        let p = g.mul(k, a);
        let sum = g.add(q, p);
        let (four, three) = (g.constant(4.0), g.constant(3.0));
        let quarter = g.bin(BinOp::Div, a, four);
        let x = g.mul(three, quarter);
        let c = g.constant(0.75);
        let y = g.mul(c, a);
        let exact = g.sub(x, y);
        let ctx = Ctx { p_spans: ParamSpans::new([(0, (0.0, 1e5))]), ..Ctx::new(0.0, vec![0.0]) };
        let mut read = Reading::new(&g, &ctx);
        let p = read.poly(sum).unwrap();
        let b = bound(&g, &ctx, &read.finite().unwrap(), &mut HashMap::new(), p).unwrap();
        assert!(interval::mag(b) < 1e-9, "{b:?}");
        assert_eq!(constant_over(&g, exact, &ctx), Some(0.0));
    }

    /// A frame slid along x by a quintic stroke over t in [0, 1] (its three
    /// terms, up to 600 mm each), a part placed in it at a float turn (30
    /// degrees about z) and (3.7, 0.1, 0): the world composes them, the
    /// question reads the part at that placement. Over the whole stroke how
    /// they differ is a few ulps, the terms cancelling (the slide's bounds
    /// over it are -1500 to 1600 mm); a world sliding 1 mm further is 1 off.
    #[test]
    fn frame_delta_cancels_a_motion_the_frame_and_world_share() {
        let mut g = Graph::new();
        let (zero, one) = (g.constant(0.0), g.constant(1.0));
        let t = g.var(Var::T);
        let mut x = zero;
        for (k, n) in [(1000.0, 3), (-1500.0, 4), (600.0, 5)] {
            let mut p = g.constant(k);
            for _ in 0..n {
                p = g.mul(p, t);
            }
            x = g.add(x, p);
        }
        let id = [[one, zero, zero], [zero, one, zero], [zero, zero, one]];
        let (sn, cs) = 30f64.to_radians().sin_cos();
        let (c, sp, sm) = (g.constant(cs), g.constant(sn), g.constant(-sn));
        let r = [[c, sm, zero], [sp, c, zero], [zero, zero, one]];
        let q = [g.constant(3.7), g.constant(0.1), zero];
        let entries = |r: [[NodeId; 3]; 3], t: [NodeId; 3]| -> [NodeId; 12] { std::array::from_fn(|n| if n < 9 { r[n / 3][n % 3] } else { t[n - 9] }) };
        // The world: the frame's slide composed with the part's placement.
        let world = |g: &mut Graph, x: NodeId| -> [NodeId; 12] {
            let w: [NodeId; 3] = std::array::from_fn(|i| if i == 0 { g.add(q[0], x) } else { q[i] });
            entries(r, w)
        };
        let ctx = Ctx { t_span: Some((0.0, 1.0)), ..Ctx::new(0.5, vec![]) };
        let (s, qe) = (entries(id, [x, zero, zero]), entries(r, q));
        let w = world(&mut g, x);
        let (m, d) = frame_delta(&g, &ctx, &s, &qe, &w).unwrap();
        assert!(m.iter().chain(&d).all(|e| interval::mag(*e) < 1e-12), "{m:?} {d:?}");
        let further = g.add(x, one);
        let w = world(&mut g, further);
        let (_, d) = frame_delta(&g, &ctx, &s, &qe, &w).unwrap();
        let length = d.iter().map(|e| e.mid() * e.mid()).sum::<f64>().sqrt();
        assert!((length - 1.0).abs() < 1e-9, "{d:?}");
    }

    /// A power past a monomial's u8 is not read: wrapped, the determinant
    /// of (cos a I)^86, cos^258 a, would read cos^2 a, and 128 such turns'
    /// r^T r, cos^256 a I, would read I.
    /// cos^255 a, the highest power, is read.
    #[test]
    fn turn_bounds_refuse_powers_past_their_type() {
        use ink_core::eval::ParamSpans;
        let mut g = Graph::new();
        let zero = g.constant(0.0);
        let a = g.var(Var::Param(0));
        let c = g.un(UnOp::Cos, a);
        let ctx = Ctx { p_spans: ParamSpans::new([(0, (0.5, 0.6))]), ..Ctx::new(0.0, vec![0.55]) };
        let scaled = |g: &mut Graph, k: usize| {
            let e = (1..k).fold(c, |e, _| g.mul(e, c));
            [[e, zero, zero], [zero, e, zero], [zero, zero, e]]
        };
        let r = scaled(&mut g, 85);
        let b = turn_bounds(&g, &r, &ctx).unwrap();
        assert!(b.det.lo <= 0.6f64.cos().powi(255) && 0.5f64.cos().powi(255) <= b.det.hi && b.det.hi < 1e-14, "{b:?}");
        for k in [86, 128] {
            let r = scaled(&mut g, k);
            assert_eq!(turn_bounds(&g, &r, &ctx), None, "{k}");
        }
    }

    /// An isometry over a scope is a rotation exactly and a translation
    /// finite at every point of it: a turn about z by 360 t moved by t is
    /// one over t in [0, 1]; moved by ln((t - 0.5)^2) (infinite at 0.5) it
    /// is not, nor is the turn scaled by 2.
    #[test]
    fn an_isometry_is_a_rotation_with_a_finite_translation_over_the_scope() {
        let mut g = Graph::new();
        let t = g.var(Var::T);
        let k = g.constant(std::f64::consts::TAU);
        let a = g.mul(k, t);
        let (c, s, z, o) = (g.un(UnOp::Cos, a), g.un(UnOp::Sin, a), g.constant(0.0), g.constant(1.0));
        let ns = g.neg(s);
        let r = [[c, ns, z], [s, c, z], [z, z, o]];
        let half = g.constant(0.5);
        let d = g.sub(t, half);
        let sq = g.un(UnOp::Square, d);
        let far = g.un(UnOp::Ln, sq);
        let two = g.constant(2.0);
        let scaled = r.map(|row| row.map(|n| if n == z || n == o { n } else { g.mul(two, n) }));
        let scope = Ctx { t_span: Some((0.0, 1.0)), ..Ctx::new(0.5, vec![]) };
        assert!(Isometry::over(&g, &r, &[t, z, z], &scope).is_some());
        assert!(Isometry::over(&g, &r, &[far, z, z], &scope).is_none());
        assert!(Isometry::over(&g, &scaled, &[t, z, z], &scope).is_none());
    }
}
