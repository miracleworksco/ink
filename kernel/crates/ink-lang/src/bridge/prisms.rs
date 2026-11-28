//! A bridge made of upright prisms, read along the directions asked: material
//! that is a union of open convex polygons in x and y (walls at any angle,
//! a part turned about z), each between two heights, is proved at its
//! actual underside, as boxes.rs proves one of boxes square to x and y.
//! Interval cells straddling a roof's underside never settle the rows
//! themselves; here the underside is known exactly, so the rule std's rows
//! read is shown on each roof's footprint in the plane instead.
//!
//! It is decided at load, in interval arithmetic on the shape's own numbers
//! (a shape or check that reads a param or t is left to the rows). Every
//! comparison is outward, so a doubt leaves the ceiling to the rows.
use super::*;
use ink_core::interval::{I, add_o, div_o, mul_o, sub_o};
use std::collections::BTreeMap;

/// n · (x, y) < c, over the real numbers these enclose. `key` is the face it
/// came from (its leaf and level nodes and side): pieces cut by one face
/// share it, so a point of one is on its inner side in every other. Bit i
/// of `along` is set where direction i runs exactly along it: n · d is 0
/// in the real numbers (`Exact`), not only within its enclosure.
#[derive(Clone, Copy)]
struct Side {
    n: [I; 2],
    c: I,
    key: (NodeId, NodeId, bool),
    along: u64,
}

/// A number as a polynomial in the turn coefficients it reads (Node::Within:
/// one identity, one number), each term's float coefficient exact: so a
/// polynomial that cancels to none of them is 0 in the real numbers,
/// whatever the coefficients are. Pure-constant arithmetic is kept only
/// where its float is exact, so it agrees whether read in floats or reals.
type Exact = BTreeMap<Vec<NodeId>, f64>;

fn exact_sum(a: &Exact, b: &Exact, sign: f64) -> Option<Exact> {
    let mut out = a.clone();
    for (m, &v) in b {
        let u = out.get(m).copied().unwrap_or(0.0);
        let v = sign * v;
        // Two-sum: the float sum is the real one where its error is 0.
        let t = u + v;
        let back = t - u;
        if !t.is_finite() || (u - (t - back)) + (v - back) != 0.0 {
            return None;
        }
        if t == 0.0 {
            out.remove(m)
        } else {
            out.insert(m.clone(), t)
        };
    }
    Some(out)
}

fn exact_product(a: &Exact, b: &Exact) -> Option<Exact> {
    let mut out = Exact::new();
    for (m, &u) in a {
        for (n, &v) in b {
            let p = u * v;
            // Normal and with no rounding error, the float product is the real one.
            if !(p.is_finite() && p.abs() >= 2f64.powi(-900) && u.mul_add(v, -p) == 0.0) {
                return None;
            }
            let mut mono: Vec<NodeId> = m.iter().chain(n).copied().collect();
            mono.sort_unstable();
            out = exact_sum(&out, &Exact::from([(mono, p)]), 1.0)?;
        }
    }
    (out.len() <= 64).then_some(out)
}

/// An open prism: z between its two bounds, inside every side.
#[derive(Clone, Default)]
struct Piece {
    z: [Option<I>; 2],
    sides: Vec<Side>,
}

fn intersect(a: &Piece, b: &Piece) -> Piece {
    let pick = |x: Option<I>, y: Option<I>, low: bool| match (x, y) {
        (Some(x), Some(y)) if low => Some(I::new(x.lo.max(y.lo), x.hi.max(y.hi))),
        (Some(x), Some(y)) => Some(I::new(x.lo.min(y.lo), x.hi.min(y.hi))),
        (x, y) => x.or(y),
    };
    Piece { z: [pick(a.z[0], b.z[0], true), pick(a.z[1], b.z[1], false)], sides: a.sides.iter().chain(&b.sides).copied().collect() }
}

fn neg(a: I) -> I {
    I::new(-a.hi, -a.lo)
}

fn zero(a: I) -> bool {
    a == I::point(0.0)
}

fn dot(n: [I; 2], x: I, y: I) -> I {
    add_o(mul_o(n[0], x), mul_o(n[1], y))
}

/// The value of a node that reads no variable (not space, t or a param),
/// enclosed as the proofs evaluate it. None where it reads one, or is not
/// finite.
pub(super) fn constant(g: &Graph, n: NodeId) -> Option<I> {
    if g.depends_on(n, |_| true) {
        return None;
    }
    let v = eval::eval_interval(&Tape::new(g, n), &Ctx::new(0.0, vec![]), [I::point(0.0); 3], None);
    (v.lo.is_finite() && v.hi.is_finite()).then_some(v)
}

/// Arithmetic of its operands alone (a number within two constants among
/// them: a turn's cosine): no variable, sample or frame.
fn plain(node: Node) -> bool {
    matches!(node, Node::Const(_) | Node::Within(..) | Node::Un(..) | Node::Bin(..) | Node::Select(..) | Node::Tag(..) | Node::Real(..))
}

struct Reader<'a> {
    g: &'a mut Graph,
    free: HashMap<NodeId, bool>,
    values: HashMap<NodeId, Option<I>>,
    affine: HashMap<NodeId, Option<[I; 4]>>,
    exact: HashMap<NodeId, Option<Exact>>,
    plane: HashMap<NodeId, Option<[Exact; 2]>>,
    half_planes: HashMap<(NodeId, NodeId, bool), Option<[Exact; 4]>>,
    /// Each direction asked, exactly, up to a positive or negative factor.
    dirs: Vec<Option<[Exact; 2]>>,
    budget: usize,
}

impl Reader<'_> {
    /// Whether n reads no variable: not space, t or a param.
    fn free(&mut self, n: NodeId) -> bool {
        if let Some(&f) = self.free.get(&n) {
            return f;
        }
        let node = self.g.nodes[n as usize];
        let (kids, k) = node.children();
        let f = plain(node) && kids[..k].iter().all(|&c| self.free(c));
        self.free.insert(n, f);
        f
    }

    /// `constant`, where `free` says it reads no variable.
    fn value(&mut self, n: NodeId) -> Option<I> {
        if let Some(&v) = self.values.get(&n) {
            return v;
        }
        let v = if self.free(n) { constant(self.g, n) } else { None };
        self.values.insert(n, v);
        v
    }

    /// The number n exactly (`Exact`), where it is one: constants, turn
    /// coefficients and their sums, differences and products.
    fn exact(&mut self, n: NodeId) -> Option<Exact> {
        if let Some(e) = self.exact.get(&n) {
            return e.clone();
        }
        let node = self.g.nodes[n as usize];
        let e = match node {
            Node::Const(0.0) => Some(Exact::new()),
            Node::Const(v) if v.is_finite() => Some(Exact::from([(vec![], v)])),
            Node::Within(..) => Some(Exact::from([(vec![n], 1.0)])),
            Node::Tag(_, a) | Node::Real(a) => self.exact(a),
            Node::Un(UnOp::Neg, a) => self.exact(a).map(|a| a.into_iter().map(|(m, v)| (m, -v)).collect()),
            Node::Bin(op @ (BinOp::Add | BinOp::Sub | BinOp::Mul), a, b) => {
                let (a, b) = (self.exact(a), self.exact(b));
                a.zip(b).and_then(|(a, b)| match op {
                    BinOp::Mul => exact_product(&a, &b),
                    _ => exact_sum(&a, &b, if op == BinOp::Add { 1.0 } else { -1.0 }),
                })
            }
            _ => None,
        };
        self.exact.insert(n, e.clone());
        e
    }

    /// n's coefficients of x and y exactly, where n is affine in x, y and z
    /// with each of them exact (`Exact`).
    fn plane(&mut self, n: NodeId) -> Option<[Exact; 2]> {
        if let Some(p) = self.plane.get(&n) {
            return p.clone();
        }
        let none = || [Exact::new(), Exact::new()];
        let one = || Exact::from([(vec![], 1.0)]);
        // A number reads no x or y.
        let p = if self.free(n) {
            Some(none())
        } else {
            let node = self.g.nodes[n as usize];
            match node {
                Node::Var(Var::X) => Some([one(), Exact::new()]),
                Node::Var(Var::Y) => Some([Exact::new(), one()]),
                Node::Var(Var::Z) => Some(none()),
                Node::Tag(_, a) | Node::Real(a) | Node::Lipschitz(a, ..) => self.plane(a),
                Node::Un(UnOp::Neg, a) => self.plane(a).map(|p| p.map(|c| c.into_iter().map(|(m, v)| (m, -v)).collect())),
                Node::Bin(op @ (BinOp::Add | BinOp::Sub), a, b) => {
                    let (a, b) = (self.plane(a), self.plane(b));
                    let sign = if op == BinOp::Add { 1.0 } else { -1.0 };
                    a.zip(b).and_then(|(a, b)| Some([exact_sum(&a[0], &b[0], sign)?, exact_sum(&a[1], &b[1], sign)?]))
                }
                Node::Bin(BinOp::Mul, a, b) if self.free(a) || self.free(b) => {
                    let (k, v) = if self.free(a) { (a, b) } else { (b, a) };
                    let k = self.exact(k);
                    k.zip(self.plane(v)).and_then(|(k, p)| Some([exact_product(&k, &p[0])?, exact_product(&k, &p[1])?]))
                }
                _ => None,
            }
        };
        self.plane.insert(n, p.clone());
        p
    }

    /// n's constant term exactly, where n is affine in x, y and z with it
    /// exact (`Exact`).
    fn offset(&mut self, n: NodeId) -> Option<Exact> {
        if self.free(n) {
            return self.exact(n);
        }
        let node = self.g.nodes[n as usize];
        match node {
            Node::Var(Var::X | Var::Y | Var::Z) => Some(Exact::new()),
            Node::Tag(_, a) | Node::Real(a) | Node::Lipschitz(a, ..) => self.offset(a),
            Node::Un(UnOp::Neg, a) => self.offset(a).map(|a| a.into_iter().map(|(m, v)| (m, -v)).collect()),
            Node::Bin(op @ (BinOp::Add | BinOp::Sub), a, b) => {
                let (a, b) = (self.offset(a), self.offset(b));
                a.zip(b).and_then(|(a, b)| exact_sum(&a, &b, if op == BinOp::Add { 1.0 } else { -1.0 }))
            }
            Node::Bin(BinOp::Mul, a, b) if self.free(a) || self.free(b) => {
                let (k, v) = if self.free(a) { (a, b) } else { (b, a) };
                let k = self.exact(k);
                k.zip(self.offset(v)).and_then(|(k, v)| exact_product(&k, &v))
            }
            _ => None,
        }
    }

    /// A side's half-plane exactly (`Exact`): n · (x, y) < l - o, as its
    /// normal n, level l and offset o (kept apart: their difference need
    /// not be a float).
    fn half_plane(&mut self, s: &Side) -> Option<[Exact; 4]> {
        if let Some(h) = self.half_planes.get(&s.key) {
            return h.clone();
        }
        let (id, level, less) = s.key;
        let flip = |e: Exact| e.into_iter().map(|(m, v)| (m, -v)).collect::<Exact>();
        let h = self.plane(id).zip(self.exact(level)).zip(self.offset(id)).map(|(([a, b], l), o)| match less {
            true => [a, b, l, o],
            false => [flip(a), flip(b), flip(l), flip(o)],
        });
        self.half_planes.insert(s.key, h.clone());
        h
    }

    /// Whether every point inside side `r` is inside `s`: the same face, or
    /// one with exactly the same normal and a bound at least r's.
    fn implies(&mut self, r: &Side, s: &Side) -> bool {
        if r.key == s.key {
            return true;
        }
        let (Some([a, b, l, o]), Some([c, d, m, p])) = (self.half_plane(r), self.half_plane(s)) else { return false };
        // Equal bounds: (l - m) - (o - p) is 0.
        let equal = || Some(exact_sum(&exact_sum(&l, &m, -1.0)?, &exact_sum(&o, &p, -1.0)?, -1.0)?.is_empty());
        a == c && b == d && (r.c.hi <= s.c.lo || equal() == Some(true))
    }

    /// The directions asked (bit i for direction i) that a face of id runs
    /// exactly along: its normal's dot with each is 0 in the real numbers.
    fn along(&mut self, id: NodeId) -> u64 {
        let Some([a, b]) = self.plane(id) else { return 0 };
        let dirs = self.dirs.clone();
        dirs.iter().enumerate().filter(|(i, _)| *i < 64).fold(0, |bits, (i, d)| {
            let flat = d.as_ref().and_then(|[dx, dy]| exact_sum(&exact_product(&a, dx)?, &exact_product(&b, dy)?, 1.0));
            if flat.is_some_and(|f| f.is_empty()) { bits | 1 << i } else { bits }
        })
    }

    /// A direction exactly, up to a factor: (p / q, r / q) is (p, r) over a
    /// q that is not 0, as std's unit directions are.
    fn direction(&mut self, [x, y]: [NodeId; 2]) -> Option<[Exact; 2]> {
        let (nx, ny) = (self.g.nodes[x as usize], self.g.nodes[y as usize]);
        let (x, y) = match (nx, ny) {
            (Node::Bin(BinOp::Div, p, q), Node::Bin(BinOp::Div, r, s)) if q == s && self.value(q).is_some_and(|q| q.lo > 0.0 || q.hi < 0.0) => (p, r),
            _ => (x, y),
        };
        Some([self.exact(x)?, self.exact(y)?])
    }

    /// n as a x + b y + c z + e, each coefficient enclosed.
    fn affine(&mut self, n: NodeId) -> Option<[I; 4]> {
        if let Some(&a) = self.affine.get(&n) {
            return a;
        }
        let o = I::point(0.0);
        let unit = I::point(1.0);
        let a = if self.free(n) {
            self.value(n).map(|v| [o, o, o, v])
        } else {
            let node = self.g.nodes[n as usize];
            match node {
                Node::Var(Var::X) => Some([unit, o, o, o]),
                Node::Var(Var::Y) => Some([o, unit, o, o]),
                Node::Var(Var::Z) => Some([o, o, unit, o]),
                Node::Tag(_, a) | Node::Real(a) | Node::Lipschitz(a, ..) => self.affine(a),
                Node::Un(UnOp::Neg, a) => self.affine(a).map(|a| a.map(neg)),
                Node::Bin(op @ (BinOp::Add | BinOp::Sub), a, b) => {
                    let (a, b) = (self.affine(a), self.affine(b));
                    a.zip(b).map(|(a, b)| std::array::from_fn(|i| if op == BinOp::Add { add_o(a[i], b[i]) } else { sub_o(a[i], b[i]) }))
                }
                Node::Bin(BinOp::Mul, a, b) if self.free(a) || self.free(b) => {
                    let (k, v) = if self.free(a) { (a, b) } else { (b, a) };
                    let k = self.value(k);
                    k.zip(self.affine(v)).map(|(k, v)| v.map(|c| mul_o(k, c)))
                }
                Node::Bin(BinOp::Div, a, b) if self.free(b) => {
                    let k = self.value(b).filter(|k| k.lo > 0.0 || k.hi < 0.0);
                    k.zip(self.affine(a)).map(|(k, v)| v.map(|c| div_o(c, k)))
                }
                _ => None,
            }
        };
        self.affine.insert(n, a);
        a
    }

    /// The strict sublevel set (`less`, else superlevel) of id at level,
    /// as prisms, never its closure. As boxes::bricks, with affine leaves.
    fn pieces(&mut self, id: NodeId, level: NodeId, less: bool) -> Option<Vec<Piece>> {
        self.budget = self.budget.checked_sub(1)?;
        if let Some(a) = self.affine(id) {
            return self.leaf(id, a, level, less);
        }
        let node = self.g.nodes[id as usize];
        match node {
            Node::Tag(_, a) | Node::Real(a) | Node::Lipschitz(a, ..) => self.pieces(a, level, less),
            Node::Un(UnOp::Neg, a) => {
                let level = self.g.neg(level);
                self.pieces(a, level, !less)
            }
            Node::Un(UnOp::Abs, a) => {
                let hi = self.pieces(a, level, less)?;
                let minus = self.g.neg(level);
                let lo = self.pieces(a, minus, !less)?;
                self.combine(hi, lo, less)
            }
            Node::Bin(op @ (BinOp::Add | BinOp::Sub), a, b) if self.free(a) || self.free(b) => {
                let level = self.g.real(level);
                if self.free(b) {
                    let l = if op == BinOp::Add { self.g.sub(level, b) } else { self.g.add(level, b) };
                    self.pieces(a, l, less)
                } else {
                    let l = if op == BinOp::Add { self.g.sub(level, a) } else { self.g.sub(a, level) };
                    self.pieces(b, l, if op == BinOp::Add { less } else { !less })
                }
            }
            Node::Bin(BinOp::Mul, a, b) if self.free(a) || self.free(b) => {
                let (k, v) = if self.free(a) { (a, b) } else { (b, a) };
                let s = self.value(k)?;
                if !(s.lo > 0.0 || s.hi < 0.0) {
                    return None;
                }
                let level = self.g.real(level);
                let l = self.g.bin_real(BinOp::Div, level, k);
                self.pieces(v, l, less == (s.lo > 0.0))
            }
            Node::Bin(op @ (BinOp::Min | BinOp::Max | BinOp::Emax), a, b) if op != BinOp::Emax || self.g.as_const(level) == Some(0.0) => {
                let left = self.pieces(a, level, less)?;
                let right = self.pieces(b, level, less)?;
                self.combine(left, right, (op != BinOp::Min) == less)
            }
            _ => None,
        }
    }

    /// a x + b y + c z + e below level (`less`) or above it: a side, a
    /// bound on z, all of space or none. A face leaning in z is not read.
    fn leaf(&mut self, id: NodeId, [a, b, c, e]: [I; 4], level: NodeId, less: bool) -> Option<Vec<Piece>> {
        let l = self.value(level)?;
        let (n, rhs) = if less { ([a, b, c], sub_o(l, e)) } else { ([a, b, c].map(neg), sub_o(e, l)) };
        let mut piece = Piece::default();
        match n.map(zero) {
            [true, true, true] if rhs.lo > 0.0 => {}
            [true, true, true] if rhs.hi <= 0.0 => return Some(vec![]),
            [_, _, true] if !(zero(n[0]) && zero(n[1])) => {
                let along = self.along(id);
                piece.sides.push(Side { n: [n[0], n[1]], c: rhs, key: (id, level, less), along });
            }
            [true, true, false] if n[2].lo > 0.0 => piece.z[1] = Some(div_o(rhs, n[2])),
            [true, true, false] if n[2].hi < 0.0 => piece.z[0] = Some(div_o(rhs, n[2])),
            _ => return None,
        }
        Some(vec![piece])
    }

    fn combine(&mut self, mut a: Vec<Piece>, b: Vec<Piece>, both: bool) -> Option<Vec<Piece>> {
        if !both {
            a.extend(b);
            return Some(a);
        }
        self.budget = self.budget.checked_sub(a.len().checked_mul(b.len())?)?;
        Some(a.iter().flat_map(|x| b.iter().map(move |y| intersect(x, y))).collect())
    }
}

/// An upper bound on e · p over the open polygon: a sum λ_i c_i of two sides
/// with λ_i n_i summing to e, λ at least 0 (weak duality: it holds wherever
/// the polygon is, bounded or not shown bounded). None: no pair shows one.
fn reach(sides: &[Side], e: [f64; 2]) -> Option<f64> {
    let e = e.map(I::point);
    let cross = |a: [I; 2], b: [I; 2]| sub_o(mul_o(a[0], b[1]), mul_o(a[1], b[0]));
    let mut best = f64::INFINITY;
    for (i, s) in sides.iter().enumerate() {
        // One side alone: n a positive multiple of e.
        if zero(cross(s.n, e)) {
            let k = dot(s.n, e[0], e[1]);
            if k.lo > 0.0 {
                best = best.min(div_o(s.c, k).hi);
            }
        }
        for t in &sides[i + 1..] {
            let det = cross(s.n, t.n);
            if !(det.lo > 0.0 || det.hi < 0.0) {
                continue;
            }
            let (a, b) = (div_o(cross(e, t.n), det), div_o(cross(s.n, e), det));
            if a.lo >= 0.0 && b.lo >= 0.0 {
                best = best.min(add_o(mul_o(a, s.c), mul_o(b, t.c)).hi);
            }
        }
    }
    best.is_finite().then_some(best)
}

/// The rule one roof's footprint is checked by.
struct Roof<'a> {
    sides: &'a [Side],
    /// Each bearing's sides that are not the roof's own.
    bearings: Vec<Vec<Side>>,
    dirs: &'a [[I; 2]],
    /// Whether direction `i` steps strictly parallel to `dirs[i]` (an axis
    /// direction, or unrounded real expressions such as `sin-cos-deg`): a
    /// diagonal direction of plain folded floats rounds each step's x and y
    /// offsets separately, so a row beside a parallel face can drift across it.
    exact_step: &'a [bool],
    step: I,
    k: i64,
}

impl Roof<'_> {
    fn over_bearing(&self, x: I, y: I) -> bool {
        self.bearings.iter().any(|b| b.iter().all(|s| dot(s.n, x, y).hi < s.c.lo))
    }

    /// Every point of the cell in the roof has its step j along direction i
    /// in the roof (moving away from none of its sides, or along one
    /// exactly when the rows step along d without per-component rounding, or
    /// past them shown) and over a bearing: a support column the full depth.
    fn supported(&self, x: I, y: I, i: usize, j: i64) -> bool {
        let d = self.dirs[i];
        let t = mul_o(I::point(j as f64), self.step);
        let (x, y) = (add_o(x, mul_o(t, d[0])), add_o(y, mul_o(t, d[1])));
        let kept = |s: &Side| (self.exact_step[i] && ((i < 64 && s.along >> i & 1 == 1) || mul_o(dot(s.n, d[0], d[1]), t).hi <= 0.0)) || dot(s.n, x, y).hi < s.c.lo;
        self.sides.iter().all(kept) && self.over_bearing(x, y)
    }

    /// Every point of the cell in the roof holds by the rows' rule: off the
    /// roof, over a bearing, or with a supported step j+ ahead and j- behind
    /// along a direction, j+ - j- at most k. Then the row runs in the roof
    /// (convex) from its point to both, and every window of k steps about
    /// the point holds one.
    fn holds(&self, x: I, y: I) -> bool {
        if self.sides.iter().any(|s| dot(s.n, x, y).lo >= s.c.hi) || self.over_bearing(x, y) {
            return true;
        }
        (0..self.dirs.len()).any(|i| {
            let ahead = (1..self.k).find(|&j| self.supported(x, y, i, j));
            let behind = (1..self.k).find(|&j| self.supported(x, y, i, -j));
            matches!((ahead, behind), (Some(a), Some(b)) if a + b <= self.k)
        })
    }
}

/// A positive constant where every piece of the material holds by the rows'
/// rule along one of `dirs` (unit, as std's rows step): on the bed, its
/// footprint cell by cell (Roof::holds) over bearings (pieces through its
/// bottom face, solid below it the full depth, and reaching the bed when the
/// first row is less than a full depth above it). None where that is not
/// shown.
#[allow(clippy::too_many_arguments)]
pub(super) fn unsupported(g: &mut Graph, f: NodeId, bed: NodeId, depth: NodeId, height: NodeId, limit: NodeId, step: NodeId, dirs: &[[NodeId; 2]]) -> Option<NodeId> {
    let mut free = HashMap::new();
    for id in g.reachable(&[f]) {
        let node = g.nodes[id as usize];
        let (kids, k) = node.children();
        let fixed = plain(node) && kids[..k].iter().all(|c| free[c]);
        free.insert(id, fixed);
    }
    let mut r =
        Reader { g, free, values: HashMap::new(), affine: HashMap::new(), exact: HashMap::new(), plane: HashMap::new(), half_planes: HashMap::new(), dirs: vec![], budget: 512 };
    r.dirs = dirs.iter().map(|&d| r.direction(d)).collect();
    // The thresholds of boxes::unsupported, as the rows read them.
    let at_h = r.g.add(bed, height);
    let at_depth = r.g.add(bed, depth);
    let first_row = r.g.min(at_h, at_depth);
    let needed = r.g.bin_real(BinOp::Add, bed, depth);
    let [bed, depth, step, limit, first_row, needed] = [bed, depth, step, limit, first_row, needed].map(|n| r.value(n));
    let (bed, depth, step, limit) = (bed?, depth?, step?, limit?);
    let full_depth_above_bed = needed?.hi <= first_row?.lo;
    let dir_nodes = dirs;
    let dirs: Vec<[I; 2]> = dir_nodes.iter().map(|&[x, y]| Some([r.value(x)?, r.value(y)?])).collect::<Option<_>>()?;
    let folded_const = |g: &Graph, mut id: NodeId| -> bool {
        loop {
            match g.nodes[id as usize] {
                Node::Const(_) => return true,
                Node::Tag(_, a) | Node::Real(a) | Node::Lipschitz(a, ..) => id = a,
                _ => return false,
            }
        }
    };
    let exact_step: Vec<bool> = dir_nodes.iter().zip(&dirs).map(|(&[nx, ny], d)| zero(d[0]) || zero(d[1]) || (!folded_const(r.g, nx) && !folded_const(r.g, ny))).collect();
    let k = div_o(limit, step);
    if !(step.lo > 0.0 && k.lo.floor() == k.hi.floor() && (2.0..=4096.0).contains(&k.lo.floor())) || dirs.is_empty() {
        return None;
    }
    let k = k.lo.floor() as i64;
    let zero = r.g.constant(0.0);
    let mut pieces = r.pieces(f, zero, true)?;
    // Certainly empty in z: no material.
    pieces.retain(|p| !matches!(p.z, [Some(lo), Some(hi)] if hi.hi <= lo.lo));
    if pieces.len() > 64 {
        return None;
    }
    // Cap total step-probes across cells, directions and k so a large k near
    // the 4096 limit cannot spend seconds at document load before giving up.
    let mut cells = (1usize << 17).min((1usize << 22) / ((k as usize) * dirs.len()).max(1));
    for b in &pieces {
        let [Some(z0), Some(z1)] = b.z else { return None };
        // On the bed, or no z between its ends.
        if z0.hi <= bed.lo || z1.hi <= z0.lo {
            continue;
        }
        let bottom = sub_o(z0, depth);
        let mut bearings = vec![];
        for s in &pieces {
            let [Some(s0), Some(s1)] = s.z else { continue };
            if s0.hi < bottom.lo && (full_depth_above_bed || s0.hi <= bed.lo) && z0.hi < s1.lo {
                // Its sides every point of the roof is inside are not read.
                bearings.push(s.sides.iter().filter(|t| !b.sides.iter().any(|u| r.implies(u, t))).copied().collect::<Vec<_>>());
            }
        }
        // A bearing under all of it.
        if bearings.iter().any(Vec::is_empty) {
            continue;
        }
        let [x1, y1, x0, y0] = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]].map(|e| reach(&b.sides, e));
        let (x0, x1, y0, y1) = (-x0?, x1?, -y0?, y1?);
        if x1 < x0 || y1 < y0 {
            continue;
        }
        let roof = Roof { sides: &b.sides, bearings, dirs: &dirs, exact_step: &exact_step, step, k };
        let finest = (x1 - x0).max(y1 - y0) * 2f64.powi(-16);
        let mut stack = vec![(I::new(x0, x1), I::new(y0, y1))];
        while let Some((x, y)) = stack.pop() {
            cells = cells.checked_sub(1)?;
            if roof.holds(x, y) {
                continue;
            }
            let (wide, tall) = (x.hi - x.lo, y.hi - y.lo);
            if wide.max(tall) <= finest {
                return None;
            }
            if wide >= tall {
                stack.extend([(I::new(x.lo, x.mid()), y), (I::new(x.mid(), x.hi), y)]);
            } else {
                stack.extend([(x, I::new(y.lo, y.mid())), (x, I::new(y.mid(), y.hi))]);
            }
        }
    }
    Some(r.g.constant(1.0))
}
