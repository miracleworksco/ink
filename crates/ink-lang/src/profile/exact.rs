//! A flat `polygon` whose edges are all straight, built exactly: every
//! decision of its construction is an exact predicate on its points (each
//! float an exact rational), and every line it emits passes through its
//! two corners in real numbers, so its sign is the drawn profile's wherever
//! interval arithmetic decides it, at any magnitude, with no tolerance.
//! Arcs, rounds and drafts are still polygon's construction (profile.rs,
//! walls.rs), within its allowance.
//!
//! The construction is profile.rs's: the profile is its convex hull less
//! the pockets between the two, each pocket read the same way, a pocket
//! leaving out its lid. Exactly: a corner is dropped only where it lies
//! exactly on the segment between its neighbours; a point lies on a hull
//! line only if exactly on it (it then splits the pocket there); a loop
//! that doubles back, crosses or touches itself is refused. Predicates are
//! filtered: a float test with its error bound decides most, expansions
//! (exact sums of floats, Shewchuk's) the rest.
//!
//! A line through u and v is emitted as ((x - ux) dy - (y - uy) dx) s, its
//! differences dx = vx - ux and dy = vy - uy each the exact sum of two
//! floats (a two-sum, written as two terms so nothing folds them), s any
//! positive float near 1 / |v - u|: zero exactly on the line through u and
//! v, whatever s rounds to.

use super::*;

/// An expansion: floats whose exact sum is the value, nonoverlapping, in
/// increasing magnitude, no zeros.
type E = Vec<f64>;

fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    (p, a.mul_add(b, -p))
}

/// e + b, exactly (Shewchuk's grow-expansion, zeros eliminated).
fn grow(e: &[f64], b: f64) -> E {
    let mut q = b;
    let mut h = Vec::with_capacity(e.len() + 1);
    for &x in e {
        let (s, err) = two_sum(q, x);
        if err != 0.0 {
            h.push(err);
        }
        q = s;
    }
    if q != 0.0 {
        h.push(q);
    }
    h
}

fn sum(e: &[f64], f: &[f64]) -> E {
    f.iter().fold(e.to_vec(), |acc, &b| grow(&acc, b))
}

fn mul(e: &[f64], f: &[f64]) -> E {
    let mut out = E::new();
    for &a in e {
        for &b in f {
            let (p, err) = two_prod(a, b);
            out = grow(&grow(&out, err), p);
        }
    }
    out
}

fn neg(e: &[f64]) -> E {
    e.iter().map(|x| -x).collect()
}

/// a - b, exactly.
fn diff(a: f64, b: f64) -> E {
    grow(&[a], -b)
}

fn sign(e: &[f64]) -> i8 {
    match e.last() {
        None => 0,
        Some(&x) => x.signum() as i8,
    }
}

/// Whether expansions of products of these differences are exact: every
/// difference 0 or within 2^-480 to 2^480 in size, so no product or its
/// error term leaves the normal floats.
fn safe(ds: &[E]) -> bool {
    ds.iter().flatten().all(|&v| v == 0.0 || (v.abs() >= 2f64.powi(-480) && v.abs() <= 2f64.powi(480)))
}

/// The sign of (b - a) x (c - a): positive where a, b, c turn
/// counter-clockwise. A float test with its error bound decides most (where
/// its products are far from underflow); expansions the rest, or integers
/// (`big`) where an expansion's products could leave the normal floats.
fn orient(a: P, b: P, c: P) -> i8 {
    let l = (b[0] - a[0]) * (c[1] - a[1]);
    let r = (b[1] - a[1]) * (c[0] - a[0]);
    let det = l - r;
    // (3 + 16 eps) eps of the terms bounds the float determinant's error
    // (Shewchuk's ccwerrboundA), doubled: differences of floats far apart
    // round, and this test is only a filter.
    let bound = 7e-16 * (l.abs() + r.abs());
    if bound >= 2f64.powi(-900) {
        if det > bound {
            return 1;
        }
        if -det > bound {
            return -1;
        }
    }
    let ds = [diff(b[0], a[0]), diff(c[1], a[1]), diff(b[1], a[1]), diff(c[0], a[0])];
    if !safe(&ds) {
        return big::orient(a, b, c);
    }
    let lhs = mul(&ds[0], &ds[1]);
    let rhs = mul(&ds[2], &ds[3]);
    sign(&sum(&lhs, &neg(&rhs)))
}

/// The sign of (b - a) . (c - b): positive where a -> b -> c runs on.
fn onward(a: P, b: P, c: P) -> i8 {
    let ds = [diff(b[0], a[0]), diff(c[0], b[0]), diff(b[1], a[1]), diff(c[1], b[1])];
    if !safe(&ds) {
        return big::onward(a, b, c);
    }
    let x = mul(&ds[0], &ds[1]);
    let y = mul(&ds[2], &ds[3]);
    sign(&sum(&x, &y))
}

/// The predicates in integers: each float an integer times a power of two,
/// all of one call's floats read at their least power, so sums and products
/// are exact at any magnitudes.
mod big {
    use super::P;
    use std::cmp::Ordering;

    /// A signed integer: its sign and its magnitude's 32-bit limbs, least
    /// first.
    #[derive(Clone)]
    pub(super) struct Z {
        neg: bool,
        mag: Vec<u32>,
    }

    fn trim(mut m: Vec<u32>) -> Vec<u32> {
        while m.last() == Some(&0) {
            m.pop();
        }
        m
    }

    fn cmp_mag(a: &[u32], b: &[u32]) -> Ordering {
        a.len().cmp(&b.len()).then_with(|| a.iter().rev().cmp(b.iter().rev()))
    }

    fn add_mag(a: &[u32], b: &[u32]) -> Vec<u32> {
        let mut out = Vec::with_capacity(a.len().max(b.len()) + 1);
        let mut carry = 0u64;
        for i in 0..a.len().max(b.len()) {
            let s = *a.get(i).unwrap_or(&0) as u64 + *b.get(i).unwrap_or(&0) as u64 + carry;
            out.push(s as u32);
            carry = s >> 32;
        }
        out.push(carry as u32);
        trim(out)
    }

    /// a - b, |a| >= |b|.
    fn sub_mag(a: &[u32], b: &[u32]) -> Vec<u32> {
        let mut out = Vec::with_capacity(a.len());
        let mut borrow = 0i64;
        for (i, &x) in a.iter().enumerate() {
            let mut d = x as i64 - *b.get(i).unwrap_or(&0) as i64 - borrow;
            borrow = 0;
            if d < 0 {
                d += 1 << 32;
                borrow = 1;
            }
            out.push(d as u32);
        }
        trim(out)
    }

    impl Z {
        /// m 2^shift for a float's integer mantissa m.
        pub(super) fn of(m: u64, neg: bool, shift: u32) -> Z {
            let mut mag = vec![0u32; (shift / 32) as usize];
            let (q, r) = (m as u128).checked_shl(shift % 32).map(|v| (v as u64 as u32, (v >> 32) as u64)).unwrap_or((0, 0));
            mag.push(q);
            mag.push(r as u32);
            mag.push((r >> 32) as u32);
            Z { neg, mag: trim(mag) }
        }
        pub(super) fn sub(&self, o: &Z) -> Z {
            self.add(&Z { neg: !o.neg, mag: o.mag.clone() })
        }
        pub(super) fn add(&self, o: &Z) -> Z {
            if self.neg == o.neg {
                return Z { neg: self.neg, mag: add_mag(&self.mag, &o.mag) };
            }
            match cmp_mag(&self.mag, &o.mag) {
                Ordering::Less => Z { neg: o.neg, mag: sub_mag(&o.mag, &self.mag) },
                _ => Z { neg: self.neg, mag: sub_mag(&self.mag, &o.mag) },
            }
        }
        pub(super) fn mul(&self, o: &Z) -> Z {
            let mut out = vec![0u64; self.mag.len() + o.mag.len() + 1];
            for (i, &x) in self.mag.iter().enumerate() {
                let mut carry = 0u64;
                for (j, &y) in o.mag.iter().enumerate() {
                    let t = out[i + j] + x as u64 * y as u64 + carry;
                    out[i + j] = t & 0xffff_ffff;
                    carry = t >> 32;
                }
                let mut k = i + o.mag.len();
                while carry != 0 {
                    let t = out[k] + carry;
                    out[k] = t & 0xffff_ffff;
                    carry = t >> 32;
                    k += 1;
                }
            }
            Z { neg: self.neg != o.neg, mag: trim(out.into_iter().map(|v| v as u32).collect()) }
        }
        pub(super) fn sign(&self) -> i8 {
            if self.mag.is_empty() {
                0
            } else if self.neg {
                -1
            } else {
                1
            }
        }
    }

    /// A finite float as (integer mantissa, exponent, sign): v = m 2^e.
    pub(super) fn parts(v: f64) -> (u64, i32, bool) {
        let bits = v.to_bits();
        let (neg, exp, frac) = (bits >> 63 == 1, ((bits >> 52) & 0x7ff) as i32, bits & ((1 << 52) - 1));
        if exp == 0 { (frac, -1074, neg) } else { (frac | (1 << 52), exp - 1075, neg) }
    }

    /// Each coordinate as an integer at the least power of two among them.
    pub(super) fn ints<const N: usize>(vs: [f64; N]) -> [Z; N] {
        let ps = vs.map(parts);
        let least = ps.iter().filter(|p| p.0 != 0).map(|p| p.1).min().unwrap_or(0);
        ps.map(|(m, e, neg)| if m == 0 { Z { neg: false, mag: vec![] } } else { Z::of(m, neg, (e - least) as u32) })
    }

    pub(super) fn orient(a: P, b: P, c: P) -> i8 {
        let [ax, ay, bx, by, cx, cy] = ints([a[0], a[1], b[0], b[1], c[0], c[1]]);
        bx.sub(&ax).mul(&cy.sub(&ay)).sub(&by.sub(&ay).mul(&cx.sub(&ax))).sign()
    }

    pub(super) fn onward(a: P, b: P, c: P) -> i8 {
        let [ax, ay, bx, by, cx, cy] = ints([a[0], a[1], b[0], b[1], c[0], c[1]]);
        bx.sub(&ax).mul(&cx.sub(&bx)).add(&by.sub(&ay).mul(&cy.sub(&by))).sign()
    }
}

fn lex_less(a: P, b: P) -> bool {
    a[0] < b[0] || (a[0] == b[0] && a[1] < b[1])
}

/// The loop's points, repeats and corners exactly on the segment between
/// their neighbours dropped, counter-clockwise; refused where it doubles
/// back on itself or has no area.
fn points(input: &[[f64; 3]]) -> Result<Vec<P>, String> {
    let mut pts: Vec<P> = Vec::with_capacity(input.len());
    for p in input {
        let q = [p[0], p[1]];
        if pts.last() != Some(&q) {
            pts.push(q);
        }
    }
    while pts.len() > 1 && pts[0] == pts[pts.len() - 1] {
        pts.pop();
    }
    loop {
        let n = pts.len();
        if n < 3 {
            return Err(format!("a polygon needs 3 corners not in line; it has {n}"));
        }
        let mut keep: Vec<P> = Vec::with_capacity(n);
        let mut dropped = false;
        for i in 0..n {
            let a = *keep.last().unwrap_or(&pts[n - 1]);
            let (b, c) = (pts[i], pts[(i + 1) % n]);
            if orient(a, b, c) == 0 {
                if onward(a, b, c) <= 0 {
                    return Err(format!("the polygon doubles back on itself at ({}, {})", b[0], b[1]));
                }
                dropped = true;
                continue;
            }
            keep.push(b);
        }
        pts = keep;
        if !dropped {
            break;
        }
    }
    let n = pts.len();
    if let Some((i, j)) = meet(&pts) {
        let (a, b) = (pts[i], pts[j]);
        return Err(format!("the polygon crosses or touches itself: its edges from ({}, {}) and from ({}, {}) meet", a[0], a[1], b[0], b[1]));
    }
    let lo = (0..n).fold(0, |m, i| if lex_less(pts[i], pts[m]) { i } else { m });
    if orient(pts[(lo + n - 1) % n], pts[lo], pts[(lo + 1) % n]) < 0 {
        pts.reverse();
    }
    Ok(pts)
}

/// Two edges of the loop, not neighbours, that meet (cross or touch),
/// exactly: by the edges' boxes swept along x, each pair the boxes leave
/// tested by orientation.
fn meet(pts: &[P]) -> Option<(usize, usize)> {
    let n = pts.len();
    let bx = |i: usize| {
        let (a, b) = (pts[i], pts[(i + 1) % n]);
        ([a[0].min(b[0]), a[1].min(b[1])], [a[0].max(b[0]), a[1].max(b[1])])
    };
    let mut ix: Vec<usize> = (0..n).collect();
    ix.sort_by(|&i, &j| bx(i).0[0].total_cmp(&bx(j).0[0]));
    let within = |p: P, a: P, b: P| a[0].min(b[0]) <= p[0] && p[0] <= a[0].max(b[0]) && a[1].min(b[1]) <= p[1] && p[1] <= a[1].max(b[1]);
    for (k, &i) in ix.iter().enumerate() {
        let (lo, hi) = bx(i);
        for &j in &ix[k + 1..] {
            let (lo2, hi2) = bx(j);
            if lo2[0] > hi[0] {
                break;
            }
            if lo2[1] > hi[1] || hi2[1] < lo[1] || (i + 1) % n == j || (j + 1) % n == i {
                continue;
            }
            let (a, b, c, d) = (pts[i], pts[(i + 1) % n], pts[j], pts[(j + 1) % n]);
            let (o1, o2, o3, o4) = (orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b));
            let touch = (o1 == 0 && within(c, a, b)) || (o2 == 0 && within(d, a, b)) || (o3 == 0 && within(a, c, d)) || (o4 == 0 && within(b, c, d));
            if touch || (o1 * o2 < 0 && o3 * o4 < 0) {
                return Some((i.min(j), i.max(j)));
            }
        }
    }
    None
}

/// The indices of the strict convex hull of `pts` (no point on an edge),
/// in loop order.
fn hull(pts: &[P]) -> Vec<usize> {
    let mut ix: Vec<usize> = (0..pts.len()).collect();
    ix.sort_by(|&a, &b| pts[a][0].total_cmp(&pts[b][0]).then(pts[a][1].total_cmp(&pts[b][1])));
    let mut h: Vec<usize> = Vec::new();
    for pass in 0..2 {
        let start = h.len();
        let it: Box<dyn Iterator<Item = &usize>> = if pass == 0 { Box::new(ix.iter()) } else { Box::new(ix.iter().rev()) };
        for &i in it {
            while h.len() >= start + 2 && orient(pts[h[h.len() - 2]], pts[h[h.len() - 1]], pts[i]) <= 0 {
                h.pop();
            }
            h.push(i);
        }
        h.pop();
    }
    h.sort_unstable();
    h.dedup();
    h
}

/// The tree of a loop (counter-clockwise, simple): its hull's half-planes
/// (but its lid, the edge from its last point to its first, if `lid`),
/// less its pockets.
fn piece(pts: &[P], lid: bool, out: &mut Vec<Leaf>) -> Tree {
    let n = pts.len();
    let h = hull(pts);
    let mut terms = Vec::new();
    for k in 0..h.len() {
        let (i, j) = (h[k], h[(k + 1) % h.len()]);
        let (a, b) = (pts[i], pts[j]);
        if !(lid && i == n - 1 && j == 0) {
            out.push(Leaf { u: a, v: b });
            terms.push(Tree::Leaf(out.len() - 1));
        }
        // The chain from i to j, cut where it touches the hull line.
        let len = (j + n - i) % n;
        let mut seg = vec![a];
        for s in 1..=len {
            let p = pts[(i + s) % n];
            seg.push(p);
            if s == len || orient(a, b, p) == 0 {
                if seg.len() > 2 {
                    seg.reverse();
                    terms.push(Tree::Neg(Box::new(piece(&seg, true, out))));
                }
                seg = vec![p];
            }
        }
    }
    Tree::Max(terms)
}

struct Leaf {
    u: P,
    v: P,
}

enum Tree {
    Leaf(usize),
    Max(Vec<Tree>),
    Neg(Box<Tree>),
}

impl Doc {
    /// A flat profile (x, y, bulge at each point): built exactly where its
    /// edges are straight, else by polygon's construction.
    pub(super) fn flat_polygon(&mut self, input: &[[f64; 3]]) -> Result<NodeId, String> {
        if input.iter().any(|p| p[2] != 0.0) {
            let t = polygon(input)?;
            return self.lower_profile(&t, None);
        }
        // Far from 1 in size (past 2^+-400), built at a power-of-two scale
        // (exact) putting its largest coordinate near 1, so no product its
        // lines emit leaves the normal floats; read at that scale. (Its
        // predicates are exact at any size: `big`.)
        let big = input.iter().fold(0.0f64, |m, p| m.max(p[0].abs()).max(p[1].abs()));
        let far = big > 0.0 && !(2f64.powi(-400)..=2f64.powi(400)).contains(&big);
        let k = if far { (-big.log2().ceil()).clamp(-1000.0, 1000.0) as i32 } else { 0 };
        let sc = 2f64.powi(k);
        let pts: Vec<[f64; 3]> = input.iter().map(|p| [p[0] * sc, p[1] * sc, 0.0]).collect();
        points(&pts).map(|pts| {
            let mut leaves = Vec::new();
            let t = piece(&pts, false, &mut leaves);
            let (x0, y0) = (self.g.var(Var::X), self.g.var(Var::Y));
            let c = self.g.constant(sc);
            let (x, y) = (self.g.mul(c, x0), self.g.mul(c, y0));
            let f = self.lower_exact(&t, &leaves, [x, y]);
            let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
            for p in &pts {
                for k in 0..2 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
            let mut terms = vec![f];
            for (v, k) in [(x, 0), (y, 1)] {
                terms.push(self.axis_face(v, hi[k], true));
                terms.push(self.axis_face(v, lo[k], false));
            }
            let f = self.reduce(BinOp::Max, &terms);
            let back = self.g.constant(1.0 / sc);
            self.g.mul(back, f)
        })
    }

    fn lower_exact(&mut self, t: &Tree, leaves: &[Leaf], xy: [NodeId; 2]) -> NodeId {
        match t {
            Tree::Leaf(k) => self.exact_line(leaves[*k].u, leaves[*k].v, xy),
            Tree::Neg(a) => {
                let a = self.lower_exact(a, leaves, xy);
                self.g.neg(a)
            }
            Tree::Max(v) => {
                let mut ids = Vec::with_capacity(v.len());
                for a in v {
                    let id = self.lower_exact(a, leaves, xy);
                    if !ids.contains(&id) {
                        ids.push(id);
                    }
                }
                self.reduce(BinOp::Max, &ids)
            }
        }
    }

    /// A line square to an axis, at coordinate c of that axis (v): v - c,
    /// positive past c (`above`), or its negation; exact, at most
    /// 1-Lipschitz, one value for every edge, box side or slab written on
    /// that line.
    pub(super) fn axis_face(&mut self, v: NodeId, c: f64, above: bool) -> NodeId {
        let at = self.g.constant(c);
        let f = self.g.sub(v, at);
        if above { f } else { self.g.neg(f) }
    }

    /// ((x - ux) dy - (y - uy) dx) s: positive right of u -> v, zero on
    /// the line through them in real numbers; square to an axis, that
    /// axis's face (`axis_face`).
    fn exact_line(&mut self, u: P, v: P, [x, y]: [NodeId; 2]) -> NodeId {
        if u[0] == v[0] {
            return self.axis_face(x, u[0], v[1] > u[1]);
        }
        if u[1] == v[1] {
            return self.axis_face(y, u[1], v[0] < u[0]);
        }
        let (dx, dy) = (diff(v[0], u[0]), diff(v[1], u[1]));
        let (ux, uy) = (self.g.constant(u[0]), self.g.constant(u[1]));
        let (xu, yu) = (self.g.sub(x, ux), self.g.sub(y, uy));
        let mut acc: Option<NodeId> = None;
        for (w, d, sg) in [(xu, &dy, 1.0), (yu, &dx, -1.0)] {
            for &c in d.iter().rev() {
                let c = self.g.constant(sg * c);
                let term = self.g.mul(w, c);
                acc = Some(match acc {
                    None => term,
                    Some(a) => self.g.add(a, term),
                });
            }
        }
        let size = |e: &E| ink_core::interval::sum_up(e.iter().map(|v| v.abs()));
        let s = self.g.constant(ink_core::interval::inv_norm_down(&[size(&dx), size(&dy)]));
        let a = acc.expect("u and v differ");
        self.g.mul(s, a)
    }
}

/// A field's exact value at a point, for tests: a dyadic rational (an
/// integer times a power of two), the graph's operations taken exactly
/// (sums, products, min, max, negation of its constants and coordinates).
#[cfg(test)]
fn dyadic_at(g: &Graph, root: NodeId, p: [f64; 2]) -> (big::Z, i32) {
    use ink_core::graph::Node;
    type D = (big::Z, i32);
    fn of(v: f64) -> D {
        let (m, e, neg) = big::parts(v);
        (big::Z::of(m, neg, 0), e)
    }
    fn align(a: &D, b: &D) -> (big::Z, big::Z, i32) {
        let e = a.1.min(b.1);
        let up = |d: &D| d.0.mul(&big::Z::of(1, false, (d.1 - e) as u32));
        (up(a), up(b), e)
    }
    let mut vals: Vec<Option<D>> = vec![None; root as usize + 1];
    for id in 0..=root as usize {
        let get = |vals: &Vec<Option<D>>, a: NodeId| vals[a as usize].clone().expect("an exact value");
        let v = match &g.nodes[id] {
            Node::Const(c) => Some(of(*c)),
            Node::Var(Var::X) => Some(of(p[0])),
            Node::Var(Var::Y) => Some(of(p[1])),
            Node::Var(_) => Some(of(0.0)),
            Node::Tag(_, a) | Node::Lipschitz(a, _, _) => vals[*a as usize].clone(),
            Node::Un(UnOp::Neg, a) => vals[*a as usize].clone().map(|(z, e)| (big::Z::of(0, false, 0).sub(&z), e)),
            Node::Bin(op @ (BinOp::Add | BinOp::Sub | BinOp::Mul | BinOp::Min | BinOp::Max), a, b) if vals[*a as usize].is_some() && vals[*b as usize].is_some() => {
                let (x, y) = (get(&vals, *a), get(&vals, *b));
                Some(match op {
                    BinOp::Mul => (x.0.mul(&y.0), x.1 + y.1),
                    _ => {
                        let (zx, zy, e) = align(&x, &y);
                        match op {
                            BinOp::Add => (zx.add(&zy), e),
                            BinOp::Sub => (zx.sub(&zy), e),
                            BinOp::Min => (if zx.sub(&zy).sign() <= 0 { zx } else { zy }, e),
                            _ => (if zx.sub(&zy).sign() >= 0 { zx } else { zy }, e),
                        }
                    }
                })
            }
            _ => None,
        };
        vals[id] = v;
    }
    vals[root as usize].clone().expect("the field's exact value")
}

/// Whether |f(q) - f(p)| <= d exactly (f the field at root).
#[cfg(test)]
pub(super) fn differs_by_at_most(g: &Graph, root: NodeId, p: [f64; 2], q: [f64; 2], d: f64) -> bool {
    let (p, q) = (&dyadic_at(g, root, p), &dyadic_at(g, root, q));
    let (m, e, neg) = big::parts(d);
    let dd = (big::Z::of(m, neg, 0), e);
    let least = p.1.min(q.1).min(dd.1);
    let up = |z: &big::Z, e: i32| z.mul(&big::Z::of(1, false, (e - least) as u32));
    let diff = up(&q.0, q.1).sub(&up(&p.0, p.1));
    let dz = up(&dd.0, dd.1);
    diff.sub(&dz).sign() <= 0 && big::Z::of(0, false, 0).sub(&diff).sub(&dz).sign() <= 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orient_is_exact_where_floats_round() {
        // Three points on a line a float determinant reads as turning.
        let (a, b, c) = ([0.5, 0.5], [12.0, 12.0], [24.0, 24.0]);
        assert_eq!(orient(a, b, c), 0);
        let e = f64::EPSILON;
        assert_eq!(orient([0.5 + e, 0.5], b, c), -1);
        assert_eq!(orient([0.5, 0.5 + e], b, c), 1);
        // Far out, where a float determinant cancels to 0.
        let k = 6e16;
        assert_eq!(orient([k, k], [k + 8.0, k], [k, k + 8.0]), 1);
    }

    #[test]
    fn a_corner_off_its_neighbours_line_is_kept() {
        let pts = [[0.0, 0.0, 0.0], [500.0, -5e-7, 0.0], [1000.0, 0.0, 0.0], [1000.0, 1000.0, 0.0], [0.0, 1000.0, 0.0]];
        assert_eq!(points(&pts).unwrap().len(), 5);
        let on = [[0.0, 0.0, 0.0], [500.0, 0.0, 0.0], [1000.0, 0.0, 0.0], [1000.0, 1000.0, 0.0], [0.0, 1000.0, 0.0]];
        assert_eq!(points(&on).unwrap().len(), 4);
    }

    #[test]
    fn orientation_is_exact_across_any_magnitudes() {
        // A convex loop from 1e-100 to 1e70: one scale to its largest
        // coordinate leaves products of its small differences below the
        // floats, where expansions round; it read as doubling back.
        let pts = [[0.0, 0.0, 0.0], [1e-100, 0.0, 0.0], [2e-100, 1e-100, 0.0], [1e70, 1e70, 0.0], [0.0, 1e70, 0.0]];
        let big = 1e70f64;
        let sc = 2f64.powi(-(big.log2().ceil() as i32));
        let scaled: Vec<[f64; 3]> = pts.iter().map(|p| [p[0] * sc, p[1] * sc, 0.0]).collect();
        let p = points(&scaled).unwrap();
        assert_eq!(p.len(), 5);
        // The integers agree with the expansions where both are exact.
        let mut r = 0x9e3779b97f4a7c15u64;
        let mut f = || {
            r ^= r << 13;
            r ^= r >> 7;
            r ^= r << 17;
            (r >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        for _ in 0..2000 {
            let (a, b) = ([f(), f()], [f(), f()]);
            let c = [a[0] + (b[0] - a[0]) * 0.5 + f() * 1e-17, a[1] + (b[1] - a[1]) * 0.5];
            assert_eq!(big::orient(a, b, c), orient(a, b, c), "{a:?} {b:?} {c:?}");
            assert_eq!(big::onward(a, b, c), onward(a, b, c));
        }
    }

    #[test]
    fn a_unit_normal_s_scale_is_proven_at_most_one() {
        // s^2 (x^2 + y^2) <= 1 in integers, for s = inv_norm_down.
        let check = |x: f64, y: f64| {
            let s = ink_core::interval::inv_norm_down(&[x, y]);
            let [zs, zx, zy] = big::ints([s, x, y]);
            let lhs = zs.mul(&zs).mul(&zx.mul(&zx).add(&zy.mul(&zy)));
            // 1 at the same scale: s, x and y read at their least power e,
            // so s^2 (x^2 + y^2) carries 2^(4e); 1 is 2^(-4e) of that unit.
            let least = [s, x, y].iter().filter(|v| **v != 0.0).map(|v| big::parts(*v).1).min().unwrap();
            let one = big::Z::of(1, false, (-4 * least) as u32);
            assert!(lhs.sub(&one).sign() <= 0, "{x} {y}: s {s}");
        };
        check(3.0, 4.0);
        check(1.0, 0.0);
        check(0.1, 0.7);
        for k in 1..500 {
            check(k as f64 * 0.37, 1.0 / k as f64);
        }
        assert_eq!(ink_core::interval::inv_norm_down(&[0.0, 2.0]), 0.5);
    }

    #[test]
    fn a_loop_that_crosses_or_touches_itself_is_refused() {
        let bow = [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert!(points(&bow).is_err());
        // A corner touching the opposite edge exactly.
        let touch = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        assert!(points(&touch).is_err());
        // The same corner a float above it: simple.
        let near = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0], [1.0, 1e-100, 0.0], [0.0, 2.0, 0.0]];
        assert!(points(&near).is_ok());
    }
}
