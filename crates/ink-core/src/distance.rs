//! Where a field's value is the distance to its shape, read from how it is
//! built: growing a shape by d (docs/IR.md, Fields) is everything within d
//! of it, which f - d gives only where f is the distance outside it, and
//! shrinking (f + d) only where -f is the distance inside.
//!
//! The walk reads structure only, never a computed coefficient: x, y and z
//! moved by constants are planes, the distance both sides; the absolute
//! value of a distance is the distance to its 0s (a slab's faces), and the
//! root of a sum of squares of such over disjoint axes the distance to
//! where all are 0 (a point, the z axis, a circle about it); a distance
//! less a constant grows it (and a point's, an axis's or a plane's both
//! ways); a union is the distance outside where each operand is, an
//! intersection inside, and an `emax` of pieces over disjoint axes both
//! sides where they are. Samples of a signed distance (an imported mesh's,
//! `redistance`'s) read in place are walked as one, but never shown one:
//! between samples the field is their trilinear interpolant, whose slope
//! reaches sqrt 3 times theirs, so what would be exact is approximate, to
//! within the samples' error and spacing (`Grows::Sampled`). Anything
//! else (a turn, a scale, a polygon's edges) is not shown: its growth's
//! bound is the field's slope, from the kernel's gradient bounds and
//! promises.
//!
//! The shape grown is where f < 0 (sign.rs), so a 0 the shape keeps out
//! (a region read as 0, a flush mouth's sheet) must not be grown as if it
//! were inside: a growth is exact only where the inside reaches every
//! point read as 0 or less (`solid_out`).

use crate::graph::{BinOp, Graph, Node, NodeId, UnOp, Var};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

/// a / b rounded down (a > 0, b >= 0): the quotient, less one step where
/// its exact residual shows it rounded up. The residual a - q b is exact
/// only where it cannot underflow (a at least 2^-900); below, one step
/// down, since a nearest quotient is within half a step.
fn div_down(a: f64, b: f64) -> f64 {
    let q = a / b;
    if b == 0.0 || q.is_infinite() {
        f64::MAX
    } else if q < f64::MIN_POSITIVE {
        0.0
    } else if a < f64::from_bits(((1023 - 900) as u64) << 52) || q.mul_add(b, -a) > 0.0 {
        q.next_down()
    } else {
        q
    }
}

/// a + b, where floats hold it exactly.
fn exact_add(a: f64, b: f64) -> Option<f64> {
    let s = a + b;
    (s.is_finite() && s - a == b && s - b == a).then_some(s)
}

/// The axes (bits x, y, z) each node's value reads.
fn axes(g: &Graph, ids: &[NodeId]) -> HashMap<NodeId, u8> {
    let mut m: HashMap<NodeId, u8> = HashMap::default();
    for &id in ids {
        let node = g.nodes[id as usize];
        let (ch, k) = match node {
            Node::Lipschitz(a, ..) | Node::Deriv(a, ..) | Node::Read(_, _, a, _) => ([a, 0, 0], 1),
            n => n.children(),
        };
        let own = match node {
            Node::Var(Var::X) => 1,
            Node::Var(Var::Y) => 2,
            Node::Var(Var::Z) => 4,
            _ => 0,
        };
        let v = ch[..k].iter().fold(own, |acc, c| acc | m.get(c).copied().unwrap_or(0));
        m.insert(id, v);
    }
    m
}

/// Whether an operation keeps its operand's sign wherever it reads it, so
/// another field of the operand's sign leaves its own sign unchanged:
/// min, max, emax, negation, names and promises (of their field).
fn keeps_sign(node: Node, child: NodeId) -> bool {
    match node {
        Node::Bin(BinOp::Min | BinOp::Max | BinOp::Emax, ..) | Node::Un(UnOp::Neg, _) | Node::Tag(..) => true,
        Node::Lipschitz(a, ..) => a == child,
        _ => false,
    }
}

/// Per `max` in `root` read only through operations that keep its sign,
/// whose operands (through the maxes under it) fall into groups over
/// disjoint axes (a box of slabs, a profile cut to a height): the groups.
/// An `emax` of the groups' maxes has the max's sign everywhere, so `root`
/// keeps its own, and is the distance outside where the groups are.
fn square_maxes(g: &Graph, root: NodeId) -> HashMap<NodeId, Vec<Vec<NodeId>>> {
    let ids = g.reachable(&[root]);
    let ax = axes(g, &ids);
    let (mut top, mut other) = (HashSet::default(), HashSet::default());
    top.insert(root);
    for &id in ids.iter().rev() {
        let node = g.nodes[id as usize];
        let (ch, k) = node.children();
        let kept = top.contains(&id) && !other.contains(&id);
        for &c in &ch[..k] {
            if kept && keeps_sign(node, c) {
                top.insert(c);
            } else {
                other.insert(c);
            }
        }
    }
    let mut out = HashMap::default();
    for &id in &ids {
        if !matches!(g.nodes[id as usize], Node::Bin(BinOp::Max, ..)) || !top.contains(&id) || other.contains(&id) {
            continue;
        }
        let (mut ops, mut todo, mut seen) = (Vec::new(), vec![id], HashSet::default());
        while let Some(n) = todo.pop() {
            match g.nodes[n as usize] {
                Node::Bin(BinOp::Max, a, b) => todo.extend([b, a]),
                _ if seen.insert(n) => ops.push(n),
                _ => {}
            }
        }
        // Operands joined while they share an axis; what reads none joins
        // the first group.
        let mut groups: Vec<(u8, Vec<NodeId>)> = Vec::new();
        for o in ops {
            let (joined, mut rest): (Vec<_>, Vec<_>) = groups.into_iter().partition(|(a, _)| a & ax[&o] != 0);
            let grp = joined.into_iter().fold((ax[&o], vec![o]), |(a, mut m), (b, n)| {
                m.extend(n);
                (a | b, m)
            });
            rest.push(grp);
            groups = rest;
        }
        if let Some(i) = groups.iter().position(|(a, _)| *a == 0)
            && groups.len() > 1
        {
            let (_, c) = groups.remove(i);
            groups[0].1.extend(c);
        }
        if groups.len() > 1 {
            out.insert(id, groups.into_iter().map(|(_, m)| m).collect());
        }
    }
    out
}

/// `root` with each `max` of pieces over disjoint axes it reads only
/// through sign-keeping operations written as an `emax` of them: the same
/// shape everywhere, and the distance outside it where the pieces are. What
/// a growth reads, so a box of slabs grows as std's box does.
pub fn exact(g: &mut Graph, root: NodeId) -> NodeId {
    let maxes = square_maxes(g, root);
    if maxes.is_empty() { root } else { g.with_emax(root, &maxes) }
}

/// Where a field is shown to be the distance (module docs).
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Reach {
    /// Where f > 0, it is the distance to where f <= 0.
    pub out: bool,
    /// Where f < 0, -f is the distance to where f >= 0.
    pub inn: bool,
    /// Every point where f <= 0 is a limit of points where f < 0: no 0 the
    /// inside does not reach.
    pub solid_out: bool,
    /// Every point where f >= 0 is a limit of points where f > 0.
    pub solid_in: bool,
    /// f is the distance to a flat (a point, an axis, a plane) less this:
    /// a ball, a rod or a slab about it.
    flat: Option<f64>,
    /// f is an axis's plane, moved or not.
    plane: bool,
    /// What it is shown to be rests on samples of a distance taken as one
    /// (`sampled_distance`): approximately so.
    pub sampled: bool,
}

/// How a field grows or shrinks by its value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Grows {
    /// Exactly: everything within d (shrinking: everything farther than d
    /// inside).
    Exact,
    /// Its slope is at most L, so it reads at most L times the distance:
    /// a growth holds everything within d / L (this), and may reach
    /// farther; a shrinking keeps only what lies farther than |d| / L
    /// inside, and may lose more.
    Within(f64),
    /// As exactly as the samples of a distance it reads are one, which
    /// their trilinear interpolant is not (its slope reaches sqrt 3 times
    /// theirs): everything within d of what they sample, to within their
    /// error plus sqrt 3 / 2 of their spacing.
    Sampled,
    /// No bound: the value may read more than the distance.
    Unknown,
}

impl Reach {
    /// How `f - d` grows (d > 0) or shrinks (d < 0) the shape, `d` in
    /// [lo, hi] (bounds of the real amount), `slope` a bound on |grad f|
    /// everywhere. Only d's sign, so bounded, decides exactness; an amount
    /// not shown to be either side has no bound (a written 0 is no offset:
    /// the caller's).
    pub fn grows(&self, (lo, hi): (f64, f64), slope: f64) -> Grows {
        let (exact, least) = match () {
            _ if lo > 0.0 => (self.out && self.solid_out, lo),
            _ if hi < 0.0 => (self.inn, -hi),
            _ => return Grows::Unknown,
        };
        if exact && self.sampled {
            Grows::Sampled
        } else if exact {
            Grows::Exact
        } else if slope.is_finite() {
            Grows::Within(div_down(least, slope))
        } else {
            Grows::Unknown
        }
    }

    fn of_value(v: Option<f64>) -> Reach {
        // Distance to nothing is infinite: a constant reads less, and no
        // point of its own sign has another side.
        v.map_or(Reach::default(), |v| Reach { out: v <= 0.0, inn: v >= 0.0, solid_out: v != 0.0, solid_in: v != 0.0, ..Reach::default() })
    }

    const PLANE: Reach = Reach { out: true, inn: true, solid_out: true, solid_in: true, flat: None, plane: true, sampled: false };

    /// Samples of a signed distance (`sampled_distance`), walked as the
    /// distance both sides.
    const SAMPLED: Reach = Reach { out: true, inn: true, solid_out: true, solid_in: true, flat: None, plane: false, sampled: true };

    /// The distance to a flat less `off`.
    fn flat(off: f64) -> Reach {
        // off 0: the flat itself, no inside; below 0, nothing.
        Reach { out: off >= 0.0, inn: true, solid_out: off != 0.0, solid_in: true, flat: Some(off), plane: false, sampled: false }
    }

    fn neg(self) -> Reach {
        Reach { out: self.inn, inn: self.out, solid_out: self.solid_in, solid_in: self.solid_out, flat: None, plane: self.plane, sampled: false }
    }

    /// f + k.
    fn shift(self, k: f64) -> Reach {
        if k == 0.0 || self.plane {
            return self;
        }
        if let Some(off) = self.flat {
            return exact_add(off, -k).map_or(Reach::default(), Reach::flat);
        }
        if k < 0.0 {
            // Grown: the distance outside stays, and every point now read
            // as 0 or less is reached by points nearer.
            Reach { inn: false, solid_out: self.out, solid_in: false, ..self }
        } else {
            Reach { out: false, solid_in: self.inn, solid_out: false, ..self }
        }
    }

    /// |f|: where f is the distance both sides, the distance to its 0s.
    fn abs(self) -> Reach {
        if !(self.out && self.inn) {
            Reach::default()
        } else if self.plane {
            Reach::flat(0.0)
        } else {
            Reach { out: true, inn: true, solid_out: false, solid_in: self.solid_out || self.solid_in, flat: None, plane: false, sampled: false }
        }
    }
}

/// max(a, b), `disjoint` where they read disjoint axes (or one reads
/// none), `euclid` where written as emax.
fn max_of(a: Reach, b: Reach, disjoint: bool, euclid: bool) -> Reach {
    if euclid && !disjoint {
        return Reach::default();
    }
    // Over disjoint axes, the distance outside an intersection is the root
    // of its operands' squares where both are positive: emax.
    let product = disjoint && euclid;
    Reach { out: product && a.out && b.out, inn: a.inn && b.inn, solid_out: product && a.solid_out && b.solid_out, solid_in: a.solid_in && b.solid_in, ..Reach::default() }
}

/// `root`'s reach. Of what reads no coordinate it trusts only constants,
/// numbers as the model writes them (or folds them, as floats): never a
/// param, a moment or a real expression of them, whose float may not be
/// its value.
pub fn reach(g: &Graph, root: NodeId) -> Reach {
    let val = |n: NodeId| g.as_const(n);
    let ids = g.reachable(&[root]);
    let ax = axes(g, &ids);
    let mut m: HashMap<NodeId, Reach> = HashMap::default();
    // f as s h + k: the node h it offsets, its sign s and the constant k,
    // where floats hold k exactly.
    let parts = |n: NodeId| -> (NodeId, f64, f64) {
        let (mut h, mut s, mut k) = (n, 1.0, 0.0);
        let c = |x: NodeId, s: f64, k: f64| if ax[&x] == 0 { val(x).and_then(|v| exact_add(k, s * v)) } else { None };
        loop {
            match g.nodes[h as usize] {
                Node::Un(UnOp::Neg, a) if ax[&a] != 0 => (h, s) = (a, -s),
                Node::Bin(BinOp::Add, a, b) if let Some(k2) = c(a, s, k) => (h, k) = (b, k2),
                Node::Bin(BinOp::Add, a, b) if let Some(k2) = c(b, s, k) => (h, k) = (a, k2),
                Node::Bin(BinOp::Sub, a, b) if let Some(k2) = c(b, -s, k) => (h, k) = (a, k2),
                Node::Bin(BinOp::Sub, a, b) if let Some(k2) = c(a, s, k) => (h, s, k) = (b, -s, k2),
                Node::Tag(_, a) | Node::Real(a) => h = a,
                _ => return (h, s, k),
            }
        }
    };
    for &id in &ids {
        let r = |n: NodeId| m.get(&n).copied().unwrap_or_default();
        let signed = |h: NodeId, s: f64| if s < 0.0 { r(h).neg() } else { r(h) };
        let reach = if ax[&id] == 0 {
            Reach::of_value(val(id))
        } else {
            match g.nodes[id as usize] {
                Node::Var(_) => Reach::PLANE,
                Node::Read(gi, path, _, at) if sampled_distance(g, gi, path, at, &ax, &r) => Reach::SAMPLED,
                Node::Tag(_, a) | Node::Real(a) | Node::Deriv(a, ..) | Node::Lipschitz(a, ..) => r(a),
                Node::Un(UnOp::Neg, a) => r(a).neg(),
                Node::Un(UnOp::Abs, a) => r(a).abs(),
                Node::Un(UnOp::Sqrt, a) => roots(g, a, &ax, &r),
                Node::Bin(BinOp::Add | BinOp::Sub, ..) => match parts(id) {
                    (h, s, k) if h != id => signed(h, s).shift(k),
                    _ => Reach::default(),
                },
                Node::Bin(op @ (BinOp::Max | BinOp::Min | BinOp::Emax), a, b) => {
                    // min(a, b) = -max(-a, -b).
                    let sign = if op == BinOp::Min { -1.0 } else { 1.0 };
                    let (pa, pb) = (parts(a), parts(b));
                    let half = |u: f64, v: f64| exact_add(u, v).map(|w| 0.5 * w).filter(|h| 2.0 * h == exact_add(u, v).unwrap_or(f64::NAN));
                    let out = match (half(pa.2, -pb.2), half(pa.2, pb.2)) {
                        // s h + ka and -s h + kb: |s h + (ka - kb) / 2| + (ka + kb) / 2.
                        (Some(m), Some(w)) if op != BinOp::Emax && pa.0 == pb.0 && pa.1 == -pb.1 => signed(pa.0, sign * pa.1).shift(sign * m).abs().shift(sign * w),
                        _ => {
                            let (ra, rb) = if sign < 0.0 { (r(a).neg(), r(b).neg()) } else { (r(a), r(b)) };
                            max_of(ra, rb, ax[&a] & ax[&b] == 0, op == BinOp::Emax)
                        }
                    };
                    if sign < 0.0 { out.neg() } else { out }
                }
                _ => Reach::default(),
            }
        };
        m.insert(id, reach);
    }
    // Whatever it is shown to be may rest on samples it reads.
    Reach { sampled: m.values().any(|r| r.sampled), ..m[&root] }
}

/// Whether grid `gi` read along `path` at frame `at` is a sampled signed
/// distance where it was sampled: the value (no derivative) of samples of
/// the distance to a surface (an imported mesh's, `redistance`'s, which
/// continue past their box by the distance to it), read at x, y and z
/// each moved by a constant or mirrored, over disjoint axes (an isometry).
/// Taken as the distance, approximately (`Grows::Sampled`): no field ink
/// could build from them is nearer.
fn sampled_distance(g: &Graph, gi: u32, path: u32, at: NodeId, ax: &HashMap<NodeId, u8>, r: &dyn Fn(NodeId) -> Reach) -> bool {
    let grid = &g.grids[gi as usize];
    let distance = grid.source().is_some_and(|s| matches!(s.kind, "mesh" | "redistance")) && !grid.holds();
    let Node::Frame(c) = g.nodes[at as usize] else { return false };
    let mut seen = 0u8;
    let isometry = c.iter().all(|&u| {
        let a = ax.get(&u).copied().unwrap_or(0);
        let fresh = a != 0 && seen & a == 0;
        seen |= a;
        fresh && r(u).plane
    });
    distance && g.paths[path as usize].is_empty() && isometry
}

/// sqrt(a1^2 + a2^2 + ...): where each ai is the distance both sides and
/// they read disjoint axes, the distance to where all are 0.
fn roots(g: &Graph, sum: NodeId, ax: &HashMap<NodeId, u8>, r: &dyn Fn(NodeId) -> Reach) -> Reach {
    let (mut terms, mut todo) = (Vec::new(), vec![sum]);
    while let Some(n) = todo.pop() {
        match g.nodes[n as usize] {
            Node::Bin(BinOp::Add, a, b) => todo.extend([a, b]),
            Node::Un(UnOp::Square, a) if ax[&a] != 0 => terms.push(a),
            _ => return Reach::default(),
        }
    }
    let (mut seen, mut flat, mut solid_in) = (0u8, true, true);
    for &t in &terms {
        let rt = r(t);
        if !(rt.out && rt.inn) || seen & ax[&t] != 0 {
            return Reach::default();
        }
        seen |= ax[&t];
        flat &= rt.plane || rt.flat == Some(0.0);
        solid_in &= rt.solid_out || rt.solid_in;
    }
    if flat { Reach::flat(0.0) } else { Reach { out: true, inn: true, solid_out: false, solid_in, flat: None, plane: false, sampled: false } }
}

/// A bound on |grad f| everywhere at `ctx`: the kernel's gradient bounds
/// over all of space, promises included (interval::lipschitz).
pub fn slope(g: &Graph, f: NodeId, ctx: &crate::eval::Ctx) -> f64 {
    let tape = crate::tape::Tape::new(g, f);
    let l = crate::interval::lipschitz(&crate::eval::eval_interval_grad(&tape, ctx, [crate::interval::I::ENTIRE; 3]));
    if l.is_nan() { f64::INFINITY } else { l }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// How (offset-shape f d) grows, values of constants alone trusted.
    fn grows(g: &Graph, f: NodeId, d: f64) -> Grows {
        reach(g, f).grows((d, d), slope(g, f, &crate::eval::Ctx::new(0.0, Vec::new())))
    }

    /// A bound held no farther than r: true wherever the slope is r's.
    fn below(w: Grows, r: f64) -> bool {
        matches!(w, Grows::Within(s) if s <= r) || w == Grows::Unknown
    }

    fn value_at(g: &Graph, f: NodeId, p: [f64; 3]) -> f64 {
        let mut v = Vec::new();
        crate::eval::eval_points(&crate::tape::Tape::new(g, f), &crate::eval::Ctx::new(0.0, Vec::new()), &[p], &mut v);
        v[0]
    }

    fn c(g: &mut Graph, v: f64) -> NodeId {
        g.constant(v)
    }

    /// lo - q and q - hi: a slab of planes.
    fn slab(g: &mut Graph, q: Var, lo: f64, hi: f64) -> NodeId {
        let q = g.var(q);
        let (lo, hi) = (c(g, lo), c(g, hi));
        let a = g.sub(lo, q);
        let b = g.sub(q, hi);
        g.max(a, b)
    }

    fn len(g: &mut Graph, xs: &[NodeId]) -> NodeId {
        let sq: Vec<NodeId> = xs.iter().map(|&a| g.un(UnOp::Square, a)).collect();
        let s = sq[1..].iter().fold(sq[0], |s, &q| g.add(s, q));
        g.un(UnOp::Sqrt, s)
    }

    fn cube_of_slabs(g: &mut Graph, lo: [f64; 3], hi: [f64; 3]) -> NodeId {
        let s: Vec<NodeId> = [Var::X, Var::Y, Var::Z].iter().enumerate().map(|(k, &v)| slab(g, v, lo[k], hi[k])).collect();
        let m = g.max(s[0], s[1]);
        g.max(m, s[2])
    }

    #[test]
    fn slabs_and_box_grow_alike() {
        let mut g = Graph::new();
        let slabs = cube_of_slabs(&mut g, [-1.0; 3], [1.0; 3]);
        // As written, the max of slabs mitres its edges.
        assert!(matches!(grows(&g, slabs, 0.5), Grows::Within(_)));
        assert_eq!(grows(&g, slabs, -0.5), Grows::Exact);
        let e = exact(&mut g, slabs);
        assert_eq!(grows(&g, e, 0.5), Grows::Exact);
        // The same values as std's box, emax(|x| - 1, |y| - 1, |z| - 1).
        let one = c(&mut g, 1.0);
        let q: Vec<NodeId> = [Var::X, Var::Y, Var::Z]
            .iter()
            .map(|&v| {
                let x = g.var(v);
                let a = g.un(UnOp::Abs, x);
                g.sub(a, one)
            })
            .collect();
        let m = g.bin(BinOp::Emax, q[0], q[1]);
        let bx = g.bin(BinOp::Emax, m, q[2]);
        assert_eq!(grows(&g, bx, 0.5), Grows::Exact);
        for p in [[1.4, 1.4, 0.0], [1.2, -1.3, 1.1], [0.2, 0.3, -0.1], [-3.0, 0.5, 2.0]] {
            assert_eq!(value_at(&g, e, p), value_at(&g, bx, p), "at {p:?}");
        }
    }

    #[test]
    fn primitives_are_distances() {
        let mut g = Graph::new();
        let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
        let (r, h, big) = (c(&mut g, 2.0), c(&mut g, 3.0), c(&mut g, 5.0));
        let l3 = len(&mut g, &[x, y, z]);
        let sphere = g.sub(l3, r);
        let l2 = len(&mut g, &[x, y]);
        let rod = g.sub(l2, r);
        let az = g.un(UnOp::Abs, z);
        let cap = g.sub(az, h);
        let cyl = g.bin(BinOp::Emax, rod, cap);
        for (name, f) in [("sphere", sphere), ("cylinder", cyl)] {
            assert_eq!((grows(&g, f, 1.0), grows(&g, f, -1.0)), (Grows::Exact, Grows::Exact), "{name}");
        }
        // A torus: exact outside; inside, not shown.
        let core = g.sub(l2, big);
        let circle = len(&mut g, &[core, z]);
        let torus = g.sub(circle, r);
        assert_eq!(grows(&g, torus, 1.0), Grows::Exact);
        assert_ne!(grows(&g, torus, -1.0), Grows::Exact);
        // A cylinder's max is not its emax: its rim mitres; written so, it is.
        let mitred = g.max(rod, cap);
        assert_ne!(grows(&g, mitred, 1.0), Grows::Exact);
        let e = exact(&mut g, mitred);
        assert_eq!(grows(&g, e, 1.0), Grows::Exact);
    }

    /// Samples of a distance are never shown to be one, read in place or
    /// not: on |p| - 1 sampled at spacing 1, the interpolant is 0 at
    /// (0, 0, 1) but 0.0742 at (0.03, 0.03, 1.05), 0.0656 from there, so
    /// growing it by 0.07 leaves out a point within 0.07 of its shape.
    /// In place (moved or mirrored too) a growth is approximate, to within
    /// the samples; scaled, baked or held past its box, bounded by slope.
    #[test]
    fn samples_of_a_distance_are_never_shown_one() {
        let pts: Vec<f32> = (0..125).map(|l| crate::grid::ijk(l, [5; 3]).iter().map(|&i| (i as f64 - 2.0).powi(2)).sum::<f64>().sqrt() as f32 - 1.0).collect();
        let grid = |kind: &'static str| {
            let source =
                crate::grid::Source { what: "ball".into(), at: "test".into(), kind, res: 1.0, band: None, measured: None, note: None, unit: "mm", span: None, component: None };
            std::sync::Arc::new(crate::grid::GridData::new([-2.0; 3], 1.0, [5; 3], pts.clone()).with_source(source))
        };
        let mut g = Graph::new();
        let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let f = g.grid(grid("redistance"), [x, y, z]);
        let (on, p) = (value_at(&g, f, [0.0, 0.0, 1.0]), value_at(&g, f, [0.03, 0.03, 1.05]));
        let near = (2.0 * 0.03f64.powi(2) + 0.05f64.powi(2)).sqrt();
        assert!(on == 0.0 && near < 0.07 && p > 0.07, "{on} {p} {near}");
        for kind in ["mesh", "redistance"] {
            let f = g.grid(grid(kind), [x, y, z]);
            let (three, nz) = (c(&mut g, 3.0), g.neg(z));
            let xm = g.sub(x, three);
            let moved = g.grid(grid(kind), [y, xm, nz]);
            let one = c(&mut g, 1.0);
            let grown = g.sub(f, one);
            for (f, d) in [(f, 0.07), (f, -0.5), (moved, 0.07), (grown, 0.5)] {
                assert_eq!(grows(&g, f, d), Grows::Sampled, "{kind} by {d}");
            }
        }
        let two = c(&mut g, 2.0);
        let x2 = g.mul(x, two);
        let scaled = g.grid(grid("mesh"), [x2, y, z]);
        let baked = g.grid(grid("bake"), [x, y, z]);
        let held = g.grid(std::sync::Arc::new(std::sync::Arc::try_unwrap(grid("mesh")).ok().unwrap().holding()), [x, y, z]);
        for f in [scaled, baked, held] {
            assert!(matches!(grows(&g, f, 0.5), Grows::Within(_) | Grows::Unknown));
        }
    }

    #[test]
    fn a_zero_kept_out_is_not_grown() {
        // (min 0 x): air read as 0 wherever x > 0; a slab of width 0, a
        // sheet read as 0.
        let mut g = Graph::new();
        let (zero, x) = (c(&mut g, 0.0), g.var(Var::X));
        let f = g.min(zero, x);
        let sheet = slab(&mut g, Var::X, 1.0, 1.0);
        assert_ne!(grows(&g, f, 1.0), Grows::Exact);
        assert_ne!(grows(&g, sheet, 1.0), Grows::Exact);
    }

    #[test]
    fn a_subtract_grows_past_its_mouth() {
        // A pocket flush with the top: inside, the distance; outside, the
        // mouth's 0 and the rim's mitre.
        let mut g = Graph::new();
        let o = cube_of_slabs(&mut g, [0.0, 0.0, 0.0], [40.0, 40.0, 10.0]);
        let i = cube_of_slabs(&mut g, [10.0, 10.0, 5.0], [30.0, 30.0, 10.0]);
        let ni = g.neg(i);
        let pocket = g.max(o, ni);
        let e = exact(&mut g, pocket);
        assert_ne!(grows(&g, e, 1.0), Grows::Exact);
        assert_eq!(grows(&g, e, -1.0), Grows::Exact);
        // A union of separate boxes grows exactly; with the pocket, not.
        let far = cube_of_slabs(&mut g, [50.0, 0.0, 0.0], [60.0, 10.0, 10.0]);
        let near = cube_of_slabs(&mut g, [-20.0, 0.0, 0.0], [-10.0, 10.0, 10.0]);
        let boxes = g.min(near, far);
        let boxes = exact(&mut g, boxes);
        assert_eq!(grows(&g, boxes, 1.0), Grows::Exact);
        let both = g.min(pocket, far);
        let both = exact(&mut g, both);
        assert_ne!(grows(&g, both, 1.0), Grows::Exact);
    }

    #[test]
    fn a_max_under_arithmetic_is_left_as_written() {
        // max(|x|, |y|, |z|) - 1 is a cube: an emax there is a sphere.
        let mut g = Graph::new();
        let a: Vec<NodeId> = [Var::X, Var::Y, Var::Z]
            .iter()
            .map(|&v| {
                let x = g.var(v);
                g.un(UnOp::Abs, x)
            })
            .collect();
        let m = g.max(a[0], a[1]);
        let m = g.max(m, a[2]);
        let one = c(&mut g, 1.0);
        let cube = g.sub(m, one);
        let e = exact(&mut g, cube);
        assert_eq!(value_at(&g, e, [0.9; 3]), value_at(&g, cube, [0.9; 3]));
        assert_ne!(grows(&g, e, 0.1), Grows::Exact);
    }

    #[test]
    fn computed_coefficients_are_not_read() {
        // Planes nearly facing apart (a wedge whose tip is far out),
        // turned, scaled or stretched by a hair: none is shown a distance,
        // and the bound is the slope's.
        let mut g = Graph::new();
        let (x, y, t) = (g.var(Var::X), g.var(Var::Y), g.var(Var::T));
        let e = c(&mut g, 1e-4);
        let a = g.sub(x, e);
        let nx = g.neg(x);
        let wedge = |g: &mut Graph, k: NodeId| {
            let ky = g.mul(k, y);
            let s = g.add(nx, ky);
            let b = g.sub(s, e);
            g.max(a, b)
        };
        let k = c(&mut g, 5e-10);
        let w = wedge(&mut g, k);
        assert!(below(grows(&g, w, 1.0), 1.0));
        // (1 + t) - 1, read at a moment whose float of it is 0.
        let one = c(&mut g, 1.0);
        let tp = g.add(one, t);
        let k = g.sub(tp, one);
        let w = wedge(&mut g, k);
        assert_ne!(reach(&g, w).grows((1.0, 1.0), 1.0), Grows::Exact);
        // A hair over unit slope: held only within 1 / slope.
        for v in [1.0000000000009, 2.0] {
            let k = c(&mut g, v);
            let kx = g.mul(k, x);
            let m = g.max(kx, y);
            assert!(below(grows(&g, kx, 1.0), 1.0 / v) && below(grows(&g, m, 1.0), 1.0 / v), "{v}");
        }
        // max(x, a y) / b shrunk: kept only farther than 1 / its slope in.
        let (va, vb) = (1.0000000000009, 0.9999999999991);
        let (ka, kb) = (c(&mut g, va), c(&mut g, vb));
        let ay = g.mul(ka, y);
        let m = g.max(x, ay);
        let f = g.bin(BinOp::Div, m, kb);
        assert!(below(grows(&g, f, -1.0), vb / va));
    }

    #[test]
    fn a_slope_is_the_kernels() {
        // |2 min(x, y)| - 1, emax(x, x / 2) and 2e200 (|1e-200 x| - 1e-200):
        // slopes 2, up to sqrt(1.25), and 2: held within d over those.
        let mut g = Graph::new();
        let (x, y) = (g.var(Var::X), g.var(Var::Y));
        let (one, two, half) = (c(&mut g, 1.0), c(&mut g, 2.0), c(&mut g, 0.5));
        let m = g.min(x, y);
        let m2 = g.mul(two, m);
        let a = g.un(UnOp::Abs, m2);
        let f = g.sub(a, one);
        assert!(below(grows(&g, f, 0.2), 0.1));
        let hx = g.mul(half, x);
        let f = g.bin(BinOp::Emax, x, hx);
        assert!(below(grows(&g, f, 1.0), 1.0 / 1.25f64.sqrt()));
        let (k, big) = (c(&mut g, 1e-200), c(&mut g, 1e100));
        let kx = g.mul(k, x);
        let ak = g.un(UnOp::Abs, kx);
        let s = g.sub(ak, k);
        let s = g.mul(big, s);
        let s = g.mul(big, s);
        let f = g.mul(two, s);
        assert!(below(grows(&g, f, 0.5), 0.25));
        // sin x promised 1-Lipschitz: a bound; warped to sin 2x, none.
        let sn = g.un(UnOp::Sin, x);
        let p = g.lipschitz(sn, one);
        assert!(below(grows(&g, p, 0.5), 0.5) && matches!(grows(&g, p, 0.5), Grows::Within(r) if r > 0.4999));
        let x2 = g.mul(two, x);
        let w = g.subst(p, &[(Var::X, x2)]);
        assert!(below(grows(&g, w, 0.5), 0.25));
    }

    #[test]
    fn a_radius_is_rounded_down() {
        // d / L rounded down: its product with L never passes d.
        // Checked scaled by 2^1000, where the residual cannot underflow.
        let big = 2f64.powi(1000);
        for (d, l) in [(1.0, 10.0), (2.728787161065714, 2.000000000000004), (0.3, 3.0), (1.0, 3.0), (2.3e-308, 2.5e-308), (2.3e-300, 2.5e-300)] {
            let Grows::Within(r) = Reach::default().grows((d, d), l) else { panic!() };
            let (d, l) = if d < 1e-200 { (d * big, l * big) } else { (d, l) };
            assert!(r.mul_add(l, -d) <= 0.0, "{d} / {l}: {r}");
        }
        // An amount whose sign is not shown has no bound, 0 included.
        assert_eq!(Reach::default().grows((-1.0, 1.0), 1.0), Grows::Unknown);
        assert_eq!(Reach::default().grows((0.0, 0.0), 1.0), Grows::Unknown);
    }

    #[test]
    fn a_halving_that_rounds_is_not_read() {
        // max(x, a - x), a = 2^-1074: |x - a/2| + a/2, never 0; halved in
        // floats, a/2 is 0, which would read it as |x|, a flat's distance.
        let mut g = Graph::new();
        let (x, y) = (g.var(Var::X), g.var(Var::Y));
        let a = c(&mut g, f64::from_bits(1));
        let ax = g.sub(a, x);
        let f = g.max(x, ax);
        let l = len(&mut g, &[f, y]);
        let lens = g.sub(l, a);
        assert_ne!(grows(&g, lens, f64::from_bits(1)), Grows::Exact);
    }
}
