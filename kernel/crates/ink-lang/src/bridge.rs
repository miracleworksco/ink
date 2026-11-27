//! Continuous material paths for the overhang and bridge checks. All bounds are fields
//! of the row point, so they follow the shape and its parameter ranges.
use super::*;
mod boxes;
mod flank;
mod footing;
mod joints;
mod rates;
use ink_core::graph::Node;
use rates::{LineNorm, Rate, Rates};

fn spatial(v: Var) -> bool {
    matches!(v, Var::X | Var::Y | Var::Z)
}

/// Bounds on each subexpression over the box enclosing a segment.
/// Missing bounds propagate without putting infinity into field arithmetic.
#[derive(Default)]
struct Bounds {
    spans: HashMap<NodeId, Option<[NodeId; 2]>>,
    independent: HashSet<NodeId>,
}

impl Bounds {
    fn new(g: &mut Graph, field: NodeId, axes: &[(Var, [NodeId; 2])]) -> Self {
        let mut bounds = Self::default();
        for id in g.reachable(&[field]) {
            bounds.insert(g, id, axes);
        }
        bounds
    }

    fn insert(&mut self, g: &mut Graph, id: NodeId, axes: &[(Var, Range)]) {
        let Self { spans, independent } = self;
        let node = g.nodes[id as usize];
        let (kids, n) = node.children();
        let axis = axes.iter().find(|(v, _)| matches!(node, Node::Var(w) if w == *v));
        let fixed = axis.is_none() && kids[..n].iter().all(|k| independent.contains(k));
        let identity = match node {
            Node::Tag(_, a) | Node::Lipschitz(a, ..) | Node::Read(_, _, a, _) | Node::Deriv(a, ..) | Node::Real(a) => Some(a),
            _ => None,
        };
        let span = if fixed {
            independent.insert(id);
            Some([id, id])
        } else if let Some((_, range)) = axis {
            Some(*range)
        } else if let Some(a) = identity {
            if independent.contains(&a) {
                independent.insert(id);
            }
            spans[&a]
        } else {
            span(g, node, spans, independent)
        };
        spans.insert(id, span);
    }
}

type Range = [NodeId; 2];

fn product(g: &mut Graph, a: Range, b: Range) -> Range {
    let zero = g.constant(0.0);
    if a == [zero; 2] || b == [zero; 2] {
        return [zero; 2];
    }
    let [al, ah] = a;
    let [bl, bh] = b;
    let mut mul = |a, b| g.bin_real(BinOp::Mul, a, b);
    let (ll, lh, hl, hh) = (mul(al, bl), mul(al, bh), mul(ah, bl), mul(ah, bh));
    let (l, r) = (g.min(ll, lh), g.min(hl, hh));
    let (u, v) = (g.max(ll, lh), g.max(hl, hh));
    [g.min(l, r), g.max(u, v)]
}

// Interval arithmetic over the moving coordinates. The bounds remain
// fields of the row point. A spatial operation must retain real arithmetic
// even when an endpoint happens to simplify to constants or parameters.
fn span(g: &mut Graph, node: Node, spans: &HashMap<NodeId, Option<[NodeId; 2]>>, independent: &HashSet<NodeId>) -> Option<[NodeId; 2]> {
    Some(match node {
        Node::Un(UnOp::Neg, a) => {
            let [lo, hi] = spans[&a]?;
            [g.neg(hi), g.neg(lo)]
        }
        Node::Un(op @ (UnOp::Abs | UnOp::Square), a) => {
            let [lo, hi] = spans[&a]?;
            let minus_hi = g.neg(hi);
            let low = g.max(minus_hi, lo);
            let zero = g.constant(0.0);
            let minus_lo = g.neg(lo);
            let abs = [g.max(zero, low), g.max(minus_lo, hi)];
            if op == UnOp::Square {
                abs.map(|n| {
                    let n = g.real(n);
                    g.un(op, n)
                })
            } else {
                abs
            }
        }
        Node::Un(UnOp::Sqrt, a) => spans[&a]?.map(|n| {
            let n = g.real(n);
            g.un(UnOp::Sqrt, n)
        }),
        Node::Bin(op @ (BinOp::Add | BinOp::Min | BinOp::Max | BinOp::Emax), a, b) => {
            let (a, b) = (spans[&a]?, spans[&b]?);
            [0, 1].map(|i| {
                let n = g.real(a[i]);
                g.bin(op, n, b[i])
            })
        }
        Node::Bin(BinOp::Sub, a, b) => {
            let (a, b) = (spans[&a]?, spans[&b]?);
            [0, 1].map(|i| {
                let n = g.real(a[i]);
                g.sub(n, b[1 - i])
            })
        }
        Node::Select(c, a, b) if independent.contains(&c) => {
            let (a, b) = (spans[&a]?, spans[&b]?);
            [0, 1].map(|i| g.select(c, a[i], b[i]))
        }
        Node::Bin(BinOp::Div, a, b) if independent.contains(&b) => {
            let [lo, hi] = spans[&a]?.map(|a| g.bin_real(BinOp::Div, a, b));
            [g.min(lo, hi), g.max(lo, hi)]
        }
        Node::Bin(BinOp::Mul, a, b) if independent.contains(&a) || independent.contains(&b) => {
            let (k, v) = if independent.contains(&a) { (a, b) } else { (b, a) };
            let k = g.real(k);
            let [a, b] = spans[&v]?.map(|v| g.mul(k, v));
            [g.min(a, b), g.max(a, b)]
        }
        Node::Bin(BinOp::Mul, a, b) => product(g, spans[&a]?, spans[&b]?),
        _ => return None,
    })
}

#[derive(Clone, Copy)]
struct Probe<'a> {
    slopes: &'a HashMap<Var, NodeId>,
    rates: &'a Rates,
}

// Search for a level crossing as z increases in a vertical column, or
// as distance from the start increases on a descending probe. A leaf
// that cannot be inverted conservatively admits either crossing.
struct Faces<'a> {
    g: &'a mut Graph,
    field: NodeId,
    bed: Option<NodeId>,
    low: NodeId,
    high: NodeId,
    column: &'a Bounds,
    probe: Option<Probe<'a>>,
    memo: HashMap<(NodeId, NodeId, bool), NodeId>,
}

impl Faces<'_> {
    fn fixed(&self, id: NodeId) -> bool {
        self.column.independent.contains(&id)
    }

    // Keep a crossing at either endpoint. For a vertical roof search only,
    // the bed itself is excluded. On a slanted line t is distance in z.
    fn at_crossing(&mut self, t: NodeId) -> NodeId {
        let g = &mut self.g;
        let yes = g.constant(-1.0);
        let no = g.constant(1.0);
        let high = g.sub(self.high, t);
        let low = g.sub(t, self.low);
        let margin = g.min(low, high);
        let above_bed = self.bed.map_or(yes, |bed| {
            let bed = g.sub(bed, t);
            g.select(bed, yes, no)
        });
        let inside = g.select(margin, no, above_bed);
        if self.column.spans[&self.field].is_none() {
            return inside;
        }
        // All operations for which we build spans are continuous in z.
        // A zero of one leaf buried in another solid cannot be a face:
        // the original field itself must be zero at the crossing.
        let field = if let Some(probe) = self.probe {
            let at: Vec<_> = [Var::X, Var::Y, Var::Z]
                .into_iter()
                .filter_map(|v| {
                    let &rate = probe.slopes.get(&v)?;
                    let start = g.var(v);
                    let distance = g.bin_real(BinOp::Mul, rate, t);
                    Some((v, g.add(start, distance)))
                })
                .collect();
            g.subst(self.field, &at)
        } else {
            g.subst(self.field, &[(Var::Z, t)])
        };
        let value = g.un(UnOp::Abs, field);
        let minus = g.neg(value);
        g.select(minus, no, inside)
    }

    // A rise crosses the given level from below along the search. Invert
    // monotone scalar operations until the crossing coordinate is known. Across
    // abs/square, either sign can cross. At a min/max either active child
    // can cross. Anything we cannot invert remains a possible crossing.
    fn candidate(&mut self, id: NodeId, level: NodeId, rise: bool) -> NodeId {
        if self.fixed(id) {
            return self.g.constant(1.0);
        }
        if let Some(&found) = self.memo.get(&(id, level, rise)) {
            return found;
        }
        let n = self.g.nodes[id as usize];
        // An affine leaf on a probe crosses at (level - start) / slope.
        // A zero slope cannot cross; keep both endpoints for a graze.
        let affine = self.probe.and_then(|p| p.rates[&id]?.affine(self.g));
        let found = if let Some(slope) = affine {
            let level = self.g.real(level);
            let offset = self.g.sub(level, id);
            let t = self.g.bin(BinOp::Div, offset, slope);
            let at = self.at_crossing(t);
            let dir = if rise { self.g.neg(slope) } else { slope };
            let no = self.g.constant(1.0);
            self.g.select(dir, at, no)
        } else {
            match n {
                Node::Var(Var::Z) if self.probe.is_some() => self.g.constant(-1.0),
                Node::Var(Var::Z) if rise => self.at_crossing(level),
                Node::Var(Var::Z) => self.g.constant(1.0),
                Node::Tag(_, a) | Node::Lipschitz(a, ..) | Node::Read(_, _, a, _) | Node::Real(a) => self.candidate(a, level, rise),
                Node::Select(c, a, b) if self.fixed(c) => {
                    let a = self.candidate(a, level, rise);
                    let b = self.candidate(b, level, rise);
                    self.g.select(c, a, b)
                }
                Node::Bin(op @ (BinOp::Min | BinOp::Max), a, b) => self.either(a, b, level, rise, op == BinOp::Min),
                Node::Bin(BinOp::Emax, a, b) => {
                    // emax differs from max only in the strictly positive
                    // quadrant. At a nonpositive level, or with either
                    // operand nonpositive throughout, the same faces cross.
                    let planar = self.either(a, b, level, rise, false);
                    let mut threshold = level;
                    for child in [a, b] {
                        if let Some([_, hi]) = self.column.spans[&child] {
                            threshold = self.g.min(threshold, hi);
                        }
                    }
                    let positive = self.g.neg(threshold);
                    let undecided = self.g.constant(-1.0);
                    self.g.select(positive, undecided, planar)
                }
                Node::Un(UnOp::Neg, a) => {
                    let level = self.g.neg(level);
                    self.candidate(a, level, !rise)
                }
                Node::Un(UnOp::Sqrt, a) => {
                    let level = self.g.real(level);
                    let square = self.g.un(UnOp::Square, level);
                    self.candidate(a, square, rise)
                }
                Node::Un(op @ (UnOp::Abs | UnOp::Square), a) => {
                    let level = if op == UnOp::Square {
                        let level = self.g.real(level);
                        self.g.un(UnOp::Sqrt, level)
                    } else {
                        level
                    };
                    let left = self.candidate(a, level, rise);
                    let minus = self.g.neg(level);
                    let right = self.candidate(a, minus, !rise);
                    self.g.min(left, right)
                }
                Node::Bin(BinOp::Sub, a, b) if self.fixed(b) => {
                    let b = self.g.real(b);
                    let level = self.g.add(level, b);
                    self.candidate(a, level, rise)
                }
                Node::Bin(BinOp::Sub, a, b) if self.fixed(a) => {
                    let a = self.g.real(a);
                    let level = self.g.sub(a, level);
                    self.candidate(b, level, !rise)
                }
                Node::Bin(BinOp::Add, a, b) if self.fixed(a) || self.fixed(b) => {
                    let (f, k) = if self.fixed(a) { (b, a) } else { (a, b) };
                    let k = self.g.real(k);
                    let level = self.g.sub(level, k);
                    self.candidate(f, level, rise)
                }
                Node::Bin(BinOp::Mul, a, b) if self.fixed(a) || self.fixed(b) => {
                    let (k, f) = if self.fixed(a) { (a, b) } else { (b, a) };
                    let k = self.g.real(k);
                    let level = self.g.bin(BinOp::Div, level, k);
                    let pos = self.candidate(f, level, rise);
                    let neg = self.candidate(f, level, !rise);
                    let minus = self.g.neg(k);
                    let unknown = self.g.constant(-1.0);
                    let pos_or_zero = self.g.select(minus, pos, unknown);
                    self.g.select(k, neg, pos_or_zero)
                }
                _ => self.g.constant(-1.0),
            }
        };
        let Some([lo, hi]) = self.column.spans[&id] else {
            self.memo.insert((id, level, rise), found);
            return found;
        };
        let threshold = self.g.real(level);
        let a = self.g.sub(threshold, lo);
        let b = self.g.sub(hi, threshold);
        let margin = self.g.min(a, b);
        let no = self.g.constant(1.0);
        let found = self.g.select(margin, no, found);
        self.memo.insert((id, level, rise), found);
        found
    }

    fn either(&mut self, a: NodeId, b: NodeId, level: NodeId, rise: bool, union: bool) -> NodeId {
        let mut active = |leaf, other| {
            let cross = self.candidate(leaf, level, rise);
            let Some([lo, hi]) = self.column.spans[&other] else {
                return cross;
            };
            let level = self.g.real(level);
            let margin = if union { self.g.sub(hi, level) } else { self.g.sub(level, lo) };
            let no = self.g.constant(1.0);
            self.g.select(margin, no, cross)
        };
        let a_cross = active(a, b);
        let b_cross = active(b, a);
        self.g.min(a_cross, b_cross)
    }
}

pub(super) fn print_slope(g: &mut Graph, angle: NodeId) -> NodeId {
    let one = g.constant(1.0);
    // At 45 degrees the exact tangent is one. Keep all other angle
    // arithmetic real and rounded outward, including the diagonals.
    let angle = g.real(angle);
    let degrees = g.constant(45.0);
    let offset = g.sub(angle, degrees);
    let offset = g.un(UnOp::Abs, offset);
    let other = g.neg(offset);
    let radians = g.un(UnOp::Rad, angle);
    let tangent = g.un(UnOp::Tan, radians);
    g.select(other, tangent, one)
}

impl Doc {
    /// A sufficient condition for a segment to stay inside, given that its
    /// starting point is inside. At an intersection both operands start
    /// inside. At a union it suffices to retain one containing operand.
    fn print_segment(g: &mut Graph, field: NodeId, direction: [NodeId; 2], d: NodeId, bed: NodeId) -> (HashMap<NodeId, [NodeId; 2]>, Rates) {
        let zero = g.constant(0.0);
        let high = g.var(Var::Z);
        let low = g.sub(high, d);
        let direction = [direction[0], direction[1], g.constant(-1.0)];
        let mut axes = Vec::new();
        let mut slopes = HashMap::new();
        for (axis, rate) in [Var::X, Var::Y, Var::Z].into_iter().zip(direction) {
            if g.as_const(rate) == Some(0.0) {
                continue;
            }
            let start = g.var(axis);
            let travel = g.mul(rate, d);
            let end = g.add(start, travel);
            let mut lo = g.min(start, end);
            // Both callers clip the descent above the bed. Keep that
            // exact lower bound when a cell straddles the clipped end.
            if axis == Var::Z {
                lo = g.max(bed, lo);
            }
            axes.push((axis, [lo, g.max(start, end)]));
            slopes.insert(axis, rate);
        }
        let mut bounds = Bounds::default();
        let yes = g.constant(-1.0);
        let no = g.constant(1.0);
        let mut rates = Rates::new();
        let mut planes = HashMap::<NodeId, usize>::new();
        let mut norms = HashMap::new();
        let mut vertical = HashSet::new();
        let mut order = HashMap::new();
        let mut safe = HashMap::<NodeId, [NodeId; 2]>::new();
        for id in g.reachable(&[field]) {
            let n = g.nodes[id as usize];
            bounds.insert(g, id, &axes);
            let fixed = bounds.independent.contains(&id);
            if !fixed && let Some(line) = LineNorm::new(g, n, &rates, d) {
                bounds.spans.insert(id, Some([line.closest, g.max(line.start, line.end)]));
                norms.insert(id, line);
            }
            let rate = if fixed {
                Some(Rate { range: [zero; 2], guard: yes })
            } else if let Node::Un(UnOp::Sqrt, a) = n
                && let Some(line) = norms.get(&a)
            {
                let lower = g.neg(line.closest);
                Some(Rate { range: line.slopes, guard: g.select(lower, yes, no) })
            } else {
                Rate::along(g, n, &bounds, &rates, &slopes)
            };
            rates.insert(id, rate);
            // Count the branches Faces can actually invert. A nonlinear
            // leaf admits a crossing, but only where its range reaches zero.
            // It must not disable a crossing proof in a different union arm.
            let size = if fixed || rate.is_some_and(|r| r.affine(g).is_some()) {
                1
            } else {
                match n {
                    Node::Real(a) | Node::Tag(_, a) | Node::Lipschitz(a, ..) | Node::Read(_, _, a, _) | Node::Un(UnOp::Neg | UnOp::Sqrt, a) => planes[&a],
                    Node::Un(UnOp::Abs | UnOp::Square, a) => 2 * planes[&a],
                    Node::Bin(BinOp::Min | BinOp::Max | BinOp::Emax, a, b) | Node::Select(_, a, b) => planes[&a] + planes[&b],
                    Node::Bin(BinOp::Add | BinOp::Sub, a, b) if bounds.independent.contains(&a) || bounds.independent.contains(&b) => planes[&a] + planes[&b],
                    Node::Bin(BinOp::Mul, a, b) if bounds.independent.contains(&a) || bounds.independent.contains(&b) => 2 * (planes[&a] + planes[&b]),
                    _ => 1,
                }
            };
            planes.insert(id, size.min(65));
            let (children, count) = n.children();
            let z_only = !matches!(n, Node::Var(Var::X | Var::Y)) && children[..count].iter().all(|n| vertical.contains(n));
            if z_only {
                vertical.insert(id);
            }
            let monotone = rates::order(g, n, &bounds, &order, rate, fixed);
            order.insert(id, monotone);
            let flags = [0, 1].map(|side| {
                if monotone[side] == yes {
                    return yes;
                }
                let signed = |g: &mut Graph, v| if side == 0 { v } else { g.neg(v) };
                match n {
                    Node::Tag(_, a) | Node::Real(a) | Node::Lipschitz(a, ..) | Node::Deriv(a, ..) | Node::Read(_, _, a, _) => safe[&a][side],
                    Node::Select(c, a, b) if bounds.independent.contains(&c) => g.select(c, safe[&a][side], safe[&b][side]),
                    Node::Bin(op @ (BinOp::Mul | BinOp::Div), a, b) if bounds.independent.contains(&b) || (op == BinOp::Mul && bounds.independent.contains(&a)) => {
                        let (x, k) = if bounds.independent.contains(&b) { (a, b) } else { (b, a) };
                        // The strict sign of k*x or x/k fixes the sign
                        // of x; k is constant along this particular probe.
                        g.select(k, safe[&x][1 - side], safe[&x][side])
                    }
                    Node::Un(UnOp::Neg, a) => safe[&a][1 - side],
                    Node::Bin(op @ (BinOp::Min | BinOp::Max | BinOp::Emax), a, b) => {
                        let (sa, sb) = (safe[&a][side], safe[&b][side]);
                        let intersection = (op != BinOp::Min) == (side == 0);
                        if intersection {
                            g.max(sa, sb)
                        } else {
                            let (a, b) = (signed(g, a), signed(g, b));
                            let left = g.select(a, sa, no);
                            let right = g.select(b, sb, no);
                            let either = g.min(left, right);
                            // The start is in the union. Each operand that
                            // could contain it may instead be proved safe.
                            let left = g.select(a, sa, yes);
                            let right = g.select(b, sb, yes);
                            let every = g.max(left, right);
                            let single = g.min(every, either);
                            if planes[&id] > 64 {
                                return single;
                            }
                            let mut faces =
                                Faces { g, field: id, bed: None, low: zero, high: d, column: &bounds, probe: Some(Probe { slopes: &slopes, rates: &rates }), memo: HashMap::new() };
                            let crossing = faces.candidate(id, zero, side == 0);
                            let clear = g.neg(crossing);
                            g.min(single, clear)
                        }
                    }
                    _ => {
                        let mut proof = no;
                        if let Some([lo, hi]) = bounds.spans[&id] {
                            let bound = if side == 0 { hi } else { g.neg(lo) };
                            proof = g.select(bound, yes, no);
                        }
                        // Above the bed, this clipped segment stays strictly
                        // above it. A z-only constraint keeps its starting
                        // sign unless its boundary is crossed on the descent.
                        if z_only {
                            let mut faces = Faces { g, field: id, bed: Some(bed), low, high, column: &bounds, probe: None, memo: HashMap::new() };
                            let crossing = faces.candidate(id, zero, side == 1);
                            let clear = g.neg(crossing);
                            proof = g.min(proof, clear);
                        }
                        g.min(proof, monotone[side])
                    }
                }
            });
            safe.insert(id, flags);
        }
        (safe, rates)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn bridge_support(&mut self, f: NodeId, bed: NodeId, angle: NodeId, height: NodeId, depth: NodeId, limit: NodeId, step: NodeId) -> [NodeId; 3] {
        let g = &mut self.g;
        let f = joints::material(g, f);
        let high = g.var(Var::Z);
        let zero = g.constant(0.0);
        let yes = g.constant(-1.0);
        let no = g.constant(1.0);
        let low = g.sub(high, depth);
        let column = Bounds::new(g, f, &[(Var::Z, [low, high])]);
        let air = column.spans[&f].map_or(yes, |[_, hi]| g.select(hi, no, yes));
        let mut faces = Faces { g, field: f, bed: None, low, high, column: &column, probe: None, memo: HashMap::new() };
        // From a solid row point, encountering air while descending
        // requires a falling crossing when read in increasing z, even
        // when the field only touches zero. Include both endpoints.
        let enters = faces.candidate(f, zero, false);
        let bare_from_solid = faces.g.max(air, enters);
        faces.bed = Some(bed);
        faces.memo.clear();
        let roof = faces.candidate(f, zero, false);
        let boxes = boxes::unsupported(g, f, bed, depth, height, limit, step).unwrap_or(yes);
        let above = g.sub(high, bed);
        let above = g.max(zero, above);
        let half = g.constant(0.5);
        let half_clearance = g.mul(above, half);
        let unsupported = [height, depth].map(|d| {
            let d = g.un(UnOp::Abs, d);
            let d = g.min(d, half_clearance);
            Self::unsupported_paths(g, f, bed, angle, d, false)
        });
        let over = g.min(unsupported[0], unsupported[1]);
        let roof = g.max(roof, over);
        [bare_from_solid, g.max(roof, boxes), f]
    }

    pub(super) fn overhang_support(&mut self, f: NodeId, bed: NodeId, angle: NodeId, height: NodeId) -> NodeId {
        let g = &mut self.g;
        let f = joints::material(g, f);
        let z = g.var(Var::Z);
        let zero = g.constant(0.0);
        let above = g.sub(z, bed);
        let above = g.max(zero, above);
        let d = g.min(height, above);
        // A vertical path starting in material can only reach air above
        // the bed through a falling face (in increasing z), including grazes.
        let bottom = g.sub(z, height);
        let low = g.max(bed, bottom);
        let column = Bounds::new(g, f, &[(Var::Z, [low, z])]);
        let mut faces = Faces { g, field: f, bed: Some(bed), low, high: z, column: &column, probe: None, memo: HashMap::new() };
        let roof = faces.candidate(f, zero, false);
        if g.as_const(roof) == Some(1.0) {
            return roof;
        }
        let angled = Self::unsupported_paths(g, f, bed, angle, d, true);
        g.max(roof, angled)
    }

    fn directions(g: &mut Graph, angle: NodeId) -> [[NodeId; 2]; 8] {
        let zero = g.constant(0.0);
        let r = print_slope(g, angle);
        let half = g.constant(0.5);
        let half = g.real(half);
        let diagonal = g.un(UnOp::Sqrt, half);
        let q = g.mul(r, diagonal);
        let (nr, nq) = (g.neg(r), g.neg(q));
        [[r, zero], [q, q], [zero, r], [nq, q], [nr, zero], [nq, nq], [zero, nr], [q, nq]]
    }

    fn unsupported_paths(g: &mut Graph, f: NodeId, bed: NodeId, angle: NodeId, d: NodeId, vertical: bool) -> NodeId {
        let zero = g.constant(0.0);
        let yes = g.constant(-1.0);
        let no = g.constant(1.0);
        let directions = Self::directions(g, angle);
        let terms = flank::MaterialTerm::of(g, f);
        let mut pieces = if vertical { footing::members(g, f, bed) } else { vec![] };
        let mut supported = no;
        let mut covered = vec![no; terms.len()];
        for direction in vertical.then_some([zero; 2]).into_iter().chain(directions) {
            let (ray, rates) = Self::print_segment(g, f, direction, d, bed);
            supported = g.min(supported, ray[&f][0]);
            if supported == yes {
                return no;
            }
            for (piece, proof) in &mut pieces {
                *proof = g.min(*proof, ray[piece][0]);
            }
            for (term, proof) in terms.iter().zip(&mut covered) {
                let keeps = term.keeps_sign(g, &ray, &rates);
                *proof = g.min(*proof, keeps);
            }
        }
        if !pieces.is_empty() {
            let mut all = yes;
            for (piece, proof) in pieces {
                let ok = g.select(piece, proof, yes);
                all = g.max(all, ok);
            }
            supported = g.min(supported, all);
        }
        if !terms.is_empty() {
            let mut all = yes;
            for (term, proof) in terms.iter().zip(covered) {
                let ok = g.select(term.active, proof, yes);
                all = g.max(all, ok);
            }
            supported = g.min(supported, all);
        }
        g.neg(supported)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink_core::interval::I;

    fn interval(g: &Graph, f: NodeId) -> I {
        eval::eval_interval(&Tape::new(g, f), &Ctx::new(0.0, vec![]), [I::point(0.0); 3], None)
    }

    #[test]
    fn column_zero_test_keeps_real_arithmetic_after_substitution() {
        let mut g = Graph::new();
        let z = g.var(Var::Z);
        let big = g.constant(1.0e16);
        let half = g.constant(0.5);
        let high = g.constant(1.0);
        let low = g.constant(0.0);
        // At z = 1/2 the field is exactly zero. Rebuilding it with
        // nearest-float constant folding would instead give -1/2 and
        // discard a possible roof across the entire support column.
        let sum = g.add(big, z);
        let restored = g.sub(sum, big);
        let f = g.sub(restored, half);
        let column = Bounds::new(&mut g, f, &[(Var::Z, [low, high])]);
        let mut faces = Faces { g: &mut g, field: f, bed: None, low, high, column: &column, probe: None, memo: HashMap::new() };
        let possible = faces.at_crossing(half);
        assert!(interval(&g, possible).contains(-1.0));
    }

    #[test]
    fn rounded_intersection_uses_the_plane_rate_only_outside_the_round() {
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let y = g.var(Var::Y);
        let one = g.constant(1.0);
        let neg = g.constant(-1.0);
        let field = g.bin(BinOp::Emax, x, y);
        let next_x = g.add(x, one);
        let next_y = g.sub(y, one);
        let bounds = Bounds::new(&mut g, field, &[(Var::X, [x, next_x]), (Var::Y, [next_y, y])]);
        let rates = HashMap::from([(x, Some(Rate { range: [one; 2], guard: neg })), (y, Some(Rate { range: [neg; 2], guard: neg }))]);
        let node = g.nodes[field as usize];
        let rate = Rate::along(&mut g, node, &bounds, &rates, &HashMap::new()).unwrap();
        let at = |id, xy: [f64; 2]| eval::eval_interval(&Tape::new(&g, id), &Ctx::new(0.0, vec![]), [I::point(xy[0]), I::point(xy[1]), I::point(0.0)], None);
        assert!(at(rate.guard, [-3., 4.]).all_inside());
        for bound in rate.range {
            assert_eq!(at(bound, [-3., 4.]), I::point(-1.0));
        }
        assert!(at(rate.guard, [3., 4.]).none_inside());
        assert!(at(rate.guard, [-0.5, 4.]).none_inside());
    }

    #[test]
    fn moving_product_encloses_cancellation_and_both_segment_ends() {
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let zero = g.constant(0.0);
        let half = g.constant(0.5);
        let one = g.constant(1.0);
        let yes = g.constant(-1.0);
        let e = 2f64.powi(-27);
        let up = g.constant(1.0 + e);
        let down = g.constant(1.0 - e);
        let scaled = g.mul(up, x);
        let a = g.sub(scaled, one);
        let b = g.add(x, down);
        let field = g.mul(a, b);
        let bounds = Bounds::new(&mut g, field, &[(Var::X, [zero, half])]);
        let slopes = HashMap::from([(Var::X, one)]);
        let mut rates = Rates::new();
        for id in g.reachable(&[field]) {
            let node = g.nodes[id as usize];
            let rate = if bounds.independent.contains(&id) { Some(Rate { range: [zero; 2], guard: yes }) } else { Rate::along(&mut g, node, &bounds, &rates, &slopes) };
            rates.insert(id, rate);
        }
        // f'(x) = 2(1+e)x - e². Nearest rounding alone loses
        // (1+e)(1-e) - 1 at the left endpoint.
        let Rate { range: [lo, hi], guard } = rates[&field].unwrap();
        assert!(interval(&g, guard).hi < 0.0);
        assert!(interval(&g, lo).lo <= -e * e);
        assert!(interval(&g, hi).hi >= 1.0 + e);

        // A square root whose argument can hit zero supplies no rate.
        let root = g.un(UnOp::Sqrt, x);
        let rate = Rate::along(&mut g, Node::Un(UnOp::Sqrt, x), &bounds, &rates, &slopes).unwrap();
        assert!(interval(&g, rate.guard).lo > 0.0);
        let nested = g.mul(root, b);
        let bounds = Bounds::new(&mut g, nested, &[(Var::X, [zero, half])]);
        rates.insert(root, Some(rate));
        let rate = Rate::along(&mut g, Node::Bin(BinOp::Mul, root, b), &bounds, &rates, &slopes).unwrap();
        assert!(interval(&g, rate.guard).lo > 0.0);
    }
}
