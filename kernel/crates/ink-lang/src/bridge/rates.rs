//! Outward bounds on change per millimetre of descent.
use super::*;

// The endpoints bound the directional derivative everywhere on the segment
// where the guard is negative. In particular ink's total sqrt has no finite
// derivative bound at zero. Never let a parent operation drop that guard.
#[derive(Clone, Copy)]
pub(super) struct Rate {
    pub(super) range: Range,
    pub(super) guard: NodeId,
}

pub(super) type Rates = HashMap<NodeId, Option<Rate>>;
impl Rate {
    pub(super) fn affine(self, g: &Graph) -> Option<NodeId> {
        (self.range[0] == self.range[1] && g.as_const(self.guard) == Some(-1.0)).then_some(self.range[0])
    }
    pub(super) fn along(g: &mut Graph, n: Node, bounds: &Bounds, rates: &Rates, slopes: &HashMap<Var, NodeId>) -> Option<Rate> {
        let yes = g.constant(-1.0);
        let (range, guard) = match n {
            Node::Var(v) => ([*slopes.get(&v)?; 2], yes),
            Node::Tag(_, a) | Node::Real(a) | Node::Lipschitz(a, ..) | Node::Deriv(a, ..) | Node::Read(_, _, a, _) => return rates[&a],
            Node::Un(op, a) => {
                let Rate { range: [dl, dh], mut guard } = rates[&a]?;
                let range = match op {
                    UnOp::Neg => [g.neg(dh), g.neg(dl)],
                    UnOp::Abs => {
                        let [lo, hi] = bounds.spans[&a]?;
                        let nl = g.neg(dh);
                        let nh = g.neg(dl);
                        let lower = g.min(dl, nl);
                        let upper = g.max(dh, nh);
                        let positive = g.neg(lo);
                        let l = g.select(positive, dl, lower);
                        let h = g.select(positive, dh, upper);
                        [g.select(hi, nl, l), g.select(hi, nh, h)]
                    }
                    UnOp::Square => {
                        let a = bounds.spans[&a]?.map(|n| {
                            let n = g.real(n);
                            g.add(n, n)
                        });
                        product(g, a, [dl, dh])
                    }
                    UnOp::Sqrt => {
                        let [lo, hi] = bounds.spans[&a]?;
                        let positive = g.neg(lo);
                        let no = g.constant(1.0);
                        let valid = g.select(positive, yes, no);
                        guard = g.max(guard, valid);
                        let half = g.constant(0.5);
                        let half = g.real(half);
                        let inv = [hi, lo].map(|n| {
                            let n = g.real(n);
                            let root = g.un(UnOp::Sqrt, n);
                            g.bin(BinOp::Div, half, root)
                        });
                        product(g, [dl, dh], inv)
                    }
                    _ => return None,
                };
                (range, guard)
            }
            Node::Select(c, a, b) if bounds.independent.contains(&c) => {
                let (a, b) = (rates[&a]?, rates[&b]?);
                ([0, 1].map(|i| g.select(c, a.range[i], b.range[i])), g.select(c, a.guard, b.guard))
            }
            Node::Bin(BinOp::Div, a, b) if bounds.independent.contains(&b) => {
                let rate = rates[&a]?;
                let [lo, hi] = rate.range.map(|v| g.bin_real(BinOp::Div, v, b));
                ([g.min(lo, hi), g.max(lo, hi)], rate.guard)
            }
            Node::Bin(op @ (BinOp::Min | BinOp::Max | BinOp::Emax), a, b) => {
                let (ra, rb) = (rates[&a]?, rates[&b]?);
                let (va, vb) = (bounds.spans[&a]?, bounds.spans[&b]?);
                let a_below = g.bin_real(BinOp::Sub, va[1], vb[0]);
                let b_below = g.bin_real(BinOp::Sub, vb[1], va[0]);
                let (only_a, only_b) = if op == BinOp::Min { (a_below, b_below) } else { (b_below, a_below) };
                let range = [0, 1].map(|i| {
                    let hull = if i == 0 { g.min(ra.range[i], rb.range[i]) } else { g.max(ra.range[i], rb.range[i]) };
                    let v = g.select(only_b, rb.range[i], hull);
                    g.select(only_a, ra.range[i], v)
                });
                let both = g.max(ra.guard, rb.guard);
                let branch_guard = g.select(only_b, rb.guard, both);
                let mut guard = g.select(only_a, ra.guard, branch_guard);
                // emax is exactly max when either operand stays nonpositive.
                // Its rounded positive quadrant needs a different derivative.
                if op == BinOp::Emax {
                    let hi = g.min(va[1], vb[1]);
                    let positive = g.neg(hi);
                    let no = g.constant(1.0);
                    guard = g.select(positive, no, guard);
                }
                (range, guard)
            }
            Node::Bin(op @ (BinOp::Add | BinOp::Sub | BinOp::Mul), a, b) => {
                let (ra, rb) = (rates[&a]?, rates[&b]?);
                let guard = g.max(ra.guard, rb.guard);
                let (da, db) = (ra.range, rb.range);
                let range = match op {
                    BinOp::Mul => {
                        let left = product(g, da, bounds.spans[&b]?);
                        let right = product(g, db, bounds.spans[&a]?);
                        [0, 1].map(|i| {
                            let v = g.real(left[i]);
                            g.add(v, right[i])
                        })
                    }
                    _ => [0, 1].map(|i| {
                        let v = g.real(da[i]);
                        if op == BinOp::Add { g.add(v, db[i]) } else { g.sub(v, db[1 - i]) }
                    }),
                };
                (range, guard)
            }
            _ => return None,
        };

        Some(Rate { range, guard })
    }
}

/// Conditions under which a value cannot increase or cannot decrease on
/// the whole segment. A monotone composition needs no finite derivative
/// at a zero radius or where a positive part starts moving.
pub(super) fn order(g: &mut Graph, node: Node, bounds: &Bounds, order: &HashMap<NodeId, Range>, rate: Option<Rate>, fixed: bool) -> Range {
    let yes = g.constant(-1.0);
    let no = g.constant(1.0);
    if fixed {
        return [yes; 2];
    }
    [0, 1].map(|i| {
        let slope = rate.map_or(no, |r| {
            let away = if i == 0 { g.neg(r.range[1]) } else { r.range[0] };
            let valid = g.select(away, no, yes);
            g.max(r.guard, valid)
        });
        if slope == yes {
            return yes;
        }
        let composed = match node {
            Node::Real(a) | Node::Lipschitz(a, ..) | Node::Read(_, _, a, _) | Node::Tag(_, a) | Node::Deriv(a, ..) | Node::Un(UnOp::Sqrt, a) => order[&a][i],
            Node::Un(UnOp::Neg, a) => order[&a][1 - i],
            Node::Bin(BinOp::Add, a, b) => g.max(order[&a][i], order[&b][i]),
            Node::Bin(op @ (BinOp::Min | BinOp::Max | BinOp::Emax), a, b) => {
                let both = g.max(order[&a][i], order[&b][i]);
                if let (Some(va), Some(vb)) = (bounds.spans[&a], bounds.spans[&b]) {
                    let active = |g: &mut Graph, x: Range, y: Range, keeps| {
                        let below = if op == BinOp::Min { g.bin_real(BinOp::Sub, x[1], y[0]) } else { g.bin_real(BinOp::Sub, y[1], x[0]) };
                        let candidate = g.select(below, keeps, both);
                        if op == BinOp::Emax {
                            let positive = g.neg(y[1]);
                            g.select(positive, both, candidate)
                        } else {
                            candidate
                        }
                    };
                    let left = active(g, va, vb, order[&a][i]);
                    let right = active(g, vb, va, order[&b][i]);
                    g.min(left, right)
                } else {
                    both
                }
            }
            Node::Bin(BinOp::Sub, a, b) => g.max(order[&a][i], order[&b][1 - i]),
            Node::Un(UnOp::Abs | UnOp::Square, a) => {
                if let Some([lo, hi]) = bounds.spans[&a] {
                    let reverse = g.neg(hi);
                    let left = g.select(reverse, no, order[&a][1 - i]);
                    let right = g.select(lo, no, order[&a][i]);
                    g.min(left, right)
                } else {
                    no
                }
            }
            Node::Bin(op @ (BinOp::Mul | BinOp::Div), a, b) if bounds.independent.contains(&b) || (op == BinOp::Mul && bounds.independent.contains(&a)) => {
                let (x, k) = if bounds.independent.contains(&b) { (a, b) } else { (b, a) };
                g.select(k, order[&x][1 - i], order[&x][i])
            }
            _ => no,
        };
        g.min(composed, slope)
    })
}

// A sum of two squares of affine coordinates is minimized at the closest
// point on the line, clamped to the segment. This retains the correlation
// between coordinates which a box around a slanted segment loses.
#[derive(Clone, Copy)]
pub(super) struct LineNorm {
    pub(super) start: NodeId,
    pub(super) end: NodeId,
    pub(super) closest: NodeId,
    pub(super) slopes: Range,
}

// Along an affine line the squared norm is quadratic and its square root
// is convex. Its minimum is the perpendicular distance to the line when
// the projection lies inside the segment, or the nearer endpoint otherwise.
// The norm's derivative increases, so its two endpoint slopes enclose it.
impl LineNorm {
    pub(super) fn new(g: &mut Graph, n: Node, rates: &Rates, d: NodeId) -> Option<LineNorm> {
        let Node::Bin(BinOp::Add, a, b) = n else { return None };
        let component = |id: NodeId| {
            let Node::Un(UnOp::Square, v) = g.nodes[id as usize] else { return None };
            // A constants-and-params square in the source may be rounded
            // before the spatial addition. Only re-express real squares.
            if !g.depends_on(v, spatial) {
                return None;
            }
            Some((v, rates[&v]?.affine(g)?))
        };
        let ((a, u), (b, v)) = (component(a)?, component(b)?);

        let uu = g.un_real(UnOp::Square, u);
        let vv = g.un_real(UnOp::Square, v);
        let speed = g.bin_real(BinOp::Add, uu, vv);
        let ua = g.bin_real(BinOp::Mul, u, a);
        let vb = g.bin_real(BinOp::Mul, v, b);
        let dot = g.bin_real(BinOp::Add, ua, vb);
        let va = g.bin_real(BinOp::Mul, v, a);
        let ub = g.bin_real(BinOp::Mul, u, b);
        let cross = g.bin_real(BinOp::Sub, va, ub);
        let cross2 = g.un_real(UnOp::Square, cross);
        let distance = g.bin_real(BinOp::Div, cross2, speed);
        let du = g.bin_real(BinOp::Mul, u, d);
        let dv = g.bin_real(BinOp::Mul, v, d);
        let ae = g.bin_real(BinOp::Add, a, du);
        let be = g.bin_real(BinOp::Add, b, dv);
        let a2 = g.un_real(UnOp::Square, ae);
        let b2 = g.un_real(UnOp::Square, be);
        let end = g.bin_real(BinOp::Add, a2, b2);
        let a2 = g.un_real(UnOp::Square, a);
        let b2 = g.un_real(UnOp::Square, b);
        let start = g.bin_real(BinOp::Add, a2, b2);
        let sd = g.bin_real(BinOp::Mul, speed, d);
        let edot = g.bin_real(BinOp::Add, dot, sd);
        let interior = g.select(edot, end, distance);
        let closest = g.select(dot, interior, start);
        let slopes = [(dot, start), (edot, end)].map(|(v, n)| {
            let norm = g.un_real(UnOp::Sqrt, n);
            g.bin_real(BinOp::Div, v, norm)
        });
        Some(LineNorm { start, end, closest, slopes })
    }
}
