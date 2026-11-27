//! A straight flank can exclude the inward crossing of an adjoining arc.
//! Work in the circle's own two affine coordinates; the half-plane may
//! come from any other constraint of the same material intersection.
use super::*;

#[derive(Clone, Copy)]
struct Flank {
    axes: [NodeId; 2],
    radius: NodeId,
    // The row point satisfies nx*x + ny*y + offset > 0.
    normal: [NodeId; 3],
}

fn circle(g: &Graph, id: NodeId, side: usize) -> Option<([NodeId; 2], NodeId)> {
    let Node::Bin(BinOp::Sub, a, b) = g.nodes[id as usize] else { return None };
    let (n, radius) = if side == 1 { (a, b) } else { (b, a) };
    let Node::Un(UnOp::Sqrt, sum) = g.nodes[n as usize] else { return None };
    let Node::Bin(BinOp::Add, a, b) = g.nodes[sum as usize] else { return None };
    let square = |id: NodeId| match g.nodes[id as usize] {
        Node::Un(UnOp::Square, v) if g.depends_on(v, spatial) => Some(v),
        _ => None,
    };
    if g.depends_on(radius, spatial) {
        return None;
    }
    Some(([square(a)?, square(b)?], radius))
}

// Write the constraint as an affine expression of the two coordinates.
// A spatial remainder, nonlinear operation or exhausted budget declines.
fn plane(g: &mut Graph, id: NodeId, axes: [NodeId; 2], budget: &mut usize) -> Option<[NodeId; 3]> {
    *budget = budget.checked_sub(1)?;
    let (zero, one) = (g.constant(0.0), g.constant(1.0));
    if id == axes[0] {
        return Some([one, zero, zero]);
    }
    if id == axes[1] {
        return Some([zero, one, zero]);
    }
    if !g.depends_on(id, spatial) {
        return Some([zero, zero, id]);
    }
    match g.nodes[id as usize] {
        Node::Real(a) | Node::Tag(_, a) | Node::Lipschitz(a, ..) => plane(g, a, axes, budget),
        Node::Un(UnOp::Neg, a) => Some(plane(g, a, axes, budget)?.map(|v| g.neg(v))),
        Node::Bin(op @ (BinOp::Add | BinOp::Sub), a, b) => {
            let (a, b) = (plane(g, a, axes, budget)?, plane(g, b, axes, budget)?);
            Some([0, 1, 2].map(|i| g.bin_real(op, a[i], b[i])))
        }
        Node::Bin(BinOp::Mul, a, b) => {
            let (v, k) = if !g.depends_on(a, spatial) {
                (b, a)
            } else if !g.depends_on(b, spatial) {
                (a, b)
            } else {
                return None;
            };
            Some(plane(g, v, axes, budget)?.map(|c| g.bin_real(BinOp::Mul, c, k)))
        }
        _ => None,
    }
}

impl Flank {
    fn for_leaf(g: &mut Graph, leaf: (NodeId, usize), term: &Leaves) -> Vec<Self> {
        let Some((axes, radius)) = circle(g, leaf.0, leaf.1) else { return vec![] };
        term.iter()
            .filter_map(|&(id, side)| {
                let p = plane(g, id, axes, &mut 32)?;
                let normal = if side == 0 { p.map(|v| g.neg(v)) } else { p };
                Some(Self { axes, radius, normal })
            })
            .collect()
    }

    // At the first hit of the closed disk, p.v <= 0. If n.v >= 0,
    // n.p can be at most R |n x v| / |v| there (project n perpendicular
    // to v). The starting half-plane stays strict while descending.
    // If its threshold is at least that maximum, no such hit exists.
    fn clear(self, g: &mut Graph, rates: &Rates) -> Option<NodeId> {
        let affine = |id: &NodeId| rates[id]?.affine(g);
        let [x, y, c] = self.normal;
        let (u, v) = (affine(&self.axes[0])?, affine(&self.axes[1])?);
        let mut mul = |a, b| g.bin_real(BinOp::Mul, a, b);
        let (xu, yv, xv, yu, uu, vv) = (mul(x, u), mul(y, v), mul(x, v), mul(y, u), mul(u, u), mul(v, v));
        let drift = g.bin_real(BinOp::Add, xu, yv);
        let cross = g.bin_real(BinOp::Sub, xv, yu);
        let speed = g.bin_real(BinOp::Add, uu, vv);
        let cc = g.un_real(UnOp::Square, c);
        let lhs = g.bin_real(BinOp::Mul, cc, speed);
        let rhs = g.bin_real(BinOp::Mul, self.radius, cross);
        let rhs = g.un_real(UnOp::Square, rhs);
        let margin = g.bin_real(BinOp::Sub, lhs, rhs);
        let minus_c = g.neg(c);
        let ok = g.min(drift, minus_c);
        let ok = g.min(ok, margin);
        let no = g.constant(1.0);
        let yes = g.constant(-1.0);
        let moving = g.neg(speed);
        let result = g.select(ok, no, yes);
        Some(g.select(moving, result, no))
    }
}

// Enumerate the ways the row point can be strictly inside. At a pointed
// roof the left and right sides can require different probe directions
// even in an arbitrarily small cell containing the apex.
type Leaves = Vec<(NodeId, usize)>;

fn start_terms(g: &Graph, id: NodeId, side: usize, budget: &mut usize) -> Option<Vec<Leaves>> {
    *budget = budget.checked_sub(1)?;
    match g.nodes[id as usize] {
        Node::Real(a) | Node::Tag(_, a) | Node::Lipschitz(a, ..) => start_terms(g, a, side, budget),
        Node::Un(UnOp::Neg, a) => start_terms(g, a, 1 - side, budget),
        Node::Bin(op @ (BinOp::Max | BinOp::Min | BinOp::Emax), a, b) => {
            let mut left = start_terms(g, a, side, budget)?;
            let right = start_terms(g, b, side, budget)?;
            if (op != BinOp::Min) == (side == 0) {
                *budget = budget.checked_sub(left.len().checked_mul(right.len())?)?;
                Some(left.iter().flat_map(|l| right.iter().map(move |r| [l.as_slice(), r].concat())).collect())
            } else {
                left.extend(right);
                Some(left)
            }
        }
        _ => Some(vec![vec![(id, side)]]),
    }
}

pub(super) struct MaterialTerm {
    pub(super) active: NodeId,
    leaves: Vec<((NodeId, usize), Vec<Flank>)>,
}

impl MaterialTerm {
    pub(super) fn of(g: &mut Graph, field: NodeId) -> Vec<Self> {
        // Distributing a branch cut repeats its fork in each term.
        // Keep such fields intact so each copy is checked once.
        if !Tape::new(g, field).jumps.is_empty() {
            return Vec::new();
        }
        let terms = start_terms(g, field, 0, &mut 384).unwrap_or_default();
        terms
            .into_iter()
            .map(|term| {
                let mut active = None;
                let mut leaves = Vec::new();
                for &(id, side) in &term {
                    let signed = if side == 0 { id } else { g.neg(id) };
                    active = Some(active.map_or(signed, |prev| g.max(prev, signed)));
                    leaves.push(((id, side), Flank::for_leaf(g, (id, side), &term)));
                }
                Self { active: active.expect("every term has a leaf"), leaves }
            })
            .collect()
    }

    pub(super) fn keeps_sign(&self, g: &mut Graph, ray: &HashMap<NodeId, [NodeId; 2]>, rates: &Rates) -> NodeId {
        let mut keeps = g.constant(-1.0);
        for ((id, side), flanks) in &self.leaves {
            let mut safe = ray[id][*side];
            for flank in flanks {
                if let Some(clear) = flank.clear(g, rates) {
                    safe = g.min(safe, clear);
                }
            }
            keeps = g.max(keeps, safe);
        }
        keeps
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink_core::interval::I;

    #[test]
    fn certified_half_planes_exclude_every_inward_crossing() {
        let mut g = Graph::new();
        let axes = [g.var(Var::X), g.var(Var::Z)];
        let radius = 3;
        let r = g.constant(f64::from(radius));
        let yes = g.constant(-1.0);
        let mut accepted = 0;
        for (u, v) in [(1, 0), (-1, 0), (0, -1), (1, -1), (-2, 1)] {
            let rates: Rates = axes
                .into_iter()
                .zip([u, v])
                .map(|(axis, s)| {
                    let rate = g.constant(f64::from(s));
                    (axis, Some(Rate { range: [rate; 2], guard: yes }))
                })
                .collect();
            for [nx, ny] in [[1, 0], [-1, 1], [1, 2], [2, -1]] {
                for offset in -9..=1 {
                    let normal = [nx, ny, offset].map(|v| g.constant(f64::from(v)));
                    let flag = Flank { axes, radius: r, normal }.clear(&mut g, &rates).unwrap();
                    let bound = eval::eval_interval(&Tape::new(&g, flag), &Ctx::new(0.0, vec![]), [I::point(0.0); 3], None);
                    if bound.hi >= 0.0 {
                        continue;
                    }
                    for x in -6..=6 {
                        for y in -6..=6 {
                            let outside = x * x + y * y - radius * radius;
                            if outside <= 0 || nx * x + ny * y + offset <= 0 {
                                continue;
                            }
                            accepted += 1;
                            // Independent integer oracle. A ray from outside hits
                            // the disk iff its quadratic has a nonnegative
                            // discriminant and the closest point is ahead.
                            let dot = x * u + y * v;
                            assert!(dot >= 0 || dot * dot < (u * u + v * v) * outside, "normal {normal:?}, start ({x},{y}), direction ({u},{v})");
                        }
                    }
                }
            }
        }
        assert!(accepted > 0);
    }
}
