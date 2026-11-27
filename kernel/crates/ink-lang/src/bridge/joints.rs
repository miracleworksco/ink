//! A field for the material used by the printing proof. Its strict interior
//! can additionally contain a proven internal face. Never substitute this
//! sign-only representation into distance arithmetic or a complement.
use super::*;

#[derive(Clone, Copy, PartialEq, Eq)]
struct Cut {
    coordinate: NodeId,
    level: NodeId,
    upper: bool,
}

#[derive(Clone, Copy)]
struct Wall {
    field: NodeId,
    cut: Option<Cut>,
}

// The level set of u +/- v has no volume when u is a spatial
// coordinate absent from v. Moving along that axis crosses the level
// exactly once. Merely sharing an expression is insufficient: it could
// be constant on a whole region (for example max(x, 0)).
fn crossing_coordinate(g: &Graph, id: NodeId) -> bool {
    let axis = |i: NodeId| if let Node::Var(v) = g.nodes[i as usize] { spatial(v).then_some(v) } else { None };
    match g.nodes[id as usize] {
        Node::Var(v) => spatial(v),
        Node::Bin(BinOp::Add | BinOp::Sub, a, b) => [(a, b), (b, a)].iter().any(|&(u, v)| axis(u).is_some_and(|k| !g.value_depends_on(v, |x| x == k))),
        _ => false,
    }
}

// Peel only exact real operations from a spatial inequality. In particular
// centre - half_width must not round to a neighbouring slab's endpoint.
fn walls(g: &mut Graph, id: NodeId, level: NodeId, upper: bool, out: &mut Vec<Wall>, budget: &mut usize) -> Option<()> {
    *budget = budget.checked_sub(1)?;
    match g.nodes[id as usize] {
        Node::Real(a) | Node::Tag(_, a) | Node::Lipschitz(a, ..) => walls(g, a, level, upper, out, budget),
        Node::Bin(BinOp::Max | BinOp::Emax, a, b) if upper && g.as_const(level) == Some(0.0) => {
            walls(g, a, level, true, out, budget)?;
            walls(g, b, level, true, out, budget)
        }
        Node::Un(UnOp::Neg, a) => {
            let level = g.un_real(UnOp::Neg, level);
            walls(g, a, level, !upper, out, budget)
        }
        Node::Un(UnOp::Abs, a) if upper => {
            walls(g, a, level, true, out, budget)?;
            let level = g.un_real(UnOp::Neg, level);
            walls(g, a, level, false, out, budget)
        }
        Node::Bin(op @ (BinOp::Sub | BinOp::Add), a, b) if !g.depends_on(b, spatial) && g.depends_on(a, spatial) => {
            let inverse = if op == BinOp::Add { BinOp::Sub } else { BinOp::Add };
            let level = g.bin_real(inverse, level, b);
            walls(g, a, level, upper, out, budget)
        }
        Node::Bin(op @ (BinOp::Sub | BinOp::Add), a, b) if !g.depends_on(a, spatial) && g.depends_on(b, spatial) => {
            let (level, upper) = if op == BinOp::Add { (g.bin_real(BinOp::Sub, level, a), upper) } else { (g.bin_real(BinOp::Sub, a, level), !upper) };
            walls(g, b, level, upper, out, budget)
        }
        _ => {
            let cut = crossing_coordinate(g, id).then_some(Cut { coordinate: id, level, upper });
            let field = if upper { g.bin_real(BinOp::Sub, id, level) } else { g.bin_real(BinOp::Sub, level, id) };
            out.push(Wall { field, cut });
            Some(())
        }
    }
}

fn ordered(g: &Graph, a: NodeId, b: NodeId) -> bool {
    a == b || g.as_const(a).zip(g.as_const(b)).is_some_and(|(a, b)| a <= b)
}

// A narrower parallel half-space implies a wider one.
fn implies(g: &Graph, a: Wall, b: Wall) -> bool {
    a.field == b.field
        || a.cut
            .zip(b.cut)
            .is_some_and(|(a, b)| a.coordinate == b.coordinate && a.upper == b.upper && if a.upper { ordered(g, a.level, b.level) } else { ordered(g, b.level, a.level) })
}

// From (A AND x<h) OR (B AND x>h), A AND B is an additional
// description of interior material. Off the plane, at least one member
// contains it. On the plane, strict A and B hold in a neighbourhood.
// Keep every transverse constraint: an edge or a corner isn't a bearing.
fn join(g: &mut Graph, a: &[Wall], b: &[Wall]) -> Option<(NodeId, bool, bool)> {
    let (i, j) = a.iter().enumerate().find_map(|(i, w)| {
        let p = w.cut?;
        b.iter().position(|v| v.cut.is_some_and(|q| p.coordinate == q.coordinate && p.upper != q.upper && p.level == q.level)).map(|j| (i, j))
    })?;
    let mut rest = Vec::<Wall>::new();
    for &w in a[..i].iter().chain(&a[i + 1..]).chain(&b[..j]).chain(&b[j + 1..]) {
        if !rest.iter().any(|&v| implies(g, v, w)) {
            rest.retain(|&v| !implies(g, w, v));
            rest.push(w);
        }
    }
    let covered = |source: &[Wall]| rest.iter().all(|&w| source.iter().any(|&v| implies(g, v, w)));
    let (left, right) = (covered(a), covered(b));
    // When extending just one member, end its extension inside the
    // neighbour. Ending it on the neighbour's far exterior face would
    // introduce another coincident edge into the path search.
    if left != right {
        let face = if left { a[i] } else { b[j] }.cut.unwrap();
        if let Some(far) = rest.iter_mut().find(|w| w.cut.is_some_and(|c| c.coordinate == face.coordinate && c.upper == face.upper)) {
            let c = far.cut.unwrap();
            let half = g.constant(0.5);
            let a = g.bin_real(BinOp::Mul, half, face.level);
            let b = g.bin_real(BinOp::Mul, half, c.level);
            let mid = g.bin_real(BinOp::Add, a, b);
            let level = if face.upper { g.min(mid, c.level) } else { g.max(mid, c.level) };
            let axis = face.coordinate;
            far.field = if face.upper { g.bin_real(BinOp::Sub, axis, level) } else { g.bin_real(BinOp::Sub, level, axis) };
        }
    }
    let solid = rest.into_iter().map(|w| w.field).reduce(|a, b| g.max(a, b)).unwrap_or_else(|| g.constant(-1.0));
    Some((solid, left, right))
}

pub(super) fn members(g: &Graph, id: NodeId, out: &mut Vec<NodeId>, budget: &mut usize) -> Option<()> {
    *budget = budget.checked_sub(1)?;
    match g.nodes[id as usize] {
        Node::Tag(_, a) | Node::Lipschitz(a, ..) => members(g, a, out, budget),
        Node::Bin(BinOp::Min, a, b) => {
            members(g, a, out, budget)?;
            members(g, b, out, budget)
        }
        _ => {
            out.push(id);
            Some(())
        }
    }
}

fn rewrite(g: &mut Graph, id: NodeId, memo: &mut HashMap<NodeId, NodeId>) -> NodeId {
    if let Some(&v) = memo.get(&id) {
        return v;
    }
    let new = match g.nodes[id as usize] {
        Node::Tag(_, a) | Node::Lipschitz(a, ..) => {
            let x = rewrite(g, a, memo);
            if x == a { id } else { x }
        }
        Node::Bin(BinOp::Max | BinOp::Emax, a, b) => {
            let (x, y) = (rewrite(g, a, memo), rewrite(g, b, memo));
            if (x, y) == (a, b) { id } else { g.max(x, y) }
        }
        Node::Bin(BinOp::Min, ..) => {
            let mut parts = vec![];
            if members(g, id, &mut parts, &mut 40).is_none() {
                return id;
            }
            let mut changed = false;
            for x in &mut parts {
                let rebuilt = rewrite(g, *x, memo);
                changed |= *x != rebuilt;
                *x = rebuilt;
            }
            let zero = g.constant(0.0);
            let sides: Vec<_> = parts
                .iter()
                .map(|&p| {
                    let mut out = Vec::new();
                    walls(g, p, zero, true, &mut out, &mut 96).map(|()| out)
                })
                .collect();
            let mut keep = vec![true; parts.len()];
            let mut extra = vec![];
            for (j, right) in sides.iter().enumerate() {
                let Some(b) = right else {
                    continue;
                };
                for (i, left) in sides[..j].iter().enumerate() {
                    if let Some(a) = left
                        && let Some((w, ca, cb)) = join(g, a, b)
                    {
                        keep[i] &= !ca;
                        keep[j] &= !cb;
                        extra.push(w);
                        changed = true;
                    }
                }
            }
            if changed { parts.into_iter().zip(keep).filter_map(|(p, k)| k.then_some(p)).chain(extra).reduce(|a, b| g.min(a, b)).unwrap() } else { id }
        }
        _ => id,
    };
    memo.insert(id, new);
    new
}

pub(super) fn material(g: &mut Graph, f: NodeId) -> NodeId {
    // Only the spatial variables vary in the neighbourhood of a seam.
    // A parameter expression (for example an angle converted to radians)
    // is one fixed number throughout that neighbourhood.
    for id in g.reachable(&[f]) {
        match g.nodes[id as usize] {
            Node::Var(_)
            | Node::Const(_)
            | Node::Within(..)
            | Node::Frame(_)
            | Node::Tag(..)
            | Node::Lipschitz(..)
            | Node::Real(_)
            | Node::Un(UnOp::Neg | UnOp::Abs | UnOp::Square | UnOp::Sqrt | UnOp::Rad | UnOp::Deg | UnOp::Cos | UnOp::Sin, _)
            | Node::Bin(BinOp::Add | BinOp::Sub | BinOp::Mul | BinOp::Min | BinOp::Max | BinOp::Emax, ..) => {}
            _ if !g.depends_on(id, spatial) => {}
            _ => return f,
        }
    }
    rewrite(g, f, &mut HashMap::new())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_face_requires_material_on_both_sides_for_each_parameter_value() {
        let mut g = Graph::new();
        let z = g.var(Var::Z);
        let x = g.var(Var::X);
        let h = g.var(Var::Param(0));
        let lower_end = g.var(Var::Param(1));
        let upper_end = g.var(Var::Param(2));
        let zero = g.constant(0.0);
        let one = g.constant(1.0);
        let left_face = g.sub(z, h);
        let right_face = g.sub(h, z);
        let bottom = g.sub(lower_end, z);
        let top = g.sub(z, upper_end);
        let width = g.sub(x, one);
        let lower = g.max(bottom, left_face);
        let lower = g.max(lower, width);
        let upper = g.max(right_face, top);
        let side = g.sub(x, zero);
        let upper = g.max(upper, side);
        let union = g.min(lower, upper);
        let joined = material(&mut g, union);
        assert_ne!(union, joined);
        let tape = Tape::new(&g, joined);
        // Independent open-rectangle oracle, including empty and reversed
        // slabs. All inputs and intermediates are exact dyadic numbers.
        for seam in [-1.0, 0.0, 1.0] {
            for bottom in [-2.0, 0.0, 2.0] {
                for top in [-2.0, 0.0, 2.0] {
                    let ctx = Ctx::new(0.0, vec![seam, bottom, top]);
                    for i in -5..=5 {
                        let height = f64::from(i) / 2.0;
                        for side in [-1.0, 0.0, 0.5, 1.0] {
                            let left = bottom < height && height < seam && side < 1.0;
                            let right = seam < height && height < top && side < 0.0;
                            let face = height == seam && bottom < seam && seam < top && side < 0.0;
                            let value = eval::eval_point(&tape, &ctx, [side, 0.0, height]);
                            assert_eq!(value < 0.0, left || right || face, "bottom {bottom}, seam {seam}, top {top}, x {side}, z {height}");
                        }
                    }
                }
            }
        }
    }
}
