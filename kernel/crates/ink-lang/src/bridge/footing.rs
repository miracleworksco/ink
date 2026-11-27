//! A rectangular member can descend through a bearing with a different
//! expression. Certify its entire open footprint from the bed right up to
//! its bottom. Both solids then approach every point of the joint.
use super::*;

// The ordinary bounds enclose the closed box. For each endpoint also
// track a sufficient condition for its value to be unattainable in the
// open box. Then an upper bound of zero proves strict material only when
// that endpoint is open. An unhandled operation keeps both endpoints.
fn interior(g: &mut Graph, f: NodeId, axes: &[(Var, Range)]) -> NodeId {
    let mut bounds = Bounds::default();
    let yes = g.constant(-1.0);
    let no = g.constant(1.0);
    let mut open = HashMap::<NodeId, Range>::new();
    let mut sign = HashMap::<NodeId, Range>::new();
    for id in g.reachable(&[f]) {
        let n = g.nodes[id as usize];
        bounds.insert(g, id, axes);
        let strict = match n {
            Node::Var(v) if spatial(v) => [yes; 2],
            Node::Lipschitz(a, ..) | Node::Tag(_, a) | Node::Real(a) | Node::Read(_, _, a, _) | Node::Deriv(a, ..) => open[&a],
            Node::Un(UnOp::Neg, a) => [open[&a][1], open[&a][0]],
            // A sum cannot attain its bound if either contributing
            // endpoint is open. A min or max may pick either child.
            Node::Bin(op @ (BinOp::Add | BinOp::Sub), a, b) => [0, 1].map(|i| g.min(open[&a][i], open[&b][if op == BinOp::Add { i } else { 1 - i }])),
            Node::Bin(BinOp::Max | BinOp::Min | BinOp::Emax, a, b) => [0, 1].map(|i| g.max(open[&a][i], open[&b][i])),
            Node::Un(UnOp::Abs, a) if let Some([l, h]) = bounds.spans[&a] => {
                let neg = g.neg(h);
                let left = g.select(neg, no, open[&a][1]);
                let right = g.select(l, no, open[&a][0]);
                [g.min(left, right), g.max(open[&a][0], open[&a][1])]
            }
            _ => [no; 2],
        };
        open.insert(id, strict);
        let safe = [0, 1].map(|i| {
            let mut ok = no;
            if let Some([lo, hi]) = bounds.spans[&id] {
                let v = if i == 0 { hi } else { g.neg(lo) };
                let away = g.neg(v);
                let at = g.select(away, no, strict[1 - i]);
                ok = g.select(v, yes, at);
            }
            let csg = match n {
                Node::Lipschitz(a, ..) | Node::Real(a) | Node::Tag(_, a) => sign[&a][i],
                Node::Un(UnOp::Neg, a) => sign[&a][1 - i],
                Node::Bin(op @ (BinOp::Min | BinOp::Max | BinOp::Emax), a, b) => {
                    if (op != BinOp::Min) == (i == 0) {
                        g.max(sign[&a][i], sign[&b][i])
                    } else {
                        g.min(sign[&a][i], sign[&b][i])
                    }
                }
                _ => no,
            };
            g.min(ok, csg)
        });
        sign.insert(id, safe);
    }
    sign[&f][0]
}

pub(super) fn members(g: &mut Graph, f: NodeId, bed: NodeId) -> Vec<(NodeId, NodeId)> {
    let mut parts = vec![];
    if g.depends_on(bed, spatial) || joints::members(g, f, &mut parts, &mut 40).is_none() || parts.len() < 2 {
        return vec![];
    }
    let zero = g.constant(0.0);
    let yes = g.constant(-1.0);
    let no = g.constant(1.0);
    let mut any = false;
    let results = parts
        .iter()
        .map(|&p| {
            let Some(bricks) = boxes::bricks(g, p, zero, true, &mut 192) else { return (p, no) };
            if bricks.len() > 6 {
                return (p, no);
            }
            let mut all = yes;
            for brick in bricks {
                let [[Some(x0), Some(x1)], [Some(y0), Some(y1)], [Some(bottom), _]] = brick else {
                    return (p, no);
                };
                let height = g.bin_real(BinOp::Sub, bed, bottom);
                let axes = [(Var::X, [x0, x1]), (Var::Y, [y0, y1]), (Var::Z, [bed, bottom])];
                let mut carried = no;
                for &other in &parts {
                    if p != other {
                        let full = interior(g, other, &axes);
                        carried = g.min(carried, full);
                    }
                }
                let ok = g.select(height, carried, yes);
                all = g.max(all, ok);
            }
            any |= all != no;
            (p, all)
        })
        .collect();
    if any { results } else { vec![] }
}
