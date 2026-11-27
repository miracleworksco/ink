//! A bridge made of open axis-aligned boxes can be proved at its actual
//! underside, even at the outer edges where interval cells straddle air.
//! Keep the bounds as real expressions, including the distinction between
//! a positive overlap and a zero seam.
use super::*;

type Slab = [Option<NodeId>; 2];
pub(super) type Brick = [Slab; 3];

fn intersect(g: &mut Graph, a: Brick, b: Brick) -> Brick {
    std::array::from_fn(|k| {
        std::array::from_fn(|i| match (a[k][i], b[k][i]) {
            (Some(x), Some(y)) => Some(if i == 0 { g.max(x, y) } else { g.min(x, y) }),
            (x, y) => x.or(y),
        })
    })
}

/// Decompose a strict sublevel (or superlevel) set, never its closure.
/// Bound work on shared and distributed unions. Unhandled shapes use the
/// general segment proof.
pub(super) fn bricks(g: &mut Graph, id: NodeId, level: NodeId, less: bool, budget: &mut usize) -> Option<Vec<Brick>> {
    *budget = budget.checked_sub(1)?;
    let n = g.nodes[id as usize];
    match n {
        Node::Tag(_, a) | Node::Real(a) | Node::Lipschitz(a, ..) => bricks(g, a, level, less, budget),
        Node::Var(v @ (Var::X | Var::Y | Var::Z)) => {
            let k = match v {
                Var::X => 0,
                Var::Y => 1,
                _ => 2,
            };
            let mut b = [[None; 2]; 3];
            b[k][usize::from(less)] = Some(g.real(level));
            Some(vec![b])
        }
        Node::Un(UnOp::Neg, a) => {
            let level = g.neg(level);
            bricks(g, a, level, !less, budget)
        }
        Node::Un(UnOp::Abs, a) => {
            let hi = bricks(g, a, level, less, budget)?;
            let minus = g.neg(level);
            let lo = bricks(g, a, minus, !less, budget)?;
            combine(g, hi, lo, less, budget)
        }
        Node::Bin(op @ (BinOp::Add | BinOp::Sub), a, b) => {
            let level = g.real(level);
            if !g.depends_on(b, spatial) {
                let l = if op == BinOp::Add { g.sub(level, b) } else { g.add(level, b) };
                bricks(g, a, l, less, budget)
            } else if !g.depends_on(a, spatial) {
                let l = if op == BinOp::Add { g.sub(level, a) } else { g.sub(a, level) };
                bricks(g, b, l, if op == BinOp::Add { less } else { !less }, budget)
            } else {
                None
            }
        }
        Node::Bin(op @ (BinOp::Min | BinOp::Max | BinOp::Emax), a, b) if op != BinOp::Emax || g.as_const(level) == Some(0.0) => {
            let left = bricks(g, a, level, less, budget)?;
            let right = bricks(g, b, level, less, budget)?;
            combine(g, left, right, (op != BinOp::Min) == less, budget)
        }
        _ => None,
    }
}

fn combine(g: &mut Graph, mut a: Vec<Brick>, b: Vec<Brick>, both: bool, budget: &mut usize) -> Option<Vec<Brick>> {
    if !both {
        a.extend(b);
        return Some(a);
    }
    *budget = budget.checked_sub(a.len().checked_mul(b.len())?)?;
    let mut out = Vec::new();
    for x in a {
        for &y in &b {
            out.push(intersect(g, x, y));
        }
    }
    Some(out)
}

/// Negative if a < b (or a <= b), positive otherwise. Equal boundaries
/// suffice in x/y because the source footprint is an *open* rectangle.
fn ordered(g: &mut Graph, a: NodeId, b: NodeId, strict: bool) -> NodeId {
    let yes = g.constant(-1.0);
    let no = g.constant(1.0);
    let a = g.real(a);
    if strict {
        let d = g.sub(a, b);
        g.select(d, yes, no)
    } else {
        let d = g.sub(b, a);
        g.select(d, no, yes)
    }
}

fn covers(g: &mut Graph, s: [NodeId; 2], b: [NodeId; 2]) -> NodeId {
    let a = ordered(g, s[0], b[0], false);
    let b = ordered(g, b[1], s[1], false);
    g.max(a, b)
}

pub(super) fn unsupported(g: &mut Graph, f: NodeId, bed: NodeId, depth: NodeId, height: NodeId, limit: NodeId, step: NodeId) -> Option<NodeId> {
    if [bed, depth, height, limit, step].iter().any(|&n| g.depends_on(n, spatial)) {
        return None;
    }
    let zero = g.constant(0.0);
    let mut list = Vec::new();
    for b in bricks(g, f, zero, true, &mut 512)? {
        let [[Some(x0), Some(x1)], [Some(y0), Some(y1)], [Some(z0), Some(z1)]] = b else { return None };
        list.push([[x0, x1], [y0, y1], [z0, z1]]);
    }
    // Pairing each bottom with every potential bearing is quadratic.
    if list.len() > 24 {
        return None;
    }
    // Both failure probes read rows strictly above bed + depth when
    // height >= depth. A bearing starting at the bed then covers every
    // required column bottom, even under a roof closer than depth to it.
    let at_h = g.add(bed, height);
    let at_depth = g.add(bed, depth);
    let first_row = g.min(at_h, at_depth);
    // The model can fold bed + depth in floats. Use the actual lean
    // threshold above and compare it to the real column length here.
    let needed = g.bin_real(BinOp::Add, bed, depth);
    let full_depth_above_bed = ordered(g, needed, first_row, false);
    let width = g.real(step);
    let two_steps = g.add(width, width);
    let room = g.sub(limit, two_steps);
    let mut all = g.constant(-1.0);
    for b in &list {
        let mut safe = ordered(g, b[2][0], bed, false);
        for &[lo, hi] in b {
            let empty = ordered(g, hi, lo, false);
            safe = g.min(safe, empty);
        }
        let bottom = g.sub(b[2][0], depth);
        // Only material extending strictly through the bottom face and
        // below the full test depth can carry this roof.
        let bearings: Vec<_> = list
            .iter()
            .map(|s| {
                let below = ordered(g, s[2][0], bottom, true);
                let on_bed = ordered(g, s[2][0], bed, false);
                let grounded = g.max(full_depth_above_bed, on_bed);
                let below = g.min(below, grounded);
                let above = ordered(g, b[2][0], s[2][1], true);
                g.max(below, above)
            })
            .collect();
        for (s, &bearing) in list.iter().zip(&bearings) {
            let x = covers(g, s[0], b[0]);
            let y = covers(g, s[1], b[1]);
            let footprint = g.max(x, y);
            let direct = g.max(bearing, footprint);
            safe = g.min(safe, direct);
        }
        // A bearing at either end must cover the full transverse extent
        // of this open roof. Their longest reaches may be chosen
        // independently. Each end needs more than one sampling step of
        // overlap, so even the least favourable row has a bearing there.
        for axis in 0..2 {
            let [lo, hi] = b[axis];
            let (mut left, mut right) = (lo, hi);
            for (s, &bearing) in list.iter().zip(&bearings) {
                let across = covers(g, s[1 - axis], b[1 - axis]);
                let full_width = g.max(bearing, across);
                let [start, end] = s[axis];
                let at_left = ordered(g, start, lo, false);
                let at_left = g.max(full_width, at_left);
                let reaches = g.select(at_left, end, lo);
                left = g.max(left, reaches);
                let at_right = ordered(g, hi, end, false);
                let at_right = g.max(full_width, at_right);
                let reaches = g.select(at_right, start, hi);
                right = g.min(right, reaches);
            }
            let left_width = g.sub(left, lo);
            let right_width = g.sub(hi, right);
            let left_ok = ordered(g, width, left_width, true);
            let right_ok = ordered(g, width, right_width, true);
            let pair = g.max(left_ok, right_ok);
            // A grid anchor can sit up to w past each edge of the gap.
            let gap = g.sub(right, left);
            let span = ordered(g, gap, room, false);
            let bridge = g.max(pair, span);
            safe = g.min(safe, bridge);
        }
        all = g.max(all, safe);
    }
    let positive = ordered(g, zero, width, true);
    let all = g.max(all, positive);
    Some(g.neg(all))
}
