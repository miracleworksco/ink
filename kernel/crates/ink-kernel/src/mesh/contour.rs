use super::brick::{B, Brick, vidx};
use crate::octree::{Grid, Local};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_points};
use ink_core::interval::I;
use ink_core::sign::{drawn, inside};
use ink_core::tape::Tape;
use std::collections::HashMap;

/// Cell corner `c` has offset bits (x = 1, y = 2, z = 4). The 12 cell edges
/// join corners differing in one bit.
pub(super) const EDGES: [(usize, usize); 12] = {
    let mut e = [(0, 0); 12];
    let mut n = 0;
    let mut c = 0;
    while c < 8 {
        let mut bit = 1;
        while bit < 8 {
            if c & bit == 0 {
                e[n] = (c, c | bit);
                n += 1;
            }
            bit <<= 1;
        }
        c += 1;
    }
    e
};

fn edge_id(a: usize, b: usize) -> usize {
    let (a, b) = (a.min(b), a.max(b));
    EDGES.iter().position(|&e| e == (a, b)).unwrap()
}

/// The six faces as rings of four corners, and the four ring edges.
pub(super) struct Face {
    axis: usize,
    pub(super) ring: [usize; 4],
    edges: [usize; 4],
}

pub(super) fn faces() -> &'static [Face; 6] {
    static F: std::sync::OnceLock<[Face; 6]> = std::sync::OnceLock::new();
    F.get_or_init(|| {
        std::array::from_fn(|f| {
            let (axis, side) = (f / 2, f % 2);
            let (a, u, v) = (1 << axis, 1 << ((axis + 1) % 3), 1 << ((axis + 2) % 3));
            let base = if side == 1 { a } else { 0 };
            let ring = [base, base | u, base | u | v, base | v];
            let edges = std::array::from_fn(|k| edge_id(ring[k], ring[(k + 1) % 4]));
            Face { axis, ring, edges }
        })
    })
}

/// Group a cell's crossing edges into patches (bit masks over the 12 edges).
pub(super) fn patches(inside_c: &[bool; 8], face_center_inside: impl Fn(usize) -> bool) -> Vec<u16> {
    let mut parent: [usize; 12] = std::array::from_fn(|i| i);
    fn find(p: &mut [usize; 12], mut x: usize) -> usize {
        while p[x] != x {
            x = p[x];
        }
        x
    }
    let union = |p: &mut [usize; 12], a: usize, b: usize| {
        let (ra, rb) = (find(p, a), find(p, b));
        p[ra] = rb;
    };
    for (fi, f) in faces().iter().enumerate() {
        let r = f.ring.map(|c| inside_c[c]);
        let crossing: Vec<usize> = (0..4).filter(|&k| r[k] != r[(k + 1) % 4]).collect();
        match crossing.len() {
            2 => union(&mut parent, f.edges[crossing[0]], f.edges[crossing[1]]),
            4 => {
                // Saddle: if the centre matches ring corner 0, corners 0 and 2
                // connect through it and the segments cut off corners 1 and 3.
                if face_center_inside(fi) == r[0] {
                    union(&mut parent, f.edges[0], f.edges[1]);
                    union(&mut parent, f.edges[2], f.edges[3]);
                } else {
                    union(&mut parent, f.edges[3], f.edges[0]);
                    union(&mut parent, f.edges[1], f.edges[2]);
                }
            }
            _ => {}
        }
    }
    let mut out: Vec<(usize, u16)> = Vec::new();
    for (e, &(a, b)) in EDGES.iter().enumerate() {
        if inside_c[a] == inside_c[b] {
            continue;
        }
        let root = find(&mut parent, e);
        match out.iter_mut().find(|(r, _)| *r == root) {
            Some((_, m)) => *m |= 1 << e,
            None => out.push((root, 1 << e)),
        }
    }
    out.into_iter().map(|(_, m)| m).collect()
}

/// Canonical key of a cell face: (global min corner of the face, normal axis).
pub(super) fn face_key(cell: [i64; 3], f: &Face) -> ([i64; 3], usize) {
    let mut g = cell;
    g[f.axis] += ((f.ring[0] >> f.axis) & 1) as i64;
    (g, f.axis)
}

/// Face centre, computed only from the canonical key so every cell sharing the
/// face gets bit-identical coordinates (and therefore the same field value).
pub(super) fn face_center(grid: &Grid, (g, axis): ([i64; 3], usize)) -> V3 {
    std::array::from_fn(|a| {
        let i = g[a] as f64 + if a == axis { 0.0 } else { 0.5 };
        grid.origin[a] + grid.h * i
    })
}

pub(super) const REFINE_STEPS: usize = 6;

/// Move each edge crossing onto the surface: batched false-position steps,
/// safeguarded to shrink the bracket every iteration. Linear interpolation
/// alone lands off the surface near sharp features, where distance fields curve.
/// A step sample on 0 is read as every sample is (`drawn`, `within` where
/// given), and those undecided are returned with the crossings, for the
/// mesh's report.
pub(super) fn refine_crossings(tape: &Tape, within: Option<([I; 3], &Tape)>, ctx: &Ctx, br: &mut [(V3, f64, V3, f64)]) -> (Vec<V3>, Vec<V3>) {
    let mut open = Vec::new();
    let lerp = |(pa, va, pb, vb): (V3, f64, V3, f64)| {
        let t = if va == vb { 0.5 } else { (va / (va - vb)).clamp(0.1, 0.9) };
        [pa[0] + (pb[0] - pa[0]) * t, pa[1] + (pb[1] - pa[1]) * t, pa[2] + (pb[2] - pa[2]) * t]
    };
    let mut mid: Vec<V3> = Vec::with_capacity(br.len());
    let mut v = Vec::new();
    for _ in 0..REFINE_STEPS {
        mid.clear();
        mid.extend(br.iter().map(|b| lerp(*b)));
        eval_points(tape, ctx, &mid, &mut v);
        open.extend(drawn(tape, ctx, &mid, &mut v, within));
        for ((b, p), &f) in br.iter_mut().zip(&mid).zip(&v) {
            if inside(f) == inside(b.1) {
                b.0 = *p;
                b.1 = f;
            } else {
                b.2 = *p;
                b.3 = f;
            }
        }
    }
    let cross = br
        .iter()
        .map(|&(pa, va, pb, vb)| {
            let t = if va == vb { 0.5 } else { (va / (va - vb)).clamp(0.0, 1.0) };
            [pa[0] + (pb[0] - pa[0]) * t, pa[1] + (pb[1] - pa[1]) * t, pa[2] + (pb[2] - pa[2]) * t]
        })
        .collect();
    (cross, open)
}

/// Minimize sum (n . (v - p))^2 + lambda |v - mass|^2, clamped near the cell;
/// and the mass point, clamped into the cell (`untangle`).
pub(super) fn qef(pts: &[V3], grads: &[[f64; 4]], cmin: V3, h: f64) -> (V3, V3) {
    let n = pts.len() as f64;
    let mut mass = [0.0; 3];
    for p in pts {
        for a in 0..3 {
            mass[a] += p[a] / n;
        }
    }
    let lambda = 0.05;
    let mut m = [[0.0; 3]; 3];
    let mut rhs = [0.0; 3];
    for a in 0..3 {
        m[a][a] = lambda;
    }
    for (p, g) in pts.iter().zip(grads) {
        let len = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
        if !(len > 0.0 && len.is_finite()) {
            continue;
        }
        let nn = [g[1] / len, g[2] / len, g[3] / len];
        let d = nn[0] * (p[0] - mass[0]) + nn[1] * (p[1] - mass[1]) + nn[2] * (p[2] - mass[2]);
        for a in 0..3 {
            for b in 0..3 {
                m[a][b] += nn[a] * nn[b];
            }
            rhs[a] += nn[a] * d;
        }
    }
    let x = solve3(m, rhs);
    // Features can sit just outside the cell whose edges detected them; allow
    // half a cell of slack so sharp edges and corners are reproduced exactly.
    let m = 0.5 * h;
    (std::array::from_fn(|a| (mass[a] + x[a]).clamp(cmin[a] - m, cmin[a] + h + m)), std::array::from_fn(|a| mass[a].clamp(cmin[a], cmin[a] + h)))
}

/// m x = b by Cramer's rule. m is lambda I plus outer products of unit
/// normals: positive definite, its determinant at least lambda^3.
fn solve3(m: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let det = |m: &[[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(&m);
    std::array::from_fn(|c| {
        let mut mc = m;
        for r in 0..3 {
            mc[r][c] = b[r];
        }
        det(&mc) / d
    })
}

/// Split each quad along the diagonal whose midpoint lies closer to the true
/// surface. At sharp creases the wrong diagonal cuts a notch into the edge.
/// Ties within floating-point evaluation error are broken deterministically
/// by the ordered grid coordinates of each diagonal's endpoints, so one-bit
/// roundoff on a flat or symmetric face cannot flip its triangles.
/// The midpoints are evaluated with the tape simplified over each block
/// they fall in (`Local`): the same values, without the whole field's cost
/// at each. In runs, each begun only until `stopped` (None then).
pub(super) fn split_quads(local: &Local, verts: &[V3], cells: &[[i64; 3]], quads: Vec<[u32; 4]>, h: f64, stopped: &dyn Fn() -> bool) -> Option<Vec<[u32; 3]>> {
    let mid = |a: u32, b: u32| {
        let (p, q) = (verts[a as usize], verts[b as usize]);
        [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0, (p[2] + q[2]) / 2.0]
    };
    let pair = |a: u32, b: u32| {
        let (ca, cb) = (cells[a as usize], cells[b as usize]);
        if ca <= cb { (ca, cb) } else { (cb, ca) }
    };
    let mut out = Vec::with_capacity(quads.len() * 2);
    for run in quads.chunks(1 << 15) {
        if stopped() {
            return None;
        }
        let pts: Vec<V3> = run.iter().flat_map(|q| [mid(q[0], q[2]), mid(q[1], q[3])]).collect();
        let f = local.eval(&pts);
        for (k, q) in run.iter().enumerate() {
            let (a, b) = (f[2 * k].abs(), f[2 * k + 1].abs());
            let reach = pts[2 * k].iter().chain(&pts[2 * k + 1]).fold(h, |m, v| m.max(v.abs()));
            let tol = 1e-12 * a.max(b).max(reach);
            let first = if (a - b).abs() <= tol { pair(q[0], q[2]) <= pair(q[1], q[3]) } else { a < b };
            if first {
                out.push([q[0], q[1], q[2]]);
                out.push([q[0], q[2], q[3]]);
            } else {
                out.push([q[0], q[1], q[3]]);
                out.push([q[1], q[2], q[3]]);
            }
        }
    }
    Some(out)
}

/// The quads brick `b`'s sign changes make, but those whose four cells lie
/// in one block drawn whole (`block_of`, Flat); and of the others, each
/// edge joining two cells of one such block, with the block.
#[allow(clippy::type_complexity)]
pub(super) fn quads(
    b: &Brick,
    cells: &HashMap<[i64; 3], (u32, u8)>,
    owners: &[(u16, u32)],
    block_of: &(dyn Fn([i64; 3]) -> Option<u32> + Sync),
) -> (Vec<[u32; 4]>, Vec<(u32, [u32; 2])>) {
    let (mut out, mut rims) = (Vec::new(), Vec::new());
    if b.patches.is_empty() {
        return (out, rims);
    }
    // Vertex of the patch in `cell` that owns the grid edge from global corner
    // `g` along axis `d`.
    let vertex = |cell: [i64; 3], g: [i64; 3], d: usize| -> Option<u32> {
        let &(start, count) = cells.get(&cell)?;
        let c0 = ((g[0] - cell[0]) | ((g[1] - cell[1]) << 1) | ((g[2] - cell[2]) << 2)) as usize;
        let e = edge_id(c0, c0 | (1 << d));
        owners[start as usize..start as usize + count as usize].iter().find(|(m, _)| m & (1 << e) != 0).map(|&(_, v)| v)
    };
    for k in 0..B as usize {
        for j in 0..B as usize {
            for i in 0..B as usize {
                let v0 = b.vals[vidx(i, j, k)];
                let g = [b.i0[0] + i as i64, b.i0[1] + j as i64, b.i0[2] + k as i64];
                for d in 0..3 {
                    let v1 = match d {
                        0 => b.vals[vidx(i + 1, j, k)],
                        1 => b.vals[vidx(i, j + 1, k)],
                        _ => b.vals[vidx(i, j, k + 1)],
                    };
                    if inside(v0) == inside(v1) {
                        continue;
                    }
                    // (d, u, v) is a right-handed cyclic triple.
                    let (u, v) = ((d + 1) % 3, (d + 2) % 3);
                    let at = |du: i64, dv: i64| {
                        let mut c = g;
                        c[u] -= du;
                        c[v] -= dv;
                        c
                    };
                    let cs = [at(1, 1), at(0, 1), at(0, 0), at(1, 0)];
                    let (Some(q00), Some(q10), Some(q11), Some(q01)) = (vertex(cs[0], g, d), vertex(cs[1], g, d), vertex(cs[2], g, d), vertex(cs[3], g, d)) else {
                        continue;
                    };
                    let blocks = cs.map(block_of);
                    if blocks[0].is_some() && blocks.iter().all(|k| *k == blocks[0]) {
                        continue;
                    }
                    // Counter-clockwise around +d when the surface normal points along +d.
                    let (q, k) = if inside(v0) { ([q00, q10, q11, q01], [0, 1, 2, 3]) } else { ([q00, q01, q11, q10], [0, 3, 2, 1]) };
                    for e in 0..4 {
                        if let (Some(x), Some(y)) = (blocks[k[e]], blocks[k[(e + 1) % 4]])
                            && x == y
                        {
                            rims.push((x, [q[e], q[(e + 1) % 4]]));
                        }
                    }
                    out.push(q);
                }
            }
        }
    }
    (out, rims)
}
