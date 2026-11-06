//! Sparse manifold dual contouring.
//!
//! The octree finds bricks of `B^3` cells that straddle the surface. Each brick
//! evaluates its corner grid in SIMD batches. Within each sign-changing cell,
//! edge crossings are grouped into surface patches (connected through the
//! cell's faces; ambiguous faces are resolved by the field at the face centre,
//! which both neighbouring cells compute identically). Each patch gets one
//! vertex placed by a quadratic error function (QEF) over surface normals, which
//! keeps sharp edges and corners. One quad is emitted per sign-changing grid
//! edge, using the patch in each of its four cells that owns that edge. One
//! vertex per patch (not per cell) keeps thin walls and sheets manifold.
//!
//! What sampling at `h` cannot see is reported, not hidden: `mesh_checked`
//! also evaluates the centre of every cell whose corners agree in sign but
//! lie within half its diagonal of the surface (by their values, as for a
//! distance), and counts the cells whose centre disagrees: material, or a
//! gap, thinner than the cell, which the mesh leaves out.

use crate::octree::{Cell, Grid, descend, root_for};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_points};
use ink_core::tape::Tape;
use ink_core::v3::{EMPTY_BOX, box_union, cross, dot, norm, normal};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

const B: i64 = 8;
const N: usize = B as usize + 1;

#[derive(Clone, Debug, Default)]
pub struct Mesh {
    pub verts: Vec<V3>,
    pub tris: Vec<[u32; 3]>,
}

impl Mesh {
    /// The sum of its triangles' areas.
    pub fn area(&self) -> f64 {
        self.tris.iter().map(|t| 0.5 * norm(normal(t.map(|i| self.verts[i as usize])))).sum()
    }
}

struct Brick {
    i0: [i64; 3],
    vals: Vec<f64>,
    /// (cell, mask of the 12 cell edges this patch owns, vertex position)
    patches: Vec<([i64; 3], u16, V3)>,
    /// Cells whose corners agree but whose centre does not, those whose
    /// centre is in material, its thickest there, and one centre.
    missed: (usize, usize, f64, Option<V3>),
}

/// Sign changes the mesh misses: cells whose eight corners agree but whose
/// centre does not (material or a gap thinner than the cell); of them, those
/// whose centre is in material, and the thickest that material is at them;
/// and the centre of one of them. A sample: features that miss every centre
/// too go uncounted.
#[derive(Clone, Copy, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct Unresolved {
    pub cells: usize,
    pub material: usize,
    pub thickest: f64,
    pub at: Option<V3>,
}

#[inline]
fn inside(v: f64) -> bool {
    v < 0.0
}

#[inline]
fn vidx(i: usize, j: usize, k: usize) -> usize {
    ink_core::grid::index([i, j, k], [N; 3])
}

/// Cell corner `c` has offset bits (x = 1, y = 2, z = 4). The 12 cell edges
/// join corners differing in one bit.
const EDGES: [(usize, usize); 12] = {
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
struct Face {
    axis: usize,
    ring: [usize; 4],
    edges: [usize; 4],
}

fn faces() -> &'static [Face; 6] {
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
fn patches(inside_c: &[bool; 8], face_center_inside: impl Fn(usize) -> bool) -> Vec<u16> {
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

pub fn mesh(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64) -> Mesh {
    sample(tape, ctx, lo, hi, h, false).0
}

/// `mesh`, and the sign changes it misses (see the module notes).
pub fn mesh_checked(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64) -> (Mesh, Unresolved) {
    sample(tape, ctx, lo, hi, h, true)
}

fn sample(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64, check: bool) -> (Mesh, Unresolved) {
    let pad = [lo[0] - h, lo[1] - h, lo[2] - h];
    let (grid, root) = root_for(pad, [hi[0] + h, hi[1] + h, hi[2] + h], h, B);
    let bricks: Vec<Brick> = descend(tape, ctx, &grid, root, B, &|t: &Arc<Tape>, c| brick(t, tape, ctx, &grid, c, check), &|_: &Arc<Tape>, _| None);
    let unresolved = bricks.iter().fold(Unresolved::default(), |u, b| Unresolved {
        cells: u.cells + b.missed.0,
        material: u.material + b.missed.1,
        thickest: u.thickest.max(b.missed.2),
        at: u.at.or(b.missed.3),
    });

    // cell -> range of its patches in `owners`; owners hold (edge mask, vertex).
    let total: usize = bricks.iter().map(|b| b.patches.len()).sum();
    let mut cells: HashMap<[i64; 3], (u32, u8)> = HashMap::with_capacity(total);
    let mut owners: Vec<(u16, u32)> = Vec::with_capacity(total);
    let mut verts = Vec::with_capacity(total);
    for b in &bricks {
        for &(c, mask, p) in &b.patches {
            let e = cells.entry(c).or_insert((owners.len() as u32, 0));
            e.1 += 1;
            owners.push((mask, verts.len() as u32));
            verts.push(p);
        }
    }
    let tris: Vec<[u32; 3]> = bricks.par_iter().flat_map_iter(|b| split_quads(tape, ctx, &verts, quads(b, &cells, &owners))).collect();
    (Mesh { verts, tris }, unresolved)
}

/// A brick's cells: `tape` simplified over it, `whole` the field anywhere.
fn brick(tape: &Tape, whole: &Arc<Tape>, ctx: &Ctx, grid: &Grid, c: Cell, check: bool) -> Option<Brick> {
    let i0 = c.i;
    let mut pts = Vec::with_capacity(N * N * N);
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                pts.push(grid.pos([i0[0] + i as i64, i0[1] + j as i64, i0[2] + k as i64]));
            }
        }
    }
    let mut vals = Vec::new();
    eval_points(tape, ctx, &pts, &mut vals);

    // Pass 1: sign-changing cells, and the centres of their saddle faces;
    // with `check`, the centres of uniform cells near the surface.
    let mut mixed: Vec<([usize; 3], [bool; 8])> = Vec::new();
    let mut face_keys: HashMap<([i64; 3], usize), usize> = HashMap::new();
    let mut centers: Vec<V3> = Vec::new();
    let mut uniform: Vec<(V3, bool)> = Vec::new();
    let reach = 0.5 * 3f64.sqrt() * grid.h;
    for k in 0..B as usize {
        for j in 0..B as usize {
            for i in 0..B as usize {
                let corner = |c: usize| vals[vidx(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1))];
                let ins: [bool; 8] = std::array::from_fn(|c| inside(corner(c)));
                if ins.iter().all(|&v| v) || ins.iter().all(|&v| !v) {
                    if check && (0..8).any(|c| corner(c).abs() <= reach) {
                        let p = grid.pos([i0[0] + i as i64, i0[1] + j as i64, i0[2] + k as i64]);
                        uniform.push((p.map(|x| x + 0.5 * grid.h), ins[0]));
                    }
                    continue;
                }
                for f in faces() {
                    let r = f.ring.map(|c| ins[c]);
                    if r[0] == r[2] && r[1] == r[3] && r[0] != r[1] {
                        let key = face_key([i0[0] + i as i64, i0[1] + j as i64, i0[2] + k as i64], f);
                        face_keys.entry(key).or_insert_with(|| {
                            centers.push(face_center(grid, key));
                            centers.len() - 1
                        });
                    }
                }
                mixed.push(([i, j, k], ins));
            }
        }
    }
    let mut missed = (0, 0, 0.0, None);
    if !uniform.is_empty() {
        let pts: Vec<V3> = uniform.iter().map(|u| u.0).collect();
        let mut v = Vec::new();
        eval_points(tape, ctx, &pts, &mut v);
        let mut material = Vec::new();
        for (&(p, ins), &f) in uniform.iter().zip(&v) {
            if inside(f) != ins {
                (missed.0, missed.3) = (missed.0 + 1, missed.3.or(Some(p)));
                if inside(f) {
                    material.push(p);
                }
            }
        }
        missed.1 = material.len();
        // Across the material, out of the brick too: the field focused on
        // the brick and a cell round it.
        if !material.is_empty() {
            let (lo, hi) = grid.cell_bounds(c);
            let near = crate::query::focus(whole, ctx, lo.map(|x| x - grid.h), hi.map(|x| x + grid.h));
            missed.2 = thickest(&near, ctx, &material, grid.h);
        }
    }
    if mixed.is_empty() {
        return (missed.0 > 0).then(|| Brick { i0, vals, patches: Vec::new(), missed });
    }
    let mut center_vals = Vec::new();
    eval_points(tape, ctx, &centers, &mut center_vals);

    // Pass 2: patches and their edge crossings.
    let mut owners: Vec<([i64; 3], u16, usize, usize)> = Vec::new(); // cell, mask, crossing start, count
    let mut brackets: Vec<(V3, f64, V3, f64)> = Vec::new();
    for ([i, j, k], ins) in &mixed {
        let cell = [i0[0] + *i as i64, i0[1] + *j as i64, i0[2] + *k as i64];
        let corner = |c: usize| vidx(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1));
        let fc = |fi: usize| inside(center_vals[face_keys[&face_key(cell, &faces()[fi])]]);
        for mask in patches(ins, fc) {
            let start = brackets.len();
            for (e, &(a, b)) in EDGES.iter().enumerate() {
                if mask & (1 << e) == 0 {
                    continue;
                }
                let (a, b) = (corner(a), corner(b));
                let (pa, va, pb, vb) = (pts[a], vals[a], pts[b], vals[b]);
                brackets.push(if inside(va) { (pa, va, pb, vb) } else { (pb, vb, pa, va) });
            }
            owners.push((cell, mask, start, brackets.len() - start));
        }
    }
    let cross = refine_crossings(tape, ctx, &mut brackets);
    let grads = eval_grads(tape, ctx, &cross);
    let placed: Vec<V3> = owners.iter().map(|&(cell, _, s, n)| qef(&cross[s..s + n], &grads[s..s + n], grid.pos(cell), grid.h)).collect();
    let patches = owners.iter().zip(placed).map(|(&(cell, mask, _, _), v)| (cell, mask, v)).collect();
    Some(Brick { i0, vals, patches, missed })
}

/// Samples a cell across, either way, on a line through the material.
const ACROSS: i32 = 16;

/// The thickest the material is at points in it: at each, its shortest
/// crossing along the field's gradient there and the three axes (one along
/// a wall reads the wall's length: the gradient at a wall's corner lies in
/// the other wall), each from where the field crosses zero either side
/// (between samples a sixteenth of `h` apart, up to `h` away).
fn thickest(tape: &Tape, ctx: &Ctx, pts: &[V3], h: f64) -> f64 {
    if pts.is_empty() {
        return 0.0;
    }
    let step = h / ACROSS as f64;
    let mut line = Vec::with_capacity(4 * pts.len() * (2 * ACROSS as usize + 1));
    for (p, g) in pts.iter().zip(eval_grads(tape, ctx, pts)) {
        let n = norm([g[1], g[2], g[3]]);
        let u = if n > 0.0 && n.is_finite() { [g[1] / n, g[2] / n, g[3] / n] } else { [0.0, 0.0, 1.0] };
        for u in [u, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            line.extend((-ACROSS..=ACROSS).map(|i| std::array::from_fn(|a| p[a] + step * i as f64 * u[a])));
        }
    }
    let mut v = Vec::new();
    eval_points(tape, ctx, &line, &mut v);
    let across = |f: &[f64]| {
        // The run of samples inside through the middle, and how far past
        // each end of it the field crosses zero.
        let (mut lo, mut hi) = (ACROSS as usize, ACROSS as usize);
        while lo > 0 && inside(f[lo - 1]) {
            lo -= 1;
        }
        while hi + 1 < f.len() && inside(f[hi + 1]) {
            hi += 1;
        }
        let past = |a: usize, b: Option<usize>| b.map_or(0.0, |b| f[a] / (f[a] - f[b]));
        step * ((hi - lo) as f64 + past(lo, lo.checked_sub(1)) + past(hi, (hi + 1 < f.len()).then_some(hi + 1)))
    };
    v.chunks(4 * (2 * ACROSS as usize + 1)).map(|f| f.chunks(2 * ACROSS as usize + 1).map(across).fold(f64::INFINITY, f64::min)).fold(0.0, f64::max)
}

/// Canonical key of a cell face: (global min corner of the face, normal axis).
fn face_key(cell: [i64; 3], f: &Face) -> ([i64; 3], usize) {
    let mut g = cell;
    g[f.axis] += ((f.ring[0] >> f.axis) & 1) as i64;
    (g, f.axis)
}

/// Face centre, computed only from the canonical key so every cell sharing the
/// face gets bit-identical coordinates (and therefore the same field value).
fn face_center(grid: &Grid, (g, axis): ([i64; 3], usize)) -> V3 {
    std::array::from_fn(|a| {
        let i = g[a] as f64 + if a == axis { 0.0 } else { 0.5 };
        grid.origin[a] + grid.h * i
    })
}

const REFINE_STEPS: usize = 6;

/// Move each edge crossing onto the surface: batched false-position steps,
/// safeguarded to shrink the bracket every iteration. Linear interpolation
/// alone lands off the surface near sharp features, where distance fields curve.
fn refine_crossings(tape: &Tape, ctx: &Ctx, br: &mut [(V3, f64, V3, f64)]) -> Vec<V3> {
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
    br.iter()
        .map(|&(pa, va, pb, vb)| {
            let t = if va == vb { 0.5 } else { (va / (va - vb)).clamp(0.0, 1.0) };
            [pa[0] + (pb[0] - pa[0]) * t, pa[1] + (pb[1] - pa[1]) * t, pa[2] + (pb[2] - pa[2]) * t]
        })
        .collect()
}

/// Minimize sum (n . (v - p))^2 + lambda |v - mass|^2, clamped near the cell.
fn qef(pts: &[V3], grads: &[[f64; 4]], cmin: V3, h: f64) -> V3 {
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
    std::array::from_fn(|a| (mass[a] + x[a]).clamp(cmin[a] - m, cmin[a] + h + m))
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
fn split_quads(tape: &Tape, ctx: &Ctx, verts: &[V3], quads: Vec<[u32; 4]>) -> Vec<[u32; 3]> {
    let mid = |a: u32, b: u32| {
        let (p, q) = (verts[a as usize], verts[b as usize]);
        [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0, (p[2] + q[2]) / 2.0]
    };
    let pts: Vec<V3> = quads.iter().flat_map(|q| [mid(q[0], q[2]), mid(q[1], q[3])]).collect();
    let mut f = Vec::new();
    eval_points(tape, ctx, &pts, &mut f);
    let mut out = Vec::with_capacity(quads.len() * 2);
    for (k, q) in quads.iter().enumerate() {
        if f[2 * k].abs() <= f[2 * k + 1].abs() {
            out.push([q[0], q[1], q[2]]);
            out.push([q[0], q[2], q[3]]);
        } else {
            out.push([q[0], q[1], q[3]]);
            out.push([q[1], q[2], q[3]]);
        }
    }
    out
}

fn quads(b: &Brick, cells: &HashMap<[i64; 3], (u32, u8)>, owners: &[(u16, u32)]) -> Vec<[u32; 4]> {
    let mut out = Vec::new();
    if b.patches.is_empty() {
        return out;
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
                    let off = |du: i64, dv: i64| {
                        let mut c = g;
                        c[u] -= du;
                        c[v] -= dv;
                        vertex(c, g, d)
                    };
                    let (Some(q00), Some(q10), Some(q11), Some(q01)) = (off(1, 1), off(0, 1), off(0, 0), off(1, 0)) else {
                        continue;
                    };
                    // Counter-clockwise around +d when the surface normal points along +d.
                    out.push(if inside(v0) { [q00, q10, q11, q01] } else { [q00, q01, q11, q10] });
                }
            }
        }
    }
    out
}

// ---------------------------------------------------------------- analysis

#[derive(Clone, Debug)]
pub struct Shell {
    pub triangles: usize,
    pub volume: f64,
    pub area: f64,
    pub lo: V3,
    pub hi: V3,
}

#[derive(Clone, Debug)]
pub struct Analysis {
    pub area: f64,
    pub volume: f64,
    pub boundary_edges: usize,
    pub nonmanifold_edges: usize,
    pub shells: Vec<Shell>,
}

fn tri_terms(m: &Mesh, t: &[u32; 3]) -> (f64, f64) {
    let (a, b, c) = (m.verts[t[0] as usize], m.verts[t[1] as usize], m.verts[t[2] as usize]);
    let area = 0.5 * norm(normal([a, b, c]));
    let vol = dot(a, cross(b, c)) / 6.0;
    (area, vol)
}

/// The root of `x` in a union-find forest of parent links, halving the
/// path to it on the way.
pub fn find(p: &mut [u32], mut x: u32) -> u32 {
    while p[x as usize] != x {
        p[x as usize] = p[p[x as usize] as usize];
        x = p[x as usize];
    }
    x
}

pub fn analyze(m: &Mesh) -> Analysis {
    let mut keys: Vec<u64> = m
        .tris
        .par_iter()
        .flat_map_iter(|t| {
            (0..3).map(move |e| {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                ((a.min(b) as u64) << 32) | a.max(b) as u64
            })
        })
        .collect();
    keys.par_sort_unstable();
    let (mut boundary_edges, mut nonmanifold_edges) = (0, 0);
    for run in keys.chunk_by(|a, b| a == b) {
        match run.len() {
            1 => boundary_edges += 1,
            2 => {}
            _ => nonmanifold_edges += 1,
        }
    }

    // Union-find over vertices to split into connected shells.
    let mut parent: Vec<u32> = (0..m.verts.len() as u32).collect();
    for t in &m.tris {
        let r0 = find(&mut parent, t[0]);
        for &v in &t[1..] {
            let r = find(&mut parent, v);
            parent[r as usize] = r0;
        }
    }
    let mut shells: HashMap<u32, Shell> = HashMap::new();
    let (mut area, mut volume) = (0.0, 0.0);
    for t in &m.tris {
        let (a, v) = tri_terms(m, t);
        area += a;
        volume += v;
        let r = find(&mut parent, t[0]);
        let (lo, hi) = EMPTY_BOX;
        let s = shells.entry(r).or_insert(Shell { triangles: 0, volume: 0.0, area: 0.0, lo, hi });
        s.triangles += 1;
        s.area += a;
        s.volume += v;
        for &i in t {
            let p = m.verts[i as usize];
            (s.lo, s.hi) = box_union((s.lo, s.hi), (p, p));
        }
    }
    let mut shells: Vec<Shell> = shells.into_values().collect();
    shells.sort_by(|a, b| b.volume.abs().total_cmp(&a.volume.abs()));
    Analysis { area, volume, boundary_edges, nonmanifold_edges, shells }
}
