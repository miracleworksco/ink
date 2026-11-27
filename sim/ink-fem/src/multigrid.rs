//! Geometric multigrid on the cell grid, used as the CG preconditioner (so the
//! solution is the same; only the iteration count changes). A coarse cell
//! merges 2x2x2 fine cells and its matrix is the exact Galerkin product
//! P^T K P of its children: trilinear coarse functions are exactly representable
//! on the fine grid, so holes, lattices, material fields and clamps all carry
//! down with no resampling. Chebyshev-Jacobi smoothing; dense Cholesky on the
//! coarsest level.
//! And the preconditioned conjugate gradients it serves.

use crate::cells::{ALL, SimCell};
use crate::numerics::dot;
use crate::system::System;
use crate::{GRAIN, NODE_GRAIN};
use ink_kernel::mesh::find;
use rayon::prelude::*;

const CHEB_DEGREE: usize = 3;
const COARSEST_DOFS: usize = 1000;

/// CSR map from the nodes of one level to the nodes of the next coarser one.
pub(crate) struct Transfer {
    pub(crate) rows: Vec<u32>,
    pub(crate) cols: Vec<u32>,
    pub(crate) w: Vec<f64>,
    /// Transpose (coarse row -> fine entries), for restriction.
    t_rows: Vec<u32>,
    t_cols: Vec<u32>,
    t_w: Vec<f64>,
}

impl Transfer {
    /// fine += P coarse
    fn prolong_add(&self, nd: usize, coarse: &[f64], fine: &mut [f64]) {
        fine.par_chunks_mut(nd).with_min_len(NODE_GRAIN).enumerate().for_each(|(n, out)| {
            for e in self.rows[n] as usize..self.rows[n + 1] as usize {
                let (c, w) = (self.cols[e] as usize, self.w[e]);
                for d in 0..nd {
                    out[d] += w * coarse[c * nd + d];
                }
            }
        });
    }

    /// coarse = P^T fine
    fn restrict(&self, nd: usize, fine: &[f64], ncoarse: usize) -> Vec<f64> {
        let mut out = vec![0.0; ncoarse * nd];
        out.par_chunks_mut(nd).with_min_len(NODE_GRAIN).enumerate().for_each(|(c, o)| {
            for e in self.t_rows[c] as usize..self.t_rows[c + 1] as usize {
                let (n, w) = (self.t_cols[e] as usize, self.t_w[e]);
                for d in 0..nd {
                    o[d] += w * fine[n * nd + d];
                }
            }
        });
        out
    }
}

pub(crate) struct Level {
    sys: System,
    inv: Vec<f64>,
    lam: f64,
}

enum Coarsest {
    Cholesky { n: usize, l: Vec<f64> },
    Smooth,
}

pub(crate) struct Multigrid {
    fine_inv: Vec<f64>,
    fine_lam: f64,
    /// Levels 1.. (level 0 is the fine system, owned by the caller).
    pub(crate) levels: Vec<Level>,
    /// transfers[l]: nodes of level l -> nodes of level l + 1.
    transfers: Vec<Transfer>,
    coarsest: Coarsest,
}

fn inverse_diagonal(sys: &System) -> Vec<f64> {
    sys.diagonal().iter().map(|d| if *d > 0.0 { 1.0 / d } else { 0.0 }).collect()
}

/// Largest eigenvalue of D^-1 K (with a safety margin): the largest of the
/// tridiagonal matrix `LANCZOS_STEPS` steps of Lanczos's iteration build
/// for D^-1/2 K D^-1/2, which has the same. As many products by K as the
/// power iteration this was took, and far nearer: a power iteration nears
/// the largest by the ratio of the next to it each step, and where a few
/// cells carry the largest (a hold's constant far over the rest, a
/// lattice's thin cut cells) or many lie close under it (any fine grid)
/// 15 steps stop well short. The smoother is then tuned low and amplifies
/// what it should damp. A plain beam's finest grid read 3.33 for 3.69 or
/// more; a bar held over a cylinder of radius 1.37 on 0.5 mm cells 3.32 for
/// 3.74, and took 32 iterations for 12.
fn max_eigen(sys: &System, inv: &[f64]) -> f64 {
    let n = sys.ndofs;
    let half: Vec<f64> = inv.iter().map(|i| i.sqrt()).collect();
    let mut seed = 0x9E37_79B9_7F4A_7C15u64;
    let mut v: Vec<f64> = (0..n)
        .map(|_| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        })
        .collect();
    let norm = dot(&v, &v).sqrt().max(1e-300);
    v.par_iter_mut().with_min_len(GRAIN).for_each(|x| *x /= norm);
    let (mut last, mut scaled, mut w) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    // The tridiagonal matrix: its diagonal, and what couples each row to the next.
    let (mut diagonal, mut beside): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
    let mut b = 0.0;
    for _ in 0..LANCZOS_STEPS {
        scaled.par_iter_mut().with_min_len(GRAIN).zip(&v).zip(&half).for_each(|((s, x), h)| *s = x * h);
        sys.apply(&scaled, &mut w);
        w.par_iter_mut().with_min_len(GRAIN).zip(&half).zip(&last).for_each(|((x, h), l)| *x = *x * h - b * l);
        let a = dot(&w, &v);
        w.par_iter_mut().with_min_len(GRAIN).zip(&v).for_each(|(x, y)| *x -= a * y);
        diagonal.push(a);
        b = dot(&w, &w).sqrt();
        // (the vectors so far span all the start reaches)
        if b <= 1e-12 * a.abs() {
            break;
        }
        beside.push(b);
        std::mem::swap(&mut last, &mut v);
        v.par_iter_mut().with_min_len(GRAIN).zip(&w).for_each(|(x, y)| *x = y / b);
    }
    beside.truncate(diagonal.len() - 1);
    1.1 * largest_of_tridiagonal(&diagonal, &beside)
}

/// The steps of Lanczos's iteration `max_eigen` takes.
const LANCZOS_STEPS: usize = 15;

/// The largest eigenvalue of the symmetric tridiagonal matrix of `diagonal`
/// and `beside` it: by bisection on how many lie under a value (the signs
/// of the pivots of the matrix less that value).
fn largest_of_tridiagonal(diagonal: &[f64], beside: &[f64]) -> f64 {
    let around = |i: usize| if i > 0 { beside[i - 1].abs() } else { 0.0 } + beside.get(i).map_or(0.0, |b| b.abs());
    let (mut lo, mut hi) = (0..diagonal.len()).fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), i| (l.min(diagonal[i] - around(i)), h.max(diagonal[i] + around(i))));
    let under = |x: f64| -> usize {
        let (mut count, mut pivot) = (0, 1.0);
        for i in 0..diagonal.len() {
            pivot = diagonal[i] - x - if i > 0 { beside[i - 1] * beside[i - 1] / pivot } else { 0.0 };
            if pivot == 0.0 {
                pivot = f64::MIN_POSITIVE;
            }
            count += usize::from(pivot < 0.0);
        }
        count
    };
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if under(mid) == diagonal.len() { hi = mid } else { lo = mid }
    }
    hi
}

/// Each of M vectors as `vecs` makes them: as slices, mutable or not.
fn refs<const M: usize>(v: &[Vec<f64>; M]) -> [&[f64]; M] {
    v.each_ref().map(|v| &v[..])
}
fn muts<const M: usize>(v: &mut [Vec<f64>; M]) -> [&mut [f64]; M] {
    v.each_mut().map(|v| &mut v[..])
}

/// Chebyshev-Jacobi smoothing of K x = b (x updated in place), for each of
/// M pairs at once (`System::applies`).
fn chebyshev<const M: usize>(sys: &System, inv: &[f64], lam: f64, degree: usize, x: &mut [Vec<f64>; M], b: [&[f64]; M], zero_start: bool) {
    let (lmax, lmin) = (lam, lam / 30.0);
    let (theta, delta) = ((lmax + lmin) / 2.0, (lmax - lmin) / 2.0);
    let sigma = theta / delta;
    let mut r: [Vec<f64>; M] = if zero_start {
        b.map(|b| b.to_vec())
    } else {
        let mut kx = std::array::from_fn(|_| vec![0.0; sys.ndofs]);
        sys.applies(refs(x), muts(&mut kx));
        std::array::from_fn(|v| b[v].iter().zip(&kx[v]).map(|(a, k)| a - k).collect())
    };
    let mut d: [Vec<f64>; M] = std::array::from_fn(|v| r[v].iter().zip(inv).map(|(r, i)| r * i / theta).collect());
    for (x, d) in x.iter_mut().zip(&d) {
        x.par_iter_mut().with_min_len(GRAIN).zip(d).for_each(|(x, d)| *x += d);
    }
    let mut rho = 1.0 / sigma;
    let mut kd = std::array::from_fn(|_| vec![0.0; sys.ndofs]);
    for _ in 1..degree {
        sys.applies(refs(&d), muts(&mut kd));
        let rho_new = 1.0 / (2.0 * sigma - rho);
        let (a, c) = (rho_new * rho, 2.0 * rho_new / delta);
        for (((x, d), r), kd) in x.iter_mut().zip(&mut d).zip(&mut r).zip(&kd) {
            r.par_iter_mut().with_min_len(GRAIN).zip(kd).for_each(|(r, k)| *r -= k);
            d.par_iter_mut().with_min_len(GRAIN).zip(&*r).zip(inv).for_each(|((d, r), i)| *d = a * *d + c * r * i);
            x.par_iter_mut().with_min_len(GRAIN).zip(&*d).for_each(|(x, d)| *x += d);
        }
        rho = rho_new;
    }
}

/// Cells in colours that share no nodes, for parallel assembly: by the
/// parity of the cell's position, and a cell's pieces (which may share
/// nodes) each in a colour of its own.
pub(crate) fn color_cells(cells: &[SimCell]) -> (Vec<Vec<u32>>, Vec<Vec<(u32, u32)>>) {
    let mut colors: Vec<Vec<u32>> = vec![Vec::new(); 8];
    let mut runs: Vec<Vec<(u32, u32)>> = vec![Vec::new(); 8];
    let parity = |c: &SimCell| ((c.idx[0] & 1) | ((c.idx[1] & 1) << 1) | ((c.idx[2] & 1) << 2)) as usize;
    let mut rank = 0;
    for (i, c) in cells.iter().enumerate() {
        if c.level > 0 {
            match runs[parity(c)].last_mut() {
                Some(r) if cells[r.1 as usize - 1].idx == c.idx && r.1 as usize == i => r.1 += 1,
                _ => runs[parity(c)].push((i as u32, i as u32 + 1)),
            }
            continue;
        }
        rank = if i > 0 && cells[i - 1].idx == c.idx { rank + 1 } else { 0 };
        let k = ((c.idx[0] & 1) | ((c.idx[1] & 1) << 1) | ((c.idx[2] & 1) << 2)) as usize + 8 * rank;
        if colors.len() <= k {
            colors.resize(k + 1, Vec::new());
        }
        colors[k].push(i as u32);
    }
    (colors, runs)
}

/// Coarsen one level: 2x2x2 cells merge, a grid cell's children (`MAX_DEPTH`)
/// with it; matrices by exact Galerkin products, each fine cell's corners
/// weighted by the coarse cell's trilinear functions where they lie.
pub(crate) fn coarsen(fine: &System, node_idx: &[[i64; 3]]) -> (System, Vec<[i64; 3]>, Transfer) {
    use std::collections::HashMap;
    let (nd, ne) = (fine.ndof, fine.ne);
    let sc = fine.scale;
    // A fine cell's nodes as it reaches them: its corners, or what a
    // hanging corner hangs on.
    let reach = |c: &SimCell| -> Vec<u32> {
        let mut v: Vec<u32> = c
            .nodes
            .iter()
            .flat_map(|&n| match fine.hang_of.get(n as usize) {
                Some(&h) if h != u32::MAX => fine.hang[h as usize].iter().map(|m| m.0).collect(),
                _ => vec![n],
            })
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let reached: Vec<Vec<u32>> = fine.cells.par_iter().map(reach).collect();
    // Coarse cells and their children (fine cells).
    let mut cell_of: HashMap<[i64; 3], usize> = HashMap::new();
    let mut groups: Vec<Vec<u32>> = Vec::new();
    let mut gidx: Vec<[i64; 3]> = Vec::new();
    for (i, c) in fine.cells.iter().enumerate() {
        let p = [c.idx[0].div_euclid(2), c.idx[1].div_euclid(2), c.idx[2].div_euclid(2)];
        let slot = *cell_of.entry(p).or_insert_with(|| {
            groups.push(Vec::new());
            gidx.push(p);
            gidx.len() - 1
        });
        groups[slot].push(i as u32);
    }
    // A coarse cell's pieces: its children joined by the fine nodes they
    // share (pieces of material apart on the fine grid stay apart, so the
    // coarse correction can move them apart).
    let mut children: Vec<Vec<u32>> = Vec::new();
    let mut cidx: Vec<[i64; 3]> = Vec::new();
    let mut piece_of = vec![0u32; fine.cells.len()];
    for (g, ch) in groups.iter().enumerate() {
        let mut local: Vec<u32> = (0..ch.len() as u32).collect();
        let mut first: HashMap<u32, u32> = HashMap::new();
        for (i, &ci) in ch.iter().enumerate() {
            for &n in &reached[ci as usize] {
                match first.get(&n) {
                    Some(&j) => {
                        let (ri, rj) = (find(&mut local, i as u32), find(&mut local, j));
                        local[ri as usize] = rj;
                    }
                    None => _ = first.insert(n, i as u32),
                }
            }
        }
        let mut roots: Vec<u32> = Vec::new();
        for (i, &ci) in ch.iter().enumerate() {
            let r = find(&mut local, i as u32);
            let k = match roots.iter().position(|&x| x == r) {
                Some(k) => k,
                None => {
                    roots.push(r);
                    children.push(Vec::new());
                    cidx.push(gidx[g]);
                    roots.len() - 1
                }
            };
            let piece = children.len() - roots.len() + k;
            children[piece].push(ci);
            piece_of[ci as usize] = piece as u32;
        }
    }
    // Coarse nodes: a coarse piece's 8 corners, shared with another piece
    // where a fine node of the one is a fine node of the other.
    let corner = |p: [i64; 3], j: usize| [p[0] + (j & 1) as i64, p[1] + ((j >> 1) & 1) as i64, p[2] + ((j >> 2) & 1) as i64];
    let mut parent: Vec<u32> = (0..8 * cidx.len() as u32).collect();
    let mut first_piece = vec![u32::MAX; fine.ndofs / nd];
    for (ci, nodes) in reached.iter().enumerate() {
        let a = piece_of[ci] as usize;
        for &n in nodes {
            let slot = &mut first_piece[n as usize];
            if *slot == u32::MAX {
                *slot = a as u32;
                continue;
            }
            let b = *slot as usize;
            if a == b {
                continue;
            }
            for j in 0..8 {
                let g = corner(cidx[a], j);
                let d: [i64; 3] = std::array::from_fn(|x| g[x] - cidx[b][x]);
                if d.iter().all(|&v| v == 0 || v == 1) {
                    let k = (d[0] | d[1] << 1 | d[2] << 2) as usize;
                    let (ra, rb) = (find(&mut parent, (8 * a + j) as u32), find(&mut parent, (8 * b + k) as u32));
                    parent[ra as usize] = rb;
                }
            }
        }
    }
    let mut node_of = vec![u32::MAX; parent.len()];
    let mut cnode_idx: Vec<[i64; 3]> = Vec::new();
    let mut cells = Vec::with_capacity(cidx.len());
    for (k, p) in cidx.iter().enumerate() {
        let mut nodes = [0u32; 8];
        for (j, n) in nodes.iter_mut().enumerate() {
            let r = find(&mut parent, (8 * k + j) as u32) as usize;
            if node_of[r] == u32::MAX {
                node_of[r] = cnode_idx.len() as u32;
                cnode_idx.push(corner(*p, j));
            }
            *n = node_of[r];
        }
        cells.push(SimCell {
            idx: *p,
            nodes,
            frac: 1.0,
            coef: 1.0,
            cut: k as u32,
            tape: None,
            mult: 1.0,
            nu: 0.0,
            mat: 0,
            material: None,
            mask: ALL,
            level: 0,
            sub: [0; 3],
            lam: None,
        });
    }
    // The coarse cell's trilinear functions at a point of the lattice:
    // per axis, (corner bit, weight) where the weight is not 0.
    let axis = |p: [i64; 3], g: [i64; 3], a: usize| -> Vec<(usize, f64)> {
        let t = (g[a] - 2 * p[a] * sc) as f64 / (2 * sc) as f64;
        [(0, 1.0 - t), (1, t)].into_iter().filter(|w| w.1 != 0.0).collect()
    };
    let weights = |p: [i64; 3], g: [i64; 3]| -> Vec<(usize, f64)> {
        let (ax, ay, az) = (axis(p, g, 0), axis(p, g, 1), axis(p, g, 2));
        let mut out = Vec::new();
        for &(z, wz) in &az {
            for &(y, wy) in &ay {
                for &(x, wx) in &ax {
                    out.push((x | y << 1 | z << 2, wx * wy * wz));
                }
            }
        }
        out
    };
    // Galerkin element matrices.
    let k_cut: Vec<Vec<f64>> = children
        .par_iter()
        .zip(&cidx)
        .map(|(ch, &p)| {
            let mut kc = vec![0.0; ne * ne];
            let mut b = vec![0.0; ne * ne];
            for &ci in ch {
                let c = &fine.cells[ci as usize];
                let (a, scale) = fine.matrix(c);
                // w[k][j]: coarse corner j's function at fine corner k.
                let side = sc >> c.level;
                let mut w = [[0.0; 8]; 8];
                for (k, wk) in w.iter_mut().enumerate() {
                    let g: [i64; 3] = std::array::from_fn(|x| c.idx[x] * sc + c.sub[x] as i64 * side + ((k >> x) & 1) as i64 * side);
                    for (j, v) in weights(p, g) {
                        wk[j] = v;
                    }
                }
                // b = A P
                for r in 0..ne {
                    for j in 0..8 {
                        for d in 0..nd {
                            let acc: f64 = (0..8).map(|k| a[(k * nd + d) * ne + r] * w[k][j]).sum();
                            b[r * ne + j * nd + d] = scale * acc;
                        }
                    }
                }
                // kc += P^T b
                for i in 0..8 {
                    for di in 0..nd {
                        for col in 0..ne {
                            let acc: f64 = (0..8).map(|k| w[k][i] * b[(k * nd + di) * ne + col]).sum();
                            kc[(i * nd + di) * ne + col] += acc;
                        }
                    }
                }
            }
            kc
        })
        .collect();
    let ncoarse = cnode_idx.len();
    let mut sys = System::new(nd, cells, Vec::new(), k_cut, ncoarse, &[], 1);
    // Prolongation: each fine node interpolates the corners of a coarse
    // piece holding it (where several do, the corners it reaches are
    // shared); a hanging node none (it is no unknown).
    let (mut rows, mut cols, mut w) = (vec![0u32], Vec::new(), Vec::new());
    for (n, &g) in node_idx.iter().enumerate() {
        let piece = first_piece[n];
        if piece != u32::MAX && fine.hang_of.get(n).is_none_or(|&h| h == u32::MAX) {
            let p = cidx[piece as usize];
            for (j, v) in weights(p, g) {
                cols.push(sys.cells[piece as usize].nodes[j]);
                w.push(v);
            }
        }
        rows.push(cols.len() as u32);
    }
    // Transpose.
    let mut count = vec![0u32; ncoarse + 1];
    for &c in &cols {
        count[c as usize + 1] += 1;
    }
    for i in 0..ncoarse {
        count[i + 1] += count[i];
    }
    let t_rows = count.clone();
    let mut fill = count;
    let (mut t_cols, mut t_w) = (vec![0u32; cols.len()], vec![0.0; cols.len()]);
    for n in 0..node_idx.len() {
        for e in rows[n] as usize..rows[n + 1] as usize {
            let c = cols[e] as usize;
            let slot = fill[c] as usize;
            t_cols[slot] = n as u32;
            t_w[slot] = w[e];
            fill[c] += 1;
        }
    }
    let t = Transfer { rows, cols, w, t_rows, t_cols, t_w };
    // A diagonal shift coarsens to its lumped restriction.
    if !fine.shift.is_empty() {
        sys.shift = t.restrict(nd, &fine.shift, ncoarse);
    }
    // Low-rank terms coarsen exactly: P^T (w g g^T) P = w (P^T g) (P^T g)^T.
    sys.lowrank = fine.lowrank.iter().map(|(w, g)| (*w, t.restrict(nd, g, ncoarse))).collect();
    sys.across = fine.across.as_ref().map(|e| e.coarsen(&t, ncoarse, nd));
    (sys, cnode_idx, t)
}

fn dense_cholesky(sys: &System) -> Option<Coarsest> {
    let n = sys.ndofs;
    let mut a = vec![0.0; n * n];
    for c in &sys.cells {
        let dofs = sys.weighted_dofs(c);
        let (k, s) = sys.matrix(c);
        for (i, di) in dofs.iter().enumerate() {
            for (j, dj) in dofs.iter().enumerate() {
                for &(g, w) in di {
                    for &(g2, w2) in dj {
                        a[g as usize * n + g2 as usize] += s * w * w2 * k[j * sys.ne + i];
                    }
                }
            }
        }
    }
    // A hanging node's dofs are no unknowns: held at 0.
    for (m, &h) in sys.hang_of.iter().enumerate() {
        if h != u32::MAX {
            for d in 0..sys.ndof {
                a[(m * sys.ndof + d) * (n + 1)] = 1.0;
            }
        }
    }
    if let Some(e) = &sys.across {
        for (r, c, v) in e.entries() {
            a[r * n + c] += v;
        }
    }
    for (i, s) in sys.shift.iter().enumerate() {
        a[i * n + i] += s;
    }
    for (w, g) in &sys.lowrank {
        for i in (0..n).filter(|&i| g[i] != 0.0) {
            for j in 0..n {
                a[i * n + j] += w * g[i] * g[j];
            }
        }
    }
    for j in 0..n {
        let mut v = a[j * n + j] - (0..j).map(|k| a[j * n + k] * a[j * n + k]).sum::<f64>();
        if v <= 0.0 {
            return None;
        }
        v = v.sqrt();
        a[j * n + j] = v;
        let (head, tail) = a.split_at_mut((j + 1) * n);
        let rowj = &head[j * n..j * n + j];
        tail.par_chunks_mut(n).for_each(|row| {
            let s: f64 = (0..j).map(|k| row[k] * rowj[k]).sum();
            row[j] = (row[j] - s) / v;
        });
    }
    Some(Coarsest::Cholesky { n, l: a })
}

impl Multigrid {
    pub(crate) fn build(fine: &System, node_idx: &[[i64; 3]]) -> Multigrid {
        let fine_inv = inverse_diagonal(fine);
        let fine_lam = max_eigen(fine, &fine_inv);
        let (mut levels, mut transfers) = (Vec::new(), Vec::new());
        let mut idx = node_idx.to_vec();
        let mut prev_dofs = fine.ndofs;
        loop {
            let last = levels.last().map(|l: &Level| &l.sys).unwrap_or(fine);
            if last.ndofs <= COARSEST_DOFS {
                break;
            }
            let (sys, cidx, t) = coarsen(last, &idx);
            let shrunk = sys.ndofs as f64 / prev_dofs as f64;
            prev_dofs = sys.ndofs;
            let inv = inverse_diagonal(&sys);
            let lam = max_eigen(&sys, &inv);
            levels.push(Level { sys, inv, lam });
            transfers.push(t);
            idx = cidx;
            if shrunk > 0.75 {
                break; // coarsening has stalled (e.g. thin sheets); stop here
            }
        }
        let bottom = levels.last().map(|l| &l.sys).unwrap_or(fine);
        let coarsest = if bottom.ndofs <= COARSEST_DOFS { dense_cholesky(bottom).unwrap_or(Coarsest::Smooth) } else { Coarsest::Smooth };
        Multigrid { fine_inv, fine_lam, levels, transfers, coarsest }
    }

    fn level<'a>(&'a self, l: usize, fine: &'a System) -> (&'a System, &'a [f64], f64) {
        if l == 0 {
            (fine, &self.fine_inv, self.fine_lam)
        } else {
            let lv = &self.levels[l - 1];
            (&lv.sys, &lv.inv, lv.lam)
        }
    }

    /// One V-cycle for K x = b at level l (a symmetric preconditioner).
    pub(crate) fn vcycle(&self, l: usize, fine: &System, b: &[f64]) -> Vec<f64> {
        let [x] = self.vcycles(l, fine, [b]);
        x
    }

    /// `vcycle` of each of `b` (the same bits), two at a time: a level's
    /// cells are read once for both (`System::applies`).
    pub(crate) fn vcycle_each(&self, fine: &System, b: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let mut out = Vec::with_capacity(b.len());
        for b in b.chunks(2) {
            match b {
                [b0, b1] => out.extend(self.vcycles(0, fine, [b0, b1])),
                _ => out.extend(self.vcycles(0, fine, [&b[0]])),
            }
        }
        out
    }

    /// One V-cycle for each of M right-hand sides at once, at level l.
    fn vcycles<const M: usize>(&self, l: usize, fine: &System, b: [&[f64]; M]) -> [Vec<f64>; M] {
        let (sys, inv, lam) = self.level(l, fine);
        let mut x = std::array::from_fn(|_| vec![0.0; sys.ndofs]);
        if l == self.levels.len() {
            match &self.coarsest {
                Coarsest::Cholesky { n, l } => {
                    let n = *n;
                    return b.map(|b| {
                        let mut y = b.to_vec();
                        for i in 0..n {
                            let s: f64 = (0..i).map(|k| l[i * n + k] * y[k]).sum();
                            y[i] = (y[i] - s) / l[i * n + i];
                        }
                        // L^T y = z by rows of L (contiguous), not its columns.
                        for i in (0..n).rev() {
                            y[i] /= l[i * n + i];
                            let yi = y[i];
                            for (yk, lk) in y[..i].iter_mut().zip(&l[i * n..i * n + i]) {
                                *yk -= lk * yi;
                            }
                        }
                        y
                    });
                }
                Coarsest::Smooth => {
                    chebyshev(sys, inv, lam, 30, &mut x, b, true);
                    return x;
                }
            }
        }
        chebyshev(sys, inv, lam, CHEB_DEGREE, &mut x, b, true);
        let mut kx = std::array::from_fn(|_| vec![0.0; sys.ndofs]);
        sys.applies(refs(&x), muts(&mut kx));
        let t = &self.transfers[l];
        let coarse_n = t.t_rows.len() - 1;
        let rc: [Vec<f64>; M] = std::array::from_fn(|v| {
            let r: Vec<f64> = b[v].iter().zip(&kx[v]).map(|(a, k)| a - k).collect();
            t.restrict(sys.ndof, &r, coarse_n)
        });
        let xc = self.vcycles(l + 1, fine, refs(&rc));
        for (x, xc) in x.iter_mut().zip(&xc) {
            t.prolong_add(sys.ndof, xc, x);
        }
        chebyshev(sys, inv, lam, CHEB_DEGREE, &mut x, b, false);
        x
    }
}

/// Preconditioned CG on (K + diag(shift)) x = b from the given start.
/// Returns (iterations, relative residual).
pub(crate) fn pcg(sys: &System, shift: &[f64], precond: &dyn Fn(&[f64]) -> Vec<f64>, b: &[f64], x: &mut [f64], tol: f64, max_iter: usize) -> (usize, f64) {
    let apply = |p: &[f64], out: &mut [f64]| {
        sys.apply(p, out);
        if !shift.is_empty() {
            out.par_iter_mut().with_min_len(GRAIN).zip(p).zip(shift).for_each(|((o, p), s)| *o += s * p);
        }
    };
    let n = b.len();
    let mut ap = vec![0.0; n];
    apply(x, &mut ap);
    let mut r: Vec<f64> = b.iter().zip(&ap).map(|(b, a)| b - a).collect();
    let bnorm = dot(b, b).sqrt().max(1e-300);
    let mut res = dot(&r, &r).sqrt() / bnorm;
    if res <= tol {
        return (0, res);
    }
    let mut z = precond(&r);
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut iters = 0;
    // Stagnation guard: stop (reporting non-convergence) rather than grind
    // when the residual has not halved in the last 500 iterations.
    let (mut checkpoint, mut checkpoint_res) = (0usize, res);
    while iters < max_iter && res > tol {
        if iters - checkpoint >= 500 {
            if res > 0.5 * checkpoint_res {
                break;
            }
            (checkpoint, checkpoint_res) = (iters, res);
        }
        apply(&p, &mut ap);
        // Breakdown: a direction of no positive curvature (a stiffness not
        // positive definite) gives no step: stop, reporting the residual
        // reached; one not finite, a residual not finite.
        let pap = dot(&p, &ap);
        if !pap.is_finite() {
            res = f64::NAN;
            break;
        }
        if pap <= 0.0 {
            break;
        }
        let alpha = rz / pap;
        x.par_iter_mut().with_min_len(GRAIN).zip(&p).for_each(|(x, p)| *x += alpha * p);
        r.par_iter_mut().with_min_len(GRAIN).zip(&ap).for_each(|(r, a)| *r -= alpha * a);
        z = precond(&r);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        p.par_iter_mut().with_min_len(GRAIN).zip(&z).for_each(|(p, z)| *p = z + beta * *p);
        iters += 1;
        res = dot(&r, &r).sqrt() / bnorm;
    }
    (iters, res)
}

#[cfg(test)]
mod tests {

    #[test]
    fn a_tridiagonal_matrixs_largest_eigenvalue_is_found() {
        // 2 on the diagonal and -1 beside it, n rows: 2 - 2 cos(n pi / (n + 1)).
        for n in [1usize, 2, 3, 15] {
            let (got, want) = (super::largest_of_tridiagonal(&vec![2.0; n], &vec![-1.0; n - 1]), 2.0 - 2.0 * (n as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos());
            assert!((got - want).abs() < 1e-12, "{n} rows: {got} for {want}");
        }
        assert!((super::largest_of_tridiagonal(&[1.0, 5.0, 3.0], &[0.0, 0.0]) - 5.0).abs() < 1e-12);
    }
}
