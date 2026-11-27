//! Step 6: final numbering and element matrices (`assemble`), and what reads
//! the assembled system (`Model`).

use crate::buckle::HeldGeometry;
use crate::cells::{ALPHA, SimCell, place};
use crate::element::{element_basis, full_moments, laminate_matrix, shape_n};
use crate::system::{DofMatrix, System};
use crate::{Physics, RegionSel};
use ink_core::V3;
use rayon::prelude::*;

use super::held::{HeldCells, Lent};
use super::holds::Holds;
use super::settle::Cells;
use super::study::Study;

/// The final cells' numbering: the kept cells and their nodes.
pub(super) struct Mesh {
    /// Each node's place on the lattice, and its position.
    pub(super) final_idx: Vec<[i64; 3]>,
    pub(super) final_pos: Vec<V3>,
    pub(super) nn: usize,
    /// Each classified cell's final cell (u32::MAX: dropped), and each final cell's classified one.
    pub(super) final_of: Vec<u32>,
    pub(super) classified_of: Vec<usize>,
    /// Each node of the cells' numbering's final node (u32::MAX: dropped).
    pub(super) remap: Vec<u32>,
    /// The first final cell at each grid cell (u32::MAX: none).
    pub(super) cell_at: Vec<u32>,
    /// The hanging nodes, by final node.
    pub(super) hang: Vec<(u32, Vec<(u32, f64)>)>,
    /// (final cell, clamp) for reactions.
    pub(super) clamp_cells: Vec<(usize, usize)>,
    /// Each final cell's held material's geometry (`HeldGeometry`).
    pub(super) held_geometry: Vec<Option<Box<HeldGeometry>>>,
    /// The nodes inside a region's held material: what is read at them is the
    /// held value, not the free field's extension past the region's edge
    /// that they carry for the cells it cuts (a region held at 100 C read
    /// 101.4 C there). Each with its fix and held value (before its free
    /// motions).
    pub(super) held_nodes: Vec<(u32, usize, [f64; 3])>,
    pub(super) n_full: usize,
    pub(super) n_cut: usize,
    /// The material's volume (mm^3).
    pub(super) volume: f64,
}

/// What reads the assembled system: the setup, the cells, their holds, the
/// lent tractions and the final numbering.
pub(super) struct Model<'a> {
    pub(super) st: &'a Study<'a>,
    pub(super) cells: &'a Cells,
    pub(super) holds: &'a Holds,
    pub(super) lent: &'a [Lent],
    pub(super) mesh: &'a Mesh,
}

impl Model<'_> {
    /// A lent traction's scale: its two cells' multipliers' geometric mean.
    pub(super) fn lent_scale(&self, sys: &System, l: &Lent) -> f64 {
        lent_scale(&self.mesh.final_of, sys, l)
    }

    pub(super) fn across_of(&self, sys: &System) -> Option<DofMatrix> {
        across_of(self.st, self.cells, self.lent, &self.mesh.final_of, &self.mesh.remap, self.mesh.nn, sys)
    }
}

/// The system (its right-hand side the holds' and lent tractions' held
/// values), and the final numbering. Lent tractions between cells not kept
/// are dropped.
pub(super) fn assemble(st: &Study, settled: &Cells, holds: &Holds, held: &mut HeldCells, keep: &[bool]) -> (System, Vec<f64>, Mesh) {
    let Study { s, ctx, ndof, ne, h, lo, total_cells, ref grid, ref distributing, .. } = *st;
    let Cells { ref classified, depth, ref cell_nodes, ref node_idx, ref hang, .. } = *settled;
    let gained = &settled.slivers.gained;
    let clamps = &holds.clamps;
    let HeldCells { ref mut lent, ref mut free_matrix, ref mut free_material, .. } = *held;
    // (By the cells' final numbering.)
    let mut held_geometry: Vec<Option<Box<HeldGeometry>>> = Vec::new();
    // (Each final cell's place among the classified.)
    let mut classified_of: Vec<usize> = Vec::new();
    // 6. Final numbering and element matrices (clamped cells get explicit ones).
    // One element basis serves every material, each cell with its own
    // coefficients; one full-cell matrix per distinct nu of the full cells (a
    // constant nu gives one, a material per body one each).
    let basis = element_basis(s.physics);
    let with_nu = |nu: f64| match s.physics {
        Physics::Elastic { .. } => Physics::Elastic { nu },
        p => p,
    };
    // (A child's, `MAX_DEPTH`, at its own size.)
    let mut mats: Vec<(f64, u8)> = classified.iter().filter(|c| c.moments.is_none() && c.lam.is_none()).map(|c| (c.nu, c.level)).collect();
    let by = |a: &(f64, u8), b: &(f64, u8)| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1));
    mats.sort_by(by);
    mats.dedup();
    let fm = full_moments();
    let half = |level: u8| 0.5 * h / (1u64 << level) as f64;
    let k_ref: Vec<Vec<f64>> = mats.iter().map(|&(nu, level)| basis.matrix(with_nu(nu), &fm, half(level))).collect();
    let mat_of = |nu: f64, level: u8| mats.binary_search_by(|m| by(m, &(nu, level))).ok();
    let mut clamp_of: Vec<Vec<usize>> = vec![Vec::new(); classified.len()];
    for (k, cl) in clamps.iter().enumerate() {
        clamp_of[cl.cell].push(k);
    }
    let mut remap = vec![u32::MAX; node_idx.len()];
    let mut final_idx: Vec<[i64; 3]> = Vec::new();
    let mut cells = Vec::new();
    let mut k_cut = Vec::new();
    let mut cell_at = vec![u32::MAX; total_cells as usize];
    let mut f = Vec::new();
    // (final cell, clamp) for reactions.
    let mut clamp_cells: Vec<(usize, usize)> = Vec::new();
    let (mut n_full, mut n_cut, mut volume) = (0, 0, 0.0);
    let mut final_of = vec![u32::MAX; classified.len()];
    for (ci, c) in classified.iter().enumerate() {
        if !keep[ci] {
            continue;
        }
        final_of[ci] = cells.len() as u32;
        let mut nodes = [0u32; 8];
        for k in 0..8 {
            let old = cell_nodes[ci][k] as usize;
            if remap[old] == u32::MAX {
                remap[old] = final_idx.len() as u32;
                final_idx.push(node_idx[old]);
                f.extend(std::iter::repeat_n(0.0, ndof));
            }
            nodes[k] = remap[old];
        }
        let size = place(grid, c).1;
        let explicit = c.moments.is_some() || c.lam.is_some() || !clamp_of[ci].is_empty() || free_matrix.contains_key(&ci) || gained.contains_key(&ci);
        let cut = if !explicit {
            n_full += 1;
            u32::MAX
        } else {
            let mut ke = match &c.moments {
                // A laminate (`laminate`), cut or full.
                _ if c.lam.is_some() => {
                    let lam = c.lam.as_deref().expect("a laminate");
                    let m = c.moments.map_or_else(full_moments, |m| m.map(|v| v / c.coef));
                    let mut ke = laminate_matrix(lam, s.physics, &m, half(c.level));
                    if c.moments.is_some() {
                        n_cut += 1;
                        for (a, b) in ke.iter_mut().zip(laminate_matrix(lam, s.physics, &fm, half(c.level))) {
                            *a += ALPHA * b;
                        }
                    } else {
                        n_full += 1;
                    }
                    ke
                }
                Some(m) => {
                    n_cut += 1;
                    // (with the slivers it took: `gained`)
                    let m: [f64; 27] = std::array::from_fn(|k| m[k] + gained.get(&ci).map_or(0.0, |g| g.0[k]));
                    let mut ke = basis.matrix(with_nu(c.nu), &m, half(c.level));
                    let full = match mat_of(c.nu, c.level) {
                        Some(k) => std::borrow::Cow::Borrowed(&k_ref[k]),
                        None => std::borrow::Cow::Owned(basis.matrix(with_nu(c.nu), &fm, half(c.level))),
                    };
                    for (a, b) in ke.iter_mut().zip(full.iter()) {
                        *a += ALPHA * c.coef * b;
                    }
                    ke
                }
                None => {
                    n_full += 1;
                    match gained.get(&ci) {
                        // (a full cell with the slivers it took: one element
                        // of both, its inner modes condensed over both, not
                        // the sum of two)
                        Some((m, _)) => basis.matrix(with_nu(c.nu), &std::array::from_fn(|k| fm[k] * c.coef + m[k]), half(c.level)),
                        None => k_ref[mat_of(c.nu, c.level).expect("a full cell's nu")].iter().map(|v| v * c.coef).collect(),
                    }
                }
            };
            // A cell holding held material: its free material's own.
            if let Some(kf) = free_matrix.remove(&ci) {
                ke = kf;
            }
            for &k in &clamp_of[ci] {
                let cl = &clamps[k];
                clamp_cells.push((cells.len(), k));
                if distributing[cl.fix] {
                    continue;
                }
                ke.iter_mut().zip(&cl.t.kh).for_each(|(k, v)| *k += v);
                for (i, g) in cl.t.gu.iter().enumerate() {
                    f[nodes[i / ndof] as usize * ndof + i % ndof] += g;
                }
            }
            k_cut.push(ke);
            (k_cut.len() - 1) as u32
        };
        volume += (c.frac + gained.get(&ci).map_or(0.0, |g| g.1)) * size * size * size;
        // The first of a cell's pieces (the rest follow it).
        if cell_at[st.lin(c.idx)] == u32::MAX {
            cell_at[st.lin(c.idx)] = cells.len() as u32;
        }
        held_geometry.push(free_material.remove(&ci).map(Box::new));
        classified_of.push(ci);
        let mat = mat_of(c.nu, c.level).unwrap_or(0);
        // (with the slivers it took, as its stiffness has them: a
        // buckling's geometric stiffness is made of these, and without
        // them a column's factor read high by the slivers' share)
        let material: Option<Box<[f64; 27]>> = match (c.moments, gained.get(&ci)) {
            (Some(m), g) => Some(Box::new(std::array::from_fn(|k| (m[k] + g.map_or(0.0, |g| g.0[k])) / c.coef))),
            (None, Some(g)) => Some(Box::new(std::array::from_fn(|k| fm[k] + g.0[k] / c.coef))),
            (None, None) => None,
        };
        cells.push(SimCell {
            idx: c.idx,
            nodes,
            frac: c.frac,
            coef: c.coef,
            cut,
            tape: c.tape.clone(),
            mult: 1.0,
            nu: c.nu,
            mat: mat as u32,
            material,
            mask: c.mask,
            level: c.level,
            sub: c.sub,
            lam: c.lam.clone(),
        });
    }
    let nn = final_idx.len();
    let hang: Vec<(u32, Vec<(u32, f64)>)> =
        hang.iter().filter(|h| remap[h.0 as usize] != u32::MAX).map(|(n, m)| (remap[*n as usize], m.iter().map(|&(m, w)| (remap[m as usize], w)).collect())).collect();
    let mut sys = System::new(ndof, cells, k_ref, k_cut, nn, &hang, 1 << depth);
    // The lent tractions (`Lent`): -w N_i (P t) between the piece's cell's
    // dofs and the lender's, and the same transposed; and the held value's
    // work on the lender's traction. Each scaled as its two cells' materials
    // are (a redesign's `mult`), by their geometric mean: the piece's own
    // cell's constant is scaled as that cell, the lender's energy as the
    // lender, and the mean keeps the two bounding the term between them.
    lent.retain(|l| final_of[l.from] != u32::MAX && final_of[l.to] != u32::MAX);
    let lent: &[Lent] = lent;
    // (called again by a redesign, at its cells' scales)
    sys.across = across_of(st, settled, lent, &final_of, &remap, nn, &sys);
    // (The held value's work is at the cells as they are: a redesign, which
    // scales them, holds its fixes still, refusing a moved one by name, as
    // the cells' own held values' are not rescaled either.)
    for l in lent {
        for col in 0..ne {
            let n = sys.cells[final_of[l.to] as usize].nodes[col / ndof];
            let to: Vec<(u32, f64)> = match sys.hang_of.get(n as usize) {
                Some(&h) if h != u32::MAX => sys.hang[h as usize].clone(),
                _ => vec![(n, 1.0)],
            };
            for a in 0..ndof {
                for &(q, wq) in &to {
                    f[q as usize * ndof + col % ndof] -= l.w * l.pf[a * ne + col] * l.t.v[a] * wq;
                }
            }
        }
    }
    // The lattice of nodes: its step, and a node's place.
    let step = h / (1i64 << depth) as f64;
    let lattice_pos = |i: [i64; 3]| -> V3 { std::array::from_fn(|a| lo[a] + step * i[a] as f64) };
    let final_pos: Vec<V3> = final_idx.par_iter().map(|&g| lattice_pos(g)).collect();
    // The nodes inside a region's held material: what is read at them is the
    // held value, not the free field's extension past the region's edge
    // that they carry for the cells it cuts (a region held at 100 C read
    // 101.4 C there). Each with its fix and held value (before its free
    // motions).
    let held_nodes: Vec<(u32, usize, [f64; 3])> = {
        let mut on: Vec<(u32, usize)> = clamp_cells
            .iter()
            .filter(|&&(_, k)| clamps[k].hn.iter().flatten().any(|v| *v > 0.0))
            .flat_map(|&(cell, k)| sys.cells[cell].nodes.map(|n| (n, clamps[k].fix)))
            .collect();
        on.sort_unstable();
        on.dedup();
        let pts: Vec<V3> = on.iter().map(|x| final_pos[x.0 as usize]).collect();
        let inside = crate::eval::eval_points_par(&s.shape, ctx, &pts);
        let mut out = Vec::new();
        for (fi, fx) in s.fixes.iter().enumerate() {
            let RegionSel::Field(t) = &fx.region else { continue };
            // (a node on the part's face reads outside it by a rounding: one
            // such, in a cell a hold's face left 5e-6 free, read its field's
            // extension, 0.10074 for the 0.1 it is held at)
            let mine: Vec<usize> = (0..on.len()).filter(|&i| on[i].1 == fi && inside[i] <= 1e-9 * h).collect();
            let at: Vec<V3> = mine.iter().map(|&i| pts[i]).collect();
            let within = crate::eval::eval_points_par(t, ctx, &at);
            let values = st.values_at(fx, &at);
            out.extend(mine.iter().zip(within).zip(values).filter(|((_, r), _)| *r < -1e-9 * h).map(|((&i, _), v)| (on[i].0, fi, v)));
        }
        out.sort_unstable_by_key(|x| x.0);
        out
    };
    let mesh = Mesh { final_idx, final_pos, nn, final_of, classified_of, remap, cell_at, hang, clamp_cells, held_geometry, held_nodes, n_full, n_cut, volume };
    (sys, f, mesh)
}

fn lent_scale(final_of: &[u32], sys: &System, l: &Lent) -> f64 {
    (sys.cells[final_of[l.from] as usize].mult * sys.cells[final_of[l.to] as usize].mult).sqrt()
}

// The lent tractions (`Lent`): -w N_i (P t) between the piece's cell's
// dofs and the lender's, and the same transposed; and the held value's
// work on the lender's traction. Each scaled as its two cells' materials
// are (a redesign's `mult`), by their geometric mean: the piece's own
// cell's constant is scaled as that cell, the lender's energy as the
// lender, and the mean keeps the two bounding the term between them.
fn across_of(st: &Study, cells: &Cells, lent: &[Lent], final_of: &[u32], remap: &[u32], nn: usize, sys: &System) -> Option<DofMatrix> {
    let Study { ndof, ne, .. } = *st;
    let cell_nodes = &cells.cell_nodes;
    if lent.is_empty() {
        return None;
    }
    let expand = |n: u32| -> Vec<(u32, f64)> {
        match sys.hang_of.get(n as usize) {
            Some(&h) if h != u32::MAX => sys.hang[h as usize].clone(),
            _ => vec![(n, 1.0)],
        }
    };
    let mut entries = std::collections::HashMap::new();
    for l in lent {
        let lists: Vec<Vec<(u32, f64)>> = cell_nodes[l.from].iter().chain(&cell_nodes[l.to]).map(|&n| expand(remap[n as usize])).collect();
        let scale = lent_scale(final_of, sys, l);
        for col in 0..ne {
            let (to, d) = (&lists[8 + col / ndof], (col % ndof) as u32);
            for a in (0..ndof).filter(|&a| l.pf[a * ne + col] != 0.0) {
                for i in 0..8 {
                    let v = -scale * l.w * shape_n(i, l.t.xi) * l.pf[a * ne + col];
                    for &(p, wp) in &lists[i] {
                        for &(q, wq) in to {
                            let (r, c) = (p * ndof as u32 + a as u32, q * ndof as u32 + d);
                            *entries.entry((r, c)).or_insert(0.0) += v * wp * wq;
                            *entries.entry((c, r)).or_insert(0.0) += v * wp * wq;
                        }
                    }
                }
            }
        }
    }
    Some(DofMatrix::from_entries(nn * ndof, &entries))
}
