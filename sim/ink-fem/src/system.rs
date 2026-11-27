//! The assembled system: the cells' matrices, applied matrix-free on all
//! cores.

use crate::cells::SimCell;
use crate::multigrid::{Transfer, color_cells};
use crate::numerics::dot;
use crate::{CELL_GRAIN, GRAIN};
use rayon::prelude::*;

pub(crate) struct Ptr(pub(crate) *mut f64);
unsafe impl Sync for Ptr {}
unsafe impl Send for Ptr {}

pub(crate) struct System {
    pub(crate) ndof: usize,
    pub(crate) ne: usize,
    pub(crate) cells: Vec<SimCell>,
    pub(crate) colors: Vec<Vec<u32>>,
    /// Full-cell matrix for unit coefficient, one per distinct nu of the full cells.
    /// It and `k_cut` hold each matrix by columns (`System::new`).
    k_ref: Vec<Vec<f64>>,
    pub(crate) k_cut: Vec<Vec<f64>>,
    pub(crate) ndofs: usize,
    /// Optional diagonal added to the operator (lumped mass terms of a
    /// time step): the system is K + diag(shift).
    pub(crate) shift: Vec<f64>,
    /// Low-rank terms: the system gains w g g^T for each (w, g), g a full
    /// vector of dofs. Supports with free rigid motions (their amplitudes
    /// eliminated) and pinned mechanisms; restricted to coarse levels
    /// exactly, as P^T g.
    pub(crate) lowrank: Vec<(f64, Vec<f64>)>,
    /// Children's cells (`MAX_DEPTH`), in colours as `colors`: runs of one
    /// grid cell's children, each applied in order by one thread (they
    /// share nodes, and the nodes their hanging corners hang on).
    pub(crate) runs: Vec<Vec<(u32, u32)>>,
    /// Per node, its index in `hang` where it hangs (u32::MAX where not;
    /// empty where none does). A hanging node is no unknown: its dofs stay
    /// 0 in the system's vectors (`fold`, `expand`).
    pub(crate) hang_of: Vec<u32>,
    /// The nodes a hanging node interpolates, weighted.
    pub(crate) hang: Vec<Vec<(u32, f64)>>,
    /// Lattice steps per cell of `cells` at level 0 (node indices' unit).
    pub(crate) scale: i64,
    /// Terms between two cells' dofs, added to the cells' (a hold's traction
    /// read in the cell beside its piece's: `Lent`).
    pub(crate) across: Option<DofMatrix>,
}

/// A sparse symmetric matrix over the dofs (CSR): terms that couple the
/// nodes of two cells, which no cell's matrix holds.
pub(crate) struct DofMatrix {
    rows: Vec<u32>,
    cols: Vec<u32>,
    vals: Vec<f64>,
}

impl DofMatrix {
    pub(crate) fn from_entries(n: usize, e: &std::collections::HashMap<(u32, u32), f64>) -> DofMatrix {
        let mut list: Vec<((u32, u32), f64)> = e.iter().map(|(k, v)| (*k, *v)).collect();
        list.sort_unstable_by_key(|x| x.0);
        let mut rows = vec![0u32; n + 1];
        for ((r, _), _) in &list {
            rows[*r as usize + 1] += 1;
        }
        for i in 0..n {
            rows[i + 1] += rows[i];
        }
        DofMatrix { rows, cols: list.iter().map(|x| x.0.1).collect(), vals: list.iter().map(|x| x.1).collect() }
    }

    pub(crate) fn entries(&self) -> impl Iterator<Item = (usize, usize, f64)> + '_ {
        (0..self.rows.len() - 1).flat_map(move |r| (self.rows[r] as usize..self.rows[r + 1] as usize).map(move |e| (r, self.cols[e] as usize, self.vals[e])))
    }

    /// y += E x.
    fn apply_add(&self, x: &[f64], y: &mut [f64]) {
        y.par_iter_mut().with_min_len(GRAIN).enumerate().for_each(|(r, out)| {
            for e in self.rows[r] as usize..self.rows[r + 1] as usize {
                *out += self.vals[e] * x[self.cols[e] as usize];
            }
        });
    }

    /// P^T E P (P a node's, each of its dofs alike): exact, as the cells'
    /// Galerkin products are.
    pub(crate) fn coarsen(&self, t: &Transfer, ncoarse: usize, nd: usize) -> DofMatrix {
        let mut out = std::collections::HashMap::new();
        for (r, c, v) in self.entries() {
            let (n, m) = (r / nd, c / nd);
            for a in t.rows[n] as usize..t.rows[n + 1] as usize {
                for b in t.rows[m] as usize..t.rows[m + 1] as usize {
                    *out.entry((t.cols[a] * nd as u32 + (r % nd) as u32, t.cols[b] * nd as u32 + (c % nd) as u32)).or_insert(0.0) += t.w[a] * v * t.w[b];
                }
            }
        }
        DofMatrix::from_entries(ncoarse * nd, &out)
    }
}

impl System {
    /// The cells' matrices `k_ref` and `k_cut` are given by rows, and kept
    /// by columns: a cell's product then adds whole columns into every row
    /// at once (`element`). They are not symmetric to the bit (a cut
    /// cell's integration, a coarse cell's Galerkin product), so the
    /// column of an entry is its own, not its row's.
    pub(crate) fn new(ndof: usize, cells: Vec<SimCell>, mut k_ref: Vec<Vec<f64>>, mut k_cut: Vec<Vec<f64>>, nodes: usize, hang: &[(u32, Vec<(u32, f64)>)], scale: i64) -> System {
        let ne = 8 * ndof;
        k_ref.par_iter_mut().chain(k_cut.par_iter_mut()).for_each(|k| transpose(k, ne));
        let (colors, runs) = color_cells(&cells);
        let mut hang_of = if hang.is_empty() { Vec::new() } else { vec![u32::MAX; nodes] };
        for (i, (n, _)) in hang.iter().enumerate() {
            hang_of[*n as usize] = i as u32;
        }
        let hang = hang.iter().map(|h| h.1.clone()).collect();
        System { ndof, ne: 8 * ndof, cells, colors, k_ref, k_cut, ndofs: nodes * ndof, shift: Vec::new(), lowrank: Vec::new(), runs, hang_of, hang, scale, across: None }
    }

    /// Whether a corner of the cell hangs.
    pub(crate) fn hangs(&self, c: &SimCell) -> bool {
        !self.hang_of.is_empty() && c.nodes.iter().any(|&n| self.hang_of[n as usize] != u32::MAX)
    }

    /// A cell's dofs with their weights: each corner's, or, where it
    /// hangs, the nodes it hangs on (`ne` lists).
    pub(crate) fn weighted_dofs(&self, c: &SimCell) -> Vec<Vec<(u32, f64)>> {
        let nd = self.ndof;
        (0..self.ne)
            .map(|i| {
                let (n, d) = (c.nodes[i / nd], (i % nd) as u32);
                match self.hang_of.get(n as usize) {
                    Some(&h) if h != u32::MAX => self.hang[h as usize].iter().map(|&(m, w)| (m * nd as u32 + d, w)).collect(),
                    _ => vec![(n * nd as u32 + d, 1.0)],
                }
            })
            .collect()
    }

    /// v += C^T v on hanging nodes' dofs (their share to the nodes they
    /// hang on), which are then 0: a vector of the nodes made one of the
    /// system's unknowns (loads, masses, low-rank terms).
    pub(crate) fn fold(&self, v: &mut [f64], nd: usize) {
        for (n, &h) in self.hang_of.iter().enumerate().filter(|(_, h)| **h != u32::MAX) {
            for d in 0..nd {
                let x = std::mem::take(&mut v[n * nd + d]);
                for &(m, w) in &self.hang[h as usize] {
                    v[m as usize * nd + d] += w * x;
                }
            }
        }
    }

    /// Low-rank terms of the nodes, folded (`fold`).
    pub(crate) fn folded(&self, terms: Vec<(f64, Vec<f64>)>) -> Vec<(f64, Vec<f64>)> {
        terms
            .into_iter()
            .map(|(w, mut g)| {
                self.fold(&mut g, self.ndof);
                (w, g)
            })
            .collect()
    }

    /// Hanging nodes' values from those they hang on (a solution's).
    pub(crate) fn expand(&self, v: &mut [f64], nd: usize) {
        for (n, &h) in self.hang_of.iter().enumerate().filter(|(_, h)| **h != u32::MAX) {
            for d in 0..nd {
                v[n * nd + d] = self.hang[h as usize].iter().map(|&(m, w)| w * v[m as usize * nd + d]).sum();
            }
        }
    }

    /// ye = s K ue for a cell.
    pub(crate) fn element(&self, c: &SimCell, ue: &[f64], ye: &mut [f64]) {
        let mut u = [[0.0; 24]];
        u[0][..ue.len()].copy_from_slice(ue);
        let mut y = [[0.0; 24]];
        self.elements(c, &u, &mut y);
        ye.copy_from_slice(&y[0][..ye.len()]);
    }

    /// ye = s K ue for a cell, for each of M vectors, its matrix read once
    /// for them all. Each row's sum is taken in order from -0.0, as the
    /// row's own `sum` would be (the same bits), but the rows side by side,
    /// a column at a time: a row alone is a chain of additions each waiting
    /// on the last, and its entries lie a row apart.
    fn elements<const M: usize>(&self, c: &SimCell, ue: &[[f64; 24]; M], ye: &mut [[f64; 24]; M]) {
        let (k, s) = self.matrix(c);
        match self.ne {
            24 => by_columns::<24, M>(k, s, ue, ye),
            8 => by_columns::<8, M>(k, s, ue, ye),
            ne => {
                for (ue, ye) in ue.iter().zip(ye) {
                    for (i, y) in ye[..ne].iter_mut().enumerate() {
                        *y = s * (0..ne).map(|j| k[j * ne + i] * ue[j]).sum::<f64>();
                    }
                }
            }
        }
    }

    pub(crate) fn dofs(&self, c: &SimCell) -> [u32; 24] {
        let mut d = [0u32; 24];
        for (i, &n) in c.nodes.iter().enumerate() {
            for k in 0..self.ndof {
                d[i * self.ndof + k] = n * self.ndof as u32 + k as u32;
            }
        }
        d
    }

    /// y = K x (all dofs, unconstrained).
    pub(crate) fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.applies([x], [y]);
    }

    /// y = K x for each of M pairs, as `apply` does each (the same bits),
    /// in one pass over the cells: each cell's matrix and nodes read once
    /// for them all.
    pub(crate) fn applies<const M: usize>(&self, x: [&[f64]; M], mut y: [&mut [f64]; M]) {
        let ne = self.ne;
        for y in y.iter_mut() {
            y.par_iter_mut().with_min_len(GRAIN).for_each(|v| *v = 0.0);
        }
        let yp: [Ptr; M] = std::array::from_fn(|v| Ptr(y[v].as_mut_ptr()));
        for color in &self.colors {
            color.par_iter().with_min_len(CELL_GRAIN).for_each(|&ci| {
                let c = &self.cells[ci as usize];
                let d = self.dofs(c);
                let ue: [[f64; 24]; M] = std::array::from_fn(|v| std::array::from_fn(|i| if i < ne { x[v][d[i] as usize] } else { 0.0 }));
                let mut ye = [[0.0; 24]; M];
                self.elements(c, &ue, &mut ye);
                for (p, ye) in yp.iter().zip(&ye) {
                    for i in 0..ne {
                        // SAFETY: cells of one colour share no nodes, so no two
                        // threads write the same dof in this phase.
                        unsafe { *p.0.add(d[i] as usize) += ye[i] };
                    }
                }
            });
        }
        for runs in &self.runs {
            runs.par_iter().for_each(|&(a, b)| {
                for c in &self.cells[a as usize..b as usize] {
                    let mut ye = [[0.0; 24]; M];
                    if !self.hangs(c) {
                        let d = self.dofs(c);
                        let ue: [[f64; 24]; M] = std::array::from_fn(|v| std::array::from_fn(|i| if i < ne { x[v][d[i] as usize] } else { 0.0 }));
                        self.elements(c, &ue, &mut ye);
                        for (p, ye) in yp.iter().zip(&ye) {
                            for i in 0..ne {
                                // SAFETY: as below.
                                unsafe { *p.0.add(d[i] as usize) += ye[i] };
                            }
                        }
                        continue;
                    }
                    let dofs = self.weighted_dofs(c);
                    let ue: [[f64; 24]; M] = std::array::from_fn(|v| std::array::from_fn(|i| dofs.get(i).map_or(0.0, |l| l.iter().map(|&(g, w)| w * x[v][g as usize]).sum())));
                    self.elements(c, &ue, &mut ye);
                    for (p, ye) in yp.iter().zip(&ye) {
                        for (l, y) in dofs.iter().zip(ye) {
                            for &(g, w) in l {
                                // SAFETY: runs of one colour are grid cells
                                // that share no nodes, nor any their children
                                // hang on (their own corners).
                                unsafe { *p.0.add(g as usize) += w * y };
                            }
                        }
                    }
                }
            });
        }
        for (x, y) in x.into_iter().zip(y) {
            if let Some(e) = &self.across {
                e.apply_add(x, y);
            }
            if !self.shift.is_empty() {
                y.par_iter_mut().with_min_len(GRAIN).zip(x).zip(&self.shift).for_each(|((y, x), s)| *y += s * x);
            }
            for (w, g) in &self.lowrank {
                let c = w * dot(g, x);
                y.par_iter_mut().with_min_len(GRAIN).zip(g).for_each(|(y, g)| *y += c * g);
            }
        }
    }

    /// The cells' diagonal alone (no shift, no low-rank terms).
    fn cell_diagonal(&self) -> Vec<f64> {
        let mut diag = vec![0.0; self.ndofs];
        let ne = self.ne;
        for c in &self.cells {
            let (k, s) = self.matrix(c);
            if !self.hangs(c) {
                let d = self.dofs(c);
                for i in 0..ne {
                    diag[d[i] as usize] += s * k[i * ne + i];
                }
                continue;
            }
            let dofs = self.weighted_dofs(c);
            for (i, a) in dofs.iter().enumerate() {
                for (j, b) in dofs.iter().enumerate() {
                    for &(g, w) in a {
                        for &(g2, w2) in b {
                            if g == g2 {
                                diag[g as usize] += s * w * w2 * k[j * ne + i];
                            }
                        }
                    }
                }
            }
        }
        if let Some(e) = &self.across {
            for (r, _, v) in e.entries().filter(|e| e.0 == e.1) {
                diag[r] += v;
            }
        }
        diag
    }

    /// Holds free rigid motions at rest (mechanisms the loads do not
    /// drive): appends to the low-rank terms one per motion g,
    /// w h h^T with h = D g, D the cells' diagonal, and w = 1 / (g^T D g).
    /// The stiffness gives g no energy, so the term only picks, of the
    /// solutions that differ by g, the one with h^T u = 0; on g the system
    /// acts as its own diagonal (A g = D g), and the term's diagonal, a
    /// share D_i g_i^2 / g^T D g of each D_i, leaves the smoother's alone
    /// however far D ranges (a redesign's void cells, at 1e-9 of solid,
    /// say). A uniform weight (the mean diagonal over g^T g) swamped the
    /// void nodes' diagonal 1e5 times over, so the smoother left their
    /// errors alone: an MBB beam free to slide out of its plane stalled at
    /// a residual of 1.3e-5. Called again when the cells' stiffness changes.
    pub(crate) fn hold(&mut self, motions: &[Vec<f64>]) {
        let d = self.cell_diagonal();
        for g in motions {
            let h: Vec<f64> = d.iter().zip(g).map(|(d, g)| d * g).collect();
            self.lowrank.push((1.0 / dot(g, &h), h));
        }
    }

    pub(crate) fn diagonal(&self) -> Vec<f64> {
        let mut diag = self.cell_diagonal();
        for (d, s) in diag.iter_mut().zip(&self.shift) {
            *d += s;
        }
        for (w, g) in &self.lowrank {
            for (d, g) in diag.iter_mut().zip(g) {
                *d += w * g * g;
            }
        }
        diag
    }
}

impl System {
    /// Element matrix of a cell (by columns) and the factor it is scaled by.
    pub(crate) fn matrix(&self, c: &SimCell) -> (&[f64], f64) {
        if c.cut == u32::MAX { (&self.k_ref[c.mat as usize], c.coef * c.mult) } else { (&self.k_cut[c.cut as usize], c.mult) }
    }
}

/// A square matrix of `n` rows, by rows, in place by columns.
fn transpose(k: &mut [f64], n: usize) {
    for i in 0..n {
        for j in i + 1..n {
            k.swap(i * n + j, j * n + i);
        }
    }
}

/// `System::elements` for a matrix of NE columns.
fn by_columns<const NE: usize, const M: usize>(k: &[f64], s: f64, ue: &[[f64; 24]; M], ye: &mut [[f64; 24]; M]) {
    let columns: &[[f64; NE]; NE] = k.as_chunks::<NE>().0[..NE].try_into().expect("a cell's matrix");
    let mut sum = [[-0.0; NE]; M];
    for (j, column) in columns.iter().enumerate() {
        for (sum, ue) in sum.iter_mut().zip(ue) {
            let u = ue[j];
            for (a, k) in sum.iter_mut().zip(column) {
                *a += k * u;
            }
        }
    }
    for (ye, sum) in ye.iter_mut().zip(sum) {
        for (y, a) in ye[..NE].iter_mut().zip(sum) {
            *y = s * a;
        }
    }
}
