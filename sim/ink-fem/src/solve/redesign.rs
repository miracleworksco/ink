//! A redesign (`solve_redesigned`): the assembled static study handed to
//! it as an `Assembly`, solved as often as it likes.

use crate::cells::{ALPHA, place};
use crate::element::shape_n;
use crate::mechanisms::rigid_mode;
use crate::multigrid::{Multigrid, pcg};
use crate::numerics::dot;
use crate::occupancy::at_xi;
use crate::octree::Grid;
use crate::system::{DofMatrix, System};
use crate::{Fix, MAX_ITER};
use ink_core::{V3, grid, v3};
use rayon::prelude::*;

use super::assemble::{Mesh, Model};
use super::freeing::{Freed, mean_material};
use super::holds::Clamp;
use super::run::{Rhs, Run};
use super::study::Study;
use super::{Assembly, PENALTY, SolveError};

// The assembled study, as `Assembly` hands it to the redesign.
struct Asm<'a> {
    sys: &'a mut System,
    x: &'a mut [f64],
    freed: Option<Vec<Freed>>,
    iters: usize,
    res: f64,
    h: f64,
    grid: &'a Grid,
    dims: [i64; 3],
    cell_at: &'a [u32],
    f: &'a [f64],
    f_loads: &'a [f64],
    holds: &'a [Vec<f64>],
    final_idx: &'a [[i64; 3]],
    final_pos: &'a [V3],
    clamps: &'a [Clamp],
    clamp_cells: &'a [(usize, usize)],
    fixes: &'a [Fix],
    freeing: &'a dyn Fn(&System, &[Clamp]) -> Result<(Vec<Freed>, Vec<(f64, Vec<f64>)>), String>,
    across: &'a dyn Fn(&System) -> Option<DofMatrix>,
    lent_energies: &'a dyn Fn(&System, &[f64], &[Option<Vec<f64>>], &mut [f64]),
}
impl Assembly for Asm<'_> {
    fn cells(&self) -> usize {
        self.sys.cells.len()
    }
    fn h(&self) -> f64 {
        self.h
    }
    fn cell(&self, i: usize) -> ([i64; 3], V3, f64) {
        let c = &self.sys.cells[i];
        let (ctr, size) = place(self.grid, c);
        (c.idx, ctr, c.frac * size * size * size)
    }
    fn cells_at(&self, idx: [i64; 3]) -> std::ops::Range<usize> {
        if (0..3).any(|a| idx[a] < 0 || idx[a] >= self.dims[a]) {
            return 0..0;
        }
        match self.cell_at[grid::index(idx.map(|v| v as usize), self.dims.map(|v| v as usize))] {
            u32::MAX => 0..0,
            k => {
                let k = k as usize;
                k..k + self.sys.cells[k..].iter().take_while(|c| c.idx == idx).count()
            }
        }
    }
    fn held(&self) -> Vec<bool> {
        let mut out = vec![false; self.sys.cells.len()];
        for &(cell, k) in self.clamp_cells {
            out[cell] |= self.clamps[k].whole;
        }
        out
    }
    fn solve(&mut self, mult: &[f64], tol: f64) -> Result<f64, String> {
        for (c, m) in self.sys.cells.iter_mut().zip(mult) {
            c.mult = *m;
        }
        // The supports' terms at this stiffness, the lent tractions too.
        self.sys.across = (self.across)(self.sys);
        let (fr, lr) = (self.freeing)(self.sys, self.clamps)?;
        self.sys.lowrank = self.sys.folded(lr);
        self.sys.hold(self.holds);
        let sys = &*self.sys;
        let mg = Multigrid::build(sys, self.final_idx);
        let pre = |r: &[f64]| mg.vcycle(0, sys, r);
        let (its, res) = pcg(sys, &[], &pre, self.f, self.x, tol, MAX_ITER);
        sys.expand(self.x, 3);
        (self.freed, self.iters, self.res) = (Some(fr), self.iters + its, res);
        Ok(dot(self.f_loads, self.x))
    }
    fn energies(&self) -> Vec<f64> {
        let (sys, x) = (&*self.sys, &*self.x);
        let mut energy: Vec<f64> = sys
            .cells
            .par_iter()
            .map(|c| {
                let d = sys.dofs(c);
                let ne = sys.ne;
                let mut ue = [0.0; 24];
                let mut ye = [0.0; 24];
                for i in 0..ne {
                    ue[i] = x[d[i] as usize];
                }
                sys.element(c, &ue[..ne], &mut ye[..ne]);
                (0..ne).map(|i| ue[i] * ye[i]).sum::<f64>() / c.mult.max(1e-300)
            })
            .collect();
        let mut turns = vec![None; self.fixes.len()];
        for fr in self.freed.iter().flatten() {
            let k = fr.b.len();
            let fx = &self.fixes[fr.fix];
            let btu: Vec<f64> = fr.b.iter().zip(&fr.moved).map(|(b, m)| dot(b, x) - m).collect();
            let q: Vec<f64> = (0..k).map(|i| (0..k).map(|j| fr.cpinv[i * k + j] * btu[j]).sum()).collect();
            if let Some(mean) = &fr.mean {
                // The distributing hold is the difference of two
                // best-fit energies: all rigid motions, then just
                // the free ones. At each optimum the change in the
                // fitted motion cancels. For either fit w its weight
                // contributes 2 u^T W w - w^T W w.
                let mut coeff = [0.0; 6];
                for ((_, g), b) in mean.all.iter().zip(&mean.basis) {
                    let a = dot(g, x);
                    for m in 0..6 {
                        coeff[m] += a * b[m];
                    }
                }
                let mut balance_grad = 0.0;
                let mine = || self.clamp_cells.iter().filter(|c| self.clamps[c.1].fix == fr.fix);
                for &(cell, kc) in mine() {
                    let cl = &self.clamps[kc];
                    let nodes = sys.cells[cell].nodes;
                    let all: [V3; 8] = nodes.map(|n| {
                        let at = self.final_pos[n as usize];
                        (0..6).fold([0.0; 3], |v, m| v3::add(v, v3::scale(rigid_mode(at, mean.center, mean.size, m), coeff[m])))
                    });
                    let free: [V3; 8] = nodes.map(|n| fx.free.iter().zip(&q).fold([0.0; 3], |v, (m, q)| v3::add(v, v3::scale(m.at(self.final_pos[n as usize]), *q))));
                    let material = mean_material(cl, mean.even);
                    let spring = cl.spring.filter(|(beta, _)| mean.even < beta * PENALTY / ALPHA);
                    for i in 0..24 {
                        let (a, d) = (i / 3, i % 3);
                        let u = x[nodes[a] as usize * 3 + d];
                        for j in 0..24 {
                            let (b, e) = (j / 3, j % 3);
                            let work = (2.0 * u - all[a][d]) * all[b][e] - (2.0 * u - free[a][d]) * free[b][e];
                            energy[cell] += cl.t.kww()[i * 24 + j] * work;
                            if d == e && fx.dofs[d] {
                                let ij = a * 8 + b;
                                energy[cell] += material[ij] * work;
                                if let Some((_, m)) = &spring {
                                    balance_grad += sys.cells[cell].mult * m[ij] * work;
                                }
                            }
                        }
                    }
                }
                if mean.spring_total > 0.0 {
                    // Differentiate the boundary-to-interior weight ratio.
                    for &(cell, kc) in mine() {
                        let cl = &self.clamps[kc];
                        let spring = cl.spring.map_or(0.0, |(_, m)| m.iter().sum());
                        energy[cell] += balance_grad * (cl.t.mass.iter().sum::<f64>() - mean.even * spring) / mean.spring_total;
                    }
                }
                continue;
            }
            // A rigid hold's energy uses its fitted turn w = Phi q:
            // u^T uu u + 2 u^T uw w + w^T ww w.
            for &(cell, kc) in self.clamp_cells.iter().filter(|c| self.clamps[c.1].fix == fr.fix) {
                let cl = &self.clamps[kc];
                let (uw, ww, nodes) = (&cl.t.kuw(), cl.t.kww(), sys.cells[cell].nodes);
                let ue: Vec<f64> = (0..24).map(|i| x[nodes[i / 3] as usize * 3 + i % 3]).collect();
                let we: Vec<f64> = (0..24).map(|i| (0..k).map(|m| q[m] * fx.free[m].at(self.final_pos[nodes[i / 3] as usize])[i % 3]).sum::<f64>()).collect();
                let form = |a: &[f64], v: &[f64], w: &[f64]| (0..24).map(|i| v[i] * (0..24).map(|j| a[i * 24 + j] * w[j]).sum::<f64>()).sum::<f64>();
                energy[cell] += 2.0 * form(uw, &ue, &we) + form(ww, &we, &we);
            }
            turns[fr.fix] = Some(q);
        }
        (self.lent_energies)(sys, x, &turns, &mut energy);
        energy
    }
}

/// Runs the redesign on the assembled study: its solves' iterations and
/// residual, and its last solve's free motions.
#[allow(clippy::too_many_arguments)]
pub(super) fn run_redesign(
    model: &Model,
    sys: &mut System,
    run: &mut Run,
    rhs: &Rhs,
    holds: &[Vec<f64>],
    freed: &mut Vec<Freed>,
    redesign: &mut dyn FnMut(&mut dyn Assembly) -> Result<(), String>,
) -> Result<(), SolveError> {
    let Model { st, cells, lent, mesh, .. } = *model;
    let Study { s, ndof, h, dims, ref grid, .. } = *st;
    let classified = &cells.classified;
    let Mesh { ref final_of, ref final_idx, ref final_pos, ref cell_at, ref clamp_cells, .. } = *mesh;
    let Run { ref mut x, ref mut iters, ref mut res, .. } = *run;
    let Rhs { ref f, ref f_loads, .. } = *rhs;
    let clamps = &model.holds.clamps;
    let freeing = |sys: &System, clamps: &[Clamp]| model.freeing(sys, clamps);
    let across_of = |sys: &System| model.across_of(sys);
    if s.modes > 0 || s.transient.is_some() || ndof != 3 {
        return Err("a redesign applies to static elastic studies".into());
    }
    let lent_energies = |sys: &System, x: &[f64], turns: &[Option<Vec<f64>>], energy: &mut [f64]| {
        for l in lent.iter() {
            let [from, to] = [l.from, l.to].map(|i| final_of[i] as usize);
            let (piece, lender) = (&sys.cells[from], &sys.cells[to]);
            let (ctr, size) = place(grid, &classified[l.from]);
            let at = at_xi(ctr, size, l.t.xi);
            let mut gap = [0.0; 3];
            for (i, &n) in piece.nodes.iter().enumerate() {
                let w = shape_n(i, l.t.xi);
                for d in 0..3 {
                    gap[d] += w * x[n as usize * 3 + d];
                }
            }
            if let Some(q) = &turns[l.fix] {
                for (m, q) in s.fixes[l.fix].free.iter().zip(q) {
                    gap = v3::sub(gap, v3::scale(m.at(at), *q));
                }
            }
            let traction = |d: usize| (0..24).map(|j| l.pf[d * 24 + j] * x[lender.nodes[j / 3] as usize * 3 + j % 3]).sum::<f64>();
            let work = -model.lent_scale(sys, l) * l.w * (0..3).map(|d| gap[d] * traction(d)).sum::<f64>();
            // The two transposed terms contribute 2 * work. Each
            // multiplier differentiates half their geometric mean.
            energy[from] += work / piece.mult;
            energy[to] += work / lender.mult;
        }
    };
    let mut a = Asm {
        sys,
        x,
        freed: None,
        iters: 0,
        res: 0.0,
        h,
        grid,
        dims,
        cell_at,
        f,
        f_loads,
        holds,
        final_idx,
        final_pos,
        clamps,
        clamp_cells,
        fixes: &s.fixes,
        freeing: &freeing,
        across: &across_of,
        lent_energies: &lent_energies,
    };
    redesign(&mut a)?;
    let (fr, it, r) = (a.freed.take(), a.iters, a.res);
    (*iters, *res) = (it, r);
    if let Some(fr) = fr {
        *freed = fr;
    }
    Ok(())
}
