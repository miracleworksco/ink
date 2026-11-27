//! Supports with free rigid motions (`free_supports`): what each holds
//! of its region's motion, as low-rank terms of the system.

use crate::GRAIN;
use crate::cells::{ALPHA, place};
use crate::eigen::dense_sym_eigen;
use crate::mechanisms::rigid_mode;
use crate::numerics::dot;
use crate::occupancy::at_xi;
use crate::system::System;
use ink_core::{V3, v3};
use rayon::prelude::*;

use super::assemble::{Mesh, Model};
use super::holds::Clamp;
use super::study::Study;
use super::{PENALTY, SolveError};

// Supports with free rigid motions. On the region what the fix holds is
// held relative to the motions' combination Phi q, its penalty
// (u - Phi q)^T A (u - Phi q) (A the held penalty). The amplitudes are
// its best fit, q = C^+ B^T u with B = A Phi and C = Phi^T A Phi (Phi
// at the nodes: trilinear cells carry a rigid motion exactly), which
// leaves A - B C^+ B^T: a low-rank term -g g^T per independent motion,
// g = B v / sqrt(l) for each eigenpair (l, v) of C.
// A fix that moves its region (a held value v) holds (u - Phi q - v):
// its right-hand side A v (the samples' sum) loses B C^+ Phi^T A v,
// and q = C^+ (B^T u - Phi^T A v).
// Each cell's penalty is scaled with its material (`mult`: a
// redesign's), so these terms are rebuilt as it changes.
//
// A distributing hold (`distributing`) holds, of its region's rigid
// motions R (three slides, three turns about its boundary's centroid),
// those its free motions leave, S (R's W-orthogonal complement of
// them, within R; W = ww, its boundary's weight): its energy is
// (u - v)^T W P (u - v) / 2, P = S (S^T W S)^+ S^T W the W-projection
// onto S, a positive term g g^T per motion of an W-orthonormal basis of
// S, g = W s, and its right-hand side W P v. Its amplitudes, reported,
// are the free motions' share of the W-projection of u - v onto R (B =
// W Phi, C = Phi^T W Phi).
pub(super) struct MeanHold {
    pub(super) held: Vec<(Vec<f64>, f64)>,
    pub(super) all: Vec<(f64, Vec<f64>)>,
    pub(super) basis: Vec<[f64; 6]>,
    pub(super) center: V3,
    pub(super) size: f64,
    pub(super) spring_total: f64,
    pub(super) even: f64,
}

// A distributing hold balances its interior weight with its boundary
// weight, capped by the penalty it would use for a rigid hold.
pub(super) fn mean_material(cl: &Clamp, even: f64) -> [f64; 64] {
    cl.spring.map_or([0.0; 64], |(beta, m)| m.map(|v| v * even.min(beta * PENALTY / ALPHA)))
}

pub(super) struct Freed {
    pub(super) fix: usize,
    pub(super) b: Vec<Vec<f64>>,
    pub(super) cpinv: Vec<f64>,
    pub(super) terms: std::ops::Range<usize>,
    /// Phi^T A v: the held values' work along each free motion.
    pub(super) moved: Vec<f64>,
    /// A distributing hold's: each of its terms' g (unfolded) and g^T v,
    /// and the terms that hold its free motions too (a mirror's odd
    /// sector, where they are zero).
    pub(super) mean: Option<MeanHold>,
}

impl Model<'_> {
    /// Each support with free motions' terms (`Freed`), at the cells' scales.
    pub(super) fn freeing(&self, sys: &System, clamps: &[Clamp]) -> Result<(Vec<Freed>, Vec<(f64, Vec<f64>)>), String> {
        let Model { st, cells, lent, mesh, .. } = *self;
        let Study { s, ne, h, ref grid, ref distributing, .. } = *st;
        let classified = &cells.classified;
        let Mesh { ref final_of, ref final_pos, ref clamp_cells, .. } = *mesh;
        let (mut freed, mut terms) = (Vec::new(), Vec::new());
        for (fi, fx) in s.fixes.iter().enumerate().filter(|(_, f)| !f.free.is_empty()) {
            let k = fx.free.len();
            if distributing[fi] {
                let mine: Vec<(usize, usize)> = clamp_cells.iter().filter(|cc| clamps[cc.1].fix == fi).copied().collect();
                // Its weight: its boundary's (ww, of mass `mass`), and its
                // material's (int N_i N_j over it, `Clamp::spring`'s, on its
                // axes), as much in all as its boundary's, and per volume no
                // more than a penalty's (PENALTY): a region's material moves
                // as it does. (Weighed by its boundary alone, a stop over a
                // lever's tip, held along y on its face across the lever,
                // could not tell its turn about a line in it from a slide.
                // Scaled to its boundary's alone, a slab 1e-6 thick on a
                // column's end weighed a million times a cell's, its solve
                // took 55 iterations for 16 and its buckling factor read 8%
                // low; at a penalty's alone, the stop's took 255 for 21.)
                let total = |f: &dyn Fn(&Clamp) -> f64| mine.iter().map(|&(cell, ki)| sys.cells[cell].mult * f(&clamps[ki])).sum::<f64>();
                let spring_total = total(&|cl| cl.spring.map_or(0.0, |x| x.1.iter().sum()));
                let even = if spring_total > 0.0 { total(&|cl| cl.t.mass.iter().sum()) / spring_total } else { 0.0 };
                let material = |cl: &Clamp| mean_material(cl, even);
                let weight = |cl: &Clamp| -> (Vec<f64>, [f64; 64]) {
                    let (m, mut ww, mut mass) = (material(cl), cl.t.kww().to_vec(), cl.t.mass);
                    for ij in 0..64 {
                        mass[ij] += m[ij];
                        for d in (0..3).filter(|&d| fx.dofs[d]) {
                            ww[((ij / 8) * 3 + d) * ne + (ij % 8) * 3 + d] += m[ij];
                        }
                    }
                    (ww, mass)
                };
                // Its boundary's centroid and spread, by its weight.
                let (mut wsum, mut cen, mut spread) = (0.0, [0.0; 3], 0.0);
                for pass in 0..2 {
                    for &(cell, ki) in &mine {
                        let (nodes, mult) = (sys.cells[cell].nodes, sys.cells[cell].mult);
                        for i in 0..8 {
                            let row = mult * weight(&clamps[ki]).1[i * 8..i * 8 + 8].iter().sum::<f64>();
                            let at = final_pos[nodes[i] as usize];
                            if pass == 0 {
                                wsum += row;
                                (0..3).for_each(|a| cen[a] += row * at[a]);
                            } else {
                                spread += row * (0..3).map(|a| (at[a] - cen[a]).powi(2)).sum::<f64>();
                            }
                        }
                    }
                    if pass == 0 {
                        if !(wsum > 0.0) {
                            return Err(format!("fix '{}' holds nothing: its region's boundary with the part's free material is all other fixes'", fx.name));
                        }
                        cen = cen.map(|v| v / wsum);
                    }
                }
                let size = (spread / wsum).sqrt().max(1e-6 * h);
                // W R, R^T W R, R^T W v; and each free motion's energy were the
                // region held on every axis.
                let (mut wr, mut g6, mut rv, mut c_all) = (vec![vec![0.0; sys.ndofs]; 6], [0.0; 36], [0.0; 6], vec![0.0; k]);
                for &(cell, ki) in &mine {
                    let t = &clamps[ki].t;
                    let (ww, mass) = weight(&clamps[ki]);
                    let (nodes, mult) = (sys.cells[cell].nodes, sys.cells[cell].mult);
                    // (its held value's work: its boundary's, and its material's at the nodes)
                    let held = st.values_at(fx, &nodes.map(|n| final_pos[n as usize]));
                    let m = material(&clamps[ki]);
                    let mut gw = t.gw.clone();
                    for ij in 0..64 {
                        for d in (0..3).filter(|&d| fx.dofs[d]) {
                            gw[(ij / 8) * 3 + d] += m[ij] * held[ij % 8][d];
                        }
                    }
                    let r: Vec<Vec<f64>> = (0..6).map(|m| (0..ne).map(|i| rigid_mode(final_pos[nodes[i / 3] as usize], cen, size, m)[i % 3]).collect()).collect();
                    for m in 0..6 {
                        let wph: Vec<f64> = (0..ne).map(|i| (0..ne).map(|j| ww[i * ne + j] * r[m][j]).sum()).collect();
                        for (i, v) in wph.iter().enumerate() {
                            wr[m][nodes[i / 3] as usize * 3 + i % 3] += mult * v;
                        }
                        rv[m] += mult * dot(&r[m], &gw);
                        for m2 in 0..6 {
                            g6[m * 6 + m2] += mult * dot(&r[m2], &wph);
                        }
                    }
                    for (q, mo) in fx.free.iter().enumerate() {
                        let ph: Vec<V3> = nodes.iter().map(|&n| mo.at(final_pos[n as usize])).collect();
                        c_all[q] += mult * (0..64).map(|ij| mass[ij] * v3::dot(ph[ij / 8], ph[ij % 8])).sum::<f64>();
                    }
                }
                // Each free motion in R's amplitudes (a turn about the
                // centroid, per `size`, and a slide).
                let a: Vec<[f64; 6]> = fx
                    .free
                    .iter()
                    .map(|mo| {
                        let t = mo.at(cen);
                        [t[0], t[1], t[2], mo.omega[0] * size, mo.omega[1] * size, mo.omega[2] * size]
                    })
                    .collect();
                let gram = |x: &[f64; 6], y: &[f64; 6]| (0..36).map(|ij| x[ij / 6] * g6[ij] * y[ij % 6]).sum::<f64>();
                for (q, mo) in fx.free.iter().enumerate() {
                    if !(gram(&a[q], &a[q]) > 1e-9 * c_all[q]) {
                        return Err(format!("fix '{}': its free {} moves nothing the fix holds (its (dofs ...) leave that motion free already)", fx.name, mo.describe()));
                    }
                }
                // R's W-orthonormal coordinates y (of what W sees of it), the
                // free motions' span in them, and its complement.
                let (vals, vecs) = dense_sym_eigen(&g6, 6);
                let top = vals.iter().fold(0.0f64, |a, &v| a.max(v));
                let kept: Vec<usize> = (0..6).filter(|&c| vals[c] > 1e-12 * top).collect();
                let n = kept.len();
                let to_y = |x: &[f64; 6]| -> Vec<f64> { kept.iter().map(|&c| vals[c].sqrt() * (0..6).map(|m| vecs[m * 6 + c] * x[m]).sum::<f64>()).collect() };
                let from_y = |y: &[f64]| -> [f64; 6] { std::array::from_fn(|m| kept.iter().zip(y).map(|(&c, y)| vecs[m * 6 + c] / vals[c].sqrt() * y).sum()) };
                let fy: Vec<Vec<f64>> = a.iter().map(to_y).collect();
                let mut cf = vec![0.0; k * k];
                for i in 0..k {
                    for j in 0..k {
                        cf[i * k + j] = dot(&fy[i], &fy[j]);
                    }
                }
                let (fv, fvecs) = dense_sym_eigen(&cf, k);
                let ftop = fv.iter().fold(0.0f64, |a, &v| a.max(v));
                let mut cpinv = vec![0.0; k * k];
                let mut perp: Vec<f64> = (0..n * n).map(|ij| if ij / n == ij % n { 1.0 } else { 0.0 }).collect();
                for (col, &l) in fv.iter().enumerate().filter(|(_, l)| **l > 1e-12 * ftop) {
                    let v: Vec<f64> = (0..k).map(|m| fvecs[m * k + col]).collect();
                    for i in 0..k {
                        for j in 0..k {
                            cpinv[i * k + j] += v[i] * v[j] / l;
                        }
                    }
                    // (its unit vector in y)
                    let e: Vec<f64> = (0..n).map(|d| (0..k).map(|m| v[m] * fy[m][d]).sum::<f64>() / l.sqrt()).collect();
                    for ij in 0..n * n {
                        perp[ij] -= e[ij / n] * e[ij % n];
                    }
                }
                let term = |c: [f64; 6]| -> ((f64, Vec<f64>), f64) {
                    let mut g = vec![0.0; sys.ndofs];
                    for (m, wm) in wr.iter().enumerate().filter(|(m, _)| c[*m] != 0.0) {
                        g.par_iter_mut().with_min_len(GRAIN).zip(wm).for_each(|(g, w)| *g += c[m] * w);
                    }
                    ((1.0, g), (0..6).map(|m| c[m] * rv[m]).sum())
                };
                let (pv, pvecs) = if n > 0 { dense_sym_eigen(&perp, n) } else { (Vec::new(), Vec::new()) };
                let first_term = terms.len();
                let mut held = Vec::new();
                for col in (0..n).filter(|&col| pv[col] > 0.5) {
                    let (t, v) = term(from_y(&(0..n).map(|d| pvecs[d * n + col]).collect::<Vec<f64>>()));
                    held.push((t.1.clone(), v));
                    terms.push(t);
                }
                let basis: Vec<[f64; 6]> = (0..n).map(|d| from_y(&(0..n).map(|e| (e == d) as u8 as f64).collect::<Vec<f64>>())).collect();
                let all = sys.folded(basis.iter().map(|&c| term(c).0).collect());
                let b: Vec<Vec<f64>> = a
                    .iter()
                    .map(|am| {
                        let mut g = vec![0.0; sys.ndofs];
                        for (m, wm) in wr.iter().enumerate().filter(|(m, _)| am[*m] != 0.0) {
                            g.par_iter_mut().with_min_len(GRAIN).zip(wm).for_each(|(g, w)| *g += am[m] * w);
                        }
                        g
                    })
                    .collect();
                let moved = a.iter().map(|am| (0..6).map(|m| am[m] * rv[m]).sum()).collect();
                let mean = MeanHold { held, all, basis, center: cen, size, spring_total, even };
                freed.push(Freed { fix: fi, b, cpinv, moved, terms: first_term..terms.len(), mean: Some(mean) });
                continue;
            }
            let mut b = vec![vec![0.0; sys.ndofs]; k];
            // C, and each motion's energy were the region held on every axis.
            let (mut c, mut c_all, mut moved) = (vec![0.0; k * k], vec![0.0; k], vec![0.0; k]);
            // B = -uw Phi and C = Phi^T ww Phi (a penalty's A, beta A Phi
            // and Phi^T beta A Phi).
            for &(cell, ki) in clamp_cells.iter().filter(|cc| clamps[cc.1].fix == fi) {
                let t = &clamps[ki].t;
                let (nodes, mult) = (sys.cells[cell].nodes, sys.cells[cell].mult);
                let (uw, ww, gw, mass) = (&t.kuw(), t.kww(), &t.gw, &t.mass);
                let phi: Vec<Vec<f64>> = fx.free.iter().map(|mo| nodes.iter().flat_map(|&n| mo.at(final_pos[n as usize])).collect()).collect();
                for (m, ph) in phi.iter().enumerate() {
                    moved[m] += mult * ph.iter().zip(gw).map(|(p, g)| p * g).sum::<f64>();
                    let bph: Vec<f64> = (0..ne).map(|r| -(0..ne).map(|q| uw[r * ne + q] * ph[q]).sum::<f64>()).collect();
                    let wph: Vec<f64> = (0..ne).map(|r| (0..ne).map(|q| ww[r * ne + q] * ph[q]).sum()).collect();
                    for (r, v) in bph.iter().enumerate() {
                        b[m][nodes[r / 3] as usize * 3 + r % 3] += mult * v;
                    }
                    for (m2, ph2) in phi.iter().enumerate() {
                        c[m * k + m2] += mult * ph2.iter().zip(&wph).map(|(p, q)| p * q).sum::<f64>();
                    }
                    c_all[m] += mult * (0..64).map(|ij| mass[ij] * (0..3).map(|d| ph[(ij / 8) * 3 + d] * ph[(ij % 8) * 3 + d]).sum::<f64>()).sum::<f64>();
                }
            }
            // (a lent traction's, against the motion at its piece: `Lent`)
            for l in lent.iter().filter(|l| l.fix == fi) {
                let (ctr, size) = place(grid, &classified[l.from]);
                let nodes = sys.cells[final_of[l.to] as usize].nodes;
                for (m, mo) in fx.free.iter().enumerate() {
                    let ph = mo.at(at_xi(ctr, size, l.t.xi));
                    for col in 0..ne {
                        b[m][nodes[col / 3] as usize * 3 + col % 3] -= self.lent_scale(sys, l) * l.w * (0..3).map(|a| l.pf[a * ne + col] * ph[a]).sum::<f64>();
                    }
                }
            }
            for (m, mo) in fx.free.iter().enumerate() {
                if !(c[m * k + m] > 1e-9 * c_all[m]) {
                    return Err(format!("fix '{}': its free {} moves nothing the fix holds (its (dofs ...) leave that motion free already)", fx.name, mo.describe()));
                }
            }
            let (vals, vecs) = dense_sym_eigen(&c, k);
            let top = vals.iter().fold(0.0f64, |a, &v| a.max(v));
            let mut cpinv = vec![0.0; k * k];
            let first_term = terms.len();
            for (col, &l) in vals.iter().enumerate().filter(|(_, l)| **l > 1e-12 * top) {
                let v: Vec<f64> = (0..k).map(|m| vecs[m * k + col]).collect();
                for i in 0..k {
                    for j in 0..k {
                        cpinv[i * k + j] += v[i] * v[j] / l;
                    }
                }
                let mut g = vec![0.0; sys.ndofs];
                for (bm, vm) in b.iter().zip(&v) {
                    g.par_iter_mut().with_min_len(GRAIN).zip(bm).for_each(|(g, b)| *g += b * vm / l.sqrt());
                }
                terms.push((-1.0, g));
            }
            freed.push(Freed { fix: fi, b, cpinv, moved, terms: first_term..terms.len(), mean: None });
        }
        Ok((freed, terms))
    }
}

/// The supports' free motions' terms, folded into the system, and their
/// held values' share taken off the right-hand side.
pub(super) fn free_supports(model: &Model, sys: &mut System, f: &mut [f64]) -> Result<Vec<Freed>, SolveError> {
    let clamps = &model.holds.clamps;
    let (freed, terms) = model.freeing(sys, clamps)?;
    sys.lowrank = sys.folded(terms);
    for fr in &freed {
        if let Some(mean) = &fr.mean {
            for (g, v) in &mean.held {
                f.par_iter_mut().with_min_len(GRAIN).zip(g).for_each(|(f, g)| *f += v * g);
            }
            continue;
        }
        let k = fr.b.len();
        for (i, bi) in fr.b.iter().enumerate() {
            let w: f64 = (0..k).map(|j| fr.cpinv[i * k + j] * fr.moved[j]).sum();
            if w != 0.0 {
                f.iter_mut().zip(bi).for_each(|(f, b)| *f -= w * b);
            }
        }
    }
    Ok(freed)
}
