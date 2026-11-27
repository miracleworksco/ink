//! Step 8, the solves: the lowest modes (`modal`), K u = f by preconditioned
//! CG from the held values lifted onto their nodes (`lift`, `solve_static`).

use crate::eigen::lobpcg;
use crate::multigrid::{Multigrid, pcg};
use crate::numerics::dot;
use crate::system::System;
use crate::{MAX_ITER, MODAL_ITER};
use ink_core::V3;
use rayon::prelude::*;

use super::SolveError;
use super::assemble::{Mesh, Model};
use super::rigid::Rigid;
use super::study::Study;

/// The right-hand side: the loads' and the holds' (`f`), the holds' alone
/// (`f_held`), and the loads' (`f_loads`, their work the compliance).
pub(super) struct Rhs {
    pub(super) f: Vec<f64>,
    pub(super) f_held: Vec<f64>,
    pub(super) f_loads: Vec<f64>,
}

/// A modal study's: its modes (frequency and participation), mass, shapes,
/// and the eigensolver's iterations and convergence.
pub(super) struct Modes {
    pub(super) modes_out: Vec<(f64, [f64; 6])>,
    pub(super) rigid_modes: usize,
    pub(super) mass_kg: f64,
    pub(super) modal_u: Option<Vec<f64>>,
    pub(super) mode_shapes: Vec<Vec<f64>>,
    pub(super) modal_iters: usize,
    pub(super) modal_ok: bool,
}

/// What the solver gives: the solution (`x`), its iterations and residual,
/// and over time, its history, peak and time levels.
pub(super) struct Run {
    pub(super) x: Vec<f64>,
    pub(super) iters: usize,
    pub(super) res: f64,
    pub(super) history: Vec<(f64, f64)>,
    pub(super) peak_time: Option<f64>,
    pub(super) peak: Option<(f64, f64, V3)>,
    pub(super) time_levels: Vec<(usize, f64, f64)>,
}

impl Run {
    pub(super) fn new(ndofs: usize) -> Self {
        let x = vec![0.0; ndofs];
        let (iters, res) = (0, 0.0);
        let (history, peak_time, peak, time_levels) = (Vec::new(), None, None, Vec::new());
        Run { x, iters, res, history, peak_time, peak, time_levels }
    }
}

/// The lowest modes, K x = w^2 M x (where the study asks for them), of
/// the motions mass-orthogonal to the rigid motions the supports leave
/// free (`rigid`: each a mode at 0 Hz, no vibration).
pub(super) fn modal(model: &Model, sys: &System, mg: &Multigrid, rigid: &Rigid) -> Result<Modes, SolveError> {
    let Study { s, ndof, .. } = *model.st;
    let Mesh { nn, ref final_pos, .. } = *model.mesh;
    let (mut modes_out, mut mass_kg) = (Vec::new(), 0.0);
    let mut modal_u: Option<Vec<f64>> = None;
    let mut mode_shapes: Vec<Vec<f64>>;
    let mut modal_iters = 0usize;
    let mut modal_ok = true;
    let mut rigid_modes = 0;
    mode_shapes = Vec::new();
    if s.modes > 0 {
        let rho = s.density.as_ref().ok_or("modal analysis needs (density ...) in (elastic ...)")?;
        // Lumped mass per node: integral of rho N_i, in tonnes (consistent with N, mm, s).
        let parts: Vec<([f64; 8], f64)> = sys.cells.par_iter().map(|c| model.sample_cell(c, None, Some(rho))).collect();
        let share = model.held_share(&sys.cells);
        let (mut m, mut md) = (vec![0.0; nn], vec![0.0; sys.ndofs]);
        for (ci, (c, (nodal, _))) in sys.cells.iter().zip(&parts).enumerate() {
            let held = share.get(&ci);
            for k in 0..8 {
                m[c.nodes[k] as usize] += nodal[k] * 1e-12;
                for d in 0..ndof {
                    md[c.nodes[k] as usize * ndof + d] += nodal[k] * 1e-12 * (1.0 - held.map_or(0.0, |h| h[k][d]));
                }
            }
        }
        sys.fold(&mut m, 1);
        sys.fold(&mut md, ndof);
        mass_kg = m.iter().sum::<f64>() * 1e3;
        let floor = 1e-8 * m.iter().sum::<f64>() / nn as f64;
        let mdiag: Vec<f64> = md.iter().map(|m| m.max(floor)).collect();
        // The free rigid motions, on the dofs that are unknowns (a hanging
        // node's are its masters'), mass-orthonormal: Y. Every vector the
        // eigensolve takes in (a preconditioned residual, a start) has its
        // share of them taken out, so its modes are those mass-orthogonal
        // to them: a free part's elastic modes, its rigid motions not among
        // them. The stiffness gives a rigid motion no energy (a support's
        // free motion under 1e-9 of its hold's), so those are the part's
        // own modes, not a constrained problem's.
        let hanging = |i: usize| !sys.hang_of.is_empty() && sys.hang_of[i / ndof] != u32::MAX;
        let mut free: Vec<Vec<f64>> = Vec::new();
        for g in &rigid.holds {
            let mut y: Vec<f64> = g.iter().enumerate().map(|(i, v)| if hanging(i) { 0.0 } else { *v }).collect();
            let own = y.iter().zip(&mdiag).map(|(y, m)| y * y * m).sum::<f64>();
            // (twice, for what rounding leaves of the earlier ones)
            for _ in 0..2 {
                for q in &free {
                    let c: f64 = y.iter().zip(q).zip(&mdiag).map(|((y, q), m)| y * q * m).sum();
                    y.iter_mut().zip(q).for_each(|(y, q)| *y -= c * q);
                }
            }
            let left = y.iter().zip(&mdiag).map(|(y, m)| y * y * m).sum::<f64>();
            if left > 1e-10 * own {
                let k = 1.0 / left.sqrt();
                free.push(y.iter().map(|v| v * k).collect());
            }
        }
        rigid_modes = free.len();
        // (a hanging node's dofs are no unknowns: none of the solve's, though
        // the coarsest grid's factor holds them at 1)
        let precond = |r: &[Vec<f64>]| {
            let mut z = mg.vcycle_each(sys, r);
            if !sys.hang_of.is_empty() {
                z.iter_mut().flatten().enumerate().filter(|(i, _)| hanging(*i % sys.ndofs)).for_each(|(_, v)| *v = 0.0);
            }
            z
        };
        // The stiffness, without the terms `rigid` holds those motions at
        // rest by (they are the preconditioner's, which they make positive
        // definite): they would stiffen every mode not orthogonal to them
        // in their own weighting.
        let holding = &sys.lowrank[rigid.support_terms..];
        let stiffness = |v: &[f64]| {
            let mut o = vec![0.0; v.len()];
            sys.apply(v, &mut o);
            for (w, h) in holding {
                let c = w * dot(h, v);
                o.par_iter_mut().zip(h).for_each(|(o, h)| *o -= c * h);
            }
            o
        };
        let mass = |v: &[f64]| v.par_iter().zip(&mdiag).map(|(x, m)| x * m).collect();
        // Each residual under 1e-4 of K x.
        let error = |_: f64, r: &[f64], _: &[f64], kx: &[f64]| dot(r, r).sqrt() / (1e-4 * dot(kx, kx).sqrt().max(1e-300));
        let (vals, vecs, iters_used, modal_converged) = lobpcg(sys.ndofs, &stiffness, &mass, &precond, (s.modes, s.modes), &error, MODAL_ITER, &free);
        modal_ok = modal_converged;
        // Effective mass fractions for the rigid motions: translations x y z
        // and rotations about the centre of mass (x y z): which way a mode
        // moves the structure (a torsion mode shows up as a rotation).
        let mass: f64 = (0..nn).map(|n| mdiag[n * ndof]).sum();
        let com: V3 = std::array::from_fn(|k| (0..nn).map(|n| mdiag[n * ndof] * final_pos[n][k]).sum::<f64>() / mass);
        let rigid: Vec<Vec<f64>> = (0..6)
            .map(|r| {
                let mut v = vec![0.0; sys.ndofs];
                for n in 0..nn {
                    let q: V3 = std::array::from_fn(|k| final_pos[n][k] - com[k]);
                    let w: V3 = if r < 3 {
                        std::array::from_fn(|k| (k == r) as u8 as f64)
                    } else {
                        let a = r - 3;
                        // e_a x q
                        let e: V3 = std::array::from_fn(|k| (k == a) as u8 as f64);
                        [e[1] * q[2] - e[2] * q[1], e[2] * q[0] - e[0] * q[2], e[0] * q[1] - e[1] * q[0]]
                    };
                    for d in 0..3 {
                        v[n * ndof + d] = w[d];
                    }
                }
                v
            })
            .collect();
        let mdot = |a: &[f64], b: &[f64]| -> f64 { a.iter().zip(b).zip(&mdiag).map(|((x, y), m)| x * y * m).sum() };
        let rmr: Vec<f64> = rigid.iter().map(|r| mdot(r, r)).collect();
        for (lam, x) in vals.iter().zip(&vecs) {
            let xmx = mdot(x, x);
            let frac: [f64; 6] = std::array::from_fn(|r| {
                let g = mdot(x, &rigid[r]);
                g * g / xmx / rmr[r].max(1e-300)
            });
            modes_out.push((lam.max(0.0).sqrt() / (2.0 * std::f64::consts::PI), frac));
        }
        // Every mode shape, scaled to 1 mm peak; the first is the displacement field.
        mode_shapes = vecs
            .into_iter()
            .map(|mut x| {
                sys.expand(&mut x, ndof);
                let peak = (0..nn).map(|n| (0..3).map(|d| x[n * 3 + d].powi(2)).sum::<f64>().sqrt()).fold(0.0, f64::max).max(1e-300);
                x.iter().map(|v| v / peak).collect()
            })
            .collect();
        modal_u = mode_shapes.first().cloned();
        modal_iters = iters_used;
    }
    Ok(Modes { modes_out, rigid_modes, mass_kg, modal_u, mode_shapes, modal_iters, modal_ok })
}

// The held values lifted onto the nodes the fixes hold (each the
// penalty-weighted mean of what its clamps hold it at; None where no
// fix holds a dof). A solve starts from them, so its residual is
// measured against the forces of the solution, not the penalty's: a
// fix's right-hand side, beta A v, is about PENALTY times the forces
// it causes.
pub(super) fn lift(model: &Model, sys: &System) -> Vec<Option<f64>> {
    let Study { s, ndof, .. } = *model.st;
    let (clamps, clamp_cells) = (&model.holds.clamps, &model.mesh.clamp_cells);
    {
        let (mut num, mut den) = (vec![0.0; sys.ndofs], vec![0.0; sys.ndofs]);
        for &(cell, k) in clamp_cells.iter().filter(|c| s.fixes[clamps[c.1].fix].film.is_none()) {
            let (cl, nodes) = (&clamps[k], sys.cells[cell].nodes);
            for i in 0..8 {
                let row: f64 = cl.t.mass[i * 8..i * 8 + 8].iter().sum();
                for d in (0..ndof).filter(|&d| s.fixes[cl.fix].dofs[d]) {
                    num[nodes[i] as usize * ndof + d] += cl.t.lift[i * ndof + d];
                    den[nodes[i] as usize * ndof + d] += row;
                }
            }
        }
        num.iter().zip(&den).map(|(n, d)| (*d > 0.0).then(|| n / d)).collect()
    }
}

/// K u = f, from the lifted held values.
pub(super) fn solve_static(sys: &System, mg: &Multigrid, f: &[f64], lift: &[Option<f64>], tol: f64, run: &mut Run) {
    let Run { ref mut x, ref mut iters, ref mut res, .. } = *run;
    // From the lifted held values: a cantilever pushed 1 mm
    // at its tip at 0.5 mm cells, solved to 1e-8 of that, read its
    // wall's reaction 1.7% off its tip's.
    let lift: Vec<f64> = lift.iter().map(|v| v.unwrap_or(0.0)).collect();
    // The correction from it, solved from zero.
    let mut ax = vec![0.0; sys.ndofs];
    sys.apply(&lift, &mut ax);
    let b: Vec<f64> = f.iter().zip(&ax).map(|(f, a)| f - a).collect();
    let precond = |r: &[f64]| mg.vcycle(0, sys, r);
    let (it, r) = pcg(sys, &[], &precond, &b, x, tol, MAX_ITER);
    (*iters, *res) = (*iters + it, r);
    x.iter_mut().zip(&lift).for_each(|(x, l)| *x += l);
}
