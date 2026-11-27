//! Rigid motions the supports leave free (`rigid`): every piece held,
//! a mechanism the loads drive refused, one at rest held there.

use crate::Physics;
use crate::cells::SimCell;
use crate::mechanisms::{driven_failure, free_rigid_motions, rigid_mode};
use crate::numerics::dot;
use crate::system::System;
use ink_core::{V3, v3};
use ink_kernel::mesh::find;
use rayon::prelude::*;

use super::SolveError;
use super::assemble::{Mesh, Model};
use super::holds::Clamp;
use super::study::Study;

/// The pieces' rigid motions: the pieces (each node's, their centroids
/// and sizes; elastic alone), the motions held at rest, those of them that
/// turn, and how many low-rank terms are the supports' own.
#[allow(clippy::type_complexity)]
pub(super) struct Rigid {
    pub(super) pieces: Option<(Vec<usize>, Vec<V3>, Vec<f64>)>,
    pub(super) holds: Vec<Vec<f64>>,
    pub(super) rotations: Vec<usize>,
    pub(super) support_terms: usize,
}

// Rigid motions the supports leave free (a mechanism: a part on one
// pin, or on a frictionless seat alone). Loads that drive one have no
// static solution, which is said instead of solved into a displacement
// set by the solver's tolerance; balanced loads leave the motion at
// rest, held there by a term in it alone (the stiffness does not see a
// rigid motion), and the study says so. Heat needs no hold (a
// temperature fix holds every piece).
// The pieces (kept, so each holds a clamp), their centroids and sizes.
// The pieces of free material, joined through cells with free material
// only (a cell all held keeps 1e-6 of its stiffness, which holds
// nothing), and whether each has a hold of its own: a hold's terms beside
// the springs keeping held material's nodes from floating, in its cells.
// One with none is refused, whatever lost its holds: a bar whose hold's
// edge was lost on a cell plane hung by its held cells on their held
// nodes and read 1,507 mm for 0.19, "converged"; in heat 188,258,757 C.
pub(super) fn rigid(model: &Model, sys: &mut System, f: &[f64], parent: &mut [u32]) -> Result<(Rigid, Vec<String>), SolveError> {
    let clamps = &model.holds.clamps;
    // (Transient heat may leave a piece insulated: its heat capacity holds it.)
    anchored(model, sys)?;
    let pieces = rigid_pieces(model, parent);
    let support_terms = sys.lowrank.len();
    let (holds, warnings, rotations) = mechanisms(model, &pieces, sys, clamps, f, &[])?;
    sys.hold(&holds);
    Ok((Rigid { pieces, holds, rotations, support_terms }, warnings))
}

/// Refused where nothing holds a piece of the part (`rigid`).
fn anchored(model: &Model, sys: &System) -> Result<(), SolveError> {
    let Model { st, holds, lent, mesh, .. } = *model;
    let Study { s, ndof, ne, ref distributing, .. } = *st;
    let clamps = &holds.clamps;
    let Mesh { ref held_geometry, nn, ref hang, ref clamp_cells, ref final_of, ref final_pos, .. } = *mesh;
    if matches!(s.physics, Physics::Elastic { .. }) || s.transient.is_none() {
        let springy = |i: usize| held_geometry[i].as_deref().is_some_and(|g| g.0[0] == 0.0 && (0..ndof).all(|d| g.2[d] && g.1[d][0] == 0.0));
        let mut up: Vec<u32> = (0..nn as u32).collect();
        for c in sys.cells.iter().enumerate().filter(|(i, _)| !springy(*i)).map(|x| x.1) {
            for &n in &c.nodes[1..] {
                let (a, b) = (find(&mut up, c.nodes[0]), find(&mut up, n));
                up[a as usize] = b;
            }
        }
        for (n, masters) in hang.iter() {
            for (m, _) in masters {
                let (a, b) = (find(&mut up, *n), find(&mut up, *m));
                up[a as usize] = b;
            }
        }
        let mut held = std::collections::HashMap::new();
        for (i, c) in sys.cells.iter().enumerate().filter(|(i, _)| !springy(*i)) {
            held.entry(find(&mut up, c.nodes[0])).or_insert((false, i));
        }
        for &(cell, k) in clamp_cells.iter().filter(|cc| !springy(cc.0)) {
            let cl = &clamps[k];
            let r = find(&mut up, sys.cells[cell].nodes[0]);
            let most = cl.t.kh.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            let mut own = cl.t.kh.clone();
            if let Some((beta, m)) = &cl.spring {
                for (ij, v) in m.iter().enumerate() {
                    for d in (0..ndof).filter(|&d| s.fixes[cl.fix].dofs[d]) {
                        own[((ij / 8) * ndof + d) * ne + (ij % 8) * ndof + d] -= beta * v;
                    }
                }
            }
            if own.iter().any(|v| v.abs() > 1e-9 * most) || distributing[cl.fix] && cl.t.mass.iter().any(|v| *v != 0.0) {
                held.entry(r).or_insert((false, cell)).0 = true;
            }
        }
        // (a piece's traction read in the cell beside it holds that cell's
        // material: its own cell may hold none free, its edge on their plane)
        for l in lent.iter() {
            let cell = final_of[l.to] as usize;
            held.entry(find(&mut up, sys.cells[cell].nodes[0])).or_insert((false, cell)).0 = true;
        }
        // (A modal study with no fix is of a free part: its rigid motions
        // are no modes, `mechanisms`.)
        let free = s.modes > 0 && s.fixes.is_empty();
        if let Some(&(_, cell)) = held.values().filter(|h| !h.0 && !free).min_by_key(|h| h.1) {
            let what = if matches!(s.physics, Physics::Heat) { "temperature" } else { "place" };
            return Err(format!(
                "nothing holds the {what} of the part's material about {:?} mm: no fix's edge lies in it (nor a face fix, a mirror or a film on it): give a fix a region whose edge lies in the part's material; where a hold leaves a sliver of the part free beyond it, extend the hold to the face; or refine (res r)",
                final_pos[sys.cells[cell].nodes[0] as usize].map(|v| (v * 1000.0).round() / 1000.0)
            )
            .into());
        }
    }
    Ok(())
}

/// The pieces (kept, so each holds a clamp), their centroids and sizes.
#[allow(clippy::type_complexity)]
fn rigid_pieces(model: &Model, parent: &mut [u32]) -> Option<(Vec<usize>, Vec<V3>, Vec<f64>)> {
    let Study { s, h, .. } = *model.st;
    let Mesh { nn, ref remap, ref final_pos, .. } = *model.mesh;
    matches!(s.physics, Physics::Elastic { .. }).then(|| {
        let mut roots = std::collections::HashMap::new();
        let mut comp_of = vec![0usize; nn];
        for (old, &new) in remap.iter().enumerate().filter(|(_, n)| **n != u32::MAX) {
            let r = find(parent, old as u32);
            let k = roots.len();
            comp_of[new as usize] = *roots.entry(r).or_insert(k);
        }
        let nc = roots.len();
        let mut cent = vec![[0.0; 3]; nc];
        let mut count = vec![0.0; nc];
        for (p, &c) in final_pos.iter().zip(&comp_of) {
            (0..3).for_each(|k| cent[c][k] += p[k]);
            count[c] += 1.0;
        }
        cent.iter_mut().zip(&count).for_each(|(c, n)| *c = c.map(|v| v / n));
        let mut size = vec![h; nc];
        for (p, &c) in final_pos.iter().zip(&comp_of) {
            size[c] = size[c].max(((0..3).map(|k| (p[k] - cent[c][k]).powi(2)).sum::<f64>()).sqrt());
        }
        (comp_of, cent, size)
    })
}

// The motions held at rest and what the study must say of them; one
// the loads drive is a failure.
#[allow(clippy::type_complexity)]
pub(super) fn mechanisms(
    model: &Model,
    pieces: &Option<(Vec<usize>, Vec<V3>, Vec<f64>)>,
    sys: &System,
    clamps: &[Clamp],
    f: &[f64],
    reversed: &[bool],
) -> Result<(Vec<Vec<f64>>, Vec<String>, Vec<usize>), String> {
    let Study { s, ref distributing, .. } = *model.st;
    let Mesh { ref final_pos, ref clamp_cells, .. } = *model.mesh;
    let Some((comp_of, cent, size)) = &pieces else { return Ok(Default::default()) };
    let nc = cent.len();
    // The supports' energy over each piece's rigid motions.
    // As with the other element assemblies, compute each cell in
    // parallel, then scatter in input order. A worker-local fold and
    // Rayon reduce regroup these sums with scheduling; the resulting
    // null-space vectors set the mechanism's stabilization operator.
    let parts: Vec<(usize, [f64; 36], f64)> = clamp_cells
        .par_iter()
        .map(|&(cell, k)| {
            let cl = &clamps[k];
            let nodes = sys.cells[cell].nodes;
            let c = comp_of[nodes[0] as usize];
            // (a distributing hold's is its low-rank terms')
            let mut a = if distributing[cl.fix] { vec![0.0; 576] } else { cl.t.kh.clone() };
            // (A symmetry plane's, which holds by its penalty alone: its
            // mass, on the axes it leaves free.)
            if reversed.get(cl.fix) == Some(&true) {
                for i in 0..8 {
                    for j in 0..8 {
                        for d in 0..3 {
                            a[(i * 3 + d) * 24 + j * 3 + d] = if s.fixes[cl.fix].dofs[d] { 0.0 } else { cl.t.mass[i * 8 + j] };
                        }
                    }
                }
            }
            // (Not the spring keeping held material's nodes from floating:
            // a body only it holds is free, and refused where the loads
            // drive it. Counted, a hold whose edge was lost read a bar
            // pulled 1 MPa 8,007 mm long for 0.1995, "converged".)
            let mut mass: f64 = cl.t.mass.iter().sum();
            if let Some((beta, m)) = cl.spring.as_ref().filter(|_| !distributing[cl.fix]) {
                for (ij, v) in m.iter().enumerate() {
                    for d in (0..3).filter(|&d| s.fixes[cl.fix].dofs[d]) {
                        a[((ij / 8) * 3 + d) * 24 + (ij % 8) * 3 + d] -= beta * v;
                    }
                    mass -= beta * v;
                }
            }
            // R (24 x 6: the rigid motions at the cell's nodes), then R^T A R.
            let r: Vec<f64> = (0..144).map(|ip| rigid_mode(final_pos[nodes[ip / 18] as usize], cent[c], size[c], ip % 6)[(ip / 6) % 3]).collect();
            let ar: Vec<f64> = (0..144).map(|ip| (0..24).map(|j| a[(ip / 6) * 24 + j] * r[j * 6 + ip % 6]).sum()).collect();
            let block = std::array::from_fn(|pq| (0..24).map(|i| r[i * 6 + pq / 6] * ar[i * 6 + pq % 6]).sum::<f64>());
            (c, block, mass.max(0.0))
        })
        .collect();
    let (mut blocks, mut scale) = (vec![[0.0; 36]; nc], vec![0.0; nc]);
    for (c, block, s) in parts {
        blocks[c].iter_mut().zip(block).for_each(|(x, y)| *x += y);
        scale[c] += s;
    }
    // (a piece nothing but springs holds, a body held through: each of
    // its motions free of the scale, which would otherwise divide by 0)
    scale.iter_mut().for_each(|s| {
        if !(*s > 0.0) {
            *s = 1.0;
        }
    });
    let (mut holds, mut warnings) = (Vec::new(), Vec::new());
    let (mut rotations, mut free) = (Vec::new(), Vec::new());
    for (what, g, turns) in free_rigid_motions(final_pos, comp_of, cent, size, &blocks, &scale, &sys.lowrank) {
        // The loads' work on it per unit: the force along a slide (N),
        // the moment about a turn's axis (N mm). Balanced loads keep a
        // remainder from their sampling on the grid (a pull through the
        // axis, the cells not quite symmetric about it: 5e-5 of the
        // loads' work on it node by node), so under a thousandth of that
        // (each node's work along it, all added as if one way) they are
        // taken as balanced, and the warning says how far: a load whose
        // line misses the motion by under about 0.057 degrees, as a pin
        // with friction 0.001 would hold it. Loads square to the motion
        // do no work on it and count for nothing (a balanced 1000 N
        // pinch through a plate, counted by its size, let a 30 N mm
        // twist about its pin pass as balanced); and rounding, a trace
        // of the motion along a load square to it (a ball in a socket
        // pushed straight down, its spin "driven" at 0 N mm), is under
        // 1e-9 of each node's load times its motion.
        if turns {
            rotations.push(holds.len());
        }
        let work = dot(f, &g);
        let along: f64 = f.chunks(3).zip(g.chunks(3)).map(|(a, b)| (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]).abs()).sum();
        let most: f64 = f.chunks(3).zip(g.chunks(3)).map(|(a, b)| v3::norm([a[0], a[1], a[2]]) * v3::norm([b[0], b[1], b[2]])).sum();
        let spread = along.max(1e-6 * most);
        let drives = work.abs() > 1e-3 * spread;
        let unit = if what.starts_with("turn") { "N mm about its axis" } else { "N along it" };
        // A modal study's free motion is a mode at 0 Hz, no vibration: it is
        // projected out of the modes (`modal`), which are the elastic ones.
        if s.modes > 0 {
            free.push(what);
            holds.push(g);
            continue;
        }
        if drives {
            return Err(driven_failure(&what, work, unit));
        }
        warnings.push(format!(
            "the supports leave the part free to {what}: a mechanism; the loads do not drive it ({:.3e} {unit}, \
                 {:.1e} of their work on it node by node), so it is held at rest: the displacement shown has none of it \
                 (its share, each node weighted by its stiffness, is zero)",
            work,
            work.abs() / spread.max(1e-300)
        ));
        holds.push(g);
    }
    if !free.is_empty() {
        let n = free.len();
        let (s_, it) = if n == 1 { ("", "it") } else { ("s", "them") };
        warnings.push(if s.fixes.is_empty() {
            format!("nothing holds the part, free to {}: {n} rigid motion{s_}, at 0 Hz, left out of the modes, which are its elastic modes alone", free.join("; "))
        } else {
            format!(
                "the supports leave the part free to {}: {n} rigid motion{s_}, at 0 Hz, left out of the modes, which are its elastic modes alone; a (min-frequency ...) claim counts {it} unless it says (excluding rigid-motions)",
                free.join("; ")
            )
        });
    }
    Ok((holds, warnings, rotations))
}

/// Each node's reading reliable: inside material, or a full cell's corner.
pub(super) fn node_reliable(model: &Model, sys: &System) -> Vec<bool> {
    let Study { s, ctx, .. } = *model.st;
    let Mesh { nn, ref final_pos, .. } = *model.mesh;
    let mut support = vec![0.0f64; nn];
    for c in &sys.cells {
        for &n in &c.nodes {
            support[n as usize] += c.frac;
        }
    }
    // Inside test per node: corners of full cells are material; others use
    // an adjacent cut cell's simplified tape.
    let mut node_cell: Vec<Option<u32>> = vec![None; nn];
    let mut on_full = vec![false; nn];
    // (A cell whole but for rounding is whole, as a sample is: on 0.1 mm
    // cells a bar's end face, on a cell plane, had its cells 0.99999999999999
    // full and its nodes at x = 6.000000000000001, where the shape reads
    // +8.9e-16, outside: none of them read, its peak displacement was taken
    // a cell short, 0.0580 for 0.059.)
    for (ci, c) in sys.cells.iter().enumerate() {
        for &n in &c.nodes {
            if c.frac > 1.0 - 1e-9 {
                on_full[n as usize] = true;
            } else {
                node_cell[n as usize].get_or_insert(ci as u32);
            }
        }
    }
    let node_reliable: Vec<bool> = (0..nn)
        .into_par_iter()
        .map(|n| {
            // (a full cell's corner is its own field's, however few cells
            // meet there: a cube's far corner, one cell's, was left out and
            // its peak read a cell short along an edge, 0.0765 for 0.0771)
            if on_full[n] {
                return true;
            }
            if support[n] < 2.0 {
                return false;
            }
            let tape = node_cell[n].and_then(|ci| sys.cells[ci as usize].tape.as_deref()).unwrap_or(&s.shape);
            crate::eval::eval_point(tape, ctx, final_pos[n]) <= 0.0
        })
        .collect();
    node_reliable
}

impl Model<'_> {
    /// The held material's share of each cell's nodes on the axes it holds
    /// (`Clamp::hn`): its mass moves with the hold and its heat capacity is the held
    /// temperature's, so neither counts in a mode or a time step (they did,
    /// and held material's nodes, held by ALPHA's spring alone, would ring).
    pub(super) fn held_share(&self, cells: &[SimCell]) -> std::collections::HashMap<usize, [[f64; 3]; 8]> {
        let model = self;
        let ndof = self.st.ndof;
        let (clamps, clamp_cells) = (&self.holds.clamps, &self.mesh.clamp_cells);
        let mut out = std::collections::HashMap::new();
        for &(cell, k) in clamp_cells.iter() {
            let cl = &clamps[k];
            if cl.hn.iter().flatten().all(|v| *v == 0.0) {
                continue;
            }
            let (vol, _) = model.sample_cell(&cells[cell], None, None);
            let e = out.entry(cell).or_insert([[0.0; 3]; 8]);
            for i in 0..8 {
                for d in 0..ndof {
                    e[i][d] = (e[i][d] + cl.hn[d][i] / vol[i].max(1e-300)).min(1.0);
                }
            }
        }
        out
    }
}
