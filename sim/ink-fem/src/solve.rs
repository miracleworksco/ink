//! A solve, from a `Setup` to a `Solution` (`solve_with`): cells, holds,
//! loads, the linear solve and what is read off it.
//!
//! `solve_with` reads as the solve's table of contents; each phase is its
//! module, in order: `study` (the grid and material fields), `classify` and
//! `settle` (the cells), `holds` and `held` (what the fixes hold), `connect`
//! (the pieces), `assemble` (the system), `freeing` (supports' free motions),
//! `loads`, `rigid` (mechanisms), `run`, `redesign` and `over_time` (the
//! solves), and `readout` (what is read off the solution).

use crate::cells::place;
use crate::element::{M3, m3_det};
use crate::eval::Ctx;
use crate::memory;
use crate::multigrid::Multigrid;
use crate::numerics::dot;
use crate::report::{Incomplete, Report, cells_of, glued_warning, material_left_warning, spread_of, thin_left_warning};
use crate::solution::{NO_FIX, Solution};
use crate::surface::Facet;
use crate::{Fix, Held, Physics, RegionSel, Setup};
use ink_core::{V3, v3};
use rayon::prelude::*;
use std::time::Instant;

mod assemble;
mod classify;
mod connect;
mod freeing;
mod held;
mod holds;
mod loads;
mod over_time;
mod readout;
mod redesign;
mod rigid;
mod run;
mod settle;
mod study;

use assemble::{Model, assemble};
use classify::{classify, refine};
use connect::connect;
use freeing::free_supports;
use held::held_cells;
use holds::{Surfaces, hold_fixes};
use loads::{Loads, loads};
use over_time::transient;
use readout::{Buckled, Reactions, Supports, buckling, field_of, fix_footprints, fully_held, held_reads, reactions, supports, surface_reads};
use redesign::run_redesign;
use rigid::{node_reliable, rigid};
use run::{Modes, Rhs, Run, lift, modal, solve_static};
use settle::{Unresolved, settle};
use study::Study;

/// Clamp penalty, relative to the material coefficient / h^2.
const PENALTY: f64 = 1e3;

/// A torque `t` (N mm) on weighted nodes at `pos`: forces w_n (omega x
/// r_n), r_n from their weighted centroid, with J omega = t (J their
/// inertia tensor about it), which sum to no force and to exactly the
/// torque. Regularised: a region along a line has no inertia about it, and
/// takes no torque about it.
fn torque_forces(nodes: &[(u32, f64)], pos: impl Fn(u32) -> V3, t: V3) -> Vec<V3> {
    let (c, j) = torque_frame(nodes, &pos);
    let omega = solve3(&j, t);
    nodes.iter().map(|&(n, w)| v3::scale(v3::cross(omega, v3::sub(pos(n), c)), w)).collect()
}

/// The weighted centroid of a torque's nodes, and their inertia tensor
/// about it (regularised).
fn torque_frame(nodes: &[(u32, f64)], pos: &impl Fn(u32) -> V3) -> (V3, M3) {
    let wsum: f64 = nodes.iter().map(|n| n.1).sum();
    let c: V3 = std::array::from_fn(|k| nodes.iter().map(|&(n, w)| w * pos(n)[k]).sum::<f64>() / wsum);
    let mut j = [[0.0f64; 3]; 3];
    for &(n, w) in nodes {
        let r = v3::sub(pos(n), c);
        let r2 = v3::dot(r, r);
        for a in 0..3 {
            for b in 0..3 {
                j[a][b] += w * (if a == b { r2 } else { 0.0 } - r[a] * r[b]);
            }
        }
    }
    let tr = (j[0][0] + j[1][1] + j[2][2]).max(1e-300);
    for (a, row) in j.iter_mut().enumerate() {
        row[a] += 1e-12 * tr;
    }
    (c, j)
}

/// J x = b by Cramer's rule.
fn solve3(j: &M3, b: V3) -> V3 {
    let dj = m3_det(j);
    std::array::from_fn(|k| {
        let mut m = *j;
        for (a, row) in m.iter_mut().enumerate() {
            row[k] = b[a];
        }
        m3_det(&m) / dj
    })
}

#[derive(Debug)]
pub enum SolveError {
    InvalidMirror(String),
    Failed(String),
}

impl From<String> for SolveError {
    fn from(message: String) -> Self {
        Self::Failed(message)
    }
}

impl From<&str> for SolveError {
    fn from(message: &str) -> Self {
        Self::Failed(message.into())
    }
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidMirror(s) | Self::Failed(s) => f.write_str(s),
        }
    }
}

pub fn solve(s: &Setup, ctx: &Ctx) -> Result<Solution, String> {
    solve_checked(s, ctx).map_err(|e| e.to_string())
}

pub fn solve_checked(s: &Setup, ctx: &Ctx) -> Result<Solution, SolveError> {
    sized(s, || solve_with(s, ctx, None))
}

/// `f`, a solve, on as many threads as its grid keeps busy, if that is
/// fewer than the caller's. A solve is thousands of short parallel loops:
/// past a point a thread added only waits, and on a busy machine it waits
/// on threads the others hold. The square root of the grid box's cells,
/// over 16; a 40 x 10 x 10 mm bar's solve on 32 threads (load 30-45), wall
/// / CPU s, before: 720 cells (now 2 threads) 0.45 / 0.61 -> 0.37 / 0.32,
/// 3,500 (4) 0.73 / 1.68 -> 0.58 / 1.01, 21,700 (10) 2.6 / 9.4 -> 2.5 /
/// 7.7, 82,700 (18) 3.2 / 34 -> 2.8 / 28.5; examples/layout.ink's topology
/// (2,880 cells, 4) 24 -> 14 CPU s, its part check among them. Its threads
/// have a command's stack (`f` runs on one); the answer is the same on any
/// number (`ordered_sum`).
fn sized<R: Send>(s: &Setup, f: impl FnOnce() -> R + Send) -> R {
    let cells: f64 = (0..3).map(|i| ((s.hi[i] - s.lo[i]) / s.h).ceil().max(1.0)).product();
    let threads = (cells.sqrt() / 16.0).ceil().max(1.0) as usize;
    if threads >= rayon::current_num_threads() {
        return f();
    }
    match rayon::ThreadPoolBuilder::new().num_threads(threads).stack_size(256 << 20).build() {
        Ok(pool) => pool.install(f),
        Err(_) => f(),
    }
}

/// The assembled static elastic study a redesign drives
/// (`solve_redesigned`): its cells (a cut cell's pieces each one), whose
/// stiffness it scales one by one and solves at, as often as it likes.
pub trait Assembly {
    /// How many cells.
    fn cells(&self) -> usize;
    /// The cell size (mm).
    fn h(&self) -> f64;
    /// Cell i's grid index, centre and material volume (mm^3).
    fn cell(&self, i: usize) -> ([i64; 3], V3, f64);
    /// The cells at a grid index: none, one, or a cut cell's pieces.
    fn cells_at(&self, idx: [i64; 3]) -> std::ops::Range<usize>;
    /// Each cell: whether a fix holds all its material (that hold is then
    /// in the matrices of the cells beside it, where its free material is).
    fn held(&self) -> Vec<bool>;
    /// Each cell's stiffness times its multiplier (supports' penalties
    /// with it), solved from the last solve's displacement to relative
    /// residual `tol`: the compliance (N mm, the loads' work).
    fn solve(&mut self, mult: &[f64], tol: f64) -> Result<f64, String>;
    /// The derivative of u^T K u with respect to each cell's multiplier,
    /// at the last solve's displacement u, including the support terms.
    /// With stationary supports their weighted sum is the compliance and
    /// their negatives are its slopes.
    fn energies(&self) -> Vec<f64>;
}

/// A static elastic solve whose cells' stiffness `redesign` sets (a package
/// built on the kernel may): it is handed the
/// assembled study in place of the one solve, and the solution is that of its
/// last solve (with its multipliers: the stresses are the scaled material's).
pub fn solve_redesigned(s: &Setup, ctx: &Ctx, redesign: &mut (dyn FnMut(&mut dyn Assembly) -> Result<(), String> + Send)) -> Result<Solution, String> {
    sized(s, || solve_with(s, ctx, Some(redesign))).map_err(|e| e.to_string())
}

#[cfg(test)]
static LAST_SOLVE_THREADS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

fn solve_with(s: &Setup, ctx: &Ctx, redesign: Option<&mut dyn FnMut(&mut dyn Assembly) -> Result<(), String>>) -> Result<Solution, SolveError> {
    #[cfg(test)]
    LAST_SOLVE_THREADS.store(rayon::current_num_threads(), std::sync::atomic::Ordering::Relaxed);
    let t0 = Instant::now();
    let redesigned = redesign.is_some();
    // The grid over the part, a cell's samples, and the material fields.
    let st = Study::new(s, ctx)?;
    let (ndof, h) = (st.ndof, st.h);

    // 1-2. Cells: full blocks from the octree, cut cells sampled.
    let (classified, said, missed_at) = classify(&st);
    // 3. Split where they hold thin features.
    let (split, short) = refine(&st, classified)?;
    // Slivers dropped or taken by the cell beside them, laminates, and nodes.
    let cells = settle(&st, split, short, said)?;

    // 4. Clamps: what each fix holds (`hold_fixes`), and each held cell's
    // free stiffness and its holds' terms (`held_cells`).
    let mut surf = Surfaces::new(&st, &cells);
    let mut holds = hold_fixes(&st, &cells, &mut surf)?;
    let mut held = held_cells(&st, &cells, &mut holds)?;
    // 5. The pieces the cells join into, each anchored, or dropped.
    let mut joined = connect(&st, &cells, &holds, &held, missed_at);
    // 6. Final numbering and element matrices.
    let (mut sys, mut f, mesh) = assemble(&st, &cells, &holds, &mut held, &joined.keep);
    let model = Model { st: &st, cells: &cells, holds: &holds, lent: &held.lent, mesh: &mesh };
    // Supports with free rigid motions (`freeing`).
    let mut freed = free_supports(&model, &mut sys, &mut f)?;

    // What the fixes hold the system at (their penalty's right-hand side),
    // apart from the loads: the compliance is the loads' work alone. With it
    // a part moved 1 mm as a whole by a fix, free of stress, read 1.9e8 N mm.
    let f_held = f.clone();
    // 7. Load vector: body force and region loads, sampled per cell.
    let mut loads = loads(&model, &mut sys, &mut f, &mut surf, &mut joined.incomplete, (joined.dropped_cells, joined.dropped_volume))?;
    // The held material's samples served the loads alone (`held_still`):
    // as many as its cells' samples, a large held region's largest part.
    holds.tables.values_mut().for_each(|t| t.held = Vec::new());
    let model = Model { st: &st, cells: &cells, holds: &holds, lent: &held.lent, mesh: &mesh };
    // What lands on hanging nodes, to the nodes they hang on.
    sys.fold(&mut f, ndof);
    let mut f_held = f_held;
    sys.fold(&mut f_held, ndof);
    let f_loads: Vec<f64> = f.iter().zip(&f_held).map(|(a, b)| a - b).collect();
    let mut rhs = Rhs { f, f_held, f_loads };
    // Rigid motions the supports leave free: refused where the loads drive
    // one, held at rest where they do not.
    let (rigid, mut warnings) = rigid(&model, &mut sys, &rhs.f, &mut joined.parent)?;
    warnings.extend(holds.notes.iter().cloned());
    let setup_ms = t0.elapsed().as_secs_f64() * 1e3;
    let node_reliable = node_reliable(&model, &sys);

    // 8. Solve K u = f (clamps make K positive definite) by preconditioned CG,
    // or, for modal analysis, K x = w^2 M x for the lowest modes.
    let t1 = Instant::now();
    let mg = Multigrid::build(&sys, &mesh.final_idx);
    let Modes { modes_out, rigid_modes, mass_kg, modal_u, mut mode_shapes, modal_iters, modal_ok } = modal(&model, &sys, &mg, &rigid)?;
    let lift = lift(&model, &sys);
    let mut run = Run::new(sys.ndofs);
    if s.modes == 0 && s.transient.is_none() && redesign.is_none() {
        solve_static(&sys, &mg, &rhs.f, &lift, s.tol, &mut run);
    }
    if let Some(redesign) = redesign {
        run_redesign(&model, &mut sys, &mut run, &rhs, &rigid.holds, &mut freed, redesign)?;
    }
    let reads = surface_reads(&model, &sys, &mut surf);
    if let (Some(tr), 0) = (s.transient, s.modes) {
        transient(&model, tr, &mut sys, &mut rhs, &mut loads, &mut joined, &mut warnings, &node_reliable, &lift, &reads, &surf, &mut run)?;
    }
    let levels = mg.levels.len() + 1;
    let mut u = modal_u.unwrap_or(std::mem::take(&mut run.x));
    sys.expand(&mut u, ndof);
    let solve_ms = t1.elapsed().as_secs_f64() * 1e3;

    // 9. Reactions: what each hold pushes its held value with.
    let Reactions { reactions, reaction_moments, free_motions, amplitudes } = reactions(&model, &sys, &u, &freed, &loads);
    let energy = dot(&u, &rhs.f_loads);
    let clamped_volume: f64 = mesh.clamp_cells.iter().map(|&(_, k)| holds.clamps[k].vol).sum();

    // 10. Cell-centre stress (or heat flux).
    let mut cell_field: Vec<[f64; 6]> = sys.cells.par_iter().map(|c| field_of(&st, &sys, c, &u)).collect();

    // What the study must say: the supports a stress peak keeps clear of, a
    // load in a fix or by it, what the cells leave unresolved.
    joined.incomplete.extend(joined.duplicate_cut.iter().enumerate().filter(|(_, duplicate)| **duplicate).map(|(fix, _)| Incomplete::DuplicateCut { fix }));
    let Supports { near_support, possible_cuts, unstated_holds, cuts, planes } = supports(&model, &sys, &mut surf);
    let fix_footprints = fix_footprints(&model, &sys, &mut surf);
    let held_in = fully_held(&model, &sys);
    let near: Vec<bool> = near_support.iter().map(|&f| f != NO_FIX).collect();
    warnings.extend(std::mem::take(&mut joined.insulated_notes));
    warnings.extend(std::mem::take(&mut joined.idle));
    warnings.extend(load_warnings(s, &loads.landed, &held_in, &near));
    let load_spread: Vec<(String, V3)> =
        s.loads.iter().zip(&loads.landed).map(|(l, cells)| (l.name.clone(), spread_of(cells.iter().map(|&(c, w)| (place(&st.grid, &sys.cells[c]).0, w))))).collect();
    let Unresolved { stepped, glued, thin_left, material_left, interface_shared } = cells.unresolved;
    let depth = cells.depth;
    warnings.extend(glued.map(|(n, at, w)| glued_warning(n, at, w, &cells_of(h, depth))));
    warnings.extend(thin_left.map(|(n, at)| thin_left_warning(n, at, &cells_of(h, depth))));
    warnings.extend(material_left.map(|(n, at)| material_left_warning(n, at, &cells_of(h, depth))));
    warnings.extend(surf.measures.iter().filter_map(|m| m.unseen()).map(|u| format!("{u}: its claims are unknown. Refine (res r), or give the region faces the cells hold")));
    warnings.extend(loads.pressure_lost.iter().map(|(name, lost, all)| {
        format!(
            "load '{name}' (a pressure): {:.3} N of its {:.3} N, over facets of the part's surface no cell holds, was not applied: the study \
             solved the part under less load than written. Refine (res r), or press faces the cells hold",
            lost, all
        )
    }));
    let Run { iters, res, history, peak_time, peak, time_levels, .. } = run;
    let incomplete = joined.incomplete;
    // Linear buckling of a static elastic solve.
    let Buckled { buckle_out, buckle_reversed, modeled_buckle, parity_sectors } =
        buckling(&model, &mut sys, mg, &u, res, redesigned, &cell_field, &rigid, &freed, &cuts, &incomplete, &surf.measures, &mut mode_shapes)?;
    if let Some(b) = &buckle_out {
        match b.answer() {
            Some((low, f)) if low <= 1.0 => warnings.push(buckled_warning(low, f, b.at, h)),
            Some(_) => {}
            None => warnings
                .push(format!("buckling not resolved: {}", b.failure.as_deref().unwrap_or("the eigensolver gave no factor. Refine (res r) or check the supported whole part"))),
        }
    }
    let converged = if s.modes > 0 {
        modal_ok
    } else if s.transient.is_some() {
        res <= s.tol.max(1e-6)
    } else {
        res <= s.tol
    };
    let finite = |v: &[f64]| v.iter().all(|x| x.is_finite());
    let modal = s.modes > 0;
    let failure = if !(finite(&u) && energy.is_finite() && (modal || res.is_finite()) && modes_out.iter().all(|m| m.0.is_finite())) {
        Some(format!(
            "the solve broke down: its residual or solution is not finite (relative residual {res:e} after {} iterations): a material \
             or load coefficient not finite, or a singular stiffness",
            if modal { modal_iters } else { iters }
        ))
    } else if !converged && modal {
        Some(format!("the eigensolver did not converge to the {} lowest modes in {modal_iters} iterations", s.modes))
    } else if !converged {
        Some(format!(
            "the solver did not reach its tolerance {:e}: relative residual {res:e} after {iters} iterations",
            if s.transient.is_some() { s.tol.max(1e-6) } else { s.tol }
        ))
    } else {
        None
    };
    let Loads { applied, applied_moments, pressure_error, footprints: load_footprints, .. } = loads;
    let report = Report {
        memory_estimate: memory::needed(s, cells.counted),
        cells_full: mesh.n_full,
        cells_cut: mesh.n_cut,
        cells_refined: sys.cells.iter().filter(|c| c.level > 0).count(),
        refine_depth: depth,
        nodes: mesh.nn,
        dofs: sys.ndofs,
        clamped_volume,
        held_edge: held.held_edge,
        modes: modes_out,
        rigid_modes,
        mass_kg,
        dropped_cells: joined.dropped_cells,
        dropped_volume: joined.dropped_volume,
        dropped_at: joined.dropped_at,
        incomplete,
        pressure_error,
        sliver_volume: cells.slivers.volume,
        volume: mesh.volume,
        // (with the slivers each took: `gained`, as the part's volume has them)
        unresolved_volume: sys
            .cells
            .iter()
            .enumerate()
            .filter(|(_, c)| c.frac < 0.5)
            .map(|(i, c)| (c.frac + cells.slivers.gained.get(&mesh.classified_of[i]).map_or(0.0, |g| g.1)) * place(&st.grid, c).1.powi(3))
            .sum(),
        iterations: if s.modes > 0 { modal_iters } else { iters },
        levels,
        residual: if s.modes > 0 { 0.0 } else { res },
        converged,
        failure,
        history,
        peak_time,
        peak,
        time_levels,
        cached: false,
        setup_ms,
        solve_ms,
        applied,
        reactions,
        applied_moments: if ndof == 3 { applied_moments } else { Vec::new() },
        reaction_moments: if ndof == 3 { reaction_moments } else { Vec::new() },
        energy,
        free_motions,
        warnings,
        buckle: buckle_out,
        buckle_reversed,
        modeled_buckle,
        parity_sectors,
        possible_cuts,
        unstated_holds,
        glued,
        thin_left,
        material_left,
        interface_shared,
        interface_stepped: stepped,
        surface_measures: surf.measures,
        load_spread,
        load_footprints,
        fix_footprints,
        centroid: {
            let (mut w, mut m) = (0.0, [0.0; 3]);
            for c in &sys.cells {
                let (ctr, size) = place(&st.grid, c);
                let v = c.frac * size.powi(3);
                w += v;
                (0..3).for_each(|a| m[a] += v * ctr[a]);
            }
            m.map(|v| v / w.max(1e-300))
        },
    };
    // What is read at the nodes inside held material.
    let held = held_reads(&model, &sys, &u, &amplitudes, &mut cell_field);
    let mut sol = Solution {
        physics: s.physics,
        h,
        bounds: (s.lo, s.hi),
        origin: st.grid.origin,
        dims: st.dims,
        cells: sys.cells,
        cell_at: mesh.cell_at,
        node_idx: mesh.final_idx,
        depth,
        node_reliable,
        u,
        held,
        cell_field,
        near_support,
        plane_holds: planes,
        surface_points: reads.points,
        surface_own: reads.own,
        surface_normals: reads.normals,
        surface_held: reads.held,
        climbed: Vec::new(),
        peaks: std::sync::OnceLock::new(),
        reads: std::sync::OnceLock::new(),
        report,
        mode_shapes,
    };
    // (A peak stress's claims read it, `climbed_scalar`; heat's none.)
    if s.modes == 0 && matches!(s.physics, Physics::Elastic { .. }) && sol.valid().is_ok() {
        sol.climb_reads(&s.shape, ctx);
    }
    Ok(sol)
}

/// A uniformly held thin layer has no interior support edge when it
/// covers every sampled point of each facet facing the held direction.
/// Include the nearly facing facets that round an edge in the mesh.
fn thin_normal_restraint(fx: &Fix, held: f64, facets: &[Facet], samples: &[V3], ctx: &Ctx, h: f64) -> bool {
    let RegionSel::Field(t) = &fx.region else { return false };
    if fx.film.is_some() || !fx.free.is_empty() || fx.value.iter().any(|v| matches!(v, Held::Field(_))) {
        return false;
    }
    let normals = if fx.normal {
        let at: Vec<_> = facets.iter().map(|f| f.at).collect();
        Some(crate::eval::eval_grads(t, ctx, &at).into_iter().map(|g| v3::direction([g[1], g[2], g[3]]).unwrap_or([0.0; 3])).collect::<Vec<_>>())
    } else {
        if fx.dofs.iter().filter(|d| **d).count() != 1 {
            return false;
        }
        None
    };
    let axis = fx.dofs.map(|d| if d { 1.0 } else { 0.0 });
    let selected = fx.region.holds_on_surface(ctx, samples, h);
    let mut area = 0.0;
    for (i, (f, inside)) in facets.iter().zip(selected.as_chunks::<4>().0).enumerate() {
        if v3::dot(f.normal, normals.as_ref().map_or(axis, |v| v[i])).abs() > 0.9 {
            match inside.iter().filter(|d| **d).count() {
                4 => area += f.area,
                0 => {}
                _ => return false,
            }
        }
    }
    area > 0.0 && held <= 0.1 * h * area
}

/// A uniform displacement normal to a complete flat exterior face.
fn flat_normal_hold(fx: &Fix, facets: &[Facet], samples: &[V3], ctx: &Ctx, h: f64) -> Option<(usize, f64)> {
    if fx.normal || fx.film.is_some() || !fx.free.is_empty() || fx.dofs.iter().filter(|d| **d).count() != 1 || fx.value.iter().any(|v| matches!(v, Held::Field(_))) {
        return None;
    }
    let axis = fx.dofs.iter().position(|d| *d)?;
    let selected = fx.region.holds_on_surface(ctx, samples, h);
    let mut plane = None;
    for (f, inside) in facets.iter().zip(selected.as_chunks::<4>().0) {
        if f.normal[axis].abs() < 1.0 - 1e-8 || !inside.iter().any(|v| *v) {
            continue;
        }
        if !inside.iter().all(|v| *v) {
            return None;
        }
        let coord = f.corners[0][axis];
        let at = plane.get_or_insert(coord);
        if (coord - *at).abs() > 1e-8 * h {
            return None;
        }
    }
    let at = plane?;
    let mut sides = [false; 2];
    for (f, inside) in facets.iter().zip(selected.as_chunks::<4>().0) {
        for p in f.corners {
            if p[axis] < at - 1e-8 * h {
                sides[0] = true;
            }
            if p[axis] > at + 1e-8 * h {
                sides[1] = true;
            }
        }
        if f.normal[axis].abs() > 1.0 - 1e-8 && (f.corners[0][axis] - at).abs() < 1e-8 * h && !inside.iter().all(|v| *v) {
            return None;
        }
    }
    (sides[0] != sides[1]).then_some((axis, at))
}

/// What a study must say of loads that land where they load nothing: more
/// than half of one inside a region a fix holds entirely goes straight
/// into its reaction; more than half within the support margin causes
/// stress the headline peak leaves out. `landed`: per load, its cells
/// and their shares; `held`: per cell, the share of its material held.
fn load_warnings(s: &Setup, landed: &[Vec<(usize, f64)>], held: &[f64], near_support: &[bool]) -> Vec<String> {
    let mut out = Vec::new();
    for (l, cells) in s.loads.iter().zip(landed) {
        let total: f64 = cells.iter().map(|c| c.1).sum();
        if total <= 0.0 {
            continue;
        }
        let share = |f: &dyn Fn(usize) -> f64| cells.iter().map(|&(c, w)| w * f(c)).sum::<f64>() / total;
        let inside = share(&|c| held[c].min(1.0));
        let near = share(&|c| if near_support[c] { 1.0 - held[c].min(1.0) } else { 0.0 });
        let what = if s.physics.ndof() == 3 { "load" } else { "source" };
        if inside > 0.5 {
            out.push(format!(
                "{what} '{}': {:.0}% of it lies inside a fix's region, which holds that material still: that share goes \
                 straight into the fix's reaction and loads nothing else. Move the {what}'s region off the fix if it is to load the part",
                l.name,
                inside * 100.0
            ));
        } else if near > 0.5 {
            out.push(format!(
                "load '{}': {:.0}% of it lies within the support margin ({} mm of a fix): the stress it causes there is \
                 left out of the peak stress and the claims on it (max_von_mises_near_supports_mpa reads it); \
                 set a smaller (support-margin mm) to count it, or model the support's contact",
                l.name,
                near * 100.0,
                s.support_margin.max(2.0 * s.h)
            ));
        }
    }
    out
}

/// What a solve says of its linear buckling factor (its band, `low` to
/// `factor`) where the loads may reach it (`low` at most 1).
pub fn buckled_warning(low: f64, factor: f64, at: V3, h: f64) -> String {
    let at = at.map(|v| (v * 100.0).round() / 100.0);
    let n = |v: f64| format!("{}", (v * 1e4).round() / 1e4);
    format!(
        "buckles: the loads may reach its linear buckling load (its factor {} to {}, by linear buckling analysis on cells of {h} mm; the \
         first mode moves most at {at:?}). A linear static study does not model the buckle: its stresses, displacements and safety \
         factors say nothing past it. Lower the load or stiffen the part where its buckling mode, (study-field s displacement (mode 1)), moves",
        n(low),
        n(factor)
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Boundary, Load};
    use ink_core::graph::{BinOp, Graph, UnOp, Var};
    use ink_core::tape::Tape;
    use std::sync::Arc;

    /// `solve_checked` caps its worker pool to the grid size: a small 4 x 2 x 2
    /// box (16 cells) asked inside a 4-thread pool runs on 1 thread.
    #[test]
    fn solve_checked_caps_threads_on_a_small_grid() {
        let mut g = Graph::new();
        let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
        let (ax, ay, az) = (g.un(UnOp::Abs, x), g.un(UnOp::Abs, y), g.un(UnOp::Abs, z));
        let (c2, c1) = (g.constant(2.0), g.constant(1.0));
        let (dx, dy, dz) = (g.bin(BinOp::Sub, ax, c2), g.bin(BinOp::Sub, ay, c1), g.bin(BinOp::Sub, az, c1));
        let dxy = g.bin(BinOp::Max, dx, dy);
        let root = g.bin(BinOp::Max, dxy, dz);
        let shape = Arc::new(Tape::new(&g, root));
        let one = Arc::new(Tape::new(&g, c1));
        let left = g.bin(BinOp::Add, x, c1);
        let right = g.bin(BinOp::Sub, c1, x);
        let setup = Setup {
            physics: Physics::Heat,
            shape,
            coef: one,
            coef_const: Some(1.0),
            nu: None,
            strength: None,
            body: None,
            fixes: vec![Fix {
                boundary: Boundary::Support,
                name: "cold".into(),
                region: RegionSel::Field(Arc::new(Tape::new(&g, left))),
                dofs: [true, false, false],
                value: [Held::Number(0.0), Held::Number(0.0), Held::Number(0.0)],
                film: None,
                normal: false,
                free: vec![],
                rigid: false,
            }],
            loads: vec![Load { name: "hot".into(), region: RegionSel::Field(Arc::new(Tape::new(&g, right))), total: [1.0, 0.0, 0.0], total_t: None, pressure: None, torque: None }],
            lo: [-2.0, -1.0, -1.0],
            hi: [2.0, 1.0, 1.0],
            h: 1.0,
            pressure_grid: 1.0,
            tol: 1e-9,
            support_margin: 0.0,
            modes: 0,
            density: None,
            transient: None,
            split_as: None,
            both_ways: false,
        };
        let ctx = Ctx::new(0.0, vec![]);
        let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        let sol = pool.install(|| solve_checked(&setup, &ctx)).unwrap();
        assert!(sol.report.converged, "{:?}", sol.report);
        assert_eq!(LAST_SOLVE_THREADS.load(std::sync::atomic::Ordering::Relaxed), 1);
    }
}
