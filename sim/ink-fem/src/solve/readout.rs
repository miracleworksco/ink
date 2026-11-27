//! Steps 9 and 10: what is read off the solution: the surface points stresses are
//! read at, reactions, each cell's stress or flux, the supports a stress
//! peak keeps clear of, buckling, and the held values at held nodes.

use crate::buckle::buckle;
use crate::cells::{SimCell, place};
use crate::element::{b_matrix, d_matrix, laminate_mean, shape_dn};
use crate::multigrid::Multigrid;
use crate::numerics::dot;
use crate::occupancy::{at_xi, surface_cell};
use crate::report::{Incomplete, SurfaceMeasure};
use crate::solution::NO_FIX;
use crate::system::System;
use crate::{Boundary, Held, Physics, RegionSel};
use ink_core::{V3, eval, v3};
use rayon::prelude::*;

use super::assemble::{Mesh, Model};
use super::freeing::Freed;
use super::holds::Surfaces;
use super::loads::{Loads, add_axis_moment};
use super::rigid::{Rigid, mechanisms};
use super::study::Study;
use super::{SolveError, flat_normal_hold, thin_normal_restraint};

/// The surface points stress peaks are read at: each facet's centroid,
/// in the cell holding it if that is at least half material, else in
/// the one half that cell's size under it if that is (read at its
/// nearest point). Facets half a cell apart, of the cells there: a split
/// grid cell's surface is meshed again at half its finest child's size.
pub(super) struct SurfaceReads {
    /// Each facet's cell and place in it.
    pub(super) points: Vec<(u32, V3)>,
    /// And the cell it lies in, where the peak value is read (`own`).
    pub(super) own: Vec<(u32, V3)>,
    pub(super) normals: Vec<V3>,
    /// A facet in a field fix's region: that fix's axes and value there.
    pub(super) held: Vec<Option<(u8, [f64; 3])>>,
    /// Each held one's fix and point (its value read again at another
    /// moment: `over_time`).
    pub(super) held_from: Vec<Option<(usize, V3)>>,
}

pub(super) fn surface_reads(model: &Model, sys: &System, surf: &mut Surfaces) -> SurfaceReads {
    let Model { st, holds, mesh, .. } = *model;
    let Study { s, ctx, ref grid, ref distributing, .. } = *st;
    let field = &holds.field;
    let held_geometry = &mesh.held_geometry;
    type Held = Option<((u8, [f64; 3]), (usize, V3))>;
    let (points, own, normals, held): (Vec<(u32, V3)>, Vec<(u32, V3)>, Vec<V3>, Vec<Held>) = {
        let facets = surf.facets.get_or_insert_with(|| st.whole_surface(&surf.split));
        let finest_level = sys.cells.iter().map(|c| c.level).max().unwrap_or(0);
        // (a facet in a field fix's region: that fix's axes and value there)
        let held_at = |at: V3| {
            field.iter().filter(|&&fi| !distributing[fi]).find_map(|&fi| {
                let fx = &s.fixes[fi];
                let RegionSel::Field(rt) = &fx.region else { return None };
                (eval::eval_point(rt, ctx, at) <= 0.0).then(|| (((0..3).fold(0u8, |m, d| if fx.dofs[d] { m | 1 << d } else { m }), st.values_at(fx, &[at])[0]), (fi, at)))
            })
        };
        let mut pts: Vec<((u32, V3), (u32, V3), V3, Held)> = facets
            .par_iter()
            .filter_map(|fc| {
                // (a cell under half full is read from the one under it, where
                // that holds free material: under a lip thinner than half a
                // cell beyond a hold's face lies a cell all held, whose field
                // is the held value, and the lip's face read 0 where its
                // compliance said it moved 0.0009 mm)
                let free = |ci: usize| !held_geometry[ci].as_deref().is_some_and(|g| g.0[0] == 0.0);
                let read = |takes: &dyn Fn(usize) -> bool| surface_cell(&sys.cells, &|i| model.cells_at(&sys.cells, i), grid, finest_level, fc.at, fc.normal, takes);
                let (under, own) = (read(&|ci| sys.cells[ci].frac >= 0.5), read(&|_| true));
                let stress = match (under, own) {
                    (Some(u), Some(o)) if u.0 != o.0 && !free(u.0) && free(o.0) => Some(o),
                    (u, _) => u,
                }?;
                // (A value, not a gradient, is read in the cell the point
                // lies in, its field the solve's own there: carried out from
                // the cell under it, a bar heated over its last mm read its
                // end 0.7 K high on 0.7 mm cells, that cell's slope, where
                // the heat curves the field flat, carried 0.3 mm on.)
                let own = own.unwrap_or(stress);
                Some(((stress.0 as u32, stress.1), (own.0 as u32, own.1), fc.normal, held_at(fc.at)))
            })
            .collect();
        pts.sort_by_key(|p| p.0.0);
        let mut out = (Vec::with_capacity(pts.len()), Vec::with_capacity(pts.len()), Vec::with_capacity(pts.len()), Vec::with_capacity(pts.len()));
        for (a, o, b, c) in pts {
            out.0.push(a);
            out.1.push(o);
            out.2.push(b);
            out.3.push(c);
        }
        out
    };
    let held_from = held.iter().map(|h| h.map(|h| h.1)).collect();
    SurfaceReads { points, own, normals, held: held.into_iter().map(|h| h.map(|h| h.0)).collect(), held_from }
}

/// What each fix pushes back with (`reactions`, and its moments about the
/// origin), and the free motions' amplitudes.
pub(super) struct Reactions {
    pub(super) reactions: Vec<(String, [f64; 3])>,
    pub(super) reaction_moments: Vec<[f64; 3]>,
    pub(super) free_motions: Vec<(String, Vec<(String, f64)>)>,
    pub(super) amplitudes: Vec<Option<Vec<f64>>>,
}

// 9. Reactions: what each hold's energy pushes its held value with,
// uw^T u + ww w (w = g + Phi q, q = C^+ (B^T u - Phi^T ww g) the
// amplitudes it moves by along its free motions, reported): a penalty's
// beta A (w - u), Nitsche's the free material's traction on its
// boundary with the held, with gamma (w - u) (`hold`). Each scaled as
// its cell's material is.
// (From what loads the material each holds still, which its field
// never saw: `held_load`.)
pub(super) fn reactions(model: &Model, sys: &System, u: &[f64], freed: &[Freed], loads: &Loads) -> Reactions {
    let Model { st, holds, lent, mesh, .. } = *model;
    let Study { s, ndof, ne, ref grid, ref distributing, .. } = *st;
    let clamps = &holds.clamps;
    let Mesh { ref final_of, ref final_pos, ref clamp_cells, .. } = *mesh;
    let Loads { ref held_load, ref held_moment, .. } = *loads;
    let mut reactions: Vec<(String, [f64; 3])> = s.fixes.iter().zip(held_load).map(|(f, w)| (f.name.clone(), w.map(|v| -v))).collect();
    let mut reaction_moments: Vec<[f64; 3]> = held_moment.iter().map(|m| m.map(|v| -v)).collect();
    let mut free_motions = Vec::new();
    let mut amplitudes: Vec<Option<Vec<f64>>> = vec![None; s.fixes.len()];
    for fr in freed.iter() {
        let k = fr.b.len();
        let btu: Vec<f64> = fr.b.iter().zip(&fr.moved).map(|(b, m)| dot(b, u) - m).collect();
        let q: Vec<f64> = (0..k).map(|i| (0..k).map(|j| fr.cpinv[i * k + j] * btu[j]).sum()).collect();
        let fx = &s.fixes[fr.fix];
        free_motions.push((fx.name.clone(), fx.free.iter().map(|m| m.describe()).zip(q.clone()).collect::<Vec<_>>()));
        amplitudes[fr.fix] = Some(q);
    }
    // (a distributing hold's: W P (v - u), its terms')
    for fr in freed.iter() {
        for (g, v) in fr.mean.iter().flat_map(|m| &m.held) {
            let c = v - dot(g, u);
            for (j, g) in g.iter().enumerate() {
                let force = c * g;
                reactions[fr.fix].1[j % ndof] += force;
                add_axis_moment(ndof, &mut reaction_moments[fr.fix], final_pos[j / ndof], j % ndof, force);
            }
        }
    }
    for &(cell, k) in clamp_cells.iter().filter(|cc| !distributing[clamps[cc.1].fix]) {
        let mult = sys.cells[cell].mult;
        let cl = &clamps[k];
        let fx = &s.fixes[cl.fix];
        let nodes = sys.cells[cell].nodes;
        let ue: Vec<f64> = (0..ne).map(|i| u[nodes[i / ndof] as usize * ndof + i % ndof]).collect();
        let (uw, ww) = (cl.t.kuw(), cl.t.kww());
        let mut push = cl.t.gw.clone();
        for (j, p) in push.iter_mut().enumerate() {
            *p += (0..ne).map(|i| uw[i * ne + j] * ue[i]).sum::<f64>();
        }
        if let Some(q) = &amplitudes[cl.fix] {
            let w: Vec<f64> = (0..ne).map(|i| fx.free.iter().zip(q).map(|(mo, q)| q * mo.at(final_pos[nodes[i / ndof] as usize])[i % ndof]).sum()).collect();
            for (j, p) in push.iter_mut().enumerate() {
                *p += (0..ne).map(|i| ww[j * ne + i] * w[i]).sum::<f64>();
            }
        }
        for (j, p) in push.iter().enumerate() {
            let force = mult * p;
            let pos = final_pos[nodes[j / ndof] as usize];
            reactions[cl.fix].1[j % ndof] += force;
            add_axis_moment(ndof, &mut reaction_moments[cl.fix], pos, j % ndof, force);
            if let Some(by) = cl.other {
                reactions[by].1[j % ndof] -= force;
                add_axis_moment(ndof, &mut reaction_moments[by], pos, j % ndof, -force);
            }
        }
    }
    // (and the traction lent to the cell beside a piece: `Lent`)
    for l in lent.iter() {
        let nodes = sys.cells[final_of[l.to] as usize].nodes;
        let (ctr, size) = place(grid, &sys.cells[final_of[l.from] as usize]);
        let pos = at_xi(ctr, size, l.t.xi);
        for a in 0..ndof {
            let force = model.lent_scale(sys, l) * l.w * (0..ne).map(|col| l.pf[a * ne + col] * u[nodes[col / ndof] as usize * ndof + col % ndof]).sum::<f64>();
            reactions[l.fix].1[a] += force;
            add_axis_moment(ndof, &mut reaction_moments[l.fix], pos, a, force);
        }
    }
    Reactions { reactions, reaction_moments, free_motions, amplitudes }
}

// 10. Cell-centre stress (or heat flux).
pub(super) fn field_of(st: &Study, sys: &System, c: &SimCell, u: &[f64]) -> [f64; 6] {
    let Study { s, ref grid, .. } = *st;
    let d = sys.dofs(c);
    let size = place(grid, c).1;
    let mut out = [0.0; 6];
    match s.physics {
        Physics::Elastic { .. } => {
            let b = b_matrix([0.0; 3]);
            let dm = d_matrix(c.nu);
            let ue: [f64; 24] = std::array::from_fn(|j| u[d[j] as usize]);
            let mut eps = [0.0; 6];
            for (k, e) in eps.iter_mut().enumerate() {
                *e = (0..24).map(|j| b[k][j] * ue[j]).sum::<f64>() * 2.0 / size;
            }
            if let Some(lam) = &c.lam {
                return laminate_mean(lam, true, &eps).map(|v| v * c.mult);
            }
            for (k, o) in out.iter_mut().enumerate() {
                *o = c.coef * c.mult * (0..6).map(|j| dm[k][j] * eps[j]).sum::<f64>();
            }
        }
        Physics::Heat => {
            let mut grad = [0.0; 3];
            for i in 0..8 {
                let g = shape_dn(i, [0.0; 3]);
                for a in 0..3 {
                    grad[a] += g[a] * 2.0 / size * u[d[i] as usize];
                }
            }
            if let Some(lam) = &c.lam {
                return laminate_mean(lam, false, &grad).map(|v| v * c.mult);
            }
            for a in 0..3 {
                out[a] = -c.coef * c.mult * grad[a];
            }
        }
    }
    out
}

/// The supports a stress peak keeps clear of: each cell's nearest fix
/// within the support margin (`NO_FIX`: none), and which fixes are planes
/// (a mirror's, or a thin flat restraint), cuts, or may be.
pub(super) struct Supports {
    pub(super) near_support: Vec<u16>,
    pub(super) possible_cuts: Vec<(usize, usize, f64)>,
    pub(super) unstated_holds: Vec<usize>,
    pub(super) cuts: Vec<usize>,
    pub(super) planes: Vec<bool>,
}

pub(super) fn supports(model: &Model, sys: &System, surf: &mut Surfaces) -> Supports {
    let st = model.st;
    let Study { s, ctx, ndof, h, ref grid, .. } = *model.st;
    let clamps = &model.holds.clamps;
    let centers: Vec<V3> = sys.cells.iter().map(|c| place(grid, c).0).collect();
    let mut near_support = vec![NO_FIX; sys.cells.len()];
    let mut possible_cuts = Vec::new();
    let mut unstated_holds = Vec::new();
    let mut cuts = Vec::new();
    let planes: Vec<bool> = {
        let facets = surf.facets.get_or_insert_with(|| st.whole_surface(&surf.split));
        let samples = surf.samples.get_or_insert_with(|| st.sample_faces(facets));
        let mut held_volume = vec![0.0; s.fixes.len()];
        for cl in clamps.iter() {
            held_volume[cl.fix] += cl.vol;
        }
        s.fixes
            .iter()
            .zip(held_volume)
            .enumerate()
            .map(|(i, (fx, volume))| {
                if matches!(fx.boundary, Boundary::Mirror { .. }) {
                    cuts.push(i);
                    return true;
                }
                if ndof != 3 {
                    return false;
                }
                if fx.boundary == Boundary::Unspecified && fx.value.iter().all(|v| matches!(v, Held::Number(0.0))) {
                    // Held at nothing along one direction only, as a cut of symmetry holds its plane's
                    // normal, and not said to be a physical support (`(support)`, or a pivot or a slide
                    // it leaves free, which only a support has): it may be a cut. A whole flat face so
                    // held is solved as one too (`possible_cuts`); any other is only said.
                    let one_way = fx.normal || fx.dofs.iter().filter(|d| **d).count() == 1;
                    match flat_normal_hold(fx, facets, samples, ctx, h) {
                        Some((axis, at)) => possible_cuts.push((i, axis, at)),
                        None if one_way && fx.free.is_empty() => unstated_holds.push(i),
                        None => {}
                    }
                }
                thin_normal_restraint(fx, volume, facets, samples, ctx, h)
            })
            .collect()
    };
    for (k, fx) in s.fixes.iter().enumerate().filter(|(k, _)| s.physics.ndof() == 3 && !planes[*k]) {
        for (flag, near) in near_support.iter_mut().zip(fx.region.near(ctx, &centers, s.support_margin.max(2.0 * h))) {
            if near && *flag == NO_FIX {
                *flag = k.min(NO_FIX as usize - 1) as u16;
            }
        }
    }
    Supports { near_support, possible_cuts, unstated_holds, cuts, planes }
}

// Each cell's material share inside regions that hold it entirely
// (every dof, no free motion: a load there moves nothing).
pub(super) fn fully_held(model: &Model, sys: &System) -> Vec<f64> {
    let Study { s, ndof, ref grid, .. } = *model.st;
    let (clamps, clamp_cells) = (&model.holds.clamps, &model.mesh.clamp_cells);
    let mut held = vec![0.0f64; sys.cells.len()];
    for &(cell, k) in clamp_cells.iter() {
        let fx = &s.fixes[clamps[k].fix];
        if fx.film.is_none() && !fx.normal && fx.free.is_empty() && fx.dofs[..ndof].iter().all(|&d| d) {
            held[cell] += clamps[k].vol / (sys.cells[cell].frac * place(grid, &sys.cells[cell]).1.powi(3)).max(1e-300);
        }
    }
    held
}

/// Linear buckling's results: the factor, of the loads reversed, of the
/// part as modeled (cut at its mirrors), and how many parity sectors.
pub(super) struct Buckled {
    pub(super) buckle_out: Option<crate::report::Buckle>,
    pub(super) buckle_reversed: Option<crate::report::Buckle>,
    pub(super) modeled_buckle: Option<crate::report::Buckle>,
    pub(super) parity_sectors: usize,
}

// Linear buckling of every static elastic solve: what its loads times
// which factor buckle it. Said where the loads may reach it, or where
// the eigensolve gives no factor (its claims then leave nothing out
// they cannot see). Its stress must be a solution's: not after a solve
// that broke down or stopped short.
#[allow(clippy::too_many_arguments)]
pub(super) fn buckling(
    model: &Model,
    sys: &mut System,
    mg: Multigrid,
    u: &[f64],
    res: f64,
    redesigned: bool,
    cell_field: &[[f64; 6]],
    rigid: &Rigid,
    freed: &[Freed],
    cuts: &[usize],
    incomplete: &[Incomplete],
    measures: &[SurfaceMeasure],
    mode_shapes: &mut Vec<Vec<f64>>,
) -> Result<Buckled, SolveError> {
    let Study { s, h, .. } = *model.st;
    let Mesh { ref final_idx, ref final_pos, ref held_geometry, ref clamp_cells, .. } = *model.mesh;
    let (clamps, shared_free) = (&model.holds.clamps, &model.holds.shared_free);
    let Rigid { ref pieces, ref holds, ref rotations, support_terms } = *rigid;
    let mut buckle_out = None;
    let mut buckle_reversed = None;
    let mut modeled_buckle = None;
    let mut parity_sectors = 0;
    let solved = res.is_finite() && res <= s.tol && u.iter().all(|v| v.is_finite());
    if s.modes == 0 && s.transient.is_none() && s.physics.ndof() == 3 && !redesigned && solved {
        let (mut b, shapes) = buckle(sys, &mg, u, cell_field, held_geometry, h, final_pos, rotations.iter().map(|&i| &holds[i]));
        *mode_shapes = shapes;
        parity_sectors = 1;
        // The same of the loads reversed: their solution is this one's, negated. (Not of a part cut at
        // a mirror, whose sectors each take their own.)
        if s.both_ways && cuts.is_empty() {
            let (back, stress): (Vec<f64>, Vec<[f64; 6]>) = (u.iter().map(|v| -v).collect(), cell_field.iter().map(|f| f.map(|v| -v)).collect());
            buckle_reversed = Some(buckle(sys, &mg, &back, &stress, held_geometry, h, final_pos, rotations.iter().map(|&i| &holds[i])).0);
        }
        if !cuts.is_empty() {
            modeled_buckle = Some(b.clone());
            // An incomplete base load cannot supply the prestress for the
            // reflected part. Keep the native factor for inspection.
            if !incomplete.is_empty() {
                b.failure = Some("the reflected part needs a resolved base study".into());
            } else if let Some(why) = measures.iter().find_map(SurfaceMeasure::unseen) {
                b.failure = Some(format!("the reflected part's prestress is unresolved: {why}"));
            }
            // Each reflected load and material distribution gives K and KG
            // commuting with each orthogonal reflection. Their common +/-1
            // subspaces exhaust all displacements, including repeated modes.
            // Freeze u and its stresses. Change only the homogeneous trace
            // penalty in K. All physical support terms stay in every sector.
            drop(mg);
            sys.lowrank.truncate(support_terms);
            let releases = std::mem::take(&mut sys.lowrank);
            let mut reversed = vec![false; s.fixes.len()];
            let zero = vec![0.0; sys.ndofs];
            for mask in 1..(1 << cuts.len()) {
                if b.failure.is_some() {
                    mode_shapes.clear();
                    break;
                }
                for (bit, &fi) in cuts.iter().enumerate() {
                    let flip = mask & (1 << bit) != 0;
                    if flip == reversed[fi] {
                        continue;
                    }
                    for &(cell, k) in clamp_cells.iter().filter(|(_, k)| clamps[*k].fix == fi) {
                        let ke = &mut sys.k_cut[sys.cells[cell].cut as usize];
                        for i in 0..8 {
                            for j in 0..8 {
                                for d in 0..3 {
                                    let sign = if flip == s.fixes[fi].dofs[d] { -1.0 } else { 1.0 };
                                    // (By columns: row i * 3 + d, column j * 3 + d.)
                                    ke[(j * 3 + d) * 24 + i * 3 + d] += sign * clamps[k].t.mass[i * 8 + j];
                                }
                            }
                        }
                    }
                    reversed[fi] = flip;
                }
                let mut odd_axes = 0;
                for &fi in cuts.iter() {
                    if reversed[fi]
                        && let Boundary::Mirror { axis, .. } = s.fixes[fi].boundary
                    {
                        odd_axes |= 1 << axis;
                    }
                }
                // A shared even motion has zero amplitude in an odd sector.
                // Leave that fix's original penalty intact by omitting its
                // low-rank release. Off-plane supports keep all releases.
                // (a distributing hold's: its terms holding them too)
                sys.lowrank = freed
                    .iter()
                    .flat_map(|fr| match (&fr.mean, shared_free[fr.fix] & odd_axes == 0) {
                        (_, true) => releases[fr.terms.clone()].to_vec(),
                        (Some(mean), false) => mean.all.clone(),
                        (None, false) => Vec::new(),
                    })
                    .collect();
                let (rigid, _, rotates) = mechanisms(model, pieces, sys, clamps, &zero, &reversed)?;
                sys.hold(&rigid);
                let mg = Multigrid::build(sys, final_idx);
                let (next, shapes) = buckle(sys, &mg, u, cell_field, held_geometry, h, final_pos, rotates.iter().map(|&i| &rigid[i]));
                parity_sectors += 1;
                if next.failure.is_some() {
                    b = next;
                    mode_shapes.clear();
                    break;
                }
                let low = b.low.min(next.low);
                let iterations = b.iterations + next.iterations;
                let ms = b.ms + next.ms;
                if next.factor < b.factor {
                    b = next;
                    for &fi in cuts.iter() {
                        if let Boundary::Mirror { axis, .. } = s.fixes[fi].boundary {
                            b.antisymmetric[axis] = reversed[fi];
                        }
                    }
                    *mode_shapes = shapes;
                }
                b.low = low;
                b.iterations = iterations;
                b.ms = ms;
            }
        }
        buckle_out = Some(b);
    }
    Ok(Buckled { buckle_out, buckle_reversed, modeled_buckle, parity_sectors })
}

// What is read at the nodes inside held material: the held value on
// the axes it holds, moved by its free motions (they carry the free
// field's extension past the region's edge, for the cells it cuts: read
// there, a region held at 100 C read 101.4 C). A cell all of whose
// nodes are held reads them (`Solution::cell_values`), its stress too;
// a cell the edge cuts reads its free field (with the held values, a
// bar's cells beside a 0.01 mm held layer read 9.87 MPa of 10).
pub(super) fn held_reads(model: &Model, sys: &System, u: &[f64], amplitudes: &[Option<Vec<f64>>], cell_field: &mut [[f64; 6]]) -> Vec<(u32, [f64; 3])> {
    let st = model.st;
    let Study { s, ndof, .. } = *st;
    let Mesh { ref final_pos, ref held_nodes, .. } = *model.mesh;
    let held: Vec<(u32, [f64; 3])> = if s.modes > 0 {
        Vec::new()
    } else {
        let mut v: Vec<(u32, [f64; 3])> = Vec::new();
        for &(n, fi, g) in held_nodes.iter() {
            let fx = &s.fixes[fi];
            let p = final_pos[n as usize];
            let moved = amplitudes[fi].as_ref().map(|q| fx.free.iter().zip(q).fold([0.0; 3], |a, (mo, q)| v3::add(a, v3::scale(mo.at(p), *q)))).unwrap_or([0.0; 3]);
            if v.last().is_none_or(|x| x.0 != n) {
                v.push((n, std::array::from_fn(|d| if d < ndof { u[n as usize * ndof + d] } else { 0.0 })));
            }
            let e = &mut v.last_mut().expect("a node").1;
            (0..ndof).filter(|&d| fx.dofs[d]).for_each(|d| e[d] = g[d] + moved[d]);
        }
        v
    };
    if !held.is_empty() {
        let mut read = u.to_vec();
        held.iter().for_each(|&(n, g)| (0..ndof).for_each(|d| read[n as usize * ndof + d] = g[d]));
        let is_held = |n: &u32| held.binary_search_by_key(n, |h| h.0).is_ok();
        cell_field.par_iter_mut().zip(&sys.cells).filter(|(_, c)| c.nodes.iter().all(is_held)).for_each(|(f, c)| *f = field_of(st, sys, c, &read));
    }
    held
}
