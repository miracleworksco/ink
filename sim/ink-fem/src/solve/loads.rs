//! Step 7: the load vector (`loads`): body force, pressures and region loads,
//! sampled per cell, less what material a fix holds still carries.

use crate::cells::{K, MIN_FRACTION, SimCell, masked, place};
use crate::element::shape_n;
use crate::eval::{eval_interval, eval_points};
use crate::occupancy::{at_xi, occupied, surface_cell};
use crate::report::{Footprint, Gather, Incomplete, PRESSURE_COVERAGE_LIMIT, SurfaceMeasure};
use crate::surface::{
    CUT_DEPTH, Facet, ON_SURFACE, RIM, clip_facets, face_strips, facet_cut, facet_in_cells, half_spaces, integrate_pressure, onto_surface, rim_pieces, six_points, surface_measure,
};
use crate::system::System;
use crate::tape::Tape;
use crate::{Boundary, Fix, Load, RegionSel};
use ink_core::{V3, eval, interval, v3};
use rayon::prelude::*;
use std::sync::Arc;

use super::assemble::{Mesh, Model};
use super::holds::Surfaces;
use super::study::Study;
use super::{SolveError, torque_forces};

/// The loads: each one's total and moment as applied, what fixes holding
/// material still carry of them, and what they leave for the study to say.
pub(super) struct Loads {
    pub(super) applied: Vec<(String, [f64; 3])>,
    pub(super) applied_moments: Vec<[f64; 3]>,
    /// Each fix's share of the loads on material it holds still (`held_still`).
    pub(super) held_load: Vec<[f64; 3]>,
    pub(super) held_moment: Vec<[f64; 3]>,
    /// Loads varying in time: their spread over the nodes (per unit total),
    /// scaled at each step.
    pub(super) timed: Vec<(usize, Vec<f64>)>,
    /// Where each load lands: per load, its cells and the share of it each
    /// takes (any measure: volume, a pressure's force), for what the
    /// study must say of a load in a fix or by it (`load_warnings`).
    pub(super) landed: Vec<Vec<(usize, f64)>>,
    /// A pressure's force its facets found no cell to land in: (load, lost, all).
    pub(super) pressure_lost: Vec<(String, f64, f64)>,
    /// Each load's pressure's error bound, where read.
    pub(super) pressure_error: Vec<Option<f64>>,
    /// Where each load landed, as `applied` (`Report::load_footprints`).
    pub(super) footprints: Vec<Option<Footprint>>,
}

impl Model<'_> {
    /// A grid cell's cells (its pieces, or its children's).
    pub(super) fn cells_at(&self, cells: &[SimCell], idx: [i64; 3]) -> std::ops::Range<usize> {
        let (st, cell_at) = (self.st, &self.mesh.cell_at);
        match st.in_grid(idx).then(|| cell_at[st.lin(idx)]) {
            Some(k) if k != u32::MAX => k as usize..k as usize + cells[k as usize..].iter().take_while(|c| c.idx == idx).count(),
            _ => 0..0,
        }
    }

    /// A load over material a fix holds still loads nothing: it is the fix's
    /// to carry (`held_load`, each fix's, in its reaction). Put on the cell's
    /// nodes with the rest, it pushed the free material's field where that
    /// runs on into the held: a block on a floor held over a band t thick
    /// sank rho g t^2 / 2E less than its free height says, and a bar heated
    /// all through rose (t / its free length)^2 less beyond a band held at a
    /// temperature, on every grid alike. So a body's weight and a force or a
    /// source over a region leave out the held material: the part of cell
    /// `i` (final) a fix with no free motions holds on axis `d`, in `region`
    /// (all of it with none) and weighed by `weight`: each piece's fix, place
    /// and measure. (A pivot's or a slide's material moves with it and keeps
    /// its load on the nodes; so do a torque's and a load varying in time,
    /// and a pressure on the held material's face.)
    pub(super) fn held_still(&self, i: usize, d: usize, region: Option<&RegionSel>, weight: Option<&Tape>) -> Vec<(usize, V3, f64)> {
        let Study { s, ctx, h, .. } = *self.st;
        let (tables, classified_of) = (&self.holds.tables, &self.mesh.classified_of);
        let Some(t) = tables.get(&classified_of[i]) else { return Vec::new() };
        let still = |mask: u32| (0..t.fixes.len()).find(|&b| mask >> b & 1 == 1 && s.fixes[t.fixes[b]].dofs[d] && s.fixes[t.fixes[b]].free.is_empty()).map(|b| t.fixes[b]);
        let on: Vec<(usize, &(u32, V3, f64, V3))> = t.held.iter().filter_map(|x| still(x.0).map(|fi| (fi, x))).collect();
        if on.is_empty() {
            return Vec::new();
        }
        let pts: Vec<V3> = on.iter().map(|x| x.1.3).collect();
        let mut w: Vec<f64> = on.iter().map(|x| x.1.2).collect();
        if let Some(r) = region {
            // (at the piece's own place: a point of it)
            w.iter_mut().zip(r.occupied(ctx, &pts, 1e-6 * h)).for_each(|(w, o)| *w *= o.0);
        }
        if let Some(g) = weight {
            let mut v = Vec::new();
            eval_points(g, ctx, &pts, &mut v);
            w.iter_mut().zip(v).for_each(|(w, g)| *w *= g);
        }
        on.iter().zip(w).filter(|x| x.1 != 0.0).map(|((fi, x), w)| (*fi, x.1, w)).collect()
    }

    /// A cell's material's share of each of its nodes (int N_i over it, in
    /// `region` and weighed by `weight`, with the slivers it took), and its
    /// measure.
    pub(super) fn sample_cell(&self, c: &SimCell, region: Option<&RegionSel>, weight: Option<&Tape>) -> ([f64; 8], f64) {
        let Study { s, ctx, ref grid, ref xis, ref full_xi, .. } = *self.st;
        let taken_by = &self.cells.slivers.taken_by;
        // Cells cut by the shape or by the region boundary need fine samples.
        let fine = c.frac < 1.0 || region.is_some();
        let xs = if fine { xis } else { full_xi };
        let (ctr, size) = place(grid, c);
        let pts: Vec<V3> = xs.iter().map(|&xi| at_xi(ctr, size, xi)).collect();
        let (wt, at): (Vec<f64>, Vec<V3>) =
            if c.frac < 1.0 { occupied(c.tape.as_deref().unwrap_or(&s.shape), ctx, &pts, xs, size).into_iter().unzip() } else { (vec![1.0; pts.len()], xs.clone()) };
        let mut wt = masked(wt, c.mask);
        let mut at = at;
        if let Some(r) = region {
            // The region's part of each sub-cell, offset as a fix's.
            for ((w, a), (o, off)) in wt.iter_mut().zip(at.iter_mut()).zip(r.occupied(ctx, &pts, size)) {
                *w *= o;
                *a = std::array::from_fn(|k| a[k] + off[k]);
            }
        }
        if let Some(g) = weight {
            let mut v = Vec::new();
            eval_points(g, ctx, &pts, &mut v);
            for (w, g) in wt.iter_mut().zip(v) {
                *w *= g;
            }
        }
        let dv = size * size * size / pts.len() as f64;
        let mut nodal = [0.0; 8];
        let mut total = 0.0;
        for (xi, w) in at.iter().zip(&wt) {
            if *w == 0.0 {
                continue;
            }
            total += w * dv;
            for (k, n) in nodal.iter_mut().enumerate() {
                *n += w * dv * shape_n(k, *xi);
            }
        }
        // The slivers it took (`taken_by`): their material's share too, on
        // its nodes, its shape functions carried on to where they lie.
        for (sl, from) in taken_by.get(&(c.idx, c.level, c.sub)).into_iter().flatten() {
            let (ctr, size) = place(grid, sl);
            let pts: Vec<V3> = xis.iter().map(|&xi| at_xi(ctr, size, xi)).collect();
            let (wt, at): (Vec<f64>, Vec<V3>) = occupied(sl.tape.as_deref().unwrap_or(&s.shape), ctx, &pts, xis, size).into_iter().unzip();
            let (mut wt, mut at) = (masked(wt, sl.mask), at);
            if let Some(r) = region {
                for ((w, a), (o, off)) in wt.iter_mut().zip(at.iter_mut()).zip(r.occupied(ctx, &pts, size)) {
                    *w *= o;
                    *a = std::array::from_fn(|k| a[k] + off[k]);
                }
            }
            if let Some(g) = weight {
                let mut v = Vec::new();
                eval_points(g, ctx, &pts, &mut v);
                wt.iter_mut().zip(v).for_each(|(w, g)| *w *= g);
            }
            let dv = size * size * size / pts.len() as f64;
            for (xi, w) in at.iter().zip(&wt).filter(|x| *x.1 != 0.0) {
                total += w * dv;
                let xi: V3 = std::array::from_fn(|a| xi[a] + 2.0 * from[a] as f64);
                for (k, n) in nodal.iter_mut().enumerate() {
                    *n += w * dv * shape_n(k, xi);
                }
            }
        }
        (nodal, total)
    }
}

/// A moment's component of a force `val` on axis `d` at `p`, added to `m` (elastic alone).
pub(super) fn add_axis_moment(ndof: usize, m: &mut [f64; 3], p: V3, d: usize, val: f64) {
    if ndof == 3 {
        let mut fv = [0.0; 3];
        fv[d] = val;
        let c = v3::cross(p, fv);
        for k in 0..3 {
            m[k] += c[k];
        }
    }
}

// A fix whose region is material it holds still: not a mirror's (the
// half of space its part is in: it holds its plane), a film's, a
// surface's normal alone, or one with motions left free.
pub(super) fn holds_still(fx: &Fix) -> bool {
    fx.free.is_empty() && fx.film.is_none() && !fx.normal && !matches!(fx.boundary, Boundary::Mirror { .. })
}

// 7. Load vector: body force and region loads, sampled per cell.
pub(super) fn loads(model: &Model, sys: &mut System, f: &mut [f64], surf: &mut Surfaces, incomplete: &mut Vec<Incomplete>, dropped: (usize, f64)) -> Result<Loads, SolveError> {
    let s = model.st.s;
    let mut out = Loads {
        applied: Vec::new(),
        applied_moments: Vec::new(),
        held_load: vec![[0.0; 3]; s.fixes.len()],
        held_moment: vec![[0.0; 3]; s.fixes.len()],
        timed: Vec::new(),
        landed: vec![Vec::new(); s.loads.len()],
        pressure_lost: Vec::new(),
        pressure_error: Vec::new(),
        footprints: Vec::new(),
    };
    if let Some((rho, g)) = &s.body {
        body_load(model, sys, f, &mut out, rho, g);
    }
    out.pressure_error = vec![None; out.applied.len() + s.loads.len()];
    for (li, l) in s.loads.iter().enumerate() {
        if let Some(pt) = &l.pressure {
            pressure_load(model, sys, f, surf, incomplete, &mut out, li, l, pt)?;
            continue;
        }
        region_load(model, sys, f, &mut out, li, l, dropped)?;
    }
    if !out.timed.is_empty() && s.transient.is_none() {
        return Err("loads that depend on t need a (transient ...) study".into());
    }
    Ok(out)
}

/// The body's weight (`Setup::body`: density `rho`, gravity `g`).
fn body_load(model: &Model, sys: &System, f: &mut [f64], out: &mut Loads, rho: &Arc<Tape>, g: &[f64; 3]) {
    let Model { st, cells, mesh, .. } = *model;
    let final_pos = &mesh.final_pos;
    let Study { s, ctx, ndof, ref grid, ref xis, .. } = *st;
    let slivers_left = &cells.slivers.left;
    let Loads { ref mut applied, ref mut applied_moments, ref mut held_load, ref mut held_moment, ref mut footprints, .. } = *out;
    let parts: Vec<([f64; 8], f64)> = sys.cells.par_iter().map(|c| model.sample_cell(c, None, Some(rho))).collect();
    let mut total = [0.0; 3];
    let mut total_m = [0.0; 3];
    // (the part, its weight at its centre of mass)
    let mut footprint = Gather::new();
    footprint.volume = mesh.volume;
    for (c, (nodal, _)) in sys.cells.iter().zip(&parts) {
        let (ctr, size) = place(grid, c);
        footprint.cell(ctr, size);
        for k in 0..8 {
            let p = final_pos[c.nodes[k] as usize];
            footprint.weigh(p, nodal[k]);
            for d in 0..ndof {
                let v = nodal[k] * g[d];
                f[c.nodes[k] as usize * ndof + d] += v;
                total[d] += v;
                add_axis_moment(ndof, &mut total_m, p, d, v);
            }
        }
    }
    applied.push(("body".to_string(), total));
    applied_moments.push(total_m);
    footprints.push(footprint.footprint((s.lo, s.hi)));
    for (i, c) in sys.cells.iter().enumerate() {
        let (ctr, size) = place(grid, c);
        for d in (0..ndof).filter(|&d| g[d] != 0.0) {
            for (fi, xi, w) in model.held_still(i, d, None, Some(rho)) {
                for k in 0..8 {
                    f[c.nodes[k] as usize * ndof + d] -= g[d] * w * shape_n(k, xi);
                }
                held_load[fi][d] += g[d] * w;
                add_axis_moment(ndof, &mut held_moment[fi], at_xi(ctr, size, xi), d, g[d] * w);
            }
        }
    }
    // A sliver still dropped (`slivers_left`) has no cell to weigh on:
    // what of it a fix holds still is that fix's to carry, as a held
    // cell's is, and counts in the body's weight. (A bar held over its
    // first three cells, its sides' cells slivers, read its fix's
    // reaction 0.11% short of its weight.)
    let mut dropped_weight = [0.0; 3];
    let mut dropped_moment = [0.0; 3];
    for sl in slivers_left {
        let (ctr, size) = place(grid, sl);
        let pts: Vec<V3> = xis.iter().map(|&xi| at_xi(ctr, size, xi)).collect();
        let mut left = masked(occupied(sl.tape.as_deref().unwrap_or(&s.shape), ctx, &pts, xis, size).into_iter().map(|o| o.0).collect(), sl.mask);
        let mut density = Vec::new();
        eval_points(rho, ctx, &pts, &mut density);
        let dv = size * size * size / pts.len() as f64;
        for (fi, fx) in s.fixes.iter().enumerate().filter(|(_, fx)| holds_still(fx)) {
            for (k, (share, _)) in fx.region.occupied(ctx, &pts, size).into_iter().enumerate() {
                if !(share > 0.0 && left[k] > 0.0) {
                    continue;
                }
                let w = left[k] * share.min(1.0) * density[k] * dv;
                for d in (0..ndof).filter(|&d| g[d] != 0.0 && fx.dofs[d]) {
                    held_load[fi][d] += g[d] * w;
                    dropped_weight[d] += g[d] * w;
                    add_axis_moment(ndof, &mut held_moment[fi], pts[k], d, g[d] * w);
                    add_axis_moment(ndof, &mut dropped_moment, pts[k], d, g[d] * w);
                }
                // (held on every axis: no other fix's)
                if (0..ndof).all(|d| fx.dofs[d]) {
                    left[k] -= left[k] * share.min(1.0);
                }
            }
        }
    }
    if let Some(body) = applied.last_mut() {
        body.1.iter_mut().zip(dropped_weight).for_each(|(t, w)| *t += w);
    }
    if let Some(body_m) = applied_moments.last_mut() {
        body_m.iter_mut().zip(dropped_moment).for_each(|(t, w)| *t += w);
    }
}

/// A pressure load `l` (its field `pt`), over the facets of the part's surface.
#[allow(clippy::too_many_arguments)]
fn pressure_load(
    model: &Model,
    sys: &System,
    f: &mut [f64],
    surf: &mut Surfaces,
    incomplete: &mut Vec<Incomplete>,
    out: &mut Loads,
    li: usize,
    l: &Load,
    pt: &Arc<Tape>,
) -> Result<(), SolveError> {
    let Model { st, cells, holds, mesh, .. } = *model;
    let Study { s, ctx, ndof, h, dims, ref grid, .. } = *st;
    let Mesh { ref final_pos, ref held_geometry, ref classified_of, .. } = *mesh;
    let taken = &cells.slivers.taken;
    let Loads {
        ref mut applied,
        ref mut applied_moments,
        ref mut held_load,
        ref mut held_moment,
        ref mut landed,
        ref mut pressure_lost,
        ref mut pressure_error,
        ref mut footprints,
        ..
    } = *out;
    // One sample per facet, distributed to its cell's nodes: every
    // sample's load lands (else load goes missing unevenly and bends
    // thin shells).
    // Each facet's share of the region, at the share's centroid.
    let (pts, samples) = if let Some((facets, points)) = &surf.pressure {
        (facets, points)
    } else {
        let facets = surf.facets.get_or_insert_with(|| st.whole_surface(&surf.split));
        let points = surf.samples.get_or_insert_with(|| st.sample_faces(facets));
        (&*facets, &*points)
    };
    let exact_surface = surf.pressure.is_some();
    // Clip a planar load region before splitting the facet into
    // cells and holds. Its partial facets then carry pressure
    // by the same quadrature as complete facets.
    let domain = std::array::from_fn(|a| interval::I::new(s.lo[a], s.hi[a]));
    let clipped = match &l.region {
        RegionSel::Field(t) if exact_surface => {
            let mid = std::array::from_fn(|a| domain[a].lo + 0.5 * (domain[a].hi - domain[a].lo));
            let bounds = domain.map(|r| interval::I::new(r.lo - h, r.hi + h));
            half_spaces(t, ctx, mid, bounds).map(|planes| clip_facets(pts, mid, &planes, 1e-9 * h))
        }
        _ => None,
    };
    let (pts, mut inside) = match &clipped {
        Some(facets) => (facets, facets.iter().map(|f| (1.0, f.at)).collect()),
        None => (pts, l.region.on_surface(ctx, pts, samples, h)),
    };
    // A facet across a rim (its four samples on faces turned more than RIM apart: an end face and
    // the wall it meets) is no plane of the surface: its pieces, each projected along its normal
    // onto the surface, pressed where the region holds that point (`rim_pieces`). The facet's own
    // share and normal put a rim's pressure half across it and short (a tube's end 0.14% short at
    // 1 mm cells; a region a micron under it, on its own, 35%).
    let (rims, strips) = if exact_surface {
        (vec![None; pts.len()], vec![Vec::new(); pts.len()])
    } else {
        let normals: Vec<Option<V3>> = crate::eval::eval_grads(&s.shape, ctx, samples).iter().map(|g| v3::direction([g[1], g[2], g[3]])).collect();
        let rims: Vec<_> = (0..pts.len())
            .into_par_iter()
            .map(|i| {
                let n = &normals[4 * i..4 * i + 4];
                let turned = n.iter().flatten().any(|a| n.iter().flatten().any(|b| v3::dot(*a, *b) < RIM));
                turned.then(|| rim_pieces(&pts[i], &s.shape, ctx, &l.region, h))
            })
            .collect();
        for (k, pieces) in rims.iter().enumerate() {
            if let Some(pieces) = pieces {
                let held: f64 = pieces.iter().map(|p| p.1).sum();
                inside[k] = if held > 0.0 { (held / pts[k].area, std::array::from_fn(|a| pieces.iter().map(|p| p.0[a] * p.1).sum::<f64>() / held)) } else { (0.0, pts[k].at) };
            }
        }
        // A facet the region holds in part has its share already: its edge is no end of the face.
        let whole: Vec<bool> = inside.iter().zip(&rims).map(|(x, r)| x.0 >= 1.0 && r.is_none()).collect();
        let strips = face_strips(pts, &whole, &s.shape, ctx, &l.region, h);
        (rims, strips)
    };
    if inside.iter().all(|x| x.0 <= 0.0) {
        return Err(format!(
            "load '{}' (a pressure) touches none of the part's surface (its region and the shape's surface do not meet): give it a region \
                     that holds the faces it presses",
            l.name
        )
        .into());
    }
    let pos: Vec<V3> = inside.iter().map(|x| x.1).collect();
    let pv = crate::eval::eval_points_par(pt, ctx, &pos);
    let mut total = [0.0; 3];
    let mut total_m = [0.0; 3];
    let finest = sys.cells.iter().map(|c| c.level).max().unwrap_or(0);
    let range = eval_interval(pt, ctx, domain, None);
    let mut worst = 0.0;
    let mut unresolved = None;
    let mut dropped_lip_load: Option<(usize, V3, f64)> = None;
    // A line crossing can fall just off an exactly selected plane.
    // Where every loaded planar facet is proved inside the region,
    // its area cannot be an overestimate from sampling the mask.
    let exact_mask = match &l.region {
        RegionSel::Field(t) if exact_surface => Some(t),
        _ => None,
    };
    let mut sampled_mask = exact_mask.is_none();
    let mut force_samples = Vec::new();
    let (mut meant, mut lost, mut covered) = (0.0, 0.0, Vec::new());
    // (each sample it put on the part: where, its pressure, and its force)
    let mut footprint = Gather::new();
    for (i, (fc, (&(reg, q), pr))) in pts.iter().zip(inside.iter().zip(&pv)).enumerate() {
        let (nrm, w) = (&fc.normal, &fc.area);
        if reg <= 0.0 {
            continue;
        }
        // Cover the pressure between the quadrature points, even
        // when all samples of a narrow ridge land on its background.
        let bounds = || {
            std::array::from_fn(|a| {
                let xs = fc.corners.map(|p| p[a]);
                interval::I::new(xs.into_iter().fold(q[a], f64::min), xs.into_iter().fold(q[a], f64::max))
            })
        };
        if !sampled_mask {
            sampled_mask = exact_mask.is_some_and(|mask| eval_interval(mask, ctx, bounds(), None).hi > 0.0);
        }
        let variation = if range.lo != range.hi {
            let p = eval_interval(pt, ctx, bounds(), None);
            (p.hi - pr).abs().max((p.lo - pr).abs()) * w * reg
        } else {
            0.0
        };
        if *pr != 0.0 && st.cut_face(fc).is_some() {
            return Err(SolveError::InvalidMirror(format!(
                "pressure '{}' loads a symmetry cut that is internal to the reflected part. Select only the exposed faces or model the whole part",
                l.name
            )));
        }
        if *pr == 0.0 && variation > 0.0 && st.cut_face(fc).is_some() {
            return Err(format!("cannot establish zero pressure '{}' on an internal mirror face. Select only the exposed faces or refine (res r)", l.name).into());
        }
        // A facet within its cells (`facet_in_cells`): a cell's shape
        // functions bend at its planes, and a facet's force put whole
        // at its centre gave the nodes beside the plane it crossed the
        // wrong shares. They even out along a face and not at its
        // edges: a 20 x 4 x 4 bar pulled 10 MPa on its end, its sides
        // 0.2 mm past the planes of 0.3 mm cells, moved the end's
        // corner 0.067% further than its middle, and three grids
        // agreed on it to 0.047%. (A facet the load's region ends
        // within is put on by its share, at the share's centre, as
        // it was.)
        let leaf = |p: V3| surface_cell(&sys.cells, &|i| model.cells_at(&sys.cells, i), grid, finest, p, fc.normal, &|_| true).map(|(ci, _)| place(grid, &sys.cells[ci]).1);
        // (A fix's region is asked on the surface, where a facet's
        // point lies once carried there: the facet's corners carried
        // onto it along its normal, and each point by its place among
        // them. A facet beside an edge lies a little inside the part,
        // out of a region drawn within a hair of the face.)
        let on = std::cell::OnceCell::new();
        let carried = |x: V3| -> V3 {
            if exact_surface {
                return x;
            }
            let on: &Vec<V3> = on.get_or_init(|| onto_surface(&s.shape, ctx, &fc.corners, fc.normal, fc.step.max(1e-9 * h)));
            let [a, b, c] = fc.corners;
            let (u, v, w) = (v3::sub(b, a), v3::sub(c, a), v3::sub(x, a));
            let (uu, uv, vv, wu, wv) = (v3::dot(u, u), v3::dot(u, v), v3::dot(v, v), v3::dot(w, u), v3::dot(w, v));
            let den = uu * vv - uv * uv;
            let (p, q) = ((vv * wu - uv * wv) / den, (uu * wv - uv * wu) / den);
            std::array::from_fn(|k| (1.0 - p - q) * on[0][k] + p * on[1][k] + q * on[2][k])
        };
        let pieces: Option<Vec<Facet>> = (rims[i].is_none() && reg == 1.0).then(|| facet_in_cells(fc, grid.origin, grid.h, &leaf));
        // A piece on the face of material a fix holds still is that
        // fix's to carry, as a weight over it is (`held_load`): each
        // piece cut by each such fix's region. On its cell's nodes
        // it moved the free material beside it: a 3 x 3 x 3 block
        // clamped over z < 0.1 and pressed on its side moved 3.2e-4
        // less than pressed over z > 0.1 alone.
        let pieces: Option<Vec<(Facet, Option<usize>)>> = pieces.map(|pieces| {
            let mut free: Vec<Facet> = pieces;
            let mut out: Vec<(Facet, Option<usize>)> = Vec::new();
            for (fi, fx) in s.fixes.iter().enumerate().filter(|(_, fx)| holds_still(fx)) {
                let RegionSel::Field(rt) = &fx.region else { continue };
                // (its region near the facet at all: by the field's bounds over the facet's box, a step wider)
                let reach: [interval::I; 3] = std::array::from_fn(|a| {
                    let xs = fc.corners.map(|p| p[a]);
                    interval::I::new(xs.into_iter().fold(f64::INFINITY, f64::min) - fc.step, xs.into_iter().fold(f64::NEG_INFINITY, f64::max) + fc.step)
                });
                if eval_interval(rt, ctx, reach, None).lo > 0.0 {
                    continue;
                }
                let (mut within, mut without) = (Vec::new(), Vec::new());
                for piece in &free {
                    facet_cut(piece, &|x| eval::eval_point(rt, ctx, carried(x)), CUT_DEPTH, &mut within, &mut without);
                }
                out.extend(within.into_iter().map(|p| (p, Some(fi))));
                free = without;
            }
            out.extend(free.into_iter().map(|p| (p, None)));
            out
        });
        let area = w * reg;
        let sampled = pr.abs() * area;
        let tol = s.tol.min(1e-8);
        let mut error = variation;
        force_samples.clear();
        // (each sample's fix, where one holds the material it is on still)
        let mut held_by: Vec<Option<usize>> = Vec::new();
        // A coverage cutoff is not a response-error bound. Even
        // a small invisible ripple can reverse a close claim. For
        // a full facet within the cutoff, integrate the pressure
        // to the solve tolerance or leave its response unknown.
        if let Some(pieces) = &rims[i] {
            // A rim's pieces, each at its own point on the surface with the pressure there.
            let at: Vec<V3> = pieces.iter().map(|p| p.0).collect();
            let pv = crate::eval::eval_points_par(pt, ctx, &at);
            force_samples.extend(pieces.iter().zip(pv).map(|(p, v)| (p.0, v, p.1)));
        } else if let Some(pieces) = &pieces {
            // Each piece by the six-point rule, exact for a cell's
            // shape function on a planar piece within it.
            let from = |piece: &(Facet, Option<usize>), samples: &mut Vec<(V3, f64, f64)>, held_by: &mut Vec<Option<usize>>, varying: bool| -> bool {
                let ok = if varying {
                    integrate_pressure(pt, ctx, &piece.0, tol * pr.abs(), error <= PRESSURE_COVERAGE_LIMIT * sampled, samples)
                } else {
                    samples.extend(six_points(&piece.0).map(|(at, weight)| (at, *pr, weight)));
                    true
                };
                held_by.resize(samples.len(), piece.1);
                ok
            };
            if error <= tol * sampled {
                pieces.iter().for_each(|piece| {
                    from(piece, &mut force_samples, &mut held_by, false);
                });
            } else if pieces.iter().all(|piece| from(piece, &mut force_samples, &mut held_by, true)) {
                error = tol * sampled;
            } else {
                force_samples.clear();
                held_by.clear();
                force_samples.push((q, *pr, area));
            }
        } else {
            force_samples.push((q, *pr, area));
        }
        // (The facet's own samples: a point of them over a cell
        // dropped as a sliver goes to the cell the facet's centre
        // has, as the whole facet did; a rim's or a strip's that no
        // cell takes is missed, and said.)
        let own = if rims[i].is_none() { force_samples.len() } else { 0 };
        // The face past the mesh's end at its edges (`face_strips`).
        if !strips[i].is_empty() {
            let at: Vec<V3> = strips[i].iter().map(|p| p.0).collect();
            force_samples.extend(strips[i].iter().zip(crate::eval::eval_points_par(pt, ctx, &at)).map(|(p, v)| (p.0, v, p.1)));
        }
        let mut missed = 0.0;
        let mut landed_area = 0.0;
        // Each sample presses along the surface's own normal there, on the surface's own area:
        // a facet's patch dA projects to n_s dA / (n_s . n_f) of the surface. A facet that chamfers
        // a rim (a curved face meeting a flat one) took its region's share along its own normal,
        // half across the rim, and on its chord's area, short of the face's: LE1's ellipse took
        // its 10 MPa 1.8% short in x and 2% in y, and put 0.16 MN along z, where the faces are
        // not pressed. Where the surface turns too far from the facet the facet's own is kept.
        let mut grads = if exact_surface { Vec::new() } else { crate::eval::eval_grads(&s.shape, ctx, &force_samples.iter().map(|x| x.0).collect::<Vec<_>>()) };
        // (A facet beside a sharp edge lies a little inside the part,
        // its corners drawn in, and a point of it there is nearer the
        // other face: the field's normal is that face's, turned too
        // far, and the facet's own, tilted a hundredth of a radian,
        // was taken, a hundredth of its force sideways. The surface's
        // normal is read a little way out along the facet's own
        // instead, where the facet's face is the nearer. A bar's
        // pulled end took 7e-4 N sideways of 40, enough to bend 18.6
        // mm of it: its end's edge read 0.035% too far.)
        let turned = |g: &[f64; 4]| !v3::direction([g[1], g[2], g[3]]).is_some_and(|n| v3::dot(n, fc.normal) > ON_SURFACE);
        for out in [0.005, 0.02, 0.08] {
            let again: Vec<usize> = (0..grads.len()).filter(|&k| turned(&grads[k])).collect();
            if again.is_empty() {
                break;
            }
            let read = crate::eval::eval_grads(&s.shape, ctx, &again.iter().map(|&k| v3::add(force_samples[k].0, v3::scale(fc.normal, out * fc.step))).collect::<Vec<_>>());
            for (&k, g) in again.iter().zip(&read) {
                if !turned(g) {
                    grads[k] = *g;
                }
            }
        }
        let mut centre: Option<Option<(usize, V3, f64)>> = None;
        for (k, &(p, pressure, patch)) in force_samples.iter().enumerate() {
            let (nrm, patch) = match grads.get(k).and_then(|g| v3::direction([g[1], g[2], g[3]])) {
                Some(n) if v3::dot(n, *nrm) > ON_SURFACE => (n, patch / v3::dot(n, *nrm)),
                _ => (*nrm, patch),
            };
            let nrm = &nrm;
            let force = pressure.abs() * patch;
            meant += force;
            if force == 0.0 {
                landed_area += patch;
                footprint.press(p, pressure, patch, force);
                continue;
            }
            // (on material a fix holds still: that fix's, on the axes it holds)
            let still = held_by.get(k).copied().flatten();
            if let Some(fi) = still {
                for d in (0..ndof.min(3)).filter(|&d| s.fixes[fi].dofs[d]) {
                    let component = -pressure * nrm[d] * patch;
                    held_load[fi][d] += component;
                    total[d] += component;
                    add_axis_moment(ndof, &mut held_moment[fi], p, d, component);
                    add_axis_moment(ndof, &mut total_m, p, d, component);
                }
                if (0..ndof.min(3)).all(|d| s.fixes[fi].dofs[d]) {
                    landed_area += patch;
                    footprint.press(p, pressure, patch, force);
                    continue;
                }
            }
            // (A point on the surface at a corner the cells round may lie in none along the surface's
            // normal: then along the facet's, as the facet's own sample was found.)
            let found = surface_cell(&sys.cells, &|i| model.cells_at(&sys.cells, i), grid, finest, p, *nrm, &|_| true)
                .or_else(|| surface_cell(&sys.cells, &|i| model.cells_at(&sys.cells, i), grid, finest, p, fc.normal, &|_| true))
                // (over a sliver a cell beside it took: that cell's, at the place carried on)
                .or_else(|| {
                    (0..=finest).rev().find_map(|level| {
                        let n = 1i64 << level;
                        let size = grid.h / n as f64;
                        // (the cell just inside the face there, at this size)
                        let g: [i64; 3] = std::array::from_fn(|a| ((p[a] - 1e-6 * size * fc.normal[a] - grid.origin[a]) / size).floor() as i64);
                        let (to, _) = taken.get(&(g.map(|v| v.div_euclid(n)), level, g.map(|v| v.rem_euclid(n) as u16)))?;
                        let ci = model.cells_at(&sys.cells, to.0).find(|&i| sys.cells[i].level == to.1 && sys.cells[i].sub == to.2)?;
                        let (ctr, size) = place(grid, &sys.cells[ci]);
                        Some((ci, std::array::from_fn(|a| (p[a] - ctr[a]) * 2.0 / size)))
                    })
                })
                .or_else(|| {
                    let (ci, xi, size) = (k < own)
                        .then(|| {
                            *centre.get_or_insert_with(|| {
                                surface_cell(&sys.cells, &|i| model.cells_at(&sys.cells, i), grid, finest, q, fc.normal, &|_| true)
                                    .map(|(ci, xi)| (ci, xi, place(grid, &sys.cells[ci]).1))
                            })
                        })
                        .flatten()?;
                    Some((ci, std::array::from_fn(|a| xi[a] + (p[a] - q[a]) * 2.0 / size)))
                });
            let Some((ci, xi)) = found else {
                missed += force;
                continue;
            };
            landed_area += patch;
            footprint.press(p, pressure, patch, force);
            landed[li].push((ci, force));
            if dropped_lip_load.is_none()
                && s.modes == 0
                && let Some(g) = held_geometry[ci].as_deref()
            {
                let at = if k < own && pieces.is_some() { carried(p) } else { p };
                for d in (0..ndof.min(3)).filter(|&d| !still.is_some_and(|fi| s.fixes[fi].dofs[d]) && nrm[d].abs() > 1e-6 && g.2[d] && g.1[d][0] == 0.0) {
                    if let Some(fi) = holds.held_free_at(st, classified_of[ci], d, at) {
                        dropped_lip_load = Some((fi, at, place(grid, &sys.cells[ci]).1));
                        break;
                    }
                }
            }
            let c = &sys.cells[ci];
            for d in (0..ndof.min(3)).filter(|&d| !still.is_some_and(|fi| s.fixes[fi].dofs[d])) {
                let component = -pressure * nrm[d] * patch;
                for k in 0..8 {
                    let v = shape_n(k, xi) * component;
                    f[c.nodes[k] as usize * ndof + d] += v;
                    add_axis_moment(ndof, &mut total_m, final_pos[c.nodes[k] as usize], d, v);
                }
                total[d] += component;
            }
        }
        lost += missed;
        if landed_area > 0.0 {
            covered.push((*fc, landed_area / w));
        }
        let fraction = if error + missed == 0.0 {
            0.0
        } else if sampled > 0.0 {
            (error + missed) / sampled
        } else {
            f64::INFINITY
        };
        if fraction > worst {
            worst = fraction;
        }
        if (missed > 0.0 || error > tol * sampled) && s.modes == 0 && unresolved.as_ref().is_none_or(|(old, _)| fraction > *old) {
            unresolved = Some((fraction, Incomplete::Pressure { applied: applied.len(), at: q, lost: missed, variation: error, meant: sampled }));
        }
    }
    pressure_error[applied.len()] = Some(worst);
    incomplete.extend(unresolved.map(|(_, why)| why));
    if let Some((fi, at, size)) = dropped_lip_load
        && !incomplete.iter().any(|e| matches!(e, Incomplete::Held(_)))
    {
        incomplete.push(Incomplete::Held(format!(
                    "fix '{}' drops or absorbs free material under {}% of a cell near {:?} mm, and load '{}' acts on it: the cells cannot carry what it puts there, so the study is not settled. To settle it, write the hold's region out to the part's face, or refine (res r) to at least {:.4} mm",
                    s.fixes[fi].name,
                    100.0 * MIN_FRACTION,
                    at.map(|v| (v * 1000.0).round() / 1000.0),
                    l.name,
                    0.5 * size
                )));
    }
    if lost > 1e-3 * meant {
        pressure_lost.push((l.name.clone(), lost, meant));
    }
    let mut measure = SurfaceMeasure::of(&l.name, "load", &covered, surface_measure(&s.shape, ctx, grid, dims, &sys.cells, &surf.split, &l.region, s.lo), grid);
    if !sampled_mask {
        measure.excess = 0.0;
    }
    measure.unused_pressure = s.modes > 0 || (range.lo == 0.0 && range.hi == 0.0);
    surf.measures.push(measure);
    applied.push((l.name.clone(), total));
    applied_moments.push(total_m);
    footprints.push(footprint.footprint((s.lo, s.hi)));
    Ok(())
}

/// A load `l` over a region: a force, a torque, a source, or one varying in time.
fn region_load(model: &Model, sys: &System, f: &mut [f64], out: &mut Loads, li: usize, l: &Load, (dropped_cells, dropped_volume): (usize, f64)) -> Result<(), SolveError> {
    let Model { st, cells, mesh, .. } = *model;
    let Mesh { ref final_pos, nn, .. } = *mesh;
    let taken_by = &cells.slivers.taken_by;
    let Study { s, ndof, h, ref grid, .. } = *st;
    let Loads { ref mut applied, ref mut applied_moments, ref mut held_load, ref mut held_moment, ref mut timed, ref mut landed, ref mut footprints, .. } = *out;
    let parts: Vec<Option<([f64; 8], f64)>> = sys
        .cells
        .par_iter()
        .map(|c| {
            // (its own cell, or a sliver's it took)
            let reached =
                st.region_overlaps(&l.region, c.idx) || taken_by.get(&(c.idx, c.level, c.sub)).is_some_and(|v| v.iter().any(|(sl, _)| st.region_overlaps(&l.region, sl.idx)));
            reached.then(|| model.sample_cell(c, Some(&l.region), None))
        })
        .collect();
    let vol: f64 = parts.iter().flatten().map(|p| p.1).sum();
    if vol <= 0.0 {
        if dropped_cells > 0 {
            return Err(format!(
                "load '{}' region holds no material joined to a support: {dropped_cells} cells ({:.4} mm3) joined to no fix were dropped as \
                     floating, and the load lies on none of the rest (a part parted from the one held: a pin in a clearance hole, a gap the \
                     cells part). Join it to the held part with material, or hold it",
                l.name, dropped_volume
            )
            .into());
        }
        // Cells its bounds cannot rule out, none of whose samples lies in it: a region smaller than
        // the samples, not one without material (a 0.1 mm cube read on 2 mm cells, as the cells fell).
        let near = parts.iter().flatten().count();
        if near > 0 {
            return Err(format!(
                "load '{}' region lies in {near} cells but on none of their samples ({K} a side, {} mm apart on cells of {} mm): it is smaller \
                     than the cells can see where they fall. Refine (res r) until a sample lies in it, or give the region the size of what presses there",
                l.name,
                h / K as f64,
                h
            )
            .into());
        }
        return Err(format!("load '{}' region contains no material", l.name).into());
    }
    landed[li] = parts.iter().enumerate().filter_map(|(ci, p)| p.as_ref().map(|p| (ci, p.1))).filter(|p| p.1 > 0.0).collect();
    // Its material: the cells it lands in, and its centroid, where its force acts (its moment's arm).
    let mut footprint = Gather::new();
    for (c, p) in sys.cells.iter().zip(&parts) {
        let Some((nodal, v)) = p.as_ref().filter(|p| p.1 > 0.0) else { continue };
        let (ctr, size) = place(grid, c);
        footprint.cell(ctr, size);
        footprint.volume += v;
        for k in 0..8 {
            footprint.weigh(final_pos[c.nodes[k] as usize], nodal[k]);
        }
    }
    footprints.push(footprint.footprint((s.lo, s.hi)));
    if let Some(tq) = l.torque {
        // The region's material as weighted nodes, and the forces of a
        // turn of it that make the torque (`torque_forces`).
        let mut wn = vec![0.0; nn];
        for (c, p) in sys.cells.iter().zip(&parts) {
            let Some((nodal, _)) = p else { continue };
            for k in 0..8 {
                wn[c.nodes[k] as usize] += nodal[k];
            }
        }
        let nodes: Vec<(u32, f64)> = (0..nn as u32).map(|n| (n, wn[n as usize])).filter(|n| n.1 != 0.0).collect();
        let mut total_m = [0.0; 3];
        for (&(n, _), fv) in nodes.iter().zip(torque_forces(&nodes, |n| final_pos[n as usize], tq)) {
            let pos = final_pos[n as usize];
            for d in 0..3 {
                f[n as usize * ndof + d] += fv[d];
                add_axis_moment(ndof, &mut total_m, pos, d, fv[d]);
            }
        }
        applied.push((l.name.clone(), [0.0; 3]));
        applied_moments.push(total_m);
        return Ok(());
    }
    if l.total_t.is_some() {
        let mut unit = vec![0.0; nn];
        for (c, p) in sys.cells.iter().zip(&parts) {
            let Some((nodal, _)) = p else { continue };
            for k in 0..8 {
                unit[c.nodes[k] as usize] += nodal[k] / vol;
            }
        }
        sys.fold(&mut unit, 1);
        timed.push((li, unit));
        applied.push((l.name.clone(), [0.0; 3]));
        applied_moments.push([0.0; 3]);
        return Ok(());
    }
    let mut total_m = [0.0; 3];
    for (i, (c, p)) in sys.cells.iter().zip(&parts).enumerate() {
        let Some((nodal, _)) = p else { continue };
        let (ctr, size) = place(grid, c);
        for k in 0..8 {
            let pos = final_pos[c.nodes[k] as usize];
            for d in 0..ndof {
                let v = l.total[d] * nodal[k] / vol;
                f[c.nodes[k] as usize * ndof + d] += v;
                add_axis_moment(ndof, &mut total_m, pos, d, v);
            }
        }
        for d in (0..ndof).filter(|&d| l.total[d] != 0.0) {
            for (fi, xi, w) in model.held_still(i, d, Some(&l.region), None) {
                for k in 0..8 {
                    f[c.nodes[k] as usize * ndof + d] -= l.total[d] * w / vol * shape_n(k, xi);
                }
                held_load[fi][d] += l.total[d] * w / vol;
                add_axis_moment(ndof, &mut held_moment[fi], at_xi(ctr, size, xi), d, l.total[d] * w / vol);
            }
        }
    }
    applied.push((l.name.clone(), l.total));
    applied_moments.push(total_m);
    Ok(())
}
