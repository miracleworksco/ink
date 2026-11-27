//! Step 4, clamps: what each fix holds (`hold_fixes`): a mirror's plane, a
//! film's and a face's surface, a contact's normal, and a field's region
//! (`field_holds`), cell by cell.

use crate::cells::{K, lattice, place};
use crate::element::{full_moments, shape_n};
use crate::eval::{Ctx, eval_interval, eval_points};
use crate::hold;
use crate::occupancy::{
    Classified, Cov, Faces, at_xi, clipped_cov, occupied_centroid, occupied_parts, planes_part, region_faces, sets_part, spread_moments, surface_cell, surface_normals,
    tangent_offsets,
};
use crate::octree::Cell;
use crate::report::SurfaceMeasure;
use crate::surface::{Facet, facet_samples, lattice_surface, polyhedral_surface, split_cells, surface_measure};
use crate::tape::Tape;
use crate::{Boundary, RegionSel};
use ink_core::{V3, eval, interval, v3};
use rayon::prelude::*;

use super::settle::Cells;
use super::study::Study;
use super::{PENALTY, SolveError};

// 4. Clamps. A fix's region holds its material by Nitsche's method
// (`hold`): the held material leaves the stiffness, and its boundary with
// the free material carries the held value, weakly and consistently, so
// nothing slips and only the held material is held (a penalty over the
// held samples stiffened each cell the region's edge cut whole: a
// cantilever turned off the grid read 3.1% stiff on 0.5 mm cells). A
// fix on a named face holds the face, and a frictionless contact its
// surface along its normal. A symmetry plane keeps a penalty over its
// plane (PENALTY: a surface's, which stiffens no cell's material; its
// buckling sectors swap its held and free axes in the matrix, as a
// penalty's terms allow), and convection is its film's.
pub(super) struct Clamp {
    pub(super) cell: usize,
    pub(super) fix: usize,
    /// Its energy (`hold::Terms`): Nitsche's where the free material
    /// meets the held, a penalty's (a symmetry plane's; the spring that
    /// keeps held material's nodes from moving freely; what one fix
    /// holds of another's free motions), or a film's.
    pub(super) t: hold::Terms,
    /// The held material's volume (mm^3), and its share of each node's
    /// (int N_i over it) on each axis it holds (material two fixes
    /// hold, the first's): it holds no mass there.
    pub(super) vol: f64,
    pub(super) hn: [[f64; 8]; 3],
    /// A fix holding this one's free motions where they share material:
    /// what this pushes, that one pushes back.
    pub(super) other: Option<usize>,
    /// It holds all its cell's material (`HeldPart::whole`).
    pub(super) whole: bool,
    /// Of its terms, the spring (its constant, and int N_i N_j over the
    /// held material, on the fix's axes) that keeps held material's
    /// nodes from floating: no hold of the part (a body only it holds
    /// is free: `mechanisms`).
    pub(super) spring: Option<(f64, [f64; 64])>,
}

pub(super) fn clamp(cell: usize, fix: usize, t: hold::Terms) -> Clamp {
    Clamp { cell, fix, t: t.tie(), vol: 0.0, hn: [[0.0; 8]; 3], other: None, whole: false, spring: None }
}

// A fix's held material in one cell (int N_i N_j, int N_i v on the axes
// it holds, its volume and `Clamp::hn`), and the free material's
// boundary with it (`hold::Trace`s).
pub(super) struct HeldPart {
    pub(super) cell: usize,
    pub(super) fix: usize,
    pub(super) m: [f64; 64],
    pub(super) g: Vec<f64>,
    pub(super) vol: f64,
    pub(super) hn: [[f64; 8]; 3],
    pub(super) traces: Vec<hold::Trace>,
    /// Each trace's lender (`Lent`; usize::MAX: none), as `traces`.
    pub(super) lenders: Vec<usize>,
    /// Every sample of the cell's material wholly in its region: the
    /// cell holds no free material (decided by its samples: its
    /// moments less each sample's held part left 0.07% of a cell at a
    /// bore, whose nodes, held by nothing else, the bore's pressure
    /// threw 0.4 mm, and its buckling factor read 0.6).
    pub(super) whole: bool,
}

pub(super) fn held_part(ne: usize, cell: usize, fix: usize, traces: Vec<hold::Trace>) -> HeldPart {
    HeldPart { cell, fix, m: [0.0; 64], g: vec![0.0; ne], vol: 0.0, hn: [[0.0; 8]; 3], traces, lenders: Vec::new(), whole: false }
}

// A field's region holds its material. Cell by cell, each sample's
// material is split by the regions it lies in (each by its tangent
// plane there) and kept as the moments of each combination of them
// (`CellHolds`): what each fix holds, and what is free of each axis or
// of every region, come from the one split, so material two fixes hold
// counts once (each fix's held material kept the stiffness of the axes
// it left free, another's or not: a clamp held again along y read 1.5%
// stiff, and as three fixes of one axis each, 1.8%). A region's boundary
// through the free material is its trace, on the axes no other region
// holds there.
pub(super) struct CellHolds {
    /// The fixes (bit k of a combination is fixes[k]).
    pub(super) fixes: Vec<usize>,
    /// The moments (coefficient-weighted) of the material in each
    /// combination of their regions (0: in none).
    pub(super) pieces: Vec<(u32, [f64; 27])>,
    /// Its held material, by samples: combination, place (xi), volume, point.
    pub(super) held: Vec<(u32, V3, f64, V3)>,
    /// A laminate's (`laminate`, one face): each combination's volume
    /// and its material a's, in samples.
    pub(super) shares: Vec<(u32, f64, f64)>,
}

/// The part's surface, for the integrals over it (convection, pressure
/// loads, frictionless holds), the test of a symmetry plane and the
/// stress reads: meshed once, when one asks, at half the cells it lies
/// in (`part_surface`).
pub(super) struct Surfaces {
    /// The cells' grid cells split, and their depth (`split_cells`).
    pub(super) split: Vec<([i64; 3], u8)>,
    /// The pressure loads' surface, if exact (lattice or polyhedral), and its facets' points.
    pub(super) pressure: Option<(Vec<Facet>, Vec<V3>)>,
    pub(super) facets: Option<Vec<Facet>>,
    /// The facets' points on the surface (`surface_samples`), for regions.
    pub(super) samples: Option<Vec<V3>>,
    /// Each surface integral's region: what its facets covered against the
    /// lines' measure (`SurfaceMeasure`).
    pub(super) measures: Vec<SurfaceMeasure>,
}

impl Surfaces {
    pub(super) fn new(st: &Study, cells: &Cells) -> Self {
        let Study { s, ctx, dims, .. } = *st;
        let classified = &cells.classified;
        let split = split_cells(classified, dims);
        let pressure_surface = s.loads.iter().any(|l| l.pressure.is_some()).then(|| lattice_surface(s, ctx).or_else(|| polyhedral_surface(s, ctx))).flatten().map(|facets| {
            let points = facet_samples(&facets);
            (facets, points)
        });
        Surfaces { split, pressure: pressure_surface, facets: None, samples: None, measures: Vec::new() }
    }
}

/// What the fixes hold: their clamps (`Clamp`), their held material and
/// edges by cell (`HeldPart`, sorted by cell, then fix), and the cells'
/// tables of their regions' material (`CellHolds`).
pub(super) struct Holds {
    pub(super) clamps: Vec<Clamp>,
    /// Each film's point on the surface: (fix, cell, point).
    pub(super) film_points: Vec<(usize, usize, V3)>,
    pub(super) held_parts: Vec<HeldPart>,
    /// The fixes that hold their region's material.
    pub(super) field: Vec<usize>,
    /// Each fix's axes free motions share across a mirror (bit per axis).
    pub(super) shared_free: Vec<u8>,
    pub(super) tables: std::collections::HashMap<usize, CellHolds>,
    /// What the study says of its holds (a contact leaving a corner).
    pub(super) notes: Vec<String>,
}

/// Two fixes holding one material at different values.
fn contradiction(a: &str, b: &str, at: V3, x: f64, y: f64) -> String {
    format!(
        "fixes '{a}' and '{b}' hold the same material near {:?} mm on one axis at different values ({x} and {y}): they contradict each other; hold it by one, or give them the same value",
        at.map(|v| (v * 1000.0).round() / 1000.0)
    )
}

pub(super) fn hold_fixes(st: &Study, cells: &Cells, surf: &mut Surfaces) -> Result<Holds, SolveError> {
    let Study { s, ctx, ndof, ne, h, dims, ref grid, ref xis, ref distributing, .. } = *st;
    let classified = &cells.classified;
    let classified_at = |idx: [i64; 3]| cells.classified_at(st, idx);
    let classified_finest = classified.iter().map(|c| c.level).max().unwrap_or(0);
    let mut clamps: Vec<Clamp> = Vec::new();
    let mut notes = Vec::new();
    let mut film_points: Vec<(usize, usize, V3)> = Vec::new();
    let mut held_parts: Vec<HeldPart> = Vec::new();
    let mut field = Vec::new();
    let mut shared_free = vec![0u8; s.fixes.len()];
    for (fi, fx) in s.fixes.iter().enumerate() {
        if let Boundary::Mirror { axis, at } = fx.boundary {
            // Require material on only one side. The author declares the
            // other half to be its reflection, including all load vectors.
            let domain = std::array::from_fn(|a| interval::I::new(s.lo[a], s.hi[a]));
            let empty_side = |lower: bool| {
                let mut b = domain;
                if lower {
                    b[axis].hi = at.min(b[axis].hi)
                } else {
                    b[axis].lo = at.max(b[axis].lo)
                };
                b[axis].lo >= b[axis].hi || eval_interval(&s.shape, ctx, b, None).lo >= 0.0
            };
            let inward = match (empty_side(true), empty_side(false)) {
                (true, false) => 1.0,
                (false, true) => -1.0,
                _ => {
                    return Err(format!("mirror '{}' needs the entire part on one side of its cut at {at} mm. Cut the shape at the mirror or use a physical fix", fx.name).into());
                }
            };
            let found: Vec<Clamp> = classified
                .par_iter()
                .enumerate()
                .filter_map(|(ci, c)| {
                    let (ctr, size) = place(grid, c);
                    let t = 2.0 * (at - ctr[axis]) / size;
                    if !(-1.0..=1.0).contains(&t) || t * inward >= 1.0 {
                        return None;
                    }
                    let mut m = [0.0; 64];
                    for (k, &xi) in xis.iter().enumerate() {
                        let layer = (k / K.pow(axis as u32)) % K;
                        let nearest = (((t + 1.0) * 0.5 * K as f64) as usize).min(K - 1);
                        if layer != nearest || c.mask & (1 << k) == 0 {
                            continue;
                        }
                        let mut trace = xi;
                        trace[axis] = t;
                        let mut inside = at_xi(ctr, size, trace);
                        inside[axis] = at + inward * 1e-8 * size;
                        if eval::eval_point(&s.shape, ctx, inside) >= 0.0 {
                            continue;
                        }
                        let n: [f64; 8] = std::array::from_fn(|i| shape_n(i, trace));
                        for i in 0..8 {
                            for j in 0..8 {
                                m[i * 8 + j] += size.powi(3) * n[i] * n[j] / (K * K) as f64;
                            }
                        }
                    }
                    m.iter().any(|v| *v > 0.0).then(|| {
                        // As stiff in every cell it touches, however little
                        // of the cell it holds: per the held volume (m's
                        // sum). Per volume, a symmetry plane (- x 0.01) at 1
                        // mm cells held at a hundredth of a full cell's
                        // strength and slid: a quarter Lame tube grew 3% too far.
                        let vol: f64 = m.iter().sum();
                        Clamp { vol, ..clamp(ci, fi, hold::penalty(ndof, &fx.dofs, &m, &[0.0; 24][..ne], PENALTY * c.coef * size / vol)) }
                    })
                })
                .collect();
            if found.is_empty() {
                return Err(format!("mirror '{}' meets no material on its plane at {at} mm. Place it on the cut face or refine (res r)", fx.name).into());
            }
            // A single rigid support touching its reflected copy has one
            // amplitude. Establish contact on the cut before removing a
            // release in any sector. An interval overlap alone can include
            // a disjoint support, which has an independent reflected copy.
            for (i, support) in s.fixes.iter().enumerate().filter(|(_, f)| !f.free.is_empty()) {
                let mut possible = false;
                let contact = found.iter().any(|cl| {
                    let (ctr, size) = place(grid, &classified[cl.cell]);
                    let mut trace = std::array::from_fn(|a| interval::I::new(ctr[a] - 0.5 * size, ctr[a] + 0.5 * size));
                    trace[axis] = interval::I::point(at);
                    if !support.region.touches_plane(ctx, trace) {
                        return false;
                    }
                    possible = true;
                    // Include the edges: end plates and pins can occupy
                    // much less than a quadrature point's width there.
                    let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
                    let offsets = [-0.5 + 1e-8, -0.25, 0.0, 0.25, 0.5 - 1e-8];
                    let points: Vec<_> = offsets
                        .iter()
                        .flat_map(|du| {
                            offsets.iter().map(move |dv| {
                                let mut p = ctr;
                                p[axis] = at;
                                p[u] += du * size;
                                p[v] += dv * size;
                                p
                            })
                        })
                        .collect();
                    support.region.holds_on_surface(ctx, &points, 0.0).into_iter().zip(points).any(|(inside, mut p)| {
                        if !inside || eval::eval_point(&s.shape, ctx, p) > 0.0 {
                            return false;
                        }
                        p[axis] += inward * 1e-8 * size;
                        eval::eval_point(&s.shape, ctx, p) < 0.0
                    })
                });
                if !contact && possible {
                    return Err(format!(
                        "cannot establish whether the free support '{}' meets the mirror. Refine (res r), model the whole part, or mark a physical roller (support)",
                        support.name
                    )
                    .into());
                }
                if contact {
                    // Reflection preserves tangential slides and turns about
                    // the plane normal. The latter's pivot can lie anywhere
                    // on its axis. No component may move across the plane.
                    if support.free.iter().any(|m| m.v[axis] != 0.0 || (0..3).any(|d| d != axis && m.omega[d] != 0.0)) {
                        return Err(SolveError::InvalidMirror(format!(
                            "fix '{}' has a free motion on a support that may touch the mirror plane and changes under reflection. Model the whole part, or mark a physical roller (support)",
                            support.name
                        )));
                    }
                    shared_free[i] |= 1 << axis;
                }
            }
            clamps.extend(found);
            continue;
        }
        if fx.film.is_some() {
            // Convection, over the surface's facets, one point each (as a
            // pressure). A band of samples under the surface, weighted by
            // its depth, counted those it caught: the count aliased with the
            // grid (a 3 mm pin's perimeter read 15% short with its axis on a
            // cell centre, 4% on average), so benchmark 16's rise moved from
            // 39.4 to 43.1 K with where the grid sat.
            // Each facet's share of the region, at the share's centroid.
            let facets = surf.facets.get_or_insert_with(|| st.whole_surface(&surf.split));
            let samples = surf.samples.get_or_insert_with(|| st.sample_faces(facets));
            let w = fx.region.on_surface(ctx, facets, samples, h);
            let ambient = st.values_at(fx, &w.iter().map(|x| x.1).collect::<Vec<V3>>());
            let (mut at, covered): (Vec<(usize, V3, f64, f64, V3)>, Vec<(Facet, f64)>) = facets
                .par_iter()
                .zip(&w)
                .zip(&ambient)
                .filter(|((_, w), _)| w.0 > 0.0)
                .filter_map(|((f, &(w, q)), v)| {
                    let (ci, xi) = surface_cell(classified, &classified_at, grid, classified_finest, q, f.normal, &|_| true)?;
                    Some(((ci, xi, f.area * w, v[0], q), (*f, w)))
                })
                .unzip();
            surf.measures.push(SurfaceMeasure::of(&fx.name, "convection", &covered, surface_measure(&s.shape, ctx, grid, dims, classified, &surf.split, &fx.region, s.lo), grid));
            if at.is_empty() {
                return Err(format!(
                    "convection '{}' touches none of the part's surface (its region and the shape's surface do not meet): give it a region that \
                     holds the faces it cools",
                    fx.name
                )
                .into());
            }
            film_points.extend(at.iter().map(|&(ci, _, _, _, q)| (fi, ci, q)));
            at.sort_by_key(|a| a.0);
            clamps.extend(at.chunk_by(|a, b| a.0 == b.0).map(|on| {
                let (mut m, mut g) = ([0.0; 64], vec![0.0; ne]);
                for &(_, xi, w, v, _) in on {
                    let n: [f64; 8] = std::array::from_fn(|k| shape_n(k, xi));
                    for i in 0..8 {
                        g[i * ndof] += w * n[i] * v;
                        for j in 0..8 {
                            m[i * 8 + j] += w * n[i] * n[j];
                        }
                    }
                }
                Clamp { vol: m.iter().sum(), ..clamp(on[0].0, fi, hold::penalty(ndof, &fx.dofs, &m, &g, fx.film.unwrap_or(0.0))) }
            }));
            continue;
        }
        // A named face, or a frictionless contact: held on the part's
        // surface, each facet's share where it lies, a contact along the
        // surface's own normal (the part's gradient there, exact where its
        // field is: a bore's radial, so a turn about its axis stays free;
        // never the region's, so the same surface holds the same way
        // however the region is written), a face on the fix's axes.
        let proj: [[f64; 3]; 3] = std::array::from_fn(|a| std::array::from_fn(|b| if a == b && a < ndof && fx.dofs[a] { 1.0 } else { 0.0 }));
        let normal = fx.normal && ndof == 3;
        if normal || matches!(fx.region, RegionSel::Face { .. }) {
            let facets = surf.facets.get_or_insert_with(|| st.whole_surface(&surf.split));
            let samples = surf.samples.get_or_insert_with(|| st.sample_faces(facets));
            let w = fx.region.on_surface(ctx, facets, samples, h);
            let at: Vec<V3> = w.iter().map(|x| x.1).collect();
            let (values, normals) = (st.values_at(fx, &at), surface_normals(&s.shape, ctx, &at));
            // The axes a field's region holds where it holds the face's
            // material (one holds what lies in it: a face held again inside
            // a clamp whose edge cut its cells read 2.6% stiff).
            // (Strictly in it: a region touching the face from outside
            // holds none of its material. Where both hold an axis at
            // values that differ, neither moving freely, the two contradict
            // each other: refused, naming both, not settled for either.)
            // Per facet: the axes other fixes hold there, the value the first
            // holding each holds it at, and the first whose region lies on
            // the facet's face (its gradient nearer the face's normal than
            // its plane: a layer on that face, as `top`'s z > 19.99 is on a
            // column's top).
            let mut taken = vec![[false; 3]; at.len()];
            let mut held_at = vec![[0.0; 3]; at.len()];
            let mut face_of: Vec<Option<usize>> = vec![None; at.len()];
            let mut holder: Vec<Option<usize>> = vec![None; at.len()];
            // (a region with no gradient there, a constant or a kink's
            // crest, neither lies on a face nor crosses it)
            let leaning = |g: &[f64; 4], n: V3, on: bool| {
                let grad = [g[1], g[2], g[3]];
                let len = v3::norm(grad);
                len > 0.0 && (v3::dot(grad, n).abs() >= std::f64::consts::FRAC_1_SQRT_2 * len) == on
            };
            for (j, f) in s.fixes.iter().enumerate().filter(|(j, f)| *j != fi && !f.normal && !distributing[*j]) {
                let RegionSel::Field(rt) = &f.region else { continue };
                let mut v = Vec::new();
                eval_points(rt, ctx, &at, &mut v);
                let theirs = st.values_at(f, &at);
                let grads = if normal { crate::eval::eval_grads(rt, ctx, &at) } else { Vec::new() };
                for (k, _) in v.iter().enumerate().filter(|(k, x)| **x < 0.0 && w[*k].0 > 0.0) {
                    if !normal {
                        // (a face fix holds its axes: two values on one, neither
                        // moving freely, contradict each other)
                        for a in (0..ndof).filter(|&a| f.dofs[a] && fx.dofs[a]) {
                            let (x, y) = (values[k][a], theirs[k][a]);
                            if fx.free.is_empty() && f.free.is_empty() && (x - y).abs() > 1e-9 * x.abs().max(y.abs()).max(1e-12) {
                                return Err(contradiction(&fx.name, &f.name, at[k], x, y).into());
                            }
                        }
                    }
                    for a in 0..3 {
                        if f.dofs[a] && !taken[k][a] {
                            (taken[k][a], held_at[k][a]) = (true, theirs[k][a]);
                        }
                    }
                    // (one holding none of the normal's axes is not what a
                    // contact meets there)
                    if let Some(n) = normals[k].filter(|_| normal)
                        && (0..3).any(|a| f.dofs[a] && n[a].abs() > 1e-9)
                    {
                        holder[k].get_or_insert(j);
                        if face_of[k].is_none() && leaning(&grads[k], n, true) {
                            face_of[k] = Some(j);
                        }
                    }
                }
            }
            // A frictionless contact holds u . n at each facet. Where other
            // fixes hold some of n's axes, it holds what is left: u . m, m
            // the rest of n, at what their values leave of its own (n.v
            // less their axes' share), the direction m / |m| held at that
            // over |m|. Where they hold all of n, it holds nothing more,
            // and a value of theirs along n not its own contradicts it,
            // refused, but at a corner: the facet lies on a face another
            // fix's region lies on, which the contact's region crosses (its
            // gradient nearer the face's plane than its normal), and the
            // contact has a face of its own elsewhere, one no other fix
            // holds along its normal (a floor's strip at y < 0.01 takes in
            // the top face's edge, its normal z, where `top` moves it: the
            // floor holds its own face, y = 0). There it leaves the facet
            // to the other fix whole. None of it reads the cells: a facet
            // is a corner's however fine they are. A pad on a moved face,
            // no face of its own, or a lid on it, its region a layer on
            // that face, is refused.
            #[derive(Clone, Copy)]
            enum Normal {
                Hold(V3, V3),
                Leave,
                Corner,
            }
            let contact: Vec<Normal> = if normal {
                let own_grads = match &fx.region {
                    RegionSel::Field(rt) => Some(crate::eval::eval_grads(rt, ctx, &at)),
                    RegionSel::Face { .. } => None,
                };
                let share = |k: usize, n: V3| -> V3 { std::array::from_fn(|a| if taken[k][a] { n[a] } else { 0.0 }) };
                let own_face = (0..at.len()).any(|k| w[k].0 > 0.0 && normals[k].is_some_and(|n| v3::dot(share(k, n), share(k, n)) <= 1e-18));
                let mut out = Vec::with_capacity(at.len());
                for k in 0..at.len() {
                    let Some(n) = normals[k].filter(|_| w[k].0 > 0.0) else {
                        out.push(Normal::Leave);
                        continue;
                    };
                    let (na, rest) = (share(k, n), std::array::from_fn::<f64, 3, _>(|a| if taken[k][a] { 0.0 } else { n[a] }));
                    if v3::dot(na, na) <= 1e-18 {
                        out.push(Normal::Hold(n, values[k]));
                        continue;
                    }
                    let (mine, theirs) = (v3::dot(n, values[k]), v3::dot(na, held_at[k]));
                    let wanted = mine - theirs;
                    let differs = fx.free.is_empty() && holder[k].is_some_and(|j| s.fixes[j].free.is_empty()) && wanted.abs() > 1e-9 * mine.abs().max(theirs.abs()).max(1e-12);
                    let corner = own_face && face_of[k].is_some() && own_grads.as_ref().is_some_and(|g| leaning(&g[k], n, false));
                    let left = v3::norm(rest);
                    if left * left <= 1e-18 {
                        if differs && !corner {
                            let other = holder[k].map_or("", |j| s.fixes[j].name.as_str());
                            return Err(contradiction(&fx.name, other, at[k], mine, theirs).into());
                        }
                        out.push(if differs { Normal::Corner } else { Normal::Leave });
                        continue;
                    }
                    if corner {
                        out.push(if differs { Normal::Corner } else { Normal::Leave });
                        continue;
                    }
                    let m = rest.map(|v| v / left);
                    let along = wanted / left - v3::dot(m, values[k]);
                    out.push(Normal::Hold(m, std::array::from_fn(|a| values[k][a] + m[a] * along)));
                }
                out
            } else {
                vec![Normal::Leave; at.len()]
            };
            let corner: f64 = facets.iter().zip(&w).zip(&contact).filter(|(_, c)| matches!(c, Normal::Corner)).map(|((f, w), _)| f.area * w.0).sum();
            if corner > 0.0 {
                notes.push(format!(
                    "contact '{}' leaves {:.4} mm^2 of the faces it takes in at a corner to the fix whose face they are ({}), holding its own face there",
                    fx.name,
                    corner,
                    {
                        let mut names: Vec<&str> =
                            contact.iter().zip(&face_of).filter(|(c, _)| matches!(c, Normal::Corner)).filter_map(|(_, j)| j.map(|j| s.fixes[j].name.as_str())).collect();
                        names.sort_unstable();
                        names.dedup();
                        names.iter().map(|n| format!("'{n}'")).collect::<Vec<_>>().join(", ")
                    }
                ));
            }
            let (on, covered): (Vec<Option<(usize, hold::Trace)>>, Vec<(Facet, f64)>) = facets
                .par_iter()
                .zip(&w)
                .zip(&values)
                .zip(&normals)
                .zip(taken.par_iter().zip(&contact))
                .filter(|((((_, w), _), _), _)| w.0 > 0.0)
                .filter_map(|((((f, &(w, q)), v), n), (taken, contact))| {
                    let (n, mut p) = if normal { (n.map(|n| (n, std::array::from_fn(|a| std::array::from_fn(|b| n[a] * n[b])))))? } else { (n.unwrap_or(f.normal), proj) };
                    let mut v = *v;
                    if normal {
                        // (a facet it leaves to another is still of the
                        // surface its region takes in, as its measure reads;
                        // the trace's normal stays the surface's, which its
                        // terms read the stress across)
                        let Normal::Hold(m, held) = *contact else { return Some((None, (*f, w))) };
                        v = held;
                        p = std::array::from_fn(|a| std::array::from_fn(|b| m[a] * m[b]));
                    } else {
                        (0..3).filter(|&a| taken[a]).for_each(|a| (0..3).for_each(|b| (p[a][b], p[b][a]) = (0.0, 0.0)));
                        if p.iter().flatten().all(|v| *v == 0.0) {
                            return None;
                        }
                    }
                    let (ci, xi) = surface_cell(classified, &classified_at, grid, classified_finest, q, f.normal, &|_| true)?;
                    Some((Some((ci, hold::Trace { hold: 0, xi, w: f.area * w, wt: f.area * w, n, v, p })), (*f, w)))
                })
                .unzip();
            let mut on: Vec<(usize, hold::Trace)> = on.into_iter().flatten().collect();
            surf.measures.push(SurfaceMeasure::of(&fx.name, "fix", &covered, surface_measure(&s.shape, ctx, grid, dims, classified, &surf.split, &fx.region, s.lo), grid));
            // (none of it left as another fix's already: it holds nothing more)
            if on.is_empty() && !w.iter().zip(&taken).any(|(x, t)| x.0 > 0.0 && (0..ndof).any(|a| t[a])) {
                return Err(if normal {
                    format!(
                        "fix '{}' holds along the surface's normal, and its region meets none of the part's surface: give it a region that takes in the faces it holds",
                        fx.name
                    )
                } else {
                    format!("fix '{}' touches no material (its region and the shape do not overlap)", fx.name)
                }
                .into());
            }
            on.sort_by_key(|a| a.0);
            held_parts.extend(on.chunk_by(|a, b| a.0 == b.0).map(|on| held_part(ne, on[0].0, fi, on.iter().map(|x| x.1).collect())));
            continue;
        }
        field.push(fi);
    }
    let tables = field_holds(st, cells, &field, &mut held_parts)?;
    // Each cell holding held material, or the boundary of the free material
    // with it: its free material's stiffness (the cell's own where it holds
    // none), and its holds' terms (`hold::nitsche`) with that. A held
    // region's material keeps its free axes' stiffness, and on its held
    // ones, ALPHA's spring to the held value, as void keeps ALPHA of the
    // material: no node is left without stiffness.
    held_parts.sort_by_key(|p| (p.cell, p.fix));
    Ok(Holds { clamps, film_points, held_parts, field, shared_free, tables, notes })
}

impl Holds {
    /// A point strictly outside every rigid field fix holding a classified
    /// cell on axis `d` (scale-invariant via the region's gradient): used
    /// where a load or film on a dropped or absorbed free lip lands on a cell
    /// with no free material on that axis.
    pub(super) fn held_free_at(&self, st: &Study, ci_cls: usize, d: usize, at: V3) -> Option<usize> {
        let Study { s, ctx, h, ref distributing, .. } = *st;
        let held_parts = &self.held_parts;
        let on = &held_parts[held_parts.partition_point(|p| p.cell < ci_cls)..held_parts.partition_point(|p| p.cell <= ci_cls)];
        let mut holding: Option<usize> = None;
        for p in on {
            let fi = p.fix;
            let fx = &s.fixes[fi];
            if !fx.dofs[d] || distributing[fi] || fx.normal || fx.film.is_some() || matches!(fx.boundary, Boundary::Mirror { .. }) {
                continue;
            }
            // A face fix holds its surface alone. It neither owns this
            // cell's material nor changes which volume fixes hold it.
            let RegionSel::Field(rt) = &fx.region else { continue };
            let g = crate::eval::eval_grads(rt, ctx, &[at])[0];
            let scale = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt().max(1e-30);
            if g[0] <= 1e-9 * h * scale {
                return None;
            }
            holding.get_or_insert(fi);
        }
        holding
    }
}

/// Whether a point `q` of a region's face, its plane `g` (held where at
/// most 0), is on the part's face or past it, facing as it does: the
/// region's field `r` and the part's field and gradient `m` (material
/// where at most 0) there, each as a distance (by its plane's slope and
/// its own), the part's no farther out than the region's but by floats'
/// rounding at q's size, and their normals within 60 degrees (on the
/// part's edge its field's gradient is its two faces' blend; the region's
/// plane, not its field's gradient, which on its own edge is its faces'
/// blend too). Free material beyond the region's face, however thin, is
/// a lip and its edge: within a millionth of a sample, a lip 1e-10 mm thick
/// pressed on its face went untraced, its pressure into the fix's reaction,
/// and the study settled.
fn on_face(q: V3, g: &[f64; 4], r: f64, m: &[f64; 4], sub: f64) -> bool {
    let slope = |p: &[f64; 4]| (p[1] * p[1] + p[2] * p[2] + p[3] * p[3]).sqrt();
    let (sg, sm) = (slope(g), slope(m));
    if !(sg > 0.0 && sm > 0.0 && sg.is_finite() && sm.is_finite()) {
        return false;
    }
    let rounding = 64.0 * f64::EPSILON * (v3::norm(q) + sub);
    r / sg - m[0] / sm <= rounding && (g[1] * m[1] + g[2] * m[2] + g[3] * m[3]) / (sg * sm) > 0.5
}

/// A cell's material (moments, coefficient-weighted).
pub(super) fn material_of(c: &Classified) -> [f64; 27] {
    let fm = full_moments();
    c.moments.unwrap_or(fm.map(|v| v * c.coef))
}

/// Adds `w` times moments `m` to a combination's (`mask`) in `to`.
pub(super) fn add(to: &mut Vec<(u32, [f64; 27])>, mask: u32, m: [f64; 27], w: f64) {
    match to.iter_mut().find(|p| p.0 == mask) {
        Some(p) => p.1.iter_mut().zip(m).for_each(|(a, b)| *a += w * b),
        None => to.push((mask, m.map(|v| w * v))),
    }
}

/// The field fixes' held material and edges, cell by cell (added to
/// `held_parts`), and each cell's table of its regions' material.
fn field_holds(st: &Study, cells: &Cells, field: &[usize], held_parts: &mut Vec<HeldPart>) -> Result<std::collections::HashMap<usize, CellHolds>, SolveError> {
    let Study { s, ctx, ndof, ne, ref grid, ref xis, ref distributing, .. } = *st;
    let Cells { ref classified, depth, .. } = *cells;
    #[allow(clippy::type_complexity)]
    let found: Vec<(usize, CellHolds, Vec<HeldPart>, Vec<usize>, Option<(usize, usize, V3, f64, f64)>)> = classified
        .par_iter()
        .enumerate()
        .filter_map(|(ci, c)| {
            let cell_box = grid.cell_box(Cell { i: c.idx, size: 1 });
            // (to a rounding past its box: a region's face on the cell's
            // plane, a rounding short of it, is this cell's samples' to trace
            // where the cell beside it cannot)
            let reach: Vec<usize> =
                field.iter().copied().filter(|&fi| matches!(&s.fixes[fi].region, RegionSel::Field(rt) if eval_interval(rt, ctx, cell_box, None).lo <= 1e-9 * grid.h)).collect();
            if reach.is_empty() {
                return None;
            }
            let (ctr, size) = place(grid, c);
            let (sub, to_xi) = (size / K as f64, 2.0 / size);
            let pts: Vec<V3> = xis.iter().map(|&xi| at_xi(ctr, size, xi)).collect();
            // Each sample's material (the cell's, as its moments took it) and
            // its tangent plane.
            let tape = c.tape.as_deref().unwrap_or(&s.shape);
            let whole = sub * sub / 12.0;
            let (mat, mg): (Vec<(f64, V3, Cov)>, Vec<[f64; 4]>) = if c.frac < 1.0 {
                (occupied_parts(tape, ctx, &pts, size), tangent_offsets(tape, ctx, &pts, size).1)
            } else {
                (vec![(1.0, [0.0; 3], [whole, whole, whole, 0.0, 0.0, 0.0]); pts.len()], vec![[-1.0, 0.0, 0.0, 0.0]; pts.len()])
            };
            // (A laminate's moments are its mean's: its stiffness is its own.)
            let co = if c.moments.is_some() && c.lam.is_none() { st.coef_at(&pts) } else { vec![c.coef; pts.len()] };
            // Each region's share of each sample, its centroid and tangent
            // plane, as they are: no share is rounded to all or none (a
            // sample left 3% free once counted all held, and its edge
            // traced elsewhere left seams: `hold::plane_piece`).
            // (A sample holding several of a region's faces carries each:
            // `region_faces`; in a laminate's cell, one.)
            let regs: Vec<(Vec<f64>, Vec<V3>, Vec<Faces>)> = reach
                .iter()
                .map(|&fi| {
                    let RegionSel::Field(rt) = &s.fixes[fi].region else { unreachable!("a field's region") };
                    let (share, faces) = region_faces(rt, ctx, &pts, size, c.lam.is_none(), false);
                    let shares = share.iter().map(|x| x.0.clamp(0.0, 1.0)).collect();
                    (shares, share.iter().map(|x| x.1).collect(), faces)
                })
                .collect();
            // (Each combination's moments summed over the samples, the free
            // part's too: taken as the cell's less the samples the regions
            // reach, a free part of 0.5% of a cell sampled otherwise than
            // the cell had moments of no material, indefinite.)
            let mut pieces = vec![(0u32, [0.0; 27])];
            let (mut held, mut bound, mut shares) = (Vec::new(), vec![Vec::new(); reach.len()], Vec::<(u32, f64, f64)>::new());
            let mut free_any = false;
            // A laminate's face, at a sample (as a region's plane): the free
            // material's share of each material is not the cell's, and its
            // stiffness is its own laminate's (held 0.05 mm off a face of k
            // 0.2 and 200 within a 0.2 mm cell, a bar read 4.4 K hot).
            let face = c.lam.as_deref().and_then(|l| Some((l.n, l.d?)));
            let mut share = |k: usize, mask: u32, v: f64, planes: &[[f64; 4]], whole_sample: bool| {
                let Some((n, d)) = face else { return };
                let g = [v3::dot(n, xis[k]) - d, n[0] * to_xi, n[1] * to_xi, n[2] * to_xi];
                let first = (!whole_sample).then_some(&mg[k]);
                let all = planes_part(first, planes, sub).0;
                let a = if g[0].abs() > 0.5 * sub * 3f64.sqrt() * v3::norm([g[1], g[2], g[3]]) || all <= 0.0 {
                    (g[0] < 0.0) as u8 as f64
                } else if first.is_none() && planes.is_empty() {
                    crate::query::occupancy(&g, [sub; 3])
                } else {
                    planes_part(first, &[planes, &[g][..]].concat(), sub).0 / all
                };
                match shares.iter_mut().find(|x| x.0 == mask) {
                    Some(x) => (x.1, x.2) = (x.1 + v, x.2 + v * a),
                    None => shares.push((mask, v, v * a)),
                }
            };
            for k in (0..pts.len()).filter(|&k| c.mask >> k & 1 == 1 && mat[k].0 > 0.0) {
                let (sv, so, scov) = mat[k];
                // (A sample within 1e-9 of whole is whole, rounded: a face on
                // the cell's plane left its samples 1 - 1e-15, and split by
                // points, a band 0.4 of a sample deep read 0.375: the bar it
                // held read 0.15 K low.)
                let whole_sample = sv > 1.0 - 1e-9;
                let inside = (0..reach.len()).filter(|&j| regs[j].0[k] >= 1.0).fold(0u32, |m, j| m | 1 << j);
                let crossing: Vec<usize> = (0..reach.len()).filter(|&j| regs[j].0[k] > 0.0 && regs[j].0[k] < 1.0).collect();
                let xi_of = |d: V3| -> V3 { std::array::from_fn(|a| xis[k][a] + d[a] * to_xi) };
                let moments = |d: V3, cov: &Cov| spread_moments(xi_of(d), &cov.map(|v| v * to_xi * to_xi));
                if inside == 0 && crossing.is_empty() {
                    free_any = true;
                    share(k, 0, sv, &[], whole_sample);
                    add(&mut pieces, 0, moments(so, &scov), 8.0 / pts.len() as f64 * co[k] * sv);
                } else {
                    // Its parts, each in its combination of regions: across one
                    // plane, exactly; across more (with the material's), by
                    // points, scaled to its material. Regions crossing it on
                    // one plane split it as one (by points, one clamp written
                    // as three fixes of one axis each read 0.045% soft).
                    let mut across: Vec<(u32, usize)> = Vec::new();
                    for &j in &crossing {
                        let unit = |g: [f64; 4]| {
                            let l = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
                            g.map(|v| v / l)
                        };
                        let one = |i: usize| (regs[i].2[k].planes.len() == 1).then(|| unit(regs[i].2[k].planes[0]));
                        let gj = one(j);
                        match across.iter_mut().find(|(_, i)| gj.zip(one(*i)).is_some_and(|(gj, gi)| gi.iter().zip(&gj).all(|(a, b)| (a - b).abs() <= 1e-12 * (1.0 + sub)))) {
                            Some(a) => a.0 |= 1 << j,
                            None => across.push((1 << j, j)),
                        }
                    }
                    let parts: Vec<(u32, f64, V3, Cov, Vec<[f64; 4]>)> = match across[..] {
                        [] => vec![(inside, sv, so, scov, Vec::new())],
                        [(m, j)] if whole_sample && regs[j].2[k].planes.len() == 1 => {
                            let (r, g) = (regs[j].0[k], regs[j].2[k].planes[0]);
                            let away = g.map(|v| -v);
                            vec![
                                (inside | m, r, regs[j].1[k], clipped_cov(&g, sub), vec![g]),
                                (inside, 1.0 - r, occupied_centroid(&away, sub), clipped_cov(&away, sub), vec![away]),
                            ]
                        }
                        _ => {
                            let mut out: Vec<(u32, f64, V3, Cov, Vec<[f64; 4]>)> = (0..1u32 << across.len())
                                .map(|t| {
                                    // (each region's part or the rest; its planes for a laminate's share)
                                    let first = (!whole_sample).then_some(&mg[k]);
                                    let single = across.iter().all(|&(_, j)| regs[j].2[k].planes.len() == 1);
                                    let planes: Vec<[f64; 4]> = if single {
                                        across.iter().enumerate().map(|(q, &(_, j))| regs[j].2[k].planes[0].map(|v| if t >> q & 1 == 1 { v } else { -v })).collect()
                                    } else {
                                        Vec::new()
                                    };
                                    let (v, d, cov) = sets_part(first, &across.iter().enumerate().map(|(q, &(_, j))| regs[j].2[k].set(t >> q & 1 == 1)).collect::<Vec<_>>(), sub);
                                    (across.iter().enumerate().filter(|(q, _)| t >> q & 1 == 1).fold(inside, |m, (_, &(a, _))| m | a), v, d, cov, planes)
                                })
                                .collect();
                            let total: f64 = out.iter().map(|p| p.1).sum();
                            if total > 0.0 {
                                out.iter_mut().for_each(|p| p.1 *= sv / total);
                            }
                            out
                        }
                    };
                    for (mask, v, d, cov, planes) in parts.into_iter().filter(|p| p.1 > 0.0) {
                        add(&mut pieces, mask, moments(d, &cov), 8.0 / pts.len() as f64 * co[k] * v);
                        share(k, mask, v, &planes, whole_sample);
                        if mask == 0 {
                            free_any = true;
                        } else {
                            held.push((mask, xi_of(d), v * sub.powi(3), v3::add(pts[k], d)));
                        }
                    }
                }
                // Each region's boundary through the sample's material, its
                // piece in the sample's cube (on a face two samples share, the
                // free side's: `hold::plane_piece`), on the axes no other region
                // holds there: one holds what lies in it, and of a boundary
                // two share, the first (else that axis was held twice).
                for (j, pi) in (0..reach.len()).flat_map(|j| (0..regs[j].2[k].planes.len()).map(move |pi| (j, pi))) {
                    // (each of its faces there, within the others)
                    let fj = &regs[j].2[k];
                    let g = fj.planes[pi];
                    let Some(n) = v3::direction([-g[1], -g[2], -g[3]]) else { continue };
                    let mut clips: Vec<[f64; 4]> = if whole_sample { Vec::new() } else { vec![mg[k]] };
                    clips.extend(fj.planes.iter().enumerate().filter(|(q, _)| *q != pi).map(|(_, o)| if fj.any { o.map(|v| -v) } else { *o }));
                    // (a plane on the sample's face: the sample on the held side
                    // traces it where the region's field at the face's centre is
                    // free, the one on the free side where it is held or 0; the
                    // centre from its place on the lattice, so both samples sharing
                    // the face evaluate one number. And where the sample so chosen
                    // does not carry this plane, the other does: a sample has one
                    // tangent plane, and at a region's corner it is the other
                    // face's. A block's side face on a sample's face at y = 0.45,
                    // which rounding gave to the held side, lost a strip a sample
                    // wide along the block's end face, and one along the edge its
                    // two side faces meet at: 11.701 mm^2 traced of 11.8825.)
                    let owns = |axis: usize, side: f64| {
                        let (cl, cs) = lattice(c, depth);
                        let at = [k % K, (k / K) % K, k / (K * K)];
                        let unit = grid.h / ((2 * K as i64) << depth) as f64;
                        let face: V3 = std::array::from_fn(|a| {
                            let l = cl[a] * 2 * K as i64 + cs * (2 * at[a] as i64 + 1) + if a == axis { side as i64 * cs } else { 0 };
                            grid.origin[a] + unit * l as f64
                        });
                        let RegionSel::Field(rt) = &s.fixes[reach[j]].region else { unreachable!("a field's region") };
                        let at_face = eval::eval_point(rt, ctx, face);
                        if (at_face > 0.0) == (g[0] < 0.0) {
                            return true;
                        }
                        // (only where the face is the region's edge: a plane
                        // carried on past the region's other faces is no edge)
                        if !(at_face.abs() <= 1e-6 * sub * (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt()) {
                            return false;
                        }
                        // The other sample's plane, as it reads it itself: in this
                        // cell, its own; in the cell beside, at its sample's centre
                        // placed as that cell places it (which the region must reach).
                        let mut o = at.map(|v| v as i64);
                        o[axis] += side as i64;
                        let other: Faces = if (0..K as i64).contains(&o[axis]) {
                            regs[j].2[(o[0] + K as i64 * (o[1] + K as i64 * o[2])) as usize].clone()
                        } else {
                            let mut low = cl;
                            low[axis] += side as i64 * cs;
                            o[axis] = o[axis].rem_euclid(K as i64);
                            let idx: [i64; 3] = std::array::from_fn(|a| low[a] >> depth);
                            if eval_interval(rt, ctx, grid.cell_box(Cell { i: idx, size: 1 }), None).lo > 1e-9 * grid.h {
                                return true;
                            }
                            let lo = grid.pos(idx);
                            let ctr: V3 = std::array::from_fn(|a| lo[a] + size * (((low[a] - (idx[a] << depth)) / cs) as f64 + 0.5));
                            region_faces(rt, ctx, &[at_xi(ctr, size, xis[(o[0] + K as i64 * (o[1] + K as i64 * o[2])) as usize])], size, c.lam.is_none(), false).1.remove(0)
                        };
                        let unit_of = |p: &[f64; 4]| p.map(|v| v / (p[1] * p[1] + p[2] * p[2] + p[3] * p[3]).sqrt());
                        let me = unit_of(&g);
                        let same = |ot: [f64; 4]| (1..4).all(|a| (me[a] - ot[a]).abs() <= 1e-6) && (ot[0] - (me[0] + me[1 + axis] * side * sub)).abs() <= 1e-6 * sub;
                        !other.planes.iter().any(|o| same(unit_of(o)))
                    };
                    let piece = hold::polygon_points(&hold::plane_piece(&g, sub, &clips, &owns));
                    // (Where the region's face is the part's own, facing as it
                    // does, beyond it is air, no free material: no edge between
                    // held and free. Traced, a plate 0.2 mm thick on 0.5 mm
                    // cells held over x < 1 and z < 0.2, its top face, took 891
                    // iterations for 104, read 0.04% stiff, and said the hold
                    // held free slivers along the face; a rod held under its
                    // own round face read 0.7% stiff. Read at the piece's
                    // points by both fields, so a curved face is one too.)
                    let on_part: Vec<bool> = if whole_sample || piece.is_empty() {
                        vec![false; piece.len()]
                    } else {
                        let RegionSel::Field(rt) = &s.fixes[reach[j]].region else { unreachable!("a field's region") };
                        let at: Vec<V3> = piece.iter().map(|p| v3::add(pts[k], p.0)).collect();
                        let (region, part) = (crate::eval::eval_grads(rt, ctx, &at), crate::eval::eval_grads(tape, ctx, &at));
                        region.iter().zip(&part).zip(&at).map(|((r, m), &q)| on_face(q, &g, r[0], m, sub)).collect()
                    };
                    for ((d, w), _) in piece.into_iter().zip(on_part).filter(|x| !x.1) {
                        let mut axes = s.fixes[reach[j]].dofs;
                        for i in (0..reach.len()).filter(|&i| i != j && !distributing[reach[i]]) {
                            let at = regs[i].2[k].at(d);
                            if regs[i].0[k] >= 1.0 || at < -1e-9 * sub || (at.abs() <= 1e-9 * sub && i < j) {
                                (0..3).for_each(|a| axes[a] &= !s.fixes[reach[i]].dofs[a]);
                            }
                        }
                        if (0..ndof).any(|a| axes[a]) {
                            bound[j].push((xi_of(d), w, n, v3::add(pts[k], d), axes));
                        }
                    }
                }
            }
            if held.is_empty() && bound.iter().all(|b| b.is_empty()) {
                return None;
            }
            // (a distributing hold's material is free material: `distributing`)
            let rigid = (0..reach.len()).filter(|&j| !distributing[reach[j]]).fold(0u32, |m, j| m | 1 << j);
            if pieces.iter().any(|p| p.0 & !rigid != 0) {
                free_any |= pieces.iter().any(|p| p.0 != 0 && p.0 & rigid == 0 && p.1[0] > 0.0);
                let mut merged = vec![(0u32, [0.0; 27])];
                for (mask, m) in pieces {
                    add(&mut merged, mask & rigid, m, 1.0);
                }
                pieces = merged;
                let mut kept: Vec<(u32, f64, f64)> = Vec::new();
                for (mask, v, a) in shares {
                    match kept.iter_mut().find(|x| x.0 == mask & rigid) {
                        Some(x) => (x.1, x.2) = (x.1 + v, x.2 + a),
                        None => kept.push((mask & rigid, v, a)),
                    }
                }
                shares = kept;
            }
            if !free_any {
                pieces[0].1 = [0.0; 27];
            }
            let mut parts = Vec::new();
            for (j, &fi) in reach.iter().enumerate() {
                let fx = &s.fixes[fi];
                let mine: Vec<&(u32, V3, f64, V3)> = held.iter().filter(|x| x.0 >> j & 1 == 1).collect();
                if mine.is_empty() && bound[j].is_empty() {
                    continue;
                }
                let values = st.values_at(fx, &mine.iter().map(|x| x.3).chain(bound[j].iter().map(|x| x.3)).collect::<Vec<V3>>());
                let mut p = held_part(ne, ci, fi, Vec::new());
                // (a distributing hold's material is free: its int N_i N_j
                // alone, its weight's (`freeing`))
                for (&&(mask, xi, dv, _), v) in mine.iter().zip(&values) {
                    let n: [f64; 8] = std::array::from_fn(|i| shape_n(i, xi));
                    if !distributing[fi] {
                        p.vol += dv;
                    }
                    for i in 0..8 {
                        for d in (0..ndof).filter(|&d| fx.dofs[d] && !distributing[fi]) {
                            p.g[i * ndof + d] += dv * n[i] * v[d];
                            if (0..j).all(|b| mask >> b & 1 == 0 || !s.fixes[reach[b]].dofs[d]) {
                                p.hn[d][i] += dv * n[i];
                            }
                        }
                        for jj in 0..8 {
                            p.m[i * 8 + jj] += dv * n[i] * n[jj];
                        }
                    }
                }
                p.whole = !distributing[fi] && (0..pts.len()).filter(|&k| c.mask >> k & 1 == 1 && mat[k].0 > 0.0).all(|k| regs[j].0[k] >= 1.0);
                p.traces = bound[j]
                    .iter()
                    .zip(&values[mine.len()..])
                    .map(|(b, v)| hold::Trace {
                        hold: 0,
                        xi: b.0,
                        w: b.1,
                        wt: b.1,
                        n: b.2,
                        v: *v,
                        p: std::array::from_fn(|a| std::array::from_fn(|q| if a == q && a < ndof && b.4[a] { 1.0 } else { 0.0 })),
                    })
                    .collect();
                parts.push(p);
            }
            // The fixes whose region and the cell's material overlap: both
            // fields at most 0 over a box, by their intervals, the cell's box
            // halved, depth first, until one field leaves it or it is a hair
            // across (sqrt(epsilon) of a sample: far under any length the
            // cells resolve, far over a coordinate's rounding). (Points half a
            // sample about each held piece fell either side of a hold 1e-6 mm
            // thick at an end face, and it was refused as touching nothing;
            // intervals bound the fields, so a region in a cell's air 0.05 mm
            // off its material is still refused.)
            fn overlap(rt: &Tape, shape: &Tape, ctx: &Ctx, lo: V3, hi: V3, finest: f64) -> bool {
                let b: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(lo[a], hi[a]));
                let (r, m) = (eval_interval(rt, ctx, b, None), eval_interval(shape, ctx, b, None));
                if r.lo > 0.0 || m.lo > 0.0 {
                    return false;
                }
                if (r.hi <= 0.0 && m.hi <= 0.0) || hi[0] - lo[0] <= finest {
                    return true;
                }
                let mid: V3 = std::array::from_fn(|a| 0.5 * (lo[a] + hi[a]));
                (0..8).any(|q| {
                    let (l, h): (V3, V3) =
                        (std::array::from_fn(|a| if q >> a & 1 == 0 { lo[a] } else { mid[a] }), std::array::from_fn(|a| if q >> a & 1 == 0 { mid[a] } else { hi[a] }));
                    overlap(rt, shape, ctx, l, h, finest)
                })
            }
            let (box_lo, box_hi): (V3, V3) = (std::array::from_fn(|a| ctr[a] - 0.5 * size), std::array::from_fn(|a| ctr[a] + 0.5 * size));
            let touching: Vec<usize> = reach
                .iter()
                .copied()
                .filter(|&fi| matches!(&s.fixes[fi].region, RegionSel::Field(rt) if overlap(rt, tape, ctx, box_lo, box_hi, sub * f64::EPSILON.sqrt())))
                .collect();
            // Two fixes holding one piece on one axis at different values,
            // neither moving it freely: they contradict each other.
            let clash = held.iter().filter(|h| h.2 > 0.0 && h.0.count_ones() > 1).find_map(|h| {
                let on: Vec<usize> = (0..reach.len()).filter(|&j| h.0 >> j & 1 == 1).map(|j| reach[j]).collect();
                let vals: Vec<[f64; 3]> = on.iter().map(|&f| st.values_at(&s.fixes[f], &[h.3])[0]).collect();
                (0..on.len()).flat_map(|a| (a + 1..on.len()).map(move |b| (a, b))).find_map(|(a, b)| {
                    let (fa, fb) = (&s.fixes[on[a]], &s.fixes[on[b]]);
                    if !(fa.free.is_empty() && fb.free.is_empty()) {
                        return None;
                    }
                    (0..ndof).filter(|&d| fa.dofs[d] && fb.dofs[d]).find_map(|d| {
                        let (x, y) = (vals[a][d], vals[b][d]);
                        ((x - y).abs() > 1e-9 * x.abs().max(y.abs()).max(1e-12)).then_some((on[a], on[b], h.3, x, y))
                    })
                })
            });
            Some((ci, CellHolds { fixes: reach, pieces, held, shares }, parts, touching, clash))
        })
        .collect();
    let mut tables: std::collections::HashMap<usize, CellHolds> = std::collections::HashMap::new();
    let mut touched = std::collections::HashSet::new();
    for (ci, table, parts, touching, clash) in found {
        if let Some((a, b, at, x, y)) = clash {
            return Err(format!(
                "fixes '{}' and '{}' hold the same material near {:?} mm on one axis at different values ({x} and {y}): they contradict each other; hold it by one, or give them the same value",
                s.fixes[a].name,
                s.fixes[b].name,
                at.map(|v| (v * 1000.0).round() / 1000.0)
            )
            .into());
        }
        held_parts.extend(parts);
        tables.insert(ci, table);
        touched.extend(touching);
    }
    if let Some(&fi) = field.iter().find(|&&fi| !touched.contains(&fi)) {
        return Err(format!("fix '{}' touches no material (its region and the shape do not overlap)", s.fixes[fi].name).into());
    }
    Ok(tables)
}
