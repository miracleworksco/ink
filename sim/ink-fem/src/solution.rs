//! A solve's answer and what is read off it: displacements, stresses and
//! heat, reactions, peaks, a surface's radius.

use crate::cells::{ALPHA, SimCell, locate, nearest_piece, place};
use crate::element::{b_matrix, d_matrix, element_basis, full_moments, laminate_mean, laminate_strains, shape_dn, shape_n, small_turn};
use crate::eval::{Ctx, eval_points};
use crate::occupancy::{at_xi, surface_cell};
use crate::octree::Grid;
use crate::report::{Report, spread_of};
use crate::tape::Tape;
use crate::{GRAIN, Physics};
use ink_core::{V3, grid, interval, v3};
use ink_kernel::octree;
use rayon::prelude::*;
use std::sync::Arc;

/// Default distance (mm) from fix regions within which stresses are reported
/// separately from the headline peak; never less than 2 cells. A fixed physical
/// distance keeps the headline peak comparable across resolutions.
pub const SUPPORT_MARGIN_MM: f64 = 2.0;

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Solution {
    pub physics: Physics,
    pub h: f64,
    /// The box analysed: the study shape's bounds.
    pub bounds: (V3, V3),
    pub(crate) origin: V3,
    pub(crate) dims: [i64; 3],
    pub(crate) cells: Vec<SimCell>,
    pub(crate) cell_at: Vec<u32>,
    /// Each node's place on the lattice of nodes, whose step is the
    /// cell over 2^depth (`MAX_DEPTH`).
    pub node_idx: Vec<[i64; 3]>,
    #[serde(default)]
    pub(crate) depth: u8,
    /// Whether each node is reliable for reporting: inside material and backed
    /// by at least a quarter of its eight surrounding cells' worth of material.
    /// Nodes on slivers of cut cells are held weakly and can exaggerate motion.
    pub(crate) node_reliable: Vec<bool>,
    /// Solution per dof: displacement (mm) or temperature.
    pub u: Vec<f64>,
    /// The nodes inside held material, sorted, and what is read at them
    /// (`nodal`): the held value where held, `u` where free.
    #[serde(default)]
    pub(crate) held: Vec<(u32, [f64; 3])>,
    /// Per cell at its centre: stress (Voigt, MPa) or heat flux (W/mm^2) in the first 3.
    pub(crate) cell_field: Vec<[f64; 6]>,
    /// Per cell, the first elastic fix (its index) within `support_margin`
    /// of it, NO_FIX where none is. Rigid clamps create singular stresses at
    /// their edges (Saint-Venant): real supports spread load.
    pub(crate) near_support: Vec<u16>,
    /// Per fix, a complete face held only in its normal direction, with
    /// no interior support edge to omit from the stress peak.
    #[serde(default)]
    pub plane_holds: Vec<bool>,
    /// Points of the part's surface (at the facets of its mesh at half a
    /// cell), each in the cell at least half material that holds it or
    /// lies under it: the cell and the point in its coordinates, sorted by
    /// cell. Where stress peaks are read (`peaks`).
    #[serde(default)]
    pub(crate) surface_points: Vec<(u32, V3)>,
    /// Each of those points in the cell it lies in (where a cell holds it),
    /// as `surface_points`: where a peak value is read (`max_nodal`).
    #[serde(default)]
    pub(crate) surface_own: Vec<(u32, V3)>,
    /// Each of those points' facet's normal, out of the material: by which
    /// a point that may lie by an inside corner is told (`crease`).
    #[serde(default)]
    pub(crate) surface_normals: Vec<V3>,
    /// Each of those points in a fix's region: the axes it holds (bits) and
    /// its held value there, which a peak reads in place of the free
    /// field's extension (`max_nodal`).
    #[serde(default)]
    pub(crate) surface_held: Vec<Option<(u8, [f64; 3])>>,
    /// Points of the surface between those (`climb_reads`), each in the
    /// cell whose stress rises to it, as `surface_points`, sorted by cell:
    /// not read, but how far the peak may lie over its reads.
    #[serde(default)]
    pub(crate) climbed: Vec<(u32, V3)>,
    /// `peaks`, found when first asked.
    #[serde(skip)]
    pub(crate) peaks: std::sync::OnceLock<Vec<(f64, V3)>>,
    /// `reads`, found when first asked.
    #[serde(skip)]
    pub(crate) reads: std::sync::OnceLock<Vec<(u32, f64, V3)>>,
    pub report: Report,
    /// Modal analysis: every mode's shape (per dof, 1 mm peak).
    pub(crate) mode_shapes: Vec<Vec<f64>>,
}

/// The least the buckling factor of a weighted sum of solves' loads may be, from each solve's own
/// analysis and its loads' reversed (`Report::buckle`, `buckle_reversed`). One over the factor is the
/// largest eigenvalue of a matrix linear in the stress, which is linear in the loads: a convex function
/// of them that scales with them, so at the sum it is at most the sum of each term's, its weight's
/// size times its own (a negative weight's, its reversed loads'). Each term is taken at its band's low
/// end. Err: a term with a weight has no analysis that way (not run, or its eigensolve gave no factor),
/// and why.
pub fn buckle_of_sum(terms: &[(f64, &Solution)]) -> Result<f64, String> {
    let mut over = 0.0;
    for &(w, t) in terms.iter().filter(|(w, _)| w.abs() > 1e-12) {
        let b = if w > 0.0 { t.report.buckle.as_ref() } else { t.report.buckle_reversed.as_ref() };
        let Some(b) = b else { return Err("a run it is summed from has no buckling analysis of its loads that way (a part cut at a mirror has none of them reversed)".into()) };
        let Some((low, _)) = b.answer() else { return Err(b.failure.clone().unwrap_or_else(|| "no factor resolved".into())) };
        over += w.abs() / low;
    }
    Ok(1.0 / over)
}

pub fn von_mises(s: &[f64; 6]) -> f64 {
    (0.5 * ((s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2)) + 3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5])).sqrt()
}

/// A scalar summary per cell: von Mises stress (MPa) or heat flux magnitude.
fn cell_scalar(p: Physics, v: &[f64; 6]) -> f64 {
    match p {
        Physics::Elastic { .. } => von_mises(v),
        Physics::Heat => v3::norm([v[0], v[1], v[2]]),
    }
}

/// A solution on its regular grid (see `Solution::grid_fields`). Arrays are
/// x-fastest. Cell `(i, j, k)` spans nodes `i..=i+1` (and so on); node
/// `(i, j, k)` sits at `origin + h * (i, j, k)`.
pub struct GridFields {
    pub origin: V3,
    pub h: f64,
    /// Cells per axis; nodes are `dims + 1` per axis.
    pub dims: [usize; 3],
    /// Per cell: von Mises stress (MPa) or heat flux magnitude; NaN without material.
    pub cell_scalar: Vec<f32>,
    /// Per cell: material volume fraction (0 without material).
    pub cell_frac: Vec<f32>,
    /// Per node: displacement (3 values, mm) or temperature (1); NaN where there is no node.
    pub node_values: Vec<f32>,
}

#[derive(Clone, Copy, Debug)]
pub struct Extreme {
    pub value: f64,
    pub at: V3,
}

/// The tightest concave radius of the surface of `shape` within `reach` of
/// `p`, curving into the material (a hole, a fillet, a notch), where stress
/// concentrates, and where it is: infinite where the surface curves nowhere
/// so (a flat, a convex edge); 0 at a crease (a sharp inside corner, or a
/// round under half of `step`, which no finer cells resolve). None: no
/// surface within reach (a peak inside the material, at an interface or
/// across a section).
///
/// Read by the field's sign alone, so any field with the same zero set
/// gives the same (CORE: a shape is where its field is negative). The
/// surface is surveyed on a lattice a quarter of `reach` apart about `p`:
/// where two neighbouring points differ in sign, the surface between them
/// (by bisection). At each such point, along directions across its normal,
/// the surface's height off its tangent plane `s` each way, (h+ + h-) /
/// s^2, its bend that way; the largest is its concave curvature. The 4
/// points bending most over the lattice's spacing (so a crease within reach
/// is straddled) are read in full: a crease where the surface's normal
/// turns, across it, by an angle that holds as the span shrinks (by over
/// 0.001 rad, at 1.5 and 0.75 spacings: a curve's turn halves); else its
/// bend at `step` on 16 directions. A round under about the span (0.4 of
/// a cell of `reach / 1.5`) reads as a crease too: unresolved either way.
pub fn surface_radius(shape: &Arc<Tape>, ctx: &Ctx, p: V3, reach: f64, step: f64) -> Option<(V3, f64)> {
    surface_radius_as(shape, ctx, p, reach, step, false)
}

/// `surface_radius`, and (`measured`) with a crease only where its turn was
/// measured: a span that runs off an edge beside the point, or onto a point
/// with no slope, has no turn, and is not taken for one of 0 (a curved
/// flank within a span of a blade's tip edge then read a crease: its short
/// span's turn against none); the next two shorter spans are compared, and
/// where neither pair is on the surface, no crease is said. So of a whole
/// surface's survey (`Solution::crease`), where a crease must be shown; a
/// peak's, which a doubt leaves unknown, takes the other.
fn surface_radius_as(shape: &Arc<Tape>, ctx: &Ctx, p: V3, reach: f64, step: f64, measured: bool) -> Option<(V3, f64)> {
    // The shape simplified over the box the survey reads: its lattice, and
    // heights up to 8 spans either side of points 1.5 lattice spacings off
    // (a spacing a quarter of the reach).
    let half = 5.0 * reach + 9.0 * step;
    let b: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(p[a] - half, p[a] + half));
    let local = octree::Local::new(shape, ctx, p.map(|v| v - half), 2.0 * half);
    let shape = &*local.tape_over(&b);
    // Material at a point: the field under 0 there, or within a millionth
    // of a step of it by its own slope. Two bodies meeting on a face are 0
    // all over it, or off 0 by rounding, and are solved as one: read by the
    // sign alone, the join was a slit, and where it met a flat face, a
    // crease there.
    let inside = |q: &[V3]| -> Vec<bool> { crate::eval::eval_grads(shape, ctx, q).into_iter().map(|g| g[0] <= 1e-6 * step * v3::norm([g[1], g[2], g[3]])).collect() };
    // The crossing on the segment a..b (their signs differ), by bisection.
    let cross = |a: V3, b: V3, ia: bool| -> V3 {
        let (mut a, mut b) = (a, b);
        for _ in 0..24 {
            let m = v3::scale(v3::add(a, b), 0.5);
            if inside(&[m])[0] == ia { a = m } else { b = m }
        }
        v3::scale(v3::add(a, b), 0.5)
    };
    // The lattice about p, and the surface between its neighbours.
    let n = 4i64;
    let d = reach / n as f64;
    let ix = |i: i64, j: i64, k: i64| (((i + n) * (2 * n + 1) + (j + n)) * (2 * n + 1) + (k + n)) as usize;
    let lattice: Vec<V3> = (-n..=n).flat_map(|i| (-n..=n).flat_map(move |j| (-n..=n).map(move |k| [p[0] + d * i as f64, p[1] + d * j as f64, p[2] + d * k as f64]))).collect();
    let sign = inside(&lattice);
    let mut found: Vec<V3> = Vec::new();
    for i in -n..=n {
        for j in -n..=n {
            for k in -n..=n {
                for o in [[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
                    let (a, b) = (ix(i, j, k), [i + o[0], j + o[1], k + o[2]]);
                    if b.iter().any(|&v| v > n) {
                        continue;
                    }
                    let b = ix(b[0], b[1], b[2]);
                    if sign[a] != sign[b] {
                        let q = cross(lattice[a], lattice[b], sign[a]);
                        if v3::norm(v3::sub(q, p)) <= reach && found.iter().all(|f| v3::norm(v3::sub(*f, q)) > 0.75 * d) {
                            found.push(q);
                        }
                    }
                }
            }
        }
    }
    // Its normal out of the material: the field's gradient's direction,
    // or, where the field is flat on its surface, a hair off it.
    let normal = |q: V3| -> Option<V3> {
        let g = crate::eval::eval_grads(shape, ctx, &[q])[0];
        v3::direction([g[1], g[2], g[3]]).or_else(|| {
            let o: Vec<V3> = (0..3).flat_map(|a| [1.0, -1.0].map(|sg| std::array::from_fn(|b| q[b] + if a == b { sg * 1e-6 * step } else { 0.0 }))).collect();
            let gs = crate::eval::eval_grads(shape, ctx, &o);
            v3::direction(gs.iter().fold([0.0; 3], |s, g| v3::add(s, v3::direction([g[1], g[2], g[3]]).unwrap_or([0.0; 3]))))
        })
    };
    // The surface's height off q along n at q + s u: the crossing nearest
    // the tangent plane within 8 s each way. None found: past that, up
    // (a wall rising: the point stepped to is in the material) or down (an
    // edge falling away).
    let height = |q: V3, nrm: V3, u: V3, s: f64| -> f64 {
        let base = v3::add(q, v3::scale(u, s));
        let m = 16;
        let ts: Vec<f64> = (-m..=m).map(|k| 8.0 * s * k as f64 / m as f64).collect();
        let pts: Vec<V3> = ts.iter().map(|&t| v3::add(base, v3::scale(nrm, t))).collect();
        let sg = inside(&pts);
        let mid = m as usize;
        let found = (0..m as usize).flat_map(|k| [mid + k, mid - k - 1]).find(|&k| sg[k] != sg[k + 1]);
        match found {
            Some(k) => v3::dot(v3::sub(cross(pts[k], pts[k + 1], sg[k]), q), nrm),
            None if sg[mid] => 8.0 * s,
            None => -8.0 * s,
        }
    };
    // The concave curvature at q on steps s: the largest bend, (h+ + h-) /
    // s^2, over 16 directions across its normal (concave positive: the
    // surface bending toward its outside normal; a smooth surface's within
    // 1% of its largest principal curvature; a convex edge's bend, negative
    // across it, does not count).
    let concave = |q: V3, nrm: V3, s: f64, dirs: usize| -> (f64, V3) {
        let a = if nrm[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let t1 = v3::unit(v3::cross(nrm, a));
        let t2 = v3::cross(nrm, t1);
        (0..dirs)
            .map(|k| {
                let th = std::f64::consts::PI * k as f64 / dirs as f64;
                let u = v3::add(v3::scale(t1, th.cos()), v3::scale(t2, th.sin()));
                ((height(q, nrm, u, s) + height(q, nrm, v3::scale(u, -1.0), s)) / (s * s), u)
            })
            .fold((0.0, t1), |m, b| if b.0 > m.0 { b } else { m })
    };
    // The turn of the surface's normal from q - r u to q + r u (on the
    // surface, as `height` finds it): a crease's holds as r shrinks, a
    // curve's shrinks with it. None: a wall rises past the heights read (a
    // notch too sharp for them), a crease.
    let turn = |q: V3, nrm: V3, u: V3, r: f64| -> Option<Option<f64>> {
        let side = |sg: f64| -> Option<Option<V3>> {
            let h = height(q, nrm, v3::scale(u, sg), r);
            if h == 8.0 * r {
                return None;
            }
            Some((h != -8.0 * r).then(|| v3::add(v3::add(q, v3::scale(u, sg * r)), v3::scale(nrm, h))).and_then(normal))
        };
        Some(side(1.0)?.zip(side(-1.0)?).map(|(a, b)| v3::dot(a, b).clamp(-1.0, 1.0).acos()))
    };
    // A normal turning by less than this is no crease (0.06 degrees).
    const KINK: f64 = 1e-3;
    let radius = |q: V3, nrm: V3, u: V3| -> f64 {
        // Of two spans' turns: a crease, or not; None where one has none.
        let held = |far: f64, near: f64| -> Option<bool> {
            match (turn(q, nrm, u, far * d), turn(q, nrm, u, near * d)) {
                (Some(far), Some(near)) if measured => far.zip(near).map(|(far, near)| near > KINK && near > 0.75 * far),
                (Some(far), Some(near)) => Some(near.unwrap_or(0.0) > KINK && near.unwrap_or(0.0) > 0.75 * far.unwrap_or(0.0)),
                _ => Some(true),
            }
        };
        let crease = held(1.5, 0.75).or_else(|| held(0.75, 0.375)).unwrap_or(false);
        if crease {
            return 0.0;
        }
        let (b, _) = concave(q, nrm, step, 16);
        // A bend under a millionth of a step's is none (the bisection's
        // rounding on a flat).
        if b <= 1e-6 / step { f64::INFINITY } else { 1.0 / b }
    };
    // Each point's bend over the lattice's spacing on 8 directions (a
    // crease within reach lies within that of a point found), then the 4
    // bending most in full.
    let mut coarse: Vec<(V3, Option<V3>, f64, V3)> = found
        .into_iter()
        .map(|q| (q, normal(q)))
        .map(|(q, n)| match n {
            Some(n) => {
                let (b, u) = concave(q, n, d, 8);
                (q, Some(n), b, u)
            }
            None => (q, None, f64::INFINITY, [1.0, 0.0, 0.0]),
        })
        .collect();
    coarse.sort_by(|a, b| b.2.total_cmp(&a.2));
    let unsloped = if measured { f64::INFINITY } else { 0.0 };
    coarse.into_iter().take(4).map(|(q, n, b, u)| (q, n.map_or(unsloped, |n| if b <= 1e-6 / d { f64::INFINITY } else { radius(q, n, u) }))).min_by(|a, b| a.1.total_cmp(&b.1))
}

/// Which cells, given as boxes, are not of one material (their E, nu or
/// yield not one value across the box): for each, None where it is, else
/// the least yield in it (NaN where unbounded). The cells a solve with one
/// stiffness per cell cannot read (`Solution::max_cell_ratio`).
pub type Mixed<'a> = &'a dyn Fn(&[(V3, V3)]) -> Vec<Option<MixedCell>>;

/// A cell not of one material (`Mixed`), over the region it is read over.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedCell {
    /// The least the scale (a yield strength) can be there.
    pub least: f64,
    /// The most its stiffness can be there: E / (1 + nu), what sets von
    /// Mises stress for a strain (heat: k). NaN or infinite where unbounded.
    pub stiffest: f64,
}

/// A cell no fix's support margin reaches (Solution::near_support).
pub(crate) const NO_FIX: u16 = u16::MAX;

/// The support margin nearest a point: the fix (its index in the setup's
/// fixes), and the distance from the point to the centre of the nearest
/// cell the margin leaves out. `beside`: that cell is a neighbour of the
/// point's own (within a cell's diagonal), the point at the margin's edge.
#[derive(Clone, Copy, Debug)]
pub struct NearestSupport {
    pub fix: usize,
    pub distance: f64,
    pub beside: bool,
}

impl Solution {
    fn grid(&self) -> Grid {
        Grid::new(self.origin, self.h)
    }

    pub fn node_pos(&self, n: usize) -> V3 {
        let step = self.h / (1i64 << self.depth) as f64;
        std::array::from_fn(|a| self.origin[a] + step * self.node_idx[n][a] as f64)
    }

    /// Makes this the solution of the weighted sum of `terms`' loads: their
    /// fields (displacement or temperature, stress or flux) summed with
    /// those weights, which is that load's solution where the solve is
    /// linear in its loads (a static one) and `terms` are solves of one
    /// setup that differ in their loads alone, as this is. Its reactions
    /// and applied totals (forces and moments, or heat) are summed too; the rest of its report (the
    /// solver's, a buckling factor) is still its own solve's, not the
    /// sum's. False, and nothing changed, where a term is on other cells.
    pub fn set_sum(&mut self, terms: &[(f64, &Solution)]) -> bool {
        if terms.iter().any(|(_, t)| {
            t.h != self.h || t.u.len() != self.u.len() || t.cell_field.len() != self.cell_field.len() || t.node_idx.len() != self.node_idx.len() || t.depth != self.depth
        }) {
            return false;
        }
        self.u.par_iter_mut().with_min_len(GRAIN).enumerate().for_each(|(i, u)| *u = terms.iter().map(|(w, t)| w * t.u[i]).sum());
        self.cell_field.par_iter_mut().with_min_len(GRAIN).enumerate().for_each(|(i, f)| *f = std::array::from_fn(|k| terms.iter().map(|(w, t)| w * t.cell_field[i][k]).sum()));
        for (i, e) in self.report.applied.iter_mut().enumerate() {
            e.1 = std::array::from_fn(|k| terms.iter().map(|(w, t)| w * t.report.applied.get(i).map_or(0.0, |a| a.1[k])).sum());
        }
        for (i, e) in self.report.reactions.iter_mut().enumerate() {
            e.1 = std::array::from_fn(|k| terms.iter().map(|(w, t)| w * t.report.reactions.get(i).map_or(0.0, |a| a.1[k])).sum());
        }
        for (i, m) in self.report.applied_moments.iter_mut().enumerate() {
            *m = std::array::from_fn(|k| terms.iter().map(|(w, t)| w * t.report.applied_moments.get(i).map_or(0.0, |a| a[k])).sum());
        }
        for (i, m) in self.report.reaction_moments.iter_mut().enumerate() {
            *m = std::array::from_fn(|k| terms.iter().map(|(w, t)| w * t.report.reaction_moments.get(i).map_or(0.0, |a| a[k])).sum());
        }
        // What is read in held material, summed as the rest (a sum read
        // the first run's there); where the runs hold other nodes, the sum.
        if terms.iter().all(|(_, t)| t.held.len() == self.held.len() && t.held.iter().zip(&self.held).all(|(a, b)| a.0 == b.0)) {
            for (i, e) in self.held.iter_mut().enumerate() {
                e.1 = std::array::from_fn(|k| terms.iter().map(|(w, t)| w * t.held[i].1[k]).sum());
            }
        } else {
            self.held.clear();
        }
        (self.peaks, self.reads) = Default::default();
        true
    }

    /// Ok if the solve gave a solution to read, else why not
    /// (`Report::failure`). Every reduction below says so first.
    pub fn valid(&self) -> Result<(), String> {
        match &self.report.failure {
            Some(why) => Err(why.clone()),
            None => Ok(()),
        }
    }

    /// Largest displacement magnitude (elastic) or temperature (heat) over
    /// reliable nodes (all nodes if none are, e.g. parts thinner than h).
    /// What a node reads (displacement, or temperature first): the
    /// solution, or, inside held material, the held value (`held`).
    pub fn nodal(&self, n: usize) -> [f64; 3] {
        let nd = self.physics.ndof();
        match self.held.binary_search_by_key(&(n as u32), |h| h.0) {
            Ok(k) => self.held[k].1,
            Err(_) => std::array::from_fn(|d| if d < nd { self.u[n * nd + d] } else { 0.0 }),
        }
    }

    /// What cell `ci`'s field reads at its nodes (per dof): the held
    /// values where all of them are held, else the solution (a cell the
    /// held region's edge cuts is its free material's field).
    fn cell_values(&self, ci: usize) -> Vec<f64> {
        let (nd, c) = (self.physics.ndof(), &self.cells[ci]);
        let held = !self.held.is_empty() && c.nodes.iter().all(|n| self.held.binary_search_by_key(n, |h| h.0).is_ok());
        (0..8 * nd).map(|j| if held { self.nodal(c.nodes[j / nd] as usize)[j % nd] } else { self.u[c.nodes[j / nd] as usize * nd + j % nd] }).collect()
    }

    /// What a reading `x` of this solve (a displacement or a temperature at
    /// a node) is known to, relative to it. The stiffness kept in void and
    /// held material (ALPHA of the material's, over every cell: so no node
    /// is left without stiffness) stiffens the solve by up to ALPHA of the
    /// cells' share of the part, and the solve stops at its residual: of
    /// what the solve produces, the field's span over its nodes (its least
    /// to its most), not of an offset it carries whole (a part moved 10 mm
    /// with 0.01 of it strain, or held at 300 C and warmed 2: of it all, a
    /// 1e-6 floor passed 1e-6 of the offset as noise, and a reading moving
    /// apart under it was settled). A bar pulled to 7.1875 mm read
    /// 7.18749875, 7.18749935 and 7.1875000 on cells of 2, 1 and 0.5 mm:
    /// 9e-8 of it apart, within this (its span the whole pull); a patch's
    /// 1.022e-4, 3.262e-4 and 7.736e-4 mm, and a pin held at 300 C read
    /// 0.07% apart, stay moving apart. (In a wall bent across cells
    /// thinner than it, the void stiffens more than its share: about ALPHA
    /// (h/t)^3; under-stated, more grids are seen moving apart, not fewer.)
    pub fn precision_at(&self, x: f64) -> f64 {
        let g = self.grid();
        let cells: f64 = self.cells.iter().map(|c| place(&g, c).1.powi(3)).sum();
        (ALPHA * cells / self.report.volume.max(1e-300) + self.report.residual.max(0.0)) * self.span(x) / x.abs().max(1e-300)
    }

    /// What the field (a displacement's size, a temperature) varies by over
    /// the part's reliable nodes: what the solve's own error is a share of,
    /// not an offset it carries whole (a held move, a temperature's 273 K).
    /// `x`'s size where it has no nodes.
    pub fn span(&self, x: f64) -> f64 {
        let nd = self.physics.ndof();
        let any = self.node_reliable.iter().any(|&b| b);
        let (lo, hi) = (0..self.node_idx.len()).filter(|&n| !any || self.node_reliable[n]).fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), n| {
            let u = self.nodal(n);
            let v = if nd == 3 { v3::norm(u) } else { u[0] };
            (a.min(v), b.max(v))
        });
        if hi >= lo { hi - lo } else { x.abs() }
    }

    pub fn max_nodal(&self) -> Result<Extreme, String> {
        self.valid()?;
        if self.node_idx.is_empty() {
            return Err("the solution has no nodes".into());
        }
        let nd = self.physics.ndof();
        let any = self.node_reliable.iter().any(|&b| b);
        let size = |v: &[f64]| if nd == 3 { v.iter().map(|x| x * x).sum::<f64>().sqrt() } else { v[0] };
        let (n, v) = (0..self.node_idx.len())
            .filter(|&n| !any || self.node_reliable[n])
            .map(|n| {
                let u = self.nodal(n);
                let v = if nd == 3 { v3::norm(u) } else { u[0] };
                (n, v)
            })
            .fold((0, f64::NEG_INFINITY), |b, x| if x.1 > b.1 { x } else { b });
        // And at the points of the part's surface, through each one's cell's
        // field, as a peak stress is read (`reads`): a cut cell's extreme is
        // on the part's face. Read at the nodes inside alone, a bar whose end
        // lay between nodes read its last node inside, up to a cell short,
        // and the grids short alike: stretched exactly 0.2 mm, 0.19255 +-
        // 0.00448 on 0.3 mm cells.
        // (In the cell the point lies in, not as a stress read from the
        // cell under it: its field's extension carried a bar's slope past
        // where its heat curved it flat, its end 0.7 K high.)
        let at_surface = self
            .surface_own
            .iter()
            .enumerate()
            .map(|(k, &(ci, xi))| {
                let ci = ci as usize;
                let xi = xi.map(|v| v.clamp(-2.0, 2.0));
                // (a face in a fix's region reads its held value on the axes
                // it holds, not the free field's extension: a cube held at
                // 100 C read 101.29 there)
                let mut v = self.interpolate(ci, xi);
                if let Some(Some((mask, held))) = self.surface_held.get(k) {
                    (0..nd).filter(|&d| mask >> d & 1 == 1).for_each(|d| v[d] = held[d]);
                }
                (size(&v), at_xi(self.cell_center(ci), self.cell_size(ci), xi))
            })
            .filter(|x| x.0.is_finite())
            .fold((f64::NEG_INFINITY, [0.0; 3]), |b, x| if x.0 > b.0 { x } else { b });
        if at_surface.0 > v {
            return Ok(Extreme { value: at_surface.0, at: at_surface.1 });
        }
        Ok(Extreme { value: v, at: self.node_pos(n) })
    }

    /// Cells used for stress and heat-flux statistics: at least half material (cut cells with
    /// slivers of material produce unreliable stresses in the finite cell
    /// method) and, for elasticity, away from supports.
    fn reliable(&self) -> impl Iterator<Item = (usize, f64)> + '_ {
        self.cells.iter().enumerate().filter(|(i, c)| c.frac >= 0.5 && self.near_support[*i] == NO_FIX).map(|(i, _)| (i, cell_scalar(self.physics, &self.cell_field[i])))
    }

    /// Peak over reliable cells away from supports, at the surface where
    /// it peaks there (`peaks`). None reliable (a part thinner than two
    /// cells, or within the support margin all over) is no measurement: an
    /// error saying why.
    pub fn max_cell_scalar(&self) -> Result<Extreme, String> {
        let peaks = self.peaks();
        let (v, at) = self.reliable_values()?.into_iter().map(|(i, _)| peaks[i]).fold((f64::NEG_INFINITY, [0.0; 3]), |b, x| if x.0 > b.0 { x } else { b });
        Ok(Extreme { value: v, at })
    }

    /// The largest of the value over `scale` at every point read in a
    /// reliable cell away from the supports (`reads`), each over the scale
    /// at its own point: over a yield strength, one over the safety factor,
    /// so each material's own sets it where it is. Read only at each cell's
    /// peak, a weak strip elsewhere in the cell was missed: a 10 MPa strip
    /// across a bar at 10 MPa held a safety factor `> 2` at 24. A cell not of
    /// one material (`mixed`) gives no value: its element takes one
    /// stiffness, so a read in a softer or stiffer part of it is not that
    /// part's (a soft, weak skin read the bulk's stress against its own
    /// yield, and failed a part that holds); it is bounded only
    /// (`ratio_between_reads`).
    pub fn max_cell_ratio(&self, scale: &dyn Fn(&[V3]) -> Vec<f64>, mixed: Mixed) -> Result<Extreme, String> {
        let cells: Vec<usize> = self.reliable_values()?.into_iter().map(|c| c.0).collect();
        let one = self.of_one_material(&cells, mixed);
        if one.is_empty() {
            return Err(format!(
                "no cell of {} mm at least half material away from the supports is of one material: the factor is only bounded there (yield-between-reads); refine res",
                self.h
            ));
        }
        Ok(ratio_peak(&self.reads_in(one.into_iter()), scale).unwrap_or(Extreme { value: f64::NEG_INFINITY, at: [0.0; 3] }))
    }

    /// `max_near_supports`, over `scale` (as `max_cell_ratio`).
    pub fn max_near_supports_ratio(&self, scale: &dyn Fn(&[V3]) -> Vec<f64>, mixed: Mixed) -> Option<Extreme> {
        let near: Vec<usize> = self.cells.iter().enumerate().filter(|(i, c)| c.frac >= 0.5 && self.near_support[*i] != NO_FIX).map(|(i, _)| i).collect();
        ratio_peak(&self.reads_in(self.of_one_material(&near, mixed).into_iter()), scale)
    }

    /// Each cell's box (a split cell's own, not its grid cell's).
    fn cell_boxes(&self, cells: &[usize]) -> Vec<(V3, V3)> {
        cells.iter().map(|&ci| self.read_region(ci)).collect()
    }

    /// The region cell `ci` is read over: its box, and every point it
    /// reads (`reads`), up to half a cell past it (a face in a cell under
    /// half material is read from the cell under it). Whether it is of one
    /// material, and its bound, are taken over this: over its box alone, a
    /// soft skin just past it, where its face read lands, went unseen, and
    /// the bulk's stress was read against the skin's yield.
    fn read_region(&self, ci: usize) -> (V3, V3) {
        // Its own box (a split cell's, not its grid cell's).
        let (c, size) = place(&self.grid(), &self.cells[ci]);
        let (mut lo, mut hi): (V3, V3) = (c.map(|v| v - 0.5 * size), c.map(|v| v + 0.5 * size));
        for (_, p) in self.reads_in(std::iter::once(ci)) {
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
        (lo, hi)
    }

    /// The cells given that are of one material (`mixed` says None).
    fn of_one_material(&self, cells: &[usize], mixed: Mixed) -> Vec<usize> {
        cells.iter().zip(mixed(&self.cell_boxes(cells))).filter(|(_, m)| m.is_none()).map(|(&ci, _)| ci).collect()
    }

    /// The reads (`reads`) in the cells given, as (value, point).
    fn reads_in(&self, cells: impl Iterator<Item = usize>) -> Vec<(f64, V3)> {
        let reads = self.reads();
        cells.flat_map(|ci| reads[reads.partition_point(|r| (r.0 as usize) < ci)..].iter().take_while(move |r| r.0 as usize == ci).map(|r| (r.1, r.2))).collect()
    }

    /// The cells not of one material (`mixed`: a strip, a gradient, an
    /// interface; they give no value, `max_cell_ratio`) bounded: at least
    /// half material, away from the supports (`near` false) or within their
    /// margins, the largest the value over the scale can be anywhere in the
    /// cell (`cell_ratio_upper`; a least scale not above 0, or NaN, bounds
    /// nothing). Reads alone miss a strip
    /// between them, and one over a read but not the stress's peak (a strip
    /// on a beam's neutral axis, read there, its faces carrying the
    /// bending). The largest such bound and the cell's centre; None where
    /// every cell is of one material. As the cells shrink it closes on the
    /// ratio: the bound's reach is a cell.
    pub fn ratio_between_reads(&self, near: bool, mixed: Mixed) -> Option<Extreme> {
        self.valid().ok()?;
        let cells: Vec<usize> = (0..self.cells.len()).filter(|&i| self.cells[i].frac >= 0.5 && (self.near_support[i] != NO_FIX) == near).collect();
        cells
            .iter()
            .zip(mixed(&self.cell_boxes(&cells)))
            // A least scale not above 0 (NaN, or an interval's overshoot
            // below 0: a smoothstep's over a cell) bounds nothing: the ratio
            // unbounded, the factor 0, so the cell leaves a claim unknown.
            // Divided by, it gave a negative ratio no maximum picked, and the
            // cell, giving no value either, was never read.
            .filter_map(|(&ci, m)| m.map(|m| (self.cell_ratio_upper(ci, m, mixed), self.cell_center(ci))))
            .filter(|x| !x.0.is_nan())
            .max_by(|a, b| a.0.total_cmp(&b.0))
            .map(|(value, at)| Extreme { value, at })
    }

    /// The largest the cell's scalar (von Mises stress, flux) over the
    /// scale (a yield) can be anywhere it is read (`read_region`), for a
    /// cell not of one material (`m`). Each strain component (the
    /// incompatible modes included) is affine along every line parallel to
    /// an axis and von Mises or a magnitude is convex in them, so over a box
    /// they peak at its corners: the element's field there, and its reads.
    /// A laminate's (`laminate`) is each material's own stress, at every
    /// corner whichever side of its face, over its own yield where each has
    /// one, on its own side of the face over the region read (a weaker or
    /// stiffer material there taken: its least yield, at its stiffest). Any other cell's element takes one stiffness and
    /// its strain is the cell's: a point of it is strained as the cell is
    /// where its materials lie side by side along the load, and stressed as
    /// the cell is where they lie across it (the softer strained more), so
    /// its stress is bounded by its stiffest material's (`MixedCell::
    /// stiffest`) at the cell's strain. Taken at the cell's own stiffness, a
    /// strip stiffer than the cell's mean along the load read under its
    /// stress.
    fn cell_ratio_upper(&self, ci: usize, m: MixedCell, mixed: Mixed) -> f64 {
        let (lo, hi) = self.read_region(ci);
        let ctr = self.cell_center(ci);
        let corners: Vec<V3> = (0..8).map(|k| std::array::from_fn(|a| ((if k >> a & 1 == 1 { hi[a] } else { lo[a] }) - ctr[a]) * 2.0 / self.cell_size(ci))).collect();
        // (A NaN bounds nothing: the bound is unbounded, not the others'
        // largest, as `max` would skip it.)
        let most = |vs: &mut dyn Iterator<Item = f64>| vs.fold(f64::NEG_INFINITY, |m, v| if v.is_nan() { f64::INFINITY } else { m.max(v) });
        let over = |v: f64, scale: f64| if scale > 0.0 { v / scale } else { f64::INFINITY };
        let c = &self.cells[ci];
        if let Some(lam) = &c.lam {
            let Some(fs) = self.laminate_fields(ci, &corners, true) else { return f64::INFINITY };
            // Each material's least yield and stiffest over its own side of
            // the face, over the region read (a cell under half material
            // beside it is read from this one: a skin there weaker or stiffer
            // than either; or a yield's step between the samples): its side
            // of one material, its own (the samples' yield); else the fields'
            // bounds over the side, short of the face by a millionth of the
            // cell (the face's own material the other's). A face not along
            // an axis, or a layer within the cell, takes the whole region's
            // for both.
            let own = |(e, nu): (f64, f64)| match self.physics {
                Physics::Elastic { .. } => e / (1.0 + nu),
                Physics::Heat => e,
            };
            let (oa, ob) = (own(lam.a), own(lam.b));
            let whole = |o: f64| (m.least, if m.stiffest.is_nan() { f64::INFINITY } else { (m.stiffest / o).max(1.0) });
            let axis = (0..3).find(|&k| (lam.n[k].abs() - 1.0).abs() < 1e-12);
            let [(ya, ka), (yb, kb)] = match (axis, lam.d) {
                (Some(k), Some(d)) => {
                    let size = self.cell_size(ci);
                    let at = ctr[k] + lam.n[k] * d * 0.5 * size;
                    let gap = 1e-6 * size;
                    // (a below the face, along n)
                    let (mut sa, mut sb) = ((lo, hi), (lo, hi));
                    if lam.n[k] > 0.0 {
                        sa.1[k] = sa.1[k].min(at - gap);
                        sb.0[k] = sb.0[k].max(at + gap);
                    } else {
                        sa.0[k] = sa.0[k].max(at + gap);
                        sb.1[k] = sb.1[k].min(at - gap);
                    }
                    let sides = mixed(&[sa, sb]);
                    let [ya0, yb0] = lam.yields.unwrap_or([m.least; 2]);
                    let side = |r: &Option<MixedCell>, y: f64, o: f64| match r {
                        None => (y, 1.0),
                        Some(r) => (r.least, if r.stiffest.is_nan() { f64::INFINITY } else { (r.stiffest / o).max(1.0) }),
                    };
                    [side(&sides[0], ya0, oa), side(&sides[1], yb0, ob)]
                }
                _ => [whole(oa), whole(ob)],
            };
            let s = |f: &Option<[f64; 6]>| f.map_or(f64::NAN, |f| cell_scalar(self.physics, &f));
            return most(&mut fs.iter().flat_map(|[a, b]| [over(s(a) * ka, ya), over(s(b) * kb, yb)]));
        }
        let at = self.fields_in_cell(ci, &corners).map_or(f64::INFINITY, |fs| most(&mut fs.iter().map(|f| cell_scalar(self.physics, f)))).max(self.peaks()[ci].0);
        // The cell's own: E / (1 + nu), or k.
        let own = match self.physics {
            Physics::Elastic { .. } => c.coef / (1.0 + c.nu),
            Physics::Heat => c.coef,
        };
        let stiffer = if m.stiffest.is_nan() { f64::INFINITY } else { (m.stiffest / own).max(1.0) };
        over(at * stiffer, m.least)
    }

    /// The reliable cells' values, or why there are none.
    fn reliable_values(&self) -> Result<Vec<(usize, f64)>, String> {
        self.valid()?;
        let v: Vec<(usize, f64)> = self.reliable().collect();
        if !v.is_empty() {
            Ok(v)
        } else if self.cells.iter().any(|c| c.frac >= 0.5) {
            Err("no cell at least half material lies beyond the support margin of the fixes (support-margin, default 2 mm and at least \
                 2 cells): nothing away from the supports to measure; refine res or set a smaller (support-margin mm)"
                .into())
        } else {
            Err("no cell is at least half material: walls about 2 cells thick or thinner are not resolved; refine res".into())
        }
    }

    /// The support margin nearest point `at` (a peak away from the
    /// supports), None where no fix has one (heat, no fixes).
    pub fn nearest_support(&self, at: V3) -> Option<NearestSupport> {
        let (i, d2) = self
            .near_support
            .iter()
            .enumerate()
            .filter(|(_, f)| **f != NO_FIX)
            .map(|(i, _)| {
                let c = self.cell_center(i);
                (i, (0..3).map(|a| (c[a] - at[a]).powi(2)).sum::<f64>())
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))?;
        let distance = d2.sqrt();
        Some(NearestSupport { fix: self.near_support[i] as usize, distance, beside: distance <= self.h * 3f64.sqrt() * (1.0 + 1e-9) })
    }

    /// Peak over reliable cells next to supports (idealised clamp edges).
    pub fn max_near_supports(&self) -> Option<Extreme> {
        let peaks = self.peaks();
        self.cells
            .iter()
            .enumerate()
            .filter(|(i, c)| c.frac >= 0.5 && self.near_support[*i] != NO_FIX)
            .map(|(i, _)| peaks[i])
            .fold(None, |b: Option<(f64, V3)>, x| if b.is_none_or(|b| x.0 > b.0) { Some(x) } else { b })
            .map(|(value, at)| Extreme { value, at })
    }

    /// Largest rotation (radians) at the centres of cells at least half
    /// material: the size of half the curl of the displacement, how far the
    /// solution is from the small rotations linear kinematics assume. None
    /// for heat.
    pub fn max_rotation(&self) -> Option<Extreme> {
        if self.physics.ndof() != 3 {
            return None;
        }
        self.cells
            .iter()
            .enumerate()
            .filter(|(_, c)| c.frac >= 0.5)
            .map(|(i, _)| (i, small_turn(&self.element_u(i), self.cell_size(i))))
            .fold(None, |b: Option<(usize, f64)>, x| if b.is_none_or(|b| x.1 > b.1) { Some(x) } else { b })
            .map(|(i, v)| Extreme { value: v, at: self.cell_center(i) })
    }

    /// Percentiles of the reliable cells' values (see `max_cell_scalar`).
    pub fn cell_scalar_percentiles(&self, qs: &[f64]) -> Result<Vec<f64>, String> {
        // Each cell weighs as its volume: a child (`MAX_DEPTH`) as its
        // share of a grid cell, so splitting cells moves no percentile.
        let mut v: Vec<(f64, f64)> = self.reliable_values()?.into_iter().map(|(i, s)| (s, 0.125f64.powi(self.cells[i].level as i32))).collect();
        v.sort_by(|a, b| a.0.total_cmp(&b.0));
        // Each cell's place: the middle of its weight, from the first's to
        // the last's (equal weights number the cells 0, 1, 2...).
        let first = v.first().map_or(0.0, |x| x.1);
        let (mut at, mut mid) = (0.0, Vec::with_capacity(v.len()));
        for &(_, w) in &v {
            mid.push(at + 0.5 * (w - first));
            at += w;
        }
        let total = mid.last().copied().unwrap_or(0.0);
        Ok(qs
            .iter()
            .map(|q| {
                let t = q * total;
                let k = mid.partition_point(|&m| m < t);
                let k = if k > 0 && (k == mid.len() || t - mid[k - 1] < mid[k] - t) { k - 1 } else { k };
                v[k].0
            })
            .collect())
    }

    /// The solution on its dense regular grid, for visualisation (e.g. upload
    /// as 3D textures).
    pub fn grid_fields(&self) -> GridFields {
        let d = self.dims.map(|v| v as usize);
        let nc = d[0] * d[1] * d[2];
        let (mut scalar, mut frac) = (vec![f32::NAN; nc], vec![0.0f32; nc]);
        for (lin, &ci) in self.cell_at.iter().enumerate() {
            if ci != u32::MAX {
                let c = ci as usize;
                scalar[lin] = cell_scalar(self.physics, &self.cell_field[c]) as f32;
                frac[lin] = self.cells[c].frac as f32;
            }
        }
        let nd = self.physics.ndof();
        let n = [d[0] + 1, d[1] + 1, d[2] + 1];
        let mut node_values = vec![f32::NAN; n[0] * n[1] * n[2] * nd];
        let sc = 1i64 << self.depth;
        for (k, g) in self.node_idx.iter().enumerate().filter(|(_, g)| g.iter().all(|v| v % sc == 0)) {
            let lin = grid::index(g.map(|v| (v / sc) as usize), n);
            let u = self.nodal(k);
            for j in 0..nd {
                node_values[lin * nd + j] = u[j] as f32;
            }
        }
        GridFields { origin: self.origin, h: self.h, dims: d, cell_scalar: scalar, cell_frac: frac, node_values }
    }

    /// Per-cell values, numbered as `Assembly` numbers the cells, on the
    /// grid of `grid_fields` (a cut cell's first piece's), NaN without
    /// material.
    pub fn cells_on_grid(&self, per_cell: &[f64]) -> Vec<f32> {
        let mut v = vec![f32::NAN; self.cell_at.len()];
        for (lin, &ci) in self.cell_at.iter().enumerate() {
            if ci != u32::MAX {
                v[lin] = per_cell[ci as usize] as f32;
            }
        }
        v
    }

    fn cell_center(&self, i: usize) -> V3 {
        place(&self.grid(), &self.cells[i]).0
    }

    fn cell_size(&self, i: usize) -> f64 {
        place(&self.grid(), &self.cells[i]).1
    }

    /// A sharp inside corner of `shape` (a crease: `surface_radius`'s 0, a
    /// round under about 0.4 of a cell too) anywhere on the surface whose
    /// stress is read (`surface_points`, in cells away from the supports):
    /// where one is, the first found in the cells' order. Its elastic
    /// stress is unbounded, whatever the cells read there. `judge`: the
    /// claim's own grid's cells (`split_as`), whose survey a coarser grid
    /// makes too, as a peak's is.
    ///
    /// Surveyed (`surface_radius`, as a peak's surface is) about each point
    /// that some other within two cells of it lies above the tangent plane
    /// of, and it above that one's, their facets' normals turned by more than a
    /// round of `resolved` cells turns over their distance: a flat's
    /// points, a convex edge's or curve's, and a bore's or a fillet's the
    /// cells resolve have none (every point of a 49 mm bore surveyed took 7
    /// minutes), nor has a dent of the mesh at a convex edge, which lifts
    /// neither point as a corner's faces would (a blade's tip edge, convex
    /// all round, was surveyed, and read a crease where its section's
    /// leading edge has no slope). So a corner turning less than such a
    /// round does between the two points nearest it, about 20 degrees, is
    /// not found here.
    pub fn crease(&self, shape: &Arc<Tape>, ctx: &Ctx, judge: f64, resolved: f64) -> Option<V3> {
        let g = self.grid();
        // (at, its facet's normal, the cells it is judged on)
        let pts: Vec<(V3, V3, f64)> = self
            .surface_points
            .iter()
            .zip(&self.surface_normals)
            .filter(|((ci, _), _)| self.near_support[*ci as usize] == NO_FIX)
            .map(|(&(ci, xi), &n)| {
                let (c, size) = place(&g, &self.cells[ci as usize]);
                (at_xi(c, size, xi), n, size * (judge / self.h).min(1.0))
            })
            .collect();
        // A facet with material beyond it is none of the surface: two
        // bodies meeting on a face are meshed there, and solved as one.
        let beyond: Vec<V3> = pts.iter().map(|&(p, n, seen)| v3::add(p, v3::scale(n, seen / 16.0))).collect();
        let mut air = Vec::new();
        eval_points(shape, ctx, &beyond, &mut air);
        let pts: Vec<(V3, V3, f64)> = pts.into_iter().zip(air).filter(|(_, v)| !(*v <= 0.0)).map(|(p, _)| p).collect();
        // Points by the block of two of the coarsest cells that holds them.
        let block = 2.0 * pts.iter().map(|p| p.2).fold(0.0, f64::max);
        let key = |p: V3| -> [i64; 3] { p.map(|v| (v / block).floor() as i64) };
        let mut blocks: std::collections::HashMap<[i64; 3], Vec<usize>> = std::collections::HashMap::new();
        for (i, p) in pts.iter().enumerate() {
            blocks.entry(key(p.0)).or_default().push(i);
        }
        let inside = |i: usize| -> bool {
            let (p, n, seen) = pts[i];
            let k = key(p);
            (0..27).flat_map(|o| blocks.get(&[k[0] + o % 3 - 1, k[1] + o / 3 % 3 - 1, k[2] + o / 9 - 1]).into_iter().flatten()).any(|&j| {
                let (d, m) = (v3::sub(pts[j].0, p), pts[j].1);
                let s = v3::norm(d);
                let turn = v3::dot(n, m).clamp(-1.0, 1.0).acos();
                // Each above the other's plane, and by what two faces
                // turned so would lift them (their heights over their
                // distance at least the turn's sine, were they flats: half
                // of it is asked): a dent of the mesh at a convex edge
                // turns its facets as far and lifts them hardly at all.
                let (up, down) = (v3::dot(d, n), -v3::dot(d, m));
                s <= 2.0 * seen
                    && up > 0.0
                    && down > 0.0
                    && up + down >= 0.5 * turn.min(std::f64::consts::FRAC_PI_2).sin() * s
                    && turn > 2.0 * (0.5 * s / (resolved * seen)).min(1.0).asin()
            })
        };
        (0..pts.len())
            .into_par_iter()
            .find_first(|&i| inside(i) && surface_radius_as(shape, ctx, pts[i].0, 1.5 * pts[i].2, pts[i].2 / 8.0, true).is_some_and(|(_, r)| r == 0.0))
            .map(|i| pts[i].0)
    }

    /// The size of the cell holding `p` (split cells' own), where one does.
    pub fn cell_size_at(&self, p: V3) -> Option<f64> {
        self.cell_of(p).map(|(i, _)| self.cell_size(i))
    }

    /// A grid cell's cells (its pieces, or its children's).
    fn cells_at(&self, i: [i64; 3]) -> std::ops::Range<usize> {
        if (0..3).any(|a| !(0..self.dims[a]).contains(&i[a])) {
            return 0..0;
        }
        match self.cell_at[grid::index(i.map(|v| v as usize), self.dims.map(|v| v as usize))] {
            u32::MAX => 0..0,
            k => k as usize..k as usize + self.cells[k as usize..].iter().take_while(|c| c.idx == i).count(),
        }
    }

    /// Cell `ci`'s nodal displacements.
    fn element_u(&self, ci: usize) -> [f64; 24] {
        let v = self.cell_values(ci);
        std::array::from_fn(|j| v[j])
    }

    /// The cell holding a point, and the point in its coordinates. A point
    /// on a plane of the grid lies in the closed cells on both sides of it
    /// (as the grid places its planes), and a face on the plane has material
    /// on one side only: the first holding material is taken.
    fn cell_of(&self, p: V3) -> Option<(usize, V3)> {
        let g = self.grid();
        let [xs, ys, zs]: [Vec<i64>; 3] = std::array::from_fn(|a| {
            let plane = |i: i64| g.origin[a] + g.h * i as f64;
            let i = ((p[a] - g.origin[a]) / g.h).floor() as i64;
            [i, i - 1, i + 1].into_iter().filter(|&j| (0..self.dims[a]).contains(&j) && plane(j) <= p[a] && p[a] <= plane(j + 1)).collect()
        });
        let (xs, ys) = (&xs, &ys);
        zs.iter().flat_map(|&z| ys.iter().flat_map(move |&y| xs.iter().map(move |&x| [x, y, z]))).find_map(|i| locate(&self.cells, self.cells_at(i), &g, p))
    }

    /// The cell holding a point (`cell_of`), else of the cells round the
    /// one it is in, the nearest holding material, with the point in its
    /// coordinates (outside -1..1): a surface's points the cells leave out
    /// (where a cell's material rounds to none) read their nearest cell's
    /// interpolant, continued.
    fn cell_near(&self, p: V3) -> Option<(usize, V3)> {
        if let Some(c) = self.cell_of(p) {
            return Some(c);
        }
        let g = self.grid();
        let at: [i64; 3] = std::array::from_fn(|a| ((p[a] - g.origin[a]) / g.h).floor() as i64);
        // The nearest cell holding material (a child's, or a grid cell's).
        let mut best: Option<(f64, usize)> = None;
        for d in 0..27 {
            let i: [i64; 3] = std::array::from_fn(|a| at[a] + [d % 3, d / 3 % 3, d / 9][a] - 1);
            for ci in self.cells_at(i) {
                let (c, size) = place(&g, &self.cells[ci]);
                let gap: f64 = (0..3).map(|a| ((p[a] - c[a]).abs() - 0.5 * size).max(0.0).powi(2)).sum();
                if best.is_none_or(|b| gap < b.0) {
                    best = Some((gap, ci));
                }
            }
        }
        let ci = best?.1;
        let (c, size) = place(&g, &self.cells[ci]);
        let xi: V3 = std::array::from_fn(|a| (p[a] - c[a]) * 2.0 / size);
        // Of its pieces, the nearest.
        let (idx, level, sub) = (self.cells[ci].idx, self.cells[ci].level, self.cells[ci].sub);
        let same: Vec<usize> = self.cells_at(idx).filter(|&k| self.cells[k].level == level && self.cells[k].sub == sub).collect();
        Some((same[nearest_piece(same.iter().map(|&k| self.cells[k].mask), xi)], xi))
    }

    /// The solution at a point near the part (within a cell of one holding
    /// material: `cell_near`), so every point of its surface has one, varying
    /// continuously: what a picture of its deformed shape moves its surface by.
    pub fn value_near(&self, p: V3) -> Option<Vec<f64>> {
        let (ci, xi) = self.cell_near(p)?;
        Some(self.interpolate(ci, xi))
    }

    /// Interpolated solution at a point: displacement (3) or temperature (1).
    pub fn value_at(&self, p: V3) -> Option<Vec<f64>> {
        let (ci, xi) = self.cell_of(p)?;
        Some(self.interpolate(ci, xi))
    }

    /// Cell `ci`'s nodal values interpolated at `xi`.
    fn interpolate(&self, ci: usize, xi: V3) -> Vec<f64> {
        let (nd, v) = (self.physics.ndof(), self.cell_values(ci));
        (0..nd).map(|d| (0..8).map(|k| shape_n(k, xi) * v[k * nd + d]).sum()).collect()
    }

    /// σ:ε at a point (twice the strain energy density), from the strain of
    /// the element's displacement field at that point (not the cell centre:
    /// on surfaces in bending the centre value is off by the stress gradient
    /// over half a cell). Elastic solutions only.
    pub fn energy_density(&self, p: V3) -> Option<f64> {
        let Physics::Elastic { .. } = self.physics else { return None };
        let (ci, xi) = self.cell_of(p)?;
        let c = &self.cells[ci];
        let xi = xi.map(|v| v.clamp(-1.0, 1.0));
        let b = b_matrix(xi);
        let dm = d_matrix(c.nu);
        let ue = self.element_u(ci);
        let eps: [f64; 6] = std::array::from_fn(|k| (0..24).map(|j| b[k][j] * ue[j]).sum::<f64>() * 2.0 / self.cell_size(ci));
        if let Some(lam) = &c.lam {
            let sigma = laminate_mean(lam, true, &eps);
            return Some(c.mult * (0..6).map(|k| eps[k] * sigma[k]).sum::<f64>());
        }
        let e = c.coef * c.mult;
        Some((0..6).map(|k| eps[k] * e * (0..6).map(|j| dm[k][j] * eps[j]).sum::<f64>()).sum())
    }

    /// Cell-centre stress tensor (Voigt) or heat flux, and its scalar summary.
    pub fn field_at(&self, p: V3) -> Option<([f64; 6], f64, f64)> {
        let (ci, _) = self.cell_of(p)?;
        let v = self.cell_field[ci];
        Some((v, cell_scalar(self.physics, &v), self.cells[ci].frac))
    }

    /// How far mode `k` (from 1; a buckling analysis's mode is 1) reaches:
    /// the root of the summed variances of its nodes' places, weighted by
    /// their displacement squared. A mode spread over the part keeps it as
    /// the cells shrink; one held to a few cells (at a point the stress has
    /// no bound) shrinks with them.
    pub fn mode_extent(&self, k: usize) -> Option<f64> {
        let shape = self.mode_shapes.get(k.checked_sub(1)?)?;
        let nd = self.physics.ndof();
        let s = spread_of((0..self.node_idx.len()).filter(|n| (n + 1) * nd <= shape.len()).map(|n| (self.node_pos(n), (0..nd).map(|d| shape[n * nd + d].powi(2)).sum::<f64>())));
        Some(s.iter().map(|v| v * v).sum::<f64>().sqrt())
    }

    /// Mode `k` (from 1) of a modal solution as the displacement field, with
    /// its stresses (per 1 mm peak displacement).
    pub fn mode(&self, k: usize) -> Option<Solution> {
        Some(self.shaped(self.mode_shapes.get(k.checked_sub(1)?)?))
    }

    /// The solution with displacement `shape` (per dof) and its stresses.
    fn shaped(&self, shape: &[f64]) -> Solution {
        let mut s = self.clone();
        s.peaks = std::sync::OnceLock::new();
        s.u = shape.to_vec();
        s.held = Vec::new();
        let b0 = b_matrix([0.0; 3]);
        let g = s.grid();
        s.cell_field = s
            .cells
            .iter()
            .map(|c| {
                let ue: Vec<f64> = (0..24).map(|j| s.u[c.nodes[j / 3] as usize * 3 + j % 3]).collect();
                let eps: [f64; 6] = std::array::from_fn(|r| (0..24).map(|j| b0[r][j] * ue[j]).sum::<f64>() * 2.0 / place(&g, c).1);
                if let Some(lam) = &c.lam {
                    return laminate_mean(lam, true, &eps).map(|v| v * c.mult);
                }
                let dm = d_matrix(c.nu);
                std::array::from_fn(|r| c.coef * c.mult * (0..6).map(|q| dm[r][q] * eps[q]).sum::<f64>())
            })
            .collect();
        s
    }

    /// `field_at` at the point itself, not the cell centre: the element's
    /// stress (or flux) varies across the cell, and in bending a surface
    /// point reads about a cell's worth of gradient more than the centre.
    /// The centre value plus the element's change from centre to point.
    pub fn field_at_point(&self, p: V3) -> Option<([f64; 6], f64, f64)> {
        let (ci, xi) = self.cell_of(p)?;
        let v = self.fields_in_cell(ci, &[xi])?[0];
        Some((v, cell_scalar(self.physics, &v), self.cells[ci].frac))
    }

    /// Cell `ci`'s stress (or flux) at each of `xis` (its coordinates,
    /// clamped to within half a cell of it): the centre value plus the
    /// element's change from centre to point (`field_at_point`). A surface
    /// point in a cell under half full is read from the one under it
    /// (`surface_points`), up to half a cell beyond that one's face: read
    /// at the face, a 0.7 mm wall on 0.25 mm cells with its faces 0.1 mm
    /// past the grid's planes peaked at the fibre 0.1 mm in, 86 MPa
    /// against 118 at its surface.
    fn fields_in_cell(&self, ci: usize, xis: &[V3]) -> Option<Vec<[f64; 6]>> {
        if self.cells[ci].lam.is_some() {
            // Each point its material's (`laminate_fields`), or the larger
            // where it may be either.
            let s = |f: &[f64; 6]| cell_scalar(self.physics, f);
            return Some(
                self.laminate_fields(ci, xis, false)?.into_iter().map(|[a, b]| a.into_iter().chain(b).max_by(|p, q| s(p).total_cmp(&s(q))).expect("a material")).collect(),
            );
        }
        let c = &self.cells[ci];
        let xis: Vec<V3> = xis.iter().map(|xi| xi.map(|v| v.clamp(-2.0, 2.0))).collect();
        let v0 = self.cell_field[ci];
        let size = self.cell_size(ci);
        let (scale, t) = (c.coef * c.mult * 2.0 / size, self.cell_values(ci));
        match self.physics {
            Physics::Elastic { .. } => {
                // The element's full strain (with its incompatible modes) at
                // the point, minus at the centre: the change to add.
                let ue = self.element_u(ci);
                let basis = element_basis(self.physics);
                let m = c.material.as_deref().copied().unwrap_or_else(full_moments);
                let mut at = xis;
                at.push([0.0; 3]);
                let e = basis.strains(c.nu, &m, &ue, &at, size)?;
                let e0 = e[at.len() - 1];
                let dm = d_matrix(c.nu);
                Some(
                    e[..at.len() - 1]
                        .iter()
                        .map(|e1| {
                            let mut v = v0;
                            for (k, o) in v.iter_mut().enumerate() {
                                *o += c.coef * c.mult * (0..6).map(|j| dm[k][j] * (e1[j] - e0[j])).sum::<f64>();
                            }
                            v
                        })
                        .collect(),
                )
            }
            Physics::Heat => Some(
                xis.iter()
                    .map(|&xi| {
                        let mut v = v0;
                        for i in 0..8 {
                            let (g1, g0) = (shape_dn(i, xi), shape_dn(i, [0.0; 3]));
                            for a in 0..3 {
                                v[a] -= scale * (g1[a] - g0[a]) * t[i];
                            }
                        }
                        v
                    })
                    .collect(),
            ),
        }
    }

    /// A laminate cell's (`laminate`) stress (or flux) at each of `xis`
    /// (clamped as `fields_in_cell`'s), each of its two materials' where
    /// the point may be read as it (its side of the face; either, where the
    /// cell holds more than one face, or `either`): its own strain,
    /// recovered from the element's, times its own stiffness.
    fn laminate_fields(&self, ci: usize, xis: &[V3], either: bool) -> Option<Vec<[Option<[f64; 6]>; 2]>> {
        let c = &self.cells[ci];
        let lam = c.lam.as_deref()?;
        let xis: Vec<V3> = xis.iter().map(|xi| xi.map(|v| v.clamp(-2.0, 2.0))).collect();
        let (size, t) = (self.cell_size(ci), self.cell_values(ci));
        let strains: Vec<Vec<f64>> = match self.physics {
            Physics::Elastic { .. } => {
                let m = c.material.as_deref().copied().unwrap_or_else(full_moments);
                laminate_strains(lam, &m, &self.element_u(ci), &xis, size)?.into_iter().map(|e| e.to_vec()).collect()
            }
            Physics::Heat => xis
                .iter()
                .map(|&xi| {
                    let mut g = vec![0.0; 3];
                    for i in 0..8 {
                        let dn = shape_dn(i, xi);
                        for a in 0..3 {
                            g[a] += dn[a] * 2.0 / size * t[i];
                        }
                    }
                    g
                })
                .collect(),
        };
        let elastic = matches!(self.physics, Physics::Elastic { .. });
        Some(
            strains
                .iter()
                .zip(&xis)
                .map(|(e, &xi)| {
                    let [a, b] = lam.stresses(elastic, e).map(|f| f.map(|v| v * c.mult));
                    let (ra, rb) = if either { (true, true) } else { lam.at(xi) };
                    [ra.then_some(a), rb.then_some(b)]
                })
                .collect(),
        )
    }

    /// Surface points where a cell's peak may lie between its facets'
    /// centroids (`climbed`): from each cell whose field could pass the
    /// part's peak, its highest surface read is moved along the
    /// surface (its field's gradient across the surface's normal, each step
    /// back onto the zero set by Newton's) while the cell's stress rises and
    /// the cell still takes the point (`surface_cell`). A mesh at half a
    /// cell cuts a convex edge's corner off, so its centroids lie inside
    /// it, and there stress in bending falls off into the part: a square
    /// bar turned 37 degrees and bent read its peak 0.035 mm (0.14 of its
    /// 0.25 mm cells) in from its edge, 4% low, alike on every grid. The
    /// climb follows the surface across an edge as the field's sign alone
    /// places it, onto the edge, where the read is the element's own. Not
    /// a read: on a fillet three cells per radius, coarser grids' climbs
    /// rose nearer the finest's read than their reads had, and the grids,
    /// agreeing, decided a claim 14% under the truth that their reads'
    /// differences had left unknown. It bounds the read from above
    /// (`climbed_scalar`): the solve's own field reaches it.
    pub(crate) fn climb_reads(&mut self, shape: &Arc<Tape>, ctx: &Ctx) {
        let g = self.grid();
        let finest = self.cells.iter().map(|c| c.level).max().unwrap_or(0);
        let peaks = self.peaks().to_vec();
        // The peak away from the supports and within their margins, each
        // read apart (`max_cell_scalar`, `max_near_supports`).
        let near = |i: usize| self.near_support[i] != NO_FIX;
        let top = |n: bool| (0..self.cells.len()).filter(|&i| self.cells[i].frac >= 0.5 && near(i) == n).map(|i| peaks[i].0).fold(f64::NEG_INFINITY, f64::max);
        let tops = [top(false), top(true)];
        // Each candidate cell's highest surface read, as its start.
        let mut start: Vec<Option<(f64, V3)>> = vec![None; self.cells.len()];
        for &(ci, v, at) in self.reads() {
            let ci = ci as usize;
            let surface = at != place(&g, &self.cells[ci]).0;
            if surface && self.cells[ci].frac >= 0.5 && start[ci].is_none_or(|s| v > s.0) {
                start[ci] = Some((v, at));
            }
        }
        // Only cells whose field could pass the part's peak somewhere in
        // their reach (the element's field at its corners bounds it: each
        // strain is affine along every line parallel to an axis, and the
        // scalar convex in them, `cell_upper`): elsewhere the climb could
        // change no peak.
        let corners: Vec<V3> = (0..8).map(|k| std::array::from_fn(|a| if k >> a & 1 == 1 { 2.0 } else { -2.0 })).collect();
        let starts: Vec<(usize, f64, V3)> = start
            .iter()
            .enumerate()
            .filter_map(|(i, s)| s.map(|s| (i, s.0, s.1)))
            .filter(|&(i, _, _)| self.fields_in_cell(i, &corners).is_some_and(|fs| fs.iter().any(|f| !(cell_scalar(self.physics, f) < tops[near(i) as usize]))))
            .collect();
        let found: Vec<(u32, V3)> = starts
            .par_iter()
            .filter_map(|&(ci, v0, p0)| {
                let (ctr, size) = place(&g, &self.cells[ci]);
                let xi = |p: V3| -> V3 { std::array::from_fn(|a| (p[a] - ctr[a]) * 2.0 / size) };
                let inside = |p: V3| xi(p).iter().all(|v| v.abs() <= 2.0);
                let value = |pts: &[V3]| -> Option<Vec<f64>> {
                    let xs: Vec<V3> = pts.iter().map(|&p| xi(p)).collect();
                    self.fields_in_cell(ci, &xs).map(|fs| fs.iter().map(|f| cell_scalar(self.physics, f)).collect())
                };
                // Onto the zero set, and its outward normal; None off it.
                let onto = |p: V3| -> Option<(V3, V3)> {
                    let mut q = p;
                    for _ in 0..4 {
                        let gr = crate::eval::eval_grads(shape, ctx, &[q])[0];
                        let n2 = gr[1] * gr[1] + gr[2] * gr[2] + gr[3] * gr[3];
                        if !(n2 > 0.0) {
                            return None;
                        }
                        q = std::array::from_fn(|a| q[a] - gr[0] * gr[1 + a] / n2);
                    }
                    let gr = crate::eval::eval_grads(shape, ctx, &[q])[0];
                    let n = v3::direction([gr[1], gr[2], gr[3]])?;
                    (gr[0].abs() <= 1e-6 * size * v3::norm([gr[1], gr[2], gr[3]])).then_some((q, n))
                };
                // Its own surface points, and those no cell takes (a convex
                // tip in cells under half material, past the half cell a
                // face is read across) within its reach.
                let takes = |q: V3, n: V3| surface_cell(&self.cells, &|i| self.cells_at(i), &g, finest, q, n, &|c| self.cells[c].frac >= 0.5).is_none_or(|c| c.0 == ci);
                let (mut p, mut n) = onto(p0)?;
                let mut v = value(&[p])?[0];
                if !(v.is_finite() && v >= v0 - 1e-9 * v0.abs() && takes(p, n)) {
                    return None;
                }
                let d = 1e-3 * size;
                let mut step = 0.25 * size;
                for _ in 0..80 {
                    if step < 1e-3 * size {
                        break;
                    }
                    let probe: Vec<V3> = (0..6).map(|k| std::array::from_fn(|a| p[a] + if a == k / 2 { if k % 2 == 0 { d } else { -d } } else { 0.0 })).collect();
                    let Some(fd) = value(&probe) else { break };
                    let grad: V3 = std::array::from_fn(|a| (fd[2 * a] - fd[2 * a + 1]) / (2.0 * d));
                    let gn = v3::dot(grad, n);
                    let Some(t) = v3::direction(std::array::from_fn(|a| grad[a] - gn * n[a])) else { break };
                    let next = onto(std::array::from_fn(|a| p[a] + step * t[a])).filter(|&(q, m)| inside(q) && takes(q, m));
                    match next.and_then(|(q, m)| value(&[q]).map(|w| (q, m, w[0]))) {
                        Some((q, m, w)) if w.is_finite() && w > v => (p, n, v) = (q, m, w),
                        _ => step *= 0.5,
                    }
                }
                (v > v0).then(|| (ci as u32, xi(p)))
            })
            .collect();
        let mut found = found;
        found.sort_by_key(|p| p.0);
        self.climbed = found;
    }

    /// The climbed points' values (`climb_reads`), as (cell, value, point).
    fn climbed_reads(&self) -> Vec<(u32, f64, V3)> {
        let g = self.grid();
        let mut out = Vec::new();
        for grp in self.climbed.chunk_by(|a, b| a.0 == b.0) {
            let ci = grp[0].0 as usize;
            let xis: Vec<V3> = grp.iter().map(|p| p.1).collect();
            let (ctr, size) = place(&g, &self.cells[ci]);
            if let Some(fs) = self.fields_in_cell(ci, &xis) {
                out.extend(
                    fs.iter().zip(&xis).map(|(f, xi)| (ci as u32, cell_scalar(self.physics, f), at_xi(ctr, size, xi.map(|v| v.clamp(-2.0, 2.0))))).filter(|x| x.1.is_finite()),
                );
            }
        }
        out
    }

    /// The peak the surface climbs to from the reads (`climb_reads`) in the
    /// cells `max_cell_scalar` reads, where it lies over that read: the
    /// stress the solve's own field reaches between the surface's reads.
    pub fn climbed_scalar(&self) -> Option<Extreme> {
        let cells: Vec<bool> = (0..self.cells.len()).map(|i| self.cells[i].frac >= 0.5 && self.near_support[i] == NO_FIX).collect();
        self.climbed_reads().into_iter().filter(|r| cells[r.0 as usize]).map(|r| Extreme { value: r.1, at: r.2 }).max_by(|a, b| a.value.total_cmp(&b.value))
    }

    /// `climbed_scalar` over `scale` (as `max_cell_ratio`): the cells it
    /// reads, of one material.
    pub fn climbed_ratio(&self, scale: &dyn Fn(&[V3]) -> Vec<f64>, mixed: Mixed) -> Option<Extreme> {
        let cells: Vec<usize> = self.reliable_values().ok()?.into_iter().map(|c| c.0).collect();
        let mut one = vec![false; self.cells.len()];
        for ci in self.of_one_material(&cells, mixed) {
            one[ci] = true;
        }
        let at: Vec<(f64, V3)> = self.climbed_reads().into_iter().filter(|r| one[r.0 as usize]).map(|r| (r.1, r.2)).collect();
        ratio_peak(&at, scale)
    }

    /// Every point a peak of stress (or heat flux) is read at, as a probe
    /// reads it: each cell's centre, and the points of the part's surface it
    /// holds (`surface_points`; stress in bending and at a concentration
    /// peaks at the surface, the centre lies half a cell under it), as
    /// (cell, value, point), sorted by cell.
    fn reads(&self) -> &[(u32, f64, V3)] {
        ink_kernel::once::get_or_compute::<Vec<_>>(&self.reads, || {
            let g = self.grid();
            let groups: Vec<&[(u32, V3)]> = self.surface_points.chunk_by(|a, b| a.0 == b.0).collect();
            let surface: Vec<Vec<(u32, f64, V3)>> = groups
                .par_iter()
                .map(|grp| {
                    let ci = grp[0].0 as usize;
                    let xis: Vec<V3> = grp.iter().map(|p| p.1).collect();
                    let ctr = self.cell_center(ci);
                    self.fields_in_cell(ci, &xis).map_or(Vec::new(), |fs| {
                        fs.iter()
                            .zip(&xis)
                            .map(|(f, xi)| (ci as u32, cell_scalar(self.physics, f), at_xi(ctr, self.cell_size(ci), xi.map(|v| v.clamp(-2.0, 2.0)))))
                            .filter(|x| x.1.is_finite())
                            .collect()
                    })
                })
                .collect();
            // (A laminate's centre as its material there: `fields_in_cell`.)
            let mut out: Vec<(u32, f64, V3)> = self
                .cells
                .iter()
                .zip(&self.cell_field)
                .enumerate()
                .map(|(ci, (c, f))| {
                    let f = if c.lam.is_some() { self.fields_in_cell(ci, &[[0.0; 3]]).map_or([f64::NAN; 6], |f| f[0]) } else { *f };
                    (ci as u32, cell_scalar(self.physics, &f), place(&g, c).0)
                })
                .collect();
            out.extend(surface.into_iter().flatten());
            out.sort_by_key(|r| r.0);
            out
        })
    }

    /// Per cell, the largest scalar (von Mises stress or flux magnitude)
    /// over its material and where: the largest of its `reads`.
    fn peaks(&self) -> &[(f64, V3)] {
        ink_kernel::once::get_or_compute::<Vec<_>>(&self.peaks, || {
            let mut out = vec![(f64::NEG_INFINITY, [0.0; 3]); self.cells.len()];
            for &(ci, v, at) in self.reads() {
                if v > out[ci as usize].0 {
                    out[ci as usize] = (v, at);
                }
            }
            out
        })
    }

    /// Scalar summary near a point, tolerant of points just outside material
    /// cells (mesh vertices sit on the surface): searches the 3x3x3 neighbourhood.
    pub fn scalar_near(&self, p: V3) -> Option<f64> {
        if let Some((_, s, f)) = self.field_at(p)
            && f >= 0.25
        {
            return Some(s);
        }
        let mut best: Option<(f64, f64)> = None;
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let q = [p[0] + dx as f64 * self.h, p[1] + dy as f64 * self.h, p[2] + dz as f64 * self.h];
                    if let Some((ci, _)) = self.cell_of(q) {
                        let c = &self.cells[ci];
                        if c.frac > best.map_or(0.0, |b| b.1) {
                            best = Some((cell_scalar(self.physics, &self.cell_field[ci]), c.frac));
                        }
                    }
                }
            }
        }
        best.map(|b| b.0)
    }
}

/// The largest of values over `scale` at their points.
fn ratio_peak(at: &[(f64, V3)], scale: &dyn Fn(&[V3]) -> Vec<f64>) -> Option<Extreme> {
    let pts: Vec<V3> = at.iter().map(|a| a.1).collect();
    let sc = scale(&pts);
    at.iter().zip(sc).map(|(&(v, p), s)| (v / s, p)).filter(|x| !x.0.is_nan()).max_by(|a, b| a.0.total_cmp(&b.0)).map(|(value, at)| Extreme { value, at })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink_core::graph;

    /// The radius a peak's surface reads on 1 mm cells (`surface_radius`
    /// as a study reads it: reach 1.5 cells, step an eighth), of the field
    /// `f` made from x, y and constants.
    fn radius_of(f: impl Fn(&mut graph::Graph, graph::NodeId, graph::NodeId) -> graph::NodeId, p: V3) -> Option<f64> {
        let mut g = graph::Graph::new();
        let (x, y) = (g.var(graph::Var::X), g.var(graph::Var::Y));
        let root = f(&mut g, x, y);
        let ctx = Ctx::new(0.0, Vec::new());
        surface_radius(&Arc::new(Tape::new(&g, root)), &ctx, p, 1.5, 0.125).map(|r| r.1)
    }

    /// 1 - |(x, y) - c| / r scaled: a bore of radius `r` about `c`, material
    /// outside.
    fn bore(g: &mut graph::Graph, x: graph::NodeId, y: graph::NodeId, c: [f64; 2], r: f64) -> graph::NodeId {
        let (cx, cy) = (g.constant(c[0]), g.constant(c[1]));
        let (dx, dy) = (g.sub(x, cx), g.sub(y, cy));
        let (dx2, dy2) = (g.un(graph::UnOp::Square, dx), g.un(graph::UnOp::Square, dy));
        let d2 = g.add(dx2, dy2);
        let d = g.un(graph::UnOp::Sqrt, d2);
        let rr = g.constant(r);
        g.sub(rr, d)
    }

    #[test]
    fn a_surface_reads_its_zero_set_s_radius_whatever_its_field() {
        // A bore of radius 1 written f, f / (1 + 1e6 f^2) and f^3: one
        // surface. The normal's turn read off it gave 1, infinity (the
        // second's gradient turns back off its zero set) and no surface
        // (the third's is 0 on it).
        let plain = radius_of(|g, x, y| bore(g, x, y, [0.0, 0.0], 1.0), [0.0, 1.0, 0.0]);
        let damped = radius_of(
            |g, x, y| {
                let f = bore(g, x, y, [0.0, 0.0], 1.0);
                let (one, k) = (g.constant(1.0), g.constant(1e6));
                let f2 = g.un(graph::UnOp::Square, f);
                let kf2 = g.mul(k, f2);
                let den = g.add(one, kf2);
                g.bin(graph::BinOp::Div, f, den)
            },
            [0.0, 1.0, 0.0],
        );
        let cubed = radius_of(
            |g, x, y| {
                let f = bore(g, x, y, [0.0, 0.0], 1.0);
                let f2 = g.mul(f, f);
                g.mul(f2, f)
            },
            [0.0, 1.0, 0.0],
        );
        for r in [plain, damped, cubed] {
            assert!(r.is_some_and(|r| (r - 1.0).abs() < 0.03), "{plain:?} {damped:?} {cubed:?}");
        }
    }

    #[test]
    fn a_crease_turning_its_normal_over_a_thousandth_of_a_radian_reads_no_radius() {
        // min(y, y - a x) at its crease: a V, its radius 0, its normal
        // turning by atan(a). The normal's turn over a step read 2.5 at a =
        // 0.1 and 25 at 0.01: resolved. Under 0.001 rad it reads as a curve,
        // its bend at an eighth of a cell: resolved, as stated.
        let read = |a: f64| {
            radius_of(
                |g, x, y| {
                    let k = g.constant(a);
                    let ax = g.mul(k, x);
                    let tilted = g.sub(y, ax);
                    g.min(y, tilted)
                },
                [0.0, 0.0, 0.0],
            )
        };
        for a in [1.0, 0.1, 0.01, 0.002] {
            assert!(read(a).is_some_and(|r| r == 0.0), "{a}: {:?}", read(a));
        }
        assert!(read(0.0005).is_some_and(|r| r > 1.45), "{:?}", read(0.0005));
    }

    #[test]
    fn the_tightest_concave_surface_within_reach_is_found() {
        // A flat face y = 0 with a bore of radius 0.6 at (1, -0.9) under it,
        // 0.75 from a peak on the face: within 1.5 cells. One projection
        // reached the face, flat: no radius.
        let r = radius_of(
            |g, x, y| {
                let b = bore(g, x, y, [1.0, -0.9], 0.6);
                g.max(y, b)
            },
            [0.0, 0.0, 0.0],
        );
        assert!(r.is_some_and(|r| (r - 0.6).abs() < 0.03), "{r:?}");
    }
}
