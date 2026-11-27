//! Field-based finite element analysis: the finite cell method on a regular grid.
//!
//! No mesh is generated. The domain is a grid of cubes (cells) of size `h`:
//! - cells the octree proves fully inside are trilinear hex elements (with
//!   Wilson's incompatible modes, for elasticity);
//! - cells the surface cuts are integrated by sampling the shape field, so holes
//!   and lattices are captured without meshing.
//!
//! Every element integrand here (B^T D B for elasticity, grad N . grad N for
//! heat) is a polynomial of degree <= 2 in each local coordinate. So a cell's
//! matrix is an exact linear combination of 27 constant matrices, weighted by
//! the moments  m_abc = integral of (coef * inside) xi^a eta^b zeta^c.
//! Material coefficients may be fields: they enter through the moments.
//!
//! The system is solved matrix-free by conjugate gradients preconditioned
//! with geometric multigrid, on all cores (cells are 8-coloured so scatters
//! never collide).
//!
//! The `study` form's solver (ink-sim): a solve reads its `Setup`
//! and context alone, so it is cached by them (`Setup::key`).

mod buckle;
mod cells;
mod eigen;
mod element;
mod hold;
mod laminate;
mod mechanisms;
mod memory;
mod multigrid;
mod numerics;
mod occupancy;
mod report;
mod slivers;
mod solution;
mod solve;
mod surface;
mod system;
mod transient;

// The crates it builds on, at the paths the solver's code uses.
use ink_core::{V3, eval, graph, interval, tape, v3};
use ink_kernel::{cache, mesh, octree, query};

pub use crate::cells::{VOID_SHARE, variation_resolved};
use crate::eval::{Ctx, eval_interval};
use crate::occupancy::tangent_offsets;
pub use crate::report::{
    Buckle, Incomplete, Report, SurfaceMeasure, cells_of, glued_warning, interface_shared_what, interface_stepped_what, material_left_warning, material_left_what, spread_of,
    thin_left_warning,
};
pub use crate::solution::{Extreme, GridFields, Mixed, MixedCell, NearestSupport, SUPPORT_MARGIN_MM, Solution, buckle_of_sum, surface_radius, von_mises};
pub use crate::solve::{Assembly, SolveError, buckled_warning, solve, solve_checked, solve_redesigned};
use crate::surface::{Facet, facet_samples, region_share};
use crate::tape::Tape;
use rayon::prelude::*;
use std::sync::Arc;

// Least work per parallel task: dof values, cells of a colour, nodes of a
// transfer. Coarse levels are small; split finer, their tasks wait on each
// other more than they compute.
const GRAIN: usize = 4096;

/// Most iterations of a CG solve.
const MAX_ITER: usize = 200_000;
/// Most iterations of the modal eigensolve (LOBPCG).
const MODAL_ITER: usize = 300;
const CELL_GRAIN: usize = 64;
const NODE_GRAIN: usize = 1024;

#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Physics {
    /// Linear elasticity. Units: mm, N, MPa. Coefficient field: Young's modulus E.
    Elastic { nu: f64 },
    /// Steady heat conduction. Units: mm, W, K. Coefficient field: conductivity in W/(mm K).
    Heat,
}

impl Physics {
    pub fn ndof(&self) -> usize {
        match self {
            Physics::Elastic { .. } => 3,
            Physics::Heat => 1,
        }
    }
}

/// Where a support or load acts.
pub enum RegionSel {
    /// Material where this field is negative.
    Field(Arc<Tape>),
    /// Material within `depth` under the surface where provenance names `tag`
    /// (and, if given, the outward normal is within 30 degrees of `normal`).
    Face { shape: Arc<Tape>, tag: u32, normal: Option<V3>, depth: f64 },
}

impl RegionSel {
    /// Can this region intersect the box? (Conservative.)
    fn may_touch(&self, ctx: &Ctx, b: [crate::interval::I; 3]) -> bool {
        match self {
            RegionSel::Field(t) => eval_interval(t, ctx, b, None).lo < 0.0,
            RegionSel::Face { shape, depth, .. } => {
                let r = eval_interval(shape, ctx, b, None);
                r.lo < 0.0 && r.hi > -depth
            }
        }
    }

    /// A support includes its boundary when deciding whether it touches
    /// a mirror. A face selector can touch at exactly zero field value.
    fn touches_plane(&self, ctx: &Ctx, plane: [interval::I; 3]) -> bool {
        match self {
            Self::Field(t) => eval_interval(t, ctx, plane, None).lo <= 0.0,
            Self::Face { shape, depth, .. } => {
                let r = eval_interval(shape, ctx, plane, None);
                r.lo <= 0.0 && r.hi >= -depth
            }
        }
    }

    /// Whether each point belongs to the face within a band of material
    /// `[-band, slack]` around its surface.
    fn face_members(shape: &Tape, tag: u32, normal: Option<V3>, ctx: &Ctx, pts: &[V3], band: f64, slack: f64) -> Vec<bool> {
        let g = crate::eval::eval_grads(shape, ctx, pts);
        pts.iter()
            .zip(&g)
            .map(|(p, g)| {
                if !(g[0] >= -band && g[0] <= slack) {
                    return false;
                }
                if let Some(n) = normal {
                    let l = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
                    if !(l > 0.0 && (g[1] * n[0] + g[2] * n[1] + g[3] * n[2]) / l > libm::cos(30f64.to_radians())) {
                        return false;
                    }
                }
                crate::eval::provenance(shape, ctx, *p).contains(&tag)
            })
            .collect()
    }

    /// Weight of each sample (0..1) as part of the region, and the offset
    /// (in xi) of the region's part of its sub-cell from the sub-cell's
    /// centre, as `occupied_parts`: a load on a thin region acts
    /// where that region is. Held at the sub-cell's centre, a symmetry
    /// plane `(- x 0.01)` held its cells' first samples, h/8 inside: an
    /// eighth of the reaction went to the next layer of nodes, and the
    /// first cell read 8.7 MPa in a bar at 10, at every cell size.
    fn occupied(&self, ctx: &Ctx, pts: &[V3], h: f64) -> Vec<(f64, V3)> {
        match self {
            RegionSel::Field(t) => tangent_offsets(t, ctx, pts, h).0.into_iter().map(|(o, off)| (o, off.map(|v| v * 2.0 / h))).collect(),
            RegionSel::Face { shape, tag, normal, depth } => {
                Self::face_members(shape, *tag, *normal, ctx, pts, *depth, 0.0).into_iter().map(|b| (if b { 1.0 } else { 0.0 }, [0.0; 3])).collect()
            }
        }
    }

    /// Each facet's share of the region, and the centroid of that share:
    /// a field region's part of the facet (of its four points on the
    /// surface, `surface_samples`, those inside it), a face's whole facet where its
    /// centroid is on the face within its depth on either side (a facet
    /// lies on the surface to the mesh's accuracy, outside it where it is
    /// concave). The share of the facet, not of a box round its centroid:
    /// a sub-cell's box at a face is half outside the part, so a region
    /// drawn as a thin slab under a face, `(- x 0.01)`, caught half the
    /// face (50 N of a 100 N pressure; a convection's rise 1.3 to 1.6
    /// times too high), and one wholly beyond the part caught half of the
    /// face it met.
    fn on_surface(&self, ctx: &Ctx, facets: &[Facet], samples: &[V3], h: f64) -> Vec<(f64, V3)> {
        match self {
            RegionSel::Face { shape, tag, normal, depth } => {
                let pts: Vec<V3> = facets.iter().map(|f| f.at).collect();
                Self::face_members(shape, *tag, *normal, ctx, &pts, *depth, *depth).into_iter().zip(&pts).map(|(b, p)| (if b { 1.0 } else { 0.0 }, *p)).collect()
            }
            RegionSel::Field(t) => {
                let tie = 1e-9 * h;
                let flat = facet_samples(facets);
                facets
                    .par_iter()
                    .enumerate()
                    .map(|(i, f)| {
                        // The facet's patch of the surface: its corners moved along its normal as its samples
                        // were, by the affine fit through its three outer ones (the inverse of their
                        // barycentric weights is 2 I - J / 3), and how far its centre's sits off that.
                        let off: [f64; 4] = std::array::from_fn(|k| v3::dot(v3::sub(samples[4 * i + k], flat[4 * i + k]), f.normal));
                        let mean = (off[0] + off[1] + off[2]) / 3.0;
                        let corner: [V3; 3] = std::array::from_fn(|k| v3::add(f.corners[k], v3::scale(f.normal, 2.0 * off[k] - mean)));
                        let slack = (off[3] - mean).abs() + tie;
                        region_share(t, ctx, corner, slack, tie)
                    })
                    .collect()
            }
        }
    }

    /// Whether each point of the surface lies in the region, as
    /// `on_surface` tests its samples.
    fn holds_on_surface(&self, ctx: &Ctx, pts: &[V3], h: f64) -> Vec<bool> {
        match self {
            RegionSel::Face { shape, tag, normal, depth } => Self::face_members(shape, *tag, *normal, ctx, pts, *depth, *depth),
            RegionSel::Field(t) => crate::eval::eval_points_par(t, ctx, pts).into_iter().map(|v| v <= 1e-9 * h).collect(),
        }
    }

    /// Whether each point lies within `margin` of the region.
    fn near(&self, ctx: &Ctx, pts: &[V3], margin: f64) -> Vec<bool> {
        match self {
            RegionSel::Field(t) => crate::eval::eval_points_par(t, ctx, pts).into_iter().map(|v| v < margin).collect(),
            RegionSel::Face { shape, tag, normal, depth } => {
                pts.par_chunks(1024).flat_map_iter(|c| Self::face_members(shape, *tag, *normal, ctx, c, depth + margin, margin)).collect()
            }
        }
    }
}

/// A value a fix holds: a number, or a field of x, y, z.
#[derive(Clone)]
pub enum Held {
    Number(f64),
    Field(Arc<Tape>),
}

/// An ordinary region hold may coincide with a symmetry cut. Only a
/// declared mirror asserts that the geometry, supports and loads reflect.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Boundary {
    Unspecified,
    Support,
    Mirror { axis: usize, at: f64 },
}

/// A Dirichlet condition: material in `region` is held at `value`.
pub struct Fix {
    pub boundary: Boundary,
    pub name: String,
    pub region: RegionSel,
    pub dofs: [bool; 3],
    /// What it holds: per axis the displacement (elastic, mm: a prescribed
    /// motion, zero for a plain support), or the temperature (heat: the
    /// first; with `film`, the ambient's).
    pub value: [Held; 3],
    /// Convection (heat): film coefficient in W/(mm^2 K). The fix is then a
    /// Robin condition on the surface inside `region`: flux h (value - T),
    /// instead of a rigid hold.
    pub film: Option<f64>,
    /// Hold only the displacement along the surface normal (frictionless).
    pub normal: bool,
    /// Rigid motions of the region left free (elastic): on the region, what
    /// the fix holds is held relative to the best fitting combination of
    /// them, whose amplitudes are unknowns of the solve.
    pub free: Vec<Motion>,
    /// With free motions: hold the region rigidly but for them, pointwise,
    /// not its mean motion (`distributing`).
    pub rigid: bool,
}

/// A rigid motion, as its displacement field per unit amplitude:
/// omega x (p - at) + v. A turn about an axis (omega a unit vector: per
/// radian) or a slide (v a unit vector: per mm).
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Motion {
    pub omega: V3,
    pub at: V3,
    pub v: V3,
}

impl Motion {
    /// A turn about the axis through `at` along the unit vector `axis`.
    pub fn turn(at: V3, axis: V3) -> Motion {
        Motion { omega: axis, at, v: [0.0; 3] }
    }

    /// A slide along the unit vector `dir`.
    pub fn slide(dir: V3) -> Motion {
        Motion { omega: [0.0; 3], at: [0.0; 3], v: dir }
    }

    /// The displacement at `p` per unit amplitude.
    pub fn at(&self, p: V3) -> V3 {
        let r: V3 = std::array::from_fn(|k| p[k] - self.at[k]);
        let w = self.omega;
        [w[1] * r[2] - w[2] * r[1] + self.v[0], w[2] * r[0] - w[0] * r[2] + self.v[1], w[0] * r[1] - w[1] * r[0] + self.v[2]]
    }

    /// What it is, in words, with its unit.
    pub fn describe(&self) -> String {
        let f = |v: V3| format!("({})", v.map(|c| format!("{}", (c * 1e6).round() / 1e6 + 0.0)).join(", "));
        if self.omega == [0.0; 3] { format!("slide along {} (mm)", f(self.v)) } else { format!("turn about {} through {} (rad)", f(self.omega), f(self.at)) }
    }
}

/// A load of total size `total` spread over the material inside `region`,
/// proportional to volume (force in N, or heat power in W).
pub struct Load {
    pub name: String,
    pub region: RegionSel,
    pub total: [f64; 3],
    /// Totals as functions of time (transient studies); `total` is unused then.
    pub total_t: Option<[Arc<Tape>; 3]>,
    /// Surface pressure (MPa): traction -p n on the surface inside `region`
    /// (`total` is then unused).
    pub pressure: Option<Arc<Tape>>,
    /// Torque (N mm) about the region's centroid, as the forces of a rigid
    /// rotation of its material (`total` is then unused).
    pub torque: Option<[f64; 3]>,
}

/// Heat over time, from t = 0 to `duration` in `steps` equal steps of
/// backward Euler.
#[derive(Clone, Copy, Debug)]
pub struct Transient {
    pub duration: f64,
    pub steps: usize,
    /// The initial temperature everywhere.
    pub initial: f64,
    /// Specific heat capacity, J/(kg K) (with the density field).
    pub capacity: f64,
}

pub struct Setup {
    pub physics: Physics,
    /// Shape field, already clipped to a bounded region.
    pub shape: Arc<Tape>,
    /// E (MPa) or conductivity (W/(mm K)); `coef_const` short-circuits sampling.
    pub coef: Arc<Tape>,
    pub coef_const: Option<f64>,
    /// Elastic: Poisson's ratio as a field, sampled per cell (None: the
    /// physics' constant nu).
    pub nu: Option<Arc<Tape>>,
    /// Elastic: the yield strength (MPa) where it is a field. The solve
    /// does not read it, but a claim reads it point by point, so a cell
    /// it varies across is split as one a coefficient varies across is
    /// (`PROPERTY`).
    pub strength: Option<Arc<Tape>>,
    /// Body force density in N/mm^3 per unit of the `body` field (e.g. rho * g).
    pub body: Option<(Arc<Tape>, [f64; 3])>,
    pub fixes: Vec<Fix>,
    pub loads: Vec<Load>,
    pub lo: V3,
    pub hi: V3,
    pub h: f64,
    /// The study's nominal grid for proving planar pressure faces, shared
    /// by every comparison grid, including time-limited reads.
    pub pressure_grid: f64,
    pub tol: f64,
    /// Stresses closer than this to a fix region are reported separately (mm).
    pub support_margin: f64,
    /// Modal analysis: number of vibration modes to find (0 = static solve).
    pub modes: usize,
    /// Density (kg/m^3) for modal mass and heat capacity.
    pub density: Option<Arc<Tape>>,
    /// Heat over time (needs `density`).
    pub transient: Option<Transient>,
    /// Cells split where a grid of cells this size (mm) would split them,
    /// as deep (`THIN`), and no thin feature left refused: a claim's other
    /// grids, whose cells then differ from its own everywhere as their
    /// sizes do (a thin feature it leaves only widens their difference).
    pub split_as: Option<f64>,
    /// Its buckling analysis run for its loads reversed too (`Report::buckle_reversed`): a solve
    /// summed with others' (`Solution::set_sum`), where a weight may be negative.
    pub both_ways: bool,
}

impl Setup {
    /// The key of this solve in the cache (`cache::key`): the question as
    /// the solver reads it: its tapes and the context they read
    /// (`cache::question`), every setting, the key of the temperatures it
    /// reads; not the tolerance (a tighter cached solve serves a looser
    /// one) nor names (a cached solve reports those of the setup that
    /// first ran it). So two models that build the same setup share the
    /// solve. With the context projected to what the tapes
    /// read, to `solve` in: its answer cannot depend on what the key leaves
    /// out. None: the question cannot be hashed (a derivation that is not
    /// data: Derive::fingerprint).
    pub fn key(&self, ctx: &Ctx) -> Option<(u128, Ctx)> {
        self.key_with(ctx, "", &[])
    }

    /// The key of this solve as a redesign of it (`solve_redesigned`)
    /// asks it: the redesign's settings as `words`, and the fields it reads
    /// (`tapes`), whose inputs the projected context keeps too.
    pub fn key_with(&self, ctx: &Ctx, words: &str, tapes: &[&Tape]) -> Option<(u128, Ctx)> {
        self.key_as(ctx, words, tapes, self.both_ways)
    }

    /// `key`, of its solve with the buckling analysis both ways (`both_ways`), which holds all its own does.
    pub fn key_both_ways(&self, ctx: &Ctx) -> Option<(u128, Ctx)> {
        self.key_as(ctx, "", &[], true)
    }

    fn key_as(&self, ctx: &Ctx, words: &str, tapes: &[&Tape], both_ways: bool) -> Option<(u128, Ctx)> {
        // Every field, named: one added must be keyed too.
        let Setup {
            physics,
            shape,
            coef,
            coef_const,
            nu,
            strength,
            body,
            fixes,
            loads,
            lo,
            hi,
            h,
            pressure_grid,
            tol: _,
            support_margin,
            modes,
            density,
            transient,
            split_as,
            both_ways: _,
        } = self;
        let mut fp = graph::Fingerprint::default();
        fp.str(words);
        if both_ways {
            fp.str("both ways");
        }
        let gravity = body.as_ref().map(|b| b.1);
        fp.str(&format!("{physics:?} {coef_const:?} {gravity:?} {lo:?} {hi:?} {h:?} {pressure_grid:?} {support_margin:?} {modes} {transient:?} {split_as:?}"));
        // Each tape in its slot, an optional one present or not.
        let mut slots: Vec<Option<&Tape>> = vec![Some(shape), Some(coef), nu.as_deref(), strength.as_deref(), density.as_deref(), body.as_ref().map(|b| &*b.0)];
        fp.word(fixes.len() as u64).word(loads.len() as u64);
        for Fix { name: _, region: _, boundary, dofs, value, film, normal, free, rigid } in fixes {
            fp.str(&format!("{boundary:?} {dofs:?} {film:?} {normal} {free:?} {rigid}"));
            for v in value {
                match v {
                    Held::Number(x) => _ = fp.word(0).f64(*x),
                    Held::Field(t) => {
                        fp.word(1);
                        slots.push(Some(t));
                    }
                }
            }
        }
        for Load { name: _, region: _, total, total_t, pressure, torque } in loads {
            fp.str(&format!("{total:?} {torque:?}"));
            slots.extend((0..3).map(|k| total_t.as_ref().map(|t| &*t[k])).chain([pressure.as_deref()]));
        }
        // A face's tag is numbered in the question, once it is hashed.
        let mut faces = Vec::new();
        for r in fixes.iter().map(|f| &f.region).chain(loads.iter().map(|l| &l.region)) {
            slots.push(Some(match r {
                RegionSel::Field(t) => {
                    fp.word(0);
                    t
                }
                RegionSel::Face { shape, tag, normal, depth } => {
                    fp.word(1).str(&format!("{normal:?} {depth:?}"));
                    faces.push(*tag);
                    shape
                }
            }));
        }
        slots.extend(tapes.iter().map(|t| Some(*t)));
        fp.word(slots.len() as u64);
        slots.iter().for_each(|s| _ = fp.word(s.is_some() as u64));
        let tapes: Vec<&Tape> = slots.into_iter().flatten().collect();
        let ctx = ctx.project(&tapes);
        let names = cache::question(&mut fp, &tapes, &ctx)?;
        faces.iter().for_each(|&t| _ = fp.word(names.tag(t) as u64));
        Some((cache::key(fp.value()), ctx))
    }
}

/// `run` of each of `jobs`, at most `n` of them at a time, on the solver's threads: each result in its
/// job's place.
pub fn at_once<J: Sync, T: Send>(n: usize, jobs: &[J], run: impl Fn(&J) -> T + Sync) -> Vec<T> {
    jobs.chunks(n.max(1)).flat_map(|some| some.par_iter().map(&run).collect::<Vec<T>>()).collect()
}

/// How many solves that each take about `need` bytes (`Report::memory_estimate`)
/// fit at once in the memory this process may still take: at least one.
/// The estimate is 0.43 to 1.35 times what a solve takes, so they are let
/// fill 0.4 of it.
pub fn solves_at_once(need: f64) -> usize {
    memory::available().map_or(1, |have| (0.4 * have / need.max(1.0)) as usize).max(1)
}
