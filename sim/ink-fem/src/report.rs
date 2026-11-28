//! What a solve says about itself: its report, what left it incomplete, and
//! the words of its warnings.

use crate::cells::THIN;
use crate::octree::Grid;
use crate::surface::{Facet, Patch, SPLIT_BLOCK};
use ink_core::V3;
use rayon::prelude::*;

/// Keep load and fix positions in cached diagnostics. A second study
/// can share the same numerical inputs under different names.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub enum Incomplete {
    Separation(String),
    Pressure {
        applied: usize,
        at: V3,
        lost: f64,
        variation: f64,
        meant: f64,
    },
    DuplicateCut {
        fix: usize,
    },
    /// Two fixes hold different axes of one material, and one moves it
    /// with free motions (`owners`).
    SharedHold {
        fixes: [usize; 2],
        at: V3,
    },
    /// Free material a fix's edge leaves that the cells cannot carry: held
    /// with it as a sliver, or a lip thinner than a sample between its
    /// edge and the part's face.
    Held(String),
}

impl Incomplete {
    pub fn reason(&self, report: &Report) -> String {
        match self {
            Self::Separation(reason) | Self::Held(reason) => reason.clone(),
            Self::Pressure { applied, at, lost, variation, meant } => format!(
                "pressure load '{}' near {at:?} mm: on this facet {lost:.4e} N missed the cells and pressure variation between samples allows another {variation:.4e} N, against {meant:.4e} N of sampled magnitude (fraction {:.5}). Refine (res r) across the loaded surface",
                report.applied[*applied].0,
                if *meant > 0.0 { (lost + variation) / meant } else { f64::INFINITY }
            ),
            Self::DuplicateCut { fix } => {
                format!("fix '{}' holds the same cut as a declared mirror. Remove the leftover fix, or mark it (support) if it is a physical restraint", report.reactions[*fix].0)
            }
            Self::SharedHold { fixes: [a, b], at } => format!(
                "fixes '{}' and '{}' hold different axes of the same material near {at:?} mm, and one of them moves it with a free motion (a pivot or a slide): \
                 the solver does not couple the two there, so the study is not settled. Hold those axes with one fix, or keep the regions apart",
                report.reactions[*a].0, report.reactions[*b].0
            ),
        }
    }
}

#[derive(Default, Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Report {
    pub cells_full: usize,
    pub cells_cut: usize,
    /// Of them, children of grid cells split where a feature thinner than
    /// two cells lay in them (`THIN`, `MAX_DEPTH`), and how many times finer
    /// the finest are (2^depth).
    #[serde(default)]
    pub cells_refined: usize,
    #[serde(default)]
    pub refine_depth: u8,
    pub nodes: usize,
    pub dofs: usize,
    /// Material volume held by fix regions (mm^3).
    pub clamped_volume: f64,
    /// The holds' boundary with the free material, as traced (mm^2).
    #[serde(default)]
    pub held_edge: f64,
    pub dropped_cells: usize,
    pub dropped_volume: f64,
    #[serde(default)]
    pub dropped_at: Option<V3>,
    /// A model discrepancy that prevents whole-study claims (also sampled fields).
    #[serde(default)]
    pub incomplete: Vec<Incomplete>,
    /// Indexed like applied (including gravity). None for a non-pressure
    /// load. Maximum facet error divided by that facet's sampled force.
    /// Unbounded when the facet can be loaded but its sample reads zero.
    #[serde(default)]
    pub pressure_error: Vec<Option<f64>>,
    /// Material in sliver cells (under MIN_FRACTION full) left out (mm^3).
    pub sliver_volume: f64,
    pub volume: f64,
    /// Material in cells under half full: features about the cell size or
    /// thinner. Their stiffness is averaged in; their peak stresses are not
    /// resolved, and they are excluded from stress statistics.
    pub unresolved_volume: f64,
    /// Iterations of the solve that gives the answer: CG, or the
    /// eigensolver's (LOBPCG for modes).
    pub iterations: usize,
    /// Multigrid levels (1 = plain Jacobi preconditioning).
    pub levels: usize,
    pub residual: f64,
    pub converged: bool,
    /// Why the solve gave nothing to read, if it failed: it broke down (a
    /// value not finite) or did not reach its tolerance. Then no metric,
    /// sampled result or claim is taken from it (`Solution::valid`).
    #[serde(default)]
    pub failure: Option<String>,
    pub setup_ms: f64,
    pub solve_ms: f64,
    /// Applied load totals (N or W) and reactions per fix (same units).
    pub applied: Vec<(String, [f64; 3])>,
    pub reactions: Vec<(String, [f64; 3])>,
    /// Applied load moments and reaction moments about the origin (N mm), elastic studies only,
    /// indexed identically to `applied` and `reactions`.
    #[serde(default)]
    pub applied_moments: Vec<[f64; 3]>,
    #[serde(default)]
    pub reaction_moments: Vec<[f64; 3]>,
    /// u^T f, f the loads alone (not what the fixes hold): compliance (N
    /// mm) for elasticity.
    pub energy: f64,
    /// Modal analysis: per mode, frequency (Hz) and effective mass fraction
    /// along x, y, z.
    /// Frequency (Hz) and effective mass fractions: translations x y z,
    /// rotations about the centre of mass x y z.
    pub modes: Vec<(f64, [f64; 6])>,
    /// Modal analysis: the rigid motions the supports leave free, each a
    /// mode at 0 Hz, taken out of `modes`.
    #[serde(default)]
    pub rigid_modes: usize,
    /// Total mass (kg) for modal analysis.
    pub mass_kg: f64,
    /// Served from the result cache (timings are those of the original solve).
    #[serde(default)]
    pub cached: bool,
    /// Transient: (t, headline value) per step: peak displacement magnitude
    /// (mm) or peak temperature over reliable nodes.
    pub history: Vec<(f64, f64)>,
    /// Transient: the time whose state is reported as the solution (the
    /// headline peak).
    pub peak_time: Option<f64>,
    /// Transient: the peak of the history, (t, value, where): the node it
    /// is at, at its moment (the fields are the final state's).
    #[serde(default)]
    pub peak: Option<(f64, f64, V3)>,
    /// Transient: its time check, the same run on steps twice and four
    /// times as long (as many as halving allows), coarsest first, then its
    /// own: (steps taken, the history's peak, its longest step, s).
    #[serde(default)]
    pub time_levels: Vec<(usize, f64, f64)>,
    /// Supports with free rigid motions: per fix, each motion and the
    /// amplitude it moved by (rad for a turn, mm for a slide).
    #[serde(default)]
    pub free_motions: Vec<(String, Vec<(String, f64)>)>,
    /// What the solve found and the reader must know (a mechanism held at rest).
    #[serde(default)]
    pub warnings: Vec<String>,
    /// Static elastic: its linear buckling analysis, run on every solve
    /// that converged (none where it did not, or in a redesign's).
    #[serde(default)]
    pub buckle: Option<Buckle>,
    #[serde(default)]
    pub modeled_buckle: Option<Buckle>,
    /// `buckle`, of the loads reversed, where the solve was asked it (`Setup::both_ways`).
    #[serde(default)]
    pub buckle_reversed: Option<Buckle>,
    #[serde(default)]
    pub parity_sectors: usize,
    #[serde(default)]
    pub possible_cuts: Vec<(usize, usize, f64)>,
    /// Fixes that hold one direction only at nothing, say neither `(support)` nor a mirror, leave no
    /// pivot or slide free, and are no whole flat face (`possible_cuts`): each may be a cut of symmetry, whose antisymmetric
    /// modes a study of the part as cut does not hold.
    #[serde(default)]
    pub unstated_holds: Vec<usize>,
    /// A gap the cells glue shut (`glued_gap`): how many cells glue it,
    /// and the narrowest's middle and width (mm). The study solved a part
    /// with the gap filled, so a claim on it is undecided.
    #[serde(default)]
    pub glued: Option<(usize, V3, f64)>,
    /// Cells left holding a thin feature when splitting stopped (`THIN`,
    /// `MAX_DEPTH`) beyond a taper's edge: how many, and one's centre. The
    /// study solved thin material there stiffer than it is and read it in
    /// part, and thin air's sides moving as one, so a claim on it is
    /// undecided.
    #[serde(default)]
    pub thin_left: Option<(usize, V3)>,
    /// Cells left not of one material when splitting stopped (`PROPERTY`),
    /// beyond where the materials they hold are resolved: how many, and
    /// one's centre. The study solved a material feature thinner than the
    /// cells (a layer, a strip) as one coefficient per cell, lost or
    /// filling the cell, so a claim on it is undecided.
    #[serde(default)]
    pub material_left: Option<(usize, V3)>,
    /// Full cells holding an interface between materials left whole
    /// (`material_held`), parallel to a wall of the cell and between its
    /// centre and that wall, where the wall lies on a plane of the grid
    /// twice as coarse (`interface_shared_what`): how many, one's
    /// centre, and that grid's cell (mm). A cell takes its centre's
    /// value, so this grid and the one twice as coarse both place the
    /// face on that plane, and a claim's ladder on them cannot see the
    /// error of its place: a claim on it is undecided.
    #[serde(default)]
    pub interface_shared: Option<(usize, V3, f64)>,
    /// Full cells holding an interface between two materials that do not
    /// lie in them as layers along one axis, alike across their slab of
    /// cells (a face oblique to the cells, or curved, or its edge): how
    /// many, and one's centre. Each cell taking its centre's value, the
    /// cells along such a face are a stair of one material or the other on
    /// every grid, and the grids' agreement does not show its error: a
    /// claim on the study is undecided (`interface_stepped_what`).
    #[serde(default)]
    pub interface_stepped: Option<(usize, V3)>,
    /// Each integral over the surface in a region (a pressure, a
    /// convection, a frictionless hold): the region's surface its facets
    /// covered against the same measure by lines (`SurfaceMeasure`).
    #[serde(default)]
    pub surface_measures: Vec<SurfaceMeasure>,
    /// Each load's reach along each axis (`spread_of` the cells it lands in,
    /// by their share of it): a load less than a cell across is concentrated.
    #[serde(default)]
    pub load_spread: Vec<(String, V3)>,
    /// Where each load landed (`Footprint`), indexed as `applied` (gravity's
    /// first): none where it landed nowhere.
    #[serde(default)]
    pub load_footprints: Vec<Option<Footprint>>,
    /// Where each fix holds: the part's surface in its region, and the
    /// material it holds; indexed as `reactions`, none for a cut of symmetry
    /// (its plane is its place) or a region the part's surface misses.
    #[serde(default)]
    pub fix_footprints: Vec<Option<Footprint>>,
    /// The centroid of the part's material as its cells hold it (mm): each
    /// cell's at its centre.
    #[serde(default)]
    pub centroid: V3,
    /// The memory the solve was estimated to take before it allocated its
    /// cells (bytes, `memory`).
    #[serde(default)]
    pub memory_estimate: f64,
}

// This coverage cutoff alone cannot bound error in a local response.
pub(crate) const PRESSURE_COVERAGE_LIMIT: f64 = 0.03;

/// A region's surface as a study's integral over it covered, checked:
/// both measures are sum |n_a| dA over the part's surface in the region
/// (mm^2: a face's area, times up to the root of 3 where it lies across
/// the axes), the integral's by its facets, the check's by counting where
/// lines through the cells cross the surface between the mesh's planes
/// (`surface_measure`). They are compared over each patch (`SPLIT_BLOCK`
/// grid cells to a side) with those beside it, and over the whole region,
/// each against the lines' noise (three times the root of their count's
/// variance) and rounding at the region's rims: a feature's surface the
/// mesh missed shows round its patch however large the region.
#[derive(Default, Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SurfaceMeasure {
    /// The load's or fix's name, and which it is (`load`, `fix`, `convection`).
    pub name: String,
    pub kind: String,
    pub covered: f64,
    pub lines: f64,
    /// How far the lines' measure may be off (noise and rim rounding).
    pub uncertainty: f64,
    /// The surface the integral missed: by how much the lines' measure
    /// passes the facets' beyond its uncertainty, over the region or over a
    /// patch and those beside it, whichever is more (0: none).
    pub shortfall: f64,
    /// Pressure applied beyond the region, past the same uncertainty.
    #[serde(default)]
    pub excess: f64,
    /// An unloaded pressure region or an unprestressed modal study.
    #[serde(default)]
    pub(crate) unused_pressure: bool,
    /// How far the region reaches along each axis: the standard deviation
    /// of its facets' places, weighted by their share (a point's is 0, a
    /// face's its width over the root of 12).
    #[serde(default)]
    pub spread: V3,
}

/// Where a load or a fix landed, as the solve took it: what a reader checks
/// a setup by (a force on the wrong end, a clamp where the bolt is not).
#[derive(Default, Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Footprint {
    /// The box of what it landed on (mm): the surface points it took, or
    /// the cells its material lies in (then to within a cell), within the
    /// part's bounds.
    pub lo: V3,
    pub hi: V3,
    /// Its centre (mm): a load's, where its force acts (its material's
    /// centroid, a pressure's centre of force — for a flat face, its moment
    /// about the origin is `at` x its total); a fix's, its surface's
    /// centroid (else its held material's).
    pub at: V3,
    /// The material it acts on or holds (mm^3; 0 for a surface's).
    pub volume: f64,
    /// The part's surface it acts on or holds (mm^2; 0 for a material's).
    pub area: f64,
    /// A pressure's least, mean (by area) and greatest value over it (MPa).
    pub pressure: Option<[f64; 3]>,
}

/// A footprint as it is gathered: its box, its weighted centre, its
/// measures and a pressure's range.
pub(crate) struct Gather {
    lo: V3,
    hi: V3,
    w: f64,
    wp: V3,
    pub(crate) volume: f64,
    pub(crate) area: f64,
    /// (least, greatest, its integral over the area)
    pressure: Option<[f64; 3]>,
}

impl Gather {
    pub(crate) fn new() -> Gather {
        Gather { lo: [f64::INFINITY; 3], hi: [f64::NEG_INFINITY; 3], w: 0.0, wp: [0.0; 3], volume: 0.0, area: 0.0, pressure: None }
    }
    /// The box takes in `p`.
    pub(crate) fn reach(&mut self, p: V3) {
        for a in 0..3 {
            self.lo[a] = self.lo[a].min(p[a]);
            self.hi[a] = self.hi[a].max(p[a]);
        }
    }
    /// A cell's box (`centre`, `size`).
    pub(crate) fn cell(&mut self, centre: V3, size: f64) {
        self.reach(centre.map(|c| c - 0.5 * size));
        self.reach(centre.map(|c| c + 0.5 * size));
    }
    /// `p` weighs `w` in the centre.
    pub(crate) fn weigh(&mut self, p: V3, w: f64) {
        self.w += w;
        for a in 0..3 {
            self.wp[a] += w * p[a];
        }
    }
    /// A pressure's sample: `pressure` (MPa) on `patch` (mm^2) at `p`, its force's magnitude `force` (N).
    pub(crate) fn press(&mut self, p: V3, pressure: f64, patch: f64, force: f64) {
        self.reach(p);
        self.weigh(p, force);
        self.area += patch;
        let [lo, hi, sum] = self.pressure.get_or_insert([f64::INFINITY, f64::NEG_INFINITY, 0.0]);
        (*lo, *hi, *sum) = (lo.min(pressure), hi.max(pressure), *sum + pressure * patch);
    }
    /// What it gathered, within the part's bounds `(lo, hi)`; none where nothing landed.
    pub(crate) fn footprint(self, (part_lo, part_hi): (V3, V3)) -> Option<Footprint> {
        if (0..3).any(|a| !(self.lo[a] <= self.hi[a] && part_lo[a] <= part_hi[a])) {
            return None;
        }
        let lo: V3 = std::array::from_fn(|a| self.lo[a].clamp(part_lo[a], part_hi[a]));
        let hi: V3 = std::array::from_fn(|a| self.hi[a].clamp(part_lo[a], part_hi[a]));
        let at = if self.w > 0.0 { self.wp.map(|v| v / self.w) } else { std::array::from_fn(|a| 0.5 * (lo[a] + hi[a])) };
        let pressure = self.pressure.map(|[l, h, sum]| [l, if self.area > 0.0 { sum / self.area } else { 0.0 }, h]);
        Some(Footprint { lo, hi, at, volume: self.volume, area: self.area, pressure })
    }
}

/// The standard deviation of `pts`' places along each axis, weighted.
pub fn spread_of(pts: impl Iterator<Item = (V3, f64)>) -> V3 {
    let (mut w, mut m, mut q) = (0.0, [0.0; 3], [0.0; 3]);
    for (p, a) in pts {
        w += a;
        for k in 0..3 {
            m[k] += a * p[k];
            q[k] += a * p[k] * p[k];
        }
    }
    if w <= 0.0 {
        return [0.0; 3];
    }
    std::array::from_fn(|k| (q[k] / w - (m[k] / w).powi(2)).max(0.0).sqrt())
}

impl SurfaceMeasure {
    /// `covered`: the facets the integral took, each with its share of the
    /// region.
    pub(crate) fn of(name: &str, kind: &str, covered: &[(Facet, f64)], lines: Vec<([i64; 3], Patch)>, grid: &Grid) -> SurfaceMeasure {
        let h = grid.h;
        let patch = |p: V3| -> [i64; 3] { std::array::from_fn(|a| ((p[a] - grid.origin[a]) / h).floor() as i64).map(|x| x.div_euclid(SPLIT_BLOCK)) };
        // A point of a facet's patch: the one holding the material just
        // under it (as `part_surface` gives a facet its cell, and
        // `surface_measure` a crossing its box).
        let under = |p: V3, f: &Facet| patch(std::array::from_fn(|a| p[a] - 1e-6 * h * f.normal[a]));
        // Each facet's measure in each patch it lies in.
        let mut by: Vec<([i64; 3], f64)> = covered
            .par_iter()
            .flat_map_iter(|(f, w)| {
                let ps: Vec<[i64; 3]> = f.corners.iter().map(|c| under(*c, f)).collect();
                if ps.iter().all(|p| *p == ps[0]) {
                    return vec![(ps[0], f.axial() * w)];
                }
                let (a, b): ([i64; 3], [i64; 3]) =
                    (std::array::from_fn(|k| ps.iter().map(|p| p[k]).min().unwrap_or(0)), std::array::from_fn(|k| ps.iter().map(|p| p[k]).max().unwrap_or(0)));
                let mut out = Vec::new();
                for x in a[0]..=b[0] {
                    for y in a[1]..=b[1] {
                        for z in a[2]..=b[2] {
                            let q = [x, y, z];
                            let (lo, hi) = (grid.pos(q.map(|v| v * SPLIT_BLOCK)), grid.pos(q.map(|v| (v + 1) * SPLIT_BLOCK)));
                            out.extend(f.clipped(lo, hi).iter().filter(|piece| under(piece.at, piece) == q).map(|piece| (q, piece.axial() * w)));
                        }
                    }
                }
                out
            })
            .collect();
        by.sort_unstable_by_key(|b| b.0);
        let mut facets: Vec<([i64; 3], f64)> = Vec::new();
        for (k, v) in by {
            match facets.last_mut() {
                Some(last) if last.0 == k => last.1 += v,
                _ => facets.push((k, v)),
            }
        }
        // Per patch, what the lines count past the facets (less where the
        // facets cover more), its variance and its rims' rounding.
        let mut per: std::collections::HashMap<[i64; 3], (f64, f64, f64)> = std::collections::HashMap::new();
        for (k, c) in &facets {
            per.entry(*k).or_default().0 -= c;
        }
        for (k, p) in &lines {
            let e = per.entry(*k).or_default();
            e.0 += p.lines;
            e.1 += p.var;
            e.2 += p.rim;
        }
        // Each patch short past its own uncertainty, with those beside it
        // (27): where the mesh crosses a patch's face off the surface, what
        // one patch lacks the next holds, and a feature missed has nothing
        // beside it to make up.
        let local = |direction: f64| {
            per.iter()
                .filter(|(_, e)| direction * e.0 > e.2 + 3.0 * e.1.sqrt())
                .map(|(k, _)| {
                    let (mut ex, mut var, mut rim) = (0.0, 0.0, 0.0);
                    for d in 0..27 {
                        if let Some(e) = per.get(&[k[0] + d % 3 - 1, k[1] + d / 3 % 3 - 1, k[2] + d / 9 - 1]) {
                            ex += e.0;
                            var += e.1;
                            rim += e.2;
                        }
                    }
                    direction * ex - rim - 3.0 * var.sqrt()
                })
                .fold(0.0, f64::max)
        };
        let (var, rim): (f64, f64) = lines.iter().fold((0.0, 0.0), |a, l| (a.0 + l.1.var, a.1 + l.1.rim));
        let total: f64 = facets.iter().map(|f| f.1).sum();
        let all: f64 = lines.iter().map(|l| l.1.lines).sum();
        let uncertainty = rim + 3.0 * var.sqrt();
        let shortfall = (all - total - uncertainty).max(local(1.0));
        let excess = if kind == "load" { (total - all - uncertainty).max(local(-1.0)) } else { 0.0 };
        let spread = spread_of(covered.iter().map(|(f, w)| (f.at, f.area * w)));
        SurfaceMeasure { name: name.into(), kind: kind.into(), covered: total, lines: all, uncertainty, shortfall, excess, unused_pressure: false, spread }
    }

    /// What a claim on a study whose integral this is cannot see, where it
    /// was short: the part loaded, cooled, heated or held by less of its
    /// region than written, by how much its claims cannot say (a load
    /// missed on a lever moves a tip far more than its share of the
    /// region, a film's missing heat either way), so they are unknown,
    /// failing as well as holding, and the claim's coarser grids, which
    /// say so, leave it to finer ones.
    pub fn unseen(&self) -> Option<String> {
        if self.unused_pressure {
            return None;
        }
        let discrepancy = if self.excess > self.shortfall {
            format!("covered {:.4} mm2 beyond the pressure region", self.excess)
        } else if self.shortfall > 0.0 {
            format!("missed {:.4} mm2 of the region", self.shortfall)
        } else {
            return None;
        };
        Some(format!(
            "{} '{}' acts on {:.4} mm2 of the surface in its region (measured along the axes), where lines through the cells cross {:.4} mm2 \
             of it, to {:.4}: its integral over the surface {discrepancy}. Refine (res r) across the region's rim",
            self.kind, self.name, self.covered, self.lines, self.uncertainty
        ))
    }
}

/// A static elastic study's linear buckling analysis: the load factors
/// lambda of K phi = -lambda K_G phi, K_G the geometric stiffness of the
/// solved stress (each cell's, at its material's centroid): the study's
/// loads times lambda buckle the part, by small-displacement theory.
#[derive(Default, Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Buckle {
    /// Axes whose reflection changes the sign of this displacement mode.
    #[serde(default)]
    pub antisymmetric: [bool; 3],
    /// The lowest factor as found: the Rayleigh quotient of its mode, so
    /// at least the lowest factor. Its mode is the solution's mode 1.
    /// Infinite where no direction it found is compressed by the loads (a
    /// bar pulled): no mode found, its band `low` and up.
    pub factor: f64,
    /// The least the lowest factor may be, by the eigensolve's residual
    /// (`buckle`): the factor's band is `low` to `factor`. Where no mode
    /// was found, the largest factor its search excluded (infinite where
    /// it excluded all).
    #[serde(default)]
    pub low: f64,
    /// The first mode's largest displacement: where it buckles.
    pub at: V3,
    pub iterations: usize,
    pub ms: f64,
    /// Why it gave no factor: the eigensolver did not resolve one (it did
    /// not converge, or found no negative direction to resolve), or broke
    /// down.
    pub failure: Option<String>,
}

impl Buckle {
    /// The lowest factor's band, low to high, if the analysis gave one.
    pub fn answer(&self) -> Option<(f64, f64)> {
        self.failure.is_none().then_some((self.low, self.factor))
    }
}

/// The cells a solve's nodes are of: its grid's, and how fine it split
/// them (`MAX_DEPTH`).
pub fn cells_of(h: f64, depth: u8) -> String {
    if depth == 0 { format!("cells of {h} mm") } else { format!("cells of {h} mm split to {} mm", h / f64::from(1u32 << depth)) }
}

/// The warning of thin features left when splitting stopped (of `cells_of`).
pub fn thin_left_warning(cells: usize, at: V3, of: &str) -> String {
    format!(
        "a feature thinner than {THIN} of its cells is left: {cells} {of}, one centred at {at:?}, hold material or air thinner than {THIN} \
         cells across beyond a taper's edge (a fin, a web or a ligament, which they hold too stiff and read only in part; a slit or a \
         clearance, whose sides they move as one). Refine (res r) to cells under a quarter of its thickness"
    )
}

/// The warning of material features left when splitting stopped (of
/// `cells_of`).
pub fn material_left_warning(cells: usize, at: V3, of: &str) -> String {
    format!("{}. Refine (res r) to cells under a quarter of its thickness", material_left_what(cells, at, of))
}

/// What an interface on a plane a claim's grids share
/// (`Report::interface_shared`) is.
pub fn interface_shared_what(cells: usize, at: V3, plane: f64, of: &str) -> String {
    format!(
        "an interface between materials less than half a cell off a plane of cells of {plane} mm: {cells} {of}, one centred at {at:?}, \
         hold a face parallel to that plane between their centres and it; a cell takes its centre's value, so these cells and those twice \
         as coarse both place the face on the plane, and their agreement does not show the error of its place (a band of another \
         stiffness or conductivity read as if it ran to the plane)"
    )
}

/// What an interface not in layers along the cells
/// (`Report::interface_stepped`) is.
pub fn interface_stepped_what(cells: usize, at: V3, of: &str) -> String {
    format!(
        "an interface between materials that is not flat along the cells: {cells} {of}, one centred at {at:?}, hold two materials not as \
         layers along one axis, alike across their row (a face oblique to the cells, or curved, or its edge); a cell takes its centre's \
         value, so the cells along the face are a stair of one material or the other on any cells, and their agreement does not show its \
         error (a soft band's face turned a degree read 36 times its stress; an insulating face turned 2, 9 K low). Turn the part so the \
         face lies along an axis"
    )
}

/// What material left at the floor (`Report::material_left`) is.
pub fn material_left_what(cells: usize, at: V3, of: &str) -> String {
    format!(
        "a material thinner than {THIN} of its cells is left: {cells} {of}, one centred at {at:?}, are not of one material \
         (a property varies across them) and what they hold is held {THIN} cells across nowhere within {THIN} of them (a layer or a \
         strip of another conductivity, stiffness, density or yield), so each takes one value for it: the feature lost, or filling the cell"
    )
}

/// The warning of a gap the cells glue shut.
pub fn glued_warning(cells: usize, at: V3, width: f64, of: &str) -> String {
    let what = if width > 0.0 {
        format!(
            "the narrowest about {width} mm wide, at {at:?}. The study solved the part with the gap filled: its answers are the filled part's. \
             Cells of {} mm or finer part it, split to an eighth (a claim's, checked on cells four times as coarse: (res {})); or close the gap \
             in the model",
            4.0 * width,
            width
        )
    } else {
        format!(
            "at {at:?} one whose sides bounds do not show joined, though no air was found between them (a hairline, or material too thin \
             for them to show solid). The study solved the part joined there. Join the pieces with overlapping material, or open the gap"
        )
    };
    format!(
        "a gap is glued shut: {cells} {of} hold material on both sides of air too thin for their samples to part, which the part \
         is not shown to join round within the cells about them (a slit, a saw cut, a clearance between parts); {what}"
    )
}
