//! Points sampled on the exact surface (not a mesh): the surface cells of an
//! interval octree, projected onto the zero set by Newton steps, with their
//! normals and areas; the names that can form the surface; wall thickness
//! by rays.

use crate::octree::{Cell, Grid, Local, descend_padded, root_for};
use crate::query::{Side, ray_to_within};
use crate::stop::Stopped;
use ink_core::V3;
use ink_core::eval::{Ctx, decided_names, eval_grads, eval_interval, names_on_paths};
use ink_core::sign::inside;
use ink_core::tape::{BOTH, Tape};
use ink_core::v3::{add, cross, direction, dot, norm, scale, sub, unit};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

/// How close to the zero set surface() puts its points, for a region of
/// largest extent `ext`: within this distance by the field's local slope.
pub fn surface_tolerance(ext: f64) -> f64 {
    1e-9 * ext
}

/// Relative tolerance on whether a plane crosses a grid cell, against
/// rounding: surface() and plane_area_in_cell() share it.
const CELL_TOL: f64 = 1e-7;

/// Points on the surface of `shape` inside `lo..hi`: one per grid cell the
/// surface crosses (about `cells` across the largest extent), projected onto
/// the zero set, with unit outward normal and area weight (sums to the
/// area). Also returns the grid step.
///
/// Work is done per block of 4^3 cells with the tape the octree simplified
/// for the block grown by one cell, so big unions cost what the nearby
/// parts cost; the Newton steps stay within that grown block.
pub fn surface(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, cells: f64) -> (Vec<(V3, V3, f64)>, f64) {
    let (blocks, h) = surface_blocks(shape, ctx, lo, hi, cells);
    (blocks.into_iter().flat_map(|b| b.0).collect(), h)
}

/// `surface`, by block of 4^3 cells, each with the length of the tape the
/// octree simplified for it: what evaluating the field there costs.
fn surface_blocks(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, cells: f64) -> (Vec<(Vec<(V3, V3, f64)>, usize)>, f64) {
    const B: i64 = 4;
    let ext = crate::query::max_extent(lo, hi);
    let h = ext / cells;
    let (pad, padhi) = ([lo[0] - h, lo[1] - h, lo[2] - h], [hi[0] + h, hi[1] + h, hi[2] + h]);
    // Only cells meeting lo..hi, padded: a field not clipped to it along an
    // axis it does not read (a profile's z) is sampled there alone.
    let (grid, top) = root_for(pad, padhi, h, B);
    let grid = crate::octree::Grid { clip: (pad, padhi), ..grid };
    let blocks: Vec<(Vec<(V3, V3, f64)>, usize)> = descend_padded(
        shape,
        ctx,
        &grid,
        top,
        B,
        1.0,
        &|t: &Arc<Tape>, blk| {
            let (blo, bhi) = grid.cell_bounds(blk);
            let mut boxes = Vec::with_capacity(64);
            for z in 0..B {
                for y in 0..B {
                    for x in 0..B {
                        boxes.push(grid.cell_bounds(Cell { i: [blk.i[0] + x, blk.i[1] + y, blk.i[2] + z], size: 1 }));
                    }
                }
            }
            let mut p: Vec<V3> = boxes.iter().map(|(a, b)| scale(add(*a, *b), 0.5)).collect();
            // Newton steps onto the zero set, up to 8, until no point that
            // can move is further from it than a thousandth of the
            // tolerance (half of a model's blocks take one step).
            let mut g = eval_grads(t, ctx, &p);
            for _ in 0..8 {
                let mut moved = false;
                for (q, g) in p.iter_mut().zip(&g) {
                    let gg = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
                    if gg > 0.0 && g[0].is_finite() && g[0].abs() > 1e-3 * surface_tolerance(ext) * gg.sqrt() {
                        *q = sub(*q, scale([g[1], g[2], g[3]], g[0] / gg));
                        moved = true;
                    }
                }
                if !moved {
                    break;
                }
                g = eval_grads(t, ctx, &p);
            }
            // The field at the block's cell corners.
            let n1 = (B + 1) as usize;
            let corners: Vec<V3> = (0..n1 * n1 * n1)
                .map(|i| {
                    let q = ink_core::grid::ijk(i, [n1; 3]);
                    std::array::from_fn(|k| blo[k] + h * q[k] as f64)
                })
                .collect();
            let mut fc = Vec::new();
            ink_core::eval::eval_points(t, ctx, &corners, &mut fc);
            let mut out = Vec::new();
            for (ci, ((q, g), (a, b))) in p.iter().zip(&g).zip(&boxes).enumerate() {
                // A cell whose corners are all outside, meeting the surface
                // only on its boundary (beside a face's edge), holds none of
                // it: its tangent plane would add a row of area the
                // neighbours already have.
                let [cx, cy, cz] = ink_core::grid::ijk(ci, [B as usize; 3]);
                let outside = (0..8).all(|k| !inside(fc[ink_core::grid::index([cx + (k & 1), cy + ((k >> 1) & 1), cz + ((k >> 2) & 1)], [n1; 3])]));
                let strictly_in = (0..3).all(|k| q[k] > a[k] + 1e-9 * h && q[k] < b[k] - 1e-9 * h);
                if outside && !strictly_in {
                    continue;
                }
                // The block's tape is valid within one cell of the block.
                if (0..3).any(|k| q[k] < blo[k] - h || q[k] > bhi[k] + h) {
                    continue;
                }
                // On the zero set by the local slope, where the field has a
                // normal: none where its gradient is zero or not finite (a
                // vertical tangent, by which any value is within tolerance).
                let grad = [g[1], g[2], g[3]];
                let Some(n) = direction(grad) else { continue };
                if !((g[0] / norm(grad)).abs() < surface_tolerance(ext)) {
                    continue;
                }
                // Keep the cell if the surface (locally its tangent plane)
                // actually crosses it: one sample per crossed cell. Half-open,
                // so a surface lying exactly on a cell boundary is counted
                // once, by the cell on the material side (as the corner test
                // above), and tolerant, so rounding cannot let both keep it.
                let d = sub(scale(add(*a, *b), 0.5), *q);
                let (s, reach) = (dot(d, n), 0.5 * h * (n[0].abs() + n[1].abs() + n[2].abs()));
                if s >= reach * (1.0 - CELL_TOL) || s < -reach * (1.0 + CELL_TOL) || norm(d) > h {
                    continue;
                }
                // Its area: the tangent plane clipped to the cell.
                out.push((*q, n, plane_area_in_cell(scale(add(*a, *b), 0.5), h, *q, n)));
            }
            Some((out, t.len()))
        },
        &|_, _| None,
    );
    (blocks, h)
}

/// A shape's surface seen coarsely, for saying what a question about it
/// will cost before it is asked (ink check): a survey's (Survey::sketch),
/// its area from the survey's samples, and the area-weighted mean length of the tape the
/// octree simplified where the surface is (the instructions an evaluation
/// there runs). One query, cached.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Sketch {
    pub area: f64,
    /// Instructions per evaluation near the surface.
    pub ops: f64,
}

impl Sketch {
    /// Of surface samples by block, each block with its tape's length.
    fn of(blocks: &[(Vec<(V3, V3, f64)>, usize)], whole: usize) -> Sketch {
        let (mut weighed, mut area) = (0.0, 0.0);
        for (b, len) in blocks {
            for (_, _, a) in b {
                (weighed, area) = (weighed + a * *len as f64, area + a);
            }
        }
        Sketch { area, ops: if area > 0.0 { weighed / area } else { whole as f64 } }
    }
}

/// What a model's lint needs to know of a shape's surface (`ink check`),
/// from about SURVEYED samples of it (`samples`): one query, so it is
/// cached whole (cache::memo_region).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Survey {
    /// The samples' spacing.
    pub h: f64,
    /// The names (tags) seen forming the surface, and those on a path it
    /// can take that none was seen to, with in how many cells, and in how
    /// many of those it was sought finer (`surface_names`), if asked.
    pub formed: Vec<u32>,
    pub unseen: Vec<(u32, usize, usize)>,
    /// Where interval bounds are much wider than the field varies
    /// (`Loose`), if they are.
    pub loose: Option<Loose>,
    /// The undefined cases at the samples and on a grid over the box
    /// (`undefined`), and at how many points.
    pub undefined: Vec<Case>,
    pub points: usize,
    /// Material lies just inside the surface at some sample: false for a
    /// shape with no interior that intervals cannot prove empty (a sheet of
    /// no thickness, as where two shapes that only touch are intersected).
    pub interior: bool,
    /// Its area and the instructions an evaluation near it runs, from its
    /// samples: what a question about it will cost (ink check).
    pub sketch: Sketch,
}

/// Interval widths over small boxes (a quarter of a cell across, of the
/// `cells` a survey is given) at up to 2000 samples, over the width the
/// gradient implies or the field's range over the box's corners and
/// centre, the larger:
/// where these ratios are large, proofs and meshing refine far deeper than
/// the geometry needs (branch cuts of atan2 or mod, products of selects,
/// repeated variables).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Loose {
    /// The ratio exceeded by a tenth of the samples (over 5).
    pub p90: f64,
    /// The share of samples whose ratio is over 5.
    pub share: f64,
    /// The largest ratio, and where.
    pub max: f64,
    pub at: V3,
}

/// An undefined case (eval::undefined_cases): its kind, the construction
/// it is in, how many points hit it, and the first.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Case {
    pub what: String,
    pub name: Option<u32>,
    pub count: usize,
    pub at: V3,
}

/// The surface samples a `Survey` reads: up to this many for undefined
/// cases, half as many for looseness, the first block for an interior.
const SURVEYED: usize = 4000;

/// The shape's `Survey` of `lo..hi`, which must hold its material: about
/// SURVEYED samples, no finer than `cells` across, and its looseness over
/// boxes a quarter of a cell `cells` across; with the names that can form
/// its surface, proven at `cells` across, if `names`.
pub fn survey(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, cells: f64, names: bool) -> Survey {
    crate::cache::memo_region(
        "survey",
        &[shape],
        ctx,
        (lo, hi),
        |fp| survey_inputs(shape, ctx, cells, names, fp),
        |ctx| {
            let (blocks, h) = samples(shape, ctx, lo, hi, cells);
            let sketch = Sketch::of(&blocks, shape.len());
            let samples: Vec<(V3, V3, f64)> = blocks.into_iter().flat_map(|b| b.0).collect();
            let interior = interior(shape, ctx, &samples, h);
            let pts: Vec<V3> = samples.into_iter().map(|p| p.0).collect();
            let (undefined, points) = undefined(shape, ctx, &pts, lo, hi);
            let (formed, unseen) = if names { surface_names(shape, ctx, lo, hi, cells) } else { Default::default() };
            let loose = loose(shape, ctx, &pts, crate::query::max_extent(lo, hi) / cells);
            Survey { h, formed, unseen, loose, undefined, points, interior, sketch }
        },
    )
}

/// What a survey's answer depends on besides its tape, context and region:
/// its names and cases give tags by the model's numbers, as the context it
/// reads (`Ctx::project`) numbers them.
fn survey_inputs(shape: &Tape, ctx: &Ctx, cells: f64, names: bool, fp: &mut ink_core::graph::Fingerprint) {
    ink_core::tape::Names::of(&[shape], &ctx.project(&[shape])).raw(fp.f64(cells).word(names as u64));
}

/// `survey`'s answer where it is cached already (ink check ran it), with
/// nothing computed where it is not.
pub fn survey_cached(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, cells: f64, names: bool) -> Option<Survey> {
    crate::cache::cached_region("survey", &[shape], ctx, (lo, hi), |fp| survey_inputs(shape, ctx, cells, names, fp))
}

/// About SURVEYED samples of the surface, one per cell crossed (`surface`),
/// no finer than `cells` across: the cells a surface crosses go as the
/// square of the cells across, so a pass 16 across says how fine to go.
fn samples(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, cells: f64) -> (Vec<(Vec<(V3, V3, f64)>, usize)>, f64) {
    let first = surface_blocks(shape, ctx, lo, hi, cells.min(16.0));
    let n: usize = first.0.iter().map(|b| b.0.len()).sum();
    if n >= SURVEYED || cells <= 16.0 {
        return first;
    }
    surface_blocks(shape, ctx, lo, hi, (16.0 * (SURVEYED as f64 / n.max(1) as f64).sqrt()).min(cells))
}

/// Whether the field is negative a thousandth of the spacing `h` inside
/// any surface sample (point, outward normal): a solid's are, at once (so
/// this stops at the first block that has one).
fn interior(shape: &Tape, ctx: &Ctx, samples: &[(V3, V3, f64)], h: f64) -> bool {
    let mut f = Vec::new();
    samples.chunks(1024).any(|c| {
        let inner: Vec<V3> = c.iter().map(|(p, n, _)| sub(*p, scale(*n, 1e-3 * h))).collect();
        ink_core::eval::eval_points(shape, ctx, &inner, &mut f);
        f.iter().any(|&v| inside(v))
    })
}

/// The names (tags) that form the surface of `shape` in `lo..hi`, and those
/// that may: (formed, unseen). A name that forms surface at a point is on
/// the path the value takes there, so on the paths (eval::names_on_paths)
/// of the tape interval arithmetic simplifies for each cell, `cells`
/// across, the surface can cross: a name on none is proven to form none.
/// One on a path is formed where it is seen deciding a patch of surface:
/// at a point projected onto it in such a cell, and at four points across
/// the surface around it (eval::decided_names at each), so not only at an
/// edge, where it ties with what covers it, nor as a side a blend weighs
/// out. A name on a path not seen at 27 points in any cell is then sought
/// so in the eighths (still on a path) of up to SOUGHT of those cells,
/// spread evenly over them, down to an eighth of a cell. Unseen is a
/// sample's answer, not a bound: with how many cells it was on a path in,
/// and in how many it was searched finer.
pub fn surface_names(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, cells: f64) -> (Vec<u32>, Vec<(u32, usize, usize)>) {
    let ext = crate::query::max_extent(lo, hi);
    let h = ext / cells / FINER as f64;
    let (grid, top) = root_for(sub(lo, [h * FINER as f64; 3]), add(hi, [h * FINER as f64; 3]), h, 1);
    let found = Found { tol: surface_tolerance(ext), names: Mutex::new((HashMap::new(), HashSet::new())), later: Mutex::new(Vec::new()) };
    names_below(shape, ctx, &grid, top, &found);
    // Up to SOUGHT cells for each name still unseen, spread over those it
    // is on a path in (in a fixed order).
    let mut later = std::mem::take(&mut *found.later.lock().unwrap());
    later.sort_by_key(|(c, _)| c.i);
    let (on, seen) = found.names.lock().unwrap().clone();
    let paths: Vec<Vec<u32>> = later.iter().map(|(_, sub)| names_on_paths(sub)).collect();
    let mut pick = vec![false; later.len()];
    let unseen_names: Vec<u32> = on.keys().filter(|t| !seen.contains(t)).copied().collect();
    for t in &unseen_names {
        let mine: Vec<usize> = (0..later.len()).filter(|&k| paths[k].contains(t)).collect();
        let m = mine.len().min(SOUGHT);
        (0..m).for_each(|k| pick[mine[k * mine.len() / m]] = true);
    }
    // How many cells each was sought finer in: those picked on its paths,
    // for it or for another name.
    let sought: HashMap<u32, usize> = unseen_names.iter().map(|t| (*t, (0..later.len()).filter(|&k| pick[k] && paths[k].contains(t)).count())).collect();
    let later: Vec<_> = later.into_iter().zip(pick).filter_map(|(c, p)| p.then_some(c)).collect();
    later.par_iter().for_each(|(cell, sub)| cell.children().par_iter().for_each(|k| names_below(sub, ctx, &grid, *k, &found)));
    let (on, seen) = found.names.into_inner().unwrap();
    let sorted = |s: HashSet<u32>| {
        let mut v: Vec<u32> = s.into_iter().collect();
        v.sort_unstable();
        v
    };
    let mut unseen: Vec<(u32, usize, usize)> = on.into_iter().filter(|(t, _)| !seen.contains(t)).map(|(t, n)| (t, n, sought.get(&t).copied().unwrap_or(0))).collect();
    unseen.sort_unstable();
    (sorted(seen), unseen)
}

/// The names (tags) that form the surface of a moving shape over its
/// motion `span`, and those that may: surface_names at moments of the
/// span, `first` the answer at its start, then the span halved again and
/// again (1/2, 1/4, 3/4, 1/8, ...: up to MOMENTS), stopping once each of
/// `names` that can form surface at some moment is seen at one. `at` gives a
/// moment's context and a box holding its material (None: it is not looked
/// at). Those that can are the names on a path over `over`, a context with
/// t ranged over the span and a box holding its material at every moment
/// (interval arithmetic bounds every moment at once): a name on none is
/// proven to form none at any moment. Without `over`, every name may.
/// (formed, unseen as surface_names gives it, with 0 cells for a name on
/// no path at a moment looked at, the moments looked at.)
pub fn surface_names_over(
    shape: &Arc<Tape>,
    cells: f64,
    span: (f64, f64),
    names: &[u32],
    first: (Vec<u32>, Vec<(u32, usize, usize)>),
    mut at: impl FnMut(f64) -> Option<(Ctx, V3, V3)>,
    over: Option<(&Ctx, V3, V3)>,
) -> (Vec<u32>, Vec<(u32, usize, usize)>, usize) {
    let (mut formed, unseen) = first;
    let mut found: HashMap<u32, (usize, usize)> = unseen.into_iter().map(|(t, n, f)| (t, (n, f))).collect();
    let asked: Vec<u32> = names.iter().filter(|t| !formed.contains(t)).copied().collect();
    let mut can: HashSet<u32> = match over {
        Some((ctx, lo, hi)) => memo_names("names-over", shape, ctx, (lo, hi), cells, &asked, |c| names_over(shape, c, lo, hi, cells, &asked)),
        None => names_on_paths(shape).into_iter().collect(),
    };
    can.extend(found.keys());
    can.retain(|t| asked.contains(t));
    let mut looked = 1;
    let len = span.1 - span.0;
    let halvings = (1..MOMENTS.ilog2() + 1).flat_map(|d| (1..1u32 << d).step_by(2).map(move |k| k as f64 / (1u32 << d) as f64));
    for t in halvings.map(|f| span.0 + f * len) {
        if can.iter().all(|t| formed.contains(t)) {
            break;
        }
        let Some((ctx, lo, hi)) = at(t) else { continue };
        looked += 1;
        let (seen, unseen) = memo_names("surface-names", shape, &ctx, (lo, hi), cells, &[], |c| surface_names(shape, c, lo, hi, cells));
        formed.extend(seen.into_iter().filter(|t| !formed.contains(t)).collect::<Vec<_>>());
        for (t, n, f) in unseen {
            let e = found.entry(t).or_default();
            *e = (e.0.max(n), e.1.max(f));
        }
    }
    formed.sort_unstable();
    let mut unseen: Vec<(u32, usize, usize)> =
        can.into_iter().filter(|t| !formed.contains(t)).map(|t| (t, found.get(&t).map_or(0, |e| e.0), found.get(&t).map_or(0, |e| e.1))).collect();
    unseen.sort_unstable();
    (formed, unseen, looked)
}

/// How many moments surface_names_over looks at, at most: a power of 2.
const MOMENTS: u32 = 16;

/// A names query of `shape` over `lo..hi`, `cells` across, cached
/// (cache::memo_region): its answer gives tags by the model's numbers, as
/// `asked` does.
fn memo_names<T: serde::Serialize + serde::de::DeserializeOwned + std::fmt::Debug>(
    query: &str,
    shape: &Arc<Tape>,
    ctx: &Ctx,
    region: (V3, V3),
    cells: f64,
    asked: &[u32],
    compute: impl FnOnce(&Ctx) -> T,
) -> T {
    let inputs = |fp: &mut ink_core::graph::Fingerprint| {
        asked.iter().for_each(|&t| _ = fp.word(t as u64));
        survey_inputs(shape, ctx, cells, true, fp)
    };
    crate::cache::memo_region(query, &[shape], ctx, region, inputs, compute)
}

/// Those of `asked` on the paths of the tapes interval arithmetic
/// simplifies for the cells of `lo..hi`, `cells` across, the surface can
/// cross over `ctx`'s scope (surface_names' proof, over every moment of a
/// motion at once). A cell whose paths hold none not found yet is not
/// descended.
fn names_over(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, cells: f64, asked: &[u32]) -> HashSet<u32> {
    let h = crate::query::max_extent(lo, hi) / cells;
    let (grid, top) = root_for(sub(lo, [h; 3]), add(hi, [h; 3]), h, 1);
    let found = Mutex::new(HashSet::new());
    fn below(tape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, cell: Cell, asked: &[u32], found: &Mutex<HashSet<u32>>) {
        let mut ch = vec![BOTH; tape.len()];
        let r = eval_interval(tape, ctx, grid.cell_box(cell), Some(&mut ch));
        if r.none_inside() || r.all_inside() {
            return;
        }
        let sub = tape.narrowed(&ch);
        let names: Vec<u32> = names_on_paths(&sub).into_iter().filter(|t| asked.contains(t)).collect();
        {
            let mut f = found.lock().unwrap();
            if names.iter().all(|t| f.contains(t)) {
                return;
            }
            if cell.size <= 1 {
                f.extend(names);
                return;
            }
        }
        cell.children().par_iter().for_each(|k| below(&sub, ctx, grid, *k, asked, found));
    }
    below(shape, ctx, &grid, top, asked, &found);
    found.into_inner().unwrap()
}

/// How much finer than its cells surface_names searches for a name it has
/// not seen, and in how many of the cells it is on a path in, at most.
const FINER: i64 = 8;
const SOUGHT: usize = 64;

/// What surface_names has found: the names on a path (with how many cells
/// of FINER fine cells) and those seen, and the cells left with a name
/// unseen, with their tapes.
struct Found {
    tol: f64,
    names: Mutex<(HashMap<u32, usize>, HashSet<u32>)>,
    later: Mutex<Vec<(Cell, Arc<Tape>)>>,
}

/// Add to `found` the names on the paths of the tapes simplified for the
/// cells under `cell` the surface can cross (as `descend` simplifies them):
/// of cells of FINER fine cells, as on a path, and in every cell, those
/// seen. A cell whose own tape's paths hold no name not seen yet is not
/// descended: its cells' tapes are simplified from that one, so their
/// paths are among its own. A cell of FINER left with one is kept for the
/// search below it, a fine cell's is searched below at once.
fn names_below(tape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, cell: Cell, found: &Found) {
    let mut ch = vec![BOTH; tape.len()];
    let r = eval_interval(tape, ctx, grid.cell_box(cell), Some(&mut ch));
    if r.none_inside() || r.all_inside() {
        return;
    }
    let sub = tape.narrowed(&ch);
    let names = names_on_paths(&sub);
    {
        let mut f = found.names.lock().unwrap();
        if names.iter().all(|t| f.1.contains(t)) {
            return;
        }
        if cell.size == FINER {
            names.iter().for_each(|&t| *f.0.entry(t).or_default() += 1);
        }
    }
    if cell.size > FINER {
        cell.children().par_iter().for_each(|k| names_below(&sub, ctx, grid, *k, found));
    } else if !seen_in(&sub, ctx, grid, cell, &names, found) {
        if cell.size == FINER {
            found.later.lock().unwrap().push((cell, sub));
        } else if cell.size > 1 {
            cell.children().par_iter().for_each(|k| names_below(&sub, ctx, grid, *k, found));
        }
    }
}

/// Look for `names` (on the paths of `sub`, the cell's tape) forming
/// surface in the cell: points on the surface, each with four around it
/// across the surface, a thousandth of the cell away. Whether all of them
/// are seen.
fn seen_in(sub: &Tape, ctx: &Ctx, grid: &Grid, cell: Cell, names: &[u32], found: &Found) -> bool {
    let tol = found.tol;
    let (lo, hi) = grid.cell_bounds(cell);
    let h = hi[0] - lo[0];
    let starts: Vec<V3> = (0..27).map(|i| std::array::from_fn(|k| lo[k] + h * ((i / [1, 3, 9][k]) % 3) as f64 / 3.0 + h / 6.0)).collect();
    let e = 1e-3 * h;
    for (p, n) in onto(sub, ctx, starts, tol).into_iter().flatten() {
        if (0..3).any(|k| p[k] < lo[k] - 0.5 * h || p[k] > hi[k] + 0.5 * h) {
            continue;
        }
        // Surface only where the sign changes across it within the
        // tolerance (a zero the field keeps one sign around, faces meeting
        // flush, is none, a surface beside it notwithstanding): its names
        // are read there.
        let Some(p) = ink_core::sign::crossing(sub, ctx, p, n, e, tol) else {
            continue;
        };
        let mut here = decided_names(sub, ctx, p, tol);
        here.retain(|t| !found.names.lock().unwrap().1.contains(t));
        if here.is_empty() {
            continue;
        }
        let u = unit(cross(n, if n[0].abs() < 0.6 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] }));
        let v = cross(n, u);
        let around: Vec<V3> = [u, v, scale(u, -1.0), scale(v, -1.0)].iter().map(|d| add(p, scale(*d, e))).collect();
        for q in onto(sub, ctx, around, tol) {
            match q {
                Some((q, _)) => {
                    let there = decided_names(sub, ctx, q, tol);
                    here.retain(|t| there.contains(t));
                }
                None => here.clear(),
            }
        }
        let mut f = found.names.lock().unwrap();
        f.1.extend(here);
        if names.iter().all(|t| f.1.contains(t)) {
            return true;
        }
    }
    false
}

/// `pts` moved onto the zero set by Newton steps, up to 8: each that ends
/// within `tol` of it by the local slope, with its unit normal.
fn onto(t: &Tape, ctx: &Ctx, mut pts: Vec<V3>, tol: f64) -> Vec<Option<(V3, V3)>> {
    let mut g = eval_grads(t, ctx, &pts);
    for _ in 0..8 {
        let mut moved = false;
        for (q, g) in pts.iter_mut().zip(&g) {
            let gg = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
            if gg > 0.0 && g[0].is_finite() && g[0].abs() > 1e-3 * tol * gg.sqrt() {
                *q = sub(*q, scale([g[1], g[2], g[3]], g[0] / gg));
                moved = true;
            }
        }
        if !moved {
            break;
        }
        g = eval_grads(t, ctx, &pts);
    }
    pts.into_iter()
        .zip(g)
        .map(|(p, g)| {
            let n = direction([g[1], g[2], g[3]])?;
            ((g[0] / norm([g[1], g[2], g[3]])).abs() < tol).then_some((p, n))
        })
        .collect()
}

fn loose(tape: &Arc<Tape>, ctx: &Ctx, pts: &[V3], h: f64) -> Option<Loose> {
    use ink_core::interval::I;
    let step = (pts.len() / 2000).max(1);
    let sample: Vec<V3> = pts.iter().step_by(step).copied().collect();
    if sample.is_empty() {
        return None;
    }
    let grads = eval_grads(tape, ctx, &sample);
    // Small boxes: a coarse box straddling several pieces of a chain
    // (sweep segments) is loose for a reason that vanishes as proofs
    // refine; looseness that persists is what hurts.
    let r = 0.25 * h;
    // The field varies over a box by at least its range over the corners
    // and centre, where the gradient says nothing: at a tie of min or max
    // it is the average of both sides, 0 on a crease |x - c| (a cut whose
    // face lies exactly on a face of what it cuts).
    let corners: Vec<V3> = sample.iter().flat_map(|p| (0..8).map(move |k| std::array::from_fn(|a| p[a] + if (k >> a) & 1 == 1 { r } else { -r }))).collect();
    let mut at = Vec::new();
    ink_core::eval::eval_points(tape, ctx, &corners, &mut at);
    let boxes: Vec<[I; 3]> = sample.iter().map(|p| p.map(|v| I::new(v - r, v + r))).collect();
    let widths = intervals_over(tape, ctx, &boxes);
    let ratios: Vec<f64> = grads
        .iter()
        .zip(at.chunks(8))
        .zip(&widths)
        .map(|((gr, at), w)| {
            let slope = (gr[1] * gr[1] + gr[2] * gr[2] + gr[3] * gr[3]).sqrt().max(1e-12);
            let (lo, hi) = at.iter().fold((gr[0], gr[0]), |(lo, hi), &v| (lo.min(v), hi.max(v)));
            (w.hi - w.lo) / (2.0 * 3f64.sqrt() * r * slope).max(hi - lo)
        })
        .collect();
    let mut sorted = ratios.clone();
    sorted.sort_by(f64::total_cmp);
    let p90 = sorted[(sorted.len() - 1) * 9 / 10];
    if p90 <= 5.0 {
        return None;
    }
    let share = ratios.iter().filter(|&&x| x > 5.0).count() as f64 / ratios.len() as f64;
    let (k, _) = ratios.iter().enumerate().fold((0, 0.0), |b, (k, &x)| if x > b.1 { (k, x) } else { b });
    Some(Loose { p90, share, max: sorted[sorted.len() - 1], at: sample[k] })
}

/// `eval_interval` of `tape` over each of `boxes`, each read on the tape
/// simplified over a cell of an octree that holds it (the same bounds:
/// what the simplification drops is beaten throughout the cell): the cells
/// halved while they hold boxes that fit in their halves, so a long tape
/// is read whole only near the root, not once a box.
pub fn intervals_over(tape: &Arc<Tape>, ctx: &Ctx, boxes: &[[ink_core::interval::I; 3]]) -> Vec<ink_core::interval::I> {
    use ink_core::interval::I;
    fn split(tape: &Arc<Tape>, ctx: &Ctx, boxes: &[[I; 3]], cell: [I; 3], idx: Vec<usize>, out: &Mutex<Vec<I>>) {
        let mut ch = vec![BOTH; tape.len()];
        counted(tape);
        eval_interval(tape, ctx, cell, Some(&mut ch));
        let sub = tape.narrowed(&ch);
        // Halved on each axis it is more than twice as wide as its boxes
        // (where it can be).
        let mid = cell.map(|x| x.mid());
        let cut: [bool; 3] = std::array::from_fn(|a| {
            let widest = idx.iter().map(|&i| boxes[i][a].hi - boxes[i][a].lo).fold(0.0, f64::max);
            cell[a].hi - cell[a].lo > 2.0 * widest && cell[a].lo < mid[a] && mid[a] < cell[a].hi
        });
        // The part of the cell a box lies in, if one.
        let octant = |b: &[I; 3]| -> Option<usize> {
            if !cut.iter().any(|&c| c) {
                return None;
            }
            (0..3).filter(|&a| cut[a]).try_fold(0, |o, a| {
                if b[a].hi <= mid[a] {
                    Some(o)
                } else if b[a].lo >= mid[a] {
                    Some(o | 1 << a)
                } else {
                    None
                }
            })
        };
        let mut parts: Vec<Vec<usize>> = vec![Vec::new(); 8];
        let mut here = Vec::new();
        for i in idx {
            match octant(&boxes[i]) {
                Some(o) => parts[o].push(i),
                None => here.push(i),
            }
        }
        // A box alone in its octant is read on this cell's tape.
        for p in parts.iter_mut().filter(|p| p.len() == 1) {
            here.append(p);
        }
        let read: Vec<(usize, I)> = here
            .iter()
            .map(|&i| {
                counted(&sub);
                (i, eval_interval(&sub, ctx, boxes[i], None))
            })
            .collect();
        {
            let mut o = out.lock().unwrap();
            read.into_iter().for_each(|(i, v)| o[i] = v);
        }
        parts.into_par_iter().enumerate().filter(|(_, p)| !p.is_empty()).for_each(|(o, p)| {
            let c: [I; 3] = std::array::from_fn(|a| match (cut[a], o >> a & 1) {
                (false, _) => cell[a],
                (true, 0) => I::new(cell[a].lo, mid[a]),
                _ => I::new(mid[a], cell[a].hi),
            });
            split(&sub, ctx, boxes, c, p, out);
        });
    }
    if boxes.is_empty() {
        return Vec::new();
    }
    let cell: [I; 3] = std::array::from_fn(|a| boxes.iter().map(|b| b[a]).reduce(|x, y| I::new(x.lo.min(y.lo), x.hi.max(y.hi))).unwrap());
    let out = Mutex::new(vec![I::point(0.0); boxes.len()]);
    split(tape, ctx, boxes, cell, (0..boxes.len()).collect(), &out);
    out.into_inner().unwrap()
}

// The instructions `intervals_over` read on this thread: its tests' measure
// of its work.
#[cfg(test)]
thread_local!(static READ: std::cell::Cell<usize> = const { std::cell::Cell::new(0) });

/// Count a read of `tape` toward `READ` (tests only).
fn counted(_tape: &Tape) {
    #[cfg(test)]
    READ.with(|r| r.set(r.get() + _tape.len()));
}

/// The undefined cases (division by zero, ln of a negative, ...) the field
/// passes through at up to 4000 of the surface samples and on a 16^3 grid
/// over `lo..hi`, one per kind and construction (the most frequent
/// instruction), most frequent first; and the number of points.
pub fn undefined(tape: &Tape, ctx: &Ctx, surface: &[V3], lo: V3, hi: V3) -> (Vec<Case>, usize) {
    const N: usize = 16;
    let grid = (0..N * N * N).map(|i| -> V3 {
        let k = ink_core::grid::ijk(i, [N; 3]);
        std::array::from_fn(|a| lo[a] + (hi[a] - lo[a]) * (k[a] as f64 + 0.5) / N as f64)
    });
    let pts: Vec<V3> = surface.iter().step_by((surface.len() / 4000).max(1)).copied().chain(grid).collect();
    let mut cases = ink_core::eval::undefined_cases(tape, ctx, &pts);
    let mut seen = std::collections::HashSet::new();
    cases.retain(|u| seen.insert((u.what, u.name)));
    (cases.into_iter().map(|u| Case { what: u.what.to_string(), name: u.name, count: u.count, at: u.at }).collect(), pts.len())
}

/// Area of the plane through `q` with unit normal `n` inside the cube of
/// side `h` centred at `c`: the polygon where the plane cuts the cube.
/// A plane on the cube's boundary, or outside it by rounding (within the
/// tolerance surface() keeps a cell by), counts as just inside: a face
/// lying on a grid plane keeps its area in the one cell that owns it.
pub fn plane_area_in_cell(c: V3, h: f64, q: V3, n: V3) -> f64 {
    let corner = |k: usize| -> V3 { std::array::from_fn(|a| c[a] + if (k >> a) & 1 == 1 { 0.5 * h } else { -0.5 * h }) };
    let reach = 0.5 * h * (n[0].abs() + n[1].abs() + n[2].abs());
    let s = dot(sub(c, q), n);
    if s.abs() > reach * (1.0 + CELL_TOL) {
        return 0.0;
    }
    let s = s.clamp(-reach * (1.0 - CELL_TOL), reach * (1.0 - CELL_TOL));
    let sd: Vec<f64> = (0..8).map(|k| dot(sub(corner(k), c), n) + s).collect();
    // The 12 edges: corner pairs differing in one axis bit.
    let mut pts: Vec<V3> = Vec::with_capacity(6);
    for k in 0..8 {
        for a in 0..3 {
            let m = k | (1 << a);
            if m == k {
                continue;
            }
            let (s0, s1) = (sd[k], sd[m]);
            if (s0 < 0.0) != (s1 < 0.0) {
                let t = s0 / (s0 - s1);
                pts.push(add(corner(k), scale(sub(corner(m), corner(k)), t)));
            }
        }
    }
    if pts.len() < 3 {
        return 0.0;
    }
    // Order around the centroid in the plane, then the fan area.
    let ctr = scale(pts.iter().fold([0.0; 3], |acc, p| add(acc, *p)), 1.0 / pts.len() as f64);
    let e1 = unit(sub(pts[0], ctr));
    let e2 = cross(n, e1);
    pts.sort_by(|x, y| {
        let (dx, dy) = (sub(*x, ctr), sub(*y, ctr));
        libm::atan2(dot(dx, e2), dot(dx, e1)).total_cmp(&libm::atan2(dot(dy, e2), dot(dy, e1)))
    });
    let mut area = 0.0;
    for i in 0..pts.len() {
        let (p0, p1) = (sub(pts[i], ctr), sub(pts[(i + 1) % pts.len()], ctr));
        area += 0.5 * norm(cross(p0, p1));
    }
    area
}

/// The thinnest wall by the ray method: from each surface sample, the
/// distance along the inward normal to where the material ends. Corners
/// take their faces' thickness; fins, ribs and skins their true width. Each
/// ray's exit is found by query::ray_to (interval arithmetic over its
/// segments, so no gap wider than a twentieth of the sample spacing is
/// stepped over), with tapes simplified per block of 4 samples across,
/// shared by the rays through it. Samples on edges, whose rays leave
/// through the neighbouring face within a quarter of the spacing, are
/// skipped (thickness finer than that is not resolved). A ray is searched
/// only as far as the thinnest wall found so far: an exit past it cannot be
/// the thinnest, so the answer is the one tracing every ray in full gives.
/// Rays spread over the samples go first (SPREAD), so how soon the thinnest
/// walls are found does not hang on where they come in the samples' order.
/// Samples whose rays see the same (`alike`: a wall's at every height of an
/// extrusion, a face's across it) share one ray's answer.
/// Returns where the thinnest sample is and its thickness (the first in
/// sample order where several tie; a ray with no exit measures the box's
/// diagonal); None: no sample measured one. Stopped: the call's time limit
/// passed before every ray was traced (crate::stop).
pub fn thickness(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, cells: f64) -> Result<Option<(V3, f64)>, Stopped> {
    if crate::stop::reached() {
        return Err(Stopped);
    }
    let (pts, h) = surface(shape, ctx, lo, hi, cells);
    let local = Local::new(shape, ctx, sub(lo, [h; 3]), 4.0 * h);
    thinnest(&local, &pts, h, norm(sub(hi, lo)) + 4.0 * h)
}

/// The order rays are traced in (`thinnest`): every 1024th sample's, then
/// every 64th's, every 4th's, and the rest; each once.
const SPREAD: [usize; 4] = [1024, 64, 4, 1];

/// `thickness` from its surface samples `pts` spaced `h`, rays traced
/// through `local` (its `ray_evals` counts them) up to `diag`; no ray is
/// begun past the call's time limit.
fn thinnest(local: &Local, pts: &[(V3, V3, f64)], h: f64, diag: f64) -> Result<Option<(V3, f64)>, Stopped> {
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
    let ctx = local.ctx();
    // The thinnest so far (positive floats order as their bits).
    let best = AtomicU64::new(diag.to_bits());
    let stopped = AtomicBool::new(false);
    // Sample i's thickness (None where it cannot be the thinnest), and how
    // far along its ray that read the field.
    let eps = 0.05 * h;
    let ray = |i: usize| -> (Option<f64>, f64) {
        if stopped.load(Relaxed) || crate::stop::reached() {
            stopped.store(true, Relaxed);
            return (None, 0.0);
        }
        let (p, n, _) = pts[i];
        let stop = f64::from_bits(best.load(Relaxed));
        let Some(t) = ray_to_within(local, p, scale(n, -1.0), (1e-3 * h, diag), Side::Outside, eps, stop).map(|x| x.t) else {
            // No exit at all, or none before the thinnest so far.
            return ((stop == diag).then_some(diag), stop);
        };
        let reach = t + eps;
        if t > stop {
            return (None, reach);
        }
        // A wall's far side faces back at the ray (normals within 60
        // degrees of opposite). An exit through a neighbouring face, as
        // from near a drafted rim's edge, measures the corner, not a
        // wall: no thickness there.
        let q = sub(p, scale(n, t));
        let g = eval_grads(&local.tape_over(&q.map(ink_core::interval::I::point)), ctx, &[q])[0];
        let gl = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
        if (gl > 0.0 && (g[1] * n[0] + g[2] * n[1] + g[3] * n[2]) / gl > -0.5) || t < 0.25 * h {
            return (None, reach);
        }
        best.fetch_min(t.to_bits(), Relaxed);
        (Some(t), reach)
    };
    // Each group's first ray, then the rest of each (shared where they see
    // the same): every group has been looked at before any is split.
    let (order, groups) = alike(pts);
    let group = |g: usize| &order[groups[g].clone()];
    let mut first = vec![(None, 0.0); groups.len()];
    for (k, &every) in SPREAD.iter().enumerate() {
        let these: Vec<usize> = (0..groups.len()).step_by(every).filter(|&i| SPREAD[..k].iter().all(|&s| i % s != 0)).collect();
        let got: Vec<(Option<f64>, f64)> = these.par_iter().map(|&g| ray(group(g)[0])).collect();
        these.into_iter().zip(got).for_each(|(g, r)| first[g] = r);
    }
    let mut each: Vec<Option<f64>> = vec![None; pts.len()];
    (0..groups.len()).for_each(|g| each[group(g)[0]] = first[g].0);
    let many: Vec<usize> = (0..groups.len()).filter(|&g| groups[g].len() > 1).collect();
    let got: Vec<Vec<(usize, Option<f64>)>> = many.par_iter().map(|&g| share(local, pts, group(g), first[g], &ray, eps, 4.0 * h)).collect();
    got.into_iter().flatten().for_each(|(i, t)| each[i] = t);
    if stopped.into_inner() {
        return Err(Stopped);
    }
    Ok(pts.iter().zip(each).filter_map(|(s, t)| Some((s.0, t?))).min_by(|a, b| a.1.total_cmp(&b.1)))
}

/// The samples whose rays may see the same, in groups: one outward normal,
/// and the same coordinates along each axis it has a component on, so the
/// rays lie apart only along axes across the normal, along which none moves
/// (a wall's at each height of an extrusion; a flat face's across it). The
/// samples in order of their group, then their own, and each group's span
/// of them, the groups in their first sample's order.
fn alike(pts: &[(V3, V3, f64)]) -> (Vec<usize>, Vec<std::ops::Range<usize>>) {
    // -0 and 0 alike.
    let key = |i: usize| {
        let (p, n, _) = pts[i];
        (n.map(|x| (x + 0.0).to_bits()), std::array::from_fn::<_, 3, _>(|k| (n[k] != 0.0).then(|| p[k].to_bits())))
    };
    let mut order: Vec<usize> = (0..pts.len()).collect();
    order.par_sort_unstable_by(|&i, &j| key(i).cmp(&key(j)).then(i.cmp(&j)));
    let mut groups: Vec<std::ops::Range<usize>> = Vec::new();
    for (k, &i) in order.iter().enumerate() {
        match groups.last_mut() {
            Some(g) if key(order[g.start]) == key(i) => g.end = k + 1,
            _ => groups.push(k..k + 1),
        }
    }
    groups.sort_unstable_by_key(|g| order[g.start]);
    (order, groups)
}

/// The answers of samples `members` (`alike`) by `ray` (an answer, and how
/// far along its ray it read the field): the first's (given) shared by all
/// where the field over the box their rays span that far, for its sign
/// alone (query::signed), reads none of the axes they lie apart along: each
/// ray then sees the field the first sees, moved along those axes, so meets
/// the surface where it does. Else they are taken a block `block` across at
/// a time, each block's first traced and shared so, or each of its rays
/// traced. Each box asked counts as an evaluation along a ray.
fn share(
    local: &Local,
    pts: &[(V3, V3, f64)],
    members: &[usize],
    first: (Option<f64>, f64),
    ray: &(dyn Fn(usize) -> (Option<f64>, f64) + Sync),
    eps: f64,
    block: f64,
) -> Vec<(usize, Option<f64>)> {
    let n = pts[members[0]].1;
    let apart: Vec<usize> = (0..3).filter(|&k| n[k] == 0.0).collect();
    // Whether the rays from `these`, read `reach` along, all see the same.
    let same = |these: &[usize], reach: f64| {
        if these.len() == 1 {
            return true;
        }
        let mut b = [ink_core::interval::I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
        for &i in these {
            let p = pts[i].0;
            for e in [p, sub(p, scale(n, reach + eps))] {
                for k in 0..3 {
                    b[k] = ink_core::interval::I::new(b[k].lo.min(e[k]).next_down(), b[k].hi.max(e[k]).next_up());
                }
            }
        }
        local.count_ray_evals(1);
        b.iter().all(|x| x.lo.is_finite() && x.hi.is_finite()) && {
            let read = crate::query::signed(&local.tape_over(&b), local.ctx(), b).space_read();
            !apart.iter().any(|&k| read[k])
        }
    };
    if same(members, first.1) {
        return members.iter().map(|&i| (i, first.0)).collect();
    }
    let mut blocks: std::collections::BTreeMap<Vec<i64>, Vec<usize>> = Default::default();
    for &i in members {
        blocks.entry(apart.iter().map(|&k| (pts[i].0[k] / block).floor() as i64).collect()).or_default().push(i);
    }
    let blocks: Vec<Vec<usize>> = blocks.into_values().collect();
    blocks
        .par_iter()
        .flat_map_iter(|these| {
            let r = if these[0] == members[0] { first } else { ray(these[0]) };
            let shared = same(these, r.1);
            let rest = these[1..].iter().map(move |&i| (i, if shared { r.0 } else { ray(i).0 }));
            std::iter::once((these[0], r.0)).chain(rest).collect::<Vec<_>>()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Bounds over many small boxes of a long field (a row of 200 balls)
    /// are each box's own, to the bit, and read far less than the whole
    /// field once a box: 2,000 boxes along the row read 315 times its
    /// length (2,000 times, a box at a time).
    #[test]
    fn bounds_over_many_boxes_read_the_field_simplified_near_each() {
        use ink_core::graph::{BinOp, UnOp, Var};
        use ink_core::interval::I;
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let mut balls = Vec::new();
        for i in 0..200 {
            let c = b.c(3.0 * i as f64);
            let dx = b.bin(BinOp::Sub, x, c);
            let (dx, dy, dz) = (b.un(UnOp::Square, dx), b.un(UnOp::Square, y), b.un(UnOp::Square, z));
            let r = b.bin(BinOp::Add, dx, dy);
            let r = b.bin(BinOp::Add, r, dz);
            let r = b.un(UnOp::Sqrt, r);
            let one = b.c(1.0);
            balls.push(b.bin(BinOp::Sub, r, one));
        }
        let f = balls.into_iter().reduce(|p, q| b.bin(BinOp::Min, p, q)).unwrap();
        let tape = Arc::new(b.tape(f));
        let ctx = Ctx::new(0.0, vec![]);
        let boxes: Vec<[I; 3]> = (0..2000).map(|k| [0.3 * k as f64, 0.9, -0.2].map(|c| I::new(c - 0.01, c + 0.01))).collect();
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        let (got, read) = pool.install(|| {
            READ.with(|r| r.set(0));
            (intervals_over(&tape, &ctx, &boxes), READ.with(|r| r.get()))
        });
        for (bx, v) in boxes.iter().zip(&got) {
            let w = eval_interval(&tape, &ctx, *bx, None);
            assert_eq!((v.lo.to_bits(), v.hi.to_bits()), (w.lo.to_bits(), w.hi.to_bits()), "over {bx:?}");
        }
        assert!(read < 400 * tape.len(), "{read} instructions read, the field {} long, over {} boxes", tape.len(), boxes.len());
    }

    #[test]
    fn a_surface_sample_has_a_normal() {
        // cbrt(x - 1/4) has an infinite gradient on its zero set x = 1/4,
        // through cell centres here: f / |grad| = 0 would pass for on the
        // surface, and grad / |grad| give a NaN normal that
        // plane_area_in_cell panics clamping by (5225K722's elbow, where
        // a cube root in its fillet's nearest ellipse point is 0 at one
        // cell centre).
        use ink_core::graph::{BinOp, Graph, Var};
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let (quarter, third) = (g.constant(0.25), g.constant(1.0 / 3.0));
        let d = g.sub(x, quarter);
        let nd = g.neg(d);
        let below = g.bin(BinOp::Pow, nd, third);
        let below = g.neg(below);
        let above = g.bin(BinOp::Pow, d, third);
        let f = g.select(d, below, above);
        let (pts, _) = surface(&Arc::new(Tape::new(&g, f)), &Ctx::new(0.0, vec![]), [-1.0; 3], [1.0; 3], 4.0);
        assert!(pts.iter().all(|(_, n, w)| (norm(*n) - 1.0).abs() < 1e-12 && w.is_finite()));
    }

    #[test]
    fn a_ray_is_traced_only_as_far_as_the_thinnest_wall_so_far() {
        // A 0.8 fin on a 40 mm block: the fin's rays end within 0.8, the
        // block's run 40 through it. Only the thinnest is wanted, so a ray
        // stops once it is past the thinnest found so far: the answer is
        // the fin's, for a fraction of the evaluations tracing every ray
        // through gives.
        use ink_core::graph::{Graph, UnOp, Var};
        let mut g = Graph::new();
        let xyz = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let mut cuboid = |c: V3, half: V3| {
            let d: Vec<_> = (0..3)
                .map(|k| {
                    let (c, r) = (g.constant(c[k]), g.constant(half[k]));
                    let o = g.sub(xyz[k], c);
                    let a = g.un(UnOp::Abs, o);
                    g.sub(a, r)
                })
                .collect();
            let m = g.max(d[0], d[1]);
            g.max(m, d[2])
        };
        let block = cuboid([0.0; 3], [20.0; 3]);
        let fin = cuboid([10.0, 0.0, 25.0], [0.4, 10.0, 5.0]);
        let f = g.min(block, fin);
        let (t, ctx) = (Arc::new(Tape::new(&g, f)), Ctx::new(0.0, vec![]));
        let (lo, hi) = ([-21.0, -21.0, -21.0], [21.0, 21.0, 31.0]);
        let (pts, h) = surface(&t, &ctx, lo, hi, 64.0);
        let diag = norm(sub(hi, lo)) + 4.0 * h;
        // Each count is its own Local's: no other test's rays reach it.
        let local = || Local::new(&t, &ctx, sub(lo, [h; 3]), 4.0 * h);
        // On one thread the rays are taken in their order, so the count
        // does not depend on how threads are scheduled. Taken in the
        // samples' own order, with the fin's samples (highest in z) first,
        // or last: rays spread over the samples go first, so the fin is
        // found early whichever it is, and the counts stay close.
        let one = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        let ordered = |by: f64| {
            let mut order = pts.clone();
            order.sort_by(|a, b| (by * a.0[2]).total_cmp(&(by * b.0[2])));
            order
        };
        let counts: Vec<((V3, f64), usize)> = [pts.clone(), ordered(-1.0), ordered(1.0)]
            .iter()
            .map(|order| {
                one.install(|| {
                    let measured = local();
                    (thinnest(&measured, order, h, diag).unwrap().unwrap(), measured.ray_evals())
                })
            })
            .collect();
        // The same wall in every order (where several tie, maybe another of them).
        for &((at, thinnest), _) in &counts {
            assert!((thinnest - 0.8).abs() < 1e-9 && (at[0] - 10.0).abs() <= 0.4 + 1e-9, "{thinnest} at {at:?}");
        }
        assert_eq!(thickness(&t, &ctx, lo, hi, 64.0), Ok(Some(counts[0].0)));
        let (least, most) = counts.iter().fold((usize::MAX, 0), |(l, m), c| (l.min(c.1), m.max(c.1)));
        assert!(2 * least > most, "{counts:?}: the order decides how soon the fin is found");
        // Every ray traced through, as thickness traced them before.
        let through = local();
        pts.par_iter().for_each(|&(p, n, _)| _ = crate::query::ray_to(&through, p, scale(n, -1.0), (1e-3 * h, diag), Side::Outside, 0.05 * h));
        let full = through.ray_evals();
        assert!(3 * most < full, "{most} evaluations, {full} tracing every ray through");
    }

    /// A stroke 0.8 wide and 10 long, on the diagonal, extruded 0 to 2 (an
    /// emax of the profile and its span), sampled 0.05 apart: its walls'
    /// rays at every height see what those at one height see (the field
    /// along them, for its sign, reads no z), and its top's and bottom's
    /// see the same at every point of the face inside the profile's walls
    /// (it reads neither x nor y there). The thinnest is the profile's,
    /// for about the evaluations the profile asked in the plane takes
    /// (every ray traced: 1.36 M, 10 times the profile's).
    #[test]
    fn an_extrusions_rays_cost_what_its_profiles_do() {
        use crate::test_fields::B;
        use ink_core::graph::{BinOp, Var};
        let thinnest = |solid: bool| {
            let mut b = B::new();
            let (x, y) = (b.v(Var::X), b.v(Var::Y));
            let k = b.c(std::f64::consts::FRAC_1_SQRT_2);
            let (u, v) = (b.bin(BinOp::Add, x, y), b.bin(BinOp::Sub, x, y));
            let (u, v) = (b.bin(BinOp::Mul, u, k), b.bin(BinOp::Mul, v, k));
            let (across, along) = (b.slab(u, 0.0, 0.4), b.slab(v, 0.0, 5.0));
            let mut f = b.bin(BinOp::Max, across, along);
            if solid {
                let z = b.v(Var::Z);
                let span = b.slab(z, 1.0, 1.0);
                f = b.bin(BinOp::Emax, f, span);
            }
            let (t, ctx) = (Arc::new(b.tape(f)), Ctx::new(0.0, vec![]));
            let (lo, hi) = ([-4.0, -4.0, 0.0], [4.0, 4.0, if solid { 2.0 } else { 0.0 }]);
            let cells = crate::query::max_extent(lo, hi) / 0.05;
            let (pts, h) = surface(&t, &ctx, lo, hi, cells);
            let local = Local::new(&t, &ctx, sub(lo, [h; 3]), 4.0 * h);
            let got = thinnest(&local, &pts, h, norm(sub(hi, lo)) + 4.0 * h).unwrap().unwrap();
            (got.1, local.ray_evals())
        };
        let ((solid, rays), (plane, plane_rays)) = (thinnest(true), thinnest(false));
        assert!((solid - 0.8).abs() < 1e-9 && (plane - 0.8).abs() < 1e-9, "{solid} {plane}");
        assert!(rays <= 2 * plane_rays, "extruded: {rays} evaluations, the profile {plane_rays}");
    }

    #[test]
    fn a_profile_is_sampled_in_the_box_asked() {
        // A disc of radius 10 that reads no z (a 2D profile: a prism), asked
        // at z = 0 alone, 64 cells across: its rim's cells in at most the
        // four layers of a block, not in every layer of the root cube (64).
        use ink_core::graph::{Graph, UnOp, Var};
        let mut g = Graph::new();
        let (x, y) = (g.var(Var::X), g.var(Var::Y));
        let (xx, yy) = (g.mul(x, x), g.mul(y, y));
        let s = g.add(xx, yy);
        let d = g.un(UnOp::Sqrt, s);
        let r = g.constant(10.0);
        let f = g.sub(d, r);
        let (pts, h) = surface(&Arc::new(Tape::new(&g, f)), &Ctx::new(0.0, vec![]), [-10.0, -10.0, 0.0], [10.0, 10.0, 0.0], 64.0);
        let per_layer = 2.0 * std::f64::consts::PI * 10.0 / h;
        assert!((pts.len() as f64) < 4.0 * 1.5 * per_layer, "{} samples, {per_layer} a layer", pts.len());
    }

    #[test]
    fn a_crease_is_as_loose_as_its_variation() {
        // max(x, -x): a cut whose face lies exactly on a face of what it
        // cuts (KQ2L06-01NS's release button in its housing, quarter turns
        // exact). On x = 0 the gradient is the average of both sides, 0, but
        // over a box the field varies as its bounds say: a bound by the
        // gradient alone is 5.8e11 times loose.
        use ink_core::graph::{Graph, Var};
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let nx = g.neg(x);
        let f = g.max(x, nx);
        let pts: Vec<V3> = (0..100).map(|k| [0.0, 0.1 * k as f64, 0.0]).collect();
        assert!(loose(&Arc::new(Tape::new(&g, f)), &Ctx::new(0.0, vec![]), &pts, 0.1).is_none());
    }

    #[test]
    fn a_survey_samples_what_it_reads() {
        // A cube's surface crosses 98,000 cells 128 across: a survey read
        // 4,000 of them and projected all (a large assembly's check: 17 M samples,
        // 70% of its time). A sphere a hundredth across crosses few at 16:
        // it is sampled finer.
        use ink_core::graph::{Graph, UnOp, Var};
        let mut g = Graph::new();
        let xyz = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let a = xyz.map(|v| g.un(UnOp::Abs, v));
        let (m, r) = (g.max(a[0], a[1]), g.constant(50.0));
        let m = g.max(m, a[2]);
        let f = g.sub(m, r);
        let ctx = Ctx::new(0.0, vec![]);
        let n: usize = samples(&Arc::new(Tape::new(&g, f)), &ctx, [-50.0; 3], [50.0; 3], 128.0).0.iter().map(|b| b.0.len()).sum();
        assert!((SURVEYED / 2..2 * SURVEYED).contains(&n), "{n}");
        let sq = xyz.map(|v| g.mul(v, v));
        let (s, r) = (g.add(sq[0], sq[1]), g.constant(0.5));
        let s = g.add(s, sq[2]);
        let d = g.un(UnOp::Sqrt, s);
        let f = g.sub(d, r);
        let n: usize = samples(&Arc::new(Tape::new(&g, f)), &ctx, [-50.0; 3], [50.0; 3], 128.0).0.iter().map(|b| b.0.len()).sum();
        assert!(n > 4, "{n}");
    }
}
