//! Points sampled on the exact surface (not a mesh): the surface cells of an
//! interval octree, projected onto the zero set by Newton steps, with their
//! normals and areas; the names that can form the surface; wall thickness
//! by rays.

use crate::octree::{Cell, Grid, Local, descend_padded, root_for};
use crate::query::{Side, ray_to};
use ink_core::V3;
use ink_core::eval::{Ctx, decided_names, eval_grads, eval_interval, names_on_paths};
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
    const B: i64 = 4;
    let ext = crate::query::max_extent(lo, hi);
    let h = ext / cells;
    let pad = [lo[0] - h, lo[1] - h, lo[2] - h];
    let (grid, top) = root_for(pad, [hi[0] + h, hi[1] + h, hi[2] + h], h, B);
    let blocks: Vec<Vec<(V3, V3, f64)>> = descend_padded(
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
                let outside = (0..8).all(|k| fc[ink_core::grid::index([cx + (k & 1), cy + ((k >> 1) & 1), cz + ((k >> 2) & 1)], [n1; 3])] >= 0.0);
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
            Some(out)
        },
        &|_, _| None,
    );
    (blocks.into_iter().flatten().collect(), h)
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
    // Its names and cases give tags by the model's numbers.
    let raw = |fp: &mut ink_core::graph::Fingerprint| ink_core::tape::Names::of(&[shape], ctx).raw(fp);
    crate::cache::memo_region(
        "survey",
        &[shape],
        ctx,
        (lo, hi),
        |fp| raw(fp.f64(cells).word(names as u64)),
        |ctx| {
            let (samples, h) = samples(shape, ctx, lo, hi, cells);
            let interior = interior(shape, ctx, &samples, h);
            let pts: Vec<V3> = samples.into_iter().map(|p| p.0).collect();
            let (undefined, points) = undefined(shape, ctx, &pts, lo, hi);
            let (formed, unseen) = if names { surface_names(shape, ctx, lo, hi, cells) } else { Default::default() };
            let loose = loose(shape, ctx, &pts, crate::query::max_extent(lo, hi) / cells);
            Survey { h, formed, unseen, loose, undefined, points, interior }
        },
    )
}

/// About SURVEYED samples of the surface, one per cell crossed (`surface`),
/// no finer than `cells` across: the cells a surface crosses go as the
/// square of the cells across, so a pass 16 across says how fine to go.
fn samples(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, cells: f64) -> (Vec<(V3, V3, f64)>, f64) {
    let first = surface(shape, ctx, lo, hi, cells.min(16.0));
    let n = first.0.len();
    if n >= SURVEYED || cells <= 16.0 {
        return first;
    }
    surface(shape, ctx, lo, hi, (16.0 * (SURVEYED as f64 / n.max(1) as f64).sqrt()).min(cells))
}

/// Whether the field is negative a thousandth of the spacing `h` inside
/// any surface sample (point, outward normal): a solid's are, at once (so
/// this stops at the first block that has one).
fn interior(shape: &Tape, ctx: &Ctx, samples: &[(V3, V3, f64)], h: f64) -> bool {
    let mut f = Vec::new();
    samples.chunks(1024).any(|c| {
        let inner: Vec<V3> = c.iter().map(|(p, n, _)| sub(*p, scale(*n, 1e-3 * h))).collect();
        ink_core::eval::eval_points(shape, ctx, &inner, &mut f);
        f.iter().any(|&v| v < 0.0)
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
    if r.lo > 0.0 || r.hi < 0.0 {
        return;
    }
    let sub = if ch.iter().any(|&c| c != BOTH) { Arc::new(tape.simplify(&ch)) } else { tape.clone() };
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

fn loose(tape: &Tape, ctx: &Ctx, pts: &[V3], h: f64) -> Option<Loose> {
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
    let ratios: Vec<f64> = sample
        .iter()
        .zip(&grads)
        .zip(at.chunks(8))
        .map(|((p, gr), at)| {
            let b = p.map(|v| I::new(v - r, v + r));
            let w = ink_core::eval::eval_interval(tape, ctx, b, None);
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

/// Wall thickness by the ray method: from each surface sample, the distance
/// along the inward normal to where the material ends. Corners take their
/// faces' thickness; fins, ribs and skins their true width. Each ray's exit
/// is found by query::ray_to (interval arithmetic over its segments, so no
/// gap wider than a twentieth of the sample spacing is stepped over), with
/// tapes simplified per block of 4 samples across, shared by the rays
/// through it. Samples on edges, whose rays leave through the neighbouring
/// face within a quarter of the spacing, are skipped (thickness finer than
/// that is not resolved). Returns (point, outward normal, area weight,
/// thickness).
pub fn thickness(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, cells: f64) -> Vec<(V3, V3, f64, f64)> {
    let (pts, h) = surface(shape, ctx, lo, hi, cells);
    let diag = norm(sub(hi, lo)) + 4.0 * h;
    let local = Local::new(shape, ctx, sub(lo, [h; 3]), 4.0 * h);
    let exits: Vec<Option<f64>> = pts.par_iter().map(|(p, n, _)| ray_to(&local, *p, scale(*n, -1.0), (1e-3 * h, diag), Side::Outside, 0.05 * h).map(|x| x.t)).collect();
    // The field's gradient where each ray left, in one batch.
    let far: Vec<V3> = pts.iter().zip(&exits).filter_map(|((p, n, _), t)| t.map(|t| sub(*p, scale(*n, t)))).collect();
    let mut g = local.grads(&far).into_iter();
    pts.iter()
        .zip(exits)
        .filter_map(|(&(p, n, w), t)| {
            let Some(t) = t else { return Some((p, n, w, diag)) };
            // A wall's far side faces back at the ray (normals within 60
            // degrees of opposite). An exit through a neighbouring face, as
            // from near a drafted rim's edge, measures the corner, not a
            // wall: no thickness there.
            let g = g.next().unwrap();
            let gl = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
            if gl > 0.0 && (g[1] * n[0] + g[2] * n[1] + g[3] * n[2]) / gl > -0.5 {
                return None;
            }
            (t >= 0.25 * h).then_some((p, n, w, t))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert!(loose(&Tape::new(&g, f), &Ctx::new(0.0, vec![]), &pts, 0.1).is_none());
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
        let n = samples(&Arc::new(Tape::new(&g, f)), &ctx, [-50.0; 3], [50.0; 3], 128.0).0.len();
        assert!((SURVEYED / 2..2 * SURVEYED).contains(&n), "{n}");
        let sq = xyz.map(|v| g.mul(v, v));
        let (s, r) = (g.add(sq[0], sq[1]), g.constant(0.5));
        let s = g.add(s, sq[2]);
        let d = g.un(UnOp::Sqrt, s);
        let f = g.sub(d, r);
        let n = samples(&Arc::new(Tape::new(&g, f)), &ctx, [-50.0; 3], [50.0; 3], 128.0).0.len();
        assert!(n > 4, "{n}");
    }
}
