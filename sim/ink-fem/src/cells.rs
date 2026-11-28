//! The grid's cells: what a cut cell's samples are, the pieces of material they
//! make, the cells' nodes, and the cell at a point.

use crate::eval::{Ctx, eval_interval, eval_points};
use crate::laminate::Laminate;
use crate::occupancy::{Classified, Cov};
use crate::octree::Grid;
use crate::tape::Tape;
use ink_core::{V3, grid, interval, tape, v3};
use ink_kernel::mesh::find;
use std::sync::Arc;

/// Sample count per axis inside cut cells (4^3 = 64 = one SIMD batch).
pub(crate) const K: usize = 4;
/// Stiffness kept in the void part of cut cells, relative to solid (finite
/// cell regularisation; keeps nearly-empty cells from making K singular).
pub(crate) const ALPHA: f64 = 1e-6;
/// Cut cells holding less material than this fraction are dropped (slivers),
/// once cells are split: a split cell's children judge all of its material.
pub(crate) const MIN_FRACTION: f64 = 0.03;
/// A cell lends its traction to a piece beside it where at least this
/// share of it is free of every hold. A lender's constant goes as the
/// inverse of its free share, so a half keeps it near twice a whole
/// cell's; the hold's positivity does not need it, and nothing measured
/// separates it from its neighbours (at 0.25 and 0.75, thirteen holds,
/// level, tilted, curved and beside a guide, read within 0.03% of their
/// readings at 0.5, in as many iterations or fewer).
pub(crate) const LEND: f64 = 0.5;
/// The most the void's kept stiffness (ALPHA) moves a reading, relative: a kept cut cell holds at least
/// MIN_FRACTION material, so its void stiffens it by at most ALPHA (1 - MIN_FRACTION) / MIN_FRACTION of
/// its material's. Grids whose readings differ by no more are alike (claims::grids_error).
pub const VOID_SHARE: f64 = ALPHA * (1.0 - MIN_FRACTION) / MIN_FRACTION;

// A cut cell's samples as bits, sample i = a + K b + K^2 c (`sample_xi`).
const _: () = assert!(K * K * K == 64);
/// Every sample of a cell: a full cell's material.
pub(crate) const ALL: u64 = u64::MAX;

/// Per sample of a cell of side `size` (at `pts`): the chord across it
/// (mm) of each phase it holds, material and air (`MATERIAL`, `AIR`): the
/// run of that phase along the field's gradient through its point nearest
/// the sample's centre (within a sample), both ways to the other phase,
/// each end at the first point past it: the longest run the points allow,
/// so a feature is thin only where it surely is (points by eighths of a
/// sample to a sample each way, then by halves to a sample past `reach`
/// (mm) each way: a run reaching an end is at least `reach` long, its own
/// end within a sample of the centre), on the whole field `shape` (a cell's
/// own tape drops the branches that do not reach into it: a wall's far face
/// beyond the cell). Where the gradient cancels (midway
/// between a thin wall's faces, or a slit's), the least such run along the
/// axes. Half `reach` or more deep in a phase by the field's value over its
/// gradient (as occupancy reads depth: a field scaled up is no deeper),
/// infinite. The
/// chord runs through the phase, not the cell: a sliver of a cell at a
/// thick body's face spans the body, and its air the open space beyond.
pub(crate) fn sample_chords(shape: &Arc<Tape>, ctx: &Ctx, pts: &[V3], size: f64, reach: f64) -> Vec<[Option<f64>; 2]> {
    if pts.is_empty() {
        return Vec::new();
    }
    let s = size / K as f64;
    let far = (reach / s).ceil() as i64 + 1;
    // The shape over the box the lines span.
    let t = &{
        let span = far as f64 * s;
        let b: [interval::I; 3] = std::array::from_fn(|a| {
            let (lo, hi) = pts.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| (lo.min(p[a]), hi.max(p[a])));
            interval::I::new(lo - span, hi + span)
        });
        let mut ch = vec![tape::BOTH; shape.len()];
        eval_interval(shape, ctx, b, Some(&mut ch));
        shape.narrowed(&ch)
    };
    // Offsets along a line, in samples, ascending.
    let offsets: Vec<f64> = {
        let fine = (-8..=8).map(|k| k as f64 / 8.0);
        let coarse = (3..=2 * far).map(|k| k as f64 / 2.0);
        let mut v: Vec<f64> = coarse.clone().map(|d| -d).chain(fine).chain(coarse).collect();
        v.sort_by(f64::total_cmp);
        v
    };
    let n = offsets.len();
    let mid = offsets.iter().position(|&d| d == 0.0).expect("the centre");
    let mut out = vec![[None, None]; pts.len()];
    let g = crate::eval::eval_grads(t, ctx, pts);
    let (mut lines, mut along): (Vec<usize>, Vec<V3>) = (Vec::new(), Vec::new());
    for (i, g) in g.iter().enumerate() {
        let slope = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
        if slope > 0.0 && g[0].abs() / slope >= 0.5 * reach {
            out[i][usize::from(g[0] > 0.0)] = Some(f64::INFINITY);
            continue;
        }
        let dirs: Vec<V3> = match v3::direction([g[1], g[2], g[3]]) {
            Some(d) => vec![d],
            None => vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        };
        for d in dirs {
            lines.push(i);
            along.extend(offsets.iter().map(|&k| v3::add(pts[i], v3::scale(d, k * s))));
        }
    }
    let mut f = Vec::new();
    eval_points(t, ctx, &along, &mut f);
    // A line along the field's gradient crosses the surface square to it only where the field is a
    // distance; another field's gradient off its zero turns where it scales (a bore written (1 - r)(1 +
    // 0.95 sin 25x) pointed along x, and its air read thinner than its diameter). On its zero a field's
    // gradient is the surface's normal, whatever its scale: where that turns from the line's by more
    // than TURNED, the line is taken again along it, through the sample.
    let mut again: Vec<(usize, V3)> = Vec::new();
    {
        let mut at: Vec<(usize, V3)> = Vec::new();
        for (j, &i) in lines.iter().enumerate() {
            let line = &f[j * n..(j + 1) * n];
            let cross = (1..n).filter(|&k| (line[k - 1] > 0.0) != (line[k] > 0.0)).min_by_key(|&k| (k as i64 - mid as i64).abs());
            if let Some(k) = cross {
                at.push((j, along[j * n + k]));
                let _ = i;
            }
        }
        let pts_at: Vec<V3> = at.iter().map(|x| x.1).collect();
        let gc = crate::eval::eval_grads(t, ctx, &pts_at);
        let slope = |g: &[f64; 4]| (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
        for ((j, _), gc) in at.iter().zip(&gc) {
            let Some(normal) = v3::direction([gc[1], gc[2], gc[3]]) else { continue };
            let d = v3::direction(v3::sub(along[j * n + mid + 1], along[j * n + mid])).unwrap_or(normal);
            // (A distance's slope is one everywhere, a scaled one's the same everywhere: where it turns
            // at a crease, a corner's other face, its line is the corner's as it should be.)
            let (here, there) = (slope(&g[lines[*j]]), slope(gc));
            let scaled = !(there <= SCALED * here && here <= SCALED * there);
            if scaled && v3::dot(normal, d).abs() < TURNED {
                again.push((*j, normal));
            }
        }
    }
    if !again.is_empty() {
        let line_of = |j: usize, d: V3| -> Vec<V3> { offsets.iter().map(|&k| v3::add(pts[lines[j]], v3::scale(d, k * s))).collect() };
        let pts2: Vec<V3> = again.iter().flat_map(|&(j, d)| line_of(j, d)).collect();
        let mut f2 = Vec::new();
        eval_points(t, ctx, &pts2, &mut f2);
        for (q, &(j, _)) in again.iter().enumerate() {
            f[j * n..(j + 1) * n].copy_from_slice(&f2[q * n..(q + 1) * n]);
        }
    }
    for (j, &i) in lines.iter().enumerate() {
        let line = &f[j * n..(j + 1) * n];
        for (phase, chord) in out[i].iter_mut().enumerate() {
            let of = |k: usize| (line[k] > 0.0) == (phase == 1);
            // The phase's point nearest the centre, within a sample, and its run.
            let Some(at) = (0..=8).flat_map(|d| [mid - d, mid + d]).find(|&k| of(k)) else { continue };
            let hi = (at..n).take_while(|&k| of(k)).last().expect("at");
            let lo = (0..=at).rev().take_while(|&k| of(k)).last().expect("at");
            let c = if lo == 0 || hi + 1 == n { f64::INFINITY } else { (offsets[hi + 1] - offsets[lo - 1]) * s };
            *chord = Some(chord.map_or(c, |v: f64| v.min(c)));
        }
    }
    out
}

/// The cosine under which the surface's normal where a chord's line crosses it turns the line along it
/// (`sample_chords`): ten degrees, past a distance field's rounding.
const TURNED: f64 = 0.985;

/// How far a field's slope may change between a sample and the surface its chord's line crosses for the
/// field to be read as a distance there, scaled or not (`sample_chords`).
const SCALED: f64 = 1.5;

/// A feature thinner than this many of its cells across is thin
/// (`sample_chords`), of either phase: a cell holding any is split until
/// it is at least that many of its (finer) cells across (`MAX_DEPTH`). A
/// cell's trilinear field cannot hold such a feature: across thin material
/// it holds it as if it filled the cell (a ligament under a slit 7 to 68
/// times too stiff, a fin's or a web's stress and motion unread in cells
/// under half full); across thin air, material on its two sides moves as
/// one (a slit narrower than the cells shut, its sides tied by the nodes
/// they share, its hinge held by one cell's field).
pub(crate) const THIN: f64 = 2.0;

/// The phases a cell holds features of (`Classified::thin`, `resolved`),
/// and its material's properties (`PROPERTY`).
pub(crate) const MATERIAL: u8 = 1;
pub(crate) const AIR: u8 = 2;
/// A cell whose material is not one material across it, its variation not
/// resolved at its size (`variation_resolved`): a layer, a strip, an
/// interface between materials, a step steeper than the cell. Its element
/// takes one coefficient, so on every grid of a claim the feature is lost
/// alike (a 0.1 mm insulating layer in 0.2 mm cells read 25.5 degrees for
/// 49.8): it is split as a thin feature is, and left at the floor beyond
/// where its materials are resolved (`material_held`) it decides no claim.
pub(crate) const PROPERTY: u8 = 4;

/// A cell's box as its material is judged over (`PROPERTY`): shrunk by
/// this share of its side, so a material's face on a cell's face (round
/// coordinates, a body's on the lattice) leaves both cells of one
/// material, as a face shared by closed boxes would not.
pub(crate) const PROPERTY_SHRINK: f64 = 1e-6;

/// A property's range over a box this small a share of its value is its
/// bounds' outward rounding, not a variation (`PROPERTY`).
pub(crate) const PROPERTY_ROUNDING: f64 = 1e-9;

/// How much steeper than its range a property may be over a box for its
/// variation to be resolved there: the box's side times the bound on the
/// field's slope over it (`interval::lipschitz`), over its range. 1 for an
/// affine field, 4 for a quadratic least or greatest at the box's centre
/// (the steepest `variation_resolved` takes, its order 2): twice that.
/// Ranges alone are blind to what lies between their boxes' faces: a jump
/// film (a `select` between values, a step) in a floor box beside two
/// slivers tuned to grow its ranges as a smooth field's grow passed
/// (`< 110 MPa` held at 96.6 against 129). A jump has no finite slope bound
/// over the box, so its cells split and, at the floor, stay not settled.
/// A steep film of finite slope passes it beside slivers that widen the
/// box's range, which any ratio to that range allows: what the cells' own
/// values see bounds it (`PROPERTY_UNSEEN`).
pub(crate) const PROPERTY_STEEP: f64 = 8.0;

/// How far a property may reach over a box beyond what the box's centre and
/// the centres of the 26 boxes of its size round it (inside the part) take,
/// as a share of that value, for its variation to be resolved there (every
/// floor box, `property`): a cell takes its centre's value (a cut cell, its
/// samples'), so what lies beyond them is missed alike on every grid, and
/// the ladder cannot see it. Neighbour centres outside the part (in air) are
/// ignored without extrapolation, so every bound comes strictly from values
/// taken inside the part (preventing both an air-clamped outer skin and an
/// inward dip from widening the face's allowed range). A field monotone along each axis reaches
/// nothing beyond them; a smooth one only near its extremes, by its
/// curvature times the square of the cell's side. A film between the centres
/// reaches its own contrast beyond them, which nothing else in the box can
/// hide: a sliver beside it moves no centre's value unless it lies at one,
/// and there the cells take it. Ranges and slopes alone are blind to this: a
/// film 2667 against 1972 on steep finite slopes (2e5 MPa/mm), beside one
/// soft dip that made its floor box's ranges grow as a smooth field's and
/// its slope ratio pass, held `< 110 MPa` at 96.9 against about 129. Each
/// side is held to the centres' on that side, a share of its own value: so
/// k's least holds 1/k's greatest to the same share (a soft film in series),
/// and nothing a centre takes on the other side widens it.
pub(crate) const PROPERTY_UNSEEN: f64 = 0.01;

/// Whether a material property's variation over a cell is resolved at the
/// cell's size: its ranges over boxes about the cell's centre a quarter, a
/// half and the whole of its side, `r`, `r2`, `r4`. A claim's ladder reads
/// a value on three sizes and takes its error from how it moves between
/// them (ink-lang `claims::grid_error`): it measures an error that shrinks
/// with the cells and is blind to one alike on all three. A cell's
/// coefficient is off its material by at most its range, which shrinks
/// with the box where the property is smooth and stays where it steps. Read
/// as such a ladder, the variation is resolved where its range grows as a
/// smooth field's does and that ladder's estimate of `r`'s distance from 0
/// covers `r`, the estimate bounded the safe way for this question:
///
/// - its order held to at least 1. Below it the grid formula's estimate,
///   1.25 d1 / (2^p - 1), grows without bound, conservative for a value's
///   error but here what lets a step pass: a step J on a slope growing
///   slower than linearly (r, r + d, r + d (1 + e)) read as resolved for J
///   up to 1.25 d / e. Held at order 1 the estimate is at most 1.25 d1, so
///   a resolved range at least 1.8 times its quarter's over the half: a
///   step on a slope passes only under a quarter of the slope's own range
///   over the quarter box (r = J + g, g linear: J + g <= 1.25 g);
/// - its order at most 2, a smooth field's at an extremum (its range
///   quadratic in the box about it): faster, the whole box holds what its
///   half does not, a step in the box's outer part. What lies there moves
///   only `r4`, so a step passes beside variation nearer the centre that
///   grows `r2` (a step, or a slope): a cell is judged over every box of
///   the floor's size too (`property`).
///
/// Affine: (r, 2r, 4r), order 1, 1.25 r. Least or greatest at the centre:
/// (r, 4r, 16r), order 2. Between (an extremum off centre, a sine): order
/// 1 to 2, its quarter box linear where its half is, r = r2 - r, at the
/// estimate's floor. A step through the quarter box: alike on all three,
/// not resolved; one narrower than the cell, barely growing, nor. A range
/// grows with its box, so a smaller one never exceeds a larger.
pub fn variation_resolved(r: f64, r2: f64, r4: f64) -> bool {
    if r4 == 0.0 {
        return true;
    }
    let (d1, d2) = (r2 - r, r4 - r2);
    let round = 1.0 + PROPERTY_ROUNDING;
    if !(r.is_finite() && r2.is_finite() && r4.is_finite()) || d1 <= 0.0 || d2 <= d1 || d2 > 4.0 * d1 * round {
        return false;
    }
    let p = (d2 / d1).log2().max(1.0);
    r <= round * d1.max(1.25 * d1 / (2f64.powf(p) - 1.0))
}

/// Of a cell's samples (their `chords`, `sample_chords`; `occ`, their
/// material), per phase (material, air): those holding a thin feature of
/// it, and those in which it is resolved: the centre in it, and the
/// feature at least THIN cells across.
pub(crate) fn thin_and_resolved(occ: &[(f64, V3, Cov)], chords: &[[Option<f64>; 2]], size: f64) -> [(u64, u64); 2] {
    let mut out = [(0, 0); 2];
    for (i, c) in chords.iter().enumerate() {
        for (phase, (thin, resolved)) in out.iter_mut().enumerate() {
            let centre = (occ[i].0 >= 0.5) == (phase == 0);
            match c[phase] {
                Some(c) if c < THIN * size => *thin |= 1 << i,
                Some(_) if centre => *resolved |= 1 << i,
                _ => {}
            }
        }
    }
    out
}

/// Refinement: a cell holding a thin feature (`THIN`) is split into 8
/// children, sampled again on their own, and a child still holding one
/// again: to this depth at most, 2^depth times finer than the grid's cells
/// (a solve's step 3).
pub(crate) const MAX_DEPTH: u8 = 3;

/// The cells' nodes (`number`).
pub(crate) struct Numbering {
    pub(crate) cell_nodes: Vec<[u32; 8]>,
    /// Each node's place on the lattice (`lattice`).
    pub(crate) node_idx: Vec<[i64; 3]>,
    /// Hanging nodes: a child's corner inside a face or an edge of a grid
    /// cell it touches, and the nodes it is the interpolant of, weighted.
    pub(crate) hang: Vec<(u32, Vec<(u32, f64)>)>,
}

/// Number the nodes of `cells` (sorted by grid cell, then child), whose
/// finest children are 2^depth finer than the grid's cells. A node is
/// shared by the cells round it whose material touches there: full cells
/// always (one copy per lattice point, which they share), a cut cell's
/// piece where its samples touch the other cell's across their common
/// face, edge or corner. So a gap narrower than a cell, which the samples
/// resolve, does not hold its two sides together. A cell's corner inside a
/// face or an edge of a larger cell (a grid cell, or a child less split)
/// whose material its samples touch hangs on that face's (edge's) corners:
/// the larger cell's field there is bilinear (linear) in them, and the
/// smaller's follows it. On the largest such, where several meet there;
/// a corner hung on may itself hang, and a hanging node is given in the
/// nodes that do not.
pub(crate) fn number(cells: &[Classified], dims: [i64; 3], depth: u8) -> Numbering {
    use std::collections::HashMap;
    let s = 1i64 << depth;
    let k = K as i64;
    let in_grid = |i: [i64; 3]| (0..3).all(|a| i[a] >= 0 && i[a] < dims[a]);
    let lin = |i: [i64; 3]| grid::index(i.map(|v| v as usize), dims.map(|v| v as usize));
    let at_cell = |idx: [i64; 3]| -> std::ops::Range<usize> {
        let l = lin(idx);
        let a = cells.partition_point(|c| lin(c.idx) < l);
        a..a + cells[a..].iter().take_while(|c| c.idx == idx).count()
    };
    let bits = |j: usize| -> [i64; 3] { [(j & 1) as i64, ((j >> 1) & 1) as i64, ((j >> 2) & 1) as i64] };
    let corner = |ci: usize, j: usize| -> [i64; 3] {
        let (lo, side) = lattice(&cells[ci], depth);
        let b = bits(j);
        std::array::from_fn(|a| lo[a] + b[a] * side)
    };
    // Children by their low corner: (first, end) of their pieces.
    let mut child: HashMap<[i64; 3], (usize, usize)> = HashMap::new();
    for (i, c) in cells.iter().enumerate().filter(|(_, c)| c.level > 0) {
        child.entry(lattice(c, depth).0).and_modify(|r| r.1 = i + 1).or_insert((i, i + 1));
    }
    // Union-find over the grid's points (full cells' copies), full
    // children's other lattice points, and each cut piece's 8 corners.
    let nd = [dims[0] + 1, dims[1] + 1, dims[2] + 1];
    let nlin = |i: [i64; 3]| grid::index(i.map(|v| (v / s) as usize), nd.map(|v| v as usize));
    let mut nslots = (nd[0] * nd[1] * nd[2]) as usize;
    let mut cut_base = vec![u32::MAX; cells.len()];
    let mut off_grid: HashMap<[i64; 3], usize> = HashMap::new();
    for (ci, c) in cells.iter().enumerate() {
        if c.moments.is_some() {
            cut_base[ci] = nslots as u32;
            nslots += 8;
        } else if c.level > 0 {
            for j in 0..8 {
                let p = corner(ci, j);
                if p.iter().any(|v| v % s != 0) {
                    off_grid.entry(p).or_insert_with(|| {
                        nslots += 1;
                        nslots - 1
                    });
                }
            }
        }
    }
    let slot = |ci: usize, j: usize| -> usize {
        if cut_base[ci] != u32::MAX {
            return cut_base[ci] as usize + j;
        }
        let p = corner(ci, j);
        if p.iter().all(|v| v % s == 0) { nlin(p) } else { off_grid[&p] }
    };
    let mut parent: Vec<u32> = (0..nslots as u32).collect();
    let mut join = |a: usize, b: usize| {
        let (ra, rb) = (find(&mut parent, a as u32), find(&mut parent, b as u32));
        parent[ra as usize] = rb;
    };
    // Join the corners of `ci` and `cj` that coincide.
    let mut hanging: Vec<(usize, i64, Vec<(usize, f64)>)> = Vec::new();
    let offsets = (0..27).map(|d| [d % 3 - 1, (d / 3) % 3 - 1, d / 9 - 1]).filter(|o: &[i64; 3]| *o != [0, 0, 0]);
    for (ci, c) in cells.iter().enumerate() {
        let cut = c.moments.is_some();
        if c.level == 0 && !cut {
            continue;
        }
        let (lo, side) = lattice(c, depth);
        for off in offsets.clone() {
            let nlo: [i64; 3] = std::array::from_fn(|a| lo[a] + off[a] * side);
            // The cell beside, of this size or larger: a child holding the
            // place, else a grid cell. Smaller cells there see this one.
            let mut found = None;
            let mut up = side;
            while up < s {
                let at: [i64; 3] = std::array::from_fn(|a| nlo[a].div_euclid(up) * up);
                // A child there: holding the place if at least this
                // size, else the place is split finer.
                if let Some(&(a, b)) = child.get(&at) {
                    if lattice(&cells[a], depth).1 >= up {
                        found = Some(a..b);
                    }
                    break;
                }
                up *= 2;
            }
            let grid_cell: [i64; 3] = std::array::from_fn(|a| nlo[a].div_euclid(s));
            let others = match found {
                Some(r) => r,
                None if up < s || !in_grid(grid_cell) => continue,
                None => {
                    let r = at_cell(grid_cell);
                    // A refined grid cell's child there holds no material.
                    if r.is_empty() || cells[r.start].level > 0 {
                        continue;
                    }
                    r
                }
            };
            for cj in others {
                let o = &cells[cj];
                let (olo, oside) = lattice(o, depth);
                if oside == side {
                    // One size: each pair once, from a cut one.
                    if !cut || o.moments.is_some() && cj < ci || !touch(c.mask, o.mask, off) {
                        continue;
                    }
                    for j in 0..8 {
                        let g = corner(ci, j);
                        let d: [i64; 3] = std::array::from_fn(|a| g[a] - olo[a]);
                        if d.iter().all(|&v| v == 0 || v == side) {
                            let j2 = (0..3).map(|a| ((d[a] / side) as usize) << a).sum();
                            join(slot(ci, j), slot(cj, j2));
                        }
                    }
                    continue;
                }
                // A cell beside a larger one: its samples against the
                // larger's in its place there.
                if !touch(c.mask, window(o.mask, olo, oside, nlo.map(|v| v * k), side), off) {
                    continue;
                }
                for j in 0..8 {
                    let g = corner(ci, j);
                    if (0..3).any(|a| g[a] != nlo[a] && g[a] != nlo[a] + side) {
                        continue;
                    }
                    let rel: [i64; 3] = std::array::from_fn(|a| g[a] - olo[a]);
                    let free: Vec<usize> = (0..3).filter(|&a| rel[a] != 0 && rel[a] != oside).collect();
                    let at_corner = |b: [i64; 3]| (0..3).map(|a| (b[a] as usize) << a).sum::<usize>();
                    if free.is_empty() {
                        join(slot(ci, j), slot(cj, at_corner(rel.map(|v| v / oside))));
                        continue;
                    }
                    // Inside the grid cell's edge (one axis free) or face
                    // (two): the corners spanning it, bilinearly.
                    let masters: Vec<(usize, f64)> = (0..1usize << free.len())
                        .map(|m| {
                            let mut b = rel.map(|v| v / oside);
                            let mut w = 1.0;
                            for (q, &a) in free.iter().enumerate() {
                                let t = rel[a] as f64 / oside as f64;
                                b[a] = (m >> q & 1) as i64;
                                w *= if b[a] == 1 { t } else { 1.0 - t };
                            }
                            (slot(cj, at_corner(b)), w)
                        })
                        .collect();
                    hanging.push((slot(ci, j), oside, masters));
                }
            }
        }
    }
    // A hanging point hangs on the largest cells it touches; on an edge
    // between such cells it hangs on each's: their corners spanning it
    // are joined.
    for h in hanging.iter_mut() {
        h.0 = find(&mut parent, h.0 as u32) as usize;
    }
    hanging.sort_by_key(|h| (h.0, std::cmp::Reverse(h.1)));
    hanging.dedup_by(|b, a| a.0 == b.0 && a.1 > b.1);
    let mut groups: Vec<(usize, Vec<(usize, f64)>)> = Vec::new();
    for (sl, _, masters) in hanging {
        let r = find(&mut parent, sl as u32) as usize;
        match groups.last() {
            Some(g) if find(&mut parent, g.0 as u32) as usize == r => {
                let first: Vec<usize> = g.1.iter().map(|m| m.0).collect();
                for (a, b) in first.into_iter().zip(masters.iter().map(|m| m.0)) {
                    let (ra, rb) = (find(&mut parent, a as u32), find(&mut parent, b as u32));
                    parent[ra as usize] = rb;
                }
            }
            _ => groups.push((r, masters)),
        }
    }
    let mut node_of = vec![u32::MAX; nslots];
    let mut node_idx: Vec<[i64; 3]> = Vec::new();
    let mut cell_nodes: Vec<[u32; 8]> = vec![[0; 8]; cells.len()];
    for ci in 0..cells.len() {
        for j in 0..8 {
            let r = find(&mut parent, slot(ci, j) as u32) as usize;
            if node_of[r] == u32::MAX {
                node_of[r] = node_idx.len() as u32;
                node_idx.push(corner(ci, j));
            }
            cell_nodes[ci][j] = node_of[r];
        }
    }
    let direct: HashMap<u32, Vec<(u32, f64)>> = groups
        .into_iter()
        .map(|(r, masters)| {
            let n = node_of[find(&mut parent, r as u32) as usize];
            (n, masters.into_iter().map(|(m, w)| (node_of[find(&mut parent, m as u32) as usize], w)).collect())
        })
        .collect();
    // In the nodes that do not hang: a corner hung on lies on a larger
    // cell's face or edge, so the chain ends.
    fn free(n: u32, w: f64, direct: &HashMap<u32, Vec<(u32, f64)>>, out: &mut Vec<(u32, f64)>) {
        match direct.get(&n) {
            Some(m) => m.iter().for_each(|&(m, v)| free(m, w * v, direct, out)),
            None => out.push((n, w)),
        }
    }
    let mut hang: Vec<(u32, Vec<(u32, f64)>)> = direct
        .keys()
        .map(|&n| {
            let mut out = Vec::new();
            free(n, 1.0, &direct, &mut out);
            out.sort_unstable_by_key(|m| m.0);
            out.dedup_by(|b, a| {
                let same = a.0 == b.0;
                if same {
                    a.1 += b.1;
                }
                same
            });
            (n, out)
        })
        .collect();
    hang.sort_unstable_by_key(|h| h.0);
    Numbering { cell_nodes, node_idx, hang }
}

/// The samples within one step of `m` (on the cell's K^3 sample grid, the
/// 26 neighbours), in the cell.
fn grow(m: u64) -> u64 {
    // Per axis: the bits not on the high (low) face of that axis, whose
    // shift along it stays in the cell.
    const NOT_HI: [u64; 3] = [0x7777_7777_7777_7777, 0x0FFF_0FFF_0FFF_0FFF, 0x0000_FFFF_FFFF_FFFF];
    const NOT_LO: [u64; 3] = [0xEEEE_EEEE_EEEE_EEEE, 0xFFF0_FFF0_FFF0_FFF0, 0xFFFF_FFFF_FFFF_0000];
    let mut m = m;
    for (a, step) in [1u32, 4, 16].into_iter().enumerate() {
        m |= ((m & NOT_HI[a]) << step) | ((m & NOT_LO[a]) >> step);
    }
    m
}

/// The pieces of a cut cell's material: its samples holding material,
/// grouped by touching (26 neighbours on the sample grid).
pub(crate) fn pieces(material: u64) -> Vec<u64> {
    let (mut left, mut out) = (material, Vec::new());
    while left != 0 {
        let mut piece = left & left.wrapping_neg();
        loop {
            let g = grow(piece) & material;
            if g == piece {
                break;
            }
            piece = g;
        }
        out.push(piece);
        left &= !piece;
    }
    out
}

/// Whether material samples `a` of a cell and `b` of its neighbour at
/// offset `off` (each axis -1, 0 or 1) touch: some sample of each within
/// one step of the other on the grid of samples the two cells share.
pub(crate) fn touch(a: u64, b: u64, off: [i64; 3]) -> bool {
    let k = K as i64;
    let at = |i: i64| [i % k, (i / k) % k, i / (k * k)];
    let mut bits = a;
    while bits != 0 {
        let i = bits.trailing_zeros() as i64;
        bits &= bits - 1;
        let p = at(i);
        // The sample must lie on the faces toward the neighbour.
        if (0..3).any(|x| off[x] == 1 && p[x] != k - 1 || off[x] == -1 && p[x] != 0) {
            continue;
        }
        // b's samples beside it: across a face the facing layer, along a
        // shared axis within one step.
        let range = |x: usize| -> (i64, i64) {
            match off[x] {
                1 => (0, 0),
                -1 => (k - 1, k - 1),
                _ => ((p[x] - 1).max(0), (p[x] + 1).min(k - 1)),
            }
        };
        let (rx, ry, rz) = (range(0), range(1), range(2));
        for z in rz.0..=rz.1 {
            for y in ry.0..=ry.1 {
                for x in rx.0..=rx.1 {
                    if b >> (x + k * y + k * k * z) & 1 == 1 {
                        return true;
                    }
                }
            }
        }
    }
    false
}

/// Of a cell's pieces, by their samples' masks, the one with the sample
/// nearest `xi`.
pub(crate) fn nearest_piece(masks: impl Iterator<Item = u64>, xi: V3) -> usize {
    let k = K as i64;
    let q: [i64; 3] = std::array::from_fn(|a| (((xi[a] + 1.0) * 0.5 * K as f64).floor() as i64).clamp(0, k - 1));
    let mut best = (0, i64::MAX);
    for (pi, mut bits) in masks.enumerate() {
        while bits != 0 {
            let i = bits.trailing_zeros() as i64;
            bits &= bits - 1;
            let p = [i % k, (i / k) % k, i / (k * k)];
            let d: i64 = (0..3).map(|a| (p[a] - q[a]).pow(2)).sum();
            if d < best.1 {
                best = (pi, d);
            }
        }
    }
    best.0
}

/// Sample weights of a cut cell kept to its piece's samples.
pub(crate) fn masked(mut w: Vec<f64>, mask: u64) -> Vec<f64> {
    if mask != ALL {
        for (i, v) in w.iter_mut().enumerate() {
            if mask >> i & 1 == 0 {
                *v = 0.0;
            }
        }
    }
    w
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub(crate) struct SimCell {
    pub(crate) idx: [i64; 3],
    pub(crate) nodes: [u32; 8],
    /// Material volume fraction.
    pub(crate) frac: f64,
    /// Mean coefficient over the material part.
    pub(crate) coef: f64,
    /// Index into the cut-cell matrices, or u32::MAX for a full cell.
    pub(crate) cut: u32,
    /// The shape's tape simplified for this cell (cut cells of the fine grid).
    #[serde(skip)]
    pub(crate) tape: Option<Arc<Tape>>,
    /// Stiffness multiplier (a redesign's, `Assembly::solve`); 1 otherwise.
    pub(crate) mult: f64,
    /// Elastic: Poisson's ratio (the cell's mean over its material).
    pub(crate) nu: f64,
    /// Index of the full-cell matrix for this nu (see System::k_ref).
    #[serde(skip)]
    pub(crate) mat: u32,
    /// Cut cells: moments of the material (not coefficient-weighted), for
    /// integrals over it other than the stiffness (strains at points).
    pub(crate) material: Option<Box<[f64; 27]>>,
    /// The samples (`sample_xi(K)`) holding this cell's piece of material:
    /// a cut cell whose material falls apart into pieces is one cell per
    /// piece, each with its own nodes. ALL for full cells.
    #[serde(default = "all_samples")]
    pub(crate) mask: u64,
    /// A child of grid cell `idx` (`MAX_DEPTH`): 2^level times finer, at
    /// `sub` (in its own cells) within it; level 0 for a grid cell.
    #[serde(default)]
    pub(crate) level: u8,
    #[serde(default)]
    pub(crate) sub: [u16; 3],
    /// Two materials across a face in it: its element their laminate's
    /// (`laminate`), each one's stress its own.
    #[serde(default)]
    pub(crate) lam: Option<Box<Laminate>>,
}

fn all_samples() -> u64 {
    ALL
}

/// A cell, a grid cell or a child of one: where it is and what it holds.
pub(crate) trait Placed {
    fn idx(&self) -> [i64; 3];
    fn level(&self) -> u8;
    fn sub(&self) -> [u16; 3];
    fn mask(&self) -> u64;
}

impl Placed for SimCell {
    fn idx(&self) -> [i64; 3] {
        self.idx
    }
    fn level(&self) -> u8 {
        self.level
    }
    fn sub(&self) -> [u16; 3] {
        self.sub
    }
    fn mask(&self) -> u64 {
        self.mask
    }
}

impl Placed for Classified {
    fn idx(&self) -> [i64; 3] {
        self.idx
    }
    fn level(&self) -> u8 {
        self.level
    }
    fn sub(&self) -> [u16; 3] {
        self.sub
    }
    fn mask(&self) -> u64 {
        self.mask
    }
}

/// A cell's centre and side (mm).
pub(crate) fn place(grid: &Grid, c: &impl Placed) -> (V3, f64) {
    let size = grid.h / (1u64 << c.level()) as f64;
    let lo = grid.pos(c.idx());
    (std::array::from_fn(|a| lo[a] + size * (c.sub()[a] as f64 + 0.5)), size)
}

/// A cell's low corner and side on the lattice of nodes, whose step is
/// the grid's cell over 2^depth.
pub(crate) fn lattice(c: &impl Placed, depth: u8) -> ([i64; 3], i64) {
    let side = 1i64 << (depth - c.level());
    (std::array::from_fn(|a| (c.idx()[a] << depth) + c.sub()[a] as i64 * side), side)
}

/// Kept leaves sharing a face of nonzero area with a dropped cell.
pub(crate) fn face_neighbours<'a>(cells: &'a [Classified], cell: &Classified, at: impl Fn([i64; 3]) -> std::ops::Range<usize> + 'a) -> impl Iterator<Item = usize> + 'a {
    let (lo, side) = lattice(cell, MAX_DEPTH);
    let hi = lo.map(|v| v + side);
    let parent_side = 1i64 << MAX_DEPTH;
    (0..6).flat_map(move |face| {
        let (axis, upper) = (face / 2, face % 2 == 1);
        let plane = if upper { hi[axis] } else { lo[axis] };
        let mut parent = lo.map(|v| v.div_euclid(parent_side));
        parent[axis] = (if upper { plane } else { plane - 1 }).div_euclid(parent_side);
        at(parent).filter(move |&i| {
            let (nb, size) = lattice(&cells[i], MAX_DEPTH);
            let on_plane = if upper { nb[axis] == plane } else { nb[axis] + size == plane };
            on_plane && (0..3).all(|a| a == axis || (nb[a] < hi[a] && lo[a] < nb[a] + size))
        })
    })
}

/// A cell's material samples (`mask`, the cell at lattice `lo`, `side`)
/// as seen on 4^3 samples of side `unit` from `at`, both in sample units
/// of the finest cells (a lattice step over K), `unit` no coarser than the
/// cell's own samples: each sample holds material where the cell's sample
/// it lies in does; none outside the cell.
fn window(mask: u64, lo: [i64; 3], side: i64, at: [i64; 3], unit: i64) -> u64 {
    let k = K as i64;
    let mut out = 0u64;
    for i in 0..64i64 {
        let p = [i % k, (i / k) % k, i / (k * k)];
        let s: [i64; 3] = std::array::from_fn(|a| (at[a] + p[a] * unit - lo[a] * k).div_euclid(side));
        if s.iter().all(|&v| (0..k).contains(&v)) && mask >> (s[0] + k * s[1] + k * k * s[2]) & 1 == 1 {
            out |= 1 << i;
        }
    }
    out
}

/// Of the cells `cells[range]` (one grid cell's: its pieces, or its
/// children's), the one holding point `q` and `q` in its coordinates: the
/// child it lies in, of its pieces the one whose samples are nearest it.
/// None where `q`'s child holds no material.
pub(crate) fn locate<T: Placed>(cells: &[T], range: std::ops::Range<usize>, grid: &Grid, q: V3) -> Option<(usize, V3)> {
    let first = cells.get(range.start)?;
    let lo = grid.pos(first.idx());
    // A grid cell's children may be of several sizes: each size's place
    // holding q, where a child of that size is.
    let mut levels: Vec<u8> = range.clone().map(|i| cells[i].level()).collect();
    levels.sort_unstable();
    levels.dedup();
    levels.into_iter().find_map(|level| {
        let size = grid.h / (1u64 << level) as f64;
        let n = 1i64 << level;
        let s0: [i64; 3] = std::array::from_fn(|a| (((q[a] - lo[a]) / size).floor() as i64).clamp(0, n - 1));
        // On a plane between children, the one below too.
        let below: [bool; 3] = std::array::from_fn(|a| s0[a] > 0 && q[a] - lo[a] == size * s0[a] as f64);
        (0..8).filter(|&m: &usize| (0..3).all(|a| m >> a & 1 == 0 || below[a])).find_map(|m| {
            let sub: [u16; 3] = std::array::from_fn(|a| (s0[a] - (m >> a & 1) as i64) as u16);
            let xi: V3 = std::array::from_fn(|a| (q[a] - lo[a] - size * (sub[a] as f64 + 0.5)) * 2.0 / size);
            let at: Vec<usize> = range.clone().filter(|&i| cells[i].level() == level && cells[i].sub() == sub).collect();
            let near = nearest_piece(at.iter().map(|&i| cells[i].mask()), xi.map(|v| v.clamp(-1.0, 1.0)));
            at.get(near).map(|&i| (i, xi))
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sample bits from (a, b, c) on the K^3 grid.
    fn bits(at: &[[usize; 3]]) -> u64 {
        at.iter().fold(0, |m, p| m | 1 << (p[0] + K * p[1] + K * K * p[2]))
    }

    #[test]
    fn pieces_are_the_touching_samples() {
        // Two layers two apart along z: two pieces; one apart: one.
        let layer = |c: usize| bits(&(0..16).map(|i| [i % 4, i / 4, c]).collect::<Vec<_>>());
        assert_eq!(pieces(layer(0) | layer(2)), vec![layer(0), layer(2)]);
        assert_eq!(pieces(layer(1) | layer(2)).len(), 1);
        // A diagonal chain touches corner to corner (26 neighbours), and
        // nothing wraps round the cell's faces.
        assert_eq!(pieces(bits(&[[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]])).len(), 1);
        assert_eq!(pieces(bits(&[[3, 0, 0], [0, 1, 0]])).len(), 2);
        assert_eq!(pieces(bits(&[[0, 3, 0], [0, 0, 1]])).len(), 2);
        assert_eq!(grow(bits(&[[0, 0, 0]])), bits(&[[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]));
    }

    #[test]
    fn neighbours_touch_where_their_samples_do() {
        // Across the +x face: a's samples at a = 3 against b's at a = 0,
        // within one step in y and z.
        let (a, b) = (bits(&[[3, 1, 1]]), bits(&[[0, 2, 2]]));
        assert!(touch(a, b, [1, 0, 0]) && touch(b, a, [-1, 0, 0]));
        assert!(!touch(a, bits(&[[0, 3, 1]]), [1, 0, 0]));
        assert!(!touch(a, bits(&[[1, 1, 1]]), [1, 0, 0]));
        // Across an edge and a corner.
        assert!(touch(bits(&[[3, 3, 2]]), bits(&[[0, 0, 1]]), [1, 1, 0]));
        assert!(!touch(bits(&[[3, 2, 2]]), bits(&[[0, 0, 1]]), [1, 1, 0]));
        assert!(touch(bits(&[[3, 3, 3]]), bits(&[[0, 0, 0]]), [1, 1, 1]));
        assert!(touch(ALL, ALL, [-1, 1, 0]));
    }

    #[test]
    fn a_window_sees_a_cell_on_finer_samples() {
        // A grid cell (side 2 lattice steps, samples of 2 units) holding
        // its low-x half: the window of a child-sized box at its low
        // corner holds it all; one at its high-x half, none; one half
        // outside it, its own half.
        let half = bits(&(0..64).map(|i| [i % 4, (i / 4) % 4, i / 16]).filter(|p| p[0] < 2).collect::<Vec<_>>());
        assert_eq!(window(half, [0; 3], 2, [0; 3], 1), ALL);
        assert_eq!(window(half, [0; 3], 2, [4, 0, 0], 1), 0);
        assert_eq!(window(ALL, [0; 3], 2, [-2, 0, 0], 1), bits(&(0..64).map(|i| [i % 4, (i / 4) % 4, i / 16]).filter(|p| p[0] >= 2).collect::<Vec<_>>()));
    }
}
