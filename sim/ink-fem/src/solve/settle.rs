//! The cells once split, settled (`settle`): slivers dropped or taken by
//! the cell beside them, laminates, the nodes, and what the cells leave
//! unresolved.

use crate::Boundary;
use crate::cells::{AIR, ALL, MATERIAL, MIN_FRACTION, Numbering, PROPERTY, PROPERTY_SHRINK, number, place};
use crate::element::moved_moments;
use crate::eval::eval_points;
use crate::laminate::Laminate;
use crate::occupancy::{Classified, at_xi, glued_gap, spread_moments};
use crate::{laminate, memory};
use ink_core::V3;
use rayon::prelude::*;

use super::SolveError;
use super::classify::{key, left_at_floor};
use super::study::{Study, same, width};

/// The solve's cells: classified, split and settled, and their nodes.
pub(super) struct Cells {
    pub(super) classified: Vec<Classified>,
    /// How many times the finest cells were split.
    pub(super) depth: u8,
    pub(super) cell_nodes: Vec<[u32; 8]>,
    pub(super) node_idx: Vec<[i64; 3]>,
    pub(super) hang: Vec<(u32, Vec<(u32, f64)>)>,
    /// The cells, as `memory` counts them.
    pub(super) counted: memory::Cells,
    pub(super) slivers: Slivers,
    pub(super) unresolved: Unresolved,
}

impl Cells {
    /// The cells at a grid index: none, one, or a cut cell's pieces.
    pub(super) fn at_cell(&self, st: &Study, idx: [i64; 3]) -> std::ops::Range<usize> {
        let classified = &self.classified;
        let l = st.lin(idx);
        let a = classified.partition_point(|c| st.lin(c.idx) < l);
        a..a + classified[a..].iter().take_while(|c| c.idx == idx).count()
    }

    /// The cell holding a point of the surface, and the point in its
    /// coordinates: the cell it lies in, or, where that holds no material (a
    /// dropped sliver, the cell beyond a face on its plane), the one half a
    /// cell inside.
    /// (`surface_cell`), of the cells so far.
    pub(super) fn classified_at(&self, st: &Study, idx: [i64; 3]) -> std::ops::Range<usize> {
        if st.in_grid(idx) { self.at_cell(st, idx) } else { 0..0 }
    }
}

/// The slivers (`absorb_slivers`): those taken by the cell beside them, what
/// each such cell gained, and those left dropped, with their volume.
#[allow(clippy::type_complexity)]
pub(super) struct Slivers {
    pub(super) taken: std::collections::HashMap<([i64; 3], u8, [u16; 3]), (([i64; 3], u8, [u16; 3]), [i64; 3])>,
    pub(super) gained: std::collections::HashMap<usize, ([f64; 27], f64)>,
    pub(super) taken_by: std::collections::HashMap<([i64; 3], u8, [u16; 3]), Vec<(Classified, [i64; 3])>>,
    pub(super) left: Vec<Classified>,
    pub(super) volume: f64,
}

/// What the cells leave unresolved, for the report: an interface stepped
/// across cells (`stepped`), a gap the samples glue shut (`glued`), thin
/// features and materials left at the floor, and an interface the grids
/// share (`interface_shared`).
pub(super) struct Unresolved {
    pub(super) stepped: Option<(usize, V3)>,
    pub(super) glued: Option<(usize, V3, f64)>,
    pub(super) thin_left: Option<(usize, V3)>,
    pub(super) material_left: Option<(usize, V3)>,
    pub(super) interface_shared: Option<(usize, V3, f64)>,
}

pub(super) fn settle(st: &Study, cur: Vec<Classified>, short: bool, said: bool) -> Result<Cells, SolveError> {
    let Study { s, h, dims, ref grid, .. } = *st;
    // Slivers (a speck of material in a cell) are dropped, once split: they
    // add nearly-free nodes that ruin conditioning (sharp edges stalled the
    // solver) and carry almost no stiffness. Their volume is reported.
    let (classified, slivers): (Vec<Classified>, Vec<Classified>) = cur.into_iter().partition(|c| c.frac >= MIN_FRACTION);
    let mut classified = classified;
    // A cell holding two materials across a face takes their laminate's
    // element (`laminate`): each one's share of the cell's material (a cut
    // cell's piece's, its samples weighed by what they hold), and the
    // face's normal, measured along lines through its samples. A cell
    // holding more (a gradient, a third material) keeps its own.
    laminates(st, &mut classified);
    // Full cells holding an interface between two materials left whole
    // (`material_held`), each of its centre's value. Where the two lie in
    // them as layers along an axis, every cell of a slab alike (a flat face
    // along the cells), the face is placed on a plane of the cells: the
    // grid's error. Where they lie otherwise (a face oblique to the cells,
    // or curved, or its edge), or a slab's cells place the face apart, the
    // cells along it are a stair of one material or the other, on every
    // grid: a soft band in a pulled bar, turned a degree, read 36 times its
    // stress there, and a heat bar's insulating face, turned 2, 9 K low
    // with a grid error of 3. `stepped`: how many such cells, and one's
    // centre; a claim on the study is undecided
    // (`Report::interface_stepped`).
    let stepped = stepped(st, &classified);
    // A sliver's material stays in the part: the cell beside it that holds
    // the most (of its size, towards the sliver's material) takes it, as
    // its own moments gain the sliver's, moved to its coordinates
    // (`moved_moments`). That cell's shape functions are carried on over
    // the sliver, which has no nodes of its own, so nothing of what
    // dropping it was for comes back; and what presses on the sliver's
    // face goes to that cell, at the place carried on (`taken`).
    // Dropped, a part lost what its faces left in cells under the floor,
    // by where they fell in the cells: a 20 x 2 x 2 bar with sides 0.1 and
    // 0.05 mm past the planes of 0.5 mm cells held 79.9 mm^3 of its 80, and
    // its pressed end's corner, taking the dropped cell's share of the
    // load, moved 1% too far. Not a sliver a fix reaches, nor into a cell
    // one reaches (a hold's own slivers are the hold's: `slivers` below),
    // nor one of several pieces of a cell, nor across two materials: those
    // are dropped as before, and counted.
    let slivers = absorb_slivers(st, &classified, slivers);
    if classified.is_empty() {
        return Err("the study shape has no material inside its bounds".into());
    }
    let depth = classified.iter().map(|c| c.level).max().unwrap_or(0);
    // Stopped short of the floor by the budget, every thin feature is left.
    let left: Vec<(V3, u8)> = if short { classified.iter().filter(|c| c.thin != 0).map(|c| (place(grid, c).0, c.thin)).collect() } else { left_at_floor(st, &classified) };
    let counted = memory::Cells::of(&classified, dims);
    memory::warn(s, h, counted, said);
    // A gap the cells' samples glue shut (`glued_gap`), on the cells the
    // solve uses (split ones part gaps their grid cell glued): how many
    // cells, and the narrowest's middle and width.
    let glued = glued(st, &classified);
    let Numbering { cell_nodes, node_idx, hang, .. } = number(&classified, dims, depth);
    // A claim's other grid (`split_as`) may leave thin features its own
    // grid resolved: its blindness there only widens the grids' difference
    // (`claims::grids_error`), so it is not refused.
    let left_of = |kind: u8| -> Option<(usize, V3)> {
        let of: Vec<V3> = left.iter().filter(|l| l.1 & kind != 0).map(|l| l.0).collect();
        of.first().filter(|_| s.split_as.is_none()).map(|at| (of.len(), at.map(|v| (v * 1000.0).round() / 1000.0)))
    };
    let (thin_left, material_left) = (left_of(MATERIAL | AIR), left_of(PROPERTY));
    // An interface left whole (`material_held`) in a full cell, parallel
    // to one of its walls and between its centre and that wall, the wall
    // on a plane of cells twice its size (a claim's grid and the next,
    // twice as coarse, both have a plane there): a cell takes its centre's
    // value, so both place the face on that plane, its error alike on
    // them, and a claim's error is read from how far those two differ (a
    // band 0.2 mm off planes 2 mm apart read 5.3% low on cells of 0.5, 1
    // and 2 mm, their grid error 1e-5; a face 0.03 mm off a plane of 0.2
    // and 0.4 mm cells, not of 0.8, read 2520 C twice and 2620, a grid
    // error of 0.02 against 7.4). Not on a claim's coarser grids
    // (`split_as` finer than these cells): their own planes are not the
    // ladder's. (A cell of two materials is their laminate, `laminate`,
    // which holds the face where it lies: this is a cell holding a third.)
    let interface_shared = interface_shared(st, &classified);
    Ok(Cells { classified, depth, cell_nodes, node_idx, hang, counted, slivers, unresolved: Unresolved { stepped, glued, thin_left, material_left, interface_shared } })
}

fn laminates(st: &Study, classified: &mut [Classified]) {
    let Study { s, ctx, w_q, ref grid, ref xis, .. } = *st;
    let laminates: Vec<Option<(Box<Laminate>, Option<[f64; 27]>)>> = classified
        .par_iter()
        .map(|c| {
            let (ctr, size) = place(grid, c);
            if !st.stiffness_varies(ctr, size) {
                return None;
            }
            let key = |at: &[V3]| -> Vec<(f64, f64)> {
                let pts: Vec<V3> = at.iter().map(|&xi| at_xi(ctr, size, xi)).collect();
                st.coef_at(&pts).into_iter().zip(st.nu_at(&pts)).collect()
            };
            // A cut cell's piece: its samples' shares, and its material's
            // moments, not weighed by a coefficient.
            let (weight, moments) = match (&c.tape, &c.moments) {
                (Some(t), Some(_)) => {
                    let (_, occ, _) = st.sampled(t, ctr, size);
                    let mut m = [0.0; 27];
                    let mut w = vec![0.0; 64];
                    for (i, &(o, xi, cov)) in occ.iter().enumerate().filter(|(i, _)| c.mask >> i & 1 == 1) {
                        w[i] = o;
                        for (m, v) in m.iter_mut().zip(spread_moments(xi, &cov)) {
                            *m += w_q * o * v;
                        }
                    }
                    (w, Some(m))
                }
                _ => (vec![1.0; 64], None),
            };
            let mut lam = laminate::measure(&key, &same, &weight)?;
            // Each material's yield, where it is one value over it.
            if let Some(y) = &s.strength {
                let pts: Vec<V3> = xis.iter().map(|&xi| at_xi(ctr, size, xi)).collect();
                let mut ys = Vec::new();
                eval_points(y, ctx, &pts, &mut ys);
                let keys = key(xis);
                let one = |a: bool| -> Option<f64> {
                    let of: Vec<f64> = keys.iter().zip(&ys).filter(|(k, _)| lam.is_a(**k, &same) == a).map(|(_, &y)| y).collect();
                    let first = *of.first()?;
                    of.iter().all(|&y| same(y, first)).then_some(first)
                };
                lam.yields = one(true).zip(one(false)).map(|(a, b)| [a, b]);
            }
            Some((Box::new(lam), moments))
        })
        .collect();
    for (c, l) in classified.iter_mut().zip(laminates) {
        if let Some((lam, m)) = l {
            (c.coef, c.nu) = lam.mean();
            if let Some(m) = m {
                c.moments = Some(m.map(|v| v * c.coef));
            }
            c.lam = Some(lam);
        }
    }
}

fn stepped(st: &Study, classified: &[Classified]) -> Option<(usize, V3)> {
    let Study { s, ref grid, ref xis, .. } = *st;
    if s.coef_const.is_some() {
        None
    } else {
        // (cell, and its layers' axis and samples, if in layers.)
        let held: Vec<(usize, Option<(usize, u8)>)> = classified
            .par_iter()
            .enumerate()
            .filter(|(_, c)| c.moments.is_none())
            .filter_map(|(i, c)| {
                let (ctr, size) = place(grid, c);
                if width(&st.hulls(ctr, size)[0]) == 0.0 || !st.material_held(ctr, size, c.level) {
                    return None;
                }
                let pts: Vec<V3> = xis.iter().map(|&xi| at_xi(ctr, size, xi)).collect();
                Some((i, st.layered(&st.coef_at(&pts), ALL)?))
            })
            .collect();
        // The layers of each slab of cells (an axis, a level, a place along
        // the axis): one, where its cells are alike.
        let slab = |i: usize, a: usize| (a, classified[i].level, classified[i].idx[a], classified[i].sub[a]);
        let mut slabs: std::collections::HashMap<(usize, u8, i64, u16), Option<u8>> = std::collections::HashMap::new();
        for (i, l) in &held {
            if let Some((a, bits)) = l {
                slabs.entry(slab(*i, *a)).and_modify(|b| *b = b.filter(|b| b == bits)).or_insert(Some(*bits));
            }
        }
        let off: Vec<usize> = held.iter().filter(|(i, l)| l.is_none_or(|(a, _)| slabs[&slab(*i, a)].is_none())).map(|h| h.0).collect();
        off.first().map(|&i| (off.len(), place(grid, &classified[i]).0.map(|v| (v * 1000.0).round() / 1000.0)))
    }
}

fn absorb_slivers(st: &Study, classified: &[Classified], mut slivers: Vec<Classified>) -> Slivers {
    let Study { s, ctx, ref grid, ref xis, .. } = *st;
    let mut taken: std::collections::HashMap<([i64; 3], u8, [u16; 3]), (([i64; 3], u8, [u16; 3]), [i64; 3])> = std::collections::HashMap::new();
    // (What each cell so gains: the slivers' moments, weighed by their
    // coefficient as a cut cell's own are, in its coordinates, and their
    // share of a cell. Kept beside the cell: a cell's own moments are what
    // its samples hold.)
    let mut gained: std::collections::HashMap<usize, ([f64; 27], f64)> = std::collections::HashMap::new();
    // (And the slivers each took, each with where it lies from it in cells:
    // what is spread through a part's material, its weight, its mass, a
    // source's heat, is sampled in them for it: `sample_cell`.)
    let mut taken_by: std::collections::HashMap<([i64; 3], u8, [u16; 3]), Vec<(Classified, [i64; 3])>> = std::collections::HashMap::new();
    {
        let mut pieces: std::collections::HashMap<([i64; 3], u8, [u16; 3]), u32> = std::collections::HashMap::new();
        for c in classified.iter().chain(&slivers) {
            *pieces.entry(key(c)).or_default() += 1;
        }
        let whole: std::collections::HashMap<([i64; 3], u8, [u16; 3]), usize> =
            classified.iter().enumerate().filter(|(_, c)| pieces[&key(c)] == 1 && c.lam.is_none()).map(|(i, c)| (key(c), i)).collect();
        let reached = |c: &Classified| -> bool {
            let (ctr, size) = place(grid, c);
            let pts: Vec<V3> = xis.iter().map(|&xi| at_xi(ctr, size, xi)).collect();
            // (a fix whose region is material it holds: a mirror's is the
            // half of space its part is in, and a film holds none)
            s.fixes.iter().filter(|f| f.film.is_none() && !matches!(f.boundary, Boundary::Mirror { .. })).any(|f| f.region.occupied(ctx, &pts, size).iter().any(|o| o.0 > 0.0))
        };
        let to: Vec<Option<(usize, [i64; 3])>> = slivers
            .par_iter()
            .map(|sl| {
                // (not a rounding's worth: a cell beyond a face on its plane holds 1e-17 of itself)
                let m = sl.moments.as_ref().filter(|m| m[0] > 1e-12 * 8.0 * sl.coef)?;
                if pieces[&key(sl)] != 1 || reached(sl) {
                    return None;
                }
                let centre = [m[1] / m[0], m[3] / m[0], m[9] / m[0]];
                let n = 1i64 << sl.level;
                let g: [i64; 3] = std::array::from_fn(|a| sl.idx[a] * n + i64::from(sl.sub[a]));
                (0..27i64)
                    .filter(|&k| k != 13)
                    .filter_map(|k| {
                        let o = [k % 3 - 1, k / 3 % 3 - 1, k / 9 - 1];
                        // (towards its material: its centroid's side of the cell)
                        if (0..3).any(|a| o[a] != 0 && o[a] as f64 * centre[a] <= 0.0) {
                            return None;
                        }
                        let at: [i64; 3] = std::array::from_fn(|a| g[a] + o[a]);
                        let i = *whole.get(&(at.map(|v| v.div_euclid(n)), sl.level, at.map(|v| v.rem_euclid(n) as u16)))?;
                        let c = &classified[i];
                        (same(c.coef, sl.coef) && same(c.nu, sl.nu)).then_some((i, o))
                    })
                    .filter(|&(i, _)| !reached(&classified[i]))
                    // (the fullest; of two alike, the one across a face before one across an edge)
                    .max_by(|a, b| {
                        classified[a.0].frac.total_cmp(&classified[b.0].frac).then(b.1.iter().filter(|v| **v != 0).count().cmp(&a.1.iter().filter(|v| **v != 0).count()))
                    })
            })
            .collect();
        let mut left = Vec::new();
        for (sl, to) in std::mem::take(&mut slivers).into_iter().zip(to) {
            let Some((i, o)) = to else {
                left.push(sl);
                continue;
            };
            taken.insert(key(&sl), (key(&classified[i]), o));
            // (the sliver lies at -o from the cell that takes it: two of its coordinates' units a cell)
            let moved = moved_moments(sl.moments.as_ref().expect("a sliver's moments"), o.map(|v| -2.0 * v as f64));
            let gain = gained.entry(i).or_insert(([0.0; 27], 0.0));
            gain.0.iter_mut().zip(moved).for_each(|(g, m)| *g += m);
            gain.1 += sl.frac;
            taken_by.entry(key(&classified[i])).or_default().push((sl, o.map(|v| -v)));
        }
        slivers = left;
    }
    // (The slivers still dropped: a fix's reach, several pieces, another material.)
    let slivers_left = slivers;
    let sliver_volume = slivers_left.iter().map(|c| c.frac * place(grid, c).1.powi(3)).sum();
    Slivers { taken, gained, taken_by, left: slivers_left, volume: sliver_volume }
}

fn glued(st: &Study, classified: &[Classified]) -> Option<(usize, V3, f64)> {
    let Study { s, ctx, ref grid, .. } = *st;
    let glued_cells: Vec<(f64, V3)> = {
        let mut cut: Vec<&Classified> = classified.iter().filter(|c| c.moments.is_some()).collect();
        cut.dedup_by_key(|c| (c.idx, c.level, c.sub));
        cut.par_iter()
            .filter_map(|c| {
                let t = c.tape.as_ref().expect("a cut cell's tape");
                let (ctr, size) = place(grid, *c);
                let (pts, occ, _) = st.sampled(t, ctr, size);
                glued_gap(t, &s.shape, ctx, &pts, &occ, size, ctr)
            })
            .collect()
    };

    glued_cells.iter().min_by(|a, b| a.0.total_cmp(&b.0)).map(|&(w, at)| (glued_cells.len(), at.map(|v| (v * 1000.0).round() / 1000.0), (w * 1000.0).round() / 1000.0))
}

fn interface_shared(st: &Study, classified: &[Classified]) -> Option<(usize, V3, f64)> {
    let Study { s, h, ref grid, ref props, .. } = *st;
    if props.is_empty() || s.split_as.is_some_and(|j| j < h * (1.0 - 1e-9)) {
        None
    } else {
        let at: Vec<(V3, f64)> = classified
            .par_iter()
            .filter(|c| c.moments.is_none() && c.lam.is_none())
            .filter_map(|c| {
                let (ctr, size) = place(grid, c);
                let cut = PROPERTY_SHRINK * size;
                // One material between t0 and t1 cells from its centre
                // along axis `a`, toward side `sg`, and across the cell
                // of the grid twice as coarse that holds it: a face both
                // place alike is flat across their cells (a ring's face,
                // tangent to a plane over one cell, is not over two, and
                // each grid places it otherwise beside).
                let block: V3 = ctr.map(|v| (v / (2.0 * size)).floor() * 2.0 * size);
                let values = |a: usize, sg: f64, t0: f64, t1: f64| -> Option<Vec<f64>> {
                    let (mut lo, mut hi): (V3, V3) = (block.map(|v| v + cut), block.map(|v| v + 2.0 * size - cut));
                    let (p, q) = (ctr[a] + sg * t0 * size, ctr[a] + sg * t1 * size);
                    (lo[a], hi[a]) = (p.min(q), p.max(q));
                    st.one_material(lo, hi)
                };
                let one = |a: usize, sg: f64, t0: f64, t1: f64| values(a, sg, t0, t1).is_some();
                const SLABS: usize = 32;
                (0..3)
                    .any(|a| {
                        [1.0f64, -1.0].into_iter().any(|sg| {
                            let wall = ctr[a] + sg * 0.5 * size;
                            let q = wall / (2.0 * size);
                            if (q - q.round()).abs() > 1e-6 || !one(a, -sg, 0.0, 0.5 - cut / size) || one(a, sg, 0.0, 0.5 - cut / size) {
                                return false;
                            }
                            // Parallel to the wall, and across all of it:
                            // of slabs across that half, the face's (one, or
                            // two beside each other) alone not of one
                            // material, and those between it and the wall
                            // one other material (a curve clipping the
                            // block's corner leaves both at the wall).
                            let slab = |k: usize| (0.5 * k as f64 / SLABS as f64, (0.5 * (k + 1) as f64 / SLABS as f64).min(0.5 - cut / size));
                            let varied: Vec<usize> = (0..SLABS).filter(|&k| !one(a, sg, slab(k).0, slab(k).1)).collect();
                            let (Some(&k0), Some(&k1)) = (varied.first(), varied.last()) else { return false };
                            varied.len() <= 2
                                && k1 - k0 < 2
                                && k1 + 1 < SLABS
                                && values(a, sg, slab(k1 + 1).0, 0.5 - cut / size).zip(values(a, -sg, 0.0, 0.5 - cut / size)).is_some_and(|(beyond, inner)| beyond != inner)
                        })
                    })
                    .then_some((ctr, 2.0 * size))
            })
            .collect();
        at.first().map(|&(p, plane)| (at.len(), p.map(|v| (v * 1000.0).round() / 1000.0), plane))
    }
}
