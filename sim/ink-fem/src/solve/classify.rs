//! The solve's cells: classified on the grid (`classify`), then split where
//! they hold thin features (`refine`).

use crate::cells::{AIR, ALL, K, MATERIAL, MAX_DEPTH, MIN_FRACTION, PROPERTY, THIN, place, sample_chords};
use crate::element::sample_xi;
use crate::eval::{eval_interval, eval_points};
use crate::occupancy::{Classified, SECANT_STEPS, at_xi, center};
use crate::octree::{Cell, descend};
use crate::tape::Tape;
use crate::{Boundary, RegionSel};
use crate::{memory, slivers};
use ink_core::{V3, eval, interval, v3};
use ink_kernel::mesh::find;
use rayon::prelude::*;
use std::sync::Arc;

use super::SolveError;
use super::study::Study;

/// The grid's cells classified: full blocks from the octree expanded into
/// unit cells, cut cells sampled into their pieces, sorted (`Study::order`);
/// whether the memory warning was said, and material the samples missed
/// that changes how the part connects (`slivers`).
#[allow(clippy::type_complexity)]
pub(super) fn classify(st: &Study) -> (Vec<Classified>, bool, Option<([i64; 3], V3, String)>) {
    let Study { s, ctx, h, ref grid, root, .. } = *st;
    // 1. Classify cells: full blocks from the octree, cut cells sampled.
    enum Out {
        Block(Cell),
        /// A cut cell's pieces.
        Cut(Vec<Classified>),
    }
    // Material between the samples (`slivers`): each cell whose samples see
    // none has its lines between them walked, and where it found some.
    // none has its lines between them walked once every cell's samples are
    // known (a line beside a sample inside the part is that sample's face).
    let to_walk: std::sync::Mutex<Vec<([i64; 3], V3, Arc<Tape>)>> = std::sync::Mutex::new(Vec::new());
    let masks: std::sync::Mutex<std::collections::HashMap<[i64; 3], u64>> = std::sync::Mutex::new(std::collections::HashMap::new());
    let out: Vec<Out> = descend(
        &s.shape,
        ctx,
        grid,
        root,
        1,
        &|t: &Arc<Tape>, c| {
            if !st.in_grid(c.i) {
                return None;
            }
            let ctr = center(grid, c.i);
            let (pts, occ, material) = st.sampled(t, ctr, h);
            if material == 0 {
                to_walk.lock().unwrap().push((c.i, ctr, t.clone()));
                return None;
            }
            masks.lock().unwrap().insert(c.i, material);
            // Every sample whole: a full cell. Intervals cannot prove one
            // whose side lies on the part's face full (the field reaches 0
            // there), and as a cut cell its samples' moments would only
            // approximate the full cell's exact ones.
            if occ.iter().all(|&(o, _, _)| o == 1.0) {
                return Some(Out::Block(c));
            }
            let mut kept = st.pieces_of(t, c.i, 0, [0; 3], h, &pts, &occ, material);
            let (thin, resolved) = st.property(ctr, h, 0);
            kept.iter_mut().for_each(|k| (k.thin, k.resolved) = (k.thin | thin, k.resolved | resolved));
            Some(Out::Cut(kept))
        },
        &|_: &Arc<Tape>, c| Some(Out::Block(c)),
    );

    // (a cell the samples see material in: one with a sample inside the
    // part, or in a full block)
    let masks = masks.into_inner().unwrap();
    let full: std::collections::HashMap<[i64; 3], i64> = out.iter().filter_map(|o| if let Out::Block(b) = o { Some((b.i, b.size)) } else { None }).collect();
    let sizes: Vec<i64> = full.values().copied().collect::<std::collections::BTreeSet<_>>().into_iter().collect();
    let in_full = |i: [i64; 3]| sizes.iter().any(|&z| full.get(&i.map(|v| v - v.rem_euclid(z))).is_some_and(|&b| b == z));
    let seen = |i: [i64; 3]| masks.contains_key(&i) || in_full(i);
    let cell_bits = |c: [i64; 3]| masks.get(&c).copied().unwrap_or(if in_full(c) { u64::MAX } else { 0 });
    let lat = slivers::Lattice { o: center(grid, [0; 3]).map(|v| v - 0.5 * h), s: h / 4.0, cell: &cell_bits };
    let to_walk = to_walk.into_inner().unwrap();
    let missed: Vec<([i64; 3], V3)> = to_walk.par_iter().filter_map(|(i, ctr, t)| slivers::missed_in_cell(t, ctx, *i, *ctr, h, &lat).map(|at| (*i, at))).collect();
    let missed_at = slivers::changes_connection(missed, &seen);
    // 2. Expand full blocks into unit cells.
    let mut classified: Vec<Classified> = Vec::new();
    let mut blocks = Vec::new();
    for o in out {
        match o {
            Out::Cut(c) => classified.extend(c),
            Out::Block(b) => blocks.push(b),
        }
    }
    // Said on stderr where its cells may not fit in memory (`memory`):
    // here, before full blocks become cells, or once cells are split.
    let cut = classified.iter().filter(|c| c.frac >= MIN_FRACTION).count() as f64;
    let said = memory::warn(s, h, memory::Cells { full: blocks.iter().map(|b| b.size.pow(3) as f64).sum(), cut, surface: cut }, false);
    let full: Vec<Classified> = blocks
        .par_iter()
        .flat_map_iter(|b| {
            let mut v = Vec::new();
            for z in 0..b.size {
                for y in 0..b.size {
                    for x in 0..b.size {
                        let idx = [b.i[0] + x, b.i[1] + y, b.i[2] + z];
                        if st.in_grid(idx) {
                            v.push(idx);
                        }
                    }
                }
            }
            let pts: Vec<V3> = v.iter().map(|&i| center(grid, i)).collect();
            let coefs: Vec<f64> = match s.coef_const {
                Some(c) => vec![c; v.len()],
                None => st.coef_at(&pts),
            };
            let nus = st.nu_at(&pts);
            let measured = st.full_cells(&pts, h);
            v.into_iter()
                .zip(coefs)
                .zip(nus)
                .zip(measured)
                .zip(&pts)
                .map(|((((idx, coef), nu), (thin, resolved, thinnest)), &ctr)| {
                    let (pt, pr) = st.property(ctr, h, 0);
                    (idx, coef, nu, thin | pt, resolved | pr, thinnest)
                })
                .map(|(idx, coef, nu, thin, resolved, thinnest)| Classified {
                    idx,
                    frac: 1.0,
                    coef,
                    nu,
                    moments: None,
                    tape: None,
                    mask: ALL,
                    thin,
                    resolved,
                    level: 0,
                    sub: [0; 3],
                    thinnest,
                    lam: None,
                })
                .collect::<Vec<_>>()
        })
        .collect();
    classified.extend(full);
    classified.par_sort_unstable_by_key(|c| st.order(c));
    (classified, said, missed_at)
}

/// A cell's key: its grid cell, level and place among its level's children.
pub(super) fn key(c: &Classified) -> ([i64; 3], u8, [u16; 3]) {
    (c.idx, c.level, c.sub)
}

// 3. Cells split where they hold thin features, and nodes (`number`).
// A cell holding material or air thinner than THIN cells across
// (`THIN`: a fin, a web, a ligament, a sheet; a slit, a clearance) is
// split into 8 children, sampled again on their own, and a child still
// holding one again, to MAX_DEPTH or the children's budget: an octree,
// fine only where the feature is thin (a taper's edge splits along the
// edge alone). Nodes are shared where the cells' samples touch, and a
// corner inside a face or an edge of a larger cell it touches hangs on
// that face's (or edge's) corners. A thin feature left at the end is
// reported.
/// The cells once split, and whether splitting stopped short of the
/// floor, by the budget. Refused where a fix holds a slab thinner than the
/// cells resolve.
pub(super) fn refine(st: &Study, classified: Vec<Classified>) -> Result<(Vec<Classified>, bool), SolveError> {
    let Study { s, ctx, h, judge, ref grid, .. } = *st;
    // A cell's children, one level finer, each judged over its own box.
    let children = |pieces: &[Classified]| -> Vec<Classified> {
        let b = &pieces[0];
        let level = b.level + 1;
        let size = h / f64::from(1u32 << level);
        let lo = grid.pos(b.idx);
        let (mut out, mut whole) = (Vec::new(), Vec::new());
        for i in 0..8usize {
            let sub: [u16; 3] = std::array::from_fn(|a| 2 * b.sub[a] + ((i >> a) & 1) as u16);
            let ctr: V3 = std::array::from_fn(|a| lo[a] + size * (sub[a] as f64 + 0.5));
            let Some(t) = &b.tape else {
                whole.push((sub, ctr));
                continue;
            };
            let (pts, occ, material) = st.sampled(t, ctr, size);
            if material == 0 {
                continue;
            }
            if occ.iter().all(|&(o, _, _)| o == 1.0) {
                whole.push((sub, ctr));
                continue;
            }
            let mut kept = st.pieces_of(t, b.idx, level, sub, size, &pts, &occ, material);
            let (thin, resolved) = st.property(ctr, size, level);
            kept.iter_mut().for_each(|k| (k.thin, k.resolved) = (k.thin | thin, k.resolved | resolved));
            out.extend(kept);
        }
        let pts: Vec<V3> = whole.iter().map(|w| w.1).collect();
        let measured = st.full_cells(&pts, size);
        for (((sub, ctr), (coef, nu)), (thin, resolved, thinnest)) in whole.iter().zip(st.coef_at(&pts).into_iter().zip(st.nu_at(&pts))).zip(measured) {
            let (pt, pr) = st.property(*ctr, size, level);
            out.push(Classified {
                idx: b.idx,
                frac: 1.0,
                coef,
                nu,
                moments: None,
                tape: None,
                mask: ALL,
                thin: thin | pt,
                resolved: resolved | pr,
                level,
                sub: *sub,
                thinnest,
                lam: None,
            });
        }
        out
    };
    // The children a solve may add: its grid cells' number, at least 2^18
    // (a solve of a few seconds). A cell is split where it holds a feature
    // thinner than THIN cells of its size, or, for a claim's other grids,
    // of the size the claim's own grid's cell there would be (`split_as`).
    let budget = classified.iter().filter(|c| c.frac >= MIN_FRACTION).count().max(1 << 18);
    // A fix's held region thinner than THIN cells across, inside material:
    // a thin feature as material or air is, split as one (`THIN`). The free
    // material on its two sides shares the nodes of the cells across it and
    // is held as one field (a bar held over a slab 0.366 mm thick at x =
    // 5.0893, on 0.5 mm cells, read 0.007351 in a band to 0.007426 against
    // 0.007458 on 0.125 mm cells). Thin at the floor, and its free material
    // still parted by it in the cells round it, refused (`parted`). Its
    // chord through the region from each held sample, as `sample_chords`
    // reads material's.
    let regions: Vec<(usize, &Arc<Tape>)> = s
        .fixes
        .iter()
        .enumerate()
        // (the fixes that hold a region's volume, as `field` takes them below)
        .filter(|(_, f)| f.film.is_none() && !f.normal && !matches!(f.boundary, Boundary::Mirror { .. }))
        .filter_map(|(fi, f)| if let RegionSel::Field(t) = &f.region { Some((fi, t)) } else { None })
        .collect();
    let sample_at = sample_xi(K);
    // Whether a fix's region parts the free material about `ctr` (a box of
    // `n` samples of side `step` a side) into two pieces: two free samples
    // beside each other joined where the region is proven absent between
    // their centres (its bound over the segment), material samples where
    // both hold material; a piece of material whose free samples fall in two
    // joined sets is parted. (A bound that cannot prove absence parts what
    // may be one: split for nothing, at worst refused at the floor.) And a
    // strip of free material between a held sample and the part's face,
    // thinner than the samples' centres see, is a free piece too: where the
    // line from a held sample's centre to one outside the part beside it
    // leaves the part, the region's field reads free. It joins the free
    // samples beside its own and the strips beside it. (Unseen, a slab 0.08
    // mm thick held 0.03 mm short of a bar's end face was solved on 0.5 mm
    // cells, the lip and the bar one field across it: 0.5% short of the
    // bar's stretch in a band of 0.14%, and 0.76% short with a lip of 0.001.)
    let parted = |fi: usize, ctr: V3, step: f64, n: usize| -> bool {
        let RegionSel::Field(rt) = &s.fixes[fi].region else { return false };
        let pts: Vec<V3> = (0..n * n * n).map(|k| std::array::from_fn(|a| ctr[a] + step * ([k % n, (k / n) % n, k / (n * n)][a] as f64 + 0.5 - 0.5 * n as f64))).collect();
        let (mut v, mut m) = (Vec::new(), Vec::new());
        eval_points(rt, ctx, &pts, &mut v);
        eval_points(&s.shape, ctx, &pts, &mut m);
        let free = |k: usize| m[k] < 0.0 && v[k] > 0.0;
        // The strips: (held sample, the axis and way to the sample outside the part).
        let mut strips: Vec<(usize, usize, bool)> = Vec::new();
        for k in 0..pts.len() {
            for (d, at) in [1, n, n * n].into_iter().enumerate() {
                for up in [true, false] {
                    if (k / at) % n == if up { n - 1 } else { 0 } {
                        continue;
                    }
                    let j = if up { k + at } else { k - at };
                    if !(m[k] < 0.0 && v[k] < 0.0 && m[j] >= 0.0) {
                        continue;
                    }
                    // Where the line leaves the part (false position on the shape's field).
                    let along = |t: f64| -> V3 { std::array::from_fn(|a| pts[k][a] + t * (pts[j][a] - pts[k][a])) };
                    let (mut lo, mut hi, mut flo, mut fhi) = (0.0, 1.0, m[k], m[j]);
                    for _ in 0..SECANT_STEPS {
                        if !(fhi > flo) || fhi == 0.0 {
                            break;
                        }
                        let t = lo - flo * (hi - lo) / (fhi - flo);
                        let f = eval::eval_point(&s.shape, ctx, along(t));
                        if f < 0.0 { (lo, flo) = (t, f) } else { (hi, fhi) = (t, f) }
                    }
                    let face = if fhi > flo { (lo - flo * (hi - lo) / (fhi - flo)).clamp(lo, hi) } else { lo };
                    // (thick enough to be a piece: under the floor of a cut
                    // cell's material, it is held with the region below)
                    if eval::eval_point(rt, ctx, along(face)) > MIN_FRACTION * K as f64 * step {
                        strips.push((k, d, up));
                    }
                }
            }
        }
        let (mut joined, mut free_joined): (Vec<u32>, Vec<u32>) = ((0..pts.len() as u32).collect(), (0..(pts.len() + strips.len()) as u32).collect());
        for k in 0..pts.len() {
            for (d, at) in [1, n, n * n].into_iter().enumerate() {
                let j = k + at;
                if (k / at) % n == n - 1 || !(m[k] < 0.0 && m[j] < 0.0) {
                    continue;
                }
                let (a, b) = (find(&mut joined, k as u32), find(&mut joined, j as u32));
                joined[a as usize] = b;
                let segment: [interval::I; 3] = std::array::from_fn(|a| if a == d { interval::I::new(pts[k][a], pts[j][a]) } else { interval::I::point(pts[k][a]) });
                if free(k) && free(j) && eval_interval(rt, ctx, segment, None).lo > 0.0 {
                    let (a, b) = (find(&mut free_joined, k as u32), find(&mut free_joined, j as u32));
                    free_joined[a as usize] = b;
                }
            }
        }
        // (a strip with the free samples beside its own, and the strips beside it)
        for (i, &(k, d, up)) in strips.iter().enumerate() {
            for at in [1, n, n * n].into_iter().enumerate().filter(|x| x.0 != d).map(|x| x.1) {
                for side in [true, false] {
                    if (k / at) % n == if side { n - 1 } else { 0 } {
                        continue;
                    }
                    let j = if side { k + at } else { k - at };
                    let other = if free(j) { Some(j) } else { strips.iter().position(|x| *x == (j, d, up)).map(|q| pts.len() + q) };
                    if let Some(o) = other {
                        let (a, b) = (find(&mut free_joined, (pts.len() + i) as u32), find(&mut free_joined, o as u32));
                        free_joined[a as usize] = b;
                    }
                }
            }
        }
        let mut piece_of: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
        let pieces: Vec<(usize, usize)> = (0..pts.len()).filter(|&k| free(k)).map(|k| (k, k)).chain(strips.iter().enumerate().map(|(i, x)| (x.0, pts.len() + i))).collect();
        pieces.into_iter().any(|(k, node)| {
            let (m, f) = (find(&mut joined, k as u32), find(&mut free_joined, node as u32));
            *piece_of.entry(m).or_insert(f) != f
        })
    };
    let thin_hold = |c: &Classified| -> Option<usize> {
        let (ctr, size) = place(grid, c);
        let cell: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(ctr[a] - 0.5 * size, ctr[a] + 0.5 * size));
        let pts: Vec<V3> = sample_at.iter().map(|&xi| at_xi(ctr, size, xi)).enumerate().filter(|(k, _)| c.mask >> k & 1 == 1).map(|x| x.1).collect();
        regions.iter().find_map(|&(fi, rt)| {
            let r = eval_interval(rt, ctx, cell, None);
            if !(r.lo < 0.0 && r.hi > 0.0) {
                return None;
            }
            let mut v = Vec::new();
            eval_points(rt, ctx, &pts, &mut v);
            let held: Vec<V3> = pts.iter().zip(&v).filter(|x| *x.1 < 0.0).map(|x| *x.0).collect();
            // (thin, or thinner than its samples, between their centres, and
            // its free material parted by it, over this cell and half of each
            // beside it, as a slab on their shared plane lies: a thin held
            // band at a face, free on one side, is no such feature)
            let thin = held.is_empty() || sample_chords(rt, ctx, &held, size, THIN * size).iter().any(|ch| ch[0].is_some_and(|l| l < THIN * size));
            // (or held at a sample's centre and free of the region, in the
            // part, half a sample either way along its gradient: a slab
            // within the sample, both its faces in it, of which its one
            // tangent plane traces one: free material beside it a sliver
            // with no sample centre, which `parted` cannot see. A bar held
            // over x 0.03 to 0.11 beside a 0.03 sliver at its end, on
            // samples of 0.125, traced the face toward the sliver alone and
            // slipped: 0.20775, past the whole bar's stretch, for 0.1989.)
            let sub = size / K as f64;
            let hair = sub * f64::EPSILON.sqrt();
            let within = held.iter().any(|&p| {
                let probe: Vec<V3> = (0..6).map(|q| std::array::from_fn(|a| p[a] + if a == q / 2 { if q % 2 == 0 { hair } else { -hair } } else { 0.0 })).collect();
                let mut g = Vec::new();
                eval_points(rt, ctx, &probe, &mut g);
                let Some(n) = v3::direction(std::array::from_fn(|a| g[2 * a] - g[2 * a + 1])) else { return false };
                let ends: Vec<V3> = [0.5, -0.5].iter().map(|&t| v3::add(p, v3::scale(n, t * sub))).collect();
                let (mut r, mut at) = (Vec::new(), Vec::new());
                eval_points(rt, ctx, &[p, ends[0], ends[1]], &mut r);
                if !(r[1] > 0.0 && r[2] > 0.0) {
                    return false;
                }
                // (the material between each face, where the field crosses
                // 0 along the line, and its end: the part's face may lie
                // on the end itself)
                let mid: Vec<V3> = (0..2)
                    .map(|e| {
                        let t = r[0] / (r[0] - r[1 + e]);
                        let face = v3::add(p, v3::scale(v3::sub(ends[e], p), t));
                        v3::scale(v3::add(face, ends[e]), 0.5)
                    })
                    .collect();
                eval_points(&s.shape, ctx, &mid, &mut at);
                at.iter().all(|m| *m < 0.0)
            });
            ((thin && parted(fi, ctr, sub, 2 * K)) || within).then_some(fi)
        })
    };
    let wants = |c: &Classified| c.level < MAX_DEPTH && (c.thin & PROPERTY != 0 || c.thinnest.is_some_and(|t| t < THIN * judge / f64::from(1u32 << c.level)));
    let mut cur = classified;
    // A cell beside one of its size (26 neighbours, across grid cells).
    let beside = |c: &Classified| {
        let (level, n) = (c.level, 1i64 << c.level);
        let g: [i64; 3] = std::array::from_fn(|a| c.idx[a] * n + i64::from(c.sub[a]));
        (0..27).filter(|&k| k != 13).map(move |k| {
            let o = [k % 3 - 1, k / 3 % 3 - 1, k / 9 - 1];
            let at: [i64; 3] = std::array::from_fn(|a| g[a] + o[a]);
            (at.map(|v| v.div_euclid(n)), level, at.map(|v| v.rem_euclid(n) as u16))
        })
    };
    // Every cell split so far for what it holds (`wants`).
    let mut wanted: std::collections::HashSet<([i64; 3], u8, [u16; 3])> = std::collections::HashSet::new();
    // Whether splitting stopped short of the floor, by the budget.
    let short = loop {
        // A material step an interface between resolved materials is
        // split no further (`material_held`).
        let held: Vec<usize> = (0..cur.len())
            .into_par_iter()
            .filter(|&i| {
                let (c, size) = place(grid, &cur[i]);
                cur[i].thin & PROPERTY != 0 && st.material_held(c, size, cur[i].level)
            })
            .collect();
        held.into_iter().for_each(|i| cur[i].thin &= !PROPERTY);
        let parts: Vec<bool> = cur.par_iter().map(|c| c.level < MAX_DEPTH && thin_hold(c).is_some()).collect();
        // The cells to split: each one's pieces (contiguous in `cur`).
        let mut split: Vec<std::ops::Range<usize>> = Vec::new();
        let mut i = 0;
        while i < cur.len() {
            let n = cur[i..].iter().take_while(|c| key(c) == key(&cur[i])).count();
            if cur[i..i + n].iter().any(wants) || parts[i..i + n].iter().any(|p| *p) {
                split.push(i..i + n);
            }
            i += n;
        }
        // A sliver is one at its own size: beside a cell of its size split
        // for what it holds, it is split as that one is, its material
        // judged by its children, as a grid of their size would judge it.
        // (Not beside a sliver split so: a run of slivers would be split
        // one more each round, though their children hold nothing.)
        wanted.extend(split.iter().map(|r| key(&cur[r.start])));
        let mut i = 0;
        while i < cur.len() {
            let n = cur[i..].iter().take_while(|c| key(c) == key(&cur[i])).count();
            let cell = &cur[i..i + n];
            if cell[0].level < MAX_DEPTH && cell.iter().any(|c| c.frac < MIN_FRACTION) && !cell.iter().any(wants) && beside(&cell[0]).any(|k| wanted.contains(&k)) {
                split.push(i..i + n);
            }
            i += n;
        }
        split.sort_by_key(|r| r.start);
        if split.is_empty() {
            break false;
        }
        if cur.iter().filter(|c| c.level > 0 && c.frac >= MIN_FRACTION).count() + 8 * split.len() > budget {
            break true;
        }
        let made: Vec<Vec<Classified>> = split.par_iter().map(|r| children(&cur[r.clone()])).collect();
        let gone: std::collections::HashSet<([i64; 3], u8, [u16; 3])> = split.iter().map(|r| key(&cur[r.start])).collect();
        let mut all: Vec<Classified> = cur.into_iter().filter(|c| !gone.contains(&key(c))).collect();
        all.extend(made.into_iter().flatten());
        all.par_sort_unstable_by_key(|c| st.order(c));
        cur = all;
    };
    if let Some((c, fi)) = cur.par_iter().find_map_first(|c| thin_hold(c).map(|fi| (c, fi))) {
        let (ctr, size) = place(grid, c);
        return Err(format!(
            "fix '{}' holds a slab of the part thinner than its cells resolve near {:?} mm: its region leaves free material on both its sides within a cell of {} mm, {}. Write the hold at least two cells thick; where it leaves a sliver of the part free beyond it, extend the hold to the face; where it stands for a plane of symmetry, model one side and declare the cut (mirror name axis coordinate); or refine (res r)",
            s.fixes[fi].name,
            ctr.map(|v| (v * 1000.0).round() / 1000.0),
            size,
            if short { "where the cells' budget stopped their splitting" } else { "the finest they split to" }
        )
        .into());
    }
    Ok((cur, short))
}

// Where a thin feature is left when splitting stops at the floor:
// in cells THIN cells or more from any where its phase is
// resolved. Within that it is a taper's edge (a wedge, or the air
// in a notch's root, is thinner than any cells next to its edge),
// inside the floor's own resolution of the resolved phase, its
// error the grid's, which the claim's coarser grids measure; beyond
// it, a member or a gap no cell this fine resolves. A material not
// one across a cell (`PROPERTY`) still flagged (`material_held`)
// holds one no cell this fine holds. Each such cell's centre, and
// which (the phases left, `PROPERTY` for material).
pub(super) fn left_at_floor(st: &Study, cur: &[Classified]) -> Vec<(V3, u8)> {
    let Study { h, ref grid, .. } = *st;
    let at_lin = |l: usize| -> std::ops::Range<usize> {
        let a = cur.partition_point(|c| st.lin(c.idx) < l);
        a..a + cur[a..].iter().take_while(|c| st.lin(c.idx) == l).count()
    };
    let resolved_at = |p: V3, phase: u8| -> bool {
        let i: [i64; 3] = std::array::from_fn(|a| ((p[a] - grid.origin[a]) / h).floor() as i64);
        st.in_grid(i)
            && at_lin(st.lin(i)).any(|k| {
                let (ctr, size) = place(grid, &cur[k]);
                cur[k].resolved & phase != 0 && (0..3).all(|a| (p[a] - ctr[a]).abs() <= 0.5 * size)
            })
    };
    let reach = THIN as i64;
    let offsets = || (-reach..=reach).flat_map(|x| (-reach..=reach).flat_map(move |y| (-reach..=reach).map(move |z| [x as f64, y as f64, z as f64])));
    cur.par_iter()
        .filter(|c| c.thin != 0)
        .filter_map(|c| {
            let (ctr, size) = place(grid, c);
            let near = |phase: u8| offsets().any(|o| resolved_at(std::array::from_fn(|a| ctr[a] + o[a] * size), phase));
            let left = [MATERIAL, AIR].into_iter().filter(|&p| c.thin & p != 0 && !near(p)).fold(0, |m, p| m | p) | c.thin & PROPERTY;
            (left != 0).then_some((ctr, left))
        })
        .collect()
}
