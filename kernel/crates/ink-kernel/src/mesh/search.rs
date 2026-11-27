use super::brick::{B, Brick, cell_of};
use super::proof::{CROSSED, FOUND, INSIDE, PROVEN, proven};
use super::report::Ran;
use crate::octree::Grid;
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_points};
use ink_core::interval::I;
use ink_core::sign::{drawn, inside};
use ink_core::tape::{Pass, Tape};
use ink_core::v3::norm;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

/// The finest length ink speaks of in `lo..hi` (docs/CORE.md, Range): a
/// nanometre, a billionth of the farthest coordinate farther out.
pub(super) fn finest_length(lo: V3, hi: V3) -> f64 {
    1e-9 * lo.iter().chain(&hi).fold(1e3f64, |m, x| m.max(x.abs()))
}

/// What a search costs per part, in samples: its 27 corners, and the
/// field's range over each of its 8 parts (about what 2 samples cost).
const PART_COST: f64 = 27.0 + 8.0 * 2.0;

/// How many open cells the surface does not cross the search reads (in
/// brick order, `search`): as many as `allowed` samples afford, none where
/// its first split would be finer than `finest`.
pub(super) fn searched_cells(allowed: f64, h: f64, finest: f64) -> usize {
    if 0.5 * h >= finest { (allowed / PART_COST).max(0.0) as usize } else { 0 }
}

/// The bricks' open cells (`Brick::open`) searched at once, level by level:
/// each part's 27 corners sampled, and its parts the field's range does
/// not prove split again, while the level costs (PART_COST a part) within
/// the `allowed` samples left, and its samples are no closer than
/// `finest`; not the cells the surface crosses, where a sample of either
/// sign is drawn. A sample of the other sign finds what the cell holds (FOUND);
/// a cell with no parts left is PROVEN; the others stay undecided. Each
/// brick's tape is the field focused on it. What a sample found: its
/// cell, the sample, whether it is material, and how thick that is there
/// (`across`); the finest spacing sampled (h where none was); and the
/// samples taken.
#[allow(clippy::type_complexity)]
pub(super) fn search(
    bricks: &mut [Brick],
    whole: &Arc<Tape>,
    ctx: &Ctx,
    grid: &Grid,
    allowed: f64,
    finest: f64,
    covered: &(dyn Fn(V3, V3) -> bool + Sync),
    ran: &Ran,
) -> (Vec<([i64; 3], V3, bool, f64)>, f64, f64, HashMap<(u32, u32), f64>) {
    let h = grid.h;
    // A cell the surface crosses is not searched: a sample of either sign
    // is drawn there, so its search could only show it drawn, never find
    // what it leaves out.
    let open: usize = bricks.iter().map(|b| b.open.iter().filter(|o| o.1 & CROSSED == 0).count()).sum();
    let box_of = |i0: [i64; 3], k: u16| {
        let lo = grid.pos(cell_of(i0, k));
        (lo, lo.map(|x| x + h))
    };
    let held_of = |b: &Brick| (grid.pos(b.i0).map(|x| x - 0.5 * h), grid.pos(b.i0).map(|x| x + (B as f64 + 0.5) * h));
    let tape_of = |b: &Brick| {
        let (lo, hi) = held_of(b);
        b.tape.clone().unwrap_or_else(|| {
            ran.tapes.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            crate::query::focus(whole, ctx, lo, hi)
        })
    };
    // The parts left (box, brick, cell's place in its open cells), in brick
    // order, as many as what is left affords; the spacing of the samples
    // the least searched cell left undecided had (h: none but its corners).
    let (mut spacing, mut left, mut samples, mut coarsest) = (h, allowed, 0.0, 0.0f64);
    let affords = |left: f64| (left / PART_COST).max(0.0) as usize;
    let mut parts: Vec<(V3, V3, u32, u32)> = Vec::new();
    let first = searched_cells(allowed, h, finest).min(open);
    'cells: for (bi, b) in bricks.iter_mut().enumerate() {
        let i0 = b.i0;
        for (oi, o) in b.open.iter_mut().enumerate().filter(|o| o.1.1 & CROSSED == 0) {
            if parts.len() == first {
                break 'cells;
            }
            o.1 |= PROVEN;
            let (lo, hi) = box_of(i0, o.0);
            parts.push((lo, hi, bi as u32, oi as u32));
        }
    }
    if first < open {
        coarsest = h;
    }
    // Parts the search stopped at, and the spacing their cell was sampled at.
    let mut stopped: Vec<(u32, u32, f64)> = Vec::new();
    let mut found: Vec<(u32, u32, V3)> = Vec::new();
    while !parts.is_empty() {
        let k = if 0.5 * spacing < finest { 0 } else { affords(left).min(parts.len()) };
        if k < parts.len() {
            coarsest = coarsest.max(spacing);
            stopped.extend(parts.drain(k..).map(|p| (p.2, p.3, spacing)));
        }
        if k == 0 {
            break;
        }
        let cost = PART_COST * k as f64;
        (left, samples, spacing) = (left - cost, samples + cost, 0.5 * spacing);
        let groups: Vec<&[(V3, V3, u32, u32)]> = parts.chunk_by(|a, b| a.2 == b.2).collect();
        let level: Vec<(Vec<(u32, u32, V3)>, Vec<(V3, V3, u32, u32)>)> = groups
            .par_iter()
            .map(|g| {
                let b = &bricks[g[0].2 as usize];
                let tape = tape_of(b);
                let lattice = |lo: V3, hi: V3| (0..27).map(move |q| -> V3 { std::array::from_fn(|a| lo[a] + (hi[a] - lo[a]) * 0.5 * ((q / [1, 3, 9][a]) % 3) as f64) });
                let at: Vec<V3> = g.iter().flat_map(|&(lo, hi, _, _)| lattice(lo, hi)).collect();
                let mut v = Vec::new();
                eval_points(&tape, ctx, &at, &mut v);
                ran.points.fetch_add(at.len() * tape.work(Pass::Point), std::sync::atomic::Ordering::Relaxed);
                ran.boxes.fetch_add(8 * g.len() * tape.work(Pass::Box), std::sync::atomic::Ordering::Relaxed);
                // A sample on 0 nothing near decides shows neither side.
                let (hl, hh) = held_of(b);
                let open = drawn(&tape, ctx, &at, &mut v, Some((std::array::from_fn(|a| I::new(hl[a], hh[a])), &**whole)));
                let shows = |i: usize| !(v[i] == 0.0 && open.contains(&at[i]));
                let (mut hits, mut next) = (Vec::new(), Vec::new());
                for (q, &(lo, hi, bi, oi)) in g.iter().enumerate() {
                    let side = b.open[oi as usize].1 & INSIDE != 0;
                    if let Some(i) = (0..27).find(|&i| inside(v[27 * q + i]) != side && shows(27 * q + i) && !covered(at[27 * q + i], at[27 * q + i])) {
                        hits.push((bi, oi, at[27 * q + i]));
                        continue;
                    }
                    let mid: V3 = std::array::from_fn(|a| 0.5 * (lo[a] + hi[a]));
                    for o in 0..8 {
                        let pick = |a: usize, x: [V3; 2]| x[(o >> a) & 1][a];
                        let (clo, chi): (V3, V3) = (std::array::from_fn(|a| pick(a, [lo, mid])), std::array::from_fn(|a| pick(a, [mid, hi])));
                        if !(proven(&tape, ctx, grid, clo, chi, side) || covered(clo, chi)) {
                            next.push((clo, chi, bi, oi));
                        }
                    }
                }
                (hits, next)
            })
            .collect();
        for (hits, _) in &level {
            for &(bi, oi, p) in hits {
                let o = &mut bricks[bi as usize].open[oi as usize];
                if o.1 & FOUND == 0 {
                    o.1 |= FOUND;
                    found.push((bi, oi, p));
                }
            }
        }
        parts = level.into_iter().flat_map(|l| l.1).filter(|p| bricks[p.2 as usize].open[p.3 as usize].1 & FOUND == 0).collect();
    }
    // Cells with parts left are undecided.
    stopped.extend(parts.iter().map(|p| (p.2, p.3, spacing)));
    for p in &stopped {
        bricks[p.0 as usize].open[p.1 as usize].1 &= !PROVEN;
    }
    // How thick what was found is: by brick, the field focused on it and
    // a cell round it.
    found.sort_by_key(|f| (f.0, f.1));
    let by_brick: Vec<&[(u32, u32, V3)]> = found.chunk_by(|a, b| a.0 == b.0).collect();
    let bricks = &*bricks;
    let out: Vec<([i64; 3], V3, bool, f64)> = by_brick
        .par_iter()
        .flat_map_iter(|g| {
            let b = &bricks[g[0].0 as usize];
            let (lo, hi) = (grid.pos(b.i0), grid.pos(b.i0).map(|x| x + B as f64 * h));
            let near = crate::query::focus(whole, ctx, lo.map(|x| x - h), hi.map(|x| x + h));
            let pts: Vec<(V3, bool)> = g.iter().map(|f| (f.2, b.open[f.1 as usize].1 & INSIDE == 0)).collect();
            let runs = across(&near, ctx, &pts, h);
            g.iter().zip(pts).zip(runs).map(|((f, (p, material)), r)| (cell_of(b.i0, b.open[f.1 as usize].0), p, material, r)).collect::<Vec<_>>()
        })
        .collect();
    let searched = if coarsest > 0.0 {
        coarsest
    } else if samples > 0.0 {
        spacing
    } else {
        h
    };
    let mut at: HashMap<(u32, u32), f64> = HashMap::new();
    for &(bi, oi, s) in &stopped {
        let e = at.entry((bi, oi)).or_insert(s);
        *e = e.max(s);
    }
    (out, searched, samples, at)
}

/// Steps a cell across, either way, on a line through the material.
const ACROSS: usize = 16;

/// How thick what holds points is (each point and whether it is in
/// material, else in a gap): at each, its run across along the field's
/// gradient there (square to its level surface: a run along an axis is cut
/// short where it leaves by a sheet's rim), to where the field crosses
/// zero either side (the first sample on the other side, a sixteenth of
/// `h` apart up to `h` away, then halving the step to it twenty times; h
/// where none is).
fn across(tape: &Tape, ctx: &Ctx, pts: &[(V3, bool)], h: f64) -> Vec<f64> {
    let step = h / ACROSS as f64;
    let at: Vec<V3> = pts.iter().map(|p| p.0).collect();
    // Each point's gradient, each way: a ray (point, unit, side).
    let mut rays: Vec<(V3, V3, bool)> = Vec::with_capacity(2 * pts.len());
    for (&(p, side), g) in pts.iter().zip(eval_grads(tape, ctx, &at)) {
        let n = norm([g[1], g[2], g[3]]);
        let u = if n > 0.0 && n.is_finite() { [g[1] / n, g[2] / n, g[3] / n] } else { [0.0, 0.0, 1.0] };
        rays.push((p, u, side));
        rays.push((p, u.map(|x| -x), side));
    }
    let along = |(p, u, _): (V3, V3, bool), s: f64| -> V3 { std::array::from_fn(|a| p[a] + s * u[a]) };
    let line: Vec<V3> = rays.iter().flat_map(|&r| (1..=ACROSS).map(move |i| along(r, step * i as f64))).collect();
    let mut v = Vec::new();
    eval_points(tape, ctx, &line, &mut v);
    // Each ray's last distance on its side and first off it (none: h).
    let mut span: Vec<(f64, f64)> = v
        .chunks(ACROSS)
        .zip(&rays)
        .map(|(f, r)| match f.iter().position(|&x| inside(x) != r.2) {
            Some(i) => (step * i as f64, step * (i + 1) as f64),
            None => (h, h),
        })
        .collect();
    for _ in 0..20 {
        let mids: Vec<V3> = rays.iter().zip(&span).map(|(&r, &(a, b))| along(r, 0.5 * (a + b))).collect();
        eval_points(tape, ctx, &mids, &mut v);
        for ((s, &f), r) in span.iter_mut().zip(&v).zip(&rays) {
            let m = 0.5 * (s.0 + s.1);
            if s.0 < s.1 {
                *s = if inside(f) == r.2 { (m, s.1) } else { (s.0, m) };
            }
        }
    }
    span.chunks(2).map(|w| w[0].1 + w[1].1).collect()
}
