use super::analysis::Mesh;
use super::brick::{B, Brick, N, brick, cell_of};
use super::contour::{quads, split_quads};
use super::flat::{Blocks, Even, FLAT_MOST, Flat, disc, flat, regions, rim, rim_cell};
use super::proof::{CROSSED, FOUND, PROVEN, SEAM, try_untried};
use super::report::{BRICK_TRIANGLES, Limit, MAX_REGIONS, Ran, Regions, Unresolved};
use super::search::{finest_length, search, searched_cells};
use super::solid::{Solid, untangle};
use crate::octree::{Cell, Grid, Local, root_for};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_points};
use ink_core::interval::I;
use ink_core::sign::inside;
use ink_core::tape::{Pass, Tape};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

pub fn mesh(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64) -> Mesh {
    sample(tape, ctx, lo, hi, h, None, 0.0, None, false).unwrap().0
}

/// `mesh`, or None where the call's time limit passes first (`stop`).
pub fn mesh_in_time(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64) -> Option<Mesh> {
    sample(tape, ctx, lo, hi, h, None, 0.0, None, true).map(|m| m.0)
}

/// `mesh`, and what it leaves out (see the module notes), within `limit`:
/// None where its grid would take more bricks than that. A sample the
/// search finds is not counted left out where the mesh draws it: material
/// inside the solid it bounds, air outside, or either within `within` of
/// it (the accuracy it is made to: a curved surface bulges past its
/// triangles by up to that, as simplifying may move it).
pub fn mesh_checked(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64, limit: Limit, within: f64) -> Option<(Mesh, Unresolved)> {
    mesh_whole(tape, ctx, lo, hi, h, limit, within, Some(&|_| true))
}

/// `mesh_checked`, drawing whole each block (`Flat`) where the surface is
/// shown flat within half of `within` (more than 0) and `even` holds over it;
/// none: sampled uniformly.
#[allow(clippy::too_many_arguments)]
pub fn mesh_whole(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64, limit: Limit, within: f64, even: Option<&Even<'_>>) -> Option<(Mesh, Unresolved)> {
    sample(tape, ctx, lo, hi, h, Some(limit), within, even.filter(|_| within > 0.0), false)
}

/// What a descent reaches: a brick, or a block drawn whole (its slabs).
enum Found {
    Brick(Brick),
    Flat(Vec<Flat>),
}

/// Where a descent goes: everywhere, a block drawn whole met going no
/// further; of a block drawn whole (its cell and axis), only its rim
/// (`rim`) and the cells between two heights along its axis (its plain
/// slabs'); or everywhere, drawing none whole.
#[derive(Clone, Copy)]
enum Part {
    All,
    Rim(Cell, usize, (i64, i64)),
    Plain,
}

/// The interval-pruned descent to bricks (as octree::descend_padded, each
/// cell's tape simplified over it grown by half a cell), trying each block
/// a brick to FLAT_MOST cells across as a block drawn whole (`flat`, only
/// where `part` is All: its slabs): one found is gone no further into. Each
/// brick is made by `leaf`, given the block drawn whole it lies in (`part`
/// Rim).
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn descend(
    tape: &Arc<Tape>,
    ctx: &Ctx,
    grid: &Grid,
    cell: Cell,
    part: Part,
    go: &(dyn Fn() -> bool + Sync),
    leaf: &(dyn Fn(&Arc<Tape>, Cell, Option<(Cell, usize)>) -> Option<Brick> + Sync),
    flat: &(dyn Fn(&Arc<Tape>, Cell) -> Option<Vec<Flat>> + Sync),
) -> Vec<Found> {
    if !grid.meets(cell) || !go() {
        return Vec::new();
    }
    if let Part::Rim(q, a, (lo, hi)) = part
        && !rim(cell, q, a)
        && !(cell.i[a] < hi && lo < cell.i[a] + cell.size)
    {
        return Vec::new();
    }
    let pad = 0.5 * grid.h;
    let b = grid.cell_box(cell).map(|i| I::new(i.lo - pad, i.hi + pad));
    let mut ch = vec![ink_core::tape::BOTH; tape.len()];
    let r = eval_interval(tape, ctx, b, Some(&mut ch));
    if r.clear_of_surface() || r.all_inside() {
        return Vec::new();
    }
    let sub = if ch.iter().any(|&c| c != ink_core::tape::BOTH) { Arc::new(tape.simplify(&ch)) } else { tape.clone() };
    if let Part::All = part
        && cell.size <= FLAT_MOST
        && let Some(f) = flat(&sub, cell)
    {
        return vec![Found::Flat(f)];
    }
    if cell.size <= B {
        let within = match part {
            Part::Rim(q, a, _) => Some((q, a)),
            _ => None,
        };
        return leaf(&sub, cell, within).map(Found::Brick).into_iter().collect();
    }
    let parts: Vec<Vec<Found>> = cell.children().par_iter().map(|k| descend(&sub, ctx, grid, *k, part, go, leaf, flat)).collect();
    parts.into_iter().flatten().collect()
}

/// `timed`: None where the call's time limit passes first. `even`: blocks
/// where it holds are drawn whole where they can be (`Flat`), within half
/// of `within`.
#[allow(clippy::too_many_arguments)]
fn sample(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64, check: Option<Limit>, within: f64, even: Option<&Even<'_>>, timed: bool) -> Option<(Mesh, Unresolved)> {
    use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
    let pad = [lo[0] - h, lo[1] - h, lo[2] - h];
    let (grid, root) = root_for(pad, [hi[0] + h, hi[1] + h, hi[2] + h], h, B);
    // The grid's points sampled, a brick's (N^3) for each it is allowed;
    // past the limit the descent stops and none is made.
    let (reached, most, ran) = (AtomicUsize::new(0), check.map_or(usize::MAX, |l| l.bricks.saturating_mul(N * N * N)), Ran::default());
    let go = || reached.load(Relaxed) <= most && !(timed && crate::stop::reached());
    let leaf = |t: &Arc<Tape>, c: Cell, _: Option<(Cell, usize)>| brick(t, tape, ctx, &grid, c, None, (&reached, most), check.is_some().then_some(&ran));
    let whole = |t: &Arc<Tape>, c| even.and_then(|e| flat(t, ctx, &grid, c, within, e, &ran));
    let none = |_: &Arc<Tape>, _: Cell| None;
    // Each brick's tape holds half a cell round it: a vertex may be placed
    // that far out of its cell (`qef`), and a cell the surface crosses is
    // shown drawn over that box (`brick`).
    let (mut bricks, mut flats): (Vec<Brick>, Vec<Flat>) = (Vec::new(), Vec::new());
    for f in descend(tape, ctx, &grid, root, Part::All, &go, &leaf, &whole) {
        match f {
            Found::Brick(b) => bricks.push(b),
            Found::Flat(f) => flats.extend(f),
        }
    }
    if !go() {
        return None;
    }
    // The slabs drawn whole joined in regions, and their blocks' rims' cells
    // meshed, a block at once.
    let mut blocks = Blocks::new(&flats);
    let slopes = regions(&flats, &mut blocks, h, within);
    {
        let (blocks, flats) = (&blocks, &flats);
        let rim_leaf = |t: &Arc<Tape>, c: Cell, q: Option<(Cell, usize)>| {
            let (q, a) = q?;
            let slabs = &blocks.at[&(q.i, q.size)];
            let plain = |n: [i64; 3]| slabs.iter().any(|&k| flats[k as usize].plain && (flats[k as usize].lo..flats[k as usize].hi).contains(&n[a]));
            let meshed = |n: [i64; 3]| plain(n) || rim_cell(flats, blocks, q, a, n);
            brick(t, tape, ctx, &grid, c, Some(&meshed), (&reached, most), check.is_some().then_some(&ran))
        };
        let mut firsts: Vec<u32> = blocks.at.values().map(|k| k[0]).collect();
        firsts.sort_unstable();
        let rims: Vec<Vec<Found>> = firsts
            .par_iter()
            .map(|&k| {
                let f = &flats[k as usize];
                let slabs = blocks.block(f.cell.i).iter().map(|&j| &flats[j as usize]).filter(|g| g.plain);
                let heights = slabs.fold((i64::MAX, i64::MIN), |(lo, hi), g| (lo.min(g.lo), hi.max(g.hi)));
                descend(&f.tape, ctx, &grid, f.cell, Part::Rim(f.cell, f.axis, heights), &go, &rim_leaf, &none)
            })
            .collect();
        bricks.extend(rims.into_iter().flatten().filter_map(|f| if let Found::Brick(b) = f { Some(b) } else { None }));
    }
    if !go() {
        return None;
    }
    // Each step after the descent stops at the time limit too (`timed`).
    let stopped = || timed && crate::stop::reached();
    // Made again, the regions whose rims make no loops to draw them whole
    // sampled whole, until each left does.
    let (mut cells, mut owners, mut verts, mut cell_of_vert, mut safe, mut quads_of, discs, stray, whole_regions);
    loop {
        // cell -> range of its patches in `owners`; owners hold (edge mask, vertex).
        let total: usize = bricks.iter().map(|b| b.patches.len()).sum();
        cells = HashMap::with_capacity(total);
        owners = Vec::with_capacity(total);
        verts = Vec::with_capacity(total);
        cell_of_vert = Vec::with_capacity(total);
        safe = Vec::with_capacity(total);
        for b in &bricks {
            for &(c, mask, p, q) in &b.patches {
                safe.push(q);
                let e = cells.entry(c).or_insert((owners.len() as u32, 0u8));
                e.1 += 1;
                owners.push((mask, verts.len() as u32));
                verts.push(p);
                cell_of_vert.push(c);
            }
        }
        if stopped() {
            return None;
        }
        let (blocks_, flats_) = (&blocks, &flats);
        let region_of = |c: [i64; 3]| if flats_.is_empty() { None } else { blocks_.region(flats_, c) };
        let made: Vec<(Vec<[u32; 4]>, Vec<(u32, [u32; 2])>)> =
            bricks.par_iter().map(|b| if stopped() { Default::default() } else { quads(b, &cells, &owners, &region_of) }).collect();
        let mut rims: Vec<Vec<[u32; 2]>> = vec![Vec::new(); slopes.len()];
        quads_of = Vec::with_capacity(made.iter().map(|m| m.0.len()).sum());
        for (q, r) in made {
            quads_of.extend(q);
            for (k, e) in r {
                rims[k as usize].push(e);
            }
        }
        let live: Vec<bool> = (0..slopes.len()).map(|r| blocks.region.contains(&Some(r as u32))).collect();
        let drawn: Vec<Option<(Vec<[u32; 3]>, f64)>> = (0..slopes.len())
            .into_par_iter()
            .map(|r| if live[r] { disc(&flats, &blocks, slopes[r], ctx, &verts, &cell_of_vert, &rims[r], within, &ran) } else { Some((Vec::new(), 0.0)) })
            .collect();
        if drawn.iter().all(Option::is_some) {
            stray = drawn.iter().flatten().fold(0.0f64, |m, d| m.max(d.1));
            whole_regions = live.iter().filter(|l| **l).count();
            discs = drawn.into_iter().flatten().flat_map(|d| d.0).collect::<Vec<_>>();
            break;
        }
        // The blocks of the regions failed, sampled whole, and so every
        // region with a slab in one (its rim there was not meshed), until
        // none is left: their rims' bricks made again with the rest.
        let mut failed: std::collections::BTreeSet<u32> = (0..slopes.len() as u32).filter(|&r| drawn[r as usize].is_none()).collect();
        let mut gone: Vec<([i64; 3], i64)> = Vec::new();
        loop {
            let more: Vec<u32> = (0..flats.len())
                .filter(|&k| blocks.region[k].is_some_and(|r| failed.contains(&r)))
                .flat_map(|k| blocks.at[&(flats[k].cell.i, flats[k].cell.size)].iter().filter_map(|&j| blocks.region[j as usize]).collect::<Vec<_>>())
                .filter(|r| !failed.contains(r))
                .collect();
            if more.is_empty() {
                break;
            }
            failed.extend(more);
        }
        for k in 0..flats.len() {
            if blocks.region[k].is_some_and(|r| failed.contains(&r)) {
                blocks.region[k] = None;
                gone.push((flats[k].cell.i, flats[k].cell.size));
            }
        }
        gone.sort_unstable();
        gone.dedup();
        let held = |b: &Brick, (i, size): ([i64; 3], i64)| (0..3).all(|k| b.i0[k] >= i[k] && b.i0[k] < i[k] + size);
        bricks.retain(|b| !gone.iter().any(|&g| held(b, g)));
        for &(i, size) in &gone {
            let f = &flats[blocks.at[&(i, size)][0] as usize];
            bricks.extend(descend(&f.tape, ctx, &grid, f.cell, Part::Plain, &go, &leaf, &none).into_iter().filter_map(|f| if let Found::Brick(b) = f { Some(b) } else { None }));
        }
        if !go() {
            return None;
        }
    }
    let local = Local::new(tape, ctx, grid.origin, B as f64 * h);
    let mut tris = split_quads(&local, &verts, &cell_of_vert, quads_of, h, &stopped)?;
    ran.points.fetch_add(local.work() + local.blocks_work(), Relaxed);
    tris.extend(discs);
    untangle(&mut verts, &safe, &tris, h);
    if stopped() {
        return None;
    }
    let unresolved = match check {
        Some(limit) => {
            // The search goes no finer than the finest length ink speaks of there.
            // Crossed cells left untried beside one the search reads are
            // tried: what their boxes hold may be drawn there.
            try_untried(&mut bricks, ctx, &grid, &ran, searched_cells(limit.search, h, finest_length(lo, hi)));
            // What lies in the box of a cell the surface crosses that is
            // shown drawn (`unshown`) is drawn; not one left untried.
            let failing: std::collections::HashSet<[i64; 3]> =
                bricks.iter().flat_map(|b| b.open.iter().filter(|o| o.1 & CROSSED != 0).map(|o| o.0).chain(b.untried.iter().map(|u| u.0)).map(move |k| cell_of(b.i0, k))).collect();
            let accepted = |c: [i64; 3]| cells.contains_key(&c) && !failing.contains(&c);
            let covered = |plo: V3, phi: V3| {
                // Cells whose box, the cell grown by h / 2, holds plo..phi.
                let range = |a: usize| {
                    let first = ((phi[a] - 1.5 * h - grid.origin[a]) / h).ceil() as i64;
                    let last = ((plo[a] + 0.5 * h - grid.origin[a]) / h).floor() as i64;
                    first..=last
                };
                range(0).any(|i| range(1).any(|j| range(2).any(|k| accepted([i, j, k]))))
            };
            let (found, searched, samples, spacing) = search(&mut bricks, tape, ctx, &grid, limit.search, finest_length(lo, hi), &covered, &ran);
            let mut regions = Regions::default();
            let around = |c: [i64; 3]| (0..27).filter(|&q| q != 13).map(move |q| [c[0] + q % 3 - 1, c[1] + q / 3 % 3 - 1, c[2] + q / 9 - 1]);
            // Of the cells in regions drawn whole (meshed only round their
            // rims) beside a cell reported, those the surface crosses: their
            // corners sampled, as their block's field reads them.
            let reported = found.iter().map(|f| f.0).chain(bricks.iter().flat_map(|b| {
                let open = b.open.iter().filter(|o| o.1 & (PROVEN | FOUND) == 0).map(|o| o.0);
                open.chain(b.untried.iter().map(|u| u.0)).map(move |k| cell_of(b.i0, k))
            }));
            let zero_cells = bricks.iter().flat_map(|b| &b.zeros).map(|p| std::array::from_fn(|a| ((p[a] - grid.origin[a]) / h).floor() as i64));
            let mut whole_beside: Vec<[i64; 3]> = match flats.is_empty() {
                true => Vec::new(),
                false => reported.chain(zero_cells).flat_map(around).filter(|n| !cells.contains_key(n) && blocks.region(&flats, *n).is_some()).collect(),
            };
            whole_beside.sort_unstable();
            whole_beside.dedup();
            let crossed_whole: std::collections::HashSet<[i64; 3]> = whole_beside
                .par_iter()
                .filter(|&&n| {
                    let f = &flats[blocks.of(&flats, n).unwrap() as usize];
                    let at: Vec<V3> = (0..8).map(|k| grid.pos([n[0] + (k & 1), n[1] + (k >> 1 & 1), n[2] + (k >> 2 & 1)])).collect();
                    let mut v = Vec::new();
                    eval_points(&f.tape, ctx, &at, &mut v);
                    ran.points.fetch_add(at.len() * f.tape.work(Pass::Point), Relaxed);
                    v.iter().any(|x| inside(*x) != inside(v[0]))
                })
                .copied()
                .collect();
            let beside = |c: [i64; 3]| around(c).any(|n| cells.contains_key(&n) || crossed_whole.contains(&n));
            let cost = (reached.load(Relaxed) as f64 + samples) * BRICK_TRIANGLES / (N * N * N) as f64;
            let mut u = Unresolved {
                searched,
                cost,
                read: ran.boxes.load(Relaxed),
                sampled: ran.points.load(Relaxed),
                tapes: ran.tapes.load(Relaxed),
                whole: whole_regions,
                whole_within: stray,
                ..Default::default()
            };
            // Material or a gap found: reported, and material made again
            // for under the thinnest of it. Not where the mesh draws it
            // (`mesh_checked`: a sheet bulging past a cell's corners, its
            // bulge drawn by the cells round it); that cell is undecided.
            let solid = if found.is_empty() { None } else { Solid::new(&verts, &tris, 4.0 * h) };
            let mut far = None;
            for &(c, p, material, across) in &found {
                let near = beside(c);
                u.cells += 1;
                u.beside += near as usize;
                u.at = u.at.or(Some(p));
                if solid.as_ref().is_some_and(|m| m.inside(&verts, &tris, p) == Some(material) || m.near(&verts, &tris, p, within)) {
                    u.undecided += 1;
                    u.away += !near as usize;
                    u.thinner_than = u.thinner_than.max(2.0 * h);
                    regions.add(&grid, c, |r| {
                        r.undecided += 1;
                        r.beside += near as usize;
                        r.thinner_than = r.thinner_than.max(2.0 * h);
                        r.searched = r.searched.max(h);
                        if r.cells == 1 {
                            r.at = p;
                        }
                    });
                    continue;
                }
                if material {
                    u.material += 1;
                    u.thinnest = if u.thinnest == 0.0 { across } else { u.thinnest.min(across) };
                }
                far = far.or((!near).then_some(p));
                regions.add(&grid, c, |r| {
                    (r.material, r.gaps) = (r.material + material as usize, r.gaps + !material as usize);
                    r.beside += near as usize;
                    if material {
                        r.thinnest = if r.thinnest == 0.0 { across } else { r.thinnest.min(across) };
                    }
                    if !r.found {
                        (r.at, r.found) = (p, true);
                    }
                });
            }
            for (bi, b) in bricks.iter().enumerate() {
                for (oi, &(k, state)) in b.open.iter().enumerate().filter(|o| o.1.1 & (PROVEN | FOUND) == 0) {
                    let c = cell_of(b.i0, k);
                    let near = beside(c);
                    u.cells += 1;
                    u.undecided += 1;
                    u.beside += near as usize;
                    u.away += !near as usize;
                    u.seams += (state & SEAM != 0) as usize;
                    let centre = grid.pos(c).map(|x| x + 0.5 * h);
                    u.at = u.at.or(Some(centre));
                    far = far.or((!near).then_some(centre));
                    let s = spacing.get(&(bi as u32, oi as u32)).copied().unwrap_or(h);
                    // What may hide there is thinner than: across a crossed
                    // cell's box, its bound; else what samples s apart
                    // miss (a ball sqrt 3 s across holds one).
                    let thinner = match state & CROSSED != 0 {
                        true => b.thin[oi],
                        false => 3f64.sqrt() * s,
                    };
                    u.thinner_than = u.thinner_than.max(thinner);
                    regions.add(&grid, c, |r| {
                        r.undecided += 1;
                        r.thinner_than = r.thinner_than.max(thinner);
                        r.seams += (state & SEAM != 0) as usize;
                        r.beside += near as usize;
                        r.searched = r.searched.max(s);
                        if r.cells == 1 {
                            r.at = centre;
                        }
                    });
                }
            }
            // Samples on 0 nothing near decides: each its cell, undecided.
            // Drawn as air, what they leave out lies in the cells round
            // them: thinner than 2 h, as across a crossed cell's box.
            for p in bricks.iter().flat_map(|b| &b.zeros) {
                let c: [i64; 3] = std::array::from_fn(|a| ((p[a] - grid.origin[a]) / h).floor() as i64);
                let near = beside(c);
                u.cells += 1;
                u.undecided += 1;
                u.thinner_than = u.thinner_than.max(2.0 * h);
                u.beside += near as usize;
                u.away += !near as usize;
                u.at = u.at.or(Some(*p));
                far = far.or((!near).then_some(*p));
                regions.add(&grid, c, |r| {
                    r.undecided += 1;
                    r.thinner_than = r.thinner_than.max(2.0 * h);
                    r.beside += near as usize;
                    if r.cells == 1 {
                        r.at = *p;
                    }
                });
            }
            // Cells left untried: in regions already reporting what may hide
            // as thick as across a crossed cell's box (2 h), as they would.
            for b in &bricks {
                for &(k, _) in &b.untried {
                    let c = cell_of(b.i0, k);
                    let near = beside(c);
                    (u.cells, u.untried, u.beside, u.away) = (u.cells + 1, u.untried + 1, u.beside + near as usize, u.away + !near as usize);
                    u.thinner_than = u.thinner_than.max(2.0 * h);
                    regions.add(&grid, c, |r| {
                        (r.untried, r.beside) = (r.untried + 1, r.beside + near as usize);
                        r.thinner_than = r.thinner_than.max(2.0 * h);
                    });
                }
            }
            // An example away from what is drawn, where there is one.
            u.at = far.or(u.at);
            u.regions = regions.merged(MAX_REGIONS);
            u
        }
        None => Default::default(),
    };
    Some((Mesh { verts, tris }, unresolved))
}
