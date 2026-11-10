//! Sparse manifold dual contouring.
//!
//! The octree finds bricks of `B^3` cells that straddle the surface. Each brick
//! evaluates its corner grid in SIMD batches. Within each sign-changing cell,
//! edge crossings are grouped into surface patches (connected through the
//! cell's faces; ambiguous faces are resolved by the field at the face centre,
//! which both neighbouring cells compute identically). Each patch gets one
//! vertex placed by a quadratic error function (QEF) over surface normals, which
//! keeps sharp edges and corners. One quad is emitted per sign-changing grid
//! edge, using the patch in each of its four cells that owns that edge. One
//! vertex per patch (not per cell) keeps thin walls and sheets manifold.
//!
//! What sampling at `h` cannot see is reported, not hidden: `mesh_checked`
//! takes every cell whose corners agree in sign and proves, by the field's
//! range over it (interval arithmetic), that it holds that sign throughout;
//! the cells it does not prove are searched, all at once level by level:
//! split in eight, the parts' corners sampled and their ranges taken, again
//! in the parts not proven, while each level costs within what the search
//! is allowed. A cell where a sample has the other sign holds material or a
//! gap the mesh leaves out (and how thin it is, is measured); one neither
//! found nor proven is undecided, down to the finest samples taken
//! (`Unresolved::searched`). Every such cell is reported. A mesh whose grid
//! would take more bricks than it is allowed is not made (`Limit`).

use crate::octree::{Cell, Grid, Local, descend_while, root_for};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_interval, eval_points};
use ink_core::interval::I;
use ink_core::sign::{drawn, inside};
use ink_core::tape::{Sheet, Tape};
use ink_core::v3::{EMPTY_BOX, box_union, cross, dot, norm, normal};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

const B: i64 = 8;
const N: usize = B as usize + 1;

#[derive(Clone, Debug, Default)]
pub struct Mesh {
    pub verts: Vec<V3>,
    pub tris: Vec<[u32; 3]>,
}

impl Mesh {
    /// The sum of its triangles' areas.
    pub fn area(&self) -> f64 {
        self.tris.iter().map(|t| 0.5 * norm(normal(t.map(|i| self.verts[i as usize])))).sum()
    }
}

struct Brick {
    i0: [i64; 3],
    /// Its corners' values, where a cell of it changes sign (else none).
    vals: Vec<f64>,
    /// (cell, mask of the 12 cell edges this patch owns, vertex position)
    patches: Vec<([i64; 3], u16, V3)>,
    /// Its cells not shown to hold nothing the mesh does not draw
    /// (`unshown`): each cell's index in the brick, and its state
    /// (`INSIDE`, `PROVEN`, `FOUND`, `CROSSED`, `SEAM`).
    open: Vec<(u16, u8)>,
    /// Its samples on 0 (corners, saddle-face centres, crossings' steps)
    /// nothing near decides (`ink_core::sign::drawn`): drawn as air,
    /// reported (a wall may be drawn where pieces meet flush).
    zeros: Vec<V3>,
    /// Beside each open cell the surface crosses, how thick what it may
    /// hide is (Tape::sheet: a seam's reach, else its box, 2 h); 0 for one
    /// it does not.
    thin: Vec<f64>,
}

/// An open cell's corners are inside; the search showed it holds nothing
/// the mesh does not (each part of it keeps one sign, or lies in the box of
/// a crossed cell shown drawn); a sample in it has the other sign; its
/// corners do not agree (the surface crosses it, drawn, but may hold more).
/// Else it is undecided.
const INSIDE: u8 = 1;
const PROVEN: u8 = 2;
const FOUND: u8 = 4;
const CROSSED: u8 = 8;
/// A crossed cell whose box holds a seam (Tape::switches).
const SEAM: u8 = 16;

/// What a mesh may sample (`mesh_checked`): its grid, in bricks of B^3
/// cells (past that it is not made), and its search for what the grid
/// leaves out, in samples (each level of it within what is left).
#[derive(Clone, Copy, Debug, Default)]
pub struct Limit {
    pub bricks: usize,
    pub search: f64,
}

/// What a brick costs, in the triangles a view mesh's budget counts: what
/// a flat surface across it draws, two triangles a cell of a face of it.
pub const BRICK_TRIANGLES: f64 = 2.0 * (B * B) as f64;

impl Limit {
    /// No limit on the grid, and no search.
    pub const GRID: Limit = Limit { bricks: usize::MAX, search: 0.0 };

    /// What `triangles` of a budget pay for: a grid of as many bricks
    /// (BRICK_TRIANGLES each), or a search of as many samples (a sample
    /// costing what one of a brick's does).
    pub fn bricks_for(triangles: f64) -> usize {
        (triangles / BRICK_TRIANGLES).max(1.0) as usize
    }
    pub fn samples_for(triangles: f64) -> f64 {
        triangles * (N * N * N) as f64 / BRICK_TRIANGLES
    }
}

/// What a mesh at `h` leaves out (see the module notes): cells whose
/// corners agree in sign that hold the other sign somewhere, found by a
/// sample, or may (undecided: the field's range over the parts left does
/// not rule it out, down to samples `searched` apart); of them, those where
/// material was found among corners outside, those undecided, and those
/// within a cell of one the surface crosses (`beside`), and the undecided
/// ones not (`away`); the thinnest the material found is across (0 where
/// there is none), a point in one (a sample of the other sign where one was
/// found), what sampling it cost (`cost`, in triangles: a brick
/// BRICK_TRIANGLES, a sample of the search as much as one of a brick's),
/// the res it was asked for where it is coarser (`refused`, set by
/// its maker: its grid there took more bricks than it was allowed; 0:
/// none), and the instructions it ran: over boxes, its cells' proofs and
/// its search's (`read`: each tape's length, a pass over it), and at
/// points, its grid's, crossings' and search's samples (`sampled`: a
/// gradient as four).
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct Unresolved {
    pub cells: usize,
    pub material: usize,
    pub thinnest: f64,
    pub undecided: usize,
    pub beside: usize,
    pub away: usize,
    pub at: Option<V3>,
    pub searched: f64,
    pub cost: f64,
    pub refused: f64,
    pub read: usize,
    pub sampled: usize,
    pub regions: Vec<Region>,
    pub seams: usize,
    pub thinner_than: f64,
}

/// Where a mesh leaves something out (`Unresolved::regions`): a box holding
/// every cell of it (a cover: much of the box may be drawn), how many cells
/// it holds, and of them those where material or a gap was found and those
/// undecided, and those beside a cell the surface crosses; the coarsest
/// spacing an undecided one was searched at, the thinnest material found
/// across, and a point in one (`found`: a sample of the other sign, else an
/// undecided cell's centre); `seams`: undecided cells the surface crosses
/// at a select's switch (Tape::switches). Bricks of B^3 cells that touch are one region;
/// `groups` counts the touching groups merged into it, where there were
/// more than MAX_REGIONS (none is left out).
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct Region {
    pub lo: V3,
    pub hi: V3,
    pub cells: usize,
    pub material: usize,
    pub gaps: usize,
    pub undecided: usize,
    pub beside: usize,
    pub seams: usize,
    pub searched: f64,
    pub thinner_than: f64,
    pub thinnest: f64,
    pub at: V3,
    pub found: bool,
    pub groups: usize,
}

/// The most regions a mesh reports: past it, touching groups are merged,
/// by the octree's blocks, into covers.
pub const MAX_REGIONS: usize = 16;

/// Unresolved cells by brick (its index in B^3 cells), each a region.
#[derive(Default)]
struct Regions(HashMap<[i64; 3], Region>);

impl Regions {
    /// Cell `c` of `grid`, counted and bounded in its brick's region, and
    /// `f` of it.
    fn add(&mut self, grid: &Grid, c: [i64; 3], f: impl FnOnce(&mut Region)) {
        let (lo, hi) = (grid.pos(c), grid.pos(c.map(|x| x + 1)));
        let r = self.0.entry(c.map(|x| x.div_euclid(B))).or_insert_with(|| Region { lo, hi, groups: 1, ..Default::default() });
        r.lo = std::array::from_fn(|a| r.lo[a].min(lo[a]));
        r.hi = std::array::from_fn(|a| r.hi[a].max(hi[a]));
        r.cells += 1;
        f(r);
    }

    /// The regions: bricks that touch joined, then, past `most`, groups
    /// merged by ever larger blocks of bricks until no more than `most`.
    fn merged(self, most: usize) -> Vec<Region> {
        let keys: Vec<[i64; 3]> = {
            let mut k: Vec<_> = self.0.keys().copied().collect();
            k.sort();
            k
        };
        let index: HashMap<[i64; 3], usize> = keys.iter().enumerate().map(|(i, k)| (*k, i)).collect();
        let mut parent: Vec<usize> = (0..keys.len()).collect();
        fn root(p: &mut [usize], mut i: usize) -> usize {
            while p[i] != i {
                p[i] = p[p[i]];
                i = p[i];
            }
            i
        }
        for (i, k) in keys.iter().enumerate() {
            for q in 0..27 {
                let n = [k[0] + q % 3 - 1, k[1] + q / 3 % 3 - 1, k[2] + q / 9 - 1];
                if let Some(&j) = index.get(&n) {
                    let (a, b) = (root(&mut parent, i), root(&mut parent, j));
                    parent[a.max(b)] = a.min(b);
                }
            }
        }
        let join = |a: &mut Region, b: &Region| {
            a.lo = std::array::from_fn(|i| a.lo[i].min(b.lo[i]));
            a.hi = std::array::from_fn(|i| a.hi[i].max(b.hi[i]));
            (a.cells, a.material, a.gaps, a.undecided, a.beside, a.seams) =
                (a.cells + b.cells, a.material + b.material, a.gaps + b.gaps, a.undecided + b.undecided, a.beside + b.beside, a.seams + b.seams);
            a.searched = a.searched.max(b.searched);
            a.thinner_than = a.thinner_than.max(b.thinner_than);
            a.thinnest = match (a.thinnest, b.thinnest) {
                (0.0, t) | (t, 0.0) => t,
                (x, y) => x.min(y),
            };
            if !a.found && b.found {
                (a.at, a.found) = (b.at, true);
            }
        };
        // Touching bricks, one region each group, keyed by its first brick.
        let mut groups: Vec<([i64; 3], Region)> = Vec::new();
        let mut at: HashMap<usize, usize> = HashMap::new();
        for (i, k) in keys.iter().enumerate() {
            let r = root(&mut parent, i);
            let reg = &self.0[k];
            match at.get(&r) {
                Some(&g) => join(&mut groups[g].1, reg),
                None => {
                    at.insert(r, groups.len());
                    groups.push((*k, reg.clone()));
                }
            }
        }
        // Past `most`, groups in one block of 2^s bricks merged, s growing.
        let mut s = 0;
        while groups.len() > most {
            s += 1;
            let mut by: Vec<([i64; 3], Region)> = Vec::new();
            let mut block: HashMap<[i64; 3], usize> = HashMap::new();
            for (k, r) in groups {
                let b = k.map(|x| x >> s);
                match block.get(&b) {
                    Some(&g) => {
                        join(&mut by[g].1, &r);
                        by[g].1.groups += r.groups;
                    }
                    None => {
                        block.insert(b, by.len());
                        by.push((k, r));
                    }
                }
            }
            groups = by;
        }
        groups.into_iter().map(|g| g.1).collect()
    }
}

#[inline]
fn vidx(i: usize, j: usize, k: usize) -> usize {
    ink_core::grid::index([i, j, k], [N; 3])
}

/// Cell corner `c` has offset bits (x = 1, y = 2, z = 4). The 12 cell edges
/// join corners differing in one bit.
const EDGES: [(usize, usize); 12] = {
    let mut e = [(0, 0); 12];
    let mut n = 0;
    let mut c = 0;
    while c < 8 {
        let mut bit = 1;
        while bit < 8 {
            if c & bit == 0 {
                e[n] = (c, c | bit);
                n += 1;
            }
            bit <<= 1;
        }
        c += 1;
    }
    e
};

fn edge_id(a: usize, b: usize) -> usize {
    let (a, b) = (a.min(b), a.max(b));
    EDGES.iter().position(|&e| e == (a, b)).unwrap()
}

/// The six faces as rings of four corners, and the four ring edges.
struct Face {
    axis: usize,
    ring: [usize; 4],
    edges: [usize; 4],
}

fn faces() -> &'static [Face; 6] {
    static F: std::sync::OnceLock<[Face; 6]> = std::sync::OnceLock::new();
    F.get_or_init(|| {
        std::array::from_fn(|f| {
            let (axis, side) = (f / 2, f % 2);
            let (a, u, v) = (1 << axis, 1 << ((axis + 1) % 3), 1 << ((axis + 2) % 3));
            let base = if side == 1 { a } else { 0 };
            let ring = [base, base | u, base | u | v, base | v];
            let edges = std::array::from_fn(|k| edge_id(ring[k], ring[(k + 1) % 4]));
            Face { axis, ring, edges }
        })
    })
}

/// Group a cell's crossing edges into patches (bit masks over the 12 edges).
fn patches(inside_c: &[bool; 8], face_center_inside: impl Fn(usize) -> bool) -> Vec<u16> {
    let mut parent: [usize; 12] = std::array::from_fn(|i| i);
    fn find(p: &mut [usize; 12], mut x: usize) -> usize {
        while p[x] != x {
            x = p[x];
        }
        x
    }
    let union = |p: &mut [usize; 12], a: usize, b: usize| {
        let (ra, rb) = (find(p, a), find(p, b));
        p[ra] = rb;
    };
    for (fi, f) in faces().iter().enumerate() {
        let r = f.ring.map(|c| inside_c[c]);
        let crossing: Vec<usize> = (0..4).filter(|&k| r[k] != r[(k + 1) % 4]).collect();
        match crossing.len() {
            2 => union(&mut parent, f.edges[crossing[0]], f.edges[crossing[1]]),
            4 => {
                // Saddle: if the centre matches ring corner 0, corners 0 and 2
                // connect through it and the segments cut off corners 1 and 3.
                if face_center_inside(fi) == r[0] {
                    union(&mut parent, f.edges[0], f.edges[1]);
                    union(&mut parent, f.edges[2], f.edges[3]);
                } else {
                    union(&mut parent, f.edges[3], f.edges[0]);
                    union(&mut parent, f.edges[1], f.edges[2]);
                }
            }
            _ => {}
        }
    }
    let mut out: Vec<(usize, u16)> = Vec::new();
    for (e, &(a, b)) in EDGES.iter().enumerate() {
        if inside_c[a] == inside_c[b] {
            continue;
        }
        let root = find(&mut parent, e);
        match out.iter_mut().find(|(r, _)| *r == root) {
            Some((_, m)) => *m |= 1 << e,
            None => out.push((root, 1 << e)),
        }
    }
    out.into_iter().map(|(_, m)| m).collect()
}

pub fn mesh(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64) -> Mesh {
    sample(tape, ctx, lo, hi, h, None, 0.0, false).unwrap().0
}

/// `mesh`, or None where the call's time limit passes first (`stop`).
pub fn mesh_in_time(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64) -> Option<Mesh> {
    sample(tape, ctx, lo, hi, h, None, 0.0, true).map(|m| m.0)
}

/// `mesh`, and what it leaves out (see the module notes), within `limit`:
/// None where its grid would take more bricks than that. A sample the
/// search finds is not counted left out where the mesh draws it: material
/// inside the solid it bounds, air outside, or either within `within` of
/// it (the accuracy it is made to: a curved surface bulges past its
/// triangles by up to that, as simplifying may move it).
pub fn mesh_checked(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64, limit: Limit, within: f64) -> Option<(Mesh, Unresolved)> {
    sample(tape, ctx, lo, hi, h, Some(limit), within, false)
}

/// `timed`: None where the call's time limit passes first.
#[allow(clippy::too_many_arguments)]
fn sample(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64, check: Option<Limit>, within: f64, timed: bool) -> Option<(Mesh, Unresolved)> {
    use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
    let pad = [lo[0] - h, lo[1] - h, lo[2] - h];
    let (grid, root) = root_for(pad, [hi[0] + h, hi[1] + h, hi[2] + h], h, B);
    // Bricks reached; past the limit the descent stops and none is made.
    let (reached, most, ran) = (AtomicUsize::new(0), check.map_or(usize::MAX, |l| l.bricks), Ran::default());
    let go = || reached.load(Relaxed) <= most && !(timed && crate::stop::reached());
    let leaf = |t: &Arc<Tape>, c| (reached.fetch_add(1, Relaxed) < most).then(|| brick(t, tape, ctx, &grid, c, check.is_some().then_some(&ran))).flatten();
    // Each brick's tape holds half a cell round it: a vertex may be placed
    // that far out of its cell (`qef`), and a cell the surface crosses is
    // shown drawn over that box (`brick`).
    let mut bricks: Vec<Brick> = descend_while(tape, ctx, &grid, root, B, 0.5, &go, &leaf);
    if !go() {
        return None;
    }

    // cell -> range of its patches in `owners`; owners hold (edge mask, vertex).
    let total: usize = bricks.iter().map(|b| b.patches.len()).sum();
    let mut cells: HashMap<[i64; 3], (u32, u8)> = HashMap::with_capacity(total);
    let mut owners: Vec<(u16, u32)> = Vec::with_capacity(total);
    let mut verts = Vec::with_capacity(total);
    for b in &bricks {
        for &(c, mask, p) in &b.patches {
            let e = cells.entry(c).or_insert((owners.len() as u32, 0));
            e.1 += 1;
            owners.push((mask, verts.len() as u32));
            verts.push(p);
        }
    }
    // Each step after the descent stops at the time limit too (`timed`).
    let stopped = || timed && crate::stop::reached();
    if stopped() {
        return None;
    }
    let quads: Vec<[u32; 4]> = bricks.par_iter().flat_map_iter(|b| if stopped() { Vec::new() } else { quads(b, &cells, &owners) }).collect();
    let tris = split_quads(&Local::new(tape, ctx, grid.origin, B as f64 * h), &verts, quads, &stopped)?;
    if stopped() {
        return None;
    }
    let unresolved = match check {
        Some(limit) => {
            // The search goes no finer than the finest length ink speaks of there.
            // What lies in the box of a cell the surface crosses that is
            // shown drawn (`unshown`) is drawn.
            let failing: std::collections::HashSet<[i64; 3]> = bricks.iter().flat_map(|b| b.open.iter().filter(|o| o.1 & CROSSED != 0).map(move |o| cell_of(b.i0, o.0))).collect();
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
            let beside = |c: [i64; 3]| (0..27).any(|q| q != 13 && cells.contains_key(&[c[0] + q % 3 - 1, c[1] + q / 3 % 3 - 1, c[2] + q / 9 - 1]));
            let cost = reached.load(Relaxed) as f64 * BRICK_TRIANGLES + samples * BRICK_TRIANGLES / (N * N * N) as f64;
            let mut u = Unresolved { searched, cost, read: ran.boxes.load(Relaxed), sampled: ran.points.load(Relaxed), ..Default::default() };
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
            // An example away from what is drawn, where there is one.
            u.at = far.or(u.at);
            u.regions = regions.merged(MAX_REGIONS);
            u
        }
        None => Default::default(),
    };
    Some((Mesh { verts, tris }, unresolved))
}

/// The instructions a mesh ran evaluating its field over boxes and at
/// points (Unresolved::read, sampled).
#[derive(Default)]
struct Ran {
    boxes: std::sync::atomic::AtomicUsize,
    points: std::sync::atomic::AtomicUsize,
}

/// A closed mesh's triangles by the squares of side `s` their extents
/// across y and z meet, to find what a ray along +x crosses.
struct Solid {
    s: f64,
    squares: HashMap<[i64; 2], Vec<u32>>,
    cubes: HashMap<[i64; 3], Vec<u32>>,
}

impl Solid {
    /// None unless the mesh is closed: each edge two triangles'.
    fn new(verts: &[V3], tris: &[[u32; 3]], s: f64) -> Option<Solid> {
        let mut edges: HashMap<(u32, u32), u32> = HashMap::with_capacity(3 * tris.len());
        for t in tris {
            for e in 0..3 {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                *edges.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        if tris.is_empty() || edges.values().any(|&n| n != 2) {
            return None;
        }
        let (mut squares, mut cubes): (HashMap<[i64; 2], Vec<u32>>, HashMap<[i64; 3], Vec<u32>>) = (HashMap::new(), HashMap::new());
        for (k, t) in tris.iter().enumerate() {
            let p = t.map(|i| verts[i as usize]);
            let span = |a: usize| {
                let (lo, hi) = (p.iter().map(|v| v[a]).fold(f64::INFINITY, f64::min), p.iter().map(|v| v[a]).fold(f64::NEG_INFINITY, f64::max));
                ((lo / s).floor() as i64)..=((hi / s).floor() as i64)
            };
            for y in span(1) {
                for z in span(2) {
                    squares.entry([y, z]).or_default().push(k as u32);
                    for x in span(0) {
                        cubes.entry([x, y, z]).or_default().push(k as u32);
                    }
                }
            }
        }
        Some(Solid { s, squares, cubes })
    }

    /// Whether `p` lies inside the solid (a ray from it along +x crosses
    /// its triangles an odd number of times); None where the ray meets an
    /// edge or vertex, or `p` lies on a triangle, as rounding can tell.
    fn inside(&self, verts: &[V3], tris: &[[u32; 3]], p: V3) -> Option<bool> {
        let key = [(p[1] / self.s).floor() as i64, (p[2] / self.s).floor() as i64];
        let mut odd = false;
        for &k in self.squares.get(&key).map_or(&[][..], Vec::as_slice) {
            let [a, b, c] = tris[k as usize].map(|i| verts[i as usize]);
            // Twice the signed areas, across y and z, of p with each edge.
            let side = |u: V3, v: V3| (u[1] - p[1]) * (v[2] - p[2]) - (u[2] - p[2]) * (v[1] - p[1]);
            let (d0, d1, d2) = (side(b, c), side(c, a), side(a, b));
            let scale = [a, b, c].iter().flat_map(|v| [(v[1] - p[1]).abs(), (v[2] - p[2]).abs()]).fold(0.0f64, f64::max);
            let eps = 1e-9 * scale * scale;
            let all = d0 + d1 + d2;
            if (d0 > eps && d1 > eps && d2 > eps) || (d0 < -eps && d1 < -eps && d2 < -eps) {
                let x = (d0 * a[0] + d1 * b[0] + d2 * c[0]) / all;
                if (x - p[0]).abs() <= 1e-9 * (x.abs() + p[0].abs() + self.s) {
                    return None;
                }
                odd ^= x > p[0];
            } else if !(d0 < -eps || d1 < -eps || d2 < -eps) || !(d0 > eps || d1 > eps || d2 > eps) {
                // On an edge or vertex across y and z (or the triangle is
                // edge-on there): only if p is not plainly outside it.
                return None;
            }
        }
        Some(odd)
    }

    /// Whether a triangle lies within `r` of `p` (r taken no more than `s`).
    fn near(&self, verts: &[V3], tris: &[[u32; 3]], p: V3, r: f64) -> bool {
        let r = r.min(self.s);
        if r <= 0.0 {
            return false;
        }
        let c = p.map(|x| (x / self.s).floor() as i64);
        (0..27).any(|q| {
            let key = [c[0] + q % 3 - 1, c[1] + q / 3 % 3 - 1, c[2] + q / 9 - 1];
            self.cubes.get(&key).is_some_and(|ks| ks.iter().any(|&k| norm(sub3(closest_on_triangle(p, tris[k as usize].map(|i| verts[i as usize])), p)) <= r))
        })
    }
}

fn sub3(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// The point of triangle `t` nearest `p` (Ericson, Real-Time Collision
/// Detection, 5.1.5).
fn closest_on_triangle(p: V3, t: [V3; 3]) -> V3 {
    let [a, b, c] = t;
    let (ab, ac, ap) = (sub3(b, a), sub3(c, a), sub3(p, a));
    let at = |u: f64, v: f64| -> V3 { std::array::from_fn(|k| a[k] + u * ab[k] + v * ac[k]) };
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub3(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return at(d1 / (d1 - d3), 0.0);
    }
    let cp = sub3(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return at(0.0, d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return std::array::from_fn(|k| b[k] + w * (c[k] - b[k]));
    }
    let denom = 1.0 / (va + vb + vc);
    at(vb * denom, vc * denom)
}

/// Cell `k` (its index in the brick) of the brick at `i0`, in the grid.
fn cell_of(i0: [i64; 3], k: u16) -> [i64; 3] {
    let (b, k) = (B as u16, k);
    [i0[0] + (k % b) as i64, i0[1] + (k / b % b) as i64, i0[2] + (k / (b * b)) as i64]
}

/// A brick's cells: `tape` simplified over it (and half a cell round it),
/// `whole` the field anywhere; with `check` (counting the instructions its
/// proofs run), those to search too, and its samples on 0 nothing near
/// decides.
fn brick(tape: &Tape, whole: &Tape, ctx: &Ctx, grid: &Grid, c: Cell, ran: Option<&Ran>) -> Option<Brick> {
    let check = ran.is_some();
    let count = |points: usize| _ = ran.map(|r| r.points.fetch_add(points * tape.len(), std::sync::atomic::Ordering::Relaxed));
    // What deciding samples on 0 runs (`drawn`), for each: its gradient and
    // two samples beside it, and at least its cube's bounds, solid's and
    // air's (more where they split it).
    let zeros_ran = |vals: &[f64]| {
        let z = vals.iter().filter(|v| **v == 0.0).count();
        count(6 * z);
        _ = ran.map(|r| r.boxes.fetch_add(4 * z * tape.len(), std::sync::atomic::Ordering::Relaxed));
    };
    let i0 = c.i;
    let mut pts = Vec::with_capacity(N * N * N);
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                pts.push(grid.pos([i0[0] + i as i64, i0[1] + j as i64, i0[2] + k as i64]));
            }
        }
    }
    let mut vals = Vec::new();
    eval_points(tape, ctx, &pts, &mut vals);
    count(pts.len());
    // Where `tape` holds: the brick and half a cell round it.
    let held: [I; 3] = std::array::from_fn(|a| I::new(pts[0][a] - 0.5 * grid.h, pts[pts.len() - 1][a] + 0.5 * grid.h));
    let within = Some((held, whole));
    // Samples on 0 as drawn (`drawn`): those nothing near decides are
    // reported.
    zeros_ran(&vals);
    let mut zeros = drawn(tape, ctx, &pts, &mut vals, within);
    if !check {
        zeros.clear();
    }

    // Pass 1: sign-changing cells, and the centres of their saddle faces.
    let mut mixed: Vec<([usize; 3], [bool; 8])> = Vec::new();
    let mut face_keys: HashMap<([i64; 3], usize), usize> = HashMap::new();
    let mut centers: Vec<V3> = Vec::new();
    for k in 0..B as usize {
        for j in 0..B as usize {
            for i in 0..B as usize {
                let corner = |c: usize| vals[vidx(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1))];
                let ins: [bool; 8] = std::array::from_fn(|c| inside(corner(c)));
                if ins.iter().all(|&v| v) || ins.iter().all(|&v| !v) {
                    continue;
                }
                for f in faces() {
                    let r = f.ring.map(|c| ins[c]);
                    if r[0] == r[2] && r[1] == r[3] && r[0] != r[1] {
                        let key = face_key([i0[0] + i as i64, i0[1] + j as i64, i0[2] + k as i64], f);
                        face_keys.entry(key).or_insert_with(|| {
                            centers.push(face_center(grid, key));
                            centers.len() - 1
                        });
                    }
                }
                mixed.push(([i, j, k], ins));
            }
        }
    }
    let (mut open, mut thin) = (Vec::new(), Vec::new());
    if let Some(ran) = ran {
        let mut crossed = vec![false; (B * B * B) as usize];
        for (c, _) in &mixed {
            crossed[c[0] + B as usize * (c[1] + B as usize * c[2])] = true;
        }
        let (mut cells, mut n) = (Vec::new(), 0);
        unshown(tape, ctx, grid, &pts, &vals, &crossed, [0; 3], B as usize, false, &mut cells, &mut n);
        cells.sort_unstable_by_key(|o| o.0);
        (open, thin) = cells.into_iter().map(|(k, state, t)| ((k, state), t)).unzip();
        ran.boxes.fetch_add(n, std::sync::atomic::Ordering::Relaxed);
    }
    if mixed.is_empty() {
        thin.resize(open.len(), 0.0);
        return (!open.is_empty() || !zeros.is_empty()).then(|| Brick { i0, vals: Vec::new(), patches: Vec::new(), open, zeros, thin });
    }
    let mut center_vals = Vec::new();
    eval_points(tape, ctx, &centers, &mut center_vals);
    count(centers.len());
    zeros_ran(&center_vals);
    let open_centres = drawn(tape, ctx, &centers, &mut center_vals, within);
    if check {
        zeros.extend(open_centres);
    }

    // Pass 2: patches and their edge crossings.
    let mut owners: Vec<([i64; 3], u16, usize, usize)> = Vec::new(); // cell, mask, crossing start, count
    let mut brackets: Vec<(V3, f64, V3, f64)> = Vec::new();
    for ([i, j, k], ins) in &mixed {
        let cell = [i0[0] + *i as i64, i0[1] + *j as i64, i0[2] + *k as i64];
        let corner = |c: usize| vidx(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1));
        let fc = |fi: usize| inside(center_vals[face_keys[&face_key(cell, &faces()[fi])]]);
        for mask in patches(ins, fc) {
            let start = brackets.len();
            for (e, &(a, b)) in EDGES.iter().enumerate() {
                if mask & (1 << e) == 0 {
                    continue;
                }
                let (a, b) = (corner(a), corner(b));
                let (pa, va, pb, vb) = (pts[a], vals[a], pts[b], vals[b]);
                brackets.push(if inside(va) { (pa, va, pb, vb) } else { (pb, vb, pa, va) });
            }
            owners.push((cell, mask, start, brackets.len() - start));
        }
    }
    let (cross, steps) = refine_crossings(tape, Some((held, whole)), ctx, &mut brackets);
    if check {
        zeros.extend(steps);
    }
    let grads = eval_grads(tape, ctx, &cross);
    count(brackets.len() * REFINE_STEPS + 4 * cross.len());
    let placed: Vec<V3> = owners.iter().map(|&(cell, _, s, n)| qef(&cross[s..s + n], &grads[s..s + n], grid.pos(cell), grid.h)).collect();
    let patches = owners.iter().zip(placed).map(|(&(cell, mask, _, _), v)| (cell, mask, v)).collect();
    thin.resize(open.len(), 0.0);
    Some(Brick { i0, vals, patches, open, zeros, thin })
}

/// The finest length ink speaks of in `lo..hi` (docs/CORE.md, Range): a
/// nanometre, a billionth of the farthest coordinate farther out.
fn finest_length(lo: V3, hi: V3) -> f64 {
    1e-9 * lo.iter().chain(&hi).fold(1e3f64, |m, x| m.max(x.abs()))
}

/// Whether the field's range over `lo..hi` (within the grid's region)
/// proves it drawn on side `inside_` throughout: solid, or air
/// (ink_core::sign: split where pieces meet flush; a region outside the
/// grid's is: nothing there is meshed).
fn proven(tape: &Tape, ctx: &Ctx, grid: &Grid, lo: V3, hi: V3, inside_: bool) -> bool {
    let b: [I; 3] = std::array::from_fn(|a| I::new(lo[a].max(grid.clip.0[a]), hi[a].min(grid.clip.1[a])));
    if b.iter().any(|i| !(i.lo <= i.hi)) {
        return true;
    }
    if inside_ { ink_core::sign::solid(tape, ctx, b) } else { ink_core::sign::air(tape, ctx, b) }
}

/// The cells of the block of n^3 at `c` of a brick (`pts` its corners,
/// `vals` the field there, `crossed` the cells whose corners do not agree)
/// not shown to hold nothing the mesh does not draw, added to `open` with
/// how thick what a crossed one may hide is (0 for the others):
///
/// - A cell whose corners agree, where the field's range does not prove it
///   keeps their sign throughout, nor split where pieces meet flush
///   (`proven`'s solid or air).
/// - A cell the surface crosses (CROSSED, and SEAM where a select may
///   switch there), where the surface is not one sheet over the box its
///   vertex may be placed in, the cell grown by half a cell (`sheet`). Not
///   where only each part of it is: that would not rule out a closed piece
///   across parts, each part one sheet along its own way (argued; no such
///   piece was found to pass).
///
/// Each is shown over a block at once (its box; for crossed cells their
/// boxes together: one sheet over it is one over each), else each eighth
/// is tried, down to a cell (`drawn`: the block's crossed cells are shown
/// already). Each block's tape is its parent's focused on its box, so a
/// cell's proof reads only what can reach it (a sweep's nearest segments,
/// not all of them). `read` counts the instructions run.
#[allow(clippy::too_many_arguments)]
fn unshown(
    tape: &Tape,
    ctx: &Ctx,
    grid: &Grid,
    pts: &[V3],
    vals: &[f64],
    crossed: &[bool],
    c: [usize; 3],
    n: usize,
    drawn: bool,
    open: &mut Vec<(u16, u8, f64)>,
    read: &mut usize,
) {
    let (b, h) = (B as usize, grid.h);
    let cell = (c[0] + b * (c[1] + b * c[2])) as u16;
    let (lo, hi) = (pts[vidx(c[0], c[1], c[2])], pts[vidx(c[0] + n, c[1] + n, c[2] + n)]);
    let any_crossed = (0..n * n * n).any(|q| crossed[(c[0] + q % n) + b * ((c[1] + q / n % n) + b * (c[2] + q / (n * n)))]);
    let (focused, drawn) = match (any_crossed, drawn) {
        (false, _) => {
            let side = inside(vals[vidx(c[0], c[1], c[2])]);
            // Within the grid's region: past it nothing is meshed.
            let bx: [I; 3] = std::array::from_fn(|a| I::new(lo[a].max(grid.clip.0[a]), hi[a].min(grid.clip.1[a])));
            if bx.iter().any(|i| !(i.lo <= i.hi)) {
                return;
            }
            let (r, focused) = focus(tape, ctx, bx, read);
            if if side { r.all_inside() } else { r.none_inside() } {
                return;
            }
            // Else split where pieces meet flush (ink_core::sign): solid
            // (a join inside the material) or air (a flush mouth).
            let t = focused.as_ref().unwrap_or(tape);
            *read += 4 * t.len();
            if if side { ink_core::sign::solid(t, ctx, bx) } else { ink_core::sign::air(t, ctx, bx) } {
                return;
            }
            if n == 1 {
                open.push((cell, if side { INSIDE } else { 0 }, 0.0));
                return;
            }
            (focused, true)
        }
        (true, true) => (None, true),
        (true, false) => {
            let (lo, hi) = (lo.map(|x| x - 0.5 * h), hi.map(|x| x + 0.5 * h));
            let focused = focus(tape, ctx, std::array::from_fn(|a| I::new(lo[a], hi[a])), read).1;
            let t = focused.as_ref().unwrap_or(tape);
            let (held, seam) = sheet(t, ctx, lo, hi, read);
            if n == 1 && held != Sheet::One {
                let thin = match held {
                    Sheet::Thin(t) => t.min(2.0 * h),
                    _ => 2.0 * h,
                };
                open.push((cell, CROSSED | if seam { SEAM } else { 0 }, thin));
            }
            (focused, held == Sheet::One)
        }
    };
    if n == 1 {
        return;
    }
    let (t, m) = (focused.as_ref().unwrap_or(tape), n / 2);
    for k in 0..8 {
        unshown(t, ctx, grid, pts, vals, crossed, [c[0] + (k & 1) * m, c[1] + (k >> 1 & 1) * m, c[2] + (k >> 2 & 1) * m], m, drawn, open, read);
    }
}

/// `tape`'s range over the box `b`, and `tape` with each min, max and
/// select whose side it decides there taken (None: it decides none): the
/// same field within `b` (`read` counts the instructions run).
fn focus(tape: &Tape, ctx: &Ctx, b: [I; 3], read: &mut usize) -> (I, Option<Tape>) {
    *read += tape.len();
    let mut ch = vec![ink_core::tape::BOTH; tape.len()];
    let r = eval_interval(tape, ctx, b, Some(&mut ch));
    (r, ch.iter().any(|&x| x != ink_core::tape::BOTH).then(|| tape.simplify(&ch)))
}

/// What a search costs per part, in samples: its 27 corners, and the
/// field's range over each of its 8 parts (about what 2 samples cost).
const PART_COST: f64 = 27.0 + 8.0 * 2.0;

/// What the field over `lo..hi`, round a cell the surface crosses (the box
/// its vertex may be placed in, the cell grown by half a cell), holds besides what the mesh draws
/// there (Tape::sheet). One sheet: it changes sign at most once along one
/// direction, so no piece lies wholly inside, and that sheet and the
/// cell's vertex both lie in the box, within its diagonal (2 sqrt 3 h) of
/// each other: drawn. Else how thick what it may hide is, where a seam's
/// reach bounds it. In a cell the surface does not cross no vertex is
/// placed: a sheet there (a spike between samples) is left out, so this is
/// no proof there. Checked along each axis (a face square to one is
/// constant along the others), then along the middle of the field's
/// gradient's bounds (between the faces at an edge) and along its
/// gradient at the box's middle (any direction would do): where a select
/// may switch (a seam), as its reach; elsewhere over the box, and else
/// over it split where pieces meet (`split`). The best found, and whether
/// there is a seam. `read` counts the instructions run.
fn sheet(tape: &Tape, ctx: &Ctx, lo: V3, hi: V3, read: &mut usize) -> (Sheet, bool) {
    *read += 4 * tape.len();
    let b: [I; 3] = std::array::from_fn(|a| I::new(lo[a], hi[a]));
    let mut ch = vec![ink_core::tape::BOTH; tape.len()];
    eval_interval(tape, ctx, b, Some(&mut ch));
    let g = ink_core::eval::eval_interval_grad_every(tape, ctx, b);
    let mid: V3 = std::array::from_fn(|a| {
        let i = g[g.len() - 1].d[a];
        0.5 * (i.lo + i.hi)
    });
    let centre: V3 = std::array::from_fn(|a| 0.5 * (lo[a] + hi[a]));
    let at = eval_grads(tape, ctx, &[centre])[0];
    let unit = |v: V3| {
        let n = norm(v);
        (n > 0.0 && n.is_finite()).then(|| v.map(|x| x / n))
    };
    let seam = tape.switches(&ch);
    let dirs = [unit(mid), unit([at[1], at[2], at[3]])];
    let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    *read += 3 * tape.len();
    if axes.iter().any(|&d| tape.sheet(&ch, &g, None, d) == Sheet::One) {
        return (Sheet::One, seam);
    }
    if seam {
        let mut best = Sheet::Unknown;
        for d in dirs.into_iter().flatten() {
            *read += tape.len();
            // A seam's reach over the box. (Over its eighths too, one
            // direction for all, it proved a third of a duct's seam cells
            // drawn, 3% of a sample's cells, for a fifth more time.)
            let held = seam_reach(tape, ctx, lo, hi, 0, d, &mut 2, read);
            if held == Sheet::One {
                return (Sheet::One, seam);
            }
            best = least(best, held);
        }
        return (best, seam);
    }
    let dirs: Vec<V3> = dirs.into_iter().flatten().chain(axes).collect();
    *read += 2 * tape.len();
    if dirs[..dirs.len() - 3].iter().any(|&d| tape.sheet(&ch, &g, None, d) == Sheet::One) || dirs.iter().any(|&d| split(tape, ctx, b, d, ink_core::sign::SPLITS, read) != 0) {
        return (Sheet::One, seam);
    }
    (Sheet::Unknown, seam)
}

/// Of two sound accounts of a region (Tape::sheet), the stronger.
fn least(a: Sheet, b: Sheet) -> Sheet {
    match (a, b) {
        (Sheet::One, _) | (_, Sheet::One) => Sheet::One,
        (Sheet::Thin(x), Sheet::Thin(y)) => Sheet::Thin(x.min(y)),
        (Sheet::Unknown, x) | (x, Sheet::Unknown) => x,
    }
}

/// Tape::sheet along `d` over `lo..hi`, and where that is unknown over its
/// eighths, `depth` times: all its parts changing sign the same way along d
/// (`way`; else unknown: parts changing sign opposite ways would not rule
/// out a closed piece across them). One sheet where every part has
/// one, else the thickest reach of a part's seam (unknown where a part's
/// is, at the last depth). `read` counts the instructions run.
#[allow(clippy::too_many_arguments)]
fn seam_reach(tape: &Tape, ctx: &Ctx, lo: V3, hi: V3, depth: u32, d: V3, way: &mut i8, read: &mut usize) -> Sheet {
    *read += 4 * tape.len();
    let b: [I; 3] = std::array::from_fn(|a| I::new(lo[a], hi[a]));
    let mut ch = vec![ink_core::tape::BOTH; tape.len()];
    eval_interval(tape, ctx, b, Some(&mut ch));
    let g = ink_core::eval::eval_interval_grad_every(tape, ctx, b);
    let centre: V3 = std::array::from_fn(|a| 0.5 * (lo[a] + hi[a]));
    let point = ink_core::eval::eval_interval_grad_every(tape, ctx, centre.map(I::point));
    let half: V3 = std::array::from_fn(|a| 0.5 * (hi[a] - lo[a]));
    let (held, w) = tape.sheet_way(&ch, &g, Some((&point, half)), d);
    // Its own account, where its way agrees with the others'.
    let mut agreed = *way;
    let own = match (held, (agreed, w)) {
        (Sheet::Unknown, _) => Sheet::Unknown,
        (h, (2, w) | (w, 2)) => {
            agreed = w;
            h
        }
        (h, (a, b)) if a == b => h,
        _ => Sheet::Unknown,
    };
    if own == Sheet::One || depth == 0 {
        *way = agreed;
        return own;
    }
    // Else its eighths', where those are all known (with their ways).
    let mut parts_way = *way;
    let mut parts = Sheet::One;
    for o in 0..8 {
        let plo: V3 = std::array::from_fn(|a| if (o >> a) & 1 == 0 { lo[a] } else { centre[a] });
        let phi: V3 = std::array::from_fn(|a| if (o >> a) & 1 == 0 { centre[a] } else { hi[a] });
        parts = match (parts, seam_reach(tape, ctx, plo, phi, depth - 1, d, &mut parts_way, read)) {
            (_, Sheet::Unknown) => Sheet::Unknown,
            (Sheet::Unknown, _) => Sheet::Unknown,
            (Sheet::One, x) | (x, Sheet::One) => x,
            (Sheet::Thin(x), Sheet::Thin(y)) => Sheet::Thin(x.max(y)),
        };
        if parts == Sheet::Unknown {
            break;
        }
    }
    match (own, parts) {
        (_, Sheet::One) => {
            *way = parts_way;
            Sheet::One
        }
        (Sheet::Thin(x), Sheet::Thin(y)) if y < x => {
            *way = parts_way;
            Sheet::Thin(y)
        }
        (Sheet::Unknown, p @ Sheet::Thin(_)) => {
            *way = parts_way;
            p
        }
        (o, _) => {
            *way = agreed;
            o
        }
    }
}

/// The field's way of changing sign along unit d over box b where it is
/// one sheet there (Tape::sheet_way: 1, -1, 2; else 0, unknown), or where
/// that is unknown,
/// over b split by a plane where pieces meet (ink_core::sign::levels, a
/// coordinate's): one half air (ink_core::sign::air), the other changing
/// sign once along d, where each line along d leaves it into the air
/// after going from inside to not, or enters it from the air before going
/// from not to inside (the air adds no change of sign); or, d along the
/// plane (no line crosses it), the halves changing sign the same way (or
/// one keeping one sign); up to `depth` splits.
fn way(tape: &Tape, ctx: &Ctx, b: [I; 3], d: V3, depth: usize, read: &mut usize) -> i8 {
    *read += 4 * tape.len();
    let mut ch = vec![ink_core::tape::BOTH; tape.len()];
    eval_interval(tape, ctx, b, Some(&mut ch));
    let g = ink_core::eval::eval_interval_grad_every(tape, ctx, b);
    let w = match tape.sheet_way(&ch, &g, None, d) {
        (Sheet::One, w) => w,
        _ => 0,
    };
    if w != 0 || depth == 0 {
        return w;
    }
    split(tape, ctx, b, d, depth, read)
}

/// `way` over b split by a plane where pieces meet, as there (0 where no
/// such plane crosses b).
fn split(tape: &Tape, ctx: &Ctx, b: [I; 3], d: V3, depth: usize, read: &mut usize) -> i8 {
    let plane = ink_core::sign::levels(tape).iter().find_map(|&(i, c)| ink_core::sign::axis(tape, i as usize).filter(|&k| b[k].lo < c && c < b[k].hi).map(|k| (k, c)));
    let Some((k, c)) = plane else { return 0 };
    let (mut lower, mut upper) = (b, b);
    lower[k] = I::new(b[k].lo, c);
    upper[k] = I::new(c, b[k].hi);
    // d along the plane: each line along d stays in one half, so the halves
    // need only agree (one changing sign throughout one way, the other too
    // or keeping one sign: the plane between them is then a face of the
    // sheet, as at an edge).
    if d[k] == 0.0 {
        let (wl, wu) = (way(tape, ctx, lower, d, depth - 1, read), way(tape, ctx, upper, d, depth - 1, read));
        match (wl, wu) {
            (0, _) | (_, 0) => {}
            (2, w) | (w, 2) => return w,
            (a, b) if a == b => return a,
            _ => {}
        }
    }
    // s: the side of the plane the air is on.
    for (air, other, s) in [(lower, upper, -1.0), (upper, lower, 1.0)] {
        *read += 2 * tape.len();
        if !ink_core::sign::air(tape, ctx, air) {
            continue;
        }
        let toward = d[k] * s;
        return match way(tape, ctx, other, d, depth - 1, read) {
            w @ (1 | -1) if toward == 0.0 || f64::from(w) * toward > 0.0 => w,
            // Inside or air throughout: its surface, if any, is the plane.
            2 if toward == 0.0 => 2,
            2 => {
                if toward > 0.0 {
                    1
                } else {
                    -1
                }
            }
            _ => 0,
        };
    }
    0
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
fn search(
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
        crate::query::focus(whole, ctx, lo, hi)
    };
    // The parts left (box, brick, cell's place in its open cells), in brick
    // order, as many as what is left affords; the spacing of the samples
    // the least searched cell left undecided had (h: none but its corners).
    let (mut spacing, mut left, mut samples, mut coarsest) = (h, allowed, 0.0, 0.0f64);
    let affords = |left: f64| (left / PART_COST).max(0.0) as usize;
    let mut parts: Vec<(V3, V3, u32, u32)> = Vec::new();
    let first = if 0.5 * h >= finest { affords(left).min(open) } else { 0 };
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
                ran.points.fetch_add(at.len() * tape.len(), std::sync::atomic::Ordering::Relaxed);
                ran.boxes.fetch_add(8 * g.len() * tape.len(), std::sync::atomic::Ordering::Relaxed);
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

/// Canonical key of a cell face: (global min corner of the face, normal axis).
fn face_key(cell: [i64; 3], f: &Face) -> ([i64; 3], usize) {
    let mut g = cell;
    g[f.axis] += ((f.ring[0] >> f.axis) & 1) as i64;
    (g, f.axis)
}

/// Face centre, computed only from the canonical key so every cell sharing the
/// face gets bit-identical coordinates (and therefore the same field value).
fn face_center(grid: &Grid, (g, axis): ([i64; 3], usize)) -> V3 {
    std::array::from_fn(|a| {
        let i = g[a] as f64 + if a == axis { 0.0 } else { 0.5 };
        grid.origin[a] + grid.h * i
    })
}

const REFINE_STEPS: usize = 6;

/// Move each edge crossing onto the surface: batched false-position steps,
/// safeguarded to shrink the bracket every iteration. Linear interpolation
/// alone lands off the surface near sharp features, where distance fields curve.
/// A step sample on 0 is read as every sample is (`drawn`, `within` where
/// given), and those undecided are returned with the crossings, for the
/// mesh's report.
fn refine_crossings(tape: &Tape, within: Option<([I; 3], &Tape)>, ctx: &Ctx, br: &mut [(V3, f64, V3, f64)]) -> (Vec<V3>, Vec<V3>) {
    let mut open = Vec::new();
    let lerp = |(pa, va, pb, vb): (V3, f64, V3, f64)| {
        let t = if va == vb { 0.5 } else { (va / (va - vb)).clamp(0.1, 0.9) };
        [pa[0] + (pb[0] - pa[0]) * t, pa[1] + (pb[1] - pa[1]) * t, pa[2] + (pb[2] - pa[2]) * t]
    };
    let mut mid: Vec<V3> = Vec::with_capacity(br.len());
    let mut v = Vec::new();
    for _ in 0..REFINE_STEPS {
        mid.clear();
        mid.extend(br.iter().map(|b| lerp(*b)));
        eval_points(tape, ctx, &mid, &mut v);
        open.extend(drawn(tape, ctx, &mid, &mut v, within));
        for ((b, p), &f) in br.iter_mut().zip(&mid).zip(&v) {
            if inside(f) == inside(b.1) {
                b.0 = *p;
                b.1 = f;
            } else {
                b.2 = *p;
                b.3 = f;
            }
        }
    }
    let cross = br
        .iter()
        .map(|&(pa, va, pb, vb)| {
            let t = if va == vb { 0.5 } else { (va / (va - vb)).clamp(0.0, 1.0) };
            [pa[0] + (pb[0] - pa[0]) * t, pa[1] + (pb[1] - pa[1]) * t, pa[2] + (pb[2] - pa[2]) * t]
        })
        .collect();
    (cross, open)
}

/// Minimize sum (n . (v - p))^2 + lambda |v - mass|^2, clamped near the cell.
fn qef(pts: &[V3], grads: &[[f64; 4]], cmin: V3, h: f64) -> V3 {
    let n = pts.len() as f64;
    let mut mass = [0.0; 3];
    for p in pts {
        for a in 0..3 {
            mass[a] += p[a] / n;
        }
    }
    let lambda = 0.05;
    let mut m = [[0.0; 3]; 3];
    let mut rhs = [0.0; 3];
    for a in 0..3 {
        m[a][a] = lambda;
    }
    for (p, g) in pts.iter().zip(grads) {
        let len = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
        if !(len > 0.0 && len.is_finite()) {
            continue;
        }
        let nn = [g[1] / len, g[2] / len, g[3] / len];
        let d = nn[0] * (p[0] - mass[0]) + nn[1] * (p[1] - mass[1]) + nn[2] * (p[2] - mass[2]);
        for a in 0..3 {
            for b in 0..3 {
                m[a][b] += nn[a] * nn[b];
            }
            rhs[a] += nn[a] * d;
        }
    }
    let x = solve3(m, rhs);
    // Features can sit just outside the cell whose edges detected them; allow
    // half a cell of slack so sharp edges and corners are reproduced exactly.
    let m = 0.5 * h;
    std::array::from_fn(|a| (mass[a] + x[a]).clamp(cmin[a] - m, cmin[a] + h + m))
}

/// m x = b by Cramer's rule. m is lambda I plus outer products of unit
/// normals: positive definite, its determinant at least lambda^3.
fn solve3(m: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let det = |m: &[[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(&m);
    std::array::from_fn(|c| {
        let mut mc = m;
        for r in 0..3 {
            mc[r][c] = b[r];
        }
        det(&mc) / d
    })
}

/// Split each quad along the diagonal whose midpoint lies closer to the true
/// surface. At sharp creases the wrong diagonal cuts a notch into the edge.
/// The midpoints are evaluated with the tape simplified over each block
/// they fall in (`Local`): the same values, without the whole field's cost
/// at each. In runs, each begun only until `stopped` (None then).
fn split_quads(local: &Local, verts: &[V3], quads: Vec<[u32; 4]>, stopped: &dyn Fn() -> bool) -> Option<Vec<[u32; 3]>> {
    let mid = |a: u32, b: u32| {
        let (p, q) = (verts[a as usize], verts[b as usize]);
        [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0, (p[2] + q[2]) / 2.0]
    };
    let mut out = Vec::with_capacity(quads.len() * 2);
    for run in quads.chunks(1 << 15) {
        if stopped() {
            return None;
        }
        let pts: Vec<V3> = run.iter().flat_map(|q| [mid(q[0], q[2]), mid(q[1], q[3])]).collect();
        let f = local.eval(&pts);
        for (k, q) in run.iter().enumerate() {
            if f[2 * k].abs() <= f[2 * k + 1].abs() {
                out.push([q[0], q[1], q[2]]);
                out.push([q[0], q[2], q[3]]);
            } else {
                out.push([q[0], q[1], q[3]]);
                out.push([q[1], q[2], q[3]]);
            }
        }
    }
    Some(out)
}

fn quads(b: &Brick, cells: &HashMap<[i64; 3], (u32, u8)>, owners: &[(u16, u32)]) -> Vec<[u32; 4]> {
    let mut out = Vec::new();
    if b.patches.is_empty() {
        return out;
    }
    // Vertex of the patch in `cell` that owns the grid edge from global corner
    // `g` along axis `d`.
    let vertex = |cell: [i64; 3], g: [i64; 3], d: usize| -> Option<u32> {
        let &(start, count) = cells.get(&cell)?;
        let c0 = ((g[0] - cell[0]) | ((g[1] - cell[1]) << 1) | ((g[2] - cell[2]) << 2)) as usize;
        let e = edge_id(c0, c0 | (1 << d));
        owners[start as usize..start as usize + count as usize].iter().find(|(m, _)| m & (1 << e) != 0).map(|&(_, v)| v)
    };
    for k in 0..B as usize {
        for j in 0..B as usize {
            for i in 0..B as usize {
                let v0 = b.vals[vidx(i, j, k)];
                let g = [b.i0[0] + i as i64, b.i0[1] + j as i64, b.i0[2] + k as i64];
                for d in 0..3 {
                    let v1 = match d {
                        0 => b.vals[vidx(i + 1, j, k)],
                        1 => b.vals[vidx(i, j + 1, k)],
                        _ => b.vals[vidx(i, j, k + 1)],
                    };
                    if inside(v0) == inside(v1) {
                        continue;
                    }
                    // (d, u, v) is a right-handed cyclic triple.
                    let (u, v) = ((d + 1) % 3, (d + 2) % 3);
                    let off = |du: i64, dv: i64| {
                        let mut c = g;
                        c[u] -= du;
                        c[v] -= dv;
                        vertex(c, g, d)
                    };
                    let (Some(q00), Some(q10), Some(q11), Some(q01)) = (off(1, 1), off(0, 1), off(0, 0), off(1, 0)) else {
                        continue;
                    };
                    // Counter-clockwise around +d when the surface normal points along +d.
                    out.push(if inside(v0) { [q00, q10, q11, q01] } else { [q00, q01, q11, q10] });
                }
            }
        }
    }
    out
}

// ---------------------------------------------------------------- analysis

#[derive(Clone, Debug)]
pub struct Shell {
    pub triangles: usize,
    pub volume: f64,
    pub area: f64,
    pub lo: V3,
    pub hi: V3,
}

#[derive(Clone, Debug)]
pub struct Analysis {
    pub area: f64,
    pub volume: f64,
    pub boundary_edges: usize,
    pub nonmanifold_edges: usize,
    pub shells: Vec<Shell>,
}

fn tri_terms(m: &Mesh, t: &[u32; 3]) -> (f64, f64) {
    let (a, b, c) = (m.verts[t[0] as usize], m.verts[t[1] as usize], m.verts[t[2] as usize]);
    let area = 0.5 * norm(normal([a, b, c]));
    let vol = dot(a, cross(b, c)) / 6.0;
    (area, vol)
}

/// The root of `x` in a union-find forest of parent links, halving the
/// path to it on the way.
pub fn find(p: &mut [u32], mut x: u32) -> u32 {
    while p[x as usize] != x {
        p[x as usize] = p[p[x as usize] as usize];
        x = p[x as usize];
    }
    x
}

pub fn analyze(m: &Mesh) -> Analysis {
    let mut keys: Vec<u64> = m
        .tris
        .par_iter()
        .flat_map_iter(|t| {
            (0..3).map(move |e| {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                ((a.min(b) as u64) << 32) | a.max(b) as u64
            })
        })
        .collect();
    keys.par_sort_unstable();
    let (mut boundary_edges, mut nonmanifold_edges) = (0, 0);
    for run in keys.chunk_by(|a, b| a == b) {
        match run.len() {
            1 => boundary_edges += 1,
            2 => {}
            _ => nonmanifold_edges += 1,
        }
    }

    // Union-find over vertices to split into connected shells.
    let mut parent: Vec<u32> = (0..m.verts.len() as u32).collect();
    for t in &m.tris {
        let r0 = find(&mut parent, t[0]);
        for &v in &t[1..] {
            let r = find(&mut parent, v);
            parent[r as usize] = r0;
        }
    }
    let mut shells: HashMap<u32, Shell> = HashMap::new();
    let (mut area, mut volume) = (0.0, 0.0);
    for t in &m.tris {
        let (a, v) = tri_terms(m, t);
        area += a;
        volume += v;
        let r = find(&mut parent, t[0]);
        let (lo, hi) = EMPTY_BOX;
        let s = shells.entry(r).or_insert(Shell { triangles: 0, volume: 0.0, area: 0.0, lo, hi });
        s.triangles += 1;
        s.area += a;
        s.volume += v;
        for &i in t {
            let p = m.verts[i as usize];
            (s.lo, s.hi) = box_union((s.lo, s.hi), (p, p));
        }
    }
    let mut shells: Vec<Shell> = shells.into_values().collect();
    shells.sort_by(|a, b| b.volume.abs().total_cmp(&a.volume.abs()));
    Analysis { area, volume, boundary_edges, nonmanifold_edges, shells }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink_core::graph::{Graph, UnOp, Var};

    /// The brick at the origin's cell 0 (`[0, 1]^3` at h = 1), and whether
    /// it is shown drawn (not left open).
    fn cell_0(g: &Graph, f: ink_core::graph::NodeId) -> (Brick, bool) {
        let (tape, ctx) = (Tape::new(g, f), Ctx::new(0.0, vec![]));
        let b = brick(&tape, &tape, &ctx, &Grid::new([0.0; 3], 1.0), Cell { i: [0; 3], size: B }, Some(&Ran::default())).unwrap();
        let drawn = !b.open.iter().any(|o| o.0 == 0);
        (b, drawn)
    }

    #[test]
    fn a_cell_is_shown_drawn_over_all_its_vertex_may_reach() {
        // z + max(|x - .5|, |y - .5|) - 1.2, a cone: one sheet round cell 0,
        // whose vertex is placed above it (z 1.19): it lies in the box
        // shown, the cell grown by half a cell.
        let mut g = Graph::new();
        let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let (half, top) = (g.constant(0.5), g.constant(1.2));
        let (dx, dy) = (g.sub(x, half), g.sub(y, half));
        let (ax, ay) = (g.un(UnOp::Abs, dx), g.un(UnOp::Abs, dy));
        let m = g.max(ax, ay);
        let s = g.add(z, m);
        let f = g.sub(s, top);
        let (b, drawn) = cell_0(&g, f);
        let p = b.patches.iter().find(|q| q.0 == [0; 3]).unwrap().2;
        assert!(drawn && p.iter().all(|&v| (-0.5..=1.5).contains(&v)), "{p:?}");
        // (z - .5) (1 - 8 max(x - 1.05, 0)): a plane in cell 0, turning over
        // past x = 1.175, within half a cell of it. One sheet in the cell,
        // not round it: not shown drawn.
        let mut g = Graph::new();
        let [x, z] = [Var::X, Var::Z].map(|v| g.var(v));
        let (half, at, zero, eight, one) = (g.constant(0.5), g.constant(1.05), g.constant(0.0), g.constant(8.0), g.constant(1.0));
        let (dz, dx) = (g.sub(z, half), g.sub(x, at));
        let past = g.max(dx, zero);
        let k = g.mul(eight, past);
        let k = g.sub(one, k);
        let f = g.mul(dz, k);
        assert!(!cell_0(&g, f).1);
    }

    #[test]
    fn a_crossed_cells_box_holding_a_piece_besides_its_sheet_is_not_shown_drawn() {
        // A plane, z = 0.25, crossing cell 0, and in its box a double cone
        // from x 0.3 to 0.7 at y 0.5, z 1.25, widest (0.05) at x 0.5: a
        // piece of its own beside the sheet. Not shown drawn. (A guard: a
        // rule accepting the box part by part is argued unsound here, a
        // closed piece across parts each one sheet along its own way, but
        // this cone is not accepted by such a rule either.)
        let mut g = Graph::new();
        let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let c = |g: &mut Graph, v: f64| g.constant(v);
        let (q, h) = (c(&mut g, 0.25), c(&mut g, 0.5));
        let plane = g.sub(z, q);
        let (dy, zc) = (g.sub(y, h), c(&mut g, 1.25));
        let dz = g.sub(z, zc);
        let (yy, zz) = (g.un(UnOp::Square, dy), g.un(UnOp::Square, dz));
        let s = g.add(yy, zz);
        let radial = g.un(UnOp::Sqrt, s);
        let dx = g.sub(x, h);
        let ax = g.un(UnOp::Abs, dx);
        let (k, w) = (c(&mut g, 0.25), c(&mut g, 0.05));
        let taper = g.mul(k, ax);
        let r = g.sub(w, taper);
        let cone = g.sub(radial, r);
        let f = g.min(plane, cone);
        assert!(!cell_0(&g, f).1);
    }

    #[test]
    fn a_comparison_that_jumps_past_equality_is_a_seam() {
        // select(a - y, x + 0.2 a^2 + z^2 - 0.05, x + 0.2 y^2 + z^2 - 0.05),
        // a = 2 floor(x) + 1: the sides agree where a = y, but a jumps from
        // -1 to 1 at x = 0, never meeting y, and the field jumps there,
        // a closed piece 0 <= x < 0.05 beside the sheet near x = -0.15.
        // Cell [-0.25, 0.25]^3 at h 0.5 is not shown drawn.
        let mut g = Graph::new();
        let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let (offset, two, one, k) = (g.constant(0.05), g.constant(2.0), g.constant(1.0), g.constant(0.2));
        let fl = g.un(UnOp::Floor, x);
        let twice = g.mul(two, fl);
        let a = g.add(twice, one);
        let cond = g.sub(a, y);
        let (aa, bb, zz) = (g.un(UnOp::Square, a), g.un(UnOp::Square, y), g.un(UnOp::Square, z));
        let (aa, bb) = (g.mul(k, aa), g.mul(k, bb));
        let (xa, xb) = (g.add(x, aa), g.add(x, bb));
        let (xa, xb) = (g.add(xa, zz), g.add(xb, zz));
        let (left, right) = (g.sub(xa, offset), g.sub(xb, offset));
        let f = g.select(cond, left, right);
        let (tape, ctx) = (Tape::new(&g, f), Ctx::new(0.0, vec![]));
        let b = brick(&tape, &tape, &ctx, &Grid::new([-0.25; 3], 0.5), Cell { i: [0; 3], size: B }, Some(&Ran::default())).unwrap();
        let mut v = Vec::new();
        eval_points(&tape, &ctx, &[[-0.2, 0.0, 0.0], [-0.1, 0.0, 0.0], [0.025, 0.0, 0.0], [0.1, 0.0, 0.0]], &mut v);
        assert_eq!(v.iter().map(|v| *v < 0.0).collect::<Vec<_>>(), [true, false, true, false]);
        assert!(b.open.iter().any(|o| o.0 == 0));
    }

    #[test]
    fn a_seams_bound_is_kept_whole() {
        // At h = 1e6 a seam's reach is about 800000.03; a ball 800000.019
        // across fits in the slab it may hide. Its bound is no less.
        let h = 1e6;
        let mut g = Graph::new();
        let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let c = |g: &mut Graph, v: f64| g.constant(v);
        let (lower, switch, upper, half, radius) = (c(&mut g, 0.1 * h), c(&mut g, 0.10000001 * h), c(&mut g, 0.90000003 * h), c(&mut g, 0.5 * h), c(&mut g, 0.45 * h));
        let (b, top) = (g.sub(y, lower), g.sub(y, upper));
        let (dx, dz) = (g.sub(x, half), g.sub(z, half));
        let (ax, az) = (g.un(UnOp::Abs, dx), g.un(UnOp::Abs, dz));
        let radial = g.max(ax, az);
        let radial = g.sub(radial, radius);
        let above = g.sub(switch, y);
        let cond = g.max(radial, above);
        let f = g.select(cond, top, b);
        let (tape, ctx) = (Tape::new(&g, f), Ctx::new(0.0, vec![]));
        let br = brick(&tape, &tape, &ctx, &Grid::new([0.0; 3], h), Cell { i: [0; 3], size: B }, Some(&Ran::default())).unwrap();
        let k = br.open.iter().position(|o| o.0 == 0).unwrap();
        let (_, u) = mesh_checked(&Arc::new(tape), &ctx, [0.0; 3], [h; 3], h, Limit::GRID, 0.0).unwrap();
        assert!(br.thin[k] >= 0.800000019 * h && u.thinner_than >= 0.800000019 * h, "{} {}", br.thin[k], u.thinner_than);
    }

    #[test]
    fn a_walls_face_is_one_sheet_in_a_box_reaching_its_middle() {
        // |z| - w: cell 0's box, z from -0.5 to 1.5, reaches the wall's
        // middle (z = 0). With w = 0.6 it holds one face only: drawn. With
        // w = 0.3 it holds both: not shown drawn.
        let wall = |w: f64| {
            let mut g = Graph::new();
            let z = g.var(Var::Z);
            let a = g.un(UnOp::Abs, z);
            let w = g.constant(w);
            let f = g.sub(a, w);
            cell_0(&g, f).1
        };
        assert!(wall(0.6));
        assert!(!wall(0.3));
    }
}
