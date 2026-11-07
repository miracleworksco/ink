//! Parallel interval-pruned octree descent.
//!
//! Cells proven empty are dropped, cells proven full are reported whole, and
//! the rest subdivide with a tape simplified for that region. Work is spread
//! across all cores with rayon. `Local` evaluates scattered points with the
//! same simplified tapes.

use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_interval, eval_points};
use ink_core::interval::I;
use ink_core::tape::{BOTH, Tape};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::{Arc, Mutex, RwLock};

#[derive(Clone, Copy, Debug)]
pub struct Grid {
    pub origin: V3,
    pub h: f64,
    /// The region its cells are clipped to (`cell_bounds`): a proof's grid
    /// is clipped to the region it was asked about (try_root_for), so every
    /// box it examines, and every witness it finds there, lies in it, and a
    /// tape simplified for that region (query::focus) is the field over all
    /// of each box. An approximate query's cells are whole cubes (root_for,
    /// Grid::new: clipped to all of space).
    pub clip: (V3, V3),
}

impl Grid {
    /// A grid of whole cubes.
    pub fn new(origin: V3, h: f64) -> Grid {
        Grid { origin, h, clip: ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]) }
    }

    pub fn pos(&self, i: [i64; 3]) -> V3 {
        [self.origin[0] + self.h * i[0] as f64, self.origin[1] + self.h * i[1] as f64, self.origin[2] + self.h * i[2] as f64]
    }

    /// A cell's corners, clipped to the grid's region (empty, lo > hi on
    /// some axis, for a cell outside it: `meets`).
    pub fn cell_bounds(&self, c: Cell) -> (V3, V3) {
        let (lo, hi) = (self.pos(c.i), self.pos([c.i[0] + c.size, c.i[1] + c.size, c.i[2] + c.size]));
        (std::array::from_fn(|k| lo[k].max(self.clip.0[k])), std::array::from_fn(|k| hi[k].min(self.clip.1[k])))
    }

    /// A cell's sides: the cube's, h times its size, on each axis it is
    /// not clipped on, else what clipping leaves (as cell_bounds).
    pub fn cell_sides(&self, c: Cell) -> V3 {
        let (lo, hi) = (self.pos(c.i), self.pos([c.i[0] + c.size, c.i[1] + c.size, c.i[2] + c.size]));
        std::array::from_fn(|k| if lo[k] >= self.clip.0[k] && hi[k] <= self.clip.1[k] { self.h * c.size as f64 } else { hi[k].min(self.clip.1[k]) - lo[k].max(self.clip.0[k]) })
    }

    /// Whether a cell holds a point of the grid's region.
    pub fn meets(&self, c: Cell) -> bool {
        let (lo, hi) = self.cell_bounds(c);
        (0..3).all(|k| lo[k] <= hi[k])
    }

    pub fn cell_box(&self, c: Cell) -> [I; 3] {
        let (lo, hi) = self.cell_bounds(c);
        [I::new(lo[0], hi[0]), I::new(lo[1], hi[1]), I::new(lo[2], hi[2])]
    }
}

/// A cube of `size` grid cells with minimum corner at grid index `i`.
#[derive(Clone, Copy, Debug)]
pub struct Cell {
    pub i: [i64; 3],
    pub size: i64,
}

impl Cell {
    pub fn children(&self) -> [Cell; 8] {
        let h = self.size / 2;
        std::array::from_fn(|k| Cell { i: [self.i[0] + h * (k as i64 & 1), self.i[1] + h * ((k as i64 >> 1) & 1), self.i[2] + h * ((k as i64 >> 2) & 1)], size: h })
    }
}

/// Grid anchored at `lo` with spacing `h`, and a root cell covering `lo..hi`
/// whose size is `leaf * 2^k`, for approximate queries (meshes, samples):
/// where no root covers the region (try_root_for), the largest one.
pub fn root_for(lo: V3, hi: V3, h: f64, leaf: i64) -> (Grid, Cell) {
    let (grid, size, _) = grow_root(lo, hi, h, leaf);
    (grid, Cell { i: [0; 3], size })
}

/// The root for a proof, its cells clipped to `lo..hi` (Grid::clip): None
/// where none covers `lo..hi` with finite
/// corners within 2^62 cells (an extent of 2^62 h or more, one that
/// overflows, an infinite one), and the proof answers unknown. The size
/// estimate rounds, so the root grows until its far corner, as the cells
/// compute it, reaches `hi`: at lo = -1, hi = 2^-54, h = 1 an estimate
/// alone gives the root [-1, 0], missing what lies above 0.
pub fn try_root_for(lo: V3, hi: V3, h: f64, leaf: i64) -> Option<(Grid, Cell)> {
    let (grid, size, covers) = grow_root(lo, hi, h, leaf);
    covers.then_some((Grid { clip: (lo, hi), ..grid }, Cell { i: [0; 3], size }))
}

fn grow_root(lo: V3, hi: V3, h: f64, leaf: i64) -> (Grid, i64, bool) {
    let ext = crate::query::max_extent(lo, hi);
    let n = ((ext / h).ceil() as i64).max(1);
    let grid = Grid::new(lo, h);
    let covers = |s: i64| (0..3).all(|a| lo[a].is_finite() && grid.pos([s; 3])[a].is_finite() && grid.pos([s; 3])[a] >= hi[a]);
    let mut size = leaf;
    while size < n || !covers(size) {
        match size.checked_mul(2) {
            Some(s) if s <= 1 << 62 => size = s,
            _ => return (grid, size, false),
        }
    }
    (grid, size, true)
}

pub fn descend<R, FL, FF>(tape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, cell: Cell, leaf: i64, on_leaf: &FL, on_full: &FF) -> Vec<R>
where
    R: Send,
    FL: Fn(&Arc<Tape>, Cell) -> Option<R> + Sync,
    FF: Fn(&Arc<Tape>, Cell) -> Option<R> + Sync,
{
    descend_inner(tape, ctx, grid, cell, leaf, 0.0, on_leaf, on_full)
}

/// `descend` where each cell's tape is simplified (and pruned) over the cell
/// grown by `pad` grid steps on every side, so the tape handed to `on_leaf`
/// stays valid a little outside the cell (for Newton steps that wander).
pub fn descend_padded<R, FL, FF>(tape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, cell: Cell, leaf: i64, pad: f64, on_leaf: &FL, on_full: &FF) -> Vec<R>
where
    R: Send,
    FL: Fn(&Arc<Tape>, Cell) -> Option<R> + Sync,
    FF: Fn(&Arc<Tape>, Cell) -> Option<R> + Sync,
{
    descend_inner(tape, ctx, grid, cell, leaf, pad, on_leaf, on_full)
}

/// `descend_padded` with no full cells kept, which stops going down once
/// `go` says so (a limit on the leaves reached).
pub fn descend_while<R, FL>(tape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, cell: Cell, leaf: i64, pad: f64, go: &(dyn Fn() -> bool + Sync), on_leaf: &FL) -> Vec<R>
where
    R: Send,
    FL: Fn(&Arc<Tape>, Cell) -> Option<R> + Sync,
{
    descend_go(tape, ctx, grid, cell, leaf, pad, go, on_leaf, &|_: &Arc<Tape>, _| None)
}

fn descend_inner<R, FL, FF>(tape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, cell: Cell, leaf: i64, pad: f64, on_leaf: &FL, on_full: &FF) -> Vec<R>
where
    R: Send,
    FL: Fn(&Arc<Tape>, Cell) -> Option<R> + Sync,
    FF: Fn(&Arc<Tape>, Cell) -> Option<R> + Sync,
{
    descend_go(tape, ctx, grid, cell, leaf, pad, &|| true, on_leaf, on_full)
}

#[allow(clippy::too_many_arguments)]
fn descend_go<R, FL, FF>(tape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, cell: Cell, leaf: i64, pad: f64, go: &(dyn Fn() -> bool + Sync), on_leaf: &FL, on_full: &FF) -> Vec<R>
where
    R: Send,
    FL: Fn(&Arc<Tape>, Cell) -> Option<R> + Sync,
    FF: Fn(&Arc<Tape>, Cell) -> Option<R> + Sync,
{
    if !grid.meets(cell) || !go() {
        return Vec::new();
    }
    let mut ch = vec![BOTH; tape.len()];
    let mut b = grid.cell_box(cell);
    if pad > 0.0 {
        let g = pad * grid.h;
        b = b.map(|i| I::new(i.lo - g, i.hi + g));
    }
    let r = eval_interval(tape, ctx, b, Some(&mut ch));
    // Skipped where no surface meets the closed cell: one whose minimum is
    // exactly 0 may hold surface on a face, which meshing and the
    // optimizers' samples need.
    if r.clear_of_surface() {
        return Vec::new();
    }
    if r.all_inside() {
        return on_full(tape, cell).into_iter().collect();
    }
    let sub = if ch.iter().any(|&c| c != BOTH) { Arc::new(tape.simplify(&ch)) } else { tape.clone() };
    if cell.size <= leaf {
        return on_leaf(&sub, cell).into_iter().collect();
    }
    let parts: Vec<Vec<R>> = cell.children().par_iter().map(|k| descend_go(&sub, ctx, grid, *k, leaf, pad, go, on_leaf, on_full)).collect();
    parts.into_iter().flatten().collect()
}

/// A tape simplified over each block of space a point is evaluated in,
/// built on first use and kept: evaluating many points near a surface costs
/// a fraction of the whole tape there (each block's tape keeps only the
/// branches its range rules cannot rule out, as in `descend`). Blocks nest:
/// each is simplified from the tape of the block twice its size around it,
/// so a block costs a range evaluation of a tape already cut down.
pub struct Local<'a> {
    tape: &'a Arc<Tape>,
    ctx: &'a Ctx,
    origin: V3,
    size: f64,
    blocks: RwLock<HashMap<(u32, [i64; 3]), Arc<Tape>>>,
    /// Each worker thread's own copy of the blocks it has used: many
    /// threads reading one map contend for it, and a block is mostly
    /// reused by the thread that first needed it.
    near: Vec<Mutex<HashMap<(u32, [i64; 3]), Arc<Tape>>>>,
    /// Instructions run on points so far (`work`).
    work: std::sync::atomic::AtomicUsize,
}

/// Levels of nesting: blocks up to 2^20 times the smallest, past any model.
const LEVELS: u32 = 20;

impl<'a> Local<'a> {
    /// The smallest blocks are of `size`, on a lattice through `origin`.
    pub fn new(tape: &'a Arc<Tape>, ctx: &'a Ctx, origin: V3, size: f64) -> Local<'a> {
        let near = (0..rayon::current_num_threads()).map(|_| Mutex::new(HashMap::new())).collect();
        Local { tape, ctx, origin, size, blocks: RwLock::new(HashMap::new()), near, work: Default::default() }
    }

    /// No blocks: every point and box gets the whole tape (for a few
    /// evaluations, where simplifying would cost more than it saves).
    pub fn whole(tape: &'a Arc<Tape>, ctx: &'a Ctx) -> Local<'a> {
        Local { tape, ctx, origin: [0.0; 3], size: f64::INFINITY, blocks: RwLock::new(HashMap::new()), near: Vec::new(), work: Default::default() }
    }

    pub fn ctx(&self) -> &'a Ctx {
        self.ctx
    }

    /// The instructions it has run on points: each point's block's tape's
    /// length, summed (what reading points here cost, not counting the
    /// blocks' simplification).
    pub fn work(&self) -> usize {
        self.work.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Whether it has blocks (else everything gets the whole tape).
    fn blocked(&self) -> bool {
        self.size > 0.0 && self.size.is_finite()
    }

    /// The tape of the smallest block holding the box `b` (the whole tape
    /// where none does): its values agree with the whole tape's in the box,
    /// and its range rules see only the branches that can matter there.
    pub fn tape_over(&self, b: &[I; 3]) -> Arc<Tape> {
        if !self.blocked() || !b.iter().all(|i| i.lo.is_finite() && i.hi.is_finite()) {
            return self.tape.clone();
        }
        let key = |v: f64, a: usize| ((v - self.origin[a]) / self.size).floor() as i64;
        let (klo, khi): ([i64; 3], [i64; 3]) = (std::array::from_fn(|a| key(b[a].lo, a)), std::array::from_fn(|a| key(b[a].hi, a)));
        for level in 0..LEVELS {
            let k = klo.map(|x| x >> level);
            if k != khi.map(|x| x >> level) {
                continue;
            }
            // The keys round: take the block only where its box, as block()
            // simplifies over it, holds b.
            let bb = self.block_box(level, k);
            if (0..3).all(|a| bb[a].lo <= b[a].lo && b[a].hi <= bb[a].hi) {
                return self.block(level, k);
            }
        }
        self.tape.clone()
    }

    /// The field's values at `pts`, as `eval_points` gives them.
    pub fn eval(&self, pts: &[V3]) -> Vec<f64> {
        self.each(pts, |tape, p, out| eval_points(tape, self.ctx, p, out))
    }

    /// Values and spatial gradients at `pts`, as `eval_grads` gives them.
    pub fn grads(&self, pts: &[V3]) -> Vec<[f64; 4]> {
        self.each(pts, |tape, p, out| *out = eval_grads(tape, self.ctx, p))
    }

    /// `run` over the points in each block, with the block's tape: points
    /// are taken in runs of 4096 (callers pass them in an order that
    /// follows space, as a mesh's are), each run grouped by block.
    fn each<R: Copy + Default + Send>(&self, pts: &[V3], run: impl Fn(&Tape, &[V3], &mut Vec<R>) + Sync) -> Vec<R> {
        let key = |p: &V3| -> [i64; 3] { std::array::from_fn(|a| ((p[a] - self.origin[a]) / self.size).floor() as i64) };
        let mut out = vec![R::default(); pts.len()];
        out.par_chunks_mut(4096).zip(pts.par_chunks(4096)).for_each(|(out, pts)| {
            let mut order: Vec<([i64; 3], u32)> = pts.iter().enumerate().map(|(i, p)| (key(p), i as u32)).collect();
            order.sort_unstable();
            let (mut p, mut v) = (Vec::new(), Vec::new());
            for g in order.chunk_by(|a, b| a.0 == b.0) {
                let tape = self.block(0, g[0].0);
                self.work.fetch_add(tape.len() * g.len(), std::sync::atomic::Ordering::Relaxed);
                p.clear();
                p.extend(g.iter().map(|&(_, i)| pts[i as usize]));
                v.clear();
                run(&tape, &p, &mut v);
                for (&(_, i), r) in g.iter().zip(&v) {
                    out[i as usize] = *r;
                }
            }
        });
        out
    }

    /// The tape for block k of size `size * 2^level`.
    fn block(&self, level: u32, k: [i64; 3]) -> Arc<Tape> {
        if level >= LEVELS || !self.blocked() {
            return self.tape.clone();
        }
        let near = rayon::current_thread_index().and_then(|i| self.near.get(i));
        if let Some(t) = near.and_then(|n| n.lock().unwrap().get(&(level, k)).cloned()) {
            return t;
        }
        let t = self.shared_block(level, k);
        if let Some(n) = near {
            n.lock().unwrap().insert((level, k), t.clone());
        }
        t
    }

    fn shared_block(&self, level: u32, k: [i64; 3]) -> Arc<Tape> {
        if let Some(t) = self.blocks.read().unwrap().get(&(level, k)) {
            return t.clone();
        }
        let parent = self.block(level + 1, k.map(|x| x.div_euclid(2)));
        let mut ch = vec![BOTH; parent.len()];
        eval_interval(&parent, self.ctx, self.block_box(level, k), Some(&mut ch));
        let t = if ch.iter().any(|&c| c != BOTH) { Arc::new(parent.simplify(&ch)) } else { parent };
        self.blocks.write().unwrap().entry((level, k)).or_insert(t).clone()
    }

    /// Block k of size `size * 2^level`, grown by a hair so points rounded
    /// onto its faces are in.
    fn block_box(&self, level: u32, k: [i64; 3]) -> [I; 3] {
        let size = self.size * (1u64 << level) as f64;
        let pad = 1e-9 * size;
        std::array::from_fn(|a| {
            let lo = self.origin[a] + size * k[a] as f64;
            I::new(lo - pad - 1e-15 * lo.abs(), lo + size + pad + 1e-15 * lo.abs())
        })
    }
}
