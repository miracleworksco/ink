use super::brick::{B, N};
use crate::octree::Grid;
use ink_core::V3;
use std::collections::HashMap;

/// What a mesh may sample (`mesh_checked`): its grid, in bricks of B^3
/// cells (the samples of as many: past that it is not made), and its search for what the grid
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
/// not rule it out, down to samples `searched` apart; crossed cells left
/// untried, `untried`, where no proof could change what is reported but
/// its counts, Proofs::untried, so `undecided` and `seams` are at least
/// the cells no proof decides); of them, those where
/// material was found among corners outside, those undecided, and those
/// within a cell of one the surface crosses (`beside`), and the undecided
/// ones not (`away`); the thinnest the material found is across (0 where
/// there is none), a point in one (a sample of the other sign where one was
/// found), what sampling it cost (`cost`, in triangles: a brick's samples
/// BRICK_TRIANGLES, a sample of the search as much as one of a brick's),
/// the res it was asked for where it is coarser (`refused`, set by
/// its maker: its grid there took more bricks than it was allowed; 0:
/// none), and the machine instructions its field's evaluations retired
/// (Tape::work): over boxes, its cells' proofs and its search's, and the
/// tapes they simplified (`read`), and at points, its grid's, crossings',
/// quads' diagonals' and search's samples and gradients (`sampled`); and
/// the tapes its proofs and search made, each the field simplified over a
/// box (`tapes`).
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
    pub tapes: usize,
    pub untried: usize,
    /// Regions drawn whole (`Flat`: not sampled inside their rims), and the
    /// farthest their triangles may lie from the surface.
    pub whole: usize,
    pub whole_within: f64,
}

/// Where a mesh leaves something out (`Unresolved::regions`): a box holding
/// every cell of it (a cover: much of the box may be drawn), how many cells
/// it holds, and of them those where material or a gap was found and those
/// undecided (and `untried`), and those beside a cell the surface crosses; the coarsest
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
    pub untried: usize,
}

/// The most regions a mesh reports: past it, touching groups are merged,
/// by the octree's blocks, into covers.
pub const MAX_REGIONS: usize = 16;

/// Unresolved cells by brick (its index in B^3 cells), each a region.
#[derive(Default)]
pub(super) struct Regions(HashMap<[i64; 3], Region>);

impl Regions {
    /// Cell `c` of `grid`, counted and bounded in its brick's region, and
    /// `f` of it.
    pub(super) fn add(&mut self, grid: &Grid, c: [i64; 3], f: impl FnOnce(&mut Region)) {
        let (lo, hi) = (grid.pos(c), grid.pos(c.map(|x| x + 1)));
        let r = self.0.entry(c.map(|x| x.div_euclid(B))).or_insert_with(|| Region { lo, hi, groups: 1, ..Default::default() });
        r.lo = std::array::from_fn(|a| r.lo[a].min(lo[a]));
        r.hi = std::array::from_fn(|a| r.hi[a].max(hi[a]));
        r.cells += 1;
        f(r);
    }

    /// The regions: bricks that touch joined, then, past `most`, groups
    /// merged by ever larger blocks of bricks until no more than `most`.
    pub(super) fn merged(self, most: usize) -> Vec<Region> {
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
            (a.cells, a.material, a.gaps, a.undecided, a.beside, a.seams, a.untried) =
                (a.cells + b.cells, a.material + b.material, a.gaps + b.gaps, a.undecided + b.undecided, a.beside + b.beside, a.seams + b.seams, a.untried + b.untried);
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

/// The machine instructions a mesh's evaluations of its field retired over
/// boxes and at points, and the tapes it made for its proofs and search
/// (Unresolved::read, sampled, tapes).
#[derive(Default)]
pub(super) struct Ran {
    pub(super) boxes: std::sync::atomic::AtomicUsize,
    pub(super) points: std::sync::atomic::AtomicUsize,
    pub(super) tapes: std::sync::atomic::AtomicUsize,
}
