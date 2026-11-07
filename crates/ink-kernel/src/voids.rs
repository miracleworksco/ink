//! A shape's enclosed voids, counted by proof: how many components of its
//! air (f > 0) no path leads out of, or bounds on how many.
//!
//! An octree of interval bounds sorts space into cells with no air (f <= 0
//! over the cell), cells all air (f > 0), and cells the surface crosses
//! where the field rises along an axis (its partial's bound has one sign):
//! there every air point slides along that axis, staying in air, to the
//! cell's face, into the cell across it. Air cells joined by a face share
//! their air; a sloped cell's air joins the air of the cell its slope leads
//! to (one at least its size, so its whole face lies in it), and has none
//! where that chain ends in a cell with none. Each void holds such a group,
//! so the groups that do not reach the box's edge bound the voids from
//! above. Each set of cells joined across faces not proven airless, which
//! does not reach the box's edge and holds an air cell, is sealed by cells
//! with no air: it holds a void, and bounds them from below.
//!
//! Where the two differ, the cells that keep them apart are halved: a cell
//! whose slope leads nowhere (round a loop, into smaller cells), one with
//! no slope, two of different groups that a face joins (more halving finds
//! the air that joins them, or the airless wall between). A wall thinner
//! than the cells is found where halving reaches it: no sample decides.
//! The halving is the one search's (crate::search): each pass searches
//! the cells to halve, settling a cell its bounds sort and setting aside
//! one at the finest cells, under one budget of cells tested and the
//! call's time limit. Where the finest cells, the budget or the time limit
//! leave the bounds apart, the count is bounds, and says which.

use crate::octree::{Cell, Grid};
use crate::query::Ended;
use crate::search::{BATCH, Node, Outcome, Precision, Search, Smear, Test};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_interval_grad};
use ink_core::interval::I;
use ink_core::tape::{BOTH, Tape};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

/// The most cells a count tests (its passes together): past it, or past
/// the search's estimate of four times as many, it answers with what it
/// has (about 1.5 GB, a few seconds).
pub const MAX_CELLS: usize = 1_000_000;

/// Proven bounds on a shape's enclosed voids: at least `at_least`, at most
/// `at_most` (None: unbounded, where cells it could not decide may each
/// hold one). `at`: where one of the cells that leave it undecided is,
/// else a point in a void found. `cells`: the cells it used.
#[derive(Clone, Copy, Debug, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Voids {
    pub at_least: usize,
    pub at_most: Option<usize>,
    pub at: Option<V3>,
    pub cells: usize,
    /// The finest cells it may halve to (the res asked, rounded down to a
    /// power of two).
    pub finest: f64,
    /// The sides of the cells that leave it undecided, smallest and
    /// largest (None: decided).
    pub undecided: Option<(f64, f64)>,
    /// Why it stopped before halving every cell it needed to as far as
    /// `finest` (Ended::Limit: MAX_CELLS; Ended::Stopped: the call's time
    /// limit). None: it did, or decided.
    pub ended: Option<Ended>,
}

impl Voids {
    /// The count, where its bounds meet.
    pub fn count(&self) -> Option<usize> {
        (self.at_most == Some(self.at_least)).then_some(self.at_least)
    }

    /// Where it is undecided: what its cells left undecided, and why, and
    /// what would decide it. None: decided.
    pub fn why_undecided(&self) -> Option<(String, String)> {
        if self.count().is_some() {
            return None;
        }
        let (a, b) = self.undecided.unwrap_or((self.finest, self.finest));
        let sides = if a == b { format!("{a} mm across") } else { format!("{a} to {b} mm across") };
        let stuck = "touching surfaces and flush faces (f = 0 across a mouth) stay undecided at any res: move them apart or overlap them";
        Some(match self.ended {
            Some(Ended::Limit) => (
                format!("voids in cells it could not decide ({sides}) when it stopped to stay within its budget of {MAX_CELLS} cells, before halving them to {} mm", self.finest),
                format!("a finer res will not help: count a smaller region, (max f (box ...)) round where it is undecided, or simplify the shape there; {stuck}"),
            ),
            Some(Ended::Stopped) => (
                format!("voids in cells it could not decide ({sides}) when the call's time limit passed, before halving them to {} mm", self.finest),
                "give the call more time, or count a smaller region".into(),
            ),
            _ => (
                format!("voids in cells at its finest ({sides}) it could not decide: a wall or gap thinner than them, surfaces that touch, or a face flush with another"),
                format!("a finer (res r), under the thinnest wall or gap there; {stuck}"),
            ),
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Kind {
    /// No air: f <= 0 over the cell.
    Airless,
    /// All air: f > 0.
    Air,
    /// f does not fall along axis .0 in direction .1 (+1 or -1).
    Slope(usize, i64),
    /// None of these, where the search left it.
    Open,
}

struct Leaf {
    cell: Cell,
    kind: Kind,
    /// The field over the cell: the tape it was tested with (its parent's,
    /// shared with its siblings: one of its own per cell would hold
    /// gigabytes).
    tape: Arc<Tape>,
}

/// Where a cell's air is: in the group of air cell .0, none, or unknown.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Lands {
    Group(usize),
    None,
    Open,
}

/// The enclosed voids of the field `tape` whose shape lies in `lo..hi`
/// (f > 0 everywhere outside it), its cells halved down to `h` (rounded
/// down to a power of two, on a grid of its multiples, so the search's
/// halves are the grid's cells exactly).
pub fn voids(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64) -> Voids {
    let h = 2f64.powi(h.log2().floor() as i32);
    // A cell's pad: the root's faces lie outside lo..hi, in air.
    let origin = lo.map(|x| ((x - h) / h).floor() * h);
    let n = (0..3).map(|k| ((hi[k] + h - origin[k]) / h).ceil() as u64).max().unwrap_or(1).max(1);
    let root = Cell { i: [0; 3], size: n.next_power_of_two() as i64 };
    let grid = Grid::new(origin, h);
    let ctx = Arc::new(ctx.clone());
    let first = Node { b: grid.cell_box(root), ctx: ctx.clone(), state: tape.clone() };
    let (test, leaf) = sort(&first, &grid);
    let mut leaves: Vec<Leaf> = leaf.into_iter().collect();
    let (mut pass, mut tested) = (vec![(first, test)], 1);
    loop {
        let (more, used, ended) = search(pass, &grid, h, MAX_CELLS.saturating_sub(tested));
        leaves.extend(more);
        tested += used;
        let (found, halve) = count(&ctx, &grid, root, &leaves);
        if halve.is_empty() || ended.is_some() || tested >= MAX_CELLS {
            let ended = ended.or((!halve.is_empty()).then_some(Ended::Limit));
            return Voids { cells: leaves.len(), finest: h, ended: ended.filter(|_| found.count().is_none()), ..found };
        }
        // The cells to halve, out of the leaves, searched again.
        let mut keep = vec![true; leaves.len()];
        halve.iter().for_each(|&i| keep[i] = false);
        let smear = |c: Cell| Smear { space: [h * c.size as f64; 3], ..Default::default() };
        pass = halve
            .iter()
            .map(|&i| {
                (
                    Node { b: grid.cell_box(leaves[i].cell), ctx: ctx.clone(), state: leaves[i].tape.clone() },
                    Test::Open { bound: 0.0, smear: smear(leaves[i].cell), state: leaves[i].tape.clone() },
                )
            })
            .collect();
        leaves = leaves.into_iter().zip(keep).filter_map(|(l, k)| k.then_some(l)).collect();
    }
}

/// One pass of the search over `roots` (tested), at most `limit` cells:
/// the leaves it sorted, those it left (Open), how many it tested, and why
/// it stopped short (its budget, the time limit).
fn search(roots: Vec<(Node<Arc<Tape>>, Test<Arc<Tape>, ()>)>, grid: &Grid, h: f64, limit: usize) -> (Vec<Leaf>, usize, Option<Ended>) {
    let mut s = Search::start(roots, Precision { finest: h, splits: 0, maps: true });
    let mut leaves = Vec::new();
    let outcome = loop {
        let mut sorted = Vec::new();
        let round = s.round_by(limit, BATCH, |nodes| {
            let tests: Vec<(Test<Arc<Tape>, ()>, Option<Leaf>)> = nodes.par_iter().map(|n| sort(n, grid)).collect();
            tests
                .into_iter()
                .map(|(t, l)| {
                    sorted.extend(l);
                    t
                })
                .collect()
        });
        leaves.extend(sorted);
        if let Some(o) = round {
            break o;
        }
    };
    let ended = match outcome {
        Outcome::Unknown(u) if u.ended != Ended::Finest => Some(u.ended),
        _ => None,
    };
    leaves.extend(s.left().map(|n| Leaf { cell: cell_of(grid, n.b), kind: Kind::Open, tape: n.state.clone() }));
    (leaves, s.tested(), ended)
}

/// A box of the search: the grid's cell.
fn cell_of(grid: &Grid, b: [I; 3]) -> Cell {
    Cell { i: std::array::from_fn(|k| ((b[k].lo - grid.origin[k]) / grid.h).round() as i64), size: ((b[0].hi - b[0].lo) / grid.h).round() as i64 }
}

/// A box tested: settled where its bounds sort it (with its leaf), else
/// open, every axis alike (it is halved into eight), with its tape
/// simplified over it for its halves.
fn sort(n: &Node<Arc<Tape>>, grid: &Grid) -> (Test<Arc<Tape>, ()>, Option<Leaf>) {
    let (tape, b) = (&n.state, n.b);
    let leaf = |kind: Kind| (Test::Settled, Some(Leaf { cell: cell_of(grid, b), kind, tape: tape.clone() }));
    let mut ch = vec![BOTH; tape.len()];
    let r = eval_interval(tape, &n.ctx, b, Some(&mut ch));
    if r.hi <= 0.0 {
        return leaf(Kind::Airless);
    }
    if r.lo > 0.0 {
        return leaf(Kind::Air);
    }
    let sub = if ch.iter().any(|&x| x != BOTH) { Arc::new(tape.simplify(&ch)) } else { tape.clone() };
    if let Some(k) = slope(&sub, &n.ctx, b) {
        return leaf(k);
    }
    let w = b[0].hi - b[0].lo;
    (Test::Open { bound: 0.0, smear: Smear { space: [w; 3], ..Default::default() }, state: sub }, None)
}

/// The axis and direction the field surely does not fall along over box
/// `b`: its steepest such, at least somewhere (not along an edge it is
/// flat along), then the one rising, so two cells flat along an axis do
/// not lead into each other.
fn slope(t: &Tape, ctx: &Ctx, b: [I; 3]) -> Option<Kind> {
    let g = eval_interval_grad(t, ctx, b);
    (0..3)
        .flat_map(|k| [(k, -1, -g.d[k].hi, -g.d[k].lo), (k, 1, g.d[k].lo, g.d[k].hi)])
        .filter(|s| s.2 >= 0.0)
        .max_by(|a, b| a.2.partial_cmp(&b.2).expect("not NaN").then(a.3.total_cmp(&b.3)).then(a.1.cmp(&b.1)))
        .map(|(k, s, ..)| Kind::Slope(k, s))
}

fn find(p: &mut [usize], mut x: usize) -> usize {
    while p[x] != x {
        p[x] = p[p[x]];
        x = p[x];
    }
    x
}

fn union(p: &mut [usize], a: usize, b: usize) {
    let (a, b) = (find(p, a), find(p, b));
    p[a] = b;
}

fn centre(grid: &Grid, c: Cell) -> V3 {
    let (lo, hi) = grid.cell_bounds(c);
    std::array::from_fn(|k| 0.5 * (lo[k] + hi[k]))
}

/// The bounds the leaves prove, and the leaves to halve to tighten them.
fn count(ctx: &Ctx, grid: &Grid, root: Cell, leaves: &[Leaf]) -> (Voids, Vec<usize>) {
    let n = leaves.len();
    let index: HashMap<([i64; 3], i64), usize> = leaves.iter().enumerate().map(|(i, l)| ((l.cell.i, l.cell.size), i)).collect();
    // The leaf holding grid cell p (None: outside the root).
    let at = |p: [i64; 3]| -> Option<usize> {
        if p.iter().any(|&x| x < 0 || x >= root.size) {
            return None;
        }
        let mut s = 1;
        while s <= root.size {
            if let Some(&i) = index.get(&(p.map(|x| x.div_euclid(s) * s), s)) {
                return Some(i);
            }
            s *= 2;
        }
        None
    };
    // Past the root is air (the shape lies in it, a cell in): OUTSIDE's.
    const OUTSIDE: usize = usize::MAX;
    let across = |c: Cell, k: usize, s: i64| {
        let mut p = c.i;
        p[k] = if s > 0 { p[k] + c.size } else { p[k] - 1 };
        p
    };
    let next: Vec<Option<usize>> = leaves
        .par_iter()
        .map(|l| match l.kind {
            Kind::Slope(k, s) => match at(across(l.cell, k, s)) {
                None => Some(OUTSIDE),
                Some(j) => (leaves[j].cell.size >= l.cell.size).then_some(j),
            },
            _ => None,
        })
        .collect();
    let mut lands: Vec<Option<Lands>> = vec![None; n];
    let mut path: Vec<usize> = Vec::new();
    let mut on_path = vec![false; n];
    for i in 0..n {
        let mut j = i;
        let end = loop {
            if let Some(e) = lands[j] {
                break e;
            }
            if on_path[j] {
                break Lands::Open;
            }
            match leaves[j].kind {
                Kind::Air => break Lands::Group(j),
                Kind::Airless => break Lands::None,
                Kind::Open => break Lands::Open,
                Kind::Slope(..) => match next[j] {
                    Some(OUTSIDE) => break Lands::Group(n),
                    Some(t) => {
                        on_path[j] = true;
                        path.push(j);
                        j = t;
                    }
                    None => break Lands::Open,
                },
            }
        };
        if lands[j].is_none() && !on_path[j] {
            lands[j] = Some(end);
        }
        for p in path.drain(..) {
            on_path[p] = false;
            lands[p] = Some(end);
        }
    }
    let lands: Vec<Lands> = lands.into_iter().map(|l| l.expect("every leaf resolved")).collect();
    // Faces between leaves where air may pass: from the smaller of each
    // two (its whole face in the other), each once.
    let carries = |i: usize| lands[i] != Lands::None;
    let faces: Vec<(usize, usize, bool, bool)> = (0..n)
        .into_par_iter()
        .flat_map_iter(|i| {
            let l = &leaves[i];
            let mut out = Vec::new();
            if !carries(i) {
                return out;
            }
            for k in 0..3 {
                for s in [-1i64, 1] {
                    let Some(j) = at(across(l.cell, k, s)) else { continue };
                    let m = &leaves[j];
                    if !carries(j) || m.cell.size < l.cell.size || (m.cell.size == l.cell.size && s < 0) {
                        continue;
                    }
                    let mut b = grid.cell_box(l.cell);
                    b[k] = I::point(if s > 0 { b[k].hi } else { b[k].lo });
                    let r = eval_interval(&l.tape, ctx, b, None);
                    if r.hi <= 0.0 {
                        continue;
                    }
                    let air = l.kind == Kind::Air || m.kind == Kind::Air || r.lo > 0.0 || {
                        let mid = b.map(|x| I::point(x.mid()));
                        eval_interval(&l.tape, ctx, mid, None).lo > 0.0
                    };
                    out.push((i, j, air, true));
                }
            }
            out
        })
        .collect();
    // Groups (by air leaf; n is the outside), and the sets faces join.
    let mut group: Vec<usize> = (0..=n).collect();
    let mut joined: Vec<usize> = (0..=n).collect();
    let edge = |c: Cell| (0..3).any(|k| c.i[k] == 0 || c.i[k] + c.size == root.size);
    for (i, l) in leaves.iter().enumerate() {
        if edge(l.cell) {
            if let Lands::Group(g) = lands[i] {
                union(&mut group, g, n);
            }
            if carries(i) {
                union(&mut joined, i, n);
            }
        }
    }
    for &(i, j, air, _) in &faces {
        union(&mut joined, i, j);
        if let (true, Lands::Group(a), Lands::Group(b)) = (air, lands[i], lands[j]) {
            union(&mut group, a, b);
        }
    }
    let mut apart: Vec<usize> = Vec::new();
    for &(i, j, _, _) in &faces {
        if let (Lands::Group(a), Lands::Group(b)) = (lands[i], lands[j])
            && find(&mut group, a) != find(&mut group, b)
        {
            apart.extend([i, j]);
        }
    }
    let open: Vec<usize> = (0..n).filter(|&i| lands[i] == Lands::Open).collect();
    let mut halve: Vec<usize> = apart.iter().chain(&open).copied().collect();
    halve.retain(|&i| leaves[i].cell.size > 1);
    halve.sort_unstable();
    halve.dedup();
    let outside = find(&mut group, n);
    let mut groups: Vec<usize> = (0..n).filter(|&i| leaves[i].kind == Kind::Air).map(|i| find(&mut group, i)).filter(|&g| g != outside).collect();
    groups.sort_unstable();
    groups.dedup();
    let sealed = find(&mut joined, n);
    let mut holding: Vec<usize> = (0..n).filter(|&i| leaves[i].kind == Kind::Air).map(|i| find(&mut joined, i)).filter(|&s| s != sealed).collect();
    holding.sort_unstable();
    holding.dedup();
    let void_at = (0..n).find(|&i| leaves[i].kind == Kind::Air && find(&mut group, i) != outside).map(|i| centre(grid, leaves[i].cell));
    let at = open.first().map(|&i| centre(grid, leaves[i].cell)).or(void_at);
    let at_most = open.is_empty().then_some(groups.len());
    let mut undecided = None;
    if at_most != Some(holding.len()) {
        // The cells that keep the bounds apart: those with no group, and
        // two of different groups a face joins.
        let side = |i: usize| grid.h * leaves[i].cell.size as f64;
        let sides = open.iter().chain(&apart).map(|&i| side(i));
        undecided = sides.fold(None, |m: Option<(f64, f64)>, s| Some(m.map_or((s, s), |(a, b)| (a.min(s), b.max(s)))));
    }
    (Voids { at_least: holding.len(), at_most, at, cells: n, finest: grid.h, undecided, ended: None }, halve)
}
