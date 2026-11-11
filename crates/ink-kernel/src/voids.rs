//! A shape's enclosed voids, counted by proof: how many components of its
//! air (where f is not below 0, `ink_core::sign`: where it is 0 too) no
//! path leads out of, or bounds on how many. To an accuracy
//! e, as an emptiness claim's solid is (query::Cover, for the air): a void
//! counts where it holds a point deeper than e; air every point of which
//! lies within e of a point shown not air (a slit, or air along a sheet,
//! thinner than 2e) is set aside (`Voids::thin`). A wall seals however thin
//! it is. At e = 0 every void counts.
//!
//! An octree of interval bounds sorts space into cells with no air (f < 0
//! over the cell), cells all air (f >= 0), and cells the surface crosses
//! whose air slides, staying air, out through faces (`slope`): along an
//! axis the field rises along (its partial's bound has one sign), else,
//! where the field is a max of terms (a cut, an intersection), each term's
//! air along a direction that term keeps it along, shown by the bounds of
//! the terms under it (a min's and a max's each), not their hull. Air cells
//! joined by a face share their air; a sloped cell's air lands where the
//! cells it slides into land (the one at least its size each face lies in,
//! else the smaller ones tiling it): none in a cell with no air, each point
//! in one of their groups where they land in several (the cell is in
//! theirs once they come to be one). Each void holds such a group,
//! so the groups that do not reach the box's edge bound the voids from
//! above, but for a group whose every cell is thin, which holds no point
//! deeper than e.
//! Each set of cells joined across faces not proven airless, which does not
//! reach the box's edge and holds an air cell, is sealed by cells with no
//! air: it holds a void, and bounds them from below; to the accuracy, where
//! the cube of half-side e about one of its air cells' centres is shown air
//! (an emptiness claim's witness).
//!
//! Where the two differ, the cells that keep them apart are halved: a cell
//! whose slope leads nowhere (round a loop, past a thin cell), one with
//! no slope, one landing in several groups, two of different groups that a
//! face joins (more halving finds the air that joins them, or the airless
//! wall between). A thin cell is not halved, unless the bounds stay apart
//! with nothing else to halve (thin air that may lead out of a void, or
//! join two): then those on a shortest way from each group that may be a
//! void to another group or the outside (one way out decides it), else
//! every thin cell, to the finest cells. A wall thinner than the cells is
//! found where halving reaches it:
//! no sample decides. The halving is the one search's (crate::search):
//! each pass searches the cells to halve, settling a cell its bounds sort
//! and setting aside one at the finest cells, under one budget of cells
//! tested and the call's time limit. Where the finest cells, the budget or
//! the time limit leave the bounds apart, the count is bounds, and says
//! which.

use crate::octree::{Cell, Grid};
use crate::query::{Cover, Ended, box_samples};
use crate::search::{BATCH, Node, Outcome, Precision, Search, Smear, Test};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_interval_grad, eval_interval_grad_every, eval_points};
use ink_core::graph::{BinOp, UnOp};
use ink_core::interval::{self, I, norm};
use ink_core::tape::{BOTH, Op, Tape};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

/// The most cells a count tests (its passes together): past it, or where a
/// pass alone is forecast past what is left of it (the search's estimate:
/// the passes after it only add to it), it answers with what it has (about
/// 1.5 GB, a few seconds).
pub const MAX_CELLS: usize = 1_000_000;

/// Proven bounds on a shape's enclosed voids, to an accuracy: at least
/// `at_least`, at most `at_most` (None: unbounded, where cells it could
/// not decide may each hold one). `at`: where one of the cells that leave
/// it undecided is, else a point in a void found. `cells`: the cells it
/// used.
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
    /// The accuracy it counted to: a void holds a point deeper than it (0:
    /// any void).
    pub accuracy: f64,
    /// The cells it set aside as holding no such point (every point of each
    /// within the accuracy of a point shown not air), and where one is.
    pub thin: usize,
    pub thin_at: Option<V3>,
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
        let stuck = "flush mouths whose pieces read no value in common there (two formulas for one face) stay undecided at any res: write the face once, or move them apart or overlap them";
        let thin = match self.accuracy > 0.0 {
            true => format!(", air no deeper than the accuracy ({} mm) that may lead out of a void, a void about as deep as it", self.accuracy),
            false => String::new(),
        };
        Some(match self.ended {
            Some(Ended::Limit) => (
                format!(
                    "voids in cells it could not decide ({sides}) when it stopped to stay within its budget of {MAX_CELLS} cells (or was forecast past it), before halving them to {} mm",
                    self.finest
                ),
                format!("a finer res will not help: count a smaller region, (max f (box ...)) round where it is undecided, or simplify the shape there; {stuck}"),
            ),
            Some(Ended::Stopped) => (
                format!("voids in cells it could not decide ({sides}) when the call's time limit passed, before halving them to {} mm", self.finest),
                "give the call more time, or count a smaller region".into(),
            ),
            _ => (
                format!("voids in cells at its finest ({sides}) it could not decide: a wall or gap thinner than them, a face flush with another{thin}"),
                format!("a finer (res r), under the thinnest wall or gap there; {stuck}"),
            ),
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Kind {
    /// No air: f < 0 over the cell.
    Airless,
    /// All air: f >= 0.
    Air,
    /// Its air slides, staying air, along the directions set: bit 2k + 1
    /// along axis k upward, 2k downward (`slope`).
    Slope(u8),
    /// None of these, every point within the accuracy of a point shown not
    /// air (Ask::thin).
    Thin,
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
    /// Whether it is thin (Ask::thin), once asked: a leaf is kept from
    /// pass to pass.
    thin: std::sync::OnceLock<bool>,
    /// Its faces, once asked (`Leaf::face`), by axis and side.
    faces: [std::sync::OnceLock<Face>; 6],
}

/// What a leaf's face holds: no air, air shown, or air it may.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Face {
    Airless,
    Air,
    May,
}

impl Leaf {
    fn new(cell: Cell, kind: Kind, tape: Arc<Tape>) -> Leaf {
        Leaf { cell, kind, tape, thin: std::sync::OnceLock::new(), faces: Default::default() }
    }

    /// Its face on axis `k`, side `s` (-1 or 1), by bounds over it (an air
    /// cell's all air), else air where its middle is shown air.
    fn face(&self, k: usize, s: i64, ctx: &Ctx, grid: &Grid) -> Face {
        *self.faces[2 * k + (s > 0) as usize].get_or_init(|| {
            if self.kind == Kind::Air {
                return Face::Air;
            }
            let mut b = grid.cell_box(self.cell);
            b[k] = I::point(if s > 0 { b[k].hi } else { b[k].lo });
            let r = eval_interval(&self.tape, ctx, b, None);
            if r.all_inside() {
                Face::Airless
            } else if r.none_inside() || eval_interval(&self.tape, ctx, b.map(|x| I::point(x.mid())), None).none_inside() {
                Face::Air
            } else {
                Face::May
            }
        })
    }

    fn thin(&self, ask: &Ask, ctx: &Ctx, grid: &Grid) -> bool {
        *self.thin.get_or_init(|| self.kind == Kind::Thin || ask.thin(&self.tape, ctx, grid.cell_box(self.cell)))
    }
}

/// A search's box state: the field over it, and whether a thin box
/// settles (not under a thin cell halved).
type State = (Arc<Tape>, bool);

/// Where a cell's air is: in the group of air cell .0, each point in one
/// of the groups of a set of air cells (`count`'s sets, by index), none,
/// unknown, unknown but no deeper than the accuracy, or, where deeper, in
/// the groups of the cells beside it (`count`).
#[derive(Clone, Copy, Debug, PartialEq)]
enum Lands {
    Group(usize),
    Groups(usize),
    None,
    Open,
    Thin,
    Beside,
}

/// What a count asks of a cell beyond its bounds: the whole field (for
/// points past the tape a cell was tested with) and the accuracy, as
/// emptiness's rule for the air (query::Cover, sign -1).
struct Ask {
    whole: Arc<Tape>,
    cover: Cover,
}

impl Ask {
    fn new(whole: &Arc<Tape>, e: f64) -> Ask {
        let all = [I::new(f64::NEG_INFINITY, f64::INFINITY); 3];
        Ask { whole: whole.clone(), cover: Cover { tol: e, region: all, asked: all, sign: -1.0, read: [true; 3], still: [false; 3] } }
    }

    /// Whether every point of box `b` (`t` the field over it) lies within
    /// the accuracy of a point shown not air (Cover::covered): its air is
    /// no deeper, as an emptiness claim's solid is set aside. None of a box
    /// whose 64 sub-boxes are each wider than the accuracy corner to corner
    /// can be (a speed rule: it is halved).
    fn thin(&self, t: &Tape, ctx: &Ctx, b: [I; 3]) -> bool {
        let side = b.map(|x| x.hi - x.lo);
        if !(self.cover.tol > 0.0) || norm(side) > 8.0 * self.cover.tol {
            return false;
        }
        let pts = box_samples(b, 4);
        let mut v = Vec::new();
        eval_points(t, ctx, &pts, &mut v);
        // The field over the box grown by the accuracy, where a point past
        // the box is read.
        let wide = std::cell::OnceCell::new();
        let grown = || {
            wide.get_or_init(|| {
                let mut ch = vec![BOTH; self.whole.len()];
                eval_interval(&self.whole, ctx, self.cover.grown(b), Some(&mut ch));
                self.whole.simplify(&ch)
            })
        };
        self.cover.covered(t, &grown, &|r: I, _: &Ctx, _: [I; 3]| r.all_inside(), ctx, b, side, &pts, &v)
    }

    /// Whether air cell `b` holds a point deeper than the accuracy, as an
    /// emptiness claim's witness is: the cube of half-side the accuracy
    /// about its centre lies in it, or is shown air (its bounds, halved up
    /// to three times where they do not show it).
    fn deep(&self, ctx: &Ctx, b: [I; 3]) -> bool {
        fn air(t: &Tape, ctx: &Ctx, b: [I; 3], depth: u32) -> bool {
            let r = eval_interval(t, ctx, b, None);
            if r.none_inside() {
                return true;
            }
            if depth == 0 || r.all_inside() {
                return false;
            }
            (0..8).all(|o: usize| air(t, ctx, std::array::from_fn(|k| if o >> k & 1 == 0 { I::new(b[k].lo, b[k].mid()) } else { I::new(b[k].mid(), b[k].hi) }), depth - 1))
        }
        let e = self.cover.tol;
        if b[0].hi - b[0].lo >= 2.0 * e {
            return true;
        }
        air(&self.whole, ctx, b.map(|x| I::new((x.mid() - e).next_down(), (x.mid() + e).next_up())), 3)
    }
}

/// The enclosed voids of the field `tape` whose shape lies in `lo..hi`
/// (f > 0 everywhere outside it), its cells halved down to `h` (rounded
/// down to a power of two, on a grid of its multiples, so the search's
/// halves are the grid's cells exactly), to accuracy `e` (a void holds a
/// point deeper than e; 0: any void).
pub fn voids(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64, e: f64) -> Voids {
    let h = 2f64.powi(h.log2().floor() as i32);
    let ask = Ask::new(tape, e);
    // A cell's pad: the root's faces lie outside lo..hi, in air.
    let origin = lo.map(|x| ((x - h) / h).floor() * h);
    let n = (0..3).map(|k| ((hi[k] + h - origin[k]) / h).ceil() as u64).max().unwrap_or(1).max(1);
    let root = Cell { i: [0; 3], size: n.next_power_of_two() as i64 };
    let grid = Grid::new(origin, h);
    let ctx = Arc::new(ctx.clone());
    let first = Node { b: grid.cell_box(root), ctx: ctx.clone(), state: (tape.clone(), true) };
    let (test, leaf) = sort(&first, &grid, &ask);
    let mut leaves: Vec<Leaf> = leaf.into_iter().collect();
    let (mut pass, mut tested) = (vec![(first, test)], 1);
    loop {
        let (more, used, ended) = search(pass, &grid, &ask, h, MAX_CELLS.saturating_sub(tested));
        leaves.extend(more);
        tested += used;
        let (found, halve, settles) = count(&ask, &ctx, &grid, root, &leaves);
        if halve.is_empty() || ended.is_some() || tested >= MAX_CELLS {
            let ended = ended.or((!halve.is_empty()).then_some(Ended::Limit));
            return Voids { cells: leaves.len(), finest: h, ended: ended.filter(|_| found.count().is_none()), accuracy: e, ..found };
        }
        // The cells to halve, out of the leaves, searched again.
        let mut keep = vec![true; leaves.len()];
        halve.iter().for_each(|&i| keep[i] = false);
        let smear = |c: Cell| Smear { space: [h * c.size as f64; 3], ..Default::default() };
        pass = halve
            .iter()
            .map(|&i| {
                let state = (leaves[i].tape.clone(), settles);
                (Node { b: grid.cell_box(leaves[i].cell), ctx: ctx.clone(), state: state.clone() }, Test::Open { bound: 0.0, smear: smear(leaves[i].cell), state })
            })
            .collect();
        leaves = leaves.into_iter().zip(keep).filter_map(|(l, k)| k.then_some(l)).collect();
    }
}

/// One pass of the search over `roots` (tested), at most `limit` cells,
/// stopped where it is forecast past them (MAX_CELLS): the leaves it
/// sorted, those it left (Open), how many it tested, and why it stopped
/// short (its budget, the time limit).
fn search(roots: Vec<(Node<State>, Test<State, ()>)>, grid: &Grid, ask: &Ask, h: f64, limit: usize) -> (Vec<Leaf>, usize, Option<Ended>) {
    let mut s = Search::start(roots, Precision { finest: h, splits: 0, maps: true });
    let mut leaves = Vec::new();
    let ended = loop {
        let mut sorted = Vec::new();
        let round = s.round_by(limit, BATCH, |nodes| {
            let tests: Vec<(Test<State, ()>, Option<Leaf>)> = nodes.par_iter().map(|n| sort(n, grid, ask)).collect();
            tests
                .into_iter()
                .map(|(t, l)| {
                    sorted.extend(l);
                    t
                })
                .collect()
        });
        leaves.extend(sorted);
        match round {
            Some(Outcome::Unknown(u)) if u.ended != Ended::Finest => break Some(u.ended),
            Some(_) => break None,
            // Forecast past what is left of the count's budget, this pass
            // alone: the passes after it only add to it.
            None if s.forecast() > limit as f64 => break Some(Ended::Limit),
            None => {}
        }
    };
    leaves.extend(s.left().map(|n| Leaf::new(cell_of(grid, n.b), Kind::Open, n.state.0.clone())));
    (leaves, s.tested(), ended)
}

/// A box of the search: the grid's cell.
fn cell_of(grid: &Grid, b: [I; 3]) -> Cell {
    Cell { i: std::array::from_fn(|k| ((b[k].lo - grid.origin[k]) / grid.h).round() as i64), size: ((b[0].hi - b[0].lo) / grid.h).round() as i64 }
}

/// A box tested: settled where its bounds sort it (with its leaf), or
/// where it is thin and thin boxes settle there, else open, every axis
/// alike (it is halved into eight), with its tape simplified over it for
/// its halves.
fn sort(n: &Node<State>, grid: &Grid, ask: &Ask) -> (Test<State, ()>, Option<Leaf>) {
    let ((tape, settles), b) = (&n.state, n.b);
    let leaf = |kind: Kind| (Test::Settled, Some(Leaf::new(cell_of(grid, b), kind, tape.clone())));
    let mut ch = vec![BOTH; tape.len()];
    let r = eval_interval(tape, &n.ctx, b, Some(&mut ch));
    if r.all_inside() {
        return leaf(Kind::Airless);
    }
    if r.none_inside() {
        return leaf(Kind::Air);
    }
    let sub = tape.narrowed(&ch);
    if let Some(k) = by_signs(&sub, &n.ctx, b) {
        return leaf(k);
    }
    if let Some(k) = slope(&sub, &n.ctx, b) {
        return leaf(k);
    }
    if *settles && ask.thin(&sub, &n.ctx, b) {
        return leaf(Kind::Thin);
    }
    let w = b[0].hi - b[0].lo;
    (Test::Open { bound: 0.0, smear: Smear { space: [w; 3], ..Default::default() }, state: (sub, *settles) }, None)
}

/// Airless or air over the box by signs (`ink_core::sign::holds_by_cuts`):
/// a flush mouth's parts each shown air one side of its face or the other.
fn by_signs(tape: &Tape, ctx: &ink_core::eval::Ctx, b: [I; 3]) -> Option<Kind> {
    if ink_core::sign::holds_by_cuts(tape, ctx, b, &|r: I| r.all_inside()) {
        return Some(Kind::Airless);
    }
    ink_core::sign::holds_by_cuts(tape, ctx, b, &|r: I| r.none_inside()).then_some(Kind::Air)
}

/// The faces the field's air slides out through over box `b`, staying air
/// (bits as Kind::Slope's): along the axis the field surely does not fall
/// along, else, where the field is a max of terms (a cut, an intersection:
/// its air is theirs together), along a direction for each term with air
/// in the box (`keeps`), out through the faces it points to.
fn slope(t: &Tape, ctx: &Ctx, b: [I; 3]) -> Option<Kind> {
    // An axis the field rises along somewhere and falls along nowhere, else
    // a direction for each term that rises along it somewhere, else an axis
    // the field is flat along, else a direction for each term.
    let flat = match rising(&eval_interval_grad(t, ctx, b)) {
        Some((d, true)) => return Some(Kind::Slope(d)),
        d => d.map(|d| d.0),
    };
    let every = eval_interval_grad_every(t, ctx, b);
    let (mut ways, mut rises, mut terms) = (0u8, true, vec![(t.ssa.len() - 1, 1.0)]);
    while let Some((i, s)) = terms.pop() {
        let ins = t.ssa[i];
        match ins.op {
            Op::Bin(BinOp::Max | BinOp::Emax) if s > 0.0 => terms.extend([(ins.a as usize, s), (ins.b as usize, s)]),
            Op::Bin(BinOp::Min) if s < 0.0 => terms.extend([(ins.a as usize, s), (ins.b as usize, s)]),
            Op::Un(UnOp::Neg) => terms.push((ins.a as usize, -s)),
            Op::Tag(_) => terms.push((ins.a as usize, s)),
            // No air from this term in the box.
            _ if (if s > 0.0 { every[i].v.hi < 0.0 } else { every[i].v.lo > 0.0 }) => {}
            _ => {
                // The axes, steepest by the term's bounds first, as `rising`
                // orders them (one face to slide out through), then its own
                // mean gradient.
                let g = &every[i];
                let mid: V3 = std::array::from_fn(|k| s * 0.5 * (g.d[k].lo + g.d[k].hi));
                let len = norm(mid);
                let axes = (0..6).map(|w| std::array::from_fn(|k| if k == w / 2 { if w % 2 == 1 { 1.0 } else { -1.0 } } else { 0.0 }));
                let mut along: Vec<(I, V3)> =
                    axes.chain((len.is_finite() && len > 0.0).then(|| mid.map(|x| x / len))).map(|d: V3| (falls(g, s, d), d)).filter(|x| x.0.lo >= 0.0).collect();
                let axis = |d: &V3| d.iter().filter(|&&x| x != 0.0).count() == 1;
                along.sort_by(|a, b| {
                    axis(&b.1).cmp(&axis(&a.1)).then(b.0.lo.total_cmp(&a.0.lo)).then(b.0.hi.total_cmp(&a.0.hi)).then(b.1.iter().sum::<f64>().total_cmp(&a.1.iter().sum::<f64>()))
                });
                let Some(&(r, d)) = along.iter().find(|&&(_, d)| keeps(t, &every, i, s, d)) else {
                    return flat.map(Kind::Slope);
                };
                rises &= r.lo > 0.0;
                (0..3).filter(|&k| d[k] != 0.0).for_each(|k| ways |= 1 << (2 * k + (d[k] > 0.0) as usize));
            }
        }
    }
    match flat {
        Some(d) if !rises || ways == 0 => Some(Kind::Slope(d)),
        _ => (ways != 0).then_some(Kind::Slope(ways)),
    }
}

/// Bounds on how value `g` read with sign `s` changes along `d`.
fn falls(g: &ink_core::interval::G, s: f64, d: V3) -> I {
    (0..3).fold(I::point(0.0), |acc, k| interval::bin(BinOp::Add, acc, interval::bin(BinOp::Mul, I::point(s * d[k]), g.d[k]).0).0)
}

/// Whether, over the box, every point where value `i` of `t` read with
/// sign `s` (-1: negated) is at least 0 stays so moved along `d` within
/// it: each term of a min or max (an emax's as its max's) keeps its own
/// such points, as then their max or min does; any other value is of one
/// sign over the box, or does not fall along d (its gradient bounds, in
/// `every`).
fn keeps(t: &Tape, every: &[ink_core::interval::G], i: usize, s: f64, d: V3) -> bool {
    let mut stack = vec![(i, s)];
    while let Some((i, s)) = stack.pop() {
        let ins = t.ssa[i];
        match ins.op {
            Op::Bin(BinOp::Max | BinOp::Min | BinOp::Emax) => stack.extend([(ins.a as usize, s), (ins.b as usize, s)]),
            Op::Un(UnOp::Neg) => stack.push((ins.a as usize, -s)),
            Op::Tag(_) => stack.push((ins.a as usize, s)),
            _ => {
                let g = &every[i];
                let kept = if s > 0.0 { g.v.lo >= 0.0 } else { g.v.hi <= 0.0 };
                if kept {
                    continue;
                }
                if !(falls(g, s, d).lo >= 0.0) {
                    return false;
                }
            }
        }
    }
    true
}

/// The direction (a bit of Kind::Slope's) a value with bounds `g` surely
/// does not fall along: its steepest such, at least somewhere (not along
/// an edge it is flat along), then the one rising, so two cells flat along
/// an axis do not lead into each other; and whether it surely rises.
fn rising(g: &ink_core::interval::G) -> Option<(u8, bool)> {
    (0..3)
        .flat_map(|k| [(k, -1, -g.d[k].hi, -g.d[k].lo), (k, 1, g.d[k].lo, g.d[k].hi)])
        .filter(|s| s.2 >= 0.0)
        .max_by(|a, b| a.2.partial_cmp(&b.2).expect("not NaN").then(a.3.total_cmp(&b.3)).then(a.1.cmp(&b.1)))
        .map(|(k, s, lo, _)| (1 << (2 * k + (s > 0) as usize), lo > 0.0))
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

/// The cells on a shortest way across `faces` (pairs of leaves air may pass
/// between) from each group in `from` to another owner (a group, the
/// outside), through cells `passes` (owned by none): the cells a pass
/// halves to show one leads out. Searched from all the groups at once,
/// each cell reached first from the nearest; where two groups' searches
/// meet, the way joins them.
fn ways(n: usize, faces: &[(usize, usize, bool)], owner: &[Option<usize>], passes: &dyn Fn(usize) -> bool, from: &[usize]) -> Vec<usize> {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(i, j, _) in faces {
        adj[i].push(j);
        adj[j].push(i);
    }
    const NONE: usize = usize::MAX;
    // Each reached cell's group and the cell it was reached from.
    let mut by: Vec<usize> = vec![NONE; n];
    let mut prev: Vec<usize> = vec![NONE; n];
    let mut queue = std::collections::VecDeque::new();
    for i in 0..n {
        if let Some(g) = owner[i].filter(|g| from.contains(g)) {
            by[i] = g;
            queue.push_back(i);
        }
    }
    let mut done: Vec<usize> = Vec::new();
    let mut out: Vec<usize> = Vec::new();
    fn back(mut i: usize, prev: &[usize], owner: &[Option<usize>], out: &mut Vec<usize>) {
        while i != usize::MAX && owner[i].is_none() {
            out.push(i);
            i = prev[i];
        }
    }
    while let Some(i) = queue.pop_front() {
        let g = by[i];
        if done.contains(&g) {
            continue;
        }
        for &j in &adj[i] {
            if done.contains(&g) {
                break;
            }
            match owner[j] {
                Some(o) if o != g => {
                    back(i, &prev, owner, &mut out);
                    done.push(g);
                }
                Some(_) => {}
                None if !passes(j) => {}
                None if by[j] == NONE => {
                    by[j] = g;
                    prev[j] = i;
                    queue.push_back(j);
                }
                None if by[j] != g => {
                    back(i, &prev, owner, &mut out);
                    back(j, &prev, owner, &mut out);
                    done.extend([g, by[j]]);
                }
                None => {}
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn centre(grid: &Grid, c: Cell) -> V3 {
    let (lo, hi) = grid.cell_bounds(c);
    std::array::from_fn(|k| 0.5 * (lo[k] + hi[k]))
}

/// The bounds the leaves prove, the leaves to halve to tighten them, and
/// whether thin boxes settle under them (not where they are thin cells,
/// halved where nothing else is left to halve).
fn count(ask: &Ask, ctx: &Ctx, grid: &Grid, root: Cell, leaves: &[Leaf]) -> (Voids, Vec<usize>, bool) {
    let n = leaves.len();
    let index: HashMap<([i64; 3], i64), usize> = leaves.iter().enumerate().map(|(i, l)| ((l.cell.i, l.cell.size), i)).collect();
    // The leaf holding grid cell p, of side `from` at least (None: outside
    // the root, or a smaller one).
    let at_least = |p: [i64; 3], from: i64| -> Option<usize> {
        if p.iter().any(|&x| x < 0 || x >= root.size) {
            return None;
        }
        let mut s = from;
        while s <= root.size {
            if let Some(&i) = index.get(&(p.map(|x| x.div_euclid(s) * s), s)) {
                return Some(i);
            }
            s *= 2;
        }
        None
    };
    let at = |p: [i64; 3]| at_least(p, 1);
    // Past the root is air (the shape lies in it, a cell in): OUTSIDE's.
    const OUTSIDE: usize = usize::MAX;
    let across = |c: Cell, k: usize, s: i64| {
        let mut p = c.i;
        p[k] = if s > 0 { p[k] + c.size } else { p[k] - 1 };
        p
    };
    // The leaves a slope cell's air passes into across each face it slides
    // toward: the one holding the cell of its size beside it, or the
    // smaller ones tiling that cell's side on the face (each point of the
    // face is in one of them); OUTSIDE past the root. None where no leaf
    // holds part of it.
    let next: Vec<Option<Vec<usize>>> = leaves
        .par_iter()
        .map(|l| {
            let Kind::Slope(ways) = l.kind else { return None };
            let mut out = Vec::new();
            for (k, s) in (0..6).filter(|w| ways >> w & 1 == 1).map(|w| (w / 2, if w % 2 == 1 { 1 } else { -1 })) {
                let p = across(l.cell, k, s);
                if p.iter().any(|&x| x < 0 || x >= root.size) {
                    out.push(OUTSIDE);
                    continue;
                }
                if let Some(i) = at_least(p, l.cell.size) {
                    out.push(i);
                    continue;
                }
                let mut beside = l.cell;
                beside.i[k] += s * l.cell.size;
                let mut cells = vec![beside];
                while let Some(c) = cells.pop() {
                    if let Some(&i) = index.get(&(c.i, c.size)) {
                        out.push(i);
                        continue;
                    }
                    if c.size == 1 {
                        return None;
                    }
                    let half = c.size / 2;
                    for o in 0..8i64 {
                        // The halves on the face: the near side along k.
                        if (o >> k & 1) != (s < 0) as i64 {
                            continue;
                        }
                        cells.push(Cell { i: std::array::from_fn(|a| c.i[a] + half * (o >> a & 1)), size: half });
                    }
                }
            }
            let mut seen = std::collections::HashSet::new();
            out.retain(|&i| seen.insert(i));
            Some(out)
        })
        .collect();
    // Where each leaf's air lands, following slopes. A slope cell's air,
    // sliding across its faces, lands where the leaves it passes into do:
    // none in an airless one; past a thin one it is open (unknown); into
    // air cells, each point in one of their groups.
    let mut lands: Vec<Option<Lands>> = leaves.iter().map(|l| (l.kind == Kind::Thin).then_some(Lands::Thin)).collect();
    let mut on_path = vec![false; n];
    let mut sets: Vec<Vec<usize>> = Vec::new();
    let mut join = |acc: Lands, e: Lands| {
        let cells = |l: Lands, sets: &Vec<Vec<usize>>| match l {
            Lands::Group(a) => Some(vec![a]),
            Lands::Groups(s) => Some(sets[s].clone()),
            _ => None,
        };
        match (acc, e) {
            (a, Lands::None) => a,
            (Lands::None, e) => e,
            (Lands::Group(a), Lands::Group(b)) if a == b => acc,
            _ => match (cells(acc, &sets), cells(e, &sets)) {
                (Some(mut a), Some(b)) => {
                    a.extend(b);
                    a.sort_unstable();
                    a.dedup();
                    sets.push(a);
                    Lands::Groups(sets.len() - 1)
                }
                _ => Lands::Open,
            },
        }
    };
    for i in 0..n {
        if lands[i].is_some() {
            continue;
        }
        // Depth first: each cell on the path with the next of its leaves
        // to follow, and where those followed land so far.
        let mut path: Vec<(usize, usize, Lands)> = vec![(i, 0, Lands::None)];
        on_path[i] = true;
        let mut done: Option<Lands> = None;
        while let Some(&(j, at, acc)) = path.last() {
            // The leaf it followed last has landed: joined in, and the
            // next one is followed.
            let (at, acc) = match done.take() {
                Some(e) => (at + 1, join(acc, if e == Lands::Thin { Lands::Open } else { e })),
                None => (at, acc),
            };
            *path.last_mut().unwrap() = (j, at, acc);
            let end = match (leaves[j].kind, &next[j]) {
                (Kind::Air, _) => Some(Lands::Group(j)),
                (Kind::Airless, _) => Some(Lands::None),
                (Kind::Open | Kind::Thin, _) | (Kind::Slope(..), None) => Some(Lands::Open),
                (Kind::Slope(..), Some(v)) if acc == Lands::Open || at == v.len() => Some(acc),
                (Kind::Slope(..), Some(_)) => None,
            };
            if let Some(e) = end {
                path.pop();
                on_path[j] = false;
                lands[j] = Some(e);
                done = Some(e);
                continue;
            }
            let t = next[j].as_ref().unwrap()[at];
            done = match t {
                OUTSIDE => Some(Lands::Group(n)),
                t if on_path[t] => Some(Lands::Open),
                t => lands[t],
            };
            if done.is_none() {
                path.push((t, 0, Lands::None));
                on_path[t] = true;
            }
        }
    }
    // A cell whose air's group is unknown holds no void's deep point where it is
    // thin.
    let lands: Vec<Lands> = lands
        .into_par_iter()
        .enumerate()
        .map(|(i, l)| match l.expect("every leaf resolved") {
            Lands::Open if leaves[i].thin(ask, ctx, grid) => Lands::Thin,
            l => l,
        })
        .collect();
    // A cell no wider corner to corner than the accuracy holding a point
    // deeper than it lies in that point's ball: all air, its air joins the
    // air of each cell across its faces. So beside a cell with no air it
    // holds none (thin), and beside cells with groups, any it holds is in
    // theirs (each counted where it may be a void's).
    let small = |i: usize| grid.h * leaves[i].cell.size as f64 * 3f64.sqrt() <= ask.cover.tol;
    let beside: Vec<(usize, Vec<usize>)> = (0..n)
        .into_par_iter()
        .filter(|&i| lands[i] == Lands::Open && small(i))
        .map(|i| {
            let near = (0..3).flat_map(|k| [-1i64, 1].map(|s| at(across(leaves[i].cell, k, s)).map_or(Lands::Group(n), |j| lands[j])));
            let near: Vec<Lands> = near.collect();
            let groups = near.iter().filter_map(|l| if let Lands::Group(a) = l { Some(*a) } else { None }).collect();
            (i, (!near.contains(&Lands::None)).then_some(groups))
        })
        .filter_map(|(i, g): (usize, Option<Vec<usize>>)| if g.as_ref().is_some_and(|g| g.is_empty()) { None } else { Some((i, g.unwrap_or_default())) })
        .collect();
    let mut lands = lands;
    for (i, g) in &beside {
        lands[*i] = if g.is_empty() { Lands::Thin } else { Lands::Beside };
    }
    // Faces between leaves where air may pass: from the smaller of each
    // two (its whole face in the other), each once.
    let carry: Vec<bool> = lands.iter().map(|&l| l != Lands::None).collect();
    let carries = |i: usize| carry[i];
    let faces: Vec<(usize, usize, bool)> = (0..n)
        .into_par_iter()
        .flat_map_iter(|i| {
            let l = &leaves[i];
            let mut out = Vec::new();
            if !carries(i) {
                return out;
            }
            for k in 0..3 {
                for s in [-1i64, 1] {
                    let Some(j) = at_least(across(l.cell, k, s), l.cell.size) else { continue };
                    let m = &leaves[j];
                    if !carries(j) || (m.cell.size == l.cell.size && s < 0) {
                        continue;
                    }
                    match l.face(k, s, ctx, grid) {
                        Face::Airless => {}
                        face => out.push((i, j, face == Face::Air || m.kind == Kind::Air)),
                    }
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
        if edge(l.cell) && carries(i) {
            union(&mut joined, i, n);
        }
    }
    for &(i, j, _) in &faces {
        union(&mut joined, i, j);
    }
    // A cell landing in several groups that come to be one lands in it.
    let mut several: Vec<usize> = (0..n).filter(|&i| matches!(lands[i], Lands::Groups(_))).collect();
    loop {
        for (i, l) in leaves.iter().enumerate() {
            if let (true, Lands::Group(g)) = (edge(l.cell), lands[i]) {
                union(&mut group, g, n);
            }
        }
        for &(i, j, air) in &faces {
            if let (true, Lands::Group(a), Lands::Group(b)) = (air, lands[i], lands[j]) {
                union(&mut group, a, b);
            }
        }
        let before = several.len();
        several.retain(|&i| {
            let Lands::Groups(s) = lands[i] else { unreachable!() };
            let r = find(&mut group, sets[s][0]);
            let one = sets[s].iter().all(|&a| find(&mut group, a) == r);
            if one {
                lands[i] = Lands::Group(r);
            }
            !one
        });
        if several.len() == before {
            break;
        }
    }
    // Those still in several: more halving may show where their air is.
    let mut apart: Vec<usize> = several.clone();
    for &(i, j, _) in &faces {
        if let (Lands::Group(a), Lands::Group(b)) = (lands[i], lands[j])
            && find(&mut group, a) != find(&mut group, b)
        {
            apart.extend([i, j]);
        }
    }
    let outside = find(&mut group, n);
    // A group that may be a void: one of its cells is not thin.
    let root_of: Vec<Option<usize>> = (0..n)
        .map(|i| match lands[i] {
            Lands::Group(a) => Some(find(&mut group, a)).filter(|&g| g != outside),
            _ => None,
        })
        .collect();
    let mut wide = vec![false; n + 1];
    for (_, near) in &beside {
        let roots: Vec<usize> = near.iter().map(|&a| find(&mut group, a)).collect();
        if !roots.contains(&outside) {
            roots.into_iter().for_each(|g| wide[g] = true);
        }
    }
    // Each other group's cells, the widest first, until one is not thin.
    let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        if let Some(g) = root_of[i].filter(|&g| !wide[g]) {
            members.entry(g).or_default().push(i);
        }
    }
    // A cell in several groups may hold a deep point of any of theirs.
    for &i in &several {
        let Lands::Groups(s) = lands[i] else { unreachable!() };
        let mut roots: Vec<usize> = sets[s].iter().map(|&a| find(&mut group, a)).filter(|&g| g != outside && !wide[g]).collect();
        roots.sort_unstable();
        roots.dedup();
        roots.into_iter().for_each(|g| members.entry(g).or_default().push(i));
    }
    let mut members: Vec<(usize, Vec<usize>)> = members.into_iter().collect();
    members.par_iter_mut().for_each(|(_, m)| m.sort_by_key(|&i| (std::cmp::Reverse(leaves[i].cell.size), i)));
    let found: Vec<(usize, bool)> = members.par_iter().map(|(g, m)| (*g, m.iter().any(|&i| !leaves[i].thin(ask, ctx, grid)))).collect();
    found.iter().for_each(|&(g, w)| wide[g] |= w);
    let thin_groups: Vec<usize> = members.iter().filter(|(g, _)| !wide[*g]).flat_map(|(_, m)| m.iter().copied()).collect();
    let mut groups: Vec<usize> = root_of.iter().flatten().copied().filter(|&g| wide[g]).collect();
    groups.sort_unstable();
    groups.dedup();
    // The sealed sets holding air, and those shown to hold a deep point.
    let sealed = find(&mut joined, n);
    let mut holding: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in (0..n).filter(|&i| leaves[i].kind == Kind::Air) {
        let s = find(&mut joined, i);
        if s != sealed {
            holding.entry(s).or_default().push(i);
        }
    }
    let mut unshown: Vec<usize> = Vec::new();
    let mut at_least = 0;
    let mut sets: Vec<(usize, Vec<usize>)> = holding.into_iter().collect();
    sets.sort_unstable_by_key(|s| s.0);
    for (_, mut air) in sets {
        // The widest air cells first: a few, each cube a few hundred bounds.
        air.sort_by_key(|&i| (std::cmp::Reverse(leaves[i].cell.size), i));
        if ask.cover.tol == 0.0 || air.iter().take(8).any(|&i| ask.deep(ctx, grid.cell_box(leaves[i].cell))) {
            at_least += 1;
        } else {
            unshown.extend(air);
        }
    }
    let open: Vec<usize> = (0..n).filter(|&i| lands[i] == Lands::Open).collect();
    let mut thin: Vec<usize> = (0..n).filter(|&i| lands[i] == Lands::Thin).chain(thin_groups).collect();
    thin.sort_unstable();
    thin.dedup();
    let at_most = open.is_empty().then_some(groups.len());
    let decided = at_most == Some(at_least);
    let mut halve: Vec<usize> = apart.iter().chain(&open).copied().filter(|&i| leaves[i].cell.size > 1).collect();
    // Nothing else left to halve, the thin cells are: the air in them may
    // join a void to the outside, or two voids.
    let settles = !halve.is_empty() || decided || !open.is_empty();
    if !settles {
        // Those on the shortest way from each group that may be a void to
        // another group or the outside, where one leads: one way out decides
        // it, and the thin air elsewhere (along every outer face) does not.
        let owner: Vec<Option<usize>> = (0..n)
            .map(|i| match lands[i] {
                _ if edge(leaves[i].cell) && carries(i) => Some(outside),
                Lands::Group(a) => Some(find(&mut group, a)).filter(|&g| wide[g] || g == outside),
                _ => None,
            })
            .collect();
        let passes = |i: usize| owner[i].is_none() && carries(i);
        halve = ways(n, &faces, &owner, &passes, &groups).into_iter().filter(|&i| leaves[i].cell.size > 1).collect();
        if halve.is_empty() {
            halve = (0..n).filter(|&i| lands[i] == Lands::Thin && leaves[i].cell.size > 1).collect();
        }
    }
    halve.sort_unstable();
    halve.dedup();
    let void_at = (0..n).find(|&i| leaves[i].kind == Kind::Air && root_of[i].is_some_and(|g| wide[g])).map(|i| centre(grid, leaves[i].cell));
    let mut undecided = None;
    let mut stuck: Vec<usize> = Vec::new();
    if !decided {
        // The cells that keep the bounds apart: those with no group, two
        // of different groups a face joins; else thin air, and air not
        // shown deeper than the accuracy.
        stuck = open.iter().chain(&apart).copied().collect();
        if stuck.is_empty() {
            stuck = thin.iter().chain(&unshown).copied().collect();
        }
        let side = |i: usize| grid.h * leaves[i].cell.size as f64;
        undecided = stuck.iter().map(|&i| side(i)).fold(None, |m: Option<(f64, f64)>, s| Some(m.map_or((s, s), |(a, b)| (a.min(s), b.max(s)))));
    }
    let at = stuck.first().map(|&i| centre(grid, leaves[i].cell)).or(void_at);
    let thin_at = thin.first().map(|&i| centre(grid, leaves[i].cell));
    let found = Voids { at_least, at_most, at, cells: n, finest: grid.h, undecided, ended: None, accuracy: ask.cover.tol, thin: thin.len(), thin_at };
    (found, halve, settles)
}
