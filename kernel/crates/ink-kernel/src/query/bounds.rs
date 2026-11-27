//! Enclosures of material in a box, over all of space, and about the z axis.

use super::proof::worked;
use super::rotation;
use crate::octree::{Cell, Grid, try_root_for};
use crate::search;
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_interval_in};
use ink_core::interval::{I, add_o, mul_o};
use ink_core::tape::{BOTH, Tape};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, OnceLock};

pub fn max_extent(lo: V3, hi: V3) -> f64 {
    (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max)
}

/// The tape simplified for the region `lo..hi`, the field there and only
/// there: one interval pass settles which min/max branches can matter in
/// it, so a union of thousands of parts costs only the parts nearby in
/// everything that follows. Searches with it examine boxes inside the
/// region (their cells clipped to it: octree::try_root_for).
pub fn focus(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3) -> Arc<Tape> {
    let b: [I; 3] = std::array::from_fn(|k| I::new(lo[k], hi[k]));
    let mut ch = vec![BOTH; tape.len()];
    eval_interval(tape, ctx, b, Some(&mut ch));
    tape.narrowed(&ch)
}

/// Proven outer bounding box of the region where the field is < 0, found in
/// two passes (coarse, then refined). Returns `(lo, hi, cell_size)`; `None`
/// means the field is proven empty.
pub fn bounds(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3) -> Option<(V3, V3, f64)> {
    crate::cache::memo_region("bounds", &[tape], ctx, (lo, hi), |_| {}, |c| bounds_passes(tape, c, lo, hi, &[64.0, 256.0]))
}

/// Where a field can have material (be negative), over all of space. Its
/// zero set without interior (a surface of a field never negative) is not
/// material: a field 0 everywhere is empty.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Extent {
    /// Interval bounds prove the field positive everywhere.
    Empty,
    /// The field is proven positive outside this box.
    Bounded(V3, V3),
    /// Not proven bounded: material may reach arbitrarily far (a half-space,
    /// an unclipped lattice), or interval bounds cannot show that it does not.
    Unbounded,
}

/// Cells one search of `extent` may evaluate before it gives up and
/// reports `Unbounded`.
const EXTENT_BUDGET: usize = 1 << 13;

/// 2^40: farther than any model reaches, in window sizes (`reaches_far`).
const TB: f64 = 1_099_511_627_776.0;

/// ink's range: the farthest coordinate (mm), 1,000 km. Material beyond it
/// is outside what ink answers about.
pub const REACH: f64 = 1e9;

/// Where the field can be negative, over all of space, proven by covering all of space with pieces on which
/// interval bounds (`eval_interval_linear`, sound with infinite ends, whose
/// reals are finite) show the field positive, except pieces that join the
/// box. Two covers, as each settles shapes the other cannot: far sectors
/// (`extent_sectors`: narrow directions out to infinity, so a part rotated
/// about a slanted axis is settled), then, if that fails, boxes in
/// compactified coordinates (`extent_boxes`: slabs out to infinity, so a
/// twist's angle narrows with its axis). So the answer cannot miss material
/// however far it lies. Where both covers give up, a rotation's cylinder
/// can still prove the whole shape fits the window (rotation::inside_window).
/// A field no proof bounds, such as a half-space or an unclipped lattice,
/// is reported as such.
///
/// The covers search the field's members, not just the whole
/// (`extent_within`): a union's material is its members' material, and an
/// intersection's lies in any member's, or in the union of a union
/// member's members' boxes, each with the rest (an assembly's box: its
/// parts', each part's clipped by what the rest leaves). So a
/// union of parts that need different covers is bounded, and a search
/// runs only on small tapes. Every search starts with a small budget,
/// doubled up to `EXTENT_BUDGET` while undecided, so no member spends more
/// than its cheapest proof needs (an intersection whose members are
/// unbounded, as an extruded profile, is settled whole at once); the
/// next budget resumes each unfinished cover at the cell that stopped it.
/// A field with material where no model reaches (`reaches_far`: a half-space,
/// a hole's outside) is unbounded at once. Members are searched cheapest first
/// (`members`), so an intersection's box is its cheapest member's proven
/// bounded at the least budget.
pub fn extent(tape: &Arc<Tape>, ctx: &Ctx, window: (V3, V3)) -> Extent {
    crate::cache::memo_region("extent", &[tape], ctx, window, |_| {}, |c| extent_now(tape, c, window))
}

pub(super) fn extent_now(tape: &Arc<Tape>, ctx: &Ctx, window: (V3, V3)) -> Extent {
    let mut found = Found::default();
    let mut budget = EXTENT_BUDGET >> 10;
    loop {
        if let Some(e) = extent_within(tape, tape.len() - 1, ctx, window, budget, &mut found) {
            return e;
        }
        if budget >= EXTENT_BUDGET {
            return Extent::Unbounded;
        }
        budget *= 2;
    }
}

/// What `extent` found of each instruction (by its index) and each set of
/// instructions searched as their max (sorted, then usize::MAX): its
/// extent, or the budget it was undecided at.
#[derive(Default)]
pub(super) struct Found {
    results: HashMap<Vec<usize>, Result<Extent, usize>>,
    /// Unfinished covers, by the instructions combined in their evaluation order.
    searches: HashMap<Vec<usize>, Covers>,
}

/// The extent of instruction `i` of `tape` within `budget` evaluations per
/// search (None: undecided at that budget; Unbounded: material shown where
/// no model reaches).
fn extent_within(tape: &Arc<Tape>, i: usize, ctx: &Ctx, window: (V3, V3), budget: usize, found: &mut Found) -> Option<Extent> {
    use ink_core::graph::BinOp;
    use ink_core::tape::Op;
    let i = value_of(tape, i);
    match found.results.get([i].as_slice()) {
        Some(Ok(e)) => return Some(*e),
        Some(Err(b)) if *b >= budget => return None,
        _ => {}
    }
    let is = |i: usize, o: BinOp| matches!(tape.ssa[i].op, Op::Bin(x) if x == o);
    let e = match tape.ssa[i].op {
        // Material in any member: the union of their boxes.
        Op::Bin(BinOp::Min) => union(members(tape, i, BinOp::Min).into_iter(), |m| extent_within(tape, m, ctx, window, budget, found)),
        // Material in every member: inside any member's box. Else, with a
        // union among them (the costliest) some of whose members are
        // bounded alone (an assembly's parts), in the union of its
        // members' boxes, each member not bounded alone found whole with
        // the rest (a part far out, with the outside of a box asked),
        // where those searches together read no more instructions than
        // the whole: never the whole union at once. Else found whole. One
        // member (others its tags or itself again) is the intersection.
        Op::Bin(BinOp::Max) => match members(tape, i, BinOp::Max)[..] {
            [m] => extent_within(tape, m, ctx, window, budget, found),
            ref ms => {
                let mut e = None;
                for &m in ms {
                    e = extent_within(tape, m, ctx, window, budget, found).filter(|e| *e != Extent::Unbounded);
                    if e.is_some() {
                        break;
                    }
                }
                let bounded = |e: &Option<Extent>| matches!(e, Some(Extent::Bounded(..) | Extent::Empty));
                // The union's members, each bounded alone or with the rest; None
                // as soon as those with the rest read more than the whole.
                let split = |u: usize, found: &mut Found| {
                    let rest: Vec<usize> = ms.iter().copied().filter(|&m| m != u).collect();
                    let mut cones = Cones::new(tape);
                    let (all, with, mut read, mut each) = (cones.len(&[i]), cones.len(&rest) + 1, 0, Vec::new());
                    for m in members(tape, u, BinOp::Min) {
                        let (s, e): (Vec<usize>, _) = (rest.iter().copied().chain([m]).collect(), extent_within(tape, m, ctx, window, budget, found));
                        if !bounded(&e) {
                            // At most the rest's and m's, apart.
                            read += with + cones.len(&[m]);
                            if read > all {
                                return None;
                            }
                        }
                        each.push((s, e));
                    }
                    each.iter().any(|(_, e)| bounded(e)).then_some(each)
                };
                match e {
                    Some(e) => Some(e),
                    None => match ms.iter().rev().find(|&&m| is(m, BinOp::Min)).and_then(|&u| split(u, found)) {
                        Some(each) => union(each.into_iter(), |(s, alone)| match alone {
                            Some(e @ (Extent::Bounded(..) | Extent::Empty)) => Some(e),
                            _ => found_whole(tape, &s, ctx, window, budget, found),
                        }),
                        None => found_whole(tape, &[i], ctx, window, budget, found),
                    },
                }
            }
        },
        _ => found.search(tape, &[i], ctx, window, budget),
    };
    found.results.insert(vec![i], e.ok_or(budget));
    e
}

/// The extent of the max of instructions `set` of `tape`, found whole
/// (slabs and profiles are unbounded alone); or none anywhere, by the
/// first-order bound of the max over all of space, which the covers' plain
/// bounds cannot see: members that never meet, as two parallel slabs, whose
/// exactly opposite slopes cancel however far they reach. That bound does
/// not depend on the budget: it is tried once, where the covers fail at the
/// first (smallest) budget, so an extrusion, which they settle there, never
/// pays for it. Kept in `found` as `extent_within` keeps its sets'.
fn found_whole(tape: &Arc<Tape>, set: &[usize], ctx: &Ctx, window: (V3, V3), budget: usize, found: &mut Found) -> Option<Extent> {
    let mut key = set.to_vec();
    key.sort_unstable();
    key.push(usize::MAX);
    match found.results.get(&key) {
        Some(Ok(e)) => return Some(*e),
        Some(Err(b)) if *b >= budget => return None,
        _ => {}
    }
    let first = !found.results.contains_key(&key);
    let never = || crate::centred::centred_at_least(&[(&whole(tape, set), 0.0)], ctx, [I::ENTIRE; 3], 0.0);
    let e = found.search(tape, set, ctx, window, budget).or_else(|| (first && never()).then_some(Extent::Empty));
    found.results.insert(key, e.ok_or(budget));
    e
}

/// The union of the extents `extent` finds of each of `items`: None at the
/// first undecided, unbounded at the first unbounded.
fn union<T>(items: impl Iterator<Item = T>, mut extent: impl FnMut(T) -> Option<Extent>) -> Option<Extent> {
    let mut u = Extent::Empty;
    for s in items {
        u = match (u, extent(s)?) {
            (_, Extent::Unbounded) => return Some(Extent::Unbounded),
            (Extent::Bounded(l0, h0), Extent::Bounded(l1, h1)) => {
                let (lo, hi) = ink_core::v3::box_union((l0, h0), (l1, h1));
                Extent::Bounded(lo, hi)
            }
            (Extent::Empty, e) => e,
            (e, _) => e,
        };
    }
    Some(u)
}

/// The instructions cones of a tape read, each walk marking those it
/// reaches with its own number (no clearing between walks).
struct Cones<'a> {
    tape: &'a Tape,
    mark: Vec<u32>,
    walks: u32,
}

impl<'a> Cones<'a> {
    fn new(tape: &'a Tape) -> Self {
        Cones { tape, mark: vec![0; tape.len()], walks: 0 }
    }

    /// The instructions the max of `set` reads (its cones, and the maxes
    /// joining them).
    fn len(&mut self, set: &[usize]) -> usize {
        self.walks += 1;
        let (mut stack, mut n) = (set.to_vec(), set.len() - 1);
        while let Some(j) = stack.pop() {
            if std::mem::replace(&mut self.mark[j], self.walks) != self.walks {
                n += 1;
                let (args, k) = self.tape.ssa[j].args();
                stack.extend(args[..k].iter().map(|&a| a as usize));
            }
        }
        n
    }
}

/// The tape computing the max of instructions `set` of `tape`: the tape
/// itself where that is its output.
fn whole(tape: &Arc<Tape>, set: &[usize]) -> Arc<Tape> {
    if set == [tape.len() - 1] { tape.clone() } else { Arc::new(tape.cone_max(set)) }
}

/// The instruction whose value instruction `i` passes on (through tags,
/// Lipschitz promises and derivatives' marks).
fn value_of(tape: &Tape, mut i: usize) -> usize {
    use ink_core::tape::Op;
    while let Op::Tag(_) | Op::Lipschitz | Op::Deriv(..) | Op::Real = tape.ssa[i].op {
        i = tape.ssa[i].a as usize;
    }
    i
}

/// The members of the `op` (min or max) at instruction `i`: its arguments,
/// with nested `op`s flattened, cheapest first: by the instructions each
/// reads (its cone: what evaluating it costs), then in tape order, which
/// is the hash order of commutative arguments (Graph::bin) and put a
/// pair's long belt before its small part. Each instruction is taken once
/// (min and max are idempotent), so a shared subgraph does not unfold.
fn members(tape: &Tape, i: usize, op: ink_core::graph::BinOp) -> Vec<usize> {
    use ink_core::tape::Op;
    // Each instruction marked with the last walk that reached it: the
    // flattening's (i), then each member's own (m).
    let (mut out, mut stack, mut mark) = (Vec::new(), vec![i], vec![usize::MAX; tape.ssa.len()]);
    while let Some(j) = stack.pop() {
        let j = value_of(tape, j);
        if std::mem::replace(&mut mark[j], i) == i {
            continue;
        }
        match tape.ssa[j].op {
            Op::Bin(o) if o == op => stack.extend([tape.ssa[j].b as usize, tape.ssa[j].a as usize]),
            _ => out.push(j),
        }
    }
    out.sort_by_cached_key(|&m| {
        let (mut n, mut stack) = (0, vec![m]);
        mark[m] = m;
        while let Some(j) = stack.pop() {
            n += 1;
            let (args, k) = tape.ssa[j].args();
            for &a in &args[..k] {
                if std::mem::replace(&mut mark[a as usize], m) != m {
                    stack.push(a as usize);
                }
            }
        }
        (n, m)
    });
    out
}

/// Each cover walks its cells depth first. On a larger budget, resume at
/// the cell that exhausted the last one, with the same count, traversal
/// order and box of material found so far. The two budgets are independent.
struct ExtentWalk<C> {
    stack: Vec<Spec<C>>,
    lo: V3,
    hi: V3,
    evals: usize,
}

/// A cell of a cover's walk, and what stepping it found where that was
/// found ahead of the walk (`ahead`).
struct Spec<C> {
    cell: C,
    step: Option<Box<Step<C>>>,
}

/// What a cover's cell comes to (a function of the cell alone, so it can be
/// found ahead of the walk, on other threads): the length of each tape it
/// evaluated, for `worked`.
enum Step<C> {
    /// Not evaluated: no real point is in it.
    Skip,
    /// Evaluated, and settled: no material, or no real point.
    Settled(usize),
    /// Evaluated, and joins the box found.
    Join(usize, [I; 3]),
    /// Evaluated (and a second tape, once more, where a far sector checks
    /// its far end), and split into `kids`, pushed in order.
    Split(usize, Option<usize>, Vec<Spec<C>>),
}

/// How deep a walk steps a cell's kids ahead (`ExtentWalk::walk`). Two,
/// three and four took the same time on Limestone's flywheel belt claim at
/// four threads (its covers' far sectors each wait on the one before).
const AHEAD_DEPTH: u32 = 2;

impl<C: Send> ExtentWalk<C> {
    fn new(stack: Vec<C>) -> Self {
        Self { stack: stack.into_iter().map(|cell| Spec { cell, step: None }).collect(), lo: [f64::INFINITY; 3], hi: [f64::NEG_INFINITY; 3], evals: 0 }
    }

    /// The walk within `budget` evaluations: each cell popped, its step
    /// taken in this order, as one thread would take them. Steps are found
    /// ahead (`ahead`) for the cell popped and the next cells on the
    /// stack without one, each with its kids a few deep, on the pool: each
    /// is a function of its cell, so the walk, its count and where a budget
    /// stops it are the same on any number of threads. A cell found ahead
    /// and not yet walked keeps its step for the next budget; what the
    /// last budget leaves (a few cells per thread) is the only work more.
    fn walk(&mut self, budget: usize, step: &(dyn Fn(&C) -> Step<C> + Sync)) -> Option<Extent> {
        let threads = rayon::current_num_threads();
        while let Some(mut x) = self.stack.pop() {
            if x.step.is_none() {
                if threads > 1 {
                    let mut picks: Vec<&mut Spec<C>> =
                        std::iter::once(&mut x).chain(self.stack.iter_mut().rev().take(4 * threads).filter(|s| s.step.is_none()).take(threads - 1)).collect();
                    picks.par_iter_mut().for_each(|s| ahead(s, AHEAD_DEPTH, step));
                } else {
                    x.step = Some(Box::new(step(&x.cell)));
                }
            }
            match *x.step.take().expect("a step found") {
                Step::Skip => {}
                Step::Settled(n) => {
                    self.evals += 1;
                    super::proof::worked_len(n);
                }
                Step::Join(n, b) => {
                    self.evals += 1;
                    super::proof::worked_len(n);
                    self.join(b);
                }
                Step::Split(n, far, kids) => {
                    self.evals += 1;
                    super::proof::worked_len(n);
                    if self.evals >= budget {
                        x.step = Some(Box::new(Step::Split(n, far, kids)));
                        self.retry(x);
                        return None;
                    }
                    if let Some(m) = far {
                        self.evals += 1;
                        super::proof::worked_len(m);
                    }
                    self.stack.extend(kids);
                }
            }
        }
        Some(self.extent())
    }

    fn join(&mut self, b: [I; 3]) {
        for k in 0..3 {
            self.lo[k] = self.lo[k].min(b[k].lo);
            self.hi[k] = self.hi[k].max(b[k].hi);
        }
    }

    fn extent(&self) -> Extent {
        if self.lo[0] > self.hi[0] { Extent::Empty } else { Extent::Bounded(self.lo, self.hi) }
    }

    fn retry(&mut self, cell: Spec<C>) {
        self.stack.push(cell);
        self.evals -= 1;
    }
}

/// `x`'s step, and its kids' to `depth` more, side by side on the pool.
fn ahead<C: Send>(x: &mut Spec<C>, depth: u32, step: &(dyn Fn(&C) -> Step<C> + Sync)) {
    if x.step.is_none() {
        x.step = Some(Box::new(step(&x.cell)));
    }
    if depth > 0
        && let Some(Step::Split(_, _, kids)) = x.step.as_deref_mut()
    {
        kids.par_iter_mut().for_each(|k| ahead(k, depth - 1, step));
    }
}

struct Covers {
    sectors: ExtentWalk<Sector>,
    boxes: ExtentWalk<CompactBox>,
}

impl Found {
    pub(super) fn search(&mut self, tape: &Arc<Tape>, set: &[usize], ctx: &Ctx, window: (V3, V3), budget: usize) -> Option<Extent> {
        let covers = match self.searches.get_mut(set) {
            Some(covers) => covers,
            None => {
                let t = whole(tape, set);
                if reaches_far(&t, ctx, window) {
                    return Some(Extent::Unbounded);
                }
                self.searches.entry(set.to_vec()).or_insert_with(|| Covers::new(&t, ctx, window))
            }
        };
        let e = extent_sectors(&mut covers.sectors, ctx, window, budget)
            .or_else(|| extent_boxes(&mut covers.boxes, ctx, window, budget))
            .or_else(|| (budget == EXTENT_BUDGET && rotation::inside_window(&whole(tape, set), ctx, window)).then_some(Extent::Bounded(window.0, window.1)));
        if e.is_some() {
            self.searches.remove(set);
        }
        e
    }
}

/// Whether the field is shown to have material where no model reaches: at
/// a point TB windows out from the window's centre, along one of the 26
/// directions to a cube's faces, edges and corners (a half-space holds
/// one, a hole's outside all of them), bounded below 0 there over the
/// scope. No cover bounds such a field, so it is not searched.
fn reaches_far(tape: &Tape, ctx: &Ctx, window: (V3, V3)) -> bool {
    let (c, r) = centre_radius(window);
    let mut slots = Vec::new();
    (0..27).filter(|&k| k != 13).any(|k| {
        let d = [k % 3, k / 3 % 3, k / 9].map(|x| x as f64 - 1.0);
        let p: [I; 3] = std::array::from_fn(|a| I::point(c[a] + TB * r[a] * d[a]));
        worked(tape);
        eval_interval_in(tape, ctx, p, None, &mut slots).all_inside()
    })
}

/// A window's centre c and half-sides r, as the covers take them: their
/// pieces are c + r (s q), whatever window they came from, the central box
/// the one they start from (s = 1), so together they cover all of space.
/// (r must not be 0, or the sectors would collapse onto c; both are halved
/// before the sum, which overflows for a window's ends near 1e308.)
fn centre_radius(window: (V3, V3)) -> (V3, V3) {
    (std::array::from_fn(|k| 0.5 * window.0[k] + 0.5 * window.1[k]), std::array::from_fn(|k| (0.5 * window.1[k] - 0.5 * window.0[k]).max(f64::MIN_POSITIVE)))
}

/// A sector: the face (axis, side), its direction ranges along the
/// other two axes (in [-1, 1]), and its range of distance (scale of the
/// window, from 1; hi infinite for the far end).
#[derive(Clone)]
struct Sector {
    axis: usize,
    side: f64,
    v: [(f64, f64); 2],
    s: (f64, f64),
    tape: Arc<Tape>,
}

struct CompactBox {
    lo: V3,
    hi: V3,
    tape: Arc<Tape>,
}

impl Covers {
    fn new(tape: &Arc<Tape>, ctx: &Ctx, window: (V3, V3)) -> Self {
        let (c, r) = centre_radius(window);
        let mut sectors = ExtentWalk::new(
            (0..6).map(|f| Sector { axis: f / 2, side: if f % 2 == 0 { -1.0 } else { 1.0 }, v: [(-1.0, 1.0); 2], s: (1.0, f64::INFINITY), tape: tape.clone() }).collect(),
        );
        // The window joins as it is, rounded outward to contain c +- r.
        let wbox = std::array::from_fn(|k| I::new((c[k] - r[k]).next_down(), (c[k] + r[k]).next_up()));
        if ink_core::eval::eval_interval_linear(tape, ctx, &ink_core::eval::Region::Box(wbox), None).may_be_inside() {
            sectors.join(wbox);
        }
        let boxes = ExtentWalk::new(vec![CompactBox { lo: [-1.0; 3], hi: [1.0; 3], tape: tape.clone() }]);
        Self { sectors, boxes }
    }
}

/// `extent` by far sectors: around the window, six cones (one per face,
/// scaled out from its centre), each split into sectors: a range of
/// directions and a range of distance, the farthest reaching infinity. Over
/// a sector a linear term's distance factor is taken once, so a rotated
/// coordinate stays one-signed on a narrow sector however far it reaches,
/// and a radius and an angle (a polar fold's) stay narrow.
/// Undecided sectors split by direction, then by distance (doubling); a
/// finite one small for its distance joins the box, a far one never does.
fn extent_sectors(cover: &mut ExtentWalk<Sector>, ctx: &Ctx, window: (V3, V3), budget: usize) -> Option<Extent> {
    use ink_core::eval::{Region, eval_interval_linear};
    let (c, r) = centre_radius(window);
    let others = |axis: usize| -> [usize; 2] { [(axis + 1) % 3, (axis + 2) % 3] };
    let region = |x: &Sector| {
        let mut q = [I::point(0.0); 3];
        q[x.axis] = I::point(x.side);
        for (j, &k) in others(x.axis).iter().enumerate() {
            q[k] = I::new(x.v[j].0, x.v[j].1);
        }
        Region::Sector { c, r, q, s: I::new(x.s.0, x.s.1) }
    };
    let step = |x: &Sector| -> Step<Sector> {
        // Doubled past the largest float: no real point is that far.
        if x.s.0.is_infinite() {
            return Step::Skip;
        }
        let reg = region(x);
        let mut ch = vec![BOTH; x.tape.len()];
        let v = eval_interval_linear(&x.tape, ctx, &reg, Some(&mut ch));
        let n = x.tape.len();
        if v.none_inside() {
            return Step::Settled(n);
        }
        let far = x.s.1.is_infinite();
        let width = (x.v[0].1 - x.v[0].0).max(x.v[1].1 - x.v[1].0);
        // A finite sector small for its distance (directions an eighth of a
        // face, or half of one twice the window's size out: the box is
        // tightened later, and a big shape far out stays cheap) is material,
        // or too fine to settle: it joins the box.
        let fine = if x.s.0 >= 2.0 { 0.5 } else { 0.125 };
        if !far && width <= fine && x.s.1 <= 2.0 * x.s.0 {
            return Step::Join(n, reg.bounds());
        }
        let sub = x.tape.narrowed(&ch);
        let kid = |f: &dyn Fn(&mut Sector)| {
            let mut y = Sector { tape: sub.clone(), ..x.clone() };
            f(&mut y);
            Spec { cell: y, step: None }
        };
        // A finite sector splits its directions. A far one splits them only
        // if they do not settle its far end beyond ink's range (REACH: there
        // only narrower directions can help); else its distance: the far end
        // doubles away (material in between stays in the finite pieces). So
        // a distance walked out is settled within the range, and a sector
        // settled only farther out (an angle about an axis off the origin,
        // whose offset fades only there) is narrowed instead.
        let checked = far && width > 1.0 / 64.0;
        let narrower = !far
            || checked && {
                let mut y = region(x);
                if let Region::Sector { s, .. } = &mut y {
                    *s = I::new(x.s.0.max((REACH + c[x.axis].abs()) / r[x.axis]), f64::INFINITY);
                }
                eval_interval_linear(&sub, ctx, &y, None).may_be_inside()
            };
        let kids = if narrower {
            let j = if x.v[0].1 - x.v[0].0 >= x.v[1].1 - x.v[1].0 { 0 } else { 1 };
            let m = 0.5 * (x.v[j].0 + x.v[j].1);
            vec![kid(&|y| y.v[j].1 = m), kid(&|y| y.v[j].0 = m)]
        } else {
            let m = 2.0 * x.s.0;
            vec![kid(&|y| y.s.1 = m), kid(&|y| y.s.0 = m)]
        };
        Step::Split(n, checked.then(|| sub.len()), kids)
    };
    cover.walk(budget, &step)
}

/// The compactified coordinate of `extent_boxes` on one axis: u in (-1, 1)
/// maps to x = c + r u / (1 - |u|), so u = +-1/2 is the window's faces c +- r,
/// and u = +-1 is infinity.
fn from_u(u: f64, c: f64, r: f64) -> f64 {
    if u >= 1.0 {
        f64::INFINITY
    } else if u <= -1.0 {
        f64::NEG_INFINITY
    } else if u == 0.5 {
        c + r
    } else if u == -0.5 {
        c - r
    } else {
        c + r * u / (1.0 - u.abs())
    }
}

/// `extent` by boxes: an octree over compactified coordinates centred on
/// the window (its root is all of space, the cell at the middle the window,
/// cells at the rim reach infinity). Only cells reaching infinity are
/// refined, widest axis first (so a slab's finite axes narrow too); finite
/// cells not proven empty join the box.
fn extent_boxes(cover: &mut ExtentWalk<CompactBox>, ctx: &Ctx, window: (V3, V3), budget: usize) -> Option<Extent> {
    use ink_core::eval::{Region, eval_interval_linear};
    let (c, r) = centre_radius(window);
    let step = |cell: &CompactBox| -> Step<CompactBox> {
        let b: [I; 3] = std::array::from_fn(|k| I::new(from_u(cell.lo[k], c[k], r[k]), from_u(cell.hi[k], c[k], r[k])));
        let mut ch = vec![BOTH; cell.tape.len()];
        let v = eval_interval_linear(&cell.tape, ctx, &Region::Box(b), Some(&mut ch));
        let n = cell.tape.len();
        // None inside; or a cell at u = +-1 alone, which holds no real point.
        if v.none_inside() || (0..3).any(|k| cell.lo[k] >= 1.0 || cell.hi[k] <= -1.0) {
            return Step::Settled(n);
        }
        let infinite: Vec<usize> = (0..3).filter(|&k| cell.lo[k] <= -1.0 || cell.hi[k] >= 1.0).collect();
        if infinite.is_empty() {
            return Step::Join(n, b);
        }
        let sub = cell.tape.narrowed(&ch);
        // Halve the widest axis (one that reaches infinity on a tie): finite
        // axes split too, so a twist's angle over a slab narrows.
        let width = |k: usize| (cell.hi[k] - cell.lo[k], infinite.contains(&k));
        let k = (0..3).max_by(|&i, &j| width(i).partial_cmp(&width(j)).unwrap()).unwrap();
        let m = 0.5 * (cell.lo[k] + cell.hi[k]);
        let (mut h1, mut l2) = (cell.hi, cell.lo);
        (h1[k], l2[k]) = (m, m);
        Step::Split(
            n,
            None,
            vec![Spec { cell: CompactBox { lo: cell.lo, hi: h1, tape: sub.clone() }, step: None }, Spec { cell: CompactBox { lo: l2, hi: cell.hi, tape: sub }, step: None }],
        )
    };
    cover.walk(budget, &step)
}

/// Two quick passes (1/32, then 1/64 of the narrowed region): enough to focus a proof's search.
pub fn bounds_quick(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3) -> Option<(V3, V3, f64)> {
    crate::cache::memo_region("bounds-quick", &[tape], ctx, (lo, hi), |_| {}, |c| bounds_passes(tape, c, lo, hi, &[32.0, 64.0]))
}

/// How far the material of a shape (where the field is < 0) inside `lo..hi`
/// lies from the z axis, at every value of ctx's scope (t a step, ranged
/// params): a lower bound (a point some material reaches: the farthest
/// corner of a box wholly material, or a box's middle proven material) and
/// an upper, within `tol` of each other where `budget` boxes allow. The one
/// search (crate::search), farthest first (a box's bound is its farthest
/// corner from the axis, negated: none of its points is farther; the upper
/// bound is the farthest of the boxes left and those kept), halving
/// across x and y, and across z where `split_z` (a field that reads z) and
/// z is the box's longest side; a box proven empty is settled, one wholly
/// material, within `tol` across (from the axis or the lower bound) or no
/// longer halved keeps its farthest corner in the upper bound. None:
/// proven empty.
pub fn radial_extent(tape: &Tape, ctx: &Ctx, lo: V3, hi: V3, split_z: bool, tol: f64, budget: usize) -> Option<(f64, f64)> {
    use std::cell::Cell;
    let far = |b: &[I; 3]| radius(b[0].lo.abs().max(b[0].hi.abs()), b[1].lo.abs().max(b[1].hi.abs())).hi;
    // The farthest corner itself, rounded down: a point of the box.
    let corner = |b: &[I; 3]| {
        let pick = |x: I| if x.lo.abs() > x.hi.abs() { x.lo } else { x.hi };
        radius(pick(b[0]), pick(b[1])).lo
    };
    // `kept`: the farthest of the boxes that may hold material and are
    // halved no further; `seen`: whether any may.
    let (lower, kept, seen) = (Cell::new(0.0f64), Cell::new(0.0f64), Cell::new(false));
    let test = |n: &search::Node<f64>| -> search::Test<f64, ()> {
        let v = eval_interval(tape, &n.ctx, n.b, None);
        if v.none_inside() {
            return search::Test::Settled;
        }
        seen.set(true);
        let gap = |x: I| {
            if x.lo > 0.0 {
                x.lo
            } else if x.hi < 0.0 {
                -x.hi
            } else {
                0.0
            }
        };
        if v.all_inside() || far(&n.b) - radius(gap(n.b[0]), gap(n.b[1])).lo.max(lower.get()) <= tol {
            if v.all_inside() {
                lower.set(lower.get().max(corner(&n.b)));
            }
            kept.set(kept.get().max(far(&n.b)));
            return search::Test::Settled;
        }
        let mid = n.b.map(|x| x.mid());
        if eval_interval(tape, &n.ctx, mid.map(I::point), None).all_inside() {
            lower.set(lower.get().max(radius(mid[0], mid[1]).lo));
        }
        let w = n.b.map(|x| x.hi - x.lo);
        let space = [w[0], w[1], if split_z && w[2] > w[0].max(w[1]) { w[2] } else { 0.0 }];
        search::Test::Open { bound: -far(&n.b), smear: search::Smear { space, ..Default::default() }, state: far(&n.b) }
    };
    let root = search::Node { b: std::array::from_fn(|k| I::new(lo[k], hi[k])), ctx: Arc::new(ctx.clone()), state: 0.0 };
    let first = test(&root);
    let mut s = search::Search::start(vec![(root, first)], search::Precision { finest: 0.0, splits: 0, maps: true, closing: false });
    s.serially(budget, 8, test);
    let upper = s.left().map(|n| n.state).fold(kept.get(), f64::max);
    (seen.get() || s.left().next().is_some()).then(|| (lower.get(), upper.max(lower.get())))
}

/// sqrt(x^2 + y^2) bounded: lo <= the real radius <= hi. Scaled by a
/// power of two near the larger term so the squares neither underflow nor
/// overflow (a radius of 1e-200 is not 0), and each operation (scale,
/// square, sum, sqrt, rescale) rounded outward: one ulp after a float norm
/// bounds nothing.
fn radius(x: f64, y: f64) -> I {
    use ink_core::graph::UnOp;
    use ink_core::interval::un;
    let m = x.abs().max(y.abs());
    if m == 0.0 || !m.is_finite() {
        return I::point(m);
    }
    let e = (m.log2().floor() as i32).clamp(-1000, 1000);
    let sq = |v: f64| un(UnOp::Square, mul_o(I::point(v), I::point(2f64.powi(-e))));
    mul_o(un(UnOp::Sqrt, add_o(sq(x), sq(y))), I::point(2f64.powi(e)))
}

/// A box narrowed to under 3 cells of the last pass (a small part far
/// inside a big region) gets more passes, while they narrow it: each such
/// pass shrinks it at least 85 times, and a box of floats cannot shrink for
/// ever (a point-like set reaches the precision of floats in a few).
fn bounds_passes(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, passes: &[f64]) -> Option<(V3, V3, f64)> {
    let (mut lo, mut hi) = (lo, hi);
    let mut h = 0.0;
    let last = passes[passes.len() - 1];
    let extra = passes.iter().copied().chain(std::iter::repeat(last));
    let mut n = 0;
    let mut before = None;
    let parts = material_parts(tape);
    for divisions in extra {
        n += 1;
        if n > passes.len() && (max_extent(lo, hi) >= 3.0 * h || before == Some((lo, hi))) {
            break;
        }
        before = Some((lo, hi));
        h = max_extent(lo, hi) / divisions;
        // No root covers the region: no bounds tighter than it.
        let Some((grid, root)) = try_root_for(lo, hi, h, 1) else { return Some((lo, hi, h)) };
        let ((nlo, nhi), ..) = bounds_cells(&parts, ctx, &grid, root);
        if nlo[0] > nhi[0] {
            return None;
        }
        for k in 0..3 {
            lo[k] = lo[k].max(nlo[k]);
            hi[k] = hi[k].min(nhi[k]);
        }
    }
    Some((lo, hi, h))
}

/// Cells already enclosed by the box found so far cannot enlarge it.
/// For a cell crossing that box, evaluate only its outside slabs first:
/// its inside material must not force subdivision of an entire face just
/// to exclude the empty strip beyond the current bound. The result is
/// an enclosure of the material, as in a full descent. Each of the
/// field's parts (`material_parts`: an assembly's) descends on its own
/// tape, all into the one box, so a part inside the box found so far
/// costs its own few slabs, and no cell reads the whole.
/// With the cells it evaluated and the longest tape it read: work a test
/// can bound, whatever the timing.
pub(super) fn bounds_cells(parts: &[Arc<Tape>], ctx: &Ctx, grid: &Grid, root: Cell) -> ((V3, V3), usize, usize) {
    let longest = parts.iter().map(|t| t.len()).max().unwrap_or(0);
    let (visited, waiting) = (AtomicUsize::new(0), AtomicUsize::new(0));
    let threads = if longest >= PARALLEL_BOUNDS { rayon::current_num_threads() } else { 1 };
    let found = if threads > 1 {
        rayon::in_place_scope(|scope| bounds_below(parts, ctx, grid, root, &visited, Some((scope, &waiting, 8 * threads))))
    } else {
        bounds_below(parts, ctx, grid, root, &visited, None)
    };
    (found, visited.into_inner(), longest)
}

const PARALLEL_BOUNDS: usize = 4096;
type CellRange = (I, Vec<u8>);

fn bounds_below<'scope>(
    parts: &[Arc<Tape>],
    ctx: &'scope Ctx,
    grid: &Grid,
    root: Cell,
    visited: &'scope AtomicUsize,
    scope: Option<(&rayon::Scope<'scope>, &'scope AtomicUsize, usize)>,
) -> (V3, V3) {
    let (mut lo, mut hi) = ink_core::v3::EMPTY_BOX;
    let mut pending: Vec<_> = parts.iter().map(|t| (t.clone(), root, None::<Arc<OnceLock<CellRange>>>)).collect();
    let parallel = |t: &Tape| t.len() >= PARALLEL_BOUNDS && scope.is_some();
    let full = move |t: &Tape, bx| {
        visited.fetch_add(1, Relaxed);
        let mut choices = vec![BOTH; t.len()];
        (eval_interval(t, ctx, bx, Some(&mut choices)), choices)
    };
    let outside = |bx: [I; 3], lo: V3, hi: V3| (0..3).any(|k| bx[k].lo < lo[k] || bx[k].hi > hi[k]);
    while let Some((t, c, prefetched)) = pending.pop() {
        let local = OnceLock::new();
        let range = prefetched.as_deref().unwrap_or(&local);
        if !grid.meets(c) {
            continue;
        }
        let bx = grid.cell_box(c);
        if !outside(bx, lo, hi) || range.get().is_some_and(|(v, _)| v.none_inside()) {
            continue;
        }
        if lo[0] <= hi[0] {
            let mut slabs = [bx; 6];
            let mut n = 0;
            let mut whole = false;
            for k in 0..3 {
                for side in
                    [(bx[k].lo < lo[k]).then(|| I::new(bx[k].lo, bx[k].hi.min(lo[k]))), (bx[k].hi > hi[k]).then(|| I::new(bx[k].lo.max(hi[k]), bx[k].hi))].into_iter().flatten()
                {
                    whole |= side.lo.to_bits() == bx[k].lo.to_bits() && side.hi.to_bits() == bx[k].hi.to_bits();
                    slabs[n][k] = side;
                    n += 1;
                }
            }
            // A slab equal to the cell answers both questions: if it may
            // hold material we must inspect the cell, and if empty there
            // is no material in the whole cell.
            let read = |&b: &[I; 3]| {
                visited.fetch_add(1, Relaxed);
                eval_interval(&t, ctx, b, None).may_be_inside()
            };
            let inside = if whole {
                crate::once::get_or_compute(range, || full(&t, bx)).0.may_be_inside()
            } else if n > 1 && parallel(&t) {
                slabs[..n].par_iter().any(read)
            } else {
                slabs[..n].iter().any(read)
            };
            if !inside {
                continue;
            }
        }
        let (r, ch) = crate::once::get_or_compute(range, || full(&t, bx));
        if r.none_inside() {
            continue;
        }
        if r.all_inside() || c.size == 1 {
            (lo, hi) = ink_core::v3::box_union((lo, hi), grid.cell_bounds(c));
        } else {
            let sub = t.narrowed(ch);
            pending.extend(c.children().into_iter().map(|child| {
                // Only the full range is independent of the current box.
                // Start those evaluations while descending the next child
                // in order, then settle its slabs against the box found so far.
                let b = grid.cell_box(child);
                let prepared =
                    scope.filter(|(_, waiting, limit)| parallel(&sub) && waiting.load(Relaxed) < *limit && grid.meets(child) && outside(b, lo, hi)).map(|(scope, waiting, _)| {
                        waiting.fetch_add(1, Relaxed);
                        let slot = Arc::new(OnceLock::new());
                        let (out, t) = (slot.clone(), sub.clone());
                        // The initializer only runs the interval evaluator, so
                        // the descent can also take an unstarted evaluation.
                        scope.spawn(move |_| {
                            crate::once::get_or_compute(&out, || full(&t, b));
                            waiting.fetch_sub(1, Relaxed);
                        });
                        slot
                    });
                (sub.clone(), child, prepared)
            }));
        }
    }
    (lo, hi)
}

/// A field as pieces whose materials' union is its material, each its own
/// tape: a union's members; an intersection's, where a union is among its
/// members, that union's members each with the rest (an assembly clipped
/// to a box: each part clipped to it); else the field whole. Split only
/// where the pieces together are no longer than twice the whole (what
/// they share is read once a piece).
pub(super) fn material_parts(tape: &Arc<Tape>) -> Vec<Arc<Tape>> {
    use ink_core::graph::BinOp;
    use ink_core::tape::Op;
    let is = |i: usize, o: BinOp| matches!(tape.ssa[i].op, Op::Bin(x) if x == o);
    let i = value_of(tape, tape.len() - 1);
    let sets: Vec<Vec<usize>> = if is(i, BinOp::Min) {
        members(tape, i, BinOp::Min).into_iter().map(|m| vec![m]).collect()
    } else if is(i, BinOp::Max) {
        let ms = members(tape, i, BinOp::Max);
        let Some(&u) = ms.iter().rev().find(|&&m| is(m, BinOp::Min)) else { return vec![tape.clone()] };
        let rest: Vec<usize> = ms.iter().copied().filter(|&m| m != u).collect();
        members(tape, u, BinOp::Min).into_iter().map(|m| rest.iter().copied().chain([m]).collect()).collect()
    } else {
        return vec![tape.clone()];
    };
    let mut cones = Cones::new(tape);
    let read: usize = sets.iter().map(|s| cones.len(s)).sum();
    if sets.len() < 2 || read > 2 * tape.len() {
        return vec![tape.clone()];
    }
    sets.iter().map(|s| whole(tape, s)).collect()
}
