//! The one search: a deterministic best-first branch and bound over boxes
//! of space, of the ranged params and of time. A kernel question is a test
//! of a box (`Question`); the search decides where to look, how to split
//! and when to stop (docs/IR.md, Resolution). Emptiness (query::prove_empty)
//! and clearance (query::clearance) ask through it over space and the
//! ranged params; a value over the ranged params and time
//! (ranged::value_span) over the scope beyond space, all of space one box
//! tested whole.
//!
//! A test finds a box settled (nothing in it can violate the claim), a
//! witness, or open: undecided, with a sound bound of what it proves there
//! (lower is nearer a violation: a clearance's is how far every point of
//! one shape in the box is from the other; an emptiness's box proves none
//! short of settling it, so its boxes tie; a motion's step, the margin it
//! still needs, negated) and each axis's smear, how far that axis spreads
//! the test's bounds over the box (its rate there times its half-width:
//! Kearfott's smear). Open boxes wait lowest bound first: the bound the
//! search has attained is the least of theirs, and only refining that box
//! can raise it (best first, as interval branch and bound does: Hansen,
//! Kearfott); ties go to the smaller box, nearer being settled, then to the
//! first made. A box is halved on its axis of largest smear, with every
//! space axis the next halving would choose anyway, and its halves start
//! from its state (its tapes simplified over it). A box is never halved
//! on a space axis along which its question is the same all through it
//! (`Smear::still`: an extrusion's span, where the walls read no z). A box
//! no axis halves (every side no wider than `finest`, the accuracy asked,
//! still, or reaching infinity, tested whole; every param halved `splits`
//! times, or spreading them less than space does; t at adjacent floats; no
//! axis spreading them)
//! is stuck: the claim is not decided there at that accuracy. Asked
//! whether a claim holds, the proof then ends unknown, the rest not
//! refined; the claim cannot be proven, but can still be disproven, and the
//! asker may go on to hunt a witness in the rest (`hunt`: largest first,
//! as many boxes as the proof tested). Asked to map where it holds over
//! the whole scope (`Precision::maps`: a motion's steps, a sweep's boxes,
//! the pieces of a scope), it sets the box aside and searches the rest, and
//! ends unknown once none is left. A stuck box is set aside either way,
//! never halved: where the boxes that seek a witness as it ends find one,
//! an asker that goes on with a nearer question (query::min_distance) ends
//! unknown at it the next round, unless that question settles it.
//!
//! Each round halves the first BATCH boxes and tests the halves side by
//! side (on a runner, a search's piece: by the runners, a share each,
//! crate::runners), and seeks a witness (a walk from samples, costlier
//! than a test) in the first box alone; outcomes are taken in their order,
//! so the answer is a function of the question, not of the threads. An
//! asker that tests boxes itself (one at a time, or reading what it cannot
//! share across threads) takes the halves of a round and gives their tests
//! back (`round_by`). Ended unknown, the stuck box and the first boxes left
//! seek a witness before it answers. One budget: the boxes it may test; a proof
//! that would test over four times as many in all, at the rate each depth
//! has opened boxes so far, ends there (Knuth's estimate of a search tree's
//! size), a hunt at as many as the proof tested. Past the call's time limit
//! (crate::stop) it ends unknown.

use crate::query::Ended;
use crate::scope;
use ink_core::V3;
use ink_core::eval::{Ctx, span_mid};
use ink_core::interval::I;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::sync::{Arc, Mutex};

/// Boxes a round halves (two to eight times as many tested), and boxes left that seek
/// a witness when the search ends.
pub const BATCH: usize = 128;

/// An axis of a box: one of space's, a ranged param's, or time's.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Axis {
    Space(usize),
    Param(u32),
    T,
}

/// A box of the search: space, the scope beyond it (`ctx`: the ranged
/// params' spans and t's, shared by the halves of space), and what the
/// question keeps for its halves (`state`).
pub struct Node<S> {
    pub b: [I; 3],
    pub ctx: Arc<Ctx>,
    pub state: S,
}

/// How far each axis spreads a test's bounds over a box, in one unit for
/// all of them (a clearance's: length): its rate there times its
/// half-width; and the space axes along which the question is the same all
/// through the box (`still`: its field there reads them not), never halved.
#[derive(Clone, Debug, Default)]
pub struct Smear {
    pub space: [f64; 3],
    pub params: Vec<(u32, f64)>,
    pub t: f64,
    pub still: [bool; 3],
}

/// A box tested.
pub enum Test<S, W> {
    /// Nothing in it violates the claim.
    Settled,
    /// A point violating it.
    Found(W),
    /// Undecided: the bound proven over it, each axis's smear, and the
    /// state its halves start from.
    Open { bound: f64, smear: Smear, state: S },
}

/// A question the search answers.
pub trait Question: Sync {
    type State: Clone + Send + Sync;
    type Witness: Send + Sync;
    /// Decide a box as far as bounds over it can.
    fn test(&self, n: &Node<Self::State>) -> Test<Self::State, Self::Witness>;
    /// Look for a witness in a box: asked of the first box of each round,
    /// and of the first boxes left where the search ends undecided.
    fn seek(&self, n: &Node<Self::State>) -> Option<Self::Witness>;
}

/// How finely a search halves, and what it does with a box no halving
/// decides: space to sides no wider than `finest` (the accuracy asked),
/// each ranged param to 2^-`splits` of its band, t to adjacent floats;
/// such a box ends the proof (a hunt sets it aside), or, where it `maps`
/// where the claim holds over the whole scope, is set aside.
#[derive(Clone, Copy, Debug)]
pub struct Precision {
    pub finest: f64,
    pub splits: u32,
    pub maps: bool,
}

impl Precision {
    /// Whether a claim holds, space halved to `finest`, the params
    /// scope::PARAM_SPLITS times.
    pub fn holds(finest: f64) -> Precision {
        Precision { finest, splits: scope::PARAM_SPLITS, maps: false }
    }
}

/// What a search ends with.
pub enum Outcome<W> {
    Proven,
    Found(W),
    Unknown(Undecided),
}

/// A search ended undecided: the bound it attained (the least over the
/// boxes left), those boxes and their volume, the box that ended it (the
/// stuck one; else the first left), and why it ended.
#[derive(Clone, Debug)]
pub struct Undecided {
    pub bound: f64,
    pub cells: usize,
    pub volume: f64,
    pub cell: Option<[V3; 2]>,
    /// Why the proof ended.
    pub ended: Ended,
    /// The hunt for a witness after it, where one ran (`Search::hunt`).
    pub hunt: Option<Hunt>,
}

impl Undecided {
    /// Whether the call's time limit stopped it, in its proof or its hunt.
    pub fn stopped(&self) -> bool {
        self.ended == Ended::Stopped || self.hunt.is_some_and(|h| h.stopped == HuntEnd::Time)
    }
}

/// A hunt for a witness once the proof ended undecided (`Search::hunt`):
/// the boxes it tested, none a certified violation, and why it stopped.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Hunt {
    pub cells: usize,
    pub stopped: HuntEnd,
}

/// Why a hunt stopped: it had tested as many boxes as the proof (`Cap`),
/// the search's budget (`Budget`), the call's time limit passed (`Time`),
/// or every box left was stuck (`Stuck`: at the refinement limit).
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum HuntEnd {
    Cap,
    Budget,
    Time,
    Stuck,
}

/// A box waiting: its bound, the axis it is halved on (None: stuck) and
/// the space axes halved with it, its volume (of space, times its share of
/// each ranged param's band and of t's), and the order it was made in.
struct Open<S> {
    node: Node<S>,
    /// Halvings from its root.
    depth: usize,
    bound: f64,
    axis: Option<Axis>,
    /// Space axes halved with it (`together`).
    with: [bool; 3],
    volume: f64,
    /// Waiting in a hunt (Search::hunt).
    hunted: bool,
    seq: u64,
}

impl<S> Open<S> {
    /// Whether it comes before `o`: proving, lower bound, then smaller;
    /// hunting, larger, then lower bound; then made first.
    fn first(&self, o: &Self) -> Ordering {
        match self.hunted {
            false => self.bound.total_cmp(&o.bound).then(self.volume.total_cmp(&o.volume)),
            true => o.volume.total_cmp(&self.volume).then(self.bound.total_cmp(&o.bound)),
        }
        .then(self.seq.cmp(&o.seq))
    }
}

impl<S> PartialEq for Open<S> {
    fn eq(&self, o: &Self) -> bool {
        self.seq == o.seq
    }
}

impl<S> Eq for Open<S> {}

impl<S> PartialOrd for Open<S> {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

/// The greatest is popped first.
impl<S> Ord for Open<S> {
    fn cmp(&self, o: &Self) -> Ordering {
        o.first(self)
    }
}

/// A box tested, as it is taken: its depth, smallest finite side, and the
/// test's outcome, open boxes placed (their axes and volume).
type Placed<S, W> = (usize, f64, Test<Open<S>, W>);

/// A search under way.
pub struct Search<S, W> {
    p: Precision,
    open: BinaryHeap<Open<S>>,
    /// Boxes stuck: none waits to be halved. Where the search maps the
    /// whole scope or hunts it goes on without them; else one ends it, at
    /// the round it is taken in or, where that round's witness let the
    /// asker go on (a nearer question, which may settle it), the next.
    aside: Vec<Open<S>>,
    made: u64,
    tested: usize,
    volume: f64,
    res: f64,
    /// t's span over the roots, which a box's share of t is of.
    t_band: Option<(f64, f64)>,
    /// An outcome the roots' tests already gave.
    done: Option<Outcome<W>>,
    /// What each depth has shown (`ahead`).
    depths: Vec<Depth>,
    /// Hunting (`hunt`): the boxes tested when it began, and why the
    /// proof ended.
    hunting: Option<(usize, Ended)>,
    /// Why the proof last ended undecided.
    unproven: Ended,
}

/// At one depth (halvings from a root): boxes tested, of them open, boxes
/// halved and the halves they made, and open boxes waiting.
#[derive(Clone, Default)]
struct Depth {
    tested: u64,
    open: u64,
    halved: u64,
    halves: u64,
    waiting: u64,
}

/// How a round's boxes are found: their witness sought by the question,
/// where it tests them.
type Seek<'a, S, W> = Option<&'a (dyn Fn(&Node<S>) -> Option<W> + Sync)>;

impl<S: Clone + Send + Sync, W: Send + Sync> Search<S, W> {
    /// A search of `roots`, tested at once by `q`, each placed as it is
    /// tested: what a box keeps goes with its test (a middle's search
    /// starts from every box another left).
    pub fn new<Q: Question<State = S, Witness = W>>(q: &Q, roots: Vec<Node<S>>, p: Precision) -> Self {
        let mut s = Search::of(&roots, p);
        let taken: Vec<_> = roots
            .into_par_iter()
            .map(|n| {
                let t = q.test(&n);
                s.place(0, n, t)
            })
            .collect();
        let seek = |n: &Node<S>| q.seek(n);
        s.done = s.take(taken, Some(&seek));
        s
    }

    /// A search of `roots` the asker has tested (in order, up to a witness,
    /// which ends it).
    pub fn start(roots: Vec<(Node<S>, Test<S, W>)>, p: Precision) -> Self {
        let mut s = Search::of(&roots.iter().map(|(n, _)| n).collect::<Vec<_>>(), p);
        let taken: Vec<_> = roots.into_par_iter().map(|(n, t)| s.place(0, n, t)).collect();
        s.done = s.take(taken, None);
        s
    }

    /// A search of `roots`, none yet taken.
    fn of<N: std::borrow::Borrow<Node<S>>>(roots: &[N], p: Precision) -> Self {
        let t_band = roots.iter().filter_map(|n| n.borrow().ctx.t_span).reduce(|a, b| (a.0.min(b.0), a.1.max(b.1)));
        Search {
            p,
            open: BinaryHeap::new(),
            aside: Vec::new(),
            made: 0,
            tested: 0,
            volume: 0.0,
            res: f64::INFINITY,
            t_band,
            done: None,
            depths: Vec::new(),
            hunting: None,
            unproven: Ended::Limit,
        }
    }

    /// Boxes tested so far.
    pub fn tested(&self) -> usize {
        self.tested
    }

    /// Boxes tested hunting a witness (`hunt`).
    fn sought(&self) -> usize {
        self.hunting.map_or(0, |(from, _)| self.tested - from)
    }

    /// The volume of the boxes left (how much is undecided).
    pub fn volume(&self) -> f64 {
        self.volume
    }

    /// The boxes left undecided (waiting, and set aside), in no order.
    pub fn left(&self) -> impl Iterator<Item = &Node<S>> {
        self.open.iter().chain(&self.aside).map(|o| &o.node)
    }

    /// How finely it looked: the smallest side of a box tested (infinite:
    /// none had a finite side).
    pub fn res(&self) -> f64 {
        self.res
    }

    /// One round: the first BATCH boxes halved and their halves tested by
    /// `q`, the first box seeking a witness meanwhile, unless the search
    /// has ended (no box left: proven, or, with boxes set aside, unknown; a
    /// witness; a box stuck; `limit` boxes tested, or over four times as
    /// many in all at the rate it has gone (`ahead`: past its estimate's
    /// errors seen, two times, so a search that will not fit ends in
    /// thousands of boxes, not millions); the time limit).
    pub fn round<Q: Question<State = S, Witness = W>>(&mut self, q: &Q, limit: usize) -> Option<Outcome<W>> {
        let seek = |n: &Node<S>| q.seek(n);
        let mut batch = match self.begin(limit, BATCH, Some(&seek)) {
            Ok(b) => b,
            Err(o) => return Some(o),
        };
        let halves = self.halves(&batch, true);
        // What a box keeps for its halves goes with them, as each is
        // tested; the first box, which seeks, is kept.
        batch.truncate(1);
        let s = &*self;
        let test = |(depth, n): (usize, Node<S>)| {
            let t = q.test(&n);
            s.place(depth, n, t)
        };
        let (taken, found) = match crate::runners::here() {
            // A search's piece (on a runner): the halves, a runner's share
            // at a time, and the first box's seek go to the runners.
            true => {
                let cells: Vec<Mutex<Option<(usize, Node<S>)>>> = halves.into_iter().map(|h| Mutex::new(Some(h))).collect();
                let share = cells.len().div_ceil(crate::runners::count()).max(1);
                let jobs: Vec<Option<&[Mutex<Option<(usize, Node<S>)>>]>> = std::iter::once(None).chain(cells.chunks(share).map(Some)).collect();
                let (taken, found): (Vec<Vec<_>>, Vec<Option<W>>) = crate::runners::each(&jobs, |j| match j {
                    None => (Vec::new(), q.seek(&batch[0].node)),
                    Some(c) => (c.iter().map(|h| test(h.lock().unwrap().take().expect("a half, tested once"))).collect(), None),
                })
                .into_iter()
                .unzip();
                (taken.into_iter().flatten().collect(), found.into_iter().flatten().next())
            }
            false => rayon::join(|| halves.into_par_iter().map(&test).collect::<Vec<_>>(), || q.seek(&batch[0].node)),
        };
        if let Some(w) = found {
            // The round's boxes wait all the same, so an asker can go on
            // with a nearer question (query::min_distance); one stuck ends
            // the next round.
            if let Some(o) = self.take(taken, None) {
                self.done = Some(o);
            }
            return Some(Outcome::Found(w));
        }
        self.take(taken, Some(&seek))
    }

    /// Settle the boxes waiting whose bound is at least `bound`: a question
    /// made nearer (a minimum distance's, once a closer pair is found)
    /// holds there by what they proved.
    pub fn settle_from(&mut self, bound: f64) {
        let open = std::mem::take(&mut self.open).into_vec();
        for o in open {
            if o.bound >= bound {
                self.volume -= o.volume;
                self.depths[o.depth].waiting -= 1;
            } else {
                self.open.push(o);
            }
        }
        self.aside.retain(|o| o.bound < bound);
    }

    /// The claim cannot be proven (the search ended undecided: a box
    /// stuck, or it would not fit its budget), but can still be disproven:
    /// go on for a witness alone. Boxes stuck are set aside, the rest
    /// halved, the largest first, then the lowest bound: breadth first, so
    /// a violation narrower than a box is reached once boxes are that small
    /// everywhere, not after every box of one region is at its finest (a
    /// region no size decides holds the lowest bounds, and would hold the
    /// hunt), until a witness, every box left stuck, the time limit, the
    /// budget, or as many boxes tested as the proof tested (its cap: a
    /// claim held only where measured, no violation to find, would
    /// otherwise spend the whole budget). A claim the proof decides never
    /// hunts. Finding no witness does not show there is none: the boxes it
    /// did not reach, and those stuck, may hold one. Deterministic as any
    /// round.
    pub fn hunt(&mut self) {
        if self.hunting.is_some() {
            return;
        }
        self.hunting = Some((self.tested, self.unproven));
        let open = std::mem::take(&mut self.open).into_vec();
        self.open = open.into_iter().map(|o| Open { hunted: true, ..o }).collect();
    }

    /// Rounds of `q` to the search's end.
    pub fn finish<Q: Question<State = S, Witness = W>>(&mut self, q: &Q, limit: usize) -> Outcome<W> {
        loop {
            if let Some(o) = self.round(q, limit) {
                return o;
            }
        }
    }

    /// Rounds of `most` boxes to the search's end, their halves tested by
    /// the asker one at a time, in order, none past the first witness: it
    /// ends the search, and what the asker would record of a box past it
    /// would be read with it.
    pub fn serially(&mut self, limit: usize, most: usize, mut test: impl FnMut(&Node<S>) -> Test<S, W>) -> Outcome<W> {
        loop {
            let tested = |nodes: &[Node<S>]| {
                let mut tests = Vec::with_capacity(nodes.len());
                for n in nodes {
                    tests.push(test(n));
                    if matches!(tests.last(), Some(Test::Found(_))) {
                        break;
                    }
                }
                tests
            };
            if let Some(o) = self.round_by(limit, most, tested) {
                return o;
            }
        }
    }

    /// A round whose halves the asker tests (`test`, given them in order:
    /// a test for each, or for each up to a witness, which ends the
    /// search), of the first `most` boxes. No witness is sought but in its
    /// tests.
    pub fn round_by(&mut self, limit: usize, most: usize, test: impl FnOnce(&[Node<S>]) -> Vec<Test<S, W>>) -> Option<Outcome<W>> {
        let batch = match self.begin(limit, most.max(1), None) {
            Ok(b) => b,
            Err(o) => return Some(o),
        };
        let (depths, nodes): (Vec<usize>, Vec<Node<S>>) = self.halves(&batch, batch.len() > 8).into_iter().unzip();
        let tests = test(&nodes);
        assert!(tests.len() == nodes.len() || matches!(tests.last(), Some(Test::Found(_))), "a test for each box up to a witness");
        let s = &*self;
        // Tiny rounds (profile-reach) cost less than waking the pool. Keep
        // larger caller-tested rounds parallel, in the same collection order.
        let taken: Vec<_> = if batch.len() <= 8 {
            depths.into_iter().zip(nodes).zip(tests).map(|((depth, n), t)| s.place(depth, n, t)).collect()
        } else {
            depths.into_par_iter().zip(nodes).zip(tests).map(|((depth, n), t)| s.place(depth, n, t)).collect()
        };
        self.take(taken, None)
    }

    /// The first `most` boxes, popped to be halved, unless the search has
    /// ended (round).
    fn begin(&mut self, limit: usize, most: usize, seek: Seek<S, W>) -> Result<Vec<Open<S>>, Outcome<W>> {
        if let Some(o) = self.done.take() {
            return Err(o);
        }
        if let Some(stuck) = self.aside.iter().filter(|_| !self.p.maps && self.hunting.is_none()).map(|o| o.seq).min() {
            return Err(self.end(Ended::Finest, Some(stuck), seek));
        }
        if self.open.is_empty() {
            return Err(match self.aside.is_empty() {
                true => Outcome::Proven,
                false => Outcome::Unknown(self.undecided(Ended::Finest, None)),
            });
        }
        if crate::stop::reached() {
            return Err(Outcome::Unknown(self.undecided(Ended::Stopped, None)));
        }
        // A hunt cannot finish the tree the projection sizes: it ends at
        // the budget, or once it has tested as many boxes as the proof.
        let ended = match self.hunting {
            Some((from, _)) => self.tested - from >= from,
            None => self.forecast() > 4.0 * limit as f64,
        };
        if self.tested >= limit || ended {
            return Err(self.end(Ended::Limit, None, seek));
        }
        let batch: Vec<Open<S>> = (0..most.min(self.open.len())).filter_map(|_| self.open.pop()).collect();
        self.volume -= batch.iter().map(|o| o.volume).sum::<f64>();
        Ok(batch)
    }

    /// The halves of a round's boxes, each with its depth, made side by
    /// side, in order.
    fn halves(&mut self, batch: &[Open<S>], parallel: bool) -> Vec<(usize, Node<S>)> {
        // Caller-tested tiny rounds stay on the caller. A parallel
        // question keeps its splitting parallel too, ready for its tests.
        let halve_all = |o: &Open<S>| {
            let depth = o.depth + 1;
            halve(o).into_iter().map(move |n| (depth, n))
        };
        let halves: Vec<(usize, Node<S>)> = if parallel { batch.par_iter().flat_map_iter(halve_all).collect() } else { batch.iter().flat_map(halve_all).collect() };
        for o in batch {
            let d = &mut self.depths[o.depth];
            (d.waiting, d.halved) = (d.waiting - 1, d.halved + 1);
        }
        for &(depth, _) in &halves {
            self.depths[depth - 1].halves += 1;
        }
        halves
    }

    /// A box tested, and where it would wait: its smallest finite side and
    /// the test's outcome, open boxes with their axes and volume (their
    /// order is the taking's).
    fn place(&self, depth: usize, n: Node<S>, t: Test<S, W>) -> Placed<S, W> {
        let side = (0..3).map(|k| n.b[k].hi - n.b[k].lo).filter(|w| w.is_finite()).fold(f64::INFINITY, f64::min);
        let t = match t {
            Test::Settled => Test::Settled,
            Test::Found(w) => Test::Found(w),
            Test::Open { bound, smear, state } => {
                let node = Node { state, ..n };
                let axis = self.axis(&node, &smear);
                let with = self.together(&node, &smear, axis);
                let volume = self.volume_of(&node);
                let bound = if bound.is_nan() { f64::NEG_INFINITY } else { bound };
                Test::Open { bound, smear, state: Open { node, depth, bound, axis, with, volume, hunted: self.hunting.is_some(), seq: 0 } }
            }
        };
        (depth, side, t)
    }

    /// Boxes tested, taken in order: a witness ends the search; open boxes
    /// wait, and one stuck is set aside, and ends it where it does not map
    /// the whole scope.
    fn take(&mut self, taken: Vec<Placed<S, W>>, seek: Seek<S, W>) -> Option<Outcome<W>> {
        self.tested += taken.len();
        let mut stuck = None;
        for (depth, side, t) in taken {
            self.res = self.res.min(side);
            if self.depths.len() <= depth {
                self.depths.resize(depth + 1, Depth::default());
            }
            self.depths[depth].tested += 1;
            match t {
                Test::Settled => {}
                Test::Found(w) => return Some(Outcome::Found(w)),
                Test::Open { state: mut o, .. } => {
                    o.seq = self.made;
                    self.made += 1;
                    if o.axis.is_none() {
                        if !self.p.maps && self.hunting.is_none() && stuck.is_none() {
                            stuck = Some(o.seq);
                        }
                        self.aside.push(o);
                        continue;
                    }
                    let d = &mut self.depths[depth];
                    (d.open, d.waiting) = (d.open + 1, d.waiting + 1);
                    self.volume += o.volume;
                    self.open.push(o);
                }
            }
        }
        stuck.map(|s| self.end(Ended::Finest, Some(s), seek))
    }

    /// The boxes it would test in all: those tested and those `ahead`
    /// (Knuth's estimate).
    pub fn forecast(&self) -> f64 {
        self.tested as f64 + self.ahead()
    }

    /// The boxes the search would still test, at the rates each depth has
    /// shown so far: an open box at depth d is halved into its share of
    /// the halves boxes there made, of which the share open at depth d + 1
    /// is halved again, and so on (Knuth's estimate of a search tree's
    /// size, by depth).
    fn ahead(&self) -> f64 {
        let (mut below, mut total) = (0.0, 0.0);
        for d in (0..self.depths.len()).rev() {
            let x = &self.depths[d];
            let open = self.depths.get(d + 1).filter(|y| y.tested > 0).map_or(0.0, |y| y.open as f64 / y.tested as f64);
            let halves = if x.halved > 0 { x.halves as f64 / x.halved as f64 } else { 0.0 };
            below = halves * (1.0 + open * below);
            total += x.waiting as f64 * below;
        }
        total
    }

    /// The axis to halve a box on: of those it can be halved on and its
    /// test's bounds spread along (smear above 0; NaN: as large as any), the
    /// one of largest smear, ties to the wider side, then space first, in
    /// order, then t, then the params. Where every space axis they spread
    /// along is at `finest`, t or a param only while its smear is above
    /// every finite side's (finer cells cannot narrow what its span adds; a
    /// side reaching infinity is tested whole, never halved, so it does not
    /// stop them), else the widest side wider than `finest`. None, stuck:
    /// every side at `finest` or infinite, and nothing else to halve.
    fn axis(&self, n: &Node<S>, s: &Smear) -> Option<Axis> {
        let key = |v: f64| if v.is_nan() { f64::INFINITY } else { v };
        let finest = self.p.finest;
        let halvable = |a: Axis| match a {
            Axis::Space(k) => {
                let x = n.b[k];
                let m = x.mid();
                !s.still[k] && x.lo.is_finite() && x.hi.is_finite() && x.hi - x.lo > finest && x.lo < m && m < x.hi
            }
            Axis::Param(p) => n.ctx.p_spans.iter().any(|s| s.k == p && !s.derived && scope::can_halve(s.span, s.band, self.p.splits)),
            Axis::T => n.ctx.t_span.is_some_and(|s| s.0 < span_mid(s) && span_mid(s) < s.1),
        };
        let best = |c: &mut dyn Iterator<Item = (Axis, f64, f64)>| {
            c.filter(|&(a, v, _)| key(v) > 0.0 && halvable(a)).fold(None, |b: Option<(Axis, f64, f64)>, (a, v, w)| match b {
                Some((_, bv, bw)) if !(key(v) > key(bv) || key(v) == key(bv) && w > bw) => b,
                _ => Some((a, v, w)),
            })
        };
        let space = best(&mut (0..3).map(|k| (Axis::Space(k), s.space[k], n.b[k].hi - n.b[k].lo)));
        let scope = best(&mut std::iter::once((Axis::T, s.t, 0.0)).chain(s.params.iter().map(|&(p, v)| (Axis::Param(p), v, 0.0))));
        // Else the widest side not yet at `finest`: what the test does not
        // see spread (the other shape's extent) may still be what keeps it
        // undecided.
        let widest = || (0..3).filter(|&k| halvable(Axis::Space(k))).max_by(|&i, &j| (n.b[i].hi - n.b[i].lo).total_cmp(&(n.b[j].hi - n.b[j].lo))).map(Axis::Space);
        match (space, scope) {
            (Some(x), Some(p)) => Some(if key(p.1) > key(x.1) { p.0 } else { x.0 }),
            (Some(x), None) => Some(x.0),
            (None, Some(p)) if key(p.1) > (0..3).filter(|&k| (n.b[k].hi - n.b[k].lo).is_finite()).map(|k| key(s.space[k])).fold(0.0, f64::max) => Some(p.0),
            _ => widest(),
        }
    }

    /// The space axes a box is halved on with space axis k: each the next
    /// halving would choose anyway, its smear at least half the largest
    /// (halving an axis at least halves its smear), so no box is tested
    /// that would only wait for its next halving. Multisection, as interval
    /// branch and bound does where tests cost more than boxes.
    fn together(&self, n: &Node<S>, s: &Smear, axis: Option<Axis>) -> [bool; 3] {
        let Some(Axis::Space(k)) = axis else { return [false; 3] };
        std::array::from_fn(|j| {
            let x = n.b[j];
            j != k
                && !s.still[j]
                && s.space[j] > 0.0
                && s.space[j] >= 0.5 * s.space[k]
                && x.lo.is_finite()
                && x.hi.is_finite()
                && x.hi - x.lo > self.p.finest
                && x.lo < x.mid()
                && x.mid() < x.hi
        })
    }

    /// The search ends undecided (`stuck`: the box that ended it), once
    /// that box and the first others left have sought a witness, side by
    /// side, where the question seeks them: the first found, in that
    /// order, is the answer.
    fn end(&mut self, ended: Ended, stuck: Option<u64>, seek: Seek<S, W>) -> Outcome<W> {
        if let Some(seek) = seek {
            let mut first: Vec<&Open<S>> = self.open.iter().chain(self.aside.iter().filter(|o| Some(o.seq) == stuck)).collect();
            first.sort_unstable_by(|x, y| (Some(x.seq) != stuck).cmp(&(Some(y.seq) != stuck)).then(x.first(y)));
            first.truncate(BATCH);
            if let Some(w) = first.par_iter().map(|o| seek(&o.node)).collect::<Vec<_>>().into_iter().flatten().next() {
                return Outcome::Found(w);
            }
        }
        Outcome::Unknown(self.undecided(ended, stuck))
    }

    /// The search ended undecided (`ended`: the proof's reason, or where
    /// it hunts the hunt's: Finest, every box left stuck; Limit, its cap or
    /// the budget; Stopped, the time limit), and what it left.
    fn undecided(&mut self, ended: Ended, stuck: Option<u64>) -> Undecided {
        let (ended, hunt) = match self.hunting {
            Some((from, why)) => {
                let stopped = match ended {
                    Ended::Stopped => HuntEnd::Time,
                    Ended::Finest => HuntEnd::Stuck,
                    Ended::Limit if self.sought() >= from => HuntEnd::Cap,
                    Ended::Limit => HuntEnd::Budget,
                };
                (why, Some(Hunt { cells: self.sought(), stopped }))
            }
            None => (ended, None),
        };
        self.unproven = ended;
        let left = || self.open.iter().chain(&self.aside);
        // A hunt's: the first box stuck, which ended the proof.
        let stuck = stuck.or(self.hunting.and(self.aside.first().map(|o| o.seq)));
        let cell = left().find(|o| Some(o.seq) == stuck).or(self.open.peek()).or(self.aside.first()).map(|o| [o.node.b.map(|x| x.lo), o.node.b.map(|x| x.hi)]);
        let bound = left().map(|o| o.bound).fold(f64::INFINITY, f64::min);
        Undecided { bound, cells: self.open.len() + self.aside.len(), volume: left().map(|o| space_volume(&o.node.b)).sum(), cell, ended, hunt }
    }

    /// A box's volume: the product over its axes of its share of each
    /// (space's sides as they are; a ranged param's span of its band; t's
    /// span of t's over the roots), along each axis of some finite width
    /// only: one of none (a point: a param narrowed to one value, as a
    /// middle's search narrows its input) or reaching infinity (all of
    /// space) is no axis of its search.
    fn volume_of(&self, n: &Node<S>) -> f64 {
        let width = |(a, b): (f64, f64)| Some(b - a).filter(|w| w.is_finite() && *w > 0.0);
        let share = |s: (f64, f64), whole: Option<(f64, f64)>| width(s).map_or(1.0, |w| w / whole.and_then(width).unwrap_or(1.0));
        let space = n.b.iter().map(|x| share((x.lo, x.hi), None));
        let params = n.ctx.p_spans.iter().filter(|s| !s.derived).map(|s| share(s.span, Some(s.band)));
        space.chain(params).chain(n.ctx.t_span.map(|s| share(s, self.t_band))).product()
    }
}

/// A box's halves on its axis (and those halved with it), each with its
/// state.
fn halve<S: Clone>(o: &Open<S>) -> Vec<Node<S>> {
    let (b, n) = (o.node.b, &o.node);
    match o.axis.expect("a stuck box is set aside, never waits to be halved") {
        Axis::Space(k) => {
            let mut boxes = vec![b];
            for j in (0..3).filter(|&j| j == k || o.with[j]) {
                let m = b[j].mid();
                boxes = boxes
                    .into_iter()
                    .flat_map(|x| {
                        let (mut lo, mut hi) = (x, x);
                        lo[j] = I::new(x[j].lo, m);
                        hi[j] = I::new(m, x[j].hi);
                        [lo, hi]
                    })
                    .collect();
            }
            boxes.into_iter().map(|b| Node { b, ctx: n.ctx.clone(), state: n.state.clone() }).collect()
        }
        Axis::Param(p) => {
            let (lo, hi) = scope::halves(&n.ctx, p);
            [lo, hi].into_iter().map(|c| Node { b, ctx: Arc::new(c), state: n.state.clone() }).collect()
        }
        Axis::T => {
            let s = n.ctx.t_span.expect("a time span");
            let m = span_mid(s);
            [(s.0, m), (m, s.1)].into_iter().map(|h| Node { b, ctx: Arc::new(Ctx { t: span_mid(h), t_span: Some(h), ..(*n.ctx).clone() }), state: n.state.clone() }).collect()
        }
    }
}

fn space_volume(b: &[I; 3]) -> f64 {
    b.iter().map(|x| x.hi - x.lo).product()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;
    use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

    /// What tests made their boxes keep, alive, and the most at once.
    static LIVE: AtomicUsize = AtomicUsize::new(0);
    static MOST: AtomicUsize = AtomicUsize::new(0);

    struct Made;

    impl Made {
        fn new() -> Made {
            MOST.fetch_max(LIVE.fetch_add(1, Relaxed) + 1, Relaxed);
            Made
        }
    }

    impl Drop for Made {
        fn drop(&mut self) {
            LIVE.fetch_sub(1, Relaxed);
        }
    }

    /// A box's state, which its test fills (an emptiness box's tape, made
    /// as it is tested).
    #[derive(Clone)]
    struct Held(Arc<Mutex<Option<Made>>>);

    struct Fill;

    impl Question for Fill {
        type State = Held;
        type Witness = ();
        fn test(&self, n: &Node<Held>) -> Test<Held, ()> {
            *n.state.0.lock().unwrap() = Some(Made::new());
            Test::Settled
        }
        fn seek(&self, _: &Node<Held>) -> Option<()> {
            None
        }
    }

    /// Small caller-tested rounds split and place on their caller.
    /// Both sides of the threshold retain exactly the serial halving order.
    #[test]
    fn small_rounds_split_on_the_caller_and_large_rounds_keep_the_order() {
        struct Tracked(Arc<Mutex<Vec<std::thread::ThreadId>>>);
        impl Clone for Tracked {
            fn clone(&self) -> Self {
                self.0.lock().unwrap().push(std::thread::current().id());
                Self(self.0.clone())
            }
        }
        impl Drop for Tracked {
            fn drop(&mut self) {
                self.0.lock().unwrap().push(std::thread::current().id());
            }
        }
        let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        pool.install(|| {
            let caller = std::thread::current().id();
            for count in [1, 8, 9, BATCH] {
                let threads = Arc::new(Mutex::new(Vec::new()));
                let roots = || {
                    (0..count)
                        .map(|i| {
                            let n = Node {
                                b: [I { lo: i as f64 * 4.0, hi: i as f64 * 4.0 + 2.0 }, I::point(0.0), I::point(0.0)],
                                ctx: Arc::new(Ctx::new(0.0, vec![])),
                                state: Tracked(threads.clone()),
                            };
                            (n, Test::<_, ()>::Open { bound: 0.0, smear: Smear { space: [1.0, 0.0, 0.0], ..Default::default() }, state: Tracked(threads.clone()) })
                        })
                        .collect()
                };
                let mut expected = Search::start(roots(), Precision::holds(0.01));
                let batch = expected.begin(10000, count, None).unwrap_or_else(|_| panic!("open roots"));
                let want: Vec<_> = batch.iter().flat_map(|o| halve(o).into_iter().map(|n| n.b)).collect();
                let mut s = Search::start(roots(), Precision::holds(0.01));
                threads.lock().unwrap().clear();
                let mut got = Vec::new();
                let outcome = s.round_by(10000, count, |nodes| {
                    got.extend(nodes.iter().map(|n| n.b));
                    nodes.iter().map(|_| Test::Settled).collect()
                });
                let outcome = outcome.or_else(|| s.round_by(10000, count, |_| panic!("all halves settled")));
                assert!(matches!(outcome, Some(Outcome::Proven)));
                assert_eq!(got, want, "{count} parents: identical ordered boxes");
                if count <= 8 {
                    let ids = threads.lock().unwrap();
                    assert!(!ids.is_empty() && ids.iter().all(|id| *id == caller), "small rounds stay on the caller: {ids:?}");
                }
            }
        });
    }

    /// A stuck box whose search ends with a witness (the boxes left seek
    /// one when a stuck box ends it) stays set aside: an asker that goes on
    /// with a nearer question (min_distance) ends unknown at that box, or
    /// goes on past it once the nearer question settles it, and never
    /// halves it (it panicked: `ink distance` of a robot arm's tool part).
    #[test]
    fn a_stuck_box_ended_by_a_witness_is_never_halved() {
        // Boxes of x halved to 0.25, each bound by where it starts; a
        // witness is found in a box at the finest only, while `finds`.
        struct Near {
            finds: bool,
        }
        impl Question for Near {
            type State = ();
            type Witness = ();
            fn test(&self, n: &Node<()>) -> Test<(), ()> {
                Test::Open { bound: n.b[0].lo, smear: Smear { space: [1.0, 0.0, 0.0], ..Default::default() }, state: () }
            }
            fn seek(&self, n: &Node<()>) -> Option<()> {
                (self.finds && n.b[0].hi - n.b[0].lo <= 0.25).then_some(())
            }
        }
        let search = || {
            let root = Node { b: [I::new(0.0, 1.0), I::point(0.0), I::point(0.0)], ctx: Arc::new(Ctx::new(0.0, vec![])), state: () };
            let mut s = Search::new(&Near { finds: true }, vec![root], Precision::holds(0.25));
            assert!(s.round(&Near { finds: true }, 1000).is_none(), "halves at 0.5");
            assert!(matches!(s.round(&Near { finds: true }, 1000), Some(Outcome::Found(()))), "quarters stuck, a witness in one");
            s
        };
        // Nearer, two of the four quarters settled: the first stuck one left ends it.
        let mut s = search();
        s.settle_from(0.5);
        match s.round(&Near { finds: false }, 1000) {
            Some(Outcome::Unknown(u)) => {
                assert_eq!((u.ended, u.bound, u.cells), (Ended::Finest, 0.0, 2));
                assert_eq!(u.cell.map(|c| (c[0][0], c[1][0])), Some((0.0, 0.25)));
            }
            _ => panic!("unknown at the stuck box"),
        }
        // Nearer still, all settled: proven.
        let mut s = search();
        s.settle_from(0.0);
        assert!(matches!(s.round(&Near { finds: false }, 1000), Some(Outcome::Proven)));
    }

    /// A claim no box decides anywhere (every box open, its bound tying,
    /// as sin^2 + cos^2 - 1 does), failing only in boxes no wider than
    /// 2^-10 that hold 0.9: depth first, its proof ends at the first box
    /// stuck at 2^-16, near 0, never reaching 0.9; the hunt goes on breadth
    /// first and finds it, or ends once it has tested as many boxes as the
    /// proof, the stuck box reported.
    #[test]
    fn a_stuck_box_ends_the_proof_and_the_hunt_finds_a_witness_elsewhere() {
        struct Fails(Option<f64>);
        impl Question for Fails {
            type State = ();
            type Witness = f64;
            fn test(&self, n: &Node<()>) -> Test<(), f64> {
                let x = n.b[0];
                match self.0 {
                    Some(a) if x.hi - x.lo <= 1.0 / 1024.0 && x.lo <= a && a <= x.hi => Test::Found(a),
                    _ => Test::Open { bound: 0.0, smear: Smear { space: [1.0, 0.0, 0.0], ..Default::default() }, state: () },
                }
            }
            fn seek(&self, _: &Node<()>) -> Option<f64> {
                None
            }
        }
        let run = |q: &Fails| {
            let root = Node { b: [I::new(0.0, 1.0), I::point(0.0), I::point(0.0)], ctx: Arc::new(Ctx::new(0.0, vec![])), state: () };
            let mut s = Search::new(q, vec![root], Precision::holds(1.0 / 65536.0));
            let Outcome::Unknown(u) = s.finish(q, 100_000) else { panic!("the proof ends at a stuck box") };
            assert_eq!((u.ended, u.hunt, u.cell.map(|c| c[0][0] < 0.5)), (Ended::Finest, None, Some(true)));
            let proof = s.tested();
            s.hunt();
            (s.finish(q, 100_000), proof, s)
        };
        let (found, proof, s) = run(&Fails(Some(0.9)));
        assert!(matches!(found, Outcome::Found(0.9)), "the hunt finds it");
        assert!(proof < 4096 && s.sought() < 10_000, "{proof} boxes, then {} hunting", s.sought());
        let (none, proof, s) = run(&Fails(None));
        let Outcome::Unknown(u) = none else { panic!("none to find") };
        // As many boxes as the proof tested, to the round's end.
        let hunt = u.hunt.expect("a hunt");
        assert_eq!((u.ended, hunt.stopped, hunt.cells), (Ended::Finest, HuntEnd::Cap, s.sought()));
        assert!(proof <= hunt.cells && hunt.cells < proof + 2 * BATCH, "{proof} boxes, then {} hunting", hunt.cells);
        assert!(u.cell.is_some_and(|c| c[1][0] - c[0][0] <= 1.0 / 65536.0), "the box that ended the proof: {:?}", u.cell);
    }

    /// A search's roots go as they are tested, with what their tests made:
    /// a middle's search starts from every box another left (thousands),
    /// and each made its tape as tested (gigabytes over a motion's steps
    /// checked side by side, when all were kept until the last was tested).
    #[test]
    fn a_root_goes_with_its_test() {
        let roots = (0..4096).map(|_| Node { b: [I { lo: f64::NEG_INFINITY, hi: f64::INFINITY }; 3], ctx: Arc::new(Ctx::new(0.0, vec![])), state: Held(Arc::default()) }).collect();
        let s = Search::new(&Fill, roots, Precision::holds(1.0));
        assert!(s.tested() == 4096 && MOST.load(Relaxed) <= rayon::current_num_threads(), "{} made at once", MOST.load(Relaxed));
    }
}
