//! Emptiness and interference proofs, with witnesses and tolerance bounds.

use super::bounds::{bounds_quick, max_extent};
use super::cover::{Cover, far};
use super::proof::{Confirm, Ended, Verdict, box_samples, confirm_inside, contract, in_span, scope_shares, share, signed, worked};
use crate::scope::{self, At};
use crate::search;
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_interval, eval_point, eval_points};
use ink_core::interval::I;
use ink_core::sign::inside;
use ink_core::tape::{BOTH, Tape};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

/// An emptiness search's answer: its verdict, the cell size it reached,
/// the tolerance it was asked to (mm; 0: exactly), and whether a proof
/// needed it (`exact` false: some cell was shown only to hold no point
/// deeper than it, so the solid may not be empty; `Proof::exact`), with
/// up to LOOSE such cells, the lowest bounds first, where a point of the
/// solid may be (`loose`).
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Proof {
    pub verdict: Verdict,
    pub res: f64,
    pub tol: f64,
    pub exact: bool,
    pub loose: Vec<Loose>,
    /// A disproof's field value at the deepest point found walking down from its witness (`deepen`):
    /// how deep what violates it goes at least, the witness left as found. None: no deeper point.
    #[serde(default)]
    pub deepest: Option<f64>,
}

/// A cell a proof settled only by its tolerance: its box (lo, hi) and the
/// piece of the scope it was settled over, each ranged param's span
/// (`params`, those derived from them left out) and time's over a step
/// (`t`).
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Loose {
    pub cell: [V3; 2],
    pub params: Vec<(u32, (f64, f64))>,
    pub t: Option<(f64, f64)>,
}

impl Loose {
    fn of(b: [I; 3], ctx: &Ctx) -> Loose {
        let params = ctx.p_spans.iter().filter(|s| !s.derived).map(|s| (s.k, s.span)).collect();
        Loose { cell: [b.map(|x| x.lo), b.map(|x| x.hi)], params, t: ctx.t_span }
    }

    /// ctx narrowed to the cell's scope.
    pub fn scope(&self, ctx: &Ctx) -> Ctx {
        let mut c = self.params.iter().fold(ctx.clone(), |c, &(k, span)| c.with_param(k, span));
        if self.t.is_some() {
            c.t_span = self.t;
        }
        c.into_scope()
    }
}

/// How many cells a proof that needed its tolerance keeps (Proof::loose).
pub const LOOSE: usize = 64;

/// Prove that a field's solid S in lo..hi (where it is below 0, as clipped
/// to lo..hi: a point past it is outside) holds no point deeper than `tol`
/// (every point of S within tol of a point outside it; 0: no point at all),
/// or find a witness (a point whose cube of half-side tol is inside: so is
/// the ball of radius tol about it, and it is deeper): the one search
/// (crate::search) of boxes of the material's box (a proven bound of S in
/// lo..hi), halved down to cells of the floor (positive: the claim's
/// precision, half its tolerance, or for an exact question a floor set by
/// the whole model's size), or of a `res` given where that is finer. Its
/// boxes lie on a lattice whose cells are the first cells (`first_cells`:
/// 1/32 of the material's box, coarsened to the floor times a power of two)
/// halved to the floor, so a claim over a whole part is decided as finely
/// as one over a feature, and a region under 32 floors keeps its own 1/32,
/// finer than the floor. Best first, the box
/// where the field can be lowest is refined first, so a witness is
/// usually found in the first few hundred boxes; a box undecided with
/// every side at the finest cells ends the proof (Ended::Finest), with the
/// accuracy it holds to (`attained`: no point of S deeper than that, from
/// the boxes left; at or under `tol`, it is proven to it). Else the search
/// hunts a point deeper than `tol` in the rest (search::Search::hunt), and
/// is unknown where it finds none. At most BOXES boxes are tested
/// (Ended::Limit).
/// Returns the answer, the cells it reached (the smallest side of a box
/// tested) and whether a proof needed the tolerance (`exact` false).
pub fn prove_empty(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, res: Option<f64>, floor: f64, tol: f64) -> Proof {
    empty_now(tape, None, ctx, lo, hi, res, floor, tol).0
}

/// `prove_empty` (of "inside a and not inside b" where `except` is given), testing at most `boxes`
/// boxes, and as many again hunting a point (Ended::Limit past them).
#[allow(clippy::too_many_arguments)]
pub fn prove_empty_within(tape: &Arc<Tape>, except: Option<&Except>, ctx: &Ctx, lo: V3, hi: V3, res: Option<f64>, floor: f64, tol: f64, boxes: usize) -> Proof {
    empty_within(tape, except, ctx, lo, hi, res, floor, tol, boxes).0
}

/// The two shapes of a question "inside a and not inside b" (a subset's
/// failures): its search reads h = max(a, -b), whose bounds and slopes are
/// a's and b's own, and a and b themselves where h alone cannot say: h < 0
/// is such a point, h = 0 one only where a is below 0 (b reads 0 there,
/// which is not inside b).
#[derive(Clone)]
pub struct Except {
    pub a: Arc<Tape>,
    pub b: Arc<Tape>,
}

/// `prove_empty` of "inside a and not inside b" (`Except`), `tape` being
/// max(a, -b).
#[allow(clippy::too_many_arguments)]
pub fn prove_empty_except(tape: &Arc<Tape>, except: &Except, ctx: &Ctx, lo: V3, hi: V3, res: Option<f64>, floor: f64, tol: f64) -> Proof {
    empty_now(tape, Some(except), ctx, lo, hi, res, floor, tol).0
}

/// `prove_empty`, and the boxes it tested.
#[allow(clippy::too_many_arguments)]
pub(super) fn empty_now(tape: &Arc<Tape>, except: Option<&Except>, ctx: &Ctx, lo: V3, hi: V3, res: Option<f64>, floor: f64, tol: f64) -> (Proof, usize) {
    empty_within(tape, except, ctx, lo, hi, res, floor, tol, BOXES)
}

/// `empty_now`, testing at most `boxes` boxes.
#[allow(clippy::too_many_arguments)]
fn empty_within(tape: &Arc<Tape>, except: Option<&Except>, ctx: &Ctx, lo: V3, hi: V3, res: Option<f64>, floor: f64, tol: f64, boxes: usize) -> (Proof, usize) {
    assert!(floor > 0.0, "an emptiness search needs a positive floor to stop at");
    let ctx = Arc::new(in_span(ctx));
    let proven = Proof { verdict: Verdict::Proven, res: max_extent(lo, hi) / 32.0, tol, exact: true, loose: Vec::new(), deepest: None };
    // Its points lie in a's material (where h reads 0 too).
    let held = except.map_or(tape, |e| &e.a);
    let Some((blo, bhi, _)) = bounds_quick(held, &ctx, lo, hi) else { return (proven, 0) };
    let (first, floor) = first_cells(max_extent(blo, bhi), res, floor);
    // Narrowed once more, not again inside the narrowed box.
    let Some((blo, bhi, _)) = bounds_quick(held, &ctx, blo, bhi) else { return (proven, 0) };
    let mut finest = first;
    while finest > floor {
        finest *= 0.5;
    }
    let Some(root) = lattice(blo, bhi, finest) else {
        return (
            Proof {
                verdict: Verdict::Unknown { cells: 1, volume: f64::INFINITY, ended: Ended::Limit, cell: None, attained: None, hunt: None },
                res: finest,
                tol,
                exact: false,
                loose: Vec::new(),
                deepest: None,
            },
            0,
        );
    };
    let q = Empty {
        tol,
        region: std::array::from_fn(|k| I::new(blo[k], bhi[k])),
        asked: std::array::from_fn(|k| I::new(lo[k], hi[k])),
        loose: Default::default(),
        t_band: ctx.t_span,
        read: tape.space_read(),
        finest,
        sections: Default::default(),
    };
    let state = (Kept::of(tape), Kept::of(tape), except.map(|e| (Kept::of(&e.a), Kept::of(&e.b))));
    let (mut searched, s) = empty_run(&q, vec![search::Node { b: root, ctx: ctx.clone(), state }], finest, boxes);
    if let Verdict::Disproven { witness, t, value, params } = &searched.verdict {
        let (mut at, mut deepest) = (*witness, *value);
        deepen(tape, &ctx, (*t, params), &q, &mut at, &mut deepest);
        searched.deepest = (deepest < *value).then_some(deepest);
    }
    let sections: usize = q.sections.lock().unwrap().values().flatten().map(|x| x.tested).sum();
    (searched, s.tested() + sections)
}

/// A point walked from a witness as deep into what violates the claim as its field leads (at its moment and
/// params, within the region asked), each step kept only where bounds show it deeper: the depth it gives is
/// the overlap's, not the first point's, and still one the overlap reaches at least.
fn deepen(tape: &Tape, ctx: &Ctx, (t, params): (f64, &[(u32, f64)]), q: &Empty, witness: &mut V3, value: &mut f64) {
    let mut at = params.iter().fold(ctx.clone(), |c, &(k, x)| c.with_param(k, (x, x)));
    (at.t, at.t_span) = (t, None);
    let at = at.point();
    let within = |p: V3| (0..3).all(|k| q.region[k].lo <= p[k] && p[k] <= q.region[k].hi && q.asked[k].lo <= p[k] && p[k] <= q.asked[k].hi);
    let scale = (0..3).map(|k| q.region[k].hi - q.region[k].lo).fold(0.0f64, f64::max);
    let mut step = value.abs().max(1e-6 * scale);
    for _ in 0..DEEPEN_STEPS {
        let g = eval_grads(tape, &at, &[*witness])[0];
        let n = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
        if !(n > 0.0 && n.is_finite()) {
            return;
        }
        // Down its slope, else along an axis (where two shapes' fields meet, the slope flips from side to side).
        let ways: [V3; 7] = [[-g[1] / n, -g[2] / n, -g[3] / n], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
        let deeper = loop {
            if step < 1e-9 * scale {
                return;
            }
            let found = ways.iter().find_map(|w| {
                let x: V3 = std::array::from_fn(|k| witness[k] + step * w[k]);
                let hi = ink_core::eval::eval_interval(tape, &at, x.map(I::point), None).hi;
                (within(x) && hi < *value).then_some((x, hi))
            });
            if let Some(d) = found {
                break d;
            }
            step *= 0.5;
        };
        (*witness, *value) = deeper;
        step *= 2.0;
    }
}

/// How many steps a witness is walked deeper at most (`deepen`).
const DEEPEN_STEPS: usize = 256;

/// The emptiness search `q` from `roots` (boxes of a lattice of cells
/// `finest`), and its answer.
fn empty_run(q: &Empty, roots: Vec<search::Node<Pieces>>, finest: f64, boxes: usize) -> (Proof, search::Search<Pieces, Witness>) {
    // The lattice's sides are finest times powers of two: a side up to half
    // as wide again is at the finest cells (its ends rounded).
    let mut s = search::Search::new(q, roots, search::Precision::holds(1.5 * finest));
    let mut outcome = s.finish(q, boxes);
    // Not proven, and not to the tolerance by the boxes left (`attained`):
    // a point of S deeper than it may be there, so the search hunts one.
    let mut attained = None;
    if let search::Outcome::Unknown(u) = &outcome
        && !u.stopped()
    {
        let a = q.attained(&s.left().collect::<Vec<_>>());
        attained = Some(a);
        if !(q.tol > 0.0 && a <= q.tol) {
            s.hunt();
            outcome = s.finish(q, boxes);
            attained = None;
        }
    }
    // The cells it reached: the lattice's, the smallest side tested (as
    // rounded) taken to the power of two of the finest it is.
    let res = finest * (s.res() / finest).log2().round().exp2();
    let loose: Vec<Loose> = q.loose.lock().unwrap().iter().map(|c| c.1.clone()).collect();
    let verdict = match outcome {
        search::Outcome::Proven => Verdict::Proven,
        search::Outcome::Found((witness, at, value)) => Verdict::Disproven { witness, t: at.t, value, params: at.params },
        search::Outcome::Unknown(u) => {
            let left: Vec<&search::Node<Pieces>> = s.left().collect();
            let attained = attained.or_else(|| (!u.stopped()).then(|| q.attained(&left)));
            if attained.is_some_and(|a| q.tol > 0.0 && a <= q.tol) {
                // Every box left is held to the tolerance only.
                let loose = left.iter().filter_map(|n| Some(Loose::of(q.clip(n.b)?, &n.ctx))).take(LOOSE).collect();
                return (Proof { verdict: Verdict::Proven, res, tol: q.tol, exact: false, loose, deepest: None }, s);
            }
            let volume = left.iter().filter_map(|n| q.clip(n.b)).map(|b| b.iter().map(|x| x.hi - x.lo).product::<f64>()).sum();
            let cell = u.cell.and_then(|[a, b]| q.clip(std::array::from_fn(|k| I::new(a[k], b[k])))).map(|b| [b.map(|x| x.lo), b.map(|x| x.hi)]);
            Verdict::Unknown { cells: u.cells, volume, ended: u.ended, cell, attained, hunt: u.hunt }
        }
    };
    (Proof { verdict, res, tol: q.tol, exact: loose.is_empty(), loose, deepest: None }, s)
}

/// Whether `section` (a field reading no y) is shown nowhere below 0 at y = 0
/// over radii x in `r` and heights z in 0..`lead`, exactly (no tolerance),
/// by an emptiness search of at most SECTION_BOXES boxes, its cells down to
/// `floor`; and the boxes it tested.
fn section_holds(section: &Arc<Tape>, ctx: &Arc<Ctx>, r: I, lead: f64, floor: f64) -> (bool, usize) {
    let (lo, hi) = ([r.lo, 0.0, 0.0], [r.hi, 0.0, lead]);
    let Some((blo, bhi, _)) = bounds_quick(section, ctx, lo, hi) else { return (true, 0) };
    let mut finest = max_extent(blo, bhi) / 32.0;
    while finest > floor {
        finest *= 0.5;
    }
    let Some(root) = lattice(blo, bhi, finest) else { return (false, 0) };
    let q = Empty {
        tol: 0.0,
        region: std::array::from_fn(|k| I::new(blo[k], bhi[k])),
        asked: std::array::from_fn(|k| I::new(lo[k], hi[k])),
        loose: Default::default(),
        t_band: ctx.t_span,
        read: section.space_read(),
        finest,
        sections: Default::default(),
    };
    let state = (Kept::of(section), Kept::of(section), None);
    let mut s = search::Search::new(&q, vec![search::Node { b: root, ctx: ctx.clone(), state }], search::Precision::holds(1.5 * finest));
    let held = matches!(s.finish(&q, SECTION_BOXES), search::Outcome::Proven);
    (held, s.tested())
}

/// The cube of a lattice of cells of side h anchored at lo that holds
/// lo..hi: its side h times the least power of two that reaches hi on
/// every axis (a cube, so halving every axis at once keeps it one: an
/// octree's cells). None where a side is not finite or needs more than
/// 2^62 cells.
pub(super) fn lattice(lo: V3, hi: V3, h: f64) -> Option<[I; 3]> {
    let mut side = h;
    for _ in 0..=62 {
        let ends: V3 = std::array::from_fn(|k| lo[k] + side);
        if !(0..3).all(|k| lo[k].is_finite() && ends[k].is_finite()) {
            return None;
        }
        if (0..3).all(|k| ends[k] >= hi[k]) {
            return Some(std::array::from_fn(|k| I::new(lo[k], ends[k])));
        }
        side *= 2.0;
    }
    None
}

/// Boxes a search may test (an emptiness's or a clearance's).
pub const BOXES: usize = 1 << 22;

/// What h reading at least 0 over a box of "inside a and not inside b"
/// leaves there (`Empty::shallow`): no point of the question's, none deeper
/// than the tolerance, or it cannot say.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Allowed {
    Exact,
    Loose,
    No,
}

/// Witness samples a box confirms against the tolerance (the deepest
/// first: each costs a bound over a cube).
const DEEP_TRIES: usize = 4;

/// A point inside, its moment and ranged params, and the field's value there.
pub(super) type Witness = (V3, At, f64);

/// What a box of an emptiness search keeps for its halves: the field
/// simplified over it, and over it grown by the tolerance (where a
/// witness's cube and the points that cover it are bounded); for "inside
/// a and not inside b", a and b simplified over that grown box, where
/// every point its test reads them at lies.
pub(super) type Pieces = (Kept, Kept, Option<(Kept, Kept)>);

/// A tape a box keeps for its halves: the tape its test read and the
/// choices made there, a step at a time, simplified when one of its halves
/// is first tested (cloned for them, they share it). A box waiting holds
/// choices, not a tape: the tapes alive are those of boxes whose halves
/// wait (an emptiness's eight, depth first; a clearance's boxes, best
/// first, wait apart from their siblings and keep their own) and of those
/// being tested.
#[derive(Clone)]
pub(super) struct Kept(Arc<std::sync::Mutex<(Arc<Tape>, Vec<Vec<u8>>)>>);

impl Kept {
    pub(super) fn of(t: &Arc<Tape>) -> Kept {
        Kept::after(t, Vec::new())
    }

    /// `from` simplified by each of `steps` in turn.
    pub(super) fn after(from: &Arc<Tape>, steps: Vec<Vec<u8>>) -> Kept {
        let steps = steps.into_iter().filter(|ch| ch.iter().any(|&c| c != BOTH)).collect();
        Kept(Arc::new(std::sync::Mutex::new((from.clone(), steps))))
    }

    /// The tape, simplified once and kept.
    pub(super) fn tape(&self) -> Arc<Tape> {
        let mut k = self.0.lock().unwrap();
        for ch in std::mem::take(&mut k.1) {
            k.0 = k.0.narrowed(&ch);
        }
        k.0.clone()
    }

    /// The tape, not kept (a box read once, as the search ends).
    fn read(&self) -> Arc<Tape> {
        let k = self.0.lock().unwrap();
        k.1.iter().fold(k.0.clone(), |t, ch| t.narrowed(ch))
    }

    /// The tape it holds: its own once made, else the one it is made from.
    #[cfg(test)]
    pub(super) fn held(&self) -> *const Tape {
        Arc::as_ptr(&self.0.lock().unwrap().0)
    }
}

/// Emptiness as a question of the one search: no point of the solid S in
/// the material's box `region` (a proven bound of S in the region asked,
/// `asked`, where the tape holds) deeper than `tol`. Boxes are a lattice's
/// (prove_empty), each tested as clipped to `region`. A box is settled
/// where a test passes, each sound for any field:
/// - the field's bounds over it are at least 0 (over the scope's axes,
///   ranged params and time over a step of a motion, also first order in
///   them, and the facing form, whose pairs cancel the rates of two faces
///   one param or a motion moves together: a shaft and its bore of one
///   diameter);
/// - its first-order bound is, where no sample of it is inside (across a
///   gap thinner than the box, or at a branch cut);
/// - every point of it lies within the tolerance of a point shown outside
///   S: by a pair of faces facing each other along an axis, at any box
///   size (centred::centred_settles: faces that touch, or overlap thinner
///   than the tolerance), or by points shown outside near each of its
///   samples (`covered`, at boxes of a few tolerances): none of its solid
///   is deeper (`loose`: a proof that needed the tolerance, which keeps
///   the LOOSE such boxes of lowest bound).
///
/// Those after the bounds ask the box narrowed to where the field can be
/// at most 0 by its terms of one axis alone (`contract`: an extrusion's
/// span), of the field's sign alone there (`signed`: a slab below 0 inside
/// its span drops out of an extrusion), and a box is never halved along an
/// axis that reads not (an extrusion's walls inside its span: the plane's
/// question, a point shown outside standing for its whole column).
///
/// A witness is its centre where the field is below 0 all over it and the
/// cube about it fits (or is shown inside), or one of its samples, the
/// deepest first, whose cube is shown inside: at the box's moment and
/// ranged params' middle values (its own: time halved inside a step finds
/// a witness at a moment of its half), or each param pushed to the end of
/// its span that lowers the field there (scope::pushed: a stack's worst
/// corner).
///
/// Open, a box proves no depth short of the tolerance (its depth is
/// bounded only for the boxes left: `attained`), so its bound is 0, and
/// open boxes tie: the smaller goes first, so the search is depth first,
/// as interval solvers prove a set empty (IBEX, RealPaver), its memory as
/// its depth, not as a contact's area (best first on the field's lower
/// bound refines a contact breadth first, every box near it bounded
/// alike), and a box no size decides is reached in as many rounds as
/// halvings. Witnesses are sought in every box tested. Its smear on each space
/// axis is its half-width on the lattice (the field's rate taken alike
/// along each: a cube is halved as an octree's cell is); where it is
/// decided at the middle of the scope's axes but not over their spans,
/// what the spans add is the cause, and space's smear is the field's
/// spread over the box (each axis's rate bound times its half-width),
/// shared out by width, each param's and time's its blur of the bound
/// (`scope_shares`): a step's moments are halved where they blur the box
/// more than its size does. Covered at the middle but not over the spans,
/// the axis halved least so far is halved first.
pub(super) struct Empty {
    pub(super) tol: f64,
    pub(super) region: [I; 3],
    pub(super) asked: [I; 3],
    pub(super) loose: std::sync::Mutex<Vec<(f64, Loose)>>,
    /// Time's span over the search (a step of a motion), which a box's
    /// share of it is of.
    pub(super) t_band: Option<(f64, f64)>,
    /// The space axes the field reads (Tape::space_read): along one it
    /// reads not, the material box ends nowhere (`to_faces`).
    pub(super) read: [bool; 3],
    /// The lattice's cells (empty_run): a box with a side over 1.5 of
    /// them can still be halved there.
    pub(super) finest: f64,
    /// The cross-sections asked so far (`by_section`), by their tapes'
    /// instructions as written.
    pub(super) sections: std::sync::Mutex<HashMap<String, Vec<Section>>>,
}

/// A cross-section's answer (`Empty::by_section`), kept for the boxes
/// after it: the scope it was asked over (the context itself, held so no
/// other takes its place: boxes halved in space share their context, a
/// param's or time's halves each have their own), the radii, whether its
/// field was shown nowhere below 0 there, and the boxes its search tested.
pub(super) struct Section {
    ctx: Arc<Ctx>,
    r: I,
    held: bool,
    tested: usize,
}

/// Cross-sections of one tape in one scope that may fail before it is
/// asked no more.
const SECTION_FAILS: usize = 4;

/// Boxes a cross-section's search may test.
const SECTION_BOXES: usize = 1 << 14;

impl search::Question for Empty {
    type State = Pieces;
    type Witness = Witness;

    fn test(&self, n: &search::Node<Pieces>) -> search::Test<Pieces, Witness> {
        use search::Test::{Found, Open, Settled};
        let Some(b) = self.clip(n.b) else { return Settled };
        let ctx = &*n.ctx;
        let (from, given) = (n.state.0.tape(), n.state.1.tape());
        let mut ch = vec![BOTH; from.len()];
        let r = eval_interval(&from, ctx, b, Some(&mut ch));
        // a and b simplified over the grown box, once read, and the choices
        // that simplified them.
        let cell = std::cell::OnceCell::new();
        let narrowed = || {
            cell.get_or_init(|| {
                n.state.2.as_ref().map(|(a, z)| {
                    let (a, z) = (a.tape(), z.tape());
                    let (mut cha, mut chz) = (vec![BOTH; a.len()], vec![BOTH; z.len()]);
                    eval_interval(&a, ctx, self.grown(b), Some(&mut cha));
                    eval_interval(&z, ctx, self.grown(b), Some(&mut chz));
                    (Except { a: a.narrowed(&cha), b: z.narrowed(&chz) }, (a, cha), (z, chz))
                })
            })
            .as_ref()
        };
        let except = || narrowed().map(|n| &n.0);
        // h at least 0 over the box settles where that leaves no point of
        // the question's there, or none deeper than the tolerance, the box
        // then kept as one held to it (`shallow`); "inside a and not inside
        // b" holds a point reading 0 only where a is below 0.
        let by_h = std::cell::OnceCell::new();
        let by_h = |x: [I; 3]| *by_h.get_or_init(|| Self::by_h(except(), ctx, x));
        let shallow = std::cell::OnceCell::new();
        let shallow = || *shallow.get_or_init(|| if n.state.2.is_none() { Allowed::Exact } else { self.shallow(except(), &by_h, ctx, b) });
        // Its value over the box holds over any part of it (a narrowed
        // box); a box it settles held to the tolerance is kept as such.
        let settled = |a: Allowed, at: [I; 3]| {
            if a == Allowed::Loose {
                self.keep_loose(r.lo, at, ctx);
            }
            Settled
        };
        if r.clear_of_surface() {
            return Settled;
        }
        if r.none_inside() && shallow() != Allowed::No {
            return settled(shallow(), b);
        }
        // a and b are narrowed only where h may be below 0 throughout.
        if r.hi <= 0.0
            && let Some(w) = self.found_except(except(), r, ctx, b)
        {
            return Found(w);
        }
        let sub = from.narrowed(&ch);
        // No point of it below 0 over the box, shown on a cross-section,
        // settles as its bounds at least 0 would.
        if self.by_section(&sub, &n.ctx, b) && shallow() != Allowed::No {
            return settled(shallow(), b);
        }
        // `wide` simplified over the box grown by the tolerance where the
        // box's own tape has halved the one it was tested with (narrowing
        // costs a bound over the grown box).
        let mut chw = vec![BOTH; given.len()];
        let wide = match self.tol > 0.0 {
            false => sub.clone(),
            true if 2 * sub.len() <= given.len() => {
                eval_interval(&given, ctx, self.grown(b), Some(&mut chw));
                given.narrowed(&chw)
            }
            true => given.clone(),
        };
        let c: V3 = std::array::from_fn(|k| b[k].mid());
        // Inside throughout: its centre is a witness where the box holds the
        // cube about it, else where that cube is inside.
        if r.all_inside() && ((0..3).all(|k| c[k] - b[k].lo >= self.tol && b[k].hi - c[k] >= self.tol) || self.deep(&sub, &wide, &mut Confirm::new(ctx), c)) {
            return Found((c, At::of(ctx), eval_point(&sub, ctx, c)));
        }
        // The rest is asked where the field can be at most 0 alone
        // (`contract`: past there it is above 0), of its sign alone there
        // (`signed`): settled and covered so, and never halved along the
        // axes that reads not.
        let Some(narrow) = contract(&sub, ctx, b) else { return Settled };
        let mut signed = signed(&sub, ctx, narrow);
        if narrow != b {
            let rs = eval_interval(&signed, ctx, narrow, None);
            if rs.clear_of_surface() {
                return Settled;
            }
            if rs.none_inside() && shallow() != Allowed::No {
                return settled(shallow(), narrow);
            }
        }
        let b = narrow;
        // The scope's axes the field reads: ranged params, and time over a
        // step of a motion.
        let (ps, ts) = (scope::splittable(&signed, ctx, scope::PARAM_SPLITS), scope::spans_t(&signed, ctx));
        // A max of groups of terms reading params apart (scope::groups) is
        // at least 0 over a box only where one group is, and halving a
        // group's params narrows that group's bound alone: each box made
        // must still be settled by one group, so the boxes multiply with
        // each group halved. Only a group whose next halving settles a half
        // by itself (the group at least 0 over it) is halved, the one
        // nearest settling (its bound at the scope's middle less what its
        // params spread it by, `shares`). Where none is, space is halved
        // first while that can narrow a group (one reads an axis of space
        // the box can still be halved along); else the nearest is halved.
        // The params of the rest are held.
        let held = |signed: &Tape, mid: &Ctx, shares: &[(u32, f64)]| -> Vec<u32> {
            let gs = scope::groups(signed, ctx, &ps);
            if gs.is_empty() {
                return Vec::new();
            }
            // The param the search would halve first: the largest share, then
            // the least halved.
            let next =
                |g: &[u32]| g.iter().copied().max_by(|&i, &j| share(shares, i).total_cmp(&share(shares, j)).then(scope::halvings(ctx, j).total_cmp(&scope::halvings(ctx, i))));
            let halves = |g: &(Vec<u32>, Tape)| {
                next(&g.0).is_some_and(|k| {
                    let (lo, hi) = scope::halves(ctx, k);
                    eval_interval(&g.1, &lo, b, None).lo >= 0.0 || eval_interval(&g.1, &hi, b, None).lo >= 0.0
                })
            };
            let near = |g: &(Vec<u32>, Tape)| eval_interval(&g.1, mid, b, None).lo - g.0.iter().map(|&k| share(shares, k)).sum::<f64>();
            let nearest = |only: &dyn Fn(&(Vec<u32>, Tape)) -> bool| {
                gs.iter()
                    .enumerate()
                    .filter(|g| only(g.1))
                    .map(|(j, g)| (j, near(g)))
                    .fold(None, |a: Option<(usize, f64)>, (j, x)| if a.is_some_and(|a| !(x > a.1)) { a } else { Some((j, x)) })
            };
            let halvable = |k: usize| n.b[k].hi - n.b[k].lo > 1.5 * self.finest;
            let space_narrows = gs.iter().any(|g| g.1.space_read().iter().enumerate().any(|(k, &r)| r && halvable(k)));
            let best = nearest(&halves).or_else(|| if space_narrows { None } else { nearest(&|_| true) });
            gs.iter().enumerate().filter(|&(j, _)| Some(j) != best.map(|x| x.0)).flat_map(|(_, g)| g.0.iter().copied()).collect()
        };
        let mid = (!ps.is_empty() || ts).then(|| scope::middle(ctx, ts));
        let (mut params, mut t, mut rates, mut centred) = (Vec::new(), 0.0, None, false);
        if let Some(mid) = &mid {
            let m = eval_interval(&signed, mid, b, None);
            if m.all_inside() && self.deep(&sub, &wide, &mut Confirm::new(mid), c) {
                return Found((c, At::of(mid), eval_point(&sub, mid, c)));
            }
            if m.none_inside() {
                let mut chp = vec![BOTH; signed.len()];
                let (r1, blur) = ink_core::eval::eval_interval_scope_blur(&signed, ctx, b, Some(&mut chp));
                if r1.none_inside() && shallow() != Allowed::No {
                    return settled(shallow(), b);
                }
                signed = signed.narrowed(&chp);
                if shallow() != Allowed::No && crate::centred::centred_at_least(&[(&signed, 0.0)], ctx, b, 0.0) {
                    return settled(shallow(), b);
                }
                centred = true;
                let (shares, tb) = scope_shares(ctx, &[(&blur, 1.0)]);
                let held = held(&signed, mid, &shares);
                params = ps.iter().filter(|k| !held.contains(k)).map(|&k| (k, share(&shares, k))).collect();
                t = if ts { tb } else { 0.0 };
                rates = Some(ink_core::eval::eval_interval_grad(&signed, ctx, b));
            }
        }
        let still = signed.space_read().map(|r| !r);
        let pts = box_samples(b, 4);
        let mut v = Vec::new();
        eval_points(&sub, ctx, &pts, &mut v);
        let within = |k: usize| (0..3).all(|a| b[a].contains(pts[k][a]));
        // The deepest samples below 0.
        let mut below: Vec<usize> = (0..pts.len()).filter(|&k| inside(v[k]) && within(k)).collect();
        below.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
        let tries = if self.tol == 0.0 { below.len() } else { DEEP_TRIES };
        let mut confirm = (!below.is_empty()).then(|| Confirm::new(ctx));
        if let Some(&k) = below.iter().take(tries).find(|&&k| confirm.as_mut().is_some_and(|c| self.deep(&sub, &wide, c, pts[k]))) {
            return Found((pts[k], At::of(ctx), v[k]));
        }
        // A sample inside leaves no lower bound over the box to try (a
        // speed rule: skipping a bound never makes a wrong answer); with a
        // tolerance, a pair of faces that touch or overlap thinner than it
        // still covers the box, whatever its size.
        use crate::centred::Settles;
        let bounded = below.is_empty() && !centred;
        // Split where pieces meet flush (ink_core::sign::air), no point of
        // the box is the question's; for "inside a and not inside b", none
        // but where b reads 0 and a is below it (`shallow`, or `by_h` over
        // the box: none there); where b is solid over the box grown by the
        // tolerance (pieces of it meeting flush, its 0 of no volume: no face
        // of b there), points of no depth, held to the tolerance.
        if bounded && ink_core::sign::air(&sub, ctx, b) {
            if shallow() != Allowed::No {
                return settled(shallow(), b);
            }
            if let Some(e) = except()
                && Self::by_h(Some(e), ctx, b)
                && let at @ (Allowed::Exact | Allowed::Loose) = self.at_edge(e, ctx, b)
            {
                return settled(at, b);
            }
            // b as narrowed over the box grown by the tolerance holds there.
            if self.tol > 0.0 && except().is_some_and(|e| ink_core::sign::solid(&e.b, ctx, self.grown(b))) {
                self.keep_loose(r.lo, b, ctx);
                return Settled;
            }
        }
        // A point within the tolerance of one where h is at least 0 is within
        // it of the question's outside only where `by_h` decides.
        let settles = match self.tol > 0.0 {
            _ if shallow() == Allowed::No => Settles::No,
            true => match crate::centred::centred_settles(&signed, ctx, b, self.tol, bounded) {
                Settles::Covered if !by_h(self.grown(b)) => Settles::No,
                s => s,
            },
            false if bounded && crate::centred::centred_at_least(&[(&signed, 0.0)], ctx, b, 0.0) => Settles::AtLeast,
            false => Settles::No,
        };
        match settles {
            Settles::AtLeast => return settled(shallow(), b),
            Settles::Covered => {
                self.keep_loose(r.lo, b, ctx);
                return Settled;
            }
            Settles::No => {}
        }
        // Each case of the switches the box leaves undecided bounded apart
        // (a slope's clamp read in several terms: each case's rate one
        // number, not every rate between its branches' in each).
        if bounded
            && sub.ssa.iter().any(|x| matches!(x.op, ink_core::tape::Op::Select))
            && ink_core::eval::eval_interval_cases(&sub, ctx, b).none_inside()
            && shallow() != Allowed::No
        {
            return settled(shallow(), b);
        }
        let sides = b.map(|x| x.hi - x.lo);
        let cover = Cover { still, ..self.cover() };
        if cover.covered(&signed, &|| &wide, &|r, x, p| Self::out_at(except(), r, x, p), ctx, b, sides, &pts, &v) {
            self.keep_loose(r.lo, b, ctx);
            return Settled;
        }
        // Covered at the scope's middle but not over it: the axis of it
        // halved least so far is halved first (time's halvings of the
        // search's span).
        if let Some(mid) = mid.as_ref().filter(|_| params.is_empty() && t == 0.0)
            && cover.covered(&signed, &|| &wide, &|r, x, p| Self::out_at(except(), r, x, p), mid, b, sides, &pts, &v)
        {
            let held = held(&signed, mid, &[]);
            let halved = ps.iter().filter(|k| !held.contains(k)).map(|&k| (Some(k), scope::halvings(ctx, k))).chain(ts.then(|| (None, self.t_halvings(ctx))));
            match halved.min_by(|x, y| x.1.total_cmp(&y.1)).map(|x| x.0) {
                Some(Some(k)) => params.push((k, f64::INFINITY)),
                Some(None) => t = f64::INFINITY,
                None => {}
            }
        }
        if !ps.is_empty() || !scope::read_params(&sub, ctx).is_empty() {
            let mut order: Vec<usize> = (0..pts.len()).filter(|&k| within(k)).collect();
            order.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
            for &k in order.iter().take(4) {
                let Some((pushed, x)) = scope::pushed(&[&sub], ctx, pts[k]) else { continue };
                if inside(x[0]) && self.deep(&sub, &wide, &mut Confirm::new(&pushed), pts[k]) {
                    return Found((pts[k], At::of(&pushed), x[0]));
                }
            }
        }
        // Space is halved as an octree halves its cells (on the lattice's
        // widths: cubes, which the tolerance's cover settles), its spread
        // shared out by width where params compete with it.
        let hw: V3 = std::array::from_fn(|k| if still[k] { 0.0 } else { 0.5 * (n.b[k].hi - n.b[k].lo) });
        let space = match &rates {
            // Settled at the scope's middle (where rates are taken), the
            // scope's spread is what keeps it open; where space's is not
            // bounded (a gear's teeth folded at its pitch) the two cannot
            // be weighed, so space waits until the scope is halved as far
            // as it goes.
            Some(g) if !(0..3).all(|k| ink_core::interval::mag(g.d[k]).is_finite()) && (t > 0.0 || params.iter().any(|p| !(p.1 <= 0.0))) => [0.0; 3],
            Some(g) => {
                let spread = (0..3).map(|k| ink_core::interval::mag(g.d[k]) * 0.5 * sides[k]).sum::<f64>();
                hw.map(|h| h * spread / hw.iter().sum::<f64>())
            }
            None => hw,
        };
        // The halves start from the field simplified over the box, and over
        // it grown by the tolerance only where that has halved the one it
        // was tested with (no tolerance reads no `wide`).
        let wide = if self.tol > 0.0 && 2 * wide.len() <= given.len() { Kept::after(&given, vec![chw]) } else { n.state.1.clone() };
        // a and b as simplified here where they were read (the halves' grown
        // boxes lie in this one's), else as they came.
        let except = match cell.into_inner() {
            Some(Some((_, (a, cha), (z, chz)))) => Some((Kept::after(&a, vec![cha]), Kept::after(&z, vec![chz]))),
            _ => n.state.2.clone(),
        };
        Open { bound: 0.0, smear: search::Smear { space, params, t, still }, state: (Kept::after(&from, vec![ch]), wide, except) }
    }

    /// From the deepest samples of the box, a step down the field's slope:
    /// asked exactly, to just under its surface (`confirm_inside`: an
    /// overlap thinner than the samples' spacing, which none lies in);
    /// with a tolerance, to where the cube about it would fit were the
    /// field a distance (its half-diagonal below 0), confirmed (`deep`).
    fn seek(&self, n: &search::Node<Pieces>) -> Option<Witness> {
        let b = self.clip(n.b)?;
        let ctx = &*n.ctx;
        let (sub, wide) = (&n.state.0.tape(), &n.state.1.tape());
        let pts = box_samples(b, 4);
        let mut v = Vec::new();
        eval_points(sub, ctx, &pts, &mut v);
        let mut order: Vec<usize> = (0..pts.len()).collect();
        order.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
        let inside = |p: &V3| (0..3).all(|a| b[a].contains(p[a]));
        let mut confirm = Confirm::new(ctx);
        order.into_iter().take(DEEP_TRIES).find_map(|k| {
            let p = match self.tol {
                0.0 => confirm_inside(sub, ctx, pts[k], max_extent(b.map(|x| x.lo), b.map(|x| x.hi))).filter(inside)?,
                tol => {
                    let [f, gx, gy, gz] = eval_grads(sub, ctx, &[pts[k]])[0];
                    let g2 = gx * gx + gy * gy + gz * gz;
                    let s = (f + 1.01 * 3f64.sqrt() * tol * g2.sqrt()) / g2;
                    let p = [pts[k][0] - s * gx, pts[k][1] - s * gy, pts[k][2] - s * gz];
                    (s.is_finite() && inside(&p) && self.deep(sub, wide, &mut confirm, p)).then_some(p)?
                }
            };
            Some((p, At::of(ctx), eval_point(sub, ctx, p)))
        })
    }
}

impl Empty {
    /// Whether box b (off the z axis) holds no point of the field `sub`
    /// (its tape over b) below 0, shown on a cross-section: where a
    /// max of some of sub's terms (`Tape::conjuncts`, a bound of it from
    /// below) is kept by every screw motion of one lead L about z
    /// (ink_core::screw: two threads of one pitch on one axis, in any
    /// phase), its value at each point of b is its value on the half-plane
    /// y = 0 at the point's radius and a height within one lead, z in 0..L.
    /// So that max nowhere below 0 there, over the radii of b (or wider: the
    /// material's box out from half b's least radius, asked first, which
    /// the boxes after it share), is sub nowhere below 0 over b: proven
    /// exactly, at the cost of one cross-section for every box of the
    /// engagement.
    pub(super) fn by_section(&self, sub: &Tape, ctx: &Arc<Ctx>, b: [I; 3]) -> bool {
        use ink_core::graph::{BinOp, Var};
        use ink_core::tape::Op;
        if !sub.ssa.iter().any(|x| matches!(x.op, Op::Bin(BinOp::Atan2))) {
            return false;
        }
        let radius = |x: f64, y: f64| x.hypot(y);
        let near = |x: I| x.lo.max(-x.hi).max(0.0);
        let far = |x: I| x.lo.abs().max(x.hi.abs());
        // Rounded outward by more than hypot's error.
        let lo = radius(near(b[0]), near(b[1])) * (1.0 - 1e-15);
        let hi = radius(far(b[0]), far(b[1])) * (1.0 + 1e-15);
        if !(lo > 0.0 && hi.is_finite()) {
            return false;
        }
        let Some(s) = ink_core::screw::screw(sub) else { return false };
        let terms: Vec<(usize, bool)> = sub.conjuncts().into_iter().filter(|t| s.invariant(t.0)).map(|t| (t.0, t.1)).collect();
        if terms.is_empty() {
            return false;
        }
        let kept = sub.cone_terms(&terms);
        let ys: Vec<(usize, f64)> = kept.ssa.iter().enumerate().filter(|x| matches!(x.1.op, Op::Var(Var::Y))).map(|x| (x.0, 0.0)).collect();
        let section = Arc::new(kept.with_consts(&ys));
        let key: String = section.ssa.iter().map(|x| format!("{:?} {} {} {};", x.op, x.a, x.b, x.c)).collect();
        let r = I::new(lo, hi);
        let wide = I::new(0.5 * lo, hi.max(radius(far(self.region[0]), far(self.region[1])) * (1.0 + 1e-15)));
        // Held over radii holding these, or failed over radii within them
        // or too often: decided without a search.
        let (known, tried) = {
            let all = self.sections.lock().unwrap();
            // Only an answer over this very scope: a param's other values,
            // or other moments, may meet where these do not.
            let mine: Vec<&Section> = all.get(&key).map_or(&[][..], |v| &v[..]).iter().filter(|x| Arc::ptr_eq(&x.ctx, ctx)).collect();
            let within = |x: &Section, w: I| x.r.lo <= w.lo && w.hi <= x.r.hi;
            let known = if mine.iter().any(|x| x.held && within(x, r)) {
                Some(true)
            } else if mine.iter().filter(|x| !x.held).count() >= SECTION_FAILS || mine.iter().any(|x| !x.held && r.lo <= x.r.lo && x.r.hi <= r.hi) {
                Some(false)
            } else {
                None
            };
            (known, mine.iter().any(|x| within(x, wide)))
        };
        if let Some(k) = known {
            return k;
        }
        let asked = |w: I| -> bool {
            let floor = self.finest.min(s.lead.lo) / 64.0;
            let (held, tested) = section_holds(&section, ctx, w, s.lead.hi, floor);
            self.sections.lock().unwrap().entry(key.clone()).or_default().push(Section { ctx: ctx.clone(), r: w, held, tested });
            held
        };
        (!tried && asked(wide)) || asked(r)
    }

    /// For "inside a and not inside b", what h reading at least 0 over box
    /// x leaves of the question's points there: 0s of b where a is below 0
    /// (h reads 0). None where a has no point inside x (Exact). Where every
    /// 0 of b over x grown by the tolerance has points above 0 beside it
    /// (Tape::zeros_rise), or `by_h` (given) decides over that box, each has
    /// points inside a and outside b beside it, h below 0 there, which no
    /// box settles outright: it settles as `at_edge` says. A 0 of b with no
    /// outside of b beside it is a point of the question's of no depth,
    /// which no test here excuses: the search refines it, or holds it to
    /// the tolerance where a cover does.
    fn shallow(&self, except: Option<&Except>, by_h: &dyn Fn([I; 3]) -> bool, ctx: &Ctx, x: [I; 3]) -> Allowed {
        let Some(e) = except else { return Allowed::Exact };
        worked(&e.a);
        if eval_interval(&e.a, ctx, x, None).none_inside() {
            return Allowed::Exact;
        }
        let grown = self.grown(x);
        worked(&e.b);
        // b for its sign alone over that box (`signed`): its 0s and their
        // sides are b's.
        let sb = signed(&e.b, ctx, grown);
        if !(sb.zeros_rise(&ink_core::eval::eval_interval_grad_every(&sb, ctx, grown)) || by_h(grown)) {
            return Allowed::No;
        }
        self.at_edge(e, ctx, x)
    }

    /// How a box x settles whose question's points (0s of b where a is
    /// below 0) each have points inside a and outside b beside them (h
    /// below 0, which no box settles outright): Exact where a has no point
    /// inside on each face of the region asked x reaches (a clipped to it
    /// reads 0 there), so those points lie in the region and in a box;
    /// else, with a tolerance, Loose (such 0s hold no ball, so none is
    /// deeper than it); else No.
    fn at_edge(&self, e: &Except, ctx: &Ctx, x: [I; 3]) -> Allowed {
        let faces = (0..3).flat_map(|k| [(k, self.asked[k].lo, x[k].lo <= self.asked[k].lo), (k, self.asked[k].hi, x[k].hi >= self.asked[k].hi)]);
        let inside = faces.filter(|f| f.2).all(|(k, at, _)| {
            let mut face = x;
            face[k] = I::point(at);
            worked(&e.a);
            eval_interval(&e.a, ctx, face, None).none_inside()
        });
        if inside {
            Allowed::Exact
        } else if self.tol > 0.0 {
            Allowed::Loose
        } else {
            Allowed::No
        }
    }

    /// Whether every point of box x where b reads 0 and a is below 0 has
    /// points inside b beside it, so lies on the edge of the question's
    /// solid (a point within the tolerance of it is within the tolerance
    /// of the outside): where every part of x (halved up to three times
    /// where needed) has no point inside a, or no 0 of b (bounds of one
    /// sign), or b rising or falling throughout along some axis (its 0 then
    /// has no volume, and is b's surface).
    fn by_h(except: Option<&Except>, ctx: &Ctx, x: [I; 3]) -> bool {
        fn parts(e: &Except, ctx: &Ctx, x: [I; 3], depth: u32) -> bool {
            // b's 0 there is that of its tape for its sign (`signed`), at
            // any size of part; b first, the container (a is often a whole
            // assembly).
            let e = &Except { a: e.a.clone(), b: signed(&e.b, ctx, x) };
            let rb = eval_interval(&e.b, ctx, x, None);
            if rb.clear_of_surface() || rb.all_inside() {
                return true;
            }
            let g = ink_core::eval::eval_interval_grad(&e.b, ctx, x);
            if (0..3).any(|k| g.d[k].lo > 0.0 || g.d[k].hi < 0.0) {
                return true;
            }
            worked(&e.a);
            let mut ch = vec![BOTH; if depth > 0 { e.a.len() } else { 0 }];
            if eval_interval(&e.a, ctx, x, (depth > 0).then_some(&mut ch)).none_inside() {
                return true;
            }
            if depth == 0 {
                return false;
            }
            // The parts read a as simplified over x.
            let e = Except { a: e.a.narrowed(&ch), b: e.b.clone() };
            (0..8).all(|o: usize| parts(&e, ctx, std::array::from_fn(|k| if o >> k & 1 == 0 { I::new(x[k].lo, x[k].mid()) } else { I::new(x[k].mid(), x[k].hi) }), depth - 1))
        }
        except.is_none_or(|e| parts(e, ctx, x, 3))
    }

    /// Whether every point of box p (a point, or one spanning a box's still
    /// axes) is shown not one of the question's, h's bounds there `r` (or
    /// those of a tape of its sign there): at least 0, and for "inside a
    /// and not inside b" above 0 or a not inside.
    fn out_at(except: Option<&Except>, r: I, ctx: &Ctx, p: [I; 3]) -> bool {
        r.none_inside() && (r.clear_of_surface() || except.is_none_or(|e| eval_interval(&e.a, ctx, p, None).none_inside()))
    }

    /// The field over box b (`sub`, simplified over it) for its sign alone
    /// (`signed`), and the space axes it reads not there: along them the
    /// question is the same all through b (an extrusion's span, inside its
    /// profile's walls), so b is never halved on them.
    fn signed(sub: &Arc<Tape>, ctx: &Ctx, b: [I; 3]) -> (Arc<Tape>, [bool; 3]) {
        let t = signed(sub, ctx, b);
        let read = t.space_read();
        (t, read.map(|r| !r))
    }

    /// For "inside a and not inside b", a witness where h is at most 0 and a
    /// below 0 over box x (b at least 0: not inside it), at x's centre,
    /// where x holds the tolerance's cube about it or that cube is so.
    fn found_except(&self, except: Option<&Except>, r: I, ctx: &Ctx, x: [I; 3]) -> Option<Witness> {
        let e = except?;
        if !(r.hi <= 0.0 && eval_interval(&e.a, ctx, x, None).all_inside()) {
            return None;
        }
        let c: V3 = std::array::from_fn(|k| x[k].mid());
        let holds = (0..3).all(|k| c[k] - x[k].lo >= self.tol && x[k].hi - c[k] >= self.tol);
        let cube: [I; 3] = std::array::from_fn(|k| I::new((c[k] - self.tol).next_down(), (c[k] + self.tol).next_up()));
        let deep = holds
            || (0..3).all(|k| self.region[k].lo <= cube[k].lo && cube[k].hi <= self.region[k].hi)
                && eval_interval(&e.a, ctx, cube, None).all_inside()
                && eval_interval(&e.b, ctx, cube, None).none_inside();
        deep.then(|| (c, At::of(ctx), eval_point(&e.a, ctx, c).max(-eval_point(&e.b, ctx, c))))
    }

    /// Keep box `b` over ctx's scope, settled only by the tolerance with
    /// `lo` its field's lower bound, among the LOOSE of lowest bound (ties
    /// by place, then scope, so the cells kept do not depend on the order
    /// threads settle them in).
    fn keep_loose(&self, lo: f64, b: [I; 3], ctx: &Ctx) {
        let key = |c: &(f64, Loose)| -> Vec<f64> {
            let [l, h] = c.1.cell;
            let scope = c.1.t.iter().chain(c.1.params.iter().map(|p| &p.1)).flat_map(|s| [s.0, s.1]);
            [c.0, l[0], l[1], l[2], h[0], h[1], h[2]].into_iter().chain(scope).collect()
        };
        let cell = (lo, Loose::of(b, ctx));
        let new = key(&cell);
        let mut kept = self.loose.lock().unwrap();
        let at = kept.partition_point(|c| key(c).iter().zip(&new).map(|(x, y)| x.total_cmp(y)).find(|o| o.is_ne()).is_none_or(|o| o.is_lt()));
        if at < LOOSE {
            kept.insert(at, cell);
            kept.truncate(LOOSE);
        }
    }

    /// How many times time's span has been halved in ctx from the search's
    /// (log2 of their widths' ratio).
    fn t_halvings(&self, ctx: &Ctx) -> f64 {
        let w = |s: Option<(f64, f64)>| s.map_or(1.0, |(a, b)| b - a);
        (w(self.t_band) / w(ctx.t_span)).log2()
    }

    /// A box of the lattice as clipped to the material's box (None: it
    /// holds none of it).
    fn clip(&self, b: [I; 3]) -> Option<[I; 3]> {
        let c: [I; 3] = std::array::from_fn(|k| I { lo: b[k].lo.max(self.region[k].lo), hi: b[k].hi.min(self.region[k].hi) });
        c.iter().all(|x| x.lo <= x.hi).then_some(c)
    }

    /// Box `b` grown by the tolerance, within the region asked (where the
    /// tape holds).
    fn grown(&self, b: [I; 3]) -> [I; 3] {
        self.cover().grown(b)
    }

    /// Whether p (in box b) is deeper than the tolerance: the cube of
    /// half-side tol about it is inside, at ctx's point of the scope
    /// (`wide`: the field over b grown by tol); with no tolerance, p itself
    /// is inside. The cube must lie in the material's box (outside it there
    /// is none).
    fn deep(&self, sub: &Tape, wide: &Tape, confirm: &mut Confirm, p: V3) -> bool {
        if self.tol == 0.0 {
            return confirm.below(sub, p, 0.0);
        }
        let cube: [I; 3] = std::array::from_fn(|k| I::new((p[k] - self.tol).next_down(), (p[k] + self.tol).next_up()));
        (0..3).all(|k| self.region[k].lo <= cube[k].lo && cube[k].hi <= self.region[k].hi) && confirm.inside(wide, cube)
    }

    /// The tolerance's rule (Cover) for this search's solid.
    fn cover(&self) -> Cover {
        Cover { tol: self.tol, region: self.region, asked: self.asked, sign: 1.0, read: self.read, still: self.read.map(|r| !r) }
    }

    /// How far a point of the material's box in box b (as clipped) can lie
    /// from the material box's faces, beyond which S has no point: per
    /// axis, the most the nearer face can be away, and of the axes the
    /// least (rounded up). Per axis that is the most of min(p - lo, hi -
    /// p) over the box's p: half the material box's width where its real
    /// middle is in the box, else the far side's distance to the face
    /// beyond it, so the least of the three, never a sample at a float
    /// middle (which may lie on a face: 1e16 + 1 is no float).
    fn to_faces(&self, b: [I; 3]) -> f64 {
        use ink_core::interval::{add_up, mul_up};
        (0..3)
            .filter(|&k| self.read[k])
            .map(|k| {
                let (r, x) = (self.region[k], b[k]);
                let half = mul_up(add_up(r.hi, -r.lo), 0.5);
                add_up(x.hi, -r.lo).min(add_up(r.hi, -x.lo)).min(half).max(0.0)
            })
            .fold(f64::INFINITY, f64::min)
    }

    /// How deep a point of S in box n can be (its distance from the
    /// outside): for each of the box's 4 x 4 x 4 sub-boxes, the farthest
    /// it reaches from the nearest of its samples shown outside (the field
    /// at least 0 there over the whole scope), or from the material box's
    /// faces where nearer; the most of those.
    fn depth(&self, n: &search::Node<Pieces>) -> f64 {
        let Some(b) = self.clip(n.b) else { return 0.0 };
        let (ctx, read) = (&*n.ctx, n.state.0.read());
        let except = n.state.2.as_ref().map(|(a, z)| Except { a: a.read(), b: z.read() });
        let Some(b) = contract(&read, ctx, b) else { return 0.0 };
        let (sub, still) = Self::signed(&read, ctx, b);
        let pts = box_samples(b, 4);
        let half: V3 = b.map(|x| 0.125 * (x.hi - x.lo));
        let mut v = Vec::new();
        eval_points(&sub, ctx, &pts, &mut v);
        let mut slots = Vec::new();
        let within = |p: &V3| (0..3).all(|k| b[k].contains(p[k]));
        // Each shown outside along the box's span of its still axes.
        let column = |p: V3| -> [I; 3] { std::array::from_fn(|k| if still[k] { b[k] } else { I::point(p[k]) }) };
        let out: Vec<V3> = (0..pts.len())
            .filter(|&k| {
                let at = column(pts[k]);
                !inside(v[k]) && within(&pts[k]) && Self::out_at(except.as_ref(), ink_core::eval::eval_interval_in(&sub, ctx, at, None, &mut slots), ctx, at)
            })
            .map(|k| pts[k])
            .collect();
        pts.iter()
            .map(|&c| {
                let around: [I; 3] = std::array::from_fn(|k| I::new((c[k] - half[k]).next_down().max(b[k].lo), (c[k] + half[k]).next_up().min(b[k].hi)));
                out.iter().map(|&q| far(q, c, half, b, still)).fold(self.to_faces(around), f64::min)
            })
            .fold(0.0, f64::max)
    }

    /// The accuracy boxes `left` hold to: the deepest a point of S in them
    /// can be (`depth`), each bounded first by the material box's faces
    /// alone, and further only where that bound is more than the most found
    /// so far (it cannot raise the most otherwise).
    fn attained(&self, left: &[&search::Node<Pieces>]) -> f64 {
        let mut order: Vec<(f64, usize)> = left.iter().enumerate().map(|(i, n)| (self.clip(n.b).map_or(0.0, |b| self.to_faces(b)), i)).collect();
        order.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));
        let mut most = 0.0f64;
        for chunk in order.chunks(search::BATCH) {
            if chunk[0].0 <= most {
                break;
            }
            let ds: Vec<f64> = chunk.par_iter().map(|&(f, i)| if f <= most { f } else { self.depth(left[i]) }).collect();
            most = ds.into_iter().fold(most, f64::max);
        }
        most
    }
}

/// The first cells and the floor of an emptiness search (prove_empty)
/// over a region `extent` across that halves down to `floor`, or to a `res`
/// given where that is finer: 1/32 of the region, coarsened (under twice)
/// to the floor times a power of two, so the halvings land on the floor
/// itself, not on cells just above it and then again just below; doubled
/// while no coarser than a res given. A region under 32 floors keeps its
/// own 1/32.
fn first_cells(extent: f64, res: Option<f64>, floor: f64) -> (f64, f64) {
    let floor = res.map_or(floor, |r| r.min(floor));
    let mut h = extent / 32.0;
    if floor > 0.0 && h > floor {
        h = floor * (h / floor).log2().ceil().exp2();
    }
    if let Some(r) = res {
        while 2.0 * h <= r {
            h *= 2.0;
        }
    }
    (h, floor)
}
