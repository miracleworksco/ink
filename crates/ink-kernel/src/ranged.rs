//! Proofs at a moment, and over ranged params (the scope beyond space,
//! `scope`). Ranged params are interval inputs (`ctx.p_spans`), so one
//! interval proof covers every param value of a scope; a piece the proof
//! cannot settle is halved where the field depends on it. These drive the
//! searches in `query` over resolution and params, and return typed
//! results: clients add names, trust labels and presentation. A claim
//! over a motion is one over ranged params, t among them (ink-lang's
//! motion).

use crate::query::{self, Verdict};
use crate::scope::{self, At};
use crate::search;
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_interval_cases, eval_interval_scope, eval_point, span_mid};
use ink_core::graph::Fingerprint;
use ink_core::interval::I;
use ink_core::tape::Tape;
use std::sync::Arc;

/// What a proof asks, and how far it may refine.
#[derive(Clone, Copy, Debug)]
pub struct ProofOptions {
    /// A cell size the search reaches (`--res`): its floor where finer than
    /// its own (`query::prove_empty`). None: the search's own.
    pub res: Option<f64>,
    /// The finest cells the search may reach where the question is exact
    /// (no tolerance): a floor set by the whole model's size, not the
    /// claim's. With a tolerance, the floor is half of it, or this where
    /// coarser.
    pub finest: f64,
    /// How many times a ranged param's span may be halved (4 times this
    /// for `value_span`).
    pub depth: u32,
    /// The question: no point of the solid deeper than this (mm; 0: none
    /// at all). See query::prove_empty.
    pub tol: f64,
}

impl ProofOptions {
    /// Hash every option into `fp` (cache::ask).
    pub fn fingerprint(&self, fp: &mut Fingerprint) {
        let ProofOptions { res, finest, depth, tol } = *self;
        fp.f64(res.unwrap_or(-1.0)).word(res.is_some() as u64).f64(finest).word(depth as u64).f64(tol);
    }

    /// The finest cells the search reaches: half the tolerance (or
    /// `finest`, where coarser), else (an exact question) `finest`.
    pub fn floor(&self) -> f64 {
        if self.tol > 0.0 { (0.5 * self.tol).max(self.finest) } else { self.finest }
    }
}

/// What a proof settled.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Status {
    /// Proven: to the tolerance asked, where one was.
    Proven,
    Disproven,
    Unknown,
}

pub use crate::query::Proof;

impl Proof {
    pub fn status(&self) -> Status {
        match self.verdict {
            Verdict::Proven => Status::Proven,
            Verdict::Disproven { .. } => Status::Disproven,
            Verdict::Unknown { .. } => Status::Unknown,
        }
    }

    pub fn disproven(&self) -> bool {
        matches!(self.verdict, Verdict::Disproven { .. })
    }

    /// Proven, and with nothing at all of the solid left (no cell needed
    /// the tolerance).
    pub fn exact(&self) -> bool {
        matches!(self.verdict, Verdict::Proven) && self.exact
    }
}

/// Prove a (domain-clipped) field's solid inside lo..hi holds no point
/// deeper than `opts.tol` (query::prove_empty), its boxes halved down to
/// `opts.floor()`, or to `opts.res` where finer, in one search: the
/// question asked is the one answered. Pure: proofs can run in parallel.
pub fn prove(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, opts: &ProofOptions) -> Proof {
    let p: Proof = crate::cache::memo_region("prove", &[tape], ctx, (lo, hi), |fp| opts.fingerprint(fp), |c| prove_now(tape, c, lo, hi, opts));
    Proof { verdict: query::at_moment(p.verdict, ctx, &[tape]), ..p }
}

/// `prove`'s key in the answer cache (None: not cached).
pub fn prove_key(tape: &Tape, ctx: &Ctx, lo: V3, hi: V3, opts: &ProofOptions) -> Option<u128> {
    crate::cache::ask_region("prove", &[tape], ctx, (lo, hi), |fp| opts.fingerprint(fp)).0
}

fn prove_now(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, opts: &ProofOptions) -> Proof {
    let tape = &query::focus(tape, ctx, lo, hi);
    query::prove_empty(tape, ctx, lo, hi, opts.res, opts.floor(), opts.tol)
}

/// A comparison against a bound: `value cmp bound`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cmp {
    Lt,
    Le,
    Gt,
    Ge,
}

impl Cmp {
    pub fn holds(self, v: f64, bound: f64) -> bool {
        match self {
            Cmp::Lt => v < bound,
            Cmp::Le => v <= bound,
            Cmp::Gt => v > bound,
            Cmp::Ge => v >= bound,
        }
    }

    /// Holds for every value in `iv` (a bound rounded outward: proven).
    pub fn holds_on(self, iv: I, bound: f64) -> bool {
        match self {
            Cmp::Lt => iv.hi < bound,
            Cmp::Le => iv.hi <= bound,
            Cmp::Gt => iv.lo > bound,
            Cmp::Ge => iv.lo >= bound,
        }
    }

    /// Fails for every value in `iv` (disproven).
    pub fn fails_on(self, iv: I, bound: f64) -> bool {
        match self {
            Cmp::Lt => iv.lo >= bound,
            Cmp::Le => iv.lo > bound,
            Cmp::Gt => iv.hi <= bound,
            Cmp::Ge => iv.hi < bound,
        }
    }
}

/// A value claim's question: `value cmp bound`, allowing a violation of no
/// more than `tol` (in the value's own units; 0: exactly). It holds where
/// the value moved `tol` toward holding compares so (`value - tol < bound`
/// for <, `value + tol >= bound` for >=), and fails where so moved it does
/// not: a violation deeper than `tol`. A comparison stays strict or not as
/// written. Enclosures are moved outward, so both answers are sound; with
/// no tolerance they are not moved at all.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Comparison {
    pub cmp: Cmp,
    pub bound: f64,
    pub tol: f64,
}

impl Comparison {
    /// `value cmp bound`, exactly.
    pub fn exact(cmp: Cmp, bound: f64) -> Comparison {
        Comparison { cmp, bound, tol: 0.0 }
    }

    /// The enclosure `iv` moved `tol` toward holding, rounded outward.
    fn eased(self, iv: I) -> I {
        if self.tol == 0.0 {
            return iv;
        }
        let by = match self.cmp {
            Cmp::Lt | Cmp::Le => -self.tol,
            Cmp::Gt | Cmp::Ge => self.tol,
        };
        ink_core::interval::add_o(iv, I::point(by))
    }

    /// Holds for every value in `iv` (a bound rounded outward).
    pub fn holds_on(self, iv: I) -> bool {
        self.cmp.holds_on(self.eased(iv), self.bound)
    }

    /// Fails, deeper than the tolerance, for every value in `iv`.
    pub fn fails_on(self, iv: I) -> bool {
        self.cmp.fails_on(self.eased(iv), self.bound)
    }
}

/// A value checked over ranged params.
#[derive(Clone, Debug, PartialEq)]
pub struct ValueSpan {
    /// Proven, Disproven or Unknown.
    pub status: Status,
    /// Its bounds over the whole scope.
    pub range: (f64, f64),
    /// A violating point of the scope and the value there.
    pub witness: Option<(At, f64)>,
    /// Pieces of the scope examined.
    pub pieces: usize,
    /// Unknown: why the proof ended (query::Ended::Finest: a piece at the
    /// refinement limit, every ranged param it reads halved `halvings`
    /// times), and the hunt for a violation after it, none found
    /// (search::Search::hunt).
    pub ended: Option<query::Ended>,
    pub hunt: Option<search::Hunt>,
    /// The most times a piece's ranged param is halved.
    pub halvings: u32,
}

/// The comparison `c` (to its tolerance) at every value of `ctx`'s ranged
/// params (`tape` read at the point `at`: a field of params there): the
/// one search over the scope (crate::search), each piece's bounds first
/// order in the ranged params (a param read twice is the same value: p - p
/// is 0), halved on the param whose narrowing to its middle narrows them
/// most, up to 4 `opts.depth` times per axis (point bounds are cheap). An
/// open piece proves nothing short of settling, so pieces tie and the
/// smaller goes first: depth first, as an emptiness search (best first on
/// how far its bounds reach past the bound refines a curve it meets breadth
/// first, all along it). A point whose outward bounds violate deeper than
/// the tolerance disproves (bounds at a piece's centre, where every ranged
/// param keeps its real meaning); a rounding across the bound does not. A
/// piece undecided at that limit (a value that meets its bound, moved by
/// the tolerance, all along a curve; or one only finer pieces decide) ends
/// the proof; the search then hunts a violation in the rest
/// (search::Search::hunt) to its cap, and is unknown where it finds none,
/// saying why each stopped. Panics without a ranged param.
pub fn value_span(tape: &Tape, ctx: &Ctx, at: V3, c: Comparison, opts: &ProofOptions) -> ValueSpan {
    assert!(!ctx.p_spans.is_empty() || ctx.t_span.is_some(), "a scope beyond space");
    let q = ValueOver { tape, at, c, splits: 4 * opts.depth };
    let p = search::Precision { finest: 0.0, splits: q.splits, maps: false };
    let mut s = search::Search::new(&q, vec![search::Node { b: at.map(I::point), ctx: Arc::new(ctx.clone()), state: () }], p);
    let mut outcome = s.finish(&q, usize::MAX);
    if matches!(&outcome, search::Outcome::Unknown(u) if !u.stopped()) {
        s.hunt();
        outcome = s.finish(&q, query::BOXES);
    }
    let iv = eval_interval_scope(tape, ctx, at.map(I::point), None);
    let (range, pieces) = ((iv.lo, iv.hi), s.tested());
    let (status, witness, ended, hunt) = match outcome {
        search::Outcome::Proven => (Status::Proven, None, None, None),
        search::Outcome::Found(w) => (Status::Disproven, Some(w), None, None),
        search::Outcome::Unknown(u) => (Status::Unknown, None, Some(u.ended), u.hunt),
    };
    ValueSpan { status, range, witness, pieces, ended, hunt, halvings: q.splits }
}

/// A value over the ranged params as a question of the one search
/// (`value_span`).
struct ValueOver<'a> {
    tape: &'a Tape,
    at: V3,
    c: Comparison,
    splits: u32,
}

impl search::Question for ValueOver<'_> {
    type State = ();
    type Witness = (At, f64);

    fn test(&self, n: &search::Node<()>) -> search::Test<(), (At, f64)> {
        let (piece, origin) = (&*n.ctx, n.b);
        let iv = eval_interval_scope(self.tape, piece, origin, None);
        // Each case of a switch it reads several times bounded apart (a
        // min's derivative read in each term of a speed's).
        if self.c.holds_on(iv) || self.c.holds_on(eval_interval_cases(self.tape, piece, origin)) {
            return search::Test::Settled;
        }
        let spans_t = scope::spans_t(self.tape, piece);
        let centre = scope::middle(piece, spans_t);
        if self.c.fails_on(eval_interval(self.tape, &centre, origin, None)) {
            return search::Test::Found((At::of(&centre), eval_point(self.tape, &centre, self.at)));
        }
        // Time, where it spans the piece, is an axis as a param is.
        let t = match piece.t_span.filter(|_| spans_t) {
            Some(s) => {
                let v = eval_interval_scope(self.tape, &Ctx { t: span_mid(s), t_span: None, ..piece.clone() }, origin, None);
                ((iv.hi - iv.lo) - (v.hi - v.lo)).max(f64::MIN_POSITIVE)
            }
            None => 0.0,
        };
        // Each param's smear: how much narrowing it to its middle narrows
        // the bounds (a solved value is derived again there), and every
        // param that can be halved may be.
        let params = scope::halvable(&scope::read_params(self.tape, piece), piece, self.splits)
            .into_iter()
            .map(|k| {
                let m = span_mid(piece.p_spans.get(k).expect("a ranged param"));
                let v = eval_interval_scope(self.tape, &piece.with_param(k, (m, m)), origin, None);
                (k, ((iv.hi - iv.lo) - (v.hi - v.lo)).max(f64::MIN_POSITIVE))
            })
            .collect();
        search::Test::Open { bound: 0.0, smear: search::Smear { t, params, ..Default::default() }, state: () }
    }

    fn seek(&self, _: &search::Node<()>) -> Option<(At, f64)> {
        None
    }
}
