//! Checking a document's claims (its `(assert ...)` forms) on a model view:
//! each claim by the cheapest method that decides it (a witness point, an
//! interval proof refined only where undecided, a measurement). Every
//! answer says whether the claim holds and how far off it could be: an
//! accuracy, checked everywhere to it, or none, an estimate (sampled,
//! simulated or solved numerically) with no guaranteed bound. Results are
//! typed; clients name and present them.

use crate::InkError;
use crate::frame::in_frame;
use crate::motion::Asked;
use crate::view::ModelView;
use crate::{Claim, Cmp, FormMetric, Measured, Placed, RunOpts};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_points};
use ink_core::graph::Var;
use ink_core::interval::{I, add_o};
use ink_core::tape::Tape;
use ink_kernel::query::{self, Verdict};
use ink_kernel::ranged::{self, Comparison, Proof, ValueSpan};
use ink_kernel::scope::At;
use std::sync::Arc;

/// Whether a claim holds. How far off that could be is the answer's
/// accuracy (Check::accuracy): guaranteed to a number, or an estimate.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Status {
    Holds,
    /// Fails: a witness violates it by more than its accuracy.
    Fails,
    Unknown,
}

impl Status {
    pub const ALL: [Status; 3] = [Status::Holds, Status::Fails, Status::Unknown];

    pub fn as_str(self) -> &'static str {
        match self {
            Status::Holds => "holds",
            Status::Fails => "fails",
            Status::Unknown => "unknown",
        }
    }

    /// A proof's status: proven, the claim holds; disproven, it fails.
    pub fn of(s: ranged::Status) -> Status {
        match s {
            ranged::Status::Proven => Status::Holds,
            ranged::Status::Disproven => Status::Fails,
            ranged::Status::Unknown => Status::Unknown,
        }
    }
}

/// What checking a claim found, by the method that decided it.
#[derive(Clone, Debug)]
pub enum Evidence {
    /// An emptiness proof of the claim's field.
    Empty(Proof),
    /// Nothing to check: the shape (or one of two) is empty, so the claim
    /// holds; the note says which.
    Vacuous(&'static str),
    /// Nothing examined: the claim needs all of a shape not proven bounded.
    Unexamined,
    /// Clearance by the one search over boxes grown by the clearance
    /// (query::clearance): the verdict, a pair of points closer than the
    /// clearance (in the second shape, in the first), and how finely it
    /// looked.
    ClearanceSearch { verdict: Verdict, pair: Option<(V3, V3)>, res: f64 },
    /// An integral, summed over the target's bodies, estimated
    /// (`estimate_parts`), for each field compared (in order), with its
    /// status, at cells of `res` (the coarsest a body settled at); over
    /// ranged params, at the values `sampled`.
    Integral { integrals: Vec<(f64, Status)>, res: f64, sampled: Option<Sampled> },
    /// A weighted mean compared with its bound: the integral of the weight,
    /// and for each field compared (in order) its mean and status, as for
    /// an integral.
    Mean { weight: f64, means: Vec<(f64, Status)>, res: f64, sampled: Option<Sampled> },
    /// A package's metric of one of its forms; `head` is the form's kind,
    /// `quantity` the metric's (its value's unit), if it has one.
    /// `undecided`: why what it measured does not decide it (what it left
    /// out would answer otherwise), and how to decide it.
    Metric { measured: Measured, head: &'static str, form: String, quantity: Option<crate::Quantity>, undecided: Option<(String, String)> },
    /// A value at the given params and t.
    Value(f64),
    /// A value at every value of the ranged params it reads.
    ValueSpan(ValueSpan),
    /// A method that answers at one value of each param (rays, a package's
    /// solve), for a claim that reads a ranged one: not run, as its answer
    /// would say nothing about the rest of the range.
    NotCovered(String),
    /// Two shapes overlapping at a moment: the emptiness proof of their
    /// intersection (a witness proves the overlap), and, where it is
    /// unknown because the intersection's proven box is thinner than twice
    /// the accuracy (a probe box, a shallow fit), that box's thinnest side:
    /// no point in it is deeper than the accuracy, so no accuracy but a
    /// smaller one can show the overlap.
    Overlaps(Proof, Option<f64>),
    /// Two shapes touching: that they overlap no deeper than the accuracy
    /// (the emptiness proof of their intersection; over a motion, the step
    /// that decided it), and the search for two points within it of each
    /// other (a clearance of the accuracy, whose witness pair is the
    /// contact).
    Touches { overlap: Box<Evidence>, gap: Box<Evidence> },
    /// A point's field value: its bounds there, rounded outward (over the
    /// ranged params' spans, where it reads them).
    Point { at: V3, bounds: (f64, f64) },
    /// Checked at steps of a motion (Check::motion), with no one step's
    /// evidence to show.
    Stepped,
    /// Not checked: it rests on a solve with no value at this moment
    /// (ModelView::refusal), whose error says why.
    Unsolved(InkError),
    /// Not checked: the call's time limit had passed (ink_kernel::stop).
    Stopped,
    /// Not checked: outside ink's range (crate::validity), saying why.
    Outside(String),
}

/// An estimate over ranged params at values of them (`estimate_scope`):
/// how many, and the values where it fails, if it does.
#[derive(Clone, Debug)]
pub struct Sampled {
    pub values: usize,
    pub failing: Option<Vec<(u32, f64)>>,
}

/// A claim's answer.
#[derive(Clone, Debug)]
pub struct Check {
    pub status: Status,
    /// How far off the answer could be, in the claim's units: it holds or
    /// fails guaranteed to this (0: exactly), under the claim's own
    /// predicate (docs/IR.md Results). None: not guaranteed: unknown, or
    /// an estimate (sampled inputs, a solve's numbers, rays, a package's
    /// solve, steps of a motion without a margin). Inside the checker Some
    /// marks a method that guarantees its answer; `settle` writes the
    /// number.
    pub accuracy: Option<f64>,
    /// The accuracy the claim was asked to, where its method guarantees
    /// its answer (None: an estimate): an unknown was undecided at it.
    pub asked: Option<f64>,
    /// What `accuracy` and `asked` are in, where known (None: the model's
    /// own units): Doc::accuracy_quantity, or the claim's (AssertDef).
    pub quantity: Option<crate::units::Quantity>,
    pub evidence: Evidence,
    /// For a shape not proven bounded, what the method settled in the
    /// window, where that is not an answer about all of it.
    pub in_window: Option<Status>,
    /// Only the window was examined, and the claim needs all of the shape
    /// (it is not proven bounded: intersect it with a containing shape).
    pub unbounded: bool,
    /// The ranged params the claim reads, and their spans: it is proven or
    /// disproven over every value of them.
    pub ranged: Vec<(u32, (f64, f64))>,
    /// How the claim was checked, where that differs from how it reads.
    pub note: Option<&'static str>,
    /// The inexact inputs the claim reads (sampled fields, solved values
    /// not certified): its answer is measured on them.
    pub inexact: crate::Inexact,
    /// Over a motion, how it was checked (crate::motion).
    pub motion: Option<crate::motion::Stepped>,
    /// The one moment it was asked at, though its question reads t
    /// (`reads_t`): no span covered it and the run asked for none, so it
    /// holds or fails there only. Recorded as it is asked (`asserted`).
    pub one_moment: Option<f64>,
    /// An emptiness claim's at a moment: how it held (`held`).
    pub held: Option<Held>,
}

/// How an emptiness claim (empty, disjoint, subset) held.
#[derive(Clone, Debug, PartialEq)]
pub enum Held {
    /// No cell needed the accuracy's allowance: nothing violates it.
    Outright,
    /// Some cells held only to the accuracy. `near`: the deepest violating
    /// point found among NEAR_SAMPLES^3 points through each (depth 0 and no
    /// point: none found), a measurement, with the moment and ranged
    /// params' values it was found at (`At`: none ranged, its params
    /// empty).
    ToAccuracy { near: Option<(f64, Option<(V3, At)>)> },
}

/// Points measured along each side of a loose cell (Held::candidates).
const NEAR_SAMPLES: usize = 8;

/// A proof's deepest points below 0 tried for confirmation (Held::candidates).
const CONFIRM_TRIES: usize = 16;

impl Held {
    /// Where a proof's loose cells may hold a violating point: each
    /// cell's NEAR_SAMPLES^3 points below 0 on `tape` (the question's, in
    /// its frame, simplified over the cell and its piece of the scope), at
    /// the middle of that piece and where each ranged param (t among them
    /// over a motion) is pushed to the end of its span that lowers the
    /// field at the cell's centre (scope::pushed: a thin overlap at the
    /// end of a motion), the CONFIRM_TRIES deepest, with where in the
    /// scope.
    pub fn scoped(p: &Proof, tape: &Tape, ctx: &Ctx) -> Vec<(V3, At)> {
        let (mut values, mut below) = (Vec::new(), Vec::new());
        for l in &p.loose {
            let [lo, hi] = l.cell;
            let over = l.scope(ctx);
            let mut ch = vec![ink_core::tape::BOTH; tape.len()];
            eval_interval(tape, &over, std::array::from_fn(|k| I::new(lo[k], hi[k])), Some(&mut ch));
            let sub = tape.simplify(&ch);
            let points: Vec<V3> = (0..NEAR_SAMPLES.pow(3))
                .map(|i| {
                    let at = [i % NEAR_SAMPLES, i / NEAR_SAMPLES % NEAR_SAMPLES, i / NEAR_SAMPLES / NEAR_SAMPLES];
                    std::array::from_fn(|k| lo[k] + (hi[k] - lo[k]) * (at[k] as f64 + 0.5) / NEAR_SAMPLES as f64)
                })
                .collect();
            let mid = ink_kernel::scope::middle(&over, over.t_span.is_some());
            let centre = std::array::from_fn(|k| 0.5 * (lo[k] + hi[k]));
            let pushed = ink_kernel::scope::pushed(&[&sub], &over, centre).map(|x| x.0);
            for at in std::iter::once(mid).chain(pushed) {
                eval_points(&sub, &at, &points, &mut values);
                below.extend(values.iter().zip(&points).filter(|x| *x.0 < 0.0).map(|(f, c)| (*f, *c, At::of(&at))));
            }
        }
        below.sort_by(|a, b| a.0.total_cmp(&b.0));
        below.into_iter().take(CONFIRM_TRIES).map(|x| (x.1, x.2)).collect()
    }

    /// `scoped`'s points (interfere's pairs, which read no scope).
    pub fn candidates(p: &Proof, tape: &Tape, ctx: &Ctx) -> Vec<V3> {
        Held::scoped(p, tape, ctx).into_iter().map(|x| x.0).collect()
    }

    /// `confirmed_at`, each point at ctx's moment.
    pub fn confirmed(exact: bool, points: impl IntoIterator<Item = V3>, world: &Tape, ctx: &Ctx) -> Held {
        Held::confirmed_at(exact, points.into_iter().map(|c| (c, At { t: ctx.t, params: Vec::new() })), world, ctx)
    }

    /// How a claim whose proof was `exact` (or not) held, its near point
    /// the deepest of `points` (in the world, each at its moment and
    /// ranged params' values: ctx's point there, Ctx::point) where
    /// `world`, the claim's field as placed (a pair's, with its zones), is
    /// below 0 by bounds rounded outward (a float value may round away
    /// what the real one is), as deep as the bound.
    pub fn confirmed_at(exact: bool, points: impl IntoIterator<Item = (V3, At)>, world: &Tape, ctx: &Ctx) -> Held {
        if exact {
            return Held::Outright;
        }
        let bound = |c: V3, at: &At| {
            let mut x = Ctx { t: at.t, ..ctx.clone() };
            for &(k, v) in &at.params {
                x.params[k as usize] = v;
            }
            eval_interval(world, &x.point(), c.map(I::point), None).hi
        };
        let deepest = points.into_iter().map(|(c, at)| (bound(c, &at), (c, at))).filter(|x| x.0 < 0.0).min_by(|a, b| a.0.total_cmp(&b.0));
        let depth = deepest.as_ref().map_or(0.0, |d| -d.0);
        Held::ToAccuracy { near: Some(Held::near(depth, deepest.map(|d| d.1))) }
    }

    /// A near point at least `depth` deep, rounded to a µm: none at 0.
    pub fn near(depth: f64, at: Option<(V3, At)>) -> (f64, Option<(V3, At)>) {
        let depth = (depth * 1e3).round().max(0.0) / 1e3;
        (depth, at.filter(|_| depth > 0.0))
    }

    /// Measured deeper than 0.
    pub fn near_depth(&self) -> bool {
        matches!(self, Held::ToAccuracy { near: Some((d, _)) } if *d > 0.0)
    }
}

impl Check {
    /// The same answer, an estimate: no accuracy is guaranteed, for it or
    /// any status in it (what it settled in the window, an integral's).
    pub(crate) fn estimate(self) -> Check {
        Check { accuracy: None, ..self }
    }

    /// Holds, to a guaranteed accuracy.
    pub fn holds_guaranteed(&self) -> bool {
        self.status == Status::Holds && self.accuracy.is_some()
    }

    /// The answer asked to `tol`: where its method guarantees it, that is
    /// the accuracy asked, and the answer's accuracy where it holds or
    /// fails (an unknown guarantees nothing).
    pub(crate) fn settle(self, tol: f64) -> Check {
        let asked = self.accuracy.map(|_| tol);
        Check { asked, accuracy: asked.filter(|_| self.status != Status::Unknown), ..self }
    }

    /// What an unknown still holds to: the accuracy its search attained
    /// (an emptiness: no point of what violates it deeper than that; a
    /// clearance d: no pair closer than d less it), and whether only
    /// the window was examined (then it holds there, and says nothing of
    /// the rest of the shapes). None unless the answer is the claim's own
    /// search at the view's moment, by a method that guarantees it: a
    /// step of a motion's or a piece of a sweep's search answers another
    /// question (a margin's, a box of the params'), one reading inexact
    /// inputs is an estimate, and a claim whose accuracy is not a length (a
    /// range, decided exactly) has none to hold to.
    pub fn attained(&self) -> Option<(f64, bool)> {
        // A length, as the claim's own accuracy must be (a range's takes
        // none).
        let own = self.status == Status::Unknown && self.asked.is_some() && self.quantity == Some(crate::units::Quantity::Length) && self.motion.is_none();
        match &self.evidence {
            Evidence::ClearanceSearch { verdict: Verdict::Unknown { attained: Some(a), .. }, .. } if own => Some((*a, self.unbounded)),
            Evidence::Empty(Proof { verdict: ink_kernel::query::Verdict::Unknown { attained: Some(a), .. }, .. }) if own => Some((*a, self.unbounded)),
            _ => None,
        }
    }
}

/// How finely claims are checked.
#[derive(Clone, Copy, Debug, Default)]
pub struct ClaimOptions {
    /// Cell size for geometry claims (None: automatic).
    pub res: Option<f64>,
    /// The accuracy the claim is checked to, in its units (None: its
    /// default, `claim_accuracy`).
    pub accuracy: Option<f64>,
    /// Asked beside other searches (interfere's moving pairs): over a
    /// motion, t ranged is not stopped early to leave time to others
    /// (motion::ranged), which would stop those beside it too.
    pub beside: bool,
}

/// One assert's answer.
#[derive(Clone, Debug)]
pub struct AssertResult {
    pub name: String,
    /// The claim as written.
    pub claim: String,
    /// Where it is written: `file:line`.
    pub source: String,
    /// The cell size it was checked at, if one was given, and whether its
    /// own `(res r)` gave it (the finer of that and the run's).
    pub res: Option<(f64, ResBy)>,
    pub check: Check,
    /// Wall-clock milliseconds it took to check.
    pub ms: f64,
    /// Answered from the cache: what it rests on was read from a cache and
    /// nothing was computed (ink_kernel::cache::HITS and MISSES).
    pub cached: bool,
    /// AssertDef::terms where the value is (its witness, or the given
    /// params and t); empty where there is no one point.
    pub terms: Vec<(String, f64)>,
}

/// What gave a claim the res it was checked at.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ResBy {
    /// The claim's own `(res r)`.
    Claim,
    /// The run's (`ClaimOptions::res`, the CLI's `--res`).
    Run,
}

/// The res a claim is checked at: the finer of its own and the run's. A
/// claim's res says how fine it must look to be decided; a run's says how
/// fine every claim looks at least, so neither makes the other coarser.
pub fn claim_res(own: Option<f64>, run: Option<f64>) -> Option<(f64, ResBy)> {
    match (own, run) {
        (Some(c), Some(r)) if c < r => Some((c, ResBy::Claim)),
        (Some(c), None) => Some((c, ResBy::Claim)),
        (_, Some(r)) => Some((r, ResBy::Run)),
        (None, None) => None,
    }
}

/// Every assert's answer.
#[derive(Clone, Debug)]
pub struct AssertReport {
    pub asserts: Vec<AssertResult>,
    /// With `own`: how many asserts written in imported files were left
    /// out, unchecked. None: none were left out for being imported.
    pub imported_not_run: Option<usize>,
}

impl AssertReport {
    /// How many asserts came out each way, in `Status::ALL` order.
    pub fn counts(&self) -> Vec<(Status, usize)> {
        Status::ALL.iter().map(|&s| (s, self.asserts.iter().filter(|a| a.check.status == s).count())).filter(|c| c.1 > 0).collect()
    }

    /// Every claim holds, guaranteed or estimated.
    pub fn all_hold(&self) -> bool {
        self.asserts.iter().all(|a| a.check.status == Status::Holds)
    }

    /// Every claim holds to a guaranteed accuracy.
    pub fn all_hold_guaranteed(&self) -> bool {
        self.asserts.iter().all(|a| a.check.holds_guaranteed())
    }

    /// The answer's note on its claims asked at one moment
    /// (Check::one_moment), from what each recorded.
    pub fn note(&self) -> Option<String> {
        one_moment_note(&self.asserts.iter().filter_map(|a| Some((format!("claim {}", a.name), a.check.one_moment?))).collect::<Vec<_>>())
    }
}

pub use ink_kernel::progress::progress;

/// Check the document's asserts (only the one named `only`, if given; only
/// those written in the file loaded, not in its imports, if `own`), at the
/// view's params and time. Each is reported as it is decided (`progress`).
pub fn run_asserts(v: &mut ModelView, only: Option<&str>, own: bool, opts: ClaimOptions) -> Result<AssertReport, InkError> {
    run_asserts_with(v, only, own, opts, |_| None)
}

/// An answer's note on what it asked at one moment though the question
/// reads t (claims, pairs: each named, with the moment it recorded): said
/// once, for the answer. None: nothing was.
pub fn one_moment_note(asked: &[(String, f64)]) -> Option<String> {
    let mut at: Vec<(f64, Vec<&str>)> = Vec::new();
    for (what, t) in asked {
        match at.iter_mut().find(|m| m.0 == *t) {
            Some(m) => m.1.push(what),
            None => at.push((*t, vec![what])),
        }
    }
    let each: Vec<String> = at.iter().map(|(t, ws)| format!("{} at t = {t}", ws.join(", "))).collect();
    (!each.is_empty()).then(|| format!("answered at one moment only, though each reads t and no (t-span) covers it: {}. Add (t-span a b) or use --t-span", each.join("; ")))
}

/// `run_asserts`, each answer reported with its row of the caller's answer,
/// `row` (so a run stopped before its answer gives the rows done as the
/// whole answer would).
pub fn run_asserts_with(
    v: &mut ModelView,
    only: Option<&str>,
    own: bool,
    opts: ClaimOptions,
    row: impl Fn(&AssertResult) -> Option<serde_json::Value>,
) -> Result<AssertReport, InkError> {
    v.scope_question()?;
    if let Some(n) = only {
        match v.doc.asserts.iter().find(|x| x.name == n) {
            None => {
                let hint = crate::choices(n, "Its asserts", v.doc.asserts.iter().map(|x| x.name.as_str()));
                return Err(InkError::msg(format!("no assert named '{n}'. {hint}")));
            }
            Some(x) if own && !v.doc.is_own(x) => {
                return Err(InkError::msg(format!("assert '{n}' is imported (from {}:{}), and --own checks only this file's", x.file, x.line)));
            }
            Some(_) => {}
        }
    }
    let imported_not_run = own.then(|| v.doc.asserts.iter().filter(|x| !v.doc.is_own(x)).count());
    type Listed = (String, String, String, Claim, Option<u32>, (Option<f64>, Option<f64>, Option<crate::units::Quantity>), Vec<(String, u32)>, Vec<usize>, std::rc::Rc<str>);
    let list: Vec<Listed> = v
        .doc
        .asserts
        .iter()
        .filter(|x| only.is_none_or(|n| n == x.name) && (!own || v.doc.is_own(x)))
        .map(|x| {
            let source = format!("{}:{}", x.file, x.line);
            (x.name.clone(), x.text.clone(), source, x.claim.clone(), x.bound, (x.res, x.accuracy, x.quantity), x.terms.clone(), x.fixed.clone(), x.file.clone())
        })
        .collect();
    let mut asserts = Vec::new();
    progress(serde_json::json!({"total": list.len(), "of": "asserts", "kept": true}));
    for (name, claim, source, c, bound, (own_res, tol, quantity), terms, fixed, file) in list {
        // Over the motion of the file it is written in, imported or not
        // (ModelView::motion_of).
        v.motion_of(&file).map_err(|e| InkError { message: format!("assert '{name}' ({source}): {}", e.message), ..e })?;
        // The bound at the view's params (solved ones included, wherever
        // the solve is in the file).
        let c = v.doc.claim_at(&c, bound, &v.ctx);
        let roots: Vec<u32> = v.doc.claim_roots(&c).into_iter().chain(bound).collect();
        // What it sweeps, said as it starts: the ranged params it reads.
        let mut running = serde_json::json!({"running": name});
        let ranged = ranged_reads(v, &roots);
        if !ranged.is_empty() {
            running["ranged"] = ranged.iter().map(|&(k, s)| (v.doc.params[k as usize].name.clone(), serde_json::json!([s.0, s.1]))).collect::<serde_json::Map<_, _>>().into();
        }
        progress(running);
        let res = claim_res(own_res, opts.res);
        let start = std::time::Instant::now();
        let count = || {
            use std::sync::atomic::Ordering::Relaxed;
            (ink_kernel::cache::MISSES.load(Relaxed), ink_kernel::cache::HITS.load(Relaxed))
        };
        let before = count();
        let (check, terms) = asserted(v, c, bound, res.map(|r| r.0), (tol, quantity), &fixed, &terms, opts)
            .map_err(|e| InkError { message: format!("assert '{name}' ({source}): {}", e.message), ..e })?;
        let ms = start.elapsed().as_secs_f64() * 1e3;
        // From the cache: read from it, nothing recomputed. One not begun
        // (past the call's time limit) is not cached either.
        let after = count();
        let cached = after.0 == before.0 && after.1 > before.1 && !matches!(check.evidence, Evidence::Stopped);
        let a = AssertResult { name, claim, source, res, check, ms, cached, terms };
        let mut done = serde_json::json!({"done": a.name, "status": a.check.status.as_str(), "ms": ms});
        if let Some(r) = row(&a) {
            done["row"] = r;
        }
        let one_moment = a.check.one_moment.is_some();
        progress(done);
        asserts.push(a);
        // The answer's note as it grows, so a run stopped before its
        // answer (ink mcp) gives it with the rows done.
        if one_moment {
            let so_far = AssertReport { asserts: std::mem::take(&mut asserts), imported_not_run: None };
            progress(serde_json::json!({"answer": {"note": so_far.note()}}));
            asserts = so_far.asserts;
        }
    }
    Ok(AssertReport { asserts, imported_not_run })
}

/// One claim asked as `ink assert` asks it (`c` at the view's params,
/// Doc::claim_at; `bound` its bound's field): what it reads provided while
/// it is checked, and no longer; unknown where it rests on a solve with no
/// value here, or was read at load with a ranged param as one number
/// (`fixed`); a bound reading a ranged param moved into the claim; at cells
/// of `res` (claim_res) and the claim's accuracy and quantity, if it gives
/// them. With the values of `terms` where its value is. Every claim asked
/// of a design is asked here: an assert's, an optimization's requirement.
#[allow(clippy::too_many_arguments)]
pub fn asserted(
    v: &mut ModelView,
    mut c: Claim,
    bound: Option<u32>,
    res: Option<f64>,
    (tol, quantity): (Option<f64>, Option<crate::units::Quantity>),
    fixed: &[usize],
    terms: &[(String, u32)],
    opts: ClaimOptions,
) -> Result<(Check, Vec<(String, f64)>), InkError> {
    let roots: Vec<u32> = v.doc.claim_roots(&c).into_iter().chain(bound).collect();
    // Resting on a solve with no value here (one depending on t, over a
    // span; one that does not solve at this moment), it is not decided:
    // unknown, with the solve's error; the rest are.
    // What is provided only where read (an effort) is provided while
    // the claim reading it is checked, and no longer.
    let wanted = v.doc.wanted.clone();
    let refused = v.want(&roots).err();
    // Over a motion, a claim that moves with t is solved at each step.
    let unsolved = refused.or_else(|| v.refusal(&roots)).filter(|_| !v.motion.is_some_and(|span| moves(v, &roots, span)));
    // Read at load with a ranged param as one number (AssertDef::fixed),
    // it is not decided over the range.
    let fixed: Vec<usize> = fixed.iter().copied().filter(|&k| v.ctx.p_spans.get(k as u32).is_some()).collect();
    let unsolved = unsolved.or_else(|| {
        (!fixed.is_empty()).then(|| InkError {
            hint: Some(v.doc.at_one_value(&fixed, &v.ctx.p_spans)),
            ..InkError::msg(format!(
                "it reads ranged {} as one number on loading (its relation's kind reads its parts' shapes). No proof covers the range",
                fixed.iter().map(|&k| v.doc.params[k].name.as_str()).collect::<Vec<_>>().join(", ")
            ))
        })
    });
    // A bound that reads a ranged param varies with it: a value claim
    // compares the difference with 0 (one expression, so a param in
    // both keeps its dependency, and a witness fails at its own
    // values); other bounds are numbers by definition. Either way the
    // bound is then in the claim, not beside it.
    let mut note = None;
    let mut beside = bound;
    if let Some(b) = bound.filter(|&b| !ranged_reads(v, &[b]).is_empty()) {
        match c {
            Claim::Value(f, cmp, _) => {
                c = Claim::Value(v.doc.g.sub(f, b), cmp, 0.0);
                beside = None;
                note = Some("its bound reads a ranged param, so the value checked is the expression minus the bound, against 0");
            }
            Claim::Mean { parts, frames, cmp, .. } => {
                let parts = parts.into_iter().map(|(s, w, hs)| (s, w, hs.into_iter().map(|h| v.doc.g.sub(h, b)).collect())).collect();
                c = Claim::Mean { parts, frames, cmp, bound: 0.0 };
                beside = None;
                note = Some("its bound reads a ranged param, so the mean checked is that of the field minus the bound, against 0");
            }
            _ => {
                let ks: Vec<usize> = ranged_reads(v, &[b]).iter().map(|(k, _)| *k as usize).collect();
                let names: Vec<&str> = ks.iter().map(|&k| v.doc.params[k].name.as_str()).collect();
                v.doc.wanted = wanted;
                return Err(InkError {
                    hint: Some(format!("a bound is one number. Write the varying part into the claim, e.g. (value (- e b) > 0), or {}", v.doc.at_one_value(&ks, &v.ctx.p_spans))),
                    ..InkError::msg(format!("its bound reads ranged param(s) {}", names.join(", ")))
                });
            }
        }
    }
    // A run's accuracy is a length, for the claims about shapes.
    let tol = tol.or(opts.accuracy.filter(|_| c.accuracy_is_length()));
    let mut check = match unsolved {
        // Past the call's time limit, the rest are not begun.
        _ if ink_kernel::stop::reached() => {
            Check { quantity: v.doc.accuracy_quantity(&c), ranged: ranged_reads(v, &roots), ..checked(Status::Unknown, Evidence::Stopped).estimate() }
        }
        Some(e) => Check { quantity: v.doc.accuracy_quantity(&c), ..checked(Status::Unknown, Evidence::Unsolved(e)).estimate() },
        None => check(v, c, beside, ClaimOptions { res, accuracy: tol, ..opts })?,
    };
    check.note = check.note.or(note);
    // Asked at one moment though it reads t, no moment asked: recorded
    // with the answer, at the moment it was asked.
    if check.motion.is_none() {
        check.one_moment = one_moment_of(v, &roots);
    }
    check.quantity = check.quantity.or(quantity);
    if let Some(b) = bound {
        for r in ranged_reads(v, &[b]) {
            if !check.ranged.contains(&r) {
                check.ranged.push(r);
            }
        }
        check.ranged.sort_by_key(|r| r.0);
    }
    // The claim's terms where its value is: at its witness, else at
    // the given params and t.
    let at = match &check.evidence {
        Evidence::ValueSpan(ValueSpan { witness: Some((at, _)), .. }) => {
            let mut c = v.ctx.clone();
            // Over a motion with t ranged, its moment is its param t's.
            c.t = check.motion.as_ref().and_then(|m| m.motion.t).unwrap_or(at.t);
            for &(k, x) in &at.params {
                c.params[k as usize] = x;
            }
            Some(c.point())
        }
        _ if check.motion.as_ref().is_some_and(|m| m.motion.t.is_none()) => None,
        Evidence::Value(_) => {
            let mut c = v.ctx.clone();
            c.t = check.motion.as_ref().and_then(|m| m.motion.t).unwrap_or(c.t);
            Some(c)
        }
        _ => None,
    };
    let terms = at.map(|c| terms_at(v, terms, c)).unwrap_or_default();
    v.doc.wanted = wanted;
    Ok((check, terms))
}

/// Check one claim (`bound`: its bound, as a field, if it reads params). A
/// claim that reads inexact inputs (Doc::inexact: sampled fields, a solved
/// value not certified as a root) is an estimate: whatever proves it proves
/// it about those samples or numbers. So is one checked by a measuring
/// method (rays, a package's solve), whatever it reads.
///
/// Every claim is checked here, at a moment, over ranged params or over a
/// motion, and here alone is its accuracy settled: it is asked once, to
/// its own (a shape no thicker than twice it holds to it: nothing violates
/// it by more). The answer's accuracy is the one asked, where its method
/// guarantees it. The view is scoped as every claim's question first
/// (ModelView::scope_question).
pub fn check(v: &mut ModelView, claim: Claim, bound: Option<u32>, opts: ClaimOptions) -> Result<Check, InkError> {
    ask(v, claim, bound, opts)?.wait(v)
}

/// `check`, its searches over a motion left to the caller to run
/// (motion::Asked).
pub fn ask(v: &mut ModelView, claim: Claim, bound: Option<u32>, opts: ClaimOptions) -> Result<Asked<Check>, InkError> {
    v.scope_question()?;
    ask_in(v, claim, bound, opts)
}

/// `ask` in the view's scope as it is: a pointwise view's, at one value
/// of each param (`distance` bisects so, at the values it is asked at).
fn ask_in(v: &mut ModelView, claim: Claim, bound: Option<u32>, opts: ClaimOptions) -> Result<Asked<Check>, InkError> {
    let roots: Vec<u32> = v.doc.claim_roots(&claim).into_iter().chain(bound).collect();
    let tol = claim_accuracy(v, &claim, opts.accuracy);
    let inexact = v.doc.inexact(&roots, &v.certified());
    let motion = v.motion.filter(|&span| moves(v, &roots, span));
    let opts = ClaimOptions { accuracy: Some(tol), ..opts };
    if let Some(why) = crate::validity::outside(v, &claim, &roots, tol, opts.res) {
        let c = Check { accuracy: None, ..checked(Status::Unknown, Evidence::Outside(why)) };
        return Ok(Asked::Done(finished(v, &claim, &roots, c, tol, inexact)));
    }
    let c = match motion {
        // Over a motion, how it held is said where t ranged decided it
        // (its near point at a moment), never of a moment measured.
        Some(span) => crate::motion::over_span(v, claim.clone(), bound, span, opts)?.then(v, |_, c| {
            let ranged = c.motion.as_ref().is_some_and(|m| m.ranged);
            Ok(Asked::Done(Check { held: c.held.clone().filter(|_| ranged), ..c }))
        })?,
        None => Asked::Done(decide(v, claim.clone(), opts)?),
    };
    c.then(v, move |v, c| Ok(Asked::Done(finished(v, &claim, &roots, c, tol, inexact))))
}

/// `ask_in`'s answer, as it is decided: `answer`, with the ranged params
/// it reads, t among them where it was ranged over a motion
/// (motion::ranged).
fn finished(v: &ModelView, claim: &Claim, roots: &[u32], c: Check, tol: f64, inexact: crate::Inexact) -> Check {
    let mut ranged = ranged_reads(v, roots);
    for &r in &c.ranged {
        if !ranged.contains(&r) {
            ranged.push(r);
        }
    }
    answer(c, tol, inexact, v.doc.accuracy_quantity(claim), ranged)
}

/// A claim's answer over a motion decided after it was first answered
/// (motion::Upgrade), finished as `ask` finishes one.
pub(crate) fn finish_later(v: &ModelView, claim: &Claim, bound: Option<u32>, c: Check, opts: ClaimOptions) -> Check {
    let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(bound).collect();
    let tol = claim_accuracy(v, claim, opts.accuracy);
    let inexact = v.doc.inexact(&roots, &v.certified());
    finished(v, claim, &roots, c, tol, inexact)
}

/// A claim's answer as it is decided (`decide`, or a still pair of `ink
/// interfere`'s proven beside the others: interfere::prove_pairs),
/// finished: settled at the accuracy asked, `tol`; an estimate where it
/// reads inexact inputs; in the claim's `quantity`, with the ranged params
/// it reads, at a moment or over a motion alike (a witness gives each
/// one's value).
pub(crate) fn answer(c: Check, tol: f64, inexact: crate::Inexact, quantity: Option<crate::units::Quantity>, ranged: Vec<(u32, (f64, f64))>) -> Check {
    let c = if inexact.any() { c.estimate() } else { c };
    Check { inexact, quantity, ranged, ..c.settle(tol) }
}

/// The accuracy a claim is checked to: the one given, else its default:
/// for a claim about shapes the model's (a length), for a comparison
/// (a value, a volume, an integral, a mean, a metric, a thickness) 0, and
/// a point or a range is decided exactly (Claim::accuracy).
pub(crate) fn claim_accuracy(v: &ModelView, claim: &Claim, given: Option<f64>) -> f64 {
    match claim.accuracy() {
        crate::Accuracy::Length => given.unwrap_or_else(|| v.accuracy()),
        crate::Accuracy::Value => given.unwrap_or(0.0),
        crate::Accuracy::Exact => 0.0,
    }
}

/// `check` at the view's moment, before the inexact inputs it reads are
/// taken into account.
pub(crate) fn decide(v: &mut ModelView, claim: Claim, opts: ClaimOptions) -> Result<Check, InkError> {
    let res = opts.res;
    let roots = v.doc.claim_roots(&claim);
    // The claim's accuracy: a shape claim's a length, a comparison's in
    // its own units (`claim_accuracy`).
    let tol = claim_accuracy(v, &claim, opts.accuracy);
    if let Some(c) = estimated(v, &claim, &roots, tol) {
        return Ok(Check { ranged: ranged_reads(v, &roots), ..c });
    }
    let mut c = match claim {
        // Values in a region are compared exactly: their accuracy is 0.
        Claim::Empty(f) => shallow(v, f, tol, res),
        Claim::Range(f) => Check { held: None, ..shallow(v, f, tol, res) },
        // A clearance d: no point of one within d of the other, exactly,
        // with cells down to half the accuracy, or a sixteenth of d where
        // finer (a clearance sets its own scale): a gap within about that
        // of d is unknown, never proven, so the accuracy can only leave an
        // answer undecided. A clearance of 0 is the separation question,
        // disjointness, to the accuracy.
        Claim::Clearance(a, b, d) if d > 0.0 => clearance(v, a, b, d, res, (if tol > 0.0 { 0.5 * tol } else { v.finest() }).min(d / 16.0)),
        Claim::Disjoint(a, b) | Claim::Clearance(a, b, _) => {
            let f = v.doc.g.max(a, b);
            shallow(v, f, tol, res)
        }
        Claim::Subset(a, b) => {
            let nb = v.doc.g.neg(b);
            let f = v.doc.g.max(a, nb);
            shallow(v, f, tol, res)
        }
        Claim::Touches(a, b) => touches(v, a, b, tol, res),
        Claim::Point { f, at, inside } => point(v, f, at, inside),
        Claim::Integral { parts, frames, cmp, bound } => {
            let c = Comparison { cmp, bound, tol };
            let parts: Vec<Part> = parts.iter().enumerate().map(|(k, (s, gs))| (*s, gs.iter().map(|&g| Some(g)).collect(), frames.get(k).copied().flatten())).collect();
            let each = |i: &[f64]| -> Vec<Status> { i.iter().map(|&x| point_status(c, x)).collect() };
            let (s, sampled) = estimate_scope(v, &parts, res, |i| all_of(each(i).into_iter()));
            let integrals: Vec<(f64, Status)> = s.sums.iter().copied().zip(each(&s.sums)).collect();
            Check { unbounded: !s.whole, ..checked(s.status, Evidence::Integral { integrals, res: s.res, sampled }).estimate() }
        }
        Claim::Mean { parts, frames, cmp, bound } => mean(v, &parts, &frames, Comparison { cmp, bound, tol }, res),
        Claim::Metric { metric, cmp, value: bound, excluding } => {
            // A package solves at one point of the scope: where its form
            // reads a ranged axis (itself, or through a form it reads, as a
            // thermal study reads its heat study), its answer says nothing
            // about the range.
            if scope_reads(v, &roots) {
                let names = &v.doc.params;
                let num = |x: f64| (x * 1e9).round() / 1e9;
                let read: Vec<String> = ranged_reads(v, &roots).iter().map(|&(k, (a, b))| format!("{} ({}..{})", names[k as usize].name, num(a), num(b))).collect();
                let e = InkError {
                    hint: Some(at_one_value(v, &roots)),
                    ..InkError::msg(format!("a package's metric is solved at one value of each param, but this form reads ranged {}. No solve covers the range", read.join(", ")))
                };
                checked(Status::Unknown, Evidence::Unsolved(e))
            } else if metric.shape.is_some_and(|f| !v.whole(f)) {
                // Of a shape seen only in the window, the window's faces
                // would be its walls.
                Check { unbounded: true, ..checked(Status::Unknown, Evidence::Unexamined) }
            } else {
                metric_claim(v, &metric, Comparison { cmp, bound, tol }, if metric.shape.is_some() { res } else { None }, &excluding)?
            }
        }
        Claim::Value(f, cmp, bound) if ranged_reads(v, &[f]).is_empty() => {
            let x = v.doc.scalar(f, &v.ctx);
            let c = Comparison { cmp, bound, tol };
            // One that reads t, at a moment, by its bounds there rounded
            // outward: its float value may round away what the real one is
            // ((1e16 + t) - 1e16 is t, 0 in floats). Undecided by them, a
            // float value holding is an estimate (a value on its bound, an
            // effort clamped to its limit: 0 >= 0); one failing is no
            // witness, and says nothing.
            if v.doc.g.depends_on(f, |n| n == Var::T) {
                let iv = ink_core::eval::eval_interval_scope(&Tape::new(&v.doc.g, f), &v.ctx, [I::point(0.0); 3], None);
                return Ok(if c.holds_on(iv) {
                    checked(Status::Holds, Evidence::Value(x))
                } else if c.fails_on(iv) {
                    checked(Status::Fails, Evidence::Value(x))
                } else if point_status(c, x) == Status::Holds {
                    Check { accuracy: None, ..checked(Status::Holds, Evidence::Value(x)) }
                } else {
                    checked(Status::Unknown, Evidence::Value(x))
                });
            }
            checked(point_status(c, x), Evidence::Value(x))
        }
        Claim::Value(f, cmp, bound) => {
            // For every value of each ranged param: interval bounds over
            // pieces, halving undecided ones; a point violating it deeper
            // than the accuracy disproves.
            let tape = Tape::new(&v.doc.g, f);
            let result = ranged::value_span(&tape, &v.ctx, [0.0; 3], Comparison { cmp, bound, tol }, &v.proof_options(res, 0.0));
            checked(Status::of(result.status), Evidence::ValueSpan(result))
        }
        Claim::Overlaps(a, b) => {
            // Overlap is proven by a witness deeper than the accuracy in
            // both: exactly when the shallowness of the intersection is
            // disproven. Disproven only by an intersection proven empty
            // outright: one shown no deeper than the accuracy may still
            // overlap.
            let f = v.doc.g.max(a, b);
            let (p, whole) = v.prove_with(f, &v.proof_options(res, tol));
            let apart = p.exact();
            let mut c = in_window(Status::of(p.status()), whole, Evidence::Overlaps(p, None));
            c.status = match c.status {
                Status::Fails => Status::Holds,
                Status::Holds if apart => Status::Fails,
                _ => Status::Unknown,
            };
            // Depth is distance from outside: every point of a box thinner
            // than twice the accuracy is within the accuracy of leaving it.
            // The intersection lies in its own proven box and in each
            // shape's.
            if c.status == Status::Unknown
                && whole
                && tol > 0.0
                && let (Some(i), Some(ea), Some(eb)) = (v.bounds(f), v.bounds(a), v.bounds(b))
                && let side = (0..3).map(|k| i.1[k].min(ea.1[k]).min(eb.1[k]) - i.0[k].max(ea.0[k]).max(eb.0[k])).fold(f64::INFINITY, f64::min)
                && side < 2.0 * tol
                && let Evidence::Overlaps(_, thin) = &mut c.evidence
            {
                *thin = Some(side.max(0.0));
            }
            // Over ranged params a witness is at one value of them, and
            // says nothing of the rest: proven apart at every value fails.
            if c.status == Status::Holds && !ranged_reads(v, &[a, b]).is_empty() {
                (c.status, c.note) = (Status::Unknown, Some(AT_ONE_VALUE));
            }
            c
        }
    };
    c.ranged = ranged_reads(v, &roots);
    c.accuracy = c.accuracy.map(|_| tol);
    Ok(c)
}

/// A claim reading params set as estimates that gauge their own error
/// (ParamProvider::error): a value claim at a point holds where it holds
/// for every value in their bands (interval bounds; its value then the
/// bands' worst), fails where it fails for every one, and is unknown
/// between, saying so; any other claim reading them is unknown. None: it
/// reads none.
fn estimated(v: &mut ModelView, claim: &Claim, roots: &[u32], tol: f64) -> Option<Check> {
    let mut band: Vec<(usize, f64)> = Vec::new();
    for i in v.doc.solves_read(roots) {
        let s = &v.doc.solves[i];
        let Some(p) = v.doc.ext[s.ext].form.provider(&v.doc) else { continue };
        if let Some(w) = p.error(&v.ctx.params) {
            band.extend(s.owns.iter().copied().zip(w));
        }
    }
    if band.is_empty() {
        return None;
    }
    let read = v.doc.params_read(roots);
    let mut at = None;
    if let Claim::Value(f, cmp, bound) = *claim
        && ranged_reads(v, &[f]).is_empty()
        && band.iter().all(|b| read.contains(&b.0) && b.1.is_finite())
    {
        let c = Comparison { cmp, bound, tol };
        let spans: Vec<(u32, (f64, f64))> = band
            .iter()
            .map(|&(k, w)| {
                let b = add_o(I::point(v.ctx.params[k]), I::new(-w, w));
                (k as u32, (b.lo, b.hi))
            })
            .collect();
        let over = Ctx { p_spans: v.ctx.p_spans.clone().with_derived(&spans), ..v.ctx.clone() };
        let tape = Tape::new(&v.doc.g, f);
        let iv = ink_core::eval::eval_interval(&tape, &over, [I::point(0.0); 3], None);
        let x = ink_core::eval::eval_point(&tape, &v.ctx, [0.0; 3]);
        if c.holds_on(iv) {
            let worst = if matches!(cmp, Cmp::Lt | Cmp::Le) { iv.hi } else { iv.lo };
            return Some(Check { accuracy: None, ..checked(Status::Holds, Evidence::Value(worst)) });
        }
        if c.fails_on(iv) || x.is_nan() {
            return Some(Check { accuracy: None, ..checked(Status::Fails, Evidence::Value(x)) });
        }
        at = Some(x);
    }
    let bands: Vec<String> = band.iter().map(|&(k, w)| format!("{} {:.6} within {w:.6}", v.doc.params[k].name, v.ctx.params[k])).collect();
    let why = format!(
        "it reads an estimate with a band either way ({}), from its own convergence measure, not a bound (its note says how). {}",
        bands.join(", "),
        match at {
            Some(x) => format!("Its value {x:.6} neither holds nor fails over all of that band"),
            None => "Only a value claim at fixed params is decided over such a band".into(),
        }
    );
    Some(checked(Status::Unknown, Evidence::Unsolved(InkError::msg(why))))
}

/// A value at a point (no ranged param) against its comparison: it holds,
/// or fails (an undefined value, NaN, too); a value within rounding of
/// the bound moved by the accuracy is unknown.
fn point_status(c: Comparison, x: f64) -> Status {
    let p = I::point(x);
    if c.holds_on(p) {
        Status::Holds
    } else if c.fails_on(p) || x.is_nan() {
        Status::Fails
    } else {
        Status::Unknown
    }
}

/// The side of the cube about a point whose sign its bounds there leave
/// undecided (`point`): every point of it in the shape, or none, puts the
/// point there.
const POINT_CUBE: f64 = 0.005;

/// How far below 0 the field must be shown over that cube for the point to
/// be inside (not on a face of no thickness inside the material, where the
/// field is 0, as where two pieces abut).
const POINT_MARGIN: f64 = 1e-9;

/// The solid of field `f` holds no point deeper than `tol` (an emptiness
/// claim, or an overlap's, a subset's, a clearance of 0's).
fn shallow(v: &mut ModelView, f: u32, tol: f64, res: Option<f64>) -> Check {
    let (p, whole) = v.prove(f, res, tol);
    let tape = v.clipped(f);
    let points = if p.exact { Vec::new() } else { Held::scoped(&p, &tape, &v.ctx) };
    let held = Some(Held::confirmed_at(p.exact, points, &tape, &v.ctx));
    Check { held, ..emptied(p, whole) }
}

/// An emptiness claim's answer from its proof, all of its field examined
/// or not (`whole`: else it holds in the window only).
pub(crate) fn emptied(p: Proof, whole: bool) -> Check {
    in_window(Status::of(p.status()), whole, Evidence::Empty(p))
}

/// A touches claim's note where its parts come no nearer each other than
/// its accuracy.
pub(crate) const NOT_TOUCHING: &str = "they come no nearer each other than the accuracy: not touching";

/// `(touches a b)`: they overlap no deeper than the accuracy `tol`
/// (`shallow`), and come within it of each other at every value of the
/// ranged params (`contact`).
fn touches(v: &mut ModelView, a: u32, b: u32, tol: f64, res: Option<f64>) -> Check {
    let f = v.doc.g.max(a, b);
    let over = shallow(v, f, tol, res);
    let overlap = Box::new(over.evidence.clone());
    if over.status != Status::Holds {
        return Check { evidence: Evidence::Touches { overlap, gap: Box::new(Evidence::Unexamined) }, ..over };
    }
    let gap = contact(v, a, b, tol, res);
    let mut c = Check { note: gap.note, ..checked(gap.status, Evidence::Touches { overlap, gap: Box::new(gap.evidence) }) };
    if c.status == Status::Fails {
        c.note = Some(NOT_TOUCHING);
    }
    if gap.unbounded {
        (c.unbounded, c.in_window) = (true, gap.in_window);
    }
    c
}

/// A point of b within d of one of a, at every value of the ranged params
/// a and b read: proven by a pair of points closer than d (the exact
/// clearance search's witness, which disproves a clearance of d), and
/// disproven where that clearance is proven. A pair found at one value of
/// ranged params says nothing of the others: there it is unknown.
pub(crate) fn contact(v: &mut ModelView, a: u32, b: u32, d: f64, res: Option<f64>) -> Check {
    let finest = if d > 0.0 { 0.5 * d } else { v.finest() };
    let read: Vec<u32> = ranged_reads(v, &[a, b]).into_iter().map(|r| r.0).collect();
    let flip = |c: Check| {
        let status = match c.status {
            Status::Fails => Status::Holds,
            Status::Holds => Status::Fails,
            s => s,
        };
        Check { status, ..c }
    };
    let mut c = flip(clearance(v, a, b, d, res, finest));
    // Over ranged params the pair found is at one value of them, and says
    // nothing of the rest: apart by d at every value fails.
    if c.status == Status::Holds && !read.is_empty() {
        (c.status, c.note) = (Status::Unknown, Some(AT_ONE_VALUE));
    }
    c
}

/// What an overlap or a contact over ranged params says found at one value
/// of them.
pub(crate) const AT_ONE_VALUE: &str = "over ranged params, a point or pair found shows it at one value of them, not every value. Only apart at every value fails";

/// A point inside the shape `f` (`inside`) or not: its field's bounds
/// there, rounded outward (real numbers), against 0; over the ranged
/// params it reads, for every value of them (value_span, halving pieces).
/// The coordinates are numbers, at their values (constants and params);
/// coordinates that read a ranged param or t move the point with it: the
/// field read there is then checked as a value of the params (at a moment,
/// t one value: its outward bounds there decide, never the field at the
/// coordinates' floats, which may round the point elsewhere). A coordinate not
/// finite somewhere in the scope (one overflowed) is no point to be inside
/// or outside of: unknown.
fn point(v: &mut ModelView, f: u32, at: [u32; 3], inside: bool) -> Check {
    let cmp = if inside { Cmp::Lt } else { Cmp::Ge };
    let p = at.map(|n| v.doc.scalar(n, &v.ctx));
    // Each coordinate's bounds over the scope, before the field reads it:
    // a moving one's values are not p's, and one the field does not read
    // still says where the point is.
    let finite = at.iter().all(|&n| {
        let c = ink_core::eval::eval_interval_scope(&Tape::new(&v.doc.g, n), &v.ctx, [I::point(0.0); 3], None);
        c.lo.is_finite() && c.hi.is_finite()
    });
    if !finite {
        return checked(Status::Unknown, Evidence::Point { at: p, bounds: (f64::NEG_INFINITY, f64::INFINITY) });
    }
    let reads_t = |v: &ModelView, ns: &[u32]| ns.iter().any(|&n| v.doc.g.depends_on(n, |x| x == Var::T));
    // The coordinates that move, the field reads as them; the others it
    // reads as x, y and z at p, spatial values with real-number bounds
    // (substituted, a number's arithmetic would be folded in floats: at
    // y = 0.5, (1e16 + y) - 1e16 is 0).
    let moves: Vec<(Var, u32)> = [Var::X, Var::Y, Var::Z].into_iter().zip(at).filter(|&(_, n)| !ranged_reads(v, &[n]).is_empty() || reads_t(v, &[n])).collect();
    let moving = !moves.is_empty();
    // Searched over a scope: a ranged param read, or t's span where it reads t.
    let ranged = !ranged_reads(v, &at).is_empty() || !ranged_reads(v, &[f]).is_empty() || (v.ctx.t_span.is_some() && (reads_t(v, &at) || reads_t(v, &[f])));
    let tape = if moving {
        let w = v.doc.g.subst(f, &moves);
        Tape::new(&v.doc.g, w)
    } else {
        Tape::new(&v.doc.g, f)
    };
    // The field reads the coordinates that do not move at p, where the
    // point is at the params' values (the evidence's).
    let iv = ink_core::eval::eval_interval_scope(&tape, &v.ctx, p.map(I::point), None);
    let evidence = Evidence::Point { at: p, bounds: (iv.lo, iv.hi) };
    let c = Comparison::exact(cmp, 0.0);
    if c.holds_on(iv) {
        return checked(Status::Holds, evidence);
    }
    if !ranged && c.fails_on(iv) {
        return checked(Status::Fails, evidence);
    }
    // Bounds at a point straddle 0 where a branch cut of the field's
    // formula runs through it (a thread's angle at its root), though its
    // value there is one number: every point of a small cube about it at
    // least POINT_MARGIN inside the shape (or not inside it), shown exactly
    // by the search, cells either side of the cut, at every value of the
    // ranged params, puts it there. Its corners are rounded outward, so the
    // point is strictly inside it however large its coordinates (at 1e16,
    // p - 0.0025 rounds to p, and a cube of no width holds nothing to
    // show: both sides would be proven).
    if !moving {
        let lo = p.map(|x| (x - 0.5 * POINT_CUBE).next_down());
        let hi = p.map(|x| (x + 0.5 * POINT_CUBE).next_up());
        let g = &mut v.doc.g;
        let side = if inside {
            let m = g.constant(POINT_MARGIN);
            let nf = g.neg(f);
            g.sub(nf, m)
        } else {
            f
        };
        let cube = v.doc.box_field(lo, hi);
        let field = v.doc.g.max(cube, side);
        let (proof, whole) = v.prove(field, None, 0.0);
        if whole && proof.status() == ranged::Status::Proven {
            return checked(Status::Holds, evidence);
        }
    }
    if !ranged {
        return checked(Status::Unknown, evidence);
    }
    let result = ranged::value_span(&tape, &v.ctx, p, c, &v.proof_options(None, 0.0));
    checked(Status::of(result.status), Evidence::ValueSpan(result))
}

/// Several comparisons, all claimed: it fails where one does, holds where
/// all do.
fn all_of(each: impl Iterator<Item = Status>) -> Status {
    each.fold(Status::Holds, |all, s| match (all, s) {
        (Status::Fails, _) | (_, Status::Fails) => Status::Fails,
        (Status::Holds, Status::Holds) => Status::Holds,
        _ => Status::Unknown,
    })
}

/// A part of an integral claim: its shape, its integrands (None meaning 1),
/// and its body's own frame, where it is a placed body.
type Part = (u32, Vec<Option<u32>>, Option<Placed>);

/// Integrals summed over parts, estimated, and the claim's status on them.
struct Summed {
    /// One per integrand: the sum over the parts.
    sums: Vec<f64>,
    /// The coarsest cells a part stopped at (0: none had material).
    res: f64,
    /// Every part is proven bounded (else the sums are of the rest).
    whole: bool,
    status: Status,
}

/// Sums over `parts` of each integrand's integral over the part's shape,
/// estimated at the view's params and moment, each part halved as
/// query::refine halves one, from a 100th of its own box (at `res`, once,
/// where given), all stopping once every part is halved and `done` takes
/// the sums with their spreads; the sums, the coarsest cells a part stopped at, and whether
/// every part is proven bounded (else the sums are of the rest, and
/// `done` takes none). A placed body is integrated in its own frame
/// (frame::in_frame), times 1 / |det R|: copies of a body are one
/// question there, answered once (the cache).
fn estimate_parts(v: &mut ModelView, parts: &[Part], res: Option<f64>, done: impl Fn(&[f64], &[f64]) -> bool) -> (Vec<f64>, f64, bool) {
    let n = parts.first().map_or(0, |p| p.1.len());
    // Each part as integrated (its frame's shape and integrands), alike
    // parts once, with how many and their scale.
    let mut alike: Vec<((u32, Vec<Option<u32>>), f64)> = Vec::new();
    for (shape, gs, frame) in parts {
        let local = frame.and_then(|f| in_frame(v, &f, gs)).filter(|f| v.whole(f.shape));
        let (key, scale) = match local {
            Some(f) => ((f.shape, f.gs), f.inv_det.mid()),
            None => ((*shape, gs.clone()), 1.0),
        };
        match alike.iter_mut().find(|a| a.0 == key) {
            Some(a) => a.1 += scale,
            None => alike.push((key, scale)),
        }
    }
    let whole = alike.iter().all(|((s, _), _)| v.whole(*s));
    let mut each = Vec::new();
    for ((shape, gs), scale) in alike {
        let Some((lo, hi, _)) = v.bounds(shape).filter(|_| v.whole(shape)) else { continue };
        let tapes: Vec<Option<Tape>> = gs.iter().map(|g| g.map(|g| Tape::new(&v.doc.g, g))).collect();
        each.push((v.clipped(shape), tapes, lo, hi, scale));
    }
    let ctx = &v.ctx;
    let at = |(tape, tapes, lo, hi, _): &(Arc<Tape>, Vec<Option<Tape>>, V3, V3, f64), h: f64| {
        let gs: Vec<Option<&Tape>> = tapes.iter().map(Option::as_ref).collect();
        let m = query::measured(tape, &gs, false, ctx, *lo, *hi, Some(h));
        query::Sums { sums: m.sums, cells: m.cells }
    };
    let mut found: Vec<query::Measured> =
        each.iter().map(|p| query::Measured { settled: res.is_some(), ..query::Measured::first(res.unwrap_or(query::max_extent(p.2, p.3) / 100.0), |h| at(p, h)) }).collect();
    loop {
        let (mut sums, mut spread) = (vec![0.0; n], vec![0.0; n]);
        for (m, p) in found.iter().zip(&each) {
            for k in 0..n {
                sums[k] += p.4 * m.sums[k];
                spread[k] += p.4 * m.spread.get(k).unwrap_or(&0.0);
            }
        }
        if !found.iter().any(query::Measured::halves) || (whole && found.iter().all(|m| m.coarse.is_some()) && done(&sums, &spread)) {
            return (sums, found.iter().map(|m| m.res).fold(0.0, f64::max), whole);
        }
        for (m, p) in found.iter_mut().zip(&each) {
            if m.halves() {
                *m = m.halved(|h| at(p, h));
            }
        }
    }
}

/// Values of each ranged param an estimate over them is taken at, where
/// it is the only one read (`estimate_scope`).
const PARAM_SAMPLES: usize = 17;

/// `estimate_parts` over the claim's scope, with the claim's status on the
/// sums (`status`): at the view's params where its fields read no ranged
/// param; else at a grid of their values, each param's span at evenly
/// spaced values, ends included (PARAM_SAMPLES of one; of n, about
/// PARAM_SAMPLES in all, at least both ends of each). It fails at the
/// first value where it fails (the sums there), holds where it holds at
/// every one (the sums at the last), else is unknown: an estimate, which a
/// violation between the values escapes. What the fields read is what is
/// sampled: of a part in its frame, its shape and integrands there, and
/// its placement.
fn estimate_scope(v: &mut ModelView, parts: &[Part], res: Option<f64>, status: impl Fn(&[f64]) -> Status) -> (Summed, Option<Sampled>) {
    let mut fields: Vec<u32> = Vec::new();
    for (s, gs, f) in parts {
        match f.and_then(|f| Some((in_frame(v, &f, gs)?, f))).filter(|(i, _)| v.whole(i.shape)) {
            // In its frame: its shape and integrands there, and its
            // placement where its determinant varies (a stretch).
            Some((i, f)) => {
                if i.inv_det.hi - i.inv_det.lo > 1e-9 * i.inv_det.hi {
                    fields.extend(f.at.r.iter().flatten());
                }
                fields.extend(std::iter::once(i.shape).chain(i.gs.into_iter().flatten()));
            }
            None => fields.extend(std::iter::once(*s).chain(gs.iter().flatten().copied())),
        }
    }
    let read = ranged_reads(v, &fields);
    let summed = |v: &mut ModelView| {
        // Settled where the status stays with the sums moved by their
        // spreads, each way (the first and the rest apart too: a mean's
        // weight and moments).
        let stays = |x: &[f64], d: &[f64]| {
            let st = status(x);
            let moved = |a: f64, b: f64| x.iter().zip(d).enumerate().map(|(k, (x, d))| x + if k == 0 { a } else { b } * d).collect::<Vec<f64>>();
            st != Status::Unknown && [(1.0, 1.0), (-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0)].iter().all(|&(a, b)| status(&moved(a, b)) == st)
        };
        let (sums, res, whole) = estimate_parts(v, parts, res, stays);
        let status = if whole { status(&sums) } else { Status::Unknown };
        Summed { sums, res, whole, status }
    };
    if read.is_empty() {
        return (summed(v), None);
    }
    let per = ((PARAM_SAMPLES as f64).powf(1.0 / read.len() as f64).floor() as usize).max(2);
    let values = per.pow(read.len() as u32);
    let mut last = None;
    for k in 0..values {
        let at: Vec<(u32, f64)> = read.iter().enumerate().map(|(i, &(p, (lo, hi)))| (p, lo + (hi - lo) * (k / per.pow(i as u32) % per) as f64 / (per - 1) as f64)).collect();
        let ctx = at.iter().fold(v.ctx.clone(), |c, &(p, x)| c.with_param(p, (x, x)));
        let s = v.at_point(ctx, |v| summed(v));
        let fails = s.status == Status::Fails;
        let unknown = s.status == Status::Unknown;
        last = Some(s);
        if fails {
            return (last.unwrap(), Some(Sampled { values: k + 1, failing: Some(at) }));
        }
        if unknown {
            break;
        }
    }
    (last.expect("a value"), Some(Sampled { values, failing: None }))
}

/// `(mean target h (weight w) cmp v)`: the integral of w h over the parts
/// over that of w, for each h (a vec), each compared with v, estimated
/// (`estimate_scope`). Where the weight's integral is not positive there
/// is no mean: unknown.
fn mean(v: &mut ModelView, parts: &[(u32, u32, Vec<u32>)], frames: &[Option<Placed>], c: Comparison, res: Option<f64>) -> Check {
    let split: Vec<Part> = parts
        .iter()
        .enumerate()
        .map(|(k, (s, w, hs))| {
            let moments = hs.iter().map(|&h| Some(v.doc.g.mul(*w, h)));
            (*s, std::iter::once(Some(*w)).chain(moments).collect(), frames.get(k).copied().flatten())
        })
        .collect();
    let each = |i: &[f64]| -> Vec<(f64, Status)> { i[1..].iter().map(|a| a / i[0]).map(|m| (m, if i[0] > 0.0 { point_status(c, m) } else { Status::Unknown })).collect() };
    let (s, sampled) = estimate_scope(v, &split, res, |i| all_of(each(i).into_iter().map(|m| m.1)));
    let means = each(&s.sums);
    let mut c = checked(s.status, Evidence::Mean { weight: s.sums[0], means, res: s.res, sampled }).estimate();
    if s.sums[0] <= 0.0 {
        c.note = Some("the weight's integral over the target is not positive (a target with no material has no mean): no mean to compare");
    }
    c.unbounded = !s.whole;
    c
}

/// An answer by a method that guarantees it (`check` writes its
/// accuracy); `Check::estimate` makes one an estimate.
pub(crate) fn checked(status: Status, evidence: Evidence) -> Check {
    Check {
        status,
        accuracy: Some(0.0),
        asked: None,
        quantity: None,
        evidence,
        in_window: None,
        unbounded: false,
        ranged: Vec::new(),
        note: None,
        inexact: Default::default(),
        motion: None,
        one_moment: None,
        held: None,
    }
}

/// A claim's terms at `c`: what a solve on t sets taken as it sets it
/// at c's moment (an actuator's effort at a witness), the rest at c.
fn terms_at(v: &ModelView, terms: &[(String, u32)], mut c: Ctx) -> Vec<(String, f64)> {
    let nodes: Vec<u32> = terms.iter().map(|t| t.1).collect();
    let moving: Vec<usize> = v.doc.solves_read(&nodes).into_iter().filter(|&i| v.doc.solves[i].on_t).collect();
    if !moving.is_empty()
        && c.t != v.ctx.t
        && let Ok((at, _)) = v.doc.ctx_at(&v.overrides, c.t)
    {
        for k in moving.iter().flat_map(|&i| v.doc.solves[i].owns.clone()) {
            c.params[k] = at.params[k];
        }
    }
    terms.iter().map(|(n, f)| (n.clone(), v.doc.scalar(*f, &c))).collect()
}

/// The ranged params (in the view's scope) these fields read, with their spans.
fn ranged_reads(v: &ModelView, nodes: &[u32]) -> Vec<(u32, (f64, f64))> {
    v.ranged_reads(nodes)
}

/// Doc::at_one_value for the ranged params these fields read.
fn at_one_value(v: &ModelView, nodes: &[u32]) -> String {
    v.doc.at_one_value(&ranged_reads(v, nodes).into_iter().map(|(k, _)| k as usize).collect::<Vec<_>>(), &v.ctx.p_spans)
}

/// Whether these fields read a param the view ranges over. A method
/// evaluated at one point cannot cover such a scope.
fn scope_reads(v: &ModelView, nodes: &[u32]) -> bool {
    !ranged_reads(v, nodes).is_empty()
}

/// The one moment a question reading `roots` is asked at, though it reads
/// t (`reads_t`): no span covers it and the run asked for none
/// (ModelView::unspanned_moment). The one rule a claim's and an interfere
/// pair's answer record (Check::one_moment, Pair::one_moment), and `ink
/// check` warns by.
pub fn one_moment_of(v: &ModelView, roots: &[u32]) -> Option<f64> {
    v.unspanned_moment().filter(|_| reads_t(v, roots)).map(|_| v.ctx.t)
}

/// Whether the model reads t anywhere (a t in its graph, a solve on t):
/// where it does not, no question does.
pub fn model_reads_t(v: &ModelView) -> bool {
    v.doc.g.nodes.iter().any(|n| matches!(n, ink_core::graph::Node::Var(Var::T))) || v.doc.solves.iter().any(|s| s.on_t)
}

/// The claims `ink assert` asks (every one, imported or not) that it would
/// ask at one moment though they read t, in a view scoped as assert scopes
/// it: each over its file's span (ModelView::motion_of), at the view's
/// params, its bound included, by `one_moment_of`. Named, for `ink check`.
pub fn claims_at_one_moment(v: &mut ModelView) -> Result<Vec<String>, InkError> {
    let list: Vec<(String, Claim, Option<u32>, std::rc::Rc<str>)> = v.doc.asserts.iter().map(|a| (a.name.clone(), a.claim.clone(), a.bound, a.file.clone())).collect();
    let mut out = Vec::new();
    for (name, claim, bound, file) in list {
        v.motion_of(&file)?;
        let c = v.doc.claim_at(&claim, bound, &v.ctx);
        let roots: Vec<u32> = v.doc.claim_roots(&c).into_iter().chain(bound).collect();
        if one_moment_of(v, &roots).is_some() {
            out.push(format!("claim {name}"));
        }
    }
    Ok(out)
}

/// Time may enter through a provider's param as well as through the graph.
pub fn reads_t(v: &ModelView, nodes: &[u32]) -> bool {
    nodes.iter().any(|&f| v.doc.g.depends_on(f, |x| x == Var::T)) || v.doc.solves_read(nodes).iter().any(|&i| v.doc.solves[i].on_t)
}

/// Whether how an `interfere` pair on `node` held goes unsaid: it reads
/// ranged params, or t over a span.
pub(crate) fn held_unsaid(v: &ModelView, node: u32) -> bool {
    !v.ranged_reads(&[node]).is_empty() || (v.motion.is_some() || v.ctx.t_span.is_some()) && reads_t(v, &[node])
}

/// Whether fields `nodes` change with t over `span`: through a solve that
/// reads t, or by a rate in t not 0 somewhere (bounded over all of space,
/// the span and every value of the ranged params). One reading t only
/// through a term that is 0 at these params (a motion times a param set
/// to 0) does not: a claim on it is the claim at any moment.
pub fn moves(v: &ModelView, nodes: &[u32], span: (f64, f64)) -> bool {
    v.doc.solves_read(nodes).iter().any(|&i| v.doc.solves[i].on_t)
        || nodes.iter().filter(|&&f| v.doc.g.depends_on(f, |x| x == Var::T)).any(|&f| ink_kernel::motion::rate(&Tape::new(&v.doc.g, f), &v.ctx, [I::ENTIRE; 3], span) != Some(0.0))
}

/// A proof's answer about a shape that may not be proven bounded: a proof
/// in the window is unknown for all of it (a witness there still counts).
pub(crate) fn in_window(status: Status, whole: bool, evidence: Evidence) -> Check {
    let mut c = checked(status, evidence);
    if !whole && status != Status::Fails {
        if status == Status::Holds {
            (c.status, c.in_window) = (Status::Unknown, Some(status));
        }
        c.unbounded = true;
    }
    c
}

/// No point of b within d of a, exactly, with boxes no finer than
/// `finest` (or a res given where finer): the one search
/// (query::clearance), covering one shape proven bounded with boxes and
/// growing the other by d around each, either way round, since what it
/// needs, the grown shape's slope within d of the covered, may be bounded
/// for one and not the other (a polygon's sign has no slope bound at its
/// surface; a field promised 1.41 grows by more than one exact near its
/// faces). Only the covered shape within d of the other's box matters.
/// With neither proven bounded, all of space is tested as one box (two
/// parallel slabs settle there); else only the window is searched, and
/// only a witness pair (a violation) is an answer.
fn clearance(v: &mut ModelView, pa: u32, pb: u32, d: f64, res: Option<f64>, finest: f64) -> Check {
    // A shape proven empty. One not proven bounded and empty in the window
    // says nothing about the rest of it: a bounded partner gives the
    // search its region.
    if [pa, pb].into_iter().any(|f| v.whole(f) && v.bounds(f).is_none()) {
        return checked(Status::Holds, Evidence::Vacuous("one of the shapes is empty"));
    }
    let finest = res.map_or(finest, |r| r.min(finest));
    let (ta, tb) = (Arc::new(Tape::new(&v.doc.g, pa)), Arc::new(Tape::new(&v.doc.g, pb)));
    let whole = v.whole(pa) || v.whole(pb);
    if !whole {
        let all = query::Way { grown: ta.clone(), covered: tb.clone(), lo: [f64::NEG_INFINITY; 3], hi: [f64::INFINITY; 3], swapped: false };
        if let (verdict @ Verdict::Proven, pair, h) = query::clearance(&[all], &v.ctx, d, finest, query::BOXES) {
            return checked(Status::Holds, Evidence::ClearanceSearch { verdict, pair, res: h });
        }
    }
    let Some(ways) = clearance_ways(v, (pa, &ta), (pb, &tb), d) else {
        // No point of the covered shape is within d of the other.
        return checked(Status::Holds, Evidence::Vacuous("the shapes' boxes are farther apart than the clearance"));
    };
    if ways.is_empty() {
        return Check { unbounded: true, ..checked(Status::Unknown, Evidence::Unexamined) };
    }
    let (verdict, pair, h) = query::clearance(&ways, &v.ctx, d, finest, query::BOXES);
    let status = match verdict {
        Verdict::Proven => Status::Holds,
        Verdict::Disproven { .. } => Status::Fails,
        Verdict::Unknown { .. } => Status::Unknown,
    };
    let mut c = checked(status, Evidence::ClearanceSearch { verdict, pair, res: h });
    if !whole && status != Status::Fails {
        if status == Status::Holds {
            (c.status, c.in_window) = (Status::Unknown, Some(Status::Holds));
        }
        c.unbounded = true;
    }
    c
}

/// The ways round a clearance d's search (`clearance`, query::Way): each
/// shape proven bounded covered with boxes (with neither, the second, in
/// the window), within d of the other's box where that is proven bounded.
/// None: the boxes are farther apart than d.
fn clearance_ways(v: &mut ModelView, (pa, ta): (u32, &Arc<Tape>), (pb, tb): (u32, &Arc<Tape>), d: f64) -> Option<Vec<query::Way>> {
    let whole = v.whole(pa) || v.whole(pb);
    let mut ways = Vec::new();
    for (cov, grown, swapped) in [(pb, pa, false), (pa, pb, true)] {
        if !(v.whole(cov) || !whole && !swapped) {
            continue;
        }
        // Neither proven bounded, and the covered one empty in the window:
        // nothing to search.
        let Some((mut lo, mut hi, _)) = v.bounds(cov) else { continue };
        if v.whole(grown) {
            let (glo, ghi, _) = v.bounds(grown).expect("not proven empty");
            for k in 0..3 {
                lo[k] = lo[k].max((glo[k] - d).next_down());
                hi[k] = hi[k].min((ghi[k] + d).next_up());
            }
        }
        if (0..3).any(|k| lo[k] > hi[k]) {
            return None;
        }
        let grown = if swapped { tb.clone() } else { ta.clone() };
        ways.push(query::Way { grown, covered: v.clipped(cov), lo, hi, swapped });
    }
    Some(ways)
}

/// The minimum distance between two solids as a bracket (`distance`).
pub struct Distance {
    /// Proven: every point of each at least this from the other.
    pub lower: f64,
    /// A pair of points confirmed inside a and b, and their distance
    /// rounded up: an upper bound. None: no pair was confirmed.
    pub upper: Option<(f64, V3, V3)>,
    /// The search ended short of the accuracy (query::Ended::Finest or
    /// Limit) with d, the clearance it was proving, this: the bracket is as
    /// tight as the search goes at its budget.
    pub undecided: Option<f64>,
    /// Boxes the search tested.
    pub boxes: usize,
    /// The call's time limit stopped the search (ink_kernel::stop).
    pub stopped: bool,
}

/// The minimum distance between a and b, bracketed to the model's
/// accuracy by one search (query::min_distance): a clearance's, covering
/// each shape proven bounded with boxes (`clearance_ways`), its d shrinking
/// to each pair closer than it that the search finds, confirmed inside
/// both, less the accuracy, from the farthest the targets' boxes are apart
/// (where it finds none, again from a pair proven inside them).
/// The lower bound is what it proved: nothing where neither shape is proven
/// bounded (only the window was searched) or the targets read inputs that
/// are not exact (samples, solved values not certified), on which nothing
/// is proven. Each bound is sound whatever the fields: a field that is not
/// a distance (a solid cut by a plane) never stands in for one.
pub fn distance(v: &mut ModelView, fa: u32, fb: u32, opts: ClaimOptions) -> Result<Distance, InkError> {
    let (Some((alo, ahi, _)), Some((blo, bhi, _))) = (v.bounds(fa), v.bounds(fb)) else { return Err(InkError::msg("a target is empty")) };
    let far = ink_core::interval::norm(std::array::from_fn(|k| (ahi[k] - blo[k]).max(bhi[k] - alo[k]).next_up()));
    let (ta, tb) = (Arc::new(Tape::new(&v.doc.g, fa)), Arc::new(Tape::new(&v.doc.g, fb)));
    let tol = v.accuracy();
    let finest = opts.res.map_or(v.finest(), |r| r.min(v.finest())).min(if tol > 0.0 { 0.5 * tol } else { f64::INFINITY });
    let ways = clearance_ways(v, (fa, &ta), (fb, &tb), far).unwrap_or_default();
    let mut found = query::min_distance(&ways, &v.ctx, far, None, tol, finest, query::BOXES);
    // Undecided with no pair found: the search's seeks walk along a
    // field's slope, which a field flat away from its solid (a union with
    // a constant) does not have there. Again from a pair of points proven
    // inside them (query::inside_pair). Proven with none, a target has no
    // point inside (its field is never below 0), and there is no pair to
    // find.
    if found.upper.is_none()
        && matches!(found.ended, Some(query::Ended::Finest | query::Ended::Limit))
        && let Some(seed) = query::inside_pair(&ta, &tb, &v.ctx, (alo, ahi), (blo, bhi), finest)
    {
        let ways = clearance_ways(v, (fa, &ta), (fb, &tb), seed.0).unwrap_or_default();
        let again = query::min_distance(&ways, &v.ctx, far, Some(seed), tol, finest, query::BOXES);
        found = query::MinDistance { lower: again.lower.max(found.lower), tested: found.tested + again.tested, ..again };
    }
    let proven = (v.whole(fa) || v.whole(fb)) && !v.doc.inexact(&[fa, fb], &v.certified()).any();
    let stopped = found.ended == Some(query::Ended::Stopped);
    let undecided = matches!(found.ended, Some(query::Ended::Finest | query::Ended::Limit)).then(|| found.upper.map_or(far, |u| u.0) - tol);
    Ok(Distance { lower: if proven { found.lower } else { 0.0 }, upper: found.upper, undecided, boxes: found.tested, stopped })
}

/// The grid error of value `x` read on cells of h, from the same read on
/// cells of 2h (`x2`) and 4h (`x4`, if it has one): Roache's grid
/// convergence index. Where the three converge monotonically, at the
/// observed order p = log2((x4 - x2) / (x2 - x)), 1.25 |x - x2| / (2^p - 1);
/// otherwise (two grids, or three that do not converge: a feature the
/// coarser grids cannot see) 3 |x - x2|, the index of an assumed first
/// order. Never under |x - x2|, the optimizer's rule
/// (ink-design `finer`). With how it was found.
pub fn grid_error(x: f64, x2: f64, x4: Option<f64>) -> (f64, &'static str) {
    let d1 = x - x2;
    match x4.map(|x4| x2 - x4) {
        Some(d2) if d1 * d2 > 0.0 && d1.abs() < d2.abs() => {
            let p = (d2 / d1).log2();
            (d1.abs().max(1.25 * d1.abs() / (2f64.powf(p) - 1.0)), "three grids converging")
        }
        Some(d2) if d1 == 0.0 && d2 == 0.0 => (0.0, "three grids alike"),
        Some(_) => (3.0 * d1.abs(), "three grids not converging"),
        None => (3.0 * d1.abs(), "two grids"),
    }
}

/// A package's metric against its comparison: an estimate (its accuracy
/// moves the comparison, and guarantees nothing; exact where the metric is
/// proven, Measured::proven), and unknown where the
/// metric has no value (a failed solve; nothing to measure over). Forms
/// keep their own (res ...): the geometry's res is not theirs (cells a
/// hundredth of a millimetre would not fit a study); `res` overrides (a
/// shape's metric: the claim's).
fn metric_claim(v: &mut ModelView, m: &FormMetric, c: Comparison, res: Option<f64>, excluding: &[&str]) -> Result<Check, InkError> {
    let mut measured = v.measure(m, RunOpts { res, tol: 1e-8 })?;
    let head = v.doc.ext_form(Some(m.package), &m.form).map_or("form", |e| e.head);
    let mut status = match measured.value {
        Ok(x) => point_status(c, x),
        Err(_) => Status::Unknown,
    };
    // What the metric leaves out, where counting it would answer
    // otherwise: the claim is undecided by what was measured, and says
    // what was left out, where, and how to decide it; unless it leaves
    // that out on purpose, when the answer stands and says so.
    let num = |x: f64| if x.is_finite() { format!("{}", (x * 1e4).round() / 1e4) } else { "unbounded".into() };
    let mut undecided: Option<(Vec<String>, Vec<String>)> = None;
    // A proven value's claim is exact, unless it leaves out on purpose
    // what would change it.
    let mut exact = measured.proven;
    if measured.value.is_ok() {
        let measured_status = status;
        for ex in &measured.excluded {
            let counted = point_status(c, ex.value);
            if counted == measured_status {
                continue;
            }
            let at = ex.at.map_or(String::new(), |p| format!(" at [{:.3}, {:.3}, {:.3}]", p[0], p[1], p[2]));
            if excluding.contains(&ex.key) {
                exact = false;
                measured.warnings.push(format!(
                    "this claim {} leaving out {} on purpose (excluding {}). Counted, that ({}{at}) would make it {}",
                    measured_status.as_str(),
                    ex.what,
                    ex.key,
                    num(ex.value),
                    counted.as_str()
                ));
            } else {
                status = Status::Unknown;
                let u = undecided.get_or_insert_with(Default::default);
                u.0.push(format!(
                    "it {} on a value that leaves out {}. Counted, that ({}{at}) would make it {}",
                    measured_status.as_str(),
                    ex.what,
                    num(ex.value),
                    counted.as_str()
                ));
                u.1.push(format!("{}, or leave it out on purpose with (excluding {})", ex.decide, ex.key));
            }
        }
    }
    // A grid's answer holds only where its grid saw the part, and its grid
    // error, estimated from coarser grids, cannot move it past its bound
    // (`grid_error`); coarser grids that do not see the part it does cannot
    // check it.
    if status == Status::Holds
        && let Some(why) = measured.unseen.clone()
    {
        status = Status::Unknown;
        let u = undecided.get_or_insert_with(Default::default);
        u.0.push(format!("it holds on a grid that cannot see the part: {why}"));
        u.1.push("refine (res r) to under half the gap's width, or close the gap in the model".into());
    }
    if status == Status::Holds
        && !measured.proven
        && let (Ok(x), Some(h)) = (measured.value.clone(), measured.res)
    {
        let mut coarse = |k: f64| -> Result<f64, String> {
            match v.measure(m, RunOpts { res: Some(k * h), tol: 1e-8 }) {
                Ok(c) => match c.unseen {
                    Some(why) => Err(format!("unseen: {why}")),
                    None => c.value,
                },
                Err(e) => Err(e.message),
            }
        };
        match coarse(2.0) {
            Err(why) if why.starts_with("unseen: ") => {
                status = Status::Unknown;
                let u = undecided.get_or_insert_with(Default::default);
                u.0.push(format!(
                    "it holds on cells of {h} mm, but cells of {} mm cannot see the part ({}), so cannot check it. The part may have a feature these cells see no better",
                    2.0 * h,
                    &why[8..]
                ));
                u.1.push(format!("refine it, (res r) at {} mm or finer, until two grids see the part and agree", h / 2.0));
            }
            Err(why) => measured.warnings.push(format!("not checked against coarser cells ({} mm): {why}", 2.0 * h)),
            Ok(x2) => {
                let x4 = coarse(4.0).ok();
                let (e, how) = grid_error(x, x2, x4);
                if point_status(c, x - e) != Status::Holds || point_status(c, x + e) != Status::Holds {
                    status = Status::Unknown;
                    let u = undecided.get_or_insert_with(Default::default);
                    u.0.push(format!(
                        "it holds at {} on cells of {h} mm, but its grid error is about {} ({how}: {} on cells of {} mm{}), more than its margin to the bound",
                        num(x),
                        num(e),
                        num(x2),
                        2.0 * h,
                        x4.map_or(String::new(), |x4| format!(", {} on cells of {} mm", num(x4), 4.0 * h)),
                    ));
                    u.1.push(format!("refine it, (res r) at {} mm or finer, until its value moves by less than its margin", h / 2.0));
                }
            }
        }
    }
    let undecided = undecided.map(|(why, how)| (why.join("; "), how.join("; ")));
    let quantity = v.doc.registry.metric_quantity(m.metric);
    let check = checked(status, Evidence::Metric { measured, head, form: m.form.clone(), quantity, undecided });
    Ok(if exact { check } else { check.estimate() })
}
