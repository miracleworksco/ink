use super::*;

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
    /// A material point and heights of a positive-thickness air layer spanning its entire support cone.
    OverhangGap { witness: V3, at: At, layer_z: [f64; 2] },
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
    /// out would answer otherwise), and how to decide it. `stopped`: the
    /// call's time limit passed before the cells asked finished.
    Metric {
        measured: Measured,
        head: &'static str,
        form: String,
        quantity: Option<crate::Quantity>,
        undecided: Option<(String, String)>,
        stopped: Option<Short>,
        /// An estimate's basis (its cells, its grid check): none where exact.
        estimate: Option<String>,
        /// Its grid error (`grid_error`) and how it was found, where its
        /// grid check ran.
        grid_error: Option<(f64, &'static str)>,
    },
    /// A value at the given params and t.
    Value(f64),
    /// A value at every value of the ranged params it reads.
    ValueSpan(ValueSpan),
    /// A method that answers at one value of each param (rays, a package's
    /// solve), for a claim that reads a ranged one: not run, as its answer
    /// would say nothing about the rest of the range. Outside what the
    /// method covers, as is a rate with no bound (a drive that may jump).
    NotCovered(String),
    /// A search that did not prove it, saying where it stopped: something
    /// could still decide it.
    Unproven(String),
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
    /// (ModelView::refusal), whose error says why, and of which kind where
    /// the solve knows (`InkError::open`: the design lacks something, or it
    /// is outside what the solve covers; else it could not decide).
    Unsolved(InkError),
    /// A declared package limit lacks an input needed to evaluate it.
    NotChecked(&'static str),
    /// Not checked: the call's time limit had passed (ink_kernel::stop).
    Stopped,
    /// Not checked: outside ink's range (crate::validity), saying why.
    Outside(String),
}

/// A grid's estimate the call's time limit stopped before the cells asked
/// (`res`) finished (metric_claim): its value is the finest coarser cells'
/// that finished, with their grid error from the next coarser (None: none
/// finished, or they were the coarsest).
#[derive(Clone, Copy, Debug)]
pub struct Short {
    pub res: f64,
    pub error: Option<f64>,
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

/// How an emptiness claim (empty, disjoint, subset), or a void count, held.
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
    /// A void count to the accuracy: `cells` of air no deeper than it
    /// (each point within it of a point shown not air) set aside, as an
    /// emptiness claim's solid is, one at `at`.
    SetAside { cells: usize, at: Option<V3> },
}

/// Points measured along each side of a loose cell (Held::scoped).
const NEAR_SAMPLES: usize = 8;

/// A proof's deepest points below 0 tried for confirmation (Held::scoped).
const CONFIRM_TRIES: usize = 16;

impl Held {
    /// How the claim held, asked of shapes up to `charge` off the world's
    /// (a pivot's drift, a drive's slack, charged to its accuracy): nothing
    /// violating those is not nothing violating the world's, which may
    /// overlap them by up to the charge, so outright is to the accuracy
    /// (its near point not measured).
    pub fn charged(self, charge: f64) -> Held {
        match self {
            Held::Outright if charge > 0.0 => Held::ToAccuracy { near: None },
            h => h,
        }
    }

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
                below.extend(values.iter().zip(&points).filter(|x| ink_core::sign::inside(*x.0)).map(|(f, c)| (*f, *c, At::of(&at))));
            }
        }
        below.sort_by(|a, b| a.0.total_cmp(&b.0));
        below.into_iter().take(CONFIRM_TRIES).map(|x| (x.1, x.2)).collect()
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
        let deepest = points.into_iter().map(|(c, at)| (bound(c, &at), (c, at))).filter(|x| ink_core::sign::inside(x.0)).min_by(|a, b| a.0.total_cmp(&b.0));
        let depth = deepest.as_ref().map_or(0.0, |d| -d.0);
        Held::ToAccuracy { near: Some(Held::near(depth, deepest.map(|d| d.1))) }
    }

    /// A near point at least `depth` deep, rounded to a µm: none at 0.
    pub fn near(depth: f64, at: Option<(V3, At)>) -> (f64, Option<(V3, At)>) {
        let depth = (depth * 1e3).round().max(0.0) / 1e3;
        (depth, at.filter(|_| depth > 0.0))
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
    /// Over a motion that repeats, the one period is asked in pieces and not whole
    /// (motion::scoped): what a whole period's search left at its limit of boxes.
    pub in_pieces: bool,
}

/// One assert's answer.
#[derive(Clone, Debug)]
pub struct AssertResult {
    pub name: String,
    /// The claim as written.
    pub claim: String,
    /// Where it is written: `file:line`.
    pub source: String,
    /// Source of the limit checked by a generated placement claim.
    pub declared_at: Option<String>,
    /// The cell size it was checked at, if one was given, and whether its
    /// own `(res r)` gave it (the finer of that and the run's).
    pub res: Option<(f64, ResBy)>,
    pub check: Check,
    /// Wall-clock milliseconds it took to check.
    pub ms: f64,
    /// Answered from the cache: what it rests on was read from a cache and
    /// nothing was computed, or the identical proof was already made in this run.
    pub cached: bool,
    /// AssertDef::terms where the value is (its witness, or the given
    /// params and t); empty where there is no one point.
    pub terms: Vec<(String, f64)>,
    /// What gave it its accuracy, where not its default.
    pub accuracy_by: Option<AccuracyBy>,
    /// Where it fails at a point, the bodies there, for each shape it compares (`bodies_at`).
    pub at_witness: Vec<Vec<String>>,
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

/// Every claim's answer, of a run.
#[derive(Clone, Debug)]
pub struct AssertReport {
    pub asserts: Vec<AssertResult>,
    /// With `own`: how many asserts written in imported files were left
    /// out, unchecked. None: none were left out for being imported.
    pub imported_not_run: Option<usize>,
}

impl AssertReport {
    /// The answer's note on its claims asked at one moment
    /// (Check::one_moment), from what each recorded.
    pub fn note(&self) -> Option<String> {
        assert_note(&self.asserts)
    }
}

pub(super) fn assert_note(asserts: &[AssertResult]) -> Option<String> {
    one_moment_note(&asserts.iter().filter_map(|a| Some((format!("claim {}", a.name), a.check.one_moment?))).collect::<Vec<_>>())
}
