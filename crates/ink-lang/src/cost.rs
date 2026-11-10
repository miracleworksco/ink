//! What a claim or an interfere pair will cost, said before it is asked
//! (`ink check`; as `ink assert` and `ink interfere` start, for a stopped
//! call's answer): its work, counted from what is known cheaply, and where
//! what that work costs was measured, its class in time on this machine's
//! threads.
//!
//! Work is cells times the instructions an evaluation runs there. Each
//! shape is sketched once (feature::Sketch: its area, and the length of
//! the tape the octree simplified where its surface is): a body's by
//! `ink check`'s survey, another's as the survey would, cached.
//!
//! - `min-thickness`: a surface sample per cell of its res (else of its
//!   box's longest side / 128) on three grids (h, 2h, 4h), each tracing a
//!   ray. Classed.
//! - A value or a point: one evaluation. Classed.
//! - An integral (mass, volume, a centroid): the cells its surface crosses
//!   at its res; else at most its passes. Counted, not classed: its rate
//!   is not measured.
//! - A proof (subset, within, disjoint, touches, overlaps, clearance,
//!   empty, range, an interfere pair): its cells are where interval
//!   bounds cannot yet tell, which no cheap look predicts (faces seated
//!   on each other settle in a few cells; a curve on a curve, at cells of
//!   about the accuracy), so counted only by its budget and its field's
//!   size. Not classed.
//!
//! Each counted once per combination of the cases it reads (a text case
//! always: text shapes the model at load), and at each of the moments or
//! values a claim over a motion or ranged params is measured at.

use crate::Claim;
use crate::claims;
use crate::view::ModelView;
use ink_core::tape::Tape;
use ink_kernel::feature::Sketch;
use std::collections::HashMap;
use std::sync::Arc;

/// What a thickness costs, as measured on real projects' claims at 8
/// threads (benchmarks/cost.md): instructions a thread runs a second over
/// its samples and rays, and its start (bounds, the first grid).
const THICKNESS_RATE: f64 = 3.0e6;
const THICKNESS_START: f64 = 0.3;

/// How an estimate bounds the cost.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Bound {
    /// Good to about its class.
    About,
    /// Its work is counted, but what that work costs is not measured, or
    /// varies too much with the model to say: no class.
    Counted,
}

/// A claim's or a pair's estimated cost.
#[derive(Clone, Debug)]
pub struct Cost {
    /// Seconds on this machine's threads, where `About`.
    pub secs: f64,
    pub bound: Bound,
    /// What it was estimated from, in words.
    pub basis: String,
}

/// The classes an estimate is said in: each about ten times the last.
pub const CLASSES: [(f64, &str); 5] = [(1.0, "under a second"), (10.0, "seconds"), (60.0, "tens of seconds"), (600.0, "minutes"), (f64::INFINITY, "over ten minutes")];

impl Cost {
    /// None: not said (`Bound::Counted`).
    pub fn class(&self) -> Option<&'static str> {
        (self.bound == Bound::About).then(|| CLASSES.iter().find(|c| self.secs < c.0).map_or("over ten minutes", |c| c.1))
    }

    /// As an answer gives it: `{"class", "basis"}`, class null where it is
    /// not said.
    pub fn json(&self) -> serde_json::Value {
        serde_json::json!({"class": self.class(), "basis": self.basis})
    }
}

/// The bodies' sketches (their surveys', feature::Survey::sketch) by
/// field: all an estimate reads of a shape's surface. A shape no body is
/// has none, and its thickness or integral is not classed; nor is a
/// body's not surveyed yet (None).
#[derive(Default)]
pub struct Sketches(HashMap<u32, Option<Arc<Sketch>>>);

impl Sketches {
    /// Body field `f`'s sketch, from its survey.
    pub fn seed(&mut self, f: u32, s: Sketch) {
        self.0.insert(f, Some(Arc::new(s)));
    }

    /// Body field `f`'s sketch; Err where it is a body not surveyed yet.
    fn of(&self, f: u32) -> Result<Option<Arc<Sketch>>, ()> {
        match self.0.get(&f) {
            Some(None) => Err(()),
            Some(Some(s)) => Ok(Some(s.clone())),
            None => Ok(None),
        }
    }

    /// The sketches of the bodies whose fields are among `fields`, where
    /// `ink check`'s survey of them is cached already: for a caller that
    /// does not survey (ink assert), which must not be slowed by an
    /// estimate. The others are marked not surveyed yet.
    pub fn surveyed(v: &mut ModelView, fields: &[u32]) -> Sketches {
        let mut out = Sketches::default();
        let bodies: Vec<(String, u32)> = v.doc.bodies.iter().filter(|b| fields.contains(&b.shape)).map(|b| (b.name.clone(), b.shape)).collect();
        for (name, shape) in bodies {
            // Its bounds: what the claim reads first, in the view after.
            // Its survey as check asks it: check's view ranges nothing, so
            // a body reading a ranged param or t is keyed there at one
            // value of it, as here without the scope's spans.
            let check = ink_core::eval::Ctx { p_spans: ink_core::eval::ParamSpans::NONE, t_span: None, ..v.ctx.clone() };
            let found = v.bounds(shape).filter(|_| v.whole(shape)).and_then(|(lo, hi, _)| {
                let named = !v.doc.named_parts(&name, shape).is_empty();
                ink_kernel::feature::survey_cached(&Arc::new(Tape::new(&v.doc.g, shape)), &check, lo, hi, 128.0, named)
            });
            out.0.insert(shape, found.map(|s| Arc::new(s.sketch)));
        }
        out
    }
}

/// What a claim on a body not surveyed yet says.
const UNSURVEYED: &str = "no survey yet; ink check computes it";

/// A count in words: 1,234 or 1.2 million.
fn count(x: f64) -> String {
    match x {
        x if x >= 1e9 => format!("{:.1} billion", x / 1e9),
        x if x >= 1e6 => format!("{:.1} million", x / 1e6),
        x if x >= 1e4 => format!("{:.0} thousand", x / 1e3),
        x => format!("{:.0}", x),
    }
}

fn mm(x: f64) -> String {
    format!("{}", (x * 1e4).round() / 1e4)
}

fn threads() -> f64 {
    rayon::current_num_threads() as f64
}

/// A proof over `fields` (a claim's, an interfere pair's): its cells are where interval bounds cannot
/// yet tell, which no cheap look predicts (faces seated on each other
/// settle in a few; a curve on a curve, at cells of about the accuracy):
/// counted only by its budget, query::BOXES and its hunt as many again,
/// of a field of the instructions they take (fewer near most boxes, which
/// the octree simplifies it for).
pub fn proof(v: &ModelView, fields: &[u32]) -> Cost {
    let ops: usize = fields.iter().map(|&f| Tape::new(&v.doc.g, f).len()).sum();
    Cost {
        secs: f64::NAN,
        bound: Bound::Counted,
        basis: format!(
            "a proof: as many boxes as its bounds leave undecided, which no cheap look predicts (at most {} with its hunt), of a field of {} instructions",
            count(2.0 * ink_kernel::query::BOXES as f64),
            count(ops as f64)
        ),
    }
}

/// The bodies assert k's estimate reads the sketch of (`Sketches`): a
/// thickness's shape, an integral's parts'.
pub fn sketched(v: &ModelView, k: usize) -> Vec<u32> {
    let x = &v.doc.asserts[k];
    match v.doc.claim_at(&x.claim, x.bound, &v.ctx) {
        Claim::Metric { metric, .. } if metric.metric == "min-thickness" => metric.shape.into_iter().collect(),
        Claim::Integral { parts, .. } => parts.iter().map(|p| p.0).collect(),
        Claim::Mean { parts, .. } => parts.iter().map(|p| p.0).collect(),
        _ => Vec::new(),
    }
}

/// How many times `ink assert` asks a claim reading `params` for the
/// declared cases: once at the params' values and once at each
/// combination of the cases of the params it reads that `--param` does not
/// set (Doc::overrides; the rest are answered from the cache), and with
/// the params whose cases it reads. At most claims::MOST_CASES
/// combinations are asked (a run refuses more), so the count stops there.
fn cased(v: &ModelView, params: &[usize]) -> (usize, Vec<String>) {
    let (mut combos, mut names) = (1usize, Vec::new());
    for (name, given) in v.doc.cases.iter().filter(|(n, _)| !v.doc.overrides.iter().any(|(o, _)| o == n)) {
        // Text shapes the model at load: every claim may read it.
        let text = v.doc.texts.iter().any(|p| &p.name == name);
        if text || v.doc.params.iter().position(|p| &p.name == name).is_some_and(|k| params.contains(&k)) {
            combos = combos.saturating_mul(given.len() + 1);
            names.push(name.clone());
        }
    }
    (combos.min(claims::MOST_CASES + 1), names)
}

/// The cost of the document's assert k as `ink assert` will ask it (at
/// `run_res`, the run's --res), read from the bodies' sketches `sk`, so
/// `ink check` and a stoppable assert say the same: None where it is not
/// estimated (a package's solve, a shape not proven bounded).
pub fn of_assert(v: &mut ModelView, sk: &Sketches, k: usize, run_res: Option<f64>) -> Option<Cost> {
    let x = &v.doc.asserts[k];
    let (claim, bound, file, own_res) = (x.claim.clone(), x.bound, x.file.clone(), x.res);
    let c = v.doc.claim_at(&claim, bound, &v.ctx);
    let roots: Vec<u32> = v.doc.claim_roots(&c).into_iter().chain(bound).collect();
    let res = claims::claim_res(own_res, run_res).map(|r| r.0);
    let by_res = |h: f64, own: Option<f64>| if own.is_some_and(|o| o == h) { "its (res r)" } else { "--res" };
    let params = v.doc.params_read(&roots);
    // What assert ranges, as it scopes the question (ModelView::
    // scope_question): the declared ranges and tolerances it reads, and t
    // over the span of the file it is written in.
    let ranged: Vec<String> = params.iter().filter(|&&k| v.doc.has_range(k)).map(|&k| v.doc.params[k].name.clone()).collect();
    let moving = v.doc.span_of(&file).is_some() && claims::reads_t(v, &roots);
    let (runs, cased) = cased(v, &params);
    let over = {
        let mut o = Vec::new();
        if moving {
            o.push("t over the motion's span".to_string());
        }
        if !ranged.is_empty() {
            o.push(format!("the ranged params {}", ranged.join(", ")));
        }
        o.join(" and ")
    };
    let mut cost = match &c {
        Claim::Subset(a, b) | Claim::Disjoint(a, b) | Claim::Touches(a, b) | Claim::Overlaps(a, b) | Claim::Clearance(a, b, _) => proof(v, &[*a, *b]),
        Claim::Empty(f) | Claim::Range(f) => proof(v, &[*f]),
        // A value over ranges is a search over them (ranged::value_span;
        // over a motion, t ranged): not classed.
        Claim::Value(..) | Claim::Point { .. } if !over.is_empty() => Cost {
            secs: f64::NAN,
            bound: Bound::Counted,
            basis: format!("a search over {over}: interval bounds over pieces of their spans, halving undecided ones, as many as it leaves undecided"),
        },
        Claim::Value(..) | Claim::Point { .. } => Cost { secs: 0.0, basis: "one evaluation".into(), bound: Bound::About },
        // A package's metric reading a ranged param is refused (claims::
        // decide).
        Claim::Metric { .. } if !ranged.is_empty() => {
            Cost { secs: 0.0, bound: Bound::About, basis: format!("refused: a package's metric is solved at one value of each param, and it reads {over}") }
        }
        Claim::Metric { metric, .. } if metric.metric == "min-thickness" => {
            let f = metric.shape?;
            let (lo, hi, _) = v.bounds(f)?;
            let h = res.unwrap_or(ink_kernel::query::max_extent(lo, hi) / 128.0);
            match sk.of(f) {
                Err(()) => Cost { secs: f64::NAN, bound: Bound::Counted, basis: format!("rays from surface samples at {} mm; {UNSURVEYED}", mm(h)) },
                Ok(None) => Cost {
                    secs: f64::NAN,
                    bound: Bound::Counted,
                    basis: format!("rays from surface samples at {} mm; its shape is no body, whose surface ink check samples, so not counted", mm(h)),
                },
                Ok(Some(s)) => {
                    // h, 2h and 4h: 1 + 1/4 + 1/16.
                    let samples = 1.3125 * s.area / (h * h);
                    Cost {
                        secs: THICKNESS_START + samples * s.ops / (THICKNESS_RATE * threads()),
                        basis: format!(
                            "rays from about {} surface samples ({} mm² at {} mm, {}, and twice and four times that to check its grid error), of {} instructions an evaluation",
                            count(samples),
                            count(s.area),
                            mm(h),
                            if res.is_some() { by_res(h, own_res) } else { "a 128th of its box" },
                            count(s.ops)
                        ),
                        bound: Bound::About,
                    }
                }
            }
        }
        Claim::Integral { parts, .. } => integral(sk, parts.iter().map(|p| (p.0, p.1.len())), res),
        Claim::Mean { parts, .. } => integral(sk, parts.iter().map(|p| (p.0, p.2.len() + 1)), res),
        _ => return None,
    };
    if runs > 1 {
        cost.secs *= runs as f64;
        let most = if runs > claims::MOST_CASES { " at most" } else { "" };
        cost.basis += &format!(", asked{most} {runs} times for the cases of {}", cased.join(", "));
    }
    // Over a motion a thickness is measured at SAMPLES moments
    // (claims::samples); a proof ranges t in its one search.
    if moving && cost.bound == Bound::About && matches!(c, Claim::Metric { .. }) {
        cost.secs *= claims::SAMPLES as f64;
        cost.basis += &format!(", at {} moments of the motion", claims::SAMPLES);
    } else if !over.is_empty() && !cost.basis.contains(&over) {
        cost.basis += &format!(", over {over}");
    }
    if cost.bound == Bound::About {
        cost.basis += &format!("; on {} threads", threads());
    }
    Some(cost)
}

/// An integral's cost over parts (shape, how many integrands): the cells
/// its surface crosses at `res`; else at most its passes, each within
/// integral::PASS_CELLS, the next four times the last (estimate_parts). Not
/// classed: what a cell costs is not measured.
fn integral(sk: &Sketches, parts: impl Iterator<Item = (u32, usize)>, res: Option<f64>) -> Cost {
    let (mut work, mut cells, mut seen, mut surveyed) = (0.0f64, 0.0f64, true, true);
    for (shape, n) in parts {
        let s = match sk.of(shape) {
            Ok(Some(s)) => s,
            Ok(None) => {
                seen = false;
                continue;
            }
            Err(()) => {
                surveyed = false;
                continue;
            }
        };
        let c = match res {
            Some(h) => s.area / (h * h),
            None => 4.0 / 3.0 * ink_kernel::integral::PASS_CELLS as f64,
        };
        cells += c;
        work += c * s.ops * n.max(1) as f64;
    }
    let basis = match (seen, res) {
        _ if !surveyed => format!("an integral; {UNSURVEYED}"),
        (false, _) => "an integral over shapes no body is, whose surface ink check samples: not counted".to_string(),
        (true, Some(h)) => format!("an integral: about {} cells on its surface at {} mm, at {} instructions an evaluation", count(cells), mm(h), count(work / cells.max(1.0))),
        (true, None) => format!(
            "an integral: passes halved until its error is within its margin, at most {} cells, at {} instructions an evaluation",
            count(cells),
            count(work / cells.max(1.0))
        ),
    };
    Cost { secs: f64::NAN, bound: Bound::Counted, basis }
}
