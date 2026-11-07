//! The design loop over a model view: an optimization's objective and
//! requirements measured at one design (`evaluate_design`), their exact
//! derivatives where the kernel has them (`design_gradient`), and the
//! optimizer driven over reloaded designs (`run_optimize`): the method of
//! moving asymptotes where every derivative is exact, Nelder-Mead
//! otherwise. A design's measures steer the search; the design returned
//! is decided by its requirements' claims, asked as `ink assert` asks them
//! (`verify`). Package metrics and their adjoints are
//! reached through the registry (`Outcome::metric`, `Outcome::adjoint`), so
//! no package is named here. Clients print the results.

use crate::InkError;
use crate::eval::{eval_dual, eval_point};
use crate::forms::{Objective, Optimize};
use crate::graph::Var;
use crate::lang::claims::{self, Check, ClaimOptions, Status};
use crate::lang::command::{Args, ParamArgs, claim_scope, param_args};
use crate::lang::view::ModelView;
use crate::lang::{Claim, Cmp, FormMetric, Measured, RunOpts};
use crate::optimize;
use crate::tape::Tape;
use ink_core::V3;
use ink_kernel::cache;
use ink_kernel::query::{self, max_extent};

/// A requirement at one design: its text; the value the search measured
/// and its violation, signed and scaled to be comparable (<= 0 within its
/// bound), which steer the search and decide nothing; and at the design
/// returned, its claim's answer asked as `ink assert` asks it (`answer`,
/// claims::asserted), which alone decides whether it is satisfied.
#[derive(Clone, Debug)]
pub struct Requirement {
    pub claim: String,
    pub measured: f64,
    pub violation: f64,
    pub answer: Option<Check>,
}

impl Requirement {
    /// Its claim holds on the design returned, as `ink assert` answers it.
    pub fn satisfied(&self) -> bool {
        self.answer.as_ref().is_some_and(|c| c.status == Status::Holds)
    }
}

/// An optimization's objective and requirements measured at one design,
/// with what the objective's measurement says of itself.
#[derive(Clone, Debug)]
pub struct Measures {
    pub objective: f64,
    pub requires: Vec<Requirement>,
    pub warnings: Vec<String>,
}

impl Measures {
    /// Every requirement within its bound by the search's measure: how the
    /// search ranks designs, not whether one is feasible (`verify`).
    pub fn within(&self) -> bool {
        self.requires.iter().all(|r| r.violation <= 0.0)
    }
}

/// Exact derivatives with respect to the params varied, per param; None
/// where there is none (finite differences are needed).
pub type Grad = Vec<Option<f64>>;

/// Requirement `(require ...)` values of package metrics are taken at this
/// solver tolerance.
const TOL: f64 = 1e-7;

/// An optimize whose objective or a requirement reads a `(ranged lo hi)`
/// (std's hidden, a view each), as a claim reads one (its scope, through
/// solves; one used at load is the scope's error): what reads one holds
/// only over every value of it, and a design is measured at one value of
/// each, which says nothing of the rest, so the run is refused, saying which.
fn quantified(v: &ModelView, o: &Optimize) -> Result<(), InkError> {
    if v.doc.ranged.is_empty() {
        return Ok(());
    }
    let scope = v.scope_spans(&[])?;
    let names = |roots: &[u32]| -> Vec<String> { v.ranged_reads_in(&scope, roots).into_iter().map(|(k, _)| v.doc.params[k as usize].name.clone()).collect() };
    let objective: Vec<u32> = match &o.objective {
        Objective::Volume(f) | Objective::Value(f) => vec![*f],
        Objective::Mass(f, rho) => vec![*f, *rho],
        Objective::Metric(m) => v.doc.metric_roots(m),
    };
    let mut reads: Vec<String> = Vec::new();
    let q = names(&objective);
    if !q.is_empty() {
        reads.push(format!("its objective reads {}", q.join(", ")));
    }
    for (text, claim, bound) in &o.requires {
        let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(*bound).collect();
        let q = names(&roots);
        if !q.is_empty() {
            reads.push(format!("(require {text}) reads {}", q.join(", ")));
        }
    }
    if reads.is_empty() {
        return Ok(());
    }
    Err(InkError {
        hint: Some("check what reads them with ink assert, which proves it over every value; optimize over what reads none".into()),
        ..InkError::msg(format!("optimize {}: {}, made by (ranged lo hi), and a design is measured at one value of each, which says nothing of the rest", o.name, reads.join("; ")))
    })
}

/// Objective value and requirements for one design (one view).
pub fn evaluate_design(v: &mut ModelView, name: &str) -> Result<Measures, InkError> {
    let (objective, requires) = requirements(v, name)?;
    measure_design(v, objective, requires)
}

/// The objective and requirements of optimization `name` at a design, with
/// what they read provided.
fn requirements(v: &mut ModelView, name: &str) -> Result<(Objective, Vec<(String, Claim)>), InkError> {
    let (objective, requires): (Objective, Vec<(String, Claim)>) = {
        let o = v.doc.pick::<Optimize>("optimize", Some(name))?;
        quantified(v, o)?;
        (o.objective.clone(), o.requires.iter().map(|(t, c, b)| (t.clone(), v.doc.claim_at(c, *b, &v.ctx))).collect())
    };
    // What is provided only where read (a reaction) is provided for what
    // the objective and requirements read, as for a claim: else it reads
    // its unset value.
    let mut roots: Vec<u32> = requires.iter().flat_map(|(_, c)| v.doc.claim_roots(c)).collect();
    if let Objective::Value(f) = objective {
        roots.push(f);
    }
    if let Objective::Metric(m) = &objective
        && m.shape.is_some()
    {
        roots.extend(v.doc.metric_roots(m));
    }
    v.want(&roots)?;
    Ok((objective, requires))
}

/// Design's metrics of a shape (Package::shape_metrics): `min-thickness`,
/// an estimate, the thinnest wall by rays from surface samples spaced
/// `res` (ModelView::thickness: 128 along its box), where it is (`at`);
/// `cavities`, the enclosed voids by proof (ModelView::voids): those it
/// proved, and, where its cells leave more possible, the most there may be
/// (`undecided-cells`, excluded). An empty shape has no wall and no void.
pub(crate) fn shape_metric(v: &mut ModelView, metric: &str, f: u32, res: Option<f64>) -> Result<Measured, InkError> {
    let measured = |value, at, res| Measured { value, at, res, warnings: Vec::new(), excluded: Vec::new(), unseen: None, proven: false };
    Ok(match metric {
        "min-thickness" => match v.thickness(f, res) {
            None => measured(Ok(f64::INFINITY), None, None),
            Some(Err(h)) => measured(Err(crate::lang::view::UNSAMPLED.into()), None, Some(h)),
            Some(Ok((samples, h))) => {
                let min = samples.iter().min_by(|a, b| a.3.total_cmp(&b.3)).expect("at least one sample");
                measured(Ok(min.3), Some(min.0), Some(h))
            }
        },
        _ => match v.voids(f, res)? {
            None => Measured { proven: true, ..measured(Ok(0.0), None, None) },
            Some(c) => {
                let mut m = Measured { proven: true, ..measured(Ok(c.at_least as f64), c.count().is_none().then_some(c.at).flatten(), None) };
                if let Some((what, decide)) = c.why_undecided() {
                    m.excluded.push(crate::lang::Excluded { key: "undecided-cells", value: c.at_most.map_or(f64::INFINITY, |n| n as f64), at: c.at, what, decide });
                }
                m
            }
        },
    })
}

/// Why the objective has no value of the design, where a claim on it would
/// be unknown for the same reason (claims::decide): an integral or a
/// shape's metric of a shape not proven bounded (the window may cut it, and
/// its faces would be the shape's walls).
fn unbounded(v: &mut ModelView, f: u32, what: &str) -> Result<(), InkError> {
    if v.whole(f) {
        return Ok(());
    }
    Err(InkError {
        hint: Some("bound the shape (intersect it with a box)".into()),
        ..InkError::msg(format!("the objective's {what} is of a shape not proven bounded: the window may cut it, so it has no value (a claim on it is unknown, unbounded)"))
    })
}

/// The objective's value at a design, and its requirements' values and
/// violations, which steer the search. The objective fails the design
/// where a claim on it would have no value: of a shape the window may cut
/// (`unbounded`), a metric with none (a failed solve), or one measured on a
/// grid that cannot see the part (`Measured::unseen`: another part's
/// value).
fn measure_design(v: &mut ModelView, objective: Objective, requires: Vec<(String, Claim)>) -> Result<Measures, InkError> {
    let opts = RunOpts { res: None, tol: TOL };
    // A proven metric that leaves values undecided (a void count's
    // bounds apart) is not known at this design: the search steers clear
    // of it, as of one with no value.
    let measured = |v: &mut ModelView, m: &FormMetric| -> Result<(f64, Measured), InkError> {
        let r = v.measure(m, opts)?;
        let no = |why: String| InkError::msg(format!("'{}' gave no '{}' at this design: {why}", m.form, m.metric));
        let x = r.value.clone().map_err(no)?;
        if let Some(ex) = r.excluded.first().filter(|_| r.proven) {
            let most = if ex.value.is_finite() { format!("{}", ex.value) } else { "unbounded".into() };
            return Err(no(format!("its value is not known, only that it is from {x} to {most}: it leaves out {}", ex.what)));
        }
        Ok((x, r))
    };
    // The integral of field g (1 where none: the volume) over shape f in
    // the window, estimated.
    let integral = |v: &mut ModelView, f: u32, g: Option<u32>| -> f64 {
        match v.bounds(f) {
            None => 0.0,
            Some((lo, hi, _)) => {
                let t = v.clipped(f);
                let g = g.map(|g| Tape::new(&v.doc.g, g));
                query::sums(&t, &[g.as_ref()], false, &v.ctx, lo, hi, max_extent(lo, hi) / 100.0).sums[0]
            }
        }
    };
    let mut warnings = Vec::new();
    let value = match objective {
        Objective::Volume(f) => {
            unbounded(v, f, "volume")?;
            integral(v, f, None)
        }
        Objective::Mass(f, rho) => {
            unbounded(v, f, "mass")?;
            integral(v, f, Some(rho))
        }
        Objective::Metric(m) => {
            if let Some(f) = m.shape {
                unbounded(v, f, m.metric)?;
            }
            let (x, r) = measured(v, &m)?;
            if let Some(why) = r.unseen {
                return Err(InkError::msg(format!(
                    "the objective's '{}' was measured on a grid that cannot see the part ({why}): refine (res r) to under half the gap's width, or close the gap",
                    m.metric
                )));
            }
            warnings.extend(r.warnings);
            x
        }
        Objective::Value(f) => v.doc.scalar(f, &v.ctx),
    };
    // The gap between two shapes: `ink distance`'s search (claims::distance),
    // its confirmed pair's distance (within the model's accuracy of the
    // minimum), else the lower bound it reached; none where one is empty.
    let gap = |v: &mut ModelView, a: u32, b: u32| -> Result<Option<f64>, InkError> {
        if v.bounds(a).is_none() || v.bounds(b).is_none() {
            return Ok(None);
        }
        let d = crate::lang::claims::distance(v, a, b, Default::default())?;
        Ok(Some(d.upper.map_or(d.lower, |u| u.0)))
    };
    let rel = |x: f64, bound: f64, cmp: Cmp| -> f64 {
        let scale = bound.abs().max(1e-9);
        match cmp {
            Cmp::Lt | Cmp::Le => (x - bound) / scale,
            Cmp::Gt | Cmp::Ge => (bound - x) / scale,
        }
    };
    let worst = |xs: &[f64], bound: f64, cmp: Cmp| -> (f64, f64) { xs.iter().map(|&x| (x, rel(x, bound, cmp))).max_by(|a, b| a.1.total_cmp(&b.1)).unwrap_or((0.0, 0.0)) };
    let mut out = Vec::new();
    for (text, claim) in requires {
        let (measured, violation) = match claim {
            Claim::Empty(f) | Claim::Range(f) => {
                let x = integral(v, f, None);
                (x, libm::cbrt(x))
            }
            Claim::Disjoint(a, b) => {
                let f = v.doc.g.max(a, b);
                let x = integral(v, f, None);
                (x, libm::cbrt(x))
            }
            Claim::Subset(a, b) => {
                let nb = v.doc.g.neg(b);
                let f = v.doc.g.max(a, nb);
                let x = integral(v, f, None);
                (x, libm::cbrt(x))
            }
            // The overlap's volume, else the gap: 0 where they touch.
            Claim::Touches(a, b) => {
                let both = v.doc.g.max(a, b);
                let overlap = integral(v, both, None);
                match gap(v, a, b)? {
                    _ if overlap > 0.0 => (overlap, libm::cbrt(overlap)),
                    Some(dist) => (dist, dist),
                    None => (0.0, 1.0),
                }
            }
            // The field's value at the point, signed so that below 0 is
            // within.
            Claim::Point { f, at, inside, .. } => {
                let p = at.map(|n| v.doc.scalar(n, &v.ctx));
                let x = ink_core::eval::eval_point(&ink_core::tape::Tape::new(&v.doc.g, f), &v.ctx, p);
                (x, if inside { x } else { -x })
            }
            Claim::Clearance(fa, fb, d) => {
                let both = v.doc.g.max(fa, fb);
                let overlap = integral(v, both, None);
                match gap(v, fa, fb)? {
                    _ if overlap > 0.0 => (0.0, 1.0 + libm::cbrt(overlap)),
                    Some(dist) => (dist, (d - dist) / d.max(1e-9)),
                    // One shape empty: nothing to come near.
                    None => (f64::INFINITY, -1.0),
                }
            }
            Claim::Metric { metric: m, cmp, value: bound, .. } => {
                let (x, _) = measured(v, &m).map_err(|e| InkError::msg(format!("{text}: {}", e.message)))?;
                // A count: < b is <= the largest whole number under b, so
                // one void more than allowed is a violation of 1, not 0 at b.
                let violation = match (m.metric, cmp) {
                    ("cavities", Cmp::Lt) => x - (bound.ceil() - 1.0),
                    ("cavities", Cmp::Le) => x - bound.floor(),
                    ("cavities", Cmp::Gt) => (bound.floor() + 1.0) - x,
                    ("cavities", Cmp::Ge) => bound.ceil() - x,
                    _ => rel(x, bound, cmp),
                };
                (x, violation)
            }
            Claim::Value(f, cmp, bound) => {
                let x = v.doc.scalar(f, &v.ctx);
                (x, rel(x, bound, cmp))
            }
            // Of several fields, the one violated most.
            Claim::Integral { parts, cmp, bound, .. } => {
                let k = parts.first().map_or(0, |p| p.1.len());
                let xs: Vec<f64> = (0..k).map(|i| parts.iter().map(|(s, gs)| integral(v, *s, Some(gs[i]))).sum()).collect();
                worst(&xs, bound, cmp)
            }
            Claim::Mean { parts, cmp, bound, .. } => {
                let k = parts.first().map_or(0, |p| p.2.len());
                let w: f64 = parts.iter().map(|&(s, wf, _)| integral(v, s, Some(wf))).sum();
                let xs: Vec<f64> = (0..k)
                    .map(|i| {
                        let a: f64 = parts
                            .iter()
                            .map(|(s, wf, hs)| {
                                let wh = v.doc.g.mul(*wf, hs[i]);
                                integral(v, *s, Some(wh))
                            })
                            .sum();
                        a / w
                    })
                    .collect();
                worst(&xs, bound, cmp)
            }
            Claim::Overlaps(a, b) => {
                // Negative overlap volume as the violation (0 once they meet).
                let f = v.doc.g.max(a, b);
                let x = integral(v, f, None);
                (x, if x > 0.0 { -libm::cbrt(x) } else { 1.0 })
            }
        };
        out.push(Requirement { claim: text, measured, violation, answer: None });
    }
    Ok(Measures { objective: value, requires: out, warnings })
}

/// Exact derivatives of an optimization's objective and signed constraints
/// with respect to params `vary` (indices), where they exist: volume and
/// mass (boundary integrals), package metrics with an adjoint (compliance:
/// Hadamard's formula, one solve), and expressions of params (forward-mode
/// duals). None marks what needs finite differences: params baked in at
/// load time, params an adjoint does not cover (that a study's material,
/// loads or supports read: any input but its shape), params a mass's
/// density reads, and non-smooth claims.
pub fn design_gradient(v: &mut ModelView, name: &str, vary: &[usize]) -> Result<(Grad, Vec<Grad>), InkError> {
    let o = v.doc.pick::<Optimize>("optimize", Some(name))?;
    let (objective, requires): (Objective, Vec<Claim>) = (o.objective.clone(), o.requires.iter().map(|(_, c, b)| v.doc.claim_at(c, *b, &v.ctx)).collect());
    // A bound written as an expression of params, (value (+ a b) <= (- 30
    // c)), moves with them: its derivative comes off the claim's.
    let bounds: Vec<_> = o.requires.iter().map(|r| r.2).collect();
    let baked = v.doc.baked_params();
    let none = || vec![None; vary.len()];
    // d/dp of a boundary-weighted integral over field f.
    let boundary = |v: &mut ModelView, f: u32, w: &(dyn Fn(V3) -> f64 + Sync)| -> Grad {
        let Some((lo, hi, _)) = v.bounds(f) else { return vec![Some(0.0); vary.len()] };
        let t = v.clipped(f);
        let (_, d) = crate::shape_derivatives(&t, &[w], &v.ctx, lo, hi, max_extent(lo, hi) / 100.0);
        vary.iter().map(|&k| (!baked.contains(&k)).then(|| d[0][k])).collect()
    };
    // A package metric's derivative through its adjoint, if it has one.
    let adjoint = |v: &mut ModelView, m: &FormMetric| -> Result<Grad, InkError> {
        // A shape's metric (a thickness, a count of voids) has none.
        if m.shape.is_some() {
            return Ok(none());
        }
        let out = v.outcome(Some(m.package), &m.form, RunOpts { res: None, tol: TOL })?;
        let Some(a) = out.adjoint(m.metric) else { return Ok(none()) };
        let g = boundary(v, a.shape, &*a.weight);
        let uncovered = |k: usize| a.uncovered.iter().any(|&f| v.doc.g.depends_on(f, |x| x == Var::Param(k as u32)));
        Ok(g.into_iter().zip(vary).map(|(g, &k)| if uncovered(k) { None } else { g }).collect())
    };
    let quantity = |v: &mut ModelView, q: &Objective| -> Result<Grad, InkError> {
        Ok(match q {
            Objective::Volume(f) => boundary(v, *f, &|_| 1.0),
            Objective::Mass(f, rho) => {
                let r = Tape::new(&v.doc.g, *rho);
                let ctx = v.ctx.clone();
                let g = boundary(v, *f, &move |p| eval_point(&r, &ctx, p));
                // The boundary term only: a density that reads the param
                // changes the mass inside too.
                g.into_iter().zip(vary).map(|(g, &k)| if v.doc.g.depends_on(*rho, |x| x == Var::Param(k as u32)) { None } else { g }).collect()
            }
            Objective::Metric(m) => adjoint(v, m)?,
            Objective::Value(f) => {
                let d = eval_dual(&Tape::new(&v.doc.g, *f), &v.ctx, [0.0; 3]);
                vary.iter().map(|&k| (!baked.contains(&k)).then(|| d[5 + k])).collect()
            }
        })
    };
    let g0 = quantity(v, &objective)?;
    let mut gs = Vec::new();
    for c in &requires {
        let rel = |bound: f64, cmp: Cmp, g: Grad| -> Grad {
            let k = match cmp {
                Cmp::Lt | Cmp::Le => 1.0,
                Cmp::Gt | Cmp::Ge => -1.0,
            } / bound.abs().max(1e-9);
            g.into_iter().map(|x| x.map(|x| k * x)).collect()
        };
        let less_bound = |g: Grad, v: &mut ModelView| -> Result<Grad, InkError> {
            let Some(b) = bounds[gs.len()] else { return Ok(g) };
            let db = quantity(v, &Objective::Value(b))?;
            Ok(g.into_iter().zip(db).map(|(g, db)| Some(g? - db?)).collect())
        };
        // A volume, `(volume f ...)`: the integral of 1 over one shape.
        let volume = match c {
            Claim::Integral { parts, frames, .. } if frames.iter().all(Option::is_none) => match &parts[..] {
                [(f, gs)] if gs.len() == 1 && v.doc.g.as_const(gs[0]) == Some(1.0) => Some(*f),
                _ => None,
            },
            _ => None,
        };
        gs.push(match c {
            Claim::Integral { cmp, bound, .. } if volume.is_some() => rel(*bound, *cmp, less_bound(quantity(v, &Objective::Volume(volume.expect("a volume")))?, v)?),
            Claim::Metric { metric, cmp, value, .. } => rel(*value, *cmp, less_bound(adjoint(v, metric)?, v)?),
            Claim::Value(f, cmp, b) => rel(*b, *cmp, less_bound(quantity(v, &Objective::Value(*f))?, v)?),
            _ => none(),
        });
    }
    Ok((g0, gs))
}

/// What an optimization found.
pub struct Optimized {
    pub name: String,
    /// The params varied, with their bounds.
    pub vary: Vec<(String, f64, f64)>,
    pub maximize: bool,
    /// The best design (values of `vary`), its objective, whether it is
    /// feasible, and its requirements: feasible where every requirement's
    /// claim holds there as `ink assert` answers it (`verify`).
    pub design: Vec<f64>,
    pub objective: f64,
    pub feasible: bool,
    pub requires: Vec<Requirement>,
    /// The first design evaluated: the model's params, clamped to the
    /// bounds; its objective, and whether its requirements were within
    /// their bounds by the search's measure (Measures::within).
    pub start: Vec<f64>,
    pub start_objective: f64,
    pub start_within: bool,
    pub trace: Trace,
    /// Why the design returned has no objective: its evaluation failed (a
    /// study that failed at it, a metric with no value there).
    pub failure: Option<String>,
    /// The names of the params the requirements' answers index (`verify`'s
    /// view's, which may hold t besides the document's).
    pub answer_params: Vec<String>,
    /// What the reader must know of the run (it ended where it started,
    /// designs it could not evaluate).
    pub warnings: Vec<String>,
    /// The designs the search tried that could not be evaluated (one the
    /// model refuses there, a study that failed), each with its error.
    pub failed: Vec<(Vec<f64>, InkError)>,
    /// The moment the search measured each design at, where its
    /// requirements were verified over a motion (`verify`): the search's
    /// measure is of that moment only.
    pub search_t: Option<f64>,
}

/// How each method got there.
pub enum Trace {
    /// Bounded Nelder-Mead with a constraint penalty: the number of designs
    /// evaluated, and each that improved on the best feasible design so far
    /// (its evaluation, counted from 1, design and objective).
    Simplex { evaluations: usize, improvements: Vec<(usize, Vec<f64>, f64)> },
    /// The method of moving asymptotes: its iterations, each design
    /// evaluated with how many evaluations its derivatives took, and how
    /// many param derivatives were exact or differenced.
    Gradient { iterations: usize, history: Vec<Iteration>, exact: usize, differenced: usize },
}

pub struct Iteration {
    pub design: Vec<f64>,
    pub objective: f64,
    pub max_violation: f64,
    /// Design evaluations: one, plus two per param differenced.
    pub evaluations: usize,
}

/// A run's designs, reported as each is evaluated (lang::claims::progress:
/// `ink mcp` says "12 of 40 designs done" of a long run, and a run stopped
/// at its time limit answers with the best design so far and every design
/// evaluated): its number, the design, whether it met the requires (or
/// failed), its objective, the requires measured, and its time. Each design
/// better than the best so far (a feasible one by its objective, else the
/// least violating) is also reported as the answer so far, in the fields a
/// finished run gives (`design`, `objective`, `feasible`, `requires`,
/// `apply`).
/// `apply`: the options that ask ink assert of a design what the run's
/// verification asks (command::claim_scope): each param's value as the run
/// gave it, in the unit it was given in (a pin, name=v:v, kept a pin: its
/// tolerance dropped), the varied ones' at the design, exactly; and the
/// run's moment, --t. One for every answer, in progress or finished.
pub struct Replay {
    given: ParamArgs,
    t: Option<String>,
}

impl Replay {
    pub fn of(a: &Args) -> Result<Replay, InkError> {
        Ok(Replay { given: param_args(a)?, t: a.get("t").map(String::from) })
    }

    /// The options for the design `x` of the params `names`.
    pub fn apply(&self, names: &[String], x: &[f64]) -> String {
        let given = &self.given;
        let pin = |k: &str| given.pinned.iter().any(|p| p == k);
        let unit = |k: &str| given.values.iter().find(|(n, _)| n == k).and_then(|(_, g)| g.unit);
        let value = |k: &str, x: f64| {
            let x = unit(k).map_or(format!("{x}"), |u| ink_lang::units::write(x, u));
            if pin(k) { format!("{k}={x}:{x}") } else { format!("{k}={x}") }
        };
        let fixed = given.values.iter().filter(|(k, _)| !names.contains(k)).map(|(k, g)| value(k, g.value));
        let varied = names.iter().zip(x).map(|(k, x)| value(k, *x));
        let flag: Vec<String> = fixed.chain(varied).collect();
        let mut apply = format!("--param {}", flag.join(","));
        if let Some(t) = &self.t {
            apply += &format!(" --t {t}");
        }
        apply
    }
}

struct Trials<'a> {
    names: &'a [String],
    maximize: bool,
    count: std::cell::Cell<usize>,
    /// The best so far: (feasible, objective or total violation, signed so
    /// that smaller is better).
    best: std::cell::Cell<Option<(bool, f64)>>,
    /// Each design's view (ModelView::reload).
    design: &'a dyn Fn(&[f64]) -> Result<ModelView, InkError>,
    replay: &'a Replay,
    /// The designs that could not be evaluated, each once, with why.
    failed: std::cell::RefCell<Vec<(Vec<f64>, InkError)>>,
}

impl<'a> Trials<'a> {
    fn new(run: &Optimized, names: &'a [String], method: &str, design: &'a dyn Fn(&[f64]) -> Result<ModelView, InkError>, replay: &'a Replay) -> Trials<'a> {
        use serde_json::json;
        crate::lang::claims::progress(json!({"answer": {"optimize": run.name, "sense": if run.maximize { "maximize" } else { "minimize" }, "method": method}}));
        Trials { names, maximize: run.maximize, count: Default::default(), best: Default::default(), design, replay, failed: Default::default() }
    }

    /// The design `x` as a view, its measures (`evaluate_design`), and the
    /// view. A design that cannot be loaded or evaluated is reported, in
    /// progress and kept for the answer (`failed`), with its error.
    fn evaluate(&self, x: &[f64], name: &str) -> Result<(ModelView, Measures), InkError> {
        use serde_json::json;
        let t0 = std::time::Instant::now();
        let r = (self.design)(x).and_then(|mut d| cache::trial(|| evaluate_design(&mut d, name)).map(|m| (d, m)));
        self.count.set(self.count.get() + 1);
        let design: serde_json::Map<String, serde_json::Value> = self.names.iter().zip(x).map(|(p, x)| (p.clone(), json!(x))).collect();
        let mut row = json!({"evaluation": self.count.get(), "design": design});
        let status = match r.as_ref().map(|(_, m)| m) {
            Ok(m) => {
                let (v, c, feasible) = (&m.objective, &m.requires, m.within());
                row["objective"] = json!(v);
                row["search_feasible"] = json!(feasible);
                row["requires"] = json!(c.iter().map(|c| json!({"claim": c.claim, "search_measured": c.measured, "search_violation": c.violation.max(0.0)})).collect::<Vec<_>>());
                let score = if feasible { if self.maximize { -v } else { *v } } else { c.iter().map(|c| c.violation.max(0.0)).sum() };
                let better = match self.best.get() {
                    None => !score.is_nan(),
                    Some((f, s)) => (feasible && !f) || (feasible == f && score < s),
                };
                if better {
                    self.best.set(Some((feasible, score)));
                    let mut answer = row.clone();
                    answer["apply"] = json!(self.replay.apply(self.names, x));
                    crate::lang::claims::progress(json!({"answer": answer}));
                }
                // By the search's measure: a design's claims are asked only
                // at the one returned (verify).
                if feasible { "within" } else { "violating" }
            }
            Err(e) => {
                row["error"] = e.to_json();
                // Each design once (a simplex may try one again).
                let mut failed = self.failed.borrow_mut();
                if !failed.iter().any(|(y, _)| y == x) {
                    failed.push((x.to_vec(), e.clone()));
                }
                // In the answer so far too (a run stopped early answers it).
                let warned: Vec<String> = failed_warning(self.names, &failed).into_iter().collect();
                crate::lang::claims::progress(json!({"answer": {"failed_designs": failed_json(self.names, &failed), "warnings": warned}}));
                "failed"
            }
        };
        crate::lang::claims::progress(json!({"done": format!("design {}", self.count.get()), "status": status, "ms": t0.elapsed().as_secs_f64() * 1e3, "row": row}));
        r
    }
}

/// The designs a run could not evaluate, as its answer gives them: each
/// design and its error (message, where, hint).
pub fn failed_json(names: &[String], failed: &[(Vec<f64>, InkError)]) -> serde_json::Value {
    use serde_json::json;
    failed.iter().map(|(x, e)| json!({"design": names.iter().zip(x).map(|(p, x)| (p.clone(), json!(x))).collect::<serde_json::Map<_, _>>(), "error": e.to_json()})).collect()
}

/// The warning for designs a run could not evaluate: the search kept clear
/// of them, so the design it returns may be limited by the model's error,
/// not by a requirement.
fn failed_warning(names: &[String], failed: &[(Vec<f64>, InkError)]) -> Option<String> {
    let (x, e) = failed.first()?;
    let at: Vec<String> = names.iter().zip(x).map(|(p, x)| format!("{p}={x}")).collect();
    let n = failed.len();
    Some(format!(
        "{n} design{} could not be evaluated (failed_designs), and the search kept clear of them: the first, at {}: {}",
        if n == 1 { "" } else { "s" },
        at.join(", "),
        e.message
    ))
}

/// Run the optimization `name` of a view loaded from a file: each design
/// is the view reloaded with its params (ModelView::reload), from the
/// view's own values, and a trial (cache::trial: its solves are not
/// stored). `budget` overrides the form's.
///
/// The method is the method of moving asymptotes where every derivative
/// the run needs (the objective's and each requirement's, by every param
/// varied) is exact at the start design (`design_gradient`), else
/// Nelder-Mead: differences of a metric a grid measures see its steps,
/// not its trend, and MMA on them stalled where it started (examples/
/// shelf.ink: 2 iterations, its start design returned as the optimum;
/// Nelder-Mead finds a design 40% lighter). Either way a run whose answer
/// is its start design says so. `options` are the run's: the design
/// returned is verified under the claim scope they give, as `ink assert`
/// with them would ask its claims (command::claim_scope).
pub fn run_optimize(v: &ModelView, name: &str, budget: Option<usize>, options: &Args) -> Result<Optimized, InkError> {
    let o = v.doc.pick::<Optimize>("optimize", Some(name))?;
    // Refused before the first design, the answer saying why.
    quantified(v, o)?;
    let replay = Replay::of(options)?;
    let name = name.to_string();
    let (vary, maximize) = (o.vary.clone(), o.maximize);
    let budget = budget.unwrap_or(o.budget);
    let lo: Vec<f64> = vary.iter().map(|v| v.1).collect();
    let hi: Vec<f64> = vary.iter().map(|v| v.2).collect();
    let x0: Vec<f64> = vary.iter().map(|(p, l, h)| v.ctx.params[v.doc.param_index(p).unwrap()].clamp(*l, *h)).collect();
    let params: Vec<String> = vary.iter().map(|p| p.0.clone()).collect();
    let design = |x: &[f64]| -> Result<ModelView, InkError> {
        let changes: Vec<(String, ink_lang::Given)> = params.iter().zip(x).map(|(p, x)| (p.clone(), (*x).into())).collect();
        cache::trial(|| v.reload(&changes))
    };
    let run = Optimized {
        name,
        vary,
        maximize,
        design: Vec::new(),
        objective: f64::NAN,
        feasible: false,
        requires: Vec::new(),
        start: Vec::new(),
        start_objective: f64::NAN,
        start_within: false,
        trace: Trace::Simplex { evaluations: 0, improvements: Vec::new() },
        failure: None,
        warnings: Vec::new(),
        failed: Vec::new(),
        answer_params: Vec::new(),
        search_t: None,
    };
    let idx: Vec<usize> = run.vary.iter().map(|p| v.doc.param_index(&p.0).unwrap()).collect();
    let exact =
        design(&x0).and_then(|mut d| cache::trial(|| design_gradient(&mut d, &run.name, &idx))).is_ok_and(|(g0, gs)| g0.iter().chain(gs.iter().flatten()).all(Option::is_some));
    if exact {
        crate::lang::claims::progress(serde_json::json!({"of": "designs", "kept": false}));
        return optimize_gradient(run, &idx, &params, budget, &lo, &hi, &x0, &design, &replay).and_then(|o| verify(o, &design, options)).map(|o| ended(o, &lo, &hi));
    }
    let name = &run.name;
    crate::lang::claims::progress(serde_json::json!({"total": budget, "of": "designs", "kept": false}));
    let trials = Trials::new(&run, &params, "derivative-free", &design, &replay);
    let mut details: Vec<Option<Measures>> = Vec::new();
    let mut failures: Vec<Option<String>> = Vec::new();
    let outcome = optimize::minimize(&lo, &hi, &x0, budget, |x| {
        let r = trials.evaluate(x, name).map(|(_, m)| m);
        // One record per evaluation, failed or not: they are looked up by
        // their position in the optimizer's history (so a failed evaluation
        // cannot shift later designs' requirements).
        let m = match r {
            Ok(m) => m,
            Err(e) => {
                details.push(None);
                failures.push(Some(e.message.clone()));
                return Err(e.message);
            }
        };
        let eval = optimize::Eval { objective: if maximize { -m.objective } else { m.objective }, violations: m.requires.iter().map(|c| c.violation.max(0.0)).collect() };
        details.push(Some(m));
        failures.push(None);
        Ok(eval)
    });
    let best = &outcome.best;
    let idx = outcome.history.iter().position(|h| h.x == best.x).unwrap_or(0);
    let start = &outcome.history[0];
    let mut best_so_far = f64::INFINITY;
    let mut improvements = Vec::new();
    for (k, h) in outcome.history.iter().enumerate() {
        if h.eval.feasible() && h.eval.objective < best_so_far {
            best_so_far = h.eval.objective;
            improvements.push((k + 1, h.x.clone(), details[k].as_ref().map_or(f64::NAN, |m| m.objective)));
        }
    }
    let found = details.get(idx).cloned().flatten();
    let failed = trials.failed.take();
    let mut warnings = found.as_ref().map(|m| m.warnings.clone()).unwrap_or_default();
    warnings.extend(failed_warning(&params, &failed));
    Ok(Optimized {
        design: best.x.clone(),
        objective: found.as_ref().map_or(f64::NAN, |m| m.objective),
        feasible: best.eval.feasible(),
        start: start.x.clone(),
        start_objective: details[0].as_ref().map_or(f64::NAN, |m| m.objective),
        start_within: start.eval.feasible(),
        warnings,
        failed,
        requires: found.map(|m| m.requires).unwrap_or_default(),
        trace: Trace::Simplex { evaluations: outcome.history.len(), improvements },
        failure: failures.get(idx).cloned().flatten(),
        ..run
    })
    .and_then(|o| verify(o, &design, options))
    .map(|o| ended(o, &lo, &hi))
}

/// The design returned, each requirement asked as `ink assert` asks its
/// claim there (claims::asserted: at its own resolution, with its own
/// checks of a grid's error, of what a metric leaves out, of a grid that
/// cannot see the part, of a shape the window may cut): it is satisfied
/// exactly where the claim holds, and the design feasible where each is.
/// An unknown is not satisfied, and its answer says why. The design's view
/// is scoped as `ink assert` scopes its session with the run's options
/// (command::claim_scope: the motion, or each file's own span; declared
/// tolerances around the design's values), and each requirement over the
/// span of the file it is written in, as each assert is.
fn verify(mut o: Optimized, design: &dyn Fn(&[f64]) -> Result<ModelView, InkError>, options: &Args) -> Result<Optimized, InkError> {
    o.feasible = false;
    if o.failure.is_some() || o.design.is_empty() {
        return Ok(o);
    }
    let mut v = design(&o.design)?;
    let searched_at = v.ctx.t;
    let scope = claim_scope(&v, options)?;
    v.set_claim_scope(&scope)?;
    // This optimize as the document keeps it, by its type: its file is the
    // one its requirements are written in (a solve may share its name).
    let file = v.doc.picked::<Optimize>("optimize", Some(&o.name))?.file.clone();
    v.motion_of(&file)?;
    o.search_t = v.motion.map(|_| searched_at);
    let asked: Vec<(Claim, Option<u32>)> = {
        let opt = v.doc.pick::<Optimize>("optimize", Some(&o.name))?;
        opt.requires.iter().map(|(_, c, b)| (v.doc.claim_at(c, *b, &v.ctx), *b)).collect()
    };
    for (r, (c, bound)) in o.requires.iter_mut().zip(asked) {
        // A trial, as every design of the run is: not stored.
        let (check, _) = cache::trial(|| claims::asserted(&mut v, c, bound, None, (None, None), &[], &[], ClaimOptions::default()))
            .map_err(|e| InkError { message: format!("(require {}): {}", r.claim, e.message), ..e })?;
        r.answer = Some(check);
    }
    o.feasible = o.requires.iter().all(Requirement::satisfied);
    // The answers' params by the names of the view that asked them: t is
    // one where a claim was asked over a motion with t ranged.
    o.answer_params = v.doc.params.iter().map(|p| p.name.clone()).collect();
    Ok(o)
}

/// The run's warnings: a design returned that is not feasible names each
/// requirement it does not meet. Over a motion, those whose search penalty
/// was 0 at the one moment the search measures are named together, saying
/// so: the penalty decides nothing (a strict bound at its bound has none),
/// so it says only that the search was not steered by them. One that is its
/// start design says so, and that a feasible one may be the optimum, or the
/// method found no way from it.
fn ended(mut o: Optimized, lo: &[f64], hi: &[f64]) -> Optimized {
    let unmet: Vec<&Requirement> = o.requires.iter().filter(|r| r.answer.is_some() && !r.satisfied()).collect();
    let (unseen, seen): (Vec<&Requirement>, Vec<&Requirement>) = unmet.iter().partition(|r| o.search_t.is_some() && r.violation <= 0.0);
    let named = |rs: &[&Requirement]| rs.iter().map(|r| format!("(require {})", r.claim)).collect::<Vec<_>>().join(", ");
    if let (Some(t), false) = (o.search_t, unseen.is_empty()) {
        o.warnings.push(format!(
            "the design returned does not meet {}: over the span its `answer` does not hold, but the search's penalty for it was 0 at the one moment \
             the search measures (t {t}), so the search was not steered by it",
            named(&unseen)
        ));
    }
    if !seen.is_empty() {
        o.warnings.push(format!("the design returned does not meet {}: its `answer` says why", named(&seen)));
    }
    let moved = o.design.iter().zip(&o.start).zip(lo.iter().zip(hi)).any(|((d, s), (l, h))| (d - s).abs() > 1e-9 * (h - l));
    if !moved && !o.start.is_empty() {
        let how = match o.trace {
            Trace::Simplex { evaluations, .. } => format!("none of the {evaluations} designs Nelder-Mead evaluated improved on it"),
            Trace::Gradient { iterations, .. } => format!("in {iterations} iterations the method of moving asymptotes found no step from it that improves it"),
        };
        o.warnings.push(if o.feasible {
            format!(
                "the design returned is the start design (the params' values): {how} within the requirements. It may be the optimum, or the run found no way from it: \
                 try another start or a larger budget"
            )
        } else if !unseen.is_empty() && seen.is_empty() {
            format!("the design returned is the start design (the params' values): {how}")
        } else {
            format!("the design returned is the start design (the params' values): {how}. The run found no way from it: try another start or a larger budget")
        });
    }
    o
}

/// The method of moving asymptotes with exact derivatives where the kernel
/// has them (design_gradient), as it has at the start, and central
/// differences for any it has not at a later design; each iteration is one
/// design evaluation plus two per param differenced.
fn optimize_gradient(
    run: Optimized,
    idx: &[usize],
    params: &[String],
    budget: usize,
    lo: &[f64],
    hi: &[f64],
    x0: &[f64],
    design: &dyn Fn(&[f64]) -> Result<ModelView, InkError>,
    replay: &Replay,
) -> Result<Optimized, InkError> {
    let (name, sign) = (&run.name, if run.maximize { -1.0 } else { 1.0 });
    let mut log: Vec<(Vec<f64>, Measures, usize)> = Vec::new();
    let (mut exact, mut differenced) = (0usize, 0usize);
    // MMA treats a constraint within MARGIN of 0 as met (minimize_mma), and
    // converges onto an active one from either side: it is given each
    // requirement tightened by MARGIN, so a design it accepts is within the
    // bound by the search's measure (then `verify` asks each claim).
    const MARGIN: f64 = 1e-6;
    let trials = Trials::new(&run, params, "gradient", design, replay);
    let out = optimize::minimize_mma(lo, hi, x0, budget, |x| {
        let (mut d, m) = trials.evaluate(x, name).map_err(|e| e.message)?;
        let (mut g0, mut gs) = cache::trial(|| design_gradient(&mut d, name, idx)).map_err(|e| e.message)?;
        // Central differences for whatever has no exact derivative.
        let mut solves = 1;
        for j in 0..idx.len() {
            if g0[j].is_some() && gs.iter().all(|g| g[j].is_some()) {
                exact += 1;
                continue;
            }
            differenced += 1;
            let h = 0.01 * (hi[j] - lo[j]);
            let at = |step: f64| -> Result<(f64, Vec<f64>), String> {
                let mut xx = x.to_vec();
                xx[j] = (x[j] + step).clamp(lo[j], hi[j]);
                let (_, m) = trials.evaluate(&xx, name).map_err(|e| e.message)?;
                Ok((m.objective, m.requires.iter().map(|c| c.violation).collect()))
            };
            let (xp, xm) = ((x[j] + h).min(hi[j]), (x[j] - h).max(lo[j]));
            let ((vp, cp), (vm, cm)) = (at(xp - x[j])?, at(xm - x[j])?);
            solves += 2;
            let dx = (xp - xm).max(1e-300);
            g0[j].get_or_insert((vp - vm) / dx);
            for (i, g) in gs.iter_mut().enumerate() {
                g[j].get_or_insert((cp[i] - cm[i]) / dx);
            }
        }
        let constraints = m.requires.iter().map(|c| c.violation + MARGIN).collect();
        let objective = sign * m.objective;
        log.push((x.to_vec(), m, solves));
        Ok(optimize::Smooth {
            objective,
            gradient: g0.iter().map(|g| sign * g.unwrap()).collect(),
            constraints,
            jacobian: gs.iter().map(|g| g.iter().map(|x| x.unwrap()).collect()).collect(),
        })
    })
    // A design that could not be evaluated stops the method: its error, as
    // the model gave it (where, hint).
    .map_err(|m| trials.failed.borrow().last().map_or_else(|| InkError::msg(m), |(_, e)| e.clone()))?;
    let (best_x, _, iterations) = out;
    // The best design that meets every bound as written, else MMA's (the
    // least violating).
    let pick = log
        .iter()
        .filter(|l| l.1.within())
        .min_by(|a, b| (sign * a.1.objective).total_cmp(&(sign * b.1.objective)))
        .or_else(|| log.iter().find(|l| l.0 == best_x))
        .unwrap_or(&log[0]);
    let (design, found) = (pick.0.clone(), pick.1.clone());
    let (start, start_objective, start_within) = (log[0].0.clone(), log[0].1.objective, log[0].1.within());
    let history = log
        .into_iter()
        .map(|(x, m, evaluations)| Iteration { design: x, objective: m.objective, max_violation: m.requires.iter().fold(0.0f64, |w, c| w.max(c.violation)), evaluations })
        .collect();
    let failed = trials.failed.take();
    let mut warnings = found.warnings;
    warnings.extend(failed_warning(params, &failed));
    Ok(Optimized {
        design,
        objective: found.objective,
        feasible: false,
        warnings,
        failed,
        requires: found.requires,
        start,
        start_objective,
        start_within,
        trace: Trace::Gradient { iterations: iterations.len(), history, exact, differenced },
        ..run
    })
}
