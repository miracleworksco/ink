//! `(solve ...)` and `(optimize ...)`: the constraint solver (LM)
//! and design objectives.

use crate::InkError;
use crate::eval::Ctx;
use crate::graph::{Graph, NodeId, Var};
use crate::lang::*;
use crate::sexp::Sexp;
use crate::tape::Tape;

type R<T> = Result<T, InkError>;

/// A system of equations solved for some params (constraints, linkage
/// kinematics). Residuals may use t: they are re-solved per moment.
pub struct SolveSys {
    pub name: String,
    pub unknowns: Vec<usize>,
    pub residuals: Vec<NodeId>,
    pub uses_t: bool,
}

/// What an optimization measures.
#[derive(Clone, Debug)]
pub enum Objective {
    Volume(NodeId),
    Mass(NodeId, NodeId),
    /// A package's metric of one of its forms.
    Metric(FormMetric),
    /// A field that does not vary in space: any expression of parameters.
    Value(NodeId),
}

pub struct Optimize {
    pub name: String,
    pub vary: Vec<(String, f64, f64)>,
    pub maximize: bool,
    pub objective: Objective,
    /// Text, claim, and its bound as a field (see Doc::claim_at).
    pub requires: Vec<(String, Claim, Option<NodeId>)>,
    pub budget: usize,
}

/// The design package: `(solve ...)` systems and `(optimize ...)`.
pub struct Design;

impl Package for Design {
    fn name(&self) -> &'static str {
        "design"
    }

    fn forms(&self) -> &'static [&'static str] {
        &["optimize", "solve"]
    }

    fn commands(&self) -> &'static [ink_lang::command::Command] {
        &crate::command::COMMANDS
    }

    fn metrics(&self) -> &'static [(&'static str, &'static str)] {
        &[
            ("min-thickness", "the thinnest wall of a shape (mm), by rays from surface samples: (min-thickness part > t); of (air-around part m), the narrowest gap"),
            ("cavities", "how many enclosed voids a shape has (air no path leads out of) with air deeper than the accuracy, by proof: (cavities part < 1)"),
        ]
    }

    fn metric_quantity(&self, metric: &str) -> Option<ink_lang::units::Quantity> {
        (metric == "min-thickness").then_some(ink_lang::units::Quantity::Length)
    }

    fn shape_metrics(&self) -> &'static [&'static str] {
        &["min-thickness", "cavities"]
    }

    fn exclusions(&self) -> &'static [(&'static str, &'static str)] {
        &[("undecided-cells", "voids a count's finest cells could not rule out")]
    }

    fn shape_accuracy(&self, metric: &str) -> bool {
        metric == "cavities"
    }

    fn shape_metric(&self, v: &mut crate::lang::view::ModelView, metric: &str, f: NodeId, res: Option<f64>, accuracy: f64) -> R<Measured> {
        crate::driver::shape_metric(v, metric, f, res, accuracy)
    }

    fn read(&self, r: &mut Reader, form: &Sexp, head: &str, args: &[Sexp]) -> R<std::sync::Arc<dyn Form>> {
        Ok(match head {
            "solve" => std::sync::Arc::new(read_solve(r, form, args)?),
            _ => std::sync::Arc::new(read_optimize(r, form, args)?),
        })
    }
}

impl Form for SolveSys {
    fn name(&self) -> &str {
        &self.name
    }

    fn provider(&self, doc: &Doc) -> Option<std::sync::Arc<dyn ParamProvider>> {
        Some(std::sync::Arc::new(SolverSystem::new(&doc.g, self)))
    }

    fn equations(&self) -> Option<&[NodeId]> {
        Some(&self.residuals)
    }
}

impl Form for Optimize {
    fn name(&self) -> &str {
        &self.name
    }
}

/// A (solve ...) system compiled to tapes: a param provider that solves
/// at any time and param values without the document.
struct SolverSystem {
    name: String,
    unknowns: Vec<usize>,
    tapes: Vec<Tape>,
    uses_t: bool,
}

impl SolverSystem {
    fn new(g: &Graph, s: &SolveSys) -> SolverSystem {
        let tapes: Vec<Tape> = s.residuals.iter().map(|&r| Tape::new(g, r)).collect();
        SolverSystem { name: s.name.clone(), unknowns: s.unknowns.clone(), tapes, uses_t: s.uses_t }
    }
}

impl ParamProvider for SolverSystem {
    fn name(&self) -> &str {
        &self.name
    }

    fn owns(&self) -> &[usize] {
        &self.unknowns
    }

    fn uses_t(&self) -> bool {
        self.uses_t
    }

    /// Levenberg-Marquardt with exact derivatives, warm-started from the
    /// current values.
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, String> {
        let sys = self;
        let tapes = &sys.tapes;
        let eval = |params: &[f64]| -> (Vec<f64>, Vec<Vec<f64>>) {
            let ctx = Ctx::new(t, params.to_vec());
            let mut r = Vec::new();
            let mut jac = Vec::new();
            for tp in tapes {
                let d = crate::eval::eval_dual(tp, &ctx, [0.0; 3]);
                r.push(d[0]);
                jac.push(sys.unknowns.iter().map(|&k| d[5 + k]).collect());
            }
            (r, jac)
        };
        let n = sys.unknowns.len();
        let cost = |r: &[f64]| r.iter().map(|v| v * v).sum::<f64>();
        let (mut r, mut jac) = eval(params);
        let mut lambda = 1e-3;
        for _ in 0..200 {
            if r.iter().all(|v| v.abs() < 1e-12) {
                break;
            }
            // (J^T J + lambda diag) delta = -J^T r
            let mut a = vec![vec![0.0; n]; n];
            let mut g = vec![0.0; n];
            for (ri, row) in r.iter().zip(&jac) {
                for i in 0..n {
                    g[i] -= row[i] * ri;
                    for j in 0..n {
                        a[i][j] += row[i] * row[j];
                    }
                }
            }
            let mut stepped = false;
            // Damping floor: unknowns the equations leave free stay put instead of making the system singular.
            let floor = 1e-9 * (0..n).map(|i| a[i][i]).fold(0.0, f64::max).max(1e-300);
            for _ in 0..30 {
                let mut m = a.clone();
                for i in 0..n {
                    m[i][i] += lambda * a[i][i].max(floor);
                }
                let Some(delta) = solve_dense(m, g.clone()) else {
                    lambda *= 10.0;
                    continue;
                };
                let mut trial = params.to_vec();
                for (i, &k) in sys.unknowns.iter().enumerate() {
                    trial[k] += delta[i];
                }
                let (r2, j2) = eval(&trial);
                if cost(&r2) < cost(&r) {
                    params.copy_from_slice(&trial);
                    (r, jac) = (r2, j2);
                    lambda = (lambda / 3.0).max(1e-12);
                    stepped = true;
                    break;
                }
                lambda *= 4.0;
            }
            if !stepped {
                break;
            }
        }
        let worst = r.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        if !(worst < 1e-6) {
            return Err(format!("solve '{}' did not converge at t = {t} (largest residual {worst:.3e}); check the equations or the guesses", sys.name));
        }
        Ok(None)
    }

    fn root(&self, params: &[f64], t: f64) -> Result<Vec<f64>, String> {
        crate::lang::roots::Roots::regular(&self.unknowns, self.tapes.clone(), &Ctx::new(t, params.to_vec()))
    }

    fn follows(&self) -> bool {
        true
    }

    /// Its unknowns over ranged params: the one root of its equations near
    /// the solution at the ranges' middles, certified over all of them by
    /// Krawczyk's test (roots.rs), and narrowed on each piece of them.
    fn over_scope(&self, ctx: &Ctx) -> Result<OverScope, String> {
        let roots = crate::lang::roots::Roots::certify(&self.unknowns, self.tapes.clone(), self.uses_t, ctx, |p, t| self.provide(p, t).map(drop))?;
        let bands = roots.band.iter().zip(&self.unknowns).map(|(b, &k)| (k as u32, (b.lo, b.hi))).collect();
        Ok(OverScope { derive: std::sync::Arc::new(roots), bands })
    }
}

/// Solve a small dense linear system (Gaussian elimination, partial pivoting).
fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
        if a[p][c].abs() < 1e-300 {
            return None;
        }
        a.swap(c, p);
        b.swap(c, p);
        for r in c + 1..n {
            let f = a[r][c] / a[c][c];
            for k in c..n {
                a[r][k] -= f * a[c][k];
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        x[r] = (b[r] - (r + 1..n).map(|k| a[r][k] * x[k]).sum::<f64>()) / a[r][r];
    }
    Some(x)
}

/// (solve name (unknowns (a guess) ...) equation ...): each equation is
/// (== lhs rhs) or an expression that must be 0. Unknowns become params.
fn read_solve(r: &mut Reader, form: &Sexp, args: &[Sexp]) -> R<SolveSys> {
    let usage = "usage: (solve name (unknowns (a guess) (b guess) ...) (== lhs rhs) ... ); equations may use params and t, not x y z";
    let Some((name, rest)) = args.split_first() else { return Err(r.err_hint(form, "bad solve", Some(usage.into()))) };
    let n = r.name(name)?;
    let Some((unk, eqs)) = rest.split_first() else { return Err(r.err_hint(form, "solve needs unknowns and equations", Some(usage.into()))) };
    let ul = unk.list().filter(|l| l.first().and_then(|h| h.sym()) == Some("unknowns")).ok_or_else(|| r.err_hint(unk, "expected (unknowns (a guess) ...)", Some(usage.into())))?;
    let mut unknowns = Vec::new();
    let solved: Vec<usize> = r.providers().owned().collect();
    for u in &ul[1..] {
        // An existing param (declared earlier, e.g. used in placements)
        // becomes an unknown, its value the initial guess.
        if let Some(name) = u.sym() {
            let Some(k) = r.param_index(name) else {
                return Err(r.err_hint(u, format!("no param named '{name}'"), Some("unknowns are (name guess) for new values, or names of existing params".into())));
            };
            if solved.contains(&k) {
                return Err(r.err(u, format!("'{name}' is already solved by another system")));
            }
            unknowns.push(k);
            continue;
        }
        let Some([un, guess]) = u.list() else { return Err(r.err_hint(u, "unknowns are (name guess)", Some(usage.into()))) };
        let un_name = r.name(un)?;
        let g = r.number(guess)?;
        unknowns.push(r.new_param(un, &un_name, g)?);
    }
    let mut residuals = Vec::new();
    for e in eqs {
        let q = match e.list() {
            Some([h, a, b]) if h.sym() == Some("==") => {
                let (a, b) = (r.field(a)?, r.field(b)?);
                r.graph().sub(a, b)
            }
            _ => r.field(e)?,
        };
        if r.graph().depends_on(q, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
            return Err(r.err(e, "solve equations must not depend on x, y or z"));
        }
        residuals.push(q);
    }
    if residuals.len() < unknowns.len() {
        return Err(r.err_hint(form, format!("{} equations for {} unknowns: underdetermined", residuals.len(), unknowns.len()), Some(usage.into())));
    }
    // A solution is a numerical method's answer (Levenberg-Marquardt to a
    // tolerance; at a tangency it can sit on either side of the true
    // root), so solved params are numerical inputs: claims on them hold
    // numerically, never proven. Over ranged params they are certified
    // roots instead (SolverSystem::over_scope), and claims are proven.
    for &k in &unknowns {
        r.solved(k);
    }
    let uses_t = residuals.iter().any(|&q| r.graph().depends_on(q, |v| matches!(v, Var::T)));
    // The core solves it at load time, with every solve read so far in
    // their order, so later load-time code sees the solution (Doc::solves).
    Ok(SolveSys { name: n, unknowns, residuals, uses_t })
}

/// (optimize name (vary p lo hi) ... (minimize metric) | (maximize metric)
///   (require claim) ... [(budget n)])
fn read_optimize(r: &mut Reader, form: &Sexp, args: &[Sexp]) -> R<Optimize> {
    let usage = "usage: (optimize name (vary param lo hi) ... (minimize metric) (require claim) ... [(budget evaluations)]); \
                 metrics: (volume f) (mass f density) (value expr), a shape's (min-thickness f) (cavities f), or a package form's (metric form)";
    let Some((name, clauses)) = args.split_first() else { return Err(r.err_hint(form, "bad optimize", Some(usage.into()))) };
    let n = r.name(name)?;
    let (mut vary, mut objective, mut maximize, mut requires, mut budget) = (Vec::new(), None, false, Vec::new(), 60usize);
    for c in &r.clauses(clauses)? {
        match c.head_args().ok_or_else(|| r.err_hint(c, "optimize clauses are lists", Some(usage.into())))? {
            ("vary", [p, lo, hi]) => {
                let p = r.name(p)?;
                if let Some((_, (lo, hi))) = r.ranged_named(&p) {
                    return Err(r.err(c, format!("'{p}' is (ranged {lo} {hi}), a range of the question: optimize varies params")));
                }
                // Kept by its public name (`ns.p` in a file imported under a
                // namespace): the one the document sets it by.
                let Some(k) = r.param_index(&p) else { return Err(r.err(c, format!("'{p}' is not a param"))) };
                let p = r.param_names().swap_remove(k);
                if vary.iter().any(|(q, _, _)| *q == p) {
                    return Err(r.err(c, format!("'{p}' is varied twice")));
                }
                let (lo, hi) = (r.number(lo)?, r.number(hi)?);
                if lo >= hi {
                    return Err(r.err(c, "vary needs lo < hi"));
                }
                vary.push((p, lo, hi));
            }
            (sense @ ("minimize" | "maximize"), [m]) => {
                maximize = sense == "maximize";
                objective = Some(read_objective(r, m)?);
            }
            ("require", [claim]) => {
                let (cl, bound) = r.claim(claim)?;
                requires.push((claim.to_string(), cl, bound));
            }
            ("budget", [b]) => budget = r.number(b)?.max(4.0) as usize,
            _ => return Err(r.err_hint(c, "bad optimize clause", Some(usage.into()))),
        }
    }
    if vary.is_empty() {
        return Err(r.err_hint(form, "optimize needs at least one (vary param lo hi)", Some(usage.into())));
    }
    let objective = objective.ok_or_else(|| r.err_hint(form, "optimize needs (minimize metric) or (maximize metric)", Some(usage.into())))?;
    Ok(Optimize { name: n, vary, maximize, objective, requires, budget })
}

fn read_objective(r: &mut Reader, m: &Sexp) -> R<Objective> {
    Ok(match m.head_args().ok_or_else(|| r.err(m, "a metric is a list, e.g. (volume f)"))? {
        ("volume", [f]) => Objective::Volume(r.field(f)?),
        ("mass", [f, rho]) => {
            let (f, d) = (r.field(f)?, r.field(rho)?);
            if let Some(e) = crate::lang::materials::density_field_error(r.graph(), d) {
                return Err(r.err(rho, e));
            }
            Objective::Mass(f, d)
        }
        ("value", [e]) => {
            let f = r.field(e)?;
            if r.graph().depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
                return Err(r.err(e, "value must not depend on x, y or z"));
            }
            Objective::Value(f)
        }
        (k, [form]) if r.registry().metric(k).is_some() => Objective::Metric(r.metric(k, form, "optimize")?),
        _ => {
            let all = r.registry().metrics();
            return Err(r.err(m, format!("unknown metric: (volume f) (mass f density) (value expr), or a package's metric of a shape or a form, (metric of): {all}")));
        }
    })
}
