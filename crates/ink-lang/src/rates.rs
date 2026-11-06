//! Derivatives through solves. A solve's unknowns x are a function of what
//! its equations F(x, v) = 0 read: along an input v (t, or a param), J
//! dx/dv = -dF/dv, J = dF/dx at the root (the implicit function theorem).
//! Graph::diff holds every param (a partial); Doc::derivative adds, for
//! each solved param read, its partial times its rate. Those rates are a
//! solve of their own (`Rates`): its equations are F's derivative along v
//! (Doc::derivative, x read as its rates), so linear in them with J for
//! coefficients; and its own rates are a solve of the same kind. So every
//! order is the chain rule taken exactly, solved at each moment with the
//! values. Where the root is not shown regular (a linkage at a dead
//! point), the rates are not defined: their solve fails there, and what
//! reads them is refused at that moment (ModelView::refusal), never given
//! a number.

use crate::eval::{Ctx, eval_dual};
use crate::graph::{NodeId, Var};
use crate::tape::Tape;
use crate::{Doc, ExtForm, Form, Param, ParamProvider, Reads, Solve};
use std::sync::Arc;

/// The rates of a solve's unknowns along one input: params this solves
/// for, from equations linear in them (the module's).
pub struct Rates {
    /// The solve's name, for messages: errors at it name it.
    name: String,
    /// The input, as written: "t" or a param's name.
    along: String,
    owns: Vec<usize>,
    equations: Vec<NodeId>,
    uses_t: bool,
    /// The equations whose root these rest on (index in Doc::ext): the
    /// solve the first rates are of, whose root must be regular.
    root: usize,
}

impl Rates {
    /// What it is, for messages.
    pub(crate) fn label(&self) -> String {
        format!("the rates of solve '{}' along {}", self.name, self.along)
    }
}

impl Form for Rates {
    fn name(&self) -> &str {
        &self.name
    }

    fn equations(&self) -> Option<&[NodeId]> {
        Some(&self.equations)
    }

    fn provider(&self, doc: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(Linear {
            name: self.name.clone(),
            along: self.along.clone(),
            owns: self.owns.clone(),
            tapes: self.equations.iter().map(|&e| Tape::new(&doc.g, e)).collect(),
            uses_t: self.uses_t,
            root: doc.ext[self.root].form.provider(doc).expect("a solve provides params"),
        }))
    }
}

/// Rates compiled: least squares of equations linear in them, at the root
/// `root` set, where it is regular.
struct Linear {
    name: String,
    along: String,
    owns: Vec<usize>,
    tapes: Vec<Tape>,
    uses_t: bool,
    root: Arc<dyn ParamProvider>,
}

impl ParamProvider for Linear {
    fn name(&self) -> &str {
        &self.name
    }

    fn owns(&self) -> &[usize] {
        &self.owns
    }

    fn uses_t(&self) -> bool {
        self.uses_t
    }

    fn warm_starts(&self) -> bool {
        false
    }

    /// J y = -G(0) at the root certified regular (ParamProvider::root),
    /// J and G(0) exact (forward-mode duals of the equations at y = 0), by
    /// Householder QR.
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, String> {
        let undefined = |why: String| format!("the rates of solve '{}' along {} are not defined at t = {t}: {why}", self.name, self.along);
        let root = self.root.root(params, t).map_err(undefined)?;
        let mut at = params.to_vec();
        for (&k, v) in self.root.owns().iter().zip(root) {
            at[k] = v;
        }
        for &k in &self.owns {
            at[k] = 0.0;
        }
        let ctx = Ctx::new(t, at);
        let (mut jac, mut rhs) = (Vec::new(), Vec::new());
        for tp in &self.tapes {
            let d = eval_dual(tp, &ctx, [0.0; 3]);
            rhs.push(-d[0]);
            jac.push(self.owns.iter().map(|&k| d[5 + k]).collect());
        }
        let y = least_squares(jac, rhs).ok_or_else(|| undefined("its Jacobian in its unknowns has no full rank there".into()))?;
        for (&k, v) in self.owns.iter().zip(y) {
            params[k] = v;
        }
        Ok(None)
    }
}

/// The least-squares solution of a y = b, a m x n with m >= n (Householder
/// QR); None where a has not full column rank in floats.
fn least_squares(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let (m, n) = (a.len(), a.first().map_or(0, Vec::len));
    for c in 0..n {
        let norm = (c..m).map(|r| a[r][c] * a[r][c]).sum::<f64>().sqrt();
        if !(norm > 0.0 && norm.is_finite()) {
            return None;
        }
        let alpha = if a[c][c] > 0.0 { -norm } else { norm };
        let mut v: Vec<f64> = (c..m).map(|r| a[r][c]).collect();
        v[0] -= alpha;
        let vv: f64 = v.iter().map(|x| x * x).sum();
        for k in c..n {
            let s = 2.0 * (c..m).map(|r| v[r - c] * a[r][k]).sum::<f64>() / vv;
            for r in c..m {
                a[r][k] -= s * v[r - c];
            }
        }
        let s = 2.0 * (c..m).map(|r| v[r - c] * b[r]).sum::<f64>() / vv;
        for r in c..m {
            b[r] -= s * v[r - c];
        }
    }
    let mut y = vec![0.0; n];
    for r in (0..n).rev() {
        y[r] = (b[r] - (r + 1..n).map(|k| a[r][k] * y[k]).sum::<f64>()) / a[r][r];
    }
    y.iter().all(|v| v.is_finite()).then_some(y)
}

impl Doc {
    /// The derivative of `root` along `var`: its partial (every param held,
    /// Graph::derivative) plus, for each solved param it reads that varies
    /// along `var`, its partial along that param times the param's rate
    /// (`rate`): the chain rule through solves.
    pub fn derivative(&mut self, root: NodeId, var: Var) -> NodeId {
        let mut d = self.g.derivative(root, var);
        if matches!(var, Var::X | Var::Y | Var::Z) {
            return d;
        }
        for k in self.params_read(&[root]) {
            if let Some(r) = self.rate(k, var) {
                let dk = self.g.derivative(root, Var::Param(k as u32));
                let term = self.g.mul(dk, r);
                d = self.g.add(d, term);
            }
        }
        d
    }

    /// The joints' speed-limit claims made with their rates through every
    /// solve their drives read (`rated`): what a value no solve sets as
    /// they were read, a solve read after them sets.
    pub(super) fn rate_joints(&mut self) {
        for (i, drive, lo, hi) in std::mem::take(&mut self.rated) {
            let rate = self.derivative(drive, Var::T);
            let (a, b) = (self.g.sub(rate, lo), self.g.sub(hi, rate));
            let claim = self.g.min(a, b);
            let x = &mut self.asserts[i];
            x.claim = crate::Claim::Value(claim, crate::Cmp::Ge, 0.0);
            x.terms[0].1 = rate;
        }
    }

    /// Param k's rate along `var` (a param its rates' solve sets, Rates),
    /// where a solve of equations sets k and k varies along `var`; None
    /// where it does not vary. What it rests on that no solve sets is held
    /// in the derivative: each such param is kept with where (`held`), and
    /// a solve read later that makes it vary along `var` is refused
    /// (`keep_form`).
    pub fn rate(&mut self, k: usize, var: Var) -> Option<NodeId> {
        if var == Var::Param(k as u32) {
            return None;
        }
        if let Some(&r) = self.rates_of.get(&(k, var)) {
            return r;
        }
        let Some(i) = self.solves.iter().position(|s| s.owns.contains(&k)) else {
            let site = self.form_at.clone();
            self.held.entry((k, var)).or_insert(site);
            return None;
        };
        let s = self.solves[i].clone();
        let eqs = self.ext[s.ext].form.equations().map(<[NodeId]>::to_vec)?;
        // Along one of its own unknowns the others are held; along what
        // it does not rest on, nothing varies.
        let under = self.solves_under(&s.owns);
        let varies = match var {
            Var::T => self.solves[i].on_t,
            Var::Param(p) => !s.owns.contains(&(p as usize)) && under.iter().any(|&j| self.solves[j].reads.contains(&(p as usize))),
            _ => false,
        };
        if !varies {
            let site = self.form_at.clone();
            let free: Vec<usize> = under.iter().flat_map(|&j| self.solves[j].reads.clone()).filter(|p| !self.solves.iter().any(|s| s.owns.contains(p))).collect();
            for p in free {
                self.held.entry((p, var)).or_insert_with(|| site.clone());
            }
            for &x in &s.owns {
                self.rates_of.insert((x, var), None);
            }
            return None;
        }
        let along = match var {
            Var::Param(p) => self.params[p as usize].name.clone(),
            _ => "t".to_string(),
        };
        let mut owns = Vec::new();
        for &x in &s.owns {
            owns.push(self.params.len());
            self.params.push(Param { name: format!("(diff {} {along})", self.params[x].name), default: 0.0, tol: None });
            let y = self.g.var(Var::Param(owns[owns.len() - 1] as u32));
            self.rates_of.insert((x, var), Some(y));
        }
        // F's derivative along var, its unknowns read as their rates; the
        // other solves' unknowns it reads get theirs first.
        let equations: Vec<NodeId> = eqs.iter().map(|&f| self.derivative(f, var)).collect();
        let (base, root) = (&self.ext[s.ext], s.rates.unwrap_or(s.ext));
        // They rest on the root's equations too: whether it is regular.
        let mut fields = equations.clone();
        fields.extend(self.ext[root].form.equations().expect("a solve of equations"));
        let uses_t = fields.iter().any(|&e| self.g.depends_on(e, |v| v == Var::T));
        let form = Rates { name: base.form.name().to_string(), along, owns: owns.clone(), equations, uses_t, root };
        let (package, head, name, file, line, col) = (base.package, base.head, base.name.clone(), base.file.clone(), base.line, base.col);
        let reads = Reads { fields, forms: Vec::new() };
        let ext = self.ext.len();
        let read = self.provider_reads(&reads, &owns);
        self.ext.push(ExtForm { package, head, name, form: Arc::new(form), reads, file, line, col });
        for &y in &owns {
            self.solved_params.insert(y);
        }
        let mut solves = self.solves.clone();
        solves.push(Solve { ext, owns, reads: read, on_t: uses_t, on_demand: false, warm_starts: false, rates: Some(root) });
        self.set_solves(self.order_solves(solves, None).expect("rates read only what is solved before them"));
        self.solve_at_load();
        self.rates_of[&(k, var)]
    }

    /// Whether a failure of solve `i` at a moment leaves the model's shapes
    /// without a value there: any solve but rates, and rates a body's
    /// shape or placement reads (an instance's, or one `placed` records).
    /// Other rates are read by claims, efforts and derivatives alone, each
    /// refusing where they have none.
    pub fn needed(&self, i: usize) -> bool {
        if self.solves[i].rates.is_none() {
            return true;
        }
        let entries = |x: &crate::Xf| x.r.into_iter().flatten().chain(x.t);
        let instances = self.assemblies.iter().flat_map(|a| &a.instances).flat_map(|n| entries(&n.at));
        let placed = self.placed.values().flat_map(|p| std::iter::once(p.shape).chain(entries(&p.at)));
        let shapes: Vec<NodeId> = self.bodies.iter().map(|b| b.shape).chain(instances).chain(placed).collect();
        self.solves_read(&shapes).contains(&i)
    }

    /// Refuse what the solves read so far make vary that a derivative
    /// held fixed (`held`): a held param a solve now sets, varying along
    /// the derivative's input through everything that solve rests on, so a
    /// solve that leaves it fixed passes the hold on to what it reads (u =
    /// 3 m, then m = 100 t). Before the solve just read none varied, so
    /// one that does now is that solve's doing. Err: which, and where.
    pub(super) fn varies_held(&self) -> Result<(), String> {
        let mut held: Vec<(&(usize, Var), &String)> = self.held.iter().collect();
        held.sort_by(|a, b| (a.0.0, a.1).cmp(&(b.0.0, b.1)));
        for (&(k, var), site) in held {
            let Some(i) = self.solves.iter().position(|s| s.owns.contains(&k)) else { continue };
            let varies = match var {
                Var::T => self.solves[i].on_t,
                Var::Param(p) => self.solves_under(&[k]).iter().any(|&j| self.solves[j].reads.contains(&(p as usize))),
                _ => false,
            };
            if varies {
                let along = match var {
                    Var::Param(p) => format!("'{}'", self.params[p as usize].name),
                    _ => "t".into(),
                };
                return Err(format!("it makes '{}' vary along {along}, which the derivative at {site} held fixed (no solve set it there)", self.params[k].name));
            }
        }
        Ok(())
    }
}
