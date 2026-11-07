//! A loaded model, viewed at one set of params and one moment of time
//! (claims over a motion step through it: `motion`): the document, its context, and what analysis has learned of its
//! fields so far (where each one is, its proven bounds, the forms run). The
//! claim runner (`claims`) and every command work on one.

use crate::InkError;
use crate::{Doc, FormMetric, Given, Item, Measured, Outcome, Reach, Registry, RunOpts};
use ink_core::V3;
use ink_core::eval::Ctx;
use ink_core::graph::{Node, Var};
use ink_core::tape::Tape;
use ink_kernel::query::max_extent;
use ink_kernel::ranged::{self, Proof, ProofOptions};
use ink_kernel::rigid::Rigid;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::sync::Arc;

/// The most cells across a shape `ModelView::components` meshes it at, or
/// counts its voids down to, as for any mesh a query makes at a res asked
/// for.
pub const COMPONENT_CELLS: f64 = 20000.0;
/// The most triangles `ModelView::components` predicts a mesh made again
/// to sample: `ink render`'s and `ink components`' budget (25 M, about 65
/// CPU s on 6 vCPU of a Ryzen AI Max+ 395).
pub const COMPONENT_TRIANGLES: f64 = 25e6;

/// A shape's shells as `ink components` reads them: its mesh's (an
/// estimate, at cells `res`, with what it still leaves out), and its
/// enclosed voids, counted by proof as the `cavities` claim counts them.
#[derive(Clone, Debug)]
pub struct Components {
    pub analysis: ink_kernel::mesh::Analysis,
    pub unresolved: ink_kernel::mesh::Unresolved,
    pub voids: ink_kernel::voids::Voids,
    /// The mesh's cells.
    pub res: f64,
}

/// A contact zone excuses every overlap inside it: one not proven bounded
/// over the whole scope (every moment of the motion, read as one interval
/// of t, and every value of the ranged params) could hide a clash anywhere.
/// A relation's zone is the intersection of a field in each part's own
/// coordinates: bounded where either is.
fn zones_bounded(doc: &Doc, ctx: &Ctx) -> Result<(), InkError> {
    let over = if ctx.t_span.is_some() || !ctx.p_spans.is_empty() { " over the motion and param ranges" } else { "" };
    for (asm, r) in doc.assemblies.iter().flat_map(|a| a.relations.iter().map(move |r| (&a.name, r))) {
        let Some(z) = r.zone else { continue };
        if !doc.reach(z.a, ctx).whole && !z.b.is_some_and(|b| doc.reach(b, ctx).whole) {
            return Err(InkError {
                hint: Some("clip the zone to where the parts engage, e.g. (max (shape-of screw thread) (slab-z -10 0))".into()),
                ..InkError::msg(format!("assembly '{asm}': the contact zone of '{}' and '{}' is not proven bounded{over}", r.a, r.b))
            });
        }
    }
    Ok(())
}

/// The scope claims are answered over, as `ink assert` asks them: the
/// motion claims that read t are checked over (`ModelView::motion`), and
/// the params ranged, each over its span (by
/// index: set_param_spans). The command line's rules make it (ink-cli
/// session::claim_scope); `set_claim_scope` sets it.
#[derive(Clone, Debug, Default)]
pub struct ClaimScope {
    pub motion: Option<(f64, f64)>,
    /// Each claim over the span of the file it is written in, and an
    /// assembly's pairs over its file's (Doc::span_of), not `motion`: no
    /// span or moment was given the run.
    pub by_file: bool,
    pub spans: Vec<(usize, (f64, f64))>,
}

pub struct ModelView {
    pub doc: Doc,
    pub ctx: Ctx,
    /// The file the document was loaded from, if any, and the param
    /// overrides it was loaded with.
    pub path: Option<PathBuf>,
    pub overrides: Vec<(String, Given)>,
    /// The motion claims that read t are checked over, t ranged over it
    /// (crate::motion): `--t-span`, else the file's `(t-span A B)`.
    pub motion: Option<(f64, f64)>,
    /// Where `motion` is that of the file a claim or assembly is written in
    /// (`motion_of`), the moment the view was scoped at: a claim of a file
    /// stating no span is checked there.
    pub(crate) by_file: Option<f64>,
    /// What the view found, each answer kept by the scope it was found
    /// over (`Scope`), so that any other scope finds it again: per field,
    /// where analysis looks (Doc::reach).
    reaches: HashMap<(u32, Scope), Reach>,
    /// Per field: proven bounds (None: empty), inside its reach.
    boxes: HashMap<(u32, Scope), Option<(V3, V3, f64)>>,
    /// Per field: its tape clipped to its reach (`clipped`): claims about
    /// one part of a big assembly read the assembly's tape again and again.
    tapes: HashMap<(u32, Scope), Arc<Tape>>,
    /// Per extension form (package, name): what running it gave, and how,
    /// by the scope of the fields it depends on (Doc::form_inputs).
    outcomes: HashMap<((&'static str, String), Scope), (RunOpts, Rc<dyn Outcome>)>,
    /// What a query found of a field (`per_body`), by the query, the
    /// field, the other fields it read and its res.
    found: HashMap<(Found, Scope), Rc<dyn std::any::Any>>,
    /// Per field: the inputs beyond space it reads (t: u32::MAX; params by
    /// index), sorted: what its answers' scope is made of.
    inputs: HashMap<u32, Rc<[u32]>>,
    /// The solves that do not solve at this view's moment (Doc::ctx_at),
    /// by index in `Doc::solves`, with why: nothing resting on one is
    /// decided here (`refusal`).
    pub unsolved: Vec<(usize, InkError)>,
    /// Per field: whether it is the same at every moment (reads no t,
    /// directly or through a solve that does), once found.
    timeless: HashMap<u32, bool>,
    /// The param claims over a motion range t as (motion::ranged), once
    /// made (`t_param`).
    t_param: Option<u32>,
    /// Searches that may replace answers given already, for their caller
    /// to run once every answer is given (motion::Upgrade: interfere).
    pub upgrades: Vec<crate::motion::Upgrade>,
}

impl ModelView {
    /// A view of `doc` at `ctx`.
    pub fn new(doc: Doc, ctx: Ctx) -> ModelView {
        ModelView {
            doc,
            ctx,
            path: None,
            overrides: Vec::new(),
            motion: None,
            by_file: None,
            reaches: HashMap::new(),
            boxes: HashMap::new(),
            tapes: HashMap::new(),
            outcomes: HashMap::new(),
            found: HashMap::new(),
            inputs: HashMap::new(),
            unsolved: Vec::new(),
            timeless: HashMap::new(),
            t_param: None,
            upgrades: Vec::new(),
        }
    }

    /// Load a file with the packages of `registry` (the standard
    /// distribution's: `ink::registry()`), its params overridden, at moment
    /// `t`. A solve that does not solve at t leaves the view, with what
    /// rests on it refused (`unsolved`, `refusal`). Err: a direction with
    /// no length at t (Doc::check_conditions).
    pub fn load_in(registry: &'static Registry, path: &Path, overrides: &[(String, Given)], t: f64) -> Result<ModelView, InkError> {
        let doc = Doc::load_file_in(registry, path, overrides)?;
        let (ctx, unsolved) = doc.ctx_at(overrides, t)?;
        doc.check_conditions(&ctx, None, &unsolved.iter().map(|u| u.0).collect::<Vec<_>>(), false)?;
        zones_bounded(&doc, &ctx)?;
        Ok(ModelView { path: Some(path.to_path_buf()), overrides: overrides.to_vec(), unsolved, ..ModelView::new(doc, ctx) })
    }

    /// The same file loaded again, with the same packages, moment and
    /// solver, its overrides changed by `changes` (a design of the
    /// optimizer, a param stepped for a difference). A change that says no
    /// unit keeps the one its param was given in: the number is in the
    /// language's unit of the same quantity. Analysis starts over.
    /// Err: a solve that does not solve there.
    pub fn reload(&self, changes: &[(String, Given)]) -> Result<ModelView, InkError> {
        let path = self.path.as_ref().ok_or_else(|| InkError::msg("only a model loaded from a file can be reloaded with other params"))?;
        let mut overrides = self.overrides.clone();
        for (k, g) in changes {
            match overrides.iter_mut().find(|(n, _)| n == k) {
                Some((_, was)) => *was = Given { value: g.value, unit: g.unit.or(was.unit) },
                None => overrides.push((k.clone(), *g)),
            }
        }
        let v = ModelView::load_in(self.doc.registry, path, &overrides, self.ctx.t)?;
        if let Some((_, e)) = v.unsolved.first() {
            return Err(e.clone());
        }
        Ok(v)
    }

    /// Provide, in this view's context, the solves on demand that the
    /// fields at `roots` rest on (a published mass's density a command
    /// weighs: Doc::want), and from now on at every moment. Err: one
    /// that cannot be provided, at its form.
    pub fn want(&mut self, roots: &[u32]) -> Result<(), InkError> {
        self.doc.want(roots);
        // The solves what `roots` rest on, through what each reads: its
        // provider's own account where it gives one (ParamProvider::reads),
        // else what it was read from.
        let (mut todo, mut under) = (self.doc.params_read(roots), std::collections::BTreeMap::new());
        while let Some(k) = todo.pop() {
            for (i, s) in self.doc.solves.iter().enumerate() {
                if s.owns.contains(&k) && !under.contains_key(&i) {
                    let reads = self.doc.ext[s.ext].form.provider(&self.doc).and_then(|p| p.reads()).unwrap_or_else(|| s.reads.clone());
                    todo.extend(&reads);
                    under.insert(i, reads);
                }
            }
        }
        for (i, reads) in under {
            let s = &self.doc.solves[i];
            if !s.on_demand {
                continue;
            }
            // Provided at one value of each param, unless derived over the
            // scope (ModelView::scope_spans): over a ranged param it reads,
            // that value says nothing of the rest of the range.
            if let Some(&k) = reads.iter().find(|&&k| self.ctx.p_spans.get(k as u32).is_some())
                && !s.owns.iter().all(|&o| self.ctx.p_spans.is_derived(o as u32))
            {
                let what: Vec<String> = s.owns.iter().map(|&o| format!("'{}'", self.doc.params[o].name)).collect();
                return Err(InkError {
                    hint: Some(self.doc.at_one_value(&[k], &self.ctx.p_spans)),
                    ..self.doc.solve_error(
                        i,
                        format!(
                            "{} is found at one value of each param, but reads ranged '{}'. That value says nothing of the rest of its range",
                            what.join(", "),
                            self.doc.params[k].name
                        ),
                    )
                });
            }
            // One that cannot be provided here is unsolved at this moment,
            // as one that does not solve is: a claim stepping over a motion
            // reads it at this moment too (set_moment keeps it).
            if let Err(m) = self.doc.provide_one(i, &mut self.ctx.params, self.ctx.t) {
                let e = self.doc.solve_error(i, m);
                self.unsolved.push((i, e.clone()));
                return Err(e);
            }
        }
        Ok(())
    }

    /// Every solve on demand (a published mass) that would fail in this
    /// view's context, as `ink check` asks (ParamProvider::refusal: a mass
    /// is not measured to say so), each by the params it sets, with the
    /// error at its form. What they read is provided first.
    pub fn mass_refusals(&mut self) -> Vec<(Vec<usize>, InkError)> {
        let demand: Vec<usize> = (0..self.doc.solves.len()).filter(|&i| self.doc.solves[i].on_demand).collect();
        let reads: Vec<usize> = demand.iter().flat_map(|&i| self.doc.solves[i].reads.clone()).collect();
        self.doc.want_params(&reads);
        let mut failed = Vec::new();
        for i in demand {
            // What it reads that is on demand (an inner assembly's mass),
            // provided; one that fails says so at its own form.
            let read = self.doc.solves[i].reads.clone();
            for j in self.doc.solves_under(&read).into_iter().filter(|&j| self.doc.solves[j].on_demand) {
                let _ = self.doc.provide_one(j, &mut self.ctx.params, self.ctx.t);
            }
            let p = self.doc.ext[self.doc.solves[i].ext].form.provider(&self.doc).expect("a solve provides params");
            if let Some(m) = p.refusal(&self.ctx.params, self.ctx.t) {
                failed.push((self.doc.solves[i].owns.clone(), self.doc.solve_error(i, m)));
            }
        }
        failed
    }

    /// The solves that do not solve at the view's moment (by index in
    /// Doc::solves).
    fn unsolved_solves(&self) -> Vec<usize> {
        self.unsolved.iter().map(|u| u.0).collect()
    }

    /// The params with no value at the view's moment: set by a solve that
    /// does not solve there, or by one resting on such a solve.
    pub fn unsolved_params(&self) -> Vec<usize> {
        let failed = self.unsolved_solves();
        let d = &self.doc;
        (0..d.solves.len()).filter(|&i| d.solves_under(&d.solves[i].owns).iter().any(|j| failed.contains(j))).flat_map(|i| d.solves[i].owns.clone()).collect()
    }

    /// The chain rule's terms through solves for the fields at `roots`:
    /// each solved param they read, each input of `along` it varies along,
    /// and the param its rate there is (Doc::rate), solved at this view's
    /// moment (one not defined there is unsolved: `refusal`). A derivative
    /// along an input is its partial plus each term's partial along the
    /// solved param times the rate: what eval's gradient and the shape
    /// derivatives add to theirs.
    pub fn rates(&mut self, roots: &[u32], along: &[Var]) -> Vec<(usize, Var, usize)> {
        let mut out = Vec::new();
        let solved: Vec<usize> = self.doc.params_read(roots).into_iter().filter(|k| self.doc.solved_params.contains(k)).collect();
        for k in solved {
            for &var in along {
                if let Some(r) = self.doc.rate(k, var)
                    && let crate::graph::Node::Var(Var::Param(y)) = self.doc.g.nodes[r as usize]
                {
                    out.push((k, var, y as usize));
                }
            }
        }
        // Rates made here: their solves are the last (Doc::rate), each
        // provided at this moment.
        let known = self.ctx.params.len();
        if self.doc.params.len() > known {
            self.ctx.params.extend(self.doc.params[known..].iter().map(|p| p.default));
            for i in 0..self.doc.solves.len() {
                if self.doc.solves[i].owns.iter().any(|&o| o >= known)
                    && let Err(m) = self.doc.provide_one(i, &mut self.ctx.params, self.ctx.t)
                {
                    self.unsolved.push((i, self.doc.solve_error(i, m)));
                }
            }
        }
        out
    }

    /// Why what the fields at `roots` give cannot be decided in this view,
    /// if it cannot: they rest on a solve that does not solve at the view's
    /// moment.
    pub fn refusal(&self, roots: &[u32]) -> Option<InkError> {
        if self.unsolved.is_empty() {
            return None;
        }
        let solves = self.doc.solves_read(roots);
        self.unsolved.iter().find(|(i, _)| solves.contains(i)).map(|(_, e)| e.clone())
    }

    /// Check claims that read t at every moment of `span` (None: at the
    /// view's moment), with t ranged over it (crate::motion). The view's moment is one of
    /// the span's (its start, where it lay outside it): what does not move
    /// over the span is checked there. Err: a direction not proven to have
    /// a length over the span (Doc::check_conditions).
    pub fn set_motion(&mut self, span: Option<(f64, f64)>) -> Result<(), InkError> {
        if let Some((a, b)) = span
            && !(a < b && a.is_finite() && b.is_finite())
        {
            return Err(InkError::msg("a time span needs finite a < b"));
        }
        if let Some((a, b)) = span
            && !(a <= self.ctx.t && self.ctx.t <= b)
        {
            self.set_moment(a)?;
        }
        if span.is_some() {
            self.doc.check_conditions(&self.ctx, span, &self.unsolved_solves(), false)?;
            zones_bounded(&self.doc, &Ctx { t_span: span, ..self.ctx.clone() })?;
            // Made now, so the answers' param names hold it.
            self.t_param();
        }
        self.motion = span;
        Ok(())
    }

    /// The param a claim over a motion reads for t (motion::ranged): "t",
    /// added to the document's params the first time, at the view's moment.
    pub fn t_param(&mut self) -> u32 {
        if let Some(k) = self.t_param {
            return k;
        }
        self.doc.params.push(crate::Param { name: "t".into(), default: self.ctx.t, tol: None });
        let k = self.doc.params.len() - 1;
        if self.ctx.params.len() <= k {
            self.ctx.params.resize(k + 1, f64::NAN);
        }
        self.ctx.params[k] = self.ctx.t;
        self.t_param = Some(k as u32);
        k as u32
    }

    /// Analyse moment t: solves that read t solved again there, and what
    /// was learned of fields that read t forgotten. Err: a solve that does
    /// not solve at t, where params are ranged (its value over them is not
    /// derived again).
    pub fn set_moment(&mut self, t: f64) -> Result<(), InkError> {
        if self.ctx.t == t {
            return Ok(());
        }
        // A solve on t is provided again here, where it is provided at all
        // (one on demand, while something reads it).
        if (0..self.doc.solves.len()).any(|i| self.doc.solves[i].on_t && self.doc.provided_always(i)) {
            if !self.ctx.p_spans.is_empty() {
                let ks: Vec<usize> = self.ctx.p_spans.iter().filter(|p| !p.derived).map(|p| p.k as usize).collect();
                return Err(InkError::msg(format!("a solve that reads t is not solved again over ranged params at each step: {}", self.doc.at_one_value(&ks, &self.ctx.p_spans))));
            }
            let (ctx, unsolved) = self.doc.ctx_at(&self.overrides, t)?;
            (self.ctx, self.unsolved) = (ctx, unsolved);
            self.outcomes.clear();
        }
        self.ctx.t = t;
        // What reads t, directly or through a solve that does (solved
        // again here), is forgotten.
        let on_t: Vec<u32> = (0..self.doc.solves.len()).filter(|&i| self.doc.solves[i].on_t).flat_map(|i| self.doc.solves[i].owns.iter().map(|&k| k as u32)).collect();
        let g = &self.doc.g;
        let memo = std::cell::RefCell::new(std::mem::take(&mut self.timeless));
        let still = |f: &u32| *memo.borrow_mut().entry(*f).or_insert_with(|| !g.depends_on(*f, |v| v == Var::T || matches!(v, Var::Param(k) if on_t.contains(&k))));
        // (Each answer is kept by its scope, which another moment never
        // matches: this is what may not be asked for again.)
        self.reaches.retain(|k, _| still(&k.0));
        self.boxes.retain(|k, _| still(&k.0));
        self.tapes.retain(|k, _| still(&k.0));
        self.found.retain(|k, _| still(&k.0.1) && k.0.2.iter().all(still));
        let doc = &self.doc;
        self.outcomes.retain(|k, _| doc.ext_form(Some(k.0.0), &k.0.1).is_none_or(|e| doc.form_inputs(e).iter().all(still)));
        self.timeless = memo.into_inner();
        Ok(())
    }

    /// The model's accuracy (twice Doc::finest, at most 0.05 mm): what a
    /// shape claim is checked to unless it states its own.
    pub fn accuracy(&self) -> f64 {
        (2.0 * self.finest()).min(0.05)
    }

    /// The params' declared tolerance bands around their values in this
    /// view (`(param d 8 (tol 0 0.022))`: 8.000..8.022), by index.
    pub fn declared_spans(&self) -> Result<Vec<(usize, (f64, f64))>, InkError> {
        let mut out = Vec::new();
        for (k, p) in self.doc.params.iter().enumerate().filter(|(_, p)| p.tol.is_some()) {
            let band = p.band(self.ctx.params[k]).ok_or_else(|| InkError::msg(format!("'{}': its tolerance band around {} overflows", p.name, self.ctx.params[k])))?;
            out.push((k, band));
        }
        Ok(out)
    }

    /// Why param k cannot be ranged, if it cannot (see set_param_spans).
    pub fn cannot_range(&self, k: usize) -> Option<InkError> {
        let name = &self.doc.params[k].name;
        if let Some(site) = self.doc.materialized(k) {
            let (what, or) = match self.doc.ranged.iter().any(|r| r.0 == k) {
                true => ("(ranged lo hi)", ""),
                false => ("ranged (a tolerance)", ", remove its (tol lo hi), or check at the values with --nominal"),
            };
            return Some(InkError {
                hint: Some(format!("make that use symbolic (a field of the param, not a number at load time){or}")),
                ..InkError::msg(format!("'{name}' is {what}, but its value is used at load time at {site}, where it becomes one fixed number. No proof there covers its range"))
            });
        }
        let providers = self.doc.providers();
        if let Some(p) = providers.owner(k) {
            return Some(InkError::msg(format!("'{name}' is solved by '{}', so it has no range of its own", p.name())));
        }
        None
    }

    /// The ranged params these fields read, with their spans (sorted): a
    /// solved param derived over the scope stands for the params it is
    /// derived from.
    pub fn ranged_reads(&self, nodes: &[u32]) -> Vec<(u32, (f64, f64))> {
        self.ranged_reads_in(&self.ctx.p_spans, nodes)
    }

    /// `ranged_reads` over the scope `spans` (ModelView::scope_spans).
    pub fn ranged_reads_in(&self, spans: &ink_core::eval::ParamSpans, nodes: &[u32]) -> Vec<(u32, (f64, f64))> {
        if spans.is_empty() {
            return Vec::new();
        }
        let mut ks: Vec<u32> = self
            .doc
            .g
            .reachable(nodes)
            .into_iter()
            .filter_map(|n| match self.doc.g.nodes[n as usize] {
                ink_core::graph::Node::Var(Var::Param(k)) if spans.get(k).is_some() => Some(k),
                _ => None,
            })
            .flat_map(|k| spans.sources(k))
            .collect();
        ks.sort_unstable();
        ks.dedup();
        ks.into_iter().filter_map(|k| Some((k, spans.get(k)?))).collect()
    }

    /// This view scoped as every claim's question, where it is not yet (a
    /// pointwise view: eval's, a library caller's): each (ranged lo hi)
    /// over its range (scope_spans). The claims' entry points call it, so
    /// none is decided at one value of one.
    pub fn scope_question(&mut self) -> Result<(), InkError> {
        if self.doc.ranged.iter().all(|&(k, r)| self.ctx.p_spans.get(k as u32) == Some(r)) {
            return Ok(());
        }
        let given: Vec<(usize, (f64, f64))> = self.ctx.p_spans.iter().filter(|p| !p.derived).map(|p| (p.k as usize, p.span)).collect();
        self.set_param_spans(&given)
    }

    /// The solved params this view derives over its ranged params as
    /// certified roots: claims reading them are proven, not measured, where
    /// the solve reads no other solved value that is not (Doc::inexact).
    /// A root of equations over sampled fields is the samples', so what
    /// reads it is measured anyway.
    pub fn certified(&self) -> Vec<u32> {
        self.certified_in(&self.ctx.p_spans)
    }

    /// `certified` over the scope `spans` (ModelView::scope_spans).
    pub fn certified_in(&self, spans: &ink_core::eval::ParamSpans) -> Vec<u32> {
        let mut out = Vec::new();
        for s in &self.doc.solves {
            let own: Vec<u32> = out.iter().copied().chain(s.owns.iter().map(|&u| u as u32)).collect();
            if s.owns.iter().all(|&u| spans.is_derived(u as u32)) && self.doc.inexact(&self.doc.form_inputs(&self.doc.ext[s.ext]), &own).solved.is_empty() {
                out.extend(s.owns.iter().map(|&u| u as u32));
            }
        }
        out
    }

    /// The scope of these ranged params (by index): their spans, and what
    /// each solve that reads one of them sets, derived over them as a
    /// certified root (ParamProvider::over_scope), each after the solves
    /// setting what it reads (Doc::solves: so it is certified over their
    /// enclosures, wherever they are in the file). Refused, with an error
    /// naming why (never a value at one point labelled as covering the
    /// range): a param whose value became a number at load time (a loop
    /// count, a sampled or measured field:
    /// Doc::materialized), a param a solve sets, a solve over them that
    /// cannot be certified, and a direction not proven to have a length
    /// over them (a rotation's axis: Doc::check_conditions).
    pub fn scope_spans(&self, spans: &[(usize, (f64, f64))]) -> Result<ink_core::eval::ParamSpans, InkError> {
        // Each (ranged lo hi) is in every scope, over all of its range: no
        // span given drops or narrows it.
        let spans: Vec<(usize, (f64, f64))> = spans.iter().filter(|s| !self.doc.ranged.iter().any(|r| r.0 == s.0)).chain(&self.doc.ranged).copied().collect();
        let spans = &spans[..];
        for &(k, (lo, hi)) in spans {
            if !(lo <= hi && lo.is_finite() && hi.is_finite()) {
                return Err(InkError::msg(format!("'{}': a range needs finite lo <= hi, got {lo}:{hi}", self.doc.params[k].name)));
            }
            if let Some(e) = self.cannot_range(k) {
                return Err(e);
            }
        }
        let mut ctx = Ctx { p_spans: ink_core::eval::ParamSpans::new(spans.iter().map(|&(k, s)| (k as u32, s))), ..self.ctx.clone() };
        let doc = &self.doc;
        // Over the solves provided here: not one on demand nothing reads
        // (a published mass a claim over the range does not weigh).
        // Rates have no proven range: what reads them is refused claim by
        // claim (set_param_spans).
        for (_, s) in doc.solves.iter().enumerate().filter(|&(i, s)| doc.provided_always(i) && s.rates.is_none()) {
            let e = &doc.ext[s.ext];
            let Some(p) = e.form.provider(doc) else { continue };
            // What it read (Doc::form_inputs): its equations...
            let ranged: Vec<usize> = s.reads.iter().copied().filter(|&k| ctx.p_spans.get(k as u32).is_some()).collect();
            if ranged.is_empty() {
                continue;
            }
            let names = ranged.iter().map(|&k| format!("'{}'", doc.params[k].name)).collect::<Vec<_>>().join(", ");
            let what = format!("{} '{}' depends on ranged {names}", e.head, e.form.name());
            let hint = Some(doc.at_one_value(&ranged, &ctx.p_spans));
            for &u in p.owns() {
                if let Some(site) = doc.materialized(u) {
                    return Err(InkError {
                        hint,
                        ..InkError::msg(format!(
                            "{what}, but what it solves, '{}', is used at load time at {site}, where it is one number. No proof there covers the range",
                            doc.params[u].name
                        ))
                    });
                }
            }
            let o = p.over_scope(&ctx).map_err(|m| InkError { hint: hint.clone(), ..InkError::msg(format!("{what}, and what it solves has no proven range: {m}")) })?;
            ctx.p_spans = ctx.p_spans.derived(o.derive, &o.bands);
        }
        if !ctx.p_spans.is_empty() {
            doc.check_conditions(&ctx, self.motion, &self.unsolved_solves(), false)?;
            zones_bounded(doc, &ctx)?;
        }
        Ok(ctx.p_spans)
    }

    /// Analyse every value of each of these params in its span (by index;
    /// none: every param at its value): proofs cover them all, a witness
    /// says which values fail. What a solve over them sets is derived as a
    /// certified root; what cannot be covered is refused (scope_spans).
    pub fn set_param_spans(&mut self, spans: &[(usize, (f64, f64))]) -> Result<(), InkError> {
        self.ctx.p_spans = self.scope_spans(spans)?;
        // Rates over ranged params are found at single values only: what
        // reads them is refused (`refusal`), the rest decided.
        for i in (0..self.doc.solves.len()).filter(|&i| self.doc.solves[i].rates.is_some()) {
            let s = &self.doc.solves[i];
            let ks: Vec<usize> = s.reads.iter().copied().filter(|&k| self.ctx.p_spans.get(k as u32).is_some()).collect();
            if !ks.is_empty() && !self.unsolved.iter().any(|u| u.0 == i) {
                let rates: Vec<&str> = s.owns.iter().map(|&k| self.doc.params[k].name.as_str()).collect();
                let ranged: Vec<String> = ks.iter().map(|&k| format!("'{}'", self.doc.params[k].name)).collect();
                let e = InkError {
                    hint: Some(self.doc.at_one_value(&ks, &self.ctx.p_spans)),
                    ..self.doc.solve_error(
                        i,
                        format!(
                            "the rates {} of solve '{}' read ranged {}, but are found at single values only, with nothing proven over the ranges",
                            rates.join(", "),
                            self.doc.ext[s.ext].form.name(),
                            ranged.join(", ")
                        ),
                    )
                };
                self.unsolved.push((i, e));
            }
        }
        self.reaches.clear();
        self.boxes.clear();
        self.tapes.clear();
        self.outcomes.clear();
        self.found.clear();
        Ok(())
    }

    /// Answer claims over `scope` (ClaimScope): its motion, its
    /// params' spans.
    pub fn set_claim_scope(&mut self, scope: &ClaimScope) -> Result<(), InkError> {
        self.by_file = scope.by_file.then_some(self.ctx.t);
        // Made now, where an import's span may be checked over, so the
        // answers' param names hold it (set_motion).
        if self.by_file.is_some() && !self.doc.file_spans.is_empty() {
            self.t_param();
        }
        self.set_motion(scope.motion)?;
        self.set_param_spans(&scope.spans)
    }

    /// The moment what reads t is answered at when no span covers it and
    /// the run asked for no moment or span (claims::one_moment); None:
    /// over a motion, or at a moment the run asked for.
    pub fn unspanned_moment(&self) -> Option<f64> {
        self.by_file.filter(|_| self.motion.is_none())
    }

    /// Scope the view for what is written in `file` (a claim's, an
    /// assembly's): over its file's span where the run gave none
    /// (ClaimScope::by_file), else as the run asked. A claim keeps the
    /// motion it was written over when its file is imported, and a static
    /// file's claims stay at the run's moment.
    pub fn motion_of(&mut self, file: &str) -> Result<(), InkError> {
        self.motion_over(self.doc.span_of(file))
    }

    /// Scope the view over `span` where the run gave none (`motion_of`).
    pub fn motion_over(&mut self, span: Option<(f64, f64)>) -> Result<(), InkError> {
        let Some(t) = self.by_file else { return Ok(()) };
        if span == self.motion {
            return Ok(());
        }
        // From the run's moment, as the run file's claims were (set_motion).
        self.set_moment(t)?;
        self.set_motion(span)
    }

    /// What running an extension form (of `package`, if given) gave: run
    /// once per view and options.
    pub fn outcome(&mut self, package: Option<&str>, form: &str, opts: RunOpts) -> Result<Rc<dyn Outcome>, InkError> {
        let e = self.doc.ext_form(package, form).ok_or_else(|| InkError::msg(format!("no package form named '{form}'")))?;
        let (key, f, inputs) = ((e.package, form.to_string()), e.form.clone(), self.doc.form_inputs(e));
        let key = (key, self.scope(&inputs));
        if let Some((o, out)) = self.outcomes.get(&key)
            && *o == opts
        {
            return Ok(out.clone());
        }
        let out: Rc<dyn Outcome> = f.run(&mut self.doc, &self.ctx, &opts)?.into();
        self.outcomes.insert(key, (opts, out.clone()));
        Ok(out)
    }

    /// A package's metric of one of its forms or of a shape: the one path
    /// claims, requirements and objectives on package metrics are
    /// evaluated by. A failed run's metric has no value (`Outcome::failure`).
    pub fn measure(&mut self, m: &FormMetric, opts: RunOpts) -> Result<Measured, InkError> {
        if let Some(f) = m.shape {
            let p = self.doc.registry.packages.iter().find(|p| p.name() == m.package).expect("a metric's package");
            return p.shape_metric(self, m.metric, f, opts.res);
        }
        if let Some(why) = self.doc.ext_form(Some(m.package), &m.form).and_then(|e| e.form.unavailable(m.metric)) {
            return Ok(Measured { value: Err(why), at: None, res: None, warnings: Vec::new(), excluded: Vec::new(), unseen: None, proven: false });
        }
        let out = self.outcome(Some(m.package), &m.form, opts)?;
        let mut measured = out.metric(m.metric).ok_or_else(|| InkError::msg(format!("'{}' gave no '{}'", m.form, m.metric)))?;
        if let Some(why) = out.failure() {
            measured.value = Err(why);
        }
        Ok(measured)
    }

    /// Where analysis of field `f` looks: a box holding all of its
    /// material, or the window (see Doc::reach).
    pub fn reach(&mut self, f: u32) -> Reach {
        let key = (f, self.scope(&[f]));
        if let Some(r) = self.reaches.get(&key) {
            return *r;
        }
        let r = self.doc.reach(f, &self.ctx);
        self.reaches.insert(key, r);
        r
    }

    /// `reach` of each of `fs`, found at once, in parallel.
    pub fn reach_all(&mut self, fs: &[u32]) {
        use rayon::prelude::*;
        let keys: Vec<(u32, Scope)> = fs.iter().map(|&f| (f, self.scope(&[f]))).collect();
        let mut todo: Vec<(u32, Scope)> = keys.into_iter().filter(|k| !self.reaches.contains_key(k)).collect();
        todo.sort_unstable_by_key(|k| k.0);
        todo.dedup_by_key(|k| k.0);
        let (g, window, ctx) = (&self.doc.g, self.doc.domain(), &self.ctx);
        let found: Vec<Reach> = todo.par_iter().map(|k| crate::reach_in(g, window, k.0, ctx)).collect();
        self.reaches.extend(todo.into_iter().zip(found));
    }

    /// `bounds` of each of `fs`, found at once, in parallel.
    pub fn bounds_all(&mut self, fs: &[u32]) {
        use rayon::prelude::*;
        self.reach_all(fs);
        let keys: Vec<(u32, Scope)> = fs.iter().map(|&f| (f, self.scope(&[f]))).collect();
        let mut todo: Vec<(u32, Scope)> = keys.into_iter().filter(|k| !self.boxes.contains_key(k)).collect();
        todo.sort_unstable_by_key(|k| k.0);
        todo.dedup_by_key(|k| k.0);
        let asked: Vec<(Arc<Tape>, Reach)> = todo.iter().map(|k| (self.clipped(k.0), self.reach(k.0))).collect();
        let ctx = &self.ctx;
        let found: Vec<Option<(V3, V3, f64)>> = asked.par_iter().map(|(t, r)| ink_kernel::query::bounds(t, ctx, r.lo, r.hi)).collect();
        self.boxes.extend(todo.into_iter().zip(found));
    }

    /// `f` with the view at `at` (a point of its scope, or a scope wider
    /// than it: t over a span, a box over a motion). What the view found is
    /// kept by the scope it was found over (`Scope`), so `f` reads only
    /// what holds in `at` and finds the rest again; what it finds is
    /// forgotten after, the view's own caches back as they were, and so are
    /// the solves unsolved and the motion and steps (what `f` set:
    /// scope_question, set_claim_scope).
    pub fn at_point<R>(&mut self, at: Ctx, f: impl FnOnce(&mut Self) -> R) -> R {
        let whole = std::mem::replace(&mut self.ctx, at);
        let kept = (self.reaches.clone(), self.boxes.clone(), self.tapes.clone(), self.found.clone(), self.outcomes.clone());
        let scope = (self.unsolved.clone(), self.motion);
        let out = f(self);
        (self.ctx, (self.reaches, self.boxes, self.tapes, self.found, self.outcomes)) = (whole, kept);
        (self.unsolved, self.motion) = scope;
        out
    }

    /// The scope fields `fs` read, as the view's context has it (`Scope`).
    pub(crate) fn scope(&mut self, fs: &[u32]) -> Scope {
        let mut reads: Vec<u32> = Vec::new();
        for &f in fs {
            let r = match self.inputs.get(&f) {
                Some(r) => r.clone(),
                None => {
                    let g = &self.doc.g;
                    let mut r: Vec<u32> = g
                        .reachable(&[f])
                        .into_iter()
                        .filter_map(|n| match g.nodes[n as usize] {
                            Node::Var(Var::T) => Some(u32::MAX),
                            Node::Var(Var::Param(k)) => Some(k),
                            _ => None,
                        })
                        .collect();
                    // A solved param is what its solve found from what that
                    // read, and its value and span do not say how those
                    // correlate it with the others (q = p and r = 1 - p or
                    // r = p, by a switch k: the same spans of q and r for k in
                    // [0, 1] and at 1, where q - r is 0). So its scope is its
                    // solve's too, transitively.
                    let ks: Vec<usize> = r.iter().filter(|&&k| k != u32::MAX).map(|&k| k as usize).collect();
                    for i in self.doc.solves_under(&ks) {
                        let s = &self.doc.solves[i];
                        r.extend(s.reads.iter().map(|&k| k as u32));
                        if s.on_t {
                            r.push(u32::MAX);
                        }
                    }
                    r.sort_unstable();
                    r.dedup();
                    let r: Rc<[u32]> = r.into();
                    self.inputs.insert(f, r.clone());
                    r
                }
            };
            reads.extend(r.iter());
        }
        reads.sort_unstable();
        reads.dedup();
        Scope::of(&self.ctx, &reads)
    }

    /// Whether analysis of `f` sees all of it (it is proven bounded).
    pub fn whole(&mut self, f: u32) -> bool {
        self.reach(f).whole
    }

    /// Tape of a field clipped to its reach.
    pub fn clipped(&mut self, f: u32) -> Arc<Tape> {
        let key = (f, self.scope(&[f]));
        if let Some(t) = self.tapes.get(&key) {
            return t.clone();
        }
        let r = self.reach(f);
        let c = self.doc.clip(f, r.lo, r.hi);
        let t = Arc::new(Tape::new(&self.doc.g, c));
        self.tapes.insert(key, t.clone());
        t
    }

    /// Proven bounds of a field (None: empty); of the part in the window if
    /// it is not proven bounded (see `whole`).
    pub fn bounds(&mut self, f: u32) -> Option<(V3, V3, f64)> {
        let key = (f, self.scope(&[f]));
        if let Some(b) = self.boxes.get(&key) {
            return *b;
        }
        let t = self.clipped(f);
        let Reach { lo, hi, .. } = self.reach(f);
        let b = ink_kernel::query::bounds(&t, &self.ctx, lo, hi);
        self.boxes.insert(key, b);
        b
    }

    /// A proof to accuracy `tol` (mm; 0: exactly): at `res`, else coarse
    /// to fine down to half the accuracy (an exact one: the model's floor);
    /// ranged params' spans halved up to PARAM_SPLITS times.
    pub fn proof_options(&self, res: Option<f64>, tol: f64) -> ProofOptions {
        let finest = if tol > 0.0 { 0.5 * tol } else { self.finest() };
        ProofOptions { res, finest, depth: ink_kernel::scope::PARAM_SPLITS, tol }
    }

    /// The finest cells automatic refinement halves to: the model's floor
    /// (Doc::finest).
    pub fn finest(&self) -> f64 {
        self.doc.finest()
    }

    /// Prove a field's solid holds no point deeper than `tol`
    /// (ranged::prove), and whether all of it was examined: a proof on a
    /// field not proven bounded holds in the window only (a witness there
    /// still disproves).
    pub fn prove(&mut self, f: u32, res: Option<f64>, tol: f64) -> (Proof, bool) {
        let opts = self.proof_options(res, tol);
        self.prove_with(f, &opts)
    }

    /// `prove` with these options.
    pub fn prove_with(&mut self, f: u32, opts: &ProofOptions) -> (Proof, bool) {
        let tape = self.clipped(f);
        let Reach { lo, hi, whole } = self.reach(f);
        (ranged::prove(&tape, &self.ctx, lo, hi, opts), whole)
    }

    /// Numeric placement of an instance where it is a rigid motion: it is
    /// fixed over the view's scope (one that varies over the time span or
    /// with a ranged param is not: a part placed at one value would hide the
    /// rest of its range), and its numbers are a turn and a move: they move
    /// no point of the body's box farther than `RIGID` of its size from
    /// where its world field, which reads them as R^T (p - t), has it
    /// (`moved`). Not a scale or a skew, whose moments, widths and surface
    /// no rigid motion maps. So too at each of `moments`, where the caller
    /// places the body by the numbers there as well (a GLB's motion
    /// samples): what it shows at each is then the world field's there.
    /// Else None: the world field is used.
    pub fn rigid(&mut self, it: &Item, moments: &[Ctx]) -> Option<(String, Rigid)> {
        let (body, xf) = it.instance.as_ref()?;
        self.rigid_at(body, xf, moments)
    }

    /// A placement's numbers where it is fixed over the view's scope (reads
    /// neither t over a span nor a ranged param): each entry's f64 value.
    pub fn numbers(&self, xf: &crate::Xf) -> Option<Rigid> {
        let num = |n: u32| -> Option<f64> {
            if let Some(c) = self.doc.g.as_const(n) {
                return Some(c);
            }
            if self.ctx.t_span.is_some() && self.doc.g.depends_on(n, |v| matches!(v, Var::T)) {
                return None;
            }
            if self.doc.g.depends_on(n, |v| matches!(v, Var::Param(k) if self.ctx.p_spans.get(k).is_some())) {
                return None;
            }
            Some(self.doc.scalar(n, &self.ctx))
        };
        let v: Vec<f64> = xf.r.iter().flatten().chain(&xf.t).map(|&n| num(n)).collect::<Option<_>>()?;
        Some(Rigid::from_entries(v.try_into().ok()?))
    }

    /// `rigid` for body `body` placed by `xf`.
    fn rigid_at(&mut self, body: &str, xf: &crate::Xf, moments: &[Ctx]) -> Option<(String, Rigid)> {
        let entries: Vec<u32> = xf.r.iter().flatten().copied().chain(xf.t).collect();
        let at = self.numbers(xf)?;
        let shape = self.doc.body(body)?.shape;
        let Some((lo, hi, _)) = self.bounds(shape) else { return Some((body.to_string(), at)) };
        let turn = |at: &Rigid| moved(at, lo, hi) <= RIGID * max_extent(lo, hi);
        if !turn(&at) {
            return None;
        }
        if !moments.is_empty() && crate::claims::reads_t(self, &entries) {
            let tapes: Vec<Tape> = entries.iter().map(|&n| Tape::new(&self.doc.g, n)).collect();
            for c in moments {
                let v: Vec<f64> = tapes.iter().map(|tp| ink_core::eval::eval_point(tp, c, [0.0; 3])).collect();
                if !turn(&Rigid::from_entries(v.try_into().ok()?)) {
                    return None;
                }
            }
        }
        Some((body.to_string(), at))
    }

    /// Proven world bounds of an item: once per body in its own frame, then
    /// moved (instancing: `local`); its placed field's otherwise. Of the part
    /// in the window if it is not proven bounded (see `item_whole`).
    pub fn item_bounds(&mut self, it: &Item) -> Option<(V3, V3, f64)> {
        match self.local(it) {
            Some(l) => {
                let (lo, hi, h) = self.bounds(l.shape)?;
                let (a, b) = l.world_box(lo, hi);
                Some((a, b, h))
            }
            None => {
                let f = self.doc.placed_field(it);
                self.bounds(f)
            }
        }
    }

    /// Whether an item is proven bounded (where it is modelled, if placed).
    pub fn item_whole(&mut self, it: &Item) -> bool {
        match self.local(it) {
            Some(_) => true,
            None => {
                let f = self.doc.placed_field(it);
                self.whole(f)
            }
        }
    }

    /// An instance as its body's field in the body's own frame and its
    /// placement's numbers, so that what a query finds of the body is found
    /// once (`per_body`) and moved to each instance as a rigid motion moves
    /// it: a box, samples and normals, moments, a mesh. Where the placement
    /// is a rigid motion over the view's scope (`rigid`), each of its
    /// entries is bounded there (one reading t is a real number at a
    /// moment, of which its f64 value is one rounding; a sampled one has its
    /// reading's error), it is proven invertible (det R bounded away from 0,
    /// and R^T R near I, read as rigid::turn_bounds reads a turn: otherwise
    /// even a bounded body can be unbounded in the world), and the body is
    /// proven bounded in its frame; else None: the world field is used.
    pub fn local(&mut self, it: &Item) -> Option<Local> {
        let (body, xf) = it.instance.as_ref()?;
        self.local_at(body, xf)
    }

    /// `local` for body `body` placed by `xf` (an item's placement, or
    /// where it lies in a frame of the assembly: interfere's shared frame).
    pub fn local_at(&mut self, body: &str, xf: &crate::Xf) -> Option<Local> {
        use ink_core::interval::{I, mag, sub_o};
        let (body, at) = self.rigid_at(body, xf, &[])?;
        let mut err = 0.0f64;
        for (&n, &x) in xf.r.iter().flatten().chain(&xf.t).zip(at.r.iter().flatten().chain(&at.t)) {
            if !x.is_finite() {
                return None;
            }
            if self.doc.g.as_const(n).is_none() {
                let e = ink_core::eval::eval_interval(&Tape::new(&self.doc.g, n), &self.ctx, [I::ENTIRE; 3], None);
                let d = mag(sub_o(e, I::point(x)));
                if !(e.lo.is_finite() && e.hi.is_finite() && d.is_finite()) {
                    return None;
                }
                err = err.max(d);
            }
        }
        let turn = ink_kernel::rigid::turn_bounds(&self.doc.g, &xf.r, &self.ctx)?;
        let det = ink_core::interval::un(ink_core::graph::UnOp::Abs, turn.det);
        let slack = turn.slack()?;
        let shape = self.doc.body(&body)?.shape;
        (det.lo > 0.0 && self.whole(shape)).then_some(Local { shape, at, det, slack, err })
    }

    /// Wall thickness of a field by rays at surface samples (see
    /// feature::thickness: point, normal, area, thickness; at least one),
    /// and their spacing. None: the field is empty. Err(spacing): no sample
    /// fell on its surface (UNSAMPLED), so nothing was measured.
    pub fn thickness(&mut self, f: u32, res: Option<f64>) -> Option<Result<(Vec<(V3, V3, f64, f64)>, f64), f64>> {
        let (lo, hi, _) = self.bounds(f)?;
        let tape = self.clipped(f);
        let cells = res.map_or(128.0, |h| max_extent(lo, hi) / h);
        let (samples, h) = (ink_kernel::feature::thickness(&tape, &self.ctx, lo, hi, cells), max_extent(lo, hi) / cells);
        // Not cached: the claim reading it is computed.
        ink_kernel::cache::computed();
        Some(if samples.is_empty() { Err(h) } else { Ok((samples, h)) })
    }

    /// A field's enclosed voids, by proof (voids::voids): its cells halved
    /// down to `res`, else to its box's longest side / COMPONENT_CELLS.
    /// None: the field is empty. More than COMPONENT_CELLS across at `res`
    /// is an error.
    pub fn voids(&mut self, f: u32, res: Option<f64>) -> Result<Option<ink_kernel::voids::Voids>, InkError> {
        let Some((lo, hi, _)) = self.bounds(f) else { return Ok(None) };
        let ext = max_extent(lo, hi);
        if let Some(h) = res
            && ext / h > COMPONENT_CELLS
        {
            return Err(InkError::msg(format!("res {h} gives {:.0} cells across the shape; the limit is {COMPONENT_CELLS}", ext / h)));
        }
        let tape = self.clipped(f);
        // Not cached: the claim reading it is computed.
        ink_kernel::cache::computed();
        Ok(Some(ink_kernel::voids::voids(&tape, &self.ctx, lo, hi, res.unwrap_or(ext / COMPONENT_CELLS))))
    }

    /// `ink components`: a field's shells (mesh::analyze: those whose
    /// triangles face inward are voids; an estimate), meshed at cells
    /// `res`, used as given; else at its box's longest side / 160, made
    /// again where it leaves material out (mesh::Unresolved: a wall or
    /// sheet thinner than the cells), each time at cells finer than the
    /// thickest of what it left, while at most COMPONENT_CELLS across and
    /// COMPONENT_TRIANGLES predicted; and its voids by proof (`voids`).
    /// None: the field is empty.
    pub fn components(&mut self, f: u32, res: Option<f64>) -> Result<Option<Components>, InkError> {
        let Some(voids) = self.voids(f, res)? else { return Ok(None) };
        let Some((lo, hi, _)) = self.bounds(f) else { return Ok(None) };
        let ext = max_extent(lo, hi);
        let tape = self.clipped(f);
        let mut h = res.unwrap_or(ext / 160.0);
        let (mut m, mut unresolved) = ink_kernel::mesh::mesh_checked(&tape, &self.ctx, lo, hi, h);
        while res.is_none() && unresolved.material > 0 {
            let k = ((h / unresolved.thickest).log2().floor().max(0.0) + 1.0).min(60.0);
            let finer = h / 2f64.powf(k);
            if ext / finer > COMPONENT_CELLS || (m.tris.len().max(1) as f64) * 4f64.powf(k) > COMPONENT_TRIANGLES {
                break;
            }
            h = finer;
            (m, unresolved) = ink_kernel::mesh::mesh_checked(&tape, &self.ctx, lo, hi, h);
        }
        Ok(Some(Components { analysis: ink_kernel::mesh::analyze(&m), unresolved, voids, res: h }))
    }

    /// Mass moments of a field (weighted by a density field, if any) and its
    /// volume, at the cells `res` asks for, else measured to the mass's
    /// accuracy (query::measured); found once per field and density
    /// (`per_body`). None: it has no box.
    pub fn moments(&mut self, shape: u32, rho: Option<u32>, res: Option<f64>) -> Option<ink_kernel::query::Measured> {
        let found = self.per_body("mass", shape, rho.as_slice(), res, |v| {
            let (lo, hi, _) = v.bounds(shape)?;
            let rho = rho.map(|f| Tape::new(&v.doc.g, f));
            Some(ink_kernel::query::measured(&v.clipped(shape), &[rho.as_ref(), None], true, &v.ctx, lo, hi, res))
        });
        (*found).clone()
    }

    /// What query `what` finds of field `shape` (with `with`, the other
    /// fields it reads: integrands, a density) at `res`: found once, however
    /// many instances place the field (`local`), and moved to each by the
    /// caller. One field is one node (the graph holds each expression once),
    /// so instances share only what is the same field exactly.
    pub fn per_body<T: 'static>(&mut self, what: &'static str, shape: u32, with: &[u32], res: Option<f64>, compute: impl FnOnce(&mut ModelView) -> T) -> Rc<T> {
        let fields: Vec<u32> = std::iter::once(shape).chain(with.iter().copied()).collect();
        let key = ((what, shape, with.to_vec(), res.map(f64::to_bits)), self.scope(&fields));
        if let Some(v) = self.found.get(&key).cloned().and_then(|v| v.downcast::<T>().ok()) {
            return v;
        }
        let v = Rc::new(compute(self));
        self.found.insert(key, v.clone());
        v
    }
}

/// What `ModelView::per_body` keeps a query's answer by: the query, the
/// field, the other fields it read and its res (bits).
type Found = (&'static str, u32, Vec<u32>, Option<u64>);

/// The scope an answer was found over, exactly, as the inputs beyond space
/// its fields read have it (and the inputs of the solves setting the params
/// they read, transitively: ModelView::scope): t's moment and span; each
/// param's value and, where it is ranged, its span, its band and whether it
/// is derived (a param ranged over [1, 1] reads its expressions in real
/// numbers, one at 1 as its floats: (p + 1e16) - 1e16 is 1 there, 0 at the
/// value). Two contexts with one scope for a field give it one value at
/// every point, so the view keeps each answer by it and finds it again
/// under any other.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Scope(Vec<u64>);

impl Scope {
    /// `ctx`'s scope for inputs `reads` (t: u32::MAX; params by index).
    fn of(ctx: &Ctx, reads: &[u32]) -> Scope {
        let mut key = Vec::with_capacity(4 * reads.len());
        for &r in reads {
            key.push(r as u64);
            if r == u32::MAX {
                key.push(ctx.t.to_bits());
                match ctx.t_span {
                    Some((a, b)) => key.extend([1, a.to_bits(), b.to_bits()]),
                    None => key.push(0),
                }
                continue;
            }
            key.push(ctx.params.get(r as usize).map_or(u64::MAX, |x| x.to_bits()));
            match ctx.p_spans.iter().find(|s| s.k == r) {
                Some(s) => key.extend([1, s.span.0.to_bits(), s.span.1.to_bits(), s.band.0.to_bits(), s.band.1.to_bits(), s.derived as u64]),
                None => key.push(0),
            }
        }
        Scope(key)
    }
}

/// How far, as a fraction of a body's size, an instance's numbers may move
/// a point of it from where its world field has it, for its placement to be
/// rigid (`ModelView::rigid`). A turn's numbers are off by rounding (1e-16
/// of the body's size and its distance from its origin); what queries
/// resolve is far coarser (a mesh's tol is at least a body's size over
/// 600000).
pub const RIGID: f64 = 1e-9;

/// Why a shape's thickness was not measured where no surface sample fell
/// on it (`ModelView::thickness`).
pub const UNSAMPLED: &str =
    "no surface sample fell on it at this spacing. Its walls or pieces are thinner than about the spacing, so nothing was measured. A finer res samples them";

/// Bound on |(R^T R)^-1 - I| for a placement's numbers R (as
/// TurnBounds::slack bounds it for its entries in real numbers), in floats:
/// an estimate, as the numbers are. Infinite where R^T R is not near I.
fn slack_of(at: &Rigid) -> f64 {
    let r = at.r;
    let s = (0..3).map(|i| (0..3).map(|j| ((0..3).map(|k| r[k][i] * r[k][j]).sum::<f64>() - if i == j { 1.0 } else { 0.0 }).abs()).sum::<f64>()).fold(0.0, f64::max);
    if s < 1.0 { s / (1.0 - s) } else { f64::INFINITY }
}

/// How far numbers `at` move a point y of the box [lo, hi] from where a
/// world field reading them as R^T (p - t) has it (t + R^-T y), as a
/// distance: R^-T y = R (y + d), |d| <= slack max |y| (`slack_of`), so at
/// most |R| slack max |y| in each coordinate, sqrt 3 (under 2) times that.
/// An estimate, as the numbers are (a proven output covers the entries in
/// real numbers: Local::world_box).
fn moved(at: &Rigid, lo: V3, hi: V3) -> f64 {
    let most = lo.into_iter().chain(hi).map(f64::abs).fold(0.0, f64::max);
    let r = at.r.iter().map(|row| row.iter().map(|x| x.abs()).sum::<f64>()).fold(0.0, f64::max);
    2.0 * r * slack_of(at) * most
}

/// An instance in its body's own frame (`ModelView::local`): its world
/// shape is the body's shape S at R^T (p - t) (`Doc::items`), so its
/// material is t + R^-T S.
#[derive(Clone, Copy, Debug)]
pub struct Local {
    /// The body's field, in its own frame.
    pub shape: u32,
    /// The placement's numbers (f64), by which an estimate is moved, as the
    /// world field's point values read them.
    pub at: Rigid,
    /// Holds |det R| over the view's scope (rigid::turn_bounds), bounded
    /// away from 0: a volume or integral in the world is the body's over
    /// it.
    pub det: ink_core::interval::I,
    /// Bound on |(R^T R)^-1 - I| (TurnBounds::slack): R^-T y = R (y + d),
    /// |d| <= slack max |y|.
    slack: f64,
    /// Bound on how far each entry of R and t is from its number in `at`.
    err: f64,
}

impl Local {
    /// A proven box of the instance from one of the body in its frame,
    /// rounded outward: widened for R^-T against R (`slack`), moved by the
    /// numbers, then widened for how far the entries may be from them
    /// (`err`: by |dt| + |dR| |y| <= err (1 + 3 max |y|)). An exact turn
    /// with exact entries keeps the box.
    pub fn world_box(&self, lo: V3, hi: V3) -> (V3, V3) {
        use ink_core::interval::{I, add_o, add_up, mul_up, sub_o};
        let widen = |(lo, hi): (V3, V3), w: f64| match w == 0.0 {
            true => (lo, hi),
            false => (lo.map(|x| sub_o(I::point(x), I::point(w)).lo), hi.map(|x| add_o(I::point(x), I::point(w)).hi)),
        };
        let most = |lo: V3, hi: V3| lo.into_iter().chain(hi).map(f64::abs).fold(0.0, f64::max);
        let (lo, hi) = widen((lo, hi), mul_up(self.slack, most(lo, hi)));
        let y = most(lo, hi);
        widen(self.at.aabb(lo, hi), mul_up(self.err, add_up(1.0, mul_up(3.0, y))))
    }

    /// How much farther than d a point within `reach` of the box [lo, hi],
    /// moved by the numbers, may be from the instance's surface as its world
    /// field has it, where it is d from the body's in its frame (d at most
    /// `reach`): `moved` there, and R^-T stretches d by at most
    /// sqrt(1 + slack). Rounding aside, as the world field's own samples.
    pub fn off(&self, lo: V3, hi: V3, reach: f64) -> f64 {
        let (lo, hi) = (lo.map(|x| x - reach), hi.map(|x| x + reach));
        moved(&self.at, lo, hi) + slack_of(&self.at) * reach
    }

    /// An integral's estimate over the instance from the body's in its
    /// frame: over |det R| (its middle).
    pub fn integral(&self, x: f64) -> f64 {
        x / self.det.mid()
    }
}
