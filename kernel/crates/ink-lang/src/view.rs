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

/// The scope claims are answered over, as `ink require` asks them: the
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
    /// A reach found ahead of the question asking it (`reach_ahead`),
    /// taken once, by the `reach` asking it.
    ahead: Option<((u32, Scope), Reach)>,
    /// Per field: its tape clipped to its reach (`clipped`): claims about
    /// one part of a big assembly read the assembly's tape again and again.
    tapes: HashMap<(u32, Scope), Arc<Tape>>,
    /// Per extension form (package, name): what running it gave, and how,
    /// by the scope of the fields it depends on (Doc::form_inputs).
    outcomes: HashMap<((&'static str, String), Scope), (RunOpts, Rc<dyn Outcome>)>,
    /// The moment each form that reads the model's worst is run at (`at_worst`), once found.
    worst: HashMap<(&'static str, String), Option<(Ctx, Vec<String>)>>,
    /// And every moment read for it: its context, and each field's value there, in `worst_reads`' order.
    worst_all: HashMap<(&'static str, String), Vec<(Ctx, Vec<f64>)>>,
    /// The steps of the motion those moments do not cover: where what the form reads jumps, or changes faster
    /// than the step budget follows (the providers' `unchecked` and `loose`). A value there is no answer.
    worst_gaps: HashMap<(&'static str, String), Vec<(f64, f64)>>,
    /// `measure_over`'s answers, by form, metric and grid: every moment's, found once.
    over: HashMap<(&'static str, String, &'static str, Option<u64>, Option<u64>), Rc<Vec<Measured>>>,
    /// The moment `measure` answers for (`pin`).
    pinned: Option<usize>,
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
    /// Per field and params: the way it moves with each (`param_ways`).
    ways: HashMap<(u32, Vec<u32>), Vec<Option<i8>>>,
    /// Per field and substitution: the field substituted (`subst`).
    substs: HashMap<(u32, Vec<(Var, u32)>), u32>,
    /// The param claims over a motion range t as (motion::ranged), once
    /// made (`t_param`).
    t_param: Option<u32>,
    /// Searches that may replace answers given already, for their caller
    /// to run once every answer is given (motion::Upgrade: interfere).
    pub upgrades: Vec<crate::motion::Upgrade>,
    /// The boxes a search of this view may test (query::BOXES; fewer while what asks keeps a budget:
    /// motion::between).
    pub search_boxes: usize,
    /// Fields that read no z are looked at in the plane (`in_plane`).
    plane: bool,
}

/// A provider's param derived over the ranged params it reads (ParamProvider::over_ranged): every
/// piece or point of the scope a search or a sample makes is provided again there, so a value read
/// at one value of a published mass is that mass's. Where it gives nothing there the band it had
/// stands: it encloses every narrower piece's.
struct OverRanged {
    provider: std::sync::Arc<dyn crate::ParamProvider>,
    owns: [u32; 1],
    reads: Vec<u32>,
}

impl ink_core::eval::Derive for OverRanged {
    fn owns(&self) -> &[u32] {
        &self.owns
    }
    fn reads(&self) -> &[u32] {
        &self.reads
    }
    fn derive(&self, ctx: &mut Ctx) {
        let ranged: Vec<(usize, (f64, f64))> = self.reads.iter().filter_map(|&k| ctx.p_spans.get(k).map(|s| (k as usize, s))).collect();
        if let Some(Ok((value, band))) = self.provider.over_ranged(&ctx.params, ctx.t, &ranged) {
            ctx.p_spans.set_span(self.owns[0], band);
            ctx.params[self.owns[0] as usize] = value;
        }
    }
    fn rates(&self, _: &Ctx) -> Vec<Vec<f64>> {
        vec![vec![0.0; self.reads.len()]]
    }
}

/// A row left this small beside the largest, each column over its largest entry, is a sum of the rows
/// picked (`spanning`): the rounding of numbers the model worked out, far under what a solve resolves.
const SPANNED: f64 = 1e-9;

/// Of `rows` (each moment's loads, a column a component), a few that every row is a weighted sum of:
/// the rows picked, and each row's weights on them. As many are picked as the rows have independent
/// directions: the largest first, then whichever lies furthest from the sums of those picked, until
/// none lies further than `SPANNED`. With `one`, every row's weights also sum to 1 (each row is taken
/// with a 1 after it), for what the rows leave out and every moment holds once.
pub(crate) fn spanning(rows: &[&[f64]], one: bool) -> (Vec<usize>, Vec<Vec<f64>>) {
    let n = rows.first().map_or(0, |r| r.len());
    let most: Vec<f64> = (0..n).map(|j| rows.iter().map(|r| r[j].abs()).fold(0.0, f64::max)).collect();
    // What is left of each row beside the rows picked so far, and its part along each of those.
    let mut left: Vec<Vec<f64>> = rows.iter().map(|r| (0..n).filter(|&j| most[j] > 0.0).map(|j| r[j] / most[j]).chain(one.then_some(1.0)).collect()).collect();
    let mut along: Vec<Vec<f64>> = vec![Vec::new(); rows.len()];
    let mut picked: Vec<usize> = Vec::new();
    let size = |r: &[f64]| r.iter().map(|x| x * x).sum::<f64>().sqrt();
    let mut first = 0.0;
    while picked.len() < left.first().map_or(0, |r| r.len()) {
        let (k, far) = left.iter().map(|r| size(r)).enumerate().fold((0, 0.0), |b, x| if x.1 > b.1 { x } else { b });
        if far <= SPANNED * first || far == 0.0 {
            break;
        }
        first = first.max(far);
        let q: Vec<f64> = left[k].iter().map(|x| x / far).collect();
        for (r, a) in left.iter_mut().zip(&mut along) {
            let part: f64 = r.iter().zip(&q).map(|(x, y)| x * y).sum();
            r.iter_mut().zip(&q).for_each(|(x, y)| *x -= part * y);
            a.push(part);
        }
        picked.push(k);
    }
    // No load at any moment: one run answers for all.
    if picked.is_empty() {
        return (vec![0], vec![vec![1.0]; rows.len()]);
    }
    // A picked row's parts stop at its own direction: each row's weights by back-substitution.
    let weights = along
        .iter()
        .map(|a| {
            let mut w = vec![0.0; picked.len()];
            for j in (0..picked.len()).rev() {
                let above: f64 = (j + 1..picked.len()).map(|i| w[i] * along[picked[i]][j]).sum();
                w[j] = (a[j] - above) / along[picked[j]][j];
            }
            w
        })
        .collect();
    (picked, weights)
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
            ahead: None,
            boxes: HashMap::new(),
            tapes: HashMap::new(),
            outcomes: HashMap::new(),
            worst: HashMap::new(),
            worst_all: HashMap::new(),
            worst_gaps: HashMap::new(),
            over: HashMap::new(),
            pinned: None,
            found: HashMap::new(),
            inputs: HashMap::new(),
            unsolved: Vec::new(),
            timeless: HashMap::new(),
            ways: HashMap::new(),
            substs: HashMap::new(),
            t_param: None,
            upgrades: Vec::new(),
            search_boxes: ink_kernel::query::BOXES,
            plane: false,
        }
    }

    /// Load a file with the packages of `registry` (the standard
    /// distribution's: `ink::registry()`), its params overridden, at moment
    /// `t`. A solve that does not solve at t leaves the view, with what
    /// rests on it refused (`unsolved`, `refusal`). Err: a direction with
    /// no length at t (Doc::check_conditions).
    pub fn load_in(registry: &'static Registry, path: &Path, overrides: &[(String, Given)], t: f64) -> Result<ModelView, InkError> {
        Self::load_for(registry, path, overrides, t, false, &[])
    }

    /// `load_in`, its imports read as Doc::load_file_importing reads them:
    /// from `dirs` before `$INK_PATH`'s, here and on every reload.
    pub fn load_importing(registry: &'static Registry, path: &Path, overrides: &[(String, Given)], t: f64, dirs: &[std::path::PathBuf]) -> Result<ModelView, InkError> {
        Self::load_for(registry, path, overrides, t, false, dirs)
    }

    /// Load a claim view before setting its ranges with `set_claim_scope`.
    /// Ranged material checks belong to the claim that reads the material.
    pub fn load_for_claims(registry: &'static Registry, path: &Path, overrides: &[(String, Given)], t: f64) -> Result<ModelView, InkError> {
        Self::load_for(registry, path, overrides, t, true, &[])
    }

    fn load_for(
        registry: &'static Registry,
        path: &Path,
        overrides: &[(String, Given)],
        t: f64,
        claim_materials: bool,
        dirs: &[std::path::PathBuf],
    ) -> Result<ModelView, InkError> {
        let doc = Doc::read_file(registry, path, overrides, claim_materials, dirs).map_err(|(e, _)| e)?;
        // Stopped at a declared case a text param refuses (Doc::refused):
        // a view of what was read, for its caller to set aside.
        if !doc.refused.is_empty() {
            let overrides = doc.overrides.clone();
            return Ok(ModelView { path: Some(path.to_path_buf()), overrides, ..ModelView::new(doc, Ctx::new(t, Vec::new())) });
        }
        // The values each declared case took (Given::case), as the load
        // read them.
        let overrides = &doc.overrides.clone();
        let (ctx, unsolved) = doc.ctx_at(overrides, t)?;
        doc.check_conditions(&ctx, None, &unsolved.iter().map(|u| u.0).collect::<Vec<_>>(), false, None)?;
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
        let v = self.reload_at(changes, self.ctx.t)?;
        if let Some((_, e)) = v.unsolved.first() {
            return Err(e.clone());
        }
        Ok(v)
    }

    /// Reload at an explicit moment, independent of the cursor left by a previous claim.
    pub fn reload_at(&self, changes: &[(String, Given)], t: f64) -> Result<ModelView, InkError> {
        let path = self.path.as_ref().ok_or_else(|| InkError::msg("only a model loaded from a file can be reloaded with other params"))?;
        let mut overrides = self.overrides.clone();
        for (k, g) in changes {
            let mut next = g.clone();
            // Optimizer steps are in canonical units. Preserve the quantity
            // for these internally generated numbers, never for authored cases.
            if next.unit.is_none() && next.text.is_none() && next.case.is_none() {
                next.unit = self.doc.param_unit(k);
            }
            match overrides.iter_mut().find(|(n, _)| n == k) {
                Some((_, was)) => *was = next,
                None => overrides.push((k.clone(), next)),
            }
        }
        ModelView::load_for(self.doc.registry, path, &overrides, t, self.doc.claim_materials, self.doc.given_import_dirs())
    }

    /// Provide, in this view's context, the solves on demand that the
    /// fields at `roots` rest on (a published mass's density a command
    /// weighs: Doc::want), and from now on at every moment. Err: one
    /// that cannot be provided, at its form.
    pub fn want(&mut self, roots: &[u32]) -> Result<(), InkError> {
        if self.doc.solves.iter().all(|s| !s.on_demand) {
            return self.doc.check_magnitudes(&self.ctx, self.motion, &self.unsolved_solves(), false, Some(roots));
        }
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
            self.unsolved.retain(|(j, _)| *j != i);
            let ranged = reads.iter().copied().find(|&k| self.ctx.p_spans.get(k as u32).is_some());
            // An estimate's guard, over ranges (its own or others'): shown as set_moment shows it
            // (Doc::show_guard), so a moment reached either way says the same. Only the guard is
            // constant over the ranges: its value stays in the original graph, its marker
            // uncertified. One on t is shown here at this moment, at each other by set_moment, and
            // over a span of t by what ranges it.
            let shown = (!self.ctx.p_spans.is_empty()).then(|| self.doc.show_guard(i, &self.ctx)).flatten();
            let result = if let Some(shown) = shown {
                let marker = if shown.is_ok() { 0.0 } else { f64::NAN };
                s.owns.iter().for_each(|&k| self.ctx.params[k] = marker);
                shown.map_err(|m| self.doc.solve_error(i, m))
            } else if let Some(k) = ranged.filter(|_| !s.owns.iter().all(|&o| self.ctx.p_spans.is_derived(o as u32))) {
                let what: Vec<String> = s.owns.iter().map(|&o| format!("'{}'", self.doc.params[o].name)).collect();
                Err(InkError {
                    hint: Some(self.doc.at_one_value(&[k], &self.ctx.p_spans)),
                    open: Some(crate::Open::Outside),
                    ..self.doc.solve_error(
                        i,
                        format!(
                            "{} is found at one value of each param, but reads ranged '{}'. That value says nothing of the rest of its range",
                            what.join(", "),
                            self.doc.params[k].name
                        ),
                    )
                })
            } else {
                self.doc.provide_one(i, &mut self.ctx.params, self.ctx.t).map_err(|m| self.doc.solve_error(i, m))
            };
            if let Err(e) = result {
                // Stepping a moving claim also consults this refusal.
                self.unsolved.push((i, e.clone()));
                return Err(e);
            }
        }
        self.doc.check_magnitudes(&self.ctx, self.motion, &self.unsolved_solves(), false, Some(roots))
    }

    /// Every solve on demand (a published mass) that would fail in this
    /// view's context, as `ink check` asks (ParamProvider::refusal: a mass
    /// is not measured to say so), each by the params it sets, with the
    /// error at its form. What they read is provided first.
    pub fn mass_refusals(&mut self) -> Vec<(Vec<usize>, InkError)> {
        self.mass_checks().into_iter().filter(|(_, impossible, _)| *impossible).map(|(k, _, e)| (k, e)).collect()
    }

    /// Infeasible and unresolved published masses, with whether impossibility was shown.
    pub fn mass_checks(&mut self) -> Vec<(Vec<usize>, bool, InkError)> {
        let demand: Vec<usize> = (0..self.doc.solves.len()).filter(|&i| self.doc.solves[i].on_demand).collect();
        let reads: Vec<usize> = demand.iter().flat_map(|&i| self.doc.solves[i].reads.clone()).collect();
        self.doc.want_params(&reads);
        let mut failed: Vec<(Vec<usize>, bool, InkError)> = Vec::new();
        let mut unresolved = std::collections::HashSet::new();
        for i in demand {
            let solve = &self.doc.solves[i];
            let dependencies = self.doc.solves_under(&solve.reads);
            // Provide the measured inputs in dependency order, including an
            // inner mass whose density only another provider reads.
            for j in dependencies.iter().copied().filter(|&j| self.doc.solves[j].on_demand) {
                let _ = self.doc.provide_one(j, &mut self.ctx.params, self.ctx.t);
            }
            let inner = dependencies.iter().find(|j| unresolved.contains(*j));
            let status = if let Some(&j) = inner.filter(|_| self.doc.ext[solve.ext].package == "mass") {
                let name = self.doc.ext[self.doc.solves[j].ext].form.name();
                Ok(Some(format!("published mass feasibility depends on an unresolved inner mass ({name})")))
            } else {
                let p = self.doc.ext[solve.ext].form.provider(&self.doc).expect("a solve provides params");
                let inputs = self.doc.measurement_inputs(i, &self.ctx.params, self.ctx.t);
                p.feasibility(&self.ctx.params, self.ctx.t, inputs.as_deref())
            };
            let (impossible, why) = match status {
                Err(m) => (true, m),
                Ok(Some(m)) => (false, crate::Refused::from(m)),
                Ok(None) => continue,
            };
            unresolved.insert(i);
            failed.push((solve.owns.clone(), impossible, self.doc.solve_error(i, why)));
        }
        failed
    }

    /// The solves that do not solve at the view's moment (by index in
    /// Doc::solves).
    pub(crate) fn unsolved_solves(&self) -> Vec<usize> {
        self.unsolved.iter().map(|u| u.0).collect()
    }

    /// The params with no value at the view's moment: set by a solve that
    /// does not solve there, or by one resting on such a solve.
    pub fn unsolved_params(&self) -> Vec<usize> {
        let failed = self.unsolved_solves();
        let d = &self.doc;
        (0..d.solves.len()).filter(|&i| d.solves_under(&d.solves[i].owns).iter().any(|j| failed.contains(j))).flat_map(|i| d.solves[i].owns.clone()).collect()
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
        self.doc.covers_t(self.ctx.t, span)?;
        if let Some((a, b)) = span
            && !(a <= self.ctx.t && self.ctx.t <= b)
        {
            self.set_moment(a)?;
        }
        if span.is_some() {
            self.doc.check_conditions(&self.ctx, span, &self.unsolved_solves(), false, None)?;
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
        self.doc.add_param(crate::Param { name: "t".into(), default: self.ctx.t, tol: None });
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
        self.doc.covers_t(t, None)?;
        // A solve on t is provided again here, where it is provided at all
        // (one on demand, while something reads it).
        let again: Vec<usize> = (0..self.doc.solves.len()).filter(|&i| self.doc.solves[i].on_t && self.doc.provided_always(i)).collect();
        // An estimate's guard (estimate-when) is a formula: shown again at t
        // over the ranges by interval bounds, its marker 0 where it holds;
        // where it may not, refused there (its marker no number).
        let guard = |doc: &Doc, i: usize| crate::estimate::Estimate::of(&doc.ext[doc.solves[i].ext]).is_some();
        if !again.is_empty() && !self.ctx.p_spans.is_empty() && again.iter().all(|&i| guard(&self.doc, i)) {
            let mut ctx = Ctx { t, ..self.ctx.clone() };
            let mut unsolved: Vec<(usize, InkError)> = self.unsolved.iter().filter(|(j, _)| !again.contains(j)).cloned().collect();
            for &i in &again {
                let shown = self.doc.show_guard(i, &ctx).expect("a guard");
                let marker = if shown.is_ok() { 0.0 } else { f64::NAN };
                self.doc.solves[i].owns.iter().for_each(|&k| ctx.params[k] = marker);
                if let Err(m) = shown {
                    unsolved.push((i, self.doc.solve_error(i, m)));
                }
            }
            if !self.doc.claim_materials {
                self.doc.check_magnitudes(&ctx, None, &unsolved.iter().map(|u| u.0).collect::<Vec<_>>(), false, None)?;
            }
            (self.ctx, self.unsolved) = (ctx, unsolved);
            self.outcomes.clear();
        } else if !again.is_empty() {
            if !self.ctx.p_spans.is_empty() {
                let ks: Vec<usize> = self.ctx.p_spans.iter().filter(|p| !p.derived).map(|p| p.k as usize).collect();
                return Err(InkError::msg(format!("a solve that reads t is not solved again over ranged params at each step: {}", self.doc.at_one_value(&ks, &self.ctx.p_spans))));
            }
            let (ctx, unsolved) = self.doc.ctx_at(&self.overrides, t)?;
            if !self.doc.claim_materials {
                self.doc.check_magnitudes(&ctx, None, &unsolved.iter().map(|u| u.0).collect::<Vec<_>>(), false, None)?;
            }
            (self.ctx, self.unsolved) = (ctx, unsolved);
            self.outcomes.clear();
        } else {
            let was = std::mem::replace(&mut self.ctx.t, t);
            if !self.doc.claim_materials
                && let Err(e) = self.doc.check_magnitudes(&self.ctx, None, &self.unsolved_solves(), false, None)
            {
                self.ctx.t = was;
                return Err(e);
            }
        }
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
        self.outcomes.retain(|k, _| doc.ext_form(Some(k.0.0), &k.0.1).is_none_or(|e| doc.form_motion_inputs(e).iter().all(still)));
        self.timeless = memo.into_inner();
        Ok(())
    }

    /// The model's accuracy (twice Doc::finest, at most 0.05 mm): what a
    /// shape claim is checked to unless it states its own.
    pub fn accuracy(&self) -> f64 {
        (2.0 * self.finest()).min(0.05)
    }

    /// A stated fact's number at the view's params, and the low and high
    /// it has over the ranged params it reads (a tolerance, a `(ranged lo
    /// hi)`, a range given): none where it reads none. Err, saying why,
    /// where it is no number. ink's values are total (an operation with no
    /// result gives 0), so that is asked here of every operation in it, by
    /// the one table of where each has none (graph::undefined_un_over,
    /// undefined_bin_over; and a sqrt below 0, which is 0), over its
    /// operands' ranges; of the numbers the model is built from
    /// (validity: an overflow, a number outside ink's range); and of its
    /// own range. A value ink works out at a run (a solve, an estimate) is
    /// not read. Worked out when asked, not at load.
    pub fn stated(&self, fact: &crate::Stated) -> Result<(f64, Option<(f64, f64)>), String> {
        use ink_core::graph::{undefined_bin_over, undefined_un_over};
        let g = &self.doc.g;
        let range = |n| {
            let over = ink_core::eval::eval_interval_scope(&Tape::new(g, n), &self.ctx, [ink_core::interval::I::point(0.0); 3], None);
            [over.lo, over.hi]
        };
        if let Some(why) = crate::validity::numbers(&self.doc, &self.ctx, &[fact.value]) {
            return Err(why);
        }
        // A solve's params (its provider's: Doc::keep_form), found from the solves themselves: compiling
        // the providers clones every form's tapes, and a rating's formula is asked this at every moment.
        if let Some(k) = self.doc.params_read(&[fact.value]).into_iter().find(|&k| self.doc.solves.iter().any(|s| s.owns.contains(&k))) {
            return Err(format!("it reads '{}', which ink works out at a run (a solve, an estimate): a fact does not read one yet", self.doc.params[k].name));
        }
        if g.may_lack_a_result(fact.value) {
            for n in g.reachable(&[fact.value]) {
                let lacks = match g.nodes[n as usize] {
                    Node::Un(ink_core::graph::UnOp::Sqrt, a) if range(a)[0] < 0.0 => Some("the square root of a negative"),
                    Node::Un(op, a) => undefined_un_over(op, range(a)),
                    Node::Bin(op, a, b) => undefined_bin_over(op, range(a), range(b)),
                    _ => None,
                };
                if let Some(what) = lacks {
                    return Err(format!("{what}, or what may be one, is in it"));
                }
            }
        }
        let ([low, high], value) = (range(fact.value), self.doc.scalar(fact.value, &self.ctx));
        if !(value.is_finite() && low.is_finite() && high.is_finite()) {
            return Err(format!("it is somewhere from {low} to {high}"));
        }
        Ok((value, (!self.ranged_reads(&[fact.value]).is_empty()).then_some((low, high))))
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
            if s.owns.iter().all(|&u| spans.is_derived(u as u32))
                && crate::estimate::Estimate::of(&self.doc.ext[s.ext]).is_none()
                && self.doc.inexact(&self.doc.form_inputs(&self.doc.ext[s.ext]), &own).solved.is_empty()
            {
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
            if let Some(((a, b), site)) = self.doc.covered(k)
                && !(a <= lo && hi <= b)
            {
                return Err(InkError {
                    hint: Some(format!("range it within {a}..{b}, or widen its band where it is declared")),
                    ..InkError::msg(format!("'{}' is ranged over {lo}..{hi}, beyond {a}..{b}, the band a box found at load at {site} holds over", self.doc.params[k].name))
                });
            }
        }
        let mut ctx = Ctx { p_spans: ink_core::eval::ParamSpans::new(spans.iter().map(|&(k, s)| (k as u32, s))), ..self.ctx.clone() };
        let doc = &self.doc;
        // Over the solves provided here: not one on demand nothing reads
        // (a published mass a claim over the range does not weigh). Rates
        // with no proven range are not derived: what reads them is refused
        // claim by claim (set_param_spans).
        for (i, s) in doc.solves.iter().enumerate().filter(|&(i, _)| doc.provided_always(i)) {
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
            // An estimate's guard: shown at this moment over the ranges (Doc::show_guard), as want
            // and set_moment show it. One on t over a motion is shown at each moment a claim reads
            // and over each span of t, where it is read: a moment it fails at refuses what reads it
            // there, not the scope.
            if crate::estimate::Estimate::of(e).is_some() {
                if !(s.on_t && self.motion.is_some()) {
                    doc.show_guard(i, &ctx).expect("a guard").map_err(|m| InkError { hint, ..InkError::msg(format!("{what}: {m}")) })?;
                }
                continue;
            }
            // A provider that encloses what it sets over the range (a published mass's tolerance):
            // that enclosure, its param derived so.
            if let [own] = p.owns()
                && let Some(over) = p.over_ranged(&ctx.params, ctx.t, &ranged.iter().filter_map(|&k| ctx.p_spans.get(k as u32).map(|s| (k, s))).collect::<Vec<_>>())
            {
                let (value, band) = over.map_err(|r| InkError { hint: hint.clone(), ..InkError::msg(format!("{what}: {}", r.why)) })?;
                ctx.params[*own] = value;
                let derive = OverRanged { provider: p.clone(), owns: [*own as u32], reads: ranged.iter().map(|&k| k as u32).collect() };
                ctx.p_spans = ctx.p_spans.derived(std::sync::Arc::new(derive), &[(*own as u32, band)]);
                continue;
            }
            // No form derives what it sets over a range: it sets its params at single values.
            return Err(InkError { hint, ..InkError::msg(format!("{what}, and what it solves has no proven range: it sets its params at single values only")) });
        }
        if !ctx.p_spans.is_empty() {
            doc.check_conditions(&ctx, self.motion, &self.unsolved_solves(), false, None)?;
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
        self.reshow_guards();
        self.reaches.clear();
        self.boxes.clear();
        self.tapes.clear();
        self.outcomes.clear();
        self.found.clear();
        Ok(())
    }

    /// Every estimate's guard already provided, shown again at this moment
    /// in the scope just set, as `want` shows it (Doc::show_guard over
    /// ranges, else at the params' values), so `refusal` says what this
    /// scope says without another `want` or `set_moment`.
    fn reshow_guards(&mut self) {
        for i in 0..self.doc.solves.len() {
            if !self.doc.provided_always(i) || crate::estimate::Estimate::of(&self.doc.ext[self.doc.solves[i].ext]).is_none() {
                continue;
            }
            let shown = match self.ctx.p_spans.is_empty() {
                false => self.doc.show_guard(i, &self.ctx).expect("a guard"),
                true => self.doc.provide_one(i, &mut self.ctx.params, self.ctx.t).map_err(|r| r.why),
            };
            let marker = if shown.is_ok() { 0.0 } else { f64::NAN };
            self.doc.solves[i].owns.iter().for_each(|&k| self.ctx.params[k] = marker);
            self.unsolved.retain(|(j, _)| *j != i);
            if let Err(m) = shown {
                let e = self.doc.solve_error(i, m);
                self.unsolved.push((i, e));
            }
        }
    }

    /// Answer claims over `scope` (ClaimScope): its motion, its
    /// params' spans.
    pub fn set_claim_scope(&mut self, scope: &ClaimScope) -> Result<(), InkError> {
        self.by_file = scope.by_file.then_some(self.ctx.t);
        // Made now, where an import's span may be checked over, so the
        // answers' param names hold it (set_motion).
        if self.by_file.is_some() && self.doc.spans_hull().is_some() {
            self.t_param();
        }
        self.set_motion(scope.motion)?;
        self.set_param_spans(&scope.spans)?;
        if self.doc.claim_materials && self.ctx.p_spans.is_empty() && scope.motion.is_none() && !(scope.by_file && self.doc.spans_hull().is_some()) {
            self.doc.check_magnitudes(&self.ctx, None, &self.unsolved_solves(), false, None)?;
        }
        Ok(())
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
    pub fn motion_of(&mut self, file: crate::Scope) -> Result<(), InkError> {
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
        let (key, f, inputs) = ((e.package, form.to_string()), e.form.clone(), self.doc.form_motion_inputs(e));
        self.want(&inputs)?;
        let key = (key, self.scope(&inputs));
        if let Some((o, out)) = self.outcomes.get(&key)
            && *o == opts
        {
            return Ok(out.clone());
        }
        // A form that reads what the model works out over the motion is run where that is largest.
        let ctx = match self.at_worst(package, form)? {
            Some((ctx, _)) => ctx,
            None => self.ctx.clone(),
        };
        let out: Rc<dyn Outcome> = f.run(&mut self.doc, &ctx, &opts)?.into();
        self.outcomes.insert(key, (opts, out.clone()));
        Ok(out)
    }

    /// The moments of the motion a form's `worst_reads` are read at: its span's ends and the middles of
    /// this many even steps.
    const WORST_STEPS: usize = 64;

    /// Where a form reads values the model works out over the motion (`Form::worst_reads`): the context
    /// to run it in, and what to say of it. The context is the design's at one moment of the span it is
    /// asked over (`motion`, else the file's (t-span)), the solves provided there: the moment where those
    /// values are largest. With one group (one load) that is its largest size; with several, where the
    /// sizes, each over its own largest, are largest together (the root of the sum of their squares), so
    /// that all are of one moment, a load case the design meets. What that can miss, said: a larger value
    /// between the moments read, and a smaller load in a direction the part bears worse. With no span,
    /// the view's own moment. None: the form reads no such value. Err: a value the model cannot give,
    /// with what is missing (the solve's own error).
    pub fn at_worst(&mut self, package: Option<&str>, form: &str) -> Result<Option<(Ctx, Vec<String>)>, InkError> {
        let e = self.doc.ext_form(package, form).ok_or_else(|| InkError::msg(format!("no package form named '{form}'")))?;
        let key = (e.package, form.to_string());
        if let Some(found) = self.worst.get(&key) {
            return Ok(found.clone());
        }
        let groups = e.form.worst_reads();
        if groups.is_empty() {
            self.worst.insert(key, None);
            return Ok(None);
        }
        let nodes: Vec<u32> = groups.iter().flat_map(|(_, nodes)| nodes.iter().copied()).collect();
        self.want(&nodes)?;
        if let Some((k, (low, high))) = self.ranged_reads(&nodes).first().copied() {
            let name = &self.doc.params[k as usize].name;
            return Err(InkError::msg(format!(
                "'{form}' reads what the model works out, which rests on '{name}', ranged {low} to {high}: it is found at one value of each param (--param {name}=V)"
            )));
        }
        let reads = self.doc.params_read(&nodes);
        let span = self.motion.or(self.doc.t_span);
        // The moments are those the model's own answers are read at (its providers' steps: they follow how
        // fast its coordinates and loads change, halved to a budget), so a load that swings between even
        // steps is read where it swings; with none to say, even steps.
        let mut gaps: Vec<(f64, f64)> = Vec::new();
        let mut how = format!("the middles of {} even steps", Self::WORST_STEPS);
        let times: Vec<f64> = match span {
            Some((a, b)) => {
                let chosen: Vec<crate::Steps> = (self.doc.solves_read(&nodes).into_iter())
                    .filter_map(|i| self.doc.ext[self.doc.solves[i].ext].form.provider(&self.doc)?.steps(&self.ctx.params, (a, b)))
                    .filter_map(Result::ok)
                    .collect();
                match crate::motion::merged(chosen) {
                    Some(steps) if steps.ends.len() > 2 => {
                        gaps = steps.unchecked.iter().chain(&steps.loose).copied().collect();
                        how = format!("the middles of the {} steps the model's own answers are read at, which follow how fast its motion and loads change", steps.ends.len() - 1);
                        std::iter::once(a).chain(steps.ends.windows(2).map(|w| ink_core::eval::span_mid((w[0], w[1])))).chain([b]).collect()
                    }
                    _ => {
                        let h = (b - a) / Self::WORST_STEPS as f64;
                        std::iter::once(a).chain((0..Self::WORST_STEPS).map(|i| a + h * (i as f64 + 0.5))).chain([b]).collect()
                    }
                }
            }
            None => vec![self.ctx.t],
        };
        gaps.sort_by(|x, y| x.0.total_cmp(&y.0));
        self.worst_gaps.insert(key.clone(), gaps);
        // Each group's size at each moment.
        let mut read: Vec<(Ctx, Vec<f64>)> = Vec::new();
        let mut all: Vec<(Ctx, Vec<f64>)> = Vec::new();
        for &t in &times {
            let (ctx, failed) = self.doc.ctx_at(&self.overrides, t).map_err(InkError::msg)?;
            if let Some((_, why)) = failed.into_iter().find(|(i, _)| self.doc.solves[*i].owns.iter().any(|k| reads.contains(k))) {
                return Err(InkError { message: format!("'{form}' reads what the model does not give at t = {t}: {}", why.message), ..why });
            }
            all.push((Ctx { t_span: None, ..ctx.clone() }, nodes.iter().map(|&n| self.doc.scalar(n, &ctx)).collect()));
            let sizes: Vec<f64> = groups.iter().map(|(_, nodes)| nodes.iter().map(|&n| self.doc.scalar(n, &ctx).powi(2)).sum::<f64>().sqrt()).collect();
            if let Some(g) = sizes.iter().position(|x| !x.is_finite()) {
                return Err(InkError::msg(format!("'{form}': {} is no number at t = {t}: the model does not give what it reads there", groups[g].0)));
            }
            read.push((ctx, sizes));
        }
        self.worst_all.insert(key.clone(), all);
        let most: Vec<f64> = (0..groups.len()).map(|g| read.iter().map(|(_, s)| s[g]).fold(0.0, f64::max)).collect();
        let score = |sizes: &[f64]| sizes.iter().zip(&most).filter(|(_, m)| **m > 0.0).map(|(s, m)| (s / m).powi(2)).sum::<f64>();
        let at = (0..read.len()).max_by(|&a, &b| score(&read[a].1).total_cmp(&score(&read[b].1))).expect("a moment");
        let (ctx, sizes) = read.swap_remove(at);
        let sig = |x: f64| crate::claims::sig(x);
        let each: Vec<String> = groups
            .iter()
            .zip(&sizes)
            .zip(&most)
            .map(
                |(((name, _), size), most)| {
                    if size == most { format!("{name} {}", sig(*size)) } else { format!("{name} {} (its largest at any moment: {})", sig(*size), sig(*most)) }
                },
            )
            .collect();
        let note = match span {
            Some((a, b)) => format!(
                "run with what the model works out at t = {}, the moment of t {a} to {b} where it is largest{}: {}. Read at the span's ends and {how}: a \
                 larger value between them, and a smaller one in a direction the part bears worse, are not looked at",
                sig(ctx.t),
                if groups.len() > 1 { " (each over its own largest, together)" } else { "" },
                each.join("; "),
            ),
            None => format!("run with what the model works out at t = {}, the one moment asked (the file states no (t-span A B)): {}", sig(ctx.t), each.join("; ")),
        };
        let found = Some((Ctx { t_span: None, ..ctx }, vec![note]));
        self.worst.insert(key, found.clone());
        Ok(found)
    }

    /// A package's metric of one of its forms or of a shape: the one path
    /// claims, requirements and objectives on package metrics are
    /// evaluated by. A failed run's metric has no value (`Outcome::failure`).
    pub fn measure(&mut self, m: &FormMetric, opts: RunOpts) -> Result<Measured, InkError> {
        if let Some(f) = m.shape {
            let p = self.doc.registry.packages.iter().find(|p| p.name() == m.package).expect("a metric's package");
            let accuracy = opts.accuracy.unwrap_or_else(|| self.accuracy());
            return p.shape_metric(self, m.metric, f, opts.res, accuracy);
        }
        if let Some(why) = self.doc.ext_form(Some(m.package), &m.form).and_then(|e| e.form.unavailable_in(&self.doc, &self.ctx, m.metric)) {
            return Ok(Measured::none(why));
        }
        // A form loaded from the model over a motion: its worst over every moment read, where its response
        // is linear in what it reads; else at the one moment that is largest, which is then no whole answer.
        if let Some((_, notes)) = self.at_worst(Some(m.package), &m.form)? {
            let form = self.doc.ext_form(Some(m.package), &m.form).expect("found").form.clone();
            // One moment too, where the model gives its loads within a band: carried as over a motion.
            let nodes: Vec<u32> = form.worst_reads().iter().flat_map(|(_, n)| n.iter().copied()).collect();
            let several = self.worst_all.get(&(m.package, m.form.clone())).is_some_and(|all| all.len() > 1 || all.iter().any(|(ctx, _)| !self.bands_at(&nodes, ctx).is_empty()));
            let between = self.worst_gaps.get(&(m.package, m.form.clone())).and_then(|gaps| crate::motion::not_followed(crate::motion::MODEL_GAP, gaps));
            let measured = self.measure_at_worst(m, opts, &form, several, notes)?;
            return Ok(match between {
                Some(why) if measured.value.is_ok() => Measured { value: Err(why), ..measured },
                _ => measured,
            });
        }
        let out = self.outcome(Some(m.package), &m.form, opts)?;
        self.measure_tail(m, out)
    }

    fn measure_at_worst(&mut self, m: &FormMetric, opts: RunOpts, form: &Arc<dyn crate::ext::Form>, several: bool, notes: Vec<String>) -> Result<Measured, InkError> {
        match form.linear_in_worst(m.metric) {
            Some(smaller) if several => {
                let all = self.measure_over(m, opts, form.besides_worst())?;
                // One moment's, where a claim asks moment by moment; else the worst value's.
                let worse = |a: &Measured, b: &Measured| match (&a.value, &b.value) {
                    (Ok(a), Ok(b)) => {
                        if smaller {
                            b.total_cmp(a)
                        } else {
                            a.total_cmp(b)
                        }
                    }
                    (a, b) => a.is_err().cmp(&b.is_err()),
                };
                let k = self.pinned.unwrap_or_else(|| (0..all.len()).max_by(|&a, &b| worse(&all[a], &all[b])).expect("a moment"));
                Ok(all[k].clone())
            }
            _ => {
                let out = self.outcome(Some(m.package), &m.form, opts)?;
                let mut measured = out.metric(m.metric).ok_or_else(|| InkError::msg(format!("'{}' gave no '{}'", m.form, m.metric)))?;
                if let Some(why) = out.failure() {
                    measured.value = Err(why);
                }
                measured.warnings.extend(notes);
                // One moment of several, and no bound over the rest: a value, and no answer for the motion.
                if several && measured.value.is_ok() {
                    measured.unseen = Some(format!(
                        "'{}' is not linear in its loads, so its worst over the motion is not bounded from a few runs: this is its value at one moment, and another \
                             moment's load, smaller but in a direction the part bears worse, may give a worse one",
                        m.metric
                    ));
                    measured.unseen_finer = false;
                }
                Ok(measured)
            }
        }
    }

    fn measure_tail(&mut self, m: &FormMetric, out: Rc<dyn Outcome>) -> Result<Measured, InkError> {
        let mut measured = out.metric(m.metric).ok_or_else(|| InkError::msg(format!("'{}' gave no '{}'", m.form, m.metric)))?;
        if let Some(why) = out.failure() {
            measured.value = Err(why);
        }
        Ok(measured)
    }

    /// How many moments a claim on `m` is decided at, one by one (`pin`): the moments of the motion whose
    /// loads differ, where `m`'s form reads the model's values over a motion and its metric is linear in
    /// them (`measure_over`). 0: it is decided once.
    pub fn moments_over(&mut self, m: &FormMetric) -> Result<usize, InkError> {
        if m.shape.is_some() || self.at_worst(Some(m.package), &m.form)?.is_none() {
            return Ok(0);
        }
        let form = self.doc.ext_form(Some(m.package), &m.form).expect("found").form.clone();
        let all = self.worst_all.get(&(m.package, m.form.clone())).map_or(0, |all| all.len());
        Ok(if all > 1 && form.linear_in_worst(m.metric).is_some() { self.distinct(m).len() } else { 0 })
    }

    /// `measure` answers for moment `k` of `moments_over` alone (None: for the one where its value is worst).
    pub fn pin(&mut self, k: Option<usize>) {
        self.pinned = k;
    }

    /// Whether moment `a` of `moments_over` reads worse on `m` than moment `b`, by their values on its own grid:
    /// which of two moments a claim answers alike at is the one it shows.
    pub fn worse_moment(&mut self, m: &FormMetric, opts: RunOpts, a: usize, b: usize) -> bool {
        let form = self.doc.ext_form(Some(m.package), &m.form).expect("found").form.clone();
        let (Some(smaller), Ok(all)) = (form.linear_in_worst(m.metric), self.measure_over(m, opts, form.besides_worst())) else { return false };
        match (&all[a].value, &all[b].value) {
            (Ok(a), Ok(b)) => {
                if smaller {
                    a < b
                } else {
                    a > b
                }
            }
            _ => false,
        }
    }

    /// The params `nodes` read that the model gives within a band at the moment `ctx`, each with its
    /// half-width there: what its provider's own measure leaves (ParamProvider::error), and how far
    /// the real value may be from the reading over the step (ParamProvider::band: a spinning body
    /// read round), about the value set.
    pub(crate) fn bands_at(&self, nodes: &[u32], ctx: &Ctx) -> Vec<(usize, f64)> {
        let read = self.doc.params_read(nodes);
        let mut out: Vec<(usize, f64)> = Vec::new();
        for i in self.doc.solves_read(nodes) {
            let s = &self.doc.solves[i];
            let Some(p) = self.doc.ext[s.ext].form.provider(&self.doc) else { continue };
            let mut widths = vec![0.0; s.owns.len()];
            if let Some(w) = p.error(&ctx.params, ctx.t) {
                widths.iter_mut().zip(w).for_each(|(x, w)| *x += w);
            }
            if let Some((w, centre)) = p.band(&ctx.params, ctx.t) {
                for ((x, w), (c, &k)) in widths.iter_mut().zip(w).zip(centre.iter().zip(&s.owns)) {
                    *x += w + (ctx.params[k] - c).abs();
                }
            }
            out.extend(s.owns.iter().copied().zip(widths).filter(|&(k, w)| read.contains(&k) && w != 0.0));
        }
        out
    }

    /// Each param banded at some moment (`bands`), with how fast each of `nodes` (a form's loads)
    /// changes with it: the load is the param times that plus what does not read it, exactly
    /// (period::form). None: some load is no such form of some banded param.
    fn banded(&self, nodes: &[u32], bands: &[Vec<(usize, f64)>], ctx: &Ctx) -> Option<Vec<(usize, Vec<f64>)>> {
        let mut ks: Vec<usize> = bands.iter().flatten().map(|b| b.0).collect();
        ks.sort_unstable();
        ks.dedup();
        ks.into_iter()
            .map(|k| {
                let at = ctx.clone().with_band(k as u32, (ctx.params[k], ctx.params[k]));
                let rate = nodes
                    .iter()
                    .map(|&n| match ink_core::period::form(&Tape::new(&self.doc.g, n), Var::Param(k as u32), &at)? {
                        None => Some(0.0),
                        Some(c) => {
                            let c = c.value();
                            Some(if c.lo.abs() > c.hi.abs() { c.lo } else { c.hi })
                        }
                    })
                    .collect::<Option<Vec<f64>>>()?;
                Some((k, rate))
            })
            .collect()
    }

    /// The moments read for `m`'s form whose loads differ, as places in `worst_all`: a motion that
    /// dwells, or goes and returns, reads one load at many moments. Where the model's loads are all
    /// that drive the form and every moment's are one moment's times a number, the response is that
    /// moment's times the number, and every size a claim reads of it is largest where the number is:
    /// only that moment, and the one where it is most negative (its loads reversed are a run of their
    /// own), answer for the rest.
    fn distinct(&self, m: &FormMetric) -> Vec<usize> {
        let all = self.worst_all.get(&(m.package, m.form.clone())).map_or(&[][..], |all| all.as_slice());
        let mut moments: Vec<usize> = Vec::new();
        for i in 0..all.len() {
            if moments.iter().all(|&k| all[k].1 != all[i].1) {
                moments.push(i);
            }
        }
        let besides = self.doc.ext_form(Some(m.package), &m.form).is_none_or(|e| e.form.besides_worst());
        // A moment of less load may be known less well: with a band, each is decided.
        let nodes: Vec<u32> = self.doc.ext_form(Some(m.package), &m.form).map_or(Vec::new(), |e| e.form.worst_reads().iter().flat_map(|(_, n)| n.iter().copied()).collect());
        if besides || moments.len() < 3 || moments.iter().any(|&i| !self.bands_at(&nodes, &all[i].0).is_empty()) {
            return moments;
        }
        let loads: Vec<&[f64]> = moments.iter().map(|&i| all[i].1.as_slice()).collect();
        let (picked, weights) = spanning(&loads, false);
        if picked.len() != 1 {
            return moments;
        }
        let by = |pick: fn(f64, f64) -> bool| (1..weights.len()).fold(0, |best, k| if pick(weights[k][0], weights[best][0]) { k } else { best });
        let (most, least) = (by(|a, b| a > b), by(|a, b| a < b));
        let mut kept = vec![moments[most]];
        if weights[least][0] < 0.0 && least != most {
            kept.push(moments[least]);
        }
        kept.sort_unstable();
        kept
    }

    /// `measure`, of a form that reads the model's values over a motion (`Form::worst_reads`) and whose
    /// metric is a size of a response linear in them (`Form::linear_in_worst`: the largest deflection or
    /// stress of a static study in its loads): its whole answer at every moment whose loads differ
    /// (`moments_over`), on one grid. Every moment's loads are a weighted sum of a few moments' (`spanning`:
    /// as many as the loads have directions over the motion, six at most for one force and torque), so
    /// its response is that sum of those moments' runs, with nothing left to bound, and all a run's
    /// answer holds is read off each sum (`Form::of_sums`: what it leaves out, what its grid cannot
    /// see, its read's error; a buckling factor bounded). Only the model's loads may differ between
    /// the moments: where anything else the form reads does (`Form::steady_besides`), no moment is a
    /// sum of others, and none has a value. `besides`: something else drives it too, which each sum
    /// then holds once.
    fn measure_over(&mut self, m: &FormMetric, opts: RunOpts, besides: bool) -> Result<Rc<Vec<Measured>>, InkError> {
        let key = (m.package, m.form.clone(), m.metric, opts.res.map(f64::to_bits), opts.split_as.map(f64::to_bits));
        if let Some(found) = self.over.get(&key) {
            return Ok(found.clone());
        }
        let e = self.doc.ext_form(Some(m.package), &m.form).expect("a form");
        let (form, inputs) = (e.form.clone(), self.doc.form_motion_inputs(e));
        self.want(&inputs)?;
        let all = self.worst_all.get(&(m.package, m.form.clone())).cloned().unwrap_or_default();
        let metric = m.metric;
        let moments = self.distinct(m);
        let every = |one: Measured| Rc::new(vec![one; moments.len()]);
        let measured_of = |out: &dyn Outcome| -> Result<Measured, InkError> {
            let mut measured = out.metric(metric).ok_or_else(|| InkError::msg(format!("'{}' gave no '{metric}'", form.name())))?;
            if let Some(why) = out.failure() {
                measured.value = Err(why);
            }
            Ok(measured)
        };
        let ctxs: Vec<Ctx> = all.iter().map(|(ctx, _)| ctx.clone()).collect();
        let found = if let Err(why) = form.steady_besides(&mut self.doc, &ctxs, &opts) {
            every(Measured::none(format!(
                "over the motion it reads more than the model's loads change: {why}. Its answer at one moment is no sum of other moments' runs then, and its worst over the \
                 motion is not found; ask it over a span where only the model's loads change, or hold the rest still"
            )))
        } else {
            // What the model gives only within a band (ParamProvider::error) moves the loads by its rate
            // times the band, each load a form of each such param exactly (`banded`); the response to
            // that change is a sum of runs too, one more run where no moment's loads lie along it.
            let nodes: Vec<u32> = form.worst_reads().iter().flat_map(|(_, nodes)| nodes.iter().copied()).collect();
            let bands: Vec<Vec<(usize, f64)>> = moments.iter().map(|&i| self.bands_at(&nodes, &all[i].0)).collect();
            let rates = self.banded(&nodes, &bands, &all[moments[0]].0);
            let largest_load = moments.iter().flat_map(|&i| all[i].1.iter().map(|x| x.abs())).fold(0.0, f64::max);
            // A moment as it would be with one such param moved far enough to load the part as the
            // largest load does: its loads are the first moment's and that much of the param's rate.
            let moved: Vec<(usize, Ctx, Vec<f64>, f64)> = rates
                .iter()
                .flatten()
                .filter_map(|(k, rate)| {
                    let most = rate.iter().map(|x| x.abs()).fold(0.0, f64::max);
                    let by = if largest_load > 0.0 { largest_load / most } else { 1.0 / most };
                    (most > 0.0 && by.is_finite()).then(|| {
                        let mut ctx = all[moments[0]].0.clone();
                        ctx.params[*k] += by;
                        let loads = all[moments[0]].1.iter().zip(rate).map(|(l, r)| l + by * r).collect();
                        (*k, ctx, loads, by)
                    })
                })
                .collect();
            let loads: Vec<&[f64]> = moments.iter().map(|&i| all[i].1.as_slice()).chain(moved.iter().map(|m| m.2.as_slice())).collect();
            let (picked, weights) = spanning(&loads, besides);
            let at: Vec<Ctx> = picked.iter().map(|&k| if k < moments.len() { all[moments[k]].0.clone() } else { moved[k - moments.len()].1.clone() }).collect();
            // A run some moment weights negatively is asked what its sums need of its loads reversed.
            let both_ways = |i: usize| weights.iter().any(|w| w[i] < -1e-12);
            let how: Vec<(Ctx, RunOpts)> = at.iter().enumerate().map(|(i, ctx)| (ctx.clone(), RunOpts { both_ways: both_ways(i), ..opts })).collect();
            let runs = form.run_each(&mut self.doc, &how)?;
            if let Some(failed) = runs.iter().find(|r| r.failure().is_some()) {
                every(measured_of(failed.as_ref())?)
            } else {
                let refs: Vec<&dyn Outcome> = runs.iter().map(|r| r.as_ref()).collect();
                let mut sums =
                    form.of_sums(&refs, &weights[..moments.len()], metric).ok_or_else(|| InkError::msg(format!("'{}' gave no '{metric}' of its runs' sum", form.name())))?;
                // The response to each banded param's change alone: the moved moment's less the first's, a unit
                // of the param. Each moment's answer is then known to its bands' shares of those, summed.
                if bands.iter().any(|b| !b.is_empty()) {
                    let unit: Option<Vec<Measured>> = match &rates {
                        Some(_) => {
                            let rows: Vec<Vec<f64>> =
                                (0..moved.len()).map(|j| weights[moments.len() + j].iter().zip(&weights[0]).map(|(a, b)| (a - b) / moved[j].3).collect()).collect();
                            form.of_sums(&refs, &rows, metric)
                        }
                        None => None,
                    };
                    let smaller = form.linear_in_worst(metric).unwrap_or(false);
                    for (measured, band) in sums.iter_mut().zip(&bands) {
                        if band.is_empty() {
                            continue;
                        }
                        let Ok(value) = measured.value else { continue };
                        // Its share of each param's unit response: for a size, added; for a factor that is
                        // a strength over a size, added to one over it.
                        let off = unit.as_ref().and_then(|unit| {
                            band.iter().try_fold(0.0, |sum, &(k, w)| {
                                let one = match moved.iter().position(|m| m.0 == k) {
                                    Some(j) => *unit[j].value.as_ref().ok()?,
                                    // It moves no load (`moved` leaves it out): nothing of it reaches the answer.
                                    None => return Some(sum),
                                };
                                let share = if smaller { if one.is_finite() && one != 0.0 { w / one.abs() } else { 0.0 } } else { w * one.abs() };
                                share.is_finite().then_some(sum + share)
                            })
                        });
                        match off {
                            Some(off) => {
                                let (low, high) = match smaller {
                                    // A size is never under nothing.
                                    false => ((value - off).max(0.0f64.min(value)), value + off),
                                    true => (1.0 / (1.0 / value + off), if 1.0 / value - off > 0.0 { 1.0 / (1.0 / value - off) } else { f64::INFINITY }),
                                };
                                let names: Vec<String> = band.iter().map(|&(k, w)| format!("{} within {}", self.doc.params[k].name, crate::claims::sig(w))).collect();
                                // What the providers of those values say they rest on, said here: a study's answer
                                // carries no provider's note of its own.
                                let mut rests: Vec<String> = Vec::new();
                                for s in self.doc.solves.iter().filter(|s| band.iter().any(|&(k, _)| s.owns.contains(&k))) {
                                    if let Some(n) = self.doc.ext[s.ext].form.provider(&self.doc).and_then(|p| p.note()).filter(|n| !rests.contains(n)) {
                                        rests.push(n);
                                    }
                                }
                                let rests = match rests.is_empty() {
                                    true => "what the model says its value rests on: a spinning body slowed or centred on its spin, a mass stated".to_string(),
                                    false => format!("what the model says its value rests on: {}", rests.join("; ")),
                                };
                                let (was_low, was_high) = measured.band.as_ref().map_or((value, value), |b| (b.low, b.high));
                                let mut parts = measured.band.as_ref().map_or(Vec::new(), |b| b.parts.clone());
                                parts.push(("what the model gives".into(), (high - value).max(value - low)));
                                measured.band = Some(crate::ext::Band {
                                    low: was_low - (value - low),
                                    high: was_high + (high - value),
                                    parts,
                                    basis: format!(
                                        "its loads are the model's, which it gives within a band ({}): the response is linear in them, so the answer moves by each band's share of the response to that value alone",
                                        names.join(", ")
                                    ),
                                    narrows: measured
                                        .band
                                        .as_ref()
                                        .map_or(Vec::new(), |b| b.narrows.clone())
                                        .into_iter()
                                        .chain([("what the model gives".to_string(), rests)])
                                        .collect(),
                                });
                            }
                            None => {
                                measured.value = Err(format!(
                                    "its loads are the model's, which it gives within a band, and how far that moves the answer is not found (a load that is no plain multiple of what the model gives, or a response with no value): {}",
                                    band.iter().map(|&(k, w)| format!("{} within {}", self.doc.params[k].name, crate::claims::sig(w))).collect::<Vec<_>>().join(", ")
                                ));
                            }
                        }
                    }
                }
                let sig = crate::claims::sig;
                let ran: Vec<String> =
                    picked.iter().zip(&at).map(|(&k, ctx)| if k < moments.len() { sig(ctx.t) } else { format!("{} with one value the model gives moved", sig(ctx.t)) }).collect();
                // One direction of load and nothing besides: the moments decided are where it is largest (`distinct`).
                let largest = !besides && picked.len() == 1 && moments.len() <= 2 && all.len() > 2;
                for (k, measured) in sums.iter_mut().enumerate() {
                    let decided = match largest {
                        true => "every moment's loads are those at one moment times a number, and the response is linear in the loads, so it is largest where they are: the claim is decided there".to_string(),
                        false => format!(
                            "every moment's loads are a weighted sum of those at {} of them (t = {}), and the response is linear in the loads, so each moment's is that sum of those {} runs, and the claim is decided at each",
                            picked.len(),
                            ran.join(", "),
                            picked.len()
                        ),
                    };
                    measured.warnings.push(format!(
                        "its answer at t = {}, the least favourable of the {} moments of the motion read: {decided}. The moments follow how fast the model's motion and loads change; between them nothing is bounded",
                        sig(all[moments[k]].0.t),
                        all.len(),
                    ));
                }
                Rc::new(sums)
            }
        };
        self.over.insert(key, found.clone());
        Ok(found)
    }

    /// The cells shape metric `m` is measured on at `res`, where it is a
    /// grid's estimate (Package::shape_res), nothing measured.
    /// The order metric `m`'s grid error falls at (Package::metric_order).
    pub fn metric_order(&self, m: &FormMetric) -> f64 {
        self.doc.registry.packages.iter().find(|p| p.name() == m.package).map_or(2.0, |p| p.metric_order(m.metric))
    }

    pub fn shape_res(&mut self, m: &FormMetric, res: Option<f64>) -> Option<f64> {
        let f = m.shape?;
        let p = self.doc.registry.packages.iter().find(|p| p.name() == m.package).expect("a metric's package");
        p.shape_res(self, m.metric, f, res)
    }

    /// Where analysis of field `f` looks: a box holding all of its
    /// material, or the window (see Doc::reach).
    pub fn reach(&mut self, f: u32) -> Reach {
        let key = (f, self.scope(&[f]));
        if let Some(r) = self.reaches.get(&key) {
            return *r;
        }
        let ahead = if self.plane { None } else { self.ahead.take_if(|(k, _)| *k == key) };
        let r = match (ahead, self.slab(f)) {
            (Some((_, r)), _) => r,
            (None, Some(s)) => crate::plane_reach(&self.doc.g, self.doc.domain(), s, &self.ctx),
            (None, None) => self.doc.reach(f, &self.ctx),
        };
        self.reaches.insert(key, r);
        r
    }

    /// Whether field `f` is looked at in the plane: it reads no z, and the
    /// view asks in the plane (`in_plane`).
    fn flat(&self, f: u32) -> bool {
        self.plane && !self.doc.g.value_depends_on(f, |v| v == Var::Z)
    }

    /// Where `f` is looked at in the plane (`flat`), f between z = -1 and
    /// 1, whose extent in x and y is f's.
    fn slab(&mut self, f: u32) -> Option<u32> {
        if !self.flat(f) {
            return None;
        }
        let g = &mut self.doc.g;
        let (z, one) = (g.var(Var::Z), g.constant(1.0));
        let z = g.un(ink_core::graph::UnOp::Abs, z);
        let z = g.sub(z, one);
        Some(g.max(f, z))
    }

    /// `run` with each field that reads no z (a 2D profile, the same at
    /// every height) looked at in the plane: its reach its extent in x and
    /// y, every z (crate::plane_reach), whole where that extent is. What
    /// `run` finds is forgotten after, the view's caches back as they were.
    /// For a question asked the same at every height, of such fields alone.
    pub fn in_plane<R>(&mut self, run: impl FnOnce(&mut Self) -> R) -> R {
        let take = |v: &mut Self| {
            (std::mem::take(&mut v.reaches), std::mem::take(&mut v.boxes), std::mem::take(&mut v.tapes), std::mem::take(&mut v.found), std::mem::take(&mut v.outcomes))
        };
        let kept = take(self);
        self.plane = true;
        let out = run(self);
        self.plane = false;
        (self.reaches, self.boxes, self.tapes, self.found, self.outcomes) = kept;
        out
    }

    /// The way field `f` moves with each param of `ks` (Tape::param_way),
    /// found once per field and params: a claim measured at each moment of
    /// a motion asks it again at each (claims' estimated, of its bands),
    /// and it reads neither the moment nor any value.
    pub(crate) fn param_ways(&mut self, f: u32, ks: &[u32]) -> Vec<Option<i8>> {
        let key = (f, ks.to_vec());
        if let Some(w) = self.ways.get(&key) {
            return w.clone();
        }
        let tape = Tape::new(&self.doc.g, f);
        let w: Vec<Option<i8>> = ks.iter().map(|&k| tape.param_way(k)).collect();
        self.ways.insert(key, w.clone());
        w
    }

    /// `f` with the vars of `map` replaced (Graph::subst), built once per
    /// field and map: a claim measured at each moment of a motion builds
    /// the same at each (claims' estimated: its bands at their ends).
    pub(crate) fn subst(&mut self, f: u32, map: &[(Var, u32)]) -> u32 {
        let key = (f, map.to_vec());
        if let Some(&n) = self.substs.get(&key) {
            return n;
        }
        let n = self.doc.g.subst(f, map);
        self.substs.insert(key, n);
        n
    }

    /// `reach` of each of `fs`, found at once, in parallel.
    pub fn reach_all(&mut self, fs: &[u32]) {
        use rayon::prelude::*;
        let keys: Vec<(u32, Scope)> = fs.iter().map(|&f| (f, self.scope(&[f]))).collect();
        let mut todo: Vec<(u32, Scope)> = keys.into_iter().filter(|k| !self.reaches.contains_key(k)).collect();
        todo.sort_unstable_by_key(|k| k.0);
        todo.dedup_by_key(|k| k.0);
        let slabs: Vec<Option<u32>> = todo.iter().map(|k| self.slab(k.0)).collect();
        let (g, window, ctx) = (&self.doc.g, self.doc.domain(), &self.ctx);
        let found: Vec<Reach> = todo
            .par_iter()
            .zip(slabs)
            .map(|(k, slab)| match slab {
                Some(s) => crate::plane_reach(g, window, s, ctx),
                None => crate::reach_in(g, window, k.0, ctx),
            })
            .collect();
        self.reaches.extend(todo.into_iter().zip(found));
    }

    /// `reach` of `f` over `at` (a scope the view is not at: a motion's
    /// span), to find ahead of the question that asks it there, under
    /// `at_point(at)`, which forgets what it finds after: its key, and its
    /// search, which reads its tape and context alone, so a caller with
    /// many questions (interfere's moving pairs, each judged over the motion
    /// before it is measured: validity::outside) runs theirs side by side
    /// as it asks the ones before (`found_ahead`): its key, None where the
    /// view has it, and its search (`reach_search`).
    pub(crate) fn reach_ahead(&mut self, at: &Ctx, f: u32) -> Option<(u32, Scope)> {
        let whole = std::mem::replace(&mut self.ctx, at.clone());
        let key = (f, self.scope(&[f]));
        self.ctx = whole;
        (!self.reaches.contains_key(&key)).then_some(key)
    }

    /// The search `reach_ahead` keys.
    pub(crate) fn reach_search(&self, at: &Ctx, f: u32) -> ReachSearch {
        ReachSearch { tape: Arc::new(Tape::new(&self.doc.g, f)), window: self.doc.domain(), at: at.clone() }
    }

    /// The reach `reach_ahead` keyed, found: the next `reach` asking it
    /// takes it (one found before and not taken is forgotten).
    pub(crate) fn found_ahead(&mut self, key: (u32, Scope), r: Reach) {
        self.ahead = Some((key, r));
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
        let c = self.clip_to(f, r);
        let t = Arc::new(Tape::new(&self.doc.g, c));
        self.tapes.insert(key, t.clone());
        t
    }

    /// `f` clipped to reach `r`: in x and y alone where it is looked at in
    /// the plane (crate::plane_reach).
    fn clip_to(&mut self, f: u32, r: Reach) -> u32 {
        let z = if self.flat(f) { (f64::NEG_INFINITY, f64::INFINITY) } else { (r.lo[2], r.hi[2]) };
        self.doc.clip(f, [r.lo[0], r.lo[1], z.0], [r.hi[0], r.hi[1], z.1])
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
        ProofOptions { res, finest, depth: ink_kernel::scope::PARAM_SPLITS, tol, boxes: self.search_boxes }
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

    /// `prove`, of "inside a and not inside b" where `except` is (a, b), f
    /// being max(a, -b) (query::Except); with the tape it searched, f
    /// clipped to where it looked (a's reach, which holds all of f's
    /// material: f is no lower than a). What confirms its points reads
    /// that tape: f's own reach is not searched for (a search that would
    /// find no more of f than a's does, and on a larger tape).
    pub fn prove_except(&mut self, f: u32, except: Option<(u32, u32)>, res: Option<f64>, tol: f64) -> (Proof, bool, Arc<Tape>) {
        let opts = self.proof_options(res, tol);
        let Some(except) = except else {
            let (p, whole) = self.prove_with(f, &opts);
            return (p, whole, self.clipped(f));
        };
        let (a, b) = except;
        let r = self.reach(crate::claims::subset_reach(a, b));
        let clip = self.clip_to(f, r);
        let Reach { lo, hi, whole } = r;
        let tape = Arc::new(Tape::new(&self.doc.g, clip));
        let e = ink_kernel::query::Except { a: self.clipped(a), b: Arc::new(Tape::new(&self.doc.g, b)) };
        (ranged::prove_except(&tape, Some(&e), &self.ctx, lo, hi, &opts), whole, tape)
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
            if let Some((c, _)) = self.doc.g.enclosed_constant(n) {
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
                let e =
                    self.doc.g.enclosed_constant(n).map(|(_, b)| b).unwrap_or_else(|| ink_core::eval::eval_interval(&Tape::new(&self.doc.g, n), &self.ctx, [I::ENTIRE; 3], None));
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

    /// The thinnest wall of a field by rays at surface samples
    /// (feature::thickness): where it is and how thick, and the samples'
    /// spacing; a kernel answer, cached. None: the field is empty.
    /// Err(spacing, why): nothing was measured, as no sample fell on its
    /// surface (UNSAMPLED) or the call's time limit stopped its rays
    /// (STOPPED). Looked at in the plane (`in_plane`), its walls are
    /// sampled at z = 0, as at every height.
    pub fn thickness(&mut self, f: u32, res: Option<f64>) -> Option<Result<(V3, f64, f64), (f64, &'static str)>> {
        let (lo, hi, cells) = self.thickness_grid(f, res)?;
        let tape = self.clipped(f);
        let thinnest: Result<Option<(V3, f64)>, ink_kernel::stop::Stopped> =
            ink_kernel::cache::memo_region("thickness", &[&tape], &self.ctx, (lo, hi), |fp| _ = fp.f64(cells), |c| ink_kernel::feature::thickness(&tape, c, lo, hi, cells));
        let h = max_extent(lo, hi) / cells;
        Some(match thinnest {
            Ok(Some((at, t))) => Ok((at, t, h)),
            Ok(None) => Err((h, UNSAMPLED)),
            Err(_) => Err((h, STOPPED)),
        })
    }

    /// The spacing of `thickness`'s samples at `res`, nothing measured.
    /// None: the field is empty.
    pub fn thickness_res(&mut self, f: u32, res: Option<f64>) -> Option<f64> {
        let (lo, hi, cells) = self.thickness_grid(f, res)?;
        Some(max_extent(lo, hi) / cells)
    }

    /// Where `thickness` samples and how many cells across.
    fn thickness_grid(&mut self, f: u32, res: Option<f64>) -> Option<(V3, V3, f64)> {
        let (mut lo, mut hi, _) = self.bounds(f)?;
        if self.flat(f) {
            (lo[2], hi[2]) = (0.0, 0.0);
        }
        Some((lo, hi, res.map_or(128.0, |h| max_extent(lo, hi) / h)))
    }

    /// A field's enclosed voids, by proof (voids::voids), to `accuracy` (a
    /// void holds a point deeper than it): its cells halved down to `res`, else to
    /// its box's longest side / COMPONENT_CELLS. None: the field is empty.
    /// An unproven extent or more than COMPONENT_CELLS across at `res` is an error.
    /// Looked at in the plane (`in_plane`), a field has none: its air goes
    /// on along z.
    pub fn voids(&mut self, f: u32, res: Option<f64>, accuracy: f64) -> Result<Option<ink_kernel::voids::Voids>, InkError> {
        if self.flat(f) {
            return Ok(None);
        }
        self.reach(f).require_whole("shape")?;
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
        Ok(Some(ink_kernel::voids::voids(&tape, &self.ctx, lo, hi, res.unwrap_or(ext / COMPONENT_CELLS), accuracy)))
    }

    /// `ink components`: a field's shells (mesh::analyze: those whose
    /// triangles face inward are voids; an estimate), meshed at cells
    /// `res`, used as given; else at its box's longest side / 160, made
    /// again where it leaves material out (mesh::Unresolved: a wall or
    /// sheet thinner than the cells, or a part between them), each time at
    /// cells finer than the thinnest of what it left, while at most
    /// COMPONENT_CELLS across and COMPONENT_TRIANGLES predicted; and its
    /// voids by proof (`voids`), to the claims' default accuracy. None: the
    /// field is empty.
    pub fn components(&mut self, f: u32, res: Option<f64>) -> Result<Option<Components>, InkError> {
        let accuracy = self.accuracy();
        let Some(voids) = self.voids(f, res, accuracy)? else { return Ok(None) };
        let Some((lo, hi, _)) = self.bounds(f) else { return Ok(None) };
        let ext = max_extent(lo, hi);
        let tape = self.clipped(f);
        let mut h = res.unwrap_or(ext / 160.0);
        // Sampled within COMPONENT_TRIANGLES (made at twice the res where
        // its grid would take more), and searched with as many.
        use ink_kernel::mesh::Limit;
        let limit = Limit { bricks: Limit::bricks_for(COMPONENT_TRIANGLES), search: Limit::samples_for(COMPONENT_TRIANGLES) };
        let checked = |h: f64| {
            let (mut at, mut limit) = (h, limit);
            loop {
                if let Some((m, mut u)) = ink_kernel::mesh::mesh_checked(&tape, &self.ctx, lo, hi, at, limit, 0.0) {
                    if at > h {
                        u.refused = h;
                    }
                    break (m, u, at);
                }
                at *= 2.0;
                limit.search *= 0.25;
            }
        };
        let (mut m, mut unresolved, mut made_at) = checked(h);
        while res.is_none() && unresolved.thinnest > 0.0 && made_at == h {
            let k = ((h / unresolved.thinnest).log2().floor().max(0.0) + 1.0).min(60.0);
            let finer = h / 2f64.powf(k);
            if ext / finer > COMPONENT_CELLS || (m.tris.len().max(1) as f64) * 4f64.powf(k) > COMPONENT_TRIANGLES {
                break;
            }
            h = finer;
            (m, unresolved, made_at) = checked(h);
        }
        Ok(Some(Components { analysis: ink_kernel::mesh::analyze(&m), unresolved, voids, res: made_at }))
    }

    /// Mass moments of a field (weighted by a density field, if any) and its
    /// volume, at the cells `res` asks for, else measured to the mass's
    /// accuracy (query::measured); found once per field and density
    /// (`per_body`). None: it has no box.
    pub fn moments(&mut self, shape: u32, rho: Option<u32>, res: Option<f64>) -> Option<ink_kernel::query::Measured> {
        self.moments_to(shape, rho, res, false)
    }

    /// `moments`, measured on until its centre and inertia are as well
    /// known as its mass (query::Until::Moments): what reads them (a
    /// mechanism's dynamics, an inertia asked for). A mass alone does not
    /// pay for them.
    pub fn inertia(&mut self, shape: u32, rho: Option<u32>, res: Option<f64>) -> Option<ink_kernel::query::Measured> {
        self.moments_to(shape, rho, res, true)
    }

    fn moments_to(&mut self, shape: u32, rho: Option<u32>, res: Option<f64>, inertia: bool) -> Option<ink_kernel::query::Measured> {
        let found = self.per_body(if inertia { "inertia" } else { "mass" }, shape, rho.as_slice(), res, |v| {
            let rho = rho.map(|f| Tape::new(&v.doc.g, f));
            let accuracy = ink_kernel::query::ACCURACY;
            let asked = if inertia { ink_kernel::query::Until::Moments(accuracy) } else { ink_kernel::query::Until::Accuracy(accuracy) };
            let until = res.map_or(asked, ink_kernel::query::Until::Res);
            // A bounded body on its own field, as a published mass weighs it
            // (crate::materials::weighed); else the part in the window.
            let r = v.reach(shape);
            if r.whole && !v.flat(shape) {
                let t = Arc::new(Tape::new(&v.doc.g, shape));
                let (lo, hi, _) = ink_kernel::query::bounds(&t, &v.ctx, r.lo, r.hi)?;
                return Some(crate::materials::weighed(&t, rho.as_ref(), &v.ctx, lo, hi, until));
            }
            let (lo, hi, _) = v.bounds(shape)?;
            Some(ink_kernel::query::measured(&v.clipped(shape), &[rho.as_ref(), None], true, &v.ctx, lo, hi, until))
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

/// A reach's search (ModelView::reach_ahead): what `reach` finds there.
pub(crate) struct ReachSearch {
    tape: Arc<Tape>,
    window: (V3, V3),
    at: Ctx,
}

impl ReachSearch {
    pub(crate) fn run(&self) -> Reach {
        crate::reach_of(&self.tape, self.window, &self.at)
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

/// Why a value was not measured where the call's time limit stopped it
/// (ink_kernel::stop).
pub const STOPPED: &str = "the call's time limit passed before it was measured";

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

#[cfg(test)]
mod spanning_tests {
    use super::spanning;

    fn check(rows: &[Vec<f64>], one: bool, want: usize) {
        let refs: Vec<&[f64]> = rows.iter().map(|r| r.as_slice()).collect();
        let (picked, weights) = spanning(&refs, one);
        assert_eq!(picked.len(), want, "{picked:?}");
        for (r, w) in rows.iter().zip(&weights) {
            for j in 0..r.len() {
                let sum: f64 = picked.iter().zip(w).map(|(&p, w)| w * rows[p][j]).sum();
                let most = rows.iter().map(|r| r[j].abs()).fold(0.0, f64::max);
                assert!((sum - r[j]).abs() <= 1e-8 * most, "row {r:?} column {j}: {sum}");
            }
            if one {
                assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-8, "{w:?}");
            }
        }
    }

    /// Every row is the weighted sum of the rows picked, which are as many as the rows' directions.
    #[test]
    fn every_row_is_a_sum_of_the_rows_picked() {
        // One direction, its size and sign varying: a hoist's pull as it speeds up and slows.
        check(&[vec![0.0, -10.0, 0.0, 50.0], vec![0.0, -12.0, 0.0, 60.0], vec![0.0, 8.0, 0.0, -40.0]], false, 1);
        // The same beside something constant: the weights must sum to 1, which takes a second row.
        check(&[vec![0.0, -10.0, 0.0, 50.0], vec![0.0, -12.0, 0.0, 60.0], vec![0.0, 8.0, 0.0, -40.0]], true, 2);
        // A load turning in a plane, columns of very different size (N beside N*mm), one of them nothing.
        let turn: Vec<Vec<f64>> = (0..33).map(|i| (i as f64 * 0.2).sin_cos()).map(|(s, c)| vec![3.0 * c, 0.0, 3.0 * s, 4000.0 * s, 1e-3 * c]).collect();
        check(&turn, false, 2);
        check(&turn, true, 3);
        // No load at all: one run.
        check(&[vec![0.0; 3], vec![0.0; 3]], false, 1);
        // A small part in a third direction is a direction: it is not dropped.
        let mut tilted = turn.clone();
        tilted[7][1] = 1e-5;
        check(&tilted, false, 3);
    }
}

#[cfg(test)]
mod asked_tests {
    use super::*;
    use crate::{Claim, Value};

    fn view(src: &str) -> ModelView {
        static CORE: crate::Registry = crate::Registry::new(&[]);
        let doc = crate::Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), &[]).unwrap();
        let (ctx, _) = doc.ctx_at(&[], 0.0).unwrap();
        ModelView::new(doc, ctx)
    }

    /// A subset's search looks where a's material is, and its points are
    /// confirmed on the field it searched: max(a, -b) is never clipped to
    /// a reach of its own, a search that finds no more of it and, for a
    /// part turning over the motion (a belt's tooth round its loop), runs
    /// to its limit to say it is not bounded.
    #[test]
    fn a_subset_searches_no_reach_of_its_own_field() {
        let src = "(import \"std\")\n(t-span 0 1s)\n(body pin (shape (box 1 1 1)))\n(body ring (shape (cylinder-z 10 2)))\n\
                   (assembly root (instance p pin (at (move 5 0 1) (rot-z (* 360 (/ t 1s))))) (instance r ring (at)))\n\
                   (require inside (subset (shape-of root/p) (shape-of root/r)))";
        let mut v = view(src);
        v.set_motion(v.doc.t_span).unwrap();
        let r = crate::claims::run_asserts(&mut v, None, false, crate::claims::ClaimOptions::default()).unwrap();
        assert_eq!(r.asserts.len(), 1);
        let Claim::Subset(a, b) = v.doc.asserts[0].claim else { panic!("a subset") };
        let f = crate::claims::subset_field(&mut v.doc.g, a, b);
        assert!(v.tapes.keys().any(|k| k.0 == a), "a's own field is clipped to its reach");
        assert!(!v.tapes.keys().any(|k| k.0 == f), "the subset's field was clipped to a reach of its own");
    }

    /// A fact that reads what a solve sets (an estimate's guard) is no
    /// number yet, and is refused saying so; one of numbers is read.
    #[test]
    fn a_fact_reading_a_solved_param_is_refused() {
        let v = view("(def e (estimate-when 2 (> 1 0) \"why\"))\n(def n (+ 2 3))");
        let of = |name: &str| {
            let Some(Value::Field(value, _)) = v.doc.def(name).cloned() else { panic!("{name} is no field") };
            crate::Stated { what: name.into(), under: vec![], frame: None, value, unit: None, file: "t".into(), line: 1 }
        };
        assert_eq!(v.stated(&of("n")), Ok((5.0, None)));
        let refused = v.stated(&of("e")).expect_err("e reads the estimate's guard");
        assert!(refused.contains("which ink works out at a run"), "{refused}");
    }
}
