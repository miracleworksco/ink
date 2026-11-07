//! Extension forms: the one mechanism packages plug into a document with
//! (docs/CORE.md, Layers: packages).
//!
//! The core owns a fixed set of top-level forms (`CORE_FORMS`). Every other
//! head belongs to a package in the document's `Registry`: the core hands the
//! form to that package's `read`, keeps the `Form` it returns (with its
//! package and head, in document order), and never interprets it. Code that
//! knows a form's type reads it back with `doc.forms::<Study>()` or
//! `doc.form::<Study>(name)`; code that does not (`study-field`, claims and
//! objectives on a package's metrics) goes through the `Form` trait.
//!
//! A package's metrics are its claim kinds: `(max-stress books < 30)`
//! compares a metric of a form with a bound, and `(minimize (compliance
//! books))` optimises it. The package declares them (`Package::metrics`),
//! each form says which it gives (`Form::metrics`), and running the form
//! measures them (`Outcome::metric`).

use super::*;
use std::any::Any;
use std::sync::Arc;

/// The top-level forms the core reads itself.
pub const CORE_FORMS: &[&str] = &["import", "param", "def", "fn", "body", "assembly", "assert", "t-span"];

/// A package: a named set of top-level forms (`study` for sim), and the
/// metrics they give.
pub trait Package: Sync {
    /// Its name ("sim").
    fn name(&self) -> &'static str;
    /// The heads of the forms it owns.
    fn forms(&self) -> &'static [&'static str];
    /// The metrics its forms can give, each with what it is (`("min-frequency",
    /// "lowest natural frequency (Hz) of a study with (modes n)")`).
    fn metrics(&self) -> &'static [(&'static str, &'static str)] {
        &[]
    }
    /// The quantity a metric is, if it has a unit (max-stress: a pressure):
    /// a claim's bound on it may say its unit, and answers say it.
    fn metric_quantity(&self, _metric: &str) -> Option<Quantity> {
        None
    }
    /// What its metrics may leave out (`Excluded::key`), each with what it
    /// is: a claim may leave one out on purpose, `(excluding KEY)`.
    fn exclusions(&self) -> &'static [(&'static str, &'static str)] {
        &[]
    }
    /// Of its metrics, those of a shape, `(head f op v)`, not of one of
    /// its forms.
    fn shape_metrics(&self) -> &'static [&'static str] {
        &[]
    }
    /// Metric `metric` (one of `shape_metrics`) of shape `f`, proven
    /// bounded, at cell size `res` (None: its own).
    fn shape_metric(&self, _v: &mut crate::view::ModelView, metric: &str, _f: NodeId, _res: Option<f64>) -> Result<Measured, InkError> {
        Err(InkError::msg(format!("'{metric}' is no metric of a shape")))
    }
    /// Read one of its forms, `(head name clause ...)` (`args` is everything
    /// after the head), through the reader (the document's scope).
    fn read(&self, r: &mut Reader, form: &Sexp, head: &str, args: &[Sexp]) -> R<Arc<dyn Form>>;
    /// Once the file run is read: what it adds that rests on the whole
    /// document (claims about what the root assembly places), through the
    /// reader.
    fn finish(&self, _r: &mut Reader) -> R<()> {
        Ok(())
    }
    /// The heads of the values it gives inside expressions (`(reaction
    /// PATH)`): numbers it provides (derived params), which claims read.
    fn values(&self) -> &'static [&'static str] {
        &[]
    }
    /// One of its values, `(head arg ...)` (`args` after the head), read
    /// through the reader in the expression's scope.
    fn value(&self, r: &mut Reader, e: &Sexp, head: &str, _args: &[Sexp]) -> R<Value> {
        Err(r.err(e, format!("'{head}' gives no value")))
    }
    /// The heads of the assembly clauses it owns (`(HEAD ...)` inside
    /// `(assembly ...)`): each read by [`Package::clause`] in the assembly's
    /// scope and kept with the assembly, placed with each instance of it
    /// (Doc::clauses_under).
    fn clauses(&self) -> &'static [&'static str] {
        &[]
    }
    /// The heads of the clauses it owns in an instance with a motion
    /// (`(HEAD ...)` after the child, as `(limits LO HI)` is): each read by
    /// [`Package::clause`] in the assembly's scope, the motion's kind there
    /// (Reader::motion), and kept with the joint (Joint::clauses).
    fn joint_clauses(&self) -> &'static [&'static str] {
        &[]
    }
    /// The heads of the clauses it owns in a relation (`(HEAD ...)` inside
    /// `(relation NAME ...)`): each read by [`Package::clause`] in the
    /// relation's scope (the kind's), and kept with it (Kind::clauses).
    fn relation_clauses(&self) -> &'static [&'static str] {
        &[]
    }
    /// Read one of its clauses (an assembly's, a joint's or a relation's),
    /// `(head arg ...)` (`args` after the head), through the reader in its
    /// scope (an assembly's: its arguments, and the instances read so far,
    /// Reader::instance_path): what it keeps of it.
    fn clause(&self, r: &mut Reader, c: &Sexp, head: &str, _args: &[Sexp]) -> R<Arc<dyn Any + Send + Sync>> {
        Err(r.err(c, format!("'{head}' is no clause of {}", self.name())))
    }
    /// The commands it adds to front ends (command.rs).
    fn commands(&self) -> &'static [crate::command::Command] {
        &[]
    }
    /// What it adds to `ink check`'s answer on a model (claims it would
    /// make and does not, and why), after the core's.
    fn check(&self, _s: &mut crate::view::ModelView) -> serde_json::Map<String, serde_json::Value> {
        serde_json::Map::new()
    }
}

/// A package's parsed form.
pub trait Form: Any + Send + Sync {
    /// Its name, the second item of `(head name clause ...)`.
    fn name(&self) -> &str;
    /// The results it gives, read as fields with `(study-field name
    /// result)` and shown by clients as their specs say.
    fn results(&self) -> Vec<ResultSpec> {
        Vec::new()
    }
    /// The metrics it gives (of its package's `metrics`), which claims
    /// bound and objectives optimise.
    fn metrics(&self) -> Vec<&'static str> {
        Vec::new()
    }
    /// Why a metric it gives has no value, where that is known without
    /// running it (a study a fix moves has no compliance): its claims are
    /// unknown at once.
    fn unavailable(&self, _metric: &str) -> Option<String> {
        None
    }
    /// The params it sets, if it is a param provider (design's `solve`),
    /// compiled so that it works without the document.
    fn provider(&self, _doc: &Doc) -> Option<Arc<dyn ParamProvider>> {
        None
    }
    /// The equations whose root the params it sets are (each field 0 there,
    /// the params its provider owns their unknowns), if they are: their
    /// derivatives are the implicit function's (Doc::derivative).
    fn equations(&self) -> Option<&[NodeId]> {
        None
    }
    /// Run it with its own settings, as `opts` overrides them. Only forms
    /// with results or metrics run. What it reads of the document is what
    /// its package read for it (`Reads`): the fields it holds, and the
    /// forms it looked up through the reader.
    fn run(&self, _doc: &mut Doc, _ctx: &Ctx, _opts: &RunOpts) -> R<Box<dyn Outcome>> {
        Err(InkError::msg(format!("'{}' has nothing to run", self.name())))
    }
}

/// What a form that defines params by a method gives (design's `solve`:
/// equations solved for unknowns). The core applies every provider, each
/// after every one that sets a param it reads (`Doc::solves`), whenever it
/// makes a context; solves them at load time as each is read, so later
/// load-time code sees what they set (`Doc::solve_at_load`); and refuses
/// overrides of the params they own.
pub trait ParamProvider: Send + Sync {
    /// Its name, for messages.
    fn name(&self) -> &str;
    /// The params it sets.
    fn owns(&self) -> &[usize];
    /// Whether what it sets depends on t: it is provided again at each
    /// moment (over a motion, at each step).
    fn uses_t(&self) -> bool;
    /// Set its params at time t, starting from their current values; and,
    /// where what it set falls short of the accuracy it states (a
    /// measurement that did not settle), why: answers reading its params
    /// say so (`Inexact::unsettled`). An error where it can set nothing.
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, String>;
    /// Why `provide` would fail at time t, found as cheaply as it can be
    /// (`ink check` asks it of each solve on demand): by default, by
    /// providing.
    fn refusal(&self, params: &[f64], t: f64) -> Option<String> {
        self.provide(&mut params.to_vec(), t).err()
    }
    /// What it sets over the ranged params of `ctx` (`ctx.p_spans`, at
    /// moment `ctx.t`), where it reads some: its params derived from them
    /// (ink_core::eval::Derive), certified for every value of those it
    /// reads, each with its enclosure over all of them; or why that cannot
    /// be done (then assert refuses the ranges: never what it sets at one
    /// value labelled as covering them).
    fn over_scope(&self, _ctx: &Ctx) -> Result<OverScope, String> {
        Err("it sets its params at single values only".into())
    }
    /// Whether it is provided only where something reads what it sets (a
    /// field in use: Doc::want), not at load and at every context: a
    /// measurement nothing in the model's shape rests on (a published
    /// mass's density).
    fn on_demand(&self) -> bool {
        false
    }
    /// Whether it starts from its params' current values (a numerical
    /// solve's initial guess): then they are inputs to it too.
    fn warm_starts(&self) -> bool {
        true
    }
    /// The params what it sets rests on (by index), where it knows them
    /// more narrowly than the fields it was read from (Solve::reads): a
    /// ranged param among those makes a claim reading it unknown.
    fn reads(&self) -> Option<Vec<usize>> {
        None
    }
    /// The steps a claim reading what it sets is measured at over the
    /// moments `span` (params at `params`), where it chooses them from what
    /// it rests on; None: the uniform steps measured claims take.
    fn steps(&self, _params: &[f64], _span: (f64, f64)) -> Option<Result<Steps, String>> {
        None
    }
    /// Where it has no equations to bound what it sets (a simulated motion),
    /// what it sets at the middle of `span` (moments), an enclosure of each
    /// over the span, and of each one's rate along t there, from what it
    /// is: a solve on t reading it is then followed along t over those
    /// (follow.rs). None: it gives none, and such a solve is not followed.
    fn over_span(&self, _params: &[f64], _span: (f64, f64)) -> Option<Result<(Vec<f64>, Vec<ink_core::interval::I>, Vec<ink_core::interval::I>), String>> {
        None
    }
    /// What it set at t (`provide`, at `params`) is the real value there;
    /// where the steps over a motion read it with a part at rest that really
    /// moves (a spinning group read round), per param it owns, how far the
    /// real value may be from that reading at any moment of the step holding
    /// t, and the reading (the band's centre): a measured claim reading them
    /// holds where it holds over that band, fails where it fails at the real
    /// values, and is undecided between (motion::measured_at). Every other
    /// path reads the real value. None: no band.
    fn band(&self, _params: &[f64], _t: f64) -> Option<(Vec<f64>, Vec<f64>)> {
        None
    }
    /// Where what it sets is an estimate with its own measure of convergence
    /// (a value on a grid against the half grid's), per param it owns that
    /// band's reach either way, infinite where it has none (no bound: the
    /// real value may lie outside it): a value claim reading them holds
    /// where it holds over that band, fails where it fails over all of it,
    /// and is unknown between; any other claim reading them is unknown
    /// (claims::estimated). None: no band.
    fn error(&self, _params: &[f64]) -> Option<Vec<f64>> {
        None
    }
    /// What an answer reading what it sets should say of it beyond its name
    /// (a reaction summing what bodies it does not carry exert), where it is
    /// known without building anything.
    fn note(&self) -> Option<String> {
        None
    }
    /// The root of its equations (Form::equations) that what it set at
    /// (params, t) stands for, its unknowns' values, shown regular: one
    /// root about it, with a Jacobian of full rank there, so their
    /// derivatives are defined. Err: why not.
    fn root(&self, _params: &[f64], _t: f64) -> Result<Vec<f64>, String> {
        Err("it is not shown to be a regular root of its equations".into())
    }
    /// Whether what it sets is one root among others of its equations
    /// (Form::equations; a numerical solve): where it is on t, it is
    /// followed along t from its root at t = 0 (Follow), not found afresh
    /// at each moment.
    fn follows(&self) -> bool {
        false
    }
}

/// Steps a provider chose over a span (ParamProvider::steps), what reads
/// it measured at each one's middle and at the span's two ends: their
/// ends, in order from the span's start to its end; those it could not
/// choose so, as it says: `unchecked` (what it rests on jumps in them) and
/// `loose` (the step budget ran out first); and how it chose them.
#[derive(Clone, Debug)]
pub struct Steps {
    pub ends: Vec<f64>,
    pub unchecked: Vec<(f64, f64)>,
    pub loose: Vec<(f64, f64)>,
    pub basis: String,
}

/// A provider's params derived over ranged params (ParamProvider::over_scope).
pub struct OverScope {
    pub derive: Arc<dyn ink_core::eval::Derive>,
    /// Each param it sets, with its enclosure over every value of the
    /// ranged params it reads.
    pub bands: Vec<(u32, (f64, f64))>,
}

/// A document's param providers in the order they are applied (each after
/// every one setting a param it reads: `Doc::solves`), compiled: they work
/// without the document (they are Send and Sync, for viewers and workers).
#[derive(Clone, Default)]
pub struct Providers(pub Vec<Arc<dyn ParamProvider>>);

impl Providers {
    /// The params they set.
    pub fn owned(&self) -> impl Iterator<Item = usize> + '_ {
        self.0.iter().flat_map(|p| p.owns().iter().copied())
    }

    /// The provider setting param `k`.
    pub fn owner(&self, k: usize) -> Option<&dyn ParamProvider> {
        self.0.iter().find(|p| p.owns().contains(&k)).map(|p| p.as_ref())
    }
}

/// A client's overrides of how a form runs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RunOpts {
    /// Grid spacing (mm), instead of the form's own (res ...).
    pub res: Option<f64>,
    /// Relative tolerance of iterative solvers.
    pub tol: f64,
}

impl Default for RunOpts {
    fn default() -> Self {
        RunOpts { res: None, tol: 1e-8 }
    }
}

/// What running a form gave.
pub trait Outcome {
    /// Why it gave nothing to read, if it failed (a solve that broke down or
    /// did not reach its tolerance). Then no result, metric or sample is
    /// taken from it: `study-field` refuses it, and a claim on one of its
    /// metrics is unknown, with this reason (`ModelView::measure`).
    fn failure(&self) -> Option<String> {
        None
    }
    /// One of its results as computed (`None` if it does not give it).
    fn samples(&self, _result: &str) -> Option<Samples> {
        None
    }
    /// One of its results as a field's data (`study-field`): by default its
    /// samples carried a few cells past where they end, 0 beyond
    /// (`sampled_result`).
    fn field(&self, result: &str) -> Option<crate::grid::GridData> {
        let s = self.samples(result)?;
        Some(crate::results::sampled_result(s.values, s.dims, s.origin, s.h, 0.0))
    }
    /// Its mode k (from 1), a vibration mode's shape, as an outcome of its
    /// own (its displacement the shape, 1 mm at its peak, its stresses the
    /// shape's: per 1 mm of peak displacement), and what the mode is ("mode
    /// 2: 412.3 Hz, 1 mm peak"). Where it has no mode k, why (what modes it
    /// has).
    fn mode(&self, _k: usize) -> Result<(Box<dyn Outcome>, String), String> {
        Err("it has no modes".into())
    }
    /// One of its metrics (`None` if it does not give it).
    fn metric(&self, _metric: &str) -> Option<Measured> {
        None
    }
    /// How one of its metrics responds to the shape, for exact derivatives
    /// (`None`: the design loop takes differences).
    fn adjoint(&self, _metric: &str) -> Option<Adjoint<'_>> {
        None
    }
}

/// A result a form gives, as clients show it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResultSpec {
    pub name: &'static str,
    /// Its unit, for legends ("MPa", "mm"; "" for none).
    pub unit: &'static str,
    /// It is the magnitude of a vector whose x, y and z components are these
    /// results: a client may deform the part by that vector.
    pub vector: Option<[&'static str; 3]>,
    /// Its samples are a shape's field (negative inside: a solid a package
    /// designed): as a field (`study-field`) it grows past its grid by the
    /// distance to it, as an imported mesh's signed distance does, so it
    /// stays bounded. Else it is a quantity and holds its edge value there
    /// (GridData::holding).
    pub shape: bool,
}

impl ResultSpec {
    /// A scalar result, a quantity.
    pub const fn scalar(name: &'static str, unit: &'static str) -> Self {
        ResultSpec { name, unit, vector: None, shape: false }
    }
}

/// A result's samples on a regular grid: sample (i, j, k) at
/// origin + h (i, j, k), x fastest.
#[derive(Clone, Debug)]
pub struct Samples {
    pub origin: V3,
    pub h: f64,
    pub dims: [usize; 3],
    /// NaN where the result has no value (outside the material).
    pub values: Vec<f32>,
    /// How much of each sample's cell the result covers, 0 to 1: a cell's
    /// material fraction, 1 at a node that has a value, 0 where none.
    pub cover: Vec<f32>,
}

/// A metric's value as a package measured it. Package metrics are
/// numerical: a claim on one is an estimate (no accuracy), never guaranteed.
#[derive(Clone, Debug)]
pub struct Measured {
    /// The value, or why there is none: the run failed (`Outcome::failure`),
    /// or there was nothing to measure over (no stress away from the
    /// supports). A claim on a metric with no value is unknown.
    pub value: Result<f64, String>,
    /// Where it is reached (a peak), if at a point.
    pub at: Option<V3>,
    /// The grid spacing it was computed at.
    pub res: Option<f64>,
    pub warnings: Vec<String>,
    /// Values of the same metric that the measurement leaves out, and
    /// why (a study's stress within the support margin; a buckle it does
    /// not model): a claim on it says so, and is unknown where counting
    /// one would change its answer, unless it leaves that out on purpose
    /// (`(excluding KEY)`).
    pub excluded: Vec<Excluded>,
    /// What of the part its grid could not see, where it says (a gap its
    /// cells glue shut): the value is another part's, so a claim on it
    /// that holds is unknown, and a coarser grid that says so while a finer
    /// does not cannot check the finer's answer.
    pub unseen: Option<String>,
    /// Whether the value is proven, not measured on a grid (a count of
    /// voids by proof): a claim on it is exact, not an estimate, unless it
    /// leaves out on purpose what the proof left undecided.
    pub proven: bool,
}

/// A value a metric leaves out (`Measured::excluded`).
#[derive(Clone, Debug)]
pub struct Excluded {
    /// Its name, as a claim leaves it out on purpose (`(excluding KEY)`):
    /// one of its package's `exclusions`.
    pub key: &'static str,
    /// The metric's value, were it counted (infinite: unbounded, as a
    /// buckled member's stress).
    pub value: f64,
    pub at: Option<V3>,
    /// What it is and why it is left out.
    pub what: String,
    /// How to decide a claim it leaves undecided.
    pub decide: String,
}

/// A metric's exact derivative through the shape (Hadamard's formula):
/// dM/dp = -surface integral of w (df/dp) / |grad f| over the boundary of
/// `shape`, for every param that moves only that boundary.
pub struct Adjoint<'a> {
    pub shape: NodeId,
    pub weight: Box<dyn Fn(V3) -> f64 + Sync + 'a>,
    /// Every other field the metric depends on (the form's inputs besides
    /// its shape: materials, densities, where loads and supports act): the
    /// weight does not cover them, so a param they read has no exact
    /// derivative (the design loop takes differences).
    pub uncovered: Vec<NodeId>,
}

/// A package's metric of one of its forms, `(max-stress books)`, or of a
/// shape, `(min-thickness part)`: what a package claim compares with a
/// bound and an objective optimises.
#[derive(Clone, Debug)]
pub struct FormMetric {
    pub package: &'static str,
    /// The form's name (empty for a shape's).
    pub form: String,
    pub metric: &'static str,
    /// The shape, for a metric of a shape (Package::shape_metrics).
    pub shape: Option<NodeId>,
}

/// The packages a document can use, found by the heads of their forms. A
/// document is loaded with one (`Doc::load_file_in`): `ink::registry()` is
/// the kernel's, and a distribution may build its own with fewer or more
/// packages (a registry lives as long as the program: a static).
pub struct Registry {
    pub packages: &'static [&'static dyn Package],
}

impl Registry {
    /// A registry of these packages.
    pub const fn new(packages: &'static [&'static dyn Package]) -> Registry {
        Registry { packages }
    }

    /// The package owning `head`, and its spelling of it.
    pub fn package(&self, head: &str) -> Option<(&'static dyn Package, &'static str)> {
        self.packages.iter().find_map(|p| p.forms().iter().find(|h| **h == head).map(|h| (*p, *h)))
    }

    /// The package owning the assembly clause `head`, and its spelling of it.
    pub fn clause(&self, head: &str) -> Option<(&'static dyn Package, &'static str)> {
        self.owner(head, |p| p.clauses())
    }

    /// The package owning the joint clause `head` (Package::joint_clauses).
    pub fn joint_clause(&self, head: &str) -> Option<(&'static dyn Package, &'static str)> {
        self.owner(head, |p| p.joint_clauses())
    }

    /// The package owning the relation clause `head`
    /// (Package::relation_clauses).
    pub fn relation_clause(&self, head: &str) -> Option<(&'static dyn Package, &'static str)> {
        self.owner(head, |p| p.relation_clauses())
    }

    fn owner(&self, head: &str, heads: impl Fn(&dyn Package) -> &'static [&'static str]) -> Option<(&'static dyn Package, &'static str)> {
        self.packages.iter().find_map(|p| heads(*p).iter().find(|h| **h == head).map(|h| (*p, *h)))
    }

    /// The package giving the value `head` in expressions, and its spelling
    /// of it.
    pub fn value(&self, head: &str) -> Option<(&'static dyn Package, &'static str)> {
        self.packages.iter().find_map(|p| p.values().iter().find(|h| **h == head).map(|h| (*p, *h)))
    }

    /// Every package form head.
    pub fn heads(&self) -> impl Iterator<Item = &'static str> + '_ {
        self.packages.iter().flat_map(|p| p.forms().iter().copied())
    }

    /// The package giving metric `head`, and its spelling of it.
    pub fn metric(&self, head: &str) -> Option<(&'static dyn Package, &'static str)> {
        self.packages.iter().find_map(|p| p.metrics().iter().find(|m| m.0 == head).map(|m| (*p, m.0)))
    }

    /// The quantity metric `head` is, if it has a unit.
    pub fn metric_quantity(&self, head: &str) -> Option<Quantity> {
        self.metric(head).and_then(|(p, m)| p.metric_quantity(m))
    }

    /// What metric `head`'s package says its metrics may leave out.
    pub fn exclusions(&self, head: &str) -> &'static [(&'static str, &'static str)] {
        self.metric(head).map_or(&[], |(p, _)| p.exclusions())
    }

    /// Every package metric, for messages.
    pub fn metrics(&self) -> String {
        self.packages.iter().flat_map(|p| p.metrics().iter().map(|m| m.0)).collect::<Vec<_>>().join(" ")
    }

    /// The package command `name`.
    pub fn command(&self, name: &str) -> Option<&'static crate::command::Command> {
        self.commands().find(|c| c.name == name)
    }

    /// Every package command, in package order.
    pub fn commands(&self) -> impl Iterator<Item = &'static crate::command::Command> + '_ {
        self.packages.iter().flat_map(|p| p.commands().iter())
    }

    /// Every top-level form head, core first, for messages.
    pub fn top_forms(&self) -> String {
        CORE_FORMS.iter().copied().chain(self.heads()).collect::<Vec<_>>().join(", ")
    }
}

/// A form that provides params (a solve), as the document orders them.
#[derive(Clone, Debug)]
pub struct Solve {
    /// Its index in `Doc::ext`.
    pub ext: usize,
    /// The params it sets.
    pub owns: Vec<usize>,
    /// The other params it read (its equations: `Doc::form_inputs`),
    /// sorted.
    pub reads: Vec<usize>,
    /// Whether it depends on t: its equations read t (ParamProvider::uses_t,
    /// as it is made), or it reads what a solve on t sets, transitively
    /// (`Doc::set_solves`). It is solved again at each moment, and where it
    /// is followed along t its root rests on the moment, not only on what
    /// it reads (one turn later, the same reads).
    pub on_t: bool,
    /// ParamProvider::on_demand and warm_starts.
    pub on_demand: bool,
    pub warm_starts: bool,
    /// Where its params are another solve's rates (crate::rates): the
    /// solve of equations they rest on (index in `Doc::ext`). Where they
    /// are not defined it fails, refusing what reads them, never the load.
    pub rates: Option<usize>,
}

/// The order param providers are applied in: each after every one setting a
/// param it reads (what it sets is a function of those), in document order
/// where that leaves a choice (`(owns, reads)` of each, by position). Err:
/// a cycle, as positions, each reading a param the next sets (and the last
/// one the first's): they are one system of equations.
fn solve_order(solves: &[(&[usize], &[usize])]) -> Result<Vec<usize>, Vec<usize>> {
    let n = solves.len();
    // after[i]: the providers i reads from.
    let after: Vec<Vec<usize>> = (0..n).map(|i| (0..n).filter(|&j| j != i && solves[i].1.iter().any(|k| solves[j].0.contains(k))).collect()).collect();
    let mut placed = vec![false; n];
    let mut order = Vec::with_capacity(n);
    while order.len() < n {
        match (0..n).find(|&i| !placed[i] && after[i].iter().all(|&j| placed[j])) {
            Some(i) => {
                placed[i] = true;
                order.push(i);
            }
            None => {
                // Every one left reads from another one left: walk back to a
                // repeat.
                let mut path = vec![(0..n).find(|&i| !placed[i]).expect("one left")];
                loop {
                    let last = *path.last().expect("not empty");
                    let j = *after[last].iter().find(|&&j| !placed[j]).expect("reads one left");
                    if let Some(at) = path.iter().position(|&p| p == j) {
                        return Err(path.split_off(at));
                    }
                    path.push(j);
                }
            }
        }
    }
    Ok(order)
}

/// A form a package read, as the document keeps it.
pub struct ExtForm {
    pub package: &'static str,
    pub head: &'static str,
    /// Its name in the document: the form's own, `ns.name` in a file
    /// imported under a namespace (Doc::public).
    pub name: String,
    pub form: Arc<dyn Form>,
    /// What the package read to make it (`Reads`).
    pub reads: Reads,
    /// Where it is written: the file (as named to the loader, or an
    /// import's path), line and column.
    pub file: Rc<str>,
    pub line: u32,
    pub col: u32,
}

/// What a package read of the document while reading one form, as the
/// reader recorded it: every field it was handed (expressions, numbers,
/// regions, faces, claims and their bounds, bodies and params looked up),
/// and every earlier form it looked up (by name, or a metric of it). A
/// package sees the document only through the reader, so this is all that
/// the form's results and metrics can depend on: the core takes their
/// dependencies from it (`Doc::form_inputs`), and no package lists them.
#[derive(Clone, Debug, Default)]
pub struct Reads {
    pub fields: Vec<NodeId>,
    /// Indices in `Doc::ext`.
    pub forms: Vec<usize>,
}

pub(super) struct Provided {
    inputs: Vec<u64>,
    values: Vec<f64>,
    /// Why they fall short of the accuracy their provider states, if they do.
    pub(super) unsettled: Option<String>,
}

impl Doc {
    /// Apply one provider once for unchanged inputs and time (where it is
    /// on t: Solve::on_t; one that is not gives the same at every moment of
    /// a motion). The owned params are inputs too where it warm-starts:
    /// numerical solves use them as initial guesses, so a changed guess
    /// must still be solved and checked.
    pub fn provide_one(&self, i: usize, params: &mut [f64], t: f64) -> Result<(), String> {
        let s = &self.solves[i];
        let moment = s.on_t.then_some(t.to_bits());
        let guesses = if s.warm_starts { &s.owns[..] } else { &[] };
        let inputs: Vec<u64> = moment.into_iter().chain(s.reads.iter().chain(guesses).map(|&k| params[k].to_bits())).collect();
        if let Some(last) = self.provided.borrow().get(&s.ext)
            && last.inputs == inputs
        {
            for (&k, &v) in s.owns.iter().zip(&last.values) {
                params[k] = v;
            }
            return Ok(());
        }
        let unsettled = self.follow().provide(i, params, t)?;
        let values = s.owns.iter().map(|&k| params[k]).collect();
        self.provided.borrow_mut().insert(s.ext, Provided { inputs, values, unsettled });
        Ok(())
    }

    /// Hand an extension form to its package and keep what it reads.
    pub(super) fn top_ext(&mut self, package: &'static dyn Package, head: &'static str, form: &Sexp, args: &[Sexp]) -> R<()> {
        let mut r = Reader::new(self);
        let f = package.read(&mut r, form, head, args)?;
        let reads = r.into_reads();
        self.keep_form(package.name(), head, form, f, reads)
    }

    /// Keep a parsed form and order its derived parameters with the others.
    pub(super) fn keep_form(&mut self, package: &'static str, head: &'static str, form: &Sexp, f: Arc<dyn Form>, reads: Reads) -> R<()> {
        let provider = f.provider(self);
        let name = self.public(f.name());
        self.files[self.scope.file as usize].objects.insert(f.name().to_string(), name.clone());
        self.ext.push(ExtForm { package, head, name, form: f, reads, file: self.file.clone(), line: form.line, col: form.col });
        if let Some(p) = provider {
            let ext = self.ext.len() - 1;
            let reads = self.provider_reads(&self.ext[ext].reads, p.owns());
            let mut solves = self.solves.clone();
            // What a solve provided at every context reads is provided
            // with it, on demand or not.
            if !p.on_demand() {
                self.want_params(&reads);
            }
            solves.push(Solve { ext, owns: p.owns().to_vec(), reads, on_t: p.uses_t(), on_demand: p.on_demand(), warm_starts: p.warm_starts(), rates: None });
            self.set_solves(self.order_solves(solves, None).map_err(|m| self.err(form, m))?);
            self.varies_held().map_err(|m| {
                let hint = "write the solve before what differentiates through what it sets (a diff)";
                self.err_hint(form, format!("{head} '{}': {m}", self.ext[ext].form.name()), Some(hint.into()))
            })?;
            self.changes_no_use(ext).map_err(|(m, hint)| self.err_hint(form, m, Some(hint)))?;
            self.solve_at_load();
        }
        Ok(())
    }

    /// Solve every solve read so far, in order (`solves`), so the forms
    /// read after them see the solution: their unknowns' values become it,
    /// which warm-starts every load-time context (`load_ctx`). From the
    /// values so far, then warm-started and checked. One that fails is
    /// not an error yet: a solve read later may set what it reads (each
    /// prefix of the document need not solve), so the values stay as they
    /// were, and it is the error once the document is read (`unsolved`).
    /// Load time is t = 0, where a solve depending on t (directly or
    /// through what it reads: Solve::on_t) need not close (a linkage that
    /// assembles only later in its motion): it is solved again at each
    /// moment a view is made at (`ctx_at`), and its failure here is an
    /// error only if a value resting on it is used at load time
    /// (`unsolved_t`, checked once the document is read).
    pub(super) fn solve_at_load(&mut self) {
        let mut params = self.load_params();
        let _ = self.provide_each(&mut params, 0.0);
        let failed = self.provide_each(&mut params, 0.0);
        let first = failed.iter().find(|(i, _)| !self.solves[*i].on_t && self.solves[*i].rates.is_none());
        self.unsolved = first.map(|(i, m)| self.solve_error(*i, m.clone()));
        if self.unsolved.is_none() {
            for k in self.solves.iter().flat_map(|s| s.owns.clone()).collect::<Vec<_>>() {
                self.params[k].default = params[k];
            }
            // Following starts from the values at load time.
            *self.follow.borrow_mut() = None;
        }
        self.unsolved_t = failed.into_iter().map(|(i, m)| (i, self.solve_error(i, m))).collect();
    }

    /// The solves, in the order they are applied (`order_solves`), each on
    /// t that reads what one on t before it sets.
    pub(super) fn set_solves(&mut self, mut solves: Vec<Solve>) {
        for i in 0..solves.len() {
            let (before, rest) = solves.split_at_mut(i);
            rest[0].on_t |= before.iter().any(|b| b.on_t && b.owns.iter().any(|k| rest[0].reads.contains(k)));
        }
        self.solves = solves;
        *self.follow.borrow_mut() = None;
    }

    /// The solves compiled in their order, those on t followed along t
    /// from load time (Follow): what sets their params at each moment.
    pub fn follow(&self) -> Arc<Follow> {
        if let Some(f) = self.follow.borrow().as_ref() {
            return f.clone();
        }
        let solves = (0..self.solves.len())
            .map(|i| {
                let (s, e) = (&self.solves[i], &self.ext[self.solves[i].ext]);
                let label = match (e.form.as_ref() as &dyn Any).downcast_ref::<crate::rates::Rates>() {
                    Some(r) => r.label(),
                    None => format!("{} '{}'", e.head, e.form.name()),
                };
                let equations = e.form.equations().map(|eqs| eqs.iter().map(|&r| Tape::new(&self.g, r)).collect());
                (s.ext, e.form.provider(self).expect("a solve provides params"), equations, s.reads.clone(), self.solves[i].on_t, label)
            })
            .collect();
        let f = Arc::new(Follow::new(solves, self.load_params(), self.paths.clone()));
        *self.follow.borrow_mut() = Some(f.clone());
        f
    }

    /// Provide every solve's params in order (`solves`) at time t, each
    /// from the values so far, but one on demand that nothing wants. One
    /// that fails leaves its params as they were, and the rest go on. The
    /// failures: (index in `solves`, why).
    fn provide_each(&self, params: &mut [f64], t: f64) -> Vec<(usize, String)> {
        let mut failed = Vec::new();
        for (i, s) in self.solves.iter().enumerate() {
            if !self.provided_always(i) {
                continue;
            }
            let before: Vec<f64> = s.owns.iter().map(|&k| params[k]).collect();
            if let Err(m) = self.provide_one(i, params, t) {
                for (&k, v) in s.owns.iter().zip(before) {
                    params[k] = v;
                }
                failed.push((i, m));
            }
        }
        failed
    }

    /// Solve `i`'s failure `message`, as an error at its form.
    pub fn solve_error(&self, i: usize, message: String) -> InkError {
        let e = &self.ext[self.solves[i].ext];
        InkError { message, file: e.file.to_string(), line: e.line, col: e.col, hint: None, called_from: Vec::new() }
    }

    /// The solves (indices in `solves`) setting the params `ks`, and those
    /// setting what each of them reads, transitively: what values of `ks`
    /// rest on.
    pub fn solves_under(&self, ks: &[usize]) -> Vec<usize> {
        let mut out: Vec<usize> = Vec::new();
        let mut todo = ks.to_vec();
        let mut seen = HashSet::new();
        while let Some(k) = todo.pop() {
            if !seen.insert(k) {
                continue;
            }
            for (i, s) in self.solves.iter().enumerate() {
                if s.owns.contains(&k) && !out.contains(&i) {
                    out.push(i);
                    todo.extend(&s.reads);
                }
            }
        }
        out.sort_unstable();
        out
    }

    /// The solves the fields at `roots` rest on (`solves_under` the params
    /// they read).
    pub fn solves_read(&self, roots: &[NodeId]) -> Vec<usize> {
        self.solves_under(&self.params_read(roots))
    }

    /// Provide, from now on, every solve on demand that the fields at
    /// `roots` rest on (a density read: a claim's `(attr density)`, a
    /// `by-body`, what a command weighs: ModelView::want).
    pub fn want(&mut self, roots: &[NodeId]) {
        let ks = self.params_read(roots);
        self.want_params(&ks);
    }

    /// Whether solve `i` is provided at load and at every context: it is
    /// not on demand, or something wants what it sets.
    pub(super) fn provided_always(&self, i: usize) -> bool {
        let s = &self.solves[i];
        !s.on_demand || s.owns.iter().any(|k| self.wanted.contains(k))
    }

    /// `want` the params `ks`.
    pub(super) fn want_params(&mut self, ks: &[usize]) {
        for i in self.solves_under(ks) {
            if self.solves[i].on_demand {
                self.wanted.extend(self.solves[i].owns.iter().copied());
            }
        }
    }

    /// The params of solves on demand (a published mass's density) that
    /// the fields at `roots` read and nothing has wanted: values not
    /// measured here (measured where a command reads them).
    pub fn unmeasured(&self, roots: &[NodeId]) -> Vec<usize> {
        let mut ks: Vec<usize> = self.params_read(roots).into_iter().filter(|&k| self.unmeasured_param(k)).collect();
        ks.sort_unstable();
        ks
    }

    /// Whether param `k` is set by a solve on demand nothing has wanted.
    pub fn unmeasured_param(&self, k: usize) -> bool {
        !self.wanted.contains(&k) && self.solves.iter().any(|s| s.on_demand && s.owns.contains(&k))
    }

    /// A context at moment t, with the params `overrides` sets and every
    /// solve provided at t, in order; and the solves that failed there
    /// (index in `solves`, the error at its form), whose params keep their
    /// load-time values: nothing reading them is decided at this moment
    /// (ModelView::refusal). Err: an unknown param, or an override of a
    /// solved one.
    pub fn ctx_at(&self, overrides: &[(String, crate::Given)], t: f64) -> Result<(Ctx, Vec<(usize, InkError)>), String> {
        let providers = self.providers();
        let mut params: Vec<f64> = self.params.iter().map(|p| p.default).collect();
        for (k, v) in overrides.iter().map(|(k, g)| (k, &g.value)) {
            // A (ranged lo hi), by its def, is looked at at a value of its
            // range (a command answering claims refuses it: ink-cli's session).
            if self.param_index(k).is_none()
                && let Some((i, (lo, hi))) = self.ranged_named(k)
            {
                if !(lo..=hi).contains(v) {
                    return Err(format!("'{k}' is (ranged {lo} {hi}): a value of it lies in [{lo}, {hi}], not {v}"));
                }
                params[i] = *v;
                continue;
            }
            let i = self.param_index(k).ok_or_else(|| {
                let hint = match self.def(k) {
                    Some(v) if let Some((r, _)) = self.ranged_of(v) => {
                        let n = &self.params[r].name;
                        format!("'{k}' is another name for the (ranged lo hi) '{n}'. eval, check and render take one value, --param {n}=V, and assert and interfere range it")
                    }
                    Some(_) => format!("'{k}' is a value its file defines (a def: change it there), not a param"),
                    None => crate::closest(k, self.params.iter().map(|p| p.name.as_str())).unwrap_or_default(),
                };
                format!("unknown param '{k}'. {hint}").trim().to_string()
            })?;
            if let Some(p) = providers.owner(i) {
                return Err(format!("'{k}' is solved by '{}'; override its inputs instead", p.name()));
            }
            params[i] = *v;
        }
        let failed = self.provide_each(&mut params, t);
        Ok((Ctx::new(t, params), failed.into_iter().map(|(i, m)| (i, self.solve_error(i, m))).collect()))
    }

    /// Refuse the solve `ext`, just read, if it changes a value used at
    /// load time before it was solved: a load-time number (a loop count,
    /// an `if`, a bake, a measurement: `materialized`) would hide that it
    /// was solved, just as a use after the solve (`exact_at_load`). What
    /// it changes: what it sets, and what each solve reading that sets
    /// (after it in `solves`), transitively. Err: why, and a hint.
    fn changes_no_use(&self, ext: usize) -> Result<(), (String, String)> {
        let Some(at) = self.solves.iter().position(|s| s.ext == ext) else { return Ok(()) };
        let name = |s: &Solve| format!("{} '{}'", self.ext[s.ext].head, self.ext[s.ext].form.name());
        let new = &self.solves[at];
        // Each param it changes, with the solve setting it (None: this one).
        let mut changed: Vec<(usize, Option<&Solve>)> = new.owns.iter().map(|&k| (k, None)).collect();
        for s in &self.solves[at + 1..] {
            if s.reads.iter().any(|k| changed.iter().any(|(c, _)| c == k)) {
                changed.extend(s.owns.iter().map(|&k| (k, Some(s))));
            }
        }
        for (k, by) in changed {
            let Some(site) = self.materialized(k) else { continue };
            let p = &self.params[k].name;
            let what = match by {
                None => format!("{} sets '{p}'", name(new)),
                Some(s) => format!("{} changes '{p}' ({} solves it from what this one sets)", name(new), name(s)),
            };
            return Err((
                format!("{what}, whose value is used at load time at {site}. A solved value cannot shape the model at load time, whatever the solve order"),
                inexact::USE_IN_FIELDS.into(),
            ));
        }
        Ok(())
    }

    /// The fields an extension form depends on: those it read, and those of
    /// the earlier forms it looked up, transitively (a thermal study, the
    /// heat study it takes temperatures from and everything that one read).
    pub fn form_inputs(&self, e: &ExtForm) -> Vec<NodeId> {
        self.inputs_of(&e.reads)
    }

    /// The fields what a reader handed out depends on (see form_inputs).
    pub(super) fn inputs_of(&self, reads: &Reads) -> Vec<NodeId> {
        let mut fields = reads.fields.clone();
        let mut todo = reads.forms.clone();
        let mut seen = HashSet::new();
        while let Some(i) = todo.pop() {
            if seen.insert(i) {
                fields.extend(&self.ext[i].reads.fields);
                todo.extend(&self.ext[i].reads.forms);
            }
        }
        fields
    }

    /// The params a provider that read `reads` depends on, besides those
    /// it sets (`owns`), sorted.
    pub(super) fn provider_reads(&self, reads: &Reads, owns: &[usize]) -> Vec<usize> {
        let mut ks: Vec<usize> = self.params_read(&self.inputs_of(reads)).into_iter().filter(|k| !owns.contains(k)).collect();
        ks.sort_unstable();
        ks.dedup();
        ks
    }

    /// `solves` in the order they are applied (`solve_order`; `new`: the
    /// name of one not in `ext` yet, the last). Err: a cycle, named.
    pub(super) fn order_solves(&self, solves: Vec<Solve>, new: Option<&str>) -> Result<Vec<Solve>, String> {
        let pairs: Vec<(&[usize], &[usize])> = solves.iter().map(|s| (&s.owns[..], &s.reads[..])).collect();
        match solve_order(&pairs) {
            Ok(order) => Ok(order.into_iter().map(|i| solves[i].clone()).collect()),
            Err(cycle) => {
                let name = |i: usize| match self.ext.get(solves[i].ext) {
                    Some(e) => format!("{} '{}'", e.head, e.form.name()),
                    None => format!("solve '{}'", new.unwrap_or("?")),
                };
                let steps: Vec<String> = (0..cycle.len())
                    .map(|c| {
                        let (i, j) = (cycle[c], cycle[(c + 1) % cycle.len()]);
                        let k = solves[i].reads.iter().find(|k| solves[j].owns.contains(k)).expect("an edge of the cycle");
                        format!("{} reads '{}', which {} sets", name(i), self.params[*k].name, name(j))
                    })
                    .collect();
                Err(format!(
                    "{}: what each sets depends on the other, so they are one system of equations. Write one solve with all their unknowns and equations",
                    steps.join(", and ")
                ))
            }
        }
    }

    /// The document's forms of type `T`, in document order.
    pub fn forms<T: Form>(&self) -> impl Iterator<Item = &T> + '_ {
        self.ext.iter().filter_map(|e| (e.form.as_ref() as &dyn Any).downcast_ref::<T>())
    }

    /// The form of type `T` named `name`.
    pub fn form<T: Form>(&self, name: &str) -> Option<&T> {
        self.ext.iter().filter(|e| e.name == name).find_map(|e| (e.form.as_ref() as &dyn Any).downcast_ref::<T>())
    }

    /// The form of type `T` named `name`, or with no name the document's
    /// only one. `kind` (its head, "study") names it in the errors, which
    /// say which there are.
    pub fn pick<T: Form>(&self, kind: &str, name: Option<&str>) -> Result<&T, InkError> {
        self.picked::<T>(kind, name).map(|e| (e.form.as_ref() as &dyn Any).downcast_ref::<T>().expect("a form of its type"))
    }

    /// `pick`'s form as the document keeps it: its public name (`ns.name`
    /// for one imported under a namespace) is the one name that finds it
    /// again. A name two of them share is refused, not resolved to either.
    pub fn picked<T: Form>(&self, kind: &str, name: Option<&str>) -> Result<&ExtForm, InkError> {
        let all: Vec<&ExtForm> = self.ext.iter().filter(|e| (e.form.as_ref() as &dyn Any).is::<T>()).collect();
        let names: Vec<&str> = all.iter().map(|e| e.name.as_str()).collect();
        let there = if names.is_empty() { format!("the document has no {kind}") } else { format!("the document's {kind} names: {}", names.join(", ")) };
        match name {
            Some(n) => match all.iter().filter(|e| e.name == n).collect::<Vec<_>>()[..] {
                [e] => Ok(e),
                [] => Err(InkError::msg(format!("no {kind} named '{n}' ({there})"))),
                [a, b, ..] => Err(InkError {
                    hint: Some(format!("import the file one is in ({} or {}) under a namespace, (import \"FILE\" as NAME), and name it NAME.{n}", a.file, b.file)),
                    ..InkError::msg(format!("more than one {kind} is named '{n}' ({}:{} and {}:{})", a.file, a.line, b.file, b.line))
                }),
            },
            None if all.len() == 1 => Ok(all[0]),
            None if all.is_empty() => Err(InkError::msg(there)),
            None => Err(InkError::msg(format!("the document has more than one {kind} ({}): name one", names.join(", ")))),
        }
    }

    /// The param providers of its forms, compiled, in the order they are
    /// applied (`solves`).
    pub fn providers(&self) -> Providers {
        Providers(self.solves.iter().filter_map(|s| self.ext[s.ext].form.provider(self)).collect())
    }

    /// The extension form named `name` (of `package`, if given).
    pub fn ext_form(&self, package: Option<&str>, name: &str) -> Option<&ExtForm> {
        self.ext.iter().find(|e| package.is_none_or(|p| p == e.package) && e.name == name)
    }

    /// `(head of ...)`: a registered package's metric `head` of a shape
    /// (Package::shape_metrics: `of` is a field), else of one of its forms,
    /// defined before the `user` (assert, optimize) that reads it. The one
    /// reader of package metrics, for claims and objectives.
    pub(super) fn form_metric(&mut self, head: &str, form: &Sexp, env: &mut Env, user: &str) -> R<FormMetric> {
        let Some((package, metric)) = self.registry.metric(head) else {
            return Err(self.err_hint(form, format!("'{head}' is not a metric"), Some(format!("package metrics: {}", self.registry.metrics()))));
        };
        if package.shape_metrics().contains(&metric) {
            return Ok(FormMetric { package: package.name(), form: String::new(), metric, shape: Some(self.field(form, env)?) });
        }
        let name = self.path_of(form)?;
        let Some(e) = self.ext_form(Some(package.name()), &name) else {
            let heads = package.forms().join(" or ");
            return Err(self.err(form, format!("no {heads} named '{name}' (define it before the {user})")));
        };
        let gives = e.form.metrics();
        if !gives.contains(&metric) {
            let about = package.metrics().iter().find(|m| m.0 == metric).map_or("", |m| m.1);
            let list = if gives.is_empty() { "nothing".to_string() } else { gives.join(", ") };
            return Err(self.err_hint(form, format!("{} '{name}' does not give '{metric}'", e.head), Some(format!("{metric} is {about}; {} '{name}' gives {list}", e.head))));
        }
        Ok(FormMetric { package: package.name(), form: name, metric, shape: None })
    }
}
