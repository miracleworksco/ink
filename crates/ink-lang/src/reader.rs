//! The reader: what a package reads its forms with (docs/CORE.md,
//! Layers: packages). It reads expressions, names and the document's
//! shared vocabulary (materials, tags, claims, metrics) in the
//! document's scope, and makes errors at a form's location. A package sees
//! the document only through it, read-only apart from what it may change:
//! the graph it builds fields in, and the params it defines. It records
//! what it hands out (`Reads`): that record, not the package, says what the
//! form depends on.

use super::*;

/// A package's view of the document while it reads one form.
pub struct Reader<'a> {
    doc: &'a mut Doc,
    env: Env,
    reads: Reads,
    /// Reading a joint's clause: whether its motion turns (Reader::motion).
    turn: Option<bool>,
}

impl<'a> Reader<'a> {
    pub(super) fn new(doc: &'a mut Doc) -> Self {
        Reader { doc, env: Vec::new(), reads: Reads::default(), turn: None }
    }

    /// One reading a package's value inside an expression, in its scope.
    pub(super) fn in_env(doc: &'a mut Doc, env: Env) -> Self {
        Reader { doc, env, reads: Reads::default(), turn: None }
    }

    /// One reading a joint's clause, in the assembly's scope: `turn`
    /// whether its motion turns.
    pub(super) fn in_joint(doc: &'a mut Doc, env: Env, turn: bool) -> Self {
        Reader { doc, env, reads: Reads::default(), turn: Some(turn) }
    }

    /// Reading a joint's clause (Package::joint_clauses): whether its
    /// motion is a turn (true) or a slide; None elsewhere.
    pub fn motion(&self) -> Option<bool> {
        self.turn
    }

    /// What it handed out.
    pub(super) fn into_reads(self) -> Reads {
        self.reads
    }

    /// Note that it looked up form `i`: what it makes depends on that one.
    fn read_form(&mut self, i: usize) {
        self.reads.forms.push(i);
    }

    /// The graph, to build fields in.
    pub fn graph(&mut self) -> &mut Graph {
        &mut self.doc.g
    }

    /// The derivative of `f` along `var`, through solves (Doc::derivative).
    pub fn derivative(&mut self, f: NodeId, var: Var) -> NodeId {
        self.doc.derivative(f, var)
    }

    /// The packages the document is read with (for their metrics' names).
    pub fn registry(&self) -> &'static Registry {
        self.doc.registry
    }

    // ---------------------------------------------------------------- expressions

    /// An expression as a field.
    pub fn field(&mut self, e: &Sexp) -> R<NodeId> {
        let f = self.doc.field(e, &mut self.env)?;
        self.reads.fields.push(f);
        Ok(f)
    }

    /// An expression as a number known at load time (constants and params).
    pub fn number(&mut self, e: &Sexp) -> R<f64> {
        let f = self.field(e)?;
        self.doc.load_value(f).ok_or_else(|| self.err(e, "expected a number known at load time (constants and params only)"))
    }

    /// An expression where the form takes quantity `q` (`what` names the
    /// slot in messages): a number in it written in another quantity's unit
    /// is refused (units::Said).
    pub fn field_of(&mut self, e: &Sexp, q: Quantity, what: &str) -> R<NodeId> {
        let f = self.doc.field_of(e, q, what, &mut self.env)?;
        self.reads.fields.push(f);
        Ok(f)
    }

    /// An expression's value where the form takes quantity `q` (`what`
    /// names the slot in messages), as `field_of` checks it.
    pub fn value_of(&mut self, e: &Sexp, q: Quantity, what: &str) -> R<Value> {
        let v = self.value(e)?;
        self.doc.fits(&v, e, q, what)?;
        Ok(v)
    }

    /// `(accuracy x)`: how far past its bounds a claim allows, in quantity
    /// `q`'s units (None: a plain number), as a joint's limits read it.
    pub fn accuracy(&mut self, e: &Sexp, q: Option<Quantity>) -> R<f64> {
        self.doc.accuracy_clause(e, q, &mut self.env)
    }

    /// A number known at load time where the form takes quantity `q`.
    pub fn number_of(&mut self, e: &Sexp, q: Quantity, what: &str) -> R<f64> {
        let f = self.field_of(e, q, what)?;
        self.doc.load_value(f).ok_or_else(|| self.err(e, "expected a number known at load time (constants and params only)"))
    }

    /// The cell size of a `(res h)` clause: above 0 and finite
    /// (Doc::cell_size).
    pub fn cell_size(&mut self, e: &Sexp) -> R<f64> {
        let h = self.number_of(e, Quantity::Length, "res")?;
        self.doc.check_cell(e, h)
    }

    /// An expression's value (a field, string, vector, data or function).
    pub fn value(&mut self, e: &Sexp) -> R<Value> {
        let v = self.doc.eval(e, &mut self.env)?;
        self.record_value(&v);
        Ok(v)
    }

    /// The fields in a value (a function's: those it captured).
    fn record_value(&mut self, v: &Value) {
        v.fields_into(&mut self.reads.fields, true)
    }

    /// The value at `ctx` of a field it was handed (or built from them),
    /// for checks and messages: a number the form keeps comes from `number`
    /// or `load_value`, which record that its params became one number.
    pub fn scalar(&self, f: NodeId, ctx: &Ctx) -> f64 {
        self.doc.scalar(f, ctx)
    }

    /// A name (a symbol).
    pub fn name(&self, e: &Sexp) -> R<String> {
        self.doc.name_of(e)
    }

    /// A form's clauses with `for` expanded, in order.
    pub fn clauses(&mut self, clauses: &[Sexp]) -> R<Vec<Sexp>> {
        self.doc.expand_for(clauses)
    }

    // ---------------------------------------------------------------- shared vocabulary

    /// A claim (as in `assert`), and its bound as a field of params.
    pub fn claim(&mut self, c: &Sexp) -> R<(Claim, Option<NodeId>)> {
        let (claim, bound) = self.doc.claim(c)?;
        match &claim {
            Claim::Metric { metric, .. } => self.record_metric(metric),
            c => self.reads.fields.extend(c.fields()),
        }
        self.reads.fields.extend(bound);
        Ok((claim, bound))
    }

    /// `(head of)`: a registered package's metric of a shape or of one of
    /// its forms, for the `user` (the form reading it) named in errors.
    pub fn metric(&mut self, head: &str, of: &Sexp, user: &str) -> R<FormMetric> {
        let m = self.doc.form_metric(head, of, &mut self.env, user)?;
        self.record_metric(&m);
        Ok(m)
    }

    /// Record what a metric reads: its shape, else its form.
    fn record_metric(&mut self, m: &FormMetric) {
        match m.shape {
            Some(f) => self.reads.fields.push(f),
            None => self.record_form(m),
        }
    }

    /// The earlier form of type `T` named `name` (a heat study a thermal
    /// study takes its temperatures from), with its public name (`ns.name`
    /// in a file imported under a namespace): the one name that finds it
    /// again, which is what to keep of it (its own may be a local form's
    /// too). This form then depends on everything that one read.
    pub fn form<T: Form>(&mut self, name: &str) -> Option<(String, &T)> {
        let name = self.doc.resolve(name);
        let i = self.doc.ext.iter().position(|e| e.name == name && (e.form.as_ref() as &dyn std::any::Any).is::<T>())?;
        self.read_form(i);
        (self.doc.ext[i].form.as_ref() as &dyn std::any::Any).downcast_ref::<T>().map(|f| (name, f))
    }

    /// Record a reference to the form whose metric it read.
    fn record_form(&mut self, m: &FormMetric) {
        if let Some(i) = self.doc.ext.iter().position(|e| e.package == m.package && e.name == m.form) {
            self.read_form(i);
        }
    }

    /// The (ranged lo hi) `name` is bound to (Doc::ranged_named).
    pub fn ranged_named(&self, name: &str) -> Option<(usize, (f64, f64))> {
        self.doc.ranged_named(name)
    }

    /// The index of the param named `name`.
    pub fn param_index(&mut self, name: &str) -> Option<usize> {
        let k = self.doc.param_index(name)?;
        let f = self.doc.g.var(Var::Param(k as u32));
        self.reads.fields.push(f);
        Some(k)
    }

    /// The piece name `name` means in the document's scope, as built with
    /// `(tag NAME ...)`; None if no piece is.
    pub fn tag(&self, name: &str) -> Option<String> {
        self.doc.resolve_tag(name).map(|(_, t)| t)
    }

    /// A material: a name from the library (std/materials.ink; symbol or
    /// string), or data with the property keys.
    pub fn material(&mut self, m: &Sexp) -> R<serde_json::Value> {
        self.doc.material_of(m, &mut self.env)
    }

    // ---------------------------------------------------------------- load time

    /// The context for load-time evaluation: params declared so far, at
    /// their overridden and solved values.
    pub fn load_ctx(&self) -> Ctx {
        self.doc.load_ctx()
    }

    /// The document's param providers so far (solves), compiled.
    pub fn providers(&self) -> Providers {
        self.doc.providers()
    }

    /// A field's value at load time, if it is a number there (it does not
    /// depend on x y z or t).
    pub fn load_value(&self, f: NodeId) -> Option<f64> {
        self.doc.load_value(f)
    }

    /// Record that these fields were evaluated at load time: their params
    /// are baked in, so derivatives through them are taken by differences.
    pub fn bake(&self, fields: &[NodeId]) {
        self.doc.bake(fields);
    }

    // ---------------------------------------------------------------- names

    /// Check that `name` is free for a form with results: `study-field`
    /// finds results by form name, so no two such forms share one.
    pub fn new_result_name(&self, at: &Sexp, name: &str) -> R<()> {
        let public = self.doc.public(name);
        match self.doc.ext.iter().find(|e| e.name == public && !e.form.results().is_empty()) {
            Some(e) => Err(self.err(at, format!("'{name}' already names a {}", e.head))),
            None => Ok(()),
        }
    }

    // ---------------------------------------------------------------- params

    /// Define a new param (a name not yet used), returning its index.
    pub fn new_param(&mut self, at: &Sexp, name: &str, default: f64) -> R<usize> {
        self.doc.check_new_value_name(at, name)?;
        self.doc.check_param_name(at, name, "the param")?;
        let key = self.doc.defined(name, at.line);
        self.doc.param_keys.insert(key, self.doc.params.len());
        self.doc.params.push(Param { name: name.to_string(), default, tol: None });
        Ok(self.doc.params.len() - 1)
    }

    /// Mark a param solved (a numerical method's answer): claims that read
    /// it are estimates (no accuracy), unless a view derives it over its
    /// ranged params as a certified root (ParamProvider::over_scope).
    pub fn solved(&mut self, k: usize) {
        self.doc.solved_params.insert(k);
    }

    /// A param a provider sets, named `name` as answers name it (no name a
    /// file can write): not a number until provided.
    pub fn derived(&mut self, name: &str) -> usize {
        self.doc.params.push(Param { name: name.to_string(), default: f64::NAN, tol: None });
        self.doc.params.len() - 1
    }

    /// `value` plus a param, at 0, that the sum carries as its coordinate's
    /// own: where the target's instance passes the sum itself (a load with
    /// coordinates), this param is that coordinate's, not one added to it,
    /// so whatever else reads the sum (a solve, a placement, another
    /// argument) moves along the coordinate as the argument does (a
    /// simulated motion, read through a def).
    pub fn as_coordinate(&mut self, value: NodeId, name: &str) -> NodeId {
        self.doc.params.push(Param { name: name.to_string(), default: 0.0, tol: None });
        let k = self.doc.params.len() - 1;
        let own = self.doc.g.var(Var::Param(k as u32));
        let sum = self.doc.g.add(value, own);
        self.doc.coordinate_params.push((sum, k));
        sum
    }

    // ---------------------------------------------------------------- once the file is read (Package::finish)

    /// Whether the form being read is the file run's (not an import's): a
    /// form about the world the model is in (its gravity) is the file run's
    /// alone, as `(t-span ...)` is.
    pub fn in_file_run(&self) -> bool {
        self.doc.import_depth == 0
    }

    /// The bodies `target` places, each with its own field, its placement
    /// and the instances above it (Doc::items), each of their fields
    /// recorded: shapes, placements, joints' drives and limits, attributes.
    /// `target` as the form's file writes it (Doc::resolve).
    pub fn items(&mut self, target: &str) -> Result<Vec<Item>, String> {
        let items = self.doc.items(&self.doc.resolve(target))?;
        for it in &items {
            self.reads.fields.push(it.shape);
            if let Some((_, xf)) = &it.instance {
                self.record_value(&Value::Frame(*xf));
            }
            it.attrs.iter().for_each(|(_, v)| self.record_value(v));
            for j in it.links.iter().filter_map(|l| l.joint.as_ref()) {
                self.reads.fields.extend(std::iter::once(j.drive).chain([j.limits, j.speed_limits].into_iter().flatten().flat_map(|(lo, hi)| [lo, hi])));
                j.clauses.iter().for_each(|c| self.record_clause(c));
            }
        }
        Ok(items)
    }

    /// The relations declared under `target` (Doc::relations_under), what
    /// their kinds' package clauses read and where they are placed recorded.
    pub fn relations_under(&mut self, target: &str) -> Vec<crate::PlacedRelation> {
        let rs = self.doc.relations_under(&self.doc.resolve(target));
        for r in &rs {
            r.kind.iter().flat_map(|k| &k.clauses).for_each(|c| self.record_clause(c));
            self.record_value(&Value::Frame(r.at));
        }
        rs
    }

    /// What a package's clause read reading it.
    fn record_clause(&mut self, c: &crate::OwnedClause) {
        self.reads.fields.extend(c.reads.fields.iter().copied());
        self.reads.forms.extend(c.reads.forms.iter().copied());
    }

    /// The clauses `package` owns placed under `target`
    /// (Doc::clauses_under), what each read reading it and where its
    /// assembly is placed recorded.
    pub fn clauses_under(&mut self, target: &str, package: &str) -> Vec<crate::PlacedClause> {
        let target = self.doc.resolve(target);
        let cs = self.doc.clauses_under(&target, package);
        for c in &cs {
            self.record_clause(&c.clause);
            self.record_value(&Value::Frame(c.at));
        }
        cs
    }

    // ---------------------------------------------------------------- in an assembly (Package::clause)

    /// An instance path of the assembly being read (an instance read so
    /// far, or a path into one: `arm/rod`), as written.
    pub fn instance_path(&self, e: &Sexp) -> R<String> {
        let p = self.doc.name_of(e)?;
        let asm = &self.doc.building.as_ref().ok_or_else(|| self.err(e, "an instance path is read in an assembly"))?.asm;
        self.doc.check_contact_path(asm, &p).map_err(|m| self.err(e, m))?;
        Ok(p)
    }

    /// What the package owning clause `head` kept of the clauses before
    /// this one in the assembly being read.
    pub fn kept(&self, head: &str) -> Vec<std::sync::Arc<dyn std::any::Any + Send + Sync>> {
        self.doc.building.as_ref().map_or(Vec::new(), |b| b.asm.clauses.iter().filter(|c| c.head == head).map(|c| c.kept.clone()).collect())
    }

    /// `e` as a placement (`what` names the slot in errors): a frame of
    /// params, t and local names, as an instance's WHERE is read.
    pub fn frame(&mut self, e: &Sexp, what: &str) -> R<Xf> {
        let mut env = self.env.clone();
        let f = self.doc.frame_value(e, &mut env, what)?;
        self.record_value(&Value::Frame(f));
        Ok(f)
    }

    /// `a` then `b` (Doc::compose).
    pub fn compose(&mut self, a: &Xf, b: &Xf) -> Xf {
        self.doc.compose(a, b)
    }

    /// A body, by name, as the form's file writes it.
    pub fn body(&self, name: &str) -> Option<&Body> {
        self.doc.body(&self.doc.resolve(name))
    }

    /// The earlier forms of type `T`, in document order (each recorded, as
    /// `form` records one).
    pub fn forms<T: Form>(&mut self) -> Vec<&T> {
        let found: Vec<usize> = (0..self.doc.ext.len()).filter(|&i| (self.doc.ext[i].form.as_ref() as &dyn std::any::Any).is::<T>()).collect();
        self.reads.forms.extend(&found);
        found.into_iter().filter_map(|i| (self.doc.ext[i].form.as_ref() as &dyn std::any::Any).downcast_ref::<T>()).collect()
    }

    /// The document's params' names, by index.
    pub fn param_names(&self) -> Vec<String> {
        self.doc.params.iter().map(|p| p.name.clone()).collect()
    }

    /// How to load the file run again (Doc::reloader).
    pub fn reloader(&self) -> Option<Reload> {
        self.doc.reloader()
    }

    /// Add a claim to the document's (names are unique).
    pub fn assert(&mut self, at: &Sexp, a: AssertDef) -> R<()> {
        if self.doc.asserts.iter().any(|x| x.name == a.name) {
            return Err(self.err(at, format!("assert '{}' already exists", a.name)));
        }
        self.doc.asserts.push(a);
        Ok(())
    }

    /// Keep a form the package made (a param provider), resting on what
    /// the reader has handed out so far.
    pub fn keep(&mut self, package: &'static str, head: &'static str, at: &Sexp, form: std::sync::Arc<dyn Form>) -> R<()> {
        let reads = self.reads.clone();
        self.doc.keep_form(package, head, at, form, reads)
    }

    // ---------------------------------------------------------------- errors

    /// Where `at` is written, `file:line`.
    pub fn site(&self, at: &Sexp) -> String {
        format!("{}:{}", self.doc.file, at.line)
    }

    /// An error at `at`, in the file being read.
    pub fn err(&self, at: &Sexp, msg: impl Into<String>) -> InkError {
        self.doc.err(at, msg)
    }

    pub fn err_hint(&self, at: &Sexp, msg: impl Into<String>, hint: Option<String>) -> InkError {
        self.doc.err_hint(at, msg, hint)
    }
}
