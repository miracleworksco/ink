//! Claims as written: `(assert ...)` read into a Claim, its options and
//! accuracy, and its integrands.

use super::*;

impl Doc {
    /// An assert's options, the clauses after its claim: `(res r)`, the
    /// cell size the claim needs (a constant above 0 and finite, as for
    /// `--res`), and `(accuracy x)`, the accuracy it is checked to, in its
    /// units (`accuracy_clause`), for the claims that have one. Returns them.
    /// The claim `claim` named `name`, with its clauses, as `form` states
    /// it: `(require NAME CLAIM ...)`.
    pub(super) fn state_claim(&mut self, form: &Sexp, name: &Sexp, claim: &Sexp, options: &[Sexp]) -> R<()> {
        let n = self.public(&self.name_of(name)?);
        if self.asserts.iter().any(|a| a.name == n) {
            return Err(self.err(name, format!("requirement '{n}' already exists")));
        }
        let (c, bound) = self.claim(claim)?;
        let (res, tol) = self.assert_options(&c, options)?;
        let quantity = self.accuracy_quantity(&c);
        self.asserts.push(AssertDef {
            name: n,
            claim: c,
            text: claim.to_string(),
            bound,
            res,
            accuracy: tol,
            quantity,
            file: self.file.clone(),
            scope: self.scope,
            line: form.line,
            terms: Vec::new(),
            fixed: Vec::new(),
        });
        Ok(())
    }

    fn assert_options(&mut self, c: &Claim, options: &[Sexp]) -> R<(Option<f64>, Option<f64>)> {
        const USAGE: &str = "a claim's options: (res r), the cell size the claim needs; (accuracy x), the accuracy it is checked to";
        let mut res = None;
        let mut accuracy = None;
        for o in options {
            match o.list() {
                Some([head, _]) if head.sym() == Some("accuracy") => {
                    if accuracy.is_some() {
                        return Err(self.err(o, "a requirement has more than one (accuracy ...)"));
                    }
                    if c.accuracy() == Accuracy::Exact {
                        return Err(self.err_hint(
                            o,
                            "this claim is decided exactly, and has no accuracy to state",
                            Some("a point claim is decided by its sign there, and a range compares values exactly. Write any slack into its bounds".into()),
                        ));
                    }
                    let q = self.accuracy_quantity(c);
                    accuracy = Some(self.accuracy_clause(o, q, &mut Vec::new())?);
                }
                Some([head, r]) if head.sym() == Some("res") => {
                    if res.is_some() {
                        return Err(self.err(o, "a requirement has more than one (res ...)"));
                    }
                    if matches!(c, Claim::Metric { metric, .. } if metric.shape.is_none()) {
                        return Err(self.err_hint(
                            o,
                            "(res r) is the cell size of a geometry claim, and this claim is not checked on cells",
                            Some("a package form's metric is solved at its form's own (res ...)".into()),
                        ));
                    }
                    // A constant: a claim's res is not a design parameter.
                    let f = self.field_of(r, Quantity::Length, "res", &mut Vec::new())?;
                    let h = self.g.as_const(f).ok_or_else(|| self.err(r, "expected a constant number"))?;
                    res = Some(self.check_cell(r, h)?);
                }
                Some([head, _]) if head.sym() == Some("tol") => {
                    return Err(self.err_hint(o, "(tol ...) is a param's tolerance. A claim states its accuracy: (accuracy x)", Some(USAGE.into())));
                }
                _ => return Err(self.err_hint(o, format!("unknown requirement option {o}"), Some(USAGE.into()))),
            }
        }
        Ok((res, accuracy))
    }

    /// Whether claim `c` is a shape's metric answered to a length accuracy
    /// on the shape (Package::shape_accuracy: a void count), not on its
    /// value.
    pub fn shape_accuracy(&self, c: &Claim) -> bool {
        matches!(c, Claim::Metric { metric, .. } if metric.shape.is_some() && self.registry.shape_accuracy(metric.metric))
    }

    /// What a claim's accuracy (and a comparison's value) is in, where
    /// known: a length for a claim about shapes or a thickness, a volume's
    /// volume, a metric's quantity. None: the model's own units (a value,
    /// an integral, a mean), or none at all (points and ranges).
    pub fn accuracy_quantity(&self, c: &Claim) -> Option<Quantity> {
        match c {
            _ if c.accuracy_is_length() || self.shape_accuracy(c) => Some(Quantity::Length),
            Claim::Metric { metric, .. } => self.registry.metric_quantity(metric.metric),
            _ => None,
        }
    }

    /// An `(accuracy x)` clause (an assert's, a joint's limits', a
    /// relation's claim's): how far the claim's answer may be off, a
    /// constant above 0 and finite, in quantity `q` where the claim's units
    /// are known. `(tol ...)` is a part's tolerance, a param's band.
    pub(super) fn accuracy_clause(&mut self, o: &Sexp, q: Option<Quantity>, env: &mut Env) -> R<f64> {
        const HINT: &str = "(accuracy x): how far the claim's answer may be off, in its units";
        let x = match o.list() {
            Some([head, x]) if head.sym() == Some("accuracy") => x,
            Some([head, _]) if head.sym() == Some("tol") => {
                return Err(self.err_hint(o, "(tol ...) is a param's tolerance. A claim states its accuracy instead", Some(HINT.into())));
            }
            _ => return Err(self.err_hint(o, format!("expected (accuracy x), not {o}"), Some(HINT.into()))),
        };
        let mut v = self.eval(x, env)?;
        if let Some(q) = q {
            self.fits(&mut v, x, q, "accuracy")?;
        }
        self.fits_change(&v, x, "an accuracy")?;
        let f = self.as_field(v, x)?;
        let t = self.g.as_const(f).ok_or_else(|| self.err(x, "expected a constant number"))?;
        if !(t > 0.0 && t.is_finite()) {
            return Err(self.err(x, format!("(accuracy x) expects an accuracy above 0 and finite, got {t}")));
        }
        Ok(t)
    }

    /// A claim, and its bound as a field of params (see `claim_at`).
    pub(super) fn claim(&mut self, c: &Sexp) -> R<(Claim, Option<NodeId>)> {
        let claim = self.claim_parse(c);
        let bound = self.bound.take();
        Ok((claim?, bound))
    }

    fn claim_parse(&mut self, c: &Sexp) -> R<Claim> {
        let usage = format!(
            "claims: (empty f) (disjoint a b) (subset a b) (within f (x0 y0 z0) (x1 y1 z1)) \
             (clearance a b d) (touches a b) (range f lo hi region) (volume f <|<=|>|>= v) \
             (value expr < v) (overlaps a b) (integral target g <|<=|>|>= v) \
             (mean target h [(weight w)] <|<=|>|>= v), and a package's metric of a form or a shape, (metric form <|<=|>|>= v): {}",
            self.registry.metrics()
        );
        let (head, a) = c.head_args().ok_or_else(|| self.err_hint(c, "a claim must be a list", Some(usage.clone())))?;
        let env = &mut Vec::new();
        let bad = |d: &Doc| d.err_hint(c, format!("malformed claim '{head}'"), Some(usage.clone()));
        // Bounds are numbers: constants and params (at their current values),
        // like load magnitudes and clearance distances.
        const BOUND: &str = "a claim's bound must be a number (constants and params), not a field of x y z or t";
        Ok(match (head, a) {
            ("empty", [f]) => match self.eval(f, env)? {
                Value::Region { certain, possible } => Claim::EmptyBounds { certain, possible },
                value => {
                    let mask = self.as_field(value, f)?;
                    match self.print_bounds.get(&mask) {
                        Some(region) => Claim::EmptyBounds { certain: region.certain, possible: region.possible },
                        None => Claim::Empty(mask),
                    }
                }
            },
            // A point overlaps a shape where it is inside it.
            ("overlaps", [p, q]) => match (self.shape_or_point(p, env)?, self.shape_or_point(q, env)?) {
                (Ok(p), Ok(q)) => Claim::Overlaps(p, q),
                (Err(at), Ok(f)) | (Ok(f), Err(at)) => Claim::Point { f, at, inside: true },
                (Err(_), Err(_)) => return Err(self.err(c, "overlaps takes a shape and a shape or a point, not two points")),
            },
            ("touches", [p, q]) => {
                let (p, q) = (self.field(p, env)?, self.field(q, env)?);
                Claim::Touches(p, q)
            }
            ("value", [f, op, v]) => {
                // The two sides are one quantity, as a sum's terms are: a torque is not under a force, nor a
                // pressure under a bare number that says no unit.
                let (value, bound) = (self.eval(f, env)?, self.eval(v, env)?);
                if matches!((&value, &bound), (Value::Field(..), Value::Field(..))) {
                    self.one_said("value", &[value.clone(), bound.clone()], c)?;
                }
                let f = self.as_field(value, f)?;
                if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
                    return Err(self.err(c, "value claims take numbers (params, t), not fields of x y z"));
                }
                let cmp = self.cmp_of(op)?;
                let b = self.as_field(bound, v)?;
                let bound_node = b;
                let b = self.claim_bound(b).ok_or_else(|| self.err(c, BOUND))?;
                self.value_bounds.insert((f, b.to_bits()), bound_node);
                Claim::Value(f, cmp, b)
            }
            // A point (std's probe) and a shape: the point outside it, or
            // inside it.
            ("disjoint", [p, q]) => match (self.shape_or_point(p, env)?, self.shape_or_point(q, env)?) {
                (Ok(p), Ok(q)) => Claim::Disjoint(p, q),
                (Err(at), Ok(f)) | (Ok(f), Err(at)) => Claim::Point { f, at, inside: false },
                (Err(_), Err(_)) => return Err(self.err(c, "disjoint takes a shape and a shape or a point, not two points")),
            },
            ("subset", [p, q]) => match (self.shape_or_point(p, env)?, self.field(q, env)?) {
                (Ok(p), q) => Claim::Subset(p, q),
                (Err(at), f) => Claim::Point { f, at, inside: true },
            },
            ("within", [f, lo, hi]) => {
                let f = self.field(f, env)?;
                // Corners of constants and params: the box stays a field of
                // its params, so a ranged one is proven over its band.
                let (lo, hi) = (self.corner(lo)?, self.corner(hi)?);
                let b = self.box_nodes(lo, hi);
                Claim::Subset(f, b)
            }
            ("clearance", [p, q, d]) => {
                let (p, q) = (self.field(p, env)?, self.field(q, env)?);
                let d = self.field_of(d, Quantity::Length, "clearance", env)?;
                let d = self.claim_bound(d).ok_or_else(|| self.err(c, "clearance distance must be a number (constants and params), not a field"))?;
                Claim::Clearance(p, q, d)
            }
            ("range", [f, lo, hi, region]) => {
                let (f, lo, hi, region) = (self.field(f, env)?, self.field(lo, env)?, self.field(hi, env)?, self.field(region, env)?);
                let below = self.g.bin_real(BinOp::Sub, f, lo);
                let above = self.g.bin_real(BinOp::Sub, hi, f);
                let violation = self.g.min(below, above);
                Claim::Range(self.g.max(region, violation))
            }
            ("volume", [f, op, v]) => {
                let f = self.field(f, env)?;
                let cmp = self.cmp_of(op)?;
                let v = self.field_of(v, Quantity::Volume, "volume", env)?;
                let bound = self.claim_bound(v).ok_or_else(|| self.err(c, BOUND))?;
                // The integral of 1 over the shape.
                let one = self.g.constant(1.0);
                Claim::Integral { parts: vec![(f, vec![one])], frames: vec![None], cmp, bound }
            }
            ("integral", [target, g, op, v]) => {
                let (parts, frames) = self.integrands(target, &[g], env)?;
                let cmp = self.cmp_of(op)?;
                let v = self.field(v, env)?;
                let bound = self.claim_bound(v).ok_or_else(|| self.err(c, BOUND))?;
                Claim::Integral { parts: parts.into_iter().map(|(s, mut g)| (s, g.remove(0))).collect(), frames, cmp, bound }
            }
            ("mean", [target, h, rest @ ..]) if matches!(rest.len(), 2 | 3) => {
                let (w, op, v) = match rest {
                    [op, v] => (None, op, v),
                    [w, op, v] => match w.list() {
                        Some([head, w]) if head.sym() == Some("weight") => (Some(w), op, v),
                        _ => {
                            return Err(self.err_hint(
                                w,
                                "expected (weight w)",
                                Some("(mean target h (weight w) < v): the mean of h weighted by w, e.g. (weight (attr density)) for the centre of mass".into()),
                            ));
                        }
                    },
                    _ => unreachable!(),
                };
                let one = Sexp { k: SK::Num(1.0, None, None), ..h.clone() };
                let (parts, frames) = self.integrands(target, &[w.unwrap_or(&one), h], env)?;
                if parts.iter().any(|(_, f)| f[0].len() != 1) {
                    return Err(self.err(c, "a mean's weight is one field"));
                }
                let cmp = self.cmp_of(op)?;
                let v = self.field(v, env)?;
                let bound = self.claim_bound(v).ok_or_else(|| self.err(c, BOUND))?;
                Claim::Mean { parts: parts.into_iter().map(|(s, mut f)| (s, f[0][0], f.remove(1))).collect(), frames, cmp, bound }
            }
            // After the kernel's kinds, so a package cannot shadow one.
            (m, [of, op, v, rest @ ..]) if self.registry.metric(m).is_some() => {
                // (excluding KEY ...): what the metric leaves out that the
                // claim leaves out on purpose.
                let mut excluding = Vec::new();
                let known = self.registry.exclusions(m);
                let usage = || {
                    let list: Vec<String> = known.iter().map(|(k, about)| format!("{k}: {about}")).collect();
                    if list.is_empty() { format!("'{m}' leaves nothing out") } else { format!("'{m}' may leave out {}", list.join("; ")) }
                };
                for x in rest {
                    let Some([head, keys @ ..]) = x.list() else { return Err(self.err_hint(x, "expected (excluding KEY ...)", Some(usage()))) };
                    if head.sym() != Some("excluding") || keys.is_empty() {
                        return Err(self.err_hint(x, "expected (excluding KEY ...)", Some(usage())));
                    }
                    for k in keys {
                        match known.iter().find(|e| Some(e.0) == k.sym()) {
                            Some(e) => excluding.push(e.0),
                            None => return Err(self.err_hint(k, format!("'{m}' leaves out nothing called {}", k.sym().unwrap_or("that")), Some(usage()))),
                        }
                    }
                }
                let metric = self.form_metric(m, of, env, "require")?;
                let cmp = self.cmp_of(op)?;
                let v = match self.registry.metric_quantity(m) {
                    Some(q) => self.field_of(v, q, m, env)?,
                    None => self.field(v, env)?,
                };
                let value = self.claim_bound(v).ok_or_else(|| self.err(c, BOUND))?;
                Claim::Metric { metric, cmp, value, excluding }
            }
            _ => return Err(bad(self)),
        })
    }

    /// A claim's bound: its value at load time; its field is returned with
    /// the claim (see `claim`), and is what claims are checked against
    /// (`claim_at`), so its params are not materialized, nor baked: an
    /// optimizer differentiates the bound's field with the claim.
    fn claim_bound(&mut self, f: NodeId) -> Option<f64> {
        self.bound = Some(f);
        if let Some(v) = self.g.as_const(f) {
            return Some(v);
        }
        if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z | Var::T)) {
            return None;
        }
        Some(self.load_eval(f))
    }

    /// An integral's parts: for each body of `target` (a body, an assembly,
    /// or an instance path: its placed bodies, in world space), its shape
    /// and the fields of each of `exprs` (a field, or a vec of them), read
    /// with `(attr NAME)` meaning that body's attribute, placed with it. A
    /// target that names none of those is a field, one part with no
    /// attributes; one with no bodies (an empty assembly) is one part of
    /// an empty shape, its fields read once with no attributes: each
    /// integral is 0, and how many there are is the fields'. Every body
    /// gives each expression the same number of fields (one integral per
    /// field): an error names two that differ. With each part, its body's
    /// own frame, where it is a placed body.
    fn integrands(&mut self, target: &Sexp, exprs: &[&Sexp], env: &mut Env) -> R<(Vec<(NodeId, Vec<Vec<NodeId>>)>, Vec<Option<Placed>>)> {
        let path = target
            .sym()
            .filter(|s| s.contains('/') || self.files[self.scope.file as usize].objects.contains_key(*s))
            .map(|s| self.resolve(s))
            .transpose()
            .map_err(|m| self.err(target, m))?;
        let mut items: Vec<(NodeId, String, String, Result<Attrs, &'static str>, Option<Placed>)> = match path {
            Some(path) => {
                let items = self.items(&path).map_err(|m| self.err(target, m))?;
                items
                    .into_iter()
                    .map(|it| {
                        let who = self.named(&it);
                        let (shape, attrs) = (self.placed_field(&it), self.placed_attrs(&it));
                        let frame = it.instance.map(|(_, at)| Placed { shape: it.shape, at });
                        (shape, it.label, who, Ok(attrs), frame)
                    })
                    .collect()
            }
            None => vec![(self.field(target, env)?, target.to_string(), format!("'{target}'"), Err("is a field, with no attributes"), None)],
        };
        if items.is_empty() {
            let empty = self.g.constant(1.0);
            items.push((empty, target.to_string(), format!("'{target}'"), Err("has no bodies, so no attributes"), None));
        }
        let mut out: Vec<(NodeId, Vec<Vec<NodeId>>)> = Vec::new();
        let mut frames = Vec::new();
        let mut first = String::new();
        for (shape, label, who, attrs, frame) in items {
            let outer = self.integrand_of.replace((who, attrs));
            let fields = exprs.iter().map(|e| self.fields(e, env)).collect::<R<Vec<_>>>();
            self.integrand_of = outer;
            let fields = fields?;
            match out.first() {
                None => first = label,
                Some((_, those)) => {
                    let count = |n: usize| if n == 1 { "1 field".to_string() } else { format!("{n} fields") };
                    if let Some(((e, f), f0)) = exprs.iter().zip(&fields).zip(those).find(|((_, f), f0)| f.len() != f0.len()) {
                        return Err(self.err_hint(
                            e,
                            format!("{e} is {} for '{first}' but {} for '{label}'", count(f0.len()), count(f.len())),
                            Some("each field is one integral, summed over the bodies: every body must give the same number".into()),
                        ));
                    }
                }
            }
            out.push((shape, fields));
            frames.push(frame);
        }
        Ok((out, frames))
    }

    /// A field, or each of a vec of them (at least one).
    pub(super) fn fields(&mut self, e: &Sexp, env: &mut Env) -> R<Vec<NodeId>> {
        match self.eval(e, env)? {
            Value::Field(f, _) => Ok(vec![f]),
            Value::Vec(v) if !v.is_empty() => v
                .iter()
                .map(|x| match x {
                    Value::Field(f, _) => Ok(*f),
                    other => Err(self.err(e, format!("expected fields, got a {} in the vec: {e}", other.kind()))),
                })
                .collect(),
            other => Err(self.err(e, format!("expected a field or a vec of fields, got a {}: {e}", other.kind()))),
        }
    }

    /// The claim with its bound evaluated at `ctx` (solved params at their
    /// solved values, wherever the solve is in the file): the one path for
    /// asserts and optimization requirements.
    pub fn claim_at(&self, c: &Claim, bound: Option<NodeId>, ctx: &Ctx) -> Claim {
        let at = |v: f64| bound.map_or(v, |b| self.scalar(b, ctx));
        match c {
            Claim::Empty(f) => Claim::Empty(*f),
            Claim::EmptyBounds { certain, possible } => Claim::EmptyBounds { certain: *certain, possible: *possible },
            Claim::Range(f) => Claim::Range(*f),
            Claim::Disjoint(p, q) => Claim::Disjoint(*p, *q),
            Claim::Subset(p, q) => Claim::Subset(*p, *q),
            Claim::Touches(p, q) => Claim::Touches(*p, *q),
            Claim::Point { f, at, inside } => Claim::Point { f: *f, at: *at, inside: *inside },
            Claim::Clearance(p, q, d) => Claim::Clearance(*p, *q, at(*d)),
            Claim::Metric { metric, cmp, value, excluding } => Claim::Metric { metric: metric.clone(), cmp: *cmp, value: at(*value), excluding: excluding.clone() },
            Claim::Value(f, cmp, v) => Claim::Value(*f, *cmp, at(*v)),
            Claim::Overlaps(p, q) => Claim::Overlaps(*p, *q),
            Claim::Integral { parts, frames, cmp, bound: v } => Claim::Integral { parts: parts.clone(), frames: frames.clone(), cmp: *cmp, bound: at(*v) },
            Claim::Mean { parts, frames, cmp, bound: v } => Claim::Mean { parts: parts.clone(), frames: frames.clone(), cmp: *cmp, bound: at(*v) },
        }
    }

    /// A construction name as written, matched to the names in the graph.
    /// Inside a namespaced file a name may be renamed `ns.name` (it is also
    /// a def there) while the tag that made it, in a family file, kept the
    /// plain name, or the other way round: the exact name wins, then the
    /// name without, then with, the current namespace.
    pub(super) fn construction_names(&self, name: &str) -> Vec<String> {
        let mut v = vec![name.to_string()];
        if let Some(ns) = &self.files[self.scope.file as usize].ns {
            match name.strip_prefix(&format!("{ns}.")) {
                Some(plain) => v.push(plain.to_string()),
                None => v.push(format!("{ns}.{name}")),
            }
        }
        v
    }

    pub(super) fn resolve_tag(&self, name: &str) -> Option<(u32, String)> {
        self.construction_names(name).into_iter().find_map(|n| self.g.tag_names.iter().position(|t| *t == n).map(|k| (k as u32, n)))
    }

    pub(super) fn cmp_of(&self, op: &Sexp) -> R<Cmp> {
        Ok(match op.sym() {
            Some("<") => Cmp::Lt,
            Some("<=") => Cmp::Le,
            Some(">") => Cmp::Gt,
            Some(">=") => Cmp::Ge,
            Some("==" | "=") => {
                return Err(self.err_hint(op, "comparison must be one of < <= > >=", Some("for equality, bound the difference: (value (abs (- a b)) <= tol)".into())));
            }
            _ => return Err(self.err(op, "comparison must be one of < <= > >=")),
        })
    }
}
