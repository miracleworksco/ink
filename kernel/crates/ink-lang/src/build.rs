//! Bodies, assemblies and relations built from their clauses: a body's
//! shape and attributes, an assembly's instances, and the relations between
//! their parts.

use super::*;

/// The kernel's own clauses of a relation (zone, claim, data) and of an
/// instance (its joint and limits): never `(optional ...)`, in either place.
const KERNEL_CLAUSES: &[&str] = &["zone", "claim", "data", "turn", "slide", "limits", "speed-limits"];

impl Doc {
    /// A body or an assembly named `n` from its clauses, its arguments
    /// taking `env`, in the names of the file it is written in; added to
    /// the document. `passed`: where the arguments were passed (the claims
    /// an assembly's joints add are there). A body's clauses read only its
    /// arguments: not the instances of an assembly it is built for.
    pub(super) fn build(
        &mut self,
        n: String,
        of: Option<String>,
        def: &ArgDef,
        env: Env,
        params: Vec<(String, String)>,
        passed: Option<(Rc<str>, u32, Scope)>,
        body: bool,
        at: &Sexp,
    ) -> R<()> {
        // The file that defines it: where it is written, or where the
        // arguments it is built with were passed.
        let file = passed.as_ref().map_or(def.file.clone(), |p| p.0.clone());
        let (name, defined) = (n.clone(), of.is_none());
        let building = if body {
            None
        } else {
            let asm = Assembly {
                name: n.clone(),
                params: params.clone(),
                of: of.clone(),
                instances: Vec::new(),
                coordinate_inputs: Vec::new(),
                frames: Vec::new(),
                frame_on: Vec::new(),
                relations: Vec::new(),
                clauses: Vec::new(),
                file: file.clone(),
                mass: None,
            };
            Some(Building { asm, later: Vec::new(), passed, mass: Vec::new(), read: Vec::new() })
        };
        let local = std::mem::replace(&mut self.local_assembly, building.as_ref().map(|b| Rc::from(b.asm.name.as_str())));
        let parent = std::mem::replace(&mut self.building, building);
        let args = std::mem::replace(&mut self.args, env);
        let saved = (std::mem::replace(&mut self.file, def.file.clone()), std::mem::replace(&mut self.scope, def.scope));
        let r = (|| {
            let clauses = self.expand_for(&def.clauses)?;
            if body {
                return self.body_clauses(n, of, params, &clauses, at).map(Some);
            }
            let later = clauses.iter().filter_map(|c| match c.list() {
                Some([h, i, ..]) if h.sym() == Some("instance") => i.sym().map(str::to_string),
                _ => None,
            });
            self.building.as_mut().unwrap().later = later.collect();
            self.assembly_clauses(&clauses).map(|_| None)
        })();
        (self.file, self.scope) = saved;
        self.args = args;
        self.local_assembly = local;
        let b = std::mem::replace(&mut self.building, parent);
        match r? {
            Some(body) => {
                let body = self.weighed_as_defined(body);
                self.bodies.push(Body { file, ..body })
            }
            None => {
                let b = b.unwrap();
                self.assemblies.push(b.asm);
                // Built again with other arguments, it places the bodies its
                // definition weighed.
                let get = |k: &str| b.mass.iter().find(|m| m.0 == k);
                if let (Some((_, mass, at)), true) = (get("mass"), defined) {
                    let saved = (std::mem::replace(&mut self.file, def.file.clone()), std::mem::replace(&mut self.scope, def.scope));
                    let r = self.assembly_mass(at, &name, mass, get("mass-params").map(|m| &m.1));
                    (self.file, self.scope) = saved;
                    r?;
                }
            }
        }
        Ok(())
    }

    /// A body from its clauses (expanded): one shape, attributes, frames.
    fn body_clauses(&mut self, n: String, of: Option<String>, params: Vec<(String, String)>, clauses: &[Sexp], at: &Sexp) -> R<Body> {
        let what = "body clauses are (shape expr), (attr key value), (frame name F) or (mark name (vec g r) ...)";
        let mut shape = None;
        let mut attrs: Vec<(String, Value)> = Vec::new();
        let mut density = None;
        let mut frames = Vec::new();
        let mut marks = Vec::new();
        for c in clauses {
            let parts = c.list().ok_or_else(|| self.err(c, what))?;
            match parts.first().and_then(|h| h.sym()) {
                Some("shape") => {
                    let [_, e] = parts else { return Err(self.err(c, "usage: (shape expr)")) };
                    if shape.is_some() {
                        return Err(self.err(c, "body has more than one (shape ...)"));
                    }
                    shape = Some(self.field(e, &mut self.arg_env())?);
                }
                Some("attr") => {
                    let [_, k, v] = parts else { return Err(self.err(c, "usage: (attr key value)")) };
                    let k = self.name_of(k)?;
                    if attrs.iter().any(|(a, _)| *a == k) {
                        return Err(self.err(c, format!("duplicate attr '{k}'")));
                    }
                    let e = v;
                    let mut v = self.eval(e, &mut self.arg_env())?;
                    let q = match k.as_str() {
                        "density" => Some(Quantity::Density),
                        "mass" => Some(Quantity::Mass),
                        "spin-inertia" => Some(Quantity::Inertia),
                        "tolerance" => Some(Quantity::Length),
                        _ => None,
                    };
                    if let Some(q) = q {
                        self.fits(&mut v, e, q, &format!("attr {k}"))?;
                    }
                    // A part's tolerance is the accuracy its claims are
                    // checked to (claims::finest): a constant length.
                    if k == "tolerance" {
                        let t = match &v {
                            Value::Field(f, _) => self.g.as_const(*f),
                            _ => None,
                        };
                        if !t.is_some_and(|t| t > 0.0 && t.is_finite()) {
                            return Err(self.err_hint(
                                c,
                                format!("attr tolerance expects a constant length above 0, got {e}"),
                                Some("(attr tolerance 0.2mm): how far the made part may differ from its shape".into()),
                            ));
                        }
                    }
                    // A look is read once, here: clients take look::Look,
                    // never the layout.
                    let v = look::read(&mut self.g, &k, v).map_err(|e| self.err(c, format!("attr {k}: {e}")))?;
                    // A density is g/mm^3: a number in another unit is
                    // refused here, not read a million times off.
                    match (k.as_str(), &v) {
                        ("density", Value::Field(f, _)) => density = Some((c, *f)),
                        ("density", other) => return Err(self.err(c, format!("attr density expects a number or a field, got a {}", other.kind()))),
                        ("material", Value::Data(d)) => {
                            if let Some(e) = materials::material_density_error(d.value()) {
                                return Err(self.err(c, format!("attr {k}: {e}")));
                            }
                        }
                        _ => {}
                    }
                    // A material is something physics can look up: a name
                    // in std/materials.ink (or data). Free text would make
                    // every study reading it fail far from here. Another
                    // spelling of a material is refused, naming it.
                    if let (true, Value::Str(m)) = (k == "material", &v)
                        && materials::lookup(m).is_none()
                    {
                        if materials::instead(m).is_some() {
                            let (what, hint) = materials::refused(m);
                            return Err(self.err_hint(c, format!("body '{}': {what}", of.as_deref().unwrap_or(&n)), Some(hint)));
                        }
                        let w = format!(
                            "{}:{}: body '{}': material '{m}' is not in std/materials.ink ({}). Name one of those, or give (data ...) with its keys. Free text goes in (attr material-note \"...\")",
                            self.file,
                            c.line,
                            of.as_deref().unwrap_or(&n),
                            materials::names().join(" ")
                        );
                        if !self.load_warnings.contains(&w) {
                            self.load_warnings.push(w);
                        }
                    }
                    attrs.push((k, v));
                }
                Some("frame") => frames.push(self.frame_clause(c, &frames)?),
                Some("mark") => marks.extend(self.mark_clause(c, parts)?),
                // A function's call is a shape written bare.
                head => {
                    let shape = head.filter(|h| self.function(h).is_some()).map(|h| format!("write (shape ({h} ...))"));
                    return Err(self.err_hint(c, what, self.top_form_hint(c, "a body").or(shape)));
                }
            }
        }
        let shape = shape.ok_or_else(|| self.err(at, format!("body '{n}' needs a (shape expr) clause")))?;
        if let Some((at, f)) = density {
            self.guard_magnitude(at, f, Some(shape), materials::DENSITY)?;
        }
        // The body's name is a construction name too (provenance: faces of this body).
        let shape = self.g.tag(&n, shape);
        self.published_mass(at, &n, of.as_deref(), shape, &mut attrs)?;
        Ok(Body { name: n, params, of, shape, attrs: Attrs::new(attrs), frames, marks, file: self.file.clone() })
    }

    /// `(mark NAME CURVE ...)`: curves drawn on the part, each CURVE `(vec
    /// g r)`, where field g is 0 inside field r < 0, or a vec of them.
    fn mark_clause(&mut self, c: &Sexp, parts: &[Sexp]) -> R<Vec<(String, NodeId, NodeId)>> {
        let usage = "usage: (mark NAME CURVE ...): curves drawn on the part (a sheet's bend lines: (mark bend (sheet-bend-lines S))), \
                     each CURVE (vec g r), where field g is 0 inside field r < 0, or a vec of them (none: (vec)); `ink slice --out` writes them on layer NAME";
        let [_, n, curves @ ..] = parts else { return Err(self.err_hint(c, "bad mark", Some(usage.into()))) };
        let name = self.name_of(n)?;
        let mut out = Vec::new();
        fn add(v: &Value, out: &mut Vec<(NodeId, NodeId)>) -> bool {
            let Value::Vec(l) = v else { return false };
            match l.as_slice() {
                [Value::Field(g, _), Value::Field(r, _)] => {
                    out.push((*g, *r));
                    true
                }
                l => l.iter().all(|x| add(x, out)),
            }
        }
        for e in curves {
            let v = self.eval(e, &mut self.arg_env())?;
            let mut got = Vec::new();
            if !add(&v, &mut got) {
                return Err(self.err_hint(e, format!("a mark's curve is (vec g r), two fields, or a vec of them, not a {}", v.kind()), Some(usage.into())));
            }
            out.extend(got.into_iter().map(|(g, r)| (name.clone(), g, r)));
        }
        Ok(out)
    }

    /// What an instance named `iname` places for `(CHILD ARG ...)`: the
    /// body or assembly CHILD built with the arguments' values (once per
    /// distinct values), named after this instance's path.
    fn instantiate(&mut self, c: &Sexp, iname: &str, child: &str, args: &[Sexp]) -> R<String> {
        let body = self.body(child).is_some();
        let Some(def) = self.arg_defs.get(child).cloned() else {
            return Err(self.err_hint(
                c,
                format!("'{child}' is {} without arguments: it takes none", if body { "a body" } else { "an assembly" }),
                Some("write the child as a name, or give it arguments: (body (name (p default) ...) ...), (assembly (name (p default) ...) ...)".into()),
            ));
        };
        let mut vals = Vec::new();
        for a in args {
            let mut env = self.arg_env();
            let v = self.eval(a, &mut env)?;
            if let Value::Field(f, _) = v
                && self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z))
            {
                return Err(self.err(a, "a body's or an assembly's arguments are values of params and t, not fields of x y z"));
            }
            vals.push((v, a.to_string()));
        }
        let inputs = &mut self.building.as_mut().unwrap().asm.coordinate_inputs;
        inputs.extend(vals.iter().filter_map(|(v, _)| if let Value::Field(f, _) = v { Some(*f) } else { None }));
        // A load with coordinates: each number the target's own instances
        // pass is its value plus a param at 0 (Coord); an argument left to its default is passed as that.
        let b = self.building.as_ref().unwrap();
        let coordinates = self.coordinates_of.as_deref() == Some(b.asm.name.as_str()) && b.asm.of.is_none() && self.import_depth == 0;
        if coordinates && vals.len() < def.params.len() {
            let (env, texts) = self.def_args(&def, child, vals, c)?;
            vals = env.into_iter().zip(texts).map(|((_, v), (_, text))| (v, text)).collect();
        }
        if coordinates {
            for ((v, _), (p, _)) in vals.iter_mut().zip(&def.params) {
                let Value::Field(f, ref said) = *v else { continue };
                let said = said.clone();
                let name = format!("{iname}.{p}");
                let param = |d: &mut Doc, what: &str| {
                    d.add_param(Param { name: format!("{what}{name}"), default: 0.0, tol: None });
                    (d.params.len() - 1, d.g.var(Var::Param(d.params.len() as u32 - 1)))
                };
                // A value that is its coordinate's own (a simulated motion:
                // Reader::as_coordinate), passed itself, is that coordinate
                // (once): what else reads it moves along it too. An argument
                // that otherwise reads such a value (twice it, its negative,
                // a sum), or passes it again, is no coordinate of its own:
                // passed as it is, it moves along that one.
                let own = self.coordinate_params.iter().find(|(x, _)| *x == f).map(|&(_, k)| k);
                let reads_own = self.coordinate_params.iter().any(|&(_, k)| self.g.depends_on(f, |x| x == Var::Param(k as u32)));
                if reads_own && (own.is_none() || self.coords.iter().any(|c| Some(c.param) == own)) {
                    continue;
                }
                let (k, e) = match own {
                    Some(k) => (k, self.g.var(Var::Param(k as u32))),
                    None => param(self, "coordinate "),
                };
                if self.states && self.g.depends_on(f, |x| x == Var::T) {
                    // At t0: e, at speed rate and acceleration acc.
                    let ((kr, rate), (ka, acc)) = (param(self, "rate of coordinate "), param(self, "acceleration of coordinate "));
                    let t0 = match self.state_t {
                        Some(k) => self.g.var(Var::Param(k as u32)),
                        None => {
                            self.add_param(Param { name: "moment of the state".into(), default: 0.0, tol: None });
                            self.state_t = Some(self.params.len() - 1);
                            self.g.var(Var::Param(self.params.len() as u32 - 1))
                        }
                    };
                    let t = self.g.var(Var::T);
                    let dt = self.g.sub(t, t0);
                    let (half, dt2) = (self.g.constant(0.5), self.g.mul(dt, dt));
                    let (r, a) = (self.g.mul(rate, dt), self.g.mul(acc, dt2));
                    let a = self.g.mul(half, a);
                    let ra = self.g.add(r, a);
                    let passed = self.g.add(e, ra);
                    *v = Value::Field(passed, said);
                    self.coords.push(Coord { name, param: k, value: f, jet: Some((kr, ka)), passed });
                } else {
                    let passed = if own.is_some() { f } else { self.g.add(f, e) };
                    *v = Value::Field(passed, said);
                    self.coords.push(Coord { name, param: k, value: f, jet: None, passed });
                }
            }
        }
        let (env, texts) = self.def_args(&def, child, vals, c)?;
        // An assembly also produces claims at the argument site. Two
        // equal sets of arguments in different imports may have different
        // time spans, even when the arguments equal the template defaults.
        let b = self.building.as_ref().unwrap();
        let passed = b.passed.clone().unwrap_or((self.file.clone(), c.line, self.scope));
        let key = (child.to_string(), format!("{env:?}"), (!body).then_some(passed.2));
        if let Some(n) = self.built.get(&key) {
            return Ok(n.clone());
        }
        let name = format!("{}/{iname}", b.asm.name);
        self.calls.push((child.to_string(), self.file.clone(), c.line));
        let r = self.build(name.clone(), Some(child.to_string()), &def, env, texts, Some(passed), body, c);
        self.calls.pop();
        r.map_err(|mut e| {
            e.called_from.push((child.to_string(), self.file.to_string(), c.line, c.col));
            e
        })?;
        self.built.insert(key, name.clone());
        Ok(name)
    }

    /// The values of the arguments of the body or assembly being read.
    pub(super) fn arg_env(&self) -> Env {
        self.args.clone()
    }

    fn assembly_clauses(&mut self, clauses: &[Sexp]) -> R<()> {
        let owned: Vec<String> = self.registry.packages.iter().flat_map(|p| p.clauses().iter().map(|h| format!(", ({h} ...)"))).collect();
        let what = format!("assembly clauses are (instance ...), (frame name F), (contact A B ZONE), (relate A B KIND){} or (attr mass M)", owned.concat());
        let what = what.as_str();
        for c in clauses {
            let parts = c.list().ok_or_else(|| self.err(c, what))?;
            match parts.first().and_then(|h| h.sym()) {
                Some("instance") => {
                    if parts.len() < 3 {
                        return Err(self.err_hint(c, "usage: (instance name child [WHERE [CFRAME] [MOTION] [(limits LO HI)] [(speed-limits LO HI)]])", Some(INSTANCE.into())));
                    }
                    let iname = self.name_of(&parts[1])?;
                    let asm = &self.building.as_ref().unwrap().asm;
                    if iname.contains('/') || asm.instances.iter().any(|i| i.name == iname) {
                        return Err(self.err(&parts[1], format!("instance name '{iname}' is invalid or duplicated")));
                    }
                    // (CHILD ARG ...): an assembly given arguments.
                    let (head, args) = match parts[2].list() {
                        Some([h, args @ ..]) => (h, args),
                        _ => (&parts[2], &[][..]),
                    };
                    let mut child = self.path_of(head)?;
                    if self.body(&child).is_none() && self.assembly(&child).is_none() {
                        let names = self
                            .bodies
                            .iter()
                            .filter(|b| b.of.is_none())
                            .map(|b| b.name.as_str())
                            .chain(self.assemblies.iter().filter(|a| a.of.is_none()).map(|a| a.name.as_str()));
                        return Err(self.err_hint(head, format!("no body or assembly named '{child}' (define it before the assembly)"), closest(&child, names)));
                    }
                    // In a load with coordinates, an instance the target places by its child's name alone passes that
                    // child's default arguments: they are its coordinates as written ones are (a pose held at the
                    // defaults is a pose, as `(hoist 0)` is).
                    let held = self.coordinates_of.is_some()
                        && self.import_depth == 0
                        && self.building.as_ref().is_some_and(|b| self.coordinates_of.as_deref() == Some(b.asm.name.as_str()) && b.asm.of.is_none())
                        && self.arg_defs.get(&child).is_some_and(|d| !d.params.is_empty());
                    if parts[2].list().is_some() || held {
                        child = self.instantiate(&parts[2], &iname, &child, args)?;
                    }
                    // (fact ...) clauses say nothing of where it is placed.
                    let mut facts = Vec::new();
                    let (at, joint, on, met) = self.instance_placement(c, &iname, &child, &parts[3..], &mut facts)?;
                    let b = self.building.as_mut().unwrap();
                    b.later.retain(|l| *l != iname);
                    b.asm.instances.push(Instance { name: iname, child, at, joint, on, met, facts });
                }
                Some("frame") => {
                    let frames = self.building.as_ref().unwrap().asm.frames.clone();
                    self.building.as_mut().unwrap().read.clear();
                    let f = self.frame_clause(c, &frames)?;
                    let read = std::mem::take(&mut self.building.as_mut().unwrap().read);
                    // Exactly an instance's frame: the language's frame-of (not a function of that name) and nothing more.
                    let its = match parts[2].head_args() {
                        Some(("frame-of", [_, its])) if self.scoped("frame-of").is_none() && read.len() == 1 => {
                            let mut part = self.building.as_ref().unwrap().asm.name.as_str();
                            for step in read[0].split('/') {
                                part = self
                                    .assembly(part)
                                    .or(self.building.as_ref().map(|b| &b.asm).filter(|a| a.name == part))
                                    .and_then(|a| a.instances.iter().find(|i| i.name == step))
                                    .map_or("", |i| i.child.as_str());
                            }
                            let frames = self.body(part).map(|b| &b.frames).or(self.assembly(part).map(|a| &a.frames));
                            its.sym().and_then(|its| frames?.iter().position(|(m, _)| m == its)).map(|k| k as u32)
                        }
                        _ => None,
                    };
                    let b = self.building.as_mut().unwrap();
                    b.asm.frame_on.push((f.0.clone(), read, its));
                    b.asm.frames.push(f);
                }
                Some("contact") => {
                    let usage = "(contact A B ZONE): instances A and B may overlap inside ZONE, a field in A's own coordinates (a thread in a nut: (shape-of screw-body thread))";
                    let [_, a, b, zone] = parts else { return Err(self.err_hint(c, "bad contact", Some(usage.into()))) };
                    let (pa, pb) = self.relation_ends(a, b, usage)?;
                    let zone = Zone { a: self.field(zone, &mut self.arg_env())?, b: None };
                    self.building.as_mut().unwrap().asm.relations.push(Relation { a: pa, b: pb, zone: Some(zone), kind: None });
                }
                Some("relate") => self.relate(c, parts)?,
                Some(h) if let Some((package, head)) = self.registry.clause(h) => {
                    let mut r = Reader::in_env(self, self.arg_env());
                    let kept = package.clause(&mut r, c, head, &parts[1..])?;
                    let reads = r.into_reads();
                    let (file, line, scope) = self.building.as_ref().unwrap().passed.clone().unwrap_or((self.file.clone(), c.line, self.scope));
                    self.building.as_mut().unwrap().asm.clauses.push(OwnedClause { package: package.name(), head, kept, reads, file, line, scope });
                }
                Some("attr") => {
                    let usage = "an assembly's attrs are its published (attr mass M) and (attr mass-params (read-json FILE)), FILE a JSON object of param values";
                    let [_, k, v] = parts else { return Err(self.err_hint(c, "usage: (attr key value)", Some(usage.into()))) };
                    let k = self.name_of(k)?;
                    if !matches!(k.as_str(), "mass" | "mass-params") {
                        return Err(self.err_hint(c, format!("an assembly has no attr '{k}'"), Some(usage.into())));
                    }
                    if self.building.as_ref().unwrap().mass.iter().any(|m| m.0 == k) {
                        return Err(self.err(c, format!("duplicate attr '{k}'")));
                    }
                    let e = v;
                    let mut v = self.eval(e, &mut self.arg_env())?;
                    if k == "mass" {
                        self.fits(&mut v, e, Quantity::Mass, "attr mass")?;
                    }
                    self.building.as_mut().unwrap().mass.push((k, v, c.clone()));
                }
                _ => return Err(self.err_hint(c, what, self.top_form_hint(c, "an assembly"))),
            }
        }
        Ok(())
    }

    /// The body or assembly a relation's end (an instance name, or a path
    /// into a sub-assembly, read so far) places.
    fn end_child(&self, path: &str) -> Result<String, String> {
        let mut segs = path.split('/');
        let first = segs.next().unwrap_or_default();
        let b = self.building.as_ref().unwrap();
        let mut child = b.asm.instances.iter().find(|i| i.name == first).map(|i| i.child.clone()).ok_or_else(|| format!("no instance '{first}'"))?;
        for seg in segs {
            child = self
                .assembly(&child)
                .and_then(|a| a.instances.iter().find(|i| i.name == seg))
                .map(|i| i.child.clone())
                .ok_or_else(|| format!("no instance '{seg}' in '{child}'"))?;
        }
        Ok(child)
    }

    /// The attributes of the body at a relation's end, as a record (`c`
    /// the relate form, for errors): numbers at load (params' values,
    /// which a run ranging one refuses the relation's claims over:
    /// AssertDef::fixed), strings, lists of them and data. The end is a
    /// body: an assembly's attributes are no one shape's. An attribute no
    /// record holds (a field of x y z, a function, a frame) is an error,
    /// never left out.
    fn end_attrs(&mut self, c: &Sexp, path: &str) -> R<Value> {
        let ch = self.end_child(path).map_err(|e| self.err(c, e))?;
        let Some(body) = self.body(&ch) else {
            return Err(self.err(c, format!("'{path}' is an assembly, and a relation reading its ends' attributes relates bodies: relate the one meshing by its path")));
        };
        let mut m = serde_json::Map::new();
        for (k, v) in body.attrs.clone().iter() {
            let x = self.attr_json(c, path, k, v)?;
            m.insert(k.clone(), x);
        }
        Ok(Value::Data(Data::new(serde_json::Value::Object(m))))
    }

    /// Attribute `k`'s value `v` (end_attrs) as data.
    fn attr_json(&mut self, c: &Sexp, path: &str, k: &str, v: &Value) -> R<serde_json::Value> {
        // One that is no number at load (a colour painted by x y z, a density solved from a published mass)
        // stays in the record as text saying so: a kind that reads it as a number is refused there, and one
        // that does not (a mesh reading "gear") is not stopped by it.
        let number = |me: &mut Self, f: NodeId| -> R<serde_json::Value> {
            let why = match me.load_number(c, f) {
                Ok(Some(x)) => {
                    return serde_json::Number::from_f64(x).map(serde_json::Value::Number).ok_or_else(|| me.err(c, format!("'{path}''s attribute '{k}' is not a finite number")));
                }
                Ok(None) => "a field of x y z or t",
                Err(_) => "it reads a solved value or sampled data",
            };
            Ok(serde_json::Value::String(format!("'{path}''s attribute '{k}' is no number a relation's kind can read: {why}")))
        };
        Ok(match v {
            Value::Field(f, _) => number(self, *f)?,
            Value::Point(p) => serde_json::Value::Array(p.iter().map(|&f| number(self, f)).collect::<R<_>>()?),
            Value::Str(t) => serde_json::Value::String(t.clone()),
            Value::Vec(xs) => serde_json::Value::Array(xs.iter().map(|x| self.attr_json(c, path, k, x)).collect::<R<_>>()?),
            Value::Data(d) => d.value().clone(),
            other => return Err(self.err(c, format!("'{path}''s attribute '{k}' is a {}, which a relation's kind cannot read as a record", other.kind()))),
        })
    }

    /// The shape of a relation's end (an instance name, or a path into a
    /// sub-assembly, read so far) in its own coordinates: its body's, or
    /// the union of its assembly's bodies.
    fn end_shape(&mut self, path: &str) -> Result<NodeId, String> {
        let child = self.end_child(path)?;
        // The child by its name, which is a path where it was built with
        // arguments (`pair/p`), not an instance path to resolve.
        if let Some(b) = self.body(&child) {
            return Ok(b.shape);
        }
        let id = self.identity();
        let items = self.flatten(&child, &id, &child, &[])?;
        let shapes: Vec<NodeId> = items.iter().map(|i| self.placed_field(i)).collect();
        Ok(self.reduce(BinOp::Min, &shapes))
    }

    /// The ends of a contact or relation: instance names, or paths into
    /// sub-assemblies, read so far.
    fn relation_ends(&self, a: &Sexp, b: &Sexp, usage: &str) -> R<(String, String)> {
        let (pa, pb) = (self.name_of(a)?, self.name_of(b)?);
        let asm = &self.building.as_ref().unwrap().asm;
        for (who, s) in [(&pa, a), (&pb, b)] {
            self.check_contact_path(asm, who).map_err(|m| self.err_hint(s, m, Some(usage.into())))?;
        }
        Ok((pa, pb))
    }

    /// `(relate A B KIND)`: KIND, a function, called with A's and B's
    /// placements in this assembly's coordinates, gives a `(relation ...)`:
    /// its zone is intended overlap, as a contact's; each of its claims is
    /// the value claim `relate:ASSEMBLY/A~B/NAME`, stated where the
    /// assembly's arguments were passed, else here (as limits are). An end
    /// written `(PATH FRAME)` is read at FRAME, a frame in its own
    /// coordinates: the kind gets its placement composed with FRAME and its
    /// shape in FRAME's coordinates, and its zone there is put back in the
    /// end's own.
    fn relate(&mut self, c: &Sexp, parts: &[Sexp]) -> R<()> {
        let usage = "(relate A B KIND): KIND a function of A's and B's frames (in this assembly's coordinates) and their shapes (each in \
                     its own coordinates) giving a (relation ...), as std's kinds do, e.g. (relate pinion wheel (gear-mesh 15 45 1)); an \
                     end (PATH FRAME) is read at FRAME, a frame in its own coordinates (a gear whose axis is not its z)";
        let [_, a, b, kind] = parts else { return Err(self.err_hint(c, "bad relate", Some(usage.into()))) };
        // An end's FRAME is in its part's own coordinates: a frame the part
        // declares, by name, or one written there.
        let end = |me: &mut Self, s: &Sexp| -> R<(Sexp, Option<Xf>)> {
            let Some(l) = s.list() else { return Ok((s.clone(), None)) };
            let [p, f] = l else { return Err(me.err_hint(s, "a relation's end is a path or (PATH FRAME)", Some(usage.into()))) };
            if let Some(name) = f.sym() {
                let child = me.end_child(&me.name_of(p)?).map_err(|m| me.err(p, m))?;
                let frames = me.frames_of(&child);
                let Some((_, xf)) = frames.iter().find(|(n, _)| n == name) else { return Err(me.err(f, Self::no_frame(&child, name, &frames))) };
                return Ok((p.clone(), Some(*xf)));
            }
            me.own_frame = true;
            let xf = me.frame_value(f, &mut me.arg_env(), "a relation end's frame");
            me.own_frame = false;
            Ok((p.clone(), Some(xf?)))
        };
        let ((a, la), (b, lb)) = (end(self, a)?, end(self, b)?);
        let (pa, pb) = self.relation_ends(&a, &b, usage)?;
        let f = match self.eval(kind, &mut self.arg_env())? {
            Value::Fn(f) => f,
            other => return Err(self.err_hint(kind, format!("KIND must be a function of two frames and two shapes, not a {}", other.kind()), Some(usage.into()))),
        };
        let read = |me: &mut Self, p: &str, l: &Option<Xf>| -> R<(Value, NodeId)> {
            let (f, s) = (me.frame_of(p, None).map_err(|m| me.err(c, m))?, me.end_shape(p).map_err(|m| me.err(c, m))?);
            let Some(l) = l else { return Ok((Value::Frame(f), s)) };
            let (inv, at) = (me.inverse(l), me.compose(l, &f));
            let Value::Field(s, _) = me.to_world(&Value::Field(s, None), &inv) else { unreachable!("a placed field is a field") };
            Ok((Value::Frame(at), s))
        };
        let ((fa, sa), (fb, sb)) = (read(self, &pa, &la)?, read(self, &pb, &lb)?);
        let before: HashSet<usize> = self.materialized.borrow().keys().copied().collect();
        let (moved, banded) = (self.load_t.len(), self.covered.len());
        // A kind of six arguments also reads its ends' bodies' attributes,
        // a record each (a gear's declared teeth and module, say).
        let takes = match &f {
            Callable::Closure(c) => c.params.len(),
            Callable::User(_, d) => d.params.len(),
            Callable::Builtin(_) => 4,
        };
        let mut args = vec![fa, fb, Value::Field(sa, None), Value::Field(sb, None)];
        if takes == 6 {
            args.push(self.end_attrs(c, &pa)?);
            args.push(self.end_attrs(c, &pb)?);
        }
        let r = match self.call(&f, args, kind)? {
            Value::Relation(r) => r,
            other => return Err(self.err_hint(kind, format!("KIND gave a {}, not a (relation ...)", other.kind()), Some(usage.into()))),
        };
        // What the kind found at load of a part's shape (its bounds) is the
        // part's at one value: its claims name the parts' relation at each.
        // Over a motion that is not a box's, which holds over every moment
        // at once: refused.
        if self.load_t.len() > moved {
            return Err(self.err_hint(
                kind,
                "KIND read a box of a part that moves (reads t): a box found at load holds over the whole motion, so it would relate the parts over all of it at once, not at each moment",
                Some("relate the parts where neither moves, or claim what the relation means directly".into()),
            ));
        }
        // A ranged param it fixed at load, or found a box over the band
        // of, is this relation's claims' to answer for, not the file's: a
        // run ranging it refuses them alone (AssertDef::fixed).
        let fixed: Vec<usize> = {
            let mut m = self.materialized.borrow_mut();
            let mut new: Vec<usize> = m.keys().copied().filter(|k| !before.contains(k) && self.has_range(*k)).collect();
            new.iter().for_each(|k| _ = m.remove(k));
            new.extend(self.covered.drain(banded..).map(|c| c.0));
            new.sort_unstable();
            new.dedup();
            new
        };
        let asm = self.building.as_ref().unwrap().asm.name.clone();
        let (file, line, scope) = self.building.as_ref().unwrap().passed.clone().unwrap_or((self.file.clone(), c.line, self.scope));
        let mut claims = Vec::new();
        for (n, v, lo, hi, tol) in &r.claims {
            let name = format!("relate:{asm}/{pa}~{pb}/{n}");
            if self.asserts.iter().any(|x| x.name == name) {
                return Err(self.err(c, format!("'{pa}' and '{pb}' are related twice with a claim '{n}'")));
            }
            let (l, h) = (self.g.bin_real(BinOp::Sub, *v, *lo), self.g.bin_real(BinOp::Sub, *hi, *v));
            let claim = Claim::Value(self.g.min(l, h), Cmp::Ge, 0.0);
            let text = format!("{} {pa}~{pb}: lo <= {n} <= hi, as value = min({n} - lo, hi - {n}) >= 0", r.kind);
            let terms = vec![(n.clone(), *v), ("lo".into(), *lo), ("hi".into(), *hi)];
            self.asserts.push(AssertDef {
                name: name.clone(),
                claim,
                text,
                bound: None,
                res: None,
                accuracy: *tol,
                quantity: None,
                file: file.clone(),
                scope,
                line,
                terms,
                fixed: fixed.clone(),
            });
            claims.push(name);
        }
        let site = format!("{}:{}", self.file, c.line);
        let kind = Kind { name: r.kind.clone(), skipped: r.zone.as_ref().map(|z| z.0.clone()), claims, data: r.data.clone(), clauses: r.clauses.clone(), site };
        let back = |me: &mut Self, z: NodeId, l: &Option<Xf>| {
            let Some(l) = l else { return z };
            let Value::Field(z, _) = me.to_world(&Value::Field(z, None), l) else { unreachable!("a placed field is a field") };
            z
        };
        let zone = r.zone.as_ref().map(|z| z.1).map(|z| Zone { a: back(self, z.a, &la), b: z.b.map(|zb| back(self, zb, &lb)) });
        self.building.as_mut().unwrap().asm.relations.push(Relation { a: pa, b: pb, zone, kind: Some(Rc::new(kind)) });
        Ok(())
    }

    /// CLAUSES with each `(optional CLAUSE)` replaced by what CLAUSE expands
    /// to (a `for`'s clauses) where a package this ink loads owns it, and
    /// left out where none does (a mesh's coupling, a stepper's stiffness:
    /// the kernel stands without them). A kernel clause (KERNEL_CLAUSES, of
    /// either kind, `for`'s included) in it is an error: it is never optional. One no loaded
    /// package owns warns where some loaded package reads clauses of this
    /// kind (HEADS, a package's list of them), since there it can only be a
    /// misspelling or a package this ink lacks, and dropping it would drop
    /// that package's check; where none does (the kernel's ink), nothing
    /// could read it.
    pub(super) fn optional_clauses(
        &mut self,
        clauses: &[Sexp],
        env: &mut Env,
        heads: fn(&dyn ext::Package) -> &'static [&'static str],
        usage: &str,
        instance: Option<&str>,
    ) -> R<Vec<Sexp>> {
        let mut out = Vec::with_capacity(clauses.len());
        for c in clauses {
            let inner = match c.list().unwrap_or(&[]) {
                [h, inner] if h.sym() == Some("optional") => inner,
                _ => {
                    out.push(c.clone());
                    continue;
                }
            };
            for x in self.expand_for_in(std::slice::from_ref(inner), env)? {
                let head = x.list().and_then(|l| l.first()).and_then(|h| h.sym());
                match head {
                    Some(h) if self.registry.packages.iter().any(|p| heads(*p).contains(&h)) => out.push(x.clone()),
                    Some(h) if !KERNEL_CLAUSES.contains(&h) && h != "optional" => {
                        if let Some(instance) = instance
                            && h.ends_with("-limits")
                        {
                            let assembly = &self.building.as_ref().expect("an instance belongs to an assembly").asm.name;
                            self.dropped_limits.entry(assembly.clone()).or_default().push((
                                instance.to_owned(),
                                UncheckedClaim {
                                    name: h.to_owned(),
                                    file: self.file.clone(),
                                    line: x.line,
                                    reason: "this build has no package for the optional joint limit",
                                    claim: "optional joint limit",
                                    declared_at: None,
                                },
                            ));
                        }
                        let loaded: Vec<&str> = self.registry.packages.iter().filter(|p| !heads(**p).is_empty()).map(|p| p.name()).collect();
                        if !loaded.is_empty() {
                            let w = format!(
                                "{}:{}: (optional ({h} ...)) is left out: no package this ink loads ({}) reads '{h}'. Check its spelling",
                                self.file,
                                x.line,
                                loaded.join(", ")
                            );
                            if !self.load_warnings.contains(&w) {
                                self.load_warnings.push(w);
                            }
                        }
                    }
                    _ => return Err(self.err_hint(&x, "(optional CLAUSE) is for a package's clause, not the kernel's", Some(usage.into()))),
                }
            }
        }
        Ok(out)
    }

    /// `(relation NAME [(zone WHAT ZA [ZB])] (claim NAME VALUE LO HI) ...
    /// [(data KEY VALUE ...)] [CLAUSE ...])`: what a kind gives `relate`.
    pub(super) fn eval_relation(&mut self, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        let usage = "(relation NAME [(zone WHAT ZA [ZB])] (claim NAME VALUE LO HI [(accuracy x)]) ... \
                     [(data KEY VALUE ...)] [CLAUSE ...] [(optional CLAUSE)]), the value of a kind for (relate A B KIND): overlap of A and B \
                     inside ZA (in A's coordinates) and ZB (in B's) is intended, WHAT saying what that leaves unchecked; each claim LO <= \
                     VALUE <= HI is a value claim (numbers, params, t, frames), allowing VALUE x past them with (accuracy x); data is shown with the \
                     pair; CLAUSE a package's, refused where no package this ink loads owns it, (optional CLAUSE) left out there";
        let Some((n, clauses)) = args.split_first() else { return Err(self.err_hint(e, "bad relation", Some(usage.into()))) };
        let kind = match &n.k {
            SK::Sym(s) | SK::Str(s) => s.clone(),
            _ => return Err(self.err_hint(n, "a relation's name is a name or a string", Some(usage.into()))),
        };
        let mut r = RelationValue { kind, zone: None, claims: Vec::new(), data: Vec::new(), clauses: Vec::new() };
        let mut data = false;
        let clauses = self.expand_for_in(clauses, env)?;
        let clauses = &self.optional_clauses(&clauses, env, |p| p.relation_clauses(), usage, None)?;
        for c in clauses {
            let parts = c.list().unwrap_or(&[]);
            let head = parts.first().and_then(|h| h.sym());
            if let Some((package, h)) = head.and_then(|h| self.registry.relation_clause(h)) {
                if r.clauses.iter().any(|x| x.head == h) {
                    return Err(self.err(c, format!("one ({h} ...) per relation")));
                }
                let mut rd = Reader::in_env(self, env.clone());
                let kept = package.clause(&mut rd, c, h, &parts[1..])?;
                let reads = rd.into_reads();
                r.clauses.push(OwnedClause { package: package.name(), head: h, kept, reads, scope: self.scope, file: self.file.clone(), line: c.line });
                continue;
            }
            match (head, &parts[1.min(parts.len())..]) {
                (Some("zone"), [what, za, zb @ ..]) if zb.len() <= 1 && r.zone.is_none() => {
                    let Value::Str(what) = self.eval(what, env)? else { return Err(self.err_hint(what, "a zone's WHAT is a string", Some(usage.into()))) };
                    let a = self.field(za, env)?;
                    let b = zb.first().map(|zb| self.field(zb, env)).transpose()?;
                    r.zone = Some((what, Zone { a, b }));
                }
                (Some("claim"), [name, v, lo, hi, tol @ ..]) if tol.len() <= 1 => {
                    let name = self.name_of(name)?;
                    if r.claims.iter().any(|x| x.0 == name) {
                        return Err(self.err(c, format!("two claims named '{name}'")));
                    }
                    let mut f = |x: &Sexp| -> R<NodeId> {
                        let f = self.field(x, env)?;
                        if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
                            return Err(self.err(x, "a relation's claims are about numbers (constants, params, t, frames), not fields of x y z"));
                        }
                        Ok(f)
                    };
                    let (v, lo, hi) = (f(v)?, f(lo)?, f(hi)?);
                    let tol = tol.first().map(|t| self.accuracy_clause(t, None, env)).transpose()?;
                    r.claims.push((name, v, lo, hi, tol));
                }
                (Some("data"), kvs) if kvs.len() % 2 == 0 && !data => {
                    data = true;
                    for kv in kvs.chunks(2) {
                        let k = self.name_of(&kv[0])?;
                        let v = self.eval(&kv[1], env)?;
                        r.data.push((k, v));
                    }
                }
                (Some(h), _) if !["zone", "claim", "data"].contains(&h) => {
                    return Err(self.err_hint(
                        c,
                        format!("bad relation clause: '{h}' is not the kernel's, nor from a package this ink loads. Check its name; a package's clause needs an ink that loads that package"),
                        Some(usage.into()),
                    ));
                }
                _ => return Err(self.err_hint(c, "bad relation clause (one zone, one data)", Some(usage.into()))),
            }
        }
        Ok(Value::Relation(Rc::new(r)))
    }

    /// `(fact WHAT [FRAME] VALUE)` in an instance of `child`.
    pub(super) fn stated(&mut self, c: &Sexp, child: &str) -> R<Stated> {
        let usage = "usage: (fact WHAT [WHERE] VALUE) in an instance: a number the design states about the part placed there; WHERE one of its frames, or a path to a part under it or to one of that part's frames (rotor, rotor/shaft)";
        let (what, place, value) = match &c.list().unwrap()[1..] {
            [what, value] => (what, None, value),
            [what, place, value] => (what, Some(place), value),
            _ => return Err(self.err_hint(c, "wrong arguments", Some(usage.into()))),
        };
        let what = self.name_of(what)?;
        let (mut under, mut frame, mut part) = (Vec::new(), None, child.to_string());
        let written = place.map(|p| self.name_of(p)).transpose()?.unwrap_or_default();
        let segments: Vec<&str> = written.split('/').filter(|s| !s.is_empty()).collect();
        for (k, segment) in segments.iter().enumerate() {
            let frames = self.body(&part).map(|b| &b.frames).or(self.assembly(&part).map(|a| &a.frames)).unwrap();
            let instances = self.assembly(&part).map_or(&[][..], |a| &a.instances[..]);
            let (is_frame, is_instance) = (k + 1 == segments.len() && frames.iter().any(|(m, _)| m == segment), instances.iter().any(|i| i.name == *segment));
            if is_frame && is_instance {
                return Err(self.err(place.unwrap(), format!("'{part}' has a frame and an instance named '{segment}': a fact cannot say which it is on")));
            } else if is_frame {
                frame = Some(segment.to_string());
            } else if let Some(i) = instances.iter().find(|i| i.name == *segment) {
                under.push(segment.to_string());
                part = i.child.clone();
            } else {
                let names = frames.iter().map(|(m, _)| m.as_str()).chain(instances.iter().map(|i| i.name.as_str()));
                return Err(self.err(place.unwrap(), format!("'{part}' has no frame or instance '{segment}'. {}", choices(segment, "Its frames and instances", names))));
            }
        }
        // Built so that no operation in it with no result is folded or dropped to a plain 0 (Graph::strict).
        self.g.strict = true;
        let v = self.eval(value, &mut self.arg_env());
        self.g.strict = false;
        let (field, unit) = self.fact_value(v?, value)?;
        Ok(Stated { what, under, frame, value: field, unit, file: self.file.clone(), line: c.line })
    }

    /// A fact's number and unit from value `v` of `at`.
    fn fact_value(&mut self, v: Value, value: &Sexp) -> R<(NodeId, Option<&'static str>)> {
        // No unit for what is no quantity, for a sum with a bare term, or for a temperature difference (a rating's
        // temperature is a point on a scale).
        // What is a quantity times a power of angle is that quantity, the angle in radians (units::Said::read).
        let (unit, times) = match &v {
            Value::Field(_, Some(said)) if said.explicit_when.is_none() && !std::ptr::eq(said.unit, &units::DIFFERENCE) => {
                said.read().map_or((None, 1.0), |(q, k)| (Some(q.unit()), k))
            }
            _ => (None, 1.0),
        };
        let mut field = self.as_field(v, value)?;
        if times != 1.0 {
            self.g.strict = true;
            let k = self.g.constant(times);
            field = self.g.mul(field, k);
            self.g.strict = false;
        }
        if self.g.depends_on(field, |v| matches!(v, Var::X | Var::Y | Var::Z | Var::T)) {
            return Err(self.err(value, "a fact is a number or params, not a field of x y z or t"));
        }
        Ok((field, unit))
    }

    /// What one entry of a part's `derives` says without building it (`formula` builds it): the fact it works
    /// out, and the facts its function's parameters name, which it reads.
    pub fn derived_reads(entry: &Value) -> Option<(String, Vec<String>)> {
        let Value::Vec(entry) = entry else { return None };
        let [Value::Str(what), Value::Fn(f), ..] = &entry[..] else { return None };
        let params = match f {
            Callable::User(_, d) => d.params.clone(),
            Callable::Closure(k) => k.params.clone(),
            _ => return None,
        };
        Some((what.to_string(), params.iter().map(|p| p.to_string()).collect()))
    }

    /// One entry of a part's `(attr derives (vec (vec NAME FN UNIT ...) ...))`, for a tool that works it out when
    /// asked (ink-ratings): the fact NAME is FN of the facts its parameters name, each read as the quantity of the
    /// UNIT written for it (`1N`, `1rpm`; a bare `1`: a pure number). FN is called once, here, on new variables,
    /// one a parameter (`place` names them), so its value is ink's own arithmetic with its units, ranges and numbers
    /// that are none; where FN gives `(vec VALUE OK WHY)`, VALUE holds only where OK is 1, and WHY says why not.
    pub fn formula(&mut self, entry: &Value, place: &str) -> Result<Formula, String> {
        let bad = || format!("{place}: (attr derives ...) is a vec of (vec NAME FN UNIT ...), FN a function of the facts its parameters name, with a unit for each");
        let Value::Vec(entry) = entry else { return Err(bad()) };
        let [Value::Str(what), Value::Fn(f), units @ ..] = &entry[..] else { return Err(bad()) };
        let (params, file, at) = match f {
            Callable::User(_, d) if !d.rest => (d.params.clone(), d.file.clone(), d.body.clone()),
            Callable::Closure(k) if !k.rest => (k.params.clone(), k.file.clone(), k.body.clone()),
            _ => return Err(bad()),
        };
        if params.len() != units.len() {
            return Err(bad());
        }
        let (mut reads, mut args) = (Vec::new(), Vec::new());
        for (param, unit) in params.iter().zip(units) {
            let said = match unit {
                Value::Field(_, Some(s)) if s.unit.quantity != Quantity::Other && !s.is_point() => {
                    Some(Rc::new(units::Said::worked(units::unit(s.unit.quantity.unit()).expect("a quantity's unit is a unit"), None, None, &file, at.line)))
                }
                Value::Field(_, None) => None,
                _ => return Err(bad()),
            };
            let k = self.add_param(Param { name: format!("{place}: {what} reads {param}"), default: 0.0, tol: None });
            reads.push((param.clone(), k, said.as_ref().map_or("", |s| s.unit.name)));
            args.push(Value::Field(self.g.var(Var::Param(k as u32)), said));
        }
        let said = |e: InkError| format!("{place}: '{what}' ({}:{}): {}", file, at.line, e.message);
        self.g.strict = true;
        let v = self.call(f, args, &at);
        self.g.strict = false;
        let (v, only) = match v.map_err(said)? {
            Value::Vec(parts) => match &parts[..] {
                [v, ok, Value::Str(why)] => (v.clone(), Some((self.as_field(ok.clone(), &at).map_err(said)?, why.to_string()))),
                _ => return Err(bad()),
            },
            v => (v, None),
        };
        let (value, unit) = self.fact_value(v, &at).map_err(said)?;
        // A fact's number is not known at load: a function that decided on one there decided on nothing.
        if let Some((param, _, _)) = reads.iter().find(|r| self.materialized(r.1).is_some()) {
            return Err(format!("{place}: '{what}' ({file}:{}) decides on '{param}' at load, before any fact is stated: use select, or (vec VALUE OK WHY)", at.line));
        }
        Ok(Formula { what: what.to_string(), reads, value, unit, only, file, line: at.line })
    }
}
