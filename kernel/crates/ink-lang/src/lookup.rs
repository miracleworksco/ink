//! Names: what a name written in a document means where it is written
//! (params, defs, functions, imports, std), and the checks a new name passes.

use super::*;

impl Doc {
    /// The span stated by this run of a file (possibly after the claim),
    /// falling back to the run file's motion.
    pub fn span_of(&self, scope: Scope) -> Option<(f64, f64)> {
        self.files[scope.file as usize].span.or(self.t_span)
    }

    /// Each imported motion, with the source path and namespace it belongs to.
    pub fn import_spans(&self) -> impl Iterator<Item = (&str, Option<&str>, (f64, f64))> {
        self.files.iter().skip(2).filter_map(|f| Some((&*f.name, f.ns.as_deref(), f.span?)))
    }

    /// The hull of every motion stated in the model, including its imports.
    pub fn spans_hull(&self) -> Option<(f64, f64)> {
        self.files.iter().filter_map(|f| f.span).chain(self.t_span).reduce(|(a, b), (c, d)| (a.min(c), b.max(d)))
    }

    pub fn is_own(&self, a: &AssertDef) -> bool {
        a.scope.file == 1
    }

    /// Whether `file` (where something is written, as `AssertDef::file`)
    /// is the file loaded, not one it imports.
    pub fn is_own_file(&self, file: &str) -> bool {
        file == &*self.file
    }

    /// The file a param, def or fn named `name` is written in. None: a
    /// name of std's, or none.
    pub fn written_in(&self, name: &str) -> Option<&str> {
        self.origin.get(name).map(|&(f, _)| &*self.files[f as usize].name)
    }

    /// The files the file loaded imports, directly or through others, in
    /// the order they were read; std is not one.
    pub fn imported_files(&self) -> Vec<&str> {
        let mut out: Vec<&str> = Vec::new();
        for f in self.files.iter().skip(2) {
            if !out.contains(&&*f.name) {
                out.push(&f.name);
            }
        }
        out
    }

    /// The authored unit to use when emitting a new override, including
    /// optimizer replay commands. The numeric value itself is canonical.
    pub fn param_unit(&self, name: &str) -> Option<&'static units::Unit> {
        let k = self.param_index(name)?;
        self.param_said.get(&k).map(|s| s.unit).or_else(|| self.param_requirements.get(&self.params[k].name)?.bound_unit)
    }

    pub fn param_index(&self, name: &str) -> Option<usize> {
        match self.scoped(name)?.value {
            BoundName::Param(i) => Some(i),
            _ => None,
        }
    }

    pub fn body(&self, name: &str) -> Option<&Body> {
        self.bodies.iter().find(|b| b.name == name)
    }

    pub fn assembly(&self, name: &str) -> Option<&Assembly> {
        self.assemblies.iter().find(|a| a.name == name)
    }

    pub fn def_index_of(&self, name: &str) -> Option<usize> {
        self.scoped(name)?.def()
    }

    pub fn std_def_index(&self, name: &str) -> Option<usize> {
        self.files[0].bindings.get(name)?.def()
    }

    pub fn def(&self, name: &str) -> Option<&Value> {
        self.def_index_of(name).map(|i| &self.defs[i].1)
    }

    /// Whether param k is ranged where claims are answered: a tolerance,
    /// or a (ranged lo hi).
    pub fn has_range(&self, k: usize) -> bool {
        self.params[k].tol.is_some() || self.ranged.iter().any(|r| r.0 == k)
    }

    /// What to do where an answer needs one value of each of the ranged
    /// params `ks` it reads, over `spans` (the view's): a tolerance is
    /// checked at one size of it, each named with the ends of its band and
    /// the flags that pin it there (which end is worst is not proven, so
    /// each answer is that size's); a (ranged lo hi) is a range of the
    /// question, which no claim takes at one value, so its answer is
    /// looked at at one (eval, render, a viewer, ink check's conditions).
    pub fn at_one_value(&self, ks: &[usize], spans: &ink_core::eval::ParamSpans) -> String {
        let ranged: Vec<&str> = ks.iter().filter(|&&k| self.ranged.iter().any(|r| r.0 == k)).map(|&k| self.params[k].name.as_str()).collect();
        if !ranged.is_empty() {
            return format!(
                "'{}' is (ranged lo hi), a range of the question, so no claim takes it at one value. See one with ink eval or ink render --param NAME=V",
                ranged.join("', '")
            );
        }
        // An end as it is written: the band rounded outward is within a
        // float or two of it.
        let num = |x: f64| format!("{}", (x * 1e9).round() / 1e9);
        let each: Vec<String> = ks
            .iter()
            .map(|&k| {
                let name = &self.params[k].name;
                match spans.get(k as u32) {
                    Some((a, b)) if a < b => format!("--param {name}={a}:{a} or --param {name}={b}:{b} (its ends)", a = num(a), b = num(b)),
                    _ => format!("--param {name}=V:V"),
                }
            })
            .collect();
        match each.is_empty() {
            true => "check at fixed params with --nominal or --param NAME=V:V".into(),
            false => format!(
                "check it at one size of each: {}. Or use --nominal: each param at its value, no tolerance. Which size is worst is not proven: each answer is that size's",
                each.join("; ")
            ),
        }
    }

    /// The (ranged lo hi) named `name`, what --param NAME reaches: the def
    /// that made it, or its generated name ("t-slot-reaches 0"). Its names
    /// are unique, and none is a declared param's (a def's name is not; a
    /// generated one has a space), so load_params finds it by name too.
    pub fn ranged_named(&self, name: &str) -> Option<(usize, (f64, f64))> {
        self.ranged.iter().find(|r| self.params[r.0].name == name).copied()
    }

    /// The (ranged lo hi) a value is, if it is one itself.
    pub(super) fn ranged_of(&self, v: &Value) -> Option<(usize, (f64, f64))> {
        let &Value::Field(f, _) = v else { return None };
        let graph::Node::Var(Var::Param(k)) = self.g.nodes[f as usize] else { return None };
        self.ranged.iter().find(|r| r.0 == k as usize).copied()
    }

    pub(super) fn value_names(&self) -> Vec<String> {
        let mut v: Vec<String> = self.written_params().map(|(_, p)| p.name.clone()).collect();
        v.extend(self.defs.iter().map(|d| d.0.clone()));
        v.extend(self.fns.keys().cloned());
        v
    }

    /// A name a param is set by on the command line, `--param name=v`: a
    /// param's, written or generated (a `for`'s, a solve's unknown), or an
    /// import's namespace, which prefixes its params' names with `.`. It is
    /// one plain identifier (`param_word`), so every param's name, `ns.p`
    /// included, is set by --param as it is printed, unquoted, in any shell.
    pub(crate) fn check_param_name(&self, at: &Sexp, name: &str, what: &str) -> R<()> {
        if crate::command::param_word(name) {
            return Ok(());
        }
        Err(self.err_hint(
            at,
            format!(
                "{what} {name:?} is not a plain name. A param is set on the command line (--param name=v), so its name is letters, digits, '-' and '_', starting with a letter"
            ),
            Some("rename it, e.g. with '-' for a space".into()),
        ))
    }

    pub(super) fn check_new_value_name(&mut self, at: &Sexp, name: &str) -> R<()> {
        // A file's own definition may take a form's or an operation's name
        // (Doc::shadows); the constants (x y z t pi inf) are read, not
        // called, so a file's own one would change what they mean unseen.
        if CONSTANTS.contains(&name) || (self.scope.file == 0 && is_reserved(name)) {
            return Err(self.err_hint(at, format!("'{name}' is a built-in name and cannot be redefined"), Some(reserved_hint())));
        }
        let already = self.files[self.scope.file as usize]
            .bindings
            .get(name)
            .filter(|b| b.file == self.scope.file)
            .map(|b| (b.file, b.line))
            .or_else(|| self.origin.get(&self.public(name)).copied().filter(|_| !Self::std_name(name)))
            .map(|(f, line)| (format!("{}:{line}", self.files[f as usize].name), f));
        if let Some((at_, f)) = already {
            // An imported param: given at its import.
            let hint = if f != self.scope.file && self.params.iter().any(|p| p.name == name) {
                format!("an importer gives an imported param its value at the import, (import \"{}\" ({name} V))", self.files[f as usize].name)
            } else {
                "names are immutable; pick a new name".into()
            };
            return Err(self.err_hint(at, format!("'{name}' is already defined, at {at_}"), Some(hint)));
        }
        if self.scope.file != 0 && is_reserved(name) {
            let file = &mut self.files[self.scope.file as usize];
            if let Some(&line) = file.built_in_read.get(name) {
                return Err(self.err_hint(
                    at,
                    format!("this file read '{name}' as the built-in {name} at line {line}, before its own '{name}' here"),
                    Some(format!("one name has one meaning in a file: rename this one (my-{name}, say)")),
                ));
            }
            file.own.insert(name.to_string(), (at.line, at.col));
        }
        Ok(())
    }

    /// Expression `e` as a number known at load time, with what says its
    /// quantity (units::Said); `none`: the error where it is no number.
    /// A constant number `e`, what it is said in, and its exact fraction where it has one (Graph::real_value).
    pub(super) fn constant(&mut self, e: &Sexp, none: &str) -> R<(f64, Option<Rc<units::Said>>, Option<ink_core::exact::Ratio>)> {
        let v = self.eval(e, &mut Vec::new())?;
        let said = if let Value::Field(_, s) = &v { s.clone() } else { None };
        let f = self.as_field(v, e)?;
        let x = self.g.as_const(f).ok_or_else(|| self.err(e, none))?;
        Ok((x, said, self.g.real_value(f).and_then(|r| r.ratio)))
    }

    /// A param's `(tol lo hi)`: deviations from a value said by `said`;
    /// with the units its ends are said in, each where.
    pub(super) fn tol(&mut self, name: &str, t: &Sexp, said: Option<&units::Said>, rule: Option<&units::UnitRequirement>, usage: &str) -> R<((f64, f64), Vec<TolEnd>)> {
        let Some([_, lo, hi]) = t.list().filter(|l| l.first().and_then(|h| h.sym()) == Some("tol")) else {
            return Err(self.err_hint(t, "expected (tol lo hi)", Some(usage.into())));
        };
        let mut units = Vec::new();
        let mut num = |d: &mut Doc, e: &Sexp| -> R<f64> {
            let mut v = d.eval(e, &mut Vec::new())?;
            d.fits_change(&v, e, "a tolerance's deviation")?;
            // A bare end of a param that asks its inputs for units is refused below, by the param's name.
            if let Some(s) = said
                && (rule.is_none() || !crate::said(&v).is_empty())
            {
                d.fits(&mut v, e, s.unit.quantity, "a tolerance's deviation")?;
            }
            let said = if let Value::Field(_, s) = &v { s.clone() } else { None };
            let f = d.as_field(v, e)?;
            let value = d.g.as_const(f).filter(|x| x.is_finite()).ok_or_else(|| d.err(e, "a tolerance's ends must be constant finite numbers"))?;
            let end = TolEnd { value, said, line: e.line, col: e.col };
            if let Some(rule) = rule {
                end.check(rule, name).map_err(|m| d.err(e, m))?;
            }
            units.push(end);
            Ok(value)
        };
        let (lo, hi) = (num(self, lo)?, num(self, hi)?);
        if lo > hi {
            return Err(self.err_hint(t, format!("tol needs lo <= hi, got {lo} > {hi}"), Some(usage.into())));
        }
        Ok(((lo, hi), units))
    }

    /// An import's values, `(name value [(tol lo hi)])` each, read in the
    /// importer's names, with the clause each is written in.
    pub(super) fn import_values(&mut self, values: &[Sexp], usage: &str) -> R<Vec<(ImportValue, Sexp)>> {
        let mut given: Vec<(ImportValue, Sexp)> = Vec::new();
        for c in values {
            let Some([name, value, tol @ ..]) = c.list().filter(|l| l.len() <= 3) else {
                return Err(self.err_hint(c, "expected (name value [(tol lo hi)]): a value for a param of the imported file", Some(usage.into())));
            };
            let n = self.name_of(name)?;
            if given.iter().any(|(g, _)| g.name == n) {
                return Err(self.err(name, format!("'{n}' is given twice")));
            }
            // A number known at load (a param's value, baked as an if's
            // is), or a text.
            let v = self.eval(value, &mut Vec::new())?;
            if let Value::Str(t) = v {
                if !tol.is_empty() {
                    return Err(self.err(c, format!("'{n}' is given text, which takes no tolerance")));
                }
                given.push((
                    ImportValue {
                        name: n,
                        value: f64::NAN,
                        text: Some(t.into()),
                        said: None,
                        tol: None,
                        file: self.file.clone(),
                        line: value.line,
                        col: value.col,
                        ends: Vec::new(),
                        requirement: None,
                    },
                    c.clone(),
                ));
                continue;
            }
            let said = if let Value::Field(_, s) = &v { s.clone() } else { None };
            let f = self.as_field(v, value)?;
            let x = self.const_of(value, f)?;
            let requirement = units::UnitRequirement::of(value, x, said.as_deref());
            let (tol, ends) = match tol.first() {
                Some(t) => self.tol(&n, t, said.as_deref(), requirement.as_ref(), usage).map(|(b, u)| (Some(b), u))?,
                None => (None, Vec::new()),
            };
            given.push((ImportValue { name: n, value: x, text: None, said, tol, file: self.file.clone(), line: value.line, col: value.col, ends, requirement }, c.clone()));
        }
        Ok(given)
    }

    /// Import value `g` (clause `at`) sets a param written in file `i`
    /// (imported as `path`): refused where it names no `(param ...)` of
    /// that file, or, where the file was already run (first imported at
    /// `first`), where that run gave the param another value.
    pub(super) fn import_sets(&self, i: usize, g: &ImportValue, at: &Sexp, path: &str, first: Option<&str>) -> R<()> {
        let (f, n) = (&self.files[i], &g.name);
        let Some(binding) = f.bindings.get(n) else {
            let params = f.bindings.iter().filter(|(_, b)| b.file as usize == i && matches!(b.value, BoundName::Param(_) | BoundName::Text { .. })).map(|(n, _)| n.as_str());
            return Err(self.err_hint(at, format!("{} has no param '{n}'", f.name), closest(n, params)));
        };
        let not = |what: &str| Err(self.err(at, format!("'{n}' is {what}: an import sets params")));
        // A text param: set by the run that read it (TextParam), the value
        // asked the one that run read.
        if binding.file as usize != i && matches!(binding.value, BoundName::Param(_) | BoundName::Text { .. }) {
            let from = &self.files[binding.file as usize].name;
            return Err(self.err(at, format!("{n} is {from}'s: import {from} with it, before {}", f.name)));
        }
        if let BoundName::Text { param: t, .. } = binding.value {
            let p = &self.texts[t];
            // CLI overrides and text cases may select a different effective
            // value. Compare the defaults requested at the two import sites.
            let initial = f.given.iter().find(|v| v.name == *n).and_then(|v| v.text.as_deref()).unwrap_or(&p.default);
            return match first {
                Some(first) if g.text.as_deref() != Some(initial) => Err(self.err_hint(
                    at,
                    format!("{} was read at {first} with {n} {}", f.name, serde_json::json!(initial)),
                    Some(format!("give the value at {}'s first import ({first}), or import a second copy under a namespace", f.name)),
                )),
                _ => Ok(()),
            };
        }
        let k = match binding.value {
            BoundName::Param(k) => k,
            BoundName::Function(_) => return not("a fn, not a param"),
            _ if self.ranged_named(&f.ns.as_ref().map_or(n.clone(), |ns| format!("{ns}.{n}"))).is_some() => return not("a (ranged ...), a range of the question"),
            _ => return not("a def, not a param"),
        };
        let p = &self.params[k];
        let Some(first) = first else {
            // Not set by its param form: a package form's (a solve's unknown).
            return if self.set_at.contains_key(&k) { Ok(()) } else { not("made by a package form, not by (param ...)") };
        };
        if let Some(s) = self.param_said.get(&k) {
            g.fits(s.unit.quantity, &format!("{}'s {n}", f.name))?;
        }
        let asked = &at.list().expect("a clause")[1];
        let copy = format!("import a second copy under a namespace, (import \"{path}\" as NS ({n} {asked}))");
        let unit = self.param_said.get(&k).map_or("", |s| s.unit.quantity.unit());
        match self.set_at.get(&k) {
            Some(_) if p.default == g.value && g.tol.is_none_or(|t| p.tol == Some(t)) => {
                // Compare what the two import sites wrote. A command-line
                // override can give their shared, bare default a unit.
                if let Some(rule) = self.param_requirements.get(&p.name) {
                    g.check_units(rule, n)?;
                }
                Ok(())
            }
            Some((_, by)) => Err(self.err_hint(
                at,
                format!("{}'s {n} is already {}{unit}, set by {by}; this import asks {asked}", f.name, p.default),
                Some(format!("one load gives a file's param one value: {copy}")),
            )),
            None => Err(self.err_hint(
                at,
                format!("{} was read at {first} with {n} {}{unit}", f.name, p.default),
                Some(format!("give the value at {}'s first import ({first}), or {copy}", f.name)),
            )),
        }
    }

    fn std_name(name: &str) -> bool {
        static NAMES: std::sync::OnceLock<HashSet<String>> = std::sync::OnceLock::new();
        NAMES
            .get_or_init(|| {
                std_forms()
                    .iter()
                    .filter_map(|f| match f.list() {
                        Some([h, n, ..]) if matches!(h.sym(), Some("fn" | "def")) => n.sym().map(str::to_string),
                        _ => None,
                    })
                    .collect()
            })
            .contains(name)
    }

    /// Imports copy these immutable references. A public name shadowed
    /// later cannot change what another file reads.
    pub(super) fn defined(&mut self, name: &str, line: u32, value: BoundName) {
        let file = self.scope.file;
        self.files[file as usize].bindings.insert(name.to_string(), Binding { file, line, value });
        if file != 0 {
            let public = self.public(name);
            self.origin.insert(public.clone(), (file, line));
            self.fns.remove(&public);
        }
    }

    /// A name a top-level form of the file being read defines, as the
    /// document names it: in a file imported under a namespace, `ns.name`
    /// (`p.plate`, `p.q.plate` through a namespace inside it); else as written.
    pub(super) fn public(&self, name: &str) -> String {
        match &self.files[self.scope.file as usize].ns {
            Some(ns) => format!("{ns}.{name}"),
            None => name.to_string(),
        }
    }

    /// A body or assembly, or a path from one (`arm/tip`), as the file
    /// being read writes it, as the document names it: its first segment
    /// read in the file's names (SourceFile::objects).
    pub(super) fn resolve(&self, path: &str) -> Result<String, String> {
        let (first, rest) = path.split_once('/').map_or((path, None), |(a, b)| (a, Some(b)));
        if let Some(n) = self.files[self.scope.file as usize].objects.get(first) {
            return Ok(rest.map_or(n.clone(), |r| format!("{n}/{r}")));
        }
        if self.body(first).is_some() || self.assembly(first).is_some() || self.ext.iter().any(|e| e.name == first) {
            return Err(format!("'{first}' is outside this file's scope: import its definition before using it"));
        }
        Ok(path.to_string())
    }

    /// A name (a body, an assembly, a path from one) written at `e`, as
    /// the document names it (Doc::resolve).
    pub(super) fn path_of(&self, e: &Sexp) -> R<String> {
        self.resolve(&self.name_of(e)?).map_err(|m| self.err(e, m))
    }

    pub(super) fn scoped(&self, name: &str) -> Option<&Binding> {
        self.files[self.scope.file as usize].bindings.get(name).or_else(|| self.files[0].bindings.get(name))
    }

    pub(super) fn function(&self, name: &str) -> Option<&Rc<FnDef>> {
        match &self.scoped(name)?.value {
            BoundName::Function(f) => Some(f),
            _ => None,
        }
    }

    /// Whether a built-in's name read here (`at`, where it is written in
    /// the file being read) means the file's own definition, which shadows
    /// the built-in from where it is written on, wherever the read is run
    /// from (a function's body is read where it is written). One name has
    /// one meaning in a file: the built-in read anywhere before the file's
    /// own definition of it (in a function or a `for`'s forms too), or a
    /// definition of it an import brings, is an error naming both. std and
    /// imported files keep the built-in: each reads its own names.
    pub(super) fn shadows(&mut self, name: &str, at: &Sexp) -> R<bool> {
        if self.scope.file == 0 || !is_reserved(name) || CONSTANTS.contains(&name) {
            return Ok(false);
        }
        let file = &self.files[self.scope.file as usize];
        let rename = |line| format!("rename the file's own at line {line} (my-{name}, say) to call the built-in {name} here");
        let before = |line| format!("'{name}' here is the built-in {name}, but this file defines its own '{name}' at line {line}");
        if let Some(&(line, col)) = file.own.get(name) {
            // Its own, read after it is written and once it is defined (not
            // in its own value).
            let defined = file.bindings.get(name).is_some_and(|b| b.file == self.scope.file);
            if (at.line, at.col) > (line, col) && defined {
                return Ok(true);
            }
            return Err(self.err_hint(at, before(line), Some(rename(line))));
        }
        if let Some(b) = file.bindings.get(name) {
            let (from, line) = (&self.files[b.file as usize].name, b.line);
            return Err(self.err_hint(
                at,
                format!("'{name}' here is both the built-in {name} and {from}'s def at line {line}, imported"),
                Some(format!("import it under a namespace, (import \"...\" as p), and call its p.{name}")),
            ));
        }
        if let Some(&line) = file.ahead.get(name) {
            return Err(self.err_hint(at, before(line), Some(rename(line))));
        }
        self.files[self.scope.file as usize].built_in_read.entry(name.to_string()).or_insert(at.line);
        Ok(false)
    }

    /// The hint where a file's own definition of a built-in's name, called
    /// at `at`, is no function of those arguments: the built-in was likely
    /// meant.
    pub(super) fn shadow_hint(&self, name: &str) -> String {
        let line = self.scoped(name).map_or(0, |b| b.line);
        format!("'{name}' here is this file's own, defined at line {line}, not the built-in {name}: rename it (my-{name}, say) to call the built-in")
    }

    /// A read of a name the file being read defines further on, not yet
    /// defined: an error, at a top-level form or in a function (a file's
    /// name means its own definition throughout, never std's before it).
    pub(super) fn later(&self, s: &str, at: &Sexp) -> Option<InkError> {
        let file = &self.files[self.scope.file as usize];
        let line = file.ahead.get(s).filter(|_| self.scope.file != 0 && self.scoped(s).is_none_or(|b| b.file != self.scope.file))?;
        Some(self.err_hint(at, format!("unknown name '{s}': it is defined at line {line}, after this"), Some("a file is read in order, so define a name before it is used".into())))
    }

    /// The directory of the file being read: the paths it names are
    /// relative to it (a function's, in its body).
    pub(super) fn dir(&self) -> &Path {
        &self.files[self.scope.file as usize].dir
    }

    /// Why a name defined outside std is not visible here (Doc::scoped).
    pub(super) fn hidden_hint(&self, name: &str) -> Option<String> {
        let &(file, line) = self.origin.get(name).filter(|_| self.scope.file != 0 && !self.files[self.scope.file as usize].bindings.contains_key(name))?;
        Some(format!(
            "'{name}' is defined at {}:{line}, in a file this one does not import. A file sees only its own names, its imports' and std's",
            self.files[file as usize].name
        ))
    }

    /// Where a form headed by a top-level form's name, written inside
    /// `inside`, belongs: (t-span 0 2) in an assembly.
    pub(super) fn top_form_hint(&self, c: &Sexp, inside: &str) -> Option<String> {
        let head = c.list()?.first()?.sym()?;
        self.registry
            .top_forms()
            .split(", ")
            .any(|f| f == head)
            .then(|| format!("({head} ...) is a top-level form. Write it at the file's top level, not inside {inside} (`ink spec --section \"Top-level forms\"`)"))
    }

    pub(super) fn check_new_object_name(&self, at: &Sexp, name: &str) -> R<()> {
        if name.contains('/') {
            return Err(self.err(at, format!("'{name}': body/assembly names cannot contain '/'")));
        }
        if self.body(name).is_some() || self.assembly(name).is_some() {
            return Err(self.err(at, format!("a body or assembly named '{name}' already exists")));
        }
        Ok(())
    }
}
