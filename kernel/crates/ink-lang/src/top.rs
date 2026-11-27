//! Top-level forms: each form of a file (`param`, `def`, `fn`, `body`,
//! `assembly`, `assert`, a package's), run as it is reached.

use super::*;

impl Doc {
    pub(super) fn top(&mut self, form: &Sexp, base: &Path) -> R<()> {
        let items = form.list().ok_or_else(|| self.err(form, format!("top-level forms must be lists: one of {}", self.registry.top_forms())))?;
        let head = items.first().and_then(|h| h.sym()).ok_or_else(|| self.err(form, "expected a form name"))?;
        let args = &items[1..];
        // Any form but a claim may define what shape-of reads.
        if head != "require" {
            self.shape_of.clear();
        }
        match head {
            "import" => {
                let usage = "usage: (import \"std\"), (import \"path.ink\") or (import \"path.ink\" as ns): with a namespace \
                             its names read ns.name; then (name value [(tol lo hi)]) ... sets the file's params";
                // A plain import inside a namespaced file joins its namespace,
                // so a part's family stays out of the importer's names.
                // `as`: the namespace its names read under here (`as`), and
                // under in the document (`ns`: this file's, then it).
                let (target, ns, as_, values) = match args {
                    [t, a, n, values @ ..] if a.sym() == Some("as") => {
                        let at = n;
                        let n = self.name_of(n)?;
                        self.check_param_name(at, &n, "the namespace")?;
                        (t, Some(self.public(&n)), Some(n), values)
                    }
                    [t, values @ ..] => (t, self.files[self.scope.file as usize].ns.clone(), None, values),
                    [] => return Err(self.err_hint(form, "bad import", Some(usage.into()))),
                };
                let SK::Str(p) = &target.k else { return Err(self.err(target, "import path must be a string")) };
                if p == "std" {
                    if args.len() > 1 {
                        return Err(self.err(target, "std is imported as it is, (import \"std\"): your names shadow its names"));
                    }
                    if !self.std_imported {
                        self.std_imported = true;
                        self.run_source_as(STD_SRC, "std", None, base, None, Vec::new())?;
                    }
                } else {
                    let given = self.import_values(values, usage)?;
                    let (path, absent) = import_path(base, p);
                    // Each path looked for and not there rests on its
                    // directory's index too: a member of its name there
                    // would be read instead (family::member_of).
                    self.files_read.extend(absent.into_iter().flat_map(|f| [index_read(&f.with_file_name(family::INDEX)), (f, None)]));
                    // A family's member (family::member_of): its family's
                    // file, run as it, named by its own path (canonical: its
                    // directory's real path, then its name).
                    let member = family::member_of(&path).map_err(|e| {
                        self.files_read.extend([(path.clone(), None), index_read(&path.with_file_name(family::INDEX))]);
                        self.err(target, e)
                    })?;
                    let canon = canonical(&path);
                    let source = member.as_ref().map_or(path.as_path(), |m| m.family.as_path());
                    if let Some(i) = self.running.iter().position(|&r| self.files[r as usize].path == canon) {
                        let chain: Vec<&str> = self.running[i..].iter().map(|&r| &*self.files[r as usize].name).chain([&*self.files[self.running[i] as usize].name]).collect();
                        return Err(self.err_hint(
                            target,
                            format!("import cycle: {}", chain.join(" -> ")),
                            Some("a file cannot import itself, directly or through the files it imports".into()),
                        ));
                    }
                    let key = (canon, ns.clone());
                    let dir = path.parent().unwrap_or(base);
                    let (i, first) = if let Some((i, first)) = self.imported.get(&key).cloned() {
                        let previous = &self.files[i].dir;
                        // A relative filename in the current directory has an empty parent.
                        let real_dir = |d: &Path| if d.as_os_str().is_empty() { Path::new(".") } else { d }.canonicalize();
                        if dir != previous && real_dir(dir).ok().zip(real_dir(previous).ok()).is_none_or(|(a, b)| a != b) {
                            return Err(self.err_hint(
                                target,
                                format!("{} was imported at {first} from a different directory", shown(&path)),
                                Some("import the second copy under its own namespace so its relative paths are read beside it".into()),
                            ));
                        }
                        self.read_again(&path);
                        (i, Some(first))
                    } else {
                        let src = std::fs::read_to_string(source).inspect_err(|_| self.files_read.push((source.to_path_buf(), None))).map_err(|e| {
                            let dirs = import_dirs();
                            let hint = std::iter::once(path.clone()).chain(dirs.iter().map(|d| d.join(p))).find_map(|f| sibling(&f, p)).or_else(|| {
                                let dirs: Vec<String> = dirs.iter().map(|d| d.display().to_string()).collect();
                                let dirs = if dirs.is_empty() { "none: INK_PATH is unset".into() } else { dirs.join(", ") };
                                let from = if self.given_import_dirs().is_empty() { "INK_PATH's directories" } else { "the directories given its load, then INK_PATH's" };
                                Path::new(p).is_relative().then(|| format!("an import is read beside its file, else from {from} ({dirs})"))
                            });
                            if member.is_none() {
                                self.files_read.push(index_read(&path.with_file_name(family::INDEX)));
                            }
                            self.err_hint(target, not_read(&path, &e), hint)
                        })?;
                        self.import_depth += 1;
                        self.member = member.clone();
                        let i = self.run_source_as(&src, &shown(&path), Some(&path), dir, ns.as_deref(), given.iter().map(|(v, _)| v.clone()).collect())?;
                        self.import_depth -= 1;
                        self.imported.insert(key, (i, format!("{}:{}", self.file, form.line)));
                        (i, None)
                    };
                    // Stopped at a case its text param refuses (Doc::refused):
                    // the file is read only so far, so no more is asked of it.
                    if !self.refused.is_empty() {
                        return Ok(());
                    }
                    // The file sees the names of what it imports (Scope),
                    // under its namespace here (`as`).
                    for (v, at) in &given {
                        self.import_sets(i, v, at, p, first.as_deref())?;
                    }
                    let named = |n: &String| as_.as_ref().map_or(n.clone(), |a| format!("{a}.{n}"));
                    let bindings: Vec<(String, Binding)> = self.files[i].bindings.iter().map(|(n, k)| (named(n), k.clone())).collect();
                    let objects: Vec<(String, String)> = self.files[i].objects.iter().map(|(n, o)| (named(n), o.clone())).collect();
                    let me = self.scope.file as usize;
                    for (name, key) in bindings {
                        let own = self.files[me].bindings.get(&name).is_some_and(|b| b.file == self.scope.file);
                        if !own {
                            self.files[me].bindings.insert(name, key);
                        }
                    }
                    for (name, o) in objects {
                        if self.files[me].objects.get(&name) != Some(&self.public(&name)) {
                            self.files[me].objects.insert(name, o);
                        }
                    }
                }
            }
            "param" => {
                let usage = "usage: (param name default-value [(tol lo hi)] [(cases v ...)]). tol: deviations from the value, lo <= hi. (param d 8 (tol 0 0.022)) is 8.000..8.022. \
                             cases: other values ink require checks every claim at. A text param: (param name \"Ada\" [(chars \"ABC...\")] [(len lo hi)] [(cases \"...\" ...)])";
                let (name, default, opts) = match args {
                    [name, default, opts @ ..] => (name, default, opts),
                    _ => return Err(self.err_hint(form, "bad param", Some(usage.into()))),
                };
                let n = self.name_of(name)?;
                self.check_param_name(name, &n, "the param")?;
                self.check_new_value_name(name, &n)?;
                let dv = self.eval(default, &mut Vec::new())?;
                if let Value::Str(text) = dv {
                    return self.text_param(form, name, n, text, opts, usage);
                }
                let (mut tol, mut cases) = (None, None);
                for o in opts {
                    match o.list().and_then(|l| l.first()).and_then(|h| h.sym()) {
                        Some("tol") if tol.is_none() => tol = Some(o),
                        Some("cases") if cases.is_none() => cases = Some(o),
                        _ => return Err(self.err_hint(o, "a number param takes (tol lo hi) and (cases v ...)", Some(usage.into()))),
                    }
                }
                let (w, n) = (n.clone(), self.public(&n));
                if let Some((_, g)) = self.overrides.iter().find(|(k, _)| *k == n)
                    && g.value.is_nan()
                    && g.case.is_none()
                {
                    let m = crate::command::not_a_number(&n, g.text.as_deref().unwrap_or(""));
                    return Err(self.err_hint(form, m, Some(format!("{n} is a number param: give it a number, --param {n}=6.2"))));
                }
                let (mut v, mut said, exact) = self.constant(default, "param default must be a constant number")?;
                let mut exacts: Vec<(u64, ink_core::exact::Ratio)> = exact.map(|r| (v.to_bits(), r)).into_iter().collect();
                let requirement = units::UnitRequirement::of(default, v, said.as_deref());
                let mut tol = tol.map(|t| self.tol(&n, t, said.as_deref(), requirement.as_ref(), usage).map(|(b, _)| b)).transpose()?;
                if let Some(rule) = requirement {
                    self.param_requirements.insert(n.clone(), rule);
                }
                let mut case_said = None;
                if let Some(c) = cases {
                    let selected = self.overrides.iter().find(|(k, _)| *k == n).and_then(|(_, g)| g.case);
                    let mut vs = Vec::new();
                    for e in &c.list().expect("a list")[1..] {
                        let mut x = self.eval(e, &mut Vec::new())?;
                        if let Some(s) = &said {
                            self.fits(&mut x, e, s.unit.quantity, "a case")?;
                        }
                        let unit = match &x {
                            Value::Field(_, s) => {
                                if selected == Some(vs.len()) {
                                    case_said = s.clone();
                                }
                                s.as_ref().map(|s| s.unit)
                            }
                            _ => None,
                        };
                        let f = self.as_field(x, e)?;
                        let x = self.g.as_const(f).filter(|x| x.is_finite()).ok_or_else(|| self.err(e, "a case is a constant finite number, as the default is"))?;
                        exacts.extend(self.g.real_value(f).and_then(|r| r.ratio).map(|r| (x.to_bits(), r)));
                        vs.push(Given { unit, ..Given::from(x) });
                    }
                    self.take_case(&n, &vs, form)?;
                    self.cases.push((n.clone(), vs));
                }
                // Its import's value (ImportValue), in the quantity written
                // here, else in its own.
                if let Some(g) = self.files[self.scope.file as usize].given.iter().find(|g| g.name == w).cloned() {
                    if g.text.is_some() {
                        return Err(InkError {
                            message: format!("{}'s {w} is a number param, but this import gives it text", self.file),
                            file: g.file.to_string(),
                            line: g.line,
                            col: 0,
                            hint: None,
                            called_from: Vec::new(),
                            open: None,
                        });
                    }
                    if let Some(s) = said.as_ref().or(g.said.as_ref()) {
                        g.fits(s.unit.quantity, &format!("{}'s {w}", self.file))?;
                    }
                    if let Some(rule) = self.param_requirements.get(&n) {
                        g.check_units(rule, &n)?;
                    }
                    if let Some(rule) = g.requirement.clone() {
                        self.param_requirements.entry(n.clone()).or_insert(rule);
                    }
                    self.set_at.insert(self.params.len(), (v, format!("{}:{}", g.file, g.line)));
                    (v, said, tol) = (g.value, g.said.or_else(|| said.map(|s| Rc::new(units::Said { explicit_when: Some(self.g.constant(0.0)), ..(*s).clone() }))), g.tol.or(tol));
                }
                self.defined(&w, name.line, BoundName::Param(self.params.len()));
                if let Some(rule) = self.param_requirements.get(&n)
                    && let Some(case) = &case_said
                    && let Some((_, g)) = self.overrides.iter().find(|(k, _)| k == &n)
                {
                    rule.check(&n, "every term of the selected case", case.explicit_when.is_none(), g.value == 0.0).map_err(|m| self.err(form, m))?;
                }
                let given = self.param_given(&n, v, said, form)?;
                if let Some(s) = case_said.or(given) {
                    self.param_said.insert(self.params.len(), s);
                }
                if !exacts.is_empty() {
                    self.param_exact.insert(self.params.len(), exacts);
                }
                let p = Param { name: n, default: v, tol };
                if let Some(w) = p.outside_band() {
                    let w = format!("{}:{}: {w}", self.file, form.line);
                    if !self.load_warnings.contains(&w) {
                        self.load_warnings.push(w);
                    }
                }
                self.add_param(p);
            }
            "members" => {
                let usage = "usage: (members NAME [ITS-NAME]): the file draws a family, a part for each member of its table in the families.json beside it, \
                             NAME the record of the member it is drawn as and ITS-NAME that member's name, a string (`ink spec --section Families`)";
                let (name, its_name) = match args {
                    [name] => (name, None),
                    [name, its_name] => (name, Some(its_name)),
                    _ => return Err(self.err_hint(form, "bad members", Some(usage.into()))),
                };
                let me = self.scope.file as usize;
                if self.files[me].members_read {
                    return Err(self.err_hint(form, "a file declares its members once", Some(usage.into())));
                }
                self.files[me].members_read = true;
                let n = self.name_of(name)?;
                self.check_new_value_name(name, &n)?;
                let (w, n) = (n.clone(), self.public(&n));
                // Its member: the one it is run as (its path's); else, run as
                // itself, its table's first.
                let (record, member) = match self.files[me].member.clone() {
                    // the load noted what it rests on (Doc::member_read)
                    Some(m) => (m.record, m.name),
                    None => {
                        let index = self.dir().join(family::INDEX);
                        let file = self.files[me].path.clone();
                        let stem = file.file_stem().and_then(|s| s.to_str()).unwrap_or_default().to_string();
                        let (fams, content) = match family::in_dir(self.dir()) {
                            Ok(Some(f)) => f,
                            Ok(None) => {
                                self.files_read.push((index.clone(), None));
                                return Err(self.err_hint(
                                    form,
                                    format!("no {} beside {}: a family's members are in its table there", family::INDEX, self.file),
                                    Some(usage.into()),
                                ));
                            }
                            Err(e) => {
                                self.files_read.push(index_read(&index));
                                return Err(self.err(form, e));
                            }
                        };
                        let Some((first, record)) = fams.family(&stem).and_then(|m| m.first()).cloned() else {
                            self.files_read.push((index.clone(), Some(content)));
                            return Err(self.err(form, format!("{} has no members of family {stem}", shown(&index))));
                        };
                        self.member_read(&family::Member { family: file, name: first.clone(), index, content, record: record.clone() });
                        (record, first)
                    }
                };
                let v = self.json_value(&Data::new(record)).map_err(|m| self.err(name, m))?;
                self.defined(&w, name.line, BoundName::Def(self.defs.len()));
                self.defs.push((n, v));
                // The member's own name, where the file asks for it: what it would be called as a file of its own.
                if let Some(its_name) = its_name {
                    let n = self.name_of(its_name)?;
                    self.check_new_value_name(its_name, &n)?;
                    let (w, n) = (n.clone(), self.public(&n));
                    self.defined(&w, its_name.line, BoundName::Def(self.defs.len()));
                    self.defs.push((n, Value::Str(member.as_str().into())));
                }
            }
            "def" => {
                if let Some(written) = fn_written_as_def(args) {
                    return Err(self.err_hint(form, "a function is defined by fn, not def", Some(format!("write {written}"))));
                }
                let [name, e] = args else { return Err(self.err(form, "usage: (def name expr)")) };
                let n = self.name_of(name)?;
                self.check_new_value_name(name, &n)?;
                let (w, n) = (n.clone(), self.public(&n));
                let made = self.ranged.len();
                let mut v = self.eval(e, &mut Vec::new())?;
                // A (ranged lo hi) the def makes and is takes its name, the
                // one --param reaches it by (Doc::ranged_named).
                if let Some((k, _)) = self.ranged_of(&v).filter(|r| self.ranged[made..].contains(r)) {
                    self.params[k].name = n.clone();
                    self.g.name_param(k as u32, &n);
                    if let Value::Field(f, said) = v {
                        let said = self.param_given(&n, self.params[k].default, said, e)?;
                        self.param_said.remove(&k);
                        self.param_said.extend(said.clone().map(|s| (k, s)));
                        v = Value::Field(f, said);
                    }
                }
                // Every def of a spatial field names a construction: its surface
                // points report this name (see eval::provenance), which under
                // a namespace differs by it. Constants stay
                // untagged so they keep folding.
                if let Value::Field(f, said) = &v
                    && self.g.depends_on(*f, |x| matches!(x, Var::X | Var::Y | Var::Z))
                {
                    v = Value::Field(self.g.tag(&n, *f), said.clone());
                }
                self.defined(&w, name.line, BoundName::Def(self.defs.len()));
                self.defs.push((n, v));
            }
            "fn" => {
                let [name, ps, body] = args else { return Err(self.err(form, "usage: (fn name (arg ...) body)")) };
                let n = self.name_of(name)?;
                self.check_new_value_name(name, &n)?;
                let (w, n) = (n.clone(), self.public(&n));
                let (params, rest) = self.arg_names(ps, "fn")?;
                let def = Rc::new(FnDef { params, rest, body: body.clone(), file: self.file.clone(), scope: self.scope });
                self.defined(&w, name.line, BoundName::Function(def.clone()));
                if !self.origin.contains_key(&n) || self.scope.file != 0 {
                    self.fns.insert(n, def);
                }
            }
            "body" => self.top_def(form, args, true)?,
            "assembly" => self.top_def(form, args, false)?,
            // The word before 0.2: the same form under `require`.
            "assert" => {
                let rest: Vec<String> = args.iter().map(|a| a.to_string()).collect();
                return Err(self.err_hint(form, "`assert` is now `require`", Some(format!("write (require {}): the form is the same, only the word changed", rest.join(" ")))));
            }
            "require" => {
                let usage = "usage: (require NAME (within-ratings PATH [(excluding P ...)])): asks whether every rating of every part at or under PATH (an assembly or an \
                             instance path) holds, less the parts P under it left out on purpose; (require NAME CLAIM [(accuracy x)] [(res r)]): asks whether CLAIM holds (Claims)";
                let [name, asked, options @ ..] = args else { return Err(self.err_hint(form, "wrong arguments", Some(usage.into()))) };
                let n = self.public(&self.name_of(name)?);
                if self.requires.iter().any(|r| r.name == n) {
                    return Err(self.err(name, format!("require '{n}' already exists")));
                }
                let mut excluding = Vec::new();
                let within = match asked.list() {
                    Some([head, path @ ..]) if head.sym() == Some("within-ratings") => {
                        let ([path, left @ ..], []) = (path, options) else {
                            return Err(self.err_hint(
                                asked,
                                "a requirement of ratings asks (within-ratings PATH [(excluding P ...)]), and takes no other clause",
                                Some(usage.into()),
                            ));
                        };
                        let within = self.path_of(path)?;
                        self.items(&within).map_err(|m| self.err(path, m))?;
                        for clause in left {
                            let Some([head, parts @ ..]) = clause.list().filter(|l| l.first().and_then(|h| h.sym()) == Some("excluding") && l.len() > 1) else {
                                return Err(self.err_hint(
                                    clause,
                                    "a requirement of ratings takes (excluding P ...), the parts under PATH left out on purpose",
                                    Some(usage.into()),
                                ));
                            };
                            let _ = head;
                            for part in parts {
                                let left = self.path_of(part)?;
                                self.items(&left).map_err(|m| self.err(part, m))?;
                                if !left.strip_prefix(within.as_str()).is_some_and(|rest| rest.starts_with('/')) {
                                    return Err(self.err(part, format!("{left} is not under {within}: (excluding ...) leaves out parts under the path asked")));
                                }
                                excluding.push(left);
                            }
                        }
                        Some(within)
                    }
                    // A claim: stated as an assert of this name states it, and answered by what answers that.
                    _ => {
                        self.state_claim(form, name, asked, options)?;
                        None
                    }
                };
                self.requires.push(Require { name: n, within, excluding, file: self.file.clone(), line: form.line });
            }
            "t-span" => {
                let usage = "usage: (t-span A B): claims that read t are checked at every moment from A to B";
                let [a, b] = args else { return Err(self.err_hint(form, "wrong arguments", Some(usage.into()))) };
                let end = |d: &mut Doc, e: &Sexp| {
                    let f = d.field_of(e, Quantity::Time, "t-span", &mut Vec::new())?;
                    d.const_of(e, f)
                };
                let (a, b) = (end(self, a)?, end(self, b)?);
                if !(a < b && a.is_finite() && b.is_finite()) {
                    return Err(self.err_hint(form, "a time span needs finite A < B", Some(usage.into())));
                }
                // The file's motion, which its own claims and assemblies are
                // checked over, in an importer too (span_of); the run file's
                // is the default of what it states.
                if self.import_depth == 0 {
                    self.t_span = Some((a, b));
                }
                self.files[self.scope.file as usize].span = Some((a, b));
            }
            other => {
                if let Some((package, head)) = self.registry.package(other) {
                    return self.top_ext(package, head, form, args);
                }
                let all = self.registry.top_forms();
                // A form of 0.1 that left with the constraint solver.
                if matches!(other, "solve" | "optimize") {
                    let hint = format!(
                        "({other} ...) left in 0.2: the constraint solver and optimiser return as a package. Until then give the value as a param and state what it must meet as a (require ...)"
                    );
                    return Err(self.err_hint(&items[0], format!("unknown top-level form '{other}'"), Some(hint)));
                }
                return Err(self.err_hint(&items[0], format!("unknown top-level form '{other}'"), closest(other, all.split(", ")).or(Some(format!("top-level forms: {all}")))));
            }
        }
        Ok(())
    }

    /// A load at a combination of declared cases: param `n`'s override,
    /// where it names one of its cases by place (Given::case), takes that
    /// case's value here, `vs` (its expressions' values, read with every
    /// param before it set).
    fn take_case(&mut self, n: &str, vs: &[Given], form: &Sexp) -> R<()> {
        if let Some((_, g)) = self.overrides.iter_mut().find(|(k, _)| k == n)
            && let Some(k) = g.case
        {
            let v = vs.get(k).ok_or_else(|| InkError::msg(format!("{n} has {} cases, not a case {}", vs.len(), k + 1)));
            let v = v.map_err(|e| InkError { file: self.file.to_string(), line: form.line, col: form.col, ..e })?;
            *g = Given { case: Some(k), ..v.clone() };
        }
        Ok(())
    }

    /// `(param name "Ada" [(chars S)] [(len lo hi)])`: a text param named
    /// `n` (written `name`), its default `text`. Its value (the default,
    /// an import's, `--param`'s), checked against what it allows, is bound
    /// as a def is.
    fn text_param(&mut self, form: &Sexp, name: &Sexp, n: String, text: String, opts: &[Sexp], usage: &str) -> R<()> {
        let (mut chars, mut len, mut cases) = (None, None, None::<Vec<(String, Sexp)>>);
        for o in opts {
            match o.list() {
                Some([h, v]) if h.sym() == Some("chars") => match self.eval(v, &mut Vec::new())? {
                    Value::Str(c) if chars.is_none() => chars = Some(c),
                    Value::Str(_) => return Err(self.err(o, "(chars S) is given twice")),
                    other => return Err(self.err_hint(v, format!("(chars S) takes a text, the characters allowed, not a {}", other.kind()), Some(usage.into()))),
                },
                Some([h, lo, hi]) if h.sym() == Some("len") && len.is_none() => {
                    let (lo, hi) = (self.const_num(lo, &mut Vec::new())?, self.const_num(hi, &mut Vec::new())?);
                    let whole = |x: f64| x >= 0.0 && x.fract() == 0.0 && x < 1e9;
                    if !(whole(lo) && whole(hi) && lo <= hi) {
                        return Err(self.err_hint(o, format!("(len lo hi) takes whole numbers 0 <= lo <= hi, not {lo} {hi}"), Some(usage.into())));
                    }
                    len = Some((lo as usize, hi as usize));
                }
                Some([h, vs @ ..]) if h.sym() == Some("cases") && cases.is_none() => {
                    let mut c = Vec::new();
                    for v in vs {
                        match self.eval(v, &mut Vec::new())? {
                            Value::Str(t) => c.push((t, v.clone())),
                            other => return Err(self.err_hint(v, format!("a text param's case is a text, not a {}", other.kind()), Some(usage.into()))),
                        }
                    }
                    cases = Some(c);
                }
                Some([h, ..]) if h.sym() == Some("tol") => return Err(self.err_hint(o, "a text param has no tolerance", Some(usage.into()))),
                _ => return Err(self.err_hint(o, "a text param takes (chars S), (len lo hi) and (cases \"...\" ...)", Some(usage.into()))),
            }
        }
        let (w, n) = (n.clone(), self.public(&n));
        let mut p = TextParam { name: n.clone(), default: text.clone(), value: text, chars, len, set_at: None };
        // A load at a combination of declared cases (ink require's, every
        // override a case) records the param's refusal of the text it
        // takes, whoever gave it (its default too: another param's case
        // may set its len), and goes on: that combination is not asked.
        let at_case = self.overrides.iter().any(|(_, g)| g.case.is_some());
        if !at_case && let Some(m) = p.refused(&p.default, &format!("{}'s default", n)) {
            return Err(self.err(form, m));
        }
        // Whether a case is a text the param allows is decided at each
        // combination ink require asks it in, every param set (another's
        // case may set its len): `refused`.
        if let Some(c) = cases {
            let vs: Vec<Given> = c.iter().map(|(t, _)| Given::text(t)).collect();
            self.take_case(&n, &vs, form)?;
            self.cases.push((n.clone(), vs));
        }
        if let Some(g) = self.files[self.scope.file as usize].given.iter().find(|g| g.name == w).cloned() {
            let at = format!("{}:{}", g.file, g.line);
            let bad = |m: String| InkError { message: m, file: g.file.to_string(), line: g.line, col: 0, hint: None, called_from: Vec::new(), open: None };
            let Some(t) = g.text.as_deref() else { return Err(bad(format!("{}'s {w} is a text param: give it text, ({w} \"...\")", self.file))) };
            if !at_case && let Some(m) = p.refused(t, &format!("the import at {at}")) {
                return Err(bad(m));
            }
            (p.value, p.set_at) = (t.to_string(), Some(at));
        }
        if let Some((_, g)) = self.overrides.iter().find(|(k, _)| *k == n).cloned() {
            let t = g.text.as_deref().unwrap_or_default().to_string();
            match p.refuses(&t) {
                Some(r) if !at_case => return Err(self.err(form, format!("--param {n}: {}", r.said(&p, &t)))),
                _ => p.value = t,
            }
        }
        // Refused at a combination: recorded, and the load stops here,
        // before anything reads the text (run_forms).
        if at_case && let Some(r) = p.refuses(&p.value) {
            self.refused.push((p.clone(), p.value.clone(), r));
            return Ok(());
        }
        self.defined(&w, name.line, BoundName::Text { param: self.texts.len(), def: self.defs.len() });
        self.defs.push((n, Value::Str(p.value.clone())));
        self.texts.push(p);
        Ok(())
    }

    /// `(body HEAD clauses...)` or `(assembly HEAD clauses...)`: HEAD is
    /// NAME, or `(NAME (P DEFAULT) ...)` for one taking arguments. Built
    /// here with its defaults; an instance passing other values builds it
    /// again (`instantiate`).
    fn top_def(&mut self, form: &Sexp, args: &[Sexp], body: bool) -> R<()> {
        let usage = if body {
            "usage: (body name (shape expr) (attr key value) (frame name F) ...), or with arguments (body (name (p default) ...) ...)"
        } else {
            "usage: (assembly name (instance name child [WHERE]) (frame name F) ...), or with arguments (assembly (name (p default) ...) ...)"
        };
        let Some((head, clauses)) = args.split_first() else { return Err(self.err(form, usage)) };
        let (name, ps) = match head.list() {
            Some([name, ps @ ..]) => (name, ps),
            Some([]) => return Err(self.err(head, usage)),
            None => (head, &[][..]),
        };
        let w = self.name_of(name)?;
        let n = self.public(&w);
        self.check_new_object_name(name, &n)?;
        self.files[self.scope.file as usize].objects.insert(w, n.clone());
        let mut params: Vec<(String, Sexp)> = Vec::new();
        for p in ps {
            let Some([pn, d]) = p.list() else { return Err(self.err_hint(p, "an argument is (name default)", Some(usage.into()))) };
            let pn = self.name_of(pn)?;
            if params.iter().any(|(q, _)| *q == pn) {
                return Err(self.err(p, format!("duplicate argument '{pn}'")));
            }
            params.push((pn, d.clone()));
        }
        let def = Rc::new(ArgDef { params, clauses: clauses.to_vec(), file: self.file.clone(), scope: self.scope });
        let (env, texts) = self.def_args(&def, &n, Vec::new(), form)?;
        let key = (n.clone(), format!("{env:?}"), (!body).then_some(self.scope));
        self.build(n.clone(), None, &def, env, texts, None, body, form)?;
        if !def.params.is_empty() {
            self.built.insert(key, n.clone());
            self.arg_defs.insert(n, def);
        }
        Ok(())
    }

    /// The values a definition's arguments take: those passed (as values,
    /// and as written), then the defaults, each read in the definition's
    /// file with the arguments before it.
    pub(super) fn def_args(&mut self, def: &ArgDef, name: &str, args: Vec<(Value, String)>, at: &Sexp) -> R<(Env, Vec<(String, String)>)> {
        if args.len() > def.params.len() {
            let ps: Vec<&str> = def.params.iter().map(|p| p.0.as_str()).collect();
            return Err(self.err(at, format!("'{name}' takes {} argument(s) ({}), got {}", ps.len(), ps.join(" "), args.len())));
        }
        let mut env = Env::new();
        let mut texts = Vec::new();
        let mut args = args.into_iter();
        let local = self.local_assembly.take();
        let saved = (std::mem::replace(&mut self.file, def.file.clone()), std::mem::replace(&mut self.scope, def.scope));
        let r = (|| {
            for (p, d) in &def.params {
                let (v, text) = match args.next() {
                    Some(a) => a,
                    None => (self.eval(d, &mut env)?, d.to_string()),
                };
                env.push((p.clone(), v));
                texts.push((p.clone(), text));
            }
            Ok(())
        })();
        self.local_assembly = local;
        (self.file, self.scope) = saved;
        r.map(|_| (env, texts))
    }
}
