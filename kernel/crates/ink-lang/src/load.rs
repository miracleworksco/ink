//! Loading a document: finding and reading its files and imports, running
//! their forms in order, and errors at the form being run.

use super::*;

/// The names a file's forms define at top level, as written: values,
/// functions, bodies, assemblies, claims and package forms.
fn top_names(forms: &[Sexp], registry: &Registry) -> HashSet<String> {
    // Core forms that define a name; every extension form does.
    const DEFINING: &[&str] = &["def", "fn", "param", "members", "body", "assembly", "require"];
    let mut all = HashSet::new();
    for f in forms {
        let Some([h, n, ..]) = f.list() else { continue };
        // A body or an assembly with arguments is named by its list's head.
        let n = match n.list() {
            Some([n, ..]) if matches!(h.sym(), Some("assembly" | "body")) => n,
            _ => n,
        };
        let (Some(h), Some(n)) = (h.sym(), n.sym()) else { continue };
        if DEFINING.contains(&h) || registry.package(h).is_some() {
            all.insert(n.to_string());
        }
    }
    all
}

thread_local! {
    /// The directories the load running on this thread was given
    /// (`Doc::load_file_importing`), read before `$INK_PATH`'s.
    static GIVEN: std::cell::RefCell<Vec<std::path::PathBuf>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// The directories an import not found beside its file is read from, in
/// order: those the load running on this thread was given
/// (`Doc::load_file_importing`), then `$INK_PATH`'s, separated as `$PATH`'s
/// are. A host or a distribution gives its libraries so (a parts catalog:
/// `parts/` in one).
pub fn import_dirs() -> Vec<std::path::PathBuf> {
    let dirs = std::env::var_os("INK_PATH").unwrap_or_default();
    let ink_path = std::env::split_paths(&dirs).filter(|d| !d.as_os_str().is_empty()).map(|d| d.canonicalize().unwrap_or(d));
    GIVEN.with_borrow(|given| given.iter().cloned().chain(ink_path).collect())
}

/// Run `f` with `dirs` as the directories given its load (`import_dirs`),
/// the ones given before back once it returns, however it returns: a load
/// on another thread, or after this one, sees none of them.
fn importing<T>(dirs: &[std::path::PathBuf], f: impl FnOnce() -> T) -> T {
    struct Restore(Vec<std::path::PathBuf>);
    impl Drop for Restore {
        fn drop(&mut self) {
            GIVEN.with_borrow_mut(|g| *g = std::mem::take(&mut self.0));
        }
    }
    let dirs = dirs.iter().filter(|d| !d.as_os_str().is_empty()).map(|d| d.canonicalize().unwrap_or_else(|_| d.clone())).collect();
    let _restore = Restore(GIVEN.with_borrow_mut(|g| std::mem::replace(g, dirs)));
    f()
}

/// A path as a load keys the file it names: its real path; for one that is
/// no file (a family's member: family::member_of), its directory's real
/// path, then its name.
pub(super) fn canonical(path: &Path) -> std::path::PathBuf {
    path.canonicalize().unwrap_or_else(|_| {
        let dir = path.parent().map(|d| if d.as_os_str().is_empty() { Path::new(".") } else { d });
        dir.and_then(|d| d.canonicalize().ok()).zip(path.file_name()).map_or_else(|| path.to_path_buf(), |(d, n)| d.join(n))
    })
}

/// The file `(import "p")` in a file in `base` reads (a family's member
/// counts as one: family::member_of): `p` relative to
/// `base` (or absolute); else the first of `import_dirs` holding it (none:
/// `p` relative to `base`). And the files looked for before it, absent: a
/// load rests on their absence as on what it reads.
pub fn import_path(base: &Path, p: &str) -> (std::path::PathBuf, Vec<std::path::PathBuf>) {
    // An index it cannot read is there: reading the path says why.
    let there = |f: &Path| f.exists() || !matches!(family::member_of(f), Ok(None));
    let beside = base.join(p);
    if there(&beside) || Path::new(p).is_absolute() {
        return (beside, vec![]);
    }
    let mut tried = vec![beside.clone()];
    for f in import_dirs().into_iter().map(|d| d.join(p)) {
        if there(&f) {
            return (f, tried);
        }
        tried.push(f);
    }
    tried.remove(0);
    (beside, tried)
}

/// A file as ink names it to the user (in errors, claims' sources, `ink
/// check`'s imports): one in an import directory (`import_dirs`) by its
/// path there, the path that imports it from anywhere (an agent of a
/// hosted ink opens it so); any other as given.
pub fn shown(path: &Path) -> String {
    let real = path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
    match import_dirs().iter().find_map(|d| real.strip_prefix(d).ok().filter(|r| !r.as_os_str().is_empty())) {
        // An import name goes back into a model, where backslashes start escapes.
        Some(r) => r.to_string_lossy().replace(std::path::MAIN_SEPARATOR, "/"),
        None => path.display().to_string(),
    }
}

/// Why `path` could not be read (`e`): where its directory has families
/// (family::INDEX), a member's name read there as none.
pub(super) fn not_read(path: &Path, e: &std::io::Error) -> String {
    let index = path.with_file_name(family::INDEX);
    match family::in_dir(path.parent().unwrap_or(Path::new(""))) {
        Ok(Some(_)) if e.kind() == std::io::ErrorKind::NotFound && path.extension().is_some_and(|x| x == "ink") => {
            format!("cannot read {}: no such file, and {} names no member {}", shown(path), shown(&index), path.file_stem().unwrap_or_default().to_string_lossy())
        }
        _ => format!("cannot read {}: {e}", shown(path)),
    }
}

/// A directory's index (family::INDEX) as a failed load rests on it:
/// its content's key, or absent.
pub(super) fn index_read(index: &Path) -> (std::path::PathBuf, Option<u128>) {
    (index.to_path_buf(), std::fs::read(index).ok().map(|b| Doc::content_key(&b)))
}

/// For an import `written` of a file `path` that cannot be read: the file
/// (or a family's member) in its directory a typo (or an old extension)
/// could make of its name.
pub(super) fn sibling(path: &Path, written: &str) -> Option<String> {
    let name = path.file_name()?.to_str()?;
    let mut names: Vec<String> = std::fs::read_dir(path.parent()?).ok()?.filter_map(|e| e.ok()?.file_name().into_string().ok()).collect();
    // A family's members are named as files are (family::member_of).
    if let Ok(Some((fams, _))) = family::in_dir(path.parent()?) {
        names.extend(fams.names().map(|n| format!("{n}.ink")));
    }
    let found = near(name, names.iter().map(String::as_str)).into_iter().next()?;
    let dir = written.rsplit_once('/').map_or("", |(d, _)| d);
    Some(format!("did you mean \"{dir}{}{found}\"?", if dir.is_empty() { "" } else { "/" }))
}

impl Doc {
    fn empty(registry: &'static Registry) -> Doc {
        Doc {
            solved: HashMap::new(),
            registry,
            g: Graph::new(),
            params: Vec::new(),
            texts: Vec::new(),
            refused: Vec::new(),
            fonts: Default::default(),
            cases: Vec::new(),
            defs: Vec::new(),
            fns: HashMap::new(),
            bodies: Vec::new(),
            assemblies: Vec::new(),
            arg_defs: HashMap::new(),
            built: HashMap::new(),
            args: Env::new(),
            asserts: Vec::new(),
            requires: Vec::new(),
            claim_declarations: HashMap::new(),
            ext: Vec::new(),
            solves: Vec::new(),
            published: Vec::new(),
            wanted: HashSet::new(),
            provided: Default::default(),
            follow: Default::default(),
            t_span: None,
            placed: HashMap::new(),
            placed_labels: HashMap::new(),
            imported: HashMap::new(),
            set_at: HashMap::new(),
            parsed: HashMap::new(),
            std_imported: false,
            json: HashMap::new(),
            sources_read: Vec::new(),
            files_read: Vec::new(),
            load_warnings: Vec::new(),
            dropped_limits: HashMap::new(),
            param_said: HashMap::new(),
            param_requirements: HashMap::new(),
            checks: Vec::new(),
            checks_seen: HashSet::new(),
            calls: Vec::new(),
            solved_params: HashSet::new(),
            found_from: Vec::new(),
            ranged: Vec::new(),
            allowances: Vec::new(),
            support_depths: HashSet::new(),
            print_bounds: HashMap::new(),
            value_bounds: HashMap::new(),
            param_exact: HashMap::new(),
            held: HashMap::new(),
            rated: Vec::new(),
            unsolved: None,
            unsolved_t: Vec::new(),
            conditions: Vec::new(),
            magnitudes: Vec::new(),
            claim_materials: false,
            bound: None,
            placed_fields: HashMap::new(),
            warped: Vec::new(),
            warps: HashMap::new(),
            xyz_read: HashMap::new(),
            integrand_of: None,
            file: Rc::from(""),
            scope: Scope::STD,
            files: vec![SourceFile::new(Rc::from("std"), "std".into(), None, ".".into())],
            running: Vec::new(),
            origin: HashMap::new(),
            depth: 0,
            inputs: Vec::new(),
            import_depth: 0,
            overrides: Vec::new(),
            baked: Default::default(),
            materialized: Default::default(),
            covered: Vec::new(),
            load_t: Vec::new(),
            undefined_at_load: Default::default(),
            shape_of: HashMap::new(),
            form_at: String::new(),
            wild_written: None,
            form_name: String::new(),
            building: None,
            local_assembly: None,
            own_frame: false,
            source: None,
            coordinates_of: None,
            states: false,
            state_t: None,
            coords: Vec::new(),
            coordinate_params: Vec::new(),
            member: None,
        }
    }

    /// Load a file with the packages of `registry`, and parameter
    /// overrides applied to load-time computations (parameters stay
    /// symbolic in fields, so derivatives still work). The standard
    /// distribution's loaders (`ink::load_file`) pass `ink::registry()`.
    pub fn load_file_in(registry: &'static Registry, path: &Path, overrides: &[(String, Given)]) -> R<Doc> {
        Self::load_file_reads(registry, path, overrides).map_err(|(e, _)| e)
    }

    /// `load_file_in`, and for a load that fails, what it read before
    /// failing (`files_read`, the file it could not read among them): a
    /// failed load rests on those as a load rests on its `files_read`, so
    /// it fails the same way while each reads as it did.
    pub fn load_file_reads(registry: &'static Registry, path: &Path, overrides: &[(String, Given)]) -> Result<Doc, (InkError, Vec<(std::path::PathBuf, Option<u128>)>)> {
        Self::read_file(registry, path, overrides, false, &[])
    }

    /// `load_file_in`, an import not found beside its file read from
    /// `dirs`, in order, before `$INK_PATH`'s (`import_dirs`): a library
    /// caller's catalog for this file, without setting `INK_PATH` for the
    /// whole process. The imports it imports, transitively, are found so
    /// too, and so is the file loaded again (Doc::with_coordinates,
    /// ModelView::reload). No other load sees them, on this thread or
    /// another.
    pub fn load_file_importing(registry: &'static Registry, path: &Path, overrides: &[(String, Given)], dirs: &[std::path::PathBuf]) -> R<Doc> {
        Self::read_file(registry, path, overrides, false, dirs).map_err(|(e, _)| e)
    }

    /// The directories the load was given (`load_file_importing`), read
    /// before `$INK_PATH`'s.
    pub fn given_import_dirs(&self) -> &[std::path::PathBuf] {
        self.source.as_ref().map_or(&[], |o| &o.import_dirs)
    }

    pub(super) fn read_file(
        registry: &'static Registry,
        path: &Path,
        overrides: &[(String, Given)],
        claim_materials: bool,
        dirs: &[std::path::PathBuf],
    ) -> Result<Doc, (InkError, Vec<(std::path::PathBuf, Option<u128>)>)> {
        importing(dirs, || Self::read_file_here(registry, path, overrides, claim_materials, dirs))
    }

    fn read_file_here(
        registry: &'static Registry,
        path: &Path,
        overrides: &[(String, Given)],
        claim_materials: bool,
        dirs: &[std::path::PathBuf],
    ) -> Result<Doc, (InkError, Vec<(std::path::PathBuf, Option<u128>)>)> {
        let base = path.parent().unwrap_or(Path::new("."));
        // A member of a family (family::member_of): the family's file, run
        // as that member, named by its own path. The load rests on there
        // being no file of its path, and on its directory's index.
        let index = path.parent().unwrap_or(Path::new("")).join(family::INDEX);
        let member = family::member_of(path).map_err(|e| (InkError::msg(e), vec![(path.to_path_buf(), None), index_read(&index)]))?;
        if let Some(m) = member {
            let src = std::fs::read_to_string(&m.family).map_err(|e| (InkError::msg(format!("cannot read {}: {e}", shown(&m.family))), vec![(m.family.clone(), None)]))?;
            let origin = Origin { src, name: shown(path), path: path.to_path_buf(), base: base.to_path_buf(), member: Some(m), import_dirs: dirs.to_vec() };
            return Self::load_origin(registry, origin, overrides, None, claim_materials);
        }
        let src = std::fs::read_to_string(path).map_err(|e| {
            let error = InkError { hint: sibling(path, &path.display().to_string()), ..InkError::msg(not_read(path, &e)) };
            (error, vec![(path.to_path_buf(), None), index_read(&index)])
        })?;
        let origin = Origin { src, name: shown(path), path: path.to_path_buf(), base: base.to_path_buf(), member: None, import_dirs: dirs.to_vec() };
        Self::load_origin(registry, origin, overrides, None, claim_materials)
    }

    /// Note a file the load read (a source it ran, std included, or data
    /// it read): its name as shown (`shown`: in the catalog, `parts/...`),
    /// where it is (std: none, it is the language's) and its content.
    pub(super) fn read(&mut self, name: &str, path: Option<&Path>, bytes: &[u8]) {
        let mut fp = graph::Fingerprint::default();
        fp.str(name).word(bytes.len() as u64).raw(bytes);
        self.sources_read.push(fp.value());
        self.files_read.extend(path.map(|p| (p.to_path_buf(), Some(Doc::content_key(bytes)))));
    }

    /// The record of the family member whose record `path` would be, were
    /// it a file (`DIR/NAME.json`, the member `DIR/NAME.ink`:
    /// family::member_of), as `read-json` reads it there; None for no
    /// member's. The load rests on there being no file of its path, and on
    /// its directory's index; what it was read from is the record
    /// (Doc::sources).
    pub(crate) fn member_record(&mut self, path: &Path) -> Result<Option<serde_json::Value>, String> {
        if path.extension().is_none_or(|x| x != "json") {
            return Ok(None);
        }
        let Some(m) = family::member_of(&path.with_extension("ink"))? else { return Ok(None) };
        self.member_read(&m);
        self.files_read.push((path.to_path_buf(), None));
        Ok(Some(m.record))
    }

    /// Note what a load reading member `m` rests on: no file of its path,
    /// its directory's index, and, read from, its own record.
    pub(super) fn member_read(&mut self, m: &family::Member) {
        self.files_read.extend([(m.family.with_file_name(format!("{}.ink", m.name)), None), (m.index.clone(), Some(m.content))]);
        let mut fp = graph::Fingerprint::default();
        fp.str(&format!("{}#{}", shown(&m.index), m.name)).str(&m.record.to_string());
        self.sources_read.push(fp.value());
    }

    /// Note a path the load used for a file it read already, maybe by
    /// another path (a symlink, `..`): the load rests on what each path it
    /// used names, so a path retargeted later is seen.
    pub(crate) fn read_again(&mut self, path: &Path) {
        if !self.files_read.iter().any(|(f, _)| f == path) {
            self.files_read.push((path.to_path_buf(), std::fs::read(path).ok().map(|b| Doc::content_key(&b))));
        }
    }

    /// A file's content as `files_read` keys it.
    pub fn content_key(bytes: &[u8]) -> u128 {
        graph::Fingerprint::default().word(bytes.len() as u64).raw(bytes).value()
    }

    /// The files on disk the load depended on: each it read (a source it
    /// ran or data it read) with its content's key (`content_key`), and
    /// each it found absent (None: a `parts/...` import not beside its
    /// file, read from the catalog; a file it could not read, where the
    /// load failed: `load_file_reads`). With the code, they are all a load
    /// rests on: where each reads as it did, the load is the same.
    pub fn files_read(&self) -> &[(std::path::PathBuf, Option<u128>)] {
        &self.files_read
    }

    /// The fingerprint of everything the model was read from: the name and
    /// content of its own file, of every file it imports, transitively
    /// (std's text included), and of every file it reads (`read-json`,
    /// `import-stl`), as a set. A file in the catalog is named by its path
    /// there, so the same catalog gives the same fingerprint wherever it
    /// is. With the kernel's code, it is all a load of the model rests on.
    pub fn sources(&self) -> u128 {
        let mut read = self.sources_read.clone();
        read.sort_unstable();
        read.dedup();
        let mut fp = graph::Fingerprint::default();
        fp.word(read.len() as u64);
        read.iter().for_each(|k| _ = fp.key(*k));
        fp.value()
    }

    /// Load source text with the packages of `registry`.
    pub fn load_str_in(registry: &'static Registry, src: &str, name: &str, base: &Path, overrides: &[(String, Given)]) -> R<Doc> {
        Self::load_named(registry, src, name, Path::new(name), base, overrides, false).map_err(|(e, _)| e)
    }

    /// Load source text read from `path`, named `name` to the user; when
    /// it fails, with what it read before failing (`load_file_reads`).
    fn load_named(
        registry: &'static Registry,
        src: &str,
        name: &str,
        path: &Path,
        base: &Path,
        overrides: &[(String, Given)],
        claim_materials: bool,
    ) -> Result<Doc, (InkError, Vec<(std::path::PathBuf, Option<u128>)>)> {
        let origin = Origin { src: src.to_string(), name: name.to_string(), path: path.to_path_buf(), base: base.to_path_buf(), member: None, import_dirs: Vec::new() };
        Self::load_origin(registry, origin, overrides, None, claim_materials)
    }

    /// The file run loaded again, the arguments `target`'s instances pass
    /// (a top-level assembly of the file) built as coordinates (`coords`):
    /// what a package's answer differentiates placements along (a dynamics
    /// package built on the kernel's, say). Every other
    /// load builds them as written, so nothing else reads the coordinates'
    /// params, and no package adds what it adds once a file is read
    /// (Package::finish).
    pub fn with_coordinates(&self, target: &str) -> R<Doc> {
        self.reloader().ok_or_else(|| InkError::msg("only a document loaded from a source can be loaded again"))?.with_coordinates(target)
    }

    /// Whether it is a load `with_states` (Doc::with_states): the mechanism
    /// at any state, its moving coordinates' values, rates and accelerations
    /// params.
    pub fn states(&self) -> bool {
        self.states
    }

    /// How to load the file run again (`with_coordinates`), apart from the
    /// document: for a provider, which works without it.
    pub fn reloader(&self) -> Option<Reload> {
        Some(Reload { registry: self.registry, origin: self.source.clone()?, overrides: self.overrides.clone() })
    }

    pub(super) fn load_origin(
        registry: &'static Registry,
        origin: Origin,
        overrides: &[(String, Given)],
        coordinates_of: Option<(String, bool)>,
        claim_materials: bool,
    ) -> Result<Doc, (InkError, Vec<(std::path::PathBuf, Option<u128>)>)> {
        let dirs = origin.import_dirs.clone();
        importing(&dirs, || Self::load_origin_here(registry, origin, overrides, coordinates_of, claim_materials))
    }

    fn load_origin_here(
        registry: &'static Registry,
        origin: Origin,
        overrides: &[(String, Given)],
        coordinates_of: Option<(String, bool)>,
        claim_materials: bool,
    ) -> Result<Doc, (InkError, Vec<(std::path::PathBuf, Option<u128>)>)> {
        let mut d = Doc::empty(registry);
        d.claim_materials = claim_materials;
        d.overrides = overrides.to_vec();
        d.file = Rc::from(origin.name.as_str());
        (d.coordinates_of, d.states) = coordinates_of.map_or((None, false), |(t, s)| (Some(t), s));
        // Known while reading: a provider made then (a simulation) loads it
        // again.
        d.source = Some(origin.clone());
        d.member = origin.member.clone();
        let mut read = d.run_source_as(&origin.src, &origin.name, Some(&origin.path), &origin.base, None, Vec::new()).map(|_| ());
        d.parsed.clear();
        // Stopped at a refused case (Doc::refused): read no further.
        if read.is_ok() && !d.refused.is_empty() {
            return Ok(d);
        }
        // A combination of cases (ink require's) that does not declare a
        // param whose cases it takes: one some combinations make.
        if read.is_ok()
            && let Some((k, _)) = d.overrides.iter().find(|(_, g)| g.case.is_some() && g.text.is_none() && g.value.is_nan())
        {
            read = Err(InkError::msg(crate::claims::cases_not_everywhere(k, "is declared at the params' values but not at this combination")));
        }
        if read.is_ok() {
            d.rate_joints();
        }
        // What reads names after loading (a region on the command line)
        // reads the file's.
        d.scope = Scope { file: 1 };
        if read.is_ok() && d.coordinates_of.is_none() {
            read = registry.packages.iter().try_for_each(|p| {
                let mut r = Reader::new(&mut d);
                p.finish(&mut r)
            });
        }
        // A solve that does not solve with every other one: the error, once
        // all are read (Doc::solve_at_load); one depending on t, if what
        // rests on it was used at load time, where t = 0.
        let unsolved = d.unsolved.take().or_else(|| {
            let mut used: Vec<(usize, String)> = d.materialized.borrow().iter().map(|(&k, site)| (k, site.clone())).collect();
            used.sort();
            used.into_iter().find_map(|(k, site)| {
                let (_, e) = d.unsolved_t.iter().find(|(i, _)| d.solves_under(&[k]).contains(i))?;
                Some(InkError {
                    hint: Some(format!(
                        "'{}' is used at load time at {site}, where t = 0. Use it as a field (solved again at each moment), or make the solve close at t = 0",
                        d.params[k].name
                    )),
                    ..e.clone()
                })
            })
        });
        // The directions at the params' values (Doc::unit), once every
        // solve is read; those depending on t at each moment a view is
        // made at.
        let read = read.and_then(|()| if unsolved.is_none() { d.check_conditions(&d.load_ctx(), None, &[], true, None) } else { Ok(()) });
        let failed = match (read, unsolved) {
            (Ok(()), None) => return Ok(d),
            (Ok(()), Some(u)) => u,
            (Err(e), None) => e,
            (Err(mut e), Some(u)) => {
                let note = format!("before this, at {}:{}: {}. What reads its unknowns after it gets no solution", u.file, u.line, u.message);
                e.hint = Some(e.hint.map_or(note.clone(), |h| format!("{h}; {note}")));
                e
            }
        };
        Err((failed, d.files_read))
    }

    /// Context for load-time evaluation: parameters declared so far,
    /// overridden, with solved params re-solved for the overrides.
    pub(super) fn load_ctx(&self) -> Ctx {
        let mut params = self.load_params();
        for i in (0..self.solves.len()).filter(|&i| self.provided_always(i)) {
            if self.provide_one(i, &mut params, 0.0).is_err() {
                break;
            }
        }
        Ctx::new(0.0, params)
    }

    /// A load-time calculation needs only the providers its inputs read,
    /// transitively. A geometry query unrelated to density must not measure
    /// every published-mass body read so far.
    pub(super) fn load_ctx_for(&self, roots: &[NodeId]) -> Ctx {
        let mut params = self.load_params();
        for i in self.solves_read(roots) {
            if self.provide_one(i, &mut params, 0.0).is_err() {
                break;
            }
        }
        Ctx::new(0.0, params)
    }

    /// The params' load-time values before solving: defaults and overrides.
    pub(super) fn load_params(&self) -> Vec<f64> {
        let mut params: Vec<f64> = self.params.iter().map(|p| p.default).collect();
        for (k, v) in &self.overrides {
            if let Some(i) = self.params.iter().position(|p| p.name == *k) {
                params[i] = v.value;
            }
        }
        params
    }

    /// Run a file's forms (read from `path`, named `name` to the user); with
    /// a namespace, every top-level name it defines is `ns.name` in the
    /// document (Doc::public), and `name` as the file reads it.
    pub(super) fn run_source_as(&mut self, src: &str, name: &str, path: Option<&Path>, base: &Path, ns: Option<&str>, given: Vec<ImportValue>) -> R<usize> {
        // A member is read from its family's file (family::member_of).
        let member = path.and(self.member.take());
        let id = if let Some(path) = path {
            let mut f = SourceFile::new(Rc::from(name), canonical(path), ns.map(str::to_string), base.to_path_buf());
            f.given = given;
            f.member = member.clone();
            self.files.push(f);
            self.files.len() as u32 - 1
        } else {
            0
        };
        self.read(name, member.as_ref().map_or(path, |m| Some(m.family.as_path())), src.as_bytes());
        if let Some(m) = &member {
            self.member_read(m);
        }
        let saved = std::mem::replace(&mut self.file, Rc::from(name));
        let saved_scope = std::mem::replace(&mut self.scope, Scope { file: id });
        self.running.push(id);
        let forms = if path.is_none() {
            std_forms().clone()
        } else if let Some(forms) = self.parsed.get(src) {
            forms.clone()
        } else {
            let forms: std::sync::Arc<[Sexp]> =
                parse(src).map_err(|e| InkError { message: e.msg, file: name.to_string(), line: e.line, col: e.col, hint: None, called_from: Vec::new(), open: None })?.into();
            self.parsed.insert(src.to_owned(), forms.clone());
            forms
        };
        // Where the file defines each value, to say so of a read before it.
        let ahead = forms
            .iter()
            .filter_map(|f| match f.list() {
                Some([h, n, ..]) if matches!(h.sym(), Some("def" | "fn" | "param" | "members")) => Some((n.sym()?.to_string(), f.line)),
                _ => None,
            })
            .collect();
        self.files[id as usize].ahead = ahead;
        self.files[id as usize].declared = top_names(&forms, self.registry);
        let r = self.run_forms(&forms, base);

        self.file = saved;
        self.scope = saved_scope;
        self.running.pop();
        r.map(|_| id as usize)
    }

    /// Read a file's top-level forms in order, each in the names defined
    /// before it. A `for` is expanded when it is reached, so its list reads
    /// the file's names too, and the forms it makes are read there, under
    /// what its list rested on (Doc::under), a `for` among them expanded as
    /// it is reached in turn.
    /// At a combination of cases a text param refuses (Doc::refused) the
    /// load reads no further: what follows might read the text.
    fn run_forms(&mut self, forms: &[Sexp], base: &Path) -> R<()> {
        for f in forms {
            if !self.refused.is_empty() {
                return Ok(());
            }
            if let Some((outer, var, list, body)) = author::for_parts(f) {
                let (var, values) = self.at_form(f, |d| d.for_values(var, list))?;
                for (val, text) in values {
                    let made = author::for_made(body, outer, &var, val, text);
                    self.run_forms(&made, base)?;
                }
                continue;
            }
            self.at_form(f, |d| d.top(f, base))?;
        }
        Ok(())
    }

    /// Run a top-level form, or a `for`'s list: what it makes a number at
    /// load time is used there (`form_at`: Doc::materialized), and the
    /// undefined cases those numbers pass through are reported, never
    /// silent (`load_eval`).
    fn at_form<T>(&mut self, f: &Sexp, run: impl FnOnce(&mut Self) -> R<T>) -> R<T> {
        let head: Vec<String> = f.list().unwrap_or_default().iter().take(2).filter_map(|x| x.sym().map(str::to_string)).collect();
        self.form_at = format!("{}:{} ({} ...)", self.file, f.line, head.join(" "));
        self.form_name = head.last().cloned().unwrap_or_default();
        let r = run(self)?;
        for what in self.undefined_at_load.take() {
            let w = format!("{}:{}: {what} in a number used at load time. An undefined case gives 0 (see `ink spec`)", self.file, f.line);
            if !self.load_warnings.contains(&w) {
                self.load_warnings.push(w);
            }
        }
        Ok(r)
    }

    /// `f` (no x, y, z) as a number at load time, its undefined cases kept
    /// for `at_form` to report.
    pub(super) fn load_eval(&self, f: NodeId) -> f64 {
        let (tape, ctx) = (Tape::new(&self.g, f), self.load_ctx_for(&[f]));
        let cases = ink_core::eval::undefined_cases(&tape, &ctx, &[[0.0; 3]]);
        self.undefined_at_load.borrow_mut().extend(cases.iter().map(|u| u.what));
        eval_point(&tape, &ctx, [0.0; 3])
    }

    pub(super) fn err(&self, at: &Sexp, msg: impl Into<String>) -> InkError {
        InkError { message: msg.into(), file: self.file.to_string(), line: at.line, col: at.col, hint: comment_hint(at), called_from: Vec::new(), open: None }
    }

    pub(super) fn err_hint(&self, at: &Sexp, msg: impl Into<String>, hint: Option<String>) -> InkError {
        InkError { hint: comment_hint(at).or(hint), ..self.err(at, msg) }
    }

    /// The load stopped at `at`, sampling a grid (a load's one long step):
    /// the call's time limit (ink_kernel::stop) passed first.
    pub(super) fn stopped(&self, at: &Sexp) -> InkError {
        let hint = "give the call a larger time limit (timeout_s); the grids sampled are cached for the next call";
        self.err_hint(at, "the call's time limit (timeout_s) passed while loading, sampling this", Some(hint.into()))
    }
}
