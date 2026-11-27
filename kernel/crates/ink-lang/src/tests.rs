use super::*;

#[test]
fn a_load_time_scalar_runs_only_its_provider_dependencies() {
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    #[derive(Clone)]
    struct Counted {
        name: &'static str,
        owns: [usize; 1],
        read: usize,
        calls: Arc<AtomicUsize>,
    }
    impl Form for Counted {
        fn name(&self) -> &str {
            self.name
        }
        fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
            Some(Arc::new(self.clone()))
        }
    }
    impl ParamProvider for Counted {
        fn name(&self) -> &str {
            self.name
        }
        fn owns(&self) -> &[usize] {
            &self.owns
        }
        fn uses_t(&self) -> bool {
            true
        }
        fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            params[self.owns[0]] = 2.0 * params[self.read] + t;
            Ok(None)
        }
    }
    let mut d = Doc::load_str_in(&CORE, "(param p0 3) (param p1 3) (param p2 3) (param p3 3)", "test", Path::new("."), &[]).unwrap();
    let calls: Vec<_> = (0..3).map(|_| Arc::new(AtomicUsize::new(0))).collect();
    for (j, (name, read)) in [("first", 0), ("second", 1), ("unrelated", 0)].into_iter().enumerate() {
        let field = d.g.var(Var::Param(read as u32));
        let form = parse(&format!("(probe {name})")).unwrap().remove(0);
        d.keep_form("probe", "probe", &form, Arc::new(Counted { name, owns: [j + 1], read, calls: calls[j].clone() }), Reads { fields: vec![field], ..Reads::default() }).unwrap();
    }
    calls.iter().for_each(|c| c.store(0, Ordering::Relaxed));
    d.provided.borrow_mut().clear();
    let second = d.g.var(Var::Param(2));
    assert_eq!(d.load_eval(second), 12.0);
    d.overrides = vec![("p0".into(), 5.0.into())];
    assert_eq!(d.load_eval(second), 20.0);
    assert_eq!(d.load_eval(second), 20.0); // unchanged inputs do no work
    let constant = d.g.constant(7.0);
    assert_eq!(d.load_eval(constant), 7.0);
    for (head, field) in [("bounds", "(max (- (abs x) 1) (- (abs y) 1) (- (abs z) 1))"), ("profile-reach", "(max (- (abs x) 1) (- (abs y) 1))")] {
        let f = parse(&format!("({head} {field})")).unwrap().remove(0);
        d.eval_services(head, &f, &f.list().unwrap()[1..], &mut Vec::new()).unwrap();
    }
    assert_eq!(calls.iter().map(|c| c.load(Ordering::Relaxed)).collect::<Vec<_>>(), [2, 2, 0]);
    // A solver's owned values are its warm-start inputs too.
    d.params[1].default = 8.0;
    assert_eq!(d.load_eval(second), 20.0);
    assert_eq!(calls[0].load(Ordering::Relaxed), 3);
    let at = d.ctx(&[], 1.0).unwrap();
    assert_eq!(at.params[2], 15.0);
    let counts: Vec<_> = calls.iter().map(|c| c.load(Ordering::Relaxed)).collect();
    assert_eq!(d.ctx(&[], 1.0).unwrap().params[2], 15.0);
    assert_eq!(calls.iter().map(|c| c.load(Ordering::Relaxed)).collect::<Vec<_>>(), counts);
    assert_eq!(d.ctx(&[("p0".into(), 7.0.into())], 1.0).unwrap().params[2], 31.0);
}

#[test]
fn a_provider_on_demand_runs_only_where_its_params_are_read() {
    // A published mass's density: not measured at load, nor in a
    // context nothing reading it asks for (check, render, distance);
    // wanted (a mass weighed, a claim over a density), measured once:
    // what it set before is not an input to it.
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    #[derive(Clone)]
    struct Lazy {
        calls: Arc<AtomicUsize>,
    }
    impl Form for Lazy {
        fn name(&self) -> &str {
            "lazy"
        }
        fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
            Some(Arc::new(self.clone()))
        }
    }
    impl ParamProvider for Lazy {
        fn name(&self) -> &str {
            "lazy"
        }
        fn owns(&self) -> &[usize] {
            &[1]
        }
        fn uses_t(&self) -> bool {
            false
        }
        fn on_demand(&self) -> bool {
            true
        }
        fn warm_starts(&self) -> bool {
            false
        }
        fn provide(&self, params: &mut [f64], _: f64) -> Result<Option<String>, Refused> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            params[1] = 4.0;
            Ok(None)
        }
    }
    let mut d = Doc::load_str_in(&CORE, "(param p0 3) (param density 0)", "test", Path::new("."), &[]).unwrap();
    let calls = Arc::new(AtomicUsize::new(0));
    let form = parse("(probe lazy)").unwrap().remove(0);
    d.keep_form("probe", "probe", &form, Arc::new(Lazy { calls: calls.clone() }), Reads::default()).unwrap();
    assert_eq!((d.load_ctx().params[1], d.ctx(&[], 0.0).unwrap().params[1]), (0.0, 0.0));
    assert_eq!(calls.load(Ordering::Relaxed), 0, "provided at load or for a context nothing reading it asks for");
    let density = d.g.var(Var::Param(1));
    d.want(&[density]);
    let at: Vec<f64> = [0.0, 0.5].iter().map(|&t| d.ctx(&[], t).unwrap().params[1]).collect();
    let mut again = d.load_params();
    d.provide_one(0, &mut again, 0.0).unwrap();
    d.provide_one(0, &mut again, 0.0).unwrap();
    assert_eq!((at, again[1]), (vec![4.0, 4.0], 4.0));
    assert_eq!(calls.load(Ordering::Relaxed), 1, "measured more than once");
}

#[test]
fn a_provider_that_does_not_read_t_runs_once_over_a_motion() {
    // A published mass's density, at each moment a GLB or a motion
    // claim samples: measured once, not at every moment.
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    #[derive(Clone)]
    struct Counted {
        uses_t: bool,
        owns: [usize; 1],
        calls: Arc<AtomicUsize>,
    }
    impl Form for Counted {
        fn name(&self) -> &str {
            "counted"
        }
        fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
            Some(Arc::new(self.clone()))
        }
    }
    impl ParamProvider for Counted {
        fn name(&self) -> &str {
            "counted"
        }
        fn owns(&self) -> &[usize] {
            &self.owns
        }
        fn uses_t(&self) -> bool {
            self.uses_t
        }
        fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            params[self.owns[0]] = if self.uses_t { t } else { 4.0 };
            Ok(None)
        }
    }
    for uses_t in [false, true] {
        let mut d = Doc::load_str_in(&CORE, "(param p0 3)", "test", Path::new("."), &[]).unwrap();
        let calls = Arc::new(AtomicUsize::new(0));
        let form = parse("(probe counted)").unwrap().remove(0);
        d.keep_form("probe", "probe", &form, Arc::new(Counted { uses_t, owns: [0], calls: calls.clone() }), Reads::default()).unwrap();
        calls.store(0, Ordering::Relaxed);
        d.provided.borrow_mut().clear();
        let at: Vec<f64> = [0.0, 0.25, 0.5].iter().map(|&t| d.ctx(&[], t).unwrap().params[0]).collect();
        assert_eq!(at, if uses_t { [0.0, 0.25, 0.5] } else { [4.0; 3] });
        assert_eq!(calls.load(Ordering::Relaxed), if uses_t { 3 } else { 1 }, "uses_t {uses_t}");
    }
}

/// A directory of a test's own (`ink-TAG-PID` in the system's temporary
/// directory), removed when dropped, a test that panics too.
struct Scratch(std::path::PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let dir = std::env::temp_dir().join(format!("ink-{tag}-{}", std::process::id()));
        _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Scratch(dir)
    }
}

impl std::ops::Deref for Scratch {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The language alone, with no packages.
static CORE: Registry = Registry::new(&[]);

#[test]
fn package_queries_report_a_target_outside_the_readers_scope() {
    let dir = Scratch::new("package-query-scope");
    std::fs::write(dir.join("part.ink"), "(body b (shape 1))").unwrap();
    let mut doc = Doc::load_str_in(&CORE, r#"(assembly root) (import "part.ink" as part)"#, "main", &dir, &[]).unwrap();
    doc.scope = Scope { file: doc.files.iter().position(|f| f.ns.as_deref() == Some("part")).unwrap() as u32 };
    let mut reader = Reader::new(&mut doc);
    assert!(reader.clauses_under("b", "gauge").unwrap().is_empty());
    assert!(reader.relations_under("b").unwrap().is_empty());
    assert!(reader.clauses_under("root", "gauge").err().unwrap().contains("outside this file's scope"));
    assert!(reader.relations_under("root").err().unwrap().contains("outside this file's scope"));
}

/// A catalog of `parts/m/block.ink` importing `parts/m/shared.ink` (read
/// from the catalog's root, not beside it), shared's `size` being `size`;
/// and a model elsewhere importing the block.
fn catalog(dir: &Path, name: &str, size: f64) -> std::path::PathBuf {
    let root = dir.join(name);
    std::fs::create_dir_all(root.join("parts/m")).unwrap();
    std::fs::write(root.join("parts/m/shared.ink"), format!("(def size {size})")).unwrap();
    std::fs::write(root.join("parts/m/block.ink"), "(import \"parts/m/shared.ink\" as s) (def size s.size)").unwrap();
    root
}

/// A model that imports from a catalog it is not in.
fn away(dir: &Path) -> std::path::PathBuf {
    std::fs::create_dir_all(dir.join("away")).unwrap();
    let model = dir.join("away/model.ink");
    std::fs::write(&model, "(import \"parts/m/block.ink\" as b) (def size b.size)").unwrap();
    model
}

fn size(d: &Doc) -> f64 {
    let Some((_, Value::Field(f, _))) = d.defs.iter().find(|(n, _)| n == "size") else { panic!("no size") };
    d.load_eval(*f)
}

/// The directories a load is given are where its imports, and theirs,
/// are found, named by their paths there; no other load sees them; a
/// reload is given them again.
#[test]
fn a_load_given_import_dirs_finds_imports_there_alone() {
    let dir = Scratch::new("import-dirs");
    let (one, two, model) = (catalog(&dir, "one", 1.0), catalog(&dir, "two", 2.0), away(&dir));
    // Transitively, with no INK_PATH: block, and shared, which block imports.
    let d = Doc::load_file_importing(&CORE, &model, &[], std::slice::from_ref(&one)).unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!(size(&d), 1.0);
    let names: Vec<&str> = d.files.iter().skip(1).map(|f| &*f.name).collect();
    assert!(names.ends_with(&["parts/m/block.ink", "parts/m/shared.ink"]), "a file in a given directory is named by its path there: {names:?}");
    assert_eq!(d.given_import_dirs(), std::slice::from_ref(&one));
    // Loads one after another, each with its own: each reads its own, and a load given none reads neither.
    for (dirs, is) in [(vec![two.clone()], Some(2.0)), (vec![one.clone(), two.clone()], Some(1.0)), (vec![two.clone(), one.clone()], Some(2.0)), (vec![], None)] {
        let got = Doc::load_file_importing(&CORE, &model, &[], &dirs);
        assert_eq!(got.as_ref().ok().map(size), is, "{dirs:?}: {:?}", got.err());
    }
    let e = Doc::load_file_in(&CORE, &model, &[]).err().expect("the catalog is given to no plain load");
    assert!(e.message.starts_with("cannot read") && !import_dirs().contains(&one.canonicalize().unwrap()), "{}", e.message);
    // A load that fails gives its directories back too.
    let broken = dir.join("away/broken.ink");
    std::fs::write(&broken, "(import \"parts/m/none.ink\")").unwrap();
    let e = Doc::load_file_importing(&CORE, &broken, &[], std::slice::from_ref(&one)).err().unwrap();
    assert!(e.hint.as_deref().unwrap_or("").contains("the directories given its load, then INK_PATH's"), "{e:?}");
    assert!(!import_dirs().contains(&one.canonicalize().unwrap()));
    // At once, on threads of their own: each its own.
    let sizes: Vec<f64> = std::thread::scope(|s| {
        let runs: Vec<_> = (0..8)
            .map(|k| {
                let (dirs, model) = (if k % 2 == 0 { one.clone() } else { two.clone() }, &model);
                s.spawn(move || (0..4).map(|_| size(&Doc::load_file_importing(&CORE, model, &[], std::slice::from_ref(&dirs)).unwrap())).sum::<f64>())
            })
            .collect();
        runs.into_iter().map(|r| r.join().unwrap()).collect()
    });
    assert_eq!(sizes, [4.0, 8.0, 4.0, 8.0, 4.0, 8.0, 4.0, 8.0]);
    // Loaded again (another param's value): from the same directories.
    let v = crate::view::ModelView::load_importing(&CORE, &model, &[], 0.0, std::slice::from_ref(&two)).unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!(size(&v.reload(&[]).unwrap_or_else(|e| panic!("{e:?}")).doc), 2.0);
}

#[test]
fn a_json_file_is_read_once_and_shared() {
    // A table read by several files is one document in memory: every
    // data value read from it, and `get` into it, is a place in it, not
    // a copy (each catalog part held a whole copy of its family's
    // table: 2.4 MB a part).
    let dir = Scratch::new("json-shared");
    std::fs::write(dir.join("table.json"), r#"{"sizes": {"a/b": {"d": 3, "tags": [{"k": 1}]}}}"#).unwrap();
    let reader = r#"(def table (read-json "table.json")) (def row (get (get table "sizes") "a/b")) (def d (get row "d")) (def first (nth (get row "tags") 0))"#;
    std::fs::write(dir.join("a.ink"), reader).unwrap();
    std::fs::write(dir.join("b.ink"), reader).unwrap();
    let d = Doc::load_str_in(&CORE, r#"(import "a.ink" as p) (import "b.ink" as q)"#, "t", &dir, &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let data = |n: &str| match d.def(n) {
        Some(Value::Data(x)) => x.clone(),
        v => panic!("{n}: {v:?}"),
    };
    let all = [data("p.table"), data("q.table"), data("p.row"), data("q.row"), data("q.first")];
    assert!(all.iter().all(|x| Rc::ptr_eq(&x.root, &all[0].root)));
    assert_eq!(all[3].value()["d"], 3);
    assert_eq!(all[4].value()["k"], 1);
    assert!(matches!(d.def("q.d"), Some(Value::Field(f, _)) if d.g.as_const(*f) == Some(3.0)));
}

#[test]
fn get_takes_a_default_where_the_key_is_missing() {
    // A record's optional key (a belt's working tension, per its width
    // or whole) is read with a default, evaluated only where the key is
    // missing; without one, the error lists the keys there are.
    let dir = Scratch::new("json-default");
    std::fs::write(dir.join("belt.json"), r#"{"whole": 600, "at": "x"}"#).unwrap();
    let d = Doc::load_str_in(&CORE, r#"(def b (read-json "belt.json")) (def w (get b "whole" (get b "per-mm"))) (def p (get b "per-mm" 7))"#, "t", &dir, &[])
        .unwrap_or_else(|e| panic!("{e:?}"));
    let num = |n: &str| match d.def(n) {
        Some(Value::Field(f, _)) => d.g.as_const(*f),
        v => panic!("{n}: {v:?}"),
    };
    assert_eq!((num("w"), num("p")), (Some(600.0), Some(7.0)));
    let e = Doc::load_str_in(&CORE, r#"(def w (get (read-json "belt.json") "per-mm"))"#, "t", &dir, &[]).err().expect("a missing key");
    assert!(e.message == "no key 'per-mm'" && e.hint.as_deref() == Some("its keys: whole, at"), "{e:?}");
}

#[test]
fn a_function_reads_paths_beside_its_own_file() {
    // A path a function reads is relative to the file it is written
    // in, wherever it is called from: a caller with a table.json of its
    // own still gets the library's (docs/IR.md, names are lexical).
    let dir = Scratch::new("fn-path");
    std::fs::create_dir_all(dir.join("lib")).unwrap();
    std::fs::write(dir.join("lib/lookup.ink"), r#"(fn lookup (k) (get (read-json "table.json") k))"#).unwrap();
    std::fs::write(dir.join("lib/table.json"), r#"{"member": 7}"#).unwrap();
    std::fs::write(dir.join("table.json"), r#"{"member": 99}"#).unwrap();
    let d = Doc::load_str_in(&CORE, r#"(import "lib/lookup.ink" as l) (def v (l.lookup "member"))"#, "t", &dir, &[]).unwrap_or_else(|e| panic!("{e:?}"));
    assert!(matches!(d.def("v"), Some(Value::Field(f, _)) if d.g.as_const(*f) == Some(7.0)));
}

#[test]
fn a_files_own_def_may_take_a_built_ins_name() {
    // From its definition on, a file's own def, fn or param means
    // itself where a built-in has its name (names are lexical); std,
    // and the files it imports, keep the built-in.
    let dir = Scratch::new("own-built-in-name");
    let load = |src: &str| Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "t", &dir, &[]);
    let num = |d: &Doc, n: &str| match d.def(n) {
        Some(Value::Field(f, _)) => eval_point(&Tape::new(&d.g, *f), &Ctx::new(0.0, vec![]), [0.0; 3]),
        v => panic!("{n}: {v:?}"),
    };
    // A gift tag, and a def named max that std's clamp (min and max)
    // does not read.
    let d = load("(fn tag (w h) (box w h 2))\n(body gift (shape (tag 60 30)))\n(def max 4)\n(def c (clamp 9 0 5))\n(def m (+ max 1))").unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!((num(&d, "c"), num(&d, "m")), (5.0, 5.0));
    let gift = d.body("gift").unwrap().shape;
    assert_eq!(eval_point(&Tape::new(&d.g, gift), &Ctx::new(0.0, vec![]), [0.0; 3]), -1.0);
    // A built-in passed as a value from another file is the built-in
    // where it is called, in a file whose own len is another.
    std::fs::write(dir.join("count.ink"), "(def count len)").unwrap();
    let d = load("(import \"count.ink\" as c)\n(fn len (v) 99)\n(def n (c.count (vec 1 2 3)))\n(def ns (apply + (map c.count (vec (vec 1) (vec 1 2)))))\n(def own (len 0))")
        .unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!((num(&d, "n"), num(&d, "ns"), num(&d, "own")), (3.0, 3.0, 99.0));
    // A `for` over vecs writes them with the built-in vec.
    std::fs::write(dir.join("points.ink"), "(def points (vec (vec 1 2) (vec 3 4)))").unwrap();
    let d = load("(import \"points.ink\" as q)\n(fn vec (a) 0)\n(for p q.points (def s-$p (apply + p)))").unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!(num(&d, "s-3-4"), 7.0);
    // One meaning per file: the built-in read before the file's own, or
    // an import's def of the name, is an error naming both.
    let e = load("(def a (len (vec 1)))\n(fn len (v) 0)").err().unwrap();
    assert!(e.message.contains("the built-in len") && e.message.contains("line 3"), "{e:?}");
    std::fs::write(dir.join("label.ink"), "(fn tag (w) w)").unwrap();
    let e = load("(import \"label.ink\")\n(def a (tag 1))").err().unwrap();
    assert!(e.message.contains("the built-in tag and ") && e.message.contains("label.ink's def at line 1") && e.hint.as_deref().unwrap_or("").contains("namespace"), "{e:?}");
    let d = load("(import \"label.ink\" as g)\n(def a (g.tag 3))").unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!(num(&d, "a"), 3.0);
    // The built-in called where the file's own is no such function.
    for src in ["(def tag \"gift\")\n(body b (shape (tag \"lid\" (box 1 1 1))))", "(fn tag (w) (box w w w))\n(def a (tag \"lid\" (box 1 1 1)))"] {
        let e = load(src).err().unwrap();
        assert!(e.hint.as_deref().unwrap_or("").contains("not the built-in tag: rename it"), "{src}: {e:?}");
    }
    // The constants are read, not called: no file renames them.
    assert!(load("(def x 3)").err().unwrap().message.contains("'x' is a built-in name"));
}

#[test]
fn a_built_ins_name_means_one_thing_where_each_read_is_written() {
    // A read means what is visible where it is written: a function or
    // lambda written before the file's own len reads the built-in, so
    // calling it after is an error, as reading it before is; a `for`'s
    // own len after the built-in was read is one too. A user's own
    // name never reaches the head ink writes for a built-in.
    let load = |src: &str| Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "t", Path::new("."), &[]);
    for (src, says) in [
        ("(fn count-items (v) (len v))\n(fn len (v) 99)\n(def n (count-items (vec 1 2)))", "the built-in len, but this file defines its own 'len' at line 3"),
        ("(def count-items (lambda (v) (len v)))\n(fn len (v) 99)\n(def n (count-items (vec 1 2)))", "the built-in len, but this file defines its own 'len' at line 3"),
        ("(fn count-items (v) (len v))\n(def before (count-items (vec 1 2)))\n(for k (vec 1) (fn len (v) 99))", "read 'len' as the built-in len at line 2, before its own"),
    ] {
        let e = load(src).err().unwrap_or_else(|| panic!("{src}: loaded"));
        assert!(e.message.contains(says), "{src}: {e:?}");
    }
    let d = load("(fn builtin:len (v) 5)\n(def a (builtin:len 0))").unwrap_or_else(|e| panic!("{e:?}"));
    assert!(matches!(d.def("a"), Some(Value::Field(f, _)) if d.g.as_const(*f) == Some(5.0)));
}

#[test]
fn no_text_spells_the_head_ink_writes_for_a_built_in() {
    // A name a `for` makes from text (a string, a text param) is a
    // name like any other, whatever it spells: the head ink writes for
    // a built-in is no name, so none reaches it.
    let models = [
        "(fn len (v) 99)\n(for k (vec \"built-in;len\") (def a ($k (vec 1 2))))\n(def n (+ (* 100 a) (len 0)))",
        "(for k (vec \"built-in;len\") (def a ($k (vec 1 2))))\n(fn len (v) 99)\n(def n (+ (* 100 a) (len 0)))",
        "(param head \"built-in;len\")\n(fn len (v) 99)\n(for k (vec head) (def a ($k (vec 1 2))))\n(def n (+ (* 100 a) (len 0)))",
        "(param head \"built-in;len\" (cases \"built-in;len\"))\n(for k (vec head) (def a ($k (vec 1 2))))\n(for k (vec 1) (fn len (v) 99))\n(def n (+ a (len 0)))",
    ];
    for src in models {
        let e = Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "t", Path::new("."), &[]).err().unwrap_or_else(|| panic!("{src}: loaded"));
        assert!(e.message.contains("unknown function 'built-in;len'"), "{src}: {e:?}");
    }
}

#[test]
fn constant_turns_and_their_inverse_store_their_enclosures() {
    let src = "(body b (shape x)) (assembly root (instance a b (at (rot 0 0 1 30) (move 0 1 0))))";
    let mut d = Doc::load_str_in(&CORE, src, "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let at = d.assembly("root").unwrap().instances[0].at;
    let inverse = d.inverse(&at);
    for (node, exact) in [(at.r[1][0], 0.5), (inverse.t[0], -0.5)] {
        let (_, bound) = d.g.enclosed_constant(node).expect("a constant coefficient");
        assert_eq!(bound, ink_core::interval::I::point(exact));
    }
    // R^T (1,1,0) has first entry sqrt(2). Its negative lies
    // strictly between these adjacent floats.
    let src = "(body b (shape x)) (assembly root (instance a b (at (rot-z 45) (move 1 1 0))))";
    let mut d = Doc::load_str_in(&CORE, src, "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let at = d.assembly("root").unwrap().instances[0].at;
    let inverse = d.inverse(&at);
    let (_, bound) = d.g.enclosed_constant(inverse.t[0]).unwrap();
    let lower = -std::f64::consts::SQRT_2;
    assert!(bound.contains(lower) && bound.contains(lower.next_up()), "{bound:?}");
    let other = d.g.constant(17.0);
    for coefficient in turn_sin_cos(&mut d.g, other) {
        let (point, bound) = d.g.enclosed_constant(coefficient).unwrap();
        assert!(bound.lo < bound.hi && bound.contains(point));
    }
    let quarter = d.g.constant(90.0);
    assert_eq!(turn_sin_cos(&mut d.g, quarter).map(|n| d.g.as_const(n)), [Some(1.0), Some(0.0)]);
}

#[test]
fn a_turn_keeps_its_axis_exact() {
    // A turn about z leaves z alone: its matrix's z row and column are
    // the constants 0 and 1, so a frame on a spinning part reads t only
    // where it moves. A z entry of (1 - c) + c would read t, and a
    // claim on a spinning gear's axis or centre (a relation's) would be
    // stepped through every turn of it instead of checked once.
    let src = "(body b (shape (- x 1))) (assembly root (instance a b (at (move 0 0 5)) (turn z (* 360deg/s t))))";
    let d = Doc::load_str_in(&CORE, src, "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let at = d.assembly("root").unwrap().instances[0].at;
    let c = |e: NodeId| d.g.as_const(e);
    assert_eq!([at.r[2][2], at.r[0][2], at.r[1][2], at.r[2][0], at.r[2][1]].map(c), [Some(1.0), Some(0.0), Some(0.0), Some(0.0), Some(0.0)]);
    assert_eq!(at.t.map(c), [Some(0.0), Some(0.0), Some(5.0)]);
}

#[test]
fn family_functions_share_syntax_and_keep_each_parts_scope() {
    // Each namespace keeps its own bindings, but a family function's
    // parsed body (often much larger than its environment) is shared.
    let dir = Scratch::new("family-once");
    std::fs::write(dir.join("family.ink"), "(fn twice (v) (* 2 v)) (def k (twice 4))").unwrap();
    std::fs::write(dir.join("sized.ink"), "(param size 3) (fn sized () size)").unwrap();
    for (i, n) in [(1, 1), (2, 5)] {
        std::fs::write(dir.join(format!("part{i}.ink")), format!("(import \"family.ink\") (import \"sized.ink\") (def v (+ (twice {n}) (sized) k))")).unwrap();
    }
    let d = Doc::load_str_in(&CORE, r#"(import "part1.ink" as p) (import "part2.ink" as q)"#, "t", &dir, &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let v = |n: &str| match d.def(n) {
        Some(Value::Field(f, _)) => d.g.as_const(*f),
        v => panic!("{n}: {v:?}"),
    };
    assert_eq!((v("p.v"), v("q.v"), v("q.k")), (None, None, Some(8.0)));
    let (p, q) = (d.function("p.twice").unwrap(), d.function("q.twice").unwrap());
    assert_ne!(p.scope, q.scope);
    assert!(std::ptr::eq(p.body.list().unwrap(), q.body.list().unwrap()));
    let names: Vec<&str> = d.params.iter().map(|p| p.name.as_str()).collect();
    assert_eq!(names, ["p.size", "q.size"]);
    let ctx = |q: f64| Ctx::new(0.0, vec![3.0, q]);
    let at = |n: &str, c: &Ctx| match d.def(n) {
        Some(Value::Field(f, _)) => ink_core::eval::eval_point(&ink_core::tape::Tape::new(&d.g, *f), c, [0.0; 3]),
        v => panic!("{n}: {v:?}"),
    };
    assert_eq!((at("p.v", &ctx(3.0)), at("q.v", &ctx(4.0))), (13.0, 22.0));
    // A part's own name cannot shadow its family's, as when it ran.
    std::fs::write(dir.join("part3.ink"), "(import \"family.ink\") (def k 1)").unwrap();
    let e = Doc::load_str_in(&CORE, r#"(import "part1.ink" as p) (import "part3.ink" as r)"#, "t", &dir, &[]).err().expect("k defined twice");
    assert!(e.message.contains("already defined"), "{e:?}");
}

#[test]
fn a_family_shape_is_named_under_each_part_that_imports_it() {
    // A family's spatial def names its surface under the namespace it
    // is read in: each part's own name, not the first part's.
    let dir = Scratch::new("family-names");
    std::fs::write(dir.join("family.ink"), "(def orb (- x 2)) (fn make () orb)").unwrap();
    for p in ["a", "b"] {
        std::fs::write(dir.join(format!("{p}.ink")), "(import \"family.ink\") (def part (make))").unwrap();
    }
    let d = Doc::load_str_in(&CORE, r#"(import "a.ink" as a) (import "b.ink" as b)"#, "t", &dir, &[]).unwrap_or_else(|e| panic!("{e:?}"));
    // The names over each part's surface, outermost first.
    let named = |n: &str| {
        let Some(&Value::Field(mut f, _)) = d.def(n) else { panic!("{n}") };
        let mut out = Vec::new();
        while let ink_core::graph::Node::Tag(t, a) = d.g.nodes[f as usize] {
            out.push(d.g.tag_names[t as usize].as_str());
            f = a;
        }
        out
    };
    assert_eq!((named("a.part"), named("b.part")), (vec!["a.part", "a.orb"], vec!["b.part", "b.orb"]));
}

#[test]
fn grids_over_package_results_key_by_their_values() {
    // A grid sampled from a package's result is keyed by what it reads
    // (no package is named): the result's values are part of the field's
    // fingerprint, so a solver that gives other values (here: a probe
    // package's version) misses, and one giving the same values hits.
    use std::sync::atomic::{AtomicU32, Ordering};
    static VERSION: AtomicU32 = AtomicU32::new(1);
    struct Probe;
    struct ProbeForm(String);
    struct Solved(f32);
    impl Package for Probe {
        fn name(&self) -> &'static str {
            "probe"
        }
        fn forms(&self) -> &'static [&'static str] {
            &["probe"]
        }
        fn read(&self, r: &mut Reader, _form: &Sexp, _head: &str, args: &[Sexp]) -> R<std::sync::Arc<dyn Form>> {
            Ok(std::sync::Arc::new(ProbeForm(r.name(&args[0])?)))
        }
    }
    impl Form for ProbeForm {
        fn name(&self) -> &str {
            &self.0
        }
        fn results(&self) -> Vec<ResultSpec> {
            vec![ResultSpec::scalar("v", "")]
        }
        fn run(&self, _doc: &mut Doc, _ctx: &Ctx, _opts: &RunOpts) -> R<Box<dyn Outcome>> {
            Ok(Box::new(Solved(VERSION.load(Ordering::SeqCst) as f32)))
        }
    }
    impl Outcome for Solved {
        fn field(&self, _: &str) -> Option<crate::grid::GridData> {
            Some(sampled_result(vec![self.0; 8], [2, 2, 2], [-5.0; 3], 10.0, 0.0))
        }
    }
    static PROBES: Registry = Registry::new(&[&Probe]);
    let key = |version: u32| {
        VERSION.store(version, Ordering::SeqCst);
        let src = "(probe p)\n(def f (max (- (abs x) 5) (study-field p v)))";
        let d = Doc::load_str_in(&PROBES, src, "test", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
        let Some(Value::Field(f, _)) = d.def("f") else { panic!("f is not a field") };
        d.grid_key("grid-bake", Some(1.0), None, Some(*f), &[], &d.load_ctx().params, &[])
    };
    let (a, b) = (key(1), key(2));
    assert_ne!(a, b);
    assert_eq!(key(1), a);
}

/// A `for` among an instance's clauses reads its assembly's arguments:
/// a drive is limited at one supply and not at another.
#[test]
fn an_instance_clause_for_reads_the_assembly_arguments() {
    let d = Doc::load_str_in(
        &CORE,
        "(import \"std\") (body a (shape (box 1 1 1)))
        (assembly (drive (angle 0) (volts 24))
          (instance r a (at) (turn z angle) (for lo (if (< volts 20) (vec) (vec -5)) (limits lo 5))))
        (assembly root (instance low (drive (* 10 (/ t 1s)) 12)) (instance at (drive (* 10 (/ t 1s)))))",
        "t",
        Path::new("."),
        &[],
    )
    .unwrap();
    let limited = |path: &str| d.assemblies.iter().find(|a| a.name == path).unwrap().instances[0].joint.as_ref().unwrap().limits.is_some();
    assert!(!limited("root/low") && limited("root/at"));
}

/// A package's clauses (a joint's effort limits and stiffness, a
/// relation's coupling) are refused by an ink that loads no package
/// owning them, saying so; `(optional CLAUSE)` is left out there, so a
/// gear mesh's claims stand without its coupling.
#[test]
fn a_package_clause_is_refused_where_no_package_owns_it() {
    let load = |src: &str| Doc::load_str_in(&CORE, &format!("(import \"std\") (body a (shape (box 1 1 1)))\n{src}"), "t", Path::new("."), &[]);
    for (src, head) in [
        ("(assembly root (instance r a (at) (turn z 0) (effort-limits -5N*mm 5N*mm)))", "'effort-limits'"),
        ("(assembly root (instance r a (at) (slide x 0) (stiffness 10N/mm)))", "'stiffness'"),
        ("(assembly root (instance p a) (instance q a (at (move 9 0 0))) (relate p q (belt-drive 5 5 10N)))", "'couple'"),
    ] {
        let e = load(src).err().unwrap_or_else(|| panic!("{src} loads"));
        assert!(e.message.contains(head) && e.message.contains("a package this ink loads"), "{src}: {}", e.message);
    }
    let d = load(
        "(assembly root (instance p a) (instance q a (at (move 9 0 0)))
           (relate p q (lambda (fa fb sa sb) (relation r (optional (couple (vec) 0)) (data n 1)))))",
    )
    .unwrap_or_else(|e| panic!("{e:?}"));
    let k = d.assemblies[0].relations[0].kind.clone().unwrap();
    assert!(k.clauses.is_empty() && k.data.len() == 1);
    // A joint's, as a catalog stepper writes them: the kernel places it
    // without them; a kernel clause is never optional.
    let d = load(
        "(assembly root (instance r a (at) (turn z 0) (optional (stiffness 10N*mm))
           (for c (vec 0 1) (optional (effort-limits -5N*mm 5N*mm)))))",
    )
    .unwrap_or_else(|e| panic!("{e:?}"));
    let j = d.assemblies[0].instances[0].joint.as_ref().unwrap();
    assert!(j.clauses.is_empty(), "{:?}", j.clauses.len());
    let e = load("(assembly root (instance r a (at) (optional (turn z 0))))").err().unwrap();
    assert!(e.message.contains("for a package's clause"), "{}", e.message);
    // What it expands to is checked: a kernel clause in a `for` is
    // still the kernel's, in an instance and in a relation.
    for src in [
        "(assembly root (instance r a (at) (turn z 5) (optional (for k (vec 0) (limits -1 1)))))",
        "(assembly root (instance p a) (instance q a (at (move 9 0 0)))
           (relate p q (lambda (fa fb sa sb) (relation r (optional (for k (vec 0) (limits -1 1)))))))",
        "(assembly root (instance p a) (instance q a (at (move 9 0 0)))
           (relate p q (lambda (fa fb sa sb) (relation r (optional (for k (vec 0) (data n 1)))))))",
    ] {
        let e = load(src).err().unwrap_or_else(|| panic!("{src} loads"));
        assert!(e.message.contains("for a package's clause"), "{src}: {}", e.message);
    }
    // Where a loaded package reads joint or relation clauses, one none
    // reads (a misspelt stiffness) warns: dropping it there would drop
    // that package's check. Where none does (the kernel's ink, whose
    // packages read neither), the kernel stands without it, silently.
    struct Other(bool);
    impl Package for Other {
        fn name(&self) -> &'static str {
            if self.0 { "joints" } else { "forms-only" }
        }
        fn forms(&self) -> &'static [&'static str] {
            &[]
        }
        fn read(&self, _: &mut Reader, f: &Sexp, _: &str, _: &[Sexp]) -> R<std::sync::Arc<dyn Form>> {
            unreachable!("{f}")
        }
        fn joint_clauses(&self) -> &'static [&'static str] {
            if self.0 { &["stiffness"] } else { &[] }
        }
        fn relation_clauses(&self) -> &'static [&'static str] {
            if self.0 { &["couple"] } else { &[] }
        }
    }
    static OTHER: Registry = Registry::new(&[&Other(true)]);
    static FORMS_ONLY: Registry = Registry::new(&[&Other(false)]);
    let typo = "(import \"std\") (body a (shape (box 1 1 1)))
        (assembly root (instance p a) (instance q a (at (move 9 0 0)) (turn z 0) (optional (stifness 10)))
          (relate p q (lambda (fa fb sa sb) (relation r (optional (cuple (vec) 0))))))";
    let warned = |reg: &'static Registry| Doc::load_str_in(reg, typo, "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}")).load_warnings;
    let w = warned(&OTHER);
    assert!(w.iter().any(|w| w.contains("'stifness'")) && w.iter().any(|w| w.contains("'cuple'")), "{w:?}");
    assert!(warned(&CORE).is_empty() && warned(&FORMS_ONLY).is_empty());
}

/// The papercuts agents hit, each as
/// they wrote it: the error says what to write.
#[test]
fn errors_say_what_to_write() {
    let e = Doc::load_str_in(&CORE, "(def a (sqr x))", "t", Path::new("."), &[]).err().unwrap();
    assert!(e.message.contains("unknown function 'sqr'"));
    assert!(e.hint.unwrap().contains("sqrt"));
    let load = |src: &str| Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "t", Path::new("."), &[]);
    let hint = |src: &str| load(src).err().unwrap().hint.unwrap_or_default();
    // A CAD word gets std's name for it, not a near spelling ('hole'
    // is not 'fold'); a misspelt part of a name finds the name.
    let h = hint("(def a (hole 1 2))");
    assert!(h.contains("'subtract'") && h.contains("'cylinder-z'") && !h.contains("fold"), "{h}");
    assert!(hint("(def a (fillet 1 2))").contains("'fillet-union'"));
    assert!(hint("(def a (extrude 1 2))").contains("'extrude-z'"));
    assert!(hint("(def a (cylnder 1 2))").contains("'cylinder-z'"));
    assert_eq!(hint("(def a (foo 1 2))"), "");
    // Rounding a number is floor's; a top-level form inside an
    // expression says what is written there ('fn', not 'ln').
    assert!(hint("(def a (round 1.5))").contains("(floor (+ v 0.5))"));
    let h = hint("(def f (map (fn (a) (* 2 a)) (vec 1 2)))");
    assert!(h.contains("top-level form") && h.contains("(lambda (args) body)") && !h.contains("'ln'"), "{h}");
    // A function written as a def, three ways.
    assert_eq!(hint("(def (hole z) (move (cylinder-z 3 40) 0 6 z))"), "write (fn hole (z) body)");
    assert_eq!(hint("(def peg (x) (move (sphere 2) x 25 0))"), "write (fn peg (x) body)");
    assert_eq!(hint("(def slot (fn (x d) (move (box d 2 2) x 0 0)))"), "write (fn slot (x d) body)");
    // A def holding a function is called by name, as a local one is.
    let d = load("(def mark (lambda (z) (move (cylinder-z 1 1) 0 0 z)))\n(def marks (min (mark 10) (mark 20)))").unwrap_or_else(|e| panic!("{e:?}"));
    assert!(matches!(d.def("marks"), Some(Value::Field(..))));
    assert!(load("(def w 5)\n(def a (w 3))").err().unwrap().message.contains("'w' is a def holding a field, not a function"));
    // A constant's name lists the built-in names.
    for n in ["t", "x", "pi"] {
        let h = hint(&format!("(param {n} 3)"));
        assert!(h.starts_with("the built-in names: x y z t pi inf;") && h.split([' ', ';']).any(|w| w == n), "{n}: {h}");
    }
    // Another language's comment: the error it makes says what one
    // is; the same words as data (a path, a glob) are no error.
    for src in ["// a bracket\n(def a 1)", "(def a 1) /* a */", "(def a (+ 1 //2\n 3))", "(def a (sphere 1 // r\n))"] {
        assert_eq!(hint(src), "a comment runs from ';' to the line's end", "{src}");
    }
    load("(def source (data (path //server/share/part.step)))\n(def stls (data (glob /*.stl)))").unwrap_or_else(|e| panic!("{e:?}"));
    // A name quoted or in parentheses.
    assert_eq!(hint("(body \"pen-cup\" (shape (box 1 1 1)))"), "a name is a bare word: pen-cup");
    assert_eq!(hint("(body \"Pill Box\" (shape (box 1 1 1)))"), "a name is a bare word: Pill-Box");
    // Words from other languages: a list is a vec; a loop in an
    // expression maps over one.
    assert!(hint("(def a (list 1 2))").contains("'vec'"));
    let h = hint("(def a (min (for i (range 3) (move (sphere 1) i 0 0))))");
    assert!(h.contains("'map'") && h.contains("(reduce min (map (lambda (i)"), "{h}");
    assert_eq!(hint("(body b (shape (box 1 1 1)))\n(assembly root (for i (range 2) (instance (motor-$i) b)))"), "a name is a bare word: motor-0");
    // A placement step written as a shape's turn.
    assert!(hint("(def b (rot-y 90 (cylinder-z 8 55)))").contains("(rotate-y f deg)"));
    // A shape's move or turn written as a placement or a step of one
    // a step's misspelling,
    // a step nested in a step; a function
    // misspelt there keeps its own hint.
    let plate = "(body plate (shape (box 40 40 5)))\n(fn seat (i) (at (move i 0 5)))\n";
    for (src, says) in [
        ("(body axle (shape (sphere 2)) (frame left-axle (move -24 0 0)))", "(at STEP ...)"),
        ("(assembly root (instance p plate (move 0 0 5)))", "(at STEP ...)"),
        ("(assembly root (instance p plate (at (move 0 0 5) (rotate-z 90))))", "(rot-z deg)"),
        ("(assembly root (instance p plate (at (move 0 0 5) (rotz 90))))", "(rot-z deg)"),
        ("(assembly root (instance p plate (sead 1)))", "'read-json', 'seat'?"),
    ] {
        assert!(hint(&format!("{plate}{src}")).contains(says), "{src}: {}", hint(&format!("{plate}{src}")));
    }
    // The model's own function making a frame, miscounted, is told so
    // only.
    for src in ["(assembly root (instance p plate (seat)))", "(assembly root (instance p plate (at (seat 2 3))))"] {
        let e = load(&format!("{plate}{src}")).err().unwrap();
        assert!(e.message.starts_with("'seat' takes 1 argument(s) (i)") && e.hint.is_none(), "{src}: {e:?}");
    }
    // A shape written bare in a body.
    assert_eq!(hint("(body b (box 1 1 1))"), "write (shape (box ...))");
    // A claim of equality says how to write one.
    assert!(hint("(param a 1)\n(require one (value a == 1))").contains("(value (abs (- a b)) <= tol)"));
    // A CAD word as a value is named as one called is.
    assert!(hint("(def s (reduce union (vec 1 2)))").contains("a union is the min of its shapes"));
    assert!(hint("(def c (or (> 1 0) (< 1 0)))").contains("(max a b) is a or b"));
    // A one-letter name undeclared is no typo of t.
    assert_eq!(hint("(body plate (shape (box w 40 3)))"), "declare it: (param w VALUE) or (def w EXPR)");
    // A name declared as a body or an assembly further on, or an
    // instance, is not to be declared again.
    for src in ["(body b (shape (box root 1 1)))\n(assembly root (instance x b))", "(body p (shape (box 4 4 1)))\n(assembly root (instance a p) (instance b p (at (move a 0 0))))"]
    {
        assert!(!hint(src).contains("declare it"), "{src}: {}", hint(src));
    }
    let e = load(&format!("{plate}(assembly root (instance r plate (at (move 0 0 14 (move 0 0 4) (move 0 0 8)))))")).err().unwrap();
    assert!(e.message == "the step 'move' takes 3 number(s) (dx dy dz), got 5" && e.col == 38, "{e:?}");
}
/// std's sheet refuses at load what it cannot build, saying what to
/// change: the sizes it cannot take, crossing bends, a bend written
/// the wrong way round, a frame in a bend's zone.
#[test]
fn sheet_refuses_what_it_cannot_build() {
    let load = |src: &str| Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "t", Path::new("."), &[]);
    let sheet = |th: &str, k: &str, bends: &str| {
        format!("(def S (sheet (polygon (vec (vec 0 0) (vec 60 0) (vec 60 40) (vec 0 40))) (vec) {th} {k} (vec {bends})))\n(def f (sheet-folded S))")
    };
    let says = |src: &str| load(src).err().map(|e| e.message).unwrap_or_default();
    let l = "(sheet-bend (vec 30 40) (vec 30 0) 90 2)";
    assert!(load(&sheet("2", "0.4", l)).is_ok());
    assert!(says(&sheet("0", "0.4", l)).contains("thickness th must be above 0"));
    assert!(says(&sheet("2", "-1.5", l)).contains("K-factor k") && says(&sheet("2", "1.2", l)).contains("0 to 1"));
    assert!(says(&sheet("2", "0.4", "(sheet-bend (vec 30 40) (vec 30 0) 90 0)")).contains("inside radius must be above 0"));
    assert!(says(&sheet("2", "0.4", "(sheet-bend (vec 30 40) (vec 30 0) 190 2)")).contains("at most 180"));
    // bend 2 runs across bend 1's zone; its midpoint lies outside it
    let cross = format!("{l} (sheet-bend (vec 20 35) (vec 60 35) 90 2)");
    assert!(says(&sheet("2", "0.4", &cross)).contains("crosses another bend's zone"));
    // walls on +x and +y, their lines written past the outline, crossing
    // in the notched corner's air: built
    let corner = "(def S (sheet (polygon (vec (vec 0 0) (vec 70 0) (vec 70 45) (vec 45 45) (vec 45 70) (vec 0 70))) (vec) 2 0.4 \
                  (vec (sheet-bend (vec 50 60) (vec 50 -5) 90 2) (sheet-bend (vec -5 50) (vec 60 50) 90 2))))\n(def f (sheet-folded S))";
    assert!(load(corner).is_ok(), "{:?}", load(corner).err());
    // a lip on the upright written b -> a: its left faces the base
    let lip = format!("{l} (sheet-bend (vec 50 0) (vec 50 40) 90 2)");
    assert!(says(&sheet("2", "0.4", &lip)).contains("swap a and b"));
    assert!(says(&format!("{}\n(def fr (sheet-frame S 30 20 1))", sheet("2", "0.4", l))).contains("in a bend's zone"));
    assert!(says(&sheet("2", "0.4", "(sheet-bend (vec 30 20) (vec 30 20) 90 2)")).contains("a and b must differ"));
    // two tabs on one edge and a bend beyond both: which flange is it on?
    let two = "(sheet-bend (vec 30 10) (vec 30 0) 90 2) (sheet-bend (vec 30 40) (vec 30 30) 45 2) (sheet-bend (vec 45 40) (vec 45 0) 90 2)";
    assert!(says(&sheet("2", "0.4", two)).contains("a bend's line lies beyond two bends neither of which lies beyond the other"));
}

/// Every expression of units_cases/cases.txt (an outside read's probes, and the cases it found wrong), after
/// the names of units_cases/names.ink: its number, and the unit it is said in ("-": a pure number; "?": no
/// quantity of the language's), or its error. expected.txt is what each must be.
#[test]
fn the_units_cases_come_out_as_expected() {
    static CORE: Registry = Registry::new(&[]);
    let names = include_str!("units_cases/names.ink");
    let mut got = String::new();
    for case in include_str!("units_cases/cases.txt").lines().filter(|l| !l.trim().is_empty()) {
        let answer = match Doc::load_str_in(&CORE, &format!("{names} (def v {case})"), "t", Path::new("."), &[]) {
            Err(e) => format!("ERROR {}", e.message.split(':').next().unwrap_or_default()),
            Ok(d) => match d.def("v").cloned() {
                Some(Value::Field(n, said)) => {
                    let unit = said.as_ref().map_or("-", |s| if std::ptr::eq(s.unit, &units::DIFFERENCE) { "K difference" } else { s.unit.quantity.unit() });
                    let show = |x: f64| format!("{:.6}", x).trim_end_matches('0').trim_end_matches('.').to_string();
                    // What a fact of it would be, where that is not what it is said in: a quantity, its angle in radians.
                    let read = said.as_ref().and_then(|s| s.read()).filter(|(_, times)| *times != 1.0);
                    let read =
                        read.map_or(String::new(), |(q, times)| format!(", read as {} {}", d.g.as_const(n).map_or("not constant".to_string(), |x| show(x * times)), q.unit()));
                    format!("{} {unit}{read}", d.g.as_const(n).map_or("not constant".to_string(), show))
                }
                _ => "not a number".to_string(),
            },
        };
        got += &format!("{case}  =>  {answer}\n");
    }
    if got != include_str!("units_cases/expected.txt") {
        panic!("the units cases differ from units_cases/expected.txt; they are now:\n{got}");
    }
}

/// A product is one value however it is grouped, ordered and reached: every choice of four of fourteen factors (two
/// angles, two turning rates, two lengths, a torque, a force, a time, a mass, an acceleration, a pressure, a bare number, a reciprocal), in
/// every order and every grouping, written with `*`, with each group a `def`, with each group a function's
/// argument, and through `fold`, `reduce` and `apply`, is the same number, said in the same unit with the same
/// powers, and read as the same quantity.
#[test]
fn a_product_is_one_value_however_it_is_grouped() {
    static CORE: Registry = Registry::new(&[]);
    const FACTORS: [&str; 14] = ["90deg", "2rad", "60rpm", "3rad/s", "2mm", "1in", "1.5N*m", "4N", "0.5s", "2kg", "3m/s2", "2MPa", "3", "(/ 1 2mm)"];
    // Every grouping of `items`, in their order: (the expression, the defs it needs), by `*` and by a function.
    fn groupings(items: &[&str], by: &str, defs: bool, n: &mut usize) -> Vec<(String, String)> {
        if let [one] = items {
            return vec![(one.to_string(), String::new())];
        }
        let mut all = Vec::new();
        for cut in 1..items.len() {
            for (a, da) in groupings(&items[..cut], by, defs, n) {
                for (b, db) in groupings(&items[cut..], by, defs, n) {
                    let product = format!("({by} {a} {b})");
                    if defs {
                        *n += 1;
                        all.push((format!("g{n}"), format!("{da}{db}(def g{n} {product})\n")));
                    } else {
                        all.push((product, String::new()));
                    }
                }
            }
        }
        all
    }
    fn orders(items: &mut Vec<usize>, k: usize, all: &mut Vec<Vec<usize>>) {
        if k == items.len() {
            return all.push(items.clone());
        }
        for i in k..items.len() {
            items.swap(k, i);
            orders(items, k + 1, all);
            items.swap(k, i);
        }
    }
    let (mut written, mut sets) = (0, 0);
    for a in 0..14 {
        for b in a + 1..14 {
            for c in b + 1..14 {
                for d in c + 1..14 {
                    let (mut all, mut src, mut n) = (Vec::new(), String::from("(fn mul (p q) (* p q))\n"), 0);
                    orders(&mut vec![a, b, c, d], 0, &mut all);
                    let mut names = Vec::new();
                    for order in &all {
                        let items: Vec<&str> = order.iter().map(|i| FACTORS[*i]).collect();
                        let flat = items.join(" ");
                        let mut ways = vec![(format!("(* {flat})"), String::new())];
                        ways.extend(["fold * 1", "reduce *", "apply *"].map(|by| (format!("({by} (vec {flat}))"), String::new())));
                        ways.extend(groupings(&items, "*", false, &mut n));
                        ways.extend(groupings(&items, "mul", false, &mut n));
                        ways.extend(groupings(&items, "*", true, &mut n));
                        for (value, defs) in ways {
                            names.push((format!("v{}", names.len()), value.clone()));
                            src += &format!("{defs}(def v{} {value})\n", names.len() - 1);
                        }
                    }
                    let doc = Doc::load_str_in(&CORE, &src, "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{}: {}", e.line, e.message));
                    let of = |name: &str| {
                        let Some(Value::Field(f, said)) = doc.def(name).cloned() else { panic!("{name}") };
                        let said = said.map(|s| (s.unit.quantity, s.dims().map(|d| d.0), s.read().map(|r| r.0)));
                        (doc.g.as_const(f).unwrap(), said)
                    };
                    let first = of("v0");
                    for (name, value) in &names {
                        let got = of(name);
                        assert!((got.0 - first.0).abs() <= 1e-12 * first.0.abs() && got.1 == first.1, "{value} is {got:?}, but {} is {first:?}", names[0].1);
                    }
                    written += names.len();
                    sets += 1;
                }
            }
        }
    }
    assert_eq!((sets, written), (1001, 1001 * 24 * 19));
}

/// A number's unit goes with its value as it is computed, and the slot
/// that uses the value checks it: through lets, calls, the branch an
/// `if` takes, functions passed and bound, sums, and products and
/// quotients, whose units combine into the quantity they make; not
/// through a branch not taken, nor into a product that is no quantity
/// of the language's.
#[test]
fn a_value_carries_the_unit_it_was_said_in() {
    static CORE: Registry = Registry::new(&[]);
    let fns = "(fn twice (p) (+ p p)) (fn choose (v) (if 1 v 1psi)) (fn other (v) (if 0 v 1psi)) (fn pass (f v) (f v)) (fn id (v) v)";
    let ids: String = (1..=17).map(|i| format!("(fn id{i} (v) ({} v))", if i == 1 { "id".into() } else { format!("id{}", i - 1) })).collect();
    let load = |mass: &str| {
        let src = format!("{fns} {ids} (param pp 50psi) (def r (ranged 1psi 2psi)) (body b (shape (- (max (abs x) (abs y) (abs z)) 1)) (attr mass {mass}))");
        Doc::load_str_in(&CORE, &src, "t", Path::new("."), &[]).err().map(|e| e.message)
    };
    let nested = |x: &str| (0..23).fold(x.to_string(), |a, _| format!("(twice {a})"));
    for mass in [
        "(let ((p 50psi)) p)".to_string(),
        nested("1psi"),
        "(other 5g)".into(),
        "(pass id 1psi)".into(),
        "(let ((g id)) (g 1psi))".into(),
        "(id17 1psi)".into(),
        "(nth (map id (vec 1psi)) 0)".into(),
        "(reduce + (vec 1psi 2psi))".into(),
        "(- (max pp 0) 1psi)".into(),
        "r".into(),
        "(* 2 1psi)".into(),
    ] {
        let e = load(&mass).unwrap_or_else(|| panic!("{mass}: accepted"));
        assert!(e.starts_with("attr mass takes a mass") && e.contains("psi (t:1)"), "{mass}: {e}");
    }
    for mass in ["(choose 5g)", "(* 2g/mm3 4mm3)", "(/ (* 2g 3mm) 1mm)", "(* 2 (/ 6g 3))", &nested("1g")] {
        assert_eq!(load(mass), None, "{mass}");
    }
    // What the language has no unit for, its dimension known, is no mass either.
    let e = load("(* 2s 1N)").unwrap();
    assert!(e.starts_with("attr mass takes a mass (g), but the value worked out at t:1 is g mm/s, not a mass"), "{e}");
    // A product is the quantity it makes: a length times a pressure is a stiffness, a force over an area a
    // pressure, a mass under gravity a force in newtons (ink works in grams: 1000 times fewer).
    for (mass, is) in
        [("(* 2mm 1psi)", "a stiffness"), ("(/ 10N 2mm2)", "a pressure or stress"), ("(* 10kg 9.81m/s2)", "a force"), ("(square 2mm)", "an area"), ("(pow 2mm 3)", "a volume")]
    {
        let e = load(mass).unwrap_or_else(|| panic!("{mass}: accepted"));
        assert!(e.starts_with("attr mass takes a mass") && e.contains(&format!("{is} worked out at t:1")), "{mass}: {e}");
    }
    let number = |value: &str| {
        let d = Doc::load_str_in(&CORE, &format!("(def v {value})"), "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{value}: {}", e.message));
        let Some(Value::Field(n, said)) = d.def("v").cloned() else { panic!() };
        (d.g.as_const(n).unwrap(), said.map(|s| s.unit.name))
    };
    assert_eq!(number("(* 10kg 9.81m/s2)").1, Some("N"));
    assert!((number("(* 10kg 9.81m/s2)").0 - 98.1).abs() < 1e-9);
    // A pure number is left where the units cancel. A turning rate times a time is an angle, in degrees, the
    // number it always was; a bare number over a turning rate is seconds per degree, no quantity ("?"), as is.
    assert_eq!([number("(/ 2in 4mm)"), number("(* 3000rpm 2s)"), number("(/ 360 60rpm)")], [(12.7, None), (36000.0, Some("deg")), (1.0, Some("?"))]);
    // An angle is a quantity, in degrees or radians; with its unit, a sine takes it as the angle it is.
    assert_eq!(
        [number("(/ 360deg 60rpm)"), number("(+ 90deg 0rad)"), number("(* (/ 124200deg 1397mm) 1397mm/s)")],
        [(1.0, Some("s")), (90.0, Some("deg")), (124200.0, Some("deg/s"))]
    );
    assert!((number("(sin 30deg)").0 - 0.5).abs() < 1e-12 && (number("(cos 1rad)").0 - 1f64.cos()).abs() < 1e-12 && number("(sin 30)").0 == 30f64.sin());
    // An angle left beside a quantity stays in the arithmetic, however the product is grouped, and is taken in
    // radians where the value is read: torque times a turning rate as power, a radius times one as a speed, a
    // radius times an angle as an arc.
    let read = |value: &str, is: f64, unit: &str| {
        let d = Doc::load_str_in(&CORE, &format!("(def v {value})"), "t", Path::new("."), &[]).unwrap();
        let Some(Value::Field(n, Some(said))) = d.def("v").cloned() else { panic!("{value}") };
        let (q, times) = said.read().unwrap_or_else(|| panic!("{value}: no quantity"));
        let x = d.g.as_const(n).unwrap() * times;
        assert!(said.unit.quantity == Quantity::Other && (x - is).abs() <= 1e-9 * is.abs() && q.unit() == unit, "{value}: {x} {}", q.unit());
    };
    read("(* 1N*m 60rpm)", std::f64::consts::TAU, "W");
    read("(* 2mm 60rpm)", 2.0 * std::f64::consts::TAU, "mm/s");
    read("(* 2mm 90deg)", std::f64::consts::PI, "mm");
    read("(/ 1kW 60rpm)", 1e6 / std::f64::consts::TAU, "N*mm");
    read("(* (/ 1rad 1mm) (* 2mm2 1rad/s))", 2.0, "mm/s");
    // A slot reads it the same way; a turning rate and a frequency are never read as each other; a root that
    // halves an angle's power is no quantity; and `rad->deg` of an angle with its unit is its degrees, once.
    let slot = |limits: &str| {
        let src = format!("(body b (shape (- x 1))) (assembly root (instance a b (at) {limits}))");
        Doc::load_str_in(&CORE, &src, "t", Path::new("."), &[]).map(|d| {
            let bounds = d.assemblies[0].instances[0].joint.as_ref().unwrap().speed_limits.unwrap();
            (d.g.as_const(bounds.0).unwrap(), d.g.as_const(bounds.1).unwrap())
        })
    };
    let arc = slot("(slide z 0mm) (speed-limits (* -2mm 60rpm) (* (* 2mm 2) 30rpm))").unwrap();
    assert!((arc.0 + 2.0 * std::f64::consts::TAU).abs() < 1e-9 && (arc.1 - 2.0 * std::f64::consts::TAU).abs() < 1e-9, "{arc:?}");
    assert_eq!(slot("(turn z 0deg) (speed-limits (* (* 2mm 60rpm) (/ -1 2mm)) (* 2mm 60rpm (/ 1 2mm)))").unwrap(), (-360.0, 360.0));
    let no = slot("(turn z 0deg) (speed-limits 0Hz 1Hz)").unwrap_err();
    assert!(no.message.contains("a speed over a radius is radians a second, (* VALUE 1rad); a count of turns a second, (* VALUE 360deg)"), "{}", no.message);
    assert_eq!(slot("(turn z 0deg) (speed-limits 0deg/s (* 360deg 1Hz))").unwrap(), (0.0, 360.0));
    let no = slot("(turn z 0deg) (speed-limits 0deg/s (* 60rpm 90deg))").unwrap_err();
    assert!(no.message.contains("with another power of angle"), "{}", no.message);
    assert_eq!(number("(sqrt (/ (/ 1N*m 1rad) 1kg*m2))").1, Some("?"));
    assert_eq!([number("(rad->deg 1rad)"), number("(rad->deg 1)"), number("(deg->rad 180deg)"), number("(deg->rad 180)")].map(|n| n.1), [None; 4]);
    for (value, is) in
        [("(rad->deg 1rad)", 180.0 / std::f64::consts::PI), ("(rad->deg 1)", 180.0 / std::f64::consts::PI), ("(rad->deg 90deg)", 90.0), ("(deg->rad 180deg)", std::f64::consts::PI)]
    {
        assert!((number(value).0 - is).abs() < 1e-12, "{value}");
    }
    assert!((number("(deg->rad 180)").0 - std::f64::consts::PI).abs() < 1e-12 && (number("(deg->rad 2rad)").0 - 2.0).abs() < 1e-12);
    // What has no quantity keeps its powers: it is never a pure number to be given a unit again, and is a
    // quantity again when what is left is one. What truly cancels is a pure number.
    assert_eq!([number("(* (pow 2N 2) 1N)").1, number("(+ (* 2N 3s) 1)").1, number("(* (floor (* 2N 3mm)) 1N)").1, number("(* (exp 2N) 1mm)").1], [Some("?"); 4]);
    assert_eq!([number("(* (/ (* 2N 3s) 1s) 1mm)"), number("(/ (* 2N 3s) (* 2N 3s))"), number("(* 1m (/ 1 1mm))")], [(6.0, Some("N*mm")), (1.0, None), (1000.0, None)]);
    // A quantity scaled by pure numbers keeps the unit it was said in; its root and its square are what they are.
    assert_eq!([number("(* (/ 6 3) 5lbf)").1, number("(/ 5lbf (+ 1 1))").1, number("(sqrt 4mm2)").1, number("(sqrt 4mm)").1], [Some("lbf"), Some("lbf"), Some("mm"), Some("?")]);
    // A temperature is a point on its scale, in every scale: it is not multiplied. A difference is an amount.
    for point in ["300K", "20degC", "68degF"] {
        let e = Doc::load_str_in(&CORE, &format!("(def v (* 2 {point}))"), "t", Path::new("."), &[]).err().unwrap();
        assert!(e.message.starts_with(&format!("a temperature is a point on its scale, not an amount: {point} (t:1) cannot be multiplied")), "{}", e.message);
    }
    assert_eq!([number("(* 2 (- 300K 290K))").0, number("(* 2 (- 30degC 20degC))").0, number("(* 9 (- 100degF 50degF))").0], [20.0, 20.0, 250.0]);
    // A sum's terms are one quantity.
    let e = load("(+ 1g 1psi)").unwrap();
    assert!(e.starts_with("'+' takes values of one quantity, but 1g (t:1) is a mass and 1psi (t:1) is a pressure"), "{e}");
}

/// Every way of giving a temperature a unit that a refusal hints is one
/// that loads when written as the hint writes it: a temperature is a point
/// on its scale, never multiplied, so `(* 1degC e)` is no hint for one.
#[test]
fn a_temperature_hint_loads_as_written() {
    let load = |src: &str| Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "t", Path::new("."), &[]);
    let refused = |src: &str| load(src).err().unwrap_or_else(|| panic!("{src}: accepted"));
    // Each (refused source, the form the refusal says to write, where it goes in the source).
    let sum = "(param hot 20degC) (def k 3) (def v (+ hot k))";
    let e = refused(sum);
    for (form, replaces) in [("(+ 0degC k)", "k))"), ("(/ (- X 0degC) (- 1degC 0degC))", "hot k")] {
        assert!(e.message.contains(form), "{}", e.message);
        let followed = match replaces {
            "k))" => sum.replace("k))", &format!("{form}))")),
            _ => sum.replace("hot k", &format!("{} k", form.replace('X', "hot"))),
        };
        assert!(load(&followed).is_ok(), "{followed}: {:?}", load(&followed).err());
    }
    // A slot that takes a temperature, given a bare value worked out.
    let case = "(def k 3) (param hot 20degC (cases 20degC (* 2 k)))";
    let e = refused(case);
    assert!(e.message.ends_with("or of the value worked out: (+ 0degC (* 2 k))"), "{}", e.message);
    assert!(load(&case.replace("(* 2 k)", "(+ 0degC (* 2 k))")).is_ok());
    // A body's bare temperature read with a default said in degC.
    let e = refused("(body a (shape (box 1 1 1)) (attr T 5)) (def f (by-body a T 20degC))");
    let form = "(+ 0degC (by-body a T 20))";
    assert!(e.message.ends_with(form), "{}", e.message);
    assert!(load(&format!("(body a (shape (box 1 1 1)) (attr T 5)) (def f {form})")).is_ok());
    // A temperature scaled: each form its hint writes, in its own scale, loads and is what it says: 160 of
    // that scale (in degC, as ink holds a temperature), and 80 of its steps (in K, as it holds a difference).
    let at = |src: &str| {
        let d = load(&format!("(def v {src})")).unwrap_or_else(|e| panic!("{src}: {e:?}"));
        let v = d.defs.iter().find(|(n, _)| n == "v").map(|(_, v)| v.clone()).unwrap();
        let Value::Field(f, _) = v else { panic!("{src}") };
        d.load_eval(f)
    };
    for (scaled, u, point, step) in [
        ("(* 80 1degC)", "degC", 160.0, 1.0),
        ("(* 1degC k)", "degC", 160.0, 1.0),
        ("(* 2 68degF)", "degF", (160.0 - 32.0) * 5.0 / 9.0, 5.0 / 9.0),
        ("(/ 300K 2)", "K", 160.0 - 273.15, 1.0),
    ] {
        let e = refused(&format!("(def k 3) (def v {scaled})"));
        let hint = e.hint.unwrap_or_default();
        let worked = if u == "degF" { "(+ 0degF (* (* 80 2) (- 1degF 0degF)))".to_string() } else { format!("(+ 0{u} (* 80 2))") };
        for (form, is) in [(worked, point), (format!("(* 80 (- 1{u} 0{u}))"), 80.0 * step)] {
            assert!(hint.contains(&form), "{scaled}: {hint}");
            assert!((at(&form) - is).abs() < 1e-9, "{form}: {} is not {is}", at(&form));
        }
    }
    // Beside a temperature difference, an amount: what the hint adds is a difference too, so the sum is
    // one (13 K), not a temperature.
    let diff = "(def k 3) (def v (+ (- 30degC 20degC) k))";
    let e = refused(diff);
    for form in ["(* k (- 1degC 0degC))", "(/ X (- 1degC 0degC))"] {
        assert!(e.message.contains(form), "{}", e.message);
    }
    let followed = load(&diff.replace("20degC) k)", "20degC) (* k (- 1degC 0degC)))")).unwrap();
    assert!(followed.defs.iter().any(|(n, v)| n == "v" && matches!(v, Value::Field(_, Some(s)) if s.unit.quantity == crate::units::Quantity::Temperature && !s.is_point())));
    assert!((at("(+ (- 30degC 20degC) (* 3 (- 1degC 0degC)))") - 13.0).abs() < 1e-9);
    let e = refused("(def v (+ (- 30degC 20degC) 3))");
    assert!(e.message.ends_with("say it is a difference too: (* 3 (- 1degC 0degC))"), "{}", e.message);
}

/// A temperature in degF or K reads a point on its scale, so a change
/// (a tolerance, an accuracy) in one is refused in a sum too, whichever
/// term it is.
#[test]
fn a_sum_keeps_a_unit_with_a_zero_of_its_own() {
    for dev in ["(+ 1degC 33.8degF)", "(+ 33.8degF 1degC)", "(+ 1degC 1K)"] {
        let src = format!("(param temp 20degC (tol -1degC {dev}))");
        let e = Doc::load_str_in(&CORE, &src, "t", Path::new("."), &[]).err().unwrap_or_else(|| panic!("{dev}: accepted"));
        assert!(e.message.starts_with("a tolerance's deviation is a change"), "{dev}: {}", e.message);
    }
    // A deviation with a unit is its default's quantity; a bare term
    // beside a length is in mm, beside another quantity it is refused.
    let e = Doc::load_str_in(&CORE, "(param temp 20degC (tol -1degC (+ 1 1K)))", "t", Path::new("."), &[]).err().unwrap();
    assert!(e.message.starts_with("'+' takes values of one quantity, but 1K (t:1) is a temperature and 1 is a bare number"), "{}", e.message);
    let e = Doc::load_str_in(&CORE, "(param dim 1mm (tol -1g 1g))", "t", Path::new("."), &[]).err().unwrap();
    assert!(e.message.starts_with("a tolerance's deviation takes a length"), "{}", e.message);
    assert!(Doc::load_str_in(&CORE, "(param dim 1mm (tol (- 0.1 1um) (+ 0.1 1um)))", "t", Path::new("."), &[]).is_ok());
}

/// A bare number is a length (its square an area, its cube a volume), a
/// pure number or a pose's angle. Where another quantity is read it is
/// refused, naming the unit to write; nothing (0) and no limit stand.
#[test]
fn a_bare_number_where_another_quantity_is_read_is_refused_naming_its_unit() {
    let load = |src: &str| Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "t", Path::new("."), &[]).err().map(|e| e.message);
    let ball = |attr: &str| format!("(param m 3)\n(body b (shape (sphere 1)) {attr})");
    for (src, says) in [
        (ball("(attr mass 5)"), "attr mass takes a mass, and 5 is a bare number: say its unit (5g, if it is in g)"),
        (ball("(attr density 0.0027)"), "attr density takes a density, and 0.0027 is a bare number: say its unit (0.0027g/mm3, if it is in g/mm3)"),
        (ball("(attr mass (* 2 m))"), "attr mass takes a mass, and (* 2 m) is a bare number: say it where its number is written (5g), or of the value worked out: (* 1g (* 2 m))"),
        ("(t-span 0 2)".into(), "t-span takes a time, and 2 is a bare number: say its unit (2s, if it is in s)"),
        // A turn's limits are angles: a time is not one.
        (
            "(body b (shape (- x 1)))\n(assembly a (instance i b (at) (turn z 10deg) (limits -1s 1s)))".into(),
            "limits takes an angle (deg, or a number with its unit: deg rad), but s is a time, written -1s (t:3)",
        ),
        // A bare term beside a quantity: nothing says its unit.
        (
            "(t-span 0 1s)\n(def late (- t 300))".into(),
            "'-' takes values of one quantity, but a time worked out at t:3 is a time and 300 is a bare number: say its unit (300s, if it is in s)",
        ),
        (
            "(param p 2MPa)\n(def q (+ p (* 2 0.1)))".into(),
            "'+' takes values of one quantity, but 2MPa (t:2) is a pressure or stress and (* 2 0.1) is a bare number: it says no unit, so it is a length in mm or a pure number. If it is a pressure or stress, say so: (* 1MPa (* 2 0.1)); if the other is meant as a plain number, take its unit out: (/ X 1MPa)",
        ),
        // A value claim's two sides are one quantity.
        (
            "(param p 2MPa)\n(body b (shape (sphere 1)))\n(require r (value p < 5))".into(),
            "a value claim takes values of one quantity, but 2MPa (t:2) is a pressure or stress and 5 is a bare number: say its unit (5MPa, if it is in MPa)",
        ),
        (
            "(param p 2MPa)\n(body b (shape (sphere 1)))\n(require r (value p < 5N))".into(),
            "'value' takes values of one quantity, but 2MPa (t:2) is a pressure or stress and 5N (t:4) is a force",
        ),
        // A unit lost in a function is no number of another's.
        (
            ball("(attr mass (exp 1g))"),
            "attr mass takes a mass (g), but the value worked out at t:3 lost its unit in a function that takes a pure number: divide what goes in by its unit first ((/ x 1mm), (/ t 1s))",
        ),
        ("(param hot 20degC (tol -1 1degC))".into(), "param hot is written 20degC; give tol bound -1 a unit"),
        // A count is told it is one, not to become a force (3 N, whatever the force).
        (
            "(def f (vec 1N 2N 2N))\n(body b (shape (sphere 1)))\n(require r (value (len f) < 10N))".into(),
            "a value claim takes values of one quantity, but 10N (t:4) is a force and (len f) is a bare number: len counts a list's items, a pure number; a vector's size, which keeps its unit, is std's (norm V), or (len3 X Y Z)",
        ),
    ] {
        assert_eq!(load(&src).as_deref(), Some(says), "{src}");
    }
    // A value of a known dimension that is not the slot's is no number of its unit either: a product with a
    // unit, made of a formula that had already mixed a bare number with a united one, is refused, not rescaled.
    let e = load(&ball("(attr mass (* 1g (/ m 2s)))")).unwrap();
    assert!(e.starts_with("attr mass takes a mass (g), but the value worked out at t:3 is g/s, not a mass"), "{e}");
    for src in [
        ball("(attr mass 5g)"),
        ball("(attr mass 0)"),
        // A zero written with its unit says the unit of what is added to it.
        "(def r 5)\n(def rate (+ 0deg/s r))\n(body b (shape (sphere 1)) (attr mass (+ 0g (* 2 r))))".into(),
        "(t-span 0 1s)\n(def late (max 0 (- t 0.3s)))".into(),
        ball("(attr mass (* 1g 2 m))"),
        ball("") + "(require small (volume (shape-of b) < 5))",
        ball("") + "(require small (volume (shape-of b) < (* m m m)))",
        "(t-span 0 2s)".into(),
        "(body b (shape (rotate-z (box 1 2 3) 30)))".into(),
    ] {
        assert_eq!(load(&src), None, "{src}");
    }
}

/// by-body's material numbers are bare: what is refused on the way to
/// a conductivity says the one form that works, `(* 1W/m/K (by-body
/// root k 0.2))` (SIM.md), and it loads.
#[test]
fn by_body_says_the_form_that_works() {
    let load = |k: &str| {
        let src = format!(
            "(import \"std\")\n(body a (shape (box 2 2 2)) (attr material \"petg\"))\n(body c (shape (move (box 2 2 2) 4 0 0)))\n(assembly root (instance a a) (instance c c))\n(def k {k})"
        );
        Doc::load_str_in(&CORE, &src, "t", Path::new("."), &[]).err().map(|e| e.message)
    };
    let e = load("(by-body root k)").unwrap();
    assert!(e.ends_with("has no attr 'k', nor a material that gives one: give it one, or by-body a default for such bodies, (by-body root k DEFAULT)"), "{e}");
    let e = load("(by-body root k 0.2W/m/K)").unwrap();
    assert!(e.contains("so its default is bare too, and the field takes its unit whole: (* 1W/m/K (by-body root k 0.2))"), "{e}");
    assert_eq!(load("(* 1W/m/K (by-body root k 0.2))"), None);
}

/// A `for` makes numbers of its list's values: each keeps the unit its
/// value was said in.
#[test]
fn a_for_keeps_its_values_units() {
    let load = |list: &str| Doc::load_str_in(&CORE, &format!("(for q {list} (body b (shape (- x 1)) (attr mass q)))"), "t", Path::new("."), &[]).err().map(|e| e.message);
    for list in ["(vec 50psi)", "(vec (+ 50psi 1psi))", "(let ((p 50psi)) (vec p))"] {
        let e = load(list).unwrap_or_else(|| panic!("{list}: accepted"));
        assert!(e.starts_with("attr mass takes a mass") && e.contains("written 50psi (t:1)"), "{list}: {e}");
    }
    assert_eq!(load("(vec 5g)"), None);
}

/// A --param value in a unit is checked as its default would be:
/// against the default's quantity, and where the param is used.
#[test]
fn a_param_given_in_a_unit_is_said_in_it() {
    let dir = Scratch::new("param-units");
    let load = |default: &str, given: &str| {
        let f = dir.join("m.ink");
        std::fs::write(&f, format!("(param m {default})\n(body b (shape (- x 1)) (attr mass m))")).unwrap();
        let a = crate::command::Args::parse(&["eval".into(), f.display().to_string(), "--param".into(), format!("m={given}")]);
        crate::command::param_args(&a).and_then(|p| Doc::load_file_in(&CORE, &f, &p.values)).err().map(|e| e.message)
    };
    let says = |default: &str, given: &str, text: &str| {
        let e = load(default, given).unwrap_or_else(|| panic!("{default} {given}: accepted"));
        assert!(e.contains(text), "{default} {given}: {e}");
    };
    says("5g", "50psi", "--param m: m is a mass, written 5g (");
    says("5g", "50psi", "m.ink:1), but psi is a pressure");
    says("5", "50psi", "attr mass takes a mass");
    says("5", "50psi", "written 50psi (--param m)");
    says("5", "1psi:2psi", "attr mass takes a mass");
    says("5", "1g:2psi", "--param m: a range's ends are one quantity, but g is a mass and psi is a pressure");
    // A mass is given with its unit, as it is written with one.
    says("5g", "50", "param m is written 5g; give m=50 a unit");
    says("5g", "1mg:2", "give both ends of the range a unit");
    for given in ["50g", "1g:2g", "1mg:2g"] {
        assert_eq!(load("5g", given), None, "{given}");
    }
}

/// A value given a param keeps its unit wherever it travels: to a
/// range by the name its form gives it, through a reload that changes
/// the value (an optimizer's design), and across a range's ends (degF
/// over degC, as in a sum, whichever end it is).
#[test]
fn a_given_unit_travels_with_its_value() {
    let dir = Scratch::new("given-travels");
    let run = |name: &str, src: &str, given: &str| {
        let f = dir.join(name);
        std::fs::write(&f, src).unwrap();
        let a = crate::command::Args::parse(&["eval".into(), f.display().to_string(), "--param".into(), given.into()]);
        crate::command::param_args(&a).and_then(|p| crate::view::ModelView::load_in(&CORE, &f, &p.values, 0.0))
    };
    // An unnamed range, by its form's name.
    let size = "(body b (shape (- (max (abs x) (abs y) (abs z)) (ranged 1mm 5mm))))";
    let e = run("size.ink", size, "b 0=500psi").err().expect("a pressure as a length");
    assert!(e.message.contains("--param b 0: b 0 is a length, written 1mm") && e.message.contains("but psi is a pressure"), "{}", e.message);
    assert!(run("size.ink", size, "b 0=0.1in").is_ok());
    // A reload changing the value keeps its unit: m is a pressure.
    let opt = "(param m 0)\n(def v (if (> m 1) m 2mm3))\n(body b (shape (- (max (abs x) (abs y) (abs z)) 0.5)))\n(require a (volume (shape-of b) < v))";
    let v = run("opt.ink", opt, "m=50psi").unwrap_or_else(|e| panic!("{e:?}"));
    let e = v.reload(&[("m".into(), 3.0.into())]).err().expect("a pressure as a volume");
    assert!(e.message.starts_with("volume takes a volume") && e.message.contains("(--param m)"), "{}", e.message);
    assert!(run("opt.ink", opt, "m=50").unwrap().reload(&[("m".into(), 3.0.into())]).is_ok());
    // A range's ends: degF or K is kept over degC, either way round.
    let eps = "(param eps 1) (for v (vec eps) (require a (value 1 > 0) (accuracy v)))";
    for given in ["eps=1degC:274.15K", "eps=274.15K:1degC"] {
        let e = run("eps.ink", eps, given).err().unwrap_or_else(|| panic!("{given}: accepted"));
        assert!(e.message.starts_with("an accuracy is a change, and K reads a point on its scale"), "{given}: {}", e.message);
    }
}
