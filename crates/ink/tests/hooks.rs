//! The extension points no package of the kernel's uses (a package built
//! on it may): a package's values in expressions and its claims made once
//! the file is read (Package::values, value, finish), and a provider's own
//! steps over a motion, its band about a round reading, its error band,
//! its note, what it rests on and its enclosure over a span
//! (ParamProvider::steps, band, error, note, reads, over_span); and an
//! assembly clause a package owns (Package::clauses, clause). A test
//! package, `gauge`, drives each as a third party's would, through the
//! public API; `ink check`'s (Package::check) is ink-cli's test.

use ink::InkError;
use ink::claims::{ClaimOptions, Status, run_asserts};
use ink::graph::Var;
use ink::interval::I;
use ink::lang::sexp::Sexp;
use ink::lang::{AssertDef, Claim, Cmp, Doc, Form, Package, ParamProvider, Reader, Registry, Steps, Value};
use ink::view::ModelView;
use std::path::Path;
use std::sync::Arc;

/// `(gauge MODE [FIELD])`: a number this package provides at each moment,
/// t itself, read off a gauge; MODE says which of the provider's hooks it
/// gives (`plain`: none but its note); FIELD, a field it is read from.
struct Gauge;

const SPAN_WIDTH: f64 = 0.25;

impl Package for Gauge {
    fn name(&self) -> &'static str {
        "gauge"
    }
    fn forms(&self) -> &'static [&'static str] {
        &[]
    }
    fn read(&self, r: &mut Reader, form: &Sexp, head: &str, _: &[Sexp]) -> Result<Arc<dyn Form>, InkError> {
        Err(r.err(form, format!("gauge reads no form '{head}'")))
    }
    fn values(&self) -> &'static [&'static str] {
        &["gauge"]
    }
    fn value(&self, r: &mut Reader, e: &Sexp, _head: &str, args: &[Sexp]) -> Result<Value, InkError> {
        let mode = r.name(args.first().ok_or_else(|| r.err(e, "(gauge MODE [FIELD])"))?)?;
        if let Some(f) = args.get(1) {
            r.field(f)?;
        }
        let k = r.derived(&format!("gauge:{mode}"));
        r.solved(k);
        r.keep("gauge", "gauge", e, Arc::new(Reading { mode, owns: [k] }))?;
        Ok(Value::Field(r.graph().var(Var::Param(k as u32))))
    }
    fn clauses(&self) -> &'static [&'static str] {
        &["mark"]
    }
    /// `(mark NAME PATH SIZE)` in an assembly: a mark on what its instance
    /// PATH places, of a size read in the assembly's scope.
    fn clause(&self, r: &mut Reader, c: &Sexp, _head: &str, args: &[Sexp]) -> Result<Arc<dyn std::any::Any + Send + Sync>, InkError> {
        let [name, path, size] = args else { return Err(r.err(c, "(mark NAME PATH SIZE)")) };
        let name = r.name(name)?;
        if r.kept("mark").iter().filter_map(|k| k.downcast_ref::<Mark>()).any(|m| m.name == name) {
            return Err(r.err(c, format!("two marks named '{name}'")));
        }
        Ok(Arc::new(Mark { name, path: r.instance_path(path)?, size: r.field(size)? }))
    }
    /// A claim of its own on each `(gauge finish)` read: under 2 at every
    /// moment.
    fn finish(&self, r: &mut Reader) -> Result<(), InkError> {
        let read: Vec<usize> = r.forms::<Reading>().iter().filter(|g| g.mode == "finish").map(|g| g.owns[0]).collect();
        for k in read {
            let f = r.graph().var(Var::Param(k as u32));
            let at = ink::lang::sexp::parse("(gauge finish)").unwrap().remove(0);
            let def = AssertDef {
                name: "gauge:under-2".into(),
                claim: Claim::Value(f, Cmp::Lt, 2.0),
                text: "a gauge's reading is under 2".into(),
                bound: None,
                res: None,
                accuracy: None,
                quantity: None,
                file: "gauge".into(),
                line: 0,
                terms: Vec::new(),
                fixed: Vec::new(),
            };
            r.assert(&at, def)?;
        }
        Ok(())
    }
}

struct Mark {
    name: String,
    path: String,
    size: ink::graph::NodeId,
}

#[derive(Clone)]
struct Reading {
    mode: String,
    owns: [usize; 1],
}

impl Form for Reading {
    fn name(&self) -> &str {
        &self.mode
    }
    fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(self.clone()))
    }
}

impl ParamProvider for Reading {
    fn name(&self) -> &str {
        &self.mode
    }
    fn owns(&self) -> &[usize] {
        &self.owns
    }
    fn uses_t(&self) -> bool {
        true
    }
    fn on_demand(&self) -> bool {
        self.mode.starts_with("lazy")
    }
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, String> {
        params[self.owns[0]] = t;
        Ok(None)
    }
    fn reads(&self) -> Option<Vec<usize>> {
        (self.mode == "lazy-reads-nothing").then(Vec::new)
    }
    fn steps(&self, _: &[f64], (a, b): (f64, f64)) -> Option<Result<Steps, String>> {
        // A band is read at the steps its provider chooses (motion::measured_at).
        matches!(self.mode.as_str(), "steps" | "band")
            .then(|| Ok(Steps { ends: vec![a, 0.5 * (a + b), b], unchecked: Vec::new(), loose: Vec::new(), basis: "the gauge's own steps".into() }))
    }
    fn band(&self, _: &[f64], t: f64) -> Option<(Vec<f64>, Vec<f64>)> {
        (self.mode == "band").then(|| (vec![SPAN_WIDTH], vec![t]))
    }
    fn error(&self, _: &[f64]) -> Option<Vec<f64>> {
        (self.mode == "error").then(|| vec![SPAN_WIDTH])
    }
    fn note(&self) -> Option<String> {
        Some("read off a gauge".into())
    }
    fn over_span(&self, _: &[f64], (a, b): (f64, f64)) -> Option<Result<(Vec<f64>, Vec<I>, Vec<I>), String>> {
        (self.mode == "span").then(|| Ok((vec![0.5 * (a + b)], vec![I::new(a, b)], vec![I::point(1.0)])))
    }
}

// With design, for a solve on t reading the gauge.
static GAUGE: Registry = Registry::new(&[&Gauge, &ink::design::Design]);

fn view(src: &str) -> ModelView {
    let d = Doc::load_str_in(&GAUGE, &format!("(import \"std\")\n{src}"), "gauge-test", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let span = d.t_span;
    let (ctx, unsolved) = d.ctx_at(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    v.unsolved = unsolved;
    v.set_motion(span).unwrap();
    v
}

/// Each claim's status, by name.
fn statuses(src: &str) -> Vec<(String, Status)> {
    let mut v = view(src);
    run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts.into_iter().map(|a| (a.name, a.check.status)).collect()
}

fn status(src: &str, name: &str) -> Status {
    statuses(src).into_iter().find(|(n, _)| n == name).unwrap_or_else(|| panic!("no claim {name}")).1
}

#[test]
fn a_package_value_is_a_number_claims_read_and_its_claims_join_the_files() {
    // (gauge ...) is t at each moment; a claim the package makes once the
    // file is read (finish) is checked with the file's own.
    let got = statuses("(t-span 0 1)\n(def g (gauge finish))\n(assert late (value g > 0.5))");
    let late = got.iter().find(|(n, _)| n == "late").map(|c| c.1);
    let own = got.iter().find(|(n, _)| n == "gauge:under-2").map(|c| c.1);
    assert_eq!((late, own), (Some(Status::Fails), Some(Status::Holds)), "{got:?}");
    let d = Doc::load_str_in(&GAUGE, "(def g (gauge finish))", "t", Path::new("."), &[]).unwrap();
    assert!(d.asserts.iter().any(|a| a.name == "gauge:under-2"));
    // Without the package, the value is no function.
    static CORE: Registry = Registry::new(&[]);
    let e = Doc::load_str_in(&CORE, "(def g (gauge finish))", "t", Path::new("."), &[]).err().expect("unknown");
    assert!(e.message.contains("unknown function 'gauge'"), "{e:?}");
}

#[test]
fn a_claim_over_a_motion_is_measured_at_the_steps_its_provider_chooses() {
    // steps: the gauge's two steps (read at their middles and the span's
    // ends: 4 readings), and the answer says whose steps they are; without
    // them, the uniform steps of a measured claim.
    let mut basis = Vec::new();
    for mode in ["steps", "plain"] {
        let mut v = view(&format!("(t-span 0 1)\n(def g (gauge {mode}))\n(assert under (value g < 2))"));
        let a = run_asserts(&mut v, Some("under"), false, ClaimOptions::default()).unwrap().asserts.remove(0);
        assert_eq!(a.check.status, Status::Holds, "{mode}");
        let m = a.check.motion.expect("over its moments");
        basis.push((m.basis, m.motion.steps));
    }
    assert_eq!(basis[0], (Some("the gauge's own steps".to_string()), 4), "{basis:?}");
    assert!(basis[1].0.is_none() && basis[1].1 > 4, "{basis:?}");
}

#[test]
fn a_band_about_a_round_reading_holds_only_where_all_of_it_holds() {
    // band: at each step the reading is t within 0.25. Under 1.1 the
    // real values (t up to 1) hold, but the band reaches 1.25: undecided.
    // Under 2 the whole band holds; over 0.5 fails at the real values.
    let src = |mode: &str, claim: &str| format!("(t-span 0 1)\n(def g (gauge {mode}))\n(assert c (value g {claim}))");
    assert_eq!(status(&src("band", "< 1.1"), "c"), Status::Unknown);
    assert_eq!(status(&src("plain", "< 1.1"), "c"), Status::Holds);
    assert_eq!(status(&src("band", "< 2"), "c"), Status::Holds);
    assert_eq!(status(&src("band", "> 0.5"), "c"), Status::Fails);
}

#[test]
fn an_estimate_with_an_error_band_is_decided_only_over_all_of_it() {
    // error: at t = 0 the gauge reads 0 within 0.25 either way.
    let src = |mode: &str, claim: &str| format!("(def g (gauge {mode}))\n(assert c (value g {claim}))");
    assert_eq!(status(&src("error", "< 0.1"), "c"), Status::Unknown);
    assert_eq!(status(&src("plain", "< 0.1"), "c"), Status::Holds);
    assert_eq!(status(&src("error", "< 1"), "c"), Status::Holds);
    assert_eq!(status(&src("error", "> 1"), "c"), Status::Fails);
}

#[test]
fn an_answer_reading_a_provided_number_carries_its_note() {
    let mut v = view("(def g (gauge plain))\n(assert c (value g < 1))");
    let a = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts.remove(0);
    assert!(a.check.inexact.notes.contains(&"read off a gauge".to_string()), "{:?}", a.check.inexact.notes);
}

#[test]
fn what_a_provider_rests_on_is_its_own_account_where_it_gives_one() {
    // reads: provided on demand at one value of each param, a gauge read
    // from a ranged w is refused over w's range, unless it says it rests
    // on nothing of it.
    for (mode, refused) in [("lazy", true), ("lazy-reads-nothing", false)] {
        let mut v = view(&format!("(param w 5 (tol -1 1))\n(def g (gauge {mode} w))"));
        let spans = v.declared_spans().unwrap();
        v.set_param_spans(&spans).unwrap();
        let k = v.doc.params.iter().position(|p| p.name == format!("gauge:{mode}")).unwrap();
        let g = v.doc.g.var(Var::Param(k as u32));
        let e = v.want(&[g]).err();
        assert_eq!(e.is_some(), refused, "{mode}: {e:?}");
        if let Some(e) = e {
            assert!(e.message.contains("reads ranged 'w'"), "{e:?}");
        }
    }
}

#[test]
fn a_solve_on_t_reading_a_provider_without_equations_is_followed_over_its_enclosure() {
    // over_span: the gauge encloses itself over a span (t's own range, its
    // rate 1), so a solve on t reading it is followed along t; one that
    // gives none is not.
    let src = |mode: &str| format!("(t-span 0 1)\n(def g (gauge {mode}))\n(solve fit (unknowns (a 0.5)) (== (* 2 a) g))\n(assert c (value a < 0.75))");
    let (span, plain) = (status(&src("span"), "c"), status(&src("plain"), "c"));
    assert_eq!(span, Status::Holds);
    assert_ne!(plain, Status::Holds, "followed without an enclosure");
}

#[test]
fn an_assembly_clause_a_package_owns_is_read_in_its_scope_and_placed_with_each_instance() {
    // (mark ...) in an assembly given arguments: its size reads the
    // argument, its path is checked against the instances read so far, and
    // each instance of the assembly places it: under root, twice, each with
    // its instance's path and placement.
    let src = "(import \"std\") (body b (shape (sphere 1)))
(assembly (pair (s 2)) (instance a b) (mark m a (* 3 s)))
(assembly root (instance l (pair 4) (at (move -10 0 0))) (instance r (pair 5) (at (move 10 0 0))))";
    let mut d = Doc::load_str_in(&GAUGE, src, "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let placed = d.clauses_under("root", "gauge");
    let got: Vec<(String, f64, f64)> = placed
        .iter()
        .map(|p| {
            let m = p.clause.kept.downcast_ref::<Mark>().expect("a mark");
            let at = ink::eval::eval_point(&ink::tape::Tape::new(&d.g, p.at.t[0]), &ctx, [0.0; 3]);
            (format!("{}/{}:{}", p.path, m.path, m.name), at, ink::eval::eval_point(&ink::tape::Tape::new(&d.g, m.size), &ctx, [0.0; 3]))
        })
        .collect();
    assert_eq!(got, [("root/l/a:m".to_string(), -10.0, 12.0), ("root/r/a:m".to_string(), 10.0, 15.0)]);
    assert!(d.clauses_under("root", "other").is_empty());
    // Its path is an instance read so far; a name once per assembly; and
    // without the package, no clause of an assembly.
    for (bad, says) in [("(assembly q (mark m nothing 1))", "no instance 'nothing'"), ("(assembly q (instance a b) (mark m a 1) (mark m a 2))", "two marks named 'm'")] {
        let e = Doc::load_str_in(&GAUGE, &format!("(import \"std\") (body b (shape (sphere 1)))\n{bad}"), "t", Path::new("."), &[]).err().expect("refused");
        assert!(e.message.contains(says), "{bad}: {e:?}");
    }
    static CORE: Registry = Registry::new(&[]);
    let e = Doc::load_str_in(&CORE, src, "t", Path::new("."), &[]).err().expect("unknown");
    assert!(e.message.contains("assembly clauses are") && !e.message.contains("(mark"), "{e:?}");
}
