//! The extension points no package of the kernel's uses (a package built
//! on it may): a package's values in expressions and its claims made once
//! the file is read (Package::values, value, finish), and a provider's own
//! steps over a motion, its band about a round reading, its error band,
//! its note, what it rests on and its enclosure over a span
//! (ParamProvider::steps, band, error, note, reads); and an
//! assembly clause a package owns (Package::clauses, clause). A test
//! package, `gauge`, drives each as a third party's would, through the
//! public API; `ink check`'s (Package::check) is ink-cli's test.

use ink::InkError;
use ink::claims::{ClaimOptions, Status, run_asserts};
use ink::graph::Var;
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
        // twice: twice the param FIELD is, read when it is solved.
        let mut from = None;
        if let Some(f) = args.get(1) {
            let n = r.field(f)?;
            if let ink::graph::Node::Var(Var::Param(k)) = r.graph().nodes[n as usize] {
                from = Some(k as usize);
            }
        }
        let k = r.derived(&format!("gauge:{mode}"));
        r.solved(k);
        r.keep("gauge", "gauge", e, Arc::new(Reading { mode, owns: [k], from }))?;
        Ok(Value::Field(r.graph().var(Var::Param(k as u32)), None))
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
                scope: r.scope(),
                line: 0,
                terms: Vec::new(),
                fixed: Vec::new(),
            };
            r.assert(&at, def)?;
        }
        Ok(())
    }
    fn metrics(&self) -> &'static [(&'static str, &'static str)] {
        &[("height", "how tall a shape is, top to bottom (mm)")]
    }
    fn shape_metrics(&self) -> &'static [&'static str] {
        &["height"]
    }
    /// `(height f)`: f's proven box's height, not the same of a prism in
    /// any slice (it says no plane_metric).
    fn shape_metric(&self, v: &mut ModelView, _metric: &str, f: u32, _res: Option<f64>, _accuracy: f64) -> Result<ink::lang::Measured, InkError> {
        let value = v.bounds(f).map_or(Ok(0.0), |(lo, hi, _)| Ok(hi[2] - lo[2]));
        Ok(ink::lang::Measured {
            value,
            at: None,
            res: None,
            finest: None,
            warnings: Vec::new(),
            excluded: Vec::new(),
            unseen: None,
            unseen_finer: true,
            resolves: 0.0,
            unbounded: None,
            extent: None,
            read_error: None,
            proven: true,
            held: None,
            band: None,
        })
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
    from: Option<usize>,
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
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, ink::lang::Refused> {
        // far: a billion, known to half a unit either way.
        params[self.owns[0]] = match (self.mode.as_str(), self.from) {
            ("far", _) => 1e9,
            ("twice", Some(k)) => 2.0 * params[k],
            _ => t,
        };
        Ok(None)
    }
    fn reads(&self) -> Option<Vec<usize>> {
        (self.mode == "lazy-reads-nothing").then(Vec::new)
    }
    fn steps(&self, _: &[f64], (a, b): (f64, f64)) -> Option<Result<Steps, ink::lang::Refused>> {
        // A band is read at the steps its provider chooses (motion::measured_at).
        matches!(self.mode.as_str(), "steps" | "band")
            .then(|| Ok(Steps { ends: vec![a, 0.5 * (a + b), b], unchecked: Vec::new(), loose: Vec::new(), basis: "the gauge's own steps".into() }))
    }
    fn band(&self, _: &[f64], t: f64) -> Option<(Vec<f64>, Vec<f64>)> {
        (self.mode == "band").then(|| (vec![SPAN_WIDTH], vec![t]))
    }
    fn error(&self, _: &[f64], _: f64) -> Option<Vec<f64>> {
        match self.mode.as_str() {
            "error" => Some(vec![SPAN_WIDTH]),
            "far" => Some(vec![0.5]),
            _ => None,
        }
    }
    fn note(&self) -> Option<String> {
        Some("read off a gauge".into())
    }
}

static GAUGE: Registry = Registry::new(&[&Gauge, &ink::Shape]);

fn view(src: &str) -> ModelView {
    view_in(src, Path::new("."))
}

/// `view`, its files read from `dir`.
fn view_in(src: &str, dir: &Path) -> ModelView {
    let d = Doc::load_str_in(&GAUGE, &format!("(import \"std\")\n{src}"), "gauge-test", dir, &[]).unwrap_or_else(|e| panic!("{e:?}"));
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

/// A shape metric is asked in the plane of a 2D profile only where its
/// package says the plane keeps its meaning (Package::plane_metric): a
/// prism's height, measured in the plane's band, was a guaranteed 200 mm.
#[test]
fn a_metric_is_asked_in_the_plane_only_where_its_package_says_so() {
    assert_eq!(status("(require tall (height (circle-2d 1) < 300))", "tall"), Status::Unknown);
    assert_eq!(status("(require walls (min-thickness (circle-2d 1) > 1.9))", "walls"), Status::Holds);
}

#[test]
fn a_package_value_is_a_number_claims_read_and_its_claims_join_the_files() {
    // (gauge ...) is t at each moment; a claim the package makes once the
    // file is read (finish) is checked with the file's own. A gauge that
    // chooses no steps says nothing of how it changes between the moments
    // read, so a claim that holds at each is not settled over the span.
    let got = statuses("(t-span 0 1s)\n(def g (gauge finish))\n(require late (value g > 0.5))");
    let late = got.iter().find(|(n, _)| n == "late").map(|c| c.1);
    let own = got.iter().find(|(n, _)| n == "gauge:under-2").map(|c| c.1);
    assert_eq!((late, own), (Some(Status::Fails), Some(Status::Unknown)), "{got:?}");
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
    // them, the uniform steps of a measured claim, and no yes: nothing says
    // what the gauge reads between them.
    let mut basis = Vec::new();
    for (mode, status) in [("steps", Status::Holds), ("plain", Status::Unknown)] {
        let mut v = view(&format!("(t-span 0 1s)\n(def g (gauge {mode}))\n(require under (value g < 2))"));
        let a = run_asserts(&mut v, Some("under"), false, ClaimOptions::default()).unwrap().asserts.remove(0);
        assert_eq!(a.check.status, status, "{mode}");
        assert_eq!(a.check.motion.as_ref().is_some_and(|m| m.between.is_some()), mode == "plain", "{mode}");
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
    let src = |mode: &str, claim: &str| format!("(t-span 0 1s)\n(def g (gauge {mode}))\n(require c (value g {claim}))");
    assert_eq!(status(&src("band", "< 1.1"), "c"), Status::Unknown);
    assert_eq!(status(&src("steps", "< 1.1"), "c"), Status::Holds);
    assert_eq!(status(&src("band", "< 2"), "c"), Status::Holds);
    assert_eq!(status(&src("band", "> 0.5"), "c"), Status::Fails);
}

#[test]
fn an_estimate_with_an_error_band_is_decided_only_over_all_of_it() {
    // error: at t = 0 the gauge reads 0 within 0.25 either way.
    let src = |mode: &str, claim: &str| format!("(def g (gauge {mode}))\n(require c (value g {claim}))");
    assert_eq!(status(&src("error", "< 0.1"), "c"), Status::Unknown);
    assert_eq!(status(&src("plain", "< 0.1"), "c"), Status::Holds);
    assert_eq!(status(&src("error", "< 1"), "c"), Status::Holds);
    assert_eq!(status(&src("error", "> 1"), "c"), Status::Fails);
    // Read directly and through a second solve, its band reaches the claim through that solve too, where
    // no bound is carried: g + h is 3g, up to 0.75 over g's band, not the 0.25 g alone moves it by.
    let through = "(def g (gauge error))\n(def h (gauge twice g))\n(require c (value (+ g h) < 0.5))";
    assert_eq!(status(through, "c"), Status::Unknown);
}

/// An estimate's guard (estimate-when) is a formula, not a solve: a band its condition reads is carried
/// through the value by interval arithmetic where the condition keeps its sign over all of the band (1 + 2g,
/// g 0 within 0.25, is 0.5 to 1.5), and has no bound where it may change sign there (g < 0.1 holds at g's
/// value, not over its band). A second solve still carries none (#205's case, above).
#[test]
fn a_band_through_an_estimates_condition_is_carried_where_the_condition_holds_over_it() {
    let src = |cond: &str, claim: &str| format!("(def g (gauge error))\n(require c (value (estimate-when (+ 1 (* 2 g)) {cond} \"g is small\") {claim}))");
    assert_eq!(status(&src("(< g 1)", "< 1.6"), "c"), Status::Holds);
    assert_eq!(status(&src("(< g 1)", "< 1.4"), "c"), Status::Unknown);
    assert_eq!(status(&src("(< g 1)", "> 1.6"), "c"), Status::Fails);
    assert_eq!(status(&src("(< g 0.1)", "< 1.6"), "c"), Status::Unknown);
    assert_eq!(status(&src("(< g 0.1)", "< 1.1"), "c"), Status::Unknown);
}

/// An estimate that does not apply leaves a claim reading it unknown, however its value is written: one
/// cancelled out of the value (est less est, 0 times est) is still read, the value the same.
#[test]
fn an_estimate_cancelled_out_of_a_value_still_marks_it() {
    let src = |value: &str| format!("(def est (estimate-when 5 0 \"never applies\"))\n(require c (value {value} > 0))");
    for value in ["est", "(+ 1 (- est est))", "(+ 1 (* 0 est))", "(+ 1 (* est 0))", "(- (+ 1 est) est)", "(+ 1 (* 0 (* 2 est)))", "(select est 1 1)"] {
        // (A corner-reach of it reads x, which no value claim takes: ink-core's
        // every_constructor_preserves_its_input_reads keeps its marker.)
        assert_eq!(status(&src(value), "c"), Status::Unknown, "{value}");
    }
    // Where it applies, what cancels it is what it was: 1, proven.
    let applies = |value: &str| format!("(def est (estimate-when 5 1 \"applies\"))\n(require c (value {value} > 0))");
    for value in ["(+ 1 (- est est))", "(+ 1 (* 0 est))", "(- 1 est)"] {
        assert_eq!(status(&applies(value), "c"), if value == "(- 1 est)" { Status::Fails } else { Status::Holds }, "{value}");
    }
}

/// A guarded radial-load ratio whose radial force is 100 N within 0.25 N (guarded on there being no
/// axial force): 0.09975 to 0.10025 of its 1000 N rating, so under 0.1003 it holds, over it fails,
/// and under 0.1001 it is not settled. An axial force within that band (0 within 0.25 N) may cross
/// the guard: not settled, whatever the limit.
#[test]
fn a_guarded_ratio_on_a_banded_load_is_decided_over_all_of_the_band() {
    let at = |fr: &str, fa: &str, claim: &str| {
        let src = format!(
            "(def g (gauge error))\n(def axial (gauge plain))\n\
             (fn guarded-ratio (fx fy fz)\n\
               (let ((fr (+ 0N (len2 fx fy))) (fa (+ 0N (abs fz))))\n\
                 (/ (estimate-when fr (select (neg (== fa 0)) 1 (<= fa (* 1e-12 (max (abs fx) (abs fy))))) \"axial force\") 1000N)))\n\
             (require c (value (guarded-ratio {fr} 0N {fa}) {claim}))"
        );
        let mut v = view(&src);
        run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts.remove(0).check.status
    };
    assert_eq!(at("(* (+ 100 g) 1N)", "(* axial 1N)", "< 0.1003"), Status::Holds);
    assert_eq!(at("(* (+ 100 g) 1N)", "(* axial 1N)", "> 0.1003"), Status::Fails);
    assert_eq!(at("(* (+ 100 g) 1N)", "(* axial 1N)", "< 0.1001"), Status::Unknown);
    assert_eq!(at("100N", "(* g 1N)", "< 0.5"), Status::Unknown);
}

#[test]
fn an_answer_reading_a_provided_number_carries_its_note() {
    let mut v = view("(def g (gauge plain))\n(require c (value g < 1))");
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

#[test]
fn package_clause_takes_the_span_where_the_assembly_arguments_were_passed() {
    let dir = std::env::temp_dir().join(format!("ink-mark-span-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("part.ink"), "(t-span 0 1s) (body pin (shape 1)) (assembly (pair (s 1)) (instance p pin) (mark m p s))").unwrap();
    let src = r#"(import "part.ink" as p)
(t-span 0 10s)
(assembly root (instance d (p.pair 1)))"#;
    let mut d = Doc::load_str_in(&GAUGE, src, "caller", &dir, &[]).unwrap();
    std::fs::remove_dir_all(dir).unwrap();
    let own = d.clauses_under("root", "gauge").remove(0).clause;
    let imported = d.clauses_under("p.pair", "gauge").remove(0).clause;
    assert_eq!(d.span_of(own.scope), Some((0.0, 10.0)));
    assert_eq!((&*own.file, own.line), ("caller", 3));
    assert_eq!(d.span_of(imported.scope), Some((0.0, 1.0)));
}

#[test]
fn what_a_claim_reads_of_t_besides_its_providers_steps_is_shown_over_each() {
    // The gauge's steps follow the gauge (t, read at 0, 0.25, 0.75 and 1),
    // not what the claim adds to it. A swing of 64 turns in the span is
    // zero at each of those moments and 5 between them: not settled, never
    // the yes its moments give. A slow one is bounded over each step, and
    // the claim holds.
    let src = |swing: &str| format!("(t-span 0 1s)\n(def g (gauge steps))\n(require c (value (+ g {swing}) < 3))");
    let mut v = view(&src("(* 5 (sin (* 23040deg/s t)))"));
    let a = run_asserts(&mut v, Some("c"), false, ClaimOptions::default()).unwrap().asserts.remove(0);
    assert_eq!(a.check.status, Status::Unknown);
    let said = a.check.motion.as_ref().and_then(|m| m.unsettled()).unwrap_or_default();
    assert!(said.contains("it reads t itself"), "{said}");
    assert_eq!(status(&src("(* 0.1 (sin (* 360deg/s t)))"), "c"), Status::Holds);
}

#[test]
fn a_band_small_beside_its_value_still_moves_a_claim_about_shapes() {
    // The gauge reads a billion within half a unit: a two-billionth of it,
    // and half a millimetre of where the 3 mm box it places sits. Beside a
    // box 0.2 away the claim is not settled (the band reaches it); 2 away it
    // holds at every value in the band; with the other box inside it, it
    // fails at every one. No band is dropped for being small beside the
    // number it is on.
    let src = |at: &str| format!("(def g (gauge far))\n(require clear (disjoint (move (box 3 1 1) (- g 1e9) 0 0) (move (box 1 1 1) {at} 0 0)))");
    assert_eq!(status(&src("2.2"), "clear"), Status::Unknown);
    assert_eq!(status(&src("4"), "clear"), Status::Holds);
    assert_eq!(status(&src("0.1"), "clear"), Status::Fails);
}
