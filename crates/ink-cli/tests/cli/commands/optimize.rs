//! `ink optimize`'s answers: each requirement of the design returned is
//! satisfied exactly where its claim, asked as `ink assert` asks it, holds
//! (its `answer`, which says why where it does not); an objective with no
//! value of the design (one a claim on it would leave unknown) fails it.

use super::*;

/// A model: a param p varied over [1, 2] from 2, minimized, with these
/// definitions and requirement (any, or none).
fn over_p(defs: &str, objective: &str, require: &str) -> String {
    format!("(import \"std\")\n(param p 2)\n{defs}\n(optimize o (vary p 1 2) (minimize {objective}) {require} (budget 4))\n")
}

/// `ink optimize` of `src`: whether it is feasible, and each requirement's
/// (satisfied, its claim's status).
fn optimized(s: &Scratch, src: &str) -> (Value, Vec<(bool, String)>) {
    let d = ink(["optimize", &s.model("o.ink", src)]);
    let rows = d["requires"].as_array().unwrap().iter().map(|r| (r["satisfied"] == true, r["answer"]["status"].as_str().unwrap_or("none").to_string())).collect();
    (d, rows)
}

#[test]
fn a_requirement_is_satisfied_exactly_where_its_claim_holds() {
    let s = Scratch::new("optimize-claims");
    let box2 = "(def f (box 2 2 2))";
    for (defs, require, want, why) in [
        // As written: 1 < 1 holds at no design, p > 2 at none of [1, 2]; 1 <= 1 at all.
        ("", "(require (value 1 < 1))", (false, "fails"), "a strict bound at its bound"),
        ("", "(require (value p > 2))", (false, "fails"), "a strict bound at the range's end"),
        ("", "(require (value 1 <= 1))", (true, "holds"), "a bound met at it"),
        // A count met at its bound: no voids, <= 0.
        (box2, "(require (cavities f <= 0))", (true, "holds"), "a count at a non-strict bound"),
        // A probe on the surface is not inside (field 0: not met).
        (box2, "(require (subset (probe 1 0 0) f))", (false, "fails"), "a point on the surface"),
        // A negative weight has no mean (not a mean of -2e-16, met).
        (box2, "(require (mean f x (weight -1) < 1))", (false, "unknown"), "a mean with no positive weight"),
        // Nothing near an empty shape: it holds, not "a target is empty".
        ("(def a 1) (def b -1)", "(require (clearance a b 1))", (true, "holds"), "a clearance from an empty shape"),
        // All space's volume is not the window's 1e9: unknown, not met.
        ("(def f -1)", "(require (volume f < 1e20))", (false, "unknown"), "a volume of a shape not proven bounded"),
        // The complement of a ball 2 m out has a void the window cuts off,
        // and a half-space no wall the window's faces make (neither met).
        ("(def f (neg (move (sphere 10) 2000 0 0)))", "(require (cavities f < 1))", (false, "unknown"), "a count of a shape not proven bounded"),
        ("", "(require (min-thickness x < 600))", (false, "unknown"), "a thickness of a shape not proven bounded"),
    ] {
        let (d, rows) = optimized(&s, &over_p(defs, "(value p)", require));
        assert!(rows == [(want.0, want.1.to_string())] && d["feasible"] == want.0, "{why}: {}", pretty(&d));
    }
}

#[test]
fn an_objective_with_no_value_of_the_design_fails_it() {
    // Of a shape not proven bounded the window is measured, not the shape:
    // a volume (not 1e9), a count of voids (not 0), a thickness (not the
    // half-window's 500 mm).
    let s = Scratch::new("optimize-objective-unbounded");
    for (defs, objective, says) in [
        ("(def f -1)", "(volume f)", "the objective's volume is of a shape not proven bounded"),
        ("(def f (neg (move (sphere 10) 2000 0 0)))", "(cavities f)", "the objective's cavities is of a shape not proven bounded"),
        ("", "(min-thickness x)", "the objective's min-thickness is of a shape not proven bounded"),
    ] {
        // Refused: the run (the gradient method stops at a design it cannot
        // evaluate), or its design (derivative-free, it says why).
        let d = Ink::new(["optimize", &s.model("o.ink", &over_p(defs, objective, ""))]).any();
        let msg = d["error"]["message"].as_str().or(d["evaluation_failed"].as_str()).unwrap_or("");
        assert!(msg.contains(says) && d["feasible"] != true, "{objective}: {}", pretty(&d));
    }
    // A bounded one is measured: the box's volume is p.
    let (d, _) = optimized(&s, &over_p("", "(volume (box p 1 1))", ""));
    assert!(d["feasible"] == true && (num(&d["objective"]) - num(&d["design"]["p"])).abs() < 0.05, "{}", pretty(&d));
}

#[test]
fn a_shapes_metric_is_an_objective() {
    // (min-thickness f) and (cavities f) load as objectives (not as a
    // form's metric: "no optimize or solve named 'f'"), and measure:
    // a ball of radius p is 2p thick; a box with a ball cut from it has one
    // void.
    let s = Scratch::new("shape-objective");
    for (def, objective, want) in [("(sphere p)", "(min-thickness f)", 2.0), ("(subtract (box 10 10 10) (sphere p))", "(cavities f)", 0.0)] {
        let (d, _) = optimized(&s, &over_p(&format!("(def f {def})"), objective, ""));
        let p = num(&d["design"]["p"]);
        let want = if want > 0.0 { want * p } else { 1.0 };
        assert!(d["feasible"] == true && (num(&d["objective"]) - want).abs() < 0.1, "{objective}: {}", pretty(&d));
    }
}

/// `ink assert FILE` with the answer's `apply` (its options, as written):
/// each assert's status.
fn asserted_at(file: &str, d: &Value) -> Vec<String> {
    let apply: Vec<&str> = d["apply"].as_str().unwrap().split_whitespace().collect();
    ink([&["assert", file][..], &apply].concat())["asserts"].as_array().unwrap().iter().map(|a| a["status"].as_str().unwrap().to_string()).collect()
}

#[test]
fn a_requirement_is_asked_over_its_files_motion_as_an_assert_is() {
    // Over (t-span 0 10) a box sliding through another meets it, and t
    // passes 1: each requirement fails, as its assert does (both would
    // hold asked of the design at t = 0 only). So in a file imported under a
    // namespace, over that file's span.
    let s = Scratch::new("optimize-motion");
    let slide = s.model(
        "slide.ink",
        "(import \"std\")\n(param p 0)\n(t-span 0 10)\n(def a (box 2 2 2))\n(def b (move (box 2 2 2) (- 10 t) 0 0))\n\
         (optimize o (vary p 0 1) (minimize (value p)) (require (disjoint a b)) (budget 4))\n(assert apart (disjoint a b))\n",
    );
    let d = ink(["optimize", &slide]);
    assert!(d["feasible"] == false && d["requires"][0]["answer"]["status"] == "fails" && asserted_at(&slide, &d) == ["fails"], "{}", pretty(&d));
    // The search's penalty is 0 at its one moment: the numbers are
    // labelled as the search's, a note says the answer decides, and the
    // warning names the requirement unmet, not the start.
    let warned = |d: &Value, what: &str| d["warnings"].as_array().unwrap().iter().any(|w| w.as_str().unwrap().contains(what));
    assert!(
        d["requires"][0]["search_violation"] == 0.0
            && d["requires"][0]["satisfied"] == false
            && d["note"].as_str().is_some_and(|n| n.contains("one moment (t 0)") && n.contains("is the verdict"))
            && warned(&d, "does not meet (require (disjoint a b)): over the span its `answer` does not hold, but the search's penalty for it was 0")
            && !warned(&d, "where it is met")
            && !warned(&d, "try another start"),
        "{}",
        pretty(&d)
    );
    // A strict bound at its bound, 0 > 0 at t 0: no penalty, and not met
    // there either; the warning says only what the penalty was.
    let strict = s.model("strict.ink", "(param p 0)\n(t-span 0 1)\n(optimize o (vary p 0 1) (minimize (value p)) (require (value t > 0)) (budget 4))\n");
    let d = ink(["optimize", &strict]);
    assert!(
        d["feasible"] == false
            && d["requires"][0]["search_violation"] == 0.0
            && d["requires"][0]["answer"]["status"] == "fails"
            && warned(&d, "does not meet (require (value t > 0)): over the span its `answer` does not hold, but the search's penalty for it was 0")
            && !warned(&d, " met at")
            && !warned(&d, "where it is met"),
        "{}",
        pretty(&d)
    );
    // Two unmet, one with no penalty at t 0 and one with a penalty: each
    // named with its own explanation, the second not as unsteered.
    let mixed = s.model("mixed.ink", "(param p 0)\n(t-span 0 1)\n(optimize o (vary p 0 1) (minimize (value p)) (require (value t < 0.5)) (require (value p < -1)) (budget 4))\n");
    let d = ink(["optimize", &mixed]);
    let ws: Vec<&str> = d["warnings"].as_array().unwrap().iter().filter_map(|w| w.as_str()).collect();
    let unsteered: Vec<&&str> = ws.iter().filter(|w| w.contains("was not steered")).collect();
    assert!(
        d["feasible"] == false
            && d["requires"][0]["search_violation"] == 0.0
            && d["requires"][1]["search_violation"] == 1.0
            && unsteered.len() == 1
            && unsteered[0].contains("(require (value t < 0.5))")
            && !unsteered[0].contains("(require (value p < -1))")
            && warned(&d, "does not meet (require (value p < -1)): its `answer` says why"),
        "{}",
        pretty(&d)
    );
    let moment =
        s.model("moment.ink", "(param p 0)\n(t-span 0 10)\n(optimize o (vary p 0 1) (minimize (value p)) (require (value t < 1)) (budget 4))\n(assert early (value t < 1))\n");
    let d = ink(["optimize", &moment]);
    assert!(d["feasible"] == false && d["requires"][0]["answer"]["status"] == "fails" && asserted_at(&moment, &d) == ["fails"], "{}", pretty(&d));
    let parent = s.model("parent.ink", "(import \"moment.ink\" as child)\n(param r 0)\n");
    let d = ink(["optimize", &parent, "--name", "child.o"]);
    assert!(d["feasible"] == false && d["requires"][0]["answer"]["status"] == "fails", "{}", pretty(&d));
}

#[test]
fn the_design_returned_is_given_exactly_as_it_was_asked_of() {
    // p = 0.99999999996 keeps the box at 0, where (floor p) jumps at 1:
    // rounded to ten digits the answer would give p = 1, a box 100 mm off, where
    // its assert fails. The design and `apply` are the values verified, so
    // ink assert with `apply` answers as each requirement does.
    let s = Scratch::new("optimize-exact");
    for (src, verdict) in [
        (
            "(import \"std\")\n(param p 0.99999999996)\n(def f (move (box 2 2 2) (* 100 (floor p)) 0 0))\n\
             (optimize o (vary p 0 1) (maximize (value p)) (require (within f (-2 -2 -2) (2 2 2))) (budget 4))\n(assert fit (within f (-2 -2 -2) (2 2 2)))\n",
            Some(true),
        ),
        // Which design the search returns is its own; whichever, the answer
        // and ink assert agree on it.
        ("(param p 0.99999999996)\n(optimize o (vary p 0 1) (maximize (value p)) (require (value (floor p) < 1)) (budget 4))\n(assert below (value (floor p) < 1))\n", None),
    ] {
        let m = s.model("exact.ink", src);
        let d = ink(["optimize", &m]);
        let r = &d["requires"][0];
        let p = d["design"]["p"].as_f64().unwrap();
        let apply = d["apply"].as_str().unwrap();
        let back: f64 = apply.strip_prefix("--param p=").unwrap().split_whitespace().next().unwrap().parse().unwrap();
        let want = if r["satisfied"] == true { "holds" } else { "fails" };
        assert!(back == p && asserted_at(&m, &d) == [want] && verdict.is_none_or(|v| (r["satisfied"] == true) == v), "{src}: {}", pretty(&d));
    }
    // The run's own --param values are in `apply` too.
    let m = s.model("given.ink", "(param p 0)\n(param q 1)\n(optimize o (vary p 0 1) (minimize (value p)) (require (value q > 1.5)) (budget 4))\n(assert big (value q > 1.5))\n");
    let d = ink(["optimize", &m, "--param", "q=2"]);
    assert!(d["feasible"] == true && d["apply"].as_str().is_some_and(|a| a.contains("q=2")) && asserted_at(&m, &d) == ["holds"], "{}", pretty(&d));
}

#[test]
fn an_imported_optimize_is_run_by_its_own_name_and_a_shared_one_is_refused() {
    // child.o varies q under an impossible 1 < 1; the run file's o varies p
    // under 1 <= 1. --name child.o runs child.o, not the local o (the
    // form's own name, `o`), and with no local o is not refused.
    let s = Scratch::new("optimize-names");
    s.model("child.ink", "(param q 0)\n(optimize o (vary q 0 1) (minimize (value q)) (require (value 1 < 1)) (budget 4))\n");
    let local = "(optimize o (vary p 0 1) (minimize (value p)) (require (value 1 <= 1)) (budget 4))\n";
    for parent in [format!("(import \"child.ink\" as child)\n(param p 0)\n{local}"), "(import \"child.ink\" as child)\n(param p 0)\n".to_string()] {
        let m = s.model("parent.ink", &parent);
        let d = ink(["optimize", &m, "--name", "child.o"]);
        assert!(
            d["optimize"] == "child.o" && d["feasible"] == false && d["requires"][0]["claim"] == "(value 1 < 1)" && d["design"].get("child.q").is_some(),
            "{parent}: {}",
            pretty(&d)
        );
    }
    let m = s.model("parent.ink", &format!("(import \"child.ink\" as child)\n(param p 0)\n{local}"));
    let d = ink(["optimize", &m, "--name", "o"]);
    assert!(d["optimize"] == "o" && d["feasible"] == true && d["requires"][0]["claim"] == "(value 1 <= 1)", "{}", pretty(&d));
    // Imported with no namespace, two are named o: refused, not either.
    let flat = s.model("flat.ink", &format!("(import \"child.ink\")\n(param p 0)\n{local}"));
    let e = ink_err(["optimize", &flat, "--name", "o"]);
    assert!(e["error"]["message"].as_str().is_some_and(|m| m.contains("more than one optimize is named 'o'")), "{}", pretty(&e));
}

#[test]
fn a_requirement_is_asked_over_its_own_files_span_not_a_same_named_forms() {
    // An imported file's (solve o) over (t-span 0 0.5) beside the run file's
    // optimize o over (t-span 0 10): the requirement is asked over the run
    // file's span, not the solve's file's (looked up by name alone), where
    // (floor t) is 0.
    let s = Scratch::new("optimize-span-owner");
    s.model("solver.ink", "(t-span 0 0.5)\n(solve o (unknowns (q 0)) (== q 0))\n");
    let m = s.model(
        "owner.ink",
        "(import \"solver.ink\")\n(param p 0)\n(t-span 0 10)\n(optimize o (vary p 0 1) (minimize (value p)) (require (value (floor t) < 1)) (budget 4))\n(assert r (value (floor t) < 1))\n",
    );
    let d = ink(["optimize", &m]);
    assert!(d["feasible"] == false && d["requires"][0]["answer"]["status"] == "fails" && asserted_at(&m, &d) == ["fails"], "{}", pretty(&d));
}

#[test]
fn a_requirement_that_holds_over_a_motion_is_answered() {
    // Proven over t ranged, its answer names t, a param of the view that
    // asked it and not of the run's (printing it indexed past the run's
    // params: a panic).
    let s = Scratch::new("optimize-motion-holds");
    let m = s.model("holds.ink", "(param p 0)\n(t-span 0 10)\n(optimize o (vary p 0 1) (minimize (value p)) (require (value t < 11)) (budget 4))\n(assert r (value t < 11))\n");
    let d = ink(["optimize", &m]);
    let a = &d["requires"][0]["answer"];
    assert!(d["feasible"] == true && a["status"] == "holds" && a["motion"]["span"] == serde_json::json!([0.0, 10.0]) && asserted_at(&m, &d) == ["holds"], "{}", pretty(&d));
}

#[test]
fn apply_asks_assert_what_the_design_was_verified_for() {
    // A pin (q=0:0) drops q's tolerance, a fixed or a varied param's alike;
    // apply giving q=0 would bring it back, where assert fails. And a run at
    // a moment (--t) is verified at it: apply says so.
    let s = Scratch::new("optimize-apply");
    let fixed = s.model(
        "fixed.ink",
        "(param p 0)\n(param q 0 (tol 0 1))\n(optimize o (vary p 0 1) (minimize (value p)) (require (value q < 0.5)) (budget 4))\n(assert r (value q < 0.5))\n",
    );
    let d = ink(["optimize", &fixed, "--param", "q=0:0"]);
    assert!(d["feasible"] == true && d["apply"].as_str().is_some_and(|a| a.contains("q=0:0")) && asserted_at(&fixed, &d) == ["holds"], "{}", pretty(&d));
    let varied = s.model("varied.ink", "(param p 0 (tol 0 1))\n(optimize o (vary p 0 1) (minimize (value p)) (require (value p < 0.5)) (budget 4))\n(assert r (value p < 0.5))\n");
    let d = ink(["optimize", &varied, "--param", "p=0:0"]);
    assert!(d["feasible"] == true && d["apply"].as_str().is_some_and(|a| a.contains("p=0:0")) && asserted_at(&varied, &d) == ["holds"], "{}", pretty(&d));
    let moment = s.model("moment.ink", "(param p 0)\n(t-span 0 10)\n(optimize o (vary p 0 1) (minimize (value p)) (require (value t < 1)) (budget 4))\n(assert r (value t < 1))\n");
    let d = ink(["optimize", &moment, "--t", "0"]);
    assert!(d["feasible"] == true && d["apply"].as_str().is_some_and(|a| a.ends_with("--t 0")) && asserted_at(&moment, &d) == ["holds"], "{}", pretty(&d));
    // A value for a range of the question is refused, as assert refuses it.
    let ranged = s.model("ranged.ink", "(param p 0)\n(def q (ranged 0 1))\n(optimize o (vary p 0 1) (minimize (value p)) (budget 4))\n");
    let e = ink_err(["optimize", &ranged, "--param", "q=0.5"]);
    assert!(e["error"]["message"].as_str().is_some_and(|m| m.contains("a range of the question")), "{}", pretty(&e));
}

#[test]
fn a_param_name_is_one_plain_name_the_command_line_sets_unquoted() {
    // A param's name is a plain word (letters, digits, '-' and '_',
    // starting with a letter), its namespaces joined by '.': written,
    // generated by a `for` (inline or from JSON) or a solve's unknown, or an
    // import's namespace. Every other is a name `apply` could not replay:
    // p=0,q sets p and q; " q" is trimmed to q; NUL cannot be an argument;
    // p'q' and p\q lose their quotes in a shell, all where assert would then
    // fail what the run verified. Each is refused at its line.
    let s = Scratch::new("optimize-param-names");
    s.model("child.ink", "(param p 0)\n");
    let opt = |n: &str| format!("(optimize o (vary {n} 0 1) (maximize (value {n})) (require (value {n} > 0.5)) (budget 4))\n");
    let mut cases: Vec<(String, &str)> = vec![
        (format!("(param p 0)\n(param q 0)\n(param p=0,q 0)\n{}", opt("p=0,q")), "p=0,q"),
        ("(param p,q 0)\n".into(), "p,q"),
        ("(param p=q 0)\n".into(), "p=q"),
        ("(param --p 0)\n".into(), "--p"),
        ("(param p'q' 0)\n".into(), "p'q'"),
        ("(param p\\q 0)\n".into(), "p\\q"),
        ("(param p\u{0}q 0)\n".into(), "p\u{0}q"),
        ("(solve s (unknowns (a,b 0)) (== a,b 1))\n".into(), "a,b"),
        ("(import \"child.ink\" as rev,2)\n".into(), "rev,2"),
    ];
    // Generated: by a `for`, inline and from JSON, as a param, an unknown
    // and a namespace.
    for (i, n) in [" q", "q ", "\tq", "q\n", "\u{a0}q", "p q", "p\"q", "p\u{0}q"].into_iter().enumerate() {
        // An ink string: only '"' and '\' escaped, the rest as they are.
        let lit = format!("\"{}\"", n.replace('\\', "\\\\").replace('"', "\\\""));
        cases.push((format!("(param q 0)\n(for n (vec {lit})\n  (param $n 0)\n  {})\n", opt("$n")), n));
        s.model(&format!("names-{i}.json"), &serde_json::to_string(&[n]).unwrap());
        cases.push((format!("(param q 0)\n(for n (read-json \"names-{i}.json\")\n  (param $n 0))\n"), n));
        cases.push((format!("(for n (vec {lit}) (solve s (unknowns ($n 0)) (== $n 1)))\n"), n));
    }
    cases.push(("(for n (vec \" ns\") (import \"child.ink\" as $n))\n".into(), " ns"));
    for (src, name) in &cases {
        let e = Ink::new(["check", &s.model("named.ink", src)]).any();
        let m = e["error"]["message"].as_str().unwrap_or("");
        assert!(m.contains(&format!("{name:?} is not a plain name")) && e["error"]["at"]["line"].as_u64().is_some(), "{src:?}: {}", pretty(&e));
    }
    // optimize's own --param values are such names too (apply gives them
    // back): not one a file cannot write, a claim's range "r 0".
    let m = s.model("plain.ink", &format!("(param p 0)\n(def r (+ p (ranged 0 1)))\n{}", opt("p")));
    let e = ink_err(["optimize", &m, "--param", "r 0=0.5"]);
    assert!(e["error"]["message"].as_str().is_some_and(|m| m.contains("optimize sets a param a file names")), "{}", pretty(&e));
    // Plain names, '-', '_' and digits, and a namespace's '.', replay.
    let m = s.model("dotted.ink", &format!("(import \"child.ink\" as child-2)\n(param a_b-3 0)\n{}(assert r (value child-2.p > 0.5))\n", opt("child-2.p")));
    let d = ink(["optimize", &m]);
    assert!(d["feasible"] == true && d["apply"].as_str().is_some_and(|a| a.starts_with("--param child-2.p=")) && asserted_at(&m, &d) == ["holds"], "{}", pretty(&d));
}

/// A model whose design is a pressure (m, given in psi) that the
/// requirement reads as a volume past m = 1, optimized: its answer and
/// progress.
fn refused_run(s: &Scratch) -> (Value, Vec<Value>) {
    let src = "(param m 0)\n(body b (shape (- (max (abs x) (abs y) (abs z)) 0.5)))\n(optimize o (vary m 0 5) (maximize (value m))\n (require (volume (shape-of b) < (if (> m 1) m 2mm3))) (budget 40))\n";
    let f = s.model("o.ink", src);
    let log = s.path("o.progress");
    let d = s.ink(["optimize", &f, "--param", "m=50psi"]).env("INK_PROGRESS", &log).ok();
    (d, std::fs::read_to_string(&log).unwrap().lines().map(|l| serde_json::from_str::<Value>(l).unwrap()).collect())
}

/// A design the model refuses is no limit the answer hides: each one the
/// search tried is in its answer and its progress (the answer so far, as a
/// run stopped early gives it), with the error as the model gives it.
#[test]
fn designs_the_model_refuses_are_reported() {
    let s = Scratch::new("optimize-refused");
    let (d, lines) = refused_run(&s);
    let refused = |e: &Value| e["message"].as_str().is_some_and(|m| m.starts_with("volume takes a volume") && m.contains("(--param m)")) && e["at"]["line"] == 4;
    let failed = d["failed_designs"].as_array().cloned().unwrap_or_default();
    assert!(!failed.is_empty() && failed.iter().all(|r| r["design"]["m"].as_f64().unwrap() > 1.0 && refused(&r["error"])), "{}", pretty(&d));
    // Each design once.
    assert!(failed.iter().enumerate().all(|(i, r)| failed[..i].iter().all(|q| q["design"] != r["design"])), "{}", pretty(&d));
    assert!(d["warnings"].as_array().unwrap().iter().any(|w| w.as_str().unwrap().contains("could not be evaluated (failed_designs)")), "{}", pretty(&d));
    assert!(lines.iter().any(|l| l["status"] == "failed" && refused(&l["row"]["error"])), "no refused design in progress");
    assert!(lines.iter().any(|l| l["answer"]["failed_designs"].as_array().is_some_and(|f| !f.is_empty())), "no failed designs in the answer so far");
    // A run stopped early answers what it has so far: its failed designs warned of too.
    let warned = |l: &&Value| l["answer"]["warnings"].as_array().is_some_and(|w| w.iter().any(|w| w.as_str().is_some_and(|w| w.contains("could not be evaluated"))));
    assert!(lines.iter().filter(|l| l["answer"]["failed_designs"].is_array()).all(|l| warned(&l)), "a failed design in the answer so far with no warning");
}

/// Every answer's `apply`, the finished one's and each in progress (what a
/// run stopped early answers), gives a param in the unit it was given in:
/// reused, m is still a pressure.
#[test]
fn every_answer_applies_a_param_in_its_unit() {
    let s = Scratch::new("optimize-apply-unit");
    let (d, lines) = refused_run(&s);
    let in_unit = |a: &str| a.ends_with("MPa") || a.ends_with("psi");
    assert!(in_unit(d["apply"].as_str().unwrap()), "{}", pretty(&d));
    let applied: Vec<&str> = lines.iter().filter_map(|l| l["answer"]["apply"].as_str()).collect();
    assert!(!applied.is_empty() && applied.iter().all(|a| in_unit(a)), "{applied:?}");
}
