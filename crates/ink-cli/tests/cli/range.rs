//! ink's range (docs/CORE.md Answers): a number outside it is seen however it reaches the model (written,
//! computed, read from JSON, a param's override or band), and every command
//! that exports a guarantee answers unknown or refuses, saying why. Each
//! beside an in-range twin that is still decided.

use crate::common::*;
use crate::mcp::Server;
use serde_json::{Value, json};

/// A square polygon of side s and its centre: its hull's cross products
/// underflow at s = 1e-170 and the polygon read empty there (holds,
/// guaranteed, for both claims), with a box of planes that says otherwise.
const SQUARE: &str = "(def f (polygon (vec (vec 0 0) (vec s 0) (vec s s) (vec 0 s))))
(def p (/ s 2))
(assert reported-outside (disjoint (probe p p 0) f))
(assert drawn-inside (subset (probe p p 0) (max (- 0 x) (- x s) (- 0 y) (- y s))))
";

/// Every claim of `d` (an assert's answer) unknown, saying `what` is outside.
fn outside(d: &Value, what: &str) {
    let all = d["asserts"].as_array().unwrap();
    assert!(
        !all.is_empty() && all.iter().all(|a| a["status"] == "unknown" && a["unavailable"].as_str().is_some_and(|u| u.contains(what) && u.contains("outside ink's range"))),
        "{what}: {}",
        pretty(d)
    );
}

#[test]
fn a_number_outside_the_range_is_seen_however_it_is_made() {
    let s = Scratch::new("range-numbers");
    s.model("numbers.json", r#"{"s": 1e-170, "one": 1}"#);
    let square = |name: &str, def: &str| s.model(name, &format!("(import \"std\")\n{def}\n{SQUARE}"));
    for (name, def) in [
        ("literal.ink", "(def s 1e-170)"),
        ("product.ink", "(def s (* 1e-85 1e-85))"),
        ("exp.ink", "(def s (exp -391.4394658089878))"),
        ("json.ink", "(def s (get (read-json \"numbers.json\") \"s\"))"),
    ] {
        outside(&s.ink(["assert", &square(name, def)]).any(), "e-17");
    }
    // At 1 mm, from the same file's JSON: decided, the probe inside.
    let d = s.ink(["assert", &square("one.ink", "(def s (get (read-json \"numbers.json\") \"one\"))")]).any();
    assert!(verdict(&d["asserts"][0]) == "fails" && verdict(&d["asserts"][1]) == "holds", "{}", pretty(&d));
    // A param overridden on the command line, at a value or pinned.
    let param = s.model("param.ink", &format!("(import \"std\")\n(param s 1)\n{SQUARE}"));
    for v in ["s=1e-170", "s=1e-170:1e-170"] {
        outside(&s.ink(["assert", &param, "--param", v]).any(), "e-170");
    }
    // The square scaled back to 1 mm by 1/s: the scaled square lost its
    // interior (a false holds and a false fails at accuracy 0).
    let scaled = s.model(
        "scaled.ink",
        "(import \"std\")
(def s (* 1e-85 1e-85))
(def tiny (polygon (vec (vec 0 0) (vec s 0) (vec s s) (vec 0 s))))
(def f (max (scale tiny (/ 1 s)) (slab-z -1 1)))
(assert outside (disjoint (probe 0.5 0.5 0) f))
(assert inside (subset (probe 0.5 0.5 0) f))
(assert empty (empty f))",
    );
    outside(&s.ink(["assert", &scaled]).any(), "1e-170");
    // A band computed past 1e100.
    let band = s.model("band.ink", "(param s 0 (tol (* 1e60 1e60) (* 2 1e60 1e60)))\n(assert positive (value s > 0))");
    outside(&s.ink(["assert", &band]).any(), "e119");
    // A param no constructor reads, overridden to a value or a band.
    let plain = s.model("plain.ink", "(param s 1)\n(assert positive (value s > 0))");
    outside(&s.ink(["assert", &plain, "--param", "s=1e-170"]).any(), "param 's' is 1e-170");
    outside(&s.ink(["assert", &plain, "--param", "s=1e-170:1"]).any(), "param 's' ranges to 1e-170");
    assert!(verdict(&s.ink(["assert", &plain, "--param", "s=1e-90:2e-90"]).any()["asserts"][0]) == "holds");
    // std's own constants (norm2-up's 2^1000 scaling) are not the model's.
    let norm = s.model("norm.ink", "(import \"std\")\n(assert five (value (norm2-up 3 4) >= 5))");
    assert!(verdict(&s.ink(["assert", &norm]).any()["asserts"][0]) == "holds");
    // Through MCP, the override as the command line's.
    let mut srv = Server::start(ROOT, &s, &[]);
    let (err, c) = srv.answer(2, json!({"command": "assert", "file": param, "options": {"param": ["s=1e-170"]}}));
    assert!(!err, "{}", pretty(&c));
    outside(&c, "e-170");
}

#[test]
fn every_command_that_proves_refuses_outside_the_range() {
    let s = Scratch::new("range-commands");
    let far = |name: &str, at: &str| s.model(name, &format!("(import \"std\")\n(body a (shape (sphere 1)))\n(body b (shape (move (sphere 1) {at} 0 0)))"));
    let (f, near) = (far("far.ink", "1e12"), far("near.ink", "100"));
    // A proven lower bound 2 mm above the real gap, and a box excluding
    // its interior, at 1e12.
    for cmd in [vec!["distance", f.as_str(), "--a", "a", "--b", "b"], vec!["bounds", f.as_str()]] {
        let d = s.ink(&cmd).err();
        assert!(d["error"]["message"].as_str().is_some_and(|m| m.contains("a coordinate 1")), "{cmd:?}: {}", pretty(&d));
    }
    let d = s.ink(["distance", &near, "--a", "a", "--b", "b"]).ok();
    let (lo, hi) = (d["distance"]["lower_proven"].as_f64().unwrap(), d["distance"]["upper_witness"].as_f64().unwrap());
    assert!(lo <= 98.0 && 98.0 <= hi, "{}", pretty(&d));
    let b = s.ink(["bounds", &near]).ok();
    let t = b["targets"].as_array().unwrap().iter().find(|t| t["target"] == "b").unwrap();
    assert!(t["min"][0].as_f64().unwrap() <= 99.0 && t["max"][0].as_f64().unwrap() >= 101.0, "{}", pretty(&b));
    // Interference's box shortcut: pairs apart by boxes are answered too.
    let pair = |name: &str, head: &str, at: &str| {
        s.model(
            name,
            &format!(
                "(import \"std\")\n{head}\n(body ball (shape (sphere 1)))\n(assembly root (instance a ball (at (move {at} 0 0))) (instance b ball (at (move (+ {at} 5) 0 0))))"
            ),
        )
    };
    for (m, what) in [(pair("tiny.ink", "(def metadata 1e-200)", "0"), "1e-200"), (pair("far-pair.ink", "", "1e12"), "a coordinate 1")] {
        let d = s.ink(["interfere", &m]).err();
        assert!(d["error"]["message"].as_str().is_some_and(|e| e.contains(what)), "{what}: {}", pretty(&d));
    }
    let apart = pair("apart.ink", "", "0");
    let d = s.ink(["interfere", &apart, "--res", "0.001"]).ok();
    assert!(d["pairs_apart_by_bounding_boxes"] == 1, "{}", pretty(&d));
    // The box shortcut asked at a resolution under a nanometre.
    let d = s.ink(["interfere", &apart, "--res", "1e-7"]).err();
    assert!(d["error"]["message"].as_str().is_some_and(|e| e.contains("a length of 1e-7")), "{}", pretty(&d));
    // A range claim asked at a resolution under a nanometre.
    let r = s.model("range.ink", "(import \"std\")\n(assert r (range x -2 2 (box 2 2 2)))");
    outside(&s.ink(["assert", &r, "--res", "1e-7"]).any(), "a length of 1e-7");
}

/// Values computed from numbers in the range, alike however they are
/// written (a def folds, a param is read in floats): judged where floats
/// lose them (1e-85 to the fourth is 0), not for their size (1e-85 squared,
/// 1e-170, is a float like any other). Over a band or a span values are
/// bounded outward, sound at every scale: their ends are judged. What Rust
/// builds from (an axis it normalizes) is judged for its size.
#[test]
fn a_value_computed_from_params_or_t_is_judged_where_floats_lose_it() {
    let s = Scratch::new("range-computed");
    let ask = |name: &str, src: &str, extra: &[&str]| {
        let m = s.model(name, src);
        s.ink(["assert", m.as_str()].into_iter().chain(extra.iter().copied())).any()
    };
    let holds = |d: &Value| d["asserts"].as_array().is_some_and(|a| !a.is_empty() && a.iter().all(|c| verdict(c) == "holds"));
    for head in ["(def p 1e-85)", "(param p 1e-85)", "(param p 1e85)", "(param p 1e-85 (tol 0 1e-86))"] {
        let d = ask("square.ink", &format!("{head}\n(assert positive (value (* p p) > 0))"), &[]);
        assert!(holds(&d), "{head}: {}", pretty(&d));
    }
    for head in ["(def p 1e-85)", "(param p 1e-85)"] {
        outside(&ask("fourth.ink", &format!("{head}\n(assert positive (value (* p p p p) > 0))"), &[]), "outside ink's range");
    }
    for span in ["1e-85 2e-85", "1e85 2e85"] {
        let d = ask("t-square.ink", &format!("(t-span {span})\n(assert positive (value (* t t) > 0))"), &[]);
        assert!(holds(&d), "{span}: {}", pretty(&d));
    }
    outside(&ask("t-tiny.ink", "(t-span 1e-170 2e-170)\n(assert positive (value t > 0))", &[]), "1e-170");
    outside(&ask("moment.ink", "(assert positive (value t > 0))", &["--t", "1e-170"]), "the moment 1e-170");
    // std's arithmetic on a number written or a param: alike.
    for p in ["1e85", "p"] {
        let d = ask("std.ink", &format!("(import \"std\")\n(param p 1e85)\n(assert positive (value (distance (vec {p} 0) (vec 0 0)) > 0))"), &[]);
        assert!(holds(&d), "{p}: {}", pretty(&d));
    }
    // An axis computed from a param, which Rust normalizes at load.
    let axis = "(import \"std\")\n(param axis 1e-81)\n(t-span 0 1)\n(body a (shape (move (sphere 1) 100 0 0)))\n(body b (shape (move (sphere 1) 123.5 0 0)))
(assembly root (instance a a (at (axes (origin 0 0 0) (x (* 1.8 axis axis) 0 0) (z 0 0 1)) (move (* 0.1 t) 0 0))) (instance b b))
(assert apart (disjoint (shape-of root/a) (shape-of root/b)))";
    outside(&ask("axis.ink", axis, &[]), "outside ink's range");
}

/// std's exemption covers its own constants, not a model's arithmetic on
/// what std returns: (norm2-up 0 0) is std's smallest float, 5e-324.
#[test]
fn stds_own_numbers_lend_no_exemption() {
    let s = Scratch::new("range-owned");
    let m = s.model(
        "owned.ink",
        "(import \"std\")
(def s (* (* (norm2-up 0 0) 1e100) 1e54))
(assert positive (value s > 0))
(def f (scale (box s s s) (/ 1 s)))
(assert inside (subset (probe 0 0 0) f))
(assert outside (disjoint (probe 1 0 0) f))
(assert empty (empty f))",
    );
    outside(&s.ink(["assert", &m]).any(), "outside ink's range");
}

/// Every interval written (an enclosure, a box, a band, a span) is exact
/// (command::bound): rounded to 10 digits, a solve's enclosure of 1/3
/// read [0.3333333333, 0.3333333333], excluding it, and a band of 1e-12
/// read [1, 1]. Through check, assert, interfere and MCP.
#[test]
fn every_interval_written_holds_what_it_bounds() {
    let s = Scratch::new("range-intervals");
    let solve = "(param p 1 (tol -1e-12 1e-12))\n(solve s (unknowns (q 0.3)) (== (* 3 q) p))\n(assert q-positive (value q > 0))\n";
    let m = s.model("solve.ink", solve);
    let bodies = s.model("bodies.ink", &format!("(import \"std\")\n{solve}(body ball (shape (sphere 1)))\n(assembly root (instance a ball) (instance b ball (at (move 5 0 0))))"));
    let holds = |iv: &serde_json::Value, x: f64, what: &str| {
        let (lo, hi) = (iv[0].as_f64().unwrap(), iv[1].as_f64().unwrap());
        assert!(lo <= x && x <= hi, "{what}: [{lo:e}, {hi:e}] excludes {x:e}");
    };
    let third = 1.0 / 3.0;
    let c = s.ink(["check", &m]).ok();
    holds(&c["solves"][0]["over_bands"]["q"], third, "check's over_bands");
    let a = s.ink(["assert", &m]).ok();
    holds(&a["solved_over_ranges"]["q"], third, "assert's solved_over_ranges");
    holds(&a["asserts"][0]["bounds_over_scope"], third, "bounds_over_scope");
    let i = s.ink(["interfere", &bodies]).ok();
    holds(&i["solved_over_ranges"]["q"], third, "interfere's solved_over_ranges");
    let band = &i["param_spans"]["p"];
    let (lo, hi) = (band[0].as_f64().unwrap(), band[1].as_f64().unwrap());
    assert!(lo < 1.0 && 1.0 < hi && ((hi - lo) - 2e-12).abs() < 1e-14, "the band checked, as checked: {}", pretty(band));
    let big = s.model("big.ink", "(param p 1e9 (tol -1e-6 1e-6))\n(solve s (unknowns (q 3.33333333333e8)) (== (* 3 q) p))\n");
    holds(&s.ink(["check", &big]).ok()["solves"][0]["over_bands"]["q"], 1e9 / 3.0, "check's over_bands at 1e9");
    let mut srv = Server::start(ROOT, &s, &[]);
    let (err, c) = srv.answer(2, json!({"command": "check", "file": m}));
    assert!(!err, "{}", pretty(&c));
    holds(&c["solves"][0]["over_bands"]["q"], third, "MCP check's over_bands");
}

/// A moving probe's field reads the coordinates that do not move as
/// spatial values (real-number bounds), never as the numbers they are,
/// folded in floats: (1e16 + y) - 1e16 - 0.25 at y = 0.5 would read -0.25
/// (a false holds and fails, 0.25 mm out, over a span on a 10 mm part), and
/// the twins underflow or overflow only in the substituted field. None
/// decides falsely.
#[test]
fn a_moving_probe_reads_its_fixed_coordinates_as_space() {
    let s = Scratch::new("range-probe");
    let ask = |name: &str, src: &str| statuses(&s.ink(["assert", &s.model(name, src)]).any());
    let claims = "(assert inside (subset (probe t p 0) f))\n(assert outside (disjoint (probe t p 0) f))";
    // (real sign of f at the probe, the model)
    for (name, outside, src) in [
        ("cancel.ink", true, "(param p 0.5)\n(def f (- (- (+ 1e16 y) 1e16) 0.25))"),
        ("tiny.ink", true, "(param p 2)\n(def s (* 1e-100 1e-100))\n(def f (- (/ (sqrt (square (* s y))) s) 1))"),
        ("cancel-bounded.ink", true, "(t-span 0 1)\n(param p 0.5)\n(def f (max (- (- (+ 1e16 y) 1e16) 0.25) (box 10 10 10)))"),
        ("tiny-bounded.ink", true, "(t-span 0 1)\n(param p 2)\n(def s (* 1e-100 1e-100))\n(def f (max (- (/ (sqrt (square (* s y))) s) 1) (box 10 10 10)))"),
        ("huge-bounded.ink", false, "(t-span 0 1)\n(param p 0.5)\n(def s (* 1e100 1e100))\n(def f (max (- (/ (sqrt (square (* s y))) s) 1) (box 10 10 10)))"),
        ("computed-cancel.ink", true, "(t-span 0 1)\n(param p 0.5)\n(def s (* 1e100 1e100))\n(def f (max (- (- (+ s y) s) 0.25) (box 10 10 10)))"),
    ] {
        let v = ask(name, &format!("(import \"std\")\n{src}\n{claims}"));
        let (wrong_in, wrong_out) = if outside { ("holds", "fails") } else { ("fails", "holds") };
        assert!(v["inside"] != wrong_in && v["outside"] != wrong_out, "{name}: {v:?}");
    }
    // Written, 1e200 is outside ink's range.
    let v = ask("literal-cancel.ink", &format!("(import \"std\")\n(t-span 0 1)\n(param p 0.5)\n(def f (max (- (- (+ 1e200 y) 1e200) 0.25) (box 10 10 10)))\n{claims}"));
    assert!(v["inside"] == "unknown" && v["outside"] == "unknown", "{v:?}");
    // An ordinary probe moving through the same box is decided.
    let v = ask("ordinary.ink", &format!("(import \"std\")\n(t-span 0 1)\n(param p 0.5)\n(def f (max (- y 0.25) (box 10 10 10)))\n{claims}"));
    assert!(v["inside"] == "fails" && v["outside"] == "holds", "{v:?}");
    // A product of params floats overflow (p^4 at 1e100) is a number
    // floats lost.
    outside(&s.ink(["assert", &s.model("overflow.ink", "(param p 1e100)\n(assert positive (value (/ (* p p p p) (* p p p p)) > 0))")]).any(), "outside ink's range");
}

/// Guaranteed numbers are written exactly (command::n; only `measured`
/// rounds): an accuracy of 1.00000000004 read 1.0, below the violation it
/// passed (1.00000000002), and a sampled field's error bound read
/// 0.1732051429, below 0.17320514290073.
#[test]
fn guaranteed_numbers_are_written_exactly() {
    let s = Scratch::new("range-exact");
    let m = s.model("accuracy.ink", "(assert tolerant (value 1.00000000002 <= 0) (accuracy 1.00000000004))");
    let d = s.ink(["assert", &m]).any();
    assert!(d["asserts"][0]["accuracy"].as_f64() == Some(1.00000000004), "{}", pretty(&d));
    let m = s.model("sampled.ink", "(import \"std\")\n(def sampled (grid-bake (sphere 1) (res 0.2)))\n(body sample (shape sampled))\n(assert a (subset (probe 0 0 0) sampled))");
    let c = s.ink(["check", &m]).ok();
    let found = |v: &Value| -> Vec<f64> {
        let mut out = Vec::new();
        let mut stack = vec![v];
        while let Some(v) = stack.pop() {
            match v {
                Value::Object(o) => o.iter().for_each(|(k, x)| if k == "error_bound" { out.extend(x.as_f64()) } else { stack.push(x) }),
                Value::Array(a) => stack.extend(a),
                _ => {}
            }
        }
        out
    };
    let bounds = found(&c);
    assert!(!bounds.is_empty() && bounds.iter().all(|&b| b > 0.1732051429), "{bounds:?}: {}", pretty(&c));
}

/// interfere's guard judges the span the
/// run asks, its ends included, and the arithmetic of the target's contact
/// zones as of its parts: refused, CLI and MCP alike, beside the in-range
/// twin that is still decided.
#[test]
fn interfere_judges_the_span_asked_and_its_contact_zones() {
    let s = Scratch::new("range-interfere-guard");
    // A ball moving 1e-100 mm/s: 10 mm over 1e101 s, its coordinates safe.
    let slow = s.model(
        "time-range-guard.ink",
        "(import \"std\")\n(body ball (shape (sphere 1)))\n(assembly root (instance a ball) (instance b ball (at (move (+ 5 (* 1e-100 t)) 0 0))))\n",
    );
    let zone = "(import \"std\")\n(param u 1e-100)\n(body block (shape (box 10 10 10)))\n(assembly root (instance a block)";
    let underflow_zone = s.model("underflow-zone.ink", &format!("{zone} (instance b block (at (move 7 0 0))) (contact a b (move (box 40 40 40) (* u u u u) 0 0)))\n"));
    let underflow_body = s.model("underflow-body.ink", &format!("{zone} (instance b block (at (move (+ 7 (* u u u u)) 0 0))) (contact a b (box 40 40 40)))\n"));
    let in_range = s.model("in-range-zone.ink", &format!("{zone} (instance b block (at (move 7 0 0))) (contact a b (move (box 40 40 40) (* u 1e50) 0 0)))\n"));
    let refused: [(&str, Vec<&str>, &str); 5] = [
        ("an end of --t-span above 1e100", vec!["--t-span", "0:1e101"], "the moment 1e101 s is outside ink's range"),
        ("an end of --t-span below 1e-100", vec!["--t-span", "0:1e-170"], "the moment 1e-170 s is outside ink's range"),
        ("--t above 1e100", vec!["--t", "1e101"], "the moment 1e101 s is outside ink's range"),
        ("a zone's underflow", vec![], "an underflow"),
        ("a placement's underflow", vec![], "an underflow"),
    ];
    let files = [&slow, &slow, &slow, &underflow_zone, &underflow_body];
    let mut srv = Server::start(ROOT, &s, &[]);
    for (id, ((what, args, why), file)) in refused.iter().zip(files).enumerate() {
        let mut cli = vec!["interfere", file.as_str()];
        cli.extend(args);
        let e = s.ink(&cli).err();
        assert!(e["error"]["message"].as_str().is_some_and(|m| m.contains(why)), "{what}: {}", pretty(&e));
        let options: serde_json::Map<String, Value> = args.chunks(2).map(|kv| (kv[0].trim_start_matches("--").to_string(), json!(kv[1]))).collect();
        let (err, c) = srv.call(id as u64 + 1, json!({"command": "interfere", "file": file, "options": options}));
        let t = c[0]["text"].as_str().unwrap_or_default();
        assert!(err && t.contains(why), "{what}, over MCP: {t}");
    }
    // In range: the same ball over (0, 1), and the zone computed in range.
    let d = s.ink(["interfere", &slow, "--t-span", "0:1"]).ok();
    assert!(d["pairs_apart_by_bounding_boxes"] == 1, "{}", pretty(&d));
    let d = s.ink(["interfere", &in_range]).ok();
    assert!(d["overlapping_pairs"] == 0 && d["pairs"][0]["status"] == "holds" && d["pairs"][0]["accuracy"] == 0.05, "{}", pretty(&d));
}
