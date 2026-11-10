//! Text params and lettering end to end: a customer's text given by
//! --param (or an import, or MCP's options), checked against what the
//! param allows, and the letters made of it.

use crate::common::*;
use crate::mcp::Server;
use serde_json::{Value, json};

const SIGN: &str = "(import \"std\")\n\
(param name \"Ada\" (chars \"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz ,'\") (len 1 16))\n\
(param w 60)\n\
(def L (lettering name (font \"sans\" (wght 700)) 8))\n\
(def B (lettering-box L))\n\
(body sign (shape (extrude-z (lettering-2d L) 0 2)))\n\
(assert fits (value (- (nth B 2) (nth B 0)) < w))\n";

fn claim<'a>(d: &'a Value, name: &str) -> &'a Value {
    d["asserts"].as_array().and_then(|a| a.iter().find(|c| c["name"] == name)).unwrap_or_else(|| panic!("no claim {name}: {}", pretty(d)))
}

#[test]
fn check_lists_a_text_param_with_what_it_allows() {
    let s = Scratch::new("text-check");
    let f = s.model("sign.ink", SIGN);
    let d = s.ink(["check", &f, "--param", "name=Grace"]).ok();
    let p = d["params"].as_array().unwrap();
    let name = p.iter().find(|p| p["name"] == "name").unwrap_or_else(|| panic!("{}", pretty(&d)));
    assert!(
        name["default"] == "Ada" && name["value"] == "Grace" && name["len"] == json!([1, 16]) && name["chars"].as_str().is_some_and(|c| c.contains("xyz ,'")),
        "{}",
        pretty(&d)
    );
    assert!(d["defs"].as_array().is_some_and(|ds| !ds.iter().any(|x| x["name"] == "name")), "a text param is a param, not a def: {}", pretty(&d));
}

#[test]
fn a_text_param_takes_any_text_it_allows() {
    let s = Scratch::new("text-assert");
    let f = s.model("sign.ink", SIGN);
    let d = s.ink(["assert", &f, "--param", "name=Grace"]).any();
    assert!(d["all_hold"] == true && d["texts"] == json!({"name": "Grace"}), "{}", pretty(&d));
    // Quoted, a text keeps its commas; the option goes on past it.
    let d = s.ink(["assert", &f, "--param", "name=\"Lovelace, Ada\",w=200"]).any();
    assert!(d["all_hold"] == true && d["texts"]["name"] == "Lovelace, Ada", "{}", pretty(&d));
    // A longer name is wider: the claim at this text fails, saying so.
    let d = s.ink(["assert", &f, "--param", "name=Augusta Ada King"]).any();
    let fits = claim(&d, "fits");
    assert!(fits["status"] == "fails" && fits["value"].as_f64().unwrap_or(0.0) > 60.0 && d["texts"]["name"] == "Augusta Ada King", "{}", pretty(&d));
}

#[test]
fn a_text_not_allowed_is_refused_saying_what_is() {
    let s = Scratch::new("text-refused");
    let f = s.model("sign.ink", SIGN);
    let m = |d: Value| d["error"]["message"].as_str().unwrap_or_default().to_string();
    let e = m(s.ink(["assert", &f, "--param", "name=Grace!"]).err());
    assert!(e.contains("\"Grace!\" holds \"!\" (character 6), which name does not allow"), "{e}");
    let e = m(s.ink(["check", &f, "--param", "name=ABCDEFGHIJKLMNOPQ"]).err());
    assert!(e.contains("holds 17 characters; name holds 1 to 16"), "{e}");
    // A number param given words says so, as before text params.
    let d = s.ink(["check", &f, "--param", "w=Grace"]).err();
    assert!(m(d.clone()).contains("unknown unit 'Grace'") && d["error"]["hint"].as_str().is_some_and(|h| h.contains("w is a number param")), "{}", pretty(&d));
}

#[test]
fn an_import_gives_a_text_param_its_value() {
    let s = Scratch::new("text-import");
    s.model("sign.ink", SIGN);
    let two = s.model(
        "two.ink",
        "(import \"std\")\n(import \"sign.ink\" as s (name \"Grace\"))\n(assert five (value (len (char-codes s.name)) >= 5))\n(assert wide (value (- (nth s.B 2) (nth s.B 0)) > 20))\n",
    );
    let d = s.ink(["assert", &two, "--own"]).any();
    assert!(d["all_hold"] == true, "{}", pretty(&d));
    // --param overrides an import's value, as for numbers.
    let d = s.ink(["assert", &two, "--own", "--param", "s.name=Ann"]).any();
    assert!(claim(&d, "five")["status"] == "fails", "{}", pretty(&d));
    let bad = s.model("bad.ink", "(import \"std\")\n(import \"sign.ink\" as s (w \"wide\"))\n");
    let e = s.ink(["check", &bad]).err();
    assert!(e["error"]["message"].as_str().is_some_and(|m| m.contains("w is a number param, but this import gives it text")), "{}", pretty(&e));
}

#[test]
fn mcp_takes_a_text_param() {
    let s = Scratch::new("text-mcp");
    let f = s.model("sign.ink", SIGN);
    let mut srv = Server::start(ROOT, &s, &[]);
    let (err, d) = srv.answer(1, json!({"command": "assert", "file": f, "options": {"param": ["name=\"Lovelace, Ada\"", "w=200"]}}));
    assert!(!err && d["all_hold"] == true && d["texts"]["name"] == "Lovelace, Ada", "{}", pretty(&d));
    let (err, d) = srv.answer(2, json!({"command": "check", "file": f, "options": {"param": ["name=Grace!"]}}));
    assert!(err && d.to_string().contains("does not allow"), "{}", pretty(&d));
}

/// A param's declared cases: every claim is asked at each too, a claim
/// failing at one fails, and a case --param sets is not asked (the repro
/// of the issue model-test-cases: a count read at load).
#[test]
fn a_claim_is_asked_at_every_declared_case() {
    let s = Scratch::new("text-cases");
    let f = s.model(
        "cases.ink",
        "(import \"std\")\n(param n 1 (cases 2))\n\
         (def part (apply min (map (lambda (i) (move (box 2 2 2) (* i 10) 0 0)) (range n))))\n\
         (assert envelope (within part (-2 -2 -2) (2 2 2)))\n",
    );
    let d = s.ink(["assert", &f]).any();
    let e = claim(&d, "envelope");
    assert!(
        d["all_hold"] == false && e["status"] == "fails" && e["cases"][0]["case"] == json!({"n": 2.0}) && e["cases"][0]["status"] == "fails" && d["cases"] == json!([{"n": 2.0}]),
        "{}",
        pretty(&d)
    );
    let d = s.ink(["assert", &f, "--param", "n=1"]).any();
    assert!(d["all_hold"] == true && d.get("cases").is_none(), "{}", pretty(&d));
    // A text param's cases are texts it allows, listed by check.
    let t = s.model("sign.ink", &SIGN.replace("(len 1 16))", "(len 1 16) (cases \"Augusta Ada King\" \"I\"))"));
    let d = s.ink(["assert", &t]).any();
    let fits = claim(&d, "fits");
    assert!(
        fits["status"] == "fails" && fits["cases"].as_array().map(Vec::len) == Some(2) && fits["cases"][0]["status"] == "fails" && fits["cases"][1]["status"] == "holds",
        "{}",
        pretty(&d)
    );
    let c = s.ink(["check", &t]).ok();
    assert!(c["params"].as_array().unwrap().iter().any(|p| p["cases"] == json!(["Augusta Ada King", "I"])), "{}", pretty(&c));
}

/// A text case is asked at each combination of cases that allows it, every
/// param set there (its len may be another param's), and only there: a
/// combination that does not is listed, not asked, and the model builds
/// for every text the param allows (a case refused at load refused every
/// customer's text; one dropped there was never asked where another
/// param's case allowed it).
#[test]
fn a_text_case_is_asked_where_its_combination_allows_it() {
    let s = Scratch::new("text-case-allowed");
    let f = s.model("cap.ink", "(param maxlen 1 (cases 2))\n(param name \"A\" (len 1 maxlen) (cases \"AA\"))\n(assert fits (value (len (char-codes name)) < 2))\n");
    let d = s.ink(["check", &f, "--param", "name=B"]).ok();
    assert!(d["warnings"].as_array().is_some_and(|w| w.is_empty()), "{}", pretty(&d));
    let d = s.ink(["assert", &f]).any();
    let fits = claim(&d, "fits");
    assert!(fits["status"] == "fails" && d["all_hold_guaranteed"] == false && d["cases"] == json!([{"maxlen": 2.0}, {"maxlen": 2.0, "name": "AA"}]), "{}", pretty(&d));
    assert!(d["cases_not_allowed"] == json!([{"case": {"name": "AA"}, "why": d["cases_not_allowed"][0]["why"]}]), "{}", pretty(&d));
    assert!(d["cases_not_allowed"][0]["why"].as_str().is_some_and(|w| w.contains("holds 2 characters; name holds 1 to 1")), "{}", pretty(&d));
}

/// A run asking claims at declared cases reports each claim at each case
/// as its own (a call stopped at its time limit counts them: it said "9 of
/// 1 asserts done"), and each load of the model at a case with its time
/// (credited to nothing, a stopped run's finished claims took 0 ms).
#[test]
fn a_run_at_cases_reports_each_case_and_its_load() {
    let s = Scratch::new("text-cases-progress");
    let f = s.model("sign.ink", &SIGN.replace("(len 1 16))", "(len 1 16) (cases \"Augusta Ada King\" \"I\"))"));
    let log = s.path("cases.progress");
    s.ink(["assert", &f]).env("INK_PROGRESS", &log).any();
    let lines: Vec<Value> = std::fs::read_to_string(&log).unwrap().lines().map(|l| serde_json::from_str::<Value>(l).unwrap()).collect();
    let total = lines.iter().rev().find_map(|l| l["total"].as_u64());
    let done: Vec<&Value> = lines.iter().filter(|l| !l["done"].is_null()).collect();
    let loaded: Vec<&Value> = lines.iter().filter(|l| !l["loaded"].is_null()).collect();
    assert!(total == Some(3) && done.len() == 3, "one claim at its value and two cases: {lines:?}");
    assert!(done[1]["case"] == json!({"name": "Augusta Ada King"}) && done[2]["case"] == json!({"name": "I"}), "{lines:?}");
    assert!(
        loaded.len() == 2
            && loaded.iter().all(|l| l["ms"].as_f64().is_some_and(|ms| ms > 0.0))
            && loaded[0]["loaded"] == json!({"load": {"name": {"case": 1}}})
            && loaded[0]["values"] == json!({"name": "Augusta Ada King"}),
        "{lines:?}"
    );
}

/// Many params' cases make more combinations than a run asks: said, not an
/// overflow of their product (64 params of 2 cases: 3^64).
#[test]
fn too_many_combinations_are_refused_however_many() {
    let s = Scratch::new("text-cases-many");
    let src: String = (0..64).map(|i| format!("(param p{i} 0 (cases 1 2))\n")).collect::<String>() + "(assert ok (value p0 < 3))\n";
    let f = s.model("many.ink", &src);
    let e = s.ink(["assert", &f]).err();
    assert!(e["error"]["message"].as_str().is_some_and(|m| m.contains("more than the 256 combinations")), "{}", pretty(&e));
}

/// A case is its expression's value at each combination, every param set
/// there: one reading another param took that param's value at the file's
/// load alone, so "AA" was never asked and a claim failing there held.
#[test]
fn a_case_reading_another_param_is_its_value_at_each_combination() {
    let s = Scratch::new("text-case-reads");
    let f = s.model(
        "reads.ink",
        "(param maxlen 1 (cases 2))\n(param name \"A\" (len 1 maxlen) (cases (if (== maxlen 1) \"A\" \"AA\")))\n(assert fits (value (len (char-codes name)) < 2))\n",
    );
    let d = s.ink(["assert", &f]).any();
    assert!(claim(&d, "fits")["status"] == "fails" && d["cases"].as_array().is_some_and(|c| c.contains(&json!({"maxlen": 2.0, "name": "AA"}))), "{}", pretty(&d));
}

/// A text param's default a combination's other cases do not allow (a
/// shorter len) is that combination's refusal, not the run's error.
#[test]
fn a_default_a_combination_does_not_allow_is_that_combinations_refusal() {
    let s = Scratch::new("text-default-refused");
    let f = s.model("dflt.ink", "(param maxlen 2 (cases 1))\n(param name \"AA\" (len 1 maxlen) (cases \"A\"))\n(assert fits (value (len (char-codes name)) <= maxlen))\n");
    let d = s.ink(["assert", &f]).any();
    assert!(d["all_hold"] == true && d["cases_not_allowed"].as_array().is_some_and(|c| c.len() == 1 && c[0]["case"] == json!({"maxlen": 1.0})), "{}", pretty(&d));
    assert!(d["cases"] == json!([{"name": "A"}, {"maxlen": 1.0, "name": "A"}]), "{}", pretty(&d));
}

/// A combination not asked (its text param refuses the text it takes) is
/// no claim to count: a stopped run said 3 of 4, 5 of 18.
#[test]
fn a_combination_not_asked_leaves_the_progress_count() {
    let s = Scratch::new("text-cases-count");
    let f = s.model("cap.ink", "(param maxlen 1 (cases 2))\n(param name \"A\" (len 1 maxlen) (cases \"AA\"))\n(assert fits (value (len (char-codes name)) < 3))\n");
    let log = s.path("cap.progress");
    s.ink(["assert", &f]).env("INK_PROGRESS", &log).any();
    let lines: Vec<Value> = std::fs::read_to_string(&log).unwrap().lines().map(|l| serde_json::from_str::<Value>(l).unwrap()).collect();
    let total = lines.iter().rev().find_map(|l| l["total"].as_u64());
    let done = lines.iter().filter(|l| !l["done"].is_null()).count();
    assert!(total == Some(3) && done == 3, "{lines:?}");
}

/// A combination is loaded once, each case taking its value where the load
/// reads it: a load at the cases' values from the file's first load (cap
/// 2, name "A") was a step to {cap 2, name "AA"}, and failed it, though
/// "A" at cap 2 is no combination; its progress said that load.
#[test]
fn a_combination_is_one_load_at_its_cases_own_values() {
    let s = Scratch::new("text-case-one-load");
    let f = s.model(
        "stable.ink",
        "(param cap 1 (cases 2))\n(param name (if (== cap 1) \"A\" \"AA\") (len 1 2) (cases (if (== cap 1) \"A\" \"AA\")))\n\
         (def last (nth (char-codes name) (- cap 1)))\n(assert ok (value last > 0))\n",
    );
    let log = s.path("stable.progress");
    let d = s.ink(["assert", &f]).env("INK_PROGRESS", &log).any();
    assert!(d["all_hold"] == true && d["cases"].as_array().is_some_and(|c| c.contains(&json!({"cap": 2.0, "name": "AA"}))), "{}", pretty(&d));
    let lines: Vec<Value> = std::fs::read_to_string(&log).unwrap().lines().map(|l| serde_json::from_str::<Value>(l).unwrap()).collect();
    let loaded: Vec<&Value> = lines.iter().filter_map(|l| l.get("values").filter(|_| !l["loaded"].is_null())).collect();
    assert!(loaded.len() == 3 && !loaded.contains(&&json!({"cap": 2.0, "name": "A"})), "{lines:?}");
}

/// A text a combination's param refuses stops its load there, before
/// anything reads it: `(nth (char-codes name) 0)` of a refused "" failed
/// the whole run.
#[test]
fn a_refused_text_is_read_by_nothing() {
    let s = Scratch::new("text-case-refused-read");
    let f = s.model("nth.ink", "(param name \"A\" (len 1 1) (cases \"\"))\n(def first (nth (char-codes name) 0))\n(assert ok (value first > 0))\n");
    let d = s.ink(["assert", &f]).any();
    assert!(d["all_hold"] == true && d["cases_not_allowed"][0]["case"] == json!({"name": ""}) && d["cases"] == json!([]), "{}", pretty(&d));
}

/// Cases go on params every combination declares: one a `for` over a
/// param makes at some combinations alone is refused at the load that
/// shows it, saying what to do, whether it is missing at the params'
/// values (name-1 at cap 2; two params taking cases together make c-0)
/// or at a combination (name-0 at cap 0). Asked as each combination
/// loaded, a combination's own cases were missed past a refusal, and the
/// run held, guaranteed. A param a `for` makes with no cases takes its
/// value where it is declared.
#[test]
fn cases_on_a_param_some_combinations_declare_are_refused() {
    let s = Scratch::new("text-cases-made");
    let refused = |name: &str, src: &str, param: &str, at: &str| {
        let e = s.ink(["assert", &s.model(name, src)]).err();
        let m = e["error"]["message"].as_str().unwrap_or_default();
        assert!(m.contains(&format!("{param}'s (cases ...) are refused: {param} {at}")) && m.contains("as its own file"), "{}", pretty(&e));
    };
    let here = "is declared at this combination but not at the params' values";
    refused("made.ink", "(param cap 1 (cases 2))\n(for i (range cap) (param name-$i \"A\" (cases \"AA\")))\n(assert ok (value 1 > 0))\n", "name-1", here);
    // At the first combination missing it, name-0 at its value there.
    let gone = "(param cap 1 (cases 0))\n(for i (range cap) (param name-$i \"A\" (cases \"AA\")))\n(assert ok (value 1 > 0))\n";
    refused("gone.ink", gone, "name-0", "is declared at the params' values but not at this combination");
    let e = s.ink(["assert", &s.model("gone.ink", gone)]).err();
    assert!(e["error"]["message"].as_str().is_some_and(|m| m.starts_with("at the cases {\"cap\":{\"case\":1}}: ")), "{}", pretty(&e));
    refused(
        "past.ink",
        "(param a 0 (cases 1))\n(for i (range a) (param b-$i 0 (cases 1)))\n(param name \"A\" (len 0 (if (== a 0) 1 b-0)))\n(assert fits (value (if (== a 0) 0 b-0) < 1))\n",
        "b-0",
        here,
    );
    refused(
        "both.ink",
        "(param a 0 (cases 1))\n(param b 0 (cases 1))\n(for i (range (if (== (+ a b) 2) 1 0)) (param c-$i 0 (cases 1)))\n\
         (param name \"A\" (len 0 (if (== (+ a b) 2) c-0 1)))\n(assert fits (value (if (== (+ a b) 2) c-0 0) < 1))\n",
        "c-0",
        here,
    );
    let f = s.model("plain.ink", "(param cap 1 (cases 2))\n(for i (range cap) (param w-$i 1))\n(assert ok (value (if (== cap 1) 0 w-1) < 1))\n");
    let d = s.ink(["assert", &f]).any();
    assert!(claim(&d, "ok")["status"] == "fails" && d["cases"] == json!([{"cap": 2.0}]), "{}", pretty(&d));
}

/// A param's cases are as many at every combination as at the params'
/// values: b-0's second case, declared only at a=1, was never asked, and
/// a claim failing there held, guaranteed. Refused, naming both counts.
#[test]
fn a_param_with_more_cases_at_a_combination_is_refused() {
    let s = Scratch::new("text-cases-count");
    for (name, src) in [
        (
            "count.ink",
            "(param a 0 (cases 1))\n(for i (range (if (== a 0) 1 0)) (param b-$i 0 (cases 1)))\n(for i (range (if (== a 1) 1 0)) (param b-$i 0 (cases 1 2)))\n(assert ok (value b-0 < 2))\n",
        ),
        (
            "text.ink",
            "(param a 0 (cases 1))\n(for i (range (if (== a 0) 1 0)) (param b-$i \"A\" (cases \"BB\")))\n(for i (range (if (== a 1) 1 0)) (param b-$i \"A\" (cases \"BB\" \"CCC\")))\n\
             (assert ok (value (len (char-codes b-0)) < 3))\n",
        ),
    ] {
        let e = s.ink(["assert", &s.model(name, src)]).err();
        let m = e["error"]["message"].as_str().unwrap_or_default();
        assert!(m.starts_with("at the cases {\"a\":{\"case\":1}}: b-0 declares 2 cases at this combination but 1 at the params' values"), "{}", pretty(&e));
    }
    // Fewer there is refused alike.
    let fewer = "(param a 0 (cases 1))\n(for i (range (if (== a 0) 1 0)) (param b-$i 0 (cases 1 2)))\n(for i (range (if (== a 1) 1 0)) (param b-$i 0 (cases 1)))\n(assert ok (value b-0 < 3))\n";
    let e = s.ink(["assert", &s.model("fewer.ink", fewer)]).err();
    assert!(e["error"]["message"].as_str().is_some_and(|m| m.contains("b-0 declares 1 case at this combination but 2 at the params' values")), "{}", pretty(&e));
}

/// The combinations are the product of the params' cases, each param at
/// its value or a case, every one asked: a claim failing only where all
/// three take a case fails.
#[test]
fn the_combinations_are_the_product_of_the_cases() {
    let s = Scratch::new("text-cases-product");
    let f = s.model("prod.ink", "(param a 0 (cases 1 2))\n(param b 0 (cases 1))\n(param c \"A\" (cases \"BB\"))\n(assert ok (value (* a b (len (char-codes c))) < 4))\n");
    let d = s.ink(["assert", &f]).any();
    let cases = d["cases"].as_array().unwrap();
    let mut seen: Vec<String> = cases.iter().map(Value::to_string).collect();
    seen.sort();
    seen.dedup();
    assert!(cases.len() == 11 && seen.len() == 11 && cases[0..4].iter().all(|c| c.as_object().unwrap().len() == 1), "{}", pretty(&d));
    let ok = claim(&d, "ok");
    let failing: Vec<&Value> = ok["cases"].as_array().unwrap().iter().filter(|c| c["status"] == "fails").map(|c| &c["case"]).collect();
    assert!(ok["status"] == "fails" && failing == [&json!({"a": 2.0, "b": 1.0, "c": "BB"})], "{}", pretty(&d));
    // 257 combinations are more than a run asks.
    let many = format!("(param a 0 (cases {}))\n(assert ok (value a > -1))\n", (1..=257).map(|i| i.to_string()).collect::<Vec<_>>().join(" "));
    let e = s.ink(["assert", &s.model("many.ink", &many)]).err();
    assert!(e["error"]["message"].as_str().is_some_and(|m| m.contains("more than the 256 combinations")), "{}", pretty(&e));
}

/// A claim made only at some combinations (a `for` over a param) is
/// reported, of its cases alone (`at_values` false): fits-1, failing at
/// cap 2, was dropped and the run held, guaranteed. Progress counts it.
#[test]
fn a_claim_made_only_at_a_case_is_reported() {
    let s = Scratch::new("text-claims-made");
    let f = s.model("made.ink", "(param cap 1 (cases 2))\n(for i (range cap) (assert fits-$i (value i < 1)))\n");
    let log = s.path("made.progress");
    let d = s.ink(["assert", &f]).env("INK_PROGRESS", &log).any();
    let one = claim(&d, "fits-1");
    assert!(d["all_hold"] == false && one["status"] == "fails" && one["at_values"] == false && one["cases"][0]["case"] == json!({"cap": 2.0}), "{}", pretty(&d));
    let lines: Vec<Value> = std::fs::read_to_string(&log).unwrap().lines().map(|l| serde_json::from_str::<Value>(l).unwrap()).collect();
    let total = lines.iter().rev().find_map(|l| l["total"].as_u64());
    let done = lines.iter().filter(|l| !l["done"].is_null()).count();
    assert!(total == Some(3) && done == 3, "{lines:?}");
}

/// A text an imported file's param refuses at a combination is that
/// combination's refusal: the importing file stopped reading the import's
/// params ("has no param 'name'") and the run aborted.
#[test]
fn a_refusal_in_an_import_is_its_combinations() {
    let s = Scratch::new("text-cases-import");
    s.model("child.ink", "(param cap 2 (cases 1))\n(param name \"A\" (len 1 cap) (cases \"A\"))\n(assert ok (value (len (char-codes name)) <= cap))\n");
    let parent = s.model("parent.ink", "(import \"child.ink\" as c (name \"AA\"))\n");
    let d = s.ink(["assert", &parent]).any();
    assert!(d["all_hold"] == true && d["cases_not_allowed"].as_array().is_some_and(|c| c.iter().any(|x| x["case"] == json!({"c.cap": 1.0}))), "{}", pretty(&d));
}

/// A run whose only claim is made at a case does not also say no claim
/// was checked.
#[test]
fn a_claim_made_at_a_case_is_a_claim_checked() {
    let s = Scratch::new("text-cases-only-claim");
    let f = s.model("only.ink", "(param a 0 (cases 1))\n(for i (range a) (assert extra-$i (value 1 < 0)))\n");
    let d = s.ink(["assert", &f]).any();
    assert!(d["all_hold"] == false && d.get("no_claims").is_none() && claim(&d, "extra-0")["status"] == "fails", "{}", pretty(&d));
}

/// Lettering whose load is slow at its cases: four grids of its letters
/// (redistance) at the cell size `fine`, fast at its default, seconds at
/// each case (about 12 s a case on 20 cores, unstopped).
const SLOW_CASES: &str = "(import \"std\")\n(param fine 1 (cases 0.15 0.14 0.13 0.12))\n\
(def L (lettering \"INK\" (font \"sans\") 20))\n\
(def g (fold min inf (map (lambda (k) (redistance (extrude-z (- (lettering-2d L) (* k 0.01)) 0 2) (res fine))) (range 4))))\n\
(body sign (shape g))\n(assert tall (value 20 > 1))\n";

/// The call's time limit (INK_DEADLINE: searches stop 3 s before it, to
/// answer whole), `s` seconds from now.
fn limit_in(s: f64) -> String {
    (std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs_f64() + s).to_string()
}

/// A load stops at the call's time limit, inside a grid it is sampling:
/// the call answers in time, saying where the load stopped, marked stopped.
#[test]
fn a_slow_load_stops_at_the_time_limit() {
    let s = Scratch::new("text-load-limit");
    let f = s.model("slow.ink", SLOW_CASES);
    let t0 = std::time::Instant::now();
    let e = s.ink(["assert", &f, "--param", "fine=0.12"]).env("INK_NO_CACHE", "1").env("INK_DEADLINE", limit_in(3.5)).err();
    let took = t0.elapsed().as_secs_f64();
    assert!(
        took < 3.0 && e["stopped"] == true && e["error"]["message"].as_str().is_some_and(|m| m.contains("time limit")) && e["error"]["at"]["line"] == 4,
        "{took:.1} s: {}",
        pretty(&e)
    );
}

/// A grid baked past the call's time limit stops the load, whatever pass
/// of it the limit falls in (its samples, or the bound on their error),
/// never answering with a grid made after it: stopped at small shares of
/// its quickest of three whole loads, a load no machine finishes by then
/// (it would have to run several times faster), it is a stopped error at
/// the grid's form every time. No clock outside the load is compared to
/// the limit inside it.
#[test]
fn a_grid_baked_past_the_time_limit_stops_the_load() {
    let s = Scratch::new("bake-limit");
    let f = s.model("bake.ink", "(import \"std\")\n(def g (grid-bake (sphere 5) (res 0.05)))\n(body b (shape g))\n");
    let load = |limit: Option<f64>| {
        let mut c = s.ink(["check", &f, "--load"]).env("INK_NO_CACHE", "1");
        // Searches stop 3 s before the call's limit.
        if let Some(l) = limit {
            c = c.env("INK_DEADLINE", limit_in(3.0 + l));
        }
        let t0 = std::time::Instant::now();
        let d = c.any();
        (d, t0.elapsed().as_secs_f64())
    };
    let best = (0..3)
        .map(|_| {
            let (d, t) = load(None);
            assert!(d["ok"] == true, "{}", pretty(&d));
            t
        })
        .fold(f64::INFINITY, f64::min);
    for share in STOP_SHARES {
        let (d, _) = load(Some(share * best));
        assert!(d["ok"] == false && d["stopped"] == true && d["error"]["at"]["line"] == 2, "stopped at {share} of {best:.2} s: {}", pretty(&d));
    }
}

/// Where a test stops a run, as shares of its quickest whole run: late
/// enough to fall past its first steps, early enough that no run finishes
/// first.
const STOP_SHARES: [f64; 2] = [0.15, 0.25];

/// Past the call's time limit no case is loaded: the claims answered at
/// the file's values stand, each case not loaded is listed and leaves
/// every claim unknown there, and the call answers in time.
#[test]
fn cases_past_the_time_limit_are_not_loaded_and_said() {
    let s = Scratch::new("text-cases-limit");
    let f = s.model("slow.ink", SLOW_CASES);
    let t0 = std::time::Instant::now();
    let d = s.ink(["assert", &f]).env("INK_NO_CACHE", "1").env("INK_DEADLINE", limit_in(5.0)).any();
    let took = t0.elapsed().as_secs_f64();
    let tall = claim(&d, "tall");
    let cases = tall["cases"].as_array().map_or(&[][..], Vec::as_slice);
    assert!(
        took < 4.5
            && d["stopped"] == true
            && d["cases_not_loaded"].as_array().is_some_and(|c| !c.is_empty() && c.len() == cases.len() - d["cases"].as_array().map_or(0, Vec::len))
            && tall["status"] == "unknown"
            && d["all_hold"] == false
            && cases.iter().any(|c| c["hint"].as_str().is_some_and(|h| h.contains("time limit"))),
        "{took:.1} s: {}",
        pretty(&d)
    );
}
