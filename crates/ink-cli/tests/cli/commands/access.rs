//! A range of the question, `(ranged lo hi)`: claims are decided over all
//! of it, nothing narrows it, and commands showing a value look at one.

use super::*;

const WITHDRAW: &str = "(import \"std\")\n(def out (ranged 0 50))\n(def tool (move (cylinder-z 3 20) 0 0 (+ 10 out)))\n\
(def ledge (move (box 20 20 4) 6 0 43))\n(assert withdraws (disjoint tool ledge))\n";

fn msg(d: &Value) -> String {
    d["error"]["message"].as_str().unwrap_or("").to_string()
}

#[test]
fn a_ranged_value_stays_ranged_under_nominal_and_refuses_a_param_where_claims_are_answered() {
    // A tool withdrawn 0..50 mm past a ledge 41 mm up fails either way (as
    // a (tol 0 50) param it held under --nominal).
    let s = Scratch::new("withdraw");
    let w = s.model("withdraw.ink", WITHDRAW);
    let (a, b) = (ink(["assert", &w]), ink(["assert", &w, "--nominal"]));
    let (a, b) = (&a["asserts"][0], &b["asserts"][0]);
    let out = num(&b["params"]["out"]);
    assert!(
        verdict(a) == "fails"
            && verdict(b) == "fails"
            && 21.0 < out
            && out < 45.0
            && b["param_spans"].as_object().map(|o| o.len()) == Some(1)
            && nums(&b["param_spans"]["out"]) == Some(vec![0.0, 50.0]),
        "{}\n{}",
        pretty(a),
        pretty(b)
    );
    // assert and interfere take no --param of it; eval looks at its
    // witness's value, or the one given, within its range.
    for c in ["assert", "interfere"] {
        let e = &ink_err([c, &w, "--param", "out=25"])["error"];
        assert!(e["message"].as_str().unwrap_or("").contains("'out' is (ranged 0 50)") && e["hint"].as_str().unwrap_or("").contains("ink eval"), "{c}: {}", pretty(e));
    }
    let inside = |p: &[&str]| {
        let mut args = vec!["eval", &w, "--target", "tool", "--points", "[[0,0,35]]"];
        args.extend_from_slice(p);
        records(&ink(args)["targets"][0]["results"])[0]["inside"].clone()
    };
    assert!(inside(&[]) == false && inside(&["--param", "out=25"]) == true, "eval at the witness's value, then at out=25");
    let far = msg(&ink_err(["eval", &w, "--param", "out=80"]));
    assert!(far.contains("lies in [0, 50]"), "{far}");
}

#[test]
fn check_lists_a_ranges_span_warns_of_its_use_at_load_and_shows_it_at_its_at() {
    // ink check evaluates its conditions at a value, lists the range, and
    // warns of one used at load; a (ranged lo hi at) is shown at at.
    let s = Scratch::new("ranged-check");
    let c = s.model(
        "cond.ink",
        "(import \"std\")\n(def out (ranged 0 50 10))\n(def tool (check (value (- 40 out) > 0) (move (cylinder-z 3 20) 0 0 (+ 10 out))))\n\
(def ledge (move (box 20 20 4) 6 0 43))\n(assert withdraws (disjoint tool ledge))\n(def n (len (range (floor (/ out 10)))))\n",
    );
    let (at, hi) = (ink(["check", &c]), ink(["check", &c, "--param", "out=45"]));
    let ranged: Vec<Option<Vec<f64>>> = at["params"].as_array().unwrap().iter().map(|p| nums(&p["ranged"])).collect();
    let warns = |d: &Value| d["warnings"].as_array().map_or(vec![], |w| w.iter().map(|x| x.as_str().unwrap().to_string()).collect::<Vec<_>>());
    let (wa, wh) = (warns(&at), warns(&hi));
    assert!(ranged == [Some(vec![0.0, 50.0])], "its range: {}", pretty(&at["params"]));
    assert!(wa.iter().any(|x| x.contains("'out' is (ranged lo hi), but its value is used at load time")), "{wa:?}");
    assert!(wh.len() == wa.len() + 1, "its condition fails at --param out=45: {wh:?} vs {wa:?}");
    let shown = ink(["eval", &c, "--target", "tool", "--points", "[[0,0,25]]"]);
    assert!(records(&shown["targets"][0]["results"])[0]["inside"] == true, "shown at its at (10): {}", pretty(&shown));
}

#[test]
fn nothing_narrows_a_ranged_value_neither_another_name_nor_a_4th_view_element() {
    // (No import sets one either: an import gives no values.)
    let s = Scratch::new("ranged-narrow");
    // A view written with its own param ((vec dx dy dz s)) says so.
    let old = s.model(
        "old-view.ink",
        "(import \"std\")\n(param s 0 (tol 0 1))\n(def c (subtract (box 20 20 20) (box 12 12 12)))\n\
(assert a (empty (hidden c (box 20 20 20) (clear-air c (box 20 20 20) 0.5) (vec (vec 0 0 1 s)))))\n",
    );
    let eo = msg(&ink_err(["check", &old]));
    assert!(eo.contains("drop its 4th element"), "{eo}");
    // --param names a (ranged lo hi) as its witness does, and reaches its
    // load-time uses too; another def of it is refused by --param
    // (load-time uses would read it at another value).
    let al = s.model("alias.ink", "(import \"std\")\n(def u (ranged 0 1))\n(def u2 u)\n(def n (if (< u 0.5) 0 1))\n(def s (sphere (+ 1 n u)))\n");
    let r = value_at(&al, "s", &serde_json::json!([0, 0, 0]), &["--param", "u=0.9"]);
    let e = msg(&ink_err(["eval", &al, "--target", "s", "--points", "[[0,0,0]]", "--param", "u2=0.9"]));
    assert!((r + 2.9).abs() < 1e-9, "--param u=0.9 reaches n (load time): {r}");
    assert!(e.contains("another name for the (ranged lo hi) 'u'"), "{e}");
}

#[test]
fn a_claims_generated_range_name_leaves_param_to_the_param_of_its_name() {
    // A generated name has a space, which no declared name has: a claim's
    // range named travel 0 leaves --param travel to the param travel; it is
    // itself looked at by --param.
    let s = Scratch::new("ranged-names");
    let t = s.model("travel.ink", "(assert travel (value (ranged 0 1) <= 1))\n(param travel 10)\n(def reading travel)\n");
    let o = serde_json::json!([0, 0, 0]);
    let r = value_at(&t, "reading", &o, &["--param", "travel=0.75"]);
    let mut names: Vec<String> = ink(["check", &t])["params"].as_array().unwrap().iter().map(|p| p["name"].as_str().unwrap().to_string()).collect();
    names.sort();
    let generated = s.model("gen.ink", "(import \"std\")\n(def r (sphere (+ 1 (ranged 0 1))))\n");
    let looked = [value_at(&generated, "r", &o, &[]), value_at(&generated, "r", &o, &["--param", "r 0=0.5"])];
    assert!(r == 0.75 && names == ["travel", "travel 0"] && looked == [-1.0, -1.5], "{r} {names:?} {looked:?}");
}

#[test]
fn an_access_claim_fails_under_nominal_its_views_named_by_the_claim() {
    // std's hidden ranges each view over its own (ranged 0 1), named by the
    // claim: under --nominal an access claim is decided as on a normal run.
    let d = ink(["assert", &root("examples/machining-access.ink"), "--nominal", "--name", "t-slot-reaches"]);
    let a = &d["asserts"][0];
    let mut params: Vec<&str> = a["params"].as_object().unwrap().keys().map(String::as_str).collect();
    params.sort();
    assert!(verdict(a) == "fails" && params == ["t-slot-reaches 0", "t-slot-reaches 1"], "{}", pretty(a));
}

#[test]
fn distance_answers_at_the_values_asked_not_over_a_range() {
    let s = Scratch::new("ranged-distance");
    let dm = s.model(
        "dist.ink",
        "(import \"std\")\n(def out (ranged 0 50))\n(def tool (move (cylinder-z 3 20) 0 0 (+ 10 out)))\n\
(def ledge (move (box 20 20 4) 30 0 43))\n(def n (len (range (floor (/ out 10)))))\n",
    );
    for (x, want) in [(0, 27.0), (50, 17.69)] {
        let d = ink(["distance", &dm, "--a", "tool", "--b", "ledge", "--param", &format!("out={x}")]);
        let d = &d["distance"];
        assert!(num(&d["lower_proven"]) <= want && want <= num(&d["upper_witness"]), "out={x}: {want} not in {}", pretty(d));
    }
}

mod design {
    use super::*;

    const HEAD: &str = "(import \"std\")\n(param k 1)\n(def stock (move (box 30 30 10) 0 0 5))\n\
(def p (subtract stock (move (box (* 12 k) 12 20) 0 0 14)))\n(def h (hidden p stock (clear-air p stock 0.5) (vec (vec 0 0 1))))\n";

    #[test]
    fn an_optimize_reading_a_ranged_value_is_refused_saying_how_and_one_beside_it_runs() {
        // An optimize measures a design at one value of each param: one whose
        // requirement or objective reads a (ranged lo hi) is refused before
        // its first design, saying which. No (vary ...) varies one either.
        let s = Scratch::new("ranged-optimize");
        let vary = s.model("vary.ink", "(def out (ranged 0 50))\n(optimize o (vary out 0 10) (minimize (value out)) (budget 4))\n");
        let ev = msg(&ink_err(["check", &vary]));
        assert!(ev.contains("'out' is (ranged 0 50), a range of the question"), "{ev}");
        let req = s.model("opt-access.ink", &format!("{HEAD}(optimize o (vary k 0.8 1.2) (minimize (value k)) (require (empty h)) (budget 2))\n"));
        let obj = s.model("opt-access-objective.ink", &format!("{HEAD}(optimize o (vary k 0.8 1.2) (minimize (volume h)) (budget 2))\n"));
        let (r, b) = (msg(&ink_err(["optimize", &req])), msg(&ink_err(["optimize", &obj])));
        assert!(r.contains("(require (empty h)) reads h") && b.contains("its objective reads h") && (r.clone() + &b).contains("(ranged lo hi)"), "{r}\n{b}");
        // Read as a claim reads it: through a solve, or used at load (the
        // scope's error).
        let k = "(param k 1)\n(def u (ranged 0 1))\n";
        let through =
            s.model("opt-solve.ink", &format!("{k}(solve alias (unknowns (q 0)) (== q u))\n(optimize o (vary k 1 2) (minimize (value k)) (require (value q < 0.5)) (budget 4))\n"));
        let loaded = s.model("opt-load.ink", &format!("{k}(def q (if (< u 0.5) 0 1))\n(optimize o (vary k 1 2) (minimize (value k)) (require (value q < 0.5)) (budget 4))\n"));
        let (et, el) = (msg(&ink_err(["optimize", &through])), msg(&ink_err(["optimize", &loaded])));
        assert!(et.contains("(require (value q < 0.5)) reads u"), "{et}");
        assert!(el.contains("used at load time"), "{el}");
        // Beside them, one reading no ranged variable (k alone) is run, and feasible.
        let both = s.model(
            "opt-access-beside.ink",
            &format!("{HEAD}(optimize o (vary k 0.8 1.2) (minimize (value k)) (require (empty h)) (budget 2))\n(optimize plain (vary k 0.8 1.2) (minimize (value k)) (require (value k < 2)) (budget 2))\n"),
        );
        let c = ink(["optimize", &both, "--name", "plain"]);
        // Its measures are estimates: accuracy null, as a solve's.
        assert!(c["feasible"] == true && c.get("accuracy") == Some(&serde_json::Value::Null), "{}", pretty(&c));
    }
}
