//! `ink assert`: failing claims with their witness, unknowns with where
//! and what decides them, and which claims a file's own are (--own).

use super::*;

#[test]
fn failing_claims_fail_an_overlap_with_its_witness() {
    let s = Scratch::new("failing");
    let failing = s.model(
        "failing.ink",
        "(import \"std\")\n(body a (shape (sphere 10)))\n(body b (shape (move (sphere 10) 15 0 0)))\n\
(assert apart (disjoint (shape-of a) (shape-of b)))\n(assert small (volume (shape-of a) < 1000))\n",
    );
    let d = ink(["assert", &failing]);
    let (apart, small) = (&d["asserts"][0], &d["asserts"][1]);
    assert!(verdict(apart) == "fails" && apart.get("witness").is_some(), "the overlap fails with a witness: {}", pretty(apart));
    assert!(verdict(small) == "fails estimate", "the volume claim fails, an estimate: {}", pretty(small));
}

#[test]
fn a_failed_clearance_gives_its_gap_and_a_failed_overlap_its_depth() {
    // A 3 mm shank in a 6.2 mm hole: a true gap of 0.1, under the 0.2
    // claimed. The gap is the witness pair's distance, at least the true
    // one and under the clearance. Two balls of 10 at 15 apart overlap 5
    // deep: the witness lies depth into both (field values), so 0 < depth
    // <= 2.5, and value is -depth.
    let s = Scratch::new("gap");
    let m = s.model(
        "gap.ink",
        "(import \"std\")\n(def plate (subtract (box 40 40 5) (cylinder-z 3.1 10)))\n(def shank (move (cylinder-z 3 20) 0 0 -5.5))\n\
(assert fit (clearance plate shank 0.2))\n(assert apart (disjoint (sphere 10) (move (sphere 10) 15 0 0)))\n",
    );
    let d = s.ink(["assert", &m]).any();
    let (fit, apart) = (&d["asserts"][0], &d["asserts"][1]);
    let (p, q) = (&fit["witness_pair"]["in_first"], &fit["witness_pair"]["in_second"]);
    let pair = (0..3).map(|k| (num(&p[k]) - num(&q[k])).powi(2)).sum::<f64>().sqrt();
    let gap = num(&fit["gap"]);
    assert!(verdict(fit) == "fails" && (gap - pair).abs() < 1e-6 && (0.1 - 1e-6..0.2).contains(&gap) && (gap - num(&fit["value"])).abs() < 1e-9, "{}", pretty(fit));
    let depth = num(&apart["depth"]);
    assert!(verdict(apart) == "fails" && depth > 0.0 && depth <= 2.5 + 1e-9 && depth == -num(&apart["value"]), "{}", pretty(apart));
}

#[test]
fn assert_exits_3_when_a_claim_does_not_hold_and_0_when_all_do() {
    // A script reads the exit status: 0 every claim holds, 3 answered with
    // one that fails or is unknown, 1 an error (no answer).
    let s = Scratch::new("exit");
    let m = s.model("exit.ink", "(import \"std\")\n(body a (shape (sphere 10)))\n(assert big (volume (shape-of a) > 1000))\n(assert small (volume (shape-of a) < 1000))\n");
    let code = |args: &[&str]| s.ink(args).run().code;
    assert_eq!(code(&["assert", &m]), Some(3), "a failing claim");
    assert_eq!(code(&["assert", &m, "--name", "big"]), Some(0), "every claim holds");
    assert_eq!(code(&["assert", &m, "--name", "none"]), Some(1), "an error");
}

#[test]
fn an_unknown_exact_claim_gives_an_undecided_cell_and_a_res_hint() {
    // An unknown says where, as a failure does; one asked exactly (a range
    // at its bound on the face x = 0.5) is told what a res can decide, not
    // to state an accuracy it refuses.
    let s = Scratch::new("at-bound");
    let m = s.model("at-bound.ink", "(import \"std\")\n(assert at-bound (range x -1 0.5 (box 1 1 1)))\n");
    let d = ink(["assert", &m]);
    let a = &d["asserts"][0];
    let c = &a["undecided_cell"];
    let hint = a["hint"].as_str().unwrap_or("");
    assert!(verdict(a) == "unknown" && (num(&c["centre"][0]) - 0.5).abs() <= num(&c["size"][0]) && hint.contains("(res r)") && !hint.contains("accuracy"), "{}", pretty(a));
}

#[test]
fn own_runs_only_the_files_claims_and_check_imports_says_where_each_is_written() {
    // --own checks the claims written in the file, not the ones it imports
    // (here a false one), and says so. The model rule relies on it.
    let s = Scratch::new("own");
    let base = s.model("own-base.ink", "(import \"std\")\n(body b (shape (box 10 10 10)))\n(assert too-big (volume (shape-of b) > 2000))\n");
    let top = s.model("own-top.ink", "(import \"own-base.ink\" as base)\n(assert fits (within (shape-of base.b) (-6 -6 -6) (6 6 6)))\n");
    let d = ink(["check", &top, "--imports"]);
    let own: Vec<(String, bool, String)> =
        d["asserts"].as_array().unwrap().iter().map(|a| (a["name"].as_str().unwrap().into(), a["own"] == true, a["source"].as_str().unwrap().into())).collect();
    assert!(own == [("base.too-big".into(), false, format!("{base}:3")), ("fits".into(), true, format!("{top}:2"))], "{own:?}");
    let d = ink(["assert", &top]);
    let names: Vec<&str> = d["asserts"].as_array().unwrap().iter().map(|a| a["name"].as_str().unwrap()).collect();
    assert!(names == ["base.too-big", "fits"] && d["all_hold"] == false && d.get("scope").is_none(), "assert runs imported claims too: {}", pretty(&d));
    let d = ink(["assert", &top, "--own"]);
    let rows: Vec<(String, String)> = d["asserts"].as_array().unwrap().iter().map(|a| (a["name"].as_str().unwrap().into(), verdict(a))).collect();
    assert!(rows == [("fits".to_string(), "holds".to_string())] && d["all_hold_guaranteed"] == true && d["scope"] == "own" && d["imported_not_run"] == 1, "{}", pretty(&d));
    let e = &ink_err(["assert", &top, "--own", "--name", "base.too-big"])["error"];
    assert!(e["message"].as_str().unwrap().contains("imported"), "--own --name of an imported claim: {}", pretty(e));
}

#[test]
fn a_pin_raised_off_its_seat_is_not_seated() {
    // examples/frames.ink with the plate's bore-top frame 100 mm up: the
    // pin floats, and its seated probe fails (a tiny box passed here).
    let src = std::fs::read_to_string(root("examples/frames.ink")).unwrap();
    let (seat, aloft) = ("(at (move hole-x hole-y 4))", "(at (move hole-x hole-y 104))");
    assert!(src.contains(seat));
    let s = Scratch::new("pin-aloft");
    let m = s.model("frames-aloft.ink", &src.replace(seat, aloft));
    let d = ink(["assert", &m, "--name", "seated"]);
    assert!(verdict(&d["asserts"][0]) == "fails", "{}", pretty(&d["asserts"]));
}

#[test]
fn a_bearing_afloat_above_its_shoulder_does_not_rest() {
    // examples/frames.ink's 608 bearing (a 22 x 8 x 7 ring) on a
    // seat 0.4 above the shoulder: the bearing floats in a bore closer than
    // the accuracy, and its rests claim fails (the example's control).
    let s = Scratch::new("afloat");
    let m = s.model(
        "bearing-afloat.ink",
        "(import \"std\")\n(def ledge 4.5)\n\
(body bearing (shape (max (ring-2d 4 11) (slab-z -3.5 3.5))) (frame face-b (at (rot-x 180) (move 0 0 -3.5))))\n\
(body housing (shape (subtract (move (rounded-box 40 40 12 3) 0 0 6) (move (cylinder-z 11.02 20) 0 0 (+ ledge 10)) (cylinder-z 8 40)))\n\
  (frame seat (at (move 0 0 (+ ledge 0.4)))))\n\
(assembly root (instance housing housing) (instance bearing bearing (frame-of housing seat) face-b))\n\
(assert rests (touches (shape-of root/bearing) (max (shape-of root/housing) (- z ledge))))\n",
    );
    let d = ink(["assert", &m, "--own", "--name", "rests"]);
    let v: Vec<String> = d["asserts"].as_array().unwrap().iter().map(verdict).collect();
    assert!(v == ["fails"], "{}", pretty(&d["asserts"]));
}

#[test]
fn an_overlapping_clearance_fails_with_value_and_gap_zero() {
    // The pair's distance is not rounded up per coordinate where the two
    // points are one: that gives 1.5e-323 beside a gap of 0.
    let s = Scratch::new("clearance-overlap");
    let m = s.model(
        "overlap.ink",
        "(import \"std\")\n(body a (shape (box 10 10 10)))\n(body b (shape (move (box 10 10 10) 4 0 0)))\n\
        (assert clear (clearance (shape-of a) (shape-of b) 0.3))\n",
    );
    let d = &ink(["assert", &m])["asserts"][0];
    assert!(verdict(d) == "fails" && d["value"] == 0.0 && d["gap"] == 0.0, "{}", pretty(d));
}

#[test]
fn a_probe_is_decided_up_to_the_arithmetics_rounding_of_a_surface() {
    // docs/IR.md said a point within about 0.0025 of a surface is unknown;
    // it is decided, exactly, down to the field's rounding: past a face
    // fails a subset, inside holds, on it is outside.
    let s = Scratch::new("probe-near");
    let rows = |x: &str| {
        let m = s.model(
            &format!("p{x}.ink"),
            &format!("(import \"std\")\n(assert in (subset (probe {x} 0 0) (box 10 10 10)))\n(assert out (disjoint (probe {x} 0 0) (box 10 10 10)))\n"),
        );
        let d = ink(["assert", &m]);
        (verdict(&d["asserts"][0]), verdict(&d["asserts"][1]), d["asserts"][0]["accuracy"].clone())
    };
    for (x, want) in [("5.0001", ("fails", "holds")), ("4.9999", ("holds", "fails")), ("5.000000001", ("fails", "holds")), ("5", ("fails", "holds"))] {
        let (i, o, acc) = rows(x);
        assert!(i == want.0 && o == want.1 && acc == 0.0, "probe at x = {x}: subset {i}, disjoint {o}, accuracy {acc}");
    }
}

#[test]
fn a_file_with_no_claims_says_so_and_guarantees_nothing() {
    // It said all_hold_guaranteed: true, vacuously.
    let s = Scratch::new("no-claims");
    let m = s.model("none.ink", "(import \"std\")\n(body b (shape (sphere 2)))\n");
    let d = ink(["assert", &m]);
    assert!(d["all_hold"].is_null() && d["all_hold_guaranteed"].is_null() && d["no_claims"].is_string(), "{}", pretty(&d));
}

/// A drafted polygon's claims at a point: each claim's verdict, or why the
/// polygon is refused.
fn drafted_verdicts(tag: &str, defs: &str, claims: &[(&str, &str)]) -> Result<std::collections::BTreeMap<String, String>, String> {
    let s = Scratch::new(tag);
    let mut text = format!("(import \"std\")\n{defs}\n(body part (shape (max f (slab-z -1 1))))\n");
    for (name, claim) in claims {
        text += &format!("(assert {name} {claim})\n");
    }
    let m = s.model("m.ink", &text);
    let d = Ink::new(["assert", &m]).env("INK_NO_CACHE", "1").any();
    if d["ok"].as_bool() != Some(true) {
        return Err(d["error"]["message"].as_str().unwrap_or("?").to_string());
    }
    Ok(statuses(&d))
}

#[test]
fn a_drafted_reflex_join_turning_in_past_its_allowance_is_not_read_wrong() {
    // A shallow arc (radius 3.1e7) meets a
    // line at a join turned in 1e-6 rad at z = 0, the 45 degree draft
    // turning it in to 8.1e-6 at z = 0.001, two metres across. An
    // unconditional 1e-5 rad tolerance read the join as convex: the line's
    // half-plane cut the arc's material at (1100, 0.00115, 0.001), 3.9e-4
    // inside its circle (the allowance 2.2e-4), read outside at accuracy 0;
    // drafted the other way and reflected in z, the same geometry read it
    // inside. It is read inside both ways, or refused.
    let pts = "(vec (vec 0 0 7.9999999999999996e-06) (vec 1000 0) (vec 2000 0.015000000000783664) (vec 2000 1000) (vec 0 1000))";
    let defs = format!("(def f (polygon {pts} (vec 45 0 0 0 0) 0 0.0013))\n(def g (warp (polygon {pts} (vec -45 0 0 0 0) -0.0013 0) x y (neg z)))");
    let p = "(probe 1100 0.00115 0.001)";
    let claims = [("f-in", format!("(subset {p} f)")), ("f-out", format!("(disjoint {p} f)")), ("g-in", format!("(subset {p} g)")), ("g-out", format!("(disjoint {p} g)"))];
    let claims: Vec<(&str, &str)> = claims.iter().map(|(n, c)| (*n, c.as_str())).collect();
    match drafted_verdicts("turn-growing", &defs, &claims) {
        Err(e) => assert!(e.contains("polygon"), "refused, saying why: {e}"),
        Ok(v) => {
            assert!(v["f-in"] == "holds" && v["f-out"] == "fails" && v["g-in"] == "holds" && v["g-out"] == "fails", "{v:?}");
        }
    }
}

#[test]
fn a_drafted_polygon_s_box_holds_its_moving_wall() {
    // A square whose right wall leans out 45 degrees over z 0..0.1: clipped
    // at the float 10.1, it would exclude (10.1, 5, 0.1), which is inside the
    // wall in real arithmetic. The clipping box is rounded outward, so this
    // probe must not be proven disjoint at accuracy 0.
    let defs = "(def f (polygon (vec (vec 0 0) (vec 10 0) (vec 10 10) (vec 0 10)) (vec 0 45 0 0) 0 0.1))";
    let v = drafted_verdicts("box-boundary", defs, &[("out", "(disjoint (probe 10.1 5 0.1) f)")]).unwrap();
    assert!(v["out"] != "holds", "{v:?}");
}

#[test]
#[ignore = "a polygon's allowance is not yet in its answers' accuracy"]
fn a_polygon_s_allowance_is_in_its_answers() {
    // A corner 5e-7 below the line through its neighbours,
    // within polygon's allowance (1e-7 of its 1414 mm diagonal), is merged
    // into that line; the point 2.5e-7 inside the specified profile read
    // outside at accuracy 0, drafted and not. Within the allowance it is
    // unknown (neither holds nor fails), until an accuracy carries it.
    let pts = "(vec (vec 0 0) (vec 500 -0.0000005) (vec 1000 0) (vec 1000 1000) (vec 0 1000))";
    for defs in [format!("(def f (polygon {pts} 0 0 1))"), format!("(def f (polygon {pts}))")] {
        let p = "(probe 500 -2.4999999999999999e-07 0.5)";
        let (i, o) = (format!("(subset {p} f)"), format!("(disjoint {p} f)"));
        let v = drafted_verdicts("allowance-corner", &defs, &[("in", &i), ("out", &o)]).unwrap();
        assert!(v["in"] != "fails" && v["out"] != "holds", "{defs}: {v:?}");
    }
}

/// A promise is taken only for its constants as computed for the query:
/// at --param p=1e200 (fixed, or a range of one value) k = 0.001 + q / q
/// is 0.001 (q overflows, inf / inf is 0), and the promise taken as
/// written would guarantee a clearance between points 0.042 apart.
#[test]
fn a_promise_is_not_taken_where_its_params_overflow() {
    let f = fixture("promise-param-overflow.ink");
    assert_eq!(statuses(&ink(["assert", &f, "--param", "p=1"]))["gap"], "fails");
    // 1e200 is outside ink's range: unknown, never the false holds.
    for p in ["p=1e200", "p=1e200:1e200", "p=1:1e200"] {
        let s = statuses(&ink(["assert", &f, "--param", p]));
        assert_eq!(s["gap"], "unknown", "--param {p}");
    }
}
