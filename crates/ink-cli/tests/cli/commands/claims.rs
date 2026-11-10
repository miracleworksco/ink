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
    // at its bound on the face x = 0.5, the box reading x as (x + 3 y) -
    // 3 y, so its face and the bound share nothing to split at) is told
    // what a res can decide, not to state an accuracy it refuses.
    let s = Scratch::new("at-bound");
    let m = s.model("at-bound.ink", "(import \"std\")\n(assert at-bound (range x -1 0.5 (warp (box 1 1 1) (- (+ x (* 3 y)) (* 3 y)) y z)))\n");
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

/// A straight, undrafted polygon is its drawn corners exactly: a corner
/// 5e-7 below the line through its neighbours in a 1414 mm square is kept,
/// so a point 2.5e-7 inside reads inside at accuracy 0 (polygon's allowance
/// merged it into that line, and the point read outside), and one as far
/// outside, outside.
#[test]
fn a_straight_polygon_is_its_drawn_corners() {
    let defs = "(def f (polygon (vec (vec 0 0) (vec 500 -0.0000005) (vec 1000 0) (vec 1000 1000) (vec 0 1000))))";
    let (pin, pout) = ("(probe 500 -2.4999999999999999e-07 0.5)", "(probe 500 -7.5e-07 0.5)");
    let (i, o) = (format!("(subset {pin} f)"), format!("(disjoint {pout} f)"));
    let v = drafted_verdicts("straight-exact", defs, &[("in", &i), ("out", &o)]).unwrap();
    assert_eq!((v["in"].as_str(), v["out"].as_str()), ("holds", "holds"), "{v:?}");
}

/// A flanged sleeve's section, its corners rounded by quarter-turn arcs
/// 0.024 across: a low term of an exact circle's sums can fall outside
/// ink's range (-5.4e-101 here), and then every claim on the model reads
/// unknown. A circle is written as bounds on its centre less its arc's end
/// and on its radius, numbers of the arc's own size.
#[test]
fn an_exact_arc_writes_no_number_outside_ink_s_range() {
    let defs = "(def f (revolve-z (polygon (vec (vec 1.5875 0.0793751) (vec 1.666875 1e-07) (vec 3.9449375 1e-07 0.41421356) (vec 3.96875 0.0238126) (vec 3.96875 1.5636876 0.41421356) (vec 3.9449375 1.5875001) (vec 3.1988125 1.5875001 -0.41421356) (vec 3.175 1.6113126) (vec 3.175 6.1118751) (vec 2.936875 6.3500001) (vec 1.666875 6.3500001) (vec 1.5875 6.2706251)))))";
    let v = drafted_verdicts("arc-range", defs, &[("in", "(subset (probe 2.5 0 3) f)")]).unwrap();
    assert_eq!(v["in"], "holds", "{v:?}");
}

#[test]
fn a_drafted_corner_off_its_neighbours_line_is_kept() {
    // Drafted, a corner 5e-7 below the line through its neighbours (within
    // 1e-7 of the 1414 mm diagonal) is a corner: the point 2.5e-7 inside
    // the drawn profile is inside (merged into that line, it read outside
    // at accuracy 0).
    // So too with its far corners rounded.
    let p = "(probe 500 -2.4999999999999999e-07 0.5)";
    let (i, o) = (format!("(subset {p} f)"), format!("(disjoint {p} f)"));
    for (tag, defs) in [
        ("allowance-corner", "(def f (polygon (vec (vec 0 0) (vec 500 -0.0000005) (vec 1000 0) (vec 1000 1000) (vec 0 1000)) 0 0 1))"),
        ("allowance-rounded", "(def f (polygon (vec (vec 0 0) (vec 500 -0.0000005) (vec 1000 0) (vec 1000 1000 0 20) (vec 0 1000 0 20)) 0 0 1))"),
    ] {
        let v = drafted_verdicts(tag, defs, &[("in", &i), ("out", &o)]).unwrap();
        assert!(v["in"] == "holds" && v["out"] == "fails", "{defs}: {v:?}");
    }
}

/// The verdicts of a model's claims, its text as given.
fn model_verdicts(tag: &str, text: &str) -> std::collections::BTreeMap<String, String> {
    let s = Scratch::new(tag);
    let m = s.model("m.ink", text);
    let d = Ink::new(["assert", &m]).env("INK_NO_CACHE", "1").any();
    assert_eq!(d["ok"].as_bool(), Some(true), "{text}: {d}");
    statuses(&d)
}

#[test]
fn a_drafted_change_s_height_is_one_number_whatever_the_span() {
    // One profile drafted over 0..h and over 0.015..h: its construction
    // changes near z = 0.0130 and 0.0171, the second span holding only the
    // second. Each height is the number its own change defines: named by
    // its place among the span's changes, the second span's read the
    // first's bounds, and a point 3e-5 inside it at z = 0.015 read 8.8e-4
    // outside at accuracy 0.
    let pts = "(vec (vec 0.029209085768858924 0.012671456905397274 0 0) (vec -8.581149304262623e-05 0.02682450207275581 0 0) (vec -0.006701061811055517 0.01742338768871102 0 0) (vec -0.012781097462207888 0.01327952064104351 0 0) (vec -0.020837383092469172 -0.003692253093367565 0 0) (vec -0.009568482445898521 -0.031182100045671863 0 0) (vec 0.03463492335699023 -0.015127474559352197 0 0) (vec 0.01884209069273787 -0.0020143964319454728 0 0)) (vec -20 0 -20 5 -20 -20 0 0)";
    let text = format!(
        "(def a (polygon {pts} 0 0.019479422232880274))\n(def b (polygon {pts} 0.015 0.019479422232880274))\n(assert a-inside (subset (probe -0.01775 0.00325 0.015) a))\n(assert b-inside (subset (probe -0.01775 0.00325 0.015) b))\n"
    );
    let v = model_verdicts("change-identity", &text);
    assert!(v["a-inside"] == "holds" && v["b-inside"] == "holds", "{v:?}");
}

#[test]
fn a_drafted_trim_s_height_holds_its_root_whatever_the_span() {
    // A chamfer drafted 45 degrees runs out at z = 0.70710678118654760...
    // just past the first span's top. Its bracket stored there was cut at
    // that top, excluding the root; the longer span read it, trimmed the
    // chamfer early, and (-1e-18, -1e-18) outside the chamfer at
    // z = 0.7071067811865476 read inside at accuracy 0.
    let pts = "(vec (vec -10 -10) (vec 0 -10) (vec 0 -1) (vec -1 0) (vec -3 0) (vec -3 -2) (vec -7 -2) (vec -7 0) (vec -10 0)) (vec 0 0 45 0 0 0 0 0 0)";
    let text = format!("(def a (polygon {pts} 0 0.7071067811865474))\n(def b (polygon {pts} 0 1))\n(assert interior (subset (probe -1e-18 -1e-18 0.7071067811865476) b))\n");
    let v = model_verdicts("trim-identity", &text);
    assert!(v["interior"] != "holds", "{v:?}");
}

#[test]
fn a_drafted_round_is_scaled_with_its_profile() {
    // A square 1e70 across, rounded 0.1 at the origin, drafted 1 degree:
    // built at a scale putting it near 1, its round read 1.4e69 (the
    // radius left unscaled), and (1, 1), past the round, outside at
    // accuracy 0. And 1e-70 across rounded 1e-71: (1e-73, 1e-73), in the
    // corner the round cuts away, read inside.
    let big =
        model_verdicts("round-scale-big", "(def f (polygon (vec (vec 0 0 0 0.1) (vec 1e70 0) (vec 1e70 1e70) (vec 0 1e70)) 1 0 1))\n(assert interior (subset (probe 1 1 0) f))\n");
    assert_eq!(big["interior"], "holds", "{big:?}");
    let small = model_verdicts(
        "round-scale-small",
        "(def f (polygon (vec (vec 0 0 0 1e-71) (vec 1e-70 0) (vec 1e-70 1e-70) (vec 0 1e-70)) (vec 0 45 0 0) 0 1e-70))\n(assert exterior (disjoint (probe 1e-73 1e-73 0) f))\n",
    );
    assert_eq!(small["exterior"], "holds", "{small:?}");
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

/// A value over ranged params left unknown says where its proof stopped
/// (the refinement limit: each ranged param halved 32 times) and where its
/// hunt for a violation did (its cap), never that no halving decides it:
/// narrower ranges prove it. It read "stays undecided however finely
/// halved".
#[test]
fn a_value_left_unknown_says_where_its_proof_and_its_hunt_stopped() {
    let s = Scratch::new("value-limits");
    let model = |name: &str, default: &str, tol: &str| {
        s.model(
            name,
            &format!(
                "(import \"std\")\n(param p {default} (tol {tol}))\n(param q {default} (tol {tol}))\n\
(assert near (value (max (+ (- p q) 1e-11) (+ (- q p) 1e-11)) > 0) (accuracy 1e-15))\n"
            ),
        )
    };
    let d = s.ink(["assert", &model("wide.ink", "0.5", "-0.5 0.5")]).ok();
    let a = &d["asserts"][0];
    let hint = a["hint"].as_str().unwrap_or("");
    assert_eq!(verdict(a), "unknown", "{}", pretty(a));
    assert!(hint.contains("refinement limit (32 halvings of each ranged param)") && !hint.contains("however finely"), "{hint}");
    assert!(a["hunt"]["cells"].as_u64().is_some_and(|n| n > 0) && a["hunt"]["stopped_by"] == "cap", "{}", pretty(a));
    assert!(hint.contains("without finding a certified violation, stopped at its cap"), "{hint}");
    // Ranges as narrow as the 32 halvings reach prove it.
    let d = s.ink(["assert", &model("narrow.ink", "0", "0 2.3283064365386963e-10")]).ok();
    assert_eq!(verdict(&d["asserts"][0]), "holds", "{}", pretty(&d));
}

/// `ink check` says what each claim will cost before it is asked, from
/// what it knows cheaply: a thickness on coarse cells and a value are under
/// a second on any machine, and a thickness on cells of a micron over a
/// 100 mm box over ten minutes on any; a proof's cost, which no cheap look
/// predicts, is counted by its budget and given no class.
#[test]
fn check_estimates_what_each_claim_will_cost() {
    let s = Scratch::new("cost");
    let m = s.model(
        "cost.ink",
        "(import \"std\")\n(body a (shape (box 100 100 10)))\n(body far (shape (move (sphere 5) 500 0 0)))\n\
         (param w 3)\n\
         (assert apart (disjoint (shape-of a) (shape-of far)))\n\
         (assert walls (min-thickness (shape-of a) > 1) (res 5))\n\
         (assert sized (value w > 1))\n\
         (assert fine-walls (min-thickness (shape-of a) > 1) (res 0.001))\n",
    );
    let d = s.ink(["check", &m]).ok();
    let cost = |n: &str| d["asserts"].as_array().unwrap().iter().find(|a| a["name"] == n).map(|a| a["cost_estimate"].clone()).unwrap_or_default();
    assert!(
        ["walls", "sized"].iter().all(|n| cost(n)["class"] == "under a second")
            && cost("fine-walls")["class"] == "over ten minutes"
            && cost("fine-walls")["basis"].as_str().is_some_and(|b| b.contains("surface samples"))
            && cost("apart")["class"].is_null()
            && cost("apart")["basis"].as_str().is_some_and(|b| b.contains("a proof")),
        "{}",
        pretty(&d["asserts"])
    );
}

/// A cost estimate reads only what `ink check` surveys already, the
/// bodies' surfaces: thicknesses of 30 shapes no body is (distinct spheres)
/// are not counted, and say so, rather than sampled by check.
#[test]
fn a_cost_estimate_samples_nothing_check_does_not() {
    let s = Scratch::new("cost-unsurveyed");
    let claims: String = (0..30).map(|i| format!("(assert w{i} (min-thickness (sphere {}) > 0.1))\n", 1.0 + 0.1 * i as f64)).collect();
    let d = s.ink(["check", &s.model("spheres.ink", &format!("(import \"std\")\n{claims}"))]).ok();
    let rows = d["asserts"].as_array().unwrap();
    assert!(
        rows.len() == 30 && rows.iter().all(|a| a["cost_estimate"]["class"].is_null() && a["cost_estimate"]["basis"].as_str().is_some_and(|b| b.contains("no body"))),
        "{}",
        pretty(&d["asserts"])
    );
}

/// A claim's estimate is of what `ink assert` will do with it: a value
/// over a ranged param, or over the motion, is a search (not classed), not
/// one evaluation; and a stoppable assert's estimates (its progress),
/// after check, are check's.
#[test]
fn a_cost_estimate_is_of_what_assert_does_and_assert_says_the_same() {
    let s = Scratch::new("cost-ranged");
    let m = s.model(
        "ranged.ink",
        "(import \"std\")\n(def w (ranged 1 4 3))\n(t-span 0 1)\n(body a (shape (box 10 10 10)))\n\
         (assert c (value w > 0))\n(assert moved (value (+ t 1) > 0))\n(assert plain (value 3 > 2))\n\
         (assert walls (min-thickness (shape-of a) > 1) (res 5))\n(assert apart (disjoint (shape-of a) (move (shape-of a) 50 0 0)))\n",
    );
    let cache = s.path("cache");
    let d = s.ink(["check", &m]).env("XDG_CACHE_HOME", &cache).env_remove("INK_NO_CACHE").ok();
    let cost = |n: &str| d["asserts"].as_array().unwrap().iter().find(|a| a["name"] == n).map(|a| a["cost_estimate"].clone()).unwrap_or_default();
    let searched = |n: &str| cost(n)["class"].is_null() && cost(n)["basis"].as_str().is_some_and(|b| b.starts_with("a search over"));
    let progress = s.path("progress.jsonl");
    s.ink(["assert", &m]).env("INK_PROGRESS", &progress).env("XDG_CACHE_HOME", &cache).env_remove("INK_NO_CACHE").any();
    let said: Value =
        std::fs::read_to_string(&progress).unwrap().lines().filter_map(|l| serde_json::from_str::<Value>(l).ok()).find(|l| l["estimates"].is_object()).unwrap_or_default();
    let same = ["c", "moved", "plain", "walls", "apart"].iter().all(|n| said["estimates"][n] == cost(n));
    assert!(searched("c") && searched("moved") && cost("plain")["class"] == "under a second" && same, "check: {}\nassert: {}", pretty(&d["asserts"]), pretty(&said));
}

/// Cases are counted without overflow, at most as many as a run asks: 64
/// text params with a case each are 2^64 combinations, which the run
/// refuses. Set by --param, a param's cases are not asked.
#[test]
fn many_cases_are_counted_without_overflow() {
    let s = Scratch::new("cost-cases");
    let params: String = (0..64).map(|i| format!("(param n{i} \"a\" (cases \"b\"))\n")).collect();
    let m = s.model("cases.ink", &format!("(import \"std\")\n{params}(assert one (value 1 > 0))\n"));
    let d = s.ink(["check", &m]).ok();
    let refused = d["asserts"][0]["cost_estimate"]["basis"].as_str().is_some_and(|b| b.contains("asked at most 257 times"));
    let e = s.ink(["assert", &m]).err();
    let mut pinned = vec!["check".to_string(), m.clone()];
    for i in 0..64 {
        pinned.extend(["--param".to_string(), format!("n{i}=a")]);
    }
    let p = s.ink(&pinned).ok();
    assert!(
        refused
            && e["error"]["message"].as_str().is_some_and(|m| m.contains("combinations"))
            && p["asserts"][0]["cost_estimate"]["basis"].as_str().is_some_and(|b| b.starts_with("one evaluation")),
        "{}\n{}\n{}",
        pretty(&d["asserts"]),
        pretty(&e),
        pretty(&p["asserts"])
    );
}

/// A stoppable assert (INK_PROGRESS) is never slowed by its estimates: on
/// a cold cache it computes no more than one without (counted as the
/// answers each stores: the same), and its thickness, whose body no check
/// has surveyed yet, is not classed; once `ink check` has surveyed it, the
/// estimate is check's.
#[test]
fn a_stoppable_assert_surveys_nothing_for_its_estimates() {
    let s = Scratch::new("cost-cold");
    let m = s.model("walls.ink", "(import \"std\")\n(body a (shape (box 40 20 10)))\n(assert walls (min-thickness (shape-of a) > 1) (res 2))\n");
    let stored = |cache: &str| std::fs::read_dir(std::path::Path::new(cache).join("ink/answers")).map_or(0, |d| d.count());
    let (plain, cold) = (s.path("plain"), s.path("cold"));
    s.ink(["assert", &m]).env("XDG_CACHE_HOME", &plain).env_remove("INK_NO_CACHE").ok();
    let estimate = |cache: &str, log: &str| {
        s.ink(["assert", &m]).env("XDG_CACHE_HOME", cache).env("INK_PROGRESS", s.path(log)).env_remove("INK_NO_CACHE").ok();
        let lines = std::fs::read_to_string(s.path(log)).unwrap();
        lines.lines().filter_map(|l| serde_json::from_str::<Value>(l).ok()).find(|l| l["estimates"].is_object()).map(|l| l["estimates"]["walls"].clone()).unwrap_or_default()
    };
    let before = estimate(&cold, "cold.jsonl");
    let same_work = stored(&plain) == stored(&cold);
    let checked = s.ink(["check", &m]).env("XDG_CACHE_HOME", &cold).env_remove("INK_NO_CACHE").ok();
    let after = estimate(&cold, "warm.jsonl");
    assert!(
        same_work
            && before["class"].is_null()
            && before["basis"].as_str().is_some_and(|b| b.contains("no survey yet; ink check computes it"))
            && after == checked["asserts"][0]["cost_estimate"]
            && after["class"].is_string(),
        "answers stored: {} without progress, {} with; cold: {before}; after check: {after}; check: {}",
        stored(&plain),
        stored(&cold),
        checked["asserts"][0]["cost_estimate"]
    );
}
