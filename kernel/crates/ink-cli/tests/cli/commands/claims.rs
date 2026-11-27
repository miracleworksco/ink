//! `ink require`: failing claims with their witness, unknowns with where
//! and what decides them, and which claims a file's own are (--own).

use super::*;
use serde_json::json;

#[test]
fn ranged_switches_cover_poles_and_large_integers() {
    let s = Scratch::new("ranged-switches");
    for (source, witness) in [
        ("(param p 0 (tol -2 2)) (require sign (value (- (select (tan p) 1 0) (select p 1 0)) <= 0))", true),
        ("(param p 9007199254740996 (tol -2 0)) (require bound (value (floor p) >= 9007199254740996))", false),
        (
            r#"(import "std") (param q 0 (tol 0 2))
            (def count (- (floor (+ 9007199254740994 q)) 9007199254740994))
            (require clear (disjoint (probe 0 0 0) (max (box 8 8 8) (- count 1))))"#,
            false,
        ),
    ] {
        let m = s.model("switch.ink", source);
        let d = ink(["require", &m]);
        let found = &d["requirements"][0];
        if witness {
            assert_eq!(verdict(found), "fails", "{d}");
        } else {
            assert_ne!(found["status"].as_str().unwrap(), "holds", "{d}");
        }
    }
}

#[test]
fn constructed_numbers_retain_their_bounds_in_claims() {
    let s = Scratch::new("real-coefficients");
    for source in [
        r#"(import "std") (t-span 0 1s)
            (def a (* 1e16 (/ t 1s))) (def travel (- (+ a (/ t 1s)) a))
            (body carriage (shape (box 2 2 2)))
            (assembly root (instance c carriage (at) (slide x travel) (speed-limits 0 0)))"#,
        "(t-span 0 1s) (def a (* 9007199254740992 (/ t 1s))) (require stationary (value (diff (- (+ a (/ t 1s)) a) t) <= 0))",
        "(t-span 0 1s) (def a (* 49 (/ (/ t 1s) 49))) (require slow (value (diff a t) < 1)) (require backwards (value (diff (- a (/ t 1s)) t) < 0))",
        "(param p -1) (require rising (value (diff (pow p 18014398509481984) p) > 0))",
        "(require rising (value (warp (diff (pow x 18014398509481984) x) -1 0 0) > 0))",
        "(param p 2 (tol 0 2)) (def n (/ 1 18014398509481984)) (require small (value (diff (pow p n) p) <= (/ n 2)))",
        "(def station (warp (floor (* 49 (/ x 49))) 1 0 0)) (require left (subset (probe station 0 0) (- x 0.5))) (require zero (value station <= 0)) (require threshold (value 0 >= station))",
        "(def h (nth (sin-cos-deg 30) 0)) (require below (value h < 0.5)) (require above (value 0.5 > h))",
    ] {
        let m = s.model("construction.ink", source);
        let d = ink(["require", &m]);
        let claims = d["requirements"].as_array().unwrap();
        assert!(!claims.is_empty(), "{d}");
        for claim in claims {
            assert_ne!(claim["status"].as_str().unwrap(), "holds", "{claim}");
        }
    }
}

#[test]
fn a_composed_translation_encloses_the_distance_between_faces() {
    let s = Scratch::new("composed-face-gap");
    let m = s.model(
        "gap.ink",
        "(body left (shape x)) (body right (shape (- 1.5 x)))
        (assembly root (instance a left (at (move 0.1 0 0) (move 0.7 0 0))) (instance b right))
        (require gap (clearance (shape-of root/a) (shape-of root/b) 0.70000000000000007))
        (require face-gap (value (warp (shape-of root/a) 1.5 0 0) >= 0.70000000000000007))",
    );
    let d = ink(["require", &m]);
    let claims = d["requirements"].as_array().unwrap();
    assert_eq!(claims.len(), 2);
    for claim in claims {
        assert_ne!(claim["status"].as_str().unwrap(), "holds", "{claim}");
    }
}

#[test]
fn exact_turns_and_nearly_parallel_frames_keep_coincident_bodies_together() {
    let s = Scratch::new("coincident-frames");
    let check = |body: &str, instances: &str, centre: &str, accuracy: &str| {
        let m = s.model(
            "placed.ink",
            &format!(
                r#"(import "std") {body}
            (assembly root {instances})
            (def a (shape-of root/one)) (def b (shape-of root/two))
            (require apart (disjoint a b) (accuracy {accuracy}))
            (require common-empty (empty (max a b)) (accuracy {accuracy}))
            (require contained (subset b (neg a)) (accuracy {accuracy}))
            (require gap (clearance a b {accuracy}))
            (require meet (overlaps a b) (accuracy {accuracy}))
            (require centre (subset (probe {centre}) b))"#
            ),
        );
        let d = ink(["require", &m]);
        let expected = ["fails", "fails", "fails", "fails", "holds", "holds"];
        let statuses: Vec<_> = d["requirements"].as_array().unwrap().iter().map(verdict).collect();
        assert_eq!(statuses, expected, "{d}");
        let d = ink(["interfere", &m]);
        assert_eq!(d["overlapping_pairs"], 1, "{d}");
    };
    check(
        "(def part (move (box 4 4 4) 20 0 0)) (body raw (shape part))
           (body turned (shape (rotate-z part 810647932926689280)))",
        "(instance one raw) (instance two turned)",
        "20 0 0",
        "0.05",
    );
    check(
        "(body fin (shape (move (box 0.0000152587890625 0.5 0.5) 0 0 20)) (attr tolerance 0.000001))
           (def good (axes (origin 0 0 0) (x 0 0 1) (z 1 1 0)))
           (def skew (axes (origin 0 0 0) (x 1 1 0.000000000116415321826934814453125) (z 1 1 0)))",
        "(instance one fin good) (instance two fin skew)",
        "14.14213562373 14.14213562373 0",
        "0.000001",
    );
    let m = s.model("quarter.ink", r#"(import "std") (require inside (subset (probe 0 -1 0) (rotate-z y 90)))"#);
    assert_eq!(verdict(&ink(["require", &m])["requirements"][0]), "fails");
}

#[test]
fn failing_claims_fail_an_overlap_with_its_witness() {
    let s = Scratch::new("failing");
    let failing = s.model(
        "failing.ink",
        "(import \"std\")\n(body a (shape (sphere 10)))\n(body b (shape (move (sphere 10) 15 0 0)))\n\
(require apart (disjoint (shape-of a) (shape-of b)))\n(require small (volume (shape-of a) < 1000))\n",
    );
    let d = ink(["require", &failing]);
    let (apart, small) = (&d["requirements"][0], &d["requirements"][1]);
    assert!(verdict(apart) == "fails" && apart.get("witness").is_some(), "the overlap fails with a witness: {}", pretty(apart));
    assert!(verdict(small) == "fails", "the volume claim fails, proven by its bounds (4189 against 1000): {}", pretty(small));
}

#[test]
fn a_failed_clearance_gives_its_gap_and_a_failed_overlap_its_depth() {
    // A 3 mm shank in a 6.2 mm hole: a true gap of 0.1, under the 0.2
    // claimed. The gap is the witness pair's distance, at least the true
    // one and under the clearance. Two balls of 10 at 15 apart overlap 5
    // deep: the witness lies -value into both (field values), and depth is
    // the overlap's, walked down from it to its middle, 2.5.
    let s = Scratch::new("gap");
    let m = s.model(
        "gap.ink",
        "(import \"std\")\n(def plate (subtract (box 40 40 5) (cylinder-z 3.1 10)))\n(def shank (move (cylinder-z 3 20) 0 0 -5.5))\n\
(require fit (clearance plate shank 0.2))\n(require apart (disjoint (sphere 10) (move (sphere 10) 15 0 0)))\n",
    );
    let d = s.ink(["require", &m]).any();
    let (fit, apart) = (&d["requirements"][0], &d["requirements"][1]);
    let (p, q) = (&fit["witness_pair"]["in_first"], &fit["witness_pair"]["in_second"]);
    let pair = (0..3).map(|k| (num(&p[k]) - num(&q[k])).powi(2)).sum::<f64>().sqrt();
    let gap = num(&fit["gap"]);
    assert!(verdict(fit) == "fails" && (gap - pair).abs() < 1e-6 && (0.1 - 1e-6..0.2).contains(&gap) && (gap - num(&fit["value"])).abs() < 1e-9, "{}", pretty(fit));
    let depth = num(&apart["depth"]);
    assert!(verdict(apart) == "fails" && (2.5 - 1e-6..=2.5).contains(&depth) && depth >= -num(&apart["value"]), "{}", pretty(apart));
}

#[test]
fn assert_exits_3_when_a_claim_does_not_hold_and_0_when_all_do() {
    // A script reads the exit status: 0 every claim holds, 3 answered with
    // one that fails or is unknown, 1 an error (no answer).
    let s = Scratch::new("exit");
    let m = s.model("exit.ink", "(import \"std\")\n(body a (shape (sphere 10)))\n(require big (volume (shape-of a) > 1000))\n(require small (volume (shape-of a) < 1000))\n");
    let code = |args: &[&str]| s.ink(args).run().code;
    assert_eq!(code(&["require", &m]), Some(3), "a failing claim");
    assert_eq!(code(&["require", &m, "--name", "big"]), Some(0), "every claim holds");
    assert_eq!(code(&["require", &m, "--name", "none"]), Some(1), "an error");
}

#[test]
fn an_unknown_exact_claim_gives_an_undecided_cell_and_a_res_hint() {
    // An unknown says where, as a failure does; one asked exactly (a range
    // at its bound on the face x = 0.5, the box reading x as (x + 3 y) -
    // 3 y, so its face and the bound share nothing to split at) is told
    // what a res can decide, not to state an accuracy it refuses.
    let s = Scratch::new("at-bound");
    let m = s.model("at-bound.ink", "(import \"std\")\n(require at-bound (range x -1 0.5 (warp (box 1 1 1) (- (+ x (* 3 y)) (* 3 y)) y z)))\n");
    let d = ink(["require", &m]);
    let a = &d["requirements"][0];
    let c = &a["undecided_cell"];
    let hint = a["hint"].as_str().unwrap_or("");
    assert!(verdict(a) == "unknown" && (num(&c["centre"][0]) - 0.5).abs() <= num(&c["size"][0]) && hint.contains("(res r)") && !hint.contains("accuracy"), "{}", pretty(a));
}

#[test]
fn own_runs_only_the_files_claims_and_check_imports_says_where_each_is_written() {
    // --own checks the claims written in the file, not the ones it imports
    // (here a false one), and says so. The model rule relies on it.
    let s = Scratch::new("own");
    let base = s.model("own-base.ink", "(import \"std\")\n(body b (shape (box 10 10 10)))\n(require too-big (volume (shape-of b) > 2000))\n");
    let top = s.model("own-top.ink", "(import \"own-base.ink\" as base)\n(require fits (within (shape-of base.b) (-6 -6 -6) (6 6 6)))\n");
    let d = ink(["check", &top, "--imports"]);
    let own: Vec<(String, bool, String)> =
        d["requirements"].as_array().unwrap().iter().map(|a| (a["name"].as_str().unwrap().into(), a["own"] == true, a["source"].as_str().unwrap().into())).collect();
    assert!(own == [("base.too-big".into(), false, format!("{base}:3")), ("fits".into(), true, format!("{top}:2"))], "{own:?}");
    // With --imports, the claims its imports require too.
    let d = Ink::new(["require", &top, "--imports"]).any();
    let names: Vec<&str> = d["requirements"].as_array().unwrap().iter().map(|a| a["name"].as_str().unwrap()).collect();
    assert!(names == ["base.too-big", "fits"] && d["all_yes"] == false && d.get("imported_not_run").is_none(), "--imports asks imported claims too: {}", pretty(&d));
    let d = ink(["require", &top]);
    let rows: Vec<(String, String)> = d["requirements"].as_array().unwrap().iter().map(|a| (a["name"].as_str().unwrap().into(), verdict(a))).collect();
    assert!(rows == [("fits".to_string(), "holds".to_string())] && d["all_proven"] == true && d["imported_not_run"] == 1, "{}", pretty(&d));
    let e = &ink_err(["require", &top, "--name", "base.too-big"])["error"];
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
    let d = ink(["require", &m, "--name", "seated"]);
    assert!(verdict(&d["requirements"][0]) == "fails", "{}", pretty(&d["requirements"]));
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
(require rests (touches (shape-of root/bearing) (max (shape-of root/housing) (- z ledge))))\n",
    );
    let d = ink(["require", &m, "--name", "rests"]);
    let v: Vec<String> = d["requirements"].as_array().unwrap().iter().map(verdict).collect();
    assert!(v == ["fails"], "{}", pretty(&d["requirements"]));
}

#[test]
fn an_overlapping_clearance_fails_with_value_and_gap_zero() {
    // The pair's distance is not rounded up per coordinate where the two
    // points are one: that gives 1.5e-323 beside a gap of 0.
    let s = Scratch::new("clearance-overlap");
    let m = s.model(
        "overlap.ink",
        "(import \"std\")\n(body a (shape (box 10 10 10)))\n(body b (shape (move (box 10 10 10) 4 0 0)))\n\
        (require clear (clearance (shape-of a) (shape-of b) 0.3))\n",
    );
    let d = &ink(["require", &m])["requirements"][0];
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
            &format!("(import \"std\")\n(require in (subset (probe {x} 0 0) (box 10 10 10)))\n(require out (disjoint (probe {x} 0 0) (box 10 10 10)))\n"),
        );
        let d = ink(["require", &m]);
        (verdict(&d["requirements"][0]), verdict(&d["requirements"][1]), d["requirements"][0]["accuracy"].clone())
    };
    for (x, want) in [("5.0001", ("fails", "holds")), ("4.9999", ("holds", "fails")), ("5.000000001", ("fails", "holds")), ("5", ("fails", "holds"))] {
        let (i, o, acc) = rows(x);
        assert!(i == want.0 && o == want.1 && acc == 0.0, "probe at x = {x}: subset {i}, disjoint {o}, accuracy {acc}");
    }
}

#[test]
fn a_file_with_no_claims_says_so_and_guarantees_nothing() {
    // It said all_proven: true, vacuously.
    let s = Scratch::new("no-claims");
    let m = s.model("none.ink", "(import \"std\")\n(body b (shape (sphere 2)))\n");
    let d = ink(["require", &m]);
    assert!(d["all_yes"].is_null() && d["all_proven"].is_null() && d["nothing_required"].is_string(), "{}", pretty(&d));
}

/// Two polygons' allowances in one claim, each asked at its own worst end:
/// a box cut by polygon A lies in polygon B, its top `gap` below B's edge,
/// each allowance about 1.1e-5. A grown and B shrunk meet at 1.5e-5 (not
/// settled); every allowance at one end, low or high, leaves the gap open.
/// At 3e-5 it holds.
#[test]
fn a_claim_reading_two_allowances_asks_each_at_its_own_worst_end() {
    let s = Scratch::new("two-allowances");
    let model = |gap: f64| {
        let top = 40.0 - gap;
        format!(
            "(import \"std\")\n(def b (polygon (vec (vec 0 0 0) (vec 100 0 1) (vec 100 40 0) (vec 0 40 1))))\n(def a (max (move (box 1 1 1) 50 39.5 0) (polygon (vec (vec 0 {lo} 0) (vec 100 {lo} 1) (vec 100 {top} 0) (vec 0 {top} 1)))))\n(require in (subset a b) (accuracy 1e-6))\n",
            lo = -gap
        )
    };
    let verdict = |gap: f64| {
        let m = s.model(&format!("m{gap}.ink"), &model(gap));
        let d = Ink::new(["require", &m]).env("INK_NO_CACHE", "1").any();
        statuses(&d)["in"].clone()
    };
    assert_eq!(verdict(1.5e-5), "unknown");
    assert_eq!(verdict(3e-5), "holds");
}

/// A box deep in a polygon built within its allowance, moved by a ranged
/// param, fails: the point found in it is read over the allowance at its
/// own value of the param, not over the param's whole range, where the box
/// leaves the polygon. It was not settled.
#[test]
fn a_moving_box_deep_in_a_polygon_built_within_its_allowance_fails() {
    let s = Scratch::new("moving-in-allowance");
    let m = s.model(
        "m.ink",
        "(import \"std\")\n(param s 0 (tol 0 1))\n(def slot (max (polygon (vec (vec 0 0 0) (vec 100 0 1) (vec 100 40 0) (vec 0 40 1))) (- (abs z) 1)))\n(require moving-in (subset (move (box 1 1 1) (+ 30 (* s 100)) 20 0) (neg slot)))\n",
    );
    let d = Ink::new(["require", &m]).env("INK_NO_CACHE", "1").any();
    assert_eq!(statuses(&d)["moving-in"], "fails", "{}", pretty(&d));
}

/// A clearance between a polygon built within its allowance and a block
/// 0.1 away fails at 0.11 and holds at 0.09: proven both ways, the shapes
/// read at their smallest where it fails. It was not settled at 0.11 (a
/// failing clearance had no point to find over the allowance).
#[test]
fn a_clearance_to_a_polygon_built_within_its_allowance_fails_where_it_is_short() {
    let s = Scratch::new("clearance-allowance");
    let m = s.model(
        "m.ink",
        "(import \"std\")\n(def p (max (polygon (vec (vec 0 0 0) (vec 10 0 0) (vec 10 4 0) (vec 0 4 0))) (- (abs z) 1)))\n(def q (move (box 5 5 2) 12.6 2 0))\n(require short (clearance p q 0.11))\n(require clear (clearance p q 0.09))\n",
    );
    let d = Ink::new(["require", &m]).env("INK_NO_CACHE", "1").any();
    let st = statuses(&d);
    assert!(st["short"] == "fails" && st["clear"] == "holds", "{}", pretty(&d));
}

/// A record's number and the same number written are one number in a value
/// claim: a record's thread depth 3.2 against a written 3.2 is a proven yes
/// at <= 0 (it was a proven no: the record's the nearest double, the written
/// the decimal, 1e-16 apart), and 3.25 a proven no.
#[test]
fn a_records_number_is_the_number_its_digits_write() {
    let s = Scratch::new("record-digits");
    s.model("rec.json", r#"{"depth": 3.2, "pitch": 0.1}"#);
    let m = s.model(
        "m.ink",
        "(import \"std\")\n(def r (read-json \"rec.json\"))\n(require same (value (abs (- (get r \"depth\") 3.2)) <= 0))\n(require sum (value (+ (get r \"pitch\") 0.2) <= 0.3))\n(require other (value (abs (- (get r \"depth\") 3.25)) <= 0))\n",
    );
    let d = Ink::new(["require", &m]).env("INK_NO_CACHE", "1").any();
    let st = statuses(&d);
    assert!(st["same"] == "holds" && st["sum"] == "holds" && st["other"] == "fails", "{}", pretty(&d));
}

/// A drafted polygon's claims at a point: each claim's verdict, or why the
/// polygon is refused.
fn drafted_verdicts(tag: &str, defs: &str, claims: &[(&str, &str)]) -> Result<std::collections::BTreeMap<String, String>, String> {
    let s = Scratch::new(tag);
    let mut text = format!("(import \"std\")\n{defs}\n(body part (shape (max f (slab-z -1 1))))\n");
    for (name, claim) in claims {
        text += &format!("(require {name} {claim})\n");
    }
    let m = s.model("m.ink", &text);
    let d = Ink::new(["require", &m]).env("INK_NO_CACHE", "1").any();
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

/// A straight polygon is built within its allowance of its drawn corners
/// (1.4e-4 for a 1414 mm square): a point 2.5e-7 inside a corner 5e-7
/// below the line through its neighbours, and one as far outside, are not
/// settled at accuracy 0, never read the wrong way.
#[test]
fn a_straight_polygon_is_its_drawn_corners() {
    let defs = "(def f (polygon (vec (vec 0 0) (vec 500 -0.0000005) (vec 1000 0) (vec 1000 1000) (vec 0 1000))))";
    let (pin, pout) = ("(probe 500 -2.4999999999999999e-07 0.5)", "(probe 500 -7.5e-07 0.5)");
    let (i, o) = (format!("(subset {pin} f)"), format!("(disjoint {pout} f)"));
    let v = drafted_verdicts("straight-exact", defs, &[("in", &i), ("out", &o)]).unwrap();
    assert_eq!((v["in"].as_str(), v["out"].as_str()), ("unknown", "unknown"), "{v:?}");
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
    // 1e-7 of the 1414 mm diagonal) is within the profile's allowance: the
    // point 2.5e-7 inside the drawn profile is not settled either way at
    // accuracy 0, never read outside. So too with its far corners rounded.
    let p = "(probe 500 -2.4999999999999999e-07 0.5)";
    let (i, o) = (format!("(subset {p} f)"), format!("(disjoint {p} f)"));
    for (tag, defs) in [
        ("allowance-corner", "(def f (polygon (vec (vec 0 0) (vec 500 -0.0000005) (vec 1000 0) (vec 1000 1000) (vec 0 1000)) 0 0 1))"),
        ("allowance-rounded", "(def f (polygon (vec (vec 0 0) (vec 500 -0.0000005) (vec 1000 0) (vec 1000 1000 0 20) (vec 0 1000 0 20)) 0 0 1))"),
    ] {
        let v = drafted_verdicts(tag, defs, &[("in", &i), ("out", &o)]).unwrap();
        assert!(v["in"] == "unknown" && v["out"] == "unknown", "{defs}: {v:?}");
    }
}

/// A drafted section that changes within its span is refused, saying where:
/// one whose construction changes near
/// z = 0.0130 (two corners' walls crossing), a chamfer drafted 45 degrees
/// that runs out at z = 0.7071, and a square 1e70 across rounded 0.1,
/// whose round vanishes under its allowance. polygon's construction holds
/// one construction over its span; the exact drafted builder, which cut
/// its span where the construction changed, is gone.
#[test]
fn a_drafted_section_that_changes_within_its_span_is_refused() {
    let pts = "(vec (vec 0.029209085768858924 0.012671456905397274 0 0) (vec -8.581149304262623e-05 0.02682450207275581 0 0) (vec -0.006701061811055517 0.01742338768871102 0 0) (vec -0.012781097462207888 0.01327952064104351 0 0) (vec -0.020837383092469172 -0.003692253093367565 0 0) (vec -0.009568482445898521 -0.031182100045671863 0 0) (vec 0.03463492335699023 -0.015127474559352197 0 0) (vec 0.01884209069273787 -0.0020143964319454728 0 0)) (vec -20 0 -20 5 -20 -20 0 0)";
    let trim = "(vec (vec -10 -10) (vec 0 -10) (vec 0 -1) (vec -1 0) (vec -3 0) (vec -3 -2) (vec -7 -2) (vec -7 0) (vec -10 0)) (vec 0 0 45 0 0 0 0 0 0)";
    for (tag, text, said) in [
        ("change", format!("(def a (polygon {pts} 0 0.019479422232880274))\n(require a (subset (probe -0.01775 0.00325 0.015) a))\n"), "walls cross, at z = 0.013"),
        ("trim", format!("(def b (polygon {trim} 0 1))\n(require b (subset (probe -5 -5 0.5) b))\n"), "walls cross, at z = 0.7071"),
        (
            "round",
            "(def f (polygon (vec (vec 0 0 0 0.1) (vec 1e70 0) (vec 1e70 1e70) (vec 0 1e70)) 1 0 1))\n(require f (subset (probe 1 1 0) f))\n".into(),
            "round at (0, 0) (point 1) vanishes",
        ),
    ] {
        let s = Scratch::new(&format!("drafted-refused-{tag}"));
        let m = s.model("m.ink", &format!("(import \"std\")\n{text}"));
        let d = Ink::new(["require", &m]).env("INK_NO_CACHE", "1").any();
        let e = d["error"]["message"].as_str().unwrap_or("");
        assert!(d["ok"] == false && e.contains(said), "{tag}: {d}");
    }
}

/// A promise is taken only for its constants as computed for the query:
/// at --param p=1e200 (fixed, or a range of one value) k = 0.001 + q / q
/// is 0.001 (q overflows, inf / inf is 0), and the promise taken as
/// written would guarantee a clearance between points 0.042 apart.
#[test]
fn a_promise_is_not_taken_where_its_params_overflow() {
    let f = fixture("promise-param-overflow.ink");
    assert_eq!(statuses(&ink(["require", &f, "--param", "p=1"]))["gap"], "fails");
    // 1e200 is outside ink's range: unknown, never the false holds.
    for p in ["p=1e200", "p=1e200:1e200", "p=1:1e200"] {
        let s = statuses(&ink(["require", &f, "--param", p]));
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
(require near (value (max (+ (- p q) 1e-11) (+ (- q p) 1e-11)) > 0) (accuracy 1e-15))\n"
            ),
        )
    };
    let d = s.ink(["require", &model("wide.ink", "0.5", "-0.5 0.5")]).ok();
    let a = &d["requirements"][0];
    let hint = a["hint"].as_str().unwrap_or("");
    assert_eq!(verdict(a), "unknown", "{}", pretty(a));
    assert!(hint.contains("refinement limit (32 halvings of each ranged param)") && !hint.contains("however finely"), "{hint}");
    assert!(a["hunt"]["cells"].as_u64().is_some_and(|n| n > 0) && a["hunt"]["stopped_by"] == "cap", "{}", pretty(a));
    assert!(hint.contains("without finding a certified violation, stopped at its cap"), "{hint}");
    // Ranges as narrow as the 32 halvings reach prove it.
    let d = s.ink(["require", &model("narrow.ink", "0", "0 2.3283064365386963e-10")]).ok();
    assert_eq!(verdict(&d["requirements"][0]), "holds", "{}", pretty(&d));
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
         (require apart (disjoint (shape-of a) (shape-of far)))\n\
         (require walls (min-thickness (shape-of a) > 1) (res 5))\n\
         (require sized (value w > 1))\n\
         (require fine-walls (min-thickness (shape-of a) > 1) (res 0.001))\n",
    );
    let d = s.ink(["check", &m]).ok();
    let cost = |n: &str| d["requirements"].as_array().unwrap().iter().find(|a| a["name"] == n).map(|a| a["cost_estimate"].clone()).unwrap_or_default();
    assert!(
        ["walls", "sized"].iter().all(|n| cost(n)["class"] == "under a second")
            && cost("fine-walls")["class"] == "over ten minutes"
            && cost("fine-walls")["basis"].as_str().is_some_and(|b| b.contains("surface samples"))
            && cost("apart")["class"].is_null()
            && cost("apart")["basis"].as_str().is_some_and(|b| b.contains("a proof")),
        "{}",
        pretty(&d["requirements"])
    );
}

/// A cost estimate reads only what `ink check` surveys already, the
/// bodies' surfaces: thicknesses of 30 shapes no body is (distinct spheres)
/// are not counted, and say so, rather than sampled by check.
#[test]
fn a_cost_estimate_samples_nothing_check_does_not() {
    let s = Scratch::new("cost-unsurveyed");
    let claims: String = (0..30).map(|i| format!("(require w{i} (min-thickness (sphere {}) > 0.1))\n", 1.0 + 0.1 * i as f64)).collect();
    let d = s.ink(["check", &s.model("spheres.ink", &format!("(import \"std\")\n{claims}"))]).ok();
    let rows = d["requirements"].as_array().unwrap();
    assert!(
        rows.len() == 30 && rows.iter().all(|a| a["cost_estimate"]["class"].is_null() && a["cost_estimate"]["basis"].as_str().is_some_and(|b| b.contains("no body"))),
        "{}",
        pretty(&d["requirements"])
    );
}

/// A claim's estimate is of what `ink require` will do with it: a value
/// over a ranged param, or over the motion, is a search (not classed), not
/// one evaluation; and a stoppable assert's estimates (its progress),
/// after check, are check's.
#[test]
fn a_cost_estimate_is_of_what_assert_does_and_assert_says_the_same() {
    let s = Scratch::new("cost-ranged");
    let m = s.model(
        "ranged.ink",
        "(import \"std\")\n(def w (ranged 1 4 3))\n(t-span 0 1s)\n(body a (shape (box 10 10 10)))\n\
         (require c (value w > 0))\n(require moved (value (+ (/ t 1s) 1) > 0))\n(require plain (value 3 > 2))\n\
         (require walls (min-thickness (shape-of a) > 1) (res 5))\n(require apart (disjoint (shape-of a) (move (shape-of a) 50 0 0)))\n",
    );
    let cache = s.path("cache");
    let d = s.ink(["check", &m]).env("XDG_CACHE_HOME", &cache).env_remove("INK_NO_CACHE").ok();
    let cost = |n: &str| d["requirements"].as_array().unwrap().iter().find(|a| a["name"] == n).map(|a| a["cost_estimate"].clone()).unwrap_or_default();
    let searched = |n: &str| cost(n)["class"].is_null() && cost(n)["basis"].as_str().is_some_and(|b| b.starts_with("a search over"));
    let progress = s.path("progress.jsonl");
    s.ink(["require", &m]).env("INK_PROGRESS", &progress).env("XDG_CACHE_HOME", &cache).env_remove("INK_NO_CACHE").any();
    let said: Value =
        std::fs::read_to_string(&progress).unwrap().lines().filter_map(|l| serde_json::from_str::<Value>(l).ok()).find(|l| l["estimates"].is_object()).unwrap_or_default();
    let same = ["c", "moved", "plain", "walls", "apart"].iter().all(|n| said["estimates"][n] == cost(n));
    assert!(searched("c") && searched("moved") && cost("plain")["class"] == "under a second" && same, "check: {}\nassert: {}", pretty(&d["requirements"]), pretty(&said));
}

/// Cases are counted without overflow, at most as many as a run asks: 64
/// text params with a case each are 2^64 combinations, which the run
/// refuses. Set by --param, a param's cases are not asked.
#[test]
fn many_cases_are_counted_without_overflow() {
    let s = Scratch::new("cost-cases");
    let params: String = (0..64).map(|i| format!("(param n{i} \"a\" (cases \"b\"))\n")).collect();
    let m = s.model("cases.ink", &format!("(import \"std\")\n{params}(require one (value 1 > 0))\n"));
    let d = s.ink(["check", &m]).ok();
    let refused = d["requirements"][0]["cost_estimate"]["basis"].as_str().is_some_and(|b| b.contains("asked at most 257 times"));
    let e = s.ink(["require", &m]).err();
    let mut pinned = vec!["check".to_string(), m.clone()];
    for i in 0..64 {
        pinned.extend(["--param".to_string(), format!("n{i}=a")]);
    }
    let p = s.ink(&pinned).ok();
    assert!(
        refused
            && e["error"]["message"].as_str().is_some_and(|m| m.contains("combinations"))
            && p["requirements"][0]["cost_estimate"]["basis"].as_str().is_some_and(|b| b.starts_with("one evaluation")),
        "{}\n{}\n{}",
        pretty(&d["requirements"]),
        pretty(&e),
        pretty(&p["requirements"])
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
    let m = s.model("walls.ink", "(import \"std\")\n(body a (shape (box 40 20 10)))\n(require walls (min-thickness (shape-of a) > 1) (res 2))\n");
    let stored = |cache: &str| std::fs::read_dir(std::path::Path::new(cache).join("ink/answers")).map_or(0, |d| d.count());
    let (plain, cold) = (s.path("plain"), s.path("cold"));
    s.ink(["require", &m]).env("XDG_CACHE_HOME", &plain).env_remove("INK_NO_CACHE").ok();
    let estimate = |cache: &str, log: &str| {
        s.ink(["require", &m]).env("XDG_CACHE_HOME", cache).env("INK_PROGRESS", s.path(log)).env_remove("INK_NO_CACHE").ok();
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
            && after == checked["requirements"][0]["cost_estimate"]
            && after["class"].is_string(),
        "answers stored: {} without progress, {} with; cold: {before}; after check: {after}; check: {}",
        stored(&plain),
        stored(&cold),
        checked["requirements"][0]["cost_estimate"]
    );
}

#[test]
fn declared_case_starts_at_the_runs_moment_after_an_imported_span() {
    let s = Scratch::new("case-clock");
    s.model("clock.ink", "(t-span 5s 6s) (require clock-ok (value t >= 5s))");
    let m = s.model("main.ink", r#"(param offset 0s (cases 1s)) (require safe (value t >= offset)) (import "clock.ink" as clock)"#);
    let d = s.ink(["require", &m, "--imports"]).ok();
    let safe = &d["requirements"][0];
    let case = &safe["cases"][0];
    assert_eq!(case["one_moment"], 0.0);
    assert_eq!(case["status"], "fails");
    assert_eq!(safe["status"], "fails");
    assert_eq!(d["all_yes"], false);
    assert_eq!(d["all_proven"], false);
    let clock = &d["requirements"][1]["cases"][0];
    assert_eq!(clock["status"], "holds");
    assert_eq!(clock["motion"]["span"], json!([5.0, 6.0]));
}

#[test]
fn an_unloaded_case_is_incomplete_even_when_no_claim_exists_at_default() {
    let s = Scratch::new("empty-at-default");
    let f = s.model("conditional.ink", "(param n 0 (cases 1)) (for i (range n) (require limit-$i (value 1 < 0)))");
    let d = s.ink(["require", &f]).env("INK_DEADLINE", "1000000").ok();
    assert_eq!(d["cases_not_loaded"], json!([{"n": {"case": 1}}]));
    assert_eq!(d["all_yes"], false);
    assert_eq!(d["all_proven"], false);
    assert!(d.get("no_claims").is_none());
}

#[test]
fn optional_joint_limits_need_their_package_in_every_placement() {
    let s = Scratch::new("optional-limit");
    s.model("part.ink", "(body b (shape 1)) (assembly part (instance j b (at) (slide x 0) (optional (effort-limits -2N 2N))))");
    let m = s.model(
        "main.ink",
        r#"(import "part.ink" as p) (param n 1 (cases 2))
        (assembly root (for i (range n) (instance m-$i p.part))) (require local (value 1 > 0))"#,
    );
    for flags in [vec![], vec!["--imports"]] {
        let run = s.ink([&["require", &m][..], &flags].concat()).run();
        assert_eq!(run.code, Some(3));
        let d = run.json();
        assert_eq!(d["all_proven"], false);
        assert_eq!(d["all_yes"], false);
        assert_eq!(d["not_checked"], json!(["effort-limits:root/m-0/j", "effort-limits:root/m-1/j"]));
        assert!(d["not_checked_note"].as_str().unwrap().contains("no package"));
    }
}

#[test]
fn physical_density_checks_values_and_the_whole_requested_band() {
    let s = Scratch::new("density-values");
    let cube = "(body weight (shape (box 10 10 10)) (attr density rho))\n(require ballast (integral weight (attr density) > 1000g))";
    let m = s.model("cube.ink", &format!("(import \"std\") (param rho 2700kg/m3) {cube}"));
    let good = s.ink(["require", &m]).any();
    let claim = &good["requirements"][0];
    assert_eq!(claim["status"], "fails");
    assert!((num(&claim["value"]) - 2.7).abs() < 1e-6);
    for cmd in ["check", "mass", "require"] {
        let d = s.ink([cmd, &m, "--param", "rho=2700"]).err();
        assert!(
            d["error"]["message"].as_str().unwrap().contains("rho is written 2700kg/m3") && d["error"]["message"].as_str().unwrap().contains("rho=2700 a unit"),
            "{}",
            pretty(&d)
        );
    }
    let bad = s.model("bad.ink", &format!("(import \"std\") (param rho 2700g/mm3) {cube}"));
    assert!(says_density(&s.ink(["require", &bad]).err()));
    s.model("material.json", r#"{"rho":2700}"#);
    s.model("material.ink", "(fn identity (x0) x0) (def rho (identity (get (read-json \"material.json\") \"rho\")))");
    let data = s.model("data.ink", &format!("(import \"std\") (import \"material.ink\" as m) (def rho (* 1g/mm3 m.rho)) {cube}"));
    assert!(says_density(&s.ink(["require", &data]).err()));

    // Density in the hollow centre of a part is never used.
    let shell =
        |r| s.model("shell.ink", &format!("(import \"std\") (body shell (shape (max (sphere 5) (- (sphere 3)))) (attr density (* 1g/mm3 (select (sphere {r}) 2700 0.0027))))"));
    s.ink(["check", &shell(2)]).ok();
    assert!(says_density(&s.ink(["check", &shell(4)]).err()));

    let graded = s.model("fade.ink", "(import \"std\") (body part (shape (box 60 10 10)) (attr density (* 0.0012g/mm3 (clamp (/ (- 20 x) 10) 0 1))))");
    s.ink(["check", &graded]).ok();
    let band = s.model("infill.ink", &format!("(import \"std\") (param infill 0.2 (tol -0.2 0.8)) (def rho (* 0.0012g/mm3 infill)) {cube}"));
    assert_eq!(s.ink(["require", &band]).any()["requirements"][0]["status"], "fails");

    let band = s.model("band.ink", &format!("(import \"std\") (param rho 0.002g/mm3 (tol 0 1g/mm3)) {cube} (require width (value 8 < 7))"));
    let d = s.ink(["require", &band]).any();
    assert_eq!(d["requirements"][0]["status"], "unknown");
    let why = d["requirements"][0]["unavailable"].as_str().unwrap();
    assert!(["rho=0.002..1.002", "density (g/mm^3)", "0.00000005..0.023", "unit"].iter().all(|t| why.contains(t)), "{why}");
    assert_eq!(d["requirements"][1]["status"], "fails");
    // The CLI loads the midpoint before installing the requested band.
    // Its bad density may only stop the claim that reads the density.
    let d = s.ink(["require", &band, "--param", "rho=0.01g/mm3:0.09g/mm3"]).any();
    assert_eq!(d["requirements"][0]["status"], "unknown");
    assert!(d["requirements"][0]["unavailable"].as_str().unwrap().contains("rho=0.01..0.09"));
    assert_eq!(d["requirements"][1]["status"], "fails");

    // Density reaches a plausible value only outside the body.
    let m = s.model("peak.ink", r#"(import "std") (body b (shape (box 2 2 2)) (attr density (* 1g/mm3 (select (box 3 3 3) 1e-12 0.01))))"#);
    assert!(says_density(&s.ink(["check", &m]).err()));
    let m = s.model("named.ink", "(import \"std\") (param rho 2700) (body part (shape (box 10 10 10)) (attr gross rho))");
    let d = s.ink(["mass", &m, "--density", "gross"]).err();
    let why = d["error"]["message"].as_str().unwrap();
    assert!(["--density gross", "part", "rho=2700", "0.00000005..0.023", "2700kg/m3"].iter().all(|t| why.contains(t)), "{why}");
}

fn says_density(d: &serde_json::Value) -> bool {
    d["error"]["message"].as_str().is_some_and(|m| m.contains("density") && m.contains("0.00000005..0.023") && (m.contains("2700") || m.contains("field bounds")))
}

#[test]
fn ranged_mean_keeps_searching_after_an_undefined_sample() {
    let s = Scratch::new("mean-samples");
    for (bound, status, count) in [(10, "fails", 2), (30, "unknown", 17)] {
        let m = s.model("mean.ink", &format!("(import \"std\") (param p 1 (tol -1 1)) (require a (mean (box 2 2 2) 20 (weight p) < {bound}))"));
        let d = s.ink(["require", &m]).any();
        let row = &d["requirements"][0];
        assert_eq!(row["status"], status);
        assert_eq!(row["sampled_values"], count);
        if status == "fails" {
            assert!((num(&row["value"]) - 20.0).abs() < 1e-6 && num(&row["params"]["p"]) > 0.0, "{}", pretty(row));
        }
    }
}

#[test]
fn a_declared_case_with_an_invalid_claim_scope_is_refused() {
    let s = Scratch::new("overflow-case");
    let m = s.model("scope.ink", "(param width 0 (tol 0 1e308) (cases 1e308)) (require control (value 1 > 0))");
    let d = s.ink(["require", &m]).err();
    let error = d["error"]["message"].as_str().unwrap();
    assert!(error.contains("at the cases") && error.contains("width") && error.contains("overflows"), "{d}");
}

#[test]
fn a_shallow_crossing_warns_and_a_half_tolerance_shortfall_is_allowed() {
    let s = Scratch::new("path-separation");
    for (pts, radius, warns) in [
        ("(vec 0 0 0) (vec 1000000 0 0) (vec 1000000 -10 0) (vec 0 -10 0) (vec 0 -0.16 0) (vec 1000000 0.64 0)", 0.05, true),
        ("(vec 0 0 0) (vec 10 0 0) (vec 10 -10 0) (vec 5 -10 1.98) (vec 5 10 1.98)", 1.0, false),
        ("(vec 0 0 0) (vec 10 0 0) (vec 10 -10 0) (vec 5 -10 1.90) (vec 5 10 1.90)", 1.0, true),
    ] {
        let m = s.model("path.ink", &format!("(import \"std\") (def routed (tube (vec {pts}) {radius}))"));
        let d = s.ink(["check", &m]).ok();
        assert_eq!(pretty(&d["warnings"]).contains("runs through itself"), warns, "{}", pretty(&d));
    }
}

#[test]
fn an_empty_line_has_zero_lettering_bounds() {
    let s = Scratch::new("empty-lettering");
    let m = s.model("empty.ink", "(import \"std\") (def size (lettering-max \"A\" 0 (font \"sans\") 4 0)) (require width (value (nth size 0) >= 0))");
    let d = s.ink(["require", &m]).ok();
    assert_eq!(d["requirements"][0]["status"], "holds");
    assert_eq!(num(&d["requirements"][0]["value"]), 0.0);
}

#[test]
fn a_moved_density_is_checked_by_the_claim_that_weighs_it() {
    let s = Scratch::new("moved-density");
    let model = |density: &str| {
        s.model(
            "part.ink",
            &format!(
                r#"(import "std")
        (param gap 1 (tol -0.1 0.1))
        (body slab (shape (box 40 10 10)) (attr density (* {density} (clamp (/ (- 20 x) 10) 0 1))))
        (assembly root (instance s slab (at (move 50 0 0))))
        (require ballast (integral root (attr density) > 1kg))
        (require independent (value 1 < 0))
        (require moved (integral root (move (attr-of slab density) 5 0 0) > 1kg))
        (require turned (integral root (rotate-z (attr-of slab density) 90) > 1kg))"#
            ),
        )
    };
    let d = s.ink(["require", &model("2.7g/mm3")]).any();
    assert_eq!(d["requirements"][0]["status"], "unknown");
    assert!(d["requirements"][0]["unavailable"].as_str().unwrap().contains("0.00000005..0.023"));
    assert_eq!(d["requirements"][1]["status"], "fails");
    // A density the claim itself moves or turns is the same density: checked as where it was written (it
    // answered yes at 10 kg, a thousand times the part).
    for k in [2, 3] {
        assert_eq!(d["requirements"][k]["status"], "unknown", "{}", pretty(&d["requirements"][k]));
        assert!(d["requirements"][k]["unavailable"].as_str().unwrap().contains("0.00000005..0.023"));
    }
    assert_eq!(s.ink(["require", &model("2.7g/cm3")]).any()["requirements"][0]["status"], "fails");
}

#[test]
fn unit_defaults_require_units_at_each_authored_input() {
    let s = Scratch::new("param-input-units");
    for default in ["200GPa", "(* 200 1GPa)", "(+ 200GPa 0)", "(+ 0MPa 200GPa)", "(/ (* 200GPa 3mm) 3mm)"] {
        let part = |options: &str, more: &str| s.model("part.ink", &format!("(param modulus {default} {options}) (require stiff (value modulus > 1000MPa)) {more}"));
        let m = part("", "");
        for bare in ["modulus=210", "modulus=-10:10"] {
            let d = s.ink(["require", &m, "--param", bare]).err();
            assert!(d["error"]["message"].as_str().unwrap().contains(&format!("{bare} a unit")), "{d}");
        }
        assert_eq!(s.ink(["require", &m, "--param", "modulus=0"]).any()["requirements"][0]["status"], "fails");
        assert_eq!(s.ink(["require", &m, "--param", "modulus=210GPa"]).any()["requirements"][0]["status"], "holds");
        for (opts, more, why) in [("(tol -10 10GPa)", "", "tol bound -10"), ("(tol -10GPa 10)", "", "tol bound 10")] {
            let d = s.ink(["require", &part(opts, more)]).err();
            let msg = d["error"]["message"].as_str().unwrap();
            assert!(msg.contains("modulus") && msg.contains(why) && msg.contains("unit"), "{msg}");
        }
        let valid = part("(tol 0 10GPa)", "");
        assert_eq!(s.ink(["require", &valid]).any()["requirements"][0]["status"], "holds");
        let file = part("(cases 210)", "");
        assert!(s.ink(["require", &file]).err()["error"]["message"].as_str().unwrap().contains("210 is a bare number: say its unit"));
        part("", "");
        let importer = s.model("import.ink", r#"(import "part.ink" (modulus 210))"#);
        assert!(s.ink(["check", &importer]).err()["error"]["message"].as_str().unwrap().contains("import value 210 a unit"));
        let importer = s.model("import.ink", r#"(import "part.ink" (modulus 0 (tol 0 10)))"#);
        assert!(s.ink(["check", &importer]).err()["error"]["message"].as_str().unwrap().contains("tol bound 10 a unit"));
    }
    let m = s.model("native.ink", "(param wall 2mm (tol -1 1)) (param preload 5N) (require size (value wall > 2)) (require unloaded (value preload < 1N))");
    let d = s.ink(["require", &m, "--nominal", "--param", "wall=3", "--param", "preload=0"]).any();
    assert_eq!(d["counts"]["yes"], 2);
}

#[test]
fn a_density_peak_inside_an_unbounded_part_is_a_witness() {
    let s = Scratch::new("unbounded-density");
    let model = |field| s.model("slab.ink", &format!("(body slab (shape (- x 1)) (attr density {field}))"));
    s.ink(["check", &model("(* 1g/mm3 (select y 0.001 0))")]).ok();
    for field in ["(* 1g/mm3 (select (- x 2) 1e-12 0.001))", "(* 1g/mm3 (select (- x -2000) 2.7 0.001))"] {
        assert!(says_density(&s.ink(["check", &model(field)]).err()));
    }
}
