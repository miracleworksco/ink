//! Assemblies through the command line: `check`'s answer (hidden
//! constructions, undefined cases), analysis over all of space (bounded by a search that cannot miss
//! material, not clipped to the window), contacts and relations, and std's
//! gear meshes and their claims.

use crate::common::*;
use serde_json::Value;
use std::collections::BTreeMap;

fn sc(tag: &str) -> Scratch {
    Scratch::new(&format!("asm-{tag}"))
}

fn run(s: &Scratch, args: &[&str]) -> Value {
    s.ink(args).ok()
}

fn refused(s: &Scratch, args: &[&str]) -> Value {
    s.ink(args).err()
}

fn msg(d: &Value) -> &str {
    d["error"]["message"].as_str().unwrap_or("")
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"))
}

fn nums(v: &Value) -> Vec<f64> {
    v.as_array().map_or(vec![], |a| a.iter().map(f).collect())
}

fn strs(v: &Value) -> Vec<String> {
    v.as_array().map_or(vec![], |a| a.iter().map(|x| x.as_str().unwrap_or("?").to_string()).collect())
}

/// An answer's warnings containing `needle`.
fn warnings(d: &Value, needle: &str) -> Vec<String> {
    strs(&d["warnings"]).into_iter().filter(|w| w.contains(needle)).collect()
}

/// `ink assert ARGS`: each claim's verdict by the last part of its name.
fn short(s: &Scratch, args: &[&str]) -> BTreeMap<String, String> {
    let mut a = vec!["assert"];
    a.extend_from_slice(args);
    statuses(&run(s, &a)).into_iter().map(|(k, v)| (k.rsplit('/').next().unwrap().to_string(), v)).collect()
}

/// The claims of `ink assert ARGS` that do not hold, by short name.
fn failed(s: &Scratch, args: &[&str]) -> Vec<String> {
    short(s, args).into_iter().filter(|(_, v)| v != "holds").map(|(k, _)| k).collect()
}

fn map<const N: usize>(kv: [(&str, &str); N]) -> BTreeMap<String, String> {
    kv.into_iter().map(|(k, v)| (k.to_string(), v.to_string())).collect()
}

/// Each pair of an interfere answer, `A~B` by the parts' last names, and
/// its verdict.
fn pairs(d: &Value) -> BTreeMap<String, String> {
    d["pairs"].as_array().unwrap().iter().map(|p| (strs(&p["pair"]).iter().map(|q| q.rsplit('/').next().unwrap().to_string()).collect::<Vec<_>>().join("~"), verdict(p))).collect()
}

/// Each pair of an interfere answer by its full paths.
fn full_pairs(d: &Value) -> BTreeMap<String, String> {
    d["pairs"].as_array().unwrap().iter().map(|p| (strs(&p["pair"]).join("~"), verdict(p))).collect()
}

// ---------------------------------------------------------------- check

#[test]
fn check_flags_a_construction_hidden_inside_another_not_a_narrow_groove() {
    let s = sc("hidden");
    let hid = s.model(
        "hidden.ink",
        "(import \"std\")\n(def block (box 20 20 20))\n(def boss (move (sphere 3) 2 0 0))\n\
        (def hole (move (cylinder-z 2 40) 5 0 0))\n(body b (shape (subtract (min block boss) hole)))\n",
    );
    let w = warnings(&run(&s, &["check", &hid]), "forms none");
    assert!(w.len() == 1 && w[0].contains("'boss'"), "{w:?}");
    let groove = s.model(
        "groove.ink",
        "(import \"std\")\n(body bar (shape (subtract (move (box 400 12 8) 200 0 4)\n\
        \x20 (tag groove (move (cylinder-x 1.1 402) 200 6.45 5.75)))))\n",
    );
    let w = warnings(&run(&s, &["check", &groove]), "forms none");
    assert!(w.is_empty(), "{w:?}");
}

#[test]
fn undefined_cases_are_reported_by_check_and_by_eval_by_row() {
    // A ratio whose denominator is 0 wherever x <= 0 gives 0 there; only
    // cases the value depends on count (a guarded or losing branch is not
    // reported).
    let s = sc("undefined");
    let m = s.model("undefined.ink", "(import \"std\")\n(body b (shape (max (box 12 12 12) (+ (box 10 10 10) (tag ratio (/ y (max x 0)))))))\n");
    let w = warnings(&run(&s, &["check", &m]), "division by zero");
    assert!(w.len() == 1 && w[0].contains("'ratio'"), "{w:?}");
    let pts = "[[0,1,0],[1,1,0]]";
    // Picked keys keep it: a value that passed through an undefined case
    // says so however little is asked for, by row.
    let t = &run(&s, &["eval", &m, "--points", pts, "--keys", "value"])["targets"][0]["results"];
    let und = t["undefined"].as_array().cloned().unwrap_or_default();
    assert!(
        strs(&t["columns"]) == ["value"]
            && t["rows"].as_array().unwrap().len() == 2
            && und.len() == 1
            && und[0][0] == 0
            && und[0][1][0].as_str().unwrap_or("").contains("division by zero"),
        "{}",
        pretty(t)
    );
}

#[test]
fn a_rotation_axis_zero_within_a_tolerance_band_is_refused_and_warned_of() {
    // The rule: crates/ink/tests/geometry.rs a_direction_without_a_length_is_refused.
    let s = sc("zero-axis");
    let m = s.model(
        "zero-axis.ink",
        "(import \"std\")\n(param a 1 (tol -2 2))\n(body ball (shape (sphere 1)))\n\
        (assembly root (instance b ball (at (rot a 0 0 60))))\n(assert x (value 1 > 0))\n",
    );
    let w = warnings(&run(&s, &["check", &m]), "has zero length at a = 0 (squared length 0), within a in [-1, 3]");
    let e = refused(&s, &["assert", &m]);
    assert!(w.len() == 1 && msg(&e).contains("(assembly root ...)"), "{w:?}\n{}", pretty(&e));
}

#[test]
fn a_solved_value_makes_claims_measured_through_shape_of_placements_and_solves() {
    let s = sc("labels");
    let m = s.model(
        "labels.ink",
        "(import \"std\")\n(solve size (unknowns (r0 1)) (== (* r0 r0 r0) 8))\n(def r r0)\n\
        (body derived (shape (sphere r)))\n(body plain (shape (sphere 1)))\n\
        (assembly root (instance placed plain (at (move r 0 0))))\n(solve fit (unknowns (u 1)) (== u r))\n\
        (assert via-body (within (shape-of derived) (-3 -3 -3) (3 3 3)))\n\
        (assert via-placement (within (shape-of root) (-4 -4 -4) (4 4 4)))\n\
        (assert via-solve (value u > 1))\n(assert exact (within (shape-of plain) (-2 -2 -2) (2 2 2)))\n",
    );
    let r = short(&s, &[&m]);
    assert_eq!(r, map([("exact", "holds"), ("via-body", "holds estimate"), ("via-placement", "holds estimate"), ("via-solve", "holds estimate")]));
}

// ------------------------------------------ joints, contacts, relations

#[test]
fn check_lists_a_joint_and_a_limits_witness_gives_the_drive_and_limits() {
    let s = sc("joint");
    let m = s.model(
        "joint.ink",
        "(import \"std\")\n(param s 100)\n(body rail (shape (move (box 270 12 8) 135 0 4)) (frame top (at (move 0 0 8))))\n\
        (body car (shape (move (box 45 27 10) 0 0 5)) (frame base (at (rot-x 180))))\n\
        (assembly root (instance r rail) (instance c car (frame-of r top) base (slide x s) (limits 25 245)))\n",
    );
    let j = &run(&s, &["check", &m])["assemblies"][0]["instances"][1]["joint"];
    assert!(
        j["kind"] == "slide"
            && j["axis"] == "x"
            && j["drive"] == "s"
            && strs(&j["params"]) == ["s"]
            && j["reads_t"] == false
            && nums(&j["limits"]) == [25.0, 245.0]
            && j["claim"] == "limits:root/c",
        "{}",
        pretty(j)
    );
    // The witness says the drive and limits there, not the margin alone.
    let w = s.model("joint-t.ink", "(import \"std\")\n(body car (shape (box 2 4 6)))\n(assembly root (instance c car (at (move 0 0 0)) (turn z (* 40 t)) (limits -10 30)))\n");
    let a = &run(&s, &["assert", &w, "--t-span", "0:1"])["asserts"][0];
    let tm = &a["terms"];
    assert!(
        verdict(a) == "fails" && (f(&tm["drive"]) - 40.0 * f(&a["t"])).abs() < 1e-9 && f(&tm["drive"]) > 30.0 && f(&tm["lo"]) == -10.0 && f(&tm["hi"]) == 30.0,
        "{}",
        pretty(a)
    );
}

#[test]
fn an_imported_assemblys_params_are_listed_and_its_limits_are_the_driving_files_claims() {
    // Listed with the values each was built with, and its joints' claims
    // stated where it is driven.
    let s = sc("driven");
    s.model(
        "motor.ink",
        "(import \"std\")\n(body housing (shape (box 4 4 4)))\n(body rotor (shape (box 1 2 1)))\n\
        (assembly (motor (angle 0)) (instance h housing) (instance r rotor (at) (turn z angle) (speed-limits 0 6000rpm)))\n",
    );
    let w = s.model("driven.ink", "(import \"std\")\n(import \"motor.ink\" as m)\n(t-span 0 1)\n(assembly root (instance a (m.motor (* 9000rpm t))))\n");
    let d = run(&s, &["check", &w, "--imports"]);
    let asm: BTreeMap<&str, &Value> = d["assemblies"].as_array().unwrap().iter().map(|a| (a["name"].as_str().unwrap(), a)).collect();
    assert!(
        asm["m.motor"]["params"] == serde_json::json!([{"name": "angle", "value": "0"}])
            && asm["root/a"]["of"] == "m.motor"
            && asm["root/a"]["params"] == serde_json::json!([{"name": "angle", "value": "(* 9000rpm t)"}])
            && nums(&asm["root/a"]["instances"][1]["joint"]["speed_limits"]) == [0.0, 36000.0],
        "{}",
        pretty(&d["assemblies"])
    );
    let d = run(&s, &["assert", &w]);
    let a = d["asserts"].as_array().unwrap().iter().find(|a| a["name"] == "speed-limits:root/a/r").unwrap_or_else(|| panic!("{}", pretty(&d)));
    assert!(verdict(a) == "fails" && f(&a["terms"]["rate"]) == 54000.0 && a["source"].as_str().unwrap().ends_with("driven.ink:4"), "{}", pretty(a));
}

#[test]
fn a_contact_zone_excuses_the_thread_not_the_head() {
    // A contact zone excuses the thread in the nut, not the head's 5 mm
    // clash with the plate.
    let s = sc("contact");
    let parts = "(import \"std\")\n(body screw (shape (min (tag thread (move (cylinder-z 5 30) 0 0 15)) (move (cylinder-z 8 5) 0 0 -2.5))))\n\
        (body nut (shape (subtract (box 16 16 8) (cylinder-z 4.2 20))))\n(body plate (shape (subtract (box 40 40 10) (cylinder-z 5.5 20))))\n\
        (assembly root (instance s screw) (instance n nut (at (move 0 0 20))) (instance p plate)";
    let d = run(&s, &["interfere", &s.model("contact-none.ink", &format!("{parts})\n"))]);
    assert_eq!(pairs(&d), map([("s~n", "fails"), ("s~p", "fails")]));
    let d = run(&s, &["interfere", &s.model("contact.ink", &format!("{parts} (contact s n (shape-of screw thread)))\n"))]);
    let st = pairs(&d);
    assert!(st["s~n"] != "fails" && st["s~p"] == "fails" && d["overlapping_pairs"] == 1, "{st:?}");
}

/// Two blocks overlapping 3 mm, their contact declared by the assemblies
/// placing them: `root` for `carried` (turned and moved), `mid` for `p`
/// (tilted), and `root` again for `m/p` through `mid`'s placement. A zone
/// placed anywhere but where its ends are misses the overlap.
pub const ANCESTOR_CONTACT: &str = "(import \"std\")\n(body block (shape (box 10 10 10)))\n\
    (assembly pair (instance a block) (instance b block (at (move 7 0 0))))\n\
    (assembly mid (instance p pair (at (rot-x 10))) (instance q pair (at (move 0 50 0))) (contact p/a p/b (box 20 20 20)))\n\
    (assembly root (instance carried pair (at (rot-z 23) (move 100 200 0))) (instance m mid (at (rot-y 15) (move -300 0 0)))\n\
    \x20 (contact carried/a carried/b (move (box 8 20 20) 3.5 0 0)) (contact m/q/a m/q/b (move (box 8 20 20) 3.5 0 0)))\n";

#[test]
fn a_contact_declared_above_a_target_holds_for_it() {
    // Every pair is excused by the contact a parent declares, whichever
    // assembly the target names: narrowing it narrows the pairs, not the
    // contacts between them.
    let s = sc("ancestor-contact");
    let w = s.model("ancestor-contact.ink", ANCESTOR_CONTACT);
    let whole = run(&s, &["interfere", &w]);
    let zoned = |d: &Value| d["pairs"].as_array().unwrap().iter().all(|p| verdict(p) == "holds" && p["contact_zones"].as_array().is_some_and(|z| z.len() == 1));
    assert!(whole["overlapping_pairs"] == 0 && whole["pairs"].as_array().unwrap().len() == 3 && zoned(&whole), "{}", pretty(&whole));
    for (target, n) in [("root/carried", 1), ("root/m", 2), ("root/m/p", 1), ("root/m/q", 1)] {
        let d = run(&s, &["interfere", &w, "--target", target]);
        assert!(d["overlapping_pairs"] == 0 && d["pairs"].as_array().unwrap().len() == n && zoned(&d), "{target}: {}", pretty(&d));
        // The same pairs, by the same paths, as the whole root's.
        for (k, v) in full_pairs(&d) {
            assert_eq!(full_pairs(&whole).get(&k), Some(&v), "{target}: {k}");
        }
    }
    // Asked by its definition's name, the pair has no parent: its clash is
    // its own.
    let d = run(&s, &["interfere", &w, "--target", "pair"]);
    assert!(d["overlapping_pairs"] == 1, "{}", pretty(&d));
}

#[test]
fn interfere_without_a_root_is_refused_once_an_assembly_places_a_body() {
    // Unplaced, every body sits at its own origin: refused, own or imported.
    // Bodies alone, or beside a library of frames, are checked as written.
    let s = sc("no-root");
    let bodies = "(import \"std\")\n(body a (shape (box 1 1 1)))\n(body b (shape (move (box 1 1 1) 5 0 0)))\n";
    let pair = s.model("pair.ink", &format!("{bodies}(assembly pair (instance a a) (instance b b))\n"));
    for m in [pair, s.model("uses.ink", "(import \"pair.ink\" as p)\n")] {
        assert!(msg(&refused(&s, &["interfere", &m])).contains("--target NAME"), "{m}");
    }
    s.model("lib.ink", "(import \"std\")\n(assembly datum (frame top (at (move 0 0 10))))\n");
    for m in [s.model("bodies.ink", bodies), s.model("frames.ink", &format!("{bodies}(import \"lib.ink\" as lib)\n"))] {
        let d = run(&s, &["interfere", &m]);
        assert!(d["overlapping_pairs"] == 0 && d["pairs_apart_by_bounding_boxes"] == 1, "{m}: {}", pretty(&d));
    }
}

#[test]
fn checks_survey_plans_only_what_interfere_would_check() {
    // What interfere asks at one moment: pair's, placed, not its bodies
    // unplaced; fit's, not its template's at the defaults.
    let s = sc("no-root-survey");
    let bodies = "(import \"std\")\n(body a (shape (box 1 1 1)))\n(body b (shape (move (box 1 1 1) t 0 0)))\n";
    for (asm, said) in [
        ("(assembly pair (instance a a) (instance b b (at (move 5 0 0))))", "(instance pair/b)"),
        ("(assembly (unit (d 0)) (instance a a) (instance b b (at (move d 0 0))))\n(assembly fit (instance u (unit 5)))", "(instance fit/u/b)"),
    ] {
        let d = run(&s, &["check", &s.model("survey.ink", &format!("{bodies}{asm}\n"))]);
        assert!(d["warnings"].as_array().unwrap().iter().any(|w| w.as_str().unwrap().contains(&format!("t = 0 {said}"))), "{asm}: {}", pretty(&d));
    }
}

#[test]
fn a_contact_zone_not_proven_bounded_is_an_error() {
    // It would hide every clash.
    let s = sc("contact-open");
    let parts = "(import \"std\")\n(body b (shape (box 1 1 1)))\n(assembly root (instance a b) (instance b b) (contact a b ";
    for (i, zone) in ["-1", "(- z 2)"].into_iter().enumerate() {
        let e = refused(&s, &["interfere", &s.model(&format!("open-{i}.ink"), &format!("{parts}{zone}))\n"))]);
        assert!(msg(&e).contains("not proven bounded"), "{zone}: {}", pretty(&e));
    }
    // Bounded at the nominal values, but everything (or open in x) over part
    // of a param band or time span (else it would prove two coincident boxes
    // apart).
    let two = "(import \"std\")\n(param p 1 (tol -1 0))\n(body b (shape (box 1 1 1)))\n(assembly root (instance a b) (instance b b) (contact a b ";
    let cases: [(&str, &str, &[&str]); 4] = [
        ("param", "(select (- p .5) -1 (box 2 2 2))", &[][..]),
        ("param-warp", "(warp (box 2 2 2) (* p x) y z)", &[]),
        ("time", "(select (- .5 t) -1 (box 2 2 2))", &["--t-span", "0:1"]),
        ("time-warp", "(warp (box 2 2 2) (* (- 1 t) x) y z)", &["--t-span", "0:1"]),
    ];
    for (name, zone, scope) in cases {
        let m = s.model(&format!("contact-{name}.ink"), &format!("{two}{zone}))\n"));
        let mut a = vec!["interfere", m.as_str()];
        a.extend_from_slice(scope);
        let e = refused(&s, &a);
        assert!(msg(&e).contains("not proven bounded over the motion and param ranges"), "{name}: {}", pretty(&e));
    }
}

/// A pinion and wheel placed in mesh by mesh-turn and mesh-distance, the
/// assembly left open for a relate clause.
const GEARS: &str = "(import \"std\")\n(param gap 0)\n(param lift 0)\n(param tilt 0)\n(param hub 4)\n(param turn 0)\n(param phase 0)\n\
    (fn gear (n mo h) (subtract (extrude-z (involute-gear-2d n mo 20 0.1) 0 h) (cylinder-z 2 20)))\n\
    (body pinion (shape (min (gear 15 1 8) (cylinder-z hub 30))))\n(body wheel (shape (gear 45 1 8)))\n\
    (assembly root (instance p pinion (at (rot-z turn)))\n\
    \x20 (instance w wheel (at (rot-z (+ phase (mesh-turn 15 45 0 turn))) (rot-y tilt) (move (+ gap (mesh-distance 15 45 1)) 0 lift)))";

fn related() -> String {
    format!("{GEARS}\n  (relate p w (gear-mesh 15 45 1 20 0.2)))\n")
}

#[test]
fn a_related_gear_pairs_teeth_are_checked_not_zoned() {
    // std's gear-mesh has no zone (gears proven, not trusted), so the teeth
    // are checked with the rest of the pair. The wheel is in phase at the
    // moment checked; half a pitch off, its teeth clash.
    let s = sc("mesh");
    let m = s.model("mesh.ink", &related());
    let d = run(&s, &["interfere", &m, "--param", "turn=7"]);
    let p = &d["pairs"][0];
    let rel = &p["relations"][0];
    assert!(
        d["overlapping_pairs"] == 0
            && rel["kind"] == "gear-mesh"
            && rel.get("skipped").is_none()
            && nums(&rel["data"]["teeth"]) == [15.0, 45.0]
            && p.get("contact_zones").is_none(),
        "{}",
        pretty(&d["pairs"])
    );
    for extra in [&["--param", "phase=4"][..], &["--param", "hub=9"]] {
        // Out of phase the teeth clash; a hub clashing outside the mesh is still found.
        let mut a = vec!["interfere", m.as_str(), "--param", "turn=7"];
        a.extend_from_slice(extra);
        let d = run(&s, &a);
        assert!(d["overlapping_pairs"] == 1 && verdict(&d["pairs"][0]) == "fails", "{extra:?}: {}", pretty(&d["pairs"]));
    }
    // mesh-turn phases the wheel: unrelated, the teeth hold, and half a pitch off they clash.
    let free = s.model("mesh-free.ink", &format!("{GEARS})\n"));
    let d: Vec<Value> = [&[][..], &["--param", "phase=4"]]
        .iter()
        .map(|p| {
            let mut a = vec!["interfere", free.as_str(), "--param", "turn=7"];
            a.extend_from_slice(p);
            run(&s, &a)["overlapping_pairs"].clone()
        })
        .collect();
    assert_eq!(d, [0, 1]);
}

#[test]
fn gear_mesh_claims_fail_too_close_tilted_or_past_tip_clearance() {
    let s = sc("mesh-claims");
    let m = s.model("mesh-claims.ink", &related());
    let r = short(&s, &[&m]);
    assert!(r.len() == 3 && r.values().all(|v| v == "holds"), "{r:?}");
    assert_eq!(failed(&s, &[&m, "--param", "gap=-0.6"]), ["centre-distance"], "0.6 mm too close");
    assert_eq!(failed(&s, &[&m, "--param", "tilt=1"]), ["parallel-axes"], "tilted 1 degree, past what the backlash allows");
    // With teeth and module alone the backlash is unknown: the gears are
    // within the quarter module of tip clearance.
    let plain = s.model("mesh-plain.ink", &format!("{GEARS}\n  (relate p w (gear-mesh 15 45 1)))\n"));
    let r: Vec<String> = ["0", "0.2", "-0.2", "0.3"].iter().map(|g| short(&s, &[&plain, "--param", &format!("gap={g}")])["centre-distance"].clone()).collect();
    assert_eq!(r, ["holds", "holds", "holds", "fails"]);
}

#[test]
fn a_gear_on_a_datum_turned_half_round_has_exactly_parallel_axes() {
    // Met face to face with a datum turned (rot-x 180): its axis is the
    // wheel's to rounding, so parallel-axes with no backlash holds (not
    // 1.2e-16 off, sin(pi), failing).
    let s = sc("mesh-flipped");
    let m = s.model(
        "mesh-flipped.ink",
        "(import \"std\")\n(fn gear (n h) (subtract (extrude-z (involute-gear-2d n 1 20 0.1) 0 h) (cylinder-z 2 20)))\n\
        (body pinion (shape (gear 15 8)) (frame top (at (move 0 0 8))))\n(body wheel (shape (gear 45 8)))\n\
        (body datum (shape (box 4 4 2)) (frame under (at (rot-x 180) (move 0 0 -1))))\n\
        (assembly root (instance d datum (at (move 0 0 9))) (instance p pinion (frame-of d under) top)\n\
        \x20 (instance w wheel (at (rot-z (mesh-turn 15 45 0 0)) (move (mesh-distance 15 45 1) 0 0)))\n\
        \x20 (relate p w (gear-mesh 15 45 1 20 0)))\n",
    );
    let r = short(&s, &[&m]);
    assert_eq!(r.get("parallel-axes").map(String::as_str), Some("holds"), "{r:?}");
}

#[test]
fn check_lists_a_relation_and_refuses_clauses_and_zones_it_cannot_take() {
    let s = sc("relate");
    let a = &run(&s, &["check", &s.model("mesh-check.ink", &related())])["assemblies"][0]["relations"][0];
    assert!(
        strs(&a["pair"]) == ["p", "w"]
            && a["kind"] == "gear-mesh"
            && f(&a["data"]["module"]) == 1.0
            && strs(&a["claims"]) == ["relate:root/p~w/centre-distance", "relate:root/p~w/parallel-axes", "relate:root/p~w/faces-meet"],
        "{}",
        pretty(a)
    );
    // A relation is a zone, claims and data: a point probe is an assert's
    // (subset (probe ...) shape), not a clause of a kind.
    let e = refused(
        &s,
        &[
            "check",
            &s.model(
                "relate-probe.ink",
                "(import \"std\")\n(body cube (shape (box 2 2 2)))\n(assembly root (instance a cube) (instance b cube (at (move 5 0 0)))\n\
        \x20 (relate a b (lambda (fa fb sa sb) (relation probed (material inside a 0 0 0)))))\n",
            ),
        ],
    );
    assert!(msg(&e).contains("bad relation clause") && !e["error"].to_string().contains("(material"), "{}", pretty(&e));
    let e = refused(
        &s,
        &[
            "interfere",
            &s.model(
                "relate-open.ink",
                "(import \"std\")\n(body b (shape (box 1 1 1)))\n(assembly root (instance a b) (instance b b)\n\
        \x20 (relate a b (lambda (fa fb sa sb) (relation open (zone \"all\" -1 -1)))))\n",
            ),
        ],
    );
    assert!(msg(&e).contains("not proven bounded"), "a relation's zone, as a contact's: {}", pretty(&e));
    // An empty assembly's shape is the empty shape.
    let d = run(
        &s,
        &[
            "check",
            &s.model(
                "relate-empty.ink",
                "(import \"std\")\n(assembly empty)\n(body cube (shape (box 1 1 1)))\n\
        (assembly root (instance e empty) (instance c cube) (relate e c (lambda (fa fb sa sb) (relation empty))))\n",
            ),
        ],
    );
    assert!(d["assemblies"][1]["relations"][0]["kind"] == "empty", "{}", pretty(&d["assemblies"]));
}

#[test]
fn a_zone_equal_to_another_items_field_leaves_that_items_pair_alone() {
    // The pair moved by a joint takes each body from its own instance, not
    // from a map by field (which proves a/c, apart by nothing, with the b/c zone).
    let s = sc("zone-alias");
    let d = run(
        &s,
        &[
            "interfere",
            &s.model(
                "zone-alias.ink",
                "(import \"std\")\n(t-span 0 1)\n(body unit (shape (box 2 2 2)))\n(body far (shape (move (sphere 1) 100 0 0)))\n\
        (assembly root (instance a unit (at (move 10 0 0))) (instance b far) (instance c unit (at) (slide x (+ 9 t)))\n\
        \x20 (relate b c (lambda (fa fb sa sb) (relation unrelated (zone \"b/c only\" (place (shape-of unit) (at (move 10 0 0))))))))\n",
            ),
        ],
    );
    assert!(d["overlapping_pairs"] == 1 && verdict(&d["pairs"][0]) == "fails", "{}", pretty(&d["pairs"]));
}

// ----------------------------------------- analysis covers all of space
//
// A shape is bounded by a search that cannot miss material (query::extent),
// and one that is not proven bounded is examined in the window, its results
// that need all of it unknown. Clipping to the window would give false
// proofs and false disproofs on these models.

/// Spheres at the origin and 1000 mm out, in bodies and an instance.
const FAR: &str = "(import \"std\")\n(body a (shape (min (sphere 1) (move (sphere 1) 1000 0 0))))\n\
    (body b (shape (min (move (sphere 1) 10 0 0) (move (sphere 1) 1003 0 0))))\n\
    (body unit (shape (sphere 1)))\n(assembly root (instance far unit (at (move 1000 0 0))))\n\
    (assert gone (empty (shape-of root)))\n(assert near (empty (sphere 1)))\n";

#[test]
fn distance_finds_pieces_far_outside_the_window_and_a_press_fit_overlaps() {
    let s = sc("distance");
    let far = s.model("beyond.ink", FAR);
    let d = &run(&s, &["distance", &far, "--a", "a", "--b", "b", "--res", "0.1"])["distance"];
    assert!(0.85 < f(&d["lower_proven"]) && f(&d["lower_proven"]) <= 1.0 + 1e-9 && f(&d["upper_witness"]) >= 1.0 - 1e-6, "gap 1: {}", pretty(d));
    // A witness needs no accuracy: a 0.04 press fit (under the 0.05) has
    // points inside both, so it overlaps; faces that only touch are unknown
    // (the fit said false).
    let press = s.model(
        "press.ink",
        "(import \"std\")\n(body a (shape (box 2 2 2)))\n(body b (shape (move (box 2 2 2) 1.96 0 0)))\n\
        (body c (shape (move (box 2 2 2) 2 0 0)))\n(assembly root (instance a a) (instance b b) (instance c c))\n",
    );
    let o: Vec<Value> = ["root/b", "root/c"].iter().map(|x| run(&s, &["distance", &press, "--a", "root/a", "--b", x])["overlapping"].clone()).collect();
    assert_eq!(o, [Value::Bool(true), Value::Null]);
}

/// A robot arm's rear fork, six turn joints out (demos/assembly's arm at
/// t 4.3125), against a wall: the distance's search halves a box to its
/// finest undecided, and the seek that ends it there finds a closer pair,
/// so the search goes on at the nearer distance. It panicked on halving
/// that stuck box; now it answers, the box undecided.
const STUCK: &str = r#"(import "std")
(body fork (shape (min (move (extrude-z (rounded-rect-2d 45 10 2) -38 -35.2) 22.5 0 0) (emax (slab-x 33 45) (slab-y -5 5) (slab-z -38 9.2)))))
(body wall (shape (move (box 6 150 52) 0 0 26)))
(def keys-t (vec 0 1.4 2.2 3.4 6.4 7.8 8.6 9.8 12))
(fn track (ks)
  (+ (nth ks 0) (fold + 0 (map (lambda (i) (* (- (nth ks (+ i 1)) (nth ks i))
    (let ((u (clamp (/ (- t (nth keys-t i)) (- (nth keys-t (+ i 1)) (nth keys-t i))) 0 1))) (* u u (- 3 (* 2 u))))))
    (range 0 8)))))
(def shoulder (track (vec 92 70 70 98 98 74 74 98 92)))
(def elbow (track (vec -88 -80 -80 -92 -92 -82 -82 -92 -88)))
(assembly s6 (instance fork fork (at (rot-z -45) (rot-y 90))))
(assembly s5 (instance s6 s6 (at (move 70 0 6)) (turn x (track (vec 0 0 0 0 90 0 0 0 0)))))
(assembly s4 (instance s5 s5 (at (move 120 0 6)) (turn z (+ 90 shoulder elbow))))
(assembly s3 (instance s4 s4 (at (move 130 0 6)) (turn z (neg elbow))))
(assembly s2 (instance s3 s3 (at (rot-x -90) (move 0 22 45)) (turn z (neg shoulder))))
(assembly s1 (instance s2 s2 (at (move 0 0 41.45)) (turn z (track (vec -40 -40 -40 -40 62 62 62 62 -40)))))
(assembly root (instance s1 s1) (instance wall wall (at (rot-z -80) (move 167.4 29.5 0))))
"#;

#[test]
fn distance_past_a_stuck_box_a_witness_ended_answers_with_it_undecided() {
    let s = sc("stuck");
    let m = s.model("stuck.ink", STUCK);
    let d = &run(&s, &["distance", &m, "--t", "4.3125", "--a", "root/s1/s2/s3/s4/s5/s6/fork", "--b", "root/wall"])["distance"];
    let (lo, at, hi) = (f(&d["lower_proven"]), f(&d["undecided_at"]), f(&d["upper_witness"]));
    assert!(100.0 < lo && lo <= at && at < hi && hi < 110.0, "{}", pretty(d));
}

#[test]
fn shapes_far_outside_the_window_have_bounds_volume_and_claims() {
    let s = sc("far");
    let far = s.model("beyond.ink", FAR);
    let d = &run(&s, &["bounds", &far, "--target", "a"])["targets"][0];
    let (mn, mx, cell) = (nums(&d["min"]), nums(&d["max"]), f(&d["cell"]));
    assert!(d.get("unbounded").is_none() && 1001.0 <= mx[0] && mx[0] <= 1001.0 + 2.0 * cell && mn[0] <= -1.0, "{}", pretty(d));
    let v = &run(&s, &["mass", &far, "--target", "root"])["targets"][0]["volume"];
    let ball = 4.0 * std::f64::consts::PI / 3.0;
    assert!((f(&v["estimate"]) / ball - 1.0).abs() < 0.02, "an instance placed far out: {}", pretty(v));
    assert_eq!(short(&s, &[&far]), map([("gone", "fails"), ("near", "fails")]));
    // A sphere at x = 1000 has volume and overlaps itself; at a coarse
    // --res too, not unknown.
    let clipped = s.model(
        "clipped.ink",
        "(import \"std\")\n(body far (shape (move (sphere 2) 1000 0 0)))\n\
        (assert volume-positive (volume (shape-of far) > 32))\n(assert self-overlap (overlaps (shape-of far) (shape-of far)))\n",
    );
    let both = map([("self-overlap", "holds"), ("volume-positive", "holds estimate")]);
    assert_eq!(short(&s, &[&clipped]), both);
    assert_eq!(short(&s, &[&clipped, "--res", "2"]), both, "at --res 2");
    let d = &run(&s, &["bounds", &clipped])["targets"][0];
    assert!((f(&d["min"][0]) - 998.0).abs() < 0.2 && (f(&d["max"][0]) - 1002.0).abs() < 0.2, "were empty: {}", pretty(d));
    // A witness and a volume inside the window stand.
    let wit = s.model("witness.ink", "(import \"std\")\n(def f (min (sphere 1) (move (sphere 1) 1000 0 0)))\n(assert self (overlaps f f))\n(assert lower (volume f > 0.01))\n");
    assert_eq!(short(&s, &[&wit, "--res", "0.5"]), map([("lower", "holds estimate"), ("self", "holds")]));
}

#[test]
fn a_part_moving_out_of_the_window_still_overlaps_itself() {
    // It failed at t = 0.75.
    let s = sc("moving");
    let m = s.model(
        "time.ink",
        "(import \"std\")\n(body root (shape (move (sphere 1) (* 1000 t) 0 0)))\n\
        (assert missing (empty (shape-of root)))\n(assert self (overlaps (shape-of root) (shape-of root)))\n",
    );
    let r = short(&s, &[&m, "--res", "0.5", "--t-span", "0:1"]);
    assert!(r["missing"] == "fails" && r["self"] != "fails", "{r:?}");
}

#[test]
fn interferes_boxes_hold_pieces_far_outside_the_window() {
    // Both boxes include the shared sphere at x = 1000: the pair overlaps,
    // and is not dismissed as apart by its bounding boxes.
    let s = sc("broad");
    let m = s.model(
        "broad.ink",
        "(import \"std\")\n(body a (shape (min (sphere 1) (move (sphere 1) 1000 0 0))))\n\
        (body b (shape (min (move (sphere 1) 10 0 0) (move (sphere 1) 1000 0 0))))\n(assembly root (instance a a) (instance b b))\n",
    );
    let d = run(&s, &["interfere", &m, "--res", "2"]);
    assert!(d["overlapping_pairs"] == 1 && d["pairs_apart_by_bounding_boxes"] == 0, "{}", pretty(&d));
}

#[test]
fn a_body_modelled_outside_the_window_and_placed_in_it_is_analysed_without_warning() {
    let s = sc("ownframe");
    let m = s.model("ownframe.ink", "(import \"std\")\n(body outside (shape (move (sphere 2) 1000 0 0)))\n(assembly root (instance here outside (at (move -1000 0 0))))\n");
    let d = &run(&s, &["bounds", &m])["targets"][0];
    assert!((f(&d["min"][0]) + 2.0).abs() < 0.2 && (f(&d["max"][0]) - 2.0).abs() < 0.2, "was empty: {}", pretty(d));
    let d = &run(&s, &["mass", &m])["targets"][0];
    let ball = 32.0 * std::f64::consts::PI / 3.0;
    assert!((f(&d["volume"]["estimate"]) / ball - 1.0).abs() < 0.02 && d.get("unbounded").is_none(), "{}", pretty(d));
    let w = &run(&s, &["check", &m])["warnings"];
    assert!(w.as_array().is_some_and(Vec::is_empty), "{w}");
}

#[test]
fn instances_turned_by_t_are_analysed_once_in_the_bodys_frame() {
    // Each query analyses the body once in its own frame and moves the
    // answer, which still holds the true value.
    let s = sc("turned");
    let inst: String = (0..3).map(|k| format!(" (instance b{k} block (at (rot-z (* {} t)) (move {} 0 0)))", 37 * k + 11, 10 * k)).collect();
    let m = s.model("turned.ink", &format!("(import \"std\")\n(body block (shape (box 2 4 6)) (attr material \"steel\"))\n(assembly root{inst})\n"));
    let c = run(&s, &["mass", &m, "--t", "0.3"]);
    let v: Vec<&Value> = c["targets"].as_array().unwrap().iter().map(|x| &x["volume"]).collect();
    assert!(v.iter().all(|x| (f(&x["estimate"]) / 48.0 - 1.0).abs() < 0.02 && *x == v[0]), "one volume, found in the body's frame, 48: {v:?}");
    for (k, x) in c["targets"].as_array().unwrap().iter().enumerate() {
        let p = nums(&x["center_of_mass"]);
        let d = ((p[0] - 10.0 * k as f64).powi(2) + p[1].powi(2) + p[2].powi(2)).sqrt();
        assert!(d < 1e-3, "centre of mass {k} where placed: {p:?}");
    }
    let b = run(&s, &["bounds", &m, "--t", "0.3"]);
    let b = b["targets"].as_array().unwrap();
    for k in 0..3 {
        let a = ((37 * k + 11) as f64 * 0.3).to_radians();
        let half = [a.cos().abs() + 2.0 * a.sin().abs(), a.sin().abs() + 2.0 * a.cos().abs(), 3.0];
        let c = [10.0 * k as f64, 0.0, 0.0];
        let (mn, mx) = (nums(&b[k]["min"]), nums(&b[k]["max"]));
        assert!((0..3).all(|i| mn[i] <= c[i] - half[i] && c[i] + half[i] <= mx[i]), "box {k} holds it: {mn:?} {mx:?}");
    }
    let d = run(&s, &["mesh", &m, "--t", "0.3", "--res", "0.5", "--out", &s.path("turned.stl")]);
    assert!(d["instances"] == 3 && d["distinct_meshes"] == 1, "one mesh: {} {}", d["instances"], d["distinct_meshes"]);
}

/// A placement whose numbers are not a turn (an axis too small to normalise
/// in doubles: R^T R is not I), at an angle ANGLE.
fn nonrigid(angle: &str) -> String {
    format!("(import \"std\")\n(body b (shape (box 2 4 6)))\n(assembly root (instance b b (at (rot 1.6e-162 0 0 {angle}))))\n(def world (shape-of root/b))\n")
}

fn placed_by_no_turn_is_its_world_field(angle: &str) {
    // The instance is analysed in the world: its mass and mesh are its world
    // field's (moved as if turned, they would be a body of a third the
    // volume, and a mesh 1.3 inside its surface).
    let s = sc(&format!("nonrigid-{}", angle.len()));
    let m = s.model("nonrigid.ink", &nonrigid(angle));
    let m = m.as_str();
    let got: Vec<Vec<f64>> = ["root/b", "world"]
        .iter()
        .map(|target| {
            let at = ["--t", "1", "--target", *target];
            let ms = run(&s, &[&["mass", m][..], &at, &["--density", "7.85g/cm3"]].concat());
            let sh = run(&s, &[&["mesh", m][..], &at, &["--res", "0.5", "--out", &s.path("nonrigid.stl")]].concat());
            let (ms, sh) = (&ms["targets"][0], &sh["objects"][0]);
            let mut v = vec![f(&ms["mass"])];
            v.extend(nums(&ms["center_of_mass"]));
            v.extend([f(&sh["mesh_volume"]), f(&sh["deviation"]["max"])]);
            v
        })
        .collect();
    let close = |x: f64, y: f64| (x - y).abs() <= (1e-6 * x.abs().max(y.abs())).max(1e-9);
    assert!(got[0].len() == 6 && got[0].iter().zip(&got[1]).all(|(x, y)| close(*x, *y)), "root/b {:?}\nworld  {:?}", got[0], got[1]);
}

#[test]
fn an_instance_placed_by_no_turn_is_analysed_as_its_world_field() {
    placed_by_no_turn_is_its_world_field("90");
}

#[test]
fn an_instance_placed_by_no_turn_at_a_moment_is_analysed_as_its_world_field() {
    placed_by_no_turn_is_its_world_field("(* 90 t)");
}

#[test]
fn number_claims_are_not_shapes() {
    // Each holds, -1 < 0 included: a number is no shape.
    let s = sc("scalar");
    let m = s.model("scalar.ink", "(assert negative (value -1 < 0))\n(assert zero (value 0 < 1))\n(assert positive (value 1 < 2))\n");
    let r = short(&s, &[&m]);
    assert!(r.len() == 3 && r.values().all(|v| v == "holds"), "{r:?}");
}

#[test]
fn claims_read_before_the_solve_of_their_param_are_estimates_at_the_solved_value() {
    // Were guaranteed, at the unsolved value.
    let s = sc("num-before");
    let nb = s.model(
        "num-before.ink",
        "(import \"std\")\n(param u 1)\n(def old u)\n(body ball (shape (sphere u)))\n(assert before (value u < 2))\n\
        (solve fit (unknowns u) (== (square (- u 2)) 0))\n(assert via-old (value old < 2))\n\
        (assert via-body (within (shape-of ball) (-3 -3 -3) (3 3 3)))\n(assert direct (value u < 2))\n",
    );
    let r = short(&s, &[&nb, "--res", "2"]);
    assert!(r.len() == 4 && r.values().all(|v| v == "holds estimate"), "{r:?}");
    let nbb = s.model("num-bound.ink", "(param u 1)\n(assert below (value 1.5 < u))\n(solve fit (unknowns u) (== u 3))\n");
    let d = &run(&s, &["assert", &nbb])["asserts"][0];
    assert_eq!(verdict(d), "holds estimate", "a bound read before the solve: {}", pretty(d));
    // A solved placement far outside the window fails, an estimate (held).
    let nc = s.model(
        "num-clipped.ink",
        "(import \"std\")\n(solve fit (unknowns (u 1)) (== u 1000))\n(body far (shape (move (sphere 2) u 0 0)))\n(assert missing (empty (shape-of far)))\n",
    );
    assert_eq!(short(&s, &[&nc]), map([("missing", "fails estimate")]));
    // Optimization requirements read their bounds at the solved values too.
    let req = s.model(
        "req.ink",
        "(param u 1)\n(param p 0.5)\n(optimize choose (vary p 0 1) (minimize (value p)) (require (value 1.5 < u)) (budget 4))\n\
        (solve fit (unknowns u) (== u 3))\n(assert same (value 1.5 < u))\n",
    );
    let d = run(&s, &["optimize", &req]);
    assert!(d["feasible"] == true, "{}", pretty(&d));
    assert_eq!(short(&s, &[&req]), map([("same", "holds estimate")]));
}

/// A four-bar that closes only near t = 0.5: no root at load time (t = 0).
const FOURBAR: &str = "(param a 15) (param b 60) (param c 15) (param d 55)\n(def th (* 2 pi t))\n\
    (solve loop (unknowns (th3s 0.3) (th4s 2.5))\n\
    \x20 (== (+ (* a (cos th)) (* b (cos th3s))) (+ d (* c (cos th4s))))\n\
    \x20 (== (+ (* a (sin th)) (* b (sin th3s))) (* c (sin th4s))))\n\
    (assert rocker-up (value (sin th4s) > 0))\n(assert crank-x (value (* a (cos th)) < 15.1))\n";

#[test]
fn a_solve_open_at_load_time_leaves_only_the_claims_reading_it_unknown() {
    // It has no root at load time to follow along t, so at 0.5, at t = 0
    // and over a span only the claims reading it are unknown, each saying
    // why; the rest are decided (the run stopped at the first; it held, on
    // whichever root its guesses found).
    let s = sc("fourbar");
    let fb = s.model("fourbar.ink", FOURBAR);
    let why = |args: &[&str]| -> BTreeMap<String, (String, String)> {
        run(&s, &[&["assert", fb.as_str()][..], args].concat())["asserts"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| (a["name"].as_str().unwrap().to_string(), (verdict(a), a["unavailable"].as_str().unwrap_or("").to_string())))
            .collect()
    };
    let r = why(&["--t", "0.5"]);
    assert!(r["rocker-up"].0 == "unknown" && r["rocker-up"].1.contains("no root at t = 0 to start from") && r["crank-x"] == ("holds".to_string(), String::new()), "{r:?}");
    let r = why(&[]);
    assert!(r["rocker-up"].0 == "unknown" && r["rocker-up"].1.contains("did not converge at t = 0") && r["crank-x"] == ("holds".to_string(), String::new()), "{r:?}");
    // Over a motion it is solved at each step: unknown where it does not
    // close (with those t), the crank's value holds.
    let d = run(&s, &["assert", &fb, "--t-span", "0:1"]);
    let und = d["asserts"][0]["motion"]["undecided"].as_array().cloned().unwrap_or_default();
    assert!(statuses(&d) == map([("crank-x", "holds"), ("rocker-up", "unknown")]) && !und.is_empty() && und[0][0] == 0.0, "{}", pretty(&d["asserts"][0]));
    let w = &run(&s, &["check", &fb])["warnings"];
    let w: Vec<&str> = w.as_array().unwrap().iter().filter_map(|w| w.as_str()).collect();
    assert!(w.len() == 2 && w.iter().any(|w| w.contains("did not converge at t = 0")) && w.iter().any(|w| w.contains("states no t-span")), "{w:?}");
    let e = refused(&s, &["eval", &fb, "--points", "[[0,0,0]]"]);
    assert!(msg(&e).contains("did not converge"), "another command needs it solved: {}", pretty(&e));
}

#[test]
fn interfere_leaves_out_a_part_resting_on_an_open_solve_and_decides_the_rest() {
    // A pair whose zone rests on it is unknown (it stopped at the first).
    let s = sc("unsolved-parts");
    let m = s.model(
        "unsolved-parts.ink",
        "(import \"std\")\n(solve s (unknowns (q 1)) (== (square q) (- t 1)))\n(body a (shape (sphere 1)))\n\
        (assembly root (instance x a (at (move 0 0 0))) (instance y a (at (move 1.5 0 0))) (instance w a (at (move 1.5 1.5 0)))\n\
        \x20 (instance z a (at (move (* 10 q) 0 0))) (contact x y (sphere (+ 3 q))))\n",
    );
    let d = run(&s, &["interfere", &m]);
    let items: Vec<String> = d["unsolved_items"].as_array().map_or(vec![], |u| u.iter().map(|u| u["item"].as_str().unwrap_or("?").to_string()).collect());
    assert!(items == ["root/z"] && full_pairs(&d) == map([("root/x~root/y", "unknown"), ("root/x~root/w", "holds"), ("root/y~root/w", "fails")]), "{}", pretty(&d));
}

#[test]
fn an_unbounded_shape_is_examined_in_the_window_and_what_needs_all_of_it_is_unknown() {
    // A half-space: witnesses stand, emptiness and volumes are unknown.
    let s = sc("half");
    let half = s.model(
        "half.ink",
        "(import \"std\")\n(body h (shape x))\n(assert self (overlaps x x))\n(assert lower (volume x > 1))\n\
        (assert upper (volume x < 1e12))\n(assert gone (empty x))\n(assert cut (empty (max x (move (sphere 1) 2000 0 0))))\n",
    );
    assert_eq!(short(&s, &[&half, "--res", "10"]), map([("cut", "holds"), ("gone", "fails"), ("lower", "unknown"), ("self", "holds"), ("upper", "unknown")]));
    let b = run(&s, &["bounds", &half]);
    assert!(b["targets"][0].get("unbounded").is_some() && b["union"].is_null(), "{}", pretty(&b));
    let w = warnings(&run(&s, &["check", &half]), "'h' is not proven bounded");
    assert!(!w.is_empty(), "check warns about a body not proven bounded");
}

#[test]
fn interfere_checks_an_unbounded_body_against_every_other() {
    let s = sc("half-pair");
    let m = s.model(
        "half-pair.ink",
        "(import \"std\")\n(body floor (shape z))\n(body ball (shape (move (sphere 1) 1000 0 5)))\n(body hit (shape (move (sphere 1) 1000 0 0)))\n\
        (assembly root (instance floor floor) (instance ball ball) (instance hit hit))\n",
    );
    let d = run(&s, &["interfere", &m]);
    let st = full_pairs(&d);
    assert!(
        st.get("root/floor~root/ball").is_some_and(|v| v == "holds")
            && st.get("root/floor~root/hit").is_some_and(|v| v == "fails")
            && strs(&d["unbounded_items"]) == ["root/floor"],
        "{}",
        pretty(&d)
    );
    // Two half-spaces overlap without bound: no volume.
    let halves = s.model("halves.ink", "(body half (shape x))\n(assembly root (instance a half) (instance b half))\n");
    let p = &run(&s, &["interfere", &halves, "--res", "1"])["pairs"][0];
    assert!(verdict(p) == "fails" && p["overlap_volume"]["estimate"].is_null(), "{}", pretty(p));
    // Two unbounded items whose intersection is bounded, or empty: settled.
    // (Two models: a and b against rod and cap would be pairs whose
    // unbounded overlap volumes are measured over the whole window.)
    let apart = s.model("apart.ink", "(body a (shape (+ x 1)))\n(body b (shape (- 1 x)))\n(assembly root (instance a a) (instance b b))\n");
    let capped = s.model(
        "capped.ink",
        "(body rod (shape (max (- 1000 x) (- (abs y) 1) (- (abs z) 1))))\n\
        (body cap (shape (- x 1002)))\n(assembly root (instance rod rod) (instance cap cap))\n",
    );
    let mut st = full_pairs(&run(&s, &["interfere", &apart, "--res", "0.5"]));
    st.extend(full_pairs(&run(&s, &["interfere", &capped, "--res", "0.5"])));
    assert!(st.get("root/a~root/b").is_some_and(|v| v == "holds") && st.get("root/rod~root/cap").is_some_and(|v| v == "fails"), "{st:?}");
}

#[test]
fn interfere_asks_to_a_tenth_of_a_mm_and_samples_overlap_volume_every_accuracy() {
    // 0.1 mm, not 1/20,000 of the window (0.0002 in a 4 mm one, which this
    // 0.025-deep overlap would fail).
    let s = sc("interfere-accuracy");
    let d = run(
        &s,
        &[
            "interfere",
            &s.model(
                "interfere-accuracy.ink",
                "(import \"std\")\n(body b (shape (box 1 1 1)))\n\
        (assembly root (instance a b) (instance c b (at (move 0.95 0 0))))\n",
            ),
        ],
    );
    assert!(d["accuracy"] == 0.1 && verdict(&d["pairs"][0]) == "holds", "{}", pretty(&d));
    // A failing pair's overlap volume is sampled every accuracy, not at the
    // proof's cells: it stops where it finds a witness, here a 12.8 mm cell,
    // whose samples read these plates' 0.3 mm overlap (750 mm3) as 5812.
    let p = run(
        &s,
        &[
            "interfere",
            &s.model(
                "overlap-thin.ink",
                "(import \"std\")\n(body plate (shape (box 50 50 5)))\n\
        (assembly root (instance a plate (at (rot-x 30))) (instance b plate (at (move 0 0 4.7) (rot-x 30))))\n",
            ),
        ],
    )["pairs"][0]
        .clone();
    let v = &p["overlap_volume"];
    assert!(verdict(&p) == "fails" && (f(&v["estimate"]) - 750.0).abs() <= 15.0, "{}", pretty(&p));
    // Asked in a shared frame (two bodies turned alike, one moved 4.7 mm
    // along z): 40 x 40 x (5 - 4.7 cos 30 deg) mm.
    let p = run(
        &s,
        &[
            "interfere",
            &s.model(
                "overlap-framed.ink",
                "(import \"std\")\n(body p (shape (box 50 50 5)))\n(body q (shape (box 40 40 5)))\n\
        (assembly root (instance a p (at (rot-x 30))) (instance b q (at (rot-x 30) (move 0 0 4.7))))\n",
            ),
        ],
    )["pairs"][0]
        .clone();
    let v = &p["overlap_volume"];
    let truth = 1600.0 * (5.0 - 4.7 * 30f64.to_radians().cos());
    assert!(verdict(&p) == "fails" && (f(&v["estimate"]) - truth).abs() <= 0.02 * truth, "{}", pretty(&p));
}

#[test]
fn a_field_never_negative_is_empty_even_where_it_is_zero() {
    let s = sc("zero");
    let m = s.model("zero.ink", "(assert none (empty 0))\n(assert square (empty (square x)))\n(assert some (empty (- 0 (square x))))\n");
    assert_eq!(short(&s, &[&m, "--res", "10"]), map([("none", "holds"), ("some", "fails"), ("square", "holds")]));
}

#[test]
fn seated_faces_hold_to_the_accuracy_at_cells_its_size_sets() {
    // A small cube seated on a part is proven to the accuracy (it overlaps
    // nowhere deeper), at cells its own size sets; a finer --res is a size
    // the search may reach, not cells it must use.
    let s = sc("seat");
    let seat = s.model("seat.ink", "(import \"std\")\n(def f (max (box 2 2 2) (move (box 0.03 0.03 0.03) 1.015 0 0)))\n(assert seated (empty f) (accuracy 0.01))\n");
    let a = run(&s, &["assert", &seat])["asserts"][0].clone();
    let b = run(&s, &["assert", &seat, "--res", "2e-4"])["asserts"][0].clone();
    assert!(a["accuracy"] == 0.01 && b["res_given"] == 2e-4 && verdict(&a) == "holds" && verdict(&b) == "holds", "{}\n{}", pretty(&a), pretty(&b));
    // Faces seated by design are proven to the model's accuracy (0.05 mm in
    // the default window), at cells no finer than half of it.
    let flush = s.model("flush.ink", "(import \"std\")\n(def f (max (box 20 20 20) (move (box 10 10 10) 15 0 0)))\n(assert seated (empty f))\n");
    let a = &run(&s, &["assert", &flush])["asserts"][0];
    assert!(verdict(a) == "holds" && a["accuracy"] == 0.05 && f(&a["res"]) >= 0.025, "{}", pretty(a));
}

// ------------------------------------------------------------ gear meshes

/// Two gears 0.08 mm apart (each tooth thinned 0.0851 at the pitch circle,
/// the pair centred: 0.0851 cos 20 between the flanks); the assembly left
/// open.
const APART: &str = "(import \"std\")\n\
    (fn gear (n) (subtract (extrude-z (involute-gear-2d n 1 20 0.0851) 0 8) (cylinder-z 2 20)))\n\
    (body pinion (shape (gear 15)))\n(body wheel (shape (gear 45)))\n\
    (assembly root (instance p pinion)\n  (instance w wheel (at (rot-z (mesh-turn 15 45 0 0)) (move (mesh-distance 15 45 1) 0 0)))";

#[test]
fn a_mesh_holds_its_gears_whatever_their_backlash() {
    // In mesh, the mesh's claims hold.
    let s = sc("mesh-apart");
    let held = s.model("mesh-held.ink", &format!("{APART}\n  (relate p w (gear-mesh 15 45 1 20 0.1702)))\n"));
    let r = short(&s, &[&held]);
    assert!(!r.is_empty() && r.values().all(|v| v == "holds"), "{r:?}");
}

/// A wheel whose axis is its x (converted CAD), turned onto the pinion's.
const TURNED: &str = "(import \"std\")\n(fn gear (n) (subtract (extrude-z (involute-gear-2d n 1 20 0.1) 0 8) (cylinder-z 2 20)))\n\
    (body pinion (shape (gear 15)))\n(body wheel (shape (place (gear 45) (at (rot-y 90)))))\n\
    (assembly root (instance p pinion)\n\
    \x20 (instance w wheel (at (rot-y -90) (rot-z (mesh-turn 15 45 0 0)) (move (mesh-distance 15 45 1) 0 0)))\n";

fn all_hold() -> BTreeMap<String, String> {
    map([("centre-distance", "holds"), ("faces-meet", "holds"), ("parallel-axes", "holds")])
}

#[test]
fn a_gear_whose_axis_is_not_its_z_is_related_at_a_frame_on_its_axis() {
    let s = sc("mesh-framed");
    let r = short(&s, &[&s.model("mesh-framed.ink", &format!("{TURNED}  (relate p (w (at (rot-y 90))) (gear-mesh 15 45 1 20 0.2)))\n"))]);
    assert_eq!(r, all_hold(), "at (PATH FRAME)");
    let r = short(&s, &[&s.model("mesh-unframed.ink", &format!("{TURNED}  (relate p w (gear-mesh 15 45 1 20 0.2)))\n"))]);
    assert_eq!(r["parallel-axes"], "fails", "at its own frame its z is not its axis: {r:?}");
    let e = refused(&s, &["check", &s.model("mesh-end.ink", &format!("{TURNED}  (relate p (w (at) 2) (gear-mesh 15 45 1)))\n"))]);
    assert!(msg(&e).contains("(PATH FRAME)"), "an end neither a path nor (PATH FRAME): {}", pretty(&e));
}

#[test]
fn an_ends_frame_is_in_its_parts_own_coordinates() {
    // The kind's zone at such an end is in the frame's coordinates, put
    // back in the part's: a cube's corner 5 mm along its x, overlapping a
    // block, zoned by a sphere about the frame's origin.
    let s = sc("end-frame");
    let d = run(
        &s,
        &[
            "interfere",
            &s.model(
                "zone-framed.ink",
                "(import \"std\")\n(body cube (shape (box 10 10 10)))\n(body block (shape (box 2 2 2)))\n\
        (assembly root (instance a cube) (instance b block (at (move 5 0 0)))\n\
        \x20 (relate (a (at (move 5 0 0))) b (lambda (fa fb sa sb) (relation corner (zone \"the corner\" (sphere 2))))))\n",
            ),
        ],
    );
    assert!(d["overlapping_pairs"] == 0 && verdict(&d["pairs"][0]) == "holds", "{}", pretty(&d["pairs"]));
    // A frame the part declares, by name, is read where the part is;
    // (frame-of PART ...) is in the assembly's coordinates, and would be
    // applied twice (a belt 100 mm off), so it is refused, naming the way.
    let named = "(import \"std\")\n(body cube (shape (box 10 10 10)) (frame hub (at (move 5 0 0))))\n(body block (shape (box 2 2 2)))\n\
        (assembly root (instance a cube (at (move 100 0 0))) (instance b block (at (move 200 0 0)))\n\
        \x20 (relate (a END) b (lambda (fa fb sa sb) (relation at (claim x (nth (get fa \"origin\") 0) 105 105)))))\n";
    assert_eq!(short(&s, &[&s.model("end-named.ink", &named.replace("END", "hub"))]), map([("x", "holds")]));
    let e = refused(&s, &["check", &s.model("end-frame-of.ink", &named.replace("END", "(frame-of a hub)"))]);
    assert!(msg(&e).contains("(PART NAME)"), "{}", pretty(&e));
}

#[test]
fn a_bevel_pairs_claims_fail_pushed_aside_tilted_or_slid_and_its_teeth_are_checked() {
    // Straight bevel gears at right angles, apexes at the origin: the
    // pinion pushed toward its apex by d moves its apex d and closes their
    // pitch cones by d sin of its pitch angle (0.1 mm: both within the
    // 0.137 a 0.1 backlash allows; 1 mm: 0.316, past it); the wheel moved
    // sideways, its axis misses the pinion's apex; tilted a degree, they are
    // not square, and its pitch cone swings 0.5 mm through the pitch point;
    // the pinion slid 30 mm along its pitch line keeps the cones touching
    // and the axes square, its apex 30 mm off.
    let s = sc("bevel");
    let bevel = "(import \"std\")\n(param push 0)\n(param side 0)\n(param tilt 0)\n(param slide 0)\n\
        (fn bevel (n d) (max (straight-bevel-gear n 1.5 20 0 (rad->deg d) 25.576 35.576) (- (len3 x y z) 40)))\n\
        (body pinion (shape (bevel 15 (atan2 15 45))))\n(body wheel (shape (bevel 45 (atan2 45 15))))\n\
        (assembly root (instance p pinion (at (move (* slide (sin (atan2 15 45))) 0 (+ push (* slide (cos (atan2 15 45)))))))\n\
        \x20 (instance w wheel (at (rot-y (+ 90 tilt)) (move 0 side 0)))\n\
        \x20 (relate p w (bevel-mesh 15 45 1.5 10 20 0.1)))\n";
    let b = s.model("bevel.ink", bevel);
    let r: Vec<Vec<String>> = [&[][..], &["--param", "push=0.1"], &["--param", "push=1"], &["--param", "side=0.5"], &["--param", "tilt=1"], &["--param", "slide=30"]]
        .iter()
        .map(|p| failed(&s, &[&[b.as_str()][..], *p].concat()))
        .collect();
    let want: [&[&str]; 6] = [&[], &[], &["apexes", "pitch-cones"], &["apexes"], &["pitch-cones", "square-axes"], &["apexes"]];
    assert!(r.iter().zip(want).all(|(x, y)| x == y), "{r:?}");
    let e = refused(&s, &["check", &s.model("bevel-centre.ink", &bevel.replace("(bevel-mesh 15 45", "(bevel-mesh (vec 15 16) 45"))]);
    assert!(msg(&e).contains("plain tooth counts"), "(vec teeth centre-teeth): {}", pretty(&e));
    let d = run(&s, &["interfere", &b]);
    assert!(verdict(&d["pairs"][0]) == "holds" && d["pairs"][0]["relations"][0].get("skipped").is_none(), "teeth checked, no zone: {}", pretty(&d["pairs"]));
    let a = &run(&s, &["check", &b])["assemblies"][0]["relations"][0];
    assert!(a["kind"] == "bevel-mesh" && f(&a["data"]["ratio"]) == 3.0 && (f(&a["data"]["pitch-angles"][0]) - 18.43494882).abs() < 1e-6, "{}", pretty(a));
}

#[test]
fn a_pinion_on_std_s_rack_meshes_and_fails_lifted_past_the_end_or_beside() {
    // Its axis its pitch radius above the rack's pitch line; lifted 0.3 mm,
    // past the 0.137 a 0.1 backlash allows.
    let s = sc("rack");
    let rk = s.model(
        "rack.ink",
        "(import \"std\")\n(param lift 0)\n(param slide 0)\n(param across 0)\n\
        (body pinion (shape (subtract (extrude-z (involute-gear-2d 20 1 20 0.05) 0 8) (cylinder-z 2 20))))\n\
        (body rack (shape (max (extrude-z (involute-rack-2d 1 20 0.05) 0 8) (slab-x -40 40) (- -4 y))))\n\
        (assembly root (instance p pinion (at (rot-z 9) (move slide (+ 10 lift) across))) (instance r rack)\n\
        \x20 (relate p r (rack-mesh 20 1 20 0.1)))\n",
    );
    let r: Vec<Vec<String>> =
        [&[][..], &["--param", "lift=0.3"], &["--param", "slide=200"], &["--param", "across=10"]].iter().map(|p| failed(&s, &[&[rk.as_str()][..], *p].concat())).collect();
    let want: [&[&str]; 4] = [&[], &["pitch-line"], &["faces-meet"], &["faces-meet"]];
    assert!(r.iter().zip(want).all(|(x, y)| x == y), "{r:?}");
    // They mesh without clashing; the rack's lead is pi m z a turn.
    let d = run(&s, &["interfere", &rk]);
    let a = &run(&s, &["check", &rk])["assemblies"][0]["relations"][0];
    assert!(verdict(&d["pairs"][0]) == "holds" && (f(&a["data"]["lead"]) - std::f64::consts::PI * 20.0).abs() < 1e-6, "{}\n{}", pretty(&d["pairs"]), pretty(a));
}

#[test]
fn a_mesh_whose_pinion_is_built_with_its_assemblys_argument_loads_and_holds() {
    // An end built with an argument of an assembly built again (a pinion an
    // option swaps) is read as the body it is: its name, pair/p, is no path
    // (not "'root' is not an assembly").
    let s = sc("mesh-swapped");
    let m = s.model(
        "mesh-swapped.ink",
        "(import \"std\")\n(fn gear (n) (subtract (extrude-z (involute-gear-2d n 1 20 0.1) 0 8) (cylinder-z 2 20)))\n\
        (body (pinion (n 15)) (shape (gear n)))\n(body wheel (shape (gear 45)))\n\
        (assembly (pair (k 15)) (instance p (pinion k))\n\
        \x20 (instance w wheel (at (rot-z (mesh-turn k 45 0 0)) (move (mesh-distance k 45 1) 0 0)))\n\
        \x20 (relate p w (gear-mesh k 45 1 20 0.2)))\n(assembly root (instance g (pair 16)))\n",
    );
    assert_eq!(short(&s, &[&m]), all_hold());
}

#[test]
fn a_backlash_written_where_the_pressure_angle_goes_is_an_error_naming_the_order() {
    // Pressure angle first, involute-gear-2d's order: backlash first reads
    // it as a pressure angle no gear has, an error (as a warning it would
    // let a wheel 50 mm off and tilted 30 degrees pass).
    let s = sc("mesh-old");
    let e = refused(&s, &["check", &s.model("mesh-old.ink", &format!("{APART}\n  (relate p w (gear-mesh 15 45 1 0.1702 20)))\n"))]);
    assert!(msg(&e).contains("pa [backlash [efficiency]]"), "{}", pretty(&e));
}

/// Spur gears in mesh, the wheel lifted along its axis by `lift`.
const LIFTED: &str = "(import \"std\")\n(param lift 0)\n(fn gear (n) (subtract (extrude-z (involute-gear-2d n 1 20 0.1) 0 8) (cylinder-z 2 20)))\n\
    (body pinion (shape (gear 15)))\n(body wheel (shape (gear 45)))\n\
    (assembly root (instance p pinion)\n  (instance w wheel (at (rot-z (mesh-turn 15 45 0 0)) (move (mesh-distance 15 45 1) 0 lift)))\n\
    \x20 (relate p w (gear-mesh 15 45 1)))\n";

#[test]
fn faces_meet_is_read_off_both_shapes_where_the_mates_teeth_work() {
    let s = sc("mesh-faces");
    let m = s.model("mesh-lift.ink", LIFTED);
    let d = run(&s, &["assert", &m]);
    assert_eq!(statuses(&d).into_iter().map(|(k, v)| (k.rsplit('/').next().unwrap().to_string(), v)).collect::<BTreeMap<_, _>>(), all_hold());
    // Its margin is the faces' overlap, 8 mm read within the boxes' slack.
    let fm = d["asserts"].as_array().unwrap().iter().find(|x| x["name"].as_str().unwrap().ends_with("faces-meet")).unwrap();
    let (v, hi) = (f(&fm["terms"]["faces-meet"]), f(&fm["terms"]["hi"]));
    assert!(7.9 < v && v < 8.4 && hi > 2.0 * v - 1e-6, "{}", pretty(fm));
    // Spur gears 100 mm apart along the axis fail faces-meet.
    assert_eq!(short(&s, &[&m, "--param", "lift=100"]), map([("centre-distance", "holds"), ("faces-meet", "fails"), ("parallel-axes", "holds")]));
    // A body with no material in the ring of its teeth (a disc for a
    // 15-tooth pinion) is refused, saying so.
    let disc = s.model("mesh-disc.ink", &LIFTED.replace("(body pinion (shape (gear 15)))", "(body pinion (shape (extrude-z (circle-2d 5) 0 8)))"));
    let e = refused(&s, &["check", &disc]);
    assert!(msg(&e).contains("where its mate's teeth work"), "{}", pretty(&e));
    // Not down to the root: a pinion cut on a 12.4 mm shaft (its root
    // circle 12.5), the wheel slid onto the shaft, its tips 0.3 mm clear of
    // it, does not mesh.
    let shafted = s.model("mesh-shaft.ink", &LIFTED.replace("(body pinion (shape (gear 15)))", "(body pinion (shape (min (gear 15) (move (cylinder-z 6.2 30) 0 0 23))))"));
    assert_eq!(short(&s, &[&shafted, "--param", "lift=20"])["faces-meet"], "fails");
}

#[test]
fn a_mesh_of_a_param_module_loads_and_one_reading_a_toleranced_face_is_unknown_over_its_range() {
    // Its numbers may read params: its claims' accuracies are numbers, 1e-6
    // mm, not a module's fraction.
    let s = sc("mesh-params");
    let pmod = LIFTED
        .replace("(import \"std\")", "(import \"std\")\n(param module 1)")
        .replace("(involute-gear-2d n 1 20 0.1)", "(involute-gear-2d n module 20 0.1)")
        .replace("(mesh-distance 15 45 1)", "(mesh-distance 15 45 module)")
        .replace("(gear-mesh 15 45 1)", "(gear-mesh 15 45 module)");
    let pm = s.model("mesh-pmod.ink", &pmod);
    assert_eq!(short(&s, &[&pm]), all_hold());
    assert_eq!(short(&s, &[&pm, "--param", "module=1.5"]), all_hold(), "at module 1.5");
    // A ranged param the kind reads at load through a part's shape (the
    // wheel's face) refuses that mesh's claims over its range, saying so,
    // not the file: its other claims are checked.
    let tol = LIFTED
        .replace("(param lift 0)", "(param lift 0)\n(param w 8 (tol -0.1 0.1))")
        .replace("(body wheel (shape (gear 45)))", "(fn wide (n h) (subtract (extrude-z (involute-gear-2d n 1 20 0.1) 0 h) (cylinder-z 2 20)))\n(body wheel (shape (wide 45 w)))")
        + "(assert wide (value w > 7))\n";
    let t = s.model("mesh-tol.ink", &tol);
    let d = run(&s, &["assert", &t]);
    let a: BTreeMap<&str, &Value> = d["asserts"].as_array().unwrap().iter().map(|x| (x["name"].as_str().unwrap().rsplit('/').next().unwrap(), x)).collect();
    assert!(
        a["wide"]["status"] == "holds"
            && ["centre-distance", "faces-meet"].iter().all(|n| a[n]["status"] == "unknown" && a[n]["unavailable"].as_str().unwrap_or("").contains("ranged w")),
        "{}",
        pretty(&d["asserts"])
    );
    let r = short(&s, &[&t, "--nominal"]);
    assert!(r.len() == 4 && r.values().all(|v| v == "holds"), "at its values: {r:?}");
}
