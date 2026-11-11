//! Claims and pairs over a motion: `ink assert` and `ink interfere` over a
//! span of t, and the JSON that says how each was decided (`motion`:
//! `how`, steps; `t` and `witness` where one fails). The tests/motion models' verdicts are
//! geometry.rs `motion_models_hold_or_fail_as_declared`'s; these are the
//! answers' shape and the kernel answers only the CLI asks for.

use crate::common::*;
use serde_json::{Value, json};

/// A call's rows: interfere's pairs or assert's claims.
fn rows(d: &Value) -> &Vec<Value> {
    d.get("pairs").or(d.get("asserts")).and_then(Value::as_array).unwrap_or_else(|| panic!("no rows:\n{}", pretty(d)))
}

/// The first row of an ink call.
fn row(i: Ink) -> Value {
    let d = i.ok();
    rows(&d).first().cloned().unwrap_or_else(|| panic!("no rows:\n{}", pretty(&d)))
}

/// A number, NaN where absent.
fn num(v: &Value) -> f64 {
    v.as_f64().unwrap_or(f64::NAN)
}

/// An array of numbers (a span, a range, a witness).
fn nums(v: &Value) -> Vec<f64> {
    v.as_array().map_or(vec![], |a| a.iter().map(num).collect())
}

fn has(v: &Value, k: &str) -> bool {
    v.get(k).is_some_and(|x| !x.is_null())
}

fn note(v: &Value) -> &str {
    v["note"].as_str().unwrap_or("")
}

/// Each row's verdict.
fn verdicts(d: &Value) -> Vec<String> {
    rows(d).iter().map(verdict).collect()
}

/// The pair `a`~`b` of an interfere answer.
fn pair<'a>(d: &'a Value, a: &str, b: &str) -> &'a Value {
    rows(d).iter().find(|p| p["pair"] == json!([a, b])).unwrap_or_else(|| panic!("no pair {a}~{b}:\n{}", pretty(d)))
}

/// An arm 35 long at z 10 (+ lift) turned a whole turn past a post 30 out.
const ARM: &str = r#"(import "std")
(param lift 0)
(body arm (shape (move (box 35 4 4) 17.5 0 (+ 10 lift))))
(body post (shape (move (cylinder-z 3 20) 30 0 0)))
(assembly root (instance arm arm (at (rot-z (* 360 t)))) (instance post post))
(assert clear (disjoint (shape-of root/arm) (shape-of root/post)))
"#;

#[test]
fn an_arm_turned_into_a_post_fails_at_a_moment_and_lifted_clear_holds_for_all_t() {
    let s = Scratch::new("motion-arm");
    let mo = s.model("motion.ink", ARM);
    let a = row(s.ink(["assert", &mo, "--t-span", "0:1"]));
    assert!(verdict(&a) == "fails" && has(&a, "t"), "{}", pretty(&a));
    let a = row(s.ink(["assert", &mo, "--t-span", "0:1", "--param", "lift=15"]));
    assert!(verdict(&a) == "holds" && nums(&a["motion"]["span"]) == [0.0, 1.0] && a["motion"]["how"] != "measured", "{}", pretty(&a));
    // interfere steps each moving pair too, and times each.
    let d = s.ink(["interfere", &mo, "--t-span", "0:1", "--param", "lift=15"]).ok();
    let n = rows(&d).len();
    assert!(rows(&d).iter().all(|p| has(p, "ms")), "{}", pretty(&d));
    assert!(d["overlapping_pairs"] == 0 && verdicts(&d).iter().all(|v| v == "holds") && n as u64 + d["pairs_apart_by_bounding_boxes"].as_u64().unwrap() == 1, "{}", pretty(&d));
}

#[test]
fn t_span_is_refused_outside_assert_and_interfere_with_a_hint() {
    let s = Scratch::new("motion-tspan-hint");
    let mo = s.model("motion.ink", ARM);
    let d = s.ink(["bounds", &mo, "--t-span", "0:1", "--target", "root/arm"]).err();
    assert!(d["error"]["hint"].as_str().unwrap_or("").contains("--t-span is an option of"), "{}", pretty(&d));
}

#[test]
fn a_motor_riding_a_carriage_on_its_plate_touches_it_asked_once() {
    // Touching at every moment, asked once (they keep their placement
    // relative to each other); one riding 0.2 above it fails so. One lifted
    // off it but at t = 1 moves relative to it, and fails where it floats.
    let s = Scratch::new("motion-seated");
    let seat = s.model(
        "seated.ink",
        r#"(import "std")
(t-span 0 2)
(body rail (shape (box 200 20 4)))
(body plate (shape (move (box 40 30 4) 0 0 4)))
(body motor (shape (move (box 20 20 20) 0 0 16)))
(assembly carriage (instance p plate) (instance m motor) (instance f motor (at (move 0 0 0.2)))
  (instance l motor (at) (slide z (* 0.5 (abs (- t 1))))))
(assembly root (instance r rail) (instance c carriage (at) (slide x (* 40 t))))
(assert seated (touches (shape-of root/c/m) (shape-of root/c/p)))
(assert floating (touches (shape-of root/c/f) (shape-of root/c/p)))
(assert lifted (touches (shape-of root/c/l) (shape-of root/c/p)))
"#,
    );
    let d = s.ink(["assert", &seat]).ok();
    let (seated, floating, lifted) = (&rows(&d)[0], &rows(&d)[1], &rows(&d)[2]);
    assert!(verdict(seated) == "holds" && seated["motion"]["how"] == "asked once", "{}", pretty(seated));
    assert!(verdict(floating) == "fails" && floating["motion"]["how"] == "asked once" && note(floating).contains("not touching"), "{}", pretty(floating));
    assert!(verdict(lifted) == "fails" && num(&lifted["t"]) != 1.0, "{}", pretty(lifted));
}

#[test]
fn pairs_riding_one_carriage_are_still_questions_and_a_belt_reading_t_is_boxed_over_the_motion() {
    // The carriage's parts ride one rigid group, so each pair of them is a
    // still question, asked once in its own frame (no motion block), the
    // frame's rounding bounded over the span: seated holds, sunk 1 fails
    // (at every moment), and two cubes seated alike on two pads are one
    // question. A belt whose field reads t has a box over the motion, so
    // it is paired only with what that box meets, not with every item.
    let s = Scratch::new("motion-riding");
    let group = s.model(
        "riding.ink",
        r#"(import "std")
(t-span 0 2)
(body rail (shape (box 200 20 4)))
(body plate (shape (move (box 40 30 4) 0 0 4)))
(body motor (shape (move (box 20 20 20) 0 0 16)))
(body pad (shape (box 6 6 2)))
(body cube (shape (move (box 4 4 4) 0 0 3)))
(body belt (shape (move (box (+ 100 (* 20 t)) 6 1.4) 0 -100 0)))
(assembly carriage (instance p plate) (instance m motor (at (move 10 0 0))) (instance x motor (at (move -10.5 0 -1)))
  (instance pa pad (at (move -10 40 0))) (instance ca cube (at (move -10 40 0))) (instance pb pad (at (move 10 40 0))) (instance cb cube (at (move 10 40 0))))
(assembly root (instance r rail) (instance b belt) (instance c carriage (at) (slide x (* 40 t))))
"#,
    );
    let d = s.ink(["interfere", &group]).ok();
    let riding = |p: &Value| p["pair"].as_array().unwrap().iter().all(|x| x.as_str().unwrap().starts_with("root/c/"));
    assert!(
        verdict(pair(&d, "root/c/p", "root/c/m")) == "holds"
            && verdict(pair(&d, "root/c/p", "root/c/x")) == "fails"
            && pair(&d, "root/c/pb", "root/c/cb")["cached"] == true
            && rows(&d).iter().all(|p| riding(p) != has(p, "motion")),
        "{}",
        pretty(&d)
    );
    assert!(!has(&d, "unbounded_items") && !rows(&d).iter().any(|p| p["pair"].as_array().unwrap().contains(&json!("root/b"))), "{}", pretty(&d));
}

#[test]
fn a_block_sliding_on_a_grooved_rail_holds_and_a_stop_in_its_travel_is_found() {
    // The block slides along the rail's grooves, its ridges in them: the
    // pair holds with t ranged; ridges too fat for the grooves fail. A stop
    // inside the travel is hit near x = 60; one past it is never near.
    let s = Scratch::new("motion-rail");
    let rail = r#"(import "std")
(t-span 0 2)
(def section (subtract (move (rect-2d 12 8) 0 4 0) (mirror-union-x (move (circle-2d 1) 6.5 5 0))))
(body rail (shape (min (subtract (extrude-x section 0 200) (pattern-linear (move (cylinder-z 1.5 20) 12.5 0 4) 8 25 0 0)) STOP)))
(body block (shape (min (subtract (move (box 30 20 8) 15 0 6) (extrude-x (- (move (rect-2d 12 8) 0 4 0) 0.2) -1 31)) (extrude-x (mirror-union-x (move (circle-2d RIDGE) 6.5 5 0)) 0 30))))
(assembly root (instance r rail) (instance b block (at) (slide x (* 40 t))))
"#;
    let slide = |name: &str, stop: &str, ridge: &str| row(s.ink(["interfere", &s.model(name, &rail.replace("STOP", stop).replace("RIDGE", ridge))]));
    let nowhere = "(move (box 1 1 1) 0 0 -50)";
    let (clear, fat) = (slide("slide.ink", nowhere, "0.7"), slide("slide-fat.ink", nowhere, "1.3"));
    assert!(verdict(&clear) == "holds" && verdict(&fat) == "fails", "{}\n{}", pretty(&clear), pretty(&fat));
    let (hit, past) = (slide("slide-stop.ink", "(move (box 4 12 14) 60 0 7)", "0.7"), slide("slide-past.ink", "(move (box 4 12 14) 115 0 7)", "0.7"));
    assert!(verdict(&hit) == "fails" && (num(&hit["witness"][0]) - 60.0).abs() < 3.0 && verdict(&past) == "holds", "{}\n{}", pretty(&hit), pretty(&past));
}

#[test]
fn how_a_claim_held_over_a_motion_is_said_with_its_moment() {
    // A block sliding 45.03 into a wall 45 away overlaps it 0.03 at the
    // end: held with t ranged only to the accuracy, its near point inside
    // both as placed at the moment it was found at.
    let s = Scratch::new("motion-held");
    let m = s.model(
        "slide.ink",
        "(import \"std\")\n(t-span 0 1)\n(body wall (shape (move (box 10 10 10) 55 0 0)))\n(body block (shape (box 10 10 10)))\n\
         (assembly root (instance w wall) (instance b block (at) (slide x (* 45.03 t))))\n\
         (assert clear (disjoint (shape-of root/w) (shape-of root/b)) (accuracy 0.1))\n",
    );
    let d = s.ink(["assert", &m]).ok();
    let a = &d["asserts"][0];
    let (t, x) = (num(&a["near"]["params"]["t"]), num(&a["near"]["at"][0]));
    assert!(verdict(a) == "holds" && a["held"] == "to the accuracy" && num(&a["near"]["depth"]) > 0.0 && 50.0 < x && x < 5.0 + 45.03 * t, "{}", pretty(&d));
}

#[test]
fn a_near_point_says_its_moment_where_the_frame_cancels_time() {
    // Two boxes on one drive, overlapping by 0.05 in x: in the first's
    // frame t cancels, but the near point as placed is at a moment, which
    // it says.
    let s = Scratch::new("motion-held-frame");
    let m = s.model(
        "pair.ink",
        "(import \"std\")\n(t-span 2 3)\n(body a (shape (box 10 10 10)))\n(body b (shape (box 10 10 10)))\n\
         (assembly root (instance a a (at) (slide x (* 10 t))) (instance b b (at (move 9.95 0 0)) (slide x (* 10 t))))\n\
         (assert clear (disjoint (shape-of root/a) (shape-of root/b)) (accuracy 0.1))\n",
    );
    let d = s.ink(["assert", &m]).ok();
    let a = &d["asserts"][0];
    let (t, x) = (num(&a["near"]["params"]["t"]), num(&a["near"]["at"][0]));
    assert!(verdict(a) == "holds" && a["held"] == "to the accuracy" && (2.0..=3.0).contains(&t) && (10.0 * t + 4.95 < x && x < 10.0 * t + 5.0), "{}", pretty(&d));
}

/// A D-shaft of radius R turning in a bore 4.3, a collar outside the block
/// wider than the bore (so no cylinder holding the shaft clears the block).
const SHAFT: &str = r#"(import "std")
(t-span 0 1)
(body shaft (shape (min (max (cylinder-z R 40) (- x 3.5)) (move (cylinder-z 8 4) 0 0 17))))
(body block (shape (subtract (box 20 20 10) (cylinder-z 4.3 12))))
(assembly root (instance blk block) (instance s shaft (at) (turn z (* 360 t))))
"#;

#[test]
fn a_d_shaft_turning_in_its_bore_holds_and_one_too_fat_fails() {
    // Its flat turns with it: the pair holds with t ranged; 0.3 too fat, it
    // fails.
    let s = Scratch::new("motion-dshaft");
    let [fits, fat] = ["4", "4.6"].map(|r| row(s.ink(["interfere", &s.model(&format!("dshaft-{r}.ink"), &SHAFT.replace("R", r))])));
    assert!(verdict(&fits) == "holds" && fits["motion"]["how"] == "t ranged" && verdict(&fat) == "fails", "{}\n{}", pretty(&fits), pretty(&fat));
}

#[test]
fn a_d_shaft_placed_off_its_joints_axis_as_cad_places_one_holds() {
    // Placed on its line as CAD places one: its axis 2.5e-3 mm off the
    // joint's and tilted 5.2e-7 rad, or 0.02 mm off: both hold with t
    // ranged.
    let s = Scratch::new("motion-dshaft-cad");
    let cad = SHAFT
        .replace("(instance s shaft (at)", "(instance s line (at)")
        .replace("(assembly root", "(assembly line (instance s shaft (at (rot 1 0 0 0.00003) (move OFF 0 0))))\n(assembly root")
        .replace("R", "4");
    let [near, far] = ["0.0025", "0.02"].map(|o| row(s.ink(["interfere", &s.model(&format!("dshaft-cad-{o}.ink"), &cad.replace("OFF", o))])));
    assert!(
        verdict(&near) == "holds" && near["motion"]["how"] == "t ranged" && verdict(&far) == "holds" && far["motion"]["how"] == "t ranged",
        "{}\n{}",
        pretty(&near),
        pretty(&far)
    );
}

#[test]
fn a_cad_rollers_key_spinning_past_a_post_is_never_held_guaranteed() {
    // A roller as CAD writes one, spinning on a swinging arm
    // (tests/motion/cad-roller.ink), with a key pointing sideways at rest:
    // spinning past a post 0.5 mm off its tube the key hits it; still, it
    // passes. The key is in the post from t = 0.4768 to 0.4797, 0.34 mm
    // deep at most: 3 ms of the motion, where the moments measured are 31
    // ms apart (0.453 and 0.484 either side of it), so they miss it. t
    // ranged reaches it (its undecided cell is the key's, in the post) but
    // halves t to 1/256 of the span, 3.9 ms in which the key sweeps 3.5 mm,
    // and ends at that first cell it cannot decide. An estimate, never
    // held: it failed where a misalignment charge held the tube still and
    // only the key moved.
    let s = Scratch::new("motion-keyed-roller");
    let keyed = r#"(import "std")
(t-span 0 1)
(body roller (shape (place (min (subtract (move (cylinder-z 12 30) -50.4 16.1 0) (move (cylinder-z 10 31) -50.4 16.1 0)) (move (box 2 2 10) -50.4 29 0)) (at (rot 0.57735 0.57735 0.57735 120)))))
(body arm (shape (move (box 100 6 4) 50 0 -20)))
(body post (shape (cylinder-z 0.5 40)))
(assembly line (instance tube roller (at (move 0 50.4003 -16.1) (rot 0.57735 0.57735 0.577351 240))))
(assembly swing (instance a arm) (instance l line (at (move 100 0 0))SPIN))
(assembly root (instance s swing (at) (turn z (* 60 t))) (instance p post (at (move 97.860871 56.5 0))))
"#;
    let [spun, still] = [("spun", " (turn z (* 3600 t))"), ("still", "")].map(|(k, spin)| {
        let d = s.ink(["interfere", &s.model(&format!("keyed-{k}.ink"), &keyed.replace("SPIN", spin))]).ok();
        rows(&d).iter().find(|p| p["pair"][1] == "root/p").cloned().unwrap_or_else(|| panic!("{}", pretty(&d)))
    });
    assert!(verdict(&spun) == "holds estimate" && verdict(&still) == "holds", "{}\n{}", pretty(&spun), pretty(&still));
}

#[test]
fn a_boss_on_a_still_stator_fails_where_a_turning_rotors_lobes_reach_it() {
    // A boss 0.66 into a 6-fold rotor's lobes over 133..167 degrees of its
    // turn, every 60 (100..200 over the span): found at a moment there.
    let s = Scratch::new("motion-boss");
    let boss = s.model(
        "boss.ink",
        r#"(import "std")
(t-span 0 1)
(body stator (shape (min (subtract (cylinder-z 30 6) (cylinder-z 20 8)) (move (box 6 6 6) 0 17.5 0))))
(body rotor (shape (max (repeat-polar (move (box 18 3 6) 7 0 0) 6) (cylinder-z 17 6))))
(assembly root (instance st stator) (instance ro rotor (at) (turn z (+ 100 (* 100 t)))))
"#,
    );
    let p = row(s.ink(["interfere", &boss]));
    assert!(verdict(&p) == "fails" && (13.0..47.0).contains(&((100.0 + 100.0 * num(&p["t"])) % 60.0)), "{}", pretty(&p));
}

#[test]
fn a_field_reading_t_only_times_a_param_at_zero_is_checked_once() {
    // It does not change over the span, so it is the claim at any moment,
    // checked once, not measured at 17 moments as an estimate; with the
    // param at 1 it is checked over the span, an estimate (a contact at
    // every moment is no range claim).
    let s = Scratch::new("motion-still-belt");
    let belt = s.model(
        "still-belt.ink",
        r#"(import "std")
(param anim 0)
(t-span 0 10)
(body belt (shape (box (+ 100 (* anim (sin t))) 6 1.4)))
(body plate (shape (move (box 50 50 4) 0 0 -2.7)))
(assert on (touches (shape-of belt) (shape-of plate)))
"#,
    );
    let still = row(s.ink(["assert", &belt]));
    let moving = row(s.ink(["assert", &belt, "--param", "anim=1"]));
    assert!(
        verdict(&still) == "holds" && !has(&still, "motion") && verdict(&moving) == "holds estimate" && moving["motion"]["how"] == "measured",
        "{}\n{}",
        pretty(&still),
        pretty(&moving)
    );
}

#[test]
fn a_belts_wrap_turning_with_its_pulley_holds_measured_and_a_slipping_one_fails() {
    // A belt of 12 teeth round a pulley its drive turns, built as std's
    // belt-along builds one: its field reads t (its teeth run round), so the
    // pair is measured at 17 moments: an estimate. A groove filled inside
    // the wrap fails; a pulley turning twice as fast as the belt meets its
    // teeth.
    let s = Scratch::new("motion-wrap");
    let wrap = r#"(import "std")
(t-span 0 1)
(def n 12)
(def p 4)
(def r (/ (* n p) (* 2 pi)))
(def P (path-loop (vec (path-arc (vec 0 0 0) (vec 0 0 1) r (vec 0 1 0) pi) (path-line (vec 0 (neg r) 0) (vec 40 (neg r) 0))
  (path-arc (vec 40 0 0) (vec 0 0 1) r (vec 0 -1 0) pi) (path-line (vec 40 r 0) (vec 0 r 0))) (vec 0 0 1)))
(def band (move (rect-2d 0.8 4) -0.6 0 0))
(def tooth (move (rect-2d 1.2 1.2) 0 0.4 0))
(body (belt (q 0)) (shape (let ((R (path-features tooth p (* q (/ (* n p) 360)))))
  (path-pieces P (lambda (c i) (min (path-sweep-at band c) (max (path-repeat-at R P i c) (- (abs (nth c 2)) 2)))) 1.5))))
(body pulley (shape (min (subtract (cylinder-z (+ r 0.1) 6) (repeat-polar (move (box 3 1.6 8) r 0 0) n) (cylinder-z 1 10)) FILL)))
(def q (* 90 t))
(assembly root (instance b (belt q)) (instance pul pulley (at) (turn z (* RATE q))))
"#;
    let nothing = "(move (box 1 1 1) 0 0 -50)";
    let [carried, filled, slips] = [("fits", nothing, "1"), ("filled", "(move (box 3 0.8 6) (neg r) 0 0)", "1"), ("slips", nothing, "2")]
        .map(|(k, f, rate)| row(s.ink(["interfere", &s.model(&format!("wrap-{k}.ink"), &wrap.replace("FILL", f).replace("RATE", rate))])));
    assert!(verdict(&carried) == "holds estimate" && carried["motion"]["how"] == "measured", "{}", pretty(&carried));
    assert!(verdict(&filled) == "fails", "{}", pretty(&filled));
    assert!(verdict(&slips) == "fails" && slips["motion"]["how"] == "measured", "{}", pretty(&slips));
}

#[test]
fn a_part_still_over_its_span_is_checked_at_a_moment_of_it_not_at_t_zero() {
    // A block slid by 10 (min 1 t) sits at 10 over [1.5, 2], past its limit
    // of 5 and 5 into the wall: checked once at a moment of the span, its
    // limits, clearance and seat fail there, by --t-span or the file's
    // span, and interfere finds the pair; slid by 5 it sits on the wall
    // (checked at t = 0 instead, limits and clear would hold, and a seat 5
    // off fail).
    let s = Scratch::new("motion-phase");
    let phase = r#"(import "std")
HEAD(body block (shape (box 10 10 10)))
(body wall (shape (move (box 10 10 10) 15 0 0)))
(assembly root (instance w wall) (instance b block (at) (slide x (* SLID (min 1 t))) (limits 0 5)))
(assert clear (disjoint (shape-of root/b) (shape-of root/w)))
(assert seated (touches (shape-of root/b) (shape-of root/w)))
"#;
    let by_flag = s.model("phase.ink", &phase.replace("HEAD", "").replace("SLID", "10"));
    let by_file = s.model("phase-file.ink", &phase.replace("HEAD", "(t-span 1.5 2)\n").replace("SLID", "10"));
    for d in [s.ink(["assert", &by_flag, "--t-span", "1.5:2"]).ok(), s.ink(["assert", &by_file]).ok()] {
        assert!(verdicts(&d) == ["fails"; 3], "{}", pretty(&d));
    }
    let touch = s.model("phase-touch.ink", &phase.replace("HEAD", "").replace("SLID", "5"));
    let seat = row(s.ink(["assert", &touch, "--t-span", "1.5:2", "--name", "seated"]));
    assert!(verdict(&seat) == "holds", "{}", pretty(&seat));
    let d = s.ink(["interfere", &by_flag, "--t-span", "1.5:2"]).ok();
    assert!(d["overlapping_pairs"] == 1 && d["pairs_apart_by_bounding_boxes"] == 0, "{}", pretty(&d));
}

/// Two blocks moved by a solve on t, as `ink assert` and `ink interfere`
/// ask them (searched in the first's frame): a at `a`, b at `b`, q = c t.
fn solved_pair(s: &Scratch, name: &str, c: f64, a: &str, b: &str) -> String {
    s.model(
        name,
        &format!(
            "(import \"std\")\n(t-span 0 1)\n(solve drive (unknowns (q 0)) (== q (* {c} t)))\n(body block (shape (box 2 2 2)))\n\
             (assembly root (instance a block (at (move {a}))) (instance b block (at (move {b}))))\n\
             (assert apart (disjoint (shape-of root/a) (shape-of root/b)))\n"
        ),
    )
}

/// A frame is built with a solved drive ranged, never as its number at the
/// view's moment: blocks at 10 t and 10 t + q + 10, q = -20 t, coincide at
/// t = 0.5, which a frame simplified with q constant called apart outright.
#[test]
fn a_frame_reads_a_solved_drive_as_ranged() {
    let s = Scratch::new("motion-solved-frame");
    let m = solved_pair(&s, "mixed.ink", -20.0, "(* 10 t) 0 0", "(+ (* 10 t) q 10) 0 0");
    for command in ["assert", "interfere"] {
        let a = row(s.ink([command, &m]));
        let t = num(&a["t"]);
        assert!(verdict(&a) == "fails" && (10.0 * t - (10.0 - 10.0 * t)).abs() < 2.0, "{command}: {}", pretty(&a));
    }
}

/// Every witness and near point found in a frame is moved into the world
/// and confirmed there with a solved drive derived at its own moment: both
/// blocks at 10 t + q, q = 100 t, the second 0.1 t up y (a witness) or
/// 1.98 - 0.02 t (an overlap under the accuracy: a near point), were
/// reported at t = 0.5 at x near 5, where both lie about x = 55.
#[test]
fn a_frames_points_are_in_the_world_at_their_moments_with_solved_drives() {
    let s = Scratch::new("motion-solved-points");
    let witness = solved_pair(&s, "witness.ink", 100.0, "(+ (* 10 t) q) 0 0", "(+ (* 10 t) q) (* 0.1 t) 0");
    let near = solved_pair(&s, "near.ink", 100.0, "(+ (* 10 t) q) 0 0", "(+ (* 10 t) q) (- 1.98 (* 0.02 t)) 0");
    // Inside a block at x 110 t and this y: below 0.
    let inside = |p: &[f64], t: f64, y: f64| (p[0] - 110.0 * t).abs().max((p[1] - y).abs()).max(p[2].abs()) - 1.0;
    for command in ["assert", "interfere"] {
        let a = row(s.ink([command, &witness]));
        let (w, t) = (nums(&a["witness"]), num(&a["t"]));
        let depth = if w.len() == 3 { inside(&w, t, 0.0).max(inside(&w, t, 0.1 * t)) } else { f64::NAN };
        assert!(verdict(&a) == "fails" && depth < -num(&a["accuracy"]), "{command}: {depth} at the witness: {}", pretty(&a));
        let a = row(s.ink([command, &near]));
        let (p, t) = (nums(&a["near"]["at"]), num(&a["near"]["params"]["t"]));
        let depth = if p.len() == 3 { inside(&p, t, 0.0).max(inside(&p, t, 1.98 - 0.02 * t)) } else { f64::NAN };
        assert!(verdict(&a) == "holds" && a["held"] == "to the accuracy" && depth < 0.0, "{command}: {depth} at the near point: {}", pretty(&a));
    }
}

/// A claim over a motion is judged in ink's range, and by its directions,
/// with what solves on t set ranged over the span, before it is searched:
/// a probe a solve carries to 2e9 mm (followed, or with more equations than
/// unknowns, measured moment by moment) is unknown, not held at the view's
/// moment's 0; and a turn whose axis a solve takes through zero at t = 0.5
/// is unknown, saying where. Where the solve is only measured, each moment
/// measured is judged so before it is asked, and a pair it carries is not
/// asked once for the span (riding) on the view's moment's judgement: an
/// axis zero at the sampled t = 0.5, two cubes riding to 2e9 mm, unknown in
/// assert and interfere. So is a claim only whose bound reads it (judged
/// before the bound is its number at the moment): a clearance, a volume of
/// shapes moving to 2e9 mm.
/// What a solve on t sets is derived from t over the claim's whole span,
/// not read through t: a turn repeating every second beside a slide
/// solved along t is not searched over one second of t standing for the
/// ten (the slide would be read only where it is in that second).
#[test]
fn a_solve_on_t_keeps_its_claim_over_the_whole_span() {
    let s = Scratch::new("motion-solve-no-drive");
    let m = s.model("slide.ink", "(t-span 0 10)\n(solve slide (unknowns (q 0)) (== q (* 10 t)))\n(assert below (value (+ (sin (deg->rad (* 360 t))) q) <= 101.5))\n");
    let a = row(s.ink(["assert", &m]));
    assert!(verdict(&a) == "holds" && !a["accuracy"].is_null() && a["param_spans"]["t"] == json!([0.0, 10.0]), "{}", pretty(&a));
}

#[test]
fn a_claim_through_a_solve_on_t_is_judged_over_the_motion_before_it_is_searched() {
    let s = Scratch::new("motion-solved-range");
    let probe = |name: &str, solve: &str| s.model(name, &format!("(import \"std\")\n(t-span 0 1)\n{solve}\n(assert inside (subset (probe q 0 0) (move (box 2 2 2) q 0 0)))\n"));
    let followed = probe("followed.ink", "(solve A (unknowns (q 0)) (== q (* 2e9 t)))");
    let measured = probe("measured.ink", "(solve A (unknowns (q 0)) (== q (* 2e9 t)) (== (* 2 q) (* 4e9 t)))");
    for m in [&followed, &measured] {
        let a = row(s.ink(["assert", m]));
        assert!(verdict(&a) == "unknown" && a["unavailable"].as_str().is_some_and(|u| u.contains("outside ink's range")), "{}", pretty(&a));
    }
    let axis = s.model(
        "axis.ink",
        "(import \"std\")\n(t-span 0 1)\n(solve A (unknowns (q -0.5)) (== q (- t 0.5)))\n(body block (shape (box 2 2 2)))\n\
         (body post (shape (move (box 2 2 2) 20 0 0)))\n(assembly root (instance a block (at (move 5 0 0) (rot q 0 0 90))) (instance p post))\n\
         (assert apart (disjoint (shape-of root/a) (shape-of root/p)))\n",
    );
    let a = row(s.ink(["assert", &axis]));
    assert!(verdict(&a) == "unknown" && a["unavailable"].as_str().is_some_and(|u| u.contains("zero length at t = 0.5")), "{}", pretty(&a));
    let squares = "(solve A (unknowns (q 1)) (== q (+ 1 t)) (== (* 2 q) (* 2 (+ 1 t))))";
    let measured_axis = s.model(
        "measured-axis.ink",
        &format!(
            "(import \"std\")\n(t-span 0 1)\n{squares}\n(body block (shape (box 2 2 2)))\n(body post (shape (move (box 2 2 2) 20 0 0)))\n\
             (assembly root (instance a block (at (move 5 0 0) (rot (* q (- t 0.5)) 0 0 60))) (instance p post))\n\
             (assert apart (disjoint (shape-of root/a) (shape-of root/p)))\n"
        ),
    );
    let riding = s.model(
        "riding.ink",
        "(import \"std\")\n(t-span 0 1)\n(solve A (unknowns (q 0)) (== q (* 2e9 t)) (== (* 2 q) (* 4e9 t)))\n(body cube (shape (box 2 2 2)))\n\
         (assembly pair (instance a cube) (instance b cube (at (move 5 0 0))))\n(assembly root (instance p pair (at (move q 0 0))))\n\
         (assert apart (disjoint (shape-of root/p/a) (shape-of root/p/b)))\n",
    );
    for (m, why) in [(&measured_axis, "zero length"), (&riding, "outside ink's range")] {
        for command in ["assert", "interfere"] {
            let a = row(s.ink([command, m]));
            assert!(verdict(&a) == "unknown" && a["unavailable"].as_str().is_some_and(|u| u.contains(why)), "{command}: {}", pretty(&a));
        }
    }
    for claim in ["(clearance a b q)", "(volume a > (* q 0.0001))"] {
        let bound = s.model(
            "bound.ink",
            &format!(
                "(import \"std\")\n(t-span 0 1)\n{squares}\n(def a (move (box 2 2 2) (* 2e9 t) 0 0))\n(def b (move (box 2 2 2) (+ 5 (* 2e9 t)) 0 0))\n\
                 (assert c {claim})\n"
            ),
        );
        let a = row(s.ink(["assert", &bound, "--res", "0.5"]));
        assert!(verdict(&a) == "unknown" && a["unavailable"].as_str().is_some_and(|u| u.contains("outside ink's range")), "{claim}: {}", pretty(&a));
    }
}

#[test]
fn pairs_moved_by_a_solve_are_ranged_or_measured_and_reported_as_decided() {
    let s = Scratch::new("motion-solved");
    // A pair moved by a solve is never culled by a swept box; over a band
    // too it is not followed, so not ranged with t: unknown.
    let d = s.ink(["interfere", &root("tests/motion/ranged-solve.ink")]).ok();
    assert!(verdicts(&d) == ["unknown"] && d["pairs_apart_by_bounding_boxes"] == 0, "{}", pretty(&d));
    // A followed solve's overlap is found with t ranged, guaranteed, and
    // each pair is reported as decided, the last of its rows the answer's
    // (a call stopped at its limit gives them): first sampled, then t
    // ranged.
    let log = s.path("slide.progress");
    let slide = root("tests/motion/solved-slide.ink");
    let d = s.ink(["interfere", &slide]).env("INK_PROGRESS", &log).ok();
    assert!(d["overlapping_pairs"] == 1 && d.get("overlapping_pairs_estimated").is_none_or(|n| n == 0), "{}", pretty(&d));
    let reported: Vec<Value> = std::fs::read_to_string(&log).unwrap().lines().map(|l| serde_json::from_str::<Value>(l).unwrap()).filter_map(|l| l.get("row").cloned()).collect();
    let last = |p: &Value| reported.iter().rev().find(|r| r["pair"] == p["pair"]).cloned();
    assert!(rows(&d).iter().all(|p| last(p).as_ref() == Some(p)) && verdict(&rows(&d)[0]) == "fails", "{reported:?}\n{}", pretty(&d));
    // One with more equations than unknowns is not followed: measured at
    // 17 moments evenly spaced, ends included (as a band's 17 values),
    // failing early, reporting its configured step.
    let squares = s.model(
        "squares.ink",
        "(import \"std\")\n(t-span 0 1)\n(solve position (unknowns (q 0)) (== q (* 100 t)) (== (* 2 q) (* 200 t)))\n(body block (shape (box 2 2 2)))\n\
         (body wall (shape (move (box 2 2 2) 50 0 0)))\n(assembly root (instance block block (at (move q 0 0))) (instance wall wall))\n\
         (assert through-wall (disjoint (shape-of root/block) (shape-of root/wall)))\n",
    );
    let a = row(s.ink(["assert", &squares]));
    let m = &a["motion"];
    assert!(
        verdict(&a) == "fails estimate"
            && m["how"] == "measured"
            && num(&m["step"]) == 1.0 / 16.0
            && m["steps"].as_u64().unwrap() < 17
            && num(&a["t"]) * 16.0 == (num(&a["t"]) * 16.0).round(),
        "{}",
        pretty(&a)
    );
}

#[test]
fn claims_over_a_motion_are_proven_t_ranged_and_those_no_range_decides_are_estimates() {
    let s = Scratch::new("motion-measured");
    // Surfaces that nearly coincide over a span, 0.01 apart: with t ranged
    // the search proves it, not measured at 17 moments as an estimate.
    let sub = s.model(
        "subset-span.ink",
        r#"(import "std")
(body ball (shape (move (sphere 4) (* 12 t) 0 0)))
(assert inside (subset (shape-of ball) (- (move (sphere 4) (* 12 t) 0 0) 0.01)))
"#,
    );
    let a = row(s.ink(["assert", &sub, "--res", "0.05", "--t-span", "0:1"]));
    assert!(verdict(&a) == "holds", "{}", pretty(&a));
    // An estimate holding (a wall thickness over a span): all_hold, not
    // all_hold_guaranteed. (At the model's res, 0.064, it took 6.6 s.)
    let thick = s.model(
        "thick-span.ink",
        r#"(import "std")
(t-span 0 1)
(body ball (shape (move (sphere 4) (* 12 t) 0 0)))
(assert thick (min-thickness (shape-of ball) > 1) (res 0.25))
"#,
    );
    let d = s.ink(["assert", &thick]).ok();
    assert!(d["all_hold"] == true && d["all_hold_guaranteed"] == false && d["counts"] == json!({"holds": 1}), "{}", pretty(&d));
    // overlaps over a span holds at the sampled moments, an estimate (an
    // overlap at every moment is no range claim; its failure when they
    // part: geometry.rs claims_run_on_a_model_view).
    let ov = s.model(
        "overlaps-span.ink",
        r#"(import "std")
(body a (shape (box 10 10 10)))
(body b (shape (move (box 10 10 10) (* 40 t) 0 0)))
(assert meet (overlaps (shape-of a) (shape-of b)))
"#,
    );
    let a = row(s.ink(["assert", &ov, "--t-span", "0:0.2"]));
    assert!(verdict(&a) == "holds estimate", "{}", pretty(&a));
    // A failure over a motion, t ranged, is a witness at a moment:
    // guaranteed, to the claim's accuracy; holding over the motion, proven.
    let out = s.model(
        "slide-out.ink",
        r#"(import "std")
(body box (shape (box 10 10 10)))
(body ball (shape (sphere 1)))
(assembly root (instance c box (at (move 0 0 0))) (instance s ball (at (move 0 0 0)) (slide x (* 8 t))))
(assert inside (subset (shape-of root/s) (shape-of root/c)))
"#,
    );
    let a = row(s.ink(["assert", &out, "--t-span", "0:1"]));
    assert!(a["status"] == "fails" && a["accuracy"] == 0.05 && num(&a["t"]) > 0.5 && num(&a["t"]) <= 1.0, "{}", pretty(&a));
    let a = row(s.ink(["assert", &out, "--t-span", "0:0.3"]));
    assert!(a["status"] == "holds" && a["accuracy"] == 0.05, "{}", pretty(&a));
}

#[test]
fn value_claims_speed_limits_and_relation_claims_hold_to_their_accuracy_and_are_unknown_without() {
    // A value meeting its bound (10 sin(2 pi t) < 10 at t = 0.25; p (2 - p)
    // <= 1 at p = 1; a turn at 10 sin(2 pi t) deg/s against speed limits of
    // 10) is unknown asked exactly over ranged params, measured over a
    // motion where the range search is undecided (the wave's moments
    // include t = 0.25, where it meets its bound: unknown; the speed's and
    // the relation's, holding at each, an estimate), and holds to
    // (accuracy x); a violation deeper than x fails. A number in accuracy
    // only where it holds or fails, guaranteed; an unknown says the
    // accuracy it was undecided at apart.
    let s = Scratch::new("motion-value-accuracy");
    let wave = "(* 10 (sin (* 6.283185307179586 t)))";
    let turn = "(turn z (* 10 (* 0.15915494309189535 (- 1 (cos (* 6.283185307179586 t))))))";
    let tv = s.model(
        "value-accuracy.ink",
        &format!(
            r#"(import "std")
(t-span 0 1)
(param p 1 (tol 0 1))
(body ball (shape (sphere 1)))
(assembly exact (instance w ball (at (move 0 0 0)) {turn} (speed-limits -10 10)))
(assembly tolerant (instance w ball (at (move 0 0 0)) {turn} (speed-limits -10 10 (accuracy 0.01))))
(assembly related (instance a ball (at (move 0 0 0))) (instance b ball (at (move 5 0 0)))
  (relate a b (lambda (fa fb sa sb) (relation demo (claim swing {wave} -10 10) (claim swing-to {wave} -10 10 (accuracy 0.01))))))
(assert exact (value {wave} < 10))
(assert to-0.01 (value {wave} < 10) (accuracy 0.01))
(assert past-0.01 (value {wave} < 9.95) (accuracy 0.01))
(assert ranged-exact (value (* p (- 2 p)) <= 1))
(assert ranged (value (* p (- 2 p)) <= 1) (accuracy 0.001))
"#
        ),
    );
    let d = s.ink(["assert", &tv]).ok();
    let opt = |v: &Value| v.as_f64().map_or("-".to_string(), |x| x.to_string());
    let got: std::collections::BTreeMap<String, String> = rows(&d)
        .iter()
        .map(|a| (a["name"].as_str().unwrap().to_string(), format!("{} {} {}", a["status"].as_str().unwrap(), opt(&a["accuracy"]), opt(&a["accuracy_asked"]))))
        .collect();
    let want: std::collections::BTreeMap<String, String> = [
        ("speed-limits:exact/w", "holds - -"),
        ("speed-limits:tolerant/w", "holds 0.01 -"),
        ("relate:related/a~b/swing", "holds - -"),
        ("relate:related/a~b/swing-to", "holds 0.01 -"),
        ("exact", "unknown - -"),
        ("to-0.01", "holds 0.01 -"),
        ("past-0.01", "fails 0.01 -"),
        ("ranged-exact", "unknown - 0"),
        ("ranged", "holds 0.001 -"),
    ]
    .map(|(k, v)| (k.to_string(), v.to_string()))
    .into();
    assert_eq!(got, want, "{}", pretty(&d));
    let past = rows(&d).iter().find(|a| a["name"] == "past-0.01").unwrap();
    assert!(num(&past["value"]) >= 9.96 && num(&past["t"]) > 0.2 && num(&past["t"]) < 0.3, "{}", pretty(past));
    assert!(d["all_hold"] == false && d["all_hold_guaranteed"] == false && d["counts"] == json!({"holds": 6, "fails": 1, "unknown": 2}), "{}", pretty(&d["counts"]));
}

#[test]
fn a_values_accuracy_is_the_models_own_and_a_volumes_an_estimate_none() {
    let s = Scratch::new("motion-accuracy");
    let vu = s.model(
        "volume-accuracy.ink",
        r#"(import "std")
(body b (shape (box 10 10 10)))
(assert vol (volume (shape-of b) > 999) (accuracy 1mm3))
(assert v (value 2 < 3) (accuracy 0.1))
"#,
    );
    let a: std::collections::BTreeMap<String, Value> = rows(&s.ink(["assert", &vu]).ok()).iter().map(|x| (x["name"].as_str().unwrap().to_string(), x.clone())).collect();
    assert!(a["v"]["accuracy"] == 0.1 && a["vol"]["accuracy"].is_null(), "{}\n{}", pretty(&a["v"]), pretty(&a["vol"]));
}

/// An arm turned past a post 0.64 off its tip (at --res 0.1), a 0.2
/// clearance claimed; POST places the post.
const NEAR_PASS: &str = r#"(import "std")
(body arm (shape (move (box 35 4 4) 17.5 0 0)))
(body post (shape (move (cylinder-z 0.5 20) 36.2 0 0)))
(assembly root (instance arm arm (at (rot-z (* 360 t)))) POST)
(assert clear (clearance (shape-of root/arm) (shape-of root/post) 0.2))
"#;

#[test]
fn a_near_pass_holds_t_ranged_the_post_still_or_turning() {
    // A clearance over a motion is asked with t ranged (as interfere asks
    // it, over t as its input), the post still or turning too.
    let s = Scratch::new("motion-near-pass");
    let still = s.model("motion2.ink", &NEAR_PASS.replace("POST", "(instance post post)"));
    let turning = s.model("motion2-turning.ink", &NEAR_PASS.replace("POST", "(instance post post (at (rot-z (* 10 t))))"));
    let a = row(s.ink(["assert", &still, "--t-span", "0:1", "--res", "0.1"]));
    assert!(verdict(&a) == "holds" && a["motion"]["how"] == "t ranged", "{}", pretty(&a));
    let a = row(s.ink(["assert", &turning, "--t-span", "0:1", "--res", "0.1"]));
    assert!(verdict(&a) == "holds", "{}", pretty(&a));
}

#[test]
fn an_arm_turning_about_its_own_axis_is_apart_by_its_boxes_or_proven_with_t_ranged() {
    // The arm's corner sweeps a circle of radius 35.06, the post at 35.7.
    // On the x axis the arm's boxes over pieces of the turn are apart from
    // the post already; at 45 degrees they meet it, and t ranged decides.
    let s = Scratch::new("motion-envelope");
    let on_axis = s.model("motion2.ink", &NEAR_PASS.replace("POST", "(instance post post)"));
    let d = s.ink(["interfere", &on_axis, "--t-span", "0:1"]).ok();
    assert!(rows(&d).is_empty() && d["pairs_apart_by_bounding_boxes"] == 1, "{}", pretty(&d));
    let at45 = s.model("motion3.ink", &NEAR_PASS.replace("POST", "(instance post post)").replace("36.2 0 0", "25.5973 25.5973 0"));
    let d = s.ink(["interfere", &at45, "--t-span", "0:1"]).ok();
    let p = &rows(&d)[0];
    assert!(rows(&d).len() == 1 && verdict(p) == "holds" && p["motion"]["how"] == "t ranged", "{}", pretty(&d));
}

#[test]
fn parts_following_each_other_on_a_track_are_apart_by_their_boxes_over_pieces_of_the_motion() {
    // Their boxes over the whole motion meet over most of it, their boxes
    // over pieces of it at no moment, so it is no moving pair searched to
    // hold. The one behind bobbing as it goes (its own field reads t, as a
    // belt's teeth do): its box over the motion in its own frame, moved by
    // its placement over the pieces.
    let s = Scratch::new("motion-track");
    let track = s.model(
        "track.ink",
        r#"(import "std")
(t-span 0 1)
(body block (shape (box 2 2 2)))
(assembly root (instance lead block (at (move (* 100 t) 0 0))) (instance next block (at (move (- (* 100 t) 10) 0 0))))
"#,
    );
    let bob = s.model(
        "track-bob.ink",
        r#"(import "std")
(t-span 0 1)
(body block (shape (box 2 2 2)))
(body bob (shape (move (box 2 2 2) 0 0 (* 0.1 (sin (* 360 t))))))
(assembly root (instance lead block (at (move (* 100 t) 0 0))) (instance next bob (at (move (- (* 100 t) 10) 0 0))))
"#,
    );
    for m in [track, bob] {
        let d = s.ink(["interfere", &m]).ok();
        assert!(rows(&d).is_empty() && d["pairs_apart_by_bounding_boxes"] == 1, "{m}:\n{}", pretty(&d));
    }
}

#[test]
fn interfere_counts_the_pairs_holding_as_estimates() {
    // A part whose own field moves is measured at the motion's moments
    // (beside other pairs t ranged is not tried where a shape's own field
    // reads t): a 4 mm box wobbling on an arm turning 130 degrees passes
    // through a post between them and holds, an estimate, which the
    // summary counts (overlapping_pairs alone read 0).
    let s = Scratch::new("motion-swing-bob");
    let swing = s.model(
        "swing-bob.ink",
        r#"(import "std")
(t-span 0 1)
(body bob (shape (move (box 4 4 4) (+ 300 (* 0.5 (sin (* 720 t)))) 0 0)))
(body post (shape (move (cylinder-z 1 20) 0 300 0)))
(assembly root (instance arm bob (at (rot-z (* 130 t)))) (instance post post))
"#,
    );
    let d = s.ink(["interfere", &swing]).ok();
    assert!(d["overlapping_pairs"] == 0 && d["holding_pairs_estimated"] == 1 && verdicts(&d) == ["holds estimate"], "{}", pretty(&d));
}

#[test]
fn a_ball_whose_field_reads_t_is_checked_against_a_clearance_t_ranged() {
    // Its shape reads t (not its placement): passing a post 1 mm away (at
    // t 0.73) it fails a 2 mm clearance; 3 mm away it holds, t ranged.
    let s = Scratch::new("motion-sweep");
    for (gap, want) in [(4, "fails"), (6, "holds")] {
        let sw = s.model(
            &format!("sweep-{gap}.ink"),
            &format!(
                r#"(import "std")
(body ball (shape (move (sphere 2) (* 100 (- t 0.5)) 0 0)))
(body post (shape (move (cylinder-z 1 10) 23 {gap} 0)))
(assembly root (instance ball ball) (instance post post))
(assert clear (clearance (shape-of root/ball) (shape-of root/post) 2))
"#
            ),
        );
        let a = row(s.ink(["assert", &sw, "--t-span", "0:1", "--res", "0.05"]));
        let at = if want == "fails" { a["t"].as_f64().unwrap_or(0.73) } else { 0.73 };
        assert!(verdict(&a) == want && (at - 0.73).abs() < 0.05, "gap {gap}: {}", pretty(&a));
    }
}

#[test]
fn a_pin_joint_holds_over_the_eyes_motion_on_the_pin_and_an_oversized_pin_fails_there() {
    // A link whose eye (r 3.4) rides on a turning crank's pin (r 3), both
    // moving 63 mm per unit t: in the crank's frame the pin rests and the
    // eye turns about it, proven with t ranged. A quarter of the crank's
    // turn.
    let s = Scratch::new("motion-pin");
    let pj = s.model(
        "pin.ink",
        r#"(import "std")
(param hole 3.4)
(def th (* 2 pi t)) (def bx (* 10 (cos th))) (def by (* 10 (sin th)))
(body crank (shape (min (extrude-z (move (rect-2d 16 8) 5 0 0) 0 4) (move (cylinder-z 3 6) 10 0 6))))
(body link (shape (max (extrude-z (move (rect-2d 50 10) 20 0 0) 5 8) (- hole (len2 x y)))))
(assembly root (instance crank crank (at (rot-z (* 360 t))))
  (instance link link (at (rot-z (* (/ 180 pi) (atan2 (neg by) 60))) (move bx by 0))))
"#,
    );
    let d = s.ink(["interfere", &pj, "--t-span", "0:0.25"]).ok();
    assert!(d["overlapping_pairs"] == 0 && rows(&d).iter().all(|p| verdict(p) == "holds"), "{}", pretty(&d));
    // Its witness is in world space, at the pin.
    let p = row(s.ink(["interfere", &pj, "--t-span", "0:0.25", "--param", "hole=2.7"]));
    let (w, t) = (nums(&p["witness"]), num(&p["t"]));
    let pin = [10.0 * (2.0 * std::f64::consts::PI * t).cos(), 10.0 * (2.0 * std::f64::consts::PI * t).sin()];
    let r = if w.len() == 3 { (w[0] - pin[0]).hypot(w[1] - pin[1]) } else { f64::NAN };
    assert!(verdict(&p) == "fails" && (2.7..=3.0).contains(&r), "{r}: {}", pretty(&p));
}

#[test]
fn a_hole_turning_about_its_axis_round_a_pin_holds_and_off_its_axis_hits_it() {
    // The hole (radius 3.2) turns about the pin (radius 3): held with t
    // ranged; set 0.5 off its axis of turning it hits the pin by 0.3, at a
    // moment found.
    let s = Scratch::new("motion-turning-hole");
    for (ecc, want) in [("0", "holds"), ("0.5", "fails")] {
        let rh = s.model(
            &format!("rotating-hole-{ecc}.ink"),
            &format!(
                r#"(import "std")
(body plate (shape (subtract (box 10 10 4) (move (cylinder-z 3.2 100) {ecc} 0 0))))
(body pin (shape (cylinder-z 3 20)))
(assembly root (instance plate plate (at (rot-z (* 360 t)))) (instance pin pin))
(assert clear (disjoint (shape-of root/pin) (shape-of root/plate)))
"#
            ),
        );
        let a = row(s.ink(["assert", &rh, "--t-span", "0:1"]));
        let ok = verdict(&a) == want && (want == "holds" || (0.0..=1.0).contains(&num(&a["t"])));
        assert!(ok, "{ecc} off: {}", pretty(&a));
    }
}

#[test]
fn a_planet_placed_by_t_in_a_ring_is_asked_over_t_as_its_input() {
    // Orbiting and spinning with no joint: t is the pair's input, each
    // piece of it proven over all its moments at once, not stepped (66
    // steps; a planetary gearbox, 520 s and unknown). Asserted disjoint it
    // is asked as interfere asks it, not stepped over time.
    let s = Scratch::new("motion-planet");
    for (orbit, want) in [("12", "holds"), ("13.5", "fails")] {
        let pl = s.model(
            &format!("planet-{orbit}.ink"),
            &format!(
                r#"(import "std")
(t-span 0 1)
(body ring (shape (subtract (cylinder-z 30 6) (cylinder-z 20 8))))
(body planet (shape (subtract (min (cylinder-z 5 4) (move (box 4 1 4) 5 0 0)) (cylinder-z 1 8))))
(assembly root (instance r ring) (instance p planet (at (rot-z (* -720 t)) (move {orbit} 0 0) (rot-z (* 72 t)))))
(assert apart (disjoint (shape-of root/r) (shape-of root/p)))
"#
            ),
        );
        let p = row(s.ink(["interfere", &pl]));
        assert!(verdict(&p) == want && !has(&p, "over"), "{orbit} out: {}", pretty(&p));
        let a = row(s.ink(["assert", &pl]));
        let ok = verdict(&a) == want && (want == "holds" || (0.0..=1.0).contains(&num(&a["t"])));
        assert!(ok, "{orbit} out, asserted: {}", pretty(&a));
    }
}

/// A lever 120 long on a pin PIN (its top touching the lever's hub) in a
/// bore BORE, POST on the base, FRAME after the lever's shape.
fn lever(pin: &str, post: &str, bore: &str, frame: &str) -> String {
    format!(
        r#"(import "std")
(t-span 0 1)
(body base (shape (min (move (box 40 40 4) 0 0 -12) (move (cylinder-z {pin} 12) 0 0 -4){post})))
(body lever (shape (subtract (min (move (box 120 12 4) 60 0 0) (move (cylinder-z 6 4) 0 0 4)) (cylinder-z {bore} 4.0001))){frame})
"#
    )
}

const LEVER_JOINT: &str = "(assembly root (instance b base) (instance l lever (at) (turn z (* 120 t))))\n";

/// The lever on its pin, as placed: on its joint, met at a frame turned
/// over 2 mm off its origin, or in an arm placed at 30 degrees.
fn levers(tag: &str, pin: &str, post: &str, want: &str) {
    // The bore and the hub are round about the joint's axis: in the pin's
    // frame its turn moves nothing of them. A post 60 mm out on y is met by
    // the arm's side from 82.3 degrees to 97.7 of its turn (120 t), found at
    // a moment there.
    let s = Scratch::new(tag);
    let offset = ("offset", "\n (frame pivot (at (rot-x 180) (move 0 0 2)))", "(assembly root (instance b base) (instance l lever (at (move 0 0 2)) pivot (turn z (* 120 t))))\n");
    let tilted = ("tilted", "", "(assembly arm (instance b base) (instance l lever (at) (turn z (* 120 t))))\n(assembly root (instance a arm (at (rot-x 30))))\n");
    for (name, frame, asm) in [("joint", "", LEVER_JOINT), offset, tilted] {
        let p = row(s.ink(["interfere", &s.model(&format!("lever-{name}.ink"), &(lever(pin, post, "3.6", frame) + asm))]));
        let ok = match want {
            "holds" => verdict(&p) == "holds",
            "fails" => verdict(&p) == "fails",
            _ => verdict(&p) == "fails" && (82.3..97.7).contains(&(120.0 * num(&p["t"]))),
        };
        assert!(ok, "a lever on a pin {pin}{post}, {name}: want {want}: {}", pretty(&p));
    }
}

#[test]
fn a_long_lever_on_a_pin_0_6_clear_of_its_bore_holds_however_placed() {
    levers("motion-lever-holds", "3", "", "holds");
}

#[test]
fn a_long_lever_on_a_pin_cutting_its_bore_fails_however_placed() {
    // A pin 0.4 oversize cuts the bore 0.2 deep (four times interfere's accuracy).
    levers("motion-lever-fails", "4", "", "fails");
}

#[test]
fn a_long_lever_meets_a_post_where_its_turn_carries_it_however_placed() {
    levers("motion-lever-post", "3", " (move (cylinder-z 2 8) 0 60 -4)", "at the post");
}

#[test]
fn a_bore_cut_short_by_twice_the_accuracy_is_unknown_at_once_saying_where() {
    // The overlap is as deep as the accuracy (the model's, 0.05) at every
    // angle: undecided, saying where.
    // Asked over t as one of two bodies moving apart, the same.
    let s = Scratch::new("motion-press");
    let two = "(body block (shape (box 4 4 4)))\n(assembly arms (instance l lever (at) (turn z (* 120 t))) (instance k block (at (rot-z (* -90 t)) (move 0 -40 0))))\n\
               (assembly root (instance b base) (instance m arms))\n(assert apart (disjoint (shape-of root/b) (shape-of root/m)))\n";
    for (name, cmd, bore, text) in [("on its joint", "interfere", "2.9", LEVER_JOINT), ("over t", "assert", "2.9", two)] {
        let p = row(s.ink([cmd, &s.model(&format!("press-{cmd}.ink"), &(lever("3", "", bore, "") + text))]));
        assert!(verdict(&p) == "unknown" && has(&p["motion"], "undecided"), "{name}: {}", pretty(&p));
    }
}

#[test]
fn an_answer_over_a_motion_is_the_same_on_any_threads() {
    // A disc 0.2 off its axis turned by t past a bump: its pieces' searches
    // meet cells that stay blurred, and a witness after one in the search's
    // order would be kept or not as the threads came to them (14, 39, 38 or
    // 13 pieces, run to run; a 2% oval disc beside it, the same). A search's
    // outcomes are taken in its order (ink-kernel search.rs), so the answer
    // is the question's. The span is cut to the first contact (0.18..0.25,
    // three pieces: the whole turn takes 19.5 s on one thread).
    let s = Scratch::new("motion-discs");
    let discs = s.model(
        "discs.ink",
        r#"(import "std")
(t-span 0.18 0.25)
(body off (shape (move (cylinder-z 10 4) 0.2 0 0)))
(body bump (shape (move (box 2 0.2 4) 0 10.15 0)))
(assembly root (instance f off (at (rot-z (* 360 t)))) (instance b bump))
(assert off (disjoint (shape-of root/f) (shape-of root/b)))
"#,
    );
    let answer = |threads: &str| {
        let a = row(s.ink(["assert", &discs]).env("INK_NO_CACHE", "1").env("RAYON_NUM_THREADS", threads));
        assert!(a["status"] == "fails" && a["motion"]["steps"].as_u64().unwrap() > 1, "{threads} threads: {}", pretty(&a));
        json!([a["status"], a["t"], a["motion"]["steps"], a["witness"], a["value"]])
    };
    let one = answer("1");
    for threads in ["16", "16"] {
        assert_eq!(answer(threads), one, "{threads} threads against 1");
    }
}

#[test]
fn a_pulley_turned_by_t_about_its_axis_on_a_shaft_its_bore_touches_holds() {
    // Its WHERE turns it, no joint (a belt's phase): in the pulley's frame
    // the shaft turns about its own axis, round there, and its turn written
    // as two half turns is the same.
    let s = Scratch::new("motion-own-turn");
    let ot = r#"(import "std")
(t-span 0 1)
(body shaft (shape (min (move (box 40 40 4) 0 0 -12) (move (cylinder-z 2.5 30) 0 0 5))))
(body pulley (shape (subtract (min (move (cylinder-z 8 6) 0 0 3) (move (box 40 4 2) 20 0 3)) (cylinder-z 2.5 30))))
(assembly root (instance s shaft (at (move 10 0 0))) (instance p pulley (at (rot-z (* 720 t)) (move 10 0 0))))
"#;
    for (name, text) in [("own-turn.ink", ot.to_string()), ("own-turn-twice.ink", ot.replace("(rot-z (* 720 t))", "(rot-z (* 360 t)) (rot-z (* 360 t))"))] {
        let p = row(s.ink(["interfere", &s.model(name, &text)]));
        assert!(verdict(&p) == "holds", "{name}: {}", pretty(&p));
    }
}

#[test]
fn a_shaft_tilted_by_its_placement_in_a_still_bore_is_one_still_question_about_any_axis() {
    // A shaft as wide as its bore, turned about its own axis by a rot-y on
    // its carrier (no joint), in a still block with a set screw 0.4 off the
    // bore. In the block's frame the shaft's turn moves nothing of it, about
    // y as about z: the claim is asked once, reading no t (no t in its
    // spans), and proven. In the shaft's frame the screw sweeps round the
    // contact: 166 s and left to samples. Oversize, it fails.
    let s = Scratch::new("motion-tilted-shaft");
    let model = |axis: &str, r: &str| {
        let (bore, turn) = match axis {
            "y" => ("(move (cylinder-y 4 12) 0 0 15)", "(rot-x -90)"),
            _ => ("(move (cylinder-z 4 12) 0 0 15)", "(move 0 0 0)"),
        };
        let screw = if axis == "y" { "(move (cylinder-z 1.5 2) 0 0 20.4)" } else { "(move (cylinder-x 1.5 2) 5.4 0 15)" };
        let block = if axis == "y" { "(box 30 10 20)" } else { "(box 30 20 10)" };
        format!(
            "(import \"std\")\n(t-span 0 12)\n(def tilt (* 35 (cos (* 30 t))))\n\
             (body block (shape (min (subtract (move {block} 0 0 15) {bore}) {screw})))\n(body shaft (shape (cylinder-z {r} 300)))\n\
             (assembly carrier (instance s shaft (at {turn})))\n\
             (assembly root (instance b block (at (move 0 0 25))) (instance c carrier (at (rot-{axis} tilt) (move 0 0 40))))\n"
        )
    };
    for axis in ["y", "z"] {
        let p = row(s.ink(["interfere", &s.model(&format!("tilted-{axis}.ink"), &model(axis, "4"))]));
        let still = p["motion"]["how"] == "t ranged" && !has(&p["param_spans"], "t");
        assert!(verdict(&p) == "holds" && still, "about {axis}: {}", pretty(&p));
        let p = row(s.ink(["interfere", &s.model(&format!("tilted-{axis}-fat.ink"), &model(axis, "4.2"))]));
        assert!(verdict(&p) == "fails", "oversize, about {axis}: {}", pretty(&p));
    }
}

#[test]
fn a_pin_a_turning_hole_leaves_fails_whichever_way_the_turn_is_shared() {
    // The pin given the disc's own turn (not nested in it), or its carrier
    // turned so, or both on joints with the same drive: only the pin's own
    // link is at rest in the parts its turn leaves unchanged; its turn's
    // nodes are shared by the others' (they held over the motion).
    let s = Scratch::new("motion-shared-turn");
    let disc = "(body disc (shape (subtract (cylinder-z 30 4) (move (cylinder-z 3 10) 10 0 0))))\n(body pin (shape (cylinder-z 2 20)))\n";
    for (name, asm, pin) in [
        ("own-turn-shared", "(assembly root (instance d disc (at (rot-z (* 90 t)))) (instance p pin (at (rot-z (* 90 t)) (move 10 0 0))))\n", "p"),
        (
            "own-turn-carrier",
            "(assembly carrier (instance p pin (at (rot-z (* 90 t)) (move 10 0 0))))\n(assembly root (instance d disc) (instance c carrier (at (rot-z (* 90 t)))))\n",
            "c/p",
        ),
        ("joints-same-drive", "(assembly root (instance d disc (at) (turn z (* 90 t))) (instance p pin (at (move 10 0 0)) (turn z (* 90 t))))\n", "p"),
    ] {
        let sh = s.model(&format!("{name}.ink"), &format!("(import \"std\")\n(t-span 0 1)\n{disc}{asm}(assert apart (disjoint (shape-of root/d) (shape-of root/{pin})))\n"));
        let (p, a) = (row(s.ink(["interfere", &sh])), row(s.ink(["assert", &sh])));
        assert!(verdict(&p) == "fails" && verdict(&a) == "fails" && num(&p["t"]) > 0.0 && num(&p["t"]) <= 1.0, "{name}:\n{}\n{}", pretty(&p), pretty(&a));
    }
}

/// An arm on a crank's joint (held to one turn by its limits) and a post
/// at 90 degrees.
const CRANK: &str = r#"(import "std")
(body arm (shape (move (box 20 2 2) 10 0 0)))
(body post (shape (move (cylinder-z 1 10) 0 15 0)))
(assembly (crank (a 0)) (instance post post) (instance arm arm (at) (turn z a) (limits -180 180)))
"#;

/// A design placing ASM (k.crank, or another of crank.ink's) driven by
/// `drive` at `at` over [0, span], HEAD first.
fn cranked(s: &Scratch, name: &str, head: &str, asm: &str, drive: &str, at: &str, span: &str) -> String {
    s.model(name, &format!("(import \"std\")\n{head}(import \"crank.ink\" as k)\n(t-span 0 {span})\n(assembly root (instance c ({asm} {drive}) (at {at})))\n"))
}

#[test]
fn a_pair_on_one_joint_clashes_where_the_drive_reaches_the_post() {
    // The post is at 80..100 degrees of the arm's turn: driven past it, the
    // pair fails at a moment there; driven to 60 degrees, it holds.
    let s = Scratch::new("motion-crank");
    s.model("crank.ink", CRANK);
    for (speed, want) in [(360.0, "fails"), (120.0, "fails"), (60.0, "holds")] {
        let p = row(s.ink(["interfere", &cranked(&s, &format!("crank-{speed}.ink"), "", "k.crank", &format!("(* {speed} t)"), "(rot-x 20) (move 3 4 5)", "1")]));
        let ok = verdict(&p) == want && (want == "holds" || (80.0..100.0).contains(&(speed * num(&p["t"]) % 360.0)));
        assert!(ok, "driven at {speed} deg/s: want {want}: {}", pretty(&p));
    }
}

#[test]
fn a_relation_without_a_zone_leaves_a_jointed_pairs_clash_and_one_with_a_zone_takes_it_off() {
    // The row carries the relation; a zone holding the clash leaves it out.
    let s = Scratch::new("motion-crank-relation");
    for (name, rel) in [("crank-rel", "(relation near (data side 1))"), ("crank-zone", r#"(relation hit (zone "the post" (move (cylinder-z 1.5 12) 0 15 0)))"#)] {
        s.model(
            &format!("{name}-k.ink"),
            &format!(
                r#"(import "std")
(body arm (shape (move (box 20 2 2) 10 0 0)))
(body post (shape (move (cylinder-z 1 10) 0 15 0)))
(assembly (crank (a 0)) (instance post post) (instance arm arm (at) (turn z a) (limits -180 180))
  (relate post arm (lambda (fa fb sa sb) {rel})))
"#
            ),
        );
        let m = s.model(&format!("{name}.ink"), &format!("(import \"std\")\n(import \"{name}-k.ink\" as k)\n(t-span 0 1)\n(assembly root (instance c (k.crank (* 360 t))))\n"));
        let p = row(s.ink(["interfere", &m]));
        let kind = &p["relations"][0]["kind"];
        let ok = if name == "crank-rel" { verdict(&p) == "fails" && kind == "near" } else { verdict(&p) != "fails" && kind == "hit" };
        assert!(ok, "{name}: {}", pretty(&p));
    }
}

#[test]
fn moving_pairs_are_searched_side_by_side_each_row_its_own() {
    // Each pair is first answered at sampled moments, then searched with t
    // ranged beside the others as the runners free: an arm's fin turning
    // clear of a post inside its sweep holds, and a pin on the arm at the
    // start fails there, each its own row. A cache of its own, so both are
    // searched.
    let s = Scratch::new("motion-side-by-side");
    let side = s.model(
        "side.ink",
        r#"(import "std")
(t-span 0 1)
(body arm (shape (min (move (box 50 2 2) 25 0 0) (move (box 2 2 12) 49 0 5))))
(body post (shape (move (cylinder-z 0.5 4) 0 0 6)))
(body pin (shape (cylinder-z 0.5 4)))
(assembly root (instance arm arm (at) (turn z (* 360 t))) (instance in post (at (move 47.3 0 0))) (instance on pin (at (move 40 0 0))))
"#,
    );
    let log = s.path("side.progress");
    let d = s.ink(["interfere", &side]).env("RAYON_NUM_THREADS", "2").env("INK_PROGRESS", &log).env("XDG_CACHE_HOME", s.path("cache")).env_remove("INK_NO_CACHE").ok();
    let lines: Vec<Value> = std::fs::read_to_string(&log).unwrap().lines().map(|l| serde_json::from_str(l).unwrap()).collect();
    let line = |key: &str, p: &str| lines.iter().any(|l| l[key] == json!(["root/arm", p]));
    let (inside, on) = (pair(&d, "root/arm", "root/in"), pair(&d, "root/arm", "root/on"));
    assert!(["root/in", "root/on"].iter().all(|p| line("running", p) && line("done", p)), "{lines:?}");
    assert!(verdict(inside) == "holds" && verdict(on) == "fails" && (0.0..0.01).contains(&num(&on["t"])), "{}", pretty(&d));
}

#[test]
fn copies_of_a_sub_assembly_searched_side_by_side_ask_one_question() {
    // Each copy's post is asked in its arm's frame, 47.25 + 200 k less
    // 200 k: exact, so one question with t ranged, the same answer for
    // every copy, and a second call answers each from the cache. (A copy's
    // first answer is sampled in the world, its own question: so a row's
    // `cached` says nothing of whether its t-ranged search waited on
    // another copy's in flight.)
    let s = Scratch::new("motion-copies");
    let copies = s.model(
        "copies.ink",
        &format!(
            r#"(import "std")
(t-span 0 1)
(body arm (shape (min (move (box 50 2 2) 25 0 0) (move (box 2 2 12) 49 0 5))))
(body post (shape (move (cylinder-z 0.5 4) 0 0 6)))
(assembly unit (instance arm arm (at) (turn z (* 360 t))) (instance in post (at (move 47.25 0 0))))
(assembly root{})
"#,
            (0..4).map(|k| format!(" (instance u{k} unit (at (move {} 0 0)))", 200 * k)).collect::<String>()
        ),
    );
    let run = || s.ink(["interfere", &copies]).env("RAYON_NUM_THREADS", "2").env("XDG_CACHE_HOME", s.path("cache")).env_remove("INK_NO_CACHE").ok();
    let (first, again) = (run(), run());
    let moving = |d: &Value| -> Vec<Value> { rows(d).iter().filter(|p| has(p, "motion")).cloned().collect() };
    assert!(moving(&first).len() == 4 && moving(&first).iter().all(|p| verdict(p) == "holds" && p["motion"]["how"] == "t ranged"), "{}", pretty(&first));
    assert!(moving(&again).iter().all(|p| p["cached"] == true && verdict(p) == "holds"), "{}", pretty(&again));
}

#[test]
fn a_stepper_rotor_and_a_pulley_on_its_shaft_hold() {
    let s = Scratch::new("motion-stepper");
    // A pulley placed on the rotor's shaft frame: the joint turns both, so
    // in the rotor's frame it does not move, and the pair in contact is
    // proven with t ranged.
    let pul = s.model(
        "pulley-rotor.ink",
        &format!(
            r#"(import "std")
(import "{}" as stepper)
(import "{}" as pulley)
(t-span 0 1)
(assembly (drv (a 0)) (instance rotor stepper.rotor (at) (turn z a)) (frame shaft (frame-of rotor shaft)))
(assembly root (instance m (drv (* 0.5 600rpm (square t)))) (instance p pulley.part (at (move 0 0 -16) (frame-of m shaft)) hub))
"#,
            fixture("motor.ink"),
            fixture("pulley.ink")
        ),
    );
    let p = row(s.ink(["interfere", &pul]));
    assert!(verdict(&p) == "holds" && p["motion"]["how"] == "t ranged", "{}", pretty(&p));
    // The 17HS4401S's rotor in its stator, 0.12 radial at the shaft's bore:
    // the rotor turns about its own axis.
    let st = s.model(
        "rotor-stator.ink",
        &format!("(import \"std\")\n(import \"{}\" as stepper)\n(t-span 0 1)\n(assembly root (instance m (stepper.drive (* 360 t))))\n", fixture("motor.ink")),
    );
    let d = s.ink(["interfere", &st]).ok();
    let p = pair(&d, "root/m/stator", "root/m/rotor");
    assert!(verdict(p) == "holds", "{}", pretty(p));
}

#[test]
fn parts_fixed_together_riding_a_motion_are_a_still_question() {
    // A plate and the block resting on it, turned and tipped by t written
    // in their assembly's WHERE (no joint): one rigid group, the pair asked
    // once in its frame (no motion block; not one step, nor 106 by how far
    // the plate's corners travel) with the answer it has at any moment.
    // The plate passing a post moves relative to it, and is apart from it
    // by its boxes over pieces of the motion (its corner 14.14 from the
    // axis, the post from 14.3), or held over the motion; the block sunk
    // 0.5 into the plate fails.
    let s = Scratch::new("motion-rigid-group");
    for (sunk, want) in [(0.0, "holds"), (0.5, "fails")] {
        let rg = s.model(
            &format!("rigid-group-{sunk}.ink"),
            &format!(
                r#"(import "std")
(t-span 0 1)
(body plate (shape (box 20 20 2)))
(body block (shape (move (box 6 6 6) 0 0 {})))
(body post (shape (move (cylinder-z 0.3 20) 14.6 0 0)))
(assembly payload (instance plate plate) (instance block block))
(assembly root (instance post post) (instance p payload (at (rot-x (* 20 t)) (rot-z (* 90 t)))))
"#,
                4.0 - sunk
            ),
        );
        let d = s.ink(["interfere", &rg]).ok();
        let rigid = pair(&d, "root/p/plate", "root/p/block");
        let passing = rows(&d).iter().find(|p| p["pair"] == json!(["root/post", "root/p/plate"]));
        let at = s.ink(["interfere", &rg, "--t", "0.7"]).ok();
        assert!(
            verdict(rigid) == want && !has(rigid, "motion") && passing.is_none_or(|p| verdict(p) == "holds") && verdict(pair(&at, "root/p/plate", "root/p/block")) == want,
            "sunk {sunk}:\n{}\n{}",
            pretty(&d),
            pretty(&at)
        );
    }
}

#[test]
fn a_keyed_shaft_turning_in_its_bore_fails_at_an_angle_of_its_turn_at_the_moment_it_reaches_it() {
    let s = Scratch::new("motion-keyed-shaft");
    let key = s.model(
        "keyed.ink",
        r#"(import "std")
(body bore (shape (subtract (box 30 30 6) (min (cylinder-z 5 20) (move (box 3 2 20) 5 0 0)))))
(body shaft (shape (min (cylinder-z 4.9 10) (move (box 2 1.6 10) 5 0 0))))
(t-span 0 1)
(assembly root (instance b bore) (instance s shaft (at) (turn z (* 360 t))))
"#,
    );
    let p = row(s.ink(["interfere", &key]));
    let (at, w) = (360.0 * num(&p["t"]), nums(&p["witness"]));
    assert!(verdict(&p) == "fails" && at > 5.0 && at < 355.0 && w.len() == 3 && (4.9..6.1).contains(&w[0].hypot(w[1])), "{}", pretty(&p));
}

/// tests/fixtures/geared-drive.ink, its gears unrelated (no gear-mesh: a
/// relation's zone would leave them out) and its wheel turned by `wheel`.
fn geared(wheel: Option<&str>, related: bool) -> String {
    let text = std::fs::read_to_string(fixture("geared-drive.ink")).unwrap();
    let text = if related { text } else { text.replace("  (relate p w (gear-mesh zp zw module 20 (* 2 backlash))))", "  )") };
    assert!(text.contains(GEAR_TURN), "geared-drive.ink no longer turns its wheel by {GEAR_TURN}");
    wheel.map_or(text.clone(), |w| text.replace(GEAR_TURN, w))
}

const GEAR_TURN: &str = "(turn z (mesh-turn zp zw 0 (motor-turn output-rpm)))";

#[test]
fn a_pinion_against_its_wheel_holds_measured_over_the_spin_up() {
    // The wheel turns by mesh-turn of its motor's drive, its pinion on the
    // motor's shaft (the gears unrelated here: a relation's zone would leave
    // them out). Over the spin-up (18000 degrees) their teeth hold at the
    // moments measured: an estimate (a proof over one 24-degree period of
    // the pair is a joint-values proof the one search does not make). A
    // pinion with a keyway, the same over its first 16 degrees.
    let s = Scratch::new("motion-geared-period");
    let unrelated = s.model("geared-unrelated.ink", &geared(None, false));
    let keyed =
        s.model("geared-keyed.ink", &geared(None, false).replace("(body pinion (shape (gear zp 8))", "(body pinion (shape (subtract (gear zp 8) (move (box 2 3 20) 4 0 0)))"));
    for (m, span) in [(unrelated, "0:1"), (keyed, "0:0.03")] {
        let d = s.ink(["interfere", &m, "--t-span", span, "--target", "root/fast"]).ok();
        let (p, rw) = (pair(&d, "root/fast/p", "root/fast/w"), pair(&d, "root/fast/m/rotor", "root/fast/w"));
        assert!(verdict(p).starts_with("holds") && verdict(rw).starts_with("holds"), "{m}:\n{}\n{}", pretty(p), pretty(rw));
    }
}

#[test]
fn a_wheel_of_the_wrong_ratio_or_phase_fails_against_its_pinion_at_an_angle_the_motor_reaches() {
    // A wrong ratio (20 teeth for 15: by 45 degrees its teeth overlap 0.45
    // deep; e2e had 17, 0.12 deep, at 26 CPU-s a call) and half a tooth out
    // of phase fail at an angle the motion reaches, over the motor's first
    // 45 degrees (t to 0.05). Related by gear-mesh as written (no zone),
    // the wrong ratio fails too: the teeth are checked.
    let s = Scratch::new("motion-geared");
    let ratio = GEAR_TURN.replace("zp zw", "20 zw");
    let phase = format!("(turn z (+ 4 {}))", &GEAR_TURN[8..GEAR_TURN.len() - 1]);
    for (name, text) in [("geared-ratio", geared(Some(&ratio), false)), ("geared-phase", geared(Some(&phase), false)), ("geared-related-ratio", geared(Some(&ratio), true))] {
        let p = s.model(&format!("{name}.ink"), &text);
        let d = s.ink(["interfere", &p, "--t-span", "0:0.05", "--target", "root/fast"]).ok();
        let (p0, rw) = (pair(&d, "root/fast/p", "root/fast/w"), pair(&d, "root/fast/m/rotor", "root/fast/w"));
        let t = num(&p0["t"]);
        let mut ok = verdict(p0) == "fails" && (0.0..=0.05).contains(&t) && verdict(rw).starts_with("holds");
        if name == "geared-related-ratio" {
            ok &= p0["relations"][0]["kind"] == "gear-mesh" && !has(p0, "contact_zones");
        }
        assert!(ok, "{name}:\n{}\n{}", pretty(p0), pretty(rw));
    }
}

#[test]
fn drives_are_read_in_real_numbers_not_the_steps_their_floats_take() {
    let s = Scratch::new("motion-real-drives");
    // (- (+ 2^53 q) 2^53) is q, so two cubes 1 mm apart sliding by q and by
    // it never meet.
    let path = s.model(
        "float-shared-small.ink",
        r#"(import "std")
(t-span 0 1)
(def q (* 4 t))
(body cube (shape (box 0.4 0.4 0.4)))
(assembly root
 (instance a cube (at) (slide x q))
 (instance b cube (at (move 1 0 0)) (slide x (- (+ 9007199254740992 q) 9007199254740992))))
"#,
    );
    let p = row(s.ink(["interfere", &path]));
    assert!(verdict(&p) == "holds", "{}", pretty(&p));
    // A cube slid by its difference from q (0 in reals) does not move, so a
    // turning rotor 0.6 from it is apart by their boxes.
    let path = s.model(
        "envelope-float.ink",
        r#"(import "std")
(t-span 0 1)
(def q (* 4 t))
(def d (- (- (+ 9007199254740992 q) 9007199254740992) q))
(body rotor (shape (cylinder-z 2 2)))
(body cube (shape (box 0.4 0.4 0.4)))
(assembly root
 (instance a rotor (at) (turn z q))
 (instance b cube (at (move 2.8 0 0)) (slide x d)))
"#,
    );
    let d = s.ink(["interfere", &path]).ok();
    assert!(rows(&d).is_empty() && d["pairs_apart_by_bounding_boxes"] == 1 && d["overlapping_pairs"] == 0, "{}", pretty(&d));
    // Each drive kept as written: B's ((0.1 t + 0.2 t) - 0.30000000000000004
    // t) 3.6e16 is -0.9992 t in reals, so B meets A near t = 1.0003; with its
    // coefficient computed (folded to 0) B would stay at 3 and the pair be
    // proven. Not guaranteed to hold: B reaches 8e10 mm, outside ink's range.
    let az = s.model(
        "affine-zero.ink",
        r#"(import "std")
(t-span 0.996 1.004)
(def a (* 2 t))
(def b (* (- (+ (* 0.1 t) (* 0.2 t)) (* 0.30000000000000004 t)) 3.6e16))
(body cube (shape (box 0.4 0.4 0.4)))
(assembly root
 (instance a cube (at) (slide x a))
 (instance b cube (at (move 3 0 0)) (slide x b)))
"#,
    );
    let d = s.ink(["interfere", &az]).err();
    assert!(d["error"]["message"].as_str().is_some_and(|m| m.contains("outside ink's range")), "{}", pretty(&d));
}

#[test]
fn a_pair_far_out_whose_world_placements_round_together_is_checked_in_the_world() {
    // The model composes the shared links, each composition rounded: two
    // 0.4 mm cubes at x = 0 and 1 in an assembly placed at x = 2^53 are at
    // one x in the world, where floats 1 mm apart leave it not held. The
    // span is a thousandth of the turn. 2^53 is outside ink's range: refused.
    let s = Scratch::new("motion-far");
    let far = s.model(
        "shared-frame-far.ink",
        r#"(import "std")
(t-span 0 0.001)
(body cube (shape (box 0.4 0.4 0.4)))
(assembly (pair (q 0))
 (instance a cube (at) (turn z q))
 (instance b cube (at (move 1 0 0)) (turn z (* 2 q))))
(assembly root (instance p (pair (* 360 t)) (at (move 9007199254740992 0 0))))
"#,
    );
    let d = s.ink(["interfere", &far]).err();
    assert!(d["error"]["message"].as_str().is_some_and(|m| m.contains("a coordinate 9")), "{}", pretty(&d));
}

#[test]
fn a_folds_seam_is_not_proven_away() {
    // One fold whose pitch is pi / 2 + 1e-10, a 1e-10 rad window at its
    // seam. The rotor against a post inside it fails (r 12 and r 100).
    let s = Scratch::new("motion-seam");
    for (name, r, off) in [("seam-12", 12.0, "1e-10"), ("seam-100", 100.0, "2.4e-10")] {
        let path = s.model(
            &format!("{name}.ink"),
            &format!(
                r#"(import "std")
(def p (+ (/ pi 2) {off}))
(def a (atan2 y x))
(def f (mod a p))
(def gate (max (- 5e-11 f) (- f 1.5e-10)))
(body rotor (shape (max (cylinder-z {r} 4) (select gate 1 -1))))
(body post (shape (sphere 0.2)))
(t-span 0 1)
(assembly root (instance r rotor (at) (turn z (* 720 t))) (instance p post (at (move {} 0 0))))
"#,
                -0.7 * r
            ),
        );
        let p = row(s.ink(["interfere", &path]));
        assert!(verdict(&p) == "fails", "{name}: {}", pretty(&p));
    }
    // Split spokes (a rib cut by a slit at the seam, deeper filled than the
    // accuracy, and a post thicker than twice it: both sized to interfere's
    // 0.1): not guaranteed to hold.
    let spokes = s.model(
        "split-spokes.ink",
        r#"(import "std")
(def p (+ (/ pi 2) 1e-10))
(def a (atan2 y x))
(def f (mod a p))
(def gate (max (- 5e-11 f) (- f 1.5e-10)))
(def spoke (- (min f (- p f)) 0.02))
(body rotor (shape (min (cylinder-z 2 4) (max (cylinder-z 12 4) spoke (select gate 1 -1)))))
(body post (shape (box 0.4 0.4 0.4)))
(t-span 0 1)
(assembly root (instance r rotor (at) (turn z (* 720 t))) (instance p post (at (move 7 0 0))))
"#,
    );
    let p = row(s.ink(["interfere", &spokes]));
    assert!(verdict(&p) != "holds", "{}", pretty(&p));
}

#[test]
fn geared_drives_rotor_holds_against_its_wheel() {
    // Each turns about its own axis: the rotor misses the wheel at the
    // moments measured, an estimate (a proof by their envelopes, the
    // cylinders each sweeps, is a joint-values proof the one search does
    // not make).
    let d = Ink::new(["interfere", &fixture("geared-drive.ink")]).ok();
    let rotors: Vec<&Value> = rows(&d)
        .iter()
        .filter(|r| {
            let (a, b) = (r["pair"][0].as_str().unwrap(), r["pair"][1].as_str().unwrap());
            a.ends_with("m/rotor") && b == a.replace("m/rotor", "w")
        })
        .collect();
    assert!(rotors.len() == 1 && rotors.iter().all(|r| verdict(r).starts_with("holds")), "{}", pretty(&d));
}

#[test]
fn a_value_meeting_its_bound_at_every_moment_is_not_failed_and_ends_in_few_steps() {
    // No step is decided, and no moment is a confirmed violation (its f64
    // value rounds above 1 at some t). Halving each undecided step to the
    // depth would take 2^(4 depth) pieces, about 2 hours; it stops at the
    // undecided budget.
    let s = Scratch::new("motion-pyth");
    let py = s.model("pyth.ink", "(assert one (value (+ (square (cos t)) (square (sin t))) <= 1))\n");
    let a = row(s.ink(["assert", &py, "--t-span", "0:1"]));
    assert!(matches!(verdict(&a).as_str(), "unknown" | "holds" | "holds estimate") && a["motion"]["steps"].as_u64().unwrap() < 1000, "{}", pretty(&a));
}

#[test]
fn a_lever_turned_by_t_to_short_of_its_stop_holds_over_the_span() {
    // An arm turned by t over 0..44 degrees, its limits 0..45 (the agent
    // tester's lever), a stop just past 45: holds with t ranged over the
    // span asked.
    let s = Scratch::new("limits-cover");
    let m = s.model(
        "lever.ink",
        "(import \"std\")\n(t-span 0 44)\n(def arm (subtract (min (move (box 100 8 6) 50 0 0) (cylinder-z 6 6)) (cylinder-z 2.7 10)))\n\
(body arm (shape arm))\n(body stop (shape (rotate-z (move (box 10 10 10) 90 9 0) 45)))\n\
(assembly root (instance stop stop) (instance arm arm (at) (turn z t) (limits 0 45 (accuracy 0.00001))))\n\
(assert clear (disjoint (shape-of root/arm) (shape-of root/stop)) (accuracy 0.01))\n",
    );
    let d = s.ink(["assert", &m, "--name", "clear"]).ok();
    let o = &d["asserts"][0]["motion"];
    assert!(verdict(&d["asserts"][0]) == "holds" && o["span"] == serde_json::json!([0.0, 44.0]) && o["how"] == "t ranged", "{}", pretty(o));
}

/// A pair over a motion is searched in the frame its first shape rests in
/// (motion::moving_frame), and a witness found there is reported in the
/// world only where it is one of the world's own fields at its moment,
/// by bounds rounded outward (motion::Frame::confirmed). Two balls 0.1
/// apart along x overlap, both carried 1e16 out: the frame's witness,
/// moved there in floats, rounds to x = 1e16, outside both at every moment
/// It is not reported
/// as a failure: the world's own search is asked, and leaves it undecided,
/// as it did before frames.
#[test]
fn a_witness_in_a_frame_far_out_is_no_failure_where_the_world_does_not_show_it() {
    let s = Scratch::new("motion-far-frame");
    let mo = s.model(
        "far.ink",
        r#"(import "std")
(t-span 0 1)
(body a (shape (move (sphere 0.25) 1 0 0)))
(body b (shape (move (sphere 0.25) 1.1 0 0)))
(assembly root
 (instance a a (at (move (+ 1e16 (* 0.001 t)) 0 0)))
 (instance b b (at (move (+ 1e16 (* 0.001 t)) 0 0))))
(assert apart (disjoint (shape-of root/a) (shape-of root/b)))
"#,
    );
    let a = row(s.ink(["assert", &mo]));
    assert!(verdict(&a) != "fails", "{}", pretty(&a));
    // A clearance's pair, confirmed at its moment with t still ranged
    // there (a span of one value, its placement's value in real numbers
    // bounded outward): pinned as a float, 1e16 + 0.001 t + 1 rounds, and
    // a pair 0.9 off both shapes would show inside both.
    let mo = s.model(
        "far-clearance.ink",
        r#"(import "std")
(t-span 0 1)
(body a (shape (sphere 0.25)))
(body b (shape (move (sphere 0.25) 0.1 0 0)))
(assembly root
 (instance a a (at (move (+ 1e16 (* 0.001 t) 1) 0 0)))
 (instance b b (at (move (+ 1e16 (* 0.001 t) 1) 0 0))))
(assert clear (clearance (shape-of root/a) (shape-of root/b) 0.5))
"#,
    );
    let a = row(s.ink(["assert", &mo]));
    assert!(verdict(&a) != "fails", "{}", pretty(&a));
}

#[test]
fn a_lever_turning_into_a_post_fails_at_its_moment() {
    // A lever turning into a post at about 45 of its 90 degrees: found at
    // a moment there.
    let s = Scratch::new("joint-clash");
    let m = s.model(
        "lever.ink",
        "(import \"std\")\n(t-span 0 1)\n(body lever (shape (move (box 60 12 5) 24 0 0)))\n(body post (shape (move (cylinder-z 3 10) 20 34.641 0)))\n\
         (assembly root (instance p post) (instance l lever (at) (turn z (* 90 t)) (limits 0 95)))\n(assert misses (disjoint (shape-of root/l) (shape-of root/p)))\n",
    );
    let d = ink(["assert", &m, "--name", "misses"]);
    let a = &d["asserts"][0];
    assert!(verdict(a) == "fails" && (0.4..0.6).contains(&num(&a["t"])), "{}", pretty(a));
}

#[test]
fn a_pin_turning_in_a_toleranced_bore_is_proven_over_the_motion() {
    // An 8 f9 pin turning in an 8.4 H9 bore: a tolerance on the still
    // bracket's bore is one more ranged param of the one search with t,
    // proven. Its rate along the tolerance is its diameter's, bounded on
    // the bore's axis though sqrt's slope there is not; unbounded, the
    // search would halve its band to its finest everywhere and give up.
    let s = Scratch::new("motion-toleranced-bore");
    let m = s.model(
        "bore.ink",
        r#"(import "std")
(param bore-d 8.4 (tol 0 0.036))
(param pin-d 7.969 (tol -0.018 0.018))
(t-span 0 1)
(body bracket (shape (subtract (box 6 20 20) (cylinder-x (/ bore-d 2) 30))))
(body pin (shape (cylinder-z (/ pin-d 2) 10)))
(assembly root (instance b bracket) (instance p pin (at (rot-y 90)) (turn z (* 180 t))))
(assert runs-free (clearance (shape-of root/p) (shape-of root/b) 0.15))
"#,
    );
    let a = row(s.ink(["assert", &m]));
    assert!(verdict(&a) == "holds" && a["motion"]["how"] == "t ranged" && has(&a["param_spans"], "bore-d"), "{}", pretty(&a));
}

#[test]
fn a_claim_measured_over_a_motion_says_what_would_decide_it() {
    // A contact over the motion is not decided with t ranged (a contact
    // found at one moment says nothing of the others): the estimate says
    // why in its hint, with the search's own answer (motion.t_ranged).
    let s = Scratch::new("motion-measured-hint");
    let m = s.model(
        "belt.ink",
        r#"(import "std")
(t-span 0 10)
(body belt (shape (box (+ 100 (sin t)) 6 1.4)))
(body plate (shape (move (box 50 50 4) 0 0 -2.7)))
(assert on (touches (shape-of belt) (shape-of plate)))
"#,
    );
    let a = row(s.ink(["assert", &m]));
    let hint = a["hint"].as_str().unwrap_or("");
    assert!(
        verdict(&a) == "holds estimate"
            && a["motion"]["how"] == "measured"
            && has(&a["motion"]["t_ranged"], "param_spans")
            && hint.contains("at one moment says nothing")
            && hint.contains("claim what can hold at every moment")
            && !hint.contains("ranged params"),
        "{}",
        pretty(&a)
    );
}

/// Two balls riding one moving instance are asked once only where that
/// instance's placement is proven a rotation over the span (also when the
/// pair is disjoint, and when the scale is near one): turned about an axis given as a z too small to normalize
/// exactly, its matrix scales with the angle (by 1.17 at t = 0.5), so the
/// balls' gap of 1 is 0.855 there, under the 0.95 asked, balls 0.08 apart
/// meet, and a scale of 1 + 1e-12 closes a gap of 1e12 by more than its
/// margin. Asked once at the view's moment, each would hold; none does.
#[test]
fn parts_riding_a_scaling_link_are_not_asked_once() {
    let s = Scratch::new("motion-riding-scale");
    let model = |a: &str, b: &str, axis: &str, claim: &str| {
        format!(
            "(import \"std\")\n(t-span 0 1)\n(body a (shape {a}))\n(body b (shape {b}))\n(assembly pair (instance a a) (instance b b))\n\
             (assembly root (instance pair pair (at (rot 0 0 {axis} (* 180 t)))))\n(assert gap {claim})\n"
        )
    };
    let models = [
        model("(move (sphere 1) -1.5 0 0)", "(move (sphere 1) 1.5 0 0)", "2.6e-162", "(clearance (shape-of root/pair/a) (shape-of root/pair/b) .95)"),
        model("(move (sphere 1) -.96 0 0)", "(move (sphere 1) .96 0 0)", "1.6e-162", "(disjoint (shape-of root/pair/a) (shape-of root/pair/b))"),
        model("(sphere 1)", "(move (sphere 1) 1000000000004 0 0)", "1.0000150000000001e-156", "(clearance (shape-of root/pair/a) (shape-of root/pair/b) 1000000000001.5)"),
    ];
    for (k, text) in models.iter().enumerate() {
        let a = row(s.ink(["assert", &s.model(&format!("riding-{k}.ink"), text)]));
        assert!(verdict(&a) != "holds", "{text}\n{}", pretty(&a));
    }
}

/// A pair over a motion is asked in the first body's frame, its placement
/// inverted by its transpose, only where that placement is proven a
/// rotation over the span. Placed by `axes` with x 2.6e-162 (normalized to
/// a 1.17 scale in x and y), the first ball's transpose puts the second
/// where it is not, which would prove the pair apart; at t = 0.5 they
/// overlap 0.08 deep. It is not proven apart.
#[test]
fn a_pair_over_its_joint_is_framed_only_by_a_rotation() {
    let s = Scratch::new("motion-joint-scaled");
    let mo = s.model(
        "scaled.ink",
        r#"(import "std")
(t-span 0 1)
(body a (shape (move (sphere 1) 1.95 0 0)))
(body b (shape (sphere 1)))
(assembly root
 (instance a a (at (axes (origin 0 0 0) (x 2.6e-162 0 0) (z 0 0 1)) (move (* 10 t) 0 0)))
 (instance b b (at (move (* 10.01 t) 0 0))))
(assert apart (disjoint (shape-of root/a) (shape-of root/b)))
"#,
    );
    let a = row(s.ink(["assert", &mo]));
    assert!(verdict(&a) != "holds", "{}", pretty(&a));
}

/// An item's placement chain runs from the world down, the target's prefix
/// included (Doc::items), so `--target` cannot change an answer: a pair
/// under an instance scaled by `axes`, its prefix left out of its frame,
/// would be proven clear asked of `root/pair` and fail asked of `root`.
/// Both answer alike: 1.6e-162 is outside ink's range, and both refuse.
#[test]
fn a_target_path_does_not_change_an_interference_answer() {
    let s = Scratch::new("motion-target-prefix");
    let mo = s.model(
        "prefix.ink",
        r#"(import "std")
(t-span 0 1)
(body a (shape (sphere 1)))
(body b (shape (move (sphere 1) 1.82 0 0)))
(assembly pair
 (instance a a (at (move (* 10 t) 0 0)))
 (instance b b (at (move (* 10.01 t) 0 0))))
(assembly root
 (instance pair pair (at (axes (origin 0 0 0) (x 1.6e-162 0 0) (z 0 0 1)))))
"#,
    );
    let refused = |target: &str| s.ink(["interfere", &mo, "--target", target]).err()["error"]["message"].clone();
    let (inner, outer) = (refused("root/pair"), refused("root"));
    assert!(inner == outer && inner.as_str().is_some_and(|m| m.contains("1.6e-162")), "{inner} / {outer}");
}

/// A shape's parts a turn leaves unchanged are placed by the move alone
/// only in that shape's own field, by identity, never another field placed
/// alike: the first shape is a ball cut by 1 - 1000 (its axial field less another
/// part's riding the partner's placement), so its solid is the ball only
/// where the two axial fields meet, and that part's z pulses by 0.005
/// sin^2 (the pulse within the drift a pivot charges). Rewriting that part
/// with the partner would freeze its pulse, empty the first shape and hold
/// the clearance (and a disjoint, a subset of the outside) at 0.05, where
/// at t = 1/1920 the ball and the partner meet. None of the three holds.
#[test]
fn a_turn_shared_by_part_of_a_shape_moves_only_that_part() {
    let s = Scratch::new("motion-pivot-shared");
    let model = |claim: &str| {
        format!(
            r#"(import "std")
(t-span 0 1)
(def pulse (* 0.005 (square (sin (deg->rad (* 172800 t))))))
(body anchor (shape (min (tag ball (sphere 2)) (tag axial (+ z 7)))))
(body pin (shape (sphere 1.5)))
(body axial (shape (+ z 7)))
(assembly root
 (instance a anchor (at (move (* 10 t) 0 0)))
 (instance b pin (at (rot-z (* 5 t)) (move (* 10 t) 0 pulse)))
 (instance c axial (at (rot-z (* 5 t)) (move (* 10 t) 0 pulse))))
(def region (max (shape-of root/a ball) (- 1 (* 1000 (- (shape-of root/a axial) (shape-of root/c))))))
(assert hit {claim})
"#
        )
    };
    for (k, claim) in ["(clearance region (shape-of root/b) 0.05)", "(disjoint region (shape-of root/b))", "(subset region (neg (shape-of root/b)))"].iter().enumerate() {
        let a = row(s.ink(["assert", &s.model(&format!("shared-{k}.ink"), &model(claim))]));
        assert!(verdict(&a) != "holds", "{claim}: {}", pretty(&a));
    }
}

/// A turn about a round part's axis is dropped only where the turn is proven
/// a proper rotation about z (motion::split_about's gate): a partner
/// placed by `axes` with x 1.8e-162 (its square subnormal: normalized to
/// 0.81, a planar scale) has its radius-1 cylinder read as radius 1.23;
/// taken as a turn, its field would be replaced by the unscaled one,
/// radius 1, and a probe at 1.15 that meets it proven clear.
#[test]
fn a_turn_about_a_round_part_is_dropped_only_where_proven_a_rotation() {
    let s = Scratch::new("motion-pivot-scaled");
    let mo = s.model(
        "scaled.ink",
        r#"(import "std")
(t-span 0 1)
(def pulse (square (sin (deg->rad (* 172800 t)))))
(body probe (shape (move (sphere 0.1) 1.15 0 0)))
(body pin (shape (cylinder-z 1 2)))
(assembly root
 (instance a probe (at (move (* 10 t) 0 0)))
 (instance b pin (at (axes (origin 0 0 0) (x 1.8e-162 0 0) (z 0 0 1)) (rot-z (* 5 t)) (move (* 10 t) 0 0))))
(def active (max (shape-of root/a) (- 0.5 pulse)))
(assert hit (clearance active (shape-of root/b) 0.01))
"#,
    );
    let a = row(s.ink(["assert", &mo]));
    assert!(verdict(&a) != "holds", "{}", pretty(&a));
}

/// A pair moving together, the first placed by a matrix a hair from a
/// rotation (a scale of 1 + 1e-9 or 1 + 1e-6, a non-uniform scale, a
/// shear) and the second a long way out in its own coordinates: its
/// defect, times the second's lever, moves the second by far more than the
/// first's own box shows. The balls overlap at every moment; none of
/// these is proven apart.
#[test]
fn a_pair_placed_a_hair_from_a_rotation_is_not_proven_apart_by_its_partners_lever() {
    let s = Scratch::new("motion-partner-lever");
    let model = |a_axes: &str, b_local: &str, b_at: &str| {
        format!(
            "(import \"std\")\n(t-span 0 1)\n(body a (shape (sphere 1)))\n(body b (shape (move (sphere 1) {b_local})))\n(assembly root\n (instance a a (at (axes (origin 0 0 0) {a_axes}) (move (* 100000 t) 0 0)))\n (instance b b (at (move {b_at}))))\n(assert apart (disjoint (shape-of root/a) (shape-of root/b)))\n"
        )
    };
    let models = [
        model("(x 2.2227587517078364e-162 0 0) (z 0 0 1)", "-999999998.5 0 0", "(+ 1000000000 (* 100000.01 t)) 0 0"),
        model("(x 2.2227609722438268e-162 0 0) (z 0 0 1)", "-999998.5 0 0", "(+ 1000000 (* 100000.01 t)) 0 0"),
        model("(x 2.2227609722438268e-162 0 0) (z 0 0 2.2227598608644524e-162)", "-999998.5 0 0", "(+ 1000000 (* 100000.01 t)) 0 0"),
        model("(x 1 0 1) (z 0 0 2.2227609722438268e-162)", "0 0 -999998.5", "(* 100000.01 t) 0 1000000"),
    ];
    for (k, text) in models.iter().enumerate() {
        let a = row(s.ink(["assert", &s.model(&format!("lever-{k}.ink"), text)]));
        assert!(verdict(&a) != "holds", "{text}\n{}", pretty(&a));
    }
}

/// A value reading t is compared at a moment by its bounds rounded outward
/// ((1e16 + t) - 1e16 is 0 in floats, t in real numbers). One on its bound, (1e16 + t) - (1e16 + max(t, 0)), 0 in
/// both, is -2 to 2 bounded outward: its float value holding stands as an
/// estimate (an effort clamped to its limit is such a
/// value; outward bounds alone left it unknown), never a guarantee.
#[test]
fn a_value_on_its_bound_at_a_moment_holds_as_an_estimate() {
    let s = Scratch::new("motion-on-bound");
    let m = s.model("on-bound.ink", "(t-span 0.5 1)\n(assert on-bound (value (- (+ 1e16 t) (+ 1e16 (max t 0))) >= 0))\n");
    let a = row(s.ink(["assert", &m]));
    assert!(verdict(&a) == "holds estimate", "{}", pretty(&a));
}

/// Swept boxes move a body's box forward, R q + T, only where its placement
/// is an isometry over the span (rigid::Isometry); the field reads
/// R^T (p - T), whose solid R^-T moves. A sphere placed by `axes` with x 1.8e-162 (a planar
/// scale of 0.81) sits at x 123.5 and meets b's there; moved forward by R
/// its box lay near x 81, and interfere dropped the pair as apart by boxes.
#[test]
fn a_placement_no_isometry_is_not_swept_forward() {
    let s = Scratch::new("motion-swept-nonrigid");
    let m = s.model(
        "swept.ink",
        r#"(import "std")
(t-span 0 1)
(body a (shape (move (sphere 1) 100 0 0)))
(body b (shape (move (sphere 1) 123.5 0 0)))
(assembly root
 (instance a a (at (axes (origin 0 0 0) (x 1.8e-162 0 0) (z 0 0 1)) (move (* 0.1 t) 0 0)))
 (instance b b))
(assert apart (disjoint (shape-of root/a) (shape-of root/b)))
"#,
    );
    // 1.8e-162 is outside ink's range: interfere refuses, the claim is unknown.
    let d = s.ink(["interfere", &m, "--target", "root"]).err();
    assert!(d["error"]["message"].as_str().is_some_and(|e| e.contains("1.8e-162")), "{}", pretty(&d));
    assert!(verdict(&row(s.ink(["assert", &m]))) == "unknown");
}

/// Parts ride one moving instance, asked once, only where every placement
/// down to it is an isometry over the span: a rotation and a translation
/// finite there. Moved by ln((t - 0.5)^2), infinite at
/// t = 0.5 where the solids are empty, a pair welded or touching at every
/// other moment is not guaranteed so over the whole span; by ln(1 - t),
/// infinite at its end, the same; nor as two root instances moved alike.
#[test]
fn parts_riding_a_translation_infinite_at_a_moment_are_not_asked_once() {
    let s = Scratch::new("motion-nonfinite-riding");
    let pair = |claim: &str, b: &str, at: &str| {
        format!(
            "(import \"std\")\n(t-span 0 1)\n(body ball (shape (sphere 1)))\n(assembly pair (instance a ball) (instance b ball (at (move {b} 0 0))))\n(assembly root (instance pair pair (at (move {at} 0 0))))\n(assert it ({claim} (shape-of root/pair/a) (shape-of root/pair/b)))\n"
        )
    };
    let roots = "(import \"std\")\n(t-span 0 1)\n(body ball (shape (sphere 1)))\n(assembly root (instance a ball (at (move (ln (square (- t 0.5))) 0 0))) (instance b ball (at (move (+ 1 (ln (square (- t 0.5)))) 0 0))))\n(assert it (overlaps (shape-of root/a) (shape-of root/b)))\n";
    for (k, text) in [pair("overlaps", "1", "(ln (square (- t 0.5)))"), pair("touches", "2", "(ln (square (- t 0.5)))"), pair("touches", "2", "(ln (- 1 t))"), roots.to_string()]
        .iter()
        .enumerate()
    {
        let a = row(s.ink(["assert", &s.model(&format!("nonfinite-{k}.ink"), text)]));
        assert!(verdict(&a) != "holds" && a["motion"]["how"] != "asked once", "{text}\n{}", pretty(&a));
    }
}

#[test]
fn an_imported_claim_is_checked_over_the_span_of_the_file_it_is_written_in() {
    // A ball sliding 10 mm in a second leaves a 4 mm box after 0.1 s: over
    // its file's span, imported or not, it fails at a moment past 0; its
    // file's static claim stays static, checked once.
    let s = Scratch::new("motion-import-span");
    let ball = "(import \"std\")\n(t-span 0 1)\n(body ball (shape (sphere 1)))\n(assembly root (instance b ball (at) (slide x (* 10 t))))\n\
        (assert inside (subset (shape-of root/b) (box 4 4 4)))\n(assert still (subset (sphere 1) (box 4 4 4)))\n\
        (assert in-run (subset (shape-of root/b) (move (box 30 4 4) 5 0 0)))\n";
    s.model("ball.ink", ball);
    // A cube sliding 10 mm/s in y, over (0, 0.1): it stays in its box.
    s.model(
        "cube.ink",
        "(import \"std\")\n(t-span 0 0.1)\n(body cube (shape (box 1 1 1)))\n(assembly cubes (instance c cube (at) (slide y (* 10 t))))\n\
        (assert cube-in (subset (shape-of cubes/c) (box 4 4 4)))\n",
    );
    // A file of no span: its claim reads t, and takes the run file's span.
    s.model(
        "loose.ink",
        "(import \"std\")\n(body pin (shape (sphere 1)))\n(assembly pins (instance p pin (at) (slide z (* 10 t))))\n\
        (assert pin-in (subset (shape-of pins/p) (box 4 4 4)))\n",
    );
    let row_of = |d: &Value, n: &str| rows(d).iter().find(|r| r["name"] == n).cloned().unwrap_or_else(|| panic!("no {n}:\n{}", pretty(d)));
    let fails_past_0 = |r: &Value, span: [f64; 2]| verdict(r) == "fails" && num(&r["t"]) > 0.1 && nums(&r["motion"]["span"]) == span;
    for (name, text) in [("plain.ink", "(import \"ball.ink\")"), ("ns.ink", "(import \"ball.ink\" as m)")] {
        let w = s.model(name, &format!("(import \"std\")\n{text}\n"));
        let d = s.ink(["assert", &w]).ok();
        let p = if name == "ns.ink" { "m." } else { "" };
        let (inside, still) = (row_of(&d, &format!("{p}inside")), row_of(&d, &format!("{p}still")));
        assert!(fails_past_0(&inside, [0.0, 1.0]) && d["all_hold"] == false, "{name}: {}", pretty(&d));
        assert!(verdict(&still) == "holds" && !has(&still, "motion") && !note(&still).contains("t = 0 only"), "{name}: {}", pretty(&still));
        // Proven over the import's span, t ranged, t among its params.
        let run = row_of(&d, &format!("{p}in-run"));
        assert!(verdict(&run) == "holds" && nums(&run["motion"]["span"]) == [0.0, 1.0] && run["motion"]["how"] == "t ranged", "{name}: {}", pretty(&run));
        // ink check says which spans its imports' claims are checked over.
        let c = s.ink(["check", &w]).ok();
        let spans: Vec<(&String, Vec<f64>)> = c["import_t_spans"].as_object().into_iter().flatten().map(|(f, s)| (f, nums(s))).collect();
        assert!(c["t_span"].is_null() && spans.len() == 1 && spans[0].0.ends_with("ball.ink") && spans[0].1 == [0.0, 1.0], "{}", pretty(&c));
    }
    // Two imports of different spans: each claim over its own; the
    // importer's own span covers its own claims and those of a file of none.
    let w = s.model(
        "both.ink",
        "(import \"std\")\n(import \"ball.ink\")\n(import \"cube.ink\")\n(import \"loose.ink\")\n(t-span 0 0.05)\n\
        (assert mine (subset (shape-of root/b) (box 4 4 4)))\n",
    );
    let d = s.ink(["assert", &w]).ok();
    assert!(fails_past_0(&row_of(&d, "inside"), [0.0, 1.0]), "{}", pretty(&d));
    let cube = row_of(&d, "cube-in");
    assert!(verdict(&cube) == "holds" && nums(&cube["motion"]["span"]) == [0.0, 0.1], "{}", pretty(&cube));
    for n in ["pin-in", "mine"] {
        let r = row_of(&d, n);
        assert!(verdict(&r) == "holds" && nums(&r["motion"]["span"]) == [0.0, 0.05], "{n}: {}", pretty(&r));
    }
    // --t-span is the run's question: every claim over it.
    let d = s.ink(["assert", &w, "--t-span", "0:0.05"]).ok();
    assert!(rows(&d).iter().all(|r| verdict(r) == "holds"), "{}", pretty(&d));
    // Of no span anywhere, a claim reading t is asked at one moment: its
    // row records the moment, and the answer says so once, at the top,
    // naming it; check warns of what the file asks so, its imports' too.
    let w = s.model("still.ink", "(import \"std\")\n(import \"loose.ink\")\n");
    let d = s.ink(["assert", &w]).ok();
    let pin = row_of(&d, "pin-in");
    let top = d["note"].as_str().unwrap_or("");
    assert!(verdict(&pin) == "holds" && !has(&pin, "motion") && pin["one_moment"] == 0.0 && top.contains("claim pin-in at t = 0"), "{}", pretty(&d));
    let warned_by = |file: &str, what: &str, load: bool| {
        let c = s.ink(["check", file].into_iter().chain(load.then_some("--load"))).ok();
        c["warnings"].as_array().is_some_and(|w| w.iter().any(|w| w.as_str().is_some_and(|w| w.starts_with("the file reads t but states no t-span") && w.contains(what))))
    };
    let warned = |file: &str, what: &str| warned_by(file, what, false);
    assert!(warned(&w, "claim pin-in") && warned(&s.path("loose.ink"), "claim pin-in"), "an importer asks its import's claim at one moment too");
    // interfere: a pair apart by its boxes at that moment is named by its
    // instance that moves.
    let pinned = s.model(
        "pinned.ink",
        "(import \"std\")\n(body pin (shape (sphere 1)))\n(body post (shape (move (sphere 1) 10 0 0)))\n\
        (assembly pins (instance p pin (at) (slide z (* 10 t))) (instance q post))\n",
    );
    let i = s.ink(["interfere", &pinned, "--target", "pins"]).ok();
    assert!(i["note"].as_str().is_some_and(|n| n.contains("instance pins/p at t = 0")), "{}", pretty(&i));
    assert!(warned(&pinned, "instance pins/p"), "an assembly interfere would ask at one moment, as it names it");
    // A claim asked at one moment and one over its import's span: each
    // says its own moment, whatever their order (the note named the last
    // claim's, t = 2, for the first).
    s.model("spanned.ink", "(t-span 2 3)\n(assert covered (value t >= 0))\n");
    let mixed = s.model("mixed.ink", "(assert first (value t < 1))\n(import \"spanned.ink\")\n");
    let d = s.ink(["assert", &mixed]).ok();
    let top = d["note"].as_str().unwrap_or("");
    assert!(
        row_of(&d, "first")["one_moment"] == 0.0 && !has(&row_of(&d, "covered"), "one_moment") && top.contains("claim first at t = 0") && !top.contains("covered"),
        "{}",
        pretty(&d)
    );
    // A solve on t read only by a claim's bound.
    let bound = s.model("bound.ink", "(solve time-bound (unknowns (q 0)) (== q t))\n(assert bound-only (value 0 <= q))\n");
    let d = s.ink(["assert", &bound]).ok();
    assert!(row_of(&d, "bound-only")["one_moment"] == 0.0 && d["note"].as_str().is_some_and(|n| n.contains("claim bound-only")), "{}", pretty(&d));
    assert!(warned(&bound, "claim bound-only"), "a bound reading t");
    // Still parts, a contact zone that moves: the pair's question reads t.
    let zone = s.model(
        "zone.ink",
        "(import \"std\")\n(body ball (shape (sphere 1)))\n(assembly root (instance a ball) (instance b ball (at (move 1 0 0)))\n\
        (contact a b (move (sphere 10) (* 100 t) 0 0)))\n",
    );
    let i = s.ink(["interfere", &zone]).ok();
    assert!(rows(&i)[0]["one_moment"] == 0.0 && i["note"].as_str().is_some_and(|n| n.contains("pair root/a & root/b at t = 0")), "{}", pretty(&i));
    assert!(warned(&zone, "pair root/a & root/b"), "a zone reading t");
    // check's warning is what interfere asks: a moving zone of a pair its
    // boxes keep apart is asked nothing (no warning of the assembly), and
    // the default target of bodies is asked.
    let apart = s.model(
        "unused-zone.ink",
        "(import \"std\")\n(body ball (shape (sphere 1)))\n(assembly root (instance a ball) (instance b ball (at (move 10 0 0)))\n\
        (contact a b (move (sphere 10) (* 100 t) 0 0)))\n",
    );
    let i = s.ink(["interfere", &apart]).ok();
    assert!(i["pairs_apart_by_bounding_boxes"] == 1 && i.get("note").is_none() && !warned(&apart, ""), "{}", pretty(&i));
    for (name, x, what) in [("bodies.ink", 10, "instance a"), ("bodies-near.ink", 1, "pair a & b")] {
        let f = s.model(name, &format!("(import \"std\")\n(body a (shape (move (sphere 1) (* 10 t) 0 0)))\n(body b (shape (move (sphere 1) {x} 0 0)))\n"));
        let i = s.ink(["interfere", &f]).ok();
        assert!(i["note"].as_str().is_some_and(|n| n.contains(&format!("{what} at t = 0"))) && warned(&f, what), "{name}: {}", pretty(&i));
    }
    let wrap = s.model("bodies-wrap.ink", "(import \"bodies-near.ink\" as m)\n");
    assert!(warned(&wrap, "pair m.a & m.b"), "an importer's default target");
    // --load is loading only: its claims, which loading asks, but not
    // interfere's plan (the planner is not run).
    let near = s.path("bodies-near.ink");
    assert!(!warned_by(&near, "", true) && warned_by(&bound, "claim bound-only", true), "--load runs no planner");
    // A warning names a few, then how many more: an answer stays short.
    let many = s.model("many.ink", &(0..12).map(|k| format!("(assert c{k} (value t >= {k}))\n")).collect::<String>());
    let c = s.ink(["check", &many, "--load"]).ok();
    let w = c["warnings"].to_string();
    assert!(w.contains("claim c0, claim c1, claim c2 and 9 more") && !w.contains("claim c3"), "{w}");
    let sole = s.model("sole.ink", "(import \"std\")\n(body ball (shape (sphere 1)))\n(assembly root (instance b ball (at) (slide x (* 10 t))))\n");
    assert!(!warned(&sole, ""), "one instance: interfere refuses, assert has no claim");
    // A moment the run asks is its question: nothing to say.
    let d = s.ink(["assert", &bound, "--t", "0.5"]).ok();
    assert!(d.get("note").is_none() && !has(&row_of(&d, "bound-only"), "one_moment"), "a moment asked is the question: {}", pretty(&d));
    // interfere: the assembly over its file's span.
    s.model(
        "wall.ink",
        "(import \"std\")\n(t-span 0 1)\n(body ball (shape (sphere 1)))\n(body wall (shape (move (box 1 10 10) 5 0 0)))\n\
        (assembly root (instance w wall) (instance b ball (at) (slide x (* 10 t))))\n",
    );
    let d = s.ink(["interfere", &s.model("walled.ink", "(import \"std\")\n(import \"wall.ink\")\n")]).ok();
    assert!(d["overlapping_pairs"] == 1 && num(&rows(&d)[0]["t"]) > 0.0, "{}", pretty(&d));
}

/// Two balls closing over their file's (t-span 0 1), coinciding at t = 1,
/// imported into a file of no span: placed as is and inside a moved and
/// turned parent.
pub const SPANNED_PART: &str = "(import \"std\")\n(t-span 0 1)\n(body ball (shape (sphere 5)))\n\
    (assembly moving (instance a ball) (instance b ball (at (move (* 8 (- 1 t)) (* 8 (- 1 t)) 0))))\n";
pub const SPANNED_PATH: &str = "(import \"spanned-part.ink\" as m)\n\
    (assembly outer (instance inner m.moving (at (rot-z 30) (move 50 0 0))))\n\
    (assembly root (instance carried m.moving) (instance outer outer (at (rot-x 40) (move 0 100 0))))\n";
/// A ball moved 10 mm/s over (t-span 0 1), by a def, a frame and a
/// template's argument, and a wall at x 4 to 6 it reaches at t = 0.3.
const DRIVE: &str = "(import \"std\")\n(t-span 0 1)\n(body ball (shape (sphere 1)) (frame seat (at)))\n(def offset (* 10 t))\n\
    (def where (at (move offset 0 0)))\n(assembly driver (instance b ball where))\n";
const DRIVE_ARG: &str = "(import \"std\")\n(t-span 0 1)\n(body ball (shape (sphere 1)))\n(assembly (driver (d 0)) (instance b ball (at (move d 0 0))))\n";
const SOLVE_DRIVE: &str = "(import \"std\")\n(t-span 0 1)\n(solve position (unknowns (q 0)) (== q (* 10 t)))\n(body ball (shape (sphere 1)))\n\
    (assembly driver (instance b ball (at) (slide x q)))\n";
const SOLVE_BODY: &str = "(import \"std\")\n(t-span 0 1)\n(solve position (unknowns (q 0)) (== q (* 10 t)))\n(body ball (shape (move (sphere 1) q 0 0)))\n";
const SOLVE_OVERLAP: &str = "(import \"std\")\n(t-span 0 1)\n(solve position (unknowns (q 8)) (== q (* 8 (- 1 t))))\n(body mover (shape (move (sphere 5) q q 0)))\n";
const STATIC_PAIR: &str = "(import \"std\")\n(body block (shape (box 10 10 10)))\n(assembly pair (instance a block) (instance b block (at (move 7 0 0))))\n";
const CONTACT_PARENT: &str = "(import \"static-pair.ink\" as child)\n(import \"std\")\n(t-span 0 1)\n\
    (assembly parent (instance carried child.pair (at (rot-z 23) (move 100 200 0)))\n\
    \x20 (contact carried/a carried/b (move (box 40 40 40) (* 100 t) 0 0)))\n";
const SHORT_DRIVE: &str = "(import \"std\")\n(t-span 0 0.1)\n(body ball (shape (sphere 1)))\n(assembly moving (instance b ball (at (move (* 10 t) 0 0))))\n";
const LOOSE_PAIR: &str = "(import \"short-drive.ink\" as d)\n(import \"std\")\n(body wall (shape (move (box 2 4 4) 5 0 0)))\n\
    (assembly pair (instance c d.moving) (instance w wall))\n";

/// Two blocks overlapping 3 mm, carried 100 mm/s, their contact zone
/// declared against the carriage's motion: still in the world, it leaves
/// the overlap from t = 0.05. In a root, and in a parent turned twice and
/// moved.
macro_rules! counter_zone {
    () => {
        "(import \"std\")\n(t-span 0 1)\n(body block (shape (box 10 10 10)))\n(assembly pair (instance a block) (instance b block (at (move 7 0 0))))\n"
    };
}
macro_rules! counter_root {
    () => {
        "(assembly root (instance c pair (at (move (* 100 t) 0 0))) (contact c/a c/b (move (box 40 40 40) (* -100 t) 0 0)))\n"
    };
}
const COUNTER_ZONE: &str = counter_zone!();
const COUNTER_ROOT: &str = counter_root!();
const COUNTER: &[(&str, &str)] = &[("counter-zone.ink", concat!(counter_zone!(), counter_root!()))];
const COUNTER_IN_PARENT: &[(&str, &str)] = &[(
    "counter-turned.ink",
    concat!(
        counter_zone!(),
        "(assembly mover (instance c pair (at (move (* 100 t) 0 0))) (contact c/a c/b (move (box 40 40 40) (* -100 t) 0 0)))\n",
        "(assembly root (instance outer mover (at (rot-z 90) (rot-x 90) (move 100 200 300))))\n"
    ),
)];

/// A pair whose question reads t: the
/// files it is written in (the last is run), its target, how many pairs
/// fail, each at a moment past 0 of the model's spans' hull (0, 1), and
/// whether guaranteed (else measured: a shape whose own field moves).
pub struct Repro {
    pub name: &'static str,
    pub files: &'static [(&'static str, &'static str)],
    pub target: Option<&'static str>,
    pub fails: u64,
    pub guaranteed: bool,
}

const fn repro(name: &'static str, files: &'static [(&'static str, &'static str)], target: Option<&'static str>, fails: u64, guaranteed: bool) -> Repro {
    Repro { name, files, target, fails, guaranteed }
}

const HUNT: &[(&str, &str)] = &[("spanned-part.ink", SPANNED_PART), ("spanned-path.ink", "(import \"spanned-part.ink\" as m)\n(assembly root (instance carried m.moving))\n")];
const NESTED: &[(&str, &str)] = &[
    ("spanned-part.ink", SPANNED_PART),
    (
        "spanned-nested.ink",
        "(import \"std\")\n(import \"spanned-part.ink\" as m)\n(assembly holder (instance inner m.moving))\n(assembly root (instance outer holder (at (rot-z 37) (move 100 200 0))))\n",
    ),
];
const BOTH: &[(&str, &str)] = &[("spanned-part.ink", SPANNED_PART), ("both.ink", SPANNED_PATH)];
const WALL: &[(&str, &str)] = &[
    ("spanned-part.ink", SPANNED_PART),
    (
        "wall.ink",
        "(import \"std\")\n(import \"spanned-part.ink\" as m)\n(body wall (shape (move (box 1 1 1) 6 2 0)))\n(assembly root (instance carried m.moving) (instance w wall))\n",
    ),
];
const CONTACT: &[(&str, &str)] = &[
    ("static-pair.ink", STATIC_PAIR),
    ("contact-parent.ink", CONTACT_PARENT),
    ("contact-run.ink", "(import \"contact-parent.ink\" as m)\n(assembly root (instance c m.parent (at (rot-x 15) (move 0 0 50))))\n"),
];
const SOLVED: &[(&str, &str)] = &[
    ("solve-drive.ink", SOLVE_DRIVE),
    (
        "solve-placement.ink",
        "(import \"solve-drive.ink\" as m)\n(import \"std\")\n(body wall (shape (move (box 2 4 4) 5 0 0)))\n(assembly root (instance carried m.driver) (instance w wall))\n",
    ),
];
const SOLVED_BODY: &[(&str, &str)] = &[
    ("solve-body-drive.ink", SOLVE_BODY),
    (
        "solve-body.ink",
        "(import \"solve-body-drive.ink\" as m)\n(import \"std\")\n(body wall (shape (move (box 2 4 4) 5 0 0)))\n(assembly root (instance b m.ball) (instance w wall))\n",
    ),
];
const SOLVED_OVERLAP: &[(&str, &str)] = &[
    ("solve-overlap-drive.ink", SOLVE_OVERLAP),
    (
        "solve-overlap.ink",
        "(import \"solve-overlap-drive.ink\" as m)\n(import \"std\")\n(body still (shape (sphere 5)))\n(assembly root (instance a still) (instance b m.mover))\n",
    ),
];
const FALLBACK: &[(&str, &str)] = &[
    ("short-drive.ink", SHORT_DRIVE),
    ("loose-pair.ink", LOOSE_PAIR),
    ("fallback-span.ink", "(import \"loose-pair.ink\" as l)\n(t-span 0.9 1)\n(assembly root (instance c l.pair))\n"),
];
const BY_DEF: &[(&str, &str)] = &[
    ("drive.ink", DRIVE),
    (
        "by-def.ink",
        concat!("(import \"drive.ink\" as m)\n(import \"std\")\n(body wall (shape (move (box 2 4 4) 5 0 0)))\n", "(assembly root (instance b m.ball m.where) (instance w wall))\n"),
    ),
];
const BY_FRAME: &[(&str, &str)] = &[
    ("drive.ink", DRIVE),
    (
        "by-frame.ink",
        concat!(
            "(import \"drive.ink\" as m)\n(import \"std\")\n(body wall (shape (move (box 2 4 4) 5 0 0)))\n",
            "(assembly root (instance b m.ball (frame-of m.driver/b)) (instance w wall))\n"
        ),
    ),
];
const BY_BODY: &[(&str, &str)] = &[
    ("drive.ink", DRIVE),
    (
        "by-body.ink",
        concat!(
            "(import \"drive.ink\" as m)\n(import \"std\")\n(body wall (shape (move (box 2 4 4) 5 0 0)))\n",
            "(body copy (shape (shape-of m.driver/b)))\n(assembly root (instance b copy) (instance w wall))\n"
        ),
    ),
];
const BY_ASSEMBLY: &[(&str, &str)] = &[
    ("drive.ink", DRIVE),
    (
        "by-assembly.ink",
        concat!(
            "(import \"drive.ink\" as m)\n(import \"std\")\n(body wall (shape (move (box 2 4 4) 5 0 0)))\n",
            "(body copy (shape (shape-of m.driver)))\n(assembly root (instance b copy) (instance w wall))\n"
        ),
    ),
];
const FOLLOWER: &[(&str, &str)] = &[
    ("drive.ink", DRIVE),
    (
        "follower.ink",
        concat!(
            "(import \"drive.ink\" as m)\n(import \"std\")\n(body wall (shape (move (box 2 4 4) 5 0 0)))\n",
            "(assembly root (instance carried m.driver) (instance follower m.ball (frame-of carried/b)) (instance w wall)\n  (contact carried/b follower (box 10 10 10)))\n"
        ),
    ),
];
const BY_ARG: &[(&str, &str)] = &[
    ("drive-arg.ink", DRIVE_ARG),
    (
        "by-arg.ink",
        "(import \"drive-arg.ink\" as m)\n(import \"std\")\n(body wall (shape (move (box 2 4 4) 5 0 0)))\n(assembly root (instance c (m.driver (* 10 t))) (instance w wall))\n",
    ),
];

pub const REPROS: &[Repro] = &[
    repro("placed import", HUNT, None, 1, true),
    repro("placed import by path", HUNT, Some("root/carried"), 1, true),
    repro("placed import by name", HUNT, Some("m.moving"), 1, true),
    repro("nested", NESTED, None, 1, true),
    repro("nested by path", NESTED, Some("root/outer/inner"), 1, true),
    repro("placed twice", BOTH, None, 2, true),
    repro("placed in a parent by path", BOTH, Some("root/outer/inner"), 1, true),
    repro("a wall apart at t = 0", WALL, None, 2, true),
    repro("a moving contact zone", CONTACT, None, 1, true),
    repro("a moving contact zone by name", CONTACT, Some("m.parent"), 1, true),
    repro("a moving contact zone by path", CONTACT, Some("root/c"), 1, true),
    repro("a moving contact zone below its parent", CONTACT, Some("root/c/carried"), 1, true),
    repro("a joint driven by a solve on t", SOLVED, None, 1, true),
    repro("a body moved by a solve on t", SOLVED_BODY, None, 1, false),
    repro("a solve on t, boxes meeting at t = 0", SOLVED_OVERLAP, None, 1, false),
    repro("a mover's span, in a file of another", FALLBACK, None, 1, true),
    repro("a moving def", BY_DEF, None, 1, true),
    repro("a moving frame", BY_FRAME, None, 1, true),
    repro("a moving body copied", BY_BODY, None, 1, true),
    repro("a moving assembly copied", BY_ASSEMBLY, None, 1, true),
    repro("a follower", FOLLOWER, None, 2, true),
    repro("a template's argument", BY_ARG, None, 1, true),
    repro("a zone still in the world, its pair riding a carriage", COUNTER, None, 1, true),
    repro("the same, by path", COUNTER, Some("root/c"), 1, true),
    repro("the same, in a turned parent", COUNTER_IN_PARENT, None, 1, true),
    repro("the same, the parent by path", COUNTER_IN_PARENT, Some("root/outer"), 1, true),
    repro("the same, the carriage by path", COUNTER_IN_PARENT, Some("root/outer/c"), 1, true),
];

impl Repro {
    /// Its files written into `s`; the path of the one run.
    pub fn write(&self, s: &Scratch) -> String {
        let mut run = String::new();
        for (n, t) in self.files {
            run = s.model(n, t);
        }
        run
    }
    pub fn args(&self, file: &str) -> Vec<String> {
        let mut a = vec!["interfere".to_string(), file.to_string()];
        if let Some(t) = self.target {
            a.extend(["--target".to_string(), t.to_string()]);
        }
        a
    }
    /// Whether answer `d` is the one expected, said where not.
    pub fn check(&self, d: &Value) -> Result<(), String> {
        let failing: Vec<&Value> = d["pairs"].as_array().into_iter().flatten().filter(|p| p["status"] == "fails").collect();
        let ok = d["overlapping_pairs"] == self.fails
            && failing.len() as u64 == self.fails
            && failing.iter().all(|p| num(&p["t"]) > 0.0 && nums(&p["motion"]["span"]) == [0.0, 1.0] && (p["accuracy"] == 0.05) == self.guaranteed);
        if ok { Ok(()) } else { Err(format!("{}: {}", self.name, pretty(d))) }
    }
}

#[test]
fn interfere_asks_a_pair_reading_t_over_every_span_the_model_states() {
    // Whatever reads t (a placement, a joint, a shape, a contact zone, a
    // solve on t; a def, a frame or a shape copied from a file that moves
    // it), and whatever target holds it: the pair over the hull of every
    // span the model states, here (0, 1), and it fails past t = 0.
    for (k, r) in REPROS.iter().enumerate() {
        let s = Scratch::new(&format!("motion-repro-{k}"));
        let d = s.ink(r.args(&r.write(&s))).ok();
        r.check(&d).unwrap_or_else(|e| panic!("{e}"));
    }
    // --t and --t-span are the run's question: at t = 0 the wall is apart,
    // and over the mover's own (0, 0.1) it is never reached.
    let s = Scratch::new("motion-repro-asked");
    let wall = REPROS[7].write(&s);
    let d = s.ink(["interfere", &wall, "--t", "0"]).ok();
    assert!(d["overlapping_pairs"] == 0 && d["pairs_apart_by_bounding_boxes"] == 2, "{}", pretty(&d));
    let fallback = REPROS[15].write(&s);
    let d = s.ink(["interfere", &fallback, "--t-span", "0:0.1"]).ok();
    assert!(d["overlapping_pairs"] == 0, "{}", pretty(&d));
    // ink's range is judged over the same span: a ball the import moves
    // 2,000 km by t = 1 is refused, though at t = 0 it is near the wall.
    s.model("far.ink", "(import \"std\")\n(t-span 0 1)\n(body ball (shape (sphere 1)))\n(assembly far (instance b ball (at (move (* 2e9 t) 0 0))))\n");
    let far = s
        .model("far-root.ink", "(import \"far.ink\" as m)\n(import \"std\")\n(body wall (shape (move (box 2 4 4) 5 0 0)))\n(assembly root (instance c m.far) (instance w wall))\n");
    let e = s.ink(["interfere", &far]).err();
    assert!(e["error"]["message"].as_str().unwrap_or("").contains("outside ink's range"), "{}", pretty(&e));
    // The counter-moving zone over an explicit span, and the claim of the
    // same question, which fails as interfere does.
    for (k, target) in [(22, None), (23, Some("root/c")), (26, Some("root/outer/c"))] {
        let file = REPROS[k].write(&s);
        let mut args = vec!["interfere", file.as_str(), "--t-span", "0:1"];
        args.extend(target.map(|t| ["--target", t]).into_iter().flatten());
        REPROS[k].check(&s.ink(&args).ok()).unwrap_or_else(|e| panic!("--t-span 0:1: {e}"));
    }
    let claim =
        s.model("counter-claim.ink", &format!("{COUNTER_ZONE}{COUNTER_ROOT}(assert uncovered (empty (max (shape-of root/c/a) (shape-of root/c/b) (neg (box 40 40 40)))))\n"));
    let c = row(s.ink(["assert", &claim]));
    assert!(verdict(&c) == "fails" && num(&c["t"]) > 0.0, "{}", pretty(&c));
    REPROS[22].check(&s.ink(["interfere", &claim]).ok()).unwrap_or_else(|e| panic!("beside its claim: {e}"));
    // A pair reading no t is static over the span: asked once, no motion.
    let d = s.ink(["interfere", &s.model("ancestor-contact.ink", crate::assemblies::ANCESTOR_CONTACT)]).ok();
    assert!(d["overlapping_pairs"] == 0 && rows(&d).iter().all(|p| !has(p, "motion")), "{}", pretty(&d));
    // No file states a span: at t = 0, as the pair's row records, and the
    // answer says so once, naming it; the row's own note says nothing of it.
    s.model("loose-part.ink", &SPANNED_PART.replace("(t-span 0 1)\n", ""));
    let d = s.ink(["interfere", &s.model("loose-path.ink", "(import \"loose-part.ink\" as m)\n(assembly root (instance carried m.moving))\n")]).ok();
    assert!(
        d["overlapping_pairs"] == 0
            && rows(&d).len() == 1
            && note(&rows(&d)[0]).is_empty()
            && rows(&d)[0]["one_moment"] == 0.0
            && note(&d).starts_with("answered at one moment only")
            && note(&d).contains("pair root/carried/")
            && note(&d).contains("at t = 0"),
        "{}",
        pretty(&d)
    );
    // Nor where a span covers them.
    assert!(rows(&s.ink(["interfere", &REPROS[0].write(&s)]).ok()).len() == 1 && note(&s.ink(["interfere", &REPROS[0].write(&s)]).ok()).is_empty());
}

/// Following a solve along t stops at the call's time limit: the moment
/// past where it got is refused, saying so, and the call answers in time.
/// tests/motion/solve-branch.ink's four-bar with its crank sped up
/// (10 t + t^2 / 1000 rad) and its rocker's pivot creeping out (350 +
/// t / 100), so t is read outside the crank's angle and no period of it
/// is followed in its stead: seconds to t = 300 unstopped.
#[test]
fn following_a_solve_stops_at_the_time_limit() {
    let s = Scratch::new("follow-limit");
    let f = s.model("chirp.ink", "(import \"std\")\n(t-span 0 0.6283185307)\n(body crank (shape (move (box 100 4 4) 50 0 0)) (frame pin (at (move 100 0 0))))\n(body coupler (shape (move (box 400 4 4) 200 0 0)) (frame pin (at (move 400 0 0))))\n(body rocker (shape (move (box 300 4 4) 150 0 0)) (frame pin (at (move 300 0 0))))\n(param th3 0.85) (param th4 1.52)\n(assembly (fourbar (p 0))\n  (instance crank crank (at) (turn z p))\n  (instance coupler coupler (at (rot-z (rad->deg th3)) (frame-of crank pin)))\n  (instance rocker rocker (at (rot-z (rad->deg th4)) (move (+ 350 (* 0.01 t)) 0 0))))\n(assembly root (instance f (fourbar (rad->deg (+ (* 10 t) (* 0.001 t t))))))\n(def gap (- (get (frame-of root/f/coupler pin) \"origin\") (get (frame-of root/f/rocker pin) \"origin\")))\n(solve loop (unknowns th3 th4) (== (nth gap 0) 0) (== (nth gap 1) 0))\n(assert open (value (sin th4) > 0.5))\n");
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs_f64();
    let t0 = std::time::Instant::now();
    // Searches stop 3 s before the call's limit, to answer whole.
    let d = s.ink(["check", &f, "--t", "300"]).env("INK_NO_CACHE", "1").env("INK_DEADLINE", (now + 3.5).to_string()).ok();
    let took = t0.elapsed().as_secs_f64();
    let w: Vec<&str> = d["warnings"].as_array().unwrap().iter().filter_map(|w| w.as_str()).filter(|w| w.contains("is not followed")).collect();
    assert!(w.len() == 1 && w[0].contains("the call's time limit passed") && took < 3.0, "{took:.1} s: {}", pretty(&d));
}

/// A motion repeats only over the period its real formula has: the sine
/// of (x + 1) t warped to x = p = 0.1 is that of (1 + 0.1) t in real
/// numbers, not of 1.1 t (0.1 and 1.1 are floats: at t = 1e16 they
/// differ by 0.094), so over 1e17 s their difference is not guaranteed
/// under 0.001 by a search of 2 pi / 1.1.
#[test]
fn a_motion_repeats_only_over_its_real_period() {
    let s = Scratch::new("false-period");
    let f = s.model(
        "false-period.ink",
        "(param p 0.1)\n(def a (warp (sin (* (+ x 1) t)) p y z))\n(def b (sin (* 1.1 t)))\n(t-span 0 100000000000000000)\n(assert same (value (abs (- a b)) < 0.001))\n",
    );
    let d = s.ink(["assert", &f]).env("INK_NO_CACHE", "1").any();
    let c = &rows(&d)[0];
    assert!(d["all_hold_guaranteed"] == false && !(verdict(c) == "holds" && has(c, "accuracy")), "{}", pretty(&d));
}

/// Whether a coefficient is a float or a real number is the evaluator's
/// call: 1 + (select c 0.1 q) with q a tolerance is real arithmetic
/// though c picks 0.1 (1 + 0.1, not 1.1's float), so over 1e17 s its sine
/// is not guaranteed to follow sin(1.1 t) within 0.001.
#[test]
fn a_select_of_a_tolerance_is_real_in_a_period() {
    let s = Scratch::new("select-period");
    let f = s.model(
        "selected.ink",
        "(param c -1)\n(param q 0 (tol -1 1))\n(def a (sin (* (+ 1 (select c 0.1 q)) t)))\n(def b (sin (* 1.1 t)))\n(t-span 0 100000000000000000)\n(assert same (value (abs (- a b)) < 0.001))\n",
    );
    let d = s.ink(["assert", &f]).env("INK_NO_CACHE", "1").any();
    let c = &rows(&d)[0];
    assert!(d["all_hold_guaranteed"] == false && !(verdict(c) == "holds" && has(c, "accuracy")), "{}", pretty(&d));
}
