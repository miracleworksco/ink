//! ink-dynamics' commands end to end: `ink effort` (at a moment, over the
//! span, `--modes`), its reactions, and the claims it
//! answers (effort-limits, reaction, zmp), as JSON agents read.
//! ink-dynamics' own tests (ink-dynamics/src) check its numbers on
//! models whose masses they set; these check what the command line gives:
//! fields, refusals and their hints, with masses ink measures.

use crate::common::*;
use serde_json::Value;
use std::f64::consts::PI;

fn n(v: &Value) -> f64 {
    v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"))
}

/// math.isclose: within `rel` of the larger.
fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs())
}

/// JSON equal, numbers compared as numbers (ink writes 5 as 5.0).
fn same(a: &Value, b: Value) -> bool {
    match (a, &b) {
        (Value::Number(x), Value::Number(y)) => x.as_f64() == y.as_f64(),
        (Value::Array(x), Value::Array(y)) => x.len() == y.len() && x.iter().zip(y).all(|(x, y)| same(x, y.clone())),
        (Value::Object(x), Value::Object(y)) => x.len() == y.len() && x.iter().all(|(k, v)| y.get(k).is_some_and(|w| same(v, w.clone()))),
        _ => *a == b,
    }
}

fn s(v: &Value) -> &str {
    v.as_str().unwrap_or("")
}

fn strs(v: &Value) -> Vec<&str> {
    v.as_array().map_or(vec![], |a| a.iter().map(s).collect())
}

/// The row of `rows` whose `key` is `name`.
fn by<'a>(rows: &'a Value, key: &str, name: &str) -> &'a Value {
    rows.as_array().and_then(|a| a.iter().find(|r| r[key] == name)).unwrap_or_else(|| panic!("no {key} {name} in {}", pretty(rows)))
}

/// A claim of an `ink require` answer, by name.
fn claim<'a>(d: &'a Value, name: &str) -> &'a Value {
    by(&d["requirements"], "name", name)
}

fn notes_have(c: &Value, what: &[&str]) -> bool {
    c["notes"].as_array().is_some_and(|ns| ns.iter().any(|x| what.iter().all(|w| s(x).contains(w))))
}

// A pendulum on an actuator (its effort against a hand formula at the mass
// ink measures), a gate on a plain hinge (unpowered), a planet turned by t
// itself (driven by nothing), a spring whose shape moves (deforming). The
// effort-limits claim is written at line 11 (the instance p).
const EFFORT: &str = concat!(
    "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 2s)\n",
    "(body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))\n",
    "(assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -150N*mm 150N*mm)))\n",
    "(body leaf (shape (move (box 60 4 80) 30 0 -40)) (attr density 2.7g/cm3) (attr spin-inertia 5000g*mm2))\n",
    "(assembly (flap (angle 0)) (instance leaf leaf (at) (turn y angle)))\n",
    "(body planet (shape (cylinder-z 10 5)) (attr density 7.85g/cm3))\n",
    "(body (spring (len 30)) (shape (move (box 5 5 len) 0 0 (/ len 2))) (attr density 7.85g/cm3))\n",
    "(assembly root\n  (instance p (pendulum (* 60 (sin (* 180deg/s t)))))\n  (instance gate (flap (* 20 (sin (* 180deg/s t)))) (at (move 0 200 0)))\n",
    "  (instance planet planet (at (move 50 0 0) (rot-z (* 72 (/ t 1s)))))\n  (instance s (spring (+ 30 (* 5 (sin (* 180deg/s t))))) (at (move -100 0 0))))\n",
);

/// The pendulum's angle (rad) and the torque its motion needs (N mm) at t,
/// from the bar's mass, centre and inertia as `ink mass` gives them.
fn pendulum(f: &str) -> (impl Fn(f64) -> f64, impl Fn(f64) -> f64, f64) {
    let bar = ink(["mass", f, "--target", "root/p/bar", "--keys", "mass,center_of_mass,inertia"])["targets"][0].clone();
    let m = n(&bar["mass"]) / 1000.0;
    let (cx, cz) = (n(&bar["center_of_mass"][0]), n(&bar["center_of_mass"][2]));
    let d2 = cx * cx + cz * cz;
    let io = n(&bar["inertia"][1][1]) + 1000.0 * m * d2; // g mm^2
    let theta = |t: f64| 60f64.to_radians() * (PI * t).sin();
    let tau = move |t: f64| io * 1e-6 * -(PI * PI) * theta(t) + m * 9.81 * d2.sqrt() * theta(t).sin();
    (theta, tau, io)
}

#[test]
fn an_actuators_effort_at_a_moment_is_its_bodys_inertia_and_weight_about_the_joint() {
    let w = Scratch::new("dyn-effort-moment");
    let f = w.model("effort.ink", EFFORT);
    let (theta, tau, io) = pendulum(&f);
    let d = w.ink(["effort", &f, "--t", "0.25"]).ok();
    let a = &d["actuators"][0];
    let e = n(&a["effort"]);
    assert!(
        a["path"] == "root/p/bar"
            && close(e, tau(0.25), 1e-3)
            && n(&d["t"]) == 0.25
            && s(&d["method"]).contains("at moment t")
            && d["target"] == "root"
            && same(&d["gravity"], serde_json::json!([0, 0, -9.81])),
        "want {}:\n{}",
        tau(0.25),
        pretty(&d)
    );
    // With its rate, limits, margin and reflected inertia, in the
    // language's units (N mm, deg/s, g mm2 a turn), with no units map.
    assert!(
        close(n(&a["rate"]), 60.0 * PI * (PI / 4.0).cos(), 1e-9)
            && n(&a["lo"]) == -150.0
            && n(&a["hi"]) == 150.0
            && (n(&a["margin"]) - (e + 150.0).min(150.0 - e)).abs() < 1e-6
            && close(n(&a["reflected_inertia"]), io, 1e-6)
            && a.get("units").is_none(),
        "{}",
        pretty(a)
    );
    // Each coordinate: its value, rate and acceleration as passed, and the
    // effort along it per unit of it (per degree).
    let p = by(&d["coordinates"], "name", "p.angle");
    let names: Vec<&str> = d["coordinates"].as_array().unwrap().iter().map(|c| s(&c["name"])).collect();
    assert!(
        close(n(&p["value"]), theta(0.25).to_degrees(), 1e-9)
            && close(n(&p["rate"]), n(&a["rate"]), 1e-9)
            && close(n(&p["acceleration"]), -(PI * PI * theta(0.25)).to_degrees(), 1e-9)
            && close(n(&p["effort"]), e * PI / 180.0, 1e-9)
            && names.len() == 2
            && names.contains(&"p.angle")
            && names.contains(&"gate.angle"),
        "{}",
        pretty(&d["coordinates"])
    );
}

#[test]
fn effort_names_the_motion_no_actuator_supplies_and_says_what_it_rests_on() {
    let w = Scratch::new("dyn-effort-unpowered");
    let f = w.model("effort.ink", EFFORT);
    let d = w.ink(["effort", &f, "--t", "0.25"]).ok();
    // A hinge with no actuator moved by its argument is unpowered, with
    // the effort it needs.
    let up = &d["unpowered"];
    assert!(
        up.as_array().is_some_and(|u| u.len() == 1)
            && up[0]["coordinate"] == "gate.angle"
            && n(&up[0]["effort"]) != 0.0
            && n(&up[0]["effort"]) == n(&by(&d["coordinates"], "name", "gate.angle")["effort"])
            && d.get("unpowered_note").is_none(),
        "{}",
        pretty(&d)
    );
    // A body placed by t itself is driven by nothing; one whose shape moves
    // is deforming, with its mass.
    let spring = 7.85e-3 * 25.0 * (30.0 + 5.0 * (PI * 0.25).sin());
    assert!(
        strs(&d["driven_by_nothing"]) == ["root/planet"]
            && d["deforming"].as_array().is_some_and(|x| x.len() == 1)
            && d["deforming"][0]["target"] == "root/s"
            && close(n(&d["deforming"][0]["mass"]), spring, 0.02),
        "want spring {spring}: {} {}",
        d["driven_by_nothing"],
        d["deforming"]
    );
    // An estimate, with the masses it rests on (its method and what it
    // leaves out are the reference's, not in every answer).
    let ms = &d["masses"];
    assert!(
        d["accuracy"].is_null()
            && ["basis", "left_out", "units", "reactions_note"].iter().all(|k| d.get(k).is_none())
            && n(&ms["bodies"]) == 3.0
            && n(&ms["published"]) == 0.0
            && strs(&ms["spin_inertia"]) == ["root/gate/leaf"]
            && ms["res"].as_array().is_some_and(|r| r.len() == 2)
            && same(&ms["unsettled"], serde_json::json!([])),
        "{}",
        pretty(&d)
    );
}

#[test]
fn a_band_whose_carried_says_its_material_moves_is_carried_not_deforming() {
    let w = Scratch::new("dyn-effort-carried");
    let belt = w.model(
        "effort-belt.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(def r 10)\n",
            "(def P (path-loop (vec (path-arc (vec 0 0 0) (vec 0 0 1) r (vec 0 -1 0) pi) (path-line (vec 0 r 0) (vec -100 r 0))\n",
            "  (path-arc (vec -100 0 0) (vec 0 0 1) r (vec 0 1 0) pi) (path-line (vec -100 (neg r) 0) (vec 0 (neg r) 0))) (vec 0 0 1)))\n",
            "(body (band (a 0)) (shape (path-sweep (max (- (abs x) 0.5) (- (abs y) 3)) P)) (attr density 1g/cm3)\n",
            "  (attr carried (carried (max (- (abs y) 0.5) (- (abs z) 3)) P (path-length P) (* r (/ (* pi a) 180)))))\n",
            "(body hub (shape (cylinder-z 10 5)) (attr density 1g/cm3))\n",
            "(assembly (loop (a 0)) (instance b (band a)) (instance p hub (at) (turn z a)))\n(assembly root (instance l (loop (* 360 (/ t 1s) (/ t 1s)))))\n"
        ),
    );
    let c = w.ink(["effort", &belt, "--t", "0.5"]).ok();
    let band = 1e-3 * 6.0 * (2.0 * PI * 10.0 + 200.0);
    assert!(
        c["carried"].as_array().is_some_and(|x| x.len() == 1)
            && c["carried"][0]["target"] == "root/l/b"
            && close(n(&c["carried"][0]["mass"]), band, 0.02)
            && same(&c["deforming"], serde_json::json!([])),
        "want band {band}: {}",
        pretty(&c)
    );
}

#[test]
fn effort_over_the_span_gives_each_actuators_peak_rms_and_least_margin() {
    let w = Scratch::new("dyn-effort-span");
    let f = w.model("effort.ink", EFFORT);
    let (_, tau, _) = pendulum(&f);
    let d = w.ink(["effort", &f]).ok();
    let (a, mo) = (&d["actuators"][0], &d["motion"]);
    let ts: Vec<f64> = (0..20000).map(|k| k as f64 / 20000.0 * 2.0).collect();
    let rms = (ts.iter().map(|&t| tau(t).powi(2)).sum::<f64>() / ts.len() as f64).sqrt();
    let peak = ts.iter().map(|&t| tau(t).abs()).fold(0.0, f64::max);
    // Steps chosen from the motion, none unchecked.
    let steps = n(&mo["steps"]);
    assert!(s(&d["method"]).contains("each step"), "{}", pretty(&d["method"]));
    assert!(
        same(&mo["span"], serde_json::json!([0, 2]))
            && 16.0 < steps
            && steps < 4096.0
            && s(&mo["basis"]).contains("1%")
            && same(&mo["unchecked"], serde_json::json!([]))
            && same(&mo["loose"], serde_json::json!([])),
        "{}",
        pretty(mo)
    );
    let pk = &a["peak"];
    assert!(
        close(n(&pk["effort"]).abs(), peak, 0.01)
            && (n(&pk["t"]) % 1.0 - 0.5).abs() < 0.05
            && pk.get("rate").is_some()
            && pk.get("reflected_inertia").is_some()
            && close(n(&a["rms"]), rms, 0.01)
            && close(n(&a["margin"]["margin"]), 150.0 - n(&pk["effort"]).abs(), 0.01)
            && n(&a["margin"]["hi"]) == 150.0,
        "want peak {peak}, rms {rms}: {}",
        pretty(a)
    );
    // Each coordinate's peak, with its value, rate and acceleration there.
    let c = &by(&d["coordinates"], "name", "p.angle")["peak"];
    assert!(
        close(n(&c["effort"]).abs(), peak * PI / 180.0, 0.01) && (n(&c["value"]).abs() - 60.0).abs() < 1.0 && ["rate", "acceleration", "t"].iter().all(|k| c.get(k).is_some()),
        "{}",
        pretty(c)
    );
    // And the peak of what is unpowered.
    assert!(
        d["unpowered"].as_array().is_some_and(|u| u.len() == 1) && d["unpowered"][0]["coordinate"] == "gate.angle" && d["unpowered"][0]["peak"].get("effort").is_some(),
        "{}",
        pretty(&d["unpowered"])
    );
}

#[test]
fn an_actuators_effort_limits_claim_fails_where_the_effort_passes_a_limit() {
    // Each actuator under root makes a claim: its effort within its limits,
    // measured at the steps effort takes; failing, with its moment and terms.
    let w = Scratch::new("dyn-effort-claim");
    let f = w.model("effort.ink", EFFORT);
    let d = w.ink(["require", &f]).ok();
    let a = &d["requirements"][0];
    let terms: Vec<&String> = a["terms"].as_object().map_or(vec![], |t| t.keys().collect());
    assert!(
        a["name"] == "effort-limits:root/p/bar"
            && verdict(a) == "fails estimate"
            && (n(&a["t"]) % 1.0 - 0.5).abs() < 0.1
            && terms.len() == 4
            && ["effort", "rate", "lo", "hi"].iter().all(|k| terms.iter().any(|t| t == k))
            && n(&a["terms"]["effort"]).abs() > 150.0
            && strs(&a["solved_inputs"]) == ["effort:root/p/bar"]
            && s(&a["motion"]["basis"]).contains("1%")
            && s(&a["from"]).ends_with("effort.ink:11"),
        "{}",
        pretty(a)
    );
}

#[test]
fn an_actuators_effort_says_which_moving_body_it_leaves_out() {
    // A body whose shape rides a coordinate is not counted: the actuator's
    // effort leaves it out and says so. (Its claim's being unknown is
    // ink-dynamics' a_claim_says_what_its_effort_leaves_out.)
    let w = Scratch::new("dyn-effort-leaves-out");
    let lift = w.model(
        "lift.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 1s)\n(body rotor (shape (cylinder-z 3 10)) (attr density 7.8g/cm3))\n",
            "(body (blk (h 0)) (shape (move (box 50 50 50) 100 0 (+ 25 h))) (attr density 7.85g/cm3))\n",
            "(assembly (lift (z 0)) (instance r rotor (at) (turn z (* 45 z)) (effort-limits -1N*mm 1N*mm)) (instance b (blk z)))\n",
            "(assembly root (instance l (lift (* 50 (/ t 1s) (/ t 1s)))))\n"
        ),
    );
    let a = w.ink(["effort", &lift, "--t", "0.5"]).ok()["actuators"][0].clone();
    let out: Vec<&str> = strs(&a["leaves_out"]).into_iter().map(|x| x.split(" (").next().unwrap()).collect();
    assert!(out == ["root/l/b"], "{}", pretty(&a));
}

#[test]
fn without_gravity_effort_claims_are_not_made_and_effort_gives_the_motion_alone() {
    let w = Scratch::new("dyn-effort-gravity");
    let still = w.model("still.ink", &EFFORT.replace("(gravity 0 0 -9.81m/s2)\n", ""));
    let d = w.ink(["check", &still, "--load"]).ok();
    assert!(
        strs(&d["not_checked"]) == ["effort-limits:root/p/bar"]
            && s(&d["not_checked_note"]).contains("(gravity")
            && !d["requirements"].as_array().unwrap().iter().any(|x| s(&x["name"]).starts_with("effort")),
        "check --load: {}",
        pretty(&d)
    );
    // The command gives what needs no gravity (the pendulum's fastest rate, 60 pi deg/s) and no effort, saying why.
    let d = w.ink(["effort", &still]).ok();
    let a = &d["actuators"][0];
    assert!(d["gravity"].is_null() && s(&a["lacks"]).contains("(gravity gx gy gz)") && a["peak"]["effort"].is_null() && a["rms"].is_null(), "{}", pretty(&d));
    assert!((num(&a["peak"]["rate"]).abs() - 60.0 * std::f64::consts::PI).abs() < 1e-6, "{}", pretty(a));
    assert!(d["reactions"].as_array().unwrap().iter().all(|r| r["lacks"].is_string() && r["peak_force"].is_null()), "{}", pretty(&d));
    assert!(d["coordinates"].as_array().unwrap().iter().all(|c| c["lacks"].is_string() && c["peak"].is_null()) && d["unpowered"].is_null(), "{}", pretty(&d));
    let twice = w.model("twice.ink", &EFFORT.replace("(t-span", "(gravity 0 0 -9.8m/s2)\n(t-span"));
    let e = w.ink(["check", &twice]).err()["error"].clone();
    assert!(s(&e["message"]).contains("twice"), "gravity twice: {}", pretty(&e));
}

/// A marker riding a simulated block, a written wobble on top: its speed is
/// the sum of the two, and its peak need be at neither one's, so the
/// simulation's turning points do not cover it. Read over the span, its
/// speed claim is not available, saying why (it held at the moments read,
/// while it passes 125 mm/s at t = 1/128 s); at that moment it fails.
#[test]
fn a_speed_of_simulated_and_written_motion_together_is_not_read_at_the_simulations_moments() {
    let w = Scratch::new("dyn-compound-speed");
    let f = w.model(
        "compound.ink",
        "(import \"std\")\n(t-span 0 1s)\n(gravity 0 0 0)\n(body block (shape (box 100 100 100)) (attr mass 1kg))\n(body marker (shape (box 1 1 1)) (attr mass 0g))\n\
         (assembly (cart (s 0)) (instance block block (at) (slide x s))\n  (instance marker marker (at) (slide x (+ s (* (/ -50 (* 64 pi)) (cos (* 64 180deg/s t))))) (speed-limits 0 125mm/s)))\n\
         (assembly root (instance c (cart (simulate 0 100))))\n",
    );
    let name = "speed-limits:root/c/marker";
    let over = claim(&w.ink(["require", &f, "--name", name]).ok(), name).clone();
    assert!(verdict(&over) == "unknown" && s(&over["note"]).contains("the drive's rate cannot be bounded"), "{}", pretty(&over));
    let at = claim(&w.ink(["require", &f, "--name", name, "--t", "0.0078125"]).ok(), name).clone();
    assert!(verdict(&at) == "fails estimate" && close(n(&at["value"]), -25.0, 1e-6), "{}", pretty(&at));
}

#[test]
fn a_simulated_motion_says_which_coordinates_are_simulated_and_how() {
    let w = Scratch::new("dyn-effort-simulated");
    let sim = w.model(
        "sim.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 1s)\n(body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))\n",
            "(assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle)))\n(assembly root (instance p (pendulum (simulate 60))))\n"
        ),
    );
    let d = w.ink(["effort", &sim, "--t", "0.5"]).ok();
    assert!(
        strs(&d["simulated"]["coordinates"]) == ["p.angle"] && d["simulated"].get("basis").is_none() && d.get("unpowered").is_none_or(|u| u == &serde_json::json!([])),
        "{}",
        pretty(&d)
    );
}

#[test]
fn a_simulations_springs_give_each_joints_following_error_and_ringing_mm_a_slide() {
    // A 1 kg slide on (stiffness 100) commanded from 5 to 0 mm rings 5 mm
    // either way at (1/2 pi) sqrt(k / m) = 50.33 Hz, in mm.
    let w = Scratch::new("dyn-effort-springs");
    let spring = w.model(
        "spring.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 0.1s)\n(body b (shape (box 100 100 100)) (attr density 1g/cm3))\n",
            "(assembly (axis (c 0)) (instance b b (at) (slide x c) (stiffness 100N/mm)))\n(assembly root (instance a (axis (select (- (/ t 1s) 0.01) 5 0))))\n",
            "(springs) (def lag (following-error root/a/b))\n"
        ),
    );
    let j = w.ink(["effort", &spring]).ok()["simulated"]["springs"]["joints"][0].clone();
    assert!(
        j["joint"] == "root/a/b" && (n(&j["peak"]).abs() - 5.0).abs() < 1e-5 && (n(&j["ringing"]["frequency"]) - 50.329).abs() < 0.01 && j.get("units").is_none(),
        "{}",
        pretty(&j)
    );
}

#[test]
fn an_unpowered_peak_over_the_span_carries_the_bound_of_what_it_reads_round() {
    // The servo arm's gate: a leaf on a plain hinge, read round over the span
    // (its weight off the hinge averaged away). Its peak with its bound
    // covers what each moment's own reading needs.
    let w = Scratch::new("dyn-effort-unpowered-bound");
    let f = root("examples/servo-arm-effort.ink");
    let gate = |d: &Value| by(&d["unpowered"], "coordinate", "gate.angle").clone();
    let span = gate(&w.ink(["effort", &f]).ok());
    let (peak, b) = (n(&span["peak"]["effort"]).abs(), span.get("bound").map_or(0.0, n));
    for t in ["0", "0.5", "1", "3"] {
        let at = n(&gate(&w.ink(["effort", &f, "--t", t]).ok())["effort"]).abs();
        assert!(at <= peak + b, "at t = {t} the gate needs {at}, past the span's peak {peak} and bound {b}: {}", pretty(&span));
    }
}

#[test]
fn an_integrated_state_checks_cancelled_speeds_in_its_command_history() {
    let w = Scratch::new("integrated-command");
    let arm = include_str!("../../../examples/simulated-arm.ink");
    let plain = w.model("arm.ink", arm);
    let d = w.ink(["require", &plain]).ok();
    for name in ["speed-limits:root/a/s/output", "effort-limits:root/a/s/output"] {
        let a = claim(&d, name);
        assert!(a["status"] == "holds" && a["accuracy"].is_null(), "{}", pretty(a));
    }
    // The whole command becomes a plain zero in dynamics. Integration
    // past the switch used an undefined speed. Even a later point query
    // of the position depends on that history.
    let commanded = "(- (* 40 (- 45 swing)) (* 1.2 (diff swing t)))";
    assert!(arm.contains(commanded));
    let bad = w.model("bad.ink", &arm.replace(commanded, "(* 0 (diff (select (- (/ t 1s) 0.53125) 0 10) t))"));
    for (t, result) in [("0.25", "holds"), ("0.75", "unknown")] {
        let d = w.ink(["require", &bad, "--t", t, "--name", "no-overshoot"]).any();
        let a = claim(&d, "no-overshoot");
        assert_eq!(a["status"], result, "{}", pretty(a));
        if result == "unknown" {
            assert!(s(&a["note"]).contains("may jump"), "{}", pretty(a));
        }
    }
}

#[test]
fn a_stopped_simulation_says_where_and_why_and_its_efforts_are_partial() {
    // A PD command lifting 1 kg on a lossy lead screw comes to rest held by
    // friction, which a simulation does not model: it stops there, and the
    // efforts over the rest of the span read nothing.
    let w = Scratch::new("dyn-effort-stopped");
    let f = w.model(
        "lift.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 1s)\n(body rotor (shape (cylinder-z 10 30)) (attr density 7.85g/cm3))\n",
            "(body block (shape (box 50 50 50)) (attr mass 1kg))\n(def a (simulate 0))\n",
            "(assembly (lift (a 0)) (instance rotor rotor (at (move 0 0 -40)) (turn z a) (effort-limits -100N*mm 100N*mm))\n",
            "  (instance c block (at (move 0 0 100)) (slide z (/ a 360))) (relate rotor c (lead-screw 1 0.9)))\n",
            "(assembly root (instance l (lift a)))\n(command root/l/rotor (+ 1.5613 (* 0.05 (- 360 a)) (* -0.0025 (diff a t))))\n"
        ),
    );
    let d = w.ink(["effort", &f]).ok();
    let (sim, mo) = (&d["simulated"], &d["motion"]);
    let stop = n(&sim["stopped"]["t"]);
    let unchecked = mo["unchecked"].as_array().cloned().unwrap_or_default();
    assert!(
        (stop - 0.3144).abs() < 1e-3
            && s(&sim["stopped"]["why"]).contains("friction")
            && n(&sim["end"]) == stop
            && mo["partial"] == true
            && unchecked.len() == 1
            && n(&unchecked[0][0]) <= stop + 1e-3
            && n(&unchecked[0][1]) == 1.0,
        "{}",
        pretty(&serde_json::json!({"simulated": sim, "motion": mo}))
    );
    // A simulation that runs its span says so too.
    let still = w.model("still.ink", &std::fs::read_to_string(&f).unwrap().replace("(lead-screw 1 0.9)", "(lead-screw 1 1)"));
    let d = w.ink(["effort", &still]).ok();
    assert!(d["simulated"]["stopped"].is_null() && n(&d["simulated"]["end"]) == 1.0 && d["motion"]["partial"] == false, "{}", pretty(&d["simulated"]));
}

#[test]
fn a_sprung_actuators_effort_is_its_springs_as_its_claim_reads_it() {
    // 1 kg on a slide its stiffness holds, commanded x = 100 t^2 mm from
    // rest: rigidly it needs m a = 0.2 N; released, its spring overshoots to
    // 2 m a = 0.4 N. ink effort gives the spring's (the claim's), the rigid
    // demand under `rigid`.
    let w = Scratch::new("dyn-effort-sprung");
    let f = w.model(
        "sprung.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 1s)\n(body cart (shape (box 100 30 30)) (attr mass 1kg))\n",
            "(assembly (axis (s 0)) (instance c cart (at) (slide x s) (effort-limits -1000N 1000N) (stiffness 100N/mm)))\n",
            "(assembly root (instance a (axis (* 100 (/ t 1s) (/ t 1s)))))\n(springs) (require lag (value (abs (following-error root/a/c)) < 0.01))\n"
        ),
    );
    let a = w.ink(["effort", &f]).ok()["actuators"][0].clone();
    assert!(a["basis"] == "spring" && close(n(&a["peak"]["effort"]), 0.4, 1e-3) && close(n(&a["rigid"]["peak"]["effort"]), 0.2, 1e-9), "{}", pretty(&a));
    let c = claim(&w.ink(["require", &f]).ok(), "effort-limits:root/a/c").clone();
    assert!((n(&c["value"]) - (n(&a["margin"]["margin"]) - n(&a["band"]))).abs() < 1e-6, "claim {} against {}", pretty(&c), pretty(&a));
    let at = w.ink(["effort", &f, "--t", "0.5"]).ok()["actuators"][0].clone();
    assert!(at["basis"] == "spring" && n(&at["effort"]) != n(&at["rigid"]["effort"]), "{}", pretty(&at));
}

#[test]
fn a_loose_motion_is_read_by_every_aggregate_and_says_partial() {
    // A 1 um, 1 kHz vibration of 1 kg: its steps reach the budget, loose,
    // and are read at their middles by the peak, the margin and the RMS
    // alike (m a peaks at 39.48 N, its RMS that over sqrt 2), partial.
    let w = Scratch::new("dyn-effort-loose");
    let f = w.model(
        "loose.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 0)\n(t-span 0 1s)\n(body cart (shape (box 50 50 50)) (attr mass 1kg))\n",
            "(assembly (axis (s 0)) (instance c cart (at) (slide x s) (effort-limits -50N 50N)))\n",
            "(assembly root (instance a (axis (* 0.001 (sin (* 2000 180deg/s t))))))\n"
        ),
    );
    let d = w.ink(["effort", &f]).ok();
    let (a, mo) = (&d["actuators"][0], &d["motion"]);
    let (peak, rms) = (n(&a["peak"]["effort"]).abs(), n(&a["rms"]));
    let want = 0.001 * (2000.0 * PI).powi(2) * 0.001;
    assert!(
        mo["partial"] == true && !mo["loose"].as_array().unwrap().is_empty() && peak >= rms && close(peak, want, 0.01) && close(rms, want / 2f64.sqrt(), 0.01),
        "want peak {want}: {}",
        pretty(&serde_json::json!({"actuator": a, "motion": mo}))
    );
}

#[test]
fn a_sprung_drives_limit_outage_is_in_its_main_answer() {
    // A sprung carriage's upper limit drops to 0.1 N for 0.2 ms around
    // 0.123456 s: its spring's effort passes it there, in the main answer
    // as under `rigid`, and the simulation says how it is stepped.
    let w = Scratch::new("dyn-effort-sprung-dip");
    let f = w.model(
        "dip.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 1s)\n(body cart (shape (box 100 30 30)) (attr mass 1kg))\n",
            "(assembly (axis (s 0)) (instance c cart (at) (slide x s) (effort-limits -1000N (* 1N (select (- (abs (- (/ t 1s) 0.123456)) 0.0001) 0.1 1000))) (stiffness 100N/mm)))\n",
            "(assembly root (instance a (axis (* 100 (/ t 1s) (/ t 1s)))))\n(springs) (require lag (value (abs (following-error root/a/c)) < 0.01))\n"
        ),
    );
    let d = w.ink(["effort", &f]).ok();
    let a = &d["actuators"][0];
    assert!(
        a["basis"] == "spring" && n(&a["margin"]["margin"]) < 0.0 && (n(&a["margin"]["t"]) - 0.123456).abs() < 1e-4 && a["band_kind"] == "measured" && a["partial"] == true,
        "{}",
        pretty(a)
    );
    let m = s(&d["simulated"]["method"]);
    assert!(m.contains("Radau IIA throughout") && !m.contains("Dormand"), "{m}");
}

#[test]
fn an_rms_over_an_impulse_is_null_and_its_step_unchecked() {
    // A step in the written motion is an impulse: its RMS is not found, so
    // none is given (ink-dynamics' an_rms_over_an_impulse_is_refused has why).
    let w = Scratch::new("dyn-effort-impulse");
    let jump = w.model(
        "jump.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 1s)\n(body block (shape (box 50 50 50)) (attr mass 1000g))\n",
            "(assembly (cart (x 0)) (instance b block (at) (slide x x) (effort-limits -100N 100N)))\n(assembly root (instance c (cart (select (- (/ t 1s) 0.5) 0 10))))\n"
        ),
    );
    let d = w.ink(["effort", &jump]).ok();
    assert!(d["actuators"][0]["rms"].is_null() && d["motion"]["unchecked"].as_array().is_some_and(|u| !u.is_empty()), "{}", pretty(&d));
}

#[test]
fn a_moving_body_without_mass_takes_its_own_effort_and_leaves_its_claim_unknown() {
    // None is guessed: the effort it would change is not given, naming the body; its rate still is, and the
    // gate's unpowered effort, which no such body is on, is what the file with the mass gives.
    let w = Scratch::new("dyn-effort-bare");
    let bare =
        w.model("bare.ink", &EFFORT.replace("(body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))", "(body bar (shape (move (box 20 300 10) 0 150 0)))"));
    let d = w.ink(["effort", &bare]).ok();
    let a = &d["actuators"][0];
    assert!(strs(&d["without_mass"]) == ["root/p/bar"] && s(&a["lacks"]).contains("no mass: bar (root/p/bar)") && a["peak"]["effort"].is_null(), "{}", pretty(&d));
    assert!((num(&a["peak"]["rate"]).abs() - 60.0 * std::f64::consts::PI).abs() < 1e-6, "{}", pretty(a));
    // (At one moment: over the span the steps are chosen for the efforts given, and so differ.)
    let (d, whole) = (w.ink(["effort", &bare, "--t", "0.3"]).ok(), w.ink(["effort", &w.model("whole.ink", EFFORT), "--t", "0.3"]).ok());
    assert!(!d["unpowered"][0]["effort"].is_null() && same(&d["unpowered"], whole["unpowered"].clone()), "{}", pretty(&d));
    assert!(d["actuators"][0]["effort"].is_null() && !whole["actuators"][0]["effort"].is_null() && same(&d["actuators"][0]["rate"], whole["actuators"][0]["rate"].clone()));
    let a = w.ink(["require", &bare]).ok()["requirements"][0].clone();
    assert!(verdict(&a) == "unknown" && s(&a["unavailable"]).contains("no mass") && a["motion"]["steps"] == 0, "{}", pretty(&a));
}

/// A study's load may read what the model works out: it is a number when the study is asked, the model's at
/// the moment of the motion where it is largest. EFFORT's pendulum hangs from a mount whose study is loaded
/// with the pendulum's joint's reaction; nothing is typed, and the load follows the design.
#[test]
fn a_studys_load_reads_the_models_reaction_at_its_worst_moment() {
    let w = Scratch::new("dyn-study-load");
    let mount = concat!(
        "(def held (reaction root/p/bar))\n",
        "(study hold (shape (box 30 30 60)) (elastic (E 2000MPa) (nu 0.35))\n",
        "  (fix base (- z -28))\n",
        "  (load pin (force (nth held 0) (nth held 1) (nth held 2)) (region (- 28 z))))\n",
        "(require stiff (max-displacement hold < 1))\n"
    );
    let file = w.model("mount.ink", &format!("{EFFORT}{mount}"));
    let pull = |d: &Value| d["loads"][0]["total"].as_array().unwrap().iter().map(|x| num(x).powi(2)).sum::<f64>().sqrt();
    // The load is the reaction's largest size over the span, as `ink effort` reads its peak (its own steps).
    let (d, e) = (w.ink(["solve", &file, "--res", "10"]).ok(), w.ink(["effort", &file]).ok());
    let peak = e["reactions"].as_array().unwrap().iter().find(|r| r["path"] == "root/p/bar").unwrap()["peak_force"]["force"].clone();
    let peak = peak.as_array().unwrap().iter().map(|x| num(x).powi(2)).sum::<f64>().sqrt();
    assert!(peak > 1.0 && (pull(&d) - peak).abs() < 0.01 * peak, "{} against {peak}: {}", pull(&d), pretty(&d));
    let said = s(&d["run_with"]);
    assert!(said.starts_with("run with what the model works out at t = ") && said.contains("of t 0 to 2 where it is largest: load 'pin' (N) "), "{said}");
    assert!(said.contains("a larger value between them, and a smaller one in a direction the part bears worse, are not looked at"), "{said}");
    // The bar twice as dense: twice the load, with nothing retyped.
    let dense = w.model(
        "dense.ink",
        &format!("{EFFORT}{mount}").replace("(move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3)", "(move (box 20 300 10) 0 150 0)) (attr density 5.4g/cm3)"),
    );
    let twice = w.ink(["solve", &dense, "--res", "10"]).ok();
    assert!((pull(&twice) - 2.0 * pull(&d)).abs() < 1e-6 * pull(&d), "{} against {}", pull(&twice), pull(&d));
    // The requirement answers from the same load, and says where it came from.
    let r = w.ink(["require", &file, "--name", "stiff"]).ok()["requirements"][0].clone();
    assert!(r["warnings"].to_string().contains("the least favourable of the 54 moments of the motion read") && !r["value"].is_null(), "{}", pretty(&r));
    // What the model cannot give is refused with what is missing: here, a gravity.
    let weightless = w.model("weightless.ink", &format!("{EFFORT}{mount}").replace("(gravity 0 0 -9.81m/s2)\n", ""));
    let e = w.ink(["solve", &weightless, "--res", "10"]).err()["error"].clone();
    assert!(s(&e["message"]).contains("gravity"), "{}", pretty(&e));
}

/// A study loaded from the model answers for the whole motion, not for the one moment its load is largest: a
/// 40 x 4 x 20 mm aluminium beam whose tip load swings from 20 N along its stiff axis to 10 N along its weak
/// one. By hand the weak 10 N bends it 10 x 40^3 / (3 x 68900 x 106.67) = 0.029 mm, the stiff 20 N 0.0023 mm.
/// The response is linear in the load, and every moment's load is a weighted sum of two moments': each
/// moment's deflection is read off that sum of two runs and decided on its own, with no band.
#[test]
fn a_study_loaded_from_the_model_answers_for_its_worst_moment_not_its_largest_load() {
    let w = Scratch::new("dyn-study-weak");
    let beam = concat!(
        "(import \"std\")\n(gravity 0 0 0)\n(t-span 0 1s)\n(body joint (shape (box 4 4 4)) (attr material \"steel\"))\n",
        "(assembly root (instance lever joint (at) (turn z 0)))\n",
        "(load wind root/lever (force 0 (* -10N (sin (* 90deg/s t))) (* -20N (cos (* 90deg/s t)))))\n",
        "(def wind (reaction root/lever))\n",
        "(study beam (shape (move (box 40 4 20) 20 0 0)) (elastic (material aluminum-6061)) (fix wall (- x 0.25))\n",
        "  (load tip (force 0 (nth wind 1) (nth wind 2)) (region (- 39.5 x))) (res 1))\n",
        "(require stiff (max-displacement beam < LIMIT))\n",
    );
    let asked = |limit: &str| w.ink(["require", &w.model(&format!("beam-{limit}.ink"), &beam.replace("LIMIT", limit))]).any()["requirements"][0].clone();
    // 0.02 mm: the largest load passes it tenfold and the weak moment fails it. No.
    let r = asked("0.02mm");
    assert!(verdict(&r) == "fails estimate", "{}", pretty(&r));
    let said = r["warnings"].to_string();
    assert!(said.contains("the least favourable of the 18 moments of the motion read") && said.contains("a weighted sum of those at 2 of them"), "{said}");
    assert!((0.025..0.031).contains(&num(&r["value"])) && r["band"].is_null(), "{}", pretty(&r));
    // 0.05 mm: every moment is under it. Yes.
    assert!(verdict(&asked("0.05mm")) == "holds estimate");
    // The same answer at the worst moment alone: the motion's worst is that moment's own run.
    let alone = w.ink(["require", &w.model("beam-alone.ink", &beam.replace("LIMIT", "0.02mm")), "--t-span", "0.99:1"]).any()["requirements"][0].clone();
    assert!((num(&alone["value"]) - num(&r["value"])).abs() < 0.01 * num(&r["value"]), "{} against {}", pretty(&alone), pretty(&r));
}

/// A study loaded from the model answers over a motion no better than over any part of it, by whatever
/// decides the answer: every moment read is decided on its own (its value, the stress its claim leaves
/// out, a buckle), and the least favourable is the answer. Three files where the whole motion said yes and
/// a part of it did not, each beside one the same but safe, which is still a yes.
#[test]
fn a_motion_answers_no_better_than_any_part_of_it() {
    let w = Scratch::new("dyn-study-parts");
    let head = concat!(
        "(import \"std\")\n(gravity 0 0 0)\n(t-span 0 1s)\n(body joint (shape (box 4 4 4)) (attr material \"steel\"))\n",
        "(assembly root (instance lever joint (at) (turn z 0)))\n",
    );
    let ask = |name: &str, body: &str, span: Option<&str>| {
        let file = w.model(name, &format!("{head}{body}"));
        let mut args = vec!["require", file.as_str()];
        args.extend(span.iter().flat_map(|s| ["--t-span", s]));
        w.ink(args).any()["requirements"][0].clone()
    };
    // The less favourable the answer, the lower: a whole motion's is never above a part's.
    let rank = |r: &Value| ["fails estimate", "unknown", "holds estimate"].iter().position(|v| *v == verdict(r)).unwrap_or_else(|| panic!("{}", pretty(r)));
    // A pressure that follows t beside the model's load: the study is another study at each moment, and no
    // moment's answer is a sum of other moments' runs. (It bends the bar 0.57 mm at t = 0.25; summed as if
    // it stood still, the motion read 0.059 mm.) On cells of 1: on 2 its grids read 0.06278 on 4 and 2 mm (the
    // peak at their last node inside the bar, y = 4) and 0.06285 on 1 (on its surface), moving apart: a real
    // difference of 1e-3, a thousand times the solver's precision, which it rightly does not settle. This test
    // asks of a motion against its parts, not of coarse grids, so it asks on 1, where they converge: 0.06278,
    // 0.06285 and 0.06286.
    let pressed = concat!(
        "(load wind root/lever (force 0 (* 10N (sin (* 180deg/s t))) 0))\n(def wind (reaction root/lever))\n",
        "(study bar (shape (move (box 100 10 10) 50 0 0)) (elastic (material aluminum-6061)) (fix wall (- x 0.25))\n",
        "  (load tip (force 0 (nth wind 1) 0) (region (- 99.5 x)))\n",
        "  (load p (pressure (* 1MPa PRESSED)) (region (max (- 4.5 y) (- 90 x)))) (res 1))\n",
        "(require stiff (max-displacement bar < 0.3mm))\n",
    );
    let follows = pressed.replace("PRESSED", "(* (sin (* 360deg/s t)) (sin (* 360deg/s t)))");
    let (whole, part) = (ask("pressed.ink", &follows, None), ask("pressed.ink", &follows, Some("0.24:0.26")));
    assert!(verdict(&whole) == "unknown" && s(&whole["unavailable"]).contains("it reads more than the model's loads change"), "{}", pretty(&whole));
    assert!(rank(&whole) <= rank(&part), "{} against {}", pretty(&whole), pretty(&part));
    assert!(verdict(&ask("pressed-still.ink", &pressed.replace("PRESSED", "0.01"), None)) == "holds estimate");
    // A column pushed along its length only at the motion's end: it buckles there (factor 0.61), and the
    // moment that bends it most is another. (The motion said yes on that one's 1.8 mm, its buckle unread.)
    // On cells of 0.5 mm, its hold and its load a whole millimetre each: there the bend's three grids settle
    // (1.77806 +- 0.00057). On the cells the study chose for it as first written (0.502 mm, held over 0.25
    // and loaded over 0.5) they agree to 1e-4 and move apart, 1.80129, 1.80144 and 1.80163, and that leaves
    // it unsettled before any buckle. (They read 1.797, 1.799 and 1.8005 there, converging, while each
    // sample a bar's edge crosses was counted by one side's plane: `occupied_parts`.)
    let column = concat!(
        "(load wind root/lever (force (* PUSH (sin (* 90deg/s t))) (* 1N (cos (* 90deg/s t))) 0))\n(def wind (reaction root/lever))\n",
        "(study col (shape (move (box 200 4 4) 100 0 0)) (elastic (material aluminum-6061)) (fix wall (- x 1))\n",
        "  (load tip (force (nth wind 0) (nth wind 1) 0) (region (- 199 x))) (res 0.5))\n",
        "(require stiff (max-displacement col < 3mm))\n",
    );
    let pushed = column.replace("PUSH", "150N");
    let (whole, part) = (ask("column.ink", &pushed, None), ask("column.ink", &pushed, Some("0.99:1")));
    assert!(verdict(&whole) == "unknown" && whole["undecided"].to_string().contains("a buckle"), "{}", pretty(&whole));
    assert!(verdict(&part) == "unknown" && rank(&whole) <= rank(&part), "{} against {}", pretty(&whole), pretty(&part));
    assert!(verdict(&ask("column-light.ink", &column.replace("PUSH", "1.5N"), None)) == "holds estimate");
    // A bar bent by 100 N, then pulled by 5180 N: the pull's stress away from the wall is the smaller, and
    // its stress by the wall, which the claim leaves out, the larger (70.97 MPa against the bend's 66.55).
    // (The motion read what it leaves out at the bend's moment alone.)
    let bar = concat!(
        "(load wind root/lever (force (* -5180N (max 0 (- 0 (cos (* 180deg/s t))))) (* 100N (max 0 (cos (* 180deg/s t)))) 0))\n(def wind (reaction root/lever))\n",
        "(study bar (shape (move (box 100 10 10) 50 0 0)) (elastic (material aluminum-6061)) (fix wall (- x 0.25))\n",
        "  (load tip (force (nth wind 0) (nth wind 1) 0) (region (- 99.5 x))) (res 2))\n",
        "(require strong (max-stress bar < LIMIT))\n",
    );
    let tight = bar.replace("LIMIT", "72MPa");
    let (whole, part) = (ask("bar.ink", &tight, None), ask("bar.ink", &tight, Some("0.6:1")));
    assert!(
        verdict(&whole) == "unknown" && whole["undecided"].to_string().contains("the support margin") && whole["undecided"].to_string().contains("70.97"),
        "{}",
        pretty(&whole)
    );
    assert!(verdict(&part) == "unknown" && rank(&whole) <= rank(&part), "{} against {}", pretty(&whole), pretty(&part));
    assert!(verdict(&ask("bar-roomy.ink", &bar.replace("LIMIT", "90MPa"), None)) == "holds estimate");
}

/// What else loads a study (a typed load, gravity) is held once at every moment, not once with each of
/// the model's loads. The same beam with 1 N typed on its tip along its weak axis, and the model's load
/// along that axis rising to 50 N and falling back to nothing: 51 N at the worst moment with the 1 N one
/// way and 49 N with it the other, so the two answers are as 49 to 51. (Run once per load component with
/// the typed load still on, the 1 N the other way cancelled the component's own 1 N, every moment's
/// bound came out the same, and the answer was a moment off the peak: 0.1290 mm for 0.1339, with no band.)
#[test]
fn what_else_loads_a_study_is_held_once_at_every_moment() {
    let w = Scratch::new("dyn-study-besides");
    let beam = concat!(
        "(import \"std\")\n(gravity 0 0 0)\n(t-span 0 1s)\n(body joint (shape (box 4 4 4)) (attr material \"steel\"))\n",
        "(assembly root (instance lever joint (at) (turn z 0)))\n",
        "(load wind root/lever (force 0 (* -50N (sin (* 180deg/s t))) 0))\n",
        "(def wind (reaction root/lever))\n",
        "(study beam (shape (move (box 40 4 20) 20 0 0)) (elastic (material aluminum-6061)) (fix wall (- x 0.25))\n",
        "  (load tip (force 0 (nth wind 1) (nth wind 2)) (region (- 39.5 x)))\n",
        "  (load dead (force 0 DEAD 0) (region (- 39.5 x))) (res 1))\n",
        "(require stiff (max-displacement beam < 0.02mm))\n",
    );
    let [with, against] = ["1N", "-1N"].map(|dead| {
        let r = w.ink(["require", &w.model(&format!("beam{dead}.ink"), &beam.replace("DEAD", dead))]).any()["requirements"][0].clone();
        assert!(verdict(&r) == "fails estimate" && r["band"].is_null(), "{dead}: {}", pretty(&r));
        num(&r["value"])
    });
    assert!((0.13..0.15).contains(&with) && (against / with - 49.0 / 51.0).abs() < 0.005, "{with} and {against}");
}

#[test]
fn an_actuator_whose_drive_reads_no_coordinate_is_refused_saying_how_to_write_one() {
    let w = Scratch::new("dyn-effort-no-coordinate");
    let f = w.model(
        "no-coordinate.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(body bar (shape (box 10 10 10)) (attr density 1g/cm3))\n",
            "(assembly root (instance j bar (at) (turn z (* 90deg/s t)) (effort-limits -1N*mm 1N*mm)))\n"
        ),
    );
    let e = w.ink(["effort", &f]).err()["error"].clone();
    assert!(s(&e["message"]).contains("reads no coordinate") && s(&e["hint"]).contains("(assembly (arm (angle 0))"), "{}", pretty(&e));
}

// A 1 kg cube on a slide held by 100 N/mm rings at (1/2 pi) sqrt(k / m),
// 50.33 Hz; no gravity is needed.
const MODES: &str = concat!(
    "(import \"std\")\n(body b (shape (box 100 100 100)) (attr density 1g/cm3))\n",
    "(assembly (axis (s 0)) (instance b b (at) (slide x s) (stiffness 100N/mm)))\n",
    "(assembly root (instance a (axis (* 10 (/ t 1s)))))\n",
);

#[test]
fn modes_give_a_mass_on_a_spring_its_natural_frequency_and_the_spring() {
    let w = Scratch::new("dyn-modes");
    let f = w.model("modes.ink", MODES);
    let d = w.ink(["effort", &f, "--modes", "--t", "0.5"]).ok();
    let m = &d["modes"];
    let want = (1e6 * 100.0 / 1000.0f64).sqrt() / (2.0 * PI);
    let sp = &d["springs"][0];
    assert!(
        m.as_array().is_some_and(|x| x.len() == 1)
            && s(&d["method"]).contains("linearised")
            && close(n(&m[0]["frequency"]), want, 2e-2)
            && same(&m[0]["moves"], serde_json::json!([{"joint": "root/a/b", "motion": 1, "energy": 1}]))
            && same(&d["pose"], serde_json::json!([{"name": "a.s", "value": 5}]))
            && d["springs"].as_array().is_some_and(|x| x.len() == 1)
            && sp["what"] == "joint root/a/b"
            && n(&sp["stiffness"]) == 100.0
            && sp["per"].is_null()
            && s(&sp["site"]).ends_with("/modes.ink:3")
            && sp.get("units").is_none()
            && d["accuracy"].is_null()
            && ["basis", "left_out", "units"].iter().all(|k| d.get(k).is_none()),
        "want {want} Hz: {}",
        pretty(&d)
    );
}

#[test]
fn modes_where_nothing_states_a_stiffness_are_refused_saying_how_to_give_one() {
    // (The refusal itself is modes.rs a_joint_held_by_its_stiffness_rings_and_nothing_held_is_refused; here, its hint.)
    let w = Scratch::new("dyn-modes-unheld");
    let f = w.model("modes-unheld.ink", &MODES.replace(" (stiffness 100N/mm)", ""));
    let e = w.ink(["effort", &f, "--modes"]).err()["error"].clone();
    assert!(s(&e["message"]).contains("states a stiffness") && s(&e["hint"]).contains("(stiffness K)"), "{}", pretty(&e));
}

#[test]
fn a_parts_stepper_is_held_by_the_stiffness_its_file_states_where_it_states_it() {
    // 17HS19-2004S1's 0.59 N m x 50 teeth, stated in the imported fixture.
    let w = Scratch::new("dyn-modes-stepper");
    let held = w.model(
        "modes-stepper.ink",
        &format!(
            concat!(
                "(import \"std\")\n(import \"{}\" as m)\n",
                "(body arm (shape (move (box 100 10 10) 50 0 0)) (attr density 2.7g/cm3))\n",
                "(assembly root (instance d (m.drive 0)) (instance arm arm (frame-of d shaft)))\n"
            ),
            fixture("dynamics/stepper.ink")
        ),
    );
    let d = w.ink(["effort", &held, "--modes"]).ok();
    let motor: Vec<&Value> = d["springs"].as_array().unwrap().iter().filter(|x| x["what"] == "joint root/d/rotor").collect();
    assert!(
        d["modes"].as_array().is_some_and(|x| x.len() == 1)
            && motor.len() == 1
            && n(&motor[0]["stiffness"]) == 29500.0
            && motor[0]["per"] == "radian"
            && s(&motor[0]["site"]).ends_with("tests/fixtures/dynamics/stepper.ink:15"),
        "{}",
        pretty(&d)
    );
}

#[test]
fn a_reaction_claim_reads_in_a_turning_frame_i_alpha_plus_omega_cross_i_omega() {
    // P9's plate turning about its diagonal, spinning up: a claim reads its
    // reaction in its own (moving) frame, against omega x I omega + I alpha
    // at the inertia ink measures.
    let w = Scratch::new("dyn-reaction-plate");
    let f = w.model(
        "plate.ink",
        concat!(
            "(import \"std\")\n(gravity 0 0 0)\n",
            "(body plate (shape (box 300 200 5)) (attr density 2.7g/cm3)\n  (frame diag (axes (origin 0 0 0) (x 0 0 1) (z -300 -200 0))))\n",
            "(assembly (spinner (a 0)) (instance plate plate (at (rot-x 30) (rot-z 20)) diag (turn z a) (effort-limits -1e6N*mm 1e6N*mm)))\n",
            "(assembly root (instance s (spinner (rad->deg (+ (* 100 (/ t 1s)) (* 25 (/ t 1s) (/ t 1s)))))))\n",
            "(def r (reaction root/s/plate (frame-of root/s/plate)))\n",
            "(require mx (value (nth r 3) > 0))\n(require my (value (nth r 4) > 0))\n(require mz (value (nth r 5) > 0))\n"
        ),
    );
    let i = ink(["mass", &f, "--target", "plate", "--keys", "inertia"])["targets"][0]["inertia"].clone();
    let h = 300f64.hypot(200.0);
    let u = [300.0 / h, 200.0 / h, 0.0];
    let iu: Vec<f64> = (0..3).map(|a| (0..3).map(|b| n(&i[a][b]) * u[b]).sum()).collect();
    let cross = [u[1] * iu[2] - u[2] * iu[1], u[2] * iu[0] - u[0] * iu[2], u[0] * iu[1] - u[1] * iu[0]];
    let want: Vec<f64> = (0..3).map(|k| 1e-6 * (50.0 * iu[k] + 100f64.powi(2) * cross[k])).collect();
    // Each holds over the band its inertia's proven bounds give it; its value is read with the
    // claim turned (failing over all of its band, it answers the value itself).
    let d = w.ink(["require", &f]).ok();
    let mz = claim(&d, "mz");
    assert!(
        ["mx", "my", "mz"].iter().all(|c| verdict(claim(&d, c)) == "holds estimate")
            && mz["accuracy"].is_null()
            && strs(&mz["solved_inputs"]).contains(&"reaction:root/s/plate.mz"),
        "{}",
        pretty(&d)
    );
    let raw = w.model("plate-turned.ink", &std::fs::read_to_string(&f).unwrap().replace(") > 0))", ") < 0))"));
    let d = w.ink(["require", &raw]).any();
    let got: Vec<f64> = ["mx", "my", "mz"].iter().map(|c| n(&claim(&d, c)["value"])).collect();
    assert!(got.iter().zip(&want).all(|(g, w)| close(*g, *w, 1e-6)), "got {got:?}, want {want:?}");
}

/// An arm on a base swinging a block that rides its tip frame, a payload
/// load at the tip, swung through 60 degrees with this period; a claim that
/// the floor's zero-moment point stays on the foot.
fn arm(period: f64) -> String {
    format!(
        concat!(
            "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 2s)\n",
            "(body base (shape (move (box 120 80 20) 0 0 10)) (attr density 2.7g/cm3))\n",
            "(body link (shape (move (box 150 10 10) 75 0 0)) (attr density 2.7g/cm3) (frame tip (at (move 150 0 0))))\n",
            "(body block (shape (box 20 20 20)) (attr density 7.85g/cm3))\n",
            "(assembly (arm (theta 0))\n  (instance base base)\n  (instance link link (at (move 0 0 100)) (turn y theta) (effort-limits -1e5N*mm 1e5N*mm))\n  (instance block block (frame-of link tip)))\n",
            "(assembly root (instance a (arm (* 60 (sin (/ (* 2 pi (/ t 1s)) {}))))))\n",
            "(load payload root/a/link (force 0 0 -0.981N) (at (frame-of root/a/link tip)))\n",
            "(def foot (max (- (abs x) 60) (- (abs y) 40)))\n(def p (zmp (reaction root (at))))\n",
            "(require stands (value (warp foot (nth p 0) (nth p 1) 0) < 0))\n"
        ),
        period
    )
}

#[test]
fn effort_gives_each_joints_reaction_carrying_what_rides_it_and_the_load_on_it() {
    let w = Scratch::new("dyn-reaction-arm");
    let slow = w.model("slow.ink", &arm(4.0));
    let d = w.ink(["effort", &slow, "--t", "0"]).ok();
    let paths: Vec<&str> = d["reactions"].as_array().unwrap().iter().map(|x| s(&x["path"])).collect();
    let (link_r, root_r) = (by(&d["reactions"], "path", "root/a/link"), by(&d["reactions"], "path", "root"));
    let (block, link) = (7.85e-3 * 8000.0 * 1e-3, 2.7e-3 * 15000.0 * 1e-3); // kg
    let om = 60f64.to_radians() * 2.0 * PI / 4.0;
    assert!(
        paths == ["root", "root/a/link"]
            && close(n(&link_r["force"][0]), -(block * 150.0 + link * 75.0) * om * om / 1000.0, 0.02)
            && close(n(&link_r["force"][2]), (block + link) * 9.81 + 0.981, 0.02)
            && root_r.get("units").is_none()
            && d.get("reactions_note").is_none(),
        "{}",
        pretty(&d["reactions"])
    );
    // Over the span: each reaction where its force and its moment are
    // largest, with the whole wrench there.
    let fast = w.model("fast.ink", &arm(0.4));
    let d = w.ink(["effort", &fast]).ok();
    let pk = by(&d["reactions"], "path", "root/a/link");
    let whole = |p: &Value| ["t", "force", "moment"].iter().all(|k| p.get(k).is_some());
    assert!(whole(&pk["peak_force"]) && whole(&pk["peak_moment"]) && pk["peak_force"]["moment"].as_array().is_some_and(|m| m.len() == 3), "{}", pretty(pk));
}

#[test]
fn the_zero_moment_point_stays_on_the_foot_slowly_not_fast_nor_lifted() {
    // Tipping (std zmp): the floor's zero-moment point within the footprint
    // over the motion. The slow swing stands, an estimate; the fast one
    // throws it past the foot's edge, with the moment and how far; lifted
    // off its support (the floor pulls up), it is at infinity.
    let w = Scratch::new("dyn-reaction-zmp");
    let slow = w.model("slow.ink", &arm(4.0));
    let fast = w.model("fast.ink", &arm(0.4));
    let lift = w.model("lift.ink", &arm(4.0).replace("(force 0 0 -0.981N)", "(force 0 0 1000N)"));
    let st = |f: &str| claim(&w.ink(["require", f]).ok(), "stands").clone();
    let (sl, fa, li) = (st(&slow), st(&fast), st(&lift));
    assert!(verdict(&sl) == "holds estimate" && n(&sl["value"]) < 0.0 && sl["accuracy"].is_null() && n(&sl["motion"]["steps"]) > 16.0, "slow: {}", pretty(&sl));
    assert!(verdict(&fa) == "fails estimate" && n(&fa["value"]) > 0.0 && (0.0..=2.0).contains(&n(&fa["t"])), "fast: {}", pretty(&fa));
    assert!(verdict(&li) == "fails estimate" && li["value"].is_null(), "lifted: {}", pretty(&li));
}

#[test]
fn an_imported_files_reaction_claim_is_answered_in_that_file_not_here() {
    let w = Scratch::new("dyn-reaction-import");
    w.model("slow.ink", &arm(4.0));
    let top = w.model("top.ink", "(import \"std\")\n(import \"slow.ink\" as s)\n(assembly root (instance a s.root))\n");
    let c = claim(&w.ink(["require", &top, "--imports"]).ok(), "s.stands").clone();
    assert!(verdict(&c) == "unknown" && s(&c["unavailable"]).contains("check the imported file itself"), "{}", pretty(&c));
}

#[test]
fn a_reaction_its_claims_would_refuse_has_no_numbers_saying_why() {
    // An unpowered hinge moves what top carries; a ground half-space is not
    // proven bounded. Each refuses top's reaction, not the arm's joint's.
    let w = Scratch::new("dyn-reaction-refused");
    let gate = w.model(
        "gate.ink",
        &(arm(4.0)
            + concat!(
                "(body leaf (shape (move (box 60 4 80) 30 0 -40)) (attr density 2.7g/cm3))\n",
                "(assembly (flap (angle 0)) (instance leaf leaf (at (move 0 0 120)) (turn y angle)))\n",
                "(assembly top (instance a (arm 0)) (instance g (flap (* 20 (sin (* 180deg/s t)))) (at (move 0 200 0))))\n"
            )),
    );
    let r = w.ink(["effort", &gate, "--target", "top", "--t", "0.25"]).ok()["reactions"].clone();
    let (top, joint) = (by(&r, "path", "top"), by(&r, "path", "top/a/link"));
    assert!(top["force"].is_null() && s(&top["refused"]).contains("no actuator supplies") && !joint["force"].is_null(), "unpowered: {}", pretty(&r));
    let floor = w.model("floor.ink", &(arm(4.0) + "(body ground (shape (+ z 50)) (attr density 1g/cm3))\n(assembly top (instance a (arm 0)) (instance g ground))\n"));
    let r = w.ink(["effort", &floor, "--target", "top", "--t", "0"]).ok()["reactions"].clone();
    let top = by(&r, "path", "top");
    assert!(strs(&top["unbounded"]) == ["top/g"] && top["force"].is_null() && s(&top["refused"]).contains("not proven bounded"), "unbounded: {}", pretty(top));
}

#[test]
fn a_reaction_names_what_moves_with_its_joint_unseparated() {
    // A block moving with the link but not placed on its frames: the
    // joint's reaction (moves_with) and a claim on it (its notes) name it.
    let w = Scratch::new("dyn-reaction-loose");
    let loose = arm(4.0).replace("(instance block block (frame-of link tip))", "(instance block block (at (move 150 0 100) (rot-y theta)))");
    let f = w.model("loose.ink", &loose);
    let r = w.ink(["effort", &f, "--t", "0.5"]).ok()["reactions"].clone();
    assert!(strs(&by(&r, "path", "root/a/link")["moves_with"]) == ["root/a/block"] && by(&r, "path", "root").get("moves_with").is_none(), "{}", pretty(&r));
    // Its push on the link is summed in, undetermined: a claim on the
    // reaction is unknown, saying so, never a hold on the sum.
    let g = w.model("loose-claim.ink", &(loose + "(require link-fz (value (nth (reaction root/a/link) 2) > 0))\n"));
    let c = claim(&w.ink(["require", &g]).ok(), "link-fz").clone();
    assert!(verdict(&c) == "unknown" && s(&c["unavailable"]).contains("root/a/block") && s(&c["unavailable"]).contains("moves_with"), "{}", pretty(&c));
}

/// A slider-crank:
/// 600 rpm, a 200 N gas load, its motor limited to 9150 N mm. Its peak
/// demand is 9169.65 N mm (a 4000-point energy balance) between the
/// moments the steps measure (9102.27 at the best of them): the peak is
/// found where the effort turns back, and the claim fails.
const SLIDER_TIGHT: &str = r#"(import "std")
(gravity 0 0 0)
(t-span 0 0.1s)
(def r 40) (def l 120)
(fn bar (len) (move (box len 10 5) (/ len 2) 0 0))
(body crank (shape (bar r)) (attr density 7.85g/cm3))
(body rod (shape (bar l)) (attr density 7.85g/cm3))
(body piston (shape (box 20 20 80)) (attr mass 500g))
(fn bx (q) (* r (cos (deg->rad q))))
(fn by (q) (* r (sin (deg->rad q))))
(fn px (q) (+ (bx q) (sqrt (- (* l l) (square (by q))))))
(assembly (sc (q 0))
  (instance crank crank (at) (turn z q) (effort-limits -9150N*mm 9150N*mm))
  (instance rod rod (at (rot-z (rad->deg (atan2 (neg (by q)) (- (px q) (bx q))))) (move (bx q) (by q) 0)))
  (instance piston piston (at) (slide x (px q))))
(assembly root (instance m (sc (* 3600 (/ t 1s)))))
(load gas root/m/piston (force -200N 0 0))
"#;

#[test]
fn a_peak_between_the_steps_is_found_where_the_effort_turns_back() {
    let w = Scratch::new("dyn-slider-peak");
    let f = w.model("tight.ink", SLIDER_TIGHT);
    let peak = n(&w.ink(["effort", &f]).ok()["actuators"][0]["peak"]["effort"]);
    // As large pushing as pulling: the earlier of the two is shown.
    assert!(close(peak.abs(), 9169.65, 1e-6), "the peak {peak} against 9169.65");
    let c = claim(&w.ink(["require", &f]).ok(), "effort-limits:root/m/crank").clone();
    assert!(c["status"] == "fails" && n(&c["value"]) < -19.0, "the 9150 N mm limit: {}", pretty(&c));
}

/// A gearbox of two
/// spur meshes, 9:1, lifting an arm at constant speed at t = 1. At 0.5 a
/// mesh the motor gives the ideal torque over 0.5 x 0.5 (350.1 N mm, not
/// the 3 times a sum of each stage's loss gives), so a 300 N mm limit
/// fails; lowering, the load drives the motor through both, which brakes
/// 0.25 of the ideal (4.55 N mm, not 0).
const GEARBOX: &str = r#"(import "std")
(gravity 0 -9.81m/s2 0)
(t-span 0 2s)
(fn gear (z) (subtract (extrude-z (involute-gear-2d z 1 20 0) -4 4) (cylinder-z 2 20)))
(body rotor (shape (cylinder-z 5 30)) (attr density 7.85g/cm3) (attr spin-inertia 5000g*mm2))
(body pinion (shape (gear 20)) (attr density 7.85g/cm3))
(body g60 (shape (gear 60)) (attr density 7.85g/cm3))
(body g15 (shape (gear 15)) (attr density 7.85g/cm3))
(body g45 (shape (gear 45)) (attr density 7.85g/cm3))
(body arm (shape (move (box 150 10 6) 75 0 0)) (attr density 7.85g/cm3) (frame tip (at (move 150 0 0))))
(body block (shape (box 20 20 80)) (attr mass 500g))
(assembly mid (instance a g60) (instance b g15 (at (move 0 0 10))))
(assembly out (instance g g45 (at (move 0 0 10))) (instance arm arm (at (move 0 0 25))) (instance block block (at (move 0 0 25) (frame-of arm tip))))
(fn mid-angle (a) (mesh-turn 20 60 0 a))
(fn out-angle (a) (mesh-turn 15 45 0 (mid-angle a)))
(assembly (box9 (a 0))
  (instance m rotor (at (move 0 0 -30)) (turn z a) (effort-limits -300N*mm 300N*mm))
  (instance p pinion (at (move 0 0 30) (frame-of m)))
  (instance i mid (at (move 40 0 0)) (turn z (mid-angle a)))
  (instance o out (at (move 70 0 0)) (turn z (out-angle a)))
  (relate p i/a (gear-mesh 20 60 1 20 0 E))
  (relate i/b o/g (gear-mesh 15 45 1 20 0 E)))
(fn cyc (amp T) (* amp (- (/ (/ t 1s) T) (/ (sin (* 2 pi (/ (/ t 1s) T))) (* 2 pi)))))
(assembly root (instance g (box9 DRIVE)))
"#;

/// Each stage's tooth force is what passes through it, the losses of the
/// stages it drives among it: lifting, the first mesh carries the second's
/// loss, 1 / 0.5 of its ideal tooth force, and the second its ideal one;
/// lowering, the load drives the second through 0.5 and the first through
/// 0.5 again, 0.5 and 0.25 of their ideal tooth forces.
#[test]
fn losses_of_stages_in_series_multiply_either_way_power_flows() {
    let w = Scratch::new("dyn-gearbox-stages");
    // The actuator's effort, and each mesh's tooth force (its loaded flank).
    let effort = |name: &str, e: &str, drive: &str| {
        let f = w.model(name, &GEARBOX.replace(" E)", &format!(" {e})")).replace("DRIVE", drive));
        let d = w.ink(["effort", &f, "--t", "1"]).ok();
        let tooth = |k: usize| d["couplings"][k]["tensions"].as_array().unwrap().iter().map(n).fold(0.0, f64::max);
        (n(&d["actuators"][0]["effort"]), [tooth(0), tooth(1)])
    };
    let (lift, lower) = ("(- (* 9 (cyc 90 2)) 1512)", "(* 9 (cyc -90 2))");
    let ((ideal_up, ideal_teeth), (lossy_up, teeth)) = (effort("up1.ink", "1", lift), effort("up.ink", "0.5", lift));
    assert!(close(ideal_up, 787.73 / 9.0, 1e-3), "the ideal lift {ideal_up}");
    assert!(close(lossy_up, ideal_up / 0.25, 1e-9), "lifting through two 0.5 meshes: {lossy_up} against {}", ideal_up / 0.25);
    assert!(close(teeth[0], ideal_teeth[0] / 0.5, 1e-9), "lifting, the first mesh passes the second's loss: {teeth:?} against {ideal_teeth:?}");
    assert!(close(teeth[1], ideal_teeth[1], 1e-9), "and the second its ideal tooth force: {teeth:?} against {ideal_teeth:?}");
    let ((ideal_down, ideal_teeth), (lossy_down, teeth)) = (effort("down1.ink", "1", lower), effort("down.ink", "0.5", lower));
    assert!(ideal_down > 0.0 && close(lossy_down, ideal_down * 0.25, 1e-9), "lowering through two 0.5 meshes: {lossy_down} against {}", ideal_down * 0.25);
    assert!(close(teeth[1], ideal_teeth[1] * 0.5, 1e-9), "lowering, the second mesh passes 0.5: {teeth:?} against {ideal_teeth:?}");
    assert!(close(teeth[0], ideal_teeth[0] * 0.25, 1e-9), "and the first 0.25: {teeth:?} against {ideal_teeth:?}");
    let f = w.model("claim.ink", &GEARBOX.replace(" E)", " 0.5)").replace("DRIVE", lift));
    let c = claim(&w.ink(["require", &f, "--t", "1"]).ok(), "effort-limits:root/g/m").clone();
    assert!(c["status"] == "fails" && close(n(&c["value"]), 300.0 - lossy_up, 1e-9), "the 300 N mm limit: {}", pretty(&c));
}

/// The same gearbox with a mesh written output-first: each mesh 0.5 both
/// ways (a gear mesh's EB is its E), so the same gearbox. Its loss is
/// charged at the end power enters, whichever end is written first:
/// lifting 350.1 N mm, lowering 4.55, as written input-first.
#[test]
fn a_mesh_written_output_first_is_charged_where_power_enters() {
    let w = Scratch::new("dyn-gearbox-reversed");
    let first = GEARBOX.replace("(relate p i/a (gear-mesh 20 60 1 20 0 E))", "(relate i/a p (gear-mesh 60 20 1 20 0 E))");
    let second = GEARBOX.replace("(relate i/b o/g (gear-mesh 15 45 1 20 0 E)))", "(relate o/g i/b (gear-mesh 45 15 1 20 0 E)))");
    assert!(first != GEARBOX && second != GEARBOX);
    let effort = |name: &str, text: &str, e: &str, drive: &str| {
        let f = w.model(name, &text.replace(" E)", &format!(" {e})")).replace("DRIVE", drive));
        n(&w.ink(["effort", &f, "--t", "1"]).ok()["actuators"][0]["effort"])
    };
    let (lift, lower) = ("(- (* 9 (cyc 90 2)) 1512)", "(* 9 (cyc -90 2))");
    let (up, down) = (effort("up1.ink", GEARBOX, "1", lift), effort("down1.ink", GEARBOX, "1", lower));
    for (k, text) in [("first", &first), ("second", &second)] {
        let got = effort(&format!("{k}-up.ink"), text, "0.5", lift);
        assert!(close(got, up / 0.25, 1e-9), "{k} reversed, lifting: {got} against {}", up / 0.25);
        let got = effort(&format!("{k}-down.ink"), text, "0.5", lower);
        assert!(close(got, down * 0.25, 1e-9), "{k} reversed, lowering: {got} against {}", down * 0.25);
    }
}

/// A four-bar, its
/// rocker placed by a formula of the crank's angle, as the catalog's
/// linkages are, a 50 N push on the coupler. The rocker's reaction is its
/// whole balance (0.26 N), not the 23 N its pivot carries through the
/// coupler's pin, which no joint declares: the claim is unknown, naming
/// the coupler, on a formula's placement and on a joint's alike.
const FOUR_BAR_PIN: &str = r#"(import "std")
(gravity 0 -9.81m/s2 0)
(t-span 0 1s)
(def ga 100) (def la 30) (def lb 90) (def lc 70)
(fn bar (l) (move (box l 10 5) (/ l 2) 0 0))
(body crank (shape (bar la)) (attr density 7.85g/cm3))
(body coupler (shape (bar lb)) (attr density 7.85g/cm3) (frame mid (at (move (/ lb 2) 0 0))))
(body rocker (shape (bar lc)) (attr density 7.85g/cm3) (frame tip (at (move lc 0 0))))
(fn pin-b (q) (vec (* la (cos (deg->rad q))) (* la (sin (deg->rad q)))))
(fn pin-c (q)
  (let ((b (pin-b q)) (dx (- ga (nth b 0))) (dy (neg (nth b 1)))
        (d (sqrt (+ (* dx dx) (* dy dy))))
        (a2 (/ (+ (* lb lb) (neg (* lc lc)) (* d d)) (* 2 d)))
        (h (sqrt (- (* lb lb) (* a2 a2)))))
    (vec (+ (nth b 0) (* a2 (/ dx d)) (* -1 h (/ dy d)))
         (+ (nth b 1) (* a2 (/ dy d)) (* h (/ dx d))))))
(fn coupler-deg (q) (let ((b (pin-b q)) (c (pin-c q))) (rad->deg (atan2 (- (nth c 1) (nth b 1)) (- (nth c 0) (nth b 0))))))
(fn rocker-deg (q) (let ((c (pin-c q))) (rad->deg (atan2 (nth c 1) (- (nth c 0) ga)))))
(assembly (fourbar (q 0))
  (instance crank crank (at) (turn z q) (effort-limits -2000N*mm 2000N*mm))
  (instance coupler coupler (at (rot-z (coupler-deg q)) (move (nth (pin-b q) 0) (nth (pin-b q) 1) 0)))
  (instance rocker rocker (at (rot-z (rocker-deg q)) (move ga 0 0))))
(assembly root (instance m (fourbar (* 360 (/ t 1s)))))
(load push root/m/coupler (force 0 -50N 0) (at (frame-of root/m/coupler mid)))
(def d (reaction root/m/rocker))
(require pivot-d (value (len2 (nth d 0) (nth d 1)) < 5N))
"#;

/// The four-bar with its rocker placed in the root by the same formula
/// written in t, not the crank's coordinate: no joint or actuator moves or
/// holds it, so its reaction is unknown, and the crank's effort, which no
/// longer counts it, says so.
#[test]
fn a_part_placed_by_t_through_no_coordinate_is_held_by_nothing_declared() {
    let w = Scratch::new("dyn-four-bar-pin-t");
    let text = FOUR_BAR_PIN
        .replace("\n  (instance rocker rocker (at (rot-z (rocker-deg q)) (move ga 0 0))))", ")")
        .replace(
            "(assembly root (instance m (fourbar (* 360 (/ t 1s)))))",
            "(assembly root (instance m (fourbar (* 360 (/ t 1s)))) (instance rocker rocker (at (rot-z (rocker-deg (* 360 (/ t 1s)))) (move ga 0 0))))",
        )
        .replace("(def d (reaction root/m/rocker))", "(def d (reaction root/rocker))");
    assert!(text.contains("(instance rocker rocker (at (rot-z (rocker-deg (* 360 (/ t 1s))))") && !text.contains("root/m/rocker") && !text.contains("(rocker-deg q)"), "{text}");
    let f = w.model("t.ink", &text);
    for span in [&["--t", "0.25"][..], &[]] {
        let mut args = vec!["require", f.as_str()];
        args.extend(span);
        let d = w.ink(args).ok();
        let c = claim(&d, "pivot-d");
        let why = s(&c["unavailable"]);
        assert!(verdict(c) == "unknown" && why.contains("root/rocker") && why.contains("driven_by_nothing"), "{span:?}: {}", pretty(c));
        let e = claim(&d, "effort-limits:root/m/crank");
        assert!(notes_have(e, &["root/rocker", "driven_by_nothing"]), "{span:?}: {}", pretty(e));
    }
}

#[test]
fn a_linkages_pivot_reaction_summing_its_pins_push_is_unknown() {
    let w = Scratch::new("dyn-four-bar-pin");
    let joint = FOUR_BAR_PIN.replace("(instance rocker rocker (at (rot-z (rocker-deg q)) (move ga 0 0)))", "(instance rocker rocker (at (move ga 0 0)) (turn z (rocker-deg q)))");
    assert_ne!(joint, FOUR_BAR_PIN);
    for (name, text) in [("formula.ink", FOUR_BAR_PIN), ("joint.ink", joint.as_str())] {
        let f = w.model(name, text);
        for span in [&["--t", "0.25"][..], &[]] {
            let mut args = vec!["require", f.as_str()];
            args.extend(span);
            let d = w.ink(args).ok();
            let c = claim(&d, "pivot-d");
            let why = s(&c["unavailable"]);
            assert!(verdict(c) == "unknown" && why.contains("root/m/coupler") && why.contains("moves_with"), "{name} {span:?}: {}", pretty(c));
            assert!(d["all_yes"] != true, "{name} {span:?}: {}", pretty(&d));
        }
    }
}

#[test]
fn a_carriages_share_says_where_its_slides_drive_pulls_unless_a_belt_is_coupled() {
    let w = Scratch::new("dyn-reaction-carriage");
    let car = fixture("dynamics/carriage.ink");
    // Its share of its slide's reaction rests on where the slide's drive
    // pulls, which no part named says: its claim says so. (The value is
    // ink-dynamics' a_rail_carriage_under_an_offset_load_by_hand.)
    let f = w.model(
        "carriage.ink",
        &format!(
            concat!(
                "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 1s)\n(import \"{car}\" as car)\n",
                "(body block (shape (box 50 50 50)) (attr mass 1500g))\n",
                "(assembly (axis (s 0)) (instance rail car.rail) (instance c car.carriage (frame-of rail) (slide x s) (effort-limits -100N 100N))\n",
                "  (instance block block (at (move 0 100 33) (frame-of c))))\n",
                "(assembly root (instance a (axis (+ 150 (* 50 (sin (* 360deg/s t)))))))\n",
                "(require static (value (carriage-static-ratio (nth (reaction root/a/c (by root/a/c)) 0) car.spec) < 0.5))\n"
            ),
            car = car
        ),
    );
    let c = claim(&w.ink(["require", &f, "--t", "0.25"]).ok(), "static").clone();
    assert!(verdict(&c) == "holds estimate" && notes_have(&c, &["slide's drive is taken along its axis"]), "{}", pretty(&c));
    // Its value, the claim turned (it holds over its masses' band, and answers that band's edge).
    let raw = w.model("carriage-turned.ink", &std::fs::read_to_string(&f).unwrap().replace("< 0.5))", "> 0.5))"));
    let c = claim(&w.ink(["require", &raw, "--t", "0.25"]).any(), "static").clone();
    assert!(close(n(&c["value"]), 1471.5 / 38220.0, 1e-6), "{}", pretty(&c));
    // Tied by a belt related as a coupling (0.95), 20 mm above it, a 100 N
    // push driving it along, the carriage holds the belt's line, and its
    // claim no longer says where the drive pulls. Power enters the belt at
    // the carriage: the belt passes 0.95 of it, 95 N, to the drum, whose
    // joint passes the 950 N mm its effort brakes; the other 5 N is the
    // belt's friction on the carriage, paired with the drum on its axis 30
    // mm up: 95 x 20 + 5 x 30 = 2050 N mm of pitch.
    let b = w.model(
        "carriage-belt.ink",
        &format!(
            concat!(
                "(import \"std\")\n(gravity 0 0 0)\n(t-span 0 1s)\n(import \"{car}\" as car)\n",
                "(body drum (shape (cylinder-z 10 8)) (attr density 7.85g/cm3))\n",
                "(assembly (axis (s 0)) (instance rail car.rail) (instance c car.carriage (frame-of rail) (slide x s))\n",
                "  (instance d drum (at (rot-x -90) (move -200 0 30) (frame-of rail)) (turn z (rad->deg (/ s -10))) (effort-limits -1e4N*mm 1e4N*mm))\n",
                "  (relate d c (belt-span (vec (+ (get (frame-of rail) \"origin\") (vec -200 0 20)) (+ (get (frame-of c) \"origin\") (vec 0 0 20))) 0 0.95)))\n",
                "(assembly root (instance a (axis (+ 150 (* 10 (/ t 1s))))))\n(load drag root/a/c (force 100N 0 0))\n",
                "(require pitch (value (abs (nth (nth (reaction root/a/c (by root/a/c)) 0) 4)) > 0))\n",
                "(require pull (value (nth (tension root/a/d root/a/c) 0) > 0))\n"
            ),
            car = car
        ),
    );
    let d = w.ink(["require", &b, "--t", "0.25"]).ok();
    assert!(verdict(claim(&d, "pitch")) == "holds estimate" && verdict(claim(&d, "pull")) == "holds estimate", "{}", pretty(&d));
    let d = w.ink(["require", &w.model("carriage-belt-turned.ink", &std::fs::read_to_string(&b).unwrap().replace(") > 0))", ") < 0))")), "--t", "0.25"]).any();
    let (pitch, pull) = (claim(&d, "pitch"), claim(&d, "pull"));
    let e = w.ink(["effort", &b, "--t", "0.25"]).ok();
    let (drum, joint) = (by(&e["actuators"], "path", "root/a/d"), by(&e["reactions"], "path", "root/a/d"));
    assert!(
        close(n(&pitch["value"]), 2050.0, 1e-9)
            && close(n(&pull["value"]), 95.0, 1e-9)
            && !notes_have(pitch, &["slide's drive"])
            && close(n(&drum["effort"]), n(&joint["moment"][2]), 1e-6)
            && n(&drum["effort"]).abs() > 900.0,
        "{}\n{}\ndrum {}\njoint {}",
        pretty(pitch),
        pretty(pull),
        pretty(drum),
        pretty(joint)
    );
}

/// A disc spinning off-balance (a hole 30 mm out) at 100 rev/s.
const DISC: &str = concat!(
    "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 1s)\n",
    "(body disc (shape (subtract (cylinder-z 50 10) (move (cylinder-z 5 20) 30 0 0))) (attr density 7.85g/cm3))\n",
    "(assembly (spinner (a 0)) (instance disc disc (at (rot-x 90)) (turn z a) (effort-limits -1e4N*mm 1e4N*mm)))\n",
    "(assembly root (instance s (spinner (* 360 100 (/ t 1s)))))\n",
);

#[test]
fn a_reaction_over_the_span_read_with_a_group_round_gives_the_bound_of_its_departure() {
    // Read at the efforts' steps with the disc at rest, round, within its
    // bound: m e w^2 on the force.
    let w = Scratch::new("dyn-reaction-disc");
    let f = w.model("disc.ink", DISC);
    let d = w.ink(["effort", &f]).ok();
    let r = by(&d["reactions"], "path", "root/s/disc");
    let m = 7.85e-3 * PI * 25.0 * 10.0 * 30.0 * (2.0 * PI * 100.0f64).powi(2) * 1e-6; // N: the hole's m e w^2
    assert!(close(n(&r["bound"]["force"]), m, 0.02) && n(&r["bound"]["moment"]) > 0.0 && r.get("units").is_none(), "want {m}: {}", pretty(r));
}

#[test]
fn a_reaction_or_load_lacking_what_it_rests_on_is_refused_naming_it() {
    let w = Scratch::new("dyn-reaction-lacking");
    let base = "(body base (shape (move (box 120 80 20) 0 0 10)) (attr density 2.7g/cm3))";
    let bare = w.model("bare.ink", &arm(4.0).replace(base, "(body base (shape (move (box 120 80 20) 0 0 10)))"));
    let c = claim(&w.ink(["require", &bare]).ok(), "stands").clone();
    assert!(verdict(&c) == "unknown" && s(&c["unavailable"]).contains("root/a/base") && s(&c["unavailable"]).contains("no mass"), "massless base: {}", pretty(&c));
    let refusal = |name: &str, text: String| w.ink(["check", &w.model(name, &text)]).err()["error"].clone();
    let e = refusal("nowhere.ink", arm(4.0).replace("root/a/link (force", "root/a/lnk (force"));
    assert!(s(&e["message"]).contains("no instance root/a/lnk") && s(&e["hint"]).contains("root/a/link"), "a load placed nowhere: {}", pretty(&e));
    let e = refusal("weightless.ink", arm(4.0).replace("(gravity 0 0 -9.81m/s2)", ""));
    assert!(s(&e["message"]).contains("gravity") && s(&e["hint"]).contains("(gravity 0 0 -9.81m/s2)"), "without gravity: {}", pretty(&e));
    let e = refusal("forceless.ink", arm(4.0).replace("(force 0 0 -0.981N)", ""));
    assert!(s(&e["message"]).contains("force") && s(&e["hint"]).contains("(load NAME PATH"), "a load of nothing: {}", pretty(&e));
}

/// A joint is an actuator by `(actuator)` or `(actuator LO HI)` as by `(effort-limits LO HI)`: its effort is
/// worked out the same, and with LO HI the solver knows what it can give (its margin), but nothing is claimed
/// of it, and a design that asks nothing is asked nothing. (A catalog motor whose torque is a rating of the
/// part writes it: the rating is the one thing asked.)
#[test]
fn an_actuator_stating_no_limit_is_worked_out_and_claims_nothing() {
    let w = Scratch::new("dyn-actuator");
    let limited = w.ink(["effort", &w.model("limited.ink", EFFORT)]).ok();
    let f = w.model("actuator.ink", &EFFORT.replace("(effort-limits -150N*mm 150N*mm)", "(actuator)"));
    let d = w.ink(["effort", &f]).ok();
    let (a, b) = (&d["actuators"][0], &limited["actuators"][0]);
    // The steps that find a limit's crossing are the limited one's alone: the two read the peak at their own moments.
    assert!(a["path"] == "root/p/bar" && close(n(&a["peak"]["effort"]), n(&b["peak"]["effort"]), 1e-3) && close(n(&a["rms"]), n(&b["rms"]), 1e-2), "{}\n{}", pretty(a), pretty(b));
    assert!(n(&b["margin"]["margin"]) < 0.0 && a["margin"].is_null(), "{}\n{}", pretty(a), pretty(b));
    let asked = w.ink(["require", &f]).ok();
    assert!(asked["requirements"].as_array().unwrap().is_empty() && s(&asked["nothing_required"]).starts_with("nothing was asked"), "{}", pretty(&asked));
    let check = w.ink(["check", &w.model("still.ink", &EFFORT.replace("(effort-limits -150N*mm 150N*mm)", "(actuator)").replace("(gravity 0 0 -9.81m/s2)\n", "")), "--load"]).ok();
    assert!(check["not_checked"].is_null() || strs(&check["not_checked"]).is_empty(), "{}", pretty(&check));
    // With what it can give: the limited one's margin and steps, to the digit, and still no claim.
    let g = w.model("gives.ink", &EFFORT.replace("(effort-limits -150N*mm 150N*mm)", "(actuator -150N*mm 150N*mm)"));
    let d = w.ink(["effort", &g]).ok();
    assert!(d["actuators"][0]["margin"] == b["margin"] && d["actuators"][0]["peak"] == b["peak"], "{}\n{}", pretty(&d["actuators"][0]), pretty(b));
    assert!(w.ink(["require", &g]).ok()["requirements"].as_array().unwrap().is_empty());
    let e = w.ink(["check", &w.model("bad.ink", &EFFORT.replace("(effort-limits -150N*mm 150N*mm)", "(actuator 150N*mm)"))]).err()["error"].clone();
    assert!(s(&e["message"]).starts_with("usage: (actuator [LO HI])"), "{}", pretty(&e));
}

/// A load's name is its own, its point may be written as an instance's place is, and `ink check` lists what
/// the design applies: the load stated twice (by a pasted line) was summed, doubling a motor's torque with
/// nothing said.
#[test]
fn a_load_is_stated_once_and_listed() {
    let w = Scratch::new("dyn-load-once");
    let at_tip = "(load payload root/a/link (force 0 0 -0.981N) (at (frame-of root/a/link tip)))";
    let twice = w.model("twice.ink", &arm(4.0).replace(at_tip, &format!("{at_tip}\n{at_tip}")));
    let e = w.ink(["check", &twice]).err()["error"].clone();
    assert!(s(&e["message"]).contains("the load 'payload' is already stated, at line 12") && s(&e["hint"]).contains("one load, one name"), "{}", pretty(&e));
    // A still point, as the reference writes it: the same torque as the frame it coincides with at t = 0.
    let effort = |name: &str, text: String| n(&w.ink(["effort", &w.model(name, &text), "--t", "0"]).ok()["actuators"][0]["effort"]);
    let (moving, still) = (effort("tip.ink", arm(4.0)), effort("still.ink", arm(4.0).replace("(at (frame-of root/a/link tip))", "(at (move 150 0 100))")));
    assert!(close(moving, still, 1e-9) && moving.abs() > 100.0, "{moving} {still}");
    let d = w.ink(["check", &w.model("listed.ink", &arm(4.0)), "--load"]).ok();
    assert!(same(&d["gravity"], serde_json::json!([0, 0, -9.81])) && d["supports"].as_array().unwrap().is_empty(), "{}", pretty(&d));
    assert!(
        same(&d["loads"], serde_json::json!([{"name": "payload", "on": "root/a/link", "force_N": [0, 0, -0.981], "torque_Nmm": [0, 0, 0], "line": 12}])),
        "{}",
        pretty(&d["loads"])
    );
}

/// What holds a joint's load is said by a form, `(supports PATH S ...)`: at the top level with paths from root,
/// or in an assembly with its own instances'. Either is the division `(reaction PATH (by S ...))` reads, listed
/// by `ink check`; and a pose held at an assembly's default arguments is a pose, as one written `(arm 0)` is.
#[test]
fn supports_are_named_by_a_form_and_a_default_pose_is_a_pose() {
    let w = Scratch::new("dyn-supports");
    let slider = include_str!("../../../examples/slider-loads.ink");
    std::fs::write(w.join("slider-loads.json"), include_str!("../../../examples/slider-loads.json")).unwrap();
    let named = serde_json::json!([{"reaction": "root/s/c", "by": ["root/s/c"]}, {"reaction": "root/s/pan", "by": ["root/s/ba", "root/s/bb", "root/s/th"]}]);
    let listed = |name: &str, text: &str| w.ink(["check", &w.model(name, text), "--load"]).ok()["supports"].clone();
    assert!(same(&listed("as-written.ink", slider), named.clone()), "{}", pretty(&listed("as-written.ink", slider)));
    // The pan's supports said at the top level, or via (reaction ... (by ...)), instead of inside slider.
    let unread = slider.replace("\n  (supports pan ba bb th))", ")");
    assert_ne!(unread, slider);
    assert!(!unread.contains("(supports "), "{unread}");
    let top = format!("{unread}(supports root/s/pan root/s/ba root/s/bb root/s/th)\n");
    assert!(same(&listed("top.ink", &top), named.clone()), "{}", pretty(&listed("top.ink", &top)));
    let by_expr = format!("{unread}(def pan (reaction root/s/pan (by root/s/ba root/s/bb root/s/th)))\n");
    assert!(same(&listed("by-expr.ink", &by_expr), named), "{}", pretty(&listed("by-expr.ink", &by_expr)));
    let e = w.ink(["check", &w.model("alone.ink", &format!("{unread}(supports root/s/pan)\n"))]).err()["error"].clone();
    assert!(s(&e["message"]).contains("supports names the parts that hold it") && s(&e["hint"]).starts_with("(supports PATH S ...)"), "{}", pretty(&e));

    // EFFORT's pendulum takes its angle as an argument with a default: placed by name alone it is held there.
    let held = w.model("held.ink", &EFFORT.replace("(instance p (pendulum (* 60 (sin (* 180deg/s t)))))", "(instance p pendulum)"));
    let (d, written) = (
        w.ink(["effort", &held, "--t", "0"]).ok(),
        w.ink(["effort", &w.model("zero.ink", &EFFORT.replace("(pendulum (* 60 (sin (* 180deg/s t))))", "(pendulum 0)")), "--t", "0"]).ok(),
    );
    assert!(
        d["actuators"][0]["path"] == "root/p/bar" && d["actuators"][0]["effort"] == written["actuators"][0]["effort"] && d["reactions"] == written["reactions"],
        "{}",
        pretty(&d)
    );
}

/// A reaction's bound is written as computed (ink_lang::command::n): to 10
/// digits it read 72.9885294 N for 72.98852940421851 N, below the bound it
/// bounds.
#[test]
fn a_reactions_bound_is_written_exactly() {
    let w = Scratch::new("dyn-reaction-bound");
    let f = w.model(
        "reaction.ink",
        "(import \"std\")\n(gravity 0 0 -9.81m/s2)\n(t-span 0 2s)\n(body disc (shape (subtract (cylinder-z 50 10) (move (cylinder-z 5 20) 30 0 0))) (attr density 7.85g/cm3))\n\
(assembly (spinner (a 0)) (instance disc disc (at (rot-x 90)) (turn z a) (effort-limits -1e4N*mm 1e4N*mm)))\n(assembly root (instance s (spinner (* 360 100 (/ t 1s)))))",
    );
    let d = w.ink(["effort", &f]).ok();
    let bounds: Vec<f64> = d["reactions"].as_array().unwrap().iter().flat_map(|r| [&r["bound"]["force"], &r["bound"]["moment"]]).filter_map(Value::as_f64).collect();
    let rounded = |v: f64| format!("{v:.9e}").parse::<f64>().unwrap();
    assert!(!bounds.is_empty() && bounds.iter().all(|&b| b != rounded(b)), "bounds written to 10 digits: {bounds:?}");
}

#[test]
fn unchecked_effort_limits_prevent_a_guarantee() {
    let w = Scratch::new("unchecked-effort");
    let motor = w.model(
        "motor.ink",
        "(gravity 0 0 0) (body block (shape (max (- (abs x) 1) (- (abs y) 1) (- (abs z) 1))) (attr density 1g/cm3))
        (assembly (axis (position 0)) (instance block block (at) (slide x position) (effort-limits 1N 2N)))
        (assembly root (instance a (axis 0)))",
    );
    let limit = "effort-limits:root/a/block";
    assert_eq!(claim(&w.ink(["require", &motor]).ok(), limit)["status"], "fails");
    let parent = w.model("parent.ink", r#"(import "motor.ink") (require control (value 1 > 0))"#);
    for args in [vec![], vec!["--name", limit]] {
        let run = w.ink([&["require", &parent][..], &args].concat()).run();
        assert_eq!(run.code, Some(3));
        let d = run.json();
        assert_eq!(strs(&d["not_checked"]), [limit]);
        assert_eq!(claim(&d, limit)["status"], "unknown");
        assert!(s(&claim(&d, limit)["declared_at"]).ends_with("motor.ink:2"));
        assert_eq!(d["all_yes"], false);
        assert_eq!(d["all_proven"], false);
    }
    let check = w.ink(["check", &parent, "--load"]).ok();
    let d = w.ink(["require", &parent]).ok();
    assert_eq!(strs(&d["not_checked"]), [limit]);
    assert_eq!(d["not_checked"], check["not_checked"]);
    assert_eq!(d["not_checked_note"], check["not_checked_note"]);
    assert_eq!(d["all_yes"], false);
    assert_eq!(d["all_proven"], false);
    assert_eq!(claim(&d, "control")["status"], "holds");
    let only = w.ink(["require", &parent, "--name", "control"]).ok();
    assert_eq!(only["all_proven"], true);
    let missing = w.ink(["require", &parent, "--name", limit]).ok();
    assert_eq!(strs(&missing["not_checked"]), [limit]);
    assert_eq!(missing["all_yes"], false);
    assert!(missing.get("no_claims").is_none());
    assert_eq!(claim(&missing, limit)["status"], "unknown");
    assert_eq!(missing["counts"], serde_json::json!({"not settled": 1}));
    let collision = w.model("duplicate.ink", &format!(r#"(import "motor.ink") (require {limit} (value 1 > 0))"#));
    assert!(s(&w.ink(["require", &collision]).err()["error"]["message"]).contains(limit));
}

#[test]
fn own_effort_limits_include_a_motor_placed_without_arguments() {
    let w = Scratch::new("imported-effort");
    w.model(
        "motor.ink",
        r#"(import "std")
        (body weight (shape (box 100 100 100)) (attr density 1g/cm3))
        (assembly motor (instance shaft weight (at) (slide z 0) (effort-limits -2N 2N)))
        (require imported (value 0 > 1))"#,
    );
    let model = w.model(
        "caller.ink",
        r#"(import "motor.ink" as m) (gravity 0 0 -9.81m/s2)
        (param count 1 (cases 2)) (assembly root (for i (range count) (instance m-$i m.motor)))
        (require local (value 1 > 0))"#,
    );
    let inventory = w.ink(["check", &model, "--load"]).ok();
    assert_eq!(claim(&inventory, "effort-limits:root/m-0/shaft")["own"], true);
    let run = w.ink(["require", &model]).run();
    assert_eq!(run.code, Some(3));
    let d = run.json();
    let motor = claim(&d, "effort-limits:root/m-0/shaft");
    assert_eq!(motor["status"], "unknown");
    assert!(s(&motor["unavailable"]).contains("reads no coordinate"), "{motor}");
    let added = claim(&d, "effort-limits:root/m-1/shaft");
    assert_eq!(added["at_values"], false);
    assert_eq!(added["cases"][0]["status"], "unknown");
    assert_eq!(d["imported_not_run"], 1);
    assert_eq!(d["counts"], serde_json::json!({"yes": 1, "not settled": 2}));
    assert_eq!(d["all_yes"], false);
    assert_eq!(d["all_proven"], false);
}

#[test]
fn unchecked_effort_limits_at_a_declared_case_prevent_a_guarantee() {
    let w = Scratch::new("unchecked-case");
    let m = w.model(
        "conditional.ink",
        "(param count 0 (cases 1 2)) (param lower 2) (body block (shape 1))
        (assembly root (for i (range count) (instance j-$i block (at) (slide x 0) (optional (effort-limits (* 1N lower) 1N)))))
        (require control (value 1 > 0))",
    );
    assert!(w.ink(["check", &m, "--load"]).ok().get("not_checked").is_none());
    let d = w.ink(["require", &m, "--keys", "name,status"]).ok();
    assert_eq!(strs(&d["not_checked"]), ["effort-limits:root/j-0", "effort-limits:root/j-1"]);
    assert_eq!(d["all_proven"], false);
    assert_eq!(d["all_yes"], false);
    assert_eq!(w.ink(["require", &m, "--param", "count=0"]).ok()["all_proven"], true);
}

#[test]
fn unchecked_effort_limits_leave_a_previously_holding_row_unknown() {
    let w = Scratch::new("unchecked-gravity-case");
    let model = w.model(
        "main.ink",
        "(param env 1 (cases 0)) (for i (range env) (gravity 0 0 0))
        (body b (shape (max (- (abs x) 1) (- (abs y) 1) (- (abs z) 1))) (attr density 1g/cm3))
        (assembly (axis (position 0)) (instance j b (at) (slide x position) (effort-limits -1N 1N)))
        (assembly root (instance a (axis 0)))",
    );
    let d = w.ink(["require", &model]).ok();
    let row = claim(&d, "effort-limits:root/a/j");
    assert_eq!(row["status"], "unknown");
    assert_eq!(row["at_values_result"]["status"], "holds");
    assert_eq!(row["cases"][0]["status"], "unknown");
    assert!(s(&row["cases"][0]["not_checked"]).contains("gravity"));
    assert_eq!(d["counts"], serde_json::json!({"not settled": 1}));
    assert_eq!(d["all_yes"], false);
}

#[test]
fn a_simulated_coordinate_requires_density_valid_through_its_motion() {
    let w = Scratch::new("simulated-density");
    let design = |rho: &str| {
        w.model(
            "cart.ink",
            &format!(
                r#"(import "std") (gravity 0 0 0) (t-span 0 0.01s)
        (body cart (shape (box 10 10 10)) (attr density {rho}))
        (def travel (simulate 0))
        (assembly (axis (q 0)) (instance cart cart (at) (slide x q) (effort-limits -2N 2N)))
        (assembly root (instance a (axis travel))) (command root/a/cart 1)
        (require travel-limit (value travel < 20)) (require independent (value 1 < 0))"#
            ),
        )
    };
    let bad = w.ink(["require", &design("(* 1g/mm3 (select (- (/ t 1s) 0.005) 0.0027 2.7))")]).any();
    assert_eq!(bad["requirements"][0]["status"], "unknown");
    assert!(s(&bad["requirements"][0]["unavailable"]).contains("density"), "{bad}");
    assert_eq!(bad["requirements"][1]["status"], "fails");
    let good = w.ink(["require", &design("2.7g/cm3")]).any();
    assert_eq!(good["requirements"][0]["status"], "holds");
    assert!((n(&good["requirements"][0]["value"]) - 18.5185).abs() < 0.001);
}

#[test]
fn a_model_steps_the_whole_motion_so_an_estimates_condition_is_shown_over_all_of_it() {
    // A slide driven by an estimate that applies before 0.25 s only, over a
    // motion of 1 s: the model reads the drive at every step, whatever the
    // moment asked, so there is none. The whole-span RMS was a yes asked at
    // t = 0.1 and unknown at t = 0.5.
    let sc = Scratch::new("effort-estimate-span");
    let src = r#"(import "std")
(gravity 0 -9.81m/s2 0)
(t-span 0 1s)
(body w (shape (box 10 10 10)) (attr mass 2g))
(assembly (lift (h 0)) (instance s w (at) (slide y h) (effort-limits -1000N 1000N)))
(assembly root (instance l (lift (estimate-when (* 10mm/s2 t t) (< t LIMIT) "startup only"))))
(require rms (value (effort-rms root/l/s) < 100N))
"#;
    let model = sc.model("lift.ink", &src.replace("LIMIT", ".25s"));
    for t in ["0.1", "0.5"] {
        let d = sc.ink(["require", &model, "--t", t]).any();
        assert!(statuses(&d)["rms"] == "unknown" && pretty(&d).contains("startup only"), "t {t}: {}", pretty(&d));
    }
    let e = sc.ink(["effort", &model]).any();
    assert!(e["ok"] == false && pretty(&e).contains("startup only") && pretty(&e).contains("every moment of t 0 to 1"), "{}", pretty(&e));
    let model = sc.model("whole.ink", &src.replace("LIMIT", "2s"));
    let d = sc.ink(["require", &model, "--t", "0.1"]).any();
    assert_eq!(statuses(&d)["rms"], "holds estimate", "{}", pretty(&d));
}

#[test]
fn a_simulated_part_is_shown_clear_between_the_simulations_steps() {
    // A block let fall from rest, its drop simulated, passes a 0.2 mm wall
    // 600 mm down at 3.4 m/s: it is in the wall for 0.6 ms, between two of
    // the simulation's steps (16 ms apart), and every step's middle is
    // clear. Each step is asked with t ranged over it, the block moving as
    // the simulation says it does about that moment: a no. With the wall
    // out of its way every step is shown, and the estimate stands.
    let sc = Scratch::new("simulated-between");
    let src = r#"(import "std")
(gravity 0 -9.81m/s2 0)
(t-span 0 0.5s)
(body block (shape (box 2 2 2)) (attr density 2.7g/cm3))
(body wall (shape (box 20 0.2 20)))
(assembly (drop (h 0)) (instance b block (at) (slide y h)))
(assembly root (instance d (drop (simulate 0))) (instance w wall (at (move ASIDE -600 0))))
(require clear (disjoint (shape-of root/d/b) (shape-of root/w)))
"#;
    let d = sc.ink(["require", &sc.model("through.ink", &src.replace("ASIDE", "0"))]).any();
    assert_eq!(statuses(&d)["clear"], "fails estimate", "{}", pretty(&d));
    let t = d["requirements"][0]["t"].as_f64().unwrap_or(0.0);
    assert!((0.34..0.36).contains(&t), "{}", pretty(&d));
    let d = sc.ink(["require", &sc.model("aside.ink", &src.replace("ASIDE", "30"))]).any();
    assert_eq!(statuses(&d)["clear"], "holds estimate", "{}", pretty(&d));
}

#[test]
fn a_band_on_what_the_model_gives_is_carried_into_a_study_loaded_with_it() {
    // A 188 kg block spun at 100 rev/s about its middle, a box turned 17 degrees in its own frame: it
    // repeats every half turn about its spin (|x| and |y| of its own turned axes), so its centre is on
    // the axis exactly and its reaction is its weight within what its mass is known to. A post loaded
    // with it deflects 0.0006 mm under the weight, and that band is carried into the answer by the
    // post's response to each of the reaction's values alone: with room (under 1 mm) and under
    // 0.005 mm, the claim holds, an estimate. (Measured round within its measure instead, its centre
    // known to 0.8 mm, an imbalance under that would add 6e4 N at its spin, and the band up to 1.29
    // mm left both not settled.) No band is dropped for its size, and none leaves a study with no
    // answer. (On cells of 2.5:
    // on 5, its grids read 5.392e-4, 5.536e-4 and 5.541e-4 mm at the post's top, which the grids' rule
    // takes as moving apart; its coarsest grid settled only while it read the top a node short, 5.017e-4.)
    let sc = Scratch::new("study-carries-band");
    let src = r#"(import "std")
(gravity 0 0 -9.81m/s2)
(t-span 0 1s)
(body wheel (shape (rotate-z (box 400 300 200) 17)) (attr density 7.85g/cm3))
(assembly (spinner (a 0)) (instance w wheel (at (rot-x 90)) (turn z a) (effort-limits -1e9N*mm 1e9N*mm)))
(assembly root (instance s (spinner (* 360 100 (/ t 1s)))))
(def r (reaction root/s/w (at)))
(study post (shape (box 40 40 100)) (elastic (material steel)) (fix base (- z -48))
  (load top (force (nth r 0) (nth r 1) (neg (nth r 2))) (region (- 48 z))) (res 2.5))
(require room (max-displacement post < 1mm))
(require tight (max-displacement post < 0.005mm))
"#;
    let d = sc.ink(["require", &sc.model("spin.ink", src)]).any();
    let got = statuses(&d);
    assert!(got["room"] == "holds estimate" && got["tight"] == "holds estimate", "{}", pretty(&d));
    let tight = d["requirements"].as_array().unwrap().iter().find(|r| r["name"] == "tight").unwrap();
    assert!(pretty(tight).contains("reaction:root/s/w.fz within"), "{}", pretty(tight));
}

#[test]
fn a_failing_effort_limits_row_agrees_with_effort_on_the_worst_margin() {
    let sc = Scratch::new("dyn-effort-worst-witness");
    // Box pendulum (exact measured mass) and cylinder pendulum (banded measured mass):
    // both fail well before their peak at t = 0.5 s, and `ink require` must report the
    // worst margin over the span (matching `ink effort`), not the first failing moment.
    for (name, shape) in [("box", "(move (box 20 300 10) 0 150 0)"), ("cyl", "(move (rotate-x (cylinder-z 8 300) 90) 0 150 0)")] {
        let f = sc.model(
            &format!("pendulum-{name}.ink"),
            &format!(
                r#"(import "std")
(gravity 0 0 -9.81m/s2)
(t-span 0 2s)
(body bar (shape {shape}) (attr density 2.7g/cm3))
(assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -100N*mm 100N*mm)))
(assembly root (instance p (pendulum (* 60 (sin (* 180deg/s t))))))"#
            ),
        );
        let req = sc.ink(["require", &f]).any();
        let eff = sc.ink(["effort", &f]).ok();
        let row = claim(&req, "effort-limits:root/p/bar");
        let margin = &eff["actuators"][0]["margin"];
        assert_eq!(verdict(row), "fails estimate", "{name}: {}", pretty(row));
        let want = n(&margin["margin"]);
        let got = n(&row["value"]);
        let t_req = n(&row["t"]);
        let t_eff = n(&margin["t"]);
        assert!(
            close(got, want, 1e-3) && (t_req % 1.0 - 0.5).abs() < 0.02 && (t_eff % 1.0 - 0.5).abs() < 0.02,
            "{name}: require ({got} at t={t_req}) vs effort ({want} at t={t_eff}):\n{}\n{}",
            pretty(row),
            pretty(margin)
        );
    }
}
