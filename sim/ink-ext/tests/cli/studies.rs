//! Sweeps over ranged params and the call's time limit, studies (`ink
//! solve`), the queries an agent's loop reads (ray, mass, distance), claims
//! on numbers, and
//! `for`'s generated forms.

use crate::common::*;
use serde_json::{Value, json};
use std::collections::BTreeMap;

/// Numbers compared as numbers (0 is 0.0), the rest as JSON.
fn same(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::Number(x), Value::Number(y)) => x.as_f64() == y.as_f64(),
        (Value::Array(x), Value::Array(y)) => x.len() == y.len() && x.iter().zip(y).all(|(p, q)| same(p, q)),
        (Value::Object(x), Value::Object(y)) => x.len() == y.len() && x.iter().all(|(k, v)| y.get(k).is_some_and(|w| same(v, w))),
        _ => a == b,
    }
}

/// A number, NaN where it is none (so every comparison with it fails).
fn n(v: &Value) -> f64 {
    v.as_f64().unwrap_or(f64::NAN)
}

/// Whether a string contains `text`.
fn has(v: &Value, text: &str) -> bool {
    v.as_str().is_some_and(|s| s.contains(text))
}

/// Whether a string names `word` (a name the test gave: a param, a solve,
/// a part) as a word of its own. What the kernel says around it is the
/// kernel's, and its tests'.
fn names(v: &Value, word: &str) -> bool {
    v.as_str().is_some_and(|s| s.split(|c: char| !(c.is_alphanumeric() || "_-./:".contains(c))).any(|w| w == word))
}

/// Whether any string of a list contains `text`.
fn any_has(v: &Value, text: &str) -> bool {
    v.as_array().is_some_and(|a| a.iter().any(|w| has(w, text)))
}

/// A list's rows by a key of theirs (`name`, `target`).
fn by(v: &Value, key: &str) -> BTreeMap<String, Value> {
    v.as_array().unwrap_or_else(|| panic!("not a list: {}", pretty(v))).iter().map(|r| (r[key].as_str().unwrap_or("?").to_string(), r.clone())).collect()
}

/// Each target's mass, by target.
fn masses(d: &Value) -> BTreeMap<String, f64> {
    by(&d["targets"], "target").into_iter().map(|(k, v)| (k, n(&v["mass"]))).collect()
}

/// The call's time limit (INK_DEADLINE) `secs` from now; searches stop 3 s
/// before it.
fn deadline(secs: f64) -> String {
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs_f64();
    format!("{:.3}", now + secs)
}

/// Two calls made side by side: the longer one's time, not their sum.
fn both<A: Send, B: Send>(a: impl FnOnce() -> A + Send, b: impl FnOnce() -> B + Send) -> (A, B) {
    std::thread::scope(|s| {
        let b = s.spawn(b);
        let a = a();
        (a, b.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
    })
}

/// A call answering afresh (INK_NO_CACHE), its TMPDIR the scratch's.
fn fresh<S: AsRef<std::ffi::OsStr>>(sc: &Scratch, args: impl IntoIterator<Item = S>) -> Ink {
    sc.ink(args).env("INK_NO_CACHE", "1")
}

/// The first claim of an `ink require` answer.
fn first(d: &Value) -> Value {
    d["requirements"][0].clone()
}

// ---- ranged sweeps and the call's time limit ----

/// A cross turned through a ranged angle inside a ring `r` in (its arms'
/// tips at 29+). The ring's field has no slope bound (a floor term, as a
/// belt's select), so the direct clearance search decides.
const CROSS: &str = r#"(import "std")
(param a 0)
(param r 30)
(param m 0)
(t-span 0 1s)
(def arm (move (box 27 4 4) 15 0 0))
(body ring (shape (max (- r (len2 x y)) (- (len2 x y) 34) (- (abs z) 3) (- (floor x) 1000))))
(body cross (shape (min arm (rotate-z arm 90) (rotate-z arm 180) (rotate-z arm 270))))
(assembly root (instance k ring) (instance c cross (at) (turn z (+ a (* m (/ t 1s))))))
(require clear (clearance (shape-of root/k) (shape-of root/c) 1))
"#;

#[test]
fn a_clearance_over_a_wide_ranged_turn_is_proven_and_a_failure_gives_the_turn_there() {
    // A finest cell settled at the angle's middle but not over its span has
    // the angle halved where it blurs the cell more than the cell's size,
    // not space: proven (it was unknown, 12648 cells undecided at every
    // size). A ring 0.9 clear fails, the witness giving the angle there and
    // the row the ranged params.
    let sc = Scratch::new("st-cross");
    let cross = sc.model("cross.ink", CROSS);
    let (held, failed) =
        both(|| first(&fresh(&sc, ["require", &cross, "--param", "a=0:180"]).ok()), || first(&fresh(&sc, ["require", &cross, "--param", "a=0:180", "--param", "r=29.4"]).ok()));
    let a = n(&failed["params"]["a"]);
    assert!(
        verdict(&held) == "holds" && verdict(&failed) == "fails" && (0.0..=180.0).contains(&a) && same(&failed["param_spans"], &json!({"a": [0.0, 180.0]})),
        "{}\n{}",
        pretty(&held),
        pretty(&failed)
    );
}

#[test]
fn a_distance_stopped_at_the_time_limit_gives_the_bracket_it_has_saying_so() {
    let sc = Scratch::new("st-dist-limit");
    let cross = sc.model("cross.ink", CROSS);
    let d = fresh(&sc, ["distance", &cross, "--a", "root/k", "--b", "root/c", "--param", "a=30"]).env("INK_DEADLINE", deadline(-1.0)).ok();
    // The true distance is 30 - |(28.5, 2)| = 1.43. Stopped before its first
    // box, the search has no pair yet: no upper bound, a lower one still sound.
    let up = &d["distance"]["upper_witness"];
    assert!(d["stopped"] == json!(true) && n(&d["distance"]["lower_proven"]) <= 1.43 && (up.is_null() || n(up) >= 1.43), "{}", pretty(&d));
}

// ---- derivatives through solves ----

// ---- solve ----

#[test]
fn the_shelf_study_converges_in_equilibrium_under_its_weight_and_takes_probes() {
    // The self-weight is the material's volume times its density. Probes take
    // --keys and --compact; one outside the grid says so, apart.
    let shelf = fixture("shelf-study.ink");
    let d = ink(["solve", &shelf, "--res", "6"]);
    let w = n(&d["loads"][0]["total"][2]);
    let weight = n(&d["grid"]["material_volume_mm3"]) * 1270e-9 * 9.81;
    assert!(
        d["solver"]["converged"] == json!(true) && n(&d["equilibrium_imbalance"]) < 1e-6 && (w + weight).abs() < 1e-6,
        "self-weight {w} against {weight}:\n{}",
        pretty(&json!({"solver": d["solver"], "imbalance": d["equilibrium_imbalance"], "loads": d["loads"]}))
    );
    let t = ink(["solve", &shelf, "--res", "6", "--points", "[[1000,0,0],[-100,5,2]]", "--keys", "von_mises_mpa", "--compact"])["probes"].clone();
    assert!(
        t["columns"] == json!(["von_mises_mpa"])
            && t["rows"][0] == json!([null])
            && n(&t["rows"][1][0]) > 0.0
            && t["outside"][0][0] == json!(0)
            && t["outside"].as_array().map(Vec::len) == Some(1)
            && has(&t["outside"][0][1], "in no cell of the solve"),
        "{}",
        pretty(&t)
    );
}

#[test]
fn a_solve_is_an_estimate_and_says_which_stress_peak_leaves_out_the_supports() {
    // accuracy null, as a claim on a study; the peak away from the fixes
    // says what it excludes (the margin, at least 2 cells: 12 at res 6),
    // and the near-support peak its margin and that it bounds nothing.
    let d = ink(["solve", &fixture("shelf-study.ink"), "--res", "6"]);
    let (away, near) = (&d["results"]["max_von_mises_mpa"], &d["results"]["max_von_mises_near_supports_mpa"]);
    assert!(
        d.get("accuracy") == Some(&Value::Null)
            && has(&d["estimate"], "not a bound")
            && has(&away["excludes"], "12 mm of fix 'wall'")
            && n(&near["within_mm"]) == 12.0
            && has(&near["note"], "bounds nothing")
            && n(&near["value"]) > 0.0,
        "{}",
        pretty(&d["results"])
    );
}

/// A solve states its material, and an elastic study gives its Poisson's ratio: none is assumed (a claim
/// on a part held sideways answers by it, and would rest on a number its file never gave).
#[test]
fn a_solve_states_its_material_and_assumes_no_poissons_ratio() {
    let sc = Scratch::new("material-line");
    let study = "(import \"std\")\n(study s (shape (box 10 10 10)) (elastic (E 70000MPa) NU (yield 200MPa)) (fix f (+ x 4)) (load l (force 0 0 -100N) (region (- x 4))) (res 1))";
    let m = ink(["solve", &sc.model("m.ink", &study.replace("NU", "(nu 0.33)")), "--study", "s"])["material"].clone();
    assert!(n(&m["E_mpa"]) == 70000.0 && n(&m["nu"]) == 0.33 && n(&m["yield_mpa"]) == 200.0 && m.get("nu_assumed").is_none(), "{}", pretty(&m));
    message_has(&sc.ink(["solve", &sc.model("none.ink", &study.replace("NU ", "")), "--study", "s"]).err()["error"]["message"], &["elastic needs (nu ratio)"]);
    // A material brings its own.
    let named = ink(["solve", &sc.model("named.ink", &study.replace("(E 70000MPa) NU (yield 200MPa)", "(material aluminum-6061)")), "--study", "s"])["material"].clone();
    assert!(n(&named["nu"]) > 0.3 && n(&named["nu"]) < 0.36, "{}", pretty(&named));
}

/// A 10 mm slab (x 0..10) whose face x < 0 is held at 100 °C from 20 °C
/// at t = 0, its far face insulated: alpha 0.1211 mm²/s, Fourier number
/// 0.727 at 600 s.
const HELD_SLAB: &str = "(import \"std\")\n(param r 0.5)\n(study heat (shape (move (box 11 1 1) 4.5 0 0)) (heat (k 0.2W/m/K) (density 1270kg/m3) (capacity 1300J/kg/K))\n  (temperature face (- x 0) 100degC) (transient (duration 600s) (steps 60) (initial 20degC)) (res r))\n";

#[test]
fn a_transient_with_a_held_face_converges_to_the_series_solution_as_the_cells_shrink() {
    // The series solution reads T(10) = 83.051 at 600 s; backward Euler's
    // 60 steps read it about 0.7% of the rise low, at every cell size. Each
    // step's residual was relative to a right-hand side the held face's
    // penalty dominates (as 1/h² against the heat capacity's), so most
    // steps took no iterations: 0.25 mm cells read 72.9, 0.125 mm 25.8.
    let sc = Scratch::new("st-held-slab");
    let slab = sc.model("slab.ink", HELD_SLAB);
    let series = 83.051;
    for r in ["0.5", "0.25", "0.125"] {
        let d = sc.ink(["solve", &slab, "--param", &format!("r={r}"), "--points", "[[10,0,0]]"]).ok();
        let t = n(&records(&d["probes"])[0]["temperature"]);
        assert!(d["solver"]["converged"] == json!(true) && (t - series).abs() / (series - 20.0) < 0.01, "at {r} mm, T(10) {t} against {series}: {}", pretty(&d["solver"]));
    }
}

#[test]
fn a_held_temperature_is_held_as_it_is_at_each_step() {
    // The slab held at 20 °C, its initial temperature, for 120 s, then at
    // 100 °C: at 720 s it reads as HELD_SLAB at 600 s, step for step (its
    // steps end at the switch, the one ending there holding 20). Read once,
    // at its start, a held temperature held 20 throughout; read at the step
    // ending at the switch, its heat came a step early.
    let sc = Scratch::new("st-held-in-time");
    let at = |src: &str| {
        let d = sc.ink(["solve", &sc.model("slab.ink", src), "--points", "[[10,0,0],[1,-0.4,0],[1,0.4,0]]"]).ok();
        assert!(d["solver"]["converged"] == json!(true), "{}", pretty(&d["solver"]));
        records(&d["probes"]).iter().map(|r| n(&r["temperature"])).collect::<Vec<f64>>()
    };
    let switched = |held: &str| HELD_SLAB.replace("(duration 600s) (steps 60)", "(duration 720s) (steps 72)").replace("100degC", held);
    let whole = at(HELD_SLAB)[0];
    let later = at(&switched("(+ 20degC (* 80 (- 1degC 0degC) (select (- t 120s) 0 1)))"))[0];
    assert!((later - whole).abs() < 1e-3 * (whole - 20.0) && (whole - 83.051).abs() / 63.051 < 0.01, "{later} at 720 s against {whole} at 600 s");
    // So too a held temperature varying over the face as well as in time,
    // against the same held from the start.
    let over = "(* 80 (- 1degC 0degC) (+ 1 y))";
    let whole = at(&HELD_SLAB.replace("100degC", &format!("(+ 20degC {over})")));
    let later = at(&switched(&format!("(+ 20degC (* {over} (select (- t 120s) 0 1)))")));
    assert!(whole[2] > whole[1] + 1.0, "{whole:?}");
    for (a, b) in whole.iter().zip(&later) {
        assert!((a - b).abs() < 1e-3 * (a - 20.0), "{later:?} at 720 s against {whole:?} at 600 s");
    }
    // A switch between the steps' ends ends a step: one at 30 s on steps of
    // 25 s, the step ending there holding the value before it.
    let d = sc
        .ink([
            "solve",
            &sc.model("off.ink", &switched("(+ 20degC (* 80 (- 1degC 0degC) (select (- t 30s) 0 1)))").replace("(duration 720s) (steps 72)", "(duration 100s) (steps 4)")),
        ])
        .ok();
    let history: Vec<(f64, f64)> = d["results"]["transient"]["history"].as_array().unwrap().iter().map(|p| (n(&p[0]), n(&p[1]))).collect();
    let at = history.iter().position(|p| (p.0 - 30.0).abs() < 1e-6).unwrap_or_else(|| panic!("{history:?}"));
    assert!((history[at].1 - 20.0).abs() < 1e-6 && history[at + 1].1 > 50.0, "{history:?}");
    // A steady study has no time to hold one at: refused alike by check
    // and solve.
    let steady =
        sc.model("steady.ink", &HELD_SLAB.replace("(transient (duration 600s) (steps 60) (initial 20degC))", "").replace("100degC", "(+ 20degC (* (- 1degC 0degC) (/ t 1s)))"));
    for cmd in ["check", "solve"] {
        let d = sc.ink([cmd, &steady]).any();
        assert!(d["ok"] == false && pretty(&d).contains("only a transient study steps through time"), "{cmd}: {}", pretty(&d));
    }
}

#[test]
fn a_held_temperatures_switch_is_read_on_its_side() {
    let sc = Scratch::new("st-held-switch");
    let slab = |shape: &str, at: &str, held: &str, transient: &str| {
        let src = format!(
            "(import \"std\")\n(study heat (shape {shape}) (heat (k 0.2W/m/K) (density 1270kg/m3) (capacity 1300J/kg/K))\n  (temperature face {at} {held}) (transient {transient} (initial 20degC)) (res 0.5))\n"
        );
        let d = sc.ink(["solve", &sc.model("slab.ink", &src), "--points", "[[0.5,0,0],[19.75,0,0]]"]).ok();
        (records(&d["probes"]).iter().map(|r| n(&r["temperature"])).collect::<Vec<f64>>(), n(&d["results"]["transient"]["peak"]["value"]))
    };
    let bar = "(move (box 11 1 1) 4.5 0 0)";
    let switch = |at: &str| format!("(+ 20degC (* 80 (- 1degC 0degC) (select (- t {at}) 0 1)))");
    // Switching to 100 C at the end of its 10 s, its last step holds 20 C
    // (read at the end, it held 100 over the step: 53.19 C at 0.5 mm).
    let (t, peak) = slab(bar, "(- x 0)", &switch("10s"), "(duration 10s) (steps 4)");
    assert!((t[0] - 20.0).abs() < 1e-6 && (peak - 20.0).abs() < 1e-6, "{t:?} {peak}");
    // A pulse of 1e-10 s at 2 s in a run of 1e6 s: its step holds 100 C
    // (the moment before its end, stepped back by a share of the run,
    // was before its start: it held nothing).
    let pulse = "(+ 20degC (* 80 (- 1degC 0degC) (- (select (- t 2s) 0 1) (select (- t 2.0000000001s) 0 1))))";
    let (_, peak) = slab(bar, "(- x 0)", pulse, "(duration 1000000s) (steps 100)");
    assert!((peak - 100.0).abs() < 1e-6, "{peak}");
    // Its switches are where its fix holds: at the far end of a 20 mm bar,
    // switching at t = x s/mm, none in its 10 s (over the whole bar, one at
    // every moment, refused).
    let (t, _) = slab("(move (box 20 1 1) 10 0 0)", "(- 19.5 x)", "(+ 20degC (* 80 (- 1degC 0degC) (select (- (/ t 1s) x) 0 1)))", "(duration 10s) (steps 4)");
    assert!((t[1] - 20.0).abs() < 1e-6, "{t:?}");
}

fn ranged_heater(power: &str) -> String {
    format!(
        r#"(import "std")
(param gain 1 (tol -0.9 99))
(study ramp (shape (box 8 8 8))
  (heat (k 10W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))
  (convection air -1 (coefficient 1W/m2/K) (ambient 0))
  (source input (power {power}) (region -1))
  (transient (duration 1s) (steps 4) (initial 0)) (res 1))"#
    )
}

fn check_ranged_heater(power: &str) {
    let sc = Scratch::new("ranged-heater");
    let heater = sc.model("heater.ink", &format!("{}\n(require below-ten (max-temperature ramp < 10degC))", ranged_heater(power)));
    // Both heaters supply 50 J at the hot endpoint to a cube of heat
    // capacity 0.512 J/K, with less than 0.04 J lost to convection.
    let hot = sc.ink(["require", &heater, "--param", "gain=100:100"]).any();
    let a = &hot["requirements"][0];
    assert!(a["status"] == "fails" && (97.0..99.0).contains(&n(&a["value"])), "{}", pretty(a));
    let ranged = sc.ink(["require", &heater]).any();
    let a = &ranged["requirements"][0];
    assert!(
        a["status"] == "unknown" && names(&a["unavailable"], "gain") && a["param_spans"]["gain"] == ranged["param_spans"]["gain"] && a.get("one_moment").is_none(),
        "{}",
        pretty(&ranged)
    );
}

#[test]
fn a_heat_sources_speed_reads_its_own_clock() {
    let sc = Scratch::new("heater-speed-clock");
    let src = r#"(import "std")
(t-span 0 2s)
(def travel (select (- (/ t 1s) 1.5) 0 10))
(def demand (max 1 (diff travel t)))
(study ramp (shape (box 8 8 8))
  (heat (k 10W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))
  (convection air -1 (coefficient 1W/m2/K) (ambient 0))
  (source input (power (* 1W demand)) (region -1))
  (transient (duration 1s) (steps 4) (initial 0)) (res 1))
(require cool (max-temperature ramp < 10degC))
(require speed (value demand < 10))
"#;
    let model = sc.model("heater.ink", src);
    let d = sc.ink(["require", &model]).any();
    assert_eq!(statuses(&d)["cool"], "holds estimate", "{}", pretty(&d));
    assert!(d["requirements"][0].get("motion").is_none());
    assert_eq!(d["requirements"][1]["status"], "unknown", "{}", pretty(&d));
    let model = sc.model("earlier.ink", &src.replace("1.5", "0.5"));
    let d = sc.ink(["require", &model]).any();
    assert_eq!(d["requirements"][0]["status"], "unknown", "{}", pretty(&d));
    // (the jump's clock: the kernel's words for it are its, and its tests')
    assert!(d["requirements"][0]["unavailable"].as_str().is_some_and(|u| !u.is_empty()), "{}", pretty(&d));
}

#[test]
fn ranged_heater_ramping_from_zero() {
    check_ranged_heater("(* 1W (/ t 1s) gain)");
}

#[test]
fn ranged_heater_switching_on_after_a_delay() {
    // The power is off at the start, middle and end of the study.
    check_ranged_heater("(* 1W (select (- (/ t 1s) 0.125) 0 (select (- (/ t 1s) 0.375) (* 2 gain) 0)))");
}

#[test]
fn a_moving_probe_of_a_ranged_heater_cannot_use_its_nominal_temperature_field() {
    let sc = Scratch::new("ranged-heater-probe");
    let heater = sc.model(
        "probe.ink",
        &format!("{}\n(t-span 0 1s)\n(def temp (study-field ramp temperature))\n(require cool (value (warp temp t 0 0) < 10))", ranged_heater("(* 1W (/ t 1s) gain)")),
    );
    let error = sc.ink(["require", &heater]).err();
    assert!(names(&error["error"]["message"], "gain"), "{}", pretty(&error));
}

/// A 20 x 10 x 2 mm block of heat capacity 0.4 J/K, barely cooled, whose
/// heater gives a 10 J triangular pulse from t = 2 to 3 s, inside the first
/// 6 s step of the default 100 over 600 s: its peak is at least 44.97 °C
/// (25 K from the pulse, convection taking at most 0.013 J of it).
fn heat_pulse(limit: &str, transient: &str) -> String {
    format!(
        "(import \"std\")\n(def part (box 20 10 2))\n(study hot (shape part) (heat (k 200W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))\n \
         (convection air -1 (coefficient 1W/m2/K) (ambient 20degC))\n (source chip (power (* 1W (* 20 (max 0 (- 1 (* 2 (abs (- (/ t 1s) 2.5)))))))) (region -1))\n \
         {transient} (res 1))\n(require cool (max-temperature hot < {limit}degC))\n(require hot-enough (max-temperature hot > {limit}degC))\n"
    )
}

/// A 0.1 s pulse of 1 J into a 4 mm end of a 100 x 10 x 2 mm conductive bar
/// (200 W/m/K), over 10 s: every step of 1 s, 2 s or 3.3 s spreads it over
/// more of the bar than it reached by its peak, so all three read the
/// end's peak low (23.08, 22.26, 21.79 degC; 28.35 on 640 steps), their
/// differences growing as the steps shrink. Read as three step counts not
/// converging, its error in time was 3 times the last difference and
/// `< 27` held; no error can be read from them, so it is unknown, saying
/// to take more steps, and on steps short enough to follow the pulse it
/// fails.
#[test]
fn a_peak_its_step_counts_do_not_converge_on_is_not_decided() {
    let sc = Scratch::new("st-heat-spread");
    let bar = |steps: usize| {
        sc.model(
            &format!("bar-{steps}.ink"),
            &format!(
                "(import \"std\")\n(def part (move (box 100 10 2) 50 0 0))\n(study hot (shape part) (heat (k 200W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))\n \
                 (convection air -1 (coefficient 1W/m2/K) (ambient 20degC))\n \
                 (source chip (power (* 1W (* 20 (max 0 (- 1 (* 20 (abs (- (/ t 1s) 0.05)))))))) (region (- (abs (- x 2)) 2)))\n \
                 (transient (duration 10s) (steps {steps}) (initial 20degC)) (res 1))\n(require cool (max-temperature hot < 27degC))\n"
            ),
        )
    };
    let d = sc.ink(["require", &bar(10)]).any();
    let a = &d["requirements"][0];
    assert!(
        a["status"] == "unknown" && has(&a["hint"], "(steps n)") && has(&a["band"]["basis"], "no error in time can be read from them (three step counts moving apart)"),
        "{}",
        pretty(a)
    );
    let d = sc.ink(["require", &bar(640)]).any();
    assert_eq!(statuses(&d)["cool"], "fails estimate", "{}", pretty(&d));
}

#[test]
fn a_heater_pulse_inside_a_time_step_puts_its_heat_in_and_its_steps_are_checked() {
    // Each step read its source at its end: the pulse put no heat in, and
    // `< 35` held at 20.0 °C (and `> 35` failed), told to add a t-span
    // (the study's own t is not the model's, and its peak is over time).
    let sc = Scratch::new("st-heat-pulse");
    let one = |name: &str, text: &str| sc.ink(["require", &sc.model(name, text)]).any();
    let d = one("pulse.ink", &heat_pulse("35", "(transient (duration 600s) (initial 20degC))"));
    let st = statuses(&d);
    assert!(st["cool"] == "fails estimate" && st["hot-enough"] == "holds estimate", "{}", pretty(&d));
    let a = &d["requirements"][0];
    let (x, band) = (n(&a["value"]), &a["band"]);
    assert!((44.5..45.1).contains(&x) && n(&band["low"]) < 44.97 && n(&band["high"]) > 44.97 && band["dominant"] == "time steps", "{}", pretty(a));
    assert!(!d.to_string().contains("one moment") && a.get("one_moment").is_none(), "{}", pretty(&d));
    // Under a t-span it is the same, asked once: the study reads no
    // moment of the model.
    let d = one("spanned.ink", &format!("(t-span 0 600s)\n{}", heat_pulse("35", "(transient (duration 600s) (initial 20degC))")));
    assert!(statuses(&d)["cool"] == "fails estimate" && d["requirements"][0].get("motion").is_none(), "{}", pretty(&d));
    // A bound inside its error in time is unknown, saying to take more
    // steps; one step has nothing to compare with, so decides nothing.
    for (limit, tr) in [("44.7", "(transient (duration 600s) (initial 20degC))"), ("35", "(transient (duration 600s) (steps 1) (initial 20degC))")] {
        let d = one("near.ink", &heat_pulse(limit, tr));
        let a = &d["requirements"][0];
        assert!(a["status"] == "unknown" && has(&a["hint"], "(steps n)") && has(&a["undecided"], "time steps"), "{limit}, {tr}: {}", pretty(a));
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_transient_peak_is_placed_where_and_when_it_is() {
    // Two separate cubes: the left heated until t = 1, the right after 9.
    // The peak, at 0.9 s or so, was placed at the final state's hottest
    // point, on the right cube, then still at 20 °C.
    let sc = Scratch::new("st-peak-place");
    let m = sc.model(
        "two.ink",
        "(import \"std\")\n(def part (min (move (box 2 2 2) -3 0 0) (move (box 2 2 2) 3 0 0)))\n(study hot (shape part) (heat (k 200W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))\n \
         (convection air -1 (coefficient 1000W/m2/K) (ambient 20degC))\n (source early (power (* 1W (select (- (/ t 1s) 1) 2 0))) (region x))\n (source late (power (* 1W (select (- 9 (/ t 1s)) 1 0))) (region (neg x)))\n \
         (transient (duration 10s) (steps 100) (initial 20degC)))\n(require cool (max-temperature hot < 90degC))\n",
    );
    let r = sc.ink(["solve", &m, "--study", "hot"]).ok()["results"].clone();
    let (mt, peak) = (&r["max_temperature"], &r["transient"]["peak"]);
    assert!(n(&mt["at"][0]) < -2.0 && n(&mt["time_s"]) <= 1.0 && mt["at"] == peak["at"] && n(&mt["value"]) > 90.0, "{}", pretty(&r));
    // Its three grids read 97.79, 97.71 and 97.60 on cells of 0.32, 0.16 and 0.08 mm, each change larger than
    // the last: no grid error, so the claim is not settled, read at the peak's place.
    let a = sc.ink(["require", &m]).any()["requirements"][0].clone();
    assert!(a["status"] == "unknown" && has(&a["estimate"], "three grids moving apart") && n(&a["at"][0]) < -2.0, "{}", pretty(&a));
}

#[test]
fn a_solve_that_glued_a_gap_shut_decides_no_claim_either_way() {
    // A hairline cut through the bar parts its loaded end from its fix:
    // there is no finite answer. Its cells glue the cut shut and read a
    // stiff, strong bar; the claims that would hold there were unknown,
    // and those that would fail failed.
    let sc = Scratch::new("st-glued-both");
    let m = sc.model(
        "cut.ink",
        "(import \"std\")\n(def part (subtract (move (box 20 8 8) 10 0 0) (move (box 0.000001 40 40) 5.01 0 0)))\n(study s (shape part) (elastic (material steel))\n \
         (fix root (- x 1))\n (load tip (force 0 0 -1N) (region (- 19 x))) (res 1))\n(require stiff (max-displacement s > 0.1))\n(require strong (min-safety-factor s < 2))\n",
    );
    let d = sc.ink(["require", &m]).any();
    for a in d["requirements"].as_array().unwrap() {
        assert!(a["status"] == "unknown" && has(&a["undecided"], "glue shut") && a["hint"].as_str().is_some_and(|h| !h.is_empty()), "{}", pretty(a));
    }
}

#[test]
fn whole_part_transient_keeps_the_insulated_heated_cube() {
    // Half of 10.24 W heats an insulated 0.512 J/K cube for one second.
    let sc = Scratch::new("st-insulated-cube");
    let cubes = sc.model(
        "cubes.ink",
        r#"(import "std")
(def cubes (min (move (box 8 8 8) -8 0 0) (move (box 8 8 8) 8 0 0)))
(study heat-up (shape cubes) (heat (k 10W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))
  (temperature cold x 0) (source both (power 10.24W) (region -1))
  (transient (duration 1s) (steps 4) (initial 0)) (res 1))
(require cool (max-temperature heat-up < 5degC))"#,
    );
    let a = first(&fresh(&sc, ["require", &cubes]).any());
    assert_eq!(a["status"], "fails", "{}", pretty(&a));
    assert!((n(&a["value"]) - 10.0).abs() < 1e-4 && n(&a["at"][0]) > 0.0 && n(&a["band"]["low"]) > 5.0, "{}", pretty(&a));
    let s = fresh(&sc, ["solve", &cubes]).ok();
    assert_eq!(s["grid"]["dropped_floating_cells"], 0);
    assert!((n(&s["results"]["transient"]["peak"]["value"]) - 10.0).abs() < 1e-4);
}

#[test]
fn whole_part_transient_needs_heat_capacity_on_each_insulated_piece() {
    let sc = Scratch::new("st-massless-cube");
    let cubes = sc.model(
        "cubes.ink",
        r#"(import "std")
(def cubes (min (move (box 8 8 8) -8 0 0) (move (box 8 8 8) 8 0 0)))
(study hot (shape cubes) (heat (k 10W/m/K) (density (select x 1000kg/m3 0)) (capacity 1000J/kg/K))
  (temperature cold x 0) (source both (power 10.24W) (region -1))
  (transient (duration 1s) (steps 4) (initial 0)) (res 1))
(require cool (max-temperature hot < 5degC))"#,
    );
    let r = fresh(&sc, ["require", &cubes]).any();
    let a = first(&r);
    assert!(r["all_yes"] == false && a["status"] == "unknown" && has(&a["unavailable"], "insulated piece has no positive heat capacity"), "{}", pretty(&r));
    let e = fresh(&sc, ["solve", &cubes]).err();
    assert!(has(&e["error"]["message"], "insulated piece has no positive heat capacity"), "{}", pretty(&e));
}

#[test]
fn whole_part_transient_distinguishes_a_gap_from_a_missed_web() {
    let sc = Scratch::new("st-web");
    let src = r#"(import "std")
(def blocks (min (move (box 4 4 4) -3 0 0) (move (box 4 4 4) 3 0 0)))
(def body (min blocks (move (box 3 1 0.001) 0 0 0.015625)))
(study s (shape body) (heat (k 100000W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))
  (temperature sink (+ x 2) 0) (source heater (power 0.05W) (region (- 1 x)))
  (transient (duration 100s) (steps 4) (initial 0)) (res 1))
(require hot (max-temperature s > 60degC))"#;
    let web = sc.model("web.ink", src);
    let a = first(&fresh(&sc, ["require", &web]).any());
    assert_eq!(a["status"], "unknown");
    assert!(has(&a["unavailable"], "web") && has(&a["unavailable"], "Refine"), "{}", pretty(&a));
    let film = src.replace("(source heater", "(convection film (max (- 1 x) (- 1.99 z)) (coefficient 0.1W/m2/K) (ambient 0)) (source heater");
    let steady = film.replace("(transient (duration 100s) (steps 4) (initial 0))", "");
    for (name, body) in [("film", film), ("steady", steady), ("temperature", src.replace("(source heater", "(temperature heater-end (- 4 x) 78degC) (source heater"))] {
        let f = sc.model(&format!("{name}.ink"), &body);
        let a = first(&fresh(&sc, ["require", &f]).any());
        assert_eq!(a["status"], "unknown", "{name}: {a}");
        assert!(has(&a["unavailable"], "web the cells missed"), "{name}: {a}");
    }
    let shapes = src.split("(study").next().unwrap();
    for (name, additional, claim) in [("elastic", "(load pull (force 1N 0 0) (region (- 3 x)))", "(max-displacement s < 1)"), ("modal", "(modes 1)", "(min-frequency s > 1Hz)")] {
        let study =
            format!("{shapes}(study s (shape body) (elastic (E 10000MPa) (nu 0) (density 1000kg/m3)) (fix ends (- 4 (abs x))) {additional} (res 1)) (require limit {claim})");
        let f = sc.model(&format!("{name}.ink"), &study);
        let a = first(&fresh(&sc, ["require", &f]).any());
        assert_eq!(a["status"], "unknown", "{name}: {a}");
        assert!(has(&a["unavailable"], "web the cells missed"), "{name}: {a}");
    }
    let split = sc.model("split.ink", &src.replace("(shape body)", "(shape blocks)").replace("(region (- 1 x))", "(region (- 3 x))"));
    let a = first(&fresh(&sc, ["require", &split]).any());
    assert_eq!(a["status"], "holds", "{}", pretty(&a));
    assert!(n(&a["value"]) >= 78.12);
    let sol = fresh(&sc, ["solve", &split]).ok();
    assert!(any_has(&sol["warnings"], "solved as insulated"));
}

#[test]
fn whole_part_clearance_around_a_heated_pin_is_real() {
    let sc = Scratch::new("st-pin-clearance");
    let f = sc.model(
        "pin.ink",
        r#"(import "std")
(def body (min (cylinder-z 2 6) (subtract (box 18 18 6) (cylinder-z 5 12))))
(study s (shape body) (heat (k 10W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))
 (temperature sink (- 7 x) 0) (source hot (power 0.1W) (region (- (len2 x y) 2.5)))
 (transient (duration 1s) (steps 2) (initial 0)) (res 0.5))
(require warm (max-temperature s > 0.5degC))"#,
    );
    // On cells of 2, 1 and 0.5 its three grids read 1.3320, 1.3270 and 1.3266,
    // converging on the exact 1.3263 faster than their elements can, so it
    // is given three times their whole movement, 0.017. (They read 1.3276,
    // 1.3275 and 1.3266, close and not in order, while the round pin's and
    // the bore's cells under 3% full were dropped: once moving apart, not
    // settled.) On 0.25 it holds too.
    let a = first(&fresh(&sc, ["require", &f]).any());
    assert!(a["status"] == "holds" && has(&a["estimate"], "three grids converging"), "{a}");
    assert!((n(&a["value"]) - 100.0 / (24.0 * std::f64::consts::PI)).abs() <= n(&a["grid_error"]), "{a}");
    let f = sc.model("pin-finer.ink", &std::fs::read_to_string(&f).unwrap().replace("(res 0.5)", "(res 0.25)"));
    let a = first(&fresh(&sc, ["require", &f]).any());
    assert_eq!(a["status"], "holds", "{a}");
    // 0.1 J heats a cylinder of volume 24 pi mm^3.
    assert!((n(&a["value"]) - 100.0 / (24.0 * std::f64::consts::PI)).abs() < 0.01, "{a}");
}

#[test]
fn whole_part_a_web_can_leave_the_pieces_box_and_come_back() {
    let sc = Scratch::new("st-detour");
    let f = sc.model(
        "detour.ink",
        r#"(import "std")
(def blocks (min (move (box 4 4 4) -3 0 0) (move (box 4 4 4) 3 0 0)))
(def wire (min (move (box 0.001 7 1) -3.015625 3 0) (move (box 0.001 7 1) 3.015625 3 0) (move (box 6.1 0.001 1) 0 6.015625 0)))
(def body (min blocks wire))
(study s (shape body) (heat (k 100000W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))
  (temperature sink (+ x 2) 0) (source heater (power 0.05W) (region (- 1 x)))
  (transient (duration 100s) (steps 4) (initial 0)) (res 1))
(require hot (max-temperature s > 60degC))"#,
    );
    let a = first(&fresh(&sc, ["require", &f]).any());
    assert_eq!(a["status"], "unknown", "{a}");
    assert!(has(&a["unavailable"], "web the cells missed"), "{a}");
}

#[test]
fn whole_part_two_reversed_mirrors_win_on_a_two_span_beam() {
    let sc = Scratch::new("st-two-spans");
    // A middle support forces the spans to bow in opposite directions.
    let f = sc.model(
        "two-spans.ink",
        r#"(import "std")
(def bar (move (box 24 4 2) 12 0 1))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (mirror middle x 0) (mirror thickness z 0)
  (fix guide -1 (dofs y) (support))
  (fix mid (- x 0.01) (dofs z) (support))
  (fix end (- 23.99 x) (dofs z) (support))
  (load squeeze (pressure 5MPa) (region (face bar (normal 1 0 0)))) (res 1))
(require stable (min-buckling-factor s > 10))"#,
    );
    let sol = fresh(&sc, ["solve", &f]).ok();
    let b = &sol["results"]["buckling"];
    assert_eq!(b["sectors_solved"], 4);
    assert_eq!(b["antisymmetric_axes"], json!(["x", "z"]));
    assert!(n(&b["factor"]) < 10.0 && n(&b["symmetric_factor"]) > 10.0, "{b}");
    // Its factor on three grids converges, 4.270, 4.283 and 4.317 on cells
    // of 0.5, 1 and 2 mm: it fails, decided. (Held by penalties over the cells
    // its supports cut, it rose 4.425, 4.459 and 4.595, by more than its
    // eigensolve's band: no grid error, and not settled.)
    let a = first(&fresh(&sc, ["require", &f]).any());
    assert!(a["status"] == "fails" && has(&a["estimate"], "three grids converging"), "{a}");
}

#[test]
fn pressure_uncertainty_on_a_soft_bar_cannot_be_diluted_by_a_stiff_bar() {
    let sc = Scratch::new("st-soft-ridge");
    let f = sc.model(
        "bars.ink",
        r#"(import "std")
(def part (min (move (box 16 4 4) 8 -4 0) (move (box 16 4 4) 8 4 0)))
(study s (shape part) (elastic (E (* 1MPa (select y 100 1000000))) (nu 0))
  (fix root (- x 1))
  (load ends (pressure (select y (* 1MPa (- -0.001 (* 0.12 (max 0 (- 1 (* 64 (abs (+ y 4.03125)))))))) -10MPa))
        (region (face part (normal 1 0 0))))
  (res 1))
(require tiny (max-displacement s < 0.00020))"#,
    );
    // A polymer bar beside a diamond-stiff one, which carries nearly all
    // of the load: 160 N against the soft bar's 0.016 N and its ridge's.
    // Over the free 15 mm the mean elongation of the soft bar alone
    // is at least 15 * (0.001 + 0.12 / (64 * 4)) / 100 = 0.0002203125 mm.
    let a = first(&fresh(&sc, ["require", &f]).any());
    assert_eq!(a["status"], "unknown", "{a}");
    assert!(has(&a["unavailable"], "pressure load 'ends'"), "{a}");
}

#[test]
fn whole_part_pressure_checks_the_magnitude_between_samples() {
    let sc = Scratch::new("st-pressure-ridge");
    let ridge = sc.model(
        "ridge.ink",
        r#"(import "std")
(def bar (box 16 8 8))
(study pull (shape bar) (elastic (E 10000MPa) (nu 0))
  (fix base (+ x 7))
  (load tip (pressure (* 1MPa (- -0.1 (* 128 (max 0 (- 1 (* 64 (abs (- y 0.03125)))))))))
            (region (face bar (normal 1 0 0))))
  (res 1))
(require little-extension (max-displacement pull < 0.00035))
(require sample (subset (probe 6 0 0) (- (study-field pull displacement) 1)))"#,
    );
    let a = by(&fresh(&sc, ["require", &ridge]).any()["requirements"], "name");
    for row in a.values() {
        assert_eq!(row["status"], "unknown", "{}", pretty(row));
        assert!(has(&row["unavailable"], "pressure"), "{}", pretty(row));
    }
}

#[test]
fn whole_part_pressure_reports_both_sides_of_its_coverage_cutoff() {
    let sc = Scratch::new("st-pressure-fraction");
    // On each affected facet the ridge can add its full amplitude
    // to the sampled background of 0.1 MPa.
    for (amplitude, refused) in [(0.002999, false), (0.003001, true)] {
        let src = format!(
            r#"(import "std")
(def bar (box 16 8 8))
(study s (shape bar) (elastic (E 10000MPa) (nu 0))
 (fix foot (+ x 7))
 (load ridge (pressure (* 1MPa (- -0.1 (* {amplitude} (max 0 (- 1 (* 64 (abs (- y 0.03125))))))))) (region (face bar (normal 1 0 0)))) (res 1))"#
        );
        let f = sc.model("ridge.ink", &src);
        let sol = fresh(&sc, ["solve", &f]).ok();
        let fraction = n(&sol["loads"][0]["pressure_error_fraction"]);
        if refused {
            assert!((fraction - amplitude / 0.1).abs() < 1e-9, "{sol}");
        } else {
            // The narrow triangular ridge contributes 8 * amplitude/64 N.
            assert!(fraction <= 1.01e-8, "{sol}");
            assert!((n(&sol["loads"][0]["total"][0]) - (6.4 + amplitude / 8.0)).abs() < 1e-8, "{sol}");
        }
        assert_eq!(any_has(&sol["warnings"], "pressure load 'ridge'"), refused, "{sol}");
        assert_eq!(any_has(&sol["warnings"], "missed the cells and pressure variation"), refused, "{sol}");
    }
}

#[test]
fn whole_part_a_small_wave_between_samples_cannot_certify_a_close_limit() {
    let sc = Scratch::new("st-pressure-wave");
    let f = sc.model(
        "wave.ink",
        r#"(import "std")
(def bar (move (box 32 8 8) 16 0 0))
(study s (shape bar) (elastic (E 10000MPa) (nu 0))
 (fix base (- x 0.01)) (load end (pressure (* 1MPa (- -1 (* 0.028 (square (sin (* 18.84955592153876 y))))))) (region (face bar (normal 1 0 0)))) (res 1))
(require travel (max-displacement s < 0.003225))"#,
    );
    // The 8 mm face spans 24 full periods of sin^2(6 pi y). Its mean
    // tension is 1.014 MPa. The last 31.99 mm of bar therefore extend
    // at least 0.003243786 mm even before the root displacement.
    let a = first(&fresh(&sc, ["require", &f]).any());
    assert!(matches!(a["status"].as_str(), Some("fails" | "unknown")), "{a}");
    if a["status"] == "unknown" {
        assert!(has(&a["unavailable"], "Refine"), "{a}");
    }
}

#[test]
fn whole_part_a_resolved_pressure_gradient_needs_no_claim_warning() {
    let sc = Scratch::new("st-pressure-gradient");
    // An affine pressure integrates exactly even where its variation
    // across one facet exceeds the coverage cutoff.
    let f = sc.model(
        "gradient.ink",
        r#"(import "std")
(def bar (box 16 8 8))
(study s (shape bar) (elastic (E 10000MPa) (nu 0))
 (fix root (+ x 7)) (load water (pressure (* 1MPa (+ -1 (* y 0.1)))) (region (face bar (normal 1 0 0)))) (res 1))
(require travel (max-displacement s < 0.1))"#,
    );
    let a = first(&fresh(&sc, ["require", &f]).any());
    assert_eq!(a["status"], "holds", "{a}");
    assert!(!any_has(&a["warnings"], "pressure load"), "{a}");
}

#[test]
fn whole_part_physical_supports_are_not_symmetry_cuts() {
    let sc = Scratch::new("st-physical-support");
    // (Its claim a displacement: the peak stress lies at a hold's edge, where
    // it grows as the cells shrink (0.1029, 0.1031, 0.1035 MPa on cells of 4,
    // 2 and 1 mm), and is not settled; the floor is what is tested. A seat
    // held along its normal takes in the top face's edge at the corner,
    // where `top` moves it along z: it holds its own face's normal there.)
    let base = r#"(import "std")
(def c (move (box 4 4 20) 2 2 10))
(study s (shape c) (elastic (E 2000MPa) (nu 0) (yield 200MPa))
  (fix base (- z 0.5)) (fix top (- 19.99 z) (dofs z) (move 0 0 -0.001))
  (fix floor (- y 0.01) (dofs y) (support)) (res 1))
(require stress (max-displacement s < 1))"#;
    for (name, src) in [
        ("roller", base.to_string()),
        ("gravity", base.replace(" (support)", "").replace("(E 2000MPa)", "(E 2000MPa) (density 1000kg/m3)").replace("(res 1)", "(gravity 0 -9.81m/s2 0) (res 1)")),
        ("legacy-floor", base.replace(" (support)", "")),
        ("two-rollers", base.replace("(support))", ") (fix ceiling (- 3.99 y) (dofs y))")),
        ("seat", base.replace("(fix floor (- y 0.01) (dofs y) (support))", "(fix floor (max (- y 0.01) (- z 19.5)) (dofs normal) (support))")),
        ("corner-seat", base.replace("(fix floor (- y 0.01) (dofs y) (support))", "(fix floor (- y 0.01) (dofs normal) (support))")),
    ] {
        let f = sc.model(&format!("{name}.ink"), &src);
        let a = first(&fresh(&sc, ["require", &f]).any());
        assert!(a["status"] == "holds" || name == "seat" && a["status"] == "unknown", "{name}: {}", pretty(&a));
        assert!(!pretty(&a).contains("symmetry cut"));
    }
    // A seat held along its normal that does not say it is a support may be a cut: not settled, saying so.
    let unsaid = first(
        &fresh(&sc, ["require", &sc.model("seat-unsaid.ink", &base.replace("(fix floor (- y 0.01) (dofs y) (support))", "(fix floor (max (- y 0.01) (- z 19.5)) (dofs normal))"))])
            .any(),
    );
    assert!(unsaid["status"] == "unknown" && has(&unsaid["undecided"], "does not say what it is"), "{}", pretty(&unsaid));
}

#[test]
fn whole_part_reflection_checks_loads_and_free_turns() {
    let sc = Scratch::new("st-reflection-loads");
    let octant = r#"(import "std")
(def c (move (box 4 4 4) 2 2 2))
(study s (shape c) (elastic (E 2000MPa) (nu 0) (density 1000kg/m3))
  (mirror mx x 0) (mirror my y 0) (mirror mz z 0)
  (load pull (pressure -1MPa) (region (max (- 3.99 (max x y z)) (- 0.000001 x) (- 0.000001 y) (- 0.000001 z)))) (res 0.5))
(require stable (min-buckling-factor s > 2))"#;
    let model = sc.model("octant.ink", octant);
    let sol = fresh(&sc, ["solve", &model]).ok();
    assert_eq!(sol["results"]["buckling"]["sectors_solved"], 8, "{}", pretty(&sol));
    assert_eq!(first(&fresh(&sc, ["require", &model]).any())["status"], "holds");
    // With zero load the newly freed turns have neither restoring energy
    // nor coupling to any elastic displacement. The full search is valid.
    let unloaded = sc.model("unloaded.ink", &octant.replace("(pressure -1MPa)", "(pressure 0)"));
    assert_eq!(first(&fresh(&sc, ["require", &unloaded]).any())["status"], "holds");
    for (name, case, reason) in [
        ("free support", octant.replace("(res 0.5)", "(fix guide (- x 0.1) (slide 0 1 0)) (res 0.5)"), "free motion on a support that may touch the mirror"),
        ("gravity", octant.replace("(res 0.5)", "(gravity 0 0 -9.81m/s2) (res 0.5)"), "gravity acts across"),
        ("pressure", octant.replace("(max (- 3.99 (max x y z)) (- 0.000001 x) (- 0.000001 y) (- 0.000001 z))", "-1"), "loads a symmetry cut"),
        ("interior", octant.replace("(mirror mx x 0)", "(mirror mx x 2)"), "entire part on one side"),
        ("parallel", octant.replace("(mirror mx x 0)", "(mirror mx x 0) (mirror end x 4)"), "parallel mirrors"),
    ] {
        let f = sc.model(&format!("{name}.ink"), &case);
        let d = fresh(&sc, ["solve", &f]).any();
        assert!(has(&d["error"]["message"], reason), "{name}: {}", pretty(&d));
    }
}

#[test]
fn whole_part_a_pinned_end_keeps_the_rotation_in_its_mirror_plane() {
    let sc = Scratch::new("st-mirror-pin");
    let f = sc.model(
        "pinned.ink",
        r#"(import "std")
(def part (move (box 100 2 4) 50 1 0))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3))
 (fix base (- x 0.01) (pivot (at 0 0 0) (axis 0 1 0)))
 (fix top (- 99.99 x) (slide 1 0 0) (pivot (at 100 0 0) (axis 0 1 0)))
 (mirror cut y 0) (load end (force -500N 0 0) (region (- 99 x))) (res 1))
(require margin (min-buckling-factor s > 3))"#,
    );
    let a = first(&fresh(&sc, ["require", &f]).any());
    // The full 4 x 4 column carries 1000 N. Euler's pinned-end
    // solution is pi^2 E I / L^2, or 4.21 times this load.
    assert_eq!(a["status"], "holds", "{a}");
    assert!((n(&a["value"]) / 4.211 - 1.0).abs() < 0.05, "{a}");
    let sol = fresh(&sc, ["solve", &f]).ok();
    assert_eq!(sol["results"]["buckling"]["sectors_solved"], 2);
}

#[test]
fn whole_part_a_floor_roller_beside_a_hinge_carries_the_pull() {
    let sc = Scratch::new("st-floor-hinge");
    let src = r#"(import "std")
(def part (move (box 20 4 4) 10 2 0))
(study s (shape part) (elastic (E 10000MPa) (nu 0))
 (fix hinge (- x 1) (pivot (at 0 0 0) (axis 0 0 1)) (rigid))
 (fix floor (- y 0.01) (dofs y))
 (load pull (pressure -1MPa) (region (face part (normal 1 0 0)))) (res 1))
(require travel (max-displacement s < 1))"#;
    for model in [src.to_owned(), src.replace("(dofs y)", "(dofs y) (support)")] {
        let f = sc.model("hinge.ink", &model);
        let a = first(&fresh(&sc, ["require", &f]).any());
        assert_eq!(a["status"], "holds", "{a}");
        // One MPa in a 10 GPa bar extends its free 19 mm by 0.0019 mm. (The hinge's block is rigid: held by its
        // mean motion, the block turned about the pin at its corner, 1.5e-4 rad, and the bar read 0.00224.)
        assert!((n(&a["value"]) - 0.0019).abs() < 0.0001, "{a}");
    }
}

#[test]
fn a_load_region_no_sample_lies_in_says_so_not_that_it_holds_no_material() {
    // A 0.1 mm cube of a block loaded, on 2 mm cells set so no sample of
    // theirs falls in it: it was refused as holding no material.
    let sc = Scratch::new("st-unsampled-load");
    let f = sc.model(
        "cube.ink",
        r#"(import "std")
(def part (move (max (move (box 16 16 8) 0 0 4) (- x) (- y)) 0.8284 1.4642 0.4722))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3))
  (fix base (move (- z 0.5) 0.8284 1.4642 0.4722)) (mirror left x 0.8284) (mirror front y 1.4642)
  (load tip (force 0 0 -100N) (region (move (max (- x 0.1) (- y 0.1) (- 7.9 z)) 0.8284 1.4642 0.4722)))
  (res 2))"#,
    );
    let e = sc.ink(["solve", &f]).err();
    let said = e["error"]["message"].as_str().unwrap_or("");
    assert!(said.contains("on none of their samples") && said.contains("Refine (res r)"), "{e}");
}

#[test]
fn a_probe_says_it_is_one_read_and_a_cut_cells_scatters() {
    // A pressurised quarter ring's bore (Lamé 51.05 MPa): its cut cells' own reads scatter 40 to 54 at every
    // cell size, and a probe there printed one with nothing said of how far off it may be.
    let d = ink(["solve", &fixture("stress-reads/vessel.ink"), "--points", "[[34.7,34.7,20],[35.6,35.6,20]]"]);
    let probes = records(&d["probes"]);
    assert!(probes.iter().any(|p| n(&p["cell_material_fraction"]) < 1.0), "{}", pretty(&d["probes"]));
    for p in probes {
        let e = p["estimate"].as_str().unwrap_or("");
        assert!(e.contains("not a bound") && (n(&p["cell_material_fraction"]) >= 1.0 || e.contains("scatters")), "{p}");
    }
}

#[test]
fn a_bar_held_at_its_end_takes_its_hold_where_it_is() {
    // Held at x < 0.1, a 6 mm bar pulled to a strain of 0.01 moves 0.059 at its tip; on 0.1 mm cells it read
    // 0.05801 +- 4.3e-4, and `> 0.0585` failed: its end face's nodes, at x = 6.000000000000001 beside cells
    // 0.99999999999999 full, were read as outside the part, and its peak taken a cell short. It reads 0.059043.
    let sc = Scratch::new("st-hold-reach");
    let f = sc.model(
        "bar.ink",
        r#"(import "std")
(def bar (move (box 6 2 2) 3 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 0.1))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))) (res 0.1))
(require lo (max-displacement s > 0.0585mm))"#,
    );
    let a = first(&sc.ink(["require", &f]).any());
    assert_ne!(a["status"], "fails", "{}", pretty(&a));
}

#[test]
fn a_bar_whose_grids_converge_faster_than_they_can_is_not_given_their_last_difference() {
    // Held at x < 1 and pulled to a strain of 0.01, the bar's tip moves 0.19; on 0.3 mm cells it read
    // 0.18616 +- 0.00046 ("faster than their elements can"), 0.075 mm cells 0.18904, and `> 0.187` failed.
    let sc = Scratch::new("st-faster");
    let f = sc.model(
        "bar.ink",
        r#"(import "std")
(def bar (move (box 20 4 4) 10 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0)) (fix wall (- x 1.0))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))) (res 0.3))
(require reach (max-displacement s > 0.187mm (excluding buckling)))"#,
    );
    let a = first(&sc.ink(["require", &f]).any());
    assert_ne!(a["status"], "fails", "{}", pretty(&a));
}

#[test]
fn a_displacement_its_grids_move_away_from_is_not_settled() {
    // 100 N on a 0.1 mm patch: its displacement read 1.022e-4, 3.262e-4
    // and 7.736e-4 mm on cells of 4, 2 and 1 mm, and 3.556e-3 and 6.336e-3
    // on 0.25 and 0.125. Three times the last difference allowed 2.12e-3,
    // and `< 0.0022` held.
    let sc = Scratch::new("st-patch-load");
    let f = sc.model(
        "patch.ink",
        r#"(import "std")
(def block (move (box 16 16 8) 0 0 4))
(def part (max block (- x) (- y)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3))
  (fix base (- z 0.5)) (mirror left x 0) (mirror front y 0)
  (load tip (force 0 0 -100N) (region (max (- x 0.1) (- y 0.1) (- 7.9 z))))
  (res 1))
(require d (max-displacement s < 0.0022mm))"#,
    );
    let d = sc.ink(["require", &f]).any();
    let a = &d["requirements"][0];
    assert!(a["status"] == "unknown" && has(&a["undecided"], "three grids moving apart") && has(&a["hint"], "settled there if its next change is smaller than"), "{}", pretty(a));
}

#[test]
fn a_displacement_carried_whole_by_a_held_move_is_not_settled_by_the_offset() {
    // The patch above, carried down 1000 mm whole by its base: its deflection
    // still grows as the cells shrink (1000.0003426, 1000.0007679, 1000.0016142
    // on cells of 4, 2 and 1 mm), but against 3.2e-5 of 1000 mm the three
    // grids read alike (1000.0016 +- 0.032), and `< 1000.1` held, false. What
    // the solves resolve is a share of what the deflection varies by.
    let sc = Scratch::new("st-patch-carried");
    let f = sc.model(
        "patch-carried.ink",
        r#"(import "std")
(def block (move (box 16 16 8) 0 0 4))
(def part (max block (- x) (- y)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3))
  (fix base (- z 0.5) (move 0 0 -1000)) (mirror left x 0) (mirror front y 0)
  (load tip (force 0 0 -100N) (region (max (- x 0.1) (- y 0.1) (- 7.9 z))))
  (res 1))
(require d (max-displacement s < 1000.1mm))"#,
    );
    let d = sc.ink(["require", &f]).any();
    let a = &d["requirements"][0];
    assert!(a["status"] == "unknown" && has(&a["undecided"], "three grids moving apart"), "{}", pretty(a));
}

#[test]
fn a_held_temperature_varies_over_its_face() {
    // ANSYS VM108: a cylinder's side held at 80 cos(theta) degF; at x 10 ft on its axis the exact
    // temperature is 40 degF, 4.4444 degC. Its held value read as a constant, the model did not load.
    let f = format!("{ROOT}/benchmarks/standard/ansys-vm108.ink");
    let d = Ink::new(["solve", &f, "--study", "s", "--points", "[[3048,0,0]]"]).run().json();
    let t = d["probes"]["rows"][0][1].as_f64().unwrap_or_else(|| panic!("{}", pretty(&d)));
    assert!((t - 4.444444).abs() < 0.0025, "{t}");
}

#[test]
fn a_temperature_rise_on_a_hot_part_is_not_settled_by_its_kelvins() {
    // The same block held at 1000 degC on its base, 0.01 W into the 0.1 mm
    // patch: its rise grows as the cells shrink (0.041, 0.082 on cells of 2
    // and 1 mm, 0.160 and 0.254 on 0.5 and 0.25), but against 3.2e-5 of 1273
    // K the grids read alike within what the solves resolve (1000.082 +-
    // 0.124), and `< 1000.3degC` held. Held at 0 degC, the same rise moved
    // apart. A temperature carries its 273 K whole, as a held move its offset.
    let sc = Scratch::new("st-patch-hot");
    let f = sc.model(
        "patch-hot.ink",
        r#"(import "std")
(def block (move (box 16 16 8) 0 0 4))
(def part (max block (- x) (- y)))
(study s (shape part) (heat (k 200W/m/K))
  (temperature base (- z 0.5) 1000degC)
  (source tip (power 0.01W) (region (max (- x 0.1) (- y 0.1) (- 7.9 z))))
  (res 1))
(require t (max-temperature s < 1000.3degC))"#,
    );
    let d = sc.ink(["require", &f]).any();
    let a = &d["requirements"][0];
    assert!(a["status"] == "unknown" && has(&a["undecided"], "three grids moving apart"), "{}", pretty(a));
}

// The mesh's rim is chords on the arcs and ends off the cuts: 0.26% short before the face past its
// edges was pressed (face_strips), 0.002% after.
#[test]
fn quarter_tube_receives_its_share_of_axial_pressure() {
    let sc = Scratch::new("st-tube-pressure");
    let f = sc.model(
        "quarter-tube.ink",
        r#"(import "std")
(def tube (subtract (cylinder-z 10 16) (cylinder-z 6 20)))
(def part (max tube (- x) (- y)))
(study s (shape part) (elastic (E 200000MPa) (nu 0))
  (fix root (+ z 7)) (mirror left x 0) (mirror bottom y 0)
  (load end (pressure -1MPa) (region (face part (normal 0 0 1)))) (res 1))
(require stress (max-stress s < 100MPa))"#,
    );
    let sol = fresh(&sc, ["solve", &f]).ok();
    let force = n(&sol["loads"][0]["total"][2]);
    let exact = 16.0 * std::f64::consts::PI;
    assert!((force - exact).abs() < 1e-4 * exact, "{force} N vs {exact} N");
}

#[test]
fn whole_part_steady_claims_say_what_the_solve_dropped() {
    let sc = Scratch::new("st-floating-cube");
    let cubes = sc.model(
        "cubes.ink",
        r#"(import "std")
(def cubes (min (move (box 8 8 8) -8 0 0) (move (box 8 8 8) 8 0 0)))
(study hot (shape cubes) (heat (k 10W/m/K)) (temperature cold x 0)
  (source both (power 10.24W) (region -1)) (res 1))
(require cool (max-temperature hot < 5degC))
(require warm (max-temperature hot > 5degC))
(require sampled (subset (probe 8 0 0) (- (study-field hot temperature) 5)))"#,
    );
    let s = fresh(&sc, ["solve", &cubes]).ok();
    let why = s["warnings"].as_array().unwrap().iter().filter_map(Value::as_str).find(|w| w.contains("was dropped")).unwrap();
    assert!(why.contains("512") && why.contains("not connected to any fix"), "{why}");
    for a in fresh(&sc, ["require", &cubes]).ok()["requirements"].as_array().unwrap() {
        assert_eq!(a["status"], "unknown", "{}", pretty(a));
        assert!(has(&a["unavailable"], why) || has(&a["undecided"], why), "{}", pretty(a));
    }
}

#[test]
fn whole_part_elastic_and_modal_claims_count_an_unheld_piece() {
    let sc = Scratch::new("st-floating-elastic");
    let cubes = sc.model(
        "cubes.ink",
        r#"(import "std")
(def cubes (min (move (box 8 8 8) -8 0 0) (move (box 8 8 8) 8 0 0)))
(study weight (shape cubes) (elastic (E 1000MPa) (nu 0.3) (density 1000kg/m3) (yield 10MPa))
  (fix foot (+ x 11)) (gravity 0 0 -9.81m/s2) (res 1))
(study vibration (shape cubes) (elastic (E 1000MPa) (nu 0.3) (density 1000kg/m3))
  (fix foot (+ x 11)) (modes 1) (res 1))
(require stress (max-stress weight < 100MPa))
(require travel (max-displacement weight < 10))
(require work (compliance weight < 1000N*mm))
(require margin (min-safety-factor weight > 2 (excluding buckling)))
(require buckle (min-buckling-factor weight > 2))
(require frequency (min-frequency vibration > 1Hz))"#,
    );
    let answer = fresh(&sc, ["require", &cubes]).ok();
    for a in answer["requirements"].as_array().unwrap() {
        assert_eq!(a["status"], "unknown", "{}", pretty(a));
        assert!(has(&a["unavailable"], "512.0 mm^3 of material is not connected to any fix and was dropped"), "{}", pretty(a));
    }
}

#[test]
fn a_supported_modal_claim_counts_the_motions_its_supports_leave_free() {
    // A steel bar its supports should hold: its lowest mode is never its
    // free modes' 1285 Hz, read as if it were held. Cut from its clamp by a
    // 2 mm gap it is dropped, not rung free; on one roller where two were
    // meant (a hinge, and slides, at 0 Hz) the claim counts the 0 Hz modes;
    // on two, it decides once it says it leaves its slides out; with no fix
    // at all it is a free part, decided on its elastic modes.
    let sc = Scratch::new("st-modal-free-motions");
    let steel = "(elastic (E 200000MPa) (nu 0.3) (density 7.85g/cm3))";
    let left = "(fix left (+ x 99) (dofs z) (pivot (at -99.5 0 0) (axis 0 1 0)))";
    let right = "(fix right (- 99 x) (dofs z) (pivot (at 99.5 0 0) (axis 0 1 0)))";
    let model = format!(
        r#"(import "std")
(def bar (box 200 10 10))
(def gapped (min bar (move (box 10 10 10) -107 0 0)))
(study cut (shape gapped) {steel} (fix foot (+ x 105)) (modes 2) (res 2))
(study one (shape bar) {steel} {left} (modes 2) (res 2))
(study two (shape bar) {steel} {left} {right} (modes 2) (res 2))
(study free (shape bar) {steel} (modes 2) (res 2))
(require cut-ring (min-frequency cut > 500Hz))
(require one-ring (min-frequency one > 800Hz))
(require two-ring (min-frequency two > 500Hz))
(require two-elastic (min-frequency two > 500Hz (excluding rigid-motions)))
(require free-ring (min-frequency free > 1250Hz))"#
    );
    let answer = fresh(&sc, ["require", &sc.model("bars.ink", &model)]).ok();
    let claim = |name: &str| answer["requirements"].as_array().unwrap().iter().find(|a| a["name"] == name).unwrap().clone();
    let cut = claim("cut-ring");
    assert!(cut["status"] == "unknown" && has(&cut["unavailable"], "not connected to any fix and was dropped"), "{}", pretty(&cut));
    for name in ["one-ring", "two-ring"] {
        let a = claim(name);
        assert_eq!(a["status"], "unknown", "{}", pretty(&a));
        assert!(pretty(&a).contains("rigid-motions") && pretty(&a).contains("excluding rigid-motions"), "{}", pretty(&a));
    }
    for name in ["two-elastic", "free-ring"] {
        let a = claim(name);
        assert_eq!(a["status"], "holds", "{}", pretty(&a));
    }
}

#[test]
fn whole_part_static_claims_count_the_buckle_a_mirror_hides() {
    let sc = Scratch::new("st-mirror-buckle");
    let source = r#"(import "std")
(def bar (move (box 32 16 2) 16 0 1))
(study push (shape bar) (elastic (E 1000MPa) (nu 0) (yield 100MPa))
  (fix foot (- x 1)) (mirror mirror z 0)
  (load end (pressure 5MPa) (region (face bar (normal 1 0 0)))) (res 1))
(require margin (min-safety-factor push > 2))
(require stress (max-stress push < 100MPa))
(require travel (max-displacement push < 10))
(require work (compliance push < 1000N*mm))
(require buckle (min-buckling-factor push > 2))
(require yield-alone (min-safety-factor push > 2 (excluding buckling)))
(require weak (min-safety-factor push > 100))
(require sampled (subset (probe 30 0 1) (- (study-field push displacement) 10)))
(require normalized (subset (probe 30 0 1) (- (study-field push displacement (mode 1)) 2)))"#;
    let column = sc.model("column.ink", source);
    let a = by(&fresh(&sc, ["require", &column]).any()["requirements"], "name");
    for name in ["margin", "stress", "travel", "work"] {
        assert_eq!(a[name]["status"], "unknown", "{name}: {}", pretty(&a[name]));
        assert!(has(&a[name]["undecided"], "linear buckling factor"), "{}", pretty(&a[name]));
    }
    assert_eq!(a["buckle"]["status"], "fails");
    assert_eq!(a["sampled"]["status"], "unknown");
    assert_eq!(a["normalized"]["status"], "holds");
    assert_eq!(a["weak"]["status"], "fails");
    assert_eq!(a["yield-alone"]["status"], "holds");
    assert!(any_has(&a["yield-alone"]["warnings"], "buckl"), "{}", pretty(&a["yield-alone"]));

    // Same axes, stresses and load: the four searches on the quarter must
    // recover the whole column's bending, including the winning shape.
    let mut whole: Option<(f64, f64)> = None;
    let mut common = (0.0_f64, f64::INFINITY);
    for (name, shape, cuts, sectors) in [
        ("whole", "(box 32 16 4) 16 0 0", "", 1),
        ("half", "(box 32 16 2) 16 0 1", "(mirror mirror z 0)", 2),
        ("quarter", "(box 32 8 2) 16 4 1", "(mirror side y 0) (mirror mirror z 0)", 4),
    ] {
        let file = sc.model(&format!("{name}.ink"), &source.replace("(box 32 16 2) 16 0 1", shape).replace("(mirror mirror z 0)", cuts));
        let solve = fresh(&sc, ["solve", &file]).ok();
        let b = &solve["results"]["buckling"];
        let (low, high) = (n(&b["low"]), n(&b["factor"]));
        common = (common.0.max(low), common.1.min(high));
        assert!(common.0 <= common.1, "no common factor after {name}: {b}");
        // (to a rounding: a pressed facet's force is summed over its pieces within the cells)
        let (pressed, meant) = (n(&solve["loads"][0]["total"][0]), -320.0 / f64::from(sectors));
        assert!((pressed / meant - 1.0).abs() < 1e-12, "{name}: {pressed} for {meant}");
        if let Some((whole_low, whole_high)) = whole {
            assert!(low <= whole_high && whole_low <= high, "{name}: {b}, whole {whole:?}");
            assert!((high - whole_high).abs() <= (whole_high - whole_low).min(high - low), "{name}: {b}, whole {whole:?}");
        } else {
            assert!((0.66..0.70).contains(&high));
            whole = Some((low, high));
        }
        if sectors > 1 {
            assert_eq!(b["sectors_solved"], sectors);
            assert_eq!(b["antisymmetric_axes"], json!(["z"]));
            assert!(n(&b["symmetric_factor"]) > 2.0);
        }
    }
    let plane = sc.model("exact-face.ink", &source.replace("(face bar (normal 1 0 0))", "(- 32 x)"));
    let pressure = fresh(&sc, ["solve", &plane]).ok();
    assert!((n(&pressure["loads"][0]["total"][0]) / -160.0 - 1.0).abs() < 1e-12, "{}", pressure["loads"][0]);
    let claims = source.split("(require normalized").next().unwrap();
    for restraint in ["(- z 0.01) (dofs z)", "(- z 0.01) (dofs normal)", "(max (- z 0.01) (- 8 x)) (dofs z)"] {
        let leftover = sc.model("leftover.ink", &claims.replace("(mirror mirror z 0)", &format!("(mirror mirror z 0) (fix old {restraint})")));
        let rows = by(&fresh(&sc, ["require", &leftover]).any()["requirements"], "name");
        assert_eq!(rows["stress"]["status"], "unknown", "{restraint}: {}", rows["stress"]);
        assert!(has(&rows["stress"]["unavailable"], "Remove the leftover fix, or mark it (support)"));
        // In words, each row with no value says why on its own row, though a row above said the same.
        if restraint.starts_with("(- z 0.01) (dofs z)") {
            let said = fresh(&sc, ["require", &leftover, "--plain"]).run().stdout;
            let unavailable: Vec<&str> = ["margin", "stress", "travel", "work"].into_iter().filter(|r| rows[*r]["unavailable"].is_string()).collect();
            assert!(unavailable.len() > 1, "{said}");
            for name in unavailable {
                let row: String = said.lines().skip_while(|l| !l.starts_with(&format!("{name}:"))).skip(1).take_while(|l| l.starts_with(' ')).collect::<Vec<_>>().join("\n");
                assert!(row.contains(&format!("  note: {}", rows[name]["unavailable"].as_str().unwrap())), "{name}:\n{said}");
            }
        }
    }
    let deep = sc.model("deep.ink", &claims.replace("(mirror mirror z 0)", "(fix old (- z 1.75) (dofs z))"));
    let rows = by(&fresh(&sc, ["require", &deep]).any()["requirements"], "name");
    assert_eq!(rows["travel"]["status"], "unknown", "{}", rows["travel"]);
    assert!(has(&rows["travel"]["undecided"], "symmetry cut"));
    // Both the native half and the full column are unstable here. An
    // exclusion cannot make the modeled half's own linear response valid.
    let overload = sc.model(
        "overload.ink",
        &source
            .replace("(pressure 5MPa)", "(pressure 100MPa)")
            .replace("(yield 100MPa)", "(yield 10000MPa)")
            .replace("travel (max-displacement push < 10)", "travel (max-displacement push < 10 (excluding buckling))"),
    );
    let a = by(&fresh(&sc, ["require", &overload]).any()["requirements"], "name");
    for name in ["travel", "yield-alone"] {
        assert_eq!(a[name]["status"], "unknown", "{}", pretty(&a[name]));
        assert!(has(&a[name]["undecided"], "buckling factor"), "{}", pretty(&a[name]));
    }
    let pulled = sc.model("pulled.ink", &source.split("(require normalized").next().unwrap().replace("(pressure 5MPa)", "(pressure -5MPa)"));
    let a = by(&fresh(&sc, ["require", &pulled]).any()["requirements"], "name");
    assert_eq!(a["margin"]["status"], "holds");
    assert_eq!(a["buckle"]["status"], "holds");
}

#[test]
fn a_bar_in_tension_reads_its_stress_beside_a_symmetry_plane() {
    // A 4 x 4 bar under 160 N reads 10 MPa throughout, beside its symmetry
    // plane too: a held layer 0.01 thick, whose free side carries the whole
    // load (its edge held by Nitsche's method). Held at their first
    // samples' centres instead (h/8 in), the cells beside it passed an
    // eighth of the reaction to the next nodes and read 8.7 at every cell
    // size; held over the layer's samples by a penalty, all but t/2h of it.
    let sc = Scratch::new("st-sym-bar");
    let bar = sc.model(
        "bar.ink",
        "(import \"std\")\n(param r 1)\n(study sym (shape (move (box 20 4 4) 10 2 2)) (elastic (E 200000MPa) (nu 0.3))\n  (fix plane (- x 0.01) (dofs x))\n  (fix line (max (- 19.99 x) (- y 0.01)) (dofs y z)) (fix point (max (- 19.99 x) (- z 0.01) (- 3.99 y)) (dofs z))\n  (load end (force 160N 0 0) (region (- 19 x))) (res r))\n",
    );
    for h in [1.0, 0.5] {
        let d = sc.ink(["solve", &bar, "--param", &format!("r={h}"), "--points", "[[0,2,2],[0,0.5,3.5]]"]).ok();
        let expect = 10.0;
        for p in records(&d["probes"]) {
            let s = n(&p["stress_mpa"]["xx"]);
            assert!((s - expect).abs() < 0.005 * expect, "at {h} mm cells, {}: {s} MPa against {expect}", p["point"]);
        }
    }
}

#[test]
fn a_moved_studys_compliance_has_no_value_saying_to_read_the_reaction() {
    // A fix that moves its region leaves compliance no value, saying why.
    let sc = Scratch::new("st-moved");
    let col = sc.model(
        "pressed-column.ink",
        r#"(import "std")
(def c (move (box 3 2 20) 0 0 10))
(study moved (shape c) (elastic (E 2000MPa) (nu 0.3)) (fix base (- z 0.5)) (fix top (- 19.5 z) (dofs z) (move 0 0 -0.001)) (res 0.5))
(study loaded (shape c) (elastic (E 2000MPa) (nu 0.3)) (fix base (- z 0.5)) (load top (force 0 0 -1N) (region (- 19.5 z))) (res 0.5))
"#,
    );
    let m = sc.ink(["solve", &col, "--study", "moved"]).ok()["results"].clone();
    let l = sc.ink(["solve", &col, "--study", "loaded"]).ok()["results"].clone();
    assert!(
        m["compliance_nmm"]["value"].is_null() && has(&m["compliance_nmm"]["unavailable"], "reaction") && l["compliance_nmm"].is_f64(),
        "moved: {}\nloaded: {}",
        pretty(&m["compliance_nmm"]),
        pretty(&l["compliance_nmm"])
    );
}

// ---- queries ----

/// A cube, a ball 30 off, and claims on fields of different slopes.
const Q: &str = r#"(import "std")
(body cube (shape (box 20 20 20)))
(body ball (shape (move (sphere 5) 30 0 0)))
(def a (sphere 10))
(def far (move (sphere 5) 20 0 0))
(def steep (* 7 (sphere 10)))
(require c4 (clearance a far 4))
(require c6 (clearance a far 6))
(require s4 (clearance steep far 4))
"#;

#[test]
fn a_ray_gives_its_first_hit_and_many_rays_answer_in_one_call_or_as_a_table() {
    let sc = Scratch::new("st-ray");
    let q = sc.model("q.ink", Q);
    let d = sc.ink(["ray", &q, "--from", "[-50,0,0]", "--dir", "[1,0,0]"]).ok();
    assert!(d["nearest"]["target"] == "cube" && (n(&d["nearest"]["t"]) - 40.0).abs() < 1e-6, "{}", pretty(&d["nearest"]));
    let rays = "[[[-50,0,0],[1,0,0]], [[30,0,50],[0,0,-2]], [[0,50,40],[0,1,0]]]";
    let d = sc.ink(["ray", &q, "--rays", rays]).ok();
    let r = records(&d["rays"]);
    assert!(
        r.len() == 3
            && (n(&r[0]["nearest"]["t"]) - 40.0).abs() < 1e-6
            && r[1]["nearest"]["target"] == "ball"
            && (n(&r[1]["nearest"]["t"]) - 45.0).abs() < 1e-6
            && r[2]["nearest"]["target"].is_null(),
        "{}",
        pretty(&d)
    );
    let t = sc.ink(["ray", &q, "--rays", rays, "--keys", "nearest.target,nearest.t", "--compact"]).ok()["rays"].clone();
    let firsts: Vec<Value> = t["rows"].as_array().map_or(vec![], |rows| rows.iter().map(|r| r[0].clone()).collect());
    assert!(
        t["columns"] == json!(["nearest.target", "nearest.t"]) && firsts == [json!("cube"), json!("ball"), Value::Null] && t["rows"][2] == json!([null, null]),
        "{}",
        pretty(&t)
    );
}

#[test]
fn mass_gives_mass_centre_and_volume_and_the_inertia_only_when_keys_names_it() {
    // A weight loop reads mass and centre: the tensor (3-4 KB a call for an
    // assembly) is there only when --keys names it.
    let sc = Scratch::new("st-mass-keys");
    let q = sc.model("q.ink", Q);
    let t = sc.ink(["mass", &q, "--target", "cube", "--density", "1g/cm3"]).ok()["targets"][0].clone();
    let c = &t["center_of_mass"];
    assert!(
        (n(&t["mass"]) - 8.0).abs() < 1e-3
            && (0..3).all(|k| n(&c[k]).abs() < 1e-9)
            && t.get("volume").is_some()
            && ["inertia", "principal_moments", "principal_axes"].iter().all(|k| t.get(k).is_none()),
        "{}",
        pretty(&t)
    );
    let d = sc.ink(["mass", &q, "--target", "cube", "--density", "1g/cm3", "--keys", "mass,principal_moments", "--compact"]).ok();
    let t = &d["targets"][0];
    assert!(
        (n(&t["mass"]) - 8.0).abs() < 1e-3
            && (n(&t["principal_moments"][0]) - 8.0 * 800.0 / 12.0).abs() < 0.1
            && t["target"] == "cube"
            && t.get("volume").is_none()
            && d["total"].get("principal_moments").is_some()
            && d["total"].get("inertia").is_none(),
        "{}",
        pretty(&d)
    );
}

/// A stepper: a motor of published mass, and a drive whose stator and steel
/// rotor weigh its published 280 g.
const STEPPER: &str = r#"(import "std")
(body motor (shape (move (box 42 42 40) 0 0 -20)) (attr mass 280g) (frame mount (at)))
(body stator (shape (subtract (move (box 42 42 40) 0 0 -20) (cylinder-z 3 100))) (frame mount (at)))
(body rotor (shape (cylinder-z 2.5 46)) (attr density 7.85g/cm3))
(assembly root (instance motor motor))
(assembly (drive (angle 0)) (attr mass 280g) (instance stator stator) (instance rotor rotor (frame-of stator) (turn z angle)))
"#;

#[test]
fn an_assemblys_mass_over_parts_carrying_theirs_keeps_theirs_and_shares_what_remains() {
    // An assembly given its mass may place parts that carry theirs (a
    // catalog root or drive, an assembly of its own): they keep it, and the
    // rest share what remains (it failed to load: "already weighs as a share
    // of assembly 'c.root''s published mass"). A part's moving pieces weigh
    // its published mass (the drive had no mass: its stator had no density).
    let sc = Scratch::new("st-nest-mass");
    sc.model("stepper.ink", STEPPER);
    sc.model("circlip.ink", "(import \"std\")\n(body ring (shape (max (ring-2d 9.5 12) (slab-z 0 1.2))))\n(assembly root (attr mass 1.3g) (instance ring ring))\n");
    let nest = sc.model(
        "nest-mass.ink",
        r#"(import "std")
(import "stepper.ink" as m)
(import "circlip.ink" as c)
(body plate (shape (box 50 50 5)))
(body cover (shape (move (box 50 50 2) 0 0 10)))
(assembly inner (attr mass 100g) (instance c cover))
(assembly root (attr mass 500g) (instance p plate) (instance i inner) (instance d m.drive (at (move 0 0 50))) (instance r c.root (at (move 80 0 0))))
"#,
    );
    let d = sc.ink(["mass", &nest]).ok();
    let t = masses(&d);
    assert!(
        (n(&d["total"]["mass"]) - 500.0).abs() < 1e-6
            && (t["root/i/c"] - 100.0).abs() < 1e-6
            && (t["root/d/stator"] + t["root/d/rotor"] - 280.0).abs() < 1e-6
            && (t["root/r/ring"] - 1.3).abs() < 1e-6
            && (t["root/p"] - (500.0 - 100.0 - 280.0 - 1.3)).abs() < 1e-6,
        "{}",
        pretty(&d)
    );
    let drive = sc.model("drive.ink", "(import \"stepper.ink\" as m)\n(assembly root (instance a m.drive))\n");
    let d = sc.ink(["mass", &drive]).ok();
    let rotor = n(&by(&d["targets"], "target")["root/a/rotor"]["mass"]);
    assert!((n(&d["total"]["mass"]) - 280.0).abs() < 1e-6 && (rotor / (std::f64::consts::PI * 2.5 * 2.5 * 46.0 * 0.00785) - 1.0).abs() < 0.01, "{}", pretty(&d));
}

#[test]
fn a_published_mass_is_measured_only_where_a_mass_is_read() {
    // Check names the derived densities and gives them no value (it
    // measured them all, 4.4 s of a robot's 1.8 s load); eval measures one
    // only where --keys names it, its default record naming it (it measured
    // every target's and printed the name anyway); mass weighs it
    // (unmeasured, its density is no number: never 0 g).
    let sc = Scratch::new("st-derived");
    let stepper = sc.model("stepper.ink", STEPPER);
    let c = sc.ink(["check", &stepper, "--load"]).ok();
    let b: BTreeMap<String, Value> = by(&c["bodies"], "name").into_iter().map(|(k, v)| (k, v["attrs"]["density"].clone())).collect();
    let p: BTreeMap<String, Value> = by(&c["params"], "name").into_iter().filter(|(_, v)| v.get("derived").is_some()).collect();
    assert!(
        b["motor"] == json!({"derived": ["motor.mass-density"]})
            && b["stator"] == json!({"derived": ["drive.mass-density"]})
            && same(&b["rotor"], &json!(0.00785))
            && b.len() == 3
            && p.keys().eq(["drive.mass-density", "motor.mass-density"])
            && p.values().all(|x| x["value"].is_null()),
        "{}\n{}",
        pretty(&json!(b)),
        pretty(&json!(p))
    );
    let pm = sc.model("eval-mass.ink", "(import \"std\")\n(body b (shape (box 10 10 20)) (attr mass 10g))\n");
    let e0 = sc.ink(["eval", &pm, "--points", "[[0,0,0]]"]).ok();
    let e1 = sc.ink(["eval", &pm, "--points", "[[0,0,0]]", "--keys", "attrs.density"]).ok();
    assert!(
        records(&e0["targets"][0]["results"])[0]["attrs"]["density"] == json!({"derived": ["b.mass-density"]})
            && e0.get("solved_inputs").is_none()
            && (n(&records(&e1["targets"][0]["results"])[0]["attrs"]["density"]) - 0.005).abs() < 0.0001
            && e1["solved_inputs"] == json!(["b.mass-density"]),
        "{}\n{}",
        pretty(&e0),
        pretty(&e1)
    );
    let v = n(&sc.ink(["mass", &pm]).ok()["targets"][0]["mass"]);
    assert!((v - 10.0).abs() < 0.05, "{v}");
}

#[test]
fn a_published_mass_follows_its_param_and_a_part_sized_by_its_roots_arguments_weighs_its_files_mass_over_its_own_volume() {
    // A mass the published one reads follows its value (162 g at --param
    // M=172). A part sized by its root's arguments weighs its file's mass
    // at its defaults over its own volume: a 3000 mm guide (2640 g) at
    // 1000 mm weighs 880 g, and one following this file's travel, 3000 mm
    // at travel 2940, 2640 g.
    let sc = Scratch::new("st-param-mass");
    let pm = sc.model("param-mass.ink", "(import \"std\")\n(param M 162)\n(body bar (shape (box 20 300 10)) (attr mass (* 1g M)))\n");
    let m = n(&sc.ink(["mass", &pm, "--param", "M=172"]).ok()["total"]["mass"]);
    assert!((m - 172.0).abs() < 1e-6, "{m}");
    sc.model(
        "guide.ink",
        "(import \"std\")\n(param length 3000)\n(body (part (length length)) (shape (box length 32 30)))\n(assembly (root (length length)) (attr mass 2640g) (instance part (part length)))\n",
    );
    let ga = sc.model(
        "guide-args.ink",
        "(param travel 400)\n(import \"guide.ink\" as g)\n(assembly root (instance a (g.root 1000)) (instance b (g.root (+ travel 60)) (at (move 0 50 0))))\n",
    );
    let t = masses(&sc.ink(["mass", &ga, "--param", "travel=2940"]).ok());
    assert!((t["root/a/part"] / 880.0 - 1.0).abs() < 0.02 && (t["root/b/part"] / 2640.0 - 1.0).abs() < 0.02, "{t:?}");
}

#[test]
fn a_body_built_with_arguments_weighs_its_definitions_density() {
    // Measured at its defaults: a length passed as an argument scales the
    // mass as a param does ((rod 300) weighed 100 g, as at its default
    // 100). A mass that reads the argument is its own.
    let sc = Scratch::new("st-rod-mass");
    let rods = sc.model(
        "rod-mass.ink",
        r#"(import "std")
(body (rod (length 100)) (shape (box length 10 10)) (attr mass 100g))
(body (bar (length 100)) (shape (box length 10 10)) (attr mass (* length 1g)))
(assembly root (instance a rod) (instance b (rod 300) (at (move 0 50 0))) (instance c (bar 300) (at (move 0 100 0))))
"#,
    );
    let t = masses(&sc.ink(["mass", &rods]).ok());
    assert!((t["root/a"] - 100.0).abs() < 1e-6 && (t["root/b"] / 300.0 - 1.0).abs() < 0.02 && (t["root/c"] - 300.0).abs() < 1e-6, "{t:?}");
}

/// A tube 100 long with 0.1 walls, and a plate 0.25 thick.
const THIN: &str = "(import \"std\")\n(body tube (shape (subtract (cylinder-z 2 100) (cylinder-z 1.9 101))) (attr density 2.7g/cm3))\n(body plate (shape (box 200 200 0.25)) (attr density 2.7g/cm3))\n";

#[test]
fn mass_and_volume_are_measured_to_two_percent_by_halving_cells() {
    // Walls thinner than any cell are proven: the tube (0.1 walls, 100
    // long) and the plate (0.25 thick, 200 across) each settle within 2%,
    // their proven bounds round the exact mass (the plate was unsettled:
    // no cell fit inside it). Volume is mass's, with its bounds.
    let sc = Scratch::new("st-thin");
    let mat = sc.model("thin.ink", THIN);
    let exact = [("tube", std::f64::consts::PI * (4.0 - 1.9 * 1.9) * 100.0), ("plate", 200.0 * 200.0 * 0.25)];
    let m = by(&sc.ink(["mass", &mat]).ok()["targets"], "target");
    assert!(
        exact.iter().all(|(k, e)| {
            let (t, mass) = (&m[*k], 0.0027 * e);
            (n(&t["mass"]) / mass - 1.0).abs() < 0.02
                && (n(&t["volume"]["estimate"]) / e - 1.0).abs() < 0.02
                && n(&t["mass_bounds"][0]) <= mass
                && mass <= n(&t["mass_bounds"][1])
                && t.get("spread").is_none()
                && t.get("unsettled").is_none()
                && t.get("sampled").is_none()
        }),
        "{}",
        pretty(&json!(m))
    );
}

#[test]
fn a_distance_is_bracketed_and_clearance_proofs_are_sound_for_scaled_fields() {
    // A field 7 times steeper than a distance proves no more than one: c6
    // fails as for the plain sphere.
    let sc = Scratch::new("st-distance");
    let q = sc.model("q.ink", Q);
    let d = sc.ink(["distance", &q, "--a", "cube", "--b", "ball"]).ok()["distance"].clone();
    let (lo, hi) = (n(&d["lower_proven"]), n(&d["upper_witness"]));
    assert!(lo <= 15.0 && 15.0 <= hi + 1e-6 && hi - lo < 0.3, "{}", pretty(&d));
    let st = statuses(&sc.ink(["require", &q]).ok());
    assert!(st == BTreeMap::from([("c4".into(), "holds".into()), ("c6".into(), "fails".into()), ("s4".into(), "holds".into())]), "{st:?}");
}

#[test]
fn a_stress_peak_at_the_edge_of_a_support_margin_says_whose_and_how_far() {
    // Guide premise 2: max_von_mises_mpa leaves out the cells within the
    // support margin, and its peak sat in the next cell, the clamp's
    // singular stress past the margin, with nothing saying so. A
    // cantilever's peak is at its root, beside the clamp's margin: the
    // answer names the fix and the distance to its nearest left-out cell,
    // and the stress claims carry the warning.
    let sc = Scratch::new("st-beside-fix");
    let m = sc.model(
        "beam.ink",
        "(import \"std\")\n(body beam (shape (move (box 42 10 4) 19 0 0)))\n\
        (study bend (shape (shape-of beam)) (elastic (material aluminum-6061)) (fix clamp x) (load tip (force 0 0 -50N) (region (- 39 x))) (res 1))\n\
        (require strong (max-stress bend < 1000MPa))\n(require safe (min-safety-factor bend > 1))\n(require stiff (max-displacement bend < 5))\n",
    );
    let d = sc.ink(["solve", &m, "--study", "bend"]).ok();
    let peak = &d["results"]["max_von_mises_mpa"];
    let a = sc.ink(["require", &m]).ok();
    let warned = |name: &str| {
        a["requirements"]
            .as_array()
            .unwrap()
            .iter()
            .find(|x| x["name"] == name)
            .is_some_and(|x| x["warnings"].as_array().is_some_and(|w| w.iter().any(|w| has(w, "at the edge of the excluded zone around fix 'clamp'"))))
    };
    // What the stress claim leaves out says its number for what it is:
    // the claim's value were the margin counted, never a bare `value`
    // that reads as a second reading of the claim's.
    let strong = a["requirements"].as_array().unwrap().iter().find(|x| x["name"] == "strong").unwrap();
    let margin = strong["excludes"].as_array().and_then(|e| e.iter().find(|e| e["key"] == "support-margin")).unwrap_or_else(|| panic!("{}", pretty(strong)));
    assert!(margin.get("value").is_none() && n(&margin["value_if_counted"]) >= n(&strong["value"]), "{}", pretty(strong));
    assert!(
        peak["nearest_excluded"]["fix"] == "clamp"
            && n(&peak["nearest_excluded"]["distance_mm"]) <= 3f64.sqrt()
            && has(&peak["note"], "at the edge of the excluded zone around fix 'clamp'")
            && warned("strong")
            && warned("safe")
            && !warned("stiff"),
        "{}\n{}",
        pretty(peak),
        pretty(&a)
    );
}

// ---- bodies outside the window in their own frame; claims on numbers ----

#[test]
fn a_body_modelled_outside_the_window_is_analysed_where_it_is_placed() {
    let sc = Scratch::new("st-own-frame");
    let o = sc.model("own-frame.ink", "(import \"std\")\n(body low (shape (move (box 10 10 10) 0 0 -600)))\n(assembly root (instance low low (at (move 0 0 660))))\n");
    let w = sc.ink(["check", &o]).ok()["warnings"].clone();
    let d = sc.ink(["bounds", &o]).ok()["targets"][0].clone();
    assert!(w == json!([]) && (n(&d["min"][2]) - 55.0).abs() < 0.5 && (n(&d["max"][2]) - 65.0).abs() < 0.5, "{}\n{}", pretty(&w), pretty(&d));
}

#[test]
fn overlaps_and_value_claims_hold_and_fail_as_their_numbers_say() {
    let sc = Scratch::new("st-claims");
    let o = sc.model(
        "claims.ink",
        r#"(import "std")
(param hole 5) (param pin 5.05)
(def plate (subtract (box 30 30 6) (cylinder-z hole 20)))
(require press-fit (overlaps plate (cylinder-z pin 30)) (accuracy 0.005))
(require no-fit (overlaps plate (cylinder-z 4.5 30)))
(require ratio (value (/ pin hole) < 1.02))
(require ratio-bad (value (/ pin hole) < 1.005))
"#,
    );
    let st = statuses(&sc.ink(["require", &o]).ok());
    let want: BTreeMap<String, String> = [("press-fit", "holds"), ("no-fit", "fails"), ("ratio", "holds"), ("ratio-bad", "fails")].map(|(a, b)| (a.into(), b.into())).into();
    assert_eq!(st, want);
}

// ---- study claims: param bounds, linear kinematics ----

#[test]
fn a_linear_study_warns_past_5_degrees_of_rotation_and_its_claims_bounds_are_param_expressions() {
    // A PETG leaf (x -2..20, clamped at x < 0) bent by an end moment to
    // `angle` degrees: M = E I theta / a, a = 19.75; tip travel ~ a theta /
    // 2. Its claims carry the warning; their bounds read params. Bent
    // through its thickness, about its weak axis, a leaf has no lateral
    // buckle, so the claim leaves it out.
    let sc = Scratch::new("st-leaf");
    let o = sc.model(
        "leaf.ink",
        r#"(import "std")
(param angle 30) (param limit 3)
(body leaf (shape (move (box 22 1 4) 9 0 0)))
(study turn (shape (shape-of leaf)) (elastic (E 2100MPa) (nu 0.38))
  (fix clamp x) (load twist (torque 0 0 (* 1N*mm (* (/ (* 2100 (/ 1 3)) 19.75) (/ (* pi angle) 180)))) (region (- 19.5 x))) (res 0.25))
(require travel (max-displacement turn < (* 2 limit) (excluding buckling)))
(require small (value angle < (* 15 limit)))
"#,
    );
    let d = sc.ink(["solve", &o]).ok();
    let w: Vec<&Value> = d["warnings"].as_array().map_or(vec![], |w| w.iter().filter(|x| has(x, "large rotation")).collect());
    assert!(w.len() == 1 && (n(&d["results"]["max_rotation_deg"]["value"]) - 30.0).abs() < 1.5, "{}\n{}", pretty(&d["results"]["max_rotation_deg"]), pretty(&d["warnings"]));
    // At 2 degrees, no warning.
    let d = sc.ink(["solve", &o, "--param", "angle=2"]).ok();
    assert!(
        d["warnings"].as_array().is_some_and(|w| w.is_empty()) && (n(&d["results"]["max_rotation_deg"]["value"]) - 2.0).abs() < 0.1,
        "{}\n{}",
        pretty(&d["results"]["max_rotation_deg"]),
        pretty(&d["warnings"])
    );
    let a = sc.ink(["require", &o]).ok();
    let travel = by(&a["requirements"], "name")["travel"].clone();
    assert!(any_has(&travel["warnings"], "large rotation"), "{}", pretty(&travel));
    let want = |t: &str, s: &str| BTreeMap::from([("travel".to_string(), t.to_string()), ("small".to_string(), s.to_string())]);
    assert_eq!(statuses(&a), want("holds estimate", "holds"));
    assert_eq!(statuses(&sc.ink(["require", &o, "--param", "limit=1.9"]).ok()), want("fails estimate", "fails"));
}

// ---- for: generated defs, instances and loads; check's listing ----

/// Defs, instances generated by `for`.
const FOR: &str = r#"(import "std")
(for r (vec 1 2 3) (def ball-$r (sphere r)))
(body flange (shape (subtract (cylinder-z 30 6) (cylinder-z 10 10))))
(body bolt (shape (cylinder-z 2 10)))
(assembly root (instance base flange) (for i (range 6) (instance bolt-$i bolt (at (move 20 0 0) (rot-z (* 60 i))))))
"#;

#[test]
fn for_generates_defs() {
    let sc = Scratch::new("st-for");
    let o = sc.model("for.ink", FOR);
    let d = sc.ink(["check", &o]).ok();
    let names: Vec<&str> = d["defs"].as_array().map_or(vec![], |a| a.iter().filter_map(|x| x["name"].as_str()).collect());
    assert_eq!(names, ["ball-1", "ball-2", "ball-3"], "{}", pretty(&d["defs"]));
}

#[test]
fn for_generates_study_loads() {
    let sc = Scratch::new("st-for-study");
    sc.model("for.ink", FOR);
    let studied = sc.model(
        "for-study.ink",
        r#"(import "for.ink")
(study twist (shape (shape-of flange)) (elastic (E 2000MPa) (nu 0.3)) (fix bore (- (len2 x y) 11))
  (for i (range 4) (load tab-$i (force (* -10N (sin (* i (/ pi 2)))) (* 10N (cos (* i (/ pi 2)))) 0)
    (region (- (len2 (- x (* 28 (cos (* i (/ pi 2))))) (- y (* 28 (sin (* i (/ pi 2)))))) 2)))))
"#,
    );
    let d = sc.ink(["solve", &studied, "--study", "twist", "--res", "2"]).ok();
    let names: Vec<&str> = d["loads"].as_array().map_or(vec![], |a| a.iter().filter_map(|x| x["name"].as_str()).collect());
    assert!(names == ["tab-0", "tab-1", "tab-2", "tab-3"] && n(&d["equilibrium_imbalance"]) < 1e-6, "{}", pretty(&d["loads"]));
}

/// A part importing big JSON data (the kernel's tests/fixtures/import-summary.ink beside
/// its tables), and a design importing it; the design's path.
fn import_summary(sc: &Scratch) -> String {
    let part = std::fs::read_to_string(kernel_fixture("import-summary.ink")).unwrap();
    sc.model("summary/part.ink", &part);
    sc.model("summary/standards.json", &json!({"rows": (0..20000).collect::<Vec<_>>()}).to_string());
    sc.model("summary/record.json", &json!({"mechanical": {"density-g-mm3": 0.00785}}).to_string());
    sc.path("summary")
}

#[test]
fn check_summarizes_a_big_imported_value_by_its_keys_and_size() {
    // Imported data tables are summarized by keys and size, not expanded.
    let sc = Scratch::new("st-summary");
    let summary = import_summary(&sc);
    let imp = sc.model("imports-screw.ink", &format!("(import \"std\")\n(import \"{summary}/part.ink\" as screw)\n"));
    let r = sc.ink(["check", &imp, "--compact", "--imports"]).text();
    let d: Value = serde_json::from_str(&r).unwrap();
    let def = by(&d["defs"], "name")["screw.standards"].clone();
    assert!(def["value"]["data"].get("keys").is_some() && n(&def["value"]["data"]["bytes"]) > 100000.0 && r.len() < 50000, "{} bytes: {}", r.len(), pretty(&def));
}

#[test]
fn check_lists_the_files_own_names_and_an_import_in_brief_and_eval_leaves_out_a_record_unless_asked() {
    // A design importing a part reads its own names; the part is a line of
    // imports: its bodies with their frames, and counts. Eval's attrs leave
    // out a record (spec) unless --keys names it.
    let sc = Scratch::new("st-brief");
    let summary = import_summary(&sc);
    let imp = sc.model("imports-screw.ink", &format!("(import \"std\")\n(import \"{summary}/part.ink\" as screw)\n(def gap 2)\n(assembly root (instance s screw.part))\n"));
    let r = sc.ink(["check", &imp, "--compact"]).text();
    let d: Value = serde_json::from_str(&r).unwrap();
    let names = |k: &str| -> Vec<String> { d[k].as_array().map_or(vec![], |a| a.iter().filter_map(|x| x["name"].as_str().map(String::from)).collect()) };
    let part = d["imports"].as_array().and_then(|i| i.iter().find(|i| i["file"].as_str().is_some_and(|f| f.ends_with("summary/part.ink")))).cloned().unwrap_or_default();
    assert!(
        names("defs") == ["gap"]
            && names("assemblies") == ["root"]
            && d["fns"].as_array().is_some_and(|a| a.is_empty())
            && d["bodies"].as_array().is_some_and(|a| a.is_empty())
            && part["bodies"] == json!({"screw.part": ["seat"]})
            && n(&part["defs"]) > 3.0
            && !r.contains("screw.standards")
            && r.len() < 3000,
        "{} bytes: {}",
        r.len(),
        pretty(&d["imports"])
    );
    let at = |keys: &[&str]| records(&sc.ink([&["eval", imp.as_str(), "--points", "[0,0,-1]"][..], keys].concat()).ok()["targets"][0]["results"])[0]["attrs"].clone();
    let plain = at(&[]);
    let asked = at(&["--keys", "attrs.spec.mechanical"]);
    assert!(plain.get("color").is_some() && plain.get("spec").is_none() && asked["spec"]["mechanical"].get("density-g-mm3").is_some(), "{}\n{}", pretty(&plain), pretty(&asked));
}

#[test]
fn check_lists_a_def_holding_a_function_with_the_fns_by_its_params() {
    let sc = Scratch::new("st-callable");
    let f = sc.model("callable-defs.ink", "(import \"std\")\n(fn cube (a) (box a a a))\n(def sq (lambda (a) (* a a)))\n(def top max)\n(def n 2)\n");
    let d = sc.ink(["check", &f, "--load"]).ok();
    let defs: Vec<&str> = d["defs"].as_array().map_or(vec![], |a| a.iter().filter_map(|x| x["name"].as_str()).collect());
    assert!(d["fns"] == json!(["(cube a)", "(sq a)", "(top ...)"]) && defs == ["n"], "{}\n{}", pretty(&d["fns"]), pretty(&d["defs"]));
}

#[test]
fn a_claim_is_cached_only_where_its_answer_came_from_the_cache() {
    // `ink check` first cached the kernel answers a study claim reads, not
    // the study's solve, which was not counted as recomputed: the assert
    // after it solved for seconds and said cached. A claim computed with no
    // cache at all (a param's value) said cached too.
    let s = Scratch::new("claim-cached");
    let m = s.model(
        "beam.ink",
        "(import \"std\")\n(param h 4)\n(body beam (shape (move (box 42 10 h) 19 0 0)))\n\
        (study bend (shape (shape-of beam)) (elastic (material aluminum-6061)) (fix clamp x) (load tip (force 0 0 -50N) (region (- 39 x))) (res 2))\n\
        (require stiff (max-displacement bend < 5))\n(require thick (value h > 3))\n",
    );
    let cache = s.path("cache");
    let run = |cmd: &str| s.ink([cmd, &m]).env("XDG_CACHE_HOME", &cache).env_remove("INK_NO_CACHE").ok();
    run("check");
    let cached =
        |d: &Value| -> BTreeMap<String, Value> { d["requirements"].as_array().unwrap().iter().map(|a| (a["name"].as_str().unwrap().to_string(), a["cached"].clone())).collect() };
    let (first, again) = (run("require"), run("require"));
    assert_eq!(cached(&first), [("stiff".into(), json!(false)), ("thick".into(), json!(false))].into(), "{}", pretty(&first));
    assert_eq!(cached(&again), [("stiff".into(), json!(true)), ("thick".into(), json!(false))].into(), "{}", pretty(&again));
}

#[test]
fn a_study_claim_on_a_toleranced_part_says_which_band_and_the_size_to_check_at() {
    // A study solves one design: a tolerance on its shape leaves its claims
    // unknown on a plain assert, the answer naming the param, its band, and
    // the flags that check it at one end of it, which then decide it.
    let s = Scratch::new("st-toleranced");
    let m = s.model(
        "column.ink",
        r#"(import "std")
(param w 3 (tol 0 0.1))
(def c (move (box w 2 20) 0 0 10))
(study loaded (shape c) (elastic (E 2000MPa) (nu 0.3)) (fix base (- z 0.5)) (load top (force 0 0 -1N) (region (- 19.5 z))) (res 0.5))
(require stiff (max-displacement loaded < 1))
"#,
    );
    let a = s.ink(["require", &m]).ok()["requirements"][0].clone();
    let at_end = s.ink(["require", &m, "--param", "w=3.1:3.1"]).ok()["requirements"][0].clone();
    assert!(
        a["status"] == "unknown" && has(&a["unavailable"], "w (3..3.1)") && has(&a["hint"], "--param w=3:3 or --param w=3.1:3.1") && at_end["status"] == "holds",
        "{}\n{}",
        pretty(&a),
        pretty(&at_end)
    );
}

/// A sheet 0.189 mm thick lying on z = 0 off a block held at x < 1, at
/// res 0.5: its cells split to 0.0625, and a mesh of the surface at half
/// the grid's cells has its planes at z 0 and 0.25, none on the sheet.
const SHEET: &str = "(import \"std\")\n(def part (min (move (box 4 6 6) 2 0 0) (move (box 20 4 0.189) 10 0 0.0945)))\n";

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_pressure_on_a_resolved_sheet_presses_all_of_it() {
    // 0.001 MPa on every face above z 0.0945: the sheet's top (64 mm2)
    // and the block's (24), 0.088 N down. On the grid's mesh the sheet had
    // no facet: 0.024 N was applied, and the sheet's tip moved 1.3e-7 mm
    // against q L^4 / 8 E I = 0.0728, so a claim under half that held.
    let sc = Scratch::new("st-press-sheet");
    let m = sc.model(
        "sheet.ink",
        &format!(
            "{SHEET}(study s (shape part) (elastic (E 200000MPa) (nu 0)) (fix root (- x 1))\n  (load top (pressure 0.001MPa) (region (- 0.0945 z))) (res 0.5))\n(require stiff (max-displacement s < 0.0364))\n"
        ),
    );
    let d = sc.ink(["solve", &m]).ok();
    let total = &d["loads"][0]["total"];
    assert!((n(&total[2]) + 0.088).abs() < 1e-4 && n(&total[0]).abs() < 1e-3, "{}", pretty(&d["loads"]));
    // The faces in the region, 149.1 mm2, as its facets covered them and
    // as lines through the cells count them.
    let sm = &d["surface_measures"][0];
    assert!(n(&sm["shortfall_mm2"]) == 0.0 && (n(&sm["covered_mm2"]) - 149.1).abs() < 1.5 && (n(&sm["lines_mm2"]) - 149.1).abs() < 1.5, "{}", pretty(sm));
    let a = &sc.ink(["require", &m]).any()["requirements"][0];
    // It fails, but each grid splits the sheet to the same cells: no three
    // grids check it, so it is unknown.
    assert!(a["status"] == "unknown" && a["grid_error"].is_null() && (n(&a["value"]) - 0.0728).abs() < 0.1 * 0.0728, "{}", pretty(a));
}

#[test]
fn a_film_on_a_resolved_sheet_cools_all_of_it() {
    // The sheet as a heat fin: 0.01 W in at the block's back, 25 W/(m2 K)
    // to 20 C over the 149.1 mm2 above z 0.0945. The surface's mean is 20
    // + Q / (h A) = 22.68; the sheet's fin efficiency (0.82) puts the
    // block at about 22.9. On the grid's mesh the film missed the sheet's
    // top (85 mm2) and the block read 24.9. The claim's grid of 2 mm
    // leaves the sheet thin and its film misses it (83 of 149 mm2, as the
    // lines say), so that grid cannot check the claim: it holds, or is
    // unknown, never fails.
    let sc = Scratch::new("st-film-sheet");
    let m = sc.model(
        "sheet.ink",
        &format!(
            "{SHEET}(study s (shape part) (heat (k 50W/m/K))\n  (convection air (- 0.0945 z) (coefficient 25W/m2/K) (ambient 20degC))\n  (source back (power 0.01W) (region (- x 1))) (res 0.5))\n(require cool (max-temperature s < 23.5degC))\n"
        ),
    );
    let d = sc.ink(["solve", &m]).ok();
    let t = n(&d["results"]["max_temperature"]["value"]);
    assert!(t > 22.68 && t < 23.05 && n(&d["surface_measures"][0]["shortfall_mm2"]) == 0.0, "{t}: {}", pretty(&d["surface_measures"]));
    let a = &sc.ink(["require", &m]).any()["requirements"][0];
    assert!((a["status"] == "holds" || a["status"] == "unknown") && n(&a["value"]) > 22.68 && n(&a["value"]) < 23.05, "{}", pretty(a));
}

#[test]
fn excess_pressure_area_cannot_certify_a_close_displacement_limit() {
    let sc = Scratch::new("st-extra-strip");
    let f = sc.model(
        "strip-close.ink",
        r#"(import "std")
(def part (move (box 40 4 2) 20 0 0))
(param pressure-value 0.001)
(study s (shape part) (elastic (E 200000MPa) (nu 0)) (fix root (- x 1))
 (load strip (pressure (* 1MPa pressure-value)) (region (max (- 0.9 z) (- 0.05 y) (- y 0.25) (- 1 x)))) (res 0.5))
(require travel (max-displacement s > 0.000125))"#,
    );
    // The top receives 0.001 * 39 * 0.2 = 0.0078 N. Sampled at four points a facet it took 0.00975 N
    // and a spurious 0.000136 mm (left unknown by the excess area its rim check found); by the area each
    // facet holds it takes the 0.0078 N, the displacement 0.000136 x 0.0078 / 0.00975 = 0.000109 mm,
    // under the limit: no.
    let a = first(&fresh(&sc, ["require", &f]).any());
    assert!(a["status"] == "fails" && (n(&a["value"]) / 0.0001088 - 1.0).abs() < 0.01, "{a}");
    let unloaded = first(&fresh(&sc, ["require", &f, "--param", "pressure-value=0"]).any());
    assert_eq!(unloaded["status"], "fails", "{unloaded}");
}

#[test]
fn a_load_on_a_strip_narrower_than_the_lines_is_measured_and_unknown_only_where_a_grid_misses_it() {
    // A 40 x 4 x 2 bar held at x < 1, 0.001 MPa on a top strip y in
    // [b, b + w] along it, at res 0.5 (the lines a quarter of a mm apart).
    // Counted at points, the lines read such a strip as no line (0) or a
    // whole one (9.75 mm2), whatever its area, and a true claim on it went
    // unknown; counted across the strip's rims on sixteen sub-lines, they
    // measure its area, 39 w, to half a sub-line at each rim. The claim
    // (true) holds, or is unknown where a grid's facets do take less than
    // the strip (their four points a facet: 2.45 of 3.9 mm2 at b 0.1 w 0.1
    // on cells of 0.5; 1.22 of 1.95 at w 0.05 on cells of 0.25), saying so.
    let sc = Scratch::new("st-strips");
    for (b, w) in [(0.1, 0.1), (0.1, 0.05), (0.1, 0.15), (0.05, 0.2), (0.0, 0.1)] {
        let m = sc.model(
            &format!("strip-{b}-{w}.ink"),
            &format!(
                "(import \"std\")\n(def part (move (box 40 4 2) 20 0 0))\n(study s (shape part) (elastic (E 200000MPa) (nu 0)) (fix root (- x 1))\n  (load strip (pressure 0.001MPa) (region (max (- 0.9 z) (- {b} y) (- y {}) (- 1 x)))) (res 0.5))\n(require stiff (max-displacement s < 1))\n",
                b + w
            ),
        );
        let sm = &sc.ink(["solve", &m]).ok()["surface_measures"][0];
        assert!((n(&sm["lines_mm2"]) - 39.0 * w).abs() <= n(&sm["uncertainty_mm2"]) && n(&sm["uncertainty_mm2"]) < 0.25 * 39.0 * w.max(0.1), "y {b}, w {w}: {}", pretty(sm));
        let a = &sc.ink(["require", &m]).any()["requirements"][0];
        assert!(
            a["status"] == "holds"
                || a["status"] == "unknown" && (has(&a["undecided"], "its integral over the surface missed") || has(&a["undecided"], "beyond the pressure region")),
            "y {b}, w {w}: {}\n{}",
            pretty(sm),
            pretty(a)
        );
    }
}

#[test]
fn a_line_down_a_strip_s_middle_counts_the_strip_not_the_line() {
    // A strip 0.05 wide across the bar's top, x 38.1 to 38.15, a line of
    // the lines (a quarter of a mm apart) down its middle, with the root's
    // top x 1 to 1.25 (1 mm2): 1.2 mm2. Its rims tie on that line, the
    // region's gradient along the surface is none there, and the strip was
    // counted as the line's whole cross-section, 1 mm2 (5 times over).
    let sc = Scratch::new("st-strip-middle");
    let m = sc.model(
        "middle.ink",
        "(import \"std\")\n(def part (move (box 40 4 2) 20 0 0))\n(study s (shape part) (elastic (E 200000MPa) (nu 0)) (fix root (- x 1))\n  (load p (pressure 0.001MPa) (region (min (max (- 0.9 z) (- 38.1 x) (- x 38.15)) (max (- 0.9 z) (- 1 x) (- x 1.25))))) (res 0.5))\n",
    );
    let sm = &sc.ink(["solve", &m]).ok()["surface_measures"][0];
    assert!((n(&sm["lines_mm2"]) - 1.2).abs() < 0.3, "{}", pretty(sm));
}

#[test]
fn a_load_or_film_its_facets_miss_on_contact_strips_never_holds_a_claim() {
    // A 40 x 4 x 2 bar held at x < 1, loaded (0.001 MPa) or heated (a film
    // to 100 C, h 5000) on nine strips 0.05 wide across its top, 4 apart,
    // and the root's top x 1 to 1.25: at these phases the facets (four
    // points each) miss every strip, and only the root's patch is loaded
    // or heated. The tip moves 2.75e-5 mm (beam theory over the nine
    // forces), the root's stress 0.0128 MPa, the far end about 47 C. Counted
    // in proportion to the surface missed, the claims held on the root's
    // patch alone (the missed loads' lever 1,500 times the count); a short
    // integral leaves them unknown, or a finer grid fails them.
    let sc = Scratch::new("st-rollers");
    let strips = |phase: f64| {
        let one = |x0: f64| format!("(max (- 0.9 z) (- {x0} x) (- x {}))", x0 + 0.05);
        let all: Vec<String> = (0..9).map(|k| one(6.0 + phase + 4.0 * k as f64)).chain(["(max (- 0.9 z) (- 1 x) (- x 1.25))".to_string()]).collect();
        format!("(min {})", all.join(" "))
    };
    for phase in [0.1, 0.225] {
        let roller = sc.model(
            &format!("rollers-{phase}.ink"),
            &format!(
                "(import \"std\")\n(def part (move (box 40 4 2) 20 0 0))\n(study s (shape part) (elastic (E 200000MPa) (nu 0)) (fix root (- x 1))\n  (load p (pressure 0.001MPa) (region {})) (res 0.5))\n(require stiff (max-displacement s < 1e-5))\n(require strong (max-stress s < 0.01MPa))\n",
                strips(phase)
            ),
        );
        let hot = sc.model(
            &format!("hot-{phase}.ink"),
            &format!(
                "(import \"std\")\n(def part (move (box 40 4 2) 20 0 0))\n(study s (shape part) (heat (k 50W/m/K))\n  (temperature cold (- x 1) 20degC)\n  (convection hot {} (coefficient 5000W/m2/K) (ambient 100degC)) (res 0.5))\n(require cool (max-temperature s < 30degC))\n",
                strips(phase)
            ),
        );
        for m in [roller, hot] {
            for a in sc.ink(["require", &m]).any()["requirements"].as_array().unwrap() {
                assert!(a["status"] == "fails" || a["status"] == "unknown", "{m}: {}", pretty(a));
            }
        }
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_fin_its_cells_cannot_see_decides_no_claim_either_way() {
    // A fin 0.05 mm thick, left thinner than two of its finest cells (1 mm
    // split to 0.125), is held too stiff: its tip reads 2.9e-10 mm, where
    // beam theory gives 0.00156 (and cells of 0.125 mm, 0.00155). `<`
    // held and `>` failed on that; each is unknown, saying why and to refine.
    let sc = Scratch::new("st-unseen-fin");
    let fin = sc.model(
        "fin.ink",
        "(import \"std\")\n(def part (min (move (box 4 6 6) 2 0 0) (move (box 20 4 0.05) 10 0 0.13)))\n\
         (study s (shape part) (elastic (E 200000MPa) (nu 0))\n (fix root (- x 1)) (load tip (force 0 0 -1e-5N) (region (- 19 x))) (res 1))\n\
         (require moves (max-displacement s > 0.0013))\n(require stiff (max-displacement s < 0.0013))\n",
    );
    let d = sc.ink(["require", &fin]).any();
    for a in d["requirements"].as_array().unwrap() {
        assert!(a["status"] == "unknown" && has(&a["undecided"], "cannot see the part") && has(&a["hint"], "refine"), "{}", pretty(a));
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_solve_that_may_not_fit_in_memory_says_so_and_answers_as_it_would() {
    // Before a solve allocates, it estimates its memory from its cells
    // (grid.memory_estimate_mb); where that is more than is available
    // (INK_MEMORY sets it) it says so on stderr and goes on, answering
    // as it would have. Here a plate too thin for coarser cells, whose
    // claim's grids go finer. (A cached solve allocates nothing, so none
    // is served here.)
    let s = Scratch::new("st-memory");
    let plate = s.model(
        "plate.ink",
        r#"(import "std")
(def part (box 20 4 0.6))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3)) (fix root (- (+ x 10) 1)) (load tip (force 0 0 -1N) (region (- 9 x))) (res 0.5))
(require stiff (max-displacement s < 1))
(require big (volume part > 10))
"#,
    );
    let run = |args: &[&str], budget: Option<&str>| {
        let i = s.ink(args).env("INK_NO_CACHE", "1");
        match budget {
            Some(b) => i.env("INK_MEMORY", b),
            None => i.env_remove("INK_MEMORY"),
        }
        .run()
    };
    let solved = run(&["solve", &plate, "--study", "s"], None).json();
    assert!(n(&solved["grid"]["memory_estimate_mb"]) > 0.0, "{}", pretty(&solved["grid"]));
    // The same claims and values, warned of or not; the estimate's own
    // sizes and a size that is none said so.
    let same = |a: &Value, b: &Value| {
        let (x, y) = (&a["requirements"], &b["requirements"]);
        x.as_array().unwrap().len() == y.as_array().unwrap().len()
            && x.as_array()
                .unwrap()
                .iter()
                .zip(y.as_array().unwrap())
                .all(|(p, q)| p["name"] == q["name"] && p["status"] == q["status"] && (n(&p["value"]) - n(&q["value"])).abs() <= 1e-9 * n(&p["value"]).abs())
    };
    let plain = run(&["require", &plate], None);
    for (budget, says) in [("10M", "it may run out of memory"), ("10MB", "it may run out of memory"), ("lots", "INK_MEMORY 'lots' is not a size")] {
        let tight = run(&["require", &plate], Some(budget));
        assert!(
            tight.stderr.contains(says) && !plain.stderr.contains(says) && same(&plain.json(), &tight.json()),
            "{budget}: {}\n{}\n{}",
            tight.stderr,
            pretty(&plain.json()),
            pretty(&tight.json())
        );
    }
}

// ---- peak stress: claims against their references ----

/// Each named claim of a model in tests/fixtures/stress-reads, and whether
/// its reference (a closed form, or CalculiX converged) makes it true.
const STRESS_REFERENCES: [(&str, &[(&str, bool)]); 7] = [
    // CalculiX 327.12 MPa at a d10 hole, on cells of 1 and of 2 mm.
    ("hole-split-r1.ink", &[("lt-310", false), ("lt-320", false), ("lt-335", true), ("lt-345", true), ("lt-360", true), ("gt-310", true), ("gt-335", false)]),
    ("hole-split-r2.ink", &[("lt-310", false), ("lt-320", false), ("lt-335", true), ("lt-345", true), ("lt-360", true), ("gt-310", true), ("gt-335", false)]),
    // CalculiX 23.8 to 24.7 MPa at a pocket's r1.5 fillets (yield 250).
    ("pocket.ink", &[("lt-22", false), ("lt-28", true), ("lt-33", true), ("gt-22", true), ("sf-gt-9", true), ("sf-gt-11", false), ("some", true), ("nonzero-safety", true)]),
    // CalculiX 2.956 MPa at a 1 mm bore whose field is far from a distance's.
    (
        "bore-skewed.ink",
        &[("lt-2-8", false), ("lt-3-0", true), ("lt-3-05", true), ("lt-3-1", true), ("lt-3-2", true), ("lt-3-4", true), ("gt-2-8", true), ("gt-3-1", false), ("sf-gt-80", true)],
    ),
    // Lame 51.05 MPa at a ring's bore, which its cut cells read high alike on every grid.
    ("vessel.ink", &[("vessel-peak-over-52", false), ("vessel-over-51-1", false), ("vessel-under-51-1", true)]),
    // CalculiX 2.27 MPa and rising at a flexure's r0.75 root fillets, three cells across.
    ("flexure.ink", &[("root-stress-under-2-24", false)]),
    // A rod's end at 1000 s, 44.19 C.
    ("transient-heat.ink", &[("true-lower", true)]),
];

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_peak_claim_near_its_reference_never_answers_against_it() {
    // A read can sit off the true peak in ways a claim's three grids
    // share, and each of these answered against its reference on some
    // build of ink: a ring's bore its cut cells read 4% high on every grid
    // (`> 52` held), a flexure's fillet read 14% low (`< 2.24` held), split
    // holes and a pocket's fillets read high (true claims failed), a rod's
    // transient end read low. No claim its reference makes true fails, and
    // none it makes false holds; unknown may be either.
    for (file, truths) in STRESS_REFERENCES {
        let d = ink(["require", &fixture(&format!("stress-reads/{file}"))]);
        let a = by(&d["requirements"], "name");
        assert_eq!(a.len(), truths.len(), "{file}: {}", pretty(&d["requirements"]));
        for (name, true_) in truths {
            let wrong = if *true_ { "fails" } else { "holds" };
            assert!(a[*name]["status"] != wrong, "{file} {name}: {}", pretty(&a[*name]));
        }
    }
}

#[test]
fn a_sharp_inside_corner_on_or_off_the_cells_planes_leaves_stress_claims_unknown() {
    // An L plate's re-entrant edge on the cells' planes, and the same
    // corner turned 17 degrees: singular, so every stress and
    // safety-factor claim at it is false or true whatever a grid reads
    // (a peak of 0.29 MPa read for `< 1`), and none is decided.
    for f in ["sharp-L.ink", "sharp-rotated.ink"] {
        let d = ink(["require", &fixture(&format!("stress-reads/{f}"))]);
        for a in d["requirements"].as_array().unwrap() {
            assert!(a["status"] == "unknown" && has(&a["undecided"], "sharp inside corner"), "{f}: {}", pretty(a));
        }
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_sharp_inside_corner_is_what_leaves_its_claims_unknown_however_far_its_clamp() {
    // The L of sharp-L.ink with 40 mm arms, clamped 19 mm from its corner,
    // the corner on the cells' planes or turned and moved off them, on 1
    // and 0.5 mm cells. A crease's peak was read with an unbounded
    // allowance (its radius spans no cells): the support margin's stress,
    // counted with it, answered every claim but `> 0.001` (unknown, but
    // said for the clamp, not the corner).
    let sc = Scratch::new("st-crease-far");
    for (size, turn, dx, dy, res) in [(40, 0, 0.0, 0.0, 1.0), (40, 17, 0.37, 0.61, 1.0), (40, 30, 0.25, 0.5, 0.5), (20, 17, 0.25, 0.5, 0.5)] {
        let (a, reach) = (size / 2, size / 2 - 1);
        let m = sc.model(
            &format!("L{size}-{turn}-{dx}-{res}.ink"),
            &format!(
                "(import \"std\")\n(def raw (subtract (box {size} {size} 4) (move (box {size} {size} 8) {a} {a} 0)))\n\
                 (def part (move (rotate-z raw {turn}) {dx} {dy} 0))\n(fn reg (f) (move (rotate-z f {turn}) {dx} {dy} 0))\n\
                 (study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa))\n (fix wall (reg (+ x {reach})))\n \
                 (load end (force 0 1N 0) (region (reg (- {reach} x)))) (res {res}))\n\
                 (require lt1 (max-stress s < 1MPa))\n(require lt1e4 (max-stress s < 10000MPa))\n(require gt1e-3 (max-stress s > 0.001MPa))\n\
                 (require sf-gt1 (min-safety-factor s > 1))\n(require sf-gt1e-3 (min-safety-factor s > 0.001))\n(require sf-lt1e6 (min-safety-factor s < 1000000))\n"
            ),
        );
        let d = sc.ink(["require", &m]).any();
        for a in d["requirements"].as_array().unwrap() {
            assert!(a["status"] == "unknown" && has(&a["undecided"], "sharp inside corner"), "{m}: {}", pretty(a));
        }
    }
}

#[test]
fn a_sharp_inside_corner_away_from_the_read_peak_holds_no_stress_bound() {
    // A stepped shaft whose keyway's floor meets its walls in sharp inside
    // corners: its elastic stress there is unbounded, so no `max-stress <`
    // claim on it is true. Its peak is read at the shoulder fillet (160.7
    // on 1 mm cells, CalculiX 151.9), where alone the surface was surveyed
    // for a crease; the corners, some 67 mm along it, were not, and `< 1e6`
    // held. Its whole read surface surveyed, the claim is unknown, naming
    // a corner; leaving them out on purpose, it is on what the cells read.
    let d = ink(["require", &fixture("stress-reads/shaft.ink")]);
    let a = by(&d["requirements"], "name");
    assert!(a["shoulder-peak-band"]["status"] == "unknown" && has(&a["shoulder-peak-band"]["undecided"], "leaves out a sharp inside corner"), "{}", pretty(&d["requirements"]));
    assert!(a["shoulder-peak-over-155"]["status"] != "fails", "{}", pretty(&d["requirements"]));
    let sc = Scratch::new("st-shaft-corners");
    let own = std::fs::read_to_string(fixture("stress-reads/shaft.ink")).unwrap().replace("(max-stress drive < 1e6MPa)", "(max-stress drive < 1e6MPa (excluding sharp-corners))");
    let d = sc.ink(["require", &sc.model("shaft.ink", &own)]).any();
    let x = &by(&d["requirements"], "name")["shoulder-peak-band"];
    assert!(x["status"] == "holds" && pretty(x).contains("on purpose (excluding sharp-corners)"), "{}", pretty(x));
    // What counting the corner would do, said as what it is: no number, an unbounded stress (a reader took
    // `Including its value (0 at ...)` for a measured 0).
    assert!(
        pretty(x).contains("wherever load passes it the elastic stress there is unbounded") && pretty(x).contains("Counted, it is unbounded, and the verdict is fails"),
        "{}",
        pretty(x)
    );
    assert!(!pretty(x).contains("Including its value"), "{}", pretty(x));
}

#[test]
fn a_face_beside_a_thin_gap_holds_no_safety_factor_its_plates_do_not_have() {
    // Two plates across a 0.125 mm gap, each bent about its own mid-plane:
    // a safety factor of 31.25 exactly, not glued (their displacements are
    // the closed form's). The cells beside the gap read the upper plate's
    // edge 41% high, and the three grids converge away from the truth
    // (22.15, 24.72, 28.71): `< 30.625` held and `> 30.625` failed.
    let d = ink(["require", &fixture("stress-reads/gap-face.ink")]);
    let a = by(&d["requirements"], "name");
    assert!(a["gap-sf-under-30-6"]["status"] != "holds" && a["gap-sf-over-30-6"]["status"] != "fails", "{}", pretty(&d["requirements"]));
}

#[test]
fn a_heat_interface_turned_off_the_planes_holds_no_bound_its_bar_breaks() {
    // The 1-D heat bar turned 2 degrees: 2539.855 C exactly, read 2530.91
    // with a grid error of 6.55, so `< 2538.585` held, false, on the unsafe
    // side, and `> 2538.585` failed, true. Its insulating face is a stair
    // of cells on every grid: its claims are unknown, saying so.
    let d = ink(["require", &fixture("stress-reads/heat-turned.ink")]);
    let a = by(&d["requirements"], "name");
    for x in a.values() {
        assert!(x["status"] == "unknown" && has(&x["undecided"], "not flat along the cells"), "{}", pretty(x));
    }
}

#[test]
fn a_band_turned_off_the_planes_holds_no_stress_its_bar_does_not_have() {
    // A soft band in a pulled bar turned 1 degree: 2.451 MPa everywhere, a
    // safety factor of 102. Its stair of split cells read 88.3 MPa (grid
    // error 11.4): `> 2.598` held and the factor's `> 95.88` failed. Its
    // claims are unknown, saying its face is a stair of cells.
    let d = ink(["require", &fixture("stress-reads/band-turned.ink")]);
    let a = by(&d["requirements"], "name");
    for x in a.values() {
        assert!(x["status"] == "unknown" && has(&x["undecided"], "not flat along the cells"), "{}", pretty(x));
    }
}

#[test]
fn a_safety_factor_beside_a_band_reads_the_bars() {
    // A soft band in a pulled bar on 0.3 mm cells: a factor of 102 exactly.
    // A full cell on its face took its centre's stiffness and a cut cell
    // beside it, at the bar's surface, its samples' mean: the surface
    // carried the load there, and the factor read 49.8, 61.7 and 86.4 on
    // 0.3, 0.6 and 1.2 mm cells: `< 95.88` held and `> 99.96` failed. Each
    // cell on the face its two materials in series, it reads 102.0.
    let d = ink(["require", &fixture("stress-reads/band-beside.ink")]);
    let a = by(&d["requirements"], "name");
    assert!(a["band-sf-under-95-9"]["status"] != "holds" && a["band-sf-over-100"]["status"] != "fails", "{}", pretty(&d["requirements"]));
    assert!((n(&a["band-sf-over-100"]["value"]) / 102.0 - 1.0).abs() < 0.02, "{}", pretty(&d["requirements"]));
}

#[test]
fn study_material_limits_follow_the_values() {
    let s = Scratch::new("material-values");
    // Pure tension, 0.625 MPa: its end extends 0.625 x (10.5 + 0.5) / E.
    // A clamp one cell thick on the finest grid (a thicker one's penalty
    // compresses its outer cells), so no cell is compressed and buckling
    // is ruled out on every grid, whatever the E.
    let bar = "(study pull (shape (box 12 4 4)) (elastic (E modulus) (nu 0)) (fix wall (+ x 5.5)) (load tip (force 10N 0 0) (region (- 5 x))) (res 1))\n(require stiff (max-displacement pull < 0.000001))";
    let model = |param: &str| s.model("bar.ink", &format!("(import \"std\") (param modulus {param}) {bar}"));
    message_has(&s.ink(["require", &model("200e9")]).err()["error"]["message"], &["E takes a pressure or stress, and modulus is a bare number"]);
    let m = model("200e9MPa");
    let error = s.ink(["require", &m]).err();
    message_has(&error["error"]["message"], &["modulus=200000000000", "E (MPa)", "0.0001..1200000"]);
    let m = model("200GPa");
    let d = s.ink(["require", &m]).any();
    assert_eq!(d["requirements"][0]["status"], "fails");
    assert!((n(&d["requirements"][0]["value"]) - 3.4375e-5).abs() < 1e-7);
    let d = s.ink(["require", &m, "--param", "modulus=200"]).err();
    message_has(&d["error"]["message"], &["modulus", "200GPa", "modulus=200"]);
    let d = s.ink(["require", &m, "--param", "modulus=100GPa:1300GPa"]).any();
    assert_eq!(d["requirements"][0]["status"], "unknown");
    message_has(&d["requirements"][0]["unavailable"], &["modulus=100000..1300000", "E (MPa)", "0.0001..1200000"]);
    let d = s.ink(["require", &model("1e6Pa")]).any();
    assert_eq!(d["requirements"][0]["status"], "fails");
    assert!((n(&d["requirements"][0]["value"]) - 6.875).abs() < 0.01);
    let bare = s.model("soft.ink", &format!("(import \"std\") (param modulus 1e6MPa) {}", bar.replace("0.000001", "0.001")));
    let d = s.ink(["require", &bare]).any();
    assert_eq!(d["requirements"][0]["status"], "holds");
    assert!((n(&d["requirements"][0]["value"]) - 6.875e-6).abs() < 1e-8);
    let explicit = s.model("soft.ink", &format!("(import \"std\") (param modulus 1e6Pa) {}", bar.replace("0.000001", "0.001")));
    assert_eq!(s.ink(["require", &explicit]).any()["requirements"][0]["status"], "fails");
    let m = s.model(
        "yield.ink",
        &format!(
            "(import \"std\") (param modulus 200GPa) (param strength 250e6) {} (require safe (min-safety-factor pull > 2))",
            bar.replace("(nu 0)", "(nu 0)\n(yield (* 1MPa strength))")
        ),
    );
    let error = s.ink(["require", &m]).err();
    message_has(&error["error"]["message"], &["strength=250000000", "yield strength (MPa)", "0.001..100000", "250e6Pa"]);
    assert_eq!(error["error"]["at"]["line"], 2);
    for (clause, value, why) in [("(density bad)", "2700", "density"), ("(nu bad)", "0.5", "nu must be"), ("(nu bad)", "0.2mm", "dimensionless")] {
        let m = s.model("invalid.ink", &format!("(import \"std\") (param modulus 200GPa) (param bad {value}) {}", bar.replace("(nu 0)", clause)));
        assert!(has(&s.ink(["check", &m]).err()["error"]["message"], why));
    }
}

fn message_has(value: &Value, expected: &[&str]) {
    let message = value.as_str().expect("diagnostic message");
    for text in expected {
        assert!(message.contains(text), "expected {text:?} in {message}");
    }
}

#[test]
fn material_checks_weigh_a_published_mass_only_when_used() {
    let s = Scratch::new("published-material");
    for density in ["(attr-of ballast density)", "(by-body ballast density)"] {
        let m = s.model(
            "part.ink",
            &format!(
                "(import \"std\") (body ballast (shape (box 10 10 10)) (attr mass 25000g))
            (study pull (shape (shape-of ballast)) (elastic (E 200GPa) (nu 0) (density {density})) (fix wall (+ x 4)) (load tip (force 10N 0 0) (region (- 4 x))) (res 1))"
            ),
        );
        let d = s.ink(["check", &m, "--load"]).ok();
        assert!(d["params"].as_array().unwrap().iter().any(|p| p["name"] == "ballast.mass-density" && p["value"].is_null()), "{}", pretty(&d));
        for cmd in ["mass", "solve"] {
            let d = s.ink([cmd, &m]).err();
            message_has(&d["error"]["message"], &["ballast.mass-density", "25000 g"]);
        }
    }
}

#[test]
fn heat_checks_the_material_input_in_the_claims_scope() {
    let s = Scratch::new("transport-units");
    let heat = "(study h (shape (box 8 4 4)) (heat (k prop)) (temperature wall (+ x 3) 20degC) (source heater (power 0.01W) (region (- 3 x))) (res 1))";
    let m = s.model("bad.ink", &format!("(import \"std\") (param prop 1000000000W/m/K) {heat}"));
    let d = s.ink(["check", &m]).err();
    message_has(&d["error"]["message"], &["prop=1000000000", "W/m/K", "0.00001..100000"]);
    let form = heat.replace("(k prop)", "(k 10W/m/K) (capacity prop)");
    for value in ["0", "0.0005", "4184000"] {
        let m = s.model("capacity.ink", &format!("(import \"std\") (param prop {value}J/kg/K) {form}"));
        let d = s.ink(["check", &m]).err();
        message_has(&d["error"]["message"], &[&format!("prop={value}"), "0.001..20000", "J/kg/K"]);
    }
    let m = s.model("capacity.ink", &format!("(import \"std\") (param prop 4184J/kg/K) {form}"));
    s.ink(["check", &m]).ok();
    let m = s.model(
        "capacity-band.ink",
        &format!(
            "(import \"std\") (param prop 4184J/kg/K) {} (require temp (max-temperature h < 100degC)) (require width (value 8 < 7))",
            heat.replace("(k prop)", "(k 10W/m/K) (capacity prop)")
        ),
    );
    let d = s.ink(["require", &m, "--param", "prop=4000:40000"]).err();
    message_has(&d["error"]["message"], &["prop"]);
    let m = s.model("transient.ink", &std::fs::read_to_string(fixture("stress-reads/transient-heat.ink")).unwrap().replace("(capacity 500J/kg/K)", "(capacity 500000J/kg/K)"));
    message_has(&s.ink(["require", &m]).err()["error"]["message"], &["500000", "0.001..20000"]);
    let m = s.model("band.ink", &format!("(import \"std\") (param prop 0.1W/m/K) {heat} (require temp (max-temperature h < 100degC)) (require width (value 8 < 7))"));
    let d = s.ink(["require", &m, "--param", "prop=1W/m/K:1000000W/m/K"]).any();
    assert_eq!(d["requirements"][0]["status"], "unknown");
    message_has(&d["requirements"][0]["unavailable"], &["prop=1..1000000", "0.00001..100000", "W/m/K"]);
    assert_eq!(d["requirements"][1]["status"], "fails");
}

#[test]
fn a_study_loaded_from_the_model_is_read_at_the_models_steps_not_at_even_ones() {
    // A lever whose reaction swings 64 times in the span loads a bar. Read
    // at 64 even moments it was zero at every one (2.3e-14 mm, a yes); the
    // model's own steps follow its loads, and the bar's worst is 0.575 mm.
    let sc = Scratch::new("study-alias");
    let src = r#"(import "std")
(gravity 0 0 0)
(t-span 0 1s)
(body joint (shape (box 4 4 4)) (attr material "steel"))
(assembly root (instance lever joint (at) (turn z 0)))
(load wind root/lever (force 0 (* 100N (sin (* 23040deg/s t))) 0))
(def wind (reaction root/lever))
(study bar (shape (move (box 100 10 10) 50 0 0)) (elastic (material aluminum-6061)) (fix wall (- x 0.25))
  (load tip (force 0 (nth wind 1) 0) (region (- 99.5 x))) (res 2))
(require stiff (max-displacement bar < 0.3mm))
"#;
    let d = sc.ink(["require", &sc.model("alias.ink", src)]).any();
    let a = &d["requirements"][0];
    assert!(verdict(a) == "fails estimate" && (0.5..0.65).contains(&a["value"].as_f64().unwrap_or(0.0)), "{}", pretty(&d));
    // A hundred times faster the steps run out before they follow it: not
    // settled, saying so, never a yes at the moments it did read.
    let d = sc.ink(["require", &sc.model("faster.ink", &src.replace("23040deg/s", "2304000deg/s"))]).any();
    let a = &d["requirements"][0];
    assert!(verdict(a) == "unknown" && pretty(a).contains("not settled between the moments read"), "{}", pretty(&d));
    assert!(d["all_yes"] == false, "{}", pretty(&d));
}

#[test]
fn an_estimates_condition_is_shown_over_the_clock_it_is_read_on() {
    // A transient reads its power on its own clock, 0 to 1 s, and the
    // estimate applies before 0.25 s only: whatever the model's moment, the
    // study has no answer. It was a yes at t = 0 and unknown at t = 0.5.
    let sc = Scratch::new("study-estimate-clock");
    let src = r#"(import "std")
(study heater (shape (box 8 8 2))
  (heat (k 20W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))
  (convection air -1 (coefficient 10W/m2/K) (ambient 20degC))
  (source heat (power (estimate-when (* 1W (+ 1 (/ t 1s))) (< t LIMIT) "startup only")) (region -1))
  (transient (duration 1s) (steps 2) (initial 20degC)) (res 2))
(require warm (max-temperature heater < 500degC))
"#;
    let model = sc.model("heater.ink", &src.replace("LIMIT", ".25s"));
    for t in ["0", "0.1", "0.5"] {
        let d = sc.ink(["require", &model, "--t", t]).any();
        assert!(statuses(&d)["warm"] == "unknown" && pretty(&d).contains("startup only"), "t {t}: {}", pretty(&d));
    }
    // One that applies over the whole clock is run, at a moment inside it.
    let model = sc.model("whole.ink", &src.replace("LIMIT", "2s"));
    let d = sc.ink(["require", &model, "--t", "0.5"]).any();
    assert_eq!(statuses(&d)["warm"], "holds estimate", "{}", pretty(&d));
}

/// A study answers for all of its shape, so a shape not proven bounded is refused: solved in the window
/// alone, a solid no fix holds, 1000 mm off, was left out, and the held box answered yes for both.
#[test]
fn a_study_of_a_shape_not_proven_bounded_is_refused() {
    let s = Scratch::new("study-unbounded");
    let f = s.model(
        "far.ink",
        "(import \"std\")\n(def remote (move (twist-z (box 20 20 40) 2) 1000 0 0))\n(body part (shape (min (box 20 8 8) remote)))\n\
         (study beam (shape (shape-of part)) (elastic (E 2000MPa) (nu 0.3) (density 1000kg/m3)) (gravity 0 0 -9.81m/s2) (fix base (+ x 8)) (res 2))\n\
         (require stiff (max-displacement beam < 1))",
    );
    let r = s.ink(["require", &f]).any();
    let a = first(&r);
    assert_eq!(r["all_yes"], false, "{}", pretty(&r));
    assert_eq!(a["status"], "unknown", "{}", pretty(&a));
    message_has(&a["unavailable"], &["study 'beam'", "not proven bounded", "intersect it with a containing shape"]);
    message_has(&s.ink(["solve", &f]).err()["error"]["message"], &["study 'beam'", "not proven bounded", "intersect it with a containing shape"]);
}

/// A hold along some directions only may be a physical support or a cut of symmetry, and only as a
/// support does the part as cut hold every buckling mode. One that is a whole flat face is solved both
/// ways; any other that does not say what it is leaves its static claims not settled, saying so, and
/// `(support)` settles them.
#[test]
fn a_partial_hold_that_does_not_say_what_it_is_leaves_its_claim_unsettled() {
    let s = Scratch::new("study-unstated-hold");
    let bar = "(import \"std\")\n(study s (shape (box 40 6 6)) (elastic (E 70000MPa) (nu 0.33)) (fix wall (+ x 19.5))\n\
               (fix guide (max (- 15 x) (- y -2)) (dofs y) SUPPORT) (load tip (force 0 0 -20N) (region (- 19.5 x))) (res 1))\n\
               (require stiff (max-displacement s < 5mm))";
    let ask = |name: &str, support: &str| s.ink(["require", &s.model(name, &bar.replace("SUPPORT", support))]).any()["requirements"][0].clone();
    let unsaid = ask("unsaid.ink", "");
    assert!(unsaid["status"] == "unknown" && has(&unsaid["undecided"], "does not say what it is") && has(&unsaid["undecided"], "guide"), "{}", pretty(&unsaid));
    let said = ask("said.ink", "(support)");
    assert!(said["status"] == "holds" && !said.to_string().contains("does not say what it is"), "{}", pretty(&said));
}

/// A study's von Mises field is each cell's stress at its centre: it holds no peak on the surface, by a
/// support or at a sharp corner. A claim that reads it for strength is not settled, and says which
/// metrics measure those.
#[test]
fn a_claim_on_a_studys_stress_field_is_not_settled() {
    let s = Scratch::new("study-field-strength");
    let f = s.model(
        "field.ink",
        "(import \"std\")\n(def bar (move (box 20 4 4) 10 0 0))\n\
         (study s (shape bar) (elastic (E 2000MPa) (nu 0.3)) (fix wall (- x 0.5)) (load tip (force 0 0 -10N) (region (- 19.5 x))) (res 1))\n\
         (def stress (study-field s von-mises))\n(require under (subset bar (- stress 1000)))",
    );
    let a = s.ink(["require", &f]).any()["requirements"][0].clone();
    assert!(a["status"] == "unknown" && a.to_string().contains("each cell's stress at its centre") && a.to_string().contains("max-stress"), "{}", pretty(&a));
}

/// A wall thinner than half a cell is read: its cells are split until it is two of them across, so its
/// peak on coarse cells is the fine cells' within their own difference, not the block's beneath it.
/// (A wall 0.8 mm thick on a block, pushed sideways at its top.)
#[test]
fn a_wall_thinner_than_half_a_cell_has_its_stress_read() {
    let s = Scratch::new("study-thin-wall");
    let wall = "(import \"std\")\n(def part (min (move (box 20 20 10) 0 0 5) (move (box 20 0.8 10) 0 0.3 15)))\n\
                (study s (shape part) (elastic (E 2000MPa) (nu 0.3)) (fix base (- z 0.5)) (load push (force 0 2N 0) (region (- 19 z))) (res RES))\n\
                (require strong (max-stress s < 1000MPa (excluding sharp-corners)))";
    let peak = |res: &str| n(&s.ink(["require", &s.model(&format!("wall-{res}.ink"), &wall.replace("RES", res))]).any()["requirements"][0]["value"]);
    let (coarse, fine) = (peak("2"), peak("0.25"));
    assert!((coarse - fine).abs() < 0.2 * fine && fine > 9.0, "{coarse} on 2 mm cells against {fine} on 0.25 mm");
}

/// A transient's power that passes through an undefined case at a moment it reads (the language takes
/// a negative base to the power 0.5 as 0, from 0.5 s on) has no answer there: the study is refused,
/// saying which source, what and when (it answered yes, at 1 W there). A square root of a negative is no such case: it is defined, the root of 0, and a power
/// written so is that (the case this was filed with, `(sqrt -1)` from 0.25 s on, a power of 0 there).
#[test]
fn a_transient_refuses_a_power_with_no_value_at_a_moment_it_reads() {
    let sc = Scratch::new("study-nan-power");
    let src = r#"(import "std")
(study heater (shape (box 8 8 2))
  (heat (k 20W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))
  (convection air -1 (coefficient 10W/m2/K) (ambient 20degC))
  (source heat (power (* 1W (+ 1 (pow (- 0.5 (/ t 1s)) 0.5)))) (region -1))
  (transient (duration 1s) (steps 2) (initial 20degC)) (res 2))
(require warm (max-temperature heater < 500degC))
"#;
    let d = sc.ink(["require", &sc.model("heater.ink", src)]).any();
    let said = d.to_string();
    assert!(said.contains("source 'heat'") && said.contains("a negative base to a power that is not a whole number at t = 0.75"), "{}", pretty(&d));
}

/// A pressure on a strip of a face narrower than its facets takes the strip's area: each facet's share
/// is the area of it the region holds, its edge found by halving the facet, not the share of four
/// points on it (a strip 0.1 mm wide on cells of 0.5 took 6.25% short, one 0.21 wide 4.2% over, and one
/// 0.03 wide was refused). A 40 x 4 x 2 bar, 0.001 MPa on its top face where y lies in [b, b + w], along
/// the 39 mm the region spans: 39 w x 0.001 N.
#[test]
fn a_pressure_on_a_strip_narrower_than_a_facet_takes_its_area() {
    let sc = Scratch::new("strip-load");
    for (b, w) in [(0.1, 0.1), (0.0, 0.1), (0.05, 0.03), (0.37, 0.21), (1.0, 2.0)] {
        let f = sc.model(
            &format!("strip-{b}-{w}.ink"),
            &format!(
                "(import \"std\")\n(study s (shape (move (box 40 4 2) 20 2 1)) (elastic (E 2000MPa) (nu 0.3)) (fix wall (- x 0.5))\n\
                 (load p (pressure 0.001MPa) (region (max (- {b} y) (- y {}) (- 1.9 z) (- 1 x)))) (res 0.5))",
                b + w
            ),
        );
        let d = ink(["solve", &f]);
        let got = -n(&d["loads"][0]["total"][2]);
        let meant = 39.0 * w * 0.001;
        assert!((got / meant - 1.0).abs() < 0.005, "b {b} w {w}: {got} against {meant}");
    }
}

/// A shape is its sign: a bore written as the cube of a distance is the same part as the distance, and
/// is solved as it. The cells place the surface by the field's value over its slope, a third of the depth
/// for a cube, so samples in the bore's air took material; where the field's bounds over a sample's box
/// prove it all air or all material it is so, and a chord through a field whose slope changes is taken
/// along the surface's normal. It was split as thin to 0.125 mm on 1 mm cells and read 10% under the
/// distance's stress; now it splits nothing and reads within 0.1% of it. (A bore written (1 - r)(1 + 0.95
/// sin 25x) no longer reads a 0.215 mm gap it does not have.)
#[test]
fn a_bore_written_as_the_cube_of_a_distance_is_the_same_part() {
    let sc = Scratch::new("bore-cube");
    let plate = |bore: &str| {
        let f = sc.model(
            &format!("p{}.ink", bore.len()),
            &format!(
                "(import \"std\")\n(def part (max (box 30 20 4) {bore}))\n(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa))\n \
                 (fix wall (+ x 14.5)) (load pull (force 160N 0 0) (region (- 14.5 x))) (res 0.5))"
            ),
        );
        let d = ink(["solve", &f]);
        (n(&d["grid"]["cells_refined"]), n(&d["results"]["max_von_mises_mpa"]["value"]))
    };
    let (cube, distance) = (plate("(pow (- 1.5 (sqrt (+ (square x) (square y)))) 3)"), plate("(- 1.5 (sqrt (+ (square x) (square y))))"));
    assert!(cube.0 == distance.0 && (cube.1 / distance.1 - 1.0).abs() < 1e-3, "cube {cube:?}, distance {distance:?}");
    let skewed = ink(["solve", &fixture("stress-reads/bore-skewed.ink")]);
    assert!(skewed["grid"]["glued"].is_null() && !skewed.to_string().contains("glue shut"), "{}", pretty(&skewed["warnings"]));
}

/// A region drawn as a slab a micron under a face holds that face whole: the face's points are found to
/// the rounding of their place, where they were 6e-5 of a cell off it either way. NAFEMS T4 (benchmarks/
/// standard, its faces written so) reads 18.249 C at (600, 200, 0) against the exact 18.2538 (it read
/// 82.4 with almost no cooling, and since facets take a region by area it was refused as touching none);
/// a tube's end pressed through such a slab takes what its face takes (it took 35% short).
#[test]
fn a_region_a_micron_under_a_face_holds_the_face() {
    let sc = Scratch::new("region-sliver");
    let t4 =
        std::fs::read_to_string(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../benchmarks/standard/nafems-t4.ink")).unwrap() + "\n(def T (study-field s temperature))\n";
    let d = sc.ink(["eval", &sc.model("t4.ink", &t4), "--target", "T", "--points", "[[600,200,0]]"]).ok();
    let t = n(&d["targets"][0]["results"]["rows"][0][1]);
    assert!((t / 18.253756537134212 - 1.0).abs() < 1e-3, "{t}");
    let tube = |region: &str| {
        let f = sc.model(
            &format!("tube-{}.ink", region.len()),
            &format!(
                "(import \"std\")\n(def tube (subtract (cylinder-z 10 16) (cylinder-z 6 20)))\n(study s (shape tube) (elastic (E 200000MPa) (nu 0))\n  (fix root (+ z 7)) (load end (pressure -1MPa) (region {region})) (res 1))"
            ),
        );
        n(&ink(["solve", &f])["loads"][0]["total"][2])
    };
    let (face, slab) = (tube("(face tube (normal 0 0 1))"), tube("(- 7.999999 z)"));
    assert!((slab / face - 1.0).abs() < 1e-4 && (face / (64.0 * std::f64::consts::PI) - 1.0).abs() < 0.002, "{slab} and {face}");
}

#[test]
fn a_pinned_column_buckles_where_euler_says() {
    // Steel, 204 x 7.2 x 4.8, under 100 N along it, pinned at both ends (`pivot` about y, the far end also
    // `slide` along x) or clamped at one: Euler's factor with shear (Engesser's, Cowper's kappa), against Poisson's
    // ratio. Each pin held its end face rigid but for its turn, and the bar could not spread there as it was
    // squeezed: pinned-pinned read 31.256 at nu 0.3 (its band 31.108 to 31.278), 31.14 at 0.45; nu 0 read right.
    let sc = Scratch::new("st-pinned-column");
    let (e, l, b, t) = (200000.0, 204.0, 7.2, 4.8);
    let truth = |nu: f64, kl: f64| {
        let (i, a, g, kappa) = (b * t * t * t / 12.0, b * t, e / (2.0 * (1.0 + nu)), 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu));
        let pe = kl * kl * e * i / (l * l);
        pe / (1.0 + pe / (kappa * g * a)) / 100.0
    };
    let mut text = String::from("(import \"std\")\n(def part (move (box 204 7.2 4.8) 102 0 0))\n");
    let mut cases = Vec::new();
    for (end, kl) in [("pp", std::f64::consts::PI), ("fp", 4.493409457909064)] {
        for nu in [0.0, 0.3, 0.45] {
            let name = format!("{end}{}", (nu * 100.0) as u32);
            let near = if end == "pp" { "(pivot (at 0 0 0) (axis 0 1 0))" } else { "" };
            text += &format!(
                "(study {name} (shape part) (elastic (E 200000MPa) (nu {nu}))\n  (fix pin (- x 1e-06) {near})\n  (fix roller (- 203.999999 x) (pivot (at 204 0 0) (axis 0 1 0)) (slide 1 0 0))\n  (load push (force -100N 0N 0N) (region (- 203.999999 x)))\n  (res 0.6))\n(require {name} (min-buckling-factor {name} > 0))\n"
            );
            cases.push((name, truth(nu, kl)));
        }
    }
    let d = sc.ink(["require", &sc.model("columns.ink", &text)]).any();
    for (name, want) in cases {
        let a = d["requirements"].as_array().and_then(|r| r.iter().find(|r| r["name"] == name.as_str())).unwrap_or_else(|| panic!("{name}: {}", pretty(&d)));
        let (x, lo, hi) = (n(&a["value"]), n(&a["band"]["low"]), n(&a["band"]["high"]));
        eprintln!("{name}: {x} (band {lo} to {hi}) against {want}");
        // (Within 0.1%, not inside its band: the beam's formula misses the solid's by up to 0.3% (benchmarks/holds
        // solid.py); CalculiX's solid reads the pinned column at nu 0.45 31.4141, 0.03% under the beam's. Fixed-pinned
        // at nu 0.45 within 0.5%: its clamped face, which cannot contract, makes the solid stiffer than the beam.)
        let within = if name == "fp45" { 5e-3 } else { 1e-3 };
        assert!((x - want).abs() < within * want && lo < x && x < hi, "{name}: {x} (band {lo} to {hi}) against {want}: {}", pretty(a));
    }
}

#[test]
fn a_pulse_ending_between_steps_is_not_read_low() {
    // A 7.3 mm copper cube (Biot 0.009, lumped), 5 W for 3.3 s, cooled at 3000 W/m2/K to 20 C: its peak is
    // 20 + P/(hA) (1 - exp(-3.3/tau)) = 24.720 C at 3.3 s, the 3-D peak no lower. On 23 steps (3.3 s inside the
    // 8th) it read 24.52, its band 24.39 to 24.65 from 12 and 6 steps, and `< 24.68` held.
    let sc = Scratch::new("st-pulse-between-steps");
    let f = sc.model(
        "pulse.ink",
        r#"(import "std")
(def cube (move (box 7.3 7.3 7.3) 1.19 -2.03 0.77))
(study s (shape cube) (heat (k 400W/m/K) (density 8960kg/m3) (capacity 385J/kg/K))
  (convection air (face cube) (coefficient 3000W/m2/K) (ambient 20degC))
  (source heater (power (* 5W (select (- t 3.3s) 1 0))) (region -1))
  (transient (duration 10s) (steps 23) (initial 20degC))
  (res 0.7))
(require q (max-temperature s < 24.68degC))"#,
    );
    let a = first(&sc.ink(["require", &f]).any());
    assert_ne!(a["status"], "holds", "{}", pretty(&a));
}

#[test]
fn a_study_grid_error_hint_names_res_inside_the_study_form() {
    let sc = Scratch::new("st-hint-study-res");
    // A thin cantilever on coarse cells (res 4) whose displacement is close enough to its bound
    // that the three-grid error leaves it undecided: its hint must tell the author to change
    // `(res ...)` inside `(study hang ...)`, not run `ink require --res`.
    let f = sc.model(
        "hang.ink",
        r#"(import "std")
(study hang (shape (move (box 60 8 4) 30 0 0)) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (- x 1)) (load tip (force 0 0 -10N) (region (- 59 x))) (res 4))
(require stiff (max-displacement hang < 7.9mm))"#,
    );
    let a = first(&sc.ink(["require", &f]).any());
    assert_eq!(a["status"], "unknown", "{}", pretty(&a));
    assert!(has(&a["hint"], "(res ") && has(&a["hint"], "in (study hang ...)"), "{}", pretty(&a));
    // Passing `--res` to `ink require` does not change the study's own cell size.
    let with_flag = first(&sc.ink(["require", &f, "--res", "1"]).any());
    assert_eq!(with_flag["status"], "unknown", "{}", pretty(&with_flag));
    assert_eq!(with_flag["value"], a["value"], "{}", pretty(&with_flag));
}

#[test]
fn a_study_setup_failure_still_reports_independent_assertions() {
    let sc = Scratch::new("st-setup-report");
    let f = sc.model(
        "outside.ink",
        r#"(import "std")
(study h (shape (box 10 10 10)) (heat (k 20W/m/K))
  (temperature base (+ z 4) 20degC)
  (source heater (power 5W) (region (- 50 x))) (res 2))
(require hot (max-temperature h < 100degC))
(require width (value 8 < 7))"#,
    );
    let r = sc.ink(["require", &f]).any();
    assert_eq!(r["all_yes"], false, "{}", pretty(&r));
    let reqs = r["requirements"].as_array().unwrap();
    assert_eq!(reqs.len(), 2, "{}", pretty(&r));
    let hot = reqs.iter().find(|x| x["name"] == "hot").unwrap();
    let width = reqs.iter().find(|x| x["name"] == "width").unwrap();
    assert_eq!(hot["status"], "unknown", "{}", pretty(hot));
    assert!(has(&hot["unavailable"], "heater"), "{}", pretty(hot));
    assert_eq!(width["status"], "fails", "{}", pretty(width));
    let e = sc.ink(["solve", &f]).err();
    assert!(has(&e["error"]["message"], "heater"), "{}", pretty(&e));
}

#[test]
fn an_elastic_study_balances_moments_as_well_as_forces() {
    let sc = Scratch::new("st-moment-balance");
    // 1) Cantilever with an off-origin tip force: M = r x F about the origin.
    let bend = sc.model(
        "bend.ink",
        r#"(import "std")
(study s (shape (move (box 40 6 6) 20 0 0)) (elastic (E 70000MPa) (nu 0.3))
  (fix wall (- x 1)) (load tip (force 0 0 -10N) (region (- 39 x))) (res 2))"#,
    );
    let d = ink(["solve", &bend]);
    let lm = &d["loads"][0]["moment_nmm"];
    let rm = &d["reactions"][0]["moment_nmm"];
    assert!(
        n(&d["moment_imbalance"]) < 1e-6 && n(&d["equilibrium_imbalance"]) < 1e-6 && (n(&lm[1]) - 395.0).abs() < 15.0 && (n(&lm[1]) + n(&rm[1])).abs() < 1e-4,
        "{}",
        pretty(&json!({"imbalance": d["equilibrium_imbalance"], "moment_imbalance": d["moment_imbalance"], "loads": d["loads"], "reactions": d["reactions"]}))
    );
    // 2) Pure torque load: zero net force, non-zero net moment.
    let twist = sc.model(
        "twist.ink",
        r#"(import "std")
(study s (shape (move (box 40 8 8) 20 0 0)) (elastic (E 70000MPa) (nu 0.3))
  (fix wall (- x 1)) (load t (torque 500N*mm 0 0) (region (- 39 x))) (res 2))"#,
    );
    let d = ink(["solve", &twist]);
    assert!(
        n(&d["moment_imbalance"]) < 1e-6
            && n(&d["equilibrium_imbalance"]) < 1e-6
            && (n(&d["loads"][0]["moment_nmm"][0]) - 500.0).abs() < 1e-3
            && (n(&d["reactions"][0]["moment_nmm"][0]) + 500.0).abs() < 1e-3,
        "{}",
        pretty(&json!({"imbalance": d["equilibrium_imbalance"], "moment_imbalance": d["moment_imbalance"], "loads": d["loads"], "reactions": d["reactions"]}))
    );
    // 3) Symmetric axial tension through the origin: zero applied moment must not divide noise by noise.
    let pull = sc.model(
        "pull.ink",
        r#"(import "std")
(study s (shape (box 40 6 6)) (elastic (E 70000MPa) (nu 0.3))
  (fix wall (+ x 19)) (load end (force 100N 0 0) (region (- x 19))) (res 2))"#,
    );
    let d = ink(["solve", &pull]);
    assert!(
        n(&d["moment_imbalance"]) < 1e-6 && n(&d["equilibrium_imbalance"]) < 1e-6,
        "{}",
        pretty(&json!({"imbalance": d["equilibrium_imbalance"], "moment_imbalance": d["moment_imbalance"], "loads": d["loads"], "reactions": d["reactions"]}))
    );
}

#[test]
fn a_cantilevers_balance_is_independent_of_its_distance_from_the_origin() {
    let sc = Scratch::new("st-translated-balance");
    let solve_at = |x: i32, name: &str| {
        let model = format!(
            "(import \"std\")
(study s (shape (move (box 20 4 4) {x} 0 0)) (elastic (E 70000MPa) (nu 0.3))
   (fix wall (- x {})) (load tip (force 0 0 -10N) (region (- {} x))) (res 1))",
            x - 9,
            x + 9
        );
        let f = sc.model(name, &model);
        ink(["solve", &f])
    };
    let (near, far) = (solve_at(0, "near.ink"), solve_at(1000, "far.ink"));
    for key in ["equilibrium_imbalance", "moment_imbalance"] {
        let (a, b) = (n(&near[key]), n(&far[key]));
        assert!((a - b).abs() < 0.02 * a.max(b) + 2e-12, "{key}: {a:e} at the origin, {b:e} at x = 1000");
    }
    // Reported load moments still use the model's origin.
    let moment = |r: &Value| n(&r["loads"][0]["moment_nmm"][1]);
    assert!((moment(&far) - moment(&near) - 10_000.0).abs() < 1e-7);
}
