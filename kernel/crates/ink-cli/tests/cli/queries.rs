//! What ink's queries (mass, distance, interfere, clearance) and claims
//! say of the cases they cannot answer, and why: a body with no density, a
//! published mass, a dead point, a call's time limit, a band the value
//! lies outside. Kernel only; they were ext's (ink-ext tests/cli
//! studies.rs), where a change to their words here went untested.

use crate::common::*;
use serde_json::{Value, json};
use std::collections::BTreeMap;

/// An egg crate's cells, 0.16 mm across in a 60 x 60 x 1 slab: no proof
/// fits in its boxes, nor its estimate in its cells.
const CRATE: &str = "(max (box 60 60 1) (* (sin (* 20 x)) (sin (* 20 y))))";

fn message_has(value: &Value, expected: &[&str]) {
    let message = value.as_str().expect("diagnostic message");
    for text in expected {
        assert!(message.contains(text), "expected {text:?} in {message}");
    }
}

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

/// A call answering afresh (INK_NO_CACHE), its TMPDIR the scratch's.
fn fresh<S: AsRef<std::ffi::OsStr>>(sc: &Scratch, args: impl IntoIterator<Item = S>) -> Ink {
    sc.ink(args).env("INK_NO_CACHE", "1")
}

/// The first claim of an `ink require` answer.
fn first(d: &Value) -> Value {
    d["requirements"][0].clone()
}

#[test]
fn an_interfere_past_the_time_limit_begins_no_moving_pair() {
    // Interfere's moving pairs are stepped one after another: past the
    // limit the rest are not begun (each pair's setup alone ran Limestone's
    // 49,000 ten minutes past its limit).
    let sc = Scratch::new("st-interfere-limit");
    let blocks: String = (0..6).map(|i| format!(" (instance s{i} block (at (move {} 0 0)))", 3 * i)).collect();
    let slid = sc.model(
        "limit-slide.ink",
        &format!("(import \"std\")\n(t-span 0 1s)\n(body block (shape (box 2 2 2)))\n(assembly root (instance m block (at (move (* 15 (/ t 1s)) 0 1.5))){blocks})\n"),
    );
    let d = fresh(&sc, ["interfere", &slid]).env("INK_DEADLINE", deadline(-1.0)).ok();
    let pairs = d["pairs"].as_array().cloned().unwrap_or_default();
    let moving: Vec<&Value> = pairs.iter().filter(|p| p["pair"].as_array().is_some_and(|x| x.iter().any(|t| t == "root/m"))).collect();
    assert!(d["stopped"] == json!(true) && moving.len() >= 5 && moving.iter().all(|p| verdict(p) == "unknown" && has(&p["hint"], "not checked")), "{}", pretty(&d));
}

#[test]
fn an_interfere_with_cli_timeout_stops_cleanly_and_reports_what_finished() {
    let sc = Scratch::new("st-interfere-cli-timeout");
    let blocks: String = (0..6).map(|i| format!(" (instance s{i} block (at (move {} 0 0)))", 3 * i)).collect();
    let slid = sc.model(
        "timeout-slide.ink",
        &format!(
            "(import \"std\")\n(t-span 0 1s)\n(body block (shape (box 2 2 2)))\n\
             (assembly root (instance m block (at (move (* 15 (/ t 1s)) 0 1.5))){blocks})\n"
        ),
    );
    // A tiny --timeout expires before the moving-pair loop runs, stopping cleanly
    // with `stopped: true` and returning all pairs (unstarted ones marked unknown).
    let d = fresh(&sc, ["interfere", &slid, "--timeout", "0.0001"]).ok();
    let pairs = d["pairs"].as_array().cloned().unwrap_or_default();
    let moving: Vec<&Value> = pairs.iter().filter(|p| p["pair"].as_array().is_some_and(|x| x.iter().any(|t| t == "root/m"))).collect();
    assert!(d["stopped"] == json!(true) && moving.len() >= 5 && moving.iter().all(|p| verdict(p) == "unknown" && has(&p["hint"], "not checked")), "{}", pretty(&d));

    // A generous --timeout with a unit suffix completes without stopping.
    let done = fresh(&sc, ["interfere", &slid, "--timeout", "2min"]).ok();
    assert!(done.get("stopped").is_none() && done["overlapping_pairs"].as_u64().unwrap_or(0) >= 1, "{}", pretty(&done));
}

#[test]
fn mass_reads_the_materials_density_else_its_density_attr_and_a_density_named_that_a_body_lacks_is_an_error() {
    // petg's 1270 kg/m^3 is 0.00127 g/mm^3, so a 20 mm cube is 10.16 g; its
    // own density attribute wins. A density named that a body lacks
    // (misspelled) is an error, not a mass of volume times 1.
    let sc = Scratch::new("st-material");
    let mat = sc.model(
        "mat.ink",
        "(import \"std\")\n(body a (shape (box 20 20 20)) (attr material \"petg\"))\n(body b (shape (box 20 20 20)) (attr material \"petg\") (attr density 2g/cm3))\n",
    );
    let m = masses(&sc.ink(["mass", &mat]).ok());
    assert!((m["a"] - 10.16).abs() < 0.01 && (m["b"] - 16.0).abs() < 0.01, "{m:?}");
    let e = sc.ink(["mass", &mat, "--density", "densitty"]).err()["error"].clone();
    assert!(e["message"] == "'a' has no attr 'densitty'" && e["hint"] == "its attrs: material", "mass: {}", pretty(&e));
}

/// A unit cube of 0.001 g/mm^3 a hundred thousand kilometres out keeps
/// its principal moments, 1/6000 g mm^2 each: its moments are taken about
/// a point near it, not about the origin, where its distance squared
/// would round them to nothing. A density attribute below zero is refused
/// before anything is weighed.
#[test]
fn a_distant_body_keeps_its_moments_and_a_negative_density_is_refused() {
    let sc = Scratch::new("mass-far");
    let far = sc.model(
        "far.ink",
        "(import \"std\")\n(body cube (shape (move (box 1 1 1) 100000000 100000000 100000000)) (attr density 0.001g/mm3))\n(assembly root (instance cube cube))\n",
    );
    let d = sc.ink(["mass", &far, "--target", "root", "--keys", "mass,principal_moments"]).ok();
    let p = d["total"]["principal_moments"].as_array().cloned().unwrap_or_default();
    assert!(p.len() == 3 && p.iter().all(|v| (v.as_f64().unwrap() * 6000.0 - 1.0).abs() < 1e-6), "{}", pretty(&d));
    let signed = sc.model("signed.ink", "(import \"std\")\n(body b (shape (box 2 2 2)) (attr mx (- x 2)))\n");
    let e = sc.ink(["mass", &signed, "--density", "mx"]).err()["error"].clone();
    assert!(has(&e["message"], "--density mx on 'b'") && has(&e["message"], "-3.95"), "{}", pretty(&e));
}

/// A steel body, and two placings of a body with no density.
const BARE: &str = r#"(import "std")
(body steel (shape (box 10 10 10)) (attr material "steel"))
(body bare (shape (box 10 10 10)))
(def shape (box 10 10 10))
(assembly root (instance a steel) (instance b bare (at (move 20 0 0))) (instance c bare (at (move 40 0 0))))
"#;

#[test]
fn a_body_with_no_density_has_no_mass_and_the_total_names_it_once_by_body() {
    // Never its volume at 1 g/mm^3 summed with grams; the total names it
    // rather than leave it out (a stepper's stator, before it read the
    // motor's published density, put 66,724 "g" in a total of 66,734). Why
    // is said once, in the total, by body: a robot placing one bare body
    // 936 times said it on every target (422 KB of 922). A shape that is no
    // body has no mass, picked keys and all.
    let sc = Scratch::new("st-bare");
    let bare = sc.model("bare.ink", BARE);
    let d = sc.ink(["mass", &bare]).ok();
    let t = by(&d["targets"], "target");
    let b = &t["root/b"];
    let no_mass = &d["total"]["no_mass"];
    assert!(
        (n(&t["root/a"]["mass"]) - 7.85).abs() < 0.01
            && b["mass"].is_null()
            && b["center_of_mass"].is_null()
            && b["density"].is_null()
            && (n(&b["volume"]["estimate"]) / 1000.0 - 1.0).abs() < 0.02
            && t.values().all(|x| x.get("no_mass").is_none())
            && d["total"]["mass"].is_null()
            && d["total"]["without_mass"] == json!({"bare": ["root/b", "root/c"]})
            && has(no_mass, "(attr material NAME)")
            && has(no_mass, "--target T --density"),
        "{}",
        pretty(&d)
    );
    let d = sc.ink(["mass", &bare, "--target", "shape", "--keys", "mass,inertia"]).ok();
    assert!(d["targets"][0]["mass"].is_null() && d["targets"][0]["inertia"].is_null() && d["total"]["without_mass"] == json!({"shape": ["shape"]}), "{}", pretty(&d));
}

#[test]
fn a_claim_over_a_body_with_no_density_names_it_and_its_file() {
    // A claim says it in the same words; an imported body is
    // named with its file.
    let sc = Scratch::new("st-bare-claim");
    sc.model("bare.ink", BARE);
    let claim = sc.model("bare-claim.ink", "(import \"bare.ink\" as k)\n(assembly root (instance b k.bare))\n(require light (integral root (attr density) < 10))\n");
    let e = sc.ink(["check", &claim]).err()["error"].clone();
    let m = e["message"].as_str().unwrap_or("");
    assert!(
        m.starts_with("'root/b' (body 'k.bare' in ")
            && m.ends_with("bare.ink) has no attr 'density'")
            && has(&e["hint"], "Give its body (attr material NAME) or (attr density 7.85g/cm3) in the file that defines it")
            && has(&e["hint"], "an assembly placing it a published (attr mass M)"),
        "{}",
        pretty(&e)
    );
}

#[test]
fn an_assemblys_published_mass_is_shared_by_its_bodies_measured_at_the_params_defaults() {
    // The bodies it places weigh it together, each placed body counted, a
    // steel shaft at its own density, the rest sharing what remains,
    // measured at the params' defaults: a wider case weighs more.
    let sc = Scratch::new("st-drive-mass");
    let am = sc.model(
        "drive-mass.ink",
        r#"(import "std")
(param w 10)
(body case (shape (box w 10 10)))
(body out (shape (move (box 2 2 2) 0 0 6)))
(body shaft (shape (move (box 1 1 4) 3 0 0)) (attr material "steel"))
(assembly (drive (angle 0)) (attr mass 20g) (instance case case) (instance o1 out) (instance o2 out (at (move 5 0 0)))
  (instance s shaft (at) (turn z angle)))
"#,
    );
    let d = sc.ink(["mass", &am, "--target", "drive"]).ok();
    let t = by(&d["targets"], "target");
    assert!(
        (n(&d["total"]["mass"]) - 20.0).abs() < 1e-6
            && (n(&t["drive/s"]["mass"]) - 4.0 * 0.00785).abs() < 1e-4
            && t["drive/o1"]["mass"] == t["drive/o2"]["mass"]
            && t["drive/case"]["density"].as_str().is_some_and(|s| s.starts_with("published mass of 'drive', less what its bodies with a density of their own weigh")),
        "{}",
        pretty(&d)
    );
    let d = sc.ink(["mass", &am, "--target", "drive", "--param", "w=20"]).ok();
    assert!((n(&d["total"]["mass"]) - (20.0 + n(&d["targets"][0]["mass"]) / 2.0)).abs() < 1e-3, "{}", pretty(&d));
}

#[test]
fn an_assemblys_mass_that_gives_no_body_a_density_is_a_warning() {
    // Every body here has its own material, so the 5000 g is not used (it
    // read 174 g, silently): check warns, its assembly says why it is not
    // met, and the bodies weigh their materials.
    let sc = Scratch::new("st-inert");
    let inert = sc.model(
        "inert-mass.ink",
        r#"(import "std")
(body plate (shape (box 50 50 5)) (attr material "aluminum-6061"))
(body post (shape (move (cylinder-z 5 40) 0 0 20)) (attr material "steel"))
(assembly root (attr mass 5000g) (instance p plate) (instance q post))
"#,
    );
    let c = sc.ink(["check", &inert, "--load"]).ok();
    let (w, a) = (&c["warnings"], &c["assemblies"][0]["mass"]);
    let d = sc.ink(["mass", &inert]).ok();
    let want = 50.0 * 50.0 * 5.0 * 0.0027 + std::f64::consts::PI * 25.0 * 40.0 * 0.00785;
    assert!(
        w.as_array().is_some_and(|w| w.len() == 1)
            && has(&w[0], "the (attr mass ...) of assembly 'root' gives no body a density")
            && n(&a["published"]) == 5000.0
            && a["met"].is_null()
            && a["why"].as_str().is_some_and(|s| s.starts_with("gives no body a density"))
            && (n(&d["total"]["mass"]) - want).abs() < 1.0,
        "{}\n{}\n{}",
        pretty(w),
        pretty(a),
        pretty(&d["total"])
    );
}

#[test]
fn check_warns_of_a_mass_too_dense_for_its_bodies_boxes_as_mass_refuses_it() {
    // 1.2 kg in 36 cm³ (trial 9's first payload): 0.033 g/mm³, past any
    // material, over its box alone: check says what mass and effort say.
    let sc = Scratch::new("st-dense");
    let dense = sc.model("dense.ink", "(import \"std\")\n(body payload (shape (box 60 30 20)) (attr mass 1.2kg))\n(assembly root (instance p payload))\n");
    let c = sc.ink(["check", &dense]).ok();
    let m = sc.ink(["mass", &dense]).err();
    assert!(
        any_has(&c["warnings"], "payload.mass-density: 1200 g published mass implies density at least 0.03")
            && any_has(&c["warnings"], "its bodies' boxes hold at most 36")
            && has(&m["error"]["message"], "do not establish the expected 0.00000005..0.023"),
        "{}\n{}",
        pretty(&c["warnings"]),
        pretty(&m)
    );
}

#[test]
fn a_full_check_warns_of_a_published_mass_that_cannot_be_weighed_and_load_measures_none() {
    // Lighter than its own-density bodies: its assembly says why it is not met.
    let sc = Scratch::new("st-heavy");
    let heavy = sc.model(
        "heavy-mass.ink",
        "(import \"std\")\n(body s (shape (box 20 20 20)) (attr material \"steel\"))\n(body c (shape (move (box 5 5 5) 0 0 20)))\n(assembly root (attr mass 20g) (instance s s) (instance c c))\n",
    );
    let (c, cl) = (sc.ink(["check", &heavy]).ok(), sc.ink(["check", &heavy, "--load"]).ok());
    let (a, al) = (&c["assemblies"][0]["mass"], &cl["assemblies"][0]["mass"]);
    assert!(
        any_has(&c["warnings"], "more than its published 20 g")
            && cl["warnings"] == json!([])
            && n(&a["published"]) == 20.0
            && a["met"] == json!(false)
            && has(&a["why"], "more than its published 20 g")
            // Check refuses a mass without weighing it, so --load says it is checked there, not measured.
            && same(al, &json!({"published": 20, "met": null, "why": "checked by ink check, not with --load"})),
        "{}\n{}\n{}\n{}",
        pretty(&c["warnings"]),
        pretty(&cl["warnings"]),
        pretty(a),
        pretty(al)
    );
}

#[test]
fn a_published_mass_nothing_weighs_does_not_refuse_a_toleranced_claim() {
    // A claim that does not weigh it holds over the band (it was refused:
    // "mass 'b.mass-density' depends on ranged 'w'"); one that weighs it
    // over the band still is.
    let sc = Scratch::new("st-tol-mass");
    let base = "(import \"std\")\n(param w 20 (tol -0.5 0.5))\n(body b (shape (box w 10 10)) (attr mass 10g))\n(require fits (subset (shape-of b) (box 30 30 30)))\n";
    let a = sc.ink(["require", &sc.model("tol-mass.ink", base)]).ok();
    let r = sc.ink(["require", &sc.model("tol-weighed.ink", &format!("{base}(require light (integral b (attr density) < 200))\n"))]).err();
    assert!(a["counts"] == json!({"yes": 1}) && has(&r["error"]["message"], "depends on ranged 'w'"), "{}\n{}", pretty(&a), pretty(&r));
}

/// Two proofs that rested on what a field does not show. 1e16 + (x - x + 0.5) less 1e16 is 0.5 in
/// real numbers; x - x dropped to a float 0 let the sums fold in floats to 0, and "<= 0" held. Now
/// the dropped term is an exact 0 and what folds with it a band about its rounding, so the claim is
/// not a yes (a float can not resolve 0.5 beside 1e16: not settled); so of a term of t over a
/// motion. The same sum of constants as written is 0.5 in real numbers too: a value claim on constants
/// is decided in them, not by the floats it folds to (which say 0), so "<= 0" is a no. And
/// a slope promise whose k reads x, warped by a frame that jumps (a select), bounded the field along
/// a way that leaves the frame's image: the range claim held where the field is +2; now it fails.
#[test]
fn a_proof_rests_on_no_more_than_the_field_shows() {
    let sc = Scratch::new("field-shows");
    let src = "(import \"std\")\n(require zero (value (- (+ 1e16 (+ (- x x) 0.5)) 1e16) <= 0))\n(require times (value (- (+ 1e16 (+ (* 0 x) 0.5)) 1e16) <= 0))\n\
               (require written (value (- (+ 1e16 0.5) 1e16) <= 0))\n";
    let d = sc.ink(["require", &sc.model("fold.ink", src)]).any();
    let status = |name: &str| d["requirements"].as_array().unwrap().iter().find(|r| r["name"] == name).unwrap()["status"].clone();
    assert!(status("zero") != "holds" && status("times") != "holds" && status("written") == "fails", "{}", pretty(&d));
    let time = "(import \"std\")\n(t-span 0 1s)\n(require of-t (value (- (+ 1e16 (+ (* 0 (/ t 1s)) 0.5)) 1e16) <= 0))\n";
    let d = sc.ink(["require", &sc.model("time.ink", time)]).any();
    assert_ne!(d["requirements"][0]["status"], "holds", "{}", pretty(&d));
    // Evaluated at a point, such a band is no float: a shape truly 0.2 outside everywhere, its float
    // 0.3 inside, is given with its band and not said inside.
    let band = sc.model("band.ink", "(import \"std\")\n(body b (shape (- (- (+ 1e16 (+ (* 0 x) 0.5)) 1e16) 0.3)))\n");
    let e = sc.ink(["eval", &band, "--points", "[[0,0,0]]"]).any();
    let row = &e["targets"][0]["results"]["rows"][0];
    let at = |key: &str| e["targets"][0]["results"]["columns"].as_array().unwrap().iter().position(|c| c == key).map(|i| row[i].clone());
    assert!(at("inside") == Some(Value::Null) && at("band").is_some_and(|b| b[0].as_f64() < Some(0.2) && b[1].as_f64() > Some(0.2)), "{}", pretty(&e));
    let warp = "(import \"std\")\n(def s (select y -1 1))\n(def k (+ 0.001 (* 3 (abs (* (- x s) (+ x s))))))\n(def f (- (* (square x) x) (* 3 x)))\n\
                (def folded (warp (lipschitz k f) s y z))\n(require range-claim (range folded -3 -1 (box 2 2 2)))\n";
    let d = sc.ink(["require", &sc.model("warp.ink", warp)]).any();
    assert_eq!(d["requirements"][0]["status"], "fails", "{}", pretty(&d));
}

/// A claim on a material's strength: Prusament PLA's yield is the maker's typical coupon result,
/// so a claim reading it is an estimate naming it, never a proven yes; 6061-T6's is a
/// specification minimum, and the same claim on it is proven; the number written out is too.
#[test]
fn a_typical_strength_makes_a_claim_on_it_an_estimate() {
    let sc = Scratch::new("typical-strength");
    let src = "(import \"std\")\n(def stress 30)\n(require pla (value stress < (get (material \"pla\") \"yield\")))\n\
               (require al (value stress < (* 0.8 (get (material \"aluminum-6061-t6-extrusion\") \"yield\"))))\n(require written (value stress < 51))\n";
    let d = sc.ink(["require", &sc.model("strength.ink", src)]).ok();
    let row = |name: &str| d["requirements"].as_array().unwrap().iter().find(|r| r["name"] == name).unwrap().clone();
    let pla = row("pla");
    assert!(verdict(&pla) == "holds estimate" && pla["typical_strengths"] == json!(["pla-prusament-fff-horizontal yield (datasheet typical)"]), "{}", pretty(&pla));
    for name in ["al", "written"] {
        assert!(verdict(&row(name)) == "holds" && row(name).get("typical_strengths").is_none(), "{}", pretty(&row(name)));
    }
}

/// A published mass written with its tolerance, 4.19 g within 0.5 mg, nearly all of it a known 1 cm
/// ball: what is left for a 1 mm cube is 1.2 mg at 4.19 g, and 0.7 to 1.7 mg over the tolerance. A
/// claim on the cube is asked over the tolerance, its density found again at each mass: under 3 mg
/// holds, under 1.3 mg fails (it held at the mass's value alone). Exact, the mass leaves 1.2 mg. A
/// tolerance whose low end leaves the cube nothing (the ball outweighs it) says so.
#[test]
fn a_published_masss_tolerance_reaches_what_shares_it() {
    let sc = Scratch::new("mass-tolerance");
    let src = |m: &str, bound: &str| {
        format!(
            "(import \"std\")\n(param m {m})\n(body known (shape (sphere 10)) (attr density 0.001g/mm3))\n(body plug (shape (box 1 1 1)))\n\
             (assembly root (attr mass m) (instance b known) (instance p plug (at (move 12 0 0))))\n(require light (integral root/p (attr density) < {bound}))\n"
        )
    };
    let ask = |name: &str, m: &str, bound: &str| {
        let d = sc.ink(["require", &sc.model(name, &src(m, bound))]).any();
        d["requirements"][0].clone()
    };
    let tol = "4.19g (tol -0.5mg 0.5mg)";
    assert_eq!(verdict(&ask("wide.ink", tol, "3mg")), "holds estimate");
    let tight = ask("tight.ink", tol, "1.3mg");
    assert!(verdict(&tight).starts_with("fails") && tight["value"].as_f64().is_some_and(|v| v > 1.3e-3), "{}", pretty(&tight));
    assert_eq!(verdict(&ask("exact.ink", "4.19g", "1.3mg")), "holds estimate");
    let none = ask("none.ink", "4.1888g (tol -2mg 2mg)", "1.3mg");
    assert!(none["status"] == "unknown" && pretty(&none).contains("leaves these bodies nothing"), "{}", pretty(&none));
}

/// A 100 g assembly of a rigid base and a spring that shrinks over the motion: each keeps the mass
/// it has at t = 0 (71.8 and 28.2 g), where one density found again at each moment moved 8.5 g from
/// the spring to the base by t = 1.5.
#[test]
fn bodies_sharing_a_published_mass_keep_theirs_as_one_deforms() {
    let sc = Scratch::new("kept-share");
    let m = sc.model(
        "spring.ink",
        "(import \"std\")\n(t-span 0 2s)\n(body base (shape (box 20 20 10)))\n(body spring (shape (move (cylinder-z 5 (- 20 (* 5 (/ t 1s)))) 30 0 0)))\n\
         (assembly root (attr mass 100g) (instance b base) (instance s spring))\n",
    );
    let at = |t: &str| -> Vec<f64> {
        let d = ink(["mass", &m, "--t", t]);
        d["targets"].as_array().unwrap().iter().map(|x| x["mass"].as_f64().unwrap()).collect()
    };
    let (start, later) = (at("0"), at("1.5"));
    assert!((start[0] + start[1] - 100.0).abs() < 0.5 && (start[0] - 71.8).abs() < 0.3, "{start:?}");
    for (a, b) in start.iter().zip(&later) {
        assert!((a - b).abs() < 0.05 * a / 100.0 + 1e-6, "{start:?} then {later:?}");
    }
}

#[test]
fn a_published_mass_whose_volume_does_not_settle_loads_and_says_so() {
    // An egg crate: its density is the finest measurement, and what
    // reads it says unsettled, and why (it was a load error). A sheet
    // (1000:1) no longer is: its walls are proven.
    let sc = Scratch::new("st-sheet");
    let sheet = sc.model("sheet.ink", &format!("(import \"std\")\n(body sheet (shape {CRATE}) (attr mass 2.7g))\n"));
    let s = sc.ink(["mass", &sheet]).ok();
    assert!(
        (n(&s["total"]["mass"]) / 2.7 - 1.0).abs() < 0.02
            && has(&s["targets"][0]["unsettled"], "halving")
            && has(&s["unsettled"]["sheet.mass-density"], "halving")
            && s["solved_inputs"] == json!(["sheet.mass-density"]),
        "{}",
        pretty(&s)
    );
    let claims = sc.model("claims.ink", r#"(import "sheet.ink") (require weight (integral sheet (attr density) < 1000g)) (require width (value 8 < 7))"#);
    let report = sc.ink(["require", &claims]).any();
    assert_eq!(report["requirements"][0]["status"], "unknown");
    message_has(&report["requirements"][0]["unavailable"], &["sheet.mass-density", "measurement bounds", "0.00000005..0.023"]);
    assert_eq!(report["requirements"][1]["status"], "fails");
}

#[test]
fn check_lists_an_assemblys_published_mass_met_and_unsettled_where_its_density_is() {
    // Met by the density it derives (check listed no assembly's mass: a
    // part's showed its bodies' density, not the mass).
    let sc = Scratch::new("st-sheets");
    let nested = sc.model(
        "sheets.ink",
        &format!("(import \"std\")\n(body sheet (shape {CRATE}))\n(body cube (shape (move (box 10 10 10) 0 0 20)))\n(assembly inner (attr mass 2.7g) (instance s sheet))\n(assembly root (attr mass 10g) (instance i inner) (instance c cube))\n"),
    );
    let c = sc.ink(["check", &nested]).ok();
    let a: BTreeMap<String, Value> = by(&c["assemblies"], "name").into_iter().map(|(k, v)| (k, v["mass"].clone())).collect();
    let (root, inner) = (&a["root"], &a["inner"]);
    assert!(
        n(&root["published"]) == 10.0
            && root["met"] == json!(true)
            && root["density"] == "root.mass-density"
            // Check measures only the densities another mass reads (inner's): root's is measured where read.
            && root.get("unsettled").is_none()
            && inner["met"] == json!(true)
            && has(&inner["unsettled"]["inner.mass-density"], "halving"),
        "{}",
        pretty(&json!(a))
    );
}

#[test]
fn a_published_mass_is_of_the_defaults_and_sums_at_the_placed_centroids() {
    // A size param scales it; its density is measured at the defaults.
    let sc = Scratch::new("st-placed-mass");
    let mat = sc.model(
        "placed.ink",
        "(import \"std\")\n(param w 20)\n(body b (shape (box w 10 10)) (attr mass 19.9g))\n(assembly root (instance a b) (instance b b (at (rot-z 90) (move 40 0 0))))\n",
    );
    for w in [20.0, 40.0] {
        let m = sc.ink(["mass", &mat, "--param", &format!("w={w}")]).ok();
        let c = &m["total"]["center_of_mass"];
        let off = ((n(&c[0]) - 20.0).powi(2) + n(&c[1]).powi(2) + n(&c[2]).powi(2)).sqrt();
        assert!(
            (n(&m["total"]["mass"]) / (39.8 * w / 20.0) - 1.0).abs() < 0.01
                && off < 0.01
                && m["targets"]
                    .as_array()
                    .is_some_and(|t| t.iter().all(|t| t["density"].as_str().is_some_and(|d| d.starts_with("published mass of 'b' over its volume") && d.contains("centroid"))))
                && m.get("solved_inputs").is_some(),
            "w = {w}: {}",
            pretty(&m)
        );
    }
}

/// What reads a field for more than its sign reads only what the field shows. A term negative
/// everywhere (1e-10 under 0 at z = 3) adds no surface, so a 10 mm cube is thick, not 2 mm. A strip
/// swept round a 3 mm bend on edge (its 10 mm square to the bend's plane) does not fold; laid flat it
/// does.
#[test]
fn a_field_is_read_for_its_sign_and_only_what_it_shows() {
    let sc = Scratch::new("field-sign");
    let thick = sc.model("thick.ink", "(import \"std\")\n(def part (max (box 10 10 10) (- -1e-10 (abs (- z 3)))))\n(require thick (min-thickness part > 6) (res 0.5))\n");
    let d = sc.ink(["require", &thick]).ok();
    assert_eq!(verdict(&d["requirements"][0]), "holds estimate", "{}", pretty(&d));
    let warns = |name: &str, src: &str| -> Vec<String> {
        let d = sc.ink(["check", &sc.model(name, src)]).ok();
        d["warnings"].as_array().map(|w| w.iter().map(|x| x.as_str().unwrap().to_string()).collect()).unwrap_or_default()
    };
    let arc = "(def arc (map (lambda (k) (let ((a (* pi (/ k 39) (/ 2 3)))) (vec (* 3 (cos a)) (* 3 (sin a)) 0))) (range 40)))";
    assert!(warns("edge.ink", &format!("(import \"std\")\n{arc}\n(body s (shape (sweep (rect-2d 10 1) arc (vec 0 0 1))))\n")).is_empty());
    let flat = warns("flat.ink", &format!("(import \"std\")\n{arc}\n(body s (shape (sweep (rect-2d 1 10) arc (vec 0 0 1))))\n"));
    assert!(flat.len() == 1 && flat[0].contains("turns too sharply"), "{flat:?}");
}

#[test]
fn a_distance_to_an_empty_assembly_and_a_thickness_no_sample_measured_say_so() {
    // An empty assembly is the empty shape to distance, as an empty body is;
    // a wall no sample falls on is not measured: never a passing thickness.
    let sc = Scratch::new("st-empty");
    let mat = sc.model(
        "empty.ink",
        "(assembly root)\n(body b (shape (- (+ (square x) (square y) (square z)) 4)))\n(body thin (shape (max (- (abs x) 10) (- (abs y) 10) (- (abs z) 0.01))))\n(require thick (min-thickness (shape-of thin) > 5))\n",
    );
    let e = sc.ink(["distance", &mat, "--a", "root", "--b", "b"]).err()["error"]["message"].clone();
    let a = first(&sc.ink(["require", &mat]).ok());
    assert!(e == "a target is empty" && has(&a["unavailable"], "no surface sample") && a["status"] == "unknown" && a["value"].is_null(), "{e}\n{}", pretty(&a));
}

#[test]
fn a_clearance_no_cell_size_decides_ends_at_its_first_stuck_box_with_the_accuracy_it_holds_to() {
    // A box of slabs turned 45 degrees about z, extruded, reads 21.9
    // beside its corner where the other's corner is 31.02 away (boxes of
    // slabs along the axes are proven apart by the gaps between their
    // sides). Claimed clear by 31, best first, the box nearest a violation
    // is refined first; stuck at the finest cells, it ends the way: under a
    // second, where passes to the limit took 8 s per claim and 30 s per
    // distance. A 2 mm cube, turned too, 24.001 mm off a's corner, claimed
    // clear by 24: the way covering a (many times the undecided volume)
    // waits while the cube's leads. Each gives the accuracy it holds to
    // (a's field proves 21.9 of the 31 claimed), and the finest cells as
    // their narrowest side (a cell is never halved along the span).
    let sc = Scratch::new("st-edges");
    let edges = sc.model(
        "edges.ink",
        r#"(fn turned (cx h) (max (- (abs (* (+ (- x cx) y) .7071067811865476)) h) (- (abs (* (- (- x cx) y) .7071067811865476)) h) (- (abs z) h)))
(body a (shape (turned 0 50)))
(body b (shape (turned 172.44135623730952 50)))
(body s (shape (turned 96.12589168102785 1)))
(assembly root (instance a a) (instance b b))
(require c (clearance (shape-of a) (shape-of b) 31))
(require lag (clearance (shape-of a) (shape-of s) 24))
"#,
    );
    let d = fresh(&sc, ["distance", &edges, "--a", "root/a", "--b", "root/b"]).ok()["distance"].clone();
    let a = fresh(&sc, ["require", &edges]).ok()["requirements"].clone();
    let (lo, at, hi) = (n(&d["lower_proven"]), n(&d["undecided_at"]), n(&d["upper_witness"]));
    let rows = a.as_array().cloned().unwrap_or_default();
    let acc = n(&a[0]["accuracy_attained"]);
    // What a's field reads at b's corner.
    let field = 31.02 / 2f64.sqrt();
    let finest = |x: &Value| x["hint"].as_str().and_then(|h| h.split("cells, about ").nth(1)?.split(' ').next()?.parse::<f64>().ok());
    assert!(
        field - 0.1 < lo
            && lo < at
            && at < 31.02
            && (31.02..31.03).contains(&hi)
            && rows.len() == 2
            && rows.iter().all(|x| x["status"] == "unknown" && n(&x["undecided_cells"]) < 100000.0 && has(&x["hint"], "It holds to") && finest(x).is_some_and(|f| f < 5.0))
            && (31.0 - field..31.01 - field).contains(&acc),
        "{}\n{}",
        pretty(&d),
        pretty(&a)
    );
}

#[test]
fn a_cavities_claim_past_20000_cells_across_is_an_error() {
    let sc = Scratch::new("st-voids-fine");
    let fine = sc.model("voids-fine.ink", "(import \"std\")\n(require v (cavities (box 20 20 20) < 1) (res 0.0005))\n");
    let e = sc.ink(["require", &fine]).err();
    assert!(e.to_string().contains("the limit is 20000"), "{}", pretty(&e));
}

#[test]
fn an_overlap_with_a_probe_box_thinner_than_twice_the_accuracy_is_unknown_its_hint_saying_so() {
    // A probe box thinner than twice the accuracy holds no point deeper than
    // it: unknown, saying so and what decides it (claims are answered once).
    let sc = Scratch::new("st-probe-box");
    let pb = sc.model(
        "probe-box.ink",
        "(import \"std\")\n(body f (shape (box 10 10 10)))\n(require reach (overlaps (shape-of f) (box 0.02 0.02 0.05)))\n(require reach-fine (overlaps (shape-of f) (box 0.02 0.02 0.05)) (accuracy 0.005))\n",
    );
    let a = by(&sc.ink(["require", &pb]).ok()["requirements"], "name");
    assert!(
        verdict(&a["reach"]) == "unknown"
            && has(&a["reach"]["hint"], "under twice the accuracy")
            && has(&a["reach"]["hint"], "(probe x y z)")
            && verdict(&a["reach-fine"]) == "holds",
        "{}",
        pretty(&json!(a))
    );
}

/// A failure names the bodies at its witness, by the paths of the instances placing them, for
/// each shape it compares: the two that overlap, the one member of a whole assembly there, the
/// nearest of each side of a clearance; none where a probe is in no body.
#[test]
fn a_failure_names_the_bodies_at_its_witness() {
    let s = Scratch::new("at-witness");
    let f = s.model(
        "w.ink",
        "(import \"std\")\n(body a (shape (box 10 10 10)))\n(body b (shape (sphere 4)))\n(body c (shape (sphere 2)))\n\
         (assembly root (instance left a (at (move -3 0 0))) (instance right b (at (move 4 0 0))) (instance far c (at (move 0 30 0))))\n\
         (require apart (disjoint (shape-of root/left) (shape-of root/right)))\n(require whole (disjoint (shape-of root) (move (sphere 1) 0 30 0)))\n\
         (require gap (clearance (shape-of root/far) (shape-of root/left) 30))\n(require in (subset (probe 0 30 0) (shape-of root/left)))\n",
    );
    let d = s.ink(["require", &f]).any();
    let at = |name: &str| d["requirements"].as_array().unwrap().iter().find(|r| r["name"] == name).unwrap()["at_witness"].clone();
    assert!(
        at("apart") == json!([["root/left"], ["root/right"]])
            && at("whole") == json!([["root/far"], []])
            && at("gap") == json!([["root/far"], ["root/left"]])
            && at("in") == json!([[]]),
        "{}",
        pretty(&d)
    );
}

/// A shape is a regular solid: two blocks written flush are one, a point on their join is in it
/// (main read it as outside, a 0), one on the outer face is not, and a point claim, a subset
/// across the join and `ink eval` all say so; a block cut through by a gap with a width is two
/// solids, a point in the gap outside.
#[test]
fn a_shape_is_a_regular_solid() {
    let s = Scratch::new("regular-solid");
    let f = s.model(
        "flush.ink",
        "(import \"std\")\n(body two (shape (min (box 10 10 10) (move (box 10 10 10) 10 0 0))))\n\
         (body cut (shape (subtract (box 10 10 10) (box 0.2 20 20))))\n\
         (require join-in (subset (probe 5 0 0) (shape-of two)))\n(require face-out (disjoint (probe 15 0 0) (shape-of two)))\n\
         (require face-in (subset (probe 15 0 0) (shape-of two)))\n(require across (subset (move (box 2 2 2) 5 0 0) (shape-of two)))\n\
         (require gap-out (disjoint (probe 0 0 0) (shape-of cut)))\n",
    );
    let d = s.ink(["require", &f]).any();
    let status = |name: &str| d["requirements"].as_array().unwrap().iter().find(|r| r["name"] == name).unwrap()["status"].clone();
    let solids = |t: &str| s.ink(["components", &f, "--target", t]).ok()["targets"][0]["solids"].clone();
    let at = s.ink(["eval", &f, "--target", "two", "--points", "[[5,0,0],[15,0,0]]"]).ok();
    assert!(
        status("join-in") == "holds"
            && status("face-out") == "holds"
            && status("face-in") == "fails"
            && status("across") == "holds"
            && status("gap-out") == "holds"
            && solids("two") == 1
            && solids("cut") == 2
            && at["targets"][0]["results"]["rows"][0][2] == true
            && at["targets"][0]["results"]["rows"][1][2] == false,
        "{}\n{}",
        pretty(&d),
        pretty(&at)
    );
}

/// A cut that is 0 where it is empty, (subtract a (min 0 x)), leaves the part's field 0 over the
/// half where x > 0, which by what a shape is reads as outside: check says so, by its point, and a
/// claim reads it as the definition does (a point there, a box there, are not in the part).
#[test]
fn a_field_of_0_over_a_volume_is_said_not_silently_read() {
    let s = Scratch::new("zero-volume");
    let f = s.model(
        "zero.ink",
        "(import \"std\")\n(body part (shape (min (move (box 4 4 4) -8 0 0) (subtract (box 10 10 10) (min 0 x)))))\n\
         (require point-in (subset (probe 2 0 0) (shape-of part)))\n(require box-in (subset (move (box 1 1 1) 2 0 0) (shape-of part)))\n\
         (require block-in (subset (move (box 1 1 1) -8 0 0) (shape-of part)))\n",
    );
    let check = s.ink(["check", &f]).ok();
    let warned = check["warnings"].as_array().unwrap().iter().any(|w| has(w, "its field is 0 throughout a volume"));
    let d = s.ink(["require", &f]).any();
    let row = |name: &str| d["requirements"].as_array().unwrap().iter().find(|r| r["name"] == name).unwrap().clone();
    assert!(warned && row("point-in")["status"] == "fails" && row("box-in")["status"] == "fails" && row("block-in")["status"] == "holds", "{}\n{}", pretty(&check), pretty(&d));
}

/// Two plates of one thickness filleted where they meet: just above their caps each reads
/// distance to its cap's edge, and the plain fillet adds a ridge past them (material above z =
/// 2.6), which ink check warns of, naming the form; given the caps' plane the fillet fills the
/// corner below it and nothing above, and check says nothing.
#[test]
fn a_fillet_of_shapes_with_caps_in_one_plane_is_warned_of_and_capped() {
    let s = Scratch::new("fillet-caps");
    let f = s.model(
        "caps.ink",
        "(import \"std\")\n(body ridge (shape (fillet-union (box 80 40 5) (box 200 20 5) 2 90)))\n\
         (body capped (shape (fillet-union (box 80 40 5) (box 200 20 5) 2 90 (- (abs z) 2.5))))\n\
         (require ridge-under (empty (max (shape-of ridge) (- 2.6 z))))\n(require capped-under (empty (max (shape-of capped) (- 2.6 z))))\n\
         (require capped-fillets (subset (probe -40.5 10.5 0) (shape-of capped)))\n",
    );
    let check = s.ink(["check", &f]).ok();
    let warned = |body: &str| check["warnings"].as_array().unwrap().iter().any(|w| has(w, &format!("body '{body}': a fillet-union joins two shapes")));
    let d = s.ink(["require", &f]).any();
    let status = |name: &str| d["requirements"].as_array().unwrap().iter().find(|r| r["name"] == name).unwrap()["status"].clone();
    assert!(
        warned("ridge") && !warned("capped") && status("ridge-under") == "fails" && status("capped-under") == "holds" && status("capped-fillets") == "holds",
        "{}\n{}",
        pretty(&check),
        pretty(&d)
    );
}

#[test]
fn a_value_outside_its_own_band_is_warned_of_and_measured_saying_so() {
    // 8 f9 written from its basic size: the band 7.951..7.987 excludes 8.
    // check warns, naming the param and its band; distance and mass, which
    // work at the value, say they measured 8 (outside_band), unless --param
    // gives the size; a value in its band says nothing.
    let s = Scratch::new("outside-band");
    let src = |pin: &str| {
        format!(
            "(import \"std\")\n(param bore-d 8.4 (tol 0 0.036))\n{pin}\n(body bracket (shape (subtract (box 6 20 20) (cylinder-x (/ bore-d 2) 30))))\n(body pin (shape (cylinder-x (/ pin-d 2) 10)))\n(assembly root (instance b bracket) (instance p pin))\n"
        )
    };
    let basic = s.model("basic.ink", &src("(param pin-d 8 (tol -0.049 -0.013))"));
    let inside = s.model("inside.ink", &src("(param pin-d 7.969 (tol -0.018 0.018))"));
    let warned = |f: &str| s.ink(["check", f]).ok()["warnings"].as_array().cloned().unwrap_or_default().iter().any(|w| has(w, "param pin-d is 8 with its band 7.951..7.987"));
    let measured = |f: &str, extra: &[&str]| {
        let mut args = vec!["distance", f, "--a", "root/p", "--b", "root/b"];
        args.extend(extra);
        s.ink(args).ok()["outside_band"].clone()
    };
    let at_basic = measured(&basic, &[]);
    let mass = s.ink(["mass", &basic]).ok()["outside_band"].clone();
    assert!(
        warned(&basic)
            && !warned(&inside)
            && at_basic["pin-d"] == json!({"value": 8.0, "band": [7.951, 7.987]})
            && mass["pin-d"]["value"] == 8.0
            && measured(&basic, &["--param", "pin-d=7.969"]).is_null()
            && measured(&inside, &[]).is_null(),
        "{}\n{}",
        pretty(&at_basic),
        pretty(&mass)
    );
}

/// A volume, integral or mean claim is decided by its proven bounds, not its estimate, which has none: a
/// ball of radius 10 (4188.790 mm^3) was answered yes to "volume > 4189" and "mass > 4.189 g".
#[test]
fn an_integral_claim_is_not_decided_by_its_estimate_alone() {
    let sc = Scratch::new("integral-estimate");
    let m = sc.model(
        "ball.ink",
        "(import \"std\")\n(body ball (shape (sphere 10)) (attr density 1g/cm3))\n\
         (require heavier (integral ball (attr density) > 4.189))\n(require bigger (volume (sphere 10) > 4189))\n",
    );
    let d = sc.ink(["require", &m]).any();
    for r in d["requirements"].as_array().unwrap() {
        assert_ne!(r["status"], "holds", "4188.790 is not over its limit: {}", pretty(r));
    }
}

/// `ink eval` says inside by the field's outward-rounded bounds, as a probe claim does: x < 0.3 written
/// (x + 1e9) - 1e9 - 0.3 has (0.30000001, 0, 0) 1e-8 outside, where its float says inside; the bounds
/// reach across 0, so it is not said inside.
#[test]
fn eval_does_not_say_inside_where_the_bounds_reach_across_zero() {
    let sc = Scratch::new("eval-absorb");
    let m = sc.model("h9.ink", "(import \"std\")\n(body h9 (shape (max (- (- (+ x 1e9) 1e9) 0.3) (box 10 10 10))))\n");
    let e = sc.ink(["eval", &m, "--points", "[[0.30000001,0,0]]", "--keys", "inside"]).any();
    let inside = &e["targets"][0]["results"]["rows"][0][0];
    assert_ne!(*inside, Value::Bool(true), "{}", pretty(&e));
}

/// A value claim on constants is one float, checked as a point and called exact: 1 ft / 12 in > 1 is a
/// yes, sin 30 deg >= 0.5 a no, cos 90 deg > 0 a yes, sqrt 2 squared > 2 a yes; each was wrong. Now
/// decided in real numbers: the foot a proven no, the rest within rounding of their limits, not settled.
#[test]
fn a_value_claim_on_constants_carries_their_rounding() {
    let sc = Scratch::new("val-float");
    let m = sc.model(
        "val.ink",
        "(import \"std\")\n(require ft (value (/ 1ft 12in) > 1))\n(require sin30 (value (sin 30deg) >= 0.5))\n\
         (require cos90 (value (cos 90deg) > 0))\n(require root2 (value (* (sqrt 2) (sqrt 2)) > 2))\n",
    );
    let d = sc.ink(["require", &m]).any();
    let status = |name: &str| d["requirements"].as_array().unwrap().iter().find(|r| r["name"] == name).unwrap()["status"].clone();
    assert!(status("ft") != "holds" && status("sin30") != "fails" && status("cos90") != "holds" && status("root2") != "holds", "{}", pretty(&d));
    // Decimals and units written as decimals are exact fractions, so what they decide is proven: a foot is
    // twelve inches, three tenths are three tenths.
    let exact = sc.model(
        "exact.ink",
        "(import \"std\")\n(require ft-is-12in (value (- 1ft 12in) >= 0))\n(require ft-not-over (value (- 1ft 12in) > 0))\n\
         (require tenths (value (+ 0.1 0.2) <= 0.3))\n(require tenth (value (* 3 0.1) > 0.3))\n",
    );
    let d = sc.ink(["require", &exact]).any();
    let got: Vec<(String, String)> = d["requirements"].as_array().unwrap().iter().map(|r| (r["name"].as_str().unwrap().into(), r["status"].as_str().unwrap().into())).collect();
    let want = [("ft-is-12in", "holds"), ("ft-not-over", "fails"), ("tenths", "holds"), ("tenth", "fails")].map(|(a, b)| (a.to_string(), b.to_string()));
    assert_eq!(got, want, "{}", pretty(&d));
}

/// A param written in digits is, in a value claim, the fraction they write, as a literal is: with big a
/// param of 1e16, (big + 1) - big >= 1 is a proven yes (floats said 0, a proven no), and a cable carrier's
/// 334 mm in links of a record's 16.7 (ceil(334 / 16.7) of them) is exactly 334, a proven yes at >= and a
/// proven no at >; times a written 1.00000000000000001 it is more than itself (its float is 1). A param
/// computed through something inexact (1e16 times sqrt 2) is its float, its
/// arithmetic rounded outward as a range's: not settled.
#[test]
fn a_value_claim_reads_a_written_param_as_its_digits() {
    let sc = Scratch::new("param-exact");
    sc.model("chain.json", r#"{"pitch": 16.7}"#);
    let m = sc.model(
        "cancels.ink",
        "(import \"std\")\n(param big 1e16)\n(param far (* 1e16 (sqrt 2)))\n(param length 334)\n(def pitch (get (read-json \"chain.json\") \"pitch\"))\n\
         (def links (neg (floor (neg (/ length pitch)))))\n\
         (require cancels (value (- (+ big 1) big) >= 1))\n(require stock (value (* pitch links) >= length))\n\
         (require over (value (* pitch links) > length))\n(require inexact (value (- (+ far 1) far) >= 1))\n\
         (require not-one (value (* length 1.00000000000000001) <= length))\n",
    );
    let d = sc.ink(["require", &m]).env("INK_NO_CACHE", "1").any();
    let got: Vec<(String, String)> =
        d["requirements"].as_array().unwrap().iter().map(|r| (r["name"].as_str().unwrap().to_string(), r["status"].as_str().unwrap().to_string())).collect();
    let want = [("cancels", "holds"), ("stock", "holds"), ("over", "fails"), ("inexact", "unknown"), ("not-one", "fails")].map(|(a, b)| (a.to_string(), b.to_string()));
    assert_eq!(got, want, "{}", pretty(&d));
}

/// `ink mass` states no error it cannot stand behind: a ball whose field is a step (-1 in, +1 out) is
/// estimated on cells (4189.0146 mm^3, where it is 4188.790), and given with the bounds the proof reached,
/// which hold it, and no spread (two halvings that happened to agree made one of 1.6e-5 g).
#[test]
fn a_sampled_mass_states_its_proven_bounds_and_no_spread() {
    let sc = Scratch::new("mass-spread");
    let m = sc.model("step.ink", "(import \"std\")\n(body stepf (shape (select (- (len3 x y z) 10) -1 1)) (attr density 1g/cm3))\n");
    let t = &sc.ink(["mass", &m]).any()["targets"][0];
    let b = &t["mass_bounds"];
    assert!(t.get("spread").is_none() && b[0].as_f64() <= Some(4.188790205) && Some(4.188790205) <= b[1].as_f64(), "{}", pretty(t));
}
