//! Tolerance bands and per-claim res, meshes imported as sampled fields
//! (and what answers on samples or uncertified roots say), and `ink render`
//! (its PNG read here without the kernel).

use crate::common::*;
use serde_json::{Value, json};

/// Each claim of an `ink assert` answer by name.
/// Whether `v` is the band `[lo, hi]` as checked: written exactly
/// (command::bound), within float steps of the decimals that made it.
fn band(v: &Value, lo: f64, hi: f64) -> bool {
    let near = |x: &Value, y: f64| x.as_f64().is_some_and(|x| (x - y).abs() < 1e-12);
    near(&v[0], lo) && near(&v[1], hi) && v.as_array().map(Vec::len) == Some(2)
}

fn by_name(d: &Value) -> std::collections::BTreeMap<String, Value> {
    d["asserts"].as_array().unwrap().iter().map(|a| (a["name"].as_str().unwrap().to_string(), a.clone())).collect()
}

/// Whether any warning of an answer contains every one of these.
fn warned(d: &Value, all: &[&str]) -> bool {
    d["warnings"].as_array().is_some_and(|w| w.iter().any(|w| w.as_str().is_some_and(|w| all.iter().all(|s| w.contains(s)))))
}

// ---------------------------------------------------------------- tolerances

/// A shaft 7.95..8.05 in a bore 8.00..8.10 (declared on the params): 0.025
/// radial clearance at the middle sizes, 0.025 interference at shaft 8.05 in
/// bore 8.00, which an accuracy of 0.005 sees. Bands apply by default;
/// --nominal drops them.
const FIT: &str = "(import \"std\")\n(param shaft-d 8 (tol -0.05 0.05))\n(param bore-d 8.05 (tol -0.05 0.05))\n\
(body shaft (shape (cylinder-z (/ shaft-d 2) 10)))\n\
(body sleeve (shape (subtract (cylinder-z 6 14) (cylinder-z (/ bore-d 2) 16))))\n\
(assert apart (disjoint (shape-of shaft) (shape-of sleeve)) (accuracy 0.005))\n\
(assert gap (value (- bore-d shaft-d) > 0))\n";

#[test]
fn check_lists_each_tolerance_band() {
    let s = Scratch::new("tol-check");
    let d = s.ink(["check", &s.model("tolerances.ink", FIT)]).ok();
    let p = d["params"].as_array().unwrap();
    assert!(
        p.len() == 2 && p[0]["name"] == "shaft-d" && band(&p[0]["band"], 7.95, 8.05) && p[1]["name"] == "bore-d" && band(&p[1]["band"], 8.0, 8.1) && d["warnings"] == json!([]),
        "{}",
        pretty(&d)
    );
}

#[test]
fn a_fit_fails_inside_its_tolerances_with_an_in_band_witness() {
    let s = Scratch::new("tol-fit");
    let d = s.ink(["assert", &s.model("tolerances.ink", FIT)]).ok();
    let a = by_name(&d);
    let (apart, gap) = (&a["apart"], &a["gap"]);
    let spans = |s: &Value| s.as_object().map(|o| o.len()) == Some(2) && band(&s["shaft-d"], 7.95, 8.05) && band(&s["bore-d"], 8.0, 8.1);
    let p = |a: &Value, k: &str| a["params"][k].as_f64().unwrap_or(f64::NAN);
    assert!(
        verdict(apart) == "fails"
            && (7.95..=8.05).contains(&p(apart, "shaft-d"))
            && (8.0..=8.1).contains(&p(apart, "bore-d"))
            && p(apart, "shaft-d") > p(apart, "bore-d")
            && spans(&d["param_spans"])
            && spans(&apart["param_spans"]),
        "{}",
        pretty(&d)
    );
    // A value over tolerances fails with its params.
    assert!(
        verdict(gap) == "fails" && gap["value"].as_f64().unwrap() <= 0.0 && gap["params"].as_object().unwrap().keys().collect::<Vec<_>>() == ["shaft-d", "bore-d"],
        "{}",
        pretty(gap)
    );
}

#[test]
fn nominal_checks_at_the_values_only() {
    let s = Scratch::new("tol-nominal");
    let d = s.ink(["assert", &s.model("tolerances.ink", FIT), "--nominal"]).ok();
    let v: Vec<String> = d["asserts"].as_array().unwrap().iter().map(verdict).collect();
    assert!(v == ["holds", "holds"] && d["param_spans"] == json!({}), "{}", pretty(&d));
}

#[test]
fn param_moves_a_band_ranges_a_param_or_pins_one() {
    let s = Scratch::new("tol-param");
    let tol = s.model("tolerances.ink", FIT);
    let d = s.ink(["assert", &tol, "--param", "bore-d=8.15"]).ok();
    assert!(verdict(&d["asserts"][0]) == "holds" && band(&d["param_spans"]["bore-d"], 8.1, 8.2), "a band follows its param's value: {}", pretty(&d));
    let d = s.ink(["assert", &tol, "--param", "bore-d=8.15:8.25", "--param", "shaft-d=8:8"]).ok();
    assert!(
        verdict(&d["asserts"][0]) == "holds" && d["param_spans"].as_object().map(|o| o.len()) == Some(1) && band(&d["param_spans"]["bore-d"], 8.15, 8.25),
        "k=lo:hi ranges, k=v:v pins: {}",
        pretty(&d)
    );
}

#[test]
fn a_param_given_twice_or_ranged_for_a_one_value_command_is_refused() {
    let s = Scratch::new("tol-refuse");
    let tol = s.model("tolerances.ink", FIT);
    let e = s.ink(["assert", &tol, "--param", "shaft-d=7.9:8", "--param", "shaft-d=8.1:8.1"]).err();
    assert!(e["error"]["message"].as_str().unwrap().contains("shaft-d is given twice"), "{}", pretty(&e));
    let e = s.ink(["bounds", &tol, "--param", "bore-d=8:8.1"]).err();
    assert!(e["error"]["message"].as_str().unwrap().contains("assert and interfere"), "{}", pretty(&e));
}

#[test]
fn interfere_covers_a_toleranced_placement() {
    // interfere never places a part at one value of a toleranced param (its
    // box covers the band), so the pair is checked and found.
    let s = Scratch::new("tol-placed");
    let pl = s.model(
        "placed.ink",
        "(import \"std\")\n(param gap 0.5 (tol -0.8 0.6))\n(body a (shape (box 10 10 10)))\n\
        (assembly root (instance a a) (instance b a (at (move (+ 10 gap) 0 0))))\n",
    );
    let d = s.ink(["interfere", &pl]).ok();
    assert!(d["overlapping_pairs"] == 1 && d["pairs"][0]["params"]["gap"].as_f64().is_some_and(|g| g < 0.0), "{}", pretty(&d));
    let d = s.ink(["interfere", &pl, "--nominal"]).ok();
    assert!(d["overlapping_pairs"] == 0, "at its value the placement is clear: {}", pretty(&d));
}

#[test]
fn check_warns_of_a_toleranced_param_used_at_load_time() {
    // It cannot be covered: assert refuses it (geometry.rs
    // a_param_whose_value_is_used_at_load_time_is_not_ranged), check warns.
    let s = Scratch::new("tol-baked");
    let bk = s.model(
        "baked-tol.ink",
        "(import \"std\")\n(param n 3 (tol 0 1))\n(body plate (shape (if (> n 2) (box 10 10 1) (box 5 5 1))))\n\
        (assert some (overlaps (shape-of plate) (sphere 1)))\n",
    );
    let d = s.ink(["check", &bk]).ok();
    assert!(warned(&d, &["'n'"]), "{}", pretty(&d));
}

// ---------------------------------------------------------------- per-claim res

#[test]
fn an_asserts_own_res_is_used_where_finer_than_res_and_says_so() {
    // A hexagonal prism 0.0025 inside a bore.
    let s = Scratch::new("claim-res");
    let hx = s.model(
        "claim-res.ink",
        "(import \"std\")\n(def hex (extrude-z (polygon (regular-polygon-points 6 4.9975)) -2 2))\n\
        (assert tight (subset hex (cylinder-z 5 6)) (res 0.01))\n(assert plain (subset hex (cylinder-z 5 6)))\n",
    );
    let d = s.ink(["assert", &hx, "--res", "0.02"]).ok();
    let got: Vec<(String, Value, Value)> = d["asserts"].as_array().unwrap().iter().map(|a| (verdict(a), a["res_given"].clone(), a["res_given_by"].clone())).collect();
    assert!(got == [("holds".into(), json!(0.01), json!("assert")), ("holds".into(), json!(0.02), json!("--res"))], "{}", pretty(&d));
}

// ---------------------------------------------------------------- import

/// A 20 x 30 x 10 box centred on the origin as a binary STL, written here
/// (12 triangles: cheaper to sample than a meshed part), the first `n` of
/// its triangles; its path.
fn box_stl(s: &Scratch, name: &str, n: usize) -> String {
    let corner = |i: usize| [0, 1, 2].map(|k| [[-10f32, 10.0], [-15.0, 15.0], [-5.0, 5.0]][k][(i >> k) & 1]);
    let tris: Vec<[usize; 3]> =
        [[0, 2, 6, 4], [1, 5, 7, 3], [0, 4, 5, 1], [2, 3, 7, 6], [0, 1, 3, 2], [4, 6, 7, 5]].iter().flat_map(|&[a, b, c, d]| [[a, b, c], [a, c, d]]).take(n).collect();
    let mut out = vec![0u8; 80];
    out.extend((tris.len() as u32).to_le_bytes());
    for t in tris {
        out.extend([0f32; 3].into_iter().chain(t.iter().flat_map(|&i| corner(i))).flat_map(f32::to_le_bytes));
        out.extend([0, 0]);
    }
    let path = s.path(name);
    std::fs::write(&path, out).unwrap();
    path
}

/// A model importing box.stl, sampled at 0.3, then `rest`; its path.
fn imported_box(s: &Scratch, rest: &str) -> String {
    box_stl(s, "box.stl", 12);
    s.model("imp.ink", &format!("(import \"std\")\n(body q (shape (import-stl \"box.stl\" (res 0.3))))\n{rest}"))
}

#[test]
fn an_stl_ink_writes_reads_back_as_its_shape_and_queries_name_its_samples() {
    // A box with a hole, written by `ink mesh`, read back sampled at 0.5
    // (its volume is 0.24% off there; at 0.3, 0.11%, for twice the time).
    let s = Scratch::new("import-trip");
    let part = s.model("part.ink", "(import \"std\")\n(body p (shape (subtract (box 20 30 10) (cylinder-z 5 20))))\n");
    s.ink(["mesh", &part, "--res", "0.25", "--out", &s.path("part.stl")]).ok();
    let imp = s.model("imp.ink", "(import \"std\")\n(body q (shape (import-stl \"part.stl\" (res 0.5))))\n");
    let d = s.ink(["mass", &imp, "--res", "1"]).ok();
    let v = d["targets"][0]["volume"]["estimate"].as_f64().unwrap();
    let exact = 6000.0 - std::f64::consts::PI * 25.0 * 10.0;
    assert!(((v - exact) / exact).abs() < 5e-3, "STL round trip volume {v} vs {exact}");
    assert!(d["samples"][0]["source"] == "part.stl", "a query on an import names its samples: {}", pretty(&d));
    // The imported field is a signed distance.
    let d = s.ink(["eval", &imp, "--points", "[[7,0,0],[0,20,0]]"]).ok();
    // Its per-point results are a table ({columns, rows}).
    let t = &d["targets"][0]["results"];
    let col = t["columns"].as_array().unwrap().iter().position(|c| c == "value").unwrap_or_else(|| panic!("no value column: {}", pretty(t)));
    let vals: Vec<f64> = t["rows"].as_array().unwrap().iter().map(|r| r[col].as_f64().unwrap()).collect();
    assert!((vals[0] + 2.0).abs() < 0.05 && (vals[1] - 5.0).abs() < 0.05, "{vals:?}");
}

#[test]
fn an_import_reports_its_sampling_error_measured() {
    // Measured, with no bound (it is sampled in floats): under the cell's
    // half diagonal.
    let s = Scratch::new("import-error");
    let d = s.ink(["check", &imported_box(&s, "")]).ok();
    let sm = &d["samples"][0];
    let max = sm["measured_error"]["max"].as_f64().unwrap_or(f64::NAN);
    assert!(
        sm["kind"] == "mesh"
            && sm["error_bound"].is_null()
            && 0.0 < max
            && max <= 0.3 * 3f64.sqrt() / 2.0
            && d["warnings"] == json!([])
            && d["bodies"][0]["samples"] == json!(["box.stl"]),
        "{}",
        pretty(&d)
    );
}

#[test]
fn claims_and_pairs_on_an_import_are_measured_on_its_samples() {
    let s = Scratch::new("import-claims");
    // The box's top face is at z = 5: the probes clear it by 1 and by 0.1.
    let ic = imported_box(
        &s,
        "(body b (shape (box 2 2 2)))\n\
        (assert apart (disjoint (shape-of q) (move (shape-of b) 0 10 7)))\n\
        (assert apart-0.1 (disjoint (shape-of q) (move (shape-of b) 0 10 6.1)))\n\
        (assert plain (disjoint (shape-of b) (move (shape-of b) 0 10 7)))\n\
        (assembly root (instance q q) (instance b b (at (move 0 10 6.1))))\n",
    );
    let d = s.ink(["assert", &ic]).ok();
    let a = by_name(&d);
    assert!(
        verdict(&a["apart"]) == "holds estimate"
            && verdict(&a["apart-0.1"]) == "holds estimate"
            && a["apart"]["samples"][0]["source"] == "box.stl"
            && a["apart"]["samples"][0]["res"] == 0.3
            && verdict(&a["plain"]) == "holds"
            && a["plain"].get("samples").is_none(),
        "{}",
        pretty(&d)
    );
    let d = s.ink(["interfere", &ic]).ok();
    let p0 = &d["pairs"][0];
    assert!(verdict(p0) == "holds estimate" && p0["samples"][0]["source"] == "box.stl", "{}", pretty(&d));
}

#[test]
fn an_open_mesh_bounds_no_solid_and_check_says_so() {
    // A mesh that is not closed has no inside.
    let s = Scratch::new("import-open");
    box_stl(&s, "open.stl", 11);
    let io = s.model(
        "imp-open.ink",
        "(import \"std\")\n(body q (shape (import-stl \"open.stl\" (res 0.3))))\n\
        (assert apart (disjoint (shape-of q) (move (box 2 2 2) 0 10 7)))\n",
    );
    let d = s.ink(["check", &io]).ok();
    assert!(warned(&d, &["not closed"]), "{}", pretty(&d));
}

#[test]
fn an_imports_box_is_measured_and_cannot_shape_the_model_at_load_time() {
    // The box of sampled data bounds the samples: a floor sized from it is
    // measured, naming the box, and a condition on it is refused, naming
    // what was sampled.
    let s = Scratch::new("import-bounds");
    box_stl(&s, "box.stl", 12);
    let head = "(import \"std\")\n(def b (bounds (import-stl \"box.stl\" (res 0.3))))\n(def x1 (nth (nth b 1) 0))\n";
    let ib = s.model("imp-bounds.ink", &format!("{head}(body floor (shape (move (box (* 2 x1) (* 2 x1) 1) 0 0 -10)))\n(assert wide (volume (shape-of floor) > 300))\n"));
    let d = s.ink(["assert", &ib]).ok();
    let w = &by_name(&d)["wide"];
    let kinds: Vec<&Value> = w["samples"].as_array().unwrap().iter().map(|x| &x["kind"]).collect();
    assert!(verdict(w) == "holds estimate" && kinds == [&json!("bounds"), &json!("mesh")], "{}", pretty(&d));
    let ic = s.model("imp-if.ink", &format!("{head}(def n (if (> x1 0) 1 2))\n"));
    let e = s.ink(["check", &ic]).err();
    let m = e["error"]["message"].as_str().unwrap();
    assert!(m.contains("reads sampled data (bounds of (import-stl \"box.stl\" (res 0.3)) at ") && m.contains("imp-if.ink:2)"), "{}", pretty(&e));
    // Found empty, it is refused, default or not: the default would carry a
    // measured decision as exact (a ball under the grid's spacing is lost).
    let empty = s.model("empty.ink", "(import \"std\")\n(def lost (grid-bake (sphere 0.1) (res 1)))\n(def d (bounds lost 7))\n(assert leak (value (if (> d 5) 1 0) > 0))\n");
    let e = s.ink(["assert", &empty]).err();
    assert!(e["error"]["message"].as_str().unwrap().starts_with("bounds found no interior in sampled data (bake of (sphere 0.1)"), "{}", pretty(&e));
    // A form's setting (an optimizer's limit, its budget) is decided at
    // load, so it is refused too; a solve's guess, where its search starts,
    // is not (the solve's answer reads it).
    for setting in ["(vary w 1 x1)", "(vary w 1 9) (budget x1)"] {
        let m = s.model("setting.ink", &format!("{head}(param w 1)\n(optimize least {setting} (minimize (value w)))\n"));
        let e = s.ink(["check", &m]).err();
        assert!(e["error"]["message"].as_str().unwrap().contains("reads sampled data (bounds of"), "{setting}: {}", pretty(&e));
    }
    let guess = s.model("guess.ink", &format!("{head}(solve near (unknowns (a x1)) (== (* a a) 4))\n"));
    s.ink(["check", &guess]).ok();
}

#[test]
fn a_numerical_claims_integral_and_mean_are_estimates_their_parts_with_them() {
    let s = Scratch::new("part-accuracy");
    // A claim resting on a numerical root is an estimate (accuracy null),
    // and so are its integral's and mean's own statuses (b rests on a
    // double root found numerically).
    let nr = s.model(
        "numerical-parts.ink",
        "(import \"std\")\n(param p 0 (tol 0 0.00000001))\n(solve first (unknowns (a 1)) (== (square (- a 2)) 0))\n\
        (solve second (unknowns (b 2)) (== b (+ a p)))\n(def cube (box 1 1 1))\n\
        (assert false-integral (integral cube (- b 2) < 0))\n(assert false-mean (mean cube b < 2))\n",
    );
    let d = s.ink(["assert", &nr, "--res", "0.1"]).ok();
    let a = by_name(&d);
    let got = [verdict(&a["false-integral"]), verdict(&a["false-integral"]["integral"]), verdict(&a["false-mean"]), verdict(&a["false-mean"]["mean"])];
    assert!(got.iter().all(|v| v == "holds estimate"), "{got:?}: {}", pretty(&d));
}

#[test]
fn a_mass_claims_number_is_its_value_and_an_imported_params_default_is_refused_saying_what_to_do() {
    // The Results table says a claim's number is `value`: a mass claim
    // gives it there, not only under integral.estimate. And an importer
    // redefining an imported param is told to give it at the import.
    let s = Scratch::new("premise2-docs");
    s.model("a.ink", "(import \"std\")\n(param base-t 5)\n(body plate (shape (box 10 10 base-t)) (attr density 2.7g/cm3))\n");
    let m = s.model("mass.ink", "(import \"a.ink\")\n(assert light (integral plate (attr density) < 2))\n(assert low (mean plate z (weight (attr density)) < 1))\n");
    let a = by_name(&s.ink(["assert", &m]).ok());
    let n = |v: &Value| v.as_f64().unwrap_or(f64::NAN);
    let (mass, low) = (n(&a["light"]["value"]), n(&a["low"]["value"]));
    assert!((mass - 1.35).abs() < 0.03 && mass == n(&a["light"]["integral"]["estimate"]) && low.abs() < 0.01 && low == n(&a["low"]["mean"]["estimate"]), "{mass} {low}");
    let b = s.model("b.ink", "(import \"a.ink\")\n(param base-t 3.2)\n");
    let e = s.ink(["check", &b]).err();
    assert!(e["error"]["hint"].as_str().is_some_and(|h| h.ends_with("a.ink\" (base-t V))")), "{}", pretty(&e));
}

// ---------------------------------------------------------------- import values

/// A design with a claim, a body, a study and an optimization, all reading
/// its param base-t (5 mm, +-0.1), and a range of the question.
const DESIGN: &str = "(import \"std\")\n(param base-t 5mm (tol -0.1 0.1))\n(def reach (ranged 0 1))\n(def w 3)\n\
(body plate (shape (box 10 10 base-t)))\n(assert thin (value base-t < 4))\n\
(solve twice (unknowns (a 1)) (== a (* 2 base-t)))\n\
(optimize lighter (vary base-t 1 6) (minimize (value base-t)) (budget 2))\n";

#[test]
fn an_imports_value_is_its_params_default_for_the_whole_load() {
    let num = |v: &Value| v.as_f64().unwrap_or(f64::NAN);
    let s = Scratch::new("import-values");
    let ink = |args: &[&str]| s.ink(args).cwd(&s);
    s.model("a.ink", DESIGN);
    s.model("b.ink", "(import \"a.ink\" (base-t 3.2))\n");
    let b = "b.ink";
    // The claim, the body, the study and optimize's start all read 3.2;
    // check says what was written and where it was set.
    let a = by_name(&ink(&["assert", b]).ok());
    assert!(a["thin"]["status"] == "holds" && band(&a["thin"]["param_spans"]["base-t"], 3.1, 3.3), "{}", pretty(&a["thin"]));
    let top = ink(&["eval", b, "--target", "plate", "--points", "[[0,0,0]]"]).ok();
    assert!((num(&records(&top["targets"][0]["results"])[0]["value"]) + 1.6).abs() < 1e-9, "the plate is 3.2 thick: {}", pretty(&top));
    let c = ink(&["check", b, "--imports"]).ok();
    let param = |name: &str| c["params"].as_array().unwrap().iter().find(|p| p["name"] == name).cloned().unwrap_or_else(|| panic!("{name}: {}", pretty(&c)));
    let base = param("base-t");
    assert!(base["default"] == 3.2 && base["written"] == 5.0 && base["set_at"] == "b.ink:1", "{}", pretty(&base));
    assert!((num(&param("a")["value"]) - 6.4).abs() < 1e-9, "the study solves at 3.2: {}", pretty(&param("a")));
    let o = ink(&["optimize", b]).ok();
    assert!(o["start"]["design"]["base-t"] == 3.2, "{}", pretty(&o));
    // --param beats it, for one command.
    let a = by_name(&ink(&["assert", b, "--param", "base-t=4.5"]).ok());
    assert!(a["thin"]["status"] == "fails" && a["thin"]["value"] == 4.5, "{}", pretty(&a["thin"]));
    // `as ns` gives a second copy, with its own value; (tol lo hi) replaces
    // the band, given with a unit of the param's quantity.
    s.model("two.ink", "(import \"a.ink\" (base-t 0.125in (tol 0 0.2)))\n(import \"a.ink\" as thin (base-t 2))\n");
    let a = by_name(&ink(&["assert", "two.ink"]).ok());
    assert!(band(&a["thin"]["param_spans"]["base-t"], 3.175, 3.375), "{}", pretty(&a["thin"]));
    assert!(band(&a["thin.thin"]["param_spans"]["thin.base-t"], 1.9, 2.1), "{}", pretty(&a["thin.thin"]));
    // A value known at load is given as a constant is: the importer's param
    // passes on, --param on the importer reaching the import.
    s.model(
        "pass.ink",
        "(param t0 3.2)
(import \"a.ink\" as q (base-t t0))
",
    );
    let a = by_name(&ink(&["assert", "pass.ink", "--param", "t0=4.5"]).ok());
    assert!(a["q.thin"]["status"] == "fails" && a["q.thin"]["value"] == 4.5, "{}", pretty(&a["q.thin"]));
}

/// What `ink check` of `src` (as m.ink, beside a.ink: a 5 mm param p) says
/// is wrong, and the line.
fn import_refused(src: &str) -> (String, u64) {
    let s = Scratch::new("import-units");
    s.model("a.ink", "(import \"std\")\n(param p 5mm)\n");
    s.model("m.ink", src);
    let e = s.ink(["check", "m.ink"]).cwd(&s).err();
    (e["error"]["message"].as_str().unwrap_or("").to_string(), e["error"]["at"]["line"].as_u64().unwrap_or(0))
}

#[test]
fn an_imports_band_is_said_in_its_params_quantity() {
    // A band in MPa for a length, the value given bare: refused at the
    // import, as the param written so is.
    let (e, line) = import_refused("\n(import \"a.ink\" (p 3 (tol -1MPa 1MPa)))\n");
    assert!(e.contains("a.ink's p takes a length") && e.contains("MPa is a pressure") && line == 2, "{e} (line {line})");
}

#[test]
fn a_repeated_imports_value_is_said_in_its_params_quantity() {
    // 3MPa is not the 3mm the first import gave, though the numbers match.
    let (e, line) = import_refused("(import \"a.ink\" (p 3mm))\n(import \"a.ink\" (p 3MPa))\n");
    assert!(e.contains("a.ink's p takes a length") && e.contains("MPa is a pressure") && line == 2, "{e} (line {line})");
}

#[test]
fn an_imports_value_is_refused_where_the_load_cannot_give_it() {
    let s = Scratch::new("import-values-refused");
    let ink = |args: &[&str]| s.ink(args).cwd(&s);
    s.model("a.ink", DESIGN);
    s.model("b.ink", "(import \"a.ink\" (base-t 3.2))\n");
    s.model("c.ink", "(import \"std\")\n(param base-t 5)\n");
    s.model("via.ink", "(import \"c.ink\")\n");
    let refused = |src: &str| {
        s.model("m.ink", src);
        let e = ink(&["check", "m.ink"]).err();
        let e = &e["error"];
        (format!("{} | {}", e["message"].as_str().unwrap_or(""), e["hint"].as_str().unwrap_or("")), e["at"]["line"].as_u64().unwrap_or(0))
    };
    for (src, says, line) in [
        // Two importers asking two values of one load's file.
        ("(import \"b.ink\")\n(import \"a.ink\" (base-t 2.5))\n", "a.ink's base-t is already 3.2mm, set by b.ink:1; this import asks 2.5", 2),
        ("(import \"b.ink\")\n(import \"a.ink\" (base-t 2.5))\n", "(import \"a.ink\" as NS (base-t 2.5))", 2),
        // A value after a plain import ran the file at its default.
        ("(import \"a.ink\")\n(import \"a.ink\" (base-t 2.5))\n", "a.ink was read at m.ink:1 with base-t 5mm | give the value at a.ink's first import (m.ink:1)", 2),
        // A unit of another quantity, cited at the import.
        ("\n(import \"a.ink\" (base-t 3psi))\n", "a.ink's base-t takes a length", 2),
        // What is not one of the file's params.
        ("(import \"a.ink\" (reach 0.5))\n", "'reach' is a (ranged ...), a range of the question: an import sets params", 1),
        ("(import \"a.ink\" (w 2))\n", "'w' is a def, not a param", 1),
        ("(import \"a.ink\" (a 2))\n", "'a' is made by a package form", 1),
        ("(import \"a.ink\" (bse-t 2))\n", "a.ink has no param 'bse-t' | did you mean 'base-t'?", 1),
        ("(import \"via.ink\" (base-t 2))\n", "base-t is c.ink's: import c.ink with it, before via.ink", 1),
    ] {
        let (e, at) = refused(src);
        assert!(e.contains(says) && at == line, "{src}: {e} (line {at})");
    }
    // The same value twice, and a plain import after a valued one, read
    // the one run.
    s.model("same.ink", "(import \"b.ink\")\n(import \"a.ink\" (base-t 3.2))\n(import \"a.ink\")\n");
    assert!(by_name(&ink(&["assert", "same.ink"]).ok())["thin"]["status"] == "holds");
}

#[test]
fn queries_on_an_uncertified_root_name_it_and_are_measured() {
    // u stops at its guess 1; the root is 2.
    let s = Scratch::new("early-root");
    let er = s.model("early-root.ink", "(import \"std\")\n(solve first (unknowns (u 1)) (== (* 1e-13 (- u 2)) 0))\n(body q (shape (sphere u)))\n");
    let er = er.as_str();
    for args in [vec!["bounds", er], vec!["mass", er, "--res", "0.1"]] {
        let d = s.ink(&args).ok();
        assert!(d["solved_inputs"] == json!(["u"]) && d.get("measured").is_some(), "{}: {}", args[0], pretty(&d));
    }
}

#[test]
fn a_sampled_res_under_the_grids_floor_is_coarsened_with_a_warning() {
    // A grid's spacing is at least 1e-4 of its largest extent: (res .001)
    // over a 1000 mm sliver is sampled at 0.1, said beyond samples[0].res.
    let s = Scratch::new("sample-floor");
    let facet = |v: [[f64; 3]; 3]| format!("facet normal 0 0 0\nouter loop\n{}endloop\nendfacet\n", v.map(|p| format!("vertex {} {} {}\n", p[0], p[1], p[2])).concat());
    let (a, b, c, e) = ([0.0, 0.0, 0.0], [1000.0, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]);
    s.model("skinny.stl", &format!("solid t\n{}{}{}{}endsolid t\n", facet([a, c, b]), facet([a, b, e]), facet([a, e, c]), facet([b, c, e])));
    let m = s.model("floor.ink", "(body b (shape (import-stl \"skinny.stl\" (res 0.001))))\n");
    let d = s.ink(["check", &m]).ok();
    let warned = d["warnings"].as_array().is_some_and(|w| w.iter().any(|w| w.as_str().unwrap().contains("(res 0.001) is under 1e-4 of the grid's largest extent")));
    assert!(warned && d["samples"][0]["res"].as_f64().is_some_and(|r| r > 0.09), "{}", pretty(&d));
    let ok = s.model("fine.ink", "(import \"std\")\n(body b (shape (grid-bake (box 2 2 2) (res 0.25))))\n");
    let d = s.ink(["check", &ok]).ok();
    assert!(d["warnings"].as_array().is_none_or(|w| w.iter().all(|w| !w.as_str().unwrap().contains("1e-4"))), "{}", pretty(&d));
}

#[test]
fn check_says_a_root_of_sampled_equations_is_certified_for_the_samples_only() {
    let s = Scratch::new("sample-root");
    let sr = s.model(
        "sample-root.ink",
        "(import \"std\")\n(def sampled (grid-bake (sphere 1) (res 0.25)))\n(param p 0 (tol 0 0.001))\n\
        (solve first (unknowns (u 1)) (== (warp sampled u 0.07 0.13) p))\n(solve second (unknowns (v 1)) (== v (+ u p)))\n",
    );
    let d = s.ink(["check", &sr]).ok();
    let solves = d["solves"].as_array().unwrap();
    assert!(solves.len() == 2 && solves.iter().all(|x| x["samples"][0]["res"] == 0.25), "{}", pretty(&d["solves"]));
}

#[test]
fn queries_name_the_samples_and_roots_their_attributes_read() {
    // A density solved on samples; bounds reads no attribute.
    let s = Scratch::new("sampled-density");
    let dq = s.model(
        "sampled-density.ink",
        "(import \"std\")\n(def sampled (grid-bake (sphere 1) (res 0.25)))\n\
        (solve s (unknowns (u 1)) (== (warp sampled u 0.07 0.13) 0))\n(body root (shape (box 1 1 1)) (attr density u))\n",
    );
    let dq = dq.as_str();
    for args in [vec!["mass", dq, "--res", "0.1"], vec!["eval", dq, "--points", "[[0,0,0]]"]] {
        let d = s.ink(&args).ok();
        assert!(d["solved_inputs"] == json!(["u"]) && d["samples"][0]["res"] == 0.25, "{}: {}", args[0], pretty(&d));
    }
    let d = s.ink(["bounds", dq]).ok();
    assert!(d.get("samples").is_none() && d.get("solved_inputs").is_none(), "{}", pretty(&d));
}

// ---------------------------------------------------------------- render

/// A PNG ink wrote: its width, height and rows (each a filter byte, then
/// RGB pixels), read without the kernel: chunks, then zlib's inflate.
pub(crate) fn png_rows(path: &str) -> (usize, usize, Vec<Vec<u8>>) {
    let data = std::fs::read(path).unwrap();
    assert_eq!(&data[..8], b"\x89PNG\r\n\x1a\n", "{path}: not a PNG");
    let be = |o: usize| u32::from_be_bytes(data[o..o + 4].try_into().unwrap()) as usize;
    let (w, h) = (be(16), be(20));
    assert_eq!(&data[24..26], &[8, 2], "{path}: not 8-bit RGB");
    let (mut idat, mut pos) = (vec![], 8);
    while pos < data.len() {
        let size = be(pos);
        if &data[pos + 4..pos + 8] == b"IDAT" {
            idat.extend_from_slice(&data[pos + 8..pos + 8 + size]);
        }
        pos += size + 12;
    }
    let raw = inflate(&idat[2..]);
    // zlib's Adler-32 of the data, after the deflate stream: a wrong decode
    // can keep the length.
    let (mut a, mut b) = (1u32, 0u32);
    for &x in &raw {
        a = (a + x as u32) % 65521;
        b = (b + a) % 65521;
    }
    assert_eq!((b << 16 | a).to_be_bytes(), idat[idat.len() - 4..], "{path}: inflated data fails its Adler-32");
    let stride = 1 + 3 * w;
    assert_eq!(raw.len(), stride * h, "{path}: {w} x {h} needs {} bytes, inflated {}", stride * h, raw.len());
    (w, h, raw.chunks(stride).map(<[u8]>::to_vec).collect())
}

/// DEFLATE (RFC 1951): stored, fixed and dynamic Huffman blocks.
fn inflate(d: &[u8]) -> Vec<u8> {
    struct Bits<'a>(&'a [u8], usize);
    impl Bits<'_> {
        fn bit(&mut self) -> usize {
            let b = (self.0[self.1 >> 3] >> (self.1 & 7)) & 1;
            self.1 += 1;
            b as usize
        }
        fn bits(&mut self, n: usize) -> usize {
            (0..n).map(|i| self.bit() << i).sum()
        }
        /// One symbol of canonical code (count per length, symbols by code).
        fn decode(&mut self, (count, sym): &(Vec<usize>, Vec<usize>)) -> usize {
            let (mut code, mut first, mut index) = (0, 0, 0);
            for &n in &count[1..] {
                code |= self.bit();
                if code < first + n {
                    return sym[index + code - first];
                }
                index += n;
                first = (first + n) << 1;
                code <<= 1;
            }
            panic!("bad Huffman code");
        }
    }
    fn huff(lengths: &[usize]) -> (Vec<usize>, Vec<usize>) {
        let mut count = vec![0; 16];
        lengths.iter().for_each(|&l| count[l] += 1);
        count[0] = 0;
        let mut syms: Vec<usize> = (0..lengths.len()).filter(|&s| lengths[s] > 0).collect();
        syms.sort_by_key(|&s| lengths[s]);
        (count, syms)
    }
    const LBASE: [usize; 29] = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
    const LEXT: [usize; 29] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
    const DBASE: [usize; 30] = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
    const DEXT: [usize; 30] = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
    let (mut b, mut out) = (Bits(d, 0), Vec::<u8>::new());
    loop {
        let last = b.bit();
        match b.bits(2) {
            0 => {
                let at = b.1.div_ceil(8);
                let len = u16::from_le_bytes([d[at], d[at + 1]]) as usize;
                out.extend_from_slice(&d[at + 4..at + 4 + len]);
                b.1 = (at + 4 + len) * 8;
            }
            kind @ (1 | 2) => {
                let (lit, dist) = if kind == 1 {
                    let l: Vec<usize> = (0..288)
                        .map(|s| match s {
                            0..=143 => 8,
                            144..=255 => 9,
                            256..=279 => 7,
                            _ => 8,
                        })
                        .collect();
                    (huff(&l), huff(&[5; 30]))
                } else {
                    let (hlit, hdist, hclen) = (b.bits(5) + 257, b.bits(5) + 1, b.bits(4) + 4);
                    let mut cl = [0; 19];
                    for &k in &[16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15][..hclen] {
                        cl[k] = b.bits(3);
                    }
                    let cl = huff(&cl);
                    let mut l = vec![];
                    while l.len() < hlit + hdist {
                        match b.decode(&cl) {
                            s @ 0..=15 => l.push(s),
                            16 => {
                                let prev = *l.last().unwrap();
                                let n = 3 + b.bits(2);
                                l.extend(std::iter::repeat_n(prev, n));
                            }
                            17 => {
                                let n = 3 + b.bits(3);
                                l.extend(std::iter::repeat_n(0, n));
                            }
                            _ => {
                                let n = 11 + b.bits(7);
                                l.extend(std::iter::repeat_n(0, n));
                            }
                        }
                    }
                    (huff(&l[..hlit]), huff(&l[hlit..]))
                };
                loop {
                    match b.decode(&lit) {
                        s @ 0..=255 => out.push(s as u8),
                        256 => break,
                        s => {
                            let len = LBASE[s - 257] + b.bits(LEXT[s - 257]);
                            let k = b.decode(&dist);
                            let back = DBASE[k] + b.bits(DEXT[k]);
                            let from = out.len() - back;
                            for i in 0..len {
                                out.push(out[from + i]);
                            }
                        }
                    }
                }
            }
            _ => panic!("bad DEFLATE block type"),
        }
        if last == 1 {
            return out;
        }
    }
}

/// Pixel x of a PNG row.
fn px(row: &[u8], x: usize) -> [u8; 3] {
    row[1 + 3 * x..4 + 3 * x].try_into().unwrap()
}

/// How many pixels are red (k = 0) or blue (k = 2): that channel over the
/// other by more than 40.
fn count(path: &str, k: usize) -> usize {
    png_rows(path).2.iter().map(|r| r[1..].chunks(3).filter(|p| p[k] as i32 > p[2 - k] as i32 + 40).count()).sum()
}

/// A JSON array's numbers (0 and 0.0 alike).
fn floats(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect()
}

fn targets(d: &Value) -> Vec<&str> {
    d["legend"].as_array().unwrap().iter().map(|l| l["target"].as_str().unwrap()).collect()
}

#[test]
fn render_writes_a_png_size_across_and_says_what_res_it_meshed_at() {
    let s = Scratch::new("render-shelf");
    let png = s.path("shelf.png");
    let d = s.ink(["render", &fixture("shelf.ink"), "--out", &png, "--size", "300"]).ok();
    let (w, h, rows) = png_rows(&png);
    assert!((w, h) == (300, 300), "--size is the image's side: {w} x {h}");
    // ink writes RGB PNGs with unfiltered rows.
    assert!(rows.iter().all(|r| r[0] == 0), "rows are unfiltered");
    assert!(d["legend"].as_array().unwrap().iter().all(|l| l["res"].as_f64().unwrap() > 0.0 && l["unresolved"]["cells"].as_u64().is_some()), "{}", pretty(&d));
}

#[test]
fn render_without_a_target_or_root_draws_the_first_assembly_no_other_places() {
    // desk-lamp: five bodies, four lamps assembled from them, no root.
    // Drawn as the lamp, not every body at its origin; the others named.
    let s = Scratch::new("render-lamp");
    let d = s.ink(["render", &root("examples/desk-lamp.ink"), "--out", &s.path("lamp.png"), "--size", "128"]).ok();
    let legend: Vec<&str> = d["legend"].as_array().unwrap().iter().filter_map(|l| l["target"].as_str()).collect();
    assert!(d["target"] == "lamp" && legend.iter().all(|t| t.starts_with("lamp/")) && legend.len() == 4, "{}", pretty(&d));
    assert_eq!(d["other_assemblies"], serde_json::json!(["long-lamp", "swung-lamp", "swung-long-lamp"]));
}

#[test]
fn render_meshes_a_thin_long_part_finer_until_it_leaves_nothing_out() {
    // A 0.6 mm wall round a 200 x 150 mm loop (a belt's back): its size's res
    // (1.26) is coarser than the wall, which it drew dotted. Made finer, no
    // material is found left out; the cells holding both its faces (within
    // a cell of each other) are reported undecided, not shown drawn.
    let s = Scratch::new("render-band");
    let band = s.model("band.ink", "(import \"std\")\n(body band (shape (subtract (box 200 150 6) (box 198.8 148.8 8))))\n");
    let d = s.ink(["render", &band, "--out", &s.path("band.png")]).ok();
    let b = &d["legend"][0];
    assert!(b["unresolved"]["material"] == 0 && b["res"].as_f64().unwrap() < 1.0, "{}", pretty(b));
}

#[test]
fn render_meshes_a_thin_walled_tube_at_the_res_under_its_walls() {
    // A 2 x 1 in tube, 1/16 in walls (1678's drive rail, from STEP): at its
    // walls' corners the field's gradient runs along the other wall, which
    // read 4.72 mm across (it rendered blank). Its walls are 1.59 mm; no
    // material is found left out (cells holding both faces of a wall are
    // reported undecided).
    let s = Scratch::new("render-tube");
    let tube = s.model("tube.ink", "(import \"std\")\n(body tube (shape (subtract (box 50.8 25.4 635) (box 47.625 22.225 640))))\n");
    let d = s.ink(["render", &tube, "--out", &s.path("tube.png")]).ok();
    let b = &d["legend"][0];
    let res = b["res"].as_f64().unwrap();
    assert!(b["triangles"].as_u64().unwrap() > 0 && b["unresolved"]["material"] == 0 && 0.79 < res && res < 1.59, "{}", pretty(b));
}

#[test]
fn render_does_not_mesh_again_for_a_speck_left_out() {
    // A crest thinning to nothing leaves one cell out at every res, not
    // meshed again to 0.057 and 4.6 M triangles: a plate's face at 30
    // degrees to a 2.29 mm hole's axis meets the hole at a point (as on
    // 1678's id2726, from STEP).
    let s = Scratch::new("render-crest");
    let crest = s.model(
        "crest.ink",
        "(import \"std\")\n(body plate (shape (subtract (place (extrude-z (polygon (vec (vec -40 6.3374) (vec -20 17.8844) \
        (vec -18 17.8844) (vec -18 -18) (vec -40 -18))) -17 17) (at (rot -0.57735 -0.57735 -0.57735 120)))\n\
        (extrude-z (move (circle-2d 1.143) 12.7 0 0) -31 -18))))\n",
    );
    let d = s.ink(["render", &crest, "--out", &s.path("crest.png")]).ok();
    let b = &d["legend"][0];
    assert!(b["res"].as_f64().unwrap() > 0.2 && b["triangles"].as_u64().unwrap() < 500_000, "{}", pretty(b));
}

#[test]
fn render_transmission_reveals_a_part_inside_a_lid_and_follows_its_instance() {
    // A red box under a blue lid whose transmission is a field (x > 0 only),
    // the lid rotated by its instance.
    let s = Scratch::new("render-window");
    let window = s.model(
        "window.ink",
        "(import \"std\")\n(param transmission 0)\n\
        (body inside (shape (box 6 6 2)) (attr color \"#FF0000\"))\n\
        (body lid (shape (box 8 8 .5)) (attr color \"#0000FF\") (attr transmission (select x 0 transmission)))\n\
        (assembly root (instance inside inside) (instance lid lid (at (rot-z 180) (move 0 0 2.5))))\n",
    );
    let mut pixels = vec![];
    for tr in ["0", "0.78", "1"] {
        let image = s.path(&format!("window-{tr}.png"));
        s.ink(["render", &window, "--out", &image, "--size", "128", "--res", ".5", "--param", &format!("transmission={tr}")]).ok();
        let rows = png_rows(&image).2;
        pixels.push([px(&rows[32], 86), px(&rows[32], 106)]);
    }
    assert!(pixels[0][0][0] == 0 && 0 < pixels[1][0][0] && pixels[1][0][0] < pixels[2][0][0] && pixels[2][0][2] == 0, "reveals the red part: {pixels:?}");
    assert!(pixels[0][1] == pixels[1][1] && pixels[1][1] == pixels[2][1] && pixels[0][1][2] > 0, "follows the rotated instance: {pixels:?}");
}

#[test]
fn render_warns_where_walls_thinner_than_res_fall_between_samples() {
    let s = Scratch::new("render-thin");
    let thin = s.model("thin.ink", "(import \"std\")\n(body cup (shape (subtract (cylinder-z 40 60) (move (cylinder-z 39.6 60) 0 0 1))))\n");
    let d = s.ink(["render", &thin, "--out", &s.path("thin.png"), "--size", "128", "--res", "0.5"]).ok();
    assert!(d["legend"][0]["unresolved"]["cells"].as_u64().unwrap() > 0 && warned(&d, &["thinner than"]), "{}", pretty(&d));
}

#[test]
fn render_warns_of_a_ball_in_a_cell_the_surface_crosses() {
    // A ball 0.3 across just outside a cube's face, in a cell at res 1
    // the face crosses: not drawn, so warned of.
    let s = Scratch::new("render-crossed-ball");
    let m = s.model("ball.ink", "(import \"std\")\n(body b (shape (min (box 6.4 6.4 6.4) (move (sphere .15) 3.45 .2 .2))))\n");
    let d = s.ink(["render", &m, "--out", &s.path("ball.png"), "--res", "1", "--size", "128", "--view", "front"]).ok();
    assert!(warned(&d, &["'b'"]), "{}", pretty(&d));
}

#[test]
fn render_reports_and_warns_of_material_found_beside_its_surface() {
    // A ball 1.4 across between two slabs, beside their faces: the search
    // finds it, which is reported and warned of, however thick.
    let s = Scratch::new("render-thick-ball");
    let m = s.model("ball.ink", "(import \"std\")\n(body b (shape (min (move (box 2.6 6.4 6.4) -1.9 0 0) (move (box 2.2 6.4 6.4) 2.1 0 0) (move (sphere .7) .2 .2 .2))))\n");
    let d = s.ink(["render", &m, "--out", &s.path("ball.png"), "--res", "1", "--size", "128", "--view", "front"]).ok();
    assert!(d["legend"][0]["unresolved"]["material"].as_u64().unwrap() > 0 && warned(&d, &["'b'", "leaves out"]), "{}", pretty(&d));
}

#[test]
fn render_meshes_again_for_a_foil_seen_face_on() {
    // A 6 x 6 x 0.015 foil beside a 6.4 mm cube, seen from the top in a
    // small picture: far thinner than a pixel, but 6 mm across, so drawn,
    // as many blue pixels as a 0.5 mm plate there shows (a block past both
    // frames the two pictures alike).
    let s = Scratch::new("render-foil-face");
    let blue = |thick: f64| {
        let m = s.model(
            "foil.ink",
            &format!("(import \"std\")\n(body b (shape (min (box 6.4 6.4 6.4) (move (box 6 6 {thick}) 7 0 0.05) (move (box 1 1 1) 11 0 0))) (attr color \"#0000FF\"))\n"),
        );
        let out = s.path("foil.png");
        let d = s.ink(["render", &m, "--out", &out, "--size", "128", "--view", "top"]).ok();
        (count(&out, 2), d)
    };
    let ((plate, _), (foil, d)) = (blue(0.5), blue(0.015));
    assert!(foil as f64 > 0.95 * plate as f64, "{foil} blue pixels, with a plate {plate}: {}", pretty(&d));
}

#[test]
fn render_meshes_a_body_by_its_box_and_a_glb_by_its_surface() {
    // A 200 x 200 x 1 plate: by its surface (a GLB's, for a viewer's zoom)
    // its res is 0.71, by its box (a picture's, as at 0.1.2) 1: half the
    // cells.
    let s = Scratch::new("render-plate");
    let plate = s.model("plate.ink", "(import \"std\")\n(body plate (shape (box 200 200 1)))\n");
    let r = s.ink(["render", &plate, "--out", &s.path("plate.png")]).ok();
    let g = s.ink(["mesh", &plate, "--out", &s.path("plate.glb")]).ok();
    let (picture, glb) = (r["legend"][0]["res"].as_f64().unwrap(), g["objects"][0]["res"].as_f64().unwrap());
    assert!(picture >= 1.4 * glb, "render res {picture}, GLB res {glb}: {}\n{}", pretty(&r["legend"]), pretty(&g["objects"]));
}

#[test]
fn a_picture_draws_a_thin_wall_seen_face_on() {
    // A cup's 0.4 mm wall, seen face on: at 512 pixels across drawn whole;
    // at 256, about a pixel thick, drawn or warned of.
    let s = Scratch::new("render-cup");
    let cup = s.model("cup.ink", "(import \"std\")\n(body cup (shape (subtract (cylinder-z 40 60) (move (cylinder-z 39.6 60) 0 0 1))))\n");
    for size in ["256", "512"] {
        let d = s.ink(["render", &cup, "--out", &s.path("cup.png"), "--size", size]).ok();
        let drawn = d["legend"][0]["res"].as_f64().unwrap() < 0.4 && d["legend"][0]["unresolved"]["material"] == 0;
        assert!(drawn || (size == "256" && warned(&d, &["'cup'", "leaves out"])), "--size {size}: {}", pretty(&d));
    }
}

#[test]
fn shapes_that_only_touch_intersected_have_no_interior_and_render_nothing() {
    // (max a b) written for a union.
    let s = Scratch::new("render-touch");
    let touch = s.model("touch.ink", "(import \"std\")\n(body tt (shape (max (move (box 14 14 12) 0 0 6) (move (box 36 30 2.5) 0 -12 13.25))))\n");
    let d = s.ink(["check", &touch]).ok();
    let w = d["warnings"].as_array().unwrap();
    assert!(w.len() == 1 && warned(&d, &["no interior", "(min a b)"]), "{}", pretty(&d));
    let d = s.ink(["render", &touch, "--out", &s.path("touch.png"), "--size", "64"]).ok();
    assert!(warned(&d, &["no triangles"]), "{}", pretty(&d));
}

/// A blue box hiding a red ball, and a green block above it all.
const HIDDEN: &str = "(import \"std\")\n(body shell (shape (subtract (box 20 20 20) (sphere 8))) (attr color \"#0000FF\"))\n\
(body ball (shape (sphere 6)) (attr color \"#FF0000\"))\n(body lid (shape (move (box 4 4 4) 15 15 20)) (attr color \"#00FF00\"))\n";

#[test]
fn render_view_picks_one_view_by_name_or_direction() {
    let s = Scratch::new("render-view");
    let hidden = s.model("hidden.ink", HIDDEN);
    let out = s.path("view.png");
    let d = s.ink(["render", &hidden, "--out", &out, "--view", "back", "--size", "64", "--res", ".5"]).ok();
    let views: Vec<&Value> = d["views"].as_array().unwrap().iter().map(|v| &v["view"]).collect();
    assert!(d["width"] == 64 && d["height"] == 64 && views == [&json!("back")] && floats(&d["views"][0]["camera_from"]) == [0.0, 1.0, 0.0], "{}", pretty(&d));
    let d = s.ink(["render", &hidden, "--out", &out, "--view", "[0,0,-1]", "--size", "64", "--res", ".5"]).ok();
    assert!(floats(&d["views"][0]["camera_from"]) == [0.0, 0.0, -1.0] && d["views"][0]["view"] == "from -z", "{}", pretty(&d["views"]));
    let e = s.ink(["render", &hidden, "--out", &out, "--view", "side"]).err();
    assert!(e["error"]["message"].as_str().unwrap().contains("iso, top, front, right, back, left, bottom"), "an unknown view lists the names: {}", pretty(&e));
}

#[test]
fn render_cut_removes_the_half_on_the_cameras_side() {
    // The ball shows through the cut box; the lid, wholly in the removed
    // half, is left out. Uncut, the box hides the ball.
    let s = Scratch::new("render-cut");
    let hidden = s.model("hidden.ink", HIDDEN);
    let centre = |path: &str| px(&png_rows(path).2[32], 33);
    let cut = s.path("cut.png");
    let d = s.ink(["render", &hidden, "--out", &cut, "--view", "top", "--cut", "z=0", "--size", "64", "--res", ".5"]).ok();
    let c = centre(&cut);
    assert!(d["cut"]["removed"] == "z > 0" && targets(&d) == ["shell", "ball"] && c[0] as u32 > 2 * c[2] as u32, "centre {c:?}: {}", pretty(&d));
    let whole = s.path("whole.png");
    s.ink(["render", &hidden, "--out", &whole, "--view", "top", "--size", "64", "--res", ".5"]).ok();
    let c = centre(&whole);
    assert!(c[2] as u32 > 2 * c[0] as u32, "uncut, the box hides it: centre {c:?}");
}

#[test]
fn render_cut_at_a_shared_face_leaves_out_a_target_with_nothing_left() {
    // A cap standing on a block, cut at their shared face, z = 5, also by
    // default: the middle of their bounds, which are rounded outward by
    // different cells (not 4.944). Nothing of the cap is left.
    let s = Scratch::new("render-capped");
    let capped = s.model("capped.ink", "(import \"std\")\n(body block (shape (box 40 40 10)))\n(body cap (shape (move (box 20 20 10) 0 0 10)))\n");
    for cut in ["z", "z=5"] {
        let d = s.ink(["render", &capped, "--out", &s.path("capped.png"), "--view", "top", "--cut", cut, "--size", "64"]).ok();
        assert!(d["cut"]["removed"] == "z > 5" && targets(&d) == ["block"] && d["warnings"].as_array().is_none_or(Vec::is_empty), "--cut {cut}: {}", pretty(&d));
    }
}

#[test]
fn render_cut_keeps_a_crossed_target_empty_for_being_thin_and_warns() {
    // A foil thinner than res across the plane, not left out like the cap.
    let s = Scratch::new("render-foil");
    let foil = s.model("foil.ink", "(import \"std\")\n(body block (shape (box 40 40 10)))\n(body foil (shape (move (box 30 0.05 20) 0 25 0)))\n");
    let d = s.ink(["render", &foil, "--out", &s.path("foil.png"), "--view", "top", "--cut", "z=0", "--res", "1", "--size", "64"]).ok();
    assert!(targets(&d) == ["block", "foil"] && warned(&d, &["'foil'", "thinner than res"]), "{}", pretty(&d));
}

#[test]
fn render_cut_keeps_a_target_with_a_stub_on_the_side_kept() {
    // A 1 m rod with 3 mm kept: its section is drawn (a margin of a bounds
    // cell, 4 mm here, left it out).
    let s = Scratch::new("render-rod");
    let rod = s.model("rod.ink", "(import \"std\")\n(body block (shape (move (box 40 40 10) 300 0 -20)))\n(body rod (shape (move (box 40 40 1000) 0 0 497)))\n");
    let d = s.ink(["render", &rod, "--out", &s.path("rod.png"), "--view", "top", "--cut", "z=0", "--size", "64"]).ok();
    let got: Vec<(&str, bool)> = d["legend"].as_array().unwrap().iter().map(|l| (l["target"].as_str().unwrap(), l["triangles"].as_u64().unwrap() > 0)).collect();
    assert!(got == [("block", true), ("rod", true)], "{}", pretty(&d["legend"]));
}

#[test]
fn render_gives_every_target_its_own_colour_where_the_palette_allows() {
    // Eight instances, seven of one body (one mesh shared): each target
    // without a colour of its own is drawn in its own.
    let s = Scratch::new("render-colours");
    let mut src = String::from("(import \"std\")\n(body cube (shape (box 4 4 4)))\n(body plate (shape (move (box 60 6 2) 0 0 -4)))\n(assembly root\n  (instance p plate)\n");
    for k in 0..7 {
        src += &format!("  (instance c{k} cube (at (move {} 0 0)))\n", 8 * k - 24);
    }
    src += ")\n";
    let m = s.model("cubes.ink", &src);
    let d = s.ink(["render", &m, "--out", &s.path("cubes.png"), "--size", "64"]).ok();
    let colours: Vec<&str> = d["legend"].as_array().unwrap().iter().map(|l| l["color"].as_str().unwrap()).collect();
    let distinct: std::collections::HashSet<&&str> = colours.iter().collect();
    assert!(colours.len() == 8 && distinct.len() == 8, "{}", pretty(&d["legend"]));
}

#[test]
fn render_coincident_surfaces_show_the_target_listed_first() {
    // A gyroid lattice filling a block, both whole: where the lattice meets
    // the block's faces its mesh crosses them by its own error, and the
    // block, listed first, shows all over (listed second, the lattice shows
    // in 20,614 pixels). Seen face on, both at res 0.125 (a 16 mm block:
    // 30 mm takes 4.4 s alone, 21 s beside the other tests): along the edges
    // between faces the lattice's mesh can stand a pixel past the block's,
    // and coarser (0.149, the lattice's own res) its own error pokes past
    // the block's top in a few dozen pixels.
    let s = Scratch::new("render-shared");
    let shared = s.model(
        "shared.ink",
        "(import \"std\")\n(body solid (shape (box 16 16 16)) (attr color \"#FF0000\"))\n\
        (body lattice (shape (lattice-fill (box 16 16 16) (gyroid 8) 1)) (attr color \"#0000FF\"))\n",
    );
    let image = s.path("shared.png");
    s.ink(["render", &shared, "--out", &image, "--size", "256", "--view", "top", "--res", "0.125"]).ok();
    let (blue, red) = (count(&image, 2), count(&image, 0));
    assert!(blue == 0 && red > 1000, "blue {blue}, red {red}");
}

#[test]
fn render_a_cover_just_proud_of_a_block_shows_in_a_large_scene() {
    // 0.8 mm proud, in a scene 2 m across.
    let s = Scratch::new("render-cover");
    let cover = s.model(
        "cover.ink",
        "(import \"std\")\n(body block (shape (box 100 100 10)))\n(body cover (shape (box 60 60 4)) (attr color \"#FF0000\"))\n\
        (body far (shape (box 10 10 10)))\n\
        (assembly root (instance block block) (instance cover cover (at (move 0 0 3.8))) (instance far far (at (move 1900 0 0))))\n",
    );
    let image = s.path("cover.png");
    s.ink(["render", cover.as_str(), "--out", image.as_str(), "--view", "top", "--size", "256"]).ok();
    let red = count(&image, 0);
    assert!(red > 20, "{red} red pixels");
}
