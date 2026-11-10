//! The benchmarks' key answers that are not claims of their files: probes
//! of a solve against hand calculations, interference summaries, claims at
//! other params or moments, estimates of volumes and masses, meshes. A
//! benchmark's own claims (the hand numbers the language can state:
//! masses, volumes, clearances) are checked by the model rule
//! (models.rs). One test per benchmark, named `bNN_...` so
//! `benchmarks/run.sh NN` can pick it; the slow ones (models.rs's SLOW) are
//! ignored here and run by benchmarks/run.sh. Each file's
//! header gives the hand calculation and where its numbers come from. A
//! test's solve is mostly the one its file's claims ask (the same study),
//! so the model rule and the test compute it once between them (cached).

use crate::common::*;
use serde_json::Value;

fn bench(name: &str) -> String {
    root(&format!("benchmarks/{name}"))
}

fn num(v: &Value) -> f64 {
    v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"))
}

/// `v` is within relative `tol` of `want`.
#[track_caller]
fn near(what: &str, v: f64, want: f64, tol: f64) {
    assert!((v - want).abs() <= tol * want.abs(), "{what}: {v}, want {want} within {}%", tol * 100.0);
}

/// A mass answer's target by name.
fn target<'a>(d: &'a Value, name: &str) -> &'a Value {
    d["targets"].as_array().unwrap().iter().find(|t| t["target"] == name || t["name"] == name).unwrap_or_else(|| panic!("no target {name}:\n{}", pretty(d)))
}

/// An interfere answer: no pair overlaps outside its declared contacts and
/// none is left undecided.
#[track_caller]
fn no_interference(d: &Value) {
    let unknown: Vec<&Value> = d["pairs"].as_array().unwrap().iter().filter(|p| p["status"] == "unknown").collect();
    assert!(d["overlapping_pairs"] == 0 && unknown.is_empty(), "overlapping pairs {}, undecided {unknown:?}:\n{}", d["overlapping_pairs"], pretty(d));
}

/// An assert answer's claims by name.
fn claims(d: &Value) -> std::collections::BTreeMap<String, Value> {
    d["asserts"].as_array().unwrap().iter().map(|a| (a["name"].as_str().unwrap().to_string(), a.clone())).collect()
}

/// A solved mesh: watertight, one shell.
#[track_caller]
fn one_watertight_shell(d: &Value) {
    for o in d["objects"].as_array().unwrap() {
        assert!(o["watertight"] == true && o["shells"] == 1, "mesh not one watertight shell:\n{}", pretty(d));
    }
}

#[test]
fn the_readme_lists_every_benchmark_and_names_each_number_missing() {
    // Its table has a row per benchmarks/NN-*.ink, and each number up to
    // the last that has none (one gone, not reused) is named in its text.
    let readme = std::fs::read_to_string(bench("README.md")).unwrap();
    let mut have: Vec<u32> = std::fs::read_dir(bench(""))
        .unwrap()
        .filter_map(|e| {
            let name = e.ok()?.file_name().into_string().ok()?;
            name.ends_with(".ink").then(|| name.split('-').next()?.parse().ok()).flatten()
        })
        .collect();
    have.sort();
    for n in &have {
        assert!(readme.contains(&format!("\n| {n:02} |")), "benchmarks/README.md has no row for {n:02}");
    }
    let (text, _) = readme.split_once("| # |").unwrap();
    for n in (1..*have.last().unwrap()).filter(|n| !have.contains(n)) {
        assert!(text.contains(&format!(" {n} (")), "benchmarks/README.md does not say where {n} went");
    }
}

#[test]
fn b06_plate_mass_properties_match_the_polygon_moments() {
    // ink mass weighs the plate by its m (rho, --density m): its mass,
    // centre of mass and inertia about it, against the exact polygon
    // moments of benchmark 04's outline (8 thick, rho 0.0027): mass within
    // 0.1%, COM within 0.001 mm, inertia within 0.1% of the largest diagonal
    // term. (The same as claims holds, but proving the inertia that close
    // takes ten times this estimate, and the COM to 0.001 mm is left
    // undecided.)
    let d = ink(["mass", &bench("06-mass-properties.ink"), "--target", "plate-mp", "--density", "m", "--res", "0.25", "--keys", "mass,center_of_mass,inertia,volume"]);
    let t = &d["targets"][0];
    let m = num(&t["mass"]);
    let c: Vec<f64> = (0..3).map(|k| num(&t["center_of_mass"][k])).collect();
    let i: Vec<Vec<f64>> = (0..3).map(|r| (0..3).map(|s| num(&t["inertia"][r][s])).collect()).collect();
    // shoelace sums of the outline: A, Sx, Sy, Sxx, Syy, Sxy
    let p = [(0., 0.), (90., 0.), (90., 20.), (70., 20.), (70., 12.), (55., 12.), (40., 35.), (40., 70.), (22., 70.), (22., 50.), (10., 60.), (0., 40.)];
    let (mut a, mut sx, mut sy, mut sxx, mut syy, mut sxy) = (0., 0., 0., 0., 0., 0.);
    for i in 0..p.len() {
        let ((x0, y0), (x1, y1)) = (p[i], p[(i + 1) % p.len()]);
        let k = x0 * y1 - x1 * y0;
        a += k / 2.;
        sx += (x0 + x1) * k / 6.;
        sy += (y0 + y1) * k / 6.;
        sxx += (x0 * x0 + x0 * x1 + x1 * x1) * k / 12.;
        syy += (y0 * y0 + y0 * y1 + y1 * y1) * k / 12.;
        sxy += (x0 * y1 + 2. * x0 * y0 + 2. * x1 * y1 + x1 * y0) * k / 24.;
    }
    let (rho, h) = (0.0027, 8.0);
    // the raw moments m, mx, my, mz, mxx, myy, mzz, mxy, mxz, myz
    let e = [a * h, sx * h, sy * h, a * h * h / 2., sxx * h, syy * h, a * h * h * h / 3., sxy * h, sx * h * h / 2., sy * h * h / 2.].map(|v| rho * v);
    // mass, COM, inertia about the COM
    let me = e[0];
    let ce = [e[1] / me, e[2] / me, e[3] / me];
    let io = [[e[5] + e[6], -e[7], -e[8]], [-e[7], e[4] + e[6], -e[9]], [-e[8], -e[9], e[4] + e[5]]];
    let c2: f64 = ce.iter().map(|v| v * v).sum();
    let ie: Vec<Vec<f64>> = (0..3).map(|r| (0..3).map(|s| io[r][s] - me * (if r == s { c2 } else { 0. } - ce[r] * ce[s])).collect()).collect();
    let com = (0..3).map(|k| (c[k] - ce[k]).abs()).fold(0., f64::max);
    let inertia = (0..9).map(|k| (i[k / 3][k % 3] - ie[k / 3][k % 3]).abs()).fold(0., f64::max) / (0..3).map(|k| ie[k][k].abs()).fold(0., f64::max);
    assert!(
        (m - me).abs() / me < 1e-3 && com < 1e-3 && inertia < 1e-3,
        "mass {m} (exact {me}), COM {c:?} ({ce:?}), inertia {i:?} ({ie:?}): errors {}, {com} mm, {inertia}",
        (m - me).abs() / me
    );
}

#[test]
fn b07_flange_parts_neither_overlap_nor_stay_undecided() {
    // At cells of 0.2: the flanges' faces touch over 6600 mm2, past the
    // search's budget at the model's floor.
    no_interference(&ink(["interfere", &bench("07-bolted-flange.ink"), "--res", "0.2"]));
}

#[test]
fn b09_housing_overlaps_only_in_its_press_fits() {
    // At cells of 0.1: the contact zones' surfaces touch the parts' over more
    // than the budget covers at the model's floor. (The mass is a claim.)
    no_interference(&ink(["interfere", &bench("09-bearing-housing.ink"), "--res", "0.1"]));
}

#[test]
fn b19_sweeps_match_area_times_length_and_the_loft_its_ruled_sections() {
    // The duct and the taper: estimates within 1% of area (mean area) x path
    // length, the expected volumes the file computes (ink check's defs). The frustum's
    // volume is a claim. The transition's ends exact within 1% (circle r 20,
    // square 40), its middle within 3% of a ruled loft's 1419.3; the duct
    // meshes as one watertight shell.
    let f = bench("19-spline-sweep-loft.ink");
    let check = ink(["check", &f]);
    let def = |n: &str| num(&check["defs"].as_array().unwrap().iter().find(|x| x["name"] == n).unwrap()["value"]);
    for (body, want) in [("duct", def("duct-volume-expected")), ("taper", def("taper-volume-expected"))] {
        let v = ink(["mass", &f, "--target", body, "--res", "0.25"]);
        near(&format!("{body} volume"), num(&target(&v, body)["volume"]["estimate"]), want, 0.01);
    }
    for (at, want, tol) in [("0.01", 1256.64, 0.01), ("49.99", 1600.0, 0.01), ("25", 1419.3, 0.03)] {
        let s = ink(["slice", &f, "--target", "transition", "--axis", "z", "--at", at]);
        near(&format!("transition section at z {at}"), num(&s["area"]), want, tol);
    }
    let dir = Scratch::new("b19");
    one_watertight_shell(&dir.ink(["mesh", &f, "--target", "duct", "--out", &dir.path("duct.stl")]).ok());
}

#[test]
fn b31_bracket_is_the_formed_volume_and_moved_holes_break_the_rules() {
    // Volume and mass at cells of 0.1: estimates within 0.1% of the formed
    // (K = 0.5) 7277.07 mm3 and 19.648 g, and 1.143% over the K = 0.33
    // blank, 7194.82 (within 1.5%; 0.5% would fail). The ligament rules fail
    // with a hole moved toward the free edge or the bend.
    let f = bench("31-sheet-bracket.ink");
    let t = &ink(["mass", &f, "--res", "0.1"])["targets"][0];
    near("volume vs formed", num(&t["volume"]["estimate"]), 7277.0707, 0.001);
    near("volume vs blank", num(&t["volume"]["estimate"]), 7194.8238, 0.015);
    near("mass", num(&t["mass"]), 19.64809, 0.001);
    for (rule, param) in [("hole-edge-distance", "hole-z=19"), ("hole-bend-distance", "hole-z=11")] {
        let s = statuses(&ink(["assert", &f, "--res", "0.1", "--name", rule, "--param", param]));
        assert_eq!(s.get(rule).map(String::as_str), Some("fails"), "{rule} at {param}: {s:?}");
    }
}

#[test]
fn b34_fit_controls_fail_in_band_and_hold_at_the_middles() {
    // Each control's witness lies in the bands and fails it: the 10 um
    // clearance at sizes with (D - d) / 2 < 0.010, the 20 um gap at D - d <=
    // 0.020, the loose spec at d > D; at the band middles (8.011, 7.9795) the
    // 10 um clearance holds; interfere: the press fit is the only overlap.
    let f = bench("34-tolerance-fit.ink");
    let a = claims(&ink(["assert", &f]));
    let p = |c: &str, k: &str| num(&a[c]["params"][k]);
    let (bore, shaft) = (p("running-clear-10um", "bore-d"), p("running-clear-10um", "shaft-d"));
    assert!(
        (8.0..=8.022).contains(&bore) && (7.972..=7.987).contains(&shaft) && (bore - shaft) / 2.0 < 0.010 && num(&a["running-clear-10um"]["value"]) < 0.010,
        "{}",
        pretty(&a["running-clear-10um"])
    );
    assert!(p("running-gap-20um", "bore-d") - p("running-gap-20um", "shaft-d") <= 0.02, "{}", pretty(&a["running-gap-20um"]));
    assert!(p("loose-apart", "loose-shaft-d") > p("loose-apart", "loose-bore-d"), "{}", pretty(&a["loose-apart"]));
    assert!(a["running-overlaps"]["params"].is_object(), "{}", pretty(&a["running-overlaps"]));
    let m = statuses(&ink(["assert", &f, "--name", "running-clear-10um", "--param", "bore-d=8.011:8.011,shaft-d=7.9795:7.9795"]));
    assert_eq!(m["running-clear-10um"], "holds", "{m:?}");
    no_interference(&ink(["interfere", &f]));
}

#[test]
fn b35_tight_stack_fails_at_an_in_band_corner_and_holds_at_the_values() {
    // Witnesses with every length in its band and a + b + c > G2 (>= for the
    // gap), tight-apart's point between G2 and a + b + c along x; with
    // --nominal both hold: the tolerances are what fails them.
    let f = bench("35-tolerance-stack.ink");
    let a = claims(&ink(["assert", &f]));
    let band = |c: &str| {
        let ps = &a[c]["params"];
        // Written exactly: a band's end is the float it computes to.
        ["a", "b", "c"].iter().all(|k| (9.9 - 1e-9..=10.1 + 1e-9).contains(&num(&ps[*k]))) && (30.1 - 1e-9..=30.3 + 1e-9).contains(&num(&ps["G2"]))
    };
    let sum = |c: &str| ["a", "b", "c"].iter().map(|k| num(&a[c]["params"][*k])).sum::<f64>();
    let (g2, x) = (num(&a["tight-apart"]["params"]["G2"]), num(&a["tight-apart"]["witness"][0]));
    assert!(
        band("tight-apart") && band("tight-gap") && sum("tight-apart") > g2 && g2 < x && x < sum("tight-apart") && sum("tight-gap") >= num(&a["tight-gap"]["params"]["G2"]),
        "{}\n{}",
        pretty(&a["tight-apart"]),
        pretty(&a["tight-gap"])
    );
    let n = statuses(&ink(["assert", &f, "--nominal"]));
    assert!(n["tight-apart"] == "holds" && n["tight-gap"] == "holds", "{n:?}");
}

#[test]
fn b37_tight_stop_fails_at_an_in_band_corner_and_holds_at_the_values() {
    // Witnesses with L and h in their bands where L^2 - h^2 is past the hand
    // value (x-under: x > 40.15, so > 1612.02; tight-apart: > 1608.01); with
    // --nominal both hold, as estimates (the solve is numerical there).
    let f = bench("37-solved-tolerance.ink");
    let a = claims(&ink(["assert", &f]));
    for (c, past) in [("x-under", 1612.02), ("tight-apart", 1608.01)] {
        let (l, h) = (num(&a[c]["params"]["L"]), num(&a[c]["params"]["h"]));
        assert!((49.9..=50.1).contains(&l) && (29.9..=30.1).contains(&h) && l * l - h * h > past, "{c}: {}", pretty(&a[c]));
    }
    let n = statuses(&ink(["assert", &f, "--nominal"]));
    assert!(n["tight-apart"] == "holds estimate" && n["x-under"] == "holds estimate", "{n:?}");
}

// ---------------------------------------------------------------- slow: benchmarks/run.sh

/// The first hit of a ray from `from` along `dir` on `target`.
fn ray(file: &str, target: &str, from: [f64; 3], dir: [f64; 3]) -> Value {
    let d = ink(["ray", file, "--target", target, "--from", &format!("{from:?}"), "--dir", &format!("{dir:?}")]);
    let hit = if d["nearest"].is_object() { d["nearest"].clone() } else { d["hits"][0].clone() };
    assert!(hit.is_object(), "ray from {from:?} along {dir:?} misses {target}:\n{}", pretty(&d));
    hit
}

#[test]
#[ignore = "slow benchmark: benchmarks/run.sh"]
fn b10_slider_crank_parts_stay_apart_over_a_turn() {
    no_interference(&ink(["interfere", &bench("10-slider-crank.ink"), "--t-span", "0:1"]));
}

#[test]
#[ignore = "slow benchmark: benchmarks/run.sh"]
fn b11_gearbox_meshes_stay_apart_over_a_pitch() {
    // Proven continuously over t: the meshing pairs hold to the accuracy
    // (their teeth come within 0.05 by design), the pair that never meets
    // holds outright.
    let d = ink(["interfere", &bench("11-gearbox.ink"), "--t-span", "0:1"]);
    no_interference(&d);
    let io = d["pairs"].as_array().unwrap().iter().find(|p| p["pair"][0].as_str().unwrap().ends_with("input") && p["pair"][1].as_str().unwrap().ends_with("output"));
    assert!(io.is_some_and(|p| verdict(p) == "holds"), "input~output: {io:?}");
}

#[test]
#[ignore = "slow benchmark: benchmarks/run.sh"]
fn b18_ten_thousand_holes_two_ways_and_ten_thousand_pins() {
    // Both screens within 0.5% of 200 x 200 x 2 - 10000 pi 0.6^2 2 = 57 380.5
    // mm3 and watertight at 8 cells across a hole; 10 000 pins in the
    // periodic screen's holes (instances are top-level forms, so the
    // assembly is written here) collide with nothing.
    let f = bench("18-large-pattern.ink");
    let v = ink(["mass", &f, "--res", "0.2"]);
    for t in ["screen-explicit", "screen-periodic"] {
        near(t, num(&target(&v, t)["volume"]["estimate"]), 57380.5, 0.005);
    }
    let dir = Scratch::new("b18");
    one_watertight_shell(&dir.ink(["mesh", &f, "--res", "0.25", "--out", &dir.path("screen.3mf")]).ok());
    let mut pins = format!("(import \"{f}\")\n(body pin (shape (cylinder-z 0.5 4)))\n(assembly root (instance screen screen-periodic)\n");
    for i in 0..100 {
        for j in 0..100 {
            pins += &format!("  (instance p{i}-{j} pin (at (move {} {} 0)))\n", -99 + 2 * i, -99 + 2 * j);
        }
    }
    pins += ")\n";
    no_interference(&dir.ink(["interfere", &dir.model("pins.ink", &pins)]).ok());
}

#[test]
#[ignore = "slow benchmark: benchmarks/run.sh"]
fn b20_redistanced_shell_wall_is_2_mm_and_the_naive_one_is_not() {
    // The solid meshes watertight; the redistanced shell is one solid round
    // one cavity; along rays from inside the cavity its wall's normal
    // thickness is 2 within 0.1, the naive shell's off by more than 0.25.
    let f = bench("20-organic-shell.ink");
    let dir = Scratch::new("b20");
    one_watertight_shell(&dir.ink(["mesh", &f, "--target", "mouse", "--out", &dir.path("mouse.stl")]).ok());
    let c = &ink(["components", &f, "--target", "mouse-shell"])["targets"][0];
    assert!(c["solids"] == 1 && c["cavities"] == 1, "{}", pretty(c));
    // inner wall: the first hit from c; outer: the first coming back from
    // 200 along; scaled by |cos| to the outer normal
    let wall = |target: &str, c: [f64; 3]| {
        [[0.0f64, 0.0, 1.0], [1.0, 0.0, 0.3]]
            .iter()
            .map(|d| {
                let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                let d = d.map(|v| v / n);
                let inner = num(&ray(&f, target, c, d)["t"]);
                let outer = ray(&f, target, [0, 1, 2].map(|k| c[k] + 200.0 * d[k]), d.map(|v| -v));
                let cos = (0..3).map(|k| num(&outer["normal"][k]) * d[k]).sum::<f64>().abs();
                ((200.0 - num(&outer["t"]) - inner) * cos - 2.0).abs()
            })
            .fold(0.0, f64::max)
    };
    let (exact, naive) = (wall("mouse-shell", [-5.0, 120.0, 12.0]), wall("naive-shell", [-5.0, 240.0, 12.0]));
    assert!(exact < 0.1 && naive > 0.25, "wall off 2 mm by {exact} redistanced, {naive} naive");
}

#[test]
#[ignore = "slow benchmark: benchmarks/run.sh"]
fn b22_imported_plate_measures_as_its_source_offset_and_shell() {
    // The supplier's STL, inputs/22-part.stl (tools/22-part-source.ink
    // meshed). Volume within 0.3% of 20 984.07; by rays the top at z 10, the
    // counterbore floor at 6, the offset's top at 11, the naive offset of 3f
    // at 10.333 and the redistanced one at 11; by slices the bore r 8 and the
    // offset's r 7; the shell one solid round one cavity.
    let f = bench("22-import-offset-measure.ink");
    near("volume", num(&ink(["mass", &f, "--res", "0.2", "--target", "part"])["targets"][0]["volume"]["estimate"]), 20984.07, 0.003);
    for (t, from, want, tol) in [
        ("part", [20.0, 15.0, 50.0], 10.0, 0.002),
        ("part", [10.0, 0.0, 50.0], 6.0, 0.004),
        ("grown", [120.0, 15.0, 50.0], 11.0, 0.003),
        ("naive", [20.0, 115.0, 50.0], 10.3333, 0.002),
        ("exact", [20.0, 215.0, 50.0], 11.0, 0.003),
    ] {
        near(&format!("{t} from {from:?}"), num(&ray(&f, t, from, [0.0, 0.0, -1.0])["point"][2]), want, tol);
    }
    let bore = |t: &str| {
        let s = ink(["slice", &f, "--target", t, "--axis", "z", "--at", "2"]);
        (-records(&s["loops"]).iter().map(|l| num(&l["area"])).fold(f64::INFINITY, f64::min) / std::f64::consts::PI).sqrt()
    };
    near("bore", bore("part"), 8.0, 0.0025);
    near("offset bore", bore("grown"), 7.0, 0.004);
    let c = &ink(["components", &f, "--target", "hollow"])["targets"][0];
    assert!(c["solids"] == 1 && c["cavities"] == 1, "{}", pretty(c));
}

#[test]
#[ignore = "slow benchmark: benchmarks/run.sh"]
fn b26_enclosure_wall_is_2_mm_without_false_thin_spots() {
    // The walls claim's value: the ribs' 1.2 minimum (5%; the file claims
    // it above 1.1).
    let a = &ink(["assert", &bench("26-enclosure.ink"), "--name", "walls"])["asserts"][0];
    near("min thickness", num(&a["value"]), 1.2, 0.05);
}

#[test]
#[ignore = "slow benchmark: benchmarks/run.sh"]
fn b32_four_bar_solve_matches_the_closed_form() {
    // tools/32-solved.ink solves the loop at each moment (a t-dependent solve
    // blocks the file's swept claims, so it lives there); at six moments
    // (the extremes of the rocker and transmission among them) the solve and
    // the closed form agree, as an estimate.
    let f = bench("tools/32-solved.ink");
    for t in ["0", "0.11611", "0.25", "0.5", "0.66387", "0.75"] {
        let s = statuses(&ink(["assert", &f, "--t", t, "--name", "solved-matches-closed-form"]));
        assert_eq!(s.get("solved-matches-closed-form").map(String::as_str), Some("holds estimate"), "t {t}: {s:?}");
    }
}

#[test]
#[ignore = "slow benchmark: benchmarks/run.sh"]
fn b33_imported_knob_keeps_its_measured_error_volume_and_knurl() {
    // The 58k-triangle STL, inputs/33-knob.stl, is what tools/33-gen-knob.py
    // writes at n 240 (written again here and compared). The
    // import's measured deviation is below the sampling bound sqrt 3 / 2 h
    // plus each sample's certified error, 0.2165; its volume within 0.2% of
    // the mesh's 43 663.15 (the samples' own is 0.1% under it: 43 619.5 at
    // res 0.15 and 0.1, rounding the knurl's tips and the fins' edges);
    // remeshed,
    // watertight, one shell, the same volume; rays inward at the 240 knurl
    // vertices (z 10) see tips past 14.75 and roots inside 14.65 (within
    // the sampling bound), and the fin faces, rim and bore within 0.005.
    let f = bench("33-large-stl.ink");
    let dir = Scratch::new("b33");
    let knob = std::process::Command::new("python3").arg(bench("tools/33-gen-knob.py")).arg(dir.path("generated.stl")).arg("240").output().expect("python3");
    assert!(knob.status.success(), "tools/33-gen-knob.py: {}", String::from_utf8_lossy(&knob.stderr));
    assert!(
        std::fs::read(dir.path("generated.stl")).unwrap() == std::fs::read(bench("inputs/33-knob.stl")).unwrap(),
        "inputs/33-knob.stl is not what tools/33-gen-knob.py writes at n 240"
    );
    let measured = num(&ink(["check", &f])["samples"][0]["measured_error"]["max"]);
    assert!(0.0 < measured && measured <= 0.2165, "measured error {measured}");
    let v = &ink(["mass", &f, "--target", "knob", "--res", "0.3"])["targets"][0]["volume"];
    near("volume", num(&v["estimate"]), 43663.15, 0.002);
    let m = dir.ink(["mesh", &f, "--target", "knob", "--res", "0.25", "--out", &dir.path("knob.stl")]).ok();
    one_watertight_shell(&m);
    near("mesh volume", num(&m["objects"][0]["mesh_volume"]), 43666.55, 0.001);
    let mut rays: Vec<[[f64; 3]; 2]> = (0..240)
        .map(|k| {
            let a = 2.0 * std::f64::consts::PI * k as f64 / 240.0;
            [[20.0 * a.cos(), 20.0 * a.sin(), 10.0], [-a.cos(), -a.sin(), 0.0]]
        })
        .collect();
    let (c, s21) = (21f64.to_radians().cos(), 21f64.to_radians().sin());
    let radial = |p: &Value| num(&p[0]).hypot(num(&p[1]));
    let z = |p: &Value| num(&p[2]);
    let extra: [(&str, [[f64; 3]; 2], &dyn Fn(&Value) -> f64, f64); 5] = [
        ("fin top z 59", [[15.0, 0.0, 70.0], [0.0, 0.0, -1.0]], &z, 59.0),
        ("fin bottom z 24 (from the gap)", [[15.0, 0.0, 23.5], [0.0, 0.0, 1.0]], &z, 24.0),
        ("fin top z 23 (from the gap)", [[15.0, 0.0, 23.5], [0.0, 0.0, -1.0]], &z, 23.0),
        ("fin rim r 22", [[30.0 * c, 30.0 * s21, 22.5], [-c, -s21, 0.0]], &radial, 22.0),
        ("bore r 4", [[0.0, 0.0, 30.0], [c, s21, 0.0]], &radial, 4.0),
    ];
    rays.extend(extra.iter().map(|e| e.1));
    let r = records(&ink(["ray", &f, "--target", "knob", "--rays", &serde_json::to_string(&rays).unwrap()])["rays"]);
    let hit = |i: usize| {
        let p = &r[i]["nearest"]["point"];
        assert!(p.is_array(), "ray {i} misses:\n{}", pretty(&r[i]));
        p.clone()
    };
    let tips = (0..240).step_by(2).map(|i| radial(&hit(i))).fold(f64::INFINITY, f64::min);
    let roots = (1..240).step_by(2).map(|i| radial(&hit(i))).fold(0.0, f64::max);
    assert!(tips >= 14.75 && roots <= 14.65, "knurl tips from r {tips}, roots to r {roots}");
    for (k, (name, _, read, want)) in extra.iter().enumerate() {
        let got = read(&hit(240 + k));
        assert!((got - want).abs() <= 0.005, "{name}: {got}");
    }
}
