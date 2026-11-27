//! The simulation benchmarks' key answers that are not claims of their
//! files (moved from the kernel's crates/ink-cli/tests/cli/benchmarks.rs):
//! probes of a solve against hand calculations, claims at other params,
//! meshes. A benchmark's own claims are checked by the model rule
//! (models.rs). One test per benchmark, named `bNN_...`; the slow ones
//! (models.rs's SLOW) run with the slow set (`./cpu-test.sh --slow`), and those
//! giving a known wrong answer only by name (`cargo test -- --ignored NAME`).
//! Each file's header gives the hand calculation.

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

/// Probe `i` of a solve (its probes a table): its field `key`.
fn probe(d: &Value, i: usize, key: &str) -> Value {
    records(&d["probes"])[i][key].clone()
}

/// A number, or a result's `value` where it is given as one.
fn value(v: &Value) -> f64 {
    num(if v.is_object() { &v["value"] } else { v })
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn b08_spring_p99_stress_is_wahls() {
    // The deflection vs G d^4 / 8 D^3 n is the file's claims; the 99th
    // percentile von Mises within 25% of Wahl's 200 MPa.
    let d = ink(["solve", &bench("08-helical-spring.ink")]);
    near("p99 von Mises", num(&d["results"]["von_mises_percentiles_mpa"]["p99"]), 200.0, 0.25);
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn b12_truss_chords_and_midspan_match_the_pin_jointed_calc() {
    // benchmarks/tools/truss.py: top chord -31.83 MPa, bottom chord +21.22 (10%), midspan
    // deflection -0.0907 mm (20%: welded joints are stiffer).
    let d = ink(["solve", &bench("12-truss-bridge.ink"), "--points", "[[150,0,60],[150,0,0],[180,0,0]]"]);
    near("top chord von Mises", num(&probe(&d, 0, "von_mises_mpa")), 31.83, 0.10);
    near("bottom chord von Mises", num(&probe(&d, 1, "von_mises_mpa")), 21.22, 0.10);
    near("midspan dz", num(&probe(&d, 2, "displacement_mm")[2]), -0.0907, 0.20);
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn b13_vessel_ring_grows_as_lame_says_by_sectors_and_by_pressure() {
    // Lame: radial growth 0.03579 mm at mid-wall (5%), von Mises 49.5 MPa
    // (10%), at 45 deg by both studies, and at 0 deg by the true pressure.
    let pts = "[[35.355,35.355,20],[49.9,0.9,20]]";
    let radial = |d: &Value, i: usize| {
        let u = probe(d, i, "displacement_mm");
        num(&u[0]).hypot(num(&u[1]))
    };
    for study in ["pressure", "pressure-true"] {
        let d = ink(["solve", &bench("13-pressure-vessel.ink"), "--study", study, "--points", pts]);
        near(&format!("{study}: u_r at 45 deg"), radial(&d, 0), 0.03579, 0.05);
        near(&format!("{study}: von Mises at 45 deg"), num(&probe(&d, 0, "von_mises_mpa")), 49.5, 0.10);
        if study == "pressure-true" {
            near("pressure-true: u_x at 0 deg", num(&probe(&d, 1, "displacement_mm")[0]), 0.03579, 0.05);
        }
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn b15_topology_bracket_settles_and_is_stiffer_than_a_web_of_its_volume() {
    // The layout settles; the solid web's tip moves -0.0928 mm, by hand (15%); the layout, verified, is stiffer
    // (less compliant) than that web; its volume is within 20% of 30% of the envelope, 86 400 mm3; it meshes
    // watertight, as one shell. (benchmarks/tools/top88_mbb.py is the layout's independent reference.)
    let f = bench("15-topology-bracket.ink");
    let layout = ink(["topology", &f, "--name", "layout"]);
    assert!(layout["results"]["settled"] == true, "{}", pretty(&layout["results"]));
    let web = ink(["solve", &f, "--study", "baseline", "--points", "[[119,0,0]]"]);
    near("the web's tip dz", num(&probe(&web, 0, "displacement_mm")[2]), -0.0928, 0.15);
    let (verified, baseline) = (num(&ink(["solve", &f, "--study", "verify"])["results"]["compliance_nmm"]), num(&web["results"]["compliance_nmm"]));
    assert!(verified < baseline, "the layout's compliance {verified}, the web's {baseline}");
    near("volume", num(&ink(["mass", &f, "--target", "bracket"])["targets"][0]["volume"]["estimate"]), 86400.0, 0.2);
    let w = Scratch::new("b15-mesh");
    let mesh = &w.ink(["mesh", &f, "--target", "bracket", "--out", &w.path("bracket.stl")]).ok()["objects"][0];
    assert!(mesh["watertight"] == true && mesh["shells"] == 1, "{}", pretty(mesh));
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn b16_heat_sink_convection_removes_the_chips_10_w() {
    // Heat balance: the convection's reaction is the 10 W in, within 0.5%.
    // (The temperature rise vs fin theory is a claim.)
    let d = ink(["solve", &bench("16-pin-fin-heatsink.ink")]);
    // A heat load's total is a power, one number (it was [W, 0, 0]).
    near("chip", num(&d["loads"][0]["total"]), 10.0, 1e-9);
    near("convection", num(&d["reactions"][0]["total"]), -10.0, 0.005);
}

#[test]
fn b17_block_passes_80_c_when_lumped_capacitance_says() {
    // T = 80 C at -tau ln(1 - 55 / 83.33) = 496.2 s, within 3%: where the
    // headline history crosses 80. (T(600 s) is a claim.)
    let d = ink(["solve", &bench("17-transient-heatup.ink")]);
    let h: Vec<(f64, f64)> = d["results"]["transient"]["history"].as_array().unwrap().iter().map(|p| (num(&p[0]), num(&p[1]))).collect();
    let cross = h.windows(2).find(|w| w[0].1 < 80.0 && 80.0 <= w[1].1).map(|w| w[0].0 + (80.0 - w[0].1) / (w[1].1 - w[0].1) * (w[1].0 - w[0].0));
    near("80 C crossing", cross.unwrap_or_else(|| panic!("never crosses 80 C: {h:?}")), 496.2, 0.03);
}

#[test]
fn b23_bonded_bar_stretches_in_series() {
    // Springs in series 0.09757 mm (3%). (Mass and COM are claims.)
    let f = bench("23-multi-material.ink");
    let pull = ink(["solve", &f, "--study", "pull", "--points", "[[100,0,0]]"]);
    near("bar tip", num(&probe(&pull, 0, "displacement_mm")[0]), 0.09757, 0.03);
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn b24_snap_fit_bends_by_the_hand_force_and_takes_it_to_press() {
    // Under P = 3 E I y / a^3 = 2.359 N (PETG's E, std's 1500 MPa: the
    // 3.302 N this once asked was E 2100's, and read as a 28% fault) the
    // flat moves 1.0 mm (5%); pressed down 1.0 mm by a fix's move the
    // finger takes 2.359 N (5%; 2.378 on the file's 0.2 mm cells), and that
    // study's compliance is unknown, saying to read the reaction. Both on
    // 0.2 mm cells: the press study's own 0.05 resolves its root fillet's
    // peak for the stress claims, 18 GB, and a deflection needs none of it.
    let f = bench("24-snap-fit.ink");
    let press = ink(["solve", &f, "--study", "press", "--res", "0.2", "--points", "[[15.75,0,2.8]]"]);
    near("deflection under P", num(&probe(&press, 0, "displacement_mm")[2]), -1.0, 0.05);
    let pressed = ink(["solve", &f, "--study", "pressed"]);
    let finger = pressed["reactions"].as_array().unwrap().iter().find(|r| r["name"] == "finger").unwrap();
    near("force to press 1 mm", -num(&finger["total"][2]), 2.359, 0.05);
    let c = &pressed["results"]["compliance_nmm"];
    assert!(c["value"].is_null() && c["unavailable"].as_str().is_some_and(|u| u.contains("reaction")), "{}", pretty(c));
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn b25_plate_torsion_is_found_by_its_rotation() {
    // Among the plate's six modes one moves mostly about x (torsion 1, 1491.6
    // Hz by St-Venant). (Bending 1 vs Euler-Bernoulli is a claim.)
    let d = ink(["solve", &bench("25-modal-plate.ink"), "--study", "ring"]);
    let dominant: Vec<String> = d["results"]["modes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|m| m["effective_mass_fraction"].as_object().unwrap().iter().max_by(|a, b| num(a.1).total_cmp(&num(b.1))).unwrap().0.clone())
        .collect();
    assert!(dominant.iter().any(|m| m == "rx"), "no torsion mode: {dominant:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn b28_shaft_bending_torsion_shoulders_and_deflection_match_the_hand_calc() {
    // The header's criteria. At (14, 0, 60): sigma = M r / I = 26.41,
    // tau = T r / J = 17.61, von Mises 40.34 (5% each); the coupling
    // stub's torsion at r = 9, z = 0: 57.3 (5%); the deflection at z = 100:
    // 0.0493 mm (10%); the left shoulder's peak (the part's) 150 by Peterson,
    // the right shoulder's (bending alone) 29.2, 25% each. The torque passes
    // z = 60 whole only if the bearings leave the shaft free to spin: both
    // are ball joints, so none of it leaks out through them.
    let mut pts = vec![[14.0, 0.0, 60.0], [9.0, 0.0, 0.0], [0.0, -9.0, 0.0], [0.0, 0.0, 100.0]];
    // The right shoulder's fillet (r 3, centre r 13, z 183), 0.1 mm into
    // the material on the tension side.
    pts.extend((0..=9).map(|k| {
        let a = (10.0 * k as f64).to_radians();
        [-(13.0 - 3.1 * a.cos()), 0.0, 183.0 - 3.1 * a.sin()]
    }));
    let d = ink(["solve", &bench("28-stepped-shaft.ink"), "--points", &format!("{pts:?}")]);
    let stress = |i: usize, c: &str| num(&probe(&d, i, "stress_mpa")[c]);
    near("bending sigma_zz", -stress(0, "zz"), 26.41, 0.05);
    near("torsion tau_yz", stress(0, "yz"), 17.61, 0.05);
    near("von Mises", num(&probe(&d, 0, "von_mises_mpa")), 40.34, 0.05);
    near("stub torsion (9, 0)", stress(1, "yz"), 57.3, 0.05);
    near("stub torsion (0, -9)", stress(2, "zx"), 57.3, 0.05);
    near("deflection at z = 100", -num(&probe(&d, 3, "displacement_mm")[0]), 0.0493, 0.10);
    let peak = &d["results"]["max_von_mises_mpa"];
    near("left shoulder peak", num(&peak["value"]), 150.0, 0.25);
    let at = num(&peak["at"][2]);
    assert!((15.0..25.0).contains(&at), "the peak is not at the left shoulder (z 20): z {at}");
    let right = (4..pts.len()).map(|i| num(&probe(&d, i, "von_mises_mpa"))).fold(0.0, f64::max);
    near("right shoulder peak", right, 29.2, 0.25);
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn b36_supports_match_beam_theory_and_mechanisms_are_named() {
    // (a) pin and roller: midspan 0.10077 mm (3%), reactions 50 / 50 N
    // (0.1%), no axial reaction, the ends turn 1.5e-3 rad each way (3%);
    // (b) self-aligning bearings: 0.07945 mm (3%), 500 / 500 N; (c) the stop
    // takes 1000 / 98.5 = 10.152 N (1e-4); free on its pin or a frictionless
    // seat, the lever is a mechanism, named by its turn and the moment.
    let f = bench("36-pin-supports.ink");
    let a = ink(["solve", &f, "--study", "simply", "--points", "[[100,0,0]]"]);
    near("beam midspan", -num(&probe(&a, 0, "displacement_mm")[2]), 0.10077, 0.03);
    near("pin reaction", num(&a["reactions"][0]["total"][2]), 50.0, 0.001);
    near("roller reaction", num(&a["reactions"][1]["total"][2]), 50.0, 0.001);
    assert!(num(&a["reactions"][0]["total"][0]).abs() < 1e-3, "axial reaction: {}", pretty(&a["reactions"]));
    near("pin turn", num(&a["free_motions"][0]["motions"][0]["amount"]), 1.5e-3, 0.03);
    near("roller turn", -num(&a["free_motions"][1]["motions"][0]["amount"]), 1.5e-3, 0.03);
    let b = ink(["solve", &f, "--study", "bearings", "--points", "[[0,0,100]]"]);
    near("shaft midspan", -num(&probe(&b, 0, "displacement_mm")[0]), 0.07945, 0.03);
    near("left bearing", num(&b["reactions"][0]["total"][0]), 500.0, 0.001);
    near("right bearing", num(&b["reactions"][1]["total"][0]), 500.0, 0.001);
    let c = ink(["solve", &f, "--study", "hinged-stop"]);
    near("stop reaction", -num(&c["reactions"][1]["total"][1]), 10.152284, 1e-4);
    for study in ["hinged", "frictionless"] {
        let e = ink_err(["solve", &f, "--study", study]);
        let msg = e["error"]["message"].as_str().unwrap_or("");
        assert!(msg.contains("turn about the axis along (0, 0, 1) through (0, 0, 0), and the loads drive it (1000 N mm"), "{study}: {msg}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn b21_gyroid_block_is_one_open_solid_of_sheet_stiffness() {
    // One solid, no cavities (the labyrinths stay open at the sides); E*/Es
    // = (100 N / 1296 mm2) / (dz / 32) / 1700 in [0.1, 0.35] (published
    // sheet gyroids); the redistanced lattice's volume fraction within 5% of
    // 3.0937 t / a = 0.3867 of its 36 x 36 x 32 region (an estimate: the
    // claim, measured on the samples, takes 25 times as long).
    let f = bench("21-gyroid-lattice.ink");
    let v = ink(["mass", &f, "--res", "0.3", "--target", "lattice-only"]);
    near("volume fraction", num(&v["targets"][0]["volume"]["estimate"]) / 41472.0, 0.3867, 0.05);
    let c = &ink(["components", &f, "--target", "block"])["targets"][0];
    assert!(c["solids"] == 1 && c["cavities"] == 0, "{}", pretty(c));
    // Its sheets are set flush on the plates: the 0 between them is no
    // surface, so the cells glue no gap there (576 cells read as glued when
    // a box across the seam had to be proven solid by bounds alone).
    let d = ink(["solve", &f, "--points", "[[0,0,17]]"]);
    assert!(!d["warnings"].as_array().unwrap().iter().any(|w| w.as_str().is_some_and(|w| w.contains("glued shut"))), "{}", pretty(&d["warnings"]));
    let dz = num(&probe(&d, 0, "displacement_mm")[2]);
    let ratio = (100.0 / 1296.0) / (-dz / 32.0) / 1700.0;
    assert!((0.1..=0.35).contains(&ratio), "E*/Es {ratio} (dz {dz})");
}

#[test]
fn b30_flexure_guide_is_as_stiff_and_stressed_as_beam_theory() {
    // F / uy within 15% of k = 24 E I / L^3 = 0.28125 N/mm; the surface
    // stress at (4, 9.5, 0) within 10% of 1.125 MPa; converged, in
    // equilibrium, no warnings (no large rotation), peaks below yield.
    // (The first mode is a claim.) The library's PETG, E 1500 MPa: held
    // against the hand values at 2100, the stress read 27% low.
    let d = ink(["solve", &bench("30-compliant-guide.ink"), "--study", "translate", "--points", "[[43,0,0],[4,9.5,0]]"]);
    near("k", 0.140625 / num(&probe(&d, 0, "displacement_mm")[1]), 0.28125, 0.15);
    near("sigma_xx", num(&probe(&d, 1, "stress_mpa")["xx"]).abs(), 1.125, 0.10);
    let peak = ["max_von_mises_mpa", "max_von_mises_near_supports_mpa"].iter().map(|k| value(&d["results"][*k])).fold(0.0, f64::max);
    assert!(
        d["solver"]["converged"] == true && num(&d["equilibrium_imbalance"]) < 1e-5 && d["warnings"].as_array().is_none_or(|w| w.is_empty()) && peak < 45.0,
        "solver {}, imbalance {}, warnings {}, peak {peak}",
        d["solver"],
        d["equilibrium_imbalance"],
        d["warnings"]
    );
}
