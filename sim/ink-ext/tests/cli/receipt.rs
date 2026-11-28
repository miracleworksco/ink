//! `ink solve`'s setup receipt: where each fix and load landed, what each
//! holds, the material's sources and the bodies left out, so a setup that
//! is wrong reads wrong at a glance.

use crate::common::*;
use serde_json::Value;

/// A number, NaN where it is none (so every comparison with it fails).
fn n(v: &Value) -> f64 {
    v.as_f64().unwrap_or(f64::NAN)
}

/// A bar along x from 0 to 40, 4 x 4, clamped and pressed down 10 N: as
/// meant (clamp at x = 0, load on the tip), with the load put halfway, with
/// clamp and load swapped end for end, and pressed on its tip by a pressure.
const BAR: &str = r#"(import "std")
(body bar (shape (move (box 40 4 4) 20 0 0)))
(study meant (shape (shape-of bar)) (elastic (E 70000MPa) (nu 0.3))
  (fix base (- x 1)) (load tip (force 0 0 -10N) (region (- 39 x))) (res 1))
(study halfway (shape (shape-of bar)) (elastic (E 70000MPa) (nu 0.3))
  (fix base (- x 1)) (load tip (force 0 0 -10N) (region (- (abs (- x 20)) 0.5))) (res 1))
(study swapped (shape (shape-of bar)) (elastic (E 70000MPa) (nu 0.3))
  (fix base (- 39 x)) (load tip (force 0 0 -10N) (region (- x 1))) (res 1))
(study pressed (shape (shape-of bar)) (elastic (E 70000MPa) (nu 0.3))
  (fix base (- x 1)) (load end (pressure 1MPa) (region (- 39.99 x))) (res 1))
"#;

#[test]
fn a_misplaced_load_or_a_clamp_on_the_wrong_end_reads_so_in_the_receipt() {
    let sc = Scratch::new("receipt-bar");
    let f = sc.model("bar.ink", BAR);
    let solve = |study: &str| sc.ink(["solve", &f, "--study", study]).ok();
    let (meant, halfway, swapped) = (solve("meant"), solve("halfway"), solve("swapped"));
    let fix = |d: &Value| (d["reactions"][0].clone(), d["loads"][0].clone());
    let at = |v: &Value| n(&v["landed"]["at"][0]);
    // As meant: the clamp holds the end face and the sides' first mm (16 + 16 mm2, 16 mm3; its
    // box to within half a cell at the rim), all of the reaction, along x y z; the load is the
    // last mm's material, at x = 39.5, its moment r x F about the origin, and about the
    // part's centroid (20, 0, 0).
    let (base, tip) = fix(&meant);
    assert!(
        (n(&base["landed"]["area_mm2"]) - 32.0).abs() < 0.5
            && (n(&base["landed"]["volume_mm3"]) - 16.0).abs() < 0.1
            && n(&base["landed"]["lo"][0]).abs() < 1e-6
            && (0.5..=1.0 + 1e-6).contains(&n(&base["landed"]["hi"][0]))
            && (at(&base) - 0.25).abs() < 0.05
            && (n(&base["share"]) - 1.0).abs() < 1e-9
            && base["dofs"] == serde_json::json!(["x", "y", "z"]),
        "{}",
        pretty(&base)
    );
    assert!(
        tip["kind"] == "force"
            && (n(&tip["landed"]["volume_mm3"]) - 16.0).abs() < 0.1
            && (at(&tip) - 39.5).abs() < 1e-6
            && (n(&tip["landed"]["lo"][0]) - 39.0).abs() < 1e-6
            && (n(&tip["moment_nmm"][1]) - 395.0).abs() < 1e-6
            && (n(&tip["moment_about_part_centroid_nmm"][1]) - 195.0).abs() < 1e-6,
        "{}",
        pretty(&tip)
    );
    // Put halfway: the same force, but at x = 20, its moment about the origin half.
    let (_, mid) = fix(&halfway);
    assert!((at(&mid) - 20.0).abs() < 1e-6 && (n(&mid["moment_nmm"][1]) - 200.0).abs() < 1e-6, "{}", pretty(&mid));
    // Swapped: the clamp at the far end, the load at the near one.
    let (far, near) = fix(&swapped);
    assert!(
        at(&far) > 39.5 && (39.0 - 1e-6..39.5).contains(&n(&far["landed"]["lo"][0])) && at(&near) < 1.0 && (n(&near["moment_nmm"][1]) - 5.0).abs() < 1e-6,
        "{}\n{}",
        pretty(&far),
        pretty(&near)
    );
    // In words, as plainly.
    let words = sc.ink(["solve", &f, "--study", "swapped", "--plain"]).text();
    for said in ["base: holds x y z", "centred at (39.7", "tip: force (0, 0, -10) N", "acting at (0.5,", "100% of all", "bodies: solved bar"] {
        assert!(words.contains(said), "no {said:?} in:\n{words}");
    }
}

#[test]
fn a_pressures_receipt_gives_its_face_and_its_range() {
    let sc = Scratch::new("receipt-pressure");
    let f = sc.model("bar.ink", BAR);
    let d = sc.ink(["solve", &f, "--study", "pressed"]).ok();
    let end = &d["loads"][0];
    let (l, p) = (&end["landed"], &end["landed"]["pressure_mpa"]);
    assert!(
        end["kind"] == "pressure"
            && (16.0..16.5).contains(&n(&l["area_mm2"]))
            && n(&l["at"][0]) > 39.9
            && [&p["min"], &p["mean"], &p["max"]].iter().all(|v| (n(v) - 1.0).abs() < 1e-9)
            && (n(&end["total"][0]) + 16.0).abs() < 0.1,
        "{}",
        pretty(end)
    );
}

#[test]
fn the_receipt_names_the_bodies_left_out_and_where_the_material_came_from() {
    let sc = Scratch::new("receipt-bodies");
    let f = sc.model(
        "two.ink",
        r#"(import "std")
(body block (shape (move (box 10 4 4) 5 0 0)))
(body jig (shape (move (box 4 4 4) 20 0 0)))
(assembly root (instance block block) (instance jig jig))
(study one (shape (shape-of root/block)) (elastic (material petg) (E 2000MPa))
  (fix base (- x 1)) (load tip (force 0 0 -1N) (region (- 9 x))) (res 1))
(study both (shape (shape-of root)) (elastic (material petg))
  (fix base (- x 1)) (load tip (force 0 0 -1N) (region (- (abs (- x 9.5)) 0.5))) (res 1))
"#,
    );
    let one = sc.ink(["solve", &f, "--study", "one"]).ok();
    let left = &one["bodies"]["left_out"];
    assert!(
        one["bodies"]["solved"].as_array().is_some_and(|s| s.len() == 1 && s[0].as_str().is_some_and(|b| b.contains("block")))
            && left[0]["body"].as_str().is_some_and(|b| b.contains("jig"))
            && left[0]["why"] == "outside the study's shape"
            && one["material"]["from"]["E"] == "(E 2000MPa)"
            && one["material"]["from"]["nu"] == "(material petg)"
            && n(&one["material"]["E_mpa"]) == 2000.0,
        "{}",
        pretty(&serde_json::json!({"bodies": one["bodies"], "material": one["material"]}))
    );
    // The jig in the study's shape but joined to no fix: dropped, and said so.
    let both = sc.ink(["solve", &f, "--study", "both"]).ok();
    let left = &both["bodies"]["left_out"];
    assert!(
        left[0]["body"].as_str().is_some_and(|b| b.contains("jig")) && left[0]["why"].as_str().is_some_and(|w| w.contains("dropped as floating")),
        "{}",
        pretty(&both["bodies"])
    );
}
