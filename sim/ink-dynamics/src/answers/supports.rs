//! Reaction shares and friction at supports.

use super::{asserts_at, claims_at_moment, every_claim_at, every_value, load_error, near, turned, values_at};
use ink_lang::claims::Status;

/// A shaft on two 608ZZ bearings 60 mm apart (at 10 and 70 mm), pulled 100
/// N across its axis 40 mm past the second, its weight W at 60: by statics
/// over their spacing the second carries (100 x 100 + 50 W) / 60, the
/// first the rest, the other way (each bearing's y is the world's -z).
#[test]
fn a_shaft_on_two_bearings_by_hand() {
    let src = format!(
        "{SHAFT}{}",
        r#"
        (load pull root/shaft (force 0 0 -100N) (at (frame-of root/shaft tip)))
        (def s (reaction root/shaft (by root/b1 root/b2)))
        (require weight (value (nth (reaction root/shaft) 1) < 0))
        (require b1 (value (nth (nth s 0) 1) > 0))
        (require b2 (value (nth (nth s 1) 1) < 0))
        (require b2-across (value (abs (nth (nth s 1) 0)) < 1N))"#
    );
    // Each holds over its masses' band; the values by hand are read raw.
    let c = claims_at_moment(&src, 0.0);
    assert!(c.iter().all(|x| x.1 == Status::Holds), "{c:?}");
    let c = claims_at_moment(&turned(&src), 0.0);
    let of = |n: &str| c.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {c:?}"));
    let w = -of("weight").3 - 100.0;
    assert!(w > 0.45 && w < 0.48, "the shaft's weight {w}");
    // The weight acts at the measured centre, 60 mm to a cell's rounding.
    let r2 = (100.0 * 100.0 + 50.0 * w) / 60.0;
    let centre = w * CENTRE / 60.0;
    near(-of("b2").3, r2, 1e-7 * r2 + centre, "the far bearing's force");
    near(of("b1").3, r2 - 100.0 - w, 1e-7 * r2 + centre, "the near bearing's force");
    near(of("b2-across").3 + 1.0, 1.0, 1e-12, "no force along the bearing's x");
}

const SHAFT: &str = r#"(import "std") (gravity 0 0 -9.81m/s2)
    (import "../tests/fixtures/dynamics/bearing.ink" as brg)
    (body shaft (shape (move (cylinder-z 4 120) 0 0 60)) (attr density 7.85g/cm3) (frame tip (at (move 0 0 110))))
    (def axis (at (rot-x -90)))
    (assembly root
      (instance b1 brg.bearing (at (move 0 0 10) axis))
      (instance b2 brg.bearing (at (move 0 0 70) axis))
      (instance b3 brg.bearing (at (move 0 0 40) axis))
      (instance shaft shaft axis (turn z 0)))"#;

/// How far the shaft's measured centre may lie from its middle: an
/// integral over cells, the grid's own rounding (its forces move by its
/// weight times this over the bearings' span).
const CENTRE: f64 = 0.01;

/// Two radial bearings under a 20 N pull along the shaft and 1e12 N across
/// it at the tip: statics leaves the 20 N to neither, however large the
/// load across (it once hid under 1e-8 of it, 1e4 N). By hand: the
/// bearings carry 1e12 x 100 / 60 and 1e12 x 160 / 60 N across, the
/// arithmetic's rounding about 1e-16 of those; 20 N is past it by 1e14,
/// so the split is refused. The same with the shaft turned 30 degrees about
/// the world's z, its load across in the world's x and y with the pull.
#[test]
fn a_radial_load_cannot_excuse_an_unsupported_axial_load() {
    let ask = |src: &str, force: &str| {
        let src = format!(
            "{src}\n(load pull root/shaft (force {force}) (at (frame-of root/shaft tip)))\n(require axial (value (abs (nth (nth (reaction root/shaft (by root/b1 root/b2)) 0) 2)) < 1N))"
        );
        claims_at_moment(&src, 0.0).remove(0)
    };
    let c = ask(SHAFT, "0 -20N -1e12N");
    assert!(c.1 == Status::Unknown && c.2.contains("cannot carry") && c.2.contains(", 20.000000) N"), "{c:?}");
    // Across: 1e12 (cos 30, sin 30, 0); along, the axis (-sin 30, cos 30, 0), -20 N of it.
    let (c30, s30) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
    let turned = SHAFT.replace("(rot-x -90)", "(rot-x -90) (rot-z 30)");
    let c = ask(&turned, &format!("{}N {}N 0", 1e12 * c30 + 20.0 * s30, 1e12 * s30 - 20.0 * c30));
    assert!(c.1 == Status::Unknown && c.2.contains("cannot carry") && c.2.contains("name the thrust bearing"), "{c:?}");
    // With no pull along the shaft, the split stands.
    assert_eq!(ask(SHAFT, "0 0 -1e12N").1, Status::Holds);
}

/// What the split refuses: an axial force no named part carries (two radial
/// bearings), shares statics leaves open (three on one shaft), a part whose
/// record rates nothing, a rating its record lacks.
#[test]
fn a_split_the_model_does_not_determine_is_refused() {
    let at = |load: &str, by: &str| {
        let src =
            format!("{SHAFT}\n(load pull root/shaft (force {load}) (at (frame-of root/shaft tip)))\n(require a (value (nth (nth (reaction root/shaft (by {by})) 0) 1) > -1e9N))");
        claims_at_moment(&src, 0.0).into_iter().find(|x| x.0 == "a").unwrap()
    };
    let axial = at("0 -20N -100N", "root/b1 root/b2");
    assert!(axial.1 == Status::Unknown && axial.2.contains("cannot carry all of the reaction") && axial.2.contains(", 20.000000) N"), "{axial:?}");
    let three = at("0 0 -100N", "root/b1 root/b2 root/b3");
    assert!(three.1 == Status::Unknown && three.2.contains("undetermined (7 unknowns, 5 of them independent"), "{three:?}");
    assert_eq!(at("0 0 -100N", "root/b1 root/b2").1, Status::Holds);
    let none = load_error(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (import "../tests/fixtures/dynamics/screw.ink" as sc)
        (body shaft (shape (cylinder-z 4 120)) (attr density 7.85g/cm3))
        (assembly root (instance b sc.part) (instance shaft shaft (at) (turn z 0)))
        (def s (reaction root/shaft (by root/b)))"#,
    );
    assert!(none.contains("its record, sc.spec, has no load-ratings"), "{none}");
    let lacking = load_error(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (def g (get (read-json "../tests/fixtures/dynamics/ratings.json") "guide-without-C0"))
        (require c (value (carriage-static-ratio (vec 0 0 1N 0 0 0) g) < 1))"#,
    );
    assert!(lacking.contains("no key 'C0-N'"), "{lacking}");
    // A base bolted on a moving carriage, the carriage named its support:
    // nothing along the rail carries the base's pull, and the path has no
    // joint whose drive would.
    let bolted = claims_at_moment(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (import "../tests/fixtures/dynamics/carriage.ink" as car)
        (body base (shape (move (box 60 60 10) 0 0 5)) (attr mass 500g))
        (assembly (axis (s 0)) (instance rail car.rail) (instance c car.carriage (frame-of rail) (slide x s) (effort-limits -100N 100N))
          (instance base base (at (move 0 0 13) (frame-of c))))
        (assembly root (instance a (axis (* 50 (sin (* 360deg/s t))))))
        (require a (value (nth (nth (reaction root/a/base (by root/a/c)) 0) 2) > -1e9N))"#,
        0.25,
    )
    .into_iter()
    .find(|x| x.0 == "a")
    .unwrap();
    assert!(bolted.1 == Status::Unknown && bolted.2.contains("the path has no joint") && !bolted.2.contains("thrust bearing"), "{bolted:?}");
}

/// An MGN12H on its rail, moved 150 + 50 sin(2 pi t) mm by its slide's drive,
/// a 1.5 kg block on it 100 mm across the rail: at 0.25 s (a = -1.974 m/s^2)
/// it carries the 1.554 kg's weight, the block's roll m g 100 and yaw -m a
/// 100 about its frame; its slide's drive the pull along the rail; its
/// static ratio is MR's, 1471.5 / 38220.
#[test]
fn a_rail_carriage_under_an_offset_load_by_hand() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (import "../tests/fixtures/dynamics/carriage.ink" as car)
        (body block (shape (box 50 50 50)) (attr mass 1500g))
        (assembly (axis (s 0))
          (instance rail car.rail)
          (instance c car.carriage (frame-of rail) (slide x s) (effort-limits -100N 100N))
          (instance block block (at (move 0 100 33) (frame-of c))))
        (assembly root (instance a (axis (+ 150 (* 50 (sin (* 360deg/s t)))))))
        (def w (nth (reaction root/a/c (by root/a/c)) 0))
        (require fx (value (abs (nth w 0)) < 1N))
        (require fz (value (nth w 2) > 0))
        (require mx (value (nth w 3) > 0))
        (require mz (value (nth w 5) > 0))
        (require static (value (carriage-static-ratio w car.spec) < 0.5))"#;
    // The readings as the model computes them (`turned`); a holding claim answers its band's worst edge.
    let read = claims_at_moment(&turned(src), 0.25);
    let of = |n: &str| read.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {read:?}")).3;
    let a = 50.0 * (2.0 * std::f64::consts::PI).powi(2) / 1000.0;
    near(of("fx"), 0.0, 1e-12, "along the rail: the drive's");
    near(of("fz"), 1.554 * 9.81, 1e-9 * 15.24, "the weight");
    near(of("mx"), 1.5 * 9.81 * 100.0, 1e-6 * 1471.5, "roll");
    near(of("mz"), 1.5 * a * 100.0, 1e-6 * 296.1, "yaw");
    near(of("static"), 1471.5 / 38220.0, 1e-6 * 0.0385, "the static ratio");
}

/// A belt's tension from what it passes, installed and locked (Gates'
/// Formulas 5 and 6): a 2 kg carriage it moves 100 sin(2 pi t) mm pulls m a
/// = 7.896 N at 0.25 s, half of it on each of two equal spans, so 53.948 N
/// on 50 N of pretension, 0.26575 of the DOLD HTD3 9 mm belt's 203 N working
/// tension; a flywheel (spin-inertia 1e7 g mm^2) it turns 45 sin(2 pi t)
/// degrees takes J alpha = 310.06 N mm, on a 10 mm pitch radius 31.006 N,
/// all of it with no pretension (its slack span slack). The flywheel's
/// torque is read from its own driven joint: a pulley's reaction where an
/// undeclared belt turns it sums the belt's pull (moves_with) and is
/// unknown.
#[test]
fn a_belt_drive_by_hand() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (def belt (get (read-json "../tests/fixtures/dynamics/ratings.json") "belt"))
        (body block (shape (box 50 50 50)) (attr mass 2000g))
        (body wheel (shape (cylinder-z 50 10)) (attr density 7.85g/cm3) (attr spin-inertia 1e7g*mm2))
        (assembly (axis (x 0)) (instance c block (at) (slide x x) (effort-limits -100N 100N)))
        (assembly (drive (a 0)) (instance f wheel (at (move 100 300 0)) (turn z (* 0.5 a)) (effort-limits -1e5N*mm 1e5N*mm)))
        (assembly root (instance a (axis (* 100 (sin (* 360deg/s t))))) (instance d (drive (* 90 (sin (* 360deg/s t))))))
        (require pulled (value (belt-tension-ratio (belt-tension (nth (reaction root/a/c) 0) 50N) belt) < 0.5))
        (require turned (value (belt-tension (/ (nth (reaction root/d/f) 5) 10mm) 0) < 100N))"#;
    // The readings as the model computes them (`turned`); a holding claim answers its band's worst edge.
    let read = claims_at_moment(&turned(src), 0.25);
    let of = |n: &str| read.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {read:?}")).3;
    let w2 = (2.0 * std::f64::consts::PI).powi(2);
    near(of("pulled"), (2.0 * 100.0 * w2 / 1000.0 / 2.0 + 50.0) / 203.0, 1e-9, "the pulled belt's tension ratio");
    near(of("turned"), 1e7 * 1e-6 * 45.0f64.to_radians() * w2 / 10.0, 1e-9 * 31.0, "the turning belt's tension");
}

/// Split shares carry the reaction's band: two bearings 40 mm apart under
/// ROTOR_FAST's 6000 rpm disc, 73.0 N of it turning, read round over its
/// span (its weight alone, none of it along their x). Each carries half: a
/// claim its force along x stays under 30 N fails (36.49 N at 0 s), one
/// under 45 holds over the band, its worst case 36.54 N (half of m e w^2,
/// within the band's own rounding).
#[test]
fn a_split_share_carries_the_reactions_band() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (import "../tests/fixtures/dynamics/bearing.ink" as brg)
        (body disc (shape (subtract (cylinder-z 50 10) (move (cylinder-z 5 20) 30 0 0))) (attr density 7.85g/cm3))
        (assembly (spinner (a 0)) (instance disc disc (at (rot-x 90)) (turn z a) (effort-limits -1e4N*mm 1e4N*mm))
          (instance b1 brg.bearing (at (move 0 0 -20) (rot-x 90))) (instance b2 brg.bearing (at (move 0 0 20) (rot-x 90))))
        (assembly root (instance s (spinner (* 360 100 (/ t 1s)))))
        (def w (nth (reaction root/s/disc (by root/s/b1 root/s/b2)) 0))
        (require under-30 (value (nth w 0) < 30N))
        (require under-45 (value (nth w 0) < 45N))"#;
    let c = every_claim_at(src, &[]);
    let of = |n: &str| c.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {c:?}"));
    assert_eq!(of("under-30").1, Status::Fails, "{:?}", of("under-30"));
    assert_eq!(of("under-45").1, Status::Holds, "{:?}", of("under-45"));
    let x = every_value(src, "under-45");
    // (36.54 from what departs, the rest from what is not known of the disc and its masses' proven
    // bounds: 0.37. It repeats about its own line but for what is cut from it.)
    assert!((x - 36.905).abs() < 0.05, "a banded share's hold: its worst case {x}");
}

/// A stepped claim holding gives the moment it came nearest to failing:
/// under an upper bound its greatest value, over a lower its least (a
/// bearing's load under its rating was given at its lightest). A pendulum's
/// pin force along z over its swing: `< L` gives the mirror of `(neg fz) >
/// -L`.
#[test]
fn a_holding_claim_gives_its_value_nearest_to_failing() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -1000N*mm 1000N*mm)))
        (assembly root (instance p (pendulum (* 60 (sin (* 180deg/s t))))))
        (def fz (nth (reaction root/p/bar (at)) 2))
        (require under (value fz < 1e9N))
        (require over (value (neg fz) > -1e9N))"#;
    let (under, over) = (every_value(src, "under"), every_value(src, "over"));
    assert!(under > 1.8 && (under + over).abs() < 1e-9 * under, "under {under}, over {over}");
}

/// A carriage driven by a screw: its nut, on the slide's axis 20 mm above
/// the carriage's frame, carries the thrust along its own line (the slide's
/// drive then none), so the guide carriage's pitch moment is the reaction's
/// less the thrust's 20 mm arm; named alone, the guide takes the thrust as
/// through its own origin.
#[test]
fn a_screws_nut_carries_a_slides_thrust_along_its_line() {
    let c = claims_at_moment(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (import "../tests/fixtures/dynamics/carriage.ink" as car)
        (def spec (get (read-json "../tests/fixtures/dynamics/ratings.json") "nut"))
        (body nut (shape (cylinder-z 14 40)) (attr density 7.85g/cm3) (frame load (at)))
        (body block (shape (box 50 50 50)) (attr mass 2000g))
        (assembly (axis (s 0))
          (instance rail car.rail)
          (instance c car.carriage (frame-of rail) (slide x s) (effort-limits -100N 100N))
          (instance block block (at (move 0 0 33) (frame-of c)))
          (instance n nut (at (rot-y 90) (move 0 0 20) (frame-of c))))
        (assembly root (instance a (axis (+ 150 (* 50 (sin (* 360deg/s t)))))))
        (def s (reaction root/a/c (by root/a/c root/a/n)))
        (def r (reaction root/a/c (frame-of root/a/c)))
        (require thrust (value (nth (nth s 1) 2) > -1e9N))
        (require pitch (value (nth (nth s 0) 4) > -1e9N*mm))
        (require alone (value (nth (nth (reaction root/a/c (by root/a/c)) 0) 4) > -1e9N*mm))
        (require fx (value (nth r 0) > -1e9N))
        (require my (value (nth r 4) > -1e9N*mm))"#,
        0.25,
    );
    let of = |n: &str| c.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {c:?}")).3;
    assert!(of("fx").abs() > 4.0, "the thrust {}", of("fx"));
    near(of("thrust"), of("fx"), 1e-9 * of("fx").abs(), "the nut's share: the thrust");
    near(of("pitch") + 20.0 * of("thrust"), of("my"), 1e-9 * of("my").abs(), "the guide's pitch moment and the thrust's arm");
    near(of("alone"), of("my"), 1e-9 * of("my").abs(), "the guide named alone");
}

/// Two UCF201 flange units, mounting faces 100 mm apart, flanges outward:
/// each carries at its insert's centre (its load frame, 15 mm in from its
/// face), so a 100 N pull 150 mm along the shaft puts 100 x 135 / 70 =
/// 192.857 N on the far one and 92.857 N on the near (50 and 150 at their
/// mounting faces).
#[test]
fn a_flange_unit_carries_at_its_insert() {
    let of = values_at(
        r#"(import "std") (gravity 0 0 0) (import "../tests/fixtures/dynamics/units.ink" as u)
        (body shaft (shape (move (cylinder-z 6 200) 0 0 75)) (attr density 7.85g/cm3) (frame tip (at (move 0 0 150))))
        (assembly root (instance u1 u.flange (at)) (instance u2 u.flange (at (rot-x 180) (move 0 0 100))) (instance shaft shaft (at) (turn z 0)))
        (load pull root/shaft (force 0 -100N 0) (at (frame-of root/shaft tip)))
        (def s (reaction root/shaft (by root/u1 root/u2)))
        (require u1 (value (len2 (nth (nth s 0) 0) (nth (nth s 0) 1)) > 0))
        (require u2 (value (len2 (nth (nth s 1) 0) (nth (nth s 1) 1)) > 0))"#,
        0.0,
    );
    near(of("u2"), 100.0 * 135.0 / 70.0, 1e-9 * 192.9, "the far unit");
    near(of("u1"), 100.0 * 135.0 / 70.0 - 100.0, 1e-9 * 192.9, "the near unit");
}

/// Two UCP201 pillow blocks 100 mm apart along a shaft 30.2 mm up (their
/// shaft along their y): statically determinate, each carrying at its
/// insert's centre with its load frame's z along the shaft (y the world's
/// -z), so a 100 N load at 150 and the shaft's weight W at 75 put 150 +
/// 0.75 W up on the far block and 50 - 0.25 W down on the near.
#[test]
fn a_pillow_block_pair_shares_a_shaft() {
    let of = values_at(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (import "../tests/fixtures/dynamics/units.ink" as u)
        (body shaft (shape (move (cylinder-z 6 200) 0 0 75)) (attr density 7.85g/cm3) (frame tip (at (move 0 0 150))))
        (assembly root (instance u1 u.pillow (at)) (instance u2 u.pillow (at (move 0 100 0))) (instance shaft shaft (at (rot-x -90) (move 0 0 30.2)) (turn z 0)))
        (load pull root/shaft (force 0 0 -100N) (at (frame-of root/shaft tip)))
        (def s (reaction root/shaft (by root/u1 root/u2)))
        (require whole (value (nth (reaction root/shaft) 1) < 0))
        (require u1 (value (nth (nth s 0) 1) > -1e9N))
        (require u2 (value (nth (nth s 1) 1) > -1e9N))
        (require u2-axial (value (abs (nth (nth s 1) 2)) < 1N))"#,
        0.0,
    );
    let w = -of("whole") - 100.0;
    assert!(w > 1.6 && w < 1.9, "the shaft's weight {w}");
    near(of("u2"), -(150.0 + 0.75 * w), 1e-6 * 151.0, "the far block (the weight at its measured centre)");
    near(of("u1"), 50.0 - 0.25 * w, 1e-6 * 151.0, "the near block");
    near(of("u2-axial") + 1.0, 1.0, 1e-12, "nothing along the shaft");
}

/// A plate sliding on one rod on two LM8UU, an anti-rotation roller (a
/// 6321K38 v-groove track roller, given an assembly round it) 30 mm to the
/// side: no part carries along the slide's motion, so the drive takes all
/// of fx and the roller nothing along x, and the bushings, 40 mm apart,
/// carry the reaction's yaw as a couple, mz / 40 each.
#[test]
fn a_roller_carries_nothing_along_its_slide() {
    let of = values_at(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (import "../tests/fixtures/dynamics/bushing.ink" as lm)
        (import "../tests/fixtures/dynamics/bearing.ink" as gr)
        (assembly lmb (attr mass 15g) (instance bushing lm.bushing))
        (assembly roller (instance r gr.bearing))
        (body plate (shape (move (box 60 40 30) 0 15 10)) (attr mass 1000g))
        (assembly (axis (s 0))
          (instance plate plate (at) (slide x s) (effort-limits -100N 100N))
          (instance b1 lmb (at (rot-y 90) (move -20 0 0) (frame-of plate)))
          (instance b2 lmb (at (rot-y 90) (move 20 0 0) (frame-of plate)))
          (instance r roller (at (rot-x -90) (move 0 30 0) (frame-of plate))))
        (assembly root (instance a (axis (* 50 (sin (* 360deg/s t))))))
        (def s (reaction root/a/plate (by root/a/b1 root/a/b2 root/a/r)))
        (def w (reaction root/a/plate))
        (require fx (value (nth w 0) > -1e9N))
        (require mz (value (nth w 5) > -1e9N*mm))
        (require b1y (value (nth (nth s 0) 1) > -1e9N))
        (require b2y (value (nth (nth s 1) 1) > -1e9N))
        (require rx (value (nth (nth s 2) 0) > -1e9N))"#,
        0.25,
    );
    assert!(of("fx").abs() > 1.0, "the slide's pull {}", of("fx"));
    near(of("rx"), 0.0, 1e-9, "the roller along the slide");
    near(of("b2y"), of("mz") / 40.0, 1e-9 * of("mz").abs(), "the far bushing's share of the yaw");
    near(of("b1y"), -of("mz") / 40.0, 1e-9 * of("mz").abs(), "the near bushing's");
}

/// `(S FRAME)` gives the share in FRAME, carried at S's load frame: b1 of
/// SHAFT under a 100 N load between its bearings (at 25 mm), asked in the
/// world's axes at its centre, carries +R1 = 75 + W / 6 up the world's z;
/// asked at its face-a (B / 2 along its bore), the same force, with its
/// moment about that point.
#[test]
fn a_share_asked_in_a_frame_is_carried_where_the_part_carries() {
    let of = values_at(
        &format!(
            "{SHAFT}{}",
            r#"(load mid root/shaft (force 0 0 -100N) (at (at (move 0 25 0))))
        (def own (nth (reaction root/shaft (by root/b1 root/b2)) 0))
        (def world (nth (reaction root/shaft (by (root/b1 (at (move 0 10 0))) root/b2)) 0))
        (def face (nth (reaction root/shaft (by (root/b1 (frame-of root/b1 face-a)) root/b2)) 0))
        (require whole (value (nth (reaction root/shaft) 1) < 0))
        (require own-y (value (nth own 1) > -1e9N))
        (require world-z (value (nth world 2) > -1e9N))
        (require face-y (value (nth face 1) > -1e9N))
        (require face-mx (value (nth face 3) > -1e9N*mm))"#
        ),
        0.0,
    );
    let w = -of("whole") - 100.0;
    near(of("world-z"), 75.0 + w / 6.0, 1e-7 * 75.0 + w * CENTRE / 60.0, "b1 in the world's axes");
    near(of("own-y"), -of("world-z"), 1e-9 * 75.0, "b1 in its own: its y is the world's -z");
    near(of("face-y"), of("own-y"), 1e-9 * 75.0, "at its face, the same force");
    near(of("face-mx"), 3.5 * of("own-y"), 1e-9 * 263.0, "and its moment about the face");
}

/// A 51100 between two 608ZZ takes a 20 N push along the shaft, along its
/// +z (housing washer toward shaft washer); loaded the other way it would
/// have to pull, which a single-direction thrust bearing cannot: refused,
/// never checked as a push.
#[test]
fn a_single_direction_thrust_bearing_is_never_pulled() {
    let at = |fy: f64| {
        let src = format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2)
            (import "../tests/fixtures/dynamics/bearing.ink" as brg) (import "../tests/fixtures/dynamics/thrust-bearing.ink" as thr)
            (body shaft (shape (move (cylinder-z 4 120) 0 0 60)) (attr density 7.85g/cm3) (frame tip (at (move 0 0 110))))
            (def axis (at (rot-x -90)))
            (assembly root (instance b1 brg.bearing (at (move 0 0 10) axis)) (instance b2 brg.bearing (at (move 0 0 70) axis))
              (instance th thr.bearing (at (move 0 0 40) axis)) (instance shaft shaft axis (turn z 0)))
            (load pull root/shaft (force 0 {fy}N -100N) (at (frame-of root/shaft tip)))
            (require thz (value (nth (nth (reaction root/shaft (by root/b1 root/b2 root/th)) 2) 2) > -1e9N))"#
        );
        claims_at_moment(&src, 0.0).into_iter().find(|x| x.0 == "thz").unwrap()
    };
    let push = at(-20.0);
    assert!(push.1 == Status::Holds && (push.3 - 20.0).abs() < 1e-9, "{push:?}");
    let pull = at(20.0);
    assert!(pull.1 == Status::Unknown && pull.2.contains("root/th would pull what it carries by 20.000000 N") && pull.2.contains("pushes only"), "{pull:?}");
}

/// Two LMF10LUU flanged bushings on one rod, flanges 60 mm apart, sleeves
/// outward: each carries at its centre, L/2 - H = 21.5 mm along its bore
/// from its flange's face, so at -51.5 and +51.5, and 100 N at 60 puts
/// 100 x 111.5 / 103 = 108.252 N on the far one, 8.252 on the near (50 and
/// 150 at the flanges' faces).
#[test]
fn a_flanged_bushing_carries_at_its_centre() {
    let of = values_at(
        r#"(import "std") (gravity 0 0 0) (import "../tests/fixtures/dynamics/bushing.ink" as lm)
        (assembly lmf (attr mass 60g) (instance bushing lm.flanged))
        (body plate (shape (box 200 20 20)) (attr mass 100g) (frame tip (at (move 60 0 0))))
        (assembly (axis (s 0)) (instance plate plate (at) (slide x s))
          (instance b1 lmf (at (rot-y -90) (move -30 0 0) (frame-of plate))) (instance b2 lmf (at (rot-y 90) (move 30 0 0) (frame-of plate))))
        (assembly root (instance a (axis 0)))
        (load pull root/a/plate (force 0 0 -100N) (at (frame-of root/a/plate tip)))
        (def s (reaction root/a/plate (by root/a/b1 root/a/b2)))
        (require b1 (value (len3 (nth (nth s 0) 0) (nth (nth s 0) 1) (nth (nth s 0) 2)) > 0))
        (require b2 (value (len3 (nth (nth s 1) 0) (nth (nth s 1) 1) (nth (nth s 1) 2)) > 0))"#,
        0.0,
    );
    near(of("b2"), 100.0 * 111.5 / 103.0, 1e-9 * 108.3, "the far bushing");
    near(of("b1"), 100.0 * 111.5 / 103.0 - 100.0, 1e-9 * 108.3, "the near bushing");
}

/// A carriage and two screws' nuts on one slide (a dual-screw axis): their
/// shares of the drive follow the screws' stiffness, so the split is
/// refused, whichever nut is named first (the first named took it all).
#[test]
fn two_nuts_on_one_slide_are_undetermined() {
    for order in ["root/a/n1 root/a/n2", "root/a/n2 root/a/n1"] {
        let c = claims_at_moment(
            &format!(
                r#"(import "std") (gravity 0 0 0)
                (import "../tests/fixtures/dynamics/carriage.ink" as car) (import "../tests/fixtures/dynamics/nut.ink" as nt)
                (assembly nut (attr mass 50g) (instance n nt.nut))
                (body plate (shape (box 80 80 10)) (attr mass 500g) (frame p (at (move 0 10 20))))
                (assembly (axis (s 0)) (instance plate plate (at) (slide x s) (effort-limits -1000N 1000N)) (instance c car.carriage (at (frame-of plate)))
                  (instance n1 nut (at (rot-y 90) (move 0 -30 20) (frame-of plate))) (instance n2 nut (at (rot-y 90) (move 0 30 20) (frame-of plate))))
                (assembly root (instance a (axis 0)))
                (load push root/a/plate (force -100N 0 0) (at (frame-of root/a/plate p)))
                (require n (value (nth (nth (reaction root/a/plate (by root/a/c {order})) 1) 2) > -1e9N))"#
            ),
            0.0,
        );
        let n = c.iter().find(|x| x.0 == "n").unwrap();
        assert!(n.1 == Status::Unknown && n.2.contains("undetermined (7 unknowns, 6 of them independent"), "{order}: {n:?}");
    }
    // A second identical carriage: the rigid-table rule
    // shares what the carriages carry, and leaves the nuts' shares as open
    // as statics does; with one nut, statics fixes its share and the
    // carriages share the rest.
    let table = |nuts: &str, by: &str| {
        let c = claims_at_moment(
            &format!(
                r#"(import "std") (gravity 0 0 0)
                (import "../tests/fixtures/dynamics/carriage.ink" as car) (import "../tests/fixtures/dynamics/nut.ink" as nt)
                (assembly nut (attr mass 50g) (instance n nt.nut))
                (body plate (shape (box 80 80 10)) (attr mass 500g) (frame p (at (move 0 10 20))))
                (assembly (axis (s 0)) (instance plate plate (at) (slide x s) (effort-limits -1000N 1000N))
                  (instance c car.carriage (at (move -20 0 0) (frame-of plate))) (instance c2 car.carriage (at (move 20 0 0) (frame-of plate))) {nuts})
                (assembly root (instance a (axis 0)))
                (load push root/a/plate (force -100N 0 0) (at (frame-of root/a/plate p)))
                (require n (value (nth (nth (reaction root/a/plate (by root/a/c root/a/c2 {by})) 2) 2) > -1e9N))"#
            ),
            0.0,
        );
        c.into_iter().find(|x| x.0 == "n").unwrap()
    };
    let two = table("(instance n1 nut (at (rot-y 90) (move 0 -30 20) (frame-of plate))) (instance n2 nut (at (rot-y 90) (move 0 30 20) (frame-of plate)))", "root/a/n1 root/a/n2");
    assert!(two.1 == Status::Unknown && two.2.contains("undetermined"), "two nuts beside two carriages: {two:?}");
    let one = table("(instance n1 nut (at (rot-y 90) (move 0 -30 20) (frame-of plate)))", "root/a/n1");
    assert!(one.1 == Status::Holds && (one.3.abs() - 100.0).abs() < 1e-6, "one nut beside two carriages takes the push: {one:?}");
}

/// Two angular-contact bearings (a 7000AC's), both facing one way, under
/// SHAFT's 100 N at 25 mm: each carries at its pressure centre, its load
/// frame, pr tan 25 deg = 9 x 0.46631 = 4.197 mm toward its inner ring's
/// high shoulder (NSK's effective load centre), so at 5.803 and 65.803, and
/// the far one takes (100 x 19.197 + W x 54.197) / 60. Two tapered roller
/// bearings (a 30202's) likewise, at NSK's a = 8.2 mm from the outer ring's
/// back face (+z, T/2 = 5.875): 2.325 mm toward their inner rings' back
/// faces.
#[test]
fn a_bearing_carries_at_its_pressure_centre() {
    let src = |e: f64| {
        let centre = -e;
        format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (import "../tests/fixtures/dynamics/bearing.ink" as brg)
            (body shaft (shape (move (cylinder-z 4 120) 0 0 60)) (attr density 7.85g/cm3) (frame q (at (move 0 0 25))))
            (def axis (at (rot-x -90)))
            (assembly root (instance b1 (brg.bearing {centre}) (at (move 0 0 10) axis)) (instance b2 (brg.bearing {centre}) (at (move 0 0 70) axis)) (instance shaft shaft axis (turn z 0)))
            (load mid root/shaft (force 0 0 -100N) (at (frame-of root/shaft q)))
            (def s (reaction root/shaft (by root/b1 root/b2)))
            (require whole (value (nth (reaction root/shaft) 1) < 0))
            (require b2 (value (nth (nth s 1) 1) > -1e9N))"#
        )
    };
    for (part, e) in [("7000AC", 9.0 * 25f64.to_radians().tan()), ("30202", 8.2 - 5.875)] {
        let of = values_at(&src(e), 0.0);
        let w = -of("whole") - 100.0;
        near(-of("b2"), (100.0 * (25.0 - 10.0 + e) + w * (60.0 - 10.0 + e)) / 60.0, 1e-7 * 32.4 + w * CENTRE / 60.0, &format!("the far {part} at its pressure centre"));
    }
}

/// A source-defined flange unit carries at its bearing's outer ring, not
/// the retaining ring before it: 4461N18's outer ring spans 1.8 to 6.8 mm.
#[test]
fn a_source_flange_unit_carries_at_its_outer_ring() {
    let of = values_at(
        r#"(import "std") (import "../tests/fixtures/dynamics/units.ink" as u) (assembly root (instance u u.mounted))
        (require z (value (nth (get (frame-of root/u load) "origin") 2) > -1e9))"#,
        0.0,
    );
    near(of("z"), 4.3, 1e-9, "its load frame along the shaft");
}

/// On a slide, a carriage whose rail is turned 30 deg to the motion would
/// move against what it carries: refused, saying so; a bushing set across
/// the slide leaves the reaction to none, and the hint says what each part
/// carries on a slide, not to name a thrust bearing.
#[test]
fn a_part_not_set_to_its_slide_is_refused() {
    let skew = claims_at_moment(
        r#"(import "std") (gravity 0 0 0) (import "../tests/fixtures/dynamics/carriage.ink" as car)
        (body block (shape (box 50 50 50)) (attr mass 1500g))
        (assembly (axis (s 0)) (instance plate block (at) (slide x s) (effort-limits -100N 100N)) (instance c car.carriage (at (rot-z 30) (frame-of plate))))
        (assembly root (instance a (axis 150)))
        (load push root/a/plate (force 0 0 -50N))
        (require fz (value (nth (nth (reaction root/a/plate (by root/a/c)) 0) 2) > -1e9N))"#,
        0.0,
    )
    .into_iter()
    .find(|x| x.0 == "fz")
    .unwrap();
    assert!(skew.1 == Status::Unknown && skew.2.contains("root/a/c is not set to the slide's motion: its x (its rail) is at 30.000 deg"), "{skew:?}");
    let across = claims_at_moment(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s) (import "../tests/fixtures/dynamics/bushing.ink" as lm)
        (assembly lmb (attr mass 15g) (instance bushing lm.bushing))
        (body plate (shape (move (box 60 40 30) 0 15 10)) (attr mass 1000g))
        (assembly (axis (s 0)) (instance plate plate (at) (slide x s) (effort-limits -100N 100N))
          (instance b1 lmb (at (rot-y 90) (move -20 0 0) (frame-of plate))) (instance b2 lmb (at (rot-x -90) (move 20 0 0) (frame-of plate))))
        (assembly root (instance a (axis (* 50 (sin (* 360deg/s t))))))
        (require b (value (nth (nth (reaction root/a/plate (by root/a/b1 root/a/b2)) 0) 0) > -1e9N))"#,
        0.25,
    )
    .into_iter()
    .find(|x| x.0 == "b")
    .unwrap();
    assert!(across.1 == Status::Unknown && across.2.contains("on a slide a bushing carries across its bore") && !across.2.contains("thrust bearing"), "{across:?}");
}

/// std's `(slip R)`, how near a support is to sliding: a cart accelerating at half of g along its
/// floor needs the floor to push it along with m a while m g presses it
/// down, so its slip is a / g = 0.5, whatever its mass (validation P19's
/// traction limit at mu = 1 is a slip of 1). Under a stated friction of 0.6
/// it holds, at 0.4 it slides; dropped faster than g it lifts off, its slip
/// infinite.
#[test]
fn a_support_slips_past_its_friction() {
    let src = |slide: &str, x: &str, mu: f64| {
        format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
            (body block (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly (cart (x 0)) (instance b block (at) (slide {slide} x) (effort-limits -1e3N 1e3N)))
            (assembly root (instance c (cart {x})))
            (require grips (value (slip (reaction root (at))) < {mu}))"#
        )
    };
    let half_g = "(* 0.5 4905 (/ t 1s) (/ t 1s))";
    let held = claims_at_moment(&src("x", half_g, 0.6), 0.5);
    let grips = held.iter().find(|c| c.0 == "grips").unwrap();
    // It holds over its band, read at the band's edge: the slip's two
    // parts' bands taken apart (dynamics-mass-unbanded) put it within 0.0012
    // of the 0.5 it is exactly.
    assert!(grips.1 == Status::Holds && grips.3 > 0.5 && grips.3 < 0.5012, "{grips:?}");
    let raw = claims_at_moment(&turned(&src("x", half_g, 0.6)), 0.5);
    assert!((raw[0].3 - 0.5).abs() < 1e-9, "{raw:?}");
    assert_eq!(claims_at_moment(&src("x", half_g, 0.4), 0.5)[0].1, Status::Fails);
    let lifted = claims_at_moment(&src("z", "(* -9810 (/ t 1s) (/ t 1s))", 0.6), 0.5);
    assert!(lifted[0].1 == Status::Fails && lifted[0].3 > 1e300, "{lifted:?}");
}

/// Two MGN12H on one rail under a plate share as a rigid table on identical blocks, as HIWIN's
/// and THK's catalogs size them (split.rs's test holds the formulas to
/// their digits): the plate's weight off the rail's line is a roll both
/// carry alike, the shares sum to the reaction (read on the values, which
/// carry no band: each share's band is its own, so an identity between
/// them is not settled as a claim), and the answer says what it rests on.
/// An MGN12C beside an MGN12H is another block, refused.
#[test]
fn carriages_of_one_record_share_as_a_rigid_table() {
    let src = |second: &str| {
        format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
            (import "../tests/fixtures/dynamics/carriage.ink" as car) (import "../tests/fixtures/dynamics/carriage-short.ink" as cc)
            (body plate (shape (move (box 120 40 6) 0 0 16)) (attr mass 300g))
            (body block (shape (box 50 50 50)) (attr mass 1500g))
            (assembly (axis (s 0))
              (instance rail car.rail)
              (instance top plate (frame-of rail) (slide x s) (effort-limits -100N 100N))
              (instance c1 car.carriage (at (move 40 0 0) (frame-of top)))
              (instance c2 {second} (at (move -40 0 0) (frame-of top)))
              (instance block block (at (move 20 30 39) (frame-of top))))
            (assembly root (instance a (axis 150)))
            (def s (reaction root/a/top (by root/a/c1 root/a/c2)))
            (def r (reaction root/a/top (at (move 150 0 0))))
            (require c1-roll (value (nth (nth s 0) 3) > 0N*mm))
            (require c2-roll (value (nth (nth s 1) 3) > 0N*mm))
            (require c1-down (value (nth (nth s 0) 2) > 0N))
            (require c2-down (value (nth (nth s 1) 2) > 0N))
            (require weight (value (nth r 2) > 0N))
            (require near-the-block (value (- (nth (nth s 0) 2) (nth (nth s 1) 2)) > 0))"#
        )
    };
    let mine = ["c1-roll", "c2-roll", "c1-down", "c2-down", "weight", "near-the-block"];
    let a: Vec<_> = asserts_at(&src("car.carriage"), &[]).into_iter().filter(|x| mine.contains(&x.name.as_str())).collect();
    assert_eq!(a.len(), 6);
    for x in &a {
        assert_eq!(x.check.status, Status::Holds, "{}: {:?}", x.name, x.check.evidence);
    }
    let near_block = a.iter().find(|x| x.name == "near-the-block").unwrap();
    assert!(near_block.check.inexact.notes.iter().any(|n| n.contains("rigid table on identical blocks") && n.contains("G99TE24")), "{:?}", near_block.check.inexact.notes);
    for t in [0.0, 0.5, 1.0] {
        let of = values_at(&src("car.carriage"), t);
        let (r1, r2, f1, f2, w) = (of("c1-roll"), of("c2-roll"), of("c1-down"), of("c2-down"), of("weight"));
        assert!(r1.abs() > 1.0 && (r1 - r2).abs() < 1e-6, "both carry the roll alike at {t}: {r1} {r2}");
        assert!((f1 + f2 - w).abs() < 1e-6, "the shares sum to the reaction at {t}: {f1} + {f2} against {w}");
    }
    let mixed: Vec<_> = every_claim_at(&src("cc.carriage"), &[]).into_iter().filter(|c| mine.contains(&c.0.as_str()) && c.0 != "weight").collect();
    assert!(mixed.len() == 5 && mixed.iter().all(|c| c.1 == Status::Unknown && c.2.contains("undetermined")), "{mixed:?}");
}
