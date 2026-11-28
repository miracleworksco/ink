//! How drives share effort, and when they cannot supply it.

use super::{BAR, asserts_at, claim_value, claims, claims_at_moment, model};
use ink_lang::claims::Status;

/// Two actuators on one shaft, limits 10 and 1000 N mm: the split weights
/// each by its span of limits (least sum of (effort / span)^2), so the
/// small one is not given half, and the shares do not depend on units.
#[test]
fn actuators_share_by_their_limits() {
    let src = format!(
        r#"(import "std") (gravity 0 -9.81m/s2 0) {BAR}
        (body disc (shape (cylinder-z 5 5)) (attr density 1g/cm3))
        (assembly (shaft (angle 0))
          (instance bar bar (at) (turn z angle) (effort-limits -1000N*mm 1000N*mm))
          (instance d disc (at (move 0 0 -20)) (turn z angle) (effort-limits -10N*mm 10N*mm)))
        (assembly root (instance s (shaft 0)))"#
    );
    let mut m = model(&src);
    m.set_mass("root/s/bar", 162.0, [150.0, 0.0, 0.0], [[0.0; 3]; 3]);
    m.set_mass("root/s/d", 0.0, [0.0; 3], [[0.0; 3]; 3]);
    let e = &m.at(0.0).effort;
    let need = 162.0 * 9.81e-3 * 150.0;
    assert!((e[0] + e[1] - need).abs() < 1e-9 * need && (e[1] / e[0] - 1e-4).abs() < 1e-12, "{e:?}");
}

#[test]
fn a_held_motor_keeps_its_load_beside_an_independent_press() {
    let src = r#"(import "std")
        (gravity 0 -10m/s2 0)
        (body arm (shape (move (box 20 40 40) 10 0 0)) (attr mass 100g))
        (body head (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (hinge (a 0)) (instance link arm (at) (turn z a) (effort-limits -2N*mm 2N*mm)))
        (assembly (press (s 0)) (instance head head (at) (slide y s) (effort-limits -40000N 40000N)))
        (assembly root (instance a (hinge 0)) (instance p (press 0) (at (move 200 0 0))))"#;
    for rating in ["30000", "36000", "40000", "1e100"] {
        let cs = claims_at_moment(&src.replace("40000", rating), 0.0);
        let arm = cs.iter().find(|c| c.0 == "effort-limits:root/a/link").unwrap();
        assert!(arm.1 == Status::Fails && (arm.3 + 8.0).abs() < 1e-5, "{rating}: {arm:?}");
    }
}

#[test]
fn coupled_slide_and_turn_balance_each_coordinate_regardless_of_rating() {
    let src = r#"(import "std")
(gravity 0 -10m/s2 0)
(body arm (shape (move (box 20 40 40) -10 0 0)) (attr mass 100g))
(body head (shape (box 800 800 800)) (attr mass 4000100g))
(assembly (machine (s 0) (a 0))
  (instance head head (at) (slide y s) (effort-limits -40000N 40000N))
  (instance link arm (at (move 200 0 0)) (turn z (+ a (* 60 s))) (effort-limits -2N*mm 2N*mm)))
(assembly root (instance m (machine 1 -60)))"#;
    let cs = claims_at_moment(src, 0.0);
    for (name, margin) in [("head", -1.0), ("link", -8.0)] {
        let c = cs.iter().find(|c| c.0.ends_with(name)).unwrap();
        assert!(c.1 == Status::Fails && (c.3 - margin).abs() < 1e-5, "{c:?}");
    }
    // Reparameterise the actual drive matrix: s is now metres and a radians.
    // The slide margin is still in N; the same turn margin is -0.008 N m.
    let metric = r#"(import "std")
        (gravity 0 -10m/s2 0)
        (body arm (shape (move (box .020m .040m .040m) -.010m 0 0)) (attr mass .1kg))
        (body head (shape (box .800m .800m .800m)) (attr mass 4000.1kg))
        (assembly (machine (s 0) (a 0))
          (instance head head (at) (slide y (* 1m s)) (effort-limits -40000N 40000N))
          (instance link arm (at (move .200m 0 0)) (turn z (rad->deg (+ a (* (/ pi 3) 1000 s)))) (effort-limits -.002N*m .002N*m)))
        (assembly root (instance m (machine .001 (/ (neg pi) 3))))"#;
    let si = claims_at_moment(metric, 0.0);
    assert_eq!(cs.len(), 2);
    assert_eq!(si.len(), cs.len());
    for (name, status, _, margin) in &cs {
        let (_, si_status, _, si_margin) = si.iter().find(|c| &c.0 == name).unwrap();
        assert_eq!(status, si_status);
        if name.ends_with("/link") {
            assert!((si_margin / 1000.0 + 0.008).abs() < 1e-8, "{si:?}");
        }
        assert!((si_margin - margin).abs() < 1e-7, "{si:?}");
    }
    // Generalized loads, in N mm per *written* unit, change by the
    // coordinate conversion; joint efforts above are physical loads.
    let (mm, metres) = (model(src).at(0.0), model(metric).at(0.0));
    for (i, factor) in [1000.0, 180.0 / std::f64::consts::PI].into_iter().enumerate() {
        assert!((metres.f[i] / factor - mm.f[i]).abs() < 1e-7, "{:?}, {:?}", mm.f, metres.f);
    }
    // One more coordinate, with its own gravitational load, but still
    // only two drive directions. Neither motor can claim a complete balance.
    let missing = src
        .replace("(s 0) (a 0)", "(s 0) (a 0) (b 0)")
        .replace("(+ a (* 60 s))", "(+ b a (* 60 s))")
        .replace("(instance head head", "(instance passive arm (at) (turn z b)) (instance head head")
        .replace("(machine 1 -60)", "(machine 1 -60 0)");
    let cs = claims_at_moment(&missing, 0.0);
    assert_eq!(cs.len(), 2);
    assert!(cs.iter().all(|c| c.1 == Status::Unknown && c.2.contains("no actuator supplies")), "{cs:?}");
}

#[test]
fn a_held_motor_with_fixed_limits_reports_the_required_effort() {
    let src = r#"(import "std")
        (gravity 0 0 -10m/s2)
        (body payload (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (lift (height 0))
          (instance cart payload (at) (slide z height) (effort-limits 0 0)))
        (assembly root (instance a (lift 50)))"#;
    for (limits, status, margin) in [("0 0", Status::Fails, -10.0), ("5N 5N", Status::Fails, -5.0)] {
        let cs = claims_at_moment(&src.replace("effort-limits 0 0", &format!("effort-limits {limits}")), 0.0);
        assert!(cs[0].1 == status && (cs[0].3 - margin).abs() < 1e-8, "{cs:?}");
    }
    let equal = src.replace("effort-limits 0 0", "effort-limits 10N 10N");
    // Measuring the block with its published density can leave a few
    // ulps of weight beyond exactly 10. That cannot establish a failure.
    let c = &claims_at_moment(&equal, 0.0)[0];
    assert!(c.1 == Status::Holds && c.3 == 0.0 || c.1 == Status::Unknown && c.2.contains("only within rounding"), "{c:?}");
    let applied = equal.replace("(gravity 0 0 -10m/s2)", "(gravity 0 0 0)") + "\n(load weight root/a/cart (force 0 0 -10N))";
    let c = &claims_at_moment(&applied, 0.0)[0];
    assert!(c.1 == Status::Holds && c.3 == 0.0, "{c:?}");
    // One representable step above and below the fixed rating, then a
    // discrepancy larger than the load's rounding, on each side.
    for (weight, expected) in [("10.000000000000002", Status::Holds), ("9.999999999999998", Status::Holds), ("10.000001", Status::Unknown), ("9.999999", Status::Unknown)] {
        let loaded = applied.replace("(force 0 0 -10N)", &format!("(force 0 0 -{weight}N)")) + "\n(require support (value (nth (reaction root/a/cart) 2) < 11N))";
        let cs = claims_at_moment(&loaded, 0.0);
        let c = cs.iter().find(|c| c.0 == "support").unwrap();
        assert_eq!(c.1, expected, "{weight}: {cs:?}");
        if expected == Status::Holds {
            assert!((c.3 - 10.0).abs() < 1e-12, "{weight}: {cs:?}");
        } else {
            assert!(c.2.contains("fixed rating"), "{weight}: {cs:?}");
        }
    }
    // Short of the rating by 5e-9 N, half the rounding of the drive's 10 N
    // load (1e-9 of it): whether it supplies the motion is undecided, and
    // the reaction is the drive's to that half of its own rounding, 10 N's
    // 1e-8 / 2 = 5e-9 either way. A claim whose margin is within it (10.000000005
    // against 10.00000001) is undecided; one past it holds.
    let short = applied.replace("(force 0 0 -10N)", "(force 0 0 -10.000000005N)")
        + "\n(require edge (value (nth (reaction root/a/cart) 2) < 10.00000001N))\n(require far (value (nth (reaction root/a/cart) 2) < 11N))";
    let cs = claims_at_moment(&short, 0.0);
    let of = |n: &str| cs.iter().find(|c| c.0 == n).unwrap_or_else(|| panic!("no claim {n}: {cs:?}"));
    assert!(of("edge").1 == Status::Unknown && of("edge").2.contains("band either way"), "{cs:?}");
    assert!(of("far").1 == Status::Holds && (of("far").3 - 10.00000001).abs() < 1e-12, "{cs:?}");
    // A stationary passive support can carry the same load while a motor
    // moves the payload horizontally.
    let passive = r#"(import "std") (gravity 0 0 -10m/s2)
        (body payload (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (carriage (x 0))
          (instance block payload (at) (slide x x) (effort-limits -1N 1N)))
        (assembly (stage (x 0) (z 0))
          (instance c (carriage x) (at) (slide z z)))
        (assembly root (instance a (stage 1 50)))"#;
    assert_eq!(claims(passive)[0].1, Status::Holds);
}

#[test]
fn a_perpendicular_load_cannot_hide_an_unsupplied_coordinate() {
    let src = r#"(import "std") (gravity 0 0 0)
(body block (shape (box 5 5 5)) (attr mass 1g))
(assembly (m (s 0) (a 0))
 (instance motor block (at) (slide x (+ s a)) (effort-limits -20N 20N))
 (instance passive block (at) (slide y a)))
(assembly root (instance m (m 0 0)))
(load normal root/m/motor (force -10N 0 1e16N))
(load side root/m/passive (force 0 -1N 1e16N))"#;
    let cs = claims_at_moment(src, 0.0);
    assert!(cs[0].1 == Status::Unknown && cs[0].2.contains("0.500000") && cs[0].2.contains("no actuator supplies"), "{cs:?}");
}

/// A swerve module parked: the steer motor holds the heading, which no
/// load turns (gravity about a vertical axis), and the drive motor,
/// geared to both the wheel and the heading, holds the wheel's weight.
/// The steer motor cancels the drive motor's pull on the heading; the
/// sum leaves a residual of their products' rounding alone, which is not
/// a load no actuator supplies (it was, before the residual counted it).
#[test]
fn rounding_in_the_drives_sum_is_not_an_unsupplied_load() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2)
(body block (shape (box 20 20 20)) (attr mass 37g))
(assembly (module (heading 0) (roll 0))
 (instance steer block (at) (turn z heading))
 (instance wheel block (at (move 0 0 30)) (turn x roll))
 (instance arm block (at (move 0 13.7 30)) (turn x roll))
 (instance drive block (at) (turn z (+ (* 7.03125 roll) (* 4.5 heading))) (effort-limits -7090N*mm 7090N*mm))
 (instance steer-motor block (at) (turn z (* 26.090909 heading)) (effort-limits -4110N*mm 4110N*mm)))
(assembly root (instance m (module 0 0)))"#;
    let cs: Vec<_> = claims_at_moment(src, 0.0).into_iter().filter(|c| c.0.starts_with("effort-limits:")).collect();
    assert!(cs.len() == 2 && cs.iter().all(|c| c.1 == Status::Holds), "{cs:?}");
}

#[test]
fn nearly_parallel_drives_cannot_pass_as_redundant() {
    let src = r#"(import "std") (gravity 0 0 0)
(body block (shape (box 50 50 50)) (attr mass 1000g))
(assembly (m (s 0) (a 0))
 (instance shared block (at) (slide x (+ s a)) (effort-limits -10N 10N))
 (instance extra block (at) (slide x (+ s (* 1.00000000000001 a))) (effort-limits -10N 10N))
 (instance passive block (at) (slide x a)))
(assembly root (instance m (m 0 0)))
(load load root/m/shared (force -1N 0 0))
(load small root/m/passive (force -0.000000000001N 0 0))"#;
    let with_reaction = format!("{src} (require support (value (nth (reaction root) 0) < 20N))");
    let cs = claims_at_moment(&with_reaction, 0.0);
    assert_eq!(cs.len(), 3);
    assert!(cs.iter().all(|c| c.1 == Status::Unknown && c.2.contains("too close")), "{cs:?}");
}

/// The idler translates on the cart. Its turn s + |100-s| is constant
/// throughout the travel, but its calculated partial sometimes differs
/// from zero by a rounding step. The two 100 g masses need at most
/// 0.2 kg x (10 + 0.05) m/s^2 = 2.01 N, each its published mass exactly
/// (read as the density it implies times the volume the dynamics
/// measures, the cylinder made it 2.009937).
#[test]
fn a_cancelling_idler_partial_resolves_the_cart_rating() {
    let src = r#"(import "std") (gravity 0 -10m/s2 0) (t-span 0 12s)
        (body cart (shape (box 25 25 25)) (attr mass 100g))
        (body wheel (shape (cylinder-z 15 20)) (attr mass 100g))
        (assembly (axis (s 0))
          (instance cart cart (at) (slide y s) (effort-limits -4N 4N))
          (instance idler wheel (frame-of cart) (turn z (rad->deg (+ s (sqrt (square (- 100 s))))))))
        (assembly root (instance a (axis (* 50 (sin (* 1rad/s t))))))"#;
    claim_value(src, "effort-limits:root/a/cart", Status::Holds, 1.99, 1e-6);
}

/// The small partial is the only one the side load acts through. Its
/// narrow peak falls between the starting samples. Scaling its resolution
/// by 1e-9 of the unit x-partial alone gives a false hold over the span.
/// Bounded over each step, the step holding the peak is halved down to it,
/// and a moment there fails (read at its moments alone, it was undecided).
#[test]
fn a_tiny_partial_with_a_real_load_peak_is_still_checked() {
    let src = r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
        (body block (shape (box 5 5 5)) (attr mass 1g))
        (assembly (stage (q 0))
          (instance motor block (at) (slide x q) (effort-limits -50N 50N))
          (instance load block (at (move q (* q 1e-20 (exp (neg (square (* 10000 (- (/ t 1s) 0.015625)))))) 0))))
        (assembly root (instance s (stage 0)))
        (load side root/s/load (force 0 1e22N 0))"#;
    let a = asserts_at(src, &[]).remove(0).check;
    let t = a.motion.as_ref().and_then(|m| m.motion.t).unwrap_or(f64::NAN);
    assert!(a.status == Status::Fails && (t - 0.015625).abs() < 1e-4, "{:?} at {t}", a.status);
    let cs = claims_at_moment(src, 0.015625);
    assert!(cs[0].1 == Status::Fails && (cs[0].3 + 50.0).abs() < 1e-10, "{cs:?}");
}

#[test]
fn rounded_proportional_drives_leave_the_balance_undecided() {
    let a = "(instance a block (at) (slide x (+ x (* 0.1 y))) (effort-limits -10N 10N))";
    let b = "(instance b block (at) (slide x (+ (* 3 x) (* 0.3 y))) (effort-limits -10N 10N))";
    for (first, second) in [(a, b), (b, a)] {
        let src = format!(
            r#"(import "std") (gravity 0 0 0)
            (body block (shape (box 10 10 10)) (attr mass 10g))
            (assembly (stage (x 0) (y 0)) {first} {second})
            (assembly root (instance s (stage 0 0)))
            (load pull root/s/a (force 1N 0 0))"#
        );
        let cs = claims_at_moment(&src, 0.0);
        assert_eq!(cs.len(), 2);
        assert!(cs.iter().all(|c| c.1 == Status::Unknown && c.2.contains("too close")), "{cs:?}");
    }
}

#[test]
fn a_missing_moving_mass_blocks_only_drives_coupled_to_it() {
    // Two independent drives: one carries a massed body, the other a massless body.
    // The massed drive's effort-limits claim still holds; only the massless drive is unknown.
    let independent = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body weighted (shape (box 10 10 10)) (attr density 2.7g/cm3))
        (body bare (shape (box 10 10 10)))
        (assembly (slide-weighted (s 0)) (instance m weighted (at) (slide z s) (effort-limits -10N 10N)))
        (assembly (slide-bare (s 0)) (instance m bare (at) (slide z s) (effort-limits -10N 10N)))
        (assembly root
          (instance a (slide-weighted (* 10 (/ t 1s))))
          (instance b (slide-bare (* 10 (/ t 1s)))))"#;
    let cs = claims_at_moment(independent, 0.5);
    let a = cs.iter().find(|c| c.0 == "effort-limits:root/a/m").unwrap();
    let b = cs.iter().find(|c| c.0 == "effort-limits:root/b/m").unwrap();
    assert_eq!(a.1, Status::Holds, "{cs:?}");
    assert!(b.1 == Status::Unknown && b.2.contains("these moving bodies have no mass"), "{cs:?}");

    // When a drive IS coupled through a belt to a massless driven pulley/shaft, the missing
    // mass propagates across the coupling and blocks that drive, while an unrelated drive still holds.
    let coupled = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body rotor (shape (cylinder-z 4 20)) (attr density 7.85g/cm3))
        (body disc (shape (cylinder-z 10 6)) (attr density 2.7g/cm3))
        (body bare-disc (shape (cylinder-z 10 6)))
        (body weight (shape (box 10 10 10)) (attr density 2.7g/cm3))
        (assembly (motor (a 0))
          (instance rotor rotor (at) (turn z a) (effort-limits -1000N*mm 1000N*mm))
          (instance pulley disc (at (move 0 0 20) (frame-of rotor))))
        (assembly (belted (a 0))
          (instance m (motor a))
          (instance driven bare-disc (at (move 50 0 20)) (turn z a))
          (relate m/pulley driven (belt-drive 10 10 50N 0.9)))
        (assembly (other (s 0))
          (instance w weight (at (move 0 100 0)) (slide z s) (effort-limits -10N 10N)))
        (assembly root
          (instance d (belted (* 360 (/ t 1s))))
          (instance o (other (* 10 (/ t 1s)))))"#;
    let cs = claims_at_moment(coupled, 0.5);
    let drive = cs.iter().find(|c| c.0 == "effort-limits:root/d/m/rotor").unwrap();
    let other = cs.iter().find(|c| c.0 == "effort-limits:root/o/w").unwrap();
    assert!(drive.1 == Status::Unknown && drive.2.contains("these moving bodies have no mass"), "{cs:?}");
    assert_eq!(other.1, Status::Holds, "{cs:?}");
}

#[test]
fn a_crossed_limit_on_one_drive_leaves_another_drives_effort_and_unrelated_reaction_decided() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body block (shape (box 10 10 10)) (attr density 2.7g/cm3))
        (assembly (bad (s 0))
          (instance m block (at) (slide z s)
            (effort-limits (neg (- 5N (* 1N (abs s)))) (- 5N (* 1N (abs s))))))
        (assembly (good (s 0))
          (instance m block (at (move 50 0 0)) (slide z s) (effort-limits -10N 10N)))
        (assembly root
          (instance a (bad (* 20 (/ t 1s))))
          (instance b (good (* 10 (/ t 1s)))))
        (require rx (value (nth (reaction root/b/m) 2) < 10N))"#;
    let cs = claims_at_moment(src, 0.5);
    let bad = cs.iter().find(|c| c.0 == "effort-limits:root/a/m").unwrap();
    let good = cs.iter().find(|c| c.0 == "effort-limits:root/b/m").unwrap();
    let rx = cs.iter().find(|c| c.0 == "rx").unwrap();
    assert!(bad.1 == Status::Unknown && bad.2.contains("cross"), "{cs:?}");
    assert_eq!(good.1, Status::Holds, "{cs:?}");
    assert_eq!(rx.1, Status::Holds, "{cs:?}");
}
