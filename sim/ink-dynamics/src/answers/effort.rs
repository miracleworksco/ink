//! Effort limits, mass bounds, and sampling over a motion.

use super::{BAR, asserts_at, claim_value, claims, claims_at, claims_at_moment, every_claim_at, load_error, model, near, rotor, witness, words};
use crate::Model;
use ink_lang::claims::Status;

/// The peak |effort| of actuator a over the span's measured moments.
fn peak(m: &Model, a: usize, span: (f64, f64)) -> f64 {
    let s = m.steps(span);
    let middles = s.ends.windows(2).map(|w| ink_core::eval::span_mid((w[0], w[1])));
    [span.0, span.1].into_iter().chain(middles).map(|t| m.at(t).effort[a].abs()).fold(0.0, f64::max)
}

/// A bar spinning 16 times a second under gravity in its plane: its
/// effort varies with its angle at constant speed, which the coordinate's
/// speed and acceleration give no steps for; 16 steps all land at one
/// phase. Its rotation's steps find the peak, m g r.
#[test]
fn steps_follow_a_turning_body() {
    let src = format!(
        r#"(import "std") (gravity 0 -9.81m/s2 0) (t-span 0 1s) {BAR}
        (assembly (spin (angle 0)) (instance bar bar (at) (turn z angle) (effort-limits -200N*mm 200N*mm)))
        (assembly root (instance s (spin (+ 135 (* 5760 (/ t 1s))))))"#
    );
    let mut m = model(&src);
    m.set_mass("root/s/bar", 162.0, [150.0, 0.0, 0.0], [[0.0; 3]; 3]);
    let true_peak = 162.0 * 9.81e-3 * 150.0;
    let p = peak(&m, 0, (0.0, 1.0));
    assert!(p <= true_peak * (1.0 + 1e-12) && p >= 0.99 * true_peak, "peak {p} of {true_peak}");
    let c = claims(&src);
    assert_eq!(c[0].1, Status::Fails, "{c:?}");
}

/// A cart pushed harder to the end: its peak is at the span's end, which
/// a step's middle never is.
#[test]
fn the_span_ends_are_measured() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body block (shape (box 100 100 100)) (attr density 1g/cm3))
        (assembly (cart (x 0)) (instance b block (at) (slide x x) (effort-limits -5.9N 5.9N)))
        (assembly root (instance c (cart (* 1000 (/ t 1s) (/ t 1s) (/ t 1s)))))"#;
    let mut m = model(src);
    m.set_mass("root/c/b", 1000.0, [0.0; 3], [[0.0; 3]; 3]);
    assert!((peak(&m, 0, (0.0, 1.0)) - 6.0).abs() < 1e-9);
    assert_eq!(claims(src)[0].1, Status::Fails);
}

/// A drive limited to 5 N for a moment (a 0.2 ms outage, and dips of 2e-10
/// and 2e-12 s), pushing 10 N all along: each of the limit's jumps is found
/// apart from the other to the resolution of t, so the moment inside the
/// outage is measured and the claim fails there.
#[test]
fn a_limits_own_jumps_are_steps() {
    for half in ["0.0001", "1e-10", "1e-12"] {
        let src = format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body cart (shape (box 100 100 100)) (attr mass 10kg))
            (assembly (belt (s 0)) (instance c cart (at) (slide x s) (effort-limits -50N (* 1N (select (- (abs (- (/ t 1s) 0.123456)) {half}) 5 50)))))
            (assembly root (instance b (belt (* 500 (/ t 1s) (/ t 1s)))))"#
        );
        let a = &asserts_at(&src, &[])[0];
        let t = a.check.motion.as_ref().and_then(|m| m.motion.t).unwrap_or(f64::NAN);
        assert!(a.check.status == Status::Fails && (t - 0.123456).abs() <= 1e-4, "dip of {half} either way: {:?} at {t}", a.check.status);
    }
}

/// A turntable spinning about its axis, gravity along it: no effort, so
/// what the efforts differ by is rounding, which splits no step and leaves
/// none loose.
#[test]
fn rounding_splits_no_step() {
    let src = format!(
        r#"(import "std") (gravity 0 -9.81m/s2 0) (t-span 0 1s) {BAR}
        (assembly (spin (angle 0)) (instance bar bar (at (rot-x 90) (rot-y 17)) (turn z angle) (effort-limits -1N*mm 1N*mm)))
        (assembly root (instance s (spin (* 360 (/ t 1s)))))"#
    );
    let s = model(&src).steps((0.0, 1.0));
    assert!(s.ends.len() < 100 && s.loose.is_empty(), "{} steps, {} loose", s.ends.len() - 1, s.loose.len());
    assert_eq!(claims(&src)[0].1, Status::Holds);
}

/// What an effort leaves out makes its claim unknown, saying so: a block
/// whose shape rides the screw's coordinate (its weight is not counted),
/// and a passive elbow moved as if driven (the motion cannot happen).
#[test]
fn a_claim_says_what_its_effort_leaves_out() {
    let lift = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body rotor (shape (cylinder-z 3 10)) (attr density 7.8g/cm3))
        (body (blk (h 0)) (shape (move (box 50 50 50) 100 0 (+ 25 h))) (attr density 7.85g/cm3))
        (assembly (lift (z 0)) (instance r rotor (at) (turn z (* 45 z)) (effort-limits -1N*mm 1N*mm)) (instance b (blk z)))
        (assembly root (instance l (lift (* 50 (/ t 1s) (/ t 1s)))))"#;
    let c = claims(lift);
    assert!(c[0].1 == Status::Unknown && c[0].2.contains("root/l/b"), "{c:?}");
    let arm = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body link (shape (move (box 200 10 10) 100 0 0)) (attr density 2.7g/cm3) (frame tip (at (move 200 0 0))))
        (assembly (arm (sh 0) (el 0)) (instance a link (at (rot-x 90)) (turn z sh) (effort-limits -2000N*mm 2000N*mm)) (instance b link (frame-of a tip) (turn z el)))
        (assembly root (instance r (arm (* 30 (/ t 1s)) (* 90 (/ t 1s)))))"#;
    let c = claims(arm);
    assert!(c[0].1 == Status::Unknown && c[0].2.contains("r.el"), "{c:?}");
}

/// A smooth band circulating round two fixed pulleys, its pulley turned by
/// the coordinate: its shape does not move, its material does, as its
/// `carried` says, so its inertia is counted: its mass at the pulley's
/// surface speed, m r^2.
#[test]
fn a_carried_belt_counts_its_material_moving() {
    let mut m = model(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (def r 10) (def c 100)
        (def P (path-loop (vec (path-arc (vec 0 0 0) (vec 0 0 1) r (vec 0 -1 0) pi) (path-line (vec 0 r 0) (vec (neg c) r 0))
                               (path-arc (vec (neg c) 0 0) (vec 0 0 1) r (vec 0 1 0) pi) (path-line (vec (neg c) (neg r) 0) (vec 0 (neg r) 0)))
                          (vec 0 0 1)))
        (body hub (shape (cylinder-z r 6)) (attr density 1g/cm3))
        (body (band (a 0)) (shape (path-sweep (max (- (abs x) 0.5) (- (abs y) 3)) P)) (attr density 1g/cm3)
          (attr carried (carried (max (- (abs y) 0.5) (- (abs z) 3)) P (path-length P) (* r (/ (* pi a) 180)))))
        (assembly (drive (a 0)) (instance p hub (at) (turn z a) (effort-limits -100N*mm 100N*mm)) (instance b (band a)))
        (assembly root (instance d (drive (* 360 (/ t 1s) (/ t 1s)))))"#,
    );
    m.set_mass("root/d/p", 0.0, [0.0; 3], [[0.0; 3]; 3]);
    assert!(m.deforming.is_empty() && m.carried.iter().map(|(l, _)| l.as_str()).eq(["root/d/b"]), "{:?} {:?}", m.deforming, m.carried);
    let (lam, _, len) = m.strip("root/d/b");
    // N mm per rad/s^2 of the pulley: m r^2, g mm^2 as 1e-6 N mm s^2.
    let (r, alpha) = (10.0, 720.0f64.to_radians());
    let want = 1e-6 * lam * len * r * r * alpha;
    let got = m.at(0.5).effort[0];
    assert!((got - want).abs() < 1e-9 * want, "{got} against {want}");
}

/// A length passed as a constant, read by the bar's shape and by its
/// drive (a ratio): a size, not a coordinate, so the bar is rigid and
/// weighs; the claim fails as the same bar of fixed length does.
#[test]
fn a_constant_a_shape_reads_is_a_size() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body (bar (l 300)) (shape (move (box 20 l 10) 0 (/ l 2) 0)) (attr density 2.7g/cm3))
        (assembly (pendulum (angle 0) (l 300)) (instance bar (bar l) (at (rot-x -90)) (turn z (* angle (/ l 300))) (effort-limits -50N*mm 50N*mm)))
        (assembly root (instance p (pendulum (* 60 (sin (* 180deg/s t))) 300)))"#;
    let m = model(src);
    assert!(m.deforming.is_empty() && m.coords.iter().map(|c| c.name.as_str()).eq(["p.angle"]), "{:?}", m.deforming);
    assert_eq!(claims(src)[0].1, Status::Fails);
}

/// A bar whose length is toleranced: an effort found at one length says
/// nothing of the others, so the claim over the band is unknown.
#[test]
fn an_effort_over_a_toleranced_param_is_unknown() {
    let c = claims(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (param L 300 (tol -30 30))
        (body bar (shape (move (box 20 L 10) 0 (/ L 2) 0)) (attr density 2.7g/cm3))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -130N*mm 130N*mm)))
        (assembly root (instance p (pendulum (* 60 (sin (* 180deg/s t))))))"#,
    );
    assert!(c[0].1 == Status::Unknown && words(&c[0].2).contains(&"L"), "{c:?}");
}

/// A massless planet that only t turns adds nothing to any effort: it
/// needs no mass, and is driven by nothing.
#[test]
fn a_body_only_t_moves_needs_no_mass() {
    let m = model(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -1000N*mm 1000N*mm)))
        (body planet (shape (cylinder-z 10 5)))
        (assembly root (instance p (pendulum (* 60 (sin (* 180deg/s t))))) (instance planet planet (at (move 50 0 0) (rot-z (* 72 (/ t 1s))))))"#,
    );
    assert_eq!(m.driven_by_nothing, ["root/planet"]);
}

/// A rotor at 6000 rpm for 2 s, its axis across gravity: round about its
/// spin, its effort reads the spin's rate and acceleration, not its angle,
/// so its steps need not resolve the rotation (4500 would, past the
/// budget), and its measured centre and inertia are made round, so its
/// effort is the same at every angle, to rounding.
#[test]
fn a_round_rotor_needs_no_steps_for_its_spin() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body rotor (shape (cylinder-z 10 30)) (attr density 7.8g/cm3))
        (assembly (motor (angle 0)) (instance r rotor (at (rot-x 90)) (turn z angle) (effort-limits -100N*mm 100N*mm)))
        (assembly root (instance m (motor (* 36000 (/ t 1s)))))"#;
    let s = model(src).steps((0.0, 2.0));
    assert!(s.ends.len() < 100 && s.loose.is_empty(), "{} steps, {} loose", s.ends.len() - 1, s.loose.len());
    assert_eq!(claims(src)[0].1, Status::Holds);
}

/// A stepper's rotor modelled as its D-cut shaft (turn symmetry of order
/// 1, its centre 0.1 mm off its axis), given its datasheet's 54 g cm2 as
/// its spin-inertia: a balanced rotor about its own z axis, read as round.
/// Spun up to 6000 rpm over 2 s, its axis across gravity, it needs no steps
/// for its spin, and its effort is J alpha throughout.
#[test]
fn a_rotor_with_a_spin_inertia_is_balanced_about_its_axis() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body rotor (shape (max (cylinder-z 2.5 46) (- y 2))) (attr density 7.85g/cm3) (attr spin-inertia 5400g*mm2))
        (assembly (motor (angle 0)) (instance r rotor (at (rot-x 90)) (turn z angle) (effort-limits -430N*mm 430N*mm)))
        (assembly root (instance m (motor (* 9000 (/ t 1s) (/ t 1s)))))"#;
    let m = model(src);
    let s = m.steps((0.0, 2.0));
    assert!(s.ends.len() < 100 && s.loose.is_empty(), "{} steps, {} loose", s.ends.len() - 1, s.loose.len());
    let j_alpha = 5400.0 * 18000f64.to_radians() * 1e-6;
    for t in [0.3, 1.1, 2.0] {
        let e = m.at(t).effort[0];
        assert!((e - j_alpha).abs() < 1e-9 * j_alpha, "{e} N mm at {t} s, not J alpha {j_alpha}");
    }
    assert_eq!(claims(src)[0].1, Status::Holds);
}

/// The same rotor 1 mm off its axis, a balanced one (a spin-inertia)
/// carried 1 mm off it, and a rotor 1 mm off its axis spun up on a
/// carriage, its axis along gravity (its centre circling shakes the
/// carriage, m r omega^2 along its slide): each one's departure from round
/// reaches past its limits (0.72 N mm of weight against 0.7; 29 N of
/// shaking against 20), so its claim is not decided with it as a margin:
/// its rotation is resolved, past the budget, and the claim fails at a
/// step it measures past its limits, or is unknown; it never holds.
#[test]
fn an_eccentric_past_its_limits_is_still_resolved() {
    let eccentric = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body rotor (shape (move (cylinder-z 10 30) 1 0 0)) (attr density 7.8g/cm3))
        (assembly (motor (angle 0)) (instance r rotor (at (rot-x 90)) (turn z angle) (effort-limits -0.7N*mm 0.7N*mm)))
        (assembly root (instance m (motor (* 36000 (/ t 1s)))))"#;
    let carried = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body rotor (shape (cylinder-z 10 30)) (attr density 7.8g/cm3) (attr spin-inertia 30000g*mm2))
        (assembly arm (instance r rotor (at (move 1 0 0))))
        (assembly (motor (angle 0)) (instance a arm (at (rot-x 90)) (turn z angle) (effort-limits -0.7N*mm 0.7N*mm)))
        (assembly root (instance m (motor (* 36000 (/ t 1s)))))"#;
    let shaken = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body cart (shape (box 100 60 10)) (attr density 2.7g/cm3) (frame axle (at (move 0 0 5))))
        (body rotor (shape (move (cylinder-z 10 30) 1 0 15)) (attr density 7.8g/cm3))
        (assembly (bench (x 0) (spin 0))
          (instance c cart (at) (slide x x) (effort-limits -20N 20N))
          (instance r rotor (frame-of c axle) (turn z spin) (effort-limits -100N*mm 100N*mm)))
        (assembly root (instance b (bench (* 5 (sin (* 180deg/s t))) (* 9000 (/ t 1s) (/ t 1s)))))"#;
    for src in [eccentric, carried, shaken] {
        assert!(!model(src).steps((0.0, 2.0)).loose.is_empty());
        let c = claims(src);
        assert!(c[0].1 != Status::Holds, "{c:?}");
    }
}

/// An FRC-style shooter shaft: a 1/2" hex shaft (round: a hexagon's turn
/// symmetry) turned by a joint, level, carrying on its frames a hub with
/// six bolt holes and a set screw's hole, a 22-tooth sprocket on a hex
/// bore (turn symmetry of order 2) bolted to the hub's face, and a clamp
/// collar with its screw on the collar's frame; its motor's limits +-limit
/// (N mm), turned by `angle` (degrees).
fn frc_shaft(limit: f64, angle: &str) -> String {
    format!(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
    (body shaft (shape (extrude-z (hexagon-2d 12.7) -75 75)) (attr density 2.7g/cm3)
      (frame hub (at (move 0 0 20))) (frame collar (at (move 0 0 -30))))
    (body hub
      (shape (max (min (extrude-z (circle-2d 25) 0 6) (extrude-z (circle-2d 12) 0 18))
                  (neg (extrude-z (hexagon-2d 12.8) -1 19))
                  (neg (repeat-polar (move (cylinder-z 2.1 20) 19 0 3) 6))
                  (neg (move (cylinder-y 1.6 30) 0 10 12))))
      (attr density 2.7g/cm3) (frame face (at (move 0 0 -3))))
    (body sprocket (shape (max (extrude-z (involute-gear-2d 22 1.5 20 0) -3 3) (neg (extrude-z (hexagon-2d 12.8) -4 4))))
      (attr density 7.85g/cm3))
    (body collar
      (shape (max (min (extrude-z (circle-2d 14) -4 4) (move (box 10 10 8) 0 13 0))
                  (neg (extrude-z (hexagon-2d 12.8) -5 5)) (neg (move (box 1 20 10) 0 14 0))))
      (attr density 2.7g/cm3) (frame screw (at (rot-y 90) (move -6 15 0))))
    (body screw (shape (cylinder-z 2 12)) (attr density 7.85g/cm3))
    (assembly (shooter (angle 0))
      (instance shaft shaft (at (rot-y 90)) (turn z angle) (effort-limits -{limit}N*mm {limit}N*mm))
      (instance hub hub (frame-of shaft hub))
      (instance sprocket sprocket (frame-of hub face))
      (instance collar collar (frame-of shaft collar))
      (instance screw screw (frame-of collar screw)))
    (assembly root (instance s (shooter {angle})))"#
    )
}

/// A spin-up to 6000 rpm with a 0.25 s time constant, and held.
const SPIN_UP: &str = "(* 36000 (- (/ t 1s) (* 0.25 (- 1 (exp (* -4 (/ t 1s)))))))";

/// The FRC shaft spun up to 6000 rpm and held: the hub, sprocket, collar
/// and screw ride the shaft's frames, so they spin with it as one mass.
/// Its weight off the axis (29 g mm) moves the effort by at most m r g,
/// 0.29 N mm, under 1% of the spin-up's 47 N mm: the steps read them with
/// the spin at rest (4096, loose, without) and say so. Held at 6000 rpm
/// from the start, that weight is all its effort, but far inside +-400 N
/// mm: taken as the claim's margin it decides, so the rotation is not
/// resolved, and the claim says the margin. So does the eccentric rotor
/// held at 6000 rpm, 0.72 N mm against +-100.
#[test]
fn a_hub_sprocket_and_collar_riding_a_hex_shaft_decide_at_6000_rpm() {
    let eccentric = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body rotor (shape (move (cylinder-z 10 30) 1 0 0)) (attr density 7.8g/cm3))
        (assembly (motor (angle 0)) (instance r rotor (at (rot-x 90)) (turn z angle) (effort-limits -100N*mm 100N*mm)))
        (assembly root (instance m (motor (* 36000 (/ t 1s)))))"#;
    for src in [frc_shaft(400.0, SPIN_UP), frc_shaft(400.0, "(* 36000 (/ t 1s))"), eccentric.to_string()] {
        let s = model(&src).steps((0.0, 2.0));
        assert!(s.ends.len() < 100 && s.loose.is_empty(), "{} steps, {} loose", s.ends.len() - 1, s.loose.len());
        let c = claims(&src);
        assert!(c[0].1 == Status::Holds, "{c:?}");
    }
    let held = frc_shaft(400.0, "(* 36000 (/ t 1s))");
    assert!(model(&held).steps((0.0, 2.0)).basis.contains("root/s/shaft with root/s/hub, root/s/sprocket, root/s/collar, root/s/screw"));
    let c = claims(&held);
    // (0.2884, from what departs: the hub, the collar and the screw repeat about their own lines, every
    // half turn or more, which leaves nothing of their measures' changes across them.)
    assert!(c[0].3.contains("within 0.2884"), "{}", c[0].3);
}

/// The FRC shaft spun up against limits of +-47 N mm: its effort read
/// round peaks at 46.84 at the start, 0.16 inside them, and its departure
/// from round can move it 0.29 at any angle. Started at 180 degrees its
/// effort there is 46.56: within its bound of a limit the claim is
/// undecided, never holding on the angle a step happens to sample (it held
/// so before). Started at 0, its effort there is 47.10, past a limit of
/// 46.7 at a moment the steps measure by more than its bodies' measures can
/// move it (0.35): a witness, so it fails (past 47 by 0.10, it is undecided).
#[test]
fn an_effort_within_its_bound_of_a_limit_is_undecided() {
    let c = claims(&frc_shaft(47.0, &format!("(+ 180 {SPIN_UP})")));
    assert!(c[0].1 == Status::Unknown && c[0].2.contains("undecided"), "{c:?}");
    let c = claims(&frc_shaft(46.7, SPIN_UP));
    assert!(c[0].1 == Status::Fails, "{c:?}");
}

/// Round bodies turning together spun up to 6000 rpm: a roller line (a
/// sub-assembly turned by one joint, tilted by a float turn) whose roller
/// sits along the axis and whose 20-tooth gear (turn symmetry of order 20)
/// is flipped onto it, and a gear riding a round rotor's shaft frame. None
/// needs steps for its spin; the effort is their inertia times the
/// spin-up, constant.
#[test]
fn round_bodies_turning_together_need_no_steps_for_their_spin() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body roller (shape (cylinder-z 20 40)) (attr density 1.2g/cm3))
        (body gear (shape (extrude-z (involute-gear-2d 20 1 20 0) -4 4)) (attr density 2.7g/cm3))
        (body rotor (shape (cylinder-z 4 40)) (attr density 7.85g/cm3) (frame shaft (at (move 0 0 20))))
        (assembly line (instance roller roller (at (move 0 0 30))) (instance gear gear (at (rot-x 180) (move 0 0 -10))))
        (assembly (intake (rollers 0))
          (instance line line (at (rot 0.57735 0.57735 0.57735 120) (move 10 20 30)) (turn z rollers) (effort-limits -100N*mm 100N*mm))
          (instance r rotor (at (move 200 0 0)) (turn z rollers) (effort-limits -100N*mm 100N*mm))
          (instance pinion gear (at (rot-z 90) (move 0 0 5) (frame-of r shaft))))
        (assembly root (instance i (intake (* 9000 (/ t 1s) (/ t 1s)))))"#;
    let s = model(src).steps((0.0, 2.0));
    assert!(s.ends.len() < 100 && s.loose.is_empty(), "{} steps, {} loose", s.ends.len() - 1, s.loose.len());
    assert!(claims(src).iter().all(|c| c.1 == Status::Holds), "{:?}", claims(src));
}

/// A follower on a sine track behind a carriage at constant speed, 16
/// periods over the span: no body turns and no coordinate accelerates,
/// and the steps' middles all fall at one phase of the track. Its
/// translation's acceleration and its partial velocity, kept within 1%
/// inside each step, find the peak.
#[test]
fn steps_follow_a_translation_and_its_partial_velocity() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body carriage (shape (box 50 50 10)) (attr density 2.7g/cm3))
        (body follower (shape (box 20 20 20)) (attr density 7.85g/cm3))
        (assembly (track (s 0))
          (instance c carriage (at) (slide x s) (effort-limits -0.5N 0.5N))
          (instance f follower (at (move s 0 (* 5 (cos (* 0.10053096491487338 s)))))))
        (assembly root (instance r (track (* 1000 (/ t 1s)))))"#;
    let mut m = model(src);
    m.set_mass("root/r/c", 0.0, [0.0; 3], [[0.0; 3]; 3]);
    m.set_mass("root/r/f", 62.8, [0.0; 3], [[0.0; 3]; 3]);
    // Along s: dz/ds m (z'' - g), z = 5 cos(k s), s = 1000 t.
    let k = 0.10053096491487338;
    let f = |th: f64| -5.0 * k * th.sin() * 62.8 * (-5.0 * k * k * 1e6 * th.cos() + 9810.0) * 1e-6;
    let exact = (0..100000).map(|i| f(i as f64 * std::f64::consts::TAU / 100000.0).abs()).fold(0.0, f64::max);
    let p = peak(&m, 0, (0.0, 1.0));
    assert!(p <= exact * (1.0 + 1e-9) && p >= 0.99 * exact, "peak {p} of {exact}");
    assert_eq!(claims(src)[0].1, Status::Fails);
}

/// A toleranced size on a body that does not move is nothing the effort
/// reads: the claim is decided over the band.
#[test]
fn a_tolerance_the_effort_does_not_read_leaves_it_decided() {
    let c = claims(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (param T 6 (tol -0.1 0.1))
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (body plate (shape (move (box 200 200 T) 0 0 -100)) (attr density 2.7g/cm3))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -130N*mm 130N*mm)))
        (assembly root (instance p (pendulum (* 60 (sin (* 180deg/s t))))) (instance plate plate))"#,
    );
    assert_eq!(c[0].1, Status::Holds, "{c:?}");
}

/// A published mass written as a param follows its value: pinned at 172
/// g, a pendulum held level needs 253.1 N mm and fails its 245 (it held
/// at 238.4, its 162 g default's, whatever --param said); over its
/// tolerance the claim is unknown, sent to single values.
#[test]
fn a_published_mass_param_follows_its_value() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2)
        (param M 162 (tol -10 10))
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr mass (* 1g M)))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -245N*mm 245N*mm)))
        (assembly root (instance p (pendulum 90)))"#;
    let c = claims(src);
    assert!(c[0].1 == Status::Unknown && words(&c[0].2).contains(&"M"), "{c:?}");
    let at = |m: f64| claims_at(src, &[("M", m)])[0].1;
    assert_eq!((at(162.0), at(172.0)), (Status::Holds, Status::Fails));
}

/// An actuator held at a pose a body's shape follows (a boot over the
/// joint): that body is left out, and the claim says so; the rest of the
/// model answers.
#[test]
fn a_held_pose_a_shape_follows_leaves_that_body_out() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (body (boot (angle 0)) (shape (rotate-z (move (box 10 30 10) 0 15 0) angle)) (attr density 1g/cm3))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -1000N*mm 1000N*mm))
          (instance boot (boot angle) (at (rot-x -90))))
        (assembly root (instance p (pendulum 30)))"#;
    assert!(model(src).actuators[0].leaves_out[0].starts_with("root/p/boot"));
    let c = claims(src);
    assert!(c[0].1 == Status::Unknown && c[0].2.contains("root/p/boot"), "{c:?}");
}

/// Past its curve's last speed a stepper's limits cross, and plain effort
/// says so: 17HS19-2004S1 lowering 5 kg at 1 m/s on a 6 mm drum turns at
/// 1591.5 rpm, past its 24 V curve's 1200: the motion needs 5 x 9.81 x 6 =
/// 294.3 N mm of it (no coordinate unpowered), and its effort and its
/// claim are refused, naming the joint, its speed and the speed its curve
/// ends at.
#[test]
fn an_effort_past_a_steppers_curve_says_so() {
    let src = r#"(import "std") (import "../tests/fixtures/dynamics/stepper.ink" as m) (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body load (shape (box 100 100 100)) (attr density 5g/cm3))
        (assembly (hoist (z 0)) (instance l load (at (move 200 0 1500)) (slide z z)) (instance motor (m.drive (rad->deg (/ z 6)))))
        (assembly root (instance h (hoist (* -1000 (/ t 1s)))))"#;
    let m = model(src);
    let at = m.at(0.5);
    assert!((at.effort[0] - 294.3).abs() < 0.5 && at.unpowered.iter().all(Option::is_none), "{:?} {:?}", at.effort, at.unpowered);
    let why = m.refusal(0, &at).unwrap();
    assert!(why.contains("root/h/motor/rotor") && why.contains("1591.5 rpm") && why.contains("curve ends at 1200.0 rpm"), "{why}");
    let c = every_claim_at(src, &[]);
    let claim = c.iter().find(|x| x.0 == "effort-limits:root/h/motor/rotor").unwrap_or_else(|| panic!("{c:?}"));
    assert!(claim.1 == Status::Unknown && claim.2.contains("curve ends at 1200.0 rpm"), "{claim:?}");
}

/// `(effort-rms PATH)`: a 1 kg cart
/// moved 100 mm in 1 s by a cycloid, x = H (t - sin(2 pi t) / 2 pi): its
/// slide's force is m a, a = 2 pi H sin(2 pi t), so its RMS over the
/// second is m 2 pi H / sqrt 2 = 0.444288 N exactly; the midpoint rule on
/// the efforts' steps (chosen for peaks, not for a quadrature) gives
/// 0.443301, 0.22% under. Over the span it holds under a rating 1% above
/// that and fails under one 1% below; a path naming no actuator is an error
/// naming it.
#[test]
fn an_rms_effort_is_a_value_over_the_span() {
    let v = |e: &str| {
        format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
            (body block (shape (box 100 100 100)) (attr mass 1kg))
            (assembly (cart (x 0)) (instance b block (at) (slide x x) (effort-limits -1e3N 1e3N)))
            (assembly root (instance c (cart (* 100 (- (/ t 1s) (/ (sin (* 360deg/s t)) (* 2 pi)))))))
            (require cool (value (effort-rms root/c/b) < {e}N))"#
        )
    };
    let want = 1e-3 * 2.0 * std::f64::consts::PI * 100.0 / 2f64.sqrt();
    let m = model(&v("1"));
    let (rms, band) = m.rms(0, (0.0, 1.0)).unwrap();
    near(rms, want, 3e-3 * want, "m 2 pi H / sqrt 2");
    assert_eq!(band, 0.0);
    assert_eq!(every_claim_at(&v(&format!("{}", want * 1.01)), &[]).iter().find(|c| c.0 == "cool").unwrap().1, Status::Holds);
    assert_eq!(every_claim_at(&v(&format!("{}", want * 0.99)), &[]).iter().find(|c| c.0 == "cool").unwrap().1, Status::Fails);
    let why = load_error(&v("1").replace("root/c/b) <", "root/c) <"));
    assert!(why.contains("root/c is no actuator"), "{why}");
}

/// A 1 kg cart moves 100 mm in a second with equal acceleration, cruise
/// and braking times. Its clamped position has continuous speed.
#[test]
fn a_clamped_stroke_has_no_impulses() {
    let m = model(
        r#"(import "std") (gravity 0 0 0) (t-span 0 2.5s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (cart (x 0)) (instance b block (at) (slide x x) (effort-limits -100N 100N)))
        (def s (- (/ t 1s) 0.25))
        (def ta (/ 1 3))
        (def c1 (clamp s 0 ta))
        (def c2 (clamp (- s ta) 0 ta))
        (def c3 (clamp (- (- s ta) ta) 0 ta))
        (def stroke (* 450 (+ (* 0.5 (square c1)) (* ta c2) (* ta c3) (* -0.5 (square c3)))))
        (assembly root (instance c (cart stroke)))"#,
    );
    let span = (0.0, 2.5);
    let steps = m.steps(span);
    assert!(steps.unchecked.is_empty() && steps.loose.is_empty(), "{:?} unchecked, {:?} loose", steps.unchecked, steps.loose);
    near(m.rms(0, span).unwrap().0, 0.45 * (2.0f64 / 3.0 / 2.5).sqrt(), 1e-3, "RMS force in N");
}

#[test]
fn a_continuous_speed_switch_has_no_impulse() {
    for (capped, moving_for) in [("(min s 0)", 0.375f64), ("(select s s 0)", 0.375), ("(* s (floor (+ 1 s)))", 0.625)] {
        let m = model(&format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body block (shape (box 50 50 50)) (attr mass 1000g))
            (assembly (cart (x 0)) (instance b block (at) (slide x x) (effort-limits -100N 100N)))
            (def s (- (/ t 1s) 0.375))
            (assembly root (instance c (cart (* 100 (square {capped})))))"#
        ));
        let span = (0.0, 1.0);
        let unchecked = m.steps(span).unchecked;
        assert!(unchecked.is_empty(), "{capped}: {unchecked:?}");
        near(m.rms(0, span).unwrap().0, 0.2 * moving_for.sqrt(), 1e-3, "RMS force in N");
    }
}

#[test]
fn a_clamped_cycloid_needs_no_impulse_at_a_full_turn() {
    let cart = |turn| {
        model(&format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 2s)
            (body block (shape (box 50 50 50)) (attr mass 1000g))
            (assembly (cart (x 0)) (instance b block (at) (slide x x) (effort-limits -100N 100N)))
            (def cycle {turn})
            (def u (min (/ t 1s) 1))
            (def stroke (* 100 (- u (/ (sin (* cycle u)) cycle))))
            (assembly root (instance c (cart stroke)))"#
        ))
    };
    let span = (0.0, 2.0);
    let exact = cart("(deg->rad (real 360))");
    let unchecked = exact.steps(span).unchecked;
    assert!(unchecked.is_empty(), "{unchecked:?}");
    near(exact.rms(0, span).unwrap().0, 0.1 * std::f64::consts::PI, 1e-3, "RMS force in N");
    // A rounded angle leaves a nonzero speed just before the stop.
    let rounded = cart("(* 2 pi)");
    assert!(rounded.steps(span).unchecked.iter().any(|&(a, b)| a <= 1.0 && 1.0 <= b));
}

/// An arm eased -45..+45 deg in 2 s by a 3-4-5 polynomial of (clamp (/ t 2s) 0 1),
/// against a DC motor's torque line, whose corner is at no speed: at t = 0.
/// The steps halve toward that corner until no float lies between their ends
/// (down to 5e-324 s), each such step its two ends, read: nothing between the
/// moments to jump, so none is unchecked (the catalog's trial8 slider's servo).
#[test]
fn steps_halved_to_no_float_at_a_clamped_start_are_their_ends() {
    let m = model(
        r#"(import "std") (gravity 0 0 0) (t-span 0 2s)
        (body arm (shape (move (box 110 130 60) 0 -60 30)) (attr mass 1200g))
        (assembly (rig (a 0)) (instance w arm (at) (turn z a)
          (effort-limits (dc-torque-min 2.5N*m 77rpm 1 1 (+ 0deg/s (diff a t))) (dc-torque-max 2.5N*m 77rpm 1 1 (+ 0deg/s (diff a t))))))
        (def u (clamp (/ t 2s) 0 1))
        (def ease (* u u u (+ 10 (* u (- (* 6 u) 15)))))
        (assembly root (instance r (rig (* 45deg (- (* 2 ease) 1)))))"#,
    );
    let steps = m.steps((0.0, 2.0));
    assert!(steps.unchecked.is_empty(), "{:?} unchecked", steps.unchecked);
}

/// A jump in speed, or in position, against a DC motor's torque line read at
/// its speed: the limit's corner is isolated down to steps with no float inside,
/// and the jump in such a step still needs an impulse, never a quiet yes.
#[test]
fn a_jump_isolated_to_no_float_against_a_dc_line_is_unchecked() {
    for (drive, limits) in [
        ("(* 45deg (max (- (/ t 1s) 0.5) 0))", "(effort-limits (dc-torque-min 2.5N*m 77rpm 1 1 (+ 0deg/s (diff a t))) (dc-torque-max 2.5N*m 77rpm 1 1 (+ 0deg/s (diff a t))))"),
        ("(* 1deg (select (- (/ t 1s) 0.5) 0 100))", "(effort-limits (select (- (/ t 1s) 0.5) -2N*m -3N*m) (select (- (/ t 1s) 0.5) 2N*m 3N*m))"),
    ] {
        let src = format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body arm (shape (move (box 110 130 60) 0 -60 30)) (attr mass 1200g))
            (assembly (rig (a 0)) (instance w arm (at) (turn z a) {limits}))
            (assembly root (instance r (rig {drive})))"#
        );
        let m = model(&src);
        let steps = m.steps((0.0, 1.0));
        assert!(steps.unchecked.iter().any(|&(a, b)| a <= 0.5 && 0.5 <= b), "{drive}: {:?} unchecked", steps.unchecked);
        assert!(m.rms(0, (0.0, 1.0)).is_err(), "{drive}");
        let cs = claims(&src);
        assert!(cs.iter().any(|c| c.0.starts_with("effort-limits") && c.1 == Status::Unknown), "{drive}: {cs:?}");
        assert!(!cs.iter().any(|c| c.0.starts_with("effort-limits") && c.1 == Status::Holds), "{drive}: {cs:?}");
    }
}

/// A position or speed jump requires an impulse even if the smooth
/// acceleration on both sides is zero.
#[test]
fn an_rms_over_an_impulse_is_refused() {
    for x in ["(select s 0 10)", "(max s 0)", "(* 1e-9 (max s 0))", "(floor (+ 1 s))", "(select s 0 (select (neg s) 0 10))"] {
        let m = model(&format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body block (shape (box 50 50 50)) (attr mass 1000g))
            (assembly (cart (x 0)) (instance b block (at) (slide x x) (effort-limits -100N 100N)))
            (def s (- (/ t 1s) 0.375))
            (assembly root (instance c (cart {x})))"#
        ));
        let span = (0.0, 1.0);
        assert!(m.steps(span).unchecked.iter().any(|&(a, b)| a <= 0.375 && 0.375 <= b), "{x}");
        let why = m.rms(0, span).err().unwrap();
        assert!(why.contains("an impulse"), "{x}: {why}");
    }
}

#[test]
fn a_turning_arms_peak_on_a_step_boundary_is_measured() {
    let src = r#"(import "std")
        (gravity 0 0 -9.81m/s2)
        (t-span 0 1s)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr mass 162g))
        (assembly (link (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -237N*mm 237N*mm)))
        (assembly root (instance arm (link (* 180 (/ t 1s)))))"#;
    // Its peak, m g r = 0.162 x 9.81 x 150 = 238.383 N mm at t = 0.5 (its centre exactly on its own line, which it
    // repeats about every half turn), is a step's end, read there, and the failure's witness: the moment of least
    // margin, not the first past 237.
    let m = model(src);
    assert!(m.steps((0.0, 1.0)).ends.iter().any(|&t| (t - 0.5).abs() < 1e-12));
    near(m.at(0.5).effort[0], 238.383, 1e-6, "m g r");
    claim_value(src, "effort-limits:root/arm/bar", Status::Fails, -1.383, 1e-5);
}

#[test]
fn a_rating_dip_on_a_step_boundary_is_measured() {
    let src = r#"(import "std")
        (gravity 0 0 0)
        (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (slide-axis (q 0))
          (instance c block (at) (slide x q)
            (effort-limits -2N (* 1N (min 2 (* 0.0016 (square (- (diff q t) 250))))))))
        (assembly root (instance a (slide-axis (* 500 (/ t 1s) (/ t 1s)))))"#;
    // A span failure returns its first witness. The rate is 1000 t and
    // the required effort 1 N, so the margin at that reported time is known.
    let a = asserts_at(src, &[]).remove(0).check;
    let t = a.motion.unwrap().motion.t.unwrap();
    let x = witness(&a.evidence);
    assert!(a.status == Status::Fails && x < 0.0 && (x - (0.0016 * (1000.0 * t - 250.0).powi(2) - 1.0)).abs() < 1e-6, "{x} at {t}");
    let m = model(src);
    let minimum = m.steps((0.0, 1.0)).ends.windows(2).map(|w| m.at(ink_core::eval::span_mid((w[0], w[1]))).margin(0)).fold(f64::INFINITY, f64::min);
    near(minimum, -1.0, 1e-6, "margin at the bottom of the rating dip");
    // A dip 2 microseconds wide, between every moment the steps read: bounded over each step, the one holding
    // it is halved down to it, and a moment inside it fails (it was undecided, read at the moments alone).
    let narrow = src.replace("0.0016", "1000000").replace("t) 250", "t) 251.23");
    let a = asserts_at(&narrow, &[]).remove(0).check;
    let t = a.motion.as_ref().and_then(|m| m.motion.t).unwrap_or(f64::NAN);
    assert!(a.status == Status::Fails && (t - 0.25123).abs() < 1e-5, "{:?} at {t}", a.status);
    let point = &claims_at_moment(&narrow, 0.25123)[0];
    assert!(point.1 == Status::Fails && (point.3 + 1.0).abs() < 1e-6, "{point:?}");
}

/// A steel cylinder r 7.3, h 41.9 spun up at 8726.646 rad/s2 needs
/// I alpha = (pi r^2 h rho) r^2 / 2 alpha = 12.8039 N mm, inside a limit of
/// 12.806 by 0.0021. Its inertia measured 0.03% high read 12.8077, and the
/// limit failed: a measured mass with no band decided it. Within what the
/// measure is known to it is undecided; past it, decided.
#[test]
fn an_effort_within_its_masses_accuracy_of_its_limit_is_not_decided() {
    for (limit, want) in [("12.806", Status::Unknown), ("12.5", Status::Fails), ("13.5", Status::Holds)] {
        let r = asserts_at(&rotor(limit, ""), &[]).into_iter().find(|x| x.name == "effort-limits:root/j/r").unwrap();
        assert_eq!(r.check.status, want, "{limit}: {:?}", r.check.evidence);
    }
}
