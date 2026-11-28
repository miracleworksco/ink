//! An effort bounded over each step of the motion (model/bounded.rs): every
//! moment's effort inside its step's bound, whatever the mechanism's tilt,
//! gravity's direction or the limits' speed; a peak between the moments
//! read caught; a claim holding between them proven; what the bound does
//! not cover read at its moments as before.

use super::{asserts_at, claims, model};
use crate::Model;
use ink_core::eval::span_mid;
use ink_core::interval::I;
use ink_lang::Steps;
use ink_lang::claims::Status;

/// `n + 1` moments evenly across `span` and each step's middle: each
/// moment's effort (as the model reads it at a moment) inside the bound of
/// every step of the motion holding it; and with each of its bodies'
/// measured moments at either side of its bounds (Model::move_source),
/// inside that bound widened by its masses' band. Returns how far the
/// widest bound reaches past the effort's own range over its step (17
/// moments across it), over the effort's peak: what the bound adds.
fn every_moment_inside(m: &mut Model, a: usize, span: (f64, f64), n: usize) -> f64 {
    let steps: Vec<(f64, f64)> = m.steps(span).ends.windows(2).map(|w| (w[0], w[1])).collect();
    let bounds: Vec<(I, f64)> = steps.iter().map(|&s| m.effort_bound(a, s).unwrap_or_else(|why| panic!("not bounded: {why}"))).collect();
    let across = |s: (f64, f64)| (0..=16).map(move |i| if i == 16 { s.1 } else { s.0 + (s.1 - s.0) * i as f64 / 16.0 });
    let times: Vec<f64> = (0..=n).map(|i| span.0 + (span.1 - span.0) * i as f64 / n as f64).chain(steps.iter().flat_map(|&s| across(s))).collect();
    let holding = |t: f64| steps.iter().zip(&bounds).filter(move |(s, _)| s.0 <= t && t <= s.1);
    let mut peak = 0.0f64;
    for &t in &times {
        let e = m.at_once(t).effort[a];
        peak = peak.max(e.abs());
        assert!(holding(t).count() > 0, "no step holds t {t}");
        for (s, (b, _)) in holding(t) {
            assert!(b.contains(e), "actuator {a} at t {t}, in {s:?}: {e} outside {b:?}");
        }
    }
    for j in 0..m.mass_source_count() {
        for by in [1.0, -1.0] {
            m.move_source(j, by);
            for &t in times.iter().step_by(13) {
                let e = m.at_once(t).effort[a];
                for (s, (b, mass)) in holding(t) {
                    assert!(b.lo - mass <= e && e <= b.hi + mass, "source {j} at {by}: t {t} in {s:?}: {e} outside {b:?} +- {mass}");
                }
            }
            m.move_source(j, -by);
        }
    }
    past_range(m, a, &steps) / peak
}

/// How far the widest of the bounds over `steps` reaches past the effort's
/// own range over its step (17 moments across it).
fn past_range(m: &Model, a: usize, steps: &[(f64, f64)]) -> f64 {
    let across = |s: (f64, f64)| (0..=16).map(move |i| if i == 16 { s.1 } else { s.0 + (s.1 - s.0) * i as f64 / 16.0 });
    let past = |&s: &(f64, f64)| {
        let b = m.effort_bound(a, s).unwrap().0;
        let e: Vec<f64> = across(s).map(|t| m.at_once(t).effort[a]).collect();
        let (lo, hi) = (e.iter().copied().fold(f64::INFINITY, f64::min), e.iter().copied().fold(f64::NEG_INFINITY, f64::max));
        (lo - b.lo).max(b.hi - hi)
    };
    steps.iter().map(past).fold(0.0, f64::max)
}

/// A block sliding along a joint turned two ways, gravity at no axis, a
/// load whose force swings with t: 1,000 moments inside their steps'
/// bounds, each measured moment of its mass either way too.
#[test]
fn a_tilted_slide_under_a_swinging_load_is_inside_its_bounds() {
    let src = r#"(import "std") (gravity 3m/s2 -4m/s2 -8.1m/s2) (t-span 0 1s)
        (body block (shape (box 40 30 20)) (attr density 2.7g/cm3))
        (assembly (cart (s 0)) (instance b block (at (rot-z 25) (rot-x 40)) (slide x s) (effort-limits -50N 50N)))
        (assembly root (instance c (cart (* 80 (sin (* 270deg/s t))))))
        (load push root/c/b (force (* 1N (+ 2 (cos (* 400deg/s t)))) (* 0.5N (sin (* 100deg/s t))) 1N))"#;
    let mut m = model(src);
    let width = every_moment_inside(&mut m, 0, (0.0, 1.0), 1000);
    assert!(width < 0.05, "its steps' bounds reach {width} of its peak past its range");
}

/// An arm swung about a tilted shaft, gravity at no axis; a crank whose
/// offset mass spins up, gravity across its axis: inside their bounds. The
/// crank turns up to 16 degrees in a step, which the bound of its rotation
/// over the step blurs (18% of its peak past its range); quartering its
/// steps shrinks that more than threefold.
#[test]
fn a_tilted_pendulum_and_a_spinning_crank_are_inside_their_bounds() {
    let pendulum = r#"(import "std") (gravity 2m/s2 -3m/s2 -9m/s2) (t-span 0 2s)
        (body arm (shape (move (box 200 20 10) 100 0 0)) (attr material "aluminum-6061"))
        (assembly (joint (angle 0)) (instance arm arm (at (rot-x 20) (rot-z 35)) (turn y angle) (effort-limits -150N*mm 150N*mm)))
        (assembly root (instance j (joint (* 60deg (sin (* 180deg/s t))))))"#;
    let mut m = model(pendulum);
    let width = every_moment_inside(&mut m, 0, (0.0, 2.0), 1000);
    assert!(width < 0.05, "pendulum: its steps' bounds reach {width} of its peak past its range");
    let crank = r#"(import "std") (gravity 0 -9.81m/s2 0) (t-span 0 1s)
        (body crank (shape (move (cylinder-z 8 6) 40 0 0)) (attr density 7.85g/cm3))
        (assembly (spin (a 0)) (instance k crank (at (rot-y 70)) (turn z a) (effort-limits -500N*mm 500N*mm)))
        (assembly root (instance s (spin (+ 30 (* 720 (/ t 1s)) (* 300 (/ t 1s) (/ t 1s))))))"#;
    let mut m = model(crank);
    every_moment_inside(&mut m, 0, (0.0, 1.0), 1000);
    let steps: Vec<(f64, f64)> = m.steps((0.0, 1.0)).ends.windows(2).map(|w| (w[0], w[1])).collect();
    let quarters: Vec<(f64, f64)> =
        steps.iter().flat_map(|&(a, b)| (0..4).map(move |i| (a + (b - a) * i as f64 / 4.0, if i == 3 { b } else { a + (b - a) * (i + 1) as f64 / 4.0 }))).collect();
    let (whole, quartered) = (past_range(&m, 0, &steps), past_range(&m, 0, &quarters));
    assert!(quartered < whole / 3.0, "quartered, its bounds reach {quartered} past its range, whole {whole}");
}

/// Two links, the shoulder's joint turned two ways off the world's axes,
/// the elbow's about the first link's y at its tip, both driven and moving
/// at once, gravity at no axis: each actuator's effort inside its bounds,
/// each one's claim proven where it clears its limits.
#[test]
fn a_two_link_arm_on_composed_joints_is_inside_its_bounds() {
    let src = |shoulder: &str, elbow: &str| {
        format!(
            r#"(import "std") (gravity 1.5m/s2 -2m/s2 -9.5m/s2) (t-span 0 1s)
        (body link (shape (move (box 200 20 10) 100 0 0)) (attr density 2.7g/cm3) (frame tip (at (move 200 0 0))))
        (assembly (arm (sh 0) (el 0))
          (instance a link (at (rot-x 30) (rot-z 20)) (turn z sh) (effort-limits -{shoulder}N*mm {shoulder}N*mm))
          (instance b link (frame-of a tip) (turn y el) (effort-limits -{elbow}N*mm {elbow}N*mm)))
        (assembly root (instance r (arm (* 40 (sin (* 180deg/s t))) (+ 20 (* 70 (/ t 1s) (/ t 1s))))))"#
        )
    };
    let mut m = model(&src("5000", "3000"));
    for a in 0..2 {
        let width = every_moment_inside(&mut m, a, (0.0, 1.0), 1000);
        assert!(width < 0.05, "actuator {a}: its steps' bounds reach {width} of its peak past its range");
    }
    let c = claims(&src("5000", "3000"));
    assert!(c.iter().all(|x| x.1 == Status::Holds && x.3.contains("at every moment of the span")), "{c:?}");
}

/// Two motors each turning a bar by both coordinates at once, one by
/// their sum and one by their difference (a differential): the effort
/// along each coordinate is split between them by the drives' ratios
/// inverted (2 by 2, in intervals), and each actuator's effort is inside
/// its bounds at every moment.
#[test]
fn a_differential_s_split_is_inside_its_bounds() {
    let src = r#"(import "std") (gravity 0m/s2 -6m/s2 -7.5m/s2) (t-span 0 1s)
        (body bar (shape (move (box 120 15 8) 60 0 0)) (attr density 2.7g/cm3))
        (assembly (diff (a 0) (b 0))
          (instance p bar (at) (turn z (+ a b)) (effort-limits -900N*mm 900N*mm))
          (instance q bar (at (move 0 0 60) (rot-x 90)) (turn z (- a b)) (effort-limits -900N*mm 900N*mm)))
        (assembly root (instance d (diff (* 30 (sin (* 90deg/s t))) (* 20 (/ t 1s) (/ t 1s)))))"#;
    let mut m = model(src);
    for a in 0..2 {
        let width = every_moment_inside(&mut m, a, (0.0, 1.0), 1000);
        assert!(width < 0.05, "actuator {a}: its steps' bounds reach {width} of its peak past its range");
    }
    let c = claims(src);
    assert!(c.iter().all(|x| x.1 == Status::Holds && x.3.contains("at every moment of the span")), "{c:?}");
}

/// Eight two-link arms drawn at random (a fixed seed): each joint turned
/// off the world's axes, the second link offset and placed at a random
/// tilt, gravity in a random direction, each coordinate a sum of a sine and
/// a ramp, a load of random direction swinging on the second link. Each
/// actuator's effort is inside its steps' bounds at every moment, and with
/// each measured moment of the masses at either side of its bounds.
#[test]
fn arms_drawn_at_random_are_inside_their_bounds() {
    let mut seed: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut draw = move |lo: f64, hi: f64| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        lo + (hi - lo) * (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    for case in 0..8 {
        let mut n = |lo: f64, hi: f64| format!("{:.4}", draw(lo, hi));
        let src = format!(
            r#"(import "std") (gravity {}m/s2 {}m/s2 {}m/s2) (t-span 0 1s)
        (body link (shape (move (box {} 18 8) {} {} 0)) (attr density 2.7g/cm3) (frame tip (at (move {} 0 0) (rot-x {}) (rot-z {}))))
        (assembly (arm (sh 0) (el 0))
          (instance a link (at (rot-x {}) (rot-y {}) (rot-z {})) (turn z sh) (effort-limits -50000N*mm 50000N*mm))
          (instance b link (frame-of a tip) (turn y el) (effort-limits -50000N*mm 50000N*mm)))
        (assembly root (instance r (arm (+ (* {} (sin (* {}deg/s t))) (* {} (/ t 1s))) (+ (* {} (cos (* {}deg/s t))) (* {} (/ t 1s) (/ t 1s))))))
        (load push root/r/b (force (* {}N (sin (* {}deg/s t))) {}N {}N))"#,
            n(-6.0, 6.0),
            n(-6.0, 6.0),
            n(-9.8, -3.0),
            n(80.0, 220.0),
            n(30.0, 110.0),
            n(-10.0, 10.0),
            n(80.0, 220.0),
            n(-40.0, 40.0),
            n(-40.0, 40.0),
            n(-60.0, 60.0),
            n(-60.0, 60.0),
            n(-60.0, 60.0),
            n(10.0, 60.0),
            n(60.0, 400.0),
            n(-90.0, 90.0),
            n(10.0, 60.0),
            n(60.0, 400.0),
            n(-60.0, 60.0),
            n(-3.0, 3.0),
            n(60.0, 600.0),
            n(-2.0, 2.0),
            n(-2.0, 2.0),
        );
        let mut m = model(&src);
        for a in 0..2 {
            let width = every_moment_inside(&mut m, a, (0.0, 1.0), 400);
            assert!(width.is_finite(), "case {case}, actuator {a}: {width}\n{src}");
        }
    }
}

/// A cart on a trapezoid's moves (its acceleration jumping at each corner):
/// inside its bounds; and a limit 2% over its force (m a, 400 g at 400
/// mm/s2: 0.16 N), inside the 4% a peak between moments could add, is
/// proven, not undecided.
#[test]
fn a_trapezoid_s_corners_are_bounded_and_its_margin_settles() {
    let src = |limit: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1.5s)
        (body cart (shape (box 50 50 50)) (attr mass 400g))
        (assembly (axis (s 0)) (instance c cart (at) (slide x s) (effort-limits -{limit}N {limit}N)))
        (assembly root (instance m (axis (trapezoid-position t 400mm/s2 100mm/s 100mm))))"#
        )
    };
    let mut m = model(&src("1"));
    every_moment_inside(&mut m, 0, (0.0, 1.5), 1500);
    let c = claims(&src("0.1632"));
    assert!(c[0].1 == Status::Holds && c[0].3.contains("at every moment of the span"), "{c:?}");
    assert_eq!(claims(&src("0.159"))[0].1, Status::Fails);
}

/// A 1 kg block held still against a push of 10 N that swells to 15 N for
/// 0.2 ms at t = 0.3141 s. Read at 16 even steps alone, every moment holds
/// (the bump is between them); bounded over each, the step holding the bump
/// is halved down to it, and a limit of 14 N fails at a moment inside it.
/// A limit of 16 N is proven at every moment; so is it over the steps the
/// model chooses, and 14 N fails there too.
#[test]
fn a_bump_between_the_moments_read_is_caught() {
    let src = |limit: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (axis (s 0)) (instance c block (at) (slide x s) (effort-limits -{limit}N {limit}N)))
        (assembly root (instance m (axis 0)))
        (load push root/m/c (force (* 1N (+ 10 (* 5 (exp (neg (square (/ (- t 0.3141s) 0.0002s))))))) 0 0))"#
        )
    };
    let even = Steps { ends: (0..=16).map(|k| k as f64 / 16.0).collect(), unchecked: Vec::new(), loose: Vec::new(), basis: String::new() };
    let tight = model(&src("14"));
    assert!(even.ends.windows(2).all(|w| tight.at(span_mid((w[0], w[1]))).margin(0) > 3.9), "every even step's middle holds");
    let plan = tight.bounded(0).unwrap();
    let b = tight.bound_steps(0, &plan, &even);
    assert!(b.failed, "the bump fails");
    let inside = b.steps.ends.windows(2).map(|w| span_mid((w[0], w[1]))).filter(|t| (t - 0.3141).abs() < 0.0002);
    assert!(inside.clone().any(|t| tight.at(t).margin(0) < 0.0), "a moment read inside the bump: {:?}", inside.collect::<Vec<_>>());
    let wide = model(&src("16"));
    let plan = wide.bounded(0).unwrap();
    let b = wide.bound_steps(0, &plan, &even);
    assert!(!b.failed && b.steps.loose.is_empty() && b.unsettled.is_empty(), "all proven: {} halvings", b.halved);
    let c = claims(&src("16"));
    assert!(c[0].1 == Status::Holds && c[0].3.contains("at every moment of the span"), "{c:?}");
    assert_eq!(claims(&src("14"))[0].1, Status::Fails);
}

/// A 1 kg cart that a smooth pulse of acceleration (a Lorentzian, 40 m/s2
/// at its peak and 0.4 ms wide at half of it) carries 25 mm/s faster: the
/// pulse lies between 16 even steps' middles; its force, 40 N, fails a 35 N
/// limit at a moment the halved steps read inside it, and a 45 N one is
/// proven.
#[test]
fn an_acceleration_pulse_between_the_moments_read_is_caught() {
    // s'' = A / (1 + u^2), u = (t - t0) / w: s = A w^2 (u atan u - ln(1 + u^2) / 2).
    let src = |limit: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (axis (s 0)) (instance c block (at) (slide x s) (effort-limits -{limit}N {limit}N)))
        (def u (/ (- t 0.3141s) 0.0002s))
        (assembly root (instance m (axis (* 0.0016 (- (* u (atan2 u 1)) (* 0.5 (ln (+ 1 (square u)))))))))"#
        )
    };
    let even = Steps { ends: (0..=16).map(|k| k as f64 / 16.0).collect(), unchecked: Vec::new(), loose: Vec::new(), basis: String::new() };
    let tight = model(&src("35"));
    assert!(even.ends.windows(2).all(|w| tight.at(span_mid((w[0], w[1]))).margin(0) > 30.0), "every even step's middle holds");
    let plan = tight.bounded(0).unwrap();
    let b = tight.bound_steps(0, &plan, &even);
    assert!(b.failed, "the pulse fails");
    let inside = b.steps.ends.windows(2).map(|w| span_mid((w[0], w[1]))).filter(|t| (t - 0.3141).abs() < 0.0002);
    assert!(inside.clone().any(|t| tight.at(t).margin(0) < 0.0), "a moment read inside the pulse: {:?}", inside.collect::<Vec<_>>());
    assert_eq!(claims(&src("35"))[0].1, Status::Fails);
    let c = claims(&src("45"));
    assert!(c[0].1 == Status::Holds && c[0].3.contains("at every moment of the span"), "{c:?}");
}

/// examples/arm.ink: its motor's torque peaks at 106.9 N mm just below
/// level. A limit of 108, inside the 4% a peak between moments could add,
/// was undecided; bounded over each step it is proven. 106 fails.
#[test]
fn a_pendulum_s_peak_close_under_its_limit_is_proven() {
    let src = |limit: &str| {
        format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body arm (shape (move (box 200 20 10) 100 0 0)) (attr material "aluminum-6061"))
        (assembly (joint (angle 0)) (instance arm arm (at) (turn y angle) (effort-limits -{limit}N*mm {limit}N*mm)))
        (assembly root (instance j (joint (* 45deg (sin (* 180deg/s t))))))"#
        )
    };
    let c = claims(&src("108"));
    assert!(c[0].1 == Status::Holds && c[0].3.contains("at every moment of the span"), "{c:?}");
    assert_eq!(claims(&src("106"))[0].1, Status::Fails);
}

/// The same arm on a DC motor whose torque falls with its speed
/// (std's dc-torque-max and -min, 60 rpm free): its limits read over every
/// speed a step passes. A stall torque clearing the line is proven; one
/// whose line passes under the effort fails.
#[test]
fn a_motor_s_torque_speed_line_is_held_over_each_step_s_speeds() {
    let src = |stall: &str| {
        format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body arm (shape (move (box 200 20 10) 100 0 0)) (attr material "aluminum-6061"))
        (assembly (joint (angle 0))
          (instance arm arm (at) (turn y angle)
            (effort-limits (dc-torque-min {stall}N*mm 60rpm 1.4 1.4 (+ 0deg/s (diff angle t))) (dc-torque-max {stall}N*mm 60rpm 1.4 1.4 (+ 0deg/s (diff angle t))))))
        (assembly root (instance j (joint (* 45deg (sin (* 180deg/s t))))))"#
        )
    };
    let m = model(&src("300"));
    assert!(m.bounded(0).is_ok(), "{:?}", m.bounded(0).err());
    let c = claims(&src("300"));
    assert!(c[0].1 == Status::Holds && c[0].3.contains("at every moment of the span"), "{c:?}");
    assert_eq!(claims(&src("150"))[0].1, Status::Fails);
}

/// What the bound does not cover is read at its moments as before, saying
/// why: a coupling (a lead screw lifting a block), and two drives sharing
/// one coordinate (their split rests on their limits).
#[test]
fn what_the_bound_does_not_cover_is_read_at_its_moments() {
    let screw = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body rotor (shape (cylinder-z 4 30)) (attr density 7.85g/cm3))
        (body block (shape (box 40 40 20)) (attr density 2.7g/cm3))
        (assembly (lift (a 0))
          (instance rotor rotor (at) (turn z a) (effort-limits -2000N*mm 2000N*mm))
          (instance c block (at (move 0 0 100)) (slide z (* 8 (/ a 360))))
          (relate rotor c (lead-screw 8 0.4 0.3)))
        (assembly root (instance l (lift (* 360 (/ t 1s)))))"#;
    let c = claims(screw);
    assert!(c[0].1 == Status::Holds && c[0].3.contains("not bounded between its moments (couplings'"), "{c:?}");
    let shared = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (axis (s 0))
          (instance c block (at) (slide x s) (effort-limits -5N 5N))
          (instance d block (at (move 0 100 0)) (slide x s) (effort-limits -5N 5N)))
        (assembly root (instance m (axis (* 100 (/ t 1s) (/ t 1s)))))"#;
    let c = claims(shared);
    assert!(c.iter().all(|x| x.1 == Status::Holds && x.3.contains("2 drives share 1 coordinates")), "{c:?}");
}

/// A 1 kg block held against 10 N with a 1e-5 N ripple: the effort read
/// at its moments lies inside a 10.000005 N rating in some steps by less
/// than a peak between moments could add. Bounded over each step, the
/// ripple's crests fail it, and 10.3 N is proven.
#[test]
fn a_ripple_s_crest_fails_a_rating_its_moments_clear() {
    let src = |rating: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
(body block (shape (box 50 50 50)) (attr mass 1000g))
(assembly (axis (s 0)) (instance c block (at) (slide x s) (effort-limits -{rating}N {rating}N)))
(assembly root (instance m (axis 0)))
(load push root/m/c (force (* 1N (+ 10 (* 0.00001 (sin (* 32 180deg/s t))))) 0 0))"#
        )
    };
    let effort = |rating: &str| asserts_at(&src(rating), &[]).into_iter().find(|x| x.name == "effort-limits:root/m/c").unwrap().check;
    for (rating, want) in [("10.000005", Status::Fails), ("10.3", Status::Holds), ("10.5", Status::Holds), ("9.9", Status::Fails)] {
        assert_eq!(effort(rating).status, want, "{rating}");
    }
    let basis = effort("10.3").motion.unwrap().basis.unwrap();
    assert!(basis.contains("at every moment of the span"), "{basis}");
}

/// Regression test (P0-1): a smooth coordinate under a kinked placement
/// (`(turn z (max 0 s))` and a non-monotone turning-point hump
/// `(turn z (max 0 (- s 44.999)))`) has velocity jumps (impulses) where the
/// `max` switches. Symbolic `dt2` drops those impulses, so `sweep` must check
/// both zero jump AND finite first-rate derivatives and leave the step
/// crossing the kink unsettled.
#[test]
fn a_velocity_jump_in_a_placement_is_not_proven() {
    for expr in ["(max 0 s)", "(max 0 (- s 44.999))"] {
        let src = format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
(body hub (shape (box 10 10 10)) (attr mass 1g))
(body arm (shape (move (box 100 20 20) 50 0 0)) (attr mass 200g))
(assembly (axis (s 0))
  (instance d hub (at) (turn z s) (effort-limits -50000N*mm 50000N*mm))
  (instance a arm (at (move 0 100 0)) (turn z {expr})))
(assembly root (instance m (axis (- (* 90 (sin (* 180deg/s t))) 45))))"#
        );
        let c = asserts_at(&src, &[]).into_iter().find(|x| x.name == "effort-limits:root/m/d").unwrap().check;
        let basis = c.motion.as_ref().and_then(|m| m.basis.as_deref()).unwrap_or("");
        assert!(
            !basis.contains("within its limits at every moment of the span") && basis.contains("did not settle"),
            "kinked placement {expr} must leave kink step unsettled: {c:?}"
        );
    }
}

/// Regression test (P0-2 & P1): when a claim has `(accuracy x)` and a moment
/// crosses the limit inside that accuracy (10.2 N against 10 N ± 0.5 N),
/// halving must continue on other open steps so a hidden 3 N bump (12.0 to
/// 12.9 N, past 10.5 N) is still caught and fails — both when the baseline
/// after t = 0.5 is inside the 4% zone (9.9 N) and outside it (9.0 N).
#[test]
fn a_crossing_within_claim_accuracy_keeps_safeguards_on_unsettled_steps() {
    for post in ["9.9", "9.0"] {
        let src = format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
(body block (shape (box 50 50 50)) (attr mass 1000g))
(assembly (axis (s 0)) (instance c block (at) (slide x s) (effort-limits -10N 10N (accuracy 0.5N))))
(assembly root (instance m (axis 0)))
(load push root/m/c (force (* 1N (+ (select (- t 0.5s) 10.2 {post}) (* 3 (max 0 (- 1 (* 10000 (/ (abs (- t 0.7s)) 1s))))))) 0 0))"#
        );
        let c = asserts_at(&src, &[]).into_iter().find(|x| x.name == "effort-limits:root/m/c").unwrap().check;
        assert_eq!(c.status, Status::Fails, "hidden bump over {post} N baseline must fail under accuracy 0.5 N: {c:?}");
    }
}

/// Regression test (P1): the fallback 4% rule (`OFF_LINE`) still protects
/// models the interval bound refuses (here, two drives sharing one
/// coordinate).
#[test]
fn an_unbounded_effort_inside_its_rating_by_less_than_a_peak_between_moments_may_add_is_undecided() {
    let src = |rating: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
(body block (shape (box 50 50 50)) (attr mass 1000g))
(assembly (axis (s 0))
  (instance c block (at) (slide x s) (effort-limits -{rating}N {rating}N))
  (instance d block (at (move 0 100 0)) (slide x s) (effort-limits -{rating}N {rating}N)))
(assembly root (instance m (axis 0)))
(load push root/m/c (force (* 1N (+ 20 (* 0.00002 (sin (* 32 180deg/s t))))) 0 0))"#
        )
    };
    let effort = |rating: &str| asserts_at(&src(rating), &[]).into_iter().find(|x| x.name == "effort-limits:root/m/c").unwrap().check;
    assert_eq!(effort("10.000005").status, Status::Unknown, "10.000005");
    assert_eq!(effort("10.3").status, Status::Unknown, "10.3");
    assert_eq!(effort("10.5").status, Status::Holds, "10.5");
    assert_eq!(effort("9.9").status, Status::Fails, "9.9");
    let basis = effort("10.3").motion.unwrap().basis.unwrap();
    assert!(basis.contains("a peak between the moments read may lie over"), "{basis}");
}
