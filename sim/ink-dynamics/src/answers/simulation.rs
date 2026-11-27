//! Simulated motion, its claim witnesses, and when it is evaluated.

use super::{BAR, REGISTRY, asserts_at, claim_value, claims_at_moment, doc, every_claim_at, load_error, model, near, witness};
use crate::{Simulated, build, simulate};
use ink_lang::Doc;
use ink_lang::claims::{Evidence, Status};
use std::path::Path;
use std::sync::Arc;

/// A simulated motion is an ordinary motion: a pendulum let go at 60
/// degrees (no actuator) never swings past 60, but for its measured centre
/// of mass sitting a few thousandths of a mm off its axis (its far swing
/// reaches 60.0024 at its turning point, which the claim's steps find),
/// and its pin holds it up,
/// claims reading it over the span; the inverse dynamics along it finds
/// nothing unpowered (what no actuator supplies along a simulated
/// coordinate is its interpolation's error); an arm on a DC motor commanded at full
/// voltage needs, by inverse dynamics along its simulated motion, the
/// motor's line at each moment's speed. Answers reading it say how it was
/// found.
#[test]
fn a_simulated_motion_is_read_as_any_motion() {
    let pendulum = String::from(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (def swing (simulate 60))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle)))
        (assembly root (instance p (pendulum swing)))
        (require within (value (abs swing) < 60.01))
        (require held (value (nth (reaction root/p/bar (at)) 2) > 0))"#,
    );
    let a = asserts_at(&pendulum, &[]);
    for x in &a {
        assert_eq!(x.check.status, Status::Holds, "{}: {:?}", x.name, x.check.evidence);
        // How it is found, and nothing of springs where none is simulated.
        assert!(x.check.inexact.notes.iter().any(|n| n.contains("Radau IIA") && n.contains("steps") && !n.contains("spring")), "{:?}", x.check.inexact.notes);
    }
    let m = model(&pendulum);
    for t in [0.1, 0.55, 1.3, 1.9] {
        assert!(m.at(t).unpowered[0].is_none(), "at {t}: {:?}", m.at(t).unpowered);
    }
    let arm = format!(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s) {BAR}
        (assembly (arm (angle 0)) (instance bar bar (at (rot-x 90)) (turn z angle)
          (effort-limits (dc-torque-min 400N*mm 60rpm 1 1 (+ 0deg/s (diff angle t))) (dc-torque-max 400N*mm 60rpm 1 1 (+ 0deg/s (diff angle t))))))
        (assembly root (instance a (arm (simulate -90))))
        (command root/a/bar 1e9)"#
    );
    let m = model(&arm);
    for t in [0.05, 0.3, 0.8] {
        let at = m.at(t);
        near(at.effort[0], at.hi[0], 1e-5 * 400.0, &format!("the arm's effort at {t} s, its line at its speed"));
        assert!(at.rate[0] > 0.0);
    }
}

/// A simulation that cannot go on stops there and says so: a flywheel on
/// a 17HS19-2004S1 at full effort speeds up to its 24 V curve's last speed,
/// 1200 rpm, past which its limits give nothing; its motion is read up to
/// there, and past it a claim reading it is refused, naming why. A command
/// on what is no actuator, or on one moving no simulated coordinate, is an
/// error.
#[test]
fn a_simulation_stops_where_its_mechanism_cannot_go_on() {
    let src = |cmd: &str| {
        format!(
            r#"(import "std") (import "../tests/fixtures/dynamics/stepper.ink" as m) (gravity 0 0 -9.81m/s2) (t-span 0 3s)
            (body wheel (shape (cylinder-z 40 10)) (attr density 7.85g/cm3))
            (def turned (simulate 0))
            (assembly (spinner (a 0)) (instance w wheel (at (move 0 0 60)) (turn z a)) (instance motor (m.drive a)))
            (assembly root (instance s (spinner turned)) (instance still (spinner 0) (at (move 200 0 0))))
            {cmd}
            (require reachable (value turned < 1e9))
            (require slow (value (diff turned t) < 1e9))"#
        )
    };
    let early = claims_at_moment(&src("(command root/s/motor/rotor 1e9)"), 0.01);
    assert_eq!(early.iter().find(|c| c.0 == "reachable").unwrap().1, Status::Holds, "{early:?}");
    assert_eq!(early.iter().find(|c| c.0 == "slow").unwrap().1, Status::Holds, "{early:?}");
    let late = claims_at_moment(&src("(command root/s/motor/rotor 1e9)"), 2.9);
    let slow = late.iter().find(|c| c.0 == "slow").unwrap();
    assert!(slow.1 == Status::Unknown && slow.2.contains("the simulation stopped at t =") && slow.2.contains("past its curve's last speed"), "{slow:?}");
    let none = load_error(&src("(command root/s/w 1e9)"));
    assert!(none.contains("root/s/w is no actuator"), "{none}");
    let idle = load_error(&src("(command root/still/motor/rotor 1e9)"));
    assert!(idle.contains("moves no simulated coordinate"), "{idle}");
}

/// A simulation is stepped only as far as what reads it asks, and once for
/// every load asking the same question (before the fix, a file's check
/// stepped its whole span, twice): the pendulum's document, loaded, has its
/// first step for the launch derivatives; its coordinates' model uses the same
/// simulation, not another; reading 0.5 s steps it to there.
#[test]
fn a_simulation_is_stepped_once_and_only_as_far_as_it_is_read() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (def swing (simulate 57))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle)))
        (assembly root (instance p (pendulum swing)))"#;
    let d = doc(src);
    let sim = |d: &Doc| d.forms::<Simulated>().next().unwrap().run.sim().unwrap();
    let first = sim(&d);
    let reached = |s: &simulate::Simulation| s.reach(0.0, |traj| *traj.t.last().unwrap());
    let knots = |s: &simulate::Simulation| s.reach(0.0, |traj| traj.t.len());
    assert_eq!(knots(&first), 2, "loading reads the first step for the launch derivatives");
    first.at(0, 0.0).unwrap();
    assert_eq!(knots(&first), 2, "reading the launch again reuses the same step");
    // Others made meanwhile (other loads, other threads) don't take it from
    // the document holding it: under load a test's neighbours made more
    // than the 16 kept, and its coordinates' model stepped another.
    for k in 0..17 {
        sim(&doc(&src.replace("(simulate 57)", &format!("(simulate {})", 100 + k))));
    }
    let coordinates = d.reloader().unwrap().with_coordinates("root").unwrap();
    assert!(Arc::ptr_eq(&first, &sim(&coordinates)), "one simulation for the document and its coordinates' model");
    first.at(0, 0.5).unwrap();
    let r = reached(&first);
    assert!((0.5..1.0).contains(&r), "read at 0.5 s, it is stepped to {r}");
    // The question is the file and its overrides: an edited start, or a
    // param overridden, is another simulation.
    let edited = doc(&src.replace("(simulate 57)", "(simulate 58)"));
    assert!(!Arc::ptr_eq(&first, &sim(&edited)), "an edit is another question");
    let longer = Doc::load_str_in(
        &REGISTRY,
        &src.replace("(t-span 0 2s)", "(param long 300) (t-span 0 2s)").replace("box 20 300 10) 0 150 0", "box 20 long 10) 0 150 0"),
        "answers.ink",
        Path::new("."),
        &[("long".into(), 250.0.into())],
    )
    .unwrap();
    let plain = Doc::load_str_in(
        &REGISTRY,
        &src.replace("(t-span 0 2s)", "(param long 300) (t-span 0 2s)").replace("box 20 300 10) 0 150 0", "box 20 long 10) 0 150 0"),
        "answers.ink",
        Path::new("."),
        &[],
    )
    .unwrap();
    assert!(!Arc::ptr_eq(&sim(&longer), &sim(&plain)), "an override is another question");
    // Another target's model reads its own coordinates as written: the
    // pendulum assembly alone, at its argument's default.
    assert!(build(&d.reloader().unwrap(), "pendulum", [0.0, 0.0, -9.81], 0.0, &[], false).is_ok());
}

/// The simulations kept shrink back once their documents close: 20 held at
/// once, then all closed, and the next one made drops the old ones past the
/// 16 kept (before, each call dropped one, so the peak stayed: 63 of 64
/// kept with their trajectories). Others' tests share the cache, so this
/// counts only its own, freed at least 20 - 15.
#[test]
fn the_simulations_kept_shrink_once_their_documents_close() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (def swing (simulate 57))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle)))
        (assembly root (instance p (pendulum swing)))"#;
    let sim = |d: &Doc| d.forms::<Simulated>().next().unwrap().run.sim().unwrap();
    let docs: Vec<Doc> = (0..20).map(|k| doc(&src.replace("(simulate 57)", &format!("(simulate {})", 200 + k)))).collect();
    let weak: Vec<std::sync::Weak<simulate::Simulation>> = docs.iter().map(|d| Arc::downgrade(&sim(d))).collect();
    drop(docs);
    let _next = sim(&doc(&src.replace("(simulate 57)", "(simulate 300)")));
    let freed = weak.iter().filter(|w| w.upgrade().is_none()).count();
    assert!(freed >= 5, "{freed} of the 20 closed simulations freed");
}

/// A claim over a simulated motion is measured where the motion itself
/// turns back (before the fix: read at 32 even moments it missed the far
/// swing): the pendulum let go at 60 degrees reaches about -60 at its half
/// period, 0.48 s, which no even step's middle lands on, so `> -59.99`
/// fails there. A command reading a reaction is refused at load (it read
/// 0); an effort claim on a commanded actuator says it holds by
/// construction.
#[test]
fn claims_over_a_simulation_read_its_own_extremes() {
    let src = |extra: &str| {
        format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
            (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
            (def swing (simulate 60))
            (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -1e4N*mm 1e4N*mm)))
            (assembly root (instance p (pendulum swing)))
            {extra}"#
        )
    };
    let c = every_claim_at(&src("(require far (value swing > -59.99))"), &[]);
    let far = c.iter().find(|x| x.0 == "far").unwrap();
    assert_eq!(far.1, Status::Fails, "{far:?}");
    let why = load_error(&src("(def r (reaction root/p/bar (at))) (command root/p/bar (* 1000 (nth r 1)))"));
    assert!(why.contains("the command on root/p/bar reads a reaction"), "{why}");
    let a = asserts_at(&src("(command root/p/bar 10)"), &[]);
    let effort = a.iter().find(|x| x.name == "effort-limits:root/p/bar").unwrap();
    assert!(
        effort.check.status == Status::Holds && effort.check.inexact.notes.iter().any(|n| n.contains("by construction")),
        "{:?} {:?}",
        effort.check.status,
        effort.check.inexact.notes
    );
    // Its steps are the 32 even ones and its position and speed extrema, whatever the
    // integrator took (one per knot made a geometric claim 19 times
    // slower); its notes say how the motion is found, not how far it
    // had got before the claim read it ("0 steps so far").
    let a = asserts_at(&src("(require far (value swing > -61))"), &[]);
    let far = &a.iter().find(|x| x.name == "far").unwrap().check;
    let steps = far.motion.as_ref().unwrap().motion.steps;
    // Four half periods, each with a position extreme and a speed extreme, plus the span ends.
    // Up to six additional times for both ends of the acceleration, jerk and snap ranges.
    assert!(steps <= 32 + 2 * (8 + 6) + 2, "{steps} steps");
    assert!(far.inexact.notes.iter().any(|n| n.contains("Radau IIA") && !n.contains("so far")), "{:?}", far.inexact.notes);
    // A commanded pulse's effort claim holds by construction, read at even
    // steps (the efforts' own steps could not split across the pulse and
    // left it unknown).
    let pulse = asserts_at(&src("(command root/p/bar (select (- (/ t 1s) 0.4) 0 (select (- (/ t 1s) 0.5) 100 0)))"), &[]);
    assert_eq!(pulse.iter().find(|x| x.name == "effort-limits:root/p/bar").unwrap().check.status, Status::Holds);
    // A slide pushed 1 N from 0.4 to 0.5 s, whose speed jumps in its
    // acceleration: the efforts' steps leave the pulse's edges unchecked,
    // the actuator's claim needs only the steps of its simulated motion.
    let slide = asserts_at(
        r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (axis (s 0)) (instance c block (at) (slide x s) (effort-limits -100N 100N)))
        (def x-slide (simulate 0))
        (assembly root (instance a (axis x-slide)))
        (require slide-position (value x-slide >= -1))
        (command root/a/c (select (- (/ t 1s) 0.4) 0 (select (- (/ t 1s) 0.5) 1 0)))"#,
        &[],
    );
    let effort = &slide.iter().find(|x| x.name == "effort-limits:root/a/c").unwrap().check;
    let position = &slide.iter().find(|x| x.name == "slide-position").unwrap().check;
    assert_eq!(effort.status, Status::Holds, "{:?}", effort.evidence);
    assert_eq!(effort.motion.as_ref().unwrap().motion.steps, position.motion.as_ref().unwrap().motion.steps);
}

/// A simulated coordinate moving nothing with mass, pushed by nothing,
/// keeps its speed (its acceleration is not determined; it was once refused,
/// and so took the effort claim of the stepper axis beside it, whose model
/// reads it), as before springs: both claims hold. One
/// pushed is refused, naming no springs where none is simulated.
#[test]
fn a_simulation_that_cannot_start_leaves_what_does_not_read_it() {
    let src = r#"(import "std") (import "../tests/fixtures/dynamics/stepper.ink" as m) (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body ghost (shape (box 20 20 20)) (attr mass 0))
        (body cart (shape (box 100 100 100)) (attr mass 2kg))
        (assembly (free (q 0)) (instance g ghost (at) (slide x q)))
        (assembly (axis (x 0)) (instance c cart (at (move 200 0 0)) (slide x x)) (instance motor (m.drive (* x 9))))
        (def u (simulate 0))
        (assembly root (instance f (free u)) (instance a (axis (* 20 (- 1 (cos (* 360deg/s t)))))))
        (require u-small (value (abs u) < 1))"#;
    let a = asserts_at(src, &[]);
    let st = |n: &str| a.iter().find(|x| x.name == n).map(|x| (x.check.status, format!("{:?}", x.check.evidence))).unwrap();
    assert_eq!((st("effort-limits:root/a/motor/rotor").0, st("u-small").0), (Status::Holds, Status::Holds));
    let pushed = asserts_at(&src.replace("(require u-small", "(load push root/f/g (force 1N 0 0)) (require u-small"), &[]);
    let why = pushed.iter().find(|x| x.name == "u-small").map(|x| format!("{:?}", x.check.evidence)).unwrap();
    assert!(why.contains("nothing with mass moves along f.q") && !why.contains("spring"), "{why}");
}

#[test]
fn a_simulated_slides_speed_peaks_between_even_samples() {
    let src = r#"(import "std")
        (gravity 0 0 0)
        (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (def q (simulate 0 100))
        (assembly (axis (s 0)) (instance c block (at) (slide x s) (speed-limits 0 125mm/s)))
        (assembly root (instance a (axis q)))
        (load excitation root/a/c (force (* 1N (* 3.2 pi (cos (* 64 180deg/s t)))) 0 0))"#;
    let a = asserts_at(src, &[]).into_iter().find(|a| a.name == "speed-limits:root/a/c").unwrap().check;
    let t = a.motion.unwrap().motion.t.unwrap();
    let v = 100.0 + 50.0 * (64.0 * std::f64::consts::PI * t).sin();
    let x = witness(&a.evidence);
    assert!(a.status == Status::Fails && v > 125.0 && (x - (125.0 - v)).abs() < 1e-3, "{x} at {t}, speed {v}");
    let peak = format!("{src}\n(require maximum-speed (value (diff q t) < 151))");
    claim_value(&peak, "maximum-speed", Status::Holds, 150.0, 1e-3);
}

#[test]
fn a_force_ramp_has_jerk_at_launch() {
    let src = r#"(import "std")
        (gravity 0 0 0)
        (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (def q (simulate 0))
        (assembly (axis (s 0)) (instance c block (at) (slide x s)))
        (assembly root (instance a (axis q)))
        (load ramp root/a/c (force (* 1N (/ t 1s)) 0 0))
        (require jerk-limit (value (abs (diff (diff (diff q t) t) t)) < 500))"#;
    for position in ["0", "1e8"] {
        let cs = claims_at_moment(&src.replace("(simulate 0)", &format!("(simulate {position})")), 0.0);
        let c = cs.iter().find(|c| c.0 == "jerk-limit").unwrap();
        assert!(c.1 == Status::Fails && (c.3 - 1000.0).abs() < 1e-4, "{position}: {c:?}");
    }
}

#[test]
fn a_light_slide_accelerates_next_to_a_heavy_slide() {
    let src = r#"(import "std") (gravity 0 0 0) (t-span 0 0.01s)
(body small (shape (box 5 5 5)) (attr mass 1g))
(body large (shape (box 1e6 1e6 1e6)) (attr mass 1e16g))
(def light (simulate 0))
(def heavy (simulate 5))
(assembly (axis (x 0) (y 0)) (instance small small (at) (slide x x)) (instance large large (at) (slide y y)))
(assembly root (instance m (axis light heavy)))
(load push root/m/small (force 0.001N 0 0))
(require light-moves (value (diff (diff light t) t) > 999))"#;
    claim_value(src, "light-moves", Status::Holds, 1000.0, 1e-4);
}

#[test]
fn a_simulated_acceleration_peak_is_measured() {
    let src = r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
(body block (shape (box 50 50 50)) (attr mass 1000g))
(def q (simulate 0 100))
(assembly (axis (s 0)) (instance c block (at) (slide x s)))
(assembly root (instance m (axis q)))
(load push root/m/c (force (* 1N (+ 2 (sin (* 64 180deg/s t)))) 0 0))
(require acceleration (value (diff (diff q t) t) < 2999))"#;
    claim_value(src, "acceleration", Status::Fails, 3000.0, 1.1);
}

/// A force pulse of 0.1 ms: 1.0e6 mm/s2 while it lasts. The step taken
/// across each of its edges holds both accelerations, and the quintic
/// through its ends, keeping its mean, bulged to 1.113e6 there: `< 1.05e6`
/// failed. Read as its ends' in such a step, it is the pulse's.
#[test]
fn a_simulated_acceleration_at_a_force_step_is_not_read_off_the_curve_across_it() {
    let src = r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
(body block (shape (box 50 50 50)) (attr mass 1000g))
(def q (simulate 0 100))
(assembly (axis (s 0)) (instance c block (at) (slide x s)))
(assembly root (instance m (axis q)))
(load push root/m/c (force (* 1N (select (- (/ t 1s) 0.3) 1 (select (- (/ t 1s) 0.3001) 1000 1))) 0 0))
(require under (value (diff (diff q t) t) < 1000001))
(require reached (value (diff (diff q t) t) < 999999))"#;
    let a = asserts_at(src, &[]);
    assert_eq!(a.iter().find(|x| x.name == "under").unwrap().check.status, Status::Holds);
    claim_value(src, "reached", Status::Fails, 1.0e6, 0.01);
}

#[test]
fn a_simulated_jerk_peak_is_measured() {
    let src = r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
(body block (shape (box 50 50 50)) (attr mass 1000g))
(def q (simulate 0 100))
(assembly (axis (s 0)) (instance c block (at) (slide x s)))
(assembly root (instance m (axis q)))
(load push root/m/c (force (* 1N (+ 2 (sin (+ (* 64 pi (/ t 1s)) 0.7853981633974483)))) 0 0))
(require jerk (value (diff (diff (diff q t) t) t) < 198000))"#;
    claim_value(src, "jerk", Status::Fails, 201062.0, 3500.0);
}

/// A simulation's acceleration has no proven continuity (a force step makes
/// it jump), so a jerk claim every sample satisfies is unknown, never held.
#[test]
fn a_jerk_that_samples_satisfy_is_unknown_where_the_acceleration_may_jump() {
    for force in ["(* 1N (+ 2 (sin (* 64 180deg/s t))))", "(* 1N (select (- (/ t 1s) 0.5) 1 2))"] {
        let src = format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
(body block (shape (box 50 50 50)) (attr mass 1000g))
(def q (simulate 0 100))
(assembly (axis (s 0)) (instance c block (at) (slide x s)))
(assembly root (instance m (axis q)))
(load push root/m/c (force {force} 0 0))
(require below (value (abs (diff (diff (diff q t) t) t)) < 1e15))
(require above (value (diff (diff (diff q t) t) t) > -1e15))"#
        );
        let asked: Vec<_> = asserts_at(&src, &[]).into_iter().filter(|a| a.name == "below" || a.name == "above").collect();
        assert_eq!(asked.len(), 2, "{force}");
        for a in asked {
            assert_eq!(a.check.status, Status::Unknown, "{force}, {}: {:?}", a.name, a.check.evidence);
            assert!(matches!(&a.check.evidence, Evidence::NotCovered(why) if why.contains("rate")), "{force}, {}: {:?}", a.name, a.check.evidence);
        }
    }
}

/// A simulated motion reads measured masses as exact: a 62.832 g cylinder (its mass measured, not stated) pushed 1 N from rest,
/// simulated, moves 1000 / (2 x 0.062832) = 7957.75 mm in 1 s. Its measured mass reads it 7957.66: a limit of
/// 7957.71 it passes, yet the claim holds, the mass read as exact (2026-11-24); within what the mass is known
/// to, it is not settled.
#[test]
#[ignore]
fn a_simulated_claim_within_its_masses_accuracy_of_its_limit_is_not_settled() {
    let exact = 1000.0 / (2.0 * 1e-3 * std::f64::consts::PI * 400.0 * 50.0 * 1e-3);
    let src = format!(
        r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
        (body c (shape (cylinder-z 20 50)) (attr density 1g/cm3))
        (def q (simulate 0))
        (assembly (axis (s 0)) (instance c c (at) (slide x s)))
        (assembly root (instance a (axis q)))
        (load push root/a/c (force 1N 0 0))
        (require short (value q < {}))"#,
        exact * (1.0 - 5e-6)
    );
    let c = claims_at_moment(&src, 1.0);
    let short = c.iter().find(|c| c.0 == "short").unwrap();
    assert_eq!(short.1, Status::Unknown, "{short:?}");
}

/// A reaction past where the call's time limit stopped the simulation is
/// refused, not read: its moments there are unsolved, and a claim on one
/// (plain, or through (max 0N ...), which reads past a NaN as 0) was
/// decided on a NaN (it failed, and the clipped one held). Unlimited, both
/// hold (the pin carries some 2 N),
/// the simulation stepping on from where it stopped.
#[test]
fn a_reaction_past_the_time_limit_is_not_read() {
    // Its own span: the test's limit is the simulation's over it.
    let pendulum = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.737s)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (def swing (simulate 60))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle)))
        (assembly root (instance p (pendulum swing)))
        (require plain (value (nth (reaction root/p/bar (at)) 0) < 1000N))
        (require clipped (value (max 0N (nth (reaction root/p/bar (at)) 0)) < 1000N))"#;
    let statuses = || asserts_at(pendulum, &[]).into_iter().map(|a| (a.name, a.check.status, a.check.evidence)).collect::<Vec<_>>();
    // Limited first: the process keeps the simulation, stepped as far as it got.
    let limit = simulate::LimitAfter::new((0.0, 0.737), 200);
    let limited = statuses();
    // The model's own refusal says so too, not only the coordinate's.
    let m = super::reacting(pendulum, &["root/p/bar"], true).unwrap();
    let late = m.at(0.7);
    assert!(late.unsolved.is_some(), "the moment past the stop is unsolved");
    let why = m.reaction_refused(0, &late).map(|r| r.to_string()).unwrap_or_default();
    assert!(why.contains("time limit"), "{why}");
    drop(limit);
    assert_eq!(limited.len(), 2);
    for (name, status, evidence) in limited {
        assert!(!matches!(status, Status::Holds | Status::Fails), "{name}: {status:?} {evidence:?}");
        assert!(format!("{evidence:?}").contains("time limit"), "{name}: {evidence:?}");
    }
    // With time, the kept simulation steps on, and both hold.
    for (name, status, evidence) in statuses() {
        assert_eq!(status, Status::Holds, "{name}: {evidence:?}");
    }
}

/// Past where the call's time limit stopped the simulation, a claim's steps
/// are not halved to the finest (no moment there solves, however short its
/// step: 4096 steps, each read): a step solving nowhere, its ends and middle,
/// is unchecked whole, and the steps up to the stop are chosen as without
/// it.
#[test]
fn steps_past_the_time_limit_are_unchecked_not_halved() {
    // Its own span: the test's limit is the simulation's over it.
    let pendulum = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.741s)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (def swing (simulate 60))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle)))
        (assembly root (instance p (pendulum swing)))"#;
    let limit = simulate::LimitAfter::new((0.0, 0.741), 200);
    let m = model(pendulum);
    let s = m.steps((0.0, 0.741));
    drop(limit);
    let stop = s.unchecked.iter().map(|u| u.0).fold(f64::INFINITY, f64::min);
    assert!(stop > 0.1 && stop < 0.741, "unchecked from {stop}: {:?}", s.unchecked);
    // Everything past the stop is unchecked, and nothing before it.
    let after: f64 = s.unchecked.iter().map(|u| u.1 - u.0).sum();
    assert!((after - (0.741 - stop)).abs() < 1e-12, "{after} unchecked of {}", 0.741 - stop);
    assert!(s.loose.iter().all(|l| l.1 <= stop), "{:?}", s.loose);
    // Past the stop, a few steps (the stop found by halving the step across
    // it: 12; before, 2181 halved to the finest); before it, the unlimited
    // run's (565 to its 561, the stop's own among them).
    let whole = model(pendulum).steps((0.0, 0.741));
    let before = |e: &[f64]| e.iter().filter(|&&t| t <= stop).count();
    let past = s.ends.len() - before(&s.ends);
    assert!(past <= 20, "{past} steps past the stop");
    assert!(before(&s.ends) <= before(&whole.ends) + 20, "{} steps before the stop, {} unlimited", before(&s.ends), before(&whole.ends));
}
