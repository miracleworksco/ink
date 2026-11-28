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
/// reaches 60.0024 at its turning point, which the claim's steps find, and
/// what its masses are known to moves it by under its 0.01 of room there),
/// claims reading it over the span. Its pin holds it up, but a reaction
/// reads the motion through its own solve, which carries no band: along a
/// motion its masses move, it is not settled. The inverse dynamics along it finds
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
        let held = x.name == "held";
        let through = format!("{:?}", x.check.evidence).contains("through a solve");
        assert!((x.check.status, through) == if held { (Status::Unknown, true) } else { (Status::Holds, false) }, "{}: {:?}", x.name, x.check.evidence);
        // How it is found, its masses' band, and nothing of springs where none is simulated.
        assert!(
            x.check.inexact.notes.iter().any(|n| n.contains("Radau IIA") && n.contains("steps") && n.contains("first-order sensitivity") && !n.contains("spring")),
            "{:?}",
            x.check.inexact.notes
        );
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
/// 0). A commanded actuator claims nothing of its effort, its command
/// clamped to its limits (within them by construction: it held, so saying,
/// until the motion's band read through its solve left it unsettled).
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
    let claims = |extra: &str| asserts_at(&src(extra), &[]).into_iter().map(|x| x.name).collect::<Vec<_>>();
    let a = claims("(command root/p/bar 10)");
    assert!(!a.iter().any(|n| n.starts_with("effort-limits")), "{a:?}");
    assert!(claims("").iter().any(|n| n == "effort-limits:root/p/bar"), "uncommanded, it claims its limits");
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
    // A slide pushed 1 N from 0.4 to 0.5 s, whose speed jumps in its
    // acceleration: its position claim holds (its mass published, exact),
    // and its actuator claims nothing.
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
    assert!(!slide.iter().any(|x| x.name.starts_with("effort-limits")));
    let position = &slide.iter().find(|x| x.name == "slide-position").unwrap().check;
    assert_eq!(position.status, Status::Holds, "{:?}", position.evidence);
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

/// A simulated motion is known to within what its masses are: a 62.832 g
/// cylinder (its mass measured, not stated) pushed 1 N from rest moves
/// 1000 / (2 x 0.062832) = 7957.75 mm in 1 s, which its measured mass reads
/// 7957.66. Read as exact (before), a limit of 7957.71 it passes held. Its
/// mass is known only to its measure's proven bounds, which move it by x
/// dm / m (its centre and inertia move a slide by nothing): a claim within
/// that of its limit is not settled either way, one past it decided. A box
/// whose mass is published (and whose volume the measure finds exactly)
/// moves exactly, and a claim a hair either side of its motion is decided.
#[test]
fn a_simulated_claim_within_its_masses_accuracy_of_its_limit_is_not_settled() {
    let exact = 1000.0 / (2.0 * 1e-3 * std::f64::consts::PI * 400.0 * 50.0 * 1e-3);
    let src = |shape: &str, mass: &str, limit: f64| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body c (shape {shape}) {mass})
            (def q (simulate 0))
            (assembly (axis (s 0)) (instance c c (at) (slide x s)))
            (assembly root (instance a (axis q)))
            (load push root/a/c (force 1N 0 0))
            (require short (value q < {limit}))"#
        )
    };
    let status = |shape: &str, mass: &str, limit: f64| claims_at_moment(&src(shape, mass, limit), 1.0).into_iter().find(|c| c.0 == "short").unwrap().1;
    let (cylinder, measured) = ("(cylinder-z 20 50)", "(attr density 1g/cm3)");
    assert_eq!(status(cylinder, measured, exact * (1.0 - 5e-6)), Status::Unknown);
    let sim = crate::simulation(&doc(&src(cylinder, measured, 0.0))).unwrap().unwrap();
    let (x, band) = (sim.at(0, 1.0).unwrap()[0], sim.mass_band(0, 1.0));
    // Its band is its mass's share of it, at most the 2% a mass is measured to.
    assert!((x - 7957.66).abs() < 0.01 && band[0] > 1e-5 * x && band[0] < 0.02 * x, "{x} within {band:?}");
    assert!((band[1] - 2.0 * band[0]).abs() < 1e-6 * band[0] && (band[2] - 2.0 * band[0]).abs() < 1e-6 * band[0], "x = t^2 / 2m, v = t / m, a = 1 / m: {band:?}");
    for (at, want) in [(1.5, Status::Holds), (0.5, Status::Unknown), (-0.5, Status::Unknown), (-1.5, Status::Fails)] {
        assert_eq!(status(cylinder, measured, x + at * band[0]), want, "{at} of its band past it");
    }
    let (box_, published, exact) = ("(box 50 50 50)", "(attr mass 100g)", 1000.0 / (2.0 * 0.1));
    assert_eq!(crate::simulation(&doc(&src(box_, published, 0.0))).unwrap().unwrap().mass_band(0, 1.0), [0.0; 5]);
    assert_eq!((status(box_, published, exact * (1.0 + 1e-9)), status(box_, published, exact * (1.0 - 1e-9))), (Status::Holds, Status::Fails));
}

/// A body's published mass is its own, exact: the model's mass is it, its
/// moments the measured ones scaled to it. Read as the density it implies
/// (the mass over the volume ink measures) times the volume the dynamics
/// measures, which differ on a curved shape, a 62.832 g cylinder so
/// published, pushed 1 N from rest, read 7958.56 mm at 1 s against
/// 1000 / (2 x 0.062832) = 7957.75, 1e-4 light, and a claim between the two
/// was decided on the wrong side. Its centre and inertia move a slide by
/// nothing: no band, and a claim a hair either side of its motion decided.
#[test]
fn a_published_mass_moves_its_body_as_published() {
    let m = 1e-3 * std::f64::consts::PI * 400.0 * 50.0;
    let src = |limit: f64| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body c (shape (cylinder-z 20 50)) (attr mass {m}g))
            (def q (simulate 0))
            (assembly (axis (s 0)) (instance c c (at) (slide x s)))
            (assembly root (instance a (axis q)))
            (load push root/a/c (force 1N 0 0))
            (require short (value q < {limit}))"#
        )
    };
    let exact = 1e6 / (2.0 * m);
    let sim = crate::simulation(&doc(&src(0.0))).unwrap().unwrap();
    near(sim.at(0, 1.0).unwrap()[0], exact, 1e-6 * exact, "the cylinder's travel at its published mass");
    assert_eq!(sim.mass_band(0, 1.0), [0.0; 5]);
    let status = |limit: f64| claims_at_moment(&src(limit), 1.0).into_iter().find(|c| c.0 == "short").unwrap().1;
    assert_eq!((status(exact * (1.0 + 1e-9)), status(exact * (1.0 - 1e-9))), (Status::Holds, Status::Fails));
}

/// An assembly's published mass is spread over its bodies by their measured
/// volumes: a 20 by 50 mm cylinder and a 15 mm ball weighing 100 g together,
/// pushed 1 N from rest on one slide, move 1e6 / 200 = 5000 mm in 1 s. Each
/// body's share is known only as far as the volumes are (the density's
/// bounds over the measured total, Doc::provided_bounds, and its own): read
/// as exact, the two moved 5000.54 mm with no band, and a claim between
/// that and 5000 was decided on the wrong side. Now they move within 83.6
/// mm of it (each share banded alone, though the two sum to 100 g exactly):
/// a claim within the band either side is not settled, one past it decided.
#[test]
fn an_assemblys_published_mass_split_by_measured_volumes_is_banded() {
    let src = |limit: f64| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body a (shape (cylinder-z 20 50)))
            (body b (shape (sphere 15)))
            (assembly pair (attr mass 100g) (instance a a) (instance b b (at (move 100 0 0))))
            (def q (simulate 0))
            (assembly (axis (s 0)) (instance p pair (at) (slide x s)))
            (assembly root (instance a (axis q)))
            (load push root/a/p (force 1N 0 0))
            (require short (value q < {limit}))"#
        )
    };
    let exact = 5000.0;
    let sim = crate::simulation(&doc(&src(0.0))).unwrap().unwrap();
    let (x, band) = (sim.at(0, 1.0).unwrap()[0], sim.mass_band(0, 1.0)[0]);
    assert!(band > 0.0 && (x - exact).abs() <= band && band < 0.1 * exact, "{x} within {band} of {exact}");
    let status = |limit: f64| claims_at_moment(&src(limit), 1.0).into_iter().find(|c| c.0 == "short").unwrap().1;
    for (by, want) in [(1e-9, Status::Unknown), (-1e-9, Status::Unknown), (1.5 * band, Status::Holds), (-1.5 * band, Status::Fails)] {
        assert_eq!(status(x + by), want, "{by} past its motion");
    }
}

/// A belt carried along a coordinate moves with it its mass, which is
/// measured (its round section, swept along its path), known to that
/// measure's bounds: a 1 mm cord looped round two 10 mm pulleys, its
/// material moving as a 10 g cart slides, pushed 1 N from rest. The cart's
/// mass is published, exact; the cord's 0.83 g is known to 0.2%, and moves
/// the cart by x dm / m (dm what its mass is known to): 7.35 mm of its
/// 46,186 at 1 s, which its simulated motion carries (it read no band), so
/// a claim within it of the motion is not settled.
#[test]
fn a_carried_belts_mass_is_in_a_simulated_motions_band() {
    let src = |limit: f64| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (def r 10) (def c 100)
            (def P (path-loop (vec (path-arc (vec 0 0 0) (vec 0 0 1) r (vec 0 -1 0) pi) (path-line (vec 0 r 0) (vec (neg c) r 0))
                                   (path-arc (vec (neg c) 0 0) (vec 0 0 1) r (vec 0 1 0) pi) (path-line (vec (neg c) (neg r) 0) (vec 0 (neg r) 0)))
                              (vec 0 0 1)))
            (body cart (shape (box 20 20 20)) (attr mass 10g))
            (body (cord (s 0)) (shape (path-sweep (- (sqrt (+ (* x x) (* y y))) 1) P)) (attr density 1g/cm3)
              (attr carried (carried (- (sqrt (+ (* y y) (* z z))) 1) P (path-length P) s)))
            (def q (simulate 0))
            (assembly (axis (s 0)) (instance c cart (at) (slide x s)) (instance b (cord s)))
            (assembly root (instance a (axis q)))
            (load push root/a/c (force 1N 0 0))
            (require short (value q < {limit}))"#
        )
    };
    let sim = crate::simulation(&doc(&src(0.0))).unwrap().unwrap();
    let m = sim.model();
    let (lam, err, len) = m.strip("root/a/b");
    let sources = m.mass_sources(0.0, &[[0.0, 0.0, 0.0]], &[]).unwrap();
    assert_eq!(sources.len(), m.mass_source_count() + 1, "the cord's one source after the cart's");
    let (x, band) = (sim.at(0, 1.0).unwrap()[0], sim.mass_band(0, 1.0)[0]);
    let mass = 10.0 + lam * len;
    near(x, 1e6 / (2.0 * mass), 1e-6 * x, "the cart and the cord's travel");
    assert!(err > 1e-4 * lam && err < 0.05 * lam, "the cord's mass per mm {lam} known to {err}");
    near(band, x * err * len / mass, 1e-3 * band, "the cord's mass's share of the motion");
    let status = |limit: f64| claims_at_moment(&src(limit), 1.0).into_iter().find(|c| c.0 == "short").unwrap().1;
    for (by, want) in [(0.5, Status::Unknown), (-0.5, Status::Unknown), (1.5, Status::Holds), (-1.5, Status::Fails)] {
        assert_eq!(status(x + by * band), want, "{by} of its band past it");
    }
}

/// A reaction past where the call's time limit stopped the simulation is
/// refused, not read: its moments there are unsolved, and a claim on one
/// (plain, or through (max 0N ...), which reads past a NaN as 0) was
/// decided on a NaN (it failed, and the clipped one held). Unlimited, the
/// simulation steps on from where it stopped, and both are decided as far
/// as the bar's masses let them: not settled (the pin carries some 2 N,
/// but a reaction reads the motion's band through its own solve, which
/// carries none), the time limit no longer why.
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
    // With time, the kept simulation steps on and computes a finite reaction at t = 0.7.
    let m2 = super::reacting(pendulum, &["root/p/bar"], true).unwrap();
    let late2 = m2.at(0.7);
    assert!(late2.unsolved.is_none() && m2.reaction_refused(0, &late2).is_none(), "resumed at 0.7: {:?} {:?}", late2.unsolved, m2.reaction_refused(0, &late2));
    let w = &late2.reactions[0];
    assert!(w.force.iter().chain(&w.moment).all(|v| v.is_finite()) && w.force[0].abs() < 1000.0, "resumed reaction: {:?} {:?}", w.force, w.moment);
    for (name, status, evidence) in statuses() {
        let why = format!("{evidence:?}");
        assert!(status == Status::Unknown && why.contains("through a solve") && !why.contains("time limit"), "{name}: {status:?} {why}");
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

/// A spinning body read round but for what its measure leaves (a steel disc
/// but for a 1 mm hole 20 mm off its axle, the hole's imbalance measured, its
/// error 0.0153 g mm at any angle across the axle: Body::unknown) carries that
/// error in its simulated motion as two sources of its own, the centre off
/// the axle each way across it, and four of its second moment: at rest on a
/// level axle the one across gravity turns it by the error times g, the
/// other by nothing. Let go, it swings by its hole's imbalance, 4.39 deg in
/// 1 s, known to within 1.5%.
#[test]
fn a_spinning_bodys_departure_within_its_measure_is_in_its_band() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body disc (shape (subtract (cylinder-y 40 10) (move (cylinder-y 1 20) 20 0 0))) (attr density 7.8g/cm3))
        (def a (simulate 0))
        (assembly (wheel (angle 0)) (instance d disc (at) (turn y angle)))
        (assembly root (instance w (wheel a)))"#;
    let sim = crate::simulation(&doc(src)).unwrap().unwrap();
    let m = sim.model();
    let first = m.made_round.iter().find(|r| r.body == "root/w/d").expect("the disc read round but for its hole").first;
    let sources = m.mass_sources(0.0, &[[0.0; 3]], &[]).unwrap();
    assert_eq!(sources.len(), m.mass_source_count() + 6);
    let across: Vec<f64> = sources[sources.len() - 6..][..2].iter().map(|d| d[0].abs()).collect();
    let turned = first * 9.81e-3 * std::f64::consts::PI / 180.0;
    assert!(across.iter().any(|x| (x - turned).abs() < 1e-9 * turned) && across.iter().any(|x| *x < 1e-9 * turned), "{across:?} against {turned} N mm a degree");
    let (x, band) = (sim.at(0, 1.0).unwrap()[0], sim.mass_band(0, 1.0)[0]);
    assert!((x + 4.39).abs() < 0.01 && band > 0.01 * x.abs() && band < 0.02 * x.abs(), "{x} within {band}");
}
