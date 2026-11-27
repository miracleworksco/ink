//! Spring effort and claims on following error and settling.

use super::{asserts_at, claims_at_moment, doc, load_error, near};
use ink_lang::claims::{Evidence, Status};

/// The same outage on a sprung drive: its effort is its spring's force, and
/// its claim reads it at the same steps, the outage's included.
#[test]
fn a_sprung_drives_limit_jumps_are_steps_too() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body cart (shape (box 100 30 30)) (attr mass 1kg))
        (assembly (axis (s 0)) (instance c cart (at) (slide x s) (effort-limits -1000N (* 1N (select (- (abs (- (/ t 1s) 0.123456)) 0.0001) 0.1 1000))) (stiffness 100N/mm)))
        (assembly root (instance a (axis (* 100 (/ t 1s) (/ t 1s)))))
        (springs) (require lag (value (abs (following-error root/a/c)) < 0.01))"#;
    let a = asserts_at(src, &[]).into_iter().find(|a| a.name == "effort-limits:root/a/c").unwrap();
    assert_eq!(a.check.status, Status::Fails);
}

/// A file's physics is what it says, not what it reads: a 1 kg carriage on (stiffness 100N/mm) commanded
/// x = 100 t^2 mm supplies m a = 1 kg x 0.2 m/s^2 = 0.2 N followed rigidly, and its spring's force with
/// (springs), whether or not the file reads its following error (reading it once switched every effort
/// to the springs': 0.4 N with the lag asserted, 0.2 N without). A following error with no (springs) to
/// lag on is refused, saying to say it.
#[test]
fn a_files_physics_is_what_it_says_not_what_it_reads() {
    let src = |springs: &str, read: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s) {springs}
            (body block (shape (box 50 50 50)) (attr mass 1000g))
            (assembly (axis (x 0)) (instance b block (at) (slide x x) (stiffness 100N/mm) (effort-limits -100N 100N)))
            (assembly root (instance a (axis (* 100 (/ t 1s) (/ t 1s)))))
            {read}"#
        )
    };
    let effort = |src: String| {
        let e = claims_at_moment(&src, 0.5).into_iter().find(|x| x.0 == "effort-limits:root/a/b").unwrap();
        100.0 - e.3
    };
    let read = "(require lag (value (abs (following-error root/a/b)) < 1))";
    near(effort(src("", "")), 0.2, 1e-9, "rigid: m a");
    let (said, said_and_read) = (effort(src("(springs)", "")), effort(src("(springs)", read)));
    assert!((said - said_and_read).abs() < 1e-9 && (said - 0.2).abs() > 1e-3, "with (springs), its spring's whether read or not: {said} and {said_and_read}");
    let why = load_error(&src("", read));
    assert!(why.contains("says no (springs)"), "{why}");
}

/// A claim reads a following error as it reads a simulated coordinate
/// (Simulation, Springs): a 1 kg slide on (stiffness 100N/mm) commanded from 5
/// to 0 mm swings 5 mm either way, at its turning points, which the claim's
/// steps find (within 5.001 holds, within 4.99 fails, within 5.0001, inside
/// its band, is undecided); a joint no spring releases is refused, naming
/// those that are.
#[test]
fn a_claim_reads_a_following_error() {
    let src = |claim: &str| {
        format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.1s)
            (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly (axis (c 0)) (instance b b (at) (slide x c) (stiffness 100N/mm)) (instance f b (at (move 0 300 0))))
            (assembly root (instance a (axis (select (- (/ t 1s) 0.01) 5 0))))
            {claim}"#
        )
    };
    let status = |claim: &str| asserts_at(&src(claim), &[]).into_iter().find(|x| x.name == "lag").map(|x| (x.check.status, x.check.evidence.clone())).unwrap();
    assert_eq!(status("(springs) (require lag (value (abs (following-error root/a/b)) < 5.001))").0, Status::Holds);
    // Its displacement relative to the command jumps when the setpoint
    // steps. A speed read there is undefined, even if multiplied by zero.
    assert_eq!(status("(springs) (require lag (value (* 0 (diff (following-error root/a/b) t)) < 1))").0, Status::Unknown);
    let invalid = src("(springs) (require lag (value (abs (following-error root/a/b)) < 5.001))")
        .replace("(axis (select (- (/ t 1s) 0.01) 5 0))", "(axis (+ (select (- (/ t 1s) 0.01) 5 0) (* 0 (diff (select (- (/ t 1s) 0.02) 0 1) t))))");
    let a = asserts_at(&invalid, &[]);
    assert!(a.iter().any(|r| r.name == "lag" && r.check.status == Status::Unknown));
    // Within its band (100 of its tolerance, 1e-6 of the 5 mm it has
    // swung: 5e-4 mm) of the bound, undecided, either side (before
    // the fix: past it by less than the band, it failed).
    assert_eq!(status("(springs) (require lag (value (abs (following-error root/a/b)) < 5.0001))").0, Status::Unknown);
    assert_eq!(status("(springs) (require lag (value (abs (following-error root/a/b)) < 4.9999))").0, Status::Unknown);
    assert_eq!(status("(springs) (require lag (value (abs (following-error root/a/b)) < 4.99))").0, Status::Fails);
    let (st, why) = status("(springs) (require lag (value (abs (following-error root/a/f)) < 1))");
    assert!(st == Status::Unknown && format!("{why:?}").contains("is no joint released for springs (those are: root/a/b)"), "{st:?} {why:?}");
    // A command reading one is refused at load: the simulation's states
    // provide none (it read 0 there).
    let fed = src("(springs) (def lag (following-error root/a/b)) (command root/a/b (* -10 lag))").replace("(stiffness 100N/mm)", "(stiffness 100N/mm) (effort-limits -100N 100N)");
    let why = load_error(&fed.replace("(axis (select (- (/ t 1s) 0.01) 5 0))", "(axis (simulate 5))"));
    assert!(why.contains("reads a reaction, a tension, an RMS effort or a following error"), "{why}");
    // One where nothing states a stiffness says so.
    let plain = src("(springs) (require lag (value (abs (following-error root/a/b)) < 1))").replace(" (stiffness 100N/mm)", "");
    let why = asserts_at(&plain, &[]).into_iter().find(|x| x.name == "lag").map(|x| format!("{:?}", x.check.evidence)).unwrap();
    assert!(why.contains("nothing under root states a stiffness"), "{why}");
}

/// A settling claim late in a long ringing: 1 kg on 4430
/// N/mm (335 Hz) let go 0.15 um out swings by that through 3 s, but Radau
/// IIA damps it away by 2.4 s (150 tolerances of 1 nm, past the band's
/// fixed 100 then, so `< 1.02e-4` after 2.5 s held): the band widens with
/// the periods stepped, so the claim is undecided, either way.
#[test]
fn a_late_settling_claim_is_undecided_where_the_ringing_is_damped_away() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 3s)
        (body b (shape (box 100 100 100)) (attr mass 1000g))
        (assembly (axis (c 0)) (instance b b (at) (slide x c) (stiffness 4430N/mm)))
        (assembly root (instance a (axis (select (- (/ t 1s) 0.01) 1.5e-4 0))))
        (springs) (def lag (following-error root/a/b))
        (require late-under (value (select (- (/ t 1s) 2.5) 0 (abs lag)) < 1.02e-4))
        (require late-over (value (select (- (/ t 1s) 2.5) 1 (abs lag)) > 1.2e-4))"#;
    let a = asserts_at(src, &[]);
    for n in ["late-under", "late-over"] {
        let c = &a.iter().find(|x| x.name == n).unwrap().check;
        // (an estimate's band, the following error's, either way: the kernel's
        // words for it are its own)
        assert!(c.status == Status::Unknown && format!("{:?}", c.evidence).contains("following-error:root/a/b"), "{n}: {:?} {:?}", c.status, c.evidence);
    }
}

/// An effort claim in a simulation of springs: a 1 kg
/// slide on (stiffness 100N/mm) with effort limits of 15 N, commanded 0.1 (1 -
/// cos omega t) at its own frequency: its rigid demand is at most 10 N (the
/// claim held), but its spring, deviating by -(A omega t / 2) sin omega t,
/// pulls up to 149.31 N, which is what the actuator holding it supplies:
/// the claim fails (at its first witness). An actuator whose coordinate moves what the
/// springs release (an arm on a written slide carrying a sprung one) has
/// its claim unknown, naming them.
#[test]
fn a_springs_effort_claim_reads_its_spring() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.1s)
        (body b (shape (box 100 100 100)) (attr density 1g/cm3))
        (assembly (axis (c 0)) (instance b b (at) (slide x c) (stiffness 100N/mm) (effort-limits -15N 15N)))
        (assembly root (instance a (axis (* 0.1 (- 1 (cos (* (sqrt 100000) (/ t 1s))))))))
        (springs) (def lag (following-error root/a/b))"#;
    let a = asserts_at(src, &[]);
    let e = &a.iter().find(|x| x.name == "effort-limits:root/a/b").unwrap().check;
    let v = match e.evidence {
        Evidence::Value(v) => v,
        ref other => panic!("{other:?}"),
    };
    assert!(e.status == Status::Fails && v < 0.0, "{:?} {v}", e.status);
    assert!(e.inexact.notes.iter().any(|n| n.contains("its spring's force")), "{:?}", e.inexact.notes);
    let peak = crate::simulation(&doc(src)).unwrap().unwrap().springs().joints[0].effort.unwrap();
    assert!((peak - 149.31).abs() < 0.05, "{peak} N");
    // Known to within its stiffness times its following error's band at
    // the span's end (100 N/mm times 1.5e-4 mm), a claim within that of its
    // limit is undecided either side, as is its RMS
    // effort's (63.887 N: one 1e-4 under it failed).
    let status = |limits: f64, rms: &str| {
        let a = asserts_at(
            &format!("{} (require r (value (effort-rms root/a/b) < {rms}N))", src.replace("(effort-limits -15N 15N)", &format!("(effort-limits -{limits}N {limits}N)"))),
            &[],
        );
        let st = |n: &str| a.iter().find(|x| x.name == n).unwrap().check.status;
        (st("effort-limits:root/a/b"), st("r"))
    };
    assert_eq!(status(peak + 0.005, "63.8872"), (Status::Unknown, Status::Unknown));
    assert_eq!(status(peak - 0.005, "63.8874"), (Status::Unknown, Status::Unknown));
    assert_eq!(status(peak + 0.05, "63.92"), (Status::Holds, Status::Holds));
    let behind = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.1s)
        (body b (shape (box 100 100 100)) (attr density 1g/cm3))
        (assembly (axis (c 0))
          (instance arm b (at) (slide y c) (effort-limits -1000N 1000N))
          (instance s b (at (move 0 0 200) (frame-of arm)) (slide x 0) (stiffness 100N/mm)))
        (assembly root (instance a (axis (* 10 (/ t 1s) (/ t 1s)))))
        (springs) (def lag (following-error root/a/s))"#;
    let a = asserts_at(behind, &[]);
    let e = &a.iter().find(|x| x.name == "effort-limits:root/a/arm").unwrap().check;
    assert!(
        e.status == Status::Unknown && format!("{:?}", e.evidence).contains("move what the file's simulation releases for its springs (root/a/s)"),
        "{:?} {:?}",
        e.status,
        e.evidence
    );
}

/// A simulation's springs enter where the file reads a following error, no
/// sooner: a pendulum simulated beside a
/// catalog stepper's axis (its rotor states its holding stiffness) is a
/// rigid simulation, its claims decided and the motor's effort the rigid
/// model's, as before springs (they were refused: the cart tied to the
/// rotor by their drives alone); read a following error and the springs
/// are simulated, the tie refused, saying so.
#[test]
fn springs_enter_a_simulation_where_a_following_error_is_read() {
    let src = |extra: &str| {
        format!(
            r#"(import "std") (import "../tests/fixtures/dynamics/stepper.ink" as m) (gravity 0 0 -9.81m/s2) (t-span 0 1s)
            (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
            (body cart (shape (box 100 100 100)) (attr mass 2kg))
            (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle)))
            (assembly (axis (x 0)) (instance c cart (at (move 200 0 0)) (slide x x)) (instance motor (m.drive (* x 9))))
            (def swing (simulate 60))
            (assembly root (instance p (pendulum swing)) (instance a (axis (* 20 (- 1 (cos (* 360deg/s t)))))))
            (require swing-within (value (abs swing) < 60.01))
            {extra}"#
        )
    };
    let a = asserts_at(&src(""), &[]);
    let status = |n: &str| a.iter().find(|x| x.name == n).map(|x| x.check.status).unwrap();
    assert_eq!((status("swing-within"), status("effort-limits:root/a/motor/rotor")), (Status::Holds, Status::Holds));
    let sim = crate::simulation(&doc(&src(""))).unwrap().unwrap();
    assert!(sim.model().released().is_empty());
    // A file simulating nothing asks no simulation of its efforts (it would
    // be refused: no span to run over).
    let plain = r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body b (shape (box 100 100 100)) (attr density 1g/cm3))
        (assembly (axis (c 0)) (instance b b (at) (slide z c) (effort-limits -15N 15N)))
        (assembly root (instance a (axis (* 10 (/ t 1s)))))"#;
    let c = &asserts_at(plain, &[])[0].check;
    assert_eq!(c.status, Status::Holds, "{:?}", c.evidence);
    let read = src("(springs) (def lag (following-error root/a/motor/rotor))");
    let e = crate::simulation(&doc(&read)).unwrap().err().expect("refused");
    assert!(e.contains("root/a/c (with a.x) move with a coordinate that also commands joints released for springs"), "{e}");
}

/// An actuator on a sprung belt that states no stiffness keeps its
/// coordinate (before the fix: it was released, its effort 0 times its
/// deviation, and a claim the rigid model fails at -40.6 N mm held): the
/// written motion holds it, its belt's springs pull on what it drives, and
/// its effort claim, its load running through them, is unknown, saying so.
#[test]
fn an_actuator_stating_no_stiffness_keeps_its_coordinate() {
    let src = |extra: &str| {
        format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.008s)
            (def R (/ 20 pi))
            (def EA 16903N)
            (body cart (shape (box 40 30 20)) (attr density 1.27g/cm3))
            (body rotor (shape (cylinder-z 6 10)) (attr density 7.85g/cm3) (attr spin-inertia 5478.7g*mm2))
            (body idler (shape (cylinder-z 6 10)) (attr density 2.7g/cm3) (attr spin-inertia 0))
            (assembly (axis (s 300))
              (instance cart cart (at) (slide x s))
              (instance rotor rotor (at) (turn z (rad->deg (/ s R))) (effort-limits -1N*mm 1N*mm))
              (instance idler idler (at (move 600 0 0)) (turn z (rad->deg (/ s R))))
              (relate rotor idler (belt-span (vec (vec 0 R 0) (vec 600 R 0)) 50N 0.95 (tensile-stiffness EA)))
              (relate idler cart (belt-span (vec (vec 600 (neg R) 0) (vec (+ s 5) (neg R) 0)) 50N 0.95 (tensile-stiffness EA)))
              (relate cart rotor (belt-span (vec (vec (- s 5) (neg R) 0) (vec 0 (neg R) 0)) 50N 0.95 (tensile-stiffness EA))))
            (assembly root (instance a (axis (+ 300 (* 10 (- 1 (cos (* 20 180deg/s t))))))))
            {extra}"#
        )
    };
    let rigid = asserts_at(&src(""), &[]);
    let r = &rigid.iter().find(|x| x.name == "effort-limits:root/a/rotor").unwrap().check;
    assert_eq!(r.status, Status::Fails);
    let read = src("(springs) (def lag (following-error root/a/cart)) (require lag-small (value (abs lag) < 1))");
    let sim = crate::simulation(&doc(&read)).unwrap().unwrap();
    assert_eq!(sim.model().released().iter().map(|r| r.0).collect::<Vec<_>>(), ["root/a/idler", "root/a/cart"]);
    let a = asserts_at(&read, &[]);
    let e = &a.iter().find(|x| x.name == "effort-limits:root/a/rotor").unwrap().check;
    assert!(e.status == Status::Unknown && format!("{:?}", e.evidence).contains("its load runs through them"), "{:?} {:?}", e.status, e.evidence);
    // The written motion drives the axis: its carriage lags by microns.
    assert_eq!(a.iter().find(|x| x.name == "lag-small").unwrap().check.status, Status::Holds);
}
