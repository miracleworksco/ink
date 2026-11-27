//! `ink require` with every tool plugged: a requirement on a part a study
//! or a motion loads is answered by that solver, and by none where the
//! design states the number.

use crate::common::ROOT;
use std::process::Command;

/// The command's exit code and what it printed, on `file` of this folder.
fn require(file: &str) -> (i32, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_ink")).args(["require", &format!("{ROOT}/{file}"), "--plain"]).output().unwrap();
    (out.status.code().unwrap(), String::from_utf8(out.stdout).unwrap() + &String::from_utf8(out.stderr).unwrap())
}

/// The line of `text` that starts (after its indent) with `start`.
fn line<'a>(text: &'a str, start: &str) -> &'a str {
    text.lines().map(str::trim_start).find(|l| l.starts_with(start)).unwrap_or_else(|| panic!("no line '{start}' in\n{text}"))
}

#[test]
fn a_rated_bar_is_answered_by_its_study() {
    let (code, said) = require("ink-study/tests/data/rated-beam.ink");
    assert_eq!(code, 3, "{said}");
    assert!(line(&said, "strong-and-stiff:").starts_with("strong-and-stiff: NO, estimated"), "{said}");
    assert_eq!(line(&said, "yield on"), "yield on root/bar: NO, estimated  (aluminum-6061: yield 240 MPa)");
    assert_eq!(line(&said, "max-stress on"), "max-stress on root/bar: 304.896 MPa, between 248.468 and 343.041 MPa, estimated  (ink solve --study press)");
    assert_eq!(line(&said, "sag on"), "sag on root/bar: YES, estimated  (the shelf under it: 4 mm away)");
    assert_eq!(line(&said, "max-displacement on"), "max-displacement on root/bar: 2.937 mm, between 2.932 and 2.941 mm, estimated  (ink solve --study press)");
    assert!(said.contains("note: its grid error about 25.94 MPa (three grids converging faster than their elements can: 303 on cells of 1 mm, 296.2 on cells of 2 mm)"), "{said}");
}

#[test]
fn a_rated_arm_is_answered_by_its_motion() {
    let (code, said) = require("ink-effort/tests/data/rated-arm.ink");
    assert_eq!(code, 3, "{said}");
    assert!(line(&said, "motor-ok:").starts_with("motor-ok: NO, estimated"), "{said}");
    assert_eq!(line(&said, "motor-torque on"), "motor-torque on root/j/arm: YES, estimated, at 42 values of t  (arm.ink: (effort-limits -150 150))");
    assert_eq!(line(&said, "the size of peak-torque on"), "the size of peak-torque on root/j/arm: 106.904 N*mm, between 106.817 and 111.267 N*mm, estimated  (ink effort on root)");
    assert!(line(&said, "weak-motor on").starts_with("weak-motor on root/j/arm: NO, estimated"), "{said}");
}

#[test]
fn two_assemblies_are_each_answered_over_their_own_motion() {
    // The rated arm in two assemblies of one file: one swings through level, where gravity's torque on its
    // joint is all of 106 N mm; the other stays within 5 degrees of 80 up (27.4 N mm at 75, less the 1.2 its
    // slowing there takes). Each requirement reads its own assembly's motion: the 100 N mm motor fails the
    // first and holds the second.
    let (code, said) = require("ink-effort/tests/data/two-arms.ink");
    assert_eq!(code, 3, "{said}");
    assert!(line(&said, "swing-ok:").starts_with("swing-ok: NO, estimated") && line(&said, "raised-ok:").starts_with("raised-ok: YES, estimated"), "{said}");
    assert!(line(&said, "weak-motor on swing/j/arm").starts_with("weak-motor on swing/j/arm: NO, estimated"), "{said}");
    assert!(line(&said, "weak-motor on raised/j/arm").starts_with("weak-motor on raised/j/arm: YES, estimated"), "{said}");
    assert!(said.contains("the size of peak-torque on swing/j/arm: 106.904 N*mm, between 106.817 and 111.267 N*mm, estimated  (ink effort on swing)"), "{said}");
    assert!(said.contains("the size of peak-torque on raised/j/arm: 26.18 N*mm, between 26.159 and 27.249 N*mm, estimated  (ink effort on raised)"), "{said}");
}

#[test]
fn what_a_design_states_is_held_beside_the_solvers_and_one_that_asks_nothing_says_nothing() {
    // The bar's stress and sag stated on it: each is held beside the study's, and a stress typed under what it
    // works out (200 MPa where the solve gives 305) is no yes: the two disagree, each shown.
    let beam = std::fs::read_to_string(format!("{ROOT}/ink-study/tests/data/rated-beam.ink")).unwrap();
    let dir = env!("CARGO_TARGET_TMPDIR");
    std::fs::copy(format!("{ROOT}/ink-study/tests/data/rated-beam.json"), format!("{dir}/rated-beam.json")).unwrap();
    let stated = beam.replace("(instance bar bar (at))", "(instance bar bar (at) (fact max-stress 200MPa) (fact max-displacement 3.5mm))");
    assert_ne!(stated, beam);
    std::fs::write(format!("{dir}/require-stated.ink"), stated).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_ink")).args(["require", &format!("{dir}/require-stated.ink"), "--plain"]).output().unwrap();
    let said = String::from_utf8(out.stdout).unwrap();
    assert!(out.status.code() == Some(3) && line(&said, "strong-and-stiff:").starts_with("strong-and-stiff: not settled: sources disagree"), "{said}");
    assert!(said.contains("require-stated.ink:21: 200 MPa\n") && said.contains("note: ink solve --study press: 304.896 MPa, between 248.468 and 343.041 MPa\n"), "{said}");
    assert!(line(&said, "sag on root/bar:").starts_with("sag on root/bar: YES"), "{said}");
    // A study and no requirement (beam.ink without its two): nothing is asked, and nothing runs.
    let asks_nothing: String = beam.lines().filter(|l| !l.starts_with("(require ")).map(|l| format!("{l}\n")).collect();
    assert_ne!(asks_nothing.len(), beam.len());
    std::fs::write(format!("{dir}/require-nothing.ink"), asks_nothing).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_ink")).args(["require", &format!("{dir}/require-nothing.ink"), "--plain"]).output().unwrap();
    assert_eq!((out.status.code(), String::from_utf8(out.stdout).unwrap()), (Some(0), String::new()));
}

/// A file's claims, those it writes and those its limit clauses make, each with over how much of what it reads.
#[test]
fn claims_are_answered_with_how_much_they_cover() {
    let asserts = |file: &str| {
        let out = Command::new(env!("CARGO_BIN_EXE_ink")).args(["require", &format!("{ROOT}/{file}"), "--plain"]).output().unwrap();
        (out.status.code().unwrap(), String::from_utf8(out.stdout).unwrap())
    };
    // beam.ink: stiff-enough holds and strong-enough fails, both estimates (its header's expect-status).
    let (code, said) = asserts("examples/beam.ink");
    assert_eq!(code, 3, "{said}");
    assert!(line(&said, "stiff-enough:").starts_with("stiff-enough: YES, estimated  (") && line(&said, "stiff-enough:").ends_with("beam.ink:24)"), "{said}");
    assert!(line(&said, "strong-enough:").starts_with("strong-enough: NO, estimated  ("), "{said}");
    assert!(said.contains("\n  note: (max-displacement press < 4)\n") && said.contains("\n  with: value 2.937, grid error 0.004\n"), "{said}");
    // arm.ink's motor is checked at the motion's 42 moments, and says so.
    let (code, said) = asserts("examples/arm.ink");
    assert!(code == 0 && line(&said, "effort-limits:root/j/arm:").starts_with("effort-limits:root/j/arm: YES, estimated, at 42 values of t  ("), "{said}");
    // Its motor's claim is one nobody wrote (its instance's effort-limits): asked without the flag too, the same.
    assert_eq!(require("examples/arm.ink"), (0, said));
    // A requirement and the asserts beside it: rated-arm.ink's requirement first, then its generated motor claim.
    let (_, said) = asserts("ink-effort/tests/data/rated-arm.ink");
    assert!(
        said.starts_with("PLACEMENT")
            && line(&said, "motor-ok:").starts_with("motor-ok: NO, estimated")
            && line(&said, "effort-limits:root/j/arm:").starts_with("effort-limits:root/j/arm: YES, estimated, at 42 values of t"),
        "{said}"
    );
}

/// The command with these words, from this folder: its exit code, and what it printed and complained of.
fn asked(words: &[&str]) -> (i32, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_ink")).arg("require").args(words).arg("--plain").current_dir(ROOT).output().unwrap();
    (out.status.code().unwrap(), String::from_utf8(out.stdout).unwrap() + &String::from_utf8(out.stderr).unwrap())
}

/// An outside read's cases, from the command line: each was once a wrong or a hidden answer.
#[test]
fn what_an_outside_read_found() {
    // A span asked on the command is the span answered over (the file's is 0 to 2, read at 42 moments); a
    // moment is not asked for here, and saying one is refused, not ignored.
    let (_, said) = asked(&["ink-effort/tests/data/rated-arm.ink", "--t-span", "0:0.1"]);
    let top = line(&said, "motor-torque on");
    assert!(top.starts_with("motor-torque on root/j/arm: YES, estimated, at ") && !top.contains("at 42 values of t"), "{said}");
    let (code, said) = asked(&["ink-effort/tests/data/rated-arm.ink", "--t", "0.5"]);
    assert!(code == 1 && said.starts_with("ink require answers a requirement of ratings over the file's (t-span) or --t-span A:B, not at one moment"), "{said}");
    // A bar under two studies: its stress is the worse of the two, whose band is across its yield. The lighter
    // study does not settle it, and what is printed is the study that answered.
    let (code, said) = asked(&["ink-study/tests/data/two-loads.ink"]);
    assert!(code == 3 && line(&said, "yield on").starts_with("yield on root/bar: not settled: the band reaches across the limit"), "{said}");
    assert!(line(&said, "max-stress on").starts_with("max-stress on root/bar: 304.896 MPa, between 248.468 and 343.041 MPa, estimated  (ink solve --study press)"), "{said}");
    assert!(line(&said, "sag on").starts_with("sag on root/bar: YES, estimated"), "{said}");
    // An effort the solver refuses does not hide a speed over its rating.
    let (code, said) = asked(&["ink-effort/tests/data/speed-crossed.ink"]);
    assert!(code == 3 && line(&said, "free-speed on").starts_with("free-speed on root/j/arm: NO, estimated, at "), "{said}");
    // A rating on a frame the part does not have is answered by nothing; one on a frame across the joint's
    // axis, by the moment about that frame's z (the arm's weight acts about y: none), not the joint's torque.
    let (_, said) = asked(&["ink-effort/tests/data/frame-torque.ink"]);
    assert!(line(&said, "bad-frame on").starts_with("bad-frame on root/j/arm: not available: nothing states it or offers it"), "{said}");
    assert!(line(&said, "cross-axis on").starts_with("cross-axis on root/j/arm: YES, estimated"), "{said}");
    // A curve's dip at a speed the joint passes between two moments read is read: not a yes at 36 moments.
    let (code, said) = asked(&["ink-effort/tests/data/speed-notch.ink"]);
    assert!(code == 3 && line(&said, "pull-out on").starts_with("pull-out on root/j/arm: not settled: the band reaches across the limit"), "{said}");
    assert!(line(&said, "speed on").contains(" deg/s, between ") && said.contains("note: its band reaches to its speed at the moments before and after"), "{said}");
}

/// An unknown says which kind it is by the code that refused, not by its words: a moving body with no mass
/// is something the design lacks; a body the model leaves out is outside what its source covers.
#[test]
fn an_unknown_is_of_the_kind_its_solver_says() {
    let dir = env!("CARGO_TARGET_TMPDIR");
    let said = |name: &str, src: &str| {
        std::fs::write(format!("{dir}/{name}"), src).unwrap();
        let out = Command::new(env!("CARGO_BIN_EXE_ink")).args(["require", &format!("{dir}/{name}"), "--plain"]).output().unwrap();
        (out.status.code().unwrap(), String::from_utf8(out.stdout).unwrap())
    };
    let arm = std::fs::read_to_string(format!("{ROOT}/examples/arm.ink")).unwrap();
    let weightless = arm.replace(" (attr material \"aluminum-6061\")", "");
    assert_ne!(weightless, arm);
    let (code, text) = said("kind-no-mass.ink", &weightless);
    assert!(code == 3 && line(&text, "effort-limits:root/j/arm:").starts_with("effort-limits:root/j/arm: not settled: the design lacks something  ("), "{text}");
    assert!(text.contains("note: these moving bodies have no mass: arm"), "{text}");
    // A block whose shape rides the screw's coordinate: its weight is not counted (dynamics' own case).
    let lift = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body rotor (shape (cylinder-z 3 10)) (attr density 7.8g/cm3))
        (body (blk (h 0)) (shape (move (box 50 50 50) 100 0 (+ 25 h))) (attr density 7.85g/cm3))
        (assembly (lift (z 0)) (instance r rotor (at) (turn z (* 45 z)) (effort-limits -1N*mm 1N*mm)) (instance b (blk z)))
        (assembly root (instance l (lift (* 50 (/ t 1s) (/ t 1s)))))"#;
    let (code, text) = said("kind-left-out.ink", lift);
    assert!(code == 3 && line(&text, "effort-limits:root/l/r:").starts_with("effort-limits:root/l/r: not available: outside what its source covers  ("), "{text}");
    assert!(text.contains("note: its effort leaves out bodies whose shape moves with what moves it"), "{text}");
}

/// A torque rating on a frame turned a quarter about z from the joint's y axis reads the moment about
/// that frame's z: an arm turning about y has none (its weight and inertia act about y), so a 50 N mm
/// limit holds. It is answered with the joint's torque about y, -106.9 N mm: no.
#[test]
fn a_torque_on_a_frame_across_the_joint_is_not_the_joints() {
    let (_, said) = asked(&["ink-effort/tests/data/frame-torque.ink"]);
    assert!(!line(&said, "cross-axis on").starts_with("cross-axis on root/j/arm: NO"), "{said}");
}

/// A block bolted to a carriage a belt pulls: the load where it meets the carriage is its own weight and
/// inertia. The drum moves along the same coordinate and is not on it, but the belt pulls the carriage, not
/// the block, so nothing the drum does reaches it.
#[test]
fn a_part_bolted_to_a_belted_carriage_has_its_own_load() {
    let dir = env!("CARGO_TARGET_TMPDIR");
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
(body drum (shape (cylinder-z 10 10)) (attr density 2.7g/cm3))
(body cart (shape (box 40 40 20)) (attr mass 200g))
(body block (shape (move (box 30 30 30) 0 0 25)) (attr mass 200g))
(assembly (rig (angle 0))
  (instance m drum (at) (turn z angle) (effort-limits -1000N*mm 1000N*mm))
  (instance c cart (at (move 100 10 0)) (slide x (* 10mm (/ angle 1rad))))
  (instance s block (frame-of c))
  (relate m c (belt-span (vec (vec 0 10 0) (vec (+ 100 (* 10 (/ angle 1rad))) 10 0)) 10N 0.95)))
(assembly root (instance r (rig (* 1rad/s t))))
(require seat (value (abs (nth (reaction root/r/s) 2)) < 2N))
"#;
    std::fs::write(format!("{dir}/bolted.ink"), src).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_ink")).args(["require", &format!("{dir}/bolted.ink"), "--plain"]).output().unwrap();
    let said = String::from_utf8(out.stdout).unwrap() + &String::from_utf8(out.stderr).unwrap();
    // At a steady speed: 0.2 kg x 9.81 = 1.962 N, under 2.
    assert!(line(&said, "seat:").starts_with("seat: YES"), "{said}");
    assert!(said.contains("value 1.962"), "{said}");
}
