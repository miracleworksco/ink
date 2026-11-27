//! Efforts read from measured masses through the library: a light
//! member's imbalance beside a heavy one. (The loops closed by a solve
//! that were here left with the kernel's solver.)

use ink_ext as ink;

mod common;

use ink::claims::{ClaimOptions, Status, run_asserts};

/// A 196 kg plate (order 2, measured exactly) spinning at 100 rev/s with a
/// 1 mm^3 lump riding it 5 mm off the axis (m r = 0.03925 g mm): the lump's
/// imbalance swings the plate's reaction by m r w^2 = 0.0155 N. Rounding is
/// each member's own (1e-9 of its mass at its reach), so the plate's size
/// does not zero the lump's. The span ends at a quarter turn, where the
/// measured force witnesses a failure. The reported band is 0.0155 N.
#[test]
fn a_light_members_imbalance_is_kept_beside_a_heavy_one() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.1625s)
        (body plate (shape (box 1000 500 50)) (attr density 7.85g/cm3) (frame rim (at (move 0 5 0))))
        (body lump (shape (box 1 1 1)) (attr density 7.85g/cm3))
        (assembly (spinner (a 0))
          (instance plate plate (at (rot-x 90)) (turn z a) (effort-limits -1e12N*mm 1e12N*mm))
          (instance lump lump (frame-of plate rim)))
        (assembly root (instance s (spinner (* 360 100 (/ t 1s)))))
        (def r (reaction root/s/plate (at)))
        (require fx-small (value (abs (nth r 0)) < 0.006198N))"#;
    let dir = common::Scratch::new("lump");
    let f = dir.join("m.ink");
    std::fs::write(&f, src).unwrap();
    let mut v = ink::view::ModelView::load_in(ink::registry(), &f, &[], 0.0).unwrap();
    v.set_motion(v.doc.t_span).unwrap();
    let c = run_asserts(&mut v, Some("fx-small"), false, ClaimOptions::default()).unwrap().asserts.remove(0).check;
    assert_eq!(c.status, Status::Fails, "{c:?}");
    let motion = c.motion.as_ref().unwrap();
    let t = motion.motion.t.unwrap();
    let omega = 100.0 * std::f64::consts::TAU;
    let force = 1e-6 * 0.03925 * omega.powi(2) * (omega * t).sin().abs();
    assert!(force > 0.006198 && matches!(c.evidence, ink::claims::Evidence::Value(v) if (v - force).abs() < 1e-7), "{c:?}");
    assert!(motion.basis.as_ref().unwrap().contains("root/s/plate (0.0155 N"));
}
