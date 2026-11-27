//! Turn symmetry (ink_core::symmetry) of fields built from std: what a
//! motion check over a turning body can reduce to.

use ink_core::symmetry::{Order, turn_symmetry};
use ink_lang::{Doc, Registry};
use std::path::Path;

static CORE: Registry = Registry::new(&[]);

/// Each body of `src` (after std) with its order.
fn orders(src: &str) -> Vec<(String, Order)> {
    let d = Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    d.bodies.iter().map(|b| (b.name.clone(), turn_symmetry(&d.g, b.shape).order)).collect()
}

#[test]
fn std_shapes_show_their_turns() {
    let src = r#"
(fn gear (teeth) (subtract (extrude-z (involute-gear-2d teeth 1 20 0.1) 0 8) (cylinder-z 4 20)))
(body g7 (shape (gear 7)))
(body g15 (shape (gear 15)))
(body g45 (shape (gear 45)))
(body shaft (shape (chamfer-intersect (circle-2d 4) (slab-z 0 30) 0.5)))
(body keyed (shape (subtract (gear 15) (move (box 2 3 20) 4 0 0))))
(body turned (shape (rotate-z (gear 15) 7)))
(body twelve-eight (shape (min (repeat-polar (move (box 2 1 4) 10 0 0) 12) (repeat-polar (move (box 2 1 4) 12 0 0) 8))))
(body six (shape (repeat-polar (move (sphere 1) 5 0 0) 6)))
(body bolt-circle (shape (subtract (cylinder-z 20 5) (pattern-polar (move (cylinder-z 2 10) 15 0 0) 6))))
(body round-in-radians (shape (- (+ (square (deg->rad x)) (square (deg->rad y))) 1)))
"#;
    let want = [
        ("g7", Order::Of(7)),
        ("g15", Order::Of(15)),
        ("g45", Order::Of(45)),
        ("shaft", Order::Every),
        ("keyed", Order::Of(1)),
        ("turned", Order::Of(15)),
        // Two folds of other pitches: nothing (each is off by its own
        // error after a turn).
        ("twelve-eight", Order::Of(1)),
        ("six", Order::Of(6)),
        // Copies at float angles are not exactly symmetric in reals.
        ("bolt-circle", Order::Of(1)),
        // x and y scaled alike by deg->rad's exact pi / 180 (Of(1) when
        // the scaling was not a factor).
        ("round-in-radians", Order::Every),
    ];
    let got = orders(src);
    for (name, o) in want {
        assert_eq!(got.iter().find(|(n, _)| n == name).map(|g| g.1), Some(o), "{name}");
    }
}

#[test]
fn motor_rotors_show_their_flats_and_keys() {
    // Rotors built as motors' are: a stepper's core and chamfered shaft, its
    // D-flat a tagged fillet-intersect over the shaft's end; a keyed shaft,
    // its keyway a warped rounded slot under a floor, cut at a pilot.
    let src = r#"
(body stepper (shape (min (cylinder-z 10 30)
  (tag shaft (subtract (chamfer-intersect (circle-2d 2.5) (slab-z 0 24) 0.3) (tag flat (fillet-intersect (- 2 y) (- 4 z) 1 90)))))))
(body keyed (shape (max (subtract (extrude-z (circle-2d 4) -1 30)
  (tag keyway (max (warp (rounded-rect-2d 2 12 1) x (- z 15) y) (- 3 y)))) (- 0 z))))
"#;
    let got = orders(src);
    for name in ["stepper", "keyed"] {
        assert_eq!(got.iter().find(|(n, _)| n == name).map(|g| g.1), Some(Order::Of(1)), "{name}");
    }
}

#[test]
fn two_folds_whose_difference_is_selected_show_nothing() {
    // Two quarter-turn folds of pitches 1e-10 apart, each within 1e-9 of
    // pi / 2: their difference is a tiny value that changes sign within the
    // quarter, and a select on it makes the field -1 at (7, -7, 0) and +1
    // at its quarter turn (7, 7, 0). Each fold's slack is small; their
    // difference after a turn is not bounded by it.
    let src = r#"
(def angle (atan2 y x))
(def p (/ pi 2))
(def p2 (+ p 1e-10))
(def fold-gap (- (mod angle p) (mod angle p2)))
(def gate (max (+ fold-gap 5e-11) (- (- fold-gap) 1.5e-10)))
(body rotor (shape (max (cylinder-z 12 4) (select gate -1 1))))
"#;
    let got = orders(src);
    assert_eq!(got.iter().find(|(n, _)| n == "rotor").map(|g| g.1), Some(Order::Of(1)));
}
