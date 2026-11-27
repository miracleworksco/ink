//! Tension, losses, and reactions through belts, gears, and screws.

use super::{claims_at_moment, model, near, reacting, turned, values_at};
use ink_lang::claims::Status;

/// A motor turning a shaft through an open belt (std belt-drive: equal
/// pulleys of pitch radius 10 mm, 100 mm apart, the belt 110 mm up its
/// shafts, 100 N pretension), the shaft on two 608ZZ at 10 and 70 mm and
/// dragged with 1000 N mm at constant speed, 100 N of effective tension.
/// By hand, installed and locked (Gates' Formulas 5 and 6), its spans
/// carry 100 +- 50 N, pulling the shaft 200 N toward the motor at 110 mm,
/// so the far bearing pushes back 200 x 100 / 60 = 333.333 N and the near
/// one pulls 133.333 N; held by a tensioner, 200 and 100 N, 500 and 200 N
/// on the bearings; dragged with 3000 N mm, more than the locked
/// pretension carries, the slack span goes slack and the tight one takes
/// all 300 N. The motor's shaft carries the pull the other way; its effort
/// is the drag over the belt's efficiency, 0.9, and its shaft passes that
/// (without an efficiency the effort is unknown). The frame holds the two
/// torques' difference, drag (1 / 0.9 - 1), as a couple across the two
/// shafts: that over their 100 mm, across the belt's line, at the belt's
/// plane, so the far bearing takes 100 / 60 of it.
#[test]
fn a_belt_drives_pull_on_its_shafts_by_hand() {
    let src = |t0: &str, drag: f64, eff: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (import "../tests/fixtures/dynamics/bearing.ink" as brg)
            (body rotor (shape (cylinder-z 4 30)) (attr density 7.85g/cm3))
            (body disc (shape (cylinder-z 10 6)) (attr density 2.7g/cm3))
            (body shaft (shape (move (cylinder-z 4 120) 0 0 60)) (attr density 7.85g/cm3))
            (assembly (motor (a 0))
              (instance rotor rotor (at) (turn z a) (effort-limits -4000N*mm 4000N*mm))
              (instance pulley disc (at (move 0 0 110) (frame-of rotor))))
            (assembly spindle (instance shaft shaft) (instance pulley disc (at (move 0 0 110))))
            (assembly (drive (a 0))
              (instance m (motor a))
              (instance b1 brg.bearing (at (move 100 0 10)))
              (instance b2 brg.bearing (at (move 100 0 70)))
              (instance s spindle (at (move 100 0 0)) (turn z a))
              (relate m/pulley s/pulley (belt-drive 10 10 {t0} {eff})))
            (assembly root (instance d (drive (* 360 (/ t 1s)))))
            (load drag root/d/s (torque 0 0 {}N*mm))
            (def sh (reaction root/d/s (by root/d/b1 root/d/b2)))
            (def m (reaction root/d/m/rotor (at)))
            (def spans (tension root/d/m/pulley root/d/s/pulley))
            (require b1 (value (nth (nth sh 0) 0) > -1e9N))
            (require b2 (value (nth (nth sh 1) 0) > -1e9N))
            (require across (value (abs (nth (nth sh 1) 1)) > -1e9N))
            (require motor-pull (value (nth m 0) > -1e9N))
            (require motor-torque (value (nth m 5) > -1e9N*mm))
            (require tight (value (max (nth spans 0) (nth spans 1)) > -1e9N))
            (require slack (value (min (nth spans 0) (nth spans 1)) > -1e9N))"#,
            -drag
        )
    };
    let check = |t0: &str, drag: f64, tight: f64, slack: f64| {
        let of = values_at(&src(t0, drag, "0.9"), 0.25);
        let pull = tight + slack;
        near(of("tight"), tight, 1e-9, "the tight span's tension");
        near(of("slack"), slack, 1e-9, "the slack span's");
        near(of("b2"), pull * 100.0 / 60.0, 1e-9, "the far bearing");
        near(of("b1"), pull - pull * 100.0 / 60.0, 1e-9, "the near bearing");
        let couple = drag * (1.0 / 0.9 - 1.0) / 100.0;
        near(of("across"), couple * 100.0 / 60.0, 1e-9, "across the belt's line: the loss's couple");
        near(of("motor-pull"), -pull, 1e-9, "the motor's shaft holds the belt's pull");
        near(of("motor-torque"), drag / 0.9, 1e-9, "the motor's joint passes the drag over the efficiency");
    };
    check("100N", 1000.0, 150.0, 50.0);
    check("(held 100N)", 1000.0, 200.0, 100.0);
    check("100N", 3000.0, 300.0, 0.0);
    let c = claims_at_moment(&src("100N", 1000.0, "0.9"), 0.25);
    let effort = c.iter().find(|x| x.0 == "effort-limits:root/d/m/rotor").unwrap();
    near(4000.0 - effort.3, 1000.0 / 0.9, 1e-9, "the effort, the drag over the belt's efficiency");
    let c = claims_at_moment(&src("100N", 1000.0, ""), 0.25);
    let effort = c.iter().find(|x| x.0 == "effort-limits:root/d/m/rotor").unwrap();
    assert!(effort.1 == Status::Unknown && effort.2.contains("gives no efficiency, so its loss is unknown"), "{effort:?}");
    // The reactions and tensions carry the loss's friction too (the motor's
    // shaft passes it, the bearings hold its couple): unknown with it.
    assert!(c.iter().filter(|x| x.0 != effort.0).all(|x| x.1 == Status::Unknown && x.2.contains("gives no efficiency")), "{c:?}");
}

/// A 20-tooth pinion on a motor driving a 40-tooth wheel, module 2,
/// pressure angle 20, the wheel on a shaft on two 608ZZ at 10 and 70 mm,
/// its face centred 25 mm up, a 2000 N mm load torque on it at constant
/// speed: by hand (AGMA's, ISO 6336's), Ft = 2000 / 40 = 50 N along the
/// pitch circles' tangent and Fr = Ft tan 20 = 18.1985 N pushing the wheel
/// from the pinion, shared 3 : 1 by the bearings 15 and 45 mm from the
/// face. Against its motion, the pinion's motor gives 1000 N mm over the
/// mesh's efficiency, 0.98; with its motion (the load drives it back), 1000
/// times it. Driving, the frame holds the two torques' difference, 1000 (1
/// / 0.98 - 1) = 20.408 N mm, as a couple across the shafts 60 mm apart:
/// the wheel's bearings take 0.340 N more along the tangent. A 30-tooth gear stacked on the pinion, meshing nothing, and a
/// hub under a wider wheel move nothing: the mesh is where the two faces
/// overlap where the teeth work (it was at the middle of both stages, or of
/// teeth and hub, 33.3/16.7).
#[test]
fn a_spur_gear_pairs_bearings_by_hand() {
    let src_on = |drag: f64, pinion: &str, wheel: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (import "../tests/fixtures/dynamics/bearing.ink" as brg)
            (fn gear (z) (subtract (extrude-z (involute-gear-2d z 2 20 0) -5 5) (cylinder-z 4 20)))
            (body pinion (shape {pinion}) (attr density 7.85g/cm3))
            (body wheel (shape {wheel}) (attr density 7.85g/cm3))
            (body shaft (shape (move (cylinder-z 4 80) 0 0 40)) (attr density 7.85g/cm3))
            (assembly spindle (instance shaft shaft) (instance wheel wheel (at (move 0 0 25))))
            (assembly (pair (a 0))
              (instance p pinion (at (move 0 0 25)) (turn z a) (effort-limits -2000N*mm 2000N*mm))
              (instance b1 brg.bearing (at (move 60 0 10)))
              (instance b2 brg.bearing (at (move 60 0 70)))
              (instance s spindle (at (move 60 0 0)) (turn z (mesh-turn 20 40 0 a)))
              (relate p s/wheel (gear-mesh 20 40 2 20 0 0.98)))
            (assembly root (instance g (pair (* -360 (/ t 1s)))))
            (load drag root/g/s (torque 0 0 {drag}N*mm))
            (def sh (reaction root/g/s (by root/g/b1 root/g/b2)))
            (require b1x (value (nth (nth sh 0) 0) > -1e9N))
            (require b1y (value (nth (nth sh 0) 1) > -1e9N))
            (require b2x (value (nth (nth sh 1) 0) > -1e9N))
            (require b2y (value (nth (nth sh 1) 1) > -1e9N))
            (require tooth (value (max (nth (tension root/g/p root/g/s/wheel) 0) (nth (tension root/g/p root/g/s/wheel) 1)) > -1e9N))
            (require idle (value (min (nth (tension root/g/p root/g/s/wheel) 0) (nth (tension root/g/p root/g/s/wheel) 1)) > -1e9N))"#
        )
    };
    let src = |drag: f64| src_on(drag, "(gear 20)", "(gear 40)");
    let (ft, fr) = (50.0, 50.0 * 20f64.to_radians().tan());
    let of = values_at(&src(-2000.0), 0.25);
    // The forces whole, exactly; their 3 : 1 split to where the pitch point
    // is, the middle of the faces' overlap, read by bounds: an outer box,
    // each face's within a quarter mm, so 50 x 0.25 / 60 = 0.21 N here.
    let split = 50.0 * 0.25 / 60.0;
    near(of("b1x") + of("b2x"), -fr, 1e-9, "the separating force, carried by the bearings");
    let ft = ft + 1000.0 * (1.0 / 0.98 - 1.0) / 60.0;
    near(of("b1y") + of("b2y"), ft, 1e-9, "the tangential force, and the loss's couple");
    near(of("b1x"), -0.75 * fr, split, "the near bearing's share of the separating force");
    near(of("b2x"), -0.25 * fr, split, "the far bearing's");
    near(of("b1y"), 0.75 * ft, split, "the near bearing's share of the tangential force");
    near(of("b2y"), 0.25 * ft, split, "the far bearing's");
    let stacked = values_at(&src_on(-2000.0, "(min (gear 20) (move (gear 30) 0 0 10))", "(gear 40)"), 0.25);
    near(stacked("b1y"), 0.75 * ft, split, "a compound pinion: the near bearing's share, at its 20-tooth stage's face");
    near(stacked("b2y"), 0.25 * ft, split, "the far bearing's");
    // A hub as wide as the pinion's root (17.2 of 17.5) above its teeth,
    // under a wheel twice as wide: the faces are read where the wheel's
    // teeth work, so the mesh stays at the pinion's teeth (it was at the
    // middle of teeth and hub, 33.3/16.7).
    let hub = "(min (gear 20) (move (subtract (cylinder-z 17.2 20) (cylinder-z 4 30)) 0 0 15))";
    let hubbed = values_at(&src_on(-2000.0, hub, "(move (subtract (extrude-z (involute-gear-2d 40 2 20 0) -10 10) (cylinder-z 4 30)) 0 0 5)"), 0.25);
    near(hubbed("b1y"), 0.75 * ft, split, "a hubbed pinion under a wider wheel: the near bearing's share, at its teeth");
    near(hubbed("b2y"), 0.25 * ft, split, "the far bearing's");
    near(of("tooth"), 50.0 / 20f64.to_radians().cos(), 1e-9, "the tooth force along the line of action: what the wheel needs");
    near(of("idle"), 0.0, 1e-9, "the idle flank");
    let effort = |drag: f64| {
        let c = claims_at_moment(&src(drag), 0.25);
        let e = c.into_iter().find(|x| x.0 == "effort-limits:root/g/p").unwrap();
        // The margin is the nearer limit's: against +-2000.
        if e.3 < 2000.0 { 2000.0 - e.3 } else { e.3 - 2000.0 }
    };
    near(effort(-2000.0).abs(), 1000.0 / 0.98, 1e-9, "driving: the ideal 1000 N mm over the efficiency");
    near(effort(2000.0).abs(), 1000.0 * 0.98, 1e-9, "driven back: the ideal times it");
}

/// A 15-tooth straight bevel pinion on a motor (axis z) driving a 45-tooth
/// wheel (axis x), module 1.5 at the heel, face 10, their apexes at the
/// origin, 2000 N mm on the wheel: by hand (AGMA 2003, at the mean cone
/// distance Rm = 1.5 x 45 / (2 sin 71.565) - 5 = 30.5756), Ft = 2000 /
/// (Rm sin 71.565) = 68.9499 N along the tangent, and Ft tan 20 parting
/// them, which is radial Ft tan 20 cos d and axial Ft tan 20 sin d on each
/// gear (d its pitch angle): the pinion's 23.8079 and 7.9360, the
/// wheel's the other way round. The load drives the pinion back: the
/// mesh passes its efficiency, 0.98, of what the wheel gives, so the tooth
/// forces are 0.98 of those, the pinion's motor gives the 666.667 N mm
/// times it, and the wheel's friction, 2% of 2000 N mm about its axis,
/// pairs with a moment on the pinion about that axis (square to its own).
#[test]
fn a_bevel_pairs_forces_by_hand() {
    let pair = r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
        (fn bevel (n d) (max (straight-bevel-gear n 1.5 20 0 (rad->deg d) 25.576 35.576) (- (len3 x y z) 40)))
        (body pinion (shape (bevel 15 (atan2 15 45))) (attr density 7.85g/cm3))
        (body wheel (shape (bevel 45 (atan2 45 15))) (attr density 7.85g/cm3))
        (assembly (pair (a 0))
          (instance p pinion (at) (turn z a) (effort-limits -2000N*mm 2000N*mm))
          (instance w wheel (at (rot-y 90)) (turn z (/ (* -15 a) 45)))
          (relate p w (bevel-mesh 15 45 1.5 10 20 0 0.98)))
        (assembly root (instance g (pair (* 360 (/ t 1s)))))
        (load drag root/g/w (torque -2000N*mm 0 0))"#;
    let of = values_at(
        &format!(
            r#"{pair}
            (require tooth (value (max (nth (tension root/g/p root/g/w) 0) (nth (tension root/g/p root/g/w) 1)) > -1e9N))
            (require idle (value (min (nth (tension root/g/p root/g/w) 0) (nth (tension root/g/p root/g/w) 1)) > -1e9N))
            (def rp (reaction root/g/p))
            (def rw (reaction root/g/w))
            (require px (value (nth rp 0) > -1e9N)) (require py (value (nth rp 1) > -1e9N)) (require pz (value (nth rp 2) > -1e9N))
            (require wx (value (nth rw 0) > -1e9N)) (require wz (value (nth rw 2) > -1e9N))"#
        ),
        0.0,
    );
    let (dp, dw) = (15f64.atan2(45.0), 45f64.atan2(15.0));
    let rm = 1.5 * 45.0 / (2.0 * dw.sin()) - 5.0;
    let ft = 2000.0 / (rm * dw.sin());
    let fr = ft * 20f64.to_radians().tan();
    let (ft, fr) = (0.98 * ft, 0.98 * fr);
    near(of("tooth"), ft / 20f64.to_radians().cos(), 1e-9, "the tooth force along the line of action");
    near(of("idle"), 0.0, 1e-9, "the idle flank");
    // The bearings' forces are known to what the gears' imbalance may be within their measures (their
    // centres read on their axes): 1.1e-5 N at the pinion and 1.7e-5 N at the wheel at this spin, and a
    // `>` claim's value is its low end, that much under.
    near(of("py"), ft, 2e-5, "the pinion's tangential force");
    near(of("px"), -fr * dp.cos(), 2e-5, "the pinion's radial force, from the wheel");
    near(of("pz"), fr * dp.sin(), 2e-5, "the pinion's axial force, away from its apex");
    near(of("wx"), fr * dw.cos(), 2e-5, "the wheel's radial force (its frame's x is the world's -z)");
    near(of("wz"), fr * dw.sin(), 2e-5, "the wheel's axial force, away from its apex");
    let e = claims_at_moment(pair, 0.25).into_iter().find(|x| x.0 == "effort-limits:root/g/p").unwrap();
    near(2000.0 - e.3, 2000.0 / 3.0 * 0.98, 1e-9, "driven back: the ideal 666.667 N mm times the efficiency");
}

/// A 20-tooth pinion of module 1 on a motor, its axis 10 mm above a rack's
/// pitch line, 400 N on the rack driving it back: by hand, the mesh passes
/// its efficiency, 0.98, of what the rack gives, so the motor gives 400 x
/// 10 x 0.98 = 3920 N mm, the tooth force is 392 / cos 20 and parts them
/// with 392 tan 20 = 142.676 N; the rack's friction, the other 8 N, pairs
/// with 8 N on the pinion's axis along the rack, so its bearings hold the
/// whole 400 N along it. The rack slid the other way is refused: its
/// coupling would do work.
#[test]
fn a_rack_and_pinions_force_by_hand() {
    let src = |dir: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body pinion (shape (subtract (extrude-z (involute-gear-2d 20 1 20 0) -4 4) (cylinder-z 2 20))) (attr density 7.85g/cm3))
            (body rack (shape (max (extrude-z (involute-rack-2d 1 20 0) -4 4) (slab-x -60 60) (- -4 y))) (attr density 7.85g/cm3))
            (assembly (pair (a 0))
              (instance p pinion (at (move 0 10 0)) (turn z a) (effort-limits -9000N*mm 9000N*mm))
              (instance r rack (at) (slide x (* {dir} (/ (* pi 20) 360) a)))
              (relate p r (rack-mesh 20 1 20 0 0.98)))
            (assembly root (instance g (pair (* 90 (/ t 1s)))))
            (load drag root/g/r (force 400N 0 0))
            (require tooth (value (max (nth (tension root/g/p root/g/r) 0) (nth (tension root/g/p root/g/r) 1)) > -1e9N))
            (def rp (reaction root/g/p))
            (require px (value (nth rp 0) > -1e9N)) (require py (value (nth rp 1) > -1e9N))"#
        )
    };
    let of = values_at(&src("1"), 0.0);
    near(of("tooth"), 392.0 / 20f64.to_radians().cos(), 1e-9, "the tooth force along the line of action");
    near(of("px"), -400.0, 1e-9, "the pinion's force along the rack: the tooth's and the friction's");
    near(of("py"), -392.0 * 20f64.to_radians().tan(), 1e-9, "the parting force, the pinion held from the rack");
    let e = claims_at_moment(&src("1"), 0.25).into_iter().find(|x| x.0 == "effort-limits:root/g/p").unwrap();
    near(9000.0 - e.3, 4000.0 * 0.98, 1e-9, "driven back: the ideal 4000 N mm times the efficiency");
    let wrong = claims_at_moment(&src("-1"), 0.0).into_iter().find(|x| x.0 == "tooth").unwrap();
    assert!(wrong.1 == Status::Unknown && wrong.2.contains("does not match the joints' motion"), "{wrong:?}");
}

/// A screw of lead 8 mm on a motor's shaft driving a nut on a vertical
/// slide, the nut held down by 200 N, a 51100 under the screw: by hand,
/// the screw's thrust is the 200 N, the 51100 and the motor's shaft carry
/// it, and the slide passes nothing (the screw does). Raising it the motor
/// gives 200 x 8 / (2 pi 0.4) = 636.620 N mm (Rexroth's Mta, efficiency
/// 0.4), its shaft passes that, and the thread's friction acts on the nut
/// as on the screw: the nut's guide holds the same 636.620 N mm, not the
/// ideal 254.648. Lowering it, the thrust is still the 200 N, the motor
/// and the guide its back efficiency 0.3 times the ideal, 76.394 (Mta').
/// A trapezoidal screw of efficiency 0.35 given alone back-drives at 2 - 1
/// / 0.35 = -0.857: it self-locks, and lowering the load the motor drives
/// it down with -218.270 N mm, the guide holding that the other way. The
/// screw on a joint of its own, driven through nothing a coupling declares,
/// and a lead of the other sign than the nut's slide reads, are refused,
/// saying so.
#[test]
fn a_lead_screws_thrust_by_hand() {
    let src = |lift: &str, screw: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (import "../tests/fixtures/dynamics/thrust-bearing.ink" as th)
            (body rotor (shape (cylinder-z 4 30)) (attr density 7.85g/cm3))
            (body screw (shape (move (cylinder-z 4 200) 0 0 100)) (attr density 7.85g/cm3))
            (body block (shape (box 40 40 20)) (attr density 2.7g/cm3))
            (assembly (lift (a 0))
              (instance rotor rotor (at (move 0 0 -40)) (turn z a) (effort-limits -2000N*mm 2000N*mm))
              {screw}
              (instance thrust th.bearing (at (move 0 0 -20)))
              (instance c block (at (move 0 0 100)) (slide z (* 8 (/ a 360))))
              (relate screw c (lead-screw 8 0.4 0.3)))
            (assembly root (instance l (lift {lift})))
            (load weight root/l/c (force 0 0 -200N))
            (require bearing (value (nth (nth (reaction root/l/rotor (by root/l/thrust)) 0) 2) > -1e9N))
            (require shaft (value (nth (reaction root/l/rotor) 2) > -1e9N))
            (require torque (value (nth (reaction root/l/rotor) 5) > -1e9N*mm))
            (require guide (value (nth (reaction root/l/c) 5) > -1e9N*mm))
            (require slide (value (nth (reaction root/l/c) 2) > -1e9N))
            (require thrust (value (nth (tension root/l/screw root/l/c) 0) > -1e9N))
            (require pull (value (nth (tension root/l/screw root/l/c) 1) > -1e9N))"#
        )
    };
    let on_rotor = "(instance screw screw (frame-of rotor))";
    let c = claims_at_moment(&src("(* 360 (/ t 1s))", on_rotor), 0.5);
    let of = |n: &str| c.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {c:?}")).3;
    let ideal = 200.0 * 8.0 / (2.0 * std::f64::consts::PI);
    near(of("bearing"), 200.0, 1e-9, "the thrust bearing's share");
    near(of("shaft"), 200.0, 1e-9, "the screw's thrust on the motor's shaft");
    near(of("torque"), ideal / 0.4, 1e-9, "the motor's shaft passes its drive torque");
    near(of("guide"), -ideal / 0.4, 1e-9, "the nut's guide holds the thread's torque, its friction with it");
    near(of("slide"), 0.0, 1e-9, "the slide passes nothing");
    near(of("thrust"), 200.0, 1e-9, "the screw's thrust, its push line's tension");
    near(of("pull"), 0.0, 1e-9, "its pull line slack");
    near(2000.0 - of("effort-limits:root/l/rotor"), ideal / 0.4, 1e-9, "raising: Rexroth's drive torque");
    let down = claims_at_moment(&src("(* -360 (/ t 1s))", on_rotor), 0.5);
    let at = |c: &[(String, Status, String, f64)], n: &str| c.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {c:?}")).3;
    near(2000.0 - at(&down, "effort-limits:root/l/rotor"), ideal * 0.3, 1e-9, "lowering: its back efficiency times the ideal");
    near(at(&down, "thrust"), 200.0, 1e-9, "lowering: the thrust is still the load");
    near(at(&down, "bearing"), 200.0, 1e-9, "and the thrust bearing carries it");
    near(at(&down, "guide"), -ideal * 0.3, 1e-9, "lowering: the guide holds what the screw passes back");
    let locking = claims_at_moment(&src("(* -360 (/ t 1s))", on_rotor).replace("(lead-screw 8 0.4 0.3)", "(lead-screw 8 0.35)"), 0.5);
    near(at(&locking, "effort-limits:root/l/rotor") - 2000.0, ideal * (2.0 - 1.0 / 0.35), 1e-9, "lowering a self-locking screw: the motor drives it down");
    near(at(&locking, "guide"), -ideal * (2.0 - 1.0 / 0.35), 1e-9, "and the guide holds that the other way");
    let signed = claims_at_moment(&src("(* 360 (/ t 1s))", on_rotor).replace("(slide z (* 8 (/ a 360)))", "(slide z (* -8 (/ a 360)))"), 0.5);
    let b = signed.iter().find(|x| x.0 == "guide").unwrap();
    assert!(b.1 == Status::Unknown && b.2.contains("separates its parts") && b.2.contains("does not match the joints' motion"), "{b:?}");
    let own = claims_at_moment(&src("(* 360 (/ t 1s))", "(instance screw screw (at (move 0 0 -40)) (turn z a))"), 0.5);
    let b = own.iter().find(|x| x.0 == "guide").unwrap();
    assert!(b.1 == Status::Unknown && b.2.contains("the balance of root/l/screw along its motion is not met by the couplings"), "{b:?}");
}

/// A screw of lead 4 mm along an MGN12H's rail, its nut 20 mm above the
/// carriage on the screw's axis, related as a coupling: the nut, named as a
/// support, carries the screw's thrust (the coupling's tension, along its
/// axis), and the carriage the rest, the thrust's line's pitch among it, as
/// the reaction left once the coupling is taken out says.
#[test]
fn a_coupled_screws_nut_carries_its_thrust() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (import "../tests/fixtures/dynamics/carriage.ink" as car)
        (def spec (get (read-json "../tests/fixtures/dynamics/ratings.json") "nut"))
        (body nut (shape (cylinder-z 14 40)) (attr density 7.85g/cm3) (frame load (at)))
        (body screw (shape (cylinder-z 6 600)) (attr density 7.85g/cm3))
        (body block (shape (box 50 50 50)) (attr mass 2000g))
        (assembly (axis (a 0))
          (instance rail car.rail)
          (instance c car.carriage (frame-of rail) (slide x (* 4 (/ a 360))))
          (instance block block (at (move 0 0 33) (frame-of c)))
          (instance n nut (at (rot-y 90) (move 0 0 20) (frame-of c)))
          (instance sc screw (at (rot-y 90) (move 0 0 20) (frame-of rail)) (turn z a) (effort-limits -1e5N*mm 1e5N*mm))
          (relate sc n (lead-screw 4 0.9)))
        (assembly root (instance x (axis (* 90 (+ 150 (* 50 (sin (* 360deg/s t))))))))
        (def s (reaction root/x/c (by root/x/c root/x/n)))
        (def r (reaction root/x/c (frame-of root/x/c)))
        (require thrust (value (nth (nth s 1) 2) > -1e9N))
        (require pitch (value (nth (nth s 0) 4) > -1e9N*mm))
        (require pull (value (max (nth (tension root/x/sc root/x/n) 0) (nth (tension root/x/sc root/x/n) 1)) > -1e9N))
        (require fx (value (nth r 0) > -1e9N))
        (require my (value (nth r 4) > -1e9N*mm))"#;
    // The readings as the model computes them (`turned`); a holding claim answers its band's worst edge.
    let read = claims_at_moment(&turned(src), 0.25);
    let of = |n: &str| read.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {read:?}")).3;
    // At 0.25 s the carriage decelerates at 50 (2 pi)^2 mm/s2 toward -x:
    // the screw pulls the 2 kg block and the carriage back.
    assert!(of("pull") > 4.0, "the thrust {}", of("pull"));
    near(of("thrust").abs(), of("pull"), 1e-9 * of("pull"), "the nut's share: the screw's thrust");
    near(of("fx"), 0.0, 1e-9 * of("pull"), "the slide passes nothing: the screw does");
    near(of("pitch"), of("my"), 1e-9 * of("my").abs(), "the guide's pitch: the reaction left, the thrust's line in it");
}

/// A coupling whose lines do work along the motion is refused: a crossed
/// belt written as an open one (the driven pulley turning the other way)
/// separates its spans' ends. Twin
/// belts in parallel are refused as sharing by their stiffness. A tension
/// resting on a balance whose coordinate no actuator powers is refused.
#[test]
fn couplings_the_motion_does_not_fix_are_refused() {
    let drive = |turn: &str, limits: &str, twin: bool| {
        let second = if twin {
            "(instance pa2 pa (at (move 0 0 20) (frame-of m))) (instance s2 pb (at (move 0 0 20) (frame-of s))) (relate pa2 s2 (belt-drive 10 20 100N 0.95))"
        } else {
            ""
        };
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body rotor (shape (cylinder-z 4 30)) (attr density 1g/cm3))
            (body pa (shape (cylinder-z 10 6)) (attr density 1g/cm3))
            (body pb (shape (cylinder-z 20 6)) (attr density 1g/cm3))
            (assembly (drive (a 0))
              (instance m rotor (at) (turn z a) (effort-limits {limits}))
              (instance pa pa (at (frame-of m)))
              (instance s pb (at (move 100 0 0)) (turn z (* {turn} a)))
              (relate pa s (belt-drive 10 20 100N 0.95)) {second})
            (assembly root (instance d (drive (* 360 (/ t 1s)))))
            (load drag root/d/s (torque 0 0 -1000N*mm))
            (require t (value (nth (tension root/d/pa root/d/s) 0) > 0))"#
        )
    };
    let t = |src: String| claims_at_moment(&src, 0.25).into_iter().find(|x| x.0 == "t").unwrap();
    let open = t(drive("0.5", "-2000N*mm 2000N*mm", false));
    assert!(open.1 == Status::Holds && (open.3 - 125.0).abs() < 1e-9, "the open belt, locked: 100 + 50 / 2: {open:?}");
    let crossed = t(drive("-0.5", "-2000N*mm 2000N*mm", false));
    assert!(crossed.1 == Status::Unknown && crossed.2.contains("separates its parts") && crossed.2.contains("crossed belt"), "{crossed:?}");
    let twin = t(drive("0.5", "-2000N*mm 2000N*mm", true));
    assert!(twin.1 == Status::Unknown && twin.2.contains("couplings in parallel") && twin.2.contains("stiffness"), "{twin:?}");
    // A fixed drive supplying this same torque, one rounding step away: the
    // tension is known to that step's share of its rounding, and holds.
    let required = model(&drive("0.5", "-2000N*mm 2000N*mm", false)).at(0.25).effort[0];
    assert!((required - 500.0 / 0.95).abs() < 1e-9);
    let fixed = required.next_up();
    let rounded = t(drive("0.5", &format!("{fixed}N*mm {fixed}N*mm"), false));
    assert!(rounded.1 == Status::Holds && (rounded.3 - 125.0).abs() < 1e-9, "{rounded:?}");
    let unpowered = t(drive("0.5", "0 0", false));
    assert!(unpowered.1 == Status::Unknown && unpowered.2.contains("different from its fixed rating of 0"), "{unpowered:?}");
}

/// Two belt drives side by side, sharing no coupling: one's motor has a
/// fixed rating of 0, so its motion is not supplied and its tension and
/// shaft are unknown, saying so; the other's are its own: 100 + 50 / 2 N
/// on its tight span (installed and locked: 1000 N mm over the 20 mm
/// wheel, 50 N of effective tension), its shaft holding the belt's pull
/// along the shafts' line, its spans' 125 + 75 N at their angle to it, sin
/// b = (20 - 10) / 100.
#[test]
fn a_tension_rests_on_its_own_trains_balances() {
    let src = r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
        (body rotor (shape (cylinder-z 4 30)) (attr density 1g/cm3))
        (body pa (shape (cylinder-z 10 6)) (attr density 1g/cm3))
        (body pb (shape (cylinder-z 20 6)) (attr density 1g/cm3))
        (assembly (drive (a 0))
          (instance m rotor (at) (turn z a) (effort-limits -2000N*mm 2000N*mm))
          (instance pa pa (at (frame-of m)))
          (instance s pb (at (move 100 0 0)) (turn z (* 0.5 a)))
          (relate pa s (belt-drive 10 20 100N 0.95)))
        (assembly (stuck (a 0))
          (instance m rotor (at) (turn z a) (effort-limits 0N*mm 0N*mm))
          (instance pa pa (at (frame-of m)))
          (instance s pb (at (move 100 0 0)) (turn z (* 0.5 a)))
          (relate pa s (belt-drive 10 20 100N 0.95)))
        (assembly root (instance d (drive (* 360 (/ t 1s)))) (instance k (stuck (* 360 (/ t 1s))) (at (move 0 500 0))))
        (load drag root/d/s (torque 0 0 -1000N*mm))
        (load drag2 root/k/s (torque 0 0 -1000N*mm))
        (require free (value (nth (tension root/d/pa root/d/s) 0) > 0))
        (require free-shaft (value (nth (reaction root/d/m (at)) 0) < 0))
        (require stuck (value (nth (tension root/k/pa root/k/s) 0) > 0))
        (require stuck-shaft (value (nth (reaction root/k/m (at)) 0) < 0))"#;
    let c = claims_at_moment(src, 0.25);
    let of = |n: &str| c.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {c:?}"));
    let free = of("free");
    assert!(free.1 == Status::Holds && (free.3 - 125.0).abs() < 1e-9, "{free:?}");
    let shaft = of("free-shaft");
    assert!(shaft.1 == Status::Holds && (shaft.3 + 200.0 * (1.0f64 - 0.01).sqrt()).abs() < 1e-9, "{shaft:?}");
    for n in ["stuck", "stuck-shaft"] {
        let x = of(n);
        assert!(x.1 == Status::Unknown && x.2.contains("root/k/m") && x.2.contains("fixed rating of 0"), "{x:?}");
    }
}

/// What moves with a balancing joint's coordinate through nothing a joint
/// or coupling accounts for is named, and the tensions say they leave it
/// out: P3's four-bar with its crank belt-driven from a motor
/// its coupler and rocker pushing
/// on the crank's pin.
#[test]
fn a_tension_names_what_its_balances_leave_out() {
    let m = reacting(
        r#"(import "std") (gravity 0 -9.81m/s2 0) (t-span 0 1s)
        (def a 100) (def b 400) (def c 300) (def o4 350)
        (body crank (shape (move (box 100 4 4) 50 0 0)) (attr density 1g/cm3))
        (body coupler (shape (move (box 400 4 4) 200 0 0)) (attr density 1g/cm3))
        (body rocker (shape (move (box 300 4 4) 150 0 0)) (attr density 1g/cm3))
        (body rotor (shape (cylinder-z 4 10)) (attr density 1g/cm3))
        (body pul (shape (cylinder-z 10 4)) (attr density 1g/cm3))
        (fn ax (p) (* a (cos (deg->rad p))))
        (fn ay (p) (* a (sin (deg->rad p))))
        (fn ee (p) (sqrt (+ (square (- o4 (ax p))) (square (ay p)))))
        (fn th3 (p) (+ (atan2 (neg (ay p)) (- o4 (ax p))) (acos (/ (- (+ (* b b) (square (ee p))) (* c c)) (* 2 b (ee p))))))
        (fn bx (p) (+ (ax p) (* b (cos (th3 p)))))
        (fn by (p) (+ (ay p) (* b (sin (th3 p)))))
        (assembly (fourbar (p 0))
          (instance m rotor (at (move -200 0 0)) (turn z p) (effort-limits -1e6N*mm 1e6N*mm))
          (instance mp pul (at (frame-of m)))
          (instance crank crank (at) (turn z p))
          (instance cp pul (at (frame-of crank)))
          (instance coupler coupler (at (rot-z (rad->deg (th3 p))) (move (ax p) (ay p) 0)))
          (instance rocker rocker (at (rot-z (rad->deg (atan2 (by p) (- (bx p) o4)))) (move o4 0 0)))
          (relate mp cp (belt-drive 10 10 50N 0.98)))
        (assembly root (instance f (fourbar (* (rad->deg 10) (/ t 1s)))))"#,
        &["root/f/crank"],
        true,
    )
    .unwrap_or_else(|e| panic!("{}", e.message));
    assert_eq!(m.unseparated, ["root/f/coupler", "root/f/rocker"], "what the crank's balance leaves out");
    assert_eq!(m.reactions[0].moves_with, ["root/f/coupler", "root/f/rocker"], "and the crank's reaction");
}

/// The no-work check measures a line's work against each end's motion
/// before it cancels: a planet whose contact with its sun is its instant centre,
/// on a turning carrier and on a fixed sun, carries 500 N mm over its 20 mm
/// pitch radius along its line of action, 25 / cos 20 = 26.604 N; where
/// the drag turns the planet on its carrier (the sun turning faster than
/// the carrier), power enters the mesh at the planet in the carrier's
/// frame, and the tooth passes the mesh's 0.98 of it, 26.072 N; a pulley
/// rolling on a belt span fixed at its other end carries 200 N mm over 10
/// mm, 20 N; radii at a catalog's rounding (20.001) pass.
#[test]
fn couplings_resting_at_an_instant_centre_answer() {
    let planet = |motion: &str, sun: &str, planet: &str| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (fn gear (z) (subtract (extrude-z (involute-gear-2d z 2 20 0) -5 5) (cylinder-z 4 20)))
            (body sun (shape (gear 20)) (attr density 7.85g/cm3))
            (body planet (shape (gear 20)) (attr density 7.85g/cm3))
            (body arm (shape (move (box 60 8 4) 20 0 -10)) (attr density 2.7g/cm3))
            (assembly (pg (a 0) (c 0))
              (instance s sun (at) {sun})
              (instance k arm (at) (turn z c) (effort-limits -5000N*mm 5000N*mm))
              (instance p planet (at (move 40 0 0) (frame-of k)) (turn z (mesh-turn 20 20 0 {planet})))
              (relate s p (gear-mesh 20 20 2 20 0 0.98)))
            (assembly root (instance g (pg {motion})))
            (load drag root/g/p (torque 0 0 -500N*mm))
            (require tooth (value (max (nth (tension root/g/s root/g/p) 0) (nth (tension root/g/s root/g/p) 1)) > -1e9N))"#
        )
    };
    let tooth = |src: String| claims_at_moment(&src, 0.25).into_iter().find(|x| x.0 == "tooth").unwrap();
    let fn_ = 25.0 / 20f64.to_radians().cos();
    // In the carrier's frame (Willis's): with the sun turning faster than the carrier the drag drives the
    // planet there, power enters the mesh at the planet, and the tooth passes the mesh's 0.98 of it; with
    // the sun held the sun drives it, and the tooth carries the planet's whole 26.604 N.
    for (motion, sun, turn, passed) in
        [("(* 360 (/ t 1s)) (* 90 (/ t 1s))", "(turn z a) (effort-limits -5000N*mm 5000N*mm)", "(- a c)", 0.98), ("0 (* 90 (/ t 1s))", "", "(neg c)", 1.0)]
    {
        let t = tooth(planet(motion, sun, turn));
        assert!(t.1 == Status::Holds && (t.3 - passed * fn_).abs() < 1e-9, "{motion}: {t:?}");
    }
    let roll_at = |rb: &str, eff: &str| {
        claims_at_moment(
            &format!(
                r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
                (body cart (shape (box 60 40 10)) (attr density 1g/cm3))
                (body pul (shape (cylinder-z 10 6)) (attr density 1g/cm3))
                (body clamp (shape (move (box 10 10 6) -205 -10 0)) (attr density 1g/cm3))
                (assembly (axis (s 0))
                  (instance clamp clamp (at))
                  (instance cart cart (at (move 0 0 -20)) (slide x s) (effort-limits -1e4N 1e4N))
                  (instance idler pul (at (move 0 0 20) (frame-of cart)) (turn z (rad->deg (/ s (neg {rb})))))
                  (relate clamp idler (belt-span (vec (vec -200 -10 0) (vec s -10 0)) 50N {eff})))
                (assembly root (instance a (axis (* 100 (/ t 1s)))))
                (load drag root/a/idler (torque 0 0 200N*mm))
                (require span (value (nth (tension root/a/clamp root/a/idler) 0) > -1e9N))"#
            ),
            0.25,
        )
        .into_iter()
        .find(|x| x.0 == "span")
        .unwrap()
    };
    // (The span that loses nothing: with 0.95 its loss is on the pulley's rolling on the carriage, which is
    // not found, and its tension is not settled, saying so.)
    let roll = |rb: &str| roll_at(rb, "1");
    let lossy = roll_at("10", "0.95");
    assert!(lossy.1 == Status::Unknown && lossy.2.contains("rides a joint that moves"), "{lossy:?}");
    let r = roll("10");
    assert!(r.1 == Status::Holds && (r.3 - 20.0).abs() < 1e-9, "{r:?}");
    assert_eq!(roll("10.001").1, Status::Holds, "a radius at a catalog's rounding passes");
    let wrong = roll("-10");
    assert!(wrong.1 == Status::Unknown && wrong.2.contains("separates its parts at -2.000000"), "{wrong:?}");
}

/// A planet on a carrier turning at 90 deg/s round a sun held still, a
/// 500 N mm drag on the planet, the mesh 0.98: by hand, in the carrier's
/// frame (Willis's) the sun turns at -90 deg/s and gives the mesh 500 x
/// pi / 2 N mm/s, of which it loses 1 / 0.98 - 1, so the carrier gives
/// the ideal 1000 N mm (the drag's 500 x pi over its pi / 2) and that loss
/// over its pi / 2, 10.204 N mm: 1010.204. Charged on each end's own joint
/// it was charged nowhere (the sun has none), and the carrier read the
/// ideal 1000; then it was said not found. In the carrier's frame it is
/// found: the sun turns back at the carrier's rate there.
#[test]
fn a_planets_loss_is_charged_in_its_carriers_frame() {
    let src = r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
        (fn gear (z) (subtract (extrude-z (involute-gear-2d z 2 20 0) -5 5) (cylinder-z 4 20)))
        (body sun (shape (gear 20)) (attr density 7.85g/cm3))
        (body planet (shape (gear 20)) (attr density 1e-6g/mm3))
        (body arm (shape (move (box 60 8 4) 20 0 -10)) (attr density 1e-6g/mm3))
        (assembly (pg (c 0))
          (instance s sun (at))
          (instance k arm (at) (turn z c) (effort-limits -5000N*mm 5000N*mm))
          (instance p planet (at (move 40 0 0) (frame-of k)) (turn z (mesh-turn 20 20 0 (neg c))))
          (relate s p (gear-mesh 20 20 2 20 0 0.98)))
        (assembly root (instance g (pg (* 90 (/ t 1s)))))
        (load drag root/g/p (torque 0 0 -500N*mm))"#;
    let e = claims_at_moment(src, 0.25).into_iter().find(|x| x.0 == "effort-limits:root/g/k").unwrap();
    near(5000.0 - e.3, 1000.0 + 500.0 * (1.0 / 0.98 - 1.0), 1e-6, "the carrier's effort, the mesh's loss in the carrier's frame among it");
}

/// A screw (E 0.35, so 1 / E - 1 = 1.857) lifting a nut held down by 200 N
/// while a deck held down by 1000 N slides with it through nothing a
/// coupling declares: the loss charged rests on
/// the screw's 200 N, and the deck's 1000 x 8 / (2 pi) = 1273.24 N mm could
/// pass through the screw too, 2364.6 N mm more (4365.4 in all). Against
/// +-5000 its margin, 2999.2, passes that: holds, saying so; against +-4000,
/// 1999.2 does not: unknown.
#[test]
fn an_effort_whose_loss_leaves_bodies_out_is_decided_past_them() {
    let src = |lim: f64| {
        format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body rotor (shape (cylinder-z 4 30)) (attr density 7.85g/cm3))
            (body screw (shape (move (cylinder-z 4 200) 0 0 100)) (attr density 7.85g/cm3))
            (body block (shape (box 40 40 20)) (attr density 2.7g/cm3))
            (body deck (shape (box 100 100 5)) (attr density 2.7g/cm3))
            (assembly (lift (a 0))
              (instance rotor rotor (at (move 0 0 -40)) (turn z a) (effort-limits -{lim}N*mm {lim}N*mm))
              (instance screw screw (frame-of rotor))
              (instance c block (at (move 0 0 100)) (slide z (* 8 (/ a 360))))
              (instance deck deck (at (move 200 0 100)) (slide z (* 8 (/ a 360))))
              (relate screw c (lead-screw 8 0.35)))
            (assembly root (instance l (lift (* 360 (/ t 1s)))))
            (load weight root/l/c (force 0 0 -200N))
            (load payload root/l/deck (force 0 0 -1000N))"#
        )
    };
    let effort = |lim: f64| claims_at_moment(&src(lim), 0.25).into_iter().find(|x| x.0 == "effort-limits:root/l/rotor").unwrap();
    let ideal = |f: f64| f * 8.0 / (2.0 * std::f64::consts::PI);
    let held = effort(5000.0);
    near(held.3, 5000.0 - (ideal(200.0) / 0.35 + ideal(1000.0)), 1e-6, "the effort, the deck's ideal");
    assert_eq!(held.1, Status::Holds);
    let close = effort(4000.0);
    assert!(close.1 == Status::Unknown && close.2.contains("root/l/deck") && close.2.contains(&format!("{:.6}", ideal(1000.0) * (1.0 / 0.35 - 1.0))), "{close:?}");
}

/// A 2 mm lead screw on a passive turn, its nut on a slide a fixed 200 N
/// drives: the tensions come from the screw's turn balance, whose rounding
/// is a moment (N mm). Each line's tension's rounding is that over the
/// line's arm, 2 / 2 pi = 0.318 mm, in newtons; it was the moment taken as
/// newtons, a third of what the balance's rounding allows.
#[test]
fn a_tensions_rounding_is_its_balances_over_its_arm() {
    let src = r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body rotor (shape (cylinder-z 4 30)) (attr density 7.85g/cm3))
            (body screw (shape (move (cylinder-z 4 200) 0 0 100)) (attr density 7.85g/cm3))
            (body block (shape (box 40 40 20)) (attr density 2.7g/cm3))
            (assembly (lift (a 0))
              (instance rotor rotor (at (move 0 0 -40)) (turn z a))
              (instance screw screw (frame-of rotor))
              (instance c block (at (move 0 0 100)) (slide z (* 2 (/ a 360))) (effort-limits 200N 200N))
              (relate screw c (lead-screw 2 1)))
            (assembly root (instance l (lift (* 360 (/ t 1s)))))
            (load weight root/l/c (force 0 0 -200.0000001N))"#;
    let m = model(src);
    let at = m.at(0.5);
    let [(r, true)] = m.balancing()[..] else { panic!("one turn balance: {:?}", m.balancing()) };
    let (force, moment) = at.reactions[r].rounding();
    let o = &at.reactions[r].frame[9..12];
    let arm = 2.0 / (2.0 * std::f64::consts::PI);
    let want = (moment + (o[0] * o[0] + o[1] * o[1] + o[2] * o[2]).sqrt() * force) / arm;
    assert!(want > 0.0 && at.tension_noise.iter().all(|n| (n / want - 1.0).abs() < 1e-12), "{:?} against {want}", at.tension_noise);
    // The drive short of its rating within its rounding: the tensions' band is
    // one share of that, each line's (and of the tensions' own, here none).
    let k = m.coupling_index("root/l/screw", "root/l/c").expect("the screw's coupling");
    let band = m.tension_shortfall(&at, k).expect("a fixed drive short within rounding");
    let share = band[0] / want;
    assert!(share > 0.0 && band.iter().all(|b| (b / (share * want) - 1.0).abs() < 1e-12), "{band:?} against {want}");
}
