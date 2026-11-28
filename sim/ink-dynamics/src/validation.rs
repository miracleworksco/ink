//! The reference problems inverse dynamics on placements can answer
//! (P1, P2, P3, P5, P6, P8, P9, P10, P16, P17, P18; reactions: P1's pin
//! force, P8's bearing moment,
//! P9's bearing couple and its bearings' forces, P4's belt force; loads:
//! P10's tip wrench), each written as an ink
//! model whose bodies are
//! given the problem's exact mass properties (Model::set_mass; a rotor's
//! by `spin-inertia`; P16's belts at the masses ink measures), and held to
//! its tier: 1e-6 of the column's peak, or the table's printed
//! digits where coarser.

use super::*;
use ink_lang::Doc;
use std::path::Path;

static REGISTRY: ink_lang::Registry = ink_lang::Registry::new(&[&Dynamics]);

fn model(src: &str) -> Model {
    reacting(src, &[])
}

/// The model of `src` with the reactions of `paths`.
fn reacting(src: &str, paths: &[&str]) -> Model {
    let doc = Doc::load_str_in(&REGISTRY, src, "validation.ink", Path::new("."), &[]).unwrap_or_else(|e| panic!("{} at {}", e.message, e.line));
    let g = gravity(&doc).expect("a gravity");
    let paths: Vec<String> = paths.iter().map(|p| p.to_string()).collect();
    build(&doc.reloader().unwrap(), "root", g, 0.0, &paths, false).unwrap_or_else(|e| panic!("{} ({:?})", e.message, e.hint))
}

fn diag(a: f64, b: f64, c: f64) -> [[f64; 3]; 3] {
    [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]]
}

const O: [f64; 3] = [0.0; 3];

fn near(got: f64, want: f64, tol: f64, what: &str) {
    assert!((got - want).abs() <= tol, "{what}: {got} against {want} (within {tol})");
}

/// P1: a bar pinned at one end, swinging 60 degrees either way: the
/// actuator's torque, and the pin's force on the bar (its reaction, in the
/// world's axes; about its joint's axis, the torque).
#[test]
fn p1_compound_pendulum() {
    let mut m = reacting(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -1000N*mm 1000N*mm)))
        (assembly root (instance p (pendulum (* 60 (sin (* 180deg/s t))))))"#,
        &["root/p/bar"],
    );
    let (mass, l, w, h) = (162.0, 300.0, 20.0, 10.0);
    m.set_mass("root/p/bar", mass, [0.0, l / 2.0, 0.0], diag(mass * (l * l + h * h) / 12.0, mass * (w * w + h * h) / 12.0, mass * (l * l + w * w) / 12.0));
    let io = 4.865400e-3;
    let world = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    let pins = [(0.0, 0.0, 1.852225), (0.25, 0.219804, 1.566478), (0.5, 0.125575, 1.371717), (1.25, -0.219804, 1.566478), (1.5, -0.125575, 1.371717)];
    for ((t, table), (_, px, pz)) in [(0.0, 0.0), (0.25, 125.266), (0.5, 156.160), (1.25, -125.266), (1.5, -156.160)].into_iter().zip(pins) {
        let pi = std::f64::consts::PI;
        let (th, thdd) = (pi / 3.0 * (pi * t).sin(), -pi / 3.0 * pi * pi * (pi * t).sin());
        let tau = 1000.0 * (io * thdd + 0.162 * 9.81 * 0.15 * th.sin());
        let at = m.at(t);
        let e = at.effort[0];
        near(e, tau, 1e-6 * 156.16, "P1 against its formula");
        near(e, table, 5e-4 + 1e-6 * 156.16, "P1 against its table");
        let (f, _) = at.reactions[0].in_frame(&world);
        near(f[0], px, 5e-7 + 1e-6 * 1.85, &format!("P1 pin force x at {t} s"));
        near(f[1], 0.0, 1e-12, &format!("P1 pin force y at {t} s"));
        near(f[2], pz, 5e-7 + 1e-6 * 1.85, &format!("P1 pin force z at {t} s"));
        // About its joint's axis the rest exerts the actuator's torque.
        near(at.reactions[0].local().1[2], e, 1e-9 * 156.16, "P1 reaction about the joint's axis");
    }
}

/// P2: a slider-crank at 600 rpm, its rod and piston placed by the crank
/// angle in closed form (the loop closed in the placements).
#[test]
fn p2_slider_crank() {
    let mut m = model(
        r#"(import "std") (gravity 0 -9.81m/s2 0)
        (def r 50) (def l 150)
        (body crank (shape (move (box 50 4 4) 25 0 0)) (attr density 1g/cm3))
        (body rod (shape (move (box 150 4 4) 75 0 0)) (attr density 1g/cm3))
        (body piston (shape (box 10 10 10)) (attr density 1g/cm3))
        (assembly (slider (phi 0))
          (instance crank crank (at) (turn z phi) (effort-limits -1e6N*mm 1e6N*mm))
          (instance rod rod (at (rot-z (rad->deg (asin (/ (* (neg r) (sin (deg->rad phi))) l)))) (move (* r (cos (deg->rad phi))) (* r (sin (deg->rad phi))) 0)))
          (instance piston piston (at (move (+ (* r (cos (deg->rad phi))) (sqrt (- (* l l) (* r r (square (sin (deg->rad phi))))))) 0 0))))
        (assembly root (instance s (slider (* 3600 (/ t 1s)))))"#,
    );
    let bar = |mass: f64, len: f64| diag(0.0, mass * len * len / 12.0, mass * len * len / 12.0);
    m.set_mass("root/s/crank", 500.0, [25.0, 0.0, 0.0], bar(500.0, 50.0));
    m.set_mass("root/s/rod", 300.0, [75.0, 0.0, 0.0], bar(300.0, 150.0));
    m.set_mass("root/s/piston", 400.0, O, diag(0.0, 0.0, 0.0));
    let table = [196.200, 4028.046, 2166.705, -1919.188, -2910.655, -1694.855, -196.200, 1355.026, 2714.455, 1919.188, -1970.505, -3688.218];
    for (k, tau) in table.into_iter().enumerate() {
        near(m.at(30.0 * k as f64 / 3600.0).effort[0], tau, 5e-4 + 1e-6 * 4028.0, &format!("P2 at {} deg", 30 * k));
    }
}

/// P3: a crank-rocker four-bar at 10 rad/s, its coupler and rocker placed
/// by the open branch's closed form.
#[test]
fn p3_four_bar() {
    let mut m = model(
        r#"(import "std") (gravity 0 -9.81m/s2 0)
        (def a 100) (def b 400) (def c 300) (def o4 350)
        (body crank (shape (move (box 100 4 4) 50 0 0)) (attr density 1g/cm3))
        (body coupler (shape (move (box 400 4 4) 200 0 0)) (attr density 1g/cm3))
        (body rocker (shape (move (box 300 4 4) 150 0 0)) (attr density 1g/cm3))
        (fn ax (p) (* a (cos (deg->rad p))))
        (fn ay (p) (* a (sin (deg->rad p))))
        (fn ee (p) (sqrt (+ (square (- o4 (ax p))) (square (ay p)))))
        (fn th3 (p) (+ (atan2 (neg (ay p)) (- o4 (ax p))) (acos (/ (- (+ (* b b) (square (ee p))) (* c c)) (* 2 b (ee p))))))
        (fn bx (p) (+ (ax p) (* b (cos (th3 p)))))
        (fn by (p) (+ (ay p) (* b (sin (th3 p)))))
        (assembly (fourbar (p 0))
          (instance crank crank (at) (turn z p) (effort-limits -1e6N*mm 1e6N*mm))
          (instance coupler coupler (at (rot-z (rad->deg (th3 p))) (move (ax p) (ay p) 0)))
          (instance rocker rocker (at (rot-z (rad->deg (atan2 (by p) (- (bx p) o4)))) (move o4 0 0))))
        (assembly root (instance f (fourbar (* (rad->deg 10) (/ t 1s)))))"#,
    );
    let bar = |mass: f64, len: f64| diag(0.0, mass * len * len / 12.0, mass * len * len / 12.0);
    m.set_mass("root/f/crank", 200.0, [50.0, 0.0, 0.0], bar(200.0, 100.0));
    m.set_mass("root/f/coupler", 800.0, [200.0, 0.0, 0.0], bar(800.0, 400.0));
    m.set_mass("root/f/rocker", 600.0, [150.0, 0.0, 0.0], bar(600.0, 300.0));
    let table = [-596.597, 758.051, 224.960, -739.629, -1001.557, -441.039, 594.408, 1154.747];
    for (k, tau) in table.into_iter().enumerate() {
        let t = (45.0 * k as f64).to_radians() / 10.0;
        near(m.at(t).effort[0], tau, 5e-4 + 1e-6 * 1154.7, &format!("P3 at {} deg", 45 * k));
    }
}

/// P5: a lead screw's motor (68 g cm^2, its `spin-inertia`) lifting 3 kg
/// two turns on a cycloidal profile, ideal.
#[test]
fn p5_lead_screw() {
    let mut m = model(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body rotor (shape (cylinder-z 3 10)) (attr density 7.8g/cm3) (attr spin-inertia 6800g*mm2))
        (body load (shape (box 20 20 20)) (attr density 1g/cm3))
        (def H 16) (def T 0.4)
        (assembly (lift (z 0))
          (instance rotor rotor (at) (turn z (* 45 z)) (effort-limits -1000N*mm 1000N*mm))
          (instance load load (at (move 50 0 0)) (slide z z)))
        (assembly root (instance l (lift (* H (- (/ (/ t 1s) T) (/ (sin (/ (* 2 pi (/ t 1s)) T)) (* 2 pi)))))))"#,
    );
    m.set_mass("root/l/load", 3000.0, O, diag(0.0, 0.0, 0.0));
    assert_eq!(m.spin_inertia, ["root/l/rotor"]);
    for (t, tau) in [(0.0, 37.471440), (0.05, 41.541310), (0.1, 43.227105), (0.2, 37.471440), (0.3, 31.715774), (0.35, 33.401570)] {
        near(m.at(t).effort[0], tau, 5e-7 + 1e-6 * 43.2, &format!("P5 at {t} s"));
    }
}

/// P6: a 50:1 geared arm with a tip payload, swung from hanging to
/// upright: the motor's torque, its rotor's 54 g cm^2 reflected.
#[test]
fn p6_geared_arm() {
    let mut m = model(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body arm (shape (move (box 400 30 10) 200 0 0)) (attr density 2.7g/cm3))
        (body payload (shape (sphere 5)) (attr density 7.8g/cm3))
        (body rotor (shape (cylinder-z 3 10)) (attr density 7.8g/cm3) (attr spin-inertia 5400g*mm2))
        (assembly link (instance arm arm) (instance payload payload (at (move 400 0 0))))
        (assembly (geared (theta 0))
          (instance l link (at (rot-x 90)) (turn z theta))
          (instance rotor rotor (at (rot-x 90) (move 0 -100 0)) (turn z (* 50 theta)) (effort-limits -1000N*mm 1000N*mm)))
        (assembly root (instance g (geared (+ -90 (* 180 (- (/ t 1s) (/ (sin (* 360deg/s t)) (* 2 pi))))))))"#,
    );
    let (mass, l, w, h) = (324.0, 400.0, 30.0, 10.0);
    m.set_mass("root/g/l/arm", mass, [l / 2.0, 0.0, 0.0], diag(mass * (w * w + h * h) / 12.0, mass * (l * l + h * h) / 12.0, mass * (l * l + w * w) / 12.0));
    m.set_mass("root/g/l/payload", 500.0, O, diag(0.0, 0.0, 0.0));
    for (t, tau) in [(0.1, 26.764807), (0.25, 58.370822), (0.5, 51.953760), (0.75, -29.116747), (0.9, -24.659095)] {
        near(m.at(t).effort[0], tau, 5e-7 + 1e-6 * 69.1, &format!("P6 at {t} s"));
    }
    // Over the swing, the peak the steps find where the effort turns back:
    // 69.0853703 N mm at 0.3721 s, the greatest of 20001 moments over
    // 0.35..0.4 s to 1e-11.
    let s = m.steps((0.0, 1.0));
    assert!(s.unchecked.is_empty() && s.loose.is_empty());
    let peak = s.ends.windows(2).map(|w| m.at(ink_core::eval::span_mid((w[0], w[1]))).effort[0].abs()).fold(0.0, f64::max);
    assert!((peak - 69.0853703).abs() < 1e-6, "P6 peak {peak} over {} steps", s.ends.len() - 1);
}

/// P8: a wheel spinning on its axle while it steers: no effort along
/// either at constant rates; spinning up the steer, its transverse inertia
/// alone. The gyroscopic moment I_axle w_s W (z x axle) is the axle
/// bearings', the wheel's reaction.
#[test]
fn p8_steered_wheel() {
    let wheel = |steer: &str| {
        let mut m = reacting(
            &format!(
                r#"(import "std") (gravity 0 0 -9.81m/s2)
            (body wheel (shape (cylinder-z 50.8 38.1)) (attr density 1g/cm3))
            (assembly (spinning (spin 0)) (instance wheel wheel (at (rot-y 90)) (turn z spin) (effort-limits -1e6N*mm 1e6N*mm)))
            (assembly (module (steer 0) (spin 0)) (instance s (spinning spin) (at) (turn z steer) (effort-limits -1e6N*mm 1e6N*mm)))
            (assembly root (instance m (module {steer} (rad->deg (* (/ 4 0.0508) (/ t 1s))))))"#
            ),
            &["root/m/s/wheel"],
        );
        let (mass, r, w) = (200.0, 50.8, 38.1);
        let (it, ia) = (mass * (3.0 * r * r + w * w) / 12.0, 0.5 * mass * r * r);
        m.set_mass("root/m/s/wheel", mass, O, diag(it, it, ia));
        m
    };
    let steady = wheel("(* 360 (/ t 1s))");
    let world = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    for t in [0.0, 0.125, 0.3] {
        let at = steady.at(t);
        let e = &at.effort;
        assert!(e.iter().all(|x| x.abs() < 1e-9 * 127.67), "P8 at constant rates, t {t}: {e:?}");
        // The table's two rows, and the same moment at 0.3 s.
        let a = 2.0 * std::f64::consts::PI * t;
        let (f, mo) = at.reactions[0].in_frame(&world);
        let want = [-127.6743 * a.sin(), 127.6743 * a.cos(), 0.0];
        for k in 0..3 {
            near(mo[k], want[k], 5e-4 + 1e-6 * 127.67, &format!("P8 bearing moment [{k}] at {t} s"));
            near(f[k], [0.0, 0.0, 0.2 * 9.81][k], 1e-9, &format!("P8 bearing force [{k}] at {t} s"));
        }
    }
    near(steady.at(0.125).reactions[0].in_frame(&world).1[0], -90.279, 5e-4 + 1e-6 * 127.67, "P8 bearing moment x at 0.125 s, the table's");
    let spun = wheel("(rad->deg (* 10 (/ t 1s) (/ t 1s)))");
    let paths: Vec<&str> = spun.actuators.iter().map(|a| a.path.as_str()).collect();
    let (steer, spin) = (paths.iter().position(|p| *p == "root/m/s").unwrap(), paths.iter().position(|p| *p == "root/m/s/wheel").unwrap());
    let at = spun.at(0.5);
    let e = &at.effort;
    near(e[steer], 3.064510, 5e-7 + 1e-6 * 3.06, "P8 steer spin-up");
    near(e[spin], 0.0, 1e-9 * 203.2, "P8 spin torque while steering up");
    // Steered 2.5 rad: the gyroscopic 203.2 N mm along z x axle, and the
    // steer's spin-up about z.
    let mo = at.reactions[0].in_frame(&world).1;
    near(mo[0], -203.200 * 2.5f64.sin(), 5e-4 + 1e-6 * 203.2, "P8 spin-up's gyroscopic moment x");
    near(mo[1], 203.200 * 2.5f64.cos(), 5e-4 + 1e-6 * 203.2, "P8 spin-up's gyroscopic moment y");
    near(mo[2], 3.064510, 5e-7 + 1e-6 * 203.2, "P8 spin-up's moment about the steer axis");
}

/// P9: a plate turning about its diagonal, tilted in the world so its
/// inertia there has products: the drive torque along the diagonal, and
/// the bearing couple across it (its reaction in its joint's frame, whose z
/// is the diagonal).
#[test]
fn p9_plate_on_its_diagonal() {
    let plate = |angle: &str| {
        let mut m = reacting(
            &format!(
                r#"(import "std") (gravity 0 0 0)
            (body plate (shape (box 300 200 5)) (attr density 2.7g/cm3)
              (frame diag (axes (origin 0 0 0) (x 0 0 1) (z -300 -200 0))))
            (assembly (spinner (a 0)) (instance plate plate (at (rot-x 30) (rot-z 20)) diag (turn z a) (effort-limits -1e6N*mm 1e6N*mm)))
            (assembly root (instance s (spinner (rad->deg {angle}))))"#
            ),
            &["root/s/plate"],
        );
        m.set_mass("root/s/plate", 810.0, O, diag(2701688.0, 6076688.0, 8775000.0));
        m
    };
    let iu = (90000.0 * 2.701688e-3 + 40000.0 * 6.076688e-3) / 130000.0;
    for ((angle, alpha), couple) in
        [("(* 25 (/ t 1s) (/ t 1s))", 50.0), ("(* 100 (/ t 1s))", 0.0), ("(+ (* 100 (/ t 1s)) (* 25 (/ t 1s) (/ t 1s)))", 50.0)].into_iter().zip([77.885, 15576.923, 15577.118])
    {
        let at = plate(angle).at(0.0);
        let e = at.effort[0];
        near(e, 1000.0 * iu * alpha, 1e-6 * 187.007 + 1e-9 * 15576.9, &format!("P9 {angle}"));
        near(e, 1000.0 * [0.0, 0.187007][usize::from(alpha > 0.0)], 5e-4, &format!("P9 {angle} against its table"));
        let (f, mo) = at.reactions[0].local();
        near(mo[2], e, 1e-9 * 15577.1, &format!("P9 {angle}: the reaction about the diagonal is the drive's"));
        near(mo[0].hypot(mo[1]), couple, 5e-4 + 1e-6 * 15577.1, &format!("P9 {angle}: the bearing couple"));
        assert!(f.iter().all(|x| x.abs() < 1e-9), "P9 {angle}: no force on a plate turning about its centre without gravity: {f:?}");
        // Split between radial bearings 500 mm apart on the diagonal: the
        // couple over their spacing, equal and opposite, turning with it.
        let w = &at.reactions[0];
        let on = |d: f64| -> [f64; 12] { std::array::from_fn(|k| if k < 9 { w.frame[k] } else { w.frame[k] + d * w.frame[3 * (k - 9) + 2] }) };
        let part = |path, load| split::Part { path, carries: split::Carries::Radial, pulls: false, load, given: load, coupled: None, rating: "" };
        let supports = [part("b1", on(-250.0)), part("b2", on(250.0))];
        let (shares, _) = split::split(w, &supports, Some(split::Drive { turn: true, axis: 2 })).unwrap_or_else(|e| panic!("P9 {angle}: {e}"));
        let force = [31.153846, 0.155769, 31.154236][usize::from(alpha > 0.0) + usize::from(angle.contains('+'))];
        for (k, s) in shares.iter().enumerate() {
            near(s[0].hypot(s[1]), force, 5e-7 + 1e-6 * 31.15, &format!("P9 {angle}: bearing {k}'s force"));
            assert!(s[2..].iter().all(|x| x.abs() < 1e-12 * 31.15), "P9 {angle}: a radial bearing carries force across its axis alone: {s:?}");
        }
        near(shares[0][0] + shares[1][0], 0.0, 1e-9 * 31.15, &format!("P9 {angle}: opposite forces"));
        near(shares[0][1] + shares[1][1], 0.0, 1e-9 * 31.15, &format!("P9 {angle}: opposite forces"));
    }
}

/// P4's belt force by inverse dynamics: its carriage (1.5 kg on a level,
/// frictionless rail) moved as the motor's line drives it, x = v_f (t -
/// tau (1 - e^(-t / tau))): what pulls it along the rail, its reaction
/// along its slide, is m x'', the table's belt force.
#[test]
fn p4_belt_force() {
    let mut m = reacting(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (def vf 660mm/s) (def tau 0.012853596s)
        (body carriage (shape (box 40 40 20)) (attr density 1g/cm3))
        (assembly (axis (x 0)) (instance c carriage (at) (slide x x) (effort-limits -1000N 1000N)))
        (assembly root (instance a (axis (* vf (- t (* tau (- 1 (exp (neg (/ t tau))))))))))"#,
        &["root/a/c"],
    );
    m.set_mass("root/a/c", 1500.0, O, diag(1.0, 1.0, 1.0));
    for (t, belt) in [(0.0, 77.021249), (0.005, 52.200138), (0.01, 35.377957), (0.02, 16.250059), (0.05, 1.574786), (0.1, 0.032198)] {
        let at = m.at(t);
        near(at.reactions[0].local().0[0], belt, 5e-7 + 1e-6 * 77.0, &format!("P4 belt force at {t} s"));
        near(at.effort[0], belt, 5e-7 + 1e-6 * 77.0, &format!("P4 the slide's effort at {t} s"));
    }
}

/// P10: Modern Robotics' 3-link chain (its first axis skew, turning sqrt 2
/// times its angle): M(theta) theta'' + c + g + J^T F_tip, its tip wrench a
/// load on the third link at the end-effector's frame {4} (the environment
/// on it: -F_tip, given in {4}'s axes), and M, c, g and J^T F_tip alone,
/// against the library's published values.
#[test]
fn p10_modern_robotics_chain() {
    let chain = |g: f64, theta: [f64; 3], rate: [f64; 3], acc: [f64; 3], tip: f64| {
        let arg = |k: usize| format!("(+ {} (* {} (/ t 1s)) (* {} (/ t 1s) (/ t 1s)))", theta[k], rate[k], 0.5 * acc[k]);
        let mut m = model(&format!(
            r#"(import "std") (gravity 0 0 (* -1m/s2 {g}))
            (body l1 (shape (move (sphere 10) 0 0 89.159)) (attr density 1g/cm3))
            (body l2 (shape (move (sphere 10) 280 135.85 89.159)) (attr density 1g/cm3))
            (body l3 (shape (move (sphere 10) 675 16.15 89.159)) (attr density 1g/cm3))
            (assembly a3 (instance l3 l3) (frame c (axes (origin 425 0 89) (x 1 0 0) (z 0 -1 0)))
              (frame ee (axes (origin 817.25 16.15 89.159) (x 0 0 -1) (z 1 0 0))))
            (assembly (a2 (t3 0)) (instance l2 l2) (instance j3 a3 (axes (origin 425 0 89) (x 1 0 0) (z 0 1 0)) c (turn z (rad->deg t3)))
              (frame c (axes (origin 0 0 89) (x 1 0 0) (z 0 -1 0))))
            (assembly (a1 (t2 0) (t3 0)) (instance l1 l1) (instance j2 (a2 t3) (axes (origin 0 0 89) (x 1 0 0) (z 0 1 0)) c (turn z (rad->deg t2)))
              (frame c (axes (origin -500 0 500) (x 1 0 -1) (z -1 0 -1))))
            (assembly (robot (t1 0) (t2 0) (t3 0))
              (instance j1 (a1 t2 t3) (axes (origin -500 0 500) (x 1 0 -1) (z 1 0 1)) c (turn z (* (sqrt 2) (rad->deg t1)))))
            (assembly root (instance arm (robot {} {} {})))
            (def ee (frame-of root/arm/j1/j2/j3 ee))
            (def s (* -1 {tip} (+ (nth (get ee "x") 0) (nth (get ee "y") 0) (nth (get ee "z") 0))))
            (def u (* -1 {tip} (+ (nth (get ee "x") 1) (nth (get ee "y") 1) (nth (get ee "z") 1))))
            (def w (* -1 {tip} (+ (nth (get ee "x") 2) (nth (get ee "y") 2) (nth (get ee "z") 2))))
            (load tip root/arm/j1/j2/j3/l3 (force (* 1N s) (* 1N u) (* 1N w)) (torque (* 1N*mm (* 1000 s)) (* 1N*mm (* 1000 u)) (* 1N*mm (* 1000 w))) (at ee))"#,
            arg(0),
            arg(1),
            arg(2)
        ));
        m.set_mass("root/arm/j1/l1", 3700.0, [0.0, 0.0, 89.159], diag(10267000.0, 10267000.0, 6660000.0));
        m.set_mass("root/arm/j1/j2/l2", 8393.0, [280.0, 135.85, 89.159], diag(15107400.0, 226890000.0, 226890000.0));
        m.set_mass("root/arm/j1/j2/j3/l3", 2275.0, [675.0, 16.15, 89.159], diag(4095000.0, 49443300.0, 49443300.0));
        // N mm per radian: N m.
        let f = m.at(0.0).f.clone();
        assert_eq!(m.coords.iter().map(|c| c.name.as_str()).collect::<Vec<_>>(), ["arm.t1", "arm.t2", "arm.t3"]);
        f.iter().map(|x| x / 1000.0).collect::<Vec<f64>>()
    };
    let (theta, rate) = ([0.1; 3], [0.1, 0.2, 0.3]);
    // To the published digits (7 or 8 decimals).
    let close = |got: Vec<f64>, want: [f64; 3], what: &str| (0..3).for_each(|k| near(got[k], want[k], 5e-8, &format!("P10 {what} [{k}]")));
    // tau = M theta'' + c + g + J^T F_tip.
    close(chain(9.8, theta, rate, [2.0, 1.5, 1.0], 1.0), [74.69616155, -33.06766016, -3.23057314], "tau with the tip wrench");
    close(chain(9.8, theta, rate, [2.0, 1.5, 1.0], 0.0), [74.69616155 - 1.40954608, -33.06766016 - 1.85771497, -3.23057314 - 1.392409], "M theta'' + c + g");
    let mass = [[22.5433380, -0.307146754, -0.00718426391], [-0.307146754, 1.96850717, 0.432157368], [-0.00718426391, 0.432157368, 0.191630858]];
    for (k, column) in mass.iter().enumerate() {
        let mut acc = [0.0; 3];
        acc[k] = 1.0;
        close(chain(0.0, theta, [0.0; 3], acc, 0.0), *column, &format!("M column {k}"));
    }
    close(chain(0.0, theta, rate, [0.0; 3], 0.0), [0.26453118, -0.05505157, -0.00689132], "c");
    close(chain(9.8, theta, [0.0; 3], [0.0; 3], 0.0), [28.40331262, -37.64094817, -5.4415892], "g");
    close(chain(0.0, theta, [0.0; 3], [0.0; 3], 1.0), [1.40954608, 1.85771497, 1.392409], "J^T F_tip");
}

/// Tipping by hand (the robot-arm and phone-stand questions): an arm on a
/// 2 kg base swings a 500 g block, which rides its tip's frame, 60 degrees
/// either way, a 100 g payload's weight a load at the tip. The ground's
/// reaction (root's) puts its zero-moment point at x = -M_y / F_z on the
/// floor, M and F summed by hand over the point masses (m (x (a_z + g) - z
/// a_x)) and the load (its weight at its x); the shoulder's reaction is the
/// block's m (a - g) less the payload's weight.
#[test]
fn tipping_by_hand() {
    let mut m = reacting(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body base (shape (move (box 120 80 20) 0 0 10)) (attr density 1g/cm3))
        (body link (shape (move (box 150 10 10) 75 0 0)) (attr density 1g/cm3) (frame tip (at (move 150 0 0))))
        (body block (shape (box 10 10 10)) (attr density 1g/cm3))
        (assembly (arm (theta 0))
          (instance base base)
          (instance link link (at (move 0 0 100)) (turn y theta) (effort-limits -1e9N*mm 1e9N*mm))
          (instance block block (frame-of link tip)))
        (assembly root (instance a (arm (* 60 (sin (* 180deg/s t))))))
        (load payload root/a/link (force 0 0 -0.981N) (at (frame-of root/a/link tip)))"#,
        &["root", "root/a/link"],
    );
    m.set_mass("root/a/base", 2000.0, [0.0, 0.0, 10.0], diag(0.0, 0.0, 0.0));
    m.set_mass("root/a/link", 0.0, O, diag(0.0, 0.0, 0.0));
    m.set_mass("root/a/block", 500.0, O, diag(0.0, 0.0, 0.0));
    let world = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    let (pi, g, l, h) = (std::f64::consts::PI, 9810.0, 150.0, 100.0);
    for t in [0.0, 0.2, 0.5, 0.8, 1.3] {
        let a = pi / 3.0;
        let (th, w, al) = (a * (pi * t).sin(), a * pi * (pi * t).cos(), -a * pi * pi * (pi * t).sin());
        let (s, c) = th.sin_cos();
        // The block at (L c, 0, H - L s); its acceleration (mm/s^2).
        let (x, z) = (l * c, h - l * s);
        let (ax, az) = (-l * (c * w * w + s * al), l * (s * w * w - c * al));
        let fz = 1e-6 * (2000.0 * g + 500.0 * (az + g)) + 0.981;
        let my = 1e-6 * 500.0 * (z * ax - x * (az + g)) - 0.981 * x;
        let at = m.at(t);
        let (f, mo) = at.reactions[0].in_frame(&world);
        near(f[2], fz, 1e-9 * fz, &format!("the ground's force at {t} s"));
        near(-mo[1] / f[2], -my / fz, 1e-9 * l, &format!("the zero-moment point at {t} s"));
        let (f, _) = at.reactions[1].in_frame(&world);
        near(f[0], 1e-6 * 500.0 * ax, 1e-9 * fz, &format!("the shoulder's force x at {t} s"));
        near(f[2], 1e-6 * 500.0 * (az + g) + 0.981, 1e-9 * fz, &format!("the shoulder's force z at {t} s"));
    }
}

/// P16: a CoreXY (tests/fixtures/dynamics/corexy.ink) on its demo's
/// Lissajous, the head position X, Y the coordinates its root passes and
/// each motor's drive turning (X +- Y) / R, R = 20 / pi mm, its rigid
/// bodies given exact mass properties: the head m_h, the gantry m_g, each
/// idler 8 g and 150 g mm2, and on each motor J about its axis (rotor 5300,
/// pulley 180, its two still idlers, and a 0.5 g grub screw 5 mm off the
/// axis, 12.5). Each belt is open, clamped to the head at both ends, and
/// its route moves: its gantry spans and the two idlers on the gantry it
/// wraps (each turning X / R) ride Y, its clamps X, and its pieces'
/// stations change with both. Its material, by hand (belt A; B is A
/// mirrored, X -> -X), moves along the gantry at (X', Y'); round a gantry
/// idler at (0, Y') + X' u, u the path's direction, so |v|2 = X'2 + Y'2 +
/// 2 X' Y' u.y, u.y integrating to R over the wrap; and elsewhere along the
/// path at X' + Y'. The gantry spans total 2W - 2h - 2R and the side spans
/// H + F - 2R at every X, Y, so its kinetic energy's coefficients are
/// constant: its mass on each axis, and lam (L - G + 2R) coupling them, L
/// its length, G its gantry part (spans and wraps, 2W - 2h - 2R + pi R),
/// B's negative. So M_XX = m_h + m_A + m_B + 2 J / R2 + 4 x 150 / R2, M_YY
/// = m_h + m_g + 4 x 8 + m_A + m_B + 2 J / R2, M_XY = A's coupling less
/// B's (m_A - m_B were the belts' mass per mm equal), and f = M q''
/// (nothing turning weighs on a motor: their axes are along gravity). Each
/// belt's mass per mm is ink's measure (the reference is at it), its
/// section's, 1.4 x 6 mm at 1.2 g/cm3, to 1% (the section measured at
/// ink's resolution). Its effective mass along x where y'' = 0 (t = 3 s)
/// and along y where x'' = 0 (t = 2 s), and the motors' torques, (R /
/// 2)(f_x +- f_y), against M q'', X'' = -110 (pi/2)2 sin(pi t/2), Y'' =
/// -110 (pi/3)2 sin(pi t/3) mm/s2. Dropping the route's motion from
/// `strip` (the pieces' own motion from the acceleration or the coordinate
/// partials, or the pieces' starts held still) fails it. Over its 12 s
/// every step is decided: its pulleys and idlers are round, and the grub
/// screws' departure from round, their axes along gravity on a still
/// frame, moves no effort, so they need no steps for their spin. Each
/// belt's material lies on its spans and wraps, so each part carries its
/// share (a wrap's material on its wheel, a span's on its two ends, by
/// where it is along it): each motor's shaft passes its effort, its rotor,
/// pulley and grub's inertia and the belt's pull on its pulley, the belt's
/// own inertia among that pull (without it the shaft passed less: 13.373
/// N mm of motor B's 13.557 at the catalog's masses, t = 0.2 s).
#[test]
fn p16_corexy_lissajous() {
    let mut m = reacting(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (import "../tests/fixtures/dynamics/corexy.ink" as c)
        (assembly root (instance xy (c.corexy (* 110 (sin (/ (* pi (/ t 1s)) 2))) (* 110 (sin (/ (* pi (/ t 1s)) 3))))))"#,
        &["root/xy/motorA/rotor", "root/xy/motorB/rotor"],
    );
    let (mh, mg, mw, jw) = (300.0, 700.0, 8.0, 150.0);
    let j = 5300.0 + 180.0 + 2.0 * jw + 0.5 * 25.0;
    m.set_mass("root/xy/head", mh, O, diag(0.0, 0.0, 0.0));
    m.set_mass("root/xy/gantry", mg, O, diag(0.0, 0.0, 0.0));
    for k in ["A", "B"] {
        m.set_mass(&format!("root/xy/motor{k}/rotor"), 30.0, O, diag(0.0, 0.0, 5300.0));
        m.set_mass(&format!("root/xy/pulley{k}"), 8.0, O, diag(0.0, 0.0, 180.0));
        m.set_mass(&format!("root/xy/grub{k}"), 0.5, [5.0, 0.0, 0.0], diag(0.0, 0.0, 0.0));
        for w in ["corner", "foot", "out", "back"] {
            m.set_mass(&format!("root/xy/{w}{k}"), mw, O, diag(0.0, 0.0, jw));
        }
    }
    // Its belts' spans are couplings at 0.95: the reference is the ideal
    // mechanism's inertia. The tensions leave nothing out.
    m.without_losses();
    assert!(m.unseparated.is_empty(), "what the tensions leave out: {:?}", m.unseparated);
    let belts: Vec<(f64, f64)> = ["A", "B"].iter().map(|k| m.strip(&format!("root/xy/belt{k}"))).map(|(lam, _, len)| (lam, len)).collect();
    let pi = std::f64::consts::PI;
    let (r, w, foot, h) = (20.0 / pi, 350.0, 350.0, 20.0);
    // Each route's length: 4 W + 2 H + 2 F - 2 h - 4 R + 3 pi R.
    for ((lam, l), top, k) in [(belts[0], 400.0, "A"), (belts[1], 350.0, "B")] {
        let want = 1.2e-3 * 1.4 * 6.0 * (4.0 * w + 2.0 * top + 2.0 * foot - 2.0 * h - 4.0 * r + 3.0 * pi * r);
        near(lam * l, want, 1e-2 * want, &format!("P16 belt {k}'s mass"));
    }
    // Its coupling: its route less its gantry part (spans and wraps, 2 W -
    // 2 h - 2 R + pi R), plus 2 R.
    let couples = |(lam, l): (f64, f64)| lam * (l - (2.0 * w - 2.0 * h - 4.0 * r + pi * r));
    let (ma, mb) = (belts[0].0 * belts[0].1, belts[1].0 * belts[1].1);
    let spin = 2.0 * j / (r * r);
    let mxx = mh + ma + mb + spin + 4.0 * jw / (r * r);
    let myy = mh + mg + 4.0 * mw + ma + mb + spin;
    let mxy = couples(belts[0]) - couples(belts[1]);
    let mass = |t: f64, k: usize| {
        let at = m.at(t);
        1e6 * at.f[k] / at.q[k][2]
    };
    near(mass(3.0, 0), mxx, 1e-6 * mxx, "P16 effective mass along x");
    near(mass(2.0, 1), myy, 1e-6 * myy, "P16 effective mass along y");
    for t in [0.5, 1.0, 1.5, 2.0, 3.0, 4.5] {
        let (xa, ya) = (-110.0 * (pi / 2.0).powi(2) * (pi * t / 2.0).sin(), -110.0 * (pi / 3.0).powi(2) * (pi * t / 3.0).sin());
        let (fx, fy) = (1e-6 * (mxx * xa + mxy * ya), 1e-6 * (mxy * xa + myy * ya));
        let f = &m.at(t).f;
        near(r / 2.0 * (f[0] + f[1]), r / 2.0 * (fx + fy), 1e-9 + 1e-6 * 1.5, &format!("P16 motor A at {t} s"));
        near(r / 2.0 * (f[0] - f[1]), r / 2.0 * (fx - fy), 1e-9 + 1e-6 * 1.5, &format!("P16 motor B at {t} s"));
        let at = m.at(t);
        for (k, motor) in ["root/xy/motorA/rotor", "root/xy/motorB/rotor"].into_iter().enumerate() {
            let a = m.actuators.iter().position(|x| x.path == motor).unwrap();
            let noise = at.reactions[k].rounding().1 + 1e-9 * (1.0 + at.effort[a].abs());
            near(at.reactions[k].local().1[2], at.effort[a], noise, &format!("P16 {motor}'s shaft passes its effort at {t} s (to rounding)"));
        }
    }
    let s = m.steps((0.0, 12.0));
    assert!(s.loose.is_empty() && s.unchecked.is_empty(), "{} loose, {} unchecked of {}", s.loose.len(), s.unchecked.len(), s.ends.len() - 1);
    assert!(s.ends.len() < 400 && s.basis.contains("root/xy/motorA/rotor with root/xy/pulleyA, root/xy/grubA"), "{} steps: {}", s.ends.len() - 1, s.basis);
}

/// P16's CoreXY moving along y alone, its head still on the gantry (y = 60
/// sin pi t), braking at t = 0.25 s (y'' = -60 pi2 sin 45 deg, y' = 60 pi
/// cos 45 deg): what the gantry carries, its 700 g, the head's 300 g and
/// its four 8 g idlers (none turns: they turn with x), gives the belts its
/// power, M y'' y', and their spans' 0.95 loses 1 / 0.95 - 1 of it there,
/// at the gantry's slide, so the motors take back that much less than the
/// ideal mechanism gives them: their power, sum of effort times rate, is
/// the ideal one less (1 / 0.95 - 1) M y'' y' (its belts here massless and
/// its still idlers' spin free, so that every loss is the gantry's). The
/// loss was charged only along each part's own joint, so nothing was
/// charged where the gantry pulls the belts through its idlers' axles.
#[test]
fn p16_corexy_braking_loses_at_its_gantry() {
    let build = || {
        let mut m = model(
            r#"(import "std") (gravity 0 0 -9.81m/s2)
            (import "../tests/fixtures/dynamics/corexy.ink" as c)
            (assembly root (instance xy (c.corexy 0 (* 60 (sin (* 180deg/s t))))))"#,
        );
        let (mh, mg, mw) = (300.0, 700.0, 8.0);
        m.set_mass("root/xy/head", mh, O, diag(0.0, 0.0, 0.0));
        m.set_mass("root/xy/gantry", mg, O, diag(0.0, 0.0, 0.0));
        for k in ["A", "B"] {
            m.set_mass(&format!("root/xy/motor{k}/rotor"), 30.0, O, diag(0.0, 0.0, 5300.0));
            m.set_mass(&format!("root/xy/pulley{k}"), 8.0, O, diag(0.0, 0.0, 180.0));
            m.set_mass(&format!("root/xy/grub{k}"), 0.5, [5.0, 0.0, 0.0], diag(0.0, 0.0, 0.0));
            for w in ["corner", "foot", "out", "back"] {
                m.set_mass(&format!("root/xy/{w}{k}"), mw, O, diag(0.0, 0.0, 0.0));
            }
            m.set_strip_mass(&format!("root/xy/belt{k}"), 0.0);
        }
        (m, mh + mg + 4.0 * mw)
    };
    let ((lossy, carried), (mut ideal, _)) = (build(), build());
    ideal.without_losses();
    let (pi, t) = (std::f64::consts::PI, 0.25);
    let (acc, speed) = (-60.0 * pi * pi * (pi * t).sin(), 60.0 * pi * (pi * t).cos());
    let power = |m: &Model| {
        let at = m.at(t);
        (0..m.actuators.len()).map(|a| at.effort[a] * at.rate[a] * pi / 180.0).sum::<f64>()
    };
    let given = 1e-6 * carried * acc * speed;
    assert!(given < 0.0, "braking: the gantry gives the belts {given} N mm/s");
    near(power(&lossy) - power(&ideal), -(1.0 / 0.95 - 1.0) * given, 1e-9 * power(&ideal).abs(), "P16's motors' power, less the gantry's loss");
}

/// The deploy pulse of Limestone's demo.ink: a logistic speed pulse.
const PULSE: &str = r#"(fn lncosh (z) (ln (* 0.5 (+ (exp z) (exp (neg z))))))
    (fn pulse (t0 t1 w a)
      (let ((tau (/ w (* 4 a))) (h (/ (- t1 t0) (* 2 tau))) (s (/ (- (/ t 1s) (* 0.5 (+ t0 t1))) tau)))
        (* w tau (+ h (lncosh (* 0.5 (+ s h))) (neg (lncosh (* 0.5 (- s h))))))))"#;

/// P17: Limestone's intake arm deploying, its hopper sliding out with it
/// by the arm's angle (the pin in the walls' slot), with the masses the
/// reference states (Limestone's own carry no material, so the catalog
/// assembly refuses).
#[test]
fn p17_intake_deploy() {
    let mut m = model(&format!(
        r#"(import "std") (gravity 0 0 -9.81m/s2) {PULSE}
        (body block (shape (box 10 10 10)) (attr density 1g/cm3))
        (fn hopper (p) (let ((a (* 0.017453 p))) (+ 124.838713 (* 146.64208 (- (* 0.240079 (sin a)) (cos a))) (* 90.817561 (+ (* 0.240079 (cos a)) (sin a))))))
        (assembly (intake (deploy 130))
          (instance arm block (at (rot-y 90) (move 0 -304.643 168.888)) (turn z (- 130 deploy)))
          (instance hopper block (at (move 0 (neg (hopper (- 130 deploy))) 0)))
          (instance motor block (at (move 0 0 300)) (turn z (* 48.796296 (- 130 deploy))) (effort-limits -1e6N*mm 1e6N*mm)))
        (assembly root (instance robot (intake (- 130 (pulse 0.05 (+ 0.05 (/ 130 360)) 360 5400)))))"#
    ));
    let (mass, dy, dz, io) = (5449.4, 156.50, 156.27, 0.33518e9);
    // The arm's frame turns the world's x to its z: its centre (dy, dz) in
    // the world's y and z is (-dz, dy, 0) in its own.
    m.set_mass("root/robot/arm", mass, [-dz, dy, 0.0], diag(0.0, 0.0, io - mass * (dy * dy + dz * dz)));
    m.set_mass("root/robot/hopper", 2957.5, O, diag(0.0, 0.0, 0.0));
    m.set_mass("root/robot/motor", 0.0, O, diag(0.0, 0.0, 0.0));
    let table = [(0.0, 304.2098), (0.05, 918.9786), (0.1, 279.3217), (0.2, -49.3014), (0.3, -212.4634), (0.411, -893.7510), (0.6, -241.3865)];
    for (t, motor) in table {
        near(m.at(t).effort[0], motor, 5e-5 + 1e-6 * 919.08, &format!("P17 motor at {t} s"));
    }
    // Holding, stowed: m g dy; the coordinate's effort per degree of deploy.
    let hold = model(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body block (shape (box 10 10 10)) (attr density 1g/cm3))
        (assembly (intake (deploy 130)) (instance arm block (at (rot-y 90) (move 0 -304.643 168.888)) (turn z (- 130 deploy)) (effort-limits -1e6N*mm 1e6N*mm)))
        (assembly root (instance robot (intake 130)))"#,
    );
    let mut hold = hold;
    hold.set_mass("root/robot/arm", mass, [-dz, dy, 0.0], diag(0.0, 0.0, io - mass * (dy * dy + dz * dz)));
    near(hold.at(0.0).effort[0] / 1000.0, 8.36627, 5e-6, "P17 holding stowed");
}

/// P18, inverse: Limestone's climber pulling the robot (68.779 kg, the
/// reference's mass) up its bar, the robot placed by the climb; the bar's
/// load on it (its reaction), m (g - c''), which the motor's torque lifts
/// through the train: tau k (k = 1551.6814 rad per m of climb).
#[test]
fn p18_climb() {
    let mut m = reacting(
        &format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) {PULSE}
        (body block (shape (box 10 10 10)) (attr density 1g/cm3))
        (assembly (climber (climb 209.55))
          (instance robot block (at (move 0 0 (neg climb))))
          (instance motor block (at) (turn z (* 88.904796 climb)) (effort-limits -1e6N*mm 1e6N*mm)))
        (assembly root (instance c (climber (- 209.55 (pulse 9.2 (+ 9.2 (/ 158.75 300)) 300 1500)))))"#
        ),
        &["root/c/robot"],
    );
    m.set_mass("root/c/robot", 68779.0, O, diag(0.0, 0.0, 0.0));
    m.set_mass("root/c/motor", 0.0, O, diag(0.0, 0.0, 0.0));
    for (t, tau) in [(9.2, -501.3143), (9.3, -462.7064), (9.464, -434.8634), (9.6, -417.5584), (9.729, -368.3516)] {
        let at = m.at(t);
        near(at.effort[0], tau, 5e-5 + 1e-6 * 501.3, &format!("P18 at {t} s"));
        near(at.reactions[0].local().0[2], -tau * 1551.6814 / 1000.0, 5e-5 * 1.552 + 1e-6 * 777.9, &format!("P18 the bar's load at {t} s"));
    }
    near(m.at(0.0).effort[0], -434.8328, 5e-5, "P18 holding");
}

/// The model of `src` built as states (build_states), for a simulation.
fn states(src: &str) -> Model {
    let doc = Doc::load_str_in(&REGISTRY, src, "validation.ink", Path::new("."), &[]).unwrap_or_else(|e| panic!("{} at {}", e.message, e.line));
    let g = gravity(&doc).expect("a gravity");
    build_states(&doc.reloader().unwrap(), "root", g).unwrap_or_else(|e| panic!("{} ({:?})", e.message, e.hint))
}

/// The motion `m` (a model of states) simulated over `span` from `start`.
fn simulated(m: Model, span: (f64, f64), start: &[(f64, f64)]) -> simulate::Trajectory {
    let traj = simulate::simulate(m, span, start).unwrap();
    assert!(traj.stopped.is_none(), "{:?}", traj.stopped);
    traj
}

/// P4 forward: the belt carriage (1.5 kg) from rest at x = 0, its motor
/// (the Pololu 4692's line at 24 V) commanded past its limits, so at full
/// voltage: v = v_f (1 - e^(-t / tau)), x(0.1 s) = 57.520173 mm, 300 mm at
/// 0.467399050 s (to 1e-6). The motor's effort along
/// the simulated motion is its line at that speed.
#[test]
fn p4_simulated() {
    let mut m = states(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.5s)
        (def r (/ 60 pi))
        (body carriage (shape (box 40 40 20)) (attr density 1g/cm3))
        (body rotor (shape (cylinder-z 5 10)) (attr density 1g/cm3))
        (assembly (axis (x 0))
          (instance c carriage (at) (slide x x))
          (instance motor rotor (at (move 0 100 0)) (turn z (rad->deg (/ x r)))
            (effort-limits (dc-torque-min 1470.997N*mm 330rpm 1 1 (+ 0deg/s (diff (rad->deg (/ x r)) t)))
                           (dc-torque-max 1470.997N*mm 330rpm 1 1 (+ 0deg/s (diff (rad->deg (/ x r)) t))))))
        (assembly root (instance a (axis (simulate 0))))
        (command root/a/motor 1e9)"#,
    );
    m.set_mass("root/a/c", 1500.0, O, diag(1.0, 1.0, 1.0));
    m.set_mass("root/a/motor", 0.0, O, diag(0.0, 0.0, 0.0));
    let traj = simulated(m, (0.0, 0.5), &[(0.0, 0.0)]);
    let (r, ts, wf) = (60.0 / std::f64::consts::PI, 1470.997, 330.0 * 2.0 * std::f64::consts::PI / 60.0);
    let (vf, tau) = (wf * r, 1.5 * r * r * 1e-3 * wf / ts);
    for (t, x) in [(0.005, 0.56612), (0.01, 2.01327), (0.02, 6.50646), (0.05, 24.69008), (0.1, 57.520173)] {
        let got = traj.at(0, t).unwrap();
        let exact = vf * (t - tau * (1.0 - (-t / tau).exp()));
        near(got[0], exact, 1e-6 * exact.max(1.0), &format!("P4 x({t}) against the closed form"));
        near(got[0], x, 1e-5 * x, &format!("P4 x({t}) against its table"));
        near(got[1], vf * (1.0 - (-t / tau).exp()), 1e-6 * vf, &format!("P4 v({t})"));
    }
    near(traj.at(0, 0.467399050).unwrap()[0], 300.0, 1e-6 * 300.0, "P4 300 mm at 0.467399050 s");
}

/// P7 forward: the hoist (5 kg on a 25 mm drum, the Pololu 4745's line at
/// 12 V) from rest at full voltage: x(1 s) = 163.892229 mm, 200 mm at
/// 1.218151787 s.
#[test]
fn p7_simulated() {
    let mut m = states(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1.3s)
        (body load (shape (box 40 40 20)) (attr density 1g/cm3))
        (body drum (shape (cylinder-z 25 10)) (attr density 1g/cm3))
        (assembly (hoist (z 0))
          (instance load load (at) (slide z z))
          (instance drum drum (at (rot-x 90) (move 0 0 300)) (turn z (rad->deg (/ z 25)))
            (effort-limits (dc-torque-min 3334.261N*mm 100rpm 1 1 (+ 0deg/s (diff (rad->deg (/ z 25)) t)))
                           (dc-torque-max 3334.261N*mm 100rpm 1 1 (+ 0deg/s (diff (rad->deg (/ z 25)) t))))))
        (assembly root (instance h (hoist (simulate 0))))
        (command root/h/drum 1e9)"#,
    );
    m.set_mass("root/h/load", 5000.0, O, diag(1.0, 1.0, 1.0));
    m.set_mass("root/h/drum", 0.0, O, diag(0.0, 0.0, 0.0));
    let traj = simulated(m, (0.0, 1.3), &[(0.0, 0.0)]);
    for (t, x) in [(0.005, 0.17914), (0.01, 0.61711), (0.02, 1.89754), (0.05, 6.66129), (1.0, 163.892229)] {
        near(traj.at(0, t).unwrap()[0], x, 1e-5 * x.max(1.0), &format!("P7 x({t})"));
    }
    near(traj.at(0, 1.0).unwrap()[1], 165.51673, 1e-6 * 165.5, "P7 its steady speed");
    near(traj.at(0, 1.218151787).unwrap()[0], 200.0, 1e-6 * 200.0, "P7 200 mm at 1.218151787 s");
}

/// P1 released: the bar (162 g, 300 mm, pinned at an end, no actuator)
/// let go at 60 degrees swings to -60 in half its period, 2 sqrt(I_O / m g
/// d) K(sin 30), K by the arithmetic-geometric mean; its energy holds to
/// the integrator's tolerance.
#[test]
fn p1_released() {
    let mut m = states(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle)))
        (assembly root (instance p (pendulum (simulate 60))))"#,
    );
    let (mass, l, w, h) = (162.0, 300.0, 20.0, 10.0);
    m.set_mass("root/p/bar", mass, [0.0, l / 2.0, 0.0], diag(mass * (l * l + h * h) / 12.0, mass * (w * w + h * h) / 12.0, mass * (l * l + w * w) / 12.0));
    let traj = simulated(m, (0.0, 2.0), &[(60.0, 0.0)]);
    let (io, mgd): (f64, f64) = (4.865400e-3, 0.162 * 9.81 * 0.15);
    let k = {
        let (mut a, mut g) = (1.0f64, (1.0f64 - 0.25).sqrt());
        for _ in 0..10 {
            (a, g) = (0.5 * (a + g), (a * g).sqrt());
        }
        std::f64::consts::FRAC_PI_2 / a
    };
    let half = 2.0 * (io / mgd).sqrt() * k;
    let energy = |t: f64| {
        let [th, w, ..] = traj.at(0, t).unwrap();
        0.5 * io * (w.to_radians()).powi(2) + mgd * (1.0 - th.to_radians().cos())
    };
    let e0 = energy(0.0);
    for t in [0.1, 0.37, half, 1.0, 1.5, 2.0] {
        near(energy(t), e0, 1e-7 * e0, &format!("P1 energy at {t} s"));
    }
    let [th, w, ..] = traj.at(0, half).unwrap();
    near(th, -60.0, 1e-5, "P1 at the far end of its swing");
    near(w, 0.0, 1e-4, "P1 at rest there");
}

/// A wheel (published spin inertia 1e5 g mm^2, massless motor) on an equal
/// belt of efficiency E, 100 deg/s at the start and its motor commanded 10
/// N mm: it gets E x 10 N mm, alpha = E 10 / I, 2864.788976 deg/s^2 at 0.5,
/// 2291.831181 at 0.4 (before the fix: a loss turning with the power's
/// direction made the iteration stop short: 0, and at 0.4 backward).
#[test]
fn a_lossy_belt_drives_its_wheel_by_its_efficiency() {
    for e in [0.5, 0.4] {
        let traj = simulated(
            states(&format!(
                r#"(import "std") (gravity 0 0 0) (t-span 0 0.01s)
                (body rotor (shape (cylinder-z 4 10)) (attr mass 0) (frame belt (at)))
                (body wheel (shape (cylinder-z 20 20)) (attr density 7.85g/cm3) (attr spin-inertia 100000g*mm2) (frame belt (at)))
                (assembly (drive (a 0))
                  (instance motor rotor (at) (turn z a) (effort-limits -1000N*mm 1000N*mm))
                  (instance output wheel (at (move 100 0 0)) (turn z a))
                  (relate motor output (belt-drive 10 10 0 {e})))
                (assembly root (instance d (drive (simulate 0 100))))
                (command root/d/motor 10)"#
            )),
            (0.0, 0.01),
            &[(0.0, 100.0)],
        );
        let alpha = (e * 10.0 / 1e-1f64).to_degrees();
        let [x, v, a, ..] = traj.at(0, 0.005).unwrap();
        near(a, alpha, 1e-9 * alpha, &format!("its acceleration at E = {e}"));
        near(v, 100.0 + alpha * 0.005, 1e-9 * v, &format!("its speed at E = {e}"));
        near(x, 0.5 + 0.5 * alpha * 0.005 * 0.005, 1e-9, &format!("its angle at E = {e}"));
    }
}

/// An idler (3 g, its own diagonal 106 g mm^2) turned alone by a
/// simulated coordinate, its spin-inertia `spin` g mm^2, posed `at`, its
/// motor commanded `push` N mm: the model, and the source.
fn idler(at: &str, spin: &str, push: &str) -> String {
    format!(
        r#"(import "std") (gravity 0 0 0) (t-span 0 0.01s)
        (body idler (shape (cylinder-z 6 10)) (attr density 2.7g/cm3) (attr spin-inertia (* 1g*mm2 {spin})))
        (assembly (drive (a 0)) (instance idler idler {at} (turn z a) (effort-limits -1000N*mm 1000N*mm)))
        (assembly root (instance d (drive (simulate 0 100))))
        (command root/d/idler {push})"#
    )
}

/// An idler given no spin-inertia, posed where its centre measures a hair
/// off its axis, sums m r^2 there and m r^2 less in its own inertia to a
/// crumb above 0, below, or 0 exactly (each pose's, measured); a real
/// 1e-12 g mm^2 is as little beside its heft. Each carries none: pushed,
/// each is refused as having no inertia, and unpushed each keeps its speed.
/// Before the fix the crumb above 0 took the push over it, 2e22 deg/s^2,
/// and its steps shrank to nothing, where the one below was refused.
#[test]
fn an_inertia_within_rounding_of_its_heft_is_none_whatever_its_sign() {
    let poses = [
        ("(at (rot-x 37) (move 600 0 0))", "0"),
        ("(at (rot-x 33.3) (rot-y 21.7) (move 123.4 56.7 89.1))", "0"),
        ("(at (move 600 0 0))", "0"),
        ("(at (rot-x 37) (move 600 0 0))", "1e-12"),
    ];
    let gm: Vec<f64> = poses.iter().map(|(at, spin)| states(&idler(at, spin, "1")).state_near(0.0, &[[0.0, 100.0, 0.0]]).gm[0][0]).collect();
    assert!(gm[0] > 0.0 && gm[1] < 0.0 && gm[2] == 0.0 && gm[3] > 0.0, "the poses' crumbs: {gm:?}");
    for (at, spin) in poses {
        let why = simulate::simulate(states(&idler(at, spin, "1")), (0.0, 0.01), &[(0.0, 100.0)]).err().unwrap_or_default();
        assert!(why.contains("nothing with mass moves along d.a"), "{at} {spin}: {why}");
        let [x, v, a, ..] = simulated(states(&idler(at, spin, "0")), (0.0, 0.01), &[(0.0, 100.0)]).at(0, 0.005).unwrap();
        assert!((x, v, a) == (0.5, 100.0, 0.0), "{at} {spin}: {x} {v} {a}");
    }
}

/// The same idler with a small real spin-inertia, 1e-5 g mm^2 (a ten-
/// millionth of its heft), posed as the crumb above 0 was, commanded 1e-8 N
/// mm: alpha = tau / I, 1000 rad/s^2, to 1e-9.
#[test]
fn a_small_real_inertia_is_simulated() {
    let traj = simulated(states(&idler("(at (rot-x 37) (move 600 0 0))", "1e-5", "1e-8")), (0.0, 0.01), &[(0.0, 100.0)]);
    let alpha = 1000f64.to_degrees();
    let [x, v, a, ..] = traj.at(0, 0.005).unwrap();
    near(a, alpha, 1e-9 * alpha, "its acceleration");
    near(v, 100.0 + alpha * 0.005, 1e-9 * v, "its speed");
    near(x, 0.5 + 0.5 * alpha * 0.005 * 0.005, 1e-9, "its angle");
}

/// What inverse dynamics refuses stops a simulation, saying why (before
/// the fix it went on with the ideal coupling): a belt passing power that
/// gives no efficiency, and one whose output turns against an open belt.
#[test]
fn a_simulation_stops_where_inverse_dynamics_refuses() {
    let belt = |turn: &str, kind: &str| {
        let m = states(&format!(
            r#"(import "std") (gravity 0 0 0) (t-span 0 0.01s)
            (body rotor (shape (cylinder-z 4 10)) (attr mass 0) (frame belt (at)))
            (body wheel (shape (cylinder-z 20 20)) (attr density 7.85g/cm3) (attr spin-inertia 100000g*mm2) (frame belt (at)))
            (assembly (drive (a 0))
              (instance motor rotor (at) (turn z a) (effort-limits -1000N*mm 1000N*mm))
              (instance output wheel (at (move 100 0 0)) (turn z {turn}))
              (relate motor output {kind}))
            (assembly root (instance d (drive (simulate 0 100))))
            (command root/d/motor 10)"#
        ));
        simulate::simulate(m, (0.0, 0.01), &[(0.0, 100.0)]).err().unwrap()
    };
    let unknown = belt("a", "(belt-drive 10 10 0)");
    assert!(unknown.contains("gives no efficiency, so its loss is unknown"), "{unknown}");
    let crossed = belt("(neg a)", "(belt-drive 10 10 0 0.9)");
    assert!(crossed.contains("does not match the joints' motion"), "{crossed}");
}

/// A 1 N push from 0.4 to 0.5 s on a 1 kg cart, written as selects on t:
/// 100 mm/s after it, 55 mm at 1 s (before the fix: steps grown past it
/// stepped over it, 0 mm). A step over which the command may jump in t is
/// halved down to the switch.
#[test]
fn a_pulse_is_not_stepped_over() {
    let traj = simulated(
        states(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body block (shape (box 50 50 50)) (attr mass 1000g))
            (assembly (axis (s 0)) (instance c block (at) (slide x s) (effort-limits -100N 100N)))
            (assembly root (instance a (axis (simulate 0))))
            (command root/a/c (select (- (/ t 1s) 0.4) 0 (select (- (/ t 1s) 0.5) 1 0)))"#,
        ),
        (0.0, 1.0),
        &[(0.0, 0.0)],
    );
    let [x, v, ..] = traj.at(0, 1.0).unwrap();
    near(v, 100.0, 1e-6 * 100.0, "its speed after the pulse");
    near(x, 55.0, 1e-6 * 55.0, "where it is at 1 s");
}

/// A 5 kg block raised by a lead screw (8 mm lead, efficiency E, rotor and
/// screw of 0.001 g) under a constant 1500 N mm: alpha = (tau - m g L / 2 pi
/// E) / (m (L / 2 pi)^2 / E), at E = 0.49, 0.45 and 0.35, screws that
/// self-lock when driven back (before the fix: at E <= 0.5 the motion
/// started from the lossless acceleration and Newton did not converge).
/// The same screw at 0.35 holding the block with no command self-locks:
/// the simulation stops, saying friction holds it at rest.
#[test]
fn a_self_locking_screw_lifts_by_its_efficiency_and_holds_by_friction() {
    let screw = |e: f64, command: &str| {
        states(&format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.02s)
            (body rotor (shape (cylinder-z 4 30)) (attr mass 0.001g))
            (body screw (shape (move (cylinder-z 4 200) 0 0 100)) (attr mass 0.001g))
            (body block (shape (box 80 80 80)) (attr mass 5kg))
            (assembly (lift (a 0))
              (instance rotor rotor (at (move 0 0 -40)) (turn z a) (effort-limits -2000N*mm 2000N*mm))
              (instance screw screw (frame-of rotor))
              (instance c block (at (move 0 0 100)) (slide z (* 8 (/ a 360))))
              (relate screw c (lead-screw 8 {e})))
            (assembly root (instance l (lift (simulate 0))))
            {command}"#
        ))
    };
    for e in [0.49, 0.45, 0.35] {
        let traj = simulated(screw(e, "(command root/l/rotor 1500)"), (0.0, 0.02), &[(0.0, 0.0)]);
        let lead = 8.0 / (2.0 * std::f64::consts::PI);
        let alpha = 1e6 * (1500.0 - 5000.0 * 9.81e-3 * lead / e) / (5000.0 * lead * lead / e);
        let want = alpha.to_degrees();
        near(traj.at(0, 0.01).unwrap()[2], want, 1e-6 * want, &format!("its acceleration at E = {e}"));
        near(traj.at(0, 0.01).unwrap()[0], 0.5 * want * 1e-4, 1e-6 * want * 1e-4, &format!("its angle at 0.01 s at E = {e}"));
    }
    let why = simulate::simulate(screw(0.35, ""), (0.0, 0.02), &[(0.0, 0.0)]).err().unwrap();
    assert!(why.contains("at rest at t = 0 is held by friction") && why.contains("self-locks is driven back"), "{why}");
}

/// What stops a simulation says what it is: a command past every limit
/// (infinite), and a mechanism whose steps its stage budget cannot carry
/// (1 kg let go 1 mm out on a 1.6 kHz spring: its 1600 undamped periods a
/// second, each followed to the tolerance, by any accurate method; it took
/// 110022 explicit steps and 4.3 s before the budget), at
/// once; a command with a square root's
/// unbounded slope where it reaches 0 does not stall the steps (it is no
/// switch); an actuator moving a simulated and a written coordinate is
/// refused.
#[test]
fn a_simulation_says_why_it_stops() {
    let pendulum = |command: &str, limits: &str| {
        states(&format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
            (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
            (def swing (simulate 60))
            (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits {limits})))
            (assembly root (instance p (pendulum swing)))
            (command root/p/bar {command})"#
        ))
    };
    let inf = simulate::simulate(pendulum("(exp 1000)", "(neg (exp 1000)) (exp 1000)"), (0.0, 1.0), &[(60.0, 0.0)]).err().unwrap();
    assert!(inf.contains("the command on root/p/bar is inf") && inf.contains("do not bound it"), "{inf}");
    let smooth = simulate::simulate(pendulum("(sqrt (- swing 59))", "-1e4N*mm 1e4N*mm"), (0.0, 1.0), &[(60.0, 0.0)]).unwrap();
    assert!(smooth.stopped.is_none() && smooth.evaluations < 4000, "{:?} {}", smooth.stopped, smooth.evaluations);
    let stiff = states(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (def u (simulate 1))
        (assembly (axis (s 0)) (instance c block (at) (slide x s)))
        (assembly root (instance a (axis u)))
        (load spring root/a/c (force (* -1e5N u) 0 0))"#,
    );
    let t0 = std::time::Instant::now();
    let traj = simulate::simulate(stiff, (0.0, 1.0), &[(1.0, 0.0)]).unwrap();
    let (at, why) = traj.stopped.unwrap();
    assert!(at < 1.0 && why.contains("stages evaluated") && t0.elapsed().as_secs_f64() < 20.0, "{at} {why}");
    let mixed = states(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (axes (x 0) (y 0)) (instance c block (at) (slide x x)) (instance d block (at (move 0 0 50)) (slide y y))
          (instance m block (at (move 0 0 100)) (slide z (+ x y)) (effort-limits -100N 100N)))
        (assembly root (instance a (axes (simulate 0) (* 10 (/ t 1s)))))
        (command root/a/m 1)"#,
    );
    let why = simulate::simulate(mixed, (0.0, 1.0), &[(0.0, 0.0)]).err().unwrap();
    assert!(why.contains("root/a/m moves the simulated coordinates and a.y, whose motion is written"), "{why}");
}

/// The call's time limit stops the stepping where it got to (before, it
/// stepped on to the stage budget or the span's end, and ink mcp killed the
/// call before it answered): a late moment asked past it says the limit
/// stopped it, its stages a step past the limit's; asked again with time
/// left, it steps on from there to the motion the unlimited run steps.
#[test]
fn the_time_limit_stops_a_simulation_where_it_got_to() {
    // Its own span: the test's limit is the simulation's over it.
    const SPAN: (f64, f64) = (0.0, 0.913);
    let pendulum = || {
        states(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.913s)
            (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
            (def swing (simulate 60))
            (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle)))
            (assembly root (instance p (pendulum swing)))"#,
        )
    };
    let whole = simulate::simulate(pendulum(), SPAN, &[(60.0, 0.0)]).unwrap();
    assert!(whole.stopped.is_none() && whole.evaluations > 400, "{:?} {}", whole.stopped, whole.evaluations);
    let sim = simulate::Simulation::start(pendulum(), SPAN, &[(60.0, 0.0)]).unwrap();
    let begun = sim.reach(0.0, |traj| traj.evaluations);
    let limit = begun + 100;
    let limited = simulate::LimitAfter::new(SPAN, limit);
    let late = sim.at(0, 0.9);
    let r = sim.stepped();
    // A moment before where it stopped is read as it is, still stopped.
    assert!(sim.at(0, 0.01).is_ok() && sim.stepped().stopped == r.stopped);
    drop(limited);
    let (at, why) = r.stopped.clone().expect("stopped by the limit");
    assert!(at < 0.9 && r.end == at && why.contains("time limit"), "{at} {why}");
    // Stopped at the first step past the limit: a step's stages over it.
    assert!(r.evaluations >= limit && r.evaluations < limit + 20, "{} stages, the limit at {limit}", r.evaluations);
    let late = late.expect_err("no moment past where it stopped");
    assert!(late.contains("time limit"), "{late}");
    let x = sim.at(0, 0.9).unwrap();
    let r = sim.stepped();
    assert!(r.stopped.is_none() && r.end == SPAN.1, "{:?}", r.stopped);
    assert!((x[0] - whole.at(0, 0.9).unwrap()[0]).abs() < 1e-3, "{} against {}", x[0], whole.at(0, 0.9).unwrap()[0]);
}

/// Sticking: a 1 kg carriage on a 1 mm lead screw of
/// efficiency 0.9 (it does not self-lock), raised 360 degrees by a PD
/// command with its weight fed forward, comes to rest where the command
/// lies within what the screw's friction holds still: the simulation stops
/// there, saying friction holds it (it said "stiff"); the same lift at
/// efficiency 1 settles. A 100 Hz PWM command over 1 s steps each of its
/// 200 switches across at once (it ran out of stages at 0.79 s); an E =
/// 0.35 screw driven back at speed by a command names its negative
/// effective inertia (it said static friction held it).
#[test]
fn a_simulation_names_friction_switches_and_negative_inertia() {
    let lift = |e: f64| {
        states(&format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
            (body rotor (shape (cylinder-z 10 30)) (attr density 7.85g/cm3))
            (body block (shape (box 50 50 50)) (attr mass 1kg))
            (def a (simulate 0))
            (assembly (lift (a 0))
              (instance rotor rotor (at (move 0 0 -40)) (turn z a) (effort-limits -100N*mm 100N*mm))
              (instance c block (at (move 0 0 100)) (slide z (/ a 360)))
              (relate rotor c (lead-screw 1 {e})))
            (assembly root (instance l (lift a)))
            (command root/l/rotor (+ 1.5613 (* 0.05 (- 360 a)) (* -0.0025 (diff a t))))"#
        ))
    };
    let held = simulate::simulate(lift(0.9), (0.0, 1.0), &[(0.0, 0.0)]).unwrap();
    let (at, why) = held.stopped.unwrap();
    assert!(at < 1.0 && why.contains("comes to rest") && why.contains("held by friction"), "{at} {why}");
    let free = simulate::simulate(lift(1.0), (0.0, 1.0), &[(0.0, 0.0)]).unwrap();
    assert!(free.stopped.is_none() && (free.at(0, 1.0).unwrap()[0] - 360.0).abs() < 0.01, "{:?}", free.stopped);
    let pwm = simulate::simulate(
        states(
            r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
            (body block (shape (box 50 50 50)) (attr mass 1000g))
            (assembly (axis (s 0)) (instance c block (at) (slide x s) (effort-limits -100N 100N)))
            (assembly root (instance a (axis (simulate 0))))
            (command root/a/c (select (- (- (* (/ t 1s) 100) (floor (* (/ t 1s) 100))) 0.5) 1 0))"#,
        ),
        (0.0, 1.0),
        &[(0.0, 0.0)],
    )
    .unwrap();
    assert!(pwm.stopped.is_none() && pwm.evaluations < 10000, "{:?} {}", pwm.stopped, pwm.evaluations);
    // Half of each 10 ms period at 1 m/s^2: 0.5 m/s at 1 s, and 250 + 1.25
    // mm (each period's first half adds a quarter of its square's area).
    near(pwm.at(0, 1.0).unwrap()[1], 500.0, 1e-6 * 500.0, "the PWM's speed at 1 s");
    near(pwm.at(0, 1.0).unwrap()[0], 251.25, 1e-6 * 251.25, "where the PWM has moved it at 1 s");
    let lower = states(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.5s)
        (body rotor (shape (cylinder-z 4 30)) (attr density 7.85g/cm3))
        (body screw (shape (move (cylinder-z 4 200) 0 0 100)) (attr density 7.85g/cm3))
        (body block (shape (box 80 80 80)) (attr mass 5kg))
        (assembly (lift (a 0))
          (instance rotor rotor (at (move 0 0 -40)) (turn z a) (effort-limits -2000N*mm 2000N*mm))
          (instance screw screw (frame-of rotor))
          (instance c block (at (move 0 0 100)) (slide z (* 8 (/ a 360))))
          (relate screw c (lead-screw 8 0.35)))
        (assembly root (instance l (lift (simulate 0 -720))))
        (command root/l/rotor -53.53)"#,
    );
    let why = simulate::simulate(lower, (0.0, 0.5), &[(0.0, -720.0)]).err().unwrap();
    assert!(why.contains("driven back at speed") && why.contains("negative") && why.contains("write the motion"), "{why}");
}

/// A simulated def passed to several instances: itself
/// once it is a coordinate; any other argument reading it (twice it, its
/// negative, a sum with another def, or itself again) is passed as it is,
/// so its bodies move along the def's coordinate. 1 kg blocks sliding
/// along gravity (+x, 9.81 m/s^2), by Lagrange: th and 2 th, T = 5/2 m
/// th'^2 and the weights' work 3 m g th, so th'' = 3/5 g = 5886 mm/s^2;
/// th and -th, none; th, ph and th + ph, each 2/3 g = 6540; th twice, g.
/// (The second argument took the def's param as its own coordinate:
/// 6540, a false friction stop at t = 0, 4905 and 9810.)
#[test]
fn a_simulated_def_read_by_several_arguments_moves_them_all() {
    let at = |root: &str, defs: &str| {
        let m = states(&format!(
            r#"(import "std") (gravity 9.81m/s2 0 0) (t-span 0 0.1s) {defs}
            (body block (shape (box 50 50 50)) (attr mass 1kg))
            (assembly (sl (s 0)) (instance b block (at) (slide x s)))
            (assembly root {root})"#
        ));
        let n = m.simulated.len();
        (simulate::accelerations(&m, 0.0, &vec![0.0; n], &vec![0.0; n], &vec![0.0; n]).unwrap(), m.coords.len())
    };
    let th = "(def th (simulate 0))";
    let (double, _) = at("(instance p (sl th)) (instance q (sl (* 2 th)))", th);
    near(double[0], 5886.0, 1e-9 * 5886.0, "th and 2 th");
    let (mirror, _) = at("(instance p (sl th)) (instance q (sl (- 0 th)))", th);
    near(mirror[0], 0.0, 1e-9, "th and -th");
    let (sums, _) = at("(instance p (sl th)) (instance q (sl ph)) (instance r (sl (+ th ph)))", "(def th (simulate 0)) (def ph (simulate 0))");
    near(sums[0], 6540.0, 1e-9 * 6540.0, "th beside th + ph");
    near(sums[1], 6540.0, 1e-9 * 6540.0, "ph beside th + ph");
    let (twice, coords) = at("(instance p (sl th)) (instance q (sl th))", th);
    assert_eq!(coords, 1, "th passed twice is one coordinate");
    near(twice[0], 9810.0, 1e-9 * 9810.0, "th passed twice");
}

/// A stiff spring the motion does not ring on (it stopped Dormand-Prince):
/// 1 kg held by a load of -1e5 N/mm times its travel (omega = 1e4 /s, 1.6
/// kHz), pushed from rest by 500 (1 - cos 2 pi t) N. By hand, with r =
/// 2 pi / omega: x = (F0 / 2k) (1 - cos(2 pi t) / (1 - r^2)) + (F0 / 2k)
/// r^2 / (1 - r^2) cos(omega t), its ringing 2e-9 mm, below the tolerance:
/// Dormand-Prince's stiffness test finds it stiff within 0.01 s and Radau
/// IIA steps on, the second in some 760 steps (following what the explicit
/// steps left ringing near the tolerance), where Dormand-Prince alone ran
/// out of its stage budget at 0.234 s.
#[test]
fn a_stiff_spring_the_motion_does_not_ring_on_takes_few_steps() {
    let m = states(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (def u (simulate 0))
        (assembly (axis (s 0)) (instance c block (at) (slide x s)))
        (assembly root (instance a (axis u)))
        (load spring root/a/c (force (* 1N (- (* 500 (- 1 (cos (* 360deg/s t)))) (* 1e5 u))) 0 0))"#,
    );
    let traj = simulated(m, (0.0, 1.0), &[(0.0, 0.0)]);
    let (f, k, w) = (1000.0f64, 1e5f64, 1e4f64);
    let r2 = (2.0 * std::f64::consts::PI / w).powi(2);
    let exact = |t: f64| f / (2.0 * k) * (1.0 - (2.0 * std::f64::consts::PI * t).cos() / (1.0 - r2)) + f / (2.0 * k) * r2 / (1.0 - r2) * (w * t).cos();
    for t in [0.1, 0.25, 0.5, 0.77, 1.0] {
        near(traj.at(0, t).unwrap()[0], exact(t), 2e-7, &format!("x({t}) against the closed form"));
    }
    assert!(traj.t.len() < 1000, "{} steps", traj.t.len());
}

/// A 1 kg slide let go 1 mm out on a 40 Hz load spring rings for a second
/// within its stage budget (before the fix: Radau IIA stepping every
/// simulation ran out at 0.87 s, after 20003 states; Dormand-Prince took
/// 17251): a motion its stiffness test finds no stiffer than its steps is
/// stepped by Dormand-Prince, x(1) = 1 mm.
#[test]
fn a_non_stiff_ringing_costs_no_more_than_it_did() {
    let k = 0.001 * (2.0 * std::f64::consts::PI * 40.0f64).powi(2);
    let m = states(&format!(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (def u (simulate 1))
        (assembly (axis (s 0)) (instance c block (at) (slide x s)))
        (assembly root (instance a (axis u)))
        (load spring root/a/c (force (* -{k}N u) 0 0))"#
    ));
    let traj = simulated(m, (0.0, 1.0), &[(1.0, 0.0)]);
    near(traj.at(0, 1.0).unwrap()[0], 1.0, 1e-5, "x(1)");
    assert!(traj.evaluations <= 17251, "{} states", traj.evaluations);
}

/// A motion stiff for a moment (before the fix: Radau, once it took over,
/// carried on to the end): a 40 Hz ringing (1 kg let go 1 mm out),
/// damped by 100 N s/mm for 10 ms from 0.1 s (1e5 /s), is stepped by Radau
/// there and by Dormand-Prince again once its Jacobian's spectral radius
/// falls back: a second within its stage budget (it stopped at 0.908 s),
/// as many states as without the damper and a few hundred more.
#[test]
fn a_moment_of_stiffness_is_stepped_past() {
    let k = 0.001 * (2.0 * std::f64::consts::PI * 40.0f64).powi(2);
    let run = |c: &str| {
        let m = states(&format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
            (body block (shape (box 50 50 50)) (attr mass 1000g))
            (def u (simulate 1))
            (assembly (axis (s 0)) (instance c block (at) (slide x s)))
            (assembly root (instance a (axis u)))
            (load spring root/a/c (force (* 1N (- (* -{k} u) (* (select (- (/ t 1s) 0.1) 0 (select (- (/ t 1s) 0.11) {c} 0)) (diff u t)))) 0 0))"#
        ));
        simulated(m, (0.0, 1.0), &[(1.0, 0.0)])
    };
    let (plain, damped) = (run("0"), run("100"));
    assert!(damped.evaluations < plain.evaluations + 3000, "{} states against {}", damped.evaluations, plain.evaluations);
    // Damped 10 ms at c / m = 1e5 /s: its speed there (at most 0.25 m/s)
    // is gone, its position held where it was, ringing on from there.
    let x = damped.at(0, 1.0).unwrap()[0];
    assert!(x.abs() < 1.0, "{x}");
}

/// P16's CoreXY braking (y = 60 sin pi t, x still), its belts at 15 lbf and 0.95: the loss the gantry's
/// joint carries is (1/0.95 - 1) of the power it gives the belts, at t = 0.45 s 0.0526 x 1.032 kg x
/// 584 mm/s^2 = 0.0317 N along y. Each of the four idlers riding the gantry takes a share of that, by the
/// size of its pull, and none more: each moved 3.512 N, (1/0.95 - 1) of its own 66.72 N pretension,
/// opposite on its partner, and none once the gantry's power changed sign.
#[test]
fn p16_a_gantry_idler_takes_no_more_than_the_gantrys_loss() {
    let paths = ["root/xy/outA", "root/xy/backA", "root/xy/outB", "root/xy/backB"];
    let build = || {
        let mut m = reacting(
            r#"(import "std") (gravity 0 0 -9.81m/s2)
            (import "../tests/fixtures/dynamics/corexy.ink" as c)
            (assembly root (instance xy (c.corexy 0 (* 60 (sin (* 180deg/s t))))))"#,
            &paths,
        );
        let (mh, mg, mw) = (300.0, 700.0, 8.0);
        m.set_mass("root/xy/head", mh, O, diag(0.0, 0.0, 0.0));
        m.set_mass("root/xy/gantry", mg, O, diag(0.0, 0.0, 0.0));
        for k in ["A", "B"] {
            m.set_mass(&format!("root/xy/motor{k}/rotor"), 30.0, O, diag(0.0, 0.0, 5300.0));
            m.set_mass(&format!("root/xy/pulley{k}"), 8.0, O, diag(0.0, 0.0, 180.0));
            m.set_mass(&format!("root/xy/grub{k}"), 0.5, [5.0, 0.0, 0.0], diag(0.0, 0.0, 0.0));
            for w in ["corner", "foot", "out", "back"] {
                m.set_mass(&format!("root/xy/{w}{k}"), mw, O, diag(0.0, 0.0, 0.0));
            }
            m.set_strip_mass(&format!("root/xy/belt{k}"), 0.0);
        }
        m
    };
    let (lossy, mut ideal) = (build(), build());
    ideal.without_losses();
    let (pi, t) = (std::f64::consts::PI, 0.45);
    let loss = (1.0 / 0.95 - 1.0) * 1.032 * 1e-3 * 60.0 * pi * pi * (pi * t).sin();
    let (l, i) = (lossy.at(t), ideal.at(t));
    for (r, path) in paths.iter().enumerate() {
        let d = l.reactions[r].force[1] - i.reactions[r].force[1];
        assert!(d.abs() <= 1.01 * loss, "{path}: {d} N along y, the gantry's whole loss {loss}");
    }
}

/// A cord between two carts on one rail, both on one coordinate, stating 0.95: its ends move as one along
/// the rail, nothing in it turns or bends, and it loses nothing. A 50 N drag on the far cart, 40 g at 0.2
/// m/s^2: the drive supplies -50 + 0.008 = -49.992 N, and the near cart's joint and the far one's carry
/// what the rail does. (Its loss, 2.5 N, was charged to the drive and put on the near cart alone with no
/// opposite: its joint read -47.492 against the rail's -49.992.)
#[test]
fn a_cord_between_carts_on_one_line_loses_nothing() {
    let m = reacting(
        r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
(body cart (shape (box 40 40 20)) (attr mass 20g))
(assembly (pair (s 0))
  (instance a cart (at) (slide x s) (effort-limits -1000N 1000N))
  (instance b cart (at (move 100 0 0)) (slide x s))
  (relate a b (belt-span (vec (vec 20 0 0) (vec 80 0 0)) 10N 0.95)))
(assembly root (instance p (pair (* 100 (/ t 1s) (/ t 1s)))))
(load drag root/p/b (force 50N 0 0))"#,
        &["root", "root/p/a", "root/p/b"],
    );
    let at = m.at(0.5);
    near(at.effort[0], -49.992, 1e-9, "the drive, nothing lost");
    near(at.reactions[1].force[0] + at.reactions[2].force[0], at.reactions[0].force[0], 1e-9, "the carts' joints sum to the rail's");
}

/// What the masses are known to moves a simulated motion by its first-order
/// sensitivity to each measured moment (simulate/masses.rs), checked against
/// reruns with each moment moved to either side of its proven bounds (half
/// their difference): a double pendulum swung from 50 and -30 degrees (its
/// links' 20 moments, the motion nonlinear and coupled) over 1 s, and a
/// 1 kg block (its mass measured) ringing on 100 N/mm over 25 periods, read
/// at 200 moments: within 13% at each where its band is at least a
/// twentieth of its largest (where a moment's change crosses 0 the two read
/// apart), and 0.4% on the block. One-sided, a rerun at the far side reads
/// the change's second-order part too: up to 1.8 times the band on the double
/// pendulum where a change turns, at an extreme of a link's speed (2026-11-27).
#[test]
fn the_masses_band_is_what_each_measured_moment_moves_the_motion_by() {
    let springs = |src: &str| {
        let doc = Doc::load_str_in(&REGISTRY, src, "validation.ink", Path::new("."), &[]).unwrap_or_else(|e| panic!("{} at {}", e.message, e.line));
        let reload = doc.reloader().unwrap();
        crate::mechanism::springs_model(reload.with_states("root").unwrap(), &reload, "root", gravity(&doc).unwrap()).unwrap_or_else(|e| panic!("{}", e.message))
    };
    let double = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body bar (shape (move (box 20 200 10) 0 100 0)) (attr density 2.7g/cm3) (frame end (at (move 0 200 0))))
        (body tip (shape (move (box 16 150 8) 0 75 0)) (attr density 7.8g/cm3))
        (def a (simulate 50))
        (def b (simulate -30))
        (assembly (double (p 0) (q 0)) (instance bar bar (at (rot-x -90)) (turn z p)) (instance tip tip (frame-of bar end) (turn z q)))
        (assembly root (instance d (double a b)))"#;
    let ringing = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.5s)
        (body b (shape (box 100 100 100)) (attr density 1g/cm3))
        (assembly (axis (c 0)) (instance b b (at) (slide x c) (stiffness 100N/mm)))
        (assembly root (instance a (axis (select (- (/ t 1s) 0.01) 5 0))))
        (springs) (def lag (following-error root/a/b))"#;
    let cases: [(&str, (f64, f64), Vec<(f64, f64)>, usize, f64); 2] = [(double, (0.0, 1.0), vec![(50.0, 0.0), (-30.0, 0.0)], 2, 0.13), (ringing, (0.0, 0.5), Vec::new(), 1, 0.004)];
    for (src, span, start, components, within) in cases {
        let build = |by: Option<(usize, f64)>| {
            let mut m = if start.is_empty() { springs(src) } else { states(src) };
            if let Some((j, by)) = by {
                m.move_source(j, by);
            }
            simulated(m, span, &start)
        };
        let traj = build(None);
        let n = states_or(src, &start, &springs);
        let moved: Vec<[simulate::Trajectory; 2]> = (0..n).map(|j| [build(Some((j, 1.0))), build(Some((j, -1.0)))]).collect();
        let at: Vec<f64> = (1..200).map(|k| span.0 + (span.1 - span.0) * k as f64 / 200.0).collect();
        for i in 0..components {
            for d in 0..2 {
                let band: Vec<f64> = at.iter().map(|&t| traj.mass_band(i, t)[d]).collect();
                let most = band.iter().fold(0.0f64, |a, b| a.max(*b));
                for (&t, &b) in at.iter().zip(&band).filter(|x| *x.1 >= 0.05 * most) {
                    let rerun: f64 = moved.iter().map(|[p, m]| 0.5 * (p.at(i, t).unwrap()[d] - m.at(i, t).unwrap()[d]).abs()).sum();
                    assert!((rerun / b - 1.0).abs() < within, "component {i}'s derivative {d} at {t}: reruns move it {rerun}, its band {b}");
                }
            }
        }
    }
}

/// How many measured moments the model of `src` has (a springs model where it starts from none).
fn states_or(src: &str, start: &[(f64, f64)], springs: &dyn Fn(&str) -> Model) -> usize {
    if start.is_empty() { springs(src) } else { states(src) }.mass_source_count()
}

/// An actuator's limits bounded over a box of states (Model::limits_over)
/// hold every limit read at a state inside it (limits_at): P4's pulley on
/// its DC line and on one giving nothing from 300 to 301 rpm (600 to 602
/// mm/s of the carriage), over speeds either way, across the band, inside
/// it and either side, each box's ends and 400 states in it read; a
/// constant limit reads no state and is its point.
#[test]
fn limits_bounded_over_speeds_hold_every_speed_in_them() {
    use ink_core::interval::I;
    let rate = "(+ 0deg/s (diff (rad->deg (/ x r)) t))";
    let notch = format!("(let ((s (/ (+ 0deg/s (abs {rate})) 1rpm))) (select (- s 300) 400 (select (- s 301) -1 400)))");
    let model = |limits: &str| {
        states(&format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (def r (/ 60 pi))
            (body carriage (shape (box 40 40 20)) (attr density 1g/cm3))
            (assembly (axis (x 0)) (instance c carriage (at) (slide x x))
              (instance m carriage (at (move 0 100 0)) (turn z (rad->deg (/ x r))) (effort-limits {limits})))
            (assembly root (instance a (axis (* 300 (/ t 1s)))))"#
        ))
    };
    let line = format!("(dc-torque-min 1470.997N*mm 330rpm 1 1 {rate}) (dc-torque-max 1470.997N*mm 330rpm 1 1 {rate})");
    for (limits, reads) in [(line, true), (format!("(* 1N*mm (neg {notch})) (* 1N*mm {notch})"), true), ("-400N*mm 400N*mm".to_string(), false)] {
        let m = model(&limits);
        assert_eq!(m.limits_read_state(), [reads], "{limits}");
        for (v0, v1) in [(0.0, 50.0), (-120.0, 80.0), (-700.0, -650.0), (590.0, 599.0), (599.0, 603.0), (600.5, 601.5), (603.0, 700.0), (1000.0, 1100.0)] {
            let state = |v: I| [[I::new(10.0, 12.0), v, I::point(0.0)]];
            let over = m.limits_over(0.5, &state(I::new(v0, v1))).unwrap()[0];
            for k in 0..=400 {
                let v = v0 + (v1 - v0) * k as f64 / 400.0;
                for x in [10.0, 11.0, 12.0] {
                    let (lo, hi) = m.limits_at(0.5, &[[x, v, 0.0]]).unwrap().0[0];
                    assert!(over.0.lo <= lo && lo <= over.0.hi && over.1.lo <= hi && hi <= over.1.hi, "{limits} over [{v0}, {v1}] at {v}: ({lo}, {hi}) outside {over:?}");
                }
            }
            if !reads {
                assert_eq!(over, (I::point(-400.0), I::point(400.0)), "a constant is its point");
            }
        }
        // Across the band its limits cross somewhere: no effort holds over it.
        if limits.contains("select") {
            for (v0, v1) in [(599.0, 603.0), (590.0, 700.0), (-603.0, -599.0)] {
                let over = m.limits_over(0.5, &[[I::point(10.0), I::new(v0, v1), I::point(0.0)]]).unwrap()[0];
                assert!(over.0.hi > over.1.lo, "over [{v0}, {v1}] the band's limits cross: {over:?}");
            }
            let below = m.limits_over(0.5, &[[I::point(10.0), I::new(0.0, 599.0), I::point(0.0)]]).unwrap()[0];
            assert!(below.0.hi <= -399.0 && below.1.lo >= 399.0, "below the band all of it: {below:?}");
        }
    }
}
