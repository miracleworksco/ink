//! Reaction values, bounds, and the bodies and loads they account for.

use super::{asserts_at, claim_value, claims_at_moment, every_claim_at, every_value, load_error, model, near, reacting, rotor, words};
use ink_lang::claims::Status;

/// A still body only a reaction carries is read by that reaction alone: a
/// toleranced plate under root leaves the pendulum's effort claim decided,
/// and the claim on root's reaction, which weighs the plate, unknown.
#[test]
fn a_tolerance_only_a_reaction_reads_leaves_the_effort_decided() {
    let c = every_claim_at(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (param T 6 (tol -0.1 0.1))
        (body bar (shape (move (box 20 300 10) 0 150 0)) (attr density 2.7g/cm3))
        (body plate (shape (move (box 200 200 T) 0 0 -100)) (attr density 2.7g/cm3))
        (assembly (pendulum (angle 0)) (instance bar bar (at (rot-x -90)) (turn z angle) (effort-limits -130N*mm 130N*mm)))
        (assembly root (instance p (pendulum (* 60 (sin (* 180deg/s t))))) (instance plate plate))
        (require held (value (nth (reaction root (at)) 2) > 0))"#,
        &[],
    );
    let of = |n: &str| c.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {c:?}"));
    assert_eq!(of("effort-limits:root/p/bar").1, Status::Holds, "{c:?}");
    assert!(of("held").1 == Status::Unknown && words(&of("held").2).contains(&"T"), "{c:?}");
}

/// An arm swinging beside a round fan spun by t alone: the
/// fan is only a reaction's (root's), so the efforts' steps are those
/// without it, in `ink effort` (every joint's reaction and root's asked,
/// none claimed) and beside a claim on root's reaction; the round fan,
/// read at rest, costs that claim no steps for its spin.
#[test]
fn a_spin_only_a_reaction_carries_leaves_the_efforts_steps_alone() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body link (shape (move (box 150 10 10) 75 0 0)) (attr density 2.7g/cm3))
        (body fan (shape (cylinder-z 40 5)) (attr density 1.2g/cm3))
        (assembly (arm (th 0)) (instance link link (at) (turn y th) (effort-limits -1e6N*mm 1e6N*mm)))
        (assembly root (instance a (arm (* 30 (sin (* 180deg/s t))))) (instance fan fan (at (move 0 200 0)) (turn z (* 360 100 (/ t 1s)))))"#;
    let alone = model(src).steps((0.0, 2.0));
    let cli = reacting(src, &["root", "root/a/link"], false).unwrap().steps((0.0, 2.0));
    assert!(alone.ends == cli.ends && cli.loose.is_empty(), "{} steps alone, {} with the reactions ({} loose)", alone.ends.len() - 1, cli.ends.len() - 1, cli.loose.len());
    let c = every_claim_at(&format!("{src}\n(require up (value (nth (reaction root) 2) > 0))"), &[]);
    assert!(c.iter().all(|x| x.1 == Status::Holds), "{c:?}");
}

/// An unround disc (a 5 mm hole 30 mm off its axis) at 6000 rpm over 2 s
/// its effort claim holds in its own steps beside a claim on
/// its bearings' load; that claim reads the disc at rest, round, within
/// its bound m e w^2 = 73.0 N, and holds under 80 N (its true peak m g +
/// m e w^2 = 79.01 N), and the claim it stays above 0 fails at a
/// moment's real value.
#[test]
fn a_reaction_claim_leaves_an_effort_claim_alone_and_decides_with_its_band() {
    let c = every_claim_at(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body disc (shape (subtract (cylinder-z 50 10) (move (cylinder-z 5 20) 30 0 0))) (attr density 7.85g/cm3))
        (assembly (spinner (a 0)) (instance disc disc (at (rot-x 90)) (turn z a) (effort-limits -1e4N*mm 1e4N*mm)))
        (assembly root (instance s (spinner (* 360 100 (/ t 1s)))))
        (def r (reaction root/s/disc (at)))
        (require fz-under (value (nth r 2) < 80N))
        (require fz-over (value (nth r 2) > 0))"#,
        &[],
    );
    let of = |n: &str| c.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {c:?}"));
    assert_eq!(of("effort-limits:root/s/disc").1, Status::Holds, "{c:?}");
    assert_eq!(of("fz-under").1, Status::Holds, "{c:?}");
    assert_eq!(of("fz-over").1, Status::Fails, "{c:?}");
}

/// Friction written as mu times the normal force: a load reading
/// a reaction, which sums it, is refused, naming the load.
#[test]
fn a_load_reading_a_reaction_is_refused() {
    let e = load_error(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body slider (shape (box 40 40 20)) (attr density 7.85g/cm3))
        (assembly (axis (x 0)) (instance s slider (at) (slide x x) (effort-limits -1e4N 1e4N)))
        (assembly root (instance a (axis (* 100 (sin (* 1rad/s t))))))
        (def n (reaction root/a/s (at)))
        (load friction root/a/s (force (* -0.3 (nth n 2)) 0 0))"#,
    );
    assert!(e.contains("the load friction reads a reaction") && e.contains("(slip R)"), "{e}");
}

/// A payload on the arm assembly, at its forearm's tip: the
/// assembly moves inside, so the load has no one part to ride: refused,
/// naming what moves, not dropped from the efforts.
#[test]
fn a_load_on_an_assembly_moving_inside_is_refused() {
    let e = match reacting(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body upper (shape (move (box 200 10 10) 100 0 0)) (attr density 2.7g/cm3) (frame elbow (at (move 200 0 0))))
        (body fore (shape (move (box 150 10 10) 75 0 0)) (attr density 2.7g/cm3) (frame tip (at (move 150 0 0))))
        (assembly (arm (sh 0) (el 0))
          (instance upper upper (at) (turn y sh) (effort-limits -1e6N*mm 1e6N*mm))
          (instance fore fore (frame-of upper elbow) (turn y el) (effort-limits -1e6N*mm 1e6N*mm)))
        (assembly root (instance a (arm (* 40 (sin (* 2 (/ t 1s)))) (* 30 (cos (* 3 (/ t 1s)))))))
        (load push root/a (force 0 0 -10N) (at (frame-of root/a/fore tip)))"#,
        &[],
        false,
    ) {
        Ok(_) => panic!("built"),
        Err(e) => e.message,
    };
    assert!(e.contains("the load on root/a") && e.contains("root/a/upper, root/a/fore"), "{e}");
}

/// A block moving with a link but not placed on its frames, and one placed
/// on frames of the link and a post: each link's reaction names
/// it as moving with it unseparated; one riding the tip frame is carried.
#[test]
fn a_reaction_names_what_moves_with_it_unseparated() {
    let m = reacting(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body link (shape (move (box 150 10 10) 75 0 0)) (attr density 2.7g/cm3) (frame tip (at (move 150 0 0))))
        (body block (shape (box 20 20 20)) (attr density 7.85g/cm3))
        (body post (shape (move (box 10 10 100) 0 0 50)) (attr density 2.7g/cm3) (frame top (at (move 0 0 100))))
        (assembly (a (th 0)) (instance link link (at) (turn y th) (effort-limits -1e6N*mm 1e6N*mm)) (instance block block (frame-of link tip)))
        (assembly (b (th 0)) (instance link link (at) (turn y th) (effort-limits -1e6N*mm 1e6N*mm)) (instance block block (at (move 150 0 0) (rot-y th))))
        (assembly (c (th 0))
          (instance post post (at (move 150 0 -100)))
          (instance link link (at) (turn y th) (effort-limits -1e6N*mm 1e6N*mm))
          (instance block block (at (move 0 0 (* 0.5 (+ (nth (get (frame-of link tip) "origin") 2) (nth (get (frame-of post top) "origin") 2)))) (frame-of link tip))))
        (assembly root (instance a (a (* 30 (sin (* 1rad/s t))))) (instance b (b (* 30 (sin (* 1rad/s t)))) (at (move 0 300 0))) (instance c (c 0) (at (move 0 600 0))))"#,
        &["root/a/link", "root/b/link", "root/c/link"],
        false,
    )
    .unwrap();
    let with: Vec<&[String]> = m.reactions.iter().map(|r| &r.moves_with[..]).collect();
    assert!(with[0].is_empty() && with[1] == ["root/b/block"] && with[2] == ["root/c/block"], "{with:?}");
}

/// A ground half-space with a density under root: its mass
/// would be the window's, so root's reaction is refused, naming it.
#[test]
fn an_unbounded_body_leaves_a_reaction_unfound() {
    let m = reacting(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body link (shape (move (box 150 10 10) 75 0 0)) (attr density 2.7g/cm3))
        (body floor (shape (+ z 50)) (attr density 1g/cm3))
        (assembly (arm (th 0)) (instance link link (at) (turn y th) (effort-limits -1e6N*mm 1e6N*mm)))
        (assembly root (instance a (arm (* 30 (sin (* 180deg/s t))))) (instance floor floor))"#,
        &["root", "root/a/link"],
        false,
    )
    .unwrap();
    let why = m.reaction_refusal(0, &m.at(0.3)).unwrap_or_default();
    assert!(m.reactions[0].unbounded == ["root/floor"] && why.contains("not proven bounded"), "{why}");
    assert!(m.reaction_refusal(1, &m.at(0.3)).is_none());
}

/// A claim on a root reaction over 12 s, a cart carrying a wheel round but
/// not shown so (four holes at angles computed in floats: turn symmetry of
/// order 1), turned by a phase of the cart's coordinate whose bounds over
/// any step are unbounded (a / E, E = 1 + 100 (sin2 u - (1 - cos2 u)): 1,
/// though its bounds hold 0), as a belt's phase round a moving idler: the
/// wheel's measure puts its centre off the axis by about a millionth of
/// its size. The wheel is a disc, which repeats about its own line and has
/// its centre on it exactly, less four holes, which are measured: what is
/// not known of its imbalance is what is not known of the holes, 0.023 g
/// mm (0.081 where the whole wheel's measure was all that was known). The
/// phase's bounds over the steps the efforts need are wide, and over them
/// that unknown moved the reaction by 3.2 N, past the 1.14 N the claim has
/// to spare (by 1.2 N with the larger unknown, whose efforts split the
/// steps more): a smaller unknown gave a wider band. A claimed reaction's
/// steps are split until what the unknown moves it by is under 4% of its
/// peak, so the claim holds, the unknown carried (its spin is about 9 rad/s:
/// 4.8e-9 kg m times 81 / s^2, 4e-7 N, in truth).
#[test]
fn a_departure_within_its_measure_gives_a_reaction_its_band() {
    let c = every_claim_at(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 12s)
        (body cart (shape (move (box 100 100 10) 0 0 -5)) (attr density 1g/cm3))
        (body wheel (shape (subtract (cylinder-z 20 5) (pattern-polar (move (cylinder-z 2 10) 15 0 0) 4))) (attr density 2.7g/cm3))
        (fn phase (a) (* 3 (/ a (+ 1 (* 100 (- (square (sin (* 0.01 a))) (- 1 (square (cos (* 0.01 a))))))))))
        (assembly (m (a 0)) (instance c cart (at) (slide x a) (effort-limits -1e4N 1e4N)) (instance w wheel (at (move 0 0 50)) (turn z (phase a))))
        (assembly root (instance m (m (* 110 (sin (/ (* pi (/ t 1s)) 2))))))
        (require up (value (nth (reaction root (at)) 2) > 0))"#,
        &[],
    );
    let up = c.iter().find(|x| x.0 == "up").unwrap();
    // What is not known of the holes, as said last: about 0.023 g mm (as each hole's measure is known
    // to), well under the whole wheel's 0.081.
    let unknown: f64 = up.3.rsplit("root/m/w (").next().and_then(|s| s.split(" g mm").next()).and_then(|s| s.parse().ok()).unwrap_or(f64::NAN);
    assert!(up.1 == Status::Holds && (0.01..0.03).contains(&unknown) && up.3.contains("a spinning body it carries moving it by no more than that"), "{unknown} {up:?}");
}

/// A 308 kg steel wheel (250 mm radius, 200 wide: 7.85 x pi 250^2 200 =
/// 308,269 g) spun at 100 rev/s about a level axis, moved 2e-7 mm off it.
/// It is a cylinder, which repeats about its own line, so its centre is on
/// that line exactly, 2e-7 mm off the spin's, whatever its measure reads
/// (0.5 um off, within the 0.1 mm that measure is known to): its imbalance
/// is real and exact, 308.27 kg x 2e-10 m x (200 pi / s)^2 = 0.0243 N,
/// going round with it, and nothing of it is unknown (it read as none,
/// known only to 12,435 N, while the measure was all that was known). So
/// its reaction is its weight, 3024 N, and 0.0243 N across; a claim that
/// what is across stays under 1 N holds. Moved 0.3 mm off, the same:
/// 308.27 x 0.3e-3 x (200 pi)^2 = 36,510 N across at t = 0.
#[test]
fn a_rotor_that_repeats_about_its_own_line_has_its_real_imbalance() {
    let src = |off: &str| {
        format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body wheel (shape (move (cylinder-z 250 200) {off} 0 0)) (attr density 7.85g/cm3))
        (assembly (spinner (a 0)) (instance w wheel (at (rot-x 90)) (turn z a) (effort-limits -1e9N*mm 1e9N*mm)))
        (assembly root (instance s (spinner (* 360 100 (/ t 1s)))))
        (def r (reaction root/s/w (at)))
        (require fx (value (nth r 0) > -1e12N)) (require fz (value (nth r 2) > -1e12N))
        (require across (value (abs (nth r 0)) < 1N)) (require far (value (abs (nth r 0)) < 20000N))"#
        )
    };
    let weight = 308.269 * 9.81;
    let real = 308.269 * 2e-10 * (200.0 * std::f64::consts::PI).powi(2);
    for (t, turned) in [(0.0, 0.0f64), (0.00125, 45.0), (0.0025, 90.0), (0.004, 144.0)] {
        let c = claims_at_moment(&src("2e-7"), t);
        let of = |n: &str| c.iter().find(|x| x.0 == n).unwrap().3;
        near(of("fz") + real * turned.to_radians().sin() * 0.0, weight, 0.02 * weight, "the weight");
        near(of("fx"), -real * turned.to_radians().cos(), 0.02 * real, "its real imbalance, going round");
    }
    let c = every_claim_at(&src("2e-7"), &[]);
    let of = |n: &str| c.iter().find(|x| x.0 == n).unwrap();
    assert!(of("across").1 == Status::Holds && of("far").1 == Status::Holds, "{:?}", of("across"));
    // Nothing of it is unknown, on its bearings' ratings either.
    let m = reacting(&src("2e-7"), &["root/s/w"], true).unwrap_or_else(|e| panic!("{}", e.message));
    let r = m.reactions.iter().position(|r| r.path == "root/s/w").expect("the wheel's joint's reaction");
    let unknown: f64 = m.unknown_at(r, 0.0).iter().map(|s| s.0).sum();
    assert_eq!(unknown, 0.0, "nothing is unknown of a wheel that repeats about its own line");
    let c = claims_at_moment(&src("0.3"), 0.0);
    near(c.iter().find(|x| x.0 == "fx").unwrap().3.abs(), 36510.0, 0.02 * 36510.0, "a real imbalance");
}

/// A pulley that repeats about its axis but for a slot cut in from its rim (a set screw's): a 20 mm
/// disc 10 thick with a 5 mm bore, less a 14 x 3 x 3 slot centred 12 mm out, spun at 100 rev/s. The
/// disc and its bore repeat about the axis, so that part's centre is on it exactly; the slot is
/// measured: 126 mm3 of aluminium, 0.3402 g, 12 mm out, 4.082 g mm, which at this spin is 4.082e-6 kg
/// m x (200 pi / s)^2 = 1.612 N going round. What is not known of the pulley's imbalance is what is
/// not known of the slot alone, a small share of that, where its whole measure left more unknown
/// than the slot is.
#[test]
fn a_pulley_round_but_for_a_slot_has_the_slots_imbalance_and_the_slots_band() {
    let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body pulley (shape (subtract (cylinder-z 20 10) (cylinder-z 2.5 12) (move (box 14 3 3) 12 0 0))) (attr density 2.7g/cm3))
        (assembly (spinner (a 0)) (instance w pulley (at (rot-x 90)) (turn z a) (effort-limits -1e9N*mm 1e9N*mm)))
        (assembly root (instance s (spinner (* 360 100 (/ t 1s)))))
        (def r (reaction root/s/w (at)))
        (require fx (value (nth r 0) > -1e12N))"#;
    let slot = 126.0 * 2.7e-3 * 12.0 * 1e-6 * (200.0 * std::f64::consts::PI).powi(2);
    let c = claims_at_moment(src, 0.0);
    near(c.iter().find(|x| x.0 == "fx").unwrap().3.abs(), slot, 0.03 * slot, "the slot's imbalance");
    let m = reacting(src, &["root/s/w"], true).unwrap_or_else(|e| panic!("{}", e.message));
    let r = m.reactions.iter().position(|r| r.path == "root/s/w").expect("the pulley's joint's reaction");
    let unknown: f64 = m.unknown_at(r, 0.0).iter().map(|s| s.0).sum();
    assert!(unknown < 0.05 * slot, "what is not known of the slot, {unknown} N, beside the slot's {slot} N");
}

const ROTOR_FAST: &str = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
    (body disc (shape (subtract (cylinder-z 50 10) (move (cylinder-z 5 20) 30 0 0))) (attr density 7.85g/cm3))
    (assembly (spinner (a 0)) (instance disc disc (at (rot-x 90)) (turn z a) (effort-limits -1e4N*mm 1e4N*mm)))
    (assembly root (instance s (spinner (* 360 100 (/ t 1s)))))"#;

/// A bearing-load claim checked at one moment: what the
/// provider sets is the real value, so at 0.0625 s, where the hole is on
/// top, fz is m g + m e w^2 = 79.009 N (to the mass's estimate) and the
/// claim it stays under 50 fails; only the
/// stepped path reads the round value with its band.
#[test]
fn a_reaction_claim_at_one_moment_reads_the_real_value() {
    let c = claims_at_moment(&format!("{ROTOR_FAST}\n(require under (value (nth (reaction root/s/disc (at)) 2) < 50N))"), 0.0625);
    let under = c.iter().find(|x| x.0 == "under").unwrap();
    assert!(under.1 == Status::Fails && (under.3 - 79.009).abs() < 0.1, "{under:?}");
}

/// Over the span the band decides: a hold gives its worst case, 79.23 N
/// (not the round 5.988). fx + fz peaks at m g + sqrt 2 m e w^2 = 109 N, the
/// band's box at 152: a claim it stays under 120 is undecided where the
/// real samples stay under, naming each param it reads with its band's
/// half-width (fx and fz within m e w^2 = 73.02, and 0.22 its masses' proven bounds, in N).
#[test]
fn a_banded_claim_says_its_worst_case_and_its_widths() {
    let c = every_claim_at(
        &format!("{ROTOR_FAST}\n(def r (reaction root/s/disc (at)))\n(require fz-under (value (nth r 2) < 80N))\n(require fx-under (value (+ (nth r 0) (nth r 2)) < 120N))"),
        &[],
    );
    let of = |n: &str| c.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {c:?}"));
    assert_eq!(of("fz-under").1, Status::Holds, "{c:?}");
    let why = &of("fx-under").2;
    assert!(
        of("fx-under").1 == Status::Unknown && why.contains("reaction:root/s/disc.fx within 73.2") && why.contains("disc.fz within 73.2") && !why.contains(".mx"),
        "{:?}",
        of("fx-under")
    );
    let at = claims_at_moment(&format!("{ROTOR_FAST}\n(require w (value (nth (reaction root/s/disc (at)) 2) < 80N))"), 0.0625);
    assert_eq!(at[0].1, Status::Holds);
    let x = every_value(&format!("{ROTOR_FAST}\n(require fz-under (value (nth (reaction root/s/disc (at)) 2) < 80N))"), "fz-under");
    // (79.009 from its weight and what departs, and 0.22 its masses' proven bounds.)
    assert!((x - 79.226).abs() < 0.1, "a banded hold's value {x}");
}

/// The tilting moment on the spinning disc's bearings: its
/// departure in the plane across the axis moves no moment in a steady spin
/// on a still frame, so the band is its weight off the axis, h |g| = 1.81 N
/// mm (1550 before): a claim it stays under 100 N mm decides (its measure
/// alone left 105 N mm of its second moment unknown at this spin; the disc
/// repeats about its own line, which leaves that to what is cut from it).
#[test]
fn a_bearing_tilting_moment_claim_decides() {
    let c = every_claim_at(&format!("{ROTOR_FAST}\n(require mx-small (value (abs (nth (reaction root/s/disc) 3)) < 100N*mm))"), &[]);
    let mx = c.iter().find(|x| x.0 == "mx-small").unwrap();
    assert_eq!(mx.1, Status::Holds, "{mx:?}");
}

/// A bracket bolted on a swinging link's tip, its payload placed by the
/// link's angle, not on the bracket's frame: the bracket's
/// reaction (a path with no joint, riding the link's) names the payload as
/// moving with it unseparated; on the bracket's frame it is carried.
#[test]
fn a_path_with_no_joint_names_what_moves_with_it() {
    let m = reacting(
        r#"(import "std") (gravity 0 0 -9.81m/s2)
        (body link (shape (move (box 150 10 10) 75 0 0)) (attr density 2.7g/cm3) (frame tip (at (move 150 0 0))))
        (body bracket (shape (move (box 20 20 20) 10 0 0)) (attr density 2.7g/cm3) (frame end (at (move 20 0 0))))
        (body block (shape (box 20 20 20)) (attr density 7.85g/cm3))
        (assembly (g (th 0)) (instance link link (at) (turn y th) (effort-limits -1e6N*mm 1e6N*mm)) (instance bracket bracket (frame-of link tip))
          (instance block block (at (move 170 0 0) (rot-y th))))
        (assembly (h (th 0)) (instance link link (at) (turn y th) (effort-limits -1e6N*mm 1e6N*mm)) (instance bracket bracket (frame-of link tip))
          (instance block block (frame-of bracket end)))
        (assembly root (instance g (g (* 30 (sin (* 1rad/s t))))) (instance h (h (* 30 (sin (* 1rad/s t)))) (at (move 0 300 0))))"#,
        &["root/g/bracket", "root/h/bracket"],
        false,
    )
    .unwrap();
    assert!(m.reactions[0].moves_with == ["root/g/block"] && m.reactions[1].moves_with.is_empty(), "{:?} {:?}", m.reactions[0].moves_with, m.reactions[1].moves_with);
}

/// A claim reading two reactions: a still bracket's, whose steps
/// are the 16 uniform ones, and a round disc's spun at 6000 rpm read in its
/// own spinning frame, whose components need the finest steps. At the
/// bracket's 16 steps the disc's weight in that frame is sampled along its
/// x, none of it along y (each middle aliases there), and the claim would
/// hold; at every step either reaction chooses, it fails.
#[test]
fn a_claim_reading_two_reactions_takes_the_steps_of_both() {
    let c = every_claim_at(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 2s)
        (body bracket (shape (box 20 20 20)) (attr density 2.7g/cm3))
        (body disc (shape (cylinder-z 40 5)) (attr density 1.2g/cm3))
        (assembly root (instance b bracket (at (move 0 0 100))) (instance d disc (at (rot-x 90) (move 0 200 0)) (turn z (* 360 100 (/ t 1s)))))
        (def a (reaction root/b (at)))
        (def s (reaction root/d (frame-of root/d)))
        (require both (value (+ (nth a 2) (nth s 1)) > 0))"#,
        &[],
    );
    let both = c.iter().find(|x| x.0 == "both").unwrap();
    assert_eq!(both.1, Status::Fails, "{both:?}");
}

/// A door on a hinge whose return spring the door's own assembly declares,
/// a pull across the leaf at its tip growing
/// with the angle, in the door's axes. Placed twice, one turned 90 degrees
/// about z, each carries it as placed: the actuators need, to rounding, the
/// efforts the same two loads written at the top level in the world's axes
/// give, and not the efforts without them.
#[test]
fn a_load_an_assembly_declares_is_placed_with_it() {
    let base = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
        (body leaf (shape (move (box 60 4 80) 30 0 -40)) (attr density 2.7g/cm3) (frame tip (at (move 60 0 -40))))"#;
    let root = r#"(assembly root (instance a (door (* 40 (/ t 1s)))) (instance b (door (* 40 (/ t 1s))) (at (rot-z 90) (move 100 0 0))))"#;
    let door = |spring: &str| format!("(assembly (door (angle 0)) (instance leaf leaf (at (move 0 0 120)) (turn z angle) (effort-limits -1e4N*mm 1e4N*mm)) {spring})");
    let own = model(&format!("{base} {} {root}", door("(load spring leaf (force 0 (* 1N (* -0.5 angle)) 0) (at (frame-of leaf tip)))")));
    let top = model(&format!(
        "{base} {} {root}
        (load sa root/a/leaf (force 0 (* 1N (* -0.5 (* 40 (/ t 1s)))) 0) (at (frame-of root/a/leaf tip)))
        (load sb root/b/leaf (force (* 1N (* 0.5 (* 40 (/ t 1s)))) 0 0) (at (frame-of root/b/leaf tip)))",
        door("")
    ));
    let none = model(&format!("{base} {} {root}", door("")));
    for t in [0.1, 0.5, 0.9] {
        let (a, b, c) = (own.at(t), top.at(t), none.at(t));
        for k in 0..2 {
            near(a.effort[k], b.effort[k], 1e-9 * b.effort[k].abs().max(1.0), "an assembly's load against the same load at the top level");
            assert!((a.effort[k] - c.effort[k]).abs() > 100.0, "the spring moves the effort at {t}: {} {}", a.effort[k], c.effort[k]);
        }
    }
}

#[test]
fn a_slide_carries_axial_force_along_its_declared_axis() {
    let src = r#"(import "std")
        (gravity -10m/s2 0 0)
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (assembly (axis (s 0)) (instance c block (at) (slide x s) (effort-limits -20N 20N)))
        (assembly root (instance a (axis 0)))
        (require axial-rating (value (abs (nth (reaction root/a/c) 0)) < 5N))"#;
    let cs = claims_at_moment(src, 0.0);
    let c = cs.iter().find(|c| c.0 == "axial-rating").unwrap();
    assert!(c.1 == Status::Fails && (c.3 - 10.0).abs() < 1e-6, "{c:?}");
}

#[test]
fn a_bearing_claim_measures_the_component_peak_in_its_frame() {
    let src = r#"(import "std") (gravity 0 0 0) (t-span 0 1s)
(body block (shape (box 50 50 50)) (attr mass 1000g))
(assembly (axis (s 0)) (instance c block (at) (slide x s)))
(assembly root (instance m (axis 0)))
(load push root/m/c (force 0 (* 1N (+ 9 (cos (* 2 pi (- (/ t 1s) 0.25))))) 0))
(require bearing (value (abs (nth (reaction root/m/c (at)) 1)) <= 9.9999N))"#;
    claim_value(src, "bearing", Status::Fails, 10.0, 1e-7);
}

/// The rotor's reaction about its axis is its drive's torque, and reads
/// the same 12.8077: `< 12.806` failed, though it holds.
#[test]
fn a_reaction_within_its_masses_accuracy_of_its_bound_is_not_decided() {
    let r = asserts_at(&rotor("12.806", "(require held (value (nth (reaction root/j/r) 5) < 12.806N*mm))"), &[]).into_iter().find(|x| x.name == "held").unwrap();
    assert_ne!(r.check.status, Status::Fails, "{:?}", r.check.evidence);
}
