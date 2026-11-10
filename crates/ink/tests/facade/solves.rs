//! Derivatives through solves (docs/IR.md, Solving constraints): what a
//! derivative holds fixed, rates over tolerances, rates near a dead point,
//! and the chain rule of derivatives through them.

use ink::claims::{ClaimOptions, Status, run_asserts};
use ink::view::ModelView;
use std::f64::consts::PI;

fn view(src: &str, tag: &str) -> Result<ModelView, ink::InkError> {
    view_at(src, tag, 0.5)
}

/// The model `src` viewed at moment t.
fn view_at(src: &str, tag: &str, t: f64) -> Result<ModelView, ink::InkError> {
    let dir = crate::common::Scratch::new(tag);
    let f = dir.join("m.ink");
    std::fs::write(&f, format!("(import \"std\")\n{src}")).unwrap();
    ModelView::load_in(ink::registry(), &f, &[], t)
}

/// A derivative taken before a solve that makes what it read vary: held
/// fixed there, so the solve is refused, naming where (else it would give
/// 0, exactly, where dw/dt is 200, dw/dp 6, and through u = 3 m then
/// m = 100 t, 600).
/// A solve that leaves it fixed passes the hold on to what it reads. One
/// that never makes it vary loads.
#[test]
fn a_solve_may_not_make_vary_what_a_derivative_held() {
    let speed = "(t-span 0 1) (body block (shape (box 2 2 2))) (assembly root (instance block block (at) (slide x w) (speed-limits -50 50)))";
    for (src, said) in [
        (
            "(param u 0) (solve s2 (unknowns (w 0)) (== w (* 2 u))) (assert w-rate (value (diff w t) < 1)) (solve s1 (unknowns u) (== u (* 100 t)))",
            "solve 's1': it makes 'u' vary along t",
        ),
        (
            "(param p 1) (param q 0) (solve s2 (unknowns (w 0)) (== w (* 2 q))) (assert w-rate (value (diff w p) < 1)) (solve s1 (unknowns q) (== q (* 3 p)))",
            "solve 's1': it makes 'q' vary along 'p'",
        ),
        (
            "(param u 0) (param m 0) (solve s2 (unknowns (w 0)) (== w (* 2 u))) (assert w-rate (value (diff w t) < 1)) (solve s1 (unknowns u) (== u (* 3 m))) (solve s0 (unknowns m) (== m (* 100 t)))",
            "solve 's0': it makes 'u' vary along t",
        ),
        (
            "(param p 1) (param u 0) (param m 0) (solve s2 (unknowns (w 0)) (== w (* 2 u))) (assert dwdp (value (diff w p) < 1)) (solve s1 (unknowns u) (== u (* 3 m))) (solve s0 (unknowns m) (== m (* 5 p)))",
            "solve 's0': it makes 'u' vary along 'p'",
        ),
        (
            "(param u 0) (param m 0) (assert r (value (diff (* 2 u) t) < 1)) (solve s1 (unknowns u) (== u (* 3 m))) (solve s0 (unknowns m) (== m (* 100 t)))",
            "solve 's0': it makes 'u' vary along t",
        ),
        (
            &format!(
                "(param u 0) (param m 0) (solve s2 (unknowns (w 0)) (== w (* 2 u))) (assert w-rate (value (diff w t) < 1)) {speed} (solve s1 (unknowns u) (== u (* 3 m))) (solve s0 (unknowns m) (== m (* 100 t)))"
            ),
            "solve 's0': it makes 'u' vary along t",
        ),
    ] {
        let Err(e) = view(src, "held") else { panic!("{src} loads") };
        assert!(e.message.contains(said) && e.message.contains("held fixed") && e.message.contains("(assert "), "{}", e.message);
    }
    // L is solved after a derivative through it, but does not vary in t;
    // nor along k, which its solve does not read.
    for (src, value) in [
        ("(param L 10) (assert grows (value (diff (* L t) t) > 0)) (solve size (unknowns L) (== (* 2 L) 100))", 50.0),
        ("(param L 10) (param k 3) (assert dk (value (diff (* L k) k) > 0)) (solve size (unknowns L) (== (* 2 L) 100))", 50.0),
    ] {
        let mut v = view(src, "held-ok").unwrap_or_else(|e| panic!("{src}: {}", e.message));
        let c = &run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts[0].check;
        assert!(c.status == Status::Holds && matches!(c.evidence, ink::claims::Evidence::Value(x) if (x - value).abs() < 1e-9), "{src}: {c:?}");
    }
}

/// Over tolerances a solve's rates are certified as its root is: their
/// equations, linear in them, over the root's certified band (benchmark
/// 37's solve, L 50 +-5, h 30 +-0.1: d sx / d L = L / sqrt(L^2 - h^2), from
/// 1.1948 at L = 55, h = 29.9 to 1.3452 at L = 45, h = 30.1). Below 1.26
/// fails at a witness, below 1.35 and above 1.19 hold.
#[test]
fn rates_over_tolerances_are_certified() {
    let mut v = view(
        "(param L 50 (tol -5 5)) (param h 30 (tol -0.1 0.1))
         (solve place (unknowns (theta 37) (sx 40)) (== (* L (cos (* theta (/ pi 180)))) sx) (== (* L (sin (* theta (/ pi 180)))) h))
         (assert rate-mid (value (diff sx L) < 1.26)) (assert rate-lo (value (diff sx L) < 1.35)) (assert rate-hi (value (diff sx L) > 1.19))
         (assert sx-high (value sx <= 40.2))",
        "ranged",
    )
    .unwrap();
    let spans = v.declared_spans().unwrap();
    v.set_param_spans(&spans).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let c: Vec<&ink::claims::Check> = r.asserts.iter().map(|a| &a.check).collect();
    let ink::claims::Evidence::ValueSpan(s) = &c[0].evidence else { panic!("{:?}", c[0]) };
    let at = s.witness.as_ref().map_or(f64::NAN, |w| w.1);
    assert!(c[0].status == Status::Fails && at > 1.26 && at < 1.3452, "{:?}", c[0]);
    assert!(c[1..].iter().all(|c| c.accuracy.is_some()) && c[1].status == Status::Holds && c[2].status == Status::Holds && c[3].status == Status::Fails, "{c:?}");
}

/// Rates are taken at the root shown regular, not where the solver
/// stopped: 50 cos th = 30 + 20 t, followed from t = 0 to 1e-12 short of
/// its dead point at t = 1, th' = -0.4 / sin th (from a dead start, at
/// 1e-13 past it, the solver's point gives 2.3% low).
#[test]
fn rates_near_a_dead_point_are_the_roots() {
    let v = view("(solve loop (unknowns (th 1)) (== (* 50 (cos th)) (+ 30 (* 20 t)))) (assert w (value (diff th t) < 0))", "near-dead").unwrap();
    let t = 1.0 - 1e-12;
    let (ctx, failed) = v.doc.ctx_at(&[], t).unwrap();
    assert!(failed.is_empty());
    let k = v.doc.params.iter().position(|p| p.name == "(diff th t)").unwrap();
    let want = -0.4 / (1.0 - 0.4 * (1.0 - t)).acos().sin();
    assert!((ctx.params[k] / want - 1.0).abs() < 1e-3, "{} against {want}", ctx.params[k]);
}

/// A solve reading what another one on t sets is followed over that one's
/// enclosures over each step: a = 100 t, b^2 = 1 + a from b = -1 (b =
/// -sqrt(1 + 100 t)), c = 2 b. Where 1 + a reaches 0 (t = -0.01) b's root
/// ends: past it b, and c resting on it, are refused.
#[test]
fn a_solve_reading_one_on_t_is_followed_over_its_enclosures() {
    let v = view("(solve A (unknowns (a 0)) (== a (* 100 t))) (solve B (unknowns (b -1)) (== (* b b) (+ 1 a))) (solve C (unknowns (c 1)) (== c (* 2 b)))", "chain").unwrap();
    let k = |n: &str| v.doc.param_index(n).unwrap();
    for t in [0.5, 2.0, -0.005] {
        let (ctx, failed) = v.doc.ctx_at(&[], t).unwrap();
        let b = -(1.0 + 100.0 * t).sqrt();
        assert!(failed.is_empty() && (ctx.params[k("b")] - b).abs() < 1e-9 && (ctx.params[k("c")] - 2.0 * b).abs() < 1e-9, "t = {t}: {failed:?} {:?}", ctx.params);
    }
    let (_, failed) = v.doc.ctx_at(&[], -0.0101).unwrap();
    let names: Vec<&str> = failed.iter().map(|(i, _)| v.doc.ext[v.doc.solves[*i].ext].form.name()).collect();
    assert!(names == ["B", "C"] && failed[0].1.message.contains("dead point"), "{failed:?}");
}

/// A solve reading a rate (`(diff th t)` of a solve on t) is followed over
/// the rate's bounds over each step (its equations, linear in it, bounded
/// with the root's enclosure), so w = th' = -0.4 / sin th, th = acos(0.6 +
/// 0.4 t), both ways from t = 0, not only at 0;
/// past th's dead point at t = 1 it is refused with th.
#[test]
fn a_solve_reading_a_rate_is_followed_over_its_bounds() {
    let v = view("(solve A (unknowns (th 1)) (== (* 50 (cos th)) (+ 30 (* 20 t)))) (solve B (unknowns (w 0)) (== w (diff th t)))", "rate-read").unwrap();
    let w = v.doc.param_index("w").unwrap();
    for t in [0.5, 0.9, -1.0] {
        let (ctx, failed) = v.doc.ctx_at(&[], t).unwrap();
        let want = -0.4 / (0.6 + 0.4 * t).acos().sin();
        assert!(failed.is_empty() && (ctx.params[w] - want).abs() < 1e-9 * want.abs(), "t = {t}: {} against {want}, {failed:?}", ctx.params[w]);
    }
    let (_, failed) = v.doc.ctx_at(&[], 1.5).unwrap();
    assert!(failed.iter().any(|(_, e)| e.message.starts_with("solve 'B' is not followed") && e.message.contains("the rates of solve 'A' along t")), "{failed:?}");
}

/// P3's four-bar with its crank's cosine and sine each the mean of four
/// angles wrapped at one moment (t = 2 pi / 10), in degrees read through
/// `deg->rad` (pi exact): each wraps by exactly a whole turn, so the
/// equations' sides agree exactly there and it is followed across, to the
/// root a continuation of the closed form gives at t = 1. In radians
/// (`(mod (* 10 t) (* 2 pi))`, 2 pi rounded) the sides differ by
/// rounding, a jump the equations really make: refused
/// (a_switch_is_not_called_a_dead_point).
#[test]
fn a_crank_wrapped_by_whole_turns_in_degrees_is_followed() {
    let src = "(param a 100) (param b 400) (param c 300) (param d 350)
        (def ang (rad->deg (* 10 t)))
        (def w1 (mod ang 360)) (def w2 (* 0.5 (mod (* 2 ang) 720))) (def w3 (* 2 (mod (* 0.5 ang) 180))) (def w4 (* 0.25 (mod (* 4 ang) 1440)))
        (fn c4 (f) (* 0.25 (+ (f (deg->rad w1)) (f (deg->rad w2)) (f (deg->rad w3)) (f (deg->rad w4)))))
        (def cs (c4 cos)) (def sn (c4 sin))
        (solve loop (unknowns (th3 0.85) (th4 1.52)) (== (+ (* a cs) (* b (cos th3))) (+ d (* c (cos th4)))) (== (+ (* a sn) (* b (sin th3))) (* c (sin th4))))";
    let v = view(src, "crank").unwrap();
    let k = |n: &str| v.doc.param_index(n).unwrap();
    let (ctx, failed) = v.doc.ctx_at(&[], 1.0).unwrap();
    assert!(failed.is_empty(), "{failed:?}");
    let (th3, th4) = (ctx.params[k("th3")], ctx.params[k("th4")]);
    assert!((th3 - 0.852406022535).abs() < 1e-9 && (th4 - 2.175833101771).abs() < 1e-9, "th3 {th3}, th4 {th4}");
}

/// Following stops where no step is certified. Where the equations are
/// smooth that is a dead point; where they switch and their bounds do not
/// show them continuous (a switch in the unknowns: u mod 2 pi; four floors
/// of t switching at one moment), it says so, not that it is a dead
/// point. The four-bar is continuous (the four angles wrap together
/// inside cos and sin), but nothing shows its sides agree at the seam.
#[test]
fn a_switch_is_not_called_a_dead_point() {
    let four = "(param a 100) (param b 400) (param c 300) (param d 350)
        (def w1 (mod (* 10 t) (* 2 pi))) (def w2 (* 0.5 (mod (* 20 t) (* 4 pi)))) (def w3 (* 2 (mod (* 5 t) pi))) (def w4 (* 0.25 (mod (* 40 t) (* 8 pi))))
        (def cs (* 0.25 (+ (cos w1) (cos w2) (cos w3) (cos w4)))) (def sn (* 0.25 (+ (sin w1) (sin w2) (sin w3) (sin w4))))
        (solve loop (unknowns (th3 0.85) (th4 1.52)) (== (+ (* a cs) (* b (cos th3))) (+ d (* c (cos th4)))) (== (+ (* a sn) (* b (sin th3))) (* c (sin th4))))";
    for (src, t, cause) in [
        ("(solve s (unknowns (u 1)) (== (mod u (* 2 pi)) (mod (+ 1 t) (* 2 pi))))", 5.5, "its equations switch there"),
        (four, 1.0, "its equations switch there"),
        ("(solve s (unknowns (u 1)) (== (square u) (- 1 t)))", 1.5, "a dead point"),
    ] {
        let v = view(src, "switch").unwrap();
        let (_, failed) = v.doc.ctx_at(&[], t).unwrap();
        let e = &failed.first().unwrap_or_else(|| panic!("{src}: followed to {t}")).1.message;
        assert!(e.contains("is not followed") && e.contains(cause) && !(cause != "a dead point" && e.contains("dead point")), "{e}");
    }
}

/// A box holding one root at every moment of a step is the root followed
/// only where the equations do not jump: u^2 = 0.9 + mean of cos(0.75 w +
/// pi / 3) over angles w wrapped together at t = 0.01 jumps from u =
/// 1.3289 to 1.1832 there, and one box holds both, with four angles (the
/// bounds take each wrap as a choice of its sides) and with three (bounded
/// case by case), so it would cross the jump as one root. Following needs
/// the equations' rate along t over the step bounded
/// (Roots::continuous): refused past the jump, followed up to it.
#[test]
fn a_root_is_not_followed_across_a_jump() {
    let wraps = ["(mod a (* 2 pi))", "(* 0.5 (mod (* 2 a) (* 4 pi)))", "(* 2 (mod (* 0.5 a) pi))", "(* 0.25 (mod (* 4 a) (* 8 pi)))"];
    for n in [3, 4] {
        let terms: String = wraps[..n].iter().map(|w| format!(" (c {w})")).collect();
        let src = format!(
            "(def a (+ (- (* 2 pi) 0.1) (* 10 t))) (fn c (w) (cos (+ (* 0.75 w) (/ pi 3))))
            (def rhs (+ 0.9 (* (/ 1 {n}) (+{terms})))) (solve s (unknowns (u 1.31)) (== (square u) rhs))"
        );
        let v = view(&src, "jump").unwrap();
        let (ctx, failed) = v.doc.ctx_at(&[], 0.0099).unwrap();
        let u = ctx.params[v.doc.param_index("u").unwrap()];
        let left = (0.9 + (0.75 * (2.0 * PI - 0.1 + 0.099) + PI / 3.0).cos()).sqrt();
        assert!(failed.is_empty() && (u - left).abs() < 1e-9, "{n} wraps: u {u} at t 0.0099 against {left}, {failed:?}");
        let (_, failed) = v.doc.ctx_at(&[], 0.02).unwrap();
        let e = &failed.first().unwrap_or_else(|| panic!("{n} wraps: followed across the jump to 0.02")).1.message;
        assert!(e.contains("is not followed") && e.contains("its equations switch there"), "{n} wraps: {e}");
    }
}

/// A moment is a number: at t = NaN a solve on t has no root (evaluation
/// is total, so its equations have values there, but none is a root).
#[test]
fn a_solve_on_t_has_no_root_at_nan() {
    let v = view("(solve s (unknowns (u 1)) (== u (* 3 t))) (solve fixed (unknowns (k 1)) (== k 2))", "nan").unwrap();
    let (_, failed) = v.doc.ctx_at(&[], f64::NAN).unwrap();
    assert!(failed.len() == 1 && failed[0].1.message.contains("solve 's' reads t, and t = NaN is not a moment"), "{failed:?}");
}

/// Twin roots of a solve reading a rate, 1e-3 apart (v = u' or u' + 1e-3,
/// u^2 = (t - 1)^2 + 1e-4): over a step the rate is bounded along its
/// root's path, by how far it moves, so the steps are as long as that
/// allows: 84 to t = 0.5, where bounded over its root's box (t and the root
/// ranged apart) it would take 12386. v is u', the root followed.
#[test]
fn a_rate_read_over_a_step_is_bounded_by_how_far_it_moves() {
    let certified = || ink::lang::follow::CERTIFIED.with(std::cell::Cell::get);
    let before = certified();
    let v =
        view("(solve A (unknowns (u 1)) (== (square u) (+ (square (- t 1)) 1e-4))) (solve B (unknowns (v -1)) (== (* (- v (diff u t)) (- v (diff u t) 1e-3)) 0))", "twin").unwrap();
    let steps = certified() - before;
    let (u, w) = ((0.25f64 + 1e-4).sqrt(), v.ctx.params[v.doc.param_index("v").unwrap()]);
    assert!((w + 0.5 / u).abs() < 1e-9 && steps < 1000, "v {w} after {steps} steps");
}

/// A solve whose equations repeat in t is followed one period, and every
/// other moment folded into it (ink_core::period): sin(th - 10 t) = 1/2
/// from th = 0.5 is th = 10 t + pi / 6, a turn a period, at t = 1e7 and
/// -1e7 (16 M periods each way, which step by step would take hours).
/// tests/motion/solve-branch.ink's four-bar at t = 628.3, 1000 turns: the
/// rocker where following every step put it (11.7 s cold).
#[test]
fn a_solve_that_repeats_in_t_is_followed_one_period() {
    let v = view("(solve A (unknowns (th 0.5)) (== (sin (- th (* 10 t))) 0.5))", "period").unwrap();
    let th = v.doc.param_index("th").unwrap();
    for t in [1e7, -1e7, 0.3] {
        let (ctx, failed) = v.doc.ctx_at(&[], t).unwrap();
        let want = 10.0 * t + PI / 6.0;
        assert!(failed.is_empty() && (ctx.params[th] - want).abs() < 1e-6, "t = {t}: th {} against {want}, {failed:?}", ctx.params[th]);
    }
    let src = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/motion/solve-branch.ink")).unwrap();
    let v = view(&src.replace("(import \"std\")", ""), "period-branch").unwrap();
    // Its crank turns in degrees (deg->rad, pi exact): a period of exactly
    // 36 degrees of t's 10 rad / s. At 628.3 + 1e6 periods, 1e6 turns more
    // (step by step, hours), the same rocker.
    for t in [628.3, 628.3 + 1e6 * 2.0 * PI / 10.0] {
        let (ctx, failed) = v.doc.ctx_at(&[], t).unwrap();
        let th4 = ctx.params[v.doc.param_index("th4").unwrap()];
        assert!(failed.is_empty() && (th4 - 1.6022386532214294).abs() < 1e-6, "t = {t}: th4 {th4}, {failed:?}");
    }
}

/// An unknown with no period of its own is not moved over t's period: q
/// = sin t repeats with q held, shown so where its root at 0 lands in the
/// tube at 2 pi. At t = 1e12 (followed step by step, about a step per
/// 1000 s: a billion) it is sin(1e12).
#[test]
fn a_root_that_repeats_unmoved_is_folded() {
    let v = view("(solve s (unknowns (q 0)) (== q (sin t)))", "period-unmoved").unwrap();
    let q = v.doc.param_index("q").unwrap();
    let (ctx, failed) = v.doc.ctx_at(&[], 1e12).unwrap();
    assert!(failed.is_empty() && (ctx.params[q] - 1e12f64.sin()).abs() < 1e-6, "q {} against {}, {failed:?}", ctx.params[q], 1e12f64.sin());
}

/// Folding is a shortcut, never a detour: q = sin 0.3 t + sin 0.7 t
/// repeats only every 1.13e17 s (0.3 and 0.7 as floats), so a moment
/// nearer 0 than that, t = -0.1 as t = 0.1, is followed straight there.
#[test]
fn a_moment_nearer_than_a_period_is_not_folded() {
    let v = view("(solve s (unknowns (q 0)) (== q (+ (sin (* 0.3 t)) (sin (* 0.7 t)))))", "period-long").unwrap();
    let q = v.doc.param_index("q").unwrap();
    for t in [0.1, -0.1] {
        let (ctx, failed) = v.doc.ctx_at(&[], t).unwrap();
        let want = (0.3 * t).sin() + (0.7 * t).sin();
        assert!(failed.is_empty() && (ctx.params[q] - want).abs() < 1e-9, "t = {t}: q {} against {want}, {failed:?}", ctx.params[q]);
    }
}

/// A moment so far that its count of periods is uncertain (t = 1e20, the
/// period 2 pi enclosed to a float: about 1e4 periods either way) is one
/// of the whole first period's moments, its root certified at t itself:
/// q = sin t is sin(1e20), never followed over all the time between.
#[test]
fn a_moment_whose_count_of_periods_is_uncertain_is_folded() {
    let v = view("(solve s (unknowns (q 0)) (== q (sin t)))", "period-far").unwrap();
    let q = v.doc.param_index("q").unwrap();
    for t in [1e20, -1e20] {
        let (ctx, failed) = v.doc.ctx_at(&[], t).unwrap();
        assert!(failed.is_empty() && (ctx.params[q] - t.sin()).abs() < 1e-9, "t = {t}: q {} against {}, {failed:?}", ctx.params[q], t.sin());
    }
}

/// Steps grow as long as they are certified, whatever the time scale: q =
/// sin(1e-20 t) at t = 3e20 (its period 6.3e20 s) is sin 3 after about 70
/// steps, as q = sin t at t = 3 is after a few.
#[test]
fn steps_grow_with_how_slowly_the_equations_change() {
    let v = view("(solve s (unknowns (q 0)) (== q (sin (* 1e-20 t))))", "period-slow").unwrap();
    let q = v.doc.param_index("q").unwrap();
    let (ctx, failed) = v.doc.ctx_at(&[], 3e20).unwrap();
    assert!(failed.is_empty() && (ctx.params[q] - 3f64.sin()).abs() < 1e-9, "q {} against {}, {failed:?}", ctx.params[q], 3f64.sin());
}

/// A solve on t only through others sets the root followed to the moment,
/// whatever moment is asked first, not its root at load time, t = 0, which
/// its reads repeating their bits there would hand back. q = v', v = u',
/// u^2 = (t - 1)^2 + 0.05^2, symmetric about t = 1: the load-time root of
/// q(2) asked first is 1e-13 off the one followed. th, the direction of
/// (a', b'), a loop about 0 that is (3, 0.5) at t = 0 and 1 exactly: its
/// load-time root of th(1) asked first is a turn short.
#[test]
fn a_solve_on_t_through_others_is_followed_whatever_is_asked_first() {
    let chain =
        "(solve A (unknowns (u 1)) (== (square u) (+ (square (- t 1)) (square 0.05)))) (solve B (unknowns (v -1)) (== v (diff u t))) (solve C (unknowns (q 0)) (== q (diff v t)))";
    let winds = "(solve A (unknowns (a 0) (b 0)) (== a (+ (* 3 t) (* -12 t t) (* 8 t t t))) (== b (+ (* 0.5 t) (* 8 t t) (* -16 t t t) (* 8 t t t t))))
        (solve B (unknowns (th 0.2)) (== (* (cos th) (diff b t)) (* (sin th) (diff a t))))";
    let e2 = 0.05f64 * 0.05;
    for (src, t, name, want) in [(chain, 2.0, "q", e2 / (1.0 + e2).powf(1.5)), (winds, 1.0, "th", 0.5f64.atan2(3.0) + 2.0 * std::f64::consts::PI)] {
        let v = view_at(src, name, t).unwrap();
        let k = v.doc.param_index(name).unwrap();
        let first = v.ctx.params[k];
        v.doc.ctx_at(&[], 0.5).unwrap();
        let (ctx, failed) = v.doc.ctx_at(&[], t).unwrap();
        let again = ctx.params[k];
        assert!(
            failed.is_empty() && first.to_bits() == again.to_bits() && (first - want).abs() < 1e-9,
            "{name}({t}): {first} asked first, {again} after t = 0.5, {want} its closed form"
        );
    }
}

/// A claim's terms at its witness take what a solve on t sets at the
/// witness's moment, one on t only through another too: b = a, a = 100 t,
/// drives a slide limited to 75 over t in [0, 1]: the limit fails past
/// t = 0.75, by the drive the witness gives less 75, not the view's
/// moment's (50, t = 0.5), which is inside it.
#[test]
fn a_witness_takes_a_solve_on_t_through_others_at_its_moment() {
    let mut v = view("(t-span 0 1) (solve A (unknowns (a 0)) (== a (* 100 t))) (solve B (unknowns (b 0)) (== b a)) (body block (shape (box 2 2 2))) (assembly root (instance block block (at) (slide x b) (limits -1000 75)))", "witness-t").unwrap();
    v.set_motion(Some((0.0, 1.0))).unwrap();
    let a = &run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts[0];
    let margin = match &a.check.evidence {
        ink::claims::Evidence::Value(m) => *m,
        ink::claims::Evidence::ValueSpan(s) => s.witness.as_ref().expect("a witness").1,
        e => panic!("{e:?}"),
    };
    let drive = a.terms.iter().find(|t| t.0 == "drive").expect("a drive").1;
    assert!(a.check.status == Status::Fails && drive > 75.0 && (drive - (75.0 - margin)).abs() < 1e-9, "drive {drive}, margin {margin}");
}

/// A claim over a motion reads what a solve on t sets through t, its rate
/// along t bounded over each piece: a four-bar's loop closes to 0.001 over
/// a whole turn of its crank, its gap's bound second order in a piece of t
/// (its rate along the roots' path is about 0). Read as boxes of their own,
/// the two angles blur the gap by their width times the links' lengths:
/// 595,539 pieces of t.
#[test]
fn a_loop_closed_by_a_solve_on_t_is_proven_along_its_roots() {
    let src = "(t-span 0 1) (param crank 20) (param coupler 60) (param rocker 40) (param ground 50) (def a (* 2 pi t))
        (solve loop (unknowns (u 0.8) (w 1.6))
          (== (+ (* crank (cos a)) (* coupler (cos u))) (+ ground (* rocker (cos w))))
          (== (+ (* crank (sin a)) (* coupler (sin u))) (* rocker (sin w))))
        (def gap (distance (vec (+ (* crank (cos a)) (* coupler (cos u))) (+ (* crank (sin a)) (* coupler (sin u)))) (vec (+ ground (* rocker (cos w))) (* rocker (sin w)))))
        (assert closes (value gap < 0.001))";
    let mut v = view(src, "loop-along").unwrap();
    v.set_motion(Some((0.0, 1.0))).unwrap();
    let c = &run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts[0].check;
    let ink::claims::Evidence::ValueSpan(s) = &c.evidence else { panic!("{c:?}") };
    assert!(c.status == Status::Holds && c.accuracy == Some(0.0) && s.pieces < 20_000, "{c:?}");
}

/// So over a tolerance band: a solve's unknowns over a ranged param are
/// read through it at their rates bounded over each piece (-J_u^-1 J_p),
/// and the four-bar's loop closes to 0.001 over its crank's band in 7
/// pieces (2,047 read as boxes of their own; with all four links
/// toleranced, 262,143 pieces in 7 s, none in minutes before).
#[test]
fn a_loop_closed_by_a_solve_over_a_band_is_proven_along_its_roots() {
    let src = "(param crank 20 (tol -0.2 0.2)) (param coupler 60) (param rocker 40) (param ground 50) (def a 0.7)
        (solve loop (unknowns (u 0.8) (w 1.6))
          (== (+ (* crank (cos a)) (* coupler (cos u))) (+ ground (* rocker (cos w))))
          (== (+ (* crank (sin a)) (* coupler (sin u))) (* rocker (sin w))))
        (def gap (distance (vec (+ (* crank (cos a)) (* coupler (cos u))) (+ (* crank (sin a)) (* coupler (sin u)))) (vec (+ ground (* rocker (cos w))) (* rocker (sin w)))))
        (assert closes (value gap < 0.001))";
    let mut v = view(src, "loop-band").unwrap();
    let spans = v.declared_spans().unwrap();
    v.set_param_spans(&spans).unwrap();
    let c = &run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts[0].check;
    let ink::claims::Evidence::ValueSpan(s) = &c.evidence else { panic!("{c:?}") };
    assert!(c.status == Status::Holds && c.accuracy == Some(0.0) && s.pieces < 100, "{c:?}");
}

/// A rate an upstream solve leaves unbounded over a step is no rate of 0 to
/// what reads it: u = t (its twin u = t - 0.1 within 0.05 where its
/// Jacobian holds 0 over a step), w = u, q = w; not q = 0.6875 at t = 0.9,
/// its residual 0.2125, with no word.
#[test]
fn an_unbounded_rate_is_unbounded_downstream() {
    let v = view("(solve A (unknowns (u 0)) (== (* (- u t) (+ (- u t) 0.1)) 0)) (solve B (unknowns (w 0)) (== (- w u) 0)) (solve C (unknowns (q 0)) (== (- q w) 0))", "chain3")
        .unwrap();
    let q = v.doc.param_index("q").unwrap();
    for t in [0.03, 0.1, 0.9, 5.0] {
        let (ctx, failed) = v.doc.ctx_at(&[], t).unwrap();
        assert!(failed.is_empty() && (ctx.params[q] - t).abs() < 1e-9, "t = {t}: q {} {failed:?}", ctx.params[q]);
    }
}

/// u = 0 and u = 1/8 have no common root, but their sum weighted by the
/// Jacobian's left inverse (1/2, 1/2) has one, 1/16: Krawczyk's test on it
/// would certify that as a root of both (Roots::enclose, a followed
/// solve's tube). A root is certified to exist only for as many equations
/// as unknowns.
#[test]
fn more_equations_than_unknowns_have_no_certified_root() {
    use ink::graph::{Graph, Var};
    use ink::lang::roots::Roots;
    use ink::tape::Tape;
    let mut g = Graph::new();
    let u = g.var(Var::Param(0));
    let eighth = g.constant(0.125);
    let other = g.sub(u, eighth);
    let tapes = || vec![Tape::new(&g, u), Tape::new(&g, other)];
    let ctx = ink::eval::Ctx::new(0.0, vec![0.0]);
    let refused = |r: Result<(), String>| r.is_err_and(|e| e.contains("2 equations for 1 unknowns"));
    assert!(refused(Roots::enclose(&[0], tapes(), &ctx, None, &[]).map(drop)));
    assert!(refused(Roots::certify(&[0], tapes(), false, &ctx, |_, _| Ok(())).map(drop)));
    let square = Roots::enclose(&[0], vec![Tape::new(&g, other)], &ctx, None, &[]).unwrap();
    assert!(square.1[0].lo <= 0.125 && 0.125 <= square.1[0].hi, "{square:?}");
}

/// A solve on t reading one with more equations than unknowns (measured at
/// each moment, no tube) is measured too, not followed asking a's tube
/// over each step: b = max(a, 100), a = t twice. At a moment and over the
/// motion, its claim holds, an estimate.
#[test]
fn a_solve_reading_a_measured_one_is_measured() {
    let src = "(solve A (unknowns (a 0)) (== a t) (== a t)) (solve B (unknowns (b 0)) (== b (max a 100))) (assert b-ok (value (abs (- b 100)) < 0.001))";
    let mut v = view_at(src, "measured-chain", 0.5).unwrap();
    let estimate = |v: &mut ModelView| {
        let a = &run_asserts(v, None, false, ClaimOptions::default()).unwrap().asserts[0];
        assert!(a.check.status == Status::Holds && a.check.accuracy.is_none(), "{a:?}");
    };
    estimate(&mut v);
    v.set_motion(Some((0.0, 1.0))).unwrap();
    estimate(&mut v);
    let b = v.doc.param_index("b").unwrap();
    for t in [0.25, 0.5, 3.0] {
        let (ctx, failed) = v.doc.ctx_at(&[], t).unwrap();
        assert!(failed.is_empty() && (ctx.params[b] - 100.0).abs() < 1e-9, "t = {t}: b {} {failed:?}", ctx.params[b]);
    }
}
