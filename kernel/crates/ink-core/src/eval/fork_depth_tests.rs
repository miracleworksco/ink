use super::{Ctx, Inputs, eval_interval, eval_interval_grad, eval_interval_rates, eval_point, forked_count};
use crate::graph::{Graph, UnOp, Var};
use crate::interval::I;
use crate::tape::Tape;

/// Nine jumps a box straddles once each, every one read again at the output through |.| (so each
/// one's window holds every jump after it): forked each inside the last, they bounded the tape
/// 2^9 times over. Past FORK_DEPTH a jump is bounded across, as a wider bound: a handful of forks,
/// and the bound still holds every point's value.
#[test]
fn nested_jumps_fork_only_so_deep() {
    let mut g = Graph::new();
    let x = g.var(Var::X);
    let mut sum = g.constant(0.0);
    for k in 1..10 {
        let kx = g.constant(k as f64);
        let scaled = g.mul(x, kx);
        let step = g.un(UnOp::Floor, scaled);
        sum = g.add(sum, step);
    }
    let shifted = g.sub(sum, x);
    let out = g.un(UnOp::Abs, shifted);
    let tape = Tape::new(&g, out);
    let ctx = Ctx::new(0.0, vec![]);
    let b = [I::new(0.97, 1.03), I::point(0.0), I::point(0.0)];
    let before = forked_count();
    let bound = eval_interval(&tape, &ctx, b, None);
    let forks = forked_count() - before;
    assert!(forks <= 64, "{forks} forks");
    for i in 0..=60 {
        let p = 0.97 + 0.001 * i as f64;
        let v = eval_point(&tape, &ctx, [p, 0.0, 0.0]);
        assert!(bound.lo <= v && v <= bound.hi, "{v} at {p} outside {bound:?}");
    }
}

#[test]
fn corner_reach_reading_time_or_moving_frame_reports_unbounded_rate() {
    // Sub-tape reading `t` directly with constant r, guard, and identity frame:
    // a = x - t, b = y, r = 1, guard = -1 at (5, 0, 0) over t in [0, 10].
    let mut g = Graph::new();
    let (x, y, t) = (g.var(Var::X), g.var(Var::Y), g.var(Var::T));
    let a = g.sub(x, t);
    let r = g.constant(1.0);
    let guard = g.constant(-1.0);
    let cr = g.corner_reach(a, y, r, guard);
    let zero = g.constant(0.0);
    let big = g.constant(1e9);
    let sel = g.select(cr, zero, big);

    let t_cr = Tape::new(&g, cr);
    let t_sel = Tape::new(&g, sel);
    let mut ctx = Ctx::new(5.0, vec![]);
    ctx.t_span = Some((0.0, 10.0));
    let b = [I::point(5.0), I::point(0.0), I::point(0.0)];

    // Inputs::of must discover `t` inside CornerReach's sub-tapes.
    let inputs_cr = Inputs::of([&t_cr], &ctx);
    assert!(inputs_cr.t, "Inputs::of missed t inside CornerReach");
    let (v_cr, d_cr) = eval_interval_rates(&t_cr, &ctx, b, &inputs_cr);
    assert!(!(d_cr[0].lo.is_finite() && d_cr[0].hi.is_finite()), "expected unbounded rate for CornerReach, got {d_cr:?}");
    for tv in [0.0, 5.0, 10.0] {
        let pt = eval_point(&t_cr, &Ctx::new(tv, vec![]), [5.0, 0.0, 0.0]);
        assert!(v_cr.lo <= pt && pt <= v_cr.hi, "{pt} at t={tv} outside {v_cr:?}");
    }

    // Driving a select: value jumps 0 -> 1e9 across t in [0, 10].
    let inputs_sel = Inputs::of([&t_sel], &ctx);
    assert!(inputs_sel.t, "Inputs::of missed t through select(CornerReach)");
    let (v_sel, d_sel) = eval_interval_rates(&t_sel, &ctx, b, &inputs_sel);
    assert!(v_sel.lo <= 0.0 && v_sel.hi >= 1e9, "select bound missed jump: {v_sel:?}");
    assert!(!(d_sel[0].lo.is_finite() && d_sel[0].hi.is_finite()), "expected unbounded rate for select(CornerReach), got {d_sel:?}");

    // Moving frame with static sub-tapes: frame coordinate x - t moves with t.
    let mut g2 = Graph::new();
    let (x2, y2, z2, t2) = (g2.var(Var::X), g2.var(Var::Y), g2.var(Var::Z), g2.var(Var::T));
    let r2 = g2.constant(1.0);
    let guard2 = g2.constant(-1.0);
    let cr0 = g2.corner_reach(x2, y2, r2, guard2);
    let shifted = g2.sub(x2, t2);
    let cr_moved = g2.subst(cr0, &[(Var::X, shifted), (Var::Y, y2), (Var::Z, z2)]);
    let t_moved = Tape::new(&g2, cr_moved);
    let inputs_m = Inputs::of([&t_moved], &ctx);
    assert!(inputs_m.t);
    let (_vm, dm) = eval_interval_rates(&t_moved, &ctx, b, &inputs_m);
    assert!(!(dm[0].lo.is_finite() && dm[0].hi.is_finite()), "expected unbounded rate for moving-frame CornerReach, got {dm:?}");
}

#[test]
fn corner_reach_across_floor_seam_is_not_marked_alike_when_frame_jumps() {
    // Frame x = x - floor(x) jumps from 1 to 0 across the floor seam at x = 0,
    // while r and guard are constants (Alike). The output jumps across x = 0,
    // so its x-derivative over a box straddling x = 0 must be unbounded.
    let mut g = Graph::new();
    let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
    let r = g.constant(1.0);
    let guard = g.constant(-1.0);
    let cr0 = g.corner_reach(x, y, r, guard);
    let fx = g.un(UnOp::Floor, x);
    let xw = g.sub(x, fx);
    let cr_wrapped = g.subst(cr0, &[(Var::X, xw), (Var::Y, y), (Var::Z, z)]);
    let tape = Tape::new(&g, cr_wrapped);
    let ctx = Ctx::new(0.0, vec![]);
    let left = eval_point(&tape, &ctx, [-1e-6, 0.25, 0.0]);
    let right = eval_point(&tape, &ctx, [1e-6, 0.25, 0.0]);
    assert!((left - right).abs() > 0.5, "expected jump across x=0 seam: left={left}, right={right}");
    let grad = eval_interval_grad(&tape, &ctx, [I::new(-0.1, 0.1), I::point(0.25), I::point(0.0)]);
    assert!(!(grad.d[0].lo.is_finite() && grad.d[0].hi.is_finite()), "expected unbounded dx across floor seam, got {:?}", grad.d[0]);
}
