use super::bounds::{Found, bounds_cells, extent_now, material_parts};
use super::clearance::clearance_now;
use super::empty::{Empty, Kept, empty_now, lattice};
use super::proof::WORK;
use super::ray::segment_box;
use super::*;
use crate::octree::{Cell, Local, try_root_for};
use crate::search;
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_point, span_mid};
use ink_core::graph::{BinOp, Graph, NodeId, UnOp, Var};
use ink_core::interval::I;
use ink_core::tape::{BOTH, Tape};
use std::sync::Arc;

/// How a nut sits on a screw in `seated`: turned about z by `turn`
/// degrees, its thread of pitch `pitch`, left-hand where `left`.
#[derive(Clone, Copy)]
struct Nut {
    turn: f64,
    pitch: f64,
    left: bool,
}

/// The overlap of an M3 screw (its thread clipped to 8 mm under its
/// head) and a 4 mm nut whose tapped thread (std's thread-iso, both,
/// as the catalog's: 2.94 and 3.08 across) is seated 5.5 mm down its
/// axis (11 turns of 0.5: in phase), placed as `nut` says; its
/// emptiness proven at accuracy `tol`, and the boxes that took.
fn seated(nut: Nut, tol: f64) -> (Proof, usize) {
    seated_at(nut, tol, None, &Ctx::new(0.0, vec![]))
}

/// `seated`, the nut moved up its axis by `lift` (param 0 or t) more,
/// asked in `ctx`.
fn seated_at(nut: Nut, tol: f64, lift: Option<Var>, ctx: &Ctx) -> (Proof, usize) {
    let mut b = crate::test_fields::B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let (xx, yy) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y));
    let rr = b.bin(BinOp::Add, xx, yy);
    let r = b.un(UnOp::Sqrt, rr);
    let cylinder = |b: &mut crate::test_fields::B, radius: f64, lo: f64, hi: f64| {
        let k = b.c(radius);
        let wall = b.bin(BinOp::Sub, r, k);
        let ends = b.slab(z, 0.5 * (lo + hi), 0.5 * (hi - lo));
        b.bin(BinOp::Max, wall, ends)
    };
    let thread = b.thread_iso([x, y, z], 2.94, 0.5);
    let shank = b.slab(z, -4.0, 4.0);
    let shank = b.bin(BinOp::Max, thread, shank);
    let head = cylinder(&mut b, 2.75, 0.0, 3.0);
    let screw = b.bin(BinOp::Min, shank, head);
    let (c, s) = (b.c(nut.turn.to_radians().cos()), b.c(nut.turn.to_radians().sin()));
    let (cx, sy) = (b.bin(BinOp::Mul, c, x), b.bin(BinOp::Mul, s, y));
    let u = b.bin(BinOp::Add, cx, sy);
    let ms = b.un(UnOp::Neg, s);
    let (sx, cy) = (b.bin(BinOp::Mul, ms, x), b.bin(BinOp::Mul, c, y));
    let v = b.bin(BinOp::Add, sx, cy);
    let v = if nut.left { b.un(UnOp::Neg, v) } else { v };
    let h = b.c(-5.5);
    let h = match lift {
        Some(v) => {
            let l = b.v(v);
            b.bin(BinOp::Add, h, l)
        }
        None => h,
    };
    let w = b.bin(BinOp::Sub, z, h);
    let tapped = b.thread_iso([u, v, w], 3.08, nut.pitch);
    let body = cylinder(&mut b, 4.0, -6.0, -2.0);
    let tapped = b.un(UnOp::Neg, tapped);
    let nut = b.bin(BinOp::Max, body, tapped);
    let f = b.bin(BinOp::Max, screw, nut);
    let tape = Arc::new(b.tape(f));
    empty_now(&tape, None, ctx, [-6.0, -6.0, -9.0], [6.0, 6.0, 4.0], None, 0.0125, tol)
}

/// The camera slider's accuracy.
const SLIDER: f64 = 0.02575;

/// A screw seated in phase in its nut, 0.035 clear on its flanks, is
/// proven clear outright on one cross-section of the engagement, one
/// pitch tall (both threads are kept by every screw motion of the
/// pitch, so the overlap's value anywhere is its value there): 49
/// boxes, the cross-section's counted, where cells of four times the
/// accuracy, each settled by its samples, took 18,009 and held it only
/// to the accuracy.
#[test]
fn a_seated_thread_is_proven_on_its_cross_section() {
    let (p, tested) = seated(Nut { turn: 0.0, pitch: 0.5, left: false }, SLIDER);
    assert_eq!(p.verdict, Verdict::Proven);
    assert!(p.exact, "held to the accuracy only");
    assert!(tested < 200, "{tested} boxes");
}

/// The cross-section never shows clear threads that meet. The nut
/// turned off its phase by a quarter turn or more (0.125 mm along its
/// axis or more, past the flanks' 0.040: 0.037 deep or more, by
/// sampling), of another pitch, or left-hand, overlaps the screw
/// deeper than the accuracy, and a witness is found.
#[test]
fn threads_out_of_phase_or_crossed_are_caught() {
    let nuts = [90.0, 120.0, 180.0, -90.0]
        .map(|turn| Nut { turn, pitch: 0.5, left: false })
        .into_iter()
        .chain([Nut { turn: 0.0, pitch: 0.45, left: false }, Nut { turn: 0.0, pitch: 0.5, left: true }]);
    for nut in nuts {
        let (p, _) = seated(nut, SLIDER);
        assert!(matches!(p.verdict, Verdict::Disproven { .. }), "turned {}, pitch {}, left {}: {:?}", nut.turn, nut.pitch, nut.left, p.verdict);
    }
}

/// Turned 30 or 60 degrees off its phase (0.042 or 0.083 mm along its
/// axis) the nut overlaps the screw on one flank, 0.0005 or 0.019 deep
/// (by sampling), under the accuracy: never proven clear outright, only
/// to the accuracy; and asked exactly, not at all.
#[test]
fn threads_that_barely_meet_are_not_proven_clear() {
    for turn in [30.0, 60.0] {
        let nut = Nut { turn, pitch: 0.5, left: false };
        let (p, _) = seated(nut, SLIDER);
        assert!(!(p.verdict == Verdict::Proven && p.exact), "turned {turn}: {:?}", p.verdict);
    }
    let (p, _) = seated(Nut { turn: 30.0, pitch: 0.5, left: false }, 0.0);
    assert_ne!(p.verdict, Verdict::Proven);
}

/// A cross-section proven over one scope settles boxes of that scope
/// alone: the same tape asked at another value of a param it reads
/// (here a field at least 0.25 less p, p 0 then 0.4) is asked again,
/// and fails, never taken as held.
#[test]
fn a_cross_section_holds_only_in_the_scope_it_was_proven_in() {
    let mut b = crate::test_fields::B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let a = b.bin(BinOp::Atan2, y, x);
    let a = b.un(UnOp::Deg, a);
    let c360 = b.c(360.0);
    let a = b.bin(BinOp::Div, a, c360);
    let q = b.bin(BinOp::Sub, a, z);
    let fl = b.un(UnOp::Floor, q);
    let w = b.bin(BinOp::Sub, q, fl);
    let (lo, hi) = (b.c(0.25), b.c(0.75));
    let (l, r) = (b.bin(BinOp::Sub, w, lo), b.bin(BinOp::Sub, hi, w));
    let f = b.bin(BinOp::Max, l, r);
    let p = b.v(Var::Param(0));
    let f = b.bin(BinOp::Sub, f, p);
    let tape = b.tape(f);
    let region = [I::new(1.0, 2.0), I::new(0.1, 0.2), I::new(0.0, 1.0)];
    let q = Empty { tol: 0.0, region, asked: region, loose: Default::default(), t_band: None, read: [true; 3], finest: 0.01, sections: Default::default() };
    assert!(q.by_section(&tape, &Arc::new(Ctx::new(0.0, vec![0.0])), region));
    let other = Arc::new(Ctx::new(0.0, vec![0.4]));
    assert!(eval_point(&tape, &other, [1.5, 0.15, 0.5158627587152768]) < 0.0);
    assert!(!q.by_section(&tape, &other, region), "held at p = 0.4 by the cross-section proven at p = 0");
}

/// A nut lifted up the screw's axis by a ranged param or over a span of
/// time, through a quarter pitch (out of phase, overlapping 0.037 deep
/// at its end): the search halves the span, and each half's boxes are
/// settled only by cross-sections proven over that half, so it is
/// never held (it is unknown, as a range's overlap may be); over a span
/// of a hundredth of a pitch, still clear, it holds outright.
#[test]
fn a_seated_thread_moved_over_a_range_is_proven_over_each_value() {
    let nut = Nut { turn: 0.0, pitch: 0.5, left: false };
    let param = |span: (f64, f64)| Ctx { p_spans: ink_core::eval::ParamSpans::new([(0, span)]), ..Ctx::new(0.0, vec![0.5 * (span.0 + span.1)]) };
    let time = |span: (f64, f64)| Ctx { t_span: Some(span), ..Ctx::new(span.0, vec![]) };
    for (lift, ctx) in [(Var::Param(0), param((0.0, 0.125))), (Var::T, time((0.0, 0.125)))] {
        let (p, _) = seated_at(nut, SLIDER, Some(lift), &ctx);
        assert_ne!(p.verdict, Verdict::Proven, "{lift:?}");
    }
    for (lift, ctx) in [(Var::Param(0), param((0.0, 0.005))), (Var::T, time((0.0, 0.005)))] {
        let (p, _) = seated_at(nut, SLIDER, Some(lift), &ctx);
        assert!(p.verdict == Verdict::Proven && p.exact, "{lift:?}: {:?} exact {}", p.verdict, p.exact);
    }
}

/// A nut turned by a few degrees off its phase but still clear (2
/// degrees: 0.003 mm along its axis) is proven clear outright too: the
/// turn's float cosine and sine are a pair of forms that turn together,
/// so the cross-section still holds them.
#[test]
fn a_nut_turned_within_its_clearance_is_proven_clear() {
    let (p, tested) = seated(Nut { turn: 2.0, pitch: 0.5, left: false }, SLIDER);
    assert_eq!(p.verdict, Verdict::Proven);
    assert!(p.exact && tested < 200, "{tested} boxes");
}

/// `extent_now`, the instructions it evaluated and the longest tape
/// among them.
fn extent_work(tape: &Tape) -> (Extent, (usize, usize)) {
    WORK.with(|w| w.set((0, 0)));
    let e = extent_now(&Arc::new(tape.clone()), &Ctx::new(0.0, vec![]), ([-500.0; 3], [500.0; 3]));
    (e, WORK.with(|w| w.get()))
}

#[test]
fn a_larger_extent_budget_resumes_both_covers() {
    // A box well away from the window, too small to meet any of
    // the probes for far material. Compare each new budget with a cold
    // search at that same budget, including the first box that settles.
    let mut b = crate::test_fields::B::new();
    let part = b.cuboid([14.0, 9.0, 7.0], [1.0; 3]);
    let tape = Arc::new(b.tape(part));
    let root = [tape.len() - 1];
    let ctx = Ctx::new(0.0, vec![]);
    let window = ([-1.0; 3], [1.0; 3]);
    let mut kept = Found::default();
    let mut work = 0;
    let mut cold = 0;
    for budget in [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192] {
        WORK.with(|w| w.set((0, 0)));
        let answer = kept.search(&tape, &root, &ctx, window, budget);
        work += WORK.with(|w| w.get().0);
        WORK.with(|w| w.set((0, 0)));
        assert_eq!(answer, Found::default().search(&tape, &root, &ctx, window, budget));
        cold += WORK.with(|w| w.get().0);
        if let Some(answer) = answer {
            assert!(matches!(answer, Extent::Bounded(..)));
            assert!(work * 4 < cold * 3, "{work} instructions with resumption, {cold} starting over");
            return;
        }
    }
    panic!("neither cover bounded the box");
}

/// A cover's walk is the same on any number of threads (its cells
/// stepped ahead on the pool: ExtentWalk::walk): at each budget a resumed
/// search answers alike and counts the same instructions on one thread
/// and on four, through budgets that stop it and the one that settles it.
#[test]
fn an_extent_walk_is_the_same_on_any_number_of_threads() {
    let mut b = crate::test_fields::B::new();
    let part = b.cuboid([14.0, 9.0, 7.0], [1.0; 3]);
    let tape = Arc::new(b.tape(part));
    let root = [tape.len() - 1];
    let ctx = Ctx::new(0.0, vec![]);
    let window = ([-1.0; 3], [1.0; 3]);
    let walk = |threads: usize| {
        rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap().install(|| {
            let mut kept = Found::default();
            [8, 16, 32, 64, 128, 256, 512, 1024, 2048].map(|budget| {
                WORK.with(|w| w.set((0, 0)));
                (kept.search(&tape, &root, &ctx, window, budget), WORK.with(|w| w.get()))
            })
        })
    };
    let one = walk(1);
    assert!(one.iter().any(|x| x.0.is_none()) && one.iter().any(|x| x.0.is_some()), "{one:?}");
    assert_eq!(one, walk(4));
}

/// A field with material where no model reaches (a half-space, a
/// hole's outside) is unbounded at once: a union with a half-space is
/// unbounded after the tests of its members' far points, not searched
/// at every budget up to the largest (32 thousand evaluations of x).
#[test]
fn material_where_no_model_reaches_is_unbounded_at_once() {
    let mut b = crate::test_fields::B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let cube = b.cuboid([0.0; 3], [1.0; 3]);
    let u = b.bin(BinOp::Min, cube, x);
    let (e, (work, _)) = extent_work(&b.tape(u));
    assert_eq!(e, Extent::Unbounded);
    assert!(work <= 26 * b.tape(u).len(), "{work} instructions evaluated");
    // Over a hole's outside a plate's box is the plate's.
    let plate = b.cuboid([100.0, 0.0, 0.0], [4.0, 4.0, 1.0]);
    let (yy, zz) = (b.un(UnOp::Square, y), b.un(UnOp::Square, z));
    let (r2, d) = (b.c(1.0), b.bin(BinOp::Add, yy, zz));
    let hole = b.bin(BinOp::Sub, r2, d);
    let cut = b.bin(BinOp::Max, plate, hole);
    let (e, _) = extent_work(&b.tape(cut));
    assert!(matches!(e, Extent::Bounded(lo, hi) if lo[0] <= 96.0 && hi[0] >= 104.0), "{e:?}");
}

/// An intersection of a union with itself under another name is the
/// union: unbounded (a union with a half-space), not split.
#[test]
fn an_intersection_of_one_member_is_that_member() {
    let mut g = Graph::new();
    let (x, y) = (g.var(Var::X), g.var(Var::Y));
    let u = g.min(x, y);
    let (a, b) = (g.tag("left", u), g.tag("right", u));
    let f = g.max(a, b);
    assert_eq!(extent_work(&Tape::new(&g, f)).0, Extent::Unbounded);
}

/// An assembly's box asked of a within claim (its parts' union, cut by
/// the outside of the box asked) is its parts' boxes, each part's
/// searched on its own tape, the box's outside with it, never the
/// whole: one part far out, which only a larger budget bounds, leaves
/// the whole undecided at the budgets below, where the whole was
/// searched (its tape 41 parts long) at each.
#[test]
fn an_assemblys_extent_searches_its_parts_not_the_whole() {
    let mut b = crate::test_fields::B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let mut parts = Vec::new();
    for i in 0..40 {
        let c = [-195.0 + 10.0 * i as f64, 7.0, 3.0];
        let plate = b.cuboid(c, [4.0, 4.0, 1.0]);
        let (cx, cy) = (b.c(c[0]), b.c(c[1]));
        let (dx, dy) = (b.bin(BinOp::Sub, x, cx), b.bin(BinOp::Sub, y, cy));
        let (dx, dy) = (b.un(UnOp::Square, dx), b.un(UnOp::Square, dy));
        let r = b.bin(BinOp::Add, dx, dy);
        let r = b.un(UnOp::Sqrt, r);
        let one = b.c(1.0);
        let hole = b.bin(BinOp::Sub, one, r);
        parts.push(b.bin(BinOp::Max, plate, hole));
    }
    // A ball 2e7 out along x.
    let far = b.c(2e7);
    let dx = b.bin(BinOp::Sub, x, far);
    let (dx, dy, dz) = (b.un(UnOp::Square, dx), b.un(UnOp::Square, y), b.un(UnOp::Square, z));
    let r = b.bin(BinOp::Add, dx, dy);
    let r = b.bin(BinOp::Add, r, dz);
    let r = b.un(UnOp::Sqrt, r);
    let one = b.c(1.0);
    parts.push(b.bin(BinOp::Sub, r, one));
    let parts_each = parts.clone();
    while parts.len() > 1 {
        parts = parts.chunks(2).map(|p| if p.len() == 2 { b.bin(BinOp::Min, p[0], p[1]) } else { p[0] }).collect();
    }
    let asked = b.cuboid([0.0; 3], [300.0; 3]);
    let outside = b.un(UnOp::Neg, asked);
    let f = b.bin(BinOp::Max, parts[0], outside);
    let tape = b.tape(f);
    let (e, (_, longest)) = extent_work(&tape);
    assert!(matches!(e, Extent::Bounded(lo, hi) if lo[0] <= -199.0 && hi[0] >= 2e7 + 1.0), "{e:?}");
    let each: Vec<_> = parts_each.iter().map(|&p| b.bin(BinOp::Max, p, outside)).collect();
    let part = each.iter().map(|&f| b.tape(f).len()).max().unwrap();
    assert!(longest <= part, "a tape {longest} long searched, a part with the box's outside {part}, the whole {}", tape.len());
}

/// "Inside an assembly and not inside a box" (a within claim) reads
/// the assembly and the box as simplified over each box searched, as
/// the question's own tape is: never the whole assembly at every box.
/// 40 plates, each a hole cut, in a box holding them: a is read less
/// than once whole a box searched, on average (ten times, read whole
/// at each).
#[test]
fn a_subsets_search_reads_its_parts_near_each_box() {
    let mut b = crate::test_fields::B::new();
    let (x, y) = (b.v(Var::X), b.v(Var::Y));
    let mut parts = Vec::new();
    for i in 0..40 {
        let c = [-195.0 + 10.0 * i as f64, 7.0, 3.0];
        let plate = b.cuboid(c, [4.0, 4.0, 1.0]);
        let (cx, cy) = (b.c(c[0]), b.c(c[1]));
        let (dx, dy) = (b.bin(BinOp::Sub, x, cx), b.bin(BinOp::Sub, y, cy));
        let (dx, dy) = (b.un(UnOp::Square, dx), b.un(UnOp::Square, dy));
        let r = b.bin(BinOp::Add, dx, dy);
        let r = b.un(UnOp::Sqrt, r);
        let one = b.c(1.0);
        let hole = b.bin(BinOp::Sub, one, r);
        parts.push(b.bin(BinOp::Max, plate, hole));
    }
    let plate = b.tape(parts[0]).len();
    while parts.len() > 1 {
        parts = parts.chunks(2).map(|p| if p.len() == 2 { b.bin(BinOp::Min, p[0], p[1]) } else { p[0] }).collect();
    }
    let asked = b.cuboid([0.0; 3], [200.0, 20.0, 20.0]);
    let outside = b.un(UnOp::Neg, asked);
    let f = b.bin(BinOp::Max, parts[0], outside);
    let (a, z) = (Arc::new(b.tape(parts[0])), Arc::new(b.tape(asked)));
    let except = Except { a: a.clone(), b: z };
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let (p, tested, (read, _)) = pool.install(|| {
        WORK.with(|w| w.set((0, 0)));
        let (p, tested) = empty_now(&Arc::new(b.tape(f)), Some(&except), &Ctx::new(0.0, vec![]), [-250.0; 3], [250.0; 3], None, 0.025, 0.05);
        (p, tested, WORK.with(|w| w.get()))
    });
    assert_eq!(p.verdict, Verdict::Proven);
    assert!(read < tested * a.len(), "{read} instructions of a read over {tested} boxes; a plate's tape {plate}, the whole {}", a.len());
}

/// "Inside a and not inside b" where h reads at least 0 over a box is
/// settled there by one bound of a and one pass over b's slopes, not by
/// bounds over the box's parts: a pin 0.03 inside its bore's wall and
/// flush with its ends, at accuracy 0.05, reads a and b, each as long,
/// under four times a's length a box searched (where h's bound is 0 at
/// its ends, reading them over each box's parts had read a about 50
/// times a box, and left 12 times the boxes open).
#[test]
fn a_subset_flush_with_b_reads_a_and_b_seldom() {
    let mut b = crate::test_fields::B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let (xx, yy) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y));
    let rr = b.bin(BinOp::Add, xx, yy);
    let r = b.un(UnOp::Sqrt, rr);
    let ends = b.slab(z, 0.0, 1.0);
    let mut pin = |radius: f64| {
        let k = b.c(radius);
        let wall = b.bin(BinOp::Sub, r, k);
        b.bin(BinOp::Max, wall, ends)
    };
    let (inner, outer) = (pin(1.0), pin(1.03));
    let out = b.un(UnOp::Neg, outer);
    let f = b.bin(BinOp::Max, inner, out);
    let (a, z) = (Arc::new(b.tape(inner)), Arc::new(b.tape(outer)));
    let except = Except { a: a.clone(), b: z };
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let (p, tested, (read, _)) = pool.install(|| {
        WORK.with(|w| w.set((0, 0)));
        let (p, tested) = empty_now(&Arc::new(b.tape(f)), Some(&except), &Ctx::new(0.0, vec![]), [-3.0; 3], [3.0; 3], None, 0.025, 0.05);
        (p, tested, WORK.with(|w| w.get()))
    });
    assert_eq!(p.verdict, Verdict::Proven);
    assert!(read < 4 * tested * a.len(), "{read} instructions of a and b read over {tested} boxes; a's tape {}", a.len());
}

/// A subset settled where a's points reading 0 on b lie on a 0 of no
/// volume holds to the accuracy, not outright: those points are not in
/// b. b reads 0 on a sphere and is below 0 off it, the sphere inside a
/// cube: b holds the cube but for the sphere, a sheet of no depth.
#[test]
fn a_subset_across_a_sheet_of_0_in_b_holds_only_to_the_accuracy() {
    let mut b = crate::test_fields::B::new();
    let cube = b.cuboid([0.0; 3], [1.0; 3]);
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let (xx, yy, zz) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y), b.un(UnOp::Square, z));
    let xy = b.bin(BinOp::Add, xx, yy);
    let rr = b.bin(BinOp::Add, xy, zz);
    let quarter = b.c(0.25);
    let d = b.bin(BinOp::Sub, rr, quarter);
    let dd = b.un(UnOp::Square, d);
    let sheet = b.un(UnOp::Neg, dd);
    let out = b.un(UnOp::Neg, sheet);
    let f = b.bin(BinOp::Max, cube, out);
    let except = Except { a: Arc::new(b.tape(cube)), b: Arc::new(b.tape(sheet)) };
    let p = prove_empty_except(&Arc::new(b.tape(f)), &except, &Ctx::new(0.0, vec![]), [-3.0; 3], [3.0; 3], None, 0.025, 0.05);
    assert_eq!(p.verdict, Verdict::Proven);
    assert!(!p.exact && !p.loose.is_empty(), "held outright, with a's points on b's sheet of 0 outside b");
}

/// A subset's points on b's 0 at the edge of the region asked, whose
/// outside of b beside them lies past it, are not shown absent: a cube
/// about the origin and b = max(-x, -y), asked over x, y at least 0,
/// reads 0 on b at (0, 0.5, 0) inside a. At accuracy 0.05 it holds to
/// the accuracy, not outright; exactly, it is not proven.
#[test]
fn a_subsets_0s_on_b_at_the_edge_asked_are_not_held_outright() {
    let mut b = crate::test_fields::B::new();
    let cube = b.cuboid([0.0; 3], [1.0; 3]);
    let (x, y) = (b.v(Var::X), b.v(Var::Y));
    let (nx, ny) = (b.un(UnOp::Neg, x), b.un(UnOp::Neg, y));
    let quadrant = b.bin(BinOp::Max, nx, ny);
    let out = b.un(UnOp::Neg, quadrant);
    let f = b.bin(BinOp::Max, cube, out);
    let except = Except { a: Arc::new(b.tape(cube)), b: Arc::new(b.tape(quadrant)) };
    let tape = Arc::new(b.tape(f));
    let ctx = Ctx::new(0.0, vec![]);
    let held = prove_empty_except(&tape, &except, &ctx, [0.0, 0.0, -3.0], [3.0; 3], None, 0.025, 0.05);
    assert_eq!(held.verdict, Verdict::Proven);
    assert!(!held.exact && !held.loose.is_empty(), "held outright, with a's points on b's 0 at the edge asked outside b");
    let exactly = prove_empty_except(&tape, &except, &ctx, [0.0, 0.0, -3.0], [3.0; 3], None, 0.025, 0.0);
    assert_ne!(exactly.verdict, Verdict::Proven, "proven exactly, with a's points on b's 0 at the edge asked outside b");
}

/// A subset's points on b's 0 at the edge of the region asked, where
/// b falls throughout, are not shown absent by any settling of a box: a
/// cube about the origin and b = -x, asked over x at least 0, reads 0
/// on b at (0, 0.5, 0) inside a, b's outside past the region; and b
/// written with switches that cancel (-x + (y < 0 ? 1 : 2) + (y < 0 ?
/// -1 : -2)), whose boxes settle by the switches' cases. At accuracy
/// 0.05 each holds to the accuracy, not outright; exactly, neither is
/// proven.
#[test]
fn a_subsets_0s_on_a_falling_b_at_the_edge_asked_are_not_held_outright() {
    for switched in [false, true] {
        let mut g = Graph::new();
        let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
        let one = g.constant(1.0);
        let sides: Vec<_> = [x, y, z]
            .into_iter()
            .map(|v| {
                let av = g.un(UnOp::Abs, v);
                g.sub(av, one)
            })
            .collect();
        let xy = g.max(sides[0], sides[1]);
        let cube = g.max(xy, sides[2]);
        let nx = g.neg(x);
        let b = if switched {
            let (two, mone, mtwo) = (g.constant(2.0), g.constant(-1.0), g.constant(-2.0));
            let (p, n) = (g.select(y, one, two), g.select(y, mone, mtwo));
            let pn = g.add(p, n);
            g.add(nx, pn)
        } else {
            nx
        };
        let out = g.neg(b);
        let f = g.max(cube, out);
        let except = Except { a: Arc::new(Tape::new(&g, cube)), b: Arc::new(Tape::new(&g, b)) };
        let tape = Arc::new(Tape::new(&g, f));
        let ctx = Ctx::new(0.0, vec![]);
        let held = prove_empty_except(&tape, &except, &ctx, [0.0, -3.0, -3.0], [3.0; 3], None, 0.025, 0.05);
        assert_eq!(held.verdict, Verdict::Proven, "switched {switched}");
        assert!(!held.exact && !held.loose.is_empty(), "switched {switched}: held outright, with a's points on b's 0 at the edge asked outside b");
        let exactly = prove_empty_except(&tape, &except, &ctx, [0.0, -3.0, -3.0], [3.0; 3], None, 0.025, 0.0);
        assert_ne!(exactly.verdict, Verdict::Proven, "switched {switched}: proven exactly, with a's points on b's 0 at the edge asked outside b");
    }
}

/// A field rewritten (x to p, y to
/// a name of p) and then untagged is p - p = 0 and 1e16 + 0.5 above it
/// in real numbers: +0.25 everywhere, the probe outside. Untagging is
/// `Graph::rebuild`, like every rewrite: float folds of its own would
/// fold the marked arithmetic to -0.25 and prove it inside (and its
/// negation empty), exactly.
#[test]
fn an_untagged_rewrite_keeps_spatial_arithmetic_real() {
    let ctx = Ctx::new(0.0, vec![0.5]);
    for negated in [false, true] {
        let mut g = Graph::new();
        let (x, y, p) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Param(0)));
        let (half, big, quarter) = (g.constant(0.5), g.constant(1e16), g.constant(0.25));
        let d = g.sub(x, y);
        let a = g.add(d, half);
        let b = g.add(big, a);
        let c = g.sub(b, big);
        let f = g.sub(c, quarter);
        let f = if negated { g.neg(f) } else { f };
        let named = g.tag("same-input", p);
        let rewritten = g.subst(f, &[(Var::X, p), (Var::Y, named)]);
        let untagged = g.untagged(rewritten);
        let tape = Arc::new(Tape::new(&g, untagged));
        let v = eval_interval(&tape, &ctx, [I::new(-5.0, 5.0); 3], None);
        let real = if negated { -0.25 } else { 0.25 };
        assert!(v.lo <= real && real <= v.hi, "untagged, {} is {real}, bounded [{}, {}]", if negated { "its negation" } else { "the field" }, v.lo, v.hi);
        let proof = prove_empty(&tape, &ctx, [-5.0; 3], [5.0; 3], Some(0.1), 0.05, 0.05);
        let wrong = if negated { matches!(proof.verdict, Verdict::Proven) } else { matches!(proof.verdict, Verdict::Disproven { .. }) };
        assert!(!wrong, "untagged, {} is decided wrong: {:?}", if negated { "its negation" } else { "the field" }, proof.verdict);
    }
}

/// The same question with its param
/// numbered 0 or 71 (a part imported after others) has one cache key,
/// marked real or not.
#[test]
fn a_real_question_keeps_its_key_whatever_its_params_number() {
    // A question has a key only where answers are kept (INK_NO_CACHE: nothing is, and nothing is asked here).
    if std::env::var_os("INK_NO_CACHE").is_some() {
        return;
    }
    let opts = crate::ranged::ProofOptions { res: Some(0.1), finest: 0.05, depth: 12, tol: 0.05, boxes: crate::query::BOXES };
    let key = |k: u32, real: bool| {
        let mut g = Graph::new();
        let (x, p) = (g.var(Var::X), g.var(Var::Param(k)));
        let (big, quarter) = (g.constant(1e16), g.constant(0.25));
        let on = if real { x } else { p };
        let s = g.add(big, on);
        let s = g.sub(s, big);
        let f = g.sub(s, quarter);
        let f = if real { g.subst(f, &[(Var::X, p)]) } else { f };
        assert_eq!(g.reachable(&[f]).iter().any(|&n| matches!(g.nodes[n as usize], ink_core::graph::Node::Real(_))), real);
        let mut params = vec![0.0; k as usize + 1];
        params[k as usize] = 0.5;
        crate::ranged::prove_key(&Tape::new(&g, f), None, &Ctx::new(0.0, params), [-5.0; 3], [5.0; 3], &opts)
    };
    for real in [true, false] {
        let (a, b) = (key(0, real), key(71, real));
        assert!(a.is_some() && a == b, "real {real}: keys {a:?}, {b:?}");
    }
    assert!(key(0, true) != key(0, false), "a real question and its float are different questions");
}

#[test]
fn a_density_varying_inside_full_cells_is_measured() {
    // A 10 mm cube of 0.001 g/mm^3 with a 0.2 mm layer of 0.011 across
    // x = 0: 1.2 g. Its full cells' 8 points each miss the layer at
    // every cell size (1.084 g, settled); their bounds see it.
    use crate::test_fields::*;
    use ink_core::graph::{BinOp, UnOp};
    let ctx = at(0.0);
    let mut b = B::new();
    let cube = b.cuboid([0.0; 3], [5.0; 3]);
    let tape = b.clip(cube);
    let x = b.v(Var::X);
    let ax = b.un(UnOp::Abs, x);
    let edge = b.c(0.1);
    let side = b.bin(BinOp::Sub, ax, edge);
    let (dense, none) = (b.c(0.01), b.c(0.0));
    let layer = b.0.select(side, dense, none);
    let base = b.c(0.001);
    let rho = b.bin(BinOp::Add, base, layer);
    let rho = b.tape(rho);
    let (lo, hi, _) = bounds(&tape, &ctx, LO, HI).unwrap();
    let m = measured(&tape, &[Some(&rho), None], true, &ctx, lo, hi, Until::Accuracy(ACCURACY));
    assert!(m.settled && m.lo[0] <= 1.2 && 1.2 <= m.hi[0] && m.error(0) <= ACCURACY * m.sums[0], "{m:?}");
    // At a cell size asked for, the layer is resolved to it.
    let m = measured(&tape, &[Some(&rho), None], true, &ctx, lo, hi, Until::Res(0.025));
    assert!(m.lo[0] <= 1.2 && 1.2 <= m.hi[0] && (m.sums[0] / 1.2 - 1.0).abs() <= ACCURACY, "{m:?}");
}

#[test]
fn sums_with_an_integrand_in_another_slot_are_another_question() {
    // 65 integrands, one a tape and the rest None (1), the tape at slot
    // 0 or at slot 64: the presence of each slot once hashed into one
    // word, bit k % 64, so the two shared a key and one's cached sums
    // answered the other.
    let mut g = Graph::new();
    let (inside, two) = (g.constant(-1.0), g.constant(2.0));
    let shape = Arc::new(Tape::new(&g, inside));
    let two = Tape::new(&g, two);
    let ctx = Ctx::new(0.0, vec![]);
    let sums = |at: usize| {
        let mut gs = vec![None; 65];
        gs[at] = Some(&two);
        measured(&shape, &gs, false, &ctx, [-1.0; 3], [1.0; 3], Until::Res(0.5)).sums
    };
    let (a, b) = (sums(0), sums(64));
    let near = |x: f64, y: f64| (x - y).abs() <= 1e-9 * y;
    assert!(near(a[0], 16.0) && near(a[64], 8.0), "{a:?}");
    assert!(near(b[0], 8.0) && near(b[64], 16.0), "{b:?}");
}

#[test]
fn bounds_read_an_outside_cell_only_once() {
    let mut g = Graph::new();
    let x = g.var(Var::X);
    let tape = Arc::new(Tape::new(&g, x));
    let ctx = Ctx::new(0.0, vec![]);
    let (grid, root) = try_root_for([-1.0; 3], [1.0; 3], 1.0, 1).unwrap();
    let (found, reads, _) = bounds_cells(&material_parts(&tape), &ctx, &grid, root);
    assert_eq!(found, ([-1.0; 3], [0.0, 1.0, 1.0]));
    // The root and seven children. Several slabs equal the whole
    // child; neither a second slab nor the full-cell test may reread it.
    assert_eq!(reads, 8);
}

#[test]
fn bounds_visit_extrema_without_walking_the_whole_surface() {
    use ink_core::graph::BinOp;
    let mut b = crate::test_fields::B::new();
    let cube = b.cuboid([0.0; 3], [10.0; 3]);
    // A disconnected component must still enlarge the box, even
    // after the first component has supplied all six extrema.
    let remote = b.cuboid([30.0, 2.0, -1.0], [1.0; 3]);
    let both = b.bin(BinOp::Min, cube, remote);
    let ctx = Ctx::new(0.0, vec![]);
    for f in [cube, both] {
        let tape = Arc::new(b.tape(f));
        let (grid, root) = try_root_for([-40.0; 3], [40.0; 3], 0.25, 1).unwrap();
        let (found, work, _) = bounds_cells(&material_parts(&tape), &ctx, &grid, root);
        // The descent's leaves that may hold a point inside (it keeps
        // those that only touch 0, for surface on a face).
        let held = |t: &Arc<Tape>, c: Cell| eval_interval(t, &ctx, grid.cell_box(c), None).may_be_inside().then(|| grid.cell_bounds(c));
        let leaves: Vec<(V3, V3)> = crate::octree::descend(&tape, &ctx, &grid, root, 1, &held, &|_, c| Some(grid.cell_bounds(c)));
        assert_eq!(found, leaves.into_iter().fold(ink_core::v3::EMPTY_BOX, ink_core::v3::box_union));
        assert!(work < 10_000, "bounding a flat face must not visit all its leaves: {work} cells");
    }
}

/// An assembly's box is bounded part by part, each part (clipped to
/// the window, as a claim's field is) on its own tape into the one
/// box: no cell reads the whole assembly, and the box is the full
/// descent's.
#[test]
fn an_assemblys_bounds_read_each_part_on_its_own() {
    let mut b = crate::test_fields::B::new();
    let (x, y) = (b.v(Var::X), b.v(Var::Y));
    let mut parts = Vec::new();
    for i in 0..40 {
        // A plate, four holes through it.
        let c = [-195.0 + 10.0 * i as f64, 7.0 * (i % 3) as f64, 3.0];
        let mut part = b.cuboid(c, [4.0, 4.0, 1.0]);
        for h in 0..4 {
            let (cx, cy) = (b.c(c[0] + [-2.0, 2.0][h % 2]), b.c(c[1] + [-2.0, 2.0][h / 2]));
            let (dx, dy) = (b.bin(BinOp::Sub, x, cx), b.bin(BinOp::Sub, y, cy));
            let (dx, dy) = (b.un(UnOp::Square, dx), b.un(UnOp::Square, dy));
            let r = b.bin(BinOp::Add, dx, dy);
            let r = b.un(UnOp::Sqrt, r);
            let one = b.c(1.0);
            let hole = b.bin(BinOp::Sub, one, r);
            part = b.bin(BinOp::Max, part, hole);
        }
        parts.push(part);
    }
    while parts.len() > 1 {
        parts = parts.chunks(2).map(|p| if p.len() == 2 { b.bin(BinOp::Min, p[0], p[1]) } else { p[0] }).collect();
    }
    let tape = b.clip(parts[0]);
    let piece = material_parts(&tape).iter().map(|t| t.len()).max().unwrap();
    let ctx = Ctx::new(0.0, vec![]);
    let (grid, root) = try_root_for([-250.0; 3], [250.0; 3], 500.0 / 64.0, 1).unwrap();
    let (found, _, longest) = bounds_cells(&material_parts(&tape), &ctx, &grid, root);
    let held = |t: &Arc<Tape>, c: Cell| eval_interval(t, &ctx, grid.cell_box(c), None).may_be_inside().then(|| grid.cell_bounds(c));
    let leaves: Vec<(V3, V3)> = crate::octree::descend(&tape, &ctx, &grid, root, 1, &held, &|_, c| Some(grid.cell_bounds(c)));
    assert_eq!(found, leaves.into_iter().fold(ink_core::v3::EMPTY_BOX, ink_core::v3::box_union));
    assert!(longest <= piece && 4 * piece < tape.len(), "a tape {longest} long read; a part clipped {piece}, the whole {}", tape.len());
}

#[test]
fn a_witness_over_a_span_is_at_a_moment_of_its_own() {
    // Two cubes of side 2, the second at x = 10 t: over t in [0.1, 0.4]
    // they overlap until 0.2 and are apart after, at the span's middle
    // (0.25) too. The search over the span halves its moments and finds
    // the overlap at one before 0.2, which the witness carries (inside
    // both there), through a proof's cached answer too. With time not
    // halved it would be undecided over the span; one reported at the
    // span's middle would be outside them.
    let mut b = crate::test_fields::B::new();
    let a = b.cuboid([0.0; 3], [1.0; 3]);
    let (x, y, z, t) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z), b.v(Var::T));
    let ten = b.c(10.0);
    let moved = b.bin(ink_core::graph::BinOp::Mul, ten, t);
    let x = b.bin(ink_core::graph::BinOp::Sub, x, moved);
    let (sx, sy, sz) = (b.slab(x, 0.0, 1.0), b.slab(y, 0.0, 1.0), b.slab(z, 0.0, 1.0));
    let m = b.max3(sx, sy, sz);
    let both = b.bin(ink_core::graph::BinOp::Max, a, m);
    let tape = Arc::new(b.tape(both));
    let (lo, hi) = ([-2.0; 3], [6.0, 2.0, 2.0]);
    let ctx = Ctx { t_span: Some((0.1, 0.4)), ..Ctx::new(0.25, vec![]) };
    let inside = |t: f64, w: V3| (0.1..0.2).contains(&t) && eval_point(&tape, &Ctx::new(t, vec![]), w) < 0.0;
    let s = prove_empty(&tape, &ctx, lo, hi, None, 0.05, 0.0);
    assert!(matches!(s.verdict, Verdict::Disproven { witness, t, .. } if inside(t, witness)), "{s:?}");
    let opts = crate::ranged::ProofOptions { res: None, finest: 0.05, depth: 8, tol: 0.0, boxes: crate::query::BOXES };
    let p = crate::ranged::prove(&tape, &ctx, lo, hi, &opts);
    assert!(matches!(p.verdict, crate::query::Verdict::Disproven { witness, t, .. } if inside(t, witness)), "{p:?}");
}

/// A witness at a corner of the params is found in a few thousand
/// boxes: a block at the end of a stack of three lengths a, b, c (each
/// 10 +-0.1) against a wall at g (30.2 +-0.1), which overlap by up to
/// 0.2 only where all four reach their worst ends (a cube of side 0.1,
/// twice the tolerance, fits only within 0.1 of that corner). Each box
/// tries its samples with every param pushed to the end of its span
/// that lowers the field there (scope::pushed), so the corner is found
/// once a box's samples fall in the overlap, whatever order the
/// params' halves are searched in. The passes this replaced halved the
/// lower halves first and searched each to its finest splits before a
/// witness after it could count (benchmark 35: 22 s, millions of
/// cells). Nominal (no ranges), it holds.
#[test]
fn a_witness_at_a_corner_of_the_params_is_found_in_few_boxes() {
    use crate::test_fields::B;
    use ink_core::graph::BinOp;
    let mut b = B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let (pa, pb, pc, pg) = (b.v(Var::Param(0)), b.v(Var::Param(1)), b.v(Var::Param(2)), b.v(Var::Param(3)));
    let ab = b.bin(BinOp::Add, pa, pb);
    let end = b.bin(BinOp::Add, ab, pc);
    let past = b.bin(BinOp::Sub, x, end);
    let before = b.bin(BinOp::Sub, pg, x);
    let (sy, sz) = (b.slab(y, 0.0, 5.0), b.slab(z, 0.0, 5.0));
    let f = b.max3(past, before, sy);
    let f = b.bin(BinOp::Max, f, sz);
    let tape = Arc::new(b.tape(f));
    let band = |v: f64| (v - 0.1, v + 0.1);
    let ranged =
        Ctx { p_spans: ink_core::eval::ParamSpans::new([(0, band(10.0)), (1, band(10.0)), (2, band(10.0)), (3, band(30.2))]), ..Ctx::new(0.0, vec![10.0, 10.0, 10.0, 30.2]) };
    let (lo, hi) = ([29.0, -6.0, -6.0], [31.0, 6.0, 6.0]);
    let (s, n) = empty_now(&tape, None, &ranged, lo, hi, None, 0.025, 0.05);
    let Verdict::Disproven { witness, params, .. } = &s.verdict else { panic!("{s:?}") };
    let at = |k: u32| params.iter().find(|p| p.0 == k).unwrap().1;
    assert!(at(0) + at(1) + at(2) - at(3) >= 0.1 && (at(3)..=at(0) + at(1) + at(2)).contains(&witness[0]) && n < 10_000, "{s:?} after {n} boxes");
    let nominal = Ctx::new(0.0, vec![10.0, 10.0, 10.0, 30.2]);
    assert!(matches!(prove_empty(&tape, &nominal, lo, hi, None, 0.025, 0.05).verdict, Verdict::Proven));
}

/// Where each of several views settles a region apart, only the view
/// nearest settling a box, by halving its reach, has it halved: a block
/// with a round hole bored in from its +x face and one through it
/// along z, a point of their air hidden where material lies ahead of it
/// along each of +z, -z and +x, at some distance each searches (a
/// ranged param per view, std's `hidden`). Only +x sees out of the
/// bore, +z and -z out of the through hole. Halving every view's reach
/// multiplied the boxes by each view's pieces.
#[test]
fn several_views_halve_only_the_one_nearest_settling() {
    use crate::test_fields::B;
    let mut b = B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    // max(sqrt(u^2 + v^2) - r, slab of half width h about c along w): a
    // round hole's field, or with r and h shrunk by 0.5, its air clear of
    // the wall.
    let round = |b: &mut B, (u, v): (NodeId, NodeId), (cu, cv): (f64, f64), r: f64, (w, c, h): (NodeId, f64, f64)| {
        let (cu, cv) = (b.c(cu), b.c(cv));
        let (du, dv) = (b.bin(BinOp::Sub, u, cu), b.bin(BinOp::Sub, v, cv));
        let (uu, vv) = (b.bin(BinOp::Mul, du, du), b.bin(BinOp::Mul, dv, dv));
        let d = b.bin(BinOp::Add, uu, vv);
        let d = b.un(UnOp::Sqrt, d);
        let r = b.c(r);
        let d = b.bin(BinOp::Sub, d, r);
        let along = b.slab(w, c, h);
        b.bin(BinOp::Max, d, along)
    };
    // The block less the bore (radius 1.6 about y = 0, z = 8, x 7 to
    // 37) and the through hole (radius 2.2 about x = -15, y = 0).
    let part = |b: &mut B, x, y, z| {
        let stock = {
            let (a, c, d) = (b.slab(x, 0.0, 30.0), b.slab(y, 0.0, 20.0), b.slab(z, 8.0, 8.0));
            b.max3(a, c, d)
        };
        let bore = round(b, (y, z), (0.0, 8.0), 1.6, (x, 22.0, 15.0));
        let hole = round(b, (x, y), (-15.0, 0.0), 2.2, (z, 8.0, 20.0));
        let holes = b.bin(BinOp::Min, bore, hole);
        let out = b.un(UnOp::Neg, holes);
        b.bin(BinOp::Max, stock, out)
    };
    // Their air 0.5 clear of the walls, inside the block.
    let band = {
        let bore = round(&mut b, (y, z), (0.0, 8.0), 1.1, (x, 18.5, 11.0));
        let hole = round(&mut b, (x, y), (-15.0, 0.0), 1.7, (z, 8.0, 8.0));
        b.bin(BinOp::Min, bore, hole)
    };
    // The part read s ahead along each view: 16 up, 16 down, 60 along +x.
    let ahead = |b: &mut B, k: u32, axis: usize, by: f64| {
        let s = b.v(Var::Param(k));
        let d = b.c(by);
        let d = b.bin(BinOp::Mul, s, d);
        let mut p = [x, y, z];
        p[axis] = b.bin(BinOp::Add, p[axis], d);
        part(b, p[0], p[1], p[2])
    };
    let (up, down, side) = (ahead(&mut b, 0, 2, 16.0), ahead(&mut b, 1, 2, -16.0), ahead(&mut b, 2, 0, 60.0));
    let f = b.max3(band, up, down);
    let f = b.bin(BinOp::Max, f, side);
    let tape = Arc::new(b.tape(f));
    let ctx = Ctx { p_spans: ink_core::eval::ParamSpans::new((0..3).map(|k| (k, (0.0, 1.0)))), ..Ctx::new(0.0, vec![0.5; 3]) };
    let (s, n) = empty_now(&tape, None, &ctx, [-30.0, -20.0, 0.0], [30.0, 20.0, 16.0], None, 0.025, 0.05);
    assert!(matches!(s.verdict, Verdict::Proven) && n < 2_000, "{s:?} after {n} boxes");
}

/// Groups of terms that read their params alone, not space, are still
/// halved where no halving of one settles a half: a unit cube where
/// cos(p) + 0.5 and cos(q) + 0.5 are both below 0 (p, q from 0 to 12)
/// holds a point at p = q = 3, and one where each of 10 (1.6 + cos(v) +
/// sin(v)) (at least 1.8 everywhere) is below 0 is empty; one halving
/// of either param settles neither half, and halving space narrows
/// neither.
#[test]
fn groups_reading_params_alone_are_halved() {
    use crate::test_fields::B;
    let field = |g: &dyn Fn(&mut B, NodeId) -> NodeId| {
        let mut b = B::new();
        let cube = b.cuboid([0.0; 3], [0.5; 3]);
        let (p, q) = (b.v(Var::Param(0)), b.v(Var::Param(1)));
        let (gp, gq) = (g(&mut b, p), g(&mut b, q));
        let f = b.max3(cube, gp, gq);
        Arc::new(b.tape(f))
    };
    let ctx = Ctx { p_spans: ink_core::eval::ParamSpans::new((0..2).map(|k| (k, (0.0, 12.0)))), ..Ctx::new(0.0, vec![6.0; 2]) };
    let (lo, hi) = ([-1.0; 3], [1.0; 3]);
    let below = field(&|b: &mut B, v| {
        let c = b.un(UnOp::Cos, v);
        let h = b.c(0.5);
        b.bin(BinOp::Add, c, h)
    });
    let s = prove_empty(&below, &ctx, lo, hi, None, 0.025, 0.05);
    assert!(matches!(s.verdict, Verdict::Disproven { .. }), "{s:?}");
    let above = field(&|b: &mut B, v| {
        let (c, s) = (b.un(UnOp::Cos, v), b.un(UnOp::Sin, v));
        let d = b.bin(BinOp::Add, c, s);
        let h = b.c(1.6);
        let d = b.bin(BinOp::Add, d, h);
        let ten = b.c(10.0);
        b.bin(BinOp::Mul, ten, d)
    });
    let s = prove_empty(&above, &ctx, lo, hi, None, 0.025, 0.05);
    assert!(matches!(s.verdict, Verdict::Proven), "{s:?}");
}

/// A search that would test far more boxes than it may ends in a few
/// thousand: a ball exactly inside an envelope of its size written
/// otherwise ((r^2 - 400) / 40, the same sphere), asked to 0.0001,
/// whose curved faces a pair's first-order bound cancels only to the
/// square of a box's size, so every box along them settles only at
/// about 0.02 across: millions of boxes. The passes this replaced
/// stopped at their pass limit in a second; the one search tested its
/// four million boxes first (25 s), until it ends where the boxes it
/// would test, at the rate each depth has opened them, pass four times
/// its budget.
#[test]
fn a_search_that_will_not_fit_ends_early() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    let mut b = B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let sq = [x, y, z].map(|v| b.un(UnOp::Square, v));
    let rr = b.bin(BinOp::Add, sq[0], sq[1]);
    let rr = b.bin(BinOp::Add, rr, sq[2]);
    let (r, size, area, scale) = (b.un(UnOp::Sqrt, rr), b.c(20.0), b.c(1600.0), b.c(160.0));
    let ball = b.bin(BinOp::Sub, r, size);
    // The same sphere written apart, (2x)^2 + ... - 1600: shared, its
    // squared radius would split the boxes where the faces meet.
    let two = b.c(2.0);
    let sq2 = [x, y, z].map(|v| {
        let v2 = b.bin(BinOp::Mul, two, v);
        b.un(UnOp::Square, v2)
    });
    let rr2 = b.bin(BinOp::Add, sq2[0], sq2[1]);
    let rr2 = b.bin(BinOp::Add, rr2, sq2[2]);
    let envelope = b.bin(BinOp::Sub, rr2, area);
    let envelope = b.bin(BinOp::Div, envelope, scale);
    let outside = b.un(UnOp::Neg, envelope);
    let f = b.bin(BinOp::Max, ball, outside);
    let tape = Arc::new(b.tape(f));
    let (s, n) = empty_now(&tape, None, &Ctx::new(0.0, vec![]), [-21.0; 3], [21.0; 3], None, 0.00005, 0.0001);
    assert!(matches!(s.verdict, Verdict::Unknown { ended: Ended::Limit, .. }) && n < 100_000, "{s:?} after {n} boxes");
}

/// Faces that touch are settled whatever their area: a plate on its
/// envelope's floor, both faces at z = 0, is apart exactly, the boxes
/// split there. Dipping a few ulps into the floor,
/// whose bounds overlap the floor's by
/// a few ulps across the face,
/// settled box by box only by the tolerance's cover of its samples, at
/// about four tolerances across (it ended unknown after 4.2 million
/// boxes), is proven to 0.01 in one: a pair of faces facing each other
/// along an axis leaves from any point of both within the tolerance,
/// whatever the box's size. Moved 0.1 through the floor, it is a
/// witness. Where two
/// such faces meet (the plate as wide as its envelope: its edges on the
/// envelope's), the envelope's walls are a choice, through which no
/// exit is shown, and the edges are refined to the tolerance: boxes as
/// their length, not as the faces' area (127 thousand here).
#[test]
fn faces_that_touch_are_covered_whatever_their_area() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    let dip = 1.9999999999999996;
    for (side, lift, holds, most) in [(49.0, 2.0, true, 10), (49.0, dip, true, 5_000), (49.0, 1.9, false, 5_000), (50.0, dip, true, 200_000)] {
        let mut b = B::new();
        let plate = b.cuboid([0.0, 0.0, lift], [side, side, 2.0]);
        let envelope = b.cuboid([0.0, 0.0, 15.0], [50.0, 50.0, 15.0]);
        let outside = b.un(UnOp::Neg, envelope);
        let f = b.bin(BinOp::Max, plate, outside);
        let tape = Arc::new(b.tape(f));
        let (s, n) = empty_now(&tape, None, &Ctx::new(0.0, vec![]), [-51.0, -51.0, -1.0], [51.0, 51.0, 31.0], None, 0.005, 0.01);
        let ok = match s.verdict {
            Verdict::Proven => holds && s.exact == (lift == 2.0),
            Verdict::Disproven { witness, .. } => !holds && witness[2] < 0.0,
            _ => false,
        };
        assert!(ok && n < most, "plate {side} lifted {lift}: {s:?} after {n} boxes");
    }
}

/// A ring's band 1.1 wide extruded 1 (z 3 to 4) inside its outer disc
/// extruded 1.2 (2.9 to 4.1): the band's outer wall lies on the disc's
/// all round, which the plane settles in one box (the walls facing). In
/// space, a box inside both spans reads no z, and one at the band's top
/// or bottom face (its slab at most 0 there) the band's sign alone,
/// so the search is the plane's, never halved along z (it halved the
/// walls' boxes as an octree's, along z too: 4,153 boxes).
#[test]
fn an_extrusion_is_searched_as_its_profile() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    let boxes = |solid: bool| {
        let mut b = B::new();
        let (x, y) = (b.v(Var::X), b.v(Var::Y));
        let (xx, yy) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y));
        let rr = b.bin(BinOp::Add, xx, yy);
        let (d, r) = (b.un(UnOp::Sqrt, rr), b.c(18.9));
        let disc = b.bin(BinOp::Sub, d, r);
        let w = b.c(1.1);
        let inner = b.bin(BinOp::Add, disc, w);
        let hole = b.un(UnOp::Neg, inner);
        let band = b.bin(BinOp::Max, disc, hole);
        let mut extrude = |g, z0: f64, z1: f64| match solid {
            true => {
                let z = b.v(Var::Z);
                let s = b.slab(z, 0.5 * (z0 + z1), 0.5 * (z1 - z0));
                b.bin(BinOp::Emax, g, s)
            }
            false => g,
        };
        let (a, c) = (extrude(band, 3.0, 4.0), extrude(disc, 2.9, 4.1));
        let out = b.un(UnOp::Neg, c);
        let h = b.bin(BinOp::Max, a, out);
        let e = Except { a: Arc::new(b.tape(a)), b: Arc::new(b.tape(c)) };
        let (s, n) = empty_now(&Arc::new(b.tape(h)), Some(&e), &Ctx::new(0.0, vec![]), [-25.0, -25.0, 0.0], [25.0, 25.0, 5.0], None, 0.025, 0.05);
        assert!(matches!(s.verdict, Verdict::Proven), "{s:?} after {n} boxes");
        n
    };
    let (solid, plane) = (boxes(true), boxes(false));
    assert!(solid <= 4 * plane, "extruded {solid} boxes, the profile {plane}");
}

/// A block 1.5 x 4 x 4 about x = 1e16 + 1, which is no float: its
/// narrowed box is [1e16, 1e16 + 2], whose float middle is a face. How
/// deep a point can lie is the real half-width, 1, not the depth at
/// that middle (0), so no search shows it empty to 0.1 (and none can
/// show it full: no float point is inside).
#[test]
fn depth_is_half_the_width_not_the_float_middles() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    let mut b = B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let (big, one, w) = (b.c(1e16), b.c(1.0), b.c(0.75));
    let d = b.bin(BinOp::Sub, x, big);
    let d = b.bin(BinOp::Sub, d, one);
    let d = b.un(UnOp::Abs, d);
    let fx = b.bin(BinOp::Sub, d, w);
    let (fy, fz) = (b.slab(y, 0.0, 2.0), b.slab(z, 0.0, 2.0));
    let f = b.max3(fx, fy, fz);
    let tape = Arc::new(b.tape(f));
    let ctx = Ctx::new(0.0, vec![]);
    let Extent::Bounded(lo, hi) = extent(&tape, &ctx, ([-500.0; 3], [500.0; 3])) else { panic!("unbounded") };
    let (lo, hi, _) = bounds(&tape, &ctx, lo, hi).unwrap();
    for tol in [0.0, 0.1, 0.5] {
        let s = prove_empty(&tape, &ctx, lo, hi, None, 0.05, tol);
        match s.verdict {
            Verdict::Unknown { attained, .. } => assert!(attained.is_none_or(|a| a >= 0.75), "{tol}: {s:?}"),
            _ => panic!("{tol}: {s:?}"),
        }
    }
}

/// The radii bound the real ones: squares that underflow (1e-200) or
/// overflow (1e154) and a norm whose float rounds below the real one.
#[test]
fn radial_extent_bounds_the_real_radii() {
    let mut g = Graph::new();
    let n = g.constant(-1.0);
    let tape = Tape::new(&g, n);
    let ctx = Ctx::new(0.0, vec![]);
    let at = |lo: V3, hi: V3| radial_extent(&tape, &ctx, lo, hi, false, 0.0, 100).unwrap();
    let (lower, upper) = at([1e-200, 0.0, 0.0], [2e-200, 0.0, 0.0]);
    assert!(lower > 0.0 && lower <= 2e-200 && upper >= 2e-200, "{lower} {upper}");
    let (lower, upper) = at([1e154, 1e154, 0.0], [1e154, 1e154, 0.0]);
    // sqrt(2) 1e154 = 1.41421356237...e154, not (f64::MAX, inf).
    assert!(lower < 1.4142135624e154 && upper > 1.4142135623e154 && upper - lower < 1e140, "{lower} {upper}");
    // sqrt(x^2 + y^2) in floats rounds to 1.5950777814546493, below the
    // real radius (by an exact rational check), and its next float down.
    let p = [1.0506403964780062, 1.2001782726657002, 0.0];
    let (lower, upper) = at(p, p);
    assert!(lower < 1.5950777814546493 && upper > 1.5950777814546493 && upper - lower < 1e-14, "{lower} {upper}");
}

/// Faces whose slopes are so small that weight times slope underflows
/// to 0 bound no distance: the slab 0 < x < 2 scaled by 1e-323 is not
/// covered by its faces (their reach once came out negative).
#[test]
fn faces_of_a_field_scaled_to_subnormals_cover_nothing() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp, Var};
    let mut b = B::new();
    let (k, x) = (b.c(1e-323), b.v(Var::X));
    let kx = b.bin(BinOp::Mul, k, x);
    let k2 = b.c(2e-323);
    let right = b.bin(BinOp::Sub, kx, k2);
    let left = b.un(UnOp::Neg, kx);
    let f = b.bin(BinOp::Max, right, left);
    let tape = Arc::new(b.tape(f));
    let (s, _) = empty_now(&tape, None, &Ctx::new(0.0, vec![]), [-10.0; 3], [10.0; 3], None, 0.01, 0.01);
    assert!(!matches!(s.verdict, Verdict::Proven), "{s:?}");
}

/// A part outside a zone that is its own feature (a screw's thread, an
/// overlap a contact leaves out: max(part, other, -feature)) is settled
/// at any box size: the max holds the feature's value and its negation,
/// whose mean is 0 at every point however the feature curves: one box
/// (bounds of the two taken apart lose it, and the sphere's surface
/// inside the block settled only by the tolerance, at about four
/// tolerances across: 273 thousand boxes). An overlap outside the
/// feature (the head) is still found.
#[test]
fn a_part_outside_its_own_feature_is_settled_at_any_size() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    for (head, holds) in [(9.0, true), (6.5, false)] {
        let mut b = B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let sq = [x, y, z].map(|v| b.un(UnOp::Square, v));
        let rr = b.bin(BinOp::Add, sq[0], sq[1]);
        let rr = b.bin(BinOp::Add, rr, sq[2]);
        let (r, five) = (b.un(UnOp::Sqrt, rr), b.c(5.0));
        let feature = b.bin(BinOp::Sub, r, five);
        let cap = b.cuboid([0.0, 0.0, head], [1.0, 1.0, 1.0]);
        let part = b.bin(BinOp::Min, feature, cap);
        let other = b.cuboid([0.0, 0.0, 4.0], [4.5, 4.5, 2.5]);
        let out = b.un(UnOp::Neg, feature);
        let f = b.bin(BinOp::Max, part, other);
        let f = b.bin(BinOp::Max, f, out);
        let tape = Arc::new(b.tape(f));
        let (s, n) = empty_now(&tape, None, &Ctx::new(0.0, vec![]), [-6.0; 3], [6.0, 6.0, 11.0], None, 0.005, 0.01);
        let ok = match s.verdict {
            Verdict::Proven => holds,
            Verdict::Disproven { witness, .. } => !holds && witness[2] > 5.5,
            _ => false,
        };
        assert!(ok && n < 5_000, "head at {head}: {s:?} after {n} boxes");
    }
}

/// A box waiting keeps the choices its test made, not a tape: the
/// tapes the boxes left hold are their parents', each shared by its
/// halves, so fewer than the boxes (681 boxes hold 304; kept each its
/// own, 623: benchmark 02's mesh-sweep peaked at 386 MB, 135 MB so).
#[test]
fn boxes_waiting_share_their_parents_tapes() {
    use crate::test_fields::B;
    // Sixty-four cubes in rows, each touching the next, under a plate
    // pressed 0.015 into their tops (thinner than twice the tolerance,
    // too thick for the faces' own cover), tilted 1e-4 along x (so no
    // face of it reads one axis alone): the boxes along the contact
    // keep narrowing the union to the cubes near them.
    use ink_core::graph::BinOp;
    let mut b = B::new();
    let cubes: Vec<_> = (0..64).map(|k| b.cuboid([(k % 8) as f64, (k / 8) as f64, 0.0], [0.5, 0.5, 0.5])).collect();
    let union = cubes.into_iter().reduce(|x, y| b.bin(BinOp::Min, x, y)).unwrap();
    let (x, y, z, tilt) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z), b.c(1e-4));
    let rise = b.bin(BinOp::Mul, tilt, x);
    let up = b.bin(BinOp::Sub, z, rise);
    let (sx, sy, sz) = (b.slab(x, 3.5, 4.0), b.slab(y, 3.5, 4.0), b.slab(up, 0.985, 0.5));
    let plate = b.max3(sx, sy, sz);
    let f = b.bin(BinOp::Max, union, plate);
    let tape = Arc::new(b.tape(f));
    let (lo, hi) = ([-0.5, -0.5, -0.5], [7.5, 7.5, 1.5]);
    let region = std::array::from_fn(|k| I::new(lo[k], hi[k]));
    let q = Empty { tol: 0.01, region, asked: region, loose: Default::default(), t_band: None, read: [true; 3], finest: 0.01, sections: Default::default() };
    let root = search::Node { b: lattice(lo, hi, 0.01).unwrap(), ctx: Arc::new(Ctx::new(0.0, vec![])), state: (Kept::of(&tape), Kept::of(&tape), None) };
    let mut s = search::Search::new(&q, vec![root], search::Precision::holds(0.015));
    for _ in 0..20 {
        if s.round(&q, BOXES).is_some() {
            break;
        }
    }
    let left: Vec<_> = s.left().collect();
    let held: std::collections::HashSet<*const Tape> = left.iter().map(|n| n.state.0.held()).collect();
    assert!(left.len() > 500 && 3 * held.len() < 2 * left.len(), "{} boxes hold {} tapes", left.len(), held.len());
}

/// A box's tape is made when it is read, not when it is cloned.
#[test]
fn a_kept_tape_is_made_when_read() {
    use crate::test_fields::B;
    let mut b = B::new();
    let (u, v) = (b.cuboid([0.0; 3], [1.0; 3]), b.cuboid([3.0, 0.0, 0.0], [1.0; 3]));
    let f = b.bin(ink_core::graph::BinOp::Min, u, v);
    let tape = Arc::new(b.tape(f));
    let mut ch = vec![BOTH; tape.len()];
    eval_interval(&tape, &Ctx::new(0.0, vec![]), [I::new(-1.0, 1.0); 3], Some(&mut ch));
    let k = Kept::after(&tape, vec![ch]);
    let c = k.clone();
    assert!(std::ptr::eq(c.held(), Arc::as_ptr(&tape)), "made as cloned");
    assert!(c.tape().len() < tape.len() && std::ptr::eq(k.held(), Arc::as_ptr(&c.tape())), "made once, shared");
}

/// A ball about (cx, 0, 0) of radius r: its distance, or (squared) the
/// square of it less r^2, a field with no slope bound of 1.
fn ball(b: &mut crate::test_fields::B, cx: f64, r: f64, squared: bool) -> ink_core::graph::NodeId {
    use ink_core::graph::{BinOp, UnOp};
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let c = b.c(cx);
    let dx = b.bin(BinOp::Sub, x, c);
    let (xx, yy, zz) = (b.bin(BinOp::Mul, dx, dx), b.bin(BinOp::Mul, y, y), b.bin(BinOp::Mul, z, z));
    let s = b.bin(BinOp::Add, xx, yy);
    let s = b.bin(BinOp::Add, s, zz);
    match squared {
        true => {
            let r2 = b.c(r * r);
            b.bin(BinOp::Sub, s, r2)
        }
        false => {
            let l = b.un(UnOp::Sqrt, s);
            let r = b.c(r);
            b.bin(BinOp::Sub, l, r)
        }
    }
}

/// A field x + 10 y cos p - 10 y cos(p + 0.001) (plus a millionth of
/// each of eight other params, all over [-2, 2]) moves with p only
/// through its slope: a millionth's change in value, but over p's span
/// its slope bound reaches 40, so a box 0.02 across beside it at 1 to 2
/// along x is cleared by 0.5 only once p is halved. With p the ninth
/// param, past the eight a first-order bound takes at once, its blur
/// was not measured (taken as 0, then as what its span adds to the
/// plain bound, which y's two reads swamp) and the search ended unknown.
#[test]
fn a_param_past_the_first_eight_is_halved_like_them() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    for special in [0, 8] {
        let mut b = B::new();
        let (x, y) = (b.v(Var::X), b.v(Var::Y));
        let ten = b.c(10.0);
        let yy = b.bin(BinOp::Mul, y, ten);
        let p = b.v(Var::Param(special));
        let eps = b.c(0.001);
        let shifted = b.bin(BinOp::Add, p, eps);
        let (cp, cq) = (b.un(UnOp::Cos, p), b.un(UnOp::Cos, shifted));
        let (xp, xq) = (b.bin(BinOp::Mul, yy, cp), b.bin(BinOp::Mul, yy, cq));
        let mut f = b.bin(BinOp::Sub, xp, xq);
        for k in (0..9).filter(|&k| k != special) {
            let (pk, tiny) = (b.v(Var::Param(k)), b.c(1e-6));
            let t = b.bin(BinOp::Mul, pk, tiny);
            f = b.bin(BinOp::Add, f, t);
        }
        let f = b.bin(BinOp::Add, f, x);
        let bar = b.cuboid([1.5, 0.0, 0.0], [0.5, 0.01, 0.01]);
        let ctx = Ctx { p_spans: ink_core::eval::ParamSpans::new((0..9).map(|k| (k, (-2.0, 2.0)))), ..Ctx::new(0.0, vec![0.0; 9]) };
        let way = Way { grown: Arc::new(b.tape(f)), covered: Arc::new(b.tape(bar)), lo: [1.0, -0.01, -0.01], hi: [2.0, 0.01, 0.01], swapped: false };
        let (v, ..) = clearance_now(&[way], &ctx, 0.5, 0.025, 20_000);
        assert!(matches!(v, Verdict::Proven), "p{special}: {v:?}");
    }
}

/// A box reaching infinity is tested whole, never halved, so its
/// infinite spread does not stop a param from being halved: B the
/// half-space x > 1 + sin p - sin(p + 0.001) (plus a millionth of each
/// of eight other params over [-0.5, 0.5]), A x < 0, asked over all of
/// space: proven once p is halved, first or ninth. Ranked with them,
/// space's infinite smear would outrank every param, and the box would
/// be stuck at once.
#[test]
fn a_box_over_all_of_space_halves_its_params() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    for special in [0, 8] {
        let mut b = B::new();
        let x = b.v(Var::X);
        let p = b.v(Var::Param(special));
        let eps = b.c(0.001);
        let shifted = b.bin(BinOp::Add, p, eps);
        let (s0, s1) = (b.un(UnOp::Sin, p), b.un(UnOp::Sin, shifted));
        let mut f = b.bin(BinOp::Sub, s0, s1);
        for k in (0..9).filter(|&k| k != special) {
            let (pk, tiny) = (b.v(Var::Param(k)), b.c(1e-6));
            let t = b.bin(BinOp::Mul, pk, tiny);
            f = b.bin(BinOp::Add, f, t);
        }
        let one = b.c(1.0);
        let f = b.bin(BinOp::Add, one, f);
        let covered = b.bin(BinOp::Sub, f, x);
        let ctx = Ctx { p_spans: ink_core::eval::ParamSpans::new((0..9).map(|k| (k, (-0.5, 0.5)))), ..Ctx::new(0.0, vec![0.0; 9]) };
        let way = Way { grown: Arc::new(b.tape(x)), covered: Arc::new(b.tape(covered)), lo: [f64::NEG_INFINITY; 3], hi: [f64::INFINITY; 3], swapped: false };
        let (v, ..) = clearance_now(&[way], &ctx, 0.5, 0.1, 20_000);
        assert!(matches!(v, Verdict::Proven), "p{special}: {v:?}");
    }
}

/// The one search's cost, counted, both ways round. A pin of radius 5
/// in a bore of 5.5 cleared by 0.49 (a margin of 0.01) is proven by the
/// facing surfaces' first-order bound at boxes of about sqrt(margin x
/// radius), not of the margin. A square of slabs turned 45 degrees,
/// extruded, reads 30 beside its corner where the other is 42.4 away:
/// claimed clear by 31, the box nearest a violation is refined first
/// and, stuck at the finest cells, ends the search, which holds to the
/// 30 its field proves. Two boxes of slabs as far apart are proven by
/// the gaps between their sides (`Flat`). Counts far above those
/// asserted mean it refines where it need not.
#[test]
fn the_one_search_decides_in_few_boxes() {
    use crate::test_fields::B;
    use ink_core::graph::BinOp;
    let mut b = B::new();
    let (x, y) = (b.v(Var::X), b.v(Var::Y));
    let (xx, yy) = (b.bin(BinOp::Mul, x, x), b.bin(BinOp::Mul, y, y));
    let r = b.bin(BinOp::Add, xx, yy);
    let r = b.un(ink_core::graph::UnOp::Sqrt, r);
    let (five, bore) = (b.c(5.0), b.c(5.5));
    let pin = b.bin(BinOp::Sub, r, five);
    let hole = b.bin(BinOp::Sub, bore, r);
    let plate = b.cuboid([0.0; 3], [20.0, 20.0, 5.0]);
    let plate = b.bin(BinOp::Max, plate, hole);
    let (plate, pin) = (Arc::new(b.tape(plate)), Arc::new(b.tape(pin)));
    let ways = [(plate.clone(), pin.clone(), false), (pin, plate, true)].map(|(grown, covered, swapped)| Way { grown, covered, lo: [-6.0; 3], hi: [6.0; 3], swapped });
    let (v, _, _, n) = clearance_now(&ways, &Ctx::new(0.0, vec![]), 0.49, 0.0125, BOXES);
    assert!(matches!(v, Verdict::Proven) && n < 100_000, "{v:?} after {n} boxes");
    let mut b = B::new();
    let (a, other) = (b.cuboid([0.0; 3], [50.0; 3]), b.cuboid([130.0, 130.0, 0.0], [50.0; 3]));
    let (a, other) = (Arc::new(b.tape(a)), Arc::new(b.tape(other)));
    let ways = [
        Way { grown: a.clone(), covered: other.clone(), lo: [80.0, 80.0, -50.0], hi: [180.0, 180.0, 50.0], swapped: false },
        Way { grown: other, covered: a, lo: [-50.0; 3], hi: [50.0; 3], swapped: true },
    ];
    let (v, _, _, n) = clearance_now(&ways, &Ctx::new(0.0, vec![]), 31.0, 0.025, BOXES);
    assert!(matches!(v, Verdict::Proven) && n < 10, "{v:?} after {n} boxes");
    // Squares of half-side 50 turned 45 degrees, their corners along x
    // 42.43 apart: each field reads 42.43 / sqrt 2 = 30 at the other's.
    let mut b = B::new();
    let diamond = |b: &mut B, cx: f64| {
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let c = b.c(cx);
        let x = b.bin(BinOp::Sub, x, c);
        let (u, v) = (b.bin(BinOp::Add, x, y), b.bin(BinOp::Sub, x, y));
        let k = b.c(std::f64::consts::FRAC_1_SQRT_2);
        let (u, v) = (b.bin(BinOp::Mul, u, k), b.bin(BinOp::Mul, v, k));
        let (u, v, w) = (b.slab(u, 0.0, 50.0), b.slab(v, 0.0, 50.0), b.slab(z, 0.0, 50.0));
        b.max3(u, v, w)
    };
    let gap = 30.0 * std::f64::consts::SQRT_2;
    let reach = 50.0 * std::f64::consts::SQRT_2;
    let (a, other) = (diamond(&mut b, 0.0), diamond(&mut b, 2.0 * reach + gap));
    let (a, other) = (Arc::new(b.tape(a)), Arc::new(b.tape(other)));
    let ways = [
        Way { grown: a.clone(), covered: other.clone(), lo: [reach + gap, -reach, -50.0], hi: [3.0 * reach + gap, reach, 50.0], swapped: false },
        Way { grown: other, covered: a, lo: [-reach, -reach, -50.0], hi: [reach, reach, 50.0], swapped: true },
    ];
    let (v, _, _, n) = clearance_now(&ways, &Ctx::new(0.0, vec![]), 31.0, 0.025, BOXES);
    let Verdict::Unknown { ended: Ended::Finest, attained: Some(e), .. } = v else { panic!("{v:?}") };
    assert!((1.0..1.01).contains(&e) && n < 25_000, "holds to {e} after {n} boxes");
}

/// A plate with a bore, the pin in it, each a profile extruded (an
/// emax of the profile and a span, as `extrude-z` writes them), the
/// pin's span past the plate's at both ends. The question is the
/// plane's at every height both spans hold, and the gap between the
/// spans' ends where they part: the boxes stay the profiles' whole
/// height, never halved along z, and cost what the profiles' do
/// (halved as cubes: 36,282 boxes, 26 times the profiles' 1,414).
#[test]
fn an_extrusions_clearance_is_searched_as_its_profiles() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    let boxes = |solid: bool| {
        let mut b = B::new();
        let (x, y) = (b.v(Var::X), b.v(Var::Y));
        let (xx, yy) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y));
        let rr = b.bin(BinOp::Add, xx, yy);
        let r = b.un(UnOp::Sqrt, rr);
        let (five, bore) = (b.c(5.0), b.c(5.5));
        let pin = b.bin(BinOp::Sub, r, five);
        let hole = b.bin(BinOp::Sub, bore, r);
        let (sx, sy) = (b.slab(x, 0.0, 20.0), b.slab(y, 0.0, 20.0));
        let square = b.bin(BinOp::Max, sx, sy);
        let plate = b.bin(BinOp::Max, square, hole);
        let mut extrude = |g, z0: f64, z1: f64| match solid {
            true => {
                let z = b.v(Var::Z);
                let s = b.slab(z, 0.5 * (z0 + z1), 0.5 * (z1 - z0));
                b.bin(BinOp::Emax, g, s)
            }
            false => g,
        };
        let (plate, pin) = (extrude(plate, 0.0, 5.0), extrude(pin, -2.0, 8.0));
        let (plate, pin) = (Arc::new(b.tape(plate)), Arc::new(b.tape(pin)));
        let ways =
            [(plate.clone(), pin.clone(), false), (pin, plate, true)].map(|(grown, covered, swapped)| Way { grown, covered, lo: [-6.0, -6.0, -3.0], hi: [6.0, 6.0, 9.0], swapped });
        let (v, _, _, n) = clearance_now(&ways, &Ctx::new(0.0, vec![]), 0.49, 0.0125, BOXES);
        assert!(matches!(v, Verdict::Proven), "{v:?} after {n} boxes");
        n
    };
    let (solid, plane) = (boxes(true), boxes(false));
    assert!(10 * solid <= 11 * plane, "extruded {solid} boxes, the profiles {plane}");
}

/// A pin beside a plate and above it, each a profile cut to its span by
/// a plain max (an intersection, not an extrusion's emax): their
/// nearest points are the plate's top edge and the pin's bottom edge,
/// sqrt 2 apart, where each field reads 1, the distance along one axis.
/// No bound of a field there proves more than 1; the gap between the
/// spans and the profiles' in the plane, together, prove 1.4 (unknown,
/// held to 1). Over sqrt 2 a pair is found.
#[test]
fn extrusions_apart_along_their_spans_are_as_far_as_the_gap_and_their_profiles() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    let mut b = B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let eight = b.c(8.0);
    let dx = b.bin(BinOp::Sub, x, eight);
    let (xx, yy) = (b.un(UnOp::Square, dx), b.un(UnOp::Square, y));
    let rr = b.bin(BinOp::Add, xx, yy);
    let r = b.un(UnOp::Sqrt, rr);
    let two = b.c(2.0);
    let pin = b.bin(BinOp::Sub, r, two);
    let (sx, sy) = (b.slab(x, 0.0, 5.0), b.slab(y, 0.0, 5.0));
    let plate = b.bin(BinOp::Max, sx, sy);
    let (sp, sq) = (b.slab(z, 2.5, 2.5), b.slab(z, 7.5, 1.5));
    let (plate, pin) = (b.bin(BinOp::Max, plate, sp), b.bin(BinOp::Max, pin, sq));
    let (plate, pin) = (Arc::new(b.tape(plate)), Arc::new(b.tape(pin)));
    let ways =
        [(plate.clone(), pin.clone(), false), (pin, plate, true)].map(|(grown, covered, swapped)| Way { grown, covered, lo: [-6.0, -6.0, -1.0], hi: [11.0, 6.0, 10.0], swapped });
    let (v, _, _, n) = clearance_now(&ways, &Ctx::new(0.0, vec![]), 1.4, 0.0125, BOXES);
    assert!(matches!(v, Verdict::Proven), "{v:?} after {n} boxes");
    let (v, pair, _, _) = clearance_now(&ways, &Ctx::new(0.0, vec![]), 1.42, 0.0125, BOXES);
    let Some((p, q)) = pair else { panic!("{v:?}") };
    assert!(distance_up(p, q) < 1.42, "{p:?} {q:?}");
}

/// Clearances of random extrusions (discs and rectangles, each cut to
/// its span by an extrusion's emax or a plain max, one or two to a
/// shape) against their distance worked out by hand: each pair of
/// prisms is apart by their profiles' distance and their spans' gap in
/// quadrature. Proven only where the distance is at least the
/// clearance, a pair found only where it is under.
#[test]
fn clearances_of_extrusions_agree_with_their_distances() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    #[derive(Clone, Copy)]
    struct Prism {
        c: [f64; 2],
        disc: Option<f64>,
        half: [f64; 2],
        z: [f64; 2],
        emax: bool,
    }
    let mut seed = 12345u64;
    let mut r = |a: f64, b: f64| {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        a + (b - a) * ((seed >> 11) as f64 / (1u64 << 53) as f64)
    };
    let gap = |a: [f64; 2], b: [f64; 2]| (b[0] - a[1]).max(a[0] - b[1]).max(0.0);
    let apart = |p: &Prism, q: &Prism| {
        let to_rect = |c: [f64; 2], h: [f64; 2], x: [f64; 2]| ((x[0] - c[0]).abs() - h[0]).max(0.0).hypot(((x[1] - c[1]).abs() - h[1]).max(0.0));
        let plane = match (p.disc, q.disc) {
            (Some(a), Some(b)) => ((p.c[0] - q.c[0]).hypot(p.c[1] - q.c[1]) - a - b).max(0.0),
            (Some(a), None) => (to_rect(q.c, q.half, p.c) - a).max(0.0),
            (None, Some(b)) => (to_rect(p.c, p.half, q.c) - b).max(0.0),
            (None, None) => gap([p.c[0] - p.half[0], p.c[0] + p.half[0]], [q.c[0] - q.half[0], q.c[0] + q.half[0]])
                .hypot(gap([p.c[1] - p.half[1], p.c[1] + p.half[1]], [q.c[1] - q.half[1], q.c[1] + q.half[1]])),
        };
        plane.hypot(gap(p.z, q.z))
    };
    let (mut proven, mut found) = (0, 0);
    for case in 0..60 {
        let shape = |r: &mut dyn FnMut(f64, f64) -> f64| -> Vec<Prism> {
            (0..if r(0.0, 1.0) < 0.3 { 2 } else { 1 })
                .map(|_| {
                    let z0 = r(-4.0, 3.0);
                    let disc = (r(0.0, 1.0) < 0.5).then(|| r(0.5, 3.0));
                    Prism { c: [r(-6.0, 6.0), r(-6.0, 6.0)], disc, half: [r(0.3, 3.0), r(0.3, 3.0)], z: [z0, z0 + r(0.5, 5.0)], emax: r(0.0, 1.0) < 0.7 }
                })
                .collect()
        };
        let (pa, pb) = (shape(&mut r), shape(&mut r));
        let truth = pa.iter().flat_map(|p| pb.iter().map(move |q| apart(p, q))).fold(f64::INFINITY, f64::min);
        let d = if case % 2 == 0 { (truth + r(-0.3, 0.3)).max(0.05) } else { r(0.05, 4.0) };
        let mut b = B::new();
        let mut field = |ps: &[Prism]| {
            let parts: Vec<_> = ps
                .iter()
                .map(|p| {
                    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
                    let profile = match p.disc {
                        Some(radius) => {
                            let (cx, cy) = (b.c(p.c[0]), b.c(p.c[1]));
                            let (dx, dy) = (b.bin(BinOp::Sub, x, cx), b.bin(BinOp::Sub, y, cy));
                            let (dx, dy) = (b.un(UnOp::Square, dx), b.un(UnOp::Square, dy));
                            let rr = b.bin(BinOp::Add, dx, dy);
                            let (d, radius) = (b.un(UnOp::Sqrt, rr), b.c(radius));
                            b.bin(BinOp::Sub, d, radius)
                        }
                        None => {
                            let (sx, sy) = (b.slab(x, p.c[0], p.half[0]), b.slab(y, p.c[1], p.half[1]));
                            b.bin(BinOp::Max, sx, sy)
                        }
                    };
                    let span = b.slab(z, 0.5 * (p.z[0] + p.z[1]), 0.5 * (p.z[1] - p.z[0]));
                    b.bin(if p.emax { BinOp::Emax } else { BinOp::Max }, profile, span)
                })
                .collect();
            let f = parts.into_iter().reduce(|p, q| b.bin(BinOp::Min, p, q)).unwrap();
            Arc::new(b.tape(f))
        };
        let (a, c) = (field(&pa), field(&pb));
        // A box holding every prism of a shape.
        let around = |ps: &[Prism]| -> (V3, V3) {
            let reach = |p: &Prism| p.disc.map_or(p.half, |r| [r, r]);
            let lo = ps.iter().fold([f64::INFINITY; 3], |l, p| [l[0].min(p.c[0] - reach(p)[0]), l[1].min(p.c[1] - reach(p)[1]), l[2].min(p.z[0])]);
            let hi = ps.iter().fold([f64::NEG_INFINITY; 3], |h, p| [h[0].max(p.c[0] + reach(p)[0]), h[1].max(p.c[1] + reach(p)[1]), h[2].max(p.z[1])]);
            (lo.map(|v| v - 0.5), hi.map(|v| v + 0.5))
        };
        let ((alo, ahi), (blo, bhi)) = (around(&pa), around(&pb));
        let ways = [Way { grown: a.clone(), covered: c.clone(), lo: blo, hi: bhi, swapped: false }, Way { grown: c, covered: a, lo: alo, hi: ahi, swapped: true }];
        match clearance_now(&ways, &Ctx::new(0.0, vec![]), d, 0.0125, BOXES).0 {
            Verdict::Proven => {
                assert!(truth >= d - 1e-9, "case {case}: {d} proven, the shapes {truth} apart");
                proven += 1;
            }
            Verdict::Disproven { value, .. } => {
                assert!(truth < d + 1e-9 && value < d, "case {case}: {d} found under by {value}, the shapes {truth} apart");
                found += 1;
            }
            Verdict::Unknown { .. } => {}
        }
    }
    assert!(proven > 15 && found > 15, "{proven} proven, {found} found");
}

/// An unknown clearance holds to the accuracy it attained: two balls
/// 0.5 apart, their fields squared (no slope bound of 1, so a box's
/// bound proves less than its distance), asked just under and just
/// over the gap with coarse cells. Over it, never proven; where
/// unknown, the clearance less its attained accuracy is no more than
/// the gap.
#[test]
fn an_unknown_clearance_holds_to_the_accuracy_it_attained() {
    let mut b = crate::test_fields::B::new();
    let (p, q) = (ball(&mut b, 0.0, 2.0, true), ball(&mut b, 4.5, 2.0, true));
    let (p, q) = (Arc::new(b.tape(p)), Arc::new(b.tape(q)));
    let mut attained = 0;
    for d in [0.45, 0.499, 0.501, 0.55] {
        let way = Way { grown: p.clone(), covered: q.clone(), lo: [2.0, -2.0, -2.0], hi: [6.5, 2.0, 2.0], swapped: false };
        match clearance_now(&[way], &Ctx::new(0.0, vec![]), d, 0.05, BOXES).0 {
            Verdict::Proven => assert!(d <= 0.5, "{d} proven past the gap"),
            Verdict::Unknown { attained: Some(e), .. } => {
                assert!(d - e <= 0.5, "{d} holds to {e}: past the gap");
                attained += 1;
            }
            v => assert!(d > 0.5 || matches!(v, Verdict::Unknown { .. }), "{d}: {v:?}"),
        }
    }
    assert!(attained > 0, "no unknown attained an accuracy");
}

#[test]
fn a_clearance_that_does_not_read_t_has_its_witness_at_the_moment_asked() {
    // Two balls 1 apart, neither moving: one answer at every moment of
    // the motion (its key holds no t), its witness at the moment asked
    // (a cached one carried the moment it was made at).
    let mut g = Graph::new();
    let ball = |g: &mut Graph, cx: f64| {
        let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
        let c = g.constant(cx);
        let dx = g.sub(x, c);
        let (xx, yy, zz) = (g.mul(dx, dx), g.mul(y, y), g.mul(z, z));
        let s = g.add(xx, yy);
        let s = g.add(s, zz);
        let r = g.un(ink_core::graph::UnOp::Sqrt, s);
        let one = g.constant(1.0);
        g.sub(r, one)
    };
    let (a, b) = (ball(&mut g, -1.5), ball(&mut g, 1.5));
    let (a, b) = (Arc::new(Tape::new(&g, a)), Arc::new(Tape::new(&g, b)));
    for t in [0.25, 0.75] {
        let ctx = Ctx { t_span: Some((0.0, 1.0)), ..Ctx::new(t, vec![]) };
        let way = Way { grown: a.clone(), covered: b.clone(), lo: [0.0, -1.5, -1.5], hi: [3.0, 1.5, 1.5], swapped: false };
        match clearance(&[way], &ctx, 2.0, 0.25, BOXES).0 {
            Verdict::Disproven { t: at, .. } => assert_eq!(at, t),
            v => panic!("{v:?}"),
        }
    }
}

#[test]
fn witnesses_carry_finite_in_span_moments() {
    // 0.5 (t0 + t1) would overflow for t0 = 1e308, t1 = 1.1e308.
    let m = span_mid((1e308, 1.1e308));
    assert!((1e308..=1.1e308).contains(&m), "{m}");
    // 0.5 - t is below 0 over [1, 2]: a whole-cell witness, which must
    // carry a moment of the span, not ctx.t = 0 (where it is 0.5).
    let mut g = Graph::new();
    let (h, t) = (g.constant(0.5), g.var(Var::T));
    let f = g.sub(h, t);
    let tape = Arc::new(Tape::new(&g, f));
    let ctx = Ctx { t_span: Some((1.0, 2.0)), ..Ctx::new(0.0, vec![]) };
    match prove_empty(&tape, &ctx, [-1.0; 3], [1.0; 3], Some(0.5), 0.5, 0.0).verdict {
        Verdict::Disproven { t, witness, .. } => assert!((1.0..=2.0).contains(&t) && witness.iter().all(|x| x.is_finite()), "{t} {witness:?}"),
        v => panic!("{v:?}"),
    }
}

#[test]
fn an_inside_pair_with_a_target_without_interior_ends() {
    // A field 0 over its box (an agent's `(move 0 0 5 (box 4 4 1))`,
    // arguments swapped): a mesh pass finds no surface, and each finer
    // one would test every cell, to 32 768 along the box (past 30 s on
    // 28 cores). The search proves it holds no point: no pair. Its cost
    // is the CPU of a process searching alone (test_fields::alone), not
    // the time it takes, which a machine busy with other work stretches.
    if std::env::var_os("INK_INSIDE_PAIR_ALONE").is_some() {
        let mut b = crate::test_fields::B::new();
        let slab = b.cuboid([0.0; 3], [5.0, 5.0, 1.0]);
        let zero = b.c(0.0);
        let (ta, tb) = (Arc::new(b.tape(zero)), Arc::new(b.tape(slab)));
        let c = inside_pair(&ta, &tb, &Ctx::new(0.0, vec![]), ([-2.0, -2.0, 4.0], [2.0, 2.0, 6.0]), ([-5.0, -5.0, -1.0], [5.0, 5.0, 1.0]), 1e-3);
        assert!(c.is_none(), "{c:?}");
        return;
    }
    // Each finer pass testing every cell took past 840 CPU s (30 s on 28
    // cores); the proof takes about 0.03 (2026-11-17).
    let (ok, cpu) = crate::test_fields::alone("query::tests::an_inside_pair_with_a_target_without_interior_ends", "INK_INSIDE_PAIR_ALONE");
    assert!(ok && cpu < 30.0, "{cpu:.1} CPU s (ok: {ok})");
}

#[test]
fn a_proof_refines_past_a_pass_with_many_undecided_cells() {
    // A rod (radius 20, 100 long) in a tube 0.002 clear of it: cells of
    // 0.256 leave the gap undecided all round (over 200 000 of them),
    // cells of 0.128 separate the faces. Refining stopped after a pass
    // of 100 000 undecided cells, unknown; the search tests up to BOXES
    // boxes, which this is well within. The tube's radius is written
    // apart (its squares): of one radius, the gap is exact at once (Same).
    use ink_core::graph::{BinOp, UnOp};
    let mut b = crate::test_fields::B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let (xx, yy) = (b.bin(BinOp::Mul, x, x), b.bin(BinOp::Mul, y, y));
    let s = b.bin(BinOp::Add, xx, yy);
    let r = b.un(UnOp::Sqrt, s);
    let (sx, sy) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y));
    let s2 = b.bin(BinOp::Add, sx, sy);
    let r2 = b.un(UnOp::Sqrt, s2);
    let (inner, outer) = (b.c(20.0), b.c(20.002));
    let rod = b.bin(BinOp::Sub, r, inner);
    let tube = b.bin(BinOp::Sub, outer, r2);
    let (span, around) = (b.slab(z, 0.0, 50.0), b.cuboid([0.0; 3], [30.0, 30.0, 50.0]));
    let rod = b.bin(BinOp::Max, rod, span);
    let tube = b.bin(BinOp::Max, tube, around);
    let both = b.bin(BinOp::Max, rod, tube);
    let tape = Arc::new(b.tape(both));
    let s = prove_empty(&tape, &Ctx::new(0.0, vec![]), [-31.0, -31.0, -51.0], [31.0, 31.0, 51.0], None, 0.001, 0.0);
    assert!(matches!(s.verdict, Verdict::Proven) && s.res < 0.2 && s.exact, "{s:?}");
}

#[test]
fn a_witness_lies_in_the_region_asked() {
    // max(x - 0.007, -19 x - 0.0075), in a unit cube's slabs: material
    // 0.007 deep at the region's face x = 0. The first term dominates
    // over the region (x in [0, 1]), so a tape focused on it is x -
    // 0.007 alone, which below x = 0 reads ever deeper material where
    // the field is above 0 (from x = -0.0004). A root cell reaching
    // below the region, which that tape reads full, would give the
    // deepest witness (-0.044, 0.956, 0.056), outside the region and the
    // shape. Cells are clipped to the region.
    use ink_core::graph::BinOp;
    let mut b = crate::test_fields::B::new();
    let x = b.v(Var::X);
    let (k3, k19, k31) = (b.c(0.007), b.c(-19.0), b.c(0.0075));
    let a = b.bin(BinOp::Sub, x, k3);
    let steep = b.bin(BinOp::Mul, k19, x);
    let steep = b.bin(BinOp::Sub, steep, k31);
    let f = b.bin(BinOp::Max, a, steep);
    let (y, z) = (b.v(Var::Y), b.v(Var::Z));
    let (sy, sz) = (b.slab(y, 0.5, 0.5), b.slab(z, 0.5, 0.5));
    let f = b.max3(f, sy, sz);
    let tape = Arc::new(b.tape(f));
    let ctx = Ctx::new(0.0, vec![]);
    let opts = crate::ranged::ProofOptions { res: Some(0.05), finest: 0.05, depth: 0, tol: 0.0, boxes: crate::query::BOXES };
    let (lo, hi) = ([0.0; 3], [1.0; 3]);
    match crate::ranged::prove(&tape, &ctx, lo, hi, &opts).verdict {
        crate::query::Verdict::Disproven { witness, .. } => {
            assert!((0..3).all(|k| (lo[k]..=hi[k]).contains(&witness[k])), "{witness:?}");
            assert!(eval_point(&tape, &ctx, witness) < 0.0, "{witness:?}");
        }
        o => panic!("{o:?}"),
    }
}

/// A ray finds the material behind a select whose condition straddles 0
/// by rounding (cos of pi/2 y at y = 1): its bounds there are 10 wide,
/// the branches' gap, not rounding (a rule taking them for rounding
/// would drop it), and the material (x in 0.5
/// .. 1.25 at y = 1, f -0.25 at x = 1) is found where it starts,
/// unsure: only rounding decides the select there, so bounds show no
/// point of it inside.
#[test]
fn a_ray_finds_material_behind_an_undecided_select() {
    use crate::test_fields::B;
    use ink_core::graph::{BinOp, UnOp};
    let mut b = B::new();
    let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
    let (half, k) = (b.c(0.5), b.c(std::f64::consts::FRAC_PI_2));
    let ky = b.bin(BinOp::Mul, k, y);
    let cond = b.un(UnOp::Cos, ky);
    let (ten, zero) = (b.c(10.0), b.c(0.0));
    let sel = b.0.select(cond, ten, zero);
    let ramp = b.bin(BinOp::Sub, half, x);
    let a = b.bin(BinOp::Add, ramp, sel);
    let (sx, sy, sz) = (b.slab(x, 0.75, 0.5), b.slab(y, 1.0, 1.0), b.slab(z, 0.0, 1.0));
    let m = b.max3(sx, sy, sz);
    let f = b.bin(BinOp::Max, a, m);
    let tape = Arc::new(b.tape(f));
    let ctx = Ctx::new(0.0, vec![]);
    let t = ray_to(&Local::whole(&tape, &ctx), [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], (0.0, 2.0), Side::Inside, 2e-9);
    assert!(t.is_some_and(|h| (h.t - 0.5).abs() < 1e-6 && h.rounding.is_some()), "{t:?}");
}

#[test]
fn ray_segment_boxes_hold_the_real_segment() {
    // o + t dir correctly rounded (fma) lies in the box; nearest ends
    // (a rounded product, then a rounded sum) often miss it.
    let mut r = 0x9e3779b97f4a7c15u64;
    let mut f = || {
        r ^= r << 13;
        r ^= r >> 7;
        r ^= r << 17;
        (r >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    };
    for _ in 0..10000 {
        let (o, dir, t) = ([f(), f(), f()], [f(), f(), f()], f().abs() * 10.0);
        let b = segment_box(o, dir, t, t);
        for k in 0..3 {
            assert!(b[k].contains(t.mul_add(dir[k], o[k])), "{o:?} {dir:?} {t}");
        }
    }
}
