//! End-to-end kernel tests: evaluator agreement, interval soundness, meshing
//! accuracy, integrals, slices, proofs, assemblies, and error reporting.

use ink::centred::{centred_at_least, centred_lower, centred_lower_max};
use ink::eval::*;
use ink::interval::I;
use ink::lang::Doc;
use ink::mesh;
use ink::query::{self, Verdict};
use ink::tape::Tape;
use std::f64::consts::PI;
use std::path::Path;
use std::sync::Arc;

/// The repository: what the tests load (catalog parts, fixtures) is found from it.
const ROOT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

fn doc(src: &str) -> Doc {
    ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"))
}

/// A published mass is of the params' defaults: a size param scales it
/// with the volume. Over t the body moves (a pose), measured at each
/// moment, so the mass stays.
#[test]
fn published_mass_derives_density_at_the_defaults_and_the_moment() {
    let mut d = doc("(param w 20) (body b (shape (box w 10 (+ 10 t))) (attr mass .199kg))
        (assembly root (instance a b) (instance b b (at (move 40 0 0))))
        (assert light (integral root (attr density) < 410))");
    let items = d.items("root").unwrap();
    // The claim reads the density: provided at every context (Doc::want).
    for (w, time) in [(20.0, 0.0), (20.0, 5.0), (40.0, 5.0)] {
        let ctx = d.ctx(&[("w".into(), w.into())], time).unwrap();
        for it in &items {
            let (rho, why) = ink::materials::mass_density(&mut d, &format!("'{}'", it.label), &it.attrs, None).unwrap().unwrap();
            assert!(why.starts_with("published mass of 'b' over its volume") && why.contains("centroid"), "{why}");
            let density = eval_point(&Tape::new(&d.g, rho), &ctx, [0.0; 3]);
            let mass = density * w * 10.0 * (10.0 + time);
            assert!((mass / (199.0 * w / 20.0) - 1.0).abs() < 0.01, "{mass}");
            assert!(d.inexact(&[rho], &[]).any(), "measured volume must not become an exact density");
        }
    }
    for bad in ["-1", "x", "\"199\"", "2mm"] {
        assert!(ink::load_str(&format!("(import \"std\") (body b (shape (box 10 10 10)) (attr mass {bad}))"), "test", Path::new(ROOT)).is_err(), "{bad}");
    }
    for other in ["(attr density 1g/cm3)", "(attr material \"steel\")"] {
        assert!(ink::load_str(&format!("(import \"std\") (body b (shape (box 10 10 10)) (attr mass 199) {other})"), "test", Path::new(ROOT)).is_err());
    }
    // A mass over what has no bounded volume, or lighter than what its
    // bodies with their own density weigh, is refused where it is weighed,
    // and by check (ModelView::mass_refusals), not at load. Check finds it
    // without measuring the volume the density shares (most of a small
    // part's check): a mass it does not refuse is left unmeasured, and one
    // whose own-density bodies could weigh it (their density's bound over
    // their box) is measured there to decide.
    for (src, why) in [
        ("(body b (shape (box 10 10 10)) (attr mass 199))", ""),
        // A box with no point inside (the field touches 0 at the origin
        // alone): no volume, which bounds alone do not show.
        ("(body h (shape (max (box 10 10 10) (abs x) (abs y) (abs z)))) (assembly root (attr mass 199) (instance h h))", "positive finite volume"),
        (
            "(body s (shape (box 20 20 20)) (attr density 7.85g/cm3)) (body c (shape (move (box 5 5 5) 0 0 20)))
                         (assembly root (attr mass 70g) (instance s s) (instance c c))",
            "",
        ),
        ("(body b (shape x) (attr mass 199))", "needs bounded bodies"),
        (
            "(body s (shape (box 20 20 20)) (attr density 7.85g/cm3)) (body c (shape (move (box 5 5 5) 0 0 20)))
                         (assembly root (attr mass 20g) (instance s s) (instance c c))",
            "more than its published 20 g",
        ),
    ] {
        let d = doc(src);
        let ctx = d.ctx(&[], 0.0).unwrap();
        let mut v = ink::view::ModelView::new(d, ctx);
        let failed = v.mass_refusals();
        let density = v.doc.published.last().unwrap().0;
        if why.is_empty() {
            assert!(failed.is_empty() && v.doc.unmeasured_param(density), "{src}: {failed:?}");
        } else {
            assert!(failed.len() == 1 && failed[0].1.message.contains(why), "{failed:?}");
        }
    }
}

#[test]
fn published_mass_volume_refines_to_its_numerical_accuracy() {
    // A grid of extent/200 overestimates this thin annulus's volume enough
    // to put its density 14% low. A sphere checks a curved solid.
    let mut d = doc("(body ball (shape (sphere 10)) (attr mass 199g))
        (body ring (shape (subtract (cylinder-z 10 .1) (cylinder-z 9.9 1))) (attr mass 1g))");
    for (name, expected) in [("ball", 199.0 / (4.0 / 3.0 * PI * 1000.0)), ("ring", 1.0 / (PI * (100.0 - 9.9f64.powi(2)) * 0.1))] {
        let attrs = d.items(name).unwrap().remove(0).attrs;
        let (rho, why) = ink::materials::mass_density(&mut d, &format!("'{name}'"), &attrs, None).unwrap().unwrap();
        d.want(&[rho]);
        let density = d.scalar(rho, &d.ctx(&[], 0.0).unwrap());
        assert!((density / expected - 1.0).abs() < 0.01, "{name}: {density} vs {expected}");
        assert!(why.contains("measured to 2%"), "{why}");
    }
}

#[test]
fn numerical_claims_cover_only_the_scope_their_methods_evaluate() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    use ink::view::ModelView;
    let src = "(def h (- 2 t)) (def plate (box 4 4 h))
        (body p (shape plate))
        (assert thick (min-thickness plate > 1.5))
        (assert fixed (min-thickness (box 4 4 2) > 1.5))";
    let d = doc(src);
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    let opts = ClaimOptions { res: Some(0.5), ..Default::default() };
    let r = run_asserts(&mut v, None, false, opts).unwrap();
    assert!(r.all_hold(), "{r:?}");
    // Over a motion they are measured at each step: the plate thins below
    // 1.5 after t = 0.5, which the steps past it see.
    v.set_motion(Some((0.0, 1.0))).unwrap();
    let r = run_asserts(&mut v, None, false, opts).unwrap();
    // Estimates (no tolerance), a failure at a step too: rays guarantee none.
    let a = &r.asserts;
    assert!(a[0].check.status == Status::Fails && a[0].check.accuracy.is_none(), "{:?}", a[0]);
    assert!(a[1].check.status == Status::Holds && a[1].check.accuracy.is_none(), "{:?}", a[1]);
    assert!(r.asserts.iter().all(|a| a.check.motion.as_ref().is_none_or(|m| !m.bound)));
    // Runtime rays can evaluate at t = 1.
    v.set_motion(None).unwrap();
    let ctx = v.doc.ctx(&[], 1.0).unwrap();
    let mut v = ModelView::new(v.doc, ctx);
    let r = run_asserts(&mut v, None, false, opts).unwrap();
    assert!(r.asserts[0].check.status == Status::Fails && r.asserts[0].check.accuracy.is_none(), "{:?}", r.asserts[0]);
}

/// Enclosed voids are counted by proof, never from samples: a void sealed
/// by a skin a fifth of the first mesh's cells thick is found however
/// steep its field (the mesh read 0 voids for `(* 30 skinned)`, its corner
/// values past its distance test) and wherever the skin falls between cell
/// centres (`big`: the mesh read 0); the count is exact where it decides
/// (two voids; none once the hole is left open; none past a 0.2 mm slit at
/// (res 1)); and it is bounds where cells cannot part two airs: a hole
/// flush with the void's ceiling (f = 0 across its mouth) and a bearing's
/// balls seated in its races' grooves (surfaces that touch) are unknown,
/// as a count over a ranged param and of a shape not proven bounded are.
#[test]
fn cavities_are_counted_by_proof() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    use ink::view::ModelView;
    let d = doc("(param wall 3 (tol 0 0.5))
        (def outer (move (box 100 100 100) 0 0 50))
        (def opened (subtract outer (move (box 94 94 94) 0 0 50) (extrude-z (circle-2d 5) 96 101)))
        (def flush (subtract outer (move (box 94 94 94) 0 0 50) (extrude-z (circle-2d 5) 97 101)))
        (def skinned (min opened (extrude-z (circle-2d 6) 98.5 98.8)))
        (def big (min (subtract (move (box 200 200 200) 0 0 100) (move (box 190 190 190) 0 0 100) (extrude-z (circle-2d 10) 190 201))
                        (extrude-z (circle-2d 12) 196 196.3)))
        (def two (subtract (box 20 20 20) (move (box 6 16 16) -4 0 0) (move (box 6 16 16) 4 0 0)))
        (def slit (subtract (box 20 20 20) (box 12 12 12) (move (box 0.2 4 20) 0.5 0 10)))
        (assert sealed-behind-skin (cavities skinned < 1))
        (assert steep-empties (cavities (* 30 skinned) < 1))
        (assert skin-between-centres (cavities big < 1))
        (assert left-open (cavities opened < 1))
        (assert two-voids (cavities two < 3))
        (assert not-one (cavities two < 2))
        (assert slit-at-res-1 (cavities slit < 1) (res 1))
        (assert flush-mouth (cavities flush < 1))
        (def races (subtract (extrude-z (ring-2d 7 13) -2 2) (extrude-z (ring-2d 9.5 10.5) -3 3) (torus-z 10 2)))
        (assert seated-balls (cavities (min races (pattern-polar (move (sphere 2) 10 0 0) 8)) < 1))
        (assert ranged (cavities (subtract (box 20 20 20) (box (* 2 wall) 4 4)) < 1))
        (assert unbounded (cavities (subtract (- z) (sphere 2)) < 1))");
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    let spans = v.declared_spans().unwrap();
    v.set_param_spans(&spans).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let s: Vec<(&str, Status)> = r.asserts.iter().map(|a| (a.name.as_str(), a.check.status)).collect();
    use Status::{Fails, Holds, Unknown};
    let want = [Fails, Fails, Fails, Holds, Holds, Fails, Holds, Unknown, Unknown, Unknown, Unknown];
    assert_eq!(s.iter().map(|x| x.1).collect::<Vec<_>>(), want, "{s:?}");
    // A count proven is an exact answer, not an estimate.
    assert!(r.asserts[..7].iter().all(|a| a.check.accuracy.is_some()) && r.asserts[4].check.holds_guaranteed(), "{:?}", r.asserts[4].check);
    let ink::claims::Evidence::Metric { measured, .. } = &r.asserts[7].check.evidence else { panic!("{:?}", r.asserts[7]) };
    assert!(measured.value == Ok(0.0) && measured.excluded.iter().any(|e| e.key == "undecided-cells" && e.value >= 1.0), "{measured:?}");
}

/// A void count is the one search's (crate::search): what it answers is a
/// function of the question, not of the threads, where it runs out of its
/// budget too (a hole flush with the void's ceiling, f = 0 across its
/// mouth, halves to the budget), and says so: the budget, not the res, is
/// what stopped it, with the cells it left and what to do.
#[test]
fn a_void_count_is_the_same_on_any_threads_and_says_its_budget_stopped_it() {
    use ink::view::ModelView;
    let d = doc("(def flush (subtract (move (box 100 100 100) 0 0 50) (move (box 94 94 94) 0 0 50) (extrude-z (circle-2d 5) 97 101)))");
    let ink::lang::Value::Field(f, _) = *d.def("flush").unwrap() else { panic!("a field") };
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    let (lo, hi, _) = v.bounds(f).unwrap();
    let (tape, ctx) = (v.clipped(f), v.ctx.clone());
    let h = ink::query::max_extent(lo, hi) / ink::view::COMPONENT_CELLS;
    let count = |threads: usize| rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap().install(|| ink::voids::voids(&tape, &ctx, lo, hi, h));
    let (one, all) = (count(1), count(8));
    assert_eq!(format!("{one:?}"), format!("{all:?}"));
    assert!(one.count().is_none() && one.ended == Some(ink::query::Ended::Limit) && one.undecided.is_some_and(|(a, b)| a >= one.finest && b > a), "{one:?}");
    let (what, decide) = one.why_undecided().unwrap();
    assert!(what.contains("budget") && decide.contains("a finer res will not help"), "{what} / {decide}");
}

/// What a view learned of a field holds for the scope it learned it over:
/// a belt growing with t, its bounds found at the view's moment (t = 0,
/// x within 1.1), keeps its box over the motion (x within 11.1 by t = 1),
/// so the post it reaches then stays its candidate, not set apart by their
/// boxes with the moment's box read as the motion's. The view's own box is
/// the moment's after.
#[test]
fn a_box_over_a_motion_is_not_a_moments_found_before() {
    use ink::view::ModelView;
    let d = doc("(body belt (shape (box (+ 2 (* 20 t)) 2 2)))
        (body post (shape (move (box 2 2 2) 8 0 0)))
        (assembly root (instance belt belt) (instance post post))");
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    v.set_motion(Some((0.0, 1.0))).unwrap();
    let items = v.doc.items("root").unwrap();
    let f = v.doc.placed_field(&items[0]);
    let moment = v.bounds(f).unwrap();
    let [belt, post]: [ink::motion::ItemBox; 2] = ink::motion::item_boxes(&mut v, &items[..2]).try_into().unwrap();
    let (candidates, apart, _) = ink::lang::pairs::candidate_pairs(&[belt.0, post.0], &[belt.1, post.1], &[belt.2.clone(), post.2.clone()], &[]);
    assert!(moment.1[0] < 1.2 && belt.1 && belt.0.is_some_and(|b| b.1[0] >= 11.0) && candidates == vec![(0, 1)] && apart == 0, "{moment:?} {belt:?} {candidates:?}");
    assert_eq!(v.bounds(f), Some(moment));
}

/// A mass is settled only where a halving showed it settling: a 0.3 mm
/// layer of density in a 100 mm cube (29.2 g) is not settled at 42 g with
/// the last halving having moved it 76%, nor a layer every first cell
/// misses settled at zero without a halving.
#[test]
fn a_mass_is_settled_only_by_a_halving_that_shows_it() {
    for (src, mass) in [("(+ .00001 (select (max (- (abs x) .15) (- (abs y) 40) (- (abs z) 40)) .01 0))", 29.2), ("(select (- (abs x) .15) .001 0)", 3.0)] {
        let (mut v, f) = view_of(&format!("(def f (box 100 100 100)) (def rho {src})"));
        let ink::lang::Value::Field(rho, _) = *v.doc.def("rho").unwrap() else { panic!("rho is a field") };
        let m = v.moments(f, Some(rho), None).unwrap();
        assert!(m.coarse.is_some() && (!m.settled || (m.sums[0] / mass - 1.0).abs() <= 0.05), "{src}: {m:?}");
        assert!(m.settled || m.unsettled().is_some());
    }
}

/// A view and its def `f` (a field), at t = 0.
fn view_of(src: &str) -> (ink::view::ModelView, u32) {
    let d = doc(src);
    let f = match d.def("f").unwrap() {
        ink::lang::Value::Field(f, _) => *f,
        _ => panic!("f is a field"),
    };
    let ctx = d.ctx(&[], 0.0).unwrap();
    (ink::view::ModelView::new(d, ctx), f)
}

/// A view keeps what it found of a field by the scope its inputs had
/// so another finds it again, and the view's own is
/// back after: a param ranged over [1, 1] reads (p + 1e16) - 1e16 as 1, at
/// its value 1 as the float 0 (a box 22 wide, not 2); a param's value
/// moved within its span (a box 22 long, 88 mm3, not 8); a
/// solve on t, read through a solve it sets, at another moment; a param a
/// solve derives over a span; an import's param bound to one of the file's;
/// a density reading a param. Each is its own scope's answer, not the
/// first scope's.
#[test]
fn a_views_answers_are_kept_by_the_scope_they_were_found_over() {
    let width = |v: &mut ink::view::ModelView, f: u32| v.bounds(f).unwrap().1[0];
    let (mut v, f) = view_of("(param p 1) (def delta (- (+ p 10000000000000000) 10000000000000000)) (def f (box (+ 2 (* 20 delta)) 2 2))");
    let point = width(&mut v, f);
    let ranged = v.at_point(Ctx { p_spans: ParamSpans::new([(0, (1.0, 1.0))]), ..v.ctx.clone() }, |v| width(v, f));
    assert!(point < 1.2 && ranged >= 11.0 && width(&mut v, f) == point, "{point} {ranged}");
    let (mut v, f) = view_of("(param p 2) (def f (box p 2 2))");
    v.set_param_spans(&[(0, (2.0, 22.0))]).unwrap();
    let mass = |v: &mut ink::view::ModelView| v.moments(f, None, Some(0.5)).unwrap().sums[0];
    let short = mass(&mut v);
    let long = v.at_point(Ctx { params: vec![22.0], ..v.ctx.clone() }, |v| mass(v));
    assert!((short / 8.0 - 1.0).abs() < 0.05 && (long / 88.0 - 1.0).abs() < 0.05 && mass(&mut v) == short, "{short} {long}");
    let (mut v, f) = view_of("(solve first (unknowns (a 0)) (== a (* 10 t))) (solve second (unknowns (b 2)) (== b (+ 2 (* 2 a)))) (def f (box b 2 2))");
    let before = width(&mut v, f);
    let at = v.doc.ctx(&[], 1.0).unwrap();
    let after = v.at_point(at, |v| width(v, f));
    assert!(before < 1.2 && after >= 11.0 && width(&mut v, f) == before, "{before} {after}");
    let (mut v, f) = view_of("(param p 0) (solve first (unknowns (a 0)) (== a (* 10 p))) (solve second (unknowns (b 2)) (== b (+ 2 (* 2 a)))) (def f (box b 2 2))");
    let before = width(&mut v, f);
    let spans = v.scope_spans(&[(0, (0.0, 1.0))]).unwrap();
    let after = v.at_point(Ctx { p_spans: spans, ..v.ctx.clone() }, |v| width(v, f));
    assert!(before < 1.2 && after >= 11.0 && width(&mut v, f) == before, "{before} {after}");
    let (mut v, f) = view_of("(param p 0) (def f (box 2 2 2)) (def rho (+ 1 (* 9 p)))");
    let rho = match v.doc.def("rho").unwrap() {
        ink::lang::Value::Field(f, _) => *f,
        _ => panic!("rho is a field"),
    };
    let weigh = |v: &mut ink::view::ModelView| v.moments(f, Some(rho), Some(0.5)).unwrap().sums[0];
    let light = weigh(&mut v);
    let heavy = v.at_point(Ctx::new(0.0, vec![1.0]), |v| weigh(v));
    assert!((heavy / light - 10.0).abs() < 0.001 && weigh(&mut v) == light, "{light} {heavy}");
}

/// Whether a constants-and-params value is taken in floats rests on the
/// params it reads alone: with param 64 ranged,
/// param 0's (p + 1e16) - 1e16 stays the float 0, a box 2 wide, as with
/// nothing ranged; with param 0 itself ranged over [1, 1] it is 1 in real
/// numbers, 22 wide. Ranging param 64 does not take it off the float path
/// (the bit k % 64 they share stands for both), so the answer for one
/// context does not rest on what an earlier one found. Each in a view of
/// its own.
#[test]
fn a_value_is_taken_in_floats_by_the_params_it_reads_alone() {
    let mut src = "(param p 1) ".to_string();
    for k in 1..=64 {
        src += &format!("(param q{k} 0) ");
    }
    src += "(def delta (- (+ p 10000000000000000) 10000000000000000)) (def f (box (+ 2 (* 20 delta)) 2 2))";
    let width = |ranged: Option<u32>| {
        let (mut v, f) = view_of(&src);
        let spans = ParamSpans::new(ranged.map(|k| (k, (v.ctx.params[k as usize], v.ctx.params[k as usize]))));
        v.at_point(Ctx { p_spans: spans, ..v.ctx.clone() }, |v| v.bounds(f).unwrap().1[0])
    };
    let (none, other, own) = (width(None), width(Some(64)), width(Some(0)));
    assert!(none < 1.2 && other == none && own >= 11.0, "{none} {other} {own}");
}

/// Deterministic pseudo-random numbers in [0, 1).
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next()
    }
}

/// A zoo of fields exercising every op.
const ZOO: &[&str] = &[
    "(sphere 10)",
    "(rotate-z (rotate-x (box 20 30 10) 30) 20)",
    "(shell (gyroid 7) 1)",
    "(smooth-union (sphere 6) (move (box 8 8 8) 5 0 0) 2)",
    "(twist-z (box 10 4 30) 6)",
    "(repeat-x (cylinder-z 2 10) 7)",
    "(torus-z 10 3)",
    "(select (- z 1) (sphere 8) (box 12 12 12))",
    "(+ (atan2 y x) (pow (abs x) 1.5) (exp (* -0.1 (square z))) (ln (+ 1 (abs y))) (tan (* 0.1 x)) (floor (* 0.3 y)))",
    "(/ (- (sqrt (+ (square x) (square y))) 5) (+ 2 (cos z)))",
    "(pow y 3)",
    "(atan (* x y))",
    "(grid-bake (subtract (box 20 30 10) (cylinder-z 5 20)) (res 0.7))",
    "(move (rotate-z (grid-bake (torus-z 8 2) (res 0.5)) 30) 3 -2 1)",
];

/// `(check TEST value)`: a condition any builder states on what it builds,
/// which `ink check` reports where it fails, with where it was stated. The
/// same test stated twice is one condition; two tests on the same values
/// are two. Every check loading runs is reported, an imported file's too,
/// at the importer's values: a body built from an imported definition
/// through a warp (`move`) fails where the definition does.
#[test]
fn checks_report_what_fails_where_it_is_built() {
    let failures = |d: ink::lang::Doc| {
        let ctx = d.ctx(&[], 0.0).unwrap();
        let v = ink::view::ModelView::new(d, ctx);
        v.doc.check_failures(&v.ctx, v.accuracy(), &[])
    };
    let w = failures(doc("(fn plate (h) (check (value h >= 2) (box 10 10 h))) (param th 1) (body p (shape (plate th))) (body q (shape (move (plate th) 20 0 0)))"));
    assert!(w.len() == 1 && w[0].contains("check (value h >= 2) fails: 1 >= 2 is false") && w[0].contains("plate (test:"), "{w:?}");
    let w = failures(doc("(param h 1) (body p (shape (check (value h > 2) (check (value h >= 2) (box 10 10 h)))))"));
    assert!(w.len() == 2 && w.iter().any(|w| w.contains("h > 2")) && w.iter().any(|w| w.contains("h >= 2")), "{w:?}");
    let dir = crate::common::Scratch::new("check");
    std::fs::write(dir.join("lib.ink"), "(import \"std\") (param r 0.5) (def ball (check (value r <= 1) (sphere r))) (body b (shape (move ball 5 0 0)))").unwrap();
    std::fs::write(dir.join("main.ink"), "(import \"lib.ink\" as L) (body mine (shape (move L.ball 3 0 0)))").unwrap();
    let load = |f: &str, r: f64| {
        let params = [(if f == "lib.ink" { "r" } else { "L.r" }.to_string(), r.into())];
        let d = ink::load_str_with(&std::fs::read_to_string(dir.join(f)).unwrap(), f, &dir, &params).unwrap();
        let ctx = d.ctx(&params, 0.0).unwrap();
        let v = ink::view::ModelView::new(d, ctx);
        v.doc.check_failures(&v.ctx, v.accuracy(), &[])
    };
    assert_eq!(load("main.ink", 0.5), Vec::<String>::new());
    let w = load("main.ink", 2.0);
    assert!(w.len() == 1 && w[0].contains("lib.ink:1 (def ball") && w[0].contains("L.r = 2"), "{w:?}");
    assert_eq!(load("lib.ink", 2.0).len(), 1);
}

fn tape_of(src: &str) -> (Doc, Tape) {
    let d = doc(&format!("(def f {src})"));
    let ink::lang::Value::Field(f, _) = d.def("f").unwrap().clone() else { panic!() };
    let t = Tape::new(&d.g, f);
    (d, t)
}

#[test]
fn built_in_functions_are_values() {
    // A form that is a function is passed like a user fn: polygon over a
    // list of outlines is the polygon of each.
    let v = |src: &str, p: [f64; 3]| {
        let (d, t) = tape_of(src);
        eval_point(&t, &Ctx::new(0.0, d.params.iter().map(|p| p.default).collect()), p)
    };
    let squares = "(vec (vec (vec -1 -1) (vec 1 -1) (vec 1 1) (vec -1 1)) (vec (vec 4 -1) (vec 6 -1) (vec 6 1) (vec 4 1)))";
    assert_eq!(v(&format!("(reduce min (map polygon {squares}))"), [5.0, 0.0, 0.0]), -1.0);
    assert_eq!(v(&format!("(apply min (map polygon {squares}))"), [0.0, 0.0, 0.0]), -1.0);
    assert_eq!(v("(apply + (map len (vec (vec 1 2) (vec 3))))", [0.0; 3]), 3.0);
    assert_eq!(v("(fold (lambda (a f) (+ a (f 4))) 0 (vec square (lambda (u) (nth (vec 7 u) 1))))", [0.0; 3]), 20.0);
    // A form that is syntax says so, and how to pass it.
    let e = ink::load_str("(import \"std\") (def v (map if (vec 1) (vec 2) (vec 3)))", "test", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("'if' is a form, not a function") && e.hint.as_deref().unwrap_or("").contains("(lambda"), "{e:?}");
}

#[test]
fn intervals_are_sound_and_simplification_is_exact() {
    let ctx = Ctx::new(0.0, vec![]);
    let mut rng = Lcg(7);
    for src in ZOO {
        let (_d, tape) = tape_of(src);
        for _ in 0..300 {
            let c = [rng.range(-20.0, 20.0), rng.range(-20.0, 20.0), rng.range(-20.0, 20.0)];
            let s = rng.range(0.01, 8.0);
            let b = [I::new(c[0] - s, c[0] + s), I::new(c[1] - s, c[1] + s), I::new(c[2] - s, c[2] + s)];
            let mut ch = vec![0u8; tape.len()];
            let r = eval_interval(&tape, &ctx, b, Some(&mut ch));
            let simple = tape.simplify(&ch);
            let pts: Vec<[f64; 3]> = (0..64).map(|_| [rng.range(b[0].lo, b[0].hi), rng.range(b[1].lo, b[1].hi), rng.range(b[2].lo, b[2].hi)]).collect();
            let (mut full, mut simp) = (Vec::new(), Vec::new());
            eval_points(&tape, &ctx, &pts, &mut full);
            eval_points(&simple, &ctx, &pts, &mut simp);
            for ((p, f), g) in pts.iter().zip(&full).zip(&simp) {
                if f.is_nan() {
                    continue;
                }
                let tol = 1e-9 * (1.0 + f.abs());
                assert!(r.lo - tol <= *f && *f <= r.hi + tol, "{src}: f({p:?}) = {f} outside [{}, {}]", r.lo, r.hi);
                assert!(f == g || (f - g).abs() <= tol, "{src}: simplified tape differs at {p:?}: {f} vs {g}");
            }
        }
    }
}

#[test]
fn span_bounds_are_sound_and_first_order() {
    // Over a time span, interval bounds are also first order in t (the mean
    // value form), which needs each value continuous in t: steps (floor,
    // mod), jumps (a select or division in t, atan2 across its cut) and
    // infinite slopes must fall back to plain bounds. Every value at a point
    // of the box, at any moment of the span (its ends included), must lie in
    // the bound over the span, and the tape simplified by the branches
    // decided over the span must agree there. So must the leaf centred form,
    // first order in t or not (the search skips the first where it failed
    // above). Spans of zero length, and a ctx.t outside the span (bounds
    // over a span ignore it), included.
    let mut fields: Vec<String> = ZOO.iter().map(|s| format!("(move (rotate-z (rotate-x {s} (* 50 t)) (* 90 t)) (* 3 (sin (* 4 t))) 0 (* 2 t))")).collect();
    fields.extend(
        [
            "(+ (sphere 5) (floor (* 7 t)))",
            "(select (- t 0.37) (sphere 8) (move (box 6 6 6) 4 0 0))",
            "(/ (sphere 5) (- t 0.5))",
            "(atan2 (- t 0.5) (- (* 0.1 x) 1))",
            "(rotate-z (extrude-z (involute-gear-2d 18 1.5 20 0.05) -5 5) (* 20 t))",
            "(* (sqrt (abs (- t 0.5))) (sphere 4))",
            "(max (sphere 4) (- (* 3 t) x) (abs (- y (* 10 t))))",
            "(pow (- t 0.3) (+ 2 (* 0.1 z)))",
            "(+ (tan (* 3 t)) (ln (abs (- t 0.6))) (exp t) (atan (* t y)))",
            // Seams turned or moved through the box (jumps as choices).
            "(rotate-z (repeat-polar (move (box 4 2 2) 10 0 0) 12) (* 90 t))",
            "(rotate-z (repeat-polar (move (repeat-polar (move (box 2 1 1) 4 0 0) 6) 9 0 0) 8) (* 45 t))",
            "(rotate-z (extrude-z (gt-pulley-teeth-2d 20 2 0) -3 3) (* 360 t))",
            "(move (repeat-x (sphere 2) 10) (* 13 t) 0 0)",
            "(+ (atan2 (- y (* 3 t)) (- x 2)) (* 0.1 (floor (+ x (* 5 t)))))",
            // Affine forms: values read from one source twice, through
            // products, quotients, curvature, maxes and mins of the two.
            "(max (move (cylinder-z 3 10) (* 9 t) 0 0) (neg (move (cylinder-z 3.1 10) (* 9 t) 0 0)))",
            "(min (rotate-z (move (sphere 2) 4 0 0) (* 200 t)) (neg (rotate-z (move (sphere 2.2) 4 0 0) (* 200 t))))",
            "(- (* (- x (* 3 t)) (+ y (sin (* 5 t)))) (/ (+ z 1) (+ 3 (cos (* 2 t)))))",
            "(max (sqrt (+ (square (- x (* 4 (cos (* 6 t))))) (square (- y (* 4 (sin (* 6 t))))))) (- (abs z) (* t t)))",
            "(emax (- x (* 2 t)) (- y t) (- (* 0.5 z) (exp (* 0.3 t))))",
        ]
        .map(String::from),
    );
    let mut rng = Lcg(23);
    for src in &fields {
        let (_d, tape) = tape_of(src);
        for _ in 0..150 {
            let c = [rng.range(-12.0, 12.0), rng.range(-12.0, 12.0), rng.range(-12.0, 12.0)];
            let s = rng.range(0.01, 4.0);
            let b = c.map(|v| I::new(v - s, v + s));
            let t0 = rng.range(-0.2, 1.0);
            let len = match rng.next() {
                x if x < 0.1 => 0.0,
                x if x < 0.3 => 1.0,
                _ => rng.range(1e-3, 0.3),
            };
            let span = (t0, t0 + len);
            let mid = 0.5 * (span.0 + span.1);
            let ctx = Ctx { t_span: Some(span), ..Ctx::new(if rng.next() < 0.2 { span.1 + 0.4 } else { mid }, vec![]) };
            let mut ch = vec![0u8; tape.len()];
            let r = eval_interval_scope(&tape, &ctx, b, Some(&mut ch));
            let (rt, rate) = eval_interval_rates(&tape, &ctx, b, &Inputs { t: true, params: Vec::new() });
            let lower = [centred_lower(&tape, &ctx, b).unwrap_or(f64::NEG_INFINITY)];
            let mut last = [0.0; 3];
            let simple = tape.simplify(&ch);
            // A few points, each swept densely over the span (its ends and
            // middle included).
            for k in 0..96 {
                let p = if k % 24 == 0 { b.map(|x| rng.range(x.lo, x.hi)) } else { [0.0; 3] };
                let p = if k % 24 == 0 {
                    last = p;
                    p
                } else {
                    last
                };
                let t = match k % 24 {
                    0 => span.0,
                    1 => span.1,
                    2 => mid,
                    j => span.0 + (span.1 - span.0) * (j - 2) as f64 / 21.0,
                };
                let at = Ctx::new(t, vec![]);
                let (f, g) = (eval_point(&tape, &at, p), eval_point(&simple, &at, p));
                if f.is_nan() {
                    continue;
                }
                let tol = 1e-9 * (1.0 + f.abs());
                assert!(r.lo - tol <= f && f <= r.hi + tol, "{src}: f({p:?}, t {t}) = {f} outside [{}, {}] over {span:?}", r.lo, r.hi);
                assert!(rt.lo - tol <= f && f <= rt.hi + tol, "{src}: f({p:?}, t {t}) = {f} outside first order [{}, {}] over {span:?}", rt.lo, rt.hi);
                assert!(lower.iter().all(|&lo| lo <= f + tol), "{src}: centred lower bounds {lower:?} above f({p:?}, t {t}) = {f} over {span:?}");
                assert!(f == g || (f - g).abs() <= tol, "{src}: simplified tape differs at {p:?}, t {t}: {f} vs {g}");
                // The rate there within the rates' bound.
                let dt = eval_dual(&tape, &at, p)[4];
                let rtol = 1e-9 * (1.0 + dt.abs());
                assert!(!dt.is_finite() || rate[0].lo - rtol <= dt && dt <= rate[0].hi + rtol, "{src}: d/dt at {p:?}, t {t} = {dt} outside rates {:?} over {span:?}", rate[0]);
            }
        }
    }
    // First order: a hole turning about its own axis. Its field does not
    // change with t, but plain bounds blur its turned coordinates along the
    // arc they sweep; the mean value form blurs it only by its rate's
    // spread over the cell (about the cell's width, not the radius, times
    // the turning speed). A cell inside the hole at r 3..3.07 (field -0.1..-0.03): over
    // 1/100 of a turn plain bounds reach above 0, first order ones do not.
    let (_d, tape) = tape_of("(rotate-z (cylinder-z 3.1 100) (* 360 t))");
    let b = [I::new(3.0, 3.05), I::new(0.0, 0.05), I::new(0.0, 0.05)];
    let ctx = Ctx { t_span: Some((0.3, 0.31)), ..Ctx::new(0.305, vec![]) };
    let (r, plain) = (eval_interval_rates(&tape, &ctx, b, &Inputs { t: true, params: Vec::new() }).0, eval_interval(&tape, &ctx, b, None));
    assert!(plain.hi > 0.0 && r.hi < 0.0 && r.hi - r.lo < 0.1, "[{}, {}], plain [{}, {}]", r.lo, r.hi, plain.lo, plain.hi);
}

#[test]
fn scope_bounds_are_sound_and_first_order_in_params() {
    // Over ranged params p and q (and over a time span too, now and then),
    // bounds are also first order in them: each value's
    // (eval_interval_scope), and the centred forms', whose pairs take each
    // term's value at the params' middles and its rates along them. That
    // needs each value continuous in them: steps, jumps (a select, a
    // division by something that may be 0, atan2 across its cut, a sampled
    // field's slope across a cell face), infinite slopes (sqrt and ln at 0)
    // must fall back to plain bounds, where a point derivative (0 away
    // from the jump) would not; and a promise bounds slopes in space, not
    // rates (x - 1000 p, promised 1-Lipschitz, moves 1000 per unit of p).
    // Every value at a point of the box, at param values of their spans
    // (their ends included) and a moment of the span, must lie within the
    // bounds, and the tape simplified by the branches decided over them
    // must agree there.
    let mut fields: Vec<String> = ZOO.iter().map(|s| format!("(move (rotate-z (rotate-x {s} (* 50 p)) (* 90 (+ q t))) (* 3 (sin (* 4 p))) q (* 2 p q))")).collect();
    fields.extend(
        [
            "(+ (sphere 5) (floor (* 7 p)))",
            "(select (- p 0.37) (sphere 8) (move (box 6 6 6) 4 0 0))",
            "(move (box 2 2 2) (select (- p 0.37) 0 3) 0 0)",
            "(/ (sphere 5) (- p 0.5))",
            "(* (- p 0.5) (/ 1 (- p 0.5)))",
            "(atan2 (- p 0.5) (- (* 0.1 x) 1))",
            "(* (sqrt (abs (- p 0.5))) (sphere 4))",
            "(pow (- p 0.3) (+ 2 (* 0.1 z)))",
            "(+ (tan (* 3 p)) (ln (abs (- q 0.6))) (exp p) (atan (* p y)))",
            "(max (lipschitz 1 (- x (* 1000 p))) (move (box 1 1 1) 5.5 0 0))",
            "(max (move (box (+ 10 p) 10 10) (/ (+ 10 p) 2) 0 0) (move (box (+ 10 q) 10 10) (+ 10 p (/ (+ 10 q) 2)) 0 0))",
            "(max (cylinder-z (+ 4 p) 10) (neg (cylinder-z (+ 4.05 p) 12)))",
            "(max (move (sphere 3) (* 10 p) 0 0) (neg (move (sphere 3.5) (+ (* 10 p) (* 0.1 q)) 0 0)))",
            "(move (grid-bake (torus-z 8 2) (res 0.5)) (* 4 p) q 0)",
            // Seams turned or moved through the box (jumps as choices),
            // and one the params do not move (the turn is about the
            // arms' own centre).
            "(rotate-z (repeat-polar (move (box 4 2 2) 10 0 0) 12) (* 90 (+ p q)))",
            "(rotate-z (repeat-polar (move (repeat-polar (move (box 2 1 1) 4 0 0) 6) 9 0 0) 8) (* 45 p))",
            "(rotate-z (extrude-z (gt-pulley-teeth-2d 20 2 0) -3 3) (* 360 p))",
            "(move (repeat-x (sphere 2) 10) (* 13 p) q 0)",
            "(repeat-polar (move (rotate-x (box 4 2 2) (* 40 p)) 10 (* 3 q) 0) 12)",
            // A slope unbounded along one input (sqrt at 0, moved by q and
            // t) and none along another (p): p's rate stays bounded.
            "(- (sqrt (+ (square (- x (* 3 q) t)) (square y))) (* 2 p))",
            "(* (sqrt (abs (- q 0.5))) (+ (sphere 4) p))",
        ]
        .map(String::from),
    );
    let mut rng = Lcg(29);
    for src in &fields {
        let d = doc(&format!("(param p 0) (param q 0) (def f {src})"));
        let ink::lang::Value::Field(f, _) = d.def("f").unwrap().clone() else { panic!() };
        let tape = Tape::new(&d.g, f);
        let (kp, kq) = (d.param_index("p").unwrap() as u32, d.param_index("q").unwrap() as u32);
        for _ in 0..150 {
            let c = [rng.range(-12.0, 12.0), rng.range(-12.0, 12.0), rng.range(-12.0, 12.0)];
            let s = rng.range(0.01, 4.0);
            let b = c.map(|v| I::new(v - s, v + s));
            let span = |rng: &mut Lcg| {
                let a = rng.range(-0.2, 1.0);
                let len = match rng.next() {
                    x if x < 0.1 => 0.0,
                    x if x < 0.3 => 1.0,
                    _ => rng.range(1e-3, 0.3),
                };
                (a, a + len)
            };
            let (sp, sq) = (span(&mut rng), span(&mut rng));
            let st = (rng.next() < 0.3).then(|| span(&mut rng));
            let mut params = vec![0.0; d.params.len()];
            params[kp as usize] = span_mid(sp);
            params[kq as usize] = span_mid(sq);
            let ctx = Ctx { t_span: st, p_spans: ParamSpans::new([(kp, sp), (kq, sq)]), ..Ctx::new(st.map_or(0.0, span_mid), params) };
            let mut ch = vec![0u8; tape.len()];
            let r = eval_interval_scope(&tape, &ctx, b, Some(&mut ch));
            let inputs = Inputs { t: st.is_some(), ..Inputs::of([&tape], &ctx) };
            let (rt, rates) = eval_interval_rates(&tape, &ctx, b, &inputs);
            let lower = [centred_lower(&tape, &ctx, b).unwrap_or(f64::NEG_INFINITY)];
            let simple = tape.simplify(&ch);
            let at_end = |rng: &mut Lcg, s: (f64, f64), k: usize| match k % 4 {
                0 => s.0,
                1 => s.1,
                2 => span_mid(s),
                _ => rng.range(s.0, s.1),
            };
            for k in 0..96 {
                let pt = b.map(|x| rng.range(x.lo, x.hi));
                let t = st.map_or(0.0, |s| at_end(&mut rng, s, k / 4));
                let mut vals = ctx.params.clone();
                vals[kp as usize] = at_end(&mut rng, sp, k);
                vals[kq as usize] = at_end(&mut rng, sq, k / 2);
                let at = Ctx::new(t, vals);
                let (v, g) = (eval_point(&tape, &at, pt), eval_point(&simple, &at, pt));
                if v.is_nan() {
                    continue;
                }
                let tol = 1e-9 * (1.0 + v.abs());
                assert!(r.lo - tol <= v && v <= r.hi + tol, "{src}: f({pt:?}, {:?}, t {t}) = {v} outside [{}, {}] over {sp:?} {sq:?} {st:?}", at.params, r.lo, r.hi);
                assert!(
                    rt.lo - tol <= v && v <= rt.hi + tol,
                    "{src}: f({pt:?}, {:?}, t {t}) = {v} outside first order in t [{}, {}] over {sp:?} {sq:?} {st:?}",
                    at.params,
                    rt.lo,
                    rt.hi
                );
                assert!(lower.iter().all(|&lo| lo <= v + tol), "{src}: centred lower bounds {lower:?} above f({pt:?}, {:?}) = {v} over {sp:?} {sq:?} {st:?}", at.params);
                assert!(v == g || (v - g).abs() <= tol, "{src}: simplified tape differs at {pt:?}, {:?}: {v} vs {g}", at.params);
                // The rates there within the rates' bounds.
                let dual = eval_dual(&tape, &at, pt);
                let along = inputs.t.then_some(4).into_iter().chain(inputs.params.iter().map(|&k| 5 + k as usize));
                for (j, x) in along.enumerate() {
                    let (dx, rtol) = (dual[x], 1e-9 * (1.0 + dual[x].abs()));
                    assert!(
                        !dx.is_finite() || rates[j].lo - rtol <= dx && dx <= rates[j].hi + rtol,
                        "{src}: rate {j} at {pt:?}, {:?}, t {t} = {dx} outside {:?}",
                        at.params,
                        rates[j]
                    );
                }
            }
        }
    }
}

#[test]
fn a_rate_along_an_input_nothing_moves_along_is_zero_beside_an_unbounded_one() {
    // A bore of diameter p turning by q: on the axis its sqrt has no slope
    // bound, so its rate along q is unbounded, but nothing under it reads
    // p: the bore's rate along p is its diameter's, -1/2, not unbounded
    // (which would send a search to halve p's band to its finest before
    // space, and leave a pin turning in a toleranced bore undecided).
    let d = doc("(param p 8) (param q 0) (def f (- (sqrt (+ (square (- (* (cos q) y) (* (sin q) z))) (square (+ (* (sin q) y) (* (cos q) z))))) (/ p 2)))");
    let ink::lang::Value::Field(f, _) = d.def("f").unwrap().clone() else { panic!() };
    let tape = Tape::new(&d.g, f);
    let (kp, kq) = (d.param_index("p").unwrap() as u32, d.param_index("q").unwrap() as u32);
    let mut params = vec![0.0; d.params.len()];
    params[kp as usize] = 8.018;
    let ctx = Ctx { p_spans: ParamSpans::new([(kp, (8.0, 8.036)), (kq, (0.0, 1.0))]), ..Ctx::new(0.0, params) };
    let b = [I::new(-1.0, 1.0); 3];
    let (_, d) = eval_interval_rates(&tape, &ctx, b, &Inputs { t: false, params: vec![kp, kq] });
    assert!(d[0].contains(-0.5) && d[0].hi - d[0].lo < 1e-12, "along p: {d:?}");
    assert!(!(d[1].lo.is_finite() && d[1].hi.is_finite()), "along q: {d:?}");
}

#[test]
fn inputs_past_the_first_eight_stay_ranged() {
    // Ten ranged params and a time span: first-order bounds take eight of
    // the params (MAX_INPUTS). The other two and time must stay ranged in
    // each value's bound at the inputs' centre: fixed at their middles, the
    // bound would miss the values at their ends by up to 3. At every corner
    // of the eleven spans the value lies within both bounds.
    let names: Vec<String> = (0..10).map(|k| format!("p{k}")).collect();
    let params: String = names.iter().map(|n| format!("(param {n} 0)")).collect();
    let sum = format!("(+ t {})", names.join(" "));
    for src in [sum.clone(), format!("(max (+ x {sum}) (- 1 x))")] {
        let d = doc(&format!("{params} (def f {src})"));
        let ink::lang::Value::Field(f, _) = d.def("f").unwrap().clone() else { panic!() };
        let tape = Tape::new(&d.g, f);
        let ks: Vec<usize> = names.iter().map(|n| d.param_index(n).unwrap()).collect();
        let ctx = Ctx { t_span: Some((-1.0, 1.0)), p_spans: ParamSpans::new(ks.iter().map(|&k| (k as u32, (-1.0, 1.0)))), ..Ctx::new(0.0, vec![0.0; d.params.len()]) };
        assert_eq!(Inputs::of([&tape], &ctx).len(), MAX_INPUTS);
        let b = [I::new(-1.0, 1.0); 3];
        let r = eval_interval_scope(&tape, &ctx, b, None);
        let lower = centred_lower(&tape, &ctx, b).unwrap();
        for bits in 0..1u32 << 11 {
            let end = |j: u32| if bits & (1 << j) == 0 { -1.0 } else { 1.0 };
            let mut vals = ctx.params.clone();
            for (j, &k) in ks.iter().enumerate() {
                vals[k] = end(j as u32 + 1);
            }
            let v = eval_point(&tape, &Ctx::new(end(0), vals.clone()), [0.25; 3]);
            assert!(r.lo <= v && v <= r.hi && lower <= v, "{src}: f at {vals:?}, t {} = {v}, bounds {r:?}, centred {lower}", end(0));
        }
    }
}

#[test]
fn faces_one_param_moves_together_are_bounded_tight() {
    // Block a, (box a 10 10) placed at a / 2, and block b
    // placed at a + 5, touch at x = a for every a in 9.9..10.1. Plain
    // bounds read a twice in each and lose the dependency: a's fixed left
    // face is blurred by the band, and the faces that touch become a slab
    // of possible contact 0.2 thick. First order in a, the left face moves
    // by its rate (0), and the pair of touching faces by their combined
    // rate (1/2 (-1) + 1/2 (1) = 0), so the touching cells are bounded
    // within rounding of 0 wherever the contact lies in the band.
    let d = doc("(param a 10) (def blocks (max (move (box a 10 10) (/ a 2) 0 0) (move (box 10 10 10) (+ a 5) 0 0)))");
    let ink::lang::Value::Field(f, _) = d.def("blocks").unwrap().clone() else { panic!() };
    let tape = Tape::new(&d.g, f);
    let k = d.param_index("a").unwrap() as u32;
    let ctx = Ctx { p_spans: ParamSpans::new([(k, (9.9, 10.1))]), ..Ctx::new(0.0, vec![10.0]) };
    for x0 in [9.88, 9.95, 9.99, 10.05] {
        let b = [I::new(x0, x0 + 0.04), I::new(1.0, 1.04), I::new(-2.0, -1.96)];
        let plain = eval_interval(&tape, &ctx, b, None);
        let lo = centred_lower(&tape, &ctx, b).unwrap();
        assert!(plain.lo < -0.01 && lo > -1e-12 && lo <= 0.0, "x {x0}: plain {plain:?}, centred {lo}");
    }
    // Block a alone, near its left face: plain bounds 0.1 off, first order
    // within rounding of the cell's own spread.
    let d = doc("(param a 10) (def block (move (box a 10 10) (/ a 2) 0 0))");
    let ink::lang::Value::Field(f, _) = d.def("block").unwrap().clone() else { panic!() };
    let tape = Tape::new(&d.g, f);
    let b = [I::new(-0.06, -0.02), I::new(1.0, 1.04), I::new(-2.0, -1.96)];
    let (plain, r) = (eval_interval(&tape, &ctx, b, None), eval_interval_scope(&tape, &ctx, b, None));
    assert!(plain.lo < 0.0 && (r.lo - 0.02).abs() < 1e-12 && (r.hi - 0.06).abs() < 1e-12, "plain {plain:?}, first order {r:?}");
}

#[test]
fn a_grid_face_a_hair_inside_a_span_is_a_jump() {
    // A sampled field's x-slope is 0 in the cell x < 0 and 1 in x >= 0.
    // Read at x = t over t in [-2^-54, 1/4]: the low end's local coordinate
    // 1 - 2^-54 rounds onto the face, yet the slope there is 0. The span
    // crosses the jump, so its rate is unbounded and its bound holds both.
    use ink::graph::{Graph, Node, Var};
    use ink::grid::GridData;
    let layer = [0.0f32, 0.0, 1.0];
    let data = GridData::new([-1.0, 0.0, 0.0], 1.0, [3, 2, 2], layer.repeat(4));
    let mut g = Graph::new();
    let (t, c) = (g.var(Var::T), g.constant(0.5));
    g.grids.push(Arc::new(data));
    g.nodes.push(Node::Grid(0, 1, [t, c, c]));
    let tape = Tape::new(&g, (g.nodes.len() - 1) as u32);
    let span = (-(2f64.powi(-54)), 0.25);
    let ctx = Ctx { t_span: Some(span), ..Ctx::new(0.125, vec![]) };
    let r = eval_interval_rates(&tape, &ctx, [I::point(0.0); 3], &Inputs { t: true, params: Vec::new() }).0;
    assert!(r.contains(0.0) && r.contains(1.0), "[{}, {}]", r.lo, r.hi);
}

#[test]
fn gradient_is_nonzero_on_box_faces() {
    // Regression: min(d, 0) ties exactly on the surface must not zero the gradient.
    let (_d, tape) = tape_of("(box 10 10 10)");
    let ctx = Ctx::new(0.0, vec![]);
    let g = eval_grads(&tape, &ctx, &[[5.0, 1.0, 2.0], [0.0, -5.0, 0.0]]);
    assert!(g[0][1] > 0.4 && g[1][2] < -0.4, "{g:?}");
}

#[test]
fn centred_lower_bounds_hold_on_tiny_boxes() {
    // Truthful promises over boxes where a naive centre, half-width or
    // slope misses part of the box: one float wide, where the centre rounds
    // to one end and the far end is a whole width away; so small that
    // squared half-widths underflow to 0; and huge, under a warp whose
    // contraction takes the promised slope below the smallest float. The
    // bound must not pass the field's value at the point named.
    let ctx = Ctx::new(0.0, vec![]);
    let u = f64::EPSILON;
    let (o, w) = (I::point(0.0), I::new(-1e230, 3e230));
    let cases = [
        ("(lipschitz 1 (- (select (- x 1.0000000000000002) (- 1.0000000000000002 x) (- x 1.0000000000000002)) 5.551115123125783e-17))", [I::new(1.0, 1.0 + u), o, o], 1.0 + u),
        ("(lipschitz 1 (- (select x (neg x) x) 5e-201))", [I::new(-1e-200, 3e-200), o, o], 0.0),
        // A truthful tiny promise under a contracting warp: its slope
        // (1e-200 x 1e-130) underflows, and rounded to nearest read as 0.
        ("(warp (lipschitz 1e-200 (- (* 1e-200 (select x (neg x) x)) 2.5e-101)) (* 1e-130 x) (* 1e-130 y) (* 1e-130 z))", [w; 3], 0.0),
    ];
    for (src, b, at) in cases {
        let (_d, tape) = tape_of(src);
        let f = eval_point(&tape, &ctx, [at, 0.0, 0.0]);
        assert!(f < 0.0, "{src}: f = {f:e}");
        let lo = centred_lower(&tape, &ctx, b).unwrap();
        assert!(lo <= f, "{src}: lower bound {lo:e} above f({at:e}) = {f:e}");
    }
}

#[test]
fn centred_lower_rounds_its_radius_up_below_the_normal_range() {
    // F = 9/8 Q - 5/8 K (x + y) with K = 2^1023, Q = K u, over the box within
    // hw = [u, u, 0] of 0 (u the smallest subnormal): |grad F| < K, and at
    // the corner (u, u, 0) F = -Q/8. The scalar spread needs |hw| = sqrt(2)
    // u rounded up (2 u), not to nearest (u, which leaves the bound at +Q/8).
    use ink::centred::term_lower;
    use ink::interval::{G, NO_GRAM};
    let (k, u) = (2f64.powi(1023), f64::from_bits(1));
    let q = k * u;
    let d = [I::point(-0.625 * k), I::point(-0.625 * k), I::point(0.0)];
    let a = G { v: I::ENTIRE, d, l: k, s2: None, m: NO_GRAM };
    let ac = G { v: I::point(1.125 * q), ..a };
    let lo = term_lower(ac, a, [u, u, 0.0]);
    assert!(lo <= -0.125 * q, "lower bound {lo:e} above F(u, u, 0) = {:e}", -0.125 * q);
}

fn clipped(d: &mut Doc, target: &str) -> (Arc<Tape>, Ctx) {
    let ctx = d.ctx(&[], 0.0).unwrap();
    let it = d.items(target).unwrap().remove(0);
    let f = d.placed_field(&it);
    let f = d.clip_reach(f, &ctx).0;
    (Arc::new(Tape::new(&d.g, f)), ctx)
}

#[test]
fn sharp_features_are_preserved() {
    let mut d = doc("(body b (shape (rotate-z (rotate-x (box 20 30 10) 30) 20)))");
    let (tape, ctx) = clipped(&mut d, "b");
    let h = 0.5;
    let m = mesh::mesh(&tape, &ctx, [-20.0; 3], [20.0; 3], h);
    let mut fv = Vec::new();
    eval_points(&tape, &ctx, &m.verts, &mut fv);
    let max = fv.iter().fold(0.0f64, |a, v| a.max(v.abs()));
    assert!(max < 0.05 * h, "vertices drift {max} off the surface");
    let a = mesh::analyze(&m);
    assert_eq!(a.boundary_edges, 0);
    assert!((a.volume - 6000.0).abs() < 1.0, "volume {}", a.volume);
    // One vertex per cell clips each box corner by a fraction of a cell, so
    // area runs ~0.1% low at h = 0.5 (it converges as h shrinks).
    assert!((a.area - 2200.0).abs() < 2.2, "area {}", a.area);
}

/// A body's uniform mesh at h, and the same simplified to tol against its
/// field (as `ink mesh` does).
fn simplified(src: &str, h: f64, tol: f64, field: bool) -> (mesh::Mesh, ink::simplify::Simplified, Arc<Tape>, Ctx) {
    let mut d = doc(&format!("(body b (shape {src}))"));
    let (tape, ctx) = clipped(&mut d, "b");
    let (lo, hi) = d.domain();
    let (blo, bhi, _) = query::bounds(&tape, &ctx, lo, hi).unwrap();
    let fine = mesh::mesh(&tape, &ctx, blo, bhi, h);
    let local = ink::octree::Local::new(&tape, &ctx, blo, 8.0 * h);
    let f = |p: &[[f64; 3]]| local.eval(p);
    let s = ink::simplify::simplify(&fine, tol, if field { Some(&f) } else { None }, None, &[]);
    (fine, s, tape, ctx)
}

/// The sampled distance from a mesh to the field's surface, as `ink mesh`
/// measures it.
fn deviation(m: &mesh::Mesh, tape: &Arc<Tape>, ctx: &Ctx, h: f64, step: f64) -> ink::simplify::Deviation {
    let local = ink::octree::Local::new(tape, ctx, [-1e3; 3], 8.0 * h);
    ink::simplify::deviation(&|p: &[[f64; 3]]| local.eval(p), &|p: &[[f64; 3]]| local.grads(p), m, h, step, 2.0 * h)
}

/// Distance from p to triangle abc, by a different route than the
/// simplifier's: the nearest of the plane's foot (if inside) and the edges.
fn point_triangle(p: [f64; 3], t: [[f64; 3]; 3]) -> f64 {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let seg = |a: [f64; 3], b: [f64; 3]| {
        let (ab, ap) = (sub(b, a), sub(p, a));
        let s = (dot(ap, ab) / dot(ab, ab)).clamp(0.0, 1.0);
        let q = sub(ap, ab.map(|c| c * s));
        dot(q, q).sqrt()
    };
    let mut d = seg(t[0], t[1]).min(seg(t[1], t[2])).min(seg(t[2], t[0]));
    let (u, v) = (sub(t[1], t[0]), sub(t[2], t[0]));
    let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    let nn = dot(n, n).sqrt();
    if nn > 0.0 {
        let h = dot(sub(p, t[0]), n) / nn;
        let q = sub(p, n.map(|c| c * h / nn));
        let side = |a: [f64; 3], b: [f64; 3]| {
            let (e, w) = (sub(b, a), sub(q, a));
            dot([e[1] * w[2] - e[2] * w[1], e[2] * w[0] - e[0] * w[2], e[0] * w[1] - e[1] * w[0]], n) >= 0.0
        };
        if side(t[0], t[1]) && side(t[1], t[2]) && side(t[2], t[0]) {
            d = d.min(h.abs());
        }
    }
    d
}

/// The largest distance from the given points to the mesh, by brute force.
fn farthest_from(m: &mesh::Mesh, pts: &[[f64; 3]]) -> f64 {
    let nearest = |p: [f64; 3]| m.tris.iter().map(|t| point_triangle(p, t.map(|i| m.verts[i as usize]))).fold(f64::INFINITY, f64::min);
    let parts = std::thread::available_parallelism().map_or(4, |n| n.get());
    std::thread::scope(|s| {
        let jobs: Vec<_> = pts.chunks(pts.len().div_ceil(parts).max(1)).map(|c| s.spawn(move || c.iter().map(|&p| nearest(p)).fold(0.0, f64::max))).collect();
        jobs.into_iter().map(|j| j.join().unwrap()).fold(0.0, f64::max)
    })
}

#[test]
fn simplified_mesh_stays_closed_and_within_tolerance() {
    // Flat faces, a curved bore, and sharp edges and corners.
    let (h, tol) = (0.25, 0.01);
    let (fine, s, tape, ctx) = simplified("(subtract (box 20 30 10) (cylinder-z 5 20))", h, tol, true);
    let m = &s.mesh;
    let a = mesh::analyze(m);
    assert_eq!((a.boundary_edges, a.nonmanifold_edges, a.shells.len()), (0, 0, 1), "closed 2-manifold, one shell");
    assert!(m.tris.len() * 100 < fine.tris.len(), "{} of {} triangles", m.tris.len(), fine.tris.len());
    // Every vertex at res within tol of the result: tracked, and checked here
    // by brute force.
    assert!(s.moved.iter().all(|&d| d <= tol));
    let far = farthest_from(m, &fine.verts);
    assert!(far <= tol * (1.0 + 1e-9), "a vertex at res lies {far} from the simplified mesh");
    // And the result near the surface: no farther than the uniform mesh is,
    // plus tol.
    let step = tol / 4.0;
    let (dev, dev0) = (deviation(m, &tape, &ctx, h, step), deviation(&fine, &tape, &ctx, h, step));
    assert_eq!(dev.beyond, 0);
    assert!(dev.max <= dev0.max + tol, "deviation {} (uniform {})", dev.max, dev0.max);
    // As crisp as the uniform mesh: its vertex at each corner of the box is
    // kept.
    let nearest =
        |m: &mesh::Mesh, c: [f64; 3]| m.verts.iter().map(|v| ((v[0] - c[0]).powi(2) + (v[1] - c[1]).powi(2) + (v[2] - c[2]).powi(2)).sqrt()).fold(f64::INFINITY, f64::min);
    for k in 0..8 {
        let corner = [10.0, 15.0, 5.0];
        let corner = std::array::from_fn(|a| if k >> a & 1 == 0 { -corner[a] } else { corner[a] });
        let (now, was) = (nearest(m, corner), nearest(&fine, corner));
        assert!(now <= was && was < 2.0 * tol, "corner {corner:?}: nearest vertex {now} away, {was} at res");
    }
    let exact = 6000.0 - PI * 25.0 * 10.0;
    assert!((a.volume - exact).abs() / exact < 2e-3, "volume {} vs {exact}", a.volume);
}

#[test]
fn simplified_thin_walls_stay_apart() {
    // A hollow ball whose wall (0.4) is under two cells thick: two shells,
    // each closed, neither reaching through the wall to the other.
    let (h, tol) = (0.25, 0.01);
    let (fine, s, _, _) = simplified("(shell (sphere 10) 0.4)", h, tol, true);
    let a = mesh::analyze(&s.mesh);
    assert_eq!((a.boundary_edges, a.nonmanifold_edges, a.shells.len()), (0, 0, 2));
    assert!(s.mesh.tris.len() * 3 < fine.tris.len(), "{} of {} triangles", s.mesh.tris.len(), fine.tris.len());
    let exact = 4.0 / 3.0 * PI * (10.2f64.powi(3) - 9.8f64.powi(3));
    assert!((a.volume - exact).abs() / exact < 5e-3, "volume {} vs {exact}", a.volume);
    let sample: Vec<[f64; 3]> = fine.verts.iter().step_by(7).copied().collect();
    assert!(farthest_from(&s.mesh, &sample) <= tol * (1.0 + 1e-9));
}

#[test]
fn the_field_keeps_simplified_triangles_on_the_surface() {
    // Vertices at res alone (spaced 1 apart) leave room between them: a long
    // edge's midpoint can sag farther than tol from a sphere with every
    // vertex within tol. The field's check at new edges and centroids keeps
    // the result itself within tol (and the uniform mesh's own error).
    let (h, tol) = (1.0, 0.04);
    let (fine, with, tape, ctx) = simplified("(sphere 20)", h, tol, true);
    let (_, without, _, _) = simplified("(sphere 20)", h, tol, false);
    let step = tol / 16.0;
    let dev0 = deviation(&fine, &tape, &ctx, h, step).max;
    let (a, b) = (deviation(&with.mesh, &tape, &ctx, h, step).max, deviation(&without.mesh, &tape, &ctx, h, step).max);
    assert!(a <= tol + dev0 + step, "with the field: {a} (uniform {dev0})");
    assert!(b > tol + dev0, "without the field: {b}, with: {a}");
}

#[test]
fn mesh_reports_features_thinner_than_res() {
    // A plate 0.4 thick: at res 1 most of it falls between samples, and the
    // cells it passes through are counted, with one centre inside it.
    // (The grid's planes, z = -0.2 + k h, miss the plate's faces.)
    let mut d = doc("(body b (shape (move (box 20 20 0.4) 0 0 0.3)))");
    let (tape, ctx) = clipped(&mut d, "b");
    let (lo, hi) = ([-10.0, -10.0, -0.2], [10.0, 10.0, 0.6]);
    let (m, u) = mesh::mesh_checked(&tape, &ctx, lo, hi, 1.0);
    assert!(u.cells > 0, "{} triangles, no gap reported", m.tris.len());
    let at = u.at.unwrap();
    assert!(eval_point(&tape, &ctx, at) < 0.0, "{at:?} is not in the plate");
    let (_, u) = mesh::mesh_checked(&tape, &ctx, lo, hi, 0.1);
    assert_eq!(u.cells, 0);
}

#[test]
fn local_tapes_agree_with_the_whole() {
    // Points scattered over blocks of 2.5: each block's simplified tape gives
    // what the whole tape gives, values and gradients.
    let ctx = Ctx::new(0.0, vec![]);
    let mut rng = Lcg(11);
    for src in ZOO {
        let (_d, tape) = tape_of(src);
        let tape = Arc::new(tape);
        let local = ink::octree::Local::new(&tape, &ctx, [-0.3, 0.1, -0.7], 2.5);
        let pts: Vec<[f64; 3]> = (0..4000).map(|_| [rng.range(-15.0, 15.0), rng.range(-15.0, 15.0), rng.range(-15.0, 15.0)]).collect();
        let mut whole = Vec::new();
        eval_points(&tape, &ctx, &pts, &mut whole);
        let near = local.eval(&pts);
        for (k, (a, b)) in whole.iter().zip(&near).enumerate() {
            assert!(a == b || (a.is_nan() && b.is_nan()), "{src} at {:?}: whole {a}, local {b}", pts[k]);
        }
        let (gw, gl) = (eval_grads(&tape, &ctx, &pts), local.grads(&pts));
        assert!(gw.iter().zip(&gl).all(|(a, b)| a.iter().zip(b).all(|(x, y)| x == y || (x - y).abs() <= 1e-12 * x.abs().max(1.0))), "{src}: gradients differ");
        // A lone point runs the batch's code in one lane: the same values.
        for (p, w) in pts.iter().zip(&whole).take(500) {
            let v = eval_point(&tape, &ctx, *p);
            assert!(v == *w || (v.is_nan() && w.is_nan()), "{src} at {p:?}: point {v}, batch {w}");
        }
    }
}

#[test]
fn a_box_gets_the_tape_of_the_smallest_block_holding_it() {
    // Boxes of every size from 1e-6 to 10, anywhere (on blocks' faces
    // too): the tape for each gives the whole tape's values in it, and its
    // bound holds them.
    let ctx = Ctx::new(0.0, vec![]);
    let mut rng = Lcg(5);
    for src in ZOO {
        let (_d, tape) = tape_of(src);
        let tape = Arc::new(tape);
        let local = ink::octree::Local::new(&tape, &ctx, [-0.3, 0.1, -0.7], 0.5);
        for i in 0..300 {
            let c: [f64; 3] = std::array::from_fn(|_| rng.range(-15.0, 15.0));
            // Every third box starts on a block's face.
            let c = if i % 3 == 0 { [-0.3 + 0.5 * (c[0] / 0.5).round(), c[1], c[2]] } else { c };
            let w = 10f64.powf(rng.range(-6.0, 1.0));
            let b: [I; 3] = std::array::from_fn(|k| I::new(c[k], c[k] + w * rng.next()));
            let t = local.tape_over(&b);
            let r = eval_interval(&t, &ctx, b, None);
            for _ in 0..8 {
                let p: [f64; 3] = std::array::from_fn(|k| b[k].lo + (b[k].hi - b[k].lo) * rng.next());
                let (v, whole) = (eval_point(&t, &ctx, p), eval_point(&tape, &ctx, p));
                assert!(v == whole && r.lo <= v && v <= r.hi, "{src} at {p:?} in {b:?}: {v} (whole {whole}) not in {r:?}");
            }
        }
    }
    // In a row of 400 spheres 3 apart, a box in one smallest block (x 600
    // to 600.5, y and z 0 to 0.5) gets one sphere's tape; one across
    // blocks' faces, the block of 4 holding it (x 598 to 602, touching
    // three).
    let row = (0..400).map(|k| format!("(move (sphere 1) {} 0 0)", 3 * k)).collect::<Vec<_>>().join(" ");
    let (_d, tape) = tape_of(&format!("(min {row})"));
    let tape = Arc::new(tape);
    let local = ink::octree::Local::new(&tape, &ctx, [-2.0; 3], 0.5);
    let one = tape_of("(move (sphere 1) 600 0 0)").1.len();
    let t = local.tape_over(&[I::new(600.1, 600.3), I::new(0.1, 0.3), I::new(0.1, 0.2)]);
    assert!(t.len() <= one, "{} (one sphere {one}, all {})", t.len(), tape.len());
    let t = local.tape_over(&[I::new(599.5, 600.2), I::new(0.1, 0.3), I::new(-0.2, 0.2)]);
    assert!(t.len() <= 3 * one + 4, "{} (one sphere {one}, all {})", t.len(), tape.len());
}

#[test]
fn mesh_deviation_bounds_its_samples() {
    // For a sphere the distance to the surface is known: | |p| - r |. The
    // reported largest bound covers it at every vertex (each a sample), and
    // is not far above the largest over the mesh.
    let mut d = doc("(body b (shape (sphere 10)))");
    let (tape, ctx) = clipped(&mut d, "b");
    let h = 1.0;
    let m = mesh::mesh(&tape, &ctx, [-10.0; 3], [10.0; 3], h);
    let step = 1e-3;
    let dev = deviation(&m, &tape, &ctx, h, step);
    let dist = |p: &[f64; 3]| ((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() - 10.0).abs();
    let at_vertices = m.verts.iter().map(dist).fold(0.0, f64::max);
    let mut rng = Lcg(7);
    let over_mesh = m
        .tris
        .iter()
        .flat_map(|t| {
            let p = t.map(|i| m.verts[i as usize]);
            let (a, b) = (rng.next(), rng.next());
            let (a, b) = if a + b > 1.0 { (1.0 - a, 1.0 - b) } else { (a, b) };
            [std::array::from_fn::<f64, 3, _>(|c| p[0][c] + a * (p[1][c] - p[0][c]) + b * (p[2][c] - p[0][c]))]
        })
        .map(|p| dist(&p))
        .fold(at_vertices, f64::max);
    assert!(dev.max >= at_vertices, "{} below the vertices' {at_vertices}", dev.max);
    assert!(dev.max <= 1.2 * over_mesh + step, "{} well above the mesh's {over_mesh}", dev.max);
}

#[test]
fn mass_integrates_graded_density() {
    // Density 1 + z/10 over a 10 mm cube centred at the origin: mass = volume.
    let mut d = doc("(body b (shape (box 10 10 10)) (attr density (+ 1 (/ z 10))))");
    let (tape, ctx) = clipped(&mut d, "b");
    let it = d.items("b").unwrap().remove(0);
    let ink::lang::Value::Field(g, _) = d.placed_attrs(&it)[0].1 else { panic!() };
    let gt = Tape::new(&d.g, g);
    let r = query::sums(&tape, &[Some(&gt)], true, &ctx, [-5.0; 3], [5.0; 3], 0.1).sums[0];
    assert!((r - 1000.0).abs() < 1.0, "{r}");
}

#[test]
fn bounds_are_proven_and_tight() {
    let mut d = doc("(body b (shape (move (box 10 20 30) 100 -50 7)))");
    let (tape, ctx) = clipped(&mut d, "b");
    let (lo, hi) = d.domain();
    let (blo, bhi, h) = query::bounds(&tape, &ctx, lo, hi).unwrap();
    let want = ([95.0, -60.0, -8.0], [105.0, -40.0, 22.0]);
    for k in 0..3 {
        assert!(blo[k] <= want.0[k] && bhi[k] >= want.1[k], "not an outer bound");
        assert!(want.0[k] - blo[k] <= 2.0 * h && bhi[k] - want.1[k] <= 2.0 * h, "not tight");
    }
}

#[test]
fn a_flat_field_inside_is_not_covered_by_its_slope() {
    // A box whose field is flat (slope 0) through its middle: the step up
    // the slope is NaN there, which no cover may take as a point outside
    // (which would prove this 0.22 mm cube empty to 0.05 mm, inside -0.001).
    let (d, t) = tape_of("(max (box 0.22 0.22 0.22) -0.001)");
    let ctx = Ctx::new(0.0, d.params.iter().map(|p| p.default).collect());
    let s = query::prove_empty(&Arc::new(t), &ctx, [-1.0; 3], [1.0; 3], None, 0.025, 0.05);
    assert!(matches!(s.verdict, Verdict::Disproven { .. }), "{:?}", s.verdict);
}

#[test]
fn emptiness_proofs() {
    let mut d = doc("(body a (shape (sphere 10)))
         (body b (shape (move (sphere 10) 25 0 0)))
         (body c (shape (move (sphere 10) 15 0 0)))
         (def ab (max (shape-of a) (shape-of b)))
         (def ac (max (shape-of a) (shape-of c)))");
    let (lo, hi) = d.domain();
    let (t, ctx) = clipped(&mut d, "ab");
    assert!(matches!(query::prove_empty(&t, &ctx, lo, hi, None, f64::INFINITY, 0.0).verdict, Verdict::Proven));
    let (t, ctx) = clipped(&mut d, "ac");
    let Verdict::Disproven { witness, value, .. } = query::prove_empty(&t, &ctx, lo, hi, None, f64::INFINITY, 0.0).verdict else { panic!("overlap not found") };
    assert!(value < 0.0);
    let r1 = (witness[0].powi(2) + witness[1].powi(2) + witness[2].powi(2)).sqrt();
    let r2 = ((witness[0] - 15.0).powi(2) + witness[1].powi(2) + witness[2].powi(2)).sqrt();
    assert!(r1 < 10.0 && r2 < 10.0, "witness {witness:?} is not in both spheres");
}

#[test]
fn named_constructions_are_placed_with_their_instance() {
    // (shape-of path name): the peg inside a placed part, in world space.
    let d = ink::load_str(
        "(import \"std\")
         (def peg (cylinder-z 1 4))
         (body part (shape (min (box 10 10 2) (move peg 3 0 3))))
         (assembly root (instance p part (at (rot-z 90) (move 20 0 0))))
         (def peg-world (shape-of root/p peg))",
        "t",
        Path::new(ROOT),
    )
    .unwrap();
    let ink::lang::Value::Field(f, _) = d.def("peg-world").unwrap().clone() else { panic!() };
    let t = Tape::new(&d.g, f);
    let ctx = Ctx::new(0.0, vec![]);
    // The peg's axis moves from (3, 0) to (20, 3) (a quarter turn, then 20 in x).
    assert!((eval_point(&t, &ctx, [20.0, 3.0, 3.0]) + 1.0).abs() < 1e-9);
    assert!(eval_point(&t, &ctx, [20.0, 0.0, 0.0]) > 1.0, "the plate is not part of the peg");
    let e = ink::load_str("(import \"std\") (body b (shape (sphere 1))) (def q (shape-of b nope))", "t", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("has no piece built as (tag nope ...)"), "{}", e.message);
}

#[test]
fn shape_of_is_kept_by_what_it_reads() {
    // (shape-of path) is kept between reads (an assembly's claims read it
    // again and again); a construction inside the same path is not the path.
    let d = ink::load_str(
        "(import \"std\")
         (def peg (cylinder-z 1 4))
         (body part (shape (min (box 10 10 2) (move peg 3 0 3))))
         (assembly root (instance p part (at (move 20 0 0))))
         (def plate-not-peg (max (shape-of root/p) (- (shape-of root/p peg))))",
        "t",
        Path::new(ROOT),
    )
    .unwrap();
    let ink::lang::Value::Field(f, _) = d.def("plate-not-peg").unwrap().clone() else { panic!() };
    let t = Tape::new(&d.g, f);
    let ctx = Ctx::new(0.0, vec![]);
    assert!(eval_point(&t, &ctx, [16.0, -4.0, 0.0]) < 0.0, "the plate away from the peg");
    assert!(eval_point(&t, &ctx, [23.0, 0.0, 3.0]) > 0.0, "the peg is not the plate's");
}

#[test]
fn namespaced_imports_keep_files_apart() {
    // The same library twice, under two names: every top-level name it
    // defines is prefixed, its own references follow, and local bindings
    // (a let's `size`) and syntax (a body's `(shape ...)`) are left alone.
    let dir = crate::common::Scratch::new("ns");
    std::fs::write(
        dir.join("lib.ink"),
        "(import \"std\") (param size 10) (fn peg (r) (cylinder-z r size))
         (def shape (subtract (box 40 40 4) (peg 3))) (body bracket (shape shape))
         (assembly root (instance bracket bracket (at (move 5 0 0))))
         (fn scaled (s) (let ((size (* s 2)) (k size)) (box k k k)))",
    )
    .unwrap();
    let d = ink::load_str(
        "(import \"std\") (import \"lib.ink\" as a) (import \"lib.ink\" as b)
         (def both (min (shape-of a.root) (shape-of b.bracket))) (def k (a.scaled 3))",
        "t",
        &dir,
    )
    .unwrap_or_else(|e| panic!("{e:?}"));
    assert!(d.body("a.bracket").is_some() && d.body("b.bracket").is_some() && d.assembly("a.root").is_some());
    // Instance names stay local: the path reads a.root/bracket.
    let mut d = d;
    assert_eq!(d.items("a.root/bracket").unwrap().len(), 1);
    assert!(d.params.iter().any(|p| p.name == "a.size") && d.params.iter().any(|p| p.name == "b.size"));
    let ink::lang::Value::Field(k, _) = d.def("k").unwrap().clone() else { panic!() };
    let ctx = Ctx::new(0.0, vec![10.0, 10.0]);
    assert!((eval_point(&Tape::new(&d.g, k), &ctx, [2.9, 0.0, 0.0]) + 0.1).abs() < 1e-9, "the let's own size (6) is used");
}

#[test]
fn namespaced_imports_preserve_language_heads() {
    let dir = crate::common::Scratch::new("ns-syntax");
    std::fs::write(
        dir.join("inner.ink"),
        r#"
        (import "std") (def shift 40)
        (fn shifted (f) (move f shift 0 0))
    "#,
    )
    .unwrap();
    std::fs::write(
        dir.join("family.ink"),
        r#"
        (import "std")
        (fn move (f dx dy dz) (warp f (- x dx) (- y dy) (- z dz)))
        (fn rot-x () 101) (fn rot-y () 102) (fn rot-z () 103) (fn rot () 104)
        (fn shape () 105) (fn instance () 107)
        (fn limits () 108) (fn value () 109) (fn accuracy () 110)
        (fn def () 111) (fn body () 112) (fn assembly () 113) (fn assert () 114)
        (fn frame () 118)
        (fn shifted (f) (move f 4 0 0))
    "#,
    )
    .unwrap();
    std::fs::write(
        dir.join("part.ink"),
        r#"
        (import "family.ink")
        (import "inner.ink" as p)
        (param shift 3)
        (def f (shifted (sphere 1)))
        (def g (p.shifted (sphere 1)))
        (body part (shape f) (attr density 0.001)
          (frame seat (at (rot-x 90) (rot-y 90) (rot-z 90) (rot 0 0 1 90)
                          (move shift 0 0)))
          (frame datum (axes (origin shift 0 0) (x 1 0 0) (z 0 0 1))))
        (assembly root (instance p part (at (move shift 0 0))
                          (turn z 20) (limits -30 30)))
        (assert ok (value shift > 0) (accuracy 0.1))
        (def datum (frame-of part datum))
    "#,
    )
    .unwrap();
    let mut d = ink::load_str(r#"(import "part.ink" as p)"#, "t", &dir).unwrap_or_else(|e| panic!("{e:?}"));
    assert!(d.body("p.part").is_some() && d.assembly("p.root").is_some());
    let ink::lang::Value::Field(f, _) = d.def("p.f").unwrap() else { panic!() };
    let ctx = Ctx::new(0.0, vec![3.0]);
    assert_eq!(eval_point(&Tape::new(&d.g, *f), &ctx, [4.0, 0.0, 0.0]), -1.0);
    // p.shifted is a local call; p.p.shifted is the nested p's. An
    // already-qualified value such as p.shift must keep its own binding.
    let ink::lang::Value::Field(g, _) = d.def("p.g").unwrap() else { panic!() };
    assert_eq!(eval_point(&Tape::new(&d.g, *g), &ctx, [40.0, 0.0, 0.0]), -1.0);
    assert_eq!(d.items("p.root").unwrap().len(), 1);
    let ink::lang::Value::Frame(frame) = d.def("p.datum").unwrap() else { panic!() };
    assert_eq!(frame.t.map(|n| eval_point(&Tape::new(&d.g, n), &ctx, [0.0; 3])), [3.0, 0.0, 0.0]);
    // A missing member of the nested p must not fall back to the outer
    // module's own function with the same public name, p.missing.
    std::fs::write(
        dir.join("missing.ink"),
        r#"
        (import "inner.ink" as p) (fn missing () 99) (def bad (p.missing))
    "#,
    )
    .unwrap();
    let e = ink::load_str(r#"(import "missing.ink" as p)"#, "t", &dir).err().expect("no nested p.missing");
    assert!(e.message.contains("unknown function"), "{e:?}");
}

#[test]
fn a_generated_part_loads_under_a_namespace() {
    // A part as the catalog's generators write one: its record read into
    // `spec`, `g` and `f` defined from it, its shape a family function
    // mapping lambdas over the record's outlines, frames placed by
    // rotations and moves, a root placing it, a claim on it.
    let dir = crate::common::Scratch::new("ns-generated");
    std::fs::write(
        dir.join("family.ink"),
        r#"(import "std")
        (fn fam-profile (p) (subtract (polygon (get p "outline")) (if (> (len (get p "islands")) 0) (reduce min (map polygon (get p "islands"))) 1)))
        (fn fam-frame (g)
          (let ((roof (extrude-z (rounded-rect-2d (get g "length") (get g "width") 1) -3 0))
                (pockets (reduce min (map (lambda (p) (extrude-z (fam-profile p) -3.01 (- (get p "depth") 3))) (get g "pockets")))))
            (subtract roof pockets)))"#,
    )
    .unwrap();
    std::fs::write(
        dir.join("part.json"),
        r#"{"geometry": {"length": 30, "width": 20,
        "pockets": [{"outline": [[-5, -5], [5, -5], [5, 5], [-5, 5]], "islands": [[[-1, -1], [1, -1], [1, 1], [-1, 1]]], "depth": 1}]}}"#,
    )
    .unwrap();
    std::fs::write(
        dir.join("part.ink"),
        r#"(import "family.ink")
        (def spec (read-json "part.json"))
        (def g (get spec "geometry"))
        (def f (fam-frame g))
        (body part (shape f)
          (frame mount (at))
          (frame seat (at (rot-x 180) (move 0 0 -3)))
          (for i (range 2) (frame pilot-$i (at (move (* 8 i) 6 0))))
          (frame left (at (rot-x -90) (move 0 -10 -1.5))))
        (assembly root (instance driver part))
        (assert roof-material (subset (probe 10 0 -1) (shape-of part)))"#,
    )
    .unwrap();
    let d = ink::load_str(r#"(import "part.ink" as p)"#, "t", &dir).unwrap_or_else(|e| panic!("{e:?}"));
    assert!(d.body("p.part").is_some());
    assert!(d.assembly("p.root").is_some());
}

#[test]
fn namespaced_parts_keep_their_families() {
    // Two parts whose families both define `paint`, each imported under a
    // namespace: a family imported plainly by a namespaced file joins its
    // namespace (so the two never clash with each other or the importer),
    // nested `as` namespaces nest, and std stays global.
    let dir = crate::common::Scratch::new("ns-fam");
    std::fs::create_dir_all(dir.join("a")).unwrap();
    std::fs::create_dir_all(dir.join("b")).unwrap();
    // The family tags `ball`; the part also defs `ball` (renamed a.ball
    // under the namespace) as that construction: the name still resolves.
    std::fs::write(dir.join("a/fam.ink"), "(import \"std\") (fn paint (r) (min (tag ball (sphere r)) (move (sphere 1) 9 0 0)))").unwrap();
    std::fs::write(dir.join("a/part.ink"), "(import \"std\") (import \"fam.ink\") (body p (shape (paint 2))) (def ball (shape-of p ball))").unwrap();
    std::fs::write(dir.join("b/fam.ink"), "(import \"std\") (fn paint (r) (box r r r))").unwrap();
    std::fs::write(dir.join("b/inner.ink"), "(import \"std\") (def size 3) (fn twice (s) (* 2 s))").unwrap();
    // Calls through a nested namespace (in.twice) are renamed like its values.
    std::fs::write(dir.join("b/part.ink"), "(import \"std\") (import \"fam.ink\") (import \"inner.ink\" as in) (body p (shape (paint (in.twice (* 0.5 in.size)))))").unwrap();
    let d = ink::load_str(
        "(import \"std\") (import \"a/part.ink\" as a) (import \"b/part.ink\" as b) (fn paint (x) x)
         (def both (min (shape-of a.p) (move (shape-of b.p) 10 0 0)))",
        "t",
        &dir,
    )
    .unwrap_or_else(|e| panic!("{e:?}"));
    assert!(d.def("b.in.size").is_some(), "nested namespaces nest");
    assert!(d.def("a.ball").is_some(), "a construction tagged in the family resolves from the part");
    let ink::lang::Value::Field(f, _) = d.def("both").unwrap().clone() else { panic!() };
    let ctx = Ctx::new(0.0, vec![]);
    let t = Tape::new(&d.g, f);
    assert!((eval_point(&t, &ctx, [0.0, 0.0, 0.0]) + 2.0).abs() < 1e-9, "a's paint is a sphere");
    assert!((eval_point(&t, &ctx, [10.0, 0.0, 0.0]) + 1.5).abs() < 1e-9, "b's paint is a box of size 3");
}

/// A file's top-level `for` is expanded when it is reached:
/// its list reads the names the file defines before it, a `for` among the
/// forms it makes reads what those made, and so do `for`s in a body's and
/// an assembly's clauses. Under a namespace, what a `for` makes is the
/// file's from the `for` on; a read of it before the `for` is an error,
/// never the importer's name of the same text.
#[test]
fn a_top_level_for_reads_the_names_before_it() {
    let num = |d: &Doc, name: &str| match d.def(name) {
        Some(ink::lang::Value::Field(f, _)) => d.g.as_const(*f).unwrap_or_else(|| panic!("{name} is no number")),
        _ => panic!("no def {name}"),
    };
    let d = doc("(def n 3) (for i (range n) (def v-$i (* 10 i)))
         (for i (range 2) (def m-$i (+ i 1)) (for j (range m-$i) (def c-$i-$j (+ (* 10 i) j))))
         (body b (shape (sphere 1)) (for i (range n) (attr a-$i i)))
         (assembly root (for i (range n) (instance u-$i b (at (move (* 3 i) 0 0)))))");
    assert_eq!((num(&d, "v-2"), num(&d, "c-1-1")), (20.0, 11.0));
    assert!(d.def("c-0-1").is_none(), "m-0 is 1: one c-0-j");
    assert_eq!((d.body("b").unwrap().attrs.len(), d.assembly("root").unwrap().instances.len()), (3, 3));
    let e = ink::load_str("(for i (range n) (def v-$i i))\n(def n 3)", "t", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("unknown name 'n': it is defined at line 2, after this") && e.line == 1, "{e:?}");

    let dir = crate::common::Scratch::new("ns-for");
    std::fs::write(
        dir.join("part.ink"),
        "(import \"std\") (def n 2) (for i (range n) (def w-$i (* 10 i))) (def total (+ w-1 n))
         (body unit (shape (box 1 1 1))) (assembly root (for i (range n) (instance u-$i unit)))",
    )
    .unwrap();
    let d = ink::load_str("(import \"part.ink\" as p) (def n 7) (def w-1 100)", "t", &dir).unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!((num(&d, "p.total"), num(&d, "p.w-1"), num(&d, "w-1")), (12.0, 10.0, 100.0));
    assert_eq!(d.assembly("p.root").unwrap().instances.len(), 2);
}

/// Names are lexical: a function reads them where it is written (its file,
/// what that file imports, std), never where it is called. A file may
/// shadow a std name once, and its name means its own definition
/// throughout it; std keeps reading its own.
#[test]
fn names_resolve_where_they_are_written() {
    let dir = crate::common::Scratch::new("scope");
    std::fs::create_dir_all(dir.join("lib")).unwrap();
    std::fs::write(dir.join("part.ink"), "(import \"std\") (fn bracket () (box 10 10 10))").unwrap();
    std::fs::write(dir.join("lib/part.ink"), "(fn part-radius () (get (read-json \"data.json\") \"radius\")) (def part-reader (lambda () (part-radius)))").unwrap();
    std::fs::write(dir.join("lib/data.json"), r#"{"radius": 1}"#).unwrap();
    std::fs::write(dir.join("data.json"), r#"{"radius": 100}"#).unwrap();
    std::fs::write(dir.join("a.ink"), "(import \"b.ink\")").unwrap();
    std::fs::write(dir.join("b.ink"), "(import \"a.ink\")").unwrap();
    let at = |src: &str, p: [f64; 3]| {
        let d = ink::load_str(src, "t", &dir).unwrap_or_else(|e| panic!("{e:?}"));
        let ink::lang::Value::Field(f, _) = d.def("v").unwrap().clone() else { panic!() };
        eval_point(&Tape::new(&d.g, f), &Ctx::new(0.0, d.params.iter().map(|p| p.default).collect()), p)
    };
    // The importer's own box (a corner-anchored one), defined after the
    // import or before it, does not change the part's: centred, -5 at 0.
    let corner = "(fn box (w d h) (max (- x w) (neg x) (- y d) (neg y) (- z h) (neg z)))";
    for src in [format!("(import \"std\") (import \"part.ink\" as p) {corner} (def v (p.bracket))"), format!("(import \"std\") {corner} (import \"part.ink\") (def v (bracket))")] {
        assert_eq!(at(&src, [0.0; 3]), -5.0, "{src}");
    }
    // Separate files each get their own shadow, even when imported without
    // namespaces. A lambda also keeps the file whose free names it reads.
    std::fs::write(dir.join("left.ink"), r#"(import "std") (fn box () 11) (fn left () (box)) (def left-lambda (lambda () (box)))"#).unwrap();
    std::fs::write(dir.join("right.ink"), r#"(import "std") (param box 22) (fn right () box)"#).unwrap();
    assert_eq!(at("(import \"left.ink\") (import \"right.ink\") (def box 33) (def v (+ (left) (right) (apply left-lambda (vec)) box))", [0.0; 3]), 77.0);
    let e = ink::load_str("(import \"left.ink\") (import \"right.ink\") (def box 33) (def box 44)", "t", &dir).err().unwrap();
    assert!(e.message.contains("already defined"), "{e:?}");
    // Two shadows can also precede any explicit std import.
    std::fs::write(dir.join("early.ink"), "(def sphere 1)").unwrap();
    assert_eq!(at(r#"(import "early.ink") (param sphere 2) (def v sphere)"#, [0.0; 3]), 2.0);
    // std's polygon reads std's here2, not a user param of that name.
    let sq = "(def v (polygon (vec (vec -1 -1) (vec 1 -1) (vec 1 1) (vec -1 1))))";
    assert_eq!(at(&format!("(import \"std\") (param here2 5) {sq}"), [0.0; 3]), -1.0);
    // A path read in a function is its file's.
    assert_eq!(at("(import \"lib/part.ink\") (def v (apply part-reader (vec)))", [0.0; 3]), 1.0);
    let err = |src: &str| ink::load_str(src, "t", &dir).err().unwrap_or_else(|| panic!("{src} loads"));
    // One definition per name outside std (a shadow of std's included).
    for src in ["(import \"std\") (def sphere 1) (def sphere 2)", "(import \"std\") (param box 1) (def box 2)", "(import \"std\") (fn box (a) 1) (fn box (a) 2)"] {
        assert!(err(src).message.contains("is already defined, at t:1"), "{src}");
    }
    // A function reads its file's names wherever they are written, once
    // they are defined; a file's own box is never std's before it.
    assert_eq!(at("(import \"std\") (fn f () (g)) (fn g () 2) (def v (f))", [0.0; 3]), 2.0);
    for src in ["(import \"std\") (fn f () (g))\n(def v (f))\n(fn g () 1)", &format!("(import \"std\")\n(def v (box 1 1 1))\n{corner}")] {
        let e = err(src);
        assert!(e.message.contains("it is defined at line 3, after this"), "{src}: {e:?}");
    }
    // Nor does a file read its importer's names.
    std::fs::write(dir.join("reads-up.ink"), "(fn up () size)").unwrap();
    let e = err("(def size 3) (import \"reads-up.ink\") (def v (up))");
    assert!(e.message.contains("unknown name 'size'") && e.hint.as_deref().unwrap_or("").contains("in a file this one does not import"), "{e:?}");
    let e = err("(import \"a.ink\")");
    assert!(e.message.starts_with("import cycle: ") && e.message.contains("a.ink -> ") && e.message.contains("b.ink -> ") && e.message.ends_with("a.ink"), "{e:?}");
}

/// A fn or lambda whose last argument is followed by `...` takes the rest
/// as a vec, so std's union, intersect and subtract take any number of
/// shapes: (subtract a b c) is a with b and c cut away.
#[test]
fn rest_arguments_make_combinators_variadic() {
    let ctx = Ctx::new(0.0, vec![]);
    // Slabs along x: [0, 10] minus [2, 3] and [6, 7], and so on.
    let slab = |a: f64, b: f64| format!("(max (- {a} x) (- x {b}))");
    let (s0, s1, s2) = (slab(0.0, 10.0), slab(2.0, 3.0), slab(6.0, 7.0));
    for (src, inside, outside) in [
        (format!("(subtract {s0} {s1} {s2})"), vec![1.0, 4.0, 9.0], vec![2.5, 6.5, 11.0]),
        (format!("(min {s1} {s2} {})", slab(20.0, 21.0)), vec![2.5, 6.5, 20.5], vec![1.0, 4.0, 9.0]),
        (format!("(max {s0} {} {})", slab(-5.0, 4.0), slab(3.0, 30.0)), vec![3.5], vec![2.5, 4.5]),
        (format!("(subtract {s0} {s1})"), vec![1.0, 4.0, 6.5], vec![2.5, 11.0]),
    ] {
        let (_d, tape) = tape_of(&src);
        for &px in &inside {
            assert!(eval_point(&tape, &ctx, [px, 0.0, 0.0]) < 0.0, "{src}: {px} should be inside");
        }
        for &px in &outside {
            assert!(eval_point(&tape, &ctx, [px, 0.0, 0.0]) > 0.0, "{src}: {px} should be outside");
        }
    }
    let d = doc("(fn count (a rest ...) (+ a (len rest))) (def n (count 1 2 3 4)) (def m (count 1)) (def k ((lambda (xs ...) (len xs))))");
    for (name, v) in [("n", 4.0), ("m", 1.0), ("k", 0.0)] {
        let ink::lang::Value::Field(f, _) = d.def(name).unwrap().clone() else { panic!() };
        assert_eq!(d.g.as_const(f), Some(v), "{name}");
    }
    let err = |src: &str| ink::load_str(&format!("(import \"std\")\n{src}"), "t", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{src} loads")).message;
    assert!(err("(fn two (a b more ...) a) (def a (two (sphere 1)))").contains("'two' takes 2 or more argument(s) (a b more ...), got 1"));
    assert!(err("(fn f (a ... b) a)").contains("'...' goes only at the end"));
    assert!(err("(fn f (...) 1)").contains("'...' follows the argument"));
}

/// A body's material is something physics can look up: a std/materials.ink
/// name or data. Free text is warned about where it is written, not when a
/// study reads it; it belongs in material-note. Another spelling of a
/// material (a case, a space, a generic name) is refused with the name to
/// write: one name per material, not a silent mapping.
#[test]
fn a_material_that_is_not_in_the_library_is_warned_about() {
    let d = doc(concat!(
        "(body a (shape (sphere 1)) (attr material \"aluminum-6061\") (attr material-note \"6061-T6, anodised\"))\n",
        "(body b (shape (sphere 1)) (attr material (data (E 1) (nu 0.3))))\n",
        "(body c (shape (sphere 1)) (attr material \"ABS nominal; not characterized\"))"
    ));
    assert_eq!(d.load_warnings.len(), 1, "{:?}", d.load_warnings);
    let w = &d.load_warnings[0];
    assert!(w.starts_with("test:4: body 'c': material 'ABS nominal; not characterized' is not in std/materials.ink"), "{w}");
    assert!(w.contains("aluminum-6061") && w.contains("material-note"), "{w}");
    let refused = |src: &str| ink::load_str(&format!("(import \"std\")\n{src}"), "t", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{src} loads"));
    for (src, to) in [
        ("(body a (shape (sphere 1)) (attr material \"Aluminium\"))", "aluminum-6061"),
        ("(body a (shape (sphere 1)) (attr material \"PETG\"))", "petg"),
        ("(def m (material delrin))", "pom"),
    ] {
        let e = refused(src);
        assert!(e.message.contains(&format!("is spelled '{to}'")) && e.hint.as_deref().is_some_and(|h| h.contains(&format!("write {to}"))), "{src}: {e:?}");
    }
}

/// A body's look is read once, at load (ink::lang::look): "#RRGGBB"
/// becomes three numbers, fields stay fields, and a malformed one (five hex
/// digits, a colour in 0..255, a word for a number) is refused at its line,
/// not taken as the palette without a word. Clients read the look by name.
#[test]
fn looks_are_read_at_load() {
    use ink::lang::look::Look;
    let d = doc("(body a (shape (sphere 1)) (attr color \"#FF8000\") (attr metallic 1) (attr roughness (select z 0 1)))");
    let look = Look::of(&d.body("a").unwrap().attrs);
    let rgb = look.color.unwrap().map(|c| d.g.as_const(c).unwrap());
    assert!(rgb[0] == 1.0 && (rgb[1] - 128.0 / 255.0).abs() < 1e-12 && rgb[2] == 0.0, "{rgb:?}");
    assert_eq!(look.metallic.and_then(|m| d.g.as_const(m)), Some(1.0));
    assert!(look.roughness.is_some_and(|r| d.g.as_const(r).is_none()) && look.transmission.is_none());
    let err = |attr: &str| {
        let e = ink::load_str(&format!("(import \"std\")\n(body a (shape (sphere 1))\n  {attr})"), "t", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{attr} loads"));
        assert_eq!(e.line, 3, "{e:?}");
        e.message
    };
    assert!(err("(attr color \"#12345\")").contains("attr color: color \"#12345\" is not a colour: give \"#RRGGBB\""));
    assert!(err("(attr color (vec 200 30 30))").contains("color's red 200 is outside 0..1"));
    assert!(err("(attr color (vec 1 0))").contains("is not a colour"));
    assert!(err("(attr metallic \"shiny\")").contains("metallic \"shiny\" is not a number or a scalar field in 0..1"));
}

/// Inputs are never coerced or dropped without a word; these are errors:
/// counts and indices are whole numbers, checked against the limit before
/// anything is allocated; string escapes and numbers out of range; a
/// lambda's argument twice; an op given vecs of other
/// lengths, an empty one too, or the wrong number of arguments; a cell size
/// not above 0 and finite; a sample past a grid's f32s.
#[test]
fn inputs_are_never_silently_coerced() {
    let err = |src: &str| ink::load_str(&format!("(import \"std\")\n{src}"), "t", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{src} loads")).message;
    for (src, says) in [
        ("(def a (range 1e308))", "range would produce 1e308 elements (limit 1e6)"),
        ("(def a (range 0 1 1e-9))", "range would produce 1000000000 elements"),
        ("(def a (linspace 0 1 1e308))", "limit 1e6"),
        ("(def a (range 0.9 2.9))", "range counts in whole numbers, got 0.9"),
        ("(def a (nth (vec 10 20) 0.9))", "an nth index must be a whole number"),
        ("(def a (slice (vec 10 20) 0.5))", "a slice index must be a whole number"),
        ("(def a \"C:\\part.stl\")", "unknown escape '\\p'"),
        ("(def a 1e999)", "outside the range of a double"),
        ("(def a 1e-400)", "outside the range of a double"),
        ("(def a ((lambda (a a) a) 1 2))", "duplicate argument 'a'"),
        ("(def a (+ (vec) (vec 1 2)))", "on vecs of lengths 0 and 2"),
        ("(def a (select (vec) 1))", "'select' takes 3 argument(s)"),
        ("(def a (- (vec) 1 2))", "'-' takes 2 argument(s)"),
        ("(def a (range (exp 1000) 0))", "range needs finite arguments"),
        ("(def a (sqrt (vec) 2))", "'sqrt' takes 1 argument(s), got 2"),
        ("(def a (grid-bake (sphere 1) (res -1)))", "(res h) expects a cell size above 0 and finite, got -1"),
        ("(def a (grid-bake (* 1e100 (sphere 1)) (res 1)))", "past the range of the grid's numbers"),
        ("(body b (shape (sphere 1))) (def a (+ b 1))", "'b' is a body, not a value"),
    ] {
        let m = err(src);
        assert!(m.contains(says), "{src}: {m}");
    }
    let ok = ink::load_str("(import \"std\") (def a (let ((r 1) (r (+ r 1))) r)) (def s \"a\\\"b\\\\c\") (def k 1e-300)", "t", Path::new(ROOT));
    assert!(ok.is_ok(), "a sequential let rebinds; known escapes and small numbers read: {:?}", ok.err());
    // A recursion's frames are listed once, with how many times.
    let e = ink::load_str("(fn f (a) (f a)) (def b (f 1))", "t", Path::new(ROOT)).err().unwrap();
    let frames = e.to_json()["called_from"].as_array().unwrap().clone();
    assert!(frames.len() == 2 && frames[0]["times"] == 255, "{frames:?}");
}

#[test]
fn an_empty_clause_is_an_error_at_it() {
    // () as an optimize's clause, or a metric (a claim's:
    // degenerate_forms_are_errors_not_panics): a clause with no name, which
    // every reader splits the same way. (A study's: ext's sim tests.)
    for src in ["(param p 1) (optimize o (vary p 0 2) ())", "(param p 1) (optimize o (vary p 0 2) (minimize ()))"] {
        let e = ink::load_str(&format!("(import \"std\")\n{src}"), "t", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{src} loaded"));
        assert_eq!(e.line, 2, "{src}: {e:?}");
    }
}

#[test]
fn errors_point_at_source_with_hints() {
    let e = ink::load_str("(import \"std\")\n(def a (sphere 1 2))", "t", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("takes 1 argument"), "{}", e.message);
    assert_eq!(e.line, 2);
    let e = ink::load_str("(bodyy b (shape x))", "t", Path::new(ROOT)).err().unwrap();
    assert!(e.hint.unwrap().contains("body"));
    let e = ink::load_str("(def a (sphere 1))", "t", Path::new(ROOT)).err().unwrap();
    assert!(e.hint.unwrap().contains("import"));
    let e = ink::load_str("(fn f (a) (f a)) (def b (f 1))", "t", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("recursion"));
    // An error deep inside std leads back to the author's call.
    let e = ink::load_str("(import \"std\")\n\n(def b (bezier3-points (vec 0 0 0) (vec 1 0 0) (vec 1 1 0) (vec 2 2 0) \"x\"))", "t", Path::new(ROOT)).err().unwrap();
    assert_eq!(e.file, "std");
    let (f, file, line, _) = e.called_from.last().unwrap();
    assert_eq!((f.as_str(), file.as_str(), *line), ("bezier3-points", "t", 3), "{:?}", e.called_from);
}

#[test]
fn gradient_bounds_are_sound() {
    let ctx = Ctx::new(0.0, vec![]);
    let mut rng = Lcg(23);
    for src in ZOO {
        let (_d, tape) = tape_of(src);
        for _ in 0..200 {
            let c = [rng.range(-20.0, 20.0), rng.range(-20.0, 20.0), rng.range(-20.0, 20.0)];
            let s = rng.range(0.01, 6.0);
            let b = [I::new(c[0] - s, c[0] + s), I::new(c[1] - s, c[1] + s), I::new(c[2] - s, c[2] + s)];
            let g = eval_interval_grad(&tape, &ctx, b);
            let pts: Vec<[f64; 3]> = (0..64).map(|_| [rng.range(b[0].lo, b[0].hi), rng.range(b[1].lo, b[1].hi), rng.range(b[2].lo, b[2].hi)]).collect();
            for (p, gr) in pts.iter().zip(eval_grads(&tape, &ctx, &pts)) {
                if gr.iter().any(|v| !v.is_finite()) {
                    continue;
                }
                let tol = 1e-9 * (1.0 + gr[0].abs());
                assert!(g.v.lo - tol <= gr[0] && gr[0] <= g.v.hi + tol, "{src}: value at {p:?}");
                for k in 0..3 {
                    let d = gr[k + 1];
                    let t = 1e-9 * (1.0 + d.abs());
                    assert!(g.d[k].lo - t <= d && d <= g.d[k].hi + t, "{src}: d/d{k} = {d} at {p:?} outside [{}, {}]", g.d[k].lo, g.d[k].hi);
                }
            }
        }
    }
}

#[test]
fn undefined_cases_are_reported_with_their_construction() {
    // Values follow the core rules (x / 0 = 0, sqrt of a negative is 0), and
    // every undefined case the evaluation passes through is reported.
    let d = doc("(def f (+ (tag ratio (/ y x)) (tag log (ln (- 1 (square z)))) (sqrt (- 1 (square z)))))");
    let ink::lang::Value::Field(f, _) = d.def("f").unwrap().clone() else { panic!() };
    let t = Tape::new(&d.g, f);
    let ctx = Ctx::new(0.0, vec![]);
    assert_eq!(eval_point(&t, &ctx, [0.0, 3.0, 0.0]), 1.0);
    assert_eq!(eval_point(&t, &ctx, [0.0, 3.0, 2.0]), 0.0);
    // sqrt of a negative is the continuous sqrt(max(a, 0)): not reported.
    let u = undefined_cases(&t, &ctx, &[[0.0, 3.0, 0.0], [0.0, 3.0, 2.0], [1.0, 3.0, 0.5]]);
    let name = |u: &Undefined| u.name.map(|k| d.g.tag_names[k as usize].clone());
    let div = u.iter().find(|u| u.what == "division by zero").unwrap();
    assert_eq!((div.count, name(div).as_deref()), (2, Some("ratio")));
    let log = u.iter().find(|u| u.what == "ln of a negative").unwrap();
    assert_eq!((log.count, name(log).as_deref(), log.at), (1, Some("log"), [0.0, 3.0, 2.0]));
    assert_eq!(u.len(), 2);
}

#[test]
fn an_undefined_case_in_constants_is_reported_where_a_value_depends_on_it() {
    // A zero-length segment's projection is 0 / 0, multiplied by the
    // segment (0): discarded, so its distance is a circle's, unreported.
    let d = doc("(def c (- (segment-distance-at here2 (vec 1 0) (vec 1 0)) 1))");
    assert!(d.load_warnings.is_empty(), "{:?}", d.load_warnings);
    let ink::lang::Value::Field(f, _) = d.def("c").unwrap().clone() else { panic!() };
    let (t, ctx) = (Tape::new(&d.g, f), Ctx::new(0.0, vec![]));
    assert_eq!(eval_point(&t, &ctx, [1.0, 2.0, 0.0]), 1.0);
    assert!(undefined_cases(&t, &ctx, &[[1.0, 2.0, 0.0], [5.0, 0.0, 0.0]]).is_empty());
    // One a value depends on: in a field, at its samples; in a number used
    // at load time, at its line.
    let d = doc("(def f (+ (sphere 1) (ln -1)))\n(for i (range (/ 3 0)) (def k i))");
    let ink::lang::Value::Field(f, _) = d.def("f").unwrap().clone() else { panic!() };
    let u = undefined_cases(&Tape::new(&d.g, f), &Ctx::new(0.0, vec![]), &[[0.0; 3]]);
    assert_eq!(u.iter().map(|u| u.what).collect::<Vec<_>>(), ["ln of a negative"]);
    assert!(d.load_warnings.len() == 1 && d.load_warnings[0].starts_with("test:3: division by zero"), "{:?}", d.load_warnings);
}

#[test]
fn distances_are_one_lipschitz_from_the_rules_alone() {
    // The gradient bound of sqrt of a sum of squares comes from the Gram
    // matrix sum grad u_i grad u_i^T: exactly 1 for the distance to a point,
    // rotated or not, with no (lipschitz ...) promise (not sqrt 3).
    use ink::interval::{I, lipschitz};
    let ctx = Ctx::new(0.0, vec![]);
    let b = [I::new(1.0, 2.0), I::new(-3.0, -1.0), I::new(0.5, 4.0)];
    for src in
        ["(- (sqrt (+ (square x) (square y) (square z))) 3)", "(let ((u (+ (* 0.6 x) (* 0.8 y))) (v (- (* 0.8 x) (* 0.6 y)))) (- (sqrt (+ (square u) (square v) (square z))) 3))"]
    {
        let (_d, t) = tape_of(src);
        let l = lipschitz(&eval_interval_grad(&t, &ctx, b));
        assert!((1.0..1.0 + 1e-9).contains(&l), "{src}: {l}");
    }
}

/// Every path takes one derivative: the symbolic one (`diff`, a field
/// evaluated like any other), the forward gradient (eval_grads) and dual
/// numbers (eval_dual) each write an op's derivative as its partials, from
/// its operands' values, times its operands' derivatives, summed in order
/// (a tie's mean da/2 + db/2; a read's partials, the reads one axis
/// further), so they agree bit for bit, at every scale. Cases: a tie of
/// slopes 5e-324, a read warped by 1e308 x, and quotients, roots and logs.
#[test]
fn every_path_takes_one_derivative() {
    let fields = [
        "(warp baked (* 1e308 x) 0 0)",
        "(warp baked (+ (* 1.7 x) (* 0.3 y)) (* 1.3 (- y (* 0.2 z))) (+ z (* 0.4 x) t))",
        "(diff (warp baked (* 2 x) (* 3 y) (+ z (* p x))) x)",
        "(diff (diff (warp baked (* x y) (- y z) (* p z)) y) z)",
        "(warp baked (* 1e200 (sin x)) (* 1e-200 y) (/ z (+ 1 (square p))))",
        "(min (* p x) (+ (* p x) y))",
        "(max (* p x) (- (* p x) y))",
        "(* 1e300 (min (* p x) (* p y)))",
        "(min (* 1e308 x) (+ (* 1.5e308 x) y))",
        "(sqrt (+ 1e-300 (* 1e300 (square x))))",
        "(/ (* 1e300 x) (+ y (* 1e-300 z) 3))",
        "(/ (+ x 1) (* 1e-200 (+ y 2)))",
        "(ln (+ 2 (* 1e-300 x) (* 1e300 (square y))))",
        "(atan2 (* 1e200 x) (+ y (* p z)))",
        "(atan (* 1e160 (+ x y)))",
        "(emax (- (abs x) 1) (- (abs (* p y)) 1))",
        "(pow (+ 2 (abs x)) (+ y p))",
        "(smooth-union (sphere 1) (move (box 1 1 1) 1 0 0) 0.5)",
        "(len2 (- (atan2 y x) pi) z)",
    ];
    let vars = ["x", "y", "z", "t", "p"];
    let mut src = "(param p 0.5) (def baked (grid-bake (rotate-z (box 3 2 1.5) 20) (res 0.25)))\n".to_string();
    for (i, f) in fields.iter().enumerate() {
        src += &format!("(def f{i} {f})\n");
        for v in vars {
            src += &format!("(def f{i}{v} (diff f{i} {v}))\n");
        }
    }
    let d = doc(&src);
    let tape = |d: &Doc, name: &str| {
        let ink::lang::Value::Field(f, _) = d.def(name).unwrap().clone() else { panic!("{name}") };
        Tape::new(&d.g, f)
    };
    let mut rng = Lcg(5);
    let mut pts: Vec<[f64; 3]> = vec![[0.0; 3], [1e-307, 0.0, 0.0], [-1.5e-308, 2e-308, 5e-324], [1.0, 1.0, 1.0], [0.0, 0.25, -0.5], [2.7e-162, 0.0, 0.0]];
    pts.extend((0..60).map(|_| [rng.range(-3.0, 3.0), rng.range(-3.0, 3.0), rng.range(-3.0, 3.0)]));
    for (t, p) in [(0.25, 5e-324), (-0.5, 0.75), (1e-310, 1e300)] {
        let ctx = d.ctx(&[("p".into(), p.into())], t).unwrap();
        for (i, f) in fields.iter().enumerate() {
            let tf = tape(&d, &format!("f{i}"));
            let dv: Vec<Tape> = vars.iter().map(|v| tape(&d, &format!("f{i}{v}"))).collect();
            let grads = eval_grads(&tf, &ctx, &pts);
            for (at, g) in pts.iter().zip(&grads) {
                let dual = eval_dual(&tf, &ctx, *at);
                assert_eq!(g[0], eval_point(&tf, &ctx, *at), "{f} at {at:?}");
                assert_eq!(dual[0], g[0], "{f} at {at:?}");
                for (k, v) in vars.iter().enumerate() {
                    let sym = eval_point(&dv[k], &ctx, *at);
                    assert!(!sym.is_nan());
                    assert_eq!(sym, dual[1 + k], "{f}: d/d{v} at {at:?}, t = {t}, p = {p}: symbolic against dual");
                    if k < 3 {
                        assert_eq!(sym, g[1 + k], "{f}: d/d{v} at {at:?}, t = {t}, p = {p}: symbolic against gradient");
                    }
                }
            }
        }
    }
    // The repros: the tie's slope is 0 (its halves round to 0), and the read
    // warped by 1e308 x has the finite slope 1e308 (its partial, about 1,
    // times the rate).
    let ctx = d.ctx(&[("p".into(), 5e-324.into())], 0.0).unwrap();
    assert_eq!(eval_point(&tape(&d, "f5x"), &ctx, [0.0; 3]), 0.0);
    let slope = eval_grads(&tape(&d, "f0"), &ctx, &[[1e-307, 0.0, 0.0]])[0][1];
    assert!(slope.is_finite() && (0.99e308..1.01e308).contains(&slope), "{slope}");
}

#[test]
fn sampled_slope_ranges_are_finite() {
    // Proofs get gradient bounds on sampled fields: over a small box, inside
    // the grid or beyond it, a sampled field's slope ranges are finite.
    // (That they hold its gradients: gradient_bounds_are_sound, on this
    // field; that its derivative paths agree: every_path_takes_one_derivative.)
    let (_d, t) = tape_of("(move (rotate-z (grid-bake (torus-z 8 2) (res 0.5)) 30) 3 -2 1)");
    let ctx = Ctx::new(0.0, vec![]);
    let mut rng = Lcg(11);
    for _ in 0..400 {
        let p = [rng.range(-15.0, 15.0), rng.range(-15.0, 15.0), rng.range(-6.0, 8.0)];
        let gb = eval_interval_grad(&t, &ctx, p.map(|v| I::new(v - 0.2, v + 0.2)));
        assert!(gb.d.iter().all(|d| d.lo.is_finite() && d.hi.is_finite()), "unbounded slope at {p:?}: {:?}", gb.d);
    }
}

#[test]
fn ranges_hold_at_box_corners() {
    // Corner cases of the range rules: each range must hold its op's value
    // at the box's corners.
    use ink::graph::{BinOp as B, UnOp as U};
    use ink::interval::{G, I, bin, g_bin, g_select, g_un, lipschitz};
    // pow(0, -1) is +inf, not 1 / 0 (which the value rules make 0).
    assert!(bin(B::Pow, I::point(0.0), I::point(-1.0)).0.contains(f64::INFINITY));
    // A negative base over a narrow exponent range holding a whole number.
    let r = bin(B::Pow, I::point(-0.5), I::new(47.0, 47.0f64.next_up())).0;
    assert!(r.contains(-(0.5f64.powi(47))), "{r:?}");
    // atan2 jumps by 2 pi across the negative x axis: no slope bound there.
    let x = G::coordinate(I::new(-0.01, 0.01), 0);
    assert!(lipschitz(&g_bin(B::Atan2, x, G::constant(-1.0))).is_infinite());
    // min and max of the same arguments have equal bounds, yet a select
    // between them jumps where its condition changes sign.
    let (x, y, z) = (G::coordinate(I::new(-1.0, 1.0), 0), G::coordinate(I::new(-1.0, 1.0), 1), G::coordinate(I::new(-1.0, 1.0), 2));
    assert!(lipschitz(&g_select(z, g_bin(B::Min, x, y), g_bin(B::Max, x, y), false)).is_infinite());
    // Underflow is not zero: slopes of tiny scale stay above 0.
    let tiny = g_bin(B::Mul, G::coordinate(I::new(1e200, 2e200), 0), G::constant(1e-200));
    let norm = g_un(U::Sqrt, g_un(U::Square, tiny));
    assert!(lipschitz(&norm) >= 1e-200 && norm.d[0].contains(1e-200), "{norm:?}");
    let q = g_bin(B::Div, G::coordinate(I::new(1e-200, 2e-200), 0), G::constant(1e-200));
    assert!(q.d[0].contains(1e200), "{:?}", q.d[0]);
}

/// atan2's range holds the real value at every corner of boxes whose x
/// crosses 0, of either sign of y (its y > 0 and y < 0 rules subtract
/// atan(x / y) from pi/2, which cancels near the +x axis; rounded to
/// nearest, with ulps of the result added at the end, the bound of
/// atan2(1e-10, [-1, 1]) would begin at 1.0000000827e-10, above
/// atan(1e-10) < 1e-10). A point x takes the x > 0 rule and misses this,
/// so the boxes are wide in x. The range holds the library's atan2 at
/// each corner (the real value is within an ulp of it, on the side the
/// rule's directed rounding keeps: its steps are exact but for libm's
/// atan, allowed 4 ulps); for 0 < t < 1, atan t < t - t^3/3 + t^5/5 (an
/// alternating series), below each case's limit.
#[test]
fn atan2_ranges_hold_near_the_axes() {
    use ink::graph::{BinOp, atan2_f64};
    use ink::interval::bin;
    let at = |y: I, x: I| bin(BinOp::Atan2, y, x).0;
    for (t, limit) in [(1e-10, 1e-10), (0.002, 0.00199999733333975)] {
        for s in [1.0, -1.0] {
            let r = at(I::point(s * t), I::new(-1.0, 1.0));
            assert!(if s > 0.0 { r.lo < limit } else { r.hi > -limit }, "atan2({}, [-1, 1]): {r:?}", s * t);
        }
    }
    let mags = [1e-300, 1e-20, 1e-10, 2e-3, 0.5, 1.0, 3.0, 1e10, 1e300];
    let mut bad = Vec::new();
    for &y0 in &mags {
        for &y1 in mags.iter().filter(|&&m| m >= y0) {
            for &a in &mags {
                for &b in &mags {
                    for s in [1.0, -1.0] {
                        let (y, x) = (if s > 0.0 { I::new(y0, y1) } else { I::new(-y1, -y0) }, I::new(-a, b));
                        let r = at(y, x);
                        for (yc, xc) in [(y.lo, x.lo), (y.lo, x.hi), (y.hi, x.lo), (y.hi, x.hi)] {
                            let v = atan2_f64(yc, xc);
                            if !(r.lo <= v && v <= r.hi) {
                                bad.push((y, x, (yc, xc), v, r));
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(bad.is_empty(), "{} corners outside their range, e.g. {:?}", bad.len(), &bad[..bad.len().min(3)]);
    // On the cut: y from 0 up, x < 0. The one zero reads atan2(0, x) as
    // +pi, so the range ends there, not at -pi.
    for (y, x) in [(I::point(0.0), I::new(-3.0, -1.0)), (I::new(0.0, 2.0), I::new(-3.0, -1.0)), (I::point(-0.0), I::point(-1.0))] {
        let r = at(y, x);
        for (yc, xc) in [(y.lo, x.lo), (y.lo, x.hi), (y.hi, x.lo), (y.hi, x.hi)] {
            let v = atan2_f64(yc, xc);
            assert!(r.lo <= v && v <= r.hi, "atan2({yc}, {xc}) = {v} outside {r:?}");
        }
        assert!(r.lo >= atan2_f64(y.hi, x.hi) - 1e-14 && r.hi <= PI + 1e-14, "{r:?}");
    }
    // Still tight: an angle bounded to within a few ulps of pi/2 near the
    // axis (the claim (value (atan2 0.002 p) > 0.0019999973), true by
    // 3e-14, is proven: tests/geometry.rs).
    let r = at(I::point(0.002), I::new(-1.0, 1.0));
    assert!(r.lo > 0.002 - 0.002f64.powi(3) / 3.0 - 1e-15 && r.hi < PI, "{r:?}");
}

/// A sampled field's mixed derivatives at extreme spacings: the point value
/// (eval_point, through Graph::diff) lies in the interval of the point,
/// finite where the derivative is (a component divided by h^k would
/// underflow to 0 at h = 1e-170, making corner / h^2 infinite, overflow at
/// h = 1e170, making it 0, and give 0 / 0 for zero data).
/// Inside the cell the field is c x y / h^2.
#[test]
fn sampled_derivatives_hold_at_extreme_spacings() {
    use ink::graph::{Graph, Var};
    use ink::grid::GridData;
    for (h, c, want) in [(1e-170, f32::from_bits(1), 1.401298464324817e295), (1e170, f32::MAX, 3.402823466385289e-302), (1e-170, 0.0, 0.0)] {
        let data = Arc::new(GridData::new([0.0; 3], h, [2; 3], vec![0.0, 0.0, 0.0, c, 0.0, 0.0, 0.0, c]));
        let p = [0.5 * h; 3];
        let raw = data.component(p, 3);
        let mut g = Graph::new();
        let at = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let f = g.grid(data, at);
        let fx = g.diff(f, Var::X);
        let fxy = g.diff(fx, Var::Y);
        let tape = Tape::new(&g, fxy);
        let ctx = Ctx::new(0.0, vec![]);
        let (v, r) = (eval_point(&tape, &ctx, p), eval_interval(&tape, &ctx, p.map(I::point), None));
        let close = |a: f64| a == want || (a - want).abs() <= 1e-12 * want.abs();
        assert!(close(raw) && close(v) && r.contains(v), "h {h:e}, corner {c:e}: component {raw:e}, point {v:e}, range {r:?}");
    }
}

#[test]
fn grid_slope_ranges_survive_cancellation() {
    // Huge samples beside small ones: the mixed difference -1 must stay in
    // the range though f64 sums round the small terms away.
    use ink::grid::GridData;
    use ink::interval::I;
    let big = 2f64.powi(54) as f32;
    let layer = [0.0f32, 1.0, big, big];
    let g = GridData::new([0.0; 3], 1.0, [2, 2, 2], layer.iter().chain(&layer).copied().collect());
    let v = g.component([0.5, 0.5, 0.0], 3);
    let b = g.component_bound([I::new(0.4, 0.6), I::new(0.4, 0.6), I::new(0.0, 0.0)], 3);
    assert!(b.contains(v) && b.contains(-1.0), "{v} {b:?}");
}

#[test]
fn sampled_ranges_are_exact_within_a_cell() {
    // Over part of a cell a trilinear interpolant's extremes are at that
    // part's corners: a box inside one cell gets the range of the box, not
    // of the cell (proofs within a cell of a sampled surface).
    use ink::grid::GridData;
    use ink::interval::I;
    // Field z - 0.5 sampled at spacing 1: exact (linear).
    let g = GridData::new([0.0; 3], 1.0, [2, 2, 2], vec![-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5]);
    let r = g.component_bound([I::new(0.2, 0.3), I::new(0.2, 0.3), I::new(0.71, 0.76)], 0);
    assert!(r.lo >= 0.21 - 1e-9 && r.hi <= 0.26 + 1e-9 && r.contains(0.21) && r.contains(0.26), "{r:?}");
}

#[test]
fn sampled_ranges_hold_at_the_points_they_bound() {
    // The range and the value share one local coordinate, widened for its
    // rounding: at origin 100, spacing 0.1, the point (101, 101, 101) lies
    // in its range.
    use ink::grid::GridData;
    let g = GridData::new([100.0; 3], 0.1, [2; 3], vec![0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
    let p = [101.0; 3];
    assert!(g.component_bound(p.map(I::point), 0).contains(g.component(p, 0)));
    // A sweep: non-binary spacings at nonzero origins, points on
    // cell faces, in the last cell and beyond it, every component.
    let mut bad = Vec::new();
    for origin in [0.0, 1.0, 10.0, 100.0, 1000.0, 1e6, 1e10] {
        for h in [0.1, 0.3, 1.0, 1e-6, 1e-100, 1e100] {
            for dims in [[2, 2, 2], [3, 3, 3], [1, 2, 3]] {
                let values = (0..dims.iter().product::<usize>()).map(|i| ((i * 37 % 17) as f32 - 8.0) * 1.3).collect();
                let g = GridData::new([origin; 3], h, dims, values);
                for mask in 0..8 {
                    for factor in [-1.0, 0.0, 0.3, 1.0, 1.5, 2.0, 3.0, 10.0] {
                        let p = [origin + factor * h; 3];
                        let (r, v) = (g.component_bound(p.map(I::point), mask), g.component(p, mask));
                        if v.is_finite() && !r.contains(v) {
                            bad.push((origin, h, dims, mask, factor, v, r));
                        }
                    }
                }
            }
        }
    }
    assert!(bad.is_empty(), "{} points outside their range, e.g. {:?}", bad.len(), &bad[..bad.len().min(3)]);
}

/// A read of sampled data is one node only when it is the same quantity of
/// the source: the same grid, derivative (the grid's axes, in order) and
/// coordinates. The samples agree where the source need not: their second
/// derivatives along x and along y are both 0 in a cell, so are their third
/// and fourth; their mixed derivatives in either order agree (the source's
/// need not: a distance has kinks); and two points beyond the grid box that
/// clamp onto one point of it, at one distance from it, read one sample.
/// Each stays its own node through warps, so their difference is not 0 and
/// names the grid, and their fingerprints differ; the same read made twice
/// is one node.
#[test]
fn reads_of_samples_are_one_node_only_when_the_same_quantity() {
    use ink::graph::{Fingerprint, Graph, Var};
    use ink::grid::{Band, GridData, Source};
    let h = 0.5;
    let band = Band { err: 0.25, outside: Some(0.5) };
    let source = Source { what: "cube.stl".into(), at: "test".into(), kind: "mesh", res: h, band: Some(band), measured: None, note: None, unit: "", span: None, component: None };
    let values: Vec<f32> = (0..729).map(|n| ((n % 9) as f32 - 4.0).abs() + 0.1 * (n / 81) as f32).collect();
    let data = Arc::new(GridData::new([-2.0; 3], h, [9; 3], values).with_source(source));
    let mut g = Graph::new();
    let xyz = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
    let d = g.grid(data.clone(), xyz);
    let dx = g.diff(d, Var::X);
    let dy = g.diff(d, Var::Y);
    let (dxx, dyy, dxy, dyx) = (g.diff(dx, Var::X), g.diff(dy, Var::Y), g.diff(dx, Var::Y), g.diff(dy, Var::X));
    let d3 = g.diff(dxx, Var::X);
    let d4 = g.diff(d3, Var::X);
    // Beyond the grid box [-2, 2]^3: both clamp onto (-2, -2, 0), 0.5 from it.
    let at = |g: &mut Graph, p: [f64; 3]| p.map(|v| g.constant(v));
    let (p1, p2) = (at(&mut g, [-2.5, -2.0, 0.0]), at(&mut g, [-2.0, -2.5, 0.0]));
    let (r1, r2) = (g.grid(data.clone(), p1), g.grid(data.clone(), p2));
    let ctx = Ctx::new(0.0, vec![]);
    let fp = |g: &Graph, f: u32| {
        let mut fp = Fingerprint::default();
        g.fingerprint(&[f], &[], false, &mut fp);
        fp.value()
    };
    let to = at(&mut g, [1.25, 1.3, 0.0]);
    for (a, b) in [(dxx, dyy), (d4, d3), (d4, dxx), (dxy, dyx), (r1, r2)] {
        let fixed = [a, b].map(|f| g.subst(f, &[(Var::X, to[0]), (Var::Y, to[1]), (Var::Z, to[2])]));
        for (a, b) in [(a, b), (fixed[0], fixed[1])] {
            let diff = g.sub(a, b);
            assert!(a != b && g.as_const(diff).is_none() && g.sources(&[diff]).len() == 1, "{a} {b}");
            assert_ne!(fp(&g, a), fp(&g, b));
            // On the samples they agree.
            let p = [0.1, 0.2, 0.0];
            assert_eq!(eval_point(&Tape::new(&g, a), &ctx, p), eval_point(&Tape::new(&g, b), &ctx, p), "{a} {b}");
        }
    }
    let again = g.diff(dx, Var::X);
    let r1_again = g.grid(data, p1);
    assert!(again == dxx && r1_again == r1);
}

#[test]
fn linear_bounds_are_sound_over_boxes_and_far_sectors() {
    // eval_interval_linear bounds values by collected linear forms too: of
    // rotated coordinates, of abs where its argument is one-signed, of a
    // min or max member that wins over the whole region, and of terms in
    // time (constant over the region). Every value at a point of the region,
    // at any moment of the span, must lie in its bound.
    let extra = [
        "(rotate-x (max (- (abs y) (+ 4 (* 0.44 z))) (slab-z -2.5 0)) 180)",
        "(move (rotate-y (box 10 20 30) (* 40 t)) (* 5 (sin t)) 0 (* 3 t))",
        "(min (abs (- x (* 3 y))) (max (+ x z) (- y 2)) (abs (+ (* 2 z) t)))",
    ];
    let span = Ctx { t_span: Some((0.0, 1.0)), ..Ctx::new(0.5, vec![]) };
    let mut rng = Lcg(19);
    for src in ZOO.iter().chain(&extra) {
        let (_d, tape) = tape_of(src);
        for k in 0..200 {
            let mut pts = Vec::new();
            let region = if k % 2 == 0 {
                let c = [rng.range(-20.0, 20.0), rng.range(-20.0, 20.0), rng.range(-20.0, 20.0)];
                let s = rng.range(0.01, 8.0);
                for _ in 0..32 {
                    pts.push(c.map(|v| rng.range(v - s, v + s)));
                }
                Region::Box(c.map(|v| I::new(v - s, v + s)))
            } else {
                // A far sector: directions q on one face of the cube, a
                // range of q on the other two axes, distances s out from c.
                let c = [rng.range(-20.0, 20.0), rng.range(-20.0, 20.0), rng.range(-20.0, 20.0)];
                let r = [rng.range(1.0, 30.0), rng.range(1.0, 30.0), rng.range(1.0, 30.0)];
                let (axis, side) = ((k / 2) % 3, if rng.next() < 0.5 { -1.0 } else { 1.0 });
                let q: [I; 3] = std::array::from_fn(|a| {
                    if a == axis {
                        I::point(side)
                    } else {
                        let (u, w) = (rng.range(-1.0, 1.0), rng.range(0.0, 0.5));
                        I::new(u, (u + w).min(1.0))
                    }
                });
                let s0 = rng.range(1.0, 100.0);
                let s1 = if rng.next() < 0.3 { f64::INFINITY } else { s0 * rng.range(1.0, 5.0) };
                for _ in 0..32 {
                    let s = rng.range(s0, s1.min(1e3 * s0));
                    pts.push(std::array::from_fn(|a| c[a] + r[a] * s * rng.range(q[a].lo, q[a].hi)));
                }
                Region::Sector { c, r, q, s: I::new(s0, s1) }
            };
            let b = eval_interval_linear(&tape, &span, &region, None);
            for p in pts {
                let t = rng.next();
                let f = eval_point(&tape, &Ctx::new(t, vec![]), p);
                if f.is_nan() {
                    continue;
                }
                let tol = 1e-9 * (1.0 + f.abs());
                assert!(b.lo - tol <= f && f <= b.hi + tol, "{src}: f({p:?}, t {t}) = {f} outside [{}, {}] over {region:?}", b.lo, b.hi);
            }
        }
    }
}

#[test]
fn extents_find_material_anywhere_and_say_when_they_cannot() {
    // query::extent covers all of space with pieces reaching infinity.
    // Cases: a part far out, a box turned about three axes (corner cells at
    // infinity),
    // a turn that varies with t composed with another (exact zeros widened
    // to +-1 ulp, times infinity), a rod (x * x not known to be a square),
    // and a mix of shapes (cancels far out; std blend does not).
    let window = ink::lang::DEFAULT_DOMAIN;
    let extent = |src: &str| {
        let mut d = doc(&format!("(body b (shape {src}))\n(assembly q (instance i b (at (rot-z 17) (rot-y (* 30 (+ 1 t))))))"));
        let ctx = d.ctx(&[], 0.0).unwrap();
        let i = d.items("q").unwrap().remove(0);
        let f = d.placed_field(&i);
        query::extent(&Arc::new(Tape::new(&d.g, f)), &ctx, window)
    };
    let bounded = |src: &str| matches!(extent(src), query::Extent::Bounded(..));
    match extent("(move (sphere 1) 20000000 0 0)") {
        query::Extent::Bounded(lo, hi) => assert!(lo[0] < 2e7 && hi[0] > 2e7, "{lo:?} {hi:?}"),
        e => panic!("{e:?}"),
    }
    for src in [
        "(rotate-z (rotate-y (rotate-x (box 10 20 30) 30) 40) 50)",
        "(rotate-x (cylinder-z 30 20) -90)",
        "(- (segment-distance-at here (vec 0 0 7.5) (vec 0 10 7.5)) 2.3)",
        "(blend (sphere 5) (box 8 8 8) 0.3)",
        "(smooth-union (sphere 5) (move (box 8 8 8) 6 0 0) 2)",
        "(max (twist-z (box 8 8 24) 3) (cylinder-z 7 26))",
        "(max x (sphere 1e5))",
        // A wedge in a slab, turned 180 degrees: the turn's 1e-16 cross
        // terms tie its coordinates, and abs of a one-signed form is
        // linear, so far out they stay tied (the printer's Z yokes).
        "(max (rotate-x (max (- (abs y) (+ 4 (* 0.44 z))) (slab-z -2.5 0)) 180) (- (abs (- x 3)) 7.5))",
    ] {
        assert!(bounded(src), "{src}: {:?}", extent(src));
    }
    for src in ["x", "(gyroid 10)", "(twist-z (box 8 8 24) 3)"] {
        assert_eq!(extent(src), query::Extent::Unbounded, "{src}");
    }
    // A union of a part only the boxes settle (a product cut by slabs: far
    // out along x near z = 0 no sector parts z x from 0) and one only the
    // sectors settle (a box turned about three axes): each member by its
    // own cover; neither cover alone settles the union.
    let d = doc("(def u (min (max (slab-z 1 2) (- 1 y) (- y 2) (- (abs (* z x)) 5)) (move (rotate-z (rotate-y (rotate-x (box 10 20 30) 30) 40) 50) 40 0 0)))");
    let Some(ink::lang::Value::Field(u, _)) = d.def("u") else { panic!() };
    let ctx = d.ctx(&[], 0.0).unwrap();
    match query::extent(&Arc::new(Tape::new(&d.g, *u)), &ctx, window) {
        query::Extent::Bounded(lo, hi) => assert!(lo[0] <= -5.0 && hi[0] >= 5.0, "{lo:?} {hi:?}"),
        e => panic!("{e:?}"),
    }
    // A union that shares its members (min of f and a tagged f, 40 deep):
    // each member is taken once, not 2^40 times.
    let d = doc("(fn dup (f i) (min f (tag alias f)))\n(def u (fold dup (sphere 1) (range 40)))");
    let Some(ink::lang::Value::Field(u, _)) = d.def("u") else { panic!() };
    let ctx = d.ctx(&[], 0.0).unwrap();
    assert!(matches!(query::extent(&Arc::new(Tape::new(&d.g, *u)), &ctx, window), query::Extent::Bounded(..)));
    // The cover starts from the window itself, however small (a radius
    // floor would leave a gap between the window and the sectors).
    let d = doc("(def ball (move (sphere 1e-8) 5e-7 0 0))");
    let Some(ink::lang::Value::Field(ball, _)) = d.def("ball") else { panic!() };
    let ctx = d.ctx(&[], 0.0).unwrap();
    for window in [([-1e-9; 3], [1e-9; 3]), ([-1e-9, -1.0, -1.0], [1e-9, 1.0, 1.0])] {
        match query::extent(&Arc::new(Tape::new(&d.g, *ball)), &ctx, window) {
            query::Extent::Bounded(lo, hi) => assert!(lo[0] <= 4.9e-7 && hi[0] >= 5.1e-7, "{lo:?} {hi:?}"),
            e => panic!("{window:?}: {e:?}"),
        }
    }
}

#[test]
fn an_intersections_extent_is_its_cheapest_decided_members() {
    // A max's extent is the first member's its search decides, its members
    // searched cheapest first (by the instructions each reads): a row of 24
    // balls in the window (bounded: the window) cut by a ball grown
    // negative (empty) and the half-space x < 0 (cheaper than both,
    // unbounded: its search fails) is the empty ball's, empty. In tape
    // order, the hash order of max's arguments, the row came first for
    // some of the balls.
    let window = ink::lang::DEFAULT_DOMAIN;
    let row: Vec<String> = (0..24).map(|k| format!("(move (sphere 1) {} 0 0)", 10 * k)).collect();
    for r in 1..=8 {
        let d = doc(&format!("(def f (max (min {}) (+ (sphere 1) {r}) x))", row.join(" ")));
        let Some(ink::lang::Value::Field(f, _)) = d.def("f") else { panic!() };
        let ctx = d.ctx(&[], 0.0).unwrap();
        assert_eq!(query::extent(&Arc::new(Tape::new(&d.g, *f)), &ctx, window), query::Extent::Empty, "ball grown {r}");
    }
}

// ---------------------------------------------------------------- ranged params

/// p p - (2 p - 1) with p param 0 (or `p` given), as a tape: 2^-54 at
/// p = 1 + 2^-27 in real numbers, 0 in floats (p p rounds to 2 p - 1).
fn rounding_trap(g: &mut ink::graph::Graph, p: ink::graph::NodeId) -> ink::graph::NodeId {
    let pp = g.mul(p, p);
    let two = g.constant(2.0);
    let one = g.constant(1.0);
    let tp = g.mul(two, p);
    let lin = g.sub(tp, one);
    g.sub(pp, lin)
}

#[test]
fn ranged_params_bound_over_their_spans() {
    use ink::graph::{Graph, Var};
    let mut g = Graph::new();
    let p = g.var(Var::Param(0));
    let pp = g.mul(p, p);
    let tape = Tape::new(&g, pp);
    let b = [I::point(0.0); 3];
    // At its value: the point 1, as a def would be.
    assert_eq!(eval_interval(&tape, &Ctx::new(0.0, vec![1.0]), b, None), I::point(1.0));
    // Ranged over [1, 3]: p^2 reaches 9. The constants-and-params fast path
    // would give the point 1, and (value (* p p) < 4) would be proven.
    let r = Ctx { p_spans: ParamSpans::new([(0, (1.0, 3.0))]), ..Ctx::new(0.0, vec![1.0]) };
    let v = eval_interval(&tape, &r, b, None);
    assert!(v.lo <= 1.0 && v.hi >= 9.0, "{v:?}");
    // Gradient bounds: the value over the span, and no slope in space.
    let gv = eval_interval_grad(&tape, &r, b);
    assert!(gv.v.lo <= 1.0 && gv.v.hi >= 9.0 && gv.d.iter().all(|d| *d == I::point(0.0)), "{gv:?}");
}

#[test]
fn a_ranged_param_narrowed_to_a_value_keeps_its_real_meaning() {
    // Unranged, a param expression is its f64 value (0 here, a point).
    // Ranged, it means its real value, even narrowed to a single value
    // (where a search confirms its witness): a proof over the span and its
    // witness are about the same expression.
    use ink::graph::{Graph, Var};
    let (q, tiny) = (1.0 + 2f64.powi(-27), 2f64.powi(-54));
    let mut g = Graph::new();
    let p = g.var(Var::Param(0));
    let f = rounding_trap(&mut g, p);
    let tape = Tape::new(&g, f);
    let b = [I::point(0.0); 3];
    assert_eq!(eval_interval(&tape, &Ctx::new(0.0, vec![q]), b, None), I::point(0.0));
    let one = Ctx { p_spans: ParamSpans::new([(0, (q, q))]), ..Ctx::new(0.0, vec![q]) };
    assert!(eval_interval(&tape, &one, b, None).contains(tiny));
    // A point of a wider span (Ctx::point) stays ranged.
    let wide = Ctx { p_spans: ParamSpans::new([(0, (1.0, 2.0))]), ..Ctx::new(0.0, vec![q]) };
    assert!(eval_interval(&tape, &wide.point(), b, None).contains(tiny));
    assert!(eval_interval_grad(&tape, &wide.point(), b).v.contains(tiny));
}

#[test]
fn a_branch_picked_over_a_span_keeps_reading_its_param() {
    // f = s s - (2 s - 1), s = select(p, q, 2): p ranged over [-2, -1]
    // picks q, and simplifying leaves arithmetic on q alone, an unranged
    // param. It still reads p in the model, so it keeps the real meaning it
    // has over the whole span (2^-54 at q = 1 + 2^-27), in the simplified
    // tape and its cones; recomputed from the simplified tape it would be
    // a float (0).
    use ink::graph::{Graph, Var};
    use ink::tape::BOTH;
    let (q, tiny) = (1.0 + 2f64.powi(-27), 2f64.powi(-54));
    let mut g = Graph::new();
    let (p, qv, two) = (g.var(Var::Param(0)), g.var(Var::Param(1)), g.constant(2.0));
    let s = g.select(p, qv, two);
    let f = rounding_trap(&mut g, s);
    let tape = Tape::new(&g, f);
    let ctx = Ctx { p_spans: ParamSpans::new([(0, (-2.0, -1.0))]), ..Ctx::new(0.0, vec![-1.5, q]) };
    let b = [I::point(0.0); 3];
    let mut ch = vec![BOTH; tape.len()];
    assert!(eval_interval(&tape, &ctx, b, Some(&mut ch)).contains(tiny));
    let sub = tape.simplify(&ch);
    assert!(sub.len() < tape.len());
    for t in [&sub, &sub.cone(sub.len() - 1)] {
        let v = eval_interval(t, &ctx, b, None);
        assert!(v.contains(tiny), "{v:?}");
    }
}

/// u = 2 p exactly: a derivation as a solve's would be, kept exact so the
/// test can say where spans must be.
struct Twice;
impl Derive for Twice {
    fn owns(&self) -> &[u32] {
        &[1]
    }
    fn reads(&self) -> &[u32] {
        &[0]
    }
    fn derive(&self, ctx: &mut Ctx) {
        let (lo, hi) = ctx.p_spans.get(0).unwrap();
        ctx.p_spans.set_span(1, (2.0 * lo, 2.0 * hi));
        ctx.params[1] = 2.0 * ctx.params[0];
    }
    fn rates(&self, _: &Ctx) -> Vec<Vec<f64>> {
        vec![vec![2.0]]
    }
}

#[test]
fn derived_params_follow_the_spans_they_are_derived_from() {
    // A derived param (a solve's unknown) is ranged, never halved: every
    // piece or point of the scope a search makes derives it again there, so
    // a witness at p is about u's value at p, not any value of its band.
    use ink::graph::{Graph, Var};
    let mut ctx = Ctx::new(0.0, vec![0.5, 1.0]);
    ctx.p_spans = ParamSpans::new([(0, (0.0, 1.0))]).derived(Arc::new(Twice), &[(1, (0.0, 2.0))]);
    let lo = ctx.with_param(0, (0.0, 0.5));
    assert_eq!((lo.p_spans.get(1), lo.params[1]), (Some((0.0, 1.0)), 0.5));
    let mid = ctx.param_mids();
    assert_eq!(mid.p_spans.get(1), Some((1.0, 1.0)));
    let mut w = ctx.clone();
    w.params[0] = 0.8;
    assert_eq!(w.point().p_spans.get(1), Some((1.6, 1.6)));
    assert_eq!(w.into_scope().params[1], 1.6);
    // A point outside the spans starts from the derived band.
    w.params[0] = 3.0;
    assert_eq!(w.point().p_spans.get(1), Some((6.0, 6.0)));
    // Witnesses name p; searches halve p for u.
    assert_eq!(ctx.p_spans.values(&ctx.params), vec![(0, 0.5)]);
    assert_eq!(ctx.p_spans.sources(1), vec![0]);
    let mut g = Graph::new();
    let (x, u) = (g.var(Var::X), g.var(Var::Param(1)));
    let f = g.sub(x, u);
    let tape = Tape::new(&g, f);
    assert_eq!(ink::scope::read_params(&tape, &ctx), vec![0]);
    let v = eval_interval(&tape, &ctx, [I::point(0.0); 3], None);
    assert!(v.lo <= -2.0 && v.hi >= 0.0, "{v:?}");
    // Rates along u bound the derivative (-1) over the scope; a push
    // toward a witness follows u's share to p (the chain rule).
    let (_, d) = eval_interval_rates(&tape, &ctx, [I::new(0.0, 1.0), I::point(0.0), I::point(0.0)], &Inputs { t: false, params: vec![1] });
    assert!(d[0].contains(-1.0) && d[0].hi - d[0].lo < 1e-12, "{d:?}");
    let mut grad = vec![0.0, -1.0];
    ctx.p_spans.chain(&ctx, &mut grad);
    assert_eq!(grad[0], -2.0);
}

#[test]
fn centred_lower_bounds_hold_through_shifts_negations_and_choices() {
    // The first-order bound sees the terms of a max through constant shifts
    // (an offset), negations (a union subtracted: -min = max of negations)
    // and choices (a min among the terms, or an intersection subtracted: a
    // pocket's floor and wall); several fields at once, each lowered by its
    // own amount (a clearance's offset). On finite and unbounded boxes
    // (over which only exactly zero slopes add nothing), the bound must not
    // pass the field's value at any point of the box.
    let ctx = Ctx::new(0.0, vec![]);
    let fields = [
        "(- (max (sphere 6) (- 2 z)) 0.3)",
        "(max (box 10 10 10) (neg (min (cylinder-z 2 20) (move (sphere 3) 4 0 0))))",
        "(max (cylinder-z 8 6) (neg (max (- (len2 x (- z (clamp z -2 2))) 2) (- 1 y))))",
        "(min (max (sphere 5) (- z 1)) (box 3 3 12))",
        "(+ 0.7 (neg (min (max (- x 1) (- y 2)) (max (- 0 x) (- z 3)))))",
        "(max (- 3 (abs y)) (- (+ (* 0.5 x) (* 0.8660254037844386 z)) 4) (- 1 (+ (* 0.5 x) (* 0.8660254037844386 z))))",
        "(max (slab-z 0 1) (- (slab-z 1.5 3) 0.4))",
    ];
    let mut rng = Lcg(31);
    let boxes = |rng: &mut Lcg| -> [I; 3] {
        let c = [rng.range(-12.0, 12.0), rng.range(-12.0, 12.0), rng.range(-12.0, 12.0)];
        let s = rng.range(0.01, 6.0);
        std::array::from_fn(|k| match (rng.next() < 0.15, rng.next() < 0.5) {
            (true, true) => I::new(c[k] - s, f64::INFINITY),
            (true, false) => I::new(f64::NEG_INFINITY, c[k] + s),
            _ => I::new(c[k] - s, c[k] + s),
        })
    };
    // Points of the box: within 40 of its finite ends.
    let point = |rng: &mut Lcg, b: &[I; 3]| -> [f64; 3] {
        std::array::from_fn(|k| match (b[k].lo.is_finite(), b[k].hi.is_finite()) {
            (true, true) => rng.range(b[k].lo, b[k].hi),
            (true, false) => rng.range(b[k].lo, b[k].lo + 40.0),
            (false, true) => rng.range(b[k].hi - 40.0, b[k].hi),
            _ => rng.range(-40.0, 40.0),
        })
    };
    for src in fields {
        let (_d, tape) = tape_of(src);
        for _ in 0..300 {
            let b = boxes(&mut rng);
            let Some(lo) = centred_lower(&tape, &ctx, b) else { continue };
            for _ in 0..48 {
                let p = point(&mut rng, &b);
                let f = eval_point(&tape, &ctx, p);
                assert!(lo <= f + 1e-9 * (1.0 + f.abs()), "{src}: lower bound {lo} above f({p:?}) = {f} over {b:?}");
            }
        }
    }
    // Two fields at once: max(a, b - s).
    for (a, b, s) in
        [("(sphere 5)", "(neg (sphere 5.5))", 0.4), ("(slab-z 0 1)", "(slab-z 1.5 3)", 0.45), ("(box 6 6 6)", "(neg (min (move (box 6 6 6) 0.2 0 0) (cylinder-z 1 9)))", 0.1)]
    {
        let ((_da, ta), (_db, tb)) = (tape_of(a), tape_of(b));
        for _ in 0..300 {
            let bx = boxes(&mut rng);
            let Some(lo) = centred_lower_max(&[(&ta, 0.0), (&tb, s)], &ctx, bx) else { continue };
            for _ in 0..48 {
                let p = point(&mut rng, &bx);
                let f = eval_point(&ta, &ctx, p).max(eval_point(&tb, &ctx, p) - s);
                assert!(lo <= f + 1e-9 * (1.0 + f.abs()), "max({a}, {b} - {s}): lower bound {lo} above {f} at {p:?} over {bx:?}");
            }
        }
    }
}

#[test]
fn intersections_of_unbounded_shapes_are_proven_empty_everywhere() {
    // Two parallel slabs, turned or not, never meet: their
    // max's first-order bound over all of space is the gap's half, as the
    // slopes cancel exactly. Shapes that do meet, however thinly, are not
    // called empty.
    let ctx = Ctx::new(0.0, vec![]);
    let window = ([-500.0; 3], [500.0; 3]);
    let turned = |s: &str| format!("(rotate-x (rotate-z {s} 25) 30)");
    let cases = [
        ("(max (slab-z 43.99 50.01) (- (slab-z 63.99 70.01) 13.9))".to_string(), true),
        (format!("(max {} {})", turned("(slab-z 0 1)"), turned("(slab-z 1.5 3)")), true),
        ("(max (slab-z 0 1) (slab-z 1.5 3) (- 2 x))".to_string(), true),
        ("(max (slab-z 43.99 50.01) (- (slab-z 63.99 70.01) 13.981))".to_string(), false),
        ("(max (slab-z 0 1) (slab-z 0.9 3))".to_string(), false),
        ("(max (slab-z 0 1) (max (neg x) (- x 1)))".to_string(), false),
    ];
    for (src, empty) in cases {
        let (_d, tape) = tape_of(&src);
        let e = query::extent(&Arc::new(tape), &ctx, window);
        assert_eq!(e == query::Extent::Empty, empty, "{src}: {e:?}");
    }
}

#[test]
fn first_order_bounds_refuse_more_choices_than_they_can_count() {
    // A max of 16 mins of 16 members each, x + 1000 + 16 i + j: x + 1240 in
    // all. Its 16^16 = 2^64 ways through the choices would wrap a count to
    // 0, the enumeration would check none and give its starting +inf, and
    // the shape would be called empty over all of space. Too many ways, and
    // the bound takes the mins as terms: it never passes the field over the
    // box.
    let ctx = Ctx::new(0.0, vec![]);
    let (_d, tape) = tape_of("(reduce max (map (lambda (i) (reduce min (map (lambda (j) (+ x (+ 1000 (* 16 i) j))) (range 16)))) (range 16)))");
    let yz = I::new(-1.0, 1.0);
    for b in [[I::new(-2010.0, -1990.0), yz, yz], [I::new(-1250.0, -1230.0), yz, yz], [I::new(f64::NEG_INFINITY, 0.0), yz, yz], [I::ENTIRE; 3]] {
        let least = b[0].lo + 1240.0;
        let lo = centred_lower(&tape, &ctx, b).unwrap_or(f64::NEG_INFINITY);
        assert!(lo <= least, "{b:?}: {lo} above the field's least, {least}");
        assert!(!centred_at_least(&[(&tape, 0.0)], &ctx, b, (least + 1.0).max(0.0)), "{b:?}");
    }
    assert_ne!(query::extent(&Arc::new(tape), &ctx, ([-500.0; 3], [500.0; 3])), query::Extent::Empty);
    // Terms that share subexpressions are walked once for each way to them:
    // max(f, f + 1) nested 40 deep, shifted (so its root alone is one term),
    // has 2^40 ways through its shifts, more than a walk can take. It stops
    // past the terms a bound tries.
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let (_d, tape) = tape_of("(- (fold (lambda (a i) (max a (+ a 1))) x (range 40)) 0.5)");
        let b = [I::new(-50.0, -45.0), yz, yz];
        tx.send(centred_lower(&tape, &ctx, b).is_none_or(|lo| lo <= -50.0 + 39.5)).ok();
    });
    assert_eq!(rx.recv_timeout(std::time::Duration::from_secs(30)), Ok(true));
}

#[test]
fn clearance_pairs_prove_thin_margins_with_coarse_cells() {
    // A pin of radius 5 in a bore of 5.5, cleared by 0.49: a margin of 0.01.
    // Tests of the grown shape alone need cells under the margin; the
    // first-order bound of the pair (the pin's surface facing the bore's)
    // needs about sqrt(margin x radius).
    let d = doc("(def plate (subtract (box 40 40 10) (cylinder-z 5.5 100))) (def pin (cylinder-z 5 30))");
    let t = |n: &str| {
        let ink::lang::Value::Field(f, _) = d.def(n).unwrap().clone() else { panic!() };
        Arc::new(Tape::new(&d.g, f))
    };
    let ctx = Ctx::new(0.0, vec![]);
    let (lo, hi) = ([-6.0, -6.0, -6.0], [6.0, 6.0, 6.0]);
    assert!(matches!(query::clearance(&[query::Way { grown: t("plate"), covered: t("pin"), lo, hi, swapped: false }], &ctx, 0.49, 0.1, query::BOXES).0, Verdict::Proven));
    // The witness is the true gap, 0.5.
    match query::clearance(&[query::Way { grown: t("plate"), covered: t("pin"), lo, hi, swapped: false }], &ctx, 0.505, 0.1, query::BOXES).0 {
        Verdict::Disproven { value, .. } => assert!((0.5 - 1e-9..0.505).contains(&value), "{value}"),
        v => panic!("0.505 should be disproven: {v:?}"),
    }
}

/// A jump a box straddles once is bounded as a choice of its two sides
/// (eval::forked): every point's value lies in the bound, and the tape
/// pruned by the choices made over the box is the field there. Boxes about
/// the seams of polar repeats (atan2's cut on a sector edge, and inside
/// one; one nested in another), linear repeats, a floor pattern, a thread's
/// root and a path loop's station, and the far-sector evaluator's boxes.
#[test]
fn jumps_bound_as_choices_hold_the_field() {
    let cases = [
        "(repeat-polar (move (box 6 2 1) 20 0 0) 9)",
        "(repeat-polar (rotate-x (move (box 6 2 1) 20 0 0) 30) 10)",
        "(repeat-polar (move (repeat-polar (move (box 2 1 1) 4 0 0) 6) 18 0 0) 8)",
        "(repeat-x (sphere 4) 10)",
        "(warp (box 9 0.7 12) (- x (* 18 (- (clamp (floor (+ (/ x 18) 1)) 0 1) 0.5))) y z)",
        "(repeat-x (repeat-y (repeat-z (sphere 3) 8) 8) 8)",
        "(subtract (cylinder-z 30 6) (thread-iso 20 3))",
        "(max (involute-gear-2d 15 4 20 0.03) (slab-z -3 3))",
    ];
    let mut rng = 12345u64;
    let mut u = || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        (rng >> 11) as f64 / (1u64 << 53) as f64
    };
    for src in cases {
        let (d, tape) = tape_of(src);
        let ctx = d.ctx(&[], 0.0).unwrap();
        for _ in 0..3000 {
            let c = [60.0 * u() - 30.0, 60.0 * u() - 30.0, 6.0 * u() - 3.0];
            let h = 0.01 + 2.0 * u();
            let b = c.map(|v| I::new(v - h, v + h));
            let mut ch = vec![ink::tape::BOTH; tape.len()];
            let bound = eval_interval(&tape, &ctx, b, Some(&mut ch));
            let linear = eval_interval_linear(&tape, &ctx, &Region::Box(b), None);
            let pruned = tape.simplify(&ch);
            for _ in 0..8 {
                let p = c.map(|v| v - h + 2.0 * h * u());
                let (a, s) = (eval_point(&tape, &ctx, p), eval_point(&pruned, &ctx, p));
                assert!(bound.contains(a) && linear.contains(a), "{src}: {a} outside {bound:?} / {linear:?} at {p:?}");
                assert!(a == s, "{src}: pruned {s} against {a} at {p:?}");
            }
        }
    }
}

/// A jump a box and scope straddle once is a choice of its two sides in
/// bounds first order in time and ranged params too (eval::scope_run), each
/// side its branch continued over the whole region: an arm of a polar
/// repeat on the seam (atan2's cut), turned by a param or over a span of
/// time, is bounded as that arm turning, and a linear repeat's seam moved
/// by a param as the copy beside it, not with the jump's rate unbounded
/// there and the bound the whole period's: the gaps between the arms
/// (above 0) and the copies' insides (below 0). Where the params do not
/// move the seam (an arm tilted about its own axis), each point takes one
/// side for every value of them, so its rates are the sides' (+-0.09; not
/// +-8.6, from the whole period's bounds). Every point of the box at every value of the scope lies in
/// the bounds, its rates in theirs, and its change from the scope's centre
/// within what bounded rates allow (the mean value bound motion's steps
/// read): so a value that jumps where the params move the cut (half an
/// angle's sine, -1 to 1) has unbounded rates. The centred bound (each
/// value of the tape kept, every value of a window joined) holds too.
#[test]
fn first_order_bounds_take_a_seam_as_a_choice() {
    let arms = "(repeat-polar (move (box 4 2 2) 10 0 0) 12)";
    let seam = [I::new(-10.05, -9.95), I::new(-0.05, 0.05), I::new(-0.05, 0.05)];
    let cases = [
        // (field, box, p's span, t's span, the bound's ends lie in, rates bounded)
        (format!("(rotate-z {arms} p)"), seam, Some((-1.0, 1.0)), None, I::new(-1.01, -0.7), false),
        (format!("(rotate-z {arms} (* 2 t))"), seam, None, Some((0.0, 0.5)), I::new(-1.01, -0.7), false),
        ("(move (repeat-x (sphere 2) 10) p 0 0)".into(), [I::new(4.95, 5.05), I::new(-0.05, 0.05), I::new(-0.05, 0.05)], Some((-0.5, 0.5)), None, I::new(2.4, 3.6), false),
        (
            "(repeat-polar (move (rotate-x (box 4 2 2) (* 40 p)) 10 0 0) 12)".into(),
            [I::new(-10.05, -9.95), I::new(-0.05, 0.05), I::new(0.35, 0.45)],
            Some((-0.25, 0.25)),
            None,
            I::new(-0.7, -0.5),
            true,
        ),
        ("(+ (* 0.01 x) (sin (* 0.5 (atan2 (- y (* 0.1 p)) x))))".into(), seam, Some((-1.0, 1.0)), None, I::new(-1.11, 1.11), false),
        ("(sin (mod (- x (* 3 p)) 10))".into(), seam, Some((-0.1, 0.1)), None, I::new(-0.85, 0.4), false),
        ("(max (sin (* 0.5 (atan2 (* 0.1 p) -1))) (- (* 0.01 x) 2))".into(), seam, Some((-1.0, 1.0)), None, I::new(-1.01, 1.01), false),
    ];
    let mut rng = Lcg(31);
    for (src, b, sp, st, within, bounded) in cases {
        let d = doc(&format!("(param p 0) (def f {src})"));
        let ink::lang::Value::Field(f, _) = d.def("f").unwrap().clone() else { panic!() };
        let tape = Tape::new(&d.g, f);
        let kp = d.param_index("p").unwrap() as u32;
        let mut params = vec![0.0; d.params.len()];
        params[kp as usize] = sp.map_or(0.0, span_mid);
        let ctx = Ctx { t_span: st, p_spans: ParamSpans::new(sp.map(|s| (kp, s))), ..Ctx::new(st.map_or(0.0, span_mid), params) };
        let inputs = Inputs { t: st.is_some(), ..Inputs::of([&tape], &ctx) };
        let (r, rates) = eval_interval_rates(&tape, &ctx, b, &inputs);
        let lower = centred_lower(&tape, &ctx, b).unwrap_or(f64::NEG_INFINITY);
        assert!(within.contains(r.lo) && within.contains(r.hi), "{src}: {r:?}, not within {within:?}");
        // Unbounded is ENTIRE, never NaN (which every test below would pass).
        assert!(rates.iter().all(|x| !x.lo.is_nan() && !x.hi.is_nan()) && !lower.is_nan(), "{src}: rates {rates:?}, centred {lower}");
        let finite = rates.iter().all(|x| x.lo.is_finite() && x.hi.is_finite());
        assert!(finite || !bounded, "{src}: rates {rates:?}");
        let centre = |pt: [f64; 3]| eval_point(&tape, &Ctx::new(st.map_or(0.0, span_mid), ctx.params.clone()), pt);
        for _ in 0..4000 {
            let pt = b.map(|x| rng.range(x.lo, x.hi));
            let t = st.map_or(0.0, |s| rng.range(s.0, s.1));
            let mut vals = ctx.params.clone();
            vals[kp as usize] = sp.map_or(0.0, |s| rng.range(s.0, s.1));
            let at = Ctx::new(t, vals);
            let v = eval_point(&tape, &at, pt);
            assert!(r.contains(v) && lower <= v, "{src}: f({pt:?}, {:?}, t {t}) = {v} outside {r:?} or below {lower}", at.params);
            let dual = eval_dual(&tape, &at, pt);
            let along: Vec<(usize, f64)> = st.map(|s| (4, t - span_mid(s))).into_iter().chain(sp.map(|s| (5 + kp as usize, at.params[kp as usize] - span_mid(s)))).collect();
            let mut moved = 0.0;
            for (j, &(x, du)) in along.iter().enumerate() {
                assert!(rates[j].contains(dual[x]) || rates[j] == I::ENTIRE, "{src}: rate {j} at {pt:?}, {:?}, t {t} = {} outside {:?}", at.params, dual[x], rates[j]);
                moved += ink::interval::mag(rates[j]) * du.abs();
            }
            let change = (v - centre(pt)).abs();
            assert!(!finite || change <= moved * (1.0 + 1e-9) + 1e-12, "{src}: f({pt:?}, {:?}, t {t}) moved {change} from the centre, rates {rates:?} allow {moved}", at.params);
        }
    }
}

/// A part turned through a whole turn by a ranged param: over all of it
/// its far points bound as any direction, so no direction of space is
/// shown free of it (its extent unbounded, and a claim about it holds only
/// in the window); over pieces of the turn each is, and its reach is the
/// union of theirs (ink-lang's extent_in_pieces): bounded, and holding it.
#[test]
fn a_part_turned_through_a_ranged_whole_turn_has_a_bounded_reach() {
    let d = doc("(param u 0)\n(def f (rotate-z (move (cylinder-z 3 10) 15 0 0) (* 360 u)))");
    let ink::lang::Value::Field(f, _) = d.def("f").unwrap().clone() else { panic!() };
    let k = d.params.iter().position(|p| p.name == "u").unwrap() as u32;
    let ctx = Ctx { p_spans: ParamSpans::new([(k, (0.0, 1.0))]), ..Ctx::new(0.0, vec![0.0]) };
    let r = d.reach(f, &ctx);
    assert!(r.whole, "{r:?}");
    let e = ink::query::extent(&Arc::new(Tape::new(&d.g, f)), &ctx, (r.lo, r.hi));
    assert!(matches!(e, ink::query::Extent::Unbounded), "over the whole turn at once: {e:?}");
}

/// Outside ink's range (docs/CORE.md Answers) a claim is unknown, saying
/// why; its twin inside is decided.
/// Each pair is a false guarantee (a box at 1e16 whose rounded
/// midpoint is a face, so its centre reads as empty; a profile of 1e-200
/// whose radius underflows; a strip whose fixed coefficient underflows
/// to 0, proven symmetric) and the same question at engineering scale.
#[test]
fn answers_outside_the_range_of_validity_are_unknown() {
    use ink::claims::{ClaimOptions, Evidence, Status, run_asserts};
    use ink::view::ModelView;
    let ask = |src: &str| {
        let d = doc(src);
        let ctx = d.ctx(&[], 0.0).unwrap();
        let mut v = ModelView::new(d, ctx);
        run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts
    };
    let bead = |at: &str| format!("(def bead (max (- (abs (- (- x {at}) 1)) 0.75) (- (abs y) 2) (- (abs z) 2))) (assert e (empty bead) (accuracy 0.1))");
    let outside = |src: &str, what: &str| {
        let a = ask(src);
        let c = &a[0].check;
        assert!(c.status == Status::Unknown && c.accuracy.is_none(), "{src}: {c:?}");
        assert!(matches!(&c.evidence, Evidence::Outside(why) if why.contains(what)), "{src}: {:?}", c.evidence);
    };
    // Where the search looks: past 1,000 km.
    outside(&bead("1e16"), "a coordinate 1.7");
    assert_eq!(ask(&bead("100"))[0].check.status, Status::Fails);
    // A length finer than a nanometre, or a billionth of the farthest
    // coordinate.
    outside("(assert e (empty (move (sphere 1) 1e8 0 0)) (accuracy 0.01))", "finer than");
    outside("(assert e (empty (sphere 1)) (accuracy 1e-7))", "finer than");
    assert_eq!(ask("(assert e (empty (move (sphere 1) 1e6 0 0)) (accuracy 0.01))")[0].check.status, Status::Fails);
    // A point probed past the reach.
    outside("(assert p (subset (probe 1e12 0 0) (move (sphere 1) 1e12 0 0)))", "a coordinate 1e12 mm");
    // A number past 1e100 or under 1e-100 in size, whatever reads it (here
    // profile-reach, in Rust, at load).
    outside("(def tiny (max (- (abs x) 1e-200) (- (abs y) 1e-200) (- (abs z) 1e-200))) (def reach (profile-reach tiny)) (assert r (value reach < 1e-200))", "1e-200");
    let strip = |p: &str, a: &str, b: &str, c: &str| {
        format!(
            "(param p {p}) (def a {a}) (def b {b}) (def c {c})
            (def u (* (* x (* (* p a) b)) c)) (def v (* (* (* y (* p b)) a) c))
            (def strip (- (sqrt (+ (square u) (square v))) 0.01))
            (body rotor (shape (max strip (cylinder-z 2 1)))) (body post (shape (sphere 0.006)))
            (t-span 0 1)
            (assembly root (instance r rotor (at) (turn z (* 360 t))) (instance p post (at (move 0 1.5 0))))
            (assert no-collision (disjoint (shape-of root/r) (shape-of root/p)) (accuracy 0.001))"
        )
    };
    let in_motion = |src: &str| {
        let d = doc(src);
        let ctx = d.ctx(&[], 0.0).unwrap();
        let mut v = ModelView::new(d, ctx);
        v.set_motion(Some((0.0, 1.0))).unwrap();
        run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts.remove(0).check
    };
    let c = in_motion(&strip("1e-200", "1e-200", "1e300", "1e100"));
    assert!(c.status == Status::Unknown && matches!(&c.evidence, Evidence::Outside(why) if why.contains("1e-200")), "{c:?}");
    // At ordinary magnitudes a strip (u 0: |y| < 0.01) turns into the
    // post: it fails.
    assert_eq!(in_motion(&strip("1", "0", "1", "1")).status, Status::Fails);
}
