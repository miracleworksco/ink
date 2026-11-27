//! The topology package through the language: flanges, the MBB beam against the 88-line
//! code, keep, settling, the design returned, the form's checks.

use ink_ext as ink;
use std::path::Path;

/// This folder: what the tests load is found from it.
#[path = "cli/roots.rs"]
mod roots;
use roots::ROOT;

fn load(src: &str) -> ink::lang::Doc {
    ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"))
}

/// Topology `name` of `src`, at `res`, solved to 1e-9.
fn topology(src: &str, name: &str, res: Option<f64>) -> ink::topology::Layout {
    let mut d = load(src);
    let ctx = d.ctx(&[], 0.0).unwrap();
    let l = ink::topology::run(&mut d, &ctx, name, res, 1e-9).unwrap();
    assert!(l.sol.report.converged, "solver did not converge: {:?}", l.sol.report);
    l
}

/// A field a model defines, as a function of the point.
fn field(d: &ink::lang::Doc, name: &str) -> impl Fn([f64; 3]) -> f64 {
    let ink::lang::Value::Field(f, _) = d.def(name).unwrap().clone() else { panic!("{name} is no field") };
    let (t, ctx) = (ink::tape::Tape::new(&d.g, f), d.ctx(&[], 0.0).unwrap());
    move |p| ink::eval::eval_point(&t, &ctx, p)
}

#[test]
fn topology_optimization_finds_the_cantilever_flanges() {
    // Quasi-2D cantilever: material should go to the top and bottom fibres at
    // the wall (bending) and leave the neutral axis there.
    let src = "(def domain (move (box 60 2 24) 30 0 0))
(study beam (shape domain) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (- x 0.5))
  (load tip (force 0 0 -100N) (region (max (- 59 x) (- (abs z) 2)))) (res 1))
(topology layout (study beam) (volume-fraction 0.4) (iterations 15) (filter 2.5))
(def dens (study-field layout density))
(def sol (study-field layout solid))";
    let mut d = load(src);
    let ctx = d.ctx(&[], 0.0).unwrap();
    // Solved once, at load time (study-field); this is the cached result.
    let l = ink::topology::run(&mut d, &ctx, "layout", None, 1e-8).unwrap();
    let h = &l.record.history;
    let (first, last) = (h[0], *h.last().unwrap());
    assert!((last.2 - 0.4).abs() < 1e-3, "volume fraction {}", last.2);
    assert!(last.1 < 0.3 * first.1, "compliance {} from {}", last.1, first.1);
    let at = field(&d, "dens");
    assert!(at([5.0, 0.0, 11.0]) > 0.8 && at([5.0, 0.0, -11.0]) > 0.8, "flanges {} {}", at([5.0, 0.0, 11.0]), at([5.0, 0.0, -11.0]));
    assert!(at([5.0, 0.0, 0.0]) < 0.2, "neutral axis {}", at([5.0, 0.0, 0.0]));
    // `solid` is a shape's field (1/2 - density), negative where material
    // was kept: alone it is the optimized part, bounded as master proves it
    // (the domain and the grid's few cells past it). Kept material reaches
    // the domain's faces (the wall's fix, the load), so held past the grid
    // as a quantity is (ResultSpec::shape) it was unbounded.
    let ink::lang::Value::Field(sol, _) = d.def("sol").unwrap().clone() else { panic!() };
    let t = std::sync::Arc::new(ink::tape::Tape::new(&d.g, sol));
    let (lo, hi, _) = ink::query::bounds(&t, &ctx, [-500.0; 3], [500.0; 3]).expect("material");
    assert!(
        (0..3).all(|a| lo[a] >= [-4.0, -4.0, -16.0][a] && hi[a] <= [64.0, 4.0, 16.0][a]) && lo[0] < 2.0 && hi[2] > 10.0,
        "solid alone: {lo:?}..{hi:?}, the domain's 0..60, -1..1, -12..12"
    );
}

#[test]
fn topology_reports_the_held_fraction() {
    // `keep` holds 80% of the domain solid, above the 30% budget: the record
    // must say so (`ink topology` warns), since the free material then
    // vanishes.
    let src = "(def domain (move (box 40 4 10) 20 0 0))
(study beam (shape domain) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (- x 0.5))
  (load tip (force 0 0 -10N) (region (- 39.5 x))) (res 1))
(topology layout (study beam) (volume-fraction 0.3) (iterations 3) (keep (- (abs z) 4)))";
    let (held, vf) = topology(src, "layout", Some(1.0)).record.fixed;
    assert!(held > 0.78 && (vf - 0.3).abs() < 1e-12, "held {held}, target {vf}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn topology_matches_the_88_line_codes_mbb_beam() {
    // Andreassen et al. 2011's MBB half-beam (their 88-line code's example):
    // 60 x 20 cells, volume 0.5, filter radius 2.4, penalty 3, a unit load
    // at the top of the symmetry line, rolled at the far bottom corner. Here
    // a slab one cell thick: x held over the symmetry face's first quarter
    // cell, z over the corner's, the load on the top corner's. Its sideways
    // slide is a free motion held at rest: so held, the final solve once
    // stalled at a residual of 1.3e-5 (the hold swamped the void cells'
    // diagonal). Compliance over the solid domain's, which cancels the
    // elements (Q4 there, bricks here: solid 125.9 against 122.2) to first
    // order: the 88-line code ported to numpy (benchmarks/tools/top88_mbb.py)
    // gives 1.595 with its own Heaviside filter and 1.512 with this
    // projection and update (tanh at 1/2, beta 1 to 16, MMA); here 1.525. Within
    // 3%, the elements' and the supports' difference; the 88-line code's own
    // three filters span 1.60 to 1.86. And the design settles before its
    // budget: optimality criteria at beta 16 flipped cells by the whole
    // move limit to the last iteration.
    let src = "(def domain (move (box 60 1 20) 30 0.5 10))
(def sym (- x 0.25))
(def roller (max (- 59.75 x) (- z 0.25)))
(def top (max (- x 0.25) (- 19.75 z)))
(study solid (shape domain) (elastic (E 1MPa) (nu 0.3))
  (fix sym sym (dofs x)) (fix roller roller (dofs z)) (load top (force 0 0 -1N) (region top)) (res 1))
(topology layout (study solid) (volume-fraction 0.5) (iterations 100) (filter 2.4))";
    let layout = topology(src, "layout", Some(1.0));
    let mut d = load(src);
    let ctx = d.ctx(&[], 0.0).unwrap();
    let (solid, _, _) = ink::sim::run_study(&mut d, &ctx, "solid", Some(1.0), 1e-9).unwrap();
    let ratio = layout.sol.report.energy / solid.report.energy;
    assert!((ratio / 1.512 - 1.0).abs() < 0.03, "compliance {} over the solid's {}: {ratio}", layout.sol.report.energy, solid.report.energy);
    let last = *layout.record.history.last().unwrap();
    assert!(last.0 < 100 && last.3 < 0.01, "not settled: iteration {}, change {}", last.0, last.3);
    assert!((last.2 - 0.5).abs() < 1e-3, "volume fraction {}", last.2);
}

#[test]
fn topology_options_are_checked() {
    let study = "(import \"std\")\n(study s (shape (box 10 1 10)) (elastic (E 1MPa) (nu 0.3)) (fix w (+ x 4.5)) (load l (force 0 0 -1N) (region (- x 4.5))))\n";
    let moved = "(import \"std\")\n(study s (shape (box 20 4 4)) (elastic (E 1000MPa) (nu 0.3)) (fix wall (+ x 9.5)) (fix end (- 9.5 x) (dofs z) (move 0 0 -0.1)))\n";
    let modal = "(import \"std\")\n(study s (shape (box 20 4 4)) (elastic (E 1000MPa) (nu 0.3) (density 1000kg/m3)) (fix wall (+ x 9.5)) (modes 2))\n";
    for (src, msg) in [
        (format!("{study}(topology t (study s) (volume-fraction 1.5))"), "between 0 and 1"),
        (format!("{study}(topology t (study s) (volume-fraction 0.5) (graded 1.5))"), "unknown topology option 'graded'"),
        (format!("{study}(topology t (volume-fraction 0.5))"), "topology needs (study s)"),
        (format!("{study}(topology t (study nope))"), "no study named 'nope'"),
        (format!("{moved}(topology t (study s) (volume-fraction 0.5) (iterations 2))"), "no (move ...)"),
        (format!("{modal}(topology t (study s))"), "a topology lays out a static elastic study"),
        // Refused, not clamped to one cell or one update.
        (format!("{study}(topology t (study s) (filter -1))"), "filter must be a positive length"),
        (format!("{study}(topology t (study s) (filter 0))"), "filter must be a positive length"),
        (format!("{study}(topology t (study s) (iterations 0))"), "iterations must be a whole number, at least 1"),
        (format!("{study}(topology t (study s) (iterations -3))"), "iterations must be a whole number, at least 1"),
        (format!("{study}(topology t (study s) (iterations 2.5))"), "iterations must be a whole number, at least 1"),
    ] {
        let e = ink::load_str(&src, "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{src} loaded"));
        assert!(e.message.contains(msg) || e.hint.as_deref().is_some_and(|h| h.contains(msg)), "{src}: {e:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn topology_holds_no_material_under_supports() {
    // Only `keep` holds material solid (passive elements, as Sigmund's and
    // Andreassen et al.'s codes have them): a plane held along its normal
    // through every cell of a one-cell plate held them all solid, leaving
    // nothing to optimize. The optimizer empties cells under it.
    let src = "(def domain (move (box 40 1 10) 20 0.5 0))
(study plate (shape domain) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (- x 0.5)) (fix plane (- y 0.25) (dofs y))
  (load tip (force 0 0 -10N) (region (- 39.5 x))) (res 1))
(topology layout (study plate) (volume-fraction 0.3) (iterations 20) (filter 1.5))
(def dens (study-field layout density))";
    let mut d = load(src);
    let ctx = d.ctx(&[], 0.0).unwrap();
    let l = ink::topology::run(&mut d, &ctx, "layout", None, 1e-8).unwrap();
    assert_eq!(l.record.fixed.0, 0.0);
    let at = field(&d, "dens");
    let void = (0..40).filter(|&i| at([i as f64 + 0.5, 0.5, 4.5]) < 0.1).count();
    assert!(void > 5, "{void} of 40 cells along z = 4.5 emptied");
}

#[test]
fn topology_settles_only_within_its_volume() {
    // One cell, 0.99 of it to keep: it stopped as settled at iteration 6,
    // a steepness just doubled, with 0.9999 kept (the update meets the
    // volume at the steepness it ran at, not the next).
    let src = "(study s (shape (move (box 1 1 1) 0.5 0.5 0.5)) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 0.25)) (load tip (force 1N 0 0) (region (- 0.75 x))))
(topology t (study s) (volume-fraction 0.99) (iterations 100) (filter 1))";
    let r = topology(src, "t", Some(1.0)).record;
    let last = *r.history.last().unwrap();
    assert!(r.settled && last.2 <= 0.99 + 1e-3, "settled {} with {}", r.settled, last.2);
}

#[test]
fn a_topology_out_of_iterations_describes_the_design_it_returns() {
    // An MBB half-beam given 3 iterations: the summary (the history's last
    // row) is the returned design's, solved after the last update, not the
    // one before it (507.3 N mm and 0.4993 reported for 469.8 and 0.4917).
    let src = "(def domain (move (box 12 1 4) 6 0.5 2))
(study s (shape domain) (elastic (E 1MPa) (nu 0.3))
  (fix sym (- x 0.25) (dofs x)) (fix roller (max (- 11.75 x) (- z 0.25)) (dofs z))
  (load top (force 0 0 -1N) (region (max (- x 0.25) (- 3.75 z)))))
(topology t (study s) (volume-fraction 0.5) (iterations 3) (filter 1.5))";
    let l = topology(src, "t", Some(1.0));
    let last = *l.record.history.last().unwrap();
    assert!(!l.record.settled && (last.1 / l.sol.report.energy - 1.0).abs() < 1e-9, "{last:?} against {}", l.sol.report.energy);
    let gf = l.sol.grid_fields();
    let (dens, frac) = (l.sol.cells_on_grid(&l.record.density), &gf.cell_frac);
    let kept = dens.iter().zip(frac).filter(|(_, f)| **f > 0.0).map(|(d, f)| *d as f64 * *f as f64).sum::<f64>() / frac.iter().map(|f| *f as f64).sum::<f64>();
    assert!((kept - last.2).abs() < 1e-5 && (last.2 - 0.5).abs() < 1e-6, "returned {kept}, reported {}", last.2);
}

#[test]
fn topology_solve_is_independent_of_thread_count() {
    // The complete serialized layout (displacements, stresses, densities,
    // history), not just a headline value.
    let src = "(body b (shape (box 20 2 8)))
(study beam (shape (shape-of b)) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (+ x 9.5)) (load tip (force 0 0 -1N) (region (- 9 x))))
(topology layout (study beam) (volume-fraction 0.4) (iterations 3) (filter 1.5))";
    let run = |threads| {
        rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap().install(|| {
            let mut d = load(src);
            let ctx = d.ctx(&[], 0.0).unwrap();
            let (setup, _, _) = ink::sim::study_setup(&mut d, &ctx, "beam", Some(1.0), 1e-9).unwrap();
            let top = ink::topology::Settings { volume_fraction: 0.4, iterations: 3, filter: Some(1.5), keep: None };
            let mut record = None;
            let mut s = ink::sim::solve_redesigned(&setup, &ctx, &mut |a| {
                record = Some(ink::topology::optimize(a, &top, &ctx, setup.tol)?);
                Ok(())
            })
            .unwrap();
            assert!(s.report.converged && s.report.failure.is_none(), "{:?}", s.report);
            (s.report.setup_ms, s.report.solve_ms) = (0.0, 0.0);
            serde_json::to_string(&(s, record.unwrap())).unwrap()
        })
    };
    let baseline = run(1);
    for threads in [2, 4] {
        assert!(run(threads) == baseline, "layout bits differ with {threads} threads");
    }
}
