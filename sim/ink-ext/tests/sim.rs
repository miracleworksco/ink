//! Simulation validated against closed-form solutions.

mod common;

use ink::eval::{Ctx, eval_points};
use ink::lang::{Doc, Value};
use ink::sim;
use ink::tape::Tape;
use ink_ext as ink;
use std::path::Path;

/// This folder: what the tests load is found from it.
#[path = "cli/roots.rs"]
mod roots;
use roots::ROOT;

fn solve(src: &str, study: &str, res: f64) -> sim::Solution {
    let mut d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let (sol, _, _) = sim::run_study(&mut d, &ctx, study, Some(res), 1e-9).unwrap();
    assert!(sol.report.converged, "solver did not converge: {:?}", sol.report);
    sol
}

/// Sum of applied loads and reactions, relative to the largest term.
fn imbalance(s: &sim::Solution) -> f64 {
    let r = &s.report;
    let mut sum = [0.0; 3];
    let mut scale = 0.0f64;
    for (_, v) in r.applied.iter().chain(&r.reactions) {
        for d in 0..3 {
            sum[d] += v[d];
            scale = scale.max(v[d].abs());
        }
    }
    sum.iter().map(|v| v.abs()).fold(0.0, f64::max) / scale
}

const BAR: &str = "(body bar (shape (box 100 10 10)))
(study pull (shape (shape-of bar)) (elastic (E 1000MPa) (nu 0))
  (fix wall (+ x 49.5)) (load pull (force 1000N 0 0) (region (- 49 x))))";

#[test]
fn bar_in_tension() {
    // sigma = F/A = 10 MPa; u(x) = 0.01 (x + 49.75) from the clamp's mid-plane.
    let s = solve(BAR, "pull", 1.0);
    let (u, (_, vm, _)) = (s.value_at([40.0, 0.0, 0.0]).unwrap(), s.field_at([0.0, 0.0, 0.0]).unwrap());
    assert!((vm - 10.0).abs() < 0.1, "stress {vm}");
    assert!((u[0] - 0.8975).abs() < 0.015, "u(40) = {}", u[0]);
    assert!(u[1].abs() < 0.01 && u[2].abs() < 0.01, "spurious bending {u:?}");
    assert!(imbalance(&s) < 1e-6);
}

/// A 20 x 2 x 2 bar of 6061 turned 30 degrees about z off the grid, held
/// over u < 0.5 (u along it), and 10 N down at its tip, u > 19.5 (or, with
/// a density, its modes).
fn turned_cantilever(modes: bool) -> String {
    let u = "(+ (* x 0.8660254037844387) (* y 0.5))";
    let rest = if modes { "(modes 2)".to_string() } else { format!("(load tip (force 0 0 -10N) (region (- 19.5 {u})))") };
    format!(
        "(body b (shape (rotate-z (move (box 20 2 2) 10 0 0) 30)))
(study s (shape (shape-of b)) (elastic (E 68900MPa) (nu 0.33) (density 2700kg/m3))
  (fix wall (- {u} 0.5)) {rest})"
    )
}

/// A 20 x 2 x 2 bar of 6061, 10 N down at its tip, held where `wall` < 0.
fn held_bar(wall: &str) -> String {
    format!(
        "(study s (shape (move (box 20 2 2) 10 0 0)) (elastic (E 68900MPa) (nu 0.33))
  (fix wall {wall}) (load tip (force 0 0 -10N) (region (- 19.5 x))))"
    )
}

#[test]
fn the_solvers_precision_is_a_floor_only_for_what_it_is_derived_for() {
    // A displacement or a temperature at a node reads the field's span of
    // the solver's precision; a frequency and a compliance read their
    // solve's residual alone (a mode's span is its amplitude's, and a
    // compliance's floor came out 1e-6 over the load in newtons).
    let bar = held_bar("(- x 0.998)").replace("(elastic (E 68900MPa) (nu 0.33))", "(elastic (E 68900MPa) (nu 0.33) (density 2700kg/m3))");
    let precision = |claim: &str| format!("{:?}", checks(&format!("{bar}\n(require c {claim})")).remove(0).1).contains("solver's precision");
    assert!(precision("(max-displacement s < 1mm)"));
    assert!(!precision("(compliance s < 100N*mm)"));
    let modal = bar.replace("(load tip (force 0 0 -10N) (region (- 19.5 x)))", "(modes 1)");
    assert!(!format!("{:?}", checks(&format!("{modal}\n(require f (min-frequency s > 1Hz))")).remove(0).1).contains("solver's precision"));
}

#[test]
fn a_block_hold_reads_alike_either_side_of_a_cell_plane() {
    // A 20 x 4 x 4 bar (nu 0, pulled 10 MPa) held by a block x < X, y and z <
    // 0.45, its face x = X a millionth either side of the cell plane x = 1.2
    // on cells of 0.1 mm: the block's side faces lie on sample faces at 0.45,
    // a size binary fractions do not hold, and a strip one sample wide along
    // its edge was not traced below the plane (0.33543 against 0.33709
    // above; 0.3351665 and 0.3351661 with it traced).
    let tip = |x: f64| {
        let src = format!(
            "(def bar (move (box 20 4 4) 10 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (max (- x {x}) (- y 0.45) (- z 0.45)))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))"
        );
        solve(&src, "s", 0.1).max_nodal().unwrap().value
    };
    let (below, above) = (tip(1.2 - 1e-6), tip(1.2 + 1e-6));
    assert!((below / above - 1.0).abs() < 1e-5, "{below} below, {above} above");
}

/// A file loaded with import directories of its own reads this
/// distribution's forms (a study) and its imports from them, with no
/// `INK_PATH`; the kernel's own `load_file_importing` knows no study.
#[test]
fn a_file_given_import_dirs_loads_with_the_distributions_forms() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("sim-import-dirs");
    let (lib, away) = (dir.join("lib"), dir.join("away"));
    for d in [lib.join("parts"), away.clone()] {
        std::fs::create_dir_all(d).unwrap();
    }
    std::fs::write(lib.join("parts/block.ink"), "(import \"std\")\n(def block (move (box 4 4 4) 2 0 0))\n").unwrap();
    let model = away.join("model.ink");
    std::fs::write(
        &model,
        "(import \"std\")\n(import \"parts/block.ink\" as p)
(study s (shape p.block) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 0.5))
  (load tip (pressure -10MPa) (region (face p.block (normal 1 0 0)))))\n",
    )
    .unwrap();
    let mut d = ink::load_file_importing(&model, &[], std::slice::from_ref(&lib)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let (sol, _, _) = sim::run_study(&mut d, &ctx, "s", Some(0.5), 1e-9).unwrap();
    assert!(sol.report.converged && sol.max_nodal().unwrap().value > 0.0);
    assert!(ink::load_file(&model).is_err(), "the directory is given to no plain load");
}

#[test]
fn a_holds_curved_piece_reads_alike_whatever_scale_it_is_written_at() {
    // A 4 x 4 x 4 block (nu 0, pulled 10 MPa) held over x < 1.23 within 1.5
    // of its axis moves 0.031148 on cells of 0.2 mm. The same region with
    // the cylinder's field a thousandth of itself read 0.041874 while each
    // sample took the plane of the larger piece: where the plane's field
    // was the larger, the cylinder's face went untraced. Held instead over
    // x < 0.25 or within 1.3 of a point (a min: where the ball's field was
    // the smaller, the slab's face went untraced), likewise.
    let tip = |hold: &str| {
        let src = format!(
            "(def bar (move (box 4 4 4) 2 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall {hold})
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))"
        );
        solve(&src, "s", 0.2).max_nodal().unwrap().value
    };
    let bore = |k: &str| format!("(max (* {k} (- (len2 y z) 1.5)) (- x 1.23))");
    let ball = |k: &str| format!("(min (* {k} (- (len2 (len2 (- x 0.4) y) z) 1.3)) (- x 0.25))");
    for (region, plain) in [(&bore as &dyn Fn(&str) -> String, 0.031148), (&ball, 0.031770)] {
        let at_one = tip(&region("1"));
        assert!((at_one - plain).abs() < 1e-5, "{at_one}");
        for k in ["0.001", "1000"] {
            let scaled = tip(&region(k));
            assert!((scaled / at_one - 1.0).abs() < 1e-6, "{}: {scaled} against {at_one}", region(k));
        }
    }
}

/// A hold's piece of a gradient far from one is divided by it, a number
/// past ink's range made in Rust, not written: the model's claims, the
/// study's and one checked after it (which a number noted out of range
/// refused), are checked as written plainly.
#[test]
fn a_holds_piece_far_from_a_distances_scale_leaves_the_models_numbers_in_range() {
    let claims = |k: &str| {
        checks(&format!(
            "(def bar (move (box 4 4 4) 2 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (max (* {k} (- (len2 y z) 1.5)) (- x 1.23)))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))) (res 0.5))
(require stiff (max-displacement s < 1 (excluding buckling)))
(require big (volume bar > 10))"
        ))
    };
    let plain = claims("1");
    assert!(plain.iter().all(|c| c.1.status == ink::claims::Status::Holds), "{plain:?}");
    for k in ["(* 1e100 2)", "(* 1e-100 0.2)"] {
        let c = claims(k);
        assert!(c.iter().all(|c| c.1.status == ink::claims::Status::Holds), "{k}: {c:?}");
    }
}

#[test]
fn a_hold_reads_alike_whatever_scale_its_faces_fields_are_written_at() {
    // A 4 x 4 x 4 block (nu 0, pulled 10 MPa) held over its corner x < 1.23,
    // y < 0.37, z < 0.33 on cells of 0.2 mm moves 0.10327 (with exact polyhedral
    // face pressure integration; 0.10324 before, 0.10301 while each facet of
    // the pulled face went whole to its centre's cell, the corner node's share
    // short: on cells of 0.05 mm that read 0.10371, and this 0.10379). Each face is
    // traced where it is nearest, so the same region with the fields of y
    // and z written a thousand times themselves, or that of x, reads the
    // same: once 0.14004 and 0.07362, a face untraced in the samples where
    // the larger field was another's.
    let tip = |hold: &str| {
        let src = format!(
            "(def bar (move (box 4 4 4) 2 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall {hold})
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))"
        );
        solve(&src, "s", 0.2).max_nodal().unwrap().value
    };
    let plain = tip("(max (- x 1.23) (- y 0.37) (- z 0.33))");
    assert!((plain - 0.10327).abs() < 1e-5, "{plain}");
    for hold in ["(max (- x 1.23) (* 1000 (- y 0.37)) (* 1000 (- z 0.33)))", "(max (* 1000 (- x 1.23)) (- y 0.37) (- z 0.33))"] {
        let scaled = tip(hold);
        assert!((scaled / plain - 1.0).abs() < 1e-9, "{hold}: {scaled} against {plain}");
    }
}

#[test]
fn a_thin_lip_beyond_a_holds_face_reads_its_stretch() {
    // A bar 4 x 4 across (nu 0) held over x < c and pulled 10 MPa on its
    // end, a lip beyond the hold's face: the load's work is 160 N x 0.01 x
    // the lip, exactly. A planar lip 5% of a cell thick, with both faces
    // in one cell, reads within 1%. Lips under the 3% cell floor
    // (including when the lip's cell is dropped whole, as at c = 1.2
    // on 0.1 / 0.2 mm cells or c = 20.000001 on 0.5 mm cells), the study
    // marks the hold incomplete and no claim settles (neither Holds nor
    // Fails).
    let bar = |c: f64, lip: f64, res: f64| {
        let l = c + lip;
        format!(
            "(def bar (move (box {l} 4 4) {} 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x {c}))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))) (res {res}))",
            l / 2.0
        )
    };
    let (c, lip, res) = (1.234, 0.01, 0.2);
    let sol = solve(&bar(c, lip, res), "s", res);
    assert!(!sol.report.incomplete.iter().any(|e| matches!(e, ink_sim::Incomplete::Held(_))), "{:?}", sol.report.incomplete);
    let work = sol.report.energy;
    assert!((work / (1.6 * lip) - 1.0).abs() < 0.01, "{work} against {}", 1.6 * lip);
    for (c, lip, res, claim) in [
        (1.2, 0.004, 0.2, "(require work (compliance s < 0.004N*mm))"),
        (1.2, 0.004, 0.1, "(require work (compliance s < 0.004N*mm))"),
        (1.234, 0.004, 0.2, "(require work (compliance s < 0.004N*mm))"),
        (1.194, 0.01, 0.2, "(require work (compliance s > 0.02N*mm))"),
        (20.000001, 0.003999, 0.5, "(require disp (max-displacement s < 0.00003mm))"),
        (20.000001, 0.003999, 0.5, "(require work (compliance s < 0.004N*mm))"),
    ] {
        let src = bar(c, lip, res);
        let sol = solve(&src, "s", res);
        assert!(
            sol.report.incomplete.iter().any(|e| matches!(e, ink_sim::Incomplete::Held(_))),
            "c {c}, lip {lip}, res {res}: expected Incomplete::Held, got {:?}",
            sol.report.incomplete
        );
        let ck = checks(&format!("{src}\n{claim}")).remove(0).1;
        assert_eq!(ck.status, ink::claims::Status::Unknown, "c {c}, lip {lip}, res {res}: {:?} at {:?}", ck.status, measured(&ck));
    }
    // A convection film or volume force on a dropped-whole lip cell also
    // leaves its study incomplete and its claims unsettled.
    let film = "(def bar (move (box 20.004 4 4) 10.002 0 0))
(study s (shape bar) (heat (k 1W/m/K))
  (temperature cold (- x 20.000001) 20degC)
  (convection tip (face bar (normal 1 0 0)) (coefficient 100W/m2/K) (ambient 100degC)) (res 0.5))
(require hot (max-temperature s < 25degC))";
    assert!(solve(film, "s", 0.5).report.incomplete.iter().any(|e| matches!(e, ink_sim::Incomplete::Held(_))));
    assert_eq!(checks(film).remove(0).1.status, ink::claims::Status::Unknown);

    let force = "(def bar (move (box 20.004 4 4) 10.002 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 20.000001))
  (load tip (force 160N 0 0) (region (- 19.9 x))) (res 0.5))
(require disp (max-displacement s < 0.00003mm))";
    assert!(solve(force, "s", 0.5).report.incomplete.iter().any(|e| matches!(e, ink_sim::Incomplete::Held(_))));
    assert_eq!(checks(force).remove(0).1.status, ink::claims::Status::Unknown);
}

#[test]
fn a_thin_lips_face_is_read_through_a_cell_with_free_material() {
    // A 20.2 x 4 x 4 bar (nu 0) held over x < 20 and pulled 10 MPa on its
    // end, on cells of 0.5 mm: the lip beyond the hold is 0.2 mm of a cell
    // two fifths full, and its face moves 0.002 mm. A face in a cell under
    // half full is read from the cell under it, which here is all held, and
    // read the held value: 0 on the claim's two coarser grids, so 0.00205
    // +- 0.0061 (and 0 +- 0, "holds", on a lip of 0.09 mm, while the
    // compliance of the same solve said it moved). Read from its own cell:
    // +- 0.00014. (A lip thinner than a sample, 0.125 mm here, is not
    // settled at all: this one the cells carry.)
    let src = "(def bar (move (box 20.2 4 4) 10.1 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 20))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))) (res 0.5))
(require d (max-displacement s < 1000mm))";
    let (v, lo, hi) = claim_band(&checks(src).remove(0).1).unwrap();
    assert!(lo <= 0.002 && 0.002 <= hi && hi - lo < 0.001, "{v} in [{lo}, {hi}] against 0.002");
}

#[test]
fn a_volume_load_on_a_dropped_lip_is_caught_beside_a_face_fix_or_across_refinement() {
    // The end past x = 4 is dropped. Load a corner next to a cell with
    // both a face and volume fix, or put a material boundary in the lip so
    // its dropped leaves are finer than the wall cells that remain.
    for (end, center, e, fixes, region) in [
        (4.004, 2.002, "1000MPa", "(fix side (face bar (normal 0 1 0)) (dofs x y) (support)) (fix wall (- x 4.000001) (dofs x z) (support))", "(max (- 3.9 x) (- 1.9 y))"),
        (4.001, 2.0005, "(select (- x 4.0005) 1000MPa 1001MPa)", "(fix wall (- x 4.000001))", "(- 3.9 x)"),
    ] {
        let src = format!(
            "(def bar (move (box {end} 4 4) {center} 0 0))
(study s (shape bar) (elastic (E {e}) (nu 0))
  {fixes} (load tip (force 160N 0 0) (region {region})) (res 0.5))
(require disp (max-displacement s < 0.00003mm))"
        );
        let sol = solve(&src, "s", 0.5);
        assert!(sol.report.sliver_volume > 0.0 && sol.report.incomplete.iter().any(|e| matches!(e, ink_sim::Incomplete::Held(_))), "{:?}", sol.report.incomplete);
        let ck = checks(&src).remove(0).1;
        assert_eq!(ck.status, ink::claims::Status::Unknown, "{ck:?}");
    }
}

#[test]
fn a_thin_one_axis_hold_holds_under_a_pressure_across_it() {
    // A 6 mm cube (nu 0) on thin one-axis roller layers (1% to 3.1% of a 1 mm
    // cell), pressed 10 MPa on x = 6 and 7 MPa on y = 6: held pieces under
    // MIN_FRACTION keep their free axes' stiffness (without incompatible
    // modes), so the layer no longer slips under the transverse pressure.
    for t in [0.01, 0.02, 0.029, 0.031] {
        let two = format!(
            "(def blk (move (box 6 6 6) 3 3 3))
(study s (shape blk) (elastic (E 1000MPa) (nu 0))
  (mirror mz z 0) (fix fx (- x {t}) (dofs x) (support)) (fix fy (- y {t}) (dofs y) (support))
  (load px (pressure -10MPa) (region (face blk (normal 1 0 0)))) (load py (pressure -7MPa) (region (face blk (normal 0 1 0)))) (res 1))"
        );
        let truth = ((10.0f64 * (6.0 - t)).powi(2) + (7.0f64 * (6.0 - t)).powi(2)).sqrt() / 1000.0;
        let v = solve(&two, "s", 1.0).max_nodal().unwrap().value;
        assert!((v - truth).abs() < 1e-5 * truth, "two layers t={t}: {v} against {truth}");
    }
    // Single 0.02 mm layer with a y-mirror: no slip along x at x = 0 or x = 6
    // (previously u_x slipped by 7.1e-5 at x = 0 and read 0.0731312), and the
    // claim on 0.5 mm cells settles with truth inside its band.
    let src = "(def blk (move (box 6 6 6) 3 3 3))
(study s (shape blk) (elastic (E 1000MPa) (nu 0))
  (mirror mz z 0) (mirror my y 0) (fix fx (- x 0.02) (dofs x) (support))
  (load px (pressure -10MPa) (region (face blk (normal 1 0 0)))) (load py (pressure -7MPa) (region (face blk (normal 0 1 0)))) (res 0.5))
(require d (max-displacement s < 1000mm))";
    let sol = solve(src, "s", 1.0);
    let ux_far = sol.value_at([5.999, 5.999, 3.0]).unwrap()[0];
    let ux_layer = sol.value_at([0.001, 5.999, 3.0]).unwrap()[0];
    assert!((ux_far - 0.05979).abs() < 1e-5 * 0.05979, "ux_far {ux_far} against 0.05979");
    assert!((ux_layer - -0.00019).abs() < 1e-6, "ux_layer {ux_layer} against -0.00019");
    let truth = ((10.0f64 * 5.98).powi(2) + (7.0f64 * 6.0).powi(2)).sqrt() / 1000.0;
    let (v, lo, hi) = claim_band(&checks(src).remove(0).1).unwrap();
    assert!(lo <= truth && truth <= hi && (v - truth).abs() < 1e-4 * truth, "{v} in [{lo}, {hi}] against {truth}");
}

#[test]
fn a_cubes_far_corner_reads_its_peak_however_its_faces_are_held() {
    // A 6 mm cube (nu 0) pressed 10, 7 and 4 MPa on three faces meeting at a
    // corner, held on the three opposite faces by roller layers 0.02 mm thick
    // (or by mirrors): three uniaxial stretches, largest at the far corner,
    // 0.0768145 mm (0.0770714 mirrored). The field was right there; its peak
    // was read a cell short along an edge, at (6, 6, 5.5): the corner node,
    // one cell's alone, was not counted as reliable. On cells of 0.5 mm it
    // read 0.0762156 +- 0.00037 and the claim held; 0.25 mm cells halved the
    // miss, one layer or two read exact (no load along the edge read).
    let cube = |holds: &str| {
        format!(
            "(def blk (move (box 6 6 6) 3 3 3))
(study s (shape blk) (elastic (E 1000MPa) (nu 0))
  {holds}
  (load px (pressure -10MPa) (region (face blk (normal 1 0 0))))
  (load py (pressure -7MPa) (region (face blk (normal 0 1 0))))
  (load pz (pressure -4MPa) (region (face blk (normal 0 0 1)))) (res 0.5))
(require d (max-displacement s < 1000mm))"
        )
    };
    for (holds, truth) in [
        ("(fix fx (- x 0.02) (dofs x) (support)) (fix fy (- y 0.02) (dofs y) (support)) (fix fz (- z 0.02) (dofs z) (support))", 0.0768145),
        ("(mirror mx x 0) (mirror my y 0) (mirror mz z 0)", 0.0770714),
    ] {
        let (v, lo, hi) = claim_band(&checks(&cube(holds)).remove(0).1).unwrap();
        assert!(lo <= truth && truth <= hi && (v - truth).abs() < 1e-4 * truth, "{holds}: {v} in [{lo}, {hi}] against {truth}");
    }
}

#[test]
fn a_pulled_bars_compliance_holds_wherever_its_sides_fall_in_the_cells() {
    // A 20 x 2 x 2 bar (nu 0) held over x < 1.4 and pulled 10 MPa, its sides
    // at y and z = -1.95 and 0.05, a tenth of a 0.5 mm cell past the planes:
    // the load's work is 40 N x 0.186 mm = 7.44 N mm. It read 7.432144 +-
    // 0.0063 on cells of 0.25 mm (7.426 and 7.405 on 0.5 and 1 mm,
    // converging to the wrong place): each sample an edge of the bar
    // crosses was counted by its nearer side's plane alone, material past
    // the other side with it, by where the edge fell in the sample
    // (`occupied_parts`). Counted within both it reads 7.435014 +- 0.0329
    // (7.441 and 7.430 on 0.5 and 1 mm); what is left is its pressure's
    // force, 0.03% short (the test below).
    let src = "(def bar (move (box 20 2 2) 10 -0.95 -0.95))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 1.4))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))) (res 0.5))
(require c (compliance s < 1e9N*mm))";
    let (v, lo, hi) = claim_band(&checks(src).remove(0).1).unwrap();
    assert!(lo <= 7.44 && 7.44 <= hi, "{v} in [{lo}, {hi}] against 7.44");
}

#[test]
fn a_bars_edges_hold_their_material_wherever_they_fall_in_the_cells() {
    // The 20 x 2 x 2 bar on cells of 0.5 mm, 80 mm^3, its sides in y and z
    // 0.1, 0.2 and 0.25 mm past the cells' planes, and one side alone 0.05
    // past: each sample an edge crosses holds the bar's part within both
    // sides. (Counted by the nearer side's plane alone it read 80.125 at
    // 0.1, 0.1, 0.16% over, and 80.0625 on cells of 0.25 mm.)
    for (y, z) in [(0.1, 0.1), (0.2, 0.2), (0.25, 0.25), (0.05, 0.0)] {
        let src = format!(
            "(def bar (move (box 20 2 2) 10 {} {}))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 1.4))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))",
            y - 1.0,
            z - 1.0
        );
        let v = solve(&src, "s", 0.5).report.volume;
        assert!((v / 80.0 - 1.0).abs() < 1e-9, "sides {y}, {z} past: {v} mm^3 for 80");
    }
}

#[test]
fn a_pressed_faces_edge_moves_with_its_middle() {
    // A 20 x 4 x 4 bar (nu 0) held over x < 1.2 and pulled 10 MPa on its
    // end, on cells of 0.3 mm, its sides 0.2 mm past their planes: every
    // point of the end moves alike, 0.188 mm. A facet of the pressed face
    // put whole on the nodes of the cell its centre lay in, though it
    // crossed a cell's plane, gave the nodes at the face's edges the wrong
    // shares: the end's corner moved 0.030% further than its middle
    // (0.188025 for 0.187968), and moves 0.0003% further. (With the bar's
    // edges holding their own material and the facets whole, the corner
    // read 0.188126, 0.067% over, 0.188178 on cells of 0.6 mm, and a claim
    // on it held with the truth outside its band.)
    let src = "(def bar (move (box 20 4 4) 10 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 1.2))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))";
    for res in [0.3, 0.6] {
        let s = solve(src, "s", res);
        // (the middle read 0.01 mm inside the end, and stretched to it)
        let (edge, middle) = (s.max_nodal().unwrap().value, s.value_at([19.99, 0.0, 0.0]).unwrap()[0] * 18.8 / 18.79);
        assert!((edge / middle - 1.0).abs() < 1e-4 && (edge / 0.188 - 1.0).abs() < 5e-5, "cells of {res}: {edge} at the edge, {middle} in the middle, for 0.188");
    }
}

#[test]
fn the_smoothers_estimate_reaches_what_a_few_cells_carry() {
    // A bar held over a cylinder of radius 1.37 on cells of 0.5 mm: the
    // largest eigenvalue the smoother is tuned to is carried by the few
    // cells of the hold's edge, 3.74 where 15 steps of a power iteration
    // gave 3.32, and the solve took 32 iterations for 12 (`max_eigen`).
    let src = "(body b (shape (move (box 20 2 2) 10 0 0)))
(study s (shape (shape-of b)) (elastic (E 68900MPa) (nu 0.33))
  (fix wall (- (len2 x z) 1.37)) (load tip (force 0 0 -10N) (region (- 19.5 x))))";
    let n = solve(src, "s", 0.5).report.iterations;
    assert!(n <= 16, "{n} iterations");
    // A bar held in x a millionth short of a cell plane where a guide of y
    // begins, a third fix of z beside (0.3 mm cells, where the hold keeps
    // its own sliver's large constant): 15 power steps took 631 iterations;
    // Lanczos `max_eigen` solves it in <= 25.
    let sliver = "(def bar (move (box 20 2 2) 10 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0.3))
  (fix fx (- x 1.199999) (dofs x) (support))
  (fix fy (- 1.2 x) (dofs y) (support))
  (fix fz (- x 1.2) (dofs z) (support))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))";
    let sol = solve(sliver, "s", 0.3);
    assert!(sol.report.converged && sol.report.iterations <= 25, "{} iterations", sol.report.iterations);
}

#[test]
fn a_pressure_lands_whole_on_a_face_whose_edge_cells_are_thin() {
    // The end's pressure is 40 N. Its corner cell is 4% full on 0.25 mm
    // cells, where the contour's facets gave only 39.987 N.
    let src = "(def bar (move (box 20 2 2) 10 -0.95 -0.95))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 1.4))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))";
    let f = solve(src, "s", 0.25).report.applied[0].1[0].abs();
    assert!((f / 40.0 - 1.0).abs() < 5e-5, "{f} N for 40");
}

#[test]
fn a_face_crossing_a_cells_plane_reads_alike_either_side() {
    // A 20 x 2 x (2 + d) bar (nu 0) held over x < 1.5 and pulled 10 MPa on
    // its end, on cells of 0.5 mm, its top face d past the plane z = 1: the
    // end moves 0.185 mm whatever d is. With d under 0.015 the cells over
    // the plane are under 3% full, slivers, and were dropped, their
    // material and what pressed on their face with them: the end's edge
    // read 0.26%, 1.1%, 2.6% and 4.2% too far at d = 0.001, 0.005, 0.01 and
    // 0.014, and 0.05% at 0.016, where the cells are kept. The cell under
    // each now takes it (`taken`).
    for d in [-0.005, 0.0, 0.001, 0.005, 0.01, 0.014, 0.016, 0.02] {
        let src = format!(
            "(def bar (max (- (abs (- x 10)) 10) (- (abs y) 1) (- z {}) (- -1 z)))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 1.5))
  (load tip (pressure -10MPa) (region (- 19.999 x))))",
            1.0 + d
        );
        let v = solve(&src, "s", 0.5).max_nodal().unwrap().value;
        assert!((v / 0.185 - 1.0).abs() < 1e-3, "top face {d} past the plane: {v} for 0.185");
    }
}

#[test]
fn a_column_buckles_alike_as_its_side_crosses_a_cells_plane() {
    // A 20 x (2 + d) x 2 column (nu 0.3) held over x < 1.4 and pushed 1 N
    // along itself, on cells of 0.5 mm, its side d past the plane y = 1:
    // its factor grows with its width, 9.7275 at d = 0 and 9.8072 at 0.016.
    // With d under 0.015 the cells over the plane are slivers: dropped,
    // the factor stood at 9.7275 until 0.016; taken into the stiffness and
    // not into the geometric stiffness it read 9.7742 at 0.005, high, the
    // unsafe side. It reads 9.7508 (9.7524 by the width).
    let factor = |d: f64| {
        let src = format!(
            "(def bar (max (- 0 x) (- x 20) (- -1 y) (- y {}) (- -1 z) (- z 1)))
(study s (shape bar) (elastic (E 1000MPa) (nu 0.3))
  (fix wall (- x 1.4))
  (load tip (force -1N 0 0) (region (- 19.5 x))))",
            1.0 + d
        );
        solve(&src, "s", 0.5).report.buckle.expect("a buckling analysis").factor
    };
    let (on, sliver, kept) = (factor(0.0), factor(0.005), factor(0.016));
    let by_width = on + (kept - on) * 0.005 / 0.016;
    assert!((sliver / by_width - 1.0).abs() < 5e-4, "{sliver} with the side 0.005 past, {by_width} by its width ({on} on the plane, {kept} at 0.016)");
}

#[test]
fn a_pressed_end_takes_no_force_sideways() {
    // A 20 x 2.012 x 2 bar (nu 0) held over x < 1.4 and pulled 10 MPa on
    // its end, on cells of 0.5 mm: the end is flat and moves 0.186 mm. A
    // facet beside an edge lies a little inside the part, and a point of it
    // there read the other face's normal, turned too far, and took the
    // facet's own, tilted a hundredth of a radian: 7.2e-4 N of the 40.24
    // went sideways, enough to bend the bar, and the end's edge read
    // 0.035% too far.
    let src = "(def bar (max (- 0 x) (- x 20) (- -1 y) (- y 1.012) (- -1 z) (- z 1)))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 1.4))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))";
    let s = solve(src, "s", 0.5);
    let (f, v) = (s.report.applied[0].1, s.max_nodal().unwrap().value);
    assert!(f[1].abs() < 1e-9 * f[0].abs() && f[2].abs() < 1e-9 * f[0].abs(), "{f:?}");
    assert!((v / 0.186 - 1.0).abs() < 2e-4, "{v} for 0.186");
}

#[test]
fn a_bars_weight_is_its_fixs_to_carry_slivers_and_all() {
    // A 20 x 2.02 x 2.01 bar held over x < 1.5 under its own weight, on
    // cells of 0.5 mm: three of its sides 0.01 mm past the cells' planes,
    // their cells 2% full. The fix carries the bar's weight, rho g V.
    // Dropped, those cells' weight was missing: 1.48% of it; taken by the
    // cells beside them where the bar is free, 0.11% (the slivers under
    // the fix, which no cell takes).
    let src = "(def bar (max (- (abs (- x 10)) 10) (- (abs y) 1.01) (- z 1.01) (- -1 z)))
(study s (shape bar) (elastic (E 1000MPa) (nu 0) (density 2000kg/m3))
  (fix wall (- x 1.5))
  (gravity 0 0 -9.81m/s2))";
    let s = solve(src, "s", 0.5);
    let (weight, carried) = (2000e-9 * 9.81 * 20.0 * 2.02 * 2.01, s.report.reactions.iter().find(|r| r.0 == "wall").unwrap().1[2]);
    assert!((carried / weight - 1.0).abs() < 1e-8, "{carried} N carried of {weight}");
}

#[test]
fn a_pressure_on_a_wedges_faces_lands_whole() {
    // A 10 degree wedge 8 mm long and 3 wide, held over x > 7 and pressed 1
    // MPa on its slope and base: the base takes 24 N and the slope takes
    // as much across it. The contour rounded off the tip of the wedge.
    let src = "(def part (max (- 0 y) (- (* 0.984807753012208 y) (* 0.17364817766693033 x)) (- x 8) (- (abs z) 1.5)))
(study s (shape part) (elastic (E 1000MPa) (nu 0.3))
  (fix wall (- 7 x))
  (load slope (pressure 1MPa) (region (face part (normal -0.17364817766693033 0.984807753012208 0))))
  (load base (pressure 1MPa) (region (face part (normal 0 -1 0)))))";
    let sol = solve(src, "s", 0.5);
    let f = sol.report.applied.iter().find(|a| a.0 == "base").unwrap().1[1];
    assert!((f / 24.0 - 1.0).abs() < 1e-3, "{f} N for 24");
    let check = |sol: &sim::Solution, load: &str, want: [f64; 3]| {
        let f = sol.report.applied.iter().find(|a| a.0 == load).unwrap().1;
        for d in 0..3 {
            assert!((f[d] - want[d]).abs() < 1e-3, "{load}: {f:?} N for {want:?}");
        }
    };
    check(&sol, "slope", [24.0 * 10.0f64.to_radians().tan(), -24.0, 0.0]);

    // Translate a steeper wedge so the base and both sides fall off the grid.
    let src = "(def part (max (- 0.03 y) (- (* 0.8660254037844386 (- y 0.03)) (* 0.5 (- x 0.07))) (- x 8.07) (- (abs (- z 1.61)) 1.5)))
(study s (shape part) (elastic (E 1000MPa) (nu 0.3))
  (fix wall (- 7 x))
  (load slope (pressure 1MPa) (region (face part (normal -0.5 0.8660254037844386 0))))
  (load base (pressure 1MPa) (region (face part (normal 0 -1 0)))))";
    let sol = solve(src, "s", 0.5);
    check(&sol, "base", [0.0, 24.0, 0.0]);
    check(&sol, "slope", [24.0 * 30.0f64.to_radians().tan(), -24.0, 0.0]);

    // The tetrahedron's pressed face projects to a right triangle of
    // side 6 mm on each coordinate plane.
    let src = "(def part (max (- 0.03 x) (- 0.07 y) (- 0.11 z) (- (* 0.5773502691896258 (+ x y z)) 3.5853451716676)))
(study s (shape part) (elastic (E 1000MPa) (nu 0.3))
  (fix wall (- (len3 x y z) 1.5))
  (load slope (pressure 1MPa) (region (face part (normal 0.5773502691896258 0.5773502691896258 0.5773502691896258)))))";
    check(&solve(src, "s", 0.5), "slope", [-18.0; 3]);
}

#[test]
#[ignore = "a rod of radius 1.03 on 0.5 mm cells holds 0.2% too much: each sample's plane lies outside the curve"]
fn a_round_bars_material_is_its_own() {
    // A rod of radius 1.03 mm along x, 10 mm long, on cells of 0.5 mm:
    // 33.3292 mm^3. It reads 33.3955, 0.20% over (0.10% under while the
    // cells under 3% full were dropped).
    let src = "(def part (max (- (len2 (- y 0.21) (- z 0.07)) 1.03) (- (abs (- x 5.03)) 5)))
(study s (shape part) (heat (k 0.2W/m/K)) (temperature cold (- x 0.4) 20degC)
  (source q (power 0.01W) (region (- x 1000))))";
    let (v, truth) = (solve(src, "s", 0.5).report.volume, std::f64::consts::PI * 1.03 * 1.03 * 10.0);
    assert!((v / truth - 1.0).abs() < 5e-4, "{v} mm^3 for {truth}");
}

#[test]
#[ignore = "a bar's sliver cells under its hold are dropped still: 79.9925 mm^3 of its 80"]
fn a_bars_material_is_all_counted_where_its_edge_cell_is_a_sliver() {
    // A 20 x 2 x 2 bar, its sides at y = 0.1 and z = 0.05 past the planes
    // of 0.5 mm cells: the cell along that edge holds 0.1 x 0.05 mm of the
    // bar's section, 2% of itself. Where the bar is free the cell beside
    // it takes that material; where the fix reaches (x < 1.5, three cells)
    // it is dropped as it was everywhere, 79.9 mm^3 for 80 then and the
    // load's work 0.1% over.
    let src = "(def bar (move (box 20 2 2) 10 -0.9 -0.95))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 1.4))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))";
    let v = solve(src, "s", 0.5).report.volume;
    assert!((v / 80.0 - 1.0).abs() < 1e-6, "{v} mm^3 for 80");
}

#[test]
fn a_guide_over_a_walls_edge_cell_leaves_the_walls_edge_where_it_is() {
    // A 20 x 2 x 2 bar (nu 0) held over x < 1.3 and guided in y over x 1.0 to
    // 1.6, pulled 10 MPa: it stretches 0.01 x (20 - 1.3) = 0.187 mm, the
    // guide holding nothing along x. On cells of 0.5 mm the wall's edge cell
    // holds no material free of every fix, and was taken as held along x
    // too: the nodes of its far face only it reaches were held hard there,
    // and it read 0.185, as held to x = 1.5.
    let src = "(def bar (move (box 20 2 2) 10.0 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x 1.3))
  (fix g (- (abs (- x 1.3)) 0.3) (dofs y) (support))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))";
    let v = solve(src, "s", 0.5).max_nodal().unwrap().value;
    assert!((v / 0.187 - 1.0).abs() < 2e-3, "{v} against 0.187");
}

#[test]
fn a_wall_a_millionth_short_of_a_cell_plane_beside_a_guide_holds_its_bar() {
    // A 20 x 2 x 2 bar (nu 0) held over x < 0.999999 and guided in y over
    // x > 1, pulled 10 MPa: it stretches 0.01 x 19 = 0.19 mm. Previously
    // the wall's edge was dropped as untraced because the guide held y in
    // the cell beyond, reading 3765 mm on 0.5 mm cells (7530 on 1, 1883 on
    // 0.25). Also check a wider 0.004 mm sliver under MIN_FRACTION.
    for (cut, truth) in [(0.999999, 0.19), (0.996, 0.19004)] {
        let src = format!(
            "(def bar (move (box 20 2 2) 10.0 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- x {cut}))
  (fix g (- 1 x) (dofs y) (support))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))"
        );
        for res in [1.0, 0.5, 0.25] {
            let v = solve(&src, "s", res).max_nodal().unwrap().value;
            assert!((v / truth - 1.0).abs() < 2e-3, "cut {cut}, res {res}: {v} against {truth}");
        }
    }
}

#[test]
fn a_column_on_a_held_floor_band_sinks_as_its_free_height_says() {
    // A 3 x 3 x 3 block (a quarter, by mirrors) on a roller floor held over z
    // < t, E 1 MPa, nu 0, under its own weight: its top sinks rho g (3 -
    // t)^2 / 2E exactly (linear cells are exact at their nodes), and the
    // floor carries all its weight, the band's too. With the band's weight
    // on its cell's nodes it read (t / (3 - t))^2 less, alike on cells of
    // 0.3 and 0.15 mm (t 0.05: 3.3499e-4 against 3.3508e-4), its three
    // grids alike to 1e-10, so a claim's band missed the truth.
    for t in [0.01, 0.05] {
        let study = format!(
            "(body bar (shape (move (box 3.0 3.0 3.0) 1.5 1.5 1.5)))
(study s (shape (shape-of bar)) (elastic (E 1MPa) (nu 0) (density 7850kg/m3))
  (mirror mx x 0) (mirror my y 0) (fix floor (- z {t}) (dofs z) (support)) (gravity 0 0 -9.81m/s2) (res 0.3))"
        );
        let (v, lo, hi) = claim_band(&checks(&format!("{study}\n(require sinks (max-displacement s < 1mm))")).remove(0).1).unwrap();
        let truth = 7850.0 * 9.81e-9 * (3.0f64 - t).powi(2) / 2.0;
        assert!(lo <= truth && truth <= hi && (v / truth - 1.0).abs() < 1e-6, "t {t}: {v} in [{lo}, {hi}] against {truth}");
        let s = solve(&study, "s", 0.3);
        // (No stress lies across its weight, so a turn about it is a free
        // coordinate of its buckling, not a buckle: solved closer, here, the
        // turn's force came out 1e-15 of the prestress's and not zero, the
        // buckling "not resolved", and the claim above, asked again and
        // given this solve from the cache, not settled.)
        assert!(s.report.buckle.as_ref().is_some_and(|b| b.failure.is_none()), "t {t}: {:?}", s.report.buckle);
        let (weight, floor) = (7850.0 * 9.81e-9 * 27.0, s.report.reactions.iter().find(|r| r.0 == "floor").unwrap().1[2]);
        assert!((floor / weight - 1.0).abs() < 1e-6 && imbalance(&s) < 1e-6, "t {t}: the floor carries {floor} of {weight}; imbalance {}", imbalance(&s));
    }
}

#[test]
fn a_source_over_a_band_held_at_a_temperature_heats_nothing_beyond_it() {
    // A 3 x 1 x 1 bar held at 20 C over x < t and heated 0.01 W all through:
    // the band's share of the heat is the fix's, and the bar beyond rises as
    // it does with the same heat kept off the band (and the fix takes all
    // 0.01 W either way it is written). With the band's share on its cell's
    // nodes the bar rose (t / (3 - t))^2 less: 2.8e-4 of its rise at t 0.05,
    // 1.2e-3 at 0.1, alike on every grid.
    for t in [0.05, 0.1] {
        let bar = |power: f64, region: &str| {
            let src = format!(
                "(body bar (shape (move (box 3 1 1) 1.5 0 0)))
(study s (shape (shape-of bar)) (heat (k 0.5W/m/K))
  (temperature cold (- x {t}) 20degC) (source q (power {power}W) (region {region})))"
            );
            let s = solve(&src, "s", 0.3);
            (s.max_nodal().unwrap().value, s.report.reactions.iter().find(|r| r.0 == "cold").unwrap().1[0])
        };
        let (all, off) = (bar(0.01, "-1"), bar(0.01 * (3.0 - t) / 3.0, &format!("(- {t} x)")));
        let rise = 0.01 / 3.0 * (3.0 - t).powi(2) / (2.0 * 0.0005);
        assert!(((all.0 - off.0) / rise).abs() < 1e-6, "t {t}: {} heated all through, {} with the band left out", all.0, off.0);
        assert!((all.1 / -0.01 - 1.0).abs() < 1e-9, "t {t}: the fix takes {} W of 0.01", all.1);
    }
}

/// A 3 x 3 x 3 block clamped over z < t and pressed 0.001 MPa on its side x =
/// 3: how far it moves pressed over the whole side, and over z > t alone.
fn pressed_over_a_clamped_band(t: f64) -> (f64, f64) {
    let moves = |region: &str| {
        let src = format!(
            "(body bar (shape (move (box 3 3 3) 1.5 0 1.5)))
(study s (shape (shape-of bar)) (elastic (E 1MPa) (nu 0))
  (fix floor (- z {t})) (load side (pressure 0.001MPa) (region {region})))"
        );
        solve(&src, "s", 0.3).max_nodal().unwrap().value
    };
    (moves("(- 2.999 x)"), moves(&format!("(max (- 2.999 x) (- {t} z))")))
}

#[test]
fn a_pressure_on_a_clamped_bands_face_is_the_clamps_to_carry() {
    // The band's part of the side is the clamp's to carry, and the block
    // moves as it does pressed over z > t alone. On its cells' nodes with
    // the rest it moved 8.0e-5 less at t 0.05 and 3.2e-4 less at 0.1
    // (0.00880232 against 0.00880515); each facet cut by the clamp's
    // region, 1.0e-5 and 1.4e-5.
    for t in [0.05, 0.1] {
        let (all, off) = pressed_over_a_clamped_band(t);
        assert!((all / off - 1.0).abs() < 2e-5, "t {t}: {all} pressed over the whole side, {off} off the band");
    }
}

#[test]
#[ignore = "a pressure on a clamped band's face still moves the free material beside it 1e-5: where the load's region ends in the facet the clamp's region ends in, and the strips at the face's edge"]
fn a_pressure_on_a_clamped_bands_face_moves_nothing() {
    for t in [0.05, 0.1] {
        let (all, off) = pressed_over_a_clamped_band(t);
        assert!((all / off - 1.0).abs() < 1e-6, "t {t}: {all} pressed over the whole side, {off} off the band");
    }
}

#[test]
fn a_hold_slab_within_one_sample_beside_a_sliver_is_split_or_refused_not_slipped() {
    // A 20 x 4 x 4 bar (nu 0, pulled 10 MPa) held over x 0.03 to 0.11, a 0.03
    // mm sliver free beyond it at the end, on cells of 0.5 mm (samples of
    // 0.125): both the slab's faces lay in one sample, whose one tangent
    // plane is the face toward the sliver, so the face toward the bar was
    // never traced, and the hold slipped: 0.20775 +- 0.0013, past 0.2, the
    // whole bar's stretch unheld, against 0.1989. Now it is split as a thin
    // slab, and answers or is refused by name; never past the bar's stretch.
    let src = "(def bar (move (box 20 4 4) 10 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0))
  (fix wall (- (abs (- x 0.07)) 0.04))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))) (res 0.5))
(require d (max-displacement s < 1000mm))";
    let mut d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let v = match sim::run_study(&mut d, &ctx, "s", Some(0.5), 1e-9) {
        Err(e) => {
            assert!(format!("{e:?}").contains("holds a slab of the part thinner than its cells resolve"), "{e:?}");
            return;
        }
        Ok((sol, _, _)) => sol.max_nodal().unwrap().value,
    };
    // It reads 0.198848 (0.198469 before each piece of the edge read its
    // traction partly in the free cell beside it). On the claim's two
    // coarser grids the sliver has no sample's centre and the cells are not
    // split: both the slab's faces are traced there now, and the one field
    // of a cell serves the sliver and the bar 0.16% short (0.198586,
    // 0.198607), where one face traced held the sliver with the slab and
    // read the bar exactly. The three do not converge and the claim is not
    // settled; a band it states holds the truth.
    assert!((v / 0.1989 - 1.0).abs() < 1e-3, "{v} against 0.1989");
    // (Not settled, and it says so: were it settled, its band is what this
    // would hold to the truth, as it did when one face was traced.)
    let claim = checks(src).remove(0).1;
    assert!(claim.status == ink::claims::Status::Unknown && claim_band(&claim).is_err(), "{:?}: {:?}", claim.status, claim_band(&claim));
    // Thinner, x 0.055 to 0.07 beside a sliver 0.055 thick: unresolved at
    // the floor, refused by name, saying to extend the hold to the face.
    let thin = src.replace("(- (abs (- x 0.07)) 0.04)", "(- (abs (- x 0.0625)) 0.0075)");
    let mut d = ink::load_str(&format!("(import \"std\")\n{thin}"), "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let e = sim::run_study(&mut d, &ctx, "s", Some(0.5), 1e-9).err().map(|e| format!("{e:?}")).unwrap_or_default();
    assert!(e.contains("extend the hold to the face"), "{e}");
}

#[test]
fn a_hold_a_millionth_thick_at_an_end_face_holds_and_one_off_the_part_is_refused() {
    // A bar off the grid, held over its end's first 1e-6 mm (the dynamics
    // hunt's columns, 45 of them): its samples' points, half a sample about
    // the held piece, fell either side of the band, and it was refused as
    // touching no material. Its tip against the same bar held over 0.01 mm
    // (its cut sides stiffen both alike): as the free lengths, 19.999999 to
    // 19.99.
    let bar = |wall: &str| {
        format!(
            "(import \"std\")
(body b (shape (move (box 20 2 2) 23.37 0.21 -0.13)))
(study s (shape (shape-of b)) (elastic (E 1000MPa) (nu 0))
  (fix wall {wall}) (load tip (pressure -10MPa) (region (face b (normal 1 0 0)))))"
        )
    };
    let tip = |wall: &str| solve(&bar(wall), "s", 0.5).max_nodal().unwrap().value;
    let (thin, thick) = (tip("(- x 13.370001)"), tip("(- x 13.38)"));
    assert!((thin / thick / (19.999999 / 19.99) - 1.0).abs() < 1e-3, "{thin} {thick}");
    // A region a cell and more off the part holds nothing: refused.
    let mut d = ink::load_str(&bar("(- x 12)"), "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let e = sim::run_study(&mut d, &ctx, "s", Some(0.5), 1e-9).err().expect("refused");
    assert!(format!("{e:?}").contains("touches no material"), "{e:?}");
}

#[test]
fn a_lip_a_hair_thick_beyond_a_holds_face_is_a_lip() {
    // The plate held to 1e-10 mm short of its top face (and on cells a
    // thousand times as large, 1e-5 mm short), pressed on that face: the
    // lip is free material, however thin, and what presses on it is not
    // settled. (Taken for the face within a millionth of a sample, it went
    // untraced, the pressure went into the fix's reaction, and the study
    // settled: its claim held.)
    for (k, gap) in [(1.0, "1e-10"), (1000.0, "1e-5")] {
        let src = format!(
            "(def plate (move (box {} {} {}) {} 0 {}))
(study s (shape plate) (elastic (E 200000MPa) (nu 0.3))
  (fix base (max (- x {}) (- z (- {} {gap}))))
  (load top (pressure 0.00001MPa) (region (max (- x {}) (- {} z))))
  (load tip (force 0 0 -0.01N) (region (- {} x))))",
            6.0 * k,
            2.0 * k,
            0.2 * k,
            3.0 * k,
            0.1 * k,
            k,
            0.2 * k,
            k,
            0.15 * k,
            5.5 * k
        );
        let mut d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
        let ctx = d.ctx(&[], 0.0).unwrap();
        let (sol, _, _) = sim::run_study(&mut d, &ctx, "s", Some(0.5 * k), 1e-9).unwrap();
        assert!(!sol.report.incomplete.is_empty(), "gap {gap}: {:?}", sol.report.warnings);
    }
}

#[test]
fn a_holds_face_on_the_parts_own_is_no_edge() {
    // A plate 0.2 mm thick on 0.5 mm cells held over x < 1 and z under its
    // top face, the region's face on the plate's (level, off the cells'
    // planes, tilted): no free material beyond it, so no edge between held
    // and free. It reads as the plate held over x < 1 alone, in as many
    // iterations, and says nothing of it. (Traced, it took 891 iterations
    // for 104, read 0.04% stiff, and said the hold held free slivers.)
    for top in ["0.2", "0.2137", "(+ 0.2 (* 0.05 y))"] {
        let plate = |wall: &str| {
            let src = format!(
                "(def plate (max (move (box 3 1 1) 1.5 0 0.5) (- z {top})))
(study s (shape plate) (elastic (E 200000MPa) (nu 0.3))
  (fix base {wall}) (load tip (force 0 0 -0.01N) (region (- 2.5 x))))"
            );
            let sol = solve(&src, "s", 0.5);
            (sol.max_nodal().unwrap().value, sol.report.iterations, sol.report.warnings)
        };
        let (alone, flush) = (plate("(- x 1)"), plate(&format!("(max (- x 1) (- z {top}))")));
        // (tilted, 4.4e-5 apart: as a region a ten-thousandth past the
        // face reads, its samples split by points)
        assert!((flush.0 / alone.0 - 1.0).abs() < 1e-4 && flush.1 == alone.1 && flush.2 == alone.2, "top {top}: {flush:?} against {alone:?}");
    }
    // So too a round face: a rod 0.6 mm across held over x < 1 and under
    // its own face (7e-5 apart, the region's face read at its points by
    // both fields: by the samples' tangent planes, which a curve's do not
    // share, it was traced, 0.81% stiff in 54 iterations for 19, and said
    // it held 1e-16 mm^3 of free material).
    let rod = |wall: &str| {
        let src = format!(
            "(def rod (max (move (box 4 2 2) 2 0 0) (- (len2 y z) 0.3)))
(study s (shape rod) (elastic (E 200000MPa) (nu 0.3))
  (fix base {wall}) (load tip (force 0 0 -0.01N) (region (- 3.5 x))))"
        );
        let sol = solve(&src, "s", 0.5);
        (sol.max_nodal().unwrap().value, sol.report.iterations, sol.report.warnings)
    };
    let (alone, flush) = (rod("(- x 1)"), rod("(max (- x 1) (- (len2 y z) 0.3))"));
    assert!((flush.0 / alone.0 - 1.0).abs() < 2e-4 && flush.1 <= alone.1 && flush.2 == alone.2, "{flush:?} against {alone:?}");
}

#[test]
fn a_hold_edge_holds_at_any_tilt_offset_or_curve() {
    let tip = |wall: &str, res: f64| solve(&held_bar(wall), "s", res).max_nodal().unwrap().value;
    // Its plane x = 0.998 (a sample left under SNAP free, all held) tilted
    // by t through y: the free sample beside took its edge on the plane
    // moved to the cube's farthest point along its normal, for any tilt a
    // corner: at 1e-8 and 1e-4 the hold vanished, 10,852 mm against 0.2444.
    let level = tip("(- x 0.998)", 0.5);
    for t in [1e-12, 1e-10, 1e-8, 1e-6, 1e-4] {
        let v = tip(&format!("(- (+ x (* {t} y)) 0.998)"), 0.5);
        assert!((v / level - 1.0).abs() < 1e-3, "tilt {t}: {v} against {level}");
    }
    // Its edge at offsets across a 0.5 mm cell: within 0.5% of the 0.125 mm
    // cells' reading (it read 0.11 to 0.13% off; with a penalty, 0.8 to 4%).
    for (o, fine) in [(0.9, 0.248573), (1.05, 0.242746), (1.1, 0.240809), (1.2, 0.23701), (1.3, 0.233248), (1.4, 0.229526)] {
        let v = tip(&format!("(- x {o})"), 0.5);
        assert!((v / fine - 1.0).abs() < 5e-3, "offset {o}: {v} against {fine}");
    }
    // A curved edge, a cylinder about y of radius 1.37: its three grids
    // converge (0.239315, 0.239632, 0.239735; the penalty's 0.2278 on 0.5).
    let r: Vec<f64> = [0.5, 0.25].iter().map(|&h| tip("(- (len2 x z) 1.37)", h)).collect();
    assert!((r[0] / 0.239735 - 1.0).abs() < 3e-3 && (r[1] / 0.239735 - 1.0).abs() < 1e-3 && r[0] < r[1], "{r:?}");
}

#[test]
fn a_regions_faces_are_traced_whole_where_they_meet() {
    // The boundary a hold traces, against its exact area, where its faces
    // lie on the samples' faces (a sixteenth of a cell) and meet there: a
    // block's end and side faces at sizes binary holds and cannot, its end
    // face a millionth either side of a cell's plane; a box inside the
    // part (eight corners); two regions joined (an edge turned inward). A
    // sample carries one tangent plane, and at a region's corner it is the
    // other face's: a side face at y = 0.45, which rounding gave to that
    // sample, lost a strip a sample wide along the end face and one along
    // the edge its side faces meet at (11.701 mm^2 of 11.8825, the reading
    // 0.58% off across the plane), and on a cell's plane the cell beside it
    // was not asked.
    // And where they lie inside a sample: a sample holding several of a
    // region's faces carries each (one plane a sample, a face ran on past
    // the others and the second of two was not traced: a block's side faces
    // at a third of a millimetre read 1.6% of its boundary short on one side
    // of a cell's plane and 0.5% over on the other, and a box inside the
    // part, 2.7% over).
    let area = |wall: &str, res: f64| {
        let src = format!(
            "(def bar (move (box 4 4 4) 2 0 0))\n(study s (shape bar) (elastic (E 1000MPa) (nu 0)) (fix wall {wall}) (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))))"
        );
        solve(&src, "s", res).report.held_edge
    };
    let mut cases: Vec<(String, f64)> = Vec::new();
    for side in [0.25, 0.35, 0.45, 0.7, 1.1] {
        for end in [1.2 - 1e-6, 1.2 + 1e-6, 1.3] {
            cases.push((format!("(max (- x {end}) (- y {side}) (- z {side}))"), 2.0 * end * (2.0 + side) + (2.0 + side) * (2.0 + side)));
        }
    }
    cases.push(("(max (- (abs (- x 2)) 0.7) (- (abs y) 0.35) (- (abs z) 0.45))".into(), 5.74));
    cases.push(("(max (- (abs (- x 2)) 0.7) (- (abs (- y 0.05)) 0.35) (- (abs z) 0.45))".into(), 5.74));
    cases.push(("(min (- x 0.7) (max (- x 1.3) (- y 0.35)))".into(), 18.4));
    for side in [1.0 / 3.0, 0.46] {
        for end in [1.2 - 1e-6, 1.2 + 1e-6, 1.23] {
            cases.push((format!("(max (- x {end}) (- y {side}) (- z {side}))"), 2.0 * end * (2.0 + side) + (2.0 + side) * (2.0 + side)));
        }
    }
    cases.push(("(max (- (abs (- x 2)) 0.7) (- (abs y) 0.3333333333333333) (- (abs z) 0.45))".into(), 2.0 * (1.4 / 1.5 + 1.4 * 0.9 + 0.9 / 1.5)));
    cases.push(("(min (- x 0.71) (max (- x 1.33) (- y 0.37)))".into(), 4.0 * (1.63 + 0.62 + 2.37)));
    for (wall, exact) in &cases {
        let a = area(wall, 0.2);
        assert!((a / exact - 1.0).abs() < 1e-5, "{wall} on 0.2 mm cells: {a} mm^2 traced of {exact}");
    }
    for k in [8, 10, 15, 16] {
        let (wall, exact) = &cases[k];
        let a = area(wall, 0.1);
        assert!((a / exact - 1.0).abs() < 1e-5, "{wall} on 0.1 mm cells: {a} mm^2 traced of {exact}");
    }
}

#[test]
fn a_face_held_again_inside_a_clamp_holds_nothing_more() {
    // A clamp over x < 0.3 (its edge inside the first 0.5 mm cells) and the
    // end face x = 0, which it holds, held again by a face fix: the face
    // was held on its own axes beside the clamp's, and the bar read 2.6%
    // stiff (0.26546 against 0.27253 mm); on y only, 0.2%.
    let tip = |extra: &str| {
        let src = format!(
            "(def bar (move (box 20 2 2) 10 0 0))\n{}",
            held_bar("(- x 0.3)").replace("(shape (move (box 20 2 2) 10 0 0))", "(shape bar)").replace("(load tip", &format!("{extra} (load tip"))
        );
        solve(&src, "s", 0.5).max_nodal().unwrap().value
    };
    let alone = tip("");
    for extra in ["(fix face (face bar (normal -1 0 0)))", "(fix face (face bar (normal -1 0 0)) (dofs y))"] {
        let v = tip(extra);
        assert!((v / alone - 1.0).abs() < 1e-9, "{extra}: {v} against {alone}");
    }
}

#[test]
fn a_hold_region_reaching_no_material_is_refused_on_every_grid() {
    // x < -0.005 beside a part from x = 0: within a snap's reach of its face
    // on 1 mm cells, where it held that face (0.2848 mm), as no finer grid's.
    for res in [1.0, 0.5] {
        let mut d = ink::load_str(&format!("(import \"std\")\n{}", held_bar("(- x -0.005)")), "test", Path::new(ROOT)).unwrap();
        let ctx = d.ctx(&[], 0.0).unwrap();
        let e = sim::run_study(&mut d, &ctx, "s", Some(res), 1e-9).err().map(|e| format!("{e:?}")).unwrap_or_default();
        assert!(e.contains("touches no material"), "res {res}: {e}");
    }
    // Reaching into it by 0.0001 it holds.
    assert!(solve(&held_bar("(- x 0.0001)"), "s", 1.0).max_nodal().unwrap().value < 0.3);
    // A second fix's region only in a cell's air, beside its material (the
    // bar's top at y 0.9 in cells of 0.5 to 1.0; the region y > 0.95 by
    // x < 0.3): a sample holding both, and the wall's hold in its cell, let
    // it pass holding nothing.
    let src = "(import \"std\")\n(study s (shape (move (box 20 1.9 2) 10 -0.05 0)) (elastic (E 68900MPa) (nu 0.33))
  (fix wall (- x 0.3)) (fix ghost (max (- x 0.3) (- 0.95 y))) (load tip (force 0 0 -10N) (region (- 19.5 x))))";
    let mut d = ink::load_str(src, "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let e = sim::run_study(&mut d, &ctx, "s", Some(0.5), 1e-9).err().map(|e| format!("{e:?}")).unwrap_or_default();
    assert!(e.contains("fix 'ghost' touches no material"), "{e}");
    // Two fixes holding one material on one axis at different values
    // contradict each other: refused, naming both, a field's or a face's
    // (the field fix won, silently, whatever the order).
    for (extra, other) in [("(fix moved (- x 0.3) (dofs y) (move 0 0.5 0))", "moved"), ("(fix face (face bar (normal -1 0 0)) (move 1 0 0))", "face")] {
        let src = format!(
            "(import \"std\")\n(def bar (move (box 20 2 2) 10 0 0))\n{}",
            held_bar("(- x 0.3)").replace("(shape (move (box 20 2 2) 10 0 0))", "(shape bar)").replace("(load tip", &format!("{extra} (load tip"))
        );
        let mut d = ink::load_str(&src, "test", Path::new(ROOT)).unwrap();
        let ctx = d.ctx(&[], 0.0).unwrap();
        let e = sim::run_study(&mut d, &ctx, "s", Some(0.5), 1e-9).err().map(|e| format!("{e:?}")).unwrap_or_default();
        assert!(e.contains("contradict") && e.contains("'wall'") && e.contains(&format!("'{other}'")), "{extra}: {e}");
    }
}

#[test]
fn a_roller_at_a_corner_holds_its_own_face() {
    // A column 4 x 4 x 20, its base held, its top moved 0.001 down, on a
    // floor held along its normal over y < 0.01: the floor's strip takes in
    // the top face's edge, where `top` moves it along z. There the floor
    // leaves the top face to `top`, holding its own (y = 0), and says so:
    // E A d / L (L 19.49 between the holds) 1.6418 N through the base, none
    // through the floor, on every grid, and alike for a strip 0.5 or 2 mm
    // deep, or written as y^2 < 0.01.
    let column = |floor: &str| {
        format!(
            "(def c (move (box 4 4 20) 2 2 10))
(study s (shape c) (elastic (E 2000MPa) (nu 0) (yield 200MPa))
  (fix base (- z 0.5)) (fix top (- 19.99 z) (dofs z) (move 0 0 -0.001)) {floor} (res 1))"
        )
    };
    for floor in ["(- y 0.01)", "(- y 0.5)", "(- y 2)", "(- (* y y) 0.01)"] {
        for res in [4.0, 2.0, 1.0, 0.5] {
            let s = solve(&column(&format!("(fix floor {floor} (dofs normal) (support))")), "s", res);
            let r = &s.report.reactions;
            assert!((r[0].1[2] / 1.6418 - 1.0).abs() < 0.002 && r[2].1.iter().all(|v| v.abs() < 1e-6), "{floor} at {res}: {r:?}");
            // (said wherever a sample of the top face lies in the strip: on
            // cells of 4 mm none of 0.01 mm's does)
            assert!(res > 1.0 || s.report.warnings.iter().any(|w| w.contains("contact 'floor' leaves") && w.contains("('top')")), "{floor} at {res}: {:?}", s.report.warnings);
        }
    }
    // Two fixes on one face along its normal at different values still
    // contradict each other, on every grid: a lid on the top face, its
    // region a layer on it; the floor's own face pushed; and a pad on the
    // top face, its region crossing it but no face of its own. So too
    // where a region has no gradient on the face (y^2 < 0.01 at y = 0, a
    // lid's |z - 20| < 0.01 at its crest, or the whole part): it neither
    // crosses the face nor lies on it, and holds it.
    let pad = "(fix pad (max (- 19 z) (- (sqrt (+ (* (- x 2) (- x 2)) (* (- y 2) (- y 2)))) 0.8)) (dofs normal) (support))";
    for (floor, a, b) in [
        ("(fix lid (- 19.99 z) (dofs normal) (support))", "lid", "top"),
        ("(fix floor (- y 0.01) (dofs normal) (support)) (fix push (- y 0.01) (dofs y) (move 0 0.001 0))", "floor", "push"),
        (pad, "pad", "top"),
        ("(fix floor (- (* y y) 0.01) (dofs normal) (support)) (fix push (- y 0.01) (dofs y) (move 0 0.001 0))", "floor", "push"),
        ("(fix lid (- (abs (- z 20)) 0.01) (dofs normal) (support))", "lid", "top"),
        ("(fix lid -1 (dofs normal) (support))", "lid", "top"),
    ] {
        for res in [4.0, 1.0, 0.5] {
            let mut d = ink::load_str(&format!("(import \"std\")\n{}", column(floor)), "test", Path::new(ROOT)).unwrap();
            let ctx = d.ctx(&[], 0.0).unwrap();
            let e = sim::run_study(&mut d, &ctx, "s", Some(res), 1e-9).err().map(|e| format!("{e:?}")).unwrap_or_default();
            assert!(e.contains("contradict") && e.contains(&format!("'{a}'")) && e.contains(&format!("'{b}'")), "{floor} at {res}: {e}");
        }
    }
    // Turned 1 degree about x, its top's normal off z: the floor holds the
    // top's edge not at all (held along the rest of that normal, at what
    // `top` leaves of it, it would lock the corner), and reads as a floor
    // stopping short of the top does, to 1e-4.
    let (c, sn) = (1f64.to_radians().cos(), 1f64.to_radians().sin());
    let (yy, zz) = (format!("(+ (* y {c}) (* z {sn}))"), format!("(- (* z {c}) (* y {sn}))"));
    let turned = |floor: &str| {
        format!(
            "(def c (max (- (abs (- x 2)) 2) (- (abs (- {yy} 2)) 2) (- (abs (- {zz} 10)) 10)))
(study s (shape c) (elastic (E 2000MPa) (nu 0) (yield 200MPa))
  (fix base (- {zz} 0.5)) (fix top (- 19.99 {zz}) (dofs z) (move 0 0 -0.001)) {floor} (res 1))"
        )
    };
    let corner = solve(&turned(&format!("(fix floor (- {yy} 0.01) (dofs normal) (support))")), "s", 1.0).report.reactions;
    let short = solve(&turned(&format!("(fix floor (max (- {yy} 0.01) (- {zz} 19.5)) (dofs normal) (support))")), "s", 1.0).report.reactions;
    for (a, b) in corner.iter().zip(&short) {
        assert!((0..3).all(|d| (a.1[d] - b.1[d]).abs() < 1e-4 * 1.6418), "{corner:?} against {short:?}");
    }
}

#[test]
fn a_contact_on_a_sloped_face_another_holds_in_part_holds_the_rest() {
    // A wedge's 45 degree slope held along z by one fix and along its
    // normal by a contact: the contact holds the rest of the normal, x,
    // its terms reading the stress across the slope's own normal, as a face
    // fix holding x there does, to rounding. (Read across x instead, the
    // slope took 13.7716 N for 13.8119.)
    let wedge = |seat: &str| {
        format!(
            "(def w (max (box 10 10 10) (* 0.70710678 (- (+ x z) 3))))
(study s (shape w) (elastic (E 2000MPa) (nu 0.3))
  (fix base (+ z 4.5)) (fix zhold (- 2.99 (+ x z)) (dofs z)) {seat}
  (load push (force 0 0 -50N) (region (+ x 4.5))))"
        )
    };
    let contact = solve(&wedge("(fix seat (face w (normal 1 0 1)) (dofs normal) (support))"), "s", 0.5).report;
    let face = solve(&wedge("(fix seat (face w (normal 1 0 1)) (dofs x))"), "s", 0.5).report;
    assert!((contact.energy / face.energy - 1.0).abs() < 1e-9, "{} against {}", contact.energy, face.energy);
    for (a, b) in contact.reactions.iter().zip(&face.reactions) {
        assert!((0..3).all(|d| (a.1[d] - b.1[d]).abs() < 1e-7 * 50.0), "{:?} against {:?}", contact.reactions, face.reactions);
    }
}

#[test]
fn a_hold_through_cells_its_edge_cuts_holds_only_its_material() {
    // CalculiX (C3D10, three meshes converging to 0.2646 mm at the tip):
    // -0.468%, -0.158% and +0.004% on cells of 0.5, 0.25 and 0.125 mm (the
    // last in the slow set); holding every cell the region's edge cut read
    // 3.1%, 1.5% and 0.7% stiff. Its wall takes the load.
    for (res, within) in [(0.5, 0.006), (0.25, 0.0025)] {
        let s = solve(&turned_cantilever(false), "s", res);
        let tip = s.max_nodal().unwrap().value;
        assert!((tip / 0.26455 - 1.0).abs() < within, "{res}: tip {tip} mm against CalculiX's 0.26455");
        assert!(imbalance(&s) < 1e-6);
    }
    // Its first bending mode: CalculiX 4238.4, 4259.2 and 4270.2 Hz on 0.5,
    // 0.25 and 0.125 mm tets (rising: coarse tets hold fewer of the root's
    // nodes), about 4282 converged (an extrapolation of those three, not a
    // value CalculiX gave); 4284.3 and 4281.7 here on 0.5 and 0.25.
    // (Timoshenko's, 4256.5, takes the root's section free to warp, which a
    // whole held root is not: 0.6% under; Euler and Bernoulli's 4292.1.)
    let f = solve(&turned_cantilever(true), "s", 0.5).report.modes[0].0;
    assert!((f / 4282.0 - 1.0).abs() < 0.003, "{f} Hz");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_hold_through_cells_its_edge_cuts_reads_calculix_on_fine_cells() {
    // The turned cantilever on cells of 0.125 mm: +0.004% of CalculiX's 0.26455.
    let tip = solve(&turned_cantilever(false), "s", 0.125).max_nodal().unwrap().value;
    assert!((tip / 0.26455 - 1.0).abs() < 0.0005, "tip {tip} mm");
}

#[test]
fn a_hold_whose_edge_grazes_the_grid_holds_as_well_at_any_offset() {
    // A hold's edge a millionth of a cell short of a plane of the grid
    // leaves the cells it cuts a free sliver: its field's constant grew as
    // 1 / its free part, and the solve broke down. A sample it leaves under
    // 3% free is held whole, the edge on the face of the free sample beside
    // it (lent to the cell beside it, a sliver's material was cut from its
    // other neighbours: a bar pulled to 10 MPa read 10.47). At every offset
    // the tip follows the free length, (L / 19)^3 of the aligned edge's,
    // in at most 30 iterations against 13 (at 0.99 of a cell, the samples
    // beside the edge 4% free).
    let bar = |at: f64| {
        format!(
            "(body b (shape (move (box 20 2 2) 10 0 0)))
(study s (shape (shape-of b)) (elastic (E 68900MPa) (nu 0.33))
  (fix wall (- x {at})) (load tip (force 0 0 -10N) (region (- 19.5 x))))"
        )
    };
    let aligned = solve(&bar(1.0), "s", 0.5);
    for d in [0.5, 0.99, 1.0 - 1e-6] {
        let s = solve(&bar(1.0 + d * 0.5), "s", 0.5);
        let want = aligned.max_nodal().unwrap().value * ((19.0 - d * 0.5) / 19.0f64).powi(3);
        let got = s.max_nodal().unwrap().value;
        assert!((got / want - 1.0).abs() < 0.003, "{d} of a cell: {got} against {want}");
        assert!(s.report.iterations <= 2 * aligned.report.iterations + 5, "{d} of a cell: {} iterations, {} aligned", s.report.iterations, aligned.report.iterations);
    }
}

/// A metric claim's value and the band its grid and read errors give it.
fn claim_band(c: &ink::claims::Check) -> Result<(f64, f64, f64), String> {
    match &c.evidence {
        ink::claims::Evidence::Metric { measured, grid_error, .. } => {
            let v = measured.value.clone()?;
            let g = grid_error.ok_or("no grid error")?.0;
            // (a read's error is a fraction of the read, as the claim widens it)
            let (lo, hi) = measured.read_error.as_ref().map_or((v, v), |r| r.span(v));
            Ok((v, lo - g, hi + g))
        }
        e => Err(format!("{e:?}")),
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_held_claims_band_covers_its_finest_reading_wherever_the_edge_falls() {
    // Where a hold's edge falls in a cell moves a reading by the grid's
    // own error, which a claim's three
    // grids measure: on the held bar, a thin plate held along an edge and a
    // cylinder's edge, at edges from 0.4 of a 0.5 mm cell free down to on
    // its plane, level and tilted, a claim read on 0.5 and 0.25 mm cells
    // states a band (its grid error and read error) holding the same claim
    // read on 0.125 mm cells, or is unknown saying why; the bar's peak
    // stress too.
    let study = |shape: &str, wall: &str| match shape {
        "plate" => {
            format!("(study s (shape (move (box 20 10 1) 10 0 0)) (elastic (E 68900MPa) (nu 0.33) (yield 240MPa)) (fix wall {wall}) (load tip (force 0 0 -1N) (region (- 19.5 x)))")
        }
        _ => {
            format!("(study s (shape (move (box 20 2 2) 10 0 0)) (elastic (E 68900MPa) (nu 0.33) (yield 240MPa)) (fix wall {wall}) (load tip (force 0 0 -10N) (region (- 19.5 x)))")
        }
    };
    let free = [0.4, 0.2, 0.1, 0.05, 0.025, 2.5e-3, 2.5e-7, 0.0];
    let mut cases: Vec<(&str, String)> = Vec::new();
    for shape in ["bar", "plate"] {
        for tilt in [0.0, 1e-3] {
            cases.extend(free.iter().map(|f| (shape, format!("(- (+ x (* {tilt} y)) {})", 1.0 - f * 0.5))));
        }
    }
    cases.extend(free.iter().map(|f| ("cylinder", format!("(- (len2 x z) {})", 1.5 - f * 0.5))));
    let (mut unknown, mut read) = (0, 0);
    for (shape, wall) in &cases {
        let s = study(shape, wall);
        let claims = if *shape == "bar" { "(require d (max-displacement s < 1000mm))\n(require st (max-stress s < 1000MPa))" } else { "(require d (max-displacement s < 1000mm))" };
        // (the truth: the same claims read on 0.125 mm cells. Cells across a
        // part thinner than two of them are split: a claim at (res 2) on the
        // bar reads on 0.5 mm cells, at (res 0.25) on 0.25 and at (res
        // 0.125) on 0.125; on the 1 mm plate, (res 1), (res 0.5) and (res
        // 0.25) read so.)
        let (coarse, finest) = if *shape == "plate" { ([1.0, 0.5], 0.25) } else { ([2.0, 0.25], 0.125) };
        let fine: Vec<(String, f64)> = checks(&format!("{s} (res {finest}))\n{claims}")).into_iter().map(|(n, c)| (n, claim_band(&c).unwrap().0)).collect();
        for res in coarse {
            for (name, c) in checks(&format!("{s} (res {res}))\n{claims}")) {
                let truth = fine.iter().find(|f| f.0 == name).unwrap().1;
                match claim_band(&c) {
                    Ok((v, lo, hi)) => {
                        read += 1;
                        assert!(lo <= truth && truth <= hi, "{shape} {wall}, (res {res}): {name} {v} in [{lo}, {hi}] against {truth}");
                    }
                    Err(e) => {
                        unknown += 1;
                        assert!(c.status == ink::claims::Status::Unknown, "{shape} {wall}, (res {res}): {name} {:?} ({e}) and not unknown", c.status);
                    }
                }
            }
        }
    }
    eprintln!("{read} claims read, {unknown} unknown");
    // Against the exact answer: a bar of nu 0 pulled 10 MPa (E 1000) moves
    // 0.01 of its free length, held at round edges on cells binary cannot
    // hold (a hold's edge on a sample's face was lost on 0.2 mm cells, and
    // on 0.1 the end face's nodes read outside the part, its peak a cell
    // short); and its heat, 1 mW into its last 0.5
    // mm (k 0.2 W/m/K) with 20 C held: 1.25 K per mm of free length, a
    // quarter of that over the source's half millimetre.
    let (mut unknown, mut read) = (0, 0);
    let mut exact = |what: &str, src: String, truth: f64| {
        for (_, c) in checks(&src) {
            match claim_band(&c) {
                Ok((v, lo, hi)) => {
                    read += 1;
                    assert!(lo <= truth && truth <= hi, "{what}: {v} in [{lo}, {hi}] against {truth}");
                }
                Err(e) => {
                    unknown += 1;
                    assert!(c.status == ink::claims::Status::Unknown, "{what}: {:?} ({e}) and not unknown", c.status);
                }
            }
        }
    };
    for (long, edges) in [(20.0, &[0.05, 0.1, 0.13, 1.0, 1.25][..]), (6.0, &[0.05, 0.1, 0.13][..])] {
        for res in [0.1, 0.2, 0.4] {
            for &x0 in edges {
                let bar = |physics: &str, hold: &str, load: &str, claim: &str| {
                    format!("(def bar (move (box {long} 2 2) {} 0 0))\n(study s (shape bar) {physics} {hold} {load} (res {res}))\n(require c {claim})", long / 2.0)
                };
                let pull = bar(
                    "(elastic (E 1000MPa) (nu 0))",
                    &format!("(fix wall (- x {x0}))"),
                    "(load tip (pressure -10MPa) (region (face bar (normal 1 0 0))))",
                    "(max-displacement s < 1000mm)",
                );
                exact(&format!("bar {long} pulled, held at {x0}, (res {res})"), pull, 0.01 * (long - x0));
                if long == 20.0 {
                    let heat = bar(
                        "(heat (k 0.2W/m/K))",
                        &format!("(temperature cold (- x {x0}) 20degC)"),
                        &format!("(source end (power 1mW) (region (- {} x)))", long - 0.5),
                        "(max-temperature s < 1000degC)",
                    );
                    exact(&format!("bar heated, held at {x0}, (res {res})"), heat, 20.0 + 1.25 * (long - 0.5 - x0) + 1.25 * 0.25);
                }
            }
        }
    }
    eprintln!("against the exact: {read} claims read, {unknown} unknown");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_hold_slab_claims_band_holds_its_truth_or_says_why() {
    // A bar of nu 0 pulled 10 MPa (E 1000), held over a slab of half-width
    // w about x = c: its tip moves 0.01 of its free length past the slab,
    // exactly. A slab thinner than THIN cells, its free material on both
    // sides, is split as a thin feature (one field for both sides held them
    // together: x = 5.0893, w 0.1831 on 0.5 mm cells read 0.007351, its band
    // to 0.007426, against 0.007458); thin still at the floor, refused by
    // name. The holds review's 140 slabs (Python's random, seed 7: c in 5 to
    // 5.5, w 10^U(-2, -0.4)): each a band holding the truth, unknown saying
    // why, or that refusal.
    let slabs = [
        (5.1619, 0.0174),
        (5.3255, 0.0131),
        (5.2679, 0.0385),
        (5.029, 0.0648),
        (5.0187, 0.0494),
        (5.0349, 0.014),
        (5.2123, 0.2104),
        (5.0619, 0.0228),
        (5.3137, 0.3283),
        (5.2886, 0.0431),
        (5.4881, 0.0119),
        (5.4292, 0.0291),
        (5.0721, 0.0154),
        (5.1542, 0.2022),
        (5.0904, 0.0852),
        (5.3195, 0.0394),
        (5.2739, 0.0126),
        (5.0298, 0.0214),
        (5.3402, 0.0483),
        (5.1571, 0.0865),
        (5.2266, 0.0302),
        (5.3972, 0.1313),
        (5.122, 0.083),
        (5.2626, 0.2513),
        (5.3647, 0.0289),
        (5.4901, 0.0154),
        (5.2091, 0.1627),
        (5.076, 0.0606),
        (5.0196, 0.1173),
        (5.3823, 0.0826),
        (5.4377, 0.0318),
        (5.3476, 0.0893),
        (5.2899, 0.0537),
        (5.42, 0.3247),
        (5.237, 0.1155),
        (5.0303, 0.1326),
        (5.3236, 0.3881),
        (5.411, 0.0285),
        (5.1929, 0.1174),
        (5.0113, 0.0548),
        (5.084, 0.0154),
        (5.0295, 0.1695),
        (5.0647, 0.0249),
        (5.1955, 0.2479),
        (5.0403, 0.0523),
        (5.2747, 0.2591),
        (5.4096, 0.2412),
        (5.1392, 0.0462),
        (5.1794, 0.2598),
        (5.4789, 0.0174),
        (5.0881, 0.0235),
        (5.1167, 0.0597),
        (5.2946, 0.0263),
        (5.002, 0.0468),
        (5.1846, 0.0806),
        (5.4765, 0.1273),
        (5.2577, 0.0973),
        (5.3381, 0.0122),
        (5.4498, 0.177),
        (5.4373, 0.1891),
        (5.1962, 0.0435),
        (5.0518, 0.1035),
        (5.0311, 0.0128),
        (5.1044, 0.0182),
        (5.17, 0.0121),
        (5.0001, 0.0175),
        (5.0507, 0.0382),
        (5.0128, 0.2506),
        (5.307, 0.0173),
        (5.1261, 0.036),
        (5.1821, 0.0157),
        (5.4245, 0.3881),
        (5.233, 0.0594),
        (5.0429, 0.0146),
        (5.1713, 0.0265),
        (5.4144, 0.0181),
        (5.0115, 0.3323),
        (5.2641, 0.0172),
        (5.2716, 0.011),
        (5.2641, 0.3678),
        (5.4317, 0.13),
        (5.1306, 0.0386),
        (5.0835, 0.1718),
        (5.2663, 0.1764),
        (5.1648, 0.0227),
        (5.4058, 0.3766),
        (5.4263, 0.1949),
        (5.4092, 0.1527),
        (5.1134, 0.0673),
        (5.1778, 0.0111),
        (5.014, 0.028),
        (5.1296, 0.1282),
        (5.4783, 0.0519),
        (5.4685, 0.3809),
        (5.4775, 0.0383),
        (5.1102, 0.0231),
        (5.0984, 0.0212),
        (5.312, 0.2757),
        (5.4202, 0.0585),
        (5.3265, 0.1903),
        (5.0424, 0.114),
        (5.4549, 0.1785),
        (5.3751, 0.0582),
        (5.0893, 0.1831),
        (5.1663, 0.1911),
        (5.4858, 0.043),
        (5.2007, 0.3272),
        (5.3624, 0.0187),
        (5.0635, 0.0175),
        (5.4524, 0.1952),
        (5.0731, 0.2101),
        (5.4902, 0.1126),
        (5.1752, 0.0755),
        (5.0655, 0.0105),
        (5.4854, 0.1095),
        (5.2633, 0.3117),
        (5.2169, 0.2482),
        (5.4131, 0.0218),
        (5.1259, 0.0294),
        (5.1203, 0.0868),
        (5.1297, 0.0468),
        (5.0655, 0.2858),
        (5.1769, 0.0541),
        (5.2917, 0.2798),
        (5.2103, 0.294),
        (5.2508, 0.0709),
        (5.2618, 0.0107),
        (5.2201, 0.0196),
        (5.002, 0.19),
        (5.0862, 0.0572),
        (5.3626, 0.0777),
        (5.163, 0.0675),
        (5.2777, 0.1798),
        (5.0531, 0.0788),
        (5.1242, 0.0277),
        (5.3861, 0.0649),
        (5.2809, 0.1644),
        (5.4562, 0.0512),
        (5.3063, 0.0644),
        (5.2561, 0.1283),
    ];
    let (mut read, mut unknown, mut refused) = (0, 0, 0);
    for (c, w) in slabs {
        let src = format!(
            "(def bar (move (box 20 4 4) 10 0 0))\n(study s (shape bar) (elastic (E 1000MPa) (nu 0)) (fix wall (- (abs (- x {c})) {w})) (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))) (res 0.5))\n(require d (max-displacement s < 1000mm))"
        );
        let truth = 0.01 * (20.0 - c - w);
        let mut v = toleranced(&src, &[]).unwrap();
        match ink::claims::run_asserts(&mut v, None, false, ink::claims::ClaimOptions::default()) {
            Err(e) => {
                refused += 1;
                assert!(e.message.contains("holds a slab of the part thinner than its cells resolve"), "c {c}, w {w}: {}", e.message);
            }
            Ok(r) => match claim_band(&r.asserts[0].check) {
                Ok((x, lo, hi)) => {
                    read += 1;
                    assert!(lo <= truth && truth <= hi, "c {c}, w {w}: {x} in [{lo}, {hi}] against {truth}");
                }
                Err(e) => {
                    unknown += 1;
                    assert!(r.asserts[0].check.status == ink::claims::Status::Unknown, "c {c}, w {w}: {e}");
                }
            },
        }
    }
    eprintln!("{read} read, {unknown} unknown, {refused} refused");
}

#[test]
fn a_cell_two_regions_hold_is_free_outside_both() {
    // A plate one cell thick held at a wall and on a plane through it: the
    // cells both reach are free outside both regions. (Each region's held
    // part taken from the cell took what they share twice: NaN.) Between
    // plane stress's 1.29 mm and plane strain's 1.18 (the plane holds a
    // quarter of the thickness across), and the same as with the plane
    // stopped at the wall: what the wall holds, the plane holds nothing
    // more of. (The plane's held material kept its stiffness along x and z
    // where the wall held them: 0.7% stiff.)
    let src = |plane: &str| {
        format!(
            "(body p (shape (move (box 40 1 10) 20 0.5 0)))
(study plate (shape (shape-of p)) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (- x 0.5)) (fix plane {plane} (dofs y))
  (load tip (force 0 0 -10N) (region (- 39.5 x))))"
        )
    };
    let s = solve(&src("(- y 0.25)"), "plate", 1.0);
    let tip = s.max_nodal().unwrap().value;
    let apart = solve(&src("(max (- y 0.25) (- 0.5 x))"), "plate", 1.0).max_nodal().unwrap().value;
    assert!((1.18..1.29).contains(&tip) && (tip / apart - 1.0).abs() < 1e-8 && imbalance(&s) < 1e-6, "{tip}, {apart} with the plane apart");
}

#[test]
fn a_clamp_holds_the_same_held_again_or_as_three_fixes_of_one_axis() {
    // A 20 x 2 x 2 cantilever clamped below x = -8.75 (mid-cell on 0.5 mm
    // cells), 1 N down at its tip: held again along y over the same region,
    // or as three fixes of one axis each, it is the one hold. (Each fix's
    // held material kept its stiffness on the axes it left free, another
    // fix holding them or not: 1.5% and 1.8% stiff.)
    let tip = |fixes: &str| {
        let src = format!("(def bar (box 20 2 2))\n(study s (shape bar) (elastic (E 1000MPa) (nu 0.3)) {fixes} (load tip (force 0 0 -1N) (region (- 9.5 x))))");
        solve(&src, "s", 0.5).max_nodal().unwrap().value
    };
    let alone = tip("(fix clamp (+ x 8.75))");
    for fixes in ["(fix clamp (+ x 8.75)) (fix again (+ x 8.75) (dofs y))", "(fix cx (+ x 8.75) (dofs x)) (fix cy (+ x 8.75) (dofs y)) (fix cz (+ x 8.75) (dofs z))"] {
        let t = tip(fixes);
        assert!((t / alone - 1.0).abs() < 1e-8, "{fixes}: {t}, alone {alone}");
    }
}

#[test]
fn two_fixes_turning_on_one_pivot_leave_their_study_unsettled() {
    // A 20 x 2 x 2 bar held below x = -8.75 along x and along y on a pivot
    // about z through its end, along z on its underside beyond, and along y
    // beyond x = 9, 1 N down at its middle. As one fix of both axes on that
    // pivot it reads 0.0585 mm. As two fixes, one of each axis on the pivot,
    // the held material's stiffness between their axes is not where their
    // free motions move it: it read 0.0302 (their turns 0.0011 and 0.0006 rad
    // against 0.0092). So held, every claim on the study is unknown, naming
    // the two fixes. (Held rigidly, `(rigid)`: a distributing pivot holds no
    // material, only its mean motion, and two of them are two supports.)
    let src = |fixes: &str| {
        format!(
            "(def bar (box 20 2 2))\n(study s (shape bar) (elastic (E 1000MPa) (nu 0.3)) {fixes}
  (fix floor (max (+ z 0.75) (- -8.75 x)) (dofs z) (support)) (fix guide (- 9 x) (dofs y) (support)) (load mid (force 0 -1N 0) (region (- (abs x) 0.5))) (res 0.5))
(require stiff (max-displacement s < 1))"
        )
    };
    let (_, one) = &checks(&src("(fix pxy (+ x 8.75) (dofs x y) (pivot (at -10 0 0) (axis 0 0 1)) (rigid))"))[0];
    assert!(one.status == ink::claims::Status::Holds && (measured(one).unwrap() / 0.0585 - 1.0).abs() < 0.01, "{one:?}");
    let (_, two) =
        &checks(&src("(fix px (+ x 8.75) (dofs x) (pivot (at -10 0 0) (axis 0 0 1)) (rigid)) (fix py (+ x 8.75) (dofs y) (pivot (at -10 0 0) (axis 0 0 1)) (rigid))"))[0];
    let said = format!("{two:?}");
    assert!(two.status == ink::claims::Status::Unknown && said.contains("fixes 'px' and 'py' hold different axes of the same material"), "{two:?}");
}

#[test]
fn a_sum_of_runs_reads_their_sum_inside_held_material() {
    // A bar clamped at its end and held along y alone over its next 0.75
    // mm, pulled along x by one run and pushed down by another: their sum
    // (once the one, twice the other) reads the runs' sum at every node.
    // (A node inside held material read the first run's.)
    let study = |name: &str, force: &str| {
        format!("(study {name} (shape bar) (elastic (E 1000MPa) (nu 0.3)) (fix end (+ x 9.5)) (fix guide (+ x 8.75) (dofs y)) (load tip (force {force}) (region (- 9.5 x))))")
    };
    let src = format!("(def bar (box 20 2 2))\n{}\n{}", study("a", "1N 0 0"), study("b", "0 0 -1N"));
    let (a, b) = (solve(&src, "a", 0.5), solve(&src, "b", 0.5));
    let mut sum = a.clone();
    assert!(sum.set_sum(&[(1.0, &a), (2.0, &b)]));
    for n in 0..a.u.len() / 3 {
        let (s, x, y) = (sum.nodal(n), a.nodal(n), b.nodal(n));
        assert!((0..3).all(|d| (s[d] - x[d] - 2.0 * y[d]).abs() <= 1e-12 * (1.0 + s[d].abs())), "node {n}: {s:?}, the runs {x:?} and {y:?}");
    }
    for (combined, one, two) in [
        (&sum.report.applied_moments, &a.report.applied_moments, &b.report.applied_moments),
        (&sum.report.reaction_moments, &a.report.reaction_moments, &b.report.reaction_moments),
    ] {
        assert!(!combined.is_empty());
        for ((s, x), y) in combined.iter().zip(one).zip(two) {
            assert!((0..3).all(|d| (s[d] - x[d] - 2.0 * y[d]).abs() < 1e-10), "{s:?} against {x:?} + 2 * {y:?}");
        }
    }
}

#[test]
fn a_held_temperature_conducts_its_flux_and_reads_as_held() {
    // A 10 mm cube turned 20 degrees, held at 100 C below z = 2.1 and at 20
    // above 7.9: k A dT / L = 0.2759 W through it, and nothing read hotter
    // than 100 C: inside held material, the held value, not the free
    // field's extension its nodes carry. (Held over every cell the edges
    // cut, it took 10% more, and read 100.17 C.)
    let src = "(body b (shape (rotate-z (move (box 10 10 10) 0 0 5) 20)))
(study hot (shape (shape-of b)) (heat (k 0.2W/m/K))
  (temperature base (- z 2.1) 100degC) (temperature top (- 7.9 z) 20degC))";
    let s = solve(src, "hot", 0.5);
    let q = s.report.reactions.iter().find(|r| r.0 == "base").unwrap().1[0];
    let want = 2e-4 * 100.0 * 80.0 / 5.8;
    assert!((q / want - 1.0).abs() < 0.005, "{q} W against {want}");
    let peak = s.max_nodal().unwrap().value;
    assert!((peak - 100.0).abs() < 1e-9, "{peak} C");
    // A point inside it reads the same, as a probe and a picture read it
    // (beside the region's edge, 101.10 and 101.24 C).
    for p in [[0.0, 0.0, 0.25], [0.0, 0.0, 1.9], [4.0, -3.0, 1.95]] {
        let t = s.value_at(p).unwrap()[0];
        assert!((t - 100.0).abs() < 1e-9, "{t} C at {p:?}");
    }
}

#[test]
fn a_claim_keeps_its_pressure_faces_on_comparison_grids() {
    // These 6 mm sides are on the reference lattice, but cut the 2 mm
    // grid. Both solves must receive the same 36 N end load.
    let src = r#"(import "std")
(def bar (box 16 6 6))
(study pull (shape bar) (elastic (E 10000MPa) (nu 0)) (fix foot (+ x 7))
  (load end (pressure -1MPa) (region (face bar (normal 1 0 0)))) (res 1))"#;
    let mut doc = ink::load_str(src, "test", Path::new(ROOT)).unwrap();
    let ctx = doc.ctx(&[], 0.0).unwrap();
    for h in [1.0, 2.0, 0.5, 0.3] {
        let (setup, _, _) = sim::study_setup(&mut doc, &ctx, "pull", Some(h), 1e-9).unwrap();
        let sol = sim::solve(&setup, &ctx).unwrap();
        assert!((sol.report.applied[0].1[0] - 36.0).abs() < 1e-10, "cells {h}: {:?}", sol.report.applied);
    }
}

#[test]
fn cantilever_converges_to_beam_theory() {
    // Built in over x < -48, 10 N spread over x > 48: at the tip
    // q / (6 E I) int_96^98 s^2 (3 98 - s) ds = 3.707 mm bending, + 0.023
    // shear = 3.730. The clamp's and the load's edges lie on the cells'
    // planes at both sizes, and the elements bend exactly: 2 mm cells (the
    // faces at +-5 cutting theirs) have converged too, 3.737 against 3.741
    // on 1 mm. (With the clamp's edge
    // inside a cell, the cell is held as far as its samples are, moving the
    // answer by up to half a cell's share of the length: 3.95 at 2 mm, 3.89
    // at 1.)
    let src = "(body beam (shape (box 100 10 10)))
      (study tip (shape (shape-of beam)) (elastic (E 1000MPa) (nu 0))
        (fix wall (+ x 48)) (load tip (force 0 0 -10N) (region (- 48 x))))";
    let coarse = solve(src, "tip", 2.0).max_nodal().unwrap().value;
    let fine = solve(src, "tip", 1.0).max_nodal().unwrap().value;
    assert!((fine - 3.73).abs() / 3.73 < 0.01, "fine {fine}");
    assert!((coarse - 3.73).abs() / 3.73 < 0.005 && (coarse - fine).abs() < 2e-3 * fine, "coarse {coarse}, fine {fine}");
}

/// Every vertex of a part's surface mesh has a displacement near it, and
/// neighbours' agree (`solve --render` draws the part deformed by it): a
/// round rod's surface runs through cells whose material rounds to none,
/// where the solution itself has no value. One vertex left unmoved while
/// its neighbours move folded the triangles between them, drawn as specks
/// and stray outlines.
#[test]
fn every_surface_point_moves_with_its_neighbours() {
    let src = "(body rod (shape (cylinder-x 4.3 100)))
      (study tip (shape (shape-of rod)) (elastic (E 1000MPa) (nu 0.3))
        (fix wall (+ x 48)) (load tip (force 0 0 -10N) (region (- 48 x))))";
    let s = solve(src, "tip", 1.0);
    let mut d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let it = d.items("rod").unwrap().remove(0);
    let f = d.placed_field(&it);
    let f = d.clip_reach(f, &ctx).0;
    let m = ink::mesh::mesh(&std::sync::Arc::new(ink::tape::Tape::new(&d.g, f)), &ctx, [-51.0, -5.0, -5.0], [51.0, 5.0, 5.0], 0.5);
    assert!(m.verts.iter().any(|&p| s.value_at(p).is_none()), "some surface points lie in no cell");
    let u: Vec<Vec<f64>> = m.verts.iter().map(|&p| s.value_near(p).expect("every surface point has a displacement")).collect();
    let peak = s.max_nodal().unwrap().value;
    let jump = m
        .tris
        .iter()
        .flat_map(|t| [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])])
        .map(|(a, b)| (0..3).map(|k| (u[a as usize][k] - u[b as usize][k]).abs()).fold(0.0, f64::max))
        .fold(0.0, f64::max);
    assert!(jump < 0.02 * peak, "neighbouring surface points move {jump} apart, against {peak} at most");
}

#[test]
fn a_slot_narrower_than_a_cell_parts_its_sides() {
    // A 100 x 10 x 10 cantilever slotted along its mid-plane (0.6 wide, |z|
    // < 0.3) over its outer half, at 1 mm cells: the two halves bend apart.
    // Beam theory: 5.733 mm bending (I = 833 inboard, 2 x 10 x 4.7^3 / 12 =
    // 173 outboard) + 0.025 shear = 5.757; bonded across the slot, as
    // cells that held both sides of it did, it was 3.876 (the solid
    // beam's 3.881).
    let src = "(body beam (shape (max (box 100 10 10) (neg (max (- (abs z) 0.3) (neg x))))))
      (study tip (shape (shape-of beam)) (elastic (E 1000MPa) (nu 0))
        (fix wall (+ x 49.5)) (load tip (force 0 0 -10N) (region (- 49 x))))";
    let s = solve(src, "tip", 1.0);
    let (upper, lower) = (-s.value_at([49.5, 0.0, 2.5]).unwrap()[2], -s.value_at([49.5, 0.0, -2.5]).unwrap()[2]);
    assert!((upper - 5.757).abs() < 0.03 * 5.757 && (lower - 5.757).abs() < 0.03 * 5.757, "tip {upper} {lower} vs 5.757");
    // Inboard of the slot the beam is whole.
    assert!(s.value_at([-10.0, 0.0, 0.0]).is_some());
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_cut_flush_with_a_face_leaves_no_film_across_its_mouth() {
    // That slot cut to the beam's end, its end flush with the beam's: in
    // the slot the field is |x - 50|, never negative. On 0.7 mm cells, off
    // that face, the samples straddling it read a third material, a film
    // across the slot's mouth bonding the halves at the tip: 4.55 mm
    // where they bend 5.76.
    let src = "(body beam (shape (max (box 100 10 10) (neg (max (- (abs z) 0.3) (neg x) (- x 50))))))
      (study tip (shape (shape-of beam)) (elastic (E 1000MPa) (nu 0))
        (fix wall (+ x 49.5)) (load tip (force 0 0 -10N) (region (- 49 x))))";
    let s = solve(src, "tip", 0.7);
    let (upper, lower) = (-s.value_at([49.5, 0.0, 2.5]).unwrap()[2], -s.value_at([49.5, 0.0, -2.5]).unwrap()[2]);
    assert!((upper - 5.757).abs() < 0.03 * 5.757 && (lower - 5.757).abs() < 0.03 * 5.757, "tip {upper} {lower} vs 5.757");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn rotated_part_gives_the_same_answer() {
    // The grid is axis-aligned; the physics must not care. The bar's axis is
    // x rotated 25 deg about y, then 40 deg about z: (0.6942, 0.5825, -0.4226).
    let src = "(body bar (shape (rotate-z (rotate-y (box 100 10 10) 25) 40)))
      (study pull (shape (shape-of bar)) (elastic (E 1000MPa) (nu 0))
        (fix wall (+ (rotate-z (rotate-y x 25) 40) 49.5))
        (load pull (force 694.2N 582.5N -422.6N) (region (- 49 (rotate-z (rotate-y x 25) 40)))))";
    let s = solve(src, "pull", 1.0);
    let (_, vm, _) = s.field_at([0.0, 0.0, 0.0]).unwrap();
    assert!((vm - 10.0).abs() < 0.25, "stress in rotated bar {vm}");
    assert!(imbalance(&s) < 1e-6);
}

#[test]
fn faces_at_multiples_of_the_cell_size_cut_no_cells() {
    // The grid's lines lie at multiples of the cell size, wherever the
    // bounds' search stops: an 8 x 2 x 4 box at 1 mm cells is 64 whole
    // cells. Moved a quarter cell along x, only the cells along its two x
    // faces are cut (it was placed from the bounds, 0.06 mm off every face
    // of a plate: 13,361 of its 30,854 cells cut).
    for (dx, full, cut) in [(0.0, 64, 0), (0.25, 56, 16)] {
        let src = format!(
            "(body b (shape (move (box 8 2 4) {dx} 0 0)))
          (study s (shape (shape-of b)) (elastic (E 1000MPa) (nu 0.3))
            (fix wall (+ x 3)) (load tip (force 0 0 -1N) (region (- 3 x))))"
        );
        let r = solve(&src, "s", 1.0).report;
        assert_eq!((r.cells_full, r.cells_cut), (full, cut), "moved {dx}");
    }
}

#[test]
fn thermal_boundaries_do_not_exclude_heat_flux_cells() {
    // The surface region covers the entire thin plate. It exchanges heat;
    // it is not a clamp with an elastic stress singularity to exclude.
    for boundary in ["(convection air -1 (coefficient 10W/m2/K) (ambient 25degC))", "(temperature cold (+ x 3.5) 25degC)"] {
        let src = format!(
            "(body plate (shape (box 8 8 2)))
          (study hot (shape (shape-of plate)) (heat (k 20W/m/K))
            {boundary} (source chip (power 1W) (region (- 2 x))))"
        );
        let s = solve(&src, "hot", 1.0);
        assert!(s.max_cell_scalar().unwrap().value > 0.0);
        assert!(s.cell_scalar_percentiles(&[0.5, 0.9, 0.99]).is_ok());
        assert!(s.max_near_supports().is_none());
    }
}

#[test]
fn heat_conducts_linearly() {
    // Ends held at 0 and 100 across 99 mm: T(x) = 100 (x + 49.5) / 99,
    // heat flow k A dT/dx = 0.2e-3 * 100 * 100/99 W.
    let src = "(body rod (shape (box 100 10 10)))
      (study hot (shape (shape-of rod)) (heat (k 0.2W/m/K))
        (temperature cold (+ x 49.5) 0) (temperature warm (- 49.5 x) 100degC))";
    let s = solve(src, "hot", 1.0);
    for x in [-30.0, 0.0, 25.0] {
        let t = s.value_at([x, 1.0, -2.0]).unwrap()[0];
        let want = 100.0 * (x + 49.5) / 99.0;
        assert!((t - want).abs() < 0.3, "T({x}) = {t}, want {want}");
    }
    let flow = 0.2e-3 * 100.0 * 100.0 / 99.0;
    let r = &s.report.reactions;
    assert!((r[0].1[0] + flow).abs() / flow < 0.02, "heat out {:?} vs {flow}", r[0].1);
    assert!((r[1].1[0] - flow).abs() / flow < 0.02, "heat in {:?} vs {flow}", r[1].1);
}

#[test]
fn floating_pieces_are_dropped_not_fatal() {
    let src = "(body two (shape (min (box 20 10 10) (move (sphere 3) 40 0 0))))
      (study s (shape (shape-of two)) (elastic (E 1000MPa) (nu 0.3))
        (fix wall (+ x 9.5)) (load push (force 0 0 -1N) (region (- 9 x))))";
    let s = solve(src, "s", 1.0);
    assert!(s.report.dropped_volume > 100.0, "sphere should be dropped: {:?}", s.report.dropped_volume);
}

#[test]
fn named_face_supports_and_loads() {
    // The bar-in-tension check, with supports and loads found by name and
    // normal instead of hand-built regions.
    let src = "(body bar (shape (box 100 10 10)))
      (study pull (shape (shape-of bar)) (elastic (E 1000MPa) (nu 0))
        (fix wall (face bar (normal -1 0 0) (depth 0.5)))
        (load pull (force 1000N 0 0) (region (face bar (normal 1 0 0) (depth 1)))))";
    let s = solve(src, "pull", 1.0);
    let (_, vm, _) = s.field_at([0.0, 0.0, 0.0]).unwrap();
    assert!((vm - 10.0).abs() < 0.15, "stress {vm}");
    assert!(imbalance(&s) < 1e-6);
    assert!((s.report.reactions[0].1[0] + 1000.0).abs() < 1e-6);
}

#[test]
fn cantilever_natural_frequency() {
    // Steel 100 x 10 x 10 cantilever: Euler-Bernoulli first bending mode
    // 1.8751^2 / (2 pi) sqrt(E I / (rho A L^4)) = 823.6 Hz (L = 99.5 mm).
    // Its density is the body's, as a catalog part gives it (g/cm^3, read
    // as g/mm^3 like every density): a study reads the same unit.
    let src = "(body beam (shape (box 100 10 10)) (attr density 7.85g/cm3))
      (study vib (shape (shape-of beam)) (elastic (E 200000MPa) (nu 0.3) (density (by-body beam density)))
        (fix wall (+ x 49.5)) (modes 2))";
    // A density in another unit, unsaid, is refused, not read a million
    // times off.
    let bare = src.replace("(by-body beam density)", "7850");
    let e = format!("{:?}", ink::load_str(&format!("(import \"std\")\n{bare}"), "test", Path::new(ROOT)).err().expect("7850 g/mm^3 refused"));
    assert!(words(&e).contains(&"7850"), "{e}");
    let s = solve(src, "vib", 1.5);
    let f: Vec<f64> = s.report.modes.iter().map(|m| m.0).collect();
    assert!((f[0] - 823.6).abs() / 823.6 < 0.04, "first mode {f:?}");
    assert!((f[1] - f[0]).abs() / f[0] < 0.01, "square section: a degenerate pair {f:?}");
    assert!((s.report.mass_kg - 0.0785).abs() < 0.0005, "mass {}", s.report.mass_kg);
    // Bending about y or z: each mode moves the tip sideways (y or z) and
    // rotates about the other axis; no mode is a twist about x.
    for (_, m) in &s.report.modes {
        assert!(m[1] + m[2] > 0.3 && m[4] + m[5] > 0.2 && m[3] < 0.01, "fractions {m:?}");
    }
    // Probing mode 2: its own shape, 1 mm peak.
    let m2 = s.mode(2).unwrap();
    let tip = m2.value_at([49.0, 0.0, 0.0]).unwrap();
    assert!((tip[1].hypot(tip[2]) - 1.0).abs() < 0.1, "mode 2 tip {tip:?}");
    assert!(s.mode(3).is_none());
}

/// A steel bar 200 x 10 x 10 (or a study of it), its modes asked with
/// `rest` for supports.
fn steel_bar(rest: &str, modes: usize) -> String {
    format!(
        "(body bar (shape (box 200 10 10)))
(study s (shape (shape-of bar)) (elastic (E 200000MPa) (nu 0.3) (density 7.85g/cm3)) {rest} (modes {modes}))"
    )
}

fn hz(s: &sim::Solution) -> Vec<f64> {
    s.report.modes.iter().map(|m| m.0).collect()
}

#[test]
fn a_free_part_has_its_elastic_modes() {
    // A bar on nothing: its six rigid motions are no modes, said, and its
    // modes are the free-free beam's. Timoshenko's (shear and rotary
    // inertia, k = 10 (1 + nu) / (12 + 11 nu)), by the first-order
    // correction 1 / sqrt(1 + (beta r / L)^2 (1 + E / k G)) of Euler and
    // Bernoulli's (beta L 4.7300, 7.8532, 10.9956: 1297.0, 3575.4, 7009.4 Hz):
    // bending pairs 1285, 3486 and 6676 Hz; torsion, St-Venant's (J 0.1406
    // a^4), 7188 Hz. Here 1285.5, 3487.3, 6686.1 and 7147.3 on cells of 1 mm.
    let s = solve(&steel_bar("", 9), "s", 1.0);
    let f = hz(&s);
    for (k, want, within) in [(0, 1285.0, 0.003), (1, 1285.0, 0.003), (2, 3486.0, 0.003), (3, 3486.0, 0.003), (4, 6676.0, 0.003), (5, 6676.0, 0.003), (6, 7188.0, 0.01)] {
        assert!((f[k] / want - 1.0).abs() < within, "mode {}: {f:?} against {want}", k + 1);
    }
    // No mode moves it as a whole: none has a share of a rigid motion.
    assert!(s.report.modes.iter().all(|m| m.1.iter().all(|v| *v < 1e-6)), "{:?}", s.report.modes);
    assert!(s.report.warnings.iter().any(|w| w.contains("nothing holds the part") && w.contains("6 rigid motions, at 0 Hz, left out of the modes")), "{:?}", s.report.warnings);
    // Each mode at its own frequency on a coarser grid: 1284.9 on cells of 2 mm.
    assert!((hz(&solve(&steel_bar("", 2), "s", 2.0))[0] / f[0] - 1.0).abs() < 0.002);
}

#[test]
fn a_partly_held_part_has_its_free_motions_left_out() {
    // On two rollers (z held on lines through its ends, x, y and a turn
    // about z free): bending in z is a simply supported beam's,
    // Timoshenko's 575.6 and 2274 Hz over its 199 mm span (Euler and
    // Bernoulli's 578.0 and 2312); in y it is free, the free bar's 1285.5
    // and 3487.3 (above) to the hertz.
    let rollers = "(fix left (+ x 99) (dofs z) (pivot (at -99.5 0 0) (axis 0 1 0))) (fix right (- 99 x) (dofs z) (pivot (at 99.5 0 0) (axis 0 1 0)))";
    let s = solve(&steel_bar(rollers, 4), "s", 1.0);
    let (f, free) = (hz(&s), hz(&solve(&steel_bar("", 4), "s", 1.0)));
    assert!((f[0] / 575.6 - 1.0).abs() < 0.003 && (f[2] / 2274.0 - 1.0).abs() < 0.003, "{f:?}");
    assert!((f[1] / free[0] - 1.0).abs() < 1e-4 && (f[3] / free[2] - 1.0).abs() < 1e-4, "{f:?} against {free:?}");
    // A simply supported beam's first mode: 8 / pi^2 of its mass along z.
    assert!((s.report.modes[0].1[2] - 0.81).abs() < 0.02, "{:?}", s.report.modes[0]);
    assert!(s.report.warnings.iter().any(|w| w.contains("3 rigid motions")), "{:?}", s.report.warnings);
}

#[test]
fn each_free_piece_keeps_its_own_modes() {
    // Two bars apart, 200 and 100 mm: twelve rigid motions, and each bar's
    // modes as it has alone. With the long one clamped, the short one,
    // which no fix holds, is dropped (as a static study drops it, its claims
    // not settled), not rung free beside it: none of its modes is among the
    // cantilever's (204.9 Hz first), and no rigid motion is left out.
    let two = |rest: &str, modes: usize| {
        format!(
            "(body a (shape (box 200 10 10))) (body b (shape (move (box 100 10 10) 0 30 0))) (assembly ab (instance a a) (instance b b))
(study s (shape (shape-of ab)) (elastic (E 200000MPa) (nu 0.3) (density 7.85g/cm3)) {rest} (modes {modes}))"
        )
    };
    let short = hz(&solve("(body b (shape (box 100 10 10))) (study s (shape (shape-of b)) (elastic (E 200000MPa) (nu 0.3) (density 7.85g/cm3)) (modes 2))", "s", 1.0));
    let long = hz(&solve(&steel_bar("", 4), "s", 1.0));
    let s = solve(&two("", 6), "s", 1.0);
    let f = hz(&s);
    let near = |a: f64, b: f64| (a / b - 1.0).abs() < 1e-4;
    assert!(near(f[0], long[0]) && near(f[2], long[2]) && near(f[4], short[0]) && near(f[5], short[1]), "{f:?}: {long:?} and {short:?}");
    assert!(s.report.warnings.iter().any(|w| w.contains("12 rigid motions")), "{:?}", s.report.warnings);
    let s = solve(&two("(fix wall (+ x 99.5))", 4), "s", 1.0);
    let f = hz(&s);
    assert!((f[0] / 204.9 - 1.0).abs() < 0.005, "{f:?}");
    assert!(!f.iter().any(|&v| near(v, short[0])), "{f:?}: the free bar's {short:?}");
    assert!(s.report.dropped_cells > 0 && s.report.rigid_modes == 0, "{:?}", s.report.warnings);
}

#[test]
fn many_modes_of_free_pieces_keep_clear_of_their_rigid_motions() {
    // Twelve modes of a long bar and a short block apart, their lowest
    // 1285.5 and about 25,000 Hz: as the bar's lower modes converge first,
    // none of the twelve falls back toward its twelve rigid motions (every
    // one well over 0 Hz, moving neither piece as a whole), and the bar's
    // first is its own.
    let src = "(body a (shape (box 200 10 10))) (body b (shape (move (box 20 10 10) 0 30 0))) (assembly ab (instance a a) (instance b b))
(study s (shape (shape-of ab)) (elastic (E 200000MPa) (nu 0.3) (density 7.85g/cm3)) (modes 12))";
    let s = solve(src, "s", 2.0);
    let f = hz(&s);
    assert_eq!(f.len(), 12);
    assert!(f.iter().all(|&v| v > 1000.0), "{f:?}");
    assert!(s.report.modes.iter().all(|m| m.1.iter().all(|v| *v < 1e-6)), "{:?}", s.report.modes);
    let alone = hz(&solve(&steel_bar("", 2), "s", 2.0));
    assert!((f[0] / alone[0] - 1.0).abs() < 1e-4, "{f:?} against {alone:?}");
    assert_eq!(s.report.rigid_modes, 12);
}

#[test]
fn a_free_ring_rings_as_a_ring() {
    // A steel ring, mean radius 50, section 4 x 4, on nothing. Its n = 2
    // modes (thin ring theory, EI / (rho A R^4) its scale): in its plane
    // n (n^2 - 1) / sqrt(n^2 + 1), 995.6 Hz; out of it n (n^2 - 1) /
    // sqrt(n^2 + EI / GJ), 945.7 Hz. Each a pair, split here by the grid:
    // 945.4, 951.1, 999.3 and 1015.1 Hz on cells of 1 mm (942.0, 943.6,
    // 994.4 and 998.8 on 0.5).
    let src = "(study s (shape (- (max (abs (- (sqrt (+ (* x x) (* y y))) 50)) (abs z)) 2)) (elastic (E 200000MPa) (nu 0.3) (density 7.85g/cm3)) (modes 4))";
    let f = hz(&solve(src, "s", 1.0));
    for (k, want) in [(0, 945.7), (1, 945.7), (2, 995.6), (3, 995.6)] {
        assert!((f[k] / want - 1.0).abs() < 0.02, "mode {}: {f:?} against {want}", k + 1);
    }
}

#[test]
fn studies_read_numbers_in_their_units() {
    // Written with units, a study reads the same numbers: a power of ten
    // exactly (1kN, 1GPa), degF to rounding.
    let a = solve(BAR, "pull", 2.0).max_nodal().unwrap().value;
    let b = solve(&BAR.replace("(E 1000MPa)", "(E 1GPa)").replace("(force 1000N 0 0)", "(force 1kN 0 0)"), "pull", 2.0).max_nodal().unwrap().value;
    assert_eq!(a, b);
    let hot = |t: &str| solve(&HOT_BAR.replace("120degC", t), "hot", 2.0).max_nodal().unwrap().value;
    let (c, f) = (hot("120degC"), hot("248degF"));
    assert!((c - f).abs() < 1e-9 * c, "{c} {f}");
    // A number in another quantity's unit is refused, naming the form and
    // the unit written: the kernel's words around them (the quantity it
    // takes) are its own, and its tests'.
    let named = BAR.replace("(force 1000N 0 0)", "(force (neg p) 0 0)");
    for (src, from, to, says) in [
        (BAR, "(force 1000N 0 0)", "(force 1000psi 0 0)", "force psi"),
        (named.as_str(), "(body bar", "(param p 50psi)\n(body bar", "psi 50psi"),
        // Through a let, as through any binding or call: the value carries its unit.
        (BAR, "(force 1000N 0 0)", "(force 0 0 (let ((p 50psi)) p))", "force 50psi"),
        (BAR, "(E 1000MPa)", "(E (- 1000N 1N))", "E"),
        (BAR, "(fix wall", "(gravity 0 0 -9.81lbf) (fix wall", "gravity"),
        (HOT_BAR, "120degC)", "120lbf)", "temperature"),
        (BAR, "(body bar (shape (box 100 10 10)))", "(body bar (shape (box 100 10 10)) (attr density 2in))", "density"),
        (BAR, "(force 1000N 0 0) (region (- 49 x))))", "(force 1000N 0 0) (region (- 49 x))))\n(require a (max-stress pull < 10mm))", "max-stress"),
        // A change (an accuracy, a tolerance) in a unit with a zero of its
        // own would read a point on its scale: 1degF as -17.2 degC.
        (HOT_BAR, "120degC))", "120degC))\n(require a (max-temperature hot < 130degC) (accuracy 1degF))", "degF"),
        (HOT_BAR, "(body bar", "(param w 10 (tol -1K 1K))\n(body bar", "K"),
        // In a sum too, whichever term it is.
        (HOT_BAR, "120degC))", "120degC))\n(require a (max-temperature hot < 130degC) (accuracy (+ 1degC 1degF)))", "degF"),
        (HOT_BAR, "120degC))", "120degC))\n(require a (max-temperature hot < 130degC) (accuracy (+ 1degF 1degC)))", "degF"),
    ] {
        let bad = src.replace(from, to);
        assert_ne!(bad, src, "{to}");
        let e = format!("{:?}", ink::load_str(&format!("(import \"std\")\n{bad}"), "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{to} refused")));
        assert!(says.split(' ').all(|w| words(&e).contains(&w)), "{to}: {e}");
    }
}

const HOT_BAR: &str = "(body bar (shape (box 100 10 10)))
  (study hot (shape (shape-of bar)) (heat (k 50W/m/K))
    (temperature left (+ x 49.5) 120degC) (temperature right (- 49.5 x) 120degC))";

#[test]
fn transient_heat_matches_the_series_solution() {
    // Bar at 20 C, one end suddenly held at 100 C, the rest insulated:
    // T = T1 + (T0 - T1) sum 4/((2n+1) pi) sin((2n+1) pi x / 2L) exp(-((2n+1) pi / 2L)^2 a (/ t 1s)).
    let src = "(body bar (shape (move (box 100 10 10) 50 0 0)))
(study warm (shape (shape-of bar)) (heat (k 0.2W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))
  (temperature hot (- x 0.5) 100degC) (transient (duration 10000s) (steps 100) (initial 20degC)))";
    let s = solve(src, "warm", 2.0);
    let (a, l, t) = (0.2, 99.75, 10000.0);
    for x in [50.0, 99.0] {
        let xs = x - 0.25;
        let series: f64 = (0..200)
            .map(|n| {
                let k = (2 * n + 1) as f64 * std::f64::consts::PI / (2.0 * l);
                4.0 / ((2 * n + 1) as f64 * std::f64::consts::PI) * (k * xs).sin() * (-k * k * a * t).exp()
            })
            .sum();
        let want = 100.0 + (20.0 - 100.0) * series;
        let got = s.value_at([x, 0.0, 0.0]).unwrap()[0];
        assert!((got - want).abs() < 0.02 * (want - 20.0), "T({x}) = {got}, series {want}");
    }
}

#[test]
fn convecting_fin_matches_fin_theory() {
    // 10 x 10 mm steel-like fin (k 20), base held at 100 C, h = 10 W/(m^2 K)
    // to 20 C on the sides and tip: 1D fin with a convecting tip.
    let src = "(body rod (shape (move (box 100 10 10) 50 0 0)))
(study fin (shape (shape-of rod)) (heat (k 20W/m/K))
  (temperature base (- x 0.5) 100degC)
  (convection air (- 0.5 x) (coefficient 10W/m2/K) (ambient 20degC)))";
    let s = solve(src, "fin", 0.5);
    let (k, hc, a, l): (f64, f64, f64, f64) = (20.0, 10.0, 0.01, 0.0995);
    let (p, area) = (4.0 * a, a * a);
    let m = (hc * p / (k * area)).sqrt();
    let r = hc / (m * k);
    let theta = |x: f64| ((m * (l - x)).cosh() + r * (m * (l - x)).sinh()) / ((m * l).cosh() + r * (m * l).sinh());
    let q = (hc * p * k * area).sqrt() * 80.0 * ((m * l).sinh() + r * (m * l).cosh()) / ((m * l).cosh() + r * (m * l).sinh());
    for x in [25.0, 50.0, 99.5] {
        let want = 20.0 + 80.0 * theta((x - 0.5) / 1000.0);
        let got = s.value_at([x, 0.0, 0.0]).unwrap()[0];
        assert!((got - want).abs() < 0.6, "T({x}) = {got}, fin theory {want}");
    }
    let removed = -s.report.reactions.iter().find(|r| r.0 == "air").unwrap().1[0];
    assert!((removed - q).abs() < 0.015 * q, "convected {removed} W, fin theory {q} W");
}

#[test]
fn convection_takes_the_surface_s_area_wherever_the_grid_sits() {
    // A 3 mm pin, 20 long, held at 100 C, gives off h A 100 to 0 C air:
    // A = 2 pi r L + 2 pi r^2, its axis on a corner of the cells, on a
    // centre, or neither. A band of samples under the surface, counted,
    // read the side 15% short on a corner or centre (benchmark 16's
    // pins: its rise moved from 39.4 to 43.1 K with where the grid sat).
    let (r, l, hc) = (1.5, 20.0, 25e-6);
    let area = 2.0 * std::f64::consts::PI * r * (l + r);
    for (x, y) in [(0.0, 0.0), (0.25, 0.25), (0.1, 0.17)] {
        let src = format!(
            "(body pin (shape (move (cylinder-z {r} {l}) {x} {y} 0)))
(study s (shape (shape-of pin)) (heat (k 167W/m/K))
  (temperature held -1degC 100degC)
  (convection air (face pin) (coefficient 25W/m2/K) (ambient 0)))"
        );
        let s = solve(&src, "s", 0.5);
        let removed = -s.report.reactions.iter().find(|r| r.0 == "air").unwrap().1[0];
        let want = hc * area * 100.0;
        assert!((removed / want - 1.0).abs() < 0.005, "axis at ({x}, {y}): convected {removed} W, h A dT {want} W");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_tip_moved_by_a_pin_takes_the_beam_theory_force() {
    // A cantilever built in over x < -48, its end block (x > 49) pushed
    // down 1 mm by a pin through (49.75, 0, 0) along y: the block moves 1 mm
    // along z up to a free turn about the pin. Beam theory over L = 97.75:
    // F = d / (L^3 / 3 E I + L / (k G A)) = 2.6599 N (k = 5/6), the tip
    // turning 1.5 d / L less shear's share. Solved from the held values
    // lifted onto the held nodes: from zero, to 1e-9 of the penalty's forces
    // (a thousand times the elastic ones), the wall's reaction read 0.2%
    // off the pin's at 0.5 mm cells.
    let src = "(body beam (shape (box 100 10 10)))
(study s (shape (shape-of beam)) (elastic (E 1000MPa) (nu 0))
  (fix wall (+ x 48))
  (fix tip (- 49 x) (dofs z) (move 0 0 -1) (pivot (at 49.75 0 0) (axis 0 1 0))))";
    for res in [1.0, 0.5] {
        let s = solve(src, "s", res);
        let r = &s.report.reactions;
        assert!((r[1].1[2] / -2.6599 - 1.0).abs() < 0.002 && (r[0].1[2] / -r[1].1[2] - 1.0).abs() < 2e-4, "res {res}: {r:?}");
        let turn = s.report.free_motions[0].1[0].1;
        assert!((turn / (1.5 / 97.75) - 1.0).abs() < 0.01, "res {res}: turn {turn}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_move_may_be_a_field() {
    // A block's end cells held along x at u = 0.001 x: a uniform strain
    // (stress E / 1000, contraction nu / 1000), to the penalty's error
    // (1e-3 of the strain over a held cell).
    let patch = "(study s (shape (box 10 4 4)) (elastic (E 1000MPa) (nu 0.3))
  (fix ends (- 4 (abs x)) (dofs x) (move (* 0.001 x) 0 0)))";
    let s = solve(patch, "s", 1.0);
    let (u, (sigma, _, _)) = (s.value_at([2.0, 1.0, 1.0]).unwrap(), s.field_at_point([0.0, 0.0, 0.0]).unwrap());
    assert!((sigma[0] - 1.0).abs() < 2e-3 && (u[0] / 0.002 - 1.0).abs() < 2e-3 && (u[1] / -0.0003 - 1.0).abs() < 2e-3, "{sigma:?} {u:?}");
    // Submodelling: the root of a cantilever at 1 mm cells, driven over
    // its cut (x in -21..-20) by a 4 mm solve's displacement, matches the
    // whole beam at 1 mm (0.6%) where the 4 mm solve is 5.5% off.
    let src = "(body beam (shape (box 100 10 10)))
(study coarse (shape (shape-of beam)) (elastic (E 1000MPa) (nu 0.3))
  (fix wall (+ x 48)) (load tip (force 0 0 -10N) (region (- 48 x))) (res 4))
(study fine (shape (shape-of beam)) (elastic (E 1000MPa) (nu 0.3))
  (fix wall (+ x 48)) (load tip (force 0 0 -10N) (region (- 48 x))))
(study s (shape (max (shape-of beam) (+ x 20))) (elastic (E 1000MPa) (nu 0.3))
  (fix wall (+ x 48))
  (fix cut (- -21 x) (move (study-field coarse displacement-x) (study-field coarse displacement-y) (study-field coarse displacement-z))))";
    let stress = |s: &sim::Solution| s.field_at_point([-40.0, 0.0, 4.0]).unwrap().0[0];
    let (fine, sub, coarse) = (stress(&solve(src, "fine", 1.0)), stress(&solve(src, "s", 1.0)), stress(&solve(src, "coarse", 4.0)));
    assert!((sub / fine - 1.0).abs() < 0.01 && (coarse / fine - 1.0).abs() > 0.04, "sub {sub}, fine {fine}, coarse {coarse}");
}

#[test]
fn balanced_loads_square_to_a_free_motion_do_not_hide_one_that_drives_it() {
    // A plate on a pin in its bore, a balanced 1000 N pinch through its far
    // end (along z, square to the pin's turn) and 1.5 N sideways there:
    // 30 N mm about the pin, a mechanism. Measured against each node's
    // load times its motion, the pinch's size let the twist pass as
    // balanced (7.2e-4 of it) and the plate solved "at rest"; against the
    // loads' work along the turn it does nothing. The pinch alone rests.
    let study = |pin: &str, twist: bool| {
        format!(
            "(def plate (subtract (move (box 40 20 5) 5 0 0) (cylinder-z 5 10)))
(study s (shape plate) (elastic (E 70000MPa) (nu 0.3))
  (fix pin (- (len2 x y) 5.5) {pin}) {}
  (load up (force 0 0 1000N) (region (max (- 15 x) (+ z 2)))) (load down (force 0 0 -1000N) (region (max (- 15 x) (- 2 z)))) (res 1))",
            if twist { "(load side (force 0 1.5N 0) (region (- 15 x)))" } else { "" }
        )
    };
    let e = study_error(&study("(dofs normal)", true));
    assert!(e.contains("free to turn about the axis along (0, 0, 1)") && e.contains("(30 N mm"), "{e}");
    let s = solve(&study("(dofs normal)", false), "s", 1.0);
    assert!(s.report.warnings.iter().any(|w| w.contains("held at rest")), "{:?}", s.report.warnings);
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_frictionless_seat_holds_the_surface_it_selects_however_its_region_is_written() {
    // A plate on a frictionless seat round its bore, pulled along the
    // bore's axis. Four regions select the same material: a band round
    // the bore, the same field times 100, and each capped past the plate's
    // faces. Its direction was the region's gradient, so the cap held the
    // bore along z where the plain band let it slide (scaled, a step a
    // cell deep did not reach past the cap: 0.0026 mm). It is the part's
    // surface's normal, for the surface in the region: each holds the
    // bore radially and its faces' rims along z, alike. A band kept inside
    // the faces holds the bore alone, and the plate slides: a mechanism.
    let study = |seat: &str| {
        format!(
            "(def plate (subtract (box 40 20 5) (cylinder-z 5 20)))
(study s (shape plate) (elastic (material steel)) (fix pin {seat} (dofs normal))
  (load axial (force 0 0 100N) (region (- 15 x))) (res 1))"
        )
    };
    let alike: Vec<f64> = ["(- (len2 x y) 5.5)", "(* 100 (- (len2 x y) 5.5))", "(max (- (len2 x y) 5.5) (- (abs z) 2.6))", "(max (* 100 (- (len2 x y) 5.5)) (- (abs z) 2.6))"]
        .iter()
        .map(|seat| {
            let s = solve(&study(seat), "s", 1.0);
            assert!((s.report.reactions[0].1[2] + 100.0).abs() < 1e-6, "{seat}: {:?}", s.report.reactions);
            s.report.energy
        })
        .collect();
    assert!(alike.iter().all(|e| (e / alike[0] - 1.0).abs() < 1e-9), "{alike:?}");
    let e = study_error(&study("(max (- (len2 x y) 5.5) (- (abs z) 2.4))"));
    assert!(e.contains("slide along (0, 0, 1)"), "{e}");
}

#[test]
fn a_moved_study_has_no_compliance() {
    // A cube pressed 0.01 mm (111 N) by a fix: compliance, the loads' work,
    // would read 0 however stiff (and a load against the move, negative),
    // so a claim on it held for any design. It is unknown, saying to read
    // the fix's reaction, and known so before solving (the claims solved
    // the cube at 0.1 mm cells, 30 s each, to say so).
    let src = "(param w 10)
(study s (shape (box w 10 10)) (elastic (E 1000MPa) (nu 0))
  (fix base (+ z 4.5) (dofs z)) (fix press (- 4.5 z) (dofs z) (move 0 0 -0.01)) (res 2))
(study mixed (shape (box 10 10 10)) (elastic (E 1000MPa) (nu 0))
  (fix base (+ z 4.5) (dofs z)) (fix press (- 4.5 z) (dofs z) (move 0 0 -0.01)) (load up (force 0 0 10N) (region (- (abs z) 1))) (res 2))
(require stiff (compliance s < 1e-9N*mm))
(require below (compliance mixed < 0))";
    for (name, c) in checks(src) {
        assert_eq!(c.status, ink::claims::Status::Unknown, "{name}: {c:?}");
        let ink::claims::Evidence::Metric { measured, .. } = &c.evidence else { panic!("{name}: {c:?}") };
        assert!(measured.value.as_ref().is_err_and(|m| m.contains("reaction")), "{name}: {c:?}");
    }
}

#[test]
fn a_move_does_no_work_in_compliance() {
    // Compliance is the loads' work: a part moved whole by a fix, free of
    // stress, read 1.9e8 N mm, its penalty's; a 10 mm cube pressed 0.01 mm
    // (111 N, E A d / L) read 19,999.
    let s = solve("(study s (shape (box 4 4 4)) (elastic (E 1000MPa) (nu 0.3)) (fix all -1 (move 1 1 1)))", "s", 1.0);
    assert!(s.report.energy.abs() < 1e-9, "compliance {}", s.report.energy);
    let s = solve("(study s (shape (box 10 10 10)) (elastic (E 1000MPa) (nu 0)) (fix base (+ z 4.5) (dofs z)) (fix press (- 4.5 z) (dofs z) (move 0 0 -0.01)))", "s", 0.5);
    assert!(s.report.energy.abs() < 1e-9 && (s.report.reactions[1].1[2] / -111.1 - 1.0).abs() < 0.01, "compliance {}, {:?}", s.report.energy, s.report.reactions);
}

#[test]
fn moves_are_refused_where_they_mean_nothing() {
    // Each refused at load, naming why.
    let bar = |physics: &str, fix: &str, extra: &str| {
        format!("(import \"std\")\n(study s (shape (box 20 4 4)) {physics}\n  (fix wall (+ x 9.5)) {fix}\n  (load tip (force 0 0 -1N) (region (- 9 x))) {extra})")
    };
    let elastic = "(elastic (E 1000MPa) (nu 0.3) (density 1000kg/m3))";
    let heat = "(heat (k 20W/m/K))";
    let moved = "(fix end (- 9.5 x) (dofs z) (move 0 0 -0.1))";
    for (physics, fix, extra, msg) in [
        (elastic, moved, "(modes 2)", "a modal study takes no loads"),
        (elastic, "", "(transient (duration 0.01s) (steps 2))", "applies to heat studies"),
        (heat, moved, "", "apply to elastic studies"),
    ] {
        let src =
            bar(physics, fix, extra).replace(if physics == heat { "(fix wall (+ x 9.5))" } else { "" }, if physics == heat { "(temperature wall (+ x 9.5) 20degC)" } else { "" });
        let src = if physics == heat { src.replace("(load tip (force 0 0 -1N) (region (- 9 x)))", "(source tip (power 1W) (region (- 9 x)))") } else { src };
        let e = ink::load_str(&src, "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{fix} {extra} loaded"));
        assert!(e.message.contains(msg), "{fix} {extra}: {e:?}");
    }
}

/// A study's field of its mode k: `(study-field s q (mode k))` is that
/// mode's shape (1 mm peak) or its stresses, of any mode it solved, its
/// source naming the mode and its frequency; a mode it did not solve is an
/// error saying which it has.
#[test]
fn a_study_field_reads_any_mode() {
    let field = |d: &ink::lang::Doc, name: &str, p: [f64; 3]| {
        let ink::lang::Value::Field(f, _) = d.def(name).unwrap().clone() else { panic!() };
        let t = ink::tape::Tape::new(&d.g, f);
        let what = t.grids.first().and_then(|g| g.source()).map(|s| s.what.clone()).unwrap_or_default();
        (ink::eval::eval_point(&t, &d.ctx(&[], 0.0).unwrap(), p), what)
    };
    let vib = "(body beam (shape (box 100 10 10)) (attr density 7.85g/cm3))
(study vib (shape (shape-of beam)) (elastic (E 200000MPa) (nu 0.3) (density (by-body beam density)))
  (fix wall (+ x 49.5)) (modes 2) (res 1.5))
(def m1 (study-field vib displacement-y (mode 1)))
(def m2 (study-field vib displacement-y (mode 2)))
(def m2z (study-field vib displacement-z (mode 2)))
(def s2 (study-field vib von-mises (mode 2)))";
    let mut d = ink::load_str(&format!("(import \"std\")\n{vib}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let (s, _, _) = sim::run_study(&mut d, &ctx, "vib", None, 1e-8).unwrap();
    let tip = [48.0, 0.0, 0.0];
    for (k, def) in [(1, "m1"), (2, "m2")] {
        let want = s.mode(k).unwrap().value_at(tip).unwrap()[1];
        let (got, what) = field(&d, def, tip);
        let hz = s.report.modes[k - 1].0;
        assert!((got - want).abs() < 1e-3, "mode {k}'s tip, y: {got} read, {want} the solution's");
        assert!(
            what.starts_with(&format!("study vib: displacement-y, mode {k}: "))
                && what.ends_with(" Hz, 1 mm peak")
                && (what.split(' ').nth(5).unwrap().parse::<f64>().unwrap() / hz - 1.0).abs() < 1e-3,
            "{what} ({hz} Hz)"
        );
    }
    let (y, z) = (field(&d, "m2", tip).0, field(&d, "m2z", tip).0);
    assert!((y.hypot(z) - 1.0).abs() < 0.1, "mode 2's tip moves 1 mm: {y}, {z}");
    assert!(field(&d, "s2", [-45.0, 0.0, 4.5]).0 > 0.0, "mode 2's stress at the wall");
    let e =
        format!("{:?}", ink::load_str(&format!("(import \"std\")\n{vib}\n(def m3 (study-field vib displacement-y (mode 3)))"), "test", Path::new(ROOT)).err().expect("no mode 3"));
    assert!(e.contains("has no mode 3") && e.contains("gives 2 vibration modes"), "{e}");
}

#[test]
fn by_body_takes_each_bodys_own_property_or_the_default() {
    // A body without the property takes the default (or loading fails),
    // never an earlier body's value, whatever the order.
    let load = |attrs: [&str; 3], default: &str| {
        let bodies: String = attrs.iter().enumerate().map(|(i, a)| format!("(body s{i} (shape (move (sphere 1) {} 0 0)) {a})\n", 3 * i)).collect();
        let src = format!("(import \"std\")\n{bodies}(assembly root (instance a s0) (instance b s1) (instance c s2))\n(def e (by-body root E {default}))");
        ink::load_str(&src, "t", Path::new(ROOT)).map(|d| {
            let ink::lang::Value::Field(f, _) = d.def("e").unwrap().clone() else { panic!() };
            let t = ink::tape::Tape::new(&d.g, f);
            [0.0, 3.0, 6.0].map(|x| ink::eval::eval_point(&t, &ink::eval::Ctx::new(0.0, vec![]), [x, 0.0, 0.0]))
        })
    };
    assert_eq!(load(["(attr E 10)", "(attr E 20)", ""], "99").unwrap(), [10.0, 20.0, 99.0]);
    assert_eq!(load(["", "(attr E 10)", "(attr E 20)"], "99").unwrap(), [99.0, 10.0, 20.0]);
    assert!(load(["(attr E 10)", "", "(attr E 20)"], "").is_err() && load(["", "(attr E 10)", "(attr E 20)"], "").is_err());
}

#[test]
fn study_keys_follow_content_not_names() {
    // The same study written two ways shares a key; changing a param the
    // study uses changes it, an unrelated one does not.
    let a = "(param w 10) (param other 1)
(def blk (box 40 w 10)) (body b (shape blk))
(study s (shape (shape-of b)) (elastic (E 2000MPa) (nu 0.3)) (fix f (- x -19)) (load l (force 0 0 -1N) (region (- 19 x))))";
    let b = "(param w 10) (param other 1)
(body b (shape (box 40 w 10)))
(study renamed (shape (shape-of b)) (elastic (E 2000MPa) (nu 0.3)) (fix f (- x -19)) (load l (force 0 0 -1N) (region (- 19 x))))";
    let key = |src: &str, name: &str, ov: &[(String, ink::lang::Given)]| {
        let mut d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
        let ctx = d.ctx(ov, 0.0).unwrap();
        sim::study_setup(&mut d, &ctx, name, None, 1e-9).unwrap().0.key(&ctx).unwrap().0
    };
    let k = key(a, "s", &[]);
    assert_eq!(k, key(b, "renamed", &[]));
    assert_ne!(k, key(a, "s", &[("w".into(), 11.0.into())]));
    assert_eq!(k, key(a, "s", &[("other".into(), 5.0.into())]));
}

#[test]
fn cached_studies_keep_which_load_is_which() {
    // Two loads, the pressure moved from one to the
    // other, and two supports with their names swapped. Each pair solved
    // in one document (the second a memory hit, or a miss) and the second
    // again alone (a disk hit, or a miss): the pressure pushes where it is
    // now, and each reaction goes to the support now so named.
    let study = |name: &str, upper: &str, lower: &str| {
        format!(
            "(study {name} (shape (shape-of bar)) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (+ x 3)) (load upper {upper} (region (- 0.5 z))) (load lower {lower} (region (+ z 0.5))) (res 0.5))"
        )
    };
    let fixes = |name: &str, a: &str, b: &str| {
        format!(
            "(study {name} (shape (shape-of bar)) (elastic (E 2000MPa) (nu 0.3))
  (fix {a} (+ x 3)) (fix {b} (- 3 x)) (load push (force 0 0 -1N) (region (+ x 2))) (res 0.5))"
        )
    };
    let bar = "(body bar (shape (box 8 2 2)))";
    let lower = study("lower", "(force 0 0 0)", "(pressure 1MPa)");
    let after = fixes("after", "right", "left");
    let both = format!("{bar}\n{}\n{lower}\n{}\n{after}", study("upper", "(pressure 1MPa)", "(force 0 0 0)"), fixes("before", "left", "right"));
    let run = |d: &mut ink::lang::Doc, name: &str| {
        let ctx = d.ctx(&[], 0.0).unwrap();
        sim::run_study(d, &ctx, name, None, 1e-9).unwrap().0
    };
    let load = |src: &str| ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
    let check = |d: &mut ink::lang::Doc| {
        // Pressure on the bottom face pushes the free end up, and its
        // stresses are this load's: the top near the wall in compression
        // (12.6 MPa; the one on top stretches it).
        let s = run(d, "lower");
        assert!(s.value_at([3.0, 0.0, 0.0]).unwrap()[2] > 0.1 && s.field_at_point([-1.0, 0.0, 0.75]).unwrap().0[0] < -10.0, "the pressure acts on 'upper'");
        // The push is near the end fixed by 'right' now.
        let s = run(d, "after");
        let r = |n: &str| s.report.reactions.iter().find(|r| r.0 == n).unwrap().1[2];
        assert!(r("right") > 0.9 && r("left") < 0.1, "{:?}", s.report.reactions);
    };
    let mut d = load(&both);
    run(&mut d, "upper");
    run(&mut d, "before");
    check(&mut d);
    let mut alone = load(&format!("{bar}\n{lower}\n{after}"));
    check(&mut alone);
}

#[test]
fn thin_walls_bend_accurately() {
    // A 4 mm plate strip as a cantilever, only 2 cells through its
    // thickness: incompatible-mode bricks stay within 5% of beam theory
    // (plain trilinear bricks were 14% too stiff here).
    let src = "(body beam (shape (move (box 100 10 4) 50 0 0)))
(study bend (shape (shape-of beam)) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (- x 0.5)) (load tip (force 0 0 -1N) (region (- 99.5 x))))";
    let s = solve(src, "bend", 2.0);
    let (l, i) = (99.25f64, 10.0 * 64.0 / 12.0);
    let theory = l.powi(3) / (3.0 * 2000.0 * i);
    let got = -s.value_at([99.75, 0.0, 0.0]).unwrap()[2];
    assert!((got - theory).abs() < 0.05 * theory, "{got} vs {theory}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn internal_pressure_gives_lame_hoop_stress() {
    // A tube (radii 19 and 21) under 1 MPa inside, ends held axially: the
    // hoop stress at mid-wall is p ri^2 (1 + ro^2 / r^2) / (ro^2 - ri^2) = 9.49 MPa.
    let src = "(def tube (max (- (abs (- (len2 x y) 20)) 1) (- (abs z) 5)))
(body vessel (shape tube))
(study press (shape (shape-of vessel)) (elastic (E 200000MPa) (nu 0.3))
  (fix ends (- 3 (abs z)) (dofs z))
  (fix nospin (max (- (abs x) 1) (- 18 y) (- (abs z) 1)) (dofs x))
  (fix noslide (max (- (abs y) 1) (- 18 x) (- (abs z) 1)) (dofs y))
  (load p (pressure 1MPa) (region (- (len2 x y) 20))))";
    let s = solve(src, "press", 0.5);
    let (sg, _, _) = s.field_at([0.0, -20.0, 0.0]).unwrap();
    assert!((sg[0] - 9.49).abs() < 0.05 * 9.49, "hoop {}", sg[0]);
    // The tube grows evenly: the load must not ripple with the grid (a
    // cos 2 theta ripple of 1e-3 visibly ovalizes thin rings).
    let ur: Vec<f64> = [0.0f64, 15.0, 30.0, 45.0, 60.0, 75.0]
        .map(|a| {
            let (c, sn) = (a.to_radians().cos(), a.to_radians().sin());
            let u = s.value_at([20.0 * c, 20.0 * sn, 0.0]).unwrap();
            u[0] * c + u[1] * sn
        })
        .to_vec();
    let (lo, hi) = ur.iter().fold((f64::INFINITY, 0.0f64), |(a, b), &u| (a.min(u), b.max(u)));
    assert!(hi - lo < 0.02 * hi, "radial growth by angle {ur:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn poissons_ratio_is_a_field() {
    // A bar in tension whose lower half has nu 0.1 and upper half nu 0.45:
    // far from the joint each half contracts sideways exactly as a bar made
    // wholly of its material (and by about nu sigma / E).
    let study = |nu: &str| {
        format!(
            "(def bar (max (- x 5) (neg x) (- y 5) (neg y) (- z 100) (neg z)))
(body b (shape bar))
(study pull (shape (shape-of b)) (elastic (E 1000MPa) (nu {nu}))
  (fix sx (- x 0.25) (dofs x))
  (fix sy (- y 0.25) (dofs y))
  (fix sz (- z 0.25) (dofs z))
  (load end (force 0 0 250N) (region (- 99.5 z))))"
        )
    };
    let contraction = |s: &sim::Solution, z: f64| s.value_at([5.0, 2.5, z]).unwrap()[0] - s.value_at([2.5, 2.5, z]).unwrap()[0];
    let mixed = solve(&study("(select (- z 50) 0.1 0.45)"), "pull", 1.0);
    for (z, nu) in [(25.0, 0.1), (75.0, 0.45)] {
        let alone = solve(&study(&nu.to_string()), "pull", 1.0);
        let (got, want) = (contraction(&mixed, z), contraction(&alone, z));
        assert!((got - want).abs() < 0.005 * want.abs(), "z {z}: {got} vs {want} for nu {nu} alone");
        let hand = -nu * 10.0 / 1000.0 * 2.5;
        assert!((got - hand).abs() < 0.06 * hand.abs(), "z {z}: {got} vs {hand} by hand");
    }
}

#[test]
fn stress_probes_read_the_point_not_the_cell() {
    // A cantilever 4 mm deep, 4 cells through the depth: bending stress is
    // linear through the depth, sigma = M z / I. Probes near the top and
    // bottom faces must follow it, not repeat their cells' centre values.
    let src = "(def beam (max (- x 60) (neg x) (- (abs y) 2) (- (abs z) 2)))
(body b (shape beam))
(study bend (shape (shape-of b)) (elastic (E 200000MPa) (nu 0.3))
  (fix wall (- x 0.5))
  (load tip (force 0 0 -10N) (region (- 59.5 x))))";
    let s = solve(src, "bend", 1.0);
    // At x = 30: M = 10 * (60 - 30 - 0.25) (load centroid 0.25 from the end), I = 4^4 / 12.
    let m = 10.0 * 29.75;
    let i = 256.0 / 12.0;
    for z in [1.8f64, 1.2, -1.2, -1.8] {
        let (sg, _, _) = s.field_at_point([30.0, 0.3, z]).unwrap();
        let want = m * z / i;
        assert!((sg[0] - want).abs() < 0.06 * m * 2.0 / i, "z {z}: {} vs {want}", sg[0]);
    }
    // The cell-centre value cannot tell 1.2 from 1.8 (same cell at res 1).
    let (a, _, _) = s.field_at([30.0, 0.3, 1.2]).unwrap();
    let (b, _, _) = s.field_at([30.0, 0.3, 1.8]).unwrap();
    assert_eq!(a[0], b[0]);
}

#[test]
fn torque_loads_twist_a_shaft() {
    // A steel rod (r 5) held at its base, 10 N m at its far end: the twist
    // rate is T / (G J) and the shear stress T r / J.
    let src = "(def rod (max (- (len2 x y) 5) (neg z) (- z 60)))
(body b (shape rod))
(study twist (shape (shape-of b)) (elastic (E 200000MPa) (nu 0.3))
  (fix base (- z 1))
  (load drive (torque 0 0 10000N*mm) (region (- 55 z))))";
    let s = solve(src, "twist", 0.5);
    let j = std::f64::consts::PI * 625.0 / 2.0;
    let g = 200000.0 / 2.6;
    let rate = 10000.0 / (g * j);
    let du = s.value_at([5.0, 0.0, 40.0]).unwrap()[1] - s.value_at([5.0, 0.0, 20.0]).unwrap()[1];
    assert!((du - rate * 20.0 * 5.0).abs() < 0.01 * rate * 100.0, "twist {du} vs {}", rate * 100.0);
    let (sg, _, _) = s.field_at_point([3.5, 0.0, 30.0]).unwrap();
    assert!((sg[4] - 10000.0 * 3.5 / j).abs() < 0.02 * 10000.0 * 3.5 / j, "shear {}", sg[4]);
    // No net force: the base reacts only the torque.
    let r: f64 = s.report.reactions.iter().map(|(_, v)| v[0].abs() + v[1].abs() + v[2].abs()).sum();
    assert!(r < 1e-6 * 10000.0, "net reaction {r}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn normal_only_supports_slide_freely() {
    // A quarter block tilted 30 degrees off the grid, held only along the
    // normals of its bottom and symmetry faces, pressed along its axis: the
    // stress is uniaxial, -2.5 MPa, everywhere (free Poisson expansion).
    // Each region takes in its face and not the rims of the faces beside
    // it, which a frictionless hold would hold along their own normals.
    let src = "(def blk (move (box 10 10 20) 5 5 10))
(body b (shape (rotate-x blk 30)))
(study press (shape (shape-of b)) (elastic (E 1000MPa) (nu 0.3))
  (fix bottom (rotate-x (max (- z 0.3) (- (abs (- x 5)) 4.99) (- (abs (- y 5)) 4.99)) 30) (dofs normal))
  (fix symx (rotate-x (max (- x 0.3) (- (abs (- y 5)) 4.99) (- (abs (- z 10)) 9.99)) 30) (dofs normal))
  (fix symy (rotate-x (max (- y 0.3) (- (abs (- x 5)) 4.99) (- (abs (- z 10)) 9.99)) 30) (dofs normal))
  (load top (force 0 125N -216.506N) (region (rotate-x (- 19.7 z) 30))))";
    let s = solve(src, "press", 0.5);
    let (c, sn) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
    let w = [0.0, -sn, c];
    for p in [[5.0, -2.5, 12.0], [8.0, -4.0, 15.0], [2.0, 1.0, 6.0]] {
        let (sg, vm, _) = s.field_at_point(p).unwrap();
        let m = [[sg[0], sg[3], sg[5]], [sg[3], sg[1], sg[4]], [sg[5], sg[4], sg[2]]];
        let sww: f64 = (0..3).map(|i| (0..3).map(|j| w[i] * m[i][j] * w[j]).sum::<f64>()).sum();
        assert!((sww + 2.5).abs() < 0.05, "sigma along the axis at {p:?}: {sww}");
        assert!((vm - 2.5).abs() < 0.08, "uniaxial at {p:?}: von Mises {vm}");
    }
}

/// A 210 x 10 x 10 beam on its bottom edge at x = 0 and 200: a pin (held,
/// free to turn about the edge line) and a roller (z held, turning about it).
const SIMPLY_SUPPORTED: &str = "(body b (shape (move (box 210 10 10) 100 0 0)))
(study ss (shape (shape-of b)) (elastic (E 70000MPa) (nu 0.3))
  (fix pin (- (len2 x (+ z 5)) 1) (pivot (at 0 0 -5) (axis 0 1 0)))
  (fix roller (- (len2 (- x 200) (+ z 5)) 1) (dofs z) (pivot (at 200 0 -5) (axis 0 1 0)))
  (load mid (force 0 0 -100N) (region (- (abs (- x 100)) 1))))";

/// A pivot's axis, a slide's direction and a face's normal must have a
/// length: each is refused at its form.
#[test]
fn directions_without_a_length_are_refused_in_supports() {
    let load = |src: String| ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).err().map(|e| e.message);
    let named = "(body bar (shape (box 100 10 10)))\n(study pull (shape (shape-of bar)) (elastic (E 1000MPa) (nu 0))\n  (fix wall (face bar (normal 0 0 0) (depth 0.5))) (load pull (force 1000N 0 0) (region (- x 49))))";
    for (src, msg) in [
        (SIMPLY_SUPPORTED.replacen("(axis 0 1 0)", "(axis 0 0 0)", 1), "a direction must be nonzero"),
        (SIMPLY_SUPPORTED.replacen("(dofs z)", "(dofs z) (slide 0 0 0)", 1), "a direction must be nonzero"),
        (named.to_string(), "normal must be non-zero"),
    ] {
        assert!(load(src.clone()).is_some_and(|m| m.contains(msg)), "{src}: {:?}", load(src.clone()));
    }
}

fn study_error(src: &str) -> String {
    let mut d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    sim::run_study(&mut d, &ctx, "s", None, 1e-9).err().expect("the study solved").message
}

#[test]
fn balanced_loads_leave_a_mechanism_at_rest() {
    // Pulled along its length, through the pin's axis, the lever is not
    // turned: the free turn is held at rest (and the study says so), and
    // the bar carries 100 N over 12 x 6.
    let src = "(body l (shape (subtract (move (box 110 12 6) 45 0 0) (cylinder-z 3 20))))
(study s (shape (shape-of l)) (elastic (E 70000MPa) (nu 0.3))
  (fix pin (- (len2 x y) 4.5) (pivot (at 0 0 0) (axis 0 0 1)))
  (load pull (force 100N 0 0) (region (- 95 x)))
  (res 1))";
    let s = solve(src, "s", 1.0);
    let w = &s.report.warnings;
    assert!(w.len() == 1 && ["(0, 0, 1)", "(0, 0, -1)"].iter().any(|a| w[0].contains(&format!("turn about the axis along {a} through (0, 0, 0)"))), "{w:?}");
    let (_, vm, _) = s.field_at_point([50.0, 0.0, 0.0]).unwrap();
    assert!((vm - 100.0 / 72.0).abs() < 0.02 * 100.0 / 72.0, "stress {vm}");
    let u = s.value_at([50.0, 0.0, 0.0]).unwrap();
    assert!(u[1].abs() < 0.01 * u[0], "no turn: {u:?}");
}

#[test]
fn claims_on_simulated_fields_are_numerical() {
    // A solved result is a numerical input: a claim on a shape built from it
    // holds numerically at best, however exactly the shape itself is proven.
    let src = "(def bar (move (box 20 4 4) 10 0 0))
(study s (shape bar) (elastic (E 2000MPa) (nu 0.3)) (fix wall (- x 0.5)) (load tip (force 0 0 -10N) (region (- 19.5 x))) (res 1))
(def stress (study-field s von-mises))
(def graded (max bar (- stress 1e9)))
(require inside-box (subset graded (box 100 100 100)))
(require plain (subset bar (box 100 100 100)))";
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let measured = |n: &str| d.asserts.iter().find(|a| a.name == n).is_some_and(|a| d.inexact(&d.claim_roots(&a.claim), &[]).any());
    assert!(measured("inside-box") && !measured("plain"));
}

#[test]
fn a_result_holds_its_boundary_value_outside_its_grid() {
    // A package's result is a quantity, not a distance: past its grid it
    // holds the boundary value. Grown by the distance (as an imported
    // mesh's signed distance is), a cantilever's tip displacement of -9.1
    // mm read +8.9 mm 20 mm past the grid, so a part warped by it,
    // exaggerated, pulled material out into empty space.
    let src = "(def bracket (max (box 80 20 6) (neg (move (cylinder-z 4 10) 10 0 0))))
(study bend (shape bracket) (elastic (E 2000MPa) (nu 0.35))
  (fix wall (+ x 39.5)) (load tip (force 0 0 -40N) (region (- 38 x))) (res 1))
(def dz (study-field bend displacement-z))";
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let ink::lang::Value::Field(f, _) = d.def("dz").unwrap().clone() else { panic!() };
    let t = ink::tape::Tape::new(&d.g, f);
    let at = |p: [f64; 3]| ink::eval::eval_point(&t, &ctx, p);
    let (_, hi) = t.grids[0].corners();
    let edge = at([hi[0], 0.0, 0.0]);
    assert!(edge < -5.0, "the tip bends down: {edge}");
    for x in [hi[0] + 3.0, hi[0] + 20.0, hi[0] + 200.0] {
        assert_eq!(at([x, 0.0, 0.0]), edge, "at x = {x}, {} past the grid", x - hi[0]);
    }
}

#[test]
fn materials_come_from_the_ir_library() {
    // The table is ../kernel/std/materials.ink (plain IR data), not kernel code: the
    // Canonical grades and compatibility names read the same plain data.
    let petg = ink::materials::lookup("petg").unwrap();
    assert_eq!((petg["E"].as_f64(), petg["nu"].as_f64(), petg["name"].as_str()), (Some(1500.0), Some(0.38), Some("petg")));
    assert_eq!(ink::materials::lookup("stainless-304").unwrap()["yield"].as_f64(), Some(210.0));
    // Brittle materials have no yield: one limit for every stress state
    // would be wrong (concrete's compressive 30 MPa gave a tension coupon at
    // 5 MPa a safety factor of 5.9; its tensile strength is 2.9).
    for m in ["glass", "cast-iron", "tungsten-carbide", "borosilicate-glass", "polystyrene", "tpu"] {
        assert!(ink::materials::lookup(m).unwrap().get("yield").is_none(), "{m}");
    }
    assert!(ink::materials::lookup("glass").unwrap().get("tensile").is_none());
    assert_eq!(ink::materials::lookup("glass").unwrap()["bending"].as_f64(), Some(45.0));
    // Another spelling is no material: it is a hint naming the one to write
    // (a generic name, its commonest grade).
    for (spelling, name) in [
        ("PETG", "petg"),
        ("al6061", "aluminum-6061"),
        ("stainless_304", "stainless-304"),
        ("aluminium", "aluminum-6061"),
        ("Stainless Steel", "stainless-304"),
        ("acetal", "pom"),
        ("PP", "polypropylene"),
    ] {
        assert!(ink::materials::lookup(spelling).is_none(), "{spelling}");
        assert_eq!(ink::materials::instead(spelling), Some(name), "{spelling}");
    }
    assert_eq!(ink::materials::instead("unobtainium"), None);
    assert_eq!(ink::materials::names().len(), 30);
    assert_eq!(ink::materials::fluids(), ["water-liquid-20c-1bar", "air-unspecified-reference", "oil-unspecified-reference"]);
    assert_eq!(ink::materials::lookup("water").unwrap()["viscosity"].as_f64(), Some(0.00100160));
}

#[test]
fn fix_options_are_checked() {
    // An option a fix does not have is an error at its location, never
    // silently dropped along with the valid (dofs ...) (which would then
    // hold every axis instead).
    let study = |fix: &str| {
        format!(
            "(import \"std\")\n(body bar (shape (move (box 20 1 4) 10 0 0)))\n(study s (shape (shape-of bar)) (elastic (material petg))\n  \
             (fix base (- x 0.5))\n  {fix}\n  (load pull (force 1N 0 0) (region (- 19.5 x))) (res 0.5))"
        )
    };
    for (fix, msg) in [
        ("(fix drive (- 19.5 x) (dofs y) (value 0.5))", "unknown fix option 'value'"),
        ("(fix drive (- 19.5 x) (dofs y) (dofs x))", "dofs given twice"),
        ("(fix drive (- 19.5 x) (dofs normal y))", "stands alone"),
        ("(fix drive (- 19.5 x) (dofs))", "at least one"),
        ("(fix drive (- 19.5 x) y)", "must be a list"),
        ("(fix drive (- 19.5 x) (pivot (axis 0 0 1)))", "pivot needs (at x y z)"),
        ("(fix drive (- 19.5 x) (pivot (at 0 0)))", "at takes three numbers"),
        ("(fix drive (- 19.5 x) (pivot (at 0 0 0) (spin 1)))", "pivot takes (at x y z)"),
        ("(fix drive (- 19.5 x) (slide 0 0 0))", "nonzero"),
        ("(fix drive (- 19.5 x) (move 1 0))", "three displacements"),
        ("(fix drive (- 19.5 x) (move 0 0 1) (move 0 0 1))", "move given twice"),
        ("(fix drive (- 19.5 x) (move 0 0 (* 1mm/s t)))", "not t"),
        ("(fix drive (- 19.5 x) (dofs y) (move 0.1 0 0))", "along x, which the fix's (dofs ...) leave free"),
        ("(fix drive (- 19.5 x) (push-only 0 0 1))", "unknown fix option 'push-only'"),
    ] {
        let e = ink::load_str(&study(fix), "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{fix} loaded"));
        assert!(e.message.contains(msg) && e.line == 5 && e.hint.as_deref().is_some_and(|h| h.contains("(dofs x y z)")), "{fix}: {e:?}");
    }
    // The valid option keeps its meaning: the drive end holds y only and
    // slides in x under the pull (F L / E A, with the library's current E).
    let mut d = ink::load_str(&study("(fix drive (- 19.5 x) (dofs y))"), "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let (s, _, _) = sim::run_study(&mut d, &ctx, "s", None, 1e-9).unwrap();
    let ux = s.value_at([19.75, 0.0, 0.0]).unwrap()[0];
    let e_modulus = ink::materials::lookup("petg").unwrap()["E"].as_f64().unwrap();
    assert!((ux / (19.5 / (e_modulus * 4.0)) - 1.0).abs() < 0.1, "ux {ux}");
    // A free motion along what the fix leaves free already is a mistake.
    let mut d = ink::load_str(&study("(fix drive (- 19.5 x) (dofs y) (slide 1 0 0))"), "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let e = sim::run_study(&mut d, &ctx, "s", None, 1e-9).err().unwrap();
    assert!(e.message.contains("moves nothing the fix holds"), "{e:?}");
}

/// A model's view at `overrides`, with its declared tolerances applied (as
/// `ink require`).
fn toleranced(src: &str, overrides: &[(String, ink::lang::Given)]) -> Result<ink::view::ModelView, ink::InkError> {
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(overrides, 0.0).unwrap();
    let mut v = ink::view::ModelView::new(d, ctx);
    let spans = v.declared_spans()?;
    v.set_param_spans(&spans)?;
    Ok(v)
}

/// Each claim of a model, checked as `ink require` checks it.
#[test]
fn material_the_samples_miss_beside_a_part_is_not_settled_and_a_plain_part_is() {
    // A 0.001 mm foil (and a 0.005 mm web) from post to post above a bar,
    // on 1 mm cells: between the samples, the mesh still one piece, so no
    // separation shows it. The lines between the samples outside the part
    // meet it, joining the material the cells hold in two places apart:
    // not settled, naming it. The bar and posts alone: decided.
    let part = |extra: &str| {
        format!(
            "(def bar (box 10 4 4))
(def posts (min (move (box 1 1 1.4) -4.5 0 2.6) (move (box 1 1 1.4) 4.5 0 2.6)))
(def body (min bar posts {extra}))
(study s (shape body) (elastic (E 10000MPa) (nu 0) (density 1000kg/m3)) (fix wall (+ x 4.5)) (load pull (force 10N 0 0) (region (- x 4.5))) (res 1))
(require stiff (max-displacement s < 1))"
        )
    };
    for web in ["(move (box 10 1 0.001) 0 0 3.2)", "(move (box 10 4 0.005) 0 0 3.2)"] {
        let c = checks(&part(web)).remove(0).1;
        assert!(c.status == ink::claims::Status::Unknown && format!("{c:?}").contains("samples missed"), "{web}: {c:?}");
    }
    let c = checks(&part("")).remove(0).1;
    assert_eq!(c.status, ink::claims::Status::Holds, "{c:?}");
}

fn checks(src: &str) -> Vec<(String, ink::claims::Check)> {
    use ink::claims::{ClaimOptions, run_asserts};
    let mut v = toleranced(src, &[]).unwrap();
    run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts.into_iter().map(|a| (a.name, a.check)).collect()
}

/// A message's words: the names and numbers a test wrote, which it may
/// look for there without the kernel's words around them.
fn words(s: &str) -> Vec<&str> {
    s.split(|c: char| !(c.is_alphanumeric() || "-_./".contains(c))).filter(|w| !w.is_empty()).collect()
}

/// A metric claim's grid error, where three grids read it.
fn grid_error(c: &ink::claims::Check) -> Option<f64> {
    match &c.evidence {
        ink::claims::Evidence::Metric { grid_error, .. } => grid_error.map(|g| g.0),
        e => panic!("not a metric: {e:?}"),
    }
}

/// A package metric claim's measured value, or why it has none.
fn measured(c: &ink::claims::Check) -> Result<f64, String> {
    match &c.evidence {
        ink::claims::Evidence::Metric { measured, .. } => measured.value.clone(),
        e => panic!("not a metric: {e:?}"),
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_nearly_incompressible_material_keeps_its_poissons_ratio() {
    // nu = 0.4999 is admissible and solved with its own nu (not rounded
    // onto the singular 0.5): the 4 x 4 section carrying 100 N has 6.25 MPa
    // mean axial stress, as with nu = 0.49.
    for nu in ["0.4999", "0.49"] {
        let src = format!(
            "(study pull (shape (box 10 4 4)) (elastic (E 2000MPa) (nu {nu}) (yield 40MPa))
  (fix left (+ x 4)) (load right (force 100N 0 0) (region (- 4 x))) (res 1))
(require stress-safe (max-stress pull < 1MPa))
(require deflection-safe (max-displacement pull < 0.00001))"
        );
        for (name, c) in checks(&src) {
            assert!(c.status == ink::claims::Status::Fails && c.accuracy.is_none(), "nu {nu}, {name}: {c:?}");
            let v = measured(&c).unwrap();
            if name == "stress-safe" {
                assert!(v > 6.0 && v < 8.0, "nu {nu}: stress {v}");
            }
        }
    }
}

#[test]
fn stress_with_no_cell_beyond_the_support_margin_is_no_measurement() {
    // A 2 mm cube held over x < -0.5 and pulled with 40 N over x > 0.5: 10
    // MPa mean axial stress. The default 2 mm support margin covers every
    // cell, so max-stress has nothing to measure: the claim is unknown,
    // saying why. With a margin of 0, the stress is measured (von Mises 7
    // MPa between the supports, sigma_xx 10), and where the hold ends on the
    // free surface it is singular: 6.03, 9.90 and 10.06 on cells of 0.5,
    // 0.25 and 0.125 mm, 11.09 on 0.0625, growing. No claim on it holds, and
    // none is settled on grids that converge only by luck (it failed on its
    // three grids' last difference, 3 x 0.16, where they moved by 4).
    let src = |margin: &str| {
        format!(
            "(study pull (shape (box 2 2 2)) (elastic (E 2000MPa) (nu 0.3) (yield 40MPa))
  (fix left (+ x 0.5)) (load right (force 40N 0 0) (region (- 0.5 x))) (res 0.5) {margin})
(require stress-safe (max-stress pull < 1MPa))
(require yield-safe (min-safety-factor pull > 10))"
        )
    };
    for (name, c) in checks(&src("")) {
        assert_eq!(c.status, ink::claims::Status::Unknown, "{name}: {c:?}");
        assert!(measured(&c).is_err_and(|m| m.contains("support margin")), "{name}: {c:?}");
    }
    for (name, c) in checks(&src("(support-margin 0)")) {
        assert!(c.status != ink::claims::Status::Holds && c.accuracy.is_none(), "{name}: {c:?}");
        let v = measured(&c).unwrap();
        assert!(if name == "stress-safe" { v > 5.0 && v < 12.0 } else { v > 3.3 && v < 8.0 }, "{name}: {v}");
    }
}

/// A steel strip 100 x 10 x 1 clamped at x < 0, `force` (N) at its tip,
/// and the claims written after it.
fn strip(force: &str, claims: &str) -> String {
    format!(
        "(body strip (shape (move (box 105 10 1) 47.5 0 0.5)))
(study tip (shape (shape-of strip)) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa))
  (fix wall (- x 0)) (load end (force {force}) (region (- 99 x))) (res 1))
{claims}"
    )
}

#[test]
fn the_peak_stress_and_its_claims_are_read_at_the_surface() {
    // 1 N across a 1 mm strip, one cell through it: the bending stress
    // 6 F (99.5 - x) / (b t^2), 58 MPa by the root, peaks at the faces,
    // where the probes read it; the cells' centres lie on the neutral
    // plane, where it is zero. Read there, the peak was 2.8 MPa, and a
    // safety factor of 90 held where the truth is 4.3. The claims are asked
    // on cells of 0.5 mm: on 1, 2 and 4 the safety factor read 4.23, 4.35
    // and 7.39, converging faster than the elements can, and its band of
    // three times that movement reached past 10; on 0.5, 1 and 2, 4.285
    // +- 2.49, a clear no.
    let src = strip("0 0 -1N", "(require safe (min-safety-factor tip > 10))\n(require low (max-stress tip < 20MPa))");
    let s = solve(&src, "tip", 1.0);
    let p = s.max_cell_scalar().unwrap();
    let beam = 6.0 * (99.5 - p.at[0]) / 10.0;
    assert!((p.value - beam).abs() < 0.03 * beam, "peak {} at {:?}, beam theory {beam}", p.value, p.at);
    let (_, probe, _) = s.field_at_point(p.at).unwrap();
    assert!((probe - p.value).abs() < 1e-9 * p.value, "the probe there reads {probe}, the peak {}", p.value);
    for (name, c) in checks(&src.replace("(res 1)", "(res 0.5)")) {
        assert_eq!(c.status, ink::claims::Status::Fails, "{name}: {:?}", measured(&c));
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_hole_s_concentration_is_read_at_its_edge() {
    // A plate 120 x 60 x 2 with a 6 mm hole under 10 MPa: Kirsch's Kt 3
    // (3.02 at d/W 0.1, Heywood) puts 30.2 MPa at the hole's edge. Read
    // at cell centres the peak was 20.4 MPa at 1 mm cells (-33%); at the
    // surface it is the probe's, 27.5 (-9%, the cells' own error).
    let src = "(body plate (shape (subtract (move (box 120 60 2) 0 0 1) (cylinder-z 3 10))))
(study pull (shape (shape-of plate)) (elastic (E 200000MPa) (nu 0.3))
  (fix left (+ x 58) (dofs x))
  (fix ly (max (+ x 58) (- (abs y) 1)) (dofs y))
  (fix lz (+ x 58) (dofs z))
  (load right (force 1200N 0 0) (region (- 58 x))))";
    let s = solve(src, "pull", 1.0);
    let p = s.max_cell_scalar().unwrap();
    assert!((p.value - 30.23).abs() < 0.12 * 30.23, "peak {} at {:?}", p.value, p.at);
    assert!((v3_norm2(p.at) - 3.0).abs() < 0.1, "at the hole's edge: {:?}", p.at);
    let (_, probe, _) = s.field_at_point([0.0, 3.0, 1.0]).unwrap();
    assert!((probe - p.value).abs() < 0.05 * probe, "probe at the edge {probe}, peak {}", p.value);
}

fn v3_norm2(p: [f64; 3]) -> f64 {
    p[0].hypot(p[1])
}

/// A metric claim's evidence: what it measured, and why that does not
/// decide it (with how to), if it does not.
fn metric_evidence(c: &ink::claims::Check) -> (&ink::lang::Measured, Option<&(String, String)>) {
    let ink::claims::Evidence::Metric { measured, undecided, .. } = &c.evidence else { panic!("{c:?}") };
    (measured, undecided.as_ref())
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_claim_that_counting_the_support_margin_would_change_is_unknown() {
    // 1000 N on x 1..2 beside a clamp at x < 0, a 10 x 10 section: 10 MPa
    // mean shear, 15 at the neutral axis, all within the 2 mm margin. The
    // peak beyond it is a few MPa (3.8, from 1.1 on cells of 2 mm: a
    // grid error of about 8), and the margin's (singular, 19) would fail
    // the claim: it is unknown, saying what it left out, where, and how
    // to decide it (it held, with a warning). Left out on purpose,
    // (excluding support-margin), it says so; and the peak beyond the
    // margin, at its edge, reads 0.47, 1.08 and 3.82 MPa on cells of 4, 2
    // and 1 mm, growing fourfold, far past what a stress's read is known to
    // (4%): it is the singular peak creeping past the margin, and the claim
    // is not settled (it held on 3 |x - x2|). The load beside
    // the clamp, and one inside it (taken whole by the reaction), are
    // warnings of the solve.
    let src = "(body beam (shape (move (box 105 10 10) 47.5 0 0)))
(study near (shape (shape-of beam)) (elastic (E 200000MPa) (nu 0.3))
  (fix wall (- x 0)) (load close (force 0 0 -1000N) (region (max (- 1 x) (- x 2)))) (res 1))
(study inside (shape (shape-of beam)) (elastic (E 200000MPa) (nu 0.3))
  (fix wall (- x 0)) (load lost (force 0 0 -1000N) (region (+ x 1))) (load tip (force 0 0 -1N) (region (- 99 x))) (res 1))
(require near-ok (max-stress near < 15MPa))
(require near-left-out (max-stress near < 15MPa (excluding support-margin)))";
    let c = checks(src);
    let (m, why) = metric_evidence(&c[0].1);
    let ex = m.excluded.iter().find(|e| e.key == "support-margin").expect("the margin's peak");
    assert!(c[0].1.status == ink::claims::Status::Unknown && ex.value > 15.0 && ex.what.contains("within 2 mm of fix 'wall'"), "{:?}", c[0].1);
    let (why, how) = why.expect("why it is undecided");
    assert!(why.contains("the verdict is fails") && how.contains("model the support's contact") && how.contains("(excluding support-margin)"), "{why} / {how}");
    let (m, why) = metric_evidence(&c[1].1);
    assert!(c[1].1.status == ink::claims::Status::Unknown && why.is_some_and(|w| w.0.contains("three grids moving apart")), "{:?}", c[1].1);
    assert!(m.warnings.iter().any(|w| w.contains("leaving out") && w.contains("on purpose")), "{:?}", m.warnings);
    let near = solve(src, "near", 1.0);
    assert!(near.report.warnings.iter().any(|w| w.contains("load 'close': 100% of it lies within the support margin")), "{:?}", near.report.warnings);
    let inside = solve(src, "inside", 1.0);
    let w = &inside.report.warnings;
    assert!(w.len() == 1 && w[0].contains("load 'lost': 100% of it lies inside a fix's region"), "{w:?}");
    let refused = ink::load_str(&format!("(import \"std\")\n{src}\n(require typo (max-stress near < 15MPa (excluding margin)))"), "test", Path::new(ROOT)).err().expect("refused");
    assert!(refused.message.contains("leaves out nothing called margin"), "{}", refused.message);
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_symmetry_plane_has_no_support_margin() {
    // A quarter of a 120 x 60 x 4 plate with a 6 mm hole, 10 MPa across:
    // Kirsch with Heywood's width correction gives 30.2 MPa at (0, 3), on
    // the x = 0 plane. Its x and y planes are holds of one axis normal to
    // the face they hold: the whole plate's field, with no singular edge,
    // so they have no margin, and the peak is read. With them set apart
    // as a clamp's edge, the peak read 14.9 and `max-stress < 25` held at
    // half the true stress. The pad holding z on part of a face is a
    // punch's edge, and keeps its margin.
    let src = "(body plate (shape (max (subtract (move (box 120 60 4) 0 0 2) (cylinder-z 3 20)) (neg x) (neg y))))
(study pull (shape (shape-of plate)) (elastic (E 200000MPa) (nu 0.3))
  (mirror symx x 0)
  (mirror symy y 0)
  (fix pad (max (- 50 x) (- y 5) (- z 0.01)) (dofs z) (support))
  (load right (force 1200N 0 0) (region (- 58 x)))
  (res 0.5))
(require peak (max-stress pull < 25MPa))";
    let s = solve(src, "pull", 0.5);
    assert_eq!(s.plane_holds, [true, true, false]);
    let p = s.max_cell_scalar().unwrap();
    assert!((p.value - 30.23).abs() < 0.06 * 30.23 && p.at[0] < 0.5 && (v3_norm2(p.at) - 3.0).abs() < 0.3, "peak {} at {:?}", p.value, p.at);
    // `< 25` is unknown: read 29.2, its coarser grids 25.8 and 23.8, its
    // grid error past its margin.
    let c = checks(src);
    assert!(c[0].1.status == ink::claims::Status::Unknown && grid_error(&c[0].1).is_some_and(|e| e > measured(&c[0].1).unwrap() - 25.0), "{:?}", c[0].1);
    // A hold of two axes keeps its margin.
    for dofs in ["x y", "x z"] {
        let s = solve(&src.replace("(mirror symx x 0)", &format!("(fix symx (- x 0.01) (dofs {dofs}))")), "pull", 1.0);
        assert!(!s.plane_holds[0], "(dofs {dofs})");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_claim_on_a_part_its_loads_buckle_is_unknown() {
    // The strip as a column pushed at 1000 N, 24 times its Euler load
    // (41.5 N): its safety factor against yield, 2.5, held. Its linear
    // buckling factor is 0.04: a static study says nothing past the
    // buckle, so the claim is unknown, saying so; left out on purpose,
    // (excluding buckling), it holds. Pulled, nothing buckles, and it
    // holds.
    let claims = "(require ok (min-safety-factor tip > 1))\n(require checked (min-safety-factor tip > 1 (excluding buckling)))";
    // (On cells of 0.5 its factor swings 2.454, 2.454, 2.456, the last change 4.6 times the one before: not
    // settled; on 0.25 its grids converge.)
    let c = checks(&strip("-1000N 0 0", claims).replace("(res 1)", "(res 0.25)"));
    let (m, why) = metric_evidence(&c[0].1);
    assert!(c[0].1.status == ink::claims::Status::Unknown && m.excluded.iter().any(|e| e.key == "buckling"), "{:?}", c[0].1);
    assert!(why.is_some_and(|w| w.0.contains("estimated linear buckling factor") && w.1.contains("stiffen the part")), "{why:?}");
    assert_eq!(c[1].1.status, ink::claims::Status::Holds, "{:?}", c[1].1);
    let pulled = checks(&strip("1000N 0 0", claims).replace("(res 1)", "(res 0.25)"));
    assert!(pulled.iter().all(|c| c.1.status == ink::claims::Status::Holds), "{pulled:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_thin_region_takes_the_whole_face_it_lies_under() {
    // 1 MPa on a 10 x 10 end face picked by a region d deep under it: 100
    // N at every d. It took 50 N at d 0.001, 54 to 63 at
    // 0.01 and 90 at 0.1 (1 mm cells): the region's share of a box round
    // each facet's centroid, half outside the part. A region beyond the
    // part is refused.
    let src = |d: f64, r: f64| {
        format!(
            "(body bar (shape (move (box 100 10 10) 50 0 0)))
(study s (shape (shape-of bar)) (elastic (E 200000MPa) (nu 0.3))
  (fix wall (- x 1)) (load p (pressure 1MPa) (region (- (- 100 {d}) x))) (res {r}))"
        )
    };
    for (d, r) in [(0.001, 1.0), (0.01, 1.0), (0.01, 0.3), (0.1, 1.0), (0.5, 1.0)] {
        let s = solve(&src(d, r), "s", r);
        let f = s.report.applied.iter().find(|a| a.0 == "p").unwrap().1;
        // The mesh half a cell fine rounds the face's edges: 99.87 N at 1 mm
        // cells however deep the region.
        assert!((f[0] + 100.0).abs() < 0.2 && f[1].abs() < 1e-3 && f[2].abs() < 1e-3, "d {d}, cells {r}: {f:?}");
    }
    assert!(study_error(&src(-0.001, 1.0)).contains("load 'p' (a pressure) touches none of the part's surface"));
}

#[test]
fn a_thin_convection_region_cools_the_whole_face() {
    // A slab 10 thick and 20 x 20, 1 W made evenly in it, k 0.2 W/(m K)
    // and h 50 W/(m^2 K) to 20 C on both large faces (and on a strip of
    // the sides as deep as the region): 1D, its faces at 20 + q / 2 h A =
    // 45 C and its centre 15.6 above them, 60.6. The regions d deep under
    // the faces read 85 C at d 0.001 and 74 to 82 at 0.01. A region
    // beyond the part is refused, as a fix's is.
    let src = |d: f64| {
        format!(
            "(body slab (shape (move (box 10 20 20) 5 0 0)))
(study s (shape (shape-of slab)) (heat (k 0.2W/m/K))
  (convection faces (- (- 5 {d}) (abs (- x 5))) (coefficient 50W/m2/K) (ambient 20degC))
  (source gen (power 1W) (region (- (abs (- x 5)) 10))))"
        )
    };
    for (d, r) in [(0.001, 1.0), (0.01, 1.0), (0.01, 0.5)] {
        let s = solve(&src(d), "s", r);
        let (centre, face) = (s.value_at([5.0, 0.0, 0.0]).unwrap()[0], s.value_at([0.0, 0.0, 0.0]).unwrap()[0]);
        assert!((centre - 60.625).abs() < 0.5 && (face - 45.0).abs() < 0.5, "d {d}, cells {r}: centre {centre}, face {face}");
    }
    assert!(study_error(&src(-0.001)).contains("convection 'faces' touches none of the part's surface"));
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_thin_wall_off_the_grid_reads_its_surface_stress() {
    // A 0.7 mm steel strip at z 0.15..0.85, 1 N at its tip: 6 F (99.5 -
    // x) / (b t^2) at its faces. On 0.25 mm cells its outer cells are
    // 0.4 full, and the peak was read from the full cell under them at
    // that cell's face, 0.1 mm in: 86 MPa against 118.
    let src = "(body strip (shape (move (box 105 10 0.7) 47.5 0 0.5)))
(study tip (shape (shape-of strip)) (elastic (E 200000MPa) (nu 0.3))
  (fix wall (- x 0)) (load end (force 0 0 -1N) (region (- 99 x))))";
    for r in [0.5, 0.25] {
        let p = solve(src, "tip", r).max_cell_scalar().unwrap();
        let beam = 6.0 * (99.5 - p.at[0]) / (10.0 * 0.49);
        assert!((p.value - beam).abs() < 0.05 * beam, "cells {r}: peak {} at {:?}, beam theory {beam}", p.value, p.at);
    }
}

#[test]
fn a_modulus_no_material_has_is_refused_with_its_unit() {
    // E in pascals (200e9) read as MPa made the part a million times
    // stiffer than steel, silently; GPa written with its unit is right.
    let src =
        |e: &str| format!("(body b (shape (box 10 10 10)))\n(study s (shape (shape-of b)) (elastic (E {e}) (nu 0.3)) (fix f (+ x 4)) (load l (force 1N 0 0) (region (- x 4))))");
    let refused = |e: &str| ink::load_str(&format!("(import \"std\")\n{}", src(e)), "test", Path::new(ROOT)).err().expect("refused").message;
    assert!(refused("200e9").contains("is a bare number"));
    let e = refused("200e9MPa");
    assert!(e.contains("stiffer than any material"), "{e}");
    assert!(refused("0").contains("not positive"));
    assert!(refused("1e-9MPa").contains("softer than any solid"));
    solve(&src("200GPa"), "s", 1.0);
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn every_static_solve_finds_its_buckling_load() {
    // The strip as a column, clamped at x < 0 and pushed at its far end:
    // Euler's fixed-free load, pi^2 E I / (4 L^2), is 41.5 N. Every static
    // elastic solve finds its buckling load (41.5 N within 5%): at 1000 N
    // the loads buckle it, and it says so; at 30 and 15 N they do not, and
    // it says nothing. Pulled, nothing compresses it: no factor below a
    // hundred; bent by 1 N, its compressed side buckles at 760 times it.
    // A safety factor claim at 30 N holds.
    let found = |force: &str| {
        let s = solve(&strip(force, ""), "tip", 0.5);
        (s.report.buckle.as_ref().and_then(|b| b.answer()), s.report.warnings)
    };
    for (force, euler) in [("-1000N 0 0", 41.5 / 1000.0), ("-30N 0 0", 41.5 / 30.0), ("-15N 0 0", 41.5 / 15.0)] {
        let (factor, warnings) = found(force);
        assert!(factor.is_some_and(|(_, f)| (f / euler - 1.0).abs() < 0.05), "{force}: {factor:?}");
        assert_eq!(warnings.iter().any(|w| w.starts_with("buckles")), euler < 1.0, "{force}: {warnings:?}");
    }
    // (Its low end is what was checked: that 100 times the loads leave it stable.)
    for force in ["1000N 0 0", "0 0 -1N"] {
        let (factor, warnings) = found(force);
        assert!(factor.is_some_and(|(low, _)| low == 100.0) && !warnings.iter().any(|w| w.contains("buckl")), "{force}: {factor:?} {warnings:?}");
    }
    // (On 0.25 mm cells: on 0.5 its factor swings 81.82, 81.81, 81.85, growing 4.6 times, and is not settled.)
    let c = checks(&strip("-30N 0 0", "(require ok (min-safety-factor tip > 1))").replace("(res 1)", "(res 0.25)"));
    assert_eq!(c[0].1.status, ink::claims::Status::Holds, "{:?}", c[0].1);
    // A safety factor counts buckling where its band reaches below the
    // claim's: at 15 N the column buckles at 2.77 times its load, long
    // before it yields, so a margin of 2 holds and one of 5 is unknown,
    // saying so; left out on purpose, the strength alone decides it.
    let claims = "(require two (min-safety-factor tip > 2))\n(require five (min-safety-factor tip > 5))\n(require five-strength (min-safety-factor tip > 5 (excluding buckling)))";
    // (On 0.25 mm cells, where its grids converge; on 0.5 they swing growing.)
    let c = checks(&strip("-15N 0 0", claims).replace("(res 1)", "(res 0.25)"));
    assert_eq!(c[0].1.status, ink::claims::Status::Holds, "{:?}", c[0].1);
    assert!(c[1].1.status == ink::claims::Status::Unknown && undecided(&c[1].1).contains("buckle"), "{:?}", c[1].1);
    assert_eq!(c[2].1.status, ink::claims::Status::Holds, "{:?}", c[2].1);
}

/// Why a package metric claim is undecided, or "".
fn undecided(c: &ink::claims::Check) -> String {
    match &c.evidence {
        ink::claims::Evidence::Metric { undecided, .. } => undecided.as_ref().map_or(String::new(), |u| u.0.clone()),
        e => panic!("not a metric: {e:?}"),
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_claim_holds_only_where_coarser_grids_agree() {
    // A d4 hole in a 100 x 50 x 5 plate pulled at 10 MPa: Peterson's
    // finite-width Kt gives a peak of 30.1 MPa. On 1 mm cells it reads
    // 22.2, and `< 25` held; on cells of 2 and 4 mm it read 16.3 and 12.7.
    // Cells of 2 mm do not resolve its radius (one of them to it), so its
    // grids do not check it, and no finer ones are taken for it (they cost
    // a plate with a hole up to 50 times the memory): unknown, naming the
    // cells at which its grid four times as coarse resolves it,
    // 2 / (4 x 1.45) mm. The plain bar's stress, 10 MPa on every grid,
    // still holds. (The wall's margin left out: its 22.5, widened by what
    // the read may miss under the hole's climb, would decide it first.)
    let src = "(body plate (shape (subtract (box 100 50 5) (cylinder-z 2 20))) (attr material \"steel\"))
(study def (shape (shape-of plate)) (elastic (material steel))
  (fix wall (+ x 49)) (load pull (force 2500N 0 0) (region (- 49 x))) (res 1))
(require peak (max-stress def < 25MPa (excluding support-margin)))
(require sf (min-safety-factor def > 9 (excluding support-margin)))";
    for (name, c) in checks(src) {
        assert_eq!(c.status, ink::claims::Status::Unknown, "{name}: {c:?}");
        assert!(undecided(&c).contains("cells of 2 mm cannot see the part") && undecided(&c).contains("refine (res r) to 0.34"), "{name}: {}", undecided(&c));
    }
    let bar = format!("{BAR}\n(require low (max-stress pull < 15MPa))");
    let (_, c) = &checks(&bar)[0];
    assert_eq!(c.status, ink::claims::Status::Holds, "{c:?}");
}

#[test]
fn a_part_too_thin_for_the_coarser_grids_is_checked_on_finer_ones() {
    // A 100 x 3 x 3 steel cantilever on 2 mm cells bends 2.4 mm under
    // 10 N, and `< 50` held; but on 8 mm cells, four times as coarse, the
    // tip load's region holds no material, and a claim checked on 2, 4
    // and 8 mm was unknown at the cells the file asks for. Its grids are
    // now the coarsest three that read the part, 1, 2 and 4 mm, each split
    // where the 2 mm grid splits the bar (thinner than two of its cells),
    // as deep: to 0.5, 1 and 2 mm. Read on them, 2.374, 2.375 and 2.356:
    // the split finest reads 0.8% under the same bar on whole 0.5 mm cells
    // (2.3748), a swing growing twentyfold. Nothing in the solve could grow
    // without bound (its load and hold span the bar's end), so it is given
    // three times the larger change, 0.057, which covers 2.3748.
    // (With the hold's cut cell gripped whole, the coarsest had read 2.301,
    // and the three were taken as converging. On 1 mm cells, nothing split,
    // it holds: 2.3748 +- 0.0028.)
    let src = "(body bar (shape (move (box 100 3 3) 50 0 0)) (attr material \"steel\"))
(study bend (shape (shape-of bar)) (elastic (material steel))
  (fix root (- x 1)) (load tip (force 0 0 -10N) (region (- 99 x))) (res 2))
(require stiff (max-displacement bend < 50))";
    let (_, c) = &checks(src)[0];
    assert_eq!(c.status, ink::claims::Status::Holds, "{c:?}");
    let ink::claims::Evidence::Metric { measured, estimate, .. } = &c.evidence else { panic!("{c:?}") };
    assert!(measured.res == Some(1.0) && measured.finest == Some(0.5), "{c:?}");
    assert!(estimate.as_ref().is_some_and(|e| e.contains("three grids close, not in order")), "{estimate:?}");
    let (v, lo, hi) = claim_band(c).unwrap();
    assert!(lo <= 2.3748 && 2.3748 <= hi, "{v} in [{lo}, {hi}]");
    let (_, c) = &checks(&src.replace("(res 2)", "(res 1)"))[0];
    assert_eq!(c.status, ink::claims::Status::Holds, "{c:?}");
    assert!(estimate.as_ref().is_some_and(|e| ["2", "1", "4"].iter().all(|h| words(e).contains(h))), "{estimate:?}");
}

/// A 0.4 mm saw slit through 9 of a 10 mm steel cantilever's depth, a 1 mm
/// ligament under it; 100 N at the tip.
const SLIT: &str = "(body bar (shape (subtract (move (box 100 10 10) 50 0 0) (move (box 0.4 20 9) 20.2 0 0.5))) (attr material \"steel\"))
(study def (shape (shape-of bar)) (elastic (material steel))
  (fix root (- x 1)) (load tip (force 0 0 -100N) (region (- 99 x))))
(require sf (min-safety-factor def > 2))
(require stiff (max-displacement def < 0.5))
(require loose (max-displacement def < 1))";

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_slit_whose_hinge_lies_inside_a_cell_opens() {
    // On cells on the slit's planes the tip bends 4.33 mm (1 mm cells),
    // 5.29 (0.5) and 5.69 (0.25). Off them, at ink's own 1.004 mm and at
    // 0.502, the cell over the ligament's top holds the block beside the
    // slit and a 0.016 mm sliver of ligament under it, and shares its
    // face's corners with the cell across the slit: one element's field
    // spans the hinge and holds the slit's sides together. It read 1.3 mm
    // (its nodes found tying the slit) and 2.6 (unwarned). The slit is air
    // thinner than two cells: the cells along it are split until it is
    // two of theirs across, and it opens.
    for h in [1.004, 0.502] {
        let s = solve(SLIT, "def", h);
        let tip = -s.value_at([99.5, 0.0, 0.0]).unwrap()[2];
        assert!(s.report.thin_left.is_none() && s.report.cells_refined > 0, "{h} mm: {:?}", s.report.warnings);
        assert!(tip > 0.95 * 4.33 && tip < 1.05 * 5.69, "{h} mm: tip {tip}");
    }
    // Its displacement claims decide: the bound at 1 mm, or under 1 N 0.01
    // mm (elastic, its aligned answer 0.043 mm), fail. Its peak stress lies
    // at the slit's sharp inside corner, singular, which no cells read: its
    // safety factor is unknown, saying so.
    for (name, c) in checks(SLIT) {
        if name == "sf" {
            assert!(c.status == ink::claims::Status::Unknown && undecided(&c).contains("sharp inside corner"), "{name}: {c:?}");
        } else {
            assert_eq!(c.status, ink::claims::Status::Fails, "{name}: {c:?}");
        }
    }
    let elastic = SLIT.replace("-100", "-1").replace("(require sf (min-safety-factor def > 2))\n", "").replace("< 0.5", "< 0.01");
    let (name, c) = &checks(&elastic)[0];
    assert!(name == "stiff" && c.status == ink::claims::Status::Fails, "{name}: {c:?}");
}

#[test]
fn a_thin_ligament_keeps_its_bending_in_a_cut_cell() {
    // A 20 x 4 x 4 steel bar slit to a 0.05 mm ligament (z -2 to -1.95)
    // under a tip load small enough to stay linear: cells of 0.1 mm bend it
    // about 1e-4 mm. On 0.5 mm cells split round it, the ligament lay in
    // one row of samples, each of whose material was taken at its
    // centroid: a layer with no thickness, no bending stiffness, and 0.14
    // mm. Each sample's material now spreads as it lies, a slab bending as
    // t^3 / 12: split as far as the floor, the bar reads near the fine
    // cells', or says the ligament is left thinner than its cells.
    let src = "(body bar (shape (subtract (move (box 20 4 4) 10 0 0) (move (box 0.4 8 3.95) 5.2 0 0.025))) (attr material \"steel\"))
(study s (shape (shape-of bar)) (elastic (material steel)) (fix root (- x 1)) (load tip (force 0 0 -0.00001N) (region (- 19 x))))";
    let tip = |s: &sim::Solution| -s.value_at([19.5, 0.0, 0.0]).unwrap()[2];
    let split = solve(src, "s", 0.5);
    assert!(split.report.cells_refined > 0, "{:?}", split.report.warnings);
    assert!(split.report.thin_left.is_some() || tip(&split) > 5e-5 && tip(&split) < 2e-4, "{}", tip(&split));
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_sheet_across_the_axes_bends_as_one_along_them() {
    // An 8 x 1 x 0.05 mm sheet cantilevered with a tip load, along the axes
    // and turned 45 degrees about x with its load: on 0.5 mm cells it bends
    // 4.88e-4 mm along them. Turned, each sample's material had been spread
    // along each axis alone, as thick as the sample across the sheet:
    // 7.96e-5, six times too stiff. With its covariance it lies as a sheet
    // (the rest, 1.5 times, is the coarse cells' 35% too much volume).
    let sheet = |turn: f64| {
        let (s, c) = (turn.to_radians().sin(), turn.to_radians().cos());
        format!(
            "(def part (rotate-x (move (box 8 1 0.05) 4 0 0) {turn}))
(study s (shape part) (elastic (E 200000MPa) (nu 0)) (fix root (- x 1)) (load tip (force 0 {}N {}N) (region (- 7 x))))",
            1e-5 * s,
            -1e-5 * c
        )
    };
    let tip = |turn: f64| solve(&sheet(turn), "s", 0.5).max_nodal().unwrap().value;
    let (along, across) = (tip(0.0), tip(45.0));
    assert!(across > 0.5 * along && across < 1.5 * along, "along {along}, turned {across}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn grids_split_to_the_same_cells_check_no_claim() {
    // The thin ligament of `a_thin_ligament_keeps_its_bending_in_a_cut_cell`
    // under a slit 0.1 mm to the left: on 0.5 mm cells it is split to 0.25,
    // and the claim's grids, 0.25, 0.5 and 1 mm, all split to 0.25 there and
    // read it alike (6.11e-5 mm, grid error 1.3e-7): it held under 8e-5,
    // where cells of 0.05 mm bend it 1.09e-4. Grids that read the part on
    // the same cells check nothing; coarser grids split no finer than twice
    // the finer's, and the claim does not hold.
    let src = "(body bar (shape (subtract (move (box 20 4 4) 10 0 0) (move (box 0.4 8 3.95) 5.1 0 0.025))) (attr material \"steel\"))
(study s (shape (shape-of bar)) (elastic (material steel)) (fix root (- x 1)) (load tip (force 0 0 -0.00001N) (region (- 19 x))) (res 0.5))
(require stiff (max-displacement s < 0.00008))";
    let (_, c) = &checks(src)[0];
    assert_ne!(c.status, ink::claims::Status::Holds, "{c:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_coarsest_grid_that_misses_a_hinge_checks_no_claim() {
    // A 0.1 mm ligament under a slit, moved 0.09 mm along z: on 0.5 mm
    // cells nothing splits and the ligament reads all but rigid (3.4e-7
    // mm), on 0.25 and 0.125 9.3e-6 and 1.055e-5, on 0.0625 1.48e-5. The
    // three grids 0.125, 0.25, 0.5 converged at an order of 2.8, above
    // their elements' 2, and their grid error was taken as their last
    // difference: a claim under 1.2e-5 held. Faster than the elements
    // converge is not convergence: the error is the non-asymptotic 3
    // |x - x2|, and the claim does not hold. So too a stress claim on a
    // 0.8 mm ligament under 0.8 of its net section's bending stress (its
    // grids 0.5, 1 and 2 mm, the coarsest all but blind to it).
    let disp = "(body bar (shape (move (subtract (move (box 20 4 4) 10 0 0) (move (box 0.4 8 3.9) 4.9 0 0.05)) 0 0 0.09)) (attr material \"steel\"))
(study s (shape (shape-of bar)) (elastic (material steel)) (fix root (- x 1)) (load tip (force 0 0 -0.00001N) (region (- 19 x))) (res 0.5))
(require stiff (max-displacement s < 1.2e-05))";
    let stress = "(body bar (shape (subtract (move (box 20 4 4) 10 0 0) (move (box 0.4 8 3.2) 5.1 0 0.4))) (attr material \"steel\"))
(study s (shape (shape-of bar)) (elastic (material steel)) (fix root (- x 1)) (load tip (force 0 0 -0.00001N) (region (- 19 x))) (res 0.5))
(require strong (max-stress s < 0.00027MPa))";
    for src in [disp, stress] {
        let (name, c) = &checks(src)[0];
        assert_ne!(c.status, ink::claims::Status::Holds, "{name}: {c:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn thin_material_is_found_wherever_the_samples_fall() {
    // The 0.05 mm ligament of `a_thin_ligament_keeps_its_bending_in_a_cut_cell`
    // (cells of 0.1 mm bend it about 1e-4 mm), moved 0.04 mm along z, its
    // slit 0.1 mm to the right: the ligament now crosses the 0.125 mm
    // samples' centres, and a sample whose centre lies in material was
    // thick. The slit's sides joined through it, nothing was split, and
    // the 0.5 mm cells read 2.75e-7: 360 times too stiff, unwarned. With
    // the slit's centre on a grid plane (x 5.0) the cells were split once,
    // to 0.25, and still tied by a corner 0.2 mm from any material: 1.9e-6,
    // unwarned. Material is thin by the chord of material across it, and
    // air by the chord of air, wherever the samples fall, and a thin
    // feature left when splitting stops is said, its claims not decided.
    let bar = |slit: &str, dz: f64| {
        format!(
            "(body bar (shape (move (subtract (move (box 20 4 4) 10 0 0) (move (box 0.4 8 3.95) {slit} 0 0.025)) 0 0 {dz})) (attr material \"steel\"))
(study s (shape (shape-of bar)) (elastic (material steel)) (fix root (- x 1)) (load tip (force 0 0 -0.00001N) (region (- 19 x))))"
        )
    };
    for (slit, dz) in [("5.1", 0.04), ("5.0", 0.0)] {
        let s = solve(&bar(slit, dz), "s", 0.5);
        let tip = s.max_nodal().unwrap().value;
        let said = s.report.thin_left.is_some();
        assert!(said || tip > 4e-5, "slit at {slit}, moved {dz}: {tip} mm, split {}, {:?}", s.report.cells_refined, s.report.warnings);
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn material_thinner_than_the_floor_is_decided_only_beside_resolved_material() {
    // Thin material is split until it is two cells thick, to an eighth of
    // a cell. Left thinner at that floor, it is decided only within two
    // floor cells of resolved material: a taper's edge, or a stub, inside
    // the floor's own resolution, its error the grid's. On 0.5 mm cells
    // (floor 0.0625) on a 4 x 6 x 6 block held at x < 1:
    let part = |shape: &str, load: &str, claim: &str| {
        format!(
            "(body part (shape (min (move (box 4 6 6) 2 0 0) {shape})) (attr material \"steel\"))
(study s (shape (shape-of part)) (elastic (material steel)) (fix root (- x 1)) (load tip (force 0 0 -1e-5N) (region {load})) (res 0.5))
(require a {claim})"
        )
    };
    let status = |src: &str| checks(src).remove(0).1.status;
    // A fin 0.05 mm thick, 1 wide and 16 long, and one turned 45 degrees
    // about its length (a sheet across the cells, which they read up to 3
    // times too stiff): thinner than the floor beyond its reach. Unknown.
    // (Its width only multiplies the cells along it: 4 wide, the same.)
    for fin in ["(move (box 20 1 0.05) 10 0 0)", "(rotate-x (move (box 21 1 0.05) 13.5 0 0) 45)"] {
        let src = part(fin, "(- 19 x)", "(max-displacement s < 1)");
        assert_eq!(status(&src), ink::claims::Status::Unknown, "{fin}");
    }
    // A fin tapering from 1 mm at its root to a knife edge 16 mm out: it is
    // thinner than the floor's two cells over its last 2 mm, 32 floor cells
    // from resolved material. Unknown: a slender wedge is a thin member.
    let wedge = "(max (move (box 20 1 1) 10 0 0) (- (abs z) (* 0.03125 (- 20 x))))";
    assert_eq!(status(&part(wedge, "(- 15 x)", "(max-displacement s < 1)")), ink::claims::Status::Unknown);
    // A stub 0.05 mm thick standing 0.1 mm off the block, the block's own
    // tip loaded: within the floor's reach of the block. Decided: the
    // block's tip moves 1.85e-11 mm, far under the bound. (On 0.25 mm cells,
    // floor 0.03125, the stub still under two of them: on 0.5 its grids read
    // 1.696e-11, 1.746e-11 and 1.813e-11, growing, not settled; with the
    // hold's cut cell gripped whole, the coarsest had read 1.006e-11.)
    let stub = "(move (box 4.1 4 0.05) 2.05 0 3.025)";
    let src = part(stub, "(max (- 3.5 x) (- x 4))", "(max-displacement s < 1e-6)").replace("(res 0.5)", "(res 0.25)");
    assert_eq!(status(&src), ink::claims::Status::Holds);
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn split_cells_carry_a_uniform_strain_exactly() {
    // A bar slit 10 mm along its length as the cantilever is, its ends
    // held at the uniform stretch u = 0.001 x (along x alone): its cells
    // along the slit are split, and their corners inside the faces of the
    // cells beside them hang on those faces. The stretch is in what they
    // can represent, so they carry it as it is: 1 MPa everywhere, the
    // displacement exact to the solver's tolerance (1e-9 mm of 4e-3; a
    // corner hung on the wrong nodes errs by the strain over a cell, 1e-3).
    let src = "(body bar (shape (subtract (box 100 10 10) (move (box 10 0.4 9) 0 1.124 0.5))))
(study pull (shape (shape-of bar)) (elastic (E 1000MPa) (nu 0.3))
  (fix a (+ x 49.5) (dofs x) (move (* 0.001 x) 0 0)) (fix b (- x 49.5) (dofs x) (move (* 0.001 x) 0 0)))";
    let s = solve(src, "pull", 1.004);
    assert!(s.report.cells_refined > 0 && s.report.thin_left.is_none(), "{:?}", s.report.warnings);
    for p in [[3.0, 1.0, 3.0], [3.0, 1.5, 3.0], [0.3, 0.9, -3.9], [0.3, 1.4, -4.2], [-4.0, 1.3, -4.1]] {
        let (u, (sigma, _, _)) = (s.value_at(p).unwrap(), s.field_at_point(p).unwrap());
        assert!((u[0] - 0.001 * p[0]).abs() < 1e-8 && (sigma[0] - 1.0).abs() < 1e-5, "{p:?}: u {u:?}, stress {sigma:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_gap_the_cells_glue_shut_decides_no_claim() {
    // A 0.4 mm saw slit through 9 of a 10 mm cantilever's depth leaves a
    // 1 mm ligament: 6 M / (b t^2) = 4770 MPa there and a 5.7 mm tip
    // deflection. Cells of 2 and 4 mm glued it shut (it is under half a
    // cell), solving SF 4.6 and 0.19 mm. It is air thinner than two cells:
    // the cells along it are split, and their samples part it, so the gap
    // check, which reads the cells the solve uses, finds nothing glued;
    // still thinner than two of the finest, it is said left.
    for h in [2.0, 4.0] {
        let s = solve(SLIT, "def", h);
        assert!(s.report.glued.is_none() && s.report.cells_refined > 0 && s.report.thin_left.is_some(), "{h} mm: {:?}", s.report.warnings);
    }
    // A 0.05 mm slit, narrower than the split cells' samples (0.0625 at 2
    // mm), is glued on them, and the solve says so: on the slit, above its
    // ligament, about its width.
    let s = solve(&SLIT.replace("(box 0.4 20 9)", "(box 0.05 20 9)"), "def", 2.0);
    let g = s.report.glued.unwrap_or_else(|| panic!("no gap found glued: {:?}", s.report.warnings));
    assert!(g.0 >= 2 && g.2 < 0.1 && (g.1[0] - 20.2).abs() < 0.5 && g.1[2] > -4.0, "{g:?}");
    assert!(s.report.warnings.iter().any(|w| w.contains("glued shut") && w.contains("split to")), "{:?}", s.report.warnings);
    // A gap the samples part (0.6 at 1 mm cells, `a_slot_narrower_than_a_cell_parts_its_sides`) glues nothing.
    let slot = "(body beam (shape (max (box 100 10 10) (neg (max (- (abs z) 0.3) (neg x))))))
      (study tip (shape (shape-of beam)) (elastic (E 1000MPa) (nu 0))
        (fix wall (+ x 49.5)) (load tip (force 0 0 -10N) (region (- 49 x))))";
    assert!(solve(slot, "tip", 1.0).report.glued.is_none());
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_bar_cut_through_by_a_hairline_gap_decides_no_claim() {
    // A 20 x 8 x 8 steel bar cut clean through by a 0.01 mm gap at x 5.01:
    // two bodies, the tip one unsupported. Boxes straddling the gap whose
    // corners and centre all lay in material were taken as solid, joining
    // its sides, and both claims held on 1, 2 and 4 mm cells (0.00004 mm,
    // safety factor 1356). Only boxes bounds prove solid join material, so
    // the gap is glued shut, and neither claim holds.
    let src = "(body bar (shape (subtract (move (box 20 8 8) 10 0 0) (move (box 0.01 16 16) 5.01 0 0))) (attr material \"steel\"))
(study s (shape (shape-of bar)) (elastic (material steel)) (fix root (- x 1)) (load tip (force 0 0 -1N) (region (- 19 x))) (res 1))
(require stiff (max-displacement s < 0.1))
(require strong (min-safety-factor s > 2))";
    for (name, c) in checks(src) {
        assert_eq!(c.status, ink::claims::Status::Unknown, "{name}: {c:?}");
    }
    let g = solve(src, "s", 1.0).report.glued.expect("the gap is found");
    assert!((g.1[0] - 5.01).abs() < 0.3, "{g:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_bar_cut_through_by_a_hair_holds_no_claim_and_flush_pieces_still_decide() {
    // A 20 x 4 x 4 steel bar cut clean through, its tip unsupported: no
    // static answer exists. A gap no point along the step fell in (1e-4,
    // 1e-6 mm; oblique; two flush boxes moved 1e-6 apart; boxes 1e-6 to
    // 1e-3 apart under a smooth union whose blend does not close them; a
    // cut bar inside a smooth union) was taken as joined, and stiffness
    // and strength held (0.00004 mm, safety factor 1356, on a bar 8 x 8:
    // the section only multiplies the cells across the cut). The sides
    // join only where bounds show them solid across; no air found proves
    // nothing, so none holds.
    let study = |part: &str| {
        format!(
            "(def part {part})
(study s (shape part) (elastic (material steel)) (fix root (- x 1)) (load tip (force 0 0 -1N) (region (- 19 x))) (res 1))
(require stiff (max-displacement s < 0.1))
(require strong (min-safety-factor s > 2))"
        )
    };
    let bar = "(move (box 20 4 4) 10 0 0)";
    let mut parts: Vec<String> = Vec::new();
    for w in ["0.001", "0.0001", "0.000001"] {
        parts.push(format!("(subtract {bar} (move (box {w} 40 40) 5.01 0 0))"));
        parts.push(format!("(smooth-union (subtract {bar} (move (box {w} 40 40) 5.01 0 0)) (move (box 1 1 1) 1 0 0) 0.1)"));
        parts.push(format!("(smooth-union (move (box 10 4 4) 5 0 0) (move (box 10 4 4) (+ 15 {w}) 0 0) {w})"));
    }
    for half in ["0.0005", "0.00005", "5e-07"] {
        parts.push(format!("(subtract {bar} (/ (- (abs (- (+ x (* 0.3 y) (* 0.2 z)) 5.01)) {half}) 1.063014581))"));
    }
    parts.push("(min (move (box 10 4 4) 5 0 0) (move (box 10 4 4) 15.000001 0 0))".into());
    for part in &parts {
        for (name, c) in checks(&study(part)) {
            assert_ne!(c.status, ink::claims::Status::Holds, "{part} {name}: {c:?}");
        }
    }
    // Pieces meeting exactly flush are one, by the kernel's rule: the two
    // boxes, and posts set into a plate, decide (the posts on 0.25 mm cells:
    // on 0.5 their grids read 1.471e-6, 1.551e-6 and 1.634e-6 mm, growing).
    for (name, c) in checks(&study("(min (move (box 10 4 4) 5 0 0) (move (box 10 4 4) 15 0 0))")) {
        assert_eq!(c.status, ink::claims::Status::Holds, "flush {name}: {c:?}");
    }
    let posts = "(def posts (min (move (box 2 2 6) -3 -2 3) (move (box 2 2 6) -3 2 3) (move (box 2 2 6) 3 -2 3) (move (box 2 2 6) 3 2 3)))
(study s (shape (min posts (box 12 8 2))) (elastic (material steel)) (fix root (- z -0.5)) (load tip (force 0 0 -1N) (region (- 5 z))) (res 0.25))
(require stiff (max-displacement s < 0.1))";
    assert_eq!(checks(posts)[0].1.status, ink::claims::Status::Holds, "posts");
}

#[test]
fn corners_and_shallow_grooves_glue_no_gap() {
    // What the cells join near a corner the material joins too, round it:
    // an L's inside corner, a 30 degree notch's root, and a groove 0.4 wide
    // and 1 deep, each under a cell wide where its walls meet, at cells of
    // 2 and 4 mm. Their error is the grid's, which shrinks with the
    // cells (a claim's coarser grids measure it).
    let shapes = [
        "(min (move (box 60 10 10) 30 0 0) (move (box 10 10 40) 5 0 20))",
        "(max (move (box 60 10 20) 30 0 0) (neg (- (* 0.966 (abs (- x 30))) (* 0.259 z))))",
        "(max (move (box 60 10 10) 30 0 0) (neg (max (- (abs (- x 30.2)) 0.2) (- 4 z))))",
    ];
    for shape in shapes {
        let src = format!(
            "(body part (shape {shape}) (attr material \"steel\"))
(study def (shape (shape-of part)) (elastic (material steel)) (fix root (- x 1)) (load tip (force 0 0 -10N) (region (- 59 x))))"
        );
        for h in [2.0, 4.0] {
            assert!(solve(&src, "def", h).report.glued.is_none(), "{shape} at {h} mm: {:?}", solve(&src, "def", h).report.glued);
        }
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_peak_on_a_curve_its_cells_do_not_resolve_is_unknown_until_they_do() {
    use ink::claims::Status::{Holds, Unknown};
    // A bore of radius 1 (written (1 - r^2)^3, flat on its surface) in an
    // 80 x 40 x 4 plate pulled by 160 N: CalculiX 2.95 MPa. On 1 mm cells
    // it read 1.71, its grids within 0.04: `< 2.9` held and `> 2.9` failed,
    // both false. Its radius spans one cell, under the 1.45 a peak's read
    // needs: unknown either way, saying so and the cells at which its
    // coarsest grid, four times as coarse, resolves it (0.17 mm).
    let bore = "(def part (max (box 80 40 4) (pow (- 1 (len2 x y)) 3)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa))
 (fix wall (+ x 39.9)) (load pull (force 160N 0 0) (region (- 39.9 x))) (res 1))
(require under (max-stress s < 2.9MPa))
(require over (max-stress s > 2.9MPa))";
    for (name, c) in checks(bore) {
        assert!(c.status == Unknown && undecided(&c).contains("curved to a radius of 0.99") && undecided(&c).contains("refine (res r) to 0.17"), "{name}: {c:?}");
    }
    // The same bore in a 16 x 8 x 2 plate: unknown on 1 mm cells (it read
    // 17.5, under half its peak), naming 0.17 mm, at which its three grids
    // span its radius 5.8, 2.9 and 1.45 times; on 0.125 mm cells it is read
    // (32.6, against Peterson's 32.4) and decided.
    let small = |res: &str| {
        format!(
            "(def part (max (box 16 8 2) (- 1 (len2 x y))))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa))
 (fix wall (+ x 7.9)) (load pull (force 160N 0 0) (region (- 7.9 x))) (res {res}))
(require under (max-stress s < 60MPa (excluding support-margin)))
(require over (max-stress s > 20MPa (excluding support-margin)))"
        )
    };
    for (res, want) in [("1", Unknown), ("0.125", Holds)] {
        for (name, c) in checks(&small(res)) {
            assert!(c.status == want && (want == Holds || undecided(&c).contains("refine (res r) to 0.17")), "{res}, {name}: {c:?}");
        }
    }
    // A bar in plain tension, 10 MPa, its faces flat: decided both ways
    // (read 10.08, its grid error 1.36: three grids oscillating, 10.47 on
    // 4.02 mm cells, whose last holds its load in a quarter of it (as when
    // held by penalties), and 10.02 on 2.01; and its read's error, 0.41 over
    // and 1.01 under). (Its 2 mm grid read 10.15 when the clamp's penalty
    // stiffened the cells it cut: a grid error of 0.23, and bounds of 11 and
    // 8.5.) On 4 mm cells on its faces' planes it reads 10 beside its wall
    // (its corner cells' free slivers lent to the cells beside them were cut
    // from the rest: 10.04 and 10.05).
    let bar = format!("{BAR}\n(require under (max-stress pull < 12MPa (excluding support-margin)))\n(require over (max-stress pull > 7.5MPa (excluding support-margin)))");
    for (name, c) in checks(&bar) {
        assert_eq!(c.status, Holds, "{name}: {c:?}");
    }
    let s = solve(BAR, "pull", 4.0);
    for p in [[-48.0, 0.0, 0.0], [-42.0, 0.0, 0.0], [-36.6, -3.2, 3.2]] {
        let (_, vm, _) = s.field_at(p).unwrap();
        assert!((vm - 10.0).abs() < 1e-4, "{vm} MPa at {p:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_bonded_rings_safety_factor_is_not_decided_across_its_curved_interface() {
    // A quarter ring, radii 10 to 14, bonded at 12 (E 200000 and 100000,
    // nu 0, yields 250 and 50), 2 MPa in its bore: Lame gives a least
    // safety factor of 14.443. Its grids read 12.80 at the interface, and
    // agree to 1.0: `< 14.154` held, false. Its interface, curved through
    // its cells, is a stair of them on every grid: neither way is decided.
    let src = "(def r (len2 x y))\n(def part (max (extrude-z (ring-2d 10 14) 0 4) (neg x) (neg y)))
(study s (shape part) (elastic (E (* 1MPa (select (- r 12) 200000 100000))) (nu 0) (yield (* 1MPa (select (- r 12) 250 50))))
 (fix sx (- x 0.001) (dofs x)) (fix sy (- y 0.001) (dofs y)) (fix sz (- z 0.001) (dofs z))
 (load inside (pressure 2MPa) (region (max (- r 10.05) (- 0.001 x) (- 0.001 y) (- 0.001 z) (- z 3.999)))) (res 0.25))
(require under (min-safety-factor s < 14.1541349454))
(require over (min-safety-factor s > 14.1541349454))";
    for (name, c) in checks(src) {
        let sf = measured(&c).unwrap();
        assert!(c.status == ink::claims::Status::Unknown && undecided(&c).contains("not flat along the cells") && (12.0..14.0).contains(&sf), "{name}: {sf}: {c:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn the_safety_factor_takes_each_materials_own_yield() {
    // Two bodies alike in stiffness, yields 250 and 50, pulled at 40 MPa
    // throughout: the least safety factor is the weaker's, 50 / 40 = 1.25
    // (one yield for the study read 250 / 40). Steel bonded to PLA (a
    // bar of both) is weaker still at its interface's edge. Each fails
    // 1.5, but every grid splits to the same cells at the interface, so no
    // three grids check it: unknown.
    let study = |e: &str, nu: &str| {
        format!(
            "(body a (shape (move (box 50 10 10) 25 0 0)) (attr E 3500MPa) (attr nu 0.3) (attr yield 250MPa))
(body b (shape (move (box 50 10 10) 75 0 0)) (attr E 3500MPa) (attr nu 0.3) (attr yield 50MPa))
(assembly root (instance sa a) (instance pb b))
(study one (shape (shape-of root)) (elastic (E (* 1MPa {e})) (nu {nu}) (yield (by-body root yield)))
  (fix sx (- x 0.01) (dofs x)) (fix sy (- y -4.99) (dofs y)) (fix sz (- z -4.99) (dofs z))
  (load pull (force 4000N 0 0) (region (- 99 x))) (res 1))
(require sf (min-safety-factor one > 1.5))"
        )
    };
    let (_, c) = &checks(&study("3500", "0.3"))[0];
    let sf = measured(c).unwrap();
    // (no grid error: no three grids read it)
    assert!(c.status == ink::claims::Status::Unknown && grid_error(c).is_none() && (sf - 1.25).abs() < 0.02 * 1.25, "{sf}: {c:?}");
    let steel_pla = study("(select (- x 50) 200000 3500)", "(select (- x 50) 0.3 0.36)");
    let (_, c) = &checks(&steel_pla)[0];
    assert!(c.status == ink::claims::Status::Unknown && grid_error(c).is_none() && measured(c).unwrap() < 1.25, "{c:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_material_layer_thinner_than_the_cells_is_split_to_or_left_undecided() {
    // A 20 x 2 x 2 conductor (k 200) held at 20 C over its first 0.01 mm, 0.2
    // W in its last millimetre (as 1 W through a 20 x 10 x 2 one), with an
    // insulating layer (k 0.2), t mm, across it at x 9.07 (or 9.013): exactly
    // 1-D, 20 + 1 W x (t / (0.2 x 20e-6) + (18.99 - t) / (200 x 20e-6) +
    // 0.125 K/W), 49.85 C where t is 0.1. One coefficient per cell lost the
    // layer on every grid of the claim's ladder alike: at its default 0.2 mm
    // cells the 0.1 mm layer read 25.5 C, `< 35` held with a grid error of 2
    // and `> 35` failed; the 0.02 and 0.04 mm layers held `< 27` and `< 30`
    // at 24.9 (truths 29.87, 34.86). Split until the layer is cells of its
    // own, the 0.1 mm one is read (49.87); the thinner ones, thinner than
    // two of the finest cells (0.025 mm), are left: their holds unknown.
    // Its faces lie within half a cell of planes its two finer grids share:
    // a cell taking its centre's k placed them there alike (0.0125 mm of the
    // layer is 3 K), and one weighing its k over its samples read it 43.6, a
    // layer across the heat's path as stiff as one along it. Each cell on a
    // face its two materials in series, it reads 49.84, decided both ways.
    let heat = |lo: f64, hi: f64, lim: f64| {
        format!(
            "(def part (move (box 20 2 2) 10 0 0))
(def pad (max (- {lo} x) (- x {hi})))
(study hot (shape part) (heat (k (* 1W/m/K (select pad 0.2 200))))
 (temperature cold (- x 0.01) 20degC)
 (source chip (power 0.2W) (region (- 19 x))))
(require cool (max-temperature hot < {lim}degC))
(require warm (max-temperature hot > {lim}degC))"
        )
    };
    let truth = |t: f64| 20.0 + t / (0.2e-3 * 20.0) + (18.99 - t) / (0.2 * 20.0) + 0.125;
    use ink::claims::Status::{Fails, Holds, Unknown};
    let c = checks(&heat(9.07, 9.17, 35.0));
    let v = measured(&c[0].1).unwrap();
    assert!((v - truth(0.1)).abs() <= grid_error(&c[0].1).unwrap() && c[0].1.status == Fails && c[1].1.status == Holds, "{v}: {c:?}");
    for (lo, hi, lim) in [(9.07, 9.09, 27.0), (9.013, 9.053, 30.0)] {
        let (c, t) = (checks(&heat(lo, hi, lim)), hi - lo);
        assert!(c[0].1.status != Holds && c[1].1.status != Fails, "{t}: {c:?}");
        assert!(c[0].1.status != Unknown || format!("{c:?}").contains("not of one material"), "{t}: {c:?}");
    }
    // A stiffer, weaker strip (E 600000, yield 25, in E 200000, yield 250)
    // 0.2 mm thick inside a 1 mm cell of a bar pulled along it shares the
    // bar's strain: 29.07 MPa in it, a factor of 0.86. The cell's one E
    // held `> 1.5` at 24.1.
    let (_, c) = &checks(
        "(def part (box 40 10 10))
(def strip (max (- 0.1 z) (- z 0.3) (- (abs y) 4) (- (abs x) 15)))
(study s (shape part)
 (elastic (E (* 1MPa (select strip 600000 200000))) (nu 0.3) (yield (* 1MPa (select strip 25 250))))
 (fix root (+ x 19)) (load pull (force 1000N 0 0) (region (- 19 x))) (res 1))
(require safe (min-safety-factor s > 1.5))",
    )[0];
    assert!(c.status != Holds, "{c:?}");
    assert!(c.status == Unknown || (measured(c).unwrap() - 0.86).abs() < 0.1, "{c:?}");
}

#[test]
fn a_peak_on_an_edge_oblique_to_the_cells_is_read_on_the_edge() {
    // A square bar of side a turned about its axis, bent by its ends' moves
    // (nu 0): sxx = E k (z - z0) at every point, exactly, so its peak is
    // E k |z - z0| at the edge farthest from z0. The surface's mesh at half
    // a cell cuts each edge's corner off, and its facets' centroids lay
    // 0.035 mm (0.14 of 0.25 mm cells) inside the edge: 16.80 read for
    // 17.505, 4% low alike on 0.25, 0.5 and 1 mm cells; and a bar bent
    // about its top edge, its bottom edge's tip in cells under half
    // material that no read reached, 32.77 for 33.52. The read climbs the
    // surface onto the edge (`climbed_scalar`), and a claim's read counts
    // what its field reaches there.
    let bar = |deg: f64, a: f64, z0: f64| {
        let (c, s) = (deg.to_radians().cos(), deg.to_radians().sin());
        let ctr = a;
        format!(
            "(def u (+ (* {c} (- y {ctr})) (* {s} (- z {ctr}))))
(def v (- (* {c} (- z {ctr})) (* {s} (- y {ctr}))))
(study s (shape (max (- (abs u) {h}) (- (abs v) {h}) (neg x) (- x 14))) (elastic (E 200000MPa) (nu 0))
 (fix a (- x 0.001) (move (* 1e-4 x (- z {z0})) 0 (* -5e-5 x x)))
 (fix b (- 13.999 x) (move (* 1e-4 x (- z {z0})) 0 (* -5e-5 x x))))",
            h = a / 2.0
        )
    };
    for (deg, a, z0, res) in [(37.0, 1.25, 1.25, 0.25), (37.0, 1.25, 1.25, 0.5), (36.0, 1.2, 1.2 + 0.5 * 1.2 * (36f64.to_radians().cos() + 36f64.to_radians().sin()), 0.5)] {
        let (c, s) = (f64::to_radians(deg).cos(), f64::to_radians(deg).sin());
        let edge = a / 2.0 * (c + s);
        let exact = 20.0 * ((a + edge - z0).abs()).max((a - edge - z0).abs());
        let sol = solve(&bar(deg, a, z0), "s", res);
        let (read, peak) = (sol.max_cell_scalar().unwrap(), sol.climbed_scalar().unwrap());
        assert!(
            (peak.value / exact - 1.0).abs() < 0.005 && read.value < 0.99 * exact,
            "{deg} deg, a {a}, z0 {z0}, res {res}: {} for {exact} at {:?}, read {}",
            peak.value,
            peak.at,
            read.value
        );
    }
}

#[test]
fn an_interface_a_cell_off_every_grids_planes_is_held_in_series_by_its_cell() {
    // A 20 x 2 x 2 bar held at 20 C over its first 0.01 mm, 0.2 W in its last
    // mm (as 1 W through 20 x 10 x 2), k 0.2 to x 10.03 and 200 beyond:
    // exactly 1-D, 20 + 2505 + 2.2425 + 0.125 = 2527.37 C. A cell taking its
    // centre's k read it 7.4 K low on cells of 0.1, 0.2 and 0.4 mm alike (the
    // 0.03 mm of k 0.2 past x 10 as k 200), and `< 2526`, false, held. The
    // cell holding the face is its two materials in series (a laminate): on
    // every grid the bar's resistance, within its grid error (2527.367,
    // 2527.367 and 2527.345 on 0.1, 0.2 and 0.4 mm cells: 0.066, three
    // times the coarsest's 0.022 at its heated end, whose source's edge cuts
    // a cell there alone); on the two finer grids, 2527.3675 to 1e-4. (Its
    // held band, 0.4 of the first samples deep, split by points in samples a
    // face left 1 - 1e-15 whole, read 0.375: 2527.45 on every grid, 0.08 K
    // over, and the grid error under that, 0.036, saw none of it.)
    let src = "(def part (move (box 20 2 2) 10 0 0))
(study hot (shape part) (heat (k (* 1W/m/K (select (- x 10.03) 0.2 200))))
 (temperature cold (- x 0.01) 20degC)
 (source chip (power 0.2W) (region (- 19 x))) (res 0.1))
(require cool (max-temperature hot < 2526degC))";
    let (_, c) = &checks(src)[0];
    let v = measured(c).unwrap();
    let e = grid_error(c).unwrap_or(f64::INFINITY);
    assert!(c.status == ink::claims::Status::Fails && (v - 2527.3675).abs() < e && e < 0.1, "{v}: {c:?}");
    for res in [0.1, 0.2] {
        let t = solve(src, "hot", res).max_nodal().unwrap().value;
        assert!((t - 2527.3675).abs() < 1e-3, "cells {res}: {t} C");
    }
}

#[test]
fn a_hold_across_a_laminates_face_reads_its_materials_in_series() {
    // A 20 x 2 x 2 bar held at 20 C below x 10.05, k 0.2 to x 10.13 and 200
    // beyond, 0.2 W in its last mm: exactly 1-D, 20 + 20 + 2.2175 + 0.125 =
    // 42.3425 C. On 0.2 and 0.4 mm cells the hold's edge and the face share
    // a cell, a laminate: its free part's stiffness is the laminate of its
    // own share of each material (of the cell's, 4.4 K hot on 0.2 mm cells),
    // and the hold's traction that stiffness's (of the cell's mean k,
    // isotropic, 28.5 K hot). The 0.4 mm cells' heated end reads 0.02 K low,
    // as when the face holds alone.
    for (res, within) in [(0.1, 1e-3), (0.2, 1e-3), (0.4, 0.03)] {
        let src = "(def part (move (box 20 2 2) 10 0 0))
(study hot (shape part) (heat (k (* 1W/m/K (select (- x 10.13) 0.2 200))))
 (temperature cold (- x 10.05) 20degC) (source chip (power 0.2W) (region (- 19 x))))";
        let t = solve(src, "hot", res).max_nodal().unwrap().value;
        assert!((t - 42.3425).abs() < within, "cells {res}: {t} C");
    }
}

/// A bar 24 x 10 x 10 (nu 0) pulled 0.01 mm by its ends, its E, and its
/// yield where it has one, as given; its claims after.
fn pulled_bar(e: &str, yield_: &str, res: f64, claims: &str) -> String {
    format!(
        "(def part (max (neg x) (- x 24) (- (abs (- y 5)) 5) (- (abs (- z 5)) 5)))
(study s (shape part) (elastic (E (* 1MPa {e})) (nu 0) (yield (* 1MPa {yield_})))
 (fix a (- x 0.001) (move 0 0 0))
 (fix b (- 23.999 x) (move 0.01 0 0)) (res {res}))
{claims}"
    )
}

#[test]
fn bonded_bars_off_every_grids_planes_read_their_materials_in_series_and_in_parallel() {
    use ink::claims::Status::{Fails, Holds};
    // In series: a soft band (E 2000) from x 8.2 to 15.76 in E 200000, its
    // faces 0.2 and 0.24 mm off the planes of every grid. Exactly 1-D: 0.01
    // / (16.44 / 200000 + 7.56 / 2000) = 2.58920 MPa throughout. A cell
    // taking its centre's E placed the faces on the planes x 8 and 16 on
    // all three grids, 2.4510 alike (5.3% low), and the claims were refused;
    // its cell the two materials in series, every grid reads the bar's.
    // In parallel: a soft, weak layer (E 2000, yield 2) under z 3.3 under a
    // stiff one (E 200000, yield 250): each strained 0.01 / 24, 83.33 MPa in
    // the stiff and 0.8333 in the soft, the least factor the soft's, 2.4.
    // Each on a claim's three grids (a peak stress read may lie
    // up to 10% over the truth and 22% under: their limits beyond it).
    let series = 0.01 / (16.44 / 200000.0 + 7.56 / 2000.0);
    let band = "(select (max (- 8.2 x) (- x 15.76)) 2000 200000)";
    let layer = |v: &str| format!("(select (- z 3.3) {v})");
    // (On 0.5 mm cells its three grids read 2.589191, 2.589202, 2.589219: a drift of 1e-5 growing, past
    // what the solves resolve, so no grid error and not settled; on 0.25 they converge.)
    let res = 0.25;
    let c = checks(&pulled_bar(
        band,
        "250",
        res,
        "(require over (max-stress s > 2.3MPa (excluding support-margin)))\n(require under (max-stress s < 2.3MPa (excluding support-margin)))",
    ));
    let v = measured(&c[0].1).unwrap();
    assert!(c[0].1.status == Holds && c[1].1.status == Fails && (v / series - 1.0).abs() < 2e-5, "series, {res}: {v} for {series}: {c:?}");
    assert!(grid_error(&c[0].1).is_some_and(|e| e < 1e-4 * series), "series, {res}: {c:?}");
    // (The layer, 3.3 mm, is split on 2 mm cells: its ladder starts at 0.25.)
    let res = 0.25;
    let c = checks(&pulled_bar(
        &layer("2000 200000"),
        &layer("2 250"),
        res,
        "(require stress-over (max-stress s > 74MPa (excluding support-margin)))
(require stress-under (max-stress s < 102MPa (excluding support-margin)))
(require sf-over (min-safety-factor s > 2.1 (excluding support-margin)))
(require sf-under (min-safety-factor s < 2.1 (excluding support-margin)))",
    ));
    let (stress, sf) = (measured(&c[0].1).unwrap(), measured(&c[2].1).unwrap());
    assert!((stress / (200000.0 * 0.01 / 24.0) - 1.0).abs() < 1e-3 && (sf / 2.4 - 1.0).abs() < 1e-3, "parallel, {res}: {stress} {sf}: {c:?}");
    // The cells holding the face bounded by each material's own stress
    // over its own yield: the soft's, the factor's, not under it.
    let bound = excluded_value(&c[2].1, "yield-between-reads").unwrap();
    assert!((bound / 2.4 - 1.0).abs() < 1e-3, "parallel, {res}: bound {bound}");
    assert!(c.iter().all(|(name, c)| c.status == if name == "sf-under" { Fails } else { Holds }), "parallel, {res}: {c:?}");
}

#[test]
fn a_laminate_cell_is_bounded_by_each_sides_least_yield_not_its_samples() {
    // The parallel layer (E 2000 yield 2 under z 3.3, E 200000 yield 250
    // over it) on 0.25 mm cells: the face's cell a laminate, its samples
    // each side's yield. Its region read holds more than they see (a skin of
    // a cell under half material beside it, or a yield's step between its
    // samples): here the fields' bounds say the stiff side yields at 100
    // from z 3.31 to 3.34, between its samples (z 3.281 and 3.344). At 83.33
    // MPa that is a least factor of 1.2, not the soft layer's 2.4; each side
    // took its samples' yield, 2 and 250, and bounded the cell at 2.4.
    let src = pulled_bar("(select (- z 3.3) 2000 200000)", "(select (- z 3.3) 2 250)", 0.25, "");
    let s = solve(&src, "s", 0.25);
    let mixed = |b: &[([f64; 3], [f64; 3])]| {
        b.iter()
            .map(|(lo, hi)| {
                let (soft, stiff, stripe) = (lo[2] < 3.3, hi[2] > 3.3, lo[2] < 3.34 && hi[2] > 3.31);
                let least = if soft {
                    2.0
                } else if stripe {
                    100.0
                } else {
                    250.0
                };
                (soft && stiff || stripe).then_some(sim::MixedCell { least, stiffest: if stiff { 200000.0 } else { 2000.0 } })
            })
            .collect::<Vec<_>>()
    };
    let factor = 1.0 / s.ratio_between_reads(false, &mixed).unwrap().value;
    assert!(factor <= 1.2 * (1.0 + 1e-3) && factor > 1.1, "{factor}");
}

#[test]
fn a_cell_not_of_one_material_is_bounded_by_its_stiffest() {
    // E 2000 (1 + z) across a bar pulled 0.01 mm along x (nu 0, yield 250):
    // every fibre strained 0.01 / 24, the stiffest, at z 10, at 9.167 MPa: a
    // least factor of 27.27. No cell is of one material, so each is bounded,
    // and the bound took the cell's own E at its centre: the top cell's,
    // half a cell under the top, read 28.57 on 1 mm cells, over the factor
    // the part has. Bounded by the stiffest E over each cell, it is the
    // factor on every grid.
    let exact = 250.0 / (22000.0 * 0.01 / 24.0);
    for res in [0.5, 1.0, 2.0] {
        let (_, c) = &checks(&pulled_bar("(* 2000 (+ 1 z))", "250", res, "(require safe (min-safety-factor s > 20 (excluding support-margin)))"))[0];
        let bound = excluded_value(c, "yield-between-reads").unwrap();
        assert!(bound <= exact * (1.0 + 1e-4) && bound > 0.95 * exact, "{res}: {bound} for {exact}: {c:?}");
    }
}

/// A 20 x 2 x 2 bar (x 0..20) of conductivity `k` held at 20 C over its
/// first 0.01 mm, 0.2 W in its last millimetre, on 0.2 mm cells: its claim
/// `cool`, `(max-temperature hot < LIMIT)`. Exactly one-dimensional: T =
/// 20 + the integral of 0.2 W over k A from the held band to the source,
/// and half the source's millimetre. (The temperatures of 1 W through a
/// 20 x 10 x 2 bar: its section only multiplies the cells across it.)
fn bar_heat(k: &str, limit: f64) -> ink::claims::Check {
    let src = format!(
        "(def part (move (box 20 2 2) 10 0 0))
(study hot (shape part) (heat (k (* 1W/m/K {k})))
 (temperature cold (- x 0.01) 20degC)
 (source chip (power 0.2W) (region (- 19 x))) (res 0.2))
(require cool (max-temperature hot < {limit}degC))"
    );
    checks(&src).remove(0).1
}

#[test]
fn a_layer_is_not_taken_as_an_interface_by_another_body_of_its_material() {
    // An insulating block (k 0.2) at x 8.8..10.4, 0.12 mm of the conductor
    // (k 200), then a 0.07 mm film of the block's material at x
    // 10.52..10.59: 400 + 17.5 + 4.33 + 0.125 K/W, 441.96 C (424.47 with
    // the film lost). The film's cell held the block's value beside it and
    // the conductor's, so it read as one interface: `< 435` held at 424.4
    // on every grid. A material held across a step must be the one that
    // reaches the step.
    let c = bar_heat("(select (min (max (- 10.52 x) (- x 10.59)) (max (- 8.8 x) (- x 10.4))) 0.2 200)", 435.0);
    assert!(c.status != ink::claims::Status::Holds, "{c:?}");
    assert!(c.status != ink::claims::Status::Fails || (measured(&c).unwrap() - 441.96).abs() < 0.03 * 441.96, "{c:?}");
}

#[test]
fn a_layer_on_a_gentle_slope_is_not_taken_as_the_slope() {
    // A k 1 layer 0.03 mm thick at x 10.405..10.435 on a conductivity
    // rising 90 to 110 (an arctangent, 0.1 mm wide, about x 10.42): the
    // cell's range barely grew from its quarter to its half to its whole,
    // read as a slope converging slowly, and `< 30.5` held at 29.90 on
    // every grid; the layer brings it to 31.39 (by the 1-D integral).
    let c = bar_heat("(select (max (- 10.405 x) (- x 10.435)) 1 (+ 100 (* 6.366197723675814 (atan (/ (- x 10.42) 0.1)))))", 30.5);
    assert!(c.status != ink::claims::Status::Holds, "{c:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_smooth_property_with_a_least_value_inside_is_resolved() {
    // k = 100 + (x - 10.1)^2 / 10: smooth, least at a cell's centre, 29.4517
    // C by the 1-D integral. The cell's largest half varied as much as the
    // cell, read as a step, split to the floor and left: unknown.
    let c = bar_heat("(+ 100 (/ (* (- x 10.1) (- x 10.1)) 10))", 70.0);
    assert!(c.status == ink::claims::Status::Holds && (measured(&c).unwrap() - 29.4517).abs() < 0.01 * 29.4517, "{c:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_film_in_a_cells_outer_part_is_not_hidden_by_what_varies_nearer_its_centre() {
    // An insulating film (k 0.1) 0.02 mm thick in the outer part of the
    // cell x 10.0..10.2, beside variation nearer the cell's centre: the
    // ranges about the cell's centre grew from its quarter to its half
    // (the nearer variation) and to the whole (the film) as a smooth
    // field's do, so the cell was taken as resolved, never split, and the
    // film lost on every grid: `< 30` held at 25.7, 24.9 and 25.0. Truths
    // by the 1-D series (or integral).
    use ink::claims::Status::{Holds, Unknown};
    let film = "(max (- 10.01 x) (- x 10.03))";
    let v = "(min 200 (+ 2 (* 990 (abs (- x 10.2)))))";
    for (what, k, truth) in [
        ("beside a joint of two metals", format!("(select {film} 0.1 (select (- x 10.06) 200 150))"), 35.654),
        ("beside a thinner film of k 150", format!("(select {film} 0.1 (select (max (- 10.06 x) (- x 10.07)) 150 200))"), 34.868),
        ("on a graded dip", format!("(select {film} 0.1 {v})"), 35.232),
    ] {
        let c = bar_heat(&k, 30.0);
        assert!(c.status != Holds, "{what}: {c:?}");
        assert!(c.status == Unknown || (measured(&c).unwrap() - truth).abs() < 0.05 * truth, "{what}: {c:?}");
    }
    // The dip alone holds (25.238 by the integral). Smooth in k, but to k 2
    // at 990 a mm: in 1/k, what cells in series sum, its bottom is 2 / 990
    // mm wide, under the finest cells, so it may be left: not failed.
    let c = bar_heat(v, 30.0);
    assert!(c.status != ink::claims::Status::Fails && (c.status == Unknown || (measured(&c).unwrap() - 25.238).abs() < 0.03 * 25.238), "{c:?}");
    // A film of the dip's own least value, 0.04 mm, moves no range of the
    // cell's: 26.226, so `< 26` fails.
    let c = bar_heat(&format!("(select (max (- 10.005 x) (- x 10.045)) 2 {v})"), 26.0);
    assert!(c.status != Holds, "{c:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn no_placement_of_a_film_and_a_step_beside_it_in_a_cell_hides_the_film() {
    // The film (k 0.1, 0.02 mm) and a thinner film of k 150 (0.01 mm) 0.05
    // mm further on, moved across the cell x 10.0..10.2 by 0.01 mm, and
    // mirrored (the step nearer the centre on the other side): `< 30`
    // fails at every placement (34.87 by the 1-D series), never holds.
    let case = |f: f64, d: f64| bar_heat(&format!("(select (max (- {f} x) (- x {})) 0.1 (select (max (- {d} x) (- x {})) 150 200))", f + 0.02, d + 0.01), 30.0);
    let placements = (0..20).map(|i| (10.0 + 0.01 * i as f64, 10.05 + 0.01 * i as f64)).chain([10.16, 10.17, 10.175, 10.18].into_iter().flat_map(|f| [(f, 10.13), (f, 10.14)]));
    for (f, d) in placements {
        let c = case(f, d);
        assert!(c.status != ink::claims::Status::Holds, "film {f}, step {d}: {c:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_property_dipping_toward_zero_is_judged_as_its_cells_average_it() {
    // k = min(200, k0 + G |x - X0|): smooth, but in series the cells
    // average 1/k, which is near singular at its bottom (0.1 to 100 over
    // one cell where k0 is 0.01). Read as k, the bottom cell passed as
    // resolved, unsplit, and its one value overshot: `< 36` failed at 37.87
    // against 33.776 (a fail was then not checked on coarser grids). Truths by
    // 2 ln(200 / k0) / G over the dip, in series with the rest.
    let wrong = |c: &ink::claims::Check, holds: bool| c.status == if holds { ink::claims::Status::Fails } else { ink::claims::Status::Holds };
    for (name, k, limit, holds) in [
        ("dip-10.0125", "(min 200 (+ 0.01 (* 100 (abs (- x 10.0125)))))", 30.0, false),
        ("dip-10.037", "(min 200 (+ 0.01 (* 100 (abs (- x 10.037)))))", 30.0, false),
        ("dip-10.05", "(min 200 (+ 0.01 (* 100 (abs (- x 10.05)))))", 30.0, false),
        ("dip-10.1-lt32", "(min 200 (+ 0.01 (* 100 (abs (- x 10.1)))))", 32.0, false),
        ("dip-10.1-lt36", "(min 200 (+ 0.01 (* 100 (abs (- x 10.1)))))", 36.0, true),
        ("dip-10.1", "(min 200 (+ 0.01 (* 100 (abs (- x 10.1)))))", 30.0, false),
        ("dip-k1-10.3", "(min 200 (+ 1 (* 100 (abs (- x 10.3)))))", 30.05, true),
        ("dipsweep/d-0.01-100-10.1", "(min 200 (+ 0.01 (* 100 (abs (- x 10.1)))))", 34.79, true),
        ("dipsweep/d-0.01-100-10.15", "(min 200 (+ 0.01 (* 100 (abs (- x 10.15)))))", 34.79, true),
        ("dipsweep/d-0.01-100-10.3", "(min 200 (+ 0.01 (* 100 (abs (- x 10.3)))))", 34.79, true),
        ("dipsweep/d-0.1-100-10.1", "(min 200 (+ 0.1 (* 100 (abs (- x 10.1)))))", 32.42, true),
        ("dipsweep/d-0.1-100-10.15", "(min 200 (+ 0.1 (* 100 (abs (- x 10.15)))))", 32.42, true),
        ("dipsweep/d-0.1-100-10.3", "(min 200 (+ 0.1 (* 100 (abs (- x 10.3)))))", 32.42, true),
        ("dipsweep/d-1-100-10.1", "(min 200 (+ 1 (* 100 (abs (- x 10.1)))))", 30.05, true),
        ("dipsweep/d-1-100-10.15", "(min 200 (+ 1 (* 100 (abs (- x 10.15)))))", 30.05, true),
        ("dipsweep/d-1-100-10.3", "(min 200 (+ 1 (* 100 (abs (- x 10.3)))))", 30.05, true),
        ("dipsweep/d-1-1000-10.1", "(min 200 (+ 1 (* 1000 (abs (- x 10.1)))))", 26.06, true),
        ("dipsweep/d-1-1000-10.15", "(min 200 (+ 1 (* 1000 (abs (- x 10.15)))))", 26.06, true),
        ("dipsweep/d-1-1000-10.3", "(min 200 (+ 1 (* 1000 (abs (- x 10.3)))))", 26.06, true),
        ("dipsweep/d-5-100-10.1", "(min 200 (+ 5 (* 100 (abs (- x 10.1)))))", 28.41, true),
        ("dipsweep/d-5-100-10.15", "(min 200 (+ 5 (* 100 (abs (- x 10.15)))))", 28.41, true),
        ("dipsweep/d-5-100-10.3", "(min 200 (+ 5 (* 100 (abs (- x 10.3)))))", 28.41, true),
    ] {
        let c = bar_heat(k, limit);
        assert!(!wrong(&c, holds), "{name}: {c:?}");
    }
    // E = min(200000, 200 + 1e5 |x - 10.1|) in a 2 x 2 bar pulled by 200
    // N, nu 0: 0.010781 mm by 50 MPa times the integral of 1/E.
    let (_, c) = &checks(
        "(def part (move (box 20 2 2) 10 0 0))
(study pull (shape part) (elastic (E (* 1MPa (min 200000 (+ 200 (* 100000 (abs (- x 10.1))))))) (nu 0))
 (fix root (- x 0.01)) (load end (force 200N 0 0) (region (- 19 x))) (res 0.2))
(require stretch (max-displacement pull < 0.0111))",
    )[0];
    assert!(!wrong(c, true), "{c:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_film_is_judged_by_its_resistance_not_only_its_value() {
    // A k 0.01 film on a slope of k 200 to 60, with a sliver of k 120
    // beside it in its own finest cell (in each, where it straddles two):
    // read as k, the sliver grew the box's range as a slope's grows, and
    // `< 35` held at 30.3 against 50.33. And a film at the bottom of a dip
    // whose least value is near its own moved no range. In 1/k the film is
    // a jump of 100 against a few thousandths beside it. Truths by the 1-D
    // series.
    let wrong = |c: &ink::claims::Check, holds: bool| c.status == if holds { ink::claims::Status::Fails } else { ink::claims::Status::Holds };
    for (name, k, limit, holds) in [
        (
            "boundary-film",
            "(select (max (- 10.023 x) (- x 10.027)) 0.01 (select (min (max (- 10.016 x) (- x 10.0175)) (max (- 10.032 x) (- x 10.033))) 120 (max 60 (min 200 (- 200 (* 700 (- x 10)))))))",
            35.0,
            false,
        ),
        ("decoy-next-box", "(select (max (- 10.001 x) (- x 10.005)) 0.01 (select (max (- 10.03 x) (- x 10.031)) 120 (max 60 (min 200 (- 200 (* 700 (- x 10)))))))", 35.0, false),
        ("floorbox-decoy", "(select (max (- 10.001 x) (- x 10.005)) 0.01 (select (max (- 10.007 x) (- x 10.008)) 120 (max 60 (min 200 (- 200 (* 700 (- x 10)))))))", 35.0, false),
        ("slope-film-control", "(select (max (- 10.001 x) (- x 10.005)) 0.01 (max 60 (min 200 (- 200 (* 700 (- x 10))))))", 35.0, false),
        ("v-film-10.0125", "(select (max (- 10.0115 x) (- x 10.0135)) 0.01 (min 200 (+ 0.5 (* 1000 (abs (- x 10.0125))))))", 28.0, false),
        ("v-film-10.037-lt38", "(select (max (- 10.036 x) (- x 10.038)) 0.01 (min 200 (+ 0.5 (* 1000 (abs (- x 10.037))))))", 38.0, true),
        ("v-film-10.037", "(select (max (- 10.036 x) (- x 10.038)) 0.01 (min 200 (+ 0.5 (* 1000 (abs (- x 10.037))))))", 28.0, false),
        ("v-film-10.05", "(select (max (- 10.049 x) (- x 10.051)) 0.01 (min 200 (+ 0.5 (* 1000 (abs (- x 10.05))))))", 28.0, false),
        ("v-film-10.1", "(select (max (- 10.099 x) (- x 10.101)) 0.01 (min 200 (+ 0.5 (* 1000 (abs (- x 10.1))))))", 28.0, false),
    ] {
        let c = bar_heat(k, limit);
        assert!(!wrong(&c, holds), "{name}: {c:?}");
    }
}

#[test]
fn a_smooth_property_is_resolved_where_a_step_is_not() {
    // A property's ranges over boxes a quarter, a half and the whole of a
    // cell about its centre, read as a claim's ladder reads three grids,
    // its order held to 1 to 2: a smooth field's grows so, a step's does
    // not.
    use ink::sim::variation_resolved as resolved;
    assert!(resolved(0.0, 0.0, 0.0), "one material");
    assert!(resolved(1.0, 2.0, 4.0), "affine");
    assert!(resolved(1.0, 4.0, 16.0), "least at the centre");
    assert!(resolved(0.14, 0.28, 0.6084), "least off the centre, past the half box");
    assert!(resolved(1.0, 2.05, 4.3) && resolved(1.0, 1.95, 3.7), "curved");
    assert!(!resolved(1.0, 1.0, 1.0), "a step");
    assert!(!resolved(0.76, 0.96, 1.0), "a step narrower than the cell");
    assert!(!resolved(10.0, 10.1, 10.2101), "a step on a slope growing slower than linearly");
    assert!(!resolved(0.25, 0.5, 5.0), "a step in the cell's outer part");
    assert!(!resolved(0.0, 0.0, 1.0), "a step in the cell's outer part, nothing within");
    assert!(!resolved(f64::INFINITY, f64::INFINITY, f64::INFINITY));
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_modal_study_with_a_symmetry_plane_decides_no_least_frequency() {
    // A steel strip 100 x 10 x 2 clamped at one end bends first through
    // its thickness at 164.7 Hz (Euler-Bernoulli). Modelled as half its
    // thickness on a symmetry plane, that mode is gone (antisymmetric
    // about it): its least is the sideways 815 Hz, and `> 300` held.
    let src = "(body full (shape (move (box 100 10 2) 50 0 0)) (attr material \"steel\"))
(body half (shape (move (box 100 10 1) 50 0 0.5)) (attr material \"steel\"))
(study fullm (shape (shape-of full)) (elastic (material steel)) (fix root (- x 0.5)) (modes 2) (res 0.5))
(study halfm (shape (shape-of half)) (elastic (material steel)) (fix root (- x 0.5)) (mirror sym z 0) (modes 2) (res 0.5))
(require f-full (min-frequency fullm > 300Hz))
(require f-half (min-frequency halfm > 300Hz))";
    let c = checks(src);
    assert!(c[0].1.status == ink::claims::Status::Fails && (measured(&c[0].1).unwrap() / 164.7 - 1.0).abs() < 0.05, "{:?}", c[0]);
    assert_eq!(c[1].1.status, ink::claims::Status::Unknown, "{:?}", c[1]);
    assert!(undecided(&c[1].1).contains("antisymmetric about the plane of symmetry 'sym'"), "{}", undecided(&c[1].1));
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_buckling_member_decides_no_displacement_claim() {
    // A fixed-free steel strip 100 x 10 x 1 pushed at 48 times its Euler
    // load (pi^2 E I / (4 L^2) = 41 N): its linear shortening, 0.1 mm,
    // held `< 0.2`. Its buckling load is found, and the claim is unknown,
    // as a stress claim is.
    let src = "(body s (shape (move (box 100 10 1) 50 0 0)) (attr material \"steel\"))
(study push (shape (shape-of s)) (elastic (material steel))
  (fix root (- x 0.5)) (load end (force -2000N 0 0) (region (- 99 x))) (res 0.25))
(require short (max-displacement push < 0.2))";
    let (_, c) = &checks(src)[0];
    assert_eq!(c.status, ink::claims::Status::Unknown, "{c:?}");
    assert!(undecided(c).contains("buckle"), "{}", undecided(c));
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_rod_on_a_squeezed_block_has_its_buckle_found() {
    // A 1.5 mm square rod 30 mm long out of a 10 mm cube held at y < 0.5,
    // pushed along its axis at twice its cantilever Euler load (it buckles
    // at 0.49 of it), the cube pressed down y. The cube's walls, too thick
    // to flag, joined the rod's member of a buckling screen; outnumbering
    // the rod's cells, they set its direction (the squeeze's) and the rod
    // was lost in one of its sections: unscreened, its claims held. Every
    // static solve now finds its buckling load: found at every squeeze.
    for press in [0.0, 1000.0, 5000.0, 20000.0] {
        let src = format!(
            "(def part (min (move (box 10 10 10) 5 5 0) (move (box 30 1.5 1.5) 25.0 5 0)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 1000MPa)) (fix base (- y 0.5))
  (load press (force 0 -{press}N 0) (region (- 9.5 y))) (load push (force -462.638N 0 0) (region (- 39.5 x))))"
        );
        let s = solve(&src, "s", 0.5);
        let b = s.report.buckle.as_ref().unwrap_or_else(|| panic!("pressed {press}: no buckling analysis: {:?}", s.report.warnings));
        assert!(b.failure.is_none() && b.low < 0.55 && b.factor > 0.45, "pressed {press}: {b:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_rod_and_a_post_out_of_one_block_have_the_rods_buckle_found() {
    // A 1.5 mm square rod 30 mm long out of a 10 mm cube's side, held at
    // y < 0.5, pushed along its axis at twice its cantilever Euler load
    // (it buckles at 0.49 of it), and a 4 x 4 post on the cube's top
    // pressed down. A buckling screen joined rod and post through the
    // cube into one member, whose axis the post's compression set: the
    // rod went unscreened. Every static solve now finds its buckling load:
    // found at every press on the post.
    for press in [0.0, 500.0, 1000.0, 1500.0] {
        let src = format!(
            "(def part (min (move (box 10 10 10) 5 5 0) (min (move (box 30 1.5 1.5) 25.0 5 0) (move (box 4 20 4) 5 20 0))))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 1000MPa)) (fix base (- y 0.5))
  (load pb (force 0 -{press}N 0) (region (- 29.5 y))) (load push (force -462.6N 0 0) (region (- 39.5 x))))"
        );
        let s = solve(&src, "s", 0.25);
        let b = s.report.buckle.as_ref().unwrap_or_else(|| panic!("post pressed {press}: no buckling analysis: {:?}", s.report.warnings));
        assert!(b.failure.is_none() && b.low < 0.55 && b.factor > 0.45, "post pressed {press}: {b:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn rods_in_line_through_thick_blocks_buckle_as_one_column() {
    // Rods 1.5 mm square in line through small thick blocks (chains, a
    // ladder, a boss), at once and a half to twice their buckling load: a
    // buckling screen judged each rod alone, stocky, and passed it. Every
    // static solve now finds its buckling load: each found near ink's on
    // finer cells.
    let cases: [(&str, f64, &str); 5] = [
        (
            "chain3_x",
            0.5005,
            r#"(def part (min (min (min (min (move (box 30 1.5 1.5) 15.0 0 0) (move (box 10 10 10) 35.0 0 0)) (move (box 30 1.5 1.5) 55.0 0 0)) (move (box 10 10 10) 75.0 0 0)) (move (box 30 1.5 1.5) 95.0 0 0)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 1000MPa)) (fix base (- x 0.5))
  (load push (force -42.1N 0 0) (region (- 109.5 x))) (res 0.5))"#,
        ),
        (
            "chain3_x15",
            0.6668,
            r#"(def part (min (min (min (min (move (box 30 1.5 1.5) 15.0 0 0) (move (box 10 10 10) 35.0 0 0)) (move (box 30 1.5 1.5) 55.0 0 0)) (move (box 10 10 10) 75.0 0 0)) (move (box 30 1.5 1.5) 95.0 0 0)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 1000MPa)) (fix base (- x 0.5))
  (load push (force -31.6N 0 0) (region (- 109.5 x))) (res 0.5))"#,
        ),
        (
            "ladder_x",
            0.5000,
            r#"(def part (min (min (min (min (min (min (min (min (move (box 25 1.5 1.5) 12.5 -4 0) (move (box 25 1.5 1.5) 12.5 4 0)) (move (box 25 1.5 1.5) 45.5 -4 0)) (move (box 25 1.5 1.5) 45.5 4 0)) (move (box 25 1.5 1.5) 78.5 -4 0)) (move (box 25 1.5 1.5) 78.5 4 0)) (move (box 8 12 10) 29.0 0 0)) (move (box 8 12 10) 62.0 0 0)) (move (box 8 12 10) 95.0 0 0)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 1000MPa)) (fix base (- x 0.5))
  (load push (force -103.8N 0 0) (region (- 98.5 x))) (res 0.5))"#,
        ),
        (
            "boss_x",
            0.4998,
            r#"(def part (min (move (box 70 1.5 1.5) 35.0 0 0) (move (box 6 14 14) 35.0 0 0)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 1000MPa)) (fix base (- x 0.5))
  (load push (force -94.3N 0 0) (region (- 69.5 x))) (res 0.5))"#,
        ),
        (
            "chain2s_x15",
            0.6669,
            r#"(def part (min (min (move (box 30 1.5 1.5) 15.0 0 0) (move (box 10 10 10) 35.0 0 0)) (move (box 30 1.5 1.5) 55.0 0 0)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 1000MPa)) (fix base (- x 0.5))
  (load push (force -75.2N 0 0) (region (- 69.5 x))) (res 0.5))"#,
        ),
    ];
    for (name, factor, src) in cases {
        let res: f64 = src.split("(res ").nth(1).and_then(|r| r.split(')').next()).and_then(|r| r.parse().ok()).expect("a res");
        let s = solve(src, "s", res);
        let b = s.report.buckle.as_ref().unwrap_or_else(|| panic!("{name}: no buckling analysis: {:?}", s.report.warnings));
        assert!(b.failure.is_none() && b.low < 1.0 && (b.factor / factor - 1.0).abs() < 0.1, "{name}: {b:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_box_beams_flange_and_a_strip_column_have_their_buckling_loads_found() {
    // A PETG box beam, 20 x 20 with 1.2 walls, 120 long, cantilevered: its
    // compressed flange, a plate 18.8 wide held on both edges, buckles at
    // 4 pi^2 E / (12 (1 - nu^2)) (t / b)^2 = 32.9 MPa at E=2100 MPa
    // by hand (its webs hold it stiffer). At that reference E its factor
    // is 1.87 at 150 N (33 MPa at its root), and 5 times that at 30 N.
    // Buckling load scales with E: compare at the reference stiffness,
    // using the same current PETG E the study reads. Neither load buckles
    // the current material. The steel strip pushed 2000 N buckles at 0.0206
    // (Euler, fixed-free), and says so.
    let beam = |f: f64| {
        format!(
            "(body beam (shape (subtract (move (box 120 20 20) 60 0 0) (move (box 114 17.6 17.6) 60 0 0))) (attr material \"petg\"))
(study bend (shape (shape-of beam)) (elastic (material petg)) (fix root (- x 1)) (load end (force 0 0 {}N) (region (- 117 x))))",
            -f
        )
    };
    let factor = |s: &sim::Solution| s.report.buckle.as_ref().and_then(|b| b.answer()).map(|a| a.1);
    let (heavy, light) = (factor(&solve(&beam(150.0), "bend", 0.6)), factor(&solve(&beam(30.0), "bend", 0.6)));
    let stiffness_ratio = ink::materials::lookup("petg").unwrap()["E"].as_f64().unwrap() / 2100.0;
    assert!(
        heavy.zip(light).is_some_and(|(h, l)| h > 1.0 && h / stiffness_ratio > 1.5 && h / stiffness_ratio < 2.5 && (l / h / 5.0 - 1.0).abs() < 0.02),
        "{heavy:?} {light:?} (E/reference E: {stiffness_ratio})"
    );
    let strip = "(body s (shape (move (box 100 10 1) 50 0 0)) (attr material \"steel\"))
(study push (shape (shape-of s)) (elastic (material steel)) (fix root (- x 0.5)) (load end (force -2000N 0 0) (region (- 99 x))))";
    let s = solve(strip, "push", 0.25);
    assert!(factor(&s).is_some_and(|f| (f / 0.0206 - 1.0).abs() < 0.1) && s.report.warnings.iter().any(|w| w.starts_with("buckles")), "{:?}", s.report.warnings);
}

/// Each solver end to end gives the same bits with
/// one thread as with 2 and 4 (moved from the kernel's crates/ink/tests/facade/threads.rs):
/// an answer is a function of the question.
#[test]
fn every_solve_is_the_same_at_any_thread_count() {
    type Run = Box<dyn Fn() -> String + Send + Sync>;
    let mut cases: Vec<(&str, Vec<usize>, usize, Run)> = Vec::new();
    // Each solver end to end, called directly (a cached solve would pass
    // without running): the complete serialized solution, all
    // displacements, stresses, reactions, modes and history, not just a
    // headline value. A balanced pull on a pivot leaves a free turn held at
    // rest: support energy rounding feeds the null-space eigenvectors and
    // stabilization. The modal solve (15,000 dofs) and transient heat (18,000)
    // take their CG and eigen reductions over several 8192-value chunks.
    for (name, src, study, res) in [
        (
            "an elastic mechanism",
            "(body l (shape (subtract (move (box 22 4 2) 9 0 0) (cylinder-z 1 6))))
(study s (shape (shape-of l)) (elastic (E 70000MPa) (nu 0.3))
  (fix pin (- (len2 x y) 1.8) (pivot (at 0 0 0) (axis 0 0 1)))
  (load pull (force 10N 0 0) (region (- 18 x))))",
            "s",
            0.5,
        ),
        (
            "a modal solve",
            "(body b (shape (box 20 4 3)))
(study vib (shape (shape-of b)) (elastic (E 200000MPa) (nu 0.3) (density 7850kg/m3))
  (fix wall (+ x 9.5)) (modes 2))",
            "vib",
            0.4,
        ),
        (
            "transient heat",
            "(body b (shape (box 20 4 3)))
(study warm (shape (shape-of b)) (heat (k 20W/m/K) (density 1000kg/m3) (capacity 1000J/kg/K))
  (temperature hot (+ x 9.5) 100degC)
  (transient (duration 10s) (steps 5) (initial 20degC)))",
            "warm",
            0.25,
        ),
    ] {
        cases.push((
            name,
            vec![2, 4],
            1,
            Box::new(move || {
                let mut d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
                let ctx = d.ctx(&[], 0.0).unwrap();
                let (setup, _, _) = ink::sim::study_setup(&mut d, &ctx, study, Some(res), 1e-9).unwrap();
                let mut s = ink::sim::solve(&setup, &ctx).unwrap();
                assert!(s.report.converged && s.report.failure.is_none(), "{:?}", s.report);
                s.report.setup_ms = 0.0;
                s.report.solve_ms = 0.0;
                if let Some(b) = &mut s.report.buckle {
                    b.ms = 0.0;
                }
                serde_json::to_string(&s).unwrap()
            }),
        ));
    }

    std::thread::scope(|scope| {
        let runs: Vec<_> = cases
            .iter()
            .map(|(name, threads, times, run)| {
                (
                    name,
                    scope.spawn(move || {
                        let at = |n| rayon::ThreadPoolBuilder::new().num_threads(n).build().unwrap().install(run);
                        let one = at(1);
                        for &n in threads {
                            for _ in 0..*times {
                                assert!(at(n) == one, "{name}: the answer differs with {n} threads");
                            }
                        }
                    }),
                )
            })
            .collect();
        let failed: Vec<_> = runs.into_iter().filter_map(|(name, r)| r.join().is_err().then_some(*name)).collect();
        assert!(failed.is_empty(), "differs with the thread count: {failed:?}");
    });
}

fn doc(src: &str) -> ink::lang::Doc {
    ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"))
}

/// (Moved from the kernel's crates/ink/tests/facade/kernel.rs, with its study's claim.)
#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn numerical_claims_cover_only_the_scope_their_methods_evaluate() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    use ink::view::ModelView;
    let src = "(def h (- 2 (/ t 1s))) (def plate (box 4 4 h))
        (body p (shape plate))
        (study bend (shape plate) (elastic (E 1000MPa) (nu 0))
          (fix wall (+ x 1.5)) (load tip (force 0 0 -1N) (region (- 1.5 x))) (res 1))
        (require metric (max-displacement bend < 1))
        (require thick (min-thickness plate > 1.5))
        (require fixed (min-thickness (box 4 4 2) > 1.5))";
    let d = doc(src);
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    let opts = ClaimOptions { res: Some(0.5), ..Default::default() };
    let r = run_asserts(&mut v, None, false, opts).unwrap();
    assert!(!r.asserts.is_empty() && r.asserts.iter().all(|a| a.check.status == ink::claims::Status::Holds), "{r:?}");
    // Over a motion they are measured at each step: the plate thins below
    // 1.5 after t = 0.5, which the steps past it see.
    v.set_motion(Some((0.0, 1.0))).unwrap();
    let r = run_asserts(&mut v, None, false, opts).unwrap();
    // Estimates (no tolerance), a failure at a step too: rays guarantee none.
    let a = &r.asserts;
    // The study holds at each moment read, and no range of t bounds what it measures: between the
    // moments it is not settled, saying so (never a yes from the moments).
    // Or a moment is not settled itself: its tip, read on its face, reads
    // 0.007124, 0.007162 and 0.007222 mm at t 0.125, growing.
    assert!(a[0].check.status == Status::Unknown && a[0].check.motion.as_ref().is_some_and(|m| m.between.is_some() || !m.motion.undecided.is_empty()), "{:?}", a[0]);
    assert!(a[1].check.status == Status::Fails && a[1].check.accuracy.is_none(), "{:?}", a[1]);
    assert!(a[2].check.status == Status::Holds && a[2].check.accuracy.is_none(), "{:?}", a[2]);
    assert!(r.asserts.iter().all(|a| a.check.motion.as_ref().is_none_or(|m| !m.bound)));
    // Runtime rays can evaluate at t = 1.
    v.set_motion(None).unwrap();
    let ctx = v.doc.ctx(&[], 1.0).unwrap();
    let mut v = ModelView::new(v.doc, ctx);
    let r = run_asserts(&mut v, None, false, opts).unwrap();
    assert!(r.asserts[1].check.status == Status::Fails && r.asserts[1].check.accuracy.is_none(), "{:?}", r.asserts[1]);
}

/// An extension's run is kept by the scope of the fields it depends on: a
/// cantilever study of a plate thinning with t, run at t = 0 and asked
/// again at another moment (set_moment, or at_point narrowing a span to
/// it), runs there, and the view's own run is back after; so a claim on its
/// deflection over a motion, measured at 4 steps, fails as it does at the
/// last of them (it held, every step reading the first's 0.0069 mm).
#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn an_extension_run_is_kept_by_the_scope_its_inputs_had() {
    use ink::view::ModelView;
    let study = |t: f64| {
        let d = doc("(def plate (box 4 4 (- 2 (/ t 1s))))
            (study bend (shape plate) (elastic (E 1000MPa) (nu 0))
              (fix wall (+ x 1.5)) (load tip (force 0 0 -1N) (region (- 1.5 x))) (res .5))
            (require stiff (max-displacement bend < .02))");
        let ctx = d.ctx(&[], t).unwrap();
        ModelView::new(d, ctx)
    };
    let sag = |v: &mut ModelView| v.outcome(None, "bend", Default::default()).unwrap().metric("max-displacement").unwrap().value.unwrap();
    let (zero, one) = (sag(&mut study(0.0)), sag(&mut study(1.0)));
    let mut v = study(0.0);
    assert_eq!(sag(&mut v), zero);
    v.set_moment(1.0).unwrap();
    assert_eq!(sag(&mut v), one);
    let mut v = study(0.0);
    v.ctx.t_span = Some((0.0, 1.0));
    let at = ink::eval::Ctx { t: 1.0, t_span: None, ..v.ctx.clone() };
    assert_eq!(v.at_point(at, |v| sag(v)), one);
    v.ctx.t_span = None;
    assert_eq!(sag(&mut v), zero);
    let mut v = study(0.0);
    v.set_motion(Some((0.0, 1.0))).unwrap();
    let over = ink::claims::run_asserts(&mut v, None, false, Default::default()).unwrap();
    let last = ink::claims::run_asserts(&mut study(0.875), None, false, Default::default()).unwrap();
    assert!(zero < 0.02 && one > 0.02 && over.asserts[0].check.status == last.asserts[0].check.status, "{zero} {one} {:?}", over.asserts[0].check);
}

#[test]
fn an_empty_clause_is_an_error_at_it() {
    // () as a study's clause: a clause with no name, which every reader splits the same way.
    let src = "(study s ())";
    let e = ink::load_str(&format!("(import \"std\")\n{src}"), "t", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{src} loaded"));
    assert_eq!(e.line, 2, "{src}: {e:?}");
}

// Moved from the kernel's crates/ink/tests/facade/geometry.rs: results as fields, package
// forms and metrics.
const CTX: Ctx = Ctx::new(0.0, Vec::new());

fn field(src: &str, name: &str) -> (ink::lang::Doc, ink::tape::Tape) {
    let d = doc(src);
    let ink::lang::Value::Field(f, _) = d.def(name).unwrap_or_else(|| panic!("no def {name}")).clone() else { panic!("{name} is not a field") };
    let t = ink::tape::Tape::new(&d.g, f);
    (d, t)
}

#[test]
fn study_results_are_fields() {
    // A bar in tension: 10 MPa everywhere away from the ends, and displacement
    // growing linearly from the clamp. Both come back as ordinary fields.
    let src = "(body bar (shape (box 100 10 10)))
(study pull (shape (shape-of bar)) (elastic (E 1000MPa) (nu 0))
  (fix wall (+ x 49.5)) (load pull (force 1000N 0 0) (region (- 49 x))) (res 1))
(def vm (study-field pull von-mises))
(def u (study-field pull displacement))
(def thick (- (+ (abs y) 0) (* 0.1 vm)))";
    let (d, vm) = field(src, "vm");
    let mut v = vec![0.0; 2];
    eval_points(&vm, &CTX, &[[0.0, 0.0, 0.0], [20.0, 2.0, -3.0]], &mut v);
    assert!(v.iter().all(|s| (s - 10.0).abs() < 0.2), "stress {v:?}");
    let Value::Field(u, _) = d.def("u").unwrap().clone() else { panic!() };
    let u = Tape::new(&d.g, u);
    eval_points(&u, &CTX, &[[0.0, 0.0, 0.0], [40.0, 0.0, 0.0]], &mut v);
    assert!((v[0] - 0.4975).abs() < 0.015 && (v[1] - 0.8975).abs() < 0.015, "u {v:?}");
    // Fields built on them stay soundly boundable (the grid pyramid).
    let Value::Field(t, _) = d.def("thick").unwrap().clone() else { panic!() };
    let t = Tape::new(&d.g, t);
    let r = ink::eval::eval_interval(&t, &CTX, [ink::interval::I::new(-10.0, 10.0), ink::interval::I::new(1.0, 2.0), ink::interval::I::new(-1.0, 1.0)], None);
    let pts: Vec<[f64; 3]> = (0..1000).map(|i| [-10.0 + 0.02 * i as f64, 1.0 + (i % 7) as f64 / 6.0, -1.0 + (i % 11) as f64 / 5.0]).collect();
    let mut s = vec![0.0; pts.len()];
    eval_points(&t, &CTX, &pts, &mut s);
    assert!(s.iter().all(|v| r.lo <= *v && *v <= r.hi), "{r:?} does not enclose samples");
}

#[test]
fn extension_forms_go_to_their_packages() {
    // A study is read by its registered package and kept; typed access
    // reads it back.
    use ink::sim::Study;
    let src = "(param w 10)
(body bar (shape (box 100 w 10)))
(study pull (shape (shape-of bar)) (elastic (E 1000MPa) (nu 0))
  (fix wall (+ x 49.5)) (load pull (force 1000N 0 0) (region (- 49 x))))";
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let kept: Vec<(&str, &str)> = d.ext.iter().map(|e| (e.package, e.head)).collect();
    assert_eq!(kept, [("sim", "study")]);
    assert_eq!(d.forms::<Study>().count(), 1);
    assert!(d.form::<Study>("pull").is_some() && d.form::<Study>("bar").is_none());
    // Unknown heads list every package's forms; a result the form does not
    // give is an error before anything is solved.
    let e = ink::load_str("(stduy x)", "test", Path::new(ROOT)).err().unwrap();
    assert!(e.hint.as_deref().is_some_and(|h| h.contains("'study'")), "{e:?}");
    let e = ink::load_str(&format!("(import \"std\")\n{src}\n(def temp (study-field pull temperature))"), "test", Path::new(ROOT)).err().unwrap();
    assert!(words(&e.message).contains(&"temperature") && words(&e.message).contains(&"pull") && e.hint.as_deref().is_some_and(|h| h.contains("von-mises")), "{e:?}");
    // A package form's (res h) is a cell size, read as the language's.
    let bad = "(def b 1) (study s (shape (shape-of bar)) (res -1) (elastic (E 1000MPa) (nu 0)) (fix wall (+ x 49.5)))";
    let e = ink::load_str(&format!("(import \"std\")\n{src}\n{bad}"), "test", Path::new(ROOT)).err().unwrap();
    assert!(["res", "-1"].iter().all(|w| words(&e.message).contains(w)), "{bad}: {e:?}");
    // A document reads the forms of the packages it is loaded with: the
    // core alone (and std, which is plain core) knows none of these.
    static CORE: ink::lang::Registry = ink::lang::Registry::new(&[]);
    let e = Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "test", Path::new(ROOT), &[]).err().unwrap();
    assert!(words(&e.message).contains(&"study"), "{e:?}");
}

/// A study's metric claims (the package metric path itself, with a gauge
/// package: the kernel's crates/ink/tests/facade/geometry.rs).
#[test]
fn a_study_metric_it_does_not_give_or_of_no_study_is_an_error() {
    // A metric the form does not give says what it needs and what the form
    // gives; a metric of a form not yet defined says so.
    let bar = "(import \"std\")\n(body bar (shape (box 100 10 10)))
(study pull (shape (shape-of bar)) (elastic (E 1000MPa) (nu 0)) (fix wall (+ x 49.5)) (load pull (force 1000N 0 0) (region (- 49 x))))";
    let e = ink::load_str(&format!("{bar}\n(require f (min-frequency pull > 5Hz))"), "test", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("study 'pull' does not give 'min-frequency'"), "{e:?}");
    assert!(e.hint.as_deref().is_some_and(|h| h.contains("(modes n)") && h.contains("max-stress")), "{e:?}");
    let e = ink::load_str(&format!("{bar}\n(require stiff (compliance push < 1))"), "test", Path::new(ROOT)).err().unwrap();
    assert!(words(&e.message).contains(&"push"), "{e:?}");
}

/// A static study's linear buckling factor as a claim reads it
/// (`min-buckling-factor`, answered on the finest of its grids).
fn buckling_factor(src: &str, study: &str) -> f64 {
    let c = checks(&format!("{src}\n(require b (min-buckling-factor {study} > 1e6))")).remove(0).1;
    measured(&c).unwrap_or_else(|e| panic!("no factor: {e}: {c:?}"))
}

/// A steel column 100 long, 4 x 4, along x from x = 0, held at its ends
/// by `ends` and pushed along -x by `force` N over its last millimetre:
/// Euler's load pi^2 E I / (K L)^2, E I = 200000 x 4^4 / 12.
fn column(ends: &str, force: f64) -> String {
    format!(
        "(body c (shape (move (box 100 4 4) 50 0 0)))
(study push (shape (shape-of c)) (elastic (E 200000MPa) (nu 0.3))
  {ends} (load end (force {}N 0 0) (region (- 99 x))) (res 1))",
        -force
    )
}

const EULER: f64 = std::f64::consts::PI * std::f64::consts::PI * 200000.0 * 256.0 / 12.0;

#[test]
fn a_pinned_column_buckles_at_eulers_load() {
    // Pinned at both ends (each end's face a plate free to turn about y;
    // the far one also free to slide along the column): K = 1, 4211 N.
    // Pushed with 1000 N it buckles at 4.21 times that. A claim of a
    // margin of 3 holds, on cells of 1 mm checked on cells of 2 and 4.
    let src = column("(fix base (- x 0.01) (pivot (at 0 0 0) (axis 0 1 0))) (fix top (- 99.99 x) (slide 1 0 0) (pivot (at 100 0 0) (axis 0 1 0)))", 1000.0);
    let f = buckling_factor(&src, "push");
    let want = EULER / 100f64.powi(2) / 1000.0;
    eprintln!("pinned-pinned: {f} against {want}");
    assert!((f / want - 1.0).abs() < 0.05, "{f} against {want}");
    let (_, c) = &checks(&format!("{src}\n(require margin (min-buckling-factor push > 3))"))[0];
    assert_eq!(c.status, ink::claims::Status::Holds, "{c:?}");
    // Its buckling mode as a field, (mode 1) of the static study: a half
    // sine bowing across the pin's axis, 1 mm at mid-length, its factor
    // Euler's within 2% (4.163 on these 1 mm cells, 4.197 and 4.205 on 0.5
    // and 0.25 mm: converging from below at second order to 4.208; 4.20
    // read before the pins were held by Nitsche's method was their end
    // cells held whole, a shorter column).
    let d = ink::load_str(&format!("(import \"std\")\n{src}\n(def bow (study-field push displacement-z (mode 1)))"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ink::lang::Value::Field(bow, _) = d.def("bow").unwrap().clone() else { panic!() };
    let t = ink::tape::Tape::new(&d.g, bow);
    let what = t.grids.first().and_then(|g| g.source()).map(|s| s.what.clone()).unwrap_or_default();
    let factor: f64 = what.split("factor ").nth(1).and_then(|r| r.split(',').next()).and_then(|v| v.trim().parse().ok()).unwrap_or(f64::NAN);
    assert!(what.contains("buckling mode 1: factor") && (factor / want - 1.0).abs() < 0.02 && what.ends_with("1 mm peak"), "{what}");
    let ctx = d.ctx(&[], 0.0).unwrap();
    let at = |x: f64| ink::eval::eval_point(&t, &ctx, [x, 0.0, 0.0]);
    let mid = at(50.0);
    assert!((mid.abs() - 1.0).abs() < 0.05 && (at(20.0) / mid - (0.2 * std::f64::consts::PI).sin()).abs() < 0.05, "{mid}, {}", at(20.0));
    let e =
        format!("{:?}", ink::load_str(&format!("(import \"std\")\n{src}\n(def m2 (study-field push displacement-z (mode 2)))"), "test", Path::new(ROOT)).err().expect("no mode 2"));
    assert!(e.contains("first buckling mode only"), "{e}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_column_with_split_cells_buckles_as_its_plain_twin() {
    // The pinned column, and it slit lengthwise near its base (x 2 to 12,
    // 0.4 mm across y, from the top down to a 1 mm ligament whose top lies
    // 0.016 mm above a cell's plane, as the cantilever's slit), its cells
    // split there. Its slit plane carries no shear in this buckle, and it
    // takes 7.5% of the section where the bow is slight: the factor barely
    // moves. The split cells' geometric stiffness was left out (only the
    // grid's were applied): 4.39 against the plain column's 4.21.
    let plain = column("(fix base (- x 0.01) (pivot (at 0 0 0) (axis 0 1 0))) (fix top (- 99.99 x) (slide 1 0 0) (pivot (at 100 0 0) (axis 0 1 0)))", 1000.0);
    let slit = plain.replace("(shape (move (box 100 4 4) 50 0 0))", "(shape (subtract (move (box 100 4 4) 50 0 0) (move (box 10 0.4 4) 7 0.12 1.016)))");
    assert!(solve(&slit, "push", 1.0).report.cells_refined > 0);
    let (f, f0) = (buckling_factor(&slit, "push"), buckling_factor(&plain, "push"));
    assert!(f <= 1.002 * f0 && f > 0.97 * f0, "slit {f}, plain {f0}");
}

#[test]
fn a_cantilevered_column_buckles_at_eulers_load() {
    // Built in at x = 0, free at the far end: K = 2, 1053 N. Pushed with
    // 300 N it buckles at 3.51 times that: both ways equally (the section
    // is square), its first two factors alike.
    let src = column("(fix base (- x 0.01))", 300.0);
    let f = buckling_factor(&src, "push");
    let want = EULER / 200f64.powi(2) / 300.0;
    eprintln!("fixed-free: {f} against {want}");
    assert!((f / want - 1.0).abs() < 0.05, "{f} against {want}");
}

#[test]
fn a_column_built_in_at_both_ends_buckles_at_eulers_load() {
    // Built in at x = 0, its far end's face a plate free only to slide
    // along it: K = 1/2, 16845 N (shear lowers it 1.6% at this
    // stockiness). Pushed with 3000 N it buckles at 5.6 times that.
    let src = column("(fix base (- x 0.01)) (fix top (- 99.99 x) (slide 1 0 0))", 3000.0);
    let f = buckling_factor(&src, "push");
    let want = EULER / 50f64.powi(2) / 3000.0;
    eprintln!("fixed-fixed: {f} against {want}");
    assert!((f / want - 1.0).abs() < 0.05, "{f} against {want}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_plate_pressed_in_its_plane_buckles_at_its_critical_stress() {
    // A steel plate 40 x 40 x 1, its four edges held from moving out of
    // its plane (simply supported), 200 MPa of pressure on two opposite
    // edges: it buckles at k pi^2 D / (b^2 t), k = 4 for a square, D = E
    // t^3 / (12 (1 - nu^2)): 451.9 MPa, 2.26 times the load; 2.2537.
    // Each edge's face is held, not a band of material under it: a band
    // held from moving out of the plane, however thin, holds the edge's
    // turn by its shear, G times its width (0.01 mm reads 2.95, 0.005
    // 2.63, 0.02 3.46, on 0.5 mm cells as on 0.25). It read 2.26 on 0.5
    // mm cells only while a held sliver under 3% of a cell lost its free
    // axes' stiffness, and 2.94 on 0.25 mm cells then too.
    // (Its guides 0.02 mm thick read 2.197, 2.8% under: the plate's two
    // halves shared their nodes across them, held as one field.)
    let free = "(body p (shape (move (box 40 40 1) 20 20 0.5)))
(study press (shape (shape-of p)) (elastic (E 200000MPa) (nu 0.3))
  (fix left-edge (face p (normal -1 0 0)) (dofs z)) (fix right-edge (face p (normal 1 0 0)) (dofs z))
  (fix near-edge (face p (normal 0 -1 0)) (dofs z)) (fix far-edge (face p (normal 0 1 0)) (dofs z))
  (load left (pressure 200MPa) (region (- x 0.01))) (load right (pressure 200MPa) (region (- 39.99 x))) (res 0.5))";
    // Balanced dead loads leave the plate free to spin in its plane.
    // Holding that rotation by a numerical gauge would hide its instability.
    let c = checks(&format!("{free}\n(require b (min-buckling-factor press > 1))")).remove(0).1;
    assert_eq!(c.status, ink::claims::Status::Unknown, "{c:?}");
    assert!(measured(&c).unwrap_err().contains("free rigid rotation"), "{c:?}");
    // Guides on both centre lines, two cells thick (1 mm, the plate's 0.5
    // mm cells: a thinner one, its free material on both sides sharing
    // nodes, is refused), hold that rigid motion and admit the classical
    // square plate's first bending mode and Poisson strain.
    let guided = free.replace("(load left", "(fix guide-x (- (abs (- x 20)) 0.5) (dofs x) (support)) (fix guide-y (- (abs (- y 20)) 0.5) (dofs y) (support)) (load left");
    let f = buckling_factor(&guided, "press");
    let d = 200000.0 / (12.0 * (1.0 - 0.09));
    let want = 4.0 * std::f64::consts::PI.powi(2) * d / 1600.0 / 200.0;
    eprintln!("plate: {f} against {want}");
    assert!((f / want - 1.0).abs() < 0.05, "{f} against {want}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_thin_tube_pressed_along_its_axis_buckles_near_the_classical_stress() {
    // A steel tube of mean radius 20 and wall 1, 60 long, built in at
    // one end and its other end's face free only to slide along it, held
    // rigidly (an end plate: a slide alone holds the face's mean motion,
    // its rim free to go oval, and reads 1.29), pushed with 2000 MPa over
    // its section. Thin shell theory puts its
    // wall's buckle at E t / (R sqrt(3 (1 - nu^2))) = 6052 MPa, 3.03
    // times the load; it leaves out the wall's shear (2.4% lower here,
    // a half-wave of 7.5 mm across a wall of 1) and the clamped ends'
    // hold on the wall's bulge. On cells of 0.5 mm it reads 2.84, 6%
    // under.
    let area = 2.0 * std::f64::consts::PI * 20.0;
    let src = format!(
        "(body t (shape (move (subtract (cylinder-z 20.5 60) (cylinder-z 19.5 61)) 0 0 30)))
(study push (shape (shape-of t)) (elastic (E 200000MPa) (nu 0.3))
  (fix base (- z 0.01)) (fix top (- 59.99 z) (slide 0 0 1) (rigid)) (load end (force 0 0 {}N) (region (- 59 z))) (res 0.5))",
        -2000.0 * area
    );
    let f = buckling_factor(&src, "push");
    let want = 200000.0 / (20.0 * (3.0 * (1.0 - 0.09f64)).sqrt()) / 2000.0;
    eprintln!("tube: {f} against {want}");
    assert!((f / want - 1.0).abs() < 0.07, "{f} against {want}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_buckling_margin_holds_only_where_its_factor_is_resolved() {
    // The cantilevered column at 0.03 N: its factor is 10,000 times the
    // 3.55 it has at 300 N (stress, and K_G with it, scale with the load),
    // so 35,500: a margin of 36,000 fails, however large the factor, on
    // every grid. A factor resolved only to a share of a fixed scale read
    // 37,193 here, and held.
    let src = column("(fix base (- x 0.01))", 0.03).replace("(res 1)", "(res 0.25)");
    let (_, c) = &checks(&format!("{src}\n(require b (min-buckling-factor push > 36000))"))[0];
    assert_eq!(c.status, ink::claims::Status::Fails, "{c:?}");
    assert!(measured(c).is_ok_and(|f| (f / 35500.0 - 1.0).abs() < 0.02), "{c:?}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_bar_in_tension_has_its_buckling_factor_bounded_by_what_was_checked() {
    // A steel bar pulled at 100 MPa: its clamp's hold on its Poisson
    // contraction compresses its corners, which buckle at a factor of at
    // most 16,324 on 1 mm cells (a Rayleigh quotient, -6.126e-5, a longer
    // search finds). The eigensolve's 20 steps find no direction its loads
    // compress, and exclude every factor below about 2800 by a residual
    // that takes no mode as missed. What is checked is that 100 times the
    // loads leave the part stable (a run that takes none as found), so
    // the factor is "at least" 100, a band with no top. `> 3` holds, `< 10` fails, and
    // `> 1e6` is unknown (it held as "no buckling mode", and failed when
    // three infinite grids differenced to NaN); its strength claims are
    // decided by strength (a safety factor of 2.4, far under the bound).
    let src = "(body bar (shape (box 100 10 10)))
(study pull (shape (shape-of bar)) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa))
  (fix wall (+ x 49.5)) (load pull (force 10000N 0 0) (region (- 49 x))) (res 1))
(require margin (min-buckling-factor pull > 3))
(require huge (min-buckling-factor pull > 1e6))
(require small (min-buckling-factor pull < 10))
(require strong (min-safety-factor pull > 1.5))
(require weak (min-safety-factor pull > 3))";
    use ink::claims::Status::{Fails, Holds, Unknown};
    let c = checks(src);
    for ((name, c), want) in c.iter().zip([Holds, Unknown, Fails, Holds, Fails]) {
        assert_eq!(c.status, want, "{name}: {c:?}");
    }
    assert!(undecided(&c[1].1).contains("eigensolve search"), "{:?}", c[1].1);
    for (name, c) in &c[..3] {
        assert!(measured(c).is_ok_and(|v| v == f64::INFINITY), "{name}: {c:?}");
        let j = ink::lang::answer::check_json(c, &[]);
        assert!(j["value"].as_str().is_some_and(|v| v.starts_with("at least ")), "{name}: {j}");
        assert!(j["band"]["low"].as_f64() == Some(100.0), "{name}: {j}");
        assert!(j.get("at").is_none_or(|a| a.is_null()), "{name}: {j}");
    }
    assert!(measured(&c[3].1).is_ok_and(|v| (2.0..3.0).contains(&v)), "{:?}", c[3].1);
    let mut d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let (setup, _, _) = ink::sim::study_setup(&mut d, &ctx, "pull", None, 1e-9).unwrap();
    let r = ink::sim::solve(&setup, &ctx).unwrap().report;
    assert!(r.buckle.as_ref().is_some_and(|b| b.factor == f64::INFINITY && b.low == 100.0), "{:?}", r.buckle);
    assert!(r.warnings.iter().all(|w| !w.contains("buckl")), "{:?}", r.warnings);
}

/// A soft bar pulled, E 1 MPa: its search for a compressed mode is long
/// and inconclusive, as every factor is small.
fn soft_pull(clamp: f64, nu: f64) -> String {
    format!(
        "(def bar (box 12 4 4))
(study pull (shape bar) (elastic (E 1MPa) (nu {nu})) (fix wall (+ x {clamp}))
  (load tip (pressure -0.625MPa) (region (face bar (normal 1 0 0)))) (res 1))
(require short (max-displacement pull < 0.000001))
(require stable (min-buckling-factor pull > 1e6))"
    )
}

#[test]
fn a_part_with_no_compressed_cell_has_no_buckle_and_needs_no_search() {
    // A clamp one cell thick on the finest grid: no cell's stress
    // compressive past roundoff, on every grid. Nothing buckles below
    // what that roundoff leaves open (millions of times the load), so no
    // search runs; the search, run, stopped unconverged on 0.5 mm cells
    // and left its claims unknown.
    let src = soft_pull(5.5, 0.0);
    for res in [0.5, 1.0, 2.0] {
        let b = solve(&src, "pull", res).report.buckle.expect("a buckling analysis");
        assert!(b.failure.is_none() && b.iterations == 0 && b.factor == f64::INFINITY && b.low > 1e6, "{res}: {b:?}");
    }
    // So in a clamp two cells thick: its held material is no free
    // material's (a penalty over it compressed its outer cells, 1/600 of
    // the pull, and the search ran and did not resolve).
    let b = solve(&soft_pull(5.0, 0.0), "pull", 0.5).report.buckle.expect("a buckling analysis");
    assert!(b.failure.is_none() && b.iterations == 0 && b.factor == f64::INFINITY && b.low > 1e6, "thick clamp: {b:?}");
    use ink::claims::Status::Holds;
    let c = checks(&src);
    assert_eq!(c[0].1.status, ink::claims::Status::Fails, "{:?}", c[0].1);
    // 0.625 MPa over the 11.5 mm out of the clamp.
    assert!(measured(&c[0].1).is_ok_and(|v| (v - 7.1875).abs() < 0.01), "{:?}", c[0].1);
    assert_eq!(c[1].1.status, Holds, "{:?}", c[1].1);
}

#[test]
fn a_small_compressed_region_leaves_a_short_search_unresolved() {
    // The same bar, its Poisson's ratio 0.3, in a clamp two cells thick:
    // where the clamp holds the bar's contraction, by its edges, it is
    // compressed across, a small region. The search runs, does not resolve
    // in its 100 steps, and the buckle stays unresolved: the displacement
    // that fails by six orders of magnitude is unknown. (Of Poisson's ratio
    // 0, as its sibling above, the bar is compressed nowhere: the penalty
    // that held clamps once compressed their cells, 1/600 of the pull.)
    let src = soft_pull(5.0, 0.3);
    let b = solve(&src, "pull", 0.5).report.buckle.expect("a buckling analysis");
    assert!(b.iterations > 0 && b.answer().is_none() && b.failure.as_deref().is_some_and(|f| f.contains("did not resolve")), "{b:?}");
    let c = checks(&src);
    for (name, c) in &c {
        assert_eq!(c.status, ink::claims::Status::Unknown, "{name}: {c:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_weak_strip_of_yield_is_read_or_left_undecided_wherever_it_lies() {
    // A 40 x 10 x 10 bar pulled at 10 MPa, its yield 250 MPa but for a
    // strip 0.2 mm wide across it at 10: a safety factor of about 1. The
    // yield was read once per cell, at its peak stress: with the strip
    // between a cell's peak and the rest, `> 2` held at 24. Each point read
    // takes its own yield, and a strip between them all leaves the claim
    // unknown, saying so; `> 0.5`, true, never fails.
    let strip = |x: f64, w: f64, claim: &str| {
        format!(
            "(def part (box 40 10 10))\n\
             (study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield (* 1MPa (select (max (- {x} x) (- x {})) 10 250))))\n\
             (fix root (+ x 19)) (load pull (force 1000N 0 0) (region (- 19 x))) (res 2))\n{claim}",
            x + w
        )
    };
    let bar = |x: f64, claim: &str| strip(x, 0.2, claim);
    // Narrower strips, 0.05 and 0.02 mm, fall between the reads.
    for (x, w) in [(1.17, 0.05), (1.1, 0.02)] {
        let (name, c) = &checks(&strip(x, w, "(require safe (min-safety-factor s > 2))"))[0];
        assert!(c.status == ink::claims::Status::Unknown && format!("{c:?}").contains("yield-between-reads"), "strip {w} at {x}: {name}: {c:?}");
    }
    for x in [0.3, 0.7, 1.1, 1.45] {
        let (name, c) = &checks(&bar(x, "(require safe (min-safety-factor s > 2))"))[0];
        assert_ne!(c.status, ink::claims::Status::Holds, "strip at {x}: {name}: {c:?}");
        if c.status == ink::claims::Status::Unknown {
            assert!(format!("{c:?}").contains("yield-between-reads"), "strip at {x}: {c:?}");
        } else {
            assert!(measured(c).is_ok_and(|v| (0.9..1.1).contains(&v)), "strip at {x}: {c:?}");
        }
        // (True, but the grids read the strip on some cells and not others,
        // so may not agree: never failed.)
        let (name, c) = &checks(&bar(x, "(require safe (min-safety-factor s > 0.5))"))[0];
        assert_ne!(c.status, ink::claims::Status::Fails, "strip at {x}: {name}: {c:?}");
    }
    // One yield throughout: nothing between the reads.
    let (name, c) = &checks(&bar(0.3, "(require safe (min-safety-factor s > 2))").replace("(select (max (- 0.3 x) (- x 0.5)) 10 250)", "250"))[0];
    assert_eq!(c.status, ink::claims::Status::Holds, "{name}: {c:?}");
}

/// A claim's value counting what it leaves out under `key`, if it does.
fn excluded_value(c: &ink::claims::Check, key: &str) -> Option<f64> {
    match &c.evidence {
        ink::claims::Evidence::Metric { measured, .. } => measured.excluded.iter().find(|e| e.key == key).map(|e| e.value),
        e => panic!("not a metric: {e:?}"),
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_weak_strip_in_bending_or_by_a_fix_is_bounded_by_the_stress_over_its_cell() {
    // A 100 x 10 x 10 cantilever, yield 1000 MPa but for a strip 0.05 mm
    // thick at z 0.9 to 0.95, inside the cells z 0 to 1, read at their
    // centres. Bent by a couple, its neutral axis on those centres (they
    // read about 0), the strip carries 4.86 MPa against 4 MPa: a factor of
    // 0.82; bent by a tip load, 10.25 MPa against 10: 0.975. A cell's bound
    // took the stress read in it, at its centre: `> 1.5` held at 11 and 18.
    // It takes the stress's bound over the cell (its corners).
    let beam = |weak: f64, x0: &str, z0: f64, loads: &str| {
        format!(
            "(def part (move (box 100 10 10) 0 0 {z0}))\n\
             (study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield (* 1MPa (select (max (- -39.95 x) {x0} (- (abs y) 3.95) (- 0.9 z) (- z 0.95)) {weak} 1000))))\n\
             (fix root (+ x 49)) {loads} (res 1))\n(require safe (min-safety-factor s > 1.5))\n(require unsafe (min-safety-factor s < 1.5))"
        )
    };
    let couple = beam(4.0, "(- x 39.95)", 0.5, "(load top (force 1000N 0 0) (region (max (- 49 x) (- 4.5 z)))) (load bottom (force -1000N 0 0) (region (max (- 49 x) (+ z 3.5))))");
    let tip = beam(10.0, "-1", 0.0, "(load tip (force 0 0 -100N) (region (- 49 x)))");
    // A 0.02 mm strip 1.1 mm from the fix, in its support margin, of the
    // bar pulled at 10 MPa: a factor of about 1.
    let margin = "(def part (box 40 10 10))\n\
         (study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield (* 1MPa (select (max (- -17.9 x) (- x -17.88)) 10 250))))\n\
         (fix root (+ x 19)) (load pull (force 1000N 0 0) (region (- 19 x))) (res 1))\n\
         (require safe (min-safety-factor s > 2))\n(require unsafe (min-safety-factor s < 2))"
        .to_string();
    for (case, src, truth, key) in [("couple", couple, 0.82, "yield-between-reads"), ("tip", tip, 0.975, "yield-between-reads"), ("margin", margin, 1.0, "support-margin")] {
        for (name, c) in checks(&src) {
            let want = if name == "safe" { ink::claims::Status::Holds } else { ink::claims::Status::Fails };
            assert_ne!(c.status, want, "{case}: {name}: {c:?}");
            let b = excluded_value(&c, key).unwrap_or_else(|| measured(&c).unwrap());
            assert!(b > 0.0 && b <= truth * 1.02, "{case}: {name}: bound {b} against {truth}: {c:?}");
        }
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_yield_gradients_bound_closes_on_its_factor_as_the_cells_shrink() {
    // The bar pulled at 10 MPa, its yield 250 - 2 x: no strip. No cell is
    // of one material, so none gives a value, and every cell is bounded:
    // the least yield over it (210 at x 20, a factor of about 21) over its
    // stress's bound. The bound reaches a cell, and settles as the cells
    // shrink.
    let src = "(def part (box 40 10 10))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield (* 1MPa (- 250 (* 2 x))))) (fix root (+ x 19)) (load pull (force 1000N 0 0) (region (- 19 x))))";
    let yield_at = |p: &[[f64; 3]]| p.iter().map(|p| 250.0 - 2.0 * p[0]).collect::<Vec<_>>();
    let mixed = |b: &[([f64; 3], [f64; 3])]| b.iter().map(|b| Some(sim::MixedCell { least: 250.0 - 2.0 * b.1[0], stiffest: 200000.0 / 1.3 })).collect::<Vec<_>>();
    let bounds: Vec<f64> = [2.0, 1.0, 0.5]
        .into_iter()
        .map(|res| {
            let s = solve(src, "s", res);
            assert!(s.max_cell_ratio(&yield_at, &mixed).is_err(), "res {res}: a value from cells not of one material");
            1.0 / s.ratio_between_reads(false, &mixed).unwrap().value
        })
        .collect();
    assert!(bounds.iter().all(|b| (b / 21.0 - 1.0).abs() < 0.05) && (bounds[2] - bounds[1]).abs() <= (bounds[1] - bounds[0]).abs() + 1e-9, "{bounds:?}");
    // A bound is the least the factor can be: a claim it holds holds, with
    // no value; a false one is unknown, never failed. (The fix's edge, its
    // stress singular, is the support margin's: left out here.) It holds
    // past what a stress's read may miss: `> 19`, not `> 20`.
    let claims = format!(
        "{}\n(require away (min-safety-factor s > 19 (excluding support-margin)))\n(require high (min-safety-factor s > 22 (excluding support-margin)))",
        src.replace("(region (- 19 x))))", "(region (- 19 x))) (res 1))")
    );
    for (name, c) in checks(&claims) {
        let want = if name == "away" { ink::claims::Status::Holds } else { ink::claims::Status::Unknown };
        assert_eq!(c.status, want, "{name}: {c:?}");
    }
}

#[test]
fn a_soft_weak_skin_is_not_read_with_its_cells_stiffness() {
    // A beam bent to 54 MPa at its faces (yield 250: a factor of about
    // 4.6; 2.7 read on 0.5 mm cells), a 0.1 mm skin on a side face of E 200
    // and yield 25. The skin, sharing the bulk's strain, carries a
    // thousandth of its stress: a factor of about 500. Its cells take the
    // bulk's stiffness, and a read in the skin divided the bulk's stress by
    // the skin's yield: `> 1` failed at 0.48. A cell not of one material
    // is bounded only.
    let src = "(def part (box 100 10 10))
(def skin (max (- 4.9 y) (- y 5.5) (- 4 z) (- z 4.8) (- (abs x) 40)))
(study s (shape part) (elastic (E (* 1MPa (select skin 200 200000))) (nu 0.3) (yield (* 1MPa (select skin 25 250)))) (fix root (+ x 49))
 (load top (force 1000N 0 0) (region (max (- 49 x) (- 4 z)))) (load bottom (force -1000N 0 0) (region (max (- 49 x) (+ z 4)))) (res 1))
(require safe (min-safety-factor s > 1))
(require unsafe (min-safety-factor s < 1))";
    for (name, c) in checks(src) {
        let wrong = if name == "safe" { ink::claims::Status::Fails } else { ink::claims::Status::Holds };
        assert_ne!(c.status, wrong, "{name}: {c:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_cells_yield_bound_counts_where_it_varies_with_its_grid_error_and_never_unmakes_a_fail() {
    let pull = |yield_: &str, claims: &str| {
        format!(
            "(def part (box 40 10 10))\n(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield (* 1MPa {yield_})))\n\
             (fix root (+ x 19)) (load pull (force 1000N 0 0) (region (- 19 x))) (res 1))\n{claims}"
        )
    };
    // A block of yield 10 (whole cells, read) in a yield graded 250 - 2 x:
    // the factor about 1. Counting the graded cells' bound (21) alone left
    // `> 2` unknown, which fails on the value: it counts the least of value
    // and bound.
    let block = pull("(select (max (- 5 x) (- x 10)) 10 (- 250 (* 2 x)))", "(require safe (min-safety-factor s > 2))\n(require unsafe (min-safety-factor s < 2))");
    for (name, c) in checks(&block) {
        let want = if name == "safe" { ink::claims::Status::Fails } else { ink::claims::Status::Holds };
        assert_eq!(c.status, want, "{name}: {c:?}");
    }
    // A 4 MPa core z 0.05 to 0.95 over the centres of a couple-bent beam's
    // neutral-axis cells: read there, at about 0, `> 1.5` held at 11; its
    // faces carry 4.86 MPa, a factor of 0.82. A cell whose yield varies is
    // bounded, a read in it or not.
    let core = "(def part (move (box 100 10 10) 0 0 0.5))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield (* 1MPa (select (max (- -39.95 x) (- x 39.95) (- (abs y) 3.95) (- 0.05 z) (- z 0.95)) 4 1000))))
 (fix root (+ x 49)) (load top (force 1000N 0 0) (region (max (- 49 x) (- 4.5 z)))) (load bottom (force -1000N 0 0) (region (max (- 49 x) (+ z 3.5)))) (res 1))
(require safe (min-safety-factor s > 1.5))";
    // A strip of yield 60 crossing a hole's edge in a pulled plate: its
    // bound on 1 mm cells, 1.62, is above the true 1.41 (42.6 MPa on 0.25
    // and 0.125 mm cells), and `> 1.5` held; with its grid error it can't.
    let hole = "(def part (max (box 40 20 4) (- 3 (sqrt (+ (* x x) (* y y))))))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield (* 1MPa (select (max (- 0.10 x) (- x 0.12)) 60 250)))) (fix root (+ x 19)) (load pull (force 1000N 0 0) (region (- 19 x))) (res 1))
(require safe (min-safety-factor s > 1.5))";
    for (case, src) in [("core", core), ("hole", hole)] {
        let (name, c) = &checks(src)[0];
        assert_ne!(c.status, ink::claims::Status::Holds, "{case}: {name}: {c:?}");
    }
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_yield_whose_bound_reaches_zero_leaves_its_cells_unknown_not_unread() {
    // The bar pulled at 10 MPa, its yield a smoothstep dip to 10 MPa (250 -
    // 240 (3 s^2 - 2 s^3)): a factor of 1 there. Over a cell near the dip
    // the blend's interval overshoots 1, so the yield's least is below 0:
    // divided by, the cell's bound was negative, no maximum picked it, and
    // with no value either the dip was never read: `> 3` held at 4.9 on a
    // 50 MPa end, `< 2` failed at 24.
    let bar = |yield_: &str, claims: &str| {
        format!(
            "(def part (box 40 10 10))\n(def sx (max 0 (- 1 (/ (abs (- x 0.5)) 4))))\n(def blend (- (* 3 sx sx) (* 2 sx sx sx)))\n\
             (study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield (* 1MPa {yield_}))) (fix root (+ x 19)) (load pull (force 1000N 0 0) (region (- 19 x))) (res 1))\n{claims}"
        )
    };
    let dip = "(- 250 (* 240 blend))";
    for (src, safe, unsafe_) in [
        (bar(&format!("(select (- 10 x) 50 {dip})"), "(require safe (min-safety-factor s > 3))"), "safe", "(require unsafe (min-safety-factor s < 3))"),
        (bar(dip, "(require safe (min-safety-factor s > 2))"), "safe", "(require unsafe (min-safety-factor s < 2))"),
    ] {
        let src = format!("{src}\n{unsafe_}");
        for (name, c) in checks(&src) {
            let wrong = if name == safe { ink::claims::Status::Holds } else { ink::claims::Status::Fails };
            assert_ne!(c.status, wrong, "{name}: {c:?}");
        }
    }
}

/// A 20 x 4 x 4 steel bar slit from the top at `x` (0.4 mm wide, a round
/// bottom of radius `r` where given) to a ligament `t` thick at its foot,
/// moved `dz` along z, 1e-5 N at its tip, on cells of `res`.
fn slit_bar(x: f64, t: f64, dz: f64, r: Option<f64>, res: f64, claim: &str) -> String {
    let cut = match r {
        None => format!("(move (box 0.4 8 {}) {x} 0 {})", 4.0 - t, t / 2.0),
        Some(r) => format!("(min (move (box 0.4 8 4) {x} 0 {}) (- (len2 (- x {x}) (- z {})) {r}))", -2.0 + t + r + 2.0, -2.0 + t + r),
    };
    format!(
        "(body bar (shape (move (subtract (move (box 20 4 4) 10 0 0) {cut}) 0 0 {dz})) (attr material \"steel\"))
(study s (shape (shape-of bar)) (elastic (material steel)) (fix root (- x 1)) (load tip (force 0 0 -0.00001N) (region (- 19 x))) (res {res}))
{claim}"
    )
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_ligament_thinner_than_a_cell_is_split_or_said() {
    // A ligament under a slit, 0.1 to 0.2 mm, on 0.5 mm cells: a sample
    // was thick where its material spanned a quarter of a cell (counted a
    // step long), so a ligament of a fifth to two fifths of a cell joined
    // the slit's sides through thick samples: nothing split or said, the
    // cells read it 7 to 68 times too stiff, and that blind grid, the
    // coarsest of a claim's three, let claims hold short of the bar's
    // bending (0.125 mm: 6.0e-6 + 1.3e-6 under 7.5e-6, where 0.05 mm cells
    // read 8.3e-6). Material thinner than two cells is split until it is two
    // cells thick, and left thinner at the floor decides no claim.
    let claims = [
        (5.1, 0.125, 0.03, None, "(require stiff (max-displacement s < 7.5e-6))"),
        (5.1, 0.1, 0.09, None, "(require stiff (max-displacement s < 1.43e-5))"),
        (5.1, 0.05, 0.09, Some(0.2), "(require stiff (max-displacement s < 3.7e-5))"),
    ];
    for (x, t, dz, r, claim) in claims {
        let (name, c) = &checks(&slit_bar(x, t, dz, r, 0.5, claim))[0];
        assert_ne!(c.status, ink::claims::Status::Holds, "ligament {t} at {x}, moved {dz}: {name}: {c:?}");
    }
    // Solves: said, or within a third of cells far finer (0.05 mm, on the
    // slit's planes where they can be).
    let solves =
        [(5.1, 0.09, 0.06, 0.5, 2.0e-5), (5.1, 0.1, 0.09, 0.5, 1.544e-5), (4.9, 0.1, 0.09, 0.5, 1.544e-5), (5.1, 0.2, 0.03, 0.5, 2.26e-6), (5.0, 0.05, 0.0, 0.25, 9.88e-5)];
    for (x, t, dz, res, reference) in solves {
        let s = solve(&slit_bar(x, t, dz, None, res, ""), "s", res);
        let tip = s.max_nodal().unwrap().value;
        let said = s.report.thin_left.is_some();
        assert!(said || (tip / reference - 1.0).abs() < 0.33, "ligament {t} at {x}, moved {dz}, {res} mm: {tip} against {reference}, {:?}", s.report.warnings);
    }
}

/// A 4 x 6 x 6 block held at x < 1 and a 20 x 1 mm fin `t` thick along x
/// from it (x 0 to 20: free from 4), its mid-plane at z `dz`.
fn fin_at(t: f64, dz: f64, study: &str) -> String {
    format!("(def part (min (move (box 4 6 6) 2 0 0) (move (box 20 1 {t}) 10 0 {dz})))\n{study}")
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_thin_fins_displacement_stress_and_temperature_hold_no_claim_short_of_them() {
    // On 0.5 mm cells a fin 0.1 or 0.2 mm thick lay in cells under half
    // full, which the maxima did not read: claims held over it. Under 2.5e-6
    // N its tip moves 1.95e-4 mm; its stress 1 mm from its root is 6 M / (b
    // t^2), 5.4e-3 MPa at 0.2 mm; heated 0.25 mW at its tip it reaches 20.78
    // C. (As a fin 4 wide under 1e-5 N and 1 mW: its width only multiplies
    // the cells across it.)
    // Split until it is two cells thick (0.2 mm), or left thinner at the
    // floor (0.1 mm), none holds.
    let elastic = |t: f64, dz: f64, claim: &str| {
        fin_at(t, dz, &format!("(study s (shape part) (elastic (E 200000MPa) (nu 0)) (fix root (- x 1)) (load tip (force 0 0 -2.5e-06N) (region (- 19 x))) (res 0.5))\n{claim}"))
    };
    for dz in [0.0, 0.1, 0.13] {
        let heat = fin_at(
            0.1,
            dz,
            "(study s (shape part) (heat (k 50W/m/K)) (temperature base (- x 1) 20degC) (source tip (power 0.00025W) (region (- 19 x))) (res 0.5))\n(require cool (max-temperature s < 20.5degC))",
        );
        for src in [elastic(0.1, dz, "(require stiff (max-displacement s < 1.6e-4))"), elastic(0.2, dz, "(require strong (max-stress s < 0.0027187MPa))"), heat] {
            let (name, c) = &checks(&src)[0];
            assert_ne!(c.status, ink::claims::Status::Holds, "dz {dz}: {name}: {c:?}");
        }
    }
}

#[test]
fn a_wall_thinner_than_two_cells_is_split_wherever_it_lies() {
    // A plate 1.5 mm thick on 1 mm cells, its faces half-spaces (a max)
    // at z 0.4 and 1.9. A cell's own tape keeps only the branches that
    // reach into it, and the cells below z 1 measured their samples' chords
    // on one without the top face: the plate ran on past it, not thin, and
    // those 100 cells were left whole under 800 split ones. Chords run on
    // the whole field: every cell holding the plate is split.
    let src = "(body plate (shape (max (- (abs x) 5) (- (abs y) 5) (- 0.4 z) (- z 1.9))) (attr material \"steel\"))
(study s (shape (shape-of plate)) (elastic (material steel)) (fix root (+ x 4)) (load tip (force 0 0 -1N) (region (- x 4))))";
    let r = solve(src, "s", 1.0).report;
    assert_eq!(r.cells_refined, r.cells_cut + r.cells_full, "{} cut and {} full cells, {} of them split", r.cells_cut, r.cells_full, r.cells_refined);
}

#[test]
fn split_cells_read_as_native_cells_their_size() {
    // A ring 10 to 11.5 across and 1.5 long, its bore pressed (Lame: von
    // Mises 15.50 at the bore). Split once, its 1 mm cells are the 0.5 mm
    // grid's cells, and read as they do (16.14 and 16.04). A cell under 3%
    // full by its own samples (a sliver) was dropped before it could be
    // split, though its children hold more of theirs: 908 cut cells against
    // 984, the peak 17.59 (with chords on a cell's own tape, 784 and 17.71).
    let src = "(def r (len2 x y))
(body ring (shape (max (tag bore (- 10 r)) (- r 11.5) (neg z) (- z 1.5))))
(study s (shape (shape-of ring)) (elastic (E 200000MPa) (nu 0))
 (fix sz (- z 0.001) (dofs z)) (fix py (max (- 10.5 x) (- (abs y) 0.5)) (dofs y)) (fix px (max (- 10.5 y) (- (abs x) 0.5)) (dofs x))
 (load inside (pressure 2MPa) (region (face bore))))";
    let (split, native) = (solve(src, "s", 1.0), solve(src, "s", 0.5));
    let (a, b) = (&split.report, &native.report);
    assert_eq!((a.cells_cut, a.cells_full, a.nodes), (b.cells_cut, b.cells_full, b.nodes), "split {a:?}\nnative {b:?}");
    assert!((a.volume - b.volume).abs() < 1e-9 * b.volume, "{} against {}", a.volume, b.volume);
    let (x, y) = (split.max_cell_scalar().unwrap().value, native.max_cell_scalar().unwrap().value);
    assert!((x - y).abs() < 0.01 * y, "split {x}, native {y}");
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn a_sliver_beside_split_cells_is_split_as_they_are() {
    // A quarter ring 10 to 14 across, 0.189 thick, its bore pressed (Lame:
    // von Mises 7.37 at the bore, nu 0). On 0.5 mm cells the slab is split
    // to 0.0625 mm, but a cell holding only a corner of material by the
    // bore (1.2% of it) measures its chords across the bore, not the slab:
    // it was not split, and was dropped, where 0.0625 mm cells hold 37 and
    // 85% of theirs. The bore lost them and read 17.87 (0.0625 mm cells
    // read 7.80). Split as the cells beside it, its children hold them: 9.32.
    let src = "(def r (len2 x y))
(body ring (shape (max (tag bore (- 10 r)) (- r 14) (neg z) (- z 0.189) (neg x) (neg y))))
(study s (shape (shape-of ring)) (elastic (E 200000MPa) (nu 0))
 (fix sx (- x 0.001) (dofs x)) (fix sy (- y 0.001) (dofs y)) (fix sz (- z 0.001) (dofs z))
 (load inside (pressure 2MPa) (region (face bore))))";
    let peak = solve(src, "s", 0.5).max_cell_scalar().unwrap().value;
    assert!(peak < 1.3 * 7.37, "{peak} against 7.37");
}

#[test]
fn a_bar_whose_length_is_no_whole_number_of_cells_reads_its_end() {
    // A 20 x 4 x 4 bar (E 1000, nu 0) held by its end face x = 0 and pulled
    // by 10 MPa at x = 20 stretches 0.2 mm exactly. On cells of 0.3 mm its
    // end face lies between nodes (66.7 cells), and read at the nodes inside
    // alone, its peak was the last node's: 0.19255 +- 0.00448, the three
    // grids short alike (0.1881 on 0.6, 0.1759 on 1.2). Its stated value
    // with its error contains the truth. (Held at x < 1 instead, it read
    // 0.18799 +- 0.00012 against 0.19 and reads 0.18820 +- 0.00016 now, its
    // end read: the rest is the hold's own reach, its grip on the cell its edge cuts.)
    let src = "(def bar (move (box 20 4 4) 10 0 0))
(study s (shape bar) (elastic (E 1000MPa) (nu 0)) (fix wall (face bar (normal -1 0 0)))
  (load tip (pressure -10MPa) (region (face bar (normal 1 0 0)))) (res 0.3))
(require d (max-displacement s < 1000mm))";
    let c = checks(src).remove(0).1;
    let (v, e) = (measured(&c).unwrap(), grid_error(&c).unwrap_or(0.0));
    assert!((v - 0.2).abs() <= e, "{v} +- {e} against 0.2: {c:?}");
    // Heat over time alike: a slab 10.3 thick (x 1.13 to 11.43, its faces between
    // nodes on every grid: the last inside 11.2) heated at both faces by
    // convection (h 500 to 100 C, from 20 C; k 2.5, rho 2000, c 900) for
    // 15 s reaches 68.04 C at its faces (the series solution); read at its
    // nodes inside alone, `< 67.5` held at 66.48 (65.75 to 67.20).
    let src = "(def slab (move (box 10.3 4.1 3.7) 6.28 0.37 -0.41))
(study s (shape slab) (heat (k 2.5W/m/K) (density 2000kg/m3) (capacity 900J/kg/K))
  (convection hot (face slab (normal 1 0 0)) (coefficient 500W/m2/K) (ambient 100degC))
  (convection hot2 (face slab (normal -1 0 0)) (coefficient 500W/m2/K) (ambient 100degC))
  (transient (duration 15s) (steps 60) (initial 20degC)) (res 0.7))
(require q (max-temperature s < 67.5degC))";
    let c = checks(src).remove(0).1;
    let (v, e) = (measured(&c).unwrap(), grid_error(&c).unwrap_or(0.0));
    assert!(c.status != ink::claims::Status::Holds && (v - 68.04).abs() <= e, "the slab's faces reach 68.04 C: {v} +- {e}: {c:?}");
}

#[test]
fn a_peak_temperature_is_read_in_the_cell_its_point_lies_in() {
    // A 7.3 x 2.2 x 1.9 bar (k 0.5) cooled by a film (200 W/m2/K to 20 C)
    // on its x = 0 face and heated by 0.01 W over its last mm reaches
    // 64.4976 C at its far end (1-D, exact). On 0.7 and 1.4 mm cells its
    // end lies in a cell 0.43 full, and read as a stress is, from the full
    // cell under it, that cell's slope carried 0.3 mm past it where the heat
    // curves the field flat: 65.31 and 65.28, where the field at the point
    // reads 64.60 (and 0.35 mm cells 64.52), the grids moving apart.
    let src = "(def bar (move (box 7.3 2.2 1.9) 3.65 0 0))
(study s (shape bar) (heat (k 0.5W/m/K))
  (convection air (- x 0.0001) (coefficient 200W/m2/K) (ambient 20degC))
  (source q (power 0.01W) (region (- 6.3 x))))";
    for res in [0.7, 1.4] {
        let peak = solve(src, "s", res).max_nodal().unwrap();
        assert!((peak.value - 64.4976).abs() < 0.3 && peak.at[0] > 7.0, "on {res} mm cells: {peak:?} against 64.4976 at x = 7.3");
    }
}
