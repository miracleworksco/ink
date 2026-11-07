//! One membership rule (ink_core::sign, docs/CORE.md), read alike by every
//! part that asks which side of a shape a point is on: a point is inside
//! where the field is below 0; a 0 is not inside, a region the field reads
//! as 0 is air, and a 0 the field keeps one sign around is no surface.

use ink::eval::Ctx;
use ink::graph::{Graph, NodeId, UnOp, Var};
use ink::tape::Tape;
use ink::{feature, mesh, query, voids};
use std::path::Path;
use std::sync::Arc;

const ROOT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

fn at0() -> Ctx {
    Ctx::new(0.0, vec![])
}

/// A box about `c`, `h` from it each way: max over the axes of |v - c| - h.
fn cube(g: &mut Graph, c: [f64; 3], h: [f64; 3]) -> NodeId {
    let faces: Vec<NodeId> = (0..3)
        .map(|k| {
            let v = g.var([Var::X, Var::Y, Var::Z][k]);
            let (c, h) = (g.constant(c[k]), g.constant(h[k]));
            let d = g.sub(v, c);
            let a = g.un(UnOp::Abs, d);
            g.sub(a, h)
        })
        .collect();
    let xy = g.max(faces[0], faces[1]);
    g.max(xy, faces[2])
}

fn tape(g: &Graph, f: NodeId) -> Arc<Tape> {
    Arc::new(Tape::new(g, f))
}

/// A 40 mm cube hollowed to 20 mm, its hollow read as 0 (`(min 0 f)`)
/// inside the outer box: air the field reads as 0.
fn zero_hollow(g: &mut Graph) -> NodeId {
    let outer = cube(g, [0.0; 3], [20.0; 3]);
    let inner = cube(g, [0.0; 3], [10.0; 3]);
    let ni = g.neg(inner);
    let hollow = g.max(outer, ni);
    let zero = g.constant(0.0);
    let capped = g.min(zero, hollow);
    g.max(outer, capped)
}

/// Two halves of a 20 mm cube meeting flush on x = 0: 0 across it inside
/// the material.
fn split_cube(g: &mut Graph) -> NodeId {
    let b = cube(g, [0.0; 3], [10.0; 3]);
    let x = g.var(Var::X);
    let nx = g.neg(x);
    let (left, right) = (g.max(b, x), g.max(b, nx));
    g.min(left, right)
}

#[test]
fn air_the_field_reads_as_0_is_a_void() {
    let mut g = Graph::new();
    let f = zero_hollow(&mut g);
    let v = voids::voids(&tape(&g, f), &at0(), [-20.0; 3], [20.0; 3], 1.0, 0.0);
    assert_eq!(v.count(), Some(1), "{v:?}");
}

#[test]
fn a_pocket_flush_with_a_cell_face_seals_no_void() {
    // A cube with a 4 mm pocket whose mouth lies on its top face, which is
    // a face of the count's cells: the field is 0 across the mouth.
    for top in [7.0, 15.0, 31.0, 63.0] {
        let mut g = Graph::new();
        let outer = cube(&mut g, [top / 2.0; 3], [top / 2.0; 3]);
        let pocket = cube(&mut g, [top / 2.0, top / 2.0, top - 2.0], [2.0; 3]);
        let np = g.neg(pocket);
        let f = g.max(outer, np);
        let v = voids::voids(&tape(&g, f), &at0(), [0.0; 3], [top; 3], 1.0, 0.0);
        assert_eq!(v.at_least, 0, "top {top}: {v:?}");
    }
}

#[test]
fn a_zero_sheet_in_air_is_no_wall() {
    for half in [4.0, 8.0, 12.0, 16.0] {
        let mut g = Graph::new();
        let outer = cube(&mut g, [0.0; 3], [20.0; 3]);
        let inner = cube(&mut g, [0.0; 3], [15.0; 3]);
        let ni = g.neg(inner);
        let hollow = g.max(outer, ni);
        let b = cube(&mut g, [0.0; 3], [half; 3]);
        let sheet = g.un(UnOp::Abs, b);
        let f = g.min(hollow, sheet);
        let v = voids::voids(&tape(&g, f), &at0(), [-20.0; 3], [20.0; 3], 1.0, 0.0);
        assert_eq!(v.count(), Some(1), "sheet {half}: {v:?}");
    }
}

#[test]
fn a_subset_is_of_what_b_holds_not_what_it_reads_as_0() {
    use ink::claims::Status::*;
    use ink::claims::{ClaimOptions, run_asserts};
    let src = "(import \"std\")
        (def a (move (box 2 2 2) 10 0 0))
        (def b (min 0 (sphere 1)))
        (assert apart (disjoint a b))
        (assert in-b (subset a b))
        (assert in-sphere (subset a (sphere 1)))
        (assert ball-in-b (subset (sphere 0.5) b))";
    let doc = ink::load_str(src, "subset", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let s: Vec<_> = r.asserts.iter().map(|a| (a.name.as_str(), a.check.status)).collect();
    assert_eq!(s, [("apart", Holds), ("in-b", Fails), ("in-sphere", Fails), ("ball-in-b", Holds)]);
}

#[test]
fn faces_meeting_flush_inside_the_material_are_meshed_only_with_a_report() {
    // The 0 between the halves reads as air (nothing shows it has no
    // thickness): a wall drawn there is reported.
    let mut g = Graph::new();
    let f = split_cube(&mut g);
    let (m, u) = mesh::mesh_checked(&tape(&g, f), &at0(), [-10.0; 3], [10.0; 3], 1.0, mesh::Limit::GRID).unwrap();
    let inner = m.verts.iter().filter(|p| p[0].abs() < 0.01 && p[1].abs() < 9.0 && p[2].abs() < 9.0).count();
    assert!(inner == 0 || u.cells > 0, "{inner} vertices on the 0 inside the material, unreported");
}

/// Names seen and unseen on a shape's surface, at 128 cells across.
fn names(g: &Graph, f: NodeId, lo: [f64; 3], hi: [f64; 3]) -> (Vec<u32>, Vec<(u32, usize, usize)>) {
    feature::surface_names(&tape(g, f), &at0(), lo, hi, 128.0)
}

#[test]
fn a_face_however_small_is_seen() {
    // A sphere 2 um across beside a box: a sign change at that size.
    let mut g = Graph::new();
    let b = cube(&mut g, [-20.0, 0.0, 0.0], [5.0; 3]);
    let terms: Vec<NodeId> = [(Var::X, 20.12345), (Var::Y, 0.12345), (Var::Z, 0.12345)]
        .into_iter()
        .map(|(v, c)| {
            let (v, c) = (g.var(v), g.constant(c));
            let d = g.sub(v, c);
            g.un(UnOp::Square, d)
        })
        .collect();
    let s = g.add(terms[0], terms[1]);
    let s = g.add(s, terms[2]);
    let d = g.un(UnOp::Sqrt, s);
    let r = g.constant(1e-6);
    let sphere = g.sub(d, r);
    let sphere = g.tag("sphere", sphere);
    let f = g.min(b, sphere);
    let (seen, unseen) = names(&g, f, [-25.0, -5.0, -5.0], [25.0, 5.0, 5.0]);
    assert!(!seen.is_empty(), "unseen {unseen:?}");
}

#[test]
fn a_name_only_on_a_zero_inside_the_material_is_not_seen() {
    // An insert filling its pocket exactly: 0 across their faces, negative
    // either side.
    let mut g = Graph::new();
    let b = cube(&mut g, [0.0; 3], [10.0; 3]);
    let insert = cube(&mut g, [0.0; 3], [2.0; 3]);
    let ni = g.neg(insert);
    let ring = g.max(b, ni);
    let insert = g.tag("insert", insert);
    let f = g.min(ring, insert);
    let (seen, _) = names(&g, f, [-10.0; 3], [10.0; 3]);
    assert!(seen.is_empty(), "{seen:?}");
}

#[test]
fn a_slice_reads_zeros_by_the_rule() {
    // Air read as 0 is a hole; a 0 inside the material, drawn, is reported.
    let mut g = Graph::new();
    let f = zero_hollow(&mut g);
    assert_eq!(query::slice(&tape(&g, f), &at0(), 2, 0.0, [-20.0; 3], [20.0; 3], 1.0).len(), 2);
    let mut g = Graph::new();
    let f = split_cube(&mut g);
    let (loops, open) = query::slice_checked(&tape(&g, f), &at0(), 2, 0.0, [-10.0; 3], [10.0; 3], 1.0);
    assert!(loops.len() == 1 || !open.is_empty(), "{} loops, unreported", loops.len());
}

/// A model loaded from source (std imported), and its first body's tape.
fn body(src: &str) -> (ink::lang::Doc, Arc<Tape>) {
    let doc = ink::load_str(&format!("(import \"std\") {src}"), "sign", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let t = tape(&doc.g, doc.bodies[0].shape);
    (doc, t)
}

#[test]
fn air_read_as_0_wider_than_the_samples_stays_air() {
    // A slit 1.5 mm wide whose air reads exactly 0, sampled 1 mm apart:
    // its centre's neighbours are inside, but nearer they read 0.
    let mut g = Graph::new();
    let b = cube(&mut g, [0.0; 3], [10.0; 3]);
    let x = g.var(Var::X);
    let ax = g.un(UnOp::Abs, x);
    let w = g.constant(0.75);
    let cut = g.sub(w, ax);
    let zero = g.constant(0.0);
    let cut = g.min(zero, cut);
    let f = g.max(b, cut);
    let t = tape(&g, f);
    assert_eq!(query::slice(&t, &at0(), 2, 0.3, [-10.0; 3], [10.0; 3], 1.0).len(), 2, "the slit splits the slice");
    let m = mesh::mesh(&t, &at0(), [-10.0; 3], [10.0; 3], 1.0);
    assert_eq!(mesh::analyze(&m).shells.len(), 2, "the slit splits the solid");
}

#[test]
fn a_zero_in_air_beside_a_surface_is_not_seen() {
    // fake-face reads 0 on x = c, positive either side; a real face lies
    // 1e-5 mm from it.
    let (doc, t) = body("(body part (shape (max (box 20 20 20) (min (tag fake-face (max (- x 0.12345) (* -2 (- x 0.12345)))) (+ (- x 0.12345) 0.00001)))))");
    let fake = doc.g.tag_names.iter().position(|n| n == "fake-face").unwrap() as u32;
    let (seen, _) = feature::surface_names(&t, &at0(), [-10.0; 3], [10.0; 3], 128.0);
    assert!(!seen.contains(&fake), "{seen:?}");
}

#[test]
fn material_bounds_leave_out_air_read_as_0() {
    let (_, t) = body("(body capped (shape (min 0 (box 20 20 20))))");
    let (lo, hi, _) = query::bounds(&t, &at0(), [-500.0; 3], [500.0; 3]).unwrap();
    assert!(query::max_extent(lo, hi) < 21.0, "{lo:?}..{hi:?}");
    assert!(feature::survey(&t, &at0(), lo, hi, 128.0, true).interior, "the box has an interior");
}

#[test]
fn a_subset_of_what_reads_0_leaves_the_model_valid() {
    // Into a constant 0 and a param read as 0: no number enters the model
    // (the question is a's inside and b's outside, not a field shifted).
    use ink::claims::Status::*;
    use ink::claims::{ClaimOptions, run_asserts};
    let src =
        "(import \"std\") (param air 0) (def a (box 2 2 2)) (assert into-empty (subset a 0)) (assert into-air (subset a air)) (assert normal (subset (sphere 0.5) (sphere 1)))";
    let doc = ink::load_str(src, "subset-0", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let s: Vec<_> = r.asserts.iter().map(|a| (a.name.as_str(), a.check.status)).collect();
    assert_eq!(s, [("into-empty", Fails), ("into-air", Fails), ("normal", Holds)]);
}

#[test]
fn a_volume_derivative_reads_the_inside_where_outside_reads_0() {
    // A 2 x 2 x 2p box, p = 1: dV/dp = 8, its outside read as 0 or not.
    let mut g = Graph::new();
    let p = g.var(Var::Param(0));
    let one = g.constant(1.0);
    let faces: Vec<NodeId> = [(Var::X, p), (Var::Y, one), (Var::Z, one)]
        .into_iter()
        .map(|(v, h)| {
            let v = g.var(v);
            let a = g.un(UnOp::Abs, v);
            g.sub(a, h)
        })
        .collect();
    let f = g.max(faces[0], faces[1]);
    let f = g.max(f, faces[2]);
    let zero = g.constant(0.0);
    let capped = g.min(f, zero);
    let ctx = Ctx::new(0.0, vec![1.0]);
    let weight = |_: [f64; 3]| 1.0;
    for f in [f, capped] {
        let (area, d) = ink::design::shape_derivatives(&tape(&g, f), &[&weight], &ctx, [-1.5; 3], [1.5; 3], 0.25);
        assert!((area - 24.0).abs() < 0.5 && (d[0][0] - 8.0).abs() < 0.2, "area {area}, dV/dp {}", d[0][0]);
    }
}

#[test]
fn a_0_where_air_wedges_meet_is_air_and_where_blocks_meet_is_undecided() {
    // min(0, |x + y| / 2 - |x - y|): air read as 0 in two wedges meeting
    // on the z axis; -min(|x|, |y|): four blocks meeting there flush.
    let mut g = Graph::new();
    let (x, y) = (g.var(Var::X), g.var(Var::Y));
    let (s, d) = (g.add(x, y), g.sub(x, y));
    let (s, d) = (g.un(UnOp::Abs, s), g.un(UnOp::Abs, d));
    let half = g.constant(0.5);
    let s = g.mul(half, s);
    let w = g.sub(s, d);
    let zero = g.constant(0.0);
    let wedges = g.min(zero, w);
    let (ax, ay) = (g.un(UnOp::Abs, x), g.un(UnOp::Abs, y));
    let m = g.min(ax, ay);
    let blocks = g.neg(m);
    for (f, decided) in [(wedges, true), (blocks, false)] {
        assert_eq!(ink::sign::undecided(&tape(&g, f), &at0(), &[[0.0; 3]], &[0.0], None).is_empty(), decided);
    }
}

#[test]
fn four_blocks_meeting_flush_mesh_and_slice_as_one_or_say_so() {
    let (_, t) = body("(body p (shape (max (box 20 20 20) (neg (min (abs x) (abs y))))))");
    let (m, u) = mesh::mesh_checked(&t, &at0(), [-10.0; 3], [10.0; 3], 1.0, mesh::Limit::GRID).unwrap();
    let inner = m.verts.iter().filter(|p| p[0].abs() < 1.0 && p[1].abs() < 1.0 && p[2].abs() < 9.0).count();
    assert!(inner == 0 || u.cells > 0, "{inner} vertices where the blocks meet, unreported");
    let (loops, open) = query::slice_checked(&t, &at0(), 2, 0.3, [-10.0; 3], [10.0; 3], 1.0);
    assert!(loops.len() == 1 || !open.is_empty());
}

#[test]
fn a_subset_reads_a_where_its_material_lies() {
    // b reads 0 everywhere: max(a, -b) is never below 0, so its own reach
    // is the default window, which a 2 mm cube at x = 1000 lies outside.
    use ink::claims::Status::*;
    use ink::claims::{ClaimOptions, run_asserts};
    let src = "(import \"std\") (def a (move (box 2 2 2) 1000 0 0)) (assert contains (subset a 0))";
    let doc = ink::load_str(src, "far", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    assert_eq!(r.asserts[0].check.status, Fails, "{:?}", r.asserts[0].check);
}

#[test]
fn air_read_as_0_no_bound_shows_is_never_filled_unreported() {
    // A 0.4 mm slit read as 0, and a 4 mm one whose field cancels
    // ((x + y) - x - y, 0 in real numbers, loose in bounds), sampled 1 mm
    // apart: neither is shown a join, so each is air or reported.
    for (src, half) in [("(min 0 (- 0.2 (abs x)))", 4.0), ("(min 0 (+ (- 2 (abs x)) (* 16 (- (- (+ x y) x) y))))", 6.0)] {
        let (_, t) = body(&format!("(body p (shape (max (box {w} {w} {w}) {src})))", w = 2.0 * half));
        let open = ink::sign::undecided(&t, &at0(), &[[0.0; 3]], &[0.0], None);
        let (m, u) = mesh::mesh_checked(&t, &at0(), [-half; 3], [half; 3], 1.0, mesh::Limit::GRID).unwrap();
        assert!(mesh::analyze(&m).shells.len() == 2 || u.cells > 0, "{src}: one solid, nothing unresolved ({open:?})");
        let (loops, undecided) = query::slice_checked(&t, &at0(), 2, 0.3, [-half; 3], [half; 3], 1.0);
        assert!(loops.len() == 2 || !undecided.is_empty(), "{src}: one loop, nothing unresolved");
    }
}

#[test]
fn an_undecided_saddle_face_centre_reaches_the_mesh_report() {
    let (_, t) = body("(body p (shape (max (box 7 7 7) (- (* 0.01 (abs (- (+ x y) 1))) (abs (- x y))))))");
    let p = [0.5, 0.5, 0.0];
    let v = [ink::eval::eval_point(&t, &at0(), p)];
    assert_eq!(ink::sign::undecided(&t, &at0(), &[p], &v, None).len(), 1, "a sample nothing near decides");
    let (_, u) = mesh::mesh_checked(&t, &at0(), [-4.0; 3], [4.0; 3], 1.0, mesh::Limit::GRID).unwrap();
    assert!(u.cells > 0, "{u:?}");
}

#[test]
fn a_0_with_nothing_beside_it_shown_is_reported_whatever_its_width() {
    // A join (-|x|), and air read as 0 half a nanometre wide or 0.05 mm
    // wide 100 km out: nothing near shows air beside the sample, so each is
    // reported (never read as material: there is no such rule).
    for (src, p) in [("(neg (abs x))", [0.0; 3]), ("(min 0 (- 0.25e-9 (abs x)))", [0.0; 3]), ("(min 0 (- 0.025 (abs x)))", [0.0, 1e8, 0.0])] {
        let (_, t) = body(&format!("(body p (shape {src}))"));
        let v = [ink::eval::eval_point(&t, &at0(), p)];
        assert_eq!(ink::sign::undecided(&t, &at0(), &[p], &v, None), vec![p], "{src}");
    }
    let (_, t) = body("(body p (shape (max (move (box 8 8 8) 0 100000000 0) (min 0 (- 0.025 (abs x))))))");
    let (lo, hi) = ([-4.0, 1e8 - 4.0, -4.0], [4.0, 1e8 + 4.0, 4.0]);
    let (m, u) = mesh::mesh_checked(&t, &at0(), lo, hi, 1.0, mesh::Limit::GRID).unwrap();
    assert!(mesh::analyze(&m).shells.len() == 2 || u.cells > 0, "{u:?}");
    let (loops, open) = query::slice_checked(&t, &at0(), 2, 0.3, lo, hi, 1.0);
    assert!(loops.len() == 2 || !open.is_empty());
}

#[test]
fn a_0_a_crossing_s_halving_reads_reaches_the_mesh_report() {
    // No sample of the grid is on the angled join x + y = 3.25; the first
    // halving of the edge (3, 0, 0)..(4, 0, 0) reads (3.25, 0, 0).
    let (_, t) = body("(body p (shape (max (box 7.8125 7.8125 7.8125) (neg (* 0.125 (abs (- (+ x y) 3.25)))))))");
    let (_, u) = mesh::mesh_checked(&t, &at0(), [-4.0; 3], [4.0; 3], 1.0, mesh::Limit::GRID).unwrap();
    assert!(u.cells > 0, "{u:?}");
}

#[test]
fn faces_one_float_apart_leave_air_that_is_reported_not_filled() {
    // Air between x = 1 and the next float: a sample at x = 1 reads 0, inside
    // either side a few billionths away. Nothing shows its 0 has no width.
    let (_, t) = body("(body p (shape (max (box 7.5 7.5 7.5) (min (- x 1) (- 1.0000000000000002 x)))))");
    let (_, open) = query::slice_checked(&t, &at0(), 2, 0.3, [-4.0; 3], [4.0; 3], 1.0);
    assert!(open.iter().any(|p| p[0] == 1.0), "{open:?}");
}

#[test]
fn a_mark_s_report_does_not_depend_on_its_region_s_scale() {
    // The same region, its field scaled by 100: the same marks, the same
    // undecided saddle centre (0.5, 0.5, 0) reported.
    let (_, g) = body("(body p (shape (- (* 0.01 (abs (- (+ x y) 1))) (abs (- x y)))))");
    let got: Vec<(usize, usize)> = [1, 100]
        .iter()
        .map(|k| {
            let (_, region) = body(&format!("(body p (shape (* {k} (- (max (abs x) (abs y)) 0.05))))"));
            let (marks, open) = query::slice_mark(&g, &region, &at0(), 2, 0.0, [-2.0; 3], [2.0; 3], 1.0);
            (marks.len(), open.len())
        })
        .collect();
    assert!(got[0].0 > 0 && got[0] == got[1] && got[0].1 > 0, "{got:?}");
}

#[test]
fn a_sweep_meshes_in_seconds() {
    // Halving its crossings lands on exact 0s of its near-linear field:
    // each is read by its neighbours, not a proof over the field.
    let doc = ink::load_str(&std::fs::read_to_string(Path::new(ROOT).join("benchmarks/19-spline-sweep-loft.ink")).unwrap(), "duct", Path::new(ROOT)).unwrap();
    let f = doc.bodies.iter().find(|b| b.name == "duct").unwrap().shape;
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let (lo, hi, _) = v.bounds(f).unwrap();
    let t = v.clipped(f);
    let start = std::time::Instant::now();
    let (m, _) = mesh::mesh_checked(&t, &v.ctx, lo, hi, 0.6, mesh::Limit::GRID).unwrap();
    let s = start.elapsed().as_secs_f64();
    assert!(s < 30.0 && mesh::analyze(&m).shells.len() == 1, "{s} s");
}

#[test]
fn a_join_on_a_brick_s_face_is_read_past_it_by_the_whole_field() {
    // -|x|'s brick over x >= 0 holds -x, which reads air past x = 0 where
    // the whole field reads material: the join at 0 stays undecided.
    let mut g = Graph::new();
    let x = g.var(Var::X);
    let nx = g.neg(x);
    let ax = g.un(UnOp::Abs, x);
    let whole = g.neg(ax);
    let within = [ink::interval::I::new(0.0, 8.0), ink::interval::I::new(-4.0, 4.0), ink::interval::I::new(-4.0, 4.0)];
    let (brick, whole) = (tape(&g, nx), tape(&g, whole));
    assert_eq!(ink::sign::undecided(&brick, &at0(), &[[0.0; 3]], &[0.0], Some((within, &whole))), vec![[0.0; 3]]);
    // Flush joins on a brick boundary (the grid starts at -8: x = 0 is
    // one), meshed: no wall, or one reported.
    for src in ["(neg (abs x))", "(min x (neg x))", "(neg (square x))", "(max (select x x (* -2 x)) (select (neg x) (* -2 x) x))"] {
        let (_, t) = body(&format!("(body p (shape (max (box 14 14 14) {src})))"));
        let (m, u) = mesh::mesh_checked(&t, &at0(), [-7.0; 3], [7.0; 3], 1.0, mesh::Limit::GRID).unwrap();
        let inner = m.verts.iter().filter(|p| p[0].abs() < 0.2 && p[1].abs() < 6.0 && p[2].abs() < 6.0).count();
        assert!(inner == 0 || u.cells > 0, "{src}: {inner} vertices on the join, unreported");
    }
}

#[test]
fn a_subset_s_range_is_checked_where_a_lies() {
    // b reads 0 everywhere, so max(a, -b) is never below 0 and reaches
    // nothing; a, a cube 2e9 mm out, is past CORE's range: unknown, as its
    // search reads a's reach.
    use ink::claims::Status::*;
    use ink::claims::{ClaimOptions, run_asserts};
    let src = "(import \"std\") (def a (move (box 2 2 2) 2000000000 0 0)) (assert contained (subset a 0) (accuracy 0.05))";
    let doc = ink::load_str(src, "far-range", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let c = &r.asserts[0].check;
    assert!(c.status == Unknown && c.accuracy.is_none(), "{c:?}");
}
