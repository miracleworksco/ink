//! Standard-library geometry checked against exact answers.

use ink::eval::{Ctx, eval_point, eval_points};
use ink::lang::{Doc, Value};
use ink::tape::Tape;
use std::f64::consts::PI;
use std::path::Path;

/// The repository: what the tests load (catalog parts, fixtures) is found from it.
pub(crate) const ROOT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

pub(crate) fn field(src: &str, name: &str) -> (Doc, Tape) {
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(f, _) = d.def(name).unwrap_or_else(|| panic!("no def {name}")).clone() else { panic!("{name} is not a field") };
    let t = Tape::new(&d.g, f);
    (d, t)
}

const CTX: Ctx = Ctx::new(0.0, Vec::new());

/// No point of b within d of a, covering b's box lo..hi with boxes down to
/// `finest` (query::clearance, one way round): the verdict and a pair.
fn clearance_of(
    a: &std::sync::Arc<Tape>,
    b: &std::sync::Arc<Tape>,
    ctx: &Ctx,
    d: f64,
    lo: [f64; 3],
    hi: [f64; 3],
    finest: f64,
) -> (ink::query::Verdict, Option<([f64; 3], [f64; 3])>) {
    let way = ink::query::Way { grown: a.clone(), covered: b.clone(), lo, hi, swapped: false };
    let (v, pair, _) = ink::query::clearance(&[way], ctx, d, finest, ink::query::BOXES);
    (v, pair)
}

/// An exact emptiness proof at cells of h (None: a 32nd of the material's
/// box; query::prove_empty): its verdict and the cells it reached.
fn empty_at(tape: &std::sync::Arc<Tape>, ctx: &Ctx, lo: [f64; 3], hi: [f64; 3], h: Option<f64>) -> (ink::query::Verdict, f64) {
    let s = ink::query::prove_empty(tape, ctx, lo, hi, h, h.unwrap_or(f64::INFINITY), 0.0);
    (s.verdict, s.res)
}

#[test]
fn static_instance_bounds_cover_the_inverse_transpose_shape() {
    // R is frozen independently of axis normalization: diag(0.75, 1, 1).
    // The world half-extent is 5 / 0.75 = 6.66667, whereas the forward
    // map gives 5 * 0.75 = 3.75.
    let mut d = ink::load_str(
        "(import \"std\") (body cube (shape (box 10 10 10)))
         (assembly root (instance c cube))",
        "test",
        Path::new(ROOT),
    )
    .unwrap();
    d.assemblies[0].instances[0].at.r[0][0] = d.g.constant(0.75);
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(d, ctx);
    let it = v.doc.items("root/c").unwrap().remove(0);
    // A scale is no rigid motion (its moments, widths and surface are not
    // the body's moved): the instance is analysed in the world.
    assert!(v.rigid(&it, &[]).is_none() && v.local(&it).is_none());
    let (lo, hi, _) = v.item_bounds(&it).unwrap();
    let world = v.doc.placed_field(&it);
    let tape = Tape::new(&v.doc.g, world);
    for p in [[6.65, 0.0, 0.0], [-6.65, 0.0, 0.0], [0.0, 4.99, 0.0], [0.0, 0.0, -4.99]] {
        assert!(eval_point(&tape, &v.ctx, p) < 0.0);
        assert!((0..3).all(|k| lo[k] <= p[k] && p[k] <= hi[k]), "{p:?} outside {lo:?}..{hi:?}");
    }
    // Where the Gram bound cannot establish invertibility, use the
    // world field. Here R = 0: the cube's world field is negative everywhere.
    let mut singular = ink::lang::Item { label: it.label.clone(), shape: it.shape, attrs: it.attrs.clone(), instance: it.instance.clone(), links: it.links.clone() };
    singular.instance.as_mut().unwrap().1.r = [[v.doc.g.constant(0.0); 3]; 3];
    assert!(v.local(&singular).is_none());
    singular.instance.as_mut().unwrap().1.r[0][0] = v.doc.g.constant(f64::INFINITY);
    assert!(v.local(&singular).is_none());
    // Identity placements keep the body's exact box.
    let mut identity = it;
    identity.instance.as_mut().unwrap().1.r = std::array::from_fn(|i| std::array::from_fn(|j| v.doc.g.constant(if i == j { 1.0 } else { 0.0 })));
    let local = v.doc.body("cube").unwrap().shape;
    assert_eq!(v.item_bounds(&identity), v.bounds(local));
    // At a fixed moment t is still a real input. A sin/cos turn's Gram
    // matrix is exactly I there, but its f64 point matrix need not be: the
    // instance is analysed in its frame all the same, its box widened by
    // how far the entries may be from their floats, and holds the cube.
    let d = ink::load_str(
        "(import \"std\") (body cube (shape (box 10 10 10)))
        (assembly root (instance c cube (at (rot-z (* 30 t)))) (instance translated cube (at (move t 0 0))))",
        "test",
        Path::new(ROOT),
    )
    .unwrap();
    let ctx = d.ctx(&[], 1.0).unwrap();
    let mut v = ink::view::ModelView::new(d, ctx);
    let it = v.doc.items("root/c").unwrap().remove(0);
    assert!(v.local(&it).is_some());
    let (lo, hi, _) = v.item_bounds(&it).unwrap();
    let (s, c) = 30f64.to_radians().sin_cos();
    for (x, y) in [(5.0, 5.0), (5.0, -5.0), (-5.0, 5.0), (-5.0, -5.0)] {
        let p = [c * x - s * y, s * x + c * y, 5.0];
        assert!((0..3).all(|k| lo[k] <= p[k] && p[k] <= hi[k]), "{p:?} outside {lo:?}..{hi:?}");
    }
    let exact = v.doc.items("root/translated").unwrap().remove(0);
    assert!(v.local(&exact).is_some());
}

#[test]
fn shape_of_except_leaves_out_the_parts_named() {
    // An assembly but one of its parts (a tool's surroundings but the screw
    // it turns): the union of the others, no surface left where it was.
    let src = "(body cube (shape (box 2 2 2)))
               (assembly pair (instance a cube) (instance b cube (at (move 3 0 0))))
               (assembly root (instance p pair) (instance c cube (at (move 0 3 0))))
               (def rest (shape-of root (except root/p/b)))
               (def none-of-p (shape-of root (except root/p)))";
    let (_, t) = field(src, "rest");
    let at = |t: &Tape, p: [f64; 3]| eval_point(t, &CTX, p);
    assert!(at(&t, [0.0, 0.0, 0.0]) < 0.0, "a is kept");
    assert!(at(&t, [0.0, 3.0, 0.0]) < 0.0, "c is kept");
    assert!((at(&t, [3.0, 0.0, 0.0]) - 2.0).abs() < 1e-9, "b's centre is 2 from a, not inside b");
    let (_, t) = field(src, "none-of-p");
    assert!(at(&t, [0.0, 0.0, 0.0]) > 0.0 && at(&t, [0.0, 3.0, 0.0]) < 0.0, "a path leaves out all under it");
    // A path leaves out what is at or under it past a '/', not every label it begins: k-1 is not k-10.
    let (_, t) = field(
        "(body cube (shape (box 2 2 2)))
        (assembly root (for k (vec 1 10 11) (instance k-$k cube (at (move (* 3 k) 0 0)))))
        (def rest (shape-of root (except root/k-1)))",
        "rest",
    );
    assert!(at(&t, [3.0, 0.0, 0.0]) > 0.0 && at(&t, [30.0, 0.0, 0.0]) < 0.0 && at(&t, [33.0, 0.0, 0.0]) < 0.0, "k-1 out, k-10 and k-11 kept");
    // A `for` among the paths repeats them for a list known at load time: a
    // build order's parts from step k on (everything not yet placed), and a
    // list that makes none.
    let order = "(body cube (shape (box 2 2 2)))
        (assembly root (for k (range 4) (instance p-$k cube (at (move (* 3 k) 0 0)))))
        (def order (vec \"p-0\" \"p-1\" \"p-2\" \"p-3\"))";
    let (_, t) = field(&format!("{order} (def placed (shape-of root (except (for q (slice order 2) root/$q) (for q (vec) root/$q))))"), "placed");
    assert!(at(&t, [0.0, 0.0, 0.0]) < 0.0 && at(&t, [3.0, 0.0, 0.0]) < 0.0, "p-0 and p-1 placed");
    assert!(at(&t, [6.0, 0.0, 0.0]) > 0.0 && at(&t, [9.0, 0.0, 0.0]) > 0.0, "p-2 and p-3 not yet");
    // Its list is read where the shape-of is: a function's argument.
    let (_, t) = field(&format!("{order} (fn placed-before (k) (shape-of root (except (for q (slice order k) root/$q)))) (def placed (placed-before 3))"), "placed");
    assert!(at(&t, [6.0, 0.0, 0.0]) < 0.0 && at(&t, [9.0, 0.0, 0.0]) > 0.0, "p-2 placed before step 3, p-3 not");
    // Imported under a namespace, the loop's variable stays the loop's though
    // the file defines the same name, and a claim named after a part's path
    // (`(assert $p ...)` over "moving/b", a sub-assembly's part) is the
    // file's, as any claim it names.
    let dir = crate::common::Scratch::new("except-for");
    std::fs::write(
        dir.join("m.ink"),
        "(import \"std\") (def q \"unused\") (body cube (shape (box 2 2 2)))
        (assembly s (instance a cube) (instance b cube (at (move 5 0 0))))
        (def left (shape-of s (except (for q (vec \"b\") s/$q))))
        (for p (vec \"moving/b\") (assert $p (disjoint (shape-of s/b) (shape-of s/a))))",
    )
    .unwrap();
    let d = ink::load_str("(import \"m.ink\" as m) (assert moving/b (value 1 > 0))", "test", &dir).unwrap_or_else(|e| panic!("{e:?}"));
    let Some(Value::Field(f, _)) = d.def("m.left").cloned() else { panic!("m.left is not a field") };
    let t = Tape::new(&d.g, f);
    assert!(at(&t, [0.0, 0.0, 0.0]) < 0.0 && at(&t, [5.0, 0.0, 0.0]) > 0.0, "a kept, b left out");
    let names: Vec<&str> = d.asserts.iter().map(|a| a.name.as_str()).collect();
    assert!(names.contains(&"m.moving/b") && names.contains(&"moving/b"), "{names:?}");
    let e = ink::load_str(&format!("(import \"std\") {order} (def r (shape-of root (except (for q (vec) root/$q))))"), "test", Path::new(ROOT))
        .err()
        .expect("an except naming nothing loads");
    assert!(format!("{e:?}").contains("(except) names no path"), "{e:?}");
    let two = "(body cube (shape (box 2 2 2))) (assembly pair (instance a cube) (instance b cube (at (move 3 0 0))))
               (assembly root (instance p pair) (instance c cube (at (move 0 3 0))))";
    for (bad, says) in [
        ("(shape-of root (except root/q))", "'root/q' names no part of 'root'"),
        ("(shape-of root (except))", "(except) names no path"),
        ("(shape-of root (except root/c root/c))", "'root/c' is named twice"),
        ("(shape-of root (except root/p root/p/b))", "'root/p/b' is under 'root/p', already left out"),
        ("(shape-of root (except root/p/b root/p))", "'root/p' holds 'root/p/b', named too"),
        ("(shape-of root (except root))", "leaves no part of 'root'"),
        ("(shape-of root (except root/p root/c))", "leaves no part of 'root'"),
        ("(shape-of root/p (except root))", "leaves no part of 'root/p'"),
        ("(shape-of cube (except cube))", "leaves no part of 'cube'"),
    ] {
        let e = ink::load_str(&format!("(import \"std\")\n{two} (def r {bad})"), "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{bad} loads"));
        assert!(format!("{e:?}").contains(says), "{bad}: {e:?}");
    }
}

#[test]
fn attr_of_reads_the_placed_bodys_attribute() {
    // (attr-of PATH NAME): what the body placed there says of itself, at load time.
    let src = "(import \"std\") (body nut (shape (box 1 1 1)) (attr spec (read-json \"tests/fixtures/cap-screw.json\")) (attr note \"m5\"))
               (body key (shape (box 1 1 1)))
               (assembly root (instance n nut) (instance k key (at (move 3 0 0))))
               (def size (get (get (attr-of root/n spec) \"drive\") \"size\")) (def note (attr-of nut note)) (param p 3)";
    let d = ink::load_str(src, "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let Some(Value::Field(f, _)) = d.def("size").cloned() else { panic!("size is not a number") };
    assert_eq!(d.g.as_const(f), Some(2.5));
    assert!(matches!(d.def("note"), Some(Value::Str(s)) if s.as_str() == "m5"));
    // An attribute keeps its name in a file imported under a namespace, though
    // the file defines the same name (a catalog part's `(attr spec spec)`).
    let d = ink::load_str("(import \"tests/fixtures/cap-screw.ink\" as m3) (def k (get (get (attr-of m3.part spec) \"drive\") \"kind\"))", "test", Path::new(ROOT))
        .unwrap_or_else(|e| panic!("{e:?}"));
    assert!(matches!(d.def("k"), Some(Value::Str(s)) if s.as_str() == "hex-socket"));
    // A field attribute is read where the part lies, as (attr NAME) reads it.
    let d = ink::load_str(
        "(import \"std\") (body k (shape (box 1 1 1)) (attr heat (- x 0)))
        (assembly root (instance k k (at (move 30 0 0)))) (def heat (attr-of root/k heat))",
        "test",
        Path::new(ROOT),
    )
    .unwrap_or_else(|e| panic!("{e:?}"));
    let Some(Value::Field(f, _)) = d.def("heat").cloned() else { panic!("heat is not a field") };
    assert!(eval_point(&Tape::new(&d.g, f), &CTX, [30.0, 0.0, 0.0]).abs() < 1e-12, "the part's own x = 0 where it lies");
    for (bad, says) in [
        ("(attr-of root/k spec)", "body 'key' (at 'root/k') has no attribute 'spec'"),
        ("(attr-of root spec)", "'root' places 2 bodies"),
        ("(met-by root/k)", "'root/k' is placed by no frame of its own (no CFRAME)"),
        ("(met-by root/nope)", "assembly 'root' has no instance 'nope'"),
        ("(error \"p=\" p \" sum=\" (+ p 1) \" v=\" (vec 1 2.5) \" \" (get (attr-of root/n spec) \"family\"))", "p=3 sum=4 v=[1, 2.5] ISO4762"),
    ] {
        let e = ink::load_str(&format!("{src} (def bad {bad})"), "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{bad} loads"));
        assert!(format!("{e:?}").contains(says), "{bad}: {e:?}");
    }
}

#[test]
fn met_by_names_the_frame_a_part_is_placed_by() {
    // (met-by PATH) is the name of the frame the instance was placed by, as
    // a tool's fit to a fastener reads it (the end it engages): a 10x14
    // spanner placed on a nut by either end.
    for end in ["engage", "engage-large"] {
        let d = ink::load_str(
            &format!(
                "(import \"std\")
            (body nut (shape (cylinder-z 5 5)) (frame drive (at (move 0 0 2.5))))
            (body spanner (shape (move (box 120 12 4) 60 0 -2)) (frame engage (at)) (frame engage-large (at (move 120 0 0))))
            (assembly root (instance n nut) (instance w spanner (frame-of n drive) {end}))
            (def end (met-by root/w))"
            ),
            "met",
            Path::new(ROOT),
        )
        .unwrap_or_else(|e| panic!("{end}: {e:?}"));
        assert!(matches!(d.def("end"), Some(Value::Str(s)) if s.as_str() == end), "{end}: {:?}", d.def("end"));
    }
    // A record's numbers are numbers (24.0 is 24); a size written as text is
    // no number, and arithmetic on it says so.
    let d = ink::load_str("(def nut24 (read-json \"tests/fixtures/drive-24.json\")) (def same (== (get (get nut24 \"drive\") \"size\") 24))", "met", Path::new(ROOT)).unwrap();
    let Some(Value::Field(f, _)) = d.def("same").cloned() else { panic!("same is not a number") };
    assert_eq!(d.g.as_const(f), Some(1.0));
    let e = ink::load_str("(def text (read-json \"tests/fixtures/drive-text.json\")) (def off (- (get (get text \"drive\") \"size\") 2.5))", "met", Path::new(ROOT))
        .err()
        .expect("a string is no number");
    assert!(format!("{e:?}").contains("'-' expects numbers or vecs, got the string \\\"2.50"), "{e:?}");
}

#[test]
fn items_place_no_copy_of_a_body() {
    // 64 copies of a pair of one body a thousand nodes big: resolving them
    // and every one's bounds builds no copy of its field, only placements'
    // numbers; a placed field is built where read, once.
    let balls: String = (0..150).map(|i| format!(" (move (sphere 0.4) {i} 0 0)")).collect();
    let copies: String = (0..64).map(|i| format!(" (instance u{i} pair (at (move 0 0 {})))", 10 * i)).collect();
    let src = format!(
        "(import \"std\") (body b (shape (min{balls})))
         (assembly pair (instance p b) (instance q b (at (rot-z 30) (move 0 50 0))))
         (assembly root{copies})"
    );
    let d = ink::load_str(&src, "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(d, ctx);
    let body = v.doc.g.reachable(&[v.doc.body("b").unwrap().shape]).len();
    let before = v.doc.g.nodes.len();
    let items = v.doc.items("root").unwrap();
    let boxes: Vec<(_, _, _)> = items.iter().map(|it| v.item_bounds(it).unwrap()).collect();
    let grown = v.doc.g.nodes.len() - before;
    assert!(items.len() == 128 && body > 1000 && grown < 16 * items.len(), "{} items of a {body}-node body grew the graph by {grown} nodes", items.len());
    assert!((boxes[2].0[2] - boxes[0].0[2] - 10.0).abs() < 1e-6, "{boxes:?}");
    let f = v.doc.placed_field(&items[3]);
    assert!(v.doc.placed_field(&items[3]) == f && f != items[3].shape);
}

fn seg_dist(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab: Vec<f64> = (0..3).map(|i| b[i] - a[i]).collect();
    let ap: Vec<f64> = (0..3).map(|i| p[i] - a[i]).collect();
    let t = ((0..3).map(|i| ap[i] * ab[i]).sum::<f64>() / (0..3).map(|i| ab[i] * ab[i]).sum::<f64>()).clamp(0.0, 1.0);
    (0..3).map(|i| (ap[i] - t * ab[i]).powi(2)).sum::<f64>().sqrt()
}

struct Lcg(u64);
impl Lcg {
    fn f(&mut self, lo: f64, hi: f64) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        lo + (hi - lo) * ((self.0 >> 11) as f64 / (1u64 << 53) as f64)
    }
}

#[test]
fn polyline_distance_is_exact() {
    let pts = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 5.0, 3.0], [2.0, 8.0, -4.0]];
    let (_d, t) = field("(def f (polyline-distance (vec (vec 0 0 0) (vec 10 0 0) (vec 10 5 3) (vec 2 8 -4))))", "f");
    let mut r = Lcg(1);
    for _ in 0..500 {
        let p = [r.f(-5.0, 15.0), r.f(-5.0, 12.0), r.f(-8.0, 8.0)];
        let want = pts.windows(2).map(|w| seg_dist(p, w[0], w[1])).fold(f64::INFINITY, f64::min);
        assert!((eval_point(&t, &CTX, p) - want).abs() < 1e-9);
    }
}

#[test]
fn polygon_is_a_signed_distance_bound_and_polygon_distance_exact_outside() {
    // L-shape: [0,10]x[0,4] union [0,4]x[0,10], listed counter-clockwise.
    // polygon: the sign exact, the value never over the distance (a bound
    // past corners: an offset mitres them). polygon-distance: the same
    // sign, and the exact distance outside (an offset rounds them).
    let pts = "(vec (vec 0 0) (vec 10 0) (vec 10 4) (vec 4 4) (vec 4 10) (vec 0 10))";
    let (_d, t) = field(&format!("(def f (polygon {pts}))"), "f");
    let (_e, te) = field(&format!("(def f (polygon-distance {pts}))"), "f");
    let inside = |x: f64, y: f64| (0.0..=10.0).contains(&x) && (0.0..=4.0).contains(&y) || (0.0..=4.0).contains(&x) && (0.0..=10.0).contains(&y);
    let verts = [[0.0, 0.0], [10.0, 0.0], [10.0, 4.0], [4.0, 4.0], [4.0, 10.0], [0.0, 10.0]];
    let mut r = Lcg(2);
    for _ in 0..2000 {
        let (x, y) = (r.f(-3.0, 13.0), r.f(-3.0, 13.0));
        let d = (0..6)
            .map(|i| {
                let (a, b) = (verts[i], verts[(i + 1) % 6]);
                seg_dist([x, y, 0.0], [a[0], a[1], 0.0], [b[0], b[1], 0.0])
            })
            .fold(f64::INFINITY, f64::min);
        let want = if inside(x, y) { -d } else { d };
        let got = eval_point(&t, &CTX, [x, y, 7.0]);
        assert!(got.abs() <= d + 1e-9 && (got < 0.0) == (want < 0.0) || d < 1e-9, "at ({x}, {y}): {got} vs {want}");
        let exact = eval_point(&te, &CTX, [x, y, 7.0]);
        assert!((exact < 0.0) == (want < 0.0) || d < 1e-9, "at ({x}, {y}): {exact} vs {want}");
        if want > 0.0 {
            assert!((exact - want).abs() < 1e-9, "outside at ({x}, {y}): {exact} vs {want}");
        }
    }
    // Winding order does not matter.
    let (_d, t2) = field("(def f (polygon (vec (vec 0 10) (vec 4 10) (vec 4 4) (vec 10 4) (vec 10 0) (vec 0 0))))", "f");
    assert!((eval_point(&t2, &CTX, [2.0, 8.0, 0.0]) + 2.0).abs() < 1e-9);
}

/// std's slope promises hold for their constants as computed, not only in
/// real numbers: no sampled gradient is over 1. A cone of 5e7 to 1 (its
/// rim's wedge 2e-8 from flat) had slope 1.34 near its axis: its normals
/// came through sqrt(1 - cs^2). Thin cones each way at 1e3 to 1e9, an
/// ordinary cone and a frustum, the hexagon, a gyroid and a GT pulley.
#[test]
fn std_promised_slopes_hold_for_their_floats() {
    let mut fields: Vec<(String, f64)> = Vec::new();
    for n in 3..=9 {
        let big = 1.0 + 10f64.powi(n);
        fields.push((format!("(cone-z 1 {big} 0 1)"), 2.0));
        fields.push((format!("(cone-z {big} 1 0 1)"), 2.0));
        fields.push((format!("(cone-z 0 1 0 {big})"), 2.0));
    }
    for f in ["(cone-z 0 1 0 1)", "(cone-z 3 1 -2 5)", "(cone-z 2 2 0 1)", "(extrude-z (hexagon-2d 7) -1 1)", "(gyroid 10)", "(extrude-z (gt-pulley-teeth-2d 20 2 0) -1 1)"] {
        fields.push((f.to_string(), 8.0));
    }
    let mut r = Lcg(5);
    for (f, size) in &fields {
        let (_d, t) = field(&format!("(def f {f})"), "f");
        let mut pts: Vec<[f64; 3]> = (0..3000).map(|_| [r.f(-size, *size), r.f(-size, *size), r.f(-0.5 * size, 1.5 * size)]).collect();
        pts.extend((0..1000).map(|_| [r.f(-0.1, 0.1), r.f(-0.1, 0.1), r.f(-0.5, 1.5)]));
        pts.push([0.0, 0.0, 0.46]);
        for (p, g) in pts.iter().zip(ink::eval::eval_grads(&t, &CTX, &pts)) {
            let slope = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
            assert!(slope <= 1.0 + 1e-9, "{f} at {p:?}: slope {slope}");
        }
    }
}

/// A cone sized by a param keeps its promise: its constants are finite. A
/// scan for infinities in its graph found norm2-up's guard, inf, in the
/// branch its finite length never takes, and dropped it (11x slower
/// clearance). Near the side, interval
/// bounds alone give a slope over 1.
#[test]
fn a_param_sized_cone_keeps_its_promise() {
    use ink::interval::{I, lipschitz};
    for src in ["(param h 5) (def f (cone-z 3 1 0 h))", "(def h 5) (def f (cone-z 3 1 0 h))", "(def f (cone-z 3 1 0 (ranged 4 6 5)))"] {
        let (d, t) = field(src, "f");
        let ctx = d.ctx(&[], 0.0).unwrap();
        let g = ink::eval::eval_interval_grad(&t, &ctx, [I::new(1.9, 2.1), I::new(-0.1, 0.1), I::new(2.4, 2.6)]);
        assert!(lipschitz(&g) <= 1.0 + 1e-9, "{src}: slope bound {}", lipschitz(&g));
    }
}

/// norm2-up never rounds below the real length: compared exactly, as whole
/// multiples of the smallest float squared (u128), for subnormal and tiny
/// normal lengths, where a widening factor moves nothing (one would give
/// eta for (eta, eta), and inward results at 1e-320 and 1e-310).
#[test]
fn norm2_up_rounds_up_at_every_scale() {
    // x as a whole multiple of 2^-1074 (x < 1e-305: its square fits).
    let units = |x: f64| -> u128 {
        let (bits, x) = (x.to_bits(), x.abs());
        assert!(x < 1e-305);
        let (e, m) = ((bits >> 52) & 0x7ff, (bits & ((1 << 52) - 1)) as u128);
        if e == 0 { m } else { (m | 1 << 52) << (e - 1) }
    };
    let mut r = Lcg(11);
    let mut cases: Vec<(f64, f64)> = vec![(5e-324, 5e-324), (1e-320, 1e-320), (1e-310, 1e-310), (5e-324, 0.0), (3e-323, 4e-323), (2.2e-308, 2.2e-308), (9e-306, 3e-306)];
    for _ in 0..300 {
        let scale = 10f64.powf(r.f(-323.0, -306.0));
        cases.push((scale * r.f(0.0, 1.0), scale * r.f(0.0, 1.0)));
    }
    for (a, b) in cases {
        let (_d, t) = field(&format!("(def f (norm2-up {a:e} {:e}))", -b), "f");
        let n = eval_point(&t, &CTX, [0.0; 3]);
        let (ua, ub, un) = (units(a), units(b), units(n));
        assert!(un * un >= ua * ua + ub * ub, "norm2-up {a:e} {b:e} = {n:e}, below the real length");
        // And tight: within 1e-14 of it, or a few of the smallest floats.
        assert!(n <= a.hypot(b) * (1.0 + 1e-14) + 4.0 * 5e-324, "norm2-up {a:e} {b:e} = {n:e}");
    }
}

/// HTD's teeth meet exactly at their lands: the sector fold mirrors the
/// groove, which is one field read at |x| (the polygon of its mirrored
/// points differed by 2^-52 between its halves: a jump at each land under
/// a promise of 1).
#[test]
fn htd_teeth_meet_exactly_at_their_lands() {
    // The land on +x, points 2e-20 apart across it: any slope moves the
    // field by under 1e-19 there, not by the halves' jump of 2.2e-16.
    let (_d, t) = field("(def f (htd-pulley-teeth-at-tip-2d 30 5 23 0))", "f");
    let mut r = Lcg(13);
    let mut radii: Vec<f64> = (0..400).map(|_| r.f(20.0, 23.5)).collect();
    radii.push(22.126);
    for rad in radii {
        let (fp, fq) = (eval_point(&t, &CTX, [rad, 1e-20, 0.0]), eval_point(&t, &CTX, [rad, -1e-20, 0.0]));
        assert!((fp - fq).abs() <= 1e-19, "at r {rad}: {fp} vs {fq}");
    }
}

/// A load-time reach of a section growing with t holds over the motion: a
/// reach of 1.002 at t = 0 would clip radius 1 + 9 t at t = 1. With no
/// motion stated such profiles are refused, in a sweep too; a section of
/// fixed size loads.
#[test]
fn profile_reach_holds_over_the_motion_or_is_refused() {
    let load = |src: &str| ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).map(|_| ());
    for src in ["(def r (profile-reach (circle-2d (+ 1 (* 9 t)))))", "(def s (sweep (circle-2d (+ 1 (* 9 t))) (vec (vec 0 0 0) (vec 0 0 20)) (vec 1 0 0)))"] {
        let e = load(src).unwrap_err();
        let e = format!("{e:?}");
        assert!(e.contains("profile-reach of a profile that reads t"), "{src}: {e}");
    }
    load("(def s (sweep (circle-2d 10) (vec (vec 0 0 0) (vec 0 0 20)) (vec 1 0 0)))").unwrap();
    // Over a motion stated before it, its reach holds at every moment.
    let d = ink::load_str("(import \"std\")\n(t-span 0 1)\n(def r (profile-reach (circle-2d (+ 1 (* 9 t)))))", "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let Some(ink::lang::Value::Field(r, _)) = d.def("r") else { panic!("no r") };
    let r = d.g.as_const(*r).unwrap();
    assert!((10.0..10.1).contains(&r), "{r}");
}

/// A load-time box of a field that reads t holds at every moment of the
/// file's span, stated before it, and at t = 0, where load reads it: a
/// sphere moving 10 t along x over [0, 1] reaches x = 11. With no span
/// it is refused, and a view at a moment or over a span it does not cover
/// is refused, naming it. A ranged param's box holds over its band, which
/// a run may narrow (--nominal) but not widen.
#[test]
fn bounds_holds_over_the_span_and_the_bands_it_reads() {
    let load = |src: &str| ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT));
    let x = |d: &ink::lang::Doc, n: &str| match d.def(n) {
        Some(ink::lang::Value::Field(f, _)) => d.g.as_const(*f).unwrap(),
        v => panic!("{n}: {v:?}"),
    };
    let ends = "(def lx (nth (nth b 0) 0)) (def hx (nth (nth b 1) 0))";
    let d = load(&format!("(t-span 0 1)\n(def b (bounds (move (sphere 1) (* 10 t) 0 0)))\n{ends}")).unwrap_or_else(|e| panic!("{e:?}"));
    let (lx, hx) = (x(&d, "lx"), x(&d, "hx"));
    assert!((-1.2..=-1.0).contains(&lx) && (11.0..11.2).contains(&hx), "{lx} {hx}");
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(d, ctx);
    v.set_moment(0.5).unwrap();
    v.set_motion(Some((0.0, 1.0))).unwrap();
    let e = v.set_moment(2.0).unwrap_err();
    assert!(e.message.contains("t = 2") && e.message.contains("test:3 (def b ...)"), "{e:?}");
    let e = v.set_motion(Some((0.0, 3.0))).unwrap_err();
    assert!(e.message.contains("test:3 (def b ...)"), "{e:?}");
    let e = v.doc.ctx(&[], 2.0).unwrap_err();
    assert!(e.contains("t = 2"), "{e}");
    let e = format!("{:?}", load("(def b (bounds (move (sphere 1) (* 10 t) 0 0)))").err().unwrap());
    assert!(e.contains("bounds of a field that reads t") && e.contains("t-span"), "{e}");
    load("(def b (bounds (move (sphere 1) 10 0 0)))").unwrap();

    let src = format!("(param px 0 (tol -1 1))\n(def b (bounds (move (sphere 1) px 0 0)))\n{ends}");
    let v = toleranced(&src).unwrap_or_else(|e| panic!("{e:?}"));
    let (lx, hx) = (x(&v.doc, "lx"), x(&v.doc, "hx"));
    assert!((-2.2..=-2.0).contains(&lx) && (2.0..2.2).contains(&hx), "{lx} {hx}");
    let px = v.doc.param_index("px").unwrap();
    v.scope_spans(&[(px, (0.0, 0.0))]).unwrap();
    let e = v.scope_spans(&[(px, (-3.0, 3.0))]).unwrap_err();
    assert!(e.message.contains("'px'") && e.message.contains("test:3 (def b ...)"), "{e:?}");
}

/// A load-time reach holds at every value a run uses: a tolerance band
/// that excludes the param's value (r 10, band 1..2) is reached at both,
/// so the sweep holds (5, 0, 15) at r = 10, and a value beyond what it
/// holds over is refused, naming the form.
#[test]
fn a_load_time_reach_covers_every_value_a_run_uses() {
    let d = ink::load_str("(import \"std\")\n(param r 10 (tol -9 -8))\n(def s (sweep (circle-2d r) (vec (vec 0 0 0) (vec 0 0 20)) (vec 1 0 0)))", "test", Path::new(ROOT))
        .unwrap_or_else(|e| panic!("{e:?}"));
    let Some(Value::Field(s, _)) = d.def("s") else { panic!("no s") };
    let ctx = d.ctx(&[], 0.0).unwrap();
    assert!(eval_point(&Tape::new(&d.g, *s), &ctx, [5.0, 0.0, 15.0]) < 0.0, "inside at r = 10");
    let e = d.ctx(&[("r".to_string(), 20.0.into())], 0.0).unwrap_err();
    assert!(e.contains("'r'") && e.contains("test:3 (def s ...)"), "{e}");
}

/// A load-time box of a field that reads t holds over the motion of the
/// file whose form found it, imported or not: a child's (t-span 0 2) with
/// no span in its parent, and not the parent's shorter one.
#[test]
fn a_load_time_reach_holds_over_its_own_files_motion() {
    let dir = crate::common::Scratch::new("reach-span");
    std::fs::write(dir.join("child.ink"), "(import \"std\")\n(t-span 0 2)\n(def s (sweep (circle-2d (+ 1 t)) (vec (vec 0 0 0) (vec 0 0 20)) (vec 1 0 0)))").unwrap();
    for parent in ["(import \"child.ink\" as c)", "(t-span 0 1)\n(import \"child.ink\" as c)"] {
        let d = ink::load_str(parent, "parent", &dir).unwrap_or_else(|e| panic!("{parent}: {e:?}"));
        d.covers_t(1.5, None).unwrap_or_else(|e| panic!("{parent}: {e:?}"));
        assert!(d.covers_t(2.5, None).is_err(), "{parent}");
    }
}

/// polygon-distance's edges are straight: a rounded corner (polygon's
/// fourth number) or a bulge would be the sharp corner outside and the
/// round one inside, a jump its slope promise forbids (a 10 mm square of
/// 2 mm rounds would hold a clearance of 0.2 to a point 0.14 from it).
/// Refused; a corner of 0s is a corner.
#[test]
fn polygon_distance_refuses_curved_corners() {
    let load = |pts: &str| ink::load_str(&format!("(import \"std\")\n(def f (polygon-distance {pts}))"), "test", Path::new(ROOT)).map(|_| ());
    for pts in [
        "(vec (vec 0 0 0 2) (vec 10 0 0 2) (vec 10 10 0 2) (vec 0 10 0 2))",
        "(vec (vec 0 0) (vec 10 0 0.5) (vec 10 10) (vec 0 10))",
        "(vec (vec 0 0) (vec 10 0) (vec 10 10 0 1) (vec 0 10))",
    ] {
        let e = load(pts).unwrap_err();
        let e = format!("{e:?}");
        assert!(e.contains("polygon-distance takes corners"), "{pts}: {e}");
    }
    load("(vec (vec 0 0 0 0) (vec 10 0) (vec 10 10 0) (vec 0 10))").unwrap();
}

#[test]
fn polygon_sign_holds_on_vertex_rows() {
    use ink::query::Verdict;
    use std::sync::Arc;
    // A T-slot whose vertices share rows y = 4, 7 and 8 with points inside it.
    let slot = "(polygon (vec (vec -3 11) (vec 3 11) (vec 3 8) (vec 6 8) (vec 6 7) (vec 3 4) (vec -3 4) (vec -6 7) (vec -6 8) (vec -3 8)))";
    let (_d, t) = field(&format!("(def f {slot})"), "f");
    for (p, want) in [([0.0, 7.0, 0.0], -3.0), ([0.0, 8.0, 0.0], -3.0), ([-4.0, 7.0, 0.0], -1.0), ([8.0, 7.0, 0.0], 2.0), ([0.0, 4.0, 0.0], 0.0)] {
        let got = eval_point(&t, &CTX, p);
        assert!((got - want).abs() < 1e-9, "at {p:?}: {got} vs {want}");
    }
    // A box straddling the y = 7 row lies inside: proven, not unknown.
    let (_d, probe) = field(&format!("(def outside-part (max (box 1 0.5 1) (neg (move {slot} 0 -7 0))))"), "outside-part");
    assert!(matches!(empty_at(&Arc::new(probe), &CTX, [-0.6, -0.3, -0.6], [0.6, 0.3, 0.6], Some(0.02)).0, Verdict::Proven));
}

#[test]
fn polygon_prunes_to_nearby_edges() {
    use ink::eval::eval_interval;
    use ink::interval::I;
    // A 40-gon: near one edge, interval arithmetic keeps only the nearby
    // edges (a winding number sums all 40 and could not drop any).
    let (_d, t) = field("(def f (polygon (regular-polygon-points 40 10)))", "f");
    let b = [I::new(9.6, 10.1), I::new(-0.3, 0.3), I::new(-1.0, 1.0)];
    let mut ch = vec![ink::tape::BOTH; t.len()];
    eval_interval(&t, &CTX, b, Some(&mut ch));
    let near = t.simplify(&ch);
    assert!(near.len() * 5 < t.len(), "{} of {} instructions near one edge", near.len(), t.len());
    // The pruned tape agrees with the full one there.
    let mut r = Lcg(3);
    for _ in 0..200 {
        let p = [r.f(9.6, 10.1), r.f(-0.3, 0.3), 0.0];
        assert!((eval_point(&near, &CTX, p) - eval_point(&t, &CTX, p)).abs() < 1e-12);
    }
}

#[test]
fn param_expressions_bound_as_points() {
    use ink::query::Verdict;
    use std::sync::Arc;
    // An angle from a param through atan2 exactly on its branch cut (180
    // degrees). Params are fixed for a query, so it is a point, as a def
    // would be; widened by interval rounding it would be the whole circle,
    // and the arm could not be proven away from +x.
    let src = "(param a 0.431)\n(def ang (* (/ 180 pi) (atan2 (- a 0.431) -1)))\n\
               (def f (max (rotate-z (move (box 1 0.4 1) 3 0 0) ang) (move (box 2 2 2) 3 0 0)))";
    let (d, t) = field(src, "f");
    let ctx = d.ctx(&[], 0.0).unwrap();
    assert!(matches!(empty_at(&Arc::new(t), &ctx, [1.5, -1.5, -1.5], [4.5, 1.5, 1.5], Some(0.02)).0, Verdict::Proven));
}

/// Claims atan2's range must not prove: at p = 1, atan2(t, p) = atan t <
/// t - t^3/3 + t^5/5, below each strict lower bound, so neither claim
/// holds over the band, though over p in [-1, 1] a range that took the
/// rule cancelling near the +x axis would prove both. The first
/// fails for p above 0.99999996, and is disproven; the second only within
/// about 1e-14 of p = 1, finer than the search's pieces, so unknown. The
/// true claims beside them stay proven (the bounds are tight to about
/// 1e-15).
#[test]
fn atan2_near_its_axis_proves_no_false_bound() {
    use ink::claims::{ClaimOptions, Status::*, run_asserts};
    let src = "(param p 0 (tol -1 1))
               (assert false-near-axis (value (atan2 1e-10 p) > 1.00000004e-10))
               (assert false-milliscale (value (atan2 0.002 p) > 0.00199999733333975))
               (assert positive (value (atan2 1e-10 p) > 0))
               (assert milliscale (value (atan2 0.002 p) > 0.0019999973))";
    let doc = ink::load_str(src, "claims", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let spans = v.declared_spans().unwrap();
    v.set_param_spans(&spans).unwrap();
    let got = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts;
    let got: Vec<(&str, ink::claims::Status)> = got.iter().map(|a| (a.name.as_str(), a.check.status)).collect();
    assert_eq!(got, vec![("false-near-axis", Fails), ("false-milliscale", Unknown), ("positive", Holds), ("milliscale", Holds)]);
}

/// A (ranged lo hi) is ranged wherever a claim is decided, by the library
/// as by ink assert: a view made with no spans (a pointwise one), or one
/// whose caller narrows it, ranges it too, never answering `holds`,
/// guaranteed, at its lo alone.
#[test]
fn a_ranged_variable_is_ranged_by_every_claim_entry_point() {
    use ink::claims::{ClaimOptions, Status::*, check, run_asserts};
    let view = |src: &str| {
        let doc = ink::load_str(src, "api", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
        let ctx = doc.ctx(&[], 0.0).unwrap();
        ink::view::ModelView::new(doc, ctx)
    };
    let safe = "(def u (ranged 0 1)) (assert safe (value u < 0.5))";
    // check: on a view with no spans, and on one narrowed to [0, 0.25].
    for narrowed in [false, true] {
        let mut v = view(safe);
        if narrowed {
            v.set_param_spans(&[(0, (0.0, 0.25))]).unwrap();
        }
        let claim = v.doc.asserts[0].claim.clone();
        let c = check(&mut v, claim, None, ClaimOptions::default()).unwrap();
        assert!(c.status == Fails && c.ranged == vec![(0, (0.0, 1.0))], "{c:?}");
    }
    // run_asserts: its own reading before check (a bound over the range,
    // here 2u from 2 to 6, which the sphere's 4.19 exceeds at u = 1).
    let mut v = view(safe);
    let a = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts.remove(0);
    assert!(a.check.status == Fails && a.check.ranged == vec![(0, (0.0, 1.0))], "{:?}", a.check);
    let mut v = view("(import \"std\") (def u (ranged 1 3 3)) (assert v (volume (sphere 1) < (* 2 u)))");
    let e = run_asserts(&mut v, None, false, ClaimOptions::default()).map(|r| r.asserts[0].check.status);
    assert!(e.as_ref().is_err_and(|e| e.message.contains("its bound reads ranged param(s) u")), "{e:?}");
}

/// A turn at a moment by an angle on atan2's cut (y exactly 0 left of the
/// origin: +pi, as the one zero reads it) is bounded as tightly as any
/// other, not as [-pi, pi], which would analyse the instance as turned
/// any way: a GT2 path's idler entered from -x (a belt down its left side)
/// would have a box 132 mm across for its 18, and its pairs with the belt
/// and its bracket be unknown after their searches' cell limits. Here a
/// cube turned so, touching a neighbour, which [-pi, pi] would leave
/// unknown after 462 462 undecided cells (1.1 s).
#[test]
fn a_turn_on_atan2s_cut_is_bounded_at_a_moment() {
    use ink::claims::{ClaimOptions, Status::*, run_asserts};
    let src = "(def y0 (+ 5 (sin (* 2 pi t))))
               (body cube (shape (- (max (abs x) (abs y) (abs z)) 1)))
               (assembly root (instance c cube (at (rot-z (* (/ 180 pi) (atan2 (/ 0 y0) -1)))))
                              (instance n cube (at (move 2 0 0))))
               (assert apart (disjoint (shape-of root/c) (shape-of root/n)))";
    let doc = ink::load_str(src, "cut", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let it = v.doc.items("root/c").unwrap().remove(0);
    let (lo, hi, _) = v.item_bounds(&it).unwrap();
    assert!((0..3).all(|k| lo[k] > -1.2 && hi[k] < 1.2), "{lo:?}..{hi:?}");
    let got = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts;
    assert_eq!(got[0].check.status, Holds);
}

#[test]
fn spline_interpolates_its_points() {
    let d = ink::load_str("(import \"std\") (def s (spline-points (vec (vec 0 0 0) (vec 10 5 0) (vec 20 0 5) (vec 30 5 5)) 8))", "t", Path::new(ROOT)).unwrap();
    let Value::Vec(pts) = d.def("s").unwrap().clone() else { panic!() };
    assert_eq!(pts.len(), 3 * 8 + 1);
    let get = |v: &Value| -> Vec<f64> {
        let Value::Vec(c) = v else { panic!() };
        c.iter().map(|x| if let Value::Field(n, _) = x { d.g.as_const(*n).unwrap() } else { panic!() }).collect()
    };
    assert_eq!(get(&pts[8]), vec![10.0, 5.0, 0.0]);
    assert_eq!(get(&pts[16]), vec![20.0, 0.0, 5.0]);
    assert_eq!(get(&pts[24]), vec![30.0, 5.0, 5.0]);
}

fn volume(src: &str, res: f64) -> f64 {
    use std::sync::Arc;
    let mut d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let it = d.items("b").unwrap().remove(0);
    let f = d.placed_field(&it);
    let (c, r) = d.clip_reach(f, &CTX);
    let t = Arc::new(Tape::new(&d.g, c));
    let (blo, bhi, _) = ink::query::bounds(&t, &CTX, r.lo, r.hi).unwrap();
    ink::query::integrate(&t, &[None], false, &CTX, blo, bhi, ink::query::Until::Res(res)).sums[0]
}

#[test]
fn pipe_has_the_right_volume() {
    // Straight pipe: annulus area times length, ends cut flat.
    let v = volume("(body b (shape (pipe (vec (vec 0 0 0) (vec 0 0 50)) 10 2)))", 0.25);
    let exact = PI * (100.0 - 64.0) * 50.0;
    assert!((v - exact).abs() / exact < 5e-3, "{v} vs {exact}");
    // Bent pipe along a sampled quarter circle (radius 40): Pappus.
    let v = volume(
        "(def arc (map (lambda (a) (vec (* 40 (cos a)) (* 40 (sin a)) 0)) (linspace 0 (/ pi 2) 21)))
         (body b (shape (pipe arc 10 2)))",
        0.3,
    );
    let exact = PI * 36.0 * (2.0 * PI * 40.0 / 4.0);
    assert!((v - exact).abs() / exact < 1e-2, "{v} vs {exact}");
}

/// A square swept round a right angle is the mitred bar, its outer corner
/// rounded about the joint at the section's reach: outside the joint,
/// where both segments' nearest point is the joint, each side of the mitre
/// takes its own segment's piece, so neither fills past the mitre.
#[test]
fn sweep_is_mitred_outside_a_joint() {
    let (_d, t) = field("(def s (sweep (rect-2d 2 2) (vec (vec 0 0 0) (vec 10 0 0) (vec 10 10 0)) (vec 0 0 1)))", "s");
    // Past the mitre on the first segment's side, 0.2 from its face, not
    // inside; on the mitre's other side, the same point mirrored.
    for p in [[10.5, -1.2, 0.0], [11.2, -0.5, 0.0]] {
        assert!((eval_point(&t, &CTX, p) - 0.2).abs() < 1e-6, "{p:?}: {}", eval_point(&t, &CTX, p));
    }
    // Inside the outer corner, within the reach (1.416) of the joint.
    assert!((eval_point(&t, &CTX, [10.95, -0.95, 0.0]) + 0.05).abs() < 1e-6);
}

/// A loop (a path whose ends coincide) is swept with its closing joint
/// mitred like the others: round a square, the corner where it closes, at
/// the origin, is the corner at (20, 0) mirrored, for a round section, a
/// flat one (in a plane, so the frames close) and a pipe, not cut flat
/// there: the circle's -0.59 at (-1, -1), outside the corner, not +1.0.
/// The last loop closes only at the param's value, so a select chooses.
#[test]
fn a_loop_is_mitred_where_it_closes() {
    let sq = "(param s 20) (def sq (closed (vec (vec 0 0 0) (vec 20 0 0) (vec 20 20 0) (vec 0 20 0))))
              (def sq2 (vec (vec 0 0 0) (vec 20 0 0) (vec 20 20 0) (vec 0 20 0) (vec (- s 20) 0 0)))";
    for body in ["(sweep (circle-2d 2) sq (vec 0 0 1))", "(sweep (rect-2d 2 3) sq (vec 0 0 1))", "(pipe sq 2 0.5)", "(sweep (circle-2d 2) sq2 (vec 0 0 1))"] {
        let (d, t) = field(&format!("{sq} (def f {body})"), "f");
        let ctx = d.ctx(&[], 0.0).unwrap();
        for p in [[-1.0, -1.0, 0.3], [-0.5, -1.0, 0.0], [-1.0, -0.5, 0.0], [-1.2, 0.5, 1.0], [0.5, -1.2, -1.0], [0.9, 0.7, 0.2]] {
            let (seam, corner) = (eval_point(&t, &ctx, p), eval_point(&t, &ctx, [20.0 - p[0], p[1], p[2]]));
            assert!((seam - corner).abs() < 1e-9, "{body} at {p:?}: {seam} where it closes, {corner} at the next corner");
        }
        if body.starts_with("(sweep (circle") {
            // The distance to the joint less the reach (2, proven to 1/1000
            // of the section's size, so up to 2.004).
            let v = eval_point(&t, &ctx, [-1.0, -1.0, 0.0]);
            assert!((v - (2f64.sqrt() - 2.0)).abs() < 5e-3, "{body}: {v}");
        }
    }
}

/// Round a loop that is not planar, rotation-minimising frames come back
/// turned about the tangent (here 2.8 degrees); path-rmf spreads the turn
/// along the loop and the sweep turns its sections with it from mitre to
/// mitre, so they meet at every joint, the closing one too: the field is
/// continuous, no steeper than 1.1 between points 2e-4 apart near each
/// joint. Unspread, it jumped by 0.097 across the closing seam.
#[test]
fn a_nonplanar_loop_s_sweep_is_continuous() {
    let corners = [[0.0, 0.0, 0.0], [20.0, 0.0, 0.0], [20.0, 20.0, 10.0], [0.0, 20.0, 20.0], [-5.0, 10.0, 5.0]];
    let pts: Vec<String> = corners.iter().map(|p| format!("(vec {} {} {})", p[0], p[1], p[2])).collect();
    let (_d, t) = field(&format!("(def f (sweep (rect-2d 1 4) (closed (vec {})) (vec 0 0 1)))", pts.join(" ")), "f");
    let f = |p: [f64; 3]| eval_point(&t, &CTX, p);
    let (a, b) = ([-0.8159, -1.976, 1.137], [-0.8161, -1.976, 1.137]);
    assert!((f(a) - f(b)).abs() < 2e-4, "across the closing seam: {} and {}", f(a), f(b));
    let mut seed = 0x9e37_79b9_7f4a_7c15u64;
    let mut r = || {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        (seed >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    };
    for c in corners {
        for _ in 0..4000 {
            let p = [c[0] + 3.0 * r(), c[1] + 3.0 * r(), c[2] + 3.0 * r()];
            let d = [r(), r(), r()];
            let k = 2e-4 / (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let q = [p[0] + k * d[0], p[1] + k * d[1], p[2] + k * d[2]];
            assert!((f(p) - f(q)).abs() <= 1.1 * 2e-4, "near {c:?}: {} at {p:?}, {} at {q:?}", f(p), f(q));
        }
    }
}

#[test]
fn sweep_of_a_circle_matches_a_tube() {
    // A circle swept along a gently bent path encloses the same solid as a tube,
    // except near the flat versus round ends; compare mid-path values.
    let src = "(def path (map (lambda (a) (vec (* 60 (cos a)) (* 60 (sin a)) (* 5 a))) (linspace 0 1.2 41)))
               (def swept (sweep (circle-2d 6) path (vec 0 0 1)))
               (def rod (tube path 6))";
    let (_d, ts) = field(src, "swept");
    let (_d2, tt) = field(src, "rod");
    let mut r = Lcg(3);
    let mut worst = 0.0f64;
    for _ in 0..400 {
        let a = r.f(0.2, 1.0);
        let c = [60.0 * a.cos(), 60.0 * a.sin(), 5.0 * a];
        let p = [c[0] + r.f(-8.0, 8.0), c[1] + r.f(-8.0, 8.0), c[2] + r.f(-8.0, 8.0)];
        let (s, t) = (eval_point(&ts, &CTX, p), eval_point(&tt, &CTX, p));
        if t.abs() < 3.0 {
            worst = worst.max((s - t).abs());
        }
    }
    assert!(worst < 0.1, "sweep deviates from tube by {worst}");
}

/// A swept body overlaps itself where its path is thinner than its section
/// reaches: std's sweeps state it (`check`, path-thickness) and `ink check`
/// warns, with where, at the params it is checked at. A path through itself
/// fails; so does a turn whose mitres cross inside the section. A jog's
/// opposite turns cancel, coils that touch pass, and a loop (first point
/// last) has no ends.
#[test]
fn swept_paths_that_overlap_themselves_are_reported() {
    let failures = |src: &str, params: &[(String, ink::lang::Given)]| {
        let d = ink::load_str_with(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT), params).unwrap_or_else(|e| panic!("{e:?}"));
        let ctx = d.ctx(params, 0.0).unwrap();
        let v = ink::view::ModelView::new(d, ctx);
        v.doc.check_failures(&v.ctx, v.accuracy(), &[])
    };
    let none = Vec::<String>::new();
    // An L whose last leg comes back across the first at (5, 0), 2 mm
    // above it at k = 2: clear for a 0.9 mm rod, not for k = 1.
    let l = "(param k 2) (body hose (shape (tube (vec (vec 0 0 0) (vec 10 0 0) (vec 10 10 k) (vec 5 10 k) (vec 5 -5 k)) 0.9)))";
    assert_eq!(failures(l, &[]), none);
    let w = failures(l, &[("k".into(), 1.0.into())]);
    assert!(w.len() == 1 && w[0].contains("runs through itself") && w[0].contains("(5, 0, 0) and (5, 0, 1)") && w[0].contains("k = 1"), "{w:?}");
    // A hairpin: the 2 mm middle segment cannot hold two right-angle
    // mitres of a 1.5 mm section (it needs 3).
    let w = failures("(body pin (shape (sweep (circle-2d 1.5) (vec (vec 0 0 0) (vec 20 0 0) (vec 20 2 0) (vec 0 2 0)) (vec 0 0 1))))", &[]);
    assert!(w.iter().any(|w| w.contains("turns too sharply")), "{w:?}");
    // A 2 mm jog of a 3 mm rod: its turns go opposite ways, and the
    // mitres at its ends are parallel.
    assert_eq!(failures("(body jog (shape (tube (vec (vec 0 0 0) (vec 10 0 0) (vec 10 2 0) (vec 20 2 0)) 3)))", &[]), none);
    // A closed-end spring: coils one wire diameter apart touch, and pass
    // (its 1.3 mm segments, two apart, are near but do not face each other);
    // a pitch 5% tighter overlaps. Swept as a circle, it is judged by the
    // circle's radius, as the tube is.
    let spring = |how: &str, pitch: f64| format!("(body s (shape {}))", how.replace("PTS", &format!("(helix-points 10 {pitch} 4 192)")));
    assert_eq!(failures(&spring("(tube PTS 1)", 2.0), &[]), none);
    assert!(failures(&spring("(tube PTS 1)", 1.9), &[])[0].contains("runs through itself"));
    assert_eq!(failures(&spring("(sweep (circle-2d 1) PTS (vec 0 0 1))", 2.05), &[]), none);
    // A loop: the square ring closes without ends; a 1 mm wide loop of a
    // 1 mm rod runs through itself.
    let ring = |w: f64| format!("(body ring (shape (tube (closed (vec (vec 0 0 0) (vec 20 0 0) (vec 20 {w} 0) (vec 0 {w} 0))) 1)))");
    assert_eq!(failures(&ring(20.0), &[]), none);
    assert!(failures(&ring(1.0), &[]).iter().any(|w| w.contains("runs through itself")));
    // A hose's cover and bore on one path too tight for both: reported for
    // the cover, the larger, which fails wherever the bore does.
    let w = failures("(def u (vec (vec 0 0 0) (vec 10 0 0) (vec 10 1.5 0) (vec 0 1.5 0))) (body hose (shape (subtract (tube u 1) (tube u 0.9))))", &[]);
    assert!(w.len() == 2 && w.iter().all(|w| w.contains("reaches 1 mm")), "{w:?}");
}

/// A path is a loop by one rule: std's path-closes, which the sweeps build
/// their ends by, and the thickness check (check::LOOP_GAP) agree: a
/// 20 mm square whose ends are 1e-10 apart is a loop to both (a check by
/// 1e-12 of the path's largest coordinate would find it open).
#[test]
fn a_loop_is_a_loop_to_the_sweep_and_its_check() {
    use ink::lang::check::{LOOP_GAP, path_thickness};
    let square = |gap: f64| vec![[0.0, 0.0, 0.0], [20.0, 0.0, 0.0], [20.0, 20.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, gap]];
    // Every segment of the loop holds two mitres (bend radius 10); open,
    // the first holds one, so the tightest is found elsewhere.
    let bend_at = |gap: f64| path_thickness(&square(gap), 15.0).0.unwrap().at[0];
    for gap in [0.5 * LOOP_GAP, 2.0 * LOOP_GAP] {
        let src = format!("(import \"std\") (def c (path-closes (vec (vec 0 0 0) (vec 20 0 0) (vec 20 20 0) (vec 0 20 0) (vec 0 0 {gap}))))");
        let d = ink::load_str(&src, "test", Path::new(ROOT)).unwrap();
        let closes =
            d.g.as_const(match d.def("c") {
                Some(Value::Field(f, _)) => *f,
                _ => panic!(),
            })
            .unwrap()
                < 0.0;
        assert_eq!(closes, bend_at(gap) == bend_at(0.0), "gap {gap}: std says loop {closes}");
    }
}

/// profile-reach is the farthest point of a section from its origin,
/// proven, to a thousandth of its size; a section that reads z is taken
/// over z in [-1, 1].
#[test]
fn profile_reach_is_the_farthest_point() {
    for (profile, reach) in [("(circle-2d 1)", 1.0), ("(rect-2d 4 0.5)", 2.0f64.hypot(0.25)), ("(move (circle-2d 1) 3 0 0)", 4.0), ("(circle-2d (+ 1.5 (* 0.5 z)))", 2.0)] {
        let d = ink::load_str(&format!("(import \"std\") (def r (profile-reach {profile}))"), "test", Path::new(ROOT)).unwrap();
        let r =
            d.g.as_const(match d.def("r") {
                Some(Value::Field(f, _)) => *f,
                _ => panic!(),
            })
            .unwrap();
        assert!(r >= reach && r <= reach * 1.003, "{profile}: {r}, not {reach}");
    }
}

/// Signed distance in the plane to a convex polygon (corners in order).
fn convex_dist(p: [f64; 2], c: &[[f64; 2]]) -> f64 {
    let n = c.len();
    let d = (0..n).map(|i| seg_dist([p[0], p[1], 0.0], [c[i][0], c[i][1], 0.0], [c[(i + 1) % n][0], c[(i + 1) % n][1], 0.0])).fold(f64::INFINITY, f64::min);
    let side = |i: usize| {
        let (a, b) = (c[i], c[(i + 1) % n]);
        (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
    };
    if (0..n).all(|i| side(i) >= 0.0) { -d } else { d }
}

/// The shapes the catalog wrote by hand, each against its exact distance:
/// a ring and a hexagon in the plane, a cone and frustums (narrowing and
/// widening, a point at either end) about z, a box by its corners, an
/// extrusion along x, and a row of copies against pattern-linear's.
#[test]
fn std_primitives_are_exact_distances() {
    let mut r = Lcg(7);
    let near = |a: f64, b: f64, at: &dyn std::fmt::Debug| assert!((a - b).abs() < 1e-9, "{at:?}: {a} vs {b}");
    let box_dist = |q: [f64; 3]| q.iter().map(|v| v.max(0.0).powi(2)).sum::<f64>().sqrt() + q.iter().copied().fold(f64::NEG_INFINITY, f64::max).min(0.0);
    let (_d, ring) = field("(def f (extrude-z (ring-2d 3 5) -2 4))", "f");
    let (_d, hex) = field("(def f (hexagon-2d 10))", "f");
    let rh = 10.0 / 3f64.sqrt();
    let hexagon: Vec<[f64; 2]> = (0..6)
        .map(|k| {
            let a = PI / 2.0 + PI / 3.0 * k as f64;
            [rh * a.cos(), rh * a.sin()]
        })
        .collect();
    let cones: Vec<(Tape, [f64; 4])> = [[6.0, 2.0, -3.0, 5.0], [2.0, 7.0, 0.0, 4.0], [5.0, 0.0, 1.0, 9.0], [0.0, 4.0, -2.0, 2.0], [3.0, 3.0, 0.0, 10.0]]
        .into_iter()
        .map(|[r0, r1, z0, z1]| (field(&format!("(def f (cone-z {r0} {r1} {z0} {z1}))"), "f").1, [r0, r1, z0, z1]))
        .collect();
    let (_d, corners) = field("(def f (emax (slab-x -1 3) (slab-y 2 5) (slab-z -4 0)))", "f");
    let (_d, along_x) = field("(def f (extrude-x (rect-2d 4 2) 1 6))", "f");
    for _ in 0..4000 {
        let p = [r.f(-12.0, 12.0), r.f(-12.0, 12.0), r.f(-8.0, 14.0)];
        let rad = p[0].hypot(p[1]);
        let (ring_2d, dz) = ((rad - 4.0).abs() - 1.0, (p[2] - 1.0).abs() - 3.0);
        let want = if ring_2d > 0.0 && dz > 0.0 { ring_2d.hypot(dz) } else { ring_2d.max(dz) };
        near(eval_point(&ring, &CTX, p), want, &p);
        near(eval_point(&hex, &CTX, p), convex_dist([p[0], p[1]], &hexagon), &p);
        for (t, [r0, r1, z0, z1]) in &cones {
            let profile = [[-r0, *z0], [*r0, *z0], [*r1, *z1], [-r1, *z1]];
            near(eval_point(t, &CTX, p), convex_dist([rad, p[2]], &profile), &(p, r0, r1));
        }
        near(eval_point(&corners, &CTX, p), box_dist([(p[0] - 1.0).abs() - 2.0, (p[1] - 3.5).abs() - 1.5, (p[2] + 2.0).abs() - 2.0]), &p);
        near(eval_point(&along_x, &CTX, p), box_dist([(p[0] - 3.5).abs() - 2.5, p[1].abs() - 2.0, p[2].abs() - 1.0]), &p);
    }
    // A row of symmetric copies is pattern-linear's field, past its ends too.
    let (_d, row) = field("(def f (repeat-linear (cylinder-z 1 2) 5 3 1 0))", "f");
    let (_d, copies) = field("(def f (pattern-linear (cylinder-z 1 2) 5 3 1 0))", "f");
    for _ in 0..2000 {
        let p = [r.f(-6.0, 20.0), r.f(-6.0, 9.0), r.f(-3.0, 3.0)];
        near(eval_point(&row, &CTX, p), eval_point(&copies, &CTX, p), &p);
    }
}

#[test]
fn patterns_make_copies() {
    let v1 = volume("(body b (shape (move (sphere 3) 20 0 0)))", 0.2);
    let v6 = volume("(body b (shape (pattern-polar (move (sphere 3) 20 0 0) 6)))", 0.2);
    assert!((v6 / v1 - 6.0).abs() < 0.02, "{}", v6 / v1);
}

#[test]
fn swept_rectangle_has_area_times_length() {
    // A 6 x 3 rectangle along a helix (radius 30, pitch 20, one turn): volume = A L
    // when the profile centroid follows the path and curvature is gentle.
    let v = volume(
        "(def path (helix-points 30 20 1 45))
         (body b (shape (sweep (rect-2d 6 3) path (vec 0 0 1))))",
        0.4,
    );
    let len = (2.0 * PI * 30.0f64).hypot(20.0);
    let exact = 18.0 * len;
    assert!((v - exact).abs() / exact < 1e-2, "{v} vs {exact}");
}

#[test]
fn provenance_names_the_construction_under_a_point() {
    let d = ink::load_str(
        "(import \"std\")
         (def left (sphere 10))
         (def right (move (box 10 10 10) 15 0 0))
         (def both (min left (tag hole-free right)))",
        "t",
        Path::new(ROOT),
    )
    .unwrap();
    let Value::Field(f, _) = d.def("both").unwrap().clone() else { panic!() };
    let t = Tape::new(&d.g, f);
    let names = |p: [f64; 3]| -> Vec<String> { ink::eval::provenance(&t, &CTX, p).iter().map(|&i| d.g.tag_names[i as usize].clone()).collect() };
    assert_eq!(names([-10.0, 0.0, 0.0]), vec!["both", "left"]);
    assert_eq!(names([20.0, 0.0, 0.0]), vec!["both", "hole-free", "right"]);
}

#[test]
fn fillet_union_rounds_the_inside_edge() {
    // Solids z > 0 and x > 0 (fields -z and -x) joined with r = 2: the empty
    // quadrant's corner gets a quarter-circle fillet centred at (-2, 0, -2).
    let (_d, t) = field("(def f (fillet-union (- z) (- x) 2 90))", "f");
    for a in [0.2f64, 0.7, 1.2] {
        let p = [-2.0 + 2.0 * a.cos(), 0.0, -2.0 + 2.0 * a.sin()];
        assert!(eval_point(&t, &CTX, p).abs() < 1e-9, "not on the fillet at {p:?}");
    }
    // Away from the edge the field is unchanged: distance 3 to the nearer face.
    assert!((eval_point(&t, &CTX, [-5.0, 0.0, -3.0]) - 3.0).abs() < 1e-9);
    assert!(eval_point(&t, &CTX, [5.0, 0.0, 1.0]) < 0.0);
}

/// Signed distance from (s, t) to the convex wedge {n_a . q >= k, n_b . q >= k}
/// (unit normals; k > 0 for a closing's ball centres, the wedge's own apex
/// solved), by its faces and apex: the geometry, not std's formula.
fn wedge_distance(p: [f64; 2], na: [f64; 2], nb: [f64; 2], k: f64) -> f64 {
    let dot = |a: [f64; 2], b: [f64; 2]| a[0] * b[0] + a[1] * b[1];
    let (ga, gb) = (dot(na, p) - k, dot(nb, p) - k);
    if ga >= 0.0 && gb >= 0.0 {
        return -ga.min(gb);
    }
    let det = na[0] * nb[1] - na[1] * nb[0];
    let q = [k * (nb[1] - na[1]) / det, k * (na[0] - nb[0]) / det];
    let mut d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt();
    for (n, g, m) in [(na, ga, nb), (nb, gb, na)] {
        let foot = [p[0] - g * n[0], p[1] - g * n[1]];
        if g < 0.0 && dot(m, foot) - k >= 0.0 {
            d = d.min(-g);
        }
    }
    d
}

#[test]
fn fillets_roll_a_ball_at_the_stated_angle() {
    // Two half-spaces through the y axis, a = z and b leaning so the air
    // between them (fillet-union) or the material (fillet-intersect) opens
    // `angle` degrees. The union's closing by a ball of radius r is the set
    // within r of no ball centre {a >= r, b >= r}: field r - (distance to
    // that wedge); the intersection's opening is within r of the centres
    // {a <= -r, b <= -r}: (distance to that wedge) - r. Checked inside and
    // out, near the edge and away from it.
    let r = 1.5;
    let mut rng = Lcg(5);
    for angle in [20.0f64, 45.0, 60.0, 90.0, 120.0, 150.0, 170.0] {
        let phi = (180.0 - angle).to_radians();
        let (na, nb) = ([0.0, 1.0], [phi.sin(), phi.cos()]); // (x, z)
        let src = format!(
            "(def b (+ (* {} x) (* {} z)))
             (def u (fillet-union z b {r} {angle}))
             (def i (fillet-intersect z b {r} {angle}))",
            nb[0], nb[1]
        );
        let (_d, tu) = field(&src, "u");
        let (_d, ti) = field(&src, "i");
        for _ in 0..400 {
            let p = [rng.f(-8.0, 8.0), rng.f(-3.0, 3.0), rng.f(-8.0, 8.0)];
            let q = [p[0], p[2]];
            let want_u = r - wedge_distance(q, na, nb, r);
            let want_i = wedge_distance([-q[0], -q[1]], na, nb, r) - r;
            let (got_u, got_i) = (eval_point(&tu, &CTX, p), eval_point(&ti, &CTX, p));
            assert!((got_u - want_u).abs() < 1e-9, "union at {angle}: {got_u} vs {want_u} at {p:?}");
            assert!((got_i - want_i).abs() < 1e-9, "intersect at {angle}: {got_i} vs {want_i} at {p:?}");
        }
        // Interval bounds prove a box 0.004 across, 0.02 off the rolled
        // surface on either side: the ball's centre is on the air's
        // bisector, r / sin(angle / 2) from the edge, and a point r + d from
        // it toward the edge reads -d.
        let bis = [na[0] + nb[0], na[1] + nb[1]];
        let len = (bis[0] * bis[0] + bis[1] * bis[1]).sqrt();
        let c = r / (angle.to_radians() * 0.5).sin();
        let centre = [bis[0] / len * c, bis[1] / len * c];
        for d in [0.02, -0.02] {
            let at = [centre[0] - bis[0] / len * (r + d), centre[1] - bis[1] / len * (r + d)];
            let bx = [at[0], 0.0, at[1]].map(|v| ink::interval::I::new(v - 0.002, v + 0.002));
            let bound = ink::eval::eval_interval(&tu, &CTX, bx, None);
            assert!(if d < 0.0 { bound.lo > 0.0 } else { bound.hi < 0.0 }, "{angle}: {bound:?} at {at:?}");
            // ... and enclose every value in their box.
            for _ in 0..20 {
                let p = bx.map(|i| rng.f(i.lo, i.hi));
                let v = eval_point(&tu, &CTX, p);
                assert!(bound.lo <= v && v <= bound.hi, "{angle}: {v} outside {bound:?}");
            }
        }
    }
    // A boss standing on a floor (the round a torus about its axis) and a
    // cone leaning 20 degrees in (110 degrees in the air): each section
    // through the axis is the planar corner at radius 5.
    for (boss, angle) in [("(move (cylinder-z 5 60) 0 0 10)", 90.0f64), ("(* (cos (deg->rad 20)) (- (len2 x y) (- 5 (* z (tan (deg->rad 20))))))", 110.0)] {
        let (_d, t) = field(&format!("(def f (fillet-union z {boss} 2 {angle}))"), "f");
        let phi = (180.0 - angle).to_radians();
        let (na, nb) = ([0.0, 1.0], [phi.sin(), phi.cos()]);
        for _ in 0..400 {
            let (rho, th, z) = (rng.f(3.0, 12.0), rng.f(0.0, 6.3), rng.f(-3.0, 6.0));
            let want = 2.0 - wedge_distance([rho - 5.0, z], na, nb, 2.0);
            let got = eval_point(&t, &CTX, [rho * th.cos(), rho * th.sin(), z]);
            assert!((got - want).abs() < 1e-9, "{boss}: {got} vs {want} at rho {rho} z {z}");
        }
    }
    // At 90 degrees, std's square fillet's formula, value for value, except
    // where a ball holds the point (that formula's beads: under the
    // coplanar bottoms of a and the boss, say), which is min(a, b) (max).
    let src = "(def a (move (box 20 20 4) 0 0 2))
               (def b (min (move (cylinder-z 3 20) 4 0 10) (move (box 2 15 12) -5 0 8)))
               (def c (max (circle-2d 8) (- (abs z) 5)))
               (def new-u (fillet-union a b 1.5 90))
               (def old-u (select (neg (min (- 1.5 a) (- 1.5 b)))
                 (- (max 1.5 (min a b)) (len2 (max (- 1.5 a) 0) (max (- 1.5 b) 0))) (min a b)))
               (def new-i (fillet-intersect c (- z 3) 1.2 90))
               (def old-i (select (neg (min (+ 1.2 c) (+ 1.2 (- z 3))))
                 (+ (min -1.2 (max c (- z 3))) (len2 (max (+ 1.2 c) 0) (max (+ 1.2 (- z 3)) 0))) (max c (- z 3))))";
    // Interval bounds: no looser than the hull of that formula's and
    // min(a, b)'s, where it cannot tell whether a ball holds the box.
    let src = format!("{src} (def sharp-u (min a b)) (def sharp-i (max c (- z 3)))");
    let mut beads = 0;
    for (new, old, sharp) in [("new-u", "old-u", "sharp-u"), ("new-i", "old-i", "sharp-i")] {
        let [tn, to, ts] = [new, old, sharp].map(|n| field(&src, n).1);
        for _ in 0..2000 {
            let p = [rng.f(-12.0, 12.0), rng.f(-12.0, 12.0), rng.f(-4.0, 16.0)];
            let (vn, vo) = (eval_point(&tn, &CTX, p), eval_point(&to, &CTX, p));
            if vn != vo {
                assert_eq!(vn, eval_point(&ts, &CTX, p), "{new} at {p:?}: {vn} vs {vo}");
                beads += 1;
            }
            let bx = p.map(|v| ink::interval::I::new(v - 0.3, v + 0.3));
            let [bn, bo, bs] = [&tn, &to, &ts].map(|t| ink::eval::eval_interval(t, &CTX, bx, None));
            assert!(bn.lo >= bo.lo.min(bs.lo) && bn.hi <= bo.hi.max(bs.hi), "{new}: {bn:?} looser than {bo:?} and {bs:?} over {bx:?}");
        }
    }
    assert!(beads > 0, "no bead removed");
}

#[test]
fn fillets_and_chamfers_add_nothing_where_the_faces_are_not_a_corner() {
    // A lug buried 4 into a 5 mm base: under the base, the lug's buried
    // faces and the base's bottom are both within reach, and a formula
    // reading only their values would take them for a square corner: a
    // bead 1.35 deep under the base, outside both. A ball
    // touching the base's bottom straight above it is clear of the lug, so
    // nothing is filled. Where the lug meets the base's top, the fillet is
    // the 90 degree closing of the two planes, value for value.
    let src = "(def base (move (box 90 20 5) 0 0 2.5))
               (def lug (move (box 24 6 29) 0 0 15.5))
               (def f (fillet-union lug base 5 90))
               (def ch (chamfer-union lug base 5))
               (def fi (fillet-intersect (neg lug) (neg base) 5 90))
               (def ci (chamfer-intersect (neg lug) (neg base) 5))";
    let [tf, tc, tfi, tci] = ["f", "ch", "fi", "ci"].map(|n| field(src, n).1);
    let mut rng = Lcg(11);
    for _ in 0..20000 {
        // half of them close under the lug's buried edges and corners
        let p = if rng.f(0.0, 1.0) < 0.5 { [rng.f(-30.0, 30.0), rng.f(-12.0, 12.0), rng.f(-4.0, -0.0001)] } else { [rng.f(-17.0, 17.0), rng.f(-8.0, 8.0), rng.f(-0.1, -0.0001)] };
        for (t, what) in [(&tf, "fillet"), (&tc, "chamfer")] {
            assert!(eval_point(t, &CTX, p) > 0.0, "{what} adds material under the base at {p:?}");
        }
        for (t, what) in [(&tfi, "fillet-intersect"), (&tci, "chamfer-intersect")] {
            assert!(eval_point(t, &CTX, p) < 0.0, "{what} at {p:?}");
        }
    }
    let (r, na, nb) = (5.0, [1.0, 0.0], [0.0, 1.0]);
    for _ in 0..400 {
        let (x, y, z) = (rng.f(-10.0, 10.0), rng.f(3.0, 9.0), rng.f(5.0, 9.0));
        let want = r - wedge_distance([y - 3.0, z - 5.0], na, nb, r);
        let got = eval_point(&tf, &CTX, [x, y, z]);
        assert!((got - want).abs() < 1e-9, "fillet {got} vs {want} at {:?}", [x, y, z]);
    }
    // The chamfer at the square corner: the 45 degree plane 5 along each face.
    let got = eval_point(&tc, &CTX, [0.0, 5.5, 7.5]);
    assert!(got.abs() < 1e-9, "chamfer {got} on its plane");
}

#[test]
fn shapes_placed_on_a_path_frame() {
    // A small cube placed at the second segment's frame of a straight path
    // along +y sits at its start, with its local z along the path.
    let src = "(def path (vec (vec 0 0 0) (vec 0 10 0) (vec 0 20 0)))
               (def f (place (box 2 2 4) (nth (path-frames path (vec 0 0 1)) 1)))";
    let (_d, t) = field(src, "f");
    assert!(eval_point(&t, &CTX, [0.0, 10.0, 0.0]) < -0.9);
    assert!(eval_point(&t, &CTX, [0.0, 11.9, 0.0]) < 0.0 && eval_point(&t, &CTX, [0.0, 12.1, 0.0]) > 0.0, "local z should run along +y");
}

#[test]
fn solve_finds_constrained_values_and_follows_time() {
    // Circle-line intersection: x^2 + y^2 = 25, y = x - 1 -> (4, 3) from a nearby guess.
    let d = ink::load_str(
        "(solve meet (unknowns (px 3) (py 2)) (== (+ (square px) (square py)) 25) (== py (- px 1)))
         (def total (+ px py))",
        "t",
        Path::new(ROOT),
    )
    .unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let i = |n: &str| ctx.params[d.param_index(n).unwrap()];
    assert!((i("px") - 4.0).abs() < 1e-9 && (i("py") - 3.0).abs() < 1e-9, "({}, {})", i("px"), i("py"));
    // Time-dependent: a point on a unit circle tracking angle 2 pi t.
    let d = ink::load_str("(solve track (unknowns (a 0)) (== (sin a) (sin (* 2 pi t))) )", "t", Path::new(ROOT)).unwrap();
    let c = d.ctx(&[], 0.1).unwrap();
    assert!((c.params[0].sin() - (0.2 * std::f64::consts::PI).sin()).abs() < 1e-9);
    // Solved params cannot be overridden directly.
    assert!(d.ctx(&[("a".into(), 1.0.into())], 0.0).is_err());
}

#[test]
fn extension_forms_go_to_their_packages() {
    // solve and optimize are read by their registered package and kept in
    // document order; typed access reads them back. (With a study: ext's
    // sim tests.)
    use ink::design::{Optimize, SolveSys};
    let src = "(param w 10)
(body bar (shape (box 100 w 10)))
(solve fit (unknowns (a 1)) (== (* a a) 4))
(optimize light (vary w 5 20) (minimize (volume (shape-of bar))))";
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let kept: Vec<(&str, &str)> = d.ext.iter().map(|e| (e.package, e.head)).collect();
    assert_eq!(kept, [("design", "solve"), ("design", "optimize")]);
    assert!(d.form::<SolveSys>("fit").is_some() && d.form::<SolveSys>("bar").is_none());
    assert!(d.form::<Optimize>("light").is_some() && d.form::<SolveSys>("light").is_none());
    let fit = d.form::<SolveSys>("fit").unwrap();
    assert!((d.params[fit.unknowns[0]].default - 2.0).abs() < 1e-9);
    // A solve is a param provider: the core applies it to every context
    // and refuses overrides of what it owns.
    assert_eq!(d.providers().owned().collect::<Vec<_>>(), fit.unknowns);
    let e = d.ctx(&[("a".into(), 3.0.into())], 0.0).err().unwrap();
    assert!(e.contains("'a' is solved by 'fit'"), "{e}");
    assert!((d.ctx(&[("w".into(), 12.0.into())], 0.0).unwrap().params[fit.unknowns[0]] - 2.0).abs() < 1e-9);
    // Unknown heads list every package's forms; the kernel's has no study.
    let e = ink::load_str("(solv x)", "test", Path::new(ROOT)).err().unwrap();
    assert!(e.hint.as_deref().is_some_and(|h| h.contains("'solve'")), "{e:?}");
    let e = ink::load_str("(study s (shape (box 1 1 1)))", "test", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("unknown top-level form 'study'"), "{e:?}");
    // A param is varied once.
    let e =
        ink::load_str(&format!("(import \"std\")\n{src}\n(optimize twice (vary w 5 20) (vary w 6 9) (minimize (volume (shape-of bar))))"), "test", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("'w' is varied twice"), "{e:?}");
    // A document reads the forms of the packages it is loaded with: the
    // core alone (and std, which is plain core) knows none of these.
    static CORE: ink::lang::Registry = ink::lang::Registry::new(&[]);
    let e = Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "test", Path::new(ROOT), &[]).err().unwrap();
    assert!(e.message.contains("unknown top-level form 'solve'"), "{e:?}");
}

#[test]
fn packages_bring_their_own_metrics() {
    // A package registers forms and the metrics they give; claims compare a
    // metric of a form with a bound and objectives optimise it, all read by
    // one reader. Here a gauge reads an expression of params.
    use ink::design::{Objective, Optimize};
    use ink::lang::{Claim, Form, Measured, Outcome, Package, Reader, Registry, RunOpts};
    use ink::sexp::Sexp;
    use std::sync::Arc;
    struct Gauge;
    struct GaugeForm(String, u32);
    struct Reading(f64);
    impl Package for Gauge {
        fn name(&self) -> &'static str {
            "gauge"
        }
        fn forms(&self) -> &'static [&'static str] {
            &["gauge"]
        }
        fn metrics(&self) -> &'static [(&'static str, &'static str)] {
            &[("reading", "the value a gauge reads")]
        }
        fn read(&self, r: &mut Reader, form: &Sexp, _head: &str, args: &[Sexp]) -> Result<Arc<dyn Form>, ink::InkError> {
            let [name, e] = args else { return Err(r.err(form, "usage: (gauge name expr)")) };
            Ok(Arc::new(GaugeForm(r.name(name)?, r.field(e)?)))
        }
    }
    impl Form for GaugeForm {
        fn name(&self) -> &str {
            &self.0
        }
        fn metrics(&self) -> Vec<&'static str> {
            vec!["reading"]
        }
        fn run(&self, doc: &mut Doc, ctx: &Ctx, _opts: &RunOpts) -> Result<Box<dyn Outcome>, ink::InkError> {
            Ok(Box::new(Reading(doc.scalar(self.1, ctx))))
        }
    }
    impl Outcome for Reading {
        fn field(&self, _: &str) -> Option<ink::grid::GridData> {
            None
        }
        fn metric(&self, m: &str) -> Option<Measured> {
            (m == "reading").then(|| Measured {
                value: Ok(self.0),
                at: None,
                res: None,
                warnings: Vec::new(),
                excluded: Vec::new(),
                unseen: None,
                proven: false,
                held: None,
                band: None,
            })
        }
        fn failure(&self) -> Option<String> {
            (self.0 < 0.0).then(|| "the gauge reads below its range".into())
        }
    }
    static GAUGES: Registry = Registry::new(&[&Gauge, &ink::design::Design]);
    let src = "(param p 3)\n(gauge g (* 2 p))\n(assert low (reading g < 10))\n\
               (optimize most (vary p 1 4) (maximize (reading g)) (require (reading g <= 7)))";
    let mut d = Doc::load_str_in(&GAUGES, src, "test", Path::new(ROOT), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Claim::Metric { metric, value, .. } = &d.asserts[0].claim else { panic!("not a metric claim") };
    assert_eq!((metric.package, metric.form.as_str(), metric.metric, *value), ("gauge", "g", "reading", 10.0));
    let o = d.form::<Optimize>("most").unwrap();
    assert!(matches!(&o.objective, Objective::Metric(m) if m.metric == "reading") && o.requires.len() == 1);
    // Running the form measures its metric.
    let (f, ctx) = (d.ext_form(Some("gauge"), "g").unwrap().form.clone(), d.ctx(&[], 0.0).unwrap());
    assert_eq!(f.run(&mut d, &ctx, &RunOpts::default()).unwrap().metric("reading").map(|m| m.value), Some(Ok(6.0)));
    // A run that failed gives no metric: a claim on it is unknown, with the
    // reason, whatever the package (the one measure path).
    let broken = Doc::load_str_in(&GAUGES, "(param p -1)\n(gauge g (* 2 p))\n(assert low (reading g < 10))", "test", Path::new(ROOT), &[]).unwrap();
    let ctx = broken.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(broken, ctx);
    let c = ink::claims::run_asserts(&mut v, None, false, ink::claims::ClaimOptions::default()).unwrap().asserts.remove(0).check;
    assert!(
        c.status == ink::claims::Status::Unknown
            && matches!(&c.evidence, ink::claims::Evidence::Metric { measured, .. } if measured.value.as_ref().is_err_and(|w| w.contains("below its range"))),
        "{c:?}"
    );
    // The registry holds through imports, namespaced ones included.
    let dir = crate::common::Scratch::new("gauge");
    std::fs::write(dir.join("gauges.ink"), "(param p 3)\n(gauge g (* 2 p))\n").unwrap();
    let main = "(import \"gauges.ink\" as m)\n(assert low (reading m.g < 10))";
    let d2 = Doc::load_str_in(&GAUGES, main, "test", &dir, &[]).unwrap_or_else(|e| panic!("{e:?}"));
    assert!(matches!(&d2.asserts[0].claim, Claim::Metric { metric, .. } if metric.form == "m.g"));
    // The design loop measures it too, each design a reload of the view
    // with the view's own packages.
    std::fs::write(dir.join("most.ink"), src).unwrap();
    let v = ink::view::ModelView::load_in(&GAUGES, &dir.join("most.ink"), &[], 0.0).unwrap();
    let o = ink::design::run_optimize(&v, "most", None, &ink::lang::command::Args::parse(&["optimize".into()])).unwrap_or_else(|e| panic!("{e:?}"));
    assert!(o.feasible && (o.design[0] - 3.5).abs() < 0.05 && (o.objective - 7.0).abs() < 0.1, "{:?} {}", o.design, o.objective);
    assert!(o.requires.len() == 1 && (o.requires[0].measured - o.objective).abs() < 1e-12);
    // The standard distribution has no gauge package.
    let e = ink::load_str(src, "test", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("unknown top-level form 'gauge'"), "{e:?}");
}

#[test]
fn a_banded_metric_holds_only_where_its_whole_band_does() {
    // A gauge reads 2p on cells of 1, measured to lie within 1 of it (its
    // band, "placement"), and 2p + (h - 1) / 2 on cells of h: its ladder
    // (6, 6.5, 7.5 on cells of 1, 2, 4) gives a grid error of 0.625, so the
    // claim's band is 4.375 to 7.625. A claim holds where both ends hold,
    // fails where both fail, and is otherwise unknown, naming what dominates,
    // whatever its value alone says: the grid's error makes both a claim
    // that would hold and one that would fail on the band alone unknown.
    use ink::lang::{Band, Form, Measured, Outcome, Package, Reader, Registry, RunOpts};
    use ink::sexp::Sexp;
    use std::sync::Arc;
    struct Gauge;
    struct GaugeForm(String, u32);
    struct Reading(f64, f64);
    impl Package for Gauge {
        fn name(&self) -> &'static str {
            "gauge"
        }
        fn forms(&self) -> &'static [&'static str] {
            &["gauge"]
        }
        fn metrics(&self) -> &'static [(&'static str, &'static str)] {
            &[("reading", "the value a gauge reads")]
        }
        fn read(&self, r: &mut Reader, form: &Sexp, _head: &str, args: &[Sexp]) -> Result<Arc<dyn Form>, ink::InkError> {
            let [name, e] = args else { return Err(r.err(form, "usage: (gauge name expr)")) };
            Ok(Arc::new(GaugeForm(r.name(name)?, r.field(e)?)))
        }
    }
    impl Form for GaugeForm {
        fn name(&self) -> &str {
            &self.0
        }
        fn metrics(&self) -> Vec<&'static str> {
            vec!["reading"]
        }
        fn run(&self, doc: &mut Doc, ctx: &Ctx, opts: &RunOpts) -> Result<Box<dyn Outcome>, ink::InkError> {
            let h = opts.res.unwrap_or(1.0);
            Ok(Box::new(Reading(doc.scalar(self.1, ctx) + (h - 1.0) / 2.0, h)))
        }
    }
    impl Outcome for Reading {
        fn metric(&self, m: &str) -> Option<Measured> {
            let band = Band { low: self.0 - 1.0, high: self.0 + 1.0, parts: vec![("placement".into(), 1.0)], basis: "a test's".into() };
            (m == "reading").then(|| Measured {
                value: Ok(self.0),
                at: None,
                res: Some(self.1),
                warnings: Vec::new(),
                excluded: Vec::new(),
                unseen: None,
                proven: false,
                held: None,
                band: Some(band),
            })
        }
    }
    static GAUGES: Registry = Registry::new(&[&Gauge]);
    let src = "(param p 3)\n(gauge g (* 2 p))\n(assert under (reading g < 8))\n(assert near (reading g < 7.5))\n(assert over (reading g > 8))\n(assert above (reading g > 7.5))";
    let d = Doc::load_str_in(&GAUGES, src, "test", Path::new(ROOT), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(d, ctx);
    let answers = ink::claims::run_asserts(&mut v, None, false, ink::claims::ClaimOptions::default()).unwrap().asserts;
    let status: Vec<_> = answers.iter().map(|a| a.check.status).collect();
    use ink::claims::Status::*;
    assert_eq!(status, [Holds, Unknown, Fails, Unknown], "{answers:?}");
    let near = &answers[1].check;
    assert!(
        matches!(&near.evidence, ink::claims::Evidence::Metric { undecided: Some((why, _)), grid_error: Some((e, _)), .. }
            if why.contains("4.375 to 7.625") && why.contains("placement dominates") && (e - 0.625).abs() < 1e-12),
        "{near:?}"
    );
}

#[test]
fn volume_derivatives_match_geometry() {
    // A plate a x b x c with a bore of radius r, moved by d:
    // V = a b c - pi r^2 c, so dV/dd = 0 and the rest follow.
    use std::sync::Arc;
    let src = "(import \"std\")\n(param a 40) (param b 30) (param c 6) (param r 5) (param d 7)
(body p (shape (move (subtract (box a b c) (cylinder-z r 20)) d 0 0)))";
    let mut doc = ink::load_str(src, "test", Path::new(ROOT)).unwrap();
    let it = doc.items("p").unwrap().remove(0);
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let f = doc.placed_field(&it);
    let f = doc.clip_reach(f, &ctx).0;
    let t = Arc::new(Tape::new(&doc.g, f));
    let (area, d) = ink::design::shape_derivatives(&t, &[&|_| 1.0], &ctx, [-20.0, -20.0, -10.0], [40.0, 20.0, 10.0], 0.25);
    let (a, b, c, r) = (40.0, 30.0, 6.0, 5.0);
    let want = [b * c, a * c, a * b - PI * r * r, -2.0 * PI * r * c, 0.0];
    for (k, w) in want.iter().enumerate() {
        let got = d[0][k];
        assert!((got - w).abs() < 0.02 * (a * b), "param {k}: {got} vs {w}");
    }
    assert!((area - (2.0 * (a * b - PI * r * r) + 2.0 * (a + b) * c + 2.0 * PI * r * c)).abs() < 0.02 * area, "area {area}");
}

#[test]
fn redistance_makes_offsets_exact() {
    // A box whose field is scaled differently per axis: offsetting it
    // directly moves each face by a different amount; after redistance the
    // offset is a true rounded box.
    let skewed = "(def skewed (max (* 3 (- (abs x) 10)) (* 0.5 (- (abs y) 10)) (- (abs z) 10)))";
    let v = volume(&format!("{skewed}\n(body b (shape (- (redistance skewed (res 0.25)) 1)))"), 0.2);
    let want = 8000.0 + 2400.0 + 3.0 * PI * 20.0 + 4.0 / 3.0 * PI;
    assert!((v - want).abs() < 0.002 * want, "{v} vs {want}");
}

#[test]
fn iso_thread_has_the_profile_volume() {
    // M6 x 1 over 10 mm: 10 * mean over one pitch of pi R(u)^2, where R is the
    // ISO profile (crest flat p/8, root flat p/4, 60 degree flanks).
    let v = volume("(body b (shape (max (thread-iso 6 1) (slab-z 0 10))))", 0.04);
    let (p, depth): (f64, f64) = (1.0, 0.5412659);
    let n = 20000;
    let want: f64 = (0..n)
        .map(|i| {
            let u = (i as f64 + 0.5) / n as f64 * p;
            let w = (u - 0.5 * p).abs();
            let r = 3.0 - depth + depth.min(((0.375 * p - w) * 3f64.sqrt()).max(0.0));
            PI * r * r
        })
        .sum::<f64>()
        * 10.0
        / n as f64;
    assert!((v - want).abs() < 1e-3 * want, "{v} vs {want}");
}

#[test]
fn char_codes_are_code_points() {
    // char-codes: code points, one per character (not per byte).
    let (_d, t) = field("(def f (+ (len (char-codes \"h°llo\")) (nth (char-codes \"A°\") 1)))", "f");
    assert_eq!(eval_point(&t, &CTX, [0.0; 3]), 5.0 + 176.0);
}

#[test]
fn involute_gear_has_the_right_teeth() {
    // z = 20, module 2: pitch radius 20. On the pitch circle the tooth
    // (without backlash) spans half the circumference; outside the addendum
    // (22) and inside the root (17.5) the gear is empty and solid.
    let (_, t) = field("(def g (involute-gear-2d 20 2 20 0))", "g");
    let n = 3600;
    let pts: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let a = 2.0 * PI * i as f64 / n as f64;
            [20.0 * a.cos(), 20.0 * a.sin(), 0.0]
        })
        .collect();
    let mut v = vec![0.0; n];
    eval_points(&t, &CTX, &pts, &mut v);
    let inside = v.iter().filter(|x| **x < 0.0).count() as f64 / n as f64;
    assert!((inside - 0.5).abs() < 0.005, "tooth fraction on the pitch circle {inside}");
    let mut w = vec![0.0; 2];
    eval_points(&t, &CTX, &[[22.1, 0.0, 0.0], [0.0, 17.4, 0.0]], &mut w);
    assert!(w[0] > 0.0 && w[1] < 0.0, "{w:?}");
    // The flank is one expression of r: a series giving way to
    // tn - atan(tn) at tn = 0.9 would put a step of 0.11 there, 0.057
    // outside and inside 2e-6 mm apart (a 6-tooth gear's flank).
    let (_, t) = field("(def g (involute-gear-2d 6 10 20 0))", "g");
    eval_points(&t, &CTX, &[[37.70576057, 4.08885592, 0.0], [37.705762, 4.0888561, 0.0]], &mut w);
    assert!((w[0] - w[1]).abs() < 1e-5, "{w:?}");
}

#[test]
fn grad_lists_the_derivatives_along_x_y_z() {
    // diff agrees with the evaluators' derivatives for every op (conformance);
    // grad is its three spatial derivatives, in order.
    let d = ink::load_str("(import \"std\") (def g (grad (+ (* 2 x) (* 3 y) (* 5 z) (* x y))))", "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Vec(g) = d.def("g").unwrap().clone() else { panic!("grad is not a vec") };
    let got: Vec<f64> = g
        .iter()
        .map(|v| match v {
            Value::Field(f, _) => eval_point(&Tape::new(&d.g, *f), &CTX, [1.0, 7.0, 0.0]),
            _ => panic!("{v:?}"),
        })
        .collect();
    assert_eq!(got, [9.0, 4.0, 5.0]);
}

#[test]
fn thickness_is_exact_on_exact_walls() {
    // A ray ends where the material does, to the precision of floats: a
    // hollow sphere's wall is 2 from either side, a 3 mm plate 3 (its
    // edges take their faces' thickness, or the diagonal, or none), and a
    // 0.8 fin on the plate is its thinnest wall.
    use std::sync::Arc;
    let thinnest = |src: &str, lo: f64, hi: f64| {
        let (_d, t) = field(&format!("(def part {src})"), "part");
        ink::feature::thickness(&Arc::new(t), &CTX, [lo; 3], [hi; 3], 64.0).unwrap().expect("a wall")
    };
    let (_, t) = thinnest("(subtract (sphere 10) (sphere 8))", -10.0, 10.0);
    assert!((t - 2.0).abs() < 1e-7, "{t}");
    let (_, t) = thinnest("(box 40 30 3)", -20.0, 20.0);
    assert!((t - 3.0).abs() < 1e-9, "{t}");
    let (at, t) = thinnest("(min (box 40 30 3) (move (box 0.8 20 10) 10 0 6))", -20.0, 20.0);
    assert!((t - 0.8).abs() < 1e-9 && (at[0] - 10.0).abs() <= 0.4 + 1e-9, "{t} at {at:?}");
}

#[test]
fn a_ray_finds_the_first_crossing_and_ends_when_it_lies_in_a_face() {
    use ink::query::{Side, ray_to};
    use std::sync::Arc;
    let ctx = Ctx::new(0.0, Vec::new());
    // Two slabs with a 0.001 gap: leaving the first, the ray finds the gap
    // (longer than eps), however far it would step through clear material.
    let (_d, t) = field("(def two (min (box 10 10 1) (move (box 10 10 1) 0 0 1.001)))", "two");
    let t = Arc::new(t);
    let local = ink::octree::Local::new(&t, &ctx, [-6.0; 3], 0.5);
    let sure = |h: Option<ink::query::RayHit>| h.filter(|h| h.rounding.is_none()).unwrap().t;
    let out = sure(ray_to(&local, [0.3, 0.2, -0.499], [0.0, 0.0, 1.0], (0.0, 10.0), Side::Outside, 1e-4));
    assert!((out - 0.999).abs() < 1e-12, "left at {out}");
    // And back in: the second slab's face, then its far face.
    let back = sure(ray_to(&local, [0.3, 0.2, -0.499], [0.0, 0.0, 1.0], (out, 10.0), Side::Inside, 1e-4));
    assert!((back - 1.0).abs() < 1e-12, "entered at {back}");
    let out = sure(ray_to(&local, [0.3, 0.2, -0.499], [0.0, 0.0, 1.0], (back, 10.0), Side::Outside, 1e-4));
    assert!((out - 2.0).abs() < 1e-12, "left at {out}");
    // A ray in a turned box's face, where only rounding decides the sign:
    // the search ends (its cost is the length over eps), with nothing or a
    // point on the face.
    let (_d, t) = field("(def b (rotate-z (box 40 30 3) 30))", "b");
    let t = Arc::new(t);
    let local = ink::octree::Local::new(&t, &ctx, [-30.0; 3], 1.0);
    // Up the face at x = 20 in the box's frame, from its corner.
    let (c, s) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
    let (o, dir) = ([20.0 * c + 15.0 * s, 20.0 * s - 15.0 * c, 0.0], [-s, c, 0.0]);
    for side in [Side::Inside, Side::Outside] {
        if let Some(ink::query::RayHit { t: at, .. }) = ray_to(&local, o, dir, (0.0, 30.0), side, 1e-3) {
            let p = [o[0] + at * dir[0], o[1] + at * dir[1], o[2] + at * dir[2]];
            assert!(eval_point(&t, &CTX, p).abs() < 1e-9, "{side:?} at {at}");
        }
    }
    // A ray along a plate's bottom face, from inside a bore: the field is
    // exactly 0 on the face, and its interval bound over any stretch of
    // the face is >= 0 (no rounding range around an exact 0), so the
    // whole face is ruled out at once, at a precision of 1e-9 of the ray,
    // not each stretch halved down to that eps (a ray on a 1/4 in plate,
    // minutes).
    let (_d, t) = field("(def p (extrude-z (subtract (rect-2d 114.3 50.8) (circle-2d 14.2875)) -3.175 3.175))", "p");
    let t = Arc::new(t);
    let (o, dir) = ([-6.5, -4.35, -3.175], [0.944681, -0.327991, 0.0]);
    let seg = |a: f64, b: f64| {
        std::array::from_fn(|k| {
            let (p, q) = (o[k] + a * dir[k], o[k] + b * dir[k]);
            ink::interval::I::new(p.min(q), p.max(q))
        })
    };
    assert!(ink::eval::eval_interval(&t, &CTX, seg(20.0, 30.0), None).lo >= 0.0);
    let local = ink::octree::Local::new(&t, &ctx, [-60.0; 3], 120.0 / 256.0);
    assert_eq!(ray_to(&local, o, dir, (0.0, 70.0), Side::Inside, 7e-8), None);
}

#[test]
fn a_ray_grazing_a_face_into_material_is_a_sure_hit() {
    // At 5e-9 to a cube's face, the ray's stretch where it crosses is flat
    // to rounding, but bounds show a point eps (`ink ray`'s, 1e-9 of the
    // ray) after it inside: a sure hit where it crosses, not one said only
    // to rounding.
    use ink::query::{Side, ray_to};
    let ctx = Ctx::new(0.0, Vec::new());
    let (_d, t) = field("(def b (box 20 20 20))", "b");
    let t = std::sync::Arc::new(t);
    let local = ink::octree::Local::new(&t, &ctx, [-30.0; 3], 60.0 / 256.0);
    let h = ray_to(&local, [-20.0, -10.00000007871675, 0.0], [1.0, 4.903941923740832e-9, 0.0], (0.0, 520.0), Side::Inside, 5.2e-7);
    assert!(h.is_some_and(|h| h.rounding.is_none() && (h.t - 16.05172961950302).abs() < 1e-6), "{h:?}");
}

#[test]
fn a_ray_ending_just_before_material_finds_none_past_its_end() {
    // A box with a tilted cut: along the ray the field stays above 0 to
    // t = 16, its end, and material starts just after. A hit's certainty
    // probe past t1 found it and called a hit inside the ray sure.
    use ink::query::{Side, ray_to};
    let ctx = Ctx::new(0.0, Vec::new());
    let (_d, t) = field("(def b (max (- (abs x) 10) (- (abs y) 20) (- (abs z) 20) (- (- (* 0.8660254037844387 z) (* 0.49999999999999994 y)) 1.5)))", "b");
    let t = std::sync::Arc::new(t);
    let local = ink::octree::Local::new(&t, &ctx, [-30.0; 3], 60.0 / 256.0);
    let h = ray_to(&local, [-20.0, -9.61378377305627, -3.818468241736015], [1.0, 0.0, -1e-7], (0.0, 16.0), Side::Inside, 1.6e-8);
    assert!(h.is_none_or(|h| h.rounding.is_some()), "a sure hit inside the ray, where material starts after it: {h:?}");
}

#[test]
fn a_ray_in_a_tilted_face_ends_at_rounding_and_says_so() {
    // A box turned 30 degrees about x, the ray in its top face (box frame
    // (-25, 2, 1.5) turned): the face's bounds hold 0 by rounding all
    // along it, and halving to 1e-9 of the ray took minutes. Now it ends
    // where it meets the face, a hit within rounding of 0; 1e-7 inside
    // the face, a sure hit at the same t.
    use ink::query::{RayHit, Side, ray_to};
    let ray = |o: [f64; 3]| {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let ctx = Ctx::new(0.0, Vec::new());
            let (_d, t) = field("(def b (rotate-x (box 40 30 3) 30))", "b");
            let t = std::sync::Arc::new(t);
            let local = ink::octree::Local::new(&t, &ctx, [-30.0; 3], 60.0 / 256.0);
            tx.send(ray_to(&local, o, [1.0, 0.0, 0.0], (0.0, 60.0), Side::Inside, 6e-8))
        });
        // A minute: the fixed ray answers in milliseconds, the old one took
        // minutes; 5 s timed out once in a loaded gate.
        rx.recv_timeout(std::time::Duration::from_secs(60)).expect("still refining after 60 s")
    };
    let face = ray([-25.0, 2.482050807568877, -0.2990381056766581]);
    assert!(face.is_some_and(|h| (h.t - 5.0).abs() < 1e-6 && h.rounding.is_some_and(|w| w < 1e-12)), "{face:?}");
    let inset = ray([-25.0, 2.482050757568877, -0.29903801907411776]);
    assert!(matches!(inset, Some(RayHit { t, rounding: None }) if (t - 5.0).abs() < 1e-9), "{inset:?}");
}

#[test]
fn sampling_an_unbounded_field_is_refused_up_front() {
    // A bare gyroid fills the whole domain: redistance would mesh all of it
    // (meshing all of it would exhaust memory). It must fail fast with a hint.
    let e = ink::load_str("(import \"std\")\n(def g (redistance (gyroid 12) (res 0.3)))", "test", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("unbounded") && e.hint.is_some(), "{e:?}");
}

#[test]
fn gradient_bounds_are_sound_and_tight_on_exact_shapes() {
    // Lipschitz bounds over boxes must cover every sampled gradient; exact
    // primitives (emax, lipschitz promises, rigid or uniformly scaled) stay
    // near 1, and a stretching warp rescales the promise instead of breaking it.
    let src = "(def rbox (rotate-z (rotate-x (box 4 3 2) 30) 20))
(def hexnut (extrude-z (max (- (abs x) 4) (- (abs (+ (* 0.5 x) (* 0.8660254 y))) 4) (- (abs (- (* 0.5 x) (* 0.8660254 y))) 4)) 0 3))
(def cyl (move (cylinder-z 2 5) 0.5 0 0))
(def small (scale (box 4 4 4) 0.5))
(def stretched (warp (box 2 2 2) (* 3 x) y z))
(def rod (thread-iso 10 1.5))";
    let mut r = Lcg(7);
    for (name, tight) in [("rbox", true), ("hexnut", false), ("cyl", true), ("small", true), ("stretched", false), ("rod", false)] {
        let (_d, t) = field(src, name);
        for _ in 0..300 {
            let c = [r.f(-6.0, 6.0), r.f(-6.0, 6.0), r.f(-3.0, 4.0)];
            let w = r.f(0.05, 1.5);
            let bx = c.map(|v| ink::interval::I::new(v - w, v + w));
            let bound = ink::interval::lipschitz(&ink::eval::eval_interval_grad(&t, &CTX, bx));
            let pts: Vec<[f64; 3]> = (0..20).map(|_| c.map(|v| v + r.f(-w, w))).collect();
            for g in ink::eval::eval_grads(&t, &CTX, &pts) {
                let n = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
                assert!(n <= bound * (1.0 + 1e-9), "{name}: |grad| {n} above bound {bound} on {bx:?}");
            }
            if tight {
                assert!(bound < 1.0 + 1e-6, "{name}: bound {bound} on {bx:?}");
            }
        }
    }
}

/// A derivative's range also holds it within its field's
/// gradient bounds (graph::Node::Deriv), which keep what interval
/// arithmetic on the formula loses. d/dx of the distance to a point is
/// x / sqrt(x^2 + y^2), a quotient of two ranges through 0 near the point,
/// and at most sqrt(sum (du_i/dx)^2) = 1 by Cauchy-Schwarz: proven over
/// the whole box at once. std's norm is sqrt (dot a a), a sum of products
/// a a, which the gradient bounds take as squares; across atan2's cut,
/// a slope along z stays 0 in the angle (its operands do not move along
/// z). The rule bounds the derivative's real value at every point, so it
/// takes no promise (a promise bounds the field's slope where it has
/// one, and these are false or say nothing at x = 0: 3 cos x promised 0.5,
/// and x written as a select whose other branch, 1000 x, is taken only at
/// x = 0, where diff follows it; no witness lies off that plane). About a
/// sampled field's source, its slope is unknown: the gradient bounds'
/// partials are the samples'. A warp rewrites a derivative into something
/// else (y, the slope of x y along x, read at (x/2, y), is not the slope
/// of x y / 2): the mark goes. A tie's mean is da/2 + db/2 (as
/// (da + db)/2, the mean of 1e308 and 1e308 would be +inf, outside the
/// bound it is proven within; derivative_bounds_hold_every_value_the_rules_give
/// checks the values), and the gradient bounds' hull makes room for its
/// halvings in floats: with slopes p = q = 5e-324 (params, so f64 values)
/// the mean is 0 at every tie, here everywhere, so 1e300 times it has
/// slope 0, which the bare hull of p and q would exclude, proving it
/// within [4e-24, 6e-24] (1e300 p is 4.94e-24). Those two models are
/// outside ink's range: their claims are unknown, not checked.
#[test]
fn derivative_ranges_keep_what_the_formula_loses() {
    use ink::claims::Status::*;
    let got = run_claims(
        "(assert repro (range (diff (len2 x y) x) -1.01 1.01 (box 2 2 2)))
         (assert steep (range (diff (len2 x y) x) -0.9 1.01 (box 2 2 2)))
         (assert segment (range (diff (segment-distance-at here2 (vec -5 0) (vec 5 0)) y) -1.01 1.01 (box 14 4 2)))
         (assert seam (range (diff (len2 (- (atan2 y x) pi) z) z) -1.01 1.01 (move (box 2 2 2) -5 0 0)))
         (assert false-promise (range (diff (lipschitz 0.5 (* 3 (sin x))) x) -1 1 (box 2 2 2)))
         (assert not-the-slope (range (diff (lipschitz 1 (select (neg (square x)) x (* 1000 x))) x) -2 2 (box 2 2 2)))
         (def baked (grid-bake (sphere 3) (res 0.25)))
         (assert baked-slope (range (diff baked x) -1.5 1.5 (move (box 1 1 1) 0 0 2)))",
    );
    // Each status, and whether it is an estimate (no accuracy): the one
    // that reads samples.
    let got: Vec<(&str, ink::claims::Status, bool)> = got.iter().map(|(n, c)| (n.as_str(), c.status, c.asked.is_none())).collect();
    assert_eq!(
        got,
        vec![
            ("repro", Holds, false),
            ("steep", Fails, false),
            ("segment", Holds, false),
            ("seam", Holds, false),
            ("false-promise", Fails, false),
            ("not-the-slope", Unknown, false),
            ("baked-slope", Holds, true),
        ]
    );
    for src in [
        "(assert tie-overflow (range (diff (min (* 1e308 x) (+ (* 1e308 x) y)) x) 0 1.1e308 (box 2 2 2)))",
        "(param p 5e-324) (param q 5e-324)
         (assert subnormal-tie (range (diff (* 1e300 (min (* p x) (* q x))) x) 4e-24 6e-24 (box 2 0.1 0.1)))",
    ] {
        let got = run_claims(src);
        assert!(!got.is_empty() && got.iter().all(|(_, c)| c.status == Unknown && matches!(c.evidence, ink::claims::Evidence::Outside(_))), "{src}");
    }
    let marked = |src: &str| field(src, "f").1.ssa.iter().any(|i| matches!(i.op, ink::tape::Op::Deriv(..)));
    assert!(marked("(def f (diff (* x y) x))"));
    assert!(!marked("(def f (warp (diff (* x y) x) (* 0.5 x) y z))"));
}

/// The bound of a derivative over a box holds its value in real numbers at
/// every point of it, as every range does: at ties (a
/// min of equal sides, abs at 0, both its average), where a quotient's
/// divisor is 0 (sqrt at 0: 0), on a select's switch, across atan2's cut,
/// at a clamp's ends. The points drawn are at moderate scale, where no
/// intermediate underflows or overflows and f64 evaluation gives the real
/// value to within 1e-9 of it; an infinite value would have to lie in the
/// bound itself. At extreme scales floats are not the real value (at
/// (2.7e-162, 0) the squares in d/dx len2 underflow and floats give 1.21):
/// those points carry their real values, worked by hand, and where floats
/// do agree (a tie's mean, halved first) that is checked too.
#[test]
fn derivative_bounds_hold_every_value_the_rules_give() {
    let fields = [
        "(len2 x y)",
        "(- (len3 (- x 0.3) y z) 1)",
        "(segment-distance-at here2 (vec -1 -1) (vec 1 2))",
        "(polyline-distance (vec (vec -2 0 0) (vec 0 0 0) (vec 1 1 1) (vec 1 -1 2)))",
        "(rotate-z (rotate-x (box 2 3 4) 30) 20)",
        "(move (cylinder-z 1 2) 0.5 0 0)",
        "(capsule-z 2 0.5)",
        "(torus-z 2 0.5)",
        "(smooth-union (sphere 1) (move (box 1 1 1) 1 0 0) 0.5)",
        "(max (abs x) (min y (- z)))",
        "(emax (- (abs x) 1) (- (abs y) 1))",
        "(select (- x 0.25) (square y) (* 2 x))",
        "(* (atan2 y x) (+ 1 z))",
        "(/ x (+ 0.5 (square y)))",
        "(sqrt (abs x))",
        "(+ (floor (* 3 x)) (* 0.5 y))",
        "(pow (abs x) 1.5)",
        "(len2 (- (atan2 y x) pi) z)",
    ];
    // (field, point, its slope along x there in real numbers, whether f64
    // evaluation gives it).
    let extreme: [(&str, [f64; 3], f64, bool); 8] = [
        // A tie of equal slopes is that slope; of 1e308 and 1.5e308, their
        // mean, though their sum overflows. Off the tie, a side's slope.
        ("(min (* 1e308 x) (+ (* 1e308 x) y))", [0.0, 0.0, 0.0], 1e308, true),
        ("(min (* 1e308 x) (+ (* 1.5e308 x) y))", [0.0, 0.0, 0.0], 0.5 * 1e308 + 0.5 * 1.5e308, true),
        ("(min (* 1e308 x) (+ (* 1.5e308 x) y))", [0.0, -1.0, 0.0], 1.5e308, true),
        ("(max (* -1e308 x) (- (* -1.5e308 x) y))", [0.0, 0.0, 0.0], 0.5 * -1e308 + 0.5 * -1.5e308, true),
        // Squares that underflow in floats: the slope is +-1.
        ("(len2 x y)", [2.7e-162, 0.0, 0.0], 1.0, false),
        ("(len2 x y)", [-2.7e-162, 0.0, 0.0], -1.0, false),
        ("(sqrt (square x))", [2.7e-162, 0.0, 0.0], 1.0, false),
        // And that overflow: sqrt((1e200 x)^2) = 1e200 |x|.
        ("(sqrt (square (* 1e200 x)))", [1.0, 0.0, 0.0], 1e200, false),
    ];
    let src = fields.iter().enumerate().map(|(i, f)| format!("(def f{i} {f}) (def g{i} (grad f{i}))\n")).collect::<String>()
        + &extreme.iter().enumerate().map(|(i, e)| format!("(def e{i} (diff {} x))\n", e.0)).collect::<String>();
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let within = |bound: ink::interval::I, v: f64| {
        assert!(!v.is_nan());
        if v.is_finite() {
            let slack = 1e-9 * (1.0 + v.abs());
            bound.lo - slack <= v && v <= bound.hi + slack
        } else {
            bound.contains(v)
        }
    };
    let mut r = Lcg(11);
    for i in 0..fields.len() {
        let Some(Value::Vec(g)) = d.def(&format!("g{i}")) else { panic!() };
        for (axis, gk) in g.iter().enumerate() {
            let Value::Field(gk, _) = gk else { panic!() };
            let t = Tape::new(&d.g, *gk);
            for trial in 0..200 {
                // Boxes about the origin's planes, so points on them (ties)
                // are drawn, and boxes anywhere.
                let (c, w) = if trial % 2 == 0 { ([0.0; 3].map(|_| r.f(-0.3, 0.3)), r.f(0.3, 1.5)) } else { ([r.f(-3.0, 3.0), r.f(-3.0, 3.0), r.f(-3.0, 3.0)], r.f(0.01, 1.0)) };
                let bx = c.map(|v| ink::interval::I::new(v - w, v + w));
                let bound = ink::eval::eval_interval(&t, &CTX, bx, None);
                let mut pts: Vec<[f64; 3]> = (0..24).map(|_| c.map(|v| v + r.f(-w, w))).collect();
                pts.push(c);
                for k in 0..pts.len() {
                    let mut p = pts[k];
                    for a in 0..3 {
                        if bx[a].lo <= 0.0 && 0.0 <= bx[a].hi && (k + a) % 3 == 0 {
                            p[a] = 0.0;
                        }
                    }
                    pts.push(p);
                }
                for p in pts {
                    let v = eval_point(&t, &CTX, p);
                    assert!(within(bound, v), "{} d/d{}: {v} at {p:?} outside {bound:?} over {bx:?}", fields[i], ["x", "y", "z"][axis]);
                }
            }
        }
    }
    for (i, &(f, p, real, floats)) in extreme.iter().enumerate() {
        let Some(Value::Field(e, _)) = d.def(&format!("e{i}")) else { panic!() };
        let t = Tape::new(&d.g, *e);
        if floats {
            assert_eq!(eval_point(&t, &CTX, p), real, "{f} at {p:?}");
        }
        for w in [0.0, 1e-300, 0.25, 2.0] {
            let bx = p.map(|v| ink::interval::I::new(v - w, v + w));
            let bound = ink::eval::eval_interval(&t, &CTX, bx, None);
            assert!(bound.contains(real), "{f} d/dx: {real} at {p:?} outside {bound:?} over {bx:?}");
        }
    }
}

/// A proven clearance d guarantees a gap of d, whatever the accuracy: the
/// accuracy only says how fine the search goes, so a gap within about it
/// of d is unknown, never proven (asking for no point within d less the
/// accuracy, a gap of 0.06 would prove (clearance a b 0.09) at 0.05).
/// Boxes 0.04 and 0.06 apart are disproven by a pair closer than 0.09, one
/// 0.095 apart is proven or unknown, one 0.2 apart proven. A clearance no
/// larger than the accuracy is asked the same way, no note, its search
/// refining to a sixteenth of its distance where that is finer than half
/// the accuracy (it sets its own scale): two prisms from polygons (no
/// slope bound at their surfaces, so the direct search decides) 0.01 apart
/// are 0.005 clear, to an accuracy below the gap and to the model's 0.05.
#[test]
fn a_proven_clearance_guarantees_its_distance() {
    use ink::claims::Status;
    let gaps = [0.04, 0.06, 0.095, 0.2];
    let src: String = gaps.iter().enumerate().map(|(k, g)| format!("(assert gap{k} (clearance (box 10 10 10) (move (box 10 10 10) (+ 10 {g}) 0 0) 0.09))\n")).collect();
    let c = run_claims(&format!(
        "{src}(fn quad (r) (polygon (vec (vec (neg r) (neg r)) (vec r (neg r)) (vec r r) (vec (neg r) r))))
         (def a (extrude-z (quad 1) -1 1))
         (def b (move (extrude-z (quad 0.5) -0.5 0.5) 1.51 0 0))
         (assert clear (clearance a b 0.005) (accuracy 0.001))
         (assert vague (clearance a b 0.005))"
    ));
    let got: Vec<Status> = c.iter().map(|(_, c)| c.status).collect();
    assert!(got[0] == Status::Fails && got[1] == Status::Fails, "{c:?}");
    assert!(matches!(got[2], Status::Holds | Status::Unknown) && got[3] == Status::Holds, "{c:?}");
    assert!(got[4] == Status::Holds && c[4].1.accuracy == Some(0.001), "{c:?}");
    assert!(got[5] == Status::Holds && c[5].1.accuracy == Some(0.05) && c[5].1.note.is_none(), "{c:?}");
}

/// Drafted faces close together clear by their planes: a drafted polygon's
/// line wall is its plane, signed and linear, so where both faces cross a
/// cell the pair's slopes cancel and a gap is proven over cells far coarser
/// than it. A 2 degree boss in a pocket 0.05 apart is 0.04 clear to 0.002;
/// one face alone (a block leaning with the boss's +x wall, 0.05 off it
/// away from its corners) at cells of 0.512. Read as an unsigned distance
/// signed by another field (a wall's |plane| either side), each face's
/// slopes would span both signs in a cell it crossed: the pocket unknown at
/// the cell budget, the one face proven only at 0.008.
#[test]
fn drafted_faces_close_together_clear_by_their_planes() {
    use ink::claims::{ClaimOptions, Evidence, Status, run_asserts};
    let src = "(import \"std\")
        (def boss (max (polygon (vec (vec -5 -5) (vec 5 -5) (vec 5 5) (vec -5 5)) 2 0 10) (slab-z 0 10)))
        (def pocket (max (polygon (vec (vec -5.05 -5.05) (vec 5.05 -5.05) (vec 5.05 5.05) (vec -5.05 5.05)) 2 0 10) (slab-z -0.1 11)))
        (def cup (subtract (max (polygon (vec (vec -10 -10) (vec 10 -10) (vec 10 10) (vec -10 10))) (slab-z -2 10)) pocket))
        (def block (max (neg (polygon (vec (vec -5.05 -5.05) (vec 5.05 -5.05) (vec 5.05 5.05) (vec -5.05 5.05)) 2 0 10)) (- 3 x) (- x 9) (- (abs y) 3) (slab-z 0.5 9.5)))
        (body b (shape boss))
        (body c (shape cup))
        (body d (shape block))
        (assert pocket (clearance boss cup 0.04) (accuracy 0.002))
        (assert face (clearance boss block 0.04) (accuracy 0.002))";
    let doc = ink::load_str(src, "gap", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let res = |a: &ink::claims::AssertResult| match &a.check.evidence {
        Evidence::ClearanceSearch { res, .. } => *res,
        e => panic!("{}: {e:?}", a.name),
    };
    assert!(r.asserts.iter().all(|a| a.check.status == Status::Holds), "{:?}", r.asserts.iter().map(|a| (&a.name, a.check.status)).collect::<Vec<_>>());
    assert!(res(&r.asserts[1]) >= 0.1, "the one face at cells of {}", res(&r.asserts[1]));
}

/// An overlap is disproven only where the shapes do not meet at all: an
/// intersection shown no deeper than the accuracy may still be one (a
/// press fit 0.04 deep at 0.05), so it is unknown. Two boxes face to face
/// share no point inside both (the bounds split at the face), so they do
/// not overlap. One 0.2 deep is proven, boxes 0.2 apart disproven.
#[test]
fn an_overlap_within_the_accuracy_is_not_disproven() {
    use ink::claims::Status;
    let c = run_claims(
        "(def a (box 10 10 10))
         (assert shallow (overlaps a (move (box 10 10 10) 9.96 0 0)))
         (assert face (overlaps a (move (box 10 10 10) 10 0 0)))
         (assert deep (overlaps a (move (box 10 10 10) 9.8 0 0)))
         (assert apart (overlaps a (move (box 10 10 10) 10.2 0 0)))",
    );
    let got: Vec<Status> = c.iter().map(|(_, c)| c.status).collect();
    assert_eq!(got, [Status::Unknown, Status::Fails, Status::Holds, Status::Fails], "{c:?}");
}

/// A clearance of 0 is disjointness: boxes resting face to face are proven
/// to overlap no deeper than the accuracy (a gap of 0 searched for
/// directly is never resolved where they touch), and boxes that overlap
/// deeper are disproven.
#[test]
fn a_clearance_of_nothing_is_disjointness() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let src = "(import \"std\")
(def a (box 20 20 4))
(assert touching (clearance a (move (box 20 20 4) 0 0 4) 0))
(assert overlapping (clearance a (move (box 20 20 4) 0 0 3.7) 0))";
    let doc = ink::load_str(src, "claims", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    assert_eq!(r.asserts[0].check.status, Status::Holds, "{r:?}");
    assert_eq!(r.asserts[0].check.accuracy, Some(0.05), "{r:?}");
    assert_eq!(r.asserts[1].check.status, Status::Fails, "{r:?}");
}

/// A lid's skirt, cut level by a half-space just under the box's rim, 0.3
/// clear of the box's walls, cleared by 0.28. Covering the box, the cells
/// under the lid's ceiling stay undecided at any size (the skirt's field
/// there reads the cut plane: 0.1 at points 0.36 from the skirt), so their
/// volume does not shrink and their count grows 5 to 7 times per halving;
/// covering the skirt decides it at cells of 0.016. Refining both ways in
/// step, the search spends this budget covering the box and ends unknown
/// (with a 0.01 margin, 62 of its 63 s; on a 64 mm box, 104 of 110 s).
/// Refining the more promising way, it is proven within it.
#[test]
fn a_clearance_refines_the_way_round_that_decides() {
    use ink::query::{Verdict, Way, clearance};
    use std::sync::Arc;
    let src = "(def outline (rounded-rect-2d 16 12 3))
(def box (subtract (extrude-z outline 0 10) (extrude-z (- outline -2) 2 11)))
(def lid (subtract (extrude-z (- outline 2.3) 6 12) (extrude-z (- outline 0.3) 5 10)))
(def skirt (max lid (- z 9.8)))";
    let (d, _) = field(src, "box");
    let tape = |n: &str| {
        let Value::Field(f, _) = d.def(n).unwrap().clone() else { panic!() };
        Arc::new(Tape::new(&d.g, f))
    };
    let (bx, sk) = (tape("box"), tape("skirt"));
    // Each shape's box within the other's grown by the clearance; the claim
    // either way round (box first, then skirt first), so the deciding way
    // is first, then second: an order the search followed would fail one.
    let covers_skirt = |swapped| Way { grown: bx.clone(), covered: sk.clone(), lo: [-8.3, -6.3, 6.0], hi: [8.3, 6.3, 9.8], swapped };
    let covers_box = |swapped| Way { grown: sk.clone(), covered: bx.clone(), lo: [-8.0, -6.0, 5.6], hi: [8.0, 6.0, 10.0], swapped };
    for ways in [[covers_skirt(false), covers_box(true)], [covers_box(false), covers_skirt(true)]] {
        let (v, _, h) = clearance(&ways, &CTX, 0.28, 1e-3, ink::query::BOXES);
        assert!(matches!(v, Verdict::Proven) && h < 0.02, "{v:?} at {h}");
        let (v, pair, _) = clearance(&ways, &CTX, 0.31, 1e-3, ink::query::BOXES);
        assert!(matches!(v, Verdict::Disproven { .. }) && pair.is_some(), "{v:?}");
    }
}

#[test]
fn micron_gaps_prove_with_coarse_cells() {
    use ink::query::Verdict;
    use std::sync::Arc;
    // A bore of radius 10 and a shaft 5 um smaller: disjoint, proven with
    // 0.1 mm cells (first-order bounds of the facing surfaces; plain bounds
    // need cells under 5 um); with the shaft 1 um too big, a witness.
    let tape = |shaft: f64| {
        let src = format!("(def pair (max (subtract (cylinder-z 15 20) (cylinder-z 10 30)) (cylinder-z {shaft} 10)))");
        let (_d, t) = field(&src, "pair");
        Arc::new(t)
    };
    let (lo, hi) = ([-11.0, -11.0, -6.0], [11.0, 11.0, 6.0]);
    assert!(matches!(empty_at(&tape(9.995), &CTX, lo, hi, Some(0.1)).0, Verdict::Proven));
    assert!(matches!(empty_at(&tape(10.001), &CTX, lo, hi, Some(0.1)).0, Verdict::Disproven { .. }));
}

#[test]
fn varying_sweeps() {
    // A circle shrinking from r 4 to 2 along a straight 40 mm path is a
    // frustum; a twisting 10 x 2 bar keeps its area.
    let path = "(def straight (vec (vec 0 0 0) (vec 0 0 10) (vec 0 0 20) (vec 0 0 40)))";
    let frustum = PI * 40.0 * (16.0 + 8.0 + 4.0) / 3.0;
    let v = volume(&format!("{path}\n(body b (shape (sweep-varying (lambda (s) (circle-2d (mix 4 2 s))) straight (vec 1 0 0))))"), 0.1);
    assert!((v - frustum).abs() < 0.002 * frustum, "frustum {v} vs {frustum}");
    // A twisting section is bounded by a containing circle: interval bounds
    // cannot show a twist bounded far out, and the sweep's pieces are
    // bounded by the sections' proven reach (without one it is an error).
    let twisted = |s: &str| format!("{path}\n(body b (shape (sweep-varying (lambda (s) {s}) straight (vec 1 0 0))))");
    let v = volume(&twisted("(max (circle-2d 6) (rotate-z (rect-2d 10 2) (* 90 s)))"), 0.1);
    assert!((v - 800.0).abs() < 0.002 * 800.0, "twisted bar {v}");
    let e = ink::load_str(&format!("(import \"std\")\n{}", twisted("(rotate-z (rect-2d 10 2) (* 90 s))")), "t", Path::new(ROOT)).err().expect("an unbounded twist");
    assert!(format!("{e:?}").contains("not proven bounded"), "{e:?}");
}

#[test]
fn drafted_rims_are_not_thin() {
    // A 2 mm wall with 2 degree draft: rays from the rim near its edge exit
    // through the drafted face; those measure the corner, not the wall.
    let (_d, t) = field(
        "(def ta (tan (deg->rad 2)))
         (def outer (* (cos (deg->rad 2)) (max (- (abs x) 40 (* ta (- z 30))) (- (abs y) 30 (* ta (- z 30))))))
         (def cup (max (max outer (neg z) (- z 30)) (neg (max (+ outer 2) (- 2 z)))))",
        "cup",
    );
    let t = std::sync::Arc::new(t);
    let (_, min) = ink::feature::thickness(&t, &CTX, [-42.0, -32.0, -1.0], [42.0, 32.0, 31.0], 128.0).unwrap().expect("a wall");
    assert!(min > 1.9, "thinnest {min}");
}

#[test]
fn promises_never_outlive_the_field_they_were_made_for() {
    use ink::query::Verdict;
    use ink::tape::Op;
    use std::sync::Arc;
    // Promises in a field (std's primitives make their own too).
    let promises = |src: &str, name: &str| {
        let (_d, t) = field(src, name);
        t.code.iter().filter(|i| matches!(i.op, Op::Lipschitz)).count()
    };
    let kept = |with: &str, without: &str| promises(with, "f") > promises(without, "f");
    // A rigid warp keeps a promise (its frame moves with the field) ...
    assert!(kept("(def f (move (lipschitz 1 (- x 1)) 3 0 0))", "(def f (move (- x 1) 3 0 0))"));
    // ... a warp of t to an expression of x, y, z: the spatial
    // slope gains d/dt dt/dx, which no promise covers. With the promise,
    // this clearance of 1.9 would be proven against a true gap of 1.6.
    let src = "(def slab (max (warp (lipschitz 1 (- (abs (- x (* 5 t))) 1)) x y z (* -0.2 x)) (box 100 20 20)))
(def ball (move (sphere 0.1) 2.2 0 0))";
    assert!(!kept(&src.replace("(def slab", "(def f"), &src.replace("(def slab", "(def f").replace("(lipschitz 1 ", "(+ 0 ")));
    // A warp of t to t alone keeps it.
    assert!(kept("(def f (warp (lipschitz 1 (- (abs (- x (* 5 t))) 1)) x y z (* 2 t)))", "(def f (warp (- (abs (- x (* 5 t))) 1) x y z (* 2 t)))"));
    let (d, _) = field(src, "slab");
    let tape = |n: &str| {
        let Value::Field(f, _) = d.def(n).unwrap().clone() else { panic!() };
        Arc::new(Tape::new(&d.g, f))
    };
    let (lo, hi) = ([-6.0, -11.0, -11.0], [6.0, 11.0, 11.0]);
    assert!(!matches!(clearance_of(&tape("slab"), &tape("ball"), &CTX, 1.9, lo, hi, 0.02).0, Verdict::Proven));
}

#[test]
fn bounds_never_claim_a_box_the_domain_cut() {
    // A 2 m bar reaches past the window (+-500 mm): its load-time bounds
    // are all of it, not the box the window cut (which would disprove
    // (value bx > 900) "exactly"). A shape that is not bounded has none.
    let d = ink::load_str("(import \"std\")\n(def b (bounds (box 2000 10 10)))\n(def bx (nth (nth b 1) 0))\n", "t", Path::new(ROOT)).unwrap();
    let Some(ink::lang::Value::Field(bx, _)) = d.def("bx") else { panic!("no bx") };
    let bx = d.g.as_const(*bx).unwrap();
    // All of it, tight to within the last pass's cell (1/256 of the box).
    assert!((1000.0..1000.0 + 2.0 * 2000.0 / 256.0).contains(&bx), "{bx}");
    let e = ink::load_str("(import \"std\")\n(def b (bounds (max x (sphere 1e5))))\n(def c (bounds x))\n", "t", Path::new(ROOT)).err().expect("bounds of a half-space");
    assert!(format!("{e:?}").contains("not proven bounded") && format!("{e:?}").contains("line: 3"), "{e:?}");
}

#[test]
fn promised_threads_prove_clearance_across_their_cut() {
    use ink::query::Verdict;
    use std::sync::Arc;
    // A tapped hole 0.04 bigger than its screw, both from one promised
    // thread field T: T >= 0.04 in the nut, T < 0 in the screw. Cells
    // crossing the thread's atan2 cut have no slope bound from its formula;
    // the promise bounds T whole, so the 0.02 half-gap proves with 0.02
    // cells (looking through the promise would leave it unknown).
    // Turned so no cell centre lies on the cut itself, where even the
    // centre's value has wide bounds.
    let src = "(def th (rotate-z (thread-iso 3 0.5) 7))
(def pair (max th (neg (- th 0.04)) (slab-z -0.25 0.25)))";
    let (_d, t) = field(src, "pair");
    let (lo, hi) = ([-1.6, -1.6, -0.3], [1.6, 1.6, 0.3]);
    let (v, _) = empty_at(&Arc::new(t), &CTX, lo, hi, Some(0.02));
    assert!(matches!(v, Verdict::Proven), "{v:?}");
}

#[test]
fn coarse_threads_make_no_promise() {
    // M1 with pitch 1: the root radius is negative, the flank reaches the
    // axis and the field jumps there, so no slope bound holds and std
    // promises none. A promise would certify this cell empty: its centre
    // sees the root cylinder (+0.04) while a point by the axis is -0.29.
    use ink::interval::I;
    let promised = |t: &Tape| t.ssa.iter().any(|i| matches!(i.op, ink::tape::Op::Lipschitz));
    let (_d, t) = field("(def f (rotate-z (thread-iso 1 1) 40))", "f");
    let b = [I::new(-1e-8, 9.9e-7), I::new(-1e-8, 9.9e-7), I::new(-5e-7, 5e-7)];
    let q = eval_point(&t, &CTX, [-5e-9, -5e-9, 0.0]);
    assert!(q < -0.29, "{q}");
    assert!(eval_point(&t, &CTX, [4.9e-7, 4.9e-7, 0.0]) > 0.04);
    let lo = ink::centred::centred_lower(&t, &CTX, b).unwrap();
    assert!(lo <= q, "lower bound {lo} above f = {q}");
    assert!(!promised(&t));
    // A fine pitch for its size keeps the promise (the catalog's threads).
    assert!(promised(&field("(def f (thread-iso 3 0.5))", "f").1));
}

/// A promise made in a function holds for its inputs held fixed, so a
/// size that reads x, y or z makes none: every promise site in std, with
/// such a size and with a fixed one. Each would otherwise make a false
/// promise (|grad| 7.6 and 42 where 1 is promised), trusted by proofs.
#[test]
fn std_promises_hold_for_fixed_sizes_only() {
    let promised = |src: &str| field(src, "f").1.ssa.iter().any(|i| matches!(i.op, ink::tape::Op::Lipschitz));
    let w = "(+ 2 (sin (* 10 x)))";
    for site in [
        "(sphere R)",
        "(box (* 6 R) 2 2)",
        "(cylinder-z R 2)",
        "(torus-z 3 R)",
        "(capsule-z 2 R)",
        "(circle-2d R)",
        "(rect-2d (* 6 R) 2)",
        "(thread-iso (+ 8 R) 1)",
        "(sweep-varying (lambda (s) (circle-2d R)) (vec (vec 0 0 0) (vec 0 0 1)) (vec 1 0 0))",
        "(gyroid R)",
    ] {
        assert!(promised(&format!("(def f {})", site.replace('R', "2"))), "{site} fixed");
        assert!(!promised(&format!("(def f {})", site.replace('R', w))), "{site} spatial");
    }
    // A section's size at s (z within the section, warped to s along the path).
    assert!(!promised("(def f (sweep-varying (lambda (s) (circle-2d (mix 1 11 s))) (vec (vec 0 0 0) (vec 0 0 1)) (vec 1 0 0)))"));
    // t and params are not space: a promise holds for each of their values.
    assert!(promised("(param r 2) (def f (sphere (+ r t)))"));
    // A value a lambda captured is its input too.
    assert!(!promised(&format!("(fn g (r) ((lambda () (lipschitz 1 (- (len3 x y z) r))))) (def f (g {w}))")));
    // half-space promises nothing: its slope is exact without one.
    assert!(!promised("(def f (half-space 1 2 3 4))"));
}

/// The gyroid is divided by its greatest gradient, sqrt(3) k (std.ink
/// derives it), so its promise of 1 is true and a shell of one is at
/// least as thick as asked (divided by 1.5 k: slopes up to 1.15, so walls
/// up to 13% thinner). Gradients sampled over a cell stay within 1,
/// and reach it where std.ink says (the last point).
#[test]
fn lattices_are_divided_by_their_greatest_gradient() {
    let p = 7.0;
    let mut pts: Vec<[f64; 3]> = (0..20_000)
        .map(|i| {
            let f = |k: usize| ((i as f64 + 0.5) * [0.618_033_988_7, 0.414_213_562_4, 0.732_050_807_6][k]).fract() * p;
            [f(0), f(1), f(2)]
        })
        .collect();
    pts.extend([[0.75 * p, 0.25 * p, 0.25 * p]]);
    {
        let name = "gyroid";
        let (_d, t) = field(&format!("(def f ({name} {p}))"), "f");
        let g = ink::eval::eval_grads(&t, &Ctx::new(0.0, vec![]), &pts);
        let most = g.iter().map(|g| (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt()).fold(0.0, f64::max);
        assert!((1.0 - 1e-12..=1.0 + 1e-12).contains(&most), "{name}: greatest sampled gradient {most}");
    }
}

/// std's graded-lattice helpers give a sheet the relative density asked,
/// measured over two periods each way (whole cells) at res 0.1: within
/// 0.0003 here (the gyroid at 0.05 and 0.9 too), against the fit's 0.0011
/// at worst (a gyroid at 0.95). Taking the level for the density, a gyroid
/// asked 0.3 would be 0.194.
#[test]
fn graded_lattices_take_the_density_asked() {
    for (name, rho) in [("gyroid", 0.2), ("gyroid", 0.6)] {
        let src = format!("(body b (shape (lattice-fill (box 20 20 20) ({name} 10) ({name}-thickness {rho} 10))))");
        let got = volume(&src, 0.1) / 8000.0;
        assert!((got - rho).abs() < 0.002, "{name} at {rho}: {got}");
    }
}

/// The false proofs those promises would give, and true controls. A
/// sphere whose radius reads x (in [1, 19]) and a section tapering 1 to 11
/// over 1 mm (slope 10) each come 0.2 and 0.11 from a probe, and would be
/// proven clear by 0.3 and 0.2 with L = 1. A half-space with normal
/// 2.7e-162 promising |n| from f64 squares, which underflow to 2.2e-162,
/// would be proven 1.1 clear of a probe 0.99 away. (The probes are
/// spheres: their promises are the "std" the half-space control lists.)
#[test]
fn false_promises_no_longer_prove_clearances() {
    use ink::claims::Status::*;
    let got = assert_statuses(
        "(import \"std\")
         (def pts (vec (vec 0 0 0) (vec 0 0 1)))
         (assert spatial-radius (clearance (sphere (+ 10 (* 9 (sin x)))) (* 100 (move (sphere 0.01) -1.5 0 0)) 0.3))
         (assert fixed-radius (clearance (sphere 10) (move (sphere 1) 12 0 0) 0.9))
         (assert taper (clearance (sweep-varying (lambda (s) (circle-2d (mix 1 11 s))) pts (vec 1 0 0)) (* 100 (move (sphere 0.001) 5 0 0.3)) 0.2))
         (assert straight (clearance (sweep-varying (lambda (s) (circle-2d 4)) pts (vec 1 0 0)) (move (sphere 0.5) 5.5 0 0.5) 0.9))
         (assert half-space (clearance (half-space 1 0 0 0) (move (sphere 0.01) 1 0 0) 0.9))",
        Some(0.01),
    );
    let got: Vec<(&str, ink::claims::Status)> = got.iter().map(|(n, c)| (n.as_str(), c.status)).collect();
    assert_eq!(got, vec![("spatial-radius", Fails), ("fixed-radius", Holds), ("taper", Fails), ("straight", Holds), ("half-space", Holds),]);
    // 2.7e-162 is outside ink's range: unknown, not checked.
    let got = assert_statuses("(import \"std\") (assert underflow (clearance (half-space 2.7e-162 0 0 0) (move (sphere 0.01) 1 0 0) 1.1))", Some(0.01));
    assert!(got.len() == 1 && matches!(got[0].1.evidence, ink::claims::Evidence::Outside(_)), "{got:?}");
}

/// std's primitives over their whole argument domains. sweep-varying
/// bounds its pieces by the sections' reach, proven over every s (a guess
/// from 17 samples x 1.25 would cut a section bulging to radius 11 between
/// samples, radius 1 at each, to 1.9 without a word). A negative scale
/// keeps the solid's sign (it scales by abs(s)); a capsule with h < 0 spans
/// h to 0, not collapsing to its far end; a straight pipe's flat ends are not distances to
/// spheres cut by max (0.55 where the rim is 0.71 away).
#[test]
fn std_primitives_hold_over_their_domains() {
    let at = |src: &str, ps: &[[f64; 3]]| {
        let (_d, t) = field(src, "f");
        ps.iter().map(|&p| eval_point(&t, &CTX, p)).collect::<Vec<_>>()
    };
    let bump = "(fn sec (s) (circle-2d (+ 1 (* 10 (exp (neg (square (* 160 (- s 0.03125)))))))))";
    let swept = at(&format!("{bump} (def f (sweep-varying sec (vec (vec 0 0 0) (vec 0 0 100)) (vec 1 0 0)))"), &[[5.0, 0.0, 3.125], [10.9, 0.0, 3.125], [11.1, 0.0, 3.125]]);
    let reference = at(&format!("{bump} (def f (max (sec (/ z 100)) (slab-z 0 100)))"), &[[5.0, 0.0, 3.125], [10.9, 0.0, 3.125], [11.1, 0.0, 3.125]]);
    assert_eq!(swept, reference);
    assert!(reference[1] < 0.0 && reference[2] > 0.0, "{reference:?}");
    assert_eq!(at("(def f (scale (sphere 1) -2))", &[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, -1.5, 0.0]]), [-2.0, 1.0, -0.5]);
    assert_eq!(at("(def f (capsule-z -2 1))", &[[0.0, 0.0, 0.0], [0.0, 0.0, -2.0], [0.0, 0.0, 1.0], [0.0, 0.0, -3.5], [2.0, 0.0, -1.0]]), [-1.0, -1.0, 0.0, 0.5, 1.0]);
    // A straight pipe is an annulus extruded between flat ends, exactly.
    let (_d, pipe) = field("(def f (pipe (vec (vec 1 2 3) (vec 7 -1 5)) 2 1))", "f");
    let (_d, exact) = field("(def f (place (extrude-z (shell (circle-2d 1.5) 1) 0 7) (nth (path-frames (vec (vec 1 2 3) (vec 7 -1 5)) (vec 0 0 1)) 0)))", "f");
    let mut r = Lcg(7);
    for _ in 0..500 {
        let p = [r.f(-3.0, 11.0), r.f(-5.0, 6.0), r.f(-1.0, 9.0)];
        let (a, b) = (eval_point(&pipe, &CTX, p), eval_point(&exact, &CTX, p));
        assert!((a - b).abs() < 1e-9, "at {p:?}: {a} vs {b}");
    }
}

#[test]
fn coincident_faces_both_form_surface() {
    // A tagged datum plane on the rectangle's own bottom face: at a point
    // where rounding makes the rectangle's field exactly 0 and the plane's
    // -2e-16, the deciding branch is the rectangle. Over a cell of the face
    // neither branch is ruled out, so both form it.
    let src = "(def rail (max (tag bottom (neg z)) (tag section (warp (extrude-z (move (rect-2d 12 8) 0 4 0) 0 300) y z x))))";
    let (d, t) = field(src, "rail");
    let names = |ids: &[u32]| -> Vec<String> {
        let mut v: Vec<String> = ids.iter().map(|&i| d.g.tag_names[i as usize].clone()).collect();
        v.sort();
        v
    };
    let p = [150.0, 0.0, 2.220446049250313e-16];
    assert_eq!(eval_point(&t, &CTX, p), 0.0);
    assert_eq!(names(&ink::eval::provenance(&t, &CTX, p)), vec!["rail", "section"]);
    let t = std::sync::Arc::new(t);
    let formed = ink::feature::surface_names(&t, &CTX, [-1.0, -7.0, -1.0], [301.0, 7.0, 9.0], 128.0).0;
    assert_eq!(names(&formed), vec!["bottom", "rail", "section"]);
}

/// `check` warns only about a name proven to form no surface: a face
/// thinner than a cell is found (one sample per cell, its centre projected
/// onto the top or bottom face, never named the 0.5 mm sides of a 350 mm
/// sheet), and a construction proven covered is not.
#[test]
fn thin_faces_form_surface() {
    let src = "(import \"std\")\n\
        (body sheet (shape (max (tag outline (rect-2d 350 350)) (tag faces (slab-z 0 0.5)))))\n\
        (body capped (shape (min (tag plate (box 40 40 4)) (tag buried (sphere 1)))))\n";
    let d = ink::load_str(src, "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let hidden = |body: &str, lo: [f64; 3], hi: [f64; 3]| -> Vec<String> {
        let shape = d.body(body).unwrap().shape;
        let t = std::sync::Arc::new(Tape::new(&d.g, shape));
        let formed = ink::feature::surface_names(&t, &CTX, lo, hi, 128.0).0;
        d.named_parts(body, shape).into_iter().filter(|t| !formed.contains(t)).map(|i| d.g.tag_names[i as usize].clone()).collect()
    };
    assert!(hidden("sheet", [-175.0, -175.0, 0.0], [175.0, 175.0, 0.5]).is_empty());
    assert_eq!(hidden("capped", [-20.0, -20.0, -2.0], [20.0, 20.0, 2.0]), vec!["buried"]);
}

/// A name that reaches the surface only where it ties with what covers it
/// (an edge) or weighted out of a blend is not seen forming it, though
/// intervals keep it on the field's path there: `check` reports it as
/// undecided. A face it shares with another is seen (coincident_faces).
/// Scaling either side, or the whole field, by any positive factor changes
/// none of it: each side's tie is read by its own slope.
#[test]
fn edge_only_and_weighted_out_names_are_unseen() {
    let scales = ["1", "1e-9", "1e9"];
    for kb in scales {
        for kc in scales {
            let src = format!(
                "(import \"std\")\n\
                (body edged (shape (max (tag block (* {kb} (box 10 10 10))) (tag cut (* {kc} (- (+ x z) 10))))))\n\
                (body weighed (shape (let ((h (clamp (- x 50) 0 1))) (+ (* (- 1 h) (tag cube (* {kb} (box 10 10 10)))) (* h (tag ghost (* {kc} (sphere 3))))))))\n"
            );
            let d = ink::load_str(&src, "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
            for (body, unseen) in [("edged", "cut"), ("weighed", "ghost")] {
                let shape = d.body(body).unwrap().shape;
                let t = std::sync::Arc::new(Tape::new(&d.g, shape));
                let (formed, maybe) = ink::feature::surface_names(&t, &CTX, [-5.0; 3], [5.0; 3], 128.0);
                let maybe: Vec<u32> = maybe.iter().map(|u| u.0).collect();
                let names = |v: Vec<u32>| -> Vec<String> { d.named_parts(body, shape).into_iter().filter(|t| v.contains(t)).map(|i| d.g.tag_names[i as usize].clone()).collect() };
                assert_eq!(names(maybe), vec![unseen], "{body}, its sides scaled by {kb} and {kc}");
                assert_eq!(names(formed).len(), 1, "{body}, its sides scaled by {kb} and {kc}: the box is seen");
            }
        }
    }
}

/// An unseen name is sought finer in up to 64 of the cells it is on a path
/// in, spread over them, and says how many: an edge 2 mm long is on a path
/// in over 64 cells (212 here), and sought finer in 64 of them (a stride of
/// 4 picked 53).
#[test]
fn an_unseen_name_says_in_how_many_cells_it_was_sought() {
    let src = "(import \"std\")\n(body edged (shape (max (tag block (box 10 2 10)) (tag cut (- (+ x z) 10)))))\n";
    let d = ink::load_str(src, "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let shape = d.body("edged").unwrap().shape;
    let t = std::sync::Arc::new(Tape::new(&d.g, shape));
    let (_, unseen) = ink::feature::surface_names(&t, &CTX, [-5.0, -1.0, -5.0], [5.0, 1.0, 5.0], 128.0);
    let &(_, on, finer) = unseen.first().expect("cut is unseen");
    assert!(on > 64 && finer == 64, "on a path in {on} cells, sought finer in {finer}");
}

/// A name's cells searched finer are all those picked on its paths, for it
/// or for another name: a tag wrapping two unseen edge planes is searched
/// finer in its parts' cells too, more than the 64 picked for it alone.
#[test]
fn an_unseen_names_count_includes_cells_picked_for_others() {
    let src = "(import \"std\")\n(body solid (shape (max (tag block (box 10 10 10)) (tag group (max (tag a (- (+ x z) 10)) (tag b (- (+ x y) 10)))))))\n";
    let d = ink::load_str(src, "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let shape = d.body("solid").unwrap().shape;
    let t = std::sync::Arc::new(Tape::new(&d.g, shape));
    let (_, unseen) = ink::feature::surface_names(&t, &CTX, [-5.0, -5.0, -5.0], [5.0, 5.0, 5.0], 128.0);
    let group = d.g.tag_names.iter().position(|n| n == "group").unwrap() as u32;
    let &(_, _, finer) = unseen.iter().find(|u| u.0 == group).expect("group is unseen");
    let parts: usize = unseen.iter().filter(|u| u.0 != group).map(|u| u.2).max().unwrap_or(0);
    assert!(finer > 64 && finer >= parts, "group sought finer in {finer}, its parts in up to {parts}");
}

/// `check`'s hidden-construction warning does not look inside other
/// bodies: another body used in this one is one construction, under its
/// own name. A cutter that cuts nothing is still reported, the cutter's
/// own features are not.
#[test]
fn hidden_constructions_stop_at_other_bodies() {
    let src = "(import \"std\")\n\
        (body block (shape (box 30 30 4)))\n\
        (body cutter (shape (min (tag case (box 10 10 10)) (tag lug (move (sphere 1) 20 0 0)))))\n\
        (body gauge (shape (move (cylinder-z 1 10) 40 0 0)))\n\
        (body fitted (shape (subtract (shape-of block) (min (shape-of cutter) (shape-of gauge) (tag hole (move (cylinder-z 1 10) 0 40 0))))))\n";
    let d = ink::load_str(src, "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let hidden = |body: &str| -> Vec<String> {
        let shape = d.body(body).unwrap().shape;
        let t = std::sync::Arc::new(Tape::new(&d.g, shape));
        let formed = ink::feature::surface_names(&t, &CTX, [-16.0, -16.0, -6.0], [16.0, 16.0, 6.0], 128.0).0;
        let hidden = d.named_parts(body, shape).into_iter().filter(|t| !formed.contains(t));
        let mut v: Vec<String> = hidden.map(|i| d.g.tag_names[i as usize].clone()).collect();
        v.sort();
        v
    };
    // The cutter forms the pocket; its lug lies outside the block, which
    // is the cutter's business, not the fitted part's. The gauge and the
    // inline hole cut nothing.
    assert_eq!(hidden("fitted"), vec!["gauge", "hole"]);
    // On the cutter itself both of its features form surface.
    assert!(hidden("cutter").is_empty());
}

/// A claim where surfaces truly touch never resolves; refinement stops at
/// a floor set by the model's accuracy (half of it; the accuracy is
/// ACCURACY of the window), landing on it, and the claim is proven to that
/// accuracy, no finer: finer cells would cost more and say no more. Asked
/// exactly, the search ends at the first box undecided at the floor (fail
/// fast: its cells reach the floor, no finer). A given resolution still goes
/// finer, and a tiny probe's cells follow its region, below the floor.
/// (Faces that touch where both pieces read one coordinate are split there
/// and proven exactly; here one reads x as (x + y) - y, sharing nothing.)
#[test]
fn touching_surfaces_are_proven_to_the_models_accuracy() {
    use ink::query::{Verdict, prove_empty};
    use ink::ranged::{ProofOptions, prove};
    use std::sync::Arc;
    // A 1 mm cube seated on a face of a 2 mm one.
    // Searched in a 100 mm window, to that size's accuracy (5 um).
    let (_, t) = field("(def f (max (box 2 2 2) (warp (move (box 1 1 1) 1.5 0 0) (- (+ x y) y) y z)))", "f");
    let accuracy = 100.0 * ink::lang::ACCURACY;
    let finest = accuracy / 2.0;
    let (lo, hi) = ([-50.0; 3], [50.0; 3]);
    let (_, shared) = field("(def f (max (box 2 2 2) (move (box 1 1 1) 1.5 0 0)))", "f");
    assert!(matches!(prove_empty(&Arc::new(shared), &CTX, lo, hi, None, finest, 0.0).verdict, Verdict::Proven));
    let tape = Arc::new(t);
    // The cells it reached and the undecided cell's least side across the
    // face (y and z, where no clip narrows it), in floors: a cell is not
    // halved along an axis its question reads not there (one along an
    // edge of the face), so one of the two may be wider.
    let floors = |p: &ink::ranged::Proof| match p.verdict {
        Verdict::Unknown { cell: Some([a, b]), .. } => [p.res, (b[1] - a[1]).min(b[2] - a[2])].map(|w| (w / finest * 1e6).round() / 1e6),
        _ => panic!("touching faces stay undecided asked exactly: {p:?}"),
    };
    let s = prove_empty(&tape, &CTX, lo, hi, None, finest, 0.0);
    assert!(floors(&s) == [1.0; 2], "lands on the floor: {s:?}");
    // Asked to the accuracy, they are proven at once, at cells no finer
    // than half of it, not exactly.
    let opts = ProofOptions { res: None, finest, depth: 0, tol: accuracy };
    let p = prove(&tape, &CTX, lo, hi, &opts);
    assert!(p.verdict == Verdict::Proven && !p.exact && p.res >= finest && p.tol == accuracy, "{p:?}");
    // The default window's floor (0.025 mm) stops the exact search coarser.
    let h0 = prove_empty(&tape, &CTX, lo, hi, None, 1000.0 * ink::lang::ACCURACY / 2.0, 0.0);
    assert_eq!(floors(&h0), [10.0; 2]);
    // A given resolution (--res) is used as given, below the floor.
    let (v, hr) = empty_at(&tape, &CTX, lo, hi, Some(0.75 * finest));
    assert!(matches!(v, Verdict::Unknown { .. }) && hr == 0.75 * finest, "{v:?} at {hr}");
    // The cells follow the claim's region, below the floor if the region
    // is: a 0.03 mm probe seated on the part is searched to cells near 1e-3
    // mm, no finer.
    let (_, probe) = field("(def f (max (box 2 2 2) (warp (move (box 0.03 0.03 0.03) 1.015 0 0) (- (+ x y) y) y z)))", "f");
    let probe = Arc::new(probe);
    let s = prove_empty(&probe, &CTX, lo, hi, None, finest, 0.0);
    let start = prove_empty(&probe, &CTX, lo, hi, None, f64::INFINITY, 0.0).res;
    assert!(matches!(s.verdict, Verdict::Unknown { .. }) && s.res == start && s.res < finest, "{s:?}, start {start}");
}

/// A narrow gap around material that fills its window (a part modelled in
/// a snug domain) is halved down to the floor, however many halvings that
/// takes, and decided. Halving at most six times from 1/32 of the material
/// would stop at 32 floors, and the accuracy there, 64 floors, is too loose
/// to give (ranged::within_limit): such claims would be unknown on any
/// machine (a bench's working-height envelope). Upright, its walls read no
/// z inside the envelope's span, and the gap is decided as in the plane, at
/// coarser cells; tilted, they read every axis.
/// The answer is the same at any thread count.
#[test]
fn a_narrow_gap_around_a_whole_part_is_refined_to_the_floor() {
    use ink::ranged::{ProofOptions, prove};
    use std::sync::Arc;
    // A 14 mm square turned 30 degrees, its corners 0.005 mm inside the
    // sides of its envelope, in a 24 mm window, to that size's floor.
    let part = "(max (rotate-z (box 14 14 20) 30) (neg (box 19.1344 19.1344 22)))";
    let (lo, hi) = ([-12.0; 3], [12.0; 3]);
    let finest = ink::lang::ACCURACY / 2.0 * 24.0;
    for (shape, most) in [(part.to_string(), 64.0), (format!("(rotate-x {part} 10)"), 32.0)] {
        let (_, t) = field(&format!("(def f {shape})"), "f");
        let tape = Arc::new(t);
        let opts = ProofOptions { res: None, finest, depth: 0, tol: 0.0 };
        let at = |threads| rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap().install(|| prove(&tape, &CTX, lo, hi, &opts));
        let p = at(1);
        assert!(p.verdict == ink::query::Verdict::Proven && p.exact() && p.res > finest && p.res <= most * finest, "{shape}: {p:?}, floor {finest}");
        assert_eq!(format!("{:?}", at(4)), format!("{p:?}"));
    }
}

#[test]
fn face_on_a_grid_plane_keeps_its_area() {
    // The block's bottom face (z = 0) lies 1e-14 below a plane of the
    // sampler's grid, as rounding put it under the quadcopter's prop hub.
    // The cell above owns the face; if the plane missed the cell by
    // rounding, its area would come out 0 there and every sample would
    // weigh nothing. The face's area must not depend on where the grid's
    // planes fall.
    let (_d, t) = field("(def blk (max (- (abs x) 2) (- (abs y) 2) (- (abs (- z 1)) 1)))", "blk");
    let t = std::sync::Arc::new(t);
    let bottom = |z0: f64| -> (f64, usize) {
        let (pts, h) = ink::feature::surface(&t, &CTX, [-4.25, -3.25, z0], [3.75, 3.25, 3.0], 16.0);
        let pts: Vec<_> = pts.into_iter().filter(|p| p.1[2] < -0.9).collect();
        assert_eq!(h, 0.5);
        (pts.iter().map(|p| p.2).sum(), pts.len())
    };
    let (on, n_on) = bottom(-3.0 + 1e-14);
    let (off, n_off) = bottom(-3.2);
    assert_eq!(n_on, n_off);
    assert!(off > 12.0 && (on - off).abs() < 1e-9, "{on} on a grid plane, {off} off it");
    // The cell areas the samples weigh are exact: an axis plane through the
    // centre, h^2; a diagonal one, normal (1,1,0)/sqrt2, a rectangle h x h
    // sqrt2; x + y + z = 2.5 cuts the corner (1,1,1), a triangle of legs
    // 0.5, sqrt3/8. An axis plane on a face, or outside it by rounding:
    // the face's area; clearly outside, none.
    use ink::feature::plane_area_in_cell as area;
    use ink::v3::{scale, unit};
    let (c, z) = ([0.0; 3], [0.0, 0.0, 1.0]);
    assert!((area(c, 2.0, c, z) - 4.0).abs() < 1e-12);
    assert!((area(c, 2.0, c, unit([1.0, 1.0, 0.0])) - 4.0 * 2f64.sqrt()).abs() < 1e-12);
    assert!((area(c, 2.0, scale([1.0; 3], 2.5 / 3.0), unit([1.0; 3])) - 3f64.sqrt() / 8.0).abs() < 1e-12);
    assert!((area(c, 2.0, [0.0, 0.0, -1.0], z) - 4.0).abs() < 1e-6);
    assert!((area(c, 2.0, [0.0, 0.0, -1.0 - 1e-14], z) - 4.0).abs() < 1e-6);
    assert_eq!(area(c, 2.0, [0.0, 0.0, -1.1], z), 0.0);
}

#[test]
fn clearance_bounds_round_outward_at_any_scale() {
    use ink::query::Verdict;
    use std::sync::Arc;
    // A cube of side r = 2^-600 and the half-space x < r/4 cutting into it:
    // they overlap. The squared half-widths of the cell underflow to 0, so
    // a mean value test that read the cell as a point would prove clearance.
    let r = 2f64.powi(-600);
    let src = format!(
        "(def a (- x {q}))
(def b (max (max (- (abs (- x {h})) {h}) (- (abs (- y {h})) {h})) (- (abs (- z {h})) {h})))",
        q = r / 4.0,
        h = r / 2.0
    );
    let (d, _) = field(&src, "a");
    let tape = |n: &str| {
        let Value::Field(f, _) = d.def(n).unwrap().clone() else { panic!() };
        Arc::new(Tape::new(&d.g, f))
    };
    match clearance_of(&tape("a"), &tape("b"), &CTX, r / 16.0, [0.0; 3], [r; 3], r) {
        (Verdict::Disproven { value, .. }, Some(_)) => assert!(value < r / 16.0, "{value}"),
        v => panic!("they overlap: {v:?}"),
    }
}

#[test]
fn clearance_over_ranged_params_grows_its_boxes_outward() {
    use ink::eval::ParamSpans;
    use ink::query::Verdict;
    use std::sync::Arc;
    // B, the slab between x = H - p and H - 1/32 (H = 2^54, p in [1/8,
    // 1/4]); A, x > H + 1/16. At p = 3/16, B's point at H - 1/16 and A's
    // at H + 1/8 are 3/16 apart, under the 1/4 claimed. Over the param's
    // span, the cell's bound of A taken over the cell grown by 1/4
    // rounded to nearest, where H + 1/4 rounds to H, would leave A's points
    // outside it, and the bound prove the clearance. The second, at unit
    // scale: A from x = 1.2 - 2^-56 (above fl(1.2)), B up to 1 - 2^-60,
    // under 0.2 apart.
    let cases = [
        (
            "(def a (select (- (- x 18014398509481984) 0.0625) 1 -1)) (def b (max (+ (- x 18014398509481984) 0.03125) (- (- 18014398509481984 x) p) (neg y) (- y 2) (neg z) (- z 2)))",
            (0.125, 0.25),
            0.25,
            [18014398509481984.0 - 2.0, 0.0, 0.0],
            [18014398509481984.0, 2.0, 2.0],
            2.0,
        ),
        (
            "(def a (neg (+ (- (- x 1) 0.2) 1.3877787807814457e-17))) (def b (max (+ (- x 1) 8.673617379884035e-19) (- (- 1 x) p) (neg y) (- y 1) (neg z) (- z 1)))",
            (2f64.powi(-56), 2f64.powi(-55)),
            0.2,
            [0.0; 3],
            [1.0; 3],
            1.0,
        ),
    ];
    for (src, span, gap, lo, hi, h) in cases {
        let (d, _) = field(&format!("(param p {}) {src}", ink::eval::span_mid(span)), "a");
        let tape = |n: &str| {
            let Value::Field(f, _) = d.def(n).unwrap().clone() else { panic!() };
            Arc::new(Tape::new(&d.g, f))
        };
        let k = d.param_index("p").unwrap();
        let mut params = vec![0.0; d.params.len()];
        params[k] = ink::eval::span_mid(span);
        let ctx = Ctx { p_spans: ParamSpans::new([(k as u32, span)]), ..Ctx::new(0.0, params) };
        let (v, _) = clearance_of(&tape("a"), &tape("b"), &ctx, gap, lo, hi, h);
        assert!(!matches!(v, Verdict::Proven), "{src}: closer than {gap}, but {v:?}");
    }
}

#[test]
fn witnesses_are_confirmed_with_outward_bounds() {
    use ink::query::Verdict;
    use std::sync::Arc;
    // (x - c)^2 expanded, with c = 0.3 and fl(c c) above c^2: positive
    // everywhere in real numbers (by 3.3e-18), but within 1.4e-9 of c its f64
    // value is often negative. Sample points there are no witnesses.
    let (_d, t) = field("(def f (+ (- (* x x) (* 0.6 x)) 0.09))", "f");
    let (lo, hi) = ([0.3 - 1e-8, 0.0, 0.0], [0.3 + 1e-8, 1e-10, 1e-10]);
    let tape = Arc::new(t);
    let mut neg = false;
    for k in 0..1000 {
        neg |= eval_point(&tape, &CTX, [0.3 + (k as f64 - 500.0) * 2e-12, 0.0, 0.0]) < 0.0;
    }
    assert!(neg, "the f64 value is negative somewhere near c");
    let (v, _) = empty_at(&tape, &CTX, lo, hi, None);
    assert!(!matches!(v, Verdict::Disproven { .. }), "{v:?}");
}

/// A distance's upper bound needs a pair of points confirmed inside the two
/// solids, whatever the fields' scale. b's field is 1e-12 times a distance,
/// below 1e-10 near a, which projection alone would take for b's surface:
/// distance 0 across an 8 mm gap. A target with no inside has no pair.
#[test]
fn distance_upper_bounds_hold_at_any_field_scale() {
    let view = |b: &str| {
        let src = format!("(import \"std\")\n(body a (shape (sphere 1)))\n(body b (shape {b}))\n");
        let doc = ink::load_str(&src, "distance", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
        let ctx = doc.ctx(&[], 0.0).unwrap();
        let v = ink::view::ModelView::new(doc, ctx);
        let (fa, fb) = (v.doc.body("a").unwrap().shape, v.doc.body("b").unwrap().shape);
        (v, fa, fb)
    };
    let (mut v, fa, fb) = view("(* 1e-12 (move (sphere 1) 10 0 0))");
    let tol = v.accuracy();
    let d = ink::claims::distance(&mut v, fa, fb, ink::claims::ClaimOptions::default()).unwrap();
    let (upper, p, q) = d.upper.expect("a confirmed pair");
    assert!((8.0..8.0 + tol).contains(&upper), "{upper}");
    for (f, x) in [(fa, p), (fb, q)] {
        let t = Tape::new(&v.doc.g, f);
        assert!(ink::eval::eval_interval(&t, &CTX, x.map(ink::interval::I::point), None).hi < 0.0, "{x:?}");
    }
    // b's field touches 0 on the plane x = 5 and is positive elsewhere: no
    // point of it is inside, and the distance to it is no number.
    let (mut v, fa, fb) = view("(max (abs (- x 5)) (move (box 2 2 2) 5 0 0))");
    let Err(e) = ink::claims::distance(&mut v, fa, fb, ink::claims::ClaimOptions::default()) else { panic!("a distance to no material") };
    assert!(e.message.contains("empty"), "{e:?}");
}

/// A distance between two bodies has an answer whatever their fields read
/// away from them. b is a ball in a union with the constant 1 (the
/// catalog's servos write a missing piece so): 1, with no slope, beyond
/// 1 mm of it, where no walk from a toward it has a slope to follow. Both
/// so, the search finds no pair: again from a pair proven inside them
/// (query::inside_pair).
#[test]
fn distance_answers_where_a_field_is_flat_away_from_its_solid() {
    for a in ["(sphere 1)", "(min (sphere 1) 1)"] {
        let src = format!("(import \"std\")\n(body a (shape {a}))\n(body b (shape (min (move (sphere 1) 10 0 0) 1)))\n");
        let doc = ink::load_str(&src, "flat", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
        let ctx = doc.ctx(&[], 0.0).unwrap();
        let mut v = ink::view::ModelView::new(doc, ctx);
        let (fa, fb) = (v.doc.body("a").unwrap().shape, v.doc.body("b").unwrap().shape);
        let tol = v.accuracy();
        let d = ink::claims::distance(&mut v, fa, fb, ink::claims::ClaimOptions::default()).unwrap_or_else(|e| panic!("{a}: {e:?}"));
        let (upper, ..) = d.upper.unwrap_or_else(|| panic!("{a}: a confirmed pair"));
        assert!(d.lower <= 8.0 && 8.0 <= upper && upper - d.lower <= tol, "{a}: [{}, {upper}] accuracy {tol}", d.lower);
    }
}

/// A distance bracket holds where a field is not a distance.
/// b is a cap over the block a, 0.3 clear of its sides, cut by z < 4.8: on
/// a's top, where the cap's plate rested, b's field reads the cut plane
/// (0.2) at points 0.36 from b, and projection onto b there never settles.
/// That pair is not taken for b's surface (a 0.3 gap reported as 0).
#[test]
fn distance_brackets_the_gap_beside_a_cut_plane() {
    let src = "(import \"std\")\n(body a (shape (box 10 10 10)))\n\
               (def cap (subtract (move (box 12.6 12.6 10) 0 0 2) (move (box 10.6 10.6 9) 0 0 0.5)))\n\
               (body b (shape (max cap (- z 4.8))))\n";
    let doc = ink::load_str(src, "cap", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let (fa, fb) = (v.doc.body("a").unwrap().shape, v.doc.body("b").unwrap().shape);
    let tol = v.accuracy();
    let d = ink::claims::distance(&mut v, fa, fb, ink::claims::ClaimOptions::default()).unwrap();
    let (upper, p, q) = d.upper.expect("a confirmed pair");
    assert!(d.lower <= 0.3 && 0.3 <= upper && upper - d.lower <= tol, "[{}, {upper}] accuracy {tol}", d.lower);
    // The witnesses are inside their solids, at their surfaces.
    for (f, x) in [(fa, p), (fb, q)] {
        let r = ink::eval::eval_interval(&Tape::new(&v.doc.g, f), &CTX, x.map(ink::interval::I::point), None);
        assert!(r.hi < 0.0 && r.lo > -1e-6, "{x:?}: {r:?}");
    }
}

#[test]
fn simplified_branches_keep_real_arithmetic() {
    use ink::query::Verdict;
    use std::sync::Arc;
    // s = c = 1 + 2^-27 where x < 0: s^2 - (2 s - 1) - 2^-55 = 2^-54 - 2^-55
    // > 0 in reals (and 1 - 2^-55 where s = 2). Once a cell picked the
    // branch, its arithmetic read as constants in floats, rounding e^2 away:
    // -2^-55, a whole-cell witness.
    let (_d, t) = field("(def s (select x 1.0000000074505806 2))\n(def f (- (- (square s) (- (* 2 s) 1)) 2.7755575615628914e-17))", "f");
    let (v, _) = empty_at(&Arc::new(t), &CTX, [-1.0; 3], [1.0; 3], Some(0.25));
    assert!(!matches!(v, Verdict::Disproven { .. }), "{v:?}");
}

#[test]
fn proofs_over_regions_no_root_covers_fail_closed() {
    use ink::query::Verdict;
    use std::sync::Arc;
    // A root whose far corner overflows: no witness at [inf; 3].
    let (_d, all) = field("(def f (- 0 1))", "f");
    match empty_at(&Arc::new(all), &CTX, [1e308; 3], [1.1e308; 3], Some(5e307)).0 {
        Verdict::Disproven { witness, .. } => assert!(witness.iter().all(|x| x.is_finite()), "{witness:?}"),
        Verdict::Proven => panic!("the field is -1 everywhere"),
        Verdict::Unknown { .. } => {}
    }
}

/// The claim runner is a library over a model view: typed answers, the
/// window rule for shapes not proven bounded, and time spans (the CLI only
/// prints these).
#[test]
fn claims_run_on_a_model_view() {
    use ink::claims::{ClaimOptions, Evidence, Status, run_asserts};
    // a: a 10 mm cube; b: one sliding along x at 40 mm per unit t; c: a
    // ball 7.9 mm off a's edge (its box within 6 of a's: a clearance of 6
    // is not settled by the boxes alone).
    let src = "(import \"std\")\n(body a (shape (box 10 10 10)))\n\
               (body b (shape (move (box 10 10 10) (* 40 t) 0 0)))\n(body c (shape (move (sphere 2) 12 12 0)))\n\
               (assert clear (clearance (shape-of a) (shape-of c) 6))\n(assert near (clearance (shape-of a) (shape-of c) 9))\n\
               (assert meet (overlaps (shape-of a) (shape-of b)))\n(assert slow (value (* t t) <= 0.25))\n\
               (assert gone (empty (max x (- 600 (abs y)))))\n";
    let doc = ink::load_str(src, "claims", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let c: Vec<_> = r.asserts.iter().map(|a| &a.check).collect();
    assert!(c[0].status == Status::Holds && matches!(c[0].evidence, Evidence::ClearanceSearch { verdict: ink::query::Verdict::Proven, .. }), "{:?}", c[0]);
    assert!(c[1].status == Status::Fails && matches!(c[1].evidence, Evidence::ClearanceSearch { pair: Some(_), .. }), "{:?}", c[1]);
    assert!(c[2].status == Status::Holds && matches!(c[2].evidence, Evidence::Overlaps(..)), "{:?}", c[2]);
    assert!(c[3].status == Status::Holds && matches!(c[3].evidence, Evidence::Value(0.0)), "{:?}", c[3]);
    // Empty in the window, but a half-slab is not proven bounded: unknown.
    assert!(c[4].status == Status::Unknown && c[4].in_window == Some(Status::Holds) && c[4].unbounded, "{:?}", c[4]);
    assert_eq!(r.counts(), vec![(Status::Holds, 3), (Status::Fails, 1), (Status::Unknown, 1)]);
    assert!(!r.all_hold() && !r.all_hold_guaranteed());
    // Over t in [0, 1], t ranged: b leaves a at t = 0.25 (a moment found
    // apart), and t^2 passes 0.25.
    assert!(v.set_motion(Some((1.0, 0.0))).is_err());
    v.set_motion(Some((0.0, 1.0))).unwrap();
    let r = run_asserts(&mut v, Some("meet"), false, ClaimOptions::default()).unwrap();
    let c = &r.asserts[0].check;
    let t = |c: &ink::claims::Check| c.motion.as_ref().and_then(|m| m.motion.t).unwrap_or(f64::NAN);
    assert!(c.status == Status::Fails && matches!(&c.evidence, Evidence::Overlaps(..)) && t(c) > 0.25, "{c:?}");
    let c = &run_asserts(&mut v, Some("slow"), false, ClaimOptions::default()).unwrap().asserts[0].check;
    // A failure at a sampled moment is a witness (the samples come first),
    // as one with t ranged is.
    let past = |e: &Evidence| match e {
        Evidence::Value(x) => *x > 0.25,
        Evidence::ValueSpan(s) => s.witness.as_ref().is_some_and(|w| w.1 > 0.25),
        _ => false,
    };
    assert!(c.status == Status::Fails && c.accuracy.is_some() && past(&c.evidence) && t(c) > 0.5, "{c:?}");
    assert!(run_asserts(&mut v, Some("missing"), false, ClaimOptions::default()).is_err());
}

/// Statuses of a model's asserts, in order (a test model in a string), at
/// cell size `res` if given.
fn assert_statuses(src: &str, res: Option<f64>) -> Vec<(String, ink::claims::Check)> {
    use ink::claims::{ClaimOptions, run_asserts};
    let doc = ink::load_str(src, "claims", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    run_asserts(&mut v, None, false, ClaimOptions { res, ..ClaimOptions::default() }).unwrap().asserts.into_iter().map(|a| (a.name, a.check)).collect()
}

/// A copy is integrated in its body's frame only where its placement is
/// finite numbers over the scope: moved past the floats (e^1000 is inf,
/// at every value of a ranged param s in 1 .. 2 too), its world field has
/// no material, so its volume is not its body's.
#[test]
fn a_placement_past_the_floats_is_not_integrated_in_its_frame() {
    use ink::claims::Status::*;
    for at in ["(move (exp 1000) 0 0)", "(move 0 (- (exp 1000)) 0)", "(move 0 0 (* s (exp 1000)))"] {
        let c = assert_statuses(
            &format!("(import \"std\")\n(param s 1 (tol 0 1))\n(body b (shape (box 1 1 1)))\n(assembly root (instance c b (at {at})))\n(assert v (integral root 1 > 0.5))\n"),
            None,
        );
        assert!(c[0].1.status != Holds, "{at}: {c:?}");
    }
}

/// Integrals and weighted means over a target's bodies,
/// against exact answers: a uniform box's centre is proven inside a
/// slightly larger region and disproven outside it, its mass bracketed;
/// two densities (two bodies, or one body's density field, placed and
/// turned) move the centre of mass off the centroid, to 7.5 of 20.
/// Each part of an integral is measured from its own box: a 0.01 mm shim
/// beside a 25 mm block holds its volume (0.01) under 0.0105, where cells
/// sized by the block, sharing one cell budget, read it 0.0124.
#[test]
fn an_integral_over_parts_of_mixed_size_resolves_each() {
    let src = "(import \"std\")\n\
               (body block (shape (box 25 25 25)) (attr g 0))\n\
               (body shim (shape (box 1 1 .01)) (attr g 1))\n\
               (assembly root (instance b block) (instance s shim (at (move 300 0 0))))\n\
               (assert shim-volume (integral root (attr g) < .0105))\n";
    let got = assert_statuses(src, None);
    assert_eq!(got[0].1.status, ink::claims::Status::Holds, "{:?}", got[0].1);
}

/// A 0.3 mm layer of density (3 g) the first boxes' samples miss is in
/// their bounds, so the claim is not decided on the samples (0 g, under 1
/// g).
#[test]
fn an_integral_claim_is_not_decided_by_samples_that_miss_a_layer() {
    let src = "(import \"std\")\n\
               (body b (shape (box 100 100 100)) (attr density (select (- (abs x) .15) .001 0)))\n\
               (assert light (integral b (attr density) < 1))\n";
    let got = assert_statuses(src, None);
    assert_eq!(got[0].1.status, ink::claims::Status::Fails, "{:?}", got[0].1);
}

#[test]
fn integral_and_mean_claims() {
    use ink::claims::{Evidence, Status::*};
    // A 10 x 20 x 30 box turned 30 degrees about z, centred at (6 12 18),
    // density 0.002 g/mm^3: volume 6000, mass 12 g. Two 10 mm cubes side
    // by side along x, densities 3 and 1 (x 0.001): centre of mass x = (3 5 + 1 15) / 4 =
    // 7.5, centroid 10. The same as one 20 mm block with a density field,
    // and that block turned half a turn about z and moved 100 along x:
    // 100 - 7.5 = 92.5. A 10 mm cube with only (attr material "pla") has
    // pla's density, 0.00124: 1.24 g.
    let src = "(import \"std\")\n\
               (body blk (shape (move (rotate-z (box 10 20 30) 30) 6 12 18)) (attr density .002))\n\
               (assert x-in (mean blk x < 6.2))\n(assert x-in-above (mean blk x > 5.8))\n(assert x-out (mean blk x < 5.8))\n\
               (assert square-in (mean blk (vec (- x 6) (- 6 x) (- y 12) (- 12 y)) < 0.2))\n\
               (assert square-out (mean blk (vec (- x 6) (- 6 x) (- 12.4 y) (- y 12)) < 0.2))\n\
               (assert z-in (mean blk z (weight (attr density)) > 17.8))\n\
               (assert mass-in (integral blk (attr density) > 11.8))\n(assert mass-out (integral blk (attr density) < 11.8))\n\
               (body heavy (shape (move (box 10 10 10) 5 0 0)) (attr density .003))\n\
               (body light (shape (move (box 10 10 10) 15 0 0)) (attr density .001))\n\
               (assembly pair (instance h heavy) (instance l light))\n\
               (assert com-in (mean pair x (weight (attr density)) < 7.65))\n(assert com-out (mean pair x (weight (attr density)) < 7.35))\n\
               (assert centroid (mean pair x > 9.85))\n\
               (body graded (shape (move (box 20 10 10) 10 0 0)) (attr density (select (- x 10) .003 .001)))\n\
               (assembly turned (instance g graded (at (rot-z 180) (move 100 0 0))))\n\
               (assert graded-in (mean graded x (weight (attr density)) > 7.35))\n\
               (assert turned-in (mean turned/g x (weight (attr density)) > 92.35))\n(assert turned-out (mean turned x (weight (attr density)) > 92.65))\n\
               (body plastic (shape (box 10 10 10)) (attr material \"pla\"))\n\
               (assert pla-in (integral plastic (attr density) < 1.25))\n(assert pla-out (integral plastic (attr density) < 1.23))\n";
    let got = assert_statuses(src, None);
    let names: Vec<(&str, ink::claims::Status)> = got.iter().map(|(n, c)| (n.as_str(), c.status)).collect();
    assert_eq!(
        names,
        vec![
            ("x-in", Holds),
            ("x-in-above", Holds),
            ("x-out", Fails),
            ("square-in", Holds),
            ("square-out", Fails),
            ("z-in", Holds),
            ("mass-in", Holds),
            ("mass-out", Fails),
            ("com-in", Holds),
            ("com-out", Fails),
            ("centroid", Holds),
            ("graded-in", Holds),
            ("turned-in", Holds),
            ("turned-out", Fails),
            ("pla-in", Holds),
            ("pla-out", Fails),
        ]
    );
    // Every side of a square about the centre is claimed: a side moved to
    // 0.2 short of the centre fails, the others hold.
    let Evidence::Mean { means, .. } = &got[3].1.evidence else { panic!("{:?}", got[3]) };
    assert_eq!(means.iter().map(|m| m.1).collect::<Vec<_>>(), vec![Holds, Holds, Holds, Holds]);
    let Evidence::Mean { means, .. } = &got[4].1.evidence else { panic!("{:?}", got[4]) };
    assert_eq!(means.iter().map(|m| m.1).collect::<Vec<_>>(), vec![Holds, Holds, Fails, Holds]);
    // What decided them: estimates of the exact centre, and of the weight
    // (the volume, or mass), each to 1%; no accuracy is guaranteed.
    for (name, c) in &got {
        let exact = match name.as_str() {
            "x-in" | "x-in-above" | "x-out" => (6.0, 6000.0),
            "z-in" => (18.0, 12.0),
            "com-in" | "com-out" | "graded-in" => (7.5, 4.0),
            "centroid" => (10.0, 2000.0),
            "turned-in" | "turned-out" => (92.5, 4.0),
            _ => continue,
        };
        let Evidence::Mean { weight, means, .. } = &c.evidence else { panic!("{name}: {c:?}") };
        assert!(means.len() == 1 && means[0].1 == c.status && c.accuracy.is_none(), "{name}: {c:?}");
        assert!((means[0].0 - exact.0).abs() < 0.01 * exact.0 && (weight - exact.1).abs() < 0.01 * exact.1, "{name}: {means:?} {weight}");
    }
}

/// What an integral or a mean cannot answer, it says: a shape not proven
/// bounded (a slab) has an integral's lower bound where the integrand is
/// at least 0 everywhere, and else no bounds. A mean exists only where the
/// weight's integral is proven finite and positive and the moment's
/// finite: not over a slab (its weight, proven above 0, is infinite), nor
/// a target with no material (an empty assembly, whose integrals are
/// exactly 0, as many as the fields), nor with a weight unbounded near its
/// material (1 / x^2 over a 2 mm box from x = 0). (attr NAME) needs a body
/// that has it, and every body gives as many fields. At 5 mm cells (a 100th
/// of the window: the slab is 20 thick), so nothing is refined; the box at
/// 0.2 mm.
#[test]
fn integral_and_mean_claims_say_what_they_cannot_answer() {
    use ink::claims::{Check, Evidence, Status::*};
    let src = "(import \"std\")\n\
               (assert slab-volume (integral (slab-z 0 20) 1 > 5))\n(assert slab-mean (mean (slab-z 0 20) z < 30))\n\
               (assert nothing (mean (max (sphere 1) (move (sphere 1) 5 0 0)) x < 0))\n\
               (assembly none)\n(assert none-mass (integral none 1 < -1))\n(assert none-moments (integral none (vec x y) <= 0))\n\
               (assert none-mean (mean none x < 1))\n";
    let got = assert_statuses(src, Some(5.0));
    // With no material there is nothing to refine: a res of no size (or
    // below 0) must not refine without end.
    for res in [0.0, -1.0] {
        let empty = assert_statuses("(body b (shape 1))\n(assert a (mean b x < 1))\n(assert v (integral b 1 > 1))\n", Some(res));
        let s: Vec<_> = empty.iter().map(|(n, c)| (n.as_str(), c.status)).collect();
        assert_eq!(s, vec![("a", Unknown), ("v", Fails)], "res {res}");
    }
    // Of a shape not proven bounded nothing is estimated.
    let s: Vec<_> = got.iter().map(|(n, c)| (n.as_str(), c.status, c.unbounded)).collect();
    assert_eq!(
        s,
        vec![
            ("slab-volume", Unknown, true),
            ("slab-mean", Unknown, true),
            ("nothing", Unknown, false),
            ("none-mass", Fails, false),
            ("none-moments", Holds, false),
            ("none-mean", Unknown, false),
        ]
    );
    let note = |c: &Check, says: &str| assert!(c.note.is_some_and(|n| n.contains(says)), "{c:?}");
    for k in [2, 5] {
        note(&got[k].1, "not positive");
    }
    let integrals = |c: &Check| match &c.evidence {
        Evidence::Integral { integrals, .. } => integrals.iter().map(|(i, _)| *i).collect::<Vec<_>>(),
        e => panic!("{e:?}"),
    };
    assert_eq!(integrals(&got[3].1), vec![0.0]);
    assert_eq!(integrals(&got[4].1), vec![0.0; 2]);
    let err = |src: &str| match ink::load_str(&format!("(import \"std\")\n{src}"), "claims", Path::new(ROOT)) {
        Err(e) => format!("{e:?}"),
        Ok(_) => panic!("loaded: {src}"),
    };
    let e = err("(assert a (mean (sphere 3) x (weight (attr density)) < 1))");
    assert!(e.contains("is a field, with no attributes"), "{e}");
    let e = err("(body b (shape (sphere 3)) (attr density .002))\n(assert a (integral b (attr mass) > 0))");
    assert!(e.contains("'b' has no attr 'mass'") && e.contains("its attrs: density"), "{e}");
    // A density is g/mm^3: one in kg/m^3 or g/cm^3, unsaid, is refused.
    let e = err("(body b (shape (sphere 3)) (attr density 1270))");
    assert!(e.contains("write 1270kg/m3") && e.contains("1270g/cm3"), "{e}");
    // So is a material's, from a record.
    let dir = crate::common::Scratch::new("material-density");
    std::fs::write(dir.join("m.json"), r#"{"E": 2100, "density": 1270}"#).unwrap();
    let e = format!("{:?}", ink::load_str("(body b (shape (- x 1)) (attr material (read-json \"m.json\")))", "m", &dir).err().expect("refused"));
    assert!(e.contains("material's density 1270"), "{e}");
    let e = err("(body b (shape (sphere 3)) (attr material \"pla\"))\n(assert a (integral b (attr material) > 0))");
    assert!(e.contains("is a string, not a number or a field"), "{e}");
    let e = err("(def d (attr density))");
    assert!(e.contains("outside an integrand"), "{e}");
    let e = err("(assembly none)\n(assert a (mean none x (weight (attr density)) < 1))");
    assert!(e.contains("the target 'none' has no bodies"), "{e}");
    let pair = "(body a (shape (sphere 1)) (attr d 1))\n(body b (shape (move (sphere 1) 5 0 0)) (attr d 2))\n\
                (assembly pair (instance a a) (instance b b))\n";
    let e = err(&format!("{pair}(assert v (integral pair (if (< (attr d) 2) (vec 0) (vec 0 100)) < 10))"));
    assert!(e.contains("is 1 field for 'pair/a' but 2 fields for 'pair/b'"), "{e}");
    let e = err(&format!("{pair}(assert m (mean pair (if (< (attr d) 2) (vec x) (vec x y)) < 10))"));
    assert!(e.contains("is 1 field for 'pair/a' but 2 fields for 'pair/b'"), "{e}");
}

/// Over a ranged param or a time span, an integral claim's bounds hold at
/// every value; where they cannot decide it, the scope is halved, and a
/// failure found at a piece's middle is a witness: a box's centre moves
/// to x = 5 + off over off in [0, 4], another's to 10 t over t in [0, 1].
/// At 0.5 mm cells.
#[test]
fn integral_claims_split_their_scope() {
    use ink::claims::{ClaimOptions, Evidence, Status::*, run_asserts};
    let src = "(import \"std\")\n(param off 0 (tol 0 4))\n(body b (shape (move (box 10 10 10) (+ 5 off) 0 0)))\n\
               (assert fails-past-the-middle (mean b x < 7.5))\n\
               (body m (shape (move (box 10 10 10) (* 10 t) 0 0)))\n(assert moving-holds (mean m x < 10.5))\n";
    let doc = ink::load_str(src, "scope.ink", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let spans = v.declared_spans().unwrap();
    v.set_param_spans(&spans).unwrap();
    let sampled = |c: &ink::claims::Check| match &c.evidence {
        Evidence::Mean { sampled: Some(s), .. } => s.clone(),
        e => panic!("{e:?}"),
    };
    let opts = ClaimOptions { res: Some(0.5), ..ClaimOptions::default() };
    let run = |v: &mut ink::view::ModelView, name: &str| run_asserts(v, Some(name), false, opts).unwrap().asserts.remove(0).check;
    // Fails for off from 2.5 (the mean 7.5), not at the scope's middle
    // (off = 2): found at the first value sampled there.
    let c = run(&mut v, "fails-past-the-middle");
    let s = sampled(&c);
    assert!(c.status == Fails && c.ranged.len() == 1 && s.failing.as_ref().is_some_and(|w| w.len() == 1 && w[0].1 >= 2.5 && w[0].1 < 2.75), "{c:?}");
    // Holds at every step of a motion: measured, an estimate.
    v.set_motion(Some((0.0, 1.0))).unwrap();
    let c = run(&mut v, "moving-holds");
    assert!(c.status == Holds && c.accuracy.is_none() && c.motion.as_ref().is_some_and(|m| !m.bound), "{c:?}");
}

/// Integrals over bodies whose placements move over the scope are taken in
/// each body's own frame: its shape stays put, the motion
/// enters only the integrands, and the placement's determinant scales the
/// result. The desk lamp (examples/desk-lamp.ink) with its arm and head
/// swung about the post through any angle: its mass (2747.2 g by hand)
/// reads no swing, so one descent bounds it (in the world, each piece of
/// the swing smeared the head over its sweep: unknown after 64 pieces);
/// the long-armed lamp's centre of mass circles the post at radius 92.97,
/// past the octagon's apothem 69.29 along every normal it faces: disproven
/// at the scope's centre, swing 180. A turn of 30 degrees is floats, which
/// no rotation is: composed with the swing, a 10 mm cube's moment of z
/// about its lift of 50 (50000 at every swing) is bounded by one descent.
#[test]
fn integrals_over_moving_bodies_are_taken_in_their_frames() {
    use ink::claims::{ClaimOptions, Evidence, Status::*, run_asserts};
    let src = "(import \"std\")\n\
               (def steel 7.85e-3)\n(def aluminium 2.70e-3)\n(def pla 1.24e-3)\n\
               (body base (shape (cylinder-z 75 15)) (attr density steel))\n\
               (body post (shape (cylinder-z 8 300)) (attr density aluminium))\n\
               (fn arm-bar (length) (move (box length 16 10) (+ 8 (* 0.5 length)) 0 0))\n\
               (body arm (shape (arm-bar 300)) (attr density aluminium))\n\
               (body long-arm (shape (arm-bar 550)) (attr density aluminium))\n\
               (body head (shape (cylinder-z 40 60)) (attr density pla))\n\
               (param swing 0 (tol 0 360))\n\
               (assembly swung (instance base base (at (move 0 0 7.5))) (instance post post (at (move 0 0 165)))\n\
                 (instance arm arm (at (move 0 0 310) (rot-z swing))) (instance head head (at (move 280 0 275) (rot-z swing))))\n\
               (assembly long-swung (instance base base (at (move 0 0 7.5))) (instance post post (at (move 0 0 165)))\n\
                 (instance arm long-arm (at (move 0 0 310) (rot-z swing))) (instance head head (at (move 530 0 275) (rot-z swing))))\n\
               (fn along (k) (+ (* (cos (deg->rad (* 45 k))) x) (* (sin (deg->rad (* 45 k))) y)))\n\
               (assert mass (integral swung (attr density) < 3000))\n\
               (assert long-stands (mean long-swung (map along (range 8)) (weight (attr density)) < (* 75 (cos (deg->rad 22.5)))))\n\
               (body cube (shape (box 10 10 10)))\n\
               (assembly tilted (instance c cube (at (rot-z swing) (rot-x 30) (move 0 0 50))))\n\
               (assert tilted-z (integral tilted z > 45000))\n";
    let doc = ink::load_str(src, "frames.ink", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let spans = v.declared_spans().unwrap();
    v.set_param_spans(&spans).unwrap();
    let got: Vec<(String, ink::claims::Check)> = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts.into_iter().map(|a| (a.name, a.check)).collect();
    let integral = |c: &ink::claims::Check| match &c.evidence {
        Evidence::Integral { integrals, sampled, .. } => (integrals[0].0, sampled.clone()),
        e => panic!("{e:?}"),
    };
    let statuses: Vec<(&str, ink::claims::Status)> = got.iter().map(|(n, c)| (n.as_str(), c.status)).collect();
    assert_eq!(statuses, vec![("mass", Holds), ("long-stands", Fails), ("tilted-z", Holds)]);
    // Taken in the bodies' frames, what they read does not read the swing:
    // one estimate each, about the exact values.
    let (mass, sampled) = integral(&got[0].1);
    assert!(sampled.is_none() && (mass - 2747.2).abs() < 0.01 * 2747.2, "{mass:?} {sampled:?}");
    let (z, sampled) = integral(&got[2].1);
    assert!((z - 50000.0).abs() < 0.01 * 50000.0, "{z:?} {sampled:?}");
    let Evidence::Mean { sampled: Some(s), .. } = &got[1].1.evidence else { panic!("{:?}", got[1]) };
    assert!(s.failing.is_some(), "{s:?}");
}

/// A placement in a body's frame reads what the evaluator reads. An
/// expression of constants and fixed params is its f64 value: (- (+ p 1) p)
/// with p = 1e16 is 0, so a unit cube moved up by it, and swung, has its
/// centre at z = 0 (read as real arithmetic, 1). An exact turn adds no
/// error: half a cube moved along x by a ranged p lies below z = 0, and the
/// jump there (an unbounded slope) still bounds it to 0.5.
#[test]
fn integrals_in_frames_read_placements_as_the_evaluator_does() {
    use ink::claims::{ClaimOptions, Evidence, Status::*, run_asserts};
    let run = |src: &str| {
        let doc = ink::load_str(src, "placed.ink", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
        let ctx = doc.ctx(&[], 0.0).unwrap();
        let mut v = ink::view::ModelView::new(doc, ctx);
        let spans = v.declared_spans().unwrap();
        v.set_param_spans(&spans).unwrap();
        run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts.into_iter().map(|a| (a.name, a.check)).collect::<Vec<_>>()
    };
    let statuses = |got: &[(String, ink::claims::Check)]| got.iter().map(|(n, c)| (n.clone(), c.status)).collect::<Vec<_>>();
    let named = |v: Vec<(&str, ink::claims::Status)>| v.into_iter().map(|(n, s)| (n.to_string(), s)).collect::<Vec<_>>();
    let got = run("(import \"std\")\n(param p 1e16)\n(param angle 0 (tol 0 360))\n(body cube (shape (box 1 1 1)))\n\
                   (assembly root (instance c cube (at (rot-z angle) (move 0 0 (- (+ p 1) p)))))\n\
                   (assert lifted (mean root z > 0.5))\n(assert at-zero (mean root z < 0.1))\n");
    assert_eq!(statuses(&got), named(vec![("lifted", Fails), ("at-zero", Holds)]), "{got:?}");
    let got = run("(import \"std\")\n(param p 0 (tol 0 1))\n(body cube (shape (box 1 1 1)))\n\
                   (assembly root (instance c cube (at (move p 0 0))))\n(assert half (integral root (select z 1 0) < 2))\n");
    let Evidence::Integral { integrals, .. } = &got[0].1.evidence else { panic!("{got:?}") };
    assert!(got[0].1.status == Holds && (integrals[0].0 - 0.5).abs() < 0.01, "{got:?}");
}

/// An integral claim stops where the bounds proven so far decide it, and
/// reports the estimate there: a sphere's volume (4189 mm^3) claimed under
/// 10000 is decided on coarser boxes than one claimed under 4250, which
/// goes on to the mass's 2%; both hold, with the estimate as the value.
#[test]
fn an_integral_claim_stops_where_its_bounds_decide_it() {
    use ink::claims::{ClaimOptions, Evidence, Status::Holds, run_asserts};
    let src = "(import \"std\")\n(body b (shape (sphere 10)))\n(assert far (integral b 1 < 10000))\n(assert near (integral b 1 < 4250))\n";
    let doc = ink::load_str(src, "stops.ink", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let got = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts;
    let res: Vec<f64> = got
        .iter()
        .map(|a| match &a.check.evidence {
            Evidence::Integral { integrals, res, .. } => {
                assert!(a.check.status == Holds && (integrals[0].0 / 4188.79 - 1.0).abs() < 0.02, "{}: {:?}", a.name, a.check);
                *res
            }
            e => panic!("{e:?}"),
        })
        .collect();
    assert!(res[0] > res[1], "far stopped at boxes of {}, near at {}", res[0], res[1]);
}

/// Each assert knows the file and line it is written in, so the runner can
/// check a file's own claims and leave its imports' to their files.
#[test]
fn asserts_know_their_file() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let dir = crate::common::Scratch::new("own");
    std::fs::write(dir.join("base.ink"), "(assert false-there (value 1 < 0))\n").unwrap();
    let src = "(import \"base.ink\" as base)\n\n(assert true-here (value 1 > 0))\n";
    let doc = ink::load_str(src, "top.ink", &dir).unwrap_or_else(|e| panic!("{e:?}"));
    let at: Vec<_> = doc.asserts.iter().map(|a| (a.name.as_str(), doc.is_own(a), a.line)).collect();
    assert_eq!(at, vec![("base.false-there", false, 1), ("true-here", true, 3)]);
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, true, ClaimOptions::default()).unwrap();
    assert!(r.asserts.len() == 1 && r.asserts[0].name == "true-here" && r.asserts[0].check.status == Status::Holds);
    assert_eq!(r.asserts[0].source, "top.ink:3");
    assert_eq!(r.imported_not_run, Some(1));
    assert!(run_asserts(&mut v, Some("base.false-there"), true, ClaimOptions::default()).is_err());
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    assert!(r.asserts.len() == 2 && !r.all_hold() && r.imported_not_run.is_none());
}

/// An ASCII STL is read strictly: a vertex that is not three numbers, a
/// facet cut short, or a file ending inside one is an error at its line
/// (a dropped vertex would shift every triangle after it); a mesh of only
/// degenerate triangles is an error, not a panic.
#[test]
fn stl_input_is_read_strictly() {
    let facet = |v: [&str; 3]| format!("facet normal 0 0 0\n outer loop\n  vertex {}\n  vertex {}\n  vertex {}\n endloop\nendfacet\n", v[0], v[1], v[2]);
    let tet = |v: &str| {
        format!(
            "solid t\n{}{}{}{}endsolid t\n",
            facet(["0 0 0", "0 1 0", "1 0 0"]),
            facet(["0 0 0", "1 0 0", v]),
            facet(["0 0 0", "0 0 1", "0 1 0"]),
            facet(["1 0 0", "0 1 0", "0 0 1"])
        )
    };
    let read = |text: &str| ink::grid::read_stl(text.as_bytes());
    assert_eq!(read(&tet("0 0 1")).unwrap().len(), 4);
    assert_eq!(read(&tet("0 0 1,0")).unwrap_err(), "STL line 13: expected 'vertex x y z', three numbers, got 'vertex 0 0 1,0'");
    assert!(read(&tet("5mm 0 0 1")).unwrap_err().starts_with("STL line 13:"));
    let cut = tet("0 0 1").lines().take(26).collect::<Vec<_>>().join("\n");
    assert_eq!(read(&cut).unwrap_err(), "STL ends inside a facet at line 26: the file is truncated");
    for (text, line) in [(tet("0 0 1").replace("facet normal 0 0 0", "facet garbage"), "STL line 2:"), (tet("0 0 1").replace("endloop", "endloop garbage"), "STL line 7:")] {
        assert!(read(&text).unwrap_err().starts_with(line));
    }
    assert!(read("solid empty\nendsolid empty\n").unwrap().is_empty());
    let mut binary = vec![b' '; 84 + 50];
    binary[..5].copy_from_slice(b"solid");
    binary[80..84].copy_from_slice(&1_u32.to_le_bytes());
    binary[84..].fill(0);
    assert_eq!(ink::grid::read_stl(&binary).unwrap().len(), 1, "a solid binary header stays binary");
    let dir = crate::common::Scratch::new("stl-strict");
    std::fs::write(dir.join("point.stl"), format!("solid p\n{}endsolid p\n", facet(["1 1 1", "1 1 1", "1 1 1"]))).unwrap();
    let e = ink::load_str("(def m (import-stl \"point.stl\" (res 0.5)))", "t", &dir).err().unwrap();
    assert!(e.message.contains("STL has no triangles of nonzero area (1 degenerate facet)"), "{e:?}");
    // An import inside a function reads the mesh beside its definition,
    // even when the caller has a degenerate mesh of the same name.
    std::fs::create_dir_all(dir.join("lib")).unwrap();
    std::fs::write(dir.join("lib/point.stl"), tet("0 0 1")).unwrap();
    std::fs::write(dir.join("lib/part.ink"), r#"(fn part () (import-stl "point.stl" (res 0.5)))"#).unwrap();
    let d = ink::load_str(r#"(import "lib/part.ink") (def m (part))"#, "t", &dir).unwrap();
    assert_eq!(d.g.grids.len(), 1);
}

/// A binary STL of these triangles.
fn write_stl(path: &Path, tris: &[[[f32; 3]; 3]]) {
    let mut b = vec![0u8; 80];
    b.extend((tris.len() as u32).to_le_bytes());
    for t in tris {
        b.extend([0u8; 12]);
        t.iter().flatten().for_each(|c| b.extend(c.to_le_bytes()));
        b.extend([0u8; 2]);
    }
    std::fs::write(path, b).unwrap();
}

/// The cube 0..10, outward-facing.
fn cube_tris() -> Vec<[[f32; 3]; 3]> {
    box_tris([0.0; 3], [10.0; 3])
}

/// The box lo..hi, outward-facing.
fn box_tris(lo: [f32; 3], hi: [f32; 3]) -> Vec<[[f32; 3]; 3]> {
    let v = |i: usize| -> [f32; 3] { std::array::from_fn(|a| if (i >> a) & 1 == 1 { hi[a] } else { lo[a] }) };
    [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]].iter().flat_map(|q| [[v(q[0]), v(q[1]), v(q[2])], [v(q[0]), v(q[2]), v(q[3])]]).collect()
}

/// The boxes -2..0 and 0..2 by -1..1 by -1..1, outward-facing, their
/// shared face split along one diagonal in the first and the other in the
/// second.
fn adjacent_cubes() -> Vec<[[f32; 3]; 3]> {
    let mut t = box_tris([-2.0, -1.0, -1.0], [0.0, 1.0, 1.0]);
    t.extend(box_tris([0.0, -1.0, -1.0], [2.0, 1.0, 1.0]).into_iter().filter(|t| t.iter().any(|v| v[0] != 0.0)));
    let v = |y: f32, z: f32| [0.0, y, z];
    t.extend([[v(-1.0, -1.0), v(1.0, -1.0), v(-1.0, 1.0)], [v(1.0, -1.0), v(1.0, 1.0), v(-1.0, 1.0)]]);
    t
}

/// The signed distance to the box lo..hi, in floats.
fn box_distance(p: [f64; 3], lo: [f64; 3], hi: [f64; 3]) -> f64 {
    let q: [f64; 3] = std::array::from_fn(|a| (lo[a] - p[a]).max(p[a] - hi[a]));
    let out = q.iter().map(|v| v.max(0.0).powi(2)).sum::<f64>().sqrt();
    out + q[0].max(q[1]).max(q[2]).min(0.0)
}

/// A closed mesh is sampled as the signed distance to the solid it
/// encloses, sides by parity: whichever way its facets face, a turned
/// solid's slanted faces, and two bodies sharing a face (parity cancels the
/// face's two sheets, the distance still counts them).
#[test]
fn a_meshs_samples_are_its_solids_signed_distance() {
    use ink::grid::MeshDistance;
    let tris = |t: &[[[f32; 3]; 3]]| -> Vec<[[f64; 3]; 3]> { t.iter().map(|t| t.map(|v| v.map(|c| c as f64))).collect() };
    let mut flipped = box_tris([-1.0; 3], [1.0; 3]);
    flipped.iter_mut().for_each(|t| t.swap(1, 2));
    // An octahedron turned off the axes: inside where |u| + |v| + |w| < 1.
    let (c, s) = (0.8f64, 0.6f64);
    let rot = |q: [f64; 3]| [c * q[0] - s * q[1], s * c * q[0] + c * c * q[1] - s * q[2], s * s * q[0] + s * c * q[1] + c * q[2]];
    let unrot = |q: [f64; 3]| [c * q[0] + s * c * q[1] + s * s * q[2], -s * q[0] + c * c * q[1] + s * c * q[2], -s * q[1] + c * q[2]];
    let ax = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut oct = Vec::new();
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                let v = |a: usize, sg: f64| rot(ax[a].map(|x: f64| x * sg));
                oct.push([v(0, sx), v(1, sy), v(2, sz)]);
            }
        }
    }
    let inside_box = |lo: [f64; 3], hi: [f64; 3]| move |p: [f64; 3]| box_distance(p, lo, hi) < 0.0;
    let cases: Vec<(&str, Vec<[[f64; 3]; 3]>, Box<dyn Fn([f64; 3]) -> bool>)> = vec![
        ("cube", tris(&box_tris([-1.0; 3], [1.0; 3])), Box::new(inside_box([-1.0; 3], [1.0; 3]))),
        ("flipped", tris(&flipped), Box::new(inside_box([-1.0; 3], [1.0; 3]))),
        ("octahedron", oct, Box::new(move |p: [f64; 3]| unrot(p).iter().map(|x| x.abs()).sum::<f64>() < 1.0)),
        ("adjacent", tris(&adjacent_cubes()), Box::new(inside_box([-2.0, -1.0, -1.0], [2.0, 1.0, 1.0]))),
    ];
    for (name, t, inside) in cases {
        let md = MeshDistance::new(&t);
        assert_eq!(md.odd_edges(), 0, "{name}");
        let (h, origin, dims) = (0.07, [-1.31, -1.29, -1.33], [38, 38, 38]);
        let (values, unclear) = md.sample(origin, h, dims);
        assert_eq!(unclear, 0, "{name}");
        for (n, v) in values.iter().enumerate() {
            let q = ink::grid::ijk(n, dims);
            let p: [f64; 3] = std::array::from_fn(|a| origin[a] + h * q[a] as f64);
            let d = t.iter().map(|t| tri_dist(p, t)).fold(f64::INFINITY, f64::min);
            let d = if inside(p) { -d } else { d };
            assert!((v - d).abs() <= 1e-9, "{name}: sample {v} at {p:?}, distance {d}");
        }
    }
    // Closed: every edge in an even number of triangles (a doubled facet
    // keeps it so); a missing facet does not.
    let mut t = box_tris([-1.0; 3], [1.0; 3]);
    assert!(MeshDistance::new(&tris(&t[..11])).odd_edges() > 0);
    t.extend([t[0], t[0]]);
    assert_eq!(MeshDistance::new(&tris(&t)).odd_edges(), 0);
}

/// The distance from p to triangle t, from its plane and its edges.
fn tri_dist(p: [f64; 3], t: &[[f64; 3]; 3]) -> f64 {
    let sub = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] { std::array::from_fn(|k| a[k] - b[k]) };
    let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
    let inside = (0..3).all(|k| dot(n, cross(sub(t[k], p), sub(t[(k + 1) % 3], p))) >= 0.0);
    let edges = (0..3).map(|k| seg_dist(p, t[k], t[(k + 1) % 3])).fold(f64::INFINITY, f64::min);
    if inside { (dot(n, sub(t[0], p)) / dot(n, n).sqrt()).abs() } else { edges }
}

/// Claims on an imported mesh are measured on its samples: they hold or
/// fail as estimates (no accuracy), and name what was sampled, at what spacing, with the
/// error measured (no bound: sampled in floats).
#[test]
fn claims_on_an_import_are_measured_on_its_samples() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let dir = crate::common::Scratch::new("import");
    let cube = cube_tris();
    write_stl(&dir.join("cube.stl"), &cube);
    write_stl(&dir.join("open.stl"), &cube[..11]);
    let src = "(import \"std\")\n(body c (shape (import-stl \"cube.stl\" (res 0.5))))\n\
               (body f (shape (import-stl \"cube.stl\" (res 0.1) (region (8 0 0) (12 10 10)))))\n\
               (body o (shape (import-stl \"open.stl\" (res 0.5))))\n\
               (def block (box 2 2 2))\n\
               (assert apart (disjoint (shape-of c) (move block 12 5 5)))\n\
               (assert apart-0.2 (disjoint (shape-of c) (move block 11.2 5 5)))\n\
               (assert hit (disjoint (shape-of c) (move block 9 5 5)))\n\
               (assert hit-0.2 (disjoint (shape-of c) (move block 10.8 5 5)))\n\
               (assert deep (value (warp (shape-of c) 5 5 5) < -4))\n\
               (assert fine-apart-0.2 (disjoint (shape-of f) (move block 11.2 5 5)))\n\
               (assert fine-cut (disjoint (shape-of f) (move block 5 5 5)))\n\
               (assert open-apart (disjoint (shape-of o) (move block 12 5 5)))\n\
               (assert plain (disjoint block (move block 30 0 0)))\n";
    let doc = ink::load_str(src, "import.ink", &dir).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    // Each status, and whether it is an estimate (no accuracy).
    let got: Vec<(&str, Status, bool)> = r.asserts.iter().map(|a| (a.name.as_str(), a.check.status, a.check.asked.is_none())).collect();
    use Status::*;
    assert_eq!(
        got,
        vec![
            ("apart", Holds, true),
            ("apart-0.2", Holds, true),
            ("hit", Fails, true),
            ("hit-0.2", Fails, true),
            ("deep", Holds, true),
            ("fine-apart-0.2", Holds, true),
            ("fine-cut", Holds, true),
            ("open-apart", Holds, true),
            ("plain", Holds, false),
        ]
    );
    // What each read: the mesh, its spacing, and the deviation measured at
    // probes, below what a 1-Lipschitz field's samples can be off by.
    let s = &r.asserts[0].check.inexact.samples;
    assert!(s.len() == 1 && s[0].what == "cube.stl" && s[0].kind == "mesh" && s[0].res == 0.5 && s[0].band.is_none(), "{s:?}");
    let (measured, probes) = s[0].measured.unwrap();
    assert!(measured <= 0.5 * 3f64.sqrt() / 2.0 && measured > 0.05 && probes > 1000, "{s:?}");
    let o = &r.asserts[7].check.inexact.samples[0];
    assert!(o.note.as_deref().is_some_and(|n| n.contains("not closed")), "{o:?}");
    assert!(r.asserts[8].check.inexact.samples.is_empty());
}

// ---------------------------------------------------------------- ranged params (tolerances)

/// A model, viewed with these params ranged over these bands (each param's
/// value at the middle of its band).
fn ranged_view(src: &str, bands: &[(&str, f64, f64)]) -> ink::view::ModelView {
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let mut ctx = d.ctx(&[], 0.0).unwrap();
    ctx.p_spans = ink::eval::ParamSpans::new(bands.iter().map(|(n, lo, hi)| (d.param_index(n).unwrap() as u32, (*lo, *hi))));
    ink::view::ModelView::new(d, ctx.into_scope())
}

fn def_field(v: &ink::view::ModelView, name: &str) -> u32 {
    let Value::Field(f, _) = v.doc.def(name).unwrap_or_else(|| panic!("no def {name}")).clone() else { panic!("{name} is not a field") };
    f
}

/// The disproof's witness point and params (by name), from its evidence.
fn witness_of(v: &ink::view::ModelView, c: &ink::claims::Check) -> ([f64; 3], Vec<(String, f64)>) {
    use ink::claims::Evidence;
    use ink::query::Verdict;
    let (w, ps) = match &c.evidence {
        Evidence::Empty(p) | Evidence::Overlaps(p, _) => match &p.verdict {
            Verdict::Disproven { witness, params, .. } => (*witness, params.clone()),
            o => panic!("{o:?}"),
        },
        Evidence::ClearanceSearch { verdict: Verdict::Disproven { witness, params, .. }, .. } => (*witness, params.clone()),
        e => panic!("{e:?}"),
    };
    (w, ps.into_iter().map(|(k, x)| (v.doc.params[k as usize].name.clone(), x)).collect())
}

/// A field's value at a point with the params at these values (unranged).
fn value_at(v: &ink::view::ModelView, f: u32, p: [f64; 3], params: &[(String, f64)]) -> f64 {
    let mut ctx = Ctx::new(v.ctx.t, v.ctx.params.clone());
    for (n, x) in params {
        ctx.params[v.doc.param_index(n).unwrap()] = *x;
    }
    eval_point(&Tape::new(&v.doc.g, f), &ctx, p)
}

/// A claim checked exactly (accuracy 0): these tests are of the ranged
/// params' witnesses and bounds, whatever the accuracy.
fn check(v: &mut ink::view::ModelView, c: ink::lang::Claim) -> ink::claims::Check {
    ink::claims::check(v, c, None, ink::claims::ClaimOptions { accuracy: Some(0.0), ..Default::default() }).unwrap()
}

const FIT: &str = "(param shaft-d 8) (param bore-d 8.2)
(def shaft (cylinder-z (/ shaft-d 2) 10))
(def sleeve (subtract (cylinder-z 6 14) (cylinder-z (/ bore-d 2) 16)))";

#[test]
fn a_fit_is_proven_over_its_tolerances_or_fails_inside_them() {
    use ink::claims::Status;
    use ink::lang::Claim;
    // Shaft 7.95..8.05 in a bore 8.15..8.25: the radial gap is at least
    // (8.15 - 8.05) / 2 = 0.05 for every pair of sizes.
    let mut v = ranged_view(FIT, &[("shaft-d", 7.95, 8.05), ("bore-d", 8.15, 8.25)]);
    let (shaft, sleeve) = (def_field(&v, "shaft"), def_field(&v, "sleeve"));
    assert_eq!(check(&mut v, Claim::Clearance(shaft, sleeve, 0.04)).status, Status::Holds);
    let both = v.doc.g.max(shaft, sleeve);
    assert_eq!(check(&mut v, Claim::Empty(both)).status, Status::Holds);
    // A bore 8.00..8.10: the gap is 0.025 at the middle sizes, but a 8.05
    // shaft in an 8.00 bore interferes by 0.025. A witness: params in their
    // bands, the point inside both.
    let mut v = ranged_view(FIT, &[("shaft-d", 7.95, 8.05), ("bore-d", 8.0, 8.1)]);
    let both = v.doc.g.max(shaft, sleeve);
    let c = check(&mut v, Claim::Empty(both));
    assert_eq!(c.status, Status::Fails, "{c:?}");
    let (w, ps) = witness_of(&v, &c);
    let get = |n: &str| ps.iter().find(|p| p.0 == n).unwrap().1;
    assert!((7.95..=8.05).contains(&get("shaft-d")) && (8.0..=8.1).contains(&get("bore-d")), "{ps:?}");
    assert!(value_at(&v, shaft, w, &ps) < 0.0 && value_at(&v, sleeve, w, &ps) < 0.0, "{w:?} {ps:?}");
    // The same clearance claim fails too, with a pair at one pair of sizes.
    let c = check(&mut v, Claim::Clearance(shaft, sleeve, 0.01));
    assert_eq!(c.status, Status::Fails, "{c:?}");
    // At the values alone (8 in 8.2) it would have passed.
    let mut nominal = ranged_view(FIT, &[]);
    let both = nominal.doc.g.max(shaft, sleeve);
    assert_eq!(check(&mut nominal, Claim::Empty(both)).status, Status::Holds);
    // To an accuracy: an interference of 0.025 is within 0.05, beyond 0.01,
    // whose witness's cube of half-side 0.01 is inside both at its sizes.
    let both = v.doc.g.max(shaft, sleeve);
    let to = |v: &mut ink::view::ModelView, t: f64| {
        ink::claims::check(v, Claim::Disjoint(shaft, sleeve), None, ink::claims::ClaimOptions { accuracy: Some(t), ..Default::default() }).unwrap()
    };
    assert_eq!(to(&mut v, 0.05).status, Status::Holds);
    let c = to(&mut v, 0.01);
    assert_eq!(c.status, Status::Fails, "{c:?}");
    let (w, ps) = witness_of(&v, &c);
    let mut at = v.ctx.clone();
    for (n, x) in &ps {
        at.params[v.doc.param_index(n).unwrap()] = *x;
    }
    let cube: [ink::interval::I; 3] = std::array::from_fn(|k| ink::interval::I::new(w[k] - 0.01, w[k] + 0.01));
    assert!(ink::eval::eval_interval(&Tape::new(&v.doc.g, both), &at.point(), cube, None).hi < 0.0, "{w:?} {ps:?}");
}

#[test]
fn ranged_values_fail_inside_their_ranges_not_only_at_corners() {
    use ink::lang::Cmp;
    use ink::ranged::{Status, value_span};
    let v = ranged_view("(param p 0) (def sq (square p)) (def lifted (+ (square p) 0.3)) (def pp (* p p))", &[("p", -1.0, 1.0)]);
    let opts = v.proof_options(None, 0.0);
    let tape = |n: &str| Tape::new(&v.doc.g, def_field(&v, n));
    // p^2 >= 0.25 holds at both ends of [-1, 1], and fails at p = 0.
    let r = value_span(&tape("sq"), &v.ctx, [0.0; 3], ink::ranged::Comparison::exact(Cmp::Ge, 0.25), &opts);
    assert!(r.status == Status::Disproven && r.witness.as_ref().is_some_and(|(at, x)| at.params[0].1.abs() < 0.5 && *x < 0.25), "{r:?}");
    assert_eq!(value_span(&tape("lifted"), &v.ctx, [0.0; 3], ink::ranged::Comparison::exact(Cmp::Ge, 0.25), &opts).status, Status::Proven);
    // p p < 4 on [1, 3]: the same param twice is the same value (9 at p = 3).
    let w = ranged_view("(param p 2) (def pp (* p p))", &[("p", 1.0, 3.0)]);
    let r = value_span(&Tape::new(&w.doc.g, def_field(&w, "pp")), &w.ctx, [0.0; 3], ink::ranged::Comparison::exact(Cmp::Lt, 4.0), &opts);
    assert_eq!(r.status, Status::Disproven, "{r:?}");
    // 2 p - p - p is 0 at every p: first order in p its rate is 0 and the
    // whole band is one piece. Plain bounds over a piece of width w give
    // +-2 w, and halving p 8 times leaves it unknown.
    let w = ranged_view("(param p 0) (def zero (- (* 2 p) (+ p p)))", &[("p", 0.0, 1.0)]);
    let shallow = ink::ranged::ProofOptions { depth: 2, ..opts };
    let r = value_span(&Tape::new(&w.doc.g, def_field(&w, "zero")), &w.ctx, [0.0; 3], ink::ranged::Comparison::exact(Cmp::Lt, 1e-9), &shallow);
    assert!(r.status == Status::Proven && r.range.1 < 1e-12, "{r:?}");
}

#[test]
fn a_ranged_value_meeting_its_bound_along_a_curve_ends_at_its_first_undecided_piece() {
    // a a - 2 a b + b b >= 0 holds for every a, b in [-1, 1], meeting 0 all
    // along a = b, where bounds over a piece reach below 0 at any size:
    // asked exactly, unknown at the first piece no halving decides (each
    // param halved 32 times), not after refining the whole line (a million
    // pieces). To 0.001 it holds; moved down by 0.0001 it fails off the line.
    use ink::lang::Cmp;
    use ink::ranged::{Comparison, Status, value_span};
    let v = ranged_view("(param a 0) (param b 0) (def sq (+ (* a a) (* -2 a b) (* b b))) (def low (- sq 0.0001))", &[("a", -1.0, 1.0), ("b", -1.0, 1.0)]);
    let opts = v.proof_options(None, 0.0);
    let tape = |n: &str| Tape::new(&v.doc.g, def_field(&v, n));
    let r = value_span(&tape("sq"), &v.ctx, [0.0; 3], Comparison::exact(Cmp::Ge, 0.0), &opts);
    assert!(r.status == Status::Unknown && r.pieces < 100_000, "{r:?}");
    let r = value_span(&tape("sq"), &v.ctx, [0.0; 3], Comparison { cmp: Cmp::Ge, bound: 0.0, tol: 0.001 }, &opts);
    assert_eq!(r.status, Status::Proven, "{r:?}");
    let r = value_span(&tape("low"), &v.ctx, [0.0; 3], Comparison::exact(Cmp::Ge, 0.0), &opts);
    assert!(r.status == Status::Disproven && r.witness.as_ref().is_some_and(|(_, x)| *x < 0.0), "{r:?}");
    // A param the value does not read, narrowed to one value (as a
    // middle's search narrows its input), is no axis: the search is the
    // same (its share of the band, 0, made every piece tie, searched in the
    // order made, 4,163 pieces then 73,283).
    let w = ranged_view("(param a 0) (param b 0) (param c 0) (def sq (+ (* a a) (* -2 a b) (* b b)))", &[("a", -1.0, 1.0), ("b", -1.0, 1.0), ("c", 0.0, 1.0)]);
    let sq = Tape::new(&w.doc.g, def_field(&w, "sq"));
    let k = w.doc.param_index("c").unwrap() as u32;
    let exact = Comparison::exact(Cmp::Ge, 0.0);
    let (r, at) = (value_span(&sq, &w.ctx, [0.0; 3], exact, &opts), value_span(&sq, &w.ctx.with_param(k, (0.5, 0.5)), [0.0; 3], exact, &opts));
    assert!(r.status == Status::Unknown && at.status == Status::Unknown && at.pieces == r.pieces, "{r:?} {at:?}");
}

#[test]
fn a_value_solved_over_tolerances_halves_the_param_that_narrows_it() {
    // sx = sqrt(L^2 - h^2) over L 50 +-0.1, h 30 +-0.1 reaches 40.2015:
    // undecided against 40.2 over the bands, decided once L or h is halved
    // where it narrows sx most. Crediting the solved value's blur to both
    // would halve the wrong one 32 times over (9.5 million pieces).
    use ink::ranged::{Comparison, Status, value_span};
    let src = "(param L 50 (tol -0.1 0.1)) (param h 30 (tol -0.1 0.1))
(solve place (unknowns (theta 37) (sx 40))
  (== (* L (cos (* theta (/ pi 180)))) sx)
  (== (* L (sin (* theta (/ pi 180)))) h))";
    let mut v = toleranced(src).unwrap();
    let sx = v.doc.param_index("sx").unwrap();
    let f = v.doc.g.var(ink::graph::Var::Param(sx as u32));
    let r = value_span(&Tape::new(&v.doc.g, f), &v.ctx, [0.0; 3], Comparison::exact(ink::lang::Cmp::Le, 40.2), &v.proof_options(None, 0.0));
    assert!(r.status == Status::Proven && r.pieces < 100, "{r:?}");
}

/// p p - (2 p - 1) - 2^-55 for p within 2^-30 of 1 + 2^-27: in real
/// numbers (p - 1)^2 - 2^-55 > 0, but in floats at the middle p p rounds
/// to 2 p - 1, and the value is -2^-55.
const TRAP: &str = "(param p 1)
(def trap (- (- (* p p) (- (* 2 p) 1)) 2.7755575615628914e-17))
(def cube (max (box 2 2 2) trap))";

fn trap_view() -> ink::view::ModelView {
    let c = 1.0 + 2f64.powi(-27);
    ranged_view(TRAP, &[("p", c - 2f64.powi(-30), c + 2f64.powi(-30))])
}

#[test]
fn a_value_witness_means_what_its_proof_means() {
    // trap > 0 holds for every p in the band (in real numbers). Read in
    // floats at the band's middle it would be -2^-55 and "violate" it.
    use ink::ranged::{Status, value_span};
    let v = trap_view();
    let opts = ink::ranged::ProofOptions { depth: 2, ..v.proof_options(None, 0.0) };
    let r = value_span(&Tape::new(&v.doc.g, def_field(&v, "trap")), &v.ctx, [0.0; 3], ink::ranged::Comparison::exact(ink::lang::Cmp::Gt, 0.0), &opts);
    assert!(r.status != Status::Disproven && r.witness.is_none(), "{r:?}");
}

#[test]
fn a_witness_point_means_what_its_proof_means() {
    // A cube where trap < 0: empty for every p in the band (in real
    // numbers); read in floats at the band's middle the whole cube would be
    // a witness.
    use ink::claims::Status;
    let mut v = trap_view();
    let cube = def_field(&v, "cube");
    let c = ink::claims::check(&mut v, ink::lang::Claim::Empty(cube), None, ink::claims::ClaimOptions { res: Some(0.5), ..Default::default() }).unwrap();
    assert!(c.status != Status::Fails, "{c:?}");
}

#[test]
fn a_part_moved_by_a_ranged_param_is_found_wherever_it_goes() {
    // A ball moved p along x, p in [0, 1000]: its reach is the whole
    // family's, not the box at one value; a ball at 500 meets it.
    use ink::claims::Status;
    let mut v = ranged_view("(param p 0) (def ball (move (sphere 1) p 0 0)) (def far (move (sphere 1) 500 0 0))", &[("p", 0.0, 1000.0)]);
    let (ball, far) = (def_field(&v, "ball"), def_field(&v, "far"));
    let r = v.reach(ball);
    assert!(r.whole && r.hi[0] >= 1001.0 && r.lo[0] <= -1.0, "{r:?}");
    let both = v.doc.g.max(ball, far);
    let c = check(&mut v, ink::lang::Claim::Empty(both));
    assert_eq!(c.status, Status::Fails, "{c:?}");
    let (w, ps) = witness_of(&v, &c);
    assert!((ps[0].1 - 500.0).abs() < 2.0 && (w[0] - 500.0).abs() < 1.0, "{w:?} {ps:?}");
}

/// An overlap over ranged params is not proven by a point found at one
/// value of them: unknown, wherever they overlap; apart at every value, it
/// fails.
#[test]
fn overlaps_over_ranged_params_are_not_proven_at_one_value() {
    // Two unit balls, the second p along x: they overlap for p < 2 only.
    use ink::claims::Status;
    let src = "(param p 0) (def a (sphere 1)) (def b (move (sphere 1) p 0 0))";
    for ((lo, hi), want) in [((0.0, 3.0), Status::Unknown), ((0.0, 1.0), Status::Unknown), ((3.0, 4.0), Status::Fails)] {
        let mut v = ranged_view(src, &[("p", lo, hi)]);
        let (a, b) = (def_field(&v, "a"), def_field(&v, "b"));
        let c = check(&mut v, ink::lang::Claim::Overlaps(a, b));
        assert_eq!(c.status, want, "{lo}..{hi}: {c:?}");
    }
}

#[test]
fn promises_bound_slopes_in_space_not_over_params() {
    // x - 1000 p promised 1-Lipschitz: true in space, but over p in [0, 1]
    // it moves by 1000. At p = 0 it is clear of the block at x = 5..6 by 5;
    // over the band it is not.
    use ink::claims::Status;
    use ink::lang::Claim;
    let mut v = ranged_view("(param p 0) (def half (lipschitz 1 (- x (* 1000 p)))) (def block (move (box 1 1 1) 5.5 0 0))", &[("p", 0.0, 1.0)]);
    let (half, block) = (def_field(&v, "half"), def_field(&v, "block"));
    assert_eq!(check(&mut v, Claim::Clearance(half, block, 0.5)).status, Status::Fails);
    let both = v.doc.g.max(half, block);
    assert_eq!(check(&mut v, Claim::Empty(both)).status, Status::Fails);
    // Over p in [0, 0.001] its face moves from 0 to 1, against a block from
    // 0.6: at the middle value (face at 0.5) the block is clear by 0.1, and
    // the search bounds the cells near it first order in p. Its rate is the
    // formula's, -1000: the face moves by 0.5 either side. Bounded by the
    // promise instead (1), it would seem to move by 0.0005, and both claims
    // would be proven.
    let mut v = ranged_view("(param p 0) (def half (lipschitz 1 (- x (* 1000 p)))) (def block (move (box 1 1 1) 1.1 0 0))", &[("p", 0.0, 0.001)]);
    let (half, block) = (def_field(&v, "half"), def_field(&v, "block"));
    let c = check(&mut v, Claim::Clearance(half, block, 0.05));
    assert_eq!(c.status, Status::Fails, "{c:?}");
    let both = v.doc.g.max(half, block);
    let c = check(&mut v, Claim::Empty(both));
    assert_eq!(c.status, Status::Fails, "{c:?}");
    let (w, ps) = witness_of(&v, &c);
    assert!(ps[0].1 > 0.0006 && w[0] > 0.6 && w[0] < 1000.0 * ps[0].1, "{w:?} {ps:?}");
}

#[test]
fn a_part_that_jumps_with_a_param_is_found_past_the_jump() {
    // A cube at x = 0 for p < 1 and x = 3 from p = 1 (a select on p),
    // against a cube at x = 3, p in [0, 1.5]. At the middle value (0.75) it
    // is clear, and its position's rate is 0 on either side of the jump: a
    // first-order bound from a point derivative would keep it at 0 for
    // every p and prove the cubes apart. Across the jump its rate is
    // unbounded: the plain bound stands, and the search halves p where a
    // half's bound rises (one half always holds the jump), down to a piece
    // past it. A floor, the same: x = 3 floor(p). The position as a value,
    // the same.
    use ink::claims::Status;
    use ink::lang::Claim;
    use ink::ranged::value_span;
    for at in ["(select (- p 1) 0 3)", "(* 3 (floor p))"] {
        let mut v = ranged_view(&format!("(param p 0) (def x0 {at}) (def a (move (box 1 1 1) x0 0 0)) (def b (move (box 1 1 1) 3 0 0))"), &[("p", 0.0, 1.5)]);
        let (a, b) = (def_field(&v, "a"), def_field(&v, "b"));
        let both = v.doc.g.max(a, b);
        for c in [Claim::Empty(both), Claim::Clearance(a, b, 0.5)] {
            let c = check(&mut v, c);
            assert_eq!(c.status, Status::Fails, "{at}: {c:?}");
            let (_, ps) = witness_of(&v, &c);
            assert!(ps[0].1 >= 1.0, "{at}: {ps:?}");
        }
        let r = value_span(&Tape::new(&v.doc.g, def_field(&v, "x0")), &v.ctx, [0.0; 3], ink::ranged::Comparison::exact(ink::lang::Cmp::Lt, 1.0), &v.proof_options(None, 0.0));
        assert!(r.status == ink::ranged::Status::Disproven && r.witness.as_ref().is_some_and(|(at, x)| at.params[0].1 >= 1.0 && *x == 3.0), "{at}: {r:?}");
    }
}

#[test]
fn one_param_read_twice_is_proven_tight() {
    // A shaft of diameter d in a sleeve whose bore is
    // d + 0.1, d anywhere in 7..9: the radial gap is 0.05 for every d.
    // Plain bounds read d once in each part and blur each surface by the
    // band, so at cells of 0.2 a search over pieces of d leaves both claims
    // unknown (halving d down to 1/256 of its band still blurs each surface
    // by 0.004, and the cells settle nothing near the gap). First order in
    // d, a facing pair's rates cancel (1/2 (-1/2) + 1/2 (1/2) = 0), and the
    // cells along the gap are settled at once wherever d puts it.
    use ink::claims::Status;
    use ink::lang::Claim;
    let src = "(param d 8) (def shaft (cylinder-z (/ d 2) 4)) (def sleeve (subtract (cylinder-z 6 5) (cylinder-z (+ (/ d 2) 0.05) 8)))";
    let mut v = ranged_view(src, &[("d", 7.0, 9.0)]);
    let (shaft, sleeve) = (def_field(&v, "shaft"), def_field(&v, "sleeve"));
    let opts = ink::claims::ClaimOptions { res: Some(0.2), accuracy: Some(0.0), ..Default::default() };
    assert_eq!(ink::claims::check(&mut v, Claim::Clearance(shaft, sleeve, 0.02), None, opts).unwrap().status, Status::Holds);
    let both = v.doc.g.max(shaft, sleeve);
    assert_eq!(ink::claims::check(&mut v, Claim::Empty(both), None, opts).unwrap().status, Status::Holds);
    // Two blocks one param places a gap of 0.01 apart (a block (box a ...)
    // at a / 2, the next at a + 0.01 + 5), a in 9.9..10.1: proven outright
    // with cells five times the gap. Touching (no gap), proven disjoint to
    // an accuracy (not exactly), as at a single value.
    let blocks = |gap: f64| format!("(param a 10) (def pa (move (box a 10 10) (/ a 2) 0 0)) (def pb (move (box 10 10 10) (+ a {gap} 5) 0 0)) (def ab (max pa pb))");
    let mut v = ranged_view(&blocks(0.01), &[("a", 9.9, 10.1)]);
    let ab = def_field(&v, "ab");
    let opts = ink::claims::ClaimOptions { res: Some(0.05), ..opts };
    assert_eq!(ink::claims::check(&mut v, Claim::Empty(ab), None, opts).unwrap().status, Status::Holds);
    let mut v = ranged_view(&blocks(0.0), &[("a", 9.9, 10.1)]);
    let (pa, pb) = (def_field(&v, "pa"), def_field(&v, "pb"));
    let c = ink::claims::check(&mut v, Claim::Disjoint(pa, pb), None, ink::claims::ClaimOptions { accuracy: Some(0.1), ..opts }).unwrap();
    let ink::claims::Evidence::Empty(p) = &c.evidence else { panic!("{c:?}") };
    assert!(p.verdict == ink::query::Verdict::Proven && !p.exact && p.tol == 0.1, "{c:?}");
}

#[test]
fn time_and_params_compose() {
    // A unit cube sliding 10 t + p along x, t in [0, 1], p in [-0.5, 0.5]:
    // its face reaches 11 at t = 1, p = 0.5. A wall from 10.8 is hit (at
    // end, t near 1, and p > 0.3); one from 11.2 never is, proven with t
    // and p ranged in one search.
    use ink::claims::Status;
    let src = "(param p 0) (def cube (move (box 1 1 1) (+ (* 10 t) p) 0 0))
(def near (max cube (move (box 2 4 4) 11.8 0 0)))
(def far (max cube (move (box 2 4 4) 12.2 0 0)))";
    let mut v = ranged_view(src, &[("p", -0.5, 0.5)]);
    v.set_motion(Some((0.0, 1.0))).unwrap();
    let (near, far) = (def_field(&v, "near"), def_field(&v, "far"));
    let c = check(&mut v, ink::lang::Claim::Empty(near));
    assert_eq!(c.status, Status::Fails, "{c:?}");
    let ink::claims::Evidence::Empty(p) = &c.evidence else { panic!() };
    let ink::query::Verdict::Disproven { witness: w, params, .. } = &p.verdict else { panic!("{p:?}") };
    let t = c.motion.as_ref().and_then(|m| m.motion.t).unwrap();
    assert!(10.0 * t + params[0].1 + 0.5 > 10.8 && w[0] > 10.8 && (0.0..=1.0).contains(&t), "{t} {params:?} {w:?}");
    let c = check(&mut v, ink::lang::Claim::Empty(far));
    assert!(c.status == Status::Holds && c.accuracy.is_some(), "{c:?}");
}

/// A model's view with its declared tolerances applied (as `ink assert`).
fn toleranced(src: &str) -> Result<ink::view::ModelView, ink::InkError> {
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(d, ctx);
    let spans = v.declared_spans()?;
    v.set_param_spans(&spans)?;
    Ok(v)
}

/// A root of equations over sampled fields is certified for the samples
/// only, and an uncertified root only for the value found: what reads
/// either, directly or through another solve, is measured, naming them.
/// u is solved on a baked sphere (0.9793; 0.9890 on the sphere), v = u + p;
/// w stops at its guess 1 (the root is 2).
#[test]
fn answers_on_sampled_or_uncertified_roots_are_measured() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let src = "(def sampled (grid-bake (sphere 1) (res 0.25)))
(param p 0 (tol 0 0.001))
(solve first (unknowns (u 1)) (== (warp sampled u 0.07 0.13) p))
(solve second (unknowns (v 1)) (== v (+ u p)))
(solve early (unknowns (w 1)) (== (* 1e-13 (- w 2)) 0))
(assert false-radius (value u < 0.985))
(assert compare-source (value (warp (sphere 1) u 0.07 0.13) > -0.001))
(assert downstream (value v < 0.986))
(assert early-root (value w < 1.5))";
    let mut v = toleranced(src).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let w = v.doc.param_index("w").unwrap() as u32;
    assert!(r.asserts.iter().all(|a| a.check.asked.is_none()), "estimates: {r:?}");
    let got: Vec<(Status, usize, Vec<u32>)> = r.asserts.iter().map(|a| (a.check.status, a.check.inexact.samples.len(), a.check.inexact.solved.clone())).collect();
    use Status::*;
    assert_eq!(got, [(Holds, 1, vec![]), (Fails, 1, vec![]), (Holds, 1, vec![]), (Holds, 0, vec![w])], "{r:?}");
}

#[test]
fn tolerance_bands_hold_their_real_ends() {
    // Deviations from the value, the ends exact where the sums are, else
    // the nearest float outside: 8 (tol 0 0.022) is 8.000..8.022, its top
    // end the least float at least 8 + 0.022 (x - 8 is exact here).
    let d =
        ink::load_str("(param a 8 (tol 0 0.022)) (param b 8 (tol -0.5 0.25)) (param c 1 (tol 0 8.673617379884035e-19)) (param e 8 (tol -0.022 0))", "t", Path::new(ROOT)).unwrap();
    let band = |n: &str| d.params[d.param_index(n).unwrap()].band(d.params[d.param_index(n).unwrap()].default).unwrap();
    let (lo, hi) = band("a");
    assert!(lo == 8.0 && hi - 8.0 >= 0.022 && hi.next_down() - 8.0 < 0.022, "{hi}");
    let (lo, hi) = band("e");
    assert!(hi == 8.0 && 8.0 - lo >= 0.022 && 8.0 - lo.next_up() < 0.022, "{lo}");
    assert_eq!(band("b"), (7.5, 8.25));
    // 1 + 2^-60 is between 1 and the next float.
    assert_eq!(band("c"), (1.0, 1f64.next_up()));
    for bad in ["(param a 8 (tol 1 0))", "(param a 8 (tol 0))", "(param a 8 (tolerance 0 1))"] {
        assert!(ink::load_str(bad, "t", Path::new(ROOT)).is_err(), "{bad}");
    }
}

#[test]
fn a_param_whose_value_is_used_at_load_time_is_not_ranged() {
    // An if on n makes it one number at load time: ranging it is refused,
    // naming where. A claim's bound is kept as an expression (not used at
    // load time), and a field of a param stays symbolic: those range.
    let src = "(param n 3 (tol 0 1)) (param w 2 (tol 0 0.5)) (param b 1 (tol 0 1))
(def plate (if (> n 2) (box 10 10 1) (box 5 5 1)))
(def bar (box w 1 1))
(assert small (value w < b))";
    let e = toleranced(src).err().expect("refused");
    assert!(e.message.contains("'n'") && e.message.contains("test:3 (def plate"), "{e:?}");
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
    assert!(d.materialized(d.param_index("n").unwrap()).is_some());
    assert!(d.materialized(d.param_index("w").unwrap()).is_none() && d.materialized(d.param_index("b").unwrap()).is_none());
    // A solve's unknown cannot shape the model at load time: the number
    // would hide that it was solved.
    let e = ink::load_str("(param p 1 (tol 0 1)) (solve s (unknowns (a 1)) (== a p)) (def n (if (> a 0.5) 1 2))", "test", Path::new(ROOT)).err().expect("refused");
    assert!(e.message.contains("reads a solved value, which cannot shape the model at load time"), "{e:?}");
}

/// A solve over ranged params: its unknowns are the one root
/// of its equations, certified over the ranges (Krawczyk) and narrowed on
/// each piece the search makes, so claims on them are proven, and a
/// witness names the params where they fail (a solve at the values would
/// put a = 1 and prove a < 1.5).
#[test]
fn a_solve_over_ranged_params_is_proven_over_them() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let src = "(param p 1 (tol 0 1)) (solve s (unknowns (a 1)) (== a p))
(assert above (value a > 0.999)) (assert below (value a < 2.001)) (assert under (value a < 1.5))";
    let mut v = toleranced(src).unwrap();
    let a = v.doc.param_index("a").unwrap() as u32;
    let band = v.ctx.p_spans.get(a).expect("a is derived over p");
    assert!(band.0 <= 1.0 && 2.0 <= band.1 && band.1 - band.0 < 1.0 + 1e-9, "{band:?}");
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let st: Vec<Status> = r.asserts.iter().map(|a| a.check.status).collect();
    assert_eq!(st, [Status::Holds, Status::Holds, Status::Fails], "{r:?}");
    // Each claim covers p, the param a is derived from.
    assert!(r.asserts.iter().all(|a| a.check.ranged.iter().map(|r| v.doc.params[r.0 as usize].name.as_str()).eq(["p"])), "{r:?}");
    let ink::claims::Evidence::ValueSpan(result) = &r.asserts[2].check.evidence else { panic!("{:?}", r.asserts[2].check) };
    let (at, value) = result.witness.as_ref().expect("a witness");
    let p = at.params.iter().find(|q| v.doc.params[q.0 as usize].name == "p").expect("p in the witness").1;
    assert!((1.5..=2.0).contains(&p) && *value >= 1.5, "{at:?} {value}");
    // At the values the solve is a numerical method's answer: unchanged.
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(d, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    assert!(r.asserts.iter().all(|a| a.check.status == Status::Holds && a.check.accuracy.is_none()), "{r:?}");
}

/// The enclosure of a nonlinear solve holds every solution over the ranges,
/// and at a point narrows to that point's root: a circle of radius r
/// (r in [4.9, 5.1]) meets the line y = x - 1, solved for (px, py).
#[test]
fn a_solved_enclosure_holds_the_root_at_every_value_and_narrows_at_points() {
    let src = "(param r 5 (tol -0.1 0.1))
(solve meet (unknowns (px 3) (py 2)) (== (+ (square px) (square py)) (square r)) (== py (- px 1)))";
    let v = toleranced(src).unwrap();
    let [k, px, py] = ["r", "px", "py"].map(|n| v.doc.param_index(n).unwrap());
    let (bx, by) = (v.ctx.p_spans.get(px as u32).unwrap(), v.ctx.p_spans.get(py as u32).unwrap());
    for i in 0..=20 {
        let r = 4.9 + 0.01 * i as f64;
        // px = (1 + sqrt(2 r^2 - 1)) / 2, from (x - 1)^2 + x^2 = r^2.
        let x = (1.0 + (2.0 * r * r - 1.0).sqrt()) / 2.0;
        // The solver at the values: a numerical answer, near the root.
        let lx = v.doc.ctx(&[("r".into(), r.into())], 0.0).unwrap().params[px];
        assert!((lx - x).abs() < 1e-9 && bx.0 <= x && x <= bx.1 && by.0 <= x - 1.0 && x - 1.0 <= by.1, "{r}: {x} {bx:?} {by:?}");
        // At the point: an enclosure a few ulps wide around the root (x in
        // floats is within an ulp or two of it), its value inside.
        let mut c = v.ctx.clone();
        c.params[k] = r;
        let c = c.point();
        let (sx, sy) = (c.p_spans.get(px as u32).unwrap(), c.p_spans.get(py as u32).unwrap());
        assert!(sx.0 - 1e-14 <= x && x <= sx.1 + 1e-14 && sy.0 - 1e-14 <= x - 1.0 && x - 1.0 <= sy.1 + 1e-14, "{r}: {sx:?} {x}");
        assert!(sx.1 - sx.0 < 1e-12 && sy.1 - sy.0 < 1e-12, "{r}: {sx:?} {sy:?}");
        assert!(sx.0 <= c.params[px] && c.params[px] <= sx.1 && sy.0 <= c.params[py] && c.params[py] <= sy.1);
    }
    // The band's enclosure is about as wide as the solutions' spread
    // (dpx/dr = r / sqrt(2 r^2 - 1), 0.72 at r = 5: 0.14 over the band).
    assert!(bx.1 - bx.0 < 0.16, "{bx:?}");
    // A piece of the range narrows it.
    let half = v.ctx.with_param(k as u32, (4.9, 5.0));
    let hx = half.p_spans.get(px as u32).unwrap();
    assert!(hx.1 - hx.0 < 0.6 * (bx.1 - bx.0) && hx.0 >= bx.0 && hx.1 <= bx.1, "{hx:?} {bx:?}");
}

/// A four-bar's loop closure over four toleranced inputs: at random points
/// of the bands, and over random pieces of them holding those points, the
/// enclosure holds the root the solver finds there (to the solver's
/// accuracy), narrowing to it at points.
#[test]
fn a_four_bar_solved_over_tolerances_holds_its_root_everywhere() {
    let src = "(param phi 0.7 (tol -0.01 0.01)) (param b 60 (tol -0.05 0.05)) (param c 50 (tol -0.05 0.05)) (param d 55 (tol -0.05 0.05))
(solve loop (unknowns (beta 0.8) (gamma 1.6))
  (== (+ (* 15 (cos phi)) (* b (cos beta))) (+ d (* c (cos gamma))))
  (== (+ (* 15 (sin phi)) (* b (sin beta))) (* c (sin gamma))))";
    let v = toleranced(src).unwrap();
    let ks = ["phi", "b", "c", "d"].map(|n| v.doc.param_index(n).unwrap());
    let us = ["beta", "gamma"].map(|n| v.doc.param_index(n).unwrap());
    let bands: Vec<(f64, f64)> = ks.iter().map(|&k| v.ctx.p_spans.get(k as u32).unwrap()).collect();
    let mut seed = 0x9e3779b97f4a7c15u64;
    let mut rnd = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let holds = |s: (f64, f64), x: f64| s.0 - 1e-12 <= x && x <= s.1 + 1e-12;
    for _ in 0..40 {
        let at: Vec<f64> = bands.iter().map(|b| b.0 + (b.1 - b.0) * rnd()).collect();
        let overrides: Vec<(String, ink::lang::Given)> = ["phi", "b", "c", "d"].iter().zip(&at).map(|(n, &x)| (n.to_string(), x.into())).collect();
        let solved = v.doc.ctx(&overrides, 0.0).unwrap();
        let root = us.map(|u| solved.params[u]);
        // The band's enclosure, a piece of the bands around the point, and
        // the point.
        let mut piece = v.ctx.clone();
        for ((&k, b), &x) in ks.iter().zip(&bands).zip(&at) {
            let (lo, hi) = (x - (x - b.0) * rnd(), x + (b.1 - x) * rnd());
            piece = piece.with_param(k as u32, (lo, hi));
        }
        let mut point = piece.clone();
        for (&k, &x) in ks.iter().zip(&at) {
            point.params[k] = x;
        }
        let point = point.point();
        for (i, &u) in us.iter().enumerate() {
            for c in [&v.ctx, &piece, &point] {
                let s = c.p_spans.get(u as u32).unwrap();
                assert!(holds(s, root[i]), "{at:?}: {} = {} outside {s:?}", ["beta", "gamma"][i], root[i]);
            }
            let s = point.p_spans.get(u as u32).unwrap();
            assert!(s.1 - s.0 < 1e-12 && holds(s, point.params[u]), "{s:?}");
        }
    }
}

/// What cannot be certified is refused, never solved at one value and
/// labelled as covering the range: two roots that meet in the range, a
/// singular Jacobian at the one root, more equations than unknowns (least
/// squares).
#[test]
fn a_solve_that_cannot_be_certified_over_ranges_is_refused() {
    let refused = |src: &str, why: &str| {
        let e = toleranced(src).err().unwrap_or_else(|| panic!("{src}: not refused"));
        assert!(e.message.contains("solve 's' depends on ranged 'p'") && e.message.contains(why), "{src}: {e:?}");
    };
    // a^2 = p: the roots +-sqrt(p) meet at p = 0, inside the range.
    refused("(param p 0.0001 (tol -0.001 0.001)) (solve s (unknowns (a 0.01)) (== (square a) p))", "no single root");
    // (a - p)^3 = 0: one root, a = p, where the Jacobian is 0.
    refused("(param p 1 (tol 0 0.1)) (solve s (unknowns (a 1)) (== (* (- a p) (square (- a p))) 0))", "");
    refused("(param p 1 (tol 0 0.1)) (solve s (unknowns (a 1)) (== a p) (== (* 2 a) (* 2 p)))", "2 equations for 1 unknowns");
}

/// A solve may read what an earlier one sets: its enclosure is over the
/// earlier one's, which follows the ranged params.
#[test]
fn chained_solves_derive_over_ranges_in_order() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let src = "(param p 1 (tol 0 0.1))
(solve s1 (unknowns (a 1)) (== a p))
(solve s2 (unknowns (b 2)) (== (* b b) (* 4 (square a))))
(assert b-range (value b <= 2.2001)) (assert b-low (value b < 2.1))";
    let mut v = toleranced(src).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let st: Vec<Status> = r.asserts.iter().map(|a| a.check.status).collect();
    assert_eq!(st, [Status::Holds, Status::Fails], "{r:?}");
    assert!(r.asserts[1].check.ranged.iter().map(|r| v.doc.params[r.0 as usize].name.as_str()).eq(["p"]), "{r:?}");
}

/// A root whose equations have a kink (min u 2u = p: u = p for p >= 0, p /
/// 2 below) is narrowed on each piece of the range, though the enclosure
/// over all of it straddles the kink (where K over it cannot narrow it: a
/// box of the piece's own, inside it, is certified instead). True ends:
/// -0.05 and 0.1.
#[test]
fn a_root_at_a_kink_is_narrowed_on_each_piece() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let src = "(param p 0 (tol -0.1 0.1)) (solve s (unknowns (u 0)) (== (min u (* 2 u)) p))
(assert high (value u < 0.101)) (assert low (value u > -0.051))
(assert over-high (value u < 0.099)) (assert over-low (value u > -0.049))";
    let mut v = toleranced(src).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let st: Vec<Status> = r.asserts.iter().map(|a| a.check.status).collect();
    assert_eq!(st, [Status::Holds, Status::Holds, Status::Fails, Status::Fails], "{r:?}");
}

/// A certificate is about the equations as they read their inputs: a root
/// certified over p but resting on an uncertified solve's answer (a double
/// root, a numerical method's approximation of 2) is numerical. Here a = 2
/// exactly and b = 2 + p >= 2, but the solver's a is 1.99999905, and the
/// root certified for it would prove b < 2.
#[test]
fn a_root_resting_on_an_uncertified_solve_is_numerical() {
    use ink::claims::{ClaimOptions, Evidence, Status, run_asserts};
    let src = "(param p 0 (tol 0 0.00000001))
(solve first (unknowns (a 1)) (== (square (- a 2)) 0))
(solve second (unknowns (b 2)) (== b (+ a p)))
(solve third (unknowns (c 0)) (== c (* 2 b)))
(assert false-proof (value b < 2)) (assert downstream (value c < 4))";
    let mut v = toleranced(src).unwrap();
    // b and c are derived over p (they follow it), but not certified.
    let [b, c] = ["b", "c"].map(|n| v.doc.param_index(n).unwrap() as u32);
    assert!(v.ctx.p_spans.is_derived(b) && v.ctx.p_spans.is_derived(c));
    assert!(v.certified().is_empty(), "{:?}", v.certified());
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    assert!(r.asserts.iter().all(|a| a.check.status == Status::Holds && a.check.accuracy.is_none()), "{r:?}");
    // And so does every part of an answer: an integral's, a mean's.
    let src = "(param p 0 (tol 0 0.00000001))
(solve first (unknowns (a 1)) (== (square (- a 2)) 0))
(solve second (unknowns (b 2)) (== b (+ a p)))
(def cube (box 1 1 1)) (assert false-integral (integral cube (- b 2) < 0)) (assert false-mean (mean cube b < 2))";
    let mut v = toleranced(src).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions { res: Some(0.1), ..Default::default() }).unwrap();
    let parts: Vec<(Status, Status)> = r
        .asserts
        .iter()
        .map(|a| match &a.check.evidence {
            Evidence::Integral { integrals, .. } => (a.check.status, integrals[0].1),
            Evidence::Mean { means, .. } => (a.check.status, means[0].1),
            e => panic!("{e:?}"),
        })
        .collect();
    assert_eq!(parts, [(Status::Holds, Status::Holds); 2], "{r:?}");
    assert!(r.asserts.iter().all(|a| a.check.asked.is_none()), "estimates, every part: {r:?}");
}

/// Solves are applied in the order of what they read, not of the file: a
/// solve reading a param a later one sets comes after it, at single values
/// and over the ranges alike. Here a = 2 p and b = a + p = 3 p, p in [0, 1].
/// Two solves that read each other's unknowns are one system: refused,
/// naming both.
#[test]
fn a_solve_reading_what_a_later_one_sets_comes_after_it() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let src = "(param p 0 (tol 0 1)) (param a 0)
(solve early (unknowns (b 0)) (== b (+ a p)))
(solve late (unknowns a) (== a (* 2 p)))
(assert false-high (value b <= 1.01)) (assert equations (value (abs (- b (+ a p))) < 0.01))";
    let mut v = toleranced(src).unwrap();
    let [a, b] = ["a", "b"].map(|n| v.doc.param_index(n).unwrap());
    let band = v.ctx.p_spans.get(b as u32).unwrap();
    assert!(band.0 <= 0.0 && 3.0 <= band.1 && band.1 < 3.01, "{band:?}");
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let st: Vec<Status> = r.asserts.iter().map(|a| a.check.status).collect();
    assert_eq!(st, [Status::Fails, Status::Holds], "{r:?}");
    // At a single value: late first, then early.
    let ctx = v.doc.ctx(&[("p".into(), 1.0.into())], 0.0).unwrap();
    assert!((ctx.params[a] - 2.0).abs() < 1e-9 && (ctx.params[b] - 3.0).abs() < 1e-9, "{:?}", ctx.params);
    let e = ink::load_str("(param p 0) (param a 0) (solve one (unknowns (b 0)) (== b (+ a p))) (solve two (unknowns a) (== a (- 1 b)))", "test", Path::new(ROOT))
        .err()
        .expect("a cycle");
    assert!(e.message.contains("solve 'one' reads 'a', which solve 'two' sets") && e.message.contains("solve 'two' reads 'b', which solve 'one' sets"), "{e:?}");
}

/// Solving at load time waits for the solves after: each is
/// solved with those read so far, and one that fails is the error only if
/// the document as read does not solve, at its own form. A later solve that
/// would change a value already used at load time is refused, naming the
/// use: that use is one number, made without it. Neither order permits a
/// load-time use; both errors direct the author to fields.
#[test]
fn load_time_values_are_at_every_solve_that_changes_them() {
    let load = |src: &str| ink::load_str(src, "test", Path::new(ROOT));
    let solvable = "(param a -1) (solve early (unknowns (b 1)) (== (square b) a))";
    let d = load(&format!("{solvable} (solve late (unknowns a) (== a 1))")).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let b = d.param_index("b").unwrap();
    assert!((ctx.params[b] - 1.0).abs() < 1e-9 && (d.params[b].default - 1.0).abs() < 1e-9, "{:?}", ctx.params);
    let e = load(&format!("{solvable}\n(def w 1)")).err().unwrap();
    assert!(e.message.contains("solve 'early' did not converge") && e.line == 1, "{e:?}");

    // A plain param used at load time cannot be set by a later solve, and a
    // solved one cannot be used at load time at all.
    let e = load("(param px 0)\n(def n (if (> px 0.5) 1 2))\n(solve s (unknowns px) (== px 1))").err().unwrap();
    assert!(e.message.contains("solve 's' sets 'px', whose value is used at load time at test:2 (def n ...)"), "{e:?}");
    assert!(e.hint.as_deref().unwrap_or("").contains("use it in fields (placements, shapes, claims)"), "{e:?}");
    let e = load("(param px 0)\n(solve s (unknowns px) (== px 1))\n(def n (if (> px 0.5) 1 2))").err().unwrap();
    assert!(e.message.contains("reads a solved value") && e.line == 3, "{e:?}");
    assert!(e.hint.as_deref().unwrap_or("").contains("use it in fields (placements, shapes, claims)"), "{e:?}");
}

/// A claim whose bound is all that reads a root of equations over samples
/// is an estimate too, and one holding in a window as an estimate says so
/// there as well (its tolerance is every status's: Check::estimate).
#[test]
fn a_claim_rests_on_its_bound_and_is_an_estimate_in_every_part() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let src = "(def sampled (grid-bake (sphere 1) (res 0.25)))
(solve s (unknowns (u 1)) (== (warp sampled u 0.07 0.13) 0))
(assert bound-only (value 0.5 < u))";
    let mut v = toleranced(src).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let c = &r.asserts[0].check;
    assert!(c.status == Status::Holds && c.accuracy.is_none() && c.inexact.samples.len() == 1, "{c:?}");
    let src = "(solve s (unknowns (u 0.9)) (== (* 1e-13 (+ u 0.9)) 0))
(def waves (+ (sin (* 0.001 x)) u))
(assert waves-empty (empty waves))";
    let mut v = toleranced(src).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions { res: Some(0.05), ..Default::default() }).unwrap();
    let c = &r.asserts[0].check;
    assert!(c.status == Status::Unknown && c.in_window == Some(Status::Holds) && c.accuracy.is_none() && c.asked.is_none(), "{c:?}");
}

/// A part placed by a solve over a ranged param is found wherever the
/// solve puts it: a ball of radius 1 at x = a, a = 2 p (p in [0, 1]),
/// reaches x = 2 p + 1 <= 3. A wall from x = 3.5 is clear of it; one from
/// x = 2.5 is clear at the value (p = 0) but not for p > 0.75: disproven,
/// the witness inside both at its params.
#[test]
fn a_part_placed_by_a_solve_is_checked_over_the_ranges() {
    use ink::claims::Status;
    let src = "(param p 0 (tol 0 1)) (solve s (unknowns (a 0)) (== a (* 2 p)))
(def ball (move (sphere 1) a 0 0)) (def wall (move (box 1 4 4) 4 0 0)) (def near (move (box 1 4 4) 3 0 0))
(def both (max ball wall)) (def both-near (max ball near))";
    let mut v = toleranced(src).unwrap();
    let (both, near) = (def_field(&v, "both"), def_field(&v, "both-near"));
    assert_eq!(check(&mut v, ink::lang::Claim::Empty(both)).status, Status::Holds);
    let c = check(&mut v, ink::lang::Claim::Empty(near));
    assert_eq!(c.status, Status::Fails, "{c:?}");
    let (w, ps) = witness_of(&v, &c);
    let p = ps.iter().find(|q| q.0 == "p").expect("p").1;
    assert!(p > 0.75 && w[0] > 2.5 && w[0] < 2.0 * p + 1.0, "{w:?} {ps:?}");
    // The witness names the independent params only; a follows from p.
    assert!(ps.iter().all(|q| q.0 != "a"), "{ps:?}");
}

#[test]
fn a_bound_that_varies_is_part_of_the_value_and_refused_elsewhere() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    // 1 < p for p in [0, 2] fails for p <= 1; the bound sampled at the
    // value (2, or any one) would say it holds.
    let mut v = toleranced("(param p 2 (tol -2 0)) (assert above (value 1 < p)) (assert below (value p <= 2))").unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    assert_eq!(r.asserts.iter().map(|a| a.check.status).collect::<Vec<_>>(), vec![Status::Fails, Status::Holds], "{r:?}");
    assert!(r.asserts[0].check.note.is_some() && r.asserts[0].check.ranged.len() == 1);
    // A clearance distance is one number.
    let mut v = toleranced("(param d 1 (tol 0 1)) (assert gap (clearance (sphere 1) (move (sphere 1) 5 0 0) d))").unwrap();
    let e = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap_err();
    assert!(e.message.contains("assert 'gap'") && e.message.contains("d"), "{e:?}");
}

/// `within`'s corners may read params: the box stays a
/// field of them, so a ranged one is proven over its band, not at its value.
#[test]
fn within_corners_read_params_over_their_bands() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    // A ball of radius 1.5. a in [2, 2.5]: inside the box at every a. w in
    // [1.3, 1.6]: inside at the value 1.6, out for w < 1.5.
    let src = "(param a 2 (tol 0 0.5)) (param w 1.6 (tol -0.3 0))
(body s (shape (sphere 1.5)))
(assert roomy (within (shape-of s) ((- 0 a) -2 -2) (a 2 2)))
(assert tight (within (shape-of s) ((- 0 a) -2 -2) (w 2 2)))";
    let mut v = toleranced(src).unwrap();
    for k in ["a", "w"] {
        assert!(v.doc.materialized(v.doc.param_index(k).unwrap()).is_none(), "{k}");
    }
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let (roomy, tight) = (&r.asserts[0].check, &r.asserts[1].check);
    assert!(roomy.status == Status::Holds && roomy.ranged.iter().any(|(k, _)| v.doc.params[*k as usize].name == "a"), "{roomy:?}");
    assert_eq!(tight.status, Status::Fails, "{tight:?}");
    match &tight.evidence {
        ink::claims::Evidence::Empty(p) => match &p.verdict {
            ink::query::Verdict::Disproven { witness, params, .. } => {
                let k = v.doc.param_index("w").unwrap() as u32;
                let w = params.iter().find(|q| q.0 == k).map(|q| q.1).expect("w in the witness");
                assert!(w < witness[0] && witness[0] < 1.5, "{w} {witness:?}");
            }
            o => panic!("{o:?}"),
        },
        e => panic!("{e:?}"),
    }
    // At the values only, both hold.
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(d, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    assert!(r.all_hold_guaranteed(), "{r:?}");
    // A corner is a point of numbers, not of space.
    let e = ink::load_str("(import \"std\") (assert c (within (sphere 1) (-2 -2 -2) (x 2 2)))", "t", Path::new(ROOT)).err().expect("refused");
    assert!(e.message.contains("corner"), "{e:?}");
}

/// An assert states the cell size it needs, and is checked
/// at the finer of that and the run's: a hexagonal prism 0.0025 inside a
/// bore is proven (to the tolerance) at either.
#[test]
fn an_assert_is_checked_at_its_own_res_or_the_runs_if_finer() {
    use ink::claims::{ClaimOptions, ResBy, Status, run_asserts};
    let src = "(import \"std\")
(def hex (extrude-z (polygon (regular-polygon-points 6 4.9975)) -2 2))
(assert tight (subset hex (cylinder-z 5 6)) (res 0.01))
(assert plain (subset hex (cylinder-z 5 6)))";
    let doc = ink::load_str(src, "claims", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let mut run = |res: Option<f64>| {
        run_asserts(&mut v, None, false, ClaimOptions { res, ..ClaimOptions::default() }).unwrap().asserts.into_iter().map(|a| (a.check.status, a.res)).collect::<Vec<_>>()
    };
    assert_eq!(run(Some(0.02)), vec![(Status::Holds, Some((0.01, ResBy::Claim))), (Status::Holds, Some((0.02, ResBy::Run)))]);
    assert_eq!(run(Some(0.005)), vec![(Status::Holds, Some((0.005, ResBy::Run))), (Status::Holds, Some((0.005, ResBy::Run)))]);
    assert_eq!(run(None)[0], (Status::Holds, Some((0.01, ResBy::Claim))));
    assert_eq!(run(None)[1].1, None);
    // A res is a constant cell size above 0 and finite, given once; it is
    // an option of claims checked on cells.
    for (bad, says) in [
        ("(res 0)", "above 0"),
        ("(res -0.1)", "above 0"),
        ("(res (exp 1000))", "above 0"),
        ("(res x)", "constant"),
        ("(res 0.1) (res 0.2)", "more than one"),
        ("(scope 1)", "unknown assert option"),
        ("(accuracy 0)", "above 0"),
        ("(accuracy p)", "constant"),
        ("(accuracy 0.1) (accuracy 0.2)", "more than one"),
        ("(tol 0.1)", "A claim states its accuracy"),
    ] {
        let e = ink::load_str(&format!("(import \"std\") (param p 1) (assert c (empty (sphere 1)) {bad})"), "t", Path::new(ROOT)).err().expect(bad);
        assert!(e.message.contains(says), "{bad}: {e:?}");
    }
    let e = ink::load_str("(param p 1) (assert c (value p < 2) (res p))", "t", Path::new(ROOT)).err().expect("a param res");
    assert!(e.message.contains("constant"), "{e:?}");
    // An accuracy is in the claim's units: a shape claim's a length, a
    // value's its own; a range or a point is decided exactly and has none.
    for claim in ["(range x -1 1 (sphere 1))", "(subset (probe 0 0 0) (sphere 1))"] {
        let e = ink::load_str(&format!("(import \"std\") (assert c {claim} (accuracy 0.1))"), "t", Path::new(ROOT)).err().expect(claim);
        assert!(e.message.contains("decided exactly"), "{claim}: {e:?}");
    }
    let d = ink::load_str("(assert c (value 1 < 2) (accuracy 0.1))", "t", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!(d.asserts[0].accuracy, Some(0.1));
}

/// Asserts in `src` (with std) run on a model view: each one's status and
/// evidence, by name.
fn run_claims(src: &str) -> Vec<(String, ink::claims::Check)> {
    use ink::claims::{ClaimOptions, run_asserts};
    let doc = ink::load_str(&format!("(import \"std\")\n{src}"), "claims", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    r.asserts.iter().map(|a| (a.name.clone(), a.check.clone())).collect()
}

/// An unknown clearance examined only in the window holds to its
/// attained accuracy there, not for the shapes: two tiny cubes turned 45
/// degrees about z, corner to corner 4.25 apart (each field reads 3 at the
/// other), and two balls in the default window, each shape also the
/// half-space x > 2000 beyond it, where they overlap (distance 0). The
/// clearance of 3.5 is unknown, holding to some accuracy in the window
/// only, and `ink distance`'s lower bound stays 0, not the window's
/// attained bound (2.99).
#[test]
fn an_attained_bound_in_the_window_says_nothing_past_it() {
    use ink::claims::{ClaimOptions, Status};
    let src =
        "(fn tiny (cx cy) (max (- (abs (* (+ (- x cx) (- y cy)) .7071067811865476)) .005) (- (abs (* (- (- x cx) (- y cy)) .7071067811865476)) .005) (- (abs (- z .345)) .005)))
(fn ball (cx) (- (sqrt (+ (square (- x cx)) (square y) (square z))) 1))
(def a (min (tiny .123 .234) (ball 50) (- 2000 x)))
(def b (min (tiny 4.387 .234) (ball 59) (- 2000 x)))";
    let doc = ink::load_str(src, "window", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let field = |n: &str| match doc.def(n).unwrap() {
        Value::Field(f, _) => *f,
        _ => panic!("{n}"),
    };
    let (a, b) = (field("a"), field("b"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    assert!(!v.whole(a) && !v.whole(b));
    let c = ink::claims::check(&mut v, ink::lang::Claim::Clearance(a, b, 3.5), None, ClaimOptions::default()).unwrap();
    assert!(c.status == Status::Unknown && c.unbounded && matches!(c.attained(), Some((_, true))), "{c:?}");
    let d = ink::claims::distance(&mut v, a, b, ClaimOptions::default()).unwrap();
    assert_eq!(d.lower, 0.0, "they overlap past the window: {} proven", d.lower);
}

#[test]
fn clearance_between_unbounded_slabs_is_proven_over_all_of_space() {
    // Two slabs 13.98 apart against 13.9. max(a, b - 13.9) is never
    // negative, which interval bounds over pieces reaching infinity cannot
    // see, but its first-order bound can (the faces' slopes cancel exactly,
    // turned or not), so the claim covers all of space, not only the window.
    use ink::claims::{Evidence, Status};
    let c = run_claims(
        "(def ra (rotate-x (rotate-z (slab-z 0 1) 25) 30))\n(def rb (rotate-x (rotate-z (slab-z 1.5 3) 25) 30))\n\
         (assert belts (clearance (slab-z 43.99 50.01) (slab-z 63.99 70.01) 13.9))\n\
         (assert turned (clearance ra rb 0.45))\n(assert turned-apart (disjoint ra rb))\n\
         (assert crossing (disjoint (slab-z 0 1) (max (neg x) (- x 1))))\n",
    );
    let at = |n: &str| &c.iter().find(|x| x.0 == n).unwrap().1;
    assert!(at("belts").status == Status::Holds && matches!(at("belts").evidence, Evidence::ClearanceSearch { .. }) && !at("belts").unbounded, "{:?}", at("belts"));
    assert!(at("turned").status == Status::Holds && !at("turned").unbounded, "{:?}", at("turned"));
    assert!(at("turned-apart").status == Status::Holds && !at("turned-apart").unbounded, "{:?}", at("turned-apart"));
    assert_eq!(at("crossing").status, Status::Fails);
}

#[test]
fn clearance_from_a_shape_beyond_the_window_searches_around_its_bounded_partner() {
    // The half-space x > 2000 is not proven bounded and is empty in the
    // window (the 1 m cube), which says nothing about the rest of it: a
    // ball, proven bounded, gives the search its region. A ball at 2000
    // meets it; one at 1900 clears it by 99, over all of space.
    use ink::claims::Status;
    let c = run_claims("(def a (- 2000 x))\n(assert meets (clearance a (move (sphere 1) 2000 0 0) 0.5))\n(assert clear (clearance a (move (sphere 1) 1900 0 0) 10))\n");
    assert_eq!(c[0].1.status, Status::Fails, "{:?}", c[0].1);
    assert!(c[1].1.status == Status::Holds && !c[1].1.unbounded, "{:?}", c[1].1);
}

#[test]
fn clearance_at_a_coarse_res_given_still_refines() {
    // Two balls 0.01 apart, cleared by 0.009. A res given is a size the
    // search reaches, its first cells where coarser than its own: at 1, as
    // at 0.01, the cells halve down to the claim's precision, and it holds
    // (at 1, not unknown for want of a cell finer than the res).
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let at = |res: f64| {
        let src = "(import \"std\")\n(assert clear (clearance (sphere 1) (move (sphere 1) 2.01 0 0) 0.009) (accuracy 0.0005))\n";
        let doc = ink::load_str(src, "claims", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
        let ctx = doc.ctx(&[], 0.0).unwrap();
        let mut v = ink::view::ModelView::new(doc, ctx);
        run_asserts(&mut v, None, false, ClaimOptions { res: Some(res), ..ClaimOptions::default() }).unwrap().asserts[0].check.clone()
    };
    for res in [0.01, 1.0] {
        let c = at(res);
        assert!(c.status == Status::Holds, "{res}: {c:?}");
    }
}

#[test]
fn clearance_slope_is_bounded_where_the_shapes_meet() {
    // A collar 0.25 inside a cup
    // turned from a polygon. The collar, cut from a socket whose taper is
    // turned from a polygon too, has an unbounded slope bound over its box
    // grown by the clearance (the taper and collar both lie within a few mm
    // of its corners, so their min is undecided and the polygon's sign
    // enters it). Only segments from one shape to the other no longer than
    // the clearance need the slope, so cells with neither shape within the
    // clearance are left out, and those near both give 1: the offset form,
    // as written by hand (disjoint hub (- collar 0.2)).
    use ink::claims::{Evidence, Status};
    let src = "(def hub (subtract (move (cylinder-z 22 8) 0 0 4) (revolve-z (polygon (vec (vec -1 1) (vec 19.5 1) (vec 19.5 9) (vec -1 9))))))\n\
               (def taper (revolve-z (polygon (vec (vec 0 12.5) (vec 19.25 12.5) (vec 17 15) (vec 0 15)))))\n\
               (def socket (min (move (cylinder-z 19.25 12) 0 0 7) taper))\n\
               (def collar (max socket (slab-z 2 8)))\n\
               (assert clear (clearance hub collar 0.2) (accuracy 0.01))\n(assert too-near (clearance hub collar 0.3) (accuracy 0.01))\n";
    let c = run_claims(src);
    match &c[0].1 {
        ink::claims::Check { status: Status::Holds, evidence: Evidence::ClearanceSearch { .. }, .. } => {}
        x => panic!("{x:?}"),
    }
    assert_eq!(c[1].1.status, Status::Fails, "{:?}", c[1].1);
}

/// A profile drawn from x = 0, as every revolve is, has an edge on the
/// axis; turned, that edge is no surface: the axis is inside where the
/// profile touches it (here 3 from the wall, its distance), and outside
/// across a notch between two such edges.
#[test]
fn a_revolve_is_inside_along_its_axis() {
    use ink::claims::Status;
    let src = "(def rod (revolve-z (polygon (vec (vec 0 0) (vec 3 0) (vec 3 10) (vec 0 10)))))\n\
               (def notched (revolve-z (polygon (vec (vec 0 0) (vec 3 0) (vec 3 10) (vec 0 10) (vec 0 7) (vec 2 5) (vec 0 3)))))\n\
               (assert rod-axis (subset (probe 0 0 5) rod))\n\
               (assert rod-axis-depth (value (warp rod 0 0 5) < -2.999))\n\
               (assert notched-axis (subset (probe 0 0 1) notched))\n\
               (assert notch-on-axis (disjoint (probe 0 0 5) notched))\n";
    for (name, c) in run_claims(src) {
        assert_eq!(c.status, Status::Holds, "{name}: {c:?}");
    }
}

/// revolve-z leaves out a term only where its coefficients are proven: a
/// is y in real numbers, though its float coefficients round to 0, so
/// a - x is a wall, not an edge on the axis. The turned field is the
/// profile read at the radius, unchanged, and at (1, 0, 3), outside in
/// real numbers, never reads inside.
#[test]
fn a_revolve_keeps_a_wall_whose_floats_look_like_its_axis() {
    use ink::claims::Status;
    let src = "(def a (- (+ (* 1e16 y) y) (* 1e16 y)))\n\
               (def p (max (- a x) (- x 5) (- (abs y) 10)))\n\
               (def turned (revolve-z p))\n\
               (def same (- (revolve-z p) (warp p (len2 x y) z 0)))\n\
               (assert wall-kept (disjoint (probe 1 0 3) turned))\n";
    let (d, _) = field(src, "same");
    let Some(Value::Field(f, _)) = d.def("same") else { panic!("same is a field") };
    assert_eq!(d.g.as_const(*f), Some(0.0), "the turned field differs from the profile read at the radius");
    let c = &run_claims(src)[0].1;
    assert_ne!(c.status, Status::Fails, "{c:?}");
}

/// Axis terms grouped in a max of their own are left out too.
#[test]
fn a_revolve_leaves_out_grouped_axis_terms() {
    use ink::claims::Status;
    let src = "(def turned (revolve-z (max (max (neg x) (* -2 x)) (- x 3) (- (abs y) 5))))\n\
               (assert axis-in (subset (probe 0 0 0) turned))\n\
               (assert axis-depth (value (warp turned 0 0 0) < -2.999))\n";
    for (name, c) in run_claims(src) {
        assert_eq!(c.status, Status::Holds, "{name}: {c:?}");
    }
}

/// Axis terms wrapped in a name or written as neg of a min are left out
/// too: each reads 3 deep on the axis, as the bare max does.
#[test]
fn a_revolve_leaves_out_wrapped_axis_groups() {
    use ink::claims::Status;
    let src = "(def tagged (revolve-z (max (tag axis (max (neg x) (* -2 x))) (- x 3) (- (abs y) 5))))\n\
               (def negmin (revolve-z (max (neg (min x (* 2 x))) (- x 3) (- (abs y) 5))))\n\
               (assert tagged-depth (value (warp tagged 0 0 0) < -2.999))\n\
               (assert negmin-depth (value (warp negmin 0 0 0) < -2.999))\n";
    for (name, c) in run_claims(src) {
        assert_eq!(c.status, Status::Holds, "{name}: {c:?}");
    }
}

/// std's ways to write a profile reaching the axis: a polygon scaled (a
/// division under a positive factor) and polygon-distance (a select).
#[test]
fn a_revolve_of_std_profiles_is_inside_along_its_axis() {
    use ink::claims::Status;
    let src = "(def pts (vec (vec 0 0) (vec 3 0) (vec 3 10) (vec 0 10)))\n\
               (def scaled (revolve-z (scale (polygon pts) 2)))\n\
               (def exact (revolve-z (polygon-distance pts)))\n\
               (assert scaled-axis (subset (probe 0 0 10) scaled))\n\
               (assert exact-axis (subset (probe 0 0 5) exact))\n\
               (assert exact-depth (value (warp exact 0 0 5) < -2.999))\n";
    for (name, c) in run_claims(src) {
        assert_eq!(c.status, Status::Holds, "{name}: {c:?}");
    }
}

/// revolve-z reads each shared expression once: a40 = a39 + a39 ... over
/// -x is 2^40 terms written out, 41 nodes shared.
#[test]
fn a_revolve_reads_a_shared_profile_once() {
    let mut src = String::from("(def a0 (neg x))\n");
    for k in 1..=40 {
        src += &format!("(def a{k} (+ a{} a{}))\n", k - 1, k - 1);
    }
    src += "(def turned (revolve-z (max a40 (- x 3) (- (abs y) 5))))\n";
    let start = std::time::Instant::now();
    let (_, t) = field(&src, "turned");
    assert!(start.elapsed().as_secs_f64() < 5.0, "{:?}", start.elapsed());
    let v = eval_point(&t, &CTX, [0.0, 0.0, 0.0]);
    assert!(v < -2.999, "{v}");
}

/// A shared min tree is rewritten once per node, not once per path: a44,
/// the min of a43 and a42 down to x and 2x, is about 10^9 paths over 46
/// nodes.
#[test]
fn a_revolve_reads_a_shared_min_tree_once() {
    let mut binds = String::from("(a0 x) (a1 (* 2 x))");
    for k in 2..=44 {
        binds += &format!(" (a{k} (min a{} a{}))", k - 1, k - 2);
    }
    let src = format!("(def turned (let ({binds}) (revolve-z (neg (min a44 (- 3 x) (- 5 (abs y)))))))\n");
    let start = std::time::Instant::now();
    let (_, t) = field(&src, "turned");
    assert!(start.elapsed().as_secs_f64() < 2.0, "{:?}", start.elapsed());
    let v = eval_point(&t, &CTX, [0.0, 0.0, 0.0]);
    assert!(v < -2.999, "{v}");
}

#[test]
fn fits_beside_a_keyway_are_proven() {
    // A hub's bore 0.01 clear of a shaft with a
    // keyway, and a key 0.005 above the keyway's floor. The first-order
    // bound needs the facing surfaces' terms: through the offset (a shifted
    // max), the keyway subtracted (a negated max: at its floor-wall corner,
    // the least of the floor's and the wall's bounds), and the hub's union
    // of bore and slot subtracted (a negated min). Taken as one term, the
    // shaft's offset would have the hull of the cylinder's and the keyway's
    // slopes at the keyway's edge: proven only to tolerance.
    use ink::claims::Status;
    let src = "(def keyway (max (- (len2 x (- z (clamp z 89 111))) 4) (- 11 y)))\n\
               (def shaft (subtract (move (cylinder-z 15 40) 0 0 100) keyway))\n\
               (def hub (subtract (move (cylinder-z 35 24) 0 0 100) (min (- (len2 x y) 15.01) (max (- (abs x) 4.02) (neg y) (- y 18.3)))))\n\
               (def key (max (- (len2 x (- z (clamp z 90 110))) 3.99) (- 11.005 y) (- y 18.005)))\n\
               (assert hub-clears-shaft (disjoint (- shaft 0.005) hub))\n\
               (assert key-clears-shaft (disjoint (- shaft 0.004) key))\n";
    for (name, c) in run_claims(src) {
        assert_eq!(c.status, Status::Holds, "{name}: {c:?}");
    }
}

#[test]
fn clearance_search_grows_either_shape() {
    // A collar 0.25 inside a cup turned from a profile with no slope bound
    // at its surface (a sign times the root of a square), with a steep blob (slope 2) 0.3 above the cup's floor.
    // The collar's slope bound near both is 2, so its offset (0.4) passes
    // the 0.25 gap: no proof by emptiness. Covering the collar and growing
    // the cup needs the cup's slope near the collar: unknown. Covering the
    // cup and growing the collar, each cell by the collar's slope near it
    // (1 at the wall, 2 at the blob), proves it.
    use ink::claims::{Evidence, Status};
    let src = "(def bore (let ((g (move (rect-2d 2 4) 0 3 0))) (* (select g -1 1) (sqrt (square g)))))\n\
               (def cup (subtract (move (cylinder-z 2 4) 0 0 2) (revolve-z bore)))\n\
               (def collar (min (max (move (cylinder-z 0.75 4) 0 0 3) (slab-z 2 2.5)) (* 2 (move (sphere 0.2) 0 0 1.5))))\n\
               (assert clear (clearance cup collar 0.2) (accuracy 0.01))\n(assert too-near (clearance cup collar 0.27) (accuracy 0.01))\n";
    let c = run_claims(src);
    assert!(c[0].1.status == Status::Holds && matches!(c[0].1.evidence, Evidence::ClearanceSearch { .. }), "{:?}", c[0].1);
    assert_eq!(c[1].1.status, Status::Fails, "{:?}", c[1].1);
}

/// A direction must have a length: a rotation about a zero axis would be
/// 0 / 0 = 0 in Rodrigues' formula, a uniform scale by cos(angle), and a
/// unit sphere turned 60 degrees about (0 0 0) would read inside at 1.5.
/// Refused at load where the axis is constant or reads fixed params, at a
/// moment where it reads t, and over a scope unless proven nonzero there;
/// std's sweeps take a path's tangents and normals with (unit v).
#[test]
fn a_direction_without_a_length_is_refused() {
    let load = |src: &str, overrides: &[(String, ink::lang::Given)]| ink::load_str_with(&format!("(import \"std\")\n{src}"), "axis.ink", Path::new(ROOT), overrides);
    let err = |src: &str| load(src, &[]).err().unwrap_or_else(|| panic!("loads: {src}"));
    let ball = "(body ball (shape (sphere 1)))\n";
    // A constant axis, refused at the step.
    let e = err(&format!("{ball}(assembly root (instance b ball (at (rot 0 0 0 60))))"));
    assert!(e.message.contains("the axis of (rot 0 0 0 60) has zero length") && e.line == 3, "{e:?}");
    // A fixed param's value, overrides included, named with the form.
    let at_a = format!("(param a 0)\n{ball}(assembly root (instance b ball (at (rot a 0 0 60))))");
    let e = err(&at_a);
    assert!(e.message.contains("(rot a 0 0 60) at axis.ink:4 (assembly root ...) has zero length at a = 0"), "{e:?}");
    let mut d = load(&at_a, &[("a".into(), 2.0.into())]).unwrap();
    let ctx = d.ctx(&[("a".into(), 2.0.into())], 0.0).unwrap();
    let b = d.items("root/b").unwrap().remove(0);
    let world = d.placed_field(&b);
    assert!((eval_point(&Tape::new(&d.g, world), &ctx, [1.5, 0.0, 0.0]) - 0.5).abs() < 1e-12);
    // A solve's unknown at its solution, once every solve is read (one
    // read after the placement sets a here).
    let solved = |root: f64| format!("(param a 0)\n{ball}(assembly root (instance b ball (at (rot a 0 0 60))))\n(solve s (unknowns a) (== a {root}))");
    load(&solved(1.0), &[]).unwrap();
    let e = err(&solved(0.0));
    assert!(e.message.contains("has zero length at a = 0"), "{e:?}");
    // Over ranged params: refused unless proven nonzero over the band.
    let ranged = |tol: &str, axis: &str| {
        let d = load(&format!("(param a 1 (tol {tol}))\n{ball}(assembly root (instance b ball (at (rot {axis} 60))))"), &[]).unwrap();
        let ctx = d.ctx(&[], 0.0).unwrap();
        let mut v = ink::view::ModelView::new(d, ctx);
        let spans = v.declared_spans().unwrap();
        v.set_param_spans(&spans)
    };
    // The zero itself, at a = 0: its squared length bounds to exactly 0
    // there (an exact result is not widened).
    let e = ranged("-2 2", "a 0 0").unwrap_err();
    assert!(e.message.contains("(rot a 0 0 60) at axis.ink:4 (assembly root ...) has zero length at a = 0 (squared length 0), within a in [-1, 3]"), "{e:?}");
    ranged("-0.5 2", "a 0 0").unwrap();
    ranged("-1 6", "(cos a) (sin a) 0").unwrap();
    // Over time: at each moment, and over a motion (by steps).
    let d = load(&format!("{ball}(assembly root (instance b ball (at (rot (- t 1) 0 0 60))))"), &[]).unwrap();
    assert!(d.ctx(&[], 1.0).unwrap_err().contains("has zero length at t = 1"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(d, ctx);
    let e = v.set_motion(Some((0.0, 2.0))).unwrap_err();
    assert!(e.message.contains("has zero length at t = 1 (squared length 0), within t in [0, 2]"), "{e:?}");
    v.set_motion(Some((1.5, 3.0))).unwrap();
    // std's frames: two equal points, and up along the first segment.
    let e = err("(body p (shape (sweep (circle-2d 1) (vec (vec 0 0 0) (vec 0 0 10) (vec 0 0 10) (vec 10 0 10)) (vec 1 0 0))))");
    assert!(e.message.contains("(unit (- b a)) has zero length") && e.called_from.iter().any(|c| c.0 == "sweep" && c.2 == 2), "{e:?}");
    let e = err("(body p (shape (sweep (circle-2d 1) (vec (vec 0 0 0) (vec 0 0 10) (vec 10 0 10)) (vec 0 0 1))))");
    assert!(e.message.contains("(unit (cross (cross (nth tans 0) up) (nth tans 0))) has zero length"), "{e:?}");
    // Over space, a field: its zero is an undefined case, as any division.
    load("(def f (unit (vec x y z)))", &[]).unwrap();
}

/// Fail fast: best first, a box undecided at the finest cells ends the
/// search. Two plates touching over a 20 x 20 face (at z = 1/3, off the
/// cells' faces) are apart exactly: the boxes across the face split there.
/// The upper one's height read as (z + 3 y) - 3 y, the same in real numbers
/// but sharing nothing with the lower's, they stay undecided there at any
/// res, asked exactly: the search ends at its first box on the
/// face at the finest cells (Ended::Finest), the rest of the face not
/// refined (all of it is four million boxes at 0.01), and it holds to
/// about that box's size (accuracy attained; their intersection is empty).
/// To the accuracy they hold, and a bump into the
/// lower plate at a corner of the face is found, 0.3 mm deep or only 0.104
/// (just over twice the accuracy, 0.05). Two boxes overlapping by 0.01
/// asked to 0.005: unknown, holding to
/// no less than the 0.005 the overlap's middle is deep (the boxes left,
/// coarse where the search ended early, bound it by their samples); to
/// 0.05 they hold.
#[test]
fn a_search_ends_where_a_box_stays_undecided_at_its_finest_cells() {
    use ink::query::{Ended, Verdict, prove_empty};
    use std::sync::Arc;
    let at = |zc: f64, bump: f64| {
        let src = format!("(def f (max (move (box 20 20 4) 0 0 (- {zc} 2)) (min (move (box 20 20 4) 0 0 (+ {zc} 2)) (move (box 2 2 {}) 8 8 {zc}))))", 2.0 * bump);
        Arc::new(field(&src, "f").1)
    };
    let zc = 1.0 / 3.0;
    let shape = |bump: f64| at(zc, bump);
    let region = ([-11.0; 3], [11.0; 3]);
    assert!(matches!(prove_empty(&at(0.0, 0.0), &CTX, region.0, region.1, None, 0.01, 0.0).verdict, Verdict::Proven));
    assert!(matches!(prove_empty(&shape(0.0), &CTX, region.0, region.1, None, 0.01, 0.0).verdict, Verdict::Proven));
    let touch = Arc::new(field(&format!("(def f (max (move (box 20 20 4) 0 0 (- {zc} 2)) (warp (move (box 20 20 4) 0 0 (+ {zc} 2)) x y (- (+ z (* 3 y)) (* 3 y)))))"), "f").1);
    let s = prove_empty(&touch, &CTX, region.0, region.1, None, 0.01, 0.0);
    let Verdict::Unknown { cells, ended: Ended::Finest, cell: Some([lo, hi]), attained: Some(a), .. } = s.verdict else { panic!("{:?}", s.verdict) };
    assert!(cells < 100_000 && a <= 1.0, "{cells} boxes left, holds to {a}");
    assert!(lo[2] <= zc && zc <= hi[2] && hi[0] - lo[0] <= 0.0100001, "{lo:?} {hi:?}");
    assert!(matches!(prove_empty(&touch, &CTX, region.0, region.1, None, 0.025, 0.05).verdict, Verdict::Proven));
    for bump in [0.3, 0.104] {
        let v = prove_empty(&shape(bump), &CTX, region.0, region.1, None, 0.025, 0.05).verdict;
        let Verdict::Disproven { witness, .. } = v else { panic!("a bump {bump} deep: {v:?}") };
        assert!(witness[0] > 7.0 && witness[1] > 7.0 && witness[2] < zc, "{witness:?}");
    }
    let (_, both) = field("(def f (max (box 10 10 10) (move (box 10 10 10) 9.99 0 0)))", "f");
    let both = Arc::new(both);
    let v = prove_empty(&both, &CTX, [-6.0; 3], [16.0; 3], None, 0.0025, 0.005).verdict;
    let Verdict::Unknown { attained: Some(a), .. } = v else { panic!("{v:?}") };
    assert!((0.005..0.1).contains(&a), "holds to {a}");
    assert!(matches!(prove_empty(&both, &CTX, [-6.0; 3], [16.0; 3], None, 0.025, 0.05).verdict, Verdict::Proven));
}

/// A def's value at the default params: a field or a vec of fields.
fn nums(d: &Doc, name: &str) -> Vec<f64> {
    let v = |f: ink::graph::NodeId| eval_point(&Tape::new(&d.g, f), &Ctx::new(0.0, d.params.iter().map(|p| p.default).collect()), [0.0; 3]);
    match d.def(name).unwrap_or_else(|| panic!("no def {name}")) {
        Value::Field(f, _) => vec![v(*f)],
        Value::Vec(items) => items.iter().map(|i| if let Value::Field(f, _) = i { v(*f) } else { panic!("{name}: not fields") }).collect(),
        _ => panic!("{name}: not a field or vec"),
    }
}

fn near(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-9)
}

#[test]
fn frames_compose_through_assemblies() {
    // A frame declared on a body, exported by an assembly, read from the
    // world: origin and x axis by hand. (at F (move ...)) moves along the
    // outer axes, (at (move ...) F) along F's own.
    let src = "(body plate (shape (box 10 10 2)) (frame top (at (move 0 0 1))))
(fn put (d) (let ((h (* 2 d))) (at (rot-z 90) (move h 0 0))))
(assembly a (instance p plate (put 2.5)) (frame t (frame-of p top)))
(assembly root (instance q a (at (move 0 0 10))))
(def o (get (frame-of root/q t) \"origin\"))
(def xa (get (frame-of root/q t) \"x\"))
(def own (get (at (move 1 0 0) (frame-of root/q t)) \"origin\"))
(def outer (get (at (frame-of root/q t) (move 1 0 0)) \"origin\"))
(def local (get (frame-of plate top) \"origin\"))
(def placed (place (box 2 2 2) (frame-of root/q t)))";
    let (d, t) = field(src, "placed");
    assert!(near(&nums(&d, "o"), &[5.0, 0.0, 11.0]), "{:?}", nums(&d, "o"));
    assert!(near(&nums(&d, "xa"), &[0.0, 1.0, 0.0]), "{:?}", nums(&d, "xa"));
    assert!(near(&nums(&d, "own"), &[5.0, 1.0, 11.0]), "{:?}", nums(&d, "own"));
    assert!(near(&nums(&d, "outer"), &[6.0, 0.0, 11.0]), "{:?}", nums(&d, "outer"));
    assert!(near(&nums(&d, "local"), &[0.0, 0.0, 1.0]));
    // The box sits centred on the frame's origin.
    assert!(eval_point(&t, &CTX, [5.0, 0.9, 11.0]) < 0.0 && eval_point(&t, &CTX, [5.0, 1.1, 11.0]) > 0.0);
    for (bad, msg) in [
        ("(def f (frame-of root/q nope))", "'a' has no frame 'nope'. Its frames: t"),
        ("(def f (frame-of root/q/p/x top))", "'plate' is a body; 'root/q/p/x' goes through it"),
        ("(def f (get (frame-of plate top) \"w\"))", "a frame has no 'w'"),
        ("(body b (shape (box 1 1 1)) (frame f (box 1 1 1)))", "F must be a frame, not a field"),
        ("(body b (shape (box 1 1 1)) (frame f (at)) (frame f (at)))", "duplicate frame 'f'"),
        ("(body b (shape (box 1 1 1)) (frame f (at (move x 0 0))))", "placements cannot depend on x, y or z"),
        ("(assembly c (instance a plate (frame-of b top)) (instance b plate))", "'b' is placed later in 'c': move its clause above this one"),
    ] {
        let e = ink::load_str(&format!("(import \"std\")\n{src}\n{bad}"), "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{bad} loads"));
        assert!(e.message.contains(msg), "{bad}: {e:?}");
    }
}

#[test]
fn select_chooses_between_frames() {
    // An agent's robot arm: a marker rests at a point until t = 1, then
    // rides a frame of the assembly (a part picked up), its placement
    // (select c F G) as written there.
    let src = "(body plate (shape (box 80 60 6)) (frame seat (at (rot-z 90) (move 0 0 3))))
(body marker (shape (sphere 3)))
(assembly root
  (instance plate plate)
  (instance mark marker (select (- t 1) (at (move 40 0 20)) (frame-of plate seat))))
(def o (get (frame-of root/mark) \"origin\"))
(def xa (get (frame-of root/mark) \"x\"))";
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let at = |name: &str, t: f64| -> Vec<f64> {
        let Some(Value::Vec(items)) = d.def(name) else { panic!("{name}: not a vec") };
        items.iter().map(|i| if let Value::Field(f, _) = i { eval_point(&Tape::new(&d.g, *f), &Ctx::new(t, Vec::new()), [0.0; 3]) } else { panic!() }).collect()
    };
    assert!(near(&at("o", 0.5), &[40.0, 0.0, 20.0]) && near(&at("xa", 0.5), &[1.0, 0.0, 0.0]), "{:?} {:?}", at("o", 0.5), at("xa", 0.5));
    assert!(near(&at("o", 1.5), &[0.0, 0.0, 3.0]) && near(&at("xa", 1.5), &[0.0, 1.0, 0.0]), "{:?} {:?}", at("o", 1.5), at("xa", 1.5));
    let e = ink::load_str("(import \"std\")\n(def f (select 1 (at) (vec 1 2 3)))", "test", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("'select' expects numbers or vecs, or a number and two frames, got a frame"), "{e:?}");
    // A choice by x is no placement: the body would be in both places at
    // once. Refused where the frame is used, however it was made.
    let e = ink::load_str(&format!("(import \"std\")\n{src}\n(assembly two (instance a marker (select (- x 1) (at) (at (move 5 0 0)))))"), "test", Path::new(ROOT)).err().unwrap();
    assert!(e.message.contains("placements cannot depend on x, y or z") && e.line == 9, "{e:?}");
}

#[test]
fn a_top_level_form_inside_another_says_where_it_goes() {
    // As agents write them: (t-span ...) in an assembly; (def t-span
    // (t-span 0 2)).
    for (src, inside) in [
        ("(body b (shape (box 1 1 1)))\n(assembly root (t-span 0 2) (instance a b))", "an assembly"),
        ("(body b (shape (box 1 1 1)) (param w 3))", "a body"),
        ("(def t-span (t-span 0 2))", "an expression"),
    ] {
        let e = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{src} loads"));
        let head = if src.contains("param") { "param" } else { "t-span" };
        let want = format!("({head} ...) is a top-level form. Write it at the file's top level, not inside {inside}");
        assert!(e.hint.as_deref().is_some_and(|h| h.starts_with(&want)), "{src}: {e:?}");
    }
}

#[test]
fn axes_writes_a_frame_as_ink_prints_one() {
    // (axes (origin ...) (x ...) (z ...)) is the frame with those
    // directions: z kept, x made square to it, y = z cross x. Quarter
    // turns are exact, so a turn written either way is the same numbers,
    // bit for bit; a tilt read from a param is the turn to rounding.
    let src = "(param tilt 65)
(def a (axes (origin 0 6 20) (x 1 0 0) (z 0 1 0)))
(def b (at (rot-x -90) (move 0 6 20)))
(def flip (at (rot-x 180)))
(def skew (axes (z 0 0 2) (origin 1 2 3) (x 1 0 1)))
(def tilted (axes (origin 0 0 0) (x 1 0 0) (z 0 (neg (sin (deg->rad tilt))) (cos (deg->rad tilt)))))
(def turned (at (rot-x tilt)))";
    let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = Ctx::new(0.0, d.params.iter().map(|p| p.default).collect());
    let frame = |name: &str| -> Vec<f64> {
        let Some(Value::Frame(xf)) = d.def(name) else { panic!("{name}: not a frame") };
        xf.r.iter().flatten().chain(&xf.t).map(|&e| eval_point(&Tape::new(&d.g, e), &ctx, [0.0; 3])).collect()
    };
    let bits = |v: Vec<f64>| v.into_iter().map(f64::to_bits).collect::<Vec<_>>();
    assert_eq!(bits(frame("a")), bits(frame("b")));
    assert!(frame("flip").iter().all(|v| [0.0, 1.0, -1.0].contains(v)), "{:?}", frame("flip"));
    assert_eq!(frame("skew"), [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 3.0]);
    assert!(near(&frame("tilted"), &frame("turned")), "{:?} {:?}", frame("tilted"), frame("turned"));
    for (bad, msg) in [
        ("(def f (axes (origin 0 0 0) (x 0 0 1) (z 0 0 -3)))", "square to its z, has zero length"),
        ("(def f (axes (origin 0 0 0) (x 1 0 0) (z 0 0 0)))", "has zero length"),
        ("(def f (axes (origin 0 0 0) (x 1 0 0) (y 0 1 0) (z 0 0 1)))", "takes no y"),
        ("(def f (axes (origin 0 0 0) (z 0 0 1)))", "needs origin, x and z"),
        ("(def f (axes (origin 0 0) (x 1 0 0) (z 0 0 1)))", "(origin ...) takes 3 numbers"),
        ("(def f (axes (origin 0 0 0) (x 1 0 0) (x 1 0 0) (z 0 0 1)))", "(x ...) is given twice"),
        ("(def f (place (sphere 1) (axes (origin x 0 0) (x 1 0 0) (z 0 0 1))))", "placements cannot depend on x, y or z"),
    ] {
        let e = ink::load_str(&format!("(import \"std\")\n{bad}"), "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{bad} loads"));
        assert!(e.message.contains(msg), "{bad}: {e:?}");
    }
}

#[test]
fn frame_of_reads_local_instance_names_in_a_namespace() {
    // Two files define `bracket`; in each, an instance named like a
    // top-level body reads the instance, not the body, in a frame-of path,
    // and a WHERE reading a def named like an instance reads the def.
    let dir = crate::common::Scratch::new("frames-ns");
    for (f, h) in [("a.ink", 2.0), ("b.ink", 7.0)] {
        let src = format!(
            "(import \"std\")
(body bracket (shape (box 10 10 {h})) (frame top (at (move 0 0 {}))))
(body base (shape (box 1 1 1)))
(assembly arm (instance base bracket (at (move 0 0 1))) (frame tip (frame-of base top)))
(def lift 3)
(assembly arm2 (instance lift bracket (at (move 0 0 lift))) (frame tip (frame-of lift top)))
(assembly row (for i (vec 1) (for lift (vec 2) (instance bracket bracket (at (move 0 0 (* i lift)))))) (frame tip (frame-of bracket top)))
(def w 5)
(body (rod (w 1)) (shape (box w w w)) (frame top (at (move 0 0 w))))
(assembly rods (instance r (rod 3)) (frame tip (frame-of r top)))",
            h / 2.0
        );
        std::fs::write(dir.join(f), src).unwrap();
    }
    std::fs::write(
        dir.join("top.ink"),
        "(import \"a.ink\" as a) (import \"b.ink\" as b)
(def ta (get (frame-of a.arm tip) \"origin\"))
(def tb (get (frame-of b.arm tip) \"origin\"))
(def tl (get (frame-of b.arm2 tip) \"origin\"))
(def tr (get (frame-of a.row tip) \"origin\"))
(def tp (get (frame-of b.row/bracket top) \"origin\"))
(def tw (get (frame-of a.rods tip) \"origin\"))
(def tv (get (frame-of a.rod top) \"origin\"))",
    )
    .unwrap();
    let d = ink::load_file_with(&dir.join("top.ink"), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    assert!(near(&nums(&d, "ta"), &[0.0, 0.0, 2.0]), "{:?}", nums(&d, "ta"));
    assert!(near(&nums(&d, "tb"), &[0.0, 0.0, 4.5]), "{:?}", nums(&d, "tb"));
    assert!(near(&nums(&d, "tl"), &[0.0, 0.0, 6.5]), "{:?}", nums(&d, "tl"));
    // Instances declared in (nested) `for` clauses are local names too, and
    // a loop variable shadows a def (not: `b.row` has no instance 'bracket').
    assert!(near(&nums(&d, "tr"), &[0.0, 0.0, 3.0]), "{:?}", nums(&d, "tr"));
    assert!(near(&nums(&d, "tp"), &[0.0, 0.0, 5.5]), "{:?}", nums(&d, "tp"));
    // A body's arguments are local names too, not the file's def `w`.
    assert!(near(&nums(&d, "tw"), &[0.0, 0.0, 3.0]), "{:?}", nums(&d, "tw"));
    assert!(near(&nums(&d, "tv"), &[0.0, 0.0, 1.0]), "{:?}", nums(&d, "tv"));
}

/// `nums` with params set by name.
fn nums_at(d: &Doc, name: &str, set: &[(&str, f64)]) -> Vec<f64> {
    let mut params: Vec<f64> = d.params.iter().map(|p| p.default).collect();
    for (n, x) in set {
        params[d.param_index(n).unwrap()] = *x;
    }
    let v = |f: ink::graph::NodeId| eval_point(&Tape::new(&d.g, f), &Ctx::new(0.0, params.clone()), [0.0; 3]);
    match d.def(name).unwrap_or_else(|| panic!("no def {name}")) {
        Value::Field(f, _) => vec![v(*f)],
        Value::Vec(items) => items.iter().map(|i| if let Value::Field(f, _) = i { v(*f) } else { panic!() }).collect(),
        _ => panic!("{name}: not a field or vec"),
    }
}

#[test]
fn frames_place_parts() {
    // A motor seated under a turned, moved bracket by its mount frame,
    // turned 45 about the bracket's underside normal (which points down).
    // Motor x goes FLIP -> Rz(45) -> Rx(180) -> Rz(30): -15 degrees in the
    // world. Its side connector (local (23, 0, -7.5)) is the off-axis probe:
    // a missing FLIP puts the motor above the plate, a reversed turn puts
    // the connector at +75 degrees, clear of everything.
    let src = "(param T 5)
(body bracket (shape (subtract (box 60 60 T) (cylinder-z 11.5 (+ T 2))))
  (frame under (at (rot-x 180) (move 0 0 (* -0.5 T)))))
(body motor (shape (min (min (move (box 42 42 40) 0 0 -20) (cylinder-z 2.5 40)) (move (box 4 4 5) 23 0 -7.5)))
  (frame mount (at)) (frame shaft (at (move 0 0 20))))
(assembly drive
  (instance b bracket (at (rot-z 30) (move 40 10 0)))
  (instance m motor (frame-of b under) mount (turn z 45)))
(assembly root (instance d drive))
(def o (get (frame-of root/d/m shaft) \"origin\"))
(def xa (get (frame-of root/d/m shaft) \"x\"))
(def za (get (frame-of root/d/m shaft) \"z\"))
(def all (shape-of root))";
    let (d, t) = field(src, "all");
    let r = (-15f64).to_radians();
    assert!(near(&nums(&d, "o"), &[40.0, 10.0, 17.5]), "{:?}", nums(&d, "o"));
    assert!(near(&nums(&d, "xa"), &[r.cos(), r.sin(), 0.0]), "{:?}", nums(&d, "xa"));
    assert!(near(&nums(&d, "za"), &[0.0, 0.0, 1.0]), "{:?}", nums(&d, "za"));
    let at = |p: [f64; 3]| eval_point(&t, &Ctx::new(0.0, vec![5.0]), p);
    assert!(at([40.0, 10.0, 16.5]) < 0.0 && at([40.0, 10.0, 18.5]) > 0.0, "shaft tip");
    let probe = |deg: f64| [40.0 + 23.0 * deg.to_radians().cos(), 10.0 + 23.0 * deg.to_radians().sin(), -10.0];
    assert!(at(probe(-15.0)) < 0.0 && at(probe(75.0)) > 0.0, "connector");
    // The bracket's underside meets the mount face: the motor's top is at z = -2.5.
    assert!(at([40.0 + 15.0, 10.0, -2.6]) < 0.0 && at([40.0 + 15.0, 10.0, -2.4]) < 0.0 && at([40.0 + 15.0, 10.0, 2.6]) > 0.0);
}

#[test]
fn a_seat_meets_rotated_and_translated_frames() {
    // C = Rx(90) at (1,2,3) met with WHERE = Ry(90) then Rz(90) at
    // (10,20,30): (x, y, z) -> (z + 7, y + 18, 31 - x). (turn z 90) sends
    // C's x to WHERE's y; (slide x 4) puts C's origin at (10, 20, 26).
    let src = "(body c (shape (box 1 1 1)) (frame f (at (rot-x 90) (move 1 2 3))))
(def w (at (rot-y 90) (rot-z 90) (move 10 20 30)))
(assembly seat (instance i c w f))
(assembly turned (instance i c w f (turn z 90)))
(assembly slid (instance i c w f (slide x 4)))
(def p (get (at (move 5 6 7) (frame-of seat/i)) \"origin\"))
(def cx (get (frame-of turned/i f) \"x\"))
(def wy (get w \"y\"))
(def cz (get (frame-of seat/i f) \"z\"))
(def wz (get w \"z\"))
(def so (get (frame-of slid/i f) \"origin\"))
(def all (shape-of seat))";
    let (d, _) = field(src, "all");
    assert!(near(&nums(&d, "p"), &[14.0, 24.0, 26.0]), "{:?}", nums(&d, "p"));
    assert!(near(&nums(&d, "cx"), &nums(&d, "wy")), "{:?}", nums(&d, "cx"));
    assert!(near(&nums(&d, "cz"), &nums(&d, "wz").iter().map(|v| -v).collect::<Vec<_>>()));
    assert!(near(&nums(&d, "so"), &[10.0, 20.0, 26.0]), "{:?}", nums(&d, "so"));
}

#[test]
fn a_reused_sub_assembly_exports_its_tool_frame() {
    // The same segment twice, each turned on the previous one's tip:
    // tool = (200 cos s, 200 sin s, 44), its x along the arm.
    let src = "(body beam (shape (move (box 116 16 6) 50 0 0))
  (frame hub (at (rot-x 180) (move 0 0 -3))) (frame tip (at (move 100 0 3))))
(body post (shape (cylinder-z 10 40)) (frame top (at (move 0 0 20))))
(body pin (shape (cylinder-z 4 6)) (frame seat (at (rot-x 180) (move 0 0 -3))) (frame top (at (move 0 0 3))))
(assembly segment
  (instance b beam)
  (instance pin pin (frame-of b tip) seat)
  (frame base (frame-of b hub))
  (frame tip (frame-of pin top)))
(param shoulder 30)
(assembly arm
  (instance p post)
  (instance s1 segment (frame-of p top) base (turn z shoulder) (limits -120 120))
  (instance s2 segment (frame-of s1 tip) base (turn z (* -40 t)) (limits -150 150))
  (frame tool (frame-of s2 tip)))
(assembly root (instance a arm))
(def o (get (frame-of root/a tool) \"origin\"))
(def xa (get (frame-of root/a tool) \"x\"))
(def all (shape-of root))";
    let (d, t) = field(src, "all");
    for s in [0.0f64, 90.0] {
        let (c, sn) = (f64::cos(s.to_radians()), f64::sin(s.to_radians()));
        assert!(near(&nums_at(&d, "o", &[("shoulder", s)]), &[200.0 * c, 200.0 * sn, 44.0]), "{s}: {:?}", nums_at(&d, "o", &[("shoulder", s)]));
        assert!(near(&nums_at(&d, "xa", &[("shoulder", s)]), &[c, sn, 0.0]));
    }
    // Pin material just under the tool frame, at the default shoulder.
    let r = 30f64.to_radians();
    let (c, s) = (r.cos(), r.sin());
    assert!(eval_point(&t, &Ctx::new(0.0, vec![30.0]), [200.0 * c, 200.0 * s, 43.9]) < 0.0);
    assert!(eval_point(&t, &Ctx::new(0.0, vec![30.0]), [200.0 * c, 200.0 * s, 44.1]) > 0.0);
}

#[test]
fn joint_limits_are_claims() {
    use ink::claims::{ClaimOptions, Status::*, run_asserts};
    let src = "(param s 100 (tol -0.5 0.5))
(body rail (shape (move (box 270 12 8) 135 0 4)) (frame top (at (move 0 0 8))))
(body car (shape (move (box 45 27 10) 0 0 5)) (frame base (at (rot-x 180))))
(assembly axis (instance r rail) (instance c car (frame-of r top) base (slide x s) (limits 25 245)))
(assembly moving (instance r rail) (instance c car (frame-of r top) base (slide x (* 300 t)) (limits 25 245)))
(assembly root (instance a axis))
(def o (get (frame-of root/a/c base) \"origin\"))
(def all (shape-of root))";
    let (d, t) = field(src, "all");
    // The carriage's underside sits on the rail top, s along it.
    assert!(near(&nums(&d, "o"), &[100.0, 0.0, 8.0]), "{:?}", nums(&d, "o"));
    assert!(eval_point(&t, &Ctx::new(0.0, vec![100.0]), [100.0 + 20.0, 0.0, 12.0]) < 0.0);
    assert!(eval_point(&t, &Ctx::new(0.0, vec![100.0]), [100.0 + 24.0, 0.0, 12.0]) > 0.0);
    let run = |set: f64, span: Option<(f64, f64)>, ranged: bool, name: &str| {
        let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
        let mut ctx = d.ctx(&[("s".into(), set.into())], 0.0).unwrap();
        if ranged {
            ctx.p_spans = ink::eval::ParamSpans::new([(d.param_index("s").unwrap() as u32, (set - 0.5, set + 0.5))]);
        }
        let mut v = ink::view::ModelView::new(d, ctx.into_scope());
        v.set_motion(span).unwrap();
        run_asserts(&mut v, Some(name), false, ClaimOptions::default()).unwrap().asserts.remove(0).check
    };
    assert_eq!(run(100.0, None, true, "limits:axis/c").status, Holds);
    let c = run(260.0, None, false, "limits:axis/c");
    assert!(c.status == Fails && matches!(c.evidence, ink::claims::Evidence::Value(v) if (v + 15.0).abs() < 1e-9), "{c:?}");
    let c = run(100.0, Some((0.0, 1.0)), false, "limits:moving/c");
    let v = match &c.evidence {
        ink::claims::Evidence::Value(v) => *v,
        ink::claims::Evidence::ValueSpan(s) => s.witness.as_ref().unwrap().1,
        _ => panic!("{c:?}"),
    };
    let t = c.motion.as_ref().and_then(|m| m.motion.t).unwrap();
    assert!(c.status == Fails && v < 0.0 && !(25.0 / 300.0..=245.0 / 300.0).contains(&t), "{c:?}");
    for (bad, msg) in [
        ("(assembly e (instance r rail) (instance c car (frame-of r top) bse))", "'car' has no frame 'bse'. Its frames: base"),
        ("(body plain (shape (box 1 1 1))) (assembly e (instance r rail) (instance c plain (frame-of r top) bse))", "Add (frame bse F) to it"),
        ("(assembly e (instance r rail) (instance c car (frame-of r top) base (turn z 1) (slide x 2)))", "one motion per instance"),
        ("(assembly e (instance r rail) (instance c car (frame-of r top) base (limits 0 1)))", "(limits LO HI) bounds a motion"),
        ("(assembly e (instance r rail) (instance c car (frame-of r top) base (slide x s) (limits 5 1)))", "limits 5 > 1"),
        ("(assembly e (instance r rail) (instance c car (frame-of r top) base (slide w s)))", "must be x, y or z"),
        ("(assembly e (instance r rail) (instance c car (frame-of r top) base (spin z s)))", "bad instance clause"),
    ] {
        let e = ink::load_str(&format!("(import \"std\")\n{src}\n{bad}"), "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{bad} loads"));
        assert!(e.message.contains(msg), "{bad}: {e:?}");
    }
}

#[test]
fn assemblies_take_arguments() {
    use ink::claims::{ClaimOptions, Status::*, run_asserts};
    // A motor whose rotor's turn is a parameter, driven by its instances:
    // a default (still), a spin-up from a def of t, a gear pair's ratio of
    // it, and a fixed angle passed through a parent's parameter.
    let src = "(body housing (shape (cylinder-z 10 20)) (frame axis (at (move 0 0 10))))
(body rotor (shape (move (box 2 8 4) 0 0 2)) (frame base (at (rot-x 180))) (frame tip (at (move 6 0 4))))
(def top 6000)
(assembly (motor (angle 0) (n 1))
  (instance h housing)
  (for i (range n) (instance r$i rotor (at (move 0 0 (* 4 i)) (frame-of h axis)) base (turn z angle) (speed-limits 0 (* top 1rpm))))
  (frame tip (frame-of r0 tip)))
(assembly (pair (a 0)) (instance m (motor a)) (instance g (motor (* -2 a) 2) (at (move 40 0 0))))
(param top-rpm 9000)
(def angle 7)
(def spin (* 0.5 (* top-rpm 1rpm) (square t)))
(assembly root
  (instance still motor)
  (instance zero (motor 0))
  (instance m (motor spin))
  (instance q (pair 30) (at (move 0 50 0))))
(def m-tip (get (frame-of root/m tip) \"origin\"))
(def q-tip (get (frame-of root/q/g tip) \"origin\"))
(def still-tip (get (frame-of root/still tip) \"origin\"))
(def all (shape-of root))";
    let (d, _) = field(src, "all");
    // Equal values share one assembly (bodies are shared by all of them).
    let child = |a: &str, i: &str| d.assembly(a).unwrap().instances.iter().find(|x| x.name == i).unwrap().child.clone();
    assert_eq!((child("root", "still"), child("root", "zero")), ("motor".into(), "motor".into()));
    assert_eq!((child("root", "m"), child("root", "q"), child("root/q", "g")), ("root/m".into(), "root/q".into(), "root/q/g".into()));
    assert_eq!(d.assembly("root/q/g").unwrap().instances.len(), 3, "n = 2: the housing and two rotors");
    assert_eq!(d.assembly("root/q").unwrap().params, vec![("a".to_string(), "30".to_string())]);
    // The rotor turns by its argument (the file's own `angle` is not the parameter).
    let at = |name: &str, t: f64| {
        let Some(Value::Vec(v)) = d.def(name) else { panic!() };
        v.iter().map(|f| if let Value::Field(f, _) = f { eval_point(&Tape::new(&d.g, *f), &Ctx::new(t, vec![9000.0]), [0.0; 3]) } else { panic!() }).collect::<Vec<_>>()
    };
    assert!(near(&at("still-tip", 0.3), &[6.0, 0.0, 14.0]), "{:?}", at("still-tip", 0.3));
    let deg = 0.5f64 * 54000.0 * 0.01;
    let r = deg.to_radians();
    assert!(near(&at("m-tip", 0.1), &[6.0 * r.cos(), 6.0 * r.sin(), 14.0]), "{:?}", at("m-tip", 0.1));
    let r = (-60f64).to_radians();
    assert!(near(&at("q-tip", 0.0), &[40.0 + 6.0 * r.cos(), 50.0 + 6.0 * r.sin(), 14.0]), "{:?}", at("q-tip", 0.0));
    // The speed limit is a claim on the drive's rate, over the motion: the
    // spin-up passes 6000 rpm at t = 2/3 of a second (at top-rpm 9000).
    let run = |top: f64, name: &str| {
        let d = ink::load_str(&format!("(import \"std\")\n{src}"), "test", Path::new(ROOT)).unwrap();
        let ctx = d.ctx(&[("top-rpm".into(), top.into())], 0.0).unwrap();
        let mut v = ink::view::ModelView::new(d, ctx.into_scope());
        v.set_motion(Some((0.0, 1.0))).unwrap();
        run_asserts(&mut v, Some(name), false, ClaimOptions::default()).unwrap().asserts.remove(0).check
    };
    // From rest: the rate starts on its floor, 0, so it holds to tolerance there.
    assert_eq!(run(5000.0, "speed-limits:root/m/r0").status, Holds);
    assert_eq!(run(5000.0, "speed-limits:motor/r0").status, Holds);
    let c = run(9000.0, "speed-limits:root/m/r0");
    let t = c.motion.as_ref().and_then(|m| m.motion.t).unwrap();
    assert!(c.status == Fails && (2.0 / 3.0..=1.0).contains(&t), "{c:?}");
    assert_eq!(run(5000.0, "speed-limits:root/q/g/r0").status, Holds, "a fixed angle does not move");
    for (bad, msg) in [
        ("(assembly e (instance m (motor 1 2 3)))", "'motor' takes 2 argument(s) (angle n), got 3"),
        ("(assembly e (instance h (housing 1)))", "'housing' is a body without arguments"),
        ("(assembly e (instance h (root 1)))", "'root' is an assembly without arguments"),
        ("(assembly (e angle) (instance h housing))", "an argument is (name default)"),
        ("(assembly e (instance m (motor x)))", "not fields of x y z"),
        ("(assembly e (instance h housing) (instance r rotor (frame-of h axis) base (speed-limits 0 1)))", "(speed-limits LO HI) bounds a motion"),
        ("(assembly e (instance h housing) (instance r rotor (frame-of h axis) base (turn z t) (speed-limits 2 1)))", "speed-limits 2 > 1"),
    ] {
        let e = ink::load_str(&format!("(import \"std\")\n{src}\n{bad}"), "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{bad} loads"));
        assert!(e.message.contains(msg), "{bad}: {e:?}");
    }
}

#[test]
fn bodies_take_arguments() {
    // A rod whose length is an argument, as an assembly's: its shape, an
    // attribute and a frame read it. A slider's argument passes on to it,
    // so the rod's shape reads t where a parent drives the slider with t.
    let src = "(body (rod (h 10) (w 2))
  (shape (move (box w w h) 0 0 (* 0.5 h)))
  (attr reach h)
  (frame top (at (move 0 0 h))))
(assembly (slider (s 0)) (instance r (rod (+ 20 s))) (frame tip (frame-of r top)))
(def h 99)
(assembly root
  (instance a (slider (* 10 t)))
  (instance plain rod (at (move 10 0 0)))
  (instance ten (rod 10) (at (move 20 0 0)))
  (instance c (rod 30) (at (move 30 0 0)))
  (instance d (rod 30) (at (move 40 0 0))))
(def tip (get (frame-of root/a tip) \"origin\"))
(def all (shape-of root))";
    let (d, all) = field(src, "all");
    let child = |a: &str, i: &str| d.assembly(a).unwrap().instances.iter().find(|x| x.name == i).unwrap().child.clone();
    // The defaults, given or not, are the body as defined; equal values
    // share one body; others build it again, named after the instance.
    assert_eq!((child("root", "plain"), child("root", "ten")), ("rod".into(), "rod".into()));
    assert_eq!((child("root", "c"), child("root", "d"), child("root/a", "r")), ("root/c".into(), "root/c".into(), "root/a/r".into()));
    let b = d.body("root/c").unwrap();
    assert_eq!((b.of.as_deref(), b.params.clone()), (Some("rod"), vec![("h".to_string(), "30".to_string()), ("w".to_string(), "2".to_string())]));
    // The moving rod reads t (through the slider's argument); the others do not.
    let reads_t = |n: &str| d.g.depends_on(d.body(n).unwrap().shape, |v| matches!(v, ink::graph::Var::T));
    assert!(reads_t("root/a/r") && !reads_t("rod") && !reads_t("root/c"));
    // Its shape and its frame follow the slider: 20 long at t = 0, 30 at t = 1.
    let inside = |t: f64, p: [f64; 3]| eval_point(&all, &Ctx::new(t, vec![]), p) < 0.0;
    assert!(!inside(0.0, [0.0, 0.0, 25.0]) && inside(1.0, [0.0, 0.0, 25.0]) && inside(0.0, [30.0, 0.0, 25.0]));
    let Some(Value::Vec(v)) = d.def("tip") else { panic!() };
    let z = |t: f64| if let Value::Field(f, _) = v[2] { eval_point(&Tape::new(&d.g, f), &Ctx::new(t, vec![]), [0.0; 3]) } else { panic!() };
    assert!((z(0.0) - 20.0).abs() < 1e-12 && (z(1.0) - 30.0).abs() < 1e-12);
    for (bad, msg) in [
        ("(assembly e (instance r (rod 1 2 3)))", "'rod' takes 2 argument(s) (h w), got 3"),
        ("(assembly e (instance r (rod x)))", "not fields of x y z"),
        ("(body (bad (h 1) h) (shape (box 1 1 1)))", "an argument is (name default)"),
        ("(body (bad (a 1) (a 2)) (shape (box a a a)))", "duplicate argument 'a'"),
    ] {
        let e = ink::load_str(&format!("(import \"std\")\n{src}\n{bad}"), "test", Path::new(ROOT)).err().unwrap_or_else(|| panic!("{bad} loads"));
        assert!(e.message.contains(msg), "{bad}: {e:?}");
    }
}

/// Over a motion, a contact is proven over each step's moments at once
/// (bounds over the step see the parts move), or at the step's middle to
/// the tolerance less how far one part moves relative to the other. A
/// block sliding 10 mm along a wall it touches is proven at rest, no
/// pieces: the wall without its end planes reads nothing along the slide
/// (motion::Relaxed; not the 16 starting pieces, its face seen still by
/// bounds over each); a block turning on a still table it rests on, 8 mm
/// out, by its envelope (the cylinder it turns in, resting on the table)
/// with no pieces, not 16 steps over time. Resting on a turntable that
/// turns it, one instance's parts,
/// the pair is asked with t ranged: one search, no pieces.
#[test]
fn a_contact_over_a_motion_steps_by_how_its_parts_move_relative_to_each_other() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    use ink::view::ModelView;
    let slide = "(import \"std\")
(t-span 0 1)
(body block (shape (box 2 2 2)))
(body wall (shape (move (box 40 2 4) 5 2 0)))
(assembly root (instance block block (at (move 0 0 0)) (slide x (* 10 t))) (instance wall wall))
(assert along (disjoint (shape-of root/block) (shape-of root/wall)))";
    let table = "(import \"std\")
(t-span 0 1)
(body table (shape (cylinder-z 10 2)))
(body block (shape (move (box 2 2 2) 8 0 2)))";
    let turn = format!(
        "{table}
(assembly root (instance t table) (instance b block (at (move 0 0 0)) (turn z (* 90 t))))
(assert seated (disjoint (shape-of root/t) (shape-of root/b)))"
    );
    let together = format!(
        "{table}
(assembly top (instance table table) (instance block block))
(assembly root (instance t top (at (move 0 0 0)) (turn z (* 90 t))))
(assert seated (disjoint (shape-of root/t/table) (shape-of root/t/block)))"
    );
    let steps = |src: &str, tol: Option<f64>| {
        let doc = ink::load_str(src, "motion", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
        let ctx = doc.ctx(&[], 0.0).unwrap();
        let mut v = ModelView::new(doc, ctx);
        v.set_motion(v.doc.t_span).unwrap();
        let r = run_asserts(&mut v, None, false, ClaimOptions { accuracy: tol, ..ClaimOptions::default() }).unwrap();
        let c = &r.asserts[0].check;
        assert!(c.status == Status::Holds && c.accuracy == Some(tol.unwrap_or(0.05)), "{c:?}");
        c.motion.as_ref().unwrap().motion.steps
    };
    assert_eq!((steps(slide, Some(0.5)), steps(slide, None)), (0, 0));
    assert_eq!((steps(&turn, Some(0.5)), steps(&turn, None)), (0, 0));
    assert_eq!((steps(&together, Some(0.5)), steps(&together, None)), (0, 0));
}

/// tests/motion: tiny moving models, each checked over its `(t-span ...)`
/// with t ranged, every claim as its
/// `; expect-status:` lines say (else it holds), over their declared
/// tolerances. Each catches a wrong answer over the motion (a body's
/// radius, a parent's turn, a slide's or a value's rate, a ranged rate, a
/// nonlinear drive, a jump) or a measured input (solved, sampled) taken
/// for a proof.
#[test]
fn motion_models_hold_or_fail_as_declared() {
    use ink::claims::{ClaimOptions, run_asserts};
    use ink::view::ModelView;
    let dir = Path::new(ROOT).join("tests/motion");
    let mut files: Vec<_> = std::fs::read_dir(&dir).unwrap().map(|e| e.unwrap().path()).filter(|p| p.extension().is_some_and(|x| x == "ink")).collect();
    files.sort();
    assert!(files.len() >= 5);
    let load = |f: &Path| {
        let mut v = ModelView::load_in(ink::registry(), f, &[], 0.0)?;
        // Over the declared tolerances, as ink assert ranges them.
        let spans = v.declared_spans().unwrap();
        v.set_param_spans(&spans).unwrap();
        v.set_motion(v.doc.t_span).unwrap();
        Ok::<_, ink::InkError>(v)
    };
    // The models side by side: one after another, they are this binary's longest test.
    use rayon::prelude::*;
    files.par_iter().for_each(|f| {
        let text = std::fs::read_to_string(f).unwrap();
        let mut v = match load(f) {
            Ok(v) => v,
            // A model using a package this build lacks (a solve, without design).
            Err(e) if e.message.ends_with("not in this build") => return,
            Err(e) => panic!("{f:?}: {e:?}"),
        };
        assert!(v.motion.is_some(), "{f:?} declares its motion");
        let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
        for a in &r.asserts {
            let want = text.lines().find_map(|l| l.strip_prefix(&format!("; expect-status: {} ", a.name))).map(str::trim);
            // As the model rule's verdict: the status, and `estimate` after
            // one that holds or fails with no tolerance.
            let st = a.check.status;
            let estimate = st != ink::claims::Status::Unknown && a.check.accuracy.is_none();
            let got = if estimate { format!("{} estimate", st.as_str()) } else { st.as_str().to_string() };
            match want {
                Some(w) => assert_eq!(got, w, "{f:?} {}: {:?}", a.name, a.check),
                None => assert!(a.check.holds_guaranteed(), "{f:?} {}: {:?}", a.name, a.check),
            }
            // A ball 1e12 out is outside ink's range: not checked.
            if matches!(a.check.evidence, ink::claims::Evidence::Outside(_)) {
                continue;
            }
            let m = a.check.motion.as_ref().unwrap_or_else(|| panic!("{f:?} {} was stepped", a.name));
            // An estimate never comes from the proof path, but a proof
            // about inexact inputs (samples).
            if estimate {
                assert!(!m.bound || a.check.inexact.any(), "{f:?} {}: {:?}", a.name, a.check);
            }
        }
    });
    // The arm's pieces find the post where the arm starts and ends.
    let mut v = load(&dir.join("spin-arm.ink")).unwrap();
    let r = run_asserts(&mut v, Some("sweeps-post"), false, ClaimOptions::default()).unwrap();
    let t = r.asserts[0].check.motion.as_ref().unwrap().motion.t.unwrap();
    // Within 0.02 of an end (NaN fails, as it must).
    assert!((t - 0.5).abs() > 0.48, "{t}");
    // Two lines turned by one drive: the post in reach is hit at a moment
    // the arm reaches it (90 degrees of a turn).
    let mut v = load(&dir.join("shared-drive.ink")).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let (clear, hits) = (&r.asserts[0].check, &r.asserts[1].check);
    assert!(clear.status == ink::claims::Status::Holds, "{clear:?}");
    let t = hits.motion.as_ref().unwrap().motion.t.unwrap();
    assert!(hits.status == ink::claims::Status::Fails && (0.0..=1.0).contains(&t) && ((3600.0 * t) % 360.0 - 90.0).abs() < 3.0, "{hits:?}");
    // Lines at 1/3 and 2/3 of a drive of 1080 degrees a second: the arms
    // that never meet hold, and those that meet only near 700 degrees of
    // each 1080 fail at a moment there, in some second.
    let mut v = load(&dir.join("train-rates.ink")).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let (meets, misses) = (&r.asserts[0].check, &r.asserts[1].check);
    let t = meets.motion.as_ref().unwrap().motion.t.unwrap();
    assert!(meets.status == ink::claims::Status::Fails && (0.638..0.712).contains(&(t % 1.0)), "{meets:?}");
    assert!(misses.status == ink::claims::Status::Holds, "{misses:?}");
    // A speed pulse past its limit fails at a moment inside it.
    let mut v = load(&dir.join("speed-step.ink")).unwrap();
    let r = run_asserts(&mut v, Some("speed-limits:root/pulse"), false, ClaimOptions::default()).unwrap();
    let t = r.asserts[0].check.motion.as_ref().unwrap().motion.t.unwrap();
    assert!((0.53..=0.54).contains(&t), "{t}");
    // A drive's jump is unbounded speed: its speed limit is unknown at the
    // jump, never held.
    let jump = |file: &str, name: &str, at: f64| {
        let mut v = load(&dir.join(file)).unwrap();
        let r = run_asserts(&mut v, Some(name), false, ClaimOptions::default()).unwrap();
        let c = r.asserts[0].check.clone();
        let t = c.motion.as_ref().unwrap().motion.t.unwrap();
        assert!((t - at).abs() < 1e-3 && c.status == ink::claims::Status::Unknown, "{c:?}");
    };
    jump("speed-step.ink", "speed-limits:root/jump", 0.5);
    jump("turn-wraps.ink", "speed-limits:root/flipped", 0.5);
}

#[test]
fn std_blends_are_constructions_of_their_own() {
    // A can's wall names no point of its rim's fillet, so the
    // wall's face is the cylinder alone.
    let src = "(def can-shell (max (tag can-wall (circle-2d 6)) (tag can-flats (- (abs y) 5))))
               (def can (fillet-intersect can-shell (- (abs z) 10) 1.5 90))
               (def ch (chamfer-union (- z) (- x) 2))
               (def sm (smooth-union (- z) (- x) 2))
               (def balls (smooth-union (tag ball (sphere 1)) (tag other (move (sphere 1) 10 0 0)) 1))";
    let (d, t) = field(src, "can");
    let names = |t: &Tape, p: [f64; 3]| -> Vec<String> { ink::eval::provenance(t, &CTX, p).iter().map(|&i| d.g.tag_names[i as usize].clone()).collect() };
    assert!(names(&t, [6.0, 0.0, 0.0]).contains(&"can-wall".to_string()));
    // 45 degrees round the fillet's arc, centred at radius 4.5, height 8.5.
    let (s, c) = 45f64.to_radians().sin_cos();
    let p = [4.5 + 1.5 * c, 0.0, 8.5 + 1.5 * s];
    assert!(eval_point(&t, &CTX, p).abs() < 1e-9);
    let n = names(&t, p);
    assert!(n.contains(&"fillet".to_string()) && !n.contains(&"can-wall".to_string()), "{n:?}");
    // The bevel and the smooth blend likewise, where the two sides are equal.
    for (def, name, p) in [("ch", "chamfer", [-1.0, 0.0, -1.0]), ("sm", "blend", [-0.5, 0.0, -0.5])] {
        let (_d, t) = field(src, def);
        assert!(eval_point(&t, &CTX, p).abs() < 1e-12, "{def}");
        assert!(names(&t, p).contains(&name.to_string()), "{def}: {:?}", names(&t, p));
    }
    // Outside its band (the balls 8 apart, blended over 1) the blend
    // changes nothing: a ball's face is the ball's, whichever operand of
    // the min the blend is.
    let (_d, t) = field(src, "balls");
    for (p, ball) in [([-1.0, 0.0, 0.0], "ball"), ([0.0, 1.0, 0.0], "ball"), ([11.0, 0.0, 0.0], "other")] {
        let n = names(&t, p);
        assert!(!n.contains(&"blend".to_string()) && n.contains(&ball.to_string()), "{p:?}: {n:?}");
    }
}

#[test]
fn a_name_tied_with_what_it_names_forms_the_face() {
    // A flat screw head: its top, the plane z named head-top, ties in the
    // max with the slab's own z. The named side is z itself, by a name and
    // a choice, not a construction on it: the survey's walk sees head-top
    // form the top (a construction computing from z would not).
    let src = "(def head (max (tag countersink (- (len2 x y) (+ 3 z))) (tag head-top z) (slab-z -2 0)))\n\
               (def top (max (tag top-face z) z))";
    let (d, t) = field(src, "head");
    let named = |n: Vec<u32>| n.iter().map(|&i| d.g.tag_names[i as usize].clone()).collect::<Vec<_>>();
    for p in [[0.5, 0.0, 0.0], [-1.0, 1.0, 0.0]] {
        assert!(eval_point(&t, &CTX, p).abs() < 1e-12);
        let decided = named(ink::eval::decided_names(&t, &CTX, p, 0.01));
        assert!(decided.contains(&"head-top".to_string()), "{p:?}: {decided:?}");
    }
    // Tied with what it names directly, the name is the point's.
    let (d2, t2) = field(src, "top");
    let shown: Vec<String> = ink::eval::provenance(&t2, &CTX, [0.3, 0.0, 0.0]).iter().map(|&i| d2.g.tag_names[i as usize].clone()).collect();
    assert!(shown.contains(&"top-face".to_string()), "{shown:?}");
}

#[test]
fn degenerate_forms_are_errors_not_panics() {
    // Operations must not fail: a form or claim missing its parts, or with
    // an empty list where a part goes, is a load error.
    let claims = ["empty", "disjoint", "subset", "within", "clearance", "range", "volume", "min-thickness", "cavities", "value", "overlaps", "integral", "mean", "metric"];
    let mut srcs: Vec<String> = ["()", "(())", "(assert)", "(assert ())", "(assert x)", "(assert x ())", "(assert x (()))"].map(String::from).into();
    for h in ink::lang::CORE_FORMS {
        srcs.extend([format!("({h})"), format!("({h} ())"), format!("({h} () ())"), format!("({h} x ())"), format!("({h} x (()))")]);
    }
    for h in claims {
        srcs.extend([format!("(assert x ({h}))"), format!("(assert x ({h} ()))"), format!("(assert x ({h} () () () () ()))")]);
    }
    for s in &srcs {
        let src = format!("(import \"std\") (body b (shape (sphere 1))) {s}");
        let loaded = std::panic::catch_unwind(|| ink::load_str(&src, "t", Path::new(ROOT)).is_err());
        assert_eq!(loaded.ok(), Some(true), "{s}: expected a load error");
    }
}

#[test]
fn a_circulating_belts_wraps_are_undecided_only_in_their_tooth_band() {
    // A belt whose shift reads t, bounded over a span in which it goes
    // round many times (as a viewer bounds a shape that moves): each wrap's
    // tooth section is turned by its angle, so a cell inside the root
    // circle, or between the land and the back, is decided at any turn.
    // Turned about the origin instead (rotate-z), a cell's box swings
    // toward the centre, and most of the wrap inside the back is maybe
    // surface.
    use ink::interval::I;
    let span = Ctx { t_span: Some((0.0, 10.0)), ..CTX };
    // Two radii exercise the arc's local coordinates at different scales.
    let belt = Path::new(ROOT).join("tests/fixtures/belt.ink");
    for (radius, h) in [(6.0, 0.1), (14.0, 0.2)] {
        let src = format!("(import \"{}\" as b)\n(def belt (b.belt {radius} 60 (* 1000 t)))", belt.display());
        let (_d, t) = field(&src, "belt");
        let (root, land, back) = (radius - 1.25, radius, radius + 0.6);
        for deg in (100..=260).step_by(10) {
            let (s, c) = (deg as f64).to_radians().sin_cos();
            for (r, outside) in [(root - 1.0, true), (0.5 * (land + back), false)] {
                let b = [r * c, r * s, 0.0].map(|v| I::new(v - 0.5 * h, v + 0.5 * h));
                let v = ink::eval::eval_interval(&t, &span, b, None);
                assert!(if outside { v.lo > 0.0 } else { v.hi < 0.0 }, "{src}: the cell at r {r}, {deg} degrees: {v:?}");
            }
        }
    }
}

/// A shape claim's tolerance is a length, whatever the field's scale (a
/// field-value tolerance would not be: 1000 a has the solid of a, and its
/// values a thousand times a's). Two 10 mm cubes whose faces overlap by
/// 0.03, or by 0.08 (no point deeper than 0.04, though the finest cells lie
/// wholly inside the overlap: points beyond them settle them), are proven
/// disjoint to 0.05 mm, with either one's field scaled by 1000 or by
/// 1/1000; by 0.2 they are disproven, whatever the scale, the witness's
/// cube of half-side 0.05 inside both.
#[test]
fn a_tolerance_is_a_length_whatever_the_fields_scale() {
    use ink::claims::Status;
    let mut src = String::from("(def a (box 10 10 10))\n");
    for (name, gap) in [("shallow", 9.97), ("shallow-thick", 9.92), ("deep", 9.8)] {
        for (s, k) in [("1", 1.0), ("big", 1000.0), ("small", 0.001)] {
            src += &format!("(assert {name}-{s} (disjoint (* {k} a) (move (box 10 10 10) {gap} 0 0)))\n");
        }
    }
    for (name, c) in run_claims(&src) {
        let want = if name.starts_with("shallow") { Status::Holds } else { Status::Fails };
        assert!(c.status == want && c.accuracy == Some(0.05), "{name}: {c:?}");
        if let ink::claims::Evidence::Empty(ink::ranged::Proof { verdict: ink::query::Verdict::Disproven { witness, .. }, .. }) = c.evidence {
            // 0.05 inside both along x: the overlap is 4.8 .. 5.0.
            assert!((4.85..=4.95).contains(&witness[0]), "{name}: {witness:?}");
        }
    }
}

/// Parts that sit against each other: `touches` holds where they overlap
/// no deeper than the tolerance and come within it of each other (flush,
/// 0.03 into each other, 0.03 apart), and fails where they float (0.2
/// apart: not seated) or overlap deeper (0.2 into each other).
#[test]
fn touching_parts_are_seated_and_not_clashing() {
    use ink::claims::Status;
    let mut src = String::from("(def a (box 10 10 10))\n");
    for (name, x) in [("flush", 10.0), ("into", 9.97), ("near", 10.03), ("floats", 10.2), ("clashes", 9.8)] {
        src += &format!("(assert {name} (touches a (move (box 10 10 10) {x} 0 0)))\n");
    }
    let got: Vec<(String, Status)> = run_claims(&src).into_iter().map(|(n, c)| (n, c.status)).collect();
    let want = [("flush", Status::Holds), ("into", Status::Holds), ("near", Status::Holds), ("floats", Status::Fails), ("clashes", Status::Fails)];
    assert_eq!(got, want.map(|(n, s)| (n.to_string(), s)).to_vec());
}

#[test]
fn touching_over_params_fails_only_apart_at_every_value() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    // A contact found at one value says nothing of the others: seated
    // throughout, or floating at some values, is unknown; apart throughout
    // fails.
    let src = "(import \"std\")
(param near 0.01 (tol 0 0.02))
(param far 0.01 (tol 0 1))
(param off 1 (tol 0 1))
(def a (box 1 1 1))
(assert seated (touches a (move (box 1 1 1) (+ 1 near) 0 0)))
(assert floats (touches a (move (box 1 1 1) (+ 1 far) 0 0)))
(assert apart (touches a (move (box 1 1 1) (+ 1 off) 0 0)))";
    let doc = ink::load_str(src, "touch", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let spans = v.declared_spans().unwrap();
    v.set_param_spans(&spans).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let got: Vec<Status> = r.asserts.iter().map(|a| a.check.status).collect();
    assert_eq!(got, [Status::Unknown, Status::Unknown, Status::Fails], "{r:?}");
}

/// A point (std's probe) is inside or outside a shape by the field's sign
/// there, exactly: no tolerance (a tiny cube held to a tolerance larger
/// than itself could not fail), and none may be given.
#[test]
fn a_probe_is_a_point_decided_exactly() {
    use ink::claims::Status;
    let c = run_claims(
        "(def a (box 10 10 10))
         (assert in (subset (probe 4.999 0 0) a))
         (assert out (disjoint a (probe 5.001 0 0)))
         (assert not-in (subset (probe 5.001 0 0) a))
         (assert not-out (disjoint (probe 4.999 0 0) a))",
    );
    let got: Vec<Status> = c.iter().map(|(_, c)| c.status).collect();
    assert_eq!(got, [Status::Holds, Status::Holds, Status::Fails, Status::Fails]);
    // Exact: a tolerance of 0.
    assert!(c.iter().all(|(_, c)| c.accuracy == Some(0.0)), "{c:?}");
}

/// A point moves with what moves it: place puts it at a frame (R p + T,
/// forward, as std's in-frame-point at a path frame). Arithmetic on it is
/// an error: a point is no field.
#[test]
fn a_probe_moves_as_a_point() {
    use ink::claims::Status;
    let c = run_claims(
        "(def at-1 (move (sphere 0.01) 0 1 0))
         (assert turned (subset (place (probe 1 0 0) (at (rot-z 90))) at-1))
         (assert moved (subset (place (probe 1 2 3) (at (move -1 -1 -3))) at-1))
         (assert placed (subset (place (probe 1 0 0) (at (rot-z 90) (move -5 0 0))) (move (sphere 0.01) -5 1 0)))
         (assert in-axes-of-vecs (subset (place (probe 1 2 3) (axes (origin (vec 5 0 0)) (x 0 1 0) (z (vec 0 0 1)))) (move (sphere 0.01) 3 1 3)))",
    );
    assert!(c.len() == 4 && c.iter().all(|(_, c)| c.status == Status::Holds), "{c:?}");
    for bad in ["(subset (min (probe 0 0 0) (probe 1 0 0)) (sphere 2))", "(subset (probe x 0 0) (sphere 2))"] {
        let e = ink::load_str(&format!("(import \"std\") (assert c {bad})"), "t", Path::new(ROOT)).err().expect(bad);
        assert!(e.message.contains("point"), "{bad}: {e:?}");
    }
}

/// A warp says where a shape is read from, not where a point goes (that
/// is its inverse), so a warp of a point is an error. Each warp below,
/// taken as moving the point, to its inverse by cofactors over the
/// determinant, would prove each claim at a point no preimage is at: a
/// warp not affine though its derivatives are numbers (floor's); a time
/// map reading x, which makes the inner warp's 1 + t^2 a field of x (x (1 +
/// x^2) = 1 has no root at 1); a matrix singular in real numbers (row 3 =
/// row 1 + row 2) whose rounded determinant is 32; a scale of 1e200 whose
/// inverse overflows to the origin; t x over t in [0, 1], singular at t =
/// 0 (there an estimate that holds, at the steps).
#[test]
fn a_warp_does_not_move_a_point() {
    let loaded: Vec<&str> = [
        "(assert c (subset (move (probe 1 2 3) -1 -1 -3) (sphere 2)))",
        "(assert c (subset (rotate-z (probe 1 0 0) 90) (sphere 2)))",
        "(assert floor (subset (warp (probe 2.5 0 0) (+ x (floor x)) y z) (sphere 2)))",
        "(def p (warp (warp (probe 1 0 0) (* (+ 1 (square t)) x) y z) x y z x))
         (assert false-location (subset p (move (sphere 0.05) 1 0 0)))",
        "(def u (+ (* 453597 x) (* 688211 y) (* 450010 z)))
         (def v (+ (* 835775 x) (* 357072 y) (* 271732 z)))
         (def p (warp (probe 0 0 1) u v (+ u v)))
         (assert impossible (subset p (sphere 1e12)))",
        "(param s 1e200)
         (def p (warp (probe s s s) (* s x) (* s y) (* s z)))
         (assert wrong-origin (subset p (sphere 0.1)))",
        "(t-span 0 1)
         (def p (warp (probe 0 0 0) (* t x) y z))
         (assert no-preimage-at-zero (subset p (sphere 0.1)))",
    ]
    .into_iter()
    .filter(|src| ink::load_str(&format!("(import \"std\") {src}"), "t", Path::new(ROOT)).map_err(|e| assert!(e.message.contains("point"), "{src}: {e:?}")).is_ok())
    .collect();
    assert!(loaded.is_empty(), "a point moved by a warp: {loaded:#?}");
}

/// A point's sign is its bounds there, rounded outward: at x = 1e16 they
/// straddle 0, and neither side is proven. A
/// point placed past the largest number is at infinity, no point: the field
/// there (a sphere's, infinite) would prove it both outside and inside the
/// sphere's complement; so too with another coordinate ranged, where the
/// field reads the coordinates as expressions and the overflow is not
/// looked for at a placeholder origin.
#[test]
fn a_point_far_out_is_never_both_inside_and_outside() {
    use ink::claims::{ClaimOptions, Status, run_asserts};
    let src = "(param u 1e16)
         (def f (- (+ x 1) (+ u 2)))
         (assert outside (disjoint (probe u 0 0) f))
         (assert inside (subset (probe u 0 0) f))
         (def far (place (probe 1e308 0 0) (at (move 1e308 0 0))))
         (assert overflowed-outside (disjoint (sphere 1) far))
         (assert overflowed-inside (subset far (neg (sphere 1))))";
    let c = run_claims(src);
    assert!(c.len() == 4 && c.iter().all(|(_, c)| c.status == Status::Unknown), "{c:?}");
    let mut v = ranged_view(
        "(param s 0)
         (def far (place (probe 1e308 s 0) (at (move 1e308 0 0))))
         (assert overflowed-outside (disjoint (sphere 1) far))
         (assert overflowed-inside (subset far (neg (sphere 1))))",
        &[("s", 0.0, 1.0)],
    );
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    assert!(r.asserts.len() == 2 && r.asserts.iter().all(|a| a.check.status == Status::Unknown), "{:?}", r.asserts.iter().map(|a| &a.check).collect::<Vec<_>>());
}

/// A point on a branch cut of its field's formula: an ISO thread's angle
/// cuts at its roots, where the angle's bounds at the point, rounded
/// outward, reach either side of atan2's cut. Bounded as a choice of the
/// two sides (eval::forked), each near pi, the point's own bounds are one
/// number's (a probe of metal at the internal thread's minor diameter; a
/// DIN 1587 nut's).
#[test]
fn a_probe_on_a_thread_root_is_decided() {
    use ink::claims::{ClaimOptions, Evidence, Status, run_asserts};
    let f = Path::new(ROOT).join("tests/fixtures/thread-root.ink");
    let mut v = ink::view::ModelView::load_in(ink::registry(), &f, &[], 0.0).unwrap();
    for name in ["thread-minor-metal", "din-root"] {
        let r = run_asserts(&mut v, Some(name), false, ClaimOptions::default()).unwrap();
        let c = &r.asserts[0].check;
        let Evidence::Point { bounds, .. } = c.evidence else { panic!("{c:?}") };
        assert!(c.status == Status::Holds && bounds.1 < 0.0, "{name}: {c:?}");
    }
}

/// A point on a tooth's centre line, where std's involute gear folds its
/// angle (a mod of atan2 that wraps there): on the +x axis the angle is
/// exactly 0 (SUS1-15's back-frame-in: placed by a half turn, exact since
/// quarter turns are); a tooth over, turned by 24 degrees, its bounds reach
/// either side of the wrap, a choice of two sides (eval::forked), not the
/// whole tooth they span.
#[test]
fn a_probe_on_a_gear_tooth_centre_line_is_decided() {
    use ink::claims::{ClaimOptions, Evidence, Status, run_asserts};
    let src = "(import \"std\")
               (def g (max (involute-gear-2d 15 1 20 0.03) (neg z) (- z 10)))
               (assert root (subset (probe 6.25 0 5) g))
               (assert turned (subset (place (probe 6.25 0 -0.06) (at (rot-x 180))) g))
               (assert tooth-1 (subset (place (probe 6.25 0 5) (at (rot-z 24))) g))";
    let doc = ink::load_str(src, "claims", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    for a in run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts {
        let c = &a.check;
        let Evidence::Point { bounds, .. } = c.evidence else { panic!("{c:?}") };
        assert!(c.status == Status::Holds && bounds.1 < 0.0, "{}: {c:?}", a.name);
    }
}

/// An overlap whose intersection's proven box is thinner than twice the
/// accuracy cannot be shown at it (no point of it is deeper than the
/// accuracy), and says so: a probe box 0.02 across against the model's
/// 0.05. At an accuracy under half its side it holds, as does a fit 0.02
/// deep (unknown at 0.05: the boxes of 10 mm shapes are proven to cells
/// coarser than the fit, so its hint is the search's); a deep overlap
/// holds and says nothing of it.
#[test]
fn an_overlap_thinner_than_twice_the_accuracy_says_so() {
    use ink::claims::{ClaimOptions, Evidence, Status, run_asserts};
    let src = "(import \"std\")\n(body f (shape (box 10 10 10)))\n\
               (assert probe (overlaps (shape-of f) (box 0.02 0.02 0.05)))\n\
               (assert fit (overlaps (shape-of f) (move (box 2 2 2) 5.98 0 0)))\n\
               (assert probe-fine (overlaps (shape-of f) (box 0.02 0.02 0.05)) (accuracy 0.005))\n\
               (assert fit-fine (overlaps (shape-of f) (move (box 2 2 2) 5.98 0 0)) (accuracy 0.005))\n\
               (assert deep (overlaps (shape-of f) (move (box 2 2 2) 4.5 0 0)))\n";
    let doc = ink::load_str(src, "thin-overlap", Path::new(ROOT)).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut v = ink::view::ModelView::new(doc, ctx);
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    let thin = |i: usize| match &r.asserts[i].check.evidence {
        Evidence::Overlaps(_, t) => *t,
        e => panic!("{e:?}"),
    };
    let c = &r.asserts[0].check;
    assert!(c.status == Status::Unknown && thin(0).is_some_and(|s| (0.02..0.1).contains(&s)), "{c:?}");
    assert_eq!(r.asserts[1].check.status, Status::Unknown);
    for i in [2, 3, 4] {
        let c = &r.asserts[i].check;
        assert!(c.status == Status::Holds && thin(i).is_none(), "{}: {c:?}", r.asserts[i].name);
    }
}
