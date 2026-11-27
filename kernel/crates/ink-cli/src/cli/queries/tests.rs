use super::budget::{READ_WORK, TRIANGLE_WORK, VIEW_TRIANGLES, cost_warnings, ladder, view_floor, work};
use super::{SEEN, cmd_distance, make_meshes, paint_tapes, plan_meshes, tape_of};
use crate::{Args, CliView, J, ModelView, Path, Tape, V3, json, mesh};
use rayon::prelude::*;
use std::sync::Arc;

/// How far, as a factor either way, the instructions a view mesh
/// retired when last counted may be from TRIANGLE_WORK times its
/// `work`.
const WORK_TOLERANCE: f64 = 1.4;

/// Shapes the budget's ratios are held to, each with its res and the
/// machine instructions making its view mesh retired, as
/// `calibrate_view_work` printed them (one thread of a Zen 5, the test
/// build): a block with a hole (most of it its triangles), a bar swept
/// along a helix (its field's samples and proofs), a gyroid sheet (its
/// many triangles and transcendental field), a threaded rod and a ball
/// painted across (its paint).
const VIEW_CORPUS: [(&str, &str, &str, f64, f64); 5] = [
    ("block", "(subtract (box 40 30 20) (cylinder-z 5 30))", "", 0.5, 4.3324e9),
    ("sweep", "(sweep (rounded-rect-2d 6 3 0.8) (helix-points 15 10 1 30) (vec 0 0 1))", "", 0.5, 3.4740e10),
    ("gyroid", "(lattice-fill (box 20 20 20) (gyroid 10) 1)", "", 0.5, 1.3449e10),
    ("thread", "(max (thread-iso 8 1.25) (box 10 10 16))", "", 0.1, 2.0380e10),
    ("painted", "(sphere 15)", "(attr color (vec (clamp (/ (+ x 15) 30) 0 1) 0.5 0.5))", 0.3, 9.2696e9),
];

/// A view mesh of a corpus shape, made (not from the cache), and the
/// instructions making it retired, counted where `count` and the
/// system allow.
fn corpus_mesh(shape: &str, attr: &str, h: f64, count: bool) -> (ink::seen::Viewed, Option<counted::Cost>) {
    let doc = ink::load_str(&format!("(import \"std\")\n(body b (shape {shape}) {attr})"), "work", Path::new(".")).unwrap();
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut s = ModelView::new(doc, ctx);
    let it = s.doc.items("b").unwrap().remove(0);
    let (lo, hi, _) = s.item_bounds(&it).unwrap();
    let (tape, paint) = (tape_of(&mut s, &it, lo, hi, h), paint_tapes(&s, &it.attrs));
    let how = ink::seen::How {
        res: h,
        limit: mesh::Limit { bricks: usize::MAX, search: mesh::Limit::samples_for(1e6) },
        tol: None,
        tol_per_res: SEEN.tol_per_res,
        bend: None,
        export: false,
    };
    let make = || ink::seen::make(&tape, &paint, &s.ctx, lo, hi, how);
    match count {
        true => counted::instructions(make).map_or_else(|| (make(), None), |(v, c)| (v, Some(c))),
        false => (make(), None),
    }
}

/// The budget's ratios (`work`, from what the kernel counts of a mesh)
/// hold against the instructions view meshes of VIEW_CORPUS retired when
/// last counted, each within WORK_TOLERANCE. Counted by the kernel, so
/// it checks the same on any machine: a mesher or evaluator that comes
/// to cost otherwise than counted is caught only by counting again
/// (scripts/view-budget-calibrate.sh).
#[test]
fn the_view_budget_charges_what_meshes_count() {
    let off: Vec<(&str, f64)> =
        VIEW_CORPUS.par_iter().map(|&(name, shape, attr, h, counted)| (name, counted / (TRIANGLE_WORK * work(&corpus_mesh(shape, attr, h, false).0)))).collect();
    assert!(off.iter().all(|o| o.1 < WORK_TOLERANCE && o.1 > 1.0 / WORK_TOLERANCE), "counted over charged: {off:?}");
}

/// A body with no inside (move's arguments swapped: its field is 0
/// over its box) is an empty target: the distance says so at once,
/// without meshing the body first, finer each time to 32,768 cells a
/// side (past 30 s).
#[test]
fn distance_to_a_body_with_no_inside_ends() {
    let src = "(import \"std\")
(body base (shape (box 10 10 2)))
(body flat (shape (move 0 0 5 (box 4 4 1))))
(assembly root (instance base base) (instance flat flat))";
    let doc = ink::load_str(src, "distance", Path::new(".")).unwrap();
    let ctx = doc.ctx(&[], 0.0).unwrap();
    let mut s = ModelView::new(doc, ctx);
    let args = Args::parse(&["distance", "unused.ink", "--a", "root/flat", "--b", "root/base"].map(str::to_string));
    let Err(e) = cmd_distance(&mut s, &args) else { panic!("a distance to no material") };
    assert!(e.message.contains("empty"), "{e:?}");
}

#[test]
fn a_view_floor_coarsens_the_finest_meshes_first_within_the_budget() {
    // (triangles at res 1, own res, times made): two big parts and many
    // small detailed ones, each at its own res sampling about the same.
    let mut parts = vec![(2e5, 0.6, 1.0), (1.8e5, 0.55, 1.0)];
    parts.extend((0..30).map(|k| (2e5 * (0.05 + 0.002 * k as f64).powi(2), 0.05 + 0.002 * k as f64, 1.0)));
    let sampled = |f: f64| parts.iter().map(|&(c, h, _)| c / h.max(f).powi(2)).sum::<f64>();
    assert!(sampled(0.0) > VIEW_TRIANGLES);
    let f = view_floor(&parts);
    assert!(sampled(f) <= VIEW_TRIANGLES && sampled(0.99 * f) > VIEW_TRIANGLES, "the least floor within it: {f}");
    assert!(f > 0.11 && f < 0.55, "small parts coarsened, the big ones not made finer: {f}");
    // On the ladder (rounded up) it stays within the budget.
    assert!(parts.iter().map(|&(c, h, _)| c / ladder(h.max(f)).powi(2)).sum::<f64>() <= VIEW_TRIANGLES);
    assert_eq!(view_floor(&parts[..3]), 0.0, "a few parts keep their own res");
    // A flipbook of 42 meshes of a big part beside many small ones: at
    // f times the cube root of 42 its meshes sample what that many
    // copies of it (once each) sample at f, so the floor is theirs, not
    // that of 42 copies (a user's gantry's belt runs: 0.84, not 13.45).
    let small: Vec<(f64, f64, f64)> = (0..60).map(|k| (2e5 * (0.05 + 0.002 * k as f64).powi(2), 0.05 + 0.002 * k as f64, 1.0)).collect();
    let with = |p: (f64, f64, f64)| [small.clone(), vec![p]].concat();
    let n = 42.0f64;
    let flipbook = view_floor(&with((n * 1e6, 0.1, n)));
    let copies = view_floor(&with((n.cbrt() * 1e6, 0.1, 1.0)));
    assert!((flipbook / copies - 1.0).abs() < 1e-6, "{flipbook} {copies}");
    assert!(flipbook < 0.6 * view_floor(&with((n * 1e6, 0.1, 1.0))), "{flipbook}");
}

#[test]
fn a_mesh_that_alone_brings_the_others_coarser_is_named_on_the_ladder() {
    // Ten small parts (2e6 triangles each at their own res 0.1), and a
    // flipbook of 41 frames: without it the floor is 0.28 (on the
    // ladder 0.297), with it 0.47 (0.5).
    let small = vec![(2e4, 0.1, 1.0); 10];
    let all = [small.clone(), vec![(42.0 * 1e5, 0.1 * 42f64.sqrt(), 42.0)]].concat();
    let warned = |parts: &[(f64, f64, f64)]| {
        let floor = view_floor(parts);
        let boxes: Vec<(V3, V3, f64, Option<f64>)> = parts.iter().map(|p| ([0.0; 3], [100.0, 10.0, 10.0], ladder(p.1.max(floor * p.2.cbrt())), None)).collect();
        let made: Vec<(f64, f64)> = parts.iter().zip(&boxes).map(|(p, b)| (b.2, p.0 / (b.2 * b.2))).collect();
        let names: Vec<String> = (0..parts.len()).map(|k| format!("p{k}")).collect();
        cost_warnings(parts, floor, &boxes, &made, &names.iter().map(String::as_str).collect::<Vec<_>>())
    };
    let w = warned(&all);
    assert!(w.len() == 1 && w[0].starts_with("'p10': its 41 frames at res 1.681792831 over 100.0 mm") && w[0].ends_with("up to res 0.5 (without it, res 0.2973017788)"), "{w:?}");
    assert!(warned(&small).is_empty());
    // The ladder rounds up to quarter powers of two.
    for (h, want) in [(0.25, 0.25), (0.2, 0.21022410381342863), (0.42, 0.42044820762685725), (0.43, 0.5), (1.0, 1.0), (1.19, std::f64::consts::SQRT_2), (0.0568, 0.0625)] {
        assert!((ladder(h) - want).abs() < 1e-12 && ladder(h) >= h * (1.0 - 1e-9), "{h}: {}", ladder(h));
    }
    assert_eq!(ladder(ladder(0.3)), ladder(0.3), "a rung is its own");
    assert_eq!(ladder(ladder(0.3) / 2.0), ladder(0.3) / 2.0, "half a rung is a rung");
}

use crate::cli::counted;

/// What a view mesh counted of its work, and the instructions it took.
fn counted_json(v: &ink::seen::Viewed, c: counted::Cost) -> J {
    let u = &v.unresolved;
    json!({"instructions": c.instructions, "cycles": c.cycles, "seconds": c.seconds, "work": work(v), "uniform_triangles": v.uniform_triangles, "triangles": v.tris.len(),
        "sampled": u.sampled, "read": u.read, "field_work": v.field_work, "paint_work": v.paint_work})
}

/// The budget's ratios, calibrated: on one thread, with the processor
/// counting instructions, VIEW_CORPUS's meshes (printed as its rows),
/// and for each model in INK_WORK_MODELS (paths, ':' between; INK_PATH
/// at the checkout) a line of JSON: the instructions its view meshes
/// (as `ink mesh M.glb` makes them, made again finer included) retire
/// and what the budget charged them (`Budget::spent`), then each first
/// mesh's own with what it counted. Last, TRIANGLE_WORK and READ_WORK
/// fitted to them all (least squares of the ratio of predicted to
/// counted), and each mesh's ratio at the fit. Run with INK_NO_CACHE
/// set (scripts/view-budget-calibrate.sh).
#[test]
#[ignore = "a calibration's, run by hand"]
fn calibrate_view_work() {
    assert!(std::env::var_os("INK_NO_CACHE").is_some(), "meshes from the cache cost nothing");
    let none = "this system does not let a process count its instructions";
    // Each mesh's name, triangles, instructions counted by the kernel
    // (as `work` reads them), and by the processor.
    let mut fit: Vec<(String, f64, f64, f64)> = Vec::new();
    let mut row = |name: String, v: &ink::seen::Viewed, c: counted::Cost| {
        let u = &v.unresolved;
        fit.push((name, v.uniform_triangles as f64, (u.sampled + u.read + v.field_work + v.paint_work) as f64, c.instructions as f64));
    };
    for (name, shape, attr, h, _) in VIEW_CORPUS {
        let (v, c) = counted::on_one_thread(|| corpus_mesh(shape, attr, h, true));
        let c = c.expect(none);
        println!("        ({name:?}, {shape:?}, {attr:?}, {h:?}, {:.4e}),", c.instructions as f64);
        row(name.to_string(), &v, c);
    }
    let models = std::env::var("INK_WORK_MODELS").unwrap_or_default();
    for path in models.split(':').filter(|p| !p.is_empty()) {
        let meshes: Vec<(String, ink::seen::Viewed, counted::Cost)> = counted::on_one_thread(|| {
            let doc = ink::load_file(Path::new(path)).unwrap();
            let ctx = doc.ctx(&[], 0.0).unwrap();
            let mut s = ModelView::new(doc, ctx);
            s.path = Some(path.into());
            let a = Args::parse(&["mesh", path].map(str::to_string));
            let items = s.items(&a).unwrap();
            let plan = plan_meshes(&mut s, &a, &items, &SEEN).unwrap();
            // The surfaces' areas make_meshes finds first, to leave out.
            let fields: Vec<Arc<Tape>> = plan.jobs.iter().map(|(it, _)| Arc::new(Tape::new(&s.doc.g, it.shape))).collect();
            let boxes = plan.boxes.clone();
            let (_, surfaces) = counted::instructions(|| fields.iter().zip(&boxes).map(|(t, b)| ink::seen::surface(t, &s.ctx, b.0, b.1)).sum::<f64>()).expect(none);
            let ((_, budget), all) = counted::instructions(|| make_meshes(&mut s, &a, plan, &SEEN).unwrap()).unwrap();
            let budget = budget.unwrap();
            let instructions = all.instructions - surfaces.instructions;
            println!(
                "{}",
                json!({"model": path, "instructions": instructions, "spent": budget.spent, "per_triangle": instructions as f64 / budget.spent, "meshes": budget.boxes.len()})
            );
            let plan = plan_meshes(&mut s, &a, &items, &SEEN).unwrap();
            let mut meshes = Vec::new();
            for (j, (it, _)) in plan.jobs.iter().enumerate() {
                let (lo, hi, h, bend) = budget.boxes[j];
                let planned = budget.parts[j].0 / (h * h);
                let limit = mesh::Limit { bricks: mesh::Limit::bricks_for(SEEN.budget), search: mesh::Limit::samples_for(planned.min(SEEN.budget)) };
                let how = ink::seen::How { res: h, limit, tol: None, tol_per_res: SEEN.tol_per_res, bend, export: false };
                let (tape, paint) = (tape_of(&mut s, it, lo, hi, h), paint_tapes(&s, &it.attrs));
                let (v, c) = counted::instructions(|| ink::seen::make(&tape, &paint, &s.ctx, lo, hi, how)).unwrap();
                let mut line = counted_json(&v, c);
                line["model"] = json!(path);
                line["mesh"] = json!(it.label);
                println!("{line}");
                meshes.push((format!("{path} {}", it.label), v, c));
            }
            meshes
        });
        for (name, v, c) in meshes {
            row(name, &v, c);
        }
    }
    // Least squares of (a t + b r) / i - 1 over the meshes: a is
    // TRIANGLE_WORK, b READ_WORK.
    let (mut m, mut y) = ([[0.0; 2]; 2], [0.0; 2]);
    for &(_, t, r, i) in &fit {
        let x = [t / i, r / i];
        for j in 0..2 {
            y[j] += x[j];
            for k in 0..2 {
                m[j][k] += x[j] * x[k];
            }
        }
    }
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let (a, b) = ((y[0] * m[1][1] - y[1] * m[0][1]) / det, (m[0][0] * y[1] - m[1][0] * y[0]) / det);
    println!("fitted: TRIANGLE_WORK {a:.4e}, READ_WORK {b:.3} (in use: {TRIANGLE_WORK:.4e}, {READ_WORK})");
    for (name, t, r, i) in &fit {
        println!("counted over predicted at the fit {:.3}, as in use {:.3}: {name}", i / (a * t + b * r), i / (TRIANGLE_WORK * t + READ_WORK * r));
    }
}
