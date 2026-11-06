//! The kernel's answer cache (ink_kernel::cache::memo_region): a claim is answered
//! from the cache exactly when what it asks the kernel is unchanged, and a
//! cached answer is the one computing it gives; a measured mass is cached as a claim is; a
//! followed solve's steps are read by a later command asking the same. The
//! tests point the cache at directories of their own (XDG_CACHE_HOME), one
//! at a time (ENV).

mod common;

use ink::lang::claims::{ClaimOptions, run_asserts};
use ink::lang::view::ModelView;
use std::path::Path;

const MODEL: &str = r#"(import "std")
(t-span 0 1)
(param w 10 (tol -0.1 0.1))
(param h 4)
(def scan (import-stl "scan.stl" (res 0.5)))
(body plate (shape (box w h 2)))
(body pin (shape (move (sphere 1) 20 0 0)))
(body slider (shape (move (sphere 1) (* 5 t) 0 0)))
(assert inside (within (shape-of plate) (-6 -3 -2) (6 3 2)))
(assert apart (clearance (shape-of plate) (shape-of pin) 5))
(assert small (volume (shape-of plate) < 100))
(assert pin-fine (within (shape-of pin) (18 -2 -2) (22 2 2)) (res 0.5))
(assert slides-clear (disjoint (shape-of slider) (shape-of pin)))
(assert scan-inside (within scan (-1 -1 -1) (3 3 3)))
"#;

const BAR: &str = r#"(import "std")
(param depth 4)
(body bar (shape (move (box 20 4 depth) 10 0 0)))
"#;

/// A binary STL of the box lo..hi, outward-facing.
fn box_stl(path: &Path, lo: [f32; 3], hi: [f32; 3]) {
    let v = |i: usize| -> [f32; 3] { std::array::from_fn(|a| if (i >> a) & 1 == 1 { hi[a] } else { lo[a] }) };
    let tris: Vec<[[f32; 3]; 3]> = [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]]
        .iter()
        .flat_map(|q| [[v(q[0]), v(q[1]), v(q[2])], [v(q[0]), v(q[2]), v(q[3])]])
        .collect();
    let mut b = vec![0u8; 80];
    b.extend((tris.len() as u32).to_le_bytes());
    for t in &tris {
        b.extend([0u8; 12]);
        t.iter().flatten().for_each(|c| b.extend(c.to_le_bytes()));
        b.extend([0u8; 2]);
    }
    std::fs::write(path, b).unwrap();
}

/// Each claim's answer, everything but its timing, and whether it was
/// answered from the cache: the model checked as `ink assert --own` does
/// (declared tolerances, the file's motion), with the cache in `cache`.
fn answers(dir: &Path, cache: &Path, src: &str, scan: f32) -> Vec<(String, String, bool)> {
    std::fs::write(dir.join("m.ink"), src).unwrap();
    box_stl(&dir.join("scan.stl"), [0.0; 3], [scan; 3]);
    // SAFETY: under ENV, and no ink work runs meanwhile.
    unsafe { std::env::set_var("XDG_CACHE_HOME", cache) };
    let mut v = ModelView::load_in(ink::registry(), &dir.join("m.ink"), &[], 0.0).unwrap_or_else(|e| panic!("{e:?}"));
    let spans = v.declared_spans().unwrap();
    v.set_param_spans(&spans).unwrap();
    let span = v.doc.t_span;
    v.set_motion(span).unwrap();
    let rep = run_asserts(&mut v, None, true, ClaimOptions::default()).unwrap_or_else(|e| panic!("{e:?}"));
    rep.asserts.iter().map(|a| (a.name.clone(), format!("{:?} {:?} {:?} {:?}", a.check, a.terms, a.res, a.claim), a.cached)).collect()
}

/// An answer without where a disproof's witness lies: the search for one
/// runs in parallel, so computing it twice may find two (`Disproven { .. }`
/// kept whole but for its contents).
fn unwitnessed(answer: &str) -> String {
    let mut out = String::new();
    let mut rest = answer;
    while let Some(i) = rest.find("Disproven {") {
        out.push_str(&rest[..i + "Disproven {".len()]);
        rest = &rest[i + "Disproven {".len()..];
        let mut depth = 1;
        let end = rest.char_indices().find(|&(_, c)| {
            depth += match c {
                '{' => 1,
                '}' => -1,
                _ => 0,
            };
            depth == 0
        });
        let end = end.map_or(rest.len(), |e| e.0);
        out.push_str("..");
        rest = &rest[end..];
    }
    out + rest
}

/// Held by a test while it points the cache somewhere.
static ENV: std::sync::Mutex<()> = std::sync::Mutex::new(());

#[test]
fn a_claim_is_cached_exactly_when_its_question_is_unchanged() {
    let _env = ENV.lock().unwrap_or_else(|e| e.into_inner());
    let root = common::Scratch::new("cache-test");
    let dir = root.join("model");
    std::fs::create_dir_all(&dir).unwrap();
    let shared = root.join("shared");
    let mut fresh_dirs = 0;
    let mut fresh = |src: &str, scan: f32| {
        fresh_dirs += 1;
        answers(&dir, &root.join(format!("fresh-{fresh_dirs}")), src, scan)
    };
    let names = |a: &[(String, String, bool)], cached: bool| -> Vec<String> { a.iter().filter(|x| x.2 == cached).map(|x| x.0.clone()).collect() };

    let base = answers(&dir, &shared, MODEL, 2.0);
    assert_eq!(names(&base, true), Vec::<String>::new(), "a cold cache answers nothing");
    assert!(base.iter().all(|a| a.1.starts_with("Check { status: Holds")), "{base:#?}");
    let again = answers(&dir, &shared, MODEL, 2.0);
    assert_eq!(names(&again, false), Vec::<String>::new(), "unchanged, every claim is answered from the cache");
    let strip = |a: &[(String, String, bool)]| a.iter().map(|x| (x.0.clone(), unwitnessed(&x.1))).collect::<Vec<_>>();
    assert_eq!(strip(&again), strip(&base), "a cached answer is the one computed");

    // (what changed, the source, the scan's size, the claims whose question
    // changed: None where the kernel's answers it reads may be cached
    // though the claim's answer changes, as a claim's bound on an integral)
    let cases: Vec<(&str, String, f32, Option<&[&str]>)> = vec![
        // what only the front end sees: a comment, a name, an unused def
        ("front end", MODEL.replace("(def scan", "; a note\n(def unused 3)\n(def scan").replace("plate", "slab"), 2.0, Some(&[])),
        ("geometry", MODEL.replace("(box w h 2)", "(box w h 2.2)"), 2.0, Some(&["inside", "apart", "small"])),
        ("a claim's box", MODEL.replace("(-6 -3 -2) (6 3 2)", "(-6 -3 -2) (5.02 3 2)"), 2.0, Some(&["inside"])),
        ("a claim's bound", MODEL.replace("< 100", "< 70"), 2.0, None),
        // a ranged param is in every claim's scope (a disproof says where in
        // all of it), so its band, or its value, is in every question
        ("a tolerance", MODEL.replace("(tol -0.1 0.1)", "(tol -0.1 0.3)"), 2.0, Some(&["inside", "apart", "small", "pin-fine", "slides-clear", "scan-inside"])),
        ("a param's value", MODEL.replace("(param h 4)", "(param h 4.5)"), 2.0, Some(&["inside", "apart", "small"])),
        ("a claim's res", MODEL.replace("(22 2 2)) (res 0.5)", "(22 2 2)) (res 0.25)"), 2.0, Some(&["pin-fine"])),
        ("the t-span", MODEL.replace("(t-span 0 1)", "(t-span 0 4)"), 2.0, Some(&["slides-clear"])),
        ("a sampled file", MODEL.to_string(), 4.0, Some(&["scan-inside"])),
    ];
    for (what, src, scan, changed) in &cases {
        let reference = fresh(src, *scan);
        let got = answers(&dir, &shared, src, *scan);
        assert_eq!(strip(&got), strip(&reference), "{what}: the answers computing gives");
        match changed {
            Some(changed) => assert_eq!(names(&got, false), *changed, "{what}: the claims computed again"),
            None => assert!(!names(&got, false).is_empty() || what.contains("bound"), "{what}: something computed again"),
        }
        // A changed question gets another answer here (the cases are chosen
        // so), so a stale one would show.
        if *what != "front end" {
            assert_ne!(strip(&got), strip(&base), "{what}: the answers change");
        }
    }
    let last = answers(&dir, &shared, MODEL, 2.0);
    assert_eq!(names(&last, false), Vec::<String>::new());
    assert_eq!(strip(&last), strip(&base));

    // A mass as `ink mass` measures it is a kernel answer too, cached so
    // that a test run does not measure every catalog part's again.
    let file = dir.join("bar.ink");
    std::fs::write(&file, BAR).unwrap();
    // SAFETY: as in `answers`.
    unsafe { std::env::set_var("XDG_CACHE_HOME", root.join("bar")) };
    let mut v = ModelView::load_in(ink::registry(), &file, &[], 0.0).unwrap_or_else(|e| panic!("{e:?}"));
    let bar = v.doc.body("bar").unwrap().shape;
    let (lo, hi, _) = v.bounds(bar).unwrap();
    let tape = v.clipped(bar);
    let misses = || ink::cache::MISSES.load(std::sync::atomic::Ordering::Relaxed);
    let (before, first) = (misses(), ink::query::measured(&tape, None, &v.ctx, lo, hi, None));
    let (after, again) = (misses(), ink::query::measured(&tape, None, &v.ctx, lo, hi, None));
    assert_eq!((after - before, misses() - after), (1, 0), "measured once, then read");
    assert_eq!(format!("{first:?}"), format!("{again:?}"));
}

/// A solve on t followed to a far moment (docs/IR.md, Following along t):
/// a second command reads the steps the first certified, and its root is
/// the same, bit for bit; a chunk gone from disk is certified again; one
/// whose equation or load-time root (the other branch, from another guess)
/// differs certifies every step itself, as with a cold cache.
#[test]
fn a_followed_solve_reads_the_steps_another_command_certified() {
    let _env = ENV.lock().unwrap_or_else(|e| e.into_inner());
    let root = common::Scratch::new("follow-cache");
    let solve = |guess: f64, gap: f64| format!("(param V 1000) (solve s (unknowns (u {guess})) (== (* (- u (* V t)) (- u (* V t) {gap})) 0))");
    // u at t = 0.25 and the steps certified, the cache in `cache`.
    let follow = |src: &str, cache: &str| -> (f64, usize) {
        let f = root.join("m.ink");
        std::fs::write(&f, src).unwrap();
        // SAFETY: under ENV, and no ink work runs meanwhile.
        unsafe { std::env::set_var("XDG_CACHE_HOME", root.join(cache)) };
        let certified = || ink::lang::follow::CERTIFIED.with(std::cell::Cell::get);
        let before = certified();
        let v = ModelView::load_in(ink::registry(), &f, &[], 0.25).unwrap_or_else(|e| panic!("{e:?}"));
        (v.ctx.params[v.doc.param_index("u").unwrap()], certified() - before)
    };
    let (u, all) = follow(&solve(0.0, 1.0), "shared");
    assert!(u == 250.0 && all > 1000, "u {u} after {all} steps");
    let (again, read) = follow(&solve(0.0, 1.0), "shared");
    assert_eq!(again.to_bits(), u.to_bits());
    assert!(read < ink::lang::follow::CHUNK, "{read} steps certified again of {all}");
    // Its chunks evicted, the far ends kept (the smallest entry).
    let mut kept: Vec<(u64, std::path::PathBuf)> =
        std::fs::read_dir(root.join("shared/ink/follow")).unwrap().map(|e| e.unwrap().path()).map(|p| (p.metadata().unwrap().len(), p)).collect();
    kept.sort();
    kept[1..].iter().for_each(|(_, p)| std::fs::remove_file(p).unwrap());
    let (gone, n) = follow(&solve(0.0, 1.0), "shared");
    assert!(gone.to_bits() == u.to_bits() && n == all, "{gone} after {n} steps");
    for (what, src, value) in [("an equation", solve(0.0, 2.0), 250.0), ("the load-time root", solve(1.0, 1.0), 251.0)] {
        let (cold, n) = follow(&src, what);
        let (got, m) = follow(&src, "shared");
        assert!(got == cold && cold == value && m == n, "{what}: {got} after {m} steps, {cold} after {n} cold");
    }
}

/// A quick box is keyed by the same complete question as a full bound,
/// including the region and ranged scope; a hit preserves every float bit.
#[test]
fn quick_bounds_are_cached_by_the_question_and_round_trip_exactly() {
    use ink::eval::Ctx;
    use ink::tape::Tape;
    use std::sync::Arc;

    let _env = ENV.lock().unwrap_or_else(|e| e.into_inner());
    let root = common::Scratch::new("quick-bounds-cache");
    // SAFETY: under ENV, and no ink work runs meanwhile.
    unsafe { std::env::set_var("XDG_CACHE_HOME", &*root) };
    let doc = ink::load_str("(param r 2) (param unused 7) (body ball (shape (- (+ (square (- x t)) (square y) (square z)) (square r))))", "ball.ink", Path::new(".")).unwrap();
    let tape = Arc::new(Tape::new(&doc.g, doc.body("ball").unwrap().shape));
    let ctx = Ctx::new(0.0, vec![2.0, 7.0]);
    let bits = |v: Option<([f64; 3], [f64; 3], f64)>| v.map(|(lo, hi, res)| lo.into_iter().chain(hi).chain([res]).map(f64::to_bits).collect::<Vec<_>>());
    let ask = |ctx: &Ctx, lo, hi| {
        let before = ink::cache::here().0;
        let value = bits(ink::query::bounds_quick(&tape, ctx, lo, hi));
        (value, ink::cache::here().0 - before)
    };
    let (lo, hi) = ([-10.0; 3], [10.0; 3]);
    let (first, misses) = ask(&ctx, lo, hi);
    assert!(first.is_some() && misses == 1, "the cold box is computed");
    assert_eq!(ask(&ctx, lo, hi), (first.clone(), 0), "the warm box is bit-identical and not recomputed");
    let mut unused = ctx.clone();
    unused.params[1] = 99.0;
    assert_eq!(ask(&unused, lo, hi), (first.clone(), 0), "an unread param is not part of the question");
    for (what, ctx, lo, hi) in [
        ("radius", Ctx::new(0.0, vec![3.0, 7.0]), lo, hi),
        ("time", Ctx::new(1.0, ctx.params.clone()), lo, hi),
        ("time span", Ctx { t_span: Some((0.0, 1.0)), ..ctx.clone() }, lo, hi),
        ("param band", ctx.with_band(0, (2.0, 3.0)), lo, hi),
        ("param span", ctx.with_band(0, (2.0, 3.0)).with_param(0, (2.0, 2.5)), lo, hi),
        ("region", ctx.clone(), [-5.0; 3], [5.0; 3]),
        ("empty region", ctx.clone(), [8.0; 3], [10.0; 3]),
    ] {
        let (value, misses) = ask(&ctx, lo, hi);
        assert_eq!(misses, 1, "a changed {what} is computed");
        assert_eq!(ask(&ctx, lo, hi), (value, 0), "a changed {what} is then cached exactly");
    }
    // The normal and bypass paths compute the same box.
    // SAFETY: under ENV, and no ink work runs meanwhile.
    unsafe { std::env::set_var("INK_NO_CACHE", "1") };
    let uncached = ask(&ctx, lo, hi);
    unsafe { std::env::remove_var("INK_NO_CACHE") };
    assert_eq!(uncached, (first, 1));
}
