//! The kernel's answer cache (ink_kernel::cache::memo_region): a claim, a
//! wall's thickness included, is answered from the cache exactly when what
//! it asks the kernel is unchanged (over a solve on ranged params too),
//! and a cached answer is the one computing it gives; a measured mass is cached as a claim is; a followed
//! solve's steps are read by a later command asking the same. The tests
//! point the cache at directories of their own (XDG_CACHE_HOME), each in a
//! process of its own (`ran_alone`), so they run side by side.

mod common;

use ink::lang::claims::{ClaimOptions, run_asserts};
use ink::lang::view::ModelView;
use std::path::Path;

const MODEL: &str = r#"(import "std")
(t-span 0 1s)
(param w 10 (tol -0.1 0.1))
(param h 4)
(def scan (import-stl "scan.stl" (res 0.5)))
(body plate (shape (box w h 2)))
(body pin (shape (move (sphere 1) 20 0 0)))
(body slider (shape (move (sphere 1) (* 5 (/ t 1s)) 0 0)))
(require inside (within (shape-of plate) (-6 -3 -2) (6 3 2)))
(require apart (clearance (shape-of plate) (shape-of pin) 5))
(require small (volume (shape-of plate) < 100))
(require walls (min-thickness (shape-of pin) > 1))
(require pin-fine (within (shape-of pin) (18 -2 -2) (22 2 2)) (res 0.5))
(require slides-clear (disjoint (shape-of slider) (shape-of pin)))
(require scan-inside (within scan (-1 -1 -1) (3 3 3)))
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
/// answered from the cache: the model checked as `ink require --own` does
/// (declared tolerances, the file's motion), with the cache in `cache`.
fn answers(dir: &Path, cache: &Path, src: &str, scan: f32) -> Vec<(String, String, bool)> {
    std::fs::write(dir.join("m.ink"), src).unwrap();
    box_stl(&dir.join("scan.stl"), [0.0; 3], [scan; 3]);
    // SAFETY: alone in its process (`ran_alone`), and no ink work runs meanwhile.
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

/// Whether the test `name` (the caller, by its own name) was run in a
/// process of its own, and passed there: this binary given that test alone.
/// False in that process, which runs it. Where the cache is (XDG_CACHE_HOME,
/// INK_NO_CACHE, INK_CACHE_MB) and what it counts (cache::MISSES) are a
/// process's, so a test sets and reads them as it likes while the others
/// run side by side. The cache is on there whatever the caller's
/// environment says (INK_NO_CACHE): these are its tests.
fn ran_alone(name: &str) -> bool {
    if std::env::var_os("INK_CACHE_TEST").is_some_and(|t| t == name) {
        return false;
    }
    let out = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", name, "--nocapture", "--test-threads=1"])
        .env("INK_CACHE_TEST", name)
        .env_remove("INK_NO_CACHE")
        .output()
        .unwrap();
    let said = format!("{}{}", String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr));
    assert!(out.status.success() && said.contains("1 passed"), "{name}, in a process of its own:\n{said}");
    true
}

#[test]
fn a_claim_is_cached_exactly_when_its_question_is_unchanged() {
    if ran_alone("a_claim_is_cached_exactly_when_its_question_is_unchanged") {
        return;
    }
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
        ("the pin", MODEL.replace("(move (sphere 1) 20 0 0)", "(move (sphere 1.1) 20 0 0)"), 2.0, Some(&["apart", "walls", "pin-fine", "slides-clear"])),
        ("a claim's box", MODEL.replace("(-6 -3 -2) (6 3 2)", "(-6 -3 -2) (5.02 3 2)"), 2.0, Some(&["inside"])),
        ("a claim's bound", MODEL.replace("< 100", "< 70"), 2.0, None),
        // a ranged param is in the questions that read it alone (a disproof
        // names those it reads)
        ("a tolerance", MODEL.replace("(tol -0.1 0.1)", "(tol -0.1 0.3)"), 2.0, Some(&["inside", "apart", "small"])),
        ("a param's value", MODEL.replace("(param h 4)", "(param h 4.5)"), 2.0, Some(&["inside", "apart", "small"])),
        ("a claim's res", MODEL.replace("(22 2 2)) (res 0.5)", "(22 2 2)) (res 0.25)"), 2.0, Some(&["pin-fine"])),
        ("the t-span", MODEL.replace("(t-span 0 1s)", "(t-span 0 4s)"), 2.0, Some(&["slides-clear"])),
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
    let (before, first) = (misses(), ink::query::measured(&tape, &[None, None], true, &v.ctx, lo, hi, ink::query::Until::Accuracy(ink::query::ACCURACY)));
    let (after, again) = (misses(), ink::query::measured(&tape, &[None, None], true, &v.ctx, lo, hi, ink::query::Until::Accuracy(ink::query::ACCURACY)));
    assert_eq!((after - before, misses() - after), (1, 0), "measured once, then read");
    assert_eq!(format!("{first:?}"), format!("{again:?}"));
}

/// A quick box is keyed by the same complete question as a full bound,
/// including the region and ranged scope; a hit preserves every float bit.
#[test]
fn quick_bounds_are_cached_by_the_question_and_round_trip_exactly() {
    use ink::eval::Ctx;
    use ink::tape::Tape;
    use std::sync::Arc;

    if ran_alone("quick_bounds_are_cached_by_the_question_and_round_trip_exactly") {
        return;
    }
    let root = common::Scratch::new("quick-bounds-cache");
    // SAFETY: alone in its process (`ran_alone`), and no ink work runs meanwhile.
    unsafe { std::env::set_var("XDG_CACHE_HOME", &*root) };
    let doc =
        ink::load_str("(param r 2) (param unused 7) (body ball (shape (- (+ (square (- x (/ t 1s))) (square y) (square z)) (square r))))", "ball.ink", Path::new(".")).unwrap();
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
    // SAFETY: alone in its process (`ran_alone`), and no ink work runs meanwhile.
    unsafe { std::env::set_var("INK_NO_CACHE", "1") };
    let uncached = ask(&ctx, lo, hi);
    unsafe { std::env::remove_var("INK_NO_CACHE") };
    assert_eq!(uncached, (first, 1));
}

/// A question reads what its tapes' evaluation reads: a simplified tape
/// keeps a param's arithmetic in real numbers where it folded the param
/// itself away. p ranged over [-2, -1] picks 1e16, and (1e16 + 1) - 1e16
/// - 1/2 is 1/2 in real numbers, nothing inside, but -1/2 in floats.
#[test]
fn a_param_a_simplified_tape_folded_away_is_still_read() {
    use ink::eval::{Ctx, ParamSpans};
    use ink::graph::{Graph, Var};
    use ink::ranged::{ProofOptions, Status};
    use ink::tape::{Op, Tape};
    use std::sync::Arc;

    if ran_alone("a_param_a_simplified_tape_folded_away_is_still_read") {
        return;
    }
    let root = common::Scratch::new("folded-param");
    // SAFETY: alone in its process (`ran_alone`), and no ink work runs meanwhile.
    unsafe { std::env::set_var("XDG_CACHE_HOME", &*root) };
    let mut g = Graph::new();
    let p = g.var(Var::Param(0));
    let (big, bigger, one, half) = (g.constant(1e16), g.constant(2e16), g.constant(1.0), g.constant(0.5));
    let s = g.select(p, big, bigger);
    let s = g.add(s, one);
    let s = g.sub(s, big);
    let f = g.sub(s, half);
    let tape = Arc::new(Tape::new(&g, f));
    let ctx = Ctx { p_spans: ParamSpans::new([(0, (-2.0, -1.0))]), ..Ctx::new(0.0, vec![-1.5]) };
    let (lo, hi) = ([-1.0; 3], [1.0; 3]);
    let focused = ink::query::focus(&tape, &ctx, lo, hi);
    let opts = ProofOptions { res: None, finest: 0.1, depth: 4, tol: 0.0, boxes: ink::query::BOXES };
    for t in [&tape, &focused] {
        let p = ink::ranged::prove(t, &ctx, lo, hi, &opts);
        assert_ne!(p.status(), Status::Disproven, "{p:?}");
    }
    assert!(!focused.ssa.iter().any(|i| matches!(i.op, Op::Var(Var::Param(_)))), "the select is gone");
    assert_eq!(focused.params_read(), vec![0], "its arithmetic still reads p");
    assert_eq!(ctx.project(&[&focused]), ctx, "so p stays ranged");
}

/// u = w + t: a solve on t over a ranged w, which a shape that reads
/// neither does not ask.
struct Plus(Vec<ink::tape::Tape>);
impl ink::eval::Derive for Plus {
    fn owns(&self) -> &[u32] {
        &[1]
    }
    fn reads(&self) -> &[u32] {
        &[0]
    }
    fn tapes(&self) -> &[ink::tape::Tape] {
        &self.0
    }
    fn fingerprint(&self, _: &mut ink::graph::Fingerprint, _: &ink::tape::Names) -> bool {
        true
    }
    fn derive(&self, ctx: &mut ink::eval::Ctx) {
        let w = ctx.p_spans.get(0).unwrap();
        let t = ctx.t_span.unwrap_or((ctx.t, ctx.t));
        ctx.p_spans.set_span(1, (w.0 + t.0, w.1 + t.1));
        ctx.params[1] = ctx.params[0] + ctx.t;
    }
    fn rates(&self, _: &ink::eval::Ctx) -> Vec<Vec<f64>> {
        vec![vec![1.0]]
    }
}

/// Whether a question reads t is asked of what it reads: the key leaves
/// the moment out where its shape does not read t, and the witness read
/// back is then at the moment asked, though a solve beside it reads t.
#[test]
fn a_witness_read_back_is_at_the_moment_asked_beside_a_solve_on_t() {
    use ink::eval::{Ctx, ParamSpans};
    use ink::graph::{Graph, Var};
    use ink::query::Verdict;
    use ink::ranged::ProofOptions;
    use ink::tape::Tape;
    use std::sync::Arc;

    if ran_alone("a_witness_read_back_is_at_the_moment_asked_beside_a_solve_on_t") {
        return;
    }
    let root = common::Scratch::new("moment-asked");
    // SAFETY: alone in its process (`ran_alone`), and no ink work runs meanwhile.
    unsafe { std::env::set_var("XDG_CACHE_HOME", &*root) };
    let mut g = Graph::new();
    let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
    let (x2, y2, z2) = (g.mul(x, x), g.mul(y, y), g.mul(z, z));
    let (one, r2) = (g.constant(1.0), g.add(x2, y2));
    let r2 = g.add(r2, z2);
    let f = g.sub(r2, one);
    let ball = Arc::new(Tape::new(&g, f));
    let (w, u, t) = (g.var(Var::Param(0)), g.var(Var::Param(1)), g.var(Var::T));
    let wt = g.add(w, t);
    let eq = g.sub(u, wt);
    let spans = ParamSpans::new([(0, (0.9, 1.1))]).derived(Arc::new(Plus(vec![Tape::new(&g, eq)])), &[(1, (0.9, 1.1))]);
    let still = Ctx { p_spans: spans, ..Ctx::new(0.0, vec![1.0, 1.0]) };
    let mut over = Ctx { t: 10.5, t_span: Some((10.0, 11.0)), ..still.clone() };
    over.derive();
    let opts = ProofOptions { res: None, finest: 0.1, depth: 4, tol: 0.0, boxes: ink::query::BOXES };
    let (lo, hi) = ([-2.0; 3], [2.0; 3]);
    let moment = |c: &Ctx| match ink::ranged::prove(&ball, c, lo, hi, &opts).verdict {
        Verdict::Disproven { t, .. } => t,
        v => panic!("{v:?}"),
    };
    assert_eq!(moment(&still), 0.0);
    let before = ink::cache::here().0;
    let t = moment(&over);
    assert_eq!(ink::cache::here().0, before, "read back: the ball reads neither t nor u");
    assert!((10.0..=11.0).contains(&t), "the witness is at {t}, outside the span asked");
}

/// A grid is sampled once in a process: a letter met twice takes the
/// grid its first sampled (a load of "OO" computes what one of "O" does),
/// and so does a load again (a case's), the cache on disk off.
#[test]
fn a_grid_is_sampled_once_in_a_process() {
    if ran_alone("a_grid_is_sampled_once_in_a_process") {
        return;
    }
    let root = common::Scratch::new("grid-kept");
    // SAFETY: alone in its process (`ran_alone`), and no ink work runs meanwhile.
    unsafe { std::env::set_var("INK_NO_CACHE", "1") };
    let misses = || ink::cache::MISSES.load(std::sync::atomic::Ordering::Relaxed);
    // The results computed loading lettering of `text`, each letter's grid
    // sampled (letters of cap height h, so no other test's grids), and
    // loading it again.
    let load = |text: &str, h: f64| -> (usize, usize) {
        let f = root.join(format!("{text}.ink"));
        let src = format!(
            "(import \"std\")\n(def L (lettering \"{text}\" (font \"rounded\") {h}))\n\
             (body b (shape (fold min inf (map (lambda (q) (move (redistance (extrude-z (lettering-glyph L q) -1 1) (res 0.25)) (nth q 1) (nth q 2) 0)) (lettering-inked L)))))\n"
        );
        std::fs::write(&f, src).unwrap();
        let m0 = misses();
        let v = ModelView::load_in(ink::registry(), &f, &[], 0.0).unwrap_or_else(|e| panic!("{e:?}"));
        let m1 = misses();
        v.reload(&[]).unwrap_or_else(|e| panic!("{e:?}"));
        (m1 - m0, misses() - m1)
    };
    let (one, two) = (load("O", 9.5), load("OO", 9.25));
    // SAFETY: as above.
    unsafe { std::env::remove_var("INK_NO_CACHE") };
    assert!(one.0 > 0 && two == (one.0, 0) && one.1 == 0, "computed loading O, then again: {one:?}; OO: {two:?}");
}

/// INK_CACHE_MB caps the whole cache on disk (a hosted call's cache lives
/// in its memory, and the host sets it): 4 MB of answers stored under a
/// 1 MB cap leave at most 1 MB, the newest kept.
#[test]
fn the_cache_keeps_to_ink_cache_mb() {
    if ran_alone("the_cache_keeps_to_ink_cache_mb") {
        return;
    }
    let root = common::Scratch::new("cache-mb");
    // SAFETY: alone in its process (`ran_alone`), and no ink work runs meanwhile.
    unsafe {
        std::env::set_var("XDG_CACHE_HOME", root.join("cache"));
        std::env::set_var("INK_CACHE_MB", "1");
    }
    // Answers that deflate little: about 200 KB each on disk.
    let mut seed = 0x9e37_79b9_7f4a_7c15u64;
    let mut answer = || -> Vec<u64> {
        (0..25_000)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                seed
            })
            .collect()
    };
    for k in 0..20u128 {
        ink::cache::store("capped", k, &answer());
    }
    // SAFETY: as above.
    unsafe { std::env::remove_var("INK_CACHE_MB") };
    let dir = root.join("cache").join("ink").join("capped");
    let entries: Vec<std::fs::Metadata> = std::fs::read_dir(&dir).unwrap().map(|e| e.unwrap().metadata().unwrap()).collect();
    #[cfg(unix)]
    let on_disk: u64 = entries.iter().map(|m| std::os::unix::fs::MetadataExt::blocks(m) * 512).sum();
    // The cache budgets file lengths where Unix block counts are unavailable.
    #[cfg(not(unix))]
    let on_disk: u64 = entries.iter().map(std::fs::Metadata::len).sum();
    assert!(entries.iter().all(|m| m.len() > 150_000), "each answer takes about 200 KB");
    assert!(on_disk <= 1 << 20, "{} entries, {on_disk} bytes on disk, over INK_CACHE_MB=1", entries.len());
    assert!(dir.join(format!("{:032x}.json.z", 19u128)).exists(), "the newest is kept");
}
