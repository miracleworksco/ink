use super::*;

/// A view of `src` over the span, and its def `name`'s field.
fn view(src: &str, span: (f64, f64)) -> ModelView {
    static CORE: crate::Registry = crate::Registry::new(&[]);
    let doc = crate::Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), &[]).unwrap();
    let (ctx, _) = doc.ctx_at(&[], span.0).unwrap();
    ModelView::new(doc, ctx)
}

fn field(v: &ModelView, name: &str) -> u32 {
    let crate::Value::Field(f, _) = v.doc.def(name).unwrap().clone() else { panic!("{name} is no field") };
    f
}

/// The window `drive` gives the claim that the def `sa` is disjoint of
/// `sb`, in turns of the drive.
fn tooth(src: &str) -> Option<f64> {
    let mut v = view(src, (0.0, 1.0));
    let (a, b) = (field(&v, "sa"), field(&v, "sb"));
    let k = v.t_param();
    let d = drive(&mut v, &[a, b], (0.0, 1.0), k)?;
    // The drive turns 10 times over the span.
    Some(10.0 * (d.window.1 - d.window.0) / d.travel)
}

/// A claim reads a spin-up (t^2) only through the turns it drives: a
/// gear and a post repeat by one of the gear's teeth, a 15:45 pair by
/// one of the pinion's (the wheel turns a whole tooth with it), and an
/// arm fixed to the gear only by a whole turn. A gear tilted, turning
/// about its own z, repeats by a tooth; turned about its x, by a turn.
/// Waves of one drive summed repeat every turn of it; waves with no
/// common period leave the turn they drive repeating.
#[test]
fn a_drive_repeats_by_the_teeth_it_turns() {
    let head = "(import \"std\")\n(t-span 0 1s)\n(def drive (* 3600 (square (/ t 1s))))\n\
                (fn gear (n) (extrude-z (involute-gear-2d n 1 20 0.05) 0 5))\n\
                (body g12 (shape (gear 12)))\n(body g15 (shape (gear 15)))\n(body g45 (shape (gear 45)))\n(body post (shape (box 2 2 5)))\n";
    let near = |got: Option<f64>, want: f64| matches!(got, Some(w) if (w - want).abs() < 1e-9 * want);
    let post = format!("{head}(assembly root (instance g g12 (at) (turn z drive)) (instance p post (at (move 9 0 2.5))))\n(def sa (shape-of root/g))\n(def sb (shape-of root/p))");
    assert!(near(tooth(&post), 1.0 / 12.0), "{:?}", tooth(&post));
    let pair = format!(
        "{head}(assembly root (instance p g15 (at) (turn z drive)) (instance w g45 (at (move (mesh-distance 15 45 1) 0 0)) (turn z (mesh-turn 15 45 0 drive))))\n\
         (def sa (shape-of root/p))\n(def sb (shape-of root/w))"
    );
    assert!(near(tooth(&pair), 1.0 / 15.0), "{:?}", tooth(&pair));
    let arm = format!(
        "{head}(assembly geared (instance g g12) (instance a post (at (move 0 8 2.5))))\n\
         (assembly root (instance g geared (at) (turn z drive)) (instance p post (at (move 9 0 2.5))))\n\
         (def sa (shape-of root/g))\n(def sb (shape-of root/p))"
    );
    assert!(near(tooth(&arm), 1.0), "{:?}", tooth(&arm));
    let tilted = format!(
        "{head}(assembly root (instance g g12 (at (rot-x 10)) (turn z drive)) (instance p post (at (move 9 0 2.5))))\n(def sa (shape-of root/g))\n(def sb (shape-of root/p))"
    );
    assert!(near(tooth(&tilted), 1.0 / 12.0), "{:?}", tooth(&tilted));
    let off = format!("{head}(assembly root (instance g g12 (at) (turn x drive)) (instance p post (at (move 9 0 2.5))))\n(def sa (shape-of root/g))\n(def sb (shape-of root/p))");
    assert!(near(tooth(&off), 1.0), "{:?}", tooth(&off));
    // Two waves of one spin-up, summed: their branches rejoin above the
    // drive, which still repeats every turn.
    let mut v = view("(t-span 0 1s)\n(def a (* 3600000000 (square (/ t 1s))))\n(def f (+ (sin (deg->rad a)) (cos (deg->rad a))))", (0.0, 1.0));
    let (f, k) = (field(&v, "f"), v.t_param());
    let d = drive(&mut v, &[f], (0.0, 1.0), k).expect("a drive");
    assert!((1e7 * (d.window.1 - d.window.0) / d.travel - 1.0).abs() < 1e-9, "{:?}", d.window);
    // Waves of it with no common period (sin of a in radians, cos of it
    // in degrees) drive an outer turn b: tried deepest first, the
    // waves' drive is shown no period and b's turn is.
    let mut v = view(
        "(t-span 0 1s)\n(def a (* 3600000000 (square (/ t 1s))))\n(def b (* 1000000000 (+ (sin a) (cos (deg->rad a)))))\n\
         (def f (sin (deg->rad b)))",
        (0.0, 1.0),
    );
    let (f, k) = (field(&v, "f"), v.t_param());
    let d = drive(&mut v, &[f], (0.0, 1.0), k).expect("a drive");
    assert!(d.window.1 - d.window.0 < 1e-6 * d.travel, "{:?}", d.window);
}

/// Finding a drive costs a few readings of the claim, however many
/// values of t it reads in a row: a value 4,096 additions from t that
/// repeats in nothing is read through t and the top of the run (the
/// claim), not again through each addition.
#[test]
fn a_drive_is_found_in_a_few_readings() {
    let n = 4096;
    let src: String = std::iter::once("(t-span 0 1s)\n(def a0 (/ t 1s))\n".to_string()).chain((1..=n).map(|i| format!("(def a{i} (+ a{} 1))\n", i - 1))).collect();
    let mut v = view(&src, (0.0, 1.0));
    let (f, k) = (field(&v, &format!("a{n}")), v.t_param());
    let tries = tries(&v, &[f], k).expect("t is read");
    assert!(tries.len() <= 2, "{} readings", tries.len());
    assert!(drive(&mut v, &[f], (0.0, 1.0), k).is_none());
}

/// A frame's witness and pair, as the search in it gives them (frame
/// coordinates; t's param absent where the frame cancels t), checked
/// in the world by `Frame::confirmed`.
fn witnesses(v: &mut ModelView, claim: Claim, span: (f64, f64), y: V3) -> (Option<Evidence>, Option<Evidence>) {
    use ink_kernel::query::Verdict;
    use ink_kernel::ranged::Proof;
    let (_, world, _) = moving_frame(v, &claim, span, None).expect("a moving frame");
    let k = v.t_param();
    let fr = Frame::new(v, &claim, world, k, span, None, 0.0);
    let proof = Proof { verdict: Verdict::Disproven { witness: y, t: 0.0, value: -0.1, params: vec![] }, res: 0.05, tol: 0.05, exact: false, loose: Vec::new(), deepest: None };
    let pair = Evidence::ClearanceSearch { verdict: Verdict::Disproven { witness: y, t: 0.0, value: 0.0, params: vec![] }, pair: Some((y, y)), res: 0.05 };
    (fr.confirmed(v, Evidence::Empty(proof)), fr.confirmed(v, pair))
}

/// Two balls 0.1 apart along x carried 1e16 out: in the first's frame
/// their overlap's point (1.05, 0, 0) is a witness, but moved into the
/// world in floats it rounds to x = 1e16, outside both: neither it nor a
/// pair there is reported. Turned together by two instances instead, the frame
/// cancels t: the witness is given at the span's start, in the world
/// there (R(90 deg) (10.25, 0, 0)), and the pair too.
#[test]
fn a_frame_witness_is_reported_only_where_the_world_shows_it() {
    let far = "(import \"std\")\n(t-span 0 1s)\n(body a (shape (move (sphere 0.25) 1 0 0)))\n(body b (shape (move (sphere 0.25) 1.1 0 0)))\n\
               (assembly root (instance a a (at (move (+ 1e16 (* 0.001 (/ t 1s))) 0 0))) (instance b b (at (move (+ 1e16 (* 0.001 (/ t 1s))) 0 0))))\n\
               (def sa (shape-of root/a))\n(def sb (shape-of root/b))";
    let mut v = view(far, (0.0, 1.0));
    let (a, b) = (field(&v, "sa"), field(&v, "sb"));
    let (e, p) = witnesses(&mut v, Claim::Disjoint(a, b), (0.0, 1.0), [1.05, 0.0, 0.0]);
    assert!(e.is_none() && p.is_none(), "{e:?} {p:?}");
    let (_, p) = witnesses(&mut v, Claim::Clearance(a, b, 0.5), (0.0, 1.0), [1.05, 0.0, 0.0]);
    assert!(p.is_none(), "{p:?}");

    let turned = "(import \"std\")\n(t-span 0.25s 1s)\n(body a (shape (move (sphere 1) 10 0 0)))\n(body b (shape (move (sphere 1) 10.5 0 0)))\n\
                  (assembly root (instance a a (at (rot-z (* 360 (/ t 1s))))) (instance b b (at (rot-z (* 360 (/ t 1s))))))\n\
                  (def sa (shape-of root/a))\n(def sb (shape-of root/b))";
    let mut v = view(turned, (0.25, 1.0));
    let (a, b) = (field(&v, "sa"), field(&v, "sb"));
    let k = v.t_param();
    let (e, p) = witnesses(&mut v, Claim::Disjoint(a, b), (0.25, 1.0), [10.25, 0.0, 0.0]);
    let Some(Evidence::Empty(ink_kernel::ranged::Proof { verdict: ink_kernel::query::Verdict::Disproven { witness, t, params, .. }, .. })) = e else { panic!("{e:?}") };
    assert!(t == 0.25 && params.contains(&(k, 0.25)) && witness[0].abs() < 1e-9 && (witness[1] - 10.25).abs() < 1e-9, "{witness:?} at {t}, {params:?}");
    let Some(Evidence::ClearanceSearch { pair: Some((p, q)), verdict: ink_kernel::query::Verdict::Disproven { t, .. }, .. }) = p else { panic!("{p:?}") };
    assert!(t == 0.25 && p == q && (p[1] - 10.25).abs() < 1e-9, "{p:?} {q:?} at {t}");
}

/// A sampled failure stays a guarantee only where the claim asked of its
/// own nodes, t a param of the moment's one value, fails too (`certified`):
/// one a pointwise check found in floats, ((1e16 + t) - 1e16 > 0.001)
/// failing at its value 0, is an estimate (the real value is t, 0.0156
/// here); (t > 0.02) failing at t = 0.0156 is certified.
#[test]
fn a_sampled_failure_is_a_guarantee_only_once_its_own_claim_fails_at_the_moment() {
    let mut v = view("(def far (- (+ 1e16 (/ t 1s)) 1e16))\n(def near t)", (0.01, 0.02));
    v.set_moment(0.015625).unwrap();
    let sampled = |x: f64| Check { accuracy: Some(1e-6), ..checked(Status::Fails, Evidence::Value(x)) };
    let opts = ClaimOptions::default();
    let (far, near) = (field(&v, "far"), field(&v, "near"));
    let c = certified(&mut v, &Claim::Value(far, Cmp::Gt, 0.001), None, sampled(0.0), opts);
    assert!(c.status == Status::Fails && c.accuracy.is_none(), "{c:?}");
    let c = certified(&mut v, &Claim::Value(near, Cmp::Gt, 0.02), None, sampled(0.015625), opts);
    assert!(c.accuracy == Some(1e-6), "{c:?}");
}

/// tests/motion/loop-closed-pin.ink's coupler, seen from its rocker,
/// turns about the pin they share, though the eye is placed through the
/// coupler's angle and the pin through the rocker's: the eye is placed
/// about that pivot, its drift bounded far under the accuracy and
/// charged. Seen from a coupler whose eye is 0.3 off the pin, the pin
/// turns about the point the rocker shares with that coupler: asked
/// there, where nothing of the claim moves, its drift charged. A union
/// of the eye and a part at rest is not one shape moved, so it is not
/// charged one move and is placed as it is.
#[test]
fn a_loop_closed_joint_is_searched_about_its_pivot() {
    let src = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/motion/loop-closed-pin.ink")).unwrap();
    let src = format!(
        "{src}\n(def pin (max (shape-of root/rocker rocker-pin) (slab-z 11 15)))\n(def in-eye (shape-of root/coupler eye))\n(def in-far (shape-of root/long far-eye))\n\
         (def coupler-w (shape-of root/coupler))\n(def rocker-w (shape-of root/rocker))\n(def with-wall (min (shape-of root/coupler eye) (move (box 1 1 1) 100 0 0)))"
    );
    let mut v = view(&src, (0.0, 1.0));
    let [pin, in_eye, in_far, coupler, rocker, with_wall] = ["pin", "in-eye", "in-far", "coupler-w", "rocker-w", "with-wall"].map(|n| field(&v, n));
    let charge = |v: &mut ModelView, c: Claim| moving_frame(v, &c, (0.0, 1.0), Some(0.05)).expect("a moving frame").2;
    let seated = charge(&mut v, Claim::Subset(pin, in_eye));
    assert!(seated > 0.0 && seated < 1e-3, "{seated}");
    let clear = charge(&mut v, Claim::Clearance(coupler, rocker, 0.05));
    assert!(clear > 0.0 && clear < 2e-3, "{clear}");
    let far = charge(&mut v, Claim::Subset(pin, in_far));
    assert!(far > 0.0 && far < 1e-3, "{far}");
    assert_eq!(charge(&mut v, Claim::Disjoint(pin, with_wall)), 0.0);
    assert_eq!(charge(&mut v, Claim::Subset(pin, in_eye)), seated);
}

/// A hold of shapes asked off the world's (a pivot's drift, a drive's
/// slack, charged to the accuracy) is to the accuracy, never outright:
/// the world's may overlap by up to the charge. Two cylinders 0.01 into
/// each other at t = 0, one turning about the other's axis on a centre
/// that drifts 0.04 (its pivot), as asserted and as interfere's pair;
/// and a 15:45 gear pair apart, through a drive whose wheel turns at a
/// float rate (its slack).
#[test]
fn a_charged_hold_is_never_outright() {
    let src = "(import \"std\")\n(t-span 0 1s)\n(body a (shape (cylinder-z 1 2)))\n(body b (shape (cylinder-z 1 2)))\n\
               (def shift (- 2.01 (* 0.02 (cos (deg->rad (* 360 (/ t 1s)))))))\n\
               (assembly root (instance a a (at (rot-z (* 360 (/ t 1s))))) (instance b b (at (rot-z (* 720 (/ t 1s))) (move shift 0 0) (rot-z (* 360 (/ t 1s))))))\n\
               (require apart (disjoint (shape-of root/a) (shape-of root/b)))";
    let gears = "(import \"std\")\n(t-span 0 1s)\n(def drive (* 3600 (square (/ t 1s))))\n(fn gear (n) (extrude-z (involute-gear-2d n 1 20 0.05) 0 5))\n\
                 (body g15 (shape (gear 15)))\n(body g45 (shape (gear 45)))\n\
                 (assembly root (instance p g15 (at) (turn z drive)) (instance w g45 (at (move (+ 2 (mesh-distance 15 45 1)) 0 0)) (turn z (mesh-turn 15 45 0 drive))))\n\
                 (require apart (disjoint (shape-of root/p) (shape-of root/w)))";
    let held = |src: &str| {
        let mut v = view(src, (0.0, 1.0));
        v.set_motion(v.doc.t_span).unwrap();
        let r = claims::run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
        r.asserts[0].check.clone()
    };
    for src in [src, gears] {
        let c = held(src);
        assert!(c.holds_guaranteed() && matches!(c.held, Some(claims::Held::ToAccuracy { .. })), "{c:?}");
    }
    let mut v = view(src, (0.0, 1.0));
    v.set_motion(v.doc.t_span).unwrap();
    let (items, target) = (v.doc.items("root").unwrap(), v.doc.target("root"));
    let plan = crate::interfere::plan(&mut v, items, target.as_ref(), None, None).unwrap();
    let pairs = plan.run(&mut v, None, &|_| None).unwrap();
    let crate::interfere::Answer::Moving { check } = &pairs[0].answer else { panic!("a moving pair") };
    assert!(check.holds_guaranteed() && matches!(check.held, Some(claims::Held::ToAccuracy { .. })), "{check:?}");
}

/// An estimate's guard on t (estimate-when) is shown again at each
/// moment over the params' ranges, by interval bounds: a time-only one
/// beside an unrelated ranged param, at every moment; one reading t and
/// ranged q, (q + t > 0.905), where it holds at every q, and refused
/// with its reason where its corner (q = 0.9, t below 0.005) fails,
/// though it holds at q's value.
#[test]
fn an_estimates_guard_on_t_is_shown_at_each_moment_over_the_ranges() {
    let src = "(t-span 0 1s)\n(param q 1 (tol -.1 .1))\n(param w 5 (tol -.1 .1))
               (def time-only (estimate-when 2 (> (+ .5 (/ t 1s)) 0) \"positive\"))
               (def both (estimate-when 2 (> (+ q (/ t 1s)) .905) \"above 0.905\"))";
    let mut v = view(src, (0.0, 1.0));
    let (q, w) = (v.doc.param_index("q").unwrap(), v.doc.param_index("w").unwrap());
    v.set_param_spans(&[(q, (0.9, 1.1)), (w, (4.9, 5.1))]).unwrap();
    let (a, b) = (field(&v, "time-only"), field(&v, "both"));
    v.want(&[a]).unwrap();
    let refused = v.want(&[b]).expect_err("q = 0.9 at t = 0");
    assert!(refused.message.contains("above 0.905"), "{refused:?}");
    for (t, holds) in [(0.25, true), (1.0, true), (0.004, false), (0.006, true), (0.0, false)] {
        v.set_moment(t).unwrap();
        assert!(v.refusal(&[a]).is_none(), "{t}");
        match v.refusal(&[b]) {
            None => assert!(holds, "{t}: q + t > 0.905 fails at q = 0.9"),
            Some(e) => assert!(!holds && e.message.contains("above 0.905"), "{t}: {e:?}"),
        }
    }
}

/// A guard reached at a moment by `want` and by `set_moment` (stepping
/// away and back) is shown the same way, beside an unrelated ranged
/// param: (1e16 + 2 t) - 1e16 > 1.75 at t = 0.75, which floats and
/// interval bounds decide differently, says one thing every way.
#[test]
fn a_guard_says_the_same_at_a_moment_however_it_is_reached() {
    let src = "(t-span 0 1s)\n(param w 5 (tol -.1 .1))
               (def e (estimate-when 2 (> (- (+ 1e16 (* 2 (/ t 1s))) 1e16) 1.75) \"rounding\"))";
    let mut v = view(src, (0.75, 1.0));
    let w = v.doc.param_index("w").unwrap();
    v.set_param_spans(&[(w, (4.9, 5.1))]).unwrap();
    let e = field(&v, "e");
    let first = v.want(&[e]).is_ok();
    v.set_moment(1.0).unwrap();
    v.set_moment(0.75).unwrap();
    let stepped = v.refusal(&[e]).is_none();
    let again = v.want(&[e]).is_ok();
    assert!(first == stepped && stepped == again, "want {first}, set_moment {stepped}, want again {again}");
}

/// A guard on t and ranged q already wanted (at q's value, 1), then
/// scoped over q's range and a motion: the scope is taken and the guard
/// shown again in it at once, then at each moment read (a corner it
/// fails at, q = 0.9 at t = 0, refuses what reads it there, not the scope).
#[test]
fn a_wanted_guard_on_t_takes_a_scope_over_its_range_and_a_motion() {
    let src = "(t-span 0 1s)\n(param q 1 (tol -.1 .1))
               (def holds (estimate-when 2 (> (+ q (/ t 1s)) 0) \"positive\"))
               (def corner (estimate-when 2 (> (+ q (/ t 1s)) .905) \"above 0.905\"))";
    let mut v = view(src, (0.0, 1.0));
    let q = v.doc.param_index("q").unwrap();
    let (a, b) = (field(&v, "holds"), field(&v, "corner"));
    v.want(&[a, b]).unwrap();
    let scope = crate::view::ClaimScope { motion: Some((0.0, 1.0)), by_file: false, spans: vec![(q, (0.9, 1.1))] };
    v.set_claim_scope(&scope).unwrap();
    // Shown again in the scope just set: no want or step between.
    assert!(v.refusal(&[a]).is_none());
    assert!(v.refusal(&[b]).is_some_and(|e| e.message.contains("above 0.905")), "q = 0.9 at t = 0");
    v.want(&[a]).unwrap();
    assert!(v.want(&[b]).is_err_and(|e| e.message.contains("above 0.905")));
    v.set_moment(0.5).unwrap();
    assert!(v.refusal(&[a]).is_none() && v.refusal(&[b]).is_none());
}

/// The moving items' boxes are found side by side, each once, and each is
/// the box found of it alone, at one thread or four: wheels turning about
/// their axes (copies of one body, swept), one about a tilted axis (a turn
/// of floats, no exact isometry: its bounds over the span), a belt growing
/// with t, and a post at rest.
#[test]
fn moving_items_are_boxed_alike_side_by_side_and_alone() {
    let src = "(import \"std\")\n(t-span 0 1s)\n(def drive (* 90 (/ t 1s)))\n\
               (body wheel (shape (move (box 4 1 1) 2 0 0)))\n(body post (shape (box 1 1 1)))\n\
               (body belt (shape (box (+ 2 (* 20 (/ t 1s))) 2 2)))\n\
               (assembly root (instance a wheel (at) (turn z drive)) (instance b wheel (at (move 0 0 5)) (turn z drive))\n\
               (instance c wheel (at (rot 0.57735 0.57735 0.57735 120) (move 0 0 10)) (turn z (/ drive 9))) (instance d belt (at (move 0 20 0)))\n\
               (instance p post (at (move 9 9 9))))";
    let planned = || {
        let mut v = view(src, (0.0, 1.0));
        v.set_motion(Some((0.0, 1.0))).unwrap();
        let items = v.doc.items("root").unwrap();
        (v, items)
    };
    let shown = |b: &[ItemBox]| format!("{b:?}");
    let on = |threads: usize| {
        rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap().install(|| {
            let (mut v, items) = planned();
            item_boxes(&mut v, &items)
        })
    };
    let (one, four) = (on(1), on(4));
    assert_eq!(shown(&one), shown(&four));
    let (mut v, items) = planned();
    let alone: Vec<ItemBox> = (0..items.len()).flat_map(|k| item_boxes(&mut v, &items[k..=k])).collect();
    assert_eq!(shown(&one), shown(&alone));
    let label = |k: usize| items[k].label.rsplit('/').next().unwrap().to_string();
    let by: std::collections::HashMap<String, &ItemBox> = (0..items.len()).map(|k| (label(k), &one[k])).collect();
    // Both paths are taken: swept pieces for an exact turn, the span's
    // bounds for the tilted one; the belt holds its growth; the post is still.
    assert!(!by["a"].2.is_empty() && !by["b"].2.is_empty(), "{by:?}");
    assert!(by["c"].2.is_empty() && by["c"].0.is_some() && by["d"].0.is_some_and(|b| b.1[0] >= 11.0), "{by:?}");
    assert!(by["p"].2.is_empty() && by["p"].0.is_some(), "{by:?}");
}
