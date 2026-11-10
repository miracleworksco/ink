//! Fonts: what a glyph fills against an independent reference, every
//! glyph of every vendored font built, text shaped as HarfBuzz shapes it,
//! lettering-max's bound, and text params.

use super::region::{Edge, filled};
use super::*;
static CORE: crate::Registry = crate::Registry::new(&[]);

/// Character c's glyph in font s, as its cmap maps it.
fn gid(d: &mut Doc, s: &Spec, c: char) -> u16 {
    d.with_face(s, |f| Ok(f.glyph_index(c).expect("a glyph"))).unwrap()
}

fn square(x0: f64, y0: f64, s: f64, ccw: bool) -> Vec<Edge> {
    let mut p = [[x0, y0], [x0 + s, y0], [x0 + s, y0 + s], [x0, y0 + s]];
    if !ccw {
        p.reverse();
    }
    (0..4).map(|i| Edge { a: p[i], b: p[(i + 1) % 4], bu: 0.0 }).collect()
}

fn loops_area(ls: &[Vec<Edge>]) -> f64 {
    ls.iter().map(|l| area(&l.iter().map(|e| [e.a[0], e.a[1], e.bu]).collect::<Vec<_>>())).sum()
}

#[test]
fn overlapping_loops_fill_their_union() {
    // Two squares overlapping by a quarter: one outline, area 7/4.
    let u = filled(&[square(0.0, 0.0, 1.0, true), square(0.5, 0.5, 1.0, true)]);
    assert_eq!(u.len(), 1, "{u:?}");
    assert!((loops_area(&u) - 1.75).abs() < 1e-12, "{}", loops_area(&u));
    // A square with a hole wound the other way: two loops, area 1 - 1/4.
    let h = filled(&[square(0.0, 0.0, 1.0, true), square(0.25, 0.25, 0.5, false)]);
    assert!(h.len() == 2 && (loops_area(&h) - 0.75).abs() < 1e-12, "{h:?}");
    // The same hole wound as an outside is covered twice: filled.
    let f = filled(&[square(0.0, 0.0, 1.0, true), square(0.25, 0.25, 0.5, true)]);
    assert!(f.len() == 1 && (loops_area(&f) - 1.0).abs() < 1e-12, "{f:?}");
}

#[test]
fn a_loop_crossing_itself_fills_both_lobes() {
    // A bow tie: (0,0) (2,1) (2,0) (0,1), crossing at (1, 0.5). Its lobes
    // wind opposite ways; nonzero fills both, as two simple loops.
    let p = [[0.0, 0.0], [2.0, 1.0], [2.0, 0.0], [0.0, 1.0]];
    let tie: Vec<Edge> = (0..4).map(|i| Edge { a: p[i], b: p[(i + 1) % 4], bu: 0.0 }).collect();
    let f = filled(&[tie]);
    assert_eq!(f.len(), 2, "{f:?}");
    for l in &f {
        let a = area(&l.iter().map(|e| [e.a[0], e.a[1], e.bu]).collect::<Vec<_>>());
        assert!((a - 0.5).abs() < 1e-12, "each lobe counter-clockwise, area 1/2: {a}");
    }
}

#[test]
fn arcs_crossing_are_cut_where_they_meet() {
    // Two discs of radius 1 (each two half-circle arcs), centres 1 apart:
    // the union's area is 2 pi - (2 pi / 3 - sqrt(3) / 2).
    let disc = |cx: f64| vec![Edge { a: [cx + 1.0, 0.0], b: [cx - 1.0, 0.0], bu: 1.0 }, Edge { a: [cx - 1.0, 0.0], b: [cx + 1.0, 0.0], bu: 1.0 }];
    let u = filled(&[disc(0.0), disc(1.0)]);
    let want = 2.0 * std::f64::consts::PI - (2.0 * std::f64::consts::PI / 3.0 - 3f64.sqrt() / 2.0);
    assert!(u.len() == 1 && (loops_area(&u) - want).abs() < 1e-9, "{} {want}", loops_area(&u));
}

/// Glyph areas (em²) of skia-pathops' overlap removal (fontTools'
/// DecomposingRecordingPen into pathops.Path.simplify(fix_winding=True),
/// AreaPen), at the axes given: the variable fonts' glyphs overlap and
/// cross themselves.
const REFERENCE: &[(&str, &[(&str, f64)], char, f64)] = &[
    ("sans", &[("wght", 400.0)], 'A', 0.172527),
    ("sans", &[("wght", 400.0)], 'B', 0.231041),
    ("sans", &[("wght", 400.0)], '8', 0.185661),
    ("sans", &[("wght", 400.0)], 'g', 0.206453),
    ("sans", &[("wght", 400.0)], '@', 0.307637),
    ("sans", &[("wght", 900.0)], 'A', 0.276825),
    ("sans", &[("wght", 900.0)], 'B', 0.344303),
    ("sans", &[("wght", 900.0)], 'g', 0.335399),
    ("rounded", &[("wght", 700.0), ("wdth", 100.0)], 'A', 0.258757),
    ("rounded", &[("wght", 700.0), ("wdth", 100.0)], 'M', 0.355333),
    ("rounded", &[("wght", 700.0), ("wdth", 100.0)], '6', 0.245436),
    ("script", &[], 'A', 0.321367),
    ("script", &[], 'g', 0.233403),
    ("serif", &[("wght", 300.0)], 'A', 0.088361),
    ("serif", &[("wght", 300.0)], '&', 0.093636),
];

#[test]
fn a_glyph_fills_what_the_font_fills() {
    let mut d = Doc::load_str_in(&CORE, "", "t", std::path::Path::new("."), &[]).unwrap();
    for (font, axes, c, want) in REFERENCE {
        let s = Spec { key: font.to_string(), axes: axes.iter().map(|(k, v)| (k.to_string(), *v)).collect() };
        let id = gid(&mut d, &s, *c);
        let g = d.glyph_outline(&s, id).unwrap();
        let a: f64 = g.loops.iter().map(|l| area(l)).sum();
        // The reference flattens curves to lines: within 0.1%.
        assert!((a - want).abs() < 1e-3 * want, "{font} {axes:?} {c}: {a} against {want}");
    }
}

#[test]
fn every_printable_glyph_of_every_vendored_font_is_a_profile() {
    for (name, _, _) in VENDORED {
        let axes: Vec<&str> = match *name {
            "sans" | "serif" => vec!["(wght 100)", "(wght 400)", "(wght 900)"],
            "rounded" => vec!["(wght 300) (wdth 75)", "(wght 700) (wdth 125)"],
            _ => vec![""],
        };
        for a in axes {
            let src: String = format!("(def f (font \"{name}\" {a}))\n") + &(33..127).map(|c| format!("(def g{c} (glyph-2d f {c}))\n")).collect::<String>();
            if let Err(e) = Doc::load_str_in(&CORE, &src, "t", std::path::Path::new("."), &[]) {
                panic!("{name} {a}: {}", e.message);
            }
        }
    }
}

#[test]
fn kerning_and_metrics_are_the_fonts() {
    // glyph-2d and glyph-metrics take a character's code: A V kerns closer,
    // A A not; O has a hole; a space no ink. Fredoka bold condensed kerns v
    // e by 0 (its -30 units' variation delta, +30 there, read from the pair
    // set's start: from its records it was lost and v e kerned by -0.03).
    let src = "(def f (font \"sans\"))\n(def av (get (glyph-metrics f 65 86) \"kern\"))\n(def aa (get (glyph-metrics f 65 65) \"kern\"))\n\
               (def hole (len (get (glyph-metrics f 79) \"holes\")))\n(def space (get (glyph-metrics f 32) \"ink\"))\n\
               (def ve (get (glyph-metrics (font \"rounded\" (wght 700) (wdth 75)) 118 101) \"kern\"))\n\
               (def ve-default (get (glyph-metrics (font \"rounded\") 118 101) \"kern\"))\n\
               (def a-run (get (nth (nth (nth (nth (glyph-run f \"A\") 0) 0) 0) 1) \"box\"))\n(def a-code (get (glyph-metrics f 65) \"box\"))\n";
    let d = Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), &[]).unwrap();
    let num = |n: &str| match d.def(n) {
        Some(Value::Field(f, _)) => d.g.as_const(*f).unwrap(),
        v => panic!("{n}: {v:?}"),
    };
    assert!(num("av") < -0.01 && num("aa") == 0.0, "A V kerns closer, A A not: {} {}", num("av"), num("aa"));
    assert!(num("hole") == 1.0 && num("space") == 0.0);
    assert_eq!((num("ve"), num("ve-default")), (0.0, -0.03));
    let nums = |n: &str| -> Vec<f64> {
        let Some(Value::Vec(v)) = d.def(n) else { panic!("{n}") };
        v.iter().map(|x| if let Value::Field(g, _) = x { d.g.as_const(*g).unwrap() } else { panic!() }).collect()
    };
    assert_eq!(nums("a-run"), nums("a-code"), "code 65 is A, as glyph-run sets it");
}

#[test]
fn a_line_is_shaped_as_harfbuzz_shapes_it() {
    let cases: serde_json::Value = serde_json::from_str(include_str!("shaped.json")).unwrap();
    let mut d = Doc::load_str_in(&CORE, "", "t", std::path::Path::new("."), &[]).unwrap();
    for c in cases.as_array().unwrap() {
        let axes = c["axes"].as_object().unwrap().iter().map(|(k, v)| (k.clone(), v.as_f64().unwrap())).collect();
        let s = Spec { key: c["font"].as_str().unwrap().to_string(), axes };
        let text = c["text"].as_str().unwrap();
        let font = d.shaper(&s).unwrap();
        let (placed, adv, n) = shape::shape(&font, &s.key, text).unwrap();
        let want = c["glyphs"].as_array().unwrap();
        let what = format!("{} {:?} {text:?}", s.key, s.axes);
        assert_eq!(placed.len(), want.len(), "{what}: {} glyphs, HarfBuzz {}", placed.len(), want.len());
        let (mut pen, mut clusters) = (0.0, Vec::new());
        for (p, w) in placed.iter().zip(want) {
            let w: Vec<&serde_json::Value> = w.as_array().unwrap().iter().collect();
            let num = |i: usize| w[i].as_f64().unwrap();
            if clusters.last() != Some(&num(4)) {
                clusters.push(num(4));
            }
            let got = (p.id as f64, p.x, p.y, p.cluster);
            assert_eq!(got, (num(0), pen + num(2), num(3), clusters.len() - 1), "{what}: glyph {}", w[5]);
            pen += num(1);
        }
        assert_eq!((adv, n), (pen, clusters.len()), "{what}");
    }
}

#[test]
fn lettering_sets_ligatures_and_accents_as_the_font_does() {
    // office in Lexend is o, the ffi ligature, c, e: four glyphs. naive
    // with its i's diaeresis typed composed or as i + U+0308 sets the same
    // in Fredoka bold condensed, five glyphs, its v one i's advance (251
    // units) past its ï (HarfBuzz's positions: 0, 497, 998, 1249, 1753).
    let src = "(import \"std\")\n(def office (len (nth (lettering \"office\" (font \"sans\") 7) 2)))\n\
               (def F (font \"rounded\" (wght 700) (wdth 75)))\n(def k (/ 1 (get F \"cap-height\")))\n\
               (fn xs (t) (map (lambda (g) (/ (- (nth g 1) (nth (nth (nth (lettering-styled t F 1 0 1 0) 2) 0) 1)) k)) (nth (lettering-styled t F 1 0 1 0) 2)))\n\
               (def composed (xs \"na\u{ef}ve\"))\n(def combining (xs \"nai\u{308}ve\"))\n";
    let d = Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), &[]).unwrap();
    let nums = |n: &str| -> Vec<f64> {
        match d.def(n) {
            Some(Value::Field(f, _)) => vec![d.g.as_const(*f).unwrap()],
            Some(Value::Vec(v)) => v.iter().map(|x| if let Value::Field(f, _) = x { (d.g.as_const(*f).unwrap() * 1e6).round() / 1e3 } else { panic!() }).collect(),
            v => panic!("{n}: {v:?}"),
        }
    };
    assert_eq!(nums("office"), vec![4.0]);
    assert_eq!(nums("composed"), vec![0.0, 497.0, 998.0, 1249.0, 1753.0]);
    assert_eq!(nums("combining"), nums("composed"));
}

#[test]
fn lettering_max_bounds_every_line_its_characters_make() {
    // Lines of up to 6 characters drawn from each set (ligatures, accents,
    // a script's joining forms among them), at three trackings: each one's
    // box within the bound, its width and how far it reaches below and
    // above the baseline (y = -h/2 for one line of cap height h = 1).
    let sets = [("sans", "office AVTWayfl\u{ef}\u{e9}"), ("script", "Grace Lovelaffi\u{eb}"), ("serif", "affluent Zo\u{eb}"), ("rounded", "office na\u{ef}veQ")];
    let mut seed = 12345u64;
    let mut next = |n: usize| {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 33) as usize % n
    };
    for (font, chars) in sets {
        let cs: Vec<char> = chars.chars().collect();
        let mut src = format!("(import \"std\")\n(def F (font \"{font}\"))\n");
        let mut lines = Vec::new();
        for (ti, t) in [0.0, 0.1, -0.05].into_iter().enumerate() {
            src += &format!("(def m{ti} (lettering-max \"{chars}\" 6 F 1 {t}))\n");
            for i in 0..60 {
                let line: String = (0..1 + next(6)).map(|_| cs[next(cs.len())]).collect();
                if line.trim().is_empty() {
                    continue;
                }
                src += &format!("(def b{ti}x{i} (lettering-box (lettering-styled \"{line}\" F 1 {t} 1 0)))\n");
                lines.push((ti, i, line));
            }
        }
        let d = Doc::load_str_in(&CORE, &src, "t", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{font}: {}", e.message));
        let v = |n: &str| -> Vec<f64> {
            let Some(Value::Vec(v)) = d.def(n) else { panic!("{n}") };
            v.iter().map(|x| if let Value::Field(f, _) = x { d.g.as_const(*f).unwrap() } else { panic!() }).collect()
        };
        for (ti, i, line) in lines {
            let (m, b) = (v(&format!("m{ti}")), v(&format!("b{ti}x{i}")));
            // The bound is tight below and above (the same glyph's box):
            // within float rounding of moving it down 0.5.
            let ok = b[2] - b[0] <= m[0] && b[1] + 0.5 >= m[1] - 1e-12 && b[3] + 0.5 <= m[2] + 1e-12;
            assert!(ok, "{font} {line:?} tracking #{ti}: box {b:?} outside the bound {m:?}");
        }
    }
}

#[test]
fn lettering_max_refuses_what_it_cannot_bound() {
    let e = |src: &str| Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), &[]).err().map(|e| e.message).unwrap_or_default();
    let m = e("(import \"std\")\n(def m (lettering-max \"ae\u{301}\" 4 (font \"sans\") 1 0))");
    assert!(m.contains("U+0301") && m.contains("combining mark"), "{m}");
    // U+01C4 (DŽ) Lexend splits into D and Ž before it sets them.
    let m = e("(import \"std\")\n(def m (lettering-max \"a\u{1c4}\" 4 (font \"sans\") 1 0))");
    assert!(m.contains("(U+01C4) sets as 2 glyphs"), "{m}");
}

#[test]
fn a_missing_glyph_and_axis_are_errors_naming_them() {
    let e = |src: &str| Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), &[]).err().map(|e| e.message).unwrap_or_default();
    assert!(e("(def g (glyph-2d (font \"sans\") 20013))").contains("has no glyph for '中' (U+4E2D)"));
    assert!(e("(def f (font \"script\" (wght 700)))").contains("has no axis 'wght'"));
    assert!(e("(def f (font \"sans\" (wght 1000)))").contains("runs from 100 to 900"));
    assert!(e("(def f (font \"comic\"))").contains("no font 'comic'"));
}

#[test]
fn a_text_param_is_its_value_and_refuses_what_it_does_not_allow() {
    let load = |src: &str, o: &[(String, Given)]| Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), o);
    let src = "(param name \"Ada\" (chars \"AdaGrcen \") (len 1 8))\n(def n (len (char-codes name)))\n";
    let d = load(src, &[]).unwrap();
    assert!(matches!(d.def("name"), Some(Value::Str(s)) if s == "Ada") && d.texts[0].value == "Ada");
    let d = load(src, &[("name".into(), Given::text("Grace"))]).unwrap();
    assert!(matches!(d.def("name"), Some(Value::Str(s)) if s == "Grace"));
    let e = load(src, &[("name".into(), Given::text("Grace!"))]).err().unwrap().message;
    assert!(e.contains("holds \"!\" (character 6), which name does not allow"), "{e}");
    let e = load(src, &[("name".into(), Given::text("GraceGrace"))]).err().unwrap().message;
    assert!(e.contains("holds 10 characters; name holds 1 to 8"), "{e}");
    let e = load("(param name \"Ada!\" (chars \"Ad\"))", &[]).err().unwrap().message;
    assert!(e.contains("name's default"), "{e}");
    let e = load("(param name \"Ada\" (tol 0 1))", &[]).err().unwrap().message;
    assert!(e.contains("a text param has no tolerance"), "{e}");
}

#[test]
fn an_advance_is_the_fonts_at_its_axes() {
    // fontTools' instancer at (wght 700) (wdth 75): i and ï 251 units, H
    // 554 (at the default axes 215, 187 and 699). An advance changed by
    // the axes twice set ï 28 units wider than i ("nai ve").
    let src = "(def f (font \"rounded\" (wght 700) (wdth 75)))\n\
               (def i (get (glyph-metrics f 105) \"advance\"))\n(def iu (get (glyph-metrics f 239) \"advance\"))\n(def h (get (glyph-metrics f 72) \"advance\"))\n";
    let d = Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), &[]).unwrap();
    let num = |n: &str| match d.def(n) {
        Some(Value::Field(f, _)) => d.g.as_const(*f).unwrap(),
        v => panic!("{n}: {v:?}"),
    };
    assert_eq!([num("i"), num("iu"), num("h")], [0.251, 0.251, 0.554]);
}

#[test]
fn a_glyphs_box_holds_its_outline() {
    // Each glyph's box (glyph-metrics) against its glyph: its edges, walked
    // in steps of a thousandth, inside it and reaching to within 1e-6 em of
    // each side, and its field outside (not below 0) on each side, at 2001
    // points along it and where the walk came nearest it. A box of samples,
    // 16 to an arc, left an e's arc 1.4e-6 em outside it (Lexend); one of
    // the arcs' extremes in floats, drawn apart from the polygon, left
    // Lexend Black's G inside the field 2e-17 em past its right side.
    let mut d = Doc::load_str_in(&CORE, "", "t", std::path::Path::new("."), &[]).unwrap();
    let ctx = ink_core::eval::Ctx::new(0.0, Vec::new());
    for (font, axes) in [("sans", vec![]), ("sans", vec![("wght".to_string(), 900.0)]), ("script", vec![])] {
        let s = Spec { key: font.to_string(), axes };
        for c in "AGOQRSaegos@&8".chars() {
            let code = gid(&mut d, &s, c);
            let g = d.glyph_outline(&s, code).unwrap();
            let Value::Data(m) = d.glyph_metrics(&s, code, 0.0).unwrap() else { panic!() };
            let b: Vec<f64> = m.value()["box"].as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect();
            let Value::Field(f, _) = d.glyph_2d(&s, code).unwrap() else { panic!() };
            let tape = Tape::new(&d.g, f);
            // The walk's extremes per side (x low, y low, x high, y high),
            // each with where it is.
            let mut ext = [(f64::MAX, [0.0; 2]), (f64::MAX, [0.0; 2]), (f64::MIN, [0.0; 2]), (f64::MIN, [0.0; 2])];
            for l in &g.loops {
                for (i, p) in l.iter().enumerate() {
                    let q = l[(i + 1) % l.len()];
                    let (a, e) = ([p[0], p[1]], [q[0], q[1]]);
                    let th = 4.0 * p[2].atan();
                    let ch = sub(e, a);
                    let len = norm(ch);
                    for k in 0..=1000 {
                        let u = k as f64 / 1000.0;
                        let v = if p[2] == 0.0 {
                            lerp(a, e, u)
                        } else {
                            let centre = add(mul(add(a, e), 0.5), mul([-ch[1] / len, ch[0] / len], len / (2.0 * (th / 2.0).tan())));
                            let r = sub(a, centre);
                            let (sn, cs) = (th * u).sin_cos();
                            add(centre, [r[0] * cs - r[1] * sn, r[0] * sn + r[1] * cs])
                        };
                        for k in 0..2 {
                            if v[k] < ext[k].0 {
                                ext[k] = (v[k], v);
                            }
                            if v[k] > ext[k + 2].0 {
                                ext[k + 2] = (v[k], v);
                            }
                        }
                    }
                }
            }
            for k in 0..4 {
                let (axis, x) = (k % 2, ext[k].0);
                let (inside, near) = if k < 2 { (b[k] <= x, x - b[k]) } else { (x <= b[k], b[k] - x) };
                assert!(inside && near < 1e-6, "{font} {c}: box {b:?}, outline reaches {x} on side {k}");
                // On the side: the point beside the extreme, and 2001 along.
                let other = 1 - axis;
                let along = (0..=2000).map(|i| b[other] + (b[other + 2] - b[other]) * i as f64 / 2000.0).chain([ext[k].1[other]]);
                for t in along {
                    let mut pt = [0.0; 3];
                    (pt[axis], pt[other]) = (b[k], t);
                    let v = ink_core::eval::eval_point(&tape, &ctx, pt);
                    assert!(v >= 0.0, "{font} {:?} {c}: box {b:?}, field {v} at {pt:?} on its side", s.axes);
                }
            }
        }
    }
}

#[test]
fn a_case_not_allowed_is_recorded_why_not_refused() {
    // A declared case (Given::case, by its place) its param does not allow
    // at the params given is recorded as the param's refusal, typed, and
    // the load stops there, before anything reads the text; the same text
    // given from outside is refused.
    let src = "(param maxlen 1)\n(param name \"A\" (chars \"AB\") (len 1 maxlen) (cases \"AA\" \"C\"))\n(def first (nth (char-codes name) 5))\n";
    let load = |o: &[(String, Given)]| Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), o);
    let d = load(&[("name".into(), Given::case(0))]).unwrap();
    assert!(matches!(d.refused.as_slice(), [(p, t, Refusal::Length { n: 2, lo: 1, hi: 1 })] if p.name == "name" && t == "AA"), "{:?}", d.refused);
    let d = load(&[("name".into(), Given::case(1))]).unwrap();
    assert!(matches!(d.refused.as_slice(), [(_, _, Refusal::Char { at: 0, c: 'C' })]), "{:?}", d.refused);
    assert!(load(&[("name".into(), Given::text("AA"))]).is_err());
}

#[test]
fn a_drawn_box_holds_the_glyph_built() {
    // lettering-max reads each glyph's box from its contours as drawn,
    // not from the polygon built of the region they fill: that box holds
    // the built one, for every character from space to ÿ in each vendored
    // font, at its axes' ends.
    let mut d = Doc::load_str_in(&CORE, "", "t", std::path::Path::new("."), &[]).unwrap();
    let fonts = [
        ("sans", vec![("wght", 100.0)]),
        ("sans", vec![("wght", 900.0)]),
        ("serif", vec![("wght", 900.0)]),
        ("rounded", vec![("wght", 700.0), ("wdth", 75.0)]),
        ("script", vec![]),
    ];
    for (font, axes) in fonts {
        let s = Spec { key: font.to_string(), axes: axes.iter().map(|(k, v)| (k.to_string(), *v)).collect() };
        for c in (33u32..0x100).filter_map(char::from_u32) {
            let Some(g) = d.with_face(&s, |f| Ok(f.glyph_index(c))).unwrap() else { continue };
            let m = d.metrics_json(&s, g).unwrap();
            let drawn = d.with_face(&s, |f| Ok(Doc::drawn_box(f, g))).unwrap();
            if m["ink"] == 0 {
                continue;
            }
            let b: Vec<f64> = m["box"].as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect();
            let x = drawn.unwrap_or_else(|| panic!("{font} {c}: ink, and no contours"));
            assert!(x[0] <= b[0] && x[1] <= b[1] && x[2] >= b[2] && x[3] >= b[3], "{font} {axes:?} {c}: drawn {x:?}, built {b:?}");
        }
    }
}

/// Crafted fonts (testdata, from make-fonts.py, make-fonts-2.py and
/// make-fonts-3.py) that a bound on every line of some characters must
/// hold for or refuse: each font, its axes, the characters, how many, and
/// a line of them. Substitutions reaching glyphs no character maps to (a
/// required feature of any tag, feature variations, lookups across
/// stages, contexts calling contexts, splits and ligatures), positioning
/// repeated by contexts, Apple tables, axis mappings (avar 1.1 and 2),
/// advances and deltas past 16 bits, fractional deltas and coordinates,
/// GDEF and HVAR versions, layout formats with 24-bit offsets, a glyph
/// 20 em tall, scripts with their own shapers, deltas and advances past
/// f32's whole numbers. Last, what a refusal says,
/// where lettering-max refuses.
const REPRODUCERS: &[(&str, &str, &str, usize, &str, Option<&str>)] = &[
    ("arabic-context-stages", "", "\u{628}\u{62a}", 2, "\u{628}\u{62a}", Some("its own shaping rules")),
    ("arabic-context-stages2", "", "\u{628}\u{62a}", 2, "\u{628}\u{62a}", Some("its own shaping rules")),
    ("arabic-ligature-stages", "", "\u{628}\u{62a}", 2, "\u{628}\u{62a}", Some("its own shaping rules")),
    ("arabic-stages", "", "\u{628}", 1, "\u{628}", Some("its own shaping rules")),
    ("hangul", "", "\u{1100}\u{1161}", 2, "\u{1100}\u{1161}", Some("its own shaping rules")),
    ("required", "", "ab", 2, "ab", None),
    ("required-context", "", "ab", 2, "ab", None),
    ("stages", "", "ab", 2, "ab", None),
    ("variation", "(wght 900)", "a", 1, "a", None),
    ("repeated", "", "a", 1, "a", Some("positions glyphs by context")),
    ("avar2", "(wght 650)", "a", 1, "a", None),
    ("deep-context", "", "ab", 8, "aaaaaaab", Some("positions glyphs by context")),
    ("fractional", "(wght 650)", "a", 1, "a", None),
    ("gpos-variation", "(wght 650)", "a", 1, "a", None),
    ("large-glyph", "", "a", 1, "a", None),
    ("mort", "", "ab", 2, "ab", Some("has a mort table")),
    ("nested", "", "ab", 8, "aaaaaaab", Some("positions glyphs by context")),
    ("nested-skip", "", "ab", 7, "aaaaaab", Some("positions glyphs by context")),
    ("shared-stages", "", "ab", 2, "ab", None),
    ("split", "", "x", 1, "x", Some("sets as 3 glyphs")),
    ("stages-required", "", "ab", 2, "ab", None),
    ("gdef11", "(wght 900)", "ab", 2, "ab", Some("GDEF table of version 1.1")),
    ("gdef14", "(wght 900)", "ab", 2, "ab", Some("GDEF table of version 1.4")),
    ("gpos-long", "(wght 900)", "ab", 4, "abab", None),
    ("gpos-long-avar11", "(wght 650)", "ab", 4, "abab", None),
    ("gpos-normalization", "(wght 650.02)", "ab", 2, "ab", None),
    ("gpos-scalar", "(wght 650)", "a", 1, "a", None),
    ("gpos-small", "(wght 900)", "ab", 4, "abab", None),
    ("gpos-small-avar11", "(wght 650)", "ab", 4, "ab", None),
    ("gpos-wide", "", "ab", 2, "ab", Some("type 1, format 3, not read here")),
    ("gpos-wide-context", "", "ab", 2, "ab", Some("positions glyphs by context")),
    ("gsub-wide", "", "ab", 2, "ab", Some("type 2, format 2, not read here")),
    ("hvar11", "(wght 900)", "ab", 2, "ab", None),
    ("hvar-negative", "(wght 900)", "ab", 4, "abab", None),
    ("hvar-normalization", "(wght 650.02)", "ab", 2, "ab", None),
    ("hvar-overflow", "(wght 900)", "ab", 4, "abab", None),
    ("gpos-float-16777217", "(wght 900)", "a", 1, "a", None),
    ("gpos-float-16777219", "(wght 900)", "a", 1, "a", None),
    ("gpos-float--16777217", "(wght 900)", "a", 1, "a", None),
    ("gpos-float--16777219", "(wght 900)", "a", 1, "a", None),
    ("gpos-float-33554435", "(wght 900)", "a", 1, "a", None),
    ("gpos-float--33554435", "(wght 900)", "a", 1, "a", None),
    ("hvar-float-67108865", "(wght 650)", "ab", 4, "abab", None),
    ("hvar-float-67108869", "(wght 650)", "ab", 4, "abab", None),
    ("hvar-float-67108873", "(wght 650)", "ab", 4, "abab", None),
    ("hvar-float-134217739", "(wght 650)", "ab", 4, "abab", None),
];

#[test]
fn every_crafted_font_is_bounded_or_refused() {
    // Each line's box within lettering-max's bound (cap height 1, so the
    // baseline is at -1/2), or lettering-max refused, saying why and to
    // name the texts as cases. The line itself is set either way.
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/font/testdata");
    for &(f, axes, chars, n, text, refused) in REPRODUCERS {
        let font = format!("(import \"std\")\n(def F (font \"{f}.ttf\" {axes}))\n");
        let line = format!("{font}(def b (lettering-box (lettering-styled \"{text}\" F 1 0 1 0)))\n");
        let d = Doc::load_str_in(&CORE, &line, "t", &dir, &[]).unwrap_or_else(|e| panic!("{f}: the line is not set: {}", e.message));
        let v = |d: &Doc, n: &str| -> Vec<f64> {
            let Some(Value::Vec(v)) = d.def(n) else { panic!("{n}") };
            v.iter().map(|x| if let Value::Field(g, _) = x { d.g.as_const(*g).unwrap() } else { panic!() }).collect()
        };
        let b = v(&d, "b");
        match Doc::load_str_in(&CORE, &format!("{font}(def m (lettering-max \"{chars}\" {n} F 1 0))\n"), "t", &dir, &[]) {
            Ok(dm) => {
                assert!(refused.is_none(), "{f}: bounded, not refused ({refused:?})");
                let m = v(&dm, "m");
                let ok = b[2] - b[0] <= m[0] && b[1] + 0.5 >= m[1] - 1e-12 && b[3] + 0.5 <= m[2] + 1e-12;
                assert!(ok, "{f} {text:?}: box {b:?} outside the bound {m:?}");
            }
            Err(e) => {
                let why = refused.unwrap_or_else(|| panic!("{f}: refused: {}", e.message));
                assert!(e.message.contains(why) && e.message.contains("name the texts to check as cases"), "{f}: {}", e.message);
            }
        }
    }
}

/// The instances the outline snapshot covers: each vendored font at its
/// axes' corners and between them.
const SNAPSHOT: &[(&str, &[(&str, f64)])] = &[
    ("sans", &[("wght", 100.0)]),
    ("sans", &[("wght", 400.0)]),
    ("sans", &[("wght", 650.0)]),
    ("sans", &[("wght", 900.0)]),
    ("serif", &[("wght", 100.0)]),
    ("serif", &[("wght", 400.0)]),
    ("serif", &[("wght", 650.0)]),
    ("serif", &[("wght", 900.0)]),
    ("rounded", &[("wght", 300.0), ("wdth", 75.0)]),
    ("rounded", &[("wght", 300.0), ("wdth", 100.0)]),
    ("rounded", &[("wght", 300.0), ("wdth", 125.0)]),
    ("rounded", &[("wght", 700.0), ("wdth", 75.0)]),
    ("rounded", &[("wght", 700.0), ("wdth", 125.0)]),
    ("rounded", &[("wght", 550.0), ("wdth", 90.0)]),
    ("script", &[]),
];

/// FNV-1a over an outline's contours: their count, each one's segments,
/// each segment's points' exact bits.
fn outline_hash(o: &Outline) -> u32 {
    let mut h: u64 = 0xcbf29ce484222325;
    let mut eat = |v: u64| {
        for b in v.to_le_bytes() {
            h = (h ^ b as u64).wrapping_mul(0x100000001b3);
        }
    };
    eat(o.contours.len() as u64);
    for c in &o.contours {
        eat(c.len() as u64);
        for seg in c {
            eat(seg.len() as u64);
            for p in seg {
                eat(p[0].to_bits());
                eat(p[1].to_bits());
            }
        }
    }
    (h ^ (h >> 32)) as u32
}

/// Each snapshot instance: its line metrics as (font ...) gives them, and
/// every glyph's outline hash, one line each.
fn snapshot_lines() -> Vec<String> {
    let mut out = Vec::new();
    for (font, axes) in SNAPSHOT {
        let a: Vec<String> = axes.iter().map(|(k, v)| format!("({k} {v})")).collect();
        let d = Doc::load_str_in(&CORE, &format!("(def f (font \"{font}\" {}))", a.join(" ")), "t", std::path::Path::new("."), &[]).unwrap();
        let Some(Value::Data(f)) = d.def("f") else { panic!() };
        let m: Vec<String> = ["ascender", "descender", "line-gap", "cap-height", "x-height"].iter().map(|k| format!("{}", f.value()[k])).collect();
        let s = Spec { key: font.to_string(), axes: axes.iter().map(|(k, v)| (k.to_string(), *v)).collect() };
        let mut d = d;
        let hashes: Vec<String> = (0..glyph_count(&mut d, &s)).map(|g| format!("{:08x}", outline_hash(&drawn_outline(&mut d, &s, g)))).collect();
        let name = format!("{font}{}", axes.iter().map(|(k, v)| format!(" {k}={v}")).collect::<String>());
        out.push(format!("{name}|{}|{}", m.join(" "), hashes.join(" ")));
    }
    out
}

fn glyph_count(d: &mut Doc, s: &Spec) -> u16 {
    d.with_face(s, |f| Ok(f.number_of_glyphs())).unwrap()
}

fn drawn_outline(d: &mut Doc, s: &Spec, g: u16) -> Outline {
    d.with_face(s, |f| {
        let mut o = Outline::default();
        f.outline_glyph(g, &mut o);
        Ok(o)
    })
    .unwrap()
}

#[test]
fn every_vendored_outline_is_as_it_was() {
    // Every glyph of each vendored font, against how ink drew it when it
    // read fonts with ttf-parser 0.25.1 (testdata/outlines.txt): the same
    // line metrics everywhere, and the same points (exact bits) at each
    // font's default and the instances listed. Elsewhere the variations'
    // deltas are summed in another order, so a point may move by a few
    // thousandths of a unit (outline_points_are_within_rounding).
    let exact = ["sans wght=100", "sans wght=400", "sans wght=900", "serif wght=100", "rounded wght=300 wdth=100", "script"];
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/font/testdata/outlines.txt");
    let want: Vec<String> = std::fs::read_to_string(path).unwrap().lines().filter(|l| !l.starts_with('#')).map(String::from).collect();
    let got = snapshot_lines();
    assert_eq!(want.len(), got.len());
    for (w, g) in want.iter().zip(&got) {
        let (wp, gp): (Vec<&str>, Vec<&str>) = (w.split('|').collect(), g.split('|').collect());
        assert_eq!((wp[0], wp[1]), (gp[0], gp[1]), "{}: line metrics", wp[0]);
        if exact.contains(&wp[0]) {
            let (wh, gh): (Vec<&str>, Vec<&str>) = (wp[2].split(' ').collect(), gp[2].split(' ').collect());
            let differ: Vec<usize> = (0..wh.len().max(gh.len())).filter(|&i| wh.get(i) != gh.get(i)).collect();
            assert!(differ.is_empty(), "{}: {} of {} glyphs drawn otherwise, first {:?}", wp[0], differ.len(), wh.len(), &differ[..differ.len().min(8)]);
        }
    }
}

#[test]
fn outline_points_are_within_rounding() {
    // Each digit and letter of each vendored font, at each instance where
    // some glyph is drawn otherwise than ttf-parser drew it: the same
    // contours and segments, each point within 0.003 units of its
    // (testdata/outline-points.txt, kept to 1/1000 unit: 0.0035 in all).
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/font/testdata/outline-points.txt");
    let mut d = Doc::load_str_in(&CORE, "", "t", std::path::Path::new("."), &[]).unwrap();
    let (mut worst, mut n) = (0.0f64, 0);
    for line in std::fs::read_to_string(path).unwrap().lines().filter(|l| !l.starts_with('#')) {
        // "sans wght=650 A=<outline>": the font, its axes, the character.
        let mut toks = line.splitn(2, ' ');
        let font = toks.next().unwrap();
        let mut rest = toks.next().unwrap();
        let mut axes = Vec::new();
        while let Some((k, r)) = rest.split_once('=').filter(|(k, _)| *k == "wght" || *k == "wdth") {
            let (v, r) = r.split_once(' ').unwrap();
            axes.push((k.to_string(), v.parse::<f64>().unwrap()));
            rest = r;
        }
        let (c, outline) = (rest.chars().next().unwrap(), &rest[2..]);
        let s = Spec { key: font.to_string(), axes };
        let g = gid(&mut d, &s, c);
        let o = drawn_outline(&mut d, &s, g);
        let contours: Vec<&str> = outline.split('|').collect();
        assert_eq!(o.contours.len(), contours.len(), "{font} {:?} {c}: contours", s.axes);
        for (got, want) in o.contours.iter().zip(contours) {
            let tok: Vec<&str> = want.split(' ').collect();
            let num = |i: usize| tok[i].parse::<f64>().unwrap();
            let mut at = 2;
            let mut close = |a: [f64; 2], i: usize| {
                let e = (a[0] - num(i)).abs().max((a[1] - num(i + 1)).abs());
                worst = worst.max(e);
                e <= 0.0035
            };
            assert!(close(got[0][0], 0), "{font} {:?} {c}: start {:?} against {} {}", s.axes, got[0][0], tok[0], tok[1]);
            assert_eq!(got.len(), tok.iter().filter(|t| matches!(**t, "L" | "Q" | "C")).count(), "{font} {:?} {c}: segments", s.axes);
            for seg in got {
                let k = match tok[at] {
                    "L" => 1,
                    "Q" => 2,
                    _ => 3,
                };
                assert_eq!(seg.len() - 1, k, "{font} {:?} {c}: a segment's kind", s.axes);
                for (j, p) in seg[1..].iter().enumerate() {
                    assert!(close(*p, at + 1 + 2 * j), "{font} {:?} {c}: {p:?} against {} {}", s.axes, tok[at + 1 + 2 * j], tok[at + 2 + 2 * j]);
                }
                at += 1 + 2 * k;
            }
            n += 1;
        }
    }
    assert!(n > 800, "{n} contours read");
    eprintln!("outline points within {worst} units of ttf-parser's over {n} contours");
}
