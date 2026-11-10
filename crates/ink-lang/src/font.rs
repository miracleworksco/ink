//! Fonts: a TrueType or OpenType font's glyphs as profiles. `(font NAME
//! [(AXIS v) ...])` reads a font (one std vendors, or a file) at a point of
//! its variation axes; `(glyph-run FONT TEXT)` is the text shaped (shape.rs):
//! its lines' glyphs, each its outline (a field of x and y in ems, made of
//! `polygon`s) placed, with its metrics; `(glyph-2d FONT CODE)` is a
//! character's glyph's outline; `(glyph-metrics FONT CODE [NEXT])` its
//! advance, its kerning before NEXT, its box and a point in each of its
//! holes; `(glyph-reach FONT CHARS)` bounds any line of the characters
//! CHARS. std lays text out from these (std's lettering section).
//!
//! An outline's quadratic and cubic curves become arcs (biarcs: two arcs
//! meeting tangent, fitted to each curve's ends and their tangents, the
//! curve halved until they lie within `TOL` of it), so a glyph is a
//! `polygon` of lines and arcs: smooth, tangent where the font is, and
//! exactly the polygon ink proves things of. The font file itself is read
//! as it ships, never converted.

mod face;
mod region;
mod shape;

use crate::*;
use face::{Face, Pen};

/// The fonts std vendors (std/fonts, each with its licence): the name a
/// model asks for, the family, and the file.
pub(crate) const VENDORED: &[(&str, &str, &[u8])] = &[
    ("sans", "Lexend", include_bytes!("../../../std/fonts/Lexend-wght.ttf")),
    ("rounded", "Fredoka", include_bytes!("../../../std/fonts/Fredoka-wdth-wght.ttf")),
    ("serif", "Bitter", include_bytes!("../../../std/fonts/Bitter-wght.ttf")),
    ("script", "Pacifico", include_bytes!("../../../std/fonts/Pacifico-Regular.ttf")),
];

/// How far an outline's arcs may lie from the font's curves, in ems (6 µm
/// on capitals 10 mm tall): far below what a printer or a cutter holds,
/// and each halving of it doubles the arcs a letter takes.
const TOL: f64 = 4e-4;

/// What a document has read of fonts: files by path, and glyphs built.
#[derive(Default)]
pub(crate) struct Fonts {
    files: HashMap<String, Rc<[u8]>>,
    glyphs: HashMap<String, NodeId>,
    outlines: HashMap<String, Rc<Glyph>>,
    shapers: HashMap<String, harfrust::Font>,
}

type P = [f64; 2];

/// A font as `(font ...)` gives it: its key (a vendored name, or a file's
/// path) and the axes set.
struct Spec {
    key: String,
    axes: Vec<(String, f64)>,
}

/// An outline as the font draws it, in font units: closed contours of
/// curves, each its control points from its start.
#[derive(Default)]
struct Outline {
    contours: Vec<Vec<Vec<P>>>,
    at: P,
}

impl Pen for Outline {
    fn move_to(&mut self, x: f32, y: f32) {
        self.contours.push(Vec::new());
        self.at = [x as f64, y as f64];
    }
    fn line_to(&mut self, x: f32, y: f32) {
        self.seg(&[[x as f64, y as f64]]);
    }
    fn quad_to(&mut self, x1: f32, y1: f32, x: f32, y: f32) {
        self.seg(&[[x1 as f64, y1 as f64], [x as f64, y as f64]]);
    }
    fn curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, x: f32, y: f32) {
        self.seg(&[[x1 as f64, y1 as f64], [x2 as f64, y2 as f64], [x as f64, y as f64]]);
    }
    fn close(&mut self) {
        let start = self.contours.last().and_then(|c| c.first()).map(|s| s[0]);
        if let Some(s) = start
            && s != self.at
        {
            self.seg(&[s]);
        }
    }
}

impl Outline {
    fn seg(&mut self, pts: &[P]) {
        let mut c = vec![self.at];
        c.extend_from_slice(pts);
        self.at = *pts.last().expect("a point");
        if let Some(k) = self.contours.last_mut() {
            k.push(c);
        }
    }
}

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1]]
}
fn add(a: P, b: P) -> P {
    [a[0] + b[0], a[1] + b[1]]
}
fn mul(a: P, s: f64) -> P {
    [a[0] * s, a[1] * s]
}
fn dot(a: P, b: P) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn cross(a: P, b: P) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn norm(a: P) -> f64 {
    a[0].hypot(a[1])
}
fn lerp(a: P, b: P, t: f64) -> P {
    add(a, mul(sub(b, a), t))
}

/// A Bezier curve's point at t (de Casteljau).
fn bez(c: &[P], t: f64) -> P {
    let mut v = c.to_vec();
    while v.len() > 1 {
        v = v.windows(2).map(|w| lerp(w[0], w[1], t)).collect();
    }
    v[0]
}

/// A Bezier curve halved at t = 1/2.
fn halve(c: &[P]) -> (Vec<P>, Vec<P>) {
    let (mut l, mut r) = (vec![c[0]], vec![*c.last().expect("a point")]);
    let mut v = c.to_vec();
    while v.len() > 1 {
        v = v.windows(2).map(|w| lerp(w[0], w[1], 0.5)).collect();
        l.push(v[0]);
        r.push(*v.last().expect("a point"));
    }
    r.reverse();
    (l, r)
}

/// The distance from p to the edge from a to b of bulge `bu` (a line where
/// 0, else an arc of sweep 4 atan bu, positive counter-clockwise).
fn edge_distance(p: P, a: P, b: P, bu: f64) -> f64 {
    let c = sub(b, a);
    let l = norm(c);
    if bu.abs() < 1e-12 || l == 0.0 {
        let t = if l == 0.0 { 0.0 } else { (dot(sub(p, a), c) / (l * l)).clamp(0.0, 1.0) };
        return norm(sub(p, add(a, mul(c, t))));
    }
    let th = 4.0 * bu.atan();
    let n = [-c[1] / l, c[0] / l];
    let centre = add(mul(add(a, b), 0.5), mul(n, l / (2.0 * (th / 2.0).tan())));
    let r = norm(sub(a, centre));
    let (u, v) = (sub(a, centre), sub(p, centre));
    let mut phi = cross(u, v).atan2(dot(u, v));
    if th > 0.0 && phi < 0.0 {
        phi += std::f64::consts::TAU;
    } else if th < 0.0 && phi > 0.0 {
        phi -= std::f64::consts::TAU;
    }
    if phi.abs() <= th.abs() { (norm(v) - r).abs() } else { norm(sub(p, a)).min(norm(sub(p, b))) }
}

/// The point halfway along the edge from a to b of bulge bu.
fn edge_mid(a: P, b: P, bu: f64) -> P {
    let c = sub(b, a);
    let m = mul(add(a, b), 0.5);
    // The sagitta: tan(sweep / 4) of half the chord, right of it.
    add(m, mul([c[1], -c[0]], bu / 2.0))
}

/// The bulge of the arc from a to b leaving a along unit tangent t.
fn bulge_from_start(a: P, b: P, t: P) -> f64 {
    let c = sub(b, a);
    (cross(t, c).atan2(dot(t, c)) / 2.0).tan()
}

/// The bulge of the arc from a to b arriving at b along unit tangent t.
fn bulge_to_end(a: P, b: P, t: P) -> f64 {
    let c = sub(b, a);
    (cross(c, t).atan2(dot(c, t)) / 2.0).tan()
}

/// Curve `c` (its control points, from its start) as edges within `tol`
/// of it: each its start and bulge, pushed to `out` (its end is the next
/// edge's start).
fn fit(c: &[P], tol: f64, depth: u32, out: &mut Vec<[f64; 3]>) {
    let (p0, p1) = (c[0], *c.last().expect("a point"));
    let chord = sub(p1, p0);
    let l = norm(chord);
    // Straight: every control point within tol of the chord.
    let off = |q: P| if l == 0.0 { norm(sub(q, p0)) } else { cross(chord, sub(q, p0)).abs() / l };
    if c.iter().all(|&q| off(q) <= tol / 4.0) && c.iter().all(|&q| l == 0.0 || (-tol..=l * l + tol).contains(&dot(sub(q, p0), chord))) {
        if l > 0.0 {
            out.push([p0[0], p0[1], 0.0]);
        }
        return;
    }
    let split = |out: &mut Vec<[f64; 3]>| {
        let (a, b) = halve(c);
        fit(&a, tol, depth + 1, out);
        fit(&b, tol, depth + 1, out);
    };
    // The tangents at the ends, from the first and last control points
    // apart from them.
    let t0 = c.iter().skip(1).map(|&q| sub(q, p0)).find(|d| norm(*d) > 1e-12 * (1.0 + l));
    let t1 = c.iter().rev().skip(1).map(|&q| sub(p1, q)).find(|d| norm(*d) > 1e-12 * (1.0 + l));
    let (Some(t0), Some(t1)) = (t0, t1) else {
        if l > 0.0 {
            out.push([p0[0], p0[1], 0.0]);
        }
        return;
    };
    let (t0, t1) = (mul(t0, 1.0 / norm(t0)), mul(t1, 1.0 / norm(t1)));
    // The biarc whose arcs reach equally far along the end tangents
    // (d1 = d2 = d): |q1 - q2| = 2d, q1 = p0 + d t0, q2 = p1 - d t1.
    let t = add(t0, t1);
    let den = 2.0 * (1.0 - dot(t0, t1));
    let vt = dot(chord, t);
    let d = if den.abs() > 1e-12 {
        (-vt + (vt * vt + den * dot(chord, chord)).sqrt()) / den
    } else if dot(chord, t1).abs() > 1e-12 {
        dot(chord, chord) / (4.0 * dot(chord, t1))
    } else {
        f64::NAN
    };
    if !(d.is_finite() && d > 0.0) || depth > 18 {
        if depth > 18 {
            out.push([p0[0], p0[1], 0.0]);
            return;
        }
        return split(out);
    }
    let pm = mul(add(add(p0, mul(t0, d)), sub(p1, mul(t1, d))), 0.5);
    let (b0, b1) = (bulge_from_start(p0, pm, t0), bulge_to_end(pm, p1, t1));
    // Within tol: the curve's samples of the biarc, and the biarc's of the
    // curve (its polyline of 64 pieces).
    let near_arcs = |q: P| edge_distance(q, p0, pm, b0).min(edge_distance(q, pm, p1, b1));
    let poly: Vec<P> = (0..=64).map(|k| bez(c, k as f64 / 64.0)).collect();
    let near_curve = |q: P| poly.windows(2).map(|w| edge_distance(q, w[0], w[1], 0.0)).fold(f64::MAX, f64::min);
    let ok = (1..16).all(|k| near_arcs(bez(c, k as f64 / 16.0)) <= tol)
        && [pm, edge_mid(p0, pm, b0), edge_mid(pm, p1, b1)].into_iter().all(|q| near_curve(q) <= tol)
        && b0.abs() < 1.0
        && b1.abs() < 1.0;
    if !ok {
        return split(out);
    }
    out.push([p0[0], p0[1], b0]);
    out.push([pm[0], pm[1], b1]);
}

/// A loop's signed area (counter-clockwise positive), its arcs' segments
/// included.
fn area(lp: &[[f64; 3]]) -> f64 {
    let mut s = 0.0;
    for (i, p) in lp.iter().enumerate() {
        let q = lp[(i + 1) % lp.len()];
        s += cross([p[0], p[1]], [q[0], q[1]]) / 2.0;
        if p[2] != 0.0 {
            let th = 4.0 * p[2].atan();
            let l = norm(sub([q[0], q[1]], [p[0], p[1]]));
            let r = l / (2.0 * (th / 2.0).sin()).abs();
            s += th.signum() * r * r / 2.0 * (th.abs() - th.abs().sin());
        }
    }
    s
}

/// Whether p lies inside the loop's chords (even-odd).
fn inside(p: P, lp: &[[f64; 3]]) -> bool {
    let mut odd = false;
    for (i, a) in lp.iter().enumerate() {
        let b = lp[(i + 1) % lp.len()];
        if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) {
            odd = !odd;
        }
    }
    odd
}

/// A point well inside a loop: the middle of the widest span of its
/// chords across the line halfway up it.
fn interior(lp: &[[f64; 3]]) -> P {
    let (lo, hi) = lp.iter().fold((f64::MAX, f64::MIN), |(l, h), p| (l.min(p[1]), h.max(p[1])));
    let y = (lo + hi) / 2.0;
    let mut xs: Vec<f64> = Vec::new();
    for (i, a) in lp.iter().enumerate() {
        let b = lp[(i + 1) % lp.len()];
        if (a[1] > y) != (b[1] > y) {
            xs.push(a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]));
        }
    }
    xs.sort_by(f64::total_cmp);
    let best = xs.chunks(2).filter(|c| c.len() == 2).max_by(|a, b| (a[1] - a[0]).total_cmp(&(b[1] - b[0])));
    best.map_or([lp[0][0], lp[0][1]], |c| [(c[0] + c[1]) / 2.0, y])
}

/// The box of the edge from a to b of bulge bu (a line where 0, else an
/// arc of sweep 4 atan bu): its ends, and the arc's extremes along x and y
/// that its sweep passes.
fn edge_box(a: P, b: P, bu: f64) -> (P, P) {
    let (mut lo, mut hi) = ([a[0].min(b[0]), a[1].min(b[1])], [a[0].max(b[0]), a[1].max(b[1])]);
    let c = sub(b, a);
    let l = norm(c);
    if bu.abs() < 1e-12 || l == 0.0 {
        return (lo, hi);
    }
    let th = 4.0 * bu.atan();
    let n = [-c[1] / l, c[0] / l];
    let centre = add(mul(add(a, b), 0.5), mul(n, l / (2.0 * (th / 2.0).tan())));
    let r = norm(sub(a, centre));
    let u = sub(a, centre);
    for (k, d) in [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]].into_iter().enumerate() {
        // The angle from a to direction d, the way the arc turns.
        let mut phi = cross(u, d).atan2(dot(u, d));
        if th > 0.0 && phi < 0.0 {
            phi += std::f64::consts::TAU;
        } else if th < 0.0 && phi > 0.0 {
            phi -= std::f64::consts::TAU;
        }
        if phi.abs() <= th.abs() {
            let e = add(centre, mul(d, r));
            if k < 2 {
                hi[k] = hi[k].max(e[k]);
            } else {
                lo[k - 2] = lo[k - 2].min(e[k - 2]);
            }
        }
    }
    (lo, hi)
}

/// A glyph's outline in ems: loops of edges (start, bulge), as the font
/// winds them.
struct Glyph {
    loops: Vec<Vec<[f64; 3]>>,
    advance: f64,
}

impl Spec {
    fn of(d: &Data) -> Result<Spec, String> {
        let v = d.value();
        let key = v.get("font").and_then(|k| k.as_str()).ok_or("a font is what (font NAME) gives: data with its \"font\"")?;
        let axes = v.get("axes").and_then(|a| a.as_object()).map(|o| o.iter().filter_map(|(k, x)| Some((k.clone(), x.as_f64()?))).collect()).unwrap_or_default();
        Ok(Spec { key: key.to_string(), axes })
    }

    fn id(&self) -> String {
        let axes: Vec<String> = self.axes.iter().map(|(k, v)| format!("{k}={v}")).collect();
        format!("{}|{}", self.key, axes.join(","))
    }
}

/// The 4-letter tag of an axis name.
fn tag(s: &str) -> Option<harfrust::Tag> {
    let b = s.as_bytes();
    (b.len() <= 4 && !b.is_empty()).then(|| {
        let mut t = [b' '; 4];
        t[..b.len()].copy_from_slice(b);
        harfrust::Tag::new(&t)
    })
}

fn char_name(code: u32) -> String {
    match char::from_u32(code) {
        Some(c) if !c.is_control() => format!("'{c}' (U+{code:04X})"),
        _ => format!("U+{code:04X}"),
    }
}

impl Doc {
    /// The bytes of font `key`: a vendored font's, else the file's, read.
    fn font_bytes(&self, key: &str) -> Result<Rc<[u8]>, String> {
        if let Some((_, _, b)) = VENDORED.iter().find(|v| v.0 == key) {
            return Ok(Rc::from(*b));
        }
        self.fonts.files.get(key).cloned().ok_or_else(|| format!("no font '{key}' was read: make one with (font NAME)"))
    }

    /// Font `s` as ink and the shaper read it (face.rs), its axes set: one
    /// instance, read once.
    fn shaper(&mut self, s: &Spec) -> Result<harfrust::Font, String> {
        if let Some(f) = self.fonts.shapers.get(&s.id()) {
            return Ok(f.clone());
        }
        let font = match VENDORED.iter().find(|v| v.0 == s.key) {
            Some((_, _, b)) => harfrust::Font::new(*b, 0),
            None => harfrust::Font::new(self.font_bytes(&s.key)?.to_vec(), 0),
        };
        let font = font.ok_or_else(|| format!("{} is not a font ink reads", s.key))?;
        for (k, _) in &s.axes {
            tag(k).filter(|t| font.axes().get_by_tag(*t).is_some()).ok_or_else(|| format!("font {} has no axis '{k}'", s.key))?;
        }
        let font = font.instance_builder().variations(s.axes.iter().map(|(k, v)| (k.as_str(), *v as f32))).build();
        self.fonts.shapers.insert(s.id(), font.clone());
        Ok(font)
    }

    /// Font `s`'s face, its axes set.
    fn with_face<T>(&mut self, s: &Spec, f: impl FnOnce(&Face) -> Result<T, String>) -> Result<T, String> {
        let font = self.shaper(s)?;
        let bytes = self.font_bytes(&s.key)?;
        let face = Face::new(&bytes, &font).map_err(|m| format!("{} is {m}", s.key))?;
        f(&face)
    }

    /// Character `code`'s glyph in font `s`.
    fn glyph_id(face: &Face, s: &Spec, code: f64) -> Result<u16, String> {
        let c = (code >= 0.0 && code.fract() == 0.0 && code <= u32::MAX as f64)
            .then_some(code as u32)
            .ok_or_else(|| format!("a character is its code, a whole number (char-codes gives them), not {code}"))?;
        char::from_u32(c).and_then(|ch| face.glyph_index(ch)).ok_or_else(|| format!("font {} has no glyph for {}", s.key, char_name(c)))
    }

    /// How far character `b`'s glyph sits closer to `a`'s than `a`'s
    /// advance, in ems, as the shaper sets the pair: its kerning.
    fn pair_kern(&mut self, s: &Spec, a: char, b: char) -> Result<f64, String> {
        let font = self.shaper(s)?;
        let (_, alone, _) = shape::shape(&font, &s.key, &a.to_string())?;
        let (pair, _, _) = shape::shape(&font, &s.key, &format!("{a}{b}"))?;
        let pen = pair.iter().find(|g| g.cluster == 1).map_or(alone, |g| g.pen);
        Ok((pen - alone) / font.units_per_em() as f64)
    }

    /// A box round glyph `id`'s ink in ems, without building it: the box
    /// of its contours as drawn (the region they fill lies on them, cut
    /// where they cross), widened for the points the polygon reads as one
    /// and its rounding (it widens its box 1e-7 of its diagonal) by 1e-6
    /// of its diagonal, at least 1e-6 em; None where it has no ink.
    fn drawn_box(face: &Face, id: u16) -> Option<[f64; 4]> {
        let upem = face.units_per_em() as f64;
        let mut o = Outline::default();
        face.outline_glyph(id, &mut o);
        let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
        for c in &o.contours {
            let mut lp = Vec::new();
            for seg in c {
                let seg: Vec<P> = seg.iter().map(|p| mul(*p, 1.0 / upem)).collect();
                fit(&seg, TOL, 0, &mut lp);
            }
            for i in 0..lp.len() {
                let q = lp[(i + 1) % lp.len()];
                let (a, b) = edge_box([lp[i][0], lp[i][1]], [q[0], q[1]], lp[i][2]);
                for k in 0..2 {
                    lo[k] = lo[k].min(a[k]);
                    hi[k] = hi[k].max(b[k]);
                }
            }
        }
        let w = 1e-6 * norm(sub(hi, lo)).max(1.0);
        (lo[0] <= hi[0]).then(|| [lo[0] - w, lo[1] - w, hi[0] + w, hi[1] + w])
    }

    fn outline(face: &Face, id: u16) -> Result<Glyph, String> {
        let upem = face.units_per_em() as f64;
        let mut o = Outline::default();
        face.outline_glyph(id, &mut o);
        let mut drawn = Vec::new();
        for c in &o.contours {
            let mut lp = Vec::new();
            for seg in c {
                let seg: Vec<P> = seg.iter().map(|p| mul(*p, 1.0 / upem)).collect();
                fit(&seg, TOL, 0, &mut lp);
            }
            if lp.len() >= 2 {
                let edges: Vec<region::Edge> =
                    (0..lp.len()).map(|i| region::Edge { a: [lp[i][0], lp[i][1]], b: [lp[(i + 1) % lp.len()][0], lp[(i + 1) % lp.len()][1]], bu: lp[i][2] }).collect();
                drawn.push(edges);
            }
        }
        // What the font fills, as simple loops (outsides counter-clockwise).
        let loops: Vec<Vec<[f64; 3]>> =
            region::filled(&drawn).into_iter().map(|l| l.iter().map(|e| [e.a[0], e.a[1], e.bu]).collect::<Vec<_>>()).filter(|l| area(l).abs() > 1e-12).collect();
        // At the face's axes (its HVAR deltas, else its phantom points).
        let advance = face.advance(id) / upem;
        Ok(Glyph { loops, advance })
    }

    /// `(font NAME [(AXIS v) ...])`, `(glyph-run FONT TEXT)`, `(glyph-2d
    /// FONT CODE)`, `(glyph-metrics FONT CODE [NEXT])` and `(glyph-reach
    /// FONT CHARS)`.
    pub(super) fn eval_font(&mut self, name: &str, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        if name == "font" {
            return self.eval_font_spec(e, args, env);
        }
        let usage = match name {
            "glyph-run" => {
                "usage: (glyph-run FONT TEXT): TEXT shaped in FONT, (font NAME)'s: its lines (split at \"\\n\"), each its glyphs, each (vec outline metrics), placed, in ems"
            }
            "glyph-reach" => "usage: (glyph-reach FONT CHARS): bounds on any line of the characters CHARS set in FONT, in ems",
            "glyph-2d" => "usage: (glyph-2d FONT CODE): character CODE's outline in FONT, (font NAME)'s, in ems, its pen at the origin on the baseline",
            _ => "usage: (glyph-metrics FONT CODE [NEXT]): character CODE's advance, its kerning before NEXT, its box and a point in each hole, in ems",
        };
        let (f, a, next) = match args {
            [f, a] => (f, a, None),
            [f, a, n] if name == "glyph-metrics" => (f, a, Some(n)),
            _ => return Err(self.err_hint(e, "wrong arguments", Some(usage.into()))),
        };
        let Value::Data(d) = self.eval(f, env)? else { return Err(self.err_hint(f, "FONT is a font, as (font NAME) gives", Some(usage.into()))) };
        let s = Spec::of(&d).map_err(|m| self.err_hint(f, m, Some(usage.into())))?;
        if name == "glyph-run" || name == "glyph-reach" {
            let Value::Str(t) = self.eval(a, env)? else { return Err(self.err_hint(a, "a text", Some(usage.into()))) };
            let r = if name == "glyph-run" { self.glyph_run(&s, &t) } else { self.glyph_reach(&s, &t) };
            return r.map_err(|m| self.err(e, m));
        }
        let code = self.const_num(a, env)?;
        let next = next.map(|n| self.const_num(n, env)).transpose()?;
        let ids = self.with_face(&s, |face| Ok((Self::glyph_id(face, &s, code)?, next.map(|n| Self::glyph_id(face, &s, n)).transpose()?)));
        let (id, next_id) = ids.map_err(|m| self.err(e, m))?;
        if name == "glyph-2d" {
            return self.glyph_2d(&s, id).map_err(|m| self.err(e, m));
        }
        // Both have glyphs: char::from_u32 holds for each code.
        let as_char = |c: f64| char::from_u32(c as u32).unwrap_or('\0');
        let kern = match (next, next_id) {
            (Some(n), Some(_)) => self.pair_kern(&s, as_char(code), as_char(n)).map_err(|m| self.err(e, m))?,
            _ => 0.0,
        };
        self.glyph_metrics(&s, id, kern).map_err(|m| self.err(e, m))
    }

    fn eval_font_spec(&mut self, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        let names: Vec<&str> = VENDORED.iter().map(|v| v.0).collect();
        let usage = format!(
            "usage: (font NAME [(AXIS v) ...]): NAME one std vendors ({}) or a .ttf or .otf file (a path relative to this file); \
             AXIS a variation axis the font has, wght its weight (100 thin .. 400 regular .. 900 black)",
            names.join(", ")
        );
        let Some((n, axes)) = args.split_first() else { return Err(self.err_hint(e, "wrong arguments", Some(usage))) };
        let Value::Str(n_) = self.eval(n, env)? else { return Err(self.err_hint(n, "a font's NAME is a text", Some(usage))) };
        let key = if VENDORED.iter().any(|v| v.0 == n_) {
            n_.clone()
        } else if n_.ends_with(".ttf") || n_.ends_with(".otf") || n_.ends_with(".TTF") || n_.ends_with(".OTF") {
            let path = self.dir().join(&n_);
            let key = path.canonicalize().unwrap_or(path.clone()).to_string_lossy().to_string();
            if self.fonts.files.contains_key(&key) {
                self.read_again(&path);
            } else {
                let bytes =
                    std::fs::read(&path).inspect_err(|_| self.files_read.push((path.clone(), None))).map_err(|er| self.err(n, format!("cannot read {}: {er}", shown(&path))))?;
                self.read(&shown(&path), Some(&path), &bytes);
                if harfrust::Font::new(bytes.clone(), 0).is_none() || Face::new(&bytes, &harfrust::Font::new(bytes.clone(), 0).expect("a font")).is_err() {
                    return Err(self.err(n, format!("{} is not a font ink reads", shown(&path))));
                }
                self.fonts.files.insert(key.clone(), Rc::from(bytes));
            }
            key
        } else {
            return Err(self.err_hint(n, format!("no font '{n_}'"), Some(usage)));
        };
        let mut set: Vec<(String, f64)> = Vec::new();
        for a in axes {
            let Some([t, v]) = a.list() else { return Err(self.err_hint(a, "an axis is (AXIS v), (wght 700)", Some(usage))) };
            let Some(t) = t.sym().map(String::from) else { return Err(self.err_hint(a, "an axis is (AXIS v), (wght 700)", Some(usage))) };
            let v = self.const_num(v, env)?;
            set.push((t, v));
        }
        let s = Spec { key: key.clone(), axes: Vec::new() };
        let mut out = self
            .with_face(&s, |face| {
                let upem = face.units_per_em() as f64;
                let mut ranges = serde_json::Map::new();
                let mut axes = serde_json::Map::new();
                for ax in face.font.axes().iter() {
                    let t = ax.tag().to_string().trim().to_string();
                    ranges.insert(t.clone(), serde_json::json!([ax.min_value(), ax.default_value(), ax.max_value()]));
                    axes.insert(t, serde_json::json!(ax.default_value()));
                }
                for (t, v) in &set {
                    let Some(r) = ranges.get(t) else {
                        let have: Vec<&String> = ranges.keys().collect();
                        return Err(format!(
                            "font {key} has no axis '{t}'{}",
                            if have.is_empty() { ": it has none".to_string() } else { format!("; it has {}", have.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ")) }
                        ));
                    };
                    let (lo, hi) = (r[0].as_f64().unwrap_or(0.0), r[2].as_f64().unwrap_or(0.0));
                    if !(lo..=hi).contains(v) {
                        return Err(format!("font {key}'s {t} runs from {lo} to {hi}, not {v}"));
                    }
                    axes.insert(t.clone(), serde_json::json!(v));
                }
                let em = |v: Option<f64>| v.map_or(serde_json::Value::Null, |v| serde_json::json!(v / upem));
                let family =
                    VENDORED.iter().find(|v| v.0 == key).map(|v| v.1.to_string()).or_else(|| {
                        face.font.names().get(harfrust::font::name::NameId::FAMILY_NAME).and_then(|n| n.english_or_first()).map(|l| l.name().chars().collect::<String>())
                    });
                let m = face.font.metrics();
                let line = m.h_line();
                Ok(serde_json::json!({
                    "font": key,
                    "family": family,
                    "axes": axes,
                    "axis-ranges": ranges,
                    "ascender": em(Some(line.map_or(0.0, |l| l.ascender.to_f64()))),
                    "descender": em(Some(line.map_or(0.0, |l| l.descender.to_f64()))),
                    "line-gap": em(Some(line.map_or(0.0, |l| l.line_gap.to_f64()))),
                    "cap-height": em(m.cap_height.map(|v| v.to_f64())),
                    "x-height": em(m.x_height.map(|v| v.to_f64())),
                }))
            })
            .map_err(|m| self.err_hint(e, m, Some(usage.clone())))?;
        // A font without a cap height in its tables: its H's top.
        if out["cap-height"].is_null() || out["cap-height"] == 0.0 {
            let s = Spec::of(&Data::new(out.clone())).map_err(|m| self.err(e, m))?;
            let h =
                self.with_face(&s, |face| Ok(Self::outline(face, face.glyph_index('H').ok_or("no H")?)?.loops.iter().flatten().fold(0.0f64, |m, p| m.max(p[1])))).unwrap_or(0.7);
            out["cap-height"] = serde_json::json!(h);
        }
        Ok(Value::Data(Data::new(out)))
    }

    /// Glyph `glyph`'s outline in font `s`, built once.
    fn glyph_2d(&mut self, s: &Spec, glyph: u16) -> Result<Value, String> {
        let id = format!("{}|{glyph}", s.id());
        if let Some(&f) = self.fonts.glyphs.get(&id) {
            return Ok(Value::Field(f, None));
        }
        let g = self.glyph_outline(s, glyph)?;
        let what = || format!("font {}'s glyph {glyph}", s.key);
        if g.loops.is_empty() {
            let f = self.g.constant(f64::INFINITY);
            self.fonts.glyphs.insert(id, f);
            return Ok(Value::Field(f, None));
        }
        // Outsides run counter-clockwise, holes clockwise (region::filled),
        // none crossing another: a hole is cut from the smallest outside
        // round it, and the outsides joined.
        let areas: Vec<f64> = g.loops.iter().map(|l| area(l)).collect();
        let outer: Vec<bool> = areas.iter().map(|a| *a > 0.0).collect();
        let mut fields: Vec<NodeId> = Vec::new();
        for l in &g.loops {
            fields.push(self.flat_polygon(l, true).map_err(|m| format!("{}: {m}", what()))?);
        }
        let mut cut: Vec<NodeId> = fields.clone();
        for (h, l) in g.loops.iter().enumerate().filter(|(i, _)| !outer[*i]) {
            let p = interior(l);
            let parent = (0..g.loops.len()).filter(|&o| outer[o] && inside(p, &g.loops[o])).min_by(|&a, &b| areas[a].abs().total_cmp(&areas[b].abs()));
            if let Some(o) = parent {
                let n = self.g.neg(fields[h]);
                cut[o] = self.g.max(cut[o], n);
            }
        }
        let f = (0..g.loops.len()).filter(|&i| outer[i]).map(|i| cut[i]).reduce(|a, b| self.g.min(a, b)).expect("an outer loop");
        self.fonts.glyphs.insert(id, f);
        Ok(Value::Field(f, None))
    }

    /// Glyph `glyph`'s outline in font `s`, read once.
    fn glyph_outline(&mut self, s: &Spec, glyph: u16) -> Result<Rc<Glyph>, String> {
        let id = format!("{}|{glyph}", s.id());
        if let Some(g) = self.fonts.outlines.get(&id) {
            return Ok(g.clone());
        }
        let g = Rc::new(self.with_face(s, |face| Self::outline(face, glyph))?);
        self.fonts.outlines.insert(id, g.clone());
        Ok(g)
    }

    /// Glyph `glyph`'s metrics in ems: its advance, whether it has ink,
    /// its box and a point in each hole.
    fn metrics_json(&mut self, s: &Spec, glyph: u16) -> Result<serde_json::Map<String, serde_json::Value>, String> {
        let g = self.glyph_outline(s, glyph)?;
        // Its loops' boxes as the polygons they are built as bound them
        // (rounded outward), so the field is outside past it.
        let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
        for l in &g.loops {
            let (a, b) = crate::profile::polygon_box(l)?;
            for k in 0..2 {
                lo[k] = lo[k].min(a[k]);
                hi[k] = hi[k].max(b[k]);
            }
        }
        let holes: Vec<serde_json::Value> = g.loops.iter().filter(|l| area(l) < 0.0).map(|l| serde_json::json!(interior(l))).collect();
        let j = serde_json::json!({
            "advance": g.advance,
            "ink": if g.loops.is_empty() { 0 } else { 1 },
            "box": if g.loops.is_empty() { serde_json::json!([0, 0, 0, 0]) } else { serde_json::json!([lo[0], lo[1], hi[0], hi[1]]) },
            "holes": holes,
        });
        Ok(j.as_object().expect("an object").clone())
    }

    fn glyph_metrics(&mut self, s: &Spec, glyph: u16, kern: f64) -> Result<Value, String> {
        let mut m = self.metrics_json(s, glyph)?;
        m.insert("kern".into(), serde_json::json!(kern));
        Ok(Value::Data(Data::new(serde_json::Value::Object(m))))
    }

    /// Text `t` shaped in font `s`: each line (split at "\n") `(vec
    /// glyphs advance clusters)`, each glyph `(vec outline data)`: its
    /// outline (glyph_2d's), and its pen "x" and offset "y" from the line's
    /// start (ems), the "cluster" it is in (counted from 0), its id
    /// ("glyph") and its metrics; the line's advance (ems) and its
    /// clusters counted.
    fn glyph_run(&mut self, s: &Spec, t: &str) -> Result<Value, String> {
        let mut lines = Vec::new();
        for line in t.split('\n') {
            let font = self.shaper(s)?;
            let (placed, adv, clusters) = shape::shape(&font, &s.key, line)?;
            let upem = font.units_per_em() as f64;
            let mut glyphs = Vec::new();
            for p in placed {
                let outline = self.glyph_2d(s, p.id)?;
                let mut m = self.metrics_json(s, p.id)?;
                m.insert("glyph".into(), serde_json::json!(p.id));
                m.insert("x".into(), serde_json::json!(p.x / upem));
                m.insert("y".into(), serde_json::json!(p.y / upem));
                m.insert("cluster".into(), serde_json::json!(p.cluster));
                glyphs.push(Value::Vec(Rc::new(vec![outline, Value::Data(Data::new(serde_json::Value::Object(m)))])));
            }
            let (a, c) = (self.g.constant(adv / upem), self.g.constant(clusters as f64));
            lines.push(Value::Vec(Rc::new(vec![Value::Vec(Rc::new(glyphs)), Value::Field(a, None), Value::Field(c, None)])));
        }
        Ok(Value::Vec(Rc::new(lines)))
    }

    /// Bounds on any line of characters from `chars` set in font `s`, in
    /// ems, over every glyph the font has (shape::bound): "glyphs-per-char",
    /// the most glyphs a character can become; "step", the least and most
    /// any glyph can advance the pen (from 0); "box", a box round any
    /// glyph's ink, where the font may move it, from its pen.
    fn glyph_reach(&mut self, s: &Spec, chars: &str) -> Result<Value, String> {
        // Every glyph of the font: its advance at the axes (exact, and the
        // whole units the shaper may round it to), its box.
        let (b, upem, (lo, hi), (b0, b1)) = self.with_face(s, |face| {
            let upem = face.units_per_em() as f64;
            let (mut lo, mut hi) = (f64::MAX, 0.0f64);
            let (mut b0, mut b1) = ([f64::MAX; 2], [f64::MIN; 2]);
            for g in 0..face.number_of_glyphs() {
                let [a_lo, a_hi] = shape::rounded(face.advance(g));
                (lo, hi) = (lo.min(a_lo / upem), hi.max(a_hi / upem));
                if let Some(x) = Self::drawn_box(face, g) {
                    b0 = [b0[0].min(x[0]), b0[1].min(x[1])];
                    b1 = [b1[0].max(x[2]), b1[1].max(x[3])];
                }
            }
            Ok((shape::bound(face, &s.key, chars)?, upem, (lo, hi), (b0, b1)))
        })?;
        let bx = if b0[0] == f64::MAX {
            serde_json::json!([0, 0, 0, 0])
        } else {
            serde_json::json!([b0[0] + b.dx[0] / upem, b0[1] + b.dy[0] / upem, b1[0] + b.dx[1] / upem, b1[1] + b.dy[1] / upem])
        };
        let step = [(lo + b.advance[0] / upem).min(0.0), (hi + b.advance[1] / upem).max(0.0)];
        Ok(Value::Data(Data::new(serde_json::json!({ "glyphs-per-char": b.per_char, "step": step, "box": bx }))))
    }
}

#[cfg(test)]
mod tests;
