//! Text shaped into glyphs: a line set as the font sets it (HarfBuzz's
//! algorithm, harfrust), and a bound on every line a set of characters
//! can make that holds whatever the font substitutes.
//!
//! The shaper reads the font's GSUB and GPOS: its ligatures, contextual
//! forms, composed accents and marks placed on their bases, and its
//! kerning, at the face's variation axes. Lines run left to right.
//!
//! The bound does not follow which glyphs a text reaches: any glyph of
//! the font may appear. It reads the font as the shaper does (face.rs: the
//! same reader, the same instance, the same variation arithmetic), walking
//! GSUB's and GPOS's lookup lists in their bytes, and counts how often each
//! lookup can act on one glyph (from the stages the shaper runs and the
//! calls contextual substitutions make) and the most each lookup can move
//! or widen a glyph or multiply the glyphs. Any table, lookup type or
//! format it does not interpret it refuses, as it does scripts with their
//! own shapers, marks, cursive attachment, contextual positioning and
//! contextual substitutions calling contextual ones.

use super::face::Face;
use std::collections::HashMap;
use unicode_properties::{GeneralCategoryGroup, UnicodeGeneralCategory};

/// A glyph of a shaped line, in font units: its id, where it is drawn
/// (its pen, and with its offset x and y), and the cluster it is in.
pub(super) struct Placed {
    pub id: u16,
    pub x: f64,
    pub y: f64,
    pub pen: f64,
    pub cluster: usize,
}

/// Line `text` shaped in `face` (its axes set): its glyphs, its advance,
/// and the characters each cluster stands for counted (one cluster each).
/// A character the font lacks is an error naming it.
pub(super) fn shape(font: &harfrust::Font, key: &str, text: &str) -> Result<(Vec<Placed>, f64, usize), String> {
    let sf = harfrust::ShaperFont::new(font);
    let mut buf = harfrust::Buffer::new();
    buf.push_str(text);
    buf.guess_segment_properties();
    buf.set_direction(harfrust::Direction::LeftToRight);
    harfrust::shape(&sf, &mut buf, harfrust::ShapeOptions::new()).map_err(|e| format!("font {key} cannot shape {text:?}: {e:?}"))?;
    let (mut pen, mut placed, mut clusters) = (0.0, Vec::new(), Vec::<u32>::new());
    for (i, p) in buf.glyph_infos().iter().zip(buf.glyph_positions()) {
        if i.glyph_id == 0 {
            let c = text[i.cluster as usize..].chars().next().unwrap_or('?');
            return Err(format!("font {key} has no glyph for {}", super::char_name(c as u32)));
        }
        if clusters.last() != Some(&i.cluster) {
            clusters.push(i.cluster);
        }
        placed.push(Placed { id: i.glyph_id as u16, x: pen + p.x_offset as f64, y: p.y_offset as f64, pen, cluster: clusters.len() - 1 });
        pen += p.x_advance as f64;
    }
    Ok((placed, pen, clusters.len()))
}

type R<T> = Result<T, String>;

fn u16_at(d: &[u8], at: usize) -> Option<u16> {
    d.get(at..at + 2).map(|b| u16::from_be_bytes([b[0], b[1]]))
}

fn u32_at(d: &[u8], at: usize) -> Option<u32> {
    d.get(at..at + 4).map(|b| u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
}

/// A table's u16 at `at`, or the error a read past its end is.
fn rd(t: &[u8], at: usize) -> R<usize> {
    u16_at(t, at).map(usize::from).ok_or_else(|| "has a layout table ink cannot read (it ends early)".into())
}

fn rd32(t: &[u8], at: usize) -> R<usize> {
    u32_at(t, at).map(|v| v as usize).ok_or_else(|| "has a layout table ink cannot read (it ends early)".into())
}

/// The coverage table at `at`: each glyph it covers, with its index.
fn coverage(t: &[u8], at: usize) -> R<Vec<(u16, usize)>> {
    let mut out = Vec::new();
    match rd(t, at)? {
        1 => {
            for i in 0..rd(t, at + 2)? {
                out.push((rd(t, at + 4 + 2 * i)? as u16, i));
            }
        }
        2 => {
            for r in 0..rd(t, at + 2)? {
                let (start, end, first) = (rd(t, at + 4 + 6 * r)?, rd(t, at + 6 + 6 * r)?, rd(t, at + 8 + 6 * r)?);
                for g in start..=end.max(start) {
                    if g <= end {
                        out.push((g as u16, first + g - start));
                    }
                }
            }
        }
        f => return Err(format!("has a coverage table of format {f}")),
    }
    Ok(out)
}

/// The class definition at `at`: each glyph's class (those not given, 0).
fn classes(t: &[u8], at: usize) -> R<HashMap<u16, usize>> {
    let mut out = HashMap::new();
    match rd(t, at)? {
        1 => {
            let start = rd(t, at + 2)?;
            for i in 0..rd(t, at + 4)? {
                out.insert((start + i) as u16, rd(t, at + 6 + 2 * i)?);
            }
        }
        2 => {
            for r in 0..rd(t, at + 2)? {
                let (start, end, class) = (rd(t, at + 4 + 6 * r)?, rd(t, at + 6 + 6 * r)?, rd(t, at + 8 + 6 * r)?);
                for g in start..=end.max(start) {
                    if g <= end {
                        out.insert(g as u16, class);
                    }
                }
            }
        }
        f => return Err(format!("has a class definition of format {f}")),
    }
    Ok(out)
}

/// A GSUB or GPOS table's bytes, walked as they are.
struct Layout<'a> {
    t: &'a [u8],
    /// The extension lookup type: 7 in GSUB, 9 in GPOS.
    ext: usize,
    name: &'static str,
}

impl<'a> Layout<'a> {
    fn new(t: &'a [u8], gsub: bool) -> R<Layout<'a>> {
        let name = if gsub { "GSUB" } else { "GPOS" };
        let (major, minor) = (rd(t, 0)?, rd(t, 2)?);
        if major != 1 || minor > 1 {
            return Err(format!("has a {name} table of version {major}.{minor}"));
        }
        Ok(Layout { t, ext: if gsub { 7 } else { 9 }, name })
    }

    fn lookups(&self) -> R<usize> {
        rd(self.t, rd(self.t, 8)?)
    }

    /// Lookup `li`: its type and its subtables' offsets, an extension's
    /// resolved to what it wraps (all of one type, format 1).
    fn lookup(&self, li: usize) -> R<(usize, Vec<usize>)> {
        let t = self.t;
        let list = rd(t, 8)?;
        let lookup = list + rd(t, list + 2 + 2 * li)?;
        let mut kind = rd(t, lookup)?;
        let wrapped = kind == self.ext;
        let mut inner = None;
        let mut out = Vec::new();
        for si in 0..rd(t, lookup + 4)? {
            let mut sub = lookup + rd(t, lookup + 6 + 2 * si)?;
            if wrapped {
                let (f, k) = (rd(t, sub)?, rd(t, sub + 2)?);
                if f != 1 || inner.is_some_and(|i| i != k) || k == self.ext {
                    return Err(format!("has an extension in {} lookup {li} ink cannot read", self.name));
                }
                inner = Some(k);
                sub += rd32(t, sub + 4)?;
            }
            out.push(sub);
        }
        if wrapped {
            kind = inner.unwrap_or(0);
        }
        Ok((kind, out))
    }

    /// The lookups every version of each feature lists: its own, and those
    /// a FeatureVariations record puts in its place at some axes (all
    /// counted, whatever the axes). Each with the feature's index and tag.
    fn feature_versions(&self) -> R<Vec<(usize, [u8; 4], Vec<usize>)>> {
        let t = self.t;
        let list = rd(t, 6)?;
        let mut tags = Vec::new();
        let mut out = Vec::new();
        for i in 0..rd(t, list)? {
            let rec = list + 2 + 6 * i;
            let tag: [u8; 4] = t.get(rec..rec + 4).and_then(|b| b.try_into().ok()).ok_or("has a layout table ink cannot read (it ends early)")?;
            let f = list + rd(t, rec + 4)?;
            let lookups = (0..rd(t, f + 2)?).map(|j| rd(t, f + 4 + 2 * j)).collect::<R<Vec<usize>>>()?;
            tags.push(tag);
            out.push((i, tag, lookups));
        }
        if rd(t, 2)? == 1 && rd32(t, 10)? != 0 {
            let fv = rd32(t, 10)?;
            for r in 0..rd32(t, fv + 4)? {
                let subst = fv + rd32(t, fv + 8 + 8 * r + 4)?;
                for k in 0..rd(t, subst + 4)? {
                    let at = subst + 6 + 6 * k;
                    let fi = rd(t, at)?;
                    let alt = subst + rd32(t, at + 2)?;
                    let lookups = (0..rd(t, alt + 2)?).map(|j| rd(t, alt + 4 + 2 * j)).collect::<R<Vec<usize>>>()?;
                    let tag = *tags.get(fi).ok_or("has a feature variation of no feature")?;
                    out.push((fi, tag, lookups));
                }
            }
        }
        Ok(out)
    }

    /// The features some language system requires.
    fn required(&self) -> R<Vec<usize>> {
        let t = self.t;
        let list = rd(t, 4)?;
        let mut out = Vec::new();
        let mut langsys = |ls: usize| -> R<()> {
            let req = rd(t, ls + 2)?;
            if req != 0xffff {
                out.push(req);
            }
            Ok(())
        };
        for i in 0..rd(t, list)? {
            let script = list + rd(t, list + 2 + 6 * i + 4)?;
            let def = rd(t, script)?;
            if def != 0 {
                langsys(script + def)?;
            }
            for j in 0..rd(t, script + 2)? {
                langsys(script + rd(t, script + 4 + 6 * j + 4)?)?;
            }
        }
        Ok(out)
    }

    /// The lookups a contextual subtable (format `f`, chained or not, at
    /// `sub`) calls, every call of every rule.
    fn calls(&self, sub: usize, f: usize, chained: bool) -> R<Vec<usize>> {
        let t = self.t;
        let mut out = Vec::new();
        let records = |at: usize, n: usize, out: &mut Vec<usize>| -> R<()> {
            for k in 0..n {
                out.push(rd(t, at + 4 * k + 2)?);
            }
            Ok(())
        };
        // A rule: its lookup records, after its glyph or class sequences.
        let rule = |r: usize, out: &mut Vec<usize>| -> R<()> {
            if chained {
                let back = rd(t, r)?;
                let input = r + 2 + 2 * back;
                let ni = rd(t, input)?;
                let ahead = input + 2 + 2 * ni.saturating_sub(1);
                let na = rd(t, ahead)?;
                let recs = ahead + 2 + 2 * na;
                records(recs + 2, rd(t, recs)?, out)
            } else {
                let n = rd(t, r)?;
                records(r + 4 + 2 * n.saturating_sub(1), rd(t, r + 2)?, out)
            }
        };
        match (chained, f) {
            (_, 1) | (_, 2) => {
                // Rule sets, by index (a set may be NULL).
                let (count_at, first) = match (chained, f) {
                    (false, 1) => (sub + 4, sub + 6),
                    (false, _) => (sub + 6, sub + 8),
                    (true, 1) => (sub + 4, sub + 6),
                    (true, _) => (sub + 10, sub + 12),
                };
                for i in 0..rd(t, count_at)? {
                    let off = rd(t, first + 2 * i)?;
                    if off == 0 {
                        continue;
                    }
                    let set = sub + off;
                    for j in 0..rd(t, set)? {
                        let r = rd(t, set + 2 + 2 * j)?;
                        if r != 0 {
                            rule(set + r, &mut out)?;
                        }
                    }
                }
            }
            (false, 3) => {
                let n = rd(t, sub + 2)?;
                records(sub + 6 + 2 * n, rd(t, sub + 4)?, &mut out)?;
            }
            (true, 3) => {
                let back = rd(t, sub + 2)?;
                let input = sub + 4 + 2 * back;
                let ahead = input + 2 + 2 * rd(t, input)?;
                let recs = ahead + 2 + 2 * rd(t, ahead)?;
                records(recs + 2, rd(t, recs)?, &mut out)?;
            }
            _ => return Err(format!("has a contextual {} subtable of format {f}", self.name)),
        }
        Ok(out)
    }
}

/// How many times each lookup can act on one glyph of a line the default
/// shaper sets (Latin, Greek, Cyrillic). It runs each lookup once a stage:
/// GSUB has two (rvrn's features and the required one; then every other
/// feature and the required one), GPOS one. A lookup any feature lists
/// (any tag, any version) counts once for each stage it can be in, and
/// once more for each call a contextual rule makes, each time its caller
/// acts: a contextual lookup's matches do not overlap, so each call acts
/// on a glyph once a time its caller acts. `calls` is each lookup's calls;
/// a contextual lookup calling another is refused (its matches may start
/// at each glyph before the one it acts on).
fn applications(layout: &Layout, gsub: bool, calls: &[Vec<usize>]) -> R<Vec<u64>> {
    let n = calls.len();
    let required = layout.required()?;
    let (mut first, mut rest) = (vec![false; n], vec![false; n]);
    for (fi, tag, ls) in layout.feature_versions()? {
        let req = required.contains(&fi);
        for l in ls {
            if l < n {
                if !gsub || &tag != b"rvrn" || req {
                    rest[l] = true;
                }
                if gsub && (&tag == b"rvrn" || req) {
                    first[l] = true;
                }
            }
        }
    }
    let stages: Vec<u64> = (0..n).map(|l| first[l] as u64 + rest[l] as u64).collect();
    let mut m = stages.clone();
    for (caller, cs) in calls.iter().enumerate() {
        for &c in cs {
            if c < n && !calls[c].is_empty() {
                return Err(format!("has a contextual lookup ({caller}) calling another ({c})"));
            }
            if c < n {
                m[c] += stages[caller];
            }
        }
    }
    Ok(m)
}

/// The value record of format `f` at `at`: its x advance, x placement and
/// y placement, each with its variation device's delta (`dev`), a device's
/// offset counted from `base`.
fn value_at(t: &[u8], at: usize, f: usize, base: usize, dev: &dyn Fn(u16, u16) -> f64) -> R<[f64; 3]> {
    let mut v = [0.0f64; 8];
    let mut p = at;
    for (k, x) in v.iter_mut().enumerate() {
        if f & (1 << k) != 0 {
            let w = rd(t, p)?;
            *x = if k < 4 {
                w as u16 as i16 as f64
            } else if w == 0 {
                0.0
            } else {
                let d = base + w;
                if rd(t, d + 4)? == 0x8000 { dev(rd(t, d)? as u16, rd(t, d + 2)? as u16) } else { 0.0 }
            };
            p += 2;
        }
    }
    Ok([v[2] + v[6], v[0] + v[4], v[1] + v[5]])
}

/// The whole units the shaper may round `v` to. It turns an exact value
/// into f32 before it rounds it (a variation delta; an advance with a
/// fraction), and past 2^24 that f32 may already lie a unit or more from
/// it: the unit below each of the two, and the unit above each. (Past
/// i32 it saturates or wraps, toward 0, which every range here holds.)
pub(super) fn rounded(v: f64) -> [f64; 2] {
    let c = v as f32 as f64;
    [v.floor().min(c.floor()), v.ceil().max(c.ceil())]
}

/// The length of a value record of format `f`.
fn value_len(f: usize) -> usize {
    2 * (f & 0xff).count_ones() as usize
}

/// What any line of characters from a set can be, whatever glyphs the font
/// puts in it: how many glyphs it can hold per character, and what the
/// font's positioning can add to any glyph's advance and offsets (font
/// units, each a range from 0).
pub(super) struct Bound {
    pub per_char: u64,
    pub advance: [f64; 2],
    pub dx: [f64; 2],
    pub dy: [f64; 2],
}

/// The scripts the default shaper sets, the only ones counted here.
const SCRIPTS: &[harfrust::Script] = &[harfrust::Script::LATIN, harfrust::Script::GREEK, harfrust::Script::CYRILLIC, harfrust::Script::COMMON, harfrust::Script::INHERITED];

/// The bound on any line of the characters `chars` in `face`. Refused,
/// saying why, where it is not counted simply: a character the font lacks,
/// a combining mark or one set as several glyphs; a script with its own
/// shaper; a font shaped by its Apple tables or kerned by a `kern` table;
/// a table, lookup type or format not read here; a substitution that makes
/// a mark of another glyph, a mark positioned on a glyph that is not one,
/// cursive attachment, contextual positioning; a contextual substitution
/// calling another.
pub(super) fn bound(face: &Face, key: &str, chars: &str) -> R<Bound> {
    fn refuse<T>(why: String) -> R<T> {
        Err(format!("{why}, which lettering-max does not bound: name the texts to check as cases"))
    }
    let font_says = |w: String| format!("font {key} {w}");
    for t in [b"morx", b"mort", b"kerx", b"trak", b"kern"] {
        if face.table(t).is_some() {
            return refuse(format!("font {key} has a {} table", String::from_utf8_lossy(t)));
        }
    }
    // Marks: GDEF's glyph classes (versions 1.0, 1.2 and 1.3 read here).
    let mut class = HashMap::new();
    if let Some(g) = face.table(b"GDEF") {
        let v = (rd(g, 0).map_err(font_says)?, rd(g, 2).map_err(font_says)?);
        if !matches!(v, (1, 0) | (1, 2) | (1, 3)) {
            return refuse(format!("font {key} has a GDEF table of version {}.{}", v.0, v.1));
        }
        let at = rd(g, 4).map_err(font_says)?;
        if at != 0 {
            class = classes(g, at).or_else(|w| refuse(font_says(w)))?;
        }
    }
    let is_mark = |g: u16| class.get(&g) == Some(&3);
    let font = face.font;
    for c in chars.chars().filter(|&c| c != '\n') {
        let name = super::char_name(c as u32);
        let Some(g) = face.glyph_index(c) else { return Err(format!("font {key} has no glyph for {name}")) };
        if c.general_category_group() == GeneralCategoryGroup::Mark || is_mark(g) {
            return refuse(format!("{name} is a combining mark, placed on the letter before it"));
        }
        let sf = harfrust::ShaperFont::new(font);
        let mut buf = harfrust::Buffer::new();
        buf.push_str(&c.to_string());
        buf.guess_segment_properties();
        if let Some(sc) = buf.script()
            && !SCRIPTS.contains(&sc)
        {
            return refuse(format!("{name} is in a script with its own shaping rules"));
        }
        buf.set_direction(harfrust::Direction::LeftToRight);
        harfrust::shape(&sf, &mut buf, harfrust::ShapeOptions::new()).map_err(|e| format!("font {key} cannot shape {name}: {e:?}"))?;
        if buf.glyph_infos().len() > 1 {
            return refuse(format!("{name} sets as {} glyphs in font {key}", buf.glyph_infos().len()));
        }
    }
    let mut per_char = 1u64;
    if let Some(t) = face.table(b"GSUB") {
        let layout = Layout::new(t, true).or_else(|w| refuse(font_says(w)))?;
        let n = layout.lookups().or_else(|w| refuse(font_says(w)))?;
        let lookups = (0..n).map(|li| layout.lookup(li)).collect::<R<Vec<_>>>().or_else(|w| refuse(font_says(w)))?;
        let mut calls = vec![Vec::new(); n];
        for (li, (kind, subs)) in lookups.iter().enumerate() {
            for &sub in subs {
                if *kind == 5 || *kind == 6 {
                    let f = rd(t, sub).or_else(|w| refuse(font_says(w)))?;
                    calls[li].extend(layout.calls(sub, f, *kind == 6).or_else(|w| refuse(font_says(w)))?);
                }
            }
        }
        let m = applications(&layout, true, &calls).or_else(|w| refuse(font_says(w)))?;
        for (li, (kind, subs)) in lookups.iter().enumerate() {
            let mut longest = 1u64;
            for &sub in subs {
                let f = rd(t, sub).or_else(|w| refuse(font_says(w)))?;
                // Each covered glyph and what it may become.
                let outs = |sub: usize| -> R<Vec<(u16, Vec<u16>)>> {
                    let cov = coverage(t, sub + rd(t, sub + 2)?)?;
                    let list = |at: usize| -> R<Vec<u16>> { (0..rd(t, at)?).map(|k| rd(t, at + 2 + 2 * k).map(|g| g as u16)).collect() };
                    let mut out = Vec::new();
                    for (g, i) in cov {
                        let o = match (kind, f) {
                            (1, 1) => vec![(g as i32 + rd(t, sub + 4)? as u16 as i16 as i32) as u16],
                            (1, 2) => vec![rd(t, sub + 6 + 2 * i)? as u16],
                            (2, 1) | (3, 1) => list(sub + rd(t, sub + 6 + 2 * i)?)?,
                            (4, 1) => {
                                let set = sub + rd(t, sub + 6 + 2 * i)?;
                                (0..rd(t, set)?).map(|k| Ok(rd(t, set + rd(t, set + 2 + 2 * k)?)? as u16)).collect::<R<Vec<u16>>>()?
                            }
                            _ => unreachable!(),
                        };
                        out.push((g, o));
                    }
                    Ok(out)
                };
                let pairs = match (kind, f) {
                    (1, 1) | (1, 2) | (2, 1) | (3, 1) | (4, 1) => outs(sub).or_else(|w| refuse(font_says(w)))?,
                    (5, 1..=3) | (6, 1..=3) => Vec::new(),
                    (8, 1) => {
                        let rev = || -> R<Vec<(u16, Vec<u16>)>> {
                            let back = rd(t, sub + 4)?;
                            let ahead = sub + 6 + 2 * back;
                            let subs = ahead + 2 + 2 * rd(t, ahead)?;
                            coverage(t, sub + rd(t, sub + 2)?)?.into_iter().map(|(g, i)| Ok((g, vec![rd(t, subs + 2 + 2 * i)? as u16]))).collect()
                        };
                        rev().or_else(|w| refuse(font_says(w)))?
                    }
                    _ => return refuse(format!("font {key} has a GSUB lookup ({li}) of type {kind}, format {f}, not read here")),
                };
                for (g, o) in pairs {
                    if *kind == 2 {
                        longest = longest.max(o.len() as u64);
                    }
                    if m[li] > 0
                        && !is_mark(g)
                        && let Some(x) = o.into_iter().find(|&x| is_mark(x))
                    {
                        return refuse(format!("font {key} substitutes glyph {x}, a mark, for glyph {g}"));
                    }
                }
            }
            // Each time the lookup acts, a glyph may become `longest`.
            for _ in 0..m[li].min(64) {
                per_char = per_char.saturating_mul(longest);
            }
        }
    }
    let (mut advance, mut dx, mut dy) = ([0.0f64; 2], [0.0f64; 2], [0.0f64; 2]);
    if let Some(t) = face.table(b"GPOS") {
        let layout = Layout::new(t, false).or_else(|w| refuse(font_says(w)))?;
        let n = layout.lookups().or_else(|w| refuse(font_says(w)))?;
        let lookups = (0..n).map(|li| layout.lookup(li)).collect::<R<Vec<_>>>().or_else(|w| refuse(font_says(w)))?;
        let m = applications(&layout, false, &vec![Vec::new(); n]).or_else(|w| refuse(font_says(w)))?;
        let store = face.var_store().or_else(|w| refuse(font_says(w)))?;
        let coords = face.coords();
        // A variation delta as the shaper computes it, and the whole units
        // it may round it to.
        let delta = |outer: u16, inner: u16| store.as_ref().map_or(0.0, |s| Face::delta(s, coords, outer, inner));
        let below = |outer: u16, inner: u16| rounded(delta(outer, inner))[0];
        let above = |outer: u16, inner: u16| rounded(delta(outer, inner))[1];
        for (li, (kind, subs)) in lookups.iter().enumerate() {
            // Each value record, as a pair's first glyph and its second.
            let mut sides: [Vec<[f64; 3]>; 2] = [Vec::new(), Vec::new()];
            for &sub in subs {
                let f = rd(t, sub).or_else(|w| refuse(font_says(w)))?;
                let read = |dev: &dyn Fn(u16, u16) -> f64, sides: &mut [Vec<[f64; 3]>; 2]| -> R<()> {
                    match (kind, f) {
                        (1, 1) => sides[0].push(value_at(t, sub + 6, rd(t, sub + 4)?, sub, dev)?),
                        (1, 2) => {
                            let vf = rd(t, sub + 4)?;
                            for i in 0..rd(t, sub + 6)? {
                                sides[0].push(value_at(t, sub + 8 + i * value_len(vf), vf, sub, dev)?);
                            }
                        }
                        (2, 1) => {
                            let (f1, f2) = (rd(t, sub + 4)?, rd(t, sub + 6)?);
                            let rec = 2 + value_len(f1) + value_len(f2);
                            for k in 0..rd(t, sub + 8)? {
                                let set = sub + rd(t, sub + 10 + 2 * k)?;
                                for i in 0..rd(t, set)? {
                                    let at = set + 2 + i * rec + 2;
                                    sides[0].push(value_at(t, at, f1, set, dev)?);
                                    sides[1].push(value_at(t, at + value_len(f1), f2, set, dev)?);
                                }
                            }
                        }
                        (2, 2) => {
                            let (f1, f2) = (rd(t, sub + 4)?, rd(t, sub + 6)?);
                            let rec = value_len(f1) + value_len(f2);
                            for i in 0..rd(t, sub + 12)? * rd(t, sub + 14)? {
                                let at = sub + 16 + i * rec;
                                sides[0].push(value_at(t, at, f1, sub, dev)?);
                                sides[1].push(value_at(t, at + value_len(f1), f2, sub, dev)?);
                            }
                        }
                        _ => {}
                    }
                    Ok(())
                };
                match (kind, f) {
                    (1, 1) | (1, 2) | (2, 1) | (2, 2) => {
                        read(&below, &mut sides).or_else(|w| refuse(font_says(w)))?;
                        read(&above, &mut sides).or_else(|w| refuse(font_says(w)))?;
                    }
                    (3, _) => return refuse(format!("font {key} joins glyphs by cursive attachment")),
                    // Marks only: every glyph they place must be one.
                    (4, 1) | (5, 1) | (6, 1) => {
                        let cov = coverage(t, sub + rd(t, sub + 2).or_else(|w| refuse(font_says(w)))?).or_else(|w| refuse(font_says(w)))?;
                        if let Some((g, _)) = cov.into_iter().find(|&(g, _)| !is_mark(g)) {
                            return refuse(format!("font {key} positions glyph {g} as a mark, though it is not one"));
                        }
                    }
                    (7, _) | (8, _) => return refuse(format!("font {key} positions glyphs by context (GPOS lookup {li})")),
                    _ => return refuse(format!("font {key} has a GPOS lookup ({li}) of type {kind}, format {f}, not read here")),
                }
            }
            let times = m[li] as f64;
            for side in sides {
                for (k, r) in [&mut advance, &mut dx, &mut dy].into_iter().enumerate() {
                    let lo = side.iter().map(|x| x[k]).fold(0.0f64, f64::min);
                    let hi = side.iter().map(|x| x[k]).fold(0.0f64, f64::max);
                    r[0] += times * lo;
                    r[1] += times * hi;
                }
            }
        }
    }
    Ok(Bound { per_char, advance, dx, dy })
}
