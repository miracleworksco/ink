//! Support ratings, load frames, and placement checks.

use crate::model::Num;
use crate::{R, model, split};
use ink_lang::sexp::Sexp;
use ink_lang::view::ModelView;
use ink_lang::{InkError, Reader, Value, Xf};
use std::sync::Arc;

/// What support `s` carries, from its part's record, its `load-ratings`:
/// the `(attr spec RECORD)` of its body (a part written in the model), else
/// the `spec` the file defining it reads its record into (a catalog part's,
/// `b.spec` where it is imported as b).
pub(super) fn carries(r: &mut Reader, items: &[ink_lang::Item], s: &str, at: &Sexp) -> R<(split::Carries, bool, String)> {
    let below = format!("{s}/");
    let under: Vec<&ink_lang::Item> = items.iter().filter(|it| it.label == s || it.label.starts_with(&below)).collect();
    let hint = "a support is a placed part whose record says what it carries: its load-ratings' direction (radial, axial) or a guide's moment ratings \
                (MR, MP, MY), its body's (attr spec RECORD) or its file's spec (ink spec --section Dynamics/Supports)";
    let mut own: Vec<serde_json::Value> = Vec::new();
    for it in &under {
        let Some((_, Value::Data(d))) = it.attrs.iter().find(|(n, _)| n == "spec") else { continue };
        if d.value().get("load-ratings").is_some() && !own.contains(d.value()) {
            own.push(d.value().clone());
        }
    }
    let (record, spec) = match own.as_slice() {
        [one] => (one.clone(), "its (attr spec ...)".to_string()),
        [_, _, ..] => return Err(InkError { hint: Some(hint.into()), ..r.err(at, format!("the support {s} is not one part: its bodies give {} records", own.len())) }),
        [] => {
            let mut spaces: Vec<String> = Vec::new();
            for it in &under {
                let Some((name, _)) = &it.instance else { continue };
                let name = r.body(name).and_then(|b| b.of.clone()).unwrap_or_else(|| name.clone());
                let space = name.rsplit_once('.').map_or(String::new(), |x| format!("{}.", x.0));
                if !spaces.contains(&space) {
                    spaces.push(space);
                }
            }
            let [space] = spaces.as_slice() else {
                return Err(InkError { hint: Some(hint.into()), ..r.err(at, format!("the support {s} is not one part: it places parts of {} files", spaces.len())) });
            };
            let spec = format!("{space}spec");
            match r.value(&Sexp { k: ink_lang::sexp::SK::Sym(spec.clone()), line: at.line, col: at.col }) {
                Ok(Value::Data(d)) => (d.value().clone(), spec),
                _ => {
                    return Err(InkError {
                        hint: Some(hint.into()),
                        ..r.err(
                            at,
                            format!("the support {s}: its body gives no (attr spec RECORD) and its part's file defines no record, {spec} (its read-json of the part's record)"),
                        )
                    });
                }
            }
        }
    };
    let Some(ratings) = record.get("load-ratings") else {
        return Err(InkError { hint: Some(hint.into()), ..r.err(at, format!("the support {s}: its record, {spec}, has no load-ratings, so what it carries is not known")) });
    };
    let pulls = split::Carries::pulls(ratings);
    split::Carries::of(ratings).map(|c| (c, pulls, ratings.to_string())).ok_or_else(|| InkError {
        hint: Some(hint.into()),
        ..r.err(at, format!("the support {s}: its record's load-ratings ({spec}) give neither a direction (radial, axial) nor a guide's moment ratings (MR, MP, MY)"))
    })
}

/// Where support `s` carries its load: the frame its part declares,
/// `load` (the bearing's centre, z along its bore; a carriage's x along its
/// rail), as `(frame-of S load)` reads it; where S places an assembly that
/// does not declare one, the one part under it that does (a bushing given
/// a mass by an assembly around it). None, or several: an error naming
/// them.
pub(super) fn load_frame(r: &mut Reader, items: &[ink_lang::Item], s: &str, at: &Sexp) -> R<Xf> {
    let frame_of = |r: &mut Reader, p: &str| {
        let sym = |x: &str| Sexp { k: ink_lang::sexp::SK::Sym(x.into()), line: at.line, col: at.col };
        match r.value(&Sexp { k: ink_lang::sexp::SK::List(vec![sym("frame-of"), sym(p), sym("load")].into()), line: at.line, col: at.col }) {
            Ok(Value::Frame(x)) => Some(x),
            _ => None,
        }
    };
    if let Some(x) = frame_of(r, s) {
        return Ok(x);
    }
    let below = format!("{s}/");
    let mut inner: Vec<String> = Vec::new();
    for l in items.iter().flat_map(|it| &it.links).filter(|l| l.path.starts_with(&below)) {
        if !inner.contains(&l.path) && frame_of(r, &l.path).is_some() {
            inner.push(l.path.clone());
        }
    }
    if let [one] = inner.as_slice() {
        return Ok(frame_of(r, one).expect("its load frame"));
    }
    let hint = if inner.is_empty() {
        "a rated part's model declares (frame load ...) where its rating acts: a bearing's centre (an angular-contact or tapered roller \
         bearing's pressure centre), z along its bore; a carriage's on its rail, x along it (the catalog's do, as their generators write it). \
         An FF screw support declares none until its bearing's place is sourced"
            .to_string()
    } else {
        format!("name the part among those it places that do: {}", inner.join(", "))
    };
    Err(InkError { hint: Some(hint), ..r.err(at, format!("the support {s} declares no load frame, so where it carries is not known")) })
}

/// A part a reaction is split among: its path, what it carries, the frame
/// it carries at (its `load`), and the frame its share is given in (its
/// load frame where none is named), compiled with the params they read.
#[derive(Clone)]
pub(super) struct Support {
    pub(super) path: String,
    pub(super) carries: split::Carries,
    pub(super) pulls: bool,
    /// Its record's load-ratings as written (split::Part::rating).
    pub(super) rating: String,
    pub(super) load: (Arc<[Num; 12]>, Vec<usize>),
    pub(super) given: Option<(Arc<[Num; 12]>, Vec<usize>)>,
}

/// What the part at `place` most plainly holds, for an answer that says the line to write: the turning
/// joint under `target` whose axis it sits on (its centre off the axis by no more than it is across; the
/// nearest, of several), with the other parts on no joint that sit on that axis too and carry load ratings
/// (the other bearing of a pair). None where it sits on none. Nothing is assumed from this: it is said.
pub fn sits_on(v: &mut ModelView, target: &str, place: &str) -> Option<(String, Vec<String>)> {
    let items = v.doc.items(target).ok()?;
    let under = |label: &str, path: &str| label == path || label.strip_prefix(path).is_some_and(|rest| rest.starts_with('/'));
    type Box3 = ([f64; 3], [f64; 3]);
    let extent = |v: &mut ModelView, path: &str| -> Option<Box3> {
        let mut boxes = Vec::new();
        for it in items.iter().filter(|it| under(&it.label, path)) {
            boxes.extend(v.item_bounds(it));
        }
        (!boxes.is_empty()).then(|| {
            (
                std::array::from_fn(|k| boxes.iter().map(|b| b.0[k]).fold(f64::INFINITY, f64::min)),
                std::array::from_fn(|k| boxes.iter().map(|b| b.1[k]).fold(f64::NEG_INFINITY, f64::max)),
            )
        })
    };
    let off = |b: &Box3, o: [f64; 3], u: [f64; 3]| {
        let c: [f64; 3] = std::array::from_fn(|k| 0.5 * (b.0[k] + b.1[k]) - o[k]);
        let along = (0..3).map(|k| c[k] * u[k]).sum::<f64>();
        ((0..3).map(|k| (c[k] - along * u[k]).powi(2)).sum::<f64>().sqrt(), (0..3).map(|k| (b.1[k] - b.0[k]).powi(2)).sum::<f64>().sqrt())
    };
    let own = extent(v, place)?;
    let mut joints: Vec<String> = Vec::new();
    for l in items.iter().flat_map(|it| &it.links).filter(|l| l.joint.as_ref().is_some_and(|j| j.turn)) {
        if !joints.contains(&l.path) && !under(place, &l.path) {
            joints.push(l.path.clone());
        }
    }
    let mut best: Option<(f64, String, [f64; 3], [f64; 3])> = None;
    for path in joints {
        let Some((true, o, u)) = model::joint_axis(v, &items, target, &path) else { continue };
        let (d, across) = off(&own, o, u);
        if d <= across && best.as_ref().is_none_or(|b| d < b.0) {
            best = Some((d, path, o, u));
        }
    }
    let (_, joint, o, u) = best?;
    // The other rated parts on that axis that no joint moves: placed at the assembly's level, with load ratings.
    let jointed: Vec<&str> = items.iter().flat_map(|it| &it.links).filter(|l| l.joint.is_some()).map(|l| l.path.as_str()).collect();
    let mut others: Vec<String> = Vec::new();
    for it in &items {
        let rated = it.attrs.iter().any(|(n, value)| match (n.as_str(), value) {
            ("spec", Value::Data(d)) => d.value().get("load-ratings").is_some(),
            ("ratings", Value::Data(d)) => d.value().as_object().is_some_and(|o| !o.is_empty()),
            _ => false,
        });
        let top = it.label.as_str();
        if !rated || under(top, place) || under(place, top) || jointed.iter().any(|j| under(top, j)) || others.iter().any(|o| o == top) {
            continue;
        }
        if let Some(b) = extent(v, top)
            && off(&b, o, u).0 <= off(&b, o, u).1
        {
            others.push(top.to_string());
        }
    }
    Some((joint, others))
}

/// The supports that touch nothing but what they hold, by their boxes: each of `held`'s supports (a joint's
/// path with the parts named as holding it) whose box is more than a millimetre clear of every body under
/// `target` that is neither its own nor carried by that joint. A box clear of all the rest cannot touch
/// any of it, so the part is in mid-air: nothing passes its load on. (Boxes that meet prove no contact:
/// `ink interfere` and a render show the rest.)
pub fn supports_floating(v: &mut ModelView, target: &str, held: &[(String, Vec<String>)]) -> Vec<(String, String)> {
    let Ok(items) = v.doc.items(target) else { return vec![] };
    let paths: Vec<String> = held.iter().map(|(path, _)| path.clone()).collect();
    let (carried, _) = model::carried(&items, &paths);
    let under = |label: &str, path: &str| label == path || label.strip_prefix(path).is_some_and(|rest| rest.starts_with('/'));
    let mut out: Vec<(String, String)> = vec![];
    for ((path, supports), carried) in held.iter().zip(&carried) {
        for s in supports {
            if out.iter().any(|(o, _)| o == s) {
                continue;
            }
            let mut own = Vec::new();
            for it in items.iter().filter(|it| under(&it.label, s)) {
                own.extend(v.item_bounds(it));
            }
            if own.is_empty() {
                continue;
            }
            // Every other body: not its own, not what the joint carries, not another support of the same joint.
            let mut meets = false;
            for it in items.iter().filter(|it| it.instance.is_some() && !under(&it.label, s) && !carried.contains(&it.label) && !supports.iter().any(|o| under(&it.label, o))) {
                let Some(b) = v.item_bounds(it) else {
                    meets = true; // not proven bounded: it may reach anywhere
                    break;
                };
                if own.iter().any(|a| (0..3).all(|k| a.0[k] - b.1[k] <= 1.0 + a.2 + b.2 && b.0[k] - a.1[k] <= 1.0 + a.2 + b.2)) {
                    meets = true;
                    break;
                }
            }
            if !meets {
                out.push((s.clone(), path.clone()));
            }
        }
    }
    out
}

/// A support holds what its joint carries, so it is at it: Err where one of `held` (each a joint's path and
/// the supports the design names for it) is clear of every body the joint carries by more than they and it
/// are across. It holds nothing there, and its share would be of a load it is not under (bearings placed
/// 300 mm from their wheel were given the wheel's). Bodies not proven bounded say nothing.
pub fn supports_near(v: &mut ModelView, target: &str, held: &[(String, Vec<String>)]) -> R<()> {
    if held.is_empty() {
        return Ok(());
    }
    let items = v.doc.items(target).map_err(InkError::msg)?;
    let paths: Vec<String> = held.iter().map(|(path, _)| path.clone()).collect();
    let (carried, _) = model::carried(&items, &paths);
    let under = |label: &str, path: &str| label == path || label.strip_prefix(path).is_some_and(|rest| rest.starts_with('/'));
    for ((path, supports), carried) in held.iter().zip(&carried) {
        let mut boxes = Vec::new();
        for it in items.iter().filter(|it| carried.contains(&it.label) && !supports.iter().any(|s| under(&it.label, s))) {
            boxes.extend(v.item_bounds(it));
        }
        for s in supports.iter().filter(|s| !under(path, s) && !under(s, path)) {
            let mut own = Vec::new();
            for it in items.iter().filter(|it| under(&it.label, s)) {
                own.extend(v.item_bounds(it));
            }
            type Box3 = ([f64; 3], [f64; 3], f64);
            // Nowhere near: clear of every body held by more than the support and all of them together are
            // across (a sketch places a bearing beside its wheel; one a wheel's width and more away is elsewhere).
            let across = |boxes: &[Box3]| {
                let (lo, hi): ([f64; 3], [f64; 3]) = (
                    std::array::from_fn(|k| boxes.iter().map(|b| b.0[k]).fold(f64::INFINITY, f64::min)),
                    std::array::from_fn(|k| boxes.iter().map(|b| b.1[k]).fold(f64::NEG_INFINITY, f64::max)),
                );
                (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt()
            };
            if own.is_empty() || boxes.is_empty() {
                continue;
            }
            // What turns is held along its axis, wherever a sketch leaves the shaft undrawn: a support of
            // it is off the axis by no more than it is across. What slides is held where it is.
            if let Some((true, o, u)) = model::joint_axis(v, &items, target, path) {
                let (lo, hi): ([f64; 3], [f64; 3]) = (
                    std::array::from_fn(|k| own.iter().map(|b| b.0[k]).fold(f64::INFINITY, f64::min)),
                    std::array::from_fn(|k| own.iter().map(|b| b.1[k]).fold(f64::NEG_INFINITY, f64::max)),
                );
                let c: [f64; 3] = std::array::from_fn(|k| 0.5 * (lo[k] + hi[k]) - o[k]);
                let along = (0..3).map(|k| c[k] * u[k]).sum::<f64>();
                let off = (0..3).map(|k| (c[k] - along * u[k]).powi(2)).sum::<f64>().sqrt();
                if off > across(&own) {
                    return Err(InkError::msg(format!(
                        "{s} is named as holding what turns at {path}, and is {off:.0} mm off that joint's axis, more than it is across: it holds nothing there. A support \
                         of what turns is on its axis (a bearing on the shaft): place it there, or name the parts that do hold it"
                    )));
                }
                continue;
            }
            let reach = across(&own) + across(&boxes);
            let near = |a: &Box3, b: &Box3| (0..3).all(|k| a.0[k] - b.1[k] <= reach && b.0[k] - a.1[k] <= reach);
            if !own.iter().any(|a| boxes.iter().any(|b| near(a, b))) {
                return Err(InkError::msg(format!(
                    "{s} is named as holding what {path} carries, and is nowhere near it: it is clear of every body that joint carries by more than the {:.0} mm they and \
                     it are across together. A support is at what it holds (a carriage under its table): place it there, or name the parts that do hold it",
                    reach
                )));
            }
        }
    }
    Ok(())
}
