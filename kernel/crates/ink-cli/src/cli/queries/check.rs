//! The model overview and its diagnostics.

use super::target::EMPTY_CAUSE;
use crate::{Arc, Args, CliView, InkError, J, Session, Tape, Value, Xf, bounds, extend, inexact_fields, json, measured, n, param_names, sources_json, unchecked_fields};
use ink::graph::Var;
use ink::lang::Callable;
use rayon::prelude::*;

/// Named constructions (defs, tags) in a body's field not seen forming any
/// of its surface (feature::surface_names): proven to form none (hidden
/// inside another shape, or cut away entirely), or never seen deciding a
/// patch of it where it was sought (only an edge, a blend's covered side,
/// or a face the search missed). Usually an authoring slip (a mirrored
/// half that vanished, a hole outside the part). Then fields interval
/// arithmetic cannot bound tightly, and undefined cases (feature::survey).
/// `named`: the body's named parts.
fn survey_warnings(s: &Session, body: &str, mut named: Vec<u32>, survey: &ink::feature::Survey, moments: usize) -> Vec<String> {
    named.retain(|t| !survey.formed.contains(t));
    let when = if moments > 1 { format!(" at the {moments} moments of its motion it was sought at") } else { String::new() };
    let mut out: Vec<String> = named
        .into_iter()
        .map(|t| {
            let name = &s.doc.g.tag_names[t as usize];
            if let Some(&(_, n, finer)) = survey.unseen.iter().find(|u| u.0 == t) {
                let cells = match n {
                    0 => " (on the field's path in no cell at those, but not proven off it at every moment)".to_string(),
                    _ => format!(": in the {n} cells it is on the field's path in, finer in {finer} of them"),
                };
                format!(
                    "body '{body}': '{name}' was not seen forming any of its surface{} where it was sought{cells}. It may form only an edge \
                     (where it ties with what covers it) or a blend's covered side, or a face the search missed. Undecided, not proven",
                    when
                )
            } else {
                format!(
                    "body '{body}': '{name}' is in its field but forms none of its surface. It is hidden inside another shape or cut away \
                     (fine for an envelope intersected on purpose)"
                )
            }
        })
        .collect();
    if let Some(l) = &survey.loose {
        out.push(format!(
            "body '{body}': interval bounds of its field are {:.0}x wider than its variation on {:.0}% of its surface (up to {:.0}x, near {:?}), \
             so proofs and meshing refine deeper than the geometry needs. Usual causes: smooth blends and other fields that aren't distances \
             (try (redistance f)), atan2 or mod wraps, and a subexpression used twice",
            l.p90,
            100.0 * l.share,
            l.max,
            l.at.map(|v| (v * 1e3).round() / 1e3)
        ));
    }
    if !survey.interior {
        out.push(format!(
            "body '{body}' has no interior: its field is nowhere negative just inside its surface (sampled at {} mm). \
                          It is empty or a sheet of no thickness. {EMPTY_CAUSE}",
            measured(survey.h)
        ));
    }
    out.extend(undefined_warnings(s, body, &survey.undefined, survey.points));
    out
}

/// A moving body's names judged over its motion, where its survey (at the
/// session's moment) did not see one of `named` forming surface
/// (feature::surface_names_over): at moments of the file's t-span, each
/// with its own context and bounds, and proven hidden at every moment over
/// the span at once where the body reads no solve on t (one is then solved
/// moment by moment, so no one context ranges over it). The moments looked
/// at: 1 for a body that does not move.
fn over_motion(s: &mut Session, shape: u32, named: &[u32], survey: &mut ink::feature::Survey) -> usize {
    let Some(span) = s.doc.t_span.filter(|sp| sp.0 < sp.1) else { return 1 };
    let solved_on_t = s.doc.solves_read(&[shape]).iter().any(|&i| s.doc.solves[i].on_t);
    if !(solved_on_t || s.doc.g.depends_on(shape, |v| v == Var::T)) || named.iter().all(|t| survey.formed.contains(t)) {
        return 1;
    }
    let tape = Arc::new(Tape::new(&s.doc.g, shape));
    let over = match solved_on_t {
        true => None,
        false => {
            let c = ink::eval::Ctx { t_span: Some(span), ..s.ctx.clone() };
            s.at_point(c.clone(), |s| s.bounds(shape)).map(|(lo, hi, _)| (c, lo, hi))
        }
    };
    let at = |t: f64| {
        let c = s.doc.ctx(&s.overrides, t).ok()?;
        let (lo, hi, _) = s.at_point(c.clone(), |s| s.bounds(shape))?;
        Some((c, lo, hi))
    };
    let first = (std::mem::take(&mut survey.formed), std::mem::take(&mut survey.unseen));
    let (formed, unseen, moments) = ink::feature::surface_names_over(&tape, 128.0, span, named, first, at, over.as_ref().map(|(c, lo, hi)| (c, *lo, *hi)));
    (survey.formed, survey.unseen) = (formed, unseen);
    moments
}

/// Undefined cases (division by zero, ln of a negative, ...) the body's
/// field passes through (feature::undefined). Their value is 0
/// (docs/CORE.md), and in a shape 0 means "on the surface": say so, with
/// where in the model it happens.
fn undefined_warnings(s: &Session, body: &str, cases: &[ink::feature::Case], points: usize) -> Vec<String> {
    let names = &s.doc.g.tag_names;
    cases
        .iter()
        .map(|u| {
            let place = u.name.map_or(String::new(), |t| format!(" in '{}'", names[t as usize]));
            format!(
                "body '{body}': {}{place} at {} of {} sample points (e.g. at {:?}). An undefined case gives 0, which in a shape means \"on the surface\" (`ink spec`). Guard it, e.g. divide by (max d 1e-9)",
                u.what,
                u.count,
                points,
                u.at.map(|v| (v * 1e3).round() / 1e3)
            )
        })
        .collect()
}

/// What `ink require` and `ink interfere` on this file (at its --param and
/// --t) would answer at one moment though the question reads t, by the
/// rule their answers record (claims::one_moment_of): check's own view,
/// once it has answered, scoped as theirs (claim_scope), each claim asked
/// as assert asks it and, but with --load (loading only), interfere's
/// plan of each target a cheap test finds reading t (interfere::at_one_moment).
/// Nothing where a moment is asked, a span covers every question, or the
/// model reads no t.
fn at_one_moment(s: &mut Session, a: &Args, load: bool) -> Vec<String> {
    let claims = s.doc.asserts.iter().any(|x| s.doc.span_of(x.scope).is_none());
    let pairs = !load && s.doc.spans_hull().is_none();
    // The file's own requirements of ratings are answered over its span, and with none at one moment of a
    // model that moves: a motor's torque under its curve there says nothing of the rest of the motion.
    let rated: Vec<String> = match s.doc.t_span.is_none() && a.get("t-span").is_none() {
        true => s.doc.requires.iter().filter(|r| r.within.is_some() && s.doc.is_own_file(&r.file)).map(|r| format!("requirement {}", r.name)).collect(),
        false => Vec::new(),
    };
    if a.get("t").is_some() || !(claims || pairs || !rated.is_empty()) || !ink::claims::model_reads_t(s) {
        return Vec::new();
    }
    let Ok(scope) = ink::lang::command::claim_scope(s, a) else { return Vec::new() };
    if s.set_claim_scope(&scope).is_err() {
        return Vec::new();
    }
    let mut out = Vec::new();
    if claims {
        out.extend(ink::claims::claims_at_one_moment(s).unwrap_or_default());
    }
    out.extend(rated);
    if pairs {
        out.extend(ink::lang::interfere::at_one_moment(s));
    }
    out
}

/// A warning's list of names: a few, then how many more (an answer stays
/// short, however many there are).
fn few(names: &[String]) -> String {
    const SHOWN: usize = 3;
    match names.len() > SHOWN + 1 {
        true => format!("{} and {} more", names[..SHOWN].join(", "), names.len() - SHOWN),
        false => names.join(", "),
    }
}

/// The fillets and chamfers of two shapes in `f` given no cap (std's fillet-union and chamfer-union:
/// a select of the corner's construction, named "fillet" or "chamfer", and min(a, b), whose
/// condition is corner-reach's alone, not a max with a cap): each its kind and its two shapes.
fn plain_corners(g: &ink::graph::Graph, f: u32) -> Vec<(&'static str, u32, u32)> {
    use ink::graph::{BinOp, Node};
    let named = |n: u32, want: &str| matches!(g.nodes[n as usize], Node::Tag(t, _) if g.tag_names[t as usize] == want);
    let mut out = Vec::new();
    for n in g.reachable(&[f]) {
        let Node::Select(c, x, y) = g.nodes[n as usize] else { continue };
        let Node::Bin(BinOp::Min, a, b) = g.nodes[y as usize] else { continue };
        if matches!(g.nodes[c as usize], Node::Bin(BinOp::Max, ..)) {
            continue;
        }
        let kind = match g.nodes[x as usize] {
            _ if named(x, "fillet") => "fillet",
            Node::Bin(BinOp::Min, cut, sharp) if sharp == y && named(cut, "chamfer") => "chamfer",
            _ => continue,
        };
        if !out.contains(&(kind, a, b)) {
            out.push((kind, a, b));
        }
    }
    out
}

pub(crate) fn cmd_check(s: &mut Session, a: &Args) -> Result<J, InkError> {
    // --load: what loading decides, without the survey below (bounds,
    // surface sampling, solves over tolerance bands): a quick look at a
    // large model, and the catalog's check of one model per template.
    let load = a.flag("load")?;
    // What the file's imports define is listed in brief (`imports`), and
    // in full with --imports: an agent checking its design reads its own
    // names, not a catalog part's internals.
    let all = a.flag("imports")?;
    // --name: one def's value, whole (the list below shortens long ones).
    if let Some(name) = a.get("name") {
        let Some(v) = s.doc.def(name).cloned() else {
            let defs = s.doc.defs.iter().map(|(k, _)| k.as_str()).filter(|k| !k.contains(':'));
            return Err(InkError { hint: Some(ink::lang::choices(name, "its defs", defs)), ..InkError::msg(format!("no def named '{name}'")) });
        };
        // What it was found from, as an answer reading it says: sampled
        // data, solved values; of every field its value holds.
        let mut roots = Vec::new();
        v.fields_into(&mut roots, false);
        let mut out = json!({"name": name, "value": s.describe(&v)});
        inexact_fields(&mut out, &s.doc.inexact(&roots, &s.certified()), &param_names(s));
        return Ok(out);
    }
    // Fields that interval arithmetic cannot bound degrade everything downstream
    // (meshing resolution, proofs, rendering); flag them early.
    let shapes: Vec<(String, u32)> = if load { Vec::new() } else { s.doc.bodies.iter().map(|b| (b.name.clone(), b.shape)).collect() };
    let mut warnings = s.doc.load_warnings.clone();
    // Each body's bounds (the view keeps them), then their surveys, each
    // all at once.
    s.bounds_all(&shapes.iter().map(|b| b.1).collect::<Vec<_>>());
    let mut looks = Vec::new();
    for (name, shape) in shapes.iter().cloned() {
        let (whole, b) = (s.whole(shape), s.bounds(shape));
        let named = if whole && b.is_some() { s.doc.named_parts(&name, shape) } else { Vec::new() };
        looks.push((name, whole, b, Arc::new(Tape::new(&s.doc.g, shape)), named));
    }
    let ctx = &s.ctx;
    let found: Vec<_> = looks
        .par_iter()
        .map(|(_, whole, b, tape, named)| {
            b.map(|(lo, hi, _)| match whole {
                true => Ok(ink::feature::survey(tape, ctx, lo, hi, 128.0, !named.is_empty())),
                // Undefined cases are reported whatever the bounds.
                false => Err(ink::feature::undefined(tape, ctx, &[], lo, hi)),
            })
        })
        .collect();
    // A field 0 throughout a volume reads as outside there (a shape is the interior of the closure of
    // where its field is below 0): only a min or max with a 0 can leave one (a cut written as (min 0
    // g)), so only a body with one is walked for it.
    for (name, shape) in &shapes {
        let zero_operand = s.doc.g.reachable(&[*shape]).into_iter().any(|n| {
            matches!(s.doc.g.nodes[n as usize], ink::graph::Node::Bin(ink::graph::BinOp::Min | ink::graph::BinOp::Max, a, b)
                if s.doc.g.constant_value(a) == Some(0.0) || s.doc.g.constant_value(b) == Some(0.0))
        });
        let (lo, hi) = ink::lang::DEFAULT_DOMAIN;
        if zero_operand && let Some(p) = ink::feature::zero_volume(&Arc::new(Tape::new(&s.doc.g, *shape)), &s.ctx, lo, hi, 2_000_000) {
            warnings.push(format!(
                "body '{name}': its field is 0 throughout a volume about ({:.6}, {:.6}, {:.6}), which reads as outside it (a shape is the interior of the \
                 closure of where its field is below 0). A cut that is 0 where it is empty, as (subtract a (min 0 g)) where g > 0, removes what it covers: \
                 write the cut positive there",
                p[0], p[1], p[2]
            ));
        }
    }
    // A fillet or chamfer of two shapes that end at one extreme on an axis (caps in one plane: two
    // plates of one thickness, a rib flush with its wall's top) with no cap given: just beyond the
    // caps each reads distance to its cap's edge, and the corner can be filled past them.
    for (name, shape) in &shapes {
        for (kind, a, b) in plain_corners(&s.doc.g, *shape) {
            let (Some((alo, ahi, ca)), Some((blo, bhi, cb))) = (s.whole(a).then(|| s.bounds(a)).flatten(), s.whole(b).then(|| s.bounds(b)).flatten()) else {
                continue;
            };
            let near = ca.max(cb);
            let shared = (0..3).flat_map(|k| [(k, ahi[k], bhi[k]), (k, alo[k], blo[k])]).find(|&(_, p, q)| (p - q).abs() <= near);
            if let Some((k, p, q)) = shared {
                let (axis, at) = (["x", "y", "z"][k], 0.5 * (p + q));
                let side = if (p - alo[k]).abs() < (p - ahi[k]).abs() { "-" } else { "+" };
                warnings.push(format!(
                    "body '{name}': a {kind}-union joins two shapes whose proven bounds both end near {axis} = {at:.3} (within {near:.3} mm, facing {side}{axis}): where \
                     they end in caps in one plane there, just beyond them each reads distance to its cap's edge, and the {kind} can add material past the caps, up \
                     to its radius. Give the caps' plane as its last argument, below 0 on the material's side (for caps at {axis} = H facing {side}{axis}: \
                     ({kind}-union a b r ... {})), or {kind} the 2D section and extrude it",
                    if side == "+" { format!("(- {axis} H)") } else { format!("(- H {axis})") }
                ));
            }
        }
    }
    let mut surveyed = Vec::new();
    for (((name, whole, _, _, named), found), shape) in looks.into_iter().zip(found).zip(shapes.iter().map(|b| b.1)) {
        if !whole {
            warnings.push(format!(
                "body '{name}' is not proven bounded. Either it is unbounded, or interval arithmetic cannot bound it \
                 (products or selects of signs, far-reaching extrusions). ink examines it only in the 1 m cube at the origin, \
                 so results needing all of it are unknown. Intersect it with a shape containing it, e.g. (max f (box ...))"
            ));
        }
        match found {
            None if whole => warnings.push(format!("body '{name}' is empty: its field is never negative. {EMPTY_CAUSE}")),
            Some(Ok(mut survey)) => {
                let moments = over_motion(s, shape, &named, &mut survey);
                warnings.extend(survey_warnings(s, &name, named, &survey, moments));
                surveyed.push((shape, survey.sketch));
            }
            Some(Err((cases, points))) => warnings.extend(undefined_warnings(s, &name, &cases, points)),
            None => {}
        }
    }
    // Every published mass that cannot be measured (lighter than its
    // own-density bodies, over an unbounded body) is refused where it is
    // weighed, so say so here (without weighing it). Not with --load, nor
    // in other commands.
    let refused_masses = if load { Vec::new() } else { s.mass_checks() };
    warnings.extend(refused_masses.iter().filter(|(_, impossible, _)| *impossible).map(|(_, _, e)| format!("{}:{}: {}", e.file, e.line, e.message)));
    // A tolerance a proof cannot cover (its value is used at load time, or
    // a solve over it cannot be certified) is refused by assert and
    // interfere: say so here. Solves over the bands are certified once.
    let refused = warnings.len();
    for (k, p) in s.doc.params.iter().enumerate().filter(|&(k, _)| s.doc.has_range(k)) {
        if let Some(e) = s.cannot_range(k) {
            warnings.push(format!("param '{}': {}{}", p.name, e.message, e.hint.map(|h| format!(" ({h})")).unwrap_or_default()));
        }
    }
    let mut over_bands = ink::eval::ParamSpans::NONE;
    if !load && warnings.len() == refused && (0..s.doc.params.len()).any(|k| s.doc.has_range(k)) {
        match s.declared_spans().and_then(|spans| s.scope_spans(&spans)) {
            Ok(ps) => over_bands = ps,
            Err(e) => warnings.push(format!("tolerances: {}{}", e.message, e.hint.map(|h| format!(" ({h})")).unwrap_or_default())),
        }
    }
    // What each claim will cost assert, estimated (ink::cost): not with
    // --load. Own claims, and imported ones listed in full (--imports).
    let costs: Vec<Option<ink::cost::Cost>> = match load {
        true => vec![None; s.doc.asserts.len()],
        false => {
            let res = a.cell("res")?;
            // What the surveys found of the bodies' surfaces, and nothing
            // more: an estimate samples nothing check does not.
            let mut sk = ink::cost::Sketches::default();
            for (shape, sketch) in surveyed {
                sk.seed(shape, sketch);
            }
            let listed: Vec<bool> = s.doc.asserts.iter().map(|x| all || s.doc.is_own(x)).collect();
            listed.into_iter().enumerate().map(|(k, l)| if l { ink::cost::of_assert(s, &sk, k, res) } else { None }).collect()
        }
    };
    let samples: Vec<std::sync::Arc<ink::grid::Source>> = s.doc.g.grids.iter().filter_map(|g| g.source().cloned()).collect();
    let d = &s.doc;
    // Each imported file in brief: its bodies and assemblies, each with its
    // frames' names, and how many params, defs, fns and claims it has.
    let mut imports: Vec<(&str, serde_json::Map<String, J>)> = d.imported_files().into_iter().map(|f| (f, serde_json::Map::new())).collect();
    let mut own = |file: Option<&str>, kind: &str, item: Option<(&str, &[(String, Xf)])>| -> bool {
        let Some((_, m)) = file.and_then(|f| imports.iter_mut().find(|(i, _)| *i == f)) else { return true };
        match item {
            Some((name, frames)) => {
                m.entry(kind).or_insert(json!({}))[name] = json!(frames.iter().map(|f| &f.0).collect::<Vec<_>>());
            }
            None => {
                let c = m.entry(kind).or_insert(json!(0));
                *c = json!(c.as_u64().unwrap_or(0) + 1);
            }
        }
        all
    };
    // What a solve that fails here sets, or one resting on it, has no
    // value at this moment (a rate at a dead point, a linkage that closes
    // only later, and their rates).
    let failed = s.unsolved_params();
    let params: Vec<J> = d
        .params
        .iter()
        .zip(&s.ctx.params)
        .enumerate()
        // An allowance is the construction's, not a param of the model's (Doc::written_params).
        .filter(|(k, (p, _))| !d.allowances.iter().any(|a| a.0 == *k) && own(d.written_in(&p.name), "params", None))
        .map(|(k, (p, v))| {
            // A published mass's density has no value until a mass is read,
            // an actuator's effort until a claim reads it.
            let mut j = match d.unmeasured_param(k) {
                _ if failed.contains(&k) => json!({"name": p.name, "default": n(p.default), "value": null}),
                true if d.published.iter().any(|(j, _)| *j == k) => {
                    json!({"name": p.name, "default": null, "value": null, "derived": "measured where a mass is read (ink mass, a claim over (attr density))"})
                }
                true => {
                    json!({"name": p.name, "default": null, "value": null, "derived": format!("provided where it is read: {}", d.providers().owner(k).map_or("", |p| p.name()))})
                }
                false => json!({"name": p.name, "default": n(p.default), "value": n(*v)}),
            };
            if let Some((lo, hi)) = p.tol {
                j["tol"] = bounds(lo, hi);
                j["band"] = p.band(*v).map_or(J::Null, |(a, b)| bounds(a, b));
            }
            if let Some((_, (lo, hi))) = d.ranged.iter().find(|r| r.0 == k) {
                j["ranged"] = bounds(*lo, *hi);
            }
            if let Some((written, at)) = d.set_at.get(&k) {
                (j["written"], j["set_at"]) = (n(*written), json!(at));
            }
            if let Some((_, c)) = d.cases.iter().find(|c| c.0 == p.name) {
                j["cases"] = json!(c.iter().map(|g| n(g.value)).collect::<Vec<_>>());
            }
            j
        })
        .collect();
    // Text params after the number params: each value, and what values it
    // allows (a client asking a customer for one checks them first).
    let params: Vec<J> = params
        .into_iter()
        .chain(d.texts.iter().filter(|p| own(d.written_in(&p.name), "params", None)).map(|p| {
            let mut j = json!({"name": p.name, "default": p.default, "value": p.value});
            if let Some(c) = &p.chars {
                j["chars"] = json!(c);
            }
            if let Some((lo, hi)) = p.len {
                j["len"] = json!([lo, hi]);
            }
            if let Some(at) = &p.set_at {
                j["set_at"] = json!(at);
            }
            if let Some((_, c)) = d.cases.iter().find(|c| c.0 == p.name) {
                j["cases"] = json!(c.iter().map(|g| g.text.as_deref().unwrap_or_default()).collect::<Vec<_>>());
            }
            j
        }))
        .collect();
    // A def holding a function is called as a fn is: listed with them, by
    // its params (a builtin's are not written: `...`).
    let signature = |ps: &[String], rest: bool| ps.join(" ") + if rest { " ..." } else { "" };
    let called = |v: &Value| match v {
        Value::Fn(Callable::User(_, f)) => Some(signature(&f.params, f.rest)),
        Value::Fn(Callable::Closure(c)) => Some(signature(&c.params, c.rest)),
        Value::Fn(Callable::Builtin(_)) => Some("...".into()),
        _ => None,
    };
    let mut fns: Vec<(&String, String)> = Vec::new();
    let mut defs: Vec<J> = Vec::new();
    for (i, (k, v)) in d.defs.iter().enumerate() {
        if d.std_def_index(k) == Some(i) || d.texts.iter().any(|p| p.name == *k) {
            continue;
        }
        match called(v) {
            Some(p) if own(d.written_in(k), "fns", None) => fns.push((k, p)),
            None if own(d.written_in(k), "defs", None) => defs.push(json!({"name": k, "value": brief(s.describe(v))})),
            _ => {}
        }
    }
    fns.extend(d.fns.iter().filter(|(_, f)| f.scope != ink::lang::Scope::STD && own(Some(&f.file), "fns", None)).map(|(k, f)| (k, signature(&f.params, f.rest))));
    fns.sort();
    let fns: Vec<J> = fns.into_iter().map(|(k, p)| json!(format!("({k} {p})").replace(" )", ")"))).collect();
    // Declared frames, in the coordinates of what declares them.
    let frames = |fs: &[(String, Xf)]| -> J { J::Array(fs.iter().map(|(n, xf)| json!({"name": n, "frame": s.frame_json(xf)})).collect()) };
    let bodies: Vec<J> = d
        .bodies
        .iter()
        .filter(|b| own(Some(&b.file), "bodies", Some((&b.name, &b.frames))))
        .map(|b| {
            let attrs: serde_json::Map<String, J> = b.attrs.iter().map(|(k, v)| (k.clone(), brief(s.describe(v)))).collect();
            let mut j = json!({"name": b.name, "attrs": attrs});
            let read = d.inexact(&[b.shape], &[]).samples;
            if !read.is_empty() {
                j["samples"] = json!(read.iter().map(|r| &r.what).collect::<Vec<_>>());
            }
            // As an assembly's: its arguments' values as written, and what
            // it is an instance of.
            if !b.params.is_empty() {
                j["params"] = J::Array(b.params.iter().map(|(p, v)| json!({"name": p, "value": v})).collect());
            }
            if let Some(of) = &b.of {
                j["of"] = json!(of);
            }
            if !b.frames.is_empty() {
                j["frames"] = frames(&b.frames);
            }
            j
        })
        .collect();
    let assemblies: Vec<J> = d
        .assemblies
        .iter()
        .filter(|asm| own(Some(&asm.file), "assemblies", Some((&asm.name, &asm.frames))))
        .map(|asm| {
            let insts: Vec<J> = asm
                .instances
                .iter()
                .map(|i| {
                    let r: Vec<J> = i.at.r.iter().map(|row| json!(row.map(|e| n(d.scalar(e, &s.ctx))))).collect();
                    let kind = if d.body(&i.child).is_some() { "body" } else { "assembly" };
                    let mut j = json!({"name": i.name, "child": i.child, "kind": kind,
                           "rotation": r, "translation": json!(i.at.t.map(|e| n(d.scalar(e, &s.ctx))))});
                    if let Some(jt) = &i.joint {
                        let mut params: Vec<&str> = d.params_read(&[jt.drive]).into_iter().map(|k| d.params[k].name.as_str()).collect();
                        params.sort_unstable();
                        params.dedup();
                        // A turn's drive is degrees (its limits deg/s, N mm);
                        // a slide's, mm (mm/s, N).
                        let mut joint = json!({"kind": (if jt.turn { "turn" } else { "slide" }), "axis": (["x", "y", "z"][jt.axis]),
                            "anchor": s.frame_json(&jt.anchor), "drive": jt.text, "value": n(d.scalar(jt.drive, &s.ctx)),
                            "params": params, "reads_t": d.g.depends_on(jt.drive, |v| matches!(v, ink::graph::Var::T))});
                        if let Some((lo, hi)) = jt.limits {
                            joint["limits"] = bounds(d.scalar(lo, &s.ctx), d.scalar(hi, &s.ctx));
                            joint["claim"] = json!(format!("limits:{}/{}", asm.name, i.name));
                        }
                        if let Some((lo, hi)) = jt.speed_limits {
                            joint["speed_limits"] = bounds(d.scalar(lo, &s.ctx), d.scalar(hi, &s.ctx));
                            joint["speed_claim"] = json!(format!("speed-limits:{}/{}", asm.name, i.name));
                        }
                        // The heads of its packages' clauses (an actuator's
                        // effort limits), with their packages.
                        if !jt.clauses.is_empty() {
                            joint["clauses"] = jt.clauses.iter().map(|c| json!({"head": c.head, "package": c.package})).collect();
                        }
                        j["joint"] = joint;
                    }
                    j
                })
                .collect();
            let mut j = json!({"name": asm.name, "instances": insts});
            // Its parameters, with the values it was built with as written:
            // the defaults, or (built for an instance) the arguments passed.
            if !asm.params.is_empty() {
                j["params"] = J::Array(asm.params.iter().map(|(p, v)| json!({"name": p, "value": v})).collect());
            }
            if let Some(of) = &asm.of {
                j["of"] = json!(of);
            }
            // Contacts are relations with a zone alone; `relate` gives a kind.
            let (related, contacts): (Vec<_>, Vec<_>) = asm.relations.iter().partition(|r| r.kind.is_some());
            if !contacts.is_empty() {
                j["contacts"] = J::Array(contacts.iter().map(|c| json!([c.a, c.b])).collect());
            }
            if !related.is_empty() {
                j["relations"] = J::Array(
                    related
                        .iter()
                        .map(|r| {
                            let mut k = s.kind_json(r.kind.as_ref().unwrap());
                            k["pair"] = json!([r.a, r.b]);
                            k
                        })
                        .collect(),
                );
            }
            if !asm.frames.is_empty() {
                j["frames"] = frames(&asm.frames);
            }
            // Its published mass, and whether the bodies it places weigh
            // it (at the configuration it is published for), or why not.
            if let Some(m) = &asm.mass {
                let mut o = json!({"published": n(d.scalar(m.mass, &s.ctx))});
                match m.density {
                    None => (o["met"], o["why"]) = (J::Null, json!(ink::lang::materials::UNUSED_MASS)),
                    Some(_) if load => (o["met"], o["why"]) = (J::Null, json!("checked by ink check, not with --load")),
                    Some(f) => {
                        // The derived param the density is, measured where
                        // a mass is read; where it was (another mass reads
                        // it), how it settled.
                        let k = d.params_read(&[f])[0];
                        match refused_masses.iter().find(|(owns, _, _)| owns.contains(&k)) {
                            Some((_, impossible, e)) => (o["met"], o["why"]) = (if *impossible { json!(false) } else { J::Null }, json!(e.message)),
                            None => {
                                (o["met"], o["density"]) = (json!(true), json!(d.params[k].name));
                                let unsettled = if d.unmeasured(&[f]).is_empty() { d.inexact(&[f], &[]).unsettled } else { Default::default() };
                                if !unsettled.is_empty() {
                                    o["unsettled"] = J::Object(unsettled.into_iter().map(|(k, why)| (d.params[k as usize].name.clone(), json!(why))).collect());
                                }
                            }
                        }
                    }
                }
                j["mass"] = o;
            }
            j
        })
        .collect();
    // Every requirement: those of ratings, then the claims (those a limit clause or a relation makes among them).
    let rated = d
        .requires
        .iter()
        .filter_map(|r| Some((r, r.within.as_ref()?)))
        .filter(|(r, _)| d.is_own_file(&r.file) || own(Some(&r.file), "requirements", None))
        .map(|(r, within)| json!({"name": r.name, "claim": format!("(within-ratings {within})"), "source": format!("{}:{}", r.file, r.line), "own": d.is_own_file(&r.file)}))
        .collect::<Vec<J>>();
    let asserts = d.asserts.iter().zip(&costs).filter(|(x, _)| d.is_own(x) || own(Some(&x.file), "requirements", None)).map(|(x, c)| {
        let mut j = json!({"name": x.name, "claim": x.text, "source": format!("{}:{}", x.file, x.line), "own": d.is_own(x)});
        if let Some(c) = c {
            j["cost_estimate"] = c.json();
        }
        j
    });
    let asserts: Vec<J> = rated.into_iter().chain(asserts).collect();
    let imports: Vec<J> = imports
        .into_iter()
        .map(|(f, m)| {
            let mut j = serde_json::Map::from_iter([("file".to_string(), json!(f))]);
            j.extend(m);
            J::Object(j)
        })
        .collect();
    // Conditions builders stated (check): those failing here, at these
    // params and this moment.
    // In the context proofs use: over the ranges and the solves over them.
    let ranges = (0..s.doc.params.len()).any(|k| s.doc.has_range(k));
    let scope = if ranges { ink::eval::Ctx { p_spans: over_bands.clone(), ..s.ctx.clone() } } else { s.ctx.clone() };
    warnings.extend(s.doc.check_failures(&scope, s.accuracy(), &s.unsolved.iter().map(|u| u.0).collect::<Vec<_>>()));
    // A solve that does not solve at this moment (one depending on t, at a
    // moment its mechanism does not close): claims reading it are refused
    // here, so say so.
    for (_, e) in &s.unsolved {
        warnings.push(format!("{} (at {}:{}), so claims using it are refused at this moment. Check them at another with --t", e.message, e.file, e.line));
    }
    let default_target = if d.body("root").is_some() || d.assembly("root").is_some() { json!("root") } else { json!("all bodies") };
    let mut out = json!({
        "file": s.file(),
        "params": params, "defs": defs, "fns": fns,
        "bodies": bodies, "assemblies": assemblies, "requirements": asserts, "imports": imports,
        "t_span": d.t_span.map(|(a, b)| bounds(a, b)),
        "default_target": default_target,
        "samples": sources_json(&samples),
        "warnings": warnings,
    });
    // An import's span, which its claims and assemblies are checked over
    // here as in it (Doc::span_of).
    let mut spans = std::collections::BTreeMap::<_, Vec<_>>::new();
    for (file, ns, span) in d.import_spans() {
        spans.entry(file).or_default().push((ns, span));
    }
    if !spans.is_empty() {
        out["import_t_spans"] = J::Object(
            spans
                .into_iter()
                .map(|(file, imports)| {
                    let value = match imports.as_slice() {
                        [(_, (a, b))] => bounds(*a, *b),
                        _ => J::Object(imports.into_iter().map(|(ns, (a, b))| (ns.unwrap_or("").to_owned(), bounds(a, b))).collect()),
                    };
                    (file.to_owned(), value)
                })
                .collect(),
        );
    }
    if !d.cases.is_empty() {
        // Imported params can declare cases even when only local params are listed.
        out["has_cases"] = json!(true);
    }
    unchecked_fields(&mut out, s.doc.unchecked_claims()?.iter().map(|c| (c.name.as_str(), c.reason)));
    // What the design states of its parts and what each place takes: every (fact ...) under an assembly
    // nothing else places, with its value as read, and the names the ratings and formulas at each place read.
    // A fact nothing reads (a misspelled one) is a warning, in the words `ink require` says it in.
    match ink_ratings::stated_and_taken(s) {
        Ok((design, taken)) if !(design.facts.is_empty() && taken.is_empty()) => {
            let unread = ink_ratings::unread(&design, &taken);
            let facts: Vec<J> = (design.facts.iter())
                .map(
                    |f| json!({"what": f.about.what, "on": f.about.on, "value": f.value, "unit": f.unit, "from": f.from, "read": !unread.iter().any(|(u, _)| std::ptr::eq(*u, f))}),
                )
                .collect();
            let mut takes = serde_json::Map::new();
            for t in &taken {
                let names = takes.entry(t.on.clone()).or_insert(json!([])).as_array_mut().unwrap();
                if !names.contains(&json!(t.what)) {
                    names.push(json!(t.what));
                }
            }
            let said: Vec<J> = unread.iter().map(|(fact, names)| json!(ink_ratings::unread_said(fact, names))).collect();
            out["warnings"].as_array_mut().unwrap().extend(said);
            (out["facts"], out["takes"]) = (json!(facts), J::Object(takes));
        }
        Ok(_) => {}
        Err(why) => out["warnings"].as_array_mut().unwrap().push(json!(format!("its parts' ratings are not read: {why}"))),
    }
    let reg = s.doc.registry;
    for p in reg.packages {
        extend(&mut out, J::Object(p.check(s)));
    }
    // require and interfere would answer it at one moment, and holding
    // there reads as holding over the motion. Last: the view is scoped as
    // theirs to ask it.
    let unspanned = at_one_moment(s, a, load);
    if !unspanned.is_empty()
        && let Some(w) = out["warnings"].as_array_mut()
    {
        w.push(json!(format!("the file reads t but states no t-span, so require and interfere check only t = 0 ({}). Add (t-span a b)", few(&unspanned))));
    }
    Ok(out)
}

/// A value as check's overview shows it: whole when small, else what it is
/// and how big (a record's keys, a list's length): a catalog part's data
/// tables are hundreds of kB, and an assembly importing many parts made the
/// overview megabytes, past what an agent reads.
fn brief(v: J) -> J {
    const WHOLE: usize = 1024;
    let bytes = v.to_string().len();
    if bytes <= WHOLE {
        return v;
    }
    match v {
        J::Object(mut m) if m.len() == 1 && m.contains_key("data") => json!({"data": brief_of(m.remove("data").unwrap(), bytes)}),
        v => brief_of(v, bytes),
    }
}

fn brief_of(v: J, bytes: usize) -> J {
    match v {
        J::Object(m) => json!({"keys": m.keys().collect::<Vec<_>>(), "bytes": bytes}),
        J::Array(l) => json!({"items": l.len(), "bytes": bytes}),
        v => json!({"bytes": bytes, "starts": v.to_string().chars().take(80).collect::<String>()}),
    }
}
