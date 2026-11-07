//! Geometry queries: check, eval, bounds, mass, slice, mesh, components,
//! ray, distance, render.

use crate::*;
use ink::graph::Var;
use ink::lang::Callable;
use ink::lang::look::Look;
use ink::view::Local;
use rayon::prelude::*;

/// The usual cause of a body that is unexpectedly empty, for its warnings.
const EMPTY_CAUSE: &str = "(max a b) is the intersection of a and b, (min a b) their union. Shapes that only touch don't intersect";

/// Named constructions (defs, tags) in a body's field not seen forming any
/// of its surface (feature::surface_names): proven to form none (hidden
/// inside another shape, or cut away entirely), or never seen deciding a
/// patch of it where it was sought (only an edge, a blend's covered side,
/// or a face the search missed). Usually an authoring slip (a mirrored
/// half that vanished, a hole outside the part). Then fields interval
/// arithmetic cannot bound tightly, and undefined cases (feature::survey).
/// `named`: the body's named parts.
fn survey_warnings(s: &Session, body: &str, mut named: Vec<u32>, survey: &ink::feature::Survey) -> Vec<String> {
    named.retain(|t| !survey.formed.contains(t));
    let mut out: Vec<String> = named
        .into_iter()
        .map(|t| {
            let name = &s.doc.g.tag_names[t as usize];
            if let Some(&(_, n, finer)) = survey.unseen.iter().find(|u| u.0 == t) {
                format!(
                    "body '{body}': '{name}' was not seen forming any of its surface where it was sought: in the {n} cells it is on the \
                     field's path in, finer in {} of them. It may form only an edge (where it ties with what covers it) or a blend's \
                     covered side, or a face the search missed. Undecided, not proven",
                    finer
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

/// What `ink assert` and `ink interfere` on this file (at its --param and
/// --t) would answer at one moment though the question reads t, by the
/// rule their answers record (claims::one_moment_of): check's own view,
/// once it has answered, scoped as theirs (claim_scope), each claim asked
/// as assert asks it and, but with --load (loading only), interfere's
/// plan of each target a cheap test finds reading t (interfere::at_one_moment).
/// Nothing where a moment is asked, a span covers every question, or the
/// model reads no t.
fn at_one_moment(s: &mut Session, a: &Args, load: bool) -> Vec<String> {
    let claims = s.doc.asserts.iter().any(|x| s.doc.span_of(&x.file).is_none());
    let pairs = !load && s.doc.spans_hull().is_none();
    if a.get("t").is_some() || !(claims || pairs) || !ink::claims::model_reads_t(s) {
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

pub(crate) fn cmd_check(s: &mut Session, a: &Args) -> Result<J, InkError> {
    // --load: what loading decides, without the survey below (bounds,
    // surface sampling, solves over tolerance bands): a quick look at a
    // large model, and the catalog's check of one model per template.
    let load = a.flag("load")?;
    // What the file's imports define is listed in brief (`imports`), and
    // in full with --imports: an agent checking its design reads its own
    // names, not a catalog part's internals.
    let all = a.flag("imports")?;
    // Fields that interval arithmetic cannot bound degrade everything downstream
    // (meshing resolution, proofs, rendering); flag them early.
    let shapes: Vec<(String, u32)> = if load { Vec::new() } else { s.doc.bodies.iter().map(|b| (b.name.clone(), b.shape)).collect() };
    let mut warnings = s.doc.load_warnings.clone();
    // Each body's bounds (the view keeps them), then their surveys, each
    // all at once.
    s.bounds_all(&shapes.iter().map(|b| b.1).collect::<Vec<_>>());
    let mut looks = Vec::new();
    for (name, shape) in shapes {
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
    for ((name, whole, _, _, named), found) in looks.into_iter().zip(found) {
        if !whole {
            warnings.push(format!(
                "body '{name}' is not proven bounded. Either it is unbounded, or interval arithmetic cannot bound it \
                 (products or selects of signs, far-reaching extrusions). ink examines it only in the 1 m cube at the origin, \
                 so results needing all of it are unknown. Intersect it with a shape containing it, e.g. (max f (box ...))"
            ));
        }
        match found {
            None if whole => warnings.push(format!("body '{name}' is empty: its field is never negative. {EMPTY_CAUSE}")),
            Some(Ok(survey)) => warnings.extend(survey_warnings(s, &name, named, &survey)),
            Some(Err((cases, points))) => warnings.extend(undefined_warnings(s, &name, &cases, points)),
            None => {}
        }
    }
    // Every published mass that cannot be measured (lighter than its
    // own-density bodies, over an unbounded body) is refused where it is
    // weighed, so say so here (without weighing it). Not with --load, nor
    // in other commands.
    let refused_masses = if load { Vec::new() } else { s.mass_refusals() };
    warnings.extend(refused_masses.iter().map(|(_, e)| format!("{}:{}: {}", e.file, e.line, e.message)));
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
    let certified = s.certified_in(&over_bands);
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
        .filter(|(_, (p, _))| own(d.written_in(&p.name), "params", None))
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
            j
        })
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
        if d.def_index_of(&format!("std:{k}")) == Some(i) {
            continue;
        }
        match called(v) {
            Some(p) if own(d.written_in(k), "fns", None) => fns.push((k, p)),
            None if own(d.written_in(k), "defs", None) => defs.push(json!({"name": k, "value": brief(s.describe(v))})),
            _ => {}
        }
    }
    fns.extend(d.fns.iter().filter(|(k, f)| !k.starts_with("file:") && &*f.file != "std" && own(Some(&f.file), "fns", None)).map(|(k, f)| (k, signature(&f.params, f.rest))));
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
                        match refused_masses.iter().find(|(owns, _)| owns.contains(&k)) {
                            Some((_, e)) => (o["met"], o["why"]) = (json!(false), json!(e.message)),
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
    let asserts: Vec<J> = d
        .asserts
        .iter()
        .filter(|x| own(Some(&x.file), "asserts", None))
        .map(|x| json!({"name": x.name, "claim": x.text, "source": format!("{}:{}", x.file, x.line), "own": d.is_own(x)}))
        .collect();
    let imports: Vec<J> = imports
        .into_iter()
        .map(|(f, m)| {
            let mut j = serde_json::Map::from_iter([("file".to_string(), json!(f))]);
            j.extend(m);
            J::Object(j)
        })
        .collect();
    // The solves: forms whose params are a root of their equations
    // (Form::equations), in document order.
    let solves: Vec<J> = d
        .ext
        .iter()
        .filter_map(|e| Some((e.form.name(), e.form.equations()?, e.form.provider(d)?)))
        .map(|(name, residuals, p)| {
            let vals: serde_json::Map<String, J> = p.owns().iter().map(|&k| (d.params[k].name.clone(), n(s.ctx.params[k]))).collect();
            let mut j = json!({"name": name, "values": vals, "equations": residuals.len(), "uses_t": p.uses_t(), "largest_residual_at_load": n(residuals.iter().map(|&q| d.scalar(q, &s.ctx).abs()).fold(0.0, f64::max))});
            // Over the declared tolerances: each unknown's enclosure.
            let bands: serde_json::Map<String, J> =
                p.owns().iter().filter_map(|&k| Some((d.params[k].name.clone(), over_bands.get(k as u32).map(|b| bounds(b.0, b.1))?))).collect();
            if !bands.is_empty() {
                // What its equations read (Doc::inexact): a root of
                // sampled fields is certified for the samples only.
                let read = d.inexact(residuals, &certified);
                // Certified (Krawczyk) for every value of the declared
                // tolerances; for sampled fields' samples only, and not
                // where it reads an uncertified solved value (inexact_fields).
                j["over_bands"] = J::Object(bands);
                inexact_fields(&mut j, &read, &param_names(s));
            }
            j
        })
        .collect();

    // Conditions builders stated (check): those failing here, at these
    // params and this moment.
    warnings.extend(s.doc.check_failures(&s.ctx, s.accuracy(), &s.unsolved.iter().map(|u| u.0).collect::<Vec<_>>()));
    // A solve that does not solve at this moment (one depending on t, at a
    // moment its mechanism does not close): claims reading it are refused
    // here, so say so.
    for (_, e) in &s.unsolved {
        warnings.push(format!("{} (at {}:{}), so claims using it are refused at this moment. Check them at another with --t", e.message, e.file, e.line));
    }
    let default_target = if d.body("root").is_some() || d.assembly("root").is_some() { json!("root") } else { json!("all bodies") };
    let mut out = json!({
        "file": s.file(),
        "params": params, "defs": defs, "fns": fns, "solves": solves,
        "bodies": bodies, "assemblies": assemblies, "asserts": asserts, "imports": imports,
        "t_span": d.t_span.map(|(a, b)| bounds(a, b)),
        "default_target": default_target,
        "samples": sources_json(&samples),
        "warnings": warnings,
    });
    // An import's span, which its claims and assemblies are checked over
    // here as in it (Doc::span_of).
    let mut spans: Vec<(&str, (f64, f64))> = d.file_spans.iter().filter(|(f, _)| !d.is_own_file(f)).map(|(f, &s)| (&**f, s)).collect();
    if !spans.is_empty() {
        spans.sort_by(|a, b| a.0.cmp(b.0));
        out["import_t_spans"] = J::Object(spans.into_iter().map(|(f, (a, b))| (f.to_string(), bounds(a, b))).collect());
    }
    let reg = s.doc.registry;
    for p in reg.packages {
        extend(&mut out, J::Object(p.check(s)));
    }
    // assert and interfere would answer it at one moment, and holding
    // there reads as holding over the motion. Last: the view is scoped as
    // theirs to ask it.
    let unspanned = at_one_moment(s, a, load);
    if !unspanned.is_empty()
        && let Some(w) = out["warnings"].as_array_mut()
    {
        w.push(json!(format!("the file reads t but states no t-span, so assert and interfere check only t = 0 ({}). Add (t-span a b)", few(&unspanned))));
    }
    Ok(out)
}

pub(crate) fn cmd_eval(s: &mut Session, a: &Args) -> Result<J, InkError> {
    let pts = parse_points(a)?;
    let items = s.items(a)?;
    let shape = Records::new(a);
    let fields: Vec<u32> = items.iter().map(|it| s.doc.placed_field(it)).collect();
    // Along t and each param no solve sets, through solves: the partial
    // plus each solved value's partial times its rate (ModelView::rates).
    let free: Vec<usize> = (0..s.doc.params.len()).filter(|&k| !s.doc.solved_params.contains(&k)).collect();
    let along: Vec<Var> = std::iter::once(Var::T).chain(free.iter().map(|&k| Var::Param(k as u32))).collect();
    let terms: Vec<Vec<(usize, Var, usize)>> = fields.iter().map(|&f| s.rates(&[f], &along)).collect();
    // A rate not defined at this moment (a dead point) leaves no number.
    let undefined = s.unsolved_params();
    // The params listed: every one but rates.
    let listed: Vec<usize> = (0..s.doc.params.len()).filter(|&k| !s.doc.solves.iter().any(|x| x.rates.is_some() && x.owns.contains(&k))).collect();
    let mut out = Vec::new();
    for ((it, &f), terms) in items.iter().zip(&fields).zip(&terms) {
        let tape = Tape::new(&s.doc.g, f);
        let attrs = if shape.wants("attrs") { s.doc.placed_attrs(it) } else { Default::default() };
        let res: Vec<J> = pts
            .iter()
            .map(|&p| {
                let d = eval_dual(&tape, &s.ctx, p);
                let total = |var: Var, partial: f64| {
                    terms.iter().filter(|x| x.1 == var).try_fold(partial, |acc, &(k, _, y)| (!undefined.contains(&y)).then(|| acc + d[5 + k] * s.ctx.params[y]))
                };
                let mut r = json!({"point": v3(p), "value": n(d[0]), "inside": d[0] < 0.0,
                       "grad": v3([d[1], d[2], d[3]]), "d_dt": total(Var::T, d[4]).map_or(J::Null, n)});
                if shape.wants("d_dparams") {
                    r["d_dparams"] =
                        listed.iter().map(|&k| (s.doc.params[k].name.clone(), total(Var::Param(k as u32), d[5 + k]).map_or(J::Null, n))).collect::<serde_json::Map<_, _>>().into();
                }
                // Attributes evaluate their own fields: only when wanted. A
                // record (a part's datasheet, its JSON) is the same at every
                // point: only when --keys names it (ink check lists its keys).
                if shape.wants("attrs") {
                    r["attrs"] = attrs
                        .iter()
                        .filter(|(k, v)| !matches!(v, Value::Data(_)) || shape.names(&format!("attrs.{k}")))
                        .map(|(k, v)| (k.clone(), s.value_at(v, p)))
                        .collect::<serde_json::Map<_, _>>()
                        .into();
                }
                // Undefined cases this value passed through (each gave 0).
                let undefined: Vec<J> = ink::eval::undefined_cases(&tape, &s.ctx, &[p])
                    .iter()
                    .map(|u| {
                        json!(match u.name {
                            Some(t) => format!("{} in '{}'", u.what, s.doc.g.tag_names[t as usize]),
                            None => u.what.to_string(),
                        })
                    })
                    .collect();
                if !undefined.is_empty() {
                    r["undefined"] = json!(undefined);
                }
                r
            })
            .collect();
        out.push(json!({"target": it.label, "results": shape.shape(res, &["undefined"])?}));
    }
    Ok(json!({"targets": out}))
}

pub(crate) fn cmd_bounds(s: &mut Session, a: &Args) -> Result<J, InkError> {
    let items = s.items(a)?;
    let boxes: Vec<_> = items.iter().filter_map(|it| s.item_bounds(it)).map(|b| (b.0, b.1)).collect();
    let fields = s.doc.placed_reads(&items);
    ink::validity::refuse(&s.doc, &s.ctx, &fields, boxes, a.cell("res")?.unwrap_or(f64::INFINITY))?;
    let mut out = Vec::new();
    let mut all = Vec::new();
    let mut open = false;
    for it in &items {
        // A shape not proven bounded has no box: what lies in the window is
        // reported, as that.
        let whole = s.item_whole(it);
        open |= !whole;
        let mut r = match s.item_bounds(it) {
            Some((lo, hi, h)) => {
                if whole {
                    all.push((lo, hi, h));
                }
                json!({"target": it.label, "min": v3(lo), "max": v3(hi), "size": v3([hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2]]), "cell": n(h)})
            }
            None if whole => json!({"target": it.label, "empty": true}),
            None => json!({"target": it.label}),
        };
        // Not proven bounded: min and max (if any) bound the part in the
        // window only.
        if !whole {
            r["unbounded"] = json!(true);
        }
        out.push(r);
    }
    let union = if open { None } else { union_bounds(&all).map(|(lo, hi)| json!({"min": v3(lo), "max": v3(hi)})) };
    Ok(json!({"targets": out, "union": union}))
}

/// A measured answer's cells: `res`, where halved to them the mass's
/// `spread` relative to it, `coarse_res` and `relative_change`, and
/// `unsettled` where it stopped short.
fn measured_json(r: &mut J, m: &query::Measured) {
    r["res"] = n(m.res);
    if let Some(s) = m.relative_spread() {
        r["spread"] = n(s);
    }
    if let Some((coarse, change)) = m.coarse {
        r["coarse_res"] = n(coarse);
        r["relative_change"] = n(change);
    }
    if let Some(why) = m.unsettled() {
        r["unsettled"] = json!(why);
    }
}

pub(crate) fn cmd_slice(s: &mut Session, a: &Args) -> Result<J, InkError> {
    let items = s.items(a)?;
    let axis = match a.get("axis").unwrap_or("z") {
        "x" => 0,
        "y" => 1,
        "z" => 2,
        o => return Err(InkError::msg(format!("--axis must be x, y or z, got '{o}'"))),
    };
    // Each target is sliced on its own, so parts that touch (or nearly touch)
    // do not merge into one loop; the plane defaults to the middle of them all.
    let fields: Vec<u32> = items.iter().map(|it| s.doc.placed_field(it)).collect();
    let bounds: Vec<Option<(V3, V3, f64)>> = fields.iter().map(|&f| s.bounds(f)).collect();
    let all: Vec<(V3, V3, f64)> = bounds.iter().flatten().copied().collect();
    let Some((ulo, uhi)) = union_bounds(&all) else { return Ok(json!({"empty": true, "loops": []})) };
    let at = a.f64("at")?.unwrap_or((ulo[axis] + uhi[axis]) / 2.0);
    let (u, v) = query::slice_axes(axis);
    let names = ["x", "y", "z"];
    let shape = Records::new(a);
    // An outline to cut (DXF) is sliced at the model's accuracy, not a
    // picture's 1/400 of the target.
    let out = a.get("out");
    if let Some(o) = out.filter(|o| !o.ends_with(".dxf")) {
        return Err(InkError::msg(format!("--out {o}: a slice writes .dxf")));
    }
    let (mut total, mut lj, mut hmin, mut tmax) = (0.0, Vec::new(), f64::INFINITY, 0.0f64);
    // Each outline's layer, points, whether it closes and whether a mark.
    let mut cut: Vec<(String, Vec<[f64; 2]>, bool, bool)> = Vec::new();
    for ((it, &f), b) in items.iter().zip(&fields).zip(&bounds) {
        let Some((lo, hi, _)) = b else { continue };
        if at < lo[axis] || at > hi[axis] {
            continue;
        }
        let h = a.cell("res")?.unwrap_or(if out.is_some() { s.accuracy() } else { (hi[u] - lo[u]).max(hi[v] - lo[v]) / 400.0 });
        hmin = hmin.min(h);
        let tol = a.f64("accuracy")?.unwrap_or(h * 0.5);
        tmax = tmax.max(tol);
        let tape = s.clipped(f);
        for l in query::slice(&tape, &s.ctx, axis, at, *lo, *hi, h) {
            total += l.area;
            let (mut b0, mut b1) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
            for p in &l.pts {
                for k in 0..2 {
                    b0[k] = b0[k].min(p[k]);
                    b1[k] = b1[k].max(p[k]);
                }
            }
            let mut r = json!({"target": it.label, "kind": if l.area >= 0.0 {"outer"} else {"hole"}, "closed": l.closed,
                   "area": n(l.area), "perimeter": n(l.perimeter),
                   "min": json!([n(b0[0]), n(b0[1])]), "max": json!([n(b1[0]), n(b1[1])])});
            // A loop's points only when --keys names them: an outline of
            // tens of loops is tens of KB.
            if shape.names("points") || out.is_some() {
                let pts = simplify_polyline(&l.pts, l.closed, tol);
                if shape.names("points") {
                    r["points"] = pts.iter().map(|p| json!([n(p[0]), n(p[1])])).collect();
                }
                cut.push((it.label.clone(), pts, l.closed, false));
            }
            lj.push(r);
        }
        // Its body's marks (a sheet's bend lines), within its bounds:
        // open curves on their own layers, to draw, not to cut.
        let body = it.instance.as_ref().map_or(it.label.as_str(), |(b, _)| b.as_str());
        let marks = s.doc.body(body).map(|b| b.marks.clone()).unwrap_or_default();
        for (name, g, rg) in marks {
            let (g, rg) = (s.doc.place(g, it), s.doc.place(rg, it));
            let (gt, rt) = (Arc::new(Tape::new(&s.doc.g, g)), Tape::new(&s.doc.g, rg));
            for l in query::slice_mark(&gt, &rt, &s.ctx, axis, at, *lo, *hi, h) {
                let mut r = json!({"target": it.label, "kind": "mark", "mark": name, "closed": l.closed, "length": n(l.perimeter)});
                if shape.names("points") || out.is_some() {
                    let pts = simplify_polyline(&l.pts, l.closed, tol);
                    if shape.names("points") {
                        r["points"] = pts.iter().map(|p| json!([n(p[0]), n(p[1])])).collect();
                    }
                    cut.push((name.clone(), pts, l.closed, true));
                }
                lj.push(r);
            }
        }
    }
    // Each target sliced on its own; outer loops counter-clockwise
    // (positive area) in plane_axes, holes clockwise.
    let mut r = json!({"axis": names[axis], "at": n(at), "plane_axes": [names[u], names[v]], "res": n(hmin),
              "area": n(total), "loop_count": lj.len(), "loops": shape.shape(lj, &[])?});
    if let Some(o) = out {
        let outlines: Vec<export::Outline> = cut.iter().map(|(l, p, c, m)| export::Outline { layer: l, pts: p, closed: *c, mark: *m }).collect();
        let note = format!("ink {}: slice {}={} in {}{} mm, res {} accuracy {}", s.stem(), names[axis], measured(at), names[u], names[v], measured(hmin), measured(tmax));
        let hint = "DXF reads a layer's name in any case, its <>/\\\":;?*|=' and spaces as _: rename the mark, (mark NAME ...) with a NAME no target has";
        let bytes = export::dxf(&outlines, &note).map_err(|m| InkError { hint: Some(hint.into()), ..InkError::msg(m) })?;
        export::write(o, &bytes).map_err(|e| InkError::msg(format!("cannot write {o}: {e}")))?;
        r["out"] = json!(o);
        r["accuracy"] = json!(n(tmax));
    }
    Ok(r)
}

pub struct Meshed {
    pub label: String,
    /// The mesh in the body's own frame, shared by all instances of that body.
    pub mesh: Arc<Mesh>,
    pub colors: Arc<Vec<[f32; 3]>>,
    /// For seeing: the field's normal at each triangle's corners
    /// (`corner_normals`); else empty.
    pub normals: Arc<Vec<[f32; 3]>>,
    /// Where each vertex's paint (colour, metallic, roughness) is read: for
    /// seeing, under the surface (`ink::seen::Viewed::paint_at`); else the
    /// vertices.
    pub paint_at: Arc<Vec<V3>>,
    /// How it was made and how far it lies from the surface.
    pub quality: Arc<Quality>,
    /// Placement of this instance (None: the mesh is already in world space).
    pub place: Option<Rigid>,
    /// How much farther than in the body's frame the placed mesh may be
    /// from the instance's surface (ink::view::Local::off; 0 unplaced).
    pub off: f64,
}

pub struct Quality {
    /// The sampling, and the accuracy it was simplified to (0: uniform).
    pub res: f64,
    pub tol: f64,
    /// Triangles at res, before simplification.
    pub uniform_triangles: usize,
    /// The farthest a vertex at res moved from the result (<= tol).
    pub moved: f64,
    /// Measured for `ink mesh`: the sampled distance to the surface.
    pub deviation: Option<ink::simplify::Deviation>,
    /// Sign changes between samples, which the mesh leaves out.
    pub unresolved: mesh::Unresolved,
    /// How far each vertex may lie from the surface, where asked
    /// (`MeshHow::errors`); else none.
    pub error: Vec<f32>,
}

impl Meshed {
    /// The mesh in world coordinates.
    fn world(&self) -> Mesh {
        match &self.place {
            None => (*self.mesh).clone(),
            Some(rg) => Mesh { verts: self.mesh.verts.iter().map(|&v| rg.apply(v)).collect(), tris: self.mesh.tris.clone() },
        }
    }
}

/// How `mesh_items` meshes: cells across a target without --res; the
/// accuracy as a fraction of res without --accuracy (0: uniform); whether to
/// measure the distance to the surface (`measure`).
pub struct MeshHow {
    pub div: f64,
    pub tol_per_res: f64,
    pub measure: bool,
    /// What the meshes may cost in all, in triangles sampled (`work`), each
    /// as often as it is made (every one `drawn` makes of it), whether made
    /// or found in the cache, and a flipbook's frames as predicted: meshes
    /// are made again only within it, so they are the model's alone, not
    /// its host's, threads', load's or cache's.
    pub budget: f64,
    /// The finest res a mesh is made again at (0: any).
    pub finest: f64,
    /// Each vertex's error is wanted (render ties two objects' surfaces by
    /// it: render::Object).
    pub errors: bool,
    /// For seeing, not measuring: a body is meshed once and placed by its
    /// placement's numbers even where they are not proven exact (a turn by
    /// a cosine of t), so a moving part keeps one mesh.
    pub seen: bool,
}

/// `ink mesh`'s default accuracy, a fraction of res: 0.01 mm at res 0.3.
const EXPORT_TOL_PER_RES: f64 = 1.0 / 30.0;

/// Mesh every target (`plan_meshes`, then `make_meshes`).
fn mesh_items(s: &mut Session, a: &Args, items: &[Item], how: MeshHow) -> Result<(Vec<Meshed>, Option<Budget>), InkError> {
    let plan = plan_meshes(s, a, items, &how)?;
    make_meshes(s, a, plan, &how)
}

/// What `make_meshes` makes: each distinct mesh once (`jobs`: its item,
/// holding the field meshed, and its first target's place in the items),
/// its box, res and bend (`boxes`: triangles kept short along x, for a
/// strip a client bends round arcs as tight as it), the times it is made
/// (`times`: 1, or a client's frames of it, which the budget predicts), and
/// each target's mesh and placement (`targets`: its label, mesh, rigid
/// placement, and its instance's frame with its mesh's box). A client may
/// change what a mesh is (its field, box and res) before it is made.
pub struct Plan {
    pub jobs: Vec<(Item, usize)>,
    pub boxes: Vec<(V3, V3, f64, Option<f64>)>,
    pub times: Vec<f64>,
    pub targets: Vec<(String, usize, Option<Rigid>, Option<(Local, V3, V3)>)>,
    /// How many frames a client makes of the meshes it makes more than
    /// once (`times` over 1), from the triangles they sample and keep at
    /// one moment: what the budget predicts for them.
    pub frames: Option<Box<dyn Fn(usize, usize) -> usize>>,
    /// --res, and --max-res (infinite without).
    pub res: Option<f64>,
    pub coarsest: f64,
    div: f64,
    floor: f64,
    seen: bool,
}

impl Plan {
    /// The res a box of its own gives a mesh without --res.
    pub fn own(&self, lo: V3, hi: V3) -> f64 {
        let h = match self.seen {
            true => (box_area(lo, hi) / (6.0 * self.div * self.div)).sqrt(),
            false => (max_extent(lo, hi) / self.div).max(self.floor),
        };
        h.min(self.coarsest)
    }
}

/// What to mesh for the targets. Without --res, each object gets its own
/// resolution (its size / `how.div`, but no finer than the whole scene's
/// size / (8 * how.div), so small parts in huge assemblies stay
/// affordable; and no coarser than --max-res). Instances of one body share
/// a single mesh made in the body's frame (ModelView::local), at the res of
/// its size there, however each is turned. For seeing (`how.seen`) a
/// body's own res gives its box as many cells over its surface as a cube
/// `how.div` cells across (so a long part is not coarse across), and a
/// placement that reads t places the body at each motion sample too (the
/// GLB), so it must be rigid at each (ModelView::rigid).
pub fn plan_meshes(s: &mut Session, a: &Args, items: &[Item], how: &MeshHow) -> Result<Plan, InkError> {
    let res = a.cell("res")?;
    let coarsest = a.cell("max-res")?.unwrap_or(f64::INFINITY);
    let tol = a.f64("accuracy")?;
    if tol.is_some_and(|t| !(t >= 0.0 && t.is_finite())) {
        return Err(InkError::msg("--accuracy must be 0 (uniform) or positive and finite"));
    }
    let moves = |s: &Session, it: &Item| it.instance.as_ref().is_some_and(|(_, xf)| ink::claims::reads_t(s, &xf.r.iter().flatten().chain(&xf.t).copied().collect::<Vec<_>>()));
    let moments = match how.seen && items.iter().any(|it| moves(s, it)) {
        true => motion(s)?.1,
        false => vec![],
    };
    let locals: Vec<Option<(u32, Rigid, Option<Local>)>> = items
        .iter()
        .map(|it| {
            let (body, rg) = s.rigid(it, &moments)?;
            match s.local(it) {
                Some(l) => Some((l.shape, l.at, Some(l))),
                None if how.seen => Some((s.doc.body(&body)?.shape, rg, None)),
                None => None,
            }
        })
        .collect();
    let bounds: Vec<Option<(V3, V3, f64)>> = items
        .iter()
        .zip(&locals)
        .map(|(it, l)| match l {
            Some((shape, rg, _)) if how.seen => s.bounds(*shape).map(|(lo, hi, h)| {
                let (a, b) = rg.aabb(lo, hi);
                (a, b, h)
            }),
            _ => s.item_bounds(it),
        })
        .collect();
    let all: Vec<(V3, V3, f64)> = bounds.iter().flatten().copied().collect();
    let Some((ulo, uhi)) = union_bounds(&all) else { return Err(InkError::msg("nothing to mesh: every target is empty")) };
    let floor = if how.seen { 0.0 } else { max_extent(ulo, uhi) / (8.0 * how.div) };
    let mut plan = Plan { jobs: Vec::new(), boxes: Vec::new(), times: Vec::new(), targets: Vec::new(), frames: None, res, coarsest, div: how.div, floor, seen: how.seen };
    let mut shared: HashMap<(String, u64), usize> = HashMap::new();
    for (k, ((it, b), l)) in items.iter().zip(&bounds).zip(&locals).enumerate() {
        let Some((lo, hi, _)) = *b else { continue };
        let h = res.unwrap_or(plan.own(lo, hi));
        let cells = max_extent(lo, hi) / h;
        if cells > 20000.0 {
            return Err(InkError::msg(format!("--res {h} gives {cells:.0} cells across '{}'; the limit is 20000", it.label)));
        }
        // Instanced: mesh the body once in its own frame.
        if let Some((shape, rg, local)) = *l {
            let body = it.instance.as_ref().unwrap().0.clone();
            let Some((llo, lhi, _)) = s.bounds(shape) else { continue };
            let h = res.unwrap_or(plan.own(llo, lhi));
            let j = *shared.entry((body.clone(), h.to_bits())).or_insert_with(|| {
                plan.jobs.push((Item { label: body.clone(), shape, attrs: it.attrs.clone(), instance: None, links: Vec::new() }, k));
                plan.boxes.push((llo, lhi, h, None));
                plan.jobs.len() - 1
            });
            plan.targets.push((it.label.clone(), j, Some(rg), local.map(|l| (l, llo, lhi))));
            continue;
        }
        plan.jobs.push((Item { label: it.label.clone(), shape: s.doc.placed_field(it), attrs: s.doc.placed_attrs(it), instance: None, links: it.links.clone() }, k));
        plan.boxes.push((lo, hi, h, None));
        plan.targets.push((it.label.clone(), plan.jobs.len() - 1, None, None));
    }
    if plan.targets.is_empty() {
        return Err(InkError::msg("nothing to mesh: every target is empty"));
    }
    plan.times = vec![1.0; plan.jobs.len()];
    Ok(plan)
}

/// Make the meshes planned (`plan_meshes`). For seeing (`how.seen`, no
/// --res) a mesh is no finer than the model's floor times the cube root of
/// the times it is made (`view_floor`: its meshes within VIEW_TRIANGLES,
/// as their surfaces' areas predict them, ink::seen::surface), on the
/// ladder (`ladder`); one whose surface lies farther from it on average
/// than its res is made finer whatever the budget (`drawn`). For seeing, a
/// mesh that leaves material out is made again finer, and kept where it
/// draws more, within `how.budget` (not one at --max-res). Each mesh is
/// reported as it is made (ink::claims::progress); the budget, as planned
/// and made, is returned for `cost_warnings`.
pub fn make_meshes(s: &mut Session, a: &Args, plan: Plan, how: &MeshHow) -> Result<(Vec<Meshed>, Option<Budget>), InkError> {
    let (res, coarsest, tol) = (plan.res, plan.coarsest, a.f64("accuracy")?);
    let Plan { jobs, mut boxes, times, targets, frames, .. } = plan;
    // For seeing, each mesh's surface (ink::seen::surface): its area
    // predicts what it samples (a body of thin strips over a large box
    // samples its strips', not its box's), and its boundary points judge
    // the mesh made (`drawn`). Where none is found, its box's area.
    use rayon::prelude::*;
    let sized = how.seen && res.is_none();
    let surfaces: Vec<ink::seen::Surface> = match sized {
        true => {
            let fields: Vec<Arc<Tape>> = jobs.iter().map(|(it, _)| Arc::new(Tape::new(&s.doc.g, it.shape))).collect();
            fields.par_iter().zip(&boxes).map(|(t, b)| ink::seen::surface(t, &s.ctx, b.0, b.1)).collect()
        }
        false => vec![Default::default(); jobs.len()],
    };
    // Each mesh's triangles sampled at res 1 (all the times it is made), as
    // its area predicts them, and as its box does.
    let at_1: Vec<f64> = surfaces.iter().zip(&boxes).zip(&times).map(|((f, b), t)| t * SAMPLED_PER_CELL * if f.area > 0.0 { f.area } else { box_area(b.0, b.1) }).collect();
    let box_1: Vec<f64> = boxes.iter().zip(&times).map(|(b, t)| t * SAMPLED_PER_CELL * box_area(b.0, b.1)).collect();
    // For seeing, no mesh finer than the model's floor times the cube root
    // of the times it is made: its meshes' triangles sampled, as their
    // areas predict them, within VIEW_TRIANGLES.
    let parts: Vec<(f64, f64, f64)> = (0..jobs.len()).map(|j| (at_1[j], boxes[j].2, times[j])).collect();
    let floor = if sized { view_floor(&parts) } else { 0.0 };
    if sized {
        for (b, p) in boxes.iter_mut().zip(&parts) {
            b.2 = ladder(p.1.max(floor * p.2.cbrt())).min(coarsest);
        }
    }
    // Its deviation is sampled within 2 res of the mesh (mesh_tape).
    let plan: Vec<(String, usize, Option<Rigid>, f64)> =
        targets.into_iter().map(|(label, j, place, local)| (label, j, place, local.map_or(0.0, |(l, lo, hi)| l.off(lo, hi, 2.0 * boxes[j].2)))).collect();
    // Every mesh at once, each on its own thread (simplifying one is
    // serial, and most of the time).
    let tapes: Vec<Arc<Tape>> = jobs.iter().zip(&boxes).map(|((it, _), &(lo, hi, h, _))| tape_of(s, it, lo, hi, h)).collect();
    let paints: Vec<Vec<Arc<Tape>>> = jobs.iter().map(|(it, _)| paint_tapes(s, &it.attrs)).collect();
    let ctx = &s.ctx;
    let labels: Vec<&str> = jobs.iter().map(|(it, _)| it.label.as_str()).collect();
    let seen_how = |j: usize, h: f64| ink::seen::How { res: h, tol, tol_per_res: how.tol_per_res, bend: boxes[j].3 };
    // Each mesh (ink::seen::mesh, as `drawn` makes it: a view's made again
    // where it leaves its surface out, an export's once), what making it
    // cost (`work`: every mesh `drawn` made of it, the same whether made or
    // found in the cache) and what the one kept cost alone, with what
    // `how` measures of it.
    let mesh = |j: usize, h: f64| {
        let (lo, hi, _, _) = boxes[j];
        let finest = finest_drawn(how.seen && sized && tol.is_none() && h < coarsest, lo, hi, h);
        let (v, spent) = drawn(&tapes[j], &paints[j], ctx, lo, hi, seen_how(j, h), finest, &|| surfaces[j].clone(), labels[j]);
        let own = work(&v);
        let mut m = viewed(v);
        if how.measure || how.errors {
            measure(&m.0, &mut m.2, &tapes[j], ctx, lo, how);
        }
        (m, spent, own)
    };
    // For seeing, each mesh is reported as it is made: it is cached then.
    let mut total = REPORTED.load(std::sync::atomic::Ordering::Relaxed) + jobs.len();
    if how.seen {
        ink::claims::progress(json!({"total": total, "of": "meshes", "kept": true}));
    }
    let first: Vec<_> = (0..jobs.len()).into_par_iter().map(|j| mesh(j, boxes[j].2)).collect();
    // What the meshes cost (`work`), made or found in the cache alike
    // (`how.budget`), and each kept's own. A client's frames of a mesh
    // (Plan::frames) are made after: they are predicted, not counted.
    let mut spent: f64 = first.iter().map(|m| m.1).sum();
    let mut own: Vec<f64> = first.iter().map(|m| m.2).collect();
    let mut made: Vec<(Mesh, Vec<[f32; 3]>, Quality, Vec<V3>)> = first.into_iter().map(|m| m.0).collect();
    // A mesh that leaves material out (mesh::Unresolved: holes in a part, or
    // the part gone) is made again where that would draw more (`drew_more`),
    // at the res that resolves it: its res halved until finer than the
    // thickest of it, but no finer than `how.finest` (and not one at
    // --max-res: at the res asked for, as with --res). The cheapest first,
    // while what the meshes cost, a client's frames as predicted and what
    // the next is predicted to cost stay within `how.budget`. The
    // finer is kept where it drew more, and only then made again: a thin
    // long part (a belt) is resolved at once, a speck or a crest thinning to
    // nothing is left as it is. A gap thinner than res, two surfaces closed
    // into one (a ball resting in its race), is not made again: it looks as
    // it would apart, and where it thins to nothing no res resolves it.
    if how.seen && res.is_none() && tol.is_none() {
        let mut area: Vec<f64> = made.par_iter().map(|m| m.0.area()).collect();
        // The area each drew before it was last made again.
        let mut before: Vec<Option<f64>> = vec![None; jobs.len()];
        loop {
            // A mesh a client makes again per frame (Plan::frames), after
            // this: each frame predicted at its mesh's.
            let live: Vec<usize> = (0..jobs.len()).filter(|&j| times[j] > 1.0).collect();
            let (sampled, kept) = live.iter().fold((0, 0), |(a, b), &j| (a + made[j].2.uniform_triangles, b + made[j].0.tris.len()));
            let frames = match &frames {
                Some(f) if !live.is_empty() => f(sampled, kept) as f64 * live.iter().map(|&j| own[j]).sum::<f64>(),
                _ => 0.0,
            };
            let mut room = how.budget - spent - frames;
            // Each to make again, its res, and what it costs: four times what
            // it cost (at least a sheet in each cell it left out) per halving,
            // and at least its box's surface's triangles (a part left out whole
            // drew none).
            let mut again: Vec<(usize, f64, f64)> = (0..jobs.len())
                .filter(|&j| {
                    // What it leaves out, by area: a sheet as thick as the
                    // thickest of it, which only that share of the cells it
                    // passes through have their centre in.
                    let u = &made[j].2.unresolved;
                    let missing = 2.0 * u.material as f64 * made[j].2.res.powi(3) / u.thickest;
                    u.material > 0 && made[j].2.res < coarsest && before[j].is_none_or(|a| drew_more(a, area[j])) && drew_more(area[j], area[j] + missing)
                })
                .filter_map(|j| {
                    let q = &made[j].2;
                    let k = ((q.res / q.unresolved.thickest).log2().floor().max(0.0) + 1.0).min((q.res / how.finest).log2().floor()).min(60.0);
                    let h = q.res / 2f64.powf(k);
                    let cost = own[j].max(8.0 * q.unresolved.cells as f64);
                    (k >= 1.0).then(|| (j, h, (times[j] * 4f64.powf(k) * cost).max(box_1[j] / (h * h))))
                })
                .collect();
            again.sort_by(|a, b| a.2.total_cmp(&b.2));
            again.retain(|&(_, _, more)| {
                more <= room && {
                    room -= more;
                    true
                }
            });
            if again.is_empty() {
                break;
            }
            total += again.len();
            if how.seen {
                ink::claims::progress(json!({"total": total, "of": "meshes", "kept": true}));
            }
            let finer: Vec<_> = again.par_iter().map(|&(j, h, _)| (j, mesh(j, h))).collect();
            spent += finer.iter().map(|(_, m)| m.1).sum::<f64>();
            for (j, (m, _, c)) in finer {
                let a = m.0.area();
                before[j] = Some(area[j]);
                if drew_more(area[j], a) {
                    (made[j], area[j], own[j]) = (m, a, c);
                }
            }
        }
    }
    // For seeing, what each mesh was planned at and made at (a client's
    // frames counted as its first until it makes them).
    let budget = sized.then(|| Budget {
        jobs: plan.iter().map(|p| (p.0.clone(), p.1)).collect(),
        names: (0..jobs.len()).map(|j| plan.iter().find(|p| p.1 == j).map_or(jobs[j].0.label.clone(), |p| p.0.clone())).collect(),
        made: made.iter().zip(&times).map(|(m, t)| (m.2.res, t * m.2.uniform_triangles as f64)).collect(),
        parts,
        floor,
        boxes: boxes.clone(),
    });
    let made: Vec<(Arc<Mesh>, Arc<Vec<[f32; 3]>>, Arc<Vec<[f32; 3]>>, Arc<Vec<V3>>, Arc<Quality>)> = made
        .into_iter()
        .zip(&jobs)
        .map(|((m, normals, q, at), (it, k))| {
            let colors = s.colors(it, &at, *k);
            (Arc::new(m), Arc::new(colors), Arc::new(normals), Arc::new(at), Arc::new(q))
        })
        .collect();
    // A target whose body gives no colour takes the palette's by its place
    // among the targets: two sharing a mesh (two instances of a screw) are
    // told apart as two meshes are.
    let palette = |k: usize| ink::lang::look::PALETTE[k % ink::lang::look::PALETTE.len()];
    let unpainted: Vec<bool> = jobs.iter().map(|(it, _)| ink::lang::look::Look::of(&it.attrs).color.is_none()).collect();
    Ok((
        plan.into_iter()
            .enumerate()
            .map(|(k, (label, j, place, off))| {
                let (mesh, mut colors, normals, paint_at, quality) = made[j].clone();
                if unpainted[j] {
                    colors = Arc::new(vec![palette(k); colors.len()]);
                }
                Meshed { label, mesh, colors, normals, paint_at, quality, place, off }
            })
            .collect(),
        budget,
    ))
}

/// Meshes for seeing as planned and made (make_meshes), for
/// `cost_warnings`: each target's mesh, each mesh's first target, its res
/// and the triangles it sampled (all its meshes: a client's frames counted
/// as its first until it makes them and says what they sampled),
/// `view_floor`'s parts and floor, and each mesh's box and planned res.
pub struct Budget {
    pub jobs: HashMap<String, usize>,
    pub names: Vec<String>,
    pub made: Vec<(f64, f64)>,
    pub parts: Vec<(f64, f64, f64)>,
    pub floor: f64,
    pub boxes: Vec<(V3, V3, f64, Option<f64>)>,
}

impl Budget {
    /// What the meshes cost (`cost_warnings`).
    pub fn warnings(&self) -> Vec<String> {
        cost_warnings(&self.parts, self.floor, &self.boxes, &self.made, &self.names.iter().map(String::as_str).collect::<Vec<_>>())
    }
}

/// The finest a view mesh is made at where it leaves its surface out
/// (`drawn`): where `free` (its res is the plan's: no --res or --accuracy,
/// and not at --max-res), what the cells' limit allows across its box (as
/// `--res` does), else its res.
pub fn finest_drawn(free: bool, lo: V3, hi: V3, h: f64) -> f64 {
    if free { max_extent(lo, hi) / 20000.0 } else { h }
}

/// A view mesh that draws its body: made at `how.res`, and made again,
/// whatever the budget, while its surface's boundary points (`surface`,
/// asked for only then; ink::seen::Surface::off) lie farther from it on
/// average than its res, the error the floor plans each mesh at
/// (`view_floor`): at the ladder's res under the thickest material it
/// leaves out (half the res where it measures none), no finer than
/// `finest`. The best is kept, and none is made after one that does not
/// bring the mean down: a part the floor brought past its features (a belt
/// run, a pulley) is drawn, a detail within a cell of the mesh is not
/// chased, nor a face of no thickness. Each mesh is reported as made or
/// found (`label`). Also what every mesh made of it cost (`work`), as much
/// whether it was made or found in the cache.
pub fn drawn(
    tape: &Arc<Tape>,
    paint: &[Arc<Tape>],
    ctx: &ink::eval::Ctx,
    lo: V3,
    hi: V3,
    mut how: ink::seen::How,
    finest: f64,
    surface: &dyn Fn() -> ink::seen::Surface,
    label: &str,
) -> (ink::seen::Viewed, f64) {
    let (mut spent, mut asked) = (0.0, None);
    let mut best: Option<(ink::seen::Viewed, f64)> = None;
    loop {
        let t = Instant::now();
        let v = ink::seen::mesh(tape, paint, ctx, lo, hi, how);
        reported(&format!("{label} at res {}", measured(how.res)), v.cached, t);
        spent += work(&v);
        let h = how.res;
        let off = (h > finest).then(|| asked.get_or_insert_with(surface).off(&v.mesh(), h, tape, ctx));
        let mean = off.map_or(0.0, |o| o.mean);
        if best.as_ref().is_some_and(|b| mean >= b.1) {
            break;
        }
        best = Some((v, mean));
        match off {
            Some(o) if o.mean > h => how.res = if o.thickest > 0.0 && o.thickest < h { ladder(o.thickest / 2f64.powf(0.25)) } else { h / 2.0 }.max(finest),
            _ => break,
        }
    }
    (best.unwrap().0, spent)
}

/// What making a view mesh costs, in triangles sampled: its own, and the
/// instructions reading its paint ran at PAINT_PER_TRIANGLE a triangle.
fn work(v: &ink::seen::Viewed) -> f64 {
    v.uniform_triangles as f64 + v.paint_work as f64 / PAINT_PER_TRIANGLE
}

/// A viewing mesh as `mesh_tape` gives one.
fn viewed(v: ink::seen::Viewed) -> (Mesh, Vec<[f32; 3]>, Quality, Vec<V3>) {
    let q = Quality { res: v.res, tol: v.tol, uniform_triangles: v.uniform_triangles, moved: v.moved, deviation: None, unresolved: v.unresolved, error: Vec::new() };
    (Mesh { verts: v.verts, tris: v.tris }, v.normals, q, v.paint_at)
}

/// The paint values of an item's look (its colour's channels, metallic,
/// roughness, transmission) that vary over space.
pub fn paint_tapes(s: &Session, attrs: &ink::lang::Attrs) -> Vec<Arc<Tape>> {
    let fs = Look::of(attrs).channels().iter().flatten().copied().collect::<Vec<_>>();
    let spatial = |f: u32| s.doc.g.depends_on(f, |v| matches!(v, ink::graph::Var::X | ink::graph::Var::Y | ink::graph::Var::Z));
    fs.into_iter().filter(|&f| spatial(f)).map(|f| Arc::new(Tape::new(&s.doc.g, f))).collect()
}

/// Metallic, roughness and transmission at each vertex (0..1; the look's
/// defaults where it gives none).
pub fn looks(s: &Session, attrs: &ink::lang::Attrs, verts: &[V3]) -> Vec<[f32; 3]> {
    use ink::lang::look::{METALLIC, ROUGHNESS, TRANSMISSION};
    let look = Look::of(attrs);
    paint(s, [look.metallic, look.roughness, look.transmission], [METALLIC, ROUGHNESS, TRANSMISSION], verts)
}

/// An item's field clipped to its box and a cell beyond.
pub fn tape_of(s: &mut Session, it: &Item, lo: V3, hi: V3, h: f64) -> Arc<Tape> {
    let f = s.doc.clip(it.shape, [lo[0] - h, lo[1] - h, lo[2] - h], [hi[0] + h, hi[1] + h, hi[2] + h]);
    Arc::new(Tape::new(&s.doc.g, f))
}

/// What `how` asks of a mesh: for `ink mesh` (`measure`), its sampled
/// distance from the surface (ink::simplify::deviation, sampled within 2
/// res of it, every quarter of its tolerance); for `ink render`
/// (`errors`), how far each vertex may lie from the surface: the field's
/// value over its gradient's length there and at the centroids of its
/// triangles (simplified, a triangle strays between its corners).
fn measure(m: &Mesh, q: &mut Quality, tape: &Arc<Tape>, ctx: &ink::eval::Ctx, lo: V3, how: &MeshHow) {
    let h = q.res;
    // The field near the surface, a tape per block of 8 cells.
    let local = ink::octree::Local::new(tape, ctx, lo, 8.0 * h);
    if how.measure {
        let step = if q.tol > 0.0 { q.tol } else { h / 30.0 } / 4.0;
        q.deviation = Some(ink::simplify::deviation(&|p: &[V3]| local.eval(p), &|p: &[V3]| local.grads(p), m, h, step, 2.0 * h));
    }
    if how.errors {
        let error = |g: [f64; 4]| Some((g[0].abs() / ink::v3::norm([g[1], g[2], g[3]])) as f32).filter(|e| e.is_finite()).unwrap_or(0.0);
        let mut e: Vec<f32> = local.grads(&m.verts).into_iter().map(error).collect();
        let centroids: Vec<V3> = m.tris.iter().map(|t| std::array::from_fn(|k| t.iter().map(|&i| m.verts[i as usize][k]).sum::<f64>() / 3.0)).collect();
        for (t, g) in m.tris.iter().zip(local.grads(&centroids)) {
            t.iter().for_each(|&i| e[i as usize] = e[i as usize].max(error(g)));
        }
        q.error = e;
    }
}

/// Warnings, one per distinct mesh: a target meshed to nothing (empty, or
/// thinner than res everywhere), and sign changes a mesh left out.
fn mesh_warnings(meshed: &[Meshed]) -> Vec<String> {
    let mut seen = std::collections::HashSet::new();
    meshed
        .iter()
        .filter(|m| seen.insert(Arc::as_ptr(&m.mesh) as usize))
        .filter_map(|m| {
            if m.mesh.tris.is_empty() {
                return Some(format!("'{}': its mesh at res {} has no triangles. It is empty, or thinner than res everywhere. {EMPTY_CAUSE}", m.label, measured(m.quality.res)));
            }
            let (q, at) = (&m.quality, m.quality.unresolved.at?);
            let at = m.place.map_or(at, |p| p.apply(at));
            Some(format!(
                "'{}': {} cells hold material or a gap thinner than res {} that the mesh leaves out (e.g. the cell centred at [{}, {}, {}]): \
                 it has holes or misses features there; use a smaller --res",
                m.label,
                q.unresolved.cells,
                measured(q.res),
                measured(at[0]),
                measured(at[1]),
                measured(at[2])
            ))
        })
        .collect()
}

/// Warnings, one per distinct mesh measured, where an accuracy was asked
/// for (`--accuracy`, not 0) and the mesh lies farther from the surface
/// (its sampled deviation): simplifying keeps within it of the mesh at
/// res, not of the surface.
fn missed_accuracy(meshed: &[Meshed], asked: Option<f64>) -> Vec<String> {
    let Some(asked) = asked.filter(|t| *t > 0.0) else { return Vec::new() };
    let mut seen = std::collections::HashSet::new();
    meshed
        .iter()
        .filter(|m| seen.insert(Arc::as_ptr(&m.mesh) as usize))
        .filter_map(|m| {
            let d = m.quality.deviation?;
            let at = m.place.map_or(d.at, |p| p.apply(d.at));
            (d.max + m.off > asked).then(|| {
                format!(
                    "'{}': the mesh lies up to {} from the surface (sampled, at [{}, {}, {}]), more than --accuracy {}, which bounds only \
                     how far simplifying moves the mesh at res {}: that mesh is off by this much itself, where edges, steps or walls lie \
                     closer together than res or the surface curves within a cell; a smaller --res brings it closer",
                    m.label,
                    measured(d.max + m.off),
                    measured(at[0]),
                    measured(at[1]),
                    measured(at[2]),
                    measured(asked),
                    measured(m.quality.res)
                )
            })
        })
        .collect()
}

fn unresolved_json(m: &Meshed) -> J {
    let u = &m.quality.unresolved;
    json!({"cells": u.cells, "at": u.at.map(|p| v3(m.place.map_or(p, |r| r.apply(p))))})
}

fn mesh_report(m: &Meshed, cache: &mut HashMap<usize, Arc<mesh::Analysis>>) -> J {
    let key = Arc::as_ptr(&m.mesh) as usize;
    let an = cache.entry(key).or_insert_with(|| Arc::new(mesh::analyze(&m.mesh))).clone();
    let q = &m.quality;
    let mut r = json!({"target": m.label, "res": n(q.res), "vertices": m.mesh.verts.len(), "triangles": m.mesh.tris.len(),
           "uniform_triangles": q.uniform_triangles,
           "watertight": an.boundary_edges == 0 && an.nonmanifold_edges == 0,
           "boundary_edges": an.boundary_edges, "nonmanifold_edges": an.nonmanifold_edges,
           "shells": an.shells.len(), "mesh_volume": n(an.volume), "mesh_area": n(an.area)});
    if let Some(d) = q.deviation {
        let place = |p: V3| v3(m.place.map_or(p, |r| r.apply(p)));
        // In the world: the instance's surface may be `off` farther.
        r["deviation"] = json!({"max": n(d.max + m.off), "at": place(d.at), "simplified_within": n(q.tol), "simplification": n(q.moved),
            "beyond": d.beyond, "beyond_at": d.first_beyond.map(place), "reach": n(2.0 * q.res)});
    }
    r["unresolved"] = unresolved_json(m);
    if an.nonmanifold_edges > 0 {
        r["hint"] = json!(format!("Features about one cell ({}) wide pinch into shared edges. Mesh at a smaller --res", measured(q.res)));
    }
    r
}

pub(crate) fn cmd_mesh(s: &mut Session, a: &Args) -> Result<J, InkError> {
    let t0 = Instant::now();
    let out = a.get("out").map(str::to_string).unwrap_or(format!("out/{}.3mf", s.stem()));
    // GLB is for seeing: no deviation measured.
    let glb = out.ends_with(".glb");
    let how = if glb { SEEN } else { MeshHow { div: 256.0, tol_per_res: EXPORT_TOL_PER_RES, measure: true, budget: 0.0, finest: 0.0, errors: false, seen: false } };
    let items = s.items(a)?;
    let (meshed, budget) = mesh_items(s, a, &items, how)?;
    let mesh_ms = t0.elapsed().as_secs_f64() * 1e3;
    // The file says how it was made: res and the sampled deviation (mm, 4
    // significant digits).
    let short = |v: f64| format!("{v:.3e}").parse::<f64>().map_or(String::new(), |v| v.to_string());
    let note = |m: &Meshed| {
        let q = &m.quality;
        let dev = q.deviation.map_or(String::new(), |d| format!(" dev {} (sampled)", short(d.max + m.off)));
        format!("{}: res {}{dev}", m.label, short(q.res))
    };
    let mut warnings = budget.map_or_else(Vec::new, |b| b.warnings());
    let bytes = if out.ends_with(".stl") {
        // STL has no instancing: every copy is written out in full.
        let total: usize = meshed.iter().map(|m| m.mesh.tris.len()).sum();
        if total > 50_000_000 {
            return Err(InkError::msg(format!("{total} triangles is too many for STL (no instancing); write .3mf instead")));
        }
        let worlds: Vec<Mesh> = meshed.iter().map(|m| m.world()).collect();
        let objs: Vec<export::Object> = worlds.iter().zip(&meshed).map(|(w, m)| export::Object { name: &m.label, mesh: w, colors: &m.colors }).collect();
        let header = if meshed.len() == 1 { note(&meshed[0]) } else { format!("{} objects", meshed.len()) };
        export::stl(&objs, &header)
    } else if out.ends_with(".3mf") {
        // One mesh resource per distinct mesh; one build item per instance.
        let mut index: HashMap<usize, usize> = HashMap::new();
        let mut uniques: Vec<&Meshed> = Vec::new();
        let mut placed = Vec::new();
        for m in &meshed {
            let key = Arc::as_ptr(&m.mesh) as usize;
            let k = *index.entry(key).or_insert_with(|| {
                uniques.push(m);
                uniques.len() - 1
            });
            placed.push(export::Placed { object: k, name: &m.label, transform: m.place.map(|p| (p.r, p.t)) });
        }
        let objs: Vec<export::Object> = uniques.iter().map(|m| export::Object { name: &m.label, mesh: &m.mesh, colors: &m.colors }).collect();
        let description = format!("mm; {}", uniques.iter().map(|m| note(m)).collect::<Vec<_>>().join("; "));
        export::three_mf_instanced(&objs, &placed, &description)
    } else if glb {
        let (bytes, shapes) = glb_bytes(s, a, &meshed)?;
        warnings.extend(shapes);
        bytes
    } else {
        return Err(InkError::msg("--out must end in .stl, .3mf or .glb"));
    };
    export::write(&out, &bytes).map_err(|e| InkError::msg(format!("cannot write {out}: {e}")))?;
    mesh_answer(s, a, &out, bytes.len(), mesh_ms, &meshed, warnings)
}

/// `ink mesh`'s answer, a file of `bytes` written to `out` from `meshed`:
/// each target's mesh (`mesh_report`), and the warnings: `first`, then each
/// mesh's own, an accuracy missed, and a target meshed to nothing.
pub fn mesh_answer(s: &mut Session, a: &Args, out: &str, bytes: usize, mesh_ms: f64, meshed: &[Meshed], first: Vec<String>) -> Result<J, InkError> {
    let mut cache = HashMap::new();
    let distinct = meshed.iter().map(|m| Arc::as_ptr(&m.mesh) as usize).collect::<std::collections::HashSet<_>>().len();
    let objects: Vec<J> = meshed.iter().map(|m| mesh_report(m, &mut cache)).collect();
    let mut r = json!({"out": out, "bytes": bytes, "mesh_ms": n(mesh_ms), "instances": meshed.len(), "distinct_meshes": distinct, "objects": objects});
    if out.ends_with(".glb") {
        // Viewing meshes taken from the cache and made.
        use std::sync::atomic::Ordering::Relaxed;
        r["meshes_cached"] = json!(ink::seen::HITS.load(Relaxed));
        r["meshes_made"] = json!(ink::seen::MADE.load(Relaxed));
    }
    let mut warnings = first;
    warnings.extend(mesh_warnings(meshed));
    warnings.extend(missed_accuracy(meshed, a.f64("accuracy")?));
    // A target with no material meshes to nothing: say so, not just 0.
    warnings.extend(
        meshed
            .iter()
            .filter(|m| m.mesh.tris.is_empty())
            .map(|m| format!("{}: no surface at res {}, so the file has none of it. It is empty, or thinner than res", m.label, m.quality.res)),
    );
    let emitted: std::collections::HashSet<&str> = meshed.iter().map(|m| m.label.as_str()).collect();
    warnings.extend(s.items(a)?.iter().filter(|it| !emitted.contains(it.label.as_str())).map(|it| format!("{}: empty target; the file has no triangles of it", it.label)));
    if !warnings.is_empty() {
        r["warnings"] = json!(warnings);
    }
    Ok(r)
}

/// Meshes reported done in this process (`reported`), which a later total
/// counts on from (a GLB's flipbook frames after its meshes).
pub static REPORTED: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// Report a mesh done (ink::claims::progress): made, or found in the
/// cache, since `t`.
pub fn reported(what: &str, cached: bool, t: Instant) {
    REPORTED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let status = if cached { "cached" } else { "made" };
    ink::claims::progress(json!({"done": what, "status": status, "ms": n(t.elapsed().as_secs_f64() * 1e3)}));
}

/// GLB's meshes, for seeing at a screen's resolution: coarser than for
/// printing, and simplified to a quarter of res. Their budget is about 65
/// CPU s on 6 vCPU of a Ryzen AI Max+ 395 at corexy's 14 µs a triangle: its first
/// meshes cost 2.5 M and the 7 it makes again 1.9 M more as predicted;
/// A large assembly's first meshes 8.3 M, so none is made again.
pub const SEEN: MeshHow = MeshHow { div: 128.0, tol_per_res: 0.25, measure: false, budget: 4.6e6, finest: 0.0, errors: false, seen: true };

/// The triangles a model's GLB meshes are planned to sample when first
/// made, as their surfaces' areas predict them (`view_floor`): the rest of
/// SEEN's budget is for meshes made again.
const VIEW_TRIANGLES: f64 = 2.5e6;

/// The paint instructions that cost what sampling and simplifying a
/// triangle does (`work`): on 6 vCPU of a Ryzen AI Max+ 395 a paint instruction
/// costs 3.0 CPU ns (2.5 to 5.6 over six boards, each with and without its
/// paint), corexy 13.4 CPU µs a triangle besides its paint.
const PAINT_PER_TRIANGLE: f64 = 4.4e3;

/// Whether a mesh made again draws more than the one before it: a
/// thousandth more area (a speck on the part, or a crest drawn a little
/// further toward nothing, is not worth four times the triangles).
fn drew_more(before: f64, after: f64) -> bool {
    after > 1.001 * before
}

/// The triangles a surface samples per cell of it, as predicted: a flat
/// face's two (a folded surface samples more).
const SAMPLED_PER_CELL: f64 = 2.0;

fn box_area(lo: V3, hi: V3) -> f64 {
    let d = ink::v3::sub(hi, lo);
    2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0])
}

/// The model's floor for its view meshes: 0 where each at its own res
/// samples at most VIEW_TRIANGLES, else the least `f` at which they do,
/// each mesh at f times the cube root of the times it is made (a
/// flipbook's frames) where that is coarser than its own res. That is
/// the least error over what is seen for the triangles: res h over area A
/// is off by about h over A, each mesh is seen once at a time, and it
/// samples 2 A / h^2 each time it is made, so the sum of A h within the
/// budget is least where every h^3 is in proportion to the times made.
/// So small detailed parts are coarsened first, large ones are not made
/// finer, and a mesh made 40 times counts as about 3.4 of its size, not
/// 40. `parts`: each mesh's triangles predicted at res 1 (all the times it
/// is made), its own res, and the times it is made.
fn view_floor(parts: &[(f64, f64, f64)]) -> f64 {
    let sampled = |f: f64| parts.iter().map(|&(c, h, t)| c / h.max(f * t.cbrt()).powi(2)).sum::<f64>();
    if sampled(0.0) <= VIEW_TRIANGLES {
        return 0.0;
    }
    // At hi every mesh's res is hi or coarser: within VIEW_TRIANGLES.
    let mut lo = parts.iter().map(|p| p.1 / p.2.cbrt()).fold(f64::INFINITY, f64::min);
    let mut hi = (parts.iter().map(|p| p.0).sum::<f64>() / VIEW_TRIANGLES).sqrt();
    for _ in 0..60 {
        let mid = (lo * hi).sqrt();
        if sampled(mid) <= VIEW_TRIANGLES { hi = mid } else { lo = mid }
    }
    hi
}

/// Warnings for each mesh that sampled more than half of the triangles
/// the meshes sampled, where they sampled more than VIEW_TRIANGLES, or
/// whose cost alone brought another mesh's res up a step of the ladder
/// (only a mesh predicted at over 1 - 2^(-1/2) of them can: without it the
/// floor falls at most as the square root of what is left): it, its
/// frames, res and size, and what it sampled. `parts` as `view_floor`
/// takes them, `boxes` each mesh's box and planned res, `made` each mesh's
/// res and the triangles it sampled (its frames' as many as its first's),
/// `names` each mesh's first target. And, where they sampled more than
/// VIEW_TRIANGLES, how many were made finer than planned.
fn cost_warnings(parts: &[(f64, f64, f64)], floor: f64, boxes: &[(V3, V3, f64, Option<f64>)], made: &[(f64, f64)], names: &[&str]) -> Vec<String> {
    let planned: Vec<f64> = parts.iter().zip(boxes).map(|(p, b)| p.0 / (b.2 * b.2)).collect();
    let (planned_all, all) = (planned.iter().sum::<f64>(), made.iter().map(|m| m.1).sum::<f64>());
    let mega = |t: f64| format!("{:.2} M", t / 1e6);
    // A mesh's res at floor f.
    let at = |k: usize, f: f64| ladder(parts[k].1.max(f * parts[k].2.cbrt()));
    (0..parts.len())
        .filter_map(|j| {
            let most = all > VIEW_TRIANGLES && made[j].1 > 0.5 * all;
            let without = (floor > 0.0 && planned[j] > (1.0 - 0.5f64.sqrt()) * planned_all)
                .then(|| view_floor(&parts.iter().enumerate().filter(|&(k, _)| k != j).map(|(_, p)| *p).collect::<Vec<_>>()))
                .filter(|&w| (0..parts.len()).any(|k| k != j && at(k, w) < at(k, floor)));
            (most || without.is_some()).then(|| {
                let what = match parts[j].2 {
                    1.0 => "its mesh".to_string(),
                    t => format!("its {} frames", t - 1.0),
                };
                let mut w = format!(
                    "'{}': {what} at res {} over {} mm sampled {} of the {} triangles the GLB's meshes sampled ({:.0}%)",
                    names[j],
                    measured(made[j].0),
                    measured(max_extent(boxes[j].0, boxes[j].1)),
                    mega(made[j].1),
                    mega(all),
                    100.0 * made[j].1 / all
                );
                if let Some(w0) = without {
                    let others = if w0 == 0.0 { "each at its own res".to_string() } else { format!("res {}", measured(ladder(w0))) };
                    w += &format!("; it brings the other meshes up to res {} (without it, {others})", measured(ladder(floor)));
                }
                w
            })
        })
        .chain(
            // Past the budget: what drawing them took.
            {
                let finer: Vec<f64> = made.iter().zip(boxes).filter(|(m, b)| m.0 < b.2).map(|(m, _)| m.0).collect();
                (all > VIEW_TRIANGLES && !finer.is_empty()).then(|| {
                    format!(
                        "the GLB's meshes sampled {} triangles, {:.1} times the {} they are planned within: made finer than planned (to res {} \
                         at the finest), {} of its {} meshes drew what the floor (res {}) brought past their features",
                        mega(all),
                        all / VIEW_TRIANGLES,
                        mega(VIEW_TRIANGLES),
                        measured(finer.iter().fold(f64::INFINITY, |a, &r| a.min(r))),
                        finer.len(),
                        made.len(),
                        measured(ladder(floor))
                    )
                })
            },
        )
        .collect()
}

/// `h` rounded up to the ladder of view resolutions, 2^(k/4) mm: a part
/// gets the same res, so the same cached mesh, in models whose floors
/// differ a little.
fn ladder(h: f64) -> f64 {
    const STEP: [f64; 4] = [1.0, 1.189207115002721, std::f64::consts::SQRT_2, 1.681792830507429];
    let k = (4.0 * h.log2() - 1e-9).ceil() as i32;
    2f64.powi(k.div_euclid(4)) * STEP[k.rem_euclid(4) as usize]
}

/// Motion samples per second of the span, with their limits.
const MOTION_RATE: f64 = 30.0;
const MOTION_SAMPLES: (usize, usize) = (60, 1200);

/// Moments over the span: `rate` a second, between lo and hi of them but
/// at most n (and at least 1), from its start to its end.
pub fn span_times(span: (f64, f64), n: usize, (lo, hi): (usize, usize), rate: f64) -> Vec<f64> {
    let len = span.1 - span.0;
    let n = ((rate * len).ceil() as usize).clamp(lo, hi).min(n).max(1);
    (0..=n).map(|i| span.0 + len * i as f64 / n as f64).collect()
}

/// The moments a GLB samples motion at, and the context at each.
pub fn motion(s: &Session) -> Result<(Vec<f64>, Vec<ink::eval::Ctx>), InkError> {
    let times = span_times(s.doc.t_span.unwrap_or((0.0, 1.0)), usize::MAX, MOTION_SAMPLES, MOTION_RATE);
    let ctxs = times.iter().map(|&t| s.doc.ctx(&s.overrides, t).map_err(|e| InkError::msg(format!("at t = {t}: {e}")))).collect::<Result<_, _>>()?;
    Ok((times, ctxs))
}

/// The meshes as a GLB with the motion over the file's span: each instance
/// whose placement reads t has its placement sampled (rigid motion). A body
/// whose shape itself reads t is meshed at the moment asked (--t), and said
/// so (warnings).
fn glb_bytes(s: &mut Session, a: &Args, meshed: &[Meshed]) -> Result<(Vec<u8>, Vec<String>), InkError> {
    use ink::gltf::{GlbMesh, GlbNode, GlbScene, Shape};
    let items = s.items(a)?;
    let by_label: HashMap<&str, &Item> = items.iter().map(|it| (it.label.as_str(), it)).collect();
    let reads_t = |s: &Session, roots: &[u32]| roots.iter().any(|&f| s.doc.g.depends_on(f, |v| v == Var::T)) || s.doc.solves_read(roots).iter().any(|&i| s.doc.solves[i].on_t);
    // Each meshed item, and whether its placement moves.
    let parts: Vec<(&Item, bool)> = meshed
        .iter()
        .map(|m| {
            let it = by_label[m.label.as_str()];
            let entries: Vec<u32> = it.instance.as_ref().filter(|_| m.place.is_some()).map_or(vec![], |(_, xf)| xf.r.iter().flatten().chain(&xf.t).copied().collect());
            (it, reads_t(s, &entries))
        })
        .collect();
    let span = s.doc.t_span.unwrap_or((0.0, 1.0));
    let any_moves = parts.iter().any(|p| p.1);
    let (times, ctxs) = if any_moves { motion(s)? } else { (vec![], vec![]) };
    // Each body's attributes (its catalog record, material, datasheet), once.
    let mut bodies = serde_json::Map::new();
    let mut owned: Vec<(String, Arc<Mesh>, Arc<Vec<[f32; 3]>>, Arc<Vec<[f32; 3]>>, Vec<[f32; 3]>)> = Vec::new();
    let mut index: HashMap<usize, usize> = HashMap::new();
    let mut warnings = Vec::new();
    let mut nodes = Vec::new();
    for (m, &(it, moves)) in meshed.iter().zip(&parts) {
        let key = Arc::as_ptr(&m.mesh) as usize;
        let (name, attrs) = match &it.instance {
            Some((body, _)) if m.place.is_some() => (body.clone(), it.attrs.clone()),
            _ => (m.label.clone(), s.doc.placed_attrs(it)),
        };
        if let std::collections::hash_map::Entry::Vacant(e) = index.entry(key) {
            e.insert(owned.len());
            let field = match &it.instance {
                Some(_) if m.place.is_some() => it.shape,
                _ => s.doc.placed_field(it),
            };
            if reads_t(s, &[field]) {
                warnings.push(format!("'{}': its shape changes over t; the GLB has it as it is at t = {} (its placement moves)", m.label, measured(s.ctx.t)));
            }
            owned.push((name.clone(), m.mesh.clone(), m.colors.clone(), m.normals.clone(), looks(s, &attrs, &m.paint_at)));
            bodies.insert(name.clone(), glb_attrs(s, &attrs));
        }
        let motion = match (&it.instance, moves) {
            (Some((_, xf)), true) => {
                let tapes: Vec<Tape> = xf.r.iter().flatten().chain(&xf.t).map(|&e| Tape::new(&s.doc.g, e)).collect();
                Some(
                    ctxs.iter()
                        .map(|c| {
                            let v: Vec<f64> = tapes.iter().map(|tp| ink::eval::eval_point(tp, c, [0.0; 3])).collect();
                            ([[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]], [v[9], v[10], v[11]])
                        })
                        .collect(),
                )
            }
            _ => None,
        };
        // A body built with arguments says which it is an instance of (a
        // parts list counts a catalog part's variants together).
        let extras = match s.doc.body(&name).and_then(|b| b.of.clone()) {
            Some(of) => json!({"body": name, "of": of}),
            None => json!({"body": name}),
        };
        nodes.push(GlbNode { name: &m.label, shape: Shape::Mesh(index[&key]), place: m.place.map(|r| (r.r, r.t)), motion, extras: extras.to_string() });
    }
    let meshes: Vec<GlbMesh> = owned.iter().map(|(name, mesh, colors, normals, looks)| GlbMesh { name, mesh, colors, normals, looks }).collect();
    let mut info = json!({"file": s.stem(), "units": "mm", "up": "z", "bodies": bodies});
    if any_moves {
        info["span"] = bounds(span.0, span.1);
    }
    let times = times.iter().map(|t| t - span.0).collect();
    let scene = GlbScene { meshes, nodes, times, flip_times: Vec::new(), extras: json!({"ink": info}).to_string() };
    Ok((ink::gltf::glb(&scene), warnings))
}

/// A body's attributes as a GLB's scene extras give them: strings,
/// constant numbers and data.
pub fn glb_attrs(s: &Session, attrs: &ink::lang::Attrs) -> J {
    let mut out = serde_json::Map::new();
    for (k, v) in attrs.iter() {
        let j = match v {
            Value::Str(x) => json!(x),
            Value::Field(f, _) => match s.doc.g.as_const(*f) {
                Some(c) => n(c),
                None => continue,
            },
            Value::Data(d) => d.value().clone(),
            _ => continue,
        };
        out.insert(k.clone(), j);
    }
    J::Object(out)
}

pub(crate) fn cmd_components(s: &mut Session, a: &Args) -> Result<J, InkError> {
    // Each target's shells, from a mesh (an estimate), and its cavities
    // counted as the `cavities` claim counts them (ModelView::components).
    let res = a.cell("res")?;
    let mut out = Vec::new();
    let mut warnings = Vec::new();
    for it in s.items(a)? {
        let field = s.doc.placed_field(&it);
        let Some(m) = s.components(field, res).map_err(|e| InkError::msg(format!("'{}': {}", it.label, e.message)))? else {
            warnings.push(format!("{}: empty target", it.label));
            continue;
        };
        let an = &m.analysis;
        let shells: Vec<J> = an
            .shells
            .iter()
            .map(|sh| {
                json!({"kind": if sh.volume >= 0.0 {"solid"} else {"cavity"}, "volume": n(sh.volume.abs()),
                             "area": n(sh.area), "min": v3(sh.lo), "max": v3(sh.hi)})
            })
            .collect();
        let solids = an.shells.iter().filter(|s| s.volume >= 0.0).count();
        if m.unresolved.cells > 0 {
            warnings.push(format!(
                "{}: {} cells hold material or a gap thinner than res {}, which the mesh misses. A smaller --res resolves them",
                it.label,
                m.unresolved.cells,
                measured(m.res)
            ));
        }
        let c = &m.voids;
        if let Some((what, decide)) = c.why_undecided() {
            let at = c.at.map_or(String::new(), |p| format!(" (one at [{:.3}, {:.3}, {:.3}])", p[0], p[1], p[2]));
            let most = c.at_most.map_or("unbounded".into(), |n| n.to_string());
            warnings.push(format!("{}: its cavities are proven between {} and {most}: {what}{at}. To decide it: {decide}", it.label, c.at_least));
        }
        out.push(json!({"target": it.label, "res": n(m.res), "solids": solids, "cavities": c.count(),
                        "proven": {"cavities": {"at_least": c.at_least, "at_most": c.at_most, "at": c.at.map(v3), "cells": c.cells, "finest": n(c.finest),
                                               "undecided": c.undecided.map(|(a, b)| json!([n(a), n(b)])), "ended": c.ended}},
                        "watertight": an.boundary_edges == 0 && an.nonmanifold_edges == 0, "shells": shells,
                        "unresolved": {"cells": m.unresolved.cells, "at": m.unresolved.at.map(v3)}}));
    }
    let mut r = json!({"targets": out});
    if !warnings.is_empty() {
        r["warnings"] = json!(warnings);
    }
    Ok(r)
}

fn parse_v3(a: &Args, key: &str) -> Result<V3, InkError> {
    let v = json_arg(a, key)?.ok_or_else(|| InkError::msg(format!("needs --{key} '[x,y,z]'")))?;
    json_v3(&v).ok_or_else(|| InkError::msg(format!("--{key} must be a JSON array of 3 numbers [x,y,z]")))
}

/// Construction names (outermost first) at a point, by provenance.
fn names_at(s: &Session, tape: &Tape, p: V3) -> Vec<String> {
    ink::eval::provenance(tape, &s.ctx, p).iter().map(|&t| s.doc.g.tag_names[t as usize].clone()).collect()
}

pub(crate) fn cmd_ray(s: &mut Session, a: &Args) -> Result<J, InkError> {
    use rayon::prelude::*;
    // One ray (--from, --dir) or many (--rays '[[[ox,oy,oz],[dx,dy,dz]], ...]',
    // or @file / - like --points): the targets compile once for all of them.
    let rays: Vec<(V3, V3)> = match json_arg(a, "rays")? {
        Some(v) => v
            .as_array()
            .and_then(|l| {
                l.iter()
                    .map(|r| {
                        let r = r.as_array().filter(|r| r.len() == 2)?;
                        Some((json_v3(&r[0])?, json_v3(&r[1])?))
                    })
                    .collect::<Option<Vec<_>>>()
            })
            .ok_or_else(|| InkError::msg("--rays expects [[[ox,oy,oz],[dx,dy,dz]], ...]"))?,
        None => vec![(parse_v3(a, "from")?, parse_v3(a, "dir")?)],
    };
    let rays: Vec<(V3, V3)> = rays
        .into_iter()
        .map(|(o, d)| Ok((o, ink::v3::direction(d).ok_or_else(|| InkError::msg("ray directions must be non-zero and finite"))?)))
        .collect::<Result<_, InkError>>()?;
    let items = s.items(a)?;
    let fields: Vec<u32> = items.iter().map(|it| s.doc.placed_field(it)).collect();
    let targets: Vec<(String, Arc<Tape>, Tape)> = items.iter().zip(&fields).map(|(it, &f)| (it.label.clone(), s.clipped(f), Tape::new(&s.doc.g, f))).collect();
    // Where the targets can be: their reaches (the window for those not
    // proven bounded, which a ray may hit beyond).
    let reaches: Vec<Reach> = fields.iter().map(|&f| s.reach(f)).collect();
    let (lo, hi) = reaches.iter().map(|r| (r.lo, r.hi)).fold(ink::v3::EMPTY_BOX, ink::v3::box_union);
    let open = reaches.iter().any(|r| !r.whole);
    let (ctx, tag_names) = (&s.ctx, &s.doc.g.tag_names);
    // Tapes simplified per block of where the targets can be, shared by
    // the rays through each block.
    let ext = max_extent(lo, hi);
    let locals: Vec<ink::octree::Local> = targets.iter().map(|(_, tape, _)| ink::octree::Local::new(tape, ctx, lo, ext / 256.0)).collect();
    let out: Vec<J> = rays
        .par_iter()
        .map(|&(o, dir)| {
            // Default reach: until the ray leaves where the targets can be.
            let exit = (0..3).filter(|&k| dir[k] != 0.0).map(|k| ((if dir[k] > 0.0 { hi[k] } else { lo[k] }) - o[k]) / dir[k]).fold(f64::INFINITY, f64::min).max(0.0);
            let tmax = exit;
            let mut hits = Vec::new();
            for ((label, _, raw), local) in targets.iter().zip(&locals) {
                // Where the ray enters material. One starting inside says
                // so, at its origin: it has entered nowhere, so no normal.
                if let Some(query::RayHit { t, rounding }) = query::ray_to(local, o, dir, (0.0, tmax), query::Side::Inside, 1e-9 * tmax.abs().max(1.0)) {
                    let p = [o[0] + t * dir[0], o[1] + t * dir[1], o[2] + t * dir[2]];
                    let g = ink::eval::eval_grads(raw, ctx, &[p])[0];
                    let normal = match t == 0.0 {
                        true => None,
                        false => ink::v3::direction([g[1], g[2], g[3]]),
                    };
                    let names: Vec<String> = ink::eval::provenance(raw, ctx, p).iter().map(|&t| tag_names[t as usize].clone()).collect();
                    let mut hit = json!({"target": label, "t": n(t), "point": v3(p), "normal": normal.map_or(J::Null, v3),
                        "name": names.last().cloned().unwrap_or_default(), "path": names});
                    if t == 0.0 {
                        hit["starts_inside"] = json!(true);
                    }
                    if let Some(w) = rounding {
                        hit["rounding"] = json!(n(w));
                    }
                    hits.push((t, hit));
                }
            }
            hits.sort_by(|a, b| a.0.total_cmp(&b.0));
            let hits: Vec<J> = hits.into_iter().map(|h| h.1).collect();
            json!({"from": v3(o), "dir": v3(dir), "max": n(tmax), "nearest": hits.first().cloned(), "hits": hits})
        })
        .collect();
    // Hits are sign changes found with interval pruning (never skipped); t
    // is the distance along the unit direction.
    let records = Records::new(a);
    let mut r = match a.get("rays") {
        // One ray, not as a table: its record is the output.
        None => records.picked(out, &[], &[])?.remove(0),
        Some(_) => json!({"rays": records.shape(out, &[])?}),
    };
    // A target not proven bounded: hits are searched in the window only.
    if open {
        r["unbounded"] = json!(true);
    }
    Ok(r)
}

pub(crate) fn cmd_mass(s: &mut Session, a: &Args) -> Result<J, InkError> {
    let items = s.items(a)?;
    let fixed = a.get("density").and_then(ink::lang::materials::density_option).transpose().map_err(|e| InkError::msg(format!("--density: {e}")))?;
    let asked = a.get("density").filter(|_| fixed.is_none());
    let res = a.cell("res")?;
    let shape = Records::new(a);
    let mut out = Vec::new();
    let mut total = [0.0; 10];
    // Mass is grams where the density is one (g/mm^3): a number, or the
    // bodies' density; else the attribute's own unit times mm^3.
    // Targets with no density have no mass, so the targets have no total.
    let mut without: Vec<(String, Vec<String>)> = Vec::new();
    for it in &items {
        // Density: a number, or the target's attribute (a field is fine;
        // ink::lang::materials::mass_density).
        let who = s.doc.named(it);
        let density_of = |attrs: &[(String, Value)], d: &mut ink::lang::Doc| -> Result<Option<(u32, String)>, InkError> {
            match fixed {
                Some(v) => Ok(Some((d.g.constant(v), String::new()))),
                None => ink::lang::materials::mass_density(d, &who, attrs, asked),
            }
        };
        // Where it comes from, which placing its attributes does not change.
        let density = density_of(&it.attrs, &mut s.doc)?;
        let label = match (fixed, &density) {
            (Some(v), _) => json!(v),
            (None, Some((_, what))) => json!(what),
            (None, None) => J::Null,
        };
        // Instanced: moments once per body in its own frame, then moved.
        let local = s.local(it);
        let found = if let Some(l) = local {
            let rho = density_of(&it.attrs, &mut s.doc)?.map(|d| d.0);
            let body = s.moments(l.shape, rho, res);
            let none = query::Measured { sums: vec![0.0; 11], spread: Vec::new(), res: 0.0, coarse: None, cells: [0; 2], settled: true };
            body.map_or(none, |mut b| {
                let m = l.at.moments(&b.moments());
                b.sums[..10].copy_from_slice(&m);
                b
            })
        } else {
            let (attrs, f) = (s.doc.placed_attrs(it), s.doc.placed_field(it));
            let rho = density_of(&attrs, &mut s.doc)?.map(|d| d.0);
            match s.moments(f, rho, res) {
                Some(r) => r,
                // No box: empty, or not proven bounded with nothing in the
                // window.
                None => {
                    out.push(match s.item_whole(it) {
                        true => json!({"target": it.label, "empty": true, "volume": {"estimate": 0}}),
                        false => json!({"target": it.label, "volume": {"estimate": null}, "unbounded": true}),
                    });
                    continue;
                }
            }
        };
        // The volume, an estimate at the cells the mass was measured at (as
        // integrals are: query::measured); of a shape not proven bounded,
        // none.
        let whole = s.item_whole(it);
        let volume = match whole {
            true => json!({"estimate": n(found.volume())}),
            false => json!({"estimate": null}),
        };
        let mut r = match density {
            Some(_) => {
                for (t, x) in total.iter_mut().zip(found.moments()) {
                    *t += x;
                }
                mass_json(Some(&found.moments()), &shape)
            }
            None => {
                // By the body it places (a shape that is no body: itself).
                let body = it.instance.as_ref().map(|i| i.0.clone()).unwrap_or_else(|| it.label.clone());
                match without.iter_mut().find(|(b, _)| *b == body) {
                    Some((_, targets)) => targets.push(it.label.clone()),
                    None => without.push((body, vec![it.label.clone()])),
                }
                mass_json(None, &shape)
            }
        };
        r["target"] = json!(it.label);
        r["volume"] = volume;
        r["density"] = label;
        measured_json(&mut r, &found);
        if !whole {
            r["unbounded"] = json!(true);
        }
        out.push(r);
    }
    // No total leaves a target out: a sum of the others would read as the
    // whole's mass and move its centre.
    let total = match without.is_empty() {
        true => mass_json(Some(&total), &shape),
        false => {
            let mut t = mass_json(None, &shape);
            // Once, not per target: a robot places one bare body hundreds
            // of times.
            t["without_mass"] = J::Object(without.into_iter().map(|(b, ts)| (b, json!(ts))).collect());
            t["no_mass"] = json!(NO_MASS);
            t
        }
    };
    let out = shape.picked(out, &["target"], &["empty", "unsettled", "unbounded"])?;
    Ok(json!({"targets": out, "total": total}))
}

/// Why the bodies in a mass answer's `without_mass` have no mass, and
/// what gives them one.
const NO_MASS: &str = "these bodies (each with the targets placing it) have no density: none of their own, no material with one, and no \
    published mass. Give the body, in the file that defines it, (attr material NAME) or (attr density 7.85g/cm3); or give an assembly \
    that places it the published mass of what it places, (attr mass M). One target alone weighs at a density you choose with \
    --target T --density 7.85g/cm3; a shape that is no body, only so";

/// What a mass answer gives only when --keys names it: the inertia tensor
/// about the centre of mass, its principal moments and axes.
const TENSOR: [&str; 3] = ["inertia", "principal_moments", "principal_axes"];

/// None: no mass, so neither a centre nor a tensor.
fn mass_json(m: Option<&[f64; 10]>, shape: &Records) -> J {
    let Some(m) = m else {
        let mut r = json!({"mass": null, "center_of_mass": null});
        for k in TENSOR.into_iter().filter(|k| shape.names(k)) {
            r[k] = J::Null;
        }
        return r;
    };
    let Some(p) = ink::rigid::mass_props(m) else { return json!({"mass": 0}) };
    let mut r = json!({"mass": n(p.mass), "center_of_mass": v3(p.center)});
    let tensor = [
        json!(p.inertia.iter().map(|r| v3(*r)).collect::<Vec<_>>()),
        json!(p.principal_moments.iter().map(|&v| n(v)).collect::<Vec<_>>()),
        json!(p.principal_axes.iter().map(|&a| v3(a)).collect::<Vec<_>>()),
    ];
    for (k, v) in TENSOR.into_iter().zip(tensor) {
        if shape.names(k) {
            r[k] = v;
        }
    }
    r
}

pub(crate) fn cmd_distance(s: &mut Session, a: &Args) -> Result<J, InkError> {
    // A target's bodies as one shape (of none, an empty assembly's, the
    // empty shape, as shape-of reads it: the distance query refuses it).
    let union = |s: &mut Session, t: &str| -> Result<(u32, String), InkError> {
        let items = s.doc.items(t)?;
        let shapes: Vec<u32> = items.iter().map(|it| s.doc.placed_field(it)).collect();
        Ok((s.doc.reduce(ink::graph::BinOp::Min, &shapes), t.to_string()))
    };
    let (fa, la) = union(s, a.get("a").ok_or_else(|| InkError::msg("distance needs --a TARGET --b TARGET"))?)?;
    let (fb, lb) = union(s, a.get("b").ok_or_else(|| InkError::msg("distance needs --a TARGET --b TARGET"))?)?;
    let res = a.cell("res")?;
    let (ra, rb) = (s.reach(fa), s.reach(fb));
    let finest = res.unwrap_or(f64::INFINITY).min(s.accuracy());
    ink::validity::refuse(&s.doc, &s.ctx, &[fa, fb], [(ra.lo, ra.hi), (rb.lo, rb.hi)], finest)?;
    let both = s.doc.g.max(fa, fb);
    // What the targets read that is not exact (sampled fields, solved
    // values not certified): whether they overlap is then not said, as no
    // claim on such inputs is guaranteed (claims::check); the bracket's
    // lower bound is proven only by guaranteed answers already.
    let inexact = s.doc.inexact(&[fa, fb], &s.certified());
    // A point inside both, exactly: a witness needs no accuracy, however
    // shallow the overlap (one in the window is one all the same).
    let (meet, all) = s.prove_with(both, &s.proof_options(res, 0.0));
    if meet.disproven() && !inexact.any() {
        return Ok(json!({"a": la, "b": lb, "overlapping": true, "distance": 0, "overlap": proof_json(&meet, &Names { all: &param_names(s), read: &[] }, None)}));
    }
    let (ta, tb) = (Tape::new(&s.doc.g, fa), Tape::new(&s.doc.g, fb));
    let whole = s.whole(fa) || s.whole(fb);
    let d = ink::claims::distance(s, fa, fb, ink::claims::ClaimOptions { res, ..Default::default() })?;
    let mut distance = json!({"lower_proven": n(d.lower), "upper_witness": d.upper.map(|u| n(u.0))});
    if let Some(u) = d.undecided {
        distance["undecided_at"] = n(u);
    }
    let (p, q) = d.upper.map_or(([f64::NAN; 3], [f64::NAN; 3]), |x| (x.1, x.2));
    let fa_names = names_at(s, &ta, p);
    let fb_names = names_at(s, &tb, q);
    // The true minimum distance lies in [lower_proven, upper_witness]:
    // lower_proven proven (every point of each target at least that from
    // the other; where the clearance at undecided_at could not be decided
    // at the finest cells the search reaches, that clearance), and the
    // distance between witness points confirmed inside their targets, to
    // the model's accuracy by one search of `boxes` boxes (stopped: by the
    // call's time limit).
    let boxes = d.boxes;
    // Apart where their intersection is proven empty, all of it, or a
    // clearance above 0 is proven; else (surfaces that coincide, or inputs
    // that are not exact) unknown.
    let overlapping = if !inexact.any() && ((meet.exact() && all) || d.lower > 0.0) { json!(false) } else { J::Null };
    let mut r = json!({"a": la, "b": lb, "overlapping": overlapping, "distance": distance, "boxes": boxes,
              "witness": {"on_a": v3(p), "on_b": v3(q), "name_a": fa_names.last(), "name_b": fb_names.last()}});
    if !whole {
        r["unbounded"] = json!(true);
    }
    inexact_fields(&mut r, &inexact, &param_names(s));
    Ok(r)
}

pub(crate) fn cmd_render(s: &mut Session, a: &Args) -> Result<J, InkError> {
    let out = a.get("out").map(str::to_string).unwrap_or(format!("out/{}.png", s.stem()));
    if !out.to_ascii_lowercase().ends_with(".png") {
        return Err(InkError::msg(format!("--out {out}: render writes a PNG; name the file .png")));
    }
    let cameras = match a.get("view") {
        None => render::Camera::sheet(),
        Some(v) if v.trim_start().starts_with('[') => {
            let from = parse_v3(a, "view")?;
            vec![render::Camera::looking_from(from).ok_or_else(|| InkError::msg("--view [x,y,z]: the direction the camera looks from, not zero"))?]
        }
        Some(v) => vec![render::Camera::named(v).ok_or_else(|| {
            let names: Vec<&str> = render::NAMED.iter().map(|n| n.0).collect();
            InkError::msg(format!("--view {v}: a view is {} or [x,y,z], the direction the camera looks from", names.join(", ")))
        })?],
    };
    // Without a target or a root, the file's first assembly that no other
    // places: its bodies as assembled, not every body at its origin.
    let tops: Vec<String> = match a.get("target").is_none() && s.doc.body("root").is_none() && s.doc.assembly("root").is_none() {
        true => s.doc.top_assemblies().into_iter().map(String::from).collect(),
        false => Vec::new(),
    };
    let mut items = match tops.first() {
        Some(t) => s.doc.items(t)?,
        None => s.items(a)?,
    };
    let cut = match a.get("cut") {
        Some(c) => Some(cut_items(s, &mut items, c, cameras[0].eye)?),
        None => None,
    };
    let size = (a.f64("size")?.unwrap_or(1024.0) as usize).clamp(64, 4096);
    // Meshes made again no finer than half a pixel: what is a pixel thick
    // is drawn whole, and what is thinner is drawn dotted at any res (a pixel
    // shows what is at its centre).
    let all: Vec<(V3, V3, f64)> = items.iter().filter_map(|it| s.item_bounds(it)).collect();
    let finest = union_bounds(&all).map_or(0.0, |(lo, hi)| render::mm_per_px(lo, hi, &cameras, size) / 2.0);
    // The meshes a GLB shows (cached: a part rendered and shared is meshed
    // once), each vertex with how far it may be from the surface.
    let (meshed, _) = mesh_items(s, a, &items, MeshHow { finest, errors: true, ..SEEN })?;
    let worlds: Vec<Mesh> = meshed.iter().map(|m| m.world()).collect();
    // The field's normals at the corners, turned as the instance is.
    let normals: Vec<Vec<[f32; 3]>> = meshed
        .iter()
        .map(|m| match &m.place {
            None => m.normals.to_vec(),
            Some(rg) => m.normals.iter().map(|n| std::array::from_fn(|i| (0..3).map(|j| rg.r[i][j] * n[j] as f64).sum::<f64>() as f32)).collect(),
        })
        .collect();
    let transmission: Vec<Vec<f32>> = meshed
        .iter()
        .map(|m| {
            let it = items.iter().find(|it| it.label == m.label).unwrap();
            let attrs = match m.place {
                Some(_) => it.attrs.clone(),
                None => s.doc.placed_attrs(it),
            };
            Look::of(&attrs).transmission.map_or_else(Vec::new, |f| paint(s, [Some(f)], [0.0], &m.paint_at).concat())
        })
        .collect();
    let objs: Vec<render::Object> = (0..meshed.len())
        .map(|k| render::Object { mesh: &worlds[k], colors: &meshed[k].colors, normals: &normals[k], transmission: &transmission[k], error: &meshed[k].quality.error })
        .collect();
    let (img, w, hh, views) = render::render(&objs, &cameras, size, cut.as_ref().map(|c| (c.axis, c.at)));
    export::write(&out, &export::png(w, hh, &img)).map_err(|e| InkError::msg(format!("cannot write {out}: {e}")))?;
    let vj: Vec<J> = views.iter().map(|v| json!({"view": v.name, "camera_from": v3(v.eye)})).collect();
    let legend: Vec<J> = meshed
        .iter()
        .map(|m| {
            let c = m.colors.first().copied().unwrap_or([0.5; 3]);
            json!({"target": m.label, "color": format!("#{:02X}{:02X}{:02X}", (c[0]*255.0) as u8, (c[1]*255.0) as u8, (c[2]*255.0) as u8),
                   "res": n(m.quality.res), "triangles": m.mesh.tris.len(), "unresolved": unresolved_json(m)})
        })
        .collect();
    let layout = match cameras.len() {
        1 => "one view".to_string(),
        _ => "2x2: iso top-left, top top-right, front bottom-left, right bottom-right".to_string(),
    };
    // Drawn from a mesh at each target's res (legend): features thinner
    // than res may merge or vanish. render::render hatches a face it finds
    // on the cut plane, but finds none of a cut target's (a box cut through
    // its middle shows no hatching), so docs/IR.md does not promise it.
    let mut r = json!({"out": out, "width": w, "height": hh, "layout": layout, "views": vj, "legend": legend});
    if let Some(t) = tops.first() {
        r["target"] = json!(t);
        if tops.len() > 1 {
            r["other_assemblies"] = json!(tops[1..]);
            r["hint"] = json!(format!("No root, so this is '{t}', the first of the file's {} assemblies no other places. --target NAME draws another", tops.len()));
        }
    }
    if let Some(c) = cut {
        r["cut"] = json!({"axis": AXES[c.axis], "at": n(c.at), "removed": c.removed});
    }
    let warnings = mesh_warnings(&meshed);
    if !warnings.is_empty() {
        r["warnings"] = json!(warnings);
    }
    Ok(r)
}

const AXES: [&str; 3] = ["x", "y", "z"];

/// A cut (`cut_items`): its axis, where, and which half went ("z > 25").
struct Cut {
    axis: usize,
    at: f64,
    removed: String,
}

/// `--cut AXIS[=AT]`: each target cut by the plane AXIS = AT (by default
/// through the middle of the targets), the half on the camera's side
/// removed (along an axis the camera looks across, iso's: +x, -y, +z).
/// A target wholly on one side is kept as it is or left out (by its
/// bounds, which are rounded outward, or `no_material` on the side kept:
/// its face on the plane); one the plane crosses is its world field and
/// the half-space's max.
fn cut_items(s: &mut Session, items: &mut Vec<Item>, spec: &str, eye: V3) -> Result<Cut, InkError> {
    use ink::graph::Var;
    let usage = || InkError::msg(format!("--cut {spec}: give an axis x, y or z and, optionally, where: z=25 (else the middle)"));
    let (name, at) = spec.split_once('=').map_or((spec, None), |(n, v)| (n, Some(v)));
    let axis = AXES.iter().position(|&n| n == name.trim()).ok_or_else(usage)?;
    let bounds: Vec<Option<(V3, V3, f64)>> = items.iter().map(|it| s.item_bounds(it)).collect();
    let at = match at {
        Some(v) => ink::sexp::number_of(v.trim(), ink::units::Quantity::Length, "--cut").map_err(InkError::msg)?,
        None => {
            // Each box is rounded outward by up to its cell: the targets'
            // ends lie in lo..lo + cell and hi - cell..hi, their middle
            // between those ends' middles, written as shortly as that
            // allows.
            let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
            for (l, h, c) in bounds.iter().flatten() {
                lo = [lo[0].min(l[axis]), lo[1].min(l[axis] + c)];
                hi = [hi[0].max(h[axis] - c), hi[1].max(h[axis])];
            }
            if lo[0] > hi[1] {
                return Err(InkError::msg("nothing to cut: every target is empty"));
            }
            shortest((lo[0] + hi[0]) / 2.0, (lo[1] + hi[1]) / 2.0)
        }
    };
    let side = match eye[axis] {
        e if e.abs() > 1e-9 => e.signum(),
        _ => render::Camera::named("iso").unwrap().eye[axis].signum(),
    };
    let v = s.doc.g.var([Var::X, Var::Y, Var::Z][axis]);
    let c = s.doc.g.constant(at);
    // Positive on the removed side.
    let half = if side > 0.0 { s.doc.g.sub(v, c) } else { s.doc.g.sub(c, v) };
    let kept: Vec<Item> = std::mem::take(items)
        .into_iter()
        .zip(bounds)
        .filter_map(|(mut it, b)| {
            // How far into the removed half its nearest and farthest material
            // reach, and its box on the side kept.
            let reach = b.map(|(mut lo, mut hi, _)| {
                let (p, q) = (side * (lo[axis] - at), side * (hi[axis] - at));
                if side > 0.0 {
                    hi[axis] = at
                } else {
                    lo[axis] = at
                }
                (p.min(q), p.max(q), lo, hi)
            });
            if let Some((near, far, _, _)) = reach
                && (far <= 0.0 || near >= 0.0)
            {
                return (far <= 0.0).then_some(it);
            }
            let f = s.doc.placed_field(&it);
            if let Some((_, _, lo, hi)) = reach
                && no_material(&Tape::new(&s.doc.g, f), &s.ctx, lo, hi)
            {
                return None;
            }
            (it.shape, it.attrs) = (s.doc.g.max(f, half), s.doc.placed_attrs(&it));
            it.instance = None;
            Some(it)
        })
        .collect();
    if kept.is_empty() {
        return Err(InkError::msg(format!("--cut {spec} removes every target")));
    }
    *items = kept;
    Ok(Cut { axis, at, removed: format!("{} {} {at}", AXES[axis], if side > 0.0 { '>' } else { '<' }) })
}

/// The number in lo..hi with the fewest decimals, nearest its middle.
fn shortest(lo: f64, hi: f64) -> f64 {
    let mid = (lo + hi) / 2.0;
    let wide = ((hi - lo).log10().ceil() as i32).clamp(-17, 17);
    (-wide..=17)
        .map(|d| if d >= 0 { (mid * 10f64.powi(d)).round() / 10f64.powi(d) } else { (mid / 10f64.powi(-d)).round() * 10f64.powi(-d) })
        .find(|x| (lo..=hi).contains(x))
        .unwrap_or(mid)
}

/// Whether interval bounds prove `tape` >= 0, no material, over lo..hi:
/// over the box, else over its eighths where it does not settle, down to a
/// 16th of it across. A face on the box's side is 0 there, not material.
fn no_material(tape: &Tape, ctx: &ink::eval::Ctx, lo: V3, hi: V3) -> bool {
    use ink::interval::I;
    let mut boxes = vec![(lo, hi, 0)];
    while let Some((lo, hi, depth)) = boxes.pop() {
        let r = ink::eval::eval_interval(tape, ctx, std::array::from_fn(|a| I::new(lo[a], hi[a])), None);
        if r.lo >= 0.0 {
            continue;
        }
        if r.hi < 0.0 || depth == 4 {
            return false;
        }
        let mid: V3 = std::array::from_fn(|a| 0.5 * (lo[a] + hi[a]));
        for k in 0..8 {
            let pick = |a: usize, low: V3, high: V3| if k >> a & 1 == 0 { low[a] } else { high[a] };
            boxes.push((std::array::from_fn(|a| pick(a, lo, mid)), std::array::from_fn(|a| pick(a, mid, hi)), depth + 1));
        }
    }
    true
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

/// Ramer-Douglas-Peucker simplification of a (possibly closed) polyline.
fn simplify_polyline(pts: &[[f64; 2]], closed: bool, tol: f64) -> Vec<[f64; 2]> {
    if pts.len() < 4 {
        return pts.to_vec();
    }
    fn rdp(p: &[[f64; 2]], tol: f64, out: &mut Vec<[f64; 2]>) {
        let (a, b) = (p[0], p[p.len() - 1]);
        let d = [b[0] - a[0], b[1] - a[1]];
        let len = (d[0] * d[0] + d[1] * d[1]).sqrt();
        let dist = |q: [f64; 2]| {
            if len == 0.0 { ((q[0] - a[0]).powi(2) + (q[1] - a[1]).powi(2)).sqrt() } else { ((q[0] - a[0]) * d[1] - (q[1] - a[1]) * d[0]).abs() / len }
        };
        let (k, m) = (1..p.len() - 1).map(|k| (k, dist(p[k]))).fold((0, 0.0), |b, x| if x.1 > b.1 { x } else { b });
        if m > tol {
            rdp(&p[..=k], tol, out);
            out.pop();
            rdp(&p[k..], tol, out);
        } else {
            out.push(a);
            out.push(b);
        }
    }
    let mut out = Vec::new();
    if closed {
        let far = (1..pts.len())
            .max_by(|&i, &j| {
                let d = |k: usize| (pts[k][0] - pts[0][0]).powi(2) + (pts[k][1] - pts[0][1]).powi(2);
                d(i).total_cmp(&d(j))
            })
            .unwrap();
        let mut ring = pts.to_vec();
        ring.push(pts[0]);
        rdp(&ring[..=far], tol, &mut out);
        out.pop();
        rdp(&ring[far..], tol, &mut out);
        out.pop();
    } else {
        rdp(pts, tol, &mut out);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Whether two targets overlap is said only on exact inputs: unit
    /// balls at x = 0 and x = u, u solved from 1e-13 u = 0 starting at 4,
    /// where the solve stops (its exact root is 0, where they overlap).
    /// Their intersection, empty at u = 4, was taken for a proof that they
    /// do not.
    #[test]
    fn distance_says_nothing_of_overlap_on_inexact_inputs() {
        let src = "(import \"std\")
(solve s (unknowns (u 4)) (== (* 1e-13 u) 0))
(def a (sphere 1))
(def b (move (sphere 1) u 0 0))";
        let doc = ink::load_str(src, "distance", Path::new(".")).unwrap();
        let ctx = doc.ctx(&[], 0.0).unwrap();
        let mut s = ModelView::new(doc, ctx);
        let args = Args::parse(&["distance", "unused.ink", "--a", "a", "--b", "b"].map(str::to_string));
        let answer = cmd_distance(&mut s, &args).unwrap();
        assert!(answer["overlapping"].is_null() && answer["solved_inputs"] == json!(["u"]), "the solve stopped at 4; the exact root is 0, where the spheres overlap: {answer}");
    }

    /// A body with no inside (move's arguments swapped: its field is 0
    /// over its box) has no pair: the distance says so at once, without
    /// meshing the body first, finer each time to 32,768 cells a side (past
    /// 30 s).
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
        let answer = cmd_distance(&mut s, &args).unwrap();
        assert!(answer["distance"]["upper_witness"].is_null() && answer["distance"]["lower_proven"].is_number(), "{answer}");
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
        assert!(
            w.len() == 1 && w[0].starts_with("'p10': its 41 frames at res 1.681792831 over 100.0 mm") && w[0].ends_with("up to res 0.5 (without it, res 0.2973017788)"),
            "{w:?}"
        );
        assert!(warned(&small).is_empty());
        // The ladder rounds up to quarter powers of two.
        for (h, want) in [(0.25, 0.25), (0.2, 0.21022410381342863), (0.42, 0.42044820762685725), (0.43, 0.5), (1.0, 1.0), (1.19, std::f64::consts::SQRT_2), (0.0568, 0.0625)] {
            assert!((ladder(h) - want).abs() < 1e-12 && ladder(h) >= h * (1.0 - 1e-9), "{h}: {}", ladder(h));
        }
        assert_eq!(ladder(ladder(0.3)), ladder(0.3), "a rung is its own");
        assert_eq!(ladder(ladder(0.3) / 2.0), ladder(0.3) / 2.0, "half a rung is a rung");
    }
}
