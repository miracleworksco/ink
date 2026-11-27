//! Point and ray queries, bounds, and distance between targets.

use super::left_out;
use super::mesh_output::regions_json;
use super::target::whole_placed_field;
use crate::{
    Arc, Args, CliView, InkError, J, Names, Reach, Records, Session, Tape, V3, Value, eval_dual, inexact_fields, json, json_arg, json_v3, max_extent, n, param_names, parse_points,
    proof_json, query, union_bounds, v3,
};

pub(crate) fn cmd_eval(s: &mut Session, a: &Args) -> Result<J, InkError> {
    let pts = parse_points(a)?;
    let items = s.items(a)?;
    let shape = Records::new(a);
    let fields: Vec<u32> = items.iter().map(|it| s.doc.placed_field(it)).collect();
    let mut out = Vec::new();
    for (it, &f) in items.iter().zip(&fields) {
        let tape = Tape::new(&s.doc.g, f);
        let attrs = if shape.wants("attrs") { s.doc.placed_attrs(it) } else { Default::default() };
        let res: Vec<J> = pts
            .iter()
            .map(|&p| {
                let d = eval_dual(&tape, &s.ctx, p);
                // Each rate is a partial, every other param held: no form gives one through what it sets.
                // Inside by the field's outward-rounded bounds at the point, as a probe claim reads it: where
                // they reach across 0 (a float that rounds to the other side), not known, and the bounds given.
                // At a 0 the sign rule reads a cube about the point.
                let b = ink::eval::eval_interval(&tape, &s.ctx, p.map(ink::interval::I::point), None);
                let band = (d[0] != 0.0 && b.lo <= 0.0 && b.hi >= 0.0 && b.lo < b.hi).then_some(b);
                let inside = if band.is_some() { None } else { ink::sign::at(&tape, &s.ctx, p, d[0]) };
                let mut r = json!({"point": v3(p), "value": n(d[0]), "inside": inside, "grad": v3([d[1], d[2], d[3]]), "d_dt": n(d[4])});
                if let Some(b) = band {
                    r["band"] = json!([b.lo, b.hi]);
                }
                if shape.wants("d_dparams") {
                    r["d_dparams"] = (0..s.doc.params.len()).map(|k| (s.doc.params[k].name.clone(), n(d[5 + k]))).collect::<serde_json::Map<_, _>>().into();
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

pub(super) fn parse_v3(a: &Args, key: &str) -> Result<V3, InkError> {
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

pub(crate) fn cmd_components(s: &mut Session, a: &Args) -> Result<J, InkError> {
    // Each target's shells, from a mesh (an estimate), and its cavities
    // counted as the `cavities` claim counts them (ModelView::components).
    let res = a.cell("res")?;
    let mut out = Vec::new();
    let mut warnings = Vec::new();
    for it in s.items(a)? {
        let field = whole_placed_field(s, &it)?;
        let Some(m) = s.components(field, res).map_err(|e| InkError { message: format!("'{}': {}", it.label, e.message), ..e })? else {
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
        if let (true, Some(at)) = (m.unresolved.cells > 0, m.unresolved.at) {
            warnings.push(format!("{}: {}", it.label, left_out(&m.unresolved, m.res, at)));
        }
        let c = &m.voids;
        if let Some((what, decide)) = c.why_undecided() {
            let at = c.at.map_or(String::new(), |p| format!(" (one at [{:.3}, {:.3}, {:.3}])", p[0], p[1], p[2]));
            let most = c.at_most.map_or("unbounded".into(), |n| n.to_string());
            warnings.push(format!("{}: its cavities are proven between {} and {most}: {what}{at}. To decide it: {decide}", it.label, c.at_least));
        }
        out.push(json!({"target": it.label, "res": n(m.res), "solids": solids, "cavities": c.count(),
                        "proven": {"cavities": {"at_least": c.at_least, "at_most": c.at_most, "at": c.at.map(v3), "cells": c.cells, "finest": n(c.finest),
                                               "undecided": c.undecided.map(|(a, b)| json!([n(a), n(b)])), "ended": c.ended,
                                               "accuracy": n(c.accuracy), "set_aside": (c.thin > 0).then(|| json!({"cells": c.thin, "at": c.thin_at.map(v3)}))}},
                        "watertight": an.boundary_edges == 0 && an.nonmanifold_edges == 0, "shells": shells,
                        "unresolved": regions_json(&m.unresolved, None, m.res)}));
    }
    let mut r = json!({"targets": out});
    if !warnings.is_empty() {
        r["warnings"] = json!(warnings);
    }
    Ok(r)
}
