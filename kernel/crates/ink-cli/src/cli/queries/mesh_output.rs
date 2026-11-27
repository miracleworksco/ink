//! Mesh export, shell queries, and reports on mesh quality.

use super::budget::SEEN;
use super::mesh::{MeshHow, Meshed, looks, mesh_items, motion};
use super::target::EMPTY_CAUSE;
use crate::{Arc, Args, CliView, HashMap, InkError, Instant, Item, J, Mesh, Rigid, Session, Tape, V3, Value, bounds, export, json, measured, mesh, n, v3};
use ink::graph::Var;

/// `ink mesh`'s default accuracy, a fraction of res: 0.01 mm at res 0.3.
const EXPORT_TOL_PER_RES: f64 = 1.0 / 30.0;

/// Each distinct mesh that may leave material out or close a gap (its
/// unresolved regions): its target, how many regions, and what it may
/// leave out (`left_out`), for whatever shows the meshes to say so.
pub fn left_out_json(meshed: &[Meshed]) -> Vec<J> {
    let mut seen = std::collections::HashSet::new();
    meshed
        .iter()
        .filter(|m| !m.quality.unresolved.regions.is_empty() && seen.insert(Arc::as_ptr(&m.mesh) as usize))
        .map(|m| {
            let u = &m.quality.unresolved;
            let at = u.at.unwrap_or(u.regions[0].at);
            let at = m.place.map_or(at, |p| p.apply(at));
            json!({"target": m.label, "regions": u.regions.len(), "what": left_out(u, m.quality.res, at)})
        })
        .collect()
}

/// Warnings, one per distinct mesh: a target meshed to nothing (empty, or
/// thinner than res everywhere), and sign changes a mesh left out.
pub(super) fn mesh_warnings(meshed: &[Meshed]) -> Vec<String> {
    let mut seen = std::collections::HashSet::new();
    meshed
        .iter()
        .filter(|m| seen.insert(Arc::as_ptr(&m.mesh) as usize))
        .filter_map(|m| {
            if m.mesh.tris.is_empty() {
                return Some(format!("'{}': its mesh at res {} has no triangles. It is empty, or thinner than res everywhere. {EMPTY_CAUSE}", m.label, measured(m.quality.res)));
            }
            let q = &m.quality;
            let at = match q.unresolved.at {
                Some(at) => at,
                None if q.unresolved.refused > 0.0 => m.mesh.verts.first().copied()?,
                None => return None,
            };
            let at = m.place.map_or(at, |p| p.apply(at));
            Some(format!("'{}': {}", m.label, left_out(&q.unresolved, q.res, at)))
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

/// What a mesh at res `h` leaves out (mesh::Unresolved), `at` a point of
/// it in the world.
pub fn left_out(u: &mesh::Unresolved, h: f64, at: V3) -> String {
    let found = u.cells - u.undecided - u.untried;
    let mut what = Vec::new();
    if u.refused > 0.0 {
        what.push(format!("made at res {}, not {}: its field's bounds left more of its box to sample there than the budget pays for", measured(h), measured(u.refused)));
    }
    match u.thinnest {
        t if t > 0.0 && t < h => what.push(format!("{} cells hold material thinner than res that it leaves out ({} across at the thinnest)", u.material, measured(t))),
        _ if u.material > 0 => what.push(format!("{} cells hold material that it leaves out", u.material)),
        _ => {}
    }
    if found > u.material {
        what.push(format!("{} cells hold a gap thinner than res that it closes", found - u.material));
    }
    if u.undecided > 0 {
        what.push(format!(
            "in {}{} cells the field's bounds could not rule out material or a gap the mesh leaves out, thinner than {} (unresolved.regions each say how thin)",
            if u.untried > 0 { "at least " } else { "" },
            u.undecided,
            measured(u.thinner_than)
        ));
    }
    if u.away > 0 {
        what.push(format!("{} of them away from its drawn surface", u.away));
    }
    if u.seams > 0 {
        what.push(format!(
            "{}{} where a select switches between two fields (a seam, as at a sweep's mitres), which bounds cannot show continuous",
            if u.untried > 0 { "at least " } else { "" },
            u.seams
        ));
    }
    let advice = if found > 0 { "it has holes or misses features there; use a smaller --res" } else { "it may miss features there; a smaller --res searches finer" };
    format!("its mesh at res {}: {} (e.g. at [{}, {}, {}]); {advice}", measured(h), what.join("; "), measured(at[0]), measured(at[1]), measured(at[2]))
}

pub(super) fn unresolved_json(m: &Meshed) -> J {
    let u = &m.quality.unresolved;
    let mut j = regions_json(u, m.place.as_ref(), m.quality.res);
    j["at"] = json!(u.at.map(|p| v3(m.place.map_or(p, |r| r.apply(p)))));
    j
}

/// A mesh's `Unresolved` at res `h` as JSON, its regions' boxes and points
/// placed by `place` (a box turned: the box round it), and how near what it
/// shows drawn lies to its vertices (mesh: 2 sqrt 3 h, or more where a
/// region drawn whole may stray farther).
pub(super) fn regions_json(u: &mesh::Unresolved, place: Option<&Rigid>, h: f64) -> J {
    let pt = |p: V3| v3(place.map_or(p, |r| r.apply(p)));
    let regions: Vec<J> = u
        .regions
        .iter()
        .map(|r| {
            let (lo, hi) = place.map_or((r.lo, r.hi), |p| p.aabb(r.lo, r.hi));
            json!({"lo": v3(lo), "hi": v3(hi), "cells": r.cells, "material": r.material, "gaps": r.gaps, "undecided_at_least": r.undecided, "untried": r.untried, "beside": r.beside, "seams_at_least": r.seams, "thinner_than": n(r.thinner_than),
                   "searched": n(r.searched), "thinnest": n(r.thinnest), "at": pt(r.at), "found": r.found, "groups": r.groups})
        })
        .collect();
    json!({"cells": u.cells, "material": u.material, "thinnest": n(u.thinnest), "undecided_at_least": u.undecided, "untried": u.untried, "beside": u.beside, "away": u.away, "seams_at_least": u.seams, "thinner_than": n(u.thinner_than), "searched": n(u.searched), "refused": n(u.refused),
           "drawn_within": n((2.0 * 3f64.sqrt() * h).max(u.whole_within)), "drawn_whole": {"regions": u.whole, "within": n(u.whole_within)}, "at": u.at.map(pt), "regions": regions})
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
    let how = if glb { SEEN } else { MeshHow { div: 256.0, tol_per_res: EXPORT_TOL_PER_RES, measure: true, budget: 0.0, by_area: false, finest: 0.0, seen: false } };
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
    // For a viewer to say what the meshes leave out.
    let left = left_out_json(meshed);
    if !left.is_empty() {
        info["left_out"] = json!(left);
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
