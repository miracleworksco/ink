//! PNG rendering and optional section cuts.

use super::budget::SEEN;
use super::mesh::{MeshHow, mesh_items};
use super::mesh_output::{left_out_json, mesh_warnings, unresolved_json};
use super::spatial::parse_v3;
use super::target::{bounded_local, whole_placed_field};
use crate::{Arc, Args, CliView, InkError, Item, J, Mesh, Session, Tape, V3, export, json, n, paint, render, union_bounds, v3};
use ink::lang::look::Look;
use rayon::prelude::*;

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
    // Meshes for the picture, each vertex with how far it may be from the
    // surface: sized by their boxes, as a picture shows them whole (a
    // GLB's by their surfaces, for a viewer's zoom), within a GLB's budget.
    let (meshed, _) = mesh_items(s, a, &items, MeshHow { finest, by_area: false, ..SEEN })?;
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
    // Each opaque target's material, by its field: where two targets'
    // surfaces meet, the one listed first shows (render::Object).
    let fields: Vec<Option<(Arc<Tape>, V3, V3)>> = meshed
        .iter()
        .zip(&transmission)
        .map(|(m, tr)| {
            let it = items.iter().find(|it| it.label == m.label).unwrap();
            let (lo, hi, _) = s.item_bounds(it).filter(|_| tr.is_empty())?;
            let f = s.doc.placed_field(it);
            Some((Arc::new(Tape::new(&s.doc.g, f)), lo, hi))
        })
        .collect();
    let ctx = &s.ctx;
    type Field<'a> = Box<dyn Fn(&[V3]) -> Vec<[f64; 4]> + Sync + 'a>;
    let reads: Vec<Option<Field>> = fields
        .iter()
        .map(|f| {
            let (tape, _, _) = f.as_ref()?;
            let tape = tape.clone();
            let d: Field = Box::new(move |p: &[V3]| p.par_chunks(4096).flat_map_iter(|c| ink::eval::eval_grads(&tape, ctx, c)).collect());
            Some(d)
        })
        .collect();
    let objs: Vec<render::Object> = (0..meshed.len())
        .map(|k| render::Object {
            mesh: &worlds[k],
            colors: &meshed[k].colors,
            normals: &normals[k],
            transmission: &transmission[k],
            solid: fields[k].as_ref().zip(reads[k].as_deref()).map(|((_, lo, hi), field)| render::Solid {
                lo: *lo,
                hi: *hi,
                field,
                reach: 2.0 * meshed[k].quality.res + meshed[k].off,
            }),
        })
        .collect();
    // Where each mesh may leave material out or close a gap (its
    // unresolved regions), in the world: outlined in the picture.
    let marks: Vec<(V3, V3)> = meshed.iter().flat_map(|m| m.quality.unresolved.regions.iter().map(move |r| m.place.map_or((r.lo, r.hi), |p| p.aabb(r.lo, r.hi)))).collect();
    let (img, w, hh, views) = render::render(&objs, &cameras, size, cut.as_ref().map(|c| (c.axis, c.at)), &marks);
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
    // What the picture leaves out, beside the picture (its red boxes there).
    let left = left_out_json(&meshed);
    if !left.is_empty() {
        r["left_out"] = json!(left);
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
    // A cut may discard a whole item by its box before mesh planning.
    let mut bounds = Vec::with_capacity(items.len());
    for it in items.iter() {
        bounded_local(s, it)?;
        bounds.push(s.item_bounds(it));
    }
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
    let mut kept = Vec::new();
    for (mut it, b) in std::mem::take(items).into_iter().zip(bounds) {
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
            if far <= 0.0 {
                kept.push(it);
            }
            continue;
        }
        let f = whole_placed_field(s, &it)?;
        if let Some((_, _, lo, hi)) = reach
            && no_material(&Tape::new(&s.doc.g, f), &s.ctx, lo, hi)
        {
            continue;
        }
        (it.shape, it.attrs) = (s.doc.g.max(f, half), s.doc.placed_attrs(&it));
        it.instance = None;
        kept.push(it);
    }
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
        if r.none_inside() {
            continue;
        }
        if r.all_inside() || depth == 4 {
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
