//! Cross-section queries and DXF outlines.

use super::target::whole_placed_field;
use crate::{Arc, Args, CliView, InkError, J, Records, Session, Tape, V3, export, json, measured, n, query, union_bounds, v3};

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
    let fields: Vec<u32> = items.iter().map(|it| whole_placed_field(s, it)).collect::<Result<_, _>>()?;
    let bounds: Vec<_> = fields.iter().map(|&f| s.bounds(f)).collect();
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
    // Samples on 0 no bound decides, drawn as air (sign::drawn).
    let mut undecided: Vec<V3> = Vec::new();
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
        let (loops, open) = query::slice_checked(&tape, &s.ctx, axis, at, *lo, *hi, h);
        undecided.extend(open);
        for l in loops {
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
            let (loops, open) = query::slice_mark(&gt, &rt, &s.ctx, axis, at, *lo, *hi, h);
            undecided.extend(open);
            for l in loops {
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
    if let Some(p) = undecided.first() {
        r["unresolved"] = json!({"samples": undecided.len(), "at": v3(*p)});
    }
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
