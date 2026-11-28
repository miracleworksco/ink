//! A solve's setup receipt: where each support and load landed, what each
//! holds or pushes, the material and where each property came from, and
//! the bodies left out, so a setup is checked in a minute (a force on the
//! wrong end, a clamp where the bolt is not, reads at a glance). `ink
//! solve` gives it in its JSON, and in words with `--plain`.

use crate::{Boundary, Footprint, SolvedStudy, Study, StudyPhysics};
use ink_core::v3 as vec3;
use ink_lang::command::{n, v3};
use ink_lang::view::ModelView;
use serde_json::{Value as J, json};

/// At most this many bodies are each checked against the study's shape.
const MAX_BODIES: usize = 200;
/// Points a side of a body's box sampled to find whether the study holds it.
const BODY_SAMPLES: usize = 16;

/// The receipt, into `out` (`study_summary`'s): each load's and fix's place,
/// each fix's directions and share, the material's sources, and the bodies.
pub(crate) fn add(out: &mut J, s: &mut ModelView, sv: &SolvedStudy, name: &str) {
    let r = &sv.sol.report;
    let st = s.doc.form::<Study>(name).expect("a solved study");
    let heat = matches!(sv.physics, StudyPhysics::Heat { .. });
    // (the part's centroid: what moments are also given about)
    let center = r.centroid;
    let about_center = |force: [f64; 3], moment: [f64; 3]| v3(vec3::sub(moment, vec3::cross(center, force)));
    if !heat {
        out["part_centroid_mm"] = v3(center);
    }
    // Loads: gravity's first, where the study has it.
    let body = r.applied.len().saturating_sub(st.loads.len());
    for (i, (_, total)) in r.applied.iter().enumerate() {
        let Some(load) = out["loads"].get_mut(i) else { break };
        load["kind"] = json!(match i.checked_sub(body).and_then(|k| st.loads.get(k)) {
            None => "gravity",
            Some(_) if heat => "power",
            Some(l) if l.pressure.is_some() => "pressure",
            Some(l) if l.torque.is_some() => "torque",
            Some(_) => "force",
        });
        if let Some(Some(f)) = r.load_footprints.get(i) {
            load["landed"] = landed(f);
        }
        if !heat && let Some(m) = r.applied_moments.get(i) {
            load["moment_about_part_centroid_nmm"] = about_center(*total, *m);
        }
    }
    // Fixes: what each holds, where, and its share of all the reactions.
    let size = |v: &[f64; 3]| if heat { v[0].abs() } else { vec3::norm(*v) };
    let all: f64 = r.reactions.iter().map(|(_, v)| size(v)).sum();
    for (j, (fx, (_, total))) in st.fixes.iter().zip(&r.reactions).enumerate() {
        let Some(rx) = out["reactions"].get_mut(j) else { break };
        let (dofs, holds) = holds(s, fx, heat);
        rx["dofs"] = json!(dofs);
        rx["holds"] = json!(holds);
        rx["landed"] = match (fx.boundary, r.fix_footprints.get(j)) {
            (Boundary::Mirror { axis, at }, _) => json!({"plane": format!("{} = {}", ["x", "y", "z"][axis], words::g(at))}),
            (_, Some(Some(f))) => landed(f),
            _ => json!({"unavailable": "its region meets none of the part's surface or material"}),
        };
        rx["share"] = if all > 0.0 { n(size(total) / all) } else { J::Null };
        if !heat && let Some(m) = r.reaction_moments.get(j) {
            rx["moment_about_part_centroid_nmm"] = about_center(*total, *m);
        }
    }
    // The material, and where each property came from.
    let num = |f: Option<crate::graph::NodeId>| f.map(|f| s.doc.g.as_const(f).map_or(json!("a field"), n));
    let from: serde_json::Map<String, J> = st.material_from.iter().map(|(k, v)| (k.clone(), json!(v))).collect();
    match sv.physics {
        StudyPhysics::Elastic { density, .. } => {
            if let (Some(d), J::Object(m)) = (num(density), &mut out["material"]) {
                m.insert("density_kg_m3".into(), d);
            }
        }
        StudyPhysics::Heat { k, density, capacity } => {
            let mut m = json!({"k_w_mk": num(Some(k))});
            if let Some(d) = num(density) {
                m["density_kg_m3"] = d;
            }
            if let Some(c) = capacity {
                m["capacity_j_kgk"] = n(c);
            }
            out["material"] = m;
        }
    }
    if let J::Object(m) = &mut out["material"] {
        m.insert("from".into(), J::Object(from));
    }
    if let Some(b) = bodies(s, sv) {
        out["bodies"] = b;
    }
}

/// A footprint as JSON: where, how big, and a pressure's range.
fn landed(f: &Footprint) -> J {
    let mut o = json!({"at": v3(f.at), "lo": v3(f.lo), "hi": v3(f.hi)});
    if f.area > 0.0 {
        o["area_mm2"] = n(f.area);
    }
    if f.volume > 0.0 {
        o["volume_mm3"] = n(f.volume);
    }
    if let Some([lo, mean, hi]) = f.pressure {
        o["pressure_mpa"] = json!({"min": n(lo), "mean": n(mean), "max": n(hi)});
    }
    o
}

/// What a fix holds: its directions (`x`, `y`, `z`, `normal`, `temperature`,
/// `convection`), and in words, with what it leaves free.
fn holds(s: &ModelView, fx: &crate::StudyFix, heat: bool) -> (Vec<&'static str>, String) {
    let num = |f: crate::graph::NodeId| s.doc.g.as_const(f);
    if heat {
        let t = num(fx.value[0]).map_or("a field".to_string(), |t| format!("{} C", words::g(t)));
        return match fx.film {
            Some(h) => (vec!["convection"], format!("convection at {} W/(m2 K) to {t}", words::g(h))),
            None => (vec!["temperature"], format!("held at {t}")),
        };
    }
    if let Boundary::Mirror { axis, at } = fx.boundary {
        let a = ["x", "y", "z"][axis];
        return (vec![a], format!("a cut of symmetry: the plane {a} = {} held along {a}", words::g(at)));
    }
    let dofs: Vec<&'static str> = if fx.normal { vec!["normal"] } else { ["x", "y", "z"].into_iter().zip(fx.dofs).filter_map(|(a, d)| d.then_some(a)).collect() };
    let mut said = vec![if fx.normal { "along the surface's own normal (frictionless)".to_string() } else { dofs.join(" ") }];
    let moved: Vec<Option<f64>> = fx.value.iter().map(|&v| num(v)).collect();
    if moved.iter().any(|v| *v != Some(0.0)) {
        let mm = moved.iter().map(|v| v.map_or("a field".to_string(), words::g)).collect::<Vec<_>>().join(", ");
        said.push(format!("moved by ({mm}) mm"));
    }
    let turns: Vec<_> = fx.free.iter().filter(|m| m.omega != [0.0; 3]).collect();
    if turns.len() == 3 && turns.iter().all(|m| m.at == turns[0].at) {
        said.push(format!("free to turn every way about {}", words::point(turns[0].at)));
    } else {
        said.extend(turns.iter().map(|m| format!("free to turn about the axis {} through {}", words::point(m.omega), words::point(m.at))));
    }
    said.extend(fx.free.iter().filter(|m| m.omega == [0.0; 3]).map(|m| format!("free to slide along {}", words::point(m.v))));
    if !fx.free.is_empty() {
        said.push(if fx.rigid { "held rigidly as one but for that".into() } else { "its region's mean motion held, its material free to deform".into() });
    }
    if fx.boundary == Boundary::Support {
        said.push("a physical support".into());
    }
    (dofs, said.join("; "))
}

/// The model's bodies the solve holds, and those it left out and why: by
/// points through each body's box, inside it, against the study's shape and
/// the cells solved. None where the model has no bodies (a study of a field);
/// says why where there are more than `MAX_BODIES`.
fn bodies(s: &mut ModelView, sv: &SolvedStudy) -> Option<J> {
    let items = s.doc.default_items().ok().filter(|v| !v.is_empty())?;
    if items.len() > MAX_BODIES {
        return Some(json!({"unavailable": format!("more than {MAX_BODIES} bodies ({})", items.len())}));
    }
    let (study, ctx) = &sv.shape_field;
    let dropped = sv.sol.report.dropped_cells > 0;
    let (mut solved, mut left) = (Vec::new(), Vec::new());
    for it in &items {
        let Some((lo, hi, _)) = s.item_bounds(it) else { continue };
        let f = s.doc.placed_field(it);
        let tape = crate::tape::Tape::new(&s.doc.g, f);
        let k = BODY_SAMPLES;
        let pts: Vec<[f64; 3]> = (0..k * k * k).map(|i| std::array::from_fn(|a| lo[a] + (hi[a] - lo[a]) * (((i / k.pow(a as u32)) % k) as f64 + 0.5) / k as f64)).collect();
        let (mut body, mut shape) = (Vec::new(), Vec::new());
        crate::eval::eval_points(&tape, &s.ctx, &pts, &mut body);
        let inside: Vec<[f64; 3]> = pts.iter().zip(&body).filter(|(_, v)| **v < 0.0).map(|(p, _)| *p).collect();
        if inside.is_empty() {
            // (thinner than its samples: not judged)
            continue;
        }
        crate::eval::eval_points(study, ctx, &inside, &mut shape);
        let in_shape: Vec<[f64; 3]> = inside.iter().zip(&shape).filter(|(_, v)| **v < 0.0).map(|(p, _)| *p).collect();
        let in_cells = in_shape.iter().filter(|p| sv.sol.value_at(**p).is_some()).count();
        let (all, shape, cells) = (inside.len() as f64, in_shape.len() as f64, in_cells as f64);
        let why = if in_shape.is_empty() {
            Some("outside the study's shape".to_string())
        } else if in_cells == 0 {
            Some(if dropped { "joined to no fix, so dropped as floating" } else { "in no cell of the solve" }.to_string())
        } else if shape < 0.95 * all {
            Some(format!("{:.0}% of it lies outside the study's shape", 100.0 * (1.0 - shape / all)))
        } else if dropped && cells < 0.9 * shape {
            Some(format!("{:.0}% of it is joined to no fix, so dropped as floating", 100.0 * (1.0 - cells / shape)))
        } else {
            None
        };
        match why {
            None => solved.push(json!(it.label)),
            Some(why) => left.push(json!({"body": it.label, "why": why, "solved_share": n(cells / all)})),
        }
    }
    Some(json!({"solved": solved, "left_out": left}))
}

/// `ink solve --plain`: the receipt and the headline results in words.
pub(crate) mod words {
    use serde_json::Value as J;
    use std::fmt::Write;

    /// A number to 4 significant digits.
    pub(crate) fn g(v: f64) -> String {
        if !v.is_finite() {
            return "none".into();
        }
        if v == 0.0 {
            return "0".into();
        }
        let mag = v.abs().log10().floor() as i32;
        if !(-4..=6).contains(&mag) {
            return format!("{v:.3e}");
        }
        let s = format!("{:.*}", (3 - mag).max(0) as usize, v);
        let s = if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s };
        if s == "-0" { "0".into() } else { s }
    }

    /// A point or a vector, its components under a billionth of its largest shown as 0.
    pub(crate) fn point(p: [f64; 3]) -> String {
        let big = p.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let c = p.map(|v| if v.abs() <= 1e-9 * big { 0.0 } else { v });
        format!("({}, {}, {})", g(c[0]), g(c[1]), g(c[2]))
    }

    fn num(v: &J) -> Option<f64> {
        v.as_f64()
    }

    fn vec(v: &J) -> Option<String> {
        let a = v.as_array().filter(|a| a.len() == 3)?;
        Some(point([num(&a[0])?, num(&a[1])?, num(&a[2])?]))
    }

    /// A value: a number, a vector, or a string as it is.
    fn any(v: &J) -> String {
        vec(v).or_else(|| num(v).map(g)).or_else(|| v.as_str().map(str::to_string)).unwrap_or_else(|| "none".into())
    }

    /// Where a footprint lies: its measures, its box, its centre.
    fn place(l: &J, centre: &str) -> String {
        if let Some(p) = l["plane"].as_str() {
            return format!("on the plane {p}");
        }
        if let Some(why) = l["unavailable"].as_str() {
            return format!("nowhere: {why}");
        }
        let mut said = Vec::new();
        if let Some(a) = num(&l["area_mm2"]) {
            said.push(format!("{} mm2 of surface", g(a)));
        }
        if let Some(v) = num(&l["volume_mm3"]) {
            said.push(format!("{} mm3 of material", g(v)));
        }
        let range = |a: usize| format!("{} {}..{}", ["x", "y", "z"][a], g(num(&l["lo"][a]).unwrap_or(f64::NAN)), g(num(&l["hi"][a]).unwrap_or(f64::NAN)));
        let mut s = format!("on {}, {} {} {}, {centre} {} mm", said.join(" and "), range(0), range(1), range(2), any(&l["at"]));
        if let Some(p) = l.get("pressure_mpa") {
            let _ = write!(s, "; pressure {}..{} MPa, mean {}", any(&p["min"]), any(&p["max"]), any(&p["mean"]));
        }
        s
    }

    /// The receipt (`receipt::add`'s, in `out`), and what the solve found, in words.
    pub(crate) fn of(out: &J) -> String {
        let mut s = String::new();
        let heat = out["material"].get("k_w_mk").is_some();
        let _ = writeln!(s, "study {} on cells of {} mm ({})", any(&out["study"]), any(&out["res"]), any(&out["estimate"]));
        let m = &out["material"];
        let from = |k: &str| m["from"][k].as_str().map_or(String::new(), |f| format!(" from {f}"));
        let props: &[(&str, &str, &str)] = if heat {
            &[("k", "k_w_mk", " W/(m K)"), ("density", "density_kg_m3", " kg/m3"), ("capacity", "capacity_j_kgk", " J/(kg K)")]
        } else {
            &[("E", "E_mpa", " MPa"), ("nu", "nu", ""), ("yield", "yield_mpa", " MPa"), ("density", "density_kg_m3", " kg/m3")]
        };
        let said: Vec<String> = props.iter().filter(|(_, key, _)| m.get(*key).is_some()).map(|(k, key, unit)| format!("{k} {}{unit}{}", any(&m[*key]), from(k))).collect();
        let _ = writeln!(s, "material: {}", said.join("; "));
        if let Some(b) = out.get("bodies") {
            if let Some(why) = b["unavailable"].as_str() {
                let _ = writeln!(s, "bodies: not checked ({why})");
            } else {
                let names = |v: &J| v.as_array().map_or(String::new(), |a| a.iter().map(any).collect::<Vec<_>>().join(", "));
                let left: Vec<String> = b["left_out"].as_array().into_iter().flatten().map(|l| format!("{} ({})", any(&l["body"]), any(&l["why"]))).collect();
                let _ = writeln!(s, "bodies: solved {}{}", names(&b["solved"]), if left.is_empty() { String::new() } else { format!("; left out {}", left.join(", ")) });
            }
        }
        let unit = if heat { "W" } else { "N" };
        let centre = out.get("part_centroid_mm").map(|c| format!(" (moments in N mm; about the part's centroid {} mm)", any(c))).unwrap_or_default();
        let _ = writeln!(s, "fixes{centre}:");
        for rx in out["reactions"].as_array().into_iter().flatten() {
            let _ = writeln!(s, "  {}: holds {}", any(&rx["name"]), any(&rx["holds"]));
            let _ = writeln!(s, "    {}", place(&rx["landed"], "centred at"));
            let share = num(&rx["share"]).map_or(String::new(), |v| format!(", {}% of all", g(100.0 * v)));
            let moment = rx.get("moment_about_part_centroid_nmm").map_or(String::new(), |m| format!("; moment {}", any(m)));
            let _ = writeln!(s, "    takes {} {unit}{share}{moment}", any(&rx["total"]));
        }
        let _ = writeln!(s, "loads:");
        for l in out["loads"].as_array().into_iter().flatten() {
            let _ = writeln!(s, "  {}: {} {} {unit}", any(&l["name"]), any(&l["kind"]), any(&l["total"]));
            if let Some(f) = l.get("landed") {
                let _ = writeln!(s, "    {}", place(f, "acting at"));
            }
            if let Some(m) = l.get("moment_nmm") {
                let _ = writeln!(s, "    moment {} about the origin, {} about the part's centroid", any(m), any(&l["moment_about_part_centroid_nmm"]));
            }
        }
        let r = &out["results"];
        let read = |key: &str, what: &str, unit: &str| r.get(key).and_then(|v| num(&v["value"]).map(|x| format!("{what} {}{unit} at {}", g(x), any(&v["at"]))));
        let mut found: Vec<String> = [
            read("max_displacement_mm", "max displacement", " mm"),
            read("max_von_mises_mpa", "max von Mises", " MPa"),
            read("min_safety_factor", "least safety factor", ""),
            read("max_temperature", "max temperature", " C"),
        ]
        .into_iter()
        .flatten()
        .collect();
        found.extend(r["modes"].as_array().and_then(|m| m.first()).map(|m| format!("lowest mode {} Hz", any(&m["frequency_hz"]))));
        found.extend(r.get("unavailable").map(|w| format!("none: {}", any(w))));
        let _ = writeln!(s, "results: {}", found.join("; "));
        let _ = writeln!(s, "equilibrium imbalance {}", any(&out["equilibrium_imbalance"]));
        for w in out["warnings"].as_array().into_iter().flatten() {
            let _ = writeln!(s, "warning: {}", any(w));
        }
        s
    }
}
