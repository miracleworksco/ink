//! Mass, volume, and inertia measurements.

use crate::{Args, CliView, InkError, J, Records, Session, Value, json, n, query, v3};

pub(crate) fn cmd_mass(s: &mut Session, a: &Args) -> Result<J, InkError> {
    use ink::interval::{I, mul_o};
    let items = s.items(a)?;
    let fixed = a.get("density").and_then(ink::lang::materials::density_option).transpose().map_err(|e| InkError::msg(format!("--density: {e}")))?;
    let asked = a.get("density").filter(|_| fixed.is_none());
    let res = a.cell("res")?;
    let shape = Records::new(a);
    // Its inertia measured as closely as its mass only where --keys asks for it (TENSOR).
    let tensor = a.get("keys").is_some_and(|k| k.split(',').any(|x| TENSOR.contains(&x.trim())));
    let mut out = Vec::new();
    let mut total: Option<ink::rigid::Moments> = None;
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
        if let (Some(key), Some((field, _))) = (asked, &density) {
            s.doc.check_density_attribute(key, &who, &it.attrs, it.shape, *field, &s.ctx)?;
        }
        let label = match (fixed, &density) {
            (Some(v), _) => json!(v),
            (None, Some((_, what))) => json!(what),
            (None, None) => J::Null,
        };
        // Instanced: moments once per body in its own frame, then moved.
        let local = s.local(it);
        let found = if let Some(l) = local {
            let rho = density_of(&it.attrs, &mut s.doc)?.map(|d| d.0);
            let body = if tensor { s.inertia(l.shape, rho, res) } else { s.moments(l.shape, rho, res) };
            let none = query::Measured {
                sums: vec![0.0; 11],
                about: [0.0; 3],
                lo: vec![0.0; 11],
                hi: vec![0.0; 11],
                res: 0.0,
                boxes: 0,
                settled: true,
                why: None,
                sampled: None,
                moved: vec![0.0; 11],
            };
            body.map_or(none, |mut b| {
                let m = l.at.moments(&b.moments());
                b.sums[..10].copy_from_slice(&m.m);
                b.about = m.about;
                b
            })
        } else {
            let (attrs, f) = (s.doc.placed_attrs(it), s.doc.placed_field(it));
            let rho = density_of(&attrs, &mut s.doc)?.map(|d| d.0);
            match if tensor { s.inertia(f, rho, res) } else { s.moments(f, rho, res) } {
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
        // The volume, an estimate and its bounds, measured with the mass
        // (as integrals are: query::measured); of a shape not proven
        // bounded, none.
        let whole = s.item_whole(it);
        let reads = s.doc.placed_reads(std::slice::from_ref(it));
        let exact_shape = !s.doc.inexact(&reads, &[]).any();
        let vol = I::new(found.lo[10], found.hi[10]);
        let volume = match (whole, exact_shape) {
            (true, true) => json!({"estimate": n(found.volume()), "bounds": [n(vol.lo), n(vol.hi)]}),
            (true, false) => json!({"estimate": n(found.volume())}),
            (false, _) => json!({"estimate": null}),
        };
        let mut r = match density {
            Some((rho, _)) => {
                let own = found.moments();
                total = Some(total.map_or(own, |t| t.plus(&own)));
                let mut m = mass_json(Some(&own), &shape);
                // Of a shape not proven bounded, none: they would bound
                // only what lies in the window.
                if whole && exact_shape {
                    let bounds = if !s.doc.inexact(&[rho], &[]).any() { Some(I::new(found.lo[0], found.hi[0])) } else { s.doc.provided_bounds(rho, &s.ctx).map(|d| mul_o(d, vol)) };
                    if let Some(b) = bounds {
                        m["mass_bounds"] = json!([n(b.lo), n(b.hi)]);
                    }
                }
                m
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
        true => mass_json(Some(&total.unwrap_or(ink::rigid::Moments { about: [0.0; 3], m: [0.0; 10] })), &shape),
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
fn mass_json(m: Option<&ink::rigid::Moments>, shape: &Records) -> J {
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

/// A measured answer's smallest box side (`res`), `sampled` where its
/// estimate was taken on cells (its bounds as far as the proof reached),
/// and `unsettled` where it stopped short.
fn measured_json(r: &mut J, m: &query::Measured) {
    r["res"] = n(m.res);
    if m.sampled.is_some() {
        r["sampled"] = json!(true);
    }
    if let Some(why) = m.unsettled() {
        r["unsettled"] = json!(why);
    }
}
