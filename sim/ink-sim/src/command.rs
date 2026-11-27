//! `ink solve`: run a (study ...) and answer its results, probes,
//! and render.

use crate::{SolvedStudy, Study, StudyPhysics};
use ink_core::tape::Tape;
use ink_core::v3 as vec3;
use ink_io::{export, render};
use ink_kernel::mesh;
use ink_kernel::query::max_extent;
use ink_lang::InkError;
use ink_lang::RunOpts;
use ink_lang::command::{Args, Command, Opt, Records, form_named, inexact_reads, measured, n, parse_points, solved, v3};
use ink_lang::view::ModelView;
use serde_json::{Value as J, json};
use std::sync::Arc;

/// The package's commands.
pub(crate) static COMMANDS: [Command; 1] = [Command {
    name: "solve",
    options: &[Opt::reads("study"), Opt::reads("res"), Opt::reads("points"), Opt::writes("render"), Opt::reads("size"), Opt::reads("keys")],
    help: "Run a (study ...): finite elements on the implicit shape, no mesh. Elastic: displacement, stress, safety
                factor and buckling factor under fixes and loads, or the lowest natural frequencies ((modes n)); heat:
                temperatures and flux, steady or over time ((transient ...)). Every number is an estimate; ink require
                decides a claim on a study by three grids. ink spec --section Studies: the form.
                [--study NAME] (needed where the file has more than one) [--res MM] (the cell size; default the study's)
                [--points '[[x,y,z]]'] (probe displacement and stress, or temperature, at points)
                [--render FILE.png] (the part coloured by stress or temperature, deformed) [--size PX] [--keys K]
                example: ink solve examples/beam.ink",
    tool: r#"Runs a (study ...): FEA, elastic (with buckling, or modes) or heat (steady or transient). {"study": "load", "res": 1, "points": [[0, 0, 10]], "render": "stress.png", "size": 1024, "keys": "von_mises_mpa"}"#,
    run: ink_lang::command::Run::Model(cmd_solve),
}];

fn run_study(s: &mut ModelView, name: &str, res: Option<f64>) -> Result<SolvedStudy, InkError> {
    // What it reads of the model, at the motion's worst moment.
    let ctx = match s.at_worst(Some("sim"), name)? {
        Some((ctx, _)) => ctx,
        None => s.ctx.clone(),
    };
    SolvedStudy::run(&mut s.doc, &ctx, name, &RunOpts { res, ..RunOpts::default() })
}

/// A solve's answer; `margin`, the support margin it was solved with
/// (cells this close to a fix are reported apart).
fn study_summary(sv: &SolvedStudy, margin: f64) -> J {
    let r = &sv.sol.report;
    // A force's total is a vector (N); a heat load's a power (W), one
    // number, which the solver keeps in the first component.
    let heat = matches!(sv.physics, StudyPhysics::Heat { .. });
    let total = |v: [f64; 3]| if heat { n(v[0]) } else { v3(v) };
    let applied: Vec<J> = r
        .applied
        .iter()
        .enumerate()
        .map(|(i, (k, v))| {
            let mut load = json!({"name": k, "total": total(*v)});
            if !heat && let Some(m) = r.applied_moments.get(i) {
                load["moment_nmm"] = v3(*m);
            }
            if let Some(Some(f)) = r.pressure_error.get(i) {
                load["pressure_error_fraction"] = if f.is_finite() { n(*f) } else { json!("unbounded") };
            }
            load
        })
        .collect();
    let reactions: Vec<J> = r
        .reactions
        .iter()
        .enumerate()
        .map(|(i, (k, v))| {
            let mut rx = json!({"name": k, "total": total(*v)});
            if !heat && let Some(m) = r.reaction_moments.get(i) {
                rx["moment_nmm"] = v3(*m);
            }
            rx
        })
        .collect();
    let mut sum = [0.0; 3];
    let mut scale = 0.0f64;
    for (_, v) in r.applied.iter().chain(&r.reactions) {
        for d in 0..3 {
            sum[d] += v[d];
            scale = scale.max(v[d].abs());
        }
    }
    let (lo, hi) = sv.sol.bounds;
    let span = vec3::norm(vec3::sub(hi, lo));
    let (moment_imbalance, moment_scale) = if !heat && (!r.applied_moments.is_empty() || !r.reaction_moments.is_empty()) {
        // The report keeps moments about the model's origin. Measure balance
        // about the part's center so the residual and its scale move with it.
        let center = vec3::scale(vec3::add(lo, hi), 0.5);
        let mut msum = [0.0; 3];
        let mut mscale = 0.0f64;
        for ((_, force), moment) in r.applied.iter().zip(&r.applied_moments).chain(r.reactions.iter().zip(&r.reaction_moments)) {
            let m = vec3::sub(*moment, vec3::cross(center, *force));
            for d in 0..3 {
                msum[d] += m[d];
                mscale = mscale.max(m[d].abs());
            }
        }
        (Some(vec3::norm(msum) / mscale.max(scale * span).max(1e-300)), mscale)
    } else {
        (None, 0.0)
    };
    let fscale = if !heat && span > 0.0 { scale.max(moment_scale / span) } else { scale };
    let force_imbalance = vec3::norm(sum) / fscale.max(1e-300);
    let imbalance = force_imbalance.max(moment_imbalance.unwrap_or(0.0));
    let unresolved = r.unresolved_volume / r.volume.max(1e-300);
    // A peak, or why there is none (a failed solve; no cell to measure).
    let peak = |e: Result<crate::Extreme, String>, extra: J| match e {
        Ok(e) => {
            let mut o = json!({"value": n(e.value), "at": v3(e.at)});
            if let (J::Object(o), J::Object(x)) = (&mut o, extra) {
                o.extend(x);
            }
            o
        }
        Err(why) => json!({"value": null, "unavailable": why}),
    };
    let pct = match sv.sol.cell_scalar_percentiles(&[0.5, 0.9, 0.99]) {
        Ok(p) => json!({"p50": n(p[0]), "p90": n(p[1]), "p99": n(p[2])}),
        Err(why) => json!({"unavailable": why}),
    };
    let results = match sv.physics {
        _ if r.failure.is_some() => json!({"unavailable": r.failure}),
        StudyPhysics::Elastic { .. } if !r.modes.is_empty() => {
            let modes: Vec<J> = r
                .modes
                .iter()
                .enumerate()
                .map(|(i, (f, m))| {
                    json!({"mode": i + 1, "frequency_hz": n(*f),
                "effective_mass_fraction": {"x": n(m[0]), "y": n(m[1]), "z": n(m[2]), "rx": n(m[3]), "ry": n(m[4]), "rz": n(m[5])}})
                })
                .collect();
            json!({"modes": modes, "mass_g": n(r.mass_kg * 1e3)})
        }
        StudyPhysics::Elastic { .. } => {
            let m = sv.sol.max_cell_scalar();
            // The peak's nearest support margin: how far, and whose; at its
            // edge, a note saying what that means.
            // What it leaves out: the margin of each fix that holds
            // rigidly; a plane of symmetry has none.
            let (planes, held): (Vec<_>, Vec<_>) = sv.fix_names.iter().zip(sv.sol.plane_holds.iter().chain(std::iter::repeat(&false))).partition(|(_, p)| **p);
            let names = |v: &[(&String, &bool)]| v.iter().map(|(n, _)| format!("'{n}'")).collect::<Vec<_>>().join(", ");
            let mut excludes = if held.is_empty() {
                "nothing: no fix has a stress margin".to_string()
            } else {
                format!("cells within {} mm of fix {} (the support margin, max_von_mises_near_supports_mpa)", (margin * 1e4).round() / 1e4, names(&held))
            };
            if !planes.is_empty() {
                let (one, many) = (planes.len() == 1, names(&planes));
                excludes += &format!(
                    "; {} {many} {} the normal displacement across the measured face: no support margin",
                    if one { "fix" } else { "fixes" },
                    if one { "holds" } else { "hold" }
                );
            }
            let mut stress_peak = json!({"excludes": excludes});
            // The surface a peak lies on, where it does: its tightest concave radius
            // in its cells, which a peak's read needs (RESOLVED_RADIUS_CELLS).
            let surface = |at: [f64; 3]| {
                sv.peak_surface(at).map(|(r, h)| json!({"radius_mm": n(r), "cell_mm": n(h), "cells_per_radius": n(r / h), "resolved": r / h >= crate::RESOLVED_RADIUS_CELLS}))
            };
            if let Ok(p) = &m
                && let Some(o) = surface(p.at)
            {
                stress_peak["surface"] = o;
            }
            if let Ok(p) = &m
                && let Some(s) = sv.sol.nearest_support(p.at)
            {
                stress_peak["nearest_excluded"] = json!({"fix": sv.fix_names.get(s.fix), "distance_mm": n(s.distance)});
                if let Some(w) = sv.beside_support_warning(p) {
                    stress_peak["note"] = json!(w);
                }
            }
            let mut res = json!({
                "max_displacement_mm": peak(sv.sol.max_nodal(), json!({})),
                "max_von_mises_mpa": peak(m.clone(), stress_peak),
                "max_von_mises_near_supports_mpa": sv.sol.max_near_supports().map(|e| json!({"value": n(e.value), "at": v3(e.at), "within_mm": n(margin),
                    "note": "the peak within the support margin: a fix holds its region rigidly, so stress at its edge is singular (it grows as res shrinks): it bounds nothing; model the support's contact to read it"})),
                "von_mises_percentiles_mpa": pct,
                "compliance_nmm": if sv.moved { json!({"value": null, "unavailable": crate::MOVED_COMPLIANCE}) } else { n(r.energy) },
            });
            if let Some(w) = sv.max_rotation() {
                res["max_rotation_deg"] = json!({"value": n(w.value.to_degrees()), "at": v3(w.at)});
            }
            // Its linear buckling (every static elastic solve that
            // converged): the loads times the factor buckle it.
            if let Some(b) = &r.buckle {
                res["buckling"] = match b.answer() {
                    Some((low, f)) if f.is_infinite() && low.is_finite() => json!({"factor": null, "low": n(low),
                        "note": if b.iterations == 0 { "no buckling mode: no cell's stress is compressive beyond rounding, which excludes factors below `low`, not above" } else { "no buckling mode found: its search found no direction the loads compress, and excluded factors below `low`, not above" },
                        "iterations": b.iterations, "solve_ms": n(b.ms)}),
                    Some((_, f)) if f.is_infinite() => json!({"factor": "none",
                        "note": if b.iterations == 0 { "no buckling mode: no cell's stress is compressive, so nothing can buckle" } else { "no buckling mode: its search excluded every factor" },
                        "iterations": b.iterations, "solve_ms": n(b.ms)}),
                    Some((low, f)) => json!({"factor": n(f), "low": n(low), "at": v3(b.at), "mode": "(study-field study displacement (mode 1))",
                        "iterations": b.iterations, "solve_ms": n(b.ms)}),
                    None => json!({"factor": null, "unavailable": b.failure, "iterations": b.iterations}),
                };
            }
            if let Some(half) = &r.modeled_buckle {
                if let Some(b) = &r.buckle
                    && b.answer().is_some_and(|(_, f)| f.is_finite())
                {
                    res["buckling"]["antisymmetric_axes"] = json!(["x", "y", "z"].into_iter().zip(b.antisymmetric).filter_map(|(a, odd)| odd.then_some(a)).collect::<Vec<_>>());
                }
                res["buckling"]["sectors_solved"] = json!(r.parity_sectors);
                res["buckling"]["symmetric_factor"] = half.answer().map_or(J::Null, |(_, f)| n(f));
            }
            if let StudyPhysics::Elastic { yield_strength: Some(_), .. } = sv.physics {
                let y = sv.yield_strength.map_or(json!("a field"), n);
                res["min_safety_factor"] = match sv.min_safety_factor() {
                    Ok(e) => {
                        let mut o = json!({"value": n(e.value), "at": v3(e.at), "yield_mpa": y});
                        if let Some(s) = surface(e.at) {
                            o["surface"] = s;
                        }
                        // Where the yield falls inside a cell below its reads:
                        // the least the factor can be there, as a claim counts it.
                        for (key, near) in [("yield_between_reads", false), ("yield_between_reads_near_supports", true)] {
                            if let Some(b) = sv.safety_factor_between_reads(near).filter(|b| b.value < e.value) {
                                o[key] = json!({"value": n(b.value), "at": v3(b.at),
                                    "note": "the least the factor can be in the cells not of one material: each one's least yield over the stress's bound over it; refine res to close it"});
                            }
                        }
                        o
                    }
                    Err(why) => {
                        let mut o = json!({"value": null, "yield_mpa": y, "unavailable": why});
                        if let Some(b) = sv.safety_factor_between_reads(false) {
                            o["yield_between_reads"] = json!({"value": n(b.value), "at": v3(b.at),
                                "note": "the least the factor can be in the cells not of one material: each one's least yield over the stress's bound over it; refine res to close it"});
                        }
                        o
                    }
                };
            }
            res
        }
        StudyPhysics::Heat { .. } => {
            let tmin = (0..sv.sol.node_idx.len()).map(|i| sv.sol.nodal(i)[0]).fold(f64::INFINITY, f64::min);
            json!({"max_temperature": peak(sv.sol.max_nodal(), json!({})), "min_temperature": n(tmin),
                   "max_heat_flux_w_mm2": peak(sv.sol.max_cell_scalar(), json!({})),
                   "heat_flux_percentiles": pct})
        }
    };
    let mut results = results;
    if let (Some(tp), Some((pt, pv, at)), None) = (r.peak_time, r.peak, &r.failure) {
        // Transient heat: the peak temperature over time, where and when it
        // is; the fields are the final state.
        let stride = r.history.len().div_ceil(60).max(1);
        let hist: Vec<J> = r.history.iter().step_by(stride).map(|&(t, v)| json!([n(t), n(v)])).collect();
        let mut transient = json!({
            "headline": "max temperature over reliable nodes",
            "peak": {"value": n(pv), "time_s": n(pt), "at": v3(at)},
            "fields_at_time_s": n(tp),
            "history": hist,
        });
        if let Some(te) = crate::forms::time_error(&r.time_levels, r.residual) {
            transient["time_error"] = json!({"value": n(te.error), "how": te.how,
                "peaks": r.time_levels.iter().map(|(k, v, _)| json!({"steps": k, "value": n(*v)})).collect::<Vec<J>>(),
                "note": "the peak against the same run on steps twice and four times as long: an estimate of what its steps miss, not a bound; more (steps n) narrows it"});
        }
        results["transient"] = transient;
        results["max_temperature"] = json!({"value": n(pv), "at": v3(at), "time_s": n(pt)});
    }
    let mut warnings = Vec::new();
    // Some material sits in cells under half full at every surface (about
    // 10% of a wall 5 cells thick); a quarter means walls of about 2 cells.
    if unresolved > 0.25 {
        let (what, peak) = match sv.physics {
            StudyPhysics::Heat { .. } => ("conductance", "local temperature gradients and fluxes"),
            StudyPhysics::Elastic { .. } => ("stiffness", "peak stresses"),
        };
        warnings.push(format!(
            "{:.0}% of the material is in cells under half full at {} mm cells: walls about 2 cells thick or thinner, or faces just past the grid's planes. \
             Their {what} is averaged in, and {peak} there are read from the fuller cell under them, out to the surface; \
             the percentiles leave them out. If a wall is under 2 cells thick, refine --res",
            unresolved * 100.0,
            measured(sv.sol.h)
        ));
    }
    warnings.extend(sv.kinematics_warning());
    warnings.extend(r.warnings.iter().cloned());
    warnings.extend(r.incomplete.iter().map(|e| e.reason(r)));
    if let Some(why) = &r.failure {
        warnings.push(format!("the solve failed, so it has no results: {why}"));
    }
    warnings.extend(sv.dropped_material());
    warnings.extend(sv.hidden_modes_warning());
    let free: Vec<J> = r
        .free_motions
        .iter()
        .map(|(k, ms)| {
            // A slide's amount is mm; a turn's, radians.
            let motions = ms.iter().map(|(m, q)| json!({"motion": m, "amount": n(*q)}));
            json!({"name": k, "motions": motions.collect::<Vec<_>>()})
        })
        .collect();
    let mut out = json!({
        "res": n(sv.sol.h),
        "warnings": warnings,
        "grid": {"cells_full": r.cells_full, "cells_cut": r.cells_cut, "cells_refined": r.cells_refined, "memory_estimate_mb": n(r.memory_estimate / 1048576.0), "refined_res": n(sv.sol.h / f64::from(1u32 << r.refine_depth)), "nodes": r.nodes, "dofs": r.dofs, "clamped_volume_mm3": n(r.clamped_volume),
                 "material_volume_mm3": n(r.volume), "dropped_floating_cells": r.dropped_cells, "dropped_floating_volume_mm3": n(r.dropped_volume),
                 "dropped_sliver_volume_mm3": n(r.sliver_volume),
                 "unresolved_material_fraction": n(r.unresolved_volume / r.volume.max(1e-300))},
        "solver": {"multigrid_levels": r.levels, "iterations": r.iterations, "relative_residual": n(r.residual), "cached": r.cached,
                   "converged": r.converged, "failure": r.failure, "setup_ms": n(r.setup_ms), "solve_ms": n(r.solve_ms)},
        "results": results,
        "loads": applied,
        "reactions": reactions,
        // A load's total is a force (N, a vector), or a heat study's power
        // (W, a number); for elastic studies, `equilibrium_imbalance` covers
        // both force and moment balance (`moment_imbalance` reports moments).
        "equilibrium_imbalance": n(imbalance),
    });
    if let Some(m) = moment_imbalance {
        out["moment_imbalance"] = n(m);
    }
    if !free.is_empty() {
        out["free_motions"] = json!(free);
    }
    // Each integral over a region's surface: the surface its facets
    // covered and the same measure by lines through the cells (mm2, sum
    // |n_a| dA), how far the lines' may be off, and what the integral missed.
    if !r.surface_measures.is_empty() {
        out["surface_measures"] = json!(
            r.surface_measures.iter().map(|m| json!({"name": m.name, "kind": m.kind, "covered_mm2": n(m.covered), "lines_mm2": n(m.lines), "uncertainty_mm2": n(m.uncertainty), "shortfall_mm2": n(m.shortfall), "excess_mm2": n(m.excess)})).collect::<Vec<_>>()
        );
    }
    out
}

/// Turbo at t in 0..1, as std's `color-map` draws it (its published
/// degree-5 fit, sRGB).
fn turbo(t: f64) -> [f32; 3] {
    const C: [[f64; 6]; 3] = [
        [59.28637943, -152.94239396, 132.13108234, -42.66032258, 4.61539260, 0.13572138],
        [2.82956604, 4.27729857, -14.18503333, 4.84296658, 2.19418839, 0.09140261],
        [27.34824973, -89.90310912, 110.36276771, -60.58204836, 12.64194608, 0.10667330],
    ];
    let t = t.clamp(0.0, 1.0);
    C.map(|c| c.iter().fold(0.0, |a, k| a * t + k).clamp(0.0, 1.0) as f32)
}

/// Why a probe's row is null: its point lies in no cell the study solved.
const OUTSIDE: &str = "in no cell of the solve: not in the study's material (air, a bore, a gap) or in material dropped from it (floating, a sliver)";

fn cmd_solve(s: &mut ModelView, a: &Args) -> Result<J, InkError> {
    solved(s)?;
    if s.doc.forms::<Study>().next().is_none() {
        return Err(InkError { hint: Some("add (study name (shape f) (elastic ...) (fix ...) (load ...)); see ink spec".into()), ..InkError::msg("document has no studies") });
    }
    let name = form_named::<Study>(s, "study", a, "study")?;
    let inputs = s.doc.form_inputs(s.doc.ext_form(Some("sim"), &name).expect("selected study"));
    s.want(&inputs)?;
    let res = a.cell("res")?;
    if let Some(p) = a.get("render").filter(|p| !p.to_ascii_lowercase().ends_with(".png")) {
        return Err(InkError::msg(format!("--render {p}: the render is a PNG; name the file .png")));
    }
    let sv = run_study(s, &name, res)?;
    let mut out = study_summary(&sv, sv.margin);
    out["study"] = json!(name);
    // A load that reads the model: the moment it was taken at, and what that can miss.
    if let Some((_, notes)) = s.at_worst(Some("sim"), &name)? {
        out["run_with"] = json!(notes.join("; "));
    }
    // The material as solved.
    if let StudyPhysics::Elastic { e, nu, .. } = sv.physics {
        let num = |f| s.doc.g.as_const(f).map_or(json!("a field"), n);
        let mut m = json!({"E_mpa": num(e), "nu": num(nu)});
        if let StudyPhysics::Elastic { yield_strength: Some(y), .. } = sv.physics {
            m["yield_mpa"] = num(y);
        }
        out["material"] = m;
    }
    // A solve is an estimate, as a claim on it is: no guaranteed bound.
    out["accuracy"] = J::Null;
    out["estimate"] = json!(format!("a finite-cell solve on cells of {} mm, not a bound", sv.sol.h));

    if let (Some(_), Err(why)) = (a.get("points").or(a.get("render")), sv.sol.valid()) {
        return Err(InkError::msg(format!("study '{name}' failed, so it has nothing to probe or render: {why}")));
    }
    if a.get("points").is_some() {
        let pts = parse_points(a)?;
        let probes: Vec<J> = pts
            .iter()
            .map(|&p| match (sv.sol.value_at(p), sv.sol.field_at_point(p)) {
                (Some(u), Some((f, scalar, frac))) => match sv.physics {
                    StudyPhysics::Elastic { .. } => json!({"point": v3(p), "displacement_mm": v3([u[0], u[1], u[2]]),
                        "stress_mpa": {"xx": n(f[0]), "yy": n(f[1]), "zz": n(f[2]), "xy": n(f[3]), "yz": n(f[4]), "zx": n(f[5])},
                        "von_mises_mpa": n(scalar), "cell_material_fraction": n(frac), "estimate": probe_estimate(sv.sol.h, frac)}),
                    StudyPhysics::Heat { .. } => json!({"point": v3(p), "temperature": n(u[0]), "heat_flux_w_mm2": v3([f[0], f[1], f[2]]), "cell_material_fraction": n(frac), "estimate": probe_estimate(sv.sol.h, frac)}),
                },
                // Its row null, and why beside it (Records keeps it whatever --keys picks).
                _ => json!({"point": v3(p), "outside": OUTSIDE}),
            })
            .collect();
        out["probes"] = Records::new(a).shape(probes, &["outside"])?;
    }

    if let Some(path) = a.get("render") {
        let st = s.doc.form::<Study>(&name).unwrap();
        let shape = st.shape;
        let (lo, hi, _) = s.bounds(shape).unwrap();
        let h = max_extent(lo, hi) / 200.0;
        let f = s.doc.clip(shape, [lo[0] - h, lo[1] - h, lo[2] - h], [hi[0] + h, hi[1] + h, hi[2] + h]);
        let m = mesh::mesh(&Arc::new(Tape::new(&s.doc.g, f)), &s.ctx, lo, hi, h);
        let heat = matches!(sv.physics, StudyPhysics::Heat { .. });
        let scalars: Vec<Option<f64>> = m.verts.iter().map(|&p| if heat { sv.sol.value_near(p).map(|u| u[0]) } else { sv.sol.scalar_near(p) }).collect();
        let (vmin, vmax) = if heat {
            let v: Vec<f64> = scalars.iter().flatten().copied().collect();
            (v.iter().copied().fold(f64::INFINITY, f64::min), v.iter().copied().fold(f64::NEG_INFINITY, f64::max))
        } else {
            (0.0, sv.sol.max_cell_scalar().map_or(0.0, |e| e.value))
        };
        let colors: Vec<[f32; 3]> = scalars.iter().map(|v| v.map(|v| turbo((v - vmin) / (vmax - vmin).max(1e-300))).unwrap_or([0.6, 0.6, 0.6])).collect();
        // The displacement drawn as 5% of the part's size.
        let d = if heat { 0.0 } else { sv.sol.max_nodal().map_or(0.0, |e| e.value) };
        let deform = if d > 0.0 { 0.05 * max_extent(lo, hi) / d } else { 0.0 };
        // Every vertex moved by the solution near it, continued past the
        // cells (value_near): one left where it is while its neighbours move
        // folds the triangles between them.
        let mut shown = m.clone();
        if deform > 0.0 {
            for v in shown.verts.iter_mut() {
                if let Some(u) = sv.sol.value_near(*v) {
                    for k in 0..3 {
                        v[k] += deform * u[k];
                    }
                }
            }
        }
        let size = a.f64("size")?.unwrap_or(1024.0) as usize;
        let (img, w, _, _) = render::render(&[render::Object::plain(&shown, &colors)], &render::Camera::sheet(), size.clamp(64, 4096), None, &[]);
        let (label, note) = match heat {
            true => ("temperature (°C)", String::new()),
            false => ("von Mises stress (MPa)", format!("displacement drawn x{deform:.0}")),
        };
        let (img, hh) = render::with_scale(img, w, label, (vmin, vmax), &|t| turbo(t), &note);
        export::write(path, &export::png(w, hh, &img)).map_err(|e| InkError::msg(format!("cannot write {path}: {e}")))?;
        out["render"] = json!({"out": path, "colormap": "blue (low) -> cyan -> green -> yellow -> red (high)",
            "quantity": if heat { "temperature" } else { "von Mises stress" }, "range": json!([n(vmin), n(vmax)]),
            "deformation_scale": n(deform), "layout": "2x2: iso, top, front, right; under them the colour scale, with its range and units"});
    }
    if let Some(e) = s.doc.ext.iter().find(|e| e.name == name) {
        let roots = s.doc.form_inputs(e);
        inexact_reads(s, &mut out, &roots);
    }
    Ok(out)
}

/// What a probe's numbers are: one read on one grid, its grid's error not measured (a claim's ladder
/// measures one), not a bound; in a cell part material, a stress (a flux) that scatters about the true one
/// by a share finer cells do not shrink (a Lamé ring's bore read 26% low on 0.2 mm cells, 2% high on 0.25).
fn probe_estimate(h: f64, frac: f64) -> String {
    let one = format!("one read on cells of {h} mm, its grid error not measured: not a bound (a claim on it measures one)");
    if frac < 1.0 {
        format!("{one}; in a cell part material its stress or flux scatters about the true one by a share finer cells do not shrink (a ring's bore read 26% low)")
    } else {
        one
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `ink solve --render`'s colour map is std's `color-map` (a viewer
    /// can draw a result's legend from that): one map, the same at every
    /// value, its coefficients written in both.
    #[test]
    fn the_render_s_colour_map_is_std_s() {
        static REGISTRY: ink_lang::Registry = ink_lang::Registry::new(&[&crate::Sim]);
        let doc = ink_lang::Doc::load_str_in(&REGISTRY, "(import \"std\")\n(def c (color-map x 0 1))", "turbo", std::path::Path::new("."), &[]).unwrap();
        let ink_lang::Value::Vec(c) = doc.def("c").unwrap().clone() else { panic!("a colour") };
        let ctx = doc.ctx(&[], 0.0).unwrap();
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            let std = c.iter().map(|v| match v {
                ink_lang::Value::Field(f, _) => ink_core::eval::eval_point(&Tape::new(&doc.g, *f), &ctx, [t, 0.0, 0.0]).clamp(0.0, 1.0),
                v => panic!("{v:?}"),
            });
            for (k, (a, b)) in std.zip(turbo(t)).enumerate() {
                assert!((a - b as f64).abs() < 1e-6, "channel {k} at {t}: std {a}, the render {b}");
            }
        }
    }
}
