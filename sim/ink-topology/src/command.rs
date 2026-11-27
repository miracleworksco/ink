//! `ink topology`: run a (topology ...) and answer its layout's history.

use crate::TopologyDef;
use ink_lang::InkError;
use ink_lang::command::{Args, Command, Opt, form_named, inexact_reads, measured, n, solved};
use ink_lang::view::ModelView;
use serde_json::{Value as J, json};

/// The package's commands.
pub(crate) static COMMANDS: [Command; 1] = [Command {
    name: "topology",
    options: &[Opt::reads("name"), Opt::reads("res")],
    help: "Run a (topology ...): a stiff layout of a study's shape within a volume budget, found by local search.
                [--name NAME] [--res MM]
                example: ink topology examples/layout.ink",
    tool: r#"Runs a (topology ...): a stiff layout of a static elastic study's shape within a volume budget, found by local search (an estimate), and its history. {"name": "layout", "res": 1}"#,
    run: ink_lang::command::Run::Model(cmd_topology),
}];

fn cmd_topology(s: &mut ModelView, a: &Args) -> Result<J, InkError> {
    solved(s)?;
    let name = form_named::<TopologyDef>(s, "topology", a, "name")?;
    let study = s.doc.form::<TopologyDef>(&name).expect("a named topology").study.clone();
    let l = crate::run(&mut s.doc, &s.ctx.clone(), &name, a.cell("res")?, ink_lang::RunOpts::default().tol)?;
    let (r, rec) = (&l.sol.report, &l.record);
    let mut warnings = r.warnings.clone();
    let mut out = json!({"topology": name, "study": study, "res": n(l.sol.h)});
    if let Some(why) = &r.failure {
        warnings.push(format!("the solve failed, so it has no results: {why}"));
    } else {
        let (first, last) = (rec.history[0], *rec.history.last().expect("a run solved a design"));
        let hist: Vec<J> = rec.history.iter().map(|(i, c, v, ch, b)| json!([i, n(*c), n(*v), n(*ch), n(*b)])).collect();
        out["results"] = json!({
            "iterations": last.0 - 1, "settled": rec.settled, "compliance_nmm": n(last.1), "volume_fraction": n(last.2), "last_change": n(last.3),
            "compliance_start_nmm": n(first.1),
            "history": hist,
            "history_columns": ["design", "compliance", "volume fraction", "max density change", "projection steepness"],
            "use": "(study-field TOPOLOGY solid) is the optimized part as a field (negative inside); intersect it with the domain and verify it with a normal study",
        });
        // Its part's buckling, as a study of the part has it (Part): said where the loads may reach it,
        // or why it is not known.
        let p = &l.part;
        match (p.buckle.as_ref().and_then(|b| b.answer()), &p.failure) {
            (Some((low, factor)), _) => {
                out["results"]["part_buckling_factor"] = json!([n(low), n(factor)]);
                if low <= 1.0 {
                    let at = p.buckle.as_ref().map_or([0.0; 3], |b| b.at);
                    warnings.push(format!("topology: its part {}", ink_sim::buckled_warning(low, factor, at, p.h)));
                }
            }
            (None, why) => warnings.push(format!(
                "topology: its part's buckling is not known ({}): a layout settled for stiffness may buckle at its loads; verify the part with a study of its own",
                why.as_deref().unwrap_or("no analysis")
            )),
        }
        // Unsettled; unless keep fills the budget, whose warning says why
        // no number of iterations would settle it.
        let (held, vf) = rec.fixed;
        if held >= vf - 0.02 {
            warnings.push(format!(
                "topology: {:.0}% of the domain is held solid (by keep, or wholly held by a fix), at or above the volume fraction {}: \
                 nothing is left to optimize; shrink keep or the fixes' regions, or raise the volume fraction",
                held * 100.0,
                measured(vf)
            ));
        } else if !rec.settled {
            warnings.push(format!(
                "topology: the layout had not settled when its iterations ran out (the returned design's last change {}, volume fraction {}): \
                 give it more iterations",
                measured(last.3),
                measured(last.2)
            ));
        }
    }
    // What the answer is: an estimate, by a local search, of the density
    // model (grey cells' stiffness density^3), not of the part extracted.
    out["accuracy"] = J::Null;
    out["estimate"] = json!(format!("a numerical solve on cells of {} mm, not a bound", n(l.sol.h)));
    out["method"] = json!(
        "SIMP minimum compliance (a density per cell, stiffness density^3, filtered and projected), updated by the method of moving asymptotes: \
         a local search from a uniform start, so `settled` says the layout stopped changing, not that it is the stiffest. compliance_nmm is \
         the density model's, not the extracted part's: verify the part, (study-field TOPOLOGY solid), with a study of its own"
    );
    out["warnings"] = json!(warnings);
    out["solver"] = json!({"iterations": r.iterations, "relative_residual": n(r.residual), "converged": r.converged, "cached": r.cached, "failure": r.failure});
    if let Some(e) = s.doc.ext.iter().find(|e| e.name == name) {
        let roots = s.doc.form_inputs(e);
        inexact_reads(s, &mut out, &roots);
    }
    Ok(out)
}
