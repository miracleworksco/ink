//! `ink optimize`, over the design loop (driver.rs): it reads the
//! options and prints.

use crate::{Optimize, Optimized, Trace};
use ink_lang::InkError;
use ink_lang::claims::progress;
use ink_lang::command::{Args, Command, Opt, extend, form_named, n, param_args, param_name, solved};
use ink_lang::view::ModelView;
use serde_json::{Value as J, json};

/// The package's commands.
pub(crate) static COMMANDS: [Command; 1] = [Command {
    name: "optimize",
    options: &[Opt::reads("name"), Opt::reads("budget")],
    help: "Run an (optimize ...): vary params to minimize or maximize a metric while required claims hold.
                [--name NAME] [--budget N]",
    tool: r#"Runs an (optimize ...). {"name": "lightest", "budget": 40}"#,
    run: ink_lang::command::Run::Model(cmd_optimize),
}];

fn cmd_optimize(s: &mut ModelView, a: &Args) -> Result<J, InkError> {
    solved(s)?;
    let budget = a.f64("budget")?.map(|b| b as usize);
    let name = form_named::<Optimize>(s, "optimize", a, "name")?;
    // The answer's fields besides its designs, reported first: a run
    // stopped before its answer (ink mcp) gives them with the designs done.
    let context = json!({});
    progress(json!({"answer": context}));
    // Its design is verified as ink assert asks claims (verify), which
    // refuses a value for a range of the question: so does it, before the
    // first design, so that `apply` is a question assert takes.
    let given = param_args(a)?;
    // `apply` gives back the run's --param values, each a name a file
    // writes (command::param_name), so that it reads as it was given.
    if let Some((k, _)) = given.values.iter().find(|(k, _)| !param_name(k)) {
        return Err(InkError::msg(format!(
            "--param {k:?}: optimize sets a param a file names (letters, digits, '-' and '_', starting with a letter, its namespaces joined by '.')"
        )));
    }
    if let Some((k, (lo, hi))) = given.values.iter().find_map(|(k, _)| s.doc.ranged_named(k).map(|r| (k, r.1))) {
        return Err(InkError::msg(format!("--param {k}: '{k}' is (ranged {lo} {hi}), a range of the question: optimize verifies its design over every value of it")));
    }
    let o = crate::run_optimize(s, &name, budget, a)?;
    let names: Vec<String> = o.vary.iter().map(|v| v.0.clone()).collect();
    let mut out = optimized_json(&o, &o.answer_params, &crate::driver::Replay::of(a)?.apply(&names, &o.design));
    if !o.failed.is_empty() {
        out["failed_designs"] = crate::driver::failed_json(&names, &o.failed);
    }
    // Gradient (the method of moving asymptotes) where every derivative the
    // run needs is exact at the start design, else derivative-free (bounded
    // Nelder-Mead with a constraint penalty).
    out["method"] = json!(match o.trace {
        Trace::Simplex { .. } => "derivative-free",
        Trace::Gradient { .. } => "gradient",
    });
    out["warnings"] = json!(o.warnings);
    // Over a motion the search's numbers are of one moment: the answers
    // over the span decide.
    if let Some(t) = o.search_t {
        out["note"] = json!(format!(
            "the search steers by each requirement measured at one moment (t {t}): search_measured, search_violation and search_feasible; each requirement's `answer`, over the span, is the verdict (satisfied, feasible)"
        ));
    }
    // Its objective and requires are measured (estimates, as a claim on a
    // study is): ink assert on `apply` checks the design.
    out["accuracy"] = J::Null;
    extend(&mut out, context);
    Ok(out)
}

/// The answer: each requirement with the search's measure of it
/// (`search_measured`, `search_violation`) and, at the design returned, its claim's answer
/// as `ink assert` gives it (`answer`), which alone decides `satisfied`.
///
/// A design's values are given exactly (the shortest decimal that reads
/// back as the same float), never rounded: they are the design its claims
/// were asked of, and a field may jump at a value's last digit
/// ((floor p) at 0.99999999996).
fn optimized_json(o: &Optimized, names: &[String], apply: &str) -> J {
    let fmt = |x: &[f64]| -> serde_json::Map<String, J> { o.vary.iter().zip(x).map(|(v, x)| (v.0.clone(), json!(x))).collect() };
    let requires: Vec<J> = o
        .requires
        .iter()
        .map(|r| {
            let mut j = json!({"claim": r.claim, "search_measured": n(r.measured), "search_violation": n(r.violation.max(0.0)), "satisfied": r.satisfied()});
            if let Some(e) = r.grid_error {
                j["search_grid_error"] = n(e);
            }
            if let Some(c) = &r.answer {
                j["answer"] = ink_lang::answer::check_json(c, names);
            }
            j
        })
        .collect();
    let sense = if o.maximize { "maximize" } else { "minimize" };
    let apply = apply.to_string();
    match &o.trace {
        Trace::Simplex { evaluations, improvements } => {
            let improvements: Vec<J> = improvements.iter().map(|(k, x, v)| json!({"evaluation": k, "design": fmt(x), "objective": n(*v)})).collect();
            let mut out = json!({"optimize": o.name, "sense": sense,
                "feasible": o.feasible, "objective": n(o.objective),
                "design": fmt(&o.design), "requires": requires,
                "start": {"design": fmt(&o.start), "objective": n(o.start_objective), "search_feasible": o.start_within},
                "evaluations": evaluations, "improvements": improvements,
                "apply": apply});
            // The design returned could not be evaluated: say why (its
            // study failed, the objective's or a requirement's metric had no
            // value there).
            if let Some(why) = &o.failure {
                out["evaluation_failed"] = json!(why);
            }
            out
        }
        Trace::Gradient { iterations, history, exact, differenced } => {
            let evaluations: usize = history.iter().map(|h| h.evaluations).sum();
            let history: Vec<J> = history
                .iter()
                .map(|h| {
                    json!({"design": fmt(&h.design), "objective": n(h.objective),
                "max_violation": n(h.max_violation), "evaluations": h.evaluations})
                })
                .collect();
            json!({"optimize": o.name, "sense": sense,
                "feasible": o.feasible, "objective": n(o.objective), "design": fmt(&o.design),
                "requires": requires,
                "start": {"design": fmt(&o.start), "objective": n(o.start_objective)},
                "iterations": iterations, "design_evaluations": evaluations,
                "derivatives": {"exact_param_iterations": exact, "differenced_param_iterations": differenced},
                "history": history,
                "apply": apply})
        }
    }
}
