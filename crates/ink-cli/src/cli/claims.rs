//! `assert` and `interfere`: claims checked by the library (ink::claims,
//! ink::lang::interfere), printed as JSON.

use crate::*;
use ink::claims::{AccuracyBy, AssertReport, AssertResult, ClaimOptions, ResBy};
pub(crate) use ink::lang::answer::{ASK_FINER, Names, check_json, proof_json};

use ink::lang::interfere::{self, Answer, Pair};

/// The document's param names, by index (to name a witness's params).
pub(crate) fn param_names(s: &Session) -> Vec<String> {
    s.doc.params.iter().map(|p| p.name.clone()).collect()
}

/// An assert's row of the answer.
fn assert_json(a: &AssertResult, names: &[String]) -> J {
    let mut r = check_json(&a.check, names);
    r["name"] = json!(a.name);
    r["claim"] = json!(a.claim);
    r["source"] = json!(a.source);
    r["ms"] = n(a.ms);
    r["cached"] = json!(a.cached);
    if !a.terms.is_empty() {
        r["terms"] = J::Object(a.terms.iter().map(|(k, x)| (k.clone(), n(*x))).collect());
    }
    // The cell size asked for (each result's `res` is the one it
    // reached), and what asked: the claim's (res r) or --res.
    if let Some((h, by)) = a.res {
        r["res_given"] = n(h);
        r["res_given_by"] = json!(match by {
            ResBy::Claim => "assert",
            ResBy::Run => "--res",
        });
    }
    if let Some(by) = &a.accuracy_by {
        accuracy_by_json(&mut r, by);
    }
    r
}

/// What gave a claim (or a pair) its accuracy, where not its default: the
/// claim's (accuracy x), --accuracy, or its parts' tolerances, each named.
fn accuracy_by_json(r: &mut J, by: &AccuracyBy) {
    r["accuracy_given_by"] = json!(match by {
        AccuracyBy::Claim => "assert",
        AccuracyBy::Run => "--accuracy",
        AccuracyBy::Tolerance(_) => "tolerance",
    });
    if let AccuracyBy::Tolerance(parts) = by {
        r["tolerances"] = J::Object(parts.iter().map(|(k, x)| (k.clone(), n(*x))).collect());
    }
}

fn report_json(rep: &AssertReport, names: &[String]) -> J {
    let asserts: Vec<J> = rep.asserts.iter().map(|a| assert_json(a, names)).collect();
    let counts: serde_json::Map<String, J> = rep.counts().into_iter().map(|(s, k)| (s.as_str().to_string(), json!(k))).collect();
    let mut out = json!({"all_hold": rep.all_hold(), "all_hold_guaranteed": rep.all_hold_guaranteed(), "counts": counts, "asserts": asserts});
    // No claims checked: nothing holds or fails, guaranteed or not.
    if rep.asserts.is_empty() {
        out["all_hold"] = J::Null;
        out["all_hold_guaranteed"] = J::Null;
        out["no_claims"] = json!("no claim was checked: the file states none (assert forms), or none was asked for");
    }
    // Claims asked at one moment though they read t, each at the moment
    // it recorded: said once, at the top, where all_hold is read.
    if let Some(note) = rep.note() {
        out["note"] = json!(note);
    }
    // Claims left out for being imported: all_hold and all_hold_guaranteed
    // say nothing about them.
    if let Some(k) = rep.imported_not_run {
        out["scope"] = json!("own");
        out["imported_not_run"] = json!(k);
    }
    out
}

pub(crate) fn cmd_assert(s: &mut Session, a: &Args) -> Result<J, InkError> {
    let only = a.get("name");
    let res = a.cell("res")?;
    let own = match a.get("own") {
        None => false,
        Some("true") => true,
        Some(v) => return Err(InkError::msg(format!("--own takes no value, got '{v}' (write the file first: ink assert FILE --own)"))),
    };
    // The answer's fields besides its rows, reported first: a run stopped
    // before its answer (ink mcp) gives them with the rows done.
    let mut context = json!({});
    scope_fields(&mut context, s);
    ink::claims::progress(json!({"answer": context}));
    let names = param_names(s);
    let rep = ink::claims::run_asserts_with(s, only, own, ClaimOptions { res, ..Default::default() }, None, |a| Some(assert_json(a, &names)))?;
    let mut out = report_json(&rep, &names);
    extend(&mut out, context);
    at_cases(s, a, only, own, res, &mut out)?;
    // --keys: of each claim, its name and what is picked.
    if let J::Array(rows) = out["asserts"].take() {
        out["asserts"] = J::Array(Records::new(a).picked(rows, &["name"], &[])?);
    }
    Ok(out)
}

/// Every claim asked again at each combination of the cases its params
/// declare (`(cases v ...)`), each param at its value or one of its cases,
/// but those `--param` sets: a claim holds where it holds at all of them,
/// fails where it fails at one, else is unknown. Each row says how it came
/// out at each (`cases`), and the answer which values were asked
/// (`cases`), and which not, a text case its param does not allow there
/// (`cases_not_allowed`). The params with cases are those the file's load
/// declares, and the combinations their product: a combination's load
/// declaring cases on a param that load does not (a `for` over a param),
/// not declaring one it does, or declaring another number of them, is
/// refused. A claim made only at some
/// combinations is reported too, saying so (`at_values` false).
fn at_cases(s: &Session, a: &Args, only: Option<&str>, own: bool, res: Option<f64>, out: &mut J) -> Result<(), InkError> {
    type Combo = std::collections::BTreeMap<String, usize>;
    let given = ink::lang::command::param_args(a)?;
    let set = |n: &str| given.values.iter().any(|(k, _)| k == n);
    // The params with cases the file's load declares, but those --param
    // sets, each with how many.
    let root: Vec<(&str, usize)> = s.doc.cases.iter().filter(|(n, _)| !set(n)).map(|(n, c)| (n.as_str(), c.len())).collect();
    if root.is_empty() {
        return Ok(());
    }
    // Their product, counted before it is made, its first the params'
    // values (the file's load): the rest, fewest cases first.
    if root.iter().try_fold(1usize, |k, (_, c)| k.checked_mul(c + 1)).is_none_or(|k| k - 1 > ink::claims::MOST_CASES) {
        return Err(InkError::msg(format!("the params' cases make more than the {} combinations one run asks: declare fewer, or set one with --param", ink::claims::MOST_CASES)));
    }
    let mut combos = vec![Combo::new()];
    for (p, k) in &root {
        combos = combos
            .iter()
            .flat_map(|c| {
                std::iter::once(c.clone()).chain((0..*k).map(|i| {
                    let mut c = c.clone();
                    c.insert(p.to_string(), i);
                    c
                }))
            })
            .collect();
    }
    combos.remove(0);
    combos.sort_by_key(Combo::len);
    let shown = |g: &ink::lang::Given| g.text.as_deref().map_or_else(|| n(g.value), |t| json!(t));
    // Each claim asked at each combination is one more: counted with those
    // done, so a run stopped says how far it got through them all; those
    // waiting count the file's claims each.
    let rows = out["asserts"].as_array().map_or(0, Vec::len);
    let mut done = rows;
    let progress_total = |done: usize, waiting: usize| ink::claims::progress(json!({"total": done + waiting * rows, "of": "asserts", "kept": true}));
    progress_total(done, combos.len());
    let (mut asked, mut not_allowed): (Vec<J>, Vec<J>) = (Vec::new(), Vec::new());
    // The combinations not loaded: the call's time limit passed first.
    let mut not_loaded: Vec<J> = Vec::new();
    let by_place = |c: &Combo| J::Object(c.iter().map(|(k, c)| (k.clone(), json!({"case": c + 1}))).collect());
    for (i, combo) in combos.iter().enumerate() {
        if ink::stop::reached() {
            not_loaded.push(by_place(combo));
            continue;
        }
        let waiting = combos.len() - i - 1;
        // One load at it: each param given takes its case where the load
        // reads it (Given::case), the case's expression read there, every
        // param before it set (one may read another's).
        let changes: Vec<(String, ink::lang::Given)> = combo.iter().map(|(k, c)| (k.clone(), ink::lang::Given::case(*c))).collect();
        // Its load (the model built again at it), said as it runs and with
        // its time, as a claim is, by the cases it takes; the values they
        // took once it is loaded.
        let load = json!({"load": J::Object(combo.iter().map(|(k, c)| (k.clone(), json!({"case": c + 1}))).collect())});
        ink::claims::progress(json!({"running": load}));
        let start = std::time::Instant::now();
        let here = |m: &str| format!("at the cases {}: {m}", load["load"]);
        let mut v = match s.reload(&changes) {
            Err(_) if ink::stop::reached() => {
                not_loaded.push(by_place(combo));
                continue;
            }
            r => r.map_err(|e| InkError { message: here(&e.message), ..e })?,
        };
        // Its params with cases are the file's load's (one past a refusal
        // not read), each with as many.
        let refused = !v.doc.refused.is_empty();
        let declares = |p: &str| v.doc.cases.iter().any(|(c, _)| c == p);
        if let Some((p, _)) = v.doc.cases.iter().find(|(p, _)| !set(p) && !root.iter().any(|(r, _)| r == p)) {
            return Err(InkError::msg(here(&ink::claims::cases_not_everywhere(p, "is declared at this combination but not at the params' values"))));
        }
        if let Some((p, first, at_load)) = v.doc.cases.iter().find_map(|(p, c)| root.iter().find(|(r, k)| r == p && *k != c.len()).map(|(_, k)| (p, *k, c.len()))) {
            return Err(InkError::msg(here(&ink::claims::cases_count_changes(p, first, at_load))));
        }
        if let Some((p, _)) = root.iter().find(|(p, _)| !refused && !declares(p)) {
            return Err(InkError::msg(here(&ink::claims::cases_not_everywhere(p, "is declared at the params' values but not at this combination"))));
        }
        // The cases its load took (Given::case resolved; one past a
        // refusal, none).
        let took: std::collections::BTreeMap<&String, ink::lang::Given> =
            v.overrides.iter().filter(|(_, g)| g.case.is_some() && (g.text.is_some() || g.value.is_finite())).map(|(k, g)| (k, g.clone())).collect();
        // What each case took, as the load read it (one it did not reach,
        // past a refusal: by its place).
        let at: serde_json::Map<String, J> = combo.iter().map(|(k, c)| (k.clone(), took.get(k).map_or_else(|| json!({"case": c + 1}), shown))).collect();
        ink::claims::progress(json!({"loaded": load, "values": J::Object(at.clone()), "ms": start.elapsed().as_secs_f64() * 1e3}));
        // A text param that does not allow the text it takes at this
        // combination (its len may be another param's) takes no text
        // anyone can give there: the load stopped there; not asked, and
        // said.
        if refused {
            let why: Vec<String> = v.doc.refused.iter().map(|(p, t, r)| r.said(p, t)).collect();
            not_allowed.push(json!({"case": J::Object(at), "why": why.join("; ")}));
            progress_total(done, waiting);
            continue;
        }
        asked.push(J::Object(at.clone()));
        if let Ok(scope) = ink::lang::command::claim_scope(&v, a) {
            v.set_claim_scope(&scope)?;
        }
        let names = param_names(&v);
        let rep = ink::claims::run_asserts_with(&mut v, only, own, ClaimOptions { res, ..Default::default() }, Some(&J::Object(at.clone())), |_| None)?;
        done += rep.asserts.len();
        progress_total(done, waiting);
        let rows = out["asserts"].as_array_mut().expect("rows");
        for r in &rep.asserts {
            let mut c = assert_json(r, &names);
            // A claim the file's load does not make (a `for` over a
            // param): its row is of its cases alone.
            if !rows.iter().any(|x| x["name"] == r.name.as_str()) {
                let mut row = c.clone();
                for k in ["ms", "cached", "held", "near"] {
                    row.as_object_mut().expect("a row").remove(k);
                }
                row["at_values"] = json!(false);
                rows.push(row);
            }
            let row = rows.iter_mut().find(|x| x["name"] == r.name.as_str()).expect("its row");
            for k in ["name", "claim", "source", "ms", "cached"] {
                c.as_object_mut().expect("a row").remove(k);
            }
            c["case"] = J::Object(at.clone());
            let rank = |st: &J| match st.as_str() {
                Some("fails") => 2,
                Some("unknown") => 1,
                _ => 0,
            };
            if rank(&c["status"]) > rank(&row["status"]) {
                row["status"] = c["status"].clone();
                // What it held by at the values is no answer for the claim.
                for k in ["held", "near"] {
                    row.as_object_mut().expect("a row").remove(k);
                }
            }
            if c["accuracy"].is_null() || row["accuracy"].as_f64().zip(c["accuracy"].as_f64()).is_some_and(|(x, y)| y > x) {
                row["accuracy"] = c["accuracy"].clone();
            }
            row.as_object_mut().expect("a row").entry("cases").or_insert(json!([])).as_array_mut().expect("cases").push(c);
        }
    }
    // A claim is not checked at a combination not loaded.
    if !not_loaded.is_empty() {
        for row in out["asserts"].as_array_mut().expect("rows") {
            if row["status"] == "holds" {
                row["status"] = json!("unknown");
                for k in ["held", "near"] {
                    row.as_object_mut().expect("a row").remove(k);
                }
            }
            row["accuracy"] = J::Null;
            let cases = row.as_object_mut().expect("a row").entry("cases").or_insert(json!([])).as_array_mut().expect("cases");
            cases.extend(not_loaded.iter().map(|c| json!({"case": c, "status": "unknown", "accuracy": null, "hint": NOT_LOADED})));
        }
        out["cases_not_loaded"] = json!(not_loaded);
    }
    // The answer's summary over every case.
    let rows = out["asserts"].as_array().expect("rows").clone();
    let holds = |r: &J| r["status"] == "holds";
    let mut counts = serde_json::Map::new();
    for st in ["holds", "fails", "unknown"] {
        let k = rows.iter().filter(|r| r["status"] == st).count();
        if k > 0 {
            counts.insert(st.into(), json!(k));
        }
    }
    if !rows.is_empty() {
        // A claim made only at a case was checked.
        out.as_object_mut().expect("an answer").remove("no_claims");
        out["all_hold"] = json!(rows.iter().all(holds));
        out["all_hold_guaranteed"] = json!(rows.iter().all(|r| holds(r) && r["accuracy"].is_number()));
        out["counts"] = J::Object(counts);
    }
    out["cases"] = json!(asked);
    if !not_allowed.is_empty() {
        out["cases_not_allowed"] = json!(not_allowed);
    }
    Ok(())
}

/// A combination of cases not loaded: the call's time limit had passed.
const NOT_LOADED: &str = "not checked at this case: the call's time limit (timeout_s) passed before it was loaded. A call with --param set to it checks it";

/// A pair not begun: the call's time limit had passed.
const NOT_BEGUN: &str = "not checked: the call's time limit (timeout_s) had passed. A call on just the assembly holding it (target) checks it";

/// An overlap's volume as a row gives it: an estimate. Of an
/// intersection not proven bounded (`whole` false), the part in the
/// window, and over ranged params (`ranged`) at one value of them: none.
fn volume_json(estimate: f64, whole: bool, ranged: bool) -> J {
    json!({"estimate": if whole && !ranged { n(estimate) } else { J::Null }})
}

/// What a pair's row names besides its answer: each relation under the
/// target as a contact zone gives it (its ends), and its kind as a
/// relation's row gives it (CliView::kind_json), by Plan::relations'
/// index; and whether the scope has ranged params.
struct Rows<'a> {
    names: &'a [String],
    zones: Vec<J>,
    kinds: Vec<Option<J>>,
    ranged: bool,
}

/// A pair's row of interfere's answer (interfere::Pair): its proof or its
/// claim's answer, its time, what it reads (an estimate's accuracy is
/// null), what holding in the window says, its items, contact zones and
/// relations.
fn pair_json(p: &Pair, rows: &Rows) -> J {
    let names = rows.names;
    let mut r = match &p.answer {
        Answer::Still { check, volume } => {
            let mut r = check_json(check, names);
            if let Some(h) = r["hint"].as_str().filter(|_| r.get("undecided_cell").is_some()) {
                r["hint"] = json!(format!(
                    "{h}. If the pair touches there by design, declare it: (contact A B ZONE), or (relate A B KIND) with a zone: a simple shape about undecided_cell, clipped to the contact"
                ));
            }
            if let Some(v) = volume {
                r["overlap_volume"] = volume_json(v.estimate, p.whole, rows.ranged);
            }
            r
        }
        Answer::Moving { check } => check_json(check, names),
        Answer::NotBegun => json!({"status": "unknown", "accuracy": null, "hint": NOT_BEGUN}),
        Answer::Unsolved(why) => json!({"status": "unknown", "accuracy": null, "unavailable": why}),
    };
    if let Some(ms) = p.ms {
        r["ms"] = n(ms);
    }
    if let Some(c) = p.cached {
        r["cached"] = json!(c);
    }
    if let Some(t) = p.one_moment {
        r["one_moment"] = n(t);
    }
    if let Some(by) = &p.accuracy_by {
        accuracy_by_json(&mut r, by);
    }
    // A pair whose intersection is not proven bounded was examined in the
    // window only (a claim's check says so itself).
    if !p.whole && matches!(p.answer, Answer::Unsolved(_)) {
        r["unbounded"] = json!(true);
    }
    // A moving pair's claim says what it read inexactly as decided (a
    // solve on t ranged with t is certified there).
    let mut f = json!({});
    inexact_fields(&mut f, if let Answer::Moving { check } = &p.answer { &check.inexact } else { &p.inexact }, names);
    f["pair"] = json!(p.items);
    if !p.zones.is_empty() {
        f["contact_zones"] = J::Array(p.zones.iter().map(|&z| rows.zones[z].clone()).collect());
    }
    // What the pair's relations say (`relate`): their kinds, data, what
    // their zones leave unchecked, their claims.
    let kinds: Vec<J> = p.relations.iter().filter_map(|&z| rows.kinds[z].clone()).collect();
    if !kinds.is_empty() {
        f["relations"] = J::Array(kinds);
    }
    extend(&mut r, f);
    // A claim's hint offers (accuracy x): interfere's is --accuracy.
    if let Some(h) = r["hint"].as_str().filter(|h| h.contains(ASK_FINER)) {
        r["hint"] = json!(h.replace(ASK_FINER, PAIR_FINER));
    }
    r
}

/// What a pair's hint says instead of a claim's (answer::ASK_FINER).
const PAIR_FINER: &str = "A smaller --accuracy looks finer";

pub(crate) fn cmd_interfere(s: &mut Session, a: &Args) -> Result<J, InkError> {
    // What the target places, resolved once (Doc::target): a definition
    // name and a path placing it are asked alike. A pair that reads t over
    // every span the model states (Doc::spans_hull), whichever file moves
    // what it reads; one reading none is static.
    let target = s.doc.target(a.get("target").unwrap_or("root"));
    s.motion_over(s.doc.spans_hull())?;
    let all = match a.get("target") {
        Some(t) => s.doc.items(t)?,
        None => interfere::default_items(s)?,
    };
    let res = a.cell("res")?;
    let accuracy = match a.f64("accuracy")? {
        Some(x) if !(x > 0.0 && x.is_finite()) => return Err(InkError::msg(format!("--accuracy expects a depth above 0, got '{}'", a.get("accuracy").unwrap_or_default()))),
        // As a claim's (accuracy x) is: a number past ink's range is refused.
        Some(x) if ink::graph::wild(x) => return Err(InkError::msg(format!("--accuracy {x:e} is outside ink's range (0, or 1e-100 to 1e100 in size)"))),
        x => x,
    };
    let plan = interfere::plan(s, all, target.as_ref(), res, accuracy)?;
    let names = param_names(s);
    let rows = Rows {
        names: &names,
        zones: plan.relations.iter().map(|r| json!([r.a, r.b])).collect(),
        kinds: plan.relations.iter().map(|r| r.kind.as_ref().map(|k| s.kind_json(k))).collect(),
        ranged: !s.ctx.p_spans.is_empty(),
    };
    // The answer's fields besides its pairs, reported first: a run stopped
    // before its answer (ink mcp) gives them with the pairs done.
    let mut out = json!({"accuracy": n(plan.accuracy), "pairs_apart_by_bounding_boxes": plan.apart});
    if !plan.unsolved_items.is_empty() {
        out["unsolved_items"] = J::Array(plan.unsolved_items.iter().map(|(it, e)| json!({"item": it, "unavailable": e.message})).collect());
    }
    if !plan.unbounded_items.is_empty() {
        out["unbounded_items"] = json!(plan.unbounded_items);
    }
    if plan.boxed_apart > plan.apart {
        out["pairs_apart_by_estimated_boxes"] = json!(plan.boxed_apart - plan.apart);
    }
    scope_fields(&mut out, s);
    // What it asks at one moment though the question reads t, known
    // before any pair runs: said once, at the top.
    if let Some(note) = ink::claims::one_moment_note(&plan.one_moment) {
        out["note"] = json!(note);
    }
    ink::claims::progress(json!({"answer": out.clone()}));
    let checked = plan.run(s, res, &|p| Some(pair_json(p, &rows)))?;
    let mut pairs = Vec::with_capacity(checked.len());
    let (mut overlaps, mut estimated_overlaps, mut estimated_holds, mut undecided, mut near) = (0, 0, 0, 0, 0);
    for p in &checked {
        let r = pair_json(p, &rows);
        if p.overlaps() {
            overlaps += 1;
            estimated_overlaps += usize::from(r["accuracy"].is_null());
        }
        estimated_holds += usize::from(r["status"] == "holds" && r["accuracy"].is_null());
        undecided += usize::from(r["status"] == "unknown");
        near += usize::from(r["near"]["depth"].as_f64().is_some_and(|d| d > 0.0));
        pairs.push(r);
    }
    out["overlapping_pairs"] = json!(overlaps);
    out["undecided_pairs"] = json!(undecided);
    out["near_pairs"] = json!(near);
    if estimated_overlaps > 0 {
        out["overlapping_pairs_estimated"] = json!(estimated_overlaps);
    }
    // Pairs holding as estimates (measured: sampled steps, inexact inputs),
    // which a clash between their samples would not fail.
    if estimated_holds > 0 {
        out["holding_pairs_estimated"] = json!(estimated_holds);
    }
    // What needs acting on first: pairs that fail, then those undecided,
    // so an answer shortened to its first rows (ink mcp) keeps them.
    pairs.sort_by_key(|r| match r["status"].as_str() {
        Some("fails") => 0,
        Some("unknown") => 1,
        _ => 2,
    });
    out["pairs"] = J::Array(pairs);
    Ok(out)
}

// ---------------------------------------------------------------- studies

#[cfg(test)]
mod tests {
    use super::*;
    use ink::claims::Status;
    use ink::lang::answer::motion_fields;

    /// A measured claim the time limit stopped says how many of its
    /// undecided moments the limit left (of those listed), with the rest's
    /// own hint (a dead point's) where they have one; all of them the
    /// limit's, the limit's alone.
    #[test]
    fn a_stopped_claim_says_what_the_limit_left_and_keeps_the_other_hints() {
        const DEAD: &str = "a dead point: the solve has no rate at t = 0";
        let hint = |stopped: usize, undecided: usize, own: Option<&str>| {
            let motion = ink::motion::Motion {
                steps: 2,
                undecided: (0..undecided).map(|k| (k as f64, k as f64)).collect(),
                stopped,
                ..ink::motion::Motion::none(ink::ranged::Status::Unknown, (0.0, 1.0))
            };
            let mut r = own.map_or(json!({}), |h| json!({"hint": h}));
            motion_fields(
                &mut r,
                &ink::motion::Stepped { motion, bound: false, ranged: false, riding: false, step: None, basis: Some("a provider's".into()), left: None, left_answer: None },
            );
            r["hint"].as_str().unwrap_or_default().to_string()
        };
        let some = hint(2, 3, Some(DEAD));
        assert!(some.contains("stopped 2 of the 3 in motion.undecided before") && some.ends_with(&format!("The rest: {DEAD}")), "{some}");
        let bare = hint(1, 2, None);
        assert!(bare.contains("stopped 1 of the 2") && bare.contains("The others are undecided for their own reason"), "{bare}");
        let all = hint(3, 3, Some(DEAD));
        assert!(all.contains("stopped it: those in motion.undecided were not decided") && !all.contains("dead point"), "{all}");
        let many = hint(9, 10, None);
        assert!(many.contains("9 of the 10 in motion.undecided (the first 8 listed)"), "{many}");
        assert_eq!(hint(0, 3, Some(DEAD)), DEAD);
    }

    /// An unknown clearance's attained accuracy as it is printed: over all
    /// of the shapes, `accuracy_attained`; examined only in the window,
    /// `accuracy_attained_in_window`; on inexact inputs (a solved value),
    /// none.
    #[test]
    fn an_attained_accuracy_is_printed_with_its_scope() {
        for (unbounded, inexact) in [(false, false), (true, false), (false, true)] {
            let tail = if unbounded { "(- 2000 x)" } else { "1000000" };
            let solve = if inexact { "(solve s (unknowns (u 1)) (== (* 1e-13 u) 0))" } else { "(param u 1)" };
            let src = format!(
                "{solve}
(fn tiny (cx cy) (max (- (abs (- x cx)) .005) (- (abs (- y cy)) .005) (- (abs (- z .345)) .005)))
(fn ball (cx) (- (sqrt (+ (square (- x cx)) (square y) (square z))) 1))
(def a (min (tiny .123 .234) (ball 50) {tail}))
(def b (min (tiny 3.133 3.244) (ball (+ 59 (* 1e-6 u))) {tail}))"
            );
            let doc = ink::load_str(&src, "attained", Path::new(".")).unwrap();
            let field = |name: &str| match doc.def(name).unwrap() {
                Value::Field(f, _) => *f,
                _ => panic!("{name}"),
            };
            let (a, b) = (field("a"), field("b"));
            let ctx = doc.ctx(&[], 0.0).unwrap();
            let mut v = ModelView::new(doc, ctx);
            let c = ink::claims::check(&mut v, ink::lang::Claim::Clearance(a, b, 3.5), None, ClaimOptions::default()).unwrap();
            assert_eq!(c.status, Status::Unknown, "{c:?}");
            let row = check_json(&c, &["u".into()]);
            if inexact {
                assert!(c.inexact.any() && c.attained().is_none(), "{c:?}");
                assert!(row.get("accuracy_attained").is_none() && row.get("accuracy_attained_in_window").is_none(), "{row}");
            } else {
                let (key, other) = if unbounded { ("accuracy_attained_in_window", "accuracy_attained") } else { ("accuracy_attained", "accuracy_attained_in_window") };
                assert!(row[key].as_f64().is_some_and(|a| a > 0.0) && row.get(other).is_none(), "{row}");
            }
        }
    }
}
