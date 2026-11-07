//! `assert` and `interfere`: claims checked by the library (ink::claims,
//! ink::lang::interfere), printed as JSON.

use crate::*;
use ink::claims::{AssertReport, AssertResult, ClaimOptions, ResBy};
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
    r
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
    let rep = ink::claims::run_asserts_with(s, only, own, ClaimOptions { res, ..Default::default() }, |a| Some(assert_json(a, &names)))?;
    let mut out = report_json(&rep, &names);
    extend(&mut out, context);
    Ok(out)
}

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
