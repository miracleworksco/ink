//! The prover's answer to a claim as a fact of the claims core (`fact`):
//! what `ink require` says of a claim beside the prover's own evidence.
//!
//! A fact carries the answer's verdict and numbers, and says what the
//! prover's row leaves to its reader:
//!
//! - how strong: proven where the check guarantees an accuracy and needed
//!   none of it (it held outright, or the accuracy is 0), proven to its
//!   accuracy where it guarantees one otherwise, an estimate where it
//!   guarantees none;
//! - over how much of each variable it reads: t at one value where the
//!   claim reads t and was asked at one moment (`Check::one_moment`), at
//!   every value where its motion was bounded (t ranged, or asked once on
//!   one moving instance: `Stepped::bound`), at so many where it was
//!   measured at steps; a ranged param at every value where the answer is
//!   guaranteed (`Check::ranged`: proven or disproven over all of them), at
//!   so many where an estimate says how many it sampled, else at one;
//! - what it left out: stretches of the motion undecided or unchecked, and
//!   the shape beyond the window where it is not proven bounded.
//!
//! An unknown is not settled or not available as the code that could not
//! answer says, never as its words read: a limit lacking an input, a shape
//! not proven bounded, or a solve refusing for what the design lacks (a
//! moving body with no mass: `InkError::open`), is something the design
//! lacks; a number outside ink's range, a method asked over a range it
//! answers at one value of, a rate with no bound, or a solve refusing what
//! its model does not cover (a body it leaves out), is outside what the
//! source covers. Every other unknown is a source that could not decide,
//! its reason in the prover's own words: a search that ran out, a
//! measurement that does not decide, a solve that says no kind.

use ink_claims::{Covers, Fact, Strength, Unavailable, Unsettled, Verdict};
use ink_lang::Open;
use ink_lang::claims::{AssertResult, Evidence, Held, Status};

/// The prover's answer to one claim as a fact of the core: its verdict and
/// kind, how strong, over how much of what it reads, what it left out.
/// `params` are the design's param names, by index.
pub fn fact(params: &[String], a: &AssertResult) -> Fact {
    let check = &a.check;
    let mut fact = Fact { about: ink_claims::about(&a.name, ""), from: a.source.clone(), terms: a.terms.clone(), notes: vec![a.claim.clone()], ..Fact::default() };
    fact.notes.extend(a.declared_at.iter().map(|at| format!("its limit is declared at {at}")));
    fact.notes.extend(check.note.map(String::from));
    // The numbers of its statement, and its reason in the old answer's words.
    let (mut sampled, mut why) = (None, None);
    match &check.evidence {
        Evidence::Value(x) => fact.terms.push(("value".into(), *x)),
        Evidence::Metric { measured, undecided, estimate, grid_error, .. } => {
            fact.terms.extend(measured.value.iter().map(|x| ("value".into(), *x)));
            fact.terms.extend(grid_error.map(|e| ("grid error".into(), e.0)));
            fact.notes.extend(estimate.clone());
            fact.notes.extend(measured.value.as_ref().err().cloned());
            fact.notes.extend(undecided.iter().map(|(why, how)| format!("{why}. {how}")));
        }
        Evidence::Integral { sampled: s, .. } | Evidence::Mean { sampled: s, .. } => sampled = s.as_ref().map(|s| s.values),
        Evidence::NotChecked(reason) => why = Some((Verdict::NotSettled(Unsettled::Lacks), reason.to_string())),
        Evidence::Unexamined => why = Some((Verdict::NotSettled(Unsettled::Lacks), "it needs all of a shape not proven bounded: nothing was examined".into())),
        Evidence::Outside(reason) => why = Some((Verdict::NotAvailable(Unavailable::Outside), reason.clone())),
        Evidence::NotCovered(reason) => why = Some((Verdict::NotAvailable(Unavailable::Outside), reason.clone())),
        Evidence::Unproven(reason) => fact.notes.push(reason.clone()),
        Evidence::Unsolved(e) => {
            let kind = match e.open {
                Some(Open::Lacks) => Verdict::NotSettled(Unsettled::Lacks),
                Some(Open::Outside) => Verdict::NotAvailable(Unavailable::Outside),
                None => Verdict::NotSettled(Unsettled::Undecided),
            };
            why = Some((kind, e.message.clone()));
        }
        Evidence::Stopped => fact.notes.push("the call's time limit had passed".into()),
        _ => {}
    }
    fact.terms.extend(check.accuracy.map(|x| ("accuracy".into(), x)));
    let inexact = &check.inexact;
    if !inexact.samples.is_empty() || !inexact.solved.is_empty() {
        fact.notes.push(format!("it reads inexact inputs: {} sampled fields, {} solved values not certified", inexact.samples.len(), inexact.solved.len()));
    }
    if !inexact.typical.is_empty() {
        fact.notes.push(format!(
            "it reads a maker's typical strength, not a specification's minimum, so it holds of a typical coupon, not as a guarantee: {}",
            inexact.typical.join(", ")
        ));
    }
    fact.notes.extend(inexact.unsettled.iter().map(|(k, why)| format!("'{}' falls short of its accuracy: {why}", params[*k as usize])));
    fact.notes.extend(check.motion.as_ref().and_then(|s| s.unsettled()));
    fact.held = Some(match check.status {
        Status::Holds => Verdict::Yes,
        Status::Fails => Verdict::No,
        Status::Unknown => {
            let (verdict, reason) = why.unwrap_or((Verdict::NotSettled(Unsettled::Undecided), String::new()));
            fact.notes.extend(Some(reason).filter(|r| !r.is_empty()));
            return Fact { held: Some(verdict), ..fact };
        }
    });
    fact.strength = Some(match (check.accuracy, &check.held) {
        (None, _) => Strength::Estimated,
        (Some(_), Some(Held::Outright)) => Strength::Proven,
        (Some(0.0), None) => Strength::Proven,
        (Some(_), _) => Strength::ToAccuracy,
    });
    let spans = |what: &str, spans: &[(f64, f64)]| spans.iter().map(|(a, b)| format!("t {a} to {b}: {what}")).collect::<Vec<_>>();
    if let Some(t) = check.one_moment {
        fact.covers.push(("t".into(), Covers::One));
        fact.notes.push(format!("asked at t = {t} only: it reads t, and no (t-span) covers it"));
    }
    if let Some(stepped) = &check.motion {
        let motion = &stepped.motion;
        fact.covers.push(("t".into(), if stepped.bound { Covers::Every } else { Covers::These(motion.steps) }));
        fact.left_out
            .extend(spans("undecided", &motion.undecided).into_iter().chain(spans("not checked", &motion.unchecked)).chain(spans("held at its middle only", &motion.loose)));
    }
    for (k, (low, high)) in &check.ranged {
        let name = &params[*k as usize];
        if fact.covers.iter().any(|(var, _)| var == name) {
            continue;
        }
        let covers = match (check.accuracy, sampled) {
            (Some(_), _) => Covers::Every,
            (None, Some(values)) => Covers::These(values),
            (None, None) => {
                fact.left_out.push(format!("'{name}' from {low} to {high}: the check does not say at which values it was made"));
                Covers::One
            }
        };
        fact.covers.push((name.clone(), covers));
    }
    if check.unbounded {
        fact.left_out.push("the shape beyond the window examined: it is not proven bounded".into());
    }
    fact
}
