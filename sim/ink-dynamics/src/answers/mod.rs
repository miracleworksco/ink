//! Tests for dynamics answers, with shared model fixtures and claim readers.

mod couplings;
mod drive_balance;
mod effort;
mod reactions;
mod simulation;
mod springs;
mod supports;

use crate::{Dynamics, Model, build, gravity};
use ink_lang::Doc;
use ink_lang::claims::{ClaimOptions, Evidence, Status, run_asserts};
use ink_lang::view::ModelView;
use std::path::Path;

static REGISTRY: ink_lang::Registry = ink_lang::Registry::new(&[&Dynamics]);

const BAR: &str = "(body bar (shape (move (box 300 20 10) 150 0 0)) (attr density 2.7g/cm3))";

fn doc(src: &str) -> Doc {
    Doc::load_str_in(&REGISTRY, src, "answers.ink", Path::new("."), &[]).unwrap_or_else(|e| panic!("{} at {}", e.message, e.line))
}

/// The error loading `src`.
fn load_error(src: &str) -> String {
    match Doc::load_str_in(&REGISTRY, src, "answers.ink", Path::new("."), &[]) {
        Ok(_) => panic!("loaded"),
        Err(e) => format!("{} | {}", e.message, e.hint.unwrap_or_default()),
    }
}

fn model(src: &str) -> Model {
    reacting(src, &[], false).unwrap_or_else(|e| panic!("{}", e.message))
}

/// The model of `src` with the reactions of `paths` (`claimed`: as claims
/// ask for them; else as ink effort does).
fn reacting(src: &str, paths: &[&str], claimed: bool) -> Result<Model, ink_lang::InkError> {
    let d = doc(src);
    let paths: Vec<String> = paths.iter().map(|p| p.to_string()).collect();
    build(&d.reloader().unwrap(), "root", gravity(&d).unwrap(), 0.0, &paths, claimed)
}

fn rotor(limit: &str, claim: &str) -> String {
    format!(
        r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.1s)
(body rotor (shape (move (cylinder-z 7.3 41.9) 0 0 5.1)) (attr density 7850kg/m3))
(assembly (spin (a 0))
  (instance r rotor (at (move 3.3 -11.7 0)) (turn z a) (effort-limits -{limit}N*mm {limit}N*mm)))
(assembly root (instance j (spin (* 250000deg/s2 t t))))
{claim}"#
    )
}

/// Each effort claim's status, why where it is unavailable, and how its
/// steps were chosen, as `ink require` checks it: over the file's span, its
/// declared tolerances ranged.
fn claims(src: &str) -> Vec<(String, Status, String, String)> {
    claims_at(src, &[])
}

/// `claims` with params pinned at values (`--param NAME=V:V`).
fn claims_at(src: &str, pinned: &[(&str, f64)]) -> Vec<(String, Status, String, String)> {
    every_claim_at(src, pinned).into_iter().filter(|c| c.0.starts_with("effort-limits:")).collect()
}

/// Every claim's, so.
fn every_claim_at(src: &str, pinned: &[(&str, f64)]) -> Vec<(String, Status, String, String)> {
    asserts_at(src, pinned)
        .into_iter()
        .map(|a| {
            let why = match &a.check.evidence {
                Evidence::Unsolved(e) => format!("{} | {}", e.message, e.hint.clone().unwrap_or_default()),
                _ => String::new(),
            };
            let basis = a.check.motion.as_ref().and_then(|m| m.basis.clone()).unwrap_or_default();
            (a.name, a.check.status, why, basis)
        })
        .collect()
}

/// Every claim's answer, so.
fn asserts_at(src: &str, pinned: &[(&str, f64)]) -> Vec<ink_lang::claims::AssertResult> {
    let overrides: Vec<(String, ink_lang::Given)> = pinned.iter().map(|&(k, v)| (k.to_string(), v.into())).collect();
    let d = Doc::load_str_in(&REGISTRY, src, "answers.ink", Path::new("."), &overrides).unwrap_or_else(|e| panic!("{} at {}", e.message, e.line));
    let span = d.t_span;
    let (ctx, unsolved) = d.ctx_at(&overrides, 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    v.unsolved = unsolved;
    v.set_motion(span).unwrap();
    let spans: Vec<_> = v.declared_spans().unwrap().into_iter().filter(|&(k, _)| !pinned.iter().any(|p| p.0 == v.doc.params[k].name)).collect();
    v.set_param_spans(&spans).unwrap();
    run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap().asserts
}

/// Each claim's status and why where it is unavailable, checked at moment t
/// alone, as `ink require --t` does.
fn claims_at_moment(src: &str, t: f64) -> Vec<(String, Status, String, f64)> {
    let d = doc(src);
    let (ctx, unsolved) = d.ctx_at(&[], t).unwrap();
    let mut v = ModelView::new(d, ctx);
    v.unsolved = unsolved;
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    r.asserts
        .into_iter()
        .map(|a| {
            let why = match &a.check.evidence {
                Evidence::Unsolved(e) => e.message.clone(),
                _ => String::new(),
            };
            let x = match a.check.evidence {
                Evidence::Value(x) => x,
                _ => f64::NAN,
            };
            (a.name, a.check.status, why, x)
        })
        .collect()
}

/// The value a claim's answer gives, checked over the file's span.
fn every_value(src: &str, name: &str) -> f64 {
    let d = doc(src);
    let span = d.t_span;
    let (ctx, unsolved) = d.ctx_at(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    v.unsolved = unsolved;
    v.set_motion(span).unwrap();
    let r = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
    match r.asserts.into_iter().find(|a| a.name == name).unwrap().check.evidence {
        Evidence::Value(x) => x,
        e => panic!("{e:?}"),
    }
}

/// Each claim's value at moment t, as `claims_at_moment` gives them, by name.
/// `src` with each value claim's comparison turned the other way: one that holds over the band of
/// what it reads (its masses') then fails at every value in it, and so answers the reading itself.
fn turned(src: &str) -> String {
    let turn = |line: &str| -> String {
        let Some(end) = line.rfind("))").filter(|_| line.contains("(require ") && line.contains("(value ")) else { return line.to_string() };
        let Some((at, cmp)) = [" < ", " > ", " <= ", " >= "].iter().filter_map(|c| line[..end].rfind(c).map(|i| (i, *c))).max_by_key(|x| x.0) else { return line.to_string() };
        let other = match cmp {
            " < " => " > ",
            " > " => " < ",
            " <= " => " >= ",
            _ => " <= ",
        };
        format!("{}{other}{}", &line[..at], &line[at + cmp.len()..])
    };
    // Each claim of a line on its own: a line may hold several.
    let each = |line: &str| -> String {
        let starts: Vec<usize> = line.match_indices("(require ").map(|m| m.0).collect();
        let Some(&first) = starts.first() else { return line.to_string() };
        let ends = starts.iter().skip(1).copied().chain([line.len()]);
        std::iter::once(line[..first].to_string()).chain(starts.iter().zip(ends).map(|(&a, b)| turn(&line[a..b]))).collect()
    };
    src.lines().map(each).collect::<Vec<_>>().join("\n")
}

/// Each value claim's reading at t, as the model computes it (the claim `turned`): a claim that holds
/// over its band answers that band's worst edge instead.
fn values_at(src: &str, t: f64) -> impl Fn(&str) -> f64 + use<> {
    let c = claims_at_moment(&turned(src), t);
    move |n: &str| c.iter().find(|x| x.0 == n).unwrap_or_else(|| panic!("no claim {n}: {c:?}")).3
}

fn witness(e: &Evidence) -> f64 {
    match e {
        Evidence::Value(x) => *x,
        Evidence::ValueSpan(v) => v.witness.as_ref().expect("a numerical witness").1,
        _ => panic!("expected numerical evidence, got {e:?}"),
    }
}

/// Check both the verdict and the scalar or margin actually reported.
fn claim_value(src: &str, name: &str, status: Status, expected: f64, tolerance: f64) {
    let a = asserts_at(src, &[]).into_iter().find(|a| a.name == name).unwrap().check;
    assert_eq!(a.status, status, "{:?}", a.evidence);
    let x = witness(&a.evidence);
    assert!((x - expected).abs() <= tolerance, "{name}: {x} against {expected}");
}

fn near(got: f64, want: f64, tol: f64, what: &str) {
    assert!((got - want).abs() <= tol, "{what}: {got} against {want} (within {tol})");
}

/// A message's words: a name the test wrote may be looked for there
/// without the kernel's words around it.
fn words(s: &str) -> Vec<&str> {
    s.split(|c: char| !(c.is_alphanumeric() || "-_./:".contains(c))).filter(|w| !w.is_empty()).collect()
}
