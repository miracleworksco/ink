//! The model rule: every model in examples/, benchmarks/ and tests/fixtures
//! passes when `ink check` gives no warning it does not declare and every
//! claim written in it holds, at the resolution ink picks (one that needs
//! another says `(res r)`). A file declares, in comment lines:
//!
//!   ; expect-status: NAME STATUS    a claim meant to fail, or one with no
//!                                   guaranteed accuracy: holds, fails or
//!                                   unknown, `estimate` after holds or fails
//!                                   whose accuracy is null; matched exactly
//!   ; expect-warning: TEXT          a warning containing TEXT, given (one
//!                                   placed by its line starts with the path)
//!
//! A model using a package this build lacks is skipped. The models in
//! `SLOW` (over ~20 s each cold) are validation, not tests: `benchmarks/run.sh`
//! checks them by this rule (INK_MODELS=slow; INK_MODELS names the files to
//! check instead).

use crate::common::*;
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Models over ~20 s each cold (an empty kernel cache, 4 threads): published
/// validation and heavy examples, checked by benchmarks/run.sh, not per
/// change. A path, or the start of one.
const SLOW: &[&str] =
    &["benchmarks/02-", "benchmarks/10-", "benchmarks/11-", "benchmarks/18-", "benchmarks/20-", "benchmarks/22-", "benchmarks/26-", "benchmarks/32-", "benchmarks/33-"];
/// Fixtures that are not models of their own: an input the import test
/// copies beside the JSON it reads.
const NOT_MODELS: &[&str] = &["tests/fixtures/import-summary.ink"];

struct Declared {
    status: BTreeMap<String, String>,
    warnings: Vec<String>,
}

fn declarations(path: &Path) -> Declared {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    let mut d = Declared { status: BTreeMap::new(), warnings: vec![] };
    for l in text.lines() {
        if let Some(s) = l.strip_prefix("; expect-status: ") {
            let (name, st) = s.trim().split_once(' ').unwrap_or_else(|| panic!("{path:?}: expect-status needs NAME STATUS: {l}"));
            d.status.insert(name.to_string(), st.trim().to_string());
        } else if let Some(w) = l.strip_prefix("; expect-warning: ")
            && !w.trim().is_empty()
        {
            d.warnings.push(w.trim().to_string());
        }
    }
    d
}

fn message(d: &Value) -> String {
    match &d["error"] {
        Value::Object(e) => e.get("message").map_or_else(|| Value::Object(e.clone()).to_string(), |m| m.as_str().unwrap_or("").to_string()),
        Value::String(s) => s.clone(),
        e => format!("no error given, and not ok: {e}"),
    }
}

/// What is wrong with one model; None when it is skipped.
fn problems(path: &Path) -> Option<Vec<String>> {
    let dec = declarations(path);
    let file = path.to_str().unwrap();
    let check = Ink::new(["check", file]).run().json();
    if check["ok"].as_bool() != Some(true) {
        let m = message(&check);
        if not_built(&m) {
            return None;
        }
        return Some(vec![format!("check failed: {m}")]);
    }
    let warnings: Vec<&str> = check["warnings"].as_array().map_or(vec![], |w| w.iter().filter_map(Value::as_str).collect());
    let mut out: Vec<String> = warnings.iter().filter(|w| !dec.warnings.iter().any(|e| w.contains(e.as_str()))).map(|w| format!("warning: {w}")).collect();
    out.extend(dec.warnings.iter().filter(|e| !warnings.iter().any(|w| w.contains(e.as_str()))).map(|e| format!("expected warning not given: {e}")));
    let listed: Vec<&str> = check["asserts"].as_array().map_or(vec![], |a| a.iter().filter(|c| c["own"].as_bool() == Some(true)).filter_map(|c| c["name"].as_str()).collect());
    let rows = if listed.is_empty() {
        vec![]
    } else {
        let a = Ink::new(["assert", file, "--own"]).run().json();
        if a["ok"].as_bool() != Some(true) {
            out.push(format!("assert failed: {}", message(&a)));
            return Some(out);
        }
        if a["scope"].as_str() != Some("own") {
            out.push("assert does not say it checked the file's own claims (scope: own)".into());
        }
        let rows = a["asserts"].as_array().cloned().unwrap_or_default();
        let mut counts = BTreeMap::<String, u64>::new();
        for r in &rows {
            let st = r["status"].as_str().filter(|s| matches!(*s, "holds" | "fails")).unwrap_or("unknown");
            *counts.entry(st.to_string()).or_default() += 1;
        }
        let given: BTreeMap<String, u64> =
            a["counts"].as_object().map_or_else(BTreeMap::new, |c| c.iter().map(|(k, v)| (k.clone(), v.as_u64().unwrap_or(0))).filter(|(_, v)| *v > 0).collect());
        if given != counts {
            out.push(format!("counts {given:?} disagree with the claims {counts:?}"));
        }
        rows
    };
    let mut got = BTreeMap::<&str, usize>::new();
    for r in &rows {
        *got.entry(r["name"].as_str().unwrap_or("?")).or_default() += 1;
    }
    out.extend(listed.iter().filter(|n| !got.contains_key(*n)).map(|n| format!("{n}: missing from assert")));
    out.extend(got.iter().filter(|(_, c)| **c > 1).map(|(n, c)| format!("{n}: reported {c} times")));
    out.extend(got.keys().filter(|n| !listed.contains(n)).map(|n| format!("{n}: not a claim check lists")));
    for r in &rows {
        let (name, st) = (r["name"].as_str().unwrap_or("?"), verdict(r));
        match dec.status.get(name) {
            None if r["status"].as_str() != Some("holds") => out.push(format!("{name}: {st}")),
            Some(want) if *want != st => out.push(format!("{name}: {st}, expected {want}")),
            _ => {}
        }
    }
    out.extend(dec.status.iter().filter(|(n, _)| !got.contains_key(n.as_str())).map(|(n, st)| format!("{n}: expected {st}, but no such claim")));
    Some(out)
}

/// The .ink files under examples/, benchmarks/ and tests/fixtures/, as
/// paths from the root.
fn all_models() -> Vec<String> {
    let mut out = vec![];
    for dir in ["examples", "benchmarks", "tests/fixtures"] {
        let Ok(entries) = std::fs::read_dir(root(dir)) else { continue };
        for e in entries {
            let name = e.unwrap().file_name().to_string_lossy().into_owned();
            if name.ends_with(".ink") {
                out.push(format!("{dir}/{name}"));
            }
        }
    }
    out
}

#[test]
fn every_slow_model_named_is_there() {
    // Each name is a model here: cargo test filters SLOW and never opens
    // it, and benchmarks/run.sh, which runs SLOW, would panic on a missing
    // one.
    let all = all_models();
    let gone: Vec<&&str> = SLOW.iter().chain(NOT_MODELS).filter(|s| !all.iter().any(|m| m.starts_with(**s))).collect();
    assert!(gone.is_empty(), "SLOW or NOT_MODELS in models.rs names no model here: {gone:?}");
}

fn models() -> Vec<String> {
    // The slow ones, read here rather than parsed out of this file.
    if std::env::var("INK_MODELS").as_deref() == Ok("slow") {
        let mut out: Vec<String> = all_models().into_iter().filter(|rel| SLOW.iter().any(|s| rel.starts_with(s))).map(|rel| root(&rel)).collect();
        out.sort();
        return out;
    }
    if let Ok(list) = std::env::var("INK_MODELS") {
        return list.split_whitespace().map(|f| if Path::new(f).is_absolute() { f.to_string() } else { root(f) }).collect();
    }
    let mut out: Vec<String> = all_models().into_iter().filter(|rel| !SLOW.iter().chain(NOT_MODELS).any(|s| rel.starts_with(s))).map(|rel| root(&rel)).collect();
    out.sort();
    out
}

/// How long each model took on the last run here: the slowest start first,
/// so the run ends with short ones, not one long one alone.
const SECONDS: &str = concat!(env!("CARGO_TARGET_TMPDIR"), "/model-seconds.json");

#[test]
fn every_model_passes_the_model_rule() {
    let mut files = models();
    assert!(files.len() > 35 || std::env::var_os("INK_MODELS").is_some(), "found only {} models", files.len());
    let before: BTreeMap<String, f64> = std::fs::read(SECONDS).ok().and_then(|b| serde_json::from_slice(&b).ok()).unwrap_or_default();
    files.sort_by(|a, b| before.get(b).unwrap_or(&f64::MAX).total_cmp(before.get(a).unwrap_or(&f64::MAX)));
    let next = AtomicUsize::new(0);
    let failed = Mutex::new(vec![]);
    let skipped = Mutex::new(vec![]);
    let seconds = Mutex::new(BTreeMap::new());
    let workers = std::thread::available_parallelism().map_or(4, |n| n.get() / 4).max(1);
    std::thread::scope(|s| {
        for _ in 0..workers {
            s.spawn(|| {
                while let Some(f) = files.get(next.fetch_add(1, Ordering::Relaxed)) {
                    let rel = f.strip_prefix(&format!("{ROOT}/")).unwrap_or(f).to_string();
                    let t0 = std::time::Instant::now();
                    let p = problems(Path::new(f));
                    seconds.lock().unwrap().insert(f.clone(), t0.elapsed().as_secs_f64());
                    match p {
                        None => skipped.lock().unwrap().push(rel),
                        Some(p) if !p.is_empty() => failed.lock().unwrap().push(format!("{rel}:\n    {}", p.join("\n    "))),
                        Some(_) => {}
                    }
                }
            });
        }
    });
    let mut seconds = seconds.into_inner().unwrap();
    for (f, t) in before {
        seconds.entry(f).or_insert(t);
    }
    _ = std::fs::write(SECONDS, serde_json::to_vec(&seconds).unwrap());
    let skipped = skipped.into_inner().unwrap();
    if !skipped.is_empty() {
        eprintln!("skipped (a package this build lacks): {}", skipped.join(" "));
    }
    let failed = failed.into_inner().unwrap();
    assert!(failed.is_empty(), "{} of {} models fail the model rule:\n{}", failed.len(), files.len(), failed.join("\n"));
}
