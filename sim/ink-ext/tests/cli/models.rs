//! The model rule (the kernel's crates/ink-cli/tests/cli/models.rs, over
//! examples and benchmarks): every model passes when `ink check` gives no warning it does not declare and every
//! claim written in it holds, at the resolution ink picks (one that needs
//! another says `(res r)`). A file declares, in comment lines:
//!
//!   ; expect-status: NAME STATUS    a claim meant to fail, or one with no
//!                                   guaranteed accuracy: holds, fails or
//!                                   unknown, `estimate` after holds or fails
//!                                   whose accuracy is null; matched exactly
//!   ; expect-warning: TEXT          a warning starting with TEXT, given
//!   ; expect-error: TEXT            `ink check` refuses it with an error
//!                                   starting with TEXT (a benchmark of a
//!                                   form not yet supported, refused by name)
//!   ; optional-input: PATH          a file it reads, relative to it, kept
//!                                   out of the repository: absent, the
//!                                   model is skipped
//!
//! A model using a package this build lacks is skipped. The quick
//! benchmarks are checked here too; the slow ones (fine studies:
//! minutes) are validation, not tests: benchmarks/run.sh
//! checks them (INK_MODELS names the files to check instead). `SLOW` (over
//! ~20 s each cold) are checked by this rule in the slow set
//! (`./cpu-test.sh --slow`).

use crate::common::*;
// What it checks, the area's own (main.rs).
use crate::{ALONE, QUICK_AT_LEAST, SHAPED_BY_A_STUDY, SLOW};
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

struct Declared {
    status: BTreeMap<String, String>,
    warnings: Vec<String>,
    error: Option<String>,
    optional: Vec<String>,
}

fn declarations(path: &Path) -> Declared {
    let text = std::fs::read_to_string(path).unwrap();
    let mut d = Declared { status: BTreeMap::new(), warnings: vec![], error: None, optional: vec![] };
    for l in text.lines() {
        if let Some(s) = l.strip_prefix("; expect-status: ") {
            let (name, st) = s.trim().split_once(' ').unwrap_or_else(|| panic!("{path:?}: expect-status needs NAME STATUS: {l}"));
            d.status.insert(name.to_string(), st.trim().to_string());
        } else if let Some(w) = l.strip_prefix("; expect-warning: ") {
            if !w.trim().is_empty() {
                d.warnings.push(w.trim().to_string());
            }
        } else if let Some(e) = l.strip_prefix("; expect-error: ") {
            d.error = Some(e.trim().to_string());
        } else if let Some(p) = l.strip_prefix("; optional-input: ") {
            d.optional.push(path.parent().unwrap().join(p.trim()).to_string_lossy().into_owned());
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

/// What is wrong with one model, its claims too where `claims`; None when it is skipped.
fn problems(path: &Path, claims: bool) -> Option<Vec<String>> {
    let dec = declarations(path);
    let file = path.to_str().unwrap();
    let check = Ink::new(["check", file]).run().json();
    if check["ok"].as_bool() != Some(true) {
        let m = message(&check);
        let absent = m.strip_prefix("cannot read ").and_then(|r| r.split_once(": No such file")).map(|(f, _)| f.to_string());
        if absent.is_some_and(|f| dec.optional.iter().any(|o| Path::new(o) == Path::new(&f) || o.ends_with(&f))) || not_built(&m) {
            return None;
        }
        if dec.error.as_ref().is_some_and(|e| m.starts_with(e.as_str())) {
            return Some(vec![]);
        }
        return Some(vec![format!("check failed: {m}")]);
    }
    if let Some(e) = &dec.error {
        return Some(vec![format!("expected error not given: {e}")]);
    }
    // A warning naming a file in this area names it from the area's root, as its declaration does.
    let root = format!("{ROOT}/");
    let warnings: Vec<&str> = check["warnings"].as_array().map_or(vec![], |w| w.iter().filter_map(Value::as_str).map(|w| w.strip_prefix(root.as_str()).unwrap_or(w)).collect());
    let mut out: Vec<String> = warnings.iter().filter(|w| !dec.warnings.iter().any(|e| w.starts_with(e.as_str()))).map(|w| format!("warning: {w}")).collect();
    out.extend(dec.warnings.iter().filter(|e| !warnings.iter().any(|w| w.starts_with(e.as_str()))).map(|e| format!("expected warning not given: {e}")));
    if !claims {
        return Some(out);
    }
    let listed: Vec<&str> = check["requirements"].as_array().map_or(vec![], |a| a.iter().filter(|c| c["own"].as_bool() == Some(true)).filter_map(|c| c["name"].as_str()).collect());
    let rows = if listed.is_empty() {
        vec![]
    } else {
        let a = Ink::new(["require", file]).run().json();
        if a["ok"].as_bool() != Some(true) {
            out.push(format!("assert failed: {}", message(&a)));
            return Some(out);
        }
        let rows = a["requirements"].as_array().cloned().unwrap_or_default();
        // Its counts are of its rows' verdicts.
        let mut counts = BTreeMap::<String, u64>::new();
        for r in &rows {
            *counts.entry(r["verdict"].as_str().unwrap_or("no answer").to_string()).or_default() += 1;
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
        let holds = r["status"].as_str() == Some("holds") || (r.get("status").is_none() && r["verdict"].as_str() == Some("yes"));
        match dec.status.get(name) {
            None if !holds => out.push(format!("{name}: {st}")),
            Some(want) if *want != st => out.push(format!("{name}: {st}, expected {want}")),
            _ => {}
        }
    }
    out.extend(dec.status.iter().filter(|(n, _)| !got.contains_key(n.as_str())).map(|(n, st)| format!("{n}: expected {st}, but no such claim")));
    Some(out)
}

/// The models checked: examples' quick ones, else (`slow`) SLOW.
fn models(slow: bool) -> Vec<String> {
    if let Ok(list) = std::env::var("INK_MODELS") {
        return list.split_whitespace().map(|f| if Path::new(f).is_absolute() { f.to_string() } else { root(f) }).collect();
    }
    let mut out = vec![];
    for dir in ["examples", "benchmarks"] {
        let Ok(entries) = std::fs::read_dir(root(dir)) else { continue };
        for e in entries {
            let p = e.unwrap().path();
            let rel = format!("{dir}/{}", p.file_name().unwrap().to_string_lossy());
            if p.extension().is_some_and(|x| x == "ink") && !ALONE.contains(&rel.as_str()) && SLOW.iter().any(|s| rel.starts_with(s)) == slow {
                out.push(p.to_string_lossy().into_owned());
            }
        }
    }
    out.sort();
    out
}

/// How long each model took on the last run here: the slowest start first,
/// so the run ends with short ones, not one long one alone.
const SECONDS: &str = concat!(env!("CARGO_TARGET_TMPDIR"), "/model-seconds.json");

#[test]
fn every_model_passes_the_model_rule() {
    check(models(false), QUICK_AT_LEAST, 4, true);
}

#[test]
fn every_slow_model_gives_only_the_warnings_it_declares() {
    // Their claims take minutes, their check seconds (24 s cold for these 13,
    // 2026-11-23): a new warning firing on one is found the day it lands,
    // not in the slow set.
    let mut slow: Vec<String> = models(true).into_iter().filter(|f| !SHAPED_BY_A_STUDY.iter().any(|s| f.starts_with(&root(s)))).collect();
    slow.extend(ALONE.iter().map(|a| root(a)));
    check(slow, SLOW.len() + ALONE.len() - SHAPED_BY_A_STUDY.len(), 4, false);
}

#[test]
#[cfg_attr(not(feature = "slow"), ignore = "slow: ./cpu-test.sh --slow")]
fn slow_models_pass_the_model_rule() {
    // (Half as many at once as the quick ones: eight of these together
    // took the slow set past the memory one job is given.)
    check(models(true), SLOW.len(), 8, true);
}

#[test]
#[ignore = "needs 18 GB and 8 minutes: run alone, by name, with that memory"]
fn the_snap_fit_passes_the_model_rule_alone() {
    check(ALONE.iter().map(|a| root(a)).collect(), ALONE.len(), 4, true);
}

/// `share`: the cores to each model solved at once; `claims`: its claims too, not its warnings alone.
pub(super) fn check(mut files: Vec<String>, at_least: usize, share: usize, claims: bool) {
    assert!(files.len() >= at_least || std::env::var_os("INK_MODELS").is_some(), "found only {} models", files.len());
    let before: BTreeMap<String, f64> = std::fs::read(SECONDS).ok().and_then(|b| serde_json::from_slice(&b).ok()).unwrap_or_default();
    files.sort_by(|a, b| before.get(b).unwrap_or(&f64::MAX).total_cmp(before.get(a).unwrap_or(&f64::MAX)));
    let next = AtomicUsize::new(0);
    let failed = Mutex::new(vec![]);
    let skipped = Mutex::new(vec![]);
    let seconds = Mutex::new(BTreeMap::new());
    let workers = std::thread::available_parallelism().map_or(4, |n| n.get() / share).max(1);
    std::thread::scope(|s| {
        for _ in 0..workers {
            s.spawn(|| {
                while let Some(f) = files.get(next.fetch_add(1, Ordering::Relaxed)) {
                    let rel = f.strip_prefix(&format!("{ROOT}/")).unwrap_or(f).to_string();
                    let t0 = std::time::Instant::now();
                    let p = problems(Path::new(f), claims);
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
    // (a check alone is no measure of what the model takes)
    if claims {
        _ = std::fs::write(SECONDS, serde_json::to_vec(&seconds).unwrap());
    }
    let skipped = skipped.into_inner().unwrap();
    if !skipped.is_empty() {
        eprintln!("skipped (an optional input absent, or a package this build lacks): {}", skipped.join(" "));
    }
    let failed = failed.into_inner().unwrap();
    assert!(failed.is_empty(), "{} of {} models fail the model rule:\n{}", failed.len(), files.len(), failed.join("\n"));
}
