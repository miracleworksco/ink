//! ext's docs can't go stale: the README and each package's reference
//! (docs) open with a model you can paste and run, which is a file the
//! model rule checks, and every command they show running is run. A block
//! is marked by an HTML comment on the line before it, as the kernel's
//! docs are (../kernel/crates/ink-cli/tests/cli/docs.rs):
//!
//!   <!-- example PATH -->  ```lisp: PATH, word for word; PATH is in
//!                          examples/, so the model rule checks it
//!   <!-- run -->           ```sh: each `ink ...` line (or `.../ink ...`)
//!                          is run where `examples` is ext's; a ```json
//!                          block right after it is part of what the last
//!                          one answers
//!   <!-- not run: WHY -->  shown and not run
//!
//! A shown answer is matched as the kernel's docs match theirs: each key
//! shown is there with the value shown, a number to the digits shown, a
//! list of objects in order with others between. Each reference's first
//! subsection is its quick start (`### Quick start`), holding an example
//! and a run.

use crate::common::*;
use serde_json::Value;

// The docs checked, the area's own (main.rs).
use crate::{DOCS, REFERENCES};

/// A command line's words, split at spaces (the docs quote none).
fn command(line: &str) -> Option<Vec<String>> {
    let line = line.split(" #").next().unwrap_or("");
    let w: Vec<String> = line.split_whitespace().map(String::from).collect();
    w.first().is_some_and(|c| c == "ink" || c.ends_with("/ink")).then_some(w)
}

/// What is wrong with running a run block at `at`, its answer shown or none.
fn run(at: &str, body: &str, answer: Option<(usize, &str)>, rel: &str) -> Vec<String> {
    let dir = Scratch::new("docs");
    #[cfg(unix)]
    std::os::unix::fs::symlink(root("examples"), dir.join("examples")).unwrap();
    let mut problems = vec![];
    let mut last = None;
    for line in body.lines().filter(|l| !l.trim().is_empty()) {
        let Some(w) = command(line) else {
            problems.push(format!("{at}: `{line}` is run, but is not an ink command (mark the block not run)"));
            continue;
        };
        let r = dir.ink(&w[1..]).cwd(&dir).run();
        match serde_json::from_str::<Value>(&r.stdout) {
            Ok(d) if d["ok"] == true => last = Some(d),
            _ => problems.push(format!("{at}: `{line}` gave no ok answer:\n{}{}", r.stdout, r.stderr)),
        }
    }
    if let (Some((line, want)), Some(got)) = (answer, &last) {
        problems.extend(differs(&shown(want), got, "").into_iter().map(|d| format!("{rel}:{line}: {d}")));
    }
    problems
}

/// What is wrong with the Markdown file at `rel`: its examples, its runs
/// and, for a reference, its quick start.
fn check_doc(rel: &str) -> Vec<String> {
    let text = std::fs::read_to_string(root(rel)).unwrap();
    let bs = blocks(&text);
    let mut problems = vec![];
    if REFERENCES.contains(&rel) {
        // The first subsection is the quick start, and it holds an example and a run.
        let start = text.lines().position(|l| l.starts_with("### "));
        let next = start.and_then(|s| text.lines().skip(s + 1).position(|l| l.starts_with("### ")).map(|k| s + 1 + k + 1));
        match start {
            Some(s) if text.lines().nth(s) == Some("### Quick start") => {
                let inside = |b: &Block| b.line > s && next.is_none_or(|n| b.line < n);
                for kind in ["example", "run"] {
                    if !bs.iter().any(|b| inside(b) && b.marker.as_deref().is_some_and(|m| m.split(' ').next() == Some(kind))) {
                        problems.push(format!("{rel}: its quick start has no {kind} block"));
                    }
                }
            }
            _ => problems.push(format!("{rel}: its first subsection is not ### Quick start")),
        }
    }
    let mut runs = vec![];
    for (i, b) in bs.iter().enumerate() {
        let at = format!("{rel}:{}", b.line);
        let marker = b.marker.as_deref().unwrap_or("");
        let (kind, arg) = marker.split_once(' ').map_or((marker, ""), |(k, a)| (k, a.trim()));
        match (kind, b.lang.as_str()) {
            ("", _) | ("not", _) if marker.is_empty() || arg.starts_with("run:") => {}
            ("example", "lisp") => match std::fs::read_to_string(root(arg)) {
                _ if !arg.starts_with("examples/") => problems.push(format!("{at}: {arg} is not in examples, which the model rule checks")),
                Ok(file) if file == b.body => {}
                Ok(_) => problems.push(format!("{at}: the block is not {arg} word for word")),
                Err(e) => problems.push(format!("{at}: {arg}: {e}")),
            },
            ("run", "sh") => {
                let answer = bs.get(i + 1).filter(|n| n.lang == "json" && n.marker.is_none()).map(|a| (a.line, a.body.as_str()));
                runs.push((at, b.body.as_str(), answer));
            }
            _ => problems.push(format!("{at}: a ```{} block marked {marker:?}: mark it example (```lisp), run (```sh) or not run", b.lang)),
        }
    }
    // In order: a later run reuses what an earlier one solved (ink's cache).
    for (at, body, answer) in runs {
        problems.extend(run(&at, body, answer, rel));
    }
    problems
}

/// Each doc on a thread of its own: a doc's commands may take a minute.
#[test]
fn the_docs_open_with_models_that_run_and_show_what_ink_answers() {
    let docs = DOCS.iter().chain(REFERENCES);
    let problems: Vec<String> = std::thread::scope(|s| {
        let handles: Vec<_> = docs.map(|rel| s.spawn(move || check_doc(rel))).collect();
        handles.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    assert!(problems.is_empty(), "the docs:\n{}", problems.join("\n"));
}
