//! What ext's `ink` says of itself: its version, and its help, each
//! command's with what it does and each example in it one that runs.

use crate::common::*;
use serde_json::Value;

#[test]
fn version_names_the_distribution_and_its_kernel() {
    // `ink --version` was refused as an unknown command.
    let v = Ink::new(["--version"]).text();
    // Its own name and version first; the kernel's version after, in the
    // kernel's words.
    assert!(v.starts_with(&format!("ink-ext {} ", env!("CARGO_PKG_VERSION"))) && v.contains(ink_cli::VERSION) && v.lines().count() == 1, "{v}");
}

/// The commands `ink help` lists: the core's, then the packages'.
fn commands() -> Vec<String> {
    let help = Ink::new(["help"]).text();
    help.lines()
        .skip_while(|l| *l != "commands:")
        .take_while(|l| !l.starts_with("options"))
        .filter_map(|l| l.strip_prefix("  ").filter(|r| !r.starts_with(' ')).and_then(|r| r.split_whitespace().next()).map(String::from))
        .collect()
}

#[test]
fn every_command_has_help_and_its_examples_run() {
    // `ink help solve` said "FEA (elastic or heat)". Each
    // package's command now says what it answers and gives an example,
    // run here as written, from this area's folder.
    let names = commands();
    let packages = ["solve", "effort", "topology"];
    assert!(names.len() > 15 && packages.iter().all(|p| names.iter().any(|n| n == p)), "ink help lists {names:?}");
    let mut problems = vec![];
    let mut examples = vec![];
    for name in &names {
        let r = Ink::new(["help", name]).run();
        if r.code != Some(0) || r.stdout.split_whitespace().count() < 5 {
            problems.push(format!("ink help {name} gives no help (exit {:?}): {}{}", r.code, r.stdout, r.stderr));
        }
        let own: Vec<String> = r.stdout.lines().filter_map(|l| l.trim().strip_prefix("example: ")).map(String::from).collect();
        if packages.contains(&name.as_str()) && own.is_empty() {
            problems.push(format!("ink help {name} gives no example"));
        }
        examples.extend(own.into_iter().map(|e| (name.clone(), e)));
    }
    let ran: Vec<String> = std::thread::scope(|s| {
        let runs: Vec<_> = examples
            .iter()
            .map(|(name, example)| {
                s.spawn(move || {
                    let w: Vec<&str> = example.split_whitespace().collect();
                    if w.first() != Some(&"ink") || w.get(1) != Some(&name.as_str()) {
                        return Some(format!("ink help {name}'s example is not `ink {name} ...`: {example}"));
                    }
                    let r = Ink::new(&w[1..]).cwd(ROOT).run();
                    match serde_json::from_str::<Value>(&r.stdout) {
                        Ok(d) if d["ok"] == true => None,
                        _ => Some(format!("ink help {name}'s example, {example}, gave no ok answer:\n{}{}", r.stdout, r.stderr)),
                    }
                })
            })
            .collect();
        runs.into_iter().filter_map(|h| h.join().unwrap()).collect()
    });
    problems.extend(ran);
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}
