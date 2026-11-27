//! ext's `ink` (ink-ext: the kernel's packages and ext's: sim,
//! dynamics, topology) end to end: the
//! binary cargo built for this test (CARGO_BIN_EXE_ink) run on models the
//! tests write and tests/fixtures, its JSON and files checked as an
//! agent reads them. `models` checks examples and benchmarks by the model rule. One
//! binary, so the tests run side by side.

/// The binary under test.
pub const INK: &str = env!("CARGO_BIN_EXE_ink");

mod benchmarks;
mod common;
mod docs;
mod dynamics;
mod help;
mod models;
mod require;
mod studies;
mod topology;
mod withheld;

// What models.rs and docs.rs check: this area's (they are included by path
// by the areas built on this one, each naming its own).

/// The docs checked: the README, and the references, whose first
/// subsection is a quick start.
pub const DOCS: &[&str] = &["README.md"];
pub const REFERENCES: &[&str] = &["docs/SIM.md", "docs/DYNAMICS.md", "docs/TOPOLOGY.md"];

/// How many quick models there are at least (a check that finds fewer
/// found the wrong folder).
pub const QUICK_AT_LEAST: usize = 4;

/// Models over ~20 s each cold (the shelf, the field-driven example and
/// the twisted blade, 45 to 142 s in the 2026-11-18 census): the topology example
/// (its 60-iteration layout, 30 s at 4 threads) and the fine studies
/// the kernel's SLOW held before simulation moved here. A path, or the
/// start of one.
pub const SLOW: &[&str] = &[
    "examples/shelf.ink",
    "examples/field-driven.ink",
    "benchmarks/03-",
    "examples/topology.ink",
    "examples/gyroid-blade-mode.ink",
    "benchmarks/08-",
    "benchmarks/12-",
    "benchmarks/13-",
    "benchmarks/15-",
    "benchmarks/21-",
    "benchmarks/25-",
    "benchmarks/28-",
    "benchmarks/30-",
];

/// Models whose claims take more memory than one job is given here (the
/// snap fit's 5 million cells: 18 GB cold, against 16): in neither set,
/// checked alone by their own test, run by name with the memory.
pub const ALONE: &[&str] = &["benchmarks/24-snap-fit.ink"];

/// Slow models whose bodies are a study's result (a topology layout, a mode's shape): `ink check`
/// solves the study to know the shape, a minute or more cold, so their warnings are the slow set's.
pub const SHAPED_BY_A_STUDY: &[&str] = &["examples/topology.ink", "examples/gyroid-blade-mode.ink", "benchmarks/15-"];

fn standard_models(dir: &std::path::Path) -> Vec<String> {
    const AT_LEAST: usize = 40;
    let files: Vec<String> = std::fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("cannot read standard models in {}: {e}", dir.display()))
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "ink"))
        .map(|p| p.to_string_lossy().into_owned())
        .collect();
    assert!(files.len() >= AT_LEAST, "found only {} standard models in {}, expected at least {AT_LEAST}", files.len(), dir.display());
    files
}

#[test]
#[should_panic(expected = "cannot read standard models in")]
fn a_declared_standard_set_requires_its_directory() {
    let scratch = common::Scratch::new("standard-models");
    standard_models(std::path::Path::new(&scratch.path("missing")));
}

#[test]
fn every_standard_model_gives_only_the_warnings_it_declares() {
    // The standard sets are measured before a release (benchmarks/standard/run.py).
    // Check that they still load on each change.
    let dir = common::root("benchmarks/standard");
    models::check(standard_models(std::path::Path::new(&dir)), 0, 4, false);
}
