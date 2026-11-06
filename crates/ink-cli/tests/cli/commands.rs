//! The kernel's commands end to end: refusals and their hints, eval,
//! bounds, mass's volume, slice (and its DXF), mesh (STL, 3MF, GLB), components,
//! ranges of the question (access), interfere and assert. The files ink
//! writes are read here without the kernel.

mod access;
mod claims;
mod interfere;
mod mesh;
mod optimize;
mod queries;
mod zip;

use crate::common::*;
use serde_json::Value;

/// tests/fixtures/sphere.ink: a ball of (param r 20).
pub fn sphere() -> String {
    fixture("sphere.ink")
}

/// tests/fixtures/shelf.ink: brackets, screws and a board.
pub fn shelf() -> String {
    fixture("shelf.ink")
}

/// A number, or a failure naming where it was looked for.
pub fn num(v: &Value) -> f64 {
    v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"))
}

/// A JSON array of numbers (ints and floats alike), or None.
pub fn nums(v: &Value) -> Option<Vec<f64>> {
    v.as_array()?.iter().map(Value::as_f64).collect()
}

/// The field's value of `target` at `point` (JSON [x,y,z]), with these extra arguments.
pub fn value_at(model: &str, target: &str, point: &Value, extra: &[&str]) -> f64 {
    let pts = format!("[{point}]");
    let mut args = vec!["eval", model, "--target", target, "--points", &pts];
    args.extend_from_slice(extra);
    num(&records(&ink(args)["targets"][0]["results"])[0]["value"])
}

/// The field's largest value over these points of `target`.
pub fn max_value_at(model: &str, target: &str, points: &[Value]) -> f64 {
    let pts = Value::Array(points.to_vec()).to_string();
    let d = ink(["eval", model, "--target", target, "--points", &pts]);
    records(&d["targets"][0]["results"]).iter().map(|r| num(&r["value"])).fold(f64::MIN, f64::max)
}

/// INK_NO_CACHE is set: checks of what the cache keeps are skipped.
pub fn no_cache() -> bool {
    std::env::var_os("INK_NO_CACHE").is_some()
}

/// The files under a cache directory's answers (kernel answers stored).
pub fn files_under(dir: &std::path::Path) -> usize {
    let Ok(entries) = std::fs::read_dir(dir) else { return 0 };
    entries.map(|e| e.unwrap().path()).map(|p| if p.is_dir() { files_under(&p) } else { 1 }).sum()
}
