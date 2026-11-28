//! What the command-line tests share: the kernel's (its
//! crates/ink-cli/tests/cli/common.rs, beside this folder here and in the
//! release alike), run on ext's ink, with the two roots and what only
//! these tests add.
#![allow(dead_code)]

#[path = "../../../../kernel/crates/ink-cli/tests/cli/common.rs"]
mod kernel;
pub use kernel::*;

/// The two roots: ROOT, the folder holding the crate under test (its
/// fixtures, examples, benchmarks and docs), and KERNEL, the kernel's.
#[path = "roots.rs"]
mod roots;
/// The binary under test, as the test crate names it.
pub use crate::INK;
pub use roots::{KERNEL, ROOT};

/// A fixture of the kernel's tests/fixtures that these tests read too.
pub fn kernel_fixture(name: &str) -> String {
    format!("{KERNEL}/tests/fixtures/{name}")
}

/// A number in an answer (NaN where it is none).
pub fn num(v: &serde_json::Value) -> f64 {
    v.as_f64().unwrap_or(f64::NAN)
}
