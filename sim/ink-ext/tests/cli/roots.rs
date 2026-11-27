//! Where the tests find their files: two folders, side by side. Every test
//! here names a path through these (ink-ext/tests/cli/common.rs for the
//! command-line tests; the others include this file), so a crate one folder
//! below either root finds both wherever the pair is checked out.
#![allow(dead_code)]

/// This folder of packages, the one holding the crate under test: its
/// examples, benchmarks, docs and tests/fixtures.
pub const ROOT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/..");
/// The kernel's folder beside it: its tests/fixtures and docs.
pub const KERNEL: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../kernel");
