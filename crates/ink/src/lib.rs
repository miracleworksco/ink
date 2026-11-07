//! Implicit Kernel: an implicit geometry kernel built for AI agents.
//! See docs/CORE.md for the representation and its rationale.

// The kernel's packages. Nothing outside this file names one.
/// The design package (crates/ink-design): solve, optimize.
pub use ink_design as design;

/// The language (crates/ink-lang).
pub use ink_lang as lang;
pub use ink_lang::{InkError, materials, sexp, units};

// The crates, re-exported so `ink::graph::...` paths keep working.
pub use ink_core::{V3, eval, graph, interval, tape, v3};
pub use ink_io::{export, gltf, render, seen, simplify};
pub use ink_kernel::{cache, centred, feature, mesh, octree, query, ranged, rigid, runners, scope, stop, voids};

/// Sampled fields (ink-core) and meshes read into fields (ink-io).
pub mod grid {
    pub use ink_core::grid::*;
    pub use ink_io::stl::{MeshDistance, read_stl};
}

/// The kernel's packages: the extension forms every `ink` document can use
/// (see lang::Registry). A distribution with more makes its own registry
/// of these and its own.
pub fn registry() -> &'static lang::Registry {
    static REGISTRY: lang::Registry = lang::Registry::new(PACKAGES);
    &REGISTRY
}

/// The kernel's packages, for a distribution's registry to begin with.
pub const PACKAGES: &[&'static dyn lang::Package] = &[&design::Design];

/// A loaded model and what analysis has learned of it, and checking a
/// document's claims on it (ink-lang).
pub use ink_lang::{claims, motion, validity, view};

/// The language reference, `ink spec`.
pub const SPEC: &str = include_str!("../../../docs/IR.md");

/// Load a file with the standard distribution's packages.
pub fn load_file(path: &std::path::Path) -> Result<lang::Doc, InkError> {
    load_file_with(path, &[])
}

/// Load a file with the standard distribution's packages and parameter
/// overrides applied to load-time computations.
pub fn load_file_with(path: &std::path::Path, overrides: &[(String, lang::Given)]) -> Result<lang::Doc, InkError> {
    lang::Doc::load_file_in(registry(), path, overrides)
}

/// Load source text (`name` for messages, `base` the directory imports are
/// found from) with the standard distribution's packages.
pub fn load_str(src: &str, name: &str, base: &std::path::Path) -> Result<lang::Doc, InkError> {
    load_str_with(src, name, base, &[])
}

pub fn load_str_with(src: &str, name: &str, base: &std::path::Path, overrides: &[(String, lang::Given)]) -> Result<lang::Doc, InkError> {
    lang::Doc::load_str_in(registry(), src, name, base, overrides)
}
