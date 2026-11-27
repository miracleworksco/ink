//! ink with ext's packages (sim, dynamics, topology): the kernel (crates/ink, every
//! path of it re-exported) with a registry of its packages and those of the
//! crates this one depends on, and loading through it. `use ink_ext as ink`
//! reads as the kernel's facade.

pub use ink::*;
/// The solver tools `ink require` asks here, after the kernel's own (solvers.rs).
mod solvers;
/// The dynamics package (ink-dynamics).
pub use ink_dynamics as dynamics;
/// The simulation package (ink-sim).
pub use ink_sim as sim;
/// The topology package (ink-topology).
pub use ink_topology as topology;
pub use solvers::solvers;
/// What this release holds back, refused by name (withheld.rs).
pub mod withheld;

/// The packages of the crates this one depends on, as each declares itself
/// (build.rs): a package joins by its dependency line in Cargo.toml.
pub const PACKAGES: &[&'static dyn lang::Package] = include!(concat!(env!("OUT_DIR"), "/packages.rs"));

/// The kernel's packages, then this distribution's, then those it holds
/// back (withheld.rs), each refused by name.
pub fn registry() -> &'static lang::Registry {
    static REGISTRY: std::sync::OnceLock<lang::Registry> = std::sync::OnceLock::new();
    REGISTRY.get_or_init(|| registry_with(&[]))
}

/// The registry of an `ink` built on this one (the catalog's, the hosted
/// one, research's): the kernel's packages, this distribution's, `added`,
/// then those this release holds back that `added` does not add, each
/// refused by name.
pub fn registry_with(added: &'static [&'static dyn lang::Package]) -> lang::Registry {
    let adds = |name: &str| added.iter().any(|p| p.name() == name);
    let withheld = withheld::WITHHELD.iter().filter(|w| !adds(w.name)).map(|w| *w as &'static dyn lang::Package);
    lang::Registry::new(Vec::leak(ink::PACKAGES.iter().chain(PACKAGES).chain(added).copied().chain(withheld).collect()))
}

/// The language reference, `ink spec`: the kernel's (docs/IR.md), then
/// each package's (build.rs).
pub const SPEC: &str = include_str!(concat!(env!("OUT_DIR"), "/SPEC.md"));

/// Load a file with this distribution's packages.
pub fn load_file(path: &std::path::Path) -> Result<lang::Doc, InkError> {
    load_file_with(path, &[])
}

/// Load a file with this distribution's packages and parameter overrides.
pub fn load_file_with(path: &std::path::Path, overrides: &[(String, lang::Given)]) -> Result<lang::Doc, InkError> {
    lang::Doc::load_file_in(registry(), path, overrides)
}

/// Load a file with this distribution's packages and parameter overrides,
/// an import not found beside its file read from `dirs` before
/// `$INK_PATH`'s (lang::Doc::load_file_importing): the kernel's
/// `load_file_importing` through this registry, so a study, a body or a
/// topology form loads.
pub fn load_file_importing(path: &std::path::Path, overrides: &[(String, lang::Given)], dirs: &[std::path::PathBuf]) -> Result<lang::Doc, InkError> {
    lang::Doc::load_file_importing(registry(), path, overrides, dirs)
}

/// Load source text (`name` for messages, `base` the directory imports are
/// found from) with this distribution's packages.
pub fn load_str(src: &str, name: &str, base: &std::path::Path) -> Result<lang::Doc, InkError> {
    load_str_with(src, name, base, &[])
}

pub fn load_str_with(src: &str, name: &str, base: &std::path::Path, overrides: &[(String, lang::Given)]) -> Result<lang::Doc, InkError> {
    lang::Doc::load_str_in(registry(), src, name, base, overrides)
}
