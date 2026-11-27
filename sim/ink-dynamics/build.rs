// The package this crate is, for a distribution depending on it
// (ink-ext/build.rs): its value and its part of the language
// reference. Cargo hands a build script's metadata only to those of the
// crates depending on it, through `links` (Cargo.toml).
fn main() {
    println!("cargo::metadata=package=ink_dynamics::Dynamics");
    println!("cargo::metadata=spec={}/../docs/DYNAMICS.md", env!("CARGO_MANIFEST_DIR"));
}
