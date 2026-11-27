// The distribution holds the packages it depends on: each package crate's
// build script declares its package (`cargo::metadata=package=PATH`) and
// its part of the language reference (`spec=FILE`), which reach this one
// as DEP_<LINKS>_PACKAGE and DEP_<LINKS>_SPEC. Adding a package to ext's
// `ink` is its dependency line in Cargo.toml.
use std::path::PathBuf;

fn main() {
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let mut found: Vec<(String, String, Option<String>)> = std::env::vars()
        .filter_map(|(k, v)| {
            let links = k.strip_prefix("DEP_")?.strip_suffix("_PACKAGE")?.to_string();
            let spec = std::env::var(format!("DEP_{links}_SPEC")).ok();
            Some((links, v, spec))
        })
        .collect();
    found.sort();
    let list: String = found.iter().map(|(_, p, _)| format!("&{p}, ")).collect();
    std::fs::write(out.join("packages.rs"), format!("&[{list}]")).unwrap();
    let ir = dir.join("../../kernel/docs/IR.md");
    let mut spec = std::fs::read_to_string(&ir).unwrap_or_else(|e| panic!("{}: {e}", ir.display()));
    println!("cargo::rerun-if-changed={}", ir.display());
    for f in found.iter().filter_map(|(_, _, s)| s.as_ref()) {
        let text = std::fs::read_to_string(f).unwrap_or_else(|e| panic!("{f}: {e}"));
        spec = format!("{}\n\n{}\n", spec.trim_end(), text.trim_end());
        println!("cargo::rerun-if-changed={f}");
    }
    std::fs::write(out.join("SPEC.md"), spec).unwrap();
}
