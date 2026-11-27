//! INK_PATH: an import not found beside its file is read from the first of
//! its directories holding it, and named by its path there.

use crate::common::*;
use crate::mcp::{num, par};
use serde_json::json;

#[test]
fn an_import_not_beside_its_file_is_read_from_ink_path() {
    let s = Scratch::new("ink-path");
    s.model("lib/parts/m/block.ink", "(import \"std\")\n(import \"shared.ink\")\n(body block (shape (box w 1 1)))\n");
    s.model("lib/parts/m/shared.ink", "(param w 3)\n");
    s.model("own/parts/m/block.ink", "(import \"std\")\n(body block (shape (box 5 1 1)))\n");
    let away = s.model("away/a.ink", "(import \"parts/m/block.ink\" as b)\n(assembly root (instance x b.block (at (move 0 0 20))))\n");
    let own = s.model("own/a.ink", "(import \"parts/m/block.ink\" as b)\n(assembly root (instance x b.block))\n");
    let typo = s.model("away/typo.ink", "(import \"parts/m/blok.ink\" as b)\n");
    let path = std::env::join_paths([s.path("nowhere"), s.path("lib")]).unwrap();
    let r = par(5, |i| match i {
        0 => s.ink(["bounds", &away]).env("INK_PATH", &path).ok(),
        1 => s.ink(["check", &away]).env("INK_PATH", &path).ok(),
        2 => s.ink(["bounds", &own]).env("INK_PATH", &path).ok(),
        3 => s.ink(["check", &typo]).env("INK_PATH", &path).err(),
        _ => s.ink(["check", &away]).env_remove("INK_PATH").err(),
    });
    assert!((num(&r[0]["union"]["max"][0]) - 1.5).abs() < 0.1, "a model anywhere imports from INK_PATH: {}", pretty(&r[0]));
    let files: Vec<_> = r[1]["imports"].as_array().unwrap().iter().map(|i| i["file"].clone()).collect();
    assert!(files == [json!("parts/m/block.ink"), json!("parts/m/shared.ink")], "a file in INK_PATH is named by its path there: {}", pretty(&r[1]["imports"]));
    assert!((num(&r[2]["union"]["max"][0]) - 2.5).abs() < 0.1, "a file beside the importer comes before INK_PATH: {}", pretty(&r[2]));
    assert!(r[3]["error"]["hint"] == "did you mean \"parts/m/block.ink\"?", "a typo of a file in INK_PATH is named: {}", pretty(&r[3]));
    let e = &r[4]["error"];
    assert!(
        e["message"].as_str().is_some_and(|m| m.starts_with("cannot read ") && m.contains("parts/m/block.ink"))
            && e["hint"] == "an import is read beside its file, else from INK_PATH's directories (none: INK_PATH is unset)",
        "without INK_PATH, an import not beside its file says where imports are read from: {}",
        pretty(e)
    );
}
