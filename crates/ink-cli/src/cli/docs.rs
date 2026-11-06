//! `ink spec` and `ink std`: the language reference (docs/IR.md) and the
//! standard library (std/std.ink, std/materials.ink and std's libraries):
//! their contents (by default, or `--contents`), the part a name finds
//! (`--section NAME`), or whole (`--all`). Either whole is longer than an
//! agent's context takes in one read.
//!
//! A name finds, one rule for both and case aside: the headings whose
//! title or path (`Solving constraints/Over ranged params`) it is, std's sections
//! whose banner it is, std's definitions it names; else those whose title,
//! banner or name holds it (`primitives`, `involute`). The document asked
//! first, then the other; where neither has it, the sections whose text
//! holds it (else each of its words), with the lines.

use crate::Args;
use crate::mcp::MAX_TEXT;
use ink::InkError;
use ink::lang::STD_SRC;

/// std/materials.ink, read as std's last section.
fn std_text() -> String {
    format!("{}\n\n; ---------------------------------------------------------------- materials\n{}", STD_SRC.trim_end(), ink::materials::LIB.trim_end())
}

/// The reference's headings: each one's line, level and path (its own
/// title after those of the headings it is in).
fn heads(lines: &[&str]) -> Vec<(usize, usize, String)> {
    let mut heads: Vec<(usize, usize, String)> = Vec::new();
    for (i, l) in lines.iter().enumerate() {
        let Some(n) = level(l) else { continue };
        let up = heads.iter().rev().find(|h| h.1 < n).map(|h| format!("{}/", h.2)).unwrap_or_default();
        heads.push((i, n, format!("{up}{}", l.trim_start_matches('#').trim())));
    }
    heads
}

fn level(l: &str) -> Option<usize> {
    l.starts_with("##").then(|| l.chars().take_while(|&c| c == '#').count())
}

/// Where a heading's own text ends: at the next heading, of any level.
fn own_end(lines: &[&str], i: usize) -> usize {
    (i + 1..lines.len()).find(|&j| level(lines[j]).is_some()).unwrap_or(lines.len())
}

/// The reference's Tour (its section `## Tour`, to the next): what `ink
/// mcp` tells every agent that connects.
pub(crate) fn tour() -> &'static str {
    let at = crate::reference().find("\n## Tour\n").map_or(0, |i| i + 1);
    let end = crate::reference()[at + 1..].find("\n## ").map_or(crate::reference().len(), |i| at + 1 + i);
    crate::reference()[at..end].trim_end()
}

/// Whether `name` is `title` (`exact`) or in it, case aside.
fn named(title: &str, name: &str, exact: bool) -> bool {
    let (t, n) = (title.to_lowercase(), name.to_lowercase());
    if exact { t == n } else { t.contains(&n) }
}

/// The reference section a name finds (Ok(None): none): the headings whose
/// path or title it is, else whose title holds it. A heading with its own
/// subsections is the heading. A section longer than an MCP answer with
/// its subsections gives its own text and their paths, each read as a
/// section (the CLI's tests keep each within). Several sections: all of
/// them, where they fit an MCP answer, else their paths.
fn spec_find(name: &str) -> Result<Option<String>, InkError> {
    find_in(crate::reference(), name)
}

/// `spec_find` in the reference `text`.
fn find_in(text: &str, name: &str) -> Result<Option<String>, InkError> {
    let lines: Vec<&str> = text.lines().collect();
    let heads = heads(&lines);
    let title = |h: &(usize, usize, String)| h.2.rsplit('/').next().unwrap().to_string();
    let mut hits: Vec<&(usize, usize, String)> = heads.iter().filter(|h| named(&h.2, name, true) || named(&title(h), name, true)).collect();
    if hits.is_empty() {
        hits = heads.iter().filter(|h| named(&title(h), name, false)).collect();
    }
    // A heading and its own subsections: the heading (it holds them).
    let tops: Vec<String> = hits.iter().map(|t| format!("{}/", t.2)).collect();
    hits.retain(|h| !tops.iter().any(|t| h.2.starts_with(t.as_str())));
    let section = |&(i, n, ref path): &(usize, usize, String)| -> String {
        let end = (i + 1..lines.len()).find(|&j| level(lines[j]).is_some_and(|m| m <= n)).unwrap_or(lines.len());
        let text = |to: usize| -> String { lines[i..to].iter().map(|l| format!("{l}\n")).collect() };
        let whole = text(end);
        let subs: Vec<&(usize, usize, String)> = heads.iter().filter(|h| h.0 > i && h.0 < end && h.2.rsplit_once('/').is_some_and(|(up, _)| up == path)).collect();
        match subs.first() {
            Some(first) if whole.len() > MAX_TEXT => {
                let names: Vec<String> = subs.iter().map(|h| format!("\"{}\"", h.2)).collect();
                format!("{}(Its subsections are too long to include. Read each as a section, by its path: {}.)\n", text(first.0), names.join(", "))
            }
            _ => whole,
        }
    };
    let all: String = hits.iter().map(|h| section(h)).collect();
    match hits.len() {
        0 => Ok(None),
        1 => Ok(Some(all)),
        _ if all.len() <= MAX_TEXT => Ok(Some(all)),
        _ => Err(InkError {
            hint: Some(format!("name one by its path: {}", hits.iter().map(|h| h.2.as_str()).collect::<Vec<_>>().join("; "))),
            ..InkError::msg(format!("'{name}' names more than one section"))
        }),
    }
}

/// std's blocks (a definition is a block of lines between blank lines
/// naming it), and a block's banner (`; ---- NAME`) and the forms it
/// defines: its lines opening (fn NAME (ARGS) or (def NAME.
fn banner(b: &str) -> Option<String> {
    b.lines().next().and_then(|l| l.trim_start_matches(';').trim().strip_prefix("---")).map(|r| r.trim_start_matches('-').trim().to_string())
}

fn defs(b: &str) -> Vec<(String, String)> {
    b.lines()
        .filter_map(|l| {
            let rest = l.strip_prefix("(fn ").or_else(|| l.strip_prefix("(def "))?;
            let name = rest.split([' ', ')']).next()?.to_string();
            let sig = match l.starts_with("(fn ") {
                true => l.find(')').map_or(l.to_string(), |i| format!("{})", &l[..=i])),
                false => format!("(def {name})"),
            };
            Some((name, sig))
        })
        .collect()
}

/// The std a name finds (None: none): the section whose banner it is,
/// else the definitions it names, each with its comments; else the section
/// whose banner holds it, else the definitions whose names hold it.
fn std_find(name: &str) -> Option<String> {
    let text = std_text();
    let blocks: Vec<&str> = text.split("\n\n").collect();
    for exact in [true, false] {
        if let Some(i) = blocks.iter().position(|b| banner(b).is_some_and(|t| named(&t, name, exact))) {
            return Some(sections(&blocks, i));
        }
        let hits: Vec<&str> = blocks.iter().copied().filter(|b| defs(b).iter().any(|(n, _)| named(n, name, exact))).collect();
        if !hits.is_empty() {
            return Some(hits.join("\n\n").trim_start().to_string() + "\n");
        }
    }
    None
}

/// std's section starting at block i, to the next banner.
fn sections(blocks: &[&str], i: usize) -> String {
    let end = (i + 1..blocks.len()).find(|&j| banner(blocks[j]).is_some()).unwrap_or(blocks.len());
    blocks[i..end].join("\n\n") + "\n"
}

/// Where a name no title has is used: the sections of both documents
/// whose own text holds it (a reference section's without its
/// subsections, a std section's by its banner), else each of its words,
/// most uses first, each with the first lines holding one. None: no
/// section holds them.
fn uses(name: &str) -> Option<String> {
    let lines: Vec<&str> = crate::reference().lines().collect();
    let mut parts: Vec<(String, Vec<&str>)> = heads(&lines).into_iter().map(|h| (format!("spec \"{}\"", h.2), lines[h.0 + 1..own_end(&lines, h.0)].to_vec())).collect();
    let text = std_text();
    let mut banner_at = "";
    for b in text.split("\n\n") {
        if let Some(t) = banner(b) {
            banner_at = b.lines().next().unwrap();
            parts.push((format!("std \"{t}\""), Vec::new()));
        }
        if !banner_at.is_empty()
            && let Some(p) = parts.last_mut()
        {
            p.1.extend(b.lines());
        }
    }
    let lower = name.to_lowercase();
    let phrase = [lower.clone()];
    let each: Vec<String> = lower.split_whitespace().map(str::to_string).collect();
    for terms in [&phrase[..], &each[..]] {
        let hit = |l: &str| terms.iter().filter(|t| l.to_lowercase().contains(t.as_str())).count();
        let mut found: Vec<(usize, &String, Vec<&str>)> = parts
            .iter()
            .filter(|(_, ls)| terms.iter().all(|t| ls.iter().any(|l| l.to_lowercase().contains(t.as_str()))))
            .map(|(path, ls)| (ls.iter().map(|l| hit(l)).sum(), path, ls.iter().copied().filter(|l| hit(l) > 0).take(3).collect()))
            .collect();
        if found.is_empty() {
            continue;
        }
        found.sort_by_key(|f| std::cmp::Reverse(f.0));
        let mut out = format!("No heading, std section or definition is called '{name}'. The sections using it, most first (read one with its command and --section):\n");
        for (n, path, ls) in found.iter().take(6) {
            out += &format!("{path} ({n} uses):\n");
            ls.iter().for_each(|l| out += &format!("  {}\n", l.trim()));
        }
        return Some(out);
    }
    None
}

/// The reference whole (`--all`), the part a name finds (`--section
/// NAME`: a section of it, else of std, else where its words are used), or
/// its headings (by default).
pub(crate) fn spec(a: &Args) -> Result<String, InkError> {
    if a.get("all").is_some() {
        return Ok(crate::reference().to_string());
    }
    let Some(name) = a.get("section") else {
        let lines: Vec<&str> = crate::reference().lines().collect();
        return Ok(heads(&lines).iter().map(|h| format!("{} {}\n", "#".repeat(h.1), h.2)).collect());
    };
    if let Some(s) = spec_find(name)? {
        return Ok(s);
    }
    std_find(name)
        .map(|s| format!("From std (`ink std --section`):\n{s}"))
        .or_else(|| uses(name))
        .ok_or_else(|| InkError { hint: Some("ink spec --contents lists the headings".into()), ..InkError::msg(format!("no section of the reference is called '{name}'")) })
}

/// The standard library whole (`--all`), the part a name finds
/// (`--section NAME`: std's, else a section of the reference, else where
/// its words are used), or its sections and definitions' signatures (by
/// default).
pub(crate) fn std_lib(a: &Args) -> Result<String, InkError> {
    let text = std_text();
    if a.get("all").is_some() {
        return Ok(text);
    }
    let Some(name) = a.get("section") else {
        let mut out = String::new();
        for b in text.split("\n\n") {
            if let Some(t) = banner(b) {
                out.push_str(&format!("## {t}\n"));
            }
            for (_, sig) in defs(b) {
                out.push_str(&format!("{sig}\n"));
            }
        }
        return Ok(out);
    };
    if let Some(s) = std_find(name) {
        return Ok(s);
    }
    match spec_find(name)? {
        Some(s) => Ok(format!("From the reference (`ink spec --section`):\n{s}")),
        None => uses(name).ok_or_else(|| InkError {
            hint: Some("ink std --contents lists its sections and definitions".into()),
            ..InkError::msg(format!("std has no section or definition called '{name}'"))
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A name finds a heading by its title or path, else by what its
    /// title holds; std's sections and definitions the same way; where
    /// neither has it, where it is used.
    #[test]
    fn names_find_sections() {
        let spec = |n: &str| spec_find(n).unwrap().unwrap_or_default();
        assert!(spec("assemblies").starts_with("### Assemblies"));
        assert!(spec("Standard library/Paths and power transmission").starts_with("### Paths and power transmission"));
        assert!(spec("Solving constraints/Over ranged params").starts_with("### Over ranged params"));
        assert!(spec("top-level").starts_with("## Top-level forms"));
        let sampled = spec("Sampled fields");
        assert!(sampled.starts_with("### Sampled fields") && sampled.matches("### Sampled fields").count() == 1);
        let std = |n: &str| std_find(n).unwrap_or_default();
        assert!(std("primitives").contains("(fn box "));
        assert!(std("petg").contains("(def petg ") && std("aluminum").contains("(def aluminum-6061 "));
        assert!(std("involute").contains("(fn involute-gear-2d "));
        assert!(std_find("place").is_none());
        let used = |n: &str| uses(n).unwrap_or_default();
        assert!(used("interfere").contains("ink interfere"));
        assert!(uses("voronoi").is_none());
        assert!(tour().starts_with("## Tour\n") && !tour().contains("\n## ") && tour().contains("ink assert FILE"));
    }

    /// A section longer than an MCP answer with its subsections gives its
    /// own text and their paths; one that fits gives them whole (in a
    /// reference made to have both, whatever the reference's own sizes).
    #[test]
    fn a_section_too_long_gives_its_subsections_paths() {
        let filler = |n: usize| "words ".repeat(n / 6) + "\n";
        let doc = format!("## Long\n{}### One\n{}### Two\n{}## Short\n{}### Three\n{}", filler(100), filler(MAX_TEXT / 2), filler(MAX_TEXT / 2), filler(100), filler(100));
        let long = find_in(&doc, "Long").unwrap().unwrap();
        assert!(!long.contains("### One") && long.contains(r#"Read each as a section, by its path: "Long/One", "Long/Two"."#), "{}", &long[long.len() - 200..]);
        assert!(find_in(&doc, "Long/One").unwrap().unwrap().starts_with("### One"));
        assert!(find_in(&doc, "Short").unwrap().unwrap().contains("### Three"));
    }
}
