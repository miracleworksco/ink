//! `ink spec` and `ink std`: the language reference (docs/IR.md) and the
//! standard library (std/std.ink, std/materials.ink and std's libraries):
//! their contents (by default, or `--contents`), the part a name finds
//! (`--section NAME`), or whole (`--all`). Either whole is longer than an
//! agent's context takes in one read.
//!
//! A name finds, one rule for both and case aside: the headings whose
//! title or path (`Solving constraints/Over ranged params`) it is, std's sections
//! whose banner it is, the forms it names (the reference's entries for
//! one, `polygon`, and std's definitions); else those whose title,
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
fn spec_find(name: &str, exact: bool) -> Result<Option<String>, InkError> {
    find_in(crate::reference(), name, exact)
}

/// `spec_find` in the reference `text`: the headings whose path or title
/// is `name` (`exact`), else whose title holds it.
fn find_in(text: &str, name: &str, exact: bool) -> Result<Option<String>, InkError> {
    let lines: Vec<&str> = text.lines().collect();
    let heads = heads(&lines);
    let title = |h: &(usize, usize, String)| h.2.rsplit('/').next().unwrap().to_string();
    let mut hits: Vec<&(usize, usize, String)> = match exact {
        true => heads.iter().filter(|h| named(&h.2, name, true) || named(&title(h), name, true)).collect(),
        false => heads.iter().filter(|h| named(&title(h), name, false)).collect(),
    };
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

/// The std a name finds (None: none): the section whose banner it is
/// (`exact`, else holds it); else the definitions it names (else whose
/// names hold it), each with its comments.
fn std_find(name: &str, exact: bool) -> Option<String> {
    let text = std_text();
    let blocks: Vec<&str> = text.split("\n\n").collect();
    if let Some(i) = blocks.iter().position(|b| banner(b).is_some_and(|t| named(&t, name, exact))) {
        return Some(sections(&blocks, i));
    }
    let hits: Vec<&str> = blocks.iter().copied().filter(|b| defs(b).iter().any(|(n, _)| named(n, name, exact))).collect();
    (!hits.is_empty()).then(|| hits.join("\n\n").trim_start().to_string() + "\n")
}

/// The reference's entries for forms, each with the path of the heading
/// it is under: a row of a table (in a plain ``` block, a line opening
/// with a form and the indented lines after it), or a list item or
/// paragraph opening with one in backticks (- `(polygon pts)` ...), to
/// its end. Each with whether it is a row, and the names of the forms it
/// opens with.
fn entries(text: &str) -> Vec<(bool, String, Vec<String>, String)> {
    let lines: Vec<&str> = text.lines().collect();
    let heads = heads(&lines);
    let path = |i: usize| heads.iter().rev().find(|h| h.0 < i).map_or(String::new(), |h| h.2.clone());
    let indent = |l: &str| l.len() - l.trim_start().len();
    let mut out = Vec::new();
    // In a code block, and whether it is plain (a table, not lisp).
    let (mut block, mut table, mut i) = (false, false, 0);
    while i < lines.len() {
        let l = lines[i];
        if l.starts_with("```") {
            (block, table) = (!block, !block && l.trim() == "```");
            i += 1;
            continue;
        }
        if block && !table {
            i += 1;
            continue;
        }
        let item = l.trim_start().strip_prefix("- ");
        let after_blank = i == 0 || lines[i - 1].trim().is_empty() || level(lines[i - 1]).is_some();
        let first = match table {
            true => Some(l),
            false => item.or(after_blank.then_some(l)).and_then(|t| t.strip_prefix('`')),
        };
        let Some(first) = first.filter(|f| f.starts_with('(')) else {
            i += 1;
            continue;
        };
        // Its lines: a row's indented ones; an item's indented past its
        // dash, through blank lines between them (its later paragraphs); a
        // paragraph's to a blank line or a list.
        let more = |m: &str| match (table, item) {
            (true, _) => m.starts_with(' '),
            (false, Some(_)) => indent(m) > indent(l),
            (false, None) => !m.trim_start().starts_with("- "),
        };
        let text = |m: &str| level(m).is_none() && !m.starts_with("```") && more(m);
        let mut end = i + 1;
        loop {
            match lines.get(end) {
                Some(m) if !m.trim().is_empty() && text(m) => end += 1,
                Some(m) if m.trim().is_empty() && item.is_some() && !table => {
                    // Blank lines go on the item where it goes on after them.
                    let next = (end..lines.len()).find(|&j| !lines[j].trim().is_empty());
                    match next.filter(|&j| text(lines[j])) {
                        Some(j) => end = j,
                        None => break,
                    }
                }
                _ => break,
            }
        }
        // A row may open with several forms ((sin a) (cos a)); prose
        // with one, the rest its text.
        let names = opened(first).into_iter().take(if table { usize::MAX } else { 1 }).collect();
        out.push((table, path(i), names, lines[i..end].join("\n")));
        i = end;
    }
    out
}

/// The names of the forms a line opens with: `(a ...) (b ...)` gives a
/// and b; what follows them is the entry's text.
fn opened(line: &str) -> Vec<String> {
    let mut names = Vec::new();
    let mut rest = line.trim_start();
    while let Some(r) = rest.strip_prefix('(') {
        names.push(r.split([' ', ')']).next().unwrap_or_default().to_string());
        let (mut depth, mut at) = (1, r.len());
        for (k, c) in r.char_indices() {
            depth += match c {
                '(' => 1,
                ')' => -1,
                _ => 0,
            };
            if depth == 0 {
                at = k + 1;
                break;
            }
        }
        rest = r[at..].trim_start().trim_start_matches('`').trim_start();
    }
    names
}

/// A form's signature as the reference's first table row opening with it
/// gives it: `(polygon pts [drafts z0 z1])`.
fn signature(text: &str, name: &str) -> Option<String> {
    entries(text).into_iter().filter(|e| e.0).find_map(|(_, _, names, row)| {
        let k = names.iter().position(|n| n == name)?;
        // The k-th form the row opens with, to its closing parenthesis.
        let mut depth = 0;
        let (mut forms, mut from) = (0, 0);
        for (i, c) in row.char_indices() {
            match c {
                '(' if depth == 0 => {
                    from = i;
                    depth = 1;
                }
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    if depth == 0 {
                        if forms == k {
                            return Some(row[from..=i].to_string());
                        }
                        forms += 1;
                    }
                }
                _ => {}
            }
        }
        None
    })
}

/// What the form `name` is (None: no entry or definition is of it): the
/// reference's entries for it (`entries`), each under its heading's path,
/// then std's definitions of it.
fn defined(name: &str) -> Option<String> {
    let found: Vec<(bool, String, Vec<String>, String)> = entries(crate::reference()).into_iter().filter(|e| e.2.iter().any(|n| named(n, name, true))).collect();
    if found.is_empty() {
        return None;
    }
    let mut out = format!("The reference's entries for ({name} ...), each under the section it is in (read one whole with --section PATH):\n");
    let mut under = "";
    for (_, path, _, text) in &found {
        if path != under {
            out += &format!("\nspec \"{path}\":\n");
            under = path;
        }
        out += &format!("{text}\n");
    }
    if let Some(s) = std_find(name, true).filter(|_| std_text().split("\n\n").any(|b| defs(b).iter().any(|(n, _)| named(n, name, true)))) {
        out += &format!("\nstd:\n{s}");
    }
    // Where the reference writes it: an entry may be one use of the form
    // ((attr carried ...)) and the section writing it most its general
    // account (Bodies, attr's), given whole where no section writing it as
    // often has an entry.
    let written = written_in(crate::reference(), name);
    let most = written.first().map_or(0, |w| w.1);
    if let Some((top, _, own)) = written.first().filter(|_| !written.iter().any(|w| w.1 == most && found.iter().any(|e| e.1 == w.0))) {
        out += &format!("\nspec \"{top}\" (written most here, its own text):\n{own}");
    }
    if written.len() > 1 {
        let each: Vec<String> = written.iter().take(6).map(|(p, k, _)| format!("\"{p}\" ({k})")).collect();
        out += &format!("\nThe sections writing ({name} ...), most first: {}\n", each.join(", "));
    }
    Some(out)
}

/// The reference's sections whose own text (not their subsections') writes
/// the form `name`, `(name ` or `(name)`, with how many times and that
/// text, most first.
fn written_in(text: &str, name: &str) -> Vec<(String, usize, String)> {
    let lines: Vec<&str> = text.lines().collect();
    let (open, bare) = (format!("({name} "), format!("({name})"));
    let mut found: Vec<(String, usize, String)> = heads(&lines)
        .into_iter()
        .map(|h| {
            let own: String = lines[h.0 + 1..own_end(&lines, h.0)].iter().map(|l| format!("{l}\n")).collect();
            (h.2, own.matches(&open).count() + own.matches(&bare).count(), own)
        })
        .filter(|w| w.1 > 0)
        .collect();
    found.sort_by_key(|w| std::cmp::Reverse(w.1));
    found
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
    for exact in [true, false] {
        if let Some(s) = spec_find(name, exact)? {
            return Ok(s);
        }
        if let Some(s) = defined(name).filter(|_| exact) {
            return Ok(s);
        }
        if let Some(s) = std_find(name, exact) {
            return Ok(format!("From std (`ink std --section`):\n{s}"));
        }
    }
    uses(name).ok_or_else(|| InkError { hint: Some("ink spec --contents lists the headings".into()), ..InkError::msg(format!("no section of the reference is called '{name}'")) })
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
        out.push_str(&built_in());
        return Ok(out);
    };
    for exact in [true, false] {
        if let Some(s) = std_find(name, exact) {
            return Ok(s);
        }
        if let Some(s) = defined(name).filter(|_| exact) {
            return Ok(s);
        }
        if let Some(s) = spec_find(name, exact)? {
            return Ok(format!("From the reference (`ink spec --section`):\n{s}"));
        }
    }
    uses(name).ok_or_else(|| InkError {
        hint: Some("ink std --contents lists its sections and definitions".into()),
        ..InkError::msg(format!("std has no section or definition called '{name}'"))
    })
}

/// std's contents' last section: the forms built into the language, which
/// std's definitions are written in, each by its signature where a table
/// of the reference gives one. `ink std --section NAME` (or spec's) reads
/// one's entries.
fn built_in() -> String {
    let reference = crate::reference();
    let mut out = "## built in (not std's: the language's own; --section NAME reads one's entries)\n".to_string();
    let mut bare = Vec::new();
    for w in ink::lang::builtin_words() {
        match signature(reference, w) {
            Some(s) => out += &format!("{s}\n"),
            None => bare.push(w),
        }
    }
    if !bare.is_empty() {
        out += &format!("also {}\n", bare.join(" "));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A form's one use found as an entry, (attr carried ...), comes with
    /// the section writing it most, Bodies: what any attribute is, in spec
    /// and std alike; one whose entry is in a section writing it as often
    /// (fn's, Top-level forms) needs no other.
    #[test]
    fn a_forms_lookup_includes_the_section_writing_it_most() {
        let args = |a: &[&str]| crate::Args::parse(&a.iter().map(|s| s.to_string()).collect::<Vec<_>>());
        for attr in [spec(&args(&["spec", "--section", "attr"])).unwrap(), std_lib(&args(&["std", "--section", "attr"])).unwrap()] {
            assert!(attr.contains("(attr carried") && attr.contains("spec \"Top-level forms/Bodies\"") && attr.contains("An attribute has any name"), "{attr}");
        }
        assert!(!spec(&args(&["spec", "--section", "fn"])).unwrap().contains("written most here"));
    }

    /// A name finds a heading by its title or path, else by what its
    /// title holds; std's sections and definitions the same way; where
    /// neither has it, where it is used.
    #[test]
    fn names_find_sections() {
        let spec = |n: &str| spec_find(n, true).unwrap().or_else(|| spec_find(n, false).unwrap()).unwrap_or_default();
        assert!(spec("assemblies").starts_with("### Assemblies"));
        assert!(spec("Standard library/Paths and power transmission").starts_with("### Paths and power transmission"));
        assert!(spec("Solving constraints/Over ranged params").starts_with("### Over ranged params"));
        assert!(spec("top-level").starts_with("## Top-level forms"));
        let sampled = spec("Sampled fields");
        assert!(sampled.starts_with("### Sampled fields") && sampled.matches("### Sampled fields").count() == 1);
        let std = |n: &str| std_find(n, true).or_else(|| std_find(n, false)).unwrap_or_default();
        assert!(std("primitives").contains("(fn box "));
        assert!(std("petg").contains("(def petg ") && std("aluminum").contains("(def aluminum-6061 "));
        assert!(std("involute").contains("(fn involute-gear-2d "));
        assert!(std_find("place", true).is_none() && std_find("place", false).is_none());
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
        let long = find_in(&doc, "Long", true).unwrap().unwrap();
        assert!(!long.contains("### One") && long.contains(r#"Read each as a section, by its path: "Long/One", "Long/Two"."#), "{}", &long[long.len() - 200..]);
        assert!(find_in(&doc, "Long/One", true).unwrap().unwrap().starts_with("### One"));
        assert!(find_in(&doc, "Short", true).unwrap().unwrap().contains("### Three"));
    }

    /// A form's name finds the reference's entries for it, a table's row
    /// and a list item's, before any name merely holding it: polygon is
    /// the built-in profile with its point format, not std's
    /// htd-groove-polygon. std's contents list the built-in forms, each
    /// by its row's signature.
    #[test]
    fn a_forms_name_finds_its_entries() {
        let doc = "## Fns\n\n```\n(sin a) (cos a)       trig\n(warp e fx fy fz)     e with x:=fx,\n                      y:=fy\n```\n\n```lisp\n(sin 1)\n```\n\n- `(probe x y z)`: a point,\n  as a field.\n- `(other)`: no.\n";
        let rows = entries(doc);
        let names: Vec<(bool, Vec<String>)> = rows.iter().map(|e| (e.0, e.2.clone())).collect();
        let want = |t: bool, n: &[&str]| (t, n.iter().map(|s| s.to_string()).collect::<Vec<_>>());
        assert_eq!(names, vec![want(true, &["sin", "cos"]), want(true, &["warp"]), want(false, &["probe"]), want(false, &["other"])]);
        assert!(rows[1].3.ends_with("y:=fy") && rows[2].3.ends_with("as a field.") && rows.iter().all(|e| e.1 == "Fns"), "{rows:?}");
        // A list item runs on through blank lines to its indented paragraphs.
        let item = entries("## A\n\n- `(f a)`: one.\n  - a sub-item.\n\n  Its second paragraph.\n- `(g)`: two.\n\nProse after.\n");
        assert!(item[0].3.ends_with("Its second paragraph.") && item[1].3 == "- `(g)`: two.", "{item:?}");
        assert_eq!(signature(doc, "cos").as_deref(), Some("(cos a)"));
        assert_eq!(signature(doc, "probe"), None, "a list item gives no signature");
        let found = |a: &[&str]| spec(&crate::Args::parse(&a.iter().map(|s| s.to_string()).collect::<Vec<_>>())).unwrap();
        let polygon = found(&["spec", "--section", "polygon"]);
        assert!(polygon.contains("(vec u v b r)") && polygon.contains("(polygon pts [drafts z0 z1])") && !polygon.contains("htd"), "{polygon}");
        assert!(polygon.contains("drops in a row add up"), "the item's second paragraph: {polygon}");
        for (form, row) in [("lambda", "(lambda (a b) body)"), ("fn", "(fn name (a b ...) body)"), ("select", "(select c a b)")] {
            assert!(found(&["spec", "--section", form]).contains(row), "{form}");
        }
        let std = std_lib(&crate::Args::parse(&["std".to_string()])).unwrap();
        let built = &std[std.find("## built in").expect("a section of built-in forms")..];
        assert!(built.contains("\n(polygon pts [drafts z0 z1])\n") && built.contains("\n(map f list ...)\n"), "{built}");
        assert!(std_lib(&crate::Args::parse(&["std", "--section", "polygon"].map(String::from))).unwrap().contains("(vec u v b r)"));
    }
}
