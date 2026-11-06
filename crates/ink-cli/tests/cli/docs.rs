//! The docs (README.md, docs/GUIDE.md, docs/guide/walkthrough.md, docs/IR.md) can't go stale: every
//! model and answer they show is run. A block is marked by an HTML comment
//! on the line before it (GitHub doesn't render one):
//!
//!   <!-- file NAME -->    ```lisp: a whole model, written to NAME
//!   <!-- add NAME -->     ```lisp: more of NAME, added to its end
//!   <!-- run -->          ```sh: each `ink ...` line is run; a ```json block
//!                         right after it is what the last one answers
//!   <!-- mcp -->          ```json: an `ink` tool call's arguments over
//!                         `ink mcp`; the ```json block after it, its answer
//!   <!-- not run: WHY --> a block not ```lisp, shown and not run
//!
//! A ```lisp block with no marker loads (`ink check --load`), std imported
//! for it where it doesn't import std itself (IR.md: an excerpt is read with
//! std). One that needs names the text around it gives is fenced
//! ```lisp fragment, and not run. A block with no language and no marker is
//! text (a table of forms), not run: one that begins `(import "std")` is a
//! model missing its ```lisp.
//!
//! A picture shown (`![...](guide/NAME.png)`) is what a command before it
//! rendered (`--out NAME.png`), to within `PICTURE_NOISE` (unlike_pixels);
//! INK_BLESS=1 rewrites them.
//!
//! Every model, once written or added to, loads with no warning and every
//! claim it states holds, or answers as the first run block after it shows
//! it answering (a claim the guide shows unknown, say). A shown answer is a
//! part of the real one: each key shown is there with the value shown, a
//! number to the digits shown (a whole number exactly), a list of numbers
//! whole, a list of objects in order with others between. Commands
//! run in a scratch directory where `examples` is the repository's. Links
//! into the reference (IR.md#SLUG) must name one of its headings.
//!
//! examples/README.md lists every example in examples/, and nothing else, each
//! with a link to the guide's section it goes with.

use crate::common::*;
use serde_json::Value;
use std::collections::BTreeMap;

#[derive(Debug)]
struct Block {
    marker: Option<String>,
    lang: String,
    body: String,
    line: usize,
}

/// The fenced blocks of a Markdown text, each with the marker comment on the
/// last non-blank line before it.
fn blocks(text: &str) -> Vec<Block> {
    let mut out = vec![];
    let mut lines = text.lines().enumerate();
    let mut last: Option<&str> = None;
    while let Some((i, l)) = lines.next() {
        if let Some(lang) = l.strip_prefix("```") {
            let mut body = String::new();
            for (_, b) in lines.by_ref() {
                if b.starts_with("```") {
                    break;
                }
                body += b;
                body.push('\n');
            }
            let marker = last.and_then(|m| m.trim().strip_prefix("<!--")).and_then(|m| m.strip_suffix("-->")).map(|m| m.trim().to_string());
            out.push(Block { marker, lang: lang.trim().to_string(), body, line: i + 1 });
            last = None;
        } else if !l.trim().is_empty() {
            last = Some(l);
        }
    }
    out
}

/// A command line's words, as a shell splits them (quotes, no expansions),
/// up to a `#` comment.
fn words(line: &str) -> Vec<String> {
    let (mut out, mut w, mut quote, mut any) = (vec![], String::new(), None, false);
    for c in line.chars() {
        match (quote, c) {
            (Some(q), c) if c == q => quote = None,
            (Some(_), c) => w.push(c),
            (None, '\'' | '"') => {
                quote = Some(c);
                any = true;
            }
            (None, '#') if w.is_empty() && !any => break,
            (None, c) if c.is_whitespace() => {
                if !w.is_empty() || any {
                    out.push(std::mem::take(&mut w));
                }
                any = false;
            }
            (None, c) => w.push(c),
        }
    }
    if !w.is_empty() || any {
        out.push(w);
    }
    out
}

/// An answer as shown: JSON whose numbers are kept as written (a string
/// "\u{1}LITERAL"), so they are compared to the digits shown.
fn shown(text: &str) -> Value {
    let (mut s, mut chars, mut in_str) = (String::new(), text.chars().peekable(), false);
    while let Some(c) = chars.next() {
        if in_str {
            s.push(c);
            if c == '\\' {
                s.extend(chars.next());
            } else if c == '"' {
                in_str = false;
            }
        } else if c == '"' {
            in_str = true;
            s.push(c);
        } else if c == '-' || c.is_ascii_digit() {
            let mut lit = String::from(c);
            while let Some(&d) = chars.peek().filter(|d| d.is_ascii_digit() || matches!(d, '.' | 'e' | 'E' | '+' | '-')) {
                lit.push(d);
                chars.next();
            }
            s += &format!("\"\\u0001{lit}\"");
        } else {
            s.push(c);
        }
    }
    serde_json::from_str(&s).unwrap_or_else(|e| panic!("a shown answer is not JSON ({e}):\n{text}"))
}

/// Whether the number `got` reads as `lit` to the digits `lit` shows; a
/// whole number shown (no point, no exponent) is that number exactly.
fn same_number(lit: &str, got: f64) -> bool {
    if !lit.contains(['.', 'e', 'E']) {
        return lit.parse::<f64>().is_ok_and(|w| w == got);
    }
    let (mant, exp) = lit.split_once(['e', 'E']).map_or((lit, 0), |(m, e)| (m, e.parse::<i32>().unwrap()));
    let decimals = mant.split_once('.').map_or(0, |(_, d)| d.len() as i32);
    let want: f64 = lit.parse().unwrap();
    (got - want).abs() <= 0.5 * 10f64.powi(exp - decimals) * (1.0 + 1e-9)
}

/// Where the answer `got` differs from the part of it shown, `want`.
fn differs(want: &Value, got: &Value, at: &str) -> Vec<String> {
    match (want, got) {
        (Value::String(w), _) if w.starts_with('\u{1}') => match got.as_f64() {
            Some(g) if same_number(&w[1..], g) => vec![],
            _ => vec![format!("{at}: shown {}, got {got}", &w[1..])],
        },
        (Value::Object(w), Value::Object(g)) => w
            .iter()
            .flat_map(|(k, v)| match g.get(k) {
                Some(gv) => differs(v, gv, &format!("{at}.{k}")),
                None => vec![format!("{at}.{k}: shown, but not in the answer")],
            })
            .collect(),
        (Value::Array(w), Value::Array(g)) if w.iter().any(Value::is_object) => {
            let mut rest = g.iter();
            let mut out = vec![];
            for (i, wv) in w.iter().enumerate() {
                if !rest.by_ref().any(|gv| differs(wv, gv, "").is_empty()) {
                    out.push(format!("{at}[{i}]: shown, but no later item of the answer is it: {}", wv.to_string().replace("\\u0001", "")));
                    break;
                }
            }
            out
        }
        (Value::Array(w), Value::Array(g)) if w.len() == g.len() => w.iter().zip(g).enumerate().flat_map(|(i, (a, b))| differs(a, b, &format!("{at}[{i}]"))).collect(),
        _ if want == got => vec![],
        _ => vec![format!("{at}: shown {want}, got {got}")],
    }
}

/// One `ink mcp` session in `dir` making the tool call `args`: its answer.
fn mcp_call(dir: &Scratch, args: &Value) -> Value {
    let init = serde_json::json!({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "docs", "version": "0"}}});
    let done = serde_json::json!({"jsonrpc": "2.0", "method": "notifications/initialized"});
    let call = serde_json::json!({"jsonrpc": "2.0", "id": 2, "method": "tools/call", "params": {"name": "ink", "arguments": args}});
    let run = dir.ink(["mcp"]).cwd(dir).stdin(&format!("{init}\n{done}\n{call}\n"));
    let reply = run
        .stdout
        .lines()
        .filter_map(|l| serde_json::from_str::<Value>(l).ok())
        .find(|r| r["id"] == 2)
        .unwrap_or_else(|| panic!("ink mcp gave no answer to the call:\n{}", run.stderr));
    let text = reply["result"]["content"][0]["text"].as_str().unwrap_or_else(|| panic!("the call's answer has no text: {reply}"));
    serde_json::from_str(text).unwrap()
}

/// What a snippet asks: a model, as written so far, that must load clean
/// and hold its claims; an example that must load; command lines to run;
/// or an `ink mcp` tool call.
enum Job {
    /// A model, and the claims the next run shows not holding: name, status.
    Model(String, BTreeMap<String, String>),
    /// An example's text, std imported where it doesn't import it.
    Load(String),
    Run(Vec<String>),
    Mcp(Value),
}

/// A snippet to check, with the files written before it.
struct Step {
    at: String,
    job: Job,
    /// The answer shown after it: its line, and its text.
    answer: Option<(usize, String)>,
    files: BTreeMap<String, String>,
    /// The steps it waits for: a command, every model before it; a model,
    /// its own last version. What ink caches for them (a study's solve, a
    /// claim) is then reused, not computed twice at once.
    after: Vec<usize>,
}

/// The snippets of a Markdown text to check, and the problems with its
/// markers.
fn steps(rel: &str, text: &str) -> (Vec<Step>, Vec<String>) {
    let mut files = BTreeMap::<String, String>::new();
    let (mut out, mut problems) = (Vec::<Step>::new(), vec![]);
    let (mut models, mut last_of) = (vec![], BTreeMap::<String, usize>::new());
    let bs = blocks(text);
    let mut i = 0;
    while i < bs.len() {
        let b = &bs[i];
        let at = format!("{rel}:{}", b.line);
        let marker = b.marker.clone().unwrap_or_default();
        let (kind, arg) = marker.split_once(' ').map_or((marker.as_str(), ""), |(k, a)| (k, a.trim()));
        // The ```json block after a run or mcp block is its answer.
        let answer = bs.get(i + 1).filter(|n| n.lang == "json" && n.marker.is_none() && matches!(kind, "run" | "mcp"));
        let shown = answer.map(|a| (a.line, a.body.clone()));
        let job = match (kind, b.lang.as_str()) {
            ("", "lisp fragment") => None,
            ("", "lisp") => Some(Job::Load(if b.body.starts_with("(import \"std\")") { b.body.clone() } else { format!("(import \"std\")\n{}", b.body) })),
            ("", "") => {
                if b.body.starts_with("(import \"std\")") {
                    problems.push(format!("{at}: a model fenced with no language: fence it ```lisp"));
                }
                None
            }
            ("not", lang) if arg.starts_with("run:") && !lang.starts_with("lisp") => None,
            ("file" | "add", "lisp") => {
                let model = if kind == "file" {
                    b.body.clone()
                } else {
                    files.get(arg).cloned().map(|m| m + &b.body).unwrap_or_else(|| panic!("{at}: add to {arg}, which no block before writes"))
                };
                files.insert(arg.to_string(), model);
                // The first run after it, where it shows claims answering
                // otherwise than holds.
                let next = bs[i + 1..].iter().position(|n| n.marker.as_deref().is_some_and(|m| m.starts_with("run"))).map(|k| i + 1 + k);
                let mut shown_not = BTreeMap::new();
                if let Some(ans) = next.and_then(|k| bs.get(k + 1)).filter(|a| a.lang == "json" && a.marker.is_none()) {
                    for r in self::shown(&ans.body)["asserts"].as_array().into_iter().flatten() {
                        if let (Some(n), Some(st)) = (r["name"].as_str(), r["status"].as_str())
                            && st != "holds"
                        {
                            shown_not.insert(n.to_string(), st.to_string());
                        }
                    }
                }
                Some(Job::Model(arg.to_string(), shown_not))
            }
            ("run", "sh") if arg.is_empty() => {
                let lines: Vec<String> = b.body.lines().filter(|l| !l.trim().is_empty()).map(String::from).collect();
                for l in lines.iter().filter(|l| words(l).first().map(String::as_str) != Some("ink")) {
                    problems.push(format!("{at}: `{l}` is run, but is not an ink command (mark the block not run)"));
                }
                Some(Job::Run(lines))
            }
            ("mcp", "json") => {
                let args: Value = serde_json::from_str(&b.body).unwrap_or_else(|e| panic!("{at}: the call is not JSON: {e}"));
                if answer.is_none() {
                    problems.push(format!("{at}: an mcp call with no answer shown after it"));
                }
                Some(Job::Mcp(args))
            }
            _ => {
                problems.push(format!("{at}: a ```{} block marked {marker:?}: mark it file, add, run or mcp, or (not ```lisp) not run", b.lang));
                None
            }
        };
        if let Some(job) = job {
            let n = out.len();
            let after = match &job {
                Job::Model(name, _) => last_of.insert(name.clone(), n).into_iter().collect(),
                Job::Load(_) => vec![],
                _ => models.clone(),
            };
            if matches!(job, Job::Model(..)) {
                models.push(n);
            }
            out.push(Step { at, job, answer: shown, files: files.clone(), after });
        }
        i += if answer.is_some() { 2 } else { 1 };
    }
    (out, problems)
}

/// What is wrong with one step: run in its own scratch directory, holding
/// the files written before it, where `examples` is the repository's.
fn check_step(rel: &str, s: &Step, dir: &Scratch) -> Vec<String> {
    #[cfg(unix)]
    std::os::unix::fs::symlink(root("examples"), dir.join("examples")).unwrap();
    for (name, text) in &s.files {
        dir.model(name, text);
    }
    let at = &s.at;
    let mut problems = vec![];
    let got = match &s.job {
        Job::Model(name, shown_not) => {
            let path = dir.path(name);
            let c = dir.ink(["check", &path]).cwd(dir).run().json();
            if c["ok"] != true {
                problems.push(format!("{at}: {name} does not load: {}", c["error"]));
            } else if c["warnings"].as_array().is_some_and(|w| !w.is_empty()) {
                problems.push(format!("{at}: {name} warns: {}", c["warnings"]));
            }
            let a = dir.ink(["assert", &path]).cwd(dir).run().json();
            let bad: Vec<String> = a["asserts"]
                .as_array()
                .into_iter()
                .flatten()
                .filter(|r| r["status"] != "holds" && r["name"].as_str().and_then(|n| shown_not.get(n)).is_none_or(|st| r["status"] != st.as_str()))
                .map(|r| format!("{} {}", r["name"], r["status"]))
                .collect();
            if a["ok"] != true || !bad.is_empty() {
                problems.push(format!("{at}: {name}'s claims do not all hold: {bad:?} {}", a["error"]));
            }
            None
        }
        Job::Load(text) => {
            let c = dir.ink(["check", &dir.model("example.ink", text), "--load"]).cwd(dir).run().json();
            if c["ok"] != true {
                problems.push(format!("{at}: the example does not load (fence it ```lisp fragment if it needs names the text gives): {}", c["error"]));
            }
            None
        }
        Job::Run(lines) => {
            let mut last = None;
            for line in lines {
                let w = words(line);
                if w.first().map(String::as_str) != Some("ink") {
                    continue;
                }
                let r = dir.ink(&w[1..]).cwd(dir).run();
                match serde_json::from_str::<Value>(&r.stdout) {
                    Ok(d) if d["ok"] == true => last = Some(d),
                    _ => problems.push(format!("{at}: `{line}` gave no ok answer:\n{}{}", r.stdout, r.stderr)),
                }
            }
            last
        }
        Job::Mcp(args) => Some(mcp_call(dir, args)),
    };
    if let (Some((line, want)), Some(got)) = (&s.answer, &got) {
        problems.extend(differs(&shown(want), got, "").into_iter().map(|d| format!("{rel}:{line}: {d}")));
    }
    problems
}

/// What is wrong with the snippets of the Markdown file at `rel`. Each runs
/// as soon as the steps it waits for are done, in parallel: the guide's are
/// a chain of solves, minutes run one by one.
fn check_doc(rel: &str) -> Vec<String> {
    use std::sync::{Condvar, Mutex};
    let text = std::fs::read_to_string(root(rel)).unwrap();
    let (steps, mut problems) = steps(rel, &text);
    let done = (Mutex::new(vec![false; steps.len()]), Condvar::new());
    /// Marks its step done however it ends, a panic included, so no step
    /// waits for it forever.
    struct Done<'a>(&'a (Mutex<Vec<bool>>, Condvar), usize);
    impl Drop for Done<'_> {
        fn drop(&mut self) {
            let mut d = self.0.0.lock().unwrap_or_else(|e| e.into_inner());
            d[self.1] = true;
            self.0.1.notify_all();
        }
    }
    let results: Vec<(Vec<String>, Scratch)> = std::thread::scope(|scope| {
        let handles: Vec<_> = steps
            .iter()
            .enumerate()
            .map(|(n, s)| {
                let done = &done;
                scope.spawn(move || {
                    let _mark = Done(done, n);
                    let mut d = done.0.lock().unwrap_or_else(|e| e.into_inner());
                    while !s.after.iter().all(|&a| d[a]) {
                        d = done.1.wait(d).unwrap_or_else(|e| e.into_inner());
                    }
                    drop(d);
                    let dir = Scratch::new("docs");
                    (check_step(rel, s, &dir), dir)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e))).collect()
    });
    let dirs: Vec<&std::path::Path> = results.iter().map(|(_, d)| d.as_ref()).collect();
    problems.extend(results.iter().flat_map(|(p, _)| p.clone()));
    problems.extend(reference_links(rel, &text));
    problems.extend(images(rel, &text, &dirs));
    problems
}

/// Of two renders' pixels, how many may differ by more than `PIXEL_LEVELS`
/// on a channel: 0.04% of them. Measured on the guide's bracket at 640 px:
/// the same build's render on Linux x86 (two machines, one Zen 5),
/// Linux arm64 and macOS arm64 differed in 0 pixels (shaded facet by
/// facet; smooth shading, which adds no arithmetic that differs by
/// platform, measured on Linux x86 only); an honest change differs in many
/// more (smooth shading): its fillet 5 -> 4 mm in 815 (0.20%), a hole 5.5
/// -> 5.3 mm in 691 (0.17%), the lug 1 mm narrower in 5089 (1.2%).
const PICTURE_NOISE: f64 = 0.0004;
/// A channel difference under this is shading, not a change of shape.
const PIXEL_LEVELS: u8 = 24;

/// How many pixels of two PNGs ink wrote differ by more than PIXEL_LEVELS
/// on a channel, and how many there are; None where their sizes differ.
fn unlike_pixels(a: &str, b: &str) -> Option<(usize, usize)> {
    let ((wa, ha, ra), (wb, hb, rb)) = (crate::imports::png_rows(a), crate::imports::png_rows(b));
    if (wa, ha) != (wb, hb) {
        return None;
    }
    let n = ra.iter().zip(&rb).map(|(x, y)| x[1..].chunks(3).zip(y[1..].chunks(3)).filter(|(p, q)| p.iter().zip(*q).any(|(c, d)| c.abs_diff(*d) > PIXEL_LEVELS)).count()).sum();
    Some((n, wa * ha))
}

/// Whether two renders show the same thing, to within PICTURE_NOISE.
fn alike(a: &str, b: &str) -> Result<(), String> {
    match unlike_pixels(a, b) {
        None => Err("its size differs".into()),
        Some((n, all)) if n as f64 > PICTURE_NOISE * all as f64 => Err(format!("{n} of its {all} pixels differ (more than {} allowed)", (PICTURE_NOISE * all as f64) as usize)),
        Some(_) => Ok(()),
    }
}

/// The pictures a document shows (`![...](guide/NAME.png)`, beside it):
/// each must be what one of its commands rendered (`--out NAME.png`), to
/// within PICTURE_NOISE. INK_BLESS=1 writes them instead.
fn images(rel: &str, text: &str, dirs: &[&std::path::Path]) -> Vec<String> {
    let base = std::path::Path::new(&root(rel)).parent().unwrap().to_path_buf();
    let mut out = vec![];
    for shown in text.split("](").skip(1).filter_map(|s| s.split(')').next()).filter(|s| s.ends_with(".png")) {
        let name = std::path::Path::new(shown).file_name().unwrap();
        // The last step to write it rendered it.
        let Some(made_at) = dirs.iter().rev().map(|d| d.join(name)).find(|p| p.exists()) else {
            out.push(format!("{rel}: shows {shown}, which none of its commands renders (--out {})", name.to_string_lossy()));
            continue;
        };
        let made = std::fs::read(&made_at).unwrap();
        let path = base.join(shown);
        if std::env::var_os("INK_BLESS").is_some() {
            std::fs::write(&path, &made).unwrap();
        } else if !path.exists() {
            out.push(format!("{rel}: {shown} is not there: rerun with INK_BLESS=1 to write it, and look at it"));
        } else if let Err(why) = alike(path.to_str().unwrap(), made_at.to_str().unwrap()) {
            out.push(format!("{rel}: {shown} is not what its command renders now ({why}): rerun with INK_BLESS=1 to rewrite it, and look at it"));
        }
    }
    out
}

/// GitHub's anchor for a heading.
fn slug(heading: &str) -> String {
    heading.trim().to_lowercase().chars().filter(|c| c.is_alphanumeric() || *c == ' ' || *c == '-' || *c == '_').map(|c| if c == ' ' { '-' } else { c }).collect()
}

fn anchors(rel: &str) -> Vec<String> {
    std::fs::read_to_string(root(rel)).unwrap().lines().filter_map(|l| l.strip_prefix('#')).map(|h| slug(h.trim_start_matches('#'))).collect()
}

/// Links into the reference that name no heading of it.
fn reference_links(rel: &str, text: &str) -> Vec<String> {
    let ir = anchors("docs/IR.md");
    text.split("IR.md#")
        .skip(1)
        .map(|s| s.split(')').next().unwrap_or(""))
        .filter(|a| !ir.iter().any(|h| h == a))
        .map(|a| format!("{rel}: IR.md#{a} names no heading of the reference"))
        .collect()
}

#[test]
fn the_docs_show_what_ink_answers() {
    let docs = ["README.md", "docs/GUIDE.md", "docs/guide/walkthrough.md", "docs/IR.md"];
    let p: Vec<String> = std::thread::scope(|s| docs.map(|d| s.spawn(move || check_doc(d))).map(|h| h.join().unwrap()).concat());
    assert!(p.is_empty(), "the docs are stale:\n{}", p.join("\n"));
}

/// The index of examples/: one line per file, linking it and the guide's
/// section it goes with.
#[test]
fn the_examples_index_lists_every_example() {
    let index = std::fs::read_to_string(root("examples/README.md")).unwrap();
    let guide = anchors("docs/guide/walkthrough.md");
    let mut listed = vec![];
    let mut problems = vec![];
    for l in index.lines().filter(|l| l.starts_with("- ")) {
        let file = l.split("](").nth(1).and_then(|s| s.split(')').next()).unwrap_or("");
        listed.push(file.to_string());
        match l.split("walkthrough.md#").nth(1).and_then(|s| s.split(')').next()) {
            Some(a) if guide.iter().any(|h| h == a) => {}
            Some(a) => problems.push(format!("{file}: walkthrough.md#{a} names no section of the walkthrough")),
            None => problems.push(format!("{file}: no link to the guide's section it goes with")),
        }
    }
    let mut files: Vec<String> = std::fs::read_dir(root("examples")).unwrap().filter_map(|e| e.ok()?.file_name().into_string().ok()).filter(|n| n.ends_with(".ink")).collect();
    files.sort();
    problems.extend(files.iter().filter(|f| !listed.contains(f)).map(|f| format!("{f}: in examples/, not in its README")));
    problems.extend(listed.iter().filter(|f| !files.contains(f)).map(|f| format!("{f}: in examples/README.md, not in examples/")));
    assert!(problems.is_empty(), "examples/README.md:\n{}", problems.join("\n"));
}

/// The picture check catches the guide's bracket with its fillet made 4 mm
/// instead of 5: the tolerance is below an honest change. (That it passes
/// the bracket rendered again is the guide's own check of bracket.png.)
#[test]
fn a_picture_is_matched_to_its_shape() {
    let text = std::fs::read_to_string(root("docs/guide/walkthrough.md")).unwrap();
    let bracket = blocks(&text).into_iter().find(|b| b.marker.as_deref() == Some("file bracket.ink")).expect("the walkthrough writes bracket.ink").body;
    let dir = Scratch::new("docs-picture");
    let shown = root("docs/guide/bracket.png");
    assert!(bracket.contains("lug 5 90)"), "the walkthrough's bracket has its 5 mm fillet");
    let path = dir.model("fillet-4.ink", &bracket.replace("lug 5 90)", "lug 4 90)"));
    let fillet = dir.path("fillet-4.png");
    dir.ink(["render", &path, "--out", &fillet, "--size", "640", "--res", "0.15"]).ok();
    let (n, all) = unlike_pixels(&shown, &fillet).unwrap();
    assert!(n as f64 > 4.0 * PICTURE_NOISE * all as f64, "a 4 mm fillet differs in {n} of {all} pixels: not well above the noise allowed");
}

#[test]
fn a_shown_answer_is_matched_to_its_digits() {
    let got: Value = serde_json::json!({"a": 7.375, "w": [22.6875, 0.0, 0.42], "rows": [{"n": "x", "v": 1}, {"n": "y", "v": 2}], "s": "holds"});
    assert!(differs(&shown(r#"{"a": 7.38, "w": [22.69, 0, 0.4], "rows": [{"n": "y"}], "s": "holds"}"#), &got, "").is_empty());
    assert!(!differs(&shown(r#"{"a": 7.36}"#), &got, "").is_empty(), "7.375 is not 7.36");
    assert!(!differs(&shown(r#"{"w": [22.69, 0]}"#), &got, "").is_empty(), "a list of numbers is shown whole");
    assert!(!differs(&shown(r#"{"rows": [{"n": "y"}, {"n": "x"}]}"#), &got, "").is_empty(), "in order");
    assert!(!differs(&shown(r#"{"s": "fails"}"#), &got, "").is_empty());
    assert!(!differs(&shown(r#"{"a": 7}"#), &got, "").is_empty(), "a whole number shown is exact");
    assert_eq!(words(r#"ink eval b.ink --points '[[0, 0, 2]]'   # a comment"#), ["ink", "eval", "b.ink", "--points", "[[0, 0, 2]]"]);
}
