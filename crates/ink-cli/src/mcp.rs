//! `ink mcp`: the CLI as an MCP server on stdio (newline-delimited JSON-RPC),
//! for agents whose clients speak MCP (`claude mcp add ink -- ink mcp`). One
//! tool runs any command: each call runs this binary once, so a call sees
//! exactly what the command line would print, and a crash leaves the server
//! standing. A package's command that reads words and asks to be answered
//! here (Run::Words: a catalog's search, from the catalog it has read) the
//! server answers itself, as the command line prints it, on a thread of its
//! own. Calls run side by side and the client can cancel one (its process
//! is killed; one answered here is told to stop). A call that runs long
//! reports its elapsed time to a client that asked for progress, and is
//! stopped at its time limit, with what it finished by then (the claims, pairs and designs the command reports as they are
//! decided, or where a solve had got, ink::claims::progress, given as the
//! command gives its answer: main's `publish`). The command is told its
//! limit (INK_DEADLINE, Unix seconds), so one that can stop early answers
//! whole before it. A picture a command writes comes back as an
//! image. An answer longer than a client takes (MAX_TEXT)
//! comes back shortened whole, its longest lists summarized, with the full
//! answer in a file. Every file a server writes is in one directory of its
//! own, removed when it ends (its input closed and its calls answered,
//! SIGTERM or SIGINT).

use serde_json::{Value as J, json};
use std::collections::HashMap;
use std::io::{BufRead, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use crate::NO_FILE;

/// The commands the tool runs: the command line's, but `mcp` itself. `help`
/// is the command line's text, which says what each option does.
fn commands() -> Vec<&'static str> {
    crate::commands(crate::registry()).into_iter().map(|(c, _)| c).filter(|&c| c != "mcp").collect()
}

/// A call's time limit unless it gives one (`timeout_s`), and how often a
/// running call reports progress.
const TIMEOUT_S: f64 = 300.0;
const PROGRESS_EVERY: Duration = Duration::from_secs(2);

/// The longest text a call answers with: clients cap a tool's result
/// (Claude Code at 25 000 tokens), and an answer over it is lost whole.
pub const MAX_TEXT: usize = 24_000;

/// Protocol versions served (the first when a client asks for another).
const VERSIONS: &[&str] = &["2025-06-18", "2025-03-26", "2024-11-05"];

/// What every agent reads when it connects: how the tool maps onto the
/// commands, then the reference's Tour (docs/IR.md), so a newcomer reading
/// the reference and an agent connecting read one text. A host may add its
/// own notes after it, so it says only what is true wherever ink runs.
fn instructions() -> String {
    const TOOL: &str = r#"The `ink` tool runs ink's commands on .ink files you write with your file tools. {"command": "assert", "file": "plate.ink", "options": {"param": ["hole=6.2"]}} runs `ink assert plate.ink --param hole=6.2`. Options go in options, never beside it. true sets a flag ({"options": {"own": true}}); a list of strings repeats an option."#;
    // Only a package's command reads words (Run::Words), so only a build
    // with one offers query.
    let words = if reads_words() { r#" A command that reads words, not a file, takes {"query": "WORDS"}."# } else { "" };
    format!("{TOOL}{words} Below is spec section Tour.\n\n{}", crate::cli::tour())
}

/// Whether a command of this build reads words (a package's: Run::Words),
/// so takes a query.
fn reads_words() -> bool {
    commands().iter().any(|c| crate::words(crate::registry(), c).is_some())
}

/// The tool's description, from the table `ink help` is printed from
/// (main's COMMANDS): a line per command, what it answers and its options
/// as JSON, the packages' after the core's (Command::tool). What each
/// option does is `help`'s, not repeated in every session. Not given: the
/// options the server sets (render's out, compact).
fn description() -> String {
    let hidden = |cmd: &str, k: &str| k == "compact" || (cmd, k) == ("render", "out");
    let mut d = String::from(
        r#"Run an ink command. Answers are JSON (spec, std: text). Options go in "options"; any command on a file takes {"param": ["w=40", "h=5"], "t": 1.5}; assert and interfere also {"nominal": true}. Units: mm, g, N, MPa, N mm, degC, s, deg; an option may give its own ("0.01in")."#,
    );
    for (name, help, opts) in crate::COMMANDS.iter().filter(|c| commands().contains(&c.0) && c.0 != "help") {
        let ex: Vec<String> = opts.iter().filter(|k| !hidden(name, k)).map(|k| format!("\"{k}\": {}", crate::option(name, k).0)).collect();
        d += &format!("\n{name}: {help}{}", if ex.is_empty() { String::new() } else { format!(" {{{}}}", ex.join(", ")) });
    }
    for c in crate::registry().commands() {
        d += &format!("\n{}: {}", c.name, c.tool);
    }
    d + "\nhelp: what each option does."
}

/// The calls under way by request id, each with its process once started:
/// a cancelled one is removed (its process killed), and then sends no reply.
/// A call is entered when it arrives, so a cancellation never misses one.
type Running = Arc<Mutex<HashMap<String, Option<Child>>>>;

/// A new file in the server's directory (`dir`, removed when it ends).
fn file(dir: &Path, ext: &str) -> PathBuf {
    static N: AtomicU64 = AtomicU64::new(0);
    dir.join(format!("{}.{ext}", N.fetch_add(1, Ordering::Relaxed)))
}

/// The server ends: its calls' processes are killed, its files removed, and
/// it exits holding the lock a call starts its process under, so none starts
/// after (`call`).
fn end(running: &Running, dir: &Path, code: i32) -> ! {
    let mut r = running.lock().unwrap();
    for mut c in r.drain().filter_map(|(_, c)| c) {
        let _ = c.kill();
        let _ = c.wait();
    }
    let _ = std::fs::remove_dir_all(dir);
    std::process::exit(code)
}

/// SIGTERM and SIGINT end the server at once, its calls killed. Blocked
/// here, before any other thread starts (main runs `mcp` on its own
/// thread), they are blocked in every thread
/// and taken by one that waits for them, which may do anything (a signal
/// handler may not). A process inherits the mask, so the calls' processes
/// start with it cleared (`unblocked`).
#[cfg(unix)]
fn end_on_signal(running: Running, dir: PathBuf) {
    // SAFETY: the set is initialized by sigemptyset before use; the calls
    // only read it and write `sig`.
    unsafe {
        let mut set: libc::sigset_t = std::mem::zeroed();
        libc::sigemptyset(&mut set);
        libc::sigaddset(&mut set, libc::SIGTERM);
        libc::sigaddset(&mut set, libc::SIGINT);
        libc::pthread_sigmask(libc::SIG_BLOCK, &set, std::ptr::null_mut());
        std::thread::spawn(move || {
            let mut sig = 0;
            if libc::sigwait(&set, &mut sig) == 0 {
                end(&running, &dir, 128 + sig);
            }
        });
    }
}

#[cfg(not(unix))]
fn end_on_signal(_: Running, _: PathBuf) {}

/// A call's process starts with no signal blocked, not the server's mask.
#[cfg(unix)]
fn unblocked(cmd: &mut Command) -> &mut Command {
    use std::os::unix::process::CommandExt;
    // SAFETY: runs between fork and exec, where sigemptyset and sigprocmask
    // may be called (they are async-signal-safe).
    unsafe {
        cmd.pre_exec(|| {
            let mut set: libc::sigset_t = std::mem::zeroed();
            libc::sigemptyset(&mut set);
            libc::sigprocmask(libc::SIG_SETMASK, &set, std::ptr::null_mut());
            Ok(())
        })
    }
}

#[cfg(not(unix))]
fn unblocked(cmd: &mut Command) -> &mut Command {
    cmd
}

pub fn serve() -> i32 {
    let out = Arc::new(Mutex::new(std::io::stdout()));
    let running: Running = Arc::default();
    // One directory per server (by process id: one left by a killed server
    // whose id this one reuses is replaced).
    let dir = std::env::temp_dir().join(format!("ink-mcp-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    if let Err(e) = std::fs::create_dir_all(&dir) {
        eprintln!("ink mcp: cannot make {}: {e}", dir.display());
        return 1;
    }
    end_on_signal(running.clone(), dir.clone());
    let send = {
        let out = out.clone();
        move |m: J| {
            let mut o = out.lock().unwrap();
            let _ = writeln!(o, "{m}").and_then(|_| o.flush());
        }
    };
    let mut calls = Vec::new();
    for line in std::io::stdin().lock().lines() {
        let Ok(line) = line else { break };
        if line.trim().is_empty() {
            continue;
        }
        let m = match serde_json::from_str::<J>(&line) {
            Ok(m) => m,
            Err(e) => {
                send(error(J::Null, -32700, format!("parse error: {e}")));
                continue;
            }
        };
        let (Some(method), true) = (m["method"].as_str(), m["jsonrpc"] == "2.0") else {
            send(error(m.get("id").cloned().unwrap_or(J::Null), -32600, "not a JSON-RPC 2.0 request".into()));
            continue;
        };
        let p = m.get("params").cloned().unwrap_or(json!({}));
        // A notification gets no reply; a cancellation kills its call.
        let Some(id) = m.get("id").cloned() else {
            if method == "notifications/cancelled"
                && let Some(Some(mut c)) = running.lock().unwrap().remove(&p["requestId"].to_string())
            {
                let _ = c.kill();
                let _ = c.wait();
            }
            continue;
        };
        match method {
            "tools/call" if p["name"] == "ink" => {
                running.lock().unwrap().insert(id.to_string(), None);
                let (send, running, dir) = (send.clone(), running.clone(), dir.clone());
                calls.retain(|c: &std::thread::JoinHandle<()>| !c.is_finished());
                calls.push(std::thread::spawn(move || {
                    // Progress goes to a client that asked for it (MCP's
                    // progressToken): the seconds the call has run.
                    let token = p["_meta"]["progressToken"].clone();
                    let progress = |secs: u64, what: &str| {
                        if !token.is_null() {
                            send(json!({"jsonrpc": "2.0", "method": "notifications/progress", "params": {
                                "progressToken": token, "progress": secs, "message": format!("{what}: running for {secs} s"),
                            }}));
                        }
                    };
                    match call(&p["arguments"], &id.to_string(), &running, &progress, &dir) {
                        Some(Ok(r)) => send(json!({"jsonrpc": "2.0", "id": id, "result": r})),
                        // The kernel crashed: the server's failure, not an answer.
                        Some(Err(e)) => send(error(id, -32603, e)),
                        None => {}
                    }
                }));
            }
            "tools/call" => send(error(id, -32602, format!("unknown tool {}", p["name"]))),
            "initialize" => {
                let asked = p["protocolVersion"].as_str().unwrap_or("");
                let version = VERSIONS.iter().find(|v| **v == asked).unwrap_or(&VERSIONS[0]);
                send(json!({"jsonrpc": "2.0", "id": id, "result": {
                    "protocolVersion": version,
                    "capabilities": {"tools": {}},
                    "serverInfo": {"name": "ink", "version": env!("CARGO_PKG_VERSION")},
                    "instructions": instructions(),
                }}));
            }
            "ping" => send(json!({"jsonrpc": "2.0", "id": id, "result": {}})),
            "tools/list" => send(json!({"jsonrpc": "2.0", "id": id, "result": {"tools": [tool()]}})),
            other => send(error(id, -32601, format!("unknown method {other}"))),
        }
    }
    // Input closed: every call it read is answered (a script piping its
    // requests closes its input after the last), then it ends. A host that
    // will not wait ends it with SIGTERM, which ends its calls too.
    for c in calls {
        let _ = c.join();
    }
    end(&running, &dir, 0)
}

fn error(id: J, code: i32, message: String) -> J {
    json!({"jsonrpc": "2.0", "id": id, "error": {"code": code, "message": message}})
}

fn tool() -> J {
    let mut t = json!({
        "name": "ink",
        "description": description(),
        "inputSchema": {
            "type": "object",
            "properties": {
                "command": {"type": "string", "enum": commands()},
                "file": {"type": "string", "description": if reads_words() { "The .ink file (none for help, spec, std and commands that read words)." } else { "The .ink file (none for help, spec and std)." }},
                "query": {"type": "string", "description": "The words for a command that reads words, not a file."},
                "options": {"type": "object", "description": "The command's options, as in the examples. true sets a flag; a list of strings repeats an option."},
                "timeout_s": {"type": "number", "description": format!("Seconds until the call stops (default {TIMEOUT_S}). It then answers with what finished.")},
            },
            "required": ["command"],
            "additionalProperties": false,
        },
    });
    if !reads_words() {
        t["inputSchema"]["properties"].as_object_mut().unwrap().remove("query");
    }
    t
}

/// Runs the command as the command line would, in a process of its own.
/// None: the call was cancelled.
fn call(args: &J, id: &str, running: &Running, progress: &dyn Fn(u64, &str), dir: &Path) -> Option<Result<J, String>> {
    let fail = |text: String| running.lock().unwrap().remove(id).map(|_| Ok(json!({"content": [{"type": "text", "text": text}], "isError": true})));
    let cmd = args["command"].as_str().unwrap_or("");
    if !commands().contains(&cmd) {
        return fail(format!("unknown command '{cmd}'; one of: {}", commands().join(", ")));
    }
    // An option beside the arguments, not in options, would be dropped and
    // the answer given without it: an error, as an unknown option is.
    const ARGS: [&str; 5] = ["command", "file", "query", "options", "timeout_s"];
    if let Some(k) = args.as_object().and_then(|a| a.keys().find(|k| !ARGS.contains(&k.as_str()))) {
        let args: Vec<&str> = ARGS.into_iter().filter(|&a| a != "query" || reads_words()).collect();
        return fail(format!("unknown argument '{k}': the arguments are {}; a command's options go in options ({{\"options\": {{\"{k}\": ...}}}})", args.join(", ")));
    }
    let mut argv = vec![cmd.to_string()];
    let reads_words = crate::words(crate::registry(), cmd).is_some();
    match args["file"].as_str() {
        Some(f) => argv.push(f.to_string()),
        None if !crate::no_file(crate::registry(), cmd) => return fail(format!("'{cmd}' needs a file")),
        None => {}
    }
    match args["query"].as_str() {
        Some(q) if reads_words => argv.extend(q.split_whitespace().map(str::to_string)),
        Some(_) => {
            let what = if crate::no_file(crate::registry(), cmd) { "what it is asked goes in options" } else { "it reads a file, and its options go in options" };
            return fail(format!("'{cmd}' takes no query: {what} (help: {{\"command\": \"help\"}})"));
        }
        None => {}
    }
    let limit = match &args["timeout_s"] {
        J::Null => TIMEOUT_S,
        v => match v.as_f64().filter(|t| *t > 0.0) {
            Some(t) => t,
            None => return fail(format!("timeout_s is a number of seconds above 0, got {v}")),
        },
    };
    let mut opts = match &args["options"] {
        J::Object(o) => o.clone(),
        J::Null => Default::default(),
        v => return fail(format!("options must be an object of option names, got {v}")),
    };
    // A picture comes back as an image: render needs somewhere to write it.
    let temp = (cmd == "render" && !opts.contains_key("out")).then(|| file(dir, "png"));
    if let Some(p) = &temp {
        opts.insert("out".into(), json!(p.to_string_lossy()));
    }
    // Columns rather than repeated keys: the same answer in fewer tokens.
    if !NO_FILE.contains(&cmd) {
        opts.entry("compact").or_insert(json!(true));
    }
    for (k, v) in &opts {
        // A list of strings repeats the option (--param, --attr); any other
        // value is the option's own (a JSON list of points is one --points).
        let repeated = v.as_array().filter(|l| !l.is_empty() && l.iter().all(J::is_string));
        for v in repeated.cloned().unwrap_or_else(|| vec![v.clone()]) {
            match v {
                J::Bool(true) => argv.push(format!("--{k}")),
                J::Bool(false) | J::Null => {}
                J::String(s) => argv.push(format!("--{k}={s}")),
                v => argv.push(format!("--{k}={v}")),
            }
        }
    }
    // A command answered here (a catalog's search, from the catalog this
    // server has read), on a thread of its own: stopped at its limit or
    // cancelled as a process is, it is told to stop and its answer is not
    // given.
    if crate::answered_here(&argv) {
        let stop = Arc::new(AtomicBool::new(false));
        let (tx, rx) = std::sync::mpsc::channel();
        let (flag, line) = (stop.clone(), argv.clone());
        std::thread::spawn(move || tx.send(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| crate::printed(&line, &|| flag.load(Ordering::Relaxed))))));
        let (start, mut reported) = (Instant::now(), Instant::now());
        let printed = loop {
            if let Ok(p) = rx.recv_timeout(Duration::from_millis(20)) {
                break p;
            }
            if !running.lock().unwrap().contains_key(id) {
                stop.store(true, Ordering::Relaxed);
                return None;
            }
            if start.elapsed().as_secs_f64() > limit {
                stop.store(true, Ordering::Relaxed);
                running.lock().unwrap().remove(id)?;
                let text = format!("ink {cmd}: stopped after {} s, its time limit. {}", start.elapsed().as_secs(), limited(None, start.elapsed().as_secs_f64()));
                return Some(Ok(json!({"content": [{"type": "text", "text": text}], "isError": true})));
            }
            if reported.elapsed() >= PROGRESS_EVERY {
                reported = Instant::now();
                progress(start.elapsed().as_secs(), &format!("ink {cmd}"));
            }
        };
        running.lock().unwrap().remove(id)?;
        return Some(match printed {
            Ok((text, ok)) => {
                let answer = serde_json::from_str(&text).unwrap_or(J::Null);
                Ok(json!({"content": [{"type": "text", "text": shorten(&text, answer, MAX_TEXT, dir)}], "isError": !ok}))
            }
            Err(e) => Err(format!("ink {cmd} crashed (panicked): {}", e.downcast_ref::<String>().map_or_else(|| e.downcast_ref::<&str>().copied().unwrap_or(""), String::as_str))),
        });
    }
    let exe = std::env::current_exe().unwrap_or_else(|_| "ink".into());
    // What the command finishes as it goes, for a call stopped at its limit.
    let log = file(dir, "progress");
    let cleanup = || {
        let _ = std::fs::remove_file(&log);
        if let Some(p) = &temp {
            let _ = std::fs::remove_file(p);
        }
    };
    // Started and entered under the lock a cancellation and the server's end
    // take: a call cancelled or a server ending first starts no process, and
    // one started is there to be killed.
    let mut r = running.lock().unwrap();
    let slot = r.get_mut(id)?;
    let deadline = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0.0, |d| d.as_secs_f64()) + limit;
    let mut child = match unblocked(&mut Command::new(exe))
        .args(&argv)
        .env("INK_PROGRESS", &log)
        .env("INK_DEADLINE", format!("{deadline:.3}"))
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
    {
        Ok(c) => c,
        Err(e) => {
            drop(r);
            return fail(format!("cannot run ink: {e}"));
        }
    };
    let (out, err) = (child.stdout.take().unwrap(), child.stderr.take().unwrap());
    *slot = Some(child);
    drop(r);
    let drain = |mut r: Box<dyn Read + Send>| {
        std::thread::spawn(move || {
            let mut b = Vec::new();
            let _ = r.read_to_end(&mut b);
            String::from_utf8_lossy(&b).into_owned()
        })
    };
    let (stdout, stderr) = (drain(Box::new(out)), drain(Box::new(err)));
    let (start, mut reported) = (Instant::now(), Instant::now());
    let status = loop {
        let mut r = running.lock().unwrap();
        let Some(Some(c)) = r.get_mut(id) else { break None };
        if let Ok(Some(s)) = c.try_wait() {
            r.remove(id);
            break Some(s);
        }
        if start.elapsed().as_secs_f64() > limit {
            let _ = c.kill();
            let _ = c.wait();
            r.remove(id);
            let so_far = Finished::read(&log);
            cleanup();
            let mut text = format!("ink {cmd}: stopped after {} s, its time limit", start.elapsed().as_secs());
            match &so_far {
                Some(f) if f.total.is_none() && f.done.is_empty() => text += &format!(", while {}", f.count()),
                Some(f) => text += &format!(", with {} (the answer so far is below)", f.count()),
                None => {}
            }
            text += &format!(". {}", limited(so_far.as_ref(), start.elapsed().as_secs_f64()));
            // The answer so far, given as the command gives its own
            // (compact or not).
            let a = crate::Args::parse(&argv);
            if let Some(f) = so_far {
                let answer = crate::publish(f.json(), &a, start.elapsed().as_secs_f64() * 1e3);
                let j = serde_json::from_str(&answer).unwrap_or(J::Null);
                text += &format!("\n{}", shorten(&answer, j, MAX_TEXT - text.len(), dir));
            }
            return Some(Ok(json!({"content": [{"type": "text", "text": text}], "isError": true})));
        }
        drop(r);
        if reported.elapsed() >= PROGRESS_EVERY {
            reported = Instant::now();
            let what = match Finished::read(&log) {
                Some(f) => format!("ink {cmd} ({})", f.count()),
                None => format!("ink {cmd}"),
            };
            progress(start.elapsed().as_secs(), &what);
        }
        std::thread::sleep(Duration::from_millis(20));
    };
    let (mut text, err) = (stdout.join().unwrap_or_default(), stderr.join().unwrap_or_default());
    let Some(status) = status else {
        cleanup();
        return None;
    };
    // An answer with a claim that does not hold is an answer, not an error.
    let ok = status.success() || status.code() == Some(crate::NOT_ALL_HOLD);
    // A panic (exit 101) or a fatal signal is the kernel's crash, not the call's answer; but
    // not SIGKILL (its memory's limit) nor running out of address space (Rust aborts saying so).
    use std::os::unix::process::ExitStatusExt;
    let crashed = match (status.code(), status.signal()) {
        (Some(101), _) => true,
        (None, Some(sig)) => sig != 9 && !(sig == 6 && err.contains("memory allocation of")),
        _ => false,
    };
    if crashed {
        cleanup();
        let how = status.signal().map_or("panicked".to_string(), |s| format!("killed by signal {s}"));
        return Some(Err(format!("ink {cmd} crashed ({how}): {}", shorten_err(&err))));
    }
    if !ok && text.trim().is_empty() {
        text = format!("ink {} failed: {err}", argv.join(" "));
    }
    // The pictures the command wrote: render's `out`, solve's `render.out`.
    let answer: J = serde_json::from_str(&text).unwrap_or(J::Null);
    let mut content = vec![json!({"type": "text", "text": shorten(&text, answer.clone(), MAX_TEXT, dir)})];
    for p in [&answer["out"], &answer["render"]["out"]].into_iter().filter_map(J::as_str) {
        if p.to_lowercase().ends_with(".png")
            && let Ok(png) = std::fs::read(p)
        {
            content.push(json!({"type": "image", "data": base64(&png), "mimeType": "image/png"}));
        }
    }
    cleanup();
    Some(Ok(json!({"content": content, "isError": !ok})))
}

/// A crash's stderr, its last 2000 characters (a panic's message is at its end).
fn shorten_err(err: &str) -> &str {
    let e = err.trim();
    let at = e.char_indices().rev().nth(1999).map_or(0, |(i, _)| i);
    &e[at..]
}

/// `text`, the answer (`answer`, if JSON), within `max` characters: else
/// written whole to a file, and of JSON its longest lists, largest first,
/// summarized in place (how many items, the first few) until it fits; other
/// text cut at a line. Each says where the whole answer is (there while
/// the server runs).
fn shorten(text: &str, mut answer: J, max: usize, dir: &Path) -> String {
    if text.len() <= max {
        return text.to_string();
    }
    let full = file(dir, "txt");
    let saved = match std::fs::write(&full, text) {
        Ok(()) => format!("the whole answer ({} characters) is in {}: read it with your file tools (jq for JSON)", text.len(), full.display()),
        Err(e) => format!("the whole answer ({} characters) could not be saved ({e}): run the command line for it", text.len()),
    };
    if answer.is_object() {
        // Room for the note, and the summaries' own keys.
        let room = max - max / 8;
        while answer.to_string().len() > room {
            let Some((list, _)) = longest_list(&mut answer) else { break };
            let items = std::mem::take(list);
            let k = items.len();
            let mut first = Vec::new();
            let mut used = 0;
            for x in items {
                used += x.to_string().len();
                if first.len() >= 3 || (!first.is_empty() && used > room / 16) {
                    break;
                }
                first.push(x);
            }
            *list = vec![json!({"items": k, "first": first, "rest": "in the whole answer (shortened)"})];
        }
        let note = format!("over {max} characters, so its longest lists are summarized as [{{items, first, rest}}]; {saved}");
        if answer.to_string().len() > room {
            // Long strings, not lists: only the note.
            let keys: Vec<&String> = answer.as_object().map(|o| o.keys().collect()).unwrap_or_default();
            return json!({"ok": answer["ok"], "keys": keys, "shortened": note}).to_string();
        }
        answer["shortened"] = json!(note);
        return answer.to_string();
    }
    let mut cut = 0;
    for (i, _) in text.match_indices('\n') {
        if i > max - max / 8 {
            break;
        }
        cut = i;
    }
    format!("{}\n... shortened: over {max} characters; {saved}", &text[..cut])
}

/// The longest list in `v` (as JSON text) of more than one item, and its
/// length: an outer list holds any inner one, so inner lists are reached
/// once their outer one is summarized, into its first items (a summary's
/// own list of them is not summarized again: each step leaves one list
/// fewer).
fn longest_list(v: &mut J) -> Option<(&mut Vec<J>, usize)> {
    match v {
        J::Array(l) if l.len() > 1 => {
            let k = v.to_string().len();
            let J::Array(l) = v else { unreachable!() };
            Some((l, k))
        }
        J::Array(l) => l.iter_mut().filter_map(longest_list).max_by_key(|x| x.1),
        J::Object(m) => {
            let summary = m.contains_key("rest") && m.contains_key("items");
            m.iter_mut()
                .filter_map(|(k, x)| match (summary, x) {
                    (false, x) => longest_list(x),
                    (true, J::Array(first)) if k == "first" => first.iter_mut().filter_map(longest_list).max_by_key(|x| x.1),
                    (true, _) => None,
                })
                .max_by_key(|x| x.1)
        }
        _ => None,
    }
}

/// What a command reported finishing (ink::claims::progress's lines): how
/// many of how many claims, pairs or designs (and whether each is kept in
/// the cache), each one done with its status, time and result (its row of
/// the finished answer), the one running (where they run one at a time),
/// the answer's other fields (its scope; an
/// optimization's best design so far; a solve's state), and what it was
/// last doing (`status`: where a solve had got).
struct Finished {
    status: Option<String>,
    total: Option<u64>,
    of: String,
    kept: bool,
    done: Vec<J>,
    /// The one running when stopped (of several, as interfere's moving
    /// pairs run, the one running longest), and when it started (Unix
    /// seconds).
    running: Option<(J, Option<f64>)>,
    /// The ranged params the one running reads, by name.
    ranged: Vec<String>,
    answer: serde_json::Map<String, J>,
}

/// What limited a call stopped at its time limit, `secs` after it started,
/// and what follows from that: the one running when stopped, where it had
/// run more than half the time (it needs more than the call gave, and may
/// have a call of its own: its claim's `name`, or a pair's assembly as the
/// `target`), else what finished; and what the same call with more time
/// does, by what the command said of what it finished (`kept`): claims and
/// pairs are cached as they are decided (unless INK_NO_CACHE), so it takes
/// them from there; an optimization's designs are trials, not cached, so
/// its search starts again.
fn limited(f: Option<&Finished>, secs: f64) -> String {
    let Some(f) = f else { return "It reported nothing as it ran: a larger timeout_s runs it again.".into() };
    if f.total.is_none() && f.done.is_empty() {
        return "Where it had got is below: a larger timeout_s runs it again from the start.".into();
    }
    let again = match (f.kept, std::env::var_os("INK_NO_CACHE")) {
        (true, None) => "What finished is cached, so the same call with a larger timeout_s takes it from there and goes on",
        (true, Some(_)) => "Nothing is cached here (INK_NO_CACHE), so a call with a larger timeout_s checks it all again",
        (false, _) => "What finished is not cached (each design is a trial), so a call with a larger timeout_s starts the search again",
    };
    let Some((r, at)) = f.running.as_ref() else {
        return format!("The time went to what finished and what was being checked when it stopped (the slowest under finished.slowest). {again}.");
    };
    let name = match r {
        J::Array(p) => p.iter().filter_map(J::as_str).collect::<Vec<_>>().join(" ~ "),
        J::String(n) => n.clone(),
        r => r.to_string(),
    };
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0.0, |d| d.as_secs_f64());
    match at.map(|at| (now - at).max(0.0)) {
        Some(ran) if ran > 0.5 * secs => {
            // A pair's own call: the instance holding both, where that is
            // narrower than the whole assembly.
            let holder = |p: &[J]| -> Option<String> {
                let [a, b] = p else { return None };
                let (a, b): (Vec<&str>, Vec<&str>) = (a.as_str()?.split('/').collect(), b.as_str()?.split('/').collect());
                let k = a.iter().zip(&b).take_while(|(x, y)| x == y).count();
                (k > 1 && k < a.len()).then(|| a[..k].join("/"))
            };
            let own = match r {
                J::String(n) => format!(", or give it a call of its own (name {n})"),
                J::Array(p) => holder(p).map_or(String::new(), |h| format!(", or check the assembly holding it alone (target {h})")),
                _ => String::new(),
            };
            format!("Most of the time went to {name}, still running after {ran:.0} s: it needs more than the call gave. {again}{own}.")
        }
        _ => format!("The time went to what finished (the slowest under finished.slowest), and {name} was running when it stopped. {again}."),
    }
}

impl Finished {
    /// None: the command reported nothing (not one that reports).
    fn read(path: &std::path::Path) -> Option<Finished> {
        let text = std::fs::read_to_string(path).ok()?;
        let mut f = Finished { status: None, total: None, of: "results".into(), kept: false, done: Vec::new(), running: None, ranged: Vec::new(), answer: Default::default() };
        // Those running, in the order they started, each with its ranged
        // params.
        let mut running: Vec<(J, Option<f64>, Vec<String>)> = Vec::new();
        // A line cut short by the stop does not parse: it is skipped.
        for l in text.lines().filter_map(|l| serde_json::from_str::<J>(l).ok()) {
            if let Some(s) = l["status"].as_str() {
                f.status = Some(s.to_string());
            }
            if let Some(of) = l["of"].as_str() {
                (f.total, f.of, f.kept) = (l["total"].as_u64(), of.to_string(), l["kept"] == true);
            } else if let J::Object(a) = &l["answer"] {
                f.answer.extend(a.clone());
            } else if !l["running"].is_null() {
                running.push((l["running"].clone(), l["at"].as_f64(), l["ranged"].as_object().map_or(Vec::new(), |o| o.keys().cloned().collect())));
            } else if !l["done"].is_null() {
                if let Some(k) = running.iter().position(|r| r.0 == l["done"]) {
                    running.remove(k);
                }
                f.done.push(l);
            }
        }
        if let Some((r, at, ranged)) = running.into_iter().next() {
            (f.running, f.ranged) = (Some((r, at)), ranged);
        }
        (f.total.is_some() || !f.done.is_empty() || f.status.is_some()).then_some(f)
    }

    /// "12 of 40 pairs done", or what it was doing; and what the one
    /// running ranges ("hand-clear over 4 ranged params: a, b, c, d").
    fn count(&self) -> String {
        let done = match (self.total, &self.status) {
            (Some(t), _) => format!("{} of {t} {} done", self.done.len(), self.of),
            (None, Some(s)) if self.done.is_empty() => s.clone(),
            (None, _) => format!("{} {} done", self.done.len(), self.of),
        };
        match (&self.running, self.ranged.len()) {
            (Some((J::String(n), _)), k) if k > 0 => format!("{done}; {n} over {k} ranged param{}: {}", if k == 1 { "" } else { "s" }, self.ranged.join(", ")),
            _ => done,
        }
    }

    /// The answer so far, as a finished run gives it: its fields besides
    /// its rows (its scope; an optimization's best
    /// design), each row done under the finished answer's key (`asserts`,
    /// `pairs`; `designs`), marked `stopped`, and what finished: statuses
    /// counted, the slowest so far, and the one running when stopped.
    fn json(&self) -> J {
        let ms = |l: &J| l["ms"].as_f64().unwrap_or(0.0);
        let mut counts = serde_json::Map::new();
        for l in &self.done {
            let s = l["status"].as_str().unwrap_or("unknown").to_string();
            let k = counts.get(&s).and_then(J::as_u64).unwrap_or(0);
            counts.insert(s, json!(k + 1));
        }
        let mut slowest: Vec<&J> = self.done.iter().collect();
        slowest.sort_by(|a, b| ms(b).total_cmp(&ms(a)));
        let mut o = J::Object(self.answer.clone());
        o["stopped"] = json!(true);
        if self.total.is_none() && self.done.is_empty() {
            o["status_when_stopped"] = json!(self.status);
            return o;
        }
        let rows: Vec<J> = self.done.iter().filter_map(|l| l.get("row")).cloned().collect();
        if !rows.is_empty() {
            o[&self.of] = J::Array(rows);
        }
        o["finished"] = json!({"done": self.done.len(), "total": self.total, "of": self.of, "counts": counts,
            "slowest": slowest.iter().take(5).map(|l| json!({"name": l["done"], "ms": (ms(l) * 10.0).round() / 10.0})).collect::<Vec<_>>()});
        if let Some((r, _)) = &self.running {
            o["finished"]["running_when_stopped"] = r.clone();
        }
        o
    }
}

fn base64(b: &[u8]) -> String {
    const A: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut s = String::with_capacity(b.len().div_ceil(3) * 4);
    for c in b.chunks(3) {
        let n = (c[0] as u32) << 16 | (*c.get(1).unwrap_or(&0) as u32) << 8 | *c.get(2).unwrap_or(&0) as u32;
        for i in 0..4 {
            s.push(if i <= c.len() { A[(n >> (18 - 6 * i) & 63) as usize] as char } else { '=' });
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A solve stopped at its time limit before it answered (a package's
    /// still stepping) gives where it had got: its last status and
    /// answer so far, not a count of rows it never had, nor "It reported
    /// nothing as it ran".
    #[test]
    fn a_solve_stopped_gives_where_it_had_got() {
        let log = std::env::temp_dir().join(format!("ink-mcp-status-{}", std::process::id()));
        std::fs::write(&log, "{\"status\":\"classifying\"}\n{\"status\":\"stepping: step 20\",\"answer\":{\"study\":\"duct\",\"solver\":{\"steps\":20}}}\n").unwrap();
        let f = Finished::read(&log).unwrap();
        let _ = std::fs::remove_file(&log);
        let j = f.json();
        assert_eq!(f.count(), "stepping: step 20");
        assert!(j["stopped"] == true && j["solver"]["steps"] == 20 && j["status_when_stopped"] == "stepping: step 20" && j.get("finished").is_none(), "{j}");
        assert!(limited(Some(&f), 10.0).starts_with("Where it had got is below"));
    }

    /// The progress of a run names what is running: the claim's range (the
    /// ranged params its line gave as it started), and of pairs side by
    /// side, the longest.
    #[test]
    fn a_run_in_progress_names_what_is_running() {
        let log = std::env::temp_dir().join(format!("ink-mcp-sweep-{}", std::process::id()));
        std::fs::write(&log, "{\"total\":1,\"of\":\"asserts\",\"kept\":true}\n{\"running\":\"hand-clear\",\"ranged\":{\"pitch0\":[-100,100],\"roll0\":[-180,180]}}\n").unwrap();
        let f = Finished::read(&log).unwrap();
        let _ = std::fs::remove_file(&log);
        assert_eq!(f.count(), "0 of 1 asserts done; hand-clear over 2 ranged params: pitch0, roll0");
        // Pairs run side by side (interfere's moving pairs): a stopped run
        // names the one running longest, not the last begun, whichever
        // finished. It named none once a pair begun after a long one finished.
        std::fs::write(&log, "{\"total\":3,\"of\":\"pairs\",\"kept\":true}\n{\"running\":[\"a\",\"b\"],\"at\":10}\n{\"running\":[\"a\",\"c\"],\"at\":11}\n{\"running\":[\"a\",\"d\"],\"at\":12}\n{\"done\":[\"a\",\"d\"],\"status\":\"holds\",\"ms\":5}\n").unwrap();
        let f = Finished::read(&log).unwrap();
        let _ = std::fs::remove_file(&log);
        assert_eq!(f.json()["finished"]["running_when_stopped"], json!(["a", "b"]));
    }
}
