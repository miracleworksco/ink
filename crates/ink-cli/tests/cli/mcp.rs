//! `ink mcp` driven as an agent's client drives it: JSON-RPC lines on its
//! standard input and output, its calls answered as the command line
//! answers them, the reference's sections, and its calls' processes found
//! and signalled through /proc and `kill`: cancelled, stopped at their
//! limits, crashed and ended with it. (A parts search answered in the
//! server is the catalog package's: ext's tests.) The tests that find a
//! call's process read /proc, so they run on Linux only.

use crate::common::*;
use serde_json::{Value, json};
use std::io::{BufRead, Write};
use std::path::Path;
use std::process::{Child, ChildStdin, Command, ExitStatus, Stdio};
use std::sync::mpsc::Receiver;
use std::time::{Duration, Instant};

/// f(0..n) in order, the calls of a test that do not depend on each other
/// run side by side, as many at once as the machine has cores: on a small
/// one (a 2-core CI runner) each call keeps its own time, so a reply's
/// minute still means a hang.
pub fn par<T: Send>(n: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let next = std::sync::atomic::AtomicUsize::new(0);
    let workers = std::thread::available_parallelism().map_or(4, |c| c.get()).min(n);
    let mut done: Vec<(usize, T)> = std::thread::scope(|s| {
        let hs: Vec<_> = (0..workers)
            .map(|_| {
                s.spawn(|| {
                    let mut out = Vec::new();
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if i >= n {
                            break out;
                        }
                        out.push((i, f(i)));
                    }
                })
            })
            .collect();
        hs.into_iter().flat_map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e))).collect()
    });
    done.sort_by_key(|d| d.0);
    done.into_iter().map(|d| d.1).collect()
}

/// A number in an answer (NaN where it is none).
pub fn num(v: &Value) -> f64 {
    v.as_f64().unwrap_or(f64::NAN)
}

/// cond(), polled until it is something, for up to a minute (past that the
/// server hangs): what a server does at once is waited for by its state,
/// never for a fixed time.
#[cfg(target_os = "linux")]
pub fn until<T>(mut cond: impl FnMut() -> Option<T>) -> Option<T> {
    let t = Instant::now();
    loop {
        if let Some(v) = cond() {
            return Some(v);
        }
        if t.elapsed() > Duration::from_secs(60) {
            return None;
        }
        std::thread::sleep(Duration::from_millis(20));
    }
}

/// `kill -SIG pid`.
#[cfg(target_os = "linux")]
pub fn signal(sig: &str, pid: u32) {
    let ok = Command::new("kill").arg(format!("-{sig}")).arg(pid.to_string()).status().is_ok_and(|s| s.success());
    assert!(ok, "kill -{sig} {pid} failed");
}

/// A process's command line, its arguments.
#[cfg(target_os = "linux")]
pub fn cmdline(pid: u32) -> Vec<String> {
    std::fs::read(format!("/proc/{pid}/cmdline")).map_or(vec![], |b| b.split(|&c| c == 0).map(|a| String::from_utf8_lossy(a).into_owned()).collect())
}

/// Every process whose command line has an argument `arg` (a path of a
/// test's own scratch directory, so no other's).
#[cfg(target_os = "linux")]
pub fn running_with(arg: &str) -> Vec<u32> {
    std::fs::read_dir("/proc").unwrap().filter_map(|e| e.ok()?.file_name().to_str()?.parse::<u32>().ok()).filter(|&p| cmdline(p).iter().any(|a| a == arg)).collect()
}

/// `ink mcp` run as a client runs it: JSON-RPC lines in and out. Dropped,
/// it is ended by SIGTERM (as a host ends it), so its calls' processes go
/// with it.
pub struct Server {
    child: Child,
    stdin: Option<ChildStdin>,
    lines: Receiver<String>,
}

impl Server {
    /// A server in `cwd`, its files under `tmp` (TMPDIR), with `env` more.
    pub fn start(cwd: &str, tmp: &Path, env: &[(&str, &str)]) -> Server {
        Server::spawn(cwd, tmp, env, true)
    }
    /// One whose output nobody reads.
    #[cfg(target_os = "linux")]
    pub fn unread(cwd: &str, tmp: &Path) -> Server {
        Server::spawn(cwd, tmp, &[], false)
    }
    fn spawn(cwd: &str, tmp: &Path, env: &[(&str, &str)], read: bool) -> Server {
        let mut cmd = Command::new(INK);
        cmd.arg("mcp").current_dir(cwd).env("RAYON_NUM_THREADS", THREADS).env("TMPDIR", tmp);
        for (k, v) in env {
            cmd.env(k, v);
        }
        cmd.stdin(Stdio::piped()).stdout(if read { Stdio::piped() } else { Stdio::null() }).stderr(Stdio::null());
        let mut child = cmd.spawn().unwrap_or_else(|e| panic!("cannot run {INK} mcp: {e}"));
        let stdin = child.stdin.take();
        let (tx, lines) = std::sync::mpsc::channel();
        if let Some(out) = child.stdout.take() {
            std::thread::spawn(move || {
                for l in std::io::BufReader::new(out).lines() {
                    let Ok(l) = l else { break };
                    if tx.send(l).is_err() {
                        break;
                    }
                }
            });
        }
        Server { child, stdin, lines }
    }
    pub fn pid(&self) -> u32 {
        self.child.id()
    }
    /// Its directory, where everything it writes is.
    pub fn dir(&self, tmp: &Path) -> std::path::PathBuf {
        tmp.join(format!("ink-mcp-{}", self.pid()))
    }
    /// Lines to its input, as they are.
    pub fn write(&mut self, text: &str) {
        let i = self.stdin.as_mut().expect("its input is closed");
        i.write_all(text.as_bytes()).unwrap();
        i.flush().unwrap();
    }
    /// A message (jsonrpc 2.0) to it.
    pub fn send(&mut self, mut m: Value) {
        m["jsonrpc"] = json!("2.0");
        self.write(&format!("{m}\n"));
    }
    /// The next line it writes, as JSON. Five minutes at most: a guard
    /// against a server that never answers, not a limit on how long an
    /// answer may take (a 10 s claim took over a minute on a loaded 4-core
    /// CI runner).
    pub fn next(&self) -> Value {
        let l = self.lines.recv_timeout(Duration::from_secs(300)).unwrap_or_else(|e| panic!("ink mcp wrote nothing in 5 min ({e})"));
        serde_json::from_str(&l).unwrap_or_else(|e| panic!("ink mcp wrote non-JSON ({e}): {l}"))
    }
    /// A request and its reply (the next line).
    pub fn rpc(&mut self, id: u64, method: &str, params: Value) -> Value {
        self.send(json!({"id": id, "method": method, "params": params}));
        self.next()
    }
    /// A tool call: whether it is an error, and its content.
    pub fn call(&mut self, id: u64, args: Value) -> (bool, Vec<Value>) {
        self.call_on(id, "ink", args)
    }
    /// A call to the tool named.
    pub fn call_on(&mut self, id: u64, tool: &str, args: Value) -> (bool, Vec<Value>) {
        let r = self.rpc(id, "tools/call", json!({"name": tool, "arguments": args}));
        let res = &r["result"];
        assert!(res.is_object(), "call {id}: no result: {}", pretty(&r));
        (res["isError"] == true, res["content"].as_array().unwrap().clone())
    }
    /// A call's text, as JSON.
    pub fn answer(&mut self, id: u64, args: Value) -> (bool, Value) {
        let (err, c) = self.call(id, args);
        let t = c[0]["text"].as_str().unwrap();
        (err, serde_json::from_str(t).unwrap_or_else(|e| panic!("call {id}: non-JSON ({e}): {t}")))
    }
    /// Its input closed; how it exits (10 s at most).
    pub fn close(&mut self) -> ExitStatus {
        self.stdin = None;
        self.wait()
    }
    pub fn wait(&mut self) -> ExitStatus {
        let t = Instant::now();
        loop {
            if let Some(s) = self.child.try_wait().unwrap() {
                return s;
            }
            assert!(t.elapsed() < Duration::from_secs(10), "ink mcp still running 10 s after it was ended");
            std::thread::sleep(Duration::from_millis(10));
        }
    }
    /// Its calls' processes (its threads' children) with an argument `arg`.
    #[cfg(target_os = "linux")]
    pub fn calls_with(&self, arg: &str) -> Vec<u32> {
        let tasks = format!("/proc/{}/task", self.pid());
        let Ok(ts) = std::fs::read_dir(&tasks) else { return vec![] };
        ts.filter_map(|t| std::fs::read_to_string(t.ok()?.path().join("children")).ok())
            .flat_map(|c| c.split_whitespace().filter_map(|p| p.parse().ok()).collect::<Vec<u32>>())
            .filter(|&p| cmdline(p).iter().any(|a| a == arg))
            .collect()
    }
    /// Whether any of its threads has a child process.
    #[cfg(target_os = "linux")]
    pub fn has_children(&self) -> bool {
        let tasks = format!("/proc/{}/task", self.pid());
        std::fs::read_dir(&tasks).is_ok_and(|ts| ts.filter_map(|t| std::fs::read_to_string(t.ok()?.path().join("children")).ok()).any(|c| !c.trim().is_empty()))
    }
}

impl Drop for Server {
    fn drop(&mut self) {
        if self.child.try_wait().ok().flatten().is_none() {
            _ = Command::new("kill").arg("-TERM").arg(self.pid().to_string()).status();
            let t = Instant::now();
            while self.child.try_wait().ok().flatten().is_none() && t.elapsed() < Duration::from_secs(10) {
                std::thread::sleep(Duration::from_millis(10));
            }
            _ = self.child.kill();
            _ = self.child.wait();
        }
    }
}

/// A FIFO nobody writes: a model whose load never ends.
#[cfg(target_os = "linux")]
pub fn fifo(s: &Scratch) -> String {
    let p = s.path("never.ink");
    let ok = Command::new("mkfifo").arg(&p).status().is_ok_and(|s| s.success());
    assert!(ok, "mkfifo {p} failed");
    p
}

fn text(c: &[Value]) -> &str {
    c[0]["text"].as_str().unwrap_or("")
}

fn cut(s: &str) -> &str {
    &s[..s.floor_char_boundary(400)]
}

#[test]
fn mcp_initializes_lists_its_tool_and_ends_with_its_input() {
    let s = Scratch::new("mcp-protocol");
    let mut srv = Server::start(ROOT, &s, &[]);
    let init = srv.rpc(1, "initialize", json!({"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "cli-test", "version": "0"}}));
    let tools = srv.rpc(2, "tools/list", json!({}));
    let (init, tools) = (&init["result"], &tools["result"]["tools"]);
    assert!(
        init["capabilities"] == json!({"tools": {}})
            && init["instructions"].as_str().is_some_and(|i| i.contains("spec"))
            && tools.as_array().map(|t| t.len()) == Some(2)
            && tools[0]["name"] == "ink"
            && tools[0]["inputSchema"]["properties"].get("query").is_none()
            && tools[0]["inputSchema"]["properties"]["file"]["description"].as_str().is_some_and(|d| !d.contains("words"))
            && !init["instructions"].as_str().unwrap_or("").contains("query")
            && tools[0]["inputSchema"]["properties"]["command"]["enum"].as_array().is_some_and(|e| e.contains(&json!("assert"))),
        "initializes, with instructions, and lists the ink tool: {}",
        pretty(init)
    );
    // The description gives each command's options, and help what each
    // does: said once, not in every session's description.
    let (err, c) = srv.call(20, json!({"command": "help"}));
    let desc = tools[0]["description"].as_str().unwrap_or("");
    assert!(
        !err && text(&c).contains("--keys")
            && desc.contains("\"keys\": ")
            && desc.ends_with("\nhelp: what each option does.")
            && !desc.contains("keep only these keys")
            && desc.len() < 240 * desc.matches('\n').count(),
        "the description ({} characters, {} lines: within 240 a line) gives the options, help what each does: {}",
        desc.len(),
        desc.matches('\n').count(),
        &desc[desc.floor_char_boundary(desc.len().saturating_sub(300))..]
    );
    let init = srv.rpc(7, "initialize", json!({"protocolVersion": "2099-01-01"}));
    assert!(init["result"]["protocolVersion"] == "2025-06-18", "a protocol version it does not serve gets one it does: {}", pretty(&init));
    // An answer over what a client takes comes back shortened, whole JSON:
    // its longest lists summarized, the full answer in a file.
    let src = format!("(import \"std\")\n{}(body b (shape (sphere 1)))\n", (0..3000).map(|k| format!("(fn helper-number-{k} (x) (+ x {k}))\n")).collect::<String>());
    let big = s.model("many-fns.ink", &src);
    let (err, c) = srv.call(16, json!({"command": "check", "file": big}));
    let t = text(&c);
    let d: Value = serde_json::from_str(t).unwrap_or_else(|e| panic!("a shortened answer is whole JSON ({e}): {}", cut(t)));
    let fns = if d["fns"].as_array().map(Vec::len) == Some(1) { d["fns"][0].clone() } else { json!({}) };
    let full = d["shortened"].as_str().unwrap_or("").rsplit(" is in ").next().unwrap().split(':').next().unwrap().to_string();
    let saved: Value = std::fs::read_to_string(&full).ok().and_then(|f| serde_json::from_str(&f).ok()).unwrap_or_default();
    assert!(
        !err && t.len() <= 24000
            && fns["items"] == 3000
            && fns["first"].as_array().map(Vec::len) == Some(3)
            && d["bodies"][0]["name"] == "b"
            && saved["fns"].as_array().map(Vec::len) == Some(3000),
        "a long answer is shortened to whole JSON, the full one saved ({full}): {}",
        cut(t)
    );
    let st = srv.close();
    assert!(st.code() == Some(0), "exits when its input ends: {st:?}");
    // Everything a server writes is in its own directory, removed when it
    // ends (/tmp is often memory).
    let mine = srv.dir(&s);
    assert!(Path::new(&full).parent() == Some(&*mine) && !mine.exists(), "its files ({full}) are in {mine:?}, removed when its input ends");
}

/// Every file under `dir` but a server's own directory, with its bytes.
pub fn files_in(dir: &Path) -> std::collections::BTreeMap<std::path::PathBuf, Vec<u8>> {
    let mut found = std::collections::BTreeMap::new();
    for e in std::fs::read_dir(dir).unwrap().flatten() {
        let p = e.path();
        match p.is_dir() {
            true if !e.file_name().to_string_lossy().starts_with("ink-mcp-") => found.extend(files_in(&p)),
            true => {}
            false => _ = found.insert(p.clone(), std::fs::read(&p).unwrap_or_default()),
        }
    }
    found
}

/// The destinations a command may be given, one of each kind of file ink
/// writes; each is seeded with other bytes first, so an overwrite shows.
pub const DESTS: [&str; 6] = ["probe.stl", "probe.3mf", "probe.glb", "probe.dxf", "probe.png", "probe.json"];

/// What a call did: the tool that took it (the other refused it), whether
/// it failed, its text, and the files it made or changed.
pub struct Took {
    pub tool: &'static str,
    pub err: bool,
    pub text: String,
    pub changed: Vec<std::path::PathBuf>,
}

/// `args` to both tools, in `dir`, its destinations seeded first. Panics
/// unless exactly one takes it.
pub fn to_one_tool(dir: &Path, call_on: &mut dyn FnMut(&str, Value) -> (bool, Vec<Value>), args: Value) -> Took {
    for d in DESTS {
        std::fs::write(dir.join(d), "old").unwrap();
    }
    let mut took = Vec::new();
    for tool in ["ink", "ink_write"] {
        let before = files_in(dir);
        let (err, c) = call_on(tool, args.clone());
        let text = c[0]["text"].as_str().unwrap_or("").to_string();
        if err && text.contains("tool with these arguments") {
            continue;
        }
        let changed = files_in(dir).into_iter().filter(|(p, b)| before.get(p) != Some(b)).map(|(p, _)| p).collect();
        took.push(Took { tool, err, text, changed });
    }
    assert!(took.len() == 1, "{args} is one tool's, not {:?}", took.iter().map(|t| t.tool).collect::<Vec<_>>());
    took.pop().unwrap()
}

/// Each tool gives a title and all four hints, which a client decides by
/// whether to ask before a call. Every option the command table lists is
/// tried, at each kind of destination: every call is one tool's (the other
/// refuses it before it runs); `ink` (read-only) changes no file; and each
/// option that writes, called on `ink_write`, makes or overwrites the file
/// it names.
#[test]
fn mcp_every_call_is_one_tools_and_its_hints_say_what_it_does() {
    use ink::lang::command::Writes;
    let s = Scratch::new("mcp-hints");
    let mut srv = Server::start(&s.to_string_lossy(), &s, &[]);
    let tools = srv.rpc(1, "tools/list", json!({}))["result"]["tools"].as_array().cloned().unwrap_or_default();
    let hints = |t: &Value| ["readOnlyHint", "destructiveHint", "idempotentHint", "openWorldHint"].map(|h| t["annotations"][h].as_bool());
    for t in &tools {
        assert!(t["title"].is_string() && hints(t).iter().all(Option::is_some), "{} gives a title and all four hints: {}", t["name"], pretty(t));
    }
    let names: Vec<&str> = tools.iter().filter_map(|t| t["name"].as_str()).collect();
    assert!(names == ["ink", "ink_write"], "a tool that reads and one that writes: {names:?}");
    let hint = |name: &str| hints(tools.iter().find(|t| t["name"] == name).unwrap());
    assert!(hint("ink") == [Some(true), Some(false), Some(true), Some(false)] && hint("ink_write") == [Some(false), Some(true), Some(true), Some(false)], "{tools:?}");
    let meta = &tools[1]["_meta"]["ink/writes"];
    assert!(
        meta["mesh"] == json!({"always": true, "given": []})
            && meta["slice"] == json!({"always": false, "given": ["out"]})
            && meta["render"] == json!({"always": false, "given": ["out"]})
            && meta.get("check").is_none(),
        "ink_write says which calls write: {meta}"
    );
    let commands: Vec<_> = ink_cli::commands(ink::registry()).into_iter().filter(|c| c.0 != "mcp").collect();
    // Each command on a server of its own, side by side.
    let wrote = par(commands.len(), |i| {
        let (cmd, opts) = commands[i];
        let s = Scratch::new(&format!("mcp-hints-{cmd}"));
        let model = s.model("ball.ink", "(import \"std\")\n(body b (shape (sphere 5)))\n");
        let file = (!ink_cli::NO_FILE.contains(&cmd)).then_some(&model);
        let mut srv = Server::start(&s.to_string_lossy(), &s, &[]);
        let mut id = 0;
        let mut call_on = |tool: &str, args: Value| {
            id += 1;
            srv.call_on(id, tool, args)
        };
        let mut wrote = Vec::new();
        let base = to_one_tool(&s, &mut call_on, json!({"command": cmd, "file": file}));
        assert!(base.tool == "ink" || !base.changed.is_empty(), "ink {cmd} on {} changed nothing: {}", base.tool, cut(&base.text));
        for o in opts {
            let mut written = false;
            for d in DESTS {
                let t = to_one_tool(&s, &mut call_on, json!({"command": cmd, "file": file, "options": {o.name: d}}));
                assert!(t.tool == "ink_write" || t.changed.is_empty(), "ink {cmd} {{{}: {d}}} ran on the read-only tool and changed {:?}: {}", o.name, t.changed, cut(&t.text));
                written |= o.writes != Writes::No && !t.err && t.tool == "ink_write" && t.changed.contains(&s.join(d));
            }
            assert!(o.writes == Writes::No || written, "ink {cmd} --{} writes, yet wrote no destination it was given", o.name);
            if written {
                wrote.push(format!("{cmd} {}", o.name));
            }
        }
        wrote
    });
    let wrote: Vec<String> = wrote.into_iter().flatten().collect();
    assert!(wrote == ["slice out", "mesh out", "render out"], "the options that write wrote: {wrote:?}");
}

/// A call names a write only by an option the command reads: an option
/// named otherwise (`out=probe.dxf`, which the command line would read as
/// out) is refused on either tool, and a file or query is never read as an
/// option. Nothing is written.
#[test]
fn mcp_no_encoding_of_an_option_or_word_writes_a_file() {
    let s = Scratch::new("mcp-encodings");
    let model = s.model("ball.ink", "(import \"std\")\n(body b (shape (sphere 5)))\n");
    let mut srv = Server::start(&s.to_string_lossy(), &s, &[]);
    let mut id = 0;
    let mut call_on = |tool: &str, args: Value| {
        id += 1;
        srv.call_on(id, tool, args)
    };
    for opts in [json!({"out=probe.dxf": true}), json!({"out=probe.dxf": "x"}), json!({"OUT": "probe.dxf"}), json!({"out ": "probe.dxf"}), json!({"--out": "probe.dxf"})] {
        for tool in ["ink", "ink_write"] {
            let before = files_in(&s);
            let (err, c) = call_on(tool, json!({"command": "slice", "file": model, "options": opts}));
            assert!(err && text(&c).contains("is no option name") && files_in(&s) == before, "{tool} slice {opts}: {}", cut(text(&c)));
        }
    }
    // A value is the option's own, whatever it holds; a file is a file.
    for args in [
        json!({"command": "slice", "file": model, "options": {"keys": "kind --out=probe.dxf"}}),
        json!({"command": "check", "file": "--out=probe.dxf"}),
        json!({"command": "help", "query": "--out=probe.dxf"}),
    ] {
        let t = to_one_tool(&s, &mut call_on, args.clone());
        assert!(t.tool == "ink" && t.changed.is_empty(), "{args} on {}: {:?} {}", t.tool, t.changed, cut(&t.text));
    }
}

#[test]
fn every_option_example_in_the_description_runs() {
    // An agent copies them: each is checked to run, so slice's and ray's
    // keys are their own, not eval's, and slice's out is not a mesh's (.stl).
    let s = Scratch::new("mcp-examples");
    let model = s.model(
        "examples.ink",
        "(import \"std\")\n(param w 40 (tol -1 1))\n(param h 5)\n(def reach (* 2 w))\n(t-span 0 2)\n(body arm (shape (box w 10 h)) (attr density 1.24g/cm3))\n(body base (shape (move (box 60 60 5) 0 0 -10)))\n(assembly root (instance arm arm) (instance base base))\n(assert fit (clearance (shape-of arm) (shape-of base) 1))\n(optimize lightest (vary w 30 50) (minimize (volume (shape-of arm))) (budget 4))\n",
    );
    // Run in the scratch directory, where the examples' files are written.
    let mut srv = Server::start(&s.to_string_lossy(), &s, &[]);
    let tools = srv.rpc(1, "tools/list", json!({}));
    let mut ran = 0;
    for t in tools["result"]["tools"].as_array().unwrap() {
        let (tool, desc) = (t["name"].as_str().unwrap(), t["description"].as_str().unwrap_or(""));
        for line in desc.lines().skip(1) {
            let Some((cmd, rest)) = line.split_once(": ") else { continue };
            let Some(mut call) = rest.rfind(" {").and_then(|i| serde_json::from_str::<Value>(&rest[i + 1..]).ok()) else { continue };
            // Each option alone, with the others the command needs.
            let alone = call.as_object().unwrap().clone();
            for (k, v) in alone.iter().filter(|_| cmd == "spec" || cmd == "std") {
                ran += 1;
                let (err, c) = srv.call_on(ran, tool, json!({"command": cmd, "options": {k: v}}));
                assert!(!err, "{cmd} {{\"{k}\": {v}}}, as the description gives it: {}", cut(text(&c)));
            }
            if cmd == "ray" {
                call.as_object_mut().unwrap().remove("rays");
            }
            if cmd != "spec" && cmd != "std" {
                ran += 1;
                let (err, c) = srv.call_on(ran, tool, json!({"command": cmd, "file": model, "options": call}));
                assert!(!err, "{tool} {cmd} {call}, as the description gives it: {}", cut(text(&c)));
            }
        }
    }
    assert!(ran > 18, "the descriptions' examples: {ran} of them in {tools}");
}

#[test]
fn mcp_a_call_is_the_command_lines_answer() {
    let s = Scratch::new("mcp-calls");
    let shelf = fixture("shelf.ink");
    let mut srv = Server::start(ROOT, &s, &[]);
    let (err, c) = srv.answer(3, json!({"command": "eval", "file": shelf, "options": {"points": [[0, 0, 0], [1, 2, 3]], "keys": "point"}}));
    assert!(!err && c["targets"][0]["results"]["rows"].as_array().map(Vec::len) == Some(2), "a list of points is one --points: {}", pretty(&c));
    let (err, c) = srv.answer(19, json!({"command": "ray", "file": shelf, "options": {"rays": [[[0, 75, 100], [0, 0, -1]]], "keys": "nearest.point"}}));
    let hit = &c["rays"]["rows"][0][0];
    assert!(!err && [0.0, 75.0, 63.0].iter().enumerate().all(|(k, &x)| (hit[k].as_f64().unwrap_or(f64::NAN) - x).abs() < 1e-9), "ray's keys pick the hit point: {}", pretty(&c));
    let (err, c) = srv.call(4, json!({"command": "render", "file": shelf, "options": {"size": 64}}));
    assert!(
        !err && c.iter().map(|x| x["type"].as_str().unwrap_or("")).collect::<Vec<_>>() == ["text", "image"] && c[1]["mimeType"] == "image/png",
        "render comes back as an image: {:?}",
        c.iter().map(|x| &x["type"]).collect::<Vec<_>>()
    );
    // An empty list repeats an option no times: no params set.
    let (err, c) = srv.answer(21, json!({"command": "check", "file": shelf, "options": {"param": []}}));
    assert!(!err && c["ok"] == true, "an empty list sets nothing: {}", pretty(&c));
    // help asked of one command gives it and what its options do.
    let (err, c) = srv.call(22, json!({"command": "help", "query": "assert"}));
    assert!(!err && text(&c).starts_with("usage: ink assert <file.ink> [options]") && text(&c).contains("  --own "), "{}", cut(text(&c)));
    let (err, c) = srv.call(5, json!({"command": "check", "file": shelf, "options": {"bogus": 1}}));
    assert!(err && text(&c).contains("unknown option --bogus"), "an error is the command line's, flagged as one: {}", cut(text(&c)));
    // An option beside options is rejected with an error showing where it
    // belongs; dropped, it would answer a mass at the default param.
    let (err, c) = srv.call(6, json!({"command": "check", "file": shelf, "param": ["w=8"]}));
    // check warns of what interfere asks at one moment, over MCP as on the
    // command line: the default target of bodies, not a zone no pair asks.
    let moving = s.model("mcp-bodies.ink", "(import \"std\")\n(body a (shape (move (sphere 1) (* 10 t) 0 0)))\n(body b (shape (move (sphere 1) 1 0 0)))\n");
    let apart = s.model(
        "mcp-unused-zone.ink",
        "(import \"std\")\n(body ball (shape (sphere 1)))\n(assembly root (instance a ball) (instance b ball (at (move 10 0 0)))\n(contact a b (move (sphere 10) (* 100 t) 0 0)))\n",
    );
    let warns = |srv: &mut Server, id, f: &str| srv.answer(id, json!({"command": "check", "file": f})).1["warnings"].to_string().contains("states no t-span");
    assert!(warns(&mut srv, 8, &moving) && !warns(&mut srv, 9, &apart), "check's one-moment warning is interfere's question");
    // No command of the kernel's reads words: no query is offered, and one
    // given to spec says where its question goes.
    let (qerr, q) = srv.call(7, json!({"command": "spec", "query": "Motion"}));
    assert!(qerr && text(&q).contains("'spec' takes no query: what it is asked goes in options"), "{}", cut(text(&q)));
    assert!(
        err && text(&c).contains("unknown argument 'param'") && text(&c).contains(r#"{"options": {"param""#),
        "an argument beside options is an error saying where it goes: {}",
        cut(text(&c))
    );
}

/// A target path over MCP is asked as on the command line, and every pair
/// reading t over the model's spans (motion::REPROS), with the contacts
/// parents declare between a target's parts.
#[test]
fn mcp_a_target_path_is_asked_as_its_definition() {
    let s = Scratch::new("mcp-target-path");
    let mut srv = Server::start(ROOT, &s, &[]);
    for (id, r) in crate::motion::REPROS.iter().enumerate() {
        let sub = Scratch::new(&format!("mcp-target-path-{id}"));
        let file = r.write(&sub);
        let options = r.target.map_or(json!({}), |t| json!({"target": t}));
        let (err, c) = srv.answer(id as u64 + 1, json!({"command": "interfere", "file": file, "options": options}));
        assert!(!err, "{}: {}", r.name, pretty(&c));
        r.check(&c).unwrap_or_else(|e| panic!("over MCP, {e}"));
    }
    let contact = s.model("ancestor-contact.ink", crate::assemblies::ANCESTOR_CONTACT);
    for (id, target) in ["root/carried", "root/m/q"].into_iter().enumerate() {
        let (err, c) = srv.answer(100 + id as u64, json!({"command": "interfere", "file": contact, "options": {"target": target}}));
        assert!(!err && c["overlapping_pairs"] == 0 && c["pairs"].as_array().unwrap().iter().all(|p| p["accuracy"] == 0.05), "{target}: {}", pretty(&c));
    }
}

/// The reference and std by section over MCP. (Which names find which
/// sections is docs.rs's names_agents_asked_for; here, what only the
/// server's path shows.)
#[test]
fn mcp_spec_and_std_answer_by_section() {
    let s = Scratch::new("mcp-docs");
    let mut srv = Server::start(ROOT, &s, &[]);
    let mut ask = |id, cmd: &str, opts: Value| {
        let (err, c) = srv.call(id, json!({"command": cmd, "options": opts}));
        assert!(!err, "{cmd} {opts}: {}", cut(text(&c)));
        text(&c).to_string()
    };
    let t = ask(6, "spec", json!({"section": "motion"}));
    assert!(t.starts_with("## Motion") && !t.contains("\n## "), "spec by section: {}", cut(&t));
    let t = ask(12, "std", json!({"contents": true}));
    let own = &t[..t.find("## built in").expect("std's contents end with the built-in forms")];
    let sigs: Vec<&str> = own.lines().filter(|l| l.starts_with('(')).collect();
    let mut files: Vec<_> = std::fs::read_dir(root("std")).unwrap().map(|e| e.unwrap().path()).filter(|p| p.extension().is_some_and(|e| e == "ink")).collect();
    files.sort();
    let std_src: String = files.iter().map(|f| std::fs::read_to_string(f).unwrap()).collect();
    let defs = std_src.matches("\n(fn ").count() + std_src.matches("\n(def ").count();
    assert!(
        sigs.len() == defs && sigs.contains(&"(def petg)") && t.len() * 10 < std_src.len(),
        "std's contents list every definition ({} of {defs}), materials too, in a tenth of its size ({} of {})",
        sigs.len(),
        t.len(),
        std_src.len()
    );
    // A built-in form's name finds its entries, not a name holding it
    // (std's htd-groove-polygon); std's contents list it by signature.
    let t = ask(13, "spec", json!({"section": "polygon"}));
    assert!(t.contains("(vec u v b r)") && !t.contains("htd"), "a built-in form by name: {}", cut(&t));
    assert!(ask(15, "std", json!({"contents": true})).contains("\n(polygon pts [drafts z0 z1])\n"), "std's contents list the built-in forms");
    let t = ask(14, "std", json!({"section": "lattices"}));
    assert!(t.contains("(fn gyroid"), "std by part of a section's banner gives the section: {}", cut(&t));
    let t = ask(18, "std", json!({"section": "belt"}));
    assert!(t.contains("(fn path-loop") && t.contains("below (belts and chains)"), "belts finds path tools and where the belts are: {}", cut(&t));
    let t = ask(32, "std", json!({"section": "Manufacturing checks"}));
    assert!(t.starts_with("From the reference"), "std's section named by a spec title: {}", cut(&t));
    // A name in a heading and its subsections picks the heading. (One too
    // long with its subsections gives its own text and their paths:
    // docs.rs's a_section_too_long_gives_its_subsections_paths, on a
    // reference made so; this one's sections may all fit.)
    let top = ask(11, "spec", json!({"section": "top-level"}));
    assert!(top.starts_with("## Top-level forms"), "{}", cut(&top));
}

/// Every section of the reference is one whole answer, never shortened
/// (agents don't read the rest from its file): within 24,000 characters
/// (mcp.rs MAX_TEXT). One longer with its subsections gives its own text
/// and their paths, so a section without them must fit by itself.
#[test]
fn every_spec_section_fits_an_mcp_answer() {
    let heads: Vec<String> = Ink::new(["spec", "--contents"]).text().lines().map(|l| l.trim_start_matches('#').trim().to_string()).filter(|h| !h.is_empty()).collect();
    let sizes = par(heads.len(), |i| Ink::new(["spec", "--section", &heads[i]]).text().len());
    let over: Vec<_> = heads.iter().zip(&sizes).filter(|(_, n)| **n > 24000).collect();
    assert!(heads.len() > 30 && over.is_empty(), "split a section that doesn't fit into subsections: {over:?} ({} sections)", heads.len());
}

/// A call that never ends (check reading a FIFO nobody writes) leaves the
/// server answering, and the client's cancellation kills it; one that
/// crashes is the server's error, not an answer.
#[test]
#[cfg(target_os = "linux")]
fn mcp_a_call_never_ending_is_cancelled_and_one_crashing_is_an_error() {
    let s = Scratch::new("mcp-cancel");
    let never = fifo(&s);
    let mut srv = Server::start(ROOT, &s, &[]);
    let check = json!({"name": "ink", "arguments": {"command": "check", "file": never}});
    srv.send(json!({"id": 8, "method": "tools/call", "params": check}));
    assert!(srv.rpc(9, "ping", json!({}))["id"] == 9, "answers while a call runs");
    // The server blocks SIGTERM and SIGINT (one thread takes them); a call's
    // process blocks none, so a signal sent to it ends it.
    // Read until clear: glibc blocks every signal for a moment while a
    // thread starts, so one read can land in that; the server's mask, if
    // inherited, would never clear.
    let kid = until(|| srv.calls_with(&never).first().copied()).expect("no process for the call");
    let blocked = || std::fs::read_to_string(format!("/proc/{kid}/status")).ok()?.lines().find_map(|l| l.strip_prefix("SigBlk:").map(|v| v.trim().to_string()));
    let clear = until(|| (blocked().as_deref() == Some("0000000000000000")).then_some(()));
    assert!(clear.is_some(), "a call's process runs with no signal blocked: {:?}", blocked());
    srv.send(json!({"method": "notifications/cancelled", "params": {"requestId": 8}}));
    let gone = until(|| srv.calls_with(&never).is_empty().then_some(()));
    assert!(gone.is_some() && srv.rpc(10, "ping", json!({}))["id"] == 10, "a cancelled call's process is killed, and the call not answered");
    // A command that crashes (here an abort): a host running ink for others
    // knows it was ink's failure, not the caller's.
    srv.send(json!({"id": 16, "method": "tools/call", "params": check}));
    let kid = until(|| srv.calls_with(&never).first().copied()).expect("no process for the call");
    signal("ABRT", kid);
    let crash = srv.next();
    assert!(
        crash["id"] == 16 && crash["error"]["code"] == -32603 && crash["error"]["message"].as_str().is_some_and(|m| m.contains("crashed (killed by signal 6)")),
        "a crashed command is a JSON-RPC error saying so: {}",
        pretty(&crash)
    );
}

/// The process running `file`'s call, once its progress log (INK_PROGRESS)
/// has a finished item, whole, and (`next`) the next one started: a stop
/// then never lands before what the call is to report is written.
#[cfg(target_os = "linux")]
fn first_done(srv: &Server, file: &str, next: bool) -> Option<u32> {
    let pid = *srv.calls_with(file).first()?;
    let env = std::fs::read(format!("/proc/{pid}/environ")).ok()?;
    let log = env.split(|&b| b == 0).find_map(|v| v.strip_prefix(b"INK_PROGRESS="))?;
    let text = std::fs::read_to_string(String::from_utf8_lossy(log).as_ref()).ok()?;
    let whole = text.rsplit_once('\n').map_or("", |(w, _)| w);
    let lines: Vec<Value> = whole.lines().filter_map(|l| serde_json::from_str(l).ok()).collect();
    let done = lines.iter().position(|l| l.get("done").is_some())?;
    (!next || lines[done..].iter().any(|l| l.get("running").is_some())).then_some(pid)
}

/// A call stopped at its time limit answers as the finished command would
/// (ok, marked stopped), with what finished as that answer gives it. Three
/// at once, each given `LIMIT` s: the FIFO's check, which reports its
/// elapsed time to the progress token it gave; two claims, the quick one
/// finished (its row, with its source and res) and the other (uncached)
/// running, which took most of the time, so the hint names it; and an
/// optimize (each design a sphere's volume, about 0.1 s), which gives its
/// best design so far and each one evaluated, not cached. A command's
/// process is stopped (SIGSTOP) once its progress log shows a first item
/// finished, in whole lines (the claims' next one started), so what
/// finished does not depend on how fast the machine is, only that a first
/// item takes under the limit. SIGTERM then ends the
/// server, its files removed.
#[test]
#[cfg(target_os = "linux")]
fn mcp_a_call_stopped_at_its_limit_gives_what_finished() {
    // Over PROGRESS_EVERY (2 s), so a progress note comes before the stop,
    // and long enough that a first item finishes under a full suite's load
    // (at 4 s, a trivial claim can miss it on a loaded 24-core machine).
    const LIMIT: u64 = 8;
    let s = Scratch::new("mcp-limit");
    let never = fifo(&s);
    let claims = s.model(
        "slow-claims.ink",
        "(import \"std\")\n(body ball (shape (move (sphere 4) (* 12 t) 0 0)))\n(assert small (within (sphere 4) (-5 -5 -5) (5 5 5)) (accuracy 0.1))\n\
         (assert inside (subset (shape-of ball) (- (move (place (sphere 4) (at (rot-z 30))) (* 12 t) 0 0) 0.01)) (accuracy 0.002))\n",
    );
    let opt = s.model(
        "slow-optimize.ink",
        "(import \"std\")\n(param r 3)\n(body ball (shape (sphere r)))\n(optimize biggest (vary r 1 20) (maximize (volume (shape-of ball)))\n  (require (within (shape-of ball) (-10 -10 -10) (10 10 10))) (budget 100))\n",
    );
    let mut srv = Server::start(ROOT, &s, &[("INK_NO_CACHE", "1")]);
    let call = |args: Value| json!({"name": "ink", "arguments": args});
    srv.send(json!({"id": 1, "method": "tools/call", "params": call(json!({"command": "assert", "file": claims, "options": {"t-span": "0:1", "res": 0.05}, "timeout_s": LIMIT}))}));
    srv.send(json!({"id": 2, "method": "tools/call", "params": call(json!({"command": "optimize", "file": opt, "timeout_s": LIMIT}))}));
    let mut p3 = call(json!({"command": "check", "file": never, "timeout_s": LIMIT}));
    p3["_meta"] = json!({"progressToken": "long"});
    srv.send(json!({"id": 3, "method": "tools/call", "params": p3}));
    for (f, next) in [(&claims, true), (&opt, false)] {
        if let Some(pid) = until(|| first_done(&srv, f, next)) {
            signal("STOP", pid);
        }
    }
    let (mut replies, mut notes) = (std::collections::BTreeMap::new(), vec![]);
    while replies.len() < 3 {
        let m = srv.next();
        match m["id"].as_u64() {
            Some(id) => _ = replies.insert(id, m["result"].clone()),
            None => notes.push(m),
        }
    }
    let secs: Vec<u64> =
        notes.iter().map(|n| n["params"]["message"].as_str().and_then(|m| m.split("running for ").nth(1)).and_then(|r| r.split(' ').next()?.parse().ok()).unwrap_or(0)).collect();
    let mut sorted = secs.clone();
    sorted.sort();
    assert!(
        !notes.is_empty()
            && notes.iter().all(|n| n["method"] == "notifications/progress" && n["params"]["progressToken"] == "long")
            && secs[0] > 0
            && secs == sorted
            && secs[secs.len() - 1] <= LIMIT,
        "a long call reports progress, its seconds as they pass: {notes:?}"
    );
    let t3 = replies[&3]["content"][0]["text"].as_str().unwrap_or("");
    let stopped = t3.split("stopped after ").nth(1).and_then(|r| r.split(' ').next()?.parse::<u64>().ok());
    assert!(
        replies[&3]["isError"] == true && stopped.is_some_and(|s| s >= LIMIT) && srv.calls_with(&never).is_empty(),
        "a call past its timeout_s is stopped, its process killed: {}",
        pretty(&replies[&3])
    );
    let answer = |t: &str| t.find("\n{").map_or(Value::Null, |i| serde_json::from_str(&t[i..]).unwrap_or(Value::Null));
    let t1 = replies[&1]["content"][0]["text"].as_str().unwrap_or("");
    let a = answer(t1);
    let (fin, row) = (&a["finished"], &a["asserts"][0]);
    assert!(
        replies[&1]["isError"] == true
            && a["ok"] == true
            && a["stopped"] == true
            && t1.contains("1 of 2 asserts done")
            && fin["counts"] == json!({"holds": 1})
            && fin["running_when_stopped"] == "inside"
            && fin["slowest"][0]["name"] == "small"
            && row["source"].as_str().is_some_and(|s| s.ends_with(":3"))
            && row.get("res").is_some()
            && row["accuracy"] == 0.1
            && t1.contains("Most of the time went to inside")
            && t1.contains("(name inside)")
            && t1.contains("Nothing is cached here"),
        "a call stopped at its limit gives what finished as the answer gives it, and what took the time: {}",
        &t1[..t1.floor_char_boundary(1200)]
    );
    {
        let t2 = replies[&2]["content"][0]["text"].as_str().unwrap_or("");
        let a = answer(t2);
        let designs = a["designs"].as_array().cloned().unwrap_or_default();
        let best = designs.iter().filter(|x| x["search_feasible"] == true).max_by(|x, y| num(&x["objective"]).total_cmp(&num(&y["objective"]))).map(|x| x["design"].clone());
        assert!(
            replies[&2]["isError"] == true
                && a["ok"] == true
                && a["stopped"] == true
                && a["optimize"] == "biggest"
                && !designs.is_empty()
                && designs.len() < 100
                && a["finished"]["total"] == 100
                && Some(&a["design"]) == best.as_ref()
                && a["apply"].as_str().is_some_and(|s| s.starts_with("--param r="))
                && designs.iter().all(|x| x["design"].get("r").is_some() && x.get("requires").is_some())
                && a["method"] == "derivative-free"
                && t2.contains("starts the search again"),
            "an optimize stopped at its limit gives its best design so far and each design evaluated, with its method; \
             its designs are not cached, so more time starts it again (it said it went on): {}",
            &t2[..t2.floor_char_boundary(1200)]
        );
    }
    let mine = srv.dir(&s);
    let there = mine.is_dir();
    signal("TERM", srv.pid());
    let st = srv.wait();
    assert!(there && st.code() == Some(143) && !mine.exists(), "SIGTERM ends it, its files removed: {st:?}, its directory there before: {there}");
}

/// A claim asked at one moment though it reads t, finished before the call
/// is stopped: the answer the stop gives says so, its row and the note,
/// as the whole answer would. Stopped as the test above stops one.
#[test]
#[cfg(target_os = "linux")]
fn mcp_a_call_stopped_keeps_what_was_asked_at_one_moment() {
    const LIMIT: u64 = 8;
    let s = Scratch::new("mcp-limit-moment");
    let claims = s.model(
        "moment-claims.ink",
        "(import \"std\")\n(assert instant (value t >= 0))\n\
         (assert heavy (subset (sphere 4) (- (place (sphere 4) (at (rot-z 30))) 0.01)) (accuracy 0.00001))\n",
    );
    let mut srv = Server::start(ROOT, &s, &[("INK_NO_CACHE", "1")]);
    srv.send(json!({"id": 1, "method": "tools/call", "params": {"name": "ink", "arguments": {"command": "assert", "file": claims, "timeout_s": LIMIT}}}));
    if let Some(pid) = until(|| first_done(&srv, &claims, true)) {
        signal("STOP", pid);
    }
    let reply = srv.next();
    let t = reply["result"]["content"][0]["text"].as_str().unwrap_or("");
    let a: Value = t.find("\n{").map_or(Value::Null, |i| serde_json::from_str(&t[i..]).unwrap_or(Value::Null));
    let instant = a["asserts"].as_array().and_then(|r| r.iter().find(|r| r["name"] == "instant")).cloned().unwrap_or_default();
    assert!(
        a["stopped"] == true
            && instant["status"] == "holds"
            && instant["one_moment"] == 0.0
            && a["note"].as_str().is_some_and(|n| n.contains("claim instant at t = 0"))
            // The claim left says what it was estimated to cost.
            && a["finished"]["left_estimated"][0]["name"] == "heavy"
            && a["finished"]["left_estimated"][0]["cost_estimate"]["basis"].as_str().is_some_and(|b| b.contains("a proof")),
        "{}",
        &t[..t.floor_char_boundary(1500)]
    );
}

/// SIGTERM while calls are starting their processes (as the first one
/// appears, 30 asked for) leaves none of them running.
#[test]
#[cfg(target_os = "linux")]
fn mcp_ended_while_calls_start_leaves_none_of_their_processes() {
    let s = Scratch::new("mcp-term");
    let never = fifo(&s);
    let mut srv = Server::unread(ROOT, &s);
    let lines: String = (0..30)
        .map(|i| format!("{}\n", json!({"jsonrpc": "2.0", "id": i, "method": "tools/call", "params": {"name": "ink", "arguments": {"command": "check", "file": never}}})))
        .collect();
    srv.write(&lines);
    let t = Instant::now();
    while !srv.has_children() {
        assert!(t.elapsed() < Duration::from_secs(60), "no call started in 60 s");
    }
    signal("TERM", srv.pid());
    srv.wait();
    let left = running_with(&never);
    for &p in &left {
        signal("KILL", p);
    }
    assert!(left.is_empty(), "processes left running: {left:?}");
}

#[test]
fn mcp_answers_every_request_it_read_before_its_input_closed() {
    // A script piping initialize, tools/list and a call closes its input
    // after the last: the call is answered, not dropped, before the server
    // exits 0.
    let s = Scratch::new("mcp-eof");
    let mut srv = Server::start(ROOT, &s, &[]);
    let call = json!({"name": "ink", "arguments": {"command": "assert", "file": fixture("sphere.ink")}});
    srv.write(
        &[
            json!({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {"protocolVersion": "2025-06-18"}}),
            json!({"jsonrpc": "2.0", "method": "notifications/initialized"}),
            json!({"jsonrpc": "2.0", "id": 2, "method": "tools/list"}),
            json!({"jsonrpc": "2.0", "id": 3, "method": "tools/call", "params": call}),
        ]
        .iter()
        .map(|m| format!("{m}\n"))
        .collect::<String>(),
    );
    let st = srv.close();
    let ids: Vec<Value> = (0..3).map(|_| srv.next()).map(|r| r["id"].clone()).collect();
    assert!(st.code() == Some(0) && ids == [json!(1), json!(2), json!(3)], "{st:?}: {ids:?}");
}

#[test]
fn mcp_an_error_comes_back_compact() {
    let s = Scratch::new("mcp-compact-error");
    let mut srv = Server::start(ROOT, &s, &[]);
    for (id, args) in [(1, json!({"command": "check", "file": "no-such-file.ink"})), (2, json!({"command": "check", "file": fixture("sphere.ink"), "options": {"bogus": 1}}))] {
        let (err, c) = srv.call(id, args);
        let t = text(&c);
        assert!(err && t.trim_end().lines().count() == 1, "an error is one line, as an answer is: {t}");
    }
}

/// A field the kernel rewrites (a moving probe's coordinates put into it, a
/// motion frame's fixed params) keeps its spatial arithmetic's meaning in
/// real numbers where its operands become constants: along x = y, 1e16 +
/// (x - y) + 0.5 is 1e16 + 0.5, not its float 1e16. Each case is unknown or right,
/// never a decided wrong answer; its ordinary twin keeps its exact answer,
/// and arithmetic written on constants and fixed params alone its float's,
/// as a def's. The command line and MCP answer alike.
#[test]
fn a_rewrite_keeps_spatial_arithmetic_real() {
    const BOX: &str = "(box 10 10 10)";
    let probe =
        |field: &str, at: &str, head: &str| format!("(import \"std\")\n{head}\n(def f (max {field} {BOX}))\n(assert inside (subset {at} f))\n(assert outside (disjoint {at} f))\n");
    let cancel = "(- (- (+ 1e16 (+ (- x y) 0.5)) 1e16) 0.25)";
    let span = "(t-span 0 1)";
    let frame = |radius: &str, head: &str, bodies: bool| {
        let disc = format!("(max (- (sqrt (+ (square x) (square y))) {radius}) (slab-z -1 1))");
        let ball = "(move (sphere 0.1) 1.25 0 0)";
        let parts = if bodies {
            format!(
                "(body disc (shape {disc}))\n(body ball (shape {ball}))\n(assembly root (instance a ball (at (rot-z (* 90 t)))) (instance b disc (at (rot-z (* 180 t)))))\n(def a (shape-of root/a))\n(def b (shape-of root/b))"
            )
        } else {
            format!("(def a (place {ball} (at (rot-z (* 90 t)))))\n(def b (place {disc} (at (rot-z (* 180 t)))))")
        };
        format!("(import \"std\")\n{head}\n(t-span 0 1)\n{parts}\n(assert apart (disjoint a b))\n(assert gap (clearance a b 0.01))\n")
    };
    let real_radius = "(+ 1 (- (+ 1e16 (+ (- x (+ x p)) 0.5)) 1e16))";
    // Its real value outside the box's corner region: the probe is outside
    // (inside fails, outside holds) or inside; the frame's ball overlaps the
    // disc (apart and gap fail).
    let (out, inn, overlap) = (&[("inside", "fails"), ("outside", "holds")][..], &[("inside", "holds"), ("outside", "fails")][..], &[("apart", "fails"), ("gap", "fails")][..]);
    // (name, model, extra options, the right answers, whether it must be exactly those)
    let cases: Vec<(&str, String, Vec<&str>, &[(&str, &str)], bool)> = vec![
        ("equal-moving", probe(cancel, "(probe t t 0)", span), vec![], out, false),
        ("one-moving", probe("(- (- (+ 1e16 (+ (- x t) 0.5)) 1e16) 0.25)", "(probe t 0 0)", span), vec![], out, false),
        ("moving-product", probe(cancel, "(probe (* 2 t) (* 2 t) 0)", span), vec![], out, false),
        ("ranged", probe(cancel, "(probe p p 0)", "(def p (ranged 0 1))"), vec![], out, false),
        ("tolerance", probe(cancel, "(probe p p 0)", "(param p 0.5 (tol -0.1 0.1))"), vec![], out, false),
        ("certified", probe(cancel, "(probe q q 0)", "(param p 0.5 (tol -0.1 0.1))\n(solve s (unknowns (q 0.5)) (== q p))"), vec![], out, false),
        ("placed", probe(cancel, "(place (probe 0 0 0) (at (move t t 0)))", span), vec![], out, false),
        ("computed-big", probe("(- (- (+ big (+ (- x y) 0.5)) big) 0.25)", "(probe t t 0)", "(t-span 0 1)\n(def big (* 1e100 1e100))"), vec![], out, false),
        ("at-moment", probe(cancel, "(probe t t 0)", span), vec!["--t", "0.5"], out, false),
        ("one-at-moment", probe("(- (- (+ 1e16 (+ (- x t) 0.5)) 1e16) 0.25)", "(probe t 0 0)", span), vec!["--t", "0"], out, false),
        ("fixed-coordinate", probe("(- (- (+ 1e16 (+ y 0.5)) 1e16) 0.25)", "(probe t q 0)", "(t-span 0 1)\n(param q 0)"), vec![], out, false),
        ("fixed-coordinate-twin", probe("(- (+ y 0.5) 0.25)", "(probe t q 0)", "(t-span 0 1)\n(param q 0)"), vec![], out, true),
        ("tiny", probe("(- (/ (sqrt (square (* s (+ (- x y) 2)))) s) 1)", "(probe t t 0)", "(t-span 0 1)\n(def s (* 1e-100 1e-100))"), vec![], out, false),
        ("huge", probe("(- (/ (sqrt (square (* s (+ (- x y) 0.5)))) s) 1)", "(probe t t 0)", "(t-span 0 1)\n(def s (* 1e100 1e100))"), vec![], inn, false),
        ("frame-fields", frame(real_radius, "(param p 0)", false), vec![], overlap, false),
        ("frame-bodies", frame(real_radius, "(param p 0)", true), vec![], overlap, false),
        // Twins: the same answers decided exactly.
        ("equal-twin", probe("(- (+ (- x y) 0.5) 0.25)", "(probe t t 0)", span), vec![], out, true),
        ("one-twin", probe("(- (+ (- x t) 0.5) 0.25)", "(probe t 0 0)", span), vec![], out, true),
        ("ranged-twin", probe("(- (+ (- x y) 0.5) 0.25)", "(probe p p 0)", "(def p (ranged 0 1))"), vec![], out, true),
        ("moment-twin", probe("(- (+ (- x y) 0.5) 0.25)", "(probe t t 0)", span), vec!["--t", "0.5"], out, true),
        ("tiny-twin", probe("(- (sqrt (square (+ (- x y) 2))) 1)", "(probe t t 0)", span), vec![], out, true),
        ("huge-twin", probe("(- (sqrt (square (+ (- x y) 0.5))) 1)", "(probe t t 0)", span), vec![], inn, true),
        ("frame-twin", frame("1.5", "(param p 0)", false), vec![], overlap, true),
        ("frame-bodies-twin", frame("1.5", "(param p 0)", true), vec![], overlap, true),
        // Written on constants and fixed params alone: a def's float, -0.25.
        ("fixed-written", probe("fixed", "(probe t t 0)", "(t-span 0 1)\n(def fixed (- (- (+ 1e16 0.5) 1e16) 0.25))"), vec![], inn, true),
        ("param-written", probe("fixed", "(probe t t 0)", "(t-span 0 1)\n(param big 1e16)\n(def fixed (- (- (+ big 0.5) big) 0.25))"), vec![], inn, true),
    ];
    let s = Scratch::new("mcp-subst-real");
    // The cases side by side, each with its own server: one after another,
    // they are this binary's longest test.
    let wrong: Vec<String> = par(cases.len(), |i| {
        let (name, text, opts, right, exact) = &cases[i];
        let s = &s;
        let tmp = s.join(format!("srv-{name}"));
        std::fs::create_dir_all(&tmp).unwrap();
        let mut srv = Server::start(ROOT, &tmp, &[("INK_NO_CACHE", "1")]);
        let (mut id, mut wrong) = (100, Vec::new());
        let file = s.model(&format!("{name}.ink"), text);
        let check = s.ink(["check", file.as_str()]).env("INK_NO_CACHE", "1").ok();
        // The bodies' disc warns that its bounds are loose (a performance
        // note); nothing else warns.
        let loose = |w: &Value| name.contains("bodies") && w.as_str().is_some_and(|w| w.contains("interval bounds of its field are"));
        assert!(check["warnings"].as_array().is_none_or(|ws| ws.iter().all(loose)), "{name}: {}", pretty(&check));
        for (claim, answer) in right.iter() {
            let mut args: Vec<&str> = vec!["assert", file.as_str(), "--name", claim];
            args.extend(opts.iter());
            let cli = verdict(&s.ink(args).env("INK_NO_CACHE", "1").any()["asserts"][0]);
            let mut o = serde_json::Map::new();
            o.insert("name".into(), json!(claim));
            for kv in opts.chunks(2) {
                o.insert(kv[0].trim_start_matches("--").into(), kv[1].parse::<f64>().map_or(json!(kv[1]), |v| json!(v)));
            }
            id += 1;
            let (err, mcp) = srv.answer(id, json!({"command": "assert", "file": file, "options": o}));
            let mcp = verdict(&mcp["asserts"][0]);
            assert!(!err && cli == mcp, "{name} {claim}: the command line says {cli}, MCP {mcp}");
            let ok = if *exact { cli == *answer } else { cli == *answer || cli == "unknown" || cli == format!("{answer} estimate") };
            if !ok {
                wrong.push(format!("{name} {claim}: {cli}, where it is {answer}{}", if *exact { "" } else { " or unknown" }));
            }
        }
        wrong
    })
    .concat();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
