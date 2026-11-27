//! What the command-line tests share: running ext's ink, a scratch
//! directory of each test's own, and models written into it. (A copy of
//! the kernel's crates/ink-cli/tests/cli/common.rs, cut to what these
//! tests use: ext builds on the kernel as a third party does, its tests
//! included.)
#![allow(dead_code)]

use serde_json::Value;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// The two roots: ROOT, the folder holding the crate under test (its
/// fixtures, examples, benchmarks and docs), and KERNEL, the kernel's.
#[path = "roots.rs"]
mod roots;
/// The binary under test, as the test crate names it.
pub use crate::INK;
pub use roots::{KERNEL, ROOT};
/// Threads per ink call (RAYON_NUM_THREADS), unless a test sets its own.
pub const THREADS: &str = "4";

/// A path under ROOT.
pub fn root(rel: &str) -> String {
    format!("{ROOT}/{rel}")
}

/// A fixture in tests/fixtures.
pub fn fixture(name: &str) -> String {
    root(&format!("tests/fixtures/{name}"))
}

/// A fixture of the kernel's tests/fixtures that these tests read too.
pub fn kernel_fixture(name: &str) -> String {
    format!("{KERNEL}/tests/fixtures/{name}")
}

/// A number in an answer (NaN where it is none).
pub fn num(v: &Value) -> f64 {
    v.as_f64().unwrap_or(f64::NAN)
}

/// INK_NO_CACHE is set: checks of what the cache keeps are skipped.
pub fn no_cache() -> bool {
    std::env::var_os("INK_NO_CACHE").is_some()
}

/// Whether a refusal says the binary was built without a package a call
/// needs (`the sim package, not in this build`).
pub fn not_built(text: &str) -> bool {
    text.split("the ").any(|s| s.split_once(" package, not in this build").is_some_and(|(p, _)| !p.is_empty() && !p.contains(' ')))
}

/// What one ink call printed and how it exited.
#[derive(Debug)]
pub struct Run {
    pub stdout: String,
    pub stderr: String,
    pub code: Option<i32>,
}

impl Run {
    /// Its output as JSON: the call's answer.
    pub fn json(&self) -> Value {
        serde_json::from_str(&self.stdout).unwrap_or_else(|e| panic!("non-JSON output ({e}):\n{}\n{}", self.stdout, self.stderr))
    }
}

/// A call's JSON answer, ok or not, its exit status agreeing (0 ok, 1 an
/// error, 3 an assert with a claim that does not hold).
fn answered(what: &str, r: Run) -> Value {
    let d = r.json();
    let ok = d["ok"].as_bool() == Some(true);
    let want = if !ok {
        1
    } else if d["all_yes"] == Value::Bool(false) {
        3
    } else {
        0
    };
    assert_eq!(r.code, Some(want), "ink {what} -> ok={ok}, all_yes {}, exited {:?}:\n{}", d["all_yes"], r.code, r.stderr);
    d
}

/// The most CPU seconds a call spends once its searches are to stop
/// (`Ink::cpu_past`): answering what they decided, and the rounds of work
/// under way when the stop came (a round of studies.rs's slow claim reads
/// 1500 spheres at each point: 2.9 CPU s past its stop, once in 12 under
/// load). One whose searches went on would spend its threads' worth for as
/// long as they ran: that claim, about 240 CPU s.
pub const PAST_CPU_S: f64 = 10.0;

/// An ink call to make: its arguments, working directory and environment.
pub struct Ink {
    cmd: Command,
    args: Vec<String>,
}

impl Ink {
    pub fn new<S: AsRef<std::ffi::OsStr>>(args: impl IntoIterator<Item = S>) -> Ink {
        let mut cmd = Command::new(INK);
        let args: Vec<_> = args.into_iter().map(|a| a.as_ref().to_string_lossy().into_owned()).collect();
        cmd.args(&args).current_dir(ROOT).env("RAYON_NUM_THREADS", THREADS);
        Ink { cmd, args }
    }
    pub fn cwd(mut self, dir: impl AsRef<Path>) -> Ink {
        self.cmd.current_dir(dir);
        self
    }
    pub fn env(mut self, k: &str, v: impl AsRef<std::ffi::OsStr>) -> Ink {
        self.cmd.env(k, v);
        self
    }
    pub fn env_remove(mut self, k: &str) -> Ink {
        self.cmd.env_remove(k);
        self
    }
    /// Run it with this text on its standard input.
    pub fn stdin(self, text: &str) -> Run {
        self.go(Some(text))
    }
    pub fn run(self) -> Run {
        self.go(None)
    }
    fn go(mut self, input: Option<&str>) -> Run {
        use std::io::Write;
        self.cmd.stdin(if input.is_some() { Stdio::piped() } else { Stdio::null() }).stdout(Stdio::piped()).stderr(Stdio::piped());
        let mut child = self.cmd.spawn().unwrap_or_else(|e| panic!("cannot run {INK}: {e}"));
        if let Some(text) = input {
            child.stdin.take().unwrap().write_all(text.as_bytes()).unwrap();
        }
        let out = child.wait_with_output().unwrap();
        Run { stdout: String::from_utf8_lossy(&out.stdout).into_owned(), stderr: String::from_utf8_lossy(&out.stderr).into_owned(), code: out.status.code() }
    }
    /// Its JSON answer (as `any`), and the CPU seconds its process spent
    /// once `stop` had passed: the work it did after its searches were to
    /// stop. That a call stops promptly is measured so, not by the time it
    /// takes, which a machine busy with other work stretches (a machine under
    /// load: 29.5 s for a call whose searches stopped at 17 s). The time
    /// from `stop` to the reading of its CPU is not counted: a call that
    /// went on would spend far more after it.
    #[cfg(target_os = "linux")]
    pub fn cpu_past(mut self, stop: std::time::SystemTime) -> (Value, f64) {
        let what = self.args.join(" ");
        self.cmd.stdin(Stdio::null()).stdout(Stdio::piped()).stderr(Stdio::piped());
        let mut child = self.cmd.spawn().unwrap_or_else(|e| panic!("cannot run {INK}: {e}"));
        let read = |mut r: Box<dyn std::io::Read + Send>| {
            std::thread::spawn(move || {
                let mut s = String::new();
                _ = r.read_to_string(&mut s);
                s
            })
        };
        let (out, err) = (read(Box::new(child.stdout.take().unwrap())), read(Box::new(child.stderr.take().unwrap())));
        // Its state and CPU (utime and stime, in USER_HZ: 100 a second on
        // Linux), from /proc: there until it is reaped, its exit unreaped
        // until it is read so.
        let stat = || {
            let s = std::fs::read_to_string(format!("/proc/{}/stat", child.id())).ok()?;
            let f: Vec<&str> = s.rsplit_once(')')?.1.split_whitespace().collect();
            Some((f[0] == "Z", (f[11].parse::<f64>().ok()? + f[12].parse::<f64>().ok()?) / 100.0))
        };
        if let Ok(d) = stop.duration_since(std::time::SystemTime::now()) {
            std::thread::sleep(d);
        }
        let (_, before) = stat().expect("the call's process");
        let after = loop {
            match stat() {
                Some((true, cpu)) => break cpu,
                Some(_) => std::thread::sleep(std::time::Duration::from_millis(5)),
                None => panic!("the call's process was reaped unread"),
            }
        };
        let code = child.wait().unwrap().code();
        let r = Run { stdout: out.join().unwrap(), stderr: err.join().unwrap(), code };
        (answered(&what, r), after - before)
    }
    /// Elsewhere than Linux, its answer, and no CPU measured.
    #[cfg(not(target_os = "linux"))]
    pub fn cpu_past(self, _stop: std::time::SystemTime) -> (Value, f64) {
        (self.any(), 0.0)
    }
    /// Its JSON answer, ok or not, its exit status agreeing (0 ok, 1 an
    /// error, 3 an assert with a claim that does not hold).
    pub fn any(self) -> Value {
        let what = self.args.join(" ");
        answered(&what, self.run())
    }
    /// Its JSON answer, which must be ok.
    pub fn ok(self) -> Value {
        let what = self.args.join(" ");
        let d = self.any();
        assert!(d["ok"].as_bool() == Some(true), "ink {what} -> not ok:\n{}", pretty(&d));
        d
    }
    /// Its JSON answer, which must be an error (not a package missing from this build).
    pub fn err(self) -> Value {
        let what = self.args.join(" ");
        let d = self.any();
        assert!(d["ok"].as_bool() == Some(false), "ink {what} -> ok, expected an error:\n{}", pretty(&d));
        d
    }
    /// Its output, which must exit 0.
    pub fn text(self) -> String {
        let what = self.args.join(" ");
        let r = self.run();
        assert_eq!(r.code, Some(0), "ink {what} exited {:?}:\n{}\n{}", r.code, r.stdout, r.stderr);
        r.stdout
    }
}

/// `ink ARGS...`: its answer, which must be ok.
pub fn ink<S: AsRef<std::ffi::OsStr>>(args: impl IntoIterator<Item = S>) -> Value {
    Ink::new(args).ok()
}

/// `ink ARGS...`: its answer, which must be an error; returns the answer.
pub fn ink_err<S: AsRef<std::ffi::OsStr>>(args: impl IntoIterator<Item = S>) -> Value {
    Ink::new(args).err()
}

/// A {columns, rows} table (eval's results, a slice's loops, a solve's probes, rays) as records:
/// an object per row, a dotted column nested ("d_dparams.r").
pub fn records(t: &Value) -> Vec<Value> {
    let cols: Vec<&str> = t["columns"].as_array().unwrap_or_else(|| panic!("not a table: {t}")).iter().map(|c| c.as_str().unwrap()).collect();
    let rows = t["rows"].as_array().unwrap_or_else(|| panic!("not a table: {t}"));
    rows.iter()
        .map(|row| {
            let mut r = serde_json::Map::new();
            for (c, v) in cols.iter().zip(row.as_array().unwrap()) {
                let mut path: Vec<&str> = c.split('.').collect();
                let last = path.pop().unwrap();
                let mut o = &mut r;
                for u in path {
                    o = o.entry(u).or_insert_with(|| Value::Object(Default::default())).as_object_mut().unwrap();
                }
                o.insert(last.to_string(), v.clone());
            }
            Value::Object(r)
        })
        .collect()
}

/// JSON, indented and cut to a readable length, for a failure's message.
pub fn pretty(v: &Value) -> String {
    let s = serde_json::to_string_pretty(v).unwrap();
    if s.len() > 4000 { format!("{}...", &s[..s.floor_char_boundary(4000)]) } else { s }
}

/// A claim's or pair's answer as tests name it: its status, and ` estimate`
/// after one that holds or fails with no guaranteed accuracy (null).
pub fn verdict(row: &Value) -> String {
    let st = match row["status"].as_str() {
        Some(s) => s,
        None => match row["verdict"].as_str() {
            Some("yes") => "holds",
            Some("no") => "fails",
            Some(_) => "unknown",
            None => "?",
        },
    };
    if matches!(st, "holds" | "fails") && row.get("accuracy").is_none_or(Value::is_null) { format!("{st} estimate") } else { st.to_string() }
}

/// Each claim of an `ink require` answer: name -> verdict.
pub fn statuses(d: &Value) -> std::collections::BTreeMap<String, String> {
    d["requirements"].as_array().unwrap().iter().map(|a| (a["name"].as_str().unwrap().to_string(), verdict(a))).collect()
}

/// A directory of a test's own, empty when made, removed when dropped (a
/// test that panics too). ink's own temporary files go in it (TMPDIR).
pub struct Scratch(PathBuf);

impl Scratch {
    pub fn new(tag: &str) -> Scratch {
        // Numbered, so two tests with one tag never share a directory.
        static MADE: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = MADE.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("ext-cli-{tag}-{}-{n}", std::process::id()));
        _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Scratch(dir)
    }
    /// A file of this text in it; its path.
    pub fn model(&self, name: &str, text: &str) -> String {
        let p = self.0.join(name);
        if let Some(d) = p.parent() {
            std::fs::create_dir_all(d).unwrap();
        }
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }
    /// A path in it (not created).
    pub fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
    /// An ink call run with TMPDIR here.
    pub fn ink<S: AsRef<std::ffi::OsStr>>(&self, args: impl IntoIterator<Item = S>) -> Ink {
        Ink::new(args).env("TMPDIR", &self.0)
    }
}

impl std::ops::Deref for Scratch {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<Path> for Scratch {
    fn as_ref(&self) -> &Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A GLB file: its JSON chunk and its binary body.
pub struct Glb {
    pub json: Value,
    pub body: Vec<u8>,
}

impl Glb {
    pub fn read(path: &str) -> Glb {
        let data = std::fs::read(path).unwrap();
        let u = |o: usize| u32::from_le_bytes(data[o..o + 4].try_into().unwrap()) as usize;
        let jlen = u(12);
        let json = serde_json::from_slice(&data[20..20 + jlen]).unwrap();
        let blen = u(20 + jlen);
        Glb { json, body: data[28 + jlen..28 + jlen + blen].to_vec() }
    }
    /// Accessor `a`'s elements, each its floats (SCALAR, VEC3, VEC4).
    pub fn floats(&self, a: usize) -> Vec<Vec<f32>> {
        let acc = &self.json["accessors"][a];
        let view = &self.json["bufferViews"][acc["bufferView"].as_u64().unwrap() as usize];
        let w = match acc["type"].as_str().unwrap() {
            "SCALAR" => 1,
            "VEC2" => 2,
            "VEC3" => 3,
            "VEC4" => 4,
            t => panic!("accessor type {t}"),
        };
        let off = view["byteOffset"].as_u64().unwrap_or(0) as usize + acc["byteOffset"].as_u64().unwrap_or(0) as usize;
        let n = acc["count"].as_u64().unwrap() as usize;
        (0..n).map(|i| (0..w).map(|k| f32::from_le_bytes(self.body[off + 4 * (i * w + k)..][..4].try_into().unwrap())).collect()).collect()
    }
    /// Accessor `a`'s elements as unsigned integers (indices).
    pub fn indices(&self, a: usize) -> Vec<u32> {
        let acc = &self.json["accessors"][a];
        let view = &self.json["bufferViews"][acc["bufferView"].as_u64().unwrap() as usize];
        let off = view["byteOffset"].as_u64().unwrap_or(0) as usize + acc["byteOffset"].as_u64().unwrap_or(0) as usize;
        let n = acc["count"].as_u64().unwrap() as usize;
        match acc["componentType"].as_u64() {
            Some(5125) => (0..n).map(|i| u32::from_le_bytes(self.body[off + 4 * i..][..4].try_into().unwrap())).collect(),
            Some(5123) => (0..n).map(|i| u16::from_le_bytes(self.body[off + 2 * i..][..2].try_into().unwrap()) as u32).collect(),
            t => panic!("index component type {t:?}"),
        }
    }
    /// Each mesh named PREFIX...: its shells (triangles joined through
    /// vertices at one position, rounded to 1e-4).
    pub fn shells(&self, prefix: &str) -> Vec<usize> {
        let mut out = vec![];
        for m in self.json["meshes"].as_array().unwrap().iter().filter(|m| m["name"].as_str().is_some_and(|n| n.starts_with(prefix))) {
            let mut parent = std::collections::HashMap::<[i64; 3], [i64; 3]>::new();
            fn find(p: &mut std::collections::HashMap<[i64; 3], [i64; 3]>, x: [i64; 3]) -> [i64; 3] {
                let mut x = x;
                loop {
                    let up = *p.entry(x).or_insert(x);
                    if up == x {
                        return x;
                    }
                    x = up;
                }
            }
            for pr in m["primitives"].as_array().unwrap() {
                let pos: Vec<[i64; 3]> =
                    self.floats(pr["attributes"]["POSITION"].as_u64().unwrap() as usize).iter().map(|v| [0, 1, 2].map(|k| (v[k] as f64 * 1e4).round() as i64)).collect();
                let ix = self.indices(pr["indices"].as_u64().unwrap() as usize);
                for t in ix.chunks(3) {
                    let a = find(&mut parent, pos[t[0] as usize]);
                    let b = find(&mut parent, pos[t[1] as usize]);
                    parent.insert(b, a);
                    let c = find(&mut parent, pos[t[2] as usize]);
                    parent.insert(c, a);
                }
            }
            let keys: Vec<_> = parent.keys().copied().collect();
            let roots: std::collections::HashSet<_> = keys.into_iter().map(|k| find(&mut parent, k)).collect();
            out.push(roots.len());
        }
        out
    }
}
