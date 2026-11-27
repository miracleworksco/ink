//! What the integration tests share.
// A file whose tests all need a package this build lacks uses none of it.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// A directory of a test's own under the system's temporary directory
/// (`ink-TAG-PID`), empty when made and removed with what is in it when
/// dropped, a test that panics too: tests leave nothing behind.
pub struct Scratch(PathBuf);

impl Scratch {
    pub fn new(tag: &str) -> Scratch {
        let dir = std::env::temp_dir().join(format!("ink-{tag}-{}", std::process::id()));
        _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Scratch(dir)
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

/// A part handed from a rest to a turning grip to a post by selects on t
/// (`turns` turns in the grip), beside a wall it never reaches: each entry
/// of its placement is a select, the grip's nested in the rest's.
pub fn handed(turns: usize) -> String {
    let grip: String = (1..=turns).map(|k| format!(" (rot-{} (* {} (/ t 1s)))", if k % 2 == 1 { 'x' } else { 'z' }, 10 + k)).collect();
    format!(
        "(import \"std\")
(t-span 0 1s)
(body part (shape (min (sphere 2) (move (sphere 1) 3 0 0))))
(body wall (shape (move (box 2 40 40) -20 0 0)))
(assembly root
  (instance part part (select (- (/ t 1s) 0.3) (at (move 10 0 0)) (select (- (/ t 1s) 0.6) (at{grip} (move 10 0 0)) (at (move 0 10 0)))))
  (instance wall wall))
(require apart (disjoint (shape-of root/part) (shape-of root/wall)))"
    )
}

/// The model `src` viewed over its motion (its file's t-span), from a
/// file in `dir`.
pub fn over_motion(dir: &Path, src: &str) -> ink::view::ModelView {
    let f = dir.join("m.ink");
    std::fs::write(&f, src).unwrap();
    let mut v = ink::view::ModelView::load_in(ink::registry(), &f, &[], 0.0).unwrap_or_else(|e| panic!("{e:?}"));
    let span = v.doc.t_span;
    v.set_motion(span).unwrap();
    v
}

/// This test binary run again for its test `name` alone, `flag` set (the
/// test then does its work and nothing else), on 4 threads: whether it
/// passed, and the CPU seconds its process spent (utime and stime, from
/// /proc as it exits: USER_HZ is 100 a second on Linux). A cost measured so
/// is the work done, not the time a machine busy with other work takes.
/// Passed means it ran that one test and it passed: a `name` that no
/// longer matches runs none, and is not passed.
#[cfg(target_os = "linux")]
pub fn alone(name: &str, flag: &str) -> (bool, f64) {
    use std::process::{Command, Stdio};
    let mut c = Command::new(std::env::current_exe().unwrap())
        .args(["--exact", name, "--include-ignored", "--test-threads", "1"])
        .env(flag, "1")
        .env("RAYON_NUM_THREADS", "4")
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let mut out = c.stdout.take().unwrap();
    let read = std::thread::spawn(move || {
        let mut s = String::new();
        _ = std::io::Read::read_to_string(&mut out, &mut s);
        s
    });
    let stat = || {
        let s = std::fs::read_to_string(format!("/proc/{}/stat", c.id())).ok()?;
        let f: Vec<&str> = s.rsplit_once(')')?.1.split_whitespace().collect();
        Some((f[0] == "Z", (f[11].parse::<f64>().ok()? + f[12].parse::<f64>().ok()?) / 100.0))
    };
    let cpu = loop {
        match stat() {
            Some((true, cpu)) => break cpu,
            Some(_) => std::thread::sleep(std::time::Duration::from_millis(10)),
            None => panic!("{name}: its process was reaped unread"),
        }
    };
    let ok = c.wait().unwrap().success();
    (ok && read.join().unwrap().contains("test result: ok. 1 passed"), cpu)
}

/// Elsewhere than Linux: run alone, its CPU not measured.
#[cfg(not(target_os = "linux"))]
pub fn alone(name: &str, flag: &str) -> (bool, f64) {
    let out = std::process::Command::new(std::env::current_exe().unwrap()).args(["--exact", name, "--include-ignored", "--test-threads", "1"]).env(flag, "1").output().unwrap();
    (out.status.success() && String::from_utf8_lossy(&out.stdout).contains("test result: ok. 1 passed"), 0.0)
}
