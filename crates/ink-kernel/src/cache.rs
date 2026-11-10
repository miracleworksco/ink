//! Results cached on disk by content: `$XDG_CACHE_HOME/ink/<kind>/<key>` (or
//! ~/.cache/ink/...), one deflated JSON file per key. A key (`key`) is the
//! cache's schema, the hash of the running executable (`code`), and the
//! fingerprint of everything the result was computed from: the expanded
//! inputs, the options, and the keys of the upstream results it used. The
//! executable's bytes hold all the code that computed a result: every
//! crate, the dependencies as resolved, the features, the compiler and its
//! flags, and the maths library (the `libm` crate: ink_core::graph::defined).
//! Keys are 128 bits, so two never meet by chance. So a rebuilt binary or a
//! change to the inputs misses and a stale result is never served; two
//! programs (the CLI and a viewer) keep their own entries. Disabled with
//! INK_NO_CACHE=1. The whole cache, every kind together, is kept near
//! 8 GB on disk, least recently used first out, so a kind no build writes
//! any more goes too. Best effort: any failure just
//! means recomputing. Results computed for trial designs (`trial`: the
//! designs an optimizer or a central difference visits) are read but not
//! written.
//!
//! The kernel's own answers (`memo`: proofs, integrals, searches) are kept
//! by what was asked, exactly as the kernel sees it: the tapes, the context
//! and the options. So a claim whose question is unchanged is answered from
//! the cache, whatever else in the file changed, while the binary is the
//! same.

use ink_core::V3;
use ink_core::eval::Ctx;
use ink_core::graph::Fingerprint;
use ink_core::tape::{Names, Tape};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

fn dir(kind: &str) -> Option<PathBuf> {
    if std::env::var_os("INK_NO_CACHE").is_some() {
        return None;
    }
    code()?;
    let base = std::env::var_os("XDG_CACHE_HOME").map(PathBuf::from).or_else(|| std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".cache")))?;
    Some(base.join("ink").join(kind))
}

/// The layout of keys and entries: bumped when it changes.
pub const SCHEMA: u64 = 2;

/// The key of a result: the cache schema, the code that computed it
/// (`code`), and the fingerprint of its inputs (including the keys of
/// upstream results it used).
pub fn key(inputs: u128) -> u128 {
    let mut fp = ink_core::graph::Fingerprint::default();
    // No code (it could not be read): nothing is loaded or stored (`dir`).
    fp.word(SCHEMA).key(code().unwrap_or(0)).key(inputs);
    fp.value()
}

/// The hash of the running executable's bytes, read once per process: the
/// identity of all the code that computes a result. On Linux the file the
/// process runs (/proc/self/exe), even if its path was rebuilt since. None
/// where it cannot be read: then nothing is cached.
pub fn code() -> Option<u128> {
    static CODE: std::sync::OnceLock<Option<u128>> = std::sync::OnceLock::new();
    *CODE.get_or_init(|| {
        let path = if cfg!(target_os = "linux") { PathBuf::from("/proc/self/exe") } else { std::env::current_exe().ok()? };
        std::fs::read(path).ok().map(|b| hash(&b))
    })
}

/// A 128-bit hash of `bytes`, at memory speed (a 33 MB binary in a few
/// milliseconds; FNV's byte at a time takes ten times as long on one
/// thread): each 8-byte word mixed with its position into two independent
/// sums, over blocks on threads of its own, then the length and the tail.
/// Not rayon's: `code` runs it once for all threads, so a rayon worker
/// waiting on it would wait for its own blocks (every worker of a busy
/// pool deadlocked so). Not for adversaries: two builds differ in many
/// words.
fn hash(bytes: &[u8]) -> u128 {
    use ink_core::graph::mix64;
    let sum = |first: usize, words: &[u8]| {
        let (mut a, mut b) = (0u64, 0u64);
        for (i, w) in words.as_chunks::<8>().0.iter().enumerate() {
            let w = u64::from_le_bytes(*w);
            let i = (first + i) as u64;
            a = a.wrapping_add(mix64(w ^ i.wrapping_mul(0x9e3779b97f4a7c15)));
            b = b.wrapping_add(mix64(w.rotate_left(32) ^ i.wrapping_mul(0xc2b2ae3d27d4eb4f) ^ 0x2545f4914f6cdd1d));
        }
        (a, b)
    };
    let whole = bytes.len() / 8 * 8;
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get()).min(8);
    let block = (whole / 8).div_ceil(threads).max(1) * 8;
    let (a, b) = std::thread::scope(|s| {
        let parts: Vec<_> = bytes[..whole].chunks(block).enumerate().map(|(k, words)| s.spawn(move || sum(k * block / 8, words))).collect();
        parts.into_iter().map(|p| p.join().unwrap()).fold((0u64, 0u64), |x, y| (x.0.wrapping_add(y.0), x.1.wrapping_add(y.1)))
    });
    let mut tail = [0u8; 8];
    tail[..bytes.len() - whole].copy_from_slice(&bytes[whole..]);
    let mut fp = Fingerprint::default();
    fp.word(bytes.len() as u64).word(u64::from_le_bytes(tail)).word(a).word(b);
    fp.value()
}

pub fn load<T: serde::de::DeserializeOwned>(kind: &str, key: u128) -> Option<T> {
    let path = dir(kind)?.join(format!("{key:032x}.json.z"));
    let bytes = std::fs::read(&path).ok()?;
    // Mark it used (eviction goes by modification time).
    if let Ok(f) = std::fs::File::options().append(true).open(&path) {
        _ = f.set_modified(std::time::SystemTime::now());
    }
    let json = miniz_oxide::inflate::decompress_to_vec(&bytes).ok()?;
    serde_json::from_slice(&json).ok()
}

pub fn store<T: serde::Serialize>(kind: &str, key: u128, value: &T) {
    if TRIALS.load(Ordering::Relaxed) > 0 || crate::stop::reached() {
        return;
    }
    let Some(dir) = dir(kind) else { return };
    let Ok(json) = serde_json::to_vec(value) else { return };
    let bytes = miniz_oxide::deflate::compress_to_vec(&json, 3);
    // Write then rename, so concurrent readers never see a partial file
    // (each write its own file: threads may store one key at once).
    static WRITES: AtomicUsize = AtomicUsize::new(0);
    let tmp = dir.join(format!("{key:032x}.{}.{}.tmp", std::process::id(), WRITES.fetch_add(1, Ordering::Relaxed)));
    if std::fs::create_dir_all(&dir).is_ok() && std::fs::write(&tmp, &bytes).is_ok() {
        _ = std::fs::rename(&tmp, dir.join(format!("{key:032x}.json.z")));
        if let Some(root) = dir.parent() {
            prune(root, bytes.len() as u64, BUDGET);
        }
    }
}

/// The bytes the whole cache keeps on disk.
const BUDGET: u64 = 8192 * 1024 * 1024;

/// Trials running in this process (`trial`).
static TRIALS: AtomicUsize = AtomicUsize::new(0);

/// Run `f` on a trial design: nothing is stored while it runs, on any
/// thread (a trial's queries run on rayon's workers, which a flag of the
/// calling thread did not reach). A search's designs are points of a
/// continuous space that no later question asks for again, and an elastic
/// solve's entry is tens of megabytes: storing each would cost a second of
/// serialising per design and evict, within one run, the entries later
/// questions do ask for.
pub fn trial<R>(f: impl FnOnce() -> R) -> R {
    // Ends the trial on the way out, unwinding included.
    struct End;
    impl Drop for End {
        fn drop(&mut self) {
            TRIALS.fetch_sub(1, Ordering::Relaxed);
        }
    }
    TRIALS.fetch_add(1, Ordering::Relaxed);
    let _end = End;
    f()
}

/// Work computed in this process where a cached result could have served
/// it (not found, not cacheable, or caching off): kernel answers (`answer`)
/// and whatever else is cached its own way (a study's solve, a sampled
/// grid, a mesh, a wall's thickness), each of which says so (`computed`).
/// What a claim recomputed is the misses during it.
pub static MISSES: AtomicUsize = AtomicUsize::new(0);

/// Results read from a cache in this process (`served`). A claim's answer
/// came from the cache where it read at least one and recomputed nothing:
/// a claim that computes without a cache (a param's value) was not.
pub static HITS: AtomicUsize = AtomicUsize::new(0);

/// Count a result read from a cache (HITS).
pub fn served() {
    HITS.fetch_add(1, Ordering::Relaxed);
}

thread_local! {
    /// This thread's misses (as MISSES counts them), and whether the last
    /// answer `answer` gave on it was found (`here`).
    static HERE: std::cell::Cell<(usize, bool)> = const { std::cell::Cell::new((0, false)) };
}

/// This thread's misses so far, and whether the last answer `answer` gave
/// on it was found: what one question of several asked side by side
/// recomputed, which MISSES, counting every thread's, does not tell.
pub fn here() -> (usize, bool) {
    HERE.get()
}

/// Count a result just computed in place of a cached one (MISSES): every
/// cache besides `answer`'s calls it after it computes, or a claim that
/// solved a study would be reported cached.
pub fn computed() {
    MISSES.fetch_add(1, Ordering::Relaxed);
    HERE.set((HERE.get().0 + 1, false));
}

/// The key of a kernel answer: the query's name,
/// and what `inputs` hashes (everything the answer depends on: its tapes
/// with `Tape::fingerprint`, its context with `Ctx::fingerprint`, every
/// option; `ask_region` for a query over a region). None: caching is off, or `inputs` says the question cannot be
/// hashed (false).
pub fn ask(query: &str, inputs: impl FnOnce(&mut Fingerprint) -> bool) -> Option<u128> {
    dir(ANSWERS)?;
    let mut fp = Fingerprint::default();
    fp.str(query);
    inputs(&mut fp).then(|| key(fp.value()))
}

const ANSWERS: &str = "answers";

/// The answer under `key` if it is cached, else `compute`d and stored. An
/// answer is stored only where reading it back gives it exactly (JSON has
/// no infinities or NaN): a hit is what computing it again would give.
pub fn answer<T: serde::Serialize + serde::de::DeserializeOwned + std::fmt::Debug>(key: Option<u128>, compute: impl FnOnce() -> T) -> T {
    if let Some(v) = key.and_then(|k| load::<T>(ANSWERS, k)) {
        served();
        HERE.set((HERE.get().0, true));
        return v;
    }
    let v = compute();
    computed();
    if let Some(k) = key
        && let Ok(json) = serde_json::to_vec(&v)
        && serde_json::from_slice::<T>(&json).is_ok_and(|back| format!("{back:?}") == format!("{v:?}"))
    {
        store(ANSWERS, k, &v);
    }
    v
}

/// A query over the region `lo..hi` of `tapes`, cached: keyed by the
/// tapes, the region, the context as they read it (`Ctx::project`) and
/// what `inputs` adds (the options). `compute` is given that context alone,
/// so its answer cannot depend on what the key leaves out.
pub fn memo_region<T: serde::Serialize + serde::de::DeserializeOwned + std::fmt::Debug>(
    query: &str,
    tapes: &[&Tape],
    ctx: &Ctx,
    (lo, hi): (V3, V3),
    inputs: impl FnOnce(&mut Fingerprint),
    compute: impl FnOnce(&Ctx) -> T,
) -> T {
    let (key, ctx) = ask_region(query, tapes, ctx, (lo, hi), inputs);
    answer(key, || compute(&ctx))
}

/// `memo_region`'s answer where it is cached already, and nothing
/// computed where it is not: for what is worth reading only for free (a
/// cost estimate reading ink check's survey).
pub fn cached_region<T: serde::de::DeserializeOwned>(query: &str, tapes: &[&Tape], ctx: &Ctx, (lo, hi): (V3, V3), inputs: impl FnOnce(&mut Fingerprint)) -> Option<T> {
    load(ANSWERS, ask_region(query, tapes, ctx, (lo, hi), inputs).0?)
}

/// `memo_region`'s key (None: not cached), and the context its answer is
/// computed in, cached or not: the projected one (`Ctx::project`).
pub fn ask_region(query: &str, tapes: &[&Tape], ctx: &Ctx, (lo, hi): (V3, V3), inputs: impl FnOnce(&mut Fingerprint)) -> (Option<u128>, Ctx) {
    let ctx = ctx.project(tapes);
    let key = ask(query, |fp| {
        lo.iter().chain(&hi).for_each(|v| _ = fp.f64(*v));
        inputs(fp);
        question(fp, tapes, &ctx).is_some()
    });
    (key, ctx)
}

/// Hash a question's tapes and the context they read (projected:
/// `Ctx::project`) into `fp`, their tags and params numbered in the
/// question (`Names`, which it gives for the caller's own inputs that name
/// one). None: it cannot be hashed (a derivation that is not data:
/// Derive::fingerprint).
pub fn question(fp: &mut Fingerprint, tapes: &[&Tape], ctx: &Ctx) -> Option<Names> {
    let names = Names::of(tapes, ctx);
    fp.word(tapes.len() as u64);
    tapes.iter().for_each(|t| t.fingerprint(fp, &names));
    ctx.fingerprint(fp, &names).then_some(names)
}

/// After a store of `wrote` bytes under `root` (the cache: a directory per
/// kind), delete the least recently used entries of every kind beyond
/// `budget` bytes on disk (a small entry takes a 4 KiB block: a million of
/// them several times their length). Reading every entry's metadata takes
/// a second over a million entries, so it is done once a tenth of the
/// budget has been written since the last time, by any process: each store
/// appends a byte per 4 KiB it wrote to `root/written`, which pruning
/// empties. One pass prunes at a time, in any thread of any process: the
/// one holding the lock on that file (the others skip it: passes over one
/// listing, each deleting to the budget counting its own deletions alone,
/// would together empty the cache). It reads the count again once it holds
/// the lock, since a pass that just ended has emptied it.
fn prune(root: &std::path::Path, wrote: u64, budget: u64) {
    use std::io::Write;
    let Ok(mut f) = std::fs::File::options().create(true).append(true).open(root.join("written")) else { return };
    let due = |f: &std::fs::File| f.metadata().is_ok_and(|m| m.len() * 4096 > budget / 10);
    if f.write_all(&vec![0; wrote.div_ceil(4096) as usize]).is_err() || !due(&f) || f.try_lock().is_err() || !due(&f) {
        return;
    }
    _ = f.set_len(0);
    let Ok(kinds) = std::fs::read_dir(root) else { return };
    let mut files: Vec<(std::time::SystemTime, u64, PathBuf)> = kinds
        .flatten()
        .filter(|k| k.file_type().is_ok_and(|t| t.is_dir()))
        .filter_map(|k| std::fs::read_dir(k.path()).ok())
        .flatten()
        .flatten()
        .filter_map(|e| {
            let m = e.metadata().ok()?;
            #[cfg(unix)]
            let size = std::os::unix::fs::MetadataExt::blocks(&m) * 512;
            #[cfg(not(unix))]
            let size = m.len();
            Some((m.modified().ok()?, size, e.path()))
        })
        .collect();
    let mut total: u64 = files.iter().map(|f| f.1).sum();
    if total <= budget {
        return;
    }
    files.sort_by_key(|f| f.0);
    for (_, len, path) in files {
        if total <= budget * 9 / 10 {
            break;
        }
        if std::fs::remove_file(&path).is_ok() {
            total -= len;
        }
    }
}

#[cfg(test)]
mod trial_tests {
    use super::{TRIALS, trial};
    use std::sync::atomic::Ordering;

    #[test]
    fn a_trial_scope_nests_reaches_workers_and_ends_on_unwinding() {
        let on = || TRIALS.load(Ordering::Relaxed) > 0;
        assert!(!on());
        trial(|| {
            assert!(on());
            trial(|| assert!(on()));
            assert!(on(), "an inner scope leaves the outer one a trial");
            // A trial's work on rayon's workers is in it too.
            assert!(rayon::join(on, on) == (true, true));
        });
        assert!(!on());
        let caught = std::panic::catch_unwind(|| trial(|| panic!("a failed design")));
        assert!(caught.is_err());
        assert!(!on(), "a panic inside a trial leaves writes enabled");
        trial(|| _ = std::panic::catch_unwind(|| trial(|| panic!("inner"))));
        assert!(!on());
    }
}

#[cfg(test)]
mod tests {
    use super::key;
    use std::path::PathBuf;

    /// A directory at `tag` (a name, then perhaps a path in it) in the
    /// temporary directory; the directory named is removed when dropped.
    struct Fixture(PathBuf, PathBuf);

    impl Fixture {
        fn new(tag: &str) -> Fixture {
            let named = |t: &str| std::env::temp_dir().join(format!("ink-cache-{}-{t}", std::process::id()));
            let root = named(tag.split('/').next().unwrap());
            _ = std::fs::remove_dir_all(&root);
            Fixture(named(tag), root)
        }
    }

    impl Drop for Fixture {
        fn drop(&mut self) {
            _ = std::fs::remove_dir_all(&self.1);
        }
    }

    #[test]
    fn the_cache_is_pruned_once_a_tenth_of_its_budget_is_written() {
        // 64 entries of a block, over a budget of 16 blocks: a store of a
        // block (under a tenth of it) leaves them, the second takes the
        // oldest out to 90% of it, of every kind: a kind no build writes
        // any more goes as the others' entries come. The kind's directory
        // holds its entries alone (a count of them is what
        // crates/ink/tests/cache.rs reads).
        let f = Fixture::new("prune/answers");
        std::fs::create_dir_all(&f.0).unwrap();
        let t0 = std::time::SystemTime::now() - std::time::Duration::from_secs(1000);
        for i in 0..64u64 {
            let p = f.0.join(format!("{i:02}.json.z"));
            std::fs::write(&p, [1; 4096]).unwrap();
            std::fs::File::options().append(true).open(&p).unwrap().set_modified(t0 + std::time::Duration::from_secs(i)).unwrap();
        }
        let left = || std::fs::read_dir(&f.0).unwrap().count();
        // Another kind's entries, older than all of them.
        let old = f.1.join("reads");
        std::fs::create_dir_all(&old).unwrap();
        for i in 0..4u64 {
            let p = old.join(format!("{i}.json.z"));
            std::fs::write(&p, [1; 4096]).unwrap();
            std::fs::File::options().append(true).open(&p).unwrap().set_modified(t0 - std::time::Duration::from_secs(10)).unwrap();
        }
        let budget = 16 * 4096;
        super::prune(&f.1, 1, budget);
        assert_eq!(left(), 64, "a block written is under a tenth of the budget");
        super::prune(&f.1, 1, budget);
        assert_eq!(left(), 14, "the oldest go");
        assert!(f.0.join("63.json.z").exists() && !f.0.join("49.json.z").exists());
        assert_eq!(std::fs::read_dir(&old).unwrap().count(), 0, "the other kind's, older, went first");
        super::prune(&f.1, 1, budget);
        assert_eq!(left(), 14, "the count starts again");
    }

    #[cfg(unix)]
    #[test]
    fn concurrent_prunes_leave_the_cache_at_its_budget() {
        // 512 entries of a block and 4 of 64, over a budget of 256 blocks,
        // pruned by 8 threads at once (as by a test suite's many ink processes),
        // each having written past a tenth of it: one prunes, to 90% of
        // the budget. Passes over overlapping listings each deleted to the
        // budget counting their own deletions alone, and left nothing.
        let f = Fixture::new("concurrent/answers");
        std::fs::create_dir_all(&f.0).unwrap();
        let t0 = std::time::SystemTime::now() - std::time::Duration::from_secs(1000);
        for i in 0..516u64 {
            let p = f.0.join(format!("{i:03}.json.z"));
            std::fs::write(&p, vec![1; if i % 129 == 0 { 64 * 4096 } else { 4096 }]).unwrap();
            std::fs::File::options().append(true).open(&p).unwrap().set_modified(t0 + std::time::Duration::from_secs(i)).unwrap();
        }
        let budget = 256 * 4096;
        let start = std::sync::Barrier::new(8);
        std::thread::scope(|s| {
            for _ in 0..8 {
                s.spawn(|| {
                    start.wait();
                    super::prune(&f.1, budget / 8, budget);
                });
            }
        });
        let on_disk = |e: std::fs::DirEntry| std::os::unix::fs::MetadataExt::blocks(&e.metadata().unwrap()) * 512;
        let total: u64 = std::fs::read_dir(&f.0).unwrap().map(|e| on_disk(e.unwrap())).sum();
        assert!(total <= budget * 9 / 10 && total > budget * 9 / 10 - 64 * 4096, "{} blocks left of a budget of 256", total / 4096);
    }

    #[test]
    fn a_proof_key_holds_all_it_reads() {
        // Every input of a proof changes its key (cache::ask): the tape, t,
        // the motion's span, a ranged param's span, band and value, a param
        // it reads, the region, each option; a param it does not read, ranged
        // or not, does not.
        use crate::ranged::ProofOptions;
        use ink_core::eval::{Ctx, ParamSpans};
        use ink_core::graph::{Graph, Var};
        use ink_core::tape::Tape;
        let tape = |k: f64| {
            let mut g = Graph::new();
            let (x, t, p) = (g.var(Var::X), g.var(Var::T), g.var(Var::Param(0)));
            let c = g.constant(k);
            let xt = g.sub(x, t);
            let a = g.mul(xt, c);
            let f = g.sub(a, p);
            Tape::new(&g, f)
        };
        let ctx = Ctx { p_spans: ParamSpans::new([(0, (0.5, 1.5))]), ..Ctx::new(0.25, vec![1.0, 7.0]) };
        let opts = ProofOptions { res: None, finest: 0.01, depth: 8, tol: 0.05 };
        let k = |t: &Tape, c: &Ctx, lo: f64, o: ProofOptions| super::ask_region("prove", &[t], c, ([lo; 3], [1.0; 3]), |fp| o.fingerprint(fp)).0;
        let (t0, base) = (tape(2.0), k(&tape(2.0), &ctx, -1.0, opts));
        if base.is_none() {
            return; // caching is off
        }
        let with = |f: &dyn Fn(&mut Ctx)| {
            let mut c = ctx.clone();
            f(&mut c);
            k(&t0, &c, -1.0, opts)
        };
        let differ = [
            ("the tape", k(&tape(3.0), &ctx, -1.0, opts)),
            ("t", with(&|c| c.t = 0.5)),
            ("the motion's span", with(&|c| c.t_span = Some((0.0, 1.0)))),
            ("a ranged param's span", with(&|c| c.p_spans.set_span(0, (0.5, 1.0)))),
            ("a ranged param's band", with(&|c| c.p_spans = ParamSpans::new([(0, (0.5, 1.6))]))),
            ("a ranged param's value", with(&|c| c.params[0] = 1.25)),
            ("a param ranged", with(&|c| c.p_spans = ParamSpans::NONE)),
            ("the region", k(&t0, &ctx, -2.0, opts)),
            ("res", k(&t0, &ctx, -1.0, ProofOptions { res: Some(0.1), ..opts })),
            ("finest", k(&t0, &ctx, -1.0, ProofOptions { finest: 0.02, ..opts })),
            ("depth", k(&t0, &ctx, -1.0, ProofOptions { depth: 4, ..opts })),
            ("tol", k(&t0, &ctx, -1.0, ProofOptions { tol: 0.02, ..opts })),
        ];
        for (what, key) in differ {
            assert_ne!(key, base, "{what}");
        }
        assert_eq!(with(&|c| c.params[1] = 8.0), base, "a param the tape does not read");
        assert_eq!(with(&|c| c.p_spans = ParamSpans::new([(0, (0.5, 1.5)), (1, (6.0, 8.0))])), base, "a ranged param the tape does not read");
        // which the proof is not given: it is computed in the context its
        // key hashes
        let (_, given) = super::ask_region("prove", &[&t0], &ctx, ([-1.0; 3], [1.0; 3]), |_| {});
        assert_eq!((given.params.len(), given.params[0], given.t, &given.p_spans), (1, 1.0, 0.25, &ctx.p_spans));
        // The same expression built another way, in a graph holding more,
        // is the same tape.
        let mut g = Graph::new();
        let (y, p, t) = (g.var(Var::Y), g.var(Var::Param(0)), g.var(Var::T));
        let c = g.constant(2.0);
        let x = g.var(Var::X);
        let _unrelated = g.mul(y, c);
        let xt = g.sub(x, t);
        let a = g.mul(c, xt);
        let f = g.sub(a, p);
        assert_eq!(k(&Tape::new(&g, f), &ctx, -1.0, opts), base);
    }

    #[test]
    fn a_key_holds_names_by_their_place_in_the_question() {
        // A part imported into two models: its tags and params numbered
        // after what each declared before it. Its questions' keys are the
        // same in both, but where a ranged param (which a witness names by
        // the model's number) or the order of the part's names differs.
        use ink_core::eval::{Ctx, ParamSpans};
        use ink_core::graph::{Graph, Var};
        use ink_core::tape::Tape;
        // The part, after `before` tags and params of the model's own.
        let part = |before: u32, swap: bool| {
            let mut g = Graph::new();
            for i in 0..before {
                let y = g.var(Var::Y);
                g.tag(&format!("own{i}"), y);
            }
            let (x, p, q) = (g.var(Var::X), g.var(Var::Param(before)), g.var(Var::Param(before + 1)));
            let (p, q) = if swap { (q, p) } else { (p, q) };
            let a = g.sub(x, p);
            let a = g.tag("rim", a);
            let b = g.add(x, q);
            let b = g.tag("hub", b);
            let f = g.max(a, b);
            Tape::new(&g, f)
        };
        let ctx = |before: u32, ranged: Option<u32>| {
            let mut params = vec![9.0; before as usize];
            params.extend([1.0, 2.0]);
            Ctx { p_spans: ParamSpans::new(ranged.map(|k| (k, (0.5, 1.5)))), ..Ctx::new(0.0, params) }
        };
        let k = |t: &Tape, c: &Ctx| super::ask_region("prove", &[t], c, ([-1.0; 3], [1.0; 3]), |_| {}).0;
        let base = k(&part(0, false), &ctx(0, None));
        if base.is_none() {
            return; // caching is off
        }
        assert_eq!(k(&part(3, false), &ctx(3, None)), base, "the same part after three tags and params of the model's");
        assert_ne!(k(&part(3, true), &ctx(3, None)), base, "its params read the other way round");
        let ranged = k(&part(0, false), &ctx(0, Some(0)));
        assert_ne!(ranged, base, "a param ranged");
        assert_ne!(k(&part(3, false), &ctx(3, Some(3))), ranged, "a ranged param, which a witness names, numbered otherwise");
        assert_eq!(k(&part(0, false), &ctx(0, Some(0))), ranged);
        // It does not read t: its key is the same at every moment of a
        // motion; a tape reading t keys each moment (a_proof_key_holds_all_it_reads).
        let at = |t: f64, span: Option<(f64, f64)>| Ctx { t, t_span: span, ..ctx(0, None) };
        assert_eq!(k(&part(0, false), &at(2.5, Some((0.0, 8.0)))), base, "a moment and a span the part does not read");
    }

    #[test]
    fn a_grid_key_holds_its_samples() {
        // A grid's samples are fixed once made (a mutation would not
        // compile: GridData), and each of them is in the key of an answer
        // on it.
        use ink_core::eval::Ctx;
        use ink_core::graph::{Graph, Var};
        use ink_core::grid::GridData;
        use ink_core::tape::Tape;
        let key = |d: GridData| {
            let mut g = Graph::new();
            let at = [g.var(Var::X), g.var(Var::Y), g.var(Var::Z)];
            let f = g.grid(std::sync::Arc::new(d), at);
            super::ask_region("prove", &[&Tape::new(&g, f)], &Ctx::new(0.0, vec![]), ([0.0; 3], [1.0; 3]), |_| {}).0
        };
        let base = key(GridData::new([0.0; 3], 1.0, [2; 3], vec![1.0; 8]));
        if base.is_none() {
            return; // caching is off
        }
        let mut one = vec![1.0; 8];
        one[7] = -1.0;
        for (what, d) in [
            ("a value", GridData::new([0.0; 3], 1.0, [2; 3], one)),
            ("the origin", GridData::new([0.0, 0.0, 0.5], 1.0, [2; 3], vec![1.0; 8])),
            ("the spacing", GridData::new([0.0; 3], 0.5, [2; 3], vec![1.0; 8])),
            ("the dims", GridData::new([0.0; 3], 1.0, [4, 2, 1], vec![1.0; 8])),
        ] {
            assert_ne!(key(d), base, "{what}");
        }
        assert_eq!(key(GridData::new([0.0; 3], 1.0, [2; 3], vec![1.0; 8])), base);
    }

    #[test]
    fn keys_hold_the_code_and_the_inputs() {
        assert_ne!(key(1), key(2));
        assert!(super::code().is_some(), "the test binary is readable");
        // any word, its place, or the length changes the code's hash
        let bytes: Vec<u8> = (0..1000u32).map(|i| (i * 7) as u8).collect();
        let h = super::hash(&bytes);
        let mut other = bytes.clone();
        other[500] ^= 1;
        assert_ne!(super::hash(&other), h);
        other = bytes.clone();
        other.swap(0, 8);
        assert_ne!(super::hash(&other), h);
        assert_ne!(super::hash(&bytes[..999]), h);
        assert_eq!(super::hash(&bytes), h);
    }

    /// `code` hashes once for every thread, on whichever asks first: here
    /// a thread outside rayon's pool, while every worker waits for it. A
    /// hash on the pool's workers deadlocked so (a viewer's tests did).
    #[test]
    fn the_code_is_hashed_while_every_rayon_worker_waits_for_it() {
        use std::sync::{Arc, OnceLock, mpsc};
        let bytes = vec![7u8; 1 << 22];
        let once: Arc<OnceLock<u128>> = Arc::default();
        let (done, waited) = mpsc::channel();
        let workers = rayon::current_num_threads();
        let h = *once.get_or_init(|| {
            for _ in 0..workers {
                let (once, done) = (once.clone(), done.clone());
                rayon::spawn(move || done.send(*once.wait()).unwrap());
            }
            std::thread::sleep(std::time::Duration::from_millis(100));
            super::hash(&bytes)
        });
        assert!((0..workers).all(|_| waited.recv() == Ok(h)));
        assert_eq!(h, super::hash(&bytes));
    }
}
