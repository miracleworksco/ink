//! The call's time limit (`ink mcp` tells a command its limit:
//! INK_DEADLINE). Past it every search stops, its cells left undecided
//! (query::Ended::Stopped), and no answer is cached: it would hold what the
//! stop left undecided.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

/// When searches stop, Unix seconds as f64 bits (0: never).
static AT: AtomicU64 = AtomicU64::new(0);
/// The time limit has passed (kept: every search after it stops too).
static PASSED: AtomicBool = AtomicBool::new(false);

/// Stop the searches at `at`, Unix seconds.
pub fn at(at: f64) {
    AT.store(at.to_bits(), Ordering::Relaxed);
}

/// Whether the searches are to stop: the time limit has passed.
pub fn reached() -> bool {
    if PASSED.load(Ordering::Relaxed) {
        return true;
    }
    let at = f64::from_bits(AT.load(Ordering::Relaxed));
    let passed = at > 0.0 && std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).is_ok_and(|d| d.as_secs_f64() >= at);
    if passed {
        PASSED.store(true, Ordering::Relaxed);
    }
    passed
}

/// `f` with the searches stopped once `share` of the call's time left has
/// passed (none: no time limit, so none sooner), and whether that earlier
/// stop is what stopped them: past it, the call's own limit has not passed,
/// so later searches go on (a claim over a motion leaves the rest of its
/// time to measuring what the search left: ink-lang's motion::ranged).
/// Searches beside `f` stop with it.
pub fn within<T>(share: f64, f: impl FnOnce() -> T) -> (T, bool) {
    let limit = f64::from_bits(AT.load(Ordering::Relaxed));
    let now = || std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0.0, |d| d.as_secs_f64());
    if !(limit > 0.0) || reached() {
        return (f(), false);
    }
    AT.store((now() + share * (limit - now()).max(0.0)).to_bits(), Ordering::Relaxed);
    let out = f();
    AT.store(limit.to_bits(), Ordering::Relaxed);
    let own = PASSED.load(Ordering::Relaxed) && now() < limit;
    if own {
        PASSED.store(false, Ordering::Relaxed);
    }
    (out, own)
}
