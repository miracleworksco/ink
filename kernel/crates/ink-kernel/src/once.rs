//! A value computed once, outside any lock: the rule for what is shared
//! between threads and costs parallel work to make.
//!
//! A thread waiting in rayon (a `join`, a `par_iter`, a `scope`) runs other
//! queued work on its own stack meanwhile, any of it: another model's load,
//! another claim. Inside `OnceLock::get_or_init`'s initialiser or under a
//! `Mutex`, that work can ask for the same cell or lock on the thread that
//! holds it (a lock is not reentrant: it waits for itself), or a worker
//! blocked on it can hold the half the holder waits for. Either way every
//! thread ends up waiting and none runs. So nothing
//! that may wait in rayon runs under a lock or in an initialiser: it is
//! computed outside, then published (`get_or_compute`), and a lock guards
//! only the publishing.

use std::sync::OnceLock;

/// `cell`'s value, computed by `f` where it has none yet: `f` runs outside
/// the cell, so it may run rayon work (and anything that does). Threads
/// asking at once each compute it, and the first to finish publishes its
/// value, which every caller then reads: `f` must give the same answer
/// wherever it runs (a pure function of what the cell is for), and the
/// cost of a race is the work done twice, not a wait.
pub fn get_or_compute<T>(cell: &OnceLock<T>, f: impl FnOnce() -> T) -> &T {
    if let Some(v) = cell.get() {
        return v;
    }
    let v = f();
    cell.get_or_init(|| v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicBool, Ordering::SeqCst};
    use std::sync::mpsc;
    use std::time::Duration;

    /// `ask` on a thread of its own, its answer within `secs`, else None
    /// (deadlocked: the thread is left waiting).
    fn within<T: Send + 'static>(secs: u64, ask: impl FnOnce() -> T + Send + 'static) -> Option<T> {
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || _ = tx.send(ask()));
        rx.recv_timeout(Duration::from_secs(secs)).ok()
    }

    /// An initialiser waiting in rayon runs other work on its stack meanwhile; work that
    /// asks for the same cell there is asked on the thread computing it.
    /// Inside `get_or_init` it waits for itself (std: "deadlocks"); computed outside,
    /// each computes its own.
    #[test]
    fn work_run_inside_a_computation_may_ask_for_its_cell() {
        static CELL: OnceLock<u32> = OnceLock::new();
        let got = within(30, || *get_or_compute(&CELL, || get_or_compute(&CELL, || 7) + 1));
        assert_eq!(got, Some(7), "the inner computation's value is published first, and read by both");
    }

    /// The other shape: a computation waits in a join for its other half,
    /// which another worker took and which asks for the cell. Inside
    /// `get_or_init` that worker waits for the initialiser and the
    /// initialiser for it, in a pool of two, every time (the first half
    /// waits until the second has begun on the other worker). Computed
    /// outside, the other worker computes its own value.
    #[test]
    fn a_computation_may_wait_for_work_that_asks_for_its_cell() {
        fn ask(inside: bool) -> u32 {
            let cell: &'static OnceLock<u32> = Box::leak(Box::new(OnceLock::new()));
            let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();
            let begun: &'static AtomicBool = Box::leak(Box::new(AtomicBool::new(false)));
            let compute = move || {
                let (_, b) = rayon::join(
                    || {
                        while !begun.load(SeqCst) {
                            std::thread::yield_now();
                        }
                    },
                    || {
                        begun.store(true, SeqCst);
                        match inside {
                            true => *cell.get_or_init(|| 5),
                            false => *get_or_compute(cell, || 5),
                        }
                    },
                );
                b
            };
            pool.install(|| match inside {
                true => *cell.get_or_init(compute),
                false => *get_or_compute(cell, compute),
            })
        }
        assert_eq!(within(30, || ask(false)), Some(5));
        // The deadlock the rule avoids, every time (its threads are left
        // waiting), so the test is of what it says.
        assert_eq!(within(2, || ask(true)), None, "get_or_init with a join inside should deadlock in a pool of two");
    }
}
