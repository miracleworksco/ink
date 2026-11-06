//! The runners: the threads a search's pieces run on, at most as many at
//! once as
//! the pool has threads (RAYON_NUM_THREADS). A search hands its pieces over
//! together (`each`), and the runners take them one at a time, in turn
//! across the searches waiting: the first search's next piece, then the
//! next search's, so each search in flight gets its share and one alone
//! gets every runner. A piece runs to its end on its runner, its parallel
//! work there too (a pool of that thread alone), but for its cells'
//! rounds, which it hands to the runners as a search does its pieces
//! (search::Search::round): a thread waits only on its own work, and a
//! search alone keeps every runner busy. In one rayon pool a thread
//! waiting in a join runs whatever else is queued meanwhile, other
//! searches' pieces and their cells' rounds among it, and returns only
//! once that has ended; with searches side by side (ink interfere's moving
//! pairs) new work keeps arriving and the first searches' pieces never
//! return (on a large assembly a pair holding alone in 13 s sits in its
//! first 16 pieces for an hour). A search on a runner (a piece's round, or one
//! `spawn`ed) runs its own work while any is not begun, then gives way to
//! another runner while the rest end.

use std::any::Any;
use std::cell::Cell;
use std::collections::VecDeque;
use std::panic::{AssertUnwindSafe, catch_unwind, resume_unwind};
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, Condvar, Mutex, MutexGuard, OnceLock};

/// Pieces handed over together: piece i is `run(i)`; how many have begun
/// and ended (changed under the runners' lock); the first panic.
struct Batch {
    run: Box<dyn Fn(usize) + Send + Sync>,
    len: usize,
    begun: AtomicUsize,
    ended: AtomicUsize,
    panic: Mutex<Option<Box<dyn Any + Send>>>,
}

impl Batch {
    fn new(run: Box<dyn Fn(usize) + Send + Sync>, len: usize) -> Arc<Batch> {
        Arc::new(Batch { run, len, begun: AtomicUsize::new(0), ended: AtomicUsize::new(0), panic: Mutex::new(None) })
    }

    /// Piece i, its panic kept (the first; a later one's payload is
    /// leaked, not dropped, as its drop could panic): this never unwinds.
    fn run(&self, i: usize) {
        if let Err(p) = catch_unwind(AssertUnwindSafe(|| (self.run)(i))) {
            let mut first = self.panic.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
            match first.is_none() {
                true => *first = Some(p),
                false => std::mem::forget(p),
            }
        }
    }

    fn ended(&self) -> bool {
        self.ended.load(Relaxed) == self.len
    }
}

/// The runners' state: the batches with pieces not begun, in turn; the
/// threads running (at most `n`), the searches waiting to run again (they
/// go first), and the runner threads idle and starting.
#[derive(Default)]
struct State {
    queue: VecDeque<Arc<Batch>>,
    running: usize,
    wanting: usize,
    idle: usize,
    starting: usize,
}

impl State {
    /// The next piece in turn: the first batch's, which then goes last.
    fn next(&mut self) -> Option<(Arc<Batch>, usize)> {
        let b = self.queue.pop_front()?;
        let i = b.begun.fetch_add(1, Relaxed);
        if i + 1 < b.len {
            self.queue.push_back(b.clone());
        }
        Some((b, i))
    }

    /// b's next piece, if one is not begun.
    fn next_of(&mut self, b: &Arc<Batch>) -> Option<usize> {
        let i = b.begun.load(Relaxed);
        if i == b.len {
            return None;
        }
        b.begun.store(i + 1, Relaxed);
        if i + 1 == b.len {
            self.queue.retain(|x| !Arc::ptr_eq(x, b));
        }
        Some(i)
    }
}

struct Runners {
    n: usize,
    state: Mutex<State>,
    /// Idle runners wait here for a piece.
    work: Condvar,
    /// A search whose pieces have ended waits here to run again.
    resume: Condvar,
    /// A search waits here for its pieces to end.
    done: Condvar,
}

thread_local! {
    /// This thread is a runner.
    static RUNNER: Cell<bool> = const { Cell::new(false) };
}

/// The runners, as many at once as the first caller's pool has threads
/// (the process's: RAYON_NUM_THREADS).
fn runners() -> &'static Runners {
    static R: OnceLock<Runners> = OnceLock::new();
    R.get_or_init(|| Runners::new(rayon::current_num_threads()))
}

impl Runners {
    fn new(n: usize) -> Runners {
        Runners { n, state: Mutex::default(), work: Condvar::new(), resume: Condvar::new(), done: Condvar::new() }
    }

    fn lock(&self) -> MutexGuard<'_, State> {
        self.state.lock().unwrap()
    }

    /// What may run is handed out: a place among the `n` running, freed,
    /// to a search waiting to run again; else a piece waiting, with a place
    /// free, to an idle runner, or a new one where none is idle or starting.
    fn wake(&'static self, st: &mut State) {
        if st.wanting > 0 && st.running < self.n {
            self.resume.notify_one();
        }
        if st.running + st.wanting < self.n && !st.queue.is_empty() {
            if st.idle > 0 {
                self.work.notify_one();
            } else if st.starting == 0 {
                st.starting += 1;
                // No thread: `each` must not return while a piece may run.
                if let Err(e) = std::thread::Builder::new().name("ink-runner".into()).spawn(move || self.runner()) {
                    eprintln!("ink: no thread for a runner: {e}");
                    std::process::abort();
                }
            }
        }
    }

    /// A piece of b has ended (under the lock).
    fn end(&self, b: &Batch) {
        if b.ended.fetch_add(1, Relaxed) + 1 == b.len {
            self.done.notify_all();
        }
    }

    /// A runner: the next piece in turn while a place is free, else idle.
    fn runner(&'static self) {
        RUNNER.set(true);
        // A pool of this thread alone: a piece's parallel work runs here.
        let _pool = rayon::ThreadPoolBuilder::new().num_threads(1).use_current_thread().build().expect("a runner's pool");
        let mut st = self.lock();
        st.starting -= 1;
        loop {
            let next = if st.running + st.wanting < self.n { st.next() } else { None };
            match next {
                Some((b, i)) => {
                    st.running += 1;
                    self.wake(&mut st);
                    drop(st);
                    b.run(i);
                    st = self.lock();
                    st.running -= 1;
                    self.end(&b);
                    if st.wanting > 0 {
                        self.resume.notify_one();
                    }
                }
                None => {
                    st.idle += 1;
                    st = self.work.wait(st).unwrap();
                    st.idle -= 1;
                }
            }
        }
    }
}

/// This thread is a runner: what it runs is a search's piece, or a search
/// run on one (`spawn`).
pub fn here() -> bool {
    RUNNER.get()
}

/// The runners that run at once.
pub fn count() -> usize {
    runners().n
}

/// `f` of each item, run by the runners, in order. A panic in one is
/// raised here once all have ended.
pub fn each<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    each_on(runners(), items, f)
}

/// `each`, on the runners `r`.
fn each_on<T: Sync, R: Send>(r: &'static Runners, items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    if items.is_empty() {
        return Vec::new();
    }
    let out: Vec<Mutex<Option<R>>> = items.iter().map(|_| Mutex::new(None)).collect();
    let run = |i: usize| {
        let r = f(&items[i]);
        *out[i].lock().unwrap() = Some(r);
    };
    // SAFETY: `run` borrows this frame, which outlives every call of it:
    // this returns only once each piece has ended (`ended`, under the lock
    // the runners count it under), and no piece begins after all have.
    // Nothing between here and there unwinds (pieces' panics are caught);
    // were it to, the process aborts (`Unwinding`) before the frame goes.
    let run = unsafe { std::mem::transmute::<&(dyn Fn(usize) + Sync), &'static (dyn Fn(usize) + Sync)>(&run) };
    let unwinding = Unwinding;
    let b = Batch::new(Box::new(run), items.len());
    let mut st = r.lock();
    st.queue.push_back(b.clone());
    r.wake(&mut st);
    // On a runner: its own pieces while any is not begun, then its place
    // goes to another until the rest end.
    let mut gave = false;
    if RUNNER.get() {
        while let Some(i) = st.next_of(&b) {
            drop(st);
            b.run(i);
            st = r.lock();
            r.end(&b);
        }
        if !b.ended() {
            (st.running, gave) = (st.running - 1, true);
            r.wake(&mut st);
        }
    }
    while !b.ended() {
        st = r.done.wait(st).unwrap();
    }
    std::mem::forget(unwinding);
    if gave {
        st.wanting += 1;
        while st.running >= r.n {
            st = r.resume.wait(st).unwrap();
        }
        (st.wanting, st.running) = (st.wanting - 1, st.running + 1);
    }
    drop(st);
    if let Some(p) = b.panic.lock().unwrap().take() {
        resume_unwind(p);
    }
    out.into_iter().map(|x| x.into_inner().unwrap().expect("each piece ran")).collect()
}

/// Aborts the process if dropped: unwinding out of `each` while its pieces
/// may still run.
struct Unwinding;

impl Drop for Unwinding {
    fn drop(&mut self) {
        std::process::abort();
    }
}

/// `job` run on a runner, beside the caller, in turn with the searches'
/// pieces (ink interfere's moving pairs, each a search). A panic in it is
/// its own to report.
pub fn spawn(job: impl FnOnce() + Send + 'static) {
    let job = Mutex::new(Some(job));
    let b = Batch::new(
        Box::new(move |_| {
            let job = job.lock().unwrap().take();
            job.into_iter().for_each(|j| j());
        }),
        1,
    );
    let r = runners();
    let mut st = r.lock();
    st.queue.push_back(b);
    r.wake(&mut st);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicBool;
    use std::thread::sleep;
    use std::time::Duration;

    thread_local! {
        /// The search whose piece this thread is running.
        static IN: Cell<Option<u8>> = const { Cell::new(None) };
    }

    /// Search s's piece, `f`: whether it ran on a thread inside another's.
    fn piece(s: u8, f: impl FnOnce()) -> bool {
        let outer = IN.replace(Some(s));
        f();
        IN.set(outer);
        outer.is_some_and(|o| o != s)
    }

    fn ms(k: u64) -> Duration {
        Duration::from_millis(k)
    }

    /// A thread waiting inside a search's piece runs none of another's
    /// work: the first search's pieces (fewer than the threads) each wait
    /// in a join for a half another thread takes, while others hand over
    /// short pieces again and again. In one rayon pool the waiting
    /// threads would run the second's pieces inside the first's (and, with
    /// searches' pieces arriving without end, never return to them: ink
    /// interfere's moving pairs on a large assembly).
    #[test]
    fn a_waiting_piece_runs_only_its_own_work() {
        let n = runners().n;
        let crossed = AtomicBool::new(false);
        std::thread::scope(|s| {
            s.spawn(|| each(&vec![(); n.div_ceil(2)], |_| crossed.fetch_or(piece(0, || _ = rayon::join(|| sleep(ms(20)), || sleep(ms(40)))), Relaxed)));
            for _ in 0..3 {
                s.spawn(|| {
                    for _ in 0..30 {
                        each(&vec![(); 2 * n], |_| crossed.fetch_or(piece(1, || sleep(ms(1))), Relaxed));
                    }
                });
            }
        });
        assert!(!crossed.load(Relaxed), "a piece ran inside another search's");
    }

    /// Searches take turns: one whose pieces were handed over first, and
    /// hold every runner, does not keep the runners from another's: as
    /// runners come free they take the first's next piece, then the
    /// second's, in turn, and the first's alone once the second's have all
    /// begun; none begins while every runner is busy. In one rayon pool a
    /// thread would take its own search's pieces first, and the second's
    /// begin once the first's had run out. On runners of its own, freed one at a
    /// time, the order pieces begin in is the order they are handed out,
    /// whatever else the machine runs.
    #[test]
    fn searches_take_turns() {
        let n = 3;
        let r: &'static Runners = Box::leak(Box::new(Runners::new(n)));
        // Each piece notes its search as it begins, then ends when handed
        // a turn (given out one at a time).
        let (begun, began) = (Mutex::new(Vec::new()), Condvar::new());
        let (turns, turn) = (Mutex::new(0), Condvar::new());
        let piece = |s: u8| {
            begun.lock().unwrap().push(s);
            began.notify_all();
            let mut t = turn.wait_while(turns.lock().unwrap(), |t| *t == 0).unwrap();
            *t -= 1;
        };
        let until = |k: usize| drop(began.wait_while(begun.lock().unwrap(), |b| b.len() < k).unwrap());
        let give = |k: usize| {
            *turns.lock().unwrap() += k;
            turn.notify_all();
        };
        std::thread::scope(|s| {
            s.spawn(|| each_on(r, &vec![(); 4 * n], |_| piece(0)));
            until(n);
            s.spawn(|| each_on(r, &vec![(); n], |_| piece(1)));
            while r.lock().queue.len() < 2 {
                std::thread::yield_now();
            }
            assert_eq!(begun.lock().unwrap().len(), n, "a piece began while every runner was busy");
            for k in n..5 * n {
                give(1);
                until(k + 1);
            }
            give(n);
        });
        let want: Vec<u8> = [vec![0; n], [0, 1].repeat(n), vec![0; 2 * n]].concat();
        assert_eq!(begun.into_inner().unwrap(), want, "the order the searches' pieces began in");
    }
}
