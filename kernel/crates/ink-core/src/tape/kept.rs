//! Tapes derived from tapes (`Tape::narrowed`, `Tape::kept_cone`), kept
//! while they fit: a query's boxes and cells reach the same few tapes again
//! and again (a box's halves simplify as it did; a field's operands' cones
//! are read in every box its surface crosses), and building one costs more
//! than evaluating it there.
//!
//! Each is kept by the tape it was derived from (by its allocation, held
//! weakly, which keeps the address from being reused while the entry is, so
//! a key never names another tape) and what was asked of it, exactly: a
//! kept tape is the one that would be built. Tapes alike in every list
//! (`Tape::same`) are kept once and shared, so a region reached by two
//! paths (choices that differ only where they change nothing, or halves of
//! two boxes alike) is one tape, and what is derived from it is derived
//! once.
//!
//! At most `BYTES` in all (`Tape::bytes`, and the keys), none of more than
//! an eighth of it. Each entry is
//! marked when reached; one over the cap drops the oldest unmarked,
//! unmarking those it passes (a clock: what is reached again stays), and a
//! tape goes with the last entry holding it. Entries whose tape derived
//! from is gone are never reached again, and go when the hand comes to them.
//!
//! Every box of every query asks here, from every thread (millions a
//! second on 32), and most asks find what is kept. So the entries are
//! split by key among `SHARDS` locks, and an ask that finds its entry takes
//! only its shard's; one that keeps a new tape takes the store's (the
//! tapes, the clock, the bytes) and then a shard's, always in that order.
//! One lock for all had its waiting threads sleeping and waking: on 32
//! threads, weighing a goBILDA kit spent 15 s of its 8 s wall in the
//! system, 4 s with the entries split.

use super::Tape;
use std::collections::{HashMap, VecDeque};
use std::hash::{Hash, Hasher};
use std::sync::{Arc, Mutex, MutexGuard, Weak};

/// The bytes the kept tapes may hold, about.
pub const BYTES: usize = 16 << 20;

/// The locks the entries are split among, by key.
const SHARDS: usize = 64;

/// What a tape was asked: simplified by choices, or an instruction's cone.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum Ask {
    Narrowed,
    Cone,
}

struct Entry {
    from: Weak<Tape>,
    ask: Ask,
    args: Box<[u8]>,
    /// The tape's key among the tapes (`Tape::content`).
    content: u64,
    tape: Arc<Tape>,
    reached: bool,
    id: u64,
}

impl Entry {
    fn bytes(&self) -> usize {
        std::mem::size_of::<Entry>() + std::mem::size_of::<Tape>() + self.args.len()
    }
}

/// A tape kept, and how many entries hold it.
struct Held {
    tape: Arc<Tape>,
    entries: usize,
}

/// One shard's entries, by key hash.
type Entries = HashMap<u64, Vec<Entry>>;

fn find<'a>(entries: &'a mut Entries, h: u64, from: &Arc<Tape>, ask: Ask, args: &[u8]) -> Option<&'a mut Entry> {
    entries.get_mut(&h)?.iter_mut().find(|e| e.from.as_ptr() == Arc::as_ptr(from) && e.ask == ask && *e.args == *args)
}

/// What every shard shares: the tapes, the clock and the bytes.
#[derive(Default)]
struct Store {
    /// By content hash: the tapes with it.
    tapes: HashMap<u64, Vec<Held>>,
    /// Entries in the order kept, by key hash and id: the clock's hand
    /// takes the first.
    order: VecDeque<(u64, u64)>,
    bytes: usize,
    next: u64,
}

impl Store {
    /// The tape kept alike to `tape` (content hash `c`), else `tape`, now
    /// held by one more entry.
    fn hold(&mut self, c: u64, tape: Arc<Tape>) -> Arc<Tape> {
        let list = self.tapes.entry(c).or_default();
        if let Some(h) = list.iter_mut().find(|h| h.tape.same(&tape)) {
            h.entries += 1;
            return h.tape.clone();
        }
        self.bytes += tape.bytes();
        list.push(Held { tape: tape.clone(), entries: 1 });
        tape
    }

    /// One entry fewer holds `tape` (content hash `c`): dropped with the
    /// last, into `gone`.
    fn release(&mut self, c: u64, tape: &Arc<Tape>, gone: &mut Vec<Arc<Tape>>) {
        let Some(list) = self.tapes.get_mut(&c) else { return };
        let Some(at) = list.iter().position(|h| Arc::ptr_eq(&h.tape, tape)) else { return };
        list[at].entries -= 1;
        if list[at].entries == 0 {
            self.bytes -= list[at].tape.bytes();
            gone.push(list.swap_remove(at).tape);
            if list.is_empty() {
                self.tapes.remove(&c);
            }
        }
    }

    /// Drops entries, oldest first, those reached since the hand passed
    /// them a turn more, until at most `cap` bytes are held. What goes is
    /// handed back, to be freed outside every lock.
    fn trim(&mut self, cap: usize, kept: &Kept) -> Vec<Arc<Tape>> {
        let mut gone = Vec::new();
        while self.bytes > cap {
            let Some((h, id)) = self.order.pop_front() else { break };
            let mut entries = kept.shard(h);
            let Some(list) = entries.get_mut(&h) else { continue };
            let Some(at) = list.iter().position(|e| e.id == id) else { continue };
            let e = &mut list[at];
            if e.reached && e.from.strong_count() > 0 {
                e.reached = false;
                self.order.push_back((h, id));
                continue;
            }
            let e = list.swap_remove(at);
            if list.is_empty() {
                entries.remove(&h);
            }
            drop(entries);
            self.bytes -= e.bytes();
            self.release(e.content, &e.tape, &mut gone);
            gone.push(e.tape);
        }
        gone
    }
}

/// Kept tapes, at most `cap` bytes.
pub(super) struct Kept {
    shards: Box<[Mutex<Entries>]>,
    store: Mutex<Store>,
    cap: usize,
}

impl Kept {
    pub(super) fn new(cap: usize) -> Kept {
        Kept { shards: (0..SHARDS).map(|_| Mutex::default()).collect(), store: Mutex::default(), cap }
    }

    /// The process's.
    pub(super) fn global() -> &'static Kept {
        static KEPT: std::sync::OnceLock<Kept> = std::sync::OnceLock::new();
        KEPT.get_or_init(|| Kept::new(BYTES))
    }

    /// The entries of key hash `h`'s shard, locked (by its top bits: the
    /// hash's low ones pick its bucket within).
    fn shard(&self, h: u64) -> MutexGuard<'_, Entries> {
        self.shards[(h >> (64 - SHARDS.trailing_zeros())) as usize].lock().unwrap()
    }

    /// The tape `ask`ed of `from` with `args`: kept, or `build`'s, kept
    /// where it fits.
    pub(super) fn get(&self, from: &Arc<Tape>, ask: Ask, args: &[u8], build: impl FnOnce() -> Tape) -> Arc<Tape> {
        let mut hasher = rustc_hash::FxHasher::default();
        (Arc::as_ptr(from) as usize, ask, args).hash(&mut hasher);
        let h = hasher.finish();
        if let Some(e) = find(&mut self.shard(h), h, from, ask, args) {
            e.reached = true;
            return e.tape.clone();
        }
        let tape = Arc::new(build());
        // One an eighth of the cap would drop most of the rest for itself
        // (a body of 170,000 instructions, whose every region is new).
        if tape.bytes() + args.len() > self.cap / 8 {
            return tape;
        }
        let c = tape.content();
        let mut s = self.store.lock().unwrap();
        let mut entries = self.shard(h);
        // Asked meanwhile on another thread: its tape, so both share it.
        if let Some(e) = find(&mut entries, h, from, ask, args) {
            return e.tape.clone();
        }
        let tape = s.hold(c, tape);
        let id = s.next;
        s.next += 1;
        let e = Entry { from: Arc::downgrade(from), ask, args: args.into(), content: c, tape: tape.clone(), reached: false, id };
        s.bytes += e.bytes();
        entries.entry(h).or_default().push(e);
        drop(entries);
        s.order.push_back((h, id));
        let gone = s.trim(self.cap, self);
        drop(s);
        drop(gone);
        tape
    }

    /// The bytes held.
    #[cfg(test)]
    fn bytes(&self) -> usize {
        self.store.lock().unwrap().bytes
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{BinOp, Graph, UnOp, Var};

    /// A tape of `n` terms: the max of x + k over k < n.
    fn tape(n: usize) -> Arc<Tape> {
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let mut f = x;
        for k in 1..n {
            let c = g.constant(k as f64);
            let t = g.bin(BinOp::Add, x, c);
            f = g.bin(BinOp::Max, f, t);
        }
        let f = g.un(UnOp::Neg, f);
        Arc::new(Tape::new(&g, f))
    }

    /// However much is asked, the cap holds; what is reached again stays
    /// while what is not goes; a tape derived alike from two tapes is one.
    #[test]
    fn kept_tapes_stay_within_the_cap_and_the_reached_stay() {
        let root = tape(64);
        let top = root.len() - 2;
        let one = root.cone(top).bytes() + 512;
        let kept = Kept::new(8 * one);
        let ask = |i: usize| kept.get(&root, Ask::Cone, &(i as u32).to_le_bytes(), || root.cone(i));
        let hot = ask(top);
        for i in 0..top {
            assert!(Arc::ptr_eq(&ask(top), &hot), "the one reached again is kept");
            ask(i);
            assert!(kept.bytes() <= 8 * one, "{} bytes held, at most {}", kept.bytes(), 8 * one);
        }
        let twin = Arc::new((*root).clone());
        let other = kept.get(&twin, Ask::Cone, &(top as u32).to_le_bytes(), || twin.cone(top));
        assert!(Arc::ptr_eq(&other, &hot), "alike from a tape alike: the one kept");
    }

    /// A tape too large to keep beside the rest (an eighth of the cap) is
    /// built and handed over, not kept: it would drop every other.
    #[test]
    fn a_tape_too_large_is_not_kept() {
        let (small, large) = (tape(16), tape(4096));
        let one = small.cone(small.len() - 2).bytes() + 512;
        let kept = Kept::new(64 * one);
        let hot: Vec<_> = (0..8).map(|i| kept.get(&small, Ask::Cone, &[i, 0, 0, 0], || small.cone(i as usize))).collect();
        let held = kept.bytes();
        let top = large.len() - 2;
        assert!(large.cone(top).bytes() > 64 * one, "the large one is larger than the cap");
        kept.get(&large, Ask::Cone, &(top as u32).to_le_bytes(), || large.cone(top));
        assert_eq!(kept.bytes(), held, "nothing kept or dropped");
        for (i, t) in hot.iter().enumerate() {
            assert!(Arc::ptr_eq(&kept.get(&small, Ask::Cone, &[i as u8, 0, 0, 0], || small.cone(i)), t));
        }
    }

    /// Asked from many threads at once, over keys spread across every
    /// shard: each answer is the tape that would be built; where all fits,
    /// every thread is handed the one tape kept for a key; where not, the
    /// cap still holds over all the shards together.
    #[test]
    fn asked_from_many_threads_each_gets_the_tape_and_the_cap_holds() {
        let root = tape(64);
        let keys = root.len() - 1;
        let built: Vec<Tape> = (0..keys).map(|i| root.cone(i)).collect();
        let one = built.iter().map(|t| t.bytes()).max().unwrap() + 512;
        for (cap, all_fit) in [(4 * keys * one, true), (6 * one, false)] {
            let kept = Kept::new(cap);
            let got: Vec<Vec<Arc<Tape>>> = std::thread::scope(|s| {
                let workers: Vec<_> = (0..8)
                    .map(|w| {
                        let (kept, root) = (&kept, &root);
                        s.spawn(move || {
                            (0..20 * keys).map(|n| (n * 7 + w * 13) % keys).map(|i| (i, kept.get(root, Ask::Cone, &(i as u32).to_le_bytes(), || root.cone(i)))).collect::<Vec<_>>()
                        })
                    })
                    .collect();
                let mut by_key = vec![Vec::new(); keys];
                for w in workers {
                    for (i, t) in w.join().unwrap() {
                        assert!(t.same(&built[i]), "key {i}: the tape that would be built");
                        by_key[i].push(t);
                    }
                }
                by_key
            });
            assert!(kept.bytes() <= cap, "{} bytes held, at most {cap}", kept.bytes());
            if all_fit {
                // Two threads building a key's tape side by side are both
                // handed the one kept first.
                for (i, ts) in got.iter().enumerate() {
                    assert!(ts.iter().all(|t| Arc::ptr_eq(t, &ts[0])), "key {i}: one tape kept, shared");
                }
            }
        }
    }

    /// At a cap of 0 nothing is held, and each ask builds its own.
    #[test]
    fn nothing_is_held_at_a_cap_of_none() {
        let kept = Kept::new(0);
        let root = tape(16);
        let a = kept.get(&root, Ask::Cone, &[3, 0, 0, 0], || root.cone(3));
        assert_eq!(kept.bytes(), 0);
        let b = kept.get(&root, Ask::Cone, &[3, 0, 0, 0], || root.cone(3));
        assert!(!Arc::ptr_eq(&a, &b) && a.same(&b));
    }
}
