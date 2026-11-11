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

use super::Tape;
use std::collections::{HashMap, VecDeque};
use std::hash::{Hash, Hasher};
use std::sync::{Arc, Mutex, Weak};

/// The bytes the kept tapes may hold, about.
pub const BYTES: usize = 16 << 20;

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

#[derive(Default)]
struct Inner {
    /// By key hash: the entries with it.
    entries: HashMap<u64, Vec<Entry>>,
    /// By content hash: the tapes with it.
    tapes: HashMap<u64, Vec<Held>>,
    /// Entries in the order kept, by key hash and id: the clock's hand
    /// takes the first.
    order: VecDeque<(u64, u64)>,
    bytes: usize,
    next: u64,
}

impl Inner {
    fn find(&mut self, h: u64, from: &Arc<Tape>, ask: Ask, args: &[u8]) -> Option<&mut Entry> {
        self.entries.get_mut(&h)?.iter_mut().find(|e| e.from.as_ptr() == Arc::as_ptr(from) && e.ask == ask && *e.args == *args)
    }

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
    /// last.
    fn release(&mut self, c: u64, tape: &Arc<Tape>) {
        let Some(list) = self.tapes.get_mut(&c) else { return };
        let Some(at) = list.iter().position(|h| Arc::ptr_eq(&h.tape, tape)) else { return };
        list[at].entries -= 1;
        if list[at].entries == 0 {
            self.bytes -= list[at].tape.bytes();
            list.swap_remove(at);
            if list.is_empty() {
                self.tapes.remove(&c);
            }
        }
    }

    /// Drops entries, oldest first, those reached since the hand passed
    /// them a turn more, until at most `cap` bytes are held.
    fn trim(&mut self, cap: usize) {
        while self.bytes > cap {
            let Some((h, id)) = self.order.pop_front() else { break };
            let Some(list) = self.entries.get_mut(&h) else { continue };
            let Some(at) = list.iter().position(|e| e.id == id) else { continue };
            let e = &mut list[at];
            if e.reached && e.from.strong_count() > 0 {
                e.reached = false;
                self.order.push_back((h, id));
                continue;
            }
            let e = list.swap_remove(at);
            if list.is_empty() {
                self.entries.remove(&h);
            }
            self.bytes -= e.bytes();
            self.release(e.content, &e.tape);
        }
    }
}

/// Kept tapes, at most `cap` bytes.
pub(super) struct Kept {
    inner: Mutex<Inner>,
    cap: usize,
}

impl Kept {
    pub(super) fn new(cap: usize) -> Kept {
        Kept { inner: Mutex::default(), cap }
    }

    /// The process's.
    pub(super) fn global() -> &'static Kept {
        static KEPT: std::sync::OnceLock<Kept> = std::sync::OnceLock::new();
        KEPT.get_or_init(|| Kept::new(BYTES))
    }

    /// The tape `ask`ed of `from` with `args`: kept, or `build`'s, kept
    /// where it fits.
    pub(super) fn get(&self, from: &Arc<Tape>, ask: Ask, args: &[u8], build: impl FnOnce() -> Tape) -> Arc<Tape> {
        let mut hasher = rustc_hash::FxHasher::default();
        (Arc::as_ptr(from) as usize, ask, args).hash(&mut hasher);
        let h = hasher.finish();
        if let Some(e) = self.inner.lock().unwrap().find(h, from, ask, args) {
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
        let mut s = self.inner.lock().unwrap();
        // Asked meanwhile on another thread: its tape, so both share it.
        if let Some(e) = s.find(h, from, ask, args) {
            return e.tape.clone();
        }
        let tape = s.hold(c, tape);
        let id = s.next;
        s.next += 1;
        let e = Entry { from: Arc::downgrade(from), ask, args: args.into(), content: c, tape: tape.clone(), reached: false, id };
        s.bytes += e.bytes();
        s.entries.entry(h).or_default().push(e);
        s.order.push_back((h, id));
        s.trim(self.cap);
        tape
    }

    /// The bytes held.
    #[cfg(test)]
    fn bytes(&self) -> usize {
        self.inner.lock().unwrap().bytes
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
