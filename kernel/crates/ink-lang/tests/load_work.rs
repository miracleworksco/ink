//! The work a load does, counted as the bytes it allocates on the
//! loading thread: a long program's grows with what it computes, not with copies
//! of the values it reads.

use ink_lang::{Doc, Registry};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::path::Path;

struct Counting;

thread_local! {
    static BYTES: Cell<u64> = const { Cell::new(0) };
}

// SAFETY: forwards to the system allocator; counting touches only a
// thread-local cell, which needs no allocation.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        BYTES.with(|n| n.set(n.get() + l.size() as u64));
        unsafe { System.alloc(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        unsafe { System.dealloc(p, l) }
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, size: usize) -> *mut u8 {
        BYTES.with(|n| n.set(n.get() + size as u64));
        unsafe { System.realloc(p, l, size) }
    }
}

#[global_allocator]
static A: Counting = Counting;

static CORE: Registry = Registry::new(&[]);

/// The bytes allocated loading `src` after std, less those of std alone.
fn work(src: &str) -> u64 {
    let load = |s: &str| {
        let before = BYTES.with(Cell::get);
        Doc::load_str_in(&CORE, &format!("(import \"std\")\n{s}"), "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
        BYTES.with(Cell::get) - before
    };
    load(src) - load("")
}

/// A 30 mm strip bent +-90 degrees every 25 mm, n bends, folded.
fn meander(n: usize) -> String {
    let bends: String = (0..n)
        .map(|i| {
            let (u, a) = (20 + 25 * i, [90, -90, -90, 90][i % 4]);
            format!("(sheet-bend (vec {u} 30) (vec {u} 0) {a} 2) ")
        })
        .collect();
    let l = 25 * n + 10;
    format!("(def S (sheet (polygon (vec (vec 0 0) (vec {l} 0) (vec {l} 30) (vec 0 30))) (vec) 2 0.4 (vec {bends})))\n(body p (shape (sheet-folded S)))")
}

/// A sheet's bends are checked in pairs (a line across another's zone),
/// so its load grows as the square of their number: three times the
/// bends, at most 10.5 times the work. Copying each list read (a bend's
/// numbers, the table of which bends lie beyond which) made it grow as the
/// cube and more; so did ordering each bend's lineage by a search per
/// depth.
#[test]
fn a_long_chain_of_bends_loads_in_work_square_in_its_bends() {
    let (w30, w90) = (work(&meander(30)), work(&meander(90)));
    let ratio = w90 as f64 / w30 as f64;
    assert!(ratio < 10.5, "30 bends: {w30} bytes, 90: {w90} ({ratio:.2} times)");
}

/// A list read many times is shared, not copied: reading each element of
/// a 20,000-long list by index does work in its length, not its square.
#[test]
fn reading_a_list_does_not_copy_it() {
    let read = |n: usize| work(&format!("(def big (range {n}))\n(def s (fold + 0 (map (lambda (i) (nth big i)) (range {n}))))"));
    let (w1, w2) = (read(10_000), read(20_000));
    let ratio = w2 as f64 / w1 as f64;
    assert!(ratio < 2.5, "10,000: {w1} bytes, 20,000: {w2} ({ratio:.2} times)");
}
