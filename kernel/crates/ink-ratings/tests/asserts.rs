//! The prover's answer to a claim as a fact of the core (`ink_ratings::fact`): its verdict, how strong, and
//! over how much of what it reads.

use ink_claims::*;
use ink_lang::claims::{ClaimOptions, run_asserts};
use ink_lang::{Registry, command};

static CORE: Registry = Registry::new(&[]);

/// The claims of the design `src`, each answered by the prover and made a fact.
fn facts(name: &str, src: &str) -> Vec<Fact> {
    let file = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    std::fs::write(&file, src).unwrap();
    let words = ["require", file.to_str().unwrap()].map(String::from);
    let mut view = ink_ratings::load(&CORE, &command::Args::parse(&words)).unwrap_or_else(|e| panic!("{}", e.message));
    let params: Vec<String> = view.doc.params.iter().map(|p| p.name.clone()).collect();
    let report = run_asserts(&mut view, None, true, ClaimOptions::default()).unwrap_or_else(|e| panic!("{}", e.message));
    report.asserts.iter().map(|answer| ink_ratings::fact(&params, answer)).collect()
}

/// Four claims: one that reads t, one that reads a ranged size, one that reads neither, one that fails.
const CLAIMS: &str = "(def w (ranged 4 6))\n(require early (value t < 1s))\n(require wide (value w > 3))\n(require plain (value 2 > 1))\n(require narrow (value w > 5))\n";

#[test]
fn a_claim_that_reads_time_says_over_how_much_of_it_it_was_shown() {
    // No span: asked at one moment, as its row's `one_moment` records. Proven there.
    let [early, wide, plain, narrow] = &facts("moment.ink", CLAIMS)[..] else { panic!("four claims") };
    assert_eq!((early.held, early.strength), (Some(Verdict::Yes), Some(Strength::Proven)), "{early:?}");
    assert_eq!(early.covers, [("t".to_string(), Covers::One)]);
    assert!(early.from.ends_with("moment.ink:2") && early.notes[0] == "(value t < 1s)" && early.notes.iter().any(|n| n.starts_with("asked at t = 0 only")), "{early:?}");
    // A ranged size is covered at every value, by a proof or a disproof; a claim that reads no variable covers none.
    assert_eq!((wide.held, wide.strength, &wide.covers), (Some(Verdict::Yes), Some(Strength::Proven), &vec![("w".to_string(), Covers::Every)]), "{wide:?}");
    assert_eq!((plain.held, plain.covers.len()), (Some(Verdict::Yes), 0), "{plain:?}");
    assert_eq!(narrow.held, Some(Verdict::No));
    // Under a span the same claim is shown over all of it.
    let early = &facts("span.ink", &format!("(t-span 0 0.5s)\n{CLAIMS}"))[0];
    assert_eq!((early.held, &early.covers), (Some(Verdict::Yes), &vec![("t".to_string(), Covers::Every)]), "{early:?}");
    // And past where it holds, it fails.
    assert_eq!(facts("late.ink", &format!("(t-span 0 2s)\n{CLAIMS}"))[0].held, Some(Verdict::No));
}

#[test]
fn an_unknown_says_what_kind() {
    // A number outside ink's range is outside what the prover covers: not available, with the prover's reason.
    let huge = &facts("unknown.ink", "(require huge (value (* 1e200 1e200) > 0))\n")[0];
    assert_eq!((huge.held, huge.strength, huge.covers.len()), (Some(Verdict::NotAvailable(Unavailable::Outside)), None, 0), "{huge:?}");
    assert!(huge.notes.len() == 2 && !huge.notes[1].is_empty(), "{huge:?}");
}

/// A joint's `(limits ...)` makes a claim nobody writes, a fact like any other.
#[test]
fn a_limit_clauses_claim_is_a_fact_without_being_written() {
    let src = "(def box (max (- (abs x) 1) (- (abs y) 1) (- (abs z) 1)))\n(body arm (shape box))\n\
               (assembly root (instance a arm (at) (turn z (* 30deg/s t)) (limits -10 10)))\n(t-span 0 1s)\n(require fits (value 2 > 1))\n";
    let found = facts("unwritten.ink", src);
    // 30 degrees at t = 1 is past the 10 allowed: no, over the file's span.
    assert_eq!(found.iter().map(|f| (f.about.what.as_str(), f.held)).collect::<Vec<_>>(), [("limits:root/a", Some(Verdict::No)), ("fits", Some(Verdict::Yes))]);
}
