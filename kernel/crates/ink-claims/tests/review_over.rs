//! `over` against an outside read's cases: each was once a wrong answer.

use ink_claims::*;

struct Sample(Option<f64>);

impl Tool for Sample {
    fn ask(&self, a: &About, _: &mut Claims) -> Option<Fact> {
        if a.what != "load" {
            return None;
        }
        Some(Fact { covers: vec![("t".into(), Covers::One)], ..Fact::stated(a.clone(), self.0?, "N", "at the sample") })
    }
}

struct Fixed(&'static str, f64);

impl Tool for Fixed {
    fn ask(&self, a: &About, _: &mut Claims) -> Option<Fact> {
        (a.what == self.0).then(|| Fact { covers: vec![("t".into(), Covers::One)], ..Fact::stated(a.clone(), self.1, "N", "read only at t = 0") })
    }
}

fn check(samples: &[Option<f64>], fixed: &Fixed, terms: Vec<Term>) -> (Fact, usize) {
    let asked = about("ok", "");
    let design = Design { facts: vec![], rules: vec![Rule { about: asked.clone(), says: Says::Ordered(terms), from: "review input".into() }] };
    let at = |i: usize| samples.get(i).map(|v| Box::new(Sample(*v)) as Box<dyn Tool>);
    over(&design, &asked, "t", &[fixed], 0, &at)
}

#[test]
fn a_tool_that_offers_nothing_first_can_still_fail_later() {
    let terms = vec![Term::Read(about("load", "")), Term::Number(5.0, "N".into())];
    let (fact, at) = check(&[None, Some(8.0)], &Fixed("load", 2.0), terms.clone());
    // Nothing is decided from the first moment: the second's 8 N is over the 5.
    assert_eq!((fact.held, at), (Some(Verdict::No), 1));
    let (late, _) = check(&[Some(8.0)], &Fixed("load", 2.0), terms);
    assert_eq!(late.held, Some(Verdict::No));
}

#[test]
fn a_fact_known_at_one_time_does_not_cover_three() {
    let terms = vec![Term::Read(about("load", "")), Term::Read(about("limit", ""))];
    let (fact, _) = check(&[Some(2.0), Some(3.0), Some(4.0)], &Fixed("limit", 5.0), terms);
    // Three loads against a limit known at one moment: shown at one value of t, not three.
    assert_eq!(fact.held, Some(Verdict::Yes));
    assert_eq!(fact.covers, vec![("t".into(), Covers::One)]);
}

#[test]
fn an_open_moment_does_not_hide_a_definite_failure() {
    let terms = vec![Term::Read(about("load", "")), Term::Number(5.0, "N".into())];
    let (fact, at) = check(&[Some(2.0), None, Some(8.0)], &Fixed("unrelated", 0.0), terms);
    assert_eq!((fact.held, at, fact.margin), (Some(Verdict::No), 2, Some(-3.0)));
}

#[test]
fn a_three_term_comparison_uses_both_adjacent_gaps() {
    let terms = vec![Term::Number(2.0, "N".into()), Term::Size(about("load", "")), Term::Number(7.0, "N".into())];
    let (fact, at) = check(&[Some(-4.0), Some(6.0), Some(2.5)], &Fixed("unrelated", 0.0), terms);
    assert_eq!((fact.held, at, fact.margin), (Some(Verdict::Yes), 2, Some(0.5)));
}

#[test]
fn a_gap_that_is_no_number_does_not_hide_the_worst_case() {
    let terms = vec![Term::Read(about("load", "")), Term::Read(about("load", "")), Term::Number(5.0, "N".into())];
    let (late, _) = check(&[Some(f64::INFINITY)], &Fixed("unrelated", 0.0), terms.clone());
    let (fact, at) = check(&[Some(6.0), Some(f64::INFINITY)], &Fixed("unrelated", 0.0), terms);
    // A load without end is the worse of the two, and its room is no "not a number" that sorts as plenty.
    assert_eq!((late.held, late.margin), (Some(Verdict::No), Some(f64::NEG_INFINITY)));
    assert_eq!((fact.held, at, fact.margin), (Some(Verdict::No), 1, Some(f64::NEG_INFINITY)));
}

#[test]
fn escalation_is_local_to_the_moment_and_each_tool_is_asked_once() {
    use std::cell::Cell;
    struct Banded {
        band: (f64, f64),
        calls: std::rc::Rc<Cell<usize>>,
    }
    impl Tool for Banded {
        fn ask(&self, a: &About, _: &mut Claims) -> Option<Fact> {
            if a.what != "load" {
                return None;
            }
            self.calls.set(self.calls.get() + 1);
            Some(Fact { band: Some(self.band), covers: vec![("t".into(), Covers::One)], ..Fact::stated(a.clone(), self.band.0, "N", "banded sample") })
        }
    }
    struct Fallback(Cell<usize>);
    impl Tool for Fallback {
        fn ask(&self, a: &About, _: &mut Claims) -> Option<Fact> {
            if a.what != "load" {
                return None;
            }
            self.0.set(self.0.get() + 1);
            Some(Fact::stated(a.clone(), 4.0, "N", "fallback"))
        }
    }
    let asked = about("ok", "");
    let design = Design {
        facts: vec![],
        rules: vec![Rule { about: asked.clone(), says: Says::Ordered(vec![Term::Read(about("load", "")), Term::Number(5.0, "N".into())]), from: "review input".into() }],
    };
    let samples = [(4.0, 6.0), (7.0, 8.0)];
    let calls = samples.map(|_| std::rc::Rc::new(Cell::new(0)));
    let at = |i: usize| samples.get(i).map(|band| Box::new(Banded { band: *band, calls: calls[i].clone() }) as Box<dyn Tool>);
    let fallback = Fallback(Cell::new(0));
    let (fact, i) = over(&design, &asked, "t", &[&fallback], 0, &at);
    assert_eq!((fact.held, i), (Some(Verdict::No), 1));
    assert_eq!((calls[0].get(), calls[1].get(), fallback.0.get()), (1, 1, 1));
}
