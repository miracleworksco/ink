//! A requirement answered at each value of a variable a tool reads its
//! facts at (`over`): a made-up motor whose limit falls with speed, and a
//! motion of three moments, the torque highest at the first and the limit
//! lowest at the last.

use ink_claims::*;
use std::cell::Cell;

/// The motion: the torque (N) and the speed at each moment.
struct Moment<'a>(&'a [(f64, f64)], usize, &'a Cell<usize>);

impl Tool for Moment<'_> {
    fn ask(&self, asked: &About, _: &mut Claims) -> Option<Fact> {
        let (torque, speed) = self.0[self.1];
        let value = [("torque", torque), ("speed", speed)].into_iter().find(|(what, _)| asked.on == "m" && asked.what == *what)?.1;
        self.2.set(self.2.get() + 1);
        Some(Fact {
            covers: vec![("t".into(), Covers::One)],
            strength: Some(Strength::Estimated),
            notes: vec![format!("at moment {}", self.1)],
            ..Fact::stated(asked.clone(), value, "N", "the motion")
        })
    }
}

/// The motor's limit at its speed: 10 N less a tenth of it.
struct Curve;

impl Tool for Curve {
    fn ask(&self, asked: &About, claims: &mut Claims) -> Option<Fact> {
        (*asked == about("limit", "m")).then(|| {
            let speed = claims.answer(&about("speed", "m"));
            let mut fact = resting(asked, "the curve", [&speed].into_iter());
            (fact.value, fact.unit) = (speed.value.map(|s| 10.0 - s / 10.0), "N".into());
            fact
        })
    }
}

fn design() -> Design {
    let rule = |name: &str, says| Rule { about: about(name, "m"), says, from: name.into() };
    Design {
        facts: vec![Fact::stated(about("weight", "m"), 3.0, "N", "stated")],
        rules: vec![
            rule("under-the-curve", Says::Ordered(vec![Term::Size(about("torque", "m")), Term::Read(about("limit", "m"))])),
            rule("light", Says::Ordered(vec![Term::Read(about("weight", "m")), Term::Number(5.0, "N".into())])),
            Rule { about: about("ok", ""), says: Says::Within("m".into()), from: "ok".into() },
        ],
    }
}

/// `asked` over `motion`, and how many facts its moments gave.
fn over_motion(asked: &About, motion: &[(f64, f64)]) -> (Fact, usize, usize) {
    let (design, gave) = (design(), Cell::new(0));
    let at = |moment: usize| (moment < motion.len()).then(|| Box::new(Moment(motion, moment, &gave)) as Box<dyn Tool + '_>);
    let (fact, moment) = over(&design, asked, "t", &[&Curve], 0, &at);
    (fact, moment, gave.get())
}

#[test]
fn a_requirement_is_answered_at_every_moment_and_its_answer_is_the_worst() {
    // The torque is highest first (8 of the 10 the curve allows at rest); the least room is last (5 of 6).
    let motion = [(8.0, 0.0), (-6.0, 20.0), (5.0, 40.0)];
    let (fact, moment, _) = over_motion(&about("under-the-curve", "m"), &motion);
    assert_eq!((fact.held, fact.strength, moment), (Some(Verdict::Yes), Some(Strength::Estimated), 2));
    assert_eq!((fact.margin, &fact.covers), (Some(1.0), &vec![("t".to_string(), Covers::These(3))]));
    assert_eq!(fact.terms, [("torque on m".to_string(), 5.0), ("limit on m".to_string(), 6.0)]);
    // A little faster at the end, the curve is under the torque there: no, though the peak is where it was.
    let (fact, moment, _) = over_motion(&about("under-the-curve", "m"), &[(8.0, 0.0), (-6.0, 20.0), (5.0, 60.0)]);
    assert_eq!((fact.held, moment, fact.margin), (Some(Verdict::No), 2, Some(-1.0)));
    assert_eq!(fact.covers, [("t".to_string(), Covers::These(3))]);
    // The requirement over every rating: no where one rating is no at one moment, yes where each is yes at each.
    let (all, _, _) = over_motion(&about("ok", ""), &[(8.0, 0.0), (-6.0, 20.0), (5.0, 60.0)]);
    assert_eq!((all.held, &all.covers), (Some(Verdict::No), &vec![("t".to_string(), Covers::These(3))]));
    assert_eq!(over_motion(&about("ok", ""), &motion).0.held, Some(Verdict::Yes));
    // One moment is one value of t.
    assert_eq!(over_motion(&about("under-the-curve", "m"), &motion[..1]).0.covers, [("t".to_string(), Covers::One)]);
}

#[test]
fn what_reads_no_moment_is_answered_once_and_says_nothing_of_time() {
    let motion = [(8.0, 0.0), (-6.0, 20.0), (5.0, 40.0)];
    let (light, moment, gave) = over_motion(&about("light", "m"), &motion);
    assert_eq!((light.held, light.strength, light.covers.len(), moment, gave), (Some(Verdict::Yes), None, 0, 0, 0));
    // With no moments at all the tools answer as they would alone.
    let (none, _, _) = over_motion(&about("under-the-curve", "m"), &[]);
    assert_eq!(none.held, Some(Verdict::NotAvailable(Unavailable::NoSource)));
}

#[test]
fn a_comparison_says_how_much_room_it_had() {
    let design = design();
    let light = Claims::new(&design, &[]).answer(&about("light", "m"));
    assert_eq!((light.held, light.margin), (Some(Verdict::Yes), Some(2.0)));
    // Of bands, from the high end of one to the low end of the next: 3 known between 2.5 and 3.5 has 1.5 of room.
    let banded = Design { facts: vec![Fact { band: Some((2.5, 3.5)), ..design.facts[0].clone() }], ..design };
    assert_eq!(Claims::new(&banded, &[]).answer(&about("light", "m")).margin, Some(1.5));
}
