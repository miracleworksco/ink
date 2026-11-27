//! The claims core against real answers, a real rating, and fake tools.

use ink_claims::*;
use serde_json::Value;
use std::cell::Cell;

const NO_TOOLS: &[&dyn Tool] = &[];

// ---------------------------------------------------------------------
// 1. Today's answers load as facts and give back today's verdicts.
//
// data/answers.jsonl is 63 of the 1,741 rows `ink require --own` printed for
// the 141 designs with an `; expect-status:` line (2026-11-19, b0bfe6a306),
// each as printed, plus its file and, for an unknown, `why`. One row per
// kind, the first met: a kind is the status; proven or estimated; `held`;
// which generated check it is (limits, speed-limits, effort-limits, relate)
// or a plain claim; whether it prints `terms`; how it covers t; whether
// another param is ranged; and, for an unknown, the class of its sentence.
// `why` is a reading of that sentence: today's rows state no typed reason.
// ---------------------------------------------------------------------

fn rows() -> Vec<Value> {
    include_str!("data/answers.jsonl").lines().map(|l| serde_json::from_str(l).unwrap()).collect()
}

/// One of today's answer rows as a fact.
fn fact(row: &Value) -> Fact {
    let text = |k: &str| row[k].as_str().map(String::from);
    let name = text("name").unwrap();
    let (what, on) = name.split_once(':').unwrap_or((&name, ""));
    let decided = row["status"] != "unknown";
    let motion = &row["motion"];
    let mut covers = vec![];
    match motion["how"].as_str() {
        Some("t ranged") => covers.push(("t".to_string(), Covers::Every)),
        Some("measured") => covers.push(("t".to_string(), Covers::These(motion["steps"].as_u64().unwrap() as usize))),
        Some(how) => panic!("{how}"),
        None if !row["one_moment"].is_null() => covers.push(("t".to_string(), Covers::One)),
        None => {}
    }
    for (param, _) in row["param_spans"].as_object().into_iter().flatten() {
        // Today only a proven row says how much of a ranged param it covered.
        if param != "t" && row["accuracy"].is_number() && decided {
            covers.push((param.clone(), Covers::Every));
        }
    }
    let pair = |v: &Value| Some((v[0].as_f64()?, v[1].as_f64()?));
    Fact {
        about: about(what, on),
        value: row["value"].as_f64(),
        unit: String::new(), // today's rows print bare numbers
        held: Some(match (row["status"].as_str().unwrap(), row["why"].as_str()) {
            ("holds", _) => Verdict::Yes,
            ("fails", _) => Verdict::No,
            ("unknown", Some("lacks")) => Verdict::NotSettled(Unsettled::Lacks),
            ("unknown", Some("undecided")) => Verdict::NotSettled(Unsettled::Undecided),
            ("unknown", Some("straddles")) => Verdict::NotSettled(Unsettled::Straddles),
            ("unknown", Some("outside")) => Verdict::NotAvailable(Unavailable::Outside),
            other => panic!("{other:?}"),
        }),
        // The accuracy a statement was proven to is a number inside the statement.
        terms: (row["terms"].as_object().into_iter().flatten().filter_map(|(k, v)| Some((k.clone(), v.as_f64()?))))
            .chain(row["accuracy"].as_f64().map(|a| ("accuracy".to_string(), a)))
            .collect(),
        band: pair(&row["bounds_over_scope"]).or(row["band"]["low"].as_f64().zip(row["band"]["high"].as_f64())),
        margin: None, // no row says how much room its comparison had
        strength: match row["accuracy"].as_f64() {
            _ if !decided => None,
            Some(_) => Some(if row["held"] == "to the accuracy" { Strength::ToAccuracy } else { Strength::Proven }),
            None => Some(Strength::Estimated),
        },
        covers,
        left_out: (row["excludes"].as_array().into_iter().flatten().map(|e| e["what"].as_str().unwrap().to_string()))
            .chain(motion["undecided"].as_array().filter(|u| !u.is_empty()).map(|u| format!("t in {}", Value::from(u.clone()))))
            .collect(),
        from: text("source").unwrap(),
        used: row["solved_inputs"].as_array().into_iter().flatten().map(|s| about(s.as_str().unwrap(), "")).collect(),
        notes: ["unavailable", "undecided", "note", "hint"]
            .iter()
            .filter_map(|k| text(k))
            .chain(row["notes"].as_array().into_iter().flatten().map(|n| n.as_str().unwrap().to_string()))
            .collect(),
    }
}

#[test]
fn todays_answers_load_as_facts_with_todays_verdicts() {
    let (mut yes, mut no, mut open) = (0, 0, 0);
    for row in rows() {
        let stated = fact(&row);
        let design = Design { facts: vec![stated.clone()], rules: vec![] };
        let answer = Claims::new(&design, NO_TOOLS).answer(&stated.about);
        assert_eq!(answer, stated, "a stated answer comes back whole");
        let status = row["status"].as_str().unwrap();
        match answer.held.unwrap() {
            Verdict::Yes => (assert_eq!(status, "holds"), yes += 1),
            Verdict::No => (assert_eq!(status, "fails"), no += 1),
            _ => (assert_eq!(status, "unknown"), open += 1),
        };
        // Today's two strengths: an accuracy, or none on a decided row.
        assert_eq!(matches!(answer.strength, Some(Strength::Proven | Strength::ToAccuracy)), status != "unknown" && row["accuracy"].is_number(), "{row}");
        assert_eq!(answer.strength == Some(Strength::Estimated), status != "unknown" && row["accuracy"].is_null(), "{row}");
    }
    assert_eq!((yes, no, open), (29, 16, 18));
}

/// The rows that print the numbers inside their statement (`terms`: a
/// reading, `lo` and `hi`) are worked out again from those numbers alone.
#[test]
fn todays_limit_rows_are_worked_out_again_from_their_numbers() {
    let mut worked = 0;
    for row in rows() {
        let today = fact(&row);
        let Some(terms) = row["terms"].as_object() else { continue };
        // The reading is the term that is not a limit; an effort row also prints the rate it was read at.
        let mut reading: Vec<&str> = terms.keys().map(|k| k.as_str()).filter(|k| !["lo", "hi"].contains(k)).collect();
        if reading.len() > 1 {
            reading.retain(|k| *k != "rate");
        }
        assert_eq!(reading.len(), 1, "{row}");
        let on = today.about.on.as_str();
        let design = Design {
            facts: terms.iter().map(|(k, v)| Fact::stated(about(k, on), v.as_f64().unwrap(), "", &today.from)).collect(),
            rules: vec![Rule { about: today.about.clone(), says: Says::Ordered(["lo", reading[0], "hi"].map(|k| Term::Read(about(k, on))).to_vec()), from: today.from.clone() }],
        };
        let answer = Claims::new(&design, NO_TOOLS).answer(&today.about);
        assert_eq!(answer.held, today.held, "{row}");
        assert_eq!(answer.terms.len(), 3);
        assert_eq!(answer.strength, None, "it rests on stated numbers: nothing was proven");
        worked += 1;
    }
    assert_eq!(worked, 10);
}

// ---------------------------------------------------------------------
// 2. The real stepper. data/17HS19-1684S-PG27.json is the `gearbox` entry
// of the catalog's record of that gearmotor in its nema parts, copied whole.
// ---------------------------------------------------------------------

const SHAFT: &str = "opener/motor/shaft";

/// The gearmotor's two shaft ratings, placed as `opener/motor`.
fn ratings() -> Vec<Rule> {
    let record: Value = serde_json::from_str(include_str!("data/17HS19-1684S-PG27.json")).unwrap();
    let rated = |name: &str, load: &str, key: &str| Rule {
        about: about(name, "opener/motor"),
        says: Says::Ordered(vec![Term::Read(about(load, SHAFT)), Term::Number(record["gearbox"][key].as_f64().unwrap(), "N".into())]),
        from: format!("17HS19-1684S-PG27 gearbox {key}"),
    };
    vec![rated("side-load", "radial-load", "shaft-radial-load-N-max"), rated("thrust", "axial-load", "shaft-axial-load-N-max")]
}

fn motor_ok() -> Rule {
    Rule { about: about("motor-ok", ""), says: Says::Within("opener".into()), from: "garage.ink".into() }
}

/// The garage door design: a 50 lbf door hangs from a sprocket on the
/// motor's shaft, which the design puts at 223 N across it.
const DOOR_N: f64 = 223.0;

#[test]
fn the_gearmotor_is_not_within_its_ratings_under_the_garage_door() {
    let load = Fact::stated(about("radial-load", SHAFT), DOOR_N, "N", "garage.ink: the door's weight on the sprocket");
    let design = Design { facts: vec![load], rules: [ratings(), vec![motor_ok()]].concat() };
    let mut claims = Claims::new(&design, NO_TOOLS);

    let answer = claims.answer(&about("motor-ok", ""));
    assert_eq!(answer.held, Some(Verdict::No));
    assert_eq!(answer.used, [about("side-load", "opener/motor"), about("thrust", "opener/motor")]);

    // The story: which rating, the numbers, and where each came from.
    let side = claims.answer(&about("side-load", "opener/motor"));
    assert_eq!(side.held, Some(Verdict::No));
    assert_eq!(side.terms, [("radial-load on opener/motor/shaft".to_string(), 223.0)]);
    assert_eq!(design.rules[0].says, Says::Ordered(vec![Term::Read(about("radial-load", SHAFT)), Term::Number(100.0, "N".into())]));
    assert_eq!(side.from, "17HS19-1684S-PG27 gearbox shaft-radial-load-N-max");
    assert_eq!(side.strength, None, "as stated");
    assert_eq!(claims.answer(&side.used[0]).from, "garage.ink: the door's weight on the sprocket");

    // Nothing states a load along the shaft, so the other rating has no answer; the no stands.
    assert_eq!(claims.answer(&about("thrust", "opener/motor")).held, Some(Verdict::NotAvailable(Unavailable::NoSource)));

    // A part elsewhere is not asked about, and "every rating holds" where nothing is rated is not a yes.
    assert_eq!(Claims::new(&design, NO_TOOLS).answer(&about("motor-ok", "")).used.len(), 2);
    let elsewhere = Rule { about: about("ok", ""), says: Says::Within("open".into()), from: String::new() };
    let design = Design { rules: [design.rules.clone(), vec![elsewhere]].concat(), ..design.clone() };
    let answer = Claims::new(&design, NO_TOOLS).answer(&about("ok", ""));
    assert_eq!((answer.held, answer.used.len()), (Some(Verdict::NotAvailable(Unavailable::Empty)), 0));

    // A load stated in another unit is refused, not compared as if it were newtons.
    let pounds = Design { facts: vec![Fact::stated(about("radial-load", SHAFT), 50.0, "lbf", "garage.ink")], ..design.clone() };
    assert_eq!(Claims::new(&pounds, NO_TOOLS).answer(&about("side-load", "opener/motor")).held, Some(Verdict::NotAvailable(Unavailable::Units)));
}

// ---------------------------------------------------------------------
// 3. Tools, plugged in from here. Their numbers are made up, except the
// door's.
// ---------------------------------------------------------------------

/// Offers one fact, and counts how often it is asked.
struct Fake {
    fact: Fact,
    asked: Cell<usize>,
}

impl Fake {
    fn new(what: &str, value: f64, band: Option<(f64, f64)>) -> Fake {
        let fact = Fact {
            band,
            strength: Some(Strength::Estimated),
            covers: vec![("t".into(), Covers::These(64))],
            ..Fact::stated(about(what, SHAFT), value, "N", "fake motion: rigid parts")
        };
        Fake { fact, asked: Cell::new(0) }
    }
}

impl Tool for Fake {
    fn ask(&self, about: &About, _: &mut Claims) -> Option<Fact> {
        (*about == self.fact.about).then(|| {
            self.asked.set(self.asked.get() + 1);
            self.fact.clone()
        })
    }
}

/// Works the load across the shaft out from a fact it asks for.
struct Hangs;

impl Tool for Hangs {
    fn ask(&self, asked: &About, claims: &mut Claims) -> Option<Fact> {
        (*asked == about("radial-load", SHAFT)).then(|| {
            let weight = claims.answer(&about("weight", "opener/door"));
            Fact {
                about: asked.clone(),
                used: vec![weight.about.clone()],
                from: "fake statics: the whole weight across the shaft".into(),
                strength: Some(Strength::Estimated),
                ..weight
            }
        })
    }
}

#[test]
fn a_tool_answers_what_the_facts_cannot_in_the_order_given_and_stated_facts_that_disagree_unsettle_it() {
    let design = Design { facts: vec![], rules: [ratings(), vec![motor_ok()]].concat() };
    let (side, motor) = (about("side-load", "opener/motor"), about("motor-ok", ""));

    // The facts alone cannot say.
    assert_eq!(Claims::new(&design, NO_TOOLS).answer(&side).held, Some(Verdict::NotAvailable(Unavailable::NoSource)));

    // Nothing runs unless asked; then each thing is asked once.
    let motion = Fake::new("radial-load", DOOR_N, None);
    let tools: &[&dyn Tool] = &[&motion];
    let mut claims = Claims::new(&design, tools);
    assert_eq!(motion.asked.get(), 0);
    let answer = claims.answer(&side);
    assert_eq!(answer.held, Some(Verdict::No));
    assert_eq!(answer.strength, Some(Strength::Estimated), "as strong as the weakest fact it rests on");
    assert_eq!(answer.covers, [("t".to_string(), Covers::These(64))]);
    assert_eq!(claims.answer(&motor).held, Some(Verdict::No));
    assert_eq!(motion.asked.get(), 1);

    // An answer covers each variable only as far as the least of the facts it rests on.
    let limit = Fact { covers: vec![("t".into(), Covers::Every)], ..Fact::stated(about("limit", SHAFT), 300.0, "N", "made up") };
    let under = Rule { about: about("under", ""), says: Says::Ordered(vec![Term::Read(about("radial-load", SHAFT)), Term::Read(limit.about.clone())]), from: String::new() };
    let both = Design { facts: vec![limit], rules: vec![under] };
    let answer = Claims::new(&both, tools).answer(&about("under", ""));
    assert_eq!((answer.held, &answer.covers[..]), (Some(Verdict::Yes), &[("t".to_string(), Covers::These(64))][..]));

    // Tools are tried in the order given. The first that settles it answers, and no later one runs.
    let (first, other) = (Fake::new("radial-load", DOOR_N, None), Fake::new("radial-load", 80.0, None));
    let tools: &[&dyn Tool] = &[&first, &other];
    assert_eq!(Claims::new(&design, tools).answer(&side).held, Some(Verdict::No));
    assert_eq!((first.asked.get(), other.asked.get()), (1, 0));
    let tools: &[&dyn Tool] = &[&other, &first];
    assert_eq!(Claims::new(&design, tools).answer(&side).held, Some(Verdict::Yes));
    assert_eq!((first.asked.get(), other.asked.get()), (1, 1));

    // A later tool runs only when the earlier ones left it open; if all do, the first reply stands.
    let mut lacks = Fake::new("radial-load", 0.0, None);
    lacks.fact.held = Some(Verdict::NotSettled(Unsettled::Lacks));
    let mut gave_up = Fake::new("radial-load", 0.0, None);
    gave_up.fact.held = Some(Verdict::NotSettled(Unsettled::Undecided));
    let tools: &[&dyn Tool] = &[&lacks, &other];
    assert_eq!(Claims::new(&design, tools).answer(&side).held, Some(Verdict::Yes));
    assert_eq!(other.asked.get(), 2);
    let tools: &[&dyn Tool] = &[&lacks, &gave_up];
    assert_eq!(Claims::new(&design, tools).answer(&side).held, Some(Verdict::NotSettled(Unsettled::Lacks)));

    // Two stated facts on either side of the limit: sources disagree, and both are shown.
    let rig = |n: f64, from: &str| Fact::stated(about("radial-load", SHAFT), n, "N", from);
    let two = Design { facts: vec![rig(223.0, "garage.ink"), rig(80.0, "measured on the rig")], ..design.clone() };
    let mut claims = Claims::new(&two, tools);
    let answer = claims.answer(&side);
    assert_eq!(answer.held, Some(Verdict::NotSettled(Unsettled::Disagree)));
    assert_eq!(answer.held.unwrap().to_string(), "not settled: sources disagree");
    assert_eq!(claims.sources(&answer.used[0]).iter().map(|f| f.value.unwrap()).collect::<Vec<_>>(), [223.0, 80.0]);
    assert_eq!(claims.answer(&motor).held, Some(Verdict::NotSettled(Unsettled::Disagree)));
    // On the same side of it they still differ, so they still disagree; in value, unit or band.
    let two = Design { facts: vec![rig(223.0, "garage.ink"), rig(250.0, "measured on the rig")], ..design.clone() };
    assert_eq!(Claims::new(&two, tools).answer(&side).held, Some(Verdict::NotSettled(Unsettled::Disagree)));
    for other in [Fact { unit: "lbf".into(), ..rig(223.0, "") }, Fact { band: Some((200.0, 240.0)), ..rig(223.0, "") }] {
        let two = Design { facts: vec![rig(223.0, "garage.ink"), other], ..design.clone() };
        assert_eq!(Claims::new(&two, tools).answer(&side).held, Some(Verdict::NotSettled(Unsettled::Disagree)));
    }
    // The same number stated twice is one answer.
    let two = Design { facts: vec![rig(223.0, "garage.ink"), rig(223.0, "measured on the rig")], ..design.clone() };
    assert_eq!(Claims::new(&two, tools).answer(&side).held, Some(Verdict::No));

    // A band across the limit is not settled; a stated fact stands and no tool is asked.
    let wide = Fake::new("radial-load", 95.0, Some((90.0, 110.0)));
    let tools: &[&dyn Tool] = &[&wide];
    assert_eq!(Claims::new(&design, tools).answer(&side).held, Some(Verdict::NotSettled(Unsettled::Straddles)));
    let stated = Design { facts: vec![Fact::stated(about("radial-load", SHAFT), 60.0, "N", "measured on the rig")], ..design.clone() };
    assert_eq!(Claims::new(&stated, tools).answer(&side).held, Some(Verdict::Yes));
    assert_eq!(wide.asked.get(), 1);

    // A tool asks for a fact while it works, and says what it used.
    let door = Design { facts: vec![Fact::stated(about("weight", "opener/door"), DOOR_N, "N", "garage.ink")], ..design.clone() };
    let tools: &[&dyn Tool] = &[&Hangs];
    let mut claims = Claims::new(&door, tools);
    assert_eq!(claims.answer(&side).held, Some(Verdict::No));
    assert_eq!(claims.answer(&about("radial-load", SHAFT)).used, [about("weight", "opener/door")]);
    // Without the weight the tool's reply is open, and its reason passes up.
    assert_eq!(Claims::new(&design, tools).answer(&side).held, Some(Verdict::NotAvailable(Unavailable::NoSource)));
}

// ---------------------------------------------------------------------
// 4. What nothing can answer is not available, and says which way.
// ---------------------------------------------------------------------

#[test]
fn what_no_fact_or_tool_can_answer_is_not_available() {
    let reads =
        |name: &str, read: &str| Rule { about: about(name, ""), says: Says::Ordered(vec![Term::Read(about(read, "")), Term::Number(2.0, String::new())]), from: String::new() };
    let design = Design {
        facts: vec![Fact { about: about("fits", ""), held: Some(Verdict::Yes), ..Fact::default() }],
        rules: vec![reads("strong", "safety-factor"), reads("a", "b"), reads("b", "a"), reads("odd", "fits")],
    };
    let unrelated = Fake::new("radial-load", DOOR_N, None);
    let tools: &[&dyn Tool] = &[&unrelated];
    let mut claims = Claims::new(&design, tools);

    let answer = claims.answer(&about("strong", ""));
    assert_eq!(answer.held, Some(Verdict::NotAvailable(Unavailable::NoSource)));
    assert_eq!((answer.strength, &answer.used[..]), (None, &[about("safety-factor", "")][..]));
    assert_eq!(answer.held.unwrap().to_string(), "not available: nothing states it or offers it");
    assert_eq!(claims.answer(&about("never-written", "")).held, Some(Verdict::NotAvailable(Unavailable::NoSource)));
    assert_eq!(unrelated.asked.get(), 0);

    // A requirement that rests on itself, and one that reads a yes/no as a number.
    assert_eq!(claims.answer(&about("a", "")).held, Some(Verdict::NotAvailable(Unavailable::Circular)));
    assert_eq!(claims.answer(&about("odd", "")).held, Some(Verdict::NotAvailable(Unavailable::Outside)));
}

// ---------------------------------------------------------------------
// 5. The outside review's wrong answers, each as it was reproduced.
// ---------------------------------------------------------------------

fn limit(name: &str, read: &str, n: f64, unit: &str) -> Rule {
    Rule { about: about(name, "m"), says: Says::Ordered(vec![Term::Read(about(read, "s")), Term::Number(n, unit.into())]), from: "r".into() }
}

fn ok(design: &Design, tools: &[&dyn Tool]) -> Fact {
    Claims::new(design, tools).answer(&about("ok", "m"))
}

/// A number within its own band of the limit is not a yes. The band is the only accuracy a number has.
#[test]
fn accuracy_at_the_limit() {
    let x = |band| Fact { band, strength: Some(Strength::Proven), ..Fact::stated(about("x", "s"), 99.95, "N", "x") };
    let design = |band| Design { facts: vec![x(band)], rules: vec![limit("ok", "x", 100.0, "N")] };
    assert_eq!(ok(&design(Some((99.85, 100.05))), NO_TOOLS).held, Some(Verdict::NotSettled(Unsettled::Straddles)));
    assert_eq!(ok(&design(Some((99.94, 99.96))), NO_TOOLS).held, Some(Verdict::Yes));
    assert_eq!(ok(&design(Some((100.01, 100.05))), NO_TOOLS).held, Some(Verdict::No));
    assert_eq!(ok(&design(None), NO_TOOLS).held, Some(Verdict::Yes), "a number with no band is exact");
}

/// Asks the core for other facts and passes the first on under its own name.
struct Asks(About, Vec<About>);

impl Tool for Asks {
    fn ask(&self, asked: &About, claims: &mut Claims) -> Option<Fact> {
        (*asked == self.0).then(|| {
            let answers: Vec<Fact> = self.1.iter().map(|a| claims.answer(a)).collect();
            Fact { about: asked.clone(), ..answers[0].clone() }
        })
    }
}

/// A is offered by a tool that asks for B, B by a tool that asks for A, and then A by a tool that just states it.
/// B was first seen while A was still being found; that must not stay as "rests on itself".
#[test]
fn circular_taint() {
    let [a, b, c, d] = ["A", "B", "C", "D"].map(|x| about(x, ""));
    let (first, back) = (Asks(a.clone(), vec![b.clone()]), Asks(b.clone(), vec![a.clone()]));
    let gives = Fake { fact: Fact::stated(a.clone(), 5.0, "N", "stated by a tool"), asked: Cell::new(0) };
    let (further, last) = (Asks(c.clone(), vec![d.clone()]), Fake { fact: Fact::stated(d.clone(), 7.0, "N", ""), asked: Cell::new(0) });
    let tools: &[&dyn Tool] = &[&first, &back, &gives, &further, &last];
    let design = Design::default();
    let mut claims = Claims::new(&design, tools);
    assert_eq!((claims.answer(&a).held, claims.answer(&a).value), (None, Some(5.0)));
    assert_eq!((claims.answer(&b).held, claims.answer(&b).value), (None, Some(5.0)));
    assert_eq!(Claims::new(&design, tools).answer(&b).value, Some(5.0), "the same asked first");
    // Once the circle is past, what is found is kept again.
    assert_eq!((claims.answer(&c).value, claims.answer(&d).value, last.asked.get(), gives.asked.get()), (Some(7.0), Some(7.0), 1, 2));

    // B asks for A while A is being found, and then for itself: A's stand-in still counts.
    let back = Asks(b.clone(), vec![a.clone(), b.clone()]);
    let tools: &[&dyn Tool] = &[&first, &back, &gives];
    let mut claims = Claims::new(&design, tools);
    assert_eq!((claims.answer(&a).value, claims.answer(&b).value), (Some(5.0), Some(5.0)));

    // A real circle still says so, whichever end is asked first.
    let back = Asks(b.clone(), vec![a.clone()]);
    let tools: &[&dyn Tool] = &[&first, &back];
    for order in [[&a, &b], [&b, &a]] {
        let mut claims = Claims::new(&design, tools);
        assert_eq!(order.map(|x| claims.answer(x).held), [Some(Verdict::NotAvailable(Unavailable::Circular)); 2]);
    }
}

/// A plain number stated under a rating's own name is no answer to "does every rating hold".
#[test]
fn within_skips_unheld() {
    let all = Rule { about: about("all", ""), says: Says::Within("m".into()), from: String::new() };
    let design = Design {
        facts: vec![Fact::stated(about("load", "s"), 500.0, "N", "x"), Fact::stated(about("ok", "m"), 1.0, "", "typo")],
        rules: vec![limit("ok", "load", 100.0, "N"), all],
    };
    assert_eq!(Claims::new(&design, NO_TOOLS).answer(&about("all", "")).held, Some(Verdict::NotAvailable(Unavailable::Outside)));
    // Without the stray number the rating answers, and it is a no.
    let design = Design { facts: vec![design.facts[0].clone()], ..design };
    assert_eq!(Claims::new(&design, NO_TOOLS).answer(&about("all", "")).held, Some(Verdict::No));
}

/// A comparison of fewer than two terms checks nothing, and a number that is not one settles nothing.
#[test]
fn degenerate_ordered() {
    let held = |terms: Vec<Term>| {
        let design = Design { facts: vec![], rules: vec![Rule { about: about("r", ""), says: Says::Ordered(terms), from: String::new() }] };
        Claims::new(&design, NO_TOOLS).answer(&about("r", "")).held.unwrap()
    };
    let n = |n: f64| Term::Number(n, String::new());
    assert_eq!(held(vec![]), Verdict::NotAvailable(Unavailable::Empty));
    assert_eq!(held(vec![n(1.0)]), Verdict::NotAvailable(Unavailable::Empty));
    assert_eq!(held(vec![n(1.0), n(2.0)]), Verdict::Yes);
    assert_eq!(held(vec![n(f64::NAN), n(1.0)]), Verdict::NotSettled(Unsettled::Straddles));
    assert_eq!(held(vec![n(1.0), n(f64::NAN)]), Verdict::NotSettled(Unsettled::Straddles));
}

/// Sources that say the same give one answer, as strong as the weakest and covering as the least.
#[test]
fn agree_keeps_the_weakest_strength() {
    let x = |strength, covers| Fact { strength, covers: vec![("t".into(), covers)], ..Fact::stated(about("x", "s"), 10.0, "N", "x") };
    let rules = vec![limit("ok", "x", 100.0, "N")];
    for facts in
        [vec![x(None, Covers::Every), x(Some(Strength::Estimated), Covers::One)], vec![x(Some(Strength::Estimated), Covers::One), x(Some(Strength::Proven), Covers::Every)]]
    {
        let answer = ok(&Design { facts, rules: rules.clone() }, NO_TOOLS);
        assert_eq!((answer.held, answer.strength, &answer.covers[..]), (Some(Verdict::Yes), Some(Strength::Estimated), &[("t".to_string(), Covers::One)][..]));
    }
    // Two rules under one name are each found once and rested on once each.
    let all = Rule { about: about("all", ""), says: Says::Within("m".into()), from: String::new() };
    let design = Design { facts: vec![x(None, Covers::Every)], rules: vec![rules[0].clone(), limit("ok", "x", 200.0, "N"), all] };
    assert_eq!(Claims::new(&design, NO_TOOLS).answer(&about("all", "")).held, Some(Verdict::Yes));
}

// ---------------------------------------------------------------------
// 6. Too close to call, the strength of a verdict, sizes and rounding.
// ---------------------------------------------------------------------

fn gives(what: &str, value: f64, band: Option<(f64, f64)>, strength: Strength) -> Fake {
    Fake { fact: Fact { band, strength: Some(strength), ..Fact::stated(about(what, "s"), value, "N", "a tool") }, asked: Cell::new(0) }
}

/// Offers one fact as another passed on, and says it used it.
struct Through(&'static str, &'static str);

impl Tool for Through {
    fn ask(&self, asked: &About, claims: &mut Claims) -> Option<Fact> {
        (*asked == about(self.0, "s")).then(|| Fact { about: asked.clone(), used: vec![about(self.1, "s")], ..claims.answer(&about(self.1, "s")) })
    }
}

/// A cheap answer too close to the limit to call goes on to the next tool that offers the fact, and no further.
#[test]
fn too_close_to_call_goes_on_to_the_next_tool() {
    let design = Design { facts: vec![], rules: vec![limit("ok", "x", 100.0, "N")] };
    let (rough, fine, never) = (gives("x", 95.0, Some((90.0, 110.0)), Strength::Estimated), gives("x", 97.0, None, Strength::Estimated), gives("x", 1.0, None, Strength::Proven));
    let answer = ok(&design, &[&rough, &fine, &never]);
    assert_eq!((answer.held, answer.terms[0].1, answer.strength), (Some(Verdict::Yes), 97.0, Some(Strength::Estimated)));
    assert_eq!((rough.asked.get(), fine.asked.get(), never.asked.get()), (1, 1, 0));
    // Alone, or with another that cannot call it either, it stays not settled; each tool was asked once.
    assert_eq!(ok(&design, &[&rough]).held, Some(Verdict::NotSettled(Unsettled::Straddles)));
    let rough_too = gives("x", 99.0, Some((98.0, 103.0)), Strength::Estimated);
    assert_eq!(ok(&design, &[&rough, &rough_too]).held, Some(Verdict::NotSettled(Unsettled::Straddles)));
    assert_eq!((rough.asked.get(), rough_too.asked.get()), (3, 1));
    // A cheap answer that settles it is the answer: the next tool is not asked.
    let clear = gives("x", 50.0, Some((40.0, 60.0)), Strength::Estimated);
    assert_eq!((ok(&design, &[&clear, &fine]).held, fine.asked.get()), (Some(Verdict::Yes), 1));
    // A stated fact too close to call stands: no tool is asked to second-guess it.
    let stated = Design { facts: vec![Fact { band: Some((90.0, 110.0)), ..Fact::stated(about("x", "s"), 95.0, "N", "stated") }], ..design.clone() };
    assert_eq!((ok(&stated, &[&fine]).held, fine.asked.get()), (Some(Verdict::NotSettled(Unsettled::Straddles)), 1));

    // The fact too close to call may be one the requirement reads through another: y is x passed on by a tool.
    let through = Through("y", "x");
    let design = Design { facts: vec![], rules: vec![limit("ok", "y", 100.0, "N")] };
    let (rough, fine) = (gives("x", 95.0, Some((90.0, 110.0)), Strength::Estimated), gives("x", 97.0, None, Strength::Estimated));
    let tools: &[&dyn Tool] = &[&through, &rough, &fine];
    let mut claims = Claims::new(&design, tools);
    assert_eq!((claims.answer(&about("ok", "m")).held, rough.asked.get(), fine.asked.get()), (Some(Verdict::Yes), 1, 1));
    // What rested on the first answer, however far up, is forgotten with it: z is y passed on, asked before.
    let further = Through("z", "y");
    let (rough, fine) = (gives("x", 95.0, Some((90.0, 110.0)), Strength::Estimated), gives("x", 97.0, None, Strength::Estimated));
    let tools: &[&dyn Tool] = &[&further, &through, &rough, &fine];
    let mut claims = Claims::new(&design, tools);
    assert_eq!(claims.answer(&about("z", "s")).value, Some(95.0));
    assert_eq!((claims.answer(&about("ok", "m")).held, claims.answer(&about("z", "s")).value), (Some(Verdict::Yes), Some(97.0)));
    let (rough, fine) = (gives("x", 95.0, Some((90.0, 110.0)), Strength::Estimated), gives("x", 97.0, None, Strength::Estimated));
    let tools: &[&dyn Tool] = &[&through, &rough, &fine];
    let mut claims = Claims::new(&design, tools);
    assert_eq!((claims.answer(&about("ok", "m")).held, rough.asked.get(), fine.asked.get()), (Some(Verdict::Yes), 1, 1));
    // The better answer is then the fact's answer for everything else in the round.
    assert_eq!((claims.answer(&about("x", "s")).value, claims.answer(&about("y", "s")).value, fine.asked.get()), (Some(97.0), Some(97.0), 1));
}

/// A verdict is as strong as the facts that decide it; one that rests only on stated numbers is not proven.
#[test]
fn a_verdict_is_as_strong_as_what_decides_it() {
    let within = Rule { about: about("all", ""), says: Says::Within("m".into()), from: String::new() };
    let rules = vec![limit("ok", "x", 100.0, "N"), limit("also", "y", 100.0, "N"), within];
    let answers = |x: &Fake, y: &Fake| {
        let design = Design { facts: vec![], rules: rules.clone() };
        let tools: &[&dyn Tool] = &[x, y];
        let mut claims = Claims::new(&design, tools);
        ["ok", "also"]
            .map(|r| claims.answer(&about(r, "m")).strength)
            .into_iter()
            .chain([claims.answer(&about("all", "")).held.map(|_| Strength::Proven), claims.answer(&about("all", "")).strength])
            .collect::<Vec<_>>()
    };
    let (proven, estimated) = (Some(Strength::Proven), Some(Strength::Estimated));
    // A no is decided by the rating that fails, however the ones that pass are known.
    assert_eq!(answers(&gives("x", 500.0, None, Strength::Proven), &gives("y", 5.0, None, Strength::Estimated)), [proven, estimated, proven, proven]);
    assert_eq!(answers(&gives("x", 500.0, None, Strength::Estimated), &gives("y", 5.0, None, Strength::Proven)), [estimated, proven, proven, estimated]);
    // Of two that fail, the stronger decides it. A yes needs every one, so it is as strong as the weakest.
    assert_eq!(answers(&gives("x", 500.0, None, Strength::Estimated), &gives("y", 500.0, None, Strength::Proven))[3], proven);
    assert_eq!(answers(&gives("x", 5.0, None, Strength::Proven), &gives("y", 5.0, None, Strength::Estimated))[3], estimated);
    assert_eq!(answers(&gives("x", 5.0, None, Strength::Proven), &gives("y", 5.0, None, Strength::ToAccuracy))[3], Some(Strength::ToAccuracy));

    // Stated numbers prove nothing: with a proven fact the verdict is as stated, with an estimate an estimate.
    let stated = |y: Option<Strength>| {
        let design =
            Design { facts: vec![Fact::stated(about("x", "s"), 5.0, "N", "stated"), Fact { strength: y, ..Fact::stated(about("y", "s"), 5.0, "N", "") }], rules: rules.clone() };
        let mut claims = Claims::new(&design, NO_TOOLS);
        (claims.answer(&about("ok", "m")).strength, claims.answer(&about("all", "")).strength)
    };
    assert_eq!([stated(None), stated(proven), stated(estimated)], [(None, None), (None, None), (None, estimated)]);
    // A comparison of numbers alone is arithmetic.
    let design = Design {
        facts: vec![],
        rules: vec![Rule { about: about("r", ""), says: Says::Ordered(vec![Term::Number(1.0, "N".into()), Term::Number(2.0, "N".into())]), from: String::new() }],
    };
    assert_eq!(Claims::new(&design, NO_TOOLS).answer(&about("r", "")).strength, proven);
}

/// A rating on a size reads the number without its sign; one on the number itself keeps it.
#[test]
fn a_size_is_compared_without_its_sign() {
    let held = |term: fn(About) -> Term, value: f64, band| {
        let rule = Rule { about: about("ok", "m"), says: Says::Ordered(vec![term(about("x", "s")), Term::Number(100.0, "N".into())]), from: String::new() };
        ok(&Design { facts: vec![Fact { band, ..Fact::stated(about("x", "s"), value, "N", "") }], rules: vec![rule] }, NO_TOOLS).held.unwrap()
    };
    assert_eq!([held(Term::Size, -500.0, None), held(Term::Size, 500.0, None), held(Term::Size, -50.0, None)], [Verdict::No, Verdict::No, Verdict::Yes]);
    assert_eq!(held(Term::Read, -500.0, None), Verdict::Yes, "signed: -500 is under 100");
    // A band about nothing has sizes from nothing to its farther end; one clear of it, between its ends' sizes.
    assert_eq!(
        [held(Term::Size, 1.0, Some((-50.0, 30.0))), held(Term::Size, 1.0, Some((-150.0, 30.0))), held(Term::Size, -120.0, Some((-150.0, -110.0)))],
        [Verdict::Yes, Verdict::NotSettled(Unsettled::Straddles), Verdict::No]
    );
}

/// Numbers that differ only in the last place are equal: 1.001 kN said as newtons is 1000.9999999999999.
#[test]
fn rounding_is_not_a_difference() {
    let held = |x: f64, limit: f64| {
        let design = Design { facts: vec![Fact::stated(about("x", "s"), x, "N", "")], rules: vec![crate::limit("ok", "x", limit, "N")] };
        ok(&design, NO_TOOLS).held.unwrap()
    };
    assert_eq!(
        [held(1001.0, 1.001 * 1000.0), held(1001.0, 1001.0), held(1001.000001, 1001.0), held(0.1, 0.09999999999999998)],
        [Verdict::Yes, Verdict::Yes, Verdict::No, Verdict::Yes]
    );
    assert!(same(0.1 + 0.2, 0.3) && !same(1.0, 1.000000001) && same(0.0, 0.0) && !same(0.0, 1e-300));
    // Four floating-point steps, at any size, and no more: 100.00000000000009 is six above 100.
    let steps = |x: f64, n: u64| f64::from_bits(x.to_bits() + n);
    for x in [1e-9, 0.1, 100.0, 1e12] {
        assert!(same(x, steps(x, 4)) && !same(x, steps(x, 5)) && same(-x, -steps(x, 4)) && !same(-x, -steps(x, 5)) && !same(x, -x), "{x}");
    }
    assert_eq!(100.00000000000009_f64, steps(100.0, 6));
    assert_eq!([held(100.00000000000009, 100.0), held(steps(100.0, 4), 100.0), held(100.0000001, 100.0)], [Verdict::No, Verdict::Yes, Verdict::No]);
}

// ---------------------------------------------------------------------
// 7. The second outside read's three, on a fact moving to a later tool.
// ---------------------------------------------------------------------

/// Offers x and y, each between 90 and 110.
struct Passes;

impl Tool for Passes {
    fn ask(&self, asked: &About, _: &mut Claims) -> Option<Fact> {
        matches!(asked.what.as_str(), "x" | "y").then(|| Fact { band: Some((90.0, 110.0)), ..Fact::stated(asked.clone(), 100.0, "N", "tool") })
    }
}

fn part(name: &str) -> About {
    about(name, "part")
}

fn rating(name: &str, terms: Vec<Term>) -> Rule {
    Rule { about: part(name), says: Says::Ordered(terms), from: "rating".into() }
}

fn newtons(n: f64) -> Term {
    Term::Number(n, "N".into())
}

/// A later rating moves x to a tool that says -1: the earlier rating, already answered yes on the first tool's x, is
/// answered again, and the whole is no.
#[test]
fn a_rating_answered_before_a_fact_moved_is_answered_again() {
    let better = Fake { fact: Fact::stated(part("x"), -1.0, "N", "tool"), asked: Cell::new(0) };
    let all = Rule { about: about("all", ""), says: Says::Within("part".into()), from: "require".into() };
    let design = Design { facts: vec![], rules: vec![rating("min", vec![newtons(0.0), Term::Read(part("x"))]), rating("max", vec![Term::Read(part("x")), newtons(100.0)]), all] };
    let tools: &[&dyn Tool] = &[&Passes, &better];
    let mut claims = Claims::new(&design, tools);
    assert_eq!(
        (claims.answer(&about("all", "")).held, claims.answer(&part("min")).held, claims.answer(&part("max")).held),
        (Some(Verdict::No), Some(Verdict::No), Some(Verdict::Yes))
    );
    assert_eq!(better.asked.get(), 1);
    // Asked the other way round, the same.
    let mut claims = Claims::new(&design, tools);
    assert_eq!((claims.answer(&part("min")).held, claims.answer(&part("max")).held, claims.answer(&part("min")).held), (Some(Verdict::Yes), Some(Verdict::Yes), Some(Verdict::No)));
    assert_eq!(claims.answer(&about("all", "")).held, Some(Verdict::No));
}

/// Counts what it is asked for x, and offers nothing.
struct Declines(Cell<usize>);

impl Tool for Declines {
    fn ask(&self, asked: &About, _: &mut Claims) -> Option<Fact> {
        if *asked == part("x") {
            self.0.set(self.0.get() + 1);
        }
        None
    }
}

#[test]
fn a_tool_that_declines_is_asked_once() {
    let design = Design { facts: vec![], rules: vec![rating("first", vec![Term::Read(part("x")), newtons(100.0)]), rating("second", vec![Term::Read(part("x")), newtons(101.0)])] };
    let declines = Declines(Cell::new(0));
    let tools: &[&dyn Tool] = &[&Passes, &declines];
    let mut claims = Claims::new(&design, tools);
    for rating in ["first", "second", "first"] {
        assert_eq!(claims.answer(&part(rating)).held, Some(Verdict::NotSettled(Unsettled::Straddles)));
    }
    assert_eq!(declines.0.get(), 1);
}

/// Offers x as y and y as x, and counts.
struct ReadsOther(Cell<usize>);

impl Tool for ReadsOther {
    fn ask(&self, asked: &About, claims: &mut Claims) -> Option<Fact> {
        let other = match asked.what.as_str() {
            "x" => part("y"),
            "y" => part("x"),
            _ => return None,
        };
        self.0.set(self.0.get() + 1);
        assert!(self.0.get() < 50, "the round does not end");
        Some(Fact { about: asked.clone(), used: vec![other.clone()], ..claims.answer(&other) })
    }
}

/// Two facts each worked out from the other by a later tool: the round ends, within its bound.
#[test]
fn facts_that_read_each_other_end() {
    let design = Design { facts: vec![], rules: vec![rating("check", vec![Term::Read(part("x")), Term::Read(part("y"))])] };
    let takes_turns = ReadsOther(Cell::new(0));
    let tools: &[&dyn Tool] = &[&Passes, &takes_turns];
    let mut claims = Claims::new(&design, tools);
    assert_eq!(claims.answer(&part("check")).held, Some(Verdict::NotSettled(Unsettled::Straddles)));
    // x moves to it (asked once); y is not moved to it, since its y is worked out from x, which rests on y (asked once).
    assert_eq!(takes_turns.0.get(), 2);
    assert_eq!(
        (claims.answer(&part("check")).held, claims.answer(&part("x")).band, takes_turns.0.get()),
        (Some(Verdict::NotSettled(Unsettled::Straddles)), Some((90.0, 110.0)), 2)
    );
}

/// A tool that says a number and a yes, each resting on a fact it did not look at.
struct Careless(Vec<About>);

impl Tool for Careless {
    fn ask(&self, asked: &About, _: &mut Claims) -> Option<Fact> {
        let fact = Fact { used: self.0.clone(), ..Fact::stated(asked.clone(), 5.0, "N", "careless") };
        match asked.what.as_str() {
            "worked" => Some(fact),
            "fine" => Some(Fact { value: None, unit: String::new(), held: Some(Verdict::Yes), ..fact }),
            "broken" => Some(Fact { value: None, unit: String::new(), held: Some(Verdict::No), ..fact }),
            _ => None,
        }
    }
}

/// A yes or a number that rests on an open fact is open, for that fact's reason, whoever said it: no line says yes
/// beside one that is not settled or not available. A no stands; and so does what rests on a name nothing
/// states or offers (an old answer's inputs), which is no fact here.
#[test]
fn a_yes_that_rests_on_an_open_fact_is_open() {
    let load = about("load", SHAFT);
    let two = Design { facts: vec![Fact::stated(load.clone(), 1.0, "N", "one"), Fact::stated(load.clone(), 2.0, "N", "two")], rules: vec![] };
    let careless = Careless(vec![load.clone()]);
    let tools: &[&dyn Tool] = &[&careless];
    let mut claims = Claims::new(&two, tools);
    let disagree = Some(Verdict::NotSettled(Unsettled::Disagree));
    assert_eq!(claims.answer(&load).held, disagree);
    for what in ["worked", "fine"] {
        let fact = claims.answer(&about(what, SHAFT));
        assert_eq!((fact.held, fact.value, fact.unit.as_str(), fact.strength), (disagree, None, "", None), "{what}");
        assert_eq!(fact.used, std::slice::from_ref(&load), "it still says what it rests on");
    }
    assert_eq!(claims.answer(&about("broken", SHAFT)).held, Some(Verdict::No));
    // On a fact that is settled, and on a name the core knows nothing of, it is what the tool said.
    let one = Design { facts: vec![Fact::stated(load.clone(), 1.0, "N", "one")], rules: vec![] };
    assert_eq!(Claims::new(&one, tools).answer(&about("fine", SHAFT)).held, Some(Verdict::Yes));
    let cites = Careless(vec![about("an input of the old check", "")]);
    let tools: &[&dyn Tool] = &[&cites];
    assert_eq!(Claims::new(&one, tools).answer(&about("worked", SHAFT)).value, Some(5.0));
}

/// A tool that measures the design's own model, giving one fact, worse the way it says.
struct Model(Fact, Worse);

impl Tool for Model {
    fn ask(&self, about: &About, _: &mut Claims) -> Option<Fact> {
        (*about == self.0.about).then(|| self.0.clone())
    }
    fn measures(&self) -> Option<Worse> {
        Some(self.1)
    }
}

/// A typed number never makes an answer better than what the model works out. Of a demand (a load: more is
/// worse): typed at or past the model's most it stands, an allowance; typed short of the model's least the two
/// disagree, each shown; typed inside the model's band it does not narrow it: the fact is the model's, and says
/// what was typed. Of what meets a demand (a life: less is worse) the same, the other way round.
#[test]
fn a_stated_fact_is_held_beside_what_the_model_works_out() {
    let load = about("radial-load", SHAFT);
    let stated = |n: f64| Design { facts: vec![Fact::stated(load.clone(), n, "N", "hoist.ink:30")], rules: vec![] };
    let model = |n: f64, band, unit: &str, worse| {
        Model(Fact { band, unit: unit.into(), strength: Some(Strength::Estimated), ..Fact::stated(load.clone(), n, "N", "ink effort on root") }, worse)
    };
    let answer = |typed: f64, worked: &Model| Claims::new(&stated(typed), &[worked]).answer(&load);
    // A bearing rated 1370 N; the model works out 1380 N within 2%: 1352.4 to 1407.6.
    let worked = model(1380.0, Some((1352.4, 1407.6)), "N", Worse::More);
    // At or past its most: the typed one stands, as stated.
    for typed in [1407.6, 1500.0, -1500.0] {
        let a = answer(typed, &worked);
        assert_eq!((a.value, a.band, a.held, a.strength, a.from.as_str()), (Some(typed), None, None, None, "hoist.ink:30"));
    }
    // Inside its band, under the model's own number or over it: the model's fact, its band never narrowed
    // toward the better (its most stays), its least no better than what was typed; and what was typed said.
    for typed in [1360.0, 1352.4, 1390.0, 1407.0] {
        let a = answer(typed, &worked);
        assert_eq!((a.value, a.band, a.strength, a.from.as_str()), (Some(typed.max(1380.0)), Some((typed, 1407.6)), Some(Strength::Estimated), "ink effort on root"), "{typed}");
        assert_eq!(a.notes, [format!("stated as {typed} N (hoist.ink:30): inside what the model works out, which a typed number does not narrow")]);
    }
    // Short of its least: they disagree, each shown.
    let a = answer(1340.0, &worked);
    assert_eq!((a.value, a.held), (None, Some(Verdict::NotSettled(Unsettled::Disagree))));
    assert_eq!(a.notes[..2], ["hoist.ink:30: 1340 N", "ink effort on root: 1380 N, between 1352.4 and 1407.6 N"]);
    assert!(a.notes[2].starts_with("the model works out more than is stated: state no less than it, or leave the fact out"), "{:?}", a.notes);
    // A model's number with no band: the same number stands, one under it disagrees.
    let exact = model(-60.0, None, "N", Worse::More);
    assert_eq!((answer(60.0, &exact).value, answer(59.0, &exact).held), (Some(60.0), Some(Verdict::NotSettled(Unsettled::Disagree))));
    // A typed band that starts inside the model's and reaches past it: the model's, reaching as far.
    let reaching = Design { facts: vec![Fact { band: Some((1400.0, 1450.0)), ..Fact::stated(load.clone(), 1420.0, "N", "hoist.ink:30") }], rules: vec![] };
    assert_eq!(Claims::new(&reaching, &[&worked]).answer(&load).band, Some((1400.0, 1450.0)));
    // The model's sign is kept: -1380 N, typed as 1360.
    let turned = model(-1380.0, Some((-1407.6, -1352.4)), "N", Worse::More);
    assert_eq!(answer(1360.0, &turned).band, Some((-1407.6, -1360.0)));
    // In another unit nothing is compared: what reads the fact says so.
    assert_eq!(answer(60.0, &model(120.0, None, "lbf", Worse::More)).value, Some(60.0));

    // Where less is worse (a life the model works out as 100 h, 90 to 110): typed under its least stands;
    // typed over its most disagrees; typed inside is the model's.
    let life = model(100.0, Some((90.0, 110.0)), "N", Worse::Less);
    assert_eq!((answer(80.0, &life).value, answer(80.0, &life).band), (Some(80.0), None));
    let over = answer(120.0, &life);
    assert_eq!(over.held, Some(Verdict::NotSettled(Unsettled::Disagree)));
    assert!(over.notes[2].starts_with("the model works out less than is stated: state no more than it"), "{:?}", over.notes);
    assert_eq!((answer(105.0, &life).value, answer(105.0, &life).band), (Some(100.0), Some((90.0, 105.0))));

    // A tool that measures nothing is not asked beside a stated fact.
    let other = Fake::new("radial-load", 120.0, None);
    assert_eq!((Claims::new(&stated(60.0), &[&other]).answer(&load).value, other.asked.get()), (Some(60.0), 0));
    // Stated twice: both, each with its line.
    let twice = Design { facts: vec![Fact::stated(load.clone(), 60.0, "N", "hoist.ink:30"), Fact::stated(load.clone(), 90.0, "N", "hoist.ink:31")], rules: vec![] };
    assert_eq!(Claims::new(&twice, &[]).answer(&load).notes, ["hoist.ink:30: 60 N", "hoist.ink:31: 90 N"]);
}
