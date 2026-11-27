//! A made-up geared part carrying two ratings, placed twice in a design
//! that states facts on each and asks whether each is within its ratings;
//! then one part with one rating at a time (data/cases.json), for each way
//! a fact was once read wrongly.

use ink_claims::*;
use ink_lang::view::ModelView;
use ink_lang::{InkError, Registry, command};
use ink_ratings::{Derived, design};
use std::sync::atomic::{AtomicUsize, Ordering};

static CORE: Registry = Registry::new(&[]);
static FILES: AtomicUsize = AtomicUsize::new(0);

/// The design `src` loaded as a command would load it, with these `--param`s.
fn load_with(src: &str, params: &[&str]) -> Result<ModelView, InkError> {
    static DATA: std::sync::Once = std::sync::Once::new();
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR"));
    DATA.call_once(|| {
        for data in ["geared.json", "bad.json", "cases.json", "record.json"] {
            std::fs::copy(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(data), dir.join(data)).unwrap();
        }
    });
    let file = dir.join(format!("design-{}.ink", FILES.fetch_add(1, Ordering::Relaxed)));
    std::fs::write(&file, src).unwrap();
    let words: Vec<String> = ["require", file.to_str().unwrap()].into_iter().map(String::from).chain(params.iter().flat_map(|p| ["--param".to_string(), p.to_string()])).collect();
    ink_ratings::load(&CORE, &command::Args::parse(&words))
}

fn load(src: &str) -> Result<ModelView, InkError> {
    load_with(src, &[])
}

const BALL: &str = "(shape (- (+ (* x x) (* y y) (* z z)) 1))";

/// The part, and a design placing it twice; `left`, `right` and `asks` are the lines a sketch lacks.
fn opener(left: &str, right: &str, asks: &str) -> String {
    format!(
        r#"(def spec (read-json "geared.json"))
(param pull 150N)
(body housing {BALL})
(body rotor {BALL} (attr ratings (get spec "ratings")) (frame shaft (at (move 0 0 5))) (frame tip (at)))
(assembly gearmotor (instance housing housing) (instance rotor rotor) (frame shaft (frame-of rotor shaft)))
(body door {BALL})
(assembly opener
  (instance left gearmotor {left})
  (instance right gearmotor (at (move 20 0 0))
    {right})
  (instance door door))
{asks}"#
    )
}

const LEFT: &str = "(fact radial-load shaft pull)";
const RIGHT: &str = "(fact radial-load shaft 50N) (fact torque shaft (* 30 100 1N*mm)) (fact speed shaft 15rpm)";
const ASKS: &str = "(require left-ok (within-ratings opener/left))\n(require right-ok (within-ratings opener/right))";

/// Each requirement's answer, by name.
fn answers(view: &mut ModelView) -> Vec<(String, Option<Verdict>)> {
    let (design, formulas) = design(view).unwrap();
    let tools: &[&dyn Tool] = &[&formulas];
    let mut claims = Claims::new(&design, tools);
    view.doc.requires.iter().map(|r| (r.name.clone(), claims.answer(&about(&r.name, "")).held)).collect()
}

#[test]
fn a_design_states_facts_where_its_parts_are_placed_and_asks_whether_they_are_within_their_ratings() {
    let mut view = load(&opener(LEFT, RIGHT, ASKS)).unwrap();
    assert_eq!(answers(&mut view), [("left-ok".to_string(), Some(Verdict::No)), ("right-ok".to_string(), Some(Verdict::Yes))]);

    let (design, formulas) = design(&mut view).unwrap();
    // A rating written once on the part is on each instance of it, and a fact on the motor's shaft is on its rotor's.
    let names: Vec<String> = design.rules.iter().map(|r| r.about.to_string()).collect();
    assert_eq!(names, ["left-ok", "side-load on opener/left/rotor", "window on opener/left/rotor", "right-ok", "side-load on opener/right/rotor", "window on opener/right/rotor"]);
    assert!(design.rules[0].about == about("left-ok", "") && design.rules[0].says == Says::Within("opener/left".into()) && design.rules[0].from.ends_with(".ink:12"));
    assert_eq!(design.rules[1].says, Says::Ordered(vec![Term::Size(about("radial-load", "opener/left/rotor shaft")), Term::Number(100.0, "N".into())]));
    let places: Vec<String> = design.facts.iter().map(|f| f.about.to_string()).collect();
    assert_eq!(
        places,
        ["radial-load on opener/left/rotor shaft", "radial-load on opener/right/rotor shaft", "torque on opener/right/rotor shaft", "speed on opener/right/rotor shaft"]
    );

    // Every number is in ink's unit of its quantity: the param's 150 N, 15 rpm as 90 deg/s, the rating's 1 N*m as 1000 N*mm.
    let fact = |i: usize| (design.facts[i].value.unwrap(), &design.facts[i].unit[..], design.facts[i].from.rsplit(':').next().unwrap());
    assert_eq!([fact(0), fact(1), fact(2), fact(3)], [(150.0, "N", "8"), (50.0, "N", "10"), (3000.0, "N*mm", "10"), (90.0, "deg/s", "10")]);
    assert_eq!(
        design.rules[5].says,
        Says::Ordered(vec![
            Term::Number(1000.0, "N*mm".into()),
            Term::Size(about("torque", "opener/right/rotor shaft")),
            Term::Read(about("window's limit", "opener/right/rotor /1"))
        ])
    );

    // The story: the left motor's load is over its rating; the right one is under its curve between two points, an estimate.
    let tools: &[&dyn Tool] = &[&formulas];
    let mut claims = Claims::new(&design, tools);
    let side = claims.answer(&about("side-load", "opener/left/rotor"));
    assert_eq!((side.held, &side.terms[..], &side.from[..]), (Some(Verdict::No), &[("radial-load on opener/left/rotor shaft".to_string(), 150.0)][..], "a made-up gearbox table"));
    assert_eq!(side.strength, None, "as stated: nothing proved the load");
    assert_eq!(claims.answer(&about("window", "opener/left/rotor")).held, Some(Verdict::NotAvailable(Unavailable::NoSource)));
    let window = claims.answer(&about("window", "opener/right/rotor"));
    assert_eq!((window.held, window.strength), (Some(Verdict::Yes), Some(Strength::Estimated)));
    assert_eq!(claims.answer(&about("window's limit", "opener/right/rotor /1")).band, Some((4000.0, 6000.0)));
    assert_eq!(claims.answer(&about("right-ok", "")).strength, Some(Strength::Estimated));

    // A param given another value moves the fact that reads it; one the file does not have is an error.
    let mut lighter = load_with(&opener(LEFT, RIGHT, ASKS), &["pull=80N"]).unwrap();
    assert_eq!(answers(&mut lighter)[0].1, Some(Verdict::NotAvailable(Unavailable::NoSource)), "its side load is now yes; its window lacks a torque");
    let (d, _) = ink_ratings::design(&mut lighter).unwrap();
    assert_eq!(Claims::new(&d, &[]).answer(&about("side-load", "opener/left/rotor")).held, Some(Verdict::Yes));
    for misspelled in ["pul=80N", "pul=80N:90N", "pul=80N:80N"] {
        assert!(load_with(&opener(LEFT, RIGHT, ASKS), &[misspelled]).err().unwrap().message.contains("unknown param 'pul'"), "{misspelled}");
    }
}

/// A tool that must never be asked.
struct Never;

impl Tool for Never {
    fn ask(&self, about: &About, _: &mut Claims) -> Option<Fact> {
        panic!("asked for {about}")
    }
}

#[test]
fn a_design_that_asks_nothing_runs_nothing() {
    // The sketch, and the sketch with facts but no requirement: nothing is collected, nothing is asked.
    for (left, right) in [("", ""), (LEFT, RIGHT)] {
        let mut view = load(&opener(left, right, "")).unwrap();
        let (design, formulas) = design(&mut view).unwrap();
        assert_eq!((design.facts.len(), design.rules.len(), formulas.0.len()), (0, 0, 0));
        let tools: &[&dyn Tool] = &[&Never];
        let mut claims = Claims::new(&design, tools);
        assert_eq!(view.doc.requires.iter().map(|r| claims.answer(&about(&r.name, ""))).count(), 0);
    }
    // A requirement over a part with no ratings is not a yes, and still asks no tool.
    let mut view = load(&opener("", "", "(require door-ok (within-ratings opener/door))")).unwrap();
    assert_eq!(answers(&mut view), [("door-ok".to_string(), Some(Verdict::NotAvailable(Unavailable::Empty)))]);
}

#[test]
fn a_mistake_in_a_fact_or_a_requirement_is_an_error_or_a_refusal_where_it_is_written() {
    let error = |left: &str, right: &str, asks: &str| {
        let e = load(&opener(left, right, asks)).err().unwrap();
        (e.line, e.message)
    };
    // A frame the part does not have, at the fact's line, naming the one meant.
    let (line, said) = error("(fact radial-load shaf 150N)", "", ASKS);
    assert!(line == 8 && said.starts_with("'gearmotor' has no frame or instance 'shaf'") && said.contains("shaft"), "{line} {said}");
    let (line, said) = error("", "(fact radial-load 1 2 3)", ASKS);
    assert!(line == 10 && said.contains("wrong arguments"), "{line} {said}");
    let (line, said) = error("", "(fact radial-load shaft (* 2 (/ t 1s)))", ASKS);
    assert!(line == 10 && said.contains("not a field of x y z or t"), "{line} {said}");
    // A requirement of a part that is not there, or asking something else.
    let (line, said) = error("", "", "(require ok (within-ratings opener/lef))");
    assert!(line == 12 && said.contains("no instance 'lef'"), "{line} {said}");
    // What is no question of ratings is a claim, read as an assert reads it; ratings take no clause.
    let (line, said) = error("", "", "(require ok (< 1 2))");
    assert!(line == 12 && said.contains("malformed claim '<'"), "{line} {said}");
    let (line, said) = error("", "", "(require ok (within-ratings opener) (res 1))");
    assert!(line == 12 && said.contains("takes no other clause"), "{line} {said}");
    // What it leaves out on purpose is under the path asked, and is there.
    let (line, said) = error("", "", "(require ok (within-ratings opener/left (excluding opener/right)))");
    assert!(line == 12 && said.contains("opener/right is not under opener/left"), "{line} {said}");
    let (line, said) = error("", "", "(require ok (within-ratings opener (excluding opener/lef)))");
    assert!(line == 12 && said.contains("no instance 'lef'"), "{line} {said}");
    let (line, said) = error("", "", "(require ok (value 2 > 1))\n(require ok (value 2 > 1))");
    assert!(line == 13 && said.contains("require 'ok' already exists"), "{line} {said}");
    // The word before 0.2 is refused, with the form to write.
    let (line, said) = error("", "", "(assert ok (value 2 > 1))");
    assert!(line == 12 && said.contains("`assert` is now `require`"), "{line} {said}");
    let (line, said) = error("", "", &format!("{ASKS}\n(require left-ok (within-ratings opener))"));
    assert!(line == 14 && said.contains("already exists"), "{line} {said}");

    // A load said as a length is refused against a rating in newtons, and the refusal's story names the line.
    let mut view = load(&opener("(fact radial-load shaft 50mm)", "", ASKS)).unwrap();
    assert_eq!(answers(&mut view)[0].1, Some(Verdict::NotAvailable(Unavailable::Units)));
    let (d, _) = design(&mut view).unwrap();
    assert!(d.facts[0].unit == "mm" && d.facts[0].from.ends_with(".ink:8"));
    // A frame of an inner part is not the placed part's, but is reached through it; so is the inner part.
    let (line, said) = error("", "(fact radial-load tip 50N)", ASKS);
    assert!(line == 10 && said.starts_with("'gearmotor' has no frame or instance 'tip'"), "{line} {said}");
    let (line, said) = error("", "(fact radial-load rotor/shaf 50N)", ASKS);
    assert!(line == 10 && said.starts_with("'rotor' has no frame or instance 'shaf'"), "{line} {said}");
    let mut view = load(&opener("", "(fact radial-load rotor/shaft 50N) (fact mass rotor 20g) (fact hot rotor/tip 1)", ASKS)).unwrap();
    assert_eq!(answers(&mut view)[1].1, Some(Verdict::NotAvailable(Unavailable::NoSource)), "its side load is yes; its window still lacks a torque");
    let places: Vec<String> = design(&mut view).unwrap().0.facts.iter().map(|f| f.about.to_string()).collect();
    assert_eq!(places, ["radial-load on opener/right/rotor shaft", "mass on opener/right/rotor", "hot on opener/right/rotor tip"]);
    // A frame is the end of a path, not a step on it.
    let (line, said) = error("", "(fact radial-load shaft/rotor 50N)", ASKS);
    assert!(line == 10 && said.starts_with("'gearmotor' has no frame or instance 'shaft'"), "{line} {said}");
    // A fact on the whole design's side, written where the motor is placed in it (a part two levels down).
    let nested = format!("{}\n(assembly plant (instance line opener (fact radial-load left/shaft 99N)))\n(require plant-ok (within-ratings plant/line/left))", opener("", "", ""));
    let mut view = load(&nested).unwrap();
    let (d, _) = design(&mut view).unwrap();
    assert_eq!(d.facts[0].about, about("radial-load", "plant/line/left/rotor shaft"));
    assert_eq!(Claims::new(&d, &[]).answer(&about("side-load", "plant/line/left/rotor")).held, Some(Verdict::Yes));
    // A fact under another name, or on the part where the rating reads its shaft, leaves the rating unanswered: never a pass.
    for typo in ["(fact radial-lode shaft 50N)", "(fact radial-load 50N)"] {
        let mut view = load(&opener("", typo, ASKS)).unwrap();
        assert_eq!(answers(&mut view)[1].1, Some(Verdict::NotAvailable(Unavailable::NoSource)), "{typo}");
    }
    // The same fact stated twice with different numbers is not settled.
    let mut view = load(&opener("(fact radial-load shaft 50N) (fact radial-load shaft 60N)", "", ASKS)).unwrap();
    assert_eq!(answers(&mut view)[0].1, Some(Verdict::NotSettled(Unsettled::Disagree)));
}

/// One part carrying the rating `key` of data/cases.json, with `facts` stated on it, after `before`.
fn case(key: &str, before: &str, facts: &str) -> String {
    format!(
        "(def spec (read-json \"cases.json\"))\n{before}\n(body part {BALL} (attr ratings (get spec \"{key}\")) (frame shaft (at)) (frame tip (at (move 0 0 1))))\n(assembly a (instance p part\n {facts}))\n(require ok (within-ratings a))"
    )
}

/// The requirement's answer, and the first fact as the design hands it on.
fn asked(src: &str, params: &[&str]) -> (Verdict, Fact) {
    let mut view = load_with(src, params).unwrap_or_else(|e| panic!("{}", e.message));
    let (design, formulas) = design(&mut view).unwrap();
    let tools: &[&dyn Tool] = &[&formulas];
    (Claims::new(&design, tools).answer(&about("ok", "")).held.unwrap(), design.facts[0].clone())
}

const UNITS: Verdict = Verdict::NotAvailable(Unavailable::Units);

/// A fact has the unit of what its value is, by ink's one arithmetic: a number written with its unit, a name, a
/// sum, and a product or quotient, whose units make the quantity they make. With no unit, or another's, it is
/// refused against a rating with its line; it is never read as if it had the rating's.
#[test]
fn a_fact_has_the_unit_of_what_it_is_or_is_refused() {
    let before = "(param pull 3N)\n(param arm 2mm)\n(def turned (* pull arm))\n(def warm 300K)";
    let load = |value: &str, key: &str, what: &str| {
        let (held, fact) = asked(&case(key, before, &format!("(fact {what} shaft {value})")), &[]);
        (held, fact.value.unwrap(), fact.unit, fact.notes.len(), fact.from.ends_with(".ink:8"))
    };
    let (yes, no) = (Verdict::Yes, Verdict::No);
    // Against 10 N either way.
    for (value, is, held) in [
        ("pull", 3.0, yes),
        ("2.5N", 2.5, yes),
        ("(* 2 pull)", 6.0, yes),
        ("(* 4 1N)", 4.0, yes),
        ("(+ pull 5N)", 8.0, yes),
        ("(- (* 3 pull) 1lbf)", 9.0 - 4.4482216152605, yes),
        ("(neg pull)", -3.0, yes),
        ("(max pull 4N)", 4.0, yes),
        ("(/ (+ pull 5N) 2)", 4.0, yes),
        // A torque over its arm is a force again; a ratio of lengths is a pure number; kg under m/s2 is newtons.
        ("(/ (* pull arm) arm)", 3.0, yes),
        ("(/ turned 1mm)", 6.0, yes),
        ("(* (/ 2in 4mm) 1N)", 12.7, no),
        ("(* 2kN (/ 4mm 1m))", 8.0, yes),
        ("(* 1kg 9.81m/s2)", 9.81, yes),
        ("(* 10kg 9.81m/s2)", 98.1, no),
        ("(- 1000N (* 1kg 1m/s2))", 999.0, no),
        ("(max 0N (* 1kg 1m/s2))", 1.0, yes),
        ("(abs (* 10kg 9.81m/s2))", 98.1, no),
        ("(neg (* 10kg 9.81m/s2))", -98.1, no),
        ("(* 2MPa 4mm2)", 8.0, yes),
        // However the product is reached: reduced, folded, applied.
        ("(reduce * (vec 10kg 9.81m/s2))", 98.1, no),
        ("(reduce * (vec 1kg 9.81m/s2))", 9.81, yes),
        ("(fold * 1 (vec 1kg 9.81m/s2))", 9.81, yes),
        ("(apply * (vec 1kg 9.81m/s2))", 9.81, yes),
    ] {
        let (got, n, unit, notes, _) = load(value, "sized", "load");
        assert!(got == held && same(n, is) && unit == "N" && notes == 0, "{value}: {got:?} {n} {unit}");
    }
    assert_eq!(load("7162rpm", "fast", "speed").0, no);
    // A quantity that is not the rating's is refused, saying which it is.
    for (value, key, what, is) in [
        ("50mm", "sized", "load", "mm"),
        ("(* 2mm 2mm)", "sized", "load", "mm2"),
        ("(* pull arm)", "sized", "load", "N*mm"),
        ("turned", "sized", "load", "N*mm"),
        // A length scaled, then times a force: a torque, not newtons.
        ("(* (reduce * (vec 3mm 2)) 1N)", "sized", "load", "N*mm"),
        // A torque times a turning rate is a power, its angle in radians; a number over a turning rate a time.
        ("(* 1N*m 60rpm)", "sized", "load", "W"),
        ("(/ 360 60rpm)", "fast", "speed", "s"),
    ] {
        let (held, _, unit, notes, at) = load(value, key, what);
        assert!(held == UNITS && unit == is && notes == 0 && at, "{value}: {held:?} {unit:?}");
    }
    // A speed over a length is a per-second: a frequency, never read as a turning rate (a cycle is not a radian),
    // and a turning rate never as a frequency. Each is refused, and says what to write; written so, it is read.
    for (value, key, what, is) in [("(/ 3000mm/s 4mm)", "fast", "speed", "Hz"), ("(/ 1W 1N*m)", "fast", "speed", "Hz"), ("60rpm", "humming", "hum", "deg/s")] {
        let design = design(&mut crate::load(&case(key, before, &format!("(fact {what} shaft {value})"))).unwrap()).unwrap().0;
        let (held, _, unit, notes, _) = load(value, key, what);
        assert!(held == UNITS && unit == is && notes == 1 && design.facts[0].notes[0].contains("a cycle is not a radian"), "{value}: {held:?} {unit:?}");
    }
    assert_eq!([load("(* 360deg (/ 3000mm/s 4mm))", "fast", "speed").0, load("(/ 60rpm 360deg)", "humming", "hum").0, load("(/ 20rpm 360deg)", "humming", "hum").0], [no, no, yes]);
    // A root that halves an angle's power (a spring's stiffness over an inertia: 1 rad/s, 0.159 Hz) is no quantity:
    // refused, never read as 1 Hz.
    let (held, _, unit, notes, _) = load("(sqrt (/ (/ 1N*m 1rad) 1kg*m2))", "humming", "hum");
    assert!(held == UNITS && unit.is_empty() && notes == 1, "{held:?} {unit}");
    // A bare number, and a product that is no quantity of ink's (a force cubed, a torque times a turning rate, which
    // would be a power only in radians), have no unit: refused, with why.
    // Nor has what an operation dropped the unit of: it is not a pure number again, to be given newtons.
    let bare = [
        ("4", "sized", "load"),
        ("(* 2 2)", "sized", "load"),
        ("(* (pow 2N 2) 1N)", "sized", "load"),
        ("(* (floor (* 2N 3mm)) 1N)", "sized", "load"),
        ("(* (sin 30mm) 1N)", "sized", "load"),
        ("(* (exp 2mm) 1N)", "sized", "load"),
    ];
    for (value, key, what) in bare {
        let (held, _, unit, notes, at) = load(value, key, what);
        assert!(held == UNITS && unit.is_empty() && notes == 1 && at, "{value}: {held:?} {unit:?}");
    }
    // A sum of two quantities, and a temperature multiplied, are errors where they are written.
    for (value, says) in [
        ("(min 50N (* 100mm 100mm))", "'min' takes values of one quantity"),
        ("(+ 1N 1mm (* 2N 3s))", "'+' takes values of one quantity"),
        ("(+ 5N (* 60rpm 1s))", "'+' takes values of one quantity"),
        ("(* 2 20degC)", "a temperature is a point on its scale"),
        ("(* 2 300K)", "a temperature is a point on its scale"),
        ("(* 2 warm)", "a temperature is a point on its scale"),
    ] {
        let e = crate::load(&case("sized", before, &format!("(fact load shaft {value})"))).err().unwrap();
        assert!(e.line == 8 && e.message.contains(says), "{value}: {}", e.message);
    }
    // A torque times a turning rate is a power, in watts, held against a rating in watts: 1 N m at 1 rad/s is 1 W,
    // at 60 rpm 6.28 W. And an angle a second apart from a turning rate is one.
    let (held, watts, unit, ..) = load("(* 1N*m 1rad/s)", "powered", "power");
    assert!(held == yes && same(watts, 1.0) && unit == "W", "{watts} {unit}");
    let (held, watts, ..) = load("(* 1N*m 60rpm)", "powered", "power");
    assert!(held == no && (watts - std::f64::consts::TAU).abs() < 1e-9, "{watts}");
    assert_eq!(load("(/ 36000deg 1s)", "fast", "speed").0, yes);
    // However the product is grouped and reached, it is one fact: a turning rate, 360 deg/s, and none at rest.
    let before = &format!("{before}\n(def linear (* 2mm 60rpm))\n(def stopped (* 2mm 0rpm))\n(fn multiply (a b) (* a b))");
    let grouped = |value: &str, key: &str, what: &str| {
        let (held, fact) = asked(&case(key, before, &format!("(fact {what} shaft {value})")), &[]);
        (held, fact.value.unwrap(), fact.unit, fact.band)
    };
    for value in ["(* 2mm 60rpm (/ 1 2mm))", "(* (* 2mm 60rpm) (/ 1 2mm))", "(* linear (/ 1 2mm))", "(multiply linear (/ 1 2mm))", "(fold * 1 (vec 2mm 60rpm (/ 1 2mm)))"] {
        let (held, rate, unit, _) = grouped(value, "fast", "speed");
        assert!(held == yes && same(rate, 360.0) && unit == "deg/s", "{value}: {held:?} {rate} {unit}");
    }
    assert_eq!(grouped("(* stopped (/ 1 2mm))", "fast", "speed"), (yes, 0.0, "deg/s".to_string(), None));
    for value in ["(* (/ 1rad 1mm) 2mm2 1rad/s)", "(* (/ 1rad 1mm) (* 2mm2 1rad/s))"] {
        let (held, speed, unit, _) = grouped(value, "quick", "pace");
        assert!(held == yes && (speed - 2.0).abs() < 1e-12 && unit == "mm/s", "{value}: {held:?} {speed} {unit}");
    }
    // A range of angles times a radius is a range of arcs, in radians: 0 to pi mm.
    let (held, _, unit, band) = grouped("(* 2mm (ranged 0deg 1.5707963267948966rad))", "long", "reach");
    assert!(held == yes && unit == "mm" && band.is_some_and(|b| b.0 == 0.0 && (b.1 - std::f64::consts::PI).abs() < 1e-12), "{held:?} {unit} {band:?}");
    assert_eq!(load("(* (/ 124200deg 1397mm) 500mm/s)", "fast", "speed").0, no);
    // A temperature as written or named is read. A difference of two is an amount, not a point: it scales, and is
    // refused against a rating that is a point on the scale.
    assert_eq!([load("300K", "warm", "temp").0, load("warm", "warm", "temp").0, load("450K", "warm", "temp").0], [yes, yes, no]);
    let (held, kelvins, unit, notes, _) = load("(* 2 (- 300K 290K))", "warm", "temp");
    assert_eq!((held, kelvins, unit.as_str(), notes), (UNITS, 20.0, "", 1));
}

/// A rating on sizes reads a load, and a speed, without its sign; a signed one keeps it.
#[test]
fn a_rating_on_sizes_ignores_the_sign() {
    let held = |key: &str, facts: &str| asked(&case(key, "", facts), &[]).0;
    assert_eq!(
        [held("sized", "(fact load shaft -500N)"), held("sized", "(fact load shaft 500N)"), held("sized", "(fact load shaft -5N)")],
        [Verdict::No, Verdict::No, Verdict::Yes]
    );
    assert_eq!([held("fast", "(fact speed shaft -7000rpm)"), held("fast", "(fact speed shaft -5000rpm)")], [Verdict::No, Verdict::Yes]);
    assert_eq!(
        [held("signed", "(fact load shaft -15N)"), held("signed", "(fact load shaft -25N)"), held("signed", "(fact load shaft 15N)")],
        [Verdict::Yes, Verdict::No, Verdict::No]
    );
    // The line is read at the speed's size: 0.5 N left at 5 rpm either way round.
    assert_eq!([held("line", "(fact load shaft -0.4N) (fact speed shaft -5rpm)"), held("line", "(fact load shaft -0.6N) (fact speed shaft -5rpm)")], [Verdict::Yes, Verdict::No]);
}

/// A value that is no number is no fact: refused, with its line and what it was.
#[test]
fn a_fact_that_is_no_number_is_refused() {
    let values = [
        "(/ pull 0)",
        "(* 1e200 1e200 pull)",
        "(* -1e200 1e200 pull)",
        "(/ pull (- pull pull))",
        "(/ 1N (- 2 2))",
        "(* (/ 1 0) 1N)",
        "(* (ln -1) 1N)",
        "(* (pow -2 0.5) 1N)",
        "(* (sqrt -1) 1N)",
        "(* (sqrt (- 1 (/ pull 1N))) 1N)",
        // An overflow on the way, gone by the end; a quotient by what may be nothing, times nothing; and each under a name.
        "(* (/ 1e200 (* 1e200 1e200)) 1e200 1000N)",
        "(* 0 (/ 1 slack) 1N)",
        "(* (/ 1 slack) 1N)",
        "(- (/ pull slack) (/ pull slack))",
        "(* 0 bad pull)",
        "(* gone pull)",
        "(+ pull (* 0 (sqrt (- (/ pull 1N) 4)) 1N))",
    ];
    for value in values {
        let before = "(param pull 3N)\n(def slack (ranged -1 1))\n(def bad (/ 1 0))\n(def gone (/ 1e200 (* 1e200 1e200)))\n(def lost-root (sqrt -4))";
        let (held, fact) = asked(&case("sized", before, &format!("(fact load shaft {value})")), &[]);
        assert_eq!((held, fact.value, fact.held), (Verdict::NotAvailable(Unavailable::Outside), None, Some(Verdict::NotAvailable(Unavailable::Outside))), "{value}");
        assert!(fact.from.ends_with(".ink:9") && fact.notes[0].starts_with("its value is no number: "), "{value}: {fact:?}");
    }
}

/// What ink itself defines, under a name of its own: the sqrt of a negative is 0, and nothing times anything is
/// nothing. A fact that reads such a name reads ink's number; written in the fact, both are refused (above).
#[test]
fn a_name_holds_the_number_ink_gave_it() {
    let before = "(def slack (ranged -1 1))\n(def lost-root (sqrt -4))\n(def dropped (* 0 (/ 1 slack)))";
    let read = |value: &str| {
        let (held, fact) = asked(&case("sized", before, &format!("(fact load shaft {value})")), &[]);
        (held, fact.value, fact.unit)
    };
    // The root is 0 to ink, and 3 N and nothing is 3 N.
    assert_eq!(read("(+ 3N (* lost-root 1N))"), (Verdict::Yes, Some(3.0), "N".to_string()));
    assert_eq!(read("(+ 3N (* dropped 1N))"), (Verdict::Yes, Some(3.0), "N".to_string()));
}

/// A fact that reads a ranged param is a band over its range: a tolerance, a `(ranged lo hi)`, a range given.
#[test]
fn a_range_is_a_band() {
    let held = |before: &str, value: &str, params: &[&str]| {
        let (held, fact) = asked(&case("sized", before, &format!("(fact load shaft {value})")), params);
        (held, fact.band)
    };
    let straddles = Verdict::NotSettled(Unsettled::Straddles);
    assert_eq!(held("(def force (ranged 0N 18N))", "force", &[]), (straddles, Some((0.0, 18.0))));
    assert_eq!(held("(def force (ranged 0N 8N))", "force", &[]), (Verdict::Yes, Some((0.0, 8.0))));
    assert_eq!(held("(def force (ranged 11N 18N))", "(* 2 force)", &[]), (Verdict::No, Some((22.0, 36.0))));
    assert_eq!(held("(param force 9N (tol -1N 2N))", "force", &[]), (straddles, Some((8.0, 11.0))));
    assert_eq!(held("(param force 5N)", "force", &["force=0N:18N"]), (straddles, Some((0.0, 18.0))));
    assert_eq!(held("(param force 5N)", "force", &["force=-18N:-12N"]), (Verdict::No, Some((-18.0, -12.0))), "a size: 12 to 18");
    assert_eq!(held("(param force 5N)", "force", &[]), (Verdict::Yes, None));
    // Two ranges, one passed as an assembly's argument: 60 to 120 and 0 to 60.
    let two = "(def spec (read-json \"cases.json\"))\n(def first (ranged 30N 60N))\n(param second 30N (tol -30N 30N))\n(def twice (* 2 first))\n(body part (shape (- (+ (* x x) (* y y) (* z z)) 1)) (attr ratings (get spec \"sized\")) (frame shaft (at)))\n(assembly (module (f 0N)) (instance p part (fact load shaft (+ f second))))\n(assembly a (instance m (module twice)))\n(require ok (within-ratings a))";
    assert_eq!(asked(two, &[]).1.band, Some((60.0, 180.0)));
    assert_eq!((asked(two, &["second=65N:70N"]).0, asked(two, &["second=65N:70N"]).1.band), (Verdict::No, Some((125.0, 190.0))));
}

/// Two tables in one rating are two limits: 5 N at least and 10 N at most.
#[test]
fn two_tables_in_one_rating_are_two_limits() {
    let held = |load: &str| asked(&case("between", "", &format!("(fact load shaft {load}) (fact speed shaft 50rpm)")), &[]).0;
    assert_eq!([held("7N"), held("4N"), held("11N")], [Verdict::Yes, Verdict::No, Verdict::No]);
    let mut view = load(&case("between", "", "")).unwrap();
    let (_, formulas) = design(&mut view).unwrap();
    assert_eq!(formulas.0.iter().map(|f| f.about.to_string()).collect::<Vec<_>>(), ["r's limit on a/p /1", "r's limit on a/p /2"]);
}

/// A frame an assembly exports is the frame it is written as, whatever it is called.
#[test]
fn an_exported_frame_is_the_frame_it_is_written_as() {
    let design_of = |exported: &str| {
        let src = format!(
            "(def spec (read-json \"cases.json\"))\n{}\n(body part {BALL} (attr ratings (get spec \"sized\")) (frame shaft (at)) (frame tip (at (move 0 0 1))))\n(assembly unit (instance r part) {exported})\n(assembly a (instance u unit (fact load shaft 50N)))\n(require ok (within-ratings a))",
            if exported.contains("my-own") { "(fn frame-of (part frame) (at (move 200 0 0)))\n(def r 0)\n(def shaft 0)" } else { "" }
        );
        let mut view = load(&src).unwrap();
        (answers(&mut view)[0].1.unwrap(), design(&mut view).unwrap().0.facts[0].about.on.clone())
    };
    // Its own name for the part's tip: the load is on the tip, where no rating reads it. Never on r/shaft.
    assert_eq!(design_of("(frame shaft (frame-of r tip))"), (Verdict::NotAvailable(Unavailable::NoSource), "a/u/r tip".to_string()));
    assert_eq!(design_of("(frame shaft (frame-of r shaft))"), (Verdict::No, "a/u/r shaft".to_string()));
    assert_eq!(design_of("(frame out (frame-of r shaft)) (frame shaft (frame-of r tip))").1, "a/u/r tip");
    // A frame of its own, placed off the part's, is the assembly's own.
    assert_eq!(design_of("(frame shaft (at (move 0 0 3) (frame-of r shaft)))").1, "a/u shaft");
    assert_eq!(design_of("(frame shaft (at (move 0 0 3)))").1, "a/u shaft");
    // A function of the file's own called frame-of is not the language's: the frame it gives is the assembly's own.
    assert_eq!(design_of("(frame shaft (frame-of r shaft)) (frame my-own (at))"), (Verdict::NotAvailable(Unavailable::NoSource), "a/u shaft".to_string()));

    // A frame and an instance of one name do not share a place: a fact that could be on either is refused where
    // it is written, and the frame's place is not the instance's.
    let both = |fact: &str| {
        format!(
            "(def spec (read-json \"cases.json\"))\n(body part {BALL} (attr ratings (get spec \"plain\")))\n(assembly inner (instance p part) (frame p (at (move 20 0 0))) (frame q (at)))\n(assembly a (instance i inner {fact}))\n(require ok (within-ratings a/i/p))"
        )
    };
    let e = load(&both("(fact load p 25N)")).err().unwrap();
    assert!(e.line == 4 && e.message.contains("has a frame and an instance named 'p'"), "{}", e.message);
    let mut view = load(&both("(fact load q 25N)")).unwrap();
    assert_eq!((answers(&mut view)[0].1.unwrap(), design(&mut view).unwrap().0.facts[0].about.on.as_str()), (Verdict::NotAvailable(Unavailable::NoSource), "a/i q"));

    // A design that defines a `fact` of its own (a function giving a frame) places its instance by it, as before.
    let own = format!("(fn fact (distance) (at (move distance 0 0)))\n(body part {BALL})\n(assembly rig (instance p part (fact 20)))");
    let view = load(&own).unwrap_or_else(|e| panic!("{}", e.message));
    assert!(view.doc.assembly("rig").unwrap().instances.iter().all(|i| i.facts.is_empty()));
}

/// Table points nearer than rounding tells apart, a rating named like a table's limit, and a line too steep to subtract.
#[test]
fn a_table_keeps_its_band_its_name_and_its_numbers() {
    let held = |key: &str, facts: &str| asked(&case(key, "", facts), &[]).0;
    // The speed is at each of three points for all rounding tells: the limit is between the least and greatest, 10 to 100.
    let mut view = load(&case("near", "", "(fact speed shaft 1.0000000000000002deg/s) (fact load shaft 50N)")).unwrap();
    let (d, formulas) = design(&mut view).unwrap();
    let tools: &[&dyn Tool] = &[&formulas];
    let mut claims = Claims::new(&d, tools);
    assert_eq!(claims.answer(&about("r's limit", "a/p /1")).band, Some((10.0, 100.0)));
    assert_eq!(claims.answer(&about("ok", "")).held, Some(Verdict::NotSettled(Unsettled::Straddles)));
    assert_eq!([held("near", "(fact speed shaft 1deg/s) (fact load shaft 5N)"), held("near", "(fact speed shaft 1deg/s) (fact load shaft 500N)")], [Verdict::Yes, Verdict::No]);
    // Ratings named as a table's limit might be do not hide the table: 50 N is over the second table's 10 N.
    let both = |load: &str| held("named-like-a-limit", &format!("(fact speed shaft 1rpm) (fact load shaft {load})"));
    assert_eq!([both("50N"), both("5N")], [Verdict::No, Verdict::Yes]);
    // A line from -1e308 to 1e308 is 0 half way; nothing overflows on the way.
    assert_eq!([held("steep", "(fact speed shaft 1rpm) (fact load shaft -5N)"), held("steep", "(fact speed shaft 1rpm) (fact load shaft 5N)")], [Verdict::Yes, Verdict::No]);
}

/// Numbers equal but for the rounding of a unit or of a table's arithmetic are equal.
#[test]
fn rounding_is_not_a_difference() {
    let held = |key: &str, facts: &str| asked(&case(key, "", facts), &[]).0;
    assert_eq!(
        [held("kilo", "(fact load shaft 1001N)"), held("kilo", "(fact load shaft 1001.001N)"), held("kilo", "(fact load shaft 1.001kN)")],
        [Verdict::Yes, Verdict::No, Verdict::Yes]
    );
    assert_eq!([held("line", "(fact load shaft 0.1N) (fact speed shaft 9rpm)"), held("line", "(fact load shaft 0.1001N) (fact speed shaft 9rpm)")], [Verdict::Yes, Verdict::No]);
    assert_eq!(held("line", "(fact load shaft 0N) (fact speed shaft 10rpm)"), Verdict::Yes, "at the line's end");
    // Four floating-point steps and no more: 100.00000000000009 N is six above 100 N.
    assert_eq!(
        [
            held("plain-shaft", "(fact load shaft 100.00000000000009N)"),
            held("plain-shaft", "(fact load shaft 100.00000000000003N)"),
            held("plain-shaft", "(fact load shaft 100.0000001N)")
        ],
        [Verdict::No, Verdict::Yes, Verdict::No]
    );
}

#[test]
fn data_that_is_not_a_rating_is_an_error_not_a_rating_skipped() {
    let bad = |ratings: &str| {
        let part = format!("(def spec (read-json \"bad.json\"))\n(body part {BALL} (attr ratings {ratings}))\n(assembly a (instance p part))\n(require ok (within-ratings a))");
        design(&mut load(&part).unwrap()).err().unwrap()
    };
    assert_eq!(bad("7"), "a/p: its ratings are not data");
    assert_eq!(bad(r#"(data (r (from "x")))"#), "a/p: its ratings are not an object by name");
    for (key, says) in [
        ("no-ordered", "no `ordered`"),
        ("no-from", "no `from`"),
        ("no-unit", "a term is a number with its unit, a read, or a table at a read"),
        ("bad-unit", "'N.m' is no unit of ink's"),
        ("no-points", "a table has points"),
        ("bad-point", "a point is [x, y]"),
        ("no-units", "a table names its x and y units"),
        ("bad-table-unit", "'Nm' is no unit of ink's"),
        ("table-at-number", "a term is a number with its unit, a read, or a table at a read"),
        (
            "too-big",
            "100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 kN is no number in ink's unit (N)",
        ),
    ] {
        assert_eq!(bad(&format!("(get spec \"{key}\")")), format!("a/p: rating 'r': {says}"));
    }
}

/// A part works out facts from other facts about it (`(attr derives (vec (vec NAME FN UNIT ...) ...))`), when
/// they are asked for: FN of the facts its parameters name, each as the core gives it and in the unit written for
/// it, in ink's own arithmetic. Here: `double` is twice the load, `work` is that times the reach, held under 100 N mm.
#[test]
fn a_part_works_out_facts_from_the_facts_about_it() {
    let part = |key: &str, derives: &str| {
        format!("(def spec (read-json \"cases.json\"))\n(param push 10N)\n(body part {BALL} (attr ratings (get spec \"{key}\"))\n (attr derives {derives}))")
    };
    let derives = "(vec (vec \"double\" (lambda (load) (* 2 load)) 1N)\n (vec \"work\" (lambda (double reach) (* double reach)) 1N 1mm))";
    let placed = |facts: &str| format!("{}\n(assembly a (instance p part {facts}))\n(require ok (within-ratings a))", part("worked", derives));
    // What the requirement answers, and the fact `what` on `on`, of the design `src`.
    let asked = |src: &str, params: &[&str], what: &str, on: &str| {
        let mut view = load_with(src, params).unwrap_or_else(|e| panic!("{}", e.message));
        let (design, formulas) = design(&mut view).unwrap();
        let view = std::cell::RefCell::new(view);
        let derived = Derived::new(&view).unwrap();
        let tools: &[&dyn Tool] = &[&formulas, &derived];
        let mut claims = Claims::new(&design, tools);
        (claims.answer(&about("ok", "")).held.unwrap(), claims.answer(&about(what, on)))
    };
    let ask = |facts: &str, params: &[&str]| asked(&placed(facts), params, "work", "a/p");
    // 2 x 10 N x 4 mm = 80 N mm: yes. It is written where its function is, and rests on what the function read.
    let (held, work) = ask("(fact load push) (fact reach 4mm)", &[]);
    assert_eq!((held, work.value, work.unit.as_str(), work.band), (Verdict::Yes, Some(80.0), "N*mm", None));
    assert_eq!(work.used, [about("double", "a/p"), about("reach", "a/p")]);
    assert!(work.from.ends_with(".ink:5"), "{}", work.from);
    // Over the limit: no. A range of loads is a range of work, and one across the limit is not settled.
    assert_eq!(ask("(fact load push) (fact reach 6mm)", &[]).0, Verdict::No);
    let (held, work) = ask("(fact load push) (fact reach 4mm)", &["push=10N:15N"]);
    assert_eq!((held, work.band), (Verdict::NotSettled(Unsettled::Straddles), Some((80.0, 120.0))));
    // A fact it reads that is not stated is not assumed, and one stated two ways is not chosen between: it is open,
    // for that reason, and so is the answer above it.
    let none = Verdict::NotAvailable(Unavailable::NoSource);
    for (facts, why) in [("(fact load push)", none), ("(fact reach 4mm)", none), ("(fact load push) (fact load 11N) (fact reach 4mm)", Verdict::NotSettled(Unsettled::Disagree))] {
        let (held, work) = ask(facts, &[]);
        assert_eq!((held, work.held, work.value), (why, Some(why), None), "{facts}");
    }
    // Each fact is read in the unit its parameter says: a load in millimetres, a bare one, and a reach that is a
    // force are refused for their units, saying which, never read as if they were newtons.
    for (facts, says) in [
        ("(fact load 10mm) (fact reach 4mm)", "its formula reads 'load' in N, and it is stated in mm"),
        ("(fact load 10) (fact reach 4mm)", "its formula reads 'load' in N, and it is stated in no unit"),
        ("(fact load 10N) (fact reach 4N)", "its formula reads 'reach' in mm, and it is stated in N"),
    ] {
        let (held, work) = ask(facts, &[]);
        let double = asked(&placed(facts), &[], "double", "a/p").1;
        assert_eq!((held, work.held, work.value), (UNITS, Some(UNITS), None), "{facts}");
        assert!(work.notes.iter().chain(&double.notes).any(|n| n.contains(says)), "{facts}: {:?} {:?}", work.notes, double.notes);
    }
    // A fact with no number makes one with none.
    let outside = Verdict::NotAvailable(Unavailable::Outside);
    let (held, work) = ask("(fact load (/ push 0)) (fact reach 4mm)", &[]);
    assert_eq!((held, work.held), (outside, Some(outside)));
    // A choice between two values with one unit has it, whichever is chosen.
    let (held, work) = ask("(fact load (select (- push 20N) push 3N)) (fact reach 4mm)", &[]);
    assert_eq!((held, work.value, work.unit.as_str()), (Verdict::Yes, Some(80.0), "N*mm"));
    // A fact the part works out is the part's to say: a design that states it too is an error, not a second source
    // that would be read in the formula's place.
    let e = design(&mut load(&placed("(fact load push) (fact reach 4mm) (fact work 70N*mm)")).unwrap()).err().unwrap();
    assert!(e.starts_with("work on a/p (") && e.contains("the part works this out from other facts"), "{e}");

    // The facts are the core's, wherever they were stated. On a part under the instance; the load from the assembly
    // above and the reach from the one below; and a load stated at both levels, two ways, which is not settled,
    // so neither is the work nor the answer: never the yes of the one nearest the part.
    let unit = |inner: &str, outer: &str| {
        format!("{}\n(assembly unit (instance q part {inner}))\n(assembly a (instance p unit {outer}))\n(require ok (within-ratings a))", part("worked", derives))
    };
    let below = |inner: &str, outer: &str| asked(&unit(inner, outer), &[], "work", "a/p/q");
    let (held, work) = below("", "(fact load q push) (fact reach q 4mm)");
    assert_eq!((held, work.value), (Verdict::Yes, Some(80.0)));
    let (held, work) = below("(fact reach 4mm)", "(fact load q push)");
    assert_eq!((held, work.value), (Verdict::Yes, Some(80.0)));
    let disagree = Verdict::NotSettled(Unsettled::Disagree);
    let (held, work) = below("(fact load 1N) (fact reach 4mm)", "(fact load q 60N)");
    assert_eq!((held, work.held, work.value), (disagree, Some(disagree), None));
    // Facts at two places are not one place's: the load on the part and the reach on the unit make no work.
    assert_eq!(below("", "(fact load q push) (fact reach 4mm)").1.held, Some(none));

    // A formula says where it holds, and why not elsewhere: (vec VALUE OK WHY). Outside, at any value read, the fact
    // is refused with WHY.
    let limited = "(vec (vec \"work\" (lambda (load reach) (vec (* load reach) (<= load 12N) \"over 12 N the part is not rated\")) 1N 1mm))";
    let only = |facts: &str, params: &[&str]| {
        asked(&format!("{}\n(assembly a (instance p part {facts}))\n(require ok (within-ratings a))", part("worked", limited)), params, "work", "a/p")
    };
    assert_eq!(only("(fact load push) (fact reach 4mm)", &[]).1.value, Some(40.0));
    for params in [&["push=15N"][..], &["push=10N:15N"][..], &["push=5N:15N"][..]] {
        let (held, work) = only("(fact load push) (fact reach 4mm)", params);
        assert_eq!((held, work.held, &work.notes[..]), (outside, Some(outside), &["over 12 N the part is not rated".to_string()][..]), "{params:?}");
    }

    // A function of nothing is a number the part states (a record's limit): as stated, resting on nothing.
    let stated = "(vec (vec \"limit\" (lambda () (* 2 50N*mm))))";
    let (_, limit) = asked(&format!("{}\n(assembly a (instance p part))\n(require ok (within-ratings a))", part("worked", stated)), &[], "limit", "a/p");
    assert_eq!((limit.value, limit.unit.as_str(), limit.strength, limit.used.len()), (Some(100.0), "N*mm", None, 0));

    // One ask leaves nothing behind for the next: the tool asked of a load with a band, then of the same design
    // with the load at one value, gives a band and then none.
    let mut view = load_with(&placed("(fact load push) (fact reach 4mm)"), &["push=10N:15N"]).unwrap();
    let (ranged, formulas) = design(&mut view).unwrap();
    let view = std::cell::RefCell::new(view);
    let derived = Derived::new(&view).unwrap();
    let tools: &[&dyn Tool] = &[&formulas, &derived];
    let mut fixed = Design { facts: ranged.facts.clone(), rules: ranged.rules.clone() };
    fixed.facts.iter_mut().for_each(|f| f.band = None);
    let work = about("work", "a/p");
    assert_eq!(Claims::new(&ranged, tools).answer(&work).band, Some((80.0, 120.0)));
    let again = Claims::new(&fixed, tools).answer(&work);
    assert_eq!((again.value, again.band), (Some(100.0), None));

    // A function that cannot take the facts in the units it says (a sum of a force and a length), one that decides
    // at load on a fact (no fact is stated then), and `derives` that is no list of names, functions and units,
    // are errors when the tool is made, naming the part, the fact and where the function is.
    let made = |derives: &str| {
        let view =
            std::cell::RefCell::new(load(&format!("{}\n(assembly a (instance p part (fact load push)))\n(require ok (within-ratings a))", part("worked", derives))).unwrap());
        Derived::new(&view).err().unwrap_or_default()
    };
    let e = made("(vec (vec \"work\" (lambda (load reach) (+ load reach)) 1N 1mm))");
    assert!(e.starts_with("a/p: 'work' (") && e.contains(".ink:4): '+' takes values of one quantity"), "{e}");
    let e = made("(vec (vec \"work\" (lambda (load) (if (> load 5N) load 0N)) 1N))");
    assert!(e.contains("a/p: 'work' (") && e.contains("decides on 'load' at load, before any fact is stated"), "{e}");
    for bad in [
        "5",
        "(vec 5)",
        "(vec (vec \"work\" 5))",
        "(vec (vec \"work\" (lambda (load) load)))",
        "(vec (vec \"work\" (lambda (load) load) 1N 1N))",
        "(vec (vec \"work\" (lambda (load rest ...) load) 1N 1N))",
        "(vec (vec \"work\" (lambda (load) load) \"N\"))",
        "(vec (vec \"work\" (lambda (load) load) 20degC))",
    ] {
        assert!(made(bad).contains("(attr derives ...) is a vec of (vec NAME FN UNIT ...)"), "{bad}: {}", made(bad));
    }
    // A rating and a fact of one name at one place, stated or worked out: one would answer for the other, so it is an error.
    let clash = format!("{}\n(assembly a (instance p part (fact load push)))\n(require ok (within-ratings a))", part("clash", "(vec)"));
    let e = design(&mut load(&clash).unwrap()).err().unwrap();
    assert!(e.contains("load on a/p: a rating and a fact stated there have this one name"), "{e}");
    let clash = format!("{}\n(assembly a (instance p part))\n(require ok (within-ratings a))", part("clash", "(vec (vec \"load\" (lambda (push) push) 1N))"));
    let e = design(&mut load(&clash).unwrap()).err().unwrap();
    assert!(e.contains("a/p: 'load' is a rating and a fact the part works out"), "{e}");
    // A design that asks nothing makes no formula: its function is never called.
    let view = std::cell::RefCell::new(
        load(&format!("{}\n(assembly a (instance p part (fact load push)))", part("worked", "(vec (vec \"work\" (lambda (load reach) (+ load reach)) 1N 1mm))"))).unwrap(),
    );
    assert!(Derived::new(&view).is_ok());
}

/// One part whose record is its `spec` and whose ratings are `key` of data/record.json, with `facts` stated on it.
fn recorded(key: &str, spec: &str, facts: &str) -> String {
    format!(
        "(def spec (read-json \"record.json\"))\n(body part {BALL} (attr ratings (get spec \"{key}\")) {spec} (frame shaft (at)))\n(assembly a (instance p part\n {facts}))\n(require ok (within-ratings a))"
    )
}

/// A limit its record holds is said once: the rating reads the number, the curve's points or one coordinate of
/// a point from the part's record (its `spec`), where a copy beside it could differ. A curve says whether it may
/// be read below its first point.
#[test]
fn a_rating_reads_its_numbers_from_the_parts_record() {
    let spec = "(attr spec spec)";
    let held = |key: &str, facts: &str| asked(&recorded(key, spec, facts), &[]).0;
    // 100 N from the record, either side of it.
    let side = |load: &str| held("ratings", &format!("(fact radial-load shaft {load}) (fact torque shaft 0.01N*m) (fact speed shaft 100rpm)"));
    assert_eq!([side("99N"), side("101N")], [Verdict::Yes, Verdict::No]);
    // The curve from the record: 0.4 N m at 60 rpm, 0.1 at 300. At rest it is read as its first point's 0.4 N m,
    // an estimate that says what it took; without leave to, it is not available, and says why and what to do.
    let at_rest = |key: &str, torque: &str| {
        let mut view = load(&recorded(key, spec, &format!("(fact torque shaft {torque}) (fact speed shaft 0rpm) (fact radial-load shaft 0N)"))).unwrap();
        let (design, formulas) = design(&mut view).unwrap();
        let tools: &[&dyn Tool] = &[&formulas];
        let mut claims = Claims::new(&design, tools);
        (claims.answer(&about("ok", "")), claims.answer(&about("pull-out's limit", "a/p /1")))
    };
    let (ok, limit) = at_rest("ratings", "0.39N*m");
    assert_eq!((ok.held, limit.value, limit.strength), (Some(Verdict::Yes), Some(400.0), Some(Strength::Estimated)), "{limit:?}");
    assert_eq!(limit.left_out, ["the table below its first point, 360 deg/s: taken as its value there"]);
    assert_eq!(at_rest("ratings", "0.41N*m").0.held, Some(Verdict::No));
    let (ok, limit) = at_rest("uncovered", "0.39N*m");
    assert_eq!((ok.held, limit.held), (Some(Verdict::NotAvailable(Unavailable::Outside)), Some(Verdict::NotAvailable(Unavailable::Outside))));
    assert_eq!(
        limit.notes,
        ["its table begins at 360 deg/s, and 'speed' is as low as 0 deg/s: its chart's first torque is above its holding torque. Ask over a span that stays above 60 rpm"]
    );
    // A line whose two ends are record fields: 0.6 N m at rest to none at 600 rpm, so 0.3 N m at 300 rpm.
    let on_line = |torque: &str| held("line", &format!("(fact torque shaft {torque}) (fact speed shaft 300rpm)"));
    assert_eq!([on_line("0.29N*m"), on_line("0.31N*m")], [Verdict::Yes, Verdict::No]);
    // A record reference with no record, or to a key the record lacks, is an error where the design is read.
    let error = |key: &str, spec: &str| design(&mut load(&recorded(key, spec, "(fact radial-load shaft 1N)")).unwrap()).err().unwrap();
    assert!(error("ratings", "").contains("it reads its record, and its body carries none: (attr spec RECORD)"), "{}", error("ratings", ""));
    assert!(error("astray", spec).contains("its record has no"), "{}", error("astray", spec));
}

/// A fact's name is not a free word: a place takes what its ratings and formulas read. A fact stated where
/// nothing takes its name is found, with what that place does take, before anything is asked.
#[test]
fn a_fact_nothing_reads_is_found_with_what_its_place_takes() {
    let facts = "(fact radial-load shaft 50N) (fact sideload shaft 10N) (fact torque shaft 0.1N*m) (fact speed shaft 100rpm) (fact weight 3N)";
    let mut view = load(&recorded("ratings", "(attr spec spec)", facts)).unwrap();
    let (design, _) = design(&mut view).unwrap();
    let takes = ink_ratings::takes(&mut view.doc, "a").unwrap();
    assert_eq!(takes, [about("radial-load", "a/p shaft"), about("torque", "a/p shaft"), about("speed", "a/p shaft")]);
    let unread = ink_ratings::unread(&design, &takes);
    let found: Vec<(&str, &str, &[String])> = unread.iter().map(|(fact, names)| (fact.about.what.as_str(), fact.about.on.as_str(), names.as_slice())).collect();
    let names = ["radial-load".to_string(), "torque".to_string(), "speed".to_string()];
    let on_frame = ["radial-load on shaft".to_string(), "torque on shaft".to_string(), "speed on shaft".to_string()];
    // On the frame: the names that frame takes. On the part itself, which takes none: those its frames take.
    assert_eq!(found, [("sideload", "a/p shaft", &names[..]), ("weight", "a/p", &on_frame[..])]);
    let said = ink_ratings::unread_said(unread[0].0, unread[0].1.as_slice());
    assert!(
        said.starts_with("the fact 'sideload' stated on a/p shaft (")
            && said.ends_with("is read by nothing: no rating or formula there takes that name; that place takes: radial-load, torque, speed"),
        "{said}"
    );
    // The same without a requirement, as `ink check` reads it.
    let mut view = load(&recorded("ratings", "(attr spec spec)", facts).replace("(require ok (within-ratings a))", "")).unwrap();
    let (stated, taken) = ink_ratings::stated_and_taken(&mut view).unwrap();
    assert_eq!((stated.facts.len(), ink_ratings::unread(&stated, &taken).len(), taken), (5, 2, takes));
}
