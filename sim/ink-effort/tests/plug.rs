//! The effort plug on examples/arm.ink (a motor swinging an aluminium arm)
//! and examples/servo-arm-effort.ink: what it offers, when the solver
//! starts, and that its numbers are `ink effort`'s.

use ink_claims::*;
use ink_effort::Effort;
use ink_lang::view::ModelView;
use std::path::Path;

static REGISTRY: ink_lang::Registry = ink_lang::Registry::new(&[&ink_dynamics::Dynamics]);

const ARM: &str = "root/j/arm";
/// What `ink effort examples/arm.ink` printed for root/j/arm, 2026-11-20:
/// its peak's effort (N mm) and rate (deg/s), over 40 steps.
const PRINTED: (f64, f64, usize) = (-106.90387153606505, 139.27907003073264, 40);

fn plug(example: &str) -> Effort<'static> {
    load(&Path::new(env!("CARGO_MANIFEST_DIR")).join("../examples").join(example))
}

/// The plug on a design of its own, kept as long as the test runs.
fn plugged(view: ModelView, target: &str) -> Result<Effort<'static>, String> {
    Effort::new(Box::leak(Box::new(std::cell::RefCell::new(view))), target)
}

fn load(file: &Path) -> Effort<'static> {
    let view = ModelView::load_in(&REGISTRY, file, &[], 0.0).unwrap_or_else(|e| panic!("{}", e.message));
    plugged(view, "root").unwrap()
}

/// A value read off a run on measured masses (each measured to its band, which the order a measure's
/// search takes moves by parts in a million), to the digits the run printed.
fn printed(got: Option<f64>, want: &str) {
    let decimals = want.split('.').nth(1).map_or(0, str::len);
    assert!(got.is_some_and(|got| format!("{got:.decimals$}") == want), "{got:?}, not {want}");
}

fn near(got: Option<f64>, want: f64) {
    assert!(got.is_some_and(|got| (got - want).abs() <= 1e-9 * want.abs()), "{got:?}, not {want}");
}

#[test]
fn the_motor_gets_its_peak_torque_and_its_story_from_one_run() {
    let (effort, design) = (plug("arm.ink"), Design::default());
    let tools: &[&dyn Tool] = &[&effort];
    let mut claims = Claims::new(&design, tools);

    // Loading the design and plugging the tool in start nothing; nor does a question it does not offer.
    let none = Some(Verdict::NotAvailable(Unavailable::NoSource));
    assert_eq!(claims.answer(&about("peak-torque", "root/j")).held, none);
    assert_eq!(claims.answer(&about("mass", ARM)).held, none);
    assert_eq!(claims.answer(&about("peak-force", ARM)).held, none, "the joint turns: it has no force");
    assert_eq!(effort.runs(), 0);

    // The number is the printed one, as a size, in N mm. The story: an estimate, at the span's two ends and
    // 40 middles, nothing left out, no band.
    let torque = claims.answer(&about("peak-torque", ARM));
    // Its band: either way what the arm's proven mass, centre and inertia can move it (no body is read
    // round here), and upward the 4% a peak between the moments may be read low by.
    let mass = -PRINTED.0 - torque.band.unwrap().0;
    assert!(mass > 0.0 && mass < 0.01 * -PRINTED.0, "{torque:?}");
    near(torque.value, -PRINTED.0);
    assert_eq!((torque.unit.as_str(), torque.held), ("N*mm", None));
    near(torque.band.map(|b| b.1), -PRINTED.0 * (1.0 + ink_dynamics::OFF_LINE) + mass);
    assert!(torque.notes.last().unwrap().starts_with("between the moments read nothing is bounded"), "{:?}", torque.notes);
    assert_eq!(torque.strength, Some(Strength::Estimated));
    assert_eq!(torque.covers, [("t".to_string(), Covers::These(PRINTED.2 + 2))]);
    assert!(torque.left_out.is_empty(), "{:?}", torque.left_out);
    assert!(torque.from == "ink effort on root" && torque.notes[0].starts_with("inverse dynamics of the rigid mechanism"), "{torque:?}");
    assert_eq!(effort.runs(), 1);

    // A second fact of the same motion, and a second round of asking, read the run already made.
    let speed = claims.answer(&about("speed-at-peak-torque", ARM));
    near(speed.value, PRINTED.1);
    assert_eq!((speed.unit.as_str(), speed.strength, &speed.covers), ("deg/s", torque.strength, &torque.covers));
    near(Claims::new(&design, tools).answer(&about("peak-torque", ARM)).value, -PRINTED.0);
    assert_eq!(effort.runs(), 1);
}

/// The torque as arm.ink's comment works it out by hand, held level:
/// the arm's weight (N) at its lever (mm), in N mm. A cheap tool, asked
/// before the solver.
struct ByHand;

impl Tool for ByHand {
    fn ask(&self, asked: &About, claims: &mut Claims) -> Option<Fact> {
        (*asked == about("peak-torque", ARM)).then(|| {
            let (weight, lever) = (claims.answer(&about("weight", ARM)), claims.answer(&about("lever", ARM)));
            let mut fact = resting(asked, "arm.ink: m g L / 2, held level", [&weight, &lever].into_iter());
            match (weight.value, lever.value) {
                (Some(w), Some(l)) => {
                    (fact.value, fact.unit) = (Some(w * l), "N*mm".into());
                    fact.band = weight.band.map(|(low, high)| (low * l, high * l));
                }
                _ => (fact.held, fact.strength) = (weight.held.or(lever.held), None),
            }
            fact
        })
    }
}

/// The motor's limit as arm.ink writes it (150 N mm), and the arm's weight
/// and lever as its comment states them.
fn sketch(weight: Option<Fact>) -> Design {
    let lever = Fact::stated(about("lever", ARM), 100.0, "mm", "arm.ink: half the arm's 200 mm");
    let limit = Rule {
        about: about("motor-ok", ""),
        says: Says::Ordered(vec![Term::Read(about("peak-torque", ARM)), Term::Number(150.0, "N*mm".into())]),
        from: "arm.ink: (effort-limits -150 150)".into(),
    };
    Design { facts: weight.into_iter().chain([lever]).collect(), rules: vec![limit] }
}

fn weight() -> Fact {
    Fact::stated(about("weight", ARM), 0.108 * 9.81, "N", "arm.ink: 108 g")
}

/// Whether the motor is enough, with the hand formula (where `by_hand`) asked before the plug, and how often the
/// solver started.
fn motor_ok(design: &Design, by_hand: bool) -> (Fact, usize) {
    let effort = plug("arm.ink");
    let tools: &[&dyn Tool] = if by_hand { &[&ByHand, &effort] } else { &[&effort] };
    let fact = Claims::new(design, tools).answer(&about("motor-ok", ""));
    (fact, effort.runs())
}

#[test]
fn the_solver_starts_only_when_the_formula_leaves_the_torque_open() {
    // The formula settles it (105.9 N mm held level, under 150): the solver never starts. Nothing was proven:
    // the answer rests on the weight and the lever as stated.
    let design = sketch(Some(weight()));
    let (ok, runs) = motor_ok(&design, true);
    assert_eq!((ok.held, ok.strength, runs), (Some(Verdict::Yes), None, 0));
    near(Some(ok.terms[0].1), 0.108 * 9.81 * 100.0);

    // No formula: the solver answers, and the yes carries its story.
    let (ok, runs) = motor_ok(&design, false);
    assert_eq!((ok.held, ok.strength, runs), (Some(Verdict::Yes), Some(Strength::Estimated), 1));
    assert_eq!(ok.covers, [("t".to_string(), Covers::These(PRINTED.2 + 2))]);
    near(Some(ok.terms[0].1), -PRINTED.0);

    // A formula that cannot be worked out (nothing states the weight) leaves it open: the solver answers.
    let (ok, runs) = motor_ok(&sketch(None), true);
    assert_eq!((ok.held, ok.strength, runs), (Some(Verdict::Yes), Some(Strength::Estimated), 1));
}

/// A formula's answer too close to the limit to call (a weight known only
/// between 0.9 and 1.6 N: 90 to 160 N mm against 150) goes on to the
/// solver, which decides it: 106.9 N mm, under 150.
#[test]
fn a_band_across_the_limit_goes_on_to_the_solver() {
    let design = sketch(Some(Fact { band: Some((0.9, 1.6)), ..weight() }));
    let (ok, runs) = motor_ok(&design, true);
    assert_eq!((ok.held, ok.strength, runs), (Some(Verdict::Yes), Some(Strength::Estimated), 1));
    near(Some(ok.terms[0].1), -PRINTED.0);
    // A band clear of the limit is the answer, and the solver is not started.
    let design = sketch(Some(Fact { band: Some((0.9, 1.2)), ..weight() }));
    assert_eq!(motor_ok(&design, true).1, 0);
}

#[test]
fn what_the_solver_read_loosely_is_left_out() {
    let (effort, design) = (plug("servo-arm-effort.ink"), Design::default());
    let tools: &[&dyn Tool] = &[&effort];
    let mut claims = Claims::new(&design, tools);
    // `ink effort examples/servo-arm-effort.ink`, 2026-11-22: 799 steps, four stretches loose; the slow
    // shoulder's peak 310.331 N mm.
    let slow = claims.answer(&about("peak-torque", "root/slow/s/output"));
    printed(slow.value, "310.331");
    // Over the moments the run read, as many as its efforts ask for (801 on 2026-11-22, 783 since).
    assert!(matches!(&slow.covers[..], [(t, Covers::These(_))] if t == "t"), "{:?}", slow.covers);
    // Its own steps were all read by the rule: the four loose stretches the command lists are other actuators',
    // where their ratings may cross their efforts inside a step. Each is left out of that actuator's facts alone.
    assert!(slow.left_out.is_empty(), "{:?}", slow.left_out);
    let mut left = |on: &str| claims.answer(&about("peak-torque", on)).left_out;
    // Said as every answer from moments says a step that is not covered (ink_lang::motion::not_followed).
    let fast = left("root/fast/s/output");
    assert!(fast.len() == 1 && fast[0].starts_with("not settled between the moments read, over t 0.409") && fast[0].contains(ink_lang::motion::MODEL_GAP), "{fast:?}");
    let quick = left("root/quick/s/output");
    assert!(quick.len() == 1 && quick[0].contains("over 3 steps from t "), "{quick:?}");
    // The gate's hinge has no effort limits: it is no actuator, and nothing is offered on it.
    assert_eq!(claims.answer(&about("peak-torque", "root/gate/leaf")).held, Some(Verdict::NotAvailable(Unavailable::NoSource)));
    assert_eq!(effort.runs(), 1);
}

#[test]
fn a_slide_gets_a_force_and_a_turn_beside_it_a_torque_from_the_same_run() {
    let (effort, design) = (plug("slider-loads.ink"), Design::default());
    let tools: &[&dyn Tool] = &[&effort];
    let mut claims = Claims::new(&design, tools);
    // `ink effort examples/slider-loads.ink`, 2026-11-20: 86 steps; the belt's peak 3.424 N, the pan's 97.63 N mm.
    let pull = claims.answer(&about("peak-force", "root/s/c"));
    printed(pull.value, "3.424");
    assert_eq!((pull.unit.as_str(), &pull.covers), ("N", &vec![("t".to_string(), Covers::These(88))]));
    assert_eq!(claims.answer(&about("speed-at-peak-force", "root/s/c")).unit, "mm/s");
    printed(claims.answer(&about("peak-torque", "root/s/pan")).value, "97.63");
    assert_eq!(claims.answer(&about("peak-torque", "root/s/c")).held, Some(Verdict::NotAvailable(Unavailable::NoSource)));
    assert_eq!(effort.runs(), 1);
}

/// arm.ink with one line taken out.
fn arm_without(line: &str) -> Effort<'static> {
    let arm = include_str!("../../examples/arm.ink");
    assert!(arm.contains(line));
    let file = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("arm-without-{}.ink", &line[1..6]));
    std::fs::write(&file, arm.replace(line, "")).unwrap();
    load(&file)
}

#[test]
fn a_design_with_no_span_is_read_at_one_moment_and_one_with_no_gravity_lacks_it() {
    // At t = 0 the arm is level and at its steadiest: m g L / 2 = 105.9 N mm, by hand in arm.ink.
    let (effort, design) = (arm_without("(t-span 0 2s)"), Design::default());
    let level = Claims::new(&design, &[&effort]).answer(&about("peak-torque", ARM));
    assert!(level.value.is_some_and(|v| (v - 105.9).abs() < 0.1), "{level:?}");
    // No motion was asked: there is no between, and its band is what its proven masses allow either way.
    assert_eq!((&level.covers, level.strength), (&vec![("t".to_string(), Covers::One)], Some(Strength::Estimated)));
    assert!(level.band.is_some_and(|(lo, hi)| lo < 105.948 && 105.948 < hi && hi - lo < 0.02 * 105.948), "{level:?}");

    let effort = arm_without("(gravity 0 0 -9.81m/s2)");
    let lacks = Claims::new(&design, &[&effort]).answer(&about("peak-torque", ARM));
    assert_eq!((lacks.held, lacks.value, lacks.strength), (Some(Verdict::NotSettled(Unsettled::Lacks)), None, None));
    assert!(lacks.notes[0].contains("gravity"), "{lacks:?}");
}

/// A part's rating, as a catalog part carries it, is answered by the solver: tests/data/rated-arm.ink is arm.ink
/// with the arm carrying three ratings in units of their own (0.15 N m, 100 N mm, 30 rpm) and the design asking
/// whether it is within them. Nothing in the design states a torque or a speed.
#[test]
fn a_parts_torque_rating_is_answered_by_the_solver() {
    let file = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/rated-arm.ink");
    let mut view = ModelView::load_in(&REGISTRY, &file, &[], 0.0).unwrap_or_else(|e| panic!("{}", e.message));
    let (design, tables) = ink_ratings::design(&mut view).unwrap();
    let effort = plugged(view, "root").unwrap();
    assert_eq!(effort.runs(), 0, "reading the design starts nothing");
    let tools: &[&dyn Tool] = &[&tables, &effort];
    let mut claims = Claims::new(&design, tools);
    // 106.9 N mm against the 0.15 N m its record states (150 N mm): yes, as the solver knows it.
    let torque = claims.answer(&about("motor-torque", ARM));
    assert_eq!((torque.held, torque.strength), (Some(Verdict::Yes), Some(Strength::Estimated)));
    near(Some(torque.terms[0].1), -PRINTED.0);
    assert_eq!(torque.covers, [("t".to_string(), Covers::These(PRINTED.2 + 2))]);
    // Against a motor of 100 N mm: no. Its speed at that peak, 139.3 deg/s, under 30 rpm (180 deg/s): yes.
    assert_eq!(claims.answer(&about("weak-motor", ARM)).held, Some(Verdict::No));
    assert_eq!(claims.answer(&about("motor-speed", ARM)).held, Some(Verdict::Yes));
    // The requirement is no, decided by the rating that fails; and all of it came from one run.
    let all = claims.answer(&about("motor-ok", ""));
    assert_eq!((all.held, all.strength, effort.runs()), (Some(Verdict::No), Some(Strength::Estimated), 1));
}

/// rated-arm.ink with its motor's rating at `limit` N mm: the rating's answer.
fn against(limit: f64) -> Fact {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (ink, json) = (include_str!("data/rated-arm.ink"), include_str!("data/rated-arm.json"));
    assert!(json.contains(r#"{"number": 0.15, "unit": "N*m"}"#) && ink.contains("rated-arm.json"));
    let name = format!("rated-{limit}.json");
    std::fs::write(dir.join(&name), json.replace(r#"{"number": 0.15, "unit": "N*m"}"#, &format!(r#"{{"number": {limit}, "unit": "N*mm"}}"#))).unwrap();
    let file = dir.join(format!("arm-{limit}.ink"));
    std::fs::write(&file, ink.replace("rated-arm.json", &name)).unwrap();
    let mut view = ModelView::load_in(&REGISTRY, &file, &[], 0.0).unwrap_or_else(|e| panic!("{}", e.message));
    let (design, tables) = ink_ratings::design(&mut view).unwrap();
    let effort = plugged(view, "root").unwrap();
    let tools: &[&dyn Tool] = &[&tables, &effort];
    Claims::new(&design, tools).answer(&about("motor-torque", ARM))
}

/// The arm's peak is 106.9 N mm at the 42 moments read, and between them may be up to 4% more, 111.2: a motor
/// rated clear of that is enough, one rated under the peak read is not, and one rated between is too close to
/// call. `ink require` calls the last enough, on the moments alone.
#[test]
fn a_rating_within_what_a_peak_between_moments_may_add_is_not_settled() {
    let straddles = Some(Verdict::NotSettled(Unsettled::Straddles));
    for (limit, held) in [(100.0, Some(Verdict::No)), (106.0, Some(Verdict::No)), (107.0, straddles), (111.0, straddles), (112.0, Some(Verdict::Yes)), (150.0, Some(Verdict::Yes))]
    {
        let rating = against(limit);
        assert_eq!(rating.held, held, "against {limit} N mm");
        assert_eq!(rating.strength, held.filter(|v| !v.open()).map(|_| Strength::Estimated), "against {limit} N mm");
    }
}

/// Each moment the run read is a tool of its own: the arm's torque and speed there, with their sign, on the
/// joint's path and on a frame of it. The largest of them is the peak.
#[test]
fn each_moment_read_gives_its_torque_and_speed() {
    // tests/data/frame-torque.ink: arm.ink, its arm with a frame.
    let (effort, design) = (load(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/frame-torque.ink")), Design::default());
    // Before the run the moments are the written motion's own (the span's ends and 64 even middles), and
    // asking for a moment's tool starts nothing.
    assert!(effort.at(65).is_some() && effort.at(66).is_none() && effort.runs() == 0);
    let at = |moment: usize, what: &str, on: &str| {
        let tool = effort.at(moment).unwrap();
        Claims::new(&design, &[&tool]).answer(&about(what, on))
    };
    // t = 0: the arm is level and at rest, turning up at 45 pi deg/s: m g L / 2 = 105.9 N mm, by hand in arm.ink.
    // On a frame of it across the joint's axis (`other-axis`, its z the world's, the joint turning about y),
    // the moment about that z: the arm's weight and inertia act about y alone, so none.
    let across = at(0, "torque", &format!("{ARM} other-axis"));
    assert!(across.value.is_some_and(|v| v.abs() < 1e-6) && across.notes.iter().any(|n| n.contains("about its frame other-axis's z")), "{across:?}");
    let first = at(0, "torque", ARM);
    assert!(first.value.is_some_and(|v| (v + 105.9).abs() < 0.1), "{first:?}");
    assert_eq!((first.unit.as_str(), first.strength, &first.covers), ("N*mm", Some(Strength::Estimated), &vec![("t".to_string(), Covers::One)]));
    assert_eq!((first.notes[0].as_str(), effort.runs()), ("at t = 0", 1));
    // Either way what its proven masses allow, and away from nothing what a peak between moments may add, 4%
    // of the motion's peak.
    let band = first.band.unwrap();
    let mass = band.1 - first.value.unwrap();
    assert!(mass > 0.0 && mass < 0.01 * -PRINTED.0, "{first:?}");
    near(Some(band.0), first.value.unwrap() - mass + ink_dynamics::OFF_LINE * PRINTED.0);
    // Its speed, and as its band every speed from there to the next moment's.
    let speed = at(0, "speed", ARM);
    near(speed.value, 45.0 * std::f64::consts::PI);
    near(speed.band.map(|b| b.1), 45.0 * std::f64::consts::PI);
    assert!(speed.unit == "deg/s" && speed.band.is_some_and(|b| b.0 < b.1 && b.0 == at(1, "speed", ARM).value.unwrap()), "{speed:?}");
    // What it does not offer: a force on a joint that turns, anything on another part, a peak (the run's own),
    // a torque on a frame the arm does not have.
    let none = Some(Verdict::NotAvailable(Unavailable::NoSource));
    assert_eq!((at(0, "force", ARM).held, at(0, "torque", "root/j").held, at(0, "peak-torque", ARM).held), (none, none, none));
    assert_eq!((at(0, "torque", &format!("{ARM} missing")).held, at(0, "torque", &format!("{ARM} shaft")).held), (none, none));
    // The run read 42 moments; the largest torque among them is the peak `ink effort` prints, at its moment.
    // The swing is symmetric: as large at its mirror moment, 0.9452; of peaks equal within their rounding the
    // earlier is the one shown.
    assert!(effort.at(PRINTED.2 + 1).is_some() && effort.at(PRINTED.2 + 2).is_none());
    let size = |a: &Fact| a.value.unwrap().abs();
    let most = (0..PRINTED.2 + 2).map(|m| at(m, "torque", ARM)).reduce(|best, a| if size(&a) > size(&best) * (1.0 + 1e-9) { a } else { best }).unwrap();
    near(most.value, PRINTED.0);
    assert_eq!((most.notes[0].as_str(), effort.runs()), ("at t = 0.05483579498241968", 1));
}

/// A speed is the written motion's: it needs no mass and no gravity, and starts no run. An effort is given
/// where the masses it rests on are there, and refused, naming them, only where one is not
/// (tests/data/massless.ink: arm.ink's motor beside one whose arm has no mass and one carrying such a body).
#[test]
fn a_speed_needs_no_mass_and_a_body_without_one_blocks_only_the_efforts_it_would_change() {
    let design = Design::default();
    let lacks = Some(Verdict::NotSettled(Unsettled::Lacks));
    let effort = load(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/massless.ink"));
    let at = |moment: usize, what: &str, on: &str| {
        let tool = effort.at(moment).unwrap();
        Claims::new(&design, &[&tool]).answer(&about(what, on))
    };
    // Every motor's speed at every moment of the motion, by hand (the rate of A sin(pi t) is A pi cos(pi t)),
    // with no run: the middle of the last step is t = 2 - 1/64.
    for (on, turn) in [(ARM, 45.0), ("root/k/arm", 30.0), ("root/c/arm", 20.0)] {
        near(at(0, "speed", on).value, turn * std::f64::consts::PI);
        let last = at(64, "speed", on);
        assert_eq!(last.notes[0], "at t = 1.984375");
        assert!(last.value.is_some_and(|v| (v - turn * std::f64::consts::PI * (std::f64::consts::PI * 1.984375).cos()).abs() < 1e-9), "{last:?}");
        assert_eq!((last.unit.as_str(), last.strength, last.held), ("deg/s", Some(Strength::Estimated), None));
        assert!(last.band.is_some_and(|b| b.0 < b.1 && b.1 == at(65, "speed", on).value.unwrap()), "{last:?}");
    }
    assert!(effort.at(66).is_none() && effort.runs() == 0);
    // A moment of the written motion has no effort of the run's moment of its number, past the first.
    let other = at(3, "torque", ARM);
    assert_eq!((other.held, other.value), (Some(Verdict::NotSettled(Unsettled::Undecided)), None), "{other:?}");
    assert!(other.notes[0].contains("read from the written motion alone"), "{other:?}");
    // The arm with its mass: the torque `ink effort` prints for arm.ink, to the digit. The bodies without
    // one are on other joints and change nothing of it.
    near(Claims::new(&design, &[&effort]).answer(&about("peak-torque", ARM)).value, -PRINTED.0);
    // The two a massless body moves on: no torque, at a moment or as a peak, and which body.
    for on in ["root/k/arm", "root/c/arm"] {
        for fact in [at(0, "torque", on), Claims::new(&design, &[&effort]).answer(&about("peak-torque", on))] {
            assert_eq!((fact.held, fact.value, fact.strength), (lacks, None, None), "{fact:?}");
            assert!(fact.notes[0].starts_with("these moving bodies its effort rests on have no mass: bare (") && fact.notes[0].contains(on), "{fact:?}");
        }
        assert!(at(5, "speed", on).value.is_some());
    }

    // No gravity declared: every speed still, as it was; no effort.
    let effort = arm_without("(gravity 0 0 -9.81m/s2)");
    let moment = effort.at(0).unwrap();
    near(Claims::new(&design, &[&moment]).answer(&about("speed", ARM)).value, 45.0 * std::f64::consts::PI);
    assert_eq!(effort.runs(), 0);
    let torque = Claims::new(&design, &[&moment]).answer(&about("torque", ARM));
    assert!(torque.held == lacks && torque.notes[0].contains("gravity"), "{torque:?}");
    assert_eq!(effort.runs(), 0, "nothing is built to say there is no gravity");
}

/// A design loaded as a command loads it, with these options.
fn asked(file: &str, options: &[&str]) -> ModelView {
    let file = Path::new(env!("CARGO_MANIFEST_DIR")).join(file);
    let words: Vec<String> = ["require", file.to_str().unwrap()].into_iter().chain(options.iter().copied()).map(String::from).collect();
    ink_ratings::load(&REGISTRY, &ink_lang::command::Args::parse(&words)).unwrap_or_else(|e| panic!("{}", e.message))
}

/// A run is over the span the design is asked over, and at one value of each param: asked over a range of one
/// its model reads, the plug gives nothing, and says what would.
#[test]
fn a_run_is_over_the_span_asked_and_refuses_a_ranged_param_it_reads() {
    let design = Design::default();
    // arm.ink's own span is 0 to 2; asked over its first tenth, no moment read is past it.
    let effort = plugged(asked("../examples/arm.ink", &["--t-span", "0:0.1"]), "root").unwrap();
    let torque = Claims::new(&design, &[&effort]).answer(&about("peak-torque", ARM));
    let Covers::These(moments) = torque.covers[0].1 else { panic!("{torque:?}") };
    let last = effort.at(moments - 1).unwrap();
    assert_eq!(Claims::new(&design, &[&last]).answer(&about("speed", ARM)).notes[0], "at t = 0.1");
    assert!(moments != PRINTED.2 + 2 && effort.at(moments).is_none(), "{moments}");

    // The arm's length as a param, asked from 150 to 250 mm: its weight and so its torque read it.
    let arm = include_str!("../../examples/arm.ink");
    let long = arm.replace("(body arm (shape (move (box 200 20 10) 100 0 0))", "(param long 200)\n(body arm (shape (move (box long 20 10) (/ long 2) 0 0))");
    assert_ne!(long, arm);
    let file = Path::new(env!("CARGO_TARGET_TMPDIR")).join("arm-long.ink");
    std::fs::write(&file, long).unwrap();
    for (options, open) in [(vec!["--param", "long=150:250"], true), (vec!["--param", "long=250"], false), (vec![], false)] {
        let effort = plugged(asked(file.to_str().unwrap(), &options), "root").unwrap();
        let moment = effort.at(0).unwrap();
        let tools: &[&dyn Tool] = &[&moment, &effort];
        for what in ["peak-torque", "torque", "speed"] {
            let fact = Claims::new(&design, tools).answer(&about(what, ARM));
            assert_eq!(fact.held, open.then_some(Verdict::NotAvailable(Unavailable::Outside)), "{what} {options:?}: {fact:?}");
            assert!(!open || (fact.value.is_none() && fact.notes[0].contains("reads ranged 'long' (150 to 250)")), "{fact:?}");
        }
    }
}

/// tests/data/speed-crossed.ink: where the solver refuses the effort (its limits cross), the speed there is
/// still the motion's, and a rating on the speed alone is answered by it.
#[test]
fn a_refused_effort_does_not_take_the_speed_with_it() {
    let mut view = asked("tests/data/speed-crossed.ink", &[]);
    let (design, tables) = ink_ratings::design(&mut view).unwrap();
    let effort = plugged(view, "root").unwrap();
    let place = "root/j/arm shaft";
    let at = |moment: usize, what: &str| {
        let tool = effort.at(moment).unwrap();
        Claims::new(&design, &[&tool]).answer(&about(what, place))
    };
    // Some moment's torque is refused, and its speed is a number.
    let refused = (0..).map_while(|m| effort.at(m).map(|_| m)).find(|&m| at(m, "torque").open()).expect("a refused moment");
    assert_eq!(at(refused, "torque").held, Some(Verdict::NotSettled(Unsettled::Undecided)));
    assert!(at(refused, "speed").value.is_some_and(|v| v > 1.0), "{:?}", at(refused, "speed"));
    // The rating, a speed of at most 1 deg/s, over the motion: no.
    let tools: &[&dyn Tool] = &[&tables, &effort];
    let moments = |m: usize| effort.at(m).map(|tool| Box::new(tool) as Box<dyn Tool + '_>);
    let (free, _) = over(&design, &about("free-speed", "root/j/arm"), "t", tools, 1, &moments);
    assert_eq!((free.held, free.strength), (Some(Verdict::No), Some(Strength::Estimated)), "{free:?}");
}

/// tests/data/speed-notch.ink: the curve dips to nothing at 2 deg/s, a speed the joint passes between two
/// moments read. A moment's speed reaches to its neighbours', so the curve is read over the dip: not settled,
/// where reading it at the moments alone said yes at all 36.
#[test]
fn a_limit_is_read_over_the_speeds_between_moments() {
    let mut view = asked("tests/data/speed-notch.ink", &[]);
    let (design, tables) = ink_ratings::design(&mut view).unwrap();
    let effort = plugged(view, "root").unwrap();
    let tools: &[&dyn Tool] = &[&tables, &effort];
    let moments = |m: usize| effort.at(m).map(|tool| Box::new(tool) as Box<dyn Tool + '_>);
    let (rating, _) = over(&design, &about("pull-out", "root/j/arm"), "t", tools, 1, &moments);
    assert_eq!(rating.held, Some(Verdict::NotSettled(Unsettled::Straddles)), "{rating:?}");
}

/// The loads the model works out, as facts (examples/slider-loads.ink): the pan axis's reaction is divided among
/// two radial bearings and a thrust bearing, `(reaction root/s/pan (by root/s/ba root/s/bb root/s/th))`, and each
/// is offered its share at each moment; the pan's own part then carries none of it.
#[test]
fn a_supports_share_and_a_joints_load_are_offered_at_each_moment() {
    let (effort, design) = (plug("slider-loads.ink"), Design::default());
    let start = effort.at(0).unwrap();
    let tools: &[&dyn Tool] = &[&start];
    let mut claims = Claims::new(&design, tools);
    let mut fact = |what: &str, on: &str| claims.answer(&about(what, on));
    let (low, high, thrust) = (fact("radial-load", "root/s/ba"), fact("radial-load", "root/s/bb"), fact("axial-load", "root/s/th"));
    for share in [&low, &high, &thrust] {
        assert_eq!((share.unit.as_str(), share.held, share.strength, &share.covers), ("N", None, Some(Strength::Estimated), &vec![("t".to_string(), Covers::One)]));
        assert_eq!(share.notes[0], "at t = 0");
        assert_eq!(share.notes[1], "its share of what root/s/pan carries, by statics: (supports root/s/pan root/s/ba root/s/bb root/s/th)");
        // As good as the masses under it: what each body's measured moments, within their proven bounds, can
        // move it, the same either way (the pan's and the block's, under 1% here).
        let (value, band) = (share.value.unwrap(), share.band.unwrap());
        let reach = value - band.0;
        assert!(reach > 0.0 && reach < 0.01 * value && (band.1 - value - reach).abs() < 1e-9 * value, "{share:?}");
        assert!(share.notes[2].starts_with("within ") && share.notes[2].contains("its proven bounds"), "{share:?}");
        assert_eq!(share.left_out, ["a reaction's peak between the moments read: they are chosen for the efforts"]);
    }
    // `ink effort examples/slider-loads.ink --t 0`, 2026-11-20: the pan's reaction is 1.0959 N across and 15.0174 N
    // along its axis, with 523.64 N mm across it. The thrust bearing takes the 15.0174; the radial two, 20 mm
    // apart, push opposite ways to hold the moment, and differ by the 1.0959 across.
    near(thrust.value, 15.01741509137868);
    near(low.value, 24.647752844211468);
    near(high.value, 25.743605118344);
    assert!((high.value.unwrap() - low.value.unwrap() - 1.0958522741325276).abs() < 1e-9);
    // A radial bearing carries nothing along its axis, the thrust bearing nothing across.
    assert_eq!((fact("axial-load", "root/s/ba").value, fact("radial-load", "root/s/th").value), (Some(0.0), Some(0.0)));
    // The bearings turn as the pan does: 90 sin(pi t) degrees is 90 pi deg/s at the start.
    near(fact("speed", "root/s/bb").value, 90.0 * std::f64::consts::PI);
    // The pan's own part: its reaction is held by the supports the design names.
    let own = fact("radial-load", "root/s/pan");
    assert_eq!(own.value, Some(0.0));
    assert!(own.notes[1].starts_with("what root/s/pan carries is held by its supports, which the design names"), "{:?}", own.notes);
    // What is no load, and a place that is neither a support nor a joint's part, are not offered.
    assert_eq!(fact("tilt", "root/s/ba").held, Some(Verdict::NotAvailable(Unavailable::NoSource)));
    assert_eq!(fact("radial-load", "root/s/camera").held, Some(Verdict::NotAvailable(Unavailable::NoSource)));
    assert_eq!(start.measures(), Some(Worse::More));
}
