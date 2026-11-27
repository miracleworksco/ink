//! The study plug on examples/beam.ink (an aluminium bar on two supports,
//! pressed at its middle): what it offers, when the solver starts, and that
//! its numbers are `ink solve`'s and `ink require`'s.

use ink_claims::*;
use ink_lang::command::Args;
use ink_lang::view::ModelView;
use ink_study::Studied;
use std::path::{Path, PathBuf};

static REGISTRY: ink_lang::Registry = ink_lang::Registry::new(&[&ink_sim::Sim]);

const BAR: &str = "root/bar";
/// What `ink solve examples/beam.ink` printed, 2026-11-25:
/// results.max_von_mises_mpa.value and results.max_displacement_mm.value,
/// and beside them max_von_mises_near_supports_mpa.value (55.4: the support margin's peak, which bounds
/// nothing, a rigid hold's edge being singular; 23.05 when the supports were penalties over their cells).
/// (The sag 2.93645047 held by Nitsche's method, 2.93639311 by penalties.)
/// (2.93651682, and the margin's peak 23.35, with each pin holding its edge's mean motion, the edge free to deform.)
const PRINTED: (f64, f64, f64) = (304.8956203236258, 2.936516820757774, 23.354410539743593);
/// What `ink require` printed as `grid_error` for a claim on each of the two (`(max-stress press < 335)`,
/// `(max-displacement press < 4)`), 2026-11-25 (three grids converging faster than their elements can: three
/// times their whole movement, 5.546 and 0.0008 before), from cells of 0.5, 1 and 2 mm; and how far over its read it said
/// the stress may lie, "12.21", as it prints it; under it, a tenth.
/// (The sag's 0.00411227 before its pins held their edges' mean motion.)
const GRID_ERROR: (f64, f64) = (25.938175440407292, 0.004049815909390464);
const OVER: f64 = 12.21;

fn data(file: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(file)
}

fn view(file: &Path) -> ModelView {
    ModelView::load_in(&REGISTRY, file, &[], 0.0).unwrap_or_else(|e| panic!("{}", e.message))
}

/// The plug on a design of its own, kept as long as the test runs.
fn plugged(view: ModelView, study: Option<&str>, target: &str) -> Result<Studied<'static>, String> {
    Studied::new(Box::leak(Box::new(std::cell::RefCell::new(view))), study, target)
}

/// The plug for the only study of `file`, on what its assembly `root` places.
fn plug(file: &Path) -> Studied<'static> {
    plugged(view(file), None, "root").unwrap()
}

/// `text` as a file beside rated-beam.ink's data, which it may read.
fn written(name: &str, text: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    // (Copied whole, then named: copied in place, another test read it half written.)
    let part = dir.join(format!("rated-beam.json.{name}"));
    std::fs::copy(data("rated-beam.json"), &part).unwrap();
    std::fs::rename(&part, dir.join("rated-beam.json")).unwrap();
    let file = dir.join(name);
    std::fs::write(&file, text).unwrap();
    file
}

/// rated-beam.ink with `from` written as `to`.
fn beam_with(name: &str, from: &str, to: &str) -> PathBuf {
    let beam = include_str!("data/rated-beam.ink");
    assert!(beam.contains(from), "{from}");
    written(name, &beam.replace(from, to))
}

fn near(got: Option<f64>, want: f64) {
    assert!(got.is_some_and(|got| (got - want).abs() <= 1e-9 * want.abs()), "{got:?}, not {want}");
}

#[test]
fn the_bar_gets_its_stress_and_its_sag_and_their_story_from_one_solve() {
    // beam.ink itself, with the three lines that place its bar: the command's own example.
    let beam = include_str!("../../examples/beam.ink").to_string() + "\n(assembly root (instance bar bar (at)))\n";
    let (studied, design) = (plug(&written("beam-placed.ink", &beam)), Design::default());
    let tools: &[&dyn Tool] = &[&studied];
    let mut claims = Claims::new(&design, tools);

    // Loading the design and plugging the tool in start nothing; nor does a question it does not offer.
    let none = Some(Verdict::NotAvailable(Unavailable::NoSource));
    assert_eq!(claims.answer(&about("max-stress", "root")).held, none);
    assert_eq!(claims.answer(&about("max-stress", "bar")).held, none);
    assert_eq!(claims.answer(&about("mass", BAR)).held, none);
    assert_eq!(studied.runs(), 0);

    // The number is the printed one, in MPa. The story: an estimate on cells of 0.5 mm, its band the grid error
    // three grids show and what its read may miss beyond it, leaving out the stress beside the two supports.
    let stress = claims.answer(&about("max-stress", BAR));
    near(stress.value, PRINTED.0);
    assert_eq!((stress.unit.as_str(), stress.held, stress.strength), ("MPa", None, Some(Strength::Estimated)));
    let (low, high) = stress.band.unwrap();
    near(Some(low), PRINTED.0 - GRID_ERROR.0 - 0.1 * PRINTED.0);
    assert!((high - (PRINTED.0 + GRID_ERROR.0 + OVER)).abs() < 0.005, "{high}");
    assert!(stress.covers.is_empty(), "{:?}", stress.covers);
    assert_eq!(stress.left_out.len(), 1, "{:?}", stress.left_out);
    assert!(stress.left_out[0].starts_with("the stress within 2 mm of fix 'pin', 'roller'"), "{:?}", stress.left_out);
    assert!(stress.left_out[0].ends_with(": 304.9 MPa with it"), "{:?}", stress.left_out);
    assert_eq!(stress.from, "ink solve --study press");
    assert_eq!(stress.notes[0], "a linear static finite-cell solve on cells of 0.5 mm");
    assert_eq!(stress.notes[1], "its grid error about 25.94 MPa (three grids converging faster than their elements can: 303 on cells of 1 mm, 296.2 on cells of 2 mm)");
    assert!(stress.notes[2].starts_with("a peak stress may lie over its read by 22%"), "{:?}", stress.notes);
    assert_eq!(studied.runs(), 1);

    // A second fact of the same solve, and a second round of asking, read the solve already made. Nothing is
    // left out of a sag, and its band is its grid error and the solver's precision (1.0036e-6 of it: the stiffness
    // kept in void and held material, 1e-6, over its cells' share of the part, and its residual; 1.0044e-6 before the
    // supports' edges read their traction partly in the cell beside, the solve then stopping at a residual 2.3e-9
    // lower; 1.0021e-6 while the smoother's estimate was 15 steps of a power iteration, the solve stopping
    // 1.4e-9 lower).
    let sag = claims.answer(&about("max-displacement", BAR));
    near(sag.value, PRINTED.1);
    assert_eq!((sag.unit.as_str(), sag.strength), ("mm", Some(Strength::Estimated)));
    let precision = 1.0035659019e-06 * PRINTED.1;
    near(sag.band.map(|b| b.0), PRINTED.1 - GRID_ERROR.1 - precision);
    near(sag.band.map(|b| b.1), PRINTED.1 + GRID_ERROR.1 + precision);
    assert!(sag.left_out.is_empty() && sag.notes.len() == 3 && sag.notes[2].starts_with("the solver's precision"), "{sag:?}");
    near(Claims::new(&design, tools).answer(&about("max-stress", BAR)).value, PRINTED.0);
    assert_eq!(studied.runs(), 1);
}

/// A part's ratings, as a catalog part carries them, are answered by the solver: tests/data/rated-beam.ink is
/// beam.ink with the bar carrying two limits in units of their own (0.24 GPa, 0.4 cm) and the design asking
/// whether it is within them. Nothing in the design states a stress or a sag.
#[test]
fn a_parts_ratings_are_answered_by_the_solver() {
    let mut view = view(&data("rated-beam.ink"));
    let (design, tables) = ink_ratings::design(&mut view).unwrap();
    let studied = plugged(view, None, "root").unwrap();
    assert_eq!(studied.runs(), 0, "reading the design starts nothing");
    let tools: &[&dyn Tool] = &[&tables, &studied];
    let mut claims = Claims::new(&design, tools);
    // 2.94 mm against the 0.4 cm it may sag: yes, as the solver knows it.
    let sag = claims.answer(&about("sag", BAR));
    assert_eq!((sag.held, sag.strength), (Some(Verdict::Yes), Some(Strength::Estimated)));
    near(Some(sag.terms[0].1), PRINTED.1);
    // 305 MPa against a yield of 0.24 GPa: no. So the requirement is no, by the rating that fails, with what
    // that left out; and all of it came from one solve.
    let all = claims.answer(&about("strong-and-stiff", ""));
    assert_eq!((all.held, all.strength, studied.runs()), (Some(Verdict::No), Some(Strength::Estimated), 1));
    assert_eq!((all.used.len(), all.left_out.len()), (2, 1), "{all:?}");
    near(Some(claims.answer(&about("yield", BAR)).terms[0].1), PRINTED.0);
}

/// The requirement of `file`, the stress its yield rating read, and how often the solver started. `by_hand`:
/// the hand formula is asked first.
fn within(file: &Path, by_hand: Option<Fact>) -> (Fact, f64, usize) {
    let (fact, stress, runs) = within_or_open(file, by_hand);
    (fact, stress.expect("its yield rating read a stress"), runs)
}

/// `within`, where the yield rating may have read no stress (its sources disagree).
fn within_or_open(file: &Path, by_hand: Option<Fact>) -> (Fact, Option<f64>, usize) {
    let mut view = view(file);
    let (mut design, tables) = ink_ratings::design(&mut view).unwrap();
    let studied = plugged(view, None, "root").unwrap();
    let tools: &[&dyn Tool] = if by_hand.is_some() { &[&tables, &ByHand, &studied] } else { &[&tables, &studied] };
    design.facts.extend(by_hand);
    let mut claims = Claims::new(&design, tools);
    let fact = claims.answer(&about("strong-and-stiff", ""));
    (fact, claims.answer(&about("yield", BAR)).terms.first().map(|term| term.1), studied.runs())
}

/// The bar's stress and sag as beam.ink's comment works them out by hand from the push P (N):
/// (P L / 4) c / I = 0.3 P MPa, and P L^3 / 48 E I, with 0.02 mm of shear at 1000 N, = 0.00293 P mm.
struct ByHand;

impl Tool for ByHand {
    fn ask(&self, asked: &About, claims: &mut Claims) -> Option<Fact> {
        let (per_newton, unit) = [("max-stress", (0.3, "MPa")), ("max-displacement", (0.00293, "mm"))].into_iter().find(|(what, _)| asked.on == BAR && asked.what == *what)?.1;
        let push = claims.answer(&about("push", BAR));
        let mut fact = resting(asked, "beam.ink: by hand", [&push].into_iter());
        match push.value {
            Some(p) => {
                (fact.value, fact.unit) = (Some(p * per_newton), unit.into());
                fact.band = push.band.map(|(low, high)| (low * per_newton, high * per_newton));
            }
            None => (fact.held, fact.strength) = (push.held, None),
        }
        Some(fact)
    }
}

fn push(band: Option<(f64, f64)>) -> Fact {
    Fact { band, ..Fact::stated(about("push", BAR), 1000.0, "N", "beam.ink: about 100 kg") }
}

#[test]
fn the_solver_starts_only_when_what_is_stated_or_worked_by_hand_leaves_it_open() {
    // The design states both numbers on the bar. A typed number never makes an answer better than the study's:
    // the solver is asked (once), works out 305 MPa (269 to 323) where 200 was typed, and the two disagree. The
    // 200 MPa alone would have been a yes under the 240 of its yield.
    let stated = beam_with("beam-stated.ink", "(instance bar bar (at))", "(instance bar bar (at) (fact max-stress 200MPa) (fact max-displacement 3.5mm))");
    let (open, stress, runs) = within_or_open(&stated, None);
    assert_eq!((open.held, stress, runs), (Some(Verdict::NotSettled(Unsettled::Disagree)), None, 1));
    // Typed past the study's most, an allowance, it stands: no, on the number typed.
    let roomy = beam_with("beam-roomy.ink", "(instance bar bar (at))", "(instance bar bar (at) (fact max-stress 340MPa) (fact max-displacement 3.5mm))");
    let (no, stress, runs) = within(&roomy, None);
    assert_eq!((no.held, stress, runs), (Some(Verdict::No), 340.0, 1));

    // The hand formula settles both (300 MPa over 240; 2.93 mm under 4): the solver never starts.
    let rated = data("rated-beam.ink");
    let (no, stress, runs) = within(&rated, Some(push(None)));
    assert_eq!((no.held, no.strength, runs), (Some(Verdict::No), None, 0));
    near(Some(stress), 300.0);

    // A push known only between 700 and 900 N is 210 to 270 MPa by hand, across the 240: too close to call, so
    // the solver decides it, 305 MPa. Its sag by hand, at most 2.64 mm, stays the formula's.
    let (no, stress, runs) = within(&rated, Some(push(Some((700.0, 900.0)))));
    assert_eq!((no.held, no.strength, runs), (Some(Verdict::No), Some(Strength::Estimated), 1));
    near(Some(stress), PRINTED.0);
    // Between 700 and 750 N, 210 to 225 MPa: clear of it, and the solver is not started.
    let (yes, _, runs) = within(&rated, Some(push(Some((700.0, 750.0)))));
    assert_eq!((yes.held, yes.strength, runs), (Some(Verdict::Yes), None, 0));

    // A formula that cannot be worked out (nothing states the push) leaves both open: the solver answers.
    let (no, stress, runs) = within(&rated, Some(Fact::stated(about("pull", BAR), 1.0, "N", "not the push")));
    assert_eq!((no.held, no.strength, runs), (Some(Verdict::No), Some(Strength::Estimated), 1));
    near(Some(stress), PRINTED.0);
}

/// The requirement of rated-beam.ink with its yield rating at `limit` MPa.
fn against(limit: f64) -> (Fact, Fact) {
    let json = include_str!("data/rated-beam.json");
    assert!(json.contains(r#"{"number": 0.24, "unit": "GPa"}"#));
    let name = format!("rated-{limit}.json");
    std::fs::write(Path::new(env!("CARGO_TARGET_TMPDIR")).join(&name), json.replace(r#"{"number": 0.24, "unit": "GPa"}"#, &format!(r#"{{"number": {limit}, "unit": "MPa"}}"#)))
        .unwrap();
    let mut view = view(&beam_with(&format!("beam-{limit}.ink"), "rated-beam.json", &name));
    let (design, tables) = ink_ratings::design(&mut view).unwrap();
    let studied = plugged(view, None, "root").unwrap();
    let tools: &[&dyn Tool] = &[&tables, &studied];
    let mut claims = Claims::new(&design, tools);
    (claims.answer(&about("yield", BAR)), claims.answer(&about("strong-and-stiff", "")))
}

/// A stress is 305 MPa and may lie between 248.5 and 343.0: a limit clear of that is decided, and one inside it
/// is too close to call, as `ink require` leaves `(max-stress press < 310)` and `< 335` unknown. A yes or a no
/// is never given off the number alone.
#[test]
fn a_limit_inside_the_band_is_not_settled() {
    let straddles = Some(Verdict::NotSettled(Unsettled::Straddles));
    for (limit, held) in
        [(240.0, Some(Verdict::No)), (248.0, Some(Verdict::No)), (249.0, straddles), (304.0, straddles), (306.0, straddles), (343.0, straddles), (344.0, Some(Verdict::Yes))]
    {
        let (rating, all) = against(limit);
        assert_eq!((rating.held, all.held), (held, held), "against {limit} MPa");
        assert_eq!(rating.strength, held.filter(|v| !v.open()).map(|_| Strength::Estimated), "against {limit} MPa");
    }
}

#[test]
fn what_the_solver_cannot_stand_behind_is_open_and_has_no_number() {
    let design = Design::default();
    let ask = |file: &Path, what: &str| Claims::new(&design, &[&plug(file)]).answer(&about(what, BAR));
    // A bar 0.3 mm thick on cells of 1 mm: its cells cannot see it.
    let thin = beam_with("beam-thin.ink", "(box 210 10 10)", "(box 210 0.3 10)");
    let unseen = ask(&thin, "max-stress");
    assert_eq!((unseen.held, unseen.value, unseen.strength), (Some(Verdict::NotSettled(Unsettled::Undecided)), None, None), "{unseen:?}");
    assert!(unseen.notes[0].starts_with("its cells could not see the part: a feature thinner than two of its cells"), "{:?}", unseen.notes);
    // Pushed along its length with 100 kN: it may buckle, and a linear solve does not say what it does then.
    let column = beam_with("beam-column.ink", "(force 0 0 -1000N)", "(force -100000N 0 0)");
    let buckles = ask(&column, "max-displacement");
    assert_eq!((buckles.held, buckles.value), (Some(Verdict::NotSettled(Unsettled::Undecided)), None), "{buckles:?}");
    assert!(buckles.notes[0].contains("buckl"), "{:?}", buckles.notes);
}

#[test]
fn a_study_is_offered_only_where_it_is_one_and_not_over_a_ranged_size() {
    // A vibration study has no stress under loads; a study of a shape the target does not place is of no part.
    let modal = beam_with("beam-modal.ink", "(res 0.5))", "(modes 2) (res 0.5))");
    assert!(plugged(view(&modal), None, "root").is_err_and(|why| why.contains("not a static elastic study")));
    let other = beam_with("beam-other.ink", "(shape (shape-of bar))", "(shape (move (box 210 10 10) 100 0 1))");
    assert!(plugged(view(&other), None, "root").is_err_and(|why| why.contains("of no part placed under root")));
    assert!(plugged(view(&data("rated-beam.ink")), Some("lift"), "root").is_err());

    // A second study of the same bar is a load case the first one's facts leave out, and they say so.
    let two = beam_with(
        "beam-two.ink",
        "(assembly root",
        "(study lift (shape (shape-of bar)) (elastic (material aluminum-6061)) (fix pin (edge 0)) (load hand (force 0 0 5N) (region (- (abs (- x 100)) 1))) (res 2))\n(assembly root",
    );
    let (studied, design) = (plugged(view(&two), Some("press"), "root").unwrap(), Design::default());
    let sag = Claims::new(&design, &[&studied]).answer(&about("max-displacement", BAR));
    assert_eq!(sag.left_out, ["the other studies of this shape: 'lift'"]);

    // The bar's height asked over 9 to 11 mm: a solve is of one height, which says nothing of the rest, so
    // nothing is given and nothing is solved.
    let sized = beam_with("beam-sized.ink", "(body bar (shape (move (box 210 10 10) 100 0 0))", "(param tall 10)\n(body bar (shape (move (box 210 10 tall) 100 0 0))");
    let args = ["require", sized.to_str().unwrap(), "--param", "tall=9:11"].map(String::from);
    let studied = plugged(ink_ratings::load(&REGISTRY, &Args::parse(&args)).unwrap_or_else(|e| panic!("{}", e.message)), None, "root").unwrap();
    let sag = Claims::new(&design, &[&studied]).answer(&about("max-displacement", BAR));
    assert_eq!((sag.held, sag.value, studied.runs()), (Some(Verdict::NotAvailable(Unavailable::Outside)), None, 0), "{sag:?}");
    assert!(sag.notes[0].contains("reads ranged 'tall'"), "{sag:?}");
}

/// tests/data/two-loads.ink: the bar under two studies, pressed with 1000 N and with 500 N. They are its load
/// cases, and its stress is the worst of them: the heavier one's 305 MPa, whose band reaches across a yield of
/// 300. The lighter one's 152 MPa does not settle that.
#[test]
fn a_parts_studies_are_its_load_cases_and_the_worst_answers() {
    let mut view = view(&data("two-loads.ink"));
    let (design, tables) = ink_ratings::design(&mut view).unwrap();
    let view = std::cell::RefCell::new(view);
    let studies = ink_study::Studies::new(&view, "root");
    let tools: &[&dyn Tool] = &[&tables, &studies];
    let mut claims = Claims::new(&design, tools);
    let stress = claims.answer(&about("max-stress", BAR));
    near(stress.value, PRINTED.0);
    assert!(stress.from == "ink solve --study press" && stress.notes.last().unwrap() == "the worst of its 2 studies; the others: ink solve --study fallback", "{stress:?}");
    assert!(!stress.left_out.iter().any(|l| l.contains("other studies")), "{:?}", stress.left_out);
    assert_eq!(claims.answer(&about("yield", BAR)).held, Some(Verdict::NotSettled(Unsettled::Straddles)));
    assert_eq!((claims.answer(&about("sag", BAR)).held, studies.runs()), (Some(Verdict::Yes), 2));
}
