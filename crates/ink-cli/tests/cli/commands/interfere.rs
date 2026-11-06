//! `ink interfere`: each pair of an assembly's instances, its verdict and
//! how it was obtained, copies asked once, witnesses where the parts lie.

use super::*;
use serde_json::json;

fn pairs(d: &Value) -> &Vec<Value> {
    d["pairs"].as_array().unwrap()
}

fn uncached(rows: &[Value]) -> usize {
    rows.iter().filter(|p| p.get("cached") != Some(&Value::Bool(true))).count()
}

fn summary(rows: &[Value]) -> Vec<(String, String, bool)> {
    rows.iter().map(|p| (p["pair"].to_string(), verdict(p), p["cached"] == true)).collect()
}

#[test]
fn seated_screws_hold_to_the_accuracy_and_threads_are_excused_by_their_zones() {
    let shelf = shelf();
    let d = ink(["interfere", &shelf, "--target", "bracket-set"]);
    let rows = pairs(&d);
    assert!(
        d["overlapping_pairs"] == 0 && !rows.is_empty() && rows.iter().all(|p| verdict(p) == "holds" && p["accuracy"] == 0.1),
        "screw heads seated on the brackets hold to the accuracy: {:?}",
        summary(rows)
    );
    // Twice a claim's default, and the answer says what it lets pass, briefly.
    assert!(d["accuracy"] == 0.1 && d["accuracy_why"].as_str().is_some_and(|w| w.contains("up to 0.2 mm thick passes") && w.len() < 200), "{}", d["accuracy_why"]);
    let d = ink(["interfere", &shelf]);
    let zoned: Vec<&Value> = pairs(&d).iter().filter(|p| p.get("contact_zones").is_some_and(|z| !z.is_null() && z.as_array().is_none_or(|a| !a.is_empty()))).collect();
    assert!(
        d["overlapping_pairs"] == 0 && zoned.len() == 4 && zoned.iter().all(|p| verdict(p) != "fails"),
        "the shelf's screws thread into the board inside their contact zones: {}",
        pretty(&Value::Array(zoned.iter().map(|p| (*p).clone()).collect()))
    );
}

#[test]
fn parts_seated_by_frames_hold_guaranteed() {
    // A part seated by frames is placed exactly: its touching faces are
    // guaranteed to the accuracy, not estimates.
    let d = ink(["interfere", &fixture("frames.ink")]);
    assert!(!pairs(&d).is_empty() && pairs(&d).iter().all(|p| verdict(p) == "holds"), "{}", pretty(&d));
}

const UNIT: &str = "(import \"std\")\n(body plate (shape (box 20 20 2)))\n(body pin (shape (move (cylinder-z 2 10) 3 0 4)))\n\
(body cap (shape (move (box 6 6 2) -6 0 2)))\n(assembly unit (instance plate plate) (instance pin pin) (instance cap cap))\n";
const PLACES: [&str; 5] = ["(move 100 0 0)", "(move 0 250.3 7)", "(rot-z 30) (move -400 10 0)", "(rot-x 90) (move 1e5 -3e4 12.25)", "(rot-y 90) (move 3 -700 1)"];

#[test]
fn copies_of_a_subassembly_ask_its_pairs_once_each_witness_where_its_copy_lies() {
    // Copies of a sub-assembly ask its pairs once:
    // a still pair is asked in the frame of the deepest instance holding it,
    // so five copies moved and turned apart store the kernel answers one
    // copy does, and each copy's witness is moved to where that copy lies.
    let s = Scratch::new("copies");
    let copies = |n: usize| {
        let body: String = PLACES[..n].iter().enumerate().map(|(k, w)| format!(" (instance c{k} unit (at {w}))")).collect();
        s.model(&format!("copies-{n}.ink"), &format!("{UNIT}(assembly root{body})\n"))
    };
    let mut stored = vec![];
    let mut d = Value::Null;
    let five = copies(5);
    for m in [copies(1), five.clone()] {
        let cache = s.path(&format!("cache-{}", stored.len()));
        d = s.ink(["interfere", &m]).env("XDG_CACHE_HOME", &cache).ok();
        stored.push(files_under(&std::path::Path::new(&cache).join("ink/answers")));
    }
    let rows = pairs(&d);
    assert!(
        rows.len() == 10 && d["overlapping_pairs"] == 5 && uncached(rows) == 2,
        "five copies' pairs are two questions: each copy's rows but the first's cached: {:?}",
        summary(rows)
    );
    if !no_cache() {
        assert!(stored[0] == stored[1] && stored[0] > 0, "five copies store the kernel answers one copy does: {stored:?}");
    }
    let mut inside = vec![];
    for p in rows.iter().filter(|p| p.get("witness").is_some()) {
        let deepest = p["pair"].as_array().unwrap().iter().map(|t| value_at(&five, t.as_str().unwrap(), &p["witness"], &[])).fold(f64::MIN, f64::max);
        inside.push((p["pair"][0].as_str().unwrap().split('/').nth(1).unwrap().to_string(), deepest));
    }
    let mut copies_seen: Vec<&str> = inside.iter().map(|(c, _)| c.as_str()).collect();
    copies_seen.sort();
    assert!(
        copies_seen == ["c0", "c1", "c2", "c3", "c4"] && inside.iter().all(|(_, v)| *v < -0.05),
        "each copy's witness lies deeper than the accuracy inside both of its parts, where that copy is: {inside:?}"
    );
}

#[test]
fn a_witness_cube_lies_inside_both_parts_in_a_copy_turned_45_degrees() {
    // A witness's whole cube of half-side the accuracy lies in the overlap
    // in the world, in a copy turned 45 degrees (its corners stick out past
    // the ball of radius the accuracy), on an overlap slab 0.3 thick: at
    // the accuracy (0.1) the cube's width across the slab is 0.283. (At
    // 0.135 thick, as written for a finer accuracy, the pair holds with no
    // witness.)
    let s = Scratch::new("slab");
    let slab = s.model(
        "slab.ink",
        "(import \"std\")\n(body block (shape (box 10 10 10)))\n(body plate (shape (move (box 4 4 4) 6.7 0 0)))\n\
(assembly unit (instance block block) (instance plate plate))\n\
(assembly root (instance c0 unit (at (rot-z 45) (move 100 0 0))) (instance c1 unit (at (rot-z 45) (move 0 100 0))))\n",
    );
    let d = ink(["interfere", &slab]);
    let mut corners = vec![];
    for p in pairs(&d).iter().filter(|p| p.get("witness").is_some()) {
        let h = num(&p["accuracy"]);
        let w: Vec<f64> = p["witness"].as_array().unwrap().iter().map(num).collect();
        let pts: Vec<Value> = (0..8).map(|k| json!([0, 1, 2].map(|a| w[a] + if k >> a & 1 == 0 { -h } else { h }))).collect();
        for t in p["pair"].as_array().unwrap() {
            corners.push(max_value_at(&slab, t.as_str().unwrap(), &pts));
        }
    }
    assert!(corners.len() == 4 && corners.iter().all(|&v| v <= 0.0), "both copies' pairs fail with a witness whose cube is inside both parts: {corners:?}");
}

fn nuts(r: f64) -> String {
    let inst: String = (0..4)
        .map(|k| {
            format!(
                " (instance n{k} nut (at (rot-z {}) (rot-x {}) (move {} {} 8))) (instance s{k} screw (at (move 0 0 0.35) (frame-of n{k})))",
                37 * k,
                11 * k,
                10.3 + 25.1 * k as f64,
                3.7 * k as f64
            )
        })
        .collect();
    format!(
        "(import \"std\")\n(body nut (shape (subtract (box 8 8 4) (cylinder-z 1.6 10))))\n(body screw (shape (min (cylinder-z {r} 12) (move (cylinder-z 2.75 3) 0 0 4))))\n(assembly root{inst})\n"
    )
}

#[test]
fn pairs_placed_alike_are_one_question_wherever_each_lies() {
    // A pair is asked in the frame of one of its parts, the other placed
    // relative to it in numbers rounded to a grid, so a screw seated on
    // each of four nuts placed one by one, turned and moved apart, is one
    // proof, not four; a fat screw's witness (0.8 mm into the nut: a
    // shallow one's cube is not shown inside in a turned frame, so the
    // world is asked) is moved to each pair's world.
    let s = Scratch::new("nuts");
    let d = ink(["interfere", &s.model("nuts.ink", &nuts(1.5))]);
    let rows = pairs(&d);
    assert!(rows.len() == 4 && rows.iter().all(|p| verdict(p) == "holds") && uncached(rows) == 1, "{:?}", summary(rows));
    let fat = s.model("nuts-fat.ink", &nuts(2.4));
    let d = ink(["interfere", &fat]);
    let rows = pairs(&d);
    let inside: Vec<f64> = rows
        .iter()
        .filter(|p| p.get("witness").is_some())
        .map(|p| p["pair"].as_array().unwrap().iter().map(|t| value_at(&fat, t.as_str().unwrap(), &p["witness"], &[])).fold(f64::MIN, f64::max))
        .collect();
    assert!(rows.len() == 4 && uncached(rows) == 1 && inside.len() == 4 && inside.iter().all(|&v| v < -0.05), "{:?} {inside:?}", summary(rows));
}

#[test]
fn pairs_placed_by_a_solved_value_are_measured_not_guaranteed() {
    // Boxes placed by a solved value set pairs apart, measured: u stops at
    // its guess 1 (the root is 2, where m would overlap b), so m's two pairs
    // are counted apart from a and b's guaranteed one.
    let s = Scratch::new("num-boxes");
    let head = "(import \"std\")\n(solve s (unknowns (u 1)) (== (* 1e-13 (- u 2)) 0))\n(body cube (shape (box 1 1 1)))\n";
    let nb = s.model("apart.ink", &format!("{head}(assembly root (instance a cube) (instance b cube (at (move 5 0 0))) (instance m cube (at (move (* 3 u) 0 0))))\n"));
    let d = ink(["interfere", &nb]);
    assert!(pairs(&d).is_empty() && d["pairs_apart_by_bounding_boxes"] == 1 && d["pairs_apart_by_estimated_boxes"] == 2, "{}", pretty(&d));
    let nb = s.model("over.ink", &format!("{head}(assembly root (instance a cube) (instance m cube (at (move (* 0.25 u) 0 0))))\n"));
    let d = ink(["interfere", &nb]);
    assert!(d["overlapping_pairs"] == 1 && d["overlapping_pairs_estimated"] == 1, "an overlap of a measured pair is an estimate: {}", pretty(&d));
}

#[test]
fn a_contact_zone_resting_on_a_solve_over_t_is_solved_at_each_step_and_measured() {
    // A contact zone's inputs count as the pair's: a zone resting on a solve
    // that reads t is solved at each step over a span, and measured.
    let s = Scratch::new("zone-solve");
    let nb = s.model(
        "zone.ink",
        "(import \"std\")\n(solve size (unknowns (r 2)) (== r (- 2 (* 2 t))))\n(body b (shape (box 1 1 1)))\n(assembly root (instance a b) (instance b b) (contact a b (box r r r)))\n",
    );
    let d = ink(["interfere", &nb, "--t-span", "0:1"]);
    let p = &d["pairs"][0];
    assert!(
        verdict(p) == "fails estimate" && num(&p["t"]) > 0.5 && p["motion"]["how"] == "measured",
        "a zone shrinking with a solve over t fails, an estimate, once smaller than the boxes (was guaranteed): {}",
        pretty(p)
    );
    let d = ink(["interfere", &nb, "--t", "0"]);
    let v: Vec<String> = pairs(&d).iter().map(verdict).collect();
    assert!(v == ["holds estimate"] && d["pairs"][0].get("solved_inputs").is_some_and(|x| !x.is_null()), "a pair whose zone is solved is measured: {}", pretty(&d));
}

#[test]
fn a_press_fit_interfere_passes_is_one_its_answer_says_passes_and_a_claim_catches() {
    // A shaft 0.4 over its bore's diameter: an overlap 0.2 thick, 0.1 deep.
    // interfere passes it at its 0.1, and says what passes in a fit's
    // terms (not "an overlap shallower than 0.1 is a claim's to check",
    // which reads as if this one were caught); the claim it points to
    // fails it.
    let s = Scratch::new("press-fit");
    let m = s.model(
        "fit.ink",
        "(import \"std\")\n(param shaft-d 10.4)\n(body bush (shape (subtract (cylinder-z 10 12) (cylinder-z 5 20))))\n\
        (body shaft (shape (cylinder-z (/ shaft-d 2) 12)))\n(assembly fit (instance b bush) (instance s shaft))\n\
        (assert pressed (disjoint (shape-of bush) (shape-of shaft)))\n",
    );
    let d = ink(["interfere", &m]);
    assert!(d["overlapping_pairs"] == 0, "{}", pretty(&d));
    assert!(d["accuracy_why"].as_str().is_some_and(|w| w.contains("0.4 mm over its bore's diameter") && w.contains("(disjoint shaft bore)")), "{}", d["accuracy_why"]);
    assert_eq!(verdict(&ink(["assert", &m])["asserts"][0]), "fails");
}

#[test]
fn an_undecided_pair_is_counted_and_its_hint_names_what_interfere_takes() {
    // A shaft 0.4 over its bore with a collar on the bushing's face: the
    // pair is unknown at 0.1. The summary said overlapping_pairs 0 and
    // counted it nowhere, and its hint offered (accuracy x), which
    // interfere does not take.
    let s = Scratch::new("undecided-pair");
    let m = s.model(
        "fit.ink",
        "(import \"std\")\n(param shaft-d 10.4)\n(body bush (shape (subtract (cylinder-z 10 12) (cylinder-z 5 20))))\n\
        (body shaft (shape (min (cylinder-z (/ shaft-d 2) 12) (move (cylinder-z 8 2) 0 0 12))))\n(assembly fit (instance b bush) (instance s shaft))\n",
    );
    let d = ink(["interfere", &m]);
    let rows: Vec<&Value> = pairs(&d).iter().filter(|p| verdict(p) == "unknown").collect();
    assert!(d["undecided_pairs"] == rows.len() && !rows.is_empty(), "{}", pretty(&d));
    assert!(rows.iter().all(|p| p["hint"].as_str().is_some_and(|h| !h.contains("A smaller (accuracy x)") && h.contains("(disjoint a b) (accuracy x)"))), "{}", pretty(&d));
}
