//! `ink require`: what a file's requirements answer, in the kernel's own
//! `ink` (no solver handed to it): the claims it requires by the prover,
//! what it states by the core, each as a row and, with `--plain`, in words.

use super::*;
use serde_json::json;

const MODEL: &str = r#"(import "std")
(def w (ranged 4 6))
(body ball (shape (sphere 4)))
(require apart (disjoint (shape-of ball) (move 20 0 0 (sphere 4))) (accuracy 0.01))
(require wide (value w > 3))
(require narrow (value w > 5))
(require early (value t < 1s))
(require plain (value 2 > 1))
"#;

#[test]
fn a_files_requirements_are_answered_as_rows_and_in_words() {
    let s = Scratch::new("require");
    let m = s.model("m.ink", MODEL);
    let d = Ink::new(["require", &m]).any();
    assert!(d["ok"] == true && d["all_yes"] == false && d["counts"] == json!({"yes": 4, "no": 1}) && d["all_proven"] == false, "{d}");
    let rows = d["requirements"].as_array().unwrap();
    assert_eq!(rows.iter().map(|r| r["name"].as_str().unwrap()).collect::<Vec<_>>(), ["apart", "wide", "narrow", "early", "plain"]);
    // A row: the verdict and its story, then the prover's evidence by its own keys.
    let wide = &rows[1];
    let story = json!({"verdict": "yes", "strength": "proven", "covers": {"w": "every value"}, "from": format!("{m}:5"), "notes": ["(value w > 3)"]});
    assert!(story.as_object().unwrap().iter().all(|(k, v)| wide[k] == *v), "{wide}");
    assert!(wide["status"] == "holds" && wide["accuracy"] == 0.0 && wide["bounds_over_scope"] == json!([4.0, 6.0]) && wide["claim"] == "(value w > 3)", "{wide}");
    assert!(rows[2]["verdict"] == "no" && rows[2]["status"] == "fails" && rows[0]["accuracy"] == 0.01, "{d}");
    // A yes at one moment of a motion no span covers is a yes there, said, and the file is not all yes for it.
    assert!(rows[3]["verdict"] == "yes" && rows[3]["covers"] == json!({"t": "one value"}) && d["yes_at_one_value"] == json!(["early"]), "{d}");
    let only = Ink::new(["require", &m, "--name", "early"]).any();
    assert!(only["all_yes"] == false && only["counts"] == json!({"yes": 1}), "{only}");
    // One by name, and its keys.
    let wide = ink(["require", &m, "--name", "wide", "--keys", "verdict,strength"]);
    assert!(wide["all_yes"] == true && wide["requirements"] == json!([{"name": "wide", "verdict": "yes", "strength": "proven"}]), "{wide}");
    let none = ink_err(["require", &m, "--name", "wid"]);
    assert!(none["error"]["message"].as_str().unwrap().starts_with("no requirement named 'wid'"), "{none}");
    // In words: the lines `--plain` prints, and its exit status says one is not yes.
    let said = Ink::new(["require", &m, "--plain", "--name", "narrow"]).run();
    assert_eq!((said.code, said.stdout.as_str()), (Some(3), format!("narrow: NO, proven, at every value of w  ({m}:6)\n  note: (value w > 5)\n  with: accuracy 0\n").as_str()));
    // A file that requires nothing answers nothing, and says nothing of all.
    let empty = ink(["require", &s.model("none.ink", "(def a 1)")]);
    assert!(empty["all_yes"].is_null() && empty["requirements"] == json!([]) && empty["nothing_required"].is_string(), "{empty}");
}

/// A param's declared cases are asked as `ink require` asks them: each requirement at the params' values and at
/// every case, its answer the worst of those, saying which case and at how many values.
#[test]
fn a_requirement_is_asked_at_each_declared_case() {
    let s = Scratch::new("require-cases");
    let m = s.model("c.ink", "(param n 2 (cases 3 7))\n(require small (value n < 5))\n(require pos (value n > 0))\n");
    let d = Ink::new(["require", &m]).any();
    assert!(d["all_yes"] == false && d["cases"] == json!([{"n": 3.0}, {"n": 7.0}]), "{d}");
    let (small, pos) = (&d["requirements"][0], &d["requirements"][1]);
    assert!(small["verdict"] == "no" && small["status"] == "fails" && small["cases"][1]["value"] == 7.0 && small["covers"] == json!({"n": 3}), "{d}");
    assert_eq!(small["notes"][1], "the worst of its cases: at the case {\"n\":7.0}", "{d}");
    assert!(pos["verdict"] == "yes" && pos["covers"] == json!({"n": 3}) && pos["notes"].as_array().unwrap().len() == 1, "{d}");
    // How it came out at each case is in its row.
    assert_eq!(pos["cases"].as_array().unwrap().iter().map(|c| c["status"].as_str().unwrap()).collect::<Vec<_>>(), ["holds", "holds"], "{d}");
    // A case --param sets is not asked.
    let set = ink(["require", &m, "--param", "n=4"]);
    assert!(set["all_yes"] == true && set.get("cases").is_none() && set["requirements"][0].get("covers").is_none(), "{set}");
}

/// A wall box: its lid on a hinge pin with 0.2 mm to spare in each knuckle, swung shut over the motion, and a
/// heater plate on four standoffs, `lift` above them.
const WALL_BOX: &str = r#"(import "std")
(t-span 0 1s)
(param lift 0mm)
(body box
  (shape (subtract (min (subtract (box 90 90 50) (move (box 84 84 50) 0 0 3))
                        (move (cylinder-x 4 10) 30 -40 30) (move (box 10 6 10) 30 -42 25)
                        (move (cylinder-x 4 10) -30 -40 30) (move (box 10 6 10) -30 -42 25))
                   (move (cylinder-x 2.2 100) 0 -40 30))))
(body pin (shape (move (cylinder-x 2 76) 0 -40 30)))
(body standoff (shape (cylinder-z 3 8)))
(body plate (shape (move (box 70 70 2) 0 0 (+ -13 (/ lift 1mm)))))
(body lid (shape (subtract (min (cylinder-x 4 20) (move (box 40 80 3) 0 42 0)) (cylinder-x 2.2 100))))
(body handle (shape (move (box 20 4 4) 0 78 3.5)))
(assembly lid-unit (instance lid lid) (instance handle handle))
(assembly root
  (instance box box)
  (instance pin pin)
  (instance s1 standoff (at (move 25 25 -18)))
  (instance s2 standoff (at (move -25 25 -18)))
  (instance s3 standoff (at (move 25 -25 -18)))
  (instance s4 standoff (at (move -25 -25 -18)))
  (instance plate plate)
  (instance lid lid-unit (at (move 0 -40 30)) (turn x (- 90deg (* 90deg/s t)))))
(require ok (within-ratings root))
"#;

/// Where the parts' placement is checked, the parts that touch nothing (none within 1 mm, no relation) are
/// named, by group: a plate lifted off its standoffs, a lid whose pin is gone, a plate with no standoffs. A
/// moving lid held by its pin through a clearance fit, a gap under 1 mm, and a relation are touching; a part
/// resting where one of its supports was is too. A loose part changes no verdict.
#[test]
fn a_part_that_touches_nothing_is_named_where_the_parts_are_placed() {
    let s = Scratch::new("require-loose");
    let asked = |name: &str, model: &str, args: &[&str]| -> (Value, String) {
        let m = s.model(name, model);
        let d = Ink::new(["require", m.as_str()].into_iter().chain(args.iter().copied())).any();
        let said = Ink::new(["require", m.as_str(), "--plain"].into_iter().chain(args.iter().copied())).run().stdout;
        (d, said.lines().next().unwrap_or_default().to_string())
    };
    let (built, said) = asked("box.ink", WALL_BOX, &[]);
    assert!(built["placement"]["loose"] == json!([]), "{built}");
    assert!(said.starts_with("PLACEMENT: no two parts overlap at the start of the motion. None is proven more than 1 mm from all the rest."), "{said}");
    for (lift, loose) in [("4.5", json!([["root/plate"]])), ("0.5", json!([]))] {
        let (d, said) = asked("box.ink", WALL_BOX, &["--param", &format!("lift={lift}")]);
        assert!(d["placement"]["loose"] == loose && d["all_yes"] == built["all_yes"] && d["counts"] == built["counts"], "{lift}: {d}");
        assert_eq!(said.contains("LOOSE: 1 part touches no other part, none within 1 mm and no relation joining it: root/plate."), lift == "4.5", "{said}");
    }
    let (d, said) = asked("nopin.ink", &WALL_BOX.replace("  (instance pin pin)\n", ""), &[]);
    assert_eq!(d["placement"]["loose"], json!([["root/lid/lid", "root/lid/handle"]]), "{d}");
    assert!(said.starts_with("PLACEMENT: no two parts overlap") && said.contains("root/lid/lid, root/lid/handle (touching only each other)"), "{said}");
    let none = WALL_BOX.lines().filter(|l| !l.contains("standoff (at")).collect::<Vec<_>>().join("\n");
    assert_eq!(asked("none.ink", &none, &[]).0["placement"]["loose"], json!([["root/plate"]]));
    let three = WALL_BOX.replace("  (instance s1 standoff (at (move 25 25 -18)))\n", "");
    assert_eq!(asked("three.ink", &three, &[]).0["placement"]["loose"], json!([]));
    let related = WALL_BOX.replace("  (instance plate plate)\n", "  (instance plate plate)\n  (contact plate s1 (move (box 10 10 10) 25 25 -10))\n");
    assert_eq!(asked("related.ink", &related, &["--param", "lift=4.5"]).0["placement"]["loose"], json!([]));
    // Said with a clash, which is what is not covered.
    let (d, said) = asked("nopin.ink", &WALL_BOX.replace("  (instance pin pin)\n", ""), &["--param", "lift=-1"]);
    assert!(d["placement"]["clashes"].as_array().unwrap().len() == 4 && d["placement"]["loose"] == json!([["root/lid/lid", "root/lid/handle"]]), "{d}");
    assert!(said.starts_with("NOT COVERED: the parts are not placed as built: root/s1 and root/plate overlap in solid") && said.contains(". LOOSE: 2 parts"), "{said}");
}
