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
