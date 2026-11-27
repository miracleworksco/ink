//! `ink topology` (ink-topology) end to end: what a layout out of
//! iterations says. Its solver's properties are ink-topology's and
//! ink-ext's tests/topology.rs.

use crate::common::*;
use serde_json::Value;

/// Whether any string of a list contains `text`.
fn any_has(v: &Value, text: &str) -> bool {
    v.as_array().is_some_and(|a| a.iter().any(|w| w.as_str().is_some_and(|s| s.contains(text))))
}

#[test]
fn a_topology_out_of_iterations_says_it_had_not_settled_unless_keep_fills_its_budget() {
    // Its iterations are its updates (it said 4 for 3), and it says it had
    // not settled, unless keep fills the budget, which more iterations
    // cannot settle (it gave both reasons).
    let s = Scratch::new("topology-layout");
    let study = |name: &str| {
        format!(
            "(study {name} (shape domain) (elastic (E 1MPa) (nu 0.3)) (fix sym (- x 0.25) (dofs x)) (fix roller (max (- 11.75 x) (- z 0.25)) (dofs z))\n\
             \x20 (load top (force 0 0 -1N) (region (max (- x 0.25) (- 3.75 z)))) (res 1))\n"
        )
    };
    let lay = s.model(
        "layout.ink",
        &format!(
            "(import \"std\")\n(def domain (move (box 12 1 4) 6 0.5 2))\n{}\
             (topology s (study mbb) (volume-fraction 0.5) (iterations 3) (filter 1.5))\n\
             (topology full (study mbb) (volume-fraction 0.3) (iterations 3) (keep (- z 3)))\n",
            study("mbb")
        ),
    );
    let d = s.ink(["topology", &lay, "--name", "s"]).ok();
    // An estimate, by a local search, of the density model, as it says.
    assert!(
        d["accuracy"].is_null()
            && d["estimate"].as_str().is_some_and(|e| e.contains("not a bound"))
            && d["method"].as_str().is_some_and(|m| m.contains("local search") && m.contains("not the extracted part's")),
        "{}",
        pretty(&d)
    );
    let t = &d["results"];
    assert!(
        num(&t["iterations"]) == 3.0 && t["settled"] == false && any_has(&d["warnings"], "had not settled"),
        "a topology out of iterations counts its updates and says it had not settled: {}\n{}",
        pretty(t),
        pretty(&d["warnings"])
    );
    let w = s.ink(["topology", &lay, "--name", "full"]).ok()["warnings"].clone();
    assert!(
        any_has(&w, "nothing is left to optimize") && !any_has(&w, "give it more iterations"),
        "a topology whose keep fills its budget says so, not to give it more iterations: {}",
        pretty(&w)
    );
}

#[test]
fn an_imported_topology_lays_out_its_own_study() {
    // child.layout's study is child.mbb. It kept the name `mbb` and looked
    // it up again in the whole document: with no local mbb it was refused,
    // and beside a local mbb (ten times the load) it laid out that one, a
    // hundredfold compliance. Its compliance is the child's own run's.
    let s = Scratch::new("topology-imported");
    let study = |load: &str| {
        format!(
            "(def domain (move (box 12 1 4) 6 0.5 2))\n(study mbb (shape domain) (elastic (E 1MPa) (nu 0.3)) (fix sym (- x 0.25) (dofs x)) (fix roller (max (- 11.75 x) (- z 0.25)) (dofs z))\n\
             \x20 (load top (force 0 0 {load}) (region (max (- x 0.25) (- 3.75 z)))) (res 1))\n"
        )
    };
    let child = s.model("child.ink", &format!("(import \"std\")\n{}(topology layout (study mbb) (volume-fraction 0.5) (iterations 3) (filter 1.5))\n", study("-1N")));
    let own = s.ink(["topology", &child, "--name", "layout"]).ok();
    let compliance = |d: &Value| d["results"]["compliance_nmm"].as_f64().unwrap_or(f64::NAN);
    for parent in ["(import \"child.ink\" as child)\n".to_string(), format!("(import \"std\")\n(import \"child.ink\" as child)\n{}", study("-10N"))] {
        let p = s.model("parent.ink", &parent);
        let d = s.ink(["topology", &p, "--name", "child.layout"]).ok();
        assert!(compliance(&d) == compliance(&own) && compliance(&own).is_finite(), "{parent}\nown: {}\nimported: {}", pretty(&own), pretty(&d));
    }
}

#[test]
fn a_layout_says_where_its_part_buckles_at_its_loads() {
    // A slender column laid out for stiffness under its compression: the
    // redesign's density model runs no buckling analysis, so the part it
    // draws is solved as a study of its own, its buckling with it. At 800
    // N it buckles near 1.6 N (its factor 0.0020 on 1 mm cells, 0.0019 and
    // 0.0020 on the grids about them; a cantilever's Euler load,
    // pi^2 E I / 4 L^2, is 83 N for the whole 6 x 6 section), and says so
    // in the study's own words; under 0.1 N its factor is past 1, and it
    // says nothing.
    let s = Scratch::new("topology-buckle");
    let column = |force: &str| {
        format!(
            "(import \"std\")\n(def domain (move (box 80 6 6) 40 0 0))\n\
             (study press (shape domain) (elastic (E 2000MPa) (nu 0.3)) (fix base (- x 1)) (load top (force -{force} 0 0) (region (- 79 x))) (res 1))\n\
             (topology column (study press) (volume-fraction 0.4) (iterations 30))\n"
        )
    };
    let heavy = s.ink(["topology", &s.model("heavy.ink", &column("800N")), "--name", "column"]).ok();
    let factor = &heavy["results"]["part_buckling_factor"];
    assert!(num(&factor[0]) > 0.0 && num(&factor[1]) < 0.05 && any_has(&heavy["warnings"], "its part buckles: the loads may reach its linear buckling load"), "{}", pretty(&heavy));
    let light = s.ink(["topology", &s.model("light.ink", &column("0.1N")), "--name", "column"]).ok();
    assert!(num(&light["results"]["part_buckling_factor"][0]) > 1.0 && !any_has(&light["warnings"], "buckle"), "{}", pretty(&light));
}
