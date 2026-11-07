//! Core conformance (docs/CORE.md): every operation's value, derivative, and
//! range rules against the reference implementation's cases
//! (conformance/cases.json, from conformance/reference.py).
//!
//! The suite states the spec; `KNOWN` lists where this implementation does
//! not meet it yet. The test fails on any deviation not listed, and on any
//! listed one that no longer happens (so the list only shrinks).

use ink::eval::*;
use ink::interval::I;
use ink::lang::Value;
use ink::tape::Tape;
use serde_json::Value as J;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

/// The repository: what the tests load (catalog parts, fixtures) is found from it.
const ROOT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

/// Deviations from the spec: "check op". New exceptions must identify
/// a specific regression and repro; an absent exception fails as FIXED.
const KNOWN: &[&str] = &[];

/// Load the grid once, and compile all fields (including symbolic diff)
/// through the public IR. Keep load-time services out of the per-case loop.
fn tapes(cases: &J, expressions: &[String]) -> HashMap<String, Tape> {
    let mut src = cases["prelude"].as_str().unwrap().to_owned();
    for (i, e) in expressions.iter().enumerate() {
        src.push_str(&format!("\n(def conformance_{i} {e})"));
    }
    let doc = ink::load_str(&src, "conformance", Path::new(ROOT)).unwrap();
    let Value::Field(sampled, _) = doc.def("sampled").unwrap() else { panic!("sampled must be a field") };
    let fixture = Tape::new(&doc.g, *sampled);
    let g = &fixture.grids[0];
    let expected = &cases["grid"];
    for a in 0..3 {
        assert_eq!(g.origin()[a], expected["origin"][a].as_f64().unwrap(), "bake fixture origin");
        assert_eq!(g.dims()[a], expected["dims"][a].as_u64().unwrap() as usize, "bake fixture dimensions");
    }
    assert_eq!(g.h(), expected["h"].as_f64().unwrap());
    let samples: Vec<f32> = expected["values"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap() as f32).collect();
    assert_eq!(g.values(), samples, "bake samples must equal the independent polynomial samples");
    expressions
        .iter()
        .enumerate()
        .map(|(i, e)| {
            let Value::Field(f, _) = doc.def(&format!("conformance_{i}")).unwrap() else { panic!("{e} is not a field") };
            (e.clone(), Tape::new(&doc.g, *f))
        })
        .collect()
}

fn num(v: &J) -> f64 {
    match v {
        J::String(s) if s == "inf" => f64::INFINITY,
        J::String(s) if s == "-inf" => f64::NEG_INFINITY,
        _ => v.as_f64().unwrap(),
    }
}

fn nums(v: &J) -> Vec<f64> {
    v.as_array().unwrap().iter().map(num).collect()
}

/// Operations whose values come from a math library: allowed a few ulps.
fn library(op: &str) -> bool {
    matches!(op, "exp" | "ln" | "sin" | "cos" | "tan" | "atan2" | "pow" | "emax")
}

/// Equal under the spec (one zero; no NaN, so NaN never matches), within
/// `ulps` units in the last place.
fn same(got: f64, want: f64, ulps: f64) -> bool {
    if got == want {
        return true;
    }
    if !got.is_finite() || !want.is_finite() {
        return false;
    }
    let scale = want.abs().max(got.abs()).max(f64::MIN_POSITIVE);
    (got - want).abs() <= ulps * f64::EPSILON * scale
}

struct Report(BTreeMap<String, (usize, String)>);

impl Report {
    fn fail(&mut self, check: &str, op: &str, detail: String) {
        let e = self.0.entry(format!("{check} {op}")).or_insert((0, detail));
        e.0 += 1;
    }
}

#[test]
fn every_operation_meets_the_core_spec() {
    let src = std::fs::read_to_string(Path::new(ROOT).join("conformance/cases.json")).unwrap();
    let cases: J = serde_json::from_str(&src).unwrap();
    let mut expressions = BTreeSet::new();
    for c in cases["points"].as_array().unwrap() {
        let e = c["expr"].as_str().unwrap();
        expressions.insert(e.to_owned());
        for axis in ["x", "y", "z", "t", "p"].iter().take(c["d"].as_array().unwrap().len()) {
            expressions.insert(format!("(diff {e} {axis})"));
        }
    }
    for c in cases["ranges"].as_array().unwrap() {
        expressions.insert(c["expr"].as_str().unwrap().to_owned());
    }
    let tapes = tapes(&cases, &expressions.into_iter().collect::<Vec<_>>());
    let tape = |e: &str| &tapes[e];
    let mut r = Report(BTreeMap::new());

    for c in cases["points"].as_array().unwrap() {
        let op = c["op"].as_str().unwrap();
        let context = c["context"].as_array();
        let ctx = Ctx::new(context.map_or(0.0, |v| num(&v[0])), vec![context.map_or(0.0, |v| num(&v[1]))]);
        let t = tape(c["expr"].as_str().unwrap());
        let at = nums(&c["at"]);
        let p = [at[0], at[1], at[2]];
        let want = num(&c["value"]);
        let dwant = nums(&c["d"]);
        let ulps = c["ulps"].as_f64().unwrap_or(if library(op) { 4.0 } else { 0.0 });

        let v = eval_point(t, &ctx, p);
        if !same(v, want, ulps) {
            r.fail("value", op, format!("at {at:?}: got {v}, spec {want}"));
        }
        // The symbolic derivative (`diff`, `grad`): a field like any other.
        let dual = eval_dual(t, &ctx, p);
        let g = eval_grads(t, &ctx, &[p])[0];
        for (k, v) in ["x", "y", "z", "t", "p"].iter().take(dwant.len()).enumerate() {
            let e = c["expr"].as_str().unwrap();
            let got = eval_point(tape(&format!("(diff {e} {v})")), &ctx, p);
            if !same(got, dwant[k], 64.0) {
                r.fail("derivative(symbolic)", op, format!("d/d{v} at {at:?}: got {got}, spec {}", dwant[k]));
            }
            // One derivative in every path, bit for bit.
            if got != dual[1 + k] || (k < 3 && got != g[1 + k]) {
                r.fail("derivative(paths)", op, format!("d/d{v} at {at:?}: symbolic {got}, dual {}, gradient {}", dual[1 + k], if k < 3 { g[1 + k] } else { dual[1 + k] }));
            }
        }
        if !same(dual[0], want, ulps) {
            r.fail("value(dual)", op, format!("at {at:?}: got {}, spec {want}", dual[0]));
        }
        if !same(g[0], want, ulps) {
            r.fail("value(grads)", op, format!("at {at:?}: got {}, spec {want}", g[0]));
        }
    }

    for c in cases["ranges"].as_array().unwrap() {
        let op = c["op"].as_str().unwrap();
        let context = c["context"].as_array();
        let mut ctx = Ctx::new(context.map_or(0.0, |v| num(&v[0])), vec![context.map_or(0.0, |v| num(&v[1]))]);
        // Over a span of t and a range of p (a tolerance): every value of each.
        if let Some(s) = c["context_spans"].as_array() {
            let (t, p) = (nums(&s[0]), nums(&s[1]));
            ctx = Ctx { t_span: Some((t[0], t[1])), p_spans: ParamSpans::new([(0, (p[0], p[1]))]), ..Ctx::new(span_mid((t[0], t[1])), vec![span_mid((p[0], p[1]))]) };
        }
        let t = tape(c["expr"].as_str().unwrap());
        let b: Vec<I> = c["box"]
            .as_array()
            .unwrap()
            .iter()
            .map(|ab| {
                let v = nums(ab);
                I::new(v[0], v[1])
            })
            .collect();
        let b = [b[0], b[1], b[2]];
        let hull = nums(&c["hull"]);
        let got = eval_interval(t, &ctx, b, None);
        if !(got.lo <= hull[0] && got.hi >= hull[1]) {
            r.fail("range", op, format!("box {:?}: got [{}, {}], values reach [{}, {}]", c["box"], got.lo, got.hi, hull[0], hull[1]));
        }
        // First order in t and p (the rates that bound motion).
        let got = eval_interval_rates(t, &ctx, b, &Inputs { t: ctx.t_span.is_some(), ..Inputs::of([t], &ctx) }).0;
        if !(got.lo <= hull[0] && got.hi >= hull[1]) {
            r.fail("range(span)", op, format!("box {:?}: got [{}, {}], values reach [{}, {}]", c["box"], got.lo, got.hi, hull[0], hull[1]));
        }
        let g = eval_interval_grad(t, &ctx, b);
        if !(g.v.lo <= hull[0] && g.v.hi >= hull[1]) {
            r.fail("range(grads)", op, format!("box {:?}: got [{}, {}], values reach [{}, {}]", c["box"], g.v.lo, g.v.hi, hull[0], hull[1]));
        }
        for (k, dh) in c["dhull"].as_array().unwrap().iter().enumerate() {
            let dh = nums(dh);
            if !(g.d[k].lo <= dh[0] && g.d[k].hi >= dh[1]) {
                r.fail(
                    "derivative range",
                    op,
                    format!("d/d{} over box {:?}: got [{}, {}], derivatives reach [{}, {}]", ["x", "y", "z"][k], c["box"], g.d[k].lo, g.d[k].hi, dh[0], dh[1]),
                );
            }
        }
    }

    let found: Vec<&String> = r.0.keys().collect();
    let mut msg = String::new();
    for (k, (n, ex)) in &r.0 {
        let tag = if KNOWN.contains(&k.as_str()) { "known" } else { "NEW" };
        msg += &format!("{tag:5} {k}: {n} cases, e.g. {ex}\n");
    }
    let fixed: Vec<&&str> = KNOWN.iter().filter(|k| !r.0.contains_key(**k)).collect();
    for k in &fixed {
        msg += &format!("FIXED {k}: remove it from KNOWN\n");
    }
    eprintln!("{} point cases, {} range cases\n{msg}", cases["points"].as_array().unwrap().len(), cases["ranges"].as_array().unwrap().len());
    assert!(found.iter().all(|k| KNOWN.contains(&k.as_str())) && fixed.is_empty(), "conformance:\n{msg}");
}
