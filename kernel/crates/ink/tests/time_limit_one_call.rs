//! A wall's thickness claim the call's time limit stops, in one call with
//! the cache off: a binary of its own, as the limit stops every search in
//! its process. (time_limit.rs reads what finished from the cache, each
//! grid asked in a call of its own; here the cells of 4h and 2h are measured
//! in the call that is stopped, and nothing is kept between calls.)

mod common;

use ink::lang::answer::check_json;
use ink::lang::claims::{ClaimOptions, check};
use ink::lang::sexp::Sexp;
use ink::lang::view::ModelView;
use ink::lang::{Form, InkError, Measured, Package, Reader, Registry};
use serde_json::Value as J;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};

const MODEL: &str = r#"(import "std")
(body shell (shape (subtract (sphere 10) (sphere 8))))
(require walls (min-thickness (shape-of shell) > 1))
"#;

/// The cells asked (mm).
const H: f64 = 0.5;

/// Whether the limit is to pass as cells of H begin.
static ARMED: AtomicBool = AtomicBool::new(false);

/// The kernel's shape package, but the call's time limit passes at the
/// moment cells of H are begun (where ARMED): a stop by what is being
/// measured, not by a clock, so load cannot move it. The measuring is the
/// kernel's own, on every grid.
struct Stops;
impl Package for Stops {
    fn name(&self) -> &'static str {
        ink::Shape.name()
    }
    fn forms(&self) -> &'static [&'static str] {
        ink::Shape.forms()
    }
    fn metrics(&self) -> &'static [(&'static str, &'static str)] {
        ink::Shape.metrics()
    }
    fn metric_quantity(&self, metric: &str) -> Option<ink::units::Quantity> {
        ink::Shape.metric_quantity(metric)
    }
    fn shape_metrics(&self) -> &'static [&'static str] {
        ink::Shape.shape_metrics()
    }
    fn shape_accuracy(&self, metric: &str) -> bool {
        ink::Shape.shape_accuracy(metric)
    }
    fn plane_metric(&self, metric: &str) -> bool {
        ink::Shape.plane_metric(metric)
    }
    fn metric_order(&self, metric: &str) -> f64 {
        ink::Shape.metric_order(metric)
    }
    fn shape_res(&self, v: &mut ModelView, metric: &str, f: u32, res: Option<f64>) -> Option<f64> {
        ink::Shape.shape_res(v, metric, f, res)
    }
    fn shape_metric(&self, v: &mut ModelView, metric: &str, f: u32, res: Option<f64>, accuracy: f64) -> Result<Measured, InkError> {
        if ARMED.load(Relaxed) && res == Some(H) {
            ink::stop::at(1.0);
        }
        ink::Shape.shape_metric(v, metric, f, res, accuracy)
    }
    fn read(&self, r: &mut Reader, form: &Sexp, head: &str, args: &[Sexp]) -> Result<Arc<dyn Form>, InkError> {
        ink::Shape.read(r, form, head, args)
    }
}
static STOPS: Registry = Registry::new(&[&Stops]);

/// The claim's answer as printed, at cells of `res`.
fn answer(file: &Path, res: f64) -> J {
    let mut v = ModelView::load_in(&STOPS, file, &[], 0.0).unwrap_or_else(|e| panic!("{e:?}"));
    let (c, b) = (v.doc.asserts[0].claim.clone(), v.doc.asserts[0].bound);
    let c = v.doc.claim_at(&c, b, &v.ctx);
    check_json(&check(&mut v, c, b, ClaimOptions { res: Some(res), ..Default::default() }).unwrap_or_else(|e| panic!("{e:?}")), &[])
}

#[test]
fn a_thickness_stopped_in_one_call_answers_with_the_cells_it_finished() {
    // SAFETY: one test in this binary, and no ink work has begun.
    unsafe { std::env::set_var("INK_NO_CACHE", "1") };
    let dir = common::Scratch::new("time-limit-one-call");
    std::fs::write(dir.join("m.ink"), MODEL).unwrap();
    // No limit, on cells of 2h: what the stopped call should read.
    let at_2h = answer(&dir.join("m.ink"), 2.0 * H);
    assert!(at_2h["res"] == 2.0 * H && at_2h["value"].as_f64().is_some(), "{at_2h:#}");
    // A limit far off, passing as cells of h begin: cells of 4h and 2h
    // finished in this call, and nothing else holds them.
    ink::stop::at(4e9);
    ARMED.store(true, Relaxed);
    let a = answer(&dir.join("m.ink"), H);
    assert!(ink::stop::reached());
    assert_eq!((&a["value"], &a["at"], a["res"].as_f64(), a["stopped_before_res"].as_f64()), (&at_2h["value"], &at_2h["at"], Some(2.0 * H), Some(H)), "{a:#}");
    assert!(a["grid_error"].as_f64().is_some_and(|e| e < 0.5), "{a:#}");
    assert!(a["warnings"][0].as_str().unwrap().starts_with(&format!("the call's time limit passed before cells of {H} mm finished")), "{a:#}");
    assert_eq!(a["status"], "holds", "{a:#}");
}
