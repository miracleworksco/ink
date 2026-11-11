//! Claims on a metric measured with a band (`Measured::band`): one decision
//! on every path a metric claim takes, a reason not to decide winning.

use ink::claims::Status::{self, *};
use ink::eval::Ctx;
use ink::lang::view::ModelView;
use ink::lang::{Band, Doc, Excluded, Form, Measured, Outcome, Package, Reader, Registry, RunOpts};
use ink::sexp::Sexp;
use std::path::Path;
use std::sync::Arc;

const ROOT: &str = env!("CARGO_MANIFEST_DIR");

/// A gauge reading 2p, measured to lie within 1 of it ("placement"): on
/// cells of h it reads 2p + (h - 1) / 2 where `grid` (else it has no
/// cells), and leaves out a peak of 9 where `peak`.
struct Gauge {
    grid: bool,
    peak: bool,
}
struct GaugeForm(String, u32, bool, bool);
struct Reading(f64, Option<f64>, bool);

impl Package for Gauge {
    fn name(&self) -> &'static str {
        "gauge"
    }
    fn forms(&self) -> &'static [&'static str] {
        &["gauge"]
    }
    fn metrics(&self) -> &'static [(&'static str, &'static str)] {
        &[("reading", "the value a gauge reads")]
    }
    fn read(&self, r: &mut Reader, form: &Sexp, _head: &str, args: &[Sexp]) -> Result<Arc<dyn Form>, ink::InkError> {
        let [name, e] = args else { return Err(r.err(form, "usage: (gauge name expr)")) };
        Ok(Arc::new(GaugeForm(r.name(name)?, r.field(e)?, self.grid, self.peak)))
    }
}

impl Form for GaugeForm {
    fn name(&self) -> &str {
        &self.0
    }
    fn metrics(&self) -> Vec<&'static str> {
        vec!["reading"]
    }
    fn run(&self, doc: &mut Doc, ctx: &Ctx, opts: &RunOpts) -> Result<Box<dyn Outcome>, ink::InkError> {
        let h = self.2.then(|| opts.res.unwrap_or(1.0));
        Ok(Box::new(Reading(doc.scalar(self.1, ctx) + h.map_or(0.0, |h| (h - 1.0) / 2.0), h, self.3)))
    }
}

impl Outcome for Reading {
    fn metric(&self, m: &str) -> Option<Measured> {
        let excluded = if self.2 {
            vec![Excluded { key: "peak", value: 9.0, at: None, what: "a peak it leaves out".into(), decide: "count the peak".into(), singular: false }]
        } else {
            Vec::new()
        };
        (m == "reading").then(|| Measured {
            value: Ok(self.0),
            at: None,
            res: self.1,
            warnings: Vec::new(),
            excluded,
            unseen: None,
            proven: false,
            held: None,
            band: Some(Band { low: self.0 - 1.0, high: self.0 + 1.0, parts: vec![("placement".into(), 1.0)], basis: "a test's".into() }),
        })
    }
}

fn statuses(src: &str, registry: &'static Registry) -> Vec<Status> {
    let d = Doc::load_str_in(registry, src, "test", Path::new(ROOT), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    let answers = ink::claims::run_asserts(&mut v, None, false, ink::claims::ClaimOptions::default()).unwrap().asserts;
    answers.iter().map(|a| a.check.status).collect()
}

/// A gauge with no cells (its value no grid's): its band, 5 to 7, decides
/// alone. (Taking the band only with a grid's ladder, the claims near it
/// read holds and fails on the value 6.)
#[test]
fn a_band_decides_a_metric_with_no_grid() {
    static GAUGES: Registry = Registry::new(&[&Gauge { grid: false, peak: false }]);
    let src = "(param p 3)\n(gauge g (* 2 p))\n(assert under (reading g < 8))\n(assert near (reading g < 6.5))\n(assert over (reading g > 8))\n(assert above (reading g > 6.5))";
    assert_eq!(statuses(src, &GAUGES), [Holds, Unknown, Fails, Unknown]);
}

/// A gauge whose band would decide each claim, but which leaves out a peak
/// of 9 that would change each answer: each is unknown, the part left out
/// winning over the band. (Deciding on the band after, < 8 held and > 8
/// failed though the answer said the peak would make them otherwise.)
#[test]
fn a_part_left_out_wins_over_a_band() {
    static GAUGES: Registry = Registry::new(&[&Gauge { grid: true, peak: true }]);
    let src = "(param p 3)\n(gauge g (* 2 p))\n(assert under (reading g < 8))\n(assert near (reading g < 7.5))\n(assert over (reading g > 8))\n(assert above (reading g > 7.5))";
    assert_eq!(statuses(src, &GAUGES), [Unknown, Unknown, Unknown, Unknown]);
}
