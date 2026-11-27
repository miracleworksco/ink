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
            vec![Excluded { possible: false, key: "peak", value: 9.0, at: None, what: "a peak it leaves out".into(), decide: "count the peak".into(), singular: false }]
        } else {
            Vec::new()
        };
        (m == "reading").then(|| Measured {
            value: Ok(self.0),
            at: None,
            res: self.1,
            finest: None,
            warnings: Vec::new(),
            excluded,
            unseen: None,
            unseen_finer: true,
            resolves: 0.0,
            unbounded: None,
            extent: None,
            read_error: None,
            proven: false,
            held: None,
            band: Some(Band { low: self.0 - 1.0, high: self.0 + 1.0, parts: vec![("placement".into(), 1.0)], basis: "a test's".into(), narrows: Vec::new() }),
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
    let src =
        "(param p 3)\n(gauge g (* 2 p))\n(require under (reading g < 8))\n(require near (reading g < 6.5))\n(require over (reading g > 8))\n(require above (reading g > 6.5))";
    assert_eq!(statuses(src, &GAUGES), [Holds, Unknown, Fails, Unknown]);
}

/// A gauge whose band would decide each claim, but which leaves out a peak
/// of 9 that would change each answer: each is unknown, the part left out
/// winning over the band. (Deciding on the band after, < 8 held and > 8
/// failed though the answer said the peak would make them otherwise.)
#[test]
fn a_part_left_out_wins_over_a_band() {
    static GAUGES: Registry = Registry::new(&[&Gauge { grid: true, peak: true }]);
    let src =
        "(param p 3)\n(gauge g (* 2 p))\n(require under (reading g < 8))\n(require near (reading g < 7.5))\n(require over (reading g > 8))\n(require above (reading g > 7.5))";
    assert_eq!(statuses(src, &GAUGES), [Unknown, Unknown, Unknown, Unknown]);
}

/// A fin's tip displacement, p, read on cells of 1 mm too coarse to see
/// the fin (held too stiff, its motion reads low, so it says nothing of
/// the fin either way); or, where `seen`, read on cells of h as p + h - 1.
struct Fin {
    seen: bool,
}
struct FinForm(String, u32, bool);
struct FinReading(f64, Option<f64>);

impl Package for Fin {
    fn name(&self) -> &'static str {
        "fin"
    }
    fn forms(&self) -> &'static [&'static str] {
        &["fin"]
    }
    fn metrics(&self) -> &'static [(&'static str, &'static str)] {
        &[("tip", "a fin's tip displacement")]
    }
    fn read(&self, r: &mut Reader, form: &Sexp, _head: &str, args: &[Sexp]) -> Result<Arc<dyn Form>, ink::InkError> {
        let [name, e] = args else { return Err(r.err(form, "usage: (fin name expr)")) };
        Ok(Arc::new(FinForm(r.name(name)?, r.field(e)?, self.seen)))
    }
}

impl Form for FinForm {
    fn name(&self) -> &str {
        &self.0
    }
    fn metrics(&self) -> Vec<&'static str> {
        vec!["tip"]
    }
    fn run(&self, doc: &mut Doc, ctx: &Ctx, opts: &RunOpts) -> Result<Box<dyn Outcome>, ink::InkError> {
        let h = opts.res.unwrap_or(1.0);
        Ok(Box::new(match self.2 {
            true => FinReading(doc.scalar(self.1, ctx) + h - 1.0, Some(h)),
            false => FinReading(doc.scalar(self.1, ctx), None),
        }))
    }
}

impl Outcome for FinReading {
    fn metric(&self, m: &str) -> Option<Measured> {
        (m == "tip").then(|| Measured {
            value: Ok(self.0),
            at: None,
            res: Some(self.1.unwrap_or(1.0)),
            finest: None,
            warnings: Vec::new(),
            excluded: Vec::new(),
            unseen: self.1.is_none().then(|| "a feature thinner than two of its cells".into()),
            unseen_finer: true,
            resolves: 0.0,
            unbounded: None,
            extent: None,
            read_error: None,
            proven: false,
            held: None,
            band: None,
        })
    }
}

/// A grid that cannot see the part decides no claim, holding or failing:
/// each is unknown, saying why and to refine. (Refusing only a hold, the
/// fin's under-read failed a true `>`.)
#[test]
fn a_grid_that_cannot_see_the_part_decides_neither_way() {
    static FINS: Registry = Registry::new(&[&Fin { seen: false }]);
    let src = "(param p 3e-10)\n(fin f p)\n(require stiff (tip f < 0.0013))\n(require moves (tip f > 0.0013))";
    for (name, status, why, how) in reasons(src, &FINS) {
        assert_eq!(status, Unknown, "{name}");
        assert!(why.contains("cannot see the part") && how.starts_with("refine"), "{name}: {why}; {how}");
    }
}

/// Each claim's name, status, and why it is undecided and how to decide it.
fn reasons(src: &str, registry: &'static Registry) -> Vec<(String, Status, String, String)> {
    let d = Doc::load_str_in(registry, src, "test", Path::new(ROOT), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    let answers = ink::claims::run_asserts(&mut v, None, false, ink::claims::ClaimOptions::default()).unwrap().asserts;
    answers
        .into_iter()
        .map(|a| {
            let (why, how) = match &a.check.evidence {
                ink::claims::Evidence::Metric { undecided: Some((why, how)), .. } => (why.clone(), how.clone()),
                _ => Default::default(),
            };
            (a.name, a.check.status, why, how)
        })
        .collect()
}

/// A tip read 3, 4 and 6 on cells of 1, 2 and 4 mm (grid error 1.25): a
/// claim failing by less than that is unknown, as one holding by less is,
/// and one failing by more fails. (Checking only a hold, `> 3.5` failed.)
#[test]
fn a_fail_within_its_grid_error_is_unknown_as_a_hold_is() {
    static FINS: Registry = Registry::new(&[&Fin { seen: true }]);
    let src = "(param p 3)\n(fin f p)\n(require under (tip f < 3.5))\n(require over (tip f > 3.5))\n(require far (tip f > 10))";
    let got: Vec<_> = reasons(src, &FINS).into_iter().map(|(name, status, why, _)| (name, status, why.contains("grid error"))).collect();
    assert_eq!(got, [("under".into(), Unknown, true), ("over".into(), Unknown, true), ("far".into(), Fails, false)]);
}

/// A peak stress read 10 (no grid), leaving out a support's peak read 10.5,
/// both read up to 4% under the true value and 8% over (`ReadError`).
struct Peak;
struct PeakForm(String);
struct PeakReading;

impl Package for Peak {
    fn name(&self) -> &'static str {
        "peak"
    }
    fn forms(&self) -> &'static [&'static str] {
        &["peak"]
    }
    fn metrics(&self) -> &'static [(&'static str, &'static str)] {
        &[("stress", "a peak stress")]
    }
    fn read(&self, r: &mut Reader, form: &Sexp, _head: &str, args: &[Sexp]) -> Result<Arc<dyn Form>, ink::InkError> {
        let [name] = args else { return Err(r.err(form, "usage: (peak name)")) };
        Ok(Arc::new(PeakForm(r.name(name)?)))
    }
}

impl Form for PeakForm {
    fn name(&self) -> &str {
        &self.0
    }
    fn metrics(&self) -> Vec<&'static str> {
        vec!["stress"]
    }
    fn run(&self, _doc: &mut Doc, _ctx: &Ctx, _opts: &RunOpts) -> Result<Box<dyn Outcome>, ink::InkError> {
        Ok(Box::new(PeakReading))
    }
}

impl Outcome for PeakReading {
    fn metric(&self, m: &str) -> Option<Measured> {
        (m == "stress").then(|| Measured {
            value: Ok(10.0),
            at: None,
            res: None,
            finest: None,
            warnings: Vec::new(),
            excluded: vec![Excluded {
                possible: false,
                key: "support",
                value: 10.5,
                at: None,
                what: "a support's peak".into(),
                decide: "model the support".into(),
                singular: true,
            }],
            unseen: None,
            unseen_finer: true,
            resolves: 0.0,
            unbounded: None,
            extent: None,
            read_error: Some(ink::lang::ReadError { below: Some(0.08), above: Some(0.04), proven: false, how: "a test's" }),
            proven: false,
            held: None,
            band: None,
        })
    }
}

/// What a claim leaves out is read as its value is: counted within the
/// read's error too. 10.5 against `< 10.6` was counted exactly and held;
/// read up to 4% higher, it may cross, so the claim is unknown. Far from
/// the bound, `< 12` holds.
#[test]
fn a_part_left_out_is_counted_within_the_read_s_error() {
    static PEAKS: Registry = Registry::new(&[&Peak]);
    let src = "(peak p)\n(require near (stress p < 10.6))\n(require far (stress p < 12))";
    let got: Vec<_> = reasons(src, &PEAKS).into_iter().map(|(n, s, why, _)| (n, s, why.contains("leaves out a support's peak"))).collect();
    assert_eq!(got, [("near".into(), Unknown, true), ("far".into(), Holds, false)]);
}
