//! A banded metric's claim and the call's time limit (ink::stop): a binary
//! of its own, as the limit stops every search in its process.

use ink::claims::Status::Unknown;
use ink::lang::view::ModelView;
use ink::lang::{Band, Doc, Form, Measured, Package, Reader, Registry};
use ink::sexp::Sexp;
use std::path::Path;
use std::sync::Arc;

/// A shape's metric reading 6 within a band of 0 to 12 on every grid, the
/// call's time limit passing while cells of 1 mm are measured: the answer is
/// on cells of 2 mm, its band and their grid error deciding, so < 8 is
/// unknown. (Taking the band only with a ladder, which a stopped measurement
/// has none of, it held on the value alone.)
#[test]
fn a_stopped_claim_decides_on_its_coarser_band() {
    struct Shaped;
    impl Package for Shaped {
        fn name(&self) -> &'static str {
            "shaped"
        }
        fn forms(&self) -> &'static [&'static str] {
            &[]
        }
        fn metrics(&self) -> &'static [(&'static str, &'static str)] {
            &[("reading", "a reading of a shape")]
        }
        fn shape_metrics(&self) -> &'static [&'static str] {
            &["reading"]
        }
        fn shape_res(&self, _v: &mut ModelView, _m: &str, _f: u32, res: Option<f64>) -> Option<f64> {
            Some(res.unwrap_or(1.0))
        }
        fn read(&self, r: &mut Reader, f: &Sexp, _head: &str, _args: &[Sexp]) -> Result<Arc<dyn Form>, ink::InkError> {
            Err(r.err(f, "no forms"))
        }
        fn shape_metric(&self, _v: &mut ModelView, _m: &str, _f: u32, res: Option<f64>, _accuracy: f64) -> Result<Measured, ink::InkError> {
            let h = res.unwrap_or(1.0);
            if h == 1.0 {
                // the call's limit passes while cells of 1 mm are measured
                ink::stop::at(1.0);
            }
            Ok(Measured {
                value: if h == 1.0 { Err("stopped".into()) } else { Ok(6.0) },
                at: None,
                res: Some(h),
                warnings: Vec::new(),
                excluded: Vec::new(),
                unseen: None,
                proven: false,
                held: None,
                band: Some(Band { low: 0.0, high: 12.0, parts: vec![("placement".into(), 6.0)], basis: "a test's".into() }),
            })
        }
    }
    static SHAPED: Registry = Registry::new(&[&Shaped]);
    ink::stop::at(4e9);
    let src = "(def s (- (sqrt (+ (square x) (square y) (square z))) 10))\n(assert c (reading s < 8))";
    let d = Doc::load_str_in(&SHAPED, src, "test", Path::new(env!("CARGO_MANIFEST_DIR")), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    let answers = ink::claims::run_asserts(&mut v, None, false, ink::claims::ClaimOptions::default()).unwrap().asserts;
    assert_eq!(answers[0].check.status, Unknown, "{:?}", answers[0].check);
}
