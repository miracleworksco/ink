//! A claim the call's time limit stops, on grids that split cells to the
//! same size (ink::stop): a binary of its own, as the limit stops every
//! search in its process.

use ink::claims::Status::Unknown;
use ink::lang::view::ModelView;
use ink::lang::{Doc, Form, Measured, Package, Reader, Registry};
use ink::sexp::Sexp;
use std::path::Path;
use std::sync::Arc;

/// A shape's metric reading 6 on cells of 2 and 4 mm, each split to 1 mm
/// where it reads it (as a study splits cells round a hinge), the call's
/// time limit passing while cells of 1 mm are measured: the answer is on
/// cells of 2 mm, checked against 4 mm. The two read the part on the same
/// cells, so they agree whatever its answer; < 8 held, its grid error 0.
/// Grids whose finest cells are not each twice the last's give no grid
/// error, and the claim is unknown.
#[test]
fn a_stopped_claim_on_grids_split_alike_is_not_decided() {
    struct Split;
    impl Package for Split {
        fn name(&self) -> &'static str {
            "split"
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
                finest: (h > 1.0).then_some(1.0),
                warnings: Vec::new(),
                excluded: Vec::new(),
                unseen: None,
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
    static SPLIT: Registry = Registry::new(&[&Split]);
    ink::stop::at(4e9);
    let src = "(def s (- (sqrt (+ (square x) (square y) (square z))) 10))\n(require c (reading s < 8))";
    let d = Doc::load_str_in(&SPLIT, src, "test", Path::new(env!("CARGO_MANIFEST_DIR")), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let ctx = d.ctx(&[], 0.0).unwrap();
    let mut v = ModelView::new(d, ctx);
    let answers = ink::claims::run_asserts(&mut v, None, false, ink::claims::ClaimOptions::default()).unwrap().asserts;
    assert_eq!(answers[0].check.status, Unknown, "{:?}", answers[0].check);
}
