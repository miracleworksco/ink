//! A wall's thickness claim and the call's time limit (ink::stop): a
//! binary of its own, as the limit stops every search in its process. The
//! claim measures cells of 4h and 2h before h where the call has a limit,
//! so one the limit stops answers with the finest that finished, saying so;
//! one it does not stop answers as with no limit.

mod common;

use ink::lang::answer::check_json;
use ink::lang::claims::{ClaimOptions, check};
use ink::lang::view::ModelView;
use serde_json::Value as J;
use std::path::Path;

const MODEL: &str = r#"(import "std")
(body shell (shape (subtract (sphere 10) (sphere 8))))
(assert walls (min-thickness (shape-of shell) > 1))
(assert thick (min-thickness (shape-of shell) > 5))
"#;

/// The cells asked (mm).
const H: f64 = 0.5;

/// Each claim's answer as printed, at cells of `res`, the cache in
/// `cache`.
fn answers(dir: &Path, cache: &Path, res: f64) -> Vec<J> {
    // SAFETY: one test in this binary, and no ink work runs meanwhile.
    unsafe { std::env::set_var("XDG_CACHE_HOME", cache) };
    let mut v = ModelView::load_in(ink::registry(), &dir.join("m.ink"), &[], 0.0).unwrap_or_else(|e| panic!("{e:?}"));
    let asserts: Vec<_> = v.doc.asserts.iter().map(|a| (a.claim.clone(), a.bound)).collect();
    asserts
        .into_iter()
        .map(|(c, b)| {
            let c = v.doc.claim_at(&c, b, &v.ctx);
            check_json(&check(&mut v, c, b, ClaimOptions { res: Some(res), ..Default::default() }).unwrap_or_else(|e| panic!("{e:?}")), &[])
        })
        .collect()
}

#[test]
fn a_stopped_thickness_answers_with_the_finest_cells_that_finished() {
    let dir = common::Scratch::new("time-limit");
    std::fs::write(dir.join("m.ink"), MODEL).unwrap();
    let cache = |k: &str| dir.join(k);
    // No limit: as ever. A limit far off measures 4h and 2h first, and
    // answers the same, to the bit.
    let unlimited = answers(&dir, &cache("a"), H);
    ink::stop::at(4e9);
    let limited = answers(&dir, &cache("b"), H);
    assert_eq!(limited, unlimited);
    assert!(unlimited.iter().all(|a| a["res"] == H && a.get("stopped_before_res").is_none()), "{unlimited:#?}");
    // What finished before the limit is cached: here 2h and 4h (and 8h),
    // or 4h alone (and 8h, 16h).
    ink::stop::at(0.0);
    let at_2h = answers(&dir, &cache("c"), 2.0 * H);
    let at_4h = answers(&dir, &cache("d"), 4.0 * H);
    // The limit passed before cells of h finished (none are begun past
    // it): the value is the finest that did, with its grid error from the
    // next coarser, deciding the claim where that cannot move it past its
    // bound; an estimate, as at h.
    ink::stop::at(1.0);
    let stopped = answers(&dir, &cache("c"), H);
    for (a, at) in stopped.iter().zip(&at_2h) {
        assert_eq!((&a["value"], &a["at"], a["res"].as_f64(), a["stopped_before_res"].as_f64()), (&at["value"], &at["at"], Some(2.0 * H), Some(H)), "{a:#}");
        assert!(a["accuracy"].is_null() && a["grid_error"].as_f64().is_some_and(|e| e < 0.5), "{a:#}");
        assert!(a["warnings"][0].as_str().unwrap().starts_with(&format!("the call's time limit passed before cells of {H} mm finished")), "{a:#}");
    }
    assert_eq!((stopped[0]["status"].as_str(), stopped[1]["status"].as_str()), (Some("holds"), Some("fails")));
    // Only the coarsest finished: its value, its grid error unknown, so
    // the claim is.
    let coarsest = answers(&dir, &cache("d"), H);
    for (a, at) in coarsest.iter().zip(&at_4h) {
        assert_eq!((&a["value"], a["res"].as_f64(), a["status"].as_str()), (&at["value"], Some(4.0 * H), Some("unknown")), "{a:#}");
        assert!(a.get("grid_error").is_none() && a["undecided"].as_str().unwrap().contains("the coarsest") && a["hint"].as_str().unwrap().contains("larger timeout_s"), "{a:#}");
    }
}
