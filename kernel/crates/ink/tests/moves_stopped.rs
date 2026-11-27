//! Whether a claim's shapes move (claims::moves) past the call's time
//! limit (ink::stop): a binary of its own, as the limit stops every search
//! in its process.

mod common;

use ink::claims::{ClaimOptions, Status, moves, run_asserts};
use ink::eval::forked_count;

/// Past the limit, a part switched by t is taken as moving without a run
/// of its field, and the search over the motion that follows stops at
/// once: the claim is unknown, said to be stopped, having forked nothing.
#[test]
fn moves_runs_nothing_past_the_time_limit() {
    let dir = common::Scratch::new("moves-stopped");
    ink::stop::at(1.0);
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let (moving, check, n) = pool.install(|| {
        let mut v = common::over_motion(&dir, &common::handed(16));
        let span = v.motion.expect("a motion");
        let claim = v.doc.asserts[0].claim.clone();
        let roots = v.doc.claim_roots(&claim);
        let before = forked_count();
        let moving = moves(&v, &roots, span);
        let a = run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
        (moving, a.asserts[0].check.clone(), forked_count() - before)
    });
    assert!(moving && n == 0, "{n} forks");
    assert_eq!(check.status, Status::Unknown, "{check:?}");
}
