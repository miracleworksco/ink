//! The same answer at any thread count, for every parallel solver and query.

use ink::design::shape_derivatives as sensitivity;
use ink::eval::{Ctx, ParamSpans};
use ink::query::{BOXES, Verdict, Way, clearance, prove_empty, sums};
use ink::ranged::{ProofOptions, prove};
use std::sync::Arc;

/// Every parallel search, integral and solve gives the same bits with
/// one thread as with any other count: an answer is a function of the
/// question. Each case runs in its own pools (the cases side by side) and
/// is compared, as text, with its one-thread run. A cached query (prove,
/// clearance) is asked as a trial (`cache::trial`: nothing stored) and must have been computed,
/// not found: a one-thread answer served from the cache made every other
/// count pass without running.
#[test]
fn every_answer_is_the_same_at_any_thread_count() {
    type Run = Box<dyn Fn() -> String + Send + Sync>;
    // A cached query's answer, asked as a trial and computed on this thread.
    fn computed<T>(ask: impl FnOnce() -> T) -> T {
        let before = ink::cache::here().0;
        let answer = ink::cache::trial(ask);
        let (misses, found) = ink::cache::here();
        assert!(misses > before && !found, "served from the cache, not computed: a run outside a trial stored this question");
        answer
    }
    // The case's name, the thread counts compared with one, how many times
    // each, and the run.
    let mut cases: Vec<(&str, Vec<usize>, usize, Run)> = Vec::new();

    // Moments, shape derivatives, and the undecided cell an
    // emptiness and a clearance report (the largest, the first in order of
    // those; a clearance's, the box that ended it). Coincident faces stay
    // undecided, and nonbinary coordinates make their cell volumes
    // sensitive to a change in the floating addition tree; the boxes' cell
    // is on the face they share, x = 3.823. The balls are 0.543 apart: a
    // clearance of 0.53 holds, finer than the cells.
    cases.push((
        "query accumulations",
        vec![2, 4],
        1,
        Box::new(|| {
            let (d, shape) = super::geometry::field("(param r 3.7) (def f (sphere r))", "f");
            let ctx = d.ctx(&[], 0.0).unwrap();
            let shape = Arc::new(shape);
            let (_, touch) = super::geometry::field("(def f (max (move (box 7.4 5.3 2.9) 0.123 0.234 0.355) (move (box 7.4 5.3 2.9) 7.523 0.234 0.355)))", "f");
            let (_, a) = super::geometry::field("(def f (sphere 3.7))", "f");
            let (_, b) = super::geometry::field("(def f (move (sphere 2.9) 7.13 0.234 0.355))", "f");
            let (touch, a, b) = (Arc::new(touch), Arc::new(a), Arc::new(b));
            let mut bits = Vec::new();
            let m = sums(&shape, &[None, None], true, &ctx, [-4.1; 3], [4.2; 3], 0.35);
            bits.extend(m.sums.iter().map(|v| v.to_bits()));
            let (area, gradients) = sensitivity(&shape, &[&|_| 1.0], &ctx, [-4.1; 3], [4.2; 3], 0.35);
            bits.push(area.to_bits());
            bits.extend(gradients.into_iter().flatten().map(f64::to_bits));
            let ctx = Ctx::new(0.0, Vec::new());
            let empty = prove_empty(&touch, &ctx, [-5.1; 3], [12.3; 3], Some(0.37), 0.37, 0.0).verdict;
            let way = Way { grown: a, covered: b, lo: [-4.1; 3], hi: [11.1; 3], swapped: false };
            let clear = computed(|| clearance(&[way], &ctx, 0.53, 0.37, BOXES).0);
            for (verdict, x) in [(empty, Some(3.823)), (clear, None)] {
                let Verdict::Unknown { cells, volume, ended, cell: Some([lo, hi]), .. } = verdict else { panic!("{verdict:?}") };
                assert!(ended == ink::query::Ended::Finest, "exercise a complete search, not the shared proof budget");
                assert!(x.is_none_or(|x| lo[0] <= x && x <= hi[0]), "{lo:?} {hi:?}");
                bits.extend([cells as u64, volume.to_bits()]);
                bits.extend(lo.iter().chain(&hi).map(|v| v.to_bits()));
            }
            format!("{bits:?}")
        }),
    ));

    // A claim left undecided at the edge of its tolerance: a search's
    // passes run whole (none stops part way at a budget shared by threads),
    // so its answer and undecided cells are the same. Two cubes
    // overlapping by 0.1: the overlap is 0.05 deep, the tolerance.
    // (These questions are none an earlier test stored in the cache.)
    cases.push((
        "an undecided proof",
        vec![2, 7, 16],
        1,
        Box::new(|| {
            let (d, t) = super::geometry::field("(def f (max (box 10 10 10) (move (box 10 10 10) 9.9 0.5 0)))", "f");
            let (lo, hi) = d.domain();
            let opts = ProofOptions { res: None, finest: 0.025, depth: 0, tol: 0.05 };
            let p = computed(|| prove(&Arc::new(t), &Ctx::new(0.0, Vec::new()), lo, hi, &opts));
            assert!(matches!(p.verdict, Verdict::Unknown { .. }), "{p:?}");
            format!("{p:?}")
        }),
    ));

    // A witness is the question's: where one ends a search (an overlap, a
    // clearance, ranged params) the first in the search's order is
    // reported. Two plates 0.05 apart (a clearance of 0.2 asked), four
    // posts through the lower one (each anywhere 0..0.5 along x): each
    // search has a witness at four places at once, and keeping the first a
    // thread found would change the witness from run to run.
    cases.push((
        "witnesses",
        vec![2, 3, 8, 16, 32],
        4,
        Box::new(|| {
            let src = "(param p 0)
(def lower (max (- (abs x) 30) (- (abs y) 30) (- (abs z) 1)))
(def upper (max (- (abs x) 30) (- (abs y) 30) (- (abs (- z 2.05)) 1)))
(fn post (cx cy) (max (- (abs (- (- x p) cx)) 1) (- (abs (- y cy)) 1) (- (abs z) 2)))
(def both (max lower (min (post -20 -20) (post 20 -20) (post -20 20) (post 20 20))))";
            let tape = |n: &str| Arc::new(super::geometry::field(src, n).1);
            let (lower, upper, both) = (tape("lower"), tape("upper"), tape("both"));
            let at = Ctx::new(0.0, vec![0.0]);
            let ranged = Ctx { p_spans: ParamSpans::new([(0, (0.0, 0.5))]), ..Ctx::new(0.0, vec![0.25]) };
            let (lo, hi) = ([-32.0; 3], [32.0; 3]);
            let overlap = prove_empty(&both, &at, lo, hi, None, 0.025, 0.05).verdict;
            let ranged = prove_empty(&both, &ranged, lo, hi, None, 0.5, 0.0).verdict;
            let way = Way { grown: lower, covered: upper, lo: [-30.0, -30.0, 1.05], hi: [30.0, 30.0, 3.05], swapped: false };
            let clearance = computed(|| clearance(&[way], &at, 0.2, 0.1, BOXES).0);
            let all = format!("{overlap:?} {ranged:?} {clearance:?}");
            assert!(all.matches("Disproven").count() == 3, "{all}");
            all
        }),
    ));

    std::thread::scope(|scope| {
        let runs: Vec<_> = cases
            .iter()
            .map(|(name, threads, times, run)| {
                (
                    name,
                    scope.spawn(move || {
                        let at = |n| rayon::ThreadPoolBuilder::new().num_threads(n).build().unwrap().install(run);
                        let one = at(1);
                        for &n in threads {
                            for _ in 0..*times {
                                assert!(at(n) == one, "{name}: the answer differs with {n} threads");
                            }
                        }
                    }),
                )
            })
            .collect();
        let failed: Vec<_> = runs.into_iter().filter_map(|(name, r)| r.join().is_err().then_some(*name)).collect();
        assert!(failed.is_empty(), "differs with the thread count: {failed:?}");
    });
}
