//! Dynamics' own stepping constants and rate bound, kept here so it rests
//! on no stepping of the kernel's (the kernel's motion stepping may change
//! or go): a span's first steps, the most one answer may use, and how fast
//! a tape can change over a step.

use ink_core::eval::{Ctx, Inputs, eval_interval_rates, span_mid};
use ink_core::interval::{I, mag};
use ink_core::tape::Tape;

/// Steps a span starts with.
pub const START_STEPS: usize = 16;

/// Steps one answer may use at most.
pub const MAX_STEPS: usize = 4096;

/// The largest |d/dt| of a tape over box `b` and step `s` (every ranged
/// param of ctx at every value): None where the rate is not bounded (the
/// value may jump in the step).
pub fn rate(tape: &Tape, ctx: &Ctx, b: [I; 3], s: (f64, f64)) -> Option<f64> {
    let over = Ctx { t: span_mid(s), t_span: Some(s), ..ctx.clone() };
    let (_, r) = eval_interval_rates(tape, &over, b, &Inputs { t: true, params: Vec::new() });
    let r = r.first().map_or(0.0, |r| mag(*r));
    r.is_finite().then_some(r)
}

/// Where a quantity turns back between `lo` and `hi`, its rate's sign at
/// `lo` `falling`'s (below 0) and at `hi` the other: found by bisection on
/// the sign of `rate`, 60 halvings. A simulated coordinate's extremes
/// (Simulation::moments) and an effort's between its steps (Model::choose)
/// are found so.
pub fn turn_back(mut lo: f64, mut hi: f64, falling: bool, rate: impl Fn(f64) -> f64) -> f64 {
    for _ in 0..60 {
        let mid = span_mid((lo, hi));
        if (rate(mid) < 0.0) == falling {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    span_mid((lo, hi))
}
