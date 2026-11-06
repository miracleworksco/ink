//! A claim over a motion as it was checked (`Motion`): with t ranged in
//! the one search, asked once, or measured at moments (ink-lang's
//! motion), and bounds on how fast a field moves along t (`rate`).

use crate::ranged::Status;
use ink_core::eval::{Ctx, Inputs, eval_interval_rates, span_mid};
use ink_core::interval::I;
use ink_core::tape::Tape;

/// Pieces one range may be cut into at most.
pub const MAX_STEPS: usize = 4096;

/// A range checked piece by piece.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Motion {
    pub status: Status,
    pub span: (f64, f64),
    /// Pieces checked.
    pub steps: usize,
    /// The value of a disproof.
    pub t: Option<f64>,
    /// Pieces left unchecked: how far the bodies move in them not bounded,
    /// and no budget left to halve them.
    pub unchecked: Vec<(f64, f64)>,
    /// Pieces undecided at their middle too.
    pub undecided: Vec<(f64, f64)>,
    /// Pieces that hold at their middle only and were not halved further
    /// (at the precision asked, or the budget spent): undecided.
    pub loose: Vec<(f64, f64)>,
    /// How many of `undecided` the call's time limit (crate::stop) left.
    #[serde(default)]
    pub stopped: usize,
}

impl Motion {
    /// Nothing checked yet over `span`, with this status.
    pub fn none(status: Status, span: (f64, f64)) -> Motion {
        Motion { status, span, steps: 0, t: None, unchecked: Vec::new(), undecided: Vec::new(), loose: Vec::new(), stopped: 0 }
    }
}

/// ctx over span s of t, t at its middle (for rates over it).
fn over(ctx: &Ctx, s: (f64, f64)) -> Ctx {
    Ctx { t: span_mid(s), t_span: Some(s), ..ctx.clone() }
}

/// The largest |d/dt| of a tape over box `b` and step `s` (every ranged
/// param of ctx at every value): None where the rate is not bounded (the
/// value may jump in the step).
pub fn rate(tape: &Tape, ctx: &Ctx, b: [I; 3], s: (f64, f64)) -> Option<f64> {
    let (_, r) = eval_interval_rates(tape, &over(ctx, s), b, &Inputs { t: true, params: Vec::new() });
    let r = r.first().map_or(0.0, |r| ink_core::interval::mag(*r));
    r.is_finite().then_some(r)
}
