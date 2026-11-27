//! Tape evaluators: batched points, intervals, batched spatial gradients, and
//! full forward-mode derivatives (x, y, z, t, every parameter).

mod cases;
mod context;
mod corner;
mod diagnostics;
mod fork;
mod gradient;
mod join;
mod linear;
mod point;
mod promise;
mod range;
mod scope;

#[cfg(test)]
mod batch_tests;

#[cfg(test)]
mod fork_depth_tests;

#[cfg(test)]
mod gradient_tests;

pub use cases::{Moment, eval_interval_cases, moment};
pub use context::{Ctx, Derivations, Derive, Inputs, MAX_INPUTS, ParamSpan, ParamSpans, span_mid};
pub(crate) use context::{fixed, fixed_value};
pub use diagnostics::{Undefined, decided_names, names_on_paths, provenance, undefined_cases};
pub use fork::forked_count;
pub use gradient::{eval_interval_grad, eval_interval_grad_every};
pub use linear::{Region, eval_interval_linear};
pub use point::{L, eval_dual, eval_dual_ranged, eval_grads, eval_grads_at, eval_point, eval_points, eval_points_par};
pub use range::{eval_interval, eval_interval_clamped, eval_interval_in};
pub use scope::{
    Scoped, continuous_in_time, eval_interval_rates, eval_interval_rates_across, eval_interval_rates_jump, eval_interval_scope, eval_interval_scope_blur, eval_interval_scope_every,
};
