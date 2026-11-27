//! Checking a document's claims (its `(require NAME CLAIM)` forms, and those its limit clauses and relations make) on a model view:
//! each claim by the cheapest method that decides it (a witness point, an
//! interval proof refined only where undecided, a measurement). Every
//! answer says whether the claim holds and how far off it could be: an
//! accuracy, checked everywhere to it, or none, an estimate (sampled,
//! simulated or solved numerically) with no guaranteed bound. Results are
//! typed; clients name and present them.

use crate::InkError;
use crate::frame::in_frame;
use crate::motion::Asked;
use crate::view::ModelView;
use crate::{Claim, Cmp, FormMetric, Measured, Placed, ReadError, RunOpts};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_points};
use ink_core::graph::Var;
use ink_core::interval::{I, add_o};
use ink_core::tape::Tape;
use ink_kernel::query::{self, Verdict};
use ink_kernel::ranged::{self, Comparison, Proof, ValueSpan};
use ink_kernel::scope::At;
use std::sync::Arc;

pub use ink_kernel::progress::progress;

// Module map:
//   types        Status, Evidence, Check, Held, ClaimOptions and the assert report types
//   asserts      run_asserts and asserted: a document's claims, each as it is decided
//   asking       check and ask: one claim's entry, its scope, accuracy and finishing
//   deciding     decide and estimated: the method chosen for each kind of claim
//   shape        point, touch, contact, overhang, clearance and distance checks
//   estimate     sampled and summed integrals over parts, and the sample spacing
//   reads        what a claim's fields read: t, ranged params, bodies at a witness
//   convergence  grid error and ladders: how far a measurement could be off
//   metric       a package's measured metric as a claim
//   tests
mod asking;
mod asserts;
mod convergence;
mod deciding;
mod estimate;
mod metric;
mod reads;
mod shape;
mod types;

pub use asking::*;
pub use asserts::*;
pub use convergence::*;
pub(crate) use deciding::*;
pub use estimate::*;
use metric::*;
pub use reads::*;
pub use shape::*;
pub use types::*;

#[cfg(test)]
mod tests;
