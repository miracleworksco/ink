//! Claims over a motion (`(t-span A B)`, `--t-span`). Time is a range of
//! the question, as a tolerance band is: a claim that reads t is the same
//! claim with t a ranged param over the span (`ranged`), decided by the
//! one search over every moment at once, so what holds or fails is
//! guaranteed as it is over any band. What that search leaves undecided
//! (within the call's time limit), or cannot ask (a solve on t not
//! followed, a package's metric, an integral), is measured at evenly
//! spaced moments (`measured`): holding there is an estimate, failing at
//! one is a witness.

use crate::InkError;
use crate::claims::{self, Check, ClaimOptions, Evidence, Status, checked};
use crate::view::ModelView;
use crate::{Claim, Cmp, Xf};
use ink_core::V3;
use ink_core::eval::Ctx;
use ink_core::graph::{BinOp, DriveRead, MotionRead, Node, UnOp, Var};
use ink_core::interval;
use ink_core::tape::Tape;
pub use ink_kernel::motion::{MAX_STEPS, Motion};
use ink_kernel::rigid::Isometry;

// Module map:
//   asked     Stepped, Search, Asked: a claim's answer and what may upgrade it
//   boxes     item boxes swept over a span
//   span      over_span: a claim over a span, in pieces and between moments
//   clock     jumps in drives and the local clock check
//   search    the ranged search over t, its scope and its parts
//   turning   a drive's period: one turn standing for every moment
//   frames    the moving frame a pair is searched in
//   pivots    relabelled fields and pivots of a frame turning about z
//   framed    Frame: a witness in the moving frame shown again in the world's
//   sampling  measured: evenly spaced moments where the search cannot decide
//   over      shown over a span, certified and banded answers, rigid groups
//   tests
mod asked;
mod boxes;
mod clock;
mod framed;
mod frames;
mod over;
mod pivots;
mod sampling;
mod search;
mod span;
mod turning;

pub use asked::*;
pub use boxes::*;
pub use clock::*;
pub(crate) use framed::*;
use frames::*;
pub use over::*;
use pivots::*;
pub use sampling::*;
pub(crate) use search::*;
pub(crate) use span::*;
use turning::*;

#[cfg(test)]
mod tests;
