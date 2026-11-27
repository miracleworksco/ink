//! The dynamics package: the world's `(gravity gx gy gz)`, loads on what
//! the root places (`(load NAME PATH (force ..) (torque ..) (at FRAME))`),
//! the effort each actuator (a joint with `(effort-limits LO HI)`) supplies
//! along the written motion, and the reaction a path's bodies get from the
//! rest of the model (`(reaction PATH [FRAME])`, a value claims read), and
//! the tensions of relations' couplings (`(tension A B)`: couple.rs), by
//! inverse dynamics (model.rs). Each actuator placed under the root
//! assembly of a file that declares gravity is claimed,
//! `effort-limits:PATH`: its effort, a value this package provides at each
//! moment (a param no file writes), within its limits. Answers are
//! estimates: the masses are measured, and a motion is sampled at steps
//! chosen from it.
//! It depends on core, kernel, language and ink-carried (std's carried
//! description), never on another package (../kernel/docs/CORE.md).

#[cfg(test)]
mod answers;
mod clause;
mod command;
mod couple;
pub mod coupling;
mod effort;
pub mod joint;
mod loads;
mod mechanism;
mod model;
mod package;
mod providers;
mod reactions;
pub mod simulate;
mod simulation;
mod split;
mod steps;
mod supports;
mod tension;
#[cfg(test)]
mod tests;
#[cfg(test)]
mod validation;

pub use self::command::{Read, read, read_kinematic};
pub use self::loads::{Load, NotApplied, loads_not_applied, loads_to_ground, loads_written};
pub use self::mechanism::{build, build_kinematic, build_states, drive_rate, joints, rate_bounds};
pub use self::model::riding;
pub use self::model::*;
pub use self::package::{Gravity, NOT_CHECKED, gravity, no_gravity, not_checked};
pub(crate) use self::providers::{Sprung, spring_effort, spring_steps, with_springs};
pub use self::reactions::{Divided, divided, placement_at, reaction, reaction_bound, reaction_in, reaction_moves};
pub use self::simulation::{Command, Following, Simulated, Springs, simulation};
pub use self::supports::{sits_on, supports_floating, supports_near};
pub use ink_core::eval::Ctx;
use ink_core::graph::NodeId;
use ink_lang::InkError;

type R<T> = Result<T, InkError>;
type V3 = [f64; 3];

/// A field's node and an evaluation's context, for one who reads a drive's rate (`drive_rate`).
pub type Rate = NodeId;

/// The package.
pub struct Dynamics;
