//! Flat instruction tapes compiled from the DAG.
//!
//! A tape keeps two views of the same program:
//! - `ssa`: one value per instruction, arguments are instruction indices.
//!   Used for simplification.
//! - `code`: the same instructions with register-allocated slots.
//!   Used by the evaluators, so scratch memory stays small.
//!
//! The last instruction is the output.

mod compile;
mod fingerprint;
mod instruction;
mod joins;
mod kept;
mod program;
mod promise;
mod region;
mod simplify;

#[cfg(test)]
mod tests;

pub use fingerprint::Names;
pub use instruction::{Ins, Op, Pass};
pub use program::Tape;
pub use promise::Promise;
pub use region::{Sheet, Sheets};
pub use simplify::{BOTH, FIRST, OWN_TAPE_AT_MOST, SECOND};
