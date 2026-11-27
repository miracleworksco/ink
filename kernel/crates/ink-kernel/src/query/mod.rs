//! Geometric queries with proof-carrying answers: bounds, integrals (volume,
//! mass), emptiness proofs, and planar slices.

mod bounds;
mod clearance;
mod cover;
mod empty;
mod measurement;
mod proof;
mod ray;
mod rotation;
mod slice;

#[cfg(test)]
mod tests;

pub use crate::integral::{ACCURACY, Measured, Until, integrate, measured};
pub use bounds::{Extent, REACH, bounds, bounds_quick, extent, focus, max_extent, radial_extent};
pub use clearance::{MinDistance, Way, clearance, closer_pair, distance_up, inside_pair, min_distance, to_surface};
pub(crate) use cover::Cover;
pub use empty::{BOXES, Except, LOOSE, Loose, Proof, prove_empty, prove_empty_except, prove_empty_within};
pub use measurement::{occupancy, sym_eigen};
pub use proof::{Ended, Verdict, asked_moment, confirm_inside};
pub(crate) use proof::{at_moment, box_samples, signed};
pub use ray::{RayHit, Side, ray_to, ray_to_within};
pub use slice::{Loop, slice, slice_axes, slice_checked, slice_mark};
