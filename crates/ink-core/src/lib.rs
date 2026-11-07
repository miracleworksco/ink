//! The ink core: field graphs and their tapes, interval arithmetic, the
//! evaluators, and sampled fields. It depends on nothing
//! else in ink; the kernel services and the packages build on it. See
//! docs/CORE.md for the layers.

pub mod eval;
pub mod graph;
pub mod grid;
pub mod interval;
pub mod period;
pub mod sign;
pub mod symmetry;
pub mod tape;
pub mod v3;

pub type V3 = [f64; 3];
