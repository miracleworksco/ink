//! The library end to end, from source to answer: conformance to
//! docs/CORE.md, the kernel, standard-library geometry, solves against
//! closed forms, and the same answers at any thread count.

#[path = "../common/mod.rs"]
mod common;

mod conformance;
mod geometry;
mod kernel;
mod solves;
mod threads;
