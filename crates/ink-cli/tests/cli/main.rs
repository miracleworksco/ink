//! The `ink` command line end to end: the binary cargo built for this test
//! (CARGO_BIN_EXE_ink) run on owned fixtures and models written by the
//! tests, its JSON and output files checked independently of the kernel.
//! `models` checks every model in examples/, benchmarks/ and tests/fixtures by
//! the model rule (docs/IR.md: a model's claims are its tests).
//!
//! They are this one binary, so they all run side by side (cargo runs
//! binaries one after another): each ink call gets `THREADS` threads (a
//! call is mostly small and serial; every core each costs more CPU and is
//! slower with the tests beside it).

mod assemblies;
mod benchmarks;
mod commands;
mod common;
mod docs;
mod imports;
mod ink_path;
mod mcp;
mod models;
mod motion;
mod range;
