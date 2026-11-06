//! The design package: `(solve ...)` systems (param providers), and `(optimize ...)` (forms.rs), with the optimizers'
//! numerics behind `optimize` (optimize.rs) and the design loop that drives
//! them over a model view, and differentiates it (driver.rs). The numerics
//! work on plain parameter vectors and closures; the forms read through the
//! language.
//! It depends on core, kernel and language, never on another package
//! (docs/CORE.md).

// The crates it builds on, at the paths the package's code uses.
use ink_core::{eval, graph, tape};
use ink_lang as lang;
use ink_lang::InkError;
use ink_lang::sexp;

mod command;
mod driver;
mod forms;
pub mod optimize;
mod shape;

pub use self::driver::*;
pub use self::forms::*;
pub use self::shape::{Weight, shape_derivatives};
