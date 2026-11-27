//! `(study ...)`: physics studies (read through the `Reader`), their
//! metrics and results, and the `study-field`, `material` and `by-body`
//! expressions.

use crate::eval::Ctx;
use crate::graph::{NodeId, Var};
use crate::lang::*;
use crate::sexp::{SK, Sexp};
use crate::{InkError, V3};

type R<T> = Result<T, InkError>;

// Module map:
// - study:   Region, StudyPhysics, StudyFix, StudyLoad, Study, the Sim package
// - solved:  SolvedStudy, TimeError, time_error, the solve (`SolvedStudy::run`)
// - measure: a solved study's bands, mixed cells, safety factors, warnings
// - outcome: `Outcome for SolvedStudy` (metrics) and `stress_read`
// - read:    `read_study`
// - parts:   `read_load` and the material magnitudes
mod measure;
mod outcome;
mod parts;
mod read;
mod solved;
mod study;

pub(crate) use self::read::read_study;
pub use self::{outcome::*, parts::*, solved::*, study::*};
