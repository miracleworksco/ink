//! The model shared by dynamics claim providers and its spring effort inputs.

use crate::mechanism::build_as;
use crate::simulation::Run;
use crate::{Model, simulate};
use ink_lang::{Refused, Reload, Steps};
use std::sync::{Arc, OnceLock};

/// What every effort provider of a document shares: how to build its
/// model (a load with the root's coordinates), and the model once built.
pub(super) struct Shared {
    pub(super) reload: Reload,
    pub(super) gravity: [f64; 3],
    /// The actuators, in the order of their providers.
    pub(super) paths: Vec<String>,
    /// The paths whose reactions claims read.
    pub(super) reactions: Vec<String>,
    /// Those of them the design divides among supports, each with the supports named (`supports_near`).
    pub(super) held: Vec<(String, Vec<String>)>,
    /// The document's params' names, by index.
    pub(super) names: Vec<String>,
    pub(super) model: OnceLock<Result<(Model, Vec<usize>), Refused>>,
    /// The file's simulation, where it simulates (Run).
    pub(super) run: Option<Arc<Run>>,
}

/// What a file's simulated springs make of an actuator's effort: its joint
/// released, its spring's force (the trajectory's component of its
/// deviation, and its stiffness); none (the springs do not touch it, or it
/// is commanded: the rigid model's, or its command). The one reading of it
/// `ink effort` and the effort claims share.
pub(crate) enum Sprung {
    Spring { i: usize, k: f64 },
    No,
}

/// A document's simulation where it releases springs: None where it
/// simulates none (no simulation, or springs none).
pub(crate) fn with_springs(sim: Option<Result<Arc<simulate::Simulation>, String>>) -> Option<Result<Arc<simulate::Simulation>, String>> {
    match sim? {
        Ok(s) if s.model().released().is_empty() => None,
        other => Some(other),
    }
}

impl Sprung {
    /// What the simulation `sim` releasing springs makes of actuator `path`
    /// (`commanded`: it moves simulated coordinates alone, its command its
    /// effort). Err: its coordinates move what the springs release (its
    /// load runs through them, and its effort along the written motion is
    /// not found).
    pub(crate) fn of(sim: &simulate::Simulation, commanded: bool, path: &str) -> Result<Sprung, String> {
        if commanded {
            return Ok(Sprung::No);
        }
        let sm = sim.model();
        if let Some(j) = sm.released().iter().position(|r| r.0 == path) {
            let (i, _) = sim.released(path)?;
            return Ok(Sprung::Spring { i, k: sm.released_stiffness(j) });
        }
        let behind = sm.springs_behind(path);
        if behind.is_empty() {
            return Ok(Sprung::No);
        }
        Err(format!(
            "its coordinates move what the file's simulation releases for its springs ({}): its load runs through them, so its effort is not the rigid \
             model's along the written motion, and is not found",
            behind.join(", ")
        ))
    }
}

/// A sprung actuator's effort at t: its spring's force against its deviation.
pub(crate) fn spring_effort(sim: &simulate::Simulation, i: usize, k: f64, t: f64) -> Result<f64, String> {
    Ok(-k * sim.at(i, t)?[0])
}

/// The steps a sprung actuator's (actuator a of `model`, its spring
/// component i of stiffness k) effort is read at over `span`: the efforts'
/// steps (Model::steps, each moment an actuator's limit may jump a step of
/// its own) and where its deviation turns back (its extremes), each the
/// middle of a step of its own; a step a limit may jump in across its
/// spring's force is loose (Model::limit_undecided), as the rigid
/// effort's is.
pub(crate) fn spring_steps(model: &Model, a: usize, sim: &simulate::Simulation, i: usize, k: f64, span: (f64, f64)) -> Steps {
    let mut ends = model.steps(span).ends;
    for m in sim.moments(i, span) {
        if span.0 < m && m < span.1 {
            ends.extend([m.next_down(), m.next_up()]);
        }
    }
    ends.sort_by(f64::total_cmp);
    ends.dedup();
    let loose = ends
        .windows(2)
        .map(|w| (w[0], w[1]))
        .filter(|&s| {
            let t = ink_core::eval::span_mid(s);
            spring_effort(sim, i, k, t).is_ok_and(|e| {
                (model.limit_undecided(a, s, e) || model.rating_dip(a, s, e)) && {
                    let m = model.at(t);
                    m.lo[a] <= e && e <= m.hi[a]
                }
            })
        })
        .collect();
    Steps {
        ends,
        unchecked: Vec::new(),
        loose,
        basis: "its spring's force, measured at the efforts' steps (each moment an actuator's limit may jump a step of its own) and where it turns \
                back (its extremes), each the middle of a step of its own: an estimate between them; a step a limit may jump in across it is loose"
            .into(),
    }
}

impl Shared {
    /// The file's simulation where it releases springs: None where it
    /// simulates none (no simulation, or springs none).
    pub(super) fn springs(&self) -> Option<Result<Arc<simulate::Simulation>, String>> {
        with_springs(self.run.as_ref().map(|r| r.sim()))
    }

    /// What the file's simulated springs make of actuator a's effort (the
    /// index among `paths`). Err: its coordinates move what they release
    /// (its load runs through them, and its effort along the written motion
    /// is not found), or the simulation does not start.
    pub(super) fn sprung(&self, a: usize) -> Result<Sprung, Refused> {
        let Some(sim) = self.springs() else { return Ok(Sprung::No) };
        let sim = sim?;
        let (model, index) = self.model()?;
        Ok(Sprung::of(&sim, model.commanded(index[a]), &self.paths[a])?)
    }

    /// The model, built on first use, and the index in it of each actuator
    /// in `paths`.
    pub(super) fn model(&self) -> Result<&(Model, Vec<usize>), Refused> {
        let m = ink_kernel::once::get_or_compute(&self.model, || {
            let model = build_as(&self.reload, "root", self.gravity, 0.0, &self.reactions, true, true, &self.held)
                .map_err(|e| Refused { why: e.hint.as_ref().map_or(e.message.clone(), |h| format!("{}: {h}", e.message)), open: e.open })?;
            let index = self
                .paths
                .iter()
                .map(|p| model.actuators.iter().position(|a| a.path == *p).ok_or_else(|| format!("actuator {p} is not placed in the load with coordinates")))
                .collect::<Result<_, _>>()?;
            Ok((model, index))
        });
        m.as_ref().map_err(Clone::clone)
    }
}

/// What a reaction or tension resting on a fixed-effort drive's motion says
/// of it.
pub(super) const FIXED_SHORTFALL: &str = "where a fixed-effort drive's required effort differs from its rating within the rounding of its load (its effort claim \
     is unknown there), the motion is that drive's to that rounding, and this is read to the shortfall's share of its own rounding either way";

/// What moves with the joints whose balances give the couplings' tensions
/// that no accounted joint carries (Model::unseparated), said: the tensions
/// leave out what it exerts.
pub(super) fn unseparated(model: &Model) -> Option<String> {
    (!model.unseparated.is_empty()).then(|| {
        format!(
            "the couplings' tensions rest on the balances of the joints they reach, which leave out what these bodies moving with them exert, \
             carried by no actuated joint or balancing one (a linkage's coupler on a driven crank, a belt's own mass): {}",
            model.unseparated.join(", ")
        )
    })
}
