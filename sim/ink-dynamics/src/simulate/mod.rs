//! Forward dynamics: a mechanism stepped forward in time from its forces.
//! The coordinates passed as `(simulate X0 [V0])` are integrated from that
//! position and speed at the start of the file's span; the others follow
//! their written motion, which also commands the joints its springs hold
//! (springs.rs): each one's deviation from what the coordinates command is
//! integrated beside them, from where the springs hold the mechanism at rest
//! at the start. Each actuator moving a simulated coordinate supplies its
//! `(command PATH EFFORT)`, clamped to its effort limits at that state (a
//! command past them is full effort: a DC motor's line at its speed), or
//! none; each spring pulls by its stretch. At a state the model gives the
//! effort each coordinate and deviation needs, f(q''), by inverse dynamics
//! (Model::state), and its generalized inertia M (Moment::gm); the
//! accelerations are the root of r(q'') = Q - f(q''), Q the commands through
//! each actuator's drive (B^T tau) and the springs' pull. f is affine in
//! them, 1e-6 M q'' + h, but where a coupling's loss turns with the power's
//! direction (then piecewise affine): Newton's method from M's step, its
//! Jacobian r's own by differences, until r is within the rounding of the
//! terms summed into it; a state where it is not found stops the
//! simulation, saying so. So does a state inverse dynamics refuses (a
//! coupling whose loss is unknown or whose tensions are not found, a body
//! whose inertia is left out: Model::forward_refusal). That is the equations
//! of motion in the mechanism's own coordinates, constraints embedded (as
//! Kane's method and Autolev form them); a loop closed by a solve is
//! followed by the solve and its implicit-function rates, so its closure
//! holds at every state and needs no stabilization (coordinate
//! partitioning, Wehage and Haug 1982). A released joint carrying no
//! inertia along itself (an idler given none) is held where the springs
//! leave it at each state (Guyan's condensation, as the modes do).
//!
//! The integrator is Dormand and Prince's 5(4) pair with error control
//! (Hairer, Norsett and Wanner, Solving ODEs I, II.5) while the motion is
//! not stiff, and Radau IIA of order 5 once its stiffness test finds it so
//! (DOPRI5's: h ||f(y7) - f(y6)|| / ||y7 - y6|| past 3.25 on 15 steps
//! running; Hairer and Wanner, Solving ODEs II, IV.8, as their RADAU5:
//! collocation at three Radau points, simplified Newton iterations on the
//! stages, its embedded error estimate and step control, starting values
//! from the last step's collocation polynomial), and back to Dormand-
//! Prince once the stiffness has passed: its Jacobian's spectral radius
//! fallen to a tenth of where Radau took over (ink's own hysteresis) and
//! the step times it within Dormand-Prince's stability by half. Switching
//! on the step times the spectral radius against the explicit method's
//! stability, over steps running, is LSODA's (Petzold 1983) and
//! DifferentialEquations.jl's AutoSwitch's. A simulation of springs is
//! stepped by Radau throughout (an explicit method's steps amplify an
//! undamped ringing near the tolerance). Radau IIA is implicit and stiffly accurate: a stiff spring's or
//! damper's modes that decay fast or that the steps do not resolve take no
//! tiny steps, where an explicit method's steps follow its fastest mode.
//! Each position's error is held to its relative tolerance and its
//! absolute one (a deviation's a nanometre at the parts it moves), a
//! speed's to its position's times its own frequency, or, once stiff, by
//! what it moves over the step. Its Jacobian along
//! the deviations is the springs' own stiffness through the inertia, and by
//! differences along the simulated coordinates. A step over which the
//! simulated coordinates' generalized inertia changes more than twofold is
//! halved, so a dead point is closed in on (FOLD), not stepped across. A step over which a command, a limit, a load or a written
//! motion may jump in t (Model::switches) is halved down to the switch, so
//! a pulse is never stepped over; where the written motion jumps (a step in
//! a commanded angle), each released joint keeps its place and speed, its
//! deviation taking the jump. Its accepted steps' states and accelerations
//! are the motion's knots, read between them by the quintic Hermite through
//! each end's position, speed and acceleration. It runs only as far as what
//! reads it asks (Simulation::reach).

// The module, by concern:
// - trajectory: the simulated motion and its quintic Hermite reading
// - types: what a simulation reports (Stepped, Springs, ...), Fail, Layout
// - forces: the accelerations at a state (Eval), the springs' hold
// - linalg: LU, Cholesky, the springs' modes, the spectral radius
// - simulation: Simulation and its progress, started and read
// - following: what the springs did (following errors, modes, rest)
// - stepping: the step loop, switches, error scales, the Jacobian
// - integrators: Dormand-Prince 5(4) and Radau IIA steps

mod following;
mod forces;
mod integrators;
mod linalg;
mod simulation;
mod stepping;
mod trajectory;
mod types;

pub use forces::accelerations;
pub use simulation::Simulation;
pub use trajectory::Trajectory;
pub use types::{Fail, Following, Layout, Ringing, Span, Springs, Stepped};

use crate::model::Model;
#[cfg(test)]
use std::sync::Mutex;

/// The integrator's relative and absolute tolerance per step (positions in
/// mm, degrees or a released joint's radians, speeds in their units a
/// second).
pub const TOLERANCE: f64 = 1e-8;

/// The least absolute tolerance per step of a released joint's deviation,
/// as the parts it carries move by it (mm: a nanometre; Simulation::scales).
pub const DEVIATION: f64 = 1e-6;

/// A deviation's tolerance per step as a share of the largest deviation the
/// springs have made so far at the parts they move (none below DEVIATION),
/// and how many of those tolerances a following error is known to within
/// either way at the start, its band (0.1 um, or 1e-4 of the largest
/// deviation so far), widening by PER_PERIOD a period stepped: measured, no
/// bound (1 kg released, by Radau IIA over a second against its closed
/// form: within 68 at 335 Hz, 11 at 50 Hz; a ringing of 10 is smoothed
/// away), so a claim within it of its bound is undecided.
pub const RESOLVE: f64 = 1e-6;
pub const UNRESOLVED: f64 = 100.0;
/// The band's widening, in tolerances a period of the springs' fastest mode
/// (at the start) stepped: Radau IIA damps an undamped ringing near the
/// tolerance by about 0.2 tolerances a period (0.1 um at 335 Hz, 68 over
/// 335 periods; 0.15 um, gone by 804 periods, 150; at 50 Hz, 147 over
/// 750), half again that.
pub const PER_PERIOD: f64 = 0.3;

/// How the motion is found, for answers.
pub const EQUATIONS: &str = "simulated forward from each (simulate X0 V0) at the span's start: the equations of motion in the mechanism's \
    coordinates (its inverse dynamics and generalized inertia at each state), each actuator moving a simulated coordinate supplying its (command \
    PATH EFFORT) clamped to its effort limits, or none";
/// How a simulation without springs is stepped.
pub const INTEGRATOR: &str = "Dormand and Prince's 5(4) pair with error control, switching to Radau IIA of order 5 (implicit) where its \
    stiffness test finds the motion stiff, and back once it is not";

/// What a simulation of springs adds to METHOD.
pub const SPRINGS_METHOD: &str = "each joint its springs hold (the file reads a following error) stepped too, from where they hold it at rest at \
    the start, as its deviation from what the coordinates command, each spring pulling by its stretch, and the whole by Radau IIA throughout (an \
    explicit method's steps grow an undamped ringing near the tolerance)";

/// Most stages a simulation evaluates, accepted or rejected, its
/// Jacobian's differences included (each the accelerations at a state: one
/// state or more): a few seconds of work on a small mechanism.
const MAX_STAGES: usize = 20_000;

/// A test's time limit: the simulation over span `.0` passes it once `.1`
/// stages are evaluated, a stop by the work done, not a clock (the call's
/// own limit stops every search in the process, every test's beside it).
/// Keyed by its span, as claims step it on rayon's workers.
#[cfg(test)]
static LIMIT_AFTER: Mutex<Vec<((f64, f64), usize)>> = Mutex::new(Vec::new());

/// While held, the simulation over `span` passes its time limit after
/// `stages` stages (LIMIT_AFTER); dropped, even by a panic, it has none.
#[cfg(test)]
pub(crate) struct LimitAfter((f64, f64));

#[cfg(test)]
impl LimitAfter {
    pub(crate) fn new(span: (f64, f64), stages: usize) -> LimitAfter {
        LIMIT_AFTER.lock().unwrap_or_else(|e| e.into_inner()).push((span, stages));
        LimitAfter(span)
    }
}

#[cfg(test)]
impl Drop for LimitAfter {
    fn drop(&mut self) {
        LIMIT_AFTER.lock().unwrap_or_else(|e| e.into_inner()).retain(|(s, _)| *s != self.0);
    }
}

/// Whether the call's time limit has passed (ink_kernel::stop) for the
/// simulation over `span`, `evaluations` stages evaluated.
fn out_of_time(span: (f64, f64), evaluations: usize) -> bool {
    #[cfg(test)]
    if let Some(&(_, n)) = LIMIT_AFTER.lock().unwrap_or_else(|e| e.into_inner()).iter().find(|(s, _)| *s == span) {
        return evaluations >= n;
    }
    let _ = (span, evaluations);
    ink_kernel::stop::reached()
}

/// The mechanism `model` simulated over the whole of `span` from `start`
/// (Simulation::start), its motion. Err: it cannot start.
pub fn simulate(model: Model, span: (f64, f64), start: &[(f64, f64)]) -> Result<Trajectory, String> {
    let sim = Simulation::start(model, span, start)?;
    sim.reach(span.1, |_| ());
    Ok(sim.progress.into_inner().unwrap_or_else(|e| e.into_inner()).traj)
}
