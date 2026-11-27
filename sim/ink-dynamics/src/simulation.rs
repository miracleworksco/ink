//! Document simulations and the parameter providers that read their motion.

use crate::mechanism::{springs_model, states_model};
#[cfg(test)]
use crate::tests::BUILDING;
use crate::{gravity, no_gravity, simulate};
use ink_core::graph::{NodeId, Var};
use ink_lang::sexp::Sexp;
use ink_lang::{Doc, Form, InkError, ParamProvider, Reader, Refused, Reload, Steps};
use std::sync::{Arc, OnceLock};

/// `(springs)`: the file's stated joint stiffnesses (`(stiffness ...)`) simulated as springs, in the
/// file run. Without it a stiffness is no spring: its joint follows its written motion rigidly. Said, not
/// inferred: a file's physics is what it says, not what it happens to read (reading a following error
/// once switched every effort to the springs').
pub struct Springs {
    pub run: bool,
    pub(super) sim: Arc<Run>,
}

impl Form for Springs {
    fn name(&self) -> &str {
        "springs"
    }
}

pub(super) const SPRINGS_USAGE: &str = "(springs): in the file run, simulate the joints its (stiffness ...) hold as springs over its (t-span A B): \
    each one's effort, reactions and (following-error PATH) are then its spring's. Without it a stiffness is no spring, and a joint \
    follows its written motion";

/// Whether the file run says `(springs)`.
pub(super) fn springs_said(r: &mut Reader) -> bool {
    r.forms::<Springs>().iter().any(|s| s.run)
}

pub(super) const SIMULATE_USAGE: &str = "(simulate X0 [V0]): a coordinate's motion simulated forward from position X0 and speed V0 (default 0) at the start of \
    the file's (t-span A B), as the argument of an instance the root places, itself: (instance arm (arm (simulate 60))); the mechanism moves under \
    gravity, its loads and its actuators' (command PATH EFFORT), the other coordinates as written";

pub(super) const COMMAND_USAGE: &str = "(command PATH EFFORT): the effort the actuator PATH (a joint with effort-limits, by its instance path) supplies in a \
    simulation, N mm a turn or N a slide, a field of params, t and the simulated coordinates (a def of (simulate ...) and its (diff ... t)), clamped \
    to its effort limits there: a command past them is full effort, (command root/m/rotor 1e9) a DC motor at full voltage";

pub(super) const FOLLOWING_USAGE: &str = "(following-error PATH): in a simulation of the file's mechanism over its (t-span A B), how far joint PATH (one a \
    stiff coupling line's end rests on, or stating (stiffness K): a belt-driven carriage, a stepper's rotor) is past what the coordinates command \
    it to, its springs' deviation: mm a slide, deg a turn; a file reading one is simulated, its coordinates' written motion the command";

/// The run a document's simulation is shared through: its simulated
/// coordinates' and its following errors' (the first made).
pub(super) fn run_of(r: &mut Reader) -> Arc<Run> {
    if let Some(f) = r.forms::<Simulated>().first() {
        return f.run.clone();
    }
    if let Some(f) = r.forms::<Springs>().first() {
        return f.sim.clone();
    }
    if let Some(f) = r.forms::<Following>().first() {
        return f.run.clone();
    }
    Arc::new(Run { reload: r.reloader(), names: OnceLock::new(), inputs: OnceLock::new(), sim: OnceLock::new() })
}

/// Six params a simulation sets at each moment, named `what`'s position,
/// rate, acceleration, jerk, snap and moment t0, and the field of their
/// Taylor polynomial about t0: x + v dt + a dt^2 / 2 + j dt^3 / 6 + s dt^4
/// / 24, by Horner (exact at t0, its derivatives along t to the fourth too).
pub(super) fn taylor(r: &mut Reader, what: &str) -> ([usize; 6], NodeId) {
    let owns: [usize; 6] = ["position", "rate", "acceleration", "jerk", "snap", "moment"].map(|w| r.derived(&format!("{what} {w}")));
    owns.iter().for_each(|&k| r.solved(k));
    let g = r.graph();
    let [x, v, a, j, s, t0] = owns.map(|k| g.var(Var::Param(k as u32)));
    let t = g.var(Var::T);
    let dt = g.sub(t, t0);
    let mut field = g.constant(0.0);
    for (c, w) in [(s, 1.0 / 24.0), (j, 1.0 / 6.0), (a, 0.5), (v, 1.0), (x, 1.0)] {
        let cw = g.constant(w);
        let term = g.mul(cw, c);
        let carried = g.mul(field, dt);
        field = g.add(carried, term);
    }
    (owns, field)
}

/// `(following-error PATH)` in an expression: a released joint's deviation
/// in the file's simulation (springs.rs), its Taylor polynomial's params as
/// a simulated coordinate's (Simulated).
pub struct Following {
    pub path: String,
    pub owns: [usize; 6],
    pub(super) run: Arc<Run>,
}

impl Form for Following {
    fn name(&self) -> &str {
        "following-error"
    }

    fn integration_inputs(&self) -> &[NodeId] {
        self.run.inputs.get().map_or(&[], Vec::as_slice)
    }

    fn provider(&self, doc: &Doc) -> Option<Arc<dyn ParamProvider>> {
        (!doc.states()).then(|| Arc::new(Simulating { run: self.run.clone(), of: Of::Joint(self.path.clone()), owns: self.owns }) as Arc<dyn ParamProvider>)
    }
}

/// `(command PATH EFFORT)`: what actuator PATH supplies in a simulation.
pub struct Command {
    pub path: String,
    pub effort: NodeId,
    pub run: bool,
    pub(super) form: Sexp,
}

impl Form for Command {
    fn name(&self) -> &str {
        &self.path
    }
}

/// `(simulate X0 [V0])` in an expression: its Taylor polynomial about t0 to
/// the fourth order, x + v (t - t0) + a (t - t0)^2 / 2 + j (t - t0)^3 / 6 +
/// s (t - t0)^4 / 24, its six params the simulated coordinate's position,
/// rate, acceleration, jerk and snap at the moment t0 (set at each moment:
/// exact there, its derivatives along t to the fourth too), its place among
/// the file's simulated coordinates, and the simulation they share. Where
/// it is written.
pub struct Simulated {
    pub field: NodeId,
    pub(super) speed: NodeId,
    pub owns: [usize; 6],
    pub line: u32,
    /// Its start: position and speed, fields of params.
    pub(super) x0: NodeId,
    pub(super) v0: NodeId,
    pub(super) run: Arc<Run>,
    pub(super) index: usize,
}

impl Form for Simulated {
    fn name(&self) -> &str {
        "simulate"
    }

    fn integrated_on(&self, field: NodeId, span: (f64, f64)) -> Option<Result<(), String>> {
        if field != self.field && field != self.speed {
            return None;
        }
        Some(self.run.sim().and_then(|sim| {
            sim.at(self.index, span.0)?;
            sim.at(self.index, span.1)?;
            Ok(())
        }))
    }

    fn integration_inputs(&self) -> &[NodeId] {
        self.run.inputs.get().map_or(&[], Vec::as_slice)
    }

    fn integrated_coordinate(&self) -> Option<NodeId> {
        Some(self.field)
    }

    /// Provided at each moment, but in a load with states (the model the
    /// simulation steps, which sets them as its state).
    fn provider(&self, doc: &Doc) -> Option<Arc<dyn ParamProvider>> {
        (!doc.states()).then(|| Arc::new(Simulating { run: self.run.clone(), of: Of::Coordinate(self.index), owns: self.owns }) as Arc<dyn ParamProvider>)
    }
}

/// The simulation of a document's mechanism, where it passes a (simulate
/// ...) or reads a (following-error ...): stepped as far as it is read.
pub fn simulation(doc: &Doc) -> Option<Result<Arc<simulate::Simulation>, String>> {
    let run = doc
        .forms::<Simulated>()
        .next()
        .map(|s| s.run.clone())
        .or_else(|| doc.forms::<Springs>().find(|s| s.run).map(|s| s.sim.clone()))
        .or_else(|| doc.forms::<Following>().next().map(|f| f.run.clone()))?;
    Some(run.sim())
}

/// A document's simulation: how to load it again, its params' names, and
/// the simulation once found (shared by every load asking the same
/// question: SIMULATIONS).
pub(super) struct Run {
    reload: Option<Reload>,
    pub(super) names: OnceLock<Vec<String>>,
    pub(super) inputs: OnceLock<Vec<NodeId>>,
    sim: OnceLock<Result<Arc<simulate::Simulation>, String>>,
}

/// The simulations this process has made, by their question: what the file
/// was read from (its sources, imports and data included: Doc::sources) and
/// its param overrides. A document and the model of its coordinates its
/// efforts are read from ask the same one, and step one motion. Each one a
/// document holds is kept, and of the rest the last few (a long-running
/// process, a viewer reloading).
static SIMULATIONS: std::sync::Mutex<Vec<(u128, Arc<simulate::Simulation>)>> = std::sync::Mutex::new(Vec::new());

impl Run {
    /// The mechanism as states, ready to simulate over the file's span (it
    /// steps only as far as it is read: Simulation::reach), once.
    pub(super) fn sim(&self) -> Result<Arc<simulate::Simulation>, String> {
        ink_kernel::once::get_or_compute(&self.sim, || self.made()).clone()
    }

    /// `sim`'s simulation: one made already for the same question
    /// (SIMULATIONS), else made here. Made outside the lock on the list,
    /// which guards only the list: building the model measures its masses
    /// on rayon's workers, and a worker waiting there runs other work on
    /// its stack, another file's load among it, which asks for this lock
    /// (ink_kernel::once).
    fn made(&self) -> Result<Arc<simulate::Simulation>, String> {
        let reload = self.reload.as_ref().ok_or("a simulation needs the file's source, to load it again")?;
        let said = |e: InkError| e.hint.as_ref().map_or(e.message.clone(), |h| format!("{}: {h}", e.message));
        let doc = reload.with_states("root").map_err(said)?;
        let mut fp = ink_core::graph::Fingerprint::default();
        fp.key(doc.sources());
        for (name, v) in &reload.overrides {
            fp.str(name).word(v.value.to_bits()).str(v.unit.map_or("", |u| u.name));
        }
        let key = fp.value();
        let made = || SIMULATIONS.lock().unwrap_or_else(|e| e.into_inner());
        let found = |made: &[(u128, Arc<simulate::Simulation>)]| made.iter().find(|(k, _)| *k == key).map(|(_, sim)| sim.clone());
        if let Some(sim) = found(&made()) {
            return Ok(sim);
        }
        #[cfg(test)]
        BUILDING.take().into_iter().for_each(|f| f());
        let g = gravity(&doc).ok_or_else(|| said(no_gravity()))?;
        let span = doc.t_span.ok_or_else(|| format!("a simulation runs over the file's span: it declares no (t-span A B) ({SIMULATE_USAGE})"))?;
        let (ctx, _) = doc.ctx_at(&reload.overrides, span.0).map_err(|e| e.to_string())?;
        let starts: Vec<(f64, f64)> = doc.forms::<Simulated>().map(|s| (doc.scalar(s.x0, &ctx), doc.scalar(s.v0, &ctx))).collect();
        // Its springs, where the file says (springs) (a stepper's stiffness, beside a simulated arm,
        // stays its effort's rigid line otherwise).
        let springs = doc.forms::<Springs>().any(|s| s.run);
        let model = if springs { springs_model(doc, reload, "root", g) } else { states_model(doc, reload, "root", g) }.map_err(said)?;
        let sim = Arc::new(simulate::Simulation::start(model, span, &starts)?);
        let mut made = made();
        // Made meanwhile on another thread: that one, so both step one motion.
        if let Some(sim) = found(&made) {
            return Ok(sim);
        }
        // Kept: 16, with this one. Dropped first, oldest first, every one
        // no document holds past those: one held is still being asked
        // (before, under load a test's neighbours made 16 others
        // meanwhile and its coordinates' model stepped another), and
        // the cache shrinks back once documents close (before, it
        // dropped one a call and kept its peak).
        let mut k = 0;
        while made.len() >= 16 && k < made.len() {
            if Arc::strong_count(&made[k].1) == 1 {
                made.remove(k);
            } else {
                k += 1;
            }
        }
        made.push((key, sim.clone()));
        Ok(sim)
    }
}

/// What a simulation's provider reads of it: a simulated coordinate (its
/// place among the file's), or a released joint's deviation (its path).
enum Of {
    Coordinate(usize),
    Joint(String),
}

/// One simulated coordinate's or following error's provider: what it
/// reads, and its six params (Simulated, Following).
struct Simulating {
    run: Arc<Run>,
    of: Of,
    owns: [usize; 6],
}

impl Simulating {
    /// Its component of the trajectory, and the factor its units are read
    /// in (a turn's deviation, radians stepped, is read in degrees).
    fn component(&self, sim: &simulate::Simulation) -> Result<(usize, f64), String> {
        match &self.of {
            Of::Coordinate(i) => Ok((*i, 1.0)),
            Of::Joint(path) => sim.released(path).map(|(i, turn)| (i, if turn { 180.0 / std::f64::consts::PI } else { 1.0 })),
        }
    }
}

impl ParamProvider for Simulating {
    fn name(&self) -> &str {
        "a coordinate's motion simulated forward from its forces ((simulate X0 V0))"
    }

    fn owns(&self) -> &[usize] {
        &self.owns
    }

    fn uses_t(&self) -> bool {
        true
    }

    fn warm_starts(&self) -> bool {
        false
    }

    /// The motion at t, stepped to there: its position and first four
    /// derivatives, and t.
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
        let sim = self.run.sim()?;
        let (i, unit) = self.component(&sim)?;
        let x = sim.at(i, t)?;
        for (&k, v) in self.owns.iter().zip(x.into_iter().map(|x| unit * x).chain([t])) {
            params[k] = v;
        }
        let model = sim.model();
        Ok((!model.unsettled.is_empty()).then(|| format!("masses it rests on are not measured to 2%: {}", model.unsettled.join("; "))))
    }

    /// A claim reading it is measured at the 32 even steps a measured claim
    /// takes (Motion), the turning points of position and speed, and the
    /// span extrema of the higher derivatives.
    fn steps(&self, _params: &[f64], span: (f64, f64)) -> Option<Result<Steps, Refused>> {
        Some(self.run.sim().and_then(|sim| {
            let (i, _) = self.component(&sim)?;
            let even = (0..=EVEN).map(|k| span.0 + (span.1 - span.0) * k as f64 / EVEN as f64);
            let mut ends: Vec<f64> = even.collect();
            sim.reach(span.1, |traj| {
                for t in traj.moments(i, span).into_iter().chain(traj.extremes(i, span, 2)).chain(traj.derivative_peaks(i, span)) {
                    if span.0 < t && t < span.1 {
                        ends.extend([t.next_down(), t.next_up()]);
                    }
                }
            });
            ends.sort_by(f64::total_cmp);
            ends.dedup();
            Ok(Steps {
                ends,
                unchecked: Vec::new(),
                loose: Vec::new(),
                basis: format!(
                    "measured at the middles of {EVEN} even steps and at the turning points of each simulated coordinate or following error and the extrema of its speed, with the minimum and maximum of its acceleration, jerk and snap, each the middle of a step of its own: an estimate between them"
                ),
            })
        })
        .map_err(Refused::from))
    }

    /// A following error is an estimate known to within its band either
    /// way (Simulation::band, at the moment it is read at): a value claim
    /// holds where it holds over the band, fails where it fails over all of
    /// it, and is undecided between; its derivatives to within that times
    /// its own frequency's powers (a ringing there).
    fn error(&self, params: &[f64], _t: f64) -> Option<Vec<f64>> {
        let Of::Joint(_) = &self.of else { return None };
        let sim = self.run.sim().ok()?;
        let (i, unit) = self.component(&sim).ok()?;
        let (w, f) = (unit * sim.band(i, params[self.owns[5]]), sim.own(i));
        Some(vec![w, w * f, w * f * f, w * f.powi(3), w * f.powi(4), 0.0])
    }

    fn reads(&self) -> Option<Vec<usize>> {
        let sim = self.run.sim().ok()?;
        let names = self.run.names.get()?;
        let model = sim.model();
        Some(model.reads.iter().filter_map(|n| names.iter().position(|m| m == n)).collect())
    }

    fn note(&self) -> Option<String> {
        Some(match self.run.sim.get() {
            Some(Ok(sim)) => sim.note(),
            _ => format!("the motion of each (simulate ...) coordinate is an estimate, {}; {}", simulate::EQUATIONS, simulate::INTEGRATOR),
        })
    }
}

/// The even steps a claim over a simulated motion is measured at, beside
/// its coordinates' extremes (as a measured claim's: Motion).
pub(super) const EVEN: usize = 32;
