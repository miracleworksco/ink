use super::Trajectory;
use super::following::{holds, rest};
use super::forces::{Eval, forces, said};
use super::linalg::modes_of;
use super::types::{Layout, Stepped};
use super::{DEVIATION, PER_PERIOD, RESOLVE, UNRESOLVED};
use crate::model::Model;
use std::sync::Mutex;

/// Where a simulation has got to: its motion so far, and the integrator's
/// state at its last knot.
pub(super) struct Progress {
    pub(super) traj: Trajectory,
    /// Positions then speeds, their derivative there, and the state read
    /// there (Eval).
    pub(super) y: Vec<f64>,
    pub(super) f0: Vec<f64>,
    pub(super) last: Eval,
    /// The next step's length, the last accepted one's, and the Jacobian
    /// at the last knot where it was found (reused while Newton converges
    /// fast); the last step's collocation polynomial (its divided
    /// differences: RADAU5's CONT), from which the next starts; Newton's
    /// last rate and its factor.
    pub(super) h: f64,
    pub(super) hold: f64,
    pub(super) jac: Option<Vec<Vec<f64>>>,
    pub(super) cont: Option<[Vec<f64>; 3]>,
    pub(super) theta: f64,
    pub(super) faccon: f64,
    pub(super) first: bool,
    pub(super) reject: bool,
    /// Whether the motion was found stiff (stepped by Radau IIA since), and
    /// the stiffness test's counts of stiff and calm steps running.
    pub(super) stiff: bool,
    pub(super) stiffs: usize,
    pub(super) calm: usize,
    /// The largest deviation the springs have made so far at the parts they
    /// move (mm), which the deviations' tolerance is a share of.
    pub(super) scale: f64,
    /// Its Jacobian's spectral radius where Radau took over (0: not yet
    /// found): what it must fall from to be stiff no longer.
    pub(super) rho: f64,
    /// The step's length before a switch was stepped across (resumed
    /// after it), the switches crossed, the simulated coordinates'
    /// generalized inertia at the start, and the largest speed so far.
    pub(super) resume: Option<f64>,
    pub(super) switches: usize,
    pub(super) inertia: Vec<f64>,
    pub(super) fastest: Vec<f64>,
    /// Whether the call's time limit stopped the last stepping (its stop,
    /// in traj.stopped, is not the motion's: a later call steps on).
    pub(super) limited: bool,
}

/// A mechanism (a model of states: Model::state) simulated over `span`,
/// stepped only as far as what reads it asks.
pub struct Simulation {
    pub(super) model: Model,
    pub(super) lay: Layout,
    pub(super) span: (f64, f64),
    pub(super) progress: Mutex<Progress>,
}

/// y's derivative: its speeds, then the accelerations.
pub(super) fn deriv(y: &[f64], acc: &[f64]) -> Vec<f64> {
    let n = acc.len();
    y[n..].iter().copied().chain(acc.iter().copied()).collect()
}

impl Simulation {
    /// The mechanism `model` from each simulated coordinate's position and
    /// speed `start` (in Model::simulated's order) at the span's start, and
    /// each released joint where the springs hold it at rest there (the
    /// coordinates at their start, still: what gravity and the loads need
    /// along it, by Newton's method on the springs' pull), still. Err:
    /// nothing is simulated, a body tied by its drive alone to a released
    /// joint, or no state at the start.
    pub fn start(model: Model, span: (f64, f64), start: &[(f64, f64)]) -> Result<Simulation, String> {
        let (ns, nr, t0) = (model.simulated.len(), model.released().len(), span.0);
        if ns == 0 && nr == 0 {
            return Err("nothing is simulated: no coordinate is passed as (simulate X0 [V0]), and nothing under root states a stiffness (a belt's \
                 (tensile-stiffness EA), a joint's (stiffness K)) for a following error to read"
                .into());
        }
        let nq = model.coords.len();
        let at_start = |e: String| format!("at the span's start, t = {t0}: {e}");
        // The mechanism at rest at its start.
        let sim: Vec<[f64; 3]> = start.iter().map(|s| [s.0, 0.0, 0.0]).collect();
        let still: Vec<[f64; 3]> = model.beside(t0, &sim).map_err(at_start)?.into_iter().map(|q| [q[0], 0.0, 0.0]).collect();
        let (mut w, mut y, mut d0) = (Vec::new(), Vec::new(), vec![0.0; nr]);
        let mut still_state = None;
        if nr > 0 {
            let cut = model.springs_cut(t0, &still).map_err(at_start)?;
            if !cut.is_empty() {
                let names: Vec<String> = cut.iter().take(3).map(|(b, q)| format!("{b} (with {q})")).collect();
                let more = if cut.len() > 3 { format!(" and {} more", cut.len() - 3) } else { String::new() };
                return Err(format!(
                    "{}{more} move with a coordinate that also commands joints released for springs, but not through them: tied to them only by the drives' \
                     expressions, which the springs' deviations cut (released: {}): relate it to what drives it (gear-mesh, belt-drive, lead-screw), or \
                     place it on its driver's frame",
                    names.join(", "),
                    model.released().iter().map(|r| r.0).collect::<Vec<_>>().join(", ")
                ));
            }
            let m = model.state_near(t0, &still);
            if let Some(why) = &m.unsolved {
                return Err(at_start(why.clone()));
            }
            let st = m.springs.as_ref().expect("the springs at a state of a model with springs");
            for j in 0..nr {
                if m.gm[nq + j][nq + j] > 1e-9 * st.heft[j] {
                    w.push(j);
                } else {
                    y.push(j);
                }
            }
            // Where they hold it (rest); where they cannot, no rest exists
            // to start from.
            let held = rest(st, &m.fr, &m.noise_r);
            if let Err(j) = held {
                let (d0, joints) = (rest(st, &m.fr, &vec![f64::INFINITY; nr]).unwrap_or_default(), model.released());
                let fs = st.force(&d0);
                {
                    let unit = if joints[j].1 { "N mm" } else { "N" };
                    return Err(at_start(format!(
                        "the springs do not hold {} at rest: what its parts need along it there (gravity, the loads), {:.6} {unit}, is more than its \
                         springs pull there ({:.6} {unit}: a span gone slack under it, or nothing holding it), so no rest exists to start from; hold it (a stiffness, a pretension, a command) or start it moving as written",
                        joints[j].0, m.fr[j], fs[j]
                    )));
                }
            }
            d0 = held.unwrap_or_default();
            still_state = Some(m);
        }
        let slowest = 1.0 / (span.1 - span.0);
        let (mut gyration, mut own) = (vec![1.0; ns], vec![slowest; ns]);
        if let Some(m) = &still_state {
            let st = m.springs.as_ref().expect("the springs at a state of a model with springs");
            let (joints, k) = (model.released(), st.stiffness(&d0));
            for &j in &w {
                let r = (m.gm[nq + j][nq + j] / st.mass[j]).sqrt();
                gyration.push(if joints[j].1 && r > 0.0 && r.is_finite() { r } else { 1.0 });
                own.push((1e6 * k[j][j] / m.gm[nq + j][nq + j]).sqrt().max(slowest));
            }
        }
        let mut lay = Layout { ns, w, y, nr, gyration, own, fastest: 0.0 };
        let n = lay.n();
        let mut y0: Vec<f64> = start.iter().map(|s| s.0).chain(lay.w.iter().map(|&j| d0[j])).collect();
        y0.extend(start.iter().map(|s| s.1).chain(lay.w.iter().map(|_| 0.0)));
        let ev = forces(&model, &lay, t0, &y0[..n], &y0[n..], &vec![0.0; n]).map_err(|f| at_start(said(&model, t0, &y0[n..], f)))?;
        // The springs' modes there: one nothing holds (no spring, no
        // stiffness, no actuator along it: a written turn with no effort
        // limits on a sprung belt) is refused, naming its joints, as its
        // deviation would drift with the whole command; the fastest sets
        // how a following error's band widens.
        let nw = lay.w.len();
        if nw > 0 {
            let m: Vec<Vec<f64>> = (0..nw).map(|a| (0..nw).map(|b| ev.mss[ns + a][ns + b]).collect()).collect();
            if let Some((lam, shapes)) = modes_of(&m, &ev.keff) {
                let top = lam.iter().fold(0.0f64, |a, x| a.max(*x));
                let st = still_state.as_ref().and_then(|m| m.springs.as_ref()).expect("the springs at the start");
                let unit: Vec<f64> = (0..nw).map(|a| m[a][a].sqrt()).collect();
                if let Some(q) = (0..nw).find(|&q| !holds(st, &d0, &lay, &unit, &shapes[q])) {
                    let size = |w: usize| (shapes[q][w] * lay.gyration[ns + w]).abs();
                    let most = (0..nw).map(size).fold(0.0, f64::max);
                    let joints = model.released();
                    let names: Vec<&str> = (0..nw).filter(|&w| size(w) >= 0.1 * most).map(|w| joints[lay.w[w]].0).collect();
                    return Err(at_start(format!(
                        "nothing holds {} along a motion they make together: no spring, stiffness or actuator resists it, so their following error \
                         would drift with the whole command; give the joint driving them its effort limits (an actuator keeps its coordinate) or \
                         a stiffness",
                        names.join(", ")
                    )));
                }
                lay.fastest = (1e6 * top).sqrt() / (2.0 * std::f64::consts::PI);
            }
        }
        // A simulated coordinate's own frequency at the start, from how its
        // acceleration changes with it (a load spring's k / m), its speed's
        // tolerance its position's times that: a fast mode it barely moves
        // weighs as little in its speed.
        let mut evaluations = 1;
        for i in 0..ns {
            let delta = 1e-7 * y0[i].abs().max(1.0);
            let mut q = y0[..n].to_vec();
            q[i] += delta;
            evaluations += 1;
            if let Ok(e) = forces(&model, &lay, t0, &q, &y0[n..], &ev.acc) {
                let own = ((e.acc[i] - ev.acc[i]) / delta).abs().sqrt();
                if own.is_finite() {
                    lay.own[i] = own.max(slowest);
                }
            }
        }
        let f0 = deriv(&y0, &ev.acc);
        let fastest = y0[n..].iter().map(|v| v.abs()).collect();
        let scale = (ns..n).map(|i| lay.gyration[i] * y0[i].abs()).fold(0.0, f64::max);
        let traj = Trajectory {
            t: vec![t0],
            displacement: Vec::new(),
            across: Vec::new(),
            x: vec![(0..n).map(|i| [y0[i], y0[n + i], ev.acc[i]]).collect()],
            tensions: vec![ev.tensions.clone()],
            scale: vec![scale],
            gyration: lay.gyration.clone(),
            springs: nr > 0,
            stopped: None,
            rejected: 0,
            evaluations,
        };
        let h = (span.1 - span.0) / 1000.0;
        let inertia = ev.inertia.clone();
        let progress = Progress {
            traj,
            y: y0,
            f0,
            last: ev,
            h,
            hold: h,
            jac: None,
            cont: None,
            theta: 0.0,
            faccon: 1.0,
            first: true,
            reject: false,
            stiff: nr > 0,
            stiffs: 0,
            calm: 0,
            scale,
            rho: 0.0,
            resume: None,
            switches: 0,
            inertia,
            fastest,
            limited: false,
        };
        Ok(Simulation { model, lay, span, progress: Mutex::new(progress) })
    }

    pub fn model(&self) -> &Model {
        &self.model
    }

    /// The component of the trajectory a released joint's deviation is, and
    /// whether it turns (radians, read in degrees). Err: no such joint
    /// released, or one carrying no inertia (held where its springs leave
    /// it: not stepped).
    pub fn released(&self, path: &str) -> Result<(usize, bool), String> {
        let joints = self.model.released();
        let Some(j) = joints.iter().position(|r| r.0 == path) else {
            let names: Vec<&str> = joints.iter().map(|r| r.0).collect();
            return Err(if names.is_empty() {
                format!("{path} is no joint released for springs: nothing in the simulated mechanism states a stiffness")
            } else {
                format!("{path} is no joint released for springs (those are: {})", names.join(", "))
            });
        };
        match self.lay.w.iter().position(|&w| w == j) {
            Some(i) => Ok((self.lay.ns + i, joints[j].1)),
            None => Err(format!("{path} carries no inertia along its own motion: it is held where its springs leave it at each state, and its deviation is not stepped")),
        }
    }

    /// The motion stepped on to moment `until` (or the span's end, or where
    /// it stops), and read.
    pub fn reach<T>(&self, until: f64, read: impl FnOnce(&Trajectory) -> T) -> T {
        let mut p = self.progress.lock().unwrap_or_else(|e| e.into_inner());
        self.step(&mut p, until.min(self.span.1));
        read(&p.traj)
    }

    /// Component i's position and first four derivatives at t.
    pub fn at(&self, i: usize, t: f64) -> Result<[f64; 5], String> {
        self.reach(t.max(self.span.0.next_up()), |traj| traj.at(i, t))
    }

    /// The moments a claim over `span` reading component i is measured at
    /// (Trajectory::moments), stepped to there.
    pub fn moments(&self, i: usize, span: (f64, f64)) -> Vec<f64> {
        self.reach(span.1, |traj| traj.moments(i, span))
    }

    /// What answers reading it say of it.
    pub fn note(&self) -> String {
        self.progress.lock().unwrap_or_else(|e| e.into_inner()).traj.note()
    }

    /// Component i's RMS over the moments `span`: the square root of its
    /// square's mean, by five-point Gauss-Legendre on each piece of its
    /// interpolant (a quintic: its square integrated to rounding).
    pub fn rms(&self, i: usize, span: (f64, f64)) -> Result<f64, String> {
        const G5: [(f64, f64); 5] = [
            (-0.906179845938664, 0.236926885056189),
            (-0.538469310105683, 0.478628670499366),
            (0.0, 0.568888888888889),
            (0.538469310105683, 0.478628670499366),
            (0.906179845938664, 0.236926885056189),
        ];
        self.reach(span.1, |traj| {
            let mut ends: Vec<f64> = traj.t.iter().copied().filter(|&t| t > span.0 && t < span.1).collect();
            ends.insert(0, span.0);
            ends.push(span.1);
            let mut sum = 0.0;
            for w in ends.windows(2) {
                let (mid, half) = (0.5 * (w[0] + w[1]), 0.5 * (w[1] - w[0]));
                for (x, wt) in G5 {
                    let v = traj.at(i, mid + half * x)?[0];
                    sum += wt * half * v * v;
                }
            }
            Ok((sum / (span.1 - span.0)).sqrt())
        })
    }

    /// The motion stepped over its whole span, or to where it stops: how
    /// far, and why it stopped short, if it did (Stepped).
    pub fn stepped(&self) -> Stepped {
        self.reach(self.span.1, |traj| Stepped {
            end: *traj.t.last().expect("a knot"),
            steps: traj.t.len() - 1,
            rejected: traj.rejected,
            evaluations: traj.evaluations,
            stopped: traj.stopped.clone(),
        })
    }

    /// Stepped component i's band at t (radians a turn, mm a slide): how
    /// far a following error may be from what was simulated, either way:
    /// UNRESOLVED of its tolerance (a share of the largest deviation so
    /// far), and PER_PERIOD more for each period of the springs' fastest
    /// mode stepped since the start (Radau's damping of an undamped ringing
    /// grows with the periods it steps).
    pub fn band(&self, i: usize, t: f64) -> f64 {
        let periods = self.lay.fastest * (t - self.span.0).max(0.0);
        self.reach(t, |traj| {
            let k = traj.t.partition_point(|&x| x <= t).clamp(1, traj.t.len()) - 1;
            (UNRESOLVED + PER_PERIOD * periods) * DEVIATION.max(RESOLVE * traj.scale[k]) / traj.gyration[i]
        })
    }

    /// Stepped component i's own frequency (/s: Layout::own).
    pub fn own(&self, i: usize) -> f64 {
        self.lay.own[i]
    }
}
