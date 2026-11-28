use super::Simulation;
use super::forces::{Eval, forces, said};
use super::linalg::{radius, solve};
use super::simulation::{Progress, deriv};
use super::types::{Fail, Layout};
use super::{DEVIATION, MAX_STAGES, RESOLVE, TOLERANCE, out_of_time};
use crate::model::Model;

impl Simulation {
    /// Why it stops at moment t, positions q and speeds v: `why`, unless it
    /// is coming to rest (its speeds within 1% of the fastest so far) where
    /// at rest no acceleration agrees with its couplings' losses: then
    /// friction holds it.
    pub(super) fn stop(&self, p: &Progress, t: f64, q: &[f64], v: &[f64], why: String) -> String {
        let slow = v.iter().zip(&p.fastest).all(|(x, peak)| x.abs() <= 1e-2 * peak);
        let rest = vec![0.0; q.len()];
        match forces(&self.model, &self.lay, t, q, &rest, &vec![0.0; q.len()]) {
            Err(f @ Fail::Unbalanced { .. }) if slow => format!("it comes to rest (speeds {v:?}): {}", said(&self.model, t, &rest, f)),
            _ => why,
        }
    }

    /// f's Jacobian at the last knot: along the stepped deviations the
    /// springs' own through the inertia, 1e6 M^-1 (0, -K) (no spring reads
    /// their speeds), along the simulated coordinates and their speeds by
    /// differences. Err: a state the differences step to, why.
    pub(super) fn jacobian(&self, p: &mut Progress, t: f64) -> Result<Vec<Vec<f64>>, (f64, Vec<f64>, String)> {
        jacobian_at(&self.model, &self.lay, t, &p.y, &p.last, &mut p.traj.evaluations)
    }

    /// Steps on from the last knot until the motion reaches `until`: by
    /// Dormand and Prince's 5(4) pair until its stiffness test finds the
    /// motion stiff, then by Radau IIA (Simulation::radau).
    pub(super) fn step(&self, p: &mut Progress, until: f64) {
        let (model, lay, (a, b)) = (&self.model, &self.lay, self.span);
        let (ns, n) = (lay.ns, lay.n());
        let finest = 1e-9 * (b - a);
        let mut t = *p.traj.t.last().expect("a knot");
        // A moment it reached before the limit is read as it is.
        if t < until && std::mem::take(&mut p.limited) {
            p.traj.stopped = None;
        }
        while t < until && p.traj.stopped.is_none() {
            // Past the call's time limit (ink_kernel::stop) it stops where
            // it got to, in time to answer with that, and steps on from
            // there if asked again with time left.
            if out_of_time(self.span, p.traj.evaluations) {
                p.traj.stopped = Some((t, format!("the call's time limit passed at t = {t} (--timeout, or ink mcp's limit): a longer limit steps further")));
                p.limited = true;
                break;
            }
            if p.traj.evaluations > MAX_STAGES {
                let cause = if p.switches * 8 > p.traj.t.len() {
                    format!("what drives it switched {} times, each stepped across", p.switches)
                } else {
                    "a mechanism whose motion its steps must follow in fine detail (a fast oscillation they resolve, or a force far past what it moves)".into()
                };
                let why = format!("past {MAX_STAGES} stages evaluated, {} steps accepted and {} rejected, its steps {:.3e} s ({cause})", p.traj.t.len() - 1, p.traj.rejected, p.h);
                let why = self.stop(p, t, &p.y[..n], &p.y[n..], why);
                p.traj.stopped = Some((t, why));
                break;
            }
            let mut h = p.h.min(b - t);
            // A switch in what drives the motion (a command, a limit, a
            // load or a written motion jumping in t) is found by bisection
            // on the switch test alone: a step to just before it, one
            // across it at the finest, then on at the length before.
            let sim: Vec<[f64; 3]> = (0..ns).map(|i| [p.y[i], p.y[n + i], p.last.acc[i]]).collect();
            let mut across = false;
            if h > finest && model.switches((t, t + h), &sim) {
                let (mut lo, mut hi) = (t, t + h);
                while hi - lo > finest {
                    let mid = ink_core::eval::span_mid((lo, hi));
                    if mid == lo || mid == hi {
                        break;
                    }
                    if model.switches((lo, mid), &sim) {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                }
                if lo > t {
                    p.resume.get_or_insert(p.h);
                    h = lo - t;
                } else {
                    h = hi - t;
                    across = true;
                }
            }
            if t + h == t {
                let why = self.stop(p, t, &p.y[..n], &p.y[n..], "its steps shrank to nothing".into());
                p.traj.stopped = Some((t, why));
                break;
            }
            let t1 = if b - (t + h) < 1e-12 * (b - a) { b } else { t + h };
            let attempt = if p.stiff { self.radau(p, t, h, t1, across) } else { self.dopri(p, t, h, t1, across) };
            let (mut delta, mut ev, mut hnew, cont) = match attempt {
                Attempt::Failed(at, ys, why) => {
                    // A stage past where the mechanism goes on: a shorter
                    // step, down to nothing, finds where.
                    p.traj.rejected += 1;
                    if h < finest || across {
                        let why = self.stop(p, at, &ys[..n], &ys[n..], why);
                        p.traj.stopped = Some((at, why));
                        break;
                    }
                    p.h = h * 0.25;
                    p.reject = true;
                    continue;
                }
                Attempt::Reject(hnew) => {
                    p.traj.rejected += 1;
                    p.reject = true;
                    p.h = hnew;
                    continue;
                }
                Attempt::Accept { delta, knot, hnew, cont } => (delta, knot, hnew, cont),
            };
            let mut y1: Vec<_> = p.y.iter().zip(&delta).map(|(y, dy)| y + dy).collect();
            delta.truncate(n);
            // Where the written motion jumped across the step, each released
            // joint keeps its place and speed: its deviation takes all the
            // jump, what its command moved past its own rate over the step
            // (whatever the command's origin: a motor's angle after
            // thousands of turns steps alike).
            if across && !lay.w.is_empty() {
                for (j, &w) in lay.w.iter().enumerate() {
                    let (d0, d1) = (p.last.drives[w], ev.drives[w]);
                    delta[ns + j] -= d1[0] - d0[0] - d0[1] * h - 0.5 * d0[2] * h * h;
                    y1[ns + j] = p.y[ns + j] + delta[ns + j];
                    y1[n + ns + j] -= d1[1] - d0[1] - d0[2] * h;
                }
                p.traj.evaluations += 1;
                ev = match forces(model, lay, t1, &y1[..n], &y1[n..], &ev.acc) {
                    Ok(ev) => ev,
                    Err(f) => {
                        let why = said(model, t1, &y1[n..], f);
                        let why = self.stop(p, t1, &y1[..n], &y1[n..], why);
                        p.traj.stopped = Some((t1, why));
                        break;
                    }
                };
            }
            p.cont = cont;
            p.hold = h;
            p.f0 = deriv(&y1, &ev.acc);
            for (peak, v) in p.fastest.iter_mut().zip(&y1[n..]) {
                *peak = peak.max(v.abs());
            }
            p.traj.t.push(t1);
            p.traj.displacement.push(delta);
            p.traj.across.push(across);
            p.traj.x.push((0..n).map(|i| [y1[i], y1[n + i], ev.acc[i]]).collect());
            p.traj.tensions.push(ev.tensions.clone());
            p.scale = (ns..n).map(|i| lay.gyration[i] * y1[i].abs()).fold(p.scale, f64::max);
            p.traj.scale.push(p.scale);
            p.y = y1;
            p.last = ev;
            self.carry(p, t1 - t);
            if p.reject {
                hnew = hnew.min(h);
            }
            p.reject = false;
            p.first = false;
            t = t1;
            if across {
                p.jac = None;
                p.cont = None;
                p.switches += 1;
                p.h = p.resume.take().unwrap_or(hnew);
                continue;
            }
            p.h = hnew;
            // A rigid motion stiff no longer (a damper's moment past, a
            // contact let go): its Jacobian found afresh every 15 of Radau's
            // steps, and where its spectral radius has fallen to a tenth of
            // what it was where Radau took over (ink's own hysteresis) and the
            // step times it would keep Dormand-Prince stable by half (under
            // 1.6), back to it (the step times the spectral radius against
            // the explicit method's stability: LSODA's and AutoSwitch's
            // test). A simulation of springs, whose undamped modes an
            // explicit method amplifies, stays on Radau.
            if p.stiff && lay.nr == 0 {
                p.calm += 1;
                if p.calm >= 15 {
                    p.calm = 0;
                    p.jac = None;
                    if let Ok(j) = self.jacobian(p, t) {
                        let r = radius(&j);
                        if r < 0.1 * p.rho && h * r < 1.6 {
                            (p.stiff, p.stiffs, p.cont, p.rho) = (false, 0, None, 0.0);
                        }
                        p.jac = Some(j);
                    }
                }
            }
        }
    }

    /// The masses' sensitivity carried over the step h long just accepted,
    /// to its knot (the last), and its band over the step kept (none kept
    /// once lost: infinite there, Trajectory::mass_band).
    fn carry(&self, p: &mut Progress, h: f64) {
        let t = *p.traj.t.last().expect("a knot");
        match p.sens.as_mut().map(|s| s.step(&self.model, &self.lay, t, h, &p.y, &p.last)) {
            Some(Ok(groups)) => p.traj.push_mass(groups),
            Some(Err(why)) => (p.sens, p.traj.lost) = (None, Some((t, why))),
            None => {}
        }
    }

    /// Each component's error scale at a step of length h: its relative
    /// TOLERANCE, and its absolute one, a simulated coordinate's TOLERANCE,
    /// a deviation's a share RESOLVE of the largest deviation the springs
    /// have made so far at the parts they move (none smaller than
    /// DEVIATION), over its radius of gyration (a long arm's angle held
    /// finer than a pulley's); a speed's that times its own frequency, or,
    /// once the motion is stiff, by what it moves over the step (its
    /// position's over h, as RADAU5 scales a higher index's), so a fast mode
    /// the motion barely moves, its speed no small thing, is stepped over.
    pub(super) fn scales(&self, p: &Progress, h: f64) -> Vec<f64> {
        let (n, ns, lay) = (self.lay.n(), self.lay.ns, &self.lay);
        let reach = DEVIATION.max(RESOLVE * p.scale);
        (0..2 * n)
            .map(|i| {
                let k = i % n;
                let absolute = if k < ns { TOLERANCE } else { reach / lay.gyration[k] };
                match (i < n, p.stiff) {
                    (true, _) => TOLERANCE * p.y[i].abs() + absolute,
                    (false, true) => (TOLERANCE * p.y[k].abs() + absolute) / h,
                    (false, false) => TOLERANCE * p.y[i].abs() + absolute * lay.own[k],
                }
            })
            .collect()
    }

    /// The simulated coordinates' generalized inertia changing more than
    /// twofold over a step from the last knot, as its stages read it: a step
    /// half as long. What the motion rests on is followed, so a dead point,
    /// where it grows without bound, is closed in on and found (FOLD), not
    /// stepped across (the coordinate turning back there smoothly).
    pub(super) fn grown(&self, p: &Progress, stages: &[Eval], h: f64, across: bool) -> Option<Attempt> {
        let grows = stages.iter().any(|e| e.inertia.iter().zip(&p.last.inertia).any(|(v, last)| *v > 2.0 * last || *last > 2.0 * v));
        (grows && !across).then_some(Attempt::Reject(h * 0.5))
    }
}

/// f's Jacobian at state y, moment t, read there as `base` (Simulation::
/// jacobian), each state its differences step to counted in `evaluations`.
pub(super) fn jacobian_at(model: &Model, lay: &Layout, t: f64, y: &[f64], base: &Eval, evaluations: &mut usize) -> Result<Vec<Vec<f64>>, (f64, Vec<f64>, String)> {
    let n = lay.n();
    let nn = 2 * n;
    let mut jac = vec![vec![0.0; nn]; nn];
    for i in 0..n {
        jac[i][n + i] = 1.0;
    }
    for j in 0..lay.w.len() {
        let mut rhs = vec![0.0; n];
        for l in 0..lay.w.len() {
            rhs[lay.ns + l] = -1e6 * base.keff[l][j];
        }
        let col = solve(&base.mss, &rhs).unwrap_or_else(|| vec![0.0; n]);
        for i in 0..n {
            jac[n + i][lay.ns + j] = col[i];
        }
    }
    for c in (0..lay.ns).chain(n..n + lay.ns) {
        // The accelerations' own rounding is some 1e-14 of them.
        let delta = 1e-7 * y[c].abs().max(1.0);
        let mut y = y.to_vec();
        y[c] += delta;
        *evaluations += 1;
        let ev = forces(model, lay, t, &y[..n], &y[n..], &base.acc).map_err(|f| (t, y.clone(), said(model, t, &y[n..], f)))?;
        for i in 0..n {
            jac[n + i][c] = (ev.acc[i] - base.acc[i]) / delta;
        }
    }
    Ok(jac)
}

/// What a step's attempt comes to: its end state, the knot's state read
/// there, the next step's length and Radau's collocation polynomial (for
/// the next's starting values); a shorter step; or a state past where the
/// mechanism goes on (when, where, why).
pub(super) enum Attempt {
    Accept { delta: Vec<f64>, knot: Eval, hnew: f64, cont: Option<[Vec<f64>; 3]> },
    Reject(f64),
    Failed(f64, Vec<f64>, String),
}

/// Why a stage at `at` is a dead point: the simulated coordinates'
/// generalized inertia grown `fold`-fold since the start.
pub(super) fn folded(at: f64, fold: f64) -> String {
    format!(
        "a dead point at t = {at}: the simulated coordinates' generalized inertia has grown {fold:.3e}-fold since the start (a loop folding there: \
         they cannot carry the motion on)"
    )
}
