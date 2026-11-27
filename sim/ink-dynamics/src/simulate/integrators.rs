use super::Simulation;
use super::forces::{Eval, forces, said};
use super::linalg::{Lu, radius};
use super::simulation::{Progress, deriv};
use super::stepping::{Attempt, folded};

/// Dormand and Prince's 5(4) pair: the nodes, the stages, and the fifth-
/// order weights less the fourth's (the error estimate). The fifth-order
/// solution is the last stage's (first same as last).
const DP_C: [f64; 7] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
const DP_A: [[f64; 6]; 7] = [
    [0.0; 6],
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0],
    [9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0],
    [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0],
];
const DP_E: [f64; 7] = [71.0 / 57600.0, 0.0, -71.0 / 16695.0, 71.0 / 1920.0, -17253.0 / 339200.0, 22.0 / 525.0, -1.0 / 40.0];

/// Radau IIA of order 5: its nodes and stages (Hairer and Wanner, Solving
/// ODEs II, IV.5, Table 5.6), the real eigenvalue of A's inverse and the
/// stages' weights in the embedded error estimate (RADAU5's U1, DD1-3).
const SQ6: f64 = 2.449489742783178;
const C: [f64; 3] = [(4.0 - SQ6) / 10.0, (4.0 + SQ6) / 10.0, 1.0];
const A: [[f64; 3]; 3] = [
    [(88.0 - 7.0 * SQ6) / 360.0, (296.0 - 169.0 * SQ6) / 1800.0, (-2.0 + 3.0 * SQ6) / 225.0],
    [(296.0 + 169.0 * SQ6) / 1800.0, (88.0 + 7.0 * SQ6) / 360.0, (-2.0 - 3.0 * SQ6) / 225.0],
    [(16.0 - SQ6) / 36.0, (16.0 + SQ6) / 36.0, 1.0 / 9.0],
];
const U1: f64 = 3.637834252744496;
const DD: [f64; 3] = [-(13.0 + 7.0 * SQ6) / 3.0, (-13.0 + 7.0 * SQ6) / 3.0, -1.0 / 3.0];
/// Most Newton iterations a step takes; where they stop, as a share of the
/// tolerance (Hairer and Wanner's kappa, 0.01 to 0.1); the convergence rate
/// past which the Jacobian is found again (RADAU5's THET where it is
/// costly).
const NEWTON: usize = 7;
const KAPPA: f64 = 0.03;
const THET: f64 = 0.1;

/// How much the simulated coordinates' generalized inertia may grow from
/// the start before the motion is taken to have reached a dead point of a
/// loop (where it grows without bound: the coordinates cannot carry the
/// motion on). Scale-free: a ratio of the same quantity.
const FOLD: f64 = 1e8;

impl Simulation {
    /// One step of Dormand and Prince's 5(4) pair from the last knot (its
    /// last stage, at the step's end, the knot's state), with Hairer's
    /// stiffness test on its last two stages, both at the step's end: h
    /// times ||f(y7) - f(y6)|| / ||y7 - y6|| past 3.25 on 15 accepted steps
    /// running (6 calm ones resetting the count) finds the motion stiff, and
    /// Radau IIA steps it on (DOPRI5's test; back: Simulation::step).
    pub(super) fn dopri(&self, p: &mut Progress, t: f64, h: f64, t1: f64, across: bool) -> Attempt {
        let (model, lay) = (&self.model, &self.lay);
        let n = lay.n();
        let nn = 2 * n;
        let scal = self.scales(p, h);
        let mut k: Vec<Vec<f64>> = vec![p.f0.clone()];
        let mut ys_all: Vec<Vec<f64>> = vec![p.y.clone()];
        let mut evals: Vec<Eval> = Vec::new();
        let mut delta = vec![0.0; nn];
        for s in 1..7 {
            for (i, d) in delta.iter_mut().enumerate() {
                *d = h * (0..s).map(|j| DP_A[s][j] * k[j][i]).sum::<f64>();
            }
            let ys: Vec<f64> = p.y.iter().zip(&delta).map(|(y, dy)| y + dy).collect();
            p.traj.evaluations += 1;
            let at = if s == 6 { t1 } else { t + DP_C[s] * h };
            let guess = evals.last().map_or(&p.last.acc, |e| &e.acc).clone();
            match forces(model, lay, at, &ys[..n], &ys[n..], &guess) {
                Ok(ev) if let Some(r) = ev.inertia.iter().zip(&p.inertia).map(|(v, start)| v / start).find(|&r| r > FOLD) => return Attempt::Failed(at, ys, folded(at, r)),
                Ok(ev) => {
                    k.push(deriv(&ys, &ev.acc));
                    ys_all.push(ys);
                    evals.push(ev);
                }
                Err(f) => {
                    let why = said(model, at, &ys[n..], f);
                    return Attempt::Failed(at, ys, why);
                }
            }
        }
        if let Some(r) = self.grown(p, &evals, h, across) {
            return r;
        }
        let err = ((0..nn)
            .map(|i| {
                let e = h * (0..7).map(|j| DP_E[j] * k[j][i]).sum::<f64>();
                (e / scal[i]) * (e / scal[i])
            })
            .sum::<f64>()
            / nn as f64)
            .sqrt();
        let hnew = h * (0.9 * err.max(1e-10).powf(-0.2)).clamp(0.2, 5.0);
        if !(err <= 1.0 || across) {
            return Attempt::Reject(hnew);
        }
        // The stiffness test (accepted steps only).
        let num: f64 = (0..nn).map(|i| (k[6][i] - k[5][i]).powi(2)).sum();
        let den: f64 = (0..nn).map(|i| (ys_all[6][i] - ys_all[5][i]).powi(2)).sum();
        if den > 0.0 && h * (num / den).sqrt() > 3.25 {
            p.calm = 0;
            p.stiffs += 1;
            if p.stiffs >= 15 {
                p.stiff = true;
                p.jac = None;
                p.rho = 0.0;
            }
        } else {
            p.calm += 1;
            if p.calm >= 6 {
                p.stiffs = 0;
            }
        }
        let knot = evals.pop().expect("its last stage");
        Attempt::Accept { delta, knot, hnew, cont: None }
    }

    /// One step of Radau IIA of order 5 from the last knot (RADAU5): its
    /// Jacobian found where it is not kept, simplified Newton iterations on
    /// the stages from the last step's collocation polynomial, the embedded
    /// estimate, and the knot's state evaluated.
    pub(super) fn radau(&self, p: &mut Progress, t: f64, h: f64, t1: f64, across: bool) -> Attempt {
        let (model, lay) = (&self.model, &self.lay);
        let n = lay.n();
        let nn = 2 * n;
        if p.jac.is_none() {
            match self.jacobian(p, t) {
                Ok(j) => {
                    if p.rho == 0.0 {
                        p.rho = radius(&j);
                    }
                    p.jac = Some(j)
                }
                Err((at, ys, why)) => return Attempt::Failed(at, ys, why),
            }
        }
        let jac = p.jac.clone().expect("a Jacobian");
        let scal = self.scales(p, h);
        let rms = |v: &[f64]| (v.iter().zip(scal.iter().cycle()).map(|(x, s)| (x / s) * (x / s)).sum::<f64>() / v.len() as f64).sqrt();
        // E = I - h A (x) J over the three stages.
        let mut e = vec![vec![0.0; 3 * nn]; 3 * nn];
        for s in 0..3 {
            for q in 0..3 {
                for i in 0..nn {
                    for j in 0..nn {
                        e[s * nn + i][q * nn + j] = -h * A[s][q] * jac[i][j];
                    }
                }
            }
        }
        for (i, row) in e.iter_mut().enumerate() {
            row[i] += 1.0;
        }
        let Some(e) = Lu::of(e) else {
            p.jac = None;
            return Attempt::Reject(h * 0.5);
        };
        // Starting values: the last step's collocation polynomial carried
        // on, or none; Newton's rate trusted from the last step where it
        // starts from that polynomial.
        let mut z: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; nn]);
        if let Some([k1, k2, k3]) = &p.cont {
            let (c1m1, c2m1) = (C[0] - 1.0, C[1] - 1.0);
            let c3q = h / p.hold;
            for s in 0..3 {
                let cq = C[s] * c3q;
                for i in 0..nn {
                    z[s][i] = cq * (k1[i] + (cq - c2m1) * (k2[i] + (cq - c1m1) * k3[i]));
                }
            }
        }
        let mut faccon = if p.cont.is_some() { p.faccon.max(f64::EPSILON).powf(0.8) } else { 1.0 };
        let mut theta = p.theta;
        let mut stages: [Option<Eval>; 3] = [None, None, None];
        let (mut newt, mut converged, mut dynold) = (0, false, 0.0);
        while newt < NEWTON {
            newt += 1;
            let mut fz: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
            for s in 0..3 {
                let ys: Vec<f64> = (0..nn).map(|i| p.y[i] + z[s][i]).collect();
                p.traj.evaluations += 1;
                let at = t + C[s] * h;
                let guess = stages[s].as_ref().map_or(&p.last.acc, |x| &x.acc).clone();
                match forces(model, lay, at, &ys[..n], &ys[n..], &guess) {
                    Ok(ev) if let Some(r) = ev.inertia.iter().zip(&p.inertia).map(|(v, start)| v / start).find(|&r| r > FOLD) => return Attempt::Failed(at, ys, folded(at, r)),
                    Ok(ev) => {
                        fz[s] = deriv(&ys, &ev.acc);
                        stages[s] = Some(ev);
                    }
                    Err(f) => {
                        let why = said(model, at, &ys[n..], f);
                        return Attempt::Failed(at, ys, why);
                    }
                }
            }
            // Newton's step on G(Z) = -Z + h (A (x) I) F(Z).
            let g: Vec<f64> = (0..3).flat_map(|s| (0..nn).map(|i| -z[s][i] + h * (0..3).map(|q| A[s][q] * fz[q][i]).sum::<f64>()).collect::<Vec<_>>()).collect();
            let dz = e.solve(&g);
            let dyno = rms(&dz);
            if newt > 1 {
                theta = dyno / dynold;
                if theta >= 0.99 {
                    break;
                }
                faccon = theta / (1.0 - theta);
            }
            for s in 0..3 {
                for i in 0..nn {
                    z[s][i] += dz[s * nn + i];
                }
            }
            if faccon * dyno <= KAPPA {
                converged = true;
                break;
            }
            dynold = dyno.max(f64::EPSILON);
        }
        if !converged {
            // Newton diverging or slow: a shorter step, a fresh Jacobian.
            p.jac = None;
            p.faccon = 1.0;
            return Attempt::Reject(h * 0.5);
        }
        let stages: Vec<Eval> = stages.into_iter().flatten().collect();
        if let Some(r) = self.grown(p, &stages, h, across) {
            return r;
        }
        p.faccon = faccon;
        p.theta = theta;
        // The embedded estimate: (U1 / h - J)^-1 (f0 + sum DD_i Z_i / h),
        // again from f at y + it where the first is past the tolerance on a
        // first or rejected step (RADAU5's ESTRAD).
        let f2: Vec<f64> = (0..nn).map(|i| (DD[0] * z[0][i] + DD[1] * z[1][i] + DD[2] * z[2][i]) / h).collect();
        let mut e1 = vec![vec![0.0; nn]; nn];
        for i in 0..nn {
            for j in 0..nn {
                e1[i][j] = -jac[i][j];
            }
            e1[i][i] += U1 / h;
        }
        let Some(e1) = Lu::of(e1) else {
            p.jac = None;
            return Attempt::Reject(h * 0.5);
        };
        let cont: Vec<f64> = (0..nn).map(|i| p.f0[i] + f2[i]).collect();
        let est = e1.solve(&cont);
        let mut err = rms(&est).max(1e-10);
        if err >= 1.0 && (p.first || p.reject) && !across {
            let yp: Vec<f64> = (0..nn).map(|i| p.y[i] + est[i]).collect();
            p.traj.evaluations += 1;
            if let Ok(ev) = forces(model, lay, t, &yp[..n], &yp[n..], &p.last.acc) {
                let fp = deriv(&yp, &ev.acc);
                let cont: Vec<f64> = (0..nn).map(|i| fp[i] + f2[i]).collect();
                err = rms(&e1.solve(&cont)).max(1e-10);
            }
        }
        let fac = (0.9 * (2 * NEWTON + 1) as f64 / (newt + 2 * NEWTON) as f64).min(0.9);
        let hnew = h / (err.powf(0.25) / fac).clamp(1.0 / 8.0, 5.0);
        if !(err < 1.0 || across) {
            return Attempt::Reject(if p.first { h * 0.1 } else { hnew });
        }
        // The knot: its accelerations at the new state.
        let y1: Vec<f64> = (0..nn).map(|i| p.y[i] + z[2][i]).collect();
        let guess = stages.last().map_or(&p.last.acc, |x| &x.acc).clone();
        p.traj.evaluations += 1;
        let knot = match forces(model, lay, t1, &y1[..n], &y1[n..], &guess) {
            Ok(ev) => ev,
            Err(f) => {
                let why = said(model, t1, &y1[n..], f);
                return Attempt::Failed(t1, y1, why);
            }
        };
        // The Jacobian again at the next knot where Newton was slow, or
        // where it costs no state (no simulated coordinate: the springs'
        // own).
        if lay.ns == 0 || theta > THET {
            p.jac = None;
        }
        // The collocation polynomial's divided differences, from which the
        // next step starts (RADAU5's CONT).
        let (c1m1, c2m1, c1mc2) = (C[0] - 1.0, C[1] - 1.0, C[0] - C[1]);
        let mut k: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; nn]);
        for i in 0..nn {
            k[0][i] = (z[1][i] - z[2][i]) / c2m1;
            let ak = (z[0][i] - z[1][i]) / c1mc2;
            let acont3 = (ak - z[0][i] / C[0]) / C[1];
            k[1][i] = (ak - k[0][i]) / c1m1;
            k[2][i] = k[1][i] - acont3;
        }
        Attempt::Accept { delta: std::mem::take(&mut z[2]), knot, hnew, cont: (!across).then_some(k) }
    }
}
