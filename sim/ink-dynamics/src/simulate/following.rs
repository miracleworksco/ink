use super::Simulation;
use super::linalg::{Lu, modes_of};
use super::simulation::Progress;
use super::types::{Following, Layout, Ringing, Span, Springs};
use crate::model::springs::At;

impl Simulation {
    /// What its springs did over the span (stepped to its end, or where it
    /// stops): each released joint's following error (Springs).
    pub fn springs(&self) -> Springs {
        let model = &self.model;
        let joints = model.released();
        let lines = model.spring_lines();
        let mut p = self.progress.lock().unwrap_or_else(|e| e.into_inner());
        self.step(&mut p, self.span.1);
        let modal = self.modal(&p);
        let traj = &p.traj;
        {
            let end = *traj.t.last().expect("a knot");
            // Where the written motion comes to rest for good: past the last
            // knot where a written coordinate moves or accelerates.
            let moving = |t: f64| {
                model
                    .written(t)
                    .is_ok_and(|w| w.iter().enumerate().any(|(k, q)| model.coords[k].moves() && !model.simulated.iter().any(|s| s.0 == k) && (q[1] != 0.0 || q[2] != 0.0)))
            };
            let still = match traj.t.iter().rposition(|&t| moving(t)) {
                None => Some(traj.t[0]),
                Some(k) if k + 1 < traj.t.len() => Some(traj.t[k + 1]),
                Some(_) => None,
            };
            // Each stepped deviation's ringing's centre (its static rest at
            // the end, a sag under gravity), every one removed before the
            // ringing is projected onto the modes: another's offset would
            // weigh on the modes it takes part in.
            let rested = self.rest_at(&p, end);
            let centres: Vec<f64> = (0..self.lay.w.len())
                .map(|w| {
                    let i = self.lay.ns + w;
                    if let Some(d) = &rested {
                        return d[self.lay.w[w]];
                    }
                    let Some(s) = still.filter(|&s| s < end) else { return 0.0 };
                    let at = traj.t.iter().copied().filter(|&t| t >= s).chain(traj.moments(i, (s, end)));
                    let (lo, hi) = at.map(|t| traj.at(i, t).map_or(0.0, |v| v[0])).fold((f64::INFINITY, f64::NEG_INFINITY), |m, v| (m.0.min(v), m.1.max(v)));
                    0.5 * (lo + hi)
                })
                .collect();
            let mut out = Vec::new();
            for (j, (path, turn)) in joints.iter().enumerate() {
                let Some(w) = self.lay.w.iter().position(|&x| x == j) else { continue };
                let i = self.lay.ns + w;
                let unit = if *turn { 180.0 / std::f64::consts::PI } else { 1.0 };
                let x = |t: f64| unit * traj.at(i, t).map_or(0.0, |v| v[0]);
                // Its knots and its turning points, where its extremes are.
                let turns = traj.moments(i, (traj.t[0], end));
                let mut at: Vec<f64> = traj.t.iter().copied().chain(turns.iter().copied()).collect();
                at.sort_by(f64::total_cmp);
                let peak = at.iter().map(|&t| (t, x(t))).fold((traj.t[0], 0.0f64), |m, (t, v)| if v.abs() > m.1.abs() { (t, v) } else { m });
                let ringing = still.filter(|&s| s < end).map(|s| {
                    let vals: Vec<f64> = at.iter().filter(|&&t| t >= s).map(|&t| x(t)).collect();
                    let (lo, hi) = vals.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |m, &v| (m.0.min(v), m.1.max(v)));
                    let (frequency, share) = modal.as_ref().filter(|_| traj.stopped.is_none()).and_then(|m| m(w, &centres)).unzip();
                    Ringing { from: s, amplitude: 0.5 * (hi - lo), frequency, share }
                });
                let k = model.released_stiffness(j);
                let effort = (k > 0.0).then(|| k * peak.1.abs() / unit);
                out.push(Following { path: path.to_string(), turn: *turn, peak: peak.1, at: peak.0, end: x(end), ringing, effort });
            }
            let held = self.lay.y.iter().map(|&j| joints[j].0.to_string()).collect();
            let spans = lines
                .iter()
                .enumerate()
                .map(|(l, (name, site))| {
                    let tl = traj.t.iter().zip(&traj.tensions).map(|(&t, ts)| (t, ts[l]));
                    let least = tl.clone().fold((traj.t[0], f64::INFINITY), |m, x| if x.1 < m.1 { x } else { m });
                    let most = tl.fold((traj.t[0], f64::NEG_INFINITY), |m, x| if x.1 > m.1 { x } else { m });
                    Span { name: name.to_string(), site: site.to_string(), least, most }
                })
                .collect();
            Springs { joints: out, held, spans, still }
        }
    }

    /// Where the springs hold the released joints at rest at moment t (the
    /// simulated coordinates where the last knot has them, every speed
    /// none): its ringing's centre, a sag under gravity. None where none is
    /// found.
    fn rest_at(&self, p: &Progress, t: f64) -> Option<Vec<f64>> {
        let (ns, nr) = (self.lay.ns, self.lay.nr);
        if nr == 0 {
            return None;
        }
        let sim: Vec<[f64; 3]> = (0..ns).map(|i| [p.y[i], 0.0, 0.0]).collect();
        let still: Vec<[f64; 3]> = self.model.beside(t, &sim).ok()?.into_iter().map(|q| [q[0], 0.0, 0.0]).collect();
        let m = self.model.state_near(t, &still);
        let st = m.springs.as_ref()?;
        rest(st, &m.fr, &m.noise_r).ok()
    }

    /// The springs' modes at the last knot (their stiffness along the
    /// stepped deviations, those without inertia condensed, and their
    /// inertia; the simulated coordinates as they move), and the state's
    /// ringing about its centres projected onto them: for stepped deviation
    /// w and every one's ringing's centre, the frequency of the mode
    /// carrying most of w's ringing there (Hz) and its share. None where
    /// nothing is released.
    fn modal(&self, p: &Progress) -> Option<Box<dyn Fn(usize, &[f64]) -> Option<(f64, f64)> + '_>> {
        let (ns, nw) = (self.lay.ns, self.lay.w.len());
        if nw == 0 || !p.last.two_sided {
            return None;
        }
        let n = self.lay.n();
        let m: Vec<Vec<f64>> = (0..nw).map(|a| (0..nw).map(|b| p.last.mss[ns + a][ns + b]).collect()).collect();
        let (lam, shapes) = modes_of(&m, &p.last.keff)?;
        let (y, v) = (p.y[ns..n].to_vec(), p.y[n + ns..2 * n].to_vec());
        Some(Box::new(move |w: usize, centres: &[f64]| {
            // The state about the ringings' centres.
            let d: Vec<f64> = y.iter().zip(centres).map(|(a, c)| a - c).collect();
            let mut best: Option<(f64, f64)> = None;
            let mut total = 0.0;
            for q in 0..nw {
                let omega = (1e6 * lam[q]).max(0.0).sqrt();
                let md: f64 = (0..nw).map(|i| (0..nw).map(|j| shapes[q][i] * m[i][j] * d[j]).sum::<f64>()).sum();
                let mv: f64 = (0..nw).map(|i| (0..nw).map(|j| shapes[q][i] * m[i][j] * v[j]).sum::<f64>()).sum();
                let amplitude = if omega > 0.0 { md.hypot(mv / omega) } else { md.abs() };
                let part = (amplitude * shapes[q][w]).abs();
                total += part;
                if best.is_none_or(|b| part > b.1) {
                    best = Some((omega / (2.0 * std::f64::consts::PI), part));
                }
            }
            best.filter(|_| total > 0.0).map(|(f, part)| (f, part / total))
        }))
    }
}

/// Where the springs hold the released joints at rest: their deviations
/// where the springs' pull along each balances what the parts need there
/// (`fr`: gravity, the loads), by Newton's method on their piecewise linear
/// pull (a span taut or slack), their stiffness's pseudo-inverse leaving a
/// direction nothing holds where it is. Err: a joint where their pull does
/// not balance what is needed past the rounding (`noise`): a span gone
/// slack under the load, a joint nothing holds against one, or the slack
/// pattern found no rest in 16 passes.
pub(super) fn rest(st: &At, fr: &[f64], noise: &[f64]) -> Result<Vec<f64>, usize> {
    let nr = fr.len();
    let mut d = vec![0.0; nr];
    let taut = |d: &[f64]| -> Vec<bool> { st.tensions(d).iter().map(|t| *t > 0.0).collect() };
    for _ in 0..16 {
        let fs = st.force(&d);
        let rr: Vec<f64> = (0..nr).map(|j| fs[j] - fr[j]).collect();
        let stiffness = st.stiffness(&d);
        let scale: Vec<_> = (0..nr).map(|j| if stiffness[j][j] > 0.0 { stiffness[j][j].sqrt() } else { 1.0 }).collect();
        let scaled: Vec<_> = stiffness.iter().zip(&scale).map(|(row, si)| row.iter().zip(&scale).map(|(k, sj)| k / si / sj).collect()).collect();
        let kinv = crate::model::pinv(&scaled);
        let before = taut(&d);
        for (j, row) in kinv.iter().enumerate() {
            d[j] += row.iter().zip(&rr).zip(&scale).map(|((a, b), si)| a * (b / si)).sum::<f64>() / scale[j];
        }
        if taut(&d) == before {
            break;
        }
    }
    let (fs, tensions) = (st.force(&d), st.tensions(&d));
    for j in 0..nr {
        let pull: f64 = (0..st.kl.len()).map(|l| (tensions[l] + st.t0[l]) * st.g[l][j].abs()).sum::<f64>() + (st.kj[j] * d[j]).abs();
        if (fs[j] - fr[j]).abs() > 1e-6 * (fr[j].abs() + pull) + 16.0 * noise[j] {
            return Err(j);
        }
    }
    Ok(d)
}

/// Whether a spring holds a motion `phi` of the stepped deviations at
/// deviations `d`, those carrying no inertia following it (K_yy^-1 K_yw): a
/// line or a joint's stiffness stretching along it (springs::stretches),
/// each deviation scaled by `unit` (a stepped one's, the square root of its
/// inertia along itself; one carrying none, 1). Each spring is judged on its
/// own, so a stiff one elsewhere hides no soft one's hold.
pub(super) fn holds(st: &At, d: &[f64], lay: &Layout, unit: &[f64], phi: &[f64]) -> bool {
    let mut z = vec![0.0; lay.nr];
    let mut scale = vec![1.0; lay.nr];
    for (a, &j) in lay.w.iter().enumerate() {
        (z[j], scale[j]) = (phi[a], unit[a]);
    }
    let (springs, _) = st.springs(d);
    if !lay.y.is_empty() {
        let k = st.stiffness(d);
        let kyy: Vec<Vec<f64>> = lay.y.iter().map(|&a| lay.y.iter().map(|&b| k[a][b]).collect()).collect();
        let pull: Vec<f64> = lay.y.iter().map(|&a| lay.w.iter().zip(phi).map(|(&b, p)| k[a][b] * p).sum()).collect();
        // As the condensed stiffness takes them (condensed): none where K_yy is singular.
        if let Some(lu) = Lu::of(kyy) {
            for (&j, x) in lay.y.iter().zip(lu.solve(&pull)) {
                z[j] = -x;
            }
        }
    }
    let zs: Vec<f64> = (0..lay.nr).map(|j| z[j] * scale[j]).collect();
    springs.iter().any(|(g, _)| crate::springs::stretches(&(0..lay.nr).map(|j| g[j] / scale[j]).collect::<Vec<_>>(), &zs))
}
