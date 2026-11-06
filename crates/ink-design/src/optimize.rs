//! Optimization for design parameters, on parameters scaled to [0, 1].
//!
//! - `minimize`: bounded Nelder-Mead with constraints folded in as a penalty.
//!   Robust for the handful of parameters a design exposes and for noisy,
//!   non-smooth metrics (peak stress, grid effects).
//! - `minimize_mma`: the method of moving asymptotes (Svanberg; ink-mma),
//!   for smooth problems with exact gradients: many parameters in few
//!   iterations.
//!
//! Each evaluation may be an expensive solve, so budgets count evaluations.

/// Result of evaluating one candidate: objective (to minimize) and the
/// violation of each constraint (0 when satisfied, scaled to be comparable).
#[derive(Clone, Debug)]
pub struct Eval {
    pub objective: f64,
    pub violations: Vec<f64>,
}

impl Eval {
    pub fn feasible(&self) -> bool {
        self.violations.iter().all(|v| *v <= 0.0)
    }
}

pub struct Trial {
    pub x: Vec<f64>,
    pub eval: Eval,
    pub merit: f64,
}

pub struct Outcome {
    pub best: Trial,
    pub history: Vec<Trial>,
}

/// Minimize over `lo..hi` boxes. `f` may fail (e.g. an empty shape): such
/// points are treated as infeasible.
pub fn minimize(lo: &[f64], hi: &[f64], x0: &[f64], budget: usize, mut f: impl FnMut(&[f64]) -> Result<Eval, String>) -> Outcome {
    let n = lo.len();
    let to_x = |u: &[f64]| -> Vec<f64> { (0..n).map(|i| lo[i] + (hi[i] - lo[i]) * u[i].clamp(0.0, 1.0)).collect() };
    let mut history: Vec<Trial> = Vec::new();
    // Penalty weight: large relative to the objective's scale, set after the
    // first evaluation so units don't matter.
    let mut weight = 0.0f64;
    let mut eval = |u: &[f64], history: &mut Vec<Trial>, weight: &mut f64| -> f64 {
        let x = to_x(u);
        let e = f(&x).unwrap_or(Eval { objective: f64::INFINITY, violations: vec![f64::INFINITY] });
        if *weight == 0.0 {
            *weight = 1e3 * e.objective.abs().max(1e-9);
        }
        let merit = e.objective + *weight * e.violations.iter().map(|v| v.max(0.0)).sum::<f64>();
        let merit = if merit.is_finite() { merit } else { f64::MAX };
        history.push(Trial { x, eval: e, merit });
        merit
    };
    let u0: Vec<f64> = (0..n).map(|i| if hi[i] > lo[i] { (x0[i] - lo[i]) / (hi[i] - lo[i]) } else { 0.0 }).collect();
    // Initial simplex: x0 plus a step of 0.25 along each axis (inside the box).
    let mut simplex: Vec<(Vec<f64>, f64)> = Vec::with_capacity(n + 1);
    let m0 = eval(&u0, &mut history, &mut weight);
    simplex.push((u0.clone(), m0));
    for i in 0..n {
        let mut u = u0.clone();
        u[i] = if u[i] + 0.25 <= 1.0 { u[i] + 0.25 } else { u[i] - 0.25 };
        let m = eval(&u, &mut history, &mut weight);
        simplex.push((u, m));
    }
    let clamp = |u: Vec<f64>| -> Vec<f64> { u.into_iter().map(|v| v.clamp(0.0, 1.0)).collect() };
    while history.len() < budget {
        simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
        // Converged: the simplex has collapsed.
        let size = simplex.iter().skip(1).map(|(u, _)| u.iter().zip(&simplex[0].0).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max)).fold(0.0, f64::max);
        if size < 1e-3 {
            break;
        }
        let worst = simplex[n].clone();
        let centroid: Vec<f64> = (0..n).map(|i| simplex[..n].iter().map(|(u, _)| u[i]).sum::<f64>() / n as f64).collect();
        let along = |t: f64| clamp((0..n).map(|i| centroid[i] + t * (worst.0[i] - centroid[i])).collect());
        let refl = along(-1.0);
        let mr = eval(&refl, &mut history, &mut weight);
        if mr < simplex[0].1 {
            let exp = along(-2.0);
            let me = eval(&exp, &mut history, &mut weight);
            simplex[n] = if me < mr { (exp, me) } else { (refl, mr) };
        } else if mr < simplex[n - 1].1 {
            simplex[n] = (refl, mr);
        } else {
            let con = if mr < worst.1 { along(-0.5) } else { along(0.5) };
            let mc = eval(&con, &mut history, &mut weight);
            if mc < worst.1.min(mr) {
                simplex[n] = (con, mc);
            } else {
                // Shrink toward the best vertex.
                let best = simplex[0].0.clone();
                for v in simplex.iter_mut().skip(1) {
                    let u: Vec<f64> = (0..n).map(|i| best[i] + 0.5 * (v.0[i] - best[i])).collect();
                    let m = eval(&u, &mut history, &mut weight);
                    *v = (u, m);
                }
            }
        }
    }
    // Best feasible trial if any, else the least-violating one.
    let pick = history
        .iter()
        .enumerate()
        .filter(|(_, t)| t.eval.feasible())
        .min_by(|a, b| a.1.eval.objective.total_cmp(&b.1.eval.objective))
        .or_else(|| history.iter().enumerate().min_by(|a, b| a.1.merit.total_cmp(&b.1.merit)))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let best = Trial { x: history[pick].x.clone(), eval: history[pick].eval.clone(), merit: history[pick].merit };
    Outcome { best, history }
}

/// A point with gradients: objective, constraints g_i <= 0, and their
/// derivatives with respect to each parameter.
#[derive(Clone, Debug)]
pub struct Smooth {
    pub objective: f64,
    pub gradient: Vec<f64>,
    pub constraints: Vec<f64>,
    pub jacobian: Vec<Vec<f64>>,
}

/// Minimize over the box `lo..hi` with the method of moving asymptotes
/// (ink-mma), on the parameters scaled to [0, 1]: each iterate moves at
/// most 0.3 of a parameter's range. Stops when the design settles or after
/// `budget` evaluations. Returns the best feasible point seen (else the
/// least violating), with the objective history. A constraint up to 1e-6
/// counts as met: a caller that needs <= 0 passes its constraints shifted
/// by that.
pub fn minimize_mma(
    lo: &[f64],
    hi: &[f64],
    x0: &[f64],
    budget: usize,
    mut f: impl FnMut(&[f64]) -> Result<Smooth, String>,
) -> Result<(Vec<f64>, Smooth, Vec<(Vec<f64>, Smooth)>), String> {
    let n = lo.len();
    let span: Vec<f64> = (0..n).map(|i| (hi[i] - lo[i]).max(1e-300)).collect();
    // Work in u in [0, 1].
    let to_x = |u: &[f64]| -> Vec<f64> { (0..n).map(|i| lo[i] + span[i] * u[i]).collect() };
    let mut u: Vec<f64> = (0..n).map(|i| ((x0[i] - lo[i]) / span[i]).clamp(0.0, 1.0)).collect();
    let mut mma = ink_mma::Mma::new(0.3);
    let mut history: Vec<(Vec<f64>, Smooth)> = Vec::new();
    let mut scale = 0.0f64;
    for it in 0..budget.max(1) {
        let x = to_x(&u);
        let e = f(&x)?;
        // Gradients with respect to u; objective scaled to order one.
        if scale == 0.0 {
            scale = 1.0 / e.objective.abs().max(1e-12);
        }
        let g0: Vec<f64> = (0..n).map(|i| e.gradient[i] * span[i] * scale).collect();
        let constraints: Vec<(f64, Vec<f64>)> = e.constraints.iter().zip(&e.jacobian).map(|(c, g)| (*c, (0..n).map(|i| g[i] * span[i]).collect())).collect();
        history.push((x, e.clone()));
        let next = mma.update(&u, &g0, &constraints, None);
        let change = next.iter().zip(&u).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        u = next;
        if change < 1e-3 && e.constraints.iter().all(|c| *c <= 1e-6) && it > 0 {
            break;
        }
    }
    let best = history
        .iter()
        .filter(|(_, e)| e.constraints.iter().all(|c| *c <= 1e-6))
        .min_by(|a, b| a.1.objective.total_cmp(&b.1.objective))
        .or_else(|| history.iter().min_by(|a, b| a.1.constraints.iter().fold(0.0, |m: f64, c| m.max(*c)).total_cmp(&b.1.constraints.iter().fold(0.0, |m: f64, c| m.max(*c)))))
        .cloned()
        .ok_or("no evaluations")?;
    Ok((best.0, best.1, history))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_a_constrained_minimum() {
        // minimize (x-3)^2 + (y+1)^2 subject to x + y >= 4, in [0,10] x [-5,5]
        // -> optimum on x + y = 4: (4, 0) with value 2.
        let out = minimize(&[0.0, -5.0], &[10.0, 5.0], &[1.0, 1.0], 400, |x| {
            Ok(Eval { objective: (x[0] - 3.0).powi(2) + (x[1] + 1.0).powi(2), violations: vec![(4.0 - x[0] - x[1]) / 4.0] })
        });
        let b = &out.best;
        assert!(b.eval.feasible());
        assert!((b.eval.objective - 2.0).abs() < 0.05, "{:?} at {:?}", b.eval, b.x);
    }

    #[test]
    fn mma_solves_a_constrained_problem_with_many_variables() {
        // Minimize sum (x_i - 3)^2 over [0, 10]^20 subject to sum x_i <= 20:
        // optimum x_i = 1, objective 80.
        let n = 20;
        let (x, e, hist) = minimize_mma(&vec![0.0; n], &vec![10.0; n], &vec![5.0; n], 60, |x| {
            Ok(Smooth {
                objective: x.iter().map(|v| (v - 3.0).powi(2)).sum(),
                gradient: x.iter().map(|v| 2.0 * (v - 3.0)).collect(),
                constraints: vec![(x.iter().sum::<f64>() - 20.0) / 20.0],
                jacobian: vec![vec![1.0 / 20.0; x.len()]],
            })
        })
        .unwrap();
        assert!((e.objective - 80.0).abs() < 0.5 && x.iter().all(|v| (v - 1.0).abs() < 0.1), "{} at {x:?} after {}", e.objective, hist.len());
    }
}
