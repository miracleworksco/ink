//! The method of moving asymptotes (Svanberg 1987, with the defaults of his
//! 2007 notes), for problems with exact gradients on variables scaled to
//! [0, 1]: minimize f(x) subject to g_j(x) <= 0. Each iteration fits a
//! separable convex approximation to every function at the iterate (terms
//! p / (U - x) + q / (x - L), from the positive and negative parts of its
//! gradient), within asymptotes L, U that widen while a variable moves
//! steadily and close in where it oscillates, and steps to the
//! approximation's minimum, found through its dual.
//!
//! Used by the design package's gradient optimizer (`optimize`), and open
//! to any package built on the kernel (topology optimization, say): no
//! package depends on another (docs/CORE.md), so it is a crate of its own,
//! plain numerics with no other dependency.

/// The method's state across iterations: the two previous iterates and the
/// asymptotes.
pub struct Mma {
    /// Most a variable moves in one iteration.
    step: f64,
    x1: Vec<f64>,
    x2: Vec<f64>,
    low: Vec<f64>,
    upp: Vec<f64>,
    k: usize,
}

impl Mma {
    /// A run whose variables each move at most `step` per iteration.
    pub fn new(step: f64) -> Mma {
        Mma { step, x1: Vec::new(), x2: Vec::new(), low: Vec::new(), upp: Vec::new(), k: 0 }
    }

    /// The next iterate from `x` (in [0, 1]), given the objective's
    /// gradient `df` and each constraint's value and gradient
    /// (`g_j(x) <= 0`). Gradients are best scaled to order one per
    /// variable, as the approximations' regularisation (1e-5) assumes.
    ///
    /// Each constraint's multiplier is the smallest at which the
    /// approximation's minimiser meets that constraint's approximation,
    /// found one at a time by bisection (the dual is concave and its slope
    /// along a multiplier falls as it grows), in sweeps until they settle.
    /// `meet`, with one constraint: whether a design breaks the constraint
    /// itself, met in place of its approximation (a volume budget met
    /// exactly at every iterate, as optimality criteria do; the
    /// approximation's lagged it as the projection sharpened).
    pub fn update(&mut self, x: &[f64], df: &[f64], constraints: &[(f64, Vec<f64>)], meet: Option<&dyn Fn(&[f64]) -> bool>) -> Vec<f64> {
        let n = x.len();
        // Asymptotes.
        let (init, incr, decr) = (0.5, 1.2, 0.7);
        if self.k < 2 {
            self.low = x.iter().map(|v| v - init).collect();
            self.upp = x.iter().map(|v| v + init).collect();
        } else {
            for j in 0..n {
                let osc = (x[j] - self.x1[j]) * (self.x1[j] - self.x2[j]);
                let gamma = if osc < 0.0 {
                    decr
                } else if osc > 0.0 {
                    incr
                } else {
                    1.0
                };
                self.low[j] = (x[j] - gamma * (self.x1[j] - self.low[j])).clamp(x[j] - 10.0, x[j] - 0.01);
                self.upp[j] = (x[j] + gamma * (self.upp[j] - self.x1[j])).clamp(x[j] + 0.01, x[j] + 10.0);
            }
        }
        let (low, upp) = (&self.low, &self.upp);
        // Move limits.
        let albefa = 0.1;
        let alfa: Vec<f64> = (0..n).map(|j| (low[j] + albefa * (x[j] - low[j])).max(x[j] - self.step).max(0.0)).collect();
        let beta: Vec<f64> = (0..n).map(|j| (upp[j] - albefa * (upp[j] - x[j])).min(x[j] + self.step).min(1.0)).collect();
        // p / (U - x) + q / (x - L) per variable for each function.
        let pq = |d: &[f64]| -> Vec<(f64, f64)> {
            (0..n)
                .map(|j| {
                    let (p, q) = (d[j].max(0.0), (-d[j]).max(0.0));
                    let r = 0.001 * (p + q) + 1e-5;
                    ((p + r) * (upp[j] - x[j]).powi(2), (q + r) * (x[j] - low[j]).powi(2))
                })
                .collect()
        };
        let obj = pq(df);
        let cons: Vec<Vec<(f64, f64)>> = constraints.iter().map(|(_, d)| pq(d)).collect();
        // The minimiser of objective + sum lam_j constraint_j's approximations.
        let at = |lam: &[f64]| -> Vec<f64> {
            (0..n)
                .map(|j| {
                    let (mut p, mut q) = obj[j];
                    for (l, c) in lam.iter().zip(&cons) {
                        (p, q) = (p + l * c[j].0, q + l * c[j].1);
                    }
                    let (sp, sq) = (p.max(1e-300).sqrt(), q.max(1e-300).sqrt());
                    ((sp * low[j] + sq * upp[j]) / (sp + sq)).clamp(alfa[j], beta[j])
                })
                .collect()
        };
        // Constraint i's approximation at y, from its value at x.
        let approx = |i: usize, y: &[f64]| -> f64 {
            let c = &cons[i];
            constraints[i].0 + (0..n).map(|j| c[j].0 * (1.0 / (upp[j] - y[j]) - 1.0 / (upp[j] - x[j])) + c[j].1 * (1.0 / (y[j] - low[j]) - 1.0 / (x[j] - low[j]))).sum::<f64>()
        };
        let m = constraints.len();
        let mut lam = vec![0.0; m];
        for _sweep in 0..if m > 1 { 30 } else { 1 } {
            for i in 0..m {
                // Whether the constraint is broken at multiplier l.
                let mut broken = |l: f64| {
                    lam[i] = l;
                    let y = at(&lam);
                    match meet {
                        Some(over) if m == 1 => over(&y),
                        _ => approx(i, &y) > 0.0,
                    }
                };
                if !broken(0.0) {
                    continue;
                }
                let (mut lo, mut hi) = (0.0, 1.0);
                while broken(hi) && hi < 1e30 {
                    (lo, hi) = (hi, hi * 10.0);
                }
                for _ in 0..200 {
                    let mid = if lo == 0.0 { 0.5 * hi } else { (lo * hi).sqrt() };
                    if broken(mid) {
                        lo = mid
                    } else {
                        hi = mid
                    }
                    if hi - lo <= 1e-9 * hi {
                        break;
                    }
                }
                lam[i] = hi;
            }
        }
        let next = at(&lam);
        self.x2 = std::mem::replace(&mut self.x1, x.to_vec());
        self.k += 1;
        next
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_a_constrained_problem_with_many_variables() {
        // Minimize sum (x_i - 0.3)^2 over [0, 1]^20 subject to
        // sum x_i <= 2: optimum x_i = 0.1.
        let mut mma = Mma::new(0.3);
        let mut x = vec![0.5; 20];
        for _ in 0..60 {
            let df: Vec<f64> = x.iter().map(|v| 2.0 * (v - 0.3)).collect();
            let g = (x.iter().sum::<f64>() - 2.0) / 2.0;
            x = mma.update(&x, &df, &[(g, vec![0.5; 20])], None);
        }
        assert!(x.iter().all(|v| (v - 0.1).abs() < 0.01), "{x:?}");
    }

    #[test]
    fn meets_a_constraint_itself_when_given_it() {
        // A linear objective falling with every variable, a budget on
        // their sum: met exactly at each iterate when the test is given.
        let mut mma = Mma::new(0.2);
        let mut x = vec![0.3; 10];
        for _ in 0..5 {
            let over = |y: &[f64]| y.iter().sum::<f64>() > 3.0;
            x = mma.update(&x, &[-1.0; 10], &[(x.iter().sum::<f64>() - 3.0, vec![1.0; 10])], Some(&over));
            assert!(x.iter().sum::<f64>() <= 3.0 && x.iter().sum::<f64>() > 3.0 - 1e-6, "{}", x.iter().sum::<f64>());
        }
    }
}
