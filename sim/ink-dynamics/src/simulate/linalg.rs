/// A square matrix's LU factors with partial pivoting (in place, and the
/// rows' order); None where it is singular.
pub(super) struct Lu {
    factors: Vec<Vec<f64>>,
    perm: Vec<usize>,
    rows: Vec<f64>,
    cols: Vec<f64>,
}

impl Lu {
    pub(super) fn of(mut a: Vec<Vec<f64>>) -> Option<Lu> {
        let n = a.len();
        let mut perm: Vec<usize> = (0..n).collect();
        // The Newton and collocation systems need not be symmetric.
        // Equilibrate both equations and unknowns before partial pivoting.
        let rows: Vec<_> = a.iter().map(|row| row.iter().map(|x| x.abs()).fold(0.0, f64::max)).collect();
        if rows.iter().any(|s| !s.is_finite() || *s == 0.0) {
            return None;
        }
        for (row, s) in a.iter_mut().zip(&rows) {
            row.iter_mut().for_each(|v| *v /= s);
        }
        let cols: Vec<_> = (0..n).map(|j| a.iter().map(|row| row[j].abs()).fold(0.0, f64::max)).collect();
        if cols.contains(&0.0) {
            return None;
        }
        for row in &mut a {
            row.iter_mut().zip(&cols).for_each(|(v, s)| *v /= s);
        }
        for c in 0..n {
            let p = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
            if !(a[p][c].abs() > 1e-14) {
                return None;
            }
            a.swap(c, p);
            perm.swap(c, p);
            for i in c + 1..n {
                let f = a[i][c] / a[c][c];
                a[i][c] = f;
                for k in c + 1..n {
                    a[i][k] -= f * a[c][k];
                }
            }
        }
        Some(Lu { factors: a, perm, rows, cols })
    }

    pub(super) fn solve(&self, b: &[f64]) -> Vec<f64> {
        let a = &self.factors;
        let n = b.len();
        let mut y: Vec<f64> = self.perm.iter().map(|&p| b[p] / self.rows[p]).collect();
        for i in 0..n {
            y[i] -= (0..i).map(|k| a[i][k] * y[k]).sum::<f64>();
        }
        for i in (0..n).rev() {
            y[i] = (y[i] - (i + 1..n).map(|k| a[i][k] * y[k]).sum::<f64>()) / a[i][i];
        }
        y.iter_mut().zip(&self.cols).for_each(|(v, s)| *v /= s);
        y
    }
}

/// The solution of a symmetric positive definite system, by Cholesky; None
/// where it is not (a direction of no inertia).
pub(super) fn solve(m: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    // Jacobi scaling: the diagonal is one in every coordinate, so
    // only dependence of the motions can make a pivot negligible.
    let scale: Vec<_> = (0..n).map(|i| m[i][i].sqrt()).collect();
    if scale.iter().any(|s| !s.is_finite() || *s == 0.0) {
        return None;
    }
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let s: f64 = m[i][j] / scale[i] / scale[j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if !(s > 1e-12) {
                    return None;
                }
                l[i][i] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    let mut y = vec![0.0; n];
    for i in 0..n {
        y[i] = (b[i] / scale[i] - (0..i).map(|k| l[i][k] * y[k]).sum::<f64>()) / l[i][i];
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        x[i] = (y[i] - ((i + 1)..n).map(|k| l[k][i] * x[k]).sum::<f64>()) / l[i][i];
    }
    x.iter_mut().zip(scale).for_each(|(v, s)| *v /= s);
    Some(x)
}

/// The springs' modes along the stepped deviations: K phi = lambda M phi
/// (lambda times 1e6 is omega^2, M in g mm^2 per unit^2, K the condensed
/// stiffness), by M = L L^T and Jacobi's method on L^-1 K L^-T; each
/// eigenvalue and its M-normal shape. None where M is not positive.
pub(super) fn modes_of(m: &[Vec<f64>], k: &[Vec<f64>]) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let nw = m.len();
    let mut l = vec![vec![0.0; nw]; nw];
    for i in 0..nw {
        for j in 0..=i {
            let s: f64 = m[i][j] - (0..j).map(|q| l[i][q] * l[j][q]).sum::<f64>();
            if i == j {
                if !(s > 0.0) {
                    return None;
                }
                l[i][i] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    let lower = |b: &[f64]| -> Vec<f64> {
        let mut y = vec![0.0; nw];
        for i in 0..nw {
            y[i] = (b[i] - (0..i).map(|q| l[i][q] * y[q]).sum::<f64>()) / l[i][i];
        }
        y
    };
    let upper = |b: &[f64]| -> Vec<f64> {
        let mut x = vec![0.0; nw];
        for i in (0..nw).rev() {
            x[i] = (b[i] - (i + 1..nw).map(|q| l[q][i] * x[q]).sum::<f64>()) / l[i][i];
        }
        x
    };
    let cols: Vec<Vec<f64>> = (0..nw)
        .map(|j| {
            let e: Vec<f64> = (0..nw).map(|i| f64::from(u8::from(i == j))).collect();
            let u = upper(&e);
            let ku: Vec<f64> = (0..nw).map(|i| (0..nw).map(|q| k[i][q] * u[q]).sum()).collect();
            lower(&ku)
        })
        .collect();
    let a: Vec<Vec<f64>> = (0..nw).map(|i| (0..nw).map(|j| 0.5 * (cols[j][i] + cols[i][j])).collect()).collect();
    let (lam, vecs) = crate::model::eigen(&a);
    let shapes = (0..nw).map(|q| upper(&(0..nw).map(|r| vecs[r][q]).collect::<Vec<_>>())).collect();
    Some((lam, shapes))
}

/// A square matrix's spectral radius, by the power method: the mean growth
/// of a vector over 40 products after 20 (a pair of imaginary eigenvalues,
/// an undamped mode's, alternates its growth: the mean over an even count).
pub(super) fn radius(j: &[Vec<f64>]) -> f64 {
    let n = j.len();
    let mut x: Vec<f64> = (0..n).map(|i| 1.0 + 0.1 * i as f64).collect();
    let mut log = 0.0;
    for k in 0..60 {
        let y: Vec<f64> = (0..n).map(|i| (0..n).map(|l| j[i][l] * x[l]).sum()).collect();
        let (a, b) = (x.iter().map(|v| v * v).sum::<f64>().sqrt(), y.iter().map(|v| v * v).sum::<f64>().sqrt());
        if !(b > 0.0 && a > 0.0) {
            return 0.0;
        }
        if k >= 20 {
            log += (b / a).ln();
        }
        x = y.iter().map(|v| v / b).collect();
    }
    (log / 40.0).exp()
}

/// e^X and phi1(X) = (e^X - I) X^-1 (its series, sum X^k / (k + 1)!, where
/// X is singular) of a square matrix X: each's Taylor series to the 16th
/// power of Y = X / 2^s, the least s putting Y's norm at most 1/2, then
/// doubled s times, e^2Y = (e^Y)^2 and phi1(2Y) = phi1(Y) (e^Y + I) / 2
/// (scaling and squaring: Higham, Functions of Matrices, 10.3).
pub(super) fn exp_phi(x: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = x.len();
    let norm = (0..n).map(|j| (0..n).map(|i| x[i][j].abs()).sum::<f64>()).fold(0.0, f64::max);
    if !norm.is_finite() {
        return (vec![vec![f64::NAN; n]; n], vec![vec![f64::NAN; n]; n]);
    }
    let s = if norm > 0.5 { (2.0 * norm).log2().ceil() as i32 } else { 0 };
    let y: Vec<Vec<f64>> = x.iter().map(|r| r.iter().map(|v| v / 2f64.powi(s)).collect()).collect();
    let mul = |a: &[Vec<f64>], b: &[Vec<f64>]| -> Vec<Vec<f64>> { (0..n).map(|i| (0..n).map(|j| (0..n).map(|k| a[i][k] * b[k][j]).sum()).collect()).collect() };
    // I + Y / k (I + Y / (k + 1) (... (I + Y / (k + 15)))).
    let series = |k: usize| {
        let mut acc: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| f64::from(u8::from(i == j))).collect()).collect();
        for d in (k..k + 16).rev() {
            acc = mul(&y, &acc);
            for (i, row) in acc.iter_mut().enumerate() {
                row.iter_mut().for_each(|v| *v /= d as f64);
                row[i] += 1.0;
            }
        }
        acc
    };
    let (mut e, mut phi) = (series(1), series(2));
    for _ in 0..s {
        let mut plus = e.clone();
        (0..n).for_each(|i| plus[i][i] += 1.0);
        phi = mul(&phi, &plus).into_iter().map(|r| r.into_iter().map(|v| 0.5 * v).collect()).collect();
        e = mul(&e, &e);
    }
    (e, phi)
}
