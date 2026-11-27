use super::Sets;
use ink_core::interval::{I, mul_o, sub_o};

/// The minimum weighted-norm efforts solving Bᵀ τ = f. Independent
/// drives are required by balance whatever their ratings; the ratings
/// select a split only when another drive can carry the same load.
pub fn shares(bm: &[Vec<f64>], span: &[f64], nq: usize) -> Result<Vec<Vec<f64>>, String> {
    let pm = effort_split(bm, span, &vec![0.0; bm.len()], nq).0;
    if pm.iter().flatten().all(|v| v.is_finite()) { Ok(pm) } else { Err(UNCERTAIN_DRIVES.into()) }
}

pub(super) const UNCERTAIN_DRIVES: &str = "the drive directions are too close to separate reliably, and their dependence could not be established";

/// Fixed drives start at their specified effort. Deviate from it only
/// along directions no adjustable drive can supply. Crossed limits take
/// last priority; the model separately reports why they are unavailable.
pub(super) fn effort_split(bm: &[Vec<f64>], span: &[f64], fixed: &[f64], nq: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    let mut sets = Sets((0..nq).collect());
    let mut active = Vec::new();
    for (a, row) in bm.iter().enumerate() {
        let axes: Vec<_> = (0..nq).filter(|&k| row[k] != 0.0).collect();
        if let Some(&first) = axes.first() {
            sets.join(&axes);
            active.push((a, first));
        }
    }
    let mut blocks = vec![(Vec::new(), Vec::new()); nq];
    for (a, first) in active {
        blocks[sets.find(first)].0.push(a);
    }
    for k in 0..nq {
        blocks[sets.find(k)].1.push(k);
    }
    let mut pm = vec![vec![0.0; nq]; bm.len()];
    let mut offset = vec![0.0; bm.len()];
    'block: for (mut drives, axes) in blocks.into_iter().filter(|(ds, _)| !ds.is_empty()) {
        // A coordinate's units multiply its whole column by one factor.
        // Remove it before testing independence, then express each drive
        // in units of its own largest component (including a turn's deg/rad).
        let scale: Vec<_> = axes.iter().map(|&k| drives.iter().map(|&a| bm[a][k].abs()).fold(0.0, f64::max)).collect();
        let size = |a: usize| axes.iter().zip(&scale).map(|(&k, s)| (bm[a][k] / s).abs()).fold(0.0, f64::max);
        // Add strongest first so ratios in the redundant split are <= 1.
        let tier = |a: usize| {
            if span[a] > 0.0 {
                2
            } else if span[a] == 0.0 {
                1
            } else {
                0
            }
        };
        let strength = |a: usize| size(a).log2() + if span[a] > 0.0 { span[a].log2() } else { 0.0 };
        drives.sort_by(|&a, &b| tier(b).cmp(&tier(a)).then_with(|| strength(b).total_cmp(&strength(a))));
        let mut basis = Vec::new();
        let mut originals = Vec::new();
        let mut inverse: Vec<Vec<f64>> = Vec::new();
        for (i, &a) in drives.iter().enumerate() {
            let unit = size(a);
            let drive: Vec<_> = axes.iter().zip(&scale).map(|(&k, s)| bm[a][k] / s / unit).collect();
            let residual = perpendicular(&drive, &basis);
            let norm = residual.iter().copied().fold(0.0, f64::hypot);
            // d is the previous drives' efforts for one unit of the new
            // drive. Giving it x changes those efforts by -d x (Greville).
            let d: Vec<_> = inverse.iter().map(|p| dot_n(p, &drive)).collect();
            let row: Vec<_> = if basis.len() < axes.len() && norm > 64.0 * f64::EPSILON * (axes.len() as f64).sqrt() {
                let q: Vec<_> = residual.iter().map(|r| r / norm).collect();
                let row = q.iter().map(|x| x / norm).collect();
                basis.push(q);
                originals.push(bm[a].as_slice());
                row
            } else {
                // A small singular value can still be the only way of
                // supplying a load. Share a direction only if its redundancy
                // is established; an ambiguous group has no numerical answer.
                if basis.len() < axes.len() && !in_span(&bm[a], &originals) {
                    for &b in &drives {
                        for &k in &axes {
                            pm[b][k] = f64::NAN;
                        }
                    }
                    continue 'block;
                }
                // Minimise x²/w² + Σ ((old_i - d_i x)/w_i)².
                let ratios: Vec<_> = drives[..i].iter().map(|&b| if tier(a) == tier(b) { (strength(a) - strength(b)).exp2() } else { 0.0 }).collect();
                let denom = 1.0 + d.iter().zip(&ratios).map(|(d, r)| (d * r).powi(2)).sum::<f64>();
                offset[a] = (fixed[a] * unit
                    + (0..i)
                        .map(|j| {
                            let b = drives[j];
                            d[j] * ratios[j] * ratios[j] * (offset[b] - fixed[b] * size(b))
                        })
                        .sum::<f64>())
                    / denom;
                (0..axes.len()).map(|k| (0..i).map(|j| d[j] * ratios[j] * ratios[j] / denom * inverse[j][k]).sum()).collect()
            };
            for ((p, &b), d) in inverse.iter_mut().zip(&drives).zip(d) {
                offset[b] -= d * offset[a];
                for (x, y) in p.iter_mut().zip(&row) {
                    *x -= d * y;
                }
            }
            inverse.push(row);
        }
        // On an incomplete span retain the ordinary least-squares projection
        // in the original coordinates (also used for reflected inertia).
        let mut physical = Vec::new();
        if basis.len() < axes.len() {
            for q in &basis {
                let v: Vec<_> = q.iter().zip(&scale).map(|(x, s)| x * s).collect();
                let v = perpendicular(&v, &physical);
                let norm = v.iter().copied().fold(0.0, f64::hypot);
                physical.push(v.iter().map(|x| x / norm).collect());
            }
        }
        for (&a, row) in drives.iter().zip(inverse) {
            let unit = size(a);
            offset[a] /= unit;
            let row: Vec<_> = row.iter().zip(&scale).map(|(x, s)| x / unit / s).collect();
            let projection: Vec<_> = physical.iter().map(|q| dot_n(q, &row)).collect();
            for (j, &k) in axes.iter().enumerate() {
                pm[a][k] = if physical.is_empty() { row[j] } else { physical.iter().zip(&projection).map(|(q, c)| q[j] * c).sum() };
            }
        }
    }
    (pm, offset)
}

/// Certify a suspected dependence of the *original*, unweighted drives.
/// Fraction-free elimination with outward-rounded intervals: only an
/// exact zero remainder establishes that adding this row cannot add rank.
fn in_span(row: &[f64], basis: &[&[f64]]) -> bool {
    if basis.iter().any(|b| *b == row || b.iter().zip(row).all(|(a, b)| *a == -*b)) {
        return true;
    }
    let mut rows: Vec<Vec<_>> = basis.iter().copied().chain([row]).map(|b| b.iter().copied().map(I::point).collect()).collect();
    for i in 0..basis.len() {
        let (before, after) = rows.split_at_mut(i + 1);
        let pivot = &before[i];
        let Some(k) = (0..row.len()).filter(|&k| pivot[k].lo > 0.0 || pivot[k].hi < 0.0).max_by(|&k, &l| pivot[k].lo.abs().total_cmp(&pivot[l].lo.abs())) else { return false };
        for v in after {
            let coeff = v[k];
            for j in 0..v.len() {
                v[j] = if j == k { I::point(0.0) } else { sub_o(mul_o(pivot[k], v[j]), mul_o(coeff, pivot[j])) };
            }
        }
    }
    rows.last().unwrap().iter().all(|v| v.lo == 0.0 && v.hi == 0.0)
}

pub(super) fn dot_n(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Twice-modified Gram–Schmidt, retaining the part outside the basis
/// even when the input is nearly along its first direction.
fn perpendicular(v: &[f64], basis: &[Vec<f64>]) -> Vec<f64> {
    let mut r = v.to_vec();
    for _ in 0..2 {
        for q in basis {
            let d = dot_n(&r, q);
            for (x, y) in r.iter_mut().zip(q) {
                *x -= d * y;
            }
        }
    }
    r
}

/// The pseudo-inverse of a symmetric positive semidefinite matrix:
/// eigenvalues below 1e-12 of the largest (rounding) are taken as 0.
pub(crate) fn pinv(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let (e, v) = eigen(a);
    let largest = e.iter().map(|x| x.abs()).fold(0.0, f64::max);
    let keep: Vec<f64> = e.iter().map(|&x| if x > 1e-12 * largest { 1.0 / x } else { 0.0 }).collect();
    (0..n).map(|i| (0..n).map(|j| (0..n).map(|k| v[i][k] * keep[k] * v[j][k]).sum()).collect()).collect()
}

/// A symmetric matrix's eigenvalues and eigenvectors (the columns of the
/// second), by Jacobi's method.
pub(crate) fn eigen(a: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut m: Vec<Vec<f64>> = a.to_vec();
    let mut v: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect();
    for _ in 0..100 {
        let off: f64 = (0..n).flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j))).map(|(i, j)| m[i][j] * m[i][j]).sum();
        if off == 0.0 {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if m[p][q] == 0.0 {
                    continue;
                }
                let theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let (c, s) = (1.0 / (t * t + 1.0).sqrt(), t / (t * t + 1.0).sqrt());
                for k in 0..n {
                    let (mkp, mkq) = (m[k][p], m[k][q]);
                    m[k][p] = c * mkp - s * mkq;
                    m[k][q] = s * mkp + c * mkq;
                }
                for k in 0..n {
                    let (mpk, mqk) = (m[p][k], m[q][k]);
                    m[p][k] = c * mpk - s * mqk;
                    m[q][k] = s * mpk + c * mqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k][p], v[k][q]);
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    ((0..n).map(|i| m[i][i]).collect(), v)
}
