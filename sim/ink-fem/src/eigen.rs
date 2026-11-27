//! Eigenvalues: of a dense symmetric matrix (Jacobi), and the lowest few of
//! a large one (LOBPCG).

use crate::GRAIN;
use crate::numerics::dot;
use rayon::prelude::*;

/// Eigen-decomposition of a dense symmetric matrix (cyclic Jacobi):
/// eigenvalues ascending, eigenvectors as columns (row-major n x n).
pub(crate) fn dense_sym_eigen(a: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut a = a.to_vec();
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _sweep in 0..100 {
        let off: f64 = (0..n).flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j))).map(|(i, j)| a[i * n + j].powi(2)).sum();
        let diag: f64 = (0..n).map(|i| a[i * n + i].powi(2)).sum();
        if off <= 1e-30 * diag.max(1e-300) {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p * n + q].abs() < 1e-300 {
                    continue;
                }
                let theta = 0.5 * libm::atan2(2.0 * a[p * n + q], a[q * n + q] - a[p * n + p]);
                let (c, s) = (libm::cos(theta), libm::sin(theta));
                for k in 0..n {
                    let (akp, akq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[i * n + i].total_cmp(&a[j * n + j]));
    let vals = order.iter().map(|&i| a[i * n + i]).collect();
    let mut vecs = vec![0.0; n * n];
    for (c, &i) in order.iter().enumerate() {
        for r in 0..n {
            vecs[r * n + c] = v[r * n + i];
        }
    }
    (vals, vecs)
}

/// The lowest `block` eigenpairs of A x = nu B x (A symmetric, B
/// symmetric positive definite), by LOBPCG with the preconditioner `t` (an
/// approximate inverse of the stiffness: multigrid). Each step takes the
/// Rayleigh-Ritz pairs on the span of the block X, its preconditioned
/// residuals W and its last step P, kept with their products by A and B
/// (each new vector applies A and B once). Every vector is kept
/// B-orthogonal to `away` (B-orthonormal; A gives each next to no energy:
/// a free part's rigid motions), its share of them taken out twice, from
/// each start and residual and from the block and its step after each
/// combination, so rounding cannot bring them back. Converged where `error` is at
/// most 1 for each of the first `want` pairs (their residual over its
/// tolerance, from the value nu, the residual r = A x - nu B x, the
/// preconditioned residual t r, and A x); stopped, not converged, where
/// the largest has not halved in 50 steps. The values ascending, the
/// vectors B-orthonormal, the iterations, and whether it converged.
#[allow(clippy::type_complexity)]
pub(crate) fn lobpcg(
    n: usize,
    a: &dyn Fn(&[f64]) -> Vec<f64>,
    b: &dyn Fn(&[f64]) -> Vec<f64>,
    t: &dyn Fn(&[Vec<f64>]) -> Vec<Vec<f64>>,
    (block, want): (usize, usize),
    error: &dyn Fn(f64, &[f64], &[f64], &[f64]) -> f64,
    max_iter: usize,
    away: &[Vec<f64>],
) -> (Vec<f64>, Vec<Vec<f64>>, usize, bool) {
    // Each of `away` with its products, and a vector's share of them out of
    // it and its products alike.
    let (away_a, away_b): (Vec<Vec<f64>>, Vec<Vec<f64>>) = (away.iter().map(|y| a(y)).collect(), away.iter().map(|y| b(y)).collect());
    let clear = |v: &mut [f64], av: Option<&mut [f64]>, bv: Option<&mut [f64]>| {
        let mut taken = vec![0.0; away.len()];
        for _ in 0..2 {
            for (k, (y, by)) in away.iter().zip(&away_b).enumerate() {
                let c = dot(by, v);
                taken[k] += c;
                v.par_iter_mut().zip(y.par_iter()).for_each(|(v, y)| *v -= c * y);
            }
        }
        for (prod, ys) in [(av, &away_a), (bv, &away_b)] {
            if let Some(p) = prod {
                for (y, c) in ys.iter().zip(&taken) {
                    p.par_iter_mut().zip(y.par_iter()).for_each(|(p, y)| *p -= c * y);
                }
            }
        }
    };
    let t = |r: &[Vec<f64>]| -> Vec<Vec<f64>> {
        let mut z = t(r);
        z.iter_mut().for_each(|z| clear(z, None, None));
        z
    };
    let mut seed = 0x2545_F491_4F6C_DD1Du64;
    // Preconditioned random vectors: smooth, low-frequency content.
    let random: Vec<Vec<f64>> = (0..block)
        .map(|_| {
            (0..n)
                .map(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
                })
                .collect()
        })
        .collect();
    let start = t(&random);
    let products = |v: &[Vec<f64>]| -> (Vec<Vec<f64>>, Vec<Vec<f64>>) { (v.iter().map(|x| a(x)).collect(), v.iter().map(|x| b(x)).collect()) };
    let (sa, sb) = products(&start);
    // Rayleigh-Ritz on the span of `s` (with its products): the lowest
    // `block` pairs' coefficients, each a column over s.
    let ritz = |s: &[&Vec<f64>], sa: &[&Vec<f64>], sb: &[&Vec<f64>]| -> (Vec<f64>, Vec<Vec<f64>>) {
        let m = s.len();
        let gram = |l: &[&Vec<f64>]| -> Vec<f64> {
            let mut g = vec![0.0; m * m];
            for i in 0..m {
                for j in i..m {
                    let v = 0.5 * (dot(s[i], l[j]) + dot(s[j], l[i]));
                    g[i * m + j] = v;
                    g[j * m + i] = v;
                }
            }
            g
        };
        let (ga, gb) = (gram(sa), gram(sb));
        // B's Gram matrix, scaled to a unit diagonal, by its eigenpairs:
        // directions it holds under 1e-10 of its largest (dependent ones)
        // dropped, the rest scaled to B-orthonormal: Z.
        let d: Vec<f64> = (0..m).map(|i| 1.0 / gb[i * m + i].max(1e-300).sqrt()).collect();
        let scaled: Vec<f64> = (0..m * m).map(|ij| gb[ij] * d[ij / m] * d[ij % m]).collect();
        let (ev, evec) = dense_sym_eigen(&scaled, m);
        let top = ev.iter().cloned().fold(0.0, f64::max);
        let keep: Vec<usize> = (0..m).filter(|&k| ev[k] > 1e-10 * top).collect();
        let r = keep.len();
        let z: Vec<f64> = (0..m * r).map(|ik| d[ik / r] * evec[(ik / r) * m + keep[ik % r]] / ev[keep[ik % r]].sqrt()).collect();
        // Z^T G_A Z, and its eigenpairs.
        let az: Vec<f64> = (0..m * r).map(|ik| (0..m).map(|j| ga[(ik / r) * m + j] * z[j * r + ik % r]).sum()).collect();
        let h: Vec<f64> = (0..r * r).map(|kl| (0..m).map(|i| z[i * r + kl / r] * az[i * r + kl % r]).sum()).collect();
        let (vals, y) = dense_sym_eigen(&h, r);
        let k = block.min(r);
        let coefs = (0..k).map(|c| (0..m).map(|i| (0..r).map(|l| z[i * r + l] * y[l * r + c]).sum()).collect()).collect();
        (vals[..k].to_vec(), coefs)
    };
    let combine = |cols: &[&Vec<f64>], coef: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0; n];
        for (c, &w) in cols.iter().zip(coef) {
            if w != 0.0 {
                out.par_iter_mut().with_min_len(GRAIN).zip(c.par_iter()).for_each(|(o, v)| *o += w * v);
            }
        }
        out
    };
    let (s0, sa0, sb0): (Vec<&Vec<f64>>, Vec<&Vec<f64>>, Vec<&Vec<f64>>) = (start.iter().collect(), sa.iter().collect(), sb.iter().collect());
    let (mut lam, coefs) = ritz(&s0, &sa0, &sb0);
    let mut x: Vec<Vec<f64>> = coefs.iter().map(|c| combine(&s0, c)).collect();
    let mut ax: Vec<Vec<f64>> = coefs.iter().map(|c| combine(&sa0, c)).collect();
    let mut bx: Vec<Vec<f64>> = coefs.iter().map(|c| combine(&sb0, c)).collect();
    let (mut p, mut ap, mut bp): (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>) = (Vec::new(), Vec::new(), Vec::new());
    let mut it = 0;
    let mut converged = false;
    let (mut checkpoint, mut checkpoint_err) = (0, f64::INFINITY);
    while it < max_iter {
        it += 1;
        // Residuals, and the convergence of the wanted pairs.
        let r: Vec<Vec<f64>> = (0..x.len()).map(|j| ax[j].par_iter().zip(&bx[j]).map(|(a, b)| a - lam[j] * b).collect()).collect();
        let w = t(&r);
        let err = (0..want.min(x.len())).map(|j| error(lam[j], &r[j], &w[j], &ax[j])).fold(0.0, f64::max);
        converged = err <= 1.0;
        if converged {
            break;
        }
        if it - checkpoint >= 50 {
            if !(err <= 0.5 * checkpoint_err) {
                break;
            }
            (checkpoint, checkpoint_err) = (it, err);
        }
        let (aw, bw) = products(&w);
        let s: Vec<&Vec<f64>> = x.iter().chain(&w).chain(&p).collect();
        let sa: Vec<&Vec<f64>> = ax.iter().chain(&aw).chain(&ap).collect();
        let sb: Vec<&Vec<f64>> = bx.iter().chain(&bw).chain(&bp).collect();
        let (vals, coefs) = ritz(&s, &sa, &sb);
        // The new block, and its step: the part of it outside the old X.
        let nx = x.len();
        let step = |c: &[f64]| -> Vec<f64> { c.iter().enumerate().map(|(i, &v)| if i < nx { 0.0 } else { v }).collect() };
        let (mut x2, mut ax2, mut bx2, mut p2, mut ap2, mut bp2) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for c in &coefs {
            let (mut xc, mut axc, mut bxc) = (combine(&s, c), combine(&sa, c), combine(&sb, c));
            clear(&mut xc, Some(&mut axc), Some(&mut bxc));
            x2.push(xc);
            ax2.push(axc);
            bx2.push(bxc);
            let d = step(c);
            let (mut pc, mut apc, mut bpc) = (combine(&s, &d), combine(&sa, &d), combine(&sb, &d));
            clear(&mut pc, Some(&mut apc), Some(&mut bpc));
            p2.push(pc);
            ap2.push(apc);
            bp2.push(bpc);
        }
        (x, ax, bx, p, ap, bp, lam) = (x2, ax2, bx2, p2, ap2, bp2, vals);
    }
    (lam, x, it, converged)
}
