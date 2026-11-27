//! Material occupancy and principal axes for mass properties.

use ink_core::V3;
use ink_core::sign::inside;

/// Fraction of a box of sides `s` centred on a sample that lies inside,
/// the field taken as linear there: the slab fraction along the gradient.
pub fn occupancy(f: &[f64; 4], s: V3) -> f64 {
    let w = s[0] * f[1].abs() + s[1] * f[2].abs() + s[2] * f[3].abs(); // field range across the box
    if !(w > 0.0 && w.is_finite() && f[0].is_finite()) {
        return if inside(f[0]) { 1.0 } else { 0.0 };
    }
    (0.5 - f[0] / w).clamp(0.0, 1.0)
}

/// Eigenvalues and eigenvectors (columns) of a symmetric 3x3 matrix (Jacobi).
pub fn sym_eigen(m: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut a = m;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let (mut p, mut q, mut big) = (0, 1, 0.0);
        for i in 0..3 {
            for j in i + 1..3 {
                if a[i][j].abs() > big {
                    (p, q, big) = (i, j, a[i][j].abs());
                }
            }
        }
        if big < 1e-15 * (a[0][0].abs() + a[1][1].abs() + a[2][2].abs()).max(1e-300) {
            break;
        }
        let theta = 0.5 * libm::atan2(2.0 * a[p][q], a[q][q] - a[p][p]);
        let (c, s) = (libm::cos(theta), libm::sin(theta));
        for k in 0..3 {
            let (akp, akq) = (a[k][p], a[k][q]);
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
        }
        for k in 0..3 {
            let (apk, aqk) = (a[p][k], a[q][k]);
            a[p][k] = c * apk - s * aqk;
            a[q][k] = s * apk + c * aqk;
        }
        for k in 0..3 {
            let (vkp, vkq) = (v[k][p], v[k][q]);
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
        }
    }
    ([a[0][0], a[1][1], a[2][2]], v)
}
