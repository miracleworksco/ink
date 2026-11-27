//! The element matrices: trilinear hexes (with Wilson's incompatible modes,
//! for elasticity) as combinations of constant matrices weighted by a
//! cell's moments, and laminates'.

use crate::Physics;
use crate::cells::MIN_FRACTION;
use crate::hold;
use crate::laminate::Laminate;
use ink_core::{V3, v3};

pub(crate) fn signs(i: usize) -> [f64; 3] {
    [((i & 1) * 2) as f64 - 1.0, (((i >> 1) & 1) * 2) as f64 - 1.0, (((i >> 2) & 1) * 2) as f64 - 1.0]
}

/// Trilinear shape function i at local coordinates xi in [-1, 1]^3.
pub(crate) fn shape_n(i: usize, xi: V3) -> f64 {
    let s = signs(i);
    (1.0 + s[0] * xi[0]) * (1.0 + s[1] * xi[1]) * (1.0 + s[2] * xi[2]) / 8.0
}

/// d N_i / d xi.
pub(crate) fn shape_dn(i: usize, xi: V3) -> V3 {
    let s = signs(i);
    let f = [1.0 + s[0] * xi[0], 1.0 + s[1] * xi[1], 1.0 + s[2] * xi[2]];
    [s[0] * f[1] * f[2] / 8.0, f[0] * s[1] * f[2] / 8.0, f[0] * f[1] * s[2] / 8.0]
}

/// Strain-displacement matrix (6 x 24) in local coordinates (scale by 2/h for physical).
pub(crate) fn b_matrix(xi: V3) -> [[f64; 24]; 6] {
    let mut b = [[0.0; 24]; 6];
    for i in 0..8 {
        let d = shape_dn(i, xi);
        let c = 3 * i;
        b[0][c] = d[0];
        b[1][c + 1] = d[1];
        b[2][c + 2] = d[2];
        b[3][c] = d[1];
        b[3][c + 1] = d[0];
        b[4][c + 1] = d[2];
        b[4][c + 2] = d[1];
        b[5][c] = d[2];
        b[5][c + 2] = d[0];
    }
    b
}

/// Strain matrix (6 x 33) with Wilson's incompatible modes: the 24 nodal
/// columns, then 9 internal ones (1 - xi_k^2 for each axis k, each
/// displacement component). The internal modes let a brick bend without the
/// spurious shear that makes plain trilinear bricks too stiff when a wall is
/// only a few cells thick; they are condensed out per element.
fn b_ext(xi: V3) -> [[f64; 33]; 6] {
    let mut b = [[0.0; 33]; 6];
    let q = b_matrix(xi);
    for r in 0..6 {
        b[r][..24].copy_from_slice(&q[r]);
    }
    for k in 0..3 {
        let mut g = [0.0; 3];
        g[k] = -2.0 * xi[k];
        let c = 24 + 3 * k;
        b[0][c] = g[0];
        b[3][c] = g[1];
        b[5][c] = g[2];
        b[1][c + 1] = g[1];
        b[3][c + 1] = g[0];
        b[4][c + 1] = g[2];
        b[2][c + 2] = g[2];
        b[4][c + 2] = g[1];
        b[5][c + 2] = g[0];
    }
    b
}

/// Lamé's coefficients (lambda, mu) for E = 1.
fn lame(nu: f64) -> (f64, f64) {
    (nu / ((1.0 + nu) * (1.0 - 2.0 * nu)), 1.0 / (2.0 * (1.0 + nu)))
}

/// Isotropic elasticity matrix for E = 1.
pub(crate) fn d_matrix(nu: f64) -> [[f64; 6]; 6] {
    let (lam, mu) = lame(nu);
    d_lame(lam, mu)
}

/// The elasticity matrix lambda A + mu M: A is ones on the normal block,
/// M = diag(2, 2, 2, 1, 1, 1).
fn d_lame(lam: f64, mu: f64) -> [[f64; 6]; 6] {
    let mut d = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            d[i][j] = lam;
        }
        d[i][i] = lam + 2.0 * mu;
        d[i + 3][i + 3] = mu;
    }
    d
}

/// Element integrand (ne x ne, row-major) at xi, for unit coefficient.
fn integrand(p: Physics, xi: V3) -> Vec<f64> {
    match p {
        Physics::Heat => {
            let mut m = vec![0.0; 64];
            for i in 0..8 {
                let a = shape_dn(i, xi);
                for j in 0..8 {
                    let b = shape_dn(j, xi);
                    m[i * 8 + j] = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                }
            }
            m
        }
        Physics::Elastic { nu } => elastic_integrand(xi, &d_matrix(nu)),
    }
}

/// The elastic integrand B^T D B, 33 x 33 with the incompatible modes (see b_ext).
fn elastic_integrand(xi: V3, d: &[[f64; 6]; 6]) -> Vec<f64> {
    let b = b_ext(xi);
    let mut db = [[0.0; 33]; 6];
    for r in 0..6 {
        for c in 0..33 {
            db[r][c] = (0..6).map(|k| d[r][k] * b[k][c]).sum();
        }
    }
    let mut m = vec![0.0; 33 * 33];
    for i in 0..33 {
        for j in 0..33 {
            m[i * 33 + j] = (0..6).map(|k| b[k][i] * db[k][j]).sum();
        }
    }
    m
}

/// Coefficients over the 27 monomials of any integrand of degree <= 2 per axis.
pub(crate) fn poly_basis(f: impl Fn(V3) -> Vec<f64>) -> Vec<Vec<f64>> {
    let pts = [-1.0, 0.0, 1.0];
    let mut v: Vec<Vec<f64>> = (0..27).map(|k| f([pts[k % 3], pts[(k / 3) % 3], pts[k / 9]])).collect();
    let len = v[0].len();
    // Per axis: values at (-1, 0, 1) -> coefficients (1, t, t^2).
    for (axis, stride) in [(0, 1), (1, 3), (2, 9)] {
        let _ = axis;
        for k in 0..27 {
            if (k / stride) % 3 != 0 {
                continue;
            }
            let (i0, i1, i2) = (k, k + stride, k + 2 * stride);
            for e in 0..len {
                let (fm, f0, fp) = (v[i0][e], v[i1][e], v[i2][e]);
                v[i0][e] = f0;
                v[i1][e] = (fp - fm) / 2.0;
                v[i2][e] = (fp + fm) / 2.0 - f0;
            }
        }
    }
    v
}

/// Moments of the unit function over [-1, 1] for powers 0, 1, 2.
const MU: [f64; 3] = [2.0, 0.0, 2.0 / 3.0];

pub(crate) fn combine(basis: &[Vec<f64>], m: &[f64; 27], scale: f64) -> Vec<f64> {
    let mut out = vec![0.0; basis[0].len()];
    for (k, c) in basis.iter().enumerate() {
        if m[k] == 0.0 {
            continue;
        }
        let w = m[k] * scale;
        for (o, x) in out.iter_mut().zip(c) {
            *o += w * x;
        }
    }
    out
}

/// Element matrices from cell moments. Elastic elements carry Wilson's
/// incompatible modes with Taylor's correction (each mode's strain minus its
/// mean over the cell's material, so uniform stress stays exact on cut
/// cells), statically condensed to the 24 nodal dofs.
///
/// The integrand's coefficients are kept in parts, each weighted by one of
/// the material's coefficients: heat has one part; elasticity two, B^T A B
/// and B^T M B for D = lambda A + mu M (Lamé's coefficients). So one basis
/// serves every Poisson's ratio, each cell's exactly (`element_basis`).
pub(crate) struct ElementBasis {
    parts: Vec<Vec<Vec<f64>>>,
    /// Elastic: coefficients of [B_ext (6 x 33) | 1], for the mode means.
    g: Option<Vec<Vec<f64>>>,
}

/// The element basis of a physics (shared: it does not depend on the
/// material's coefficients).
pub(crate) fn element_basis(p: Physics) -> &'static ElementBasis {
    static HEAT: std::sync::OnceLock<ElementBasis> = std::sync::OnceLock::new();
    static ELASTIC: std::sync::OnceLock<ElementBasis> = std::sync::OnceLock::new();
    match p {
        Physics::Heat => HEAT.get_or_init(|| ElementBasis { parts: vec![poly_basis(|xi| integrand(Physics::Heat, xi))], g: None }),
        Physics::Elastic { .. } => ELASTIC.get_or_init(|| ElementBasis {
            parts: [d_lame(1.0, 0.0), d_lame(0.0, 1.0)].iter().map(|d| poly_basis(|xi| elastic_integrand(xi, d))).collect(),
            g: Some(poly_basis(|xi| {
                let b = b_ext(xi);
                b.iter().flatten().copied().chain(std::iter::once(1.0)).collect()
            })),
        }),
    }
}

impl ElementBasis {
    /// The integrand's moments over a cell for the material of `p`.
    fn combined(&self, p: Physics, m: &[f64; 27]) -> Vec<f64> {
        let weights = match p {
            Physics::Heat => vec![1.0],
            Physics::Elastic { nu } => {
                let (lam, mu) = lame(nu);
                vec![lam, mu]
            }
        };
        let mut out = combine(&self.parts[0], m, weights[0]);
        for (part, w) in self.parts.iter().zip(&weights).skip(1) {
            for (o, x) in out.iter_mut().zip(combine(part, m, *w)) {
                *o += x;
            }
        }
        out
    }

    /// Elastic: the 33 x 33 matrix with the incompatible modes (Taylor's
    /// correction applied), and the modes' coefficient table [B_ext | 1].
    pub(crate) fn extended(&self, nu: f64, m: &[f64; 27]) -> Option<(Vec<f64>, Vec<f64>)> {
        self.extended_d(&d_matrix(nu), self.combined(Physics::Elastic { nu }, m), m)
    }

    /// `extended` for a material of stiffness `d`, its integrand's
    /// moments `kk` (`combined`, or a laminate's: `aniso_basis`).
    fn extended_d(&self, d: &[[f64; 6]; 6], mut kk: Vec<f64>, m: &[f64; 27]) -> Option<(Vec<f64>, Vec<f64>)> {
        let gb = self.g.as_ref()?;
        let gm = combine(gb, m, 1.0);
        let vol = gm[198];
        let g = |r: usize, c: usize| gm[r * 33 + c];
        // Taylor's correction: K_a* -= G_a^T D G_* / V.
        if vol > 0.0 {
            for i in 24..33 {
                let dg: [f64; 6] = std::array::from_fn(|r| (0..6).map(|q| d[r][q] * g(q, i)).sum());
                for j in 0..33 {
                    let c: f64 = (0..6).map(|r| dg[r] * g(r, j)).sum::<f64>() / vol;
                    kk[i * 33 + j] -= c;
                    if j < 24 {
                        kk[j * 33 + i] -= c;
                    }
                }
            }
        }
        Some((kk, gm))
    }

    /// Elastic: the element's strain at each of `xis` for nodal displacements `ue`,
    /// with its incompatible modes recovered (K_aa a = -K_au u) and each
    /// mode's mean over the material removed, as in the stiffness.
    pub(crate) fn strains(&self, nu: f64, m: &[f64; 27], ue: &[f64; 24], xis: &[V3], h: f64) -> Option<Vec<[f64; 6]>> {
        Self::strains_ext(self.extended(nu, m)?, ue, xis, h)
    }

    /// `strains` from the element's `extended` matrices.
    fn strains_ext((kk, gm): (Vec<f64>, Vec<f64>), ue: &[f64; 24], xis: &[V3], h: f64) -> Option<Vec<[f64; 6]>> {
        // Solve the 9 x 9 system by Gaussian elimination with pivoting.
        let mut a = [[0.0f64; 10]; 9];
        let tr: f64 = (24..33).map(|i| kk[i * 33 + i].abs()).sum::<f64>().max(1e-300);
        for i in 0..9 {
            for j in 0..9 {
                a[i][j] = kk[(24 + i) * 33 + 24 + j] + if i == j { 1e-12 * tr } else { 0.0 };
            }
            a[i][9] = -(0..24).map(|c| kk[(24 + i) * 33 + c] * ue[c]).sum::<f64>();
        }
        for c in 0..9 {
            let p = (c..9).max_by(|&x, &y| a[x][c].abs().total_cmp(&a[y][c].abs()))?;
            a.swap(c, p);
            for r in 0..9 {
                if r != c && a[c][c] != 0.0 {
                    let f = a[r][c] / a[c][c];
                    for k in c..10 {
                        a[r][k] -= f * a[c][k];
                    }
                }
            }
        }
        let modes: [f64; 9] = std::array::from_fn(|i| if a[i][i] != 0.0 { a[i][9] / a[i][i] } else { 0.0 });
        let vol = gm[198];
        Some(
            xis.iter()
                .map(|&xi| {
                    let b = b_ext(xi);
                    std::array::from_fn(|r| {
                        let nodal: f64 = (0..24).map(|c| b[r][c] * ue[c]).sum();
                        let enh: f64 = (0..9).map(|i| (b[r][24 + i] - if vol > 0.0 { gm[r * 33 + 24 + i] / vol } else { 0.0 }) * modes[i]).sum();
                        (nodal + enh) * 2.0 / h
                    })
                })
                .collect(),
        )
    }

    /// Elastic: the strain operator (6 x 24, in xi: scale by 2/h) of an
    /// `extended` element (`kk`, `gm`), its incompatible modes condensed as
    /// `condense` condenses them, K_aa a = -K_au u, each mode's mean over
    /// the material taken out.
    fn strain_operator((kk, gm): (Vec<f64>, Vec<f64>), modes: bool) -> impl Fn(V3) -> [[f64; 24]; 6] + Sync + use<> {
        let tr: f64 = (24..33).map(|i| kk[i * 33 + i].abs()).sum::<f64>().max(1e-300);
        let kaa: Vec<f64> = (0..81).map(|ij| kk[(24 + ij / 9) * 33 + 24 + ij % 9] + if ij % 10 == 0 { 1e-12 * tr } else { 0.0 }).collect();
        // T = -K_aa^-1 K_au, column by column.
        let l = modes.then(|| hold::cholesky(&kaa, 9)).flatten();
        let cols: Vec<Vec<f64>> = (0..24).map(|c| l.as_ref().map_or(vec![0.0; 9], |l| hold::solve(l, 9, &(0..9).map(|i| -kk[(24 + i) * 33 + c]).collect::<Vec<f64>>()))).collect();
        let vol = gm[198];
        let mean: [[f64; 9]; 6] = std::array::from_fn(|r| std::array::from_fn(|i| if vol > 0.0 { gm[r * 33 + 24 + i] / vol } else { 0.0 }));
        move |xi: V3| {
            let b = b_ext(xi);
            std::array::from_fn(|r| std::array::from_fn(|c| b[r][c] + (0..9).map(|i| (b[r][24 + i] - mean[r][i]) * cols[c][i]).sum::<f64>()))
        }
    }

    /// The element matrix (ne x ne) for the material of `p` and cell
    /// moments `m`, times `scale`.
    pub(crate) fn matrix(&self, p: Physics, m: &[f64; 27], scale: f64) -> Vec<f64> {
        let Physics::Elastic { nu } = p else { return self.combined(p, m).into_iter().map(|v| v * scale).collect() };
        Self::condense(&self.extended(nu, m).expect("an elastic basis").0, scale)
    }

    /// The 24 x 24 matrix of an `extended` one, its modes condensed out,
    /// times `scale`.
    fn condense(kk: &[f64], scale: f64) -> Vec<f64> {
        // Condense: K = K_uu - K_ua K_aa^-1 K_au (Cholesky of the 9 x 9).
        let mut l = [[0.0f64; 9]; 9];
        let tr: f64 = (24..33).map(|i| kk[i * 33 + i].abs()).sum::<f64>().max(1e-300);
        for j in 0..9 {
            let v = kk[(24 + j) * 33 + 24 + j] + 1e-12 * tr - (0..j).map(|q| l[j][q] * l[j][q]).sum::<f64>();
            l[j][j] = v.max(1e-300 * tr).sqrt();
            for i in j + 1..9 {
                l[i][j] = (kk[(24 + i) * 33 + 24 + j] - (0..j).map(|q| l[i][q] * l[j][q]).sum::<f64>()) / l[j][j];
            }
        }
        // Y = L^-1 K_au (9 x 24); K -= Y^T Y.
        let mut y = [[0.0f64; 24]; 9];
        for c in 0..24 {
            for i in 0..9 {
                y[i][c] = (kk[(24 + i) * 33 + c] - (0..i).map(|q| l[i][q] * y[q][c]).sum::<f64>()) / l[i][i];
            }
        }
        let mut out = vec![0.0; 576];
        for i in 0..24 {
            for j in 0..24 {
                out[i * 24 + j] = scale * (kk[i * 33 + j] - (0..9).map(|q| y[q][i] * y[q][j]).sum::<f64>());
            }
        }
        out
    }
}

/// The traction of a cell's field on a surface through it, per dof, at a
/// point (xi) and the surface's normal: ndof x ne, for the material of
/// moments `m` (coefficient-weighted) and the stiffness its matrix has: its
/// laminate's (`laminate`: with the cell's mean modulus, isotropic, a hold
/// across a laminate's interface was a penalty of about E / h), else the
/// cell's modulus `coef`'s.
/// Elastic: the stress of its element (its incompatible modes condensed
/// as its matrix condenses them) on the normal; heat: the flux k grad T . n.
pub(crate) fn traction<'a>(basis: &'a ElementBasis, p: Physics, m: &[f64; 27], coef: f64, size: f64, lam: Option<&Laminate>) -> Box<dyn Fn(V3, V3) -> Vec<f64> + Sync + 'a> {
    let to = 2.0 / size;
    // (A sliver's, without its incompatible modes, as its stiffness.)
    let modes = m[0] >= MIN_FRACTION * 8.0 * coef;
    match p {
        Physics::Heat => {
            let k: Vec<f64> = lam.map_or_else(|| (0..9).map(|i| if i % 4 == 0 { coef } else { 0.0 }).collect(), |l| l.effective(false));
            Box::new(move |xi: V3, n: V3| (0..8).map(|i| to * (0..9).map(|ab| n[ab / 3] * k[ab] * shape_dn(i, xi)[ab % 3]).sum::<f64>()).collect())
        }
        Physics::Elastic { nu } => {
            let (strain, d) = match lam {
                Some(l) => (ElementBasis::strain_operator(laminate_extended(l, &m.map(|v| v / coef)).expect("an elastic basis"), modes), l.effective6()),
                None => (ElementBasis::strain_operator(basis.extended(nu, m).expect("an elastic basis"), modes), d_matrix(nu).map(|r| r.map(|v| v * coef))),
            };
            Box::new(move |xi: V3, n: V3| {
                let b = strain(xi);
                let mut t = vec![0.0; 72];
                for c in 0..24 {
                    let s: [f64; 6] = std::array::from_fn(|r| to * (0..6).map(|k| d[r][k] * b[k][c]).sum::<f64>());
                    t[c] = s[0] * n[0] + s[3] * n[1] + s[5] * n[2];
                    t[24 + c] = s[3] * n[0] + s[1] * n[1] + s[4] * n[2];
                    t[48 + c] = s[5] * n[0] + s[4] * n[1] + s[2] * n[2];
                }
                t
            })
        }
    }
}

/// The element basis of an anisotropic material (a laminate's,
/// `laminate`): one part per entry of its stiffness's upper triangle (6 x
/// 6, Voigt; heat's conductivity 3 x 3), each the integrand for a stiffness
/// of ones there and at its mirror.
struct AnisoBasis {
    n: usize,
    parts: Vec<((usize, usize), Vec<Vec<f64>>)>,
}

fn aniso_basis(p: Physics) -> &'static AnisoBasis {
    static HEAT: std::sync::OnceLock<AnisoBasis> = std::sync::OnceLock::new();
    static ELASTIC: std::sync::OnceLock<AnisoBasis> = std::sync::OnceLock::new();
    let make = |n: usize, f: &dyn Fn(V3, &[f64]) -> Vec<f64>| AnisoBasis {
        n,
        parts: (0..n)
            .flat_map(|k| (k..n).map(move |l| (k, l)))
            .map(|(k, l)| {
                let mut d = vec![0.0; n * n];
                d[k * n + l] = 1.0;
                d[l * n + k] = 1.0;
                ((k, l), poly_basis(|xi| f(xi, &d)))
            })
            .collect(),
    };
    match p {
        Physics::Heat => HEAT.get_or_init(|| {
            make(3, &|xi, k| {
                let mut m = vec![0.0; 64];
                for i in 0..8 {
                    let a = shape_dn(i, xi);
                    for j in 0..8 {
                        let b = shape_dn(j, xi);
                        m[i * 8 + j] = (0..3).map(|p| (0..3).map(|q| a[p] * k[p * 3 + q] * b[q]).sum::<f64>()).sum();
                    }
                }
                m
            })
        }),
        Physics::Elastic { .. } => ELASTIC.get_or_init(|| make(6, &|xi, d| elastic_integrand(xi, &std::array::from_fn(|i| std::array::from_fn(|j| d[i * 6 + j]))))),
    }
}

impl AnisoBasis {
    /// The integrand's moments over a cell (`m`, its material's: not
    /// weighed by a coefficient) for stiffness `d` (n x n, symmetric).
    fn combined(&self, d: &[f64], m: &[f64; 27]) -> Vec<f64> {
        let mut out = vec![0.0; self.parts[0].1[0].len()];
        for ((k, l), part) in &self.parts {
            let w = d[k * self.n + l];
            if w != 0.0 {
                for (o, x) in out.iter_mut().zip(combine(part, m, w)) {
                    *o += x;
                }
            }
        }
        out
    }
}

/// A laminate's (`laminate`) element: its `extended` matrices (elastic),
/// for its material's moments `m`.
pub(crate) fn laminate_extended(lam: &Laminate, m: &[f64; 27]) -> Option<(Vec<f64>, Vec<f64>)> {
    let p = Physics::Elastic { nu: 0.0 };
    let d = lam.effective6();
    let kk = aniso_basis(p).combined(&lam.effective(true), m);
    element_basis(p).extended_d(&d, kk, m)
}

/// A laminate's mean stress for its mean strain `e` (Voigt), or its mean
/// heat flux for its mean temperature gradient (the first 3).
pub(crate) fn laminate_mean(lam: &Laminate, elastic: bool, e: &[f64]) -> [f64; 6] {
    let k = lam.effective(elastic);
    let n = if elastic { 6 } else { 3 };
    let sign = if elastic { 1.0 } else { -1.0 };
    std::array::from_fn(|r| if r < n { sign * (0..n).map(|q| k[r * n + q] * e[q]).sum::<f64>() } else { 0.0 })
}

/// A laminate element's mean strain at each of `xis` (`ElementBasis::strains`).
pub(crate) fn laminate_strains(lam: &Laminate, m: &[f64; 27], ue: &[f64; 24], xis: &[V3], h: f64) -> Option<Vec<[f64; 6]>> {
    ElementBasis::strains_ext(laminate_extended(lam, m)?, ue, xis, h)
}

/// A laminate's element matrix (ne x ne), times `scale`.
pub(crate) fn laminate_matrix(lam: &Laminate, p: Physics, m: &[f64; 27], scale: f64) -> Vec<f64> {
    match p {
        Physics::Heat => aniso_basis(p).combined(&lam.effective(false), m).into_iter().map(|v| v * scale).collect(),
        Physics::Elastic { .. } => ElementBasis::condense(&laminate_extended(lam, m).expect("an elastic basis").0, scale),
    }
}

/// An `extended` element's matrix without its incompatible modes (ne x ne),
/// times `scale`.
pub(crate) fn uncondensed(kk: &[f64], scale: f64) -> Vec<f64> {
    (0..576).map(|ij| scale * kk[(ij / 24) * 33 + ij % 24]).collect()
}

pub(crate) fn full_moments() -> [f64; 27] {
    std::array::from_fn(|k| MU[k % 3] * MU[(k / 3) % 3] * MU[k / 9])
}

/// A material's moments (of xi^a eta^b zeta^c, a b c to 2 each, at a + 3 b +
/// 9 c) in coordinates moved by `d`, where a point at xi is at xi + d:
/// each a sum of the moments of no higher power, by the binomial.
pub(crate) fn moved_moments(m: &[f64; 27], d: V3) -> [f64; 27] {
    // (x + d)^a = sum over i of C(a, i) d^(a - i) x^i
    let term = |a: usize, i: usize, d: f64| [[1.0, 0.0, 0.0], [d, 1.0, 0.0], [d * d, 2.0 * d, 1.0]][a][i];
    std::array::from_fn(|k| {
        let (a, b, c) = (k % 3, k / 3 % 3, k / 9);
        (0..=a)
            .flat_map(|i| (0..=b).flat_map(move |j| (0..=c).map(move |l| (i, j, l))))
            .map(|(i, j, l)| term(a, i, d[0]) * term(b, j, d[1]) * term(c, l, d[2]) * m[i + 3 * j + 9 * l])
            .sum()
    })
}

/// Local sample coordinates for k samples per axis (sub-cell midpoints).
pub(crate) fn sample_xi(k: usize) -> Vec<V3> {
    let t = |j: usize| -1.0 + (2 * j + 1) as f64 / k as f64;
    let mut out = Vec::with_capacity(k * k * k);
    for c in 0..k {
        for b in 0..k {
            for a in 0..k {
                out.push([t(a), t(b), t(c)]);
            }
        }
    }
    out
}

pub(crate) type M3 = [[f64; 3]; 3];

pub(crate) fn m3_det(a: &M3) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// The size of a cell's small rotation at its centre: half the curl of
/// its displacement (radians).
pub(crate) fn small_turn(ue: &[f64; 24], h: f64) -> f64 {
    // g[d][a] = du_d / dx_a at the centre.
    let mut g = [[0.0; 3]; 3];
    for k in 0..8 {
        let dn = shape_dn(k, [0.0; 3]);
        for (d, row) in g.iter_mut().enumerate() {
            let u = ue[k * 3 + d] * 2.0 / h;
            for a in 0..3 {
                row[a] += dn[a] * u;
            }
        }
    }
    let w = [g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1]];
    0.5 * v3::norm(w)
}

#[cfg(test)]
mod tests {

    #[test]
    fn a_materials_moments_move_with_its_coordinates() {
        // A whole cell's, seen from the cell two over in x and one back in
        // z: the cube 3..5 by -1..1 by -3..-1, each moment a product of
        // its sides' integrals.
        let side = |lo: f64, hi: f64, a: usize| (hi.powi(a as i32 + 1) - lo.powi(a as i32 + 1)) / (a as f64 + 1.0);
        let moved = super::moved_moments(&super::full_moments(), [4.0, 0.0, -2.0]);
        for k in 0..27 {
            let want = side(3.0, 5.0, k % 3) * side(-1.0, 1.0, k / 3 % 3) * side(-3.0, -1.0, k / 9);
            assert!((moved[k] - want).abs() < 1e-12 * want.abs().max(1.0), "moment {k}: {} for {want}", moved[k]);
        }
    }
}
