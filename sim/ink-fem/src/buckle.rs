//! Linear buckling: the geometric stiffness of a solve's stress, and the
//! factor at which its loads buckle the part.

use crate::eigen::{dense_sym_eigen, lobpcg};
use crate::element::{combine, d_matrix, element_basis, full_moments, laminate_mean, laminate_strains, poly_basis, shape_dn};
use crate::multigrid::Multigrid;
use crate::numerics::dot;
use crate::report::Buckle;
use crate::system::{Ptr, System};
use crate::{CELL_GRAIN, GRAIN, MAX_ITER, Physics};
use ink_core::V3;
use rayon::prelude::*;
use std::time::Instant;

/// Whether the symmetric operator `a` (a stiffness plus a load factor times its geometric stiffness)
/// is positive definite: the loads at that factor buckle nothing. By conjugate gradients on a x = b
/// for a random b, preconditioned by `t` (an approximate inverse of the stiffness). A step along
/// which a has no positive curvature is a direction that buckles: Err, a proof, whatever `t` is (the
/// signs of a's eigenvalues are not the preconditioner's to change). Converged with none, to a
/// residual a hundred-thousandth of b's: Ok. b is random at every node, so every mode has a part in it (its
/// share about one over the square root of the unknowns, far over that residual), and the iteration
/// cannot converge past a mode without meeting its curvature: unlike a search for the lowest mode
/// from smooth starts, it takes no mode as the lowest and reads no residual as a distance. None: not
/// converged (its residual no longer halving in `STABLE_ITER` steps), and nothing learned. `own` puts the right-hand side on the
/// unknowns the operator owns (a hanging node's part to those it hangs on).
fn stable_at(n: usize, a: &dyn Fn(&[f64]) -> Vec<f64>, t: &dyn Fn(&[f64]) -> Vec<f64>, own: &dyn Fn(&mut Vec<f64>)) -> Option<Result<(), ()>> {
    let mut seed = 0x9E37_79B9_7F4A_7C15u64;
    let b: Vec<f64> = (0..n)
        .map(|_| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        })
        .collect();
    // (As a load is: what lands on a node the operator does not own, to those that own it.)
    let mut b = b;
    own(&mut b);
    let mut r = b;
    let mut now = dot(&r, &r).sqrt();
    let bnorm = now.max(1e-300);
    let mut z = t(&r);
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    // As the static solve goes on (`pcg`): while the residual halves within 500 steps.
    let (mut since, mut then) = (0usize, f64::INFINITY);
    for _ in 0..MAX_ITER {
        if since >= STABLE_ITER {
            if now > 0.5 * then {
                return None;
            }
            (since, then) = (0, now);
        }
        since += 1;
        let ap = a(&p);
        let pap = dot(&p, &ap);
        if !pap.is_finite() {
            return None;
        }
        if pap <= 0.0 {
            return Some(Err(()));
        }
        let alpha = rz / pap;
        r.par_iter_mut().with_min_len(GRAIN).zip(&ap).for_each(|(r, a)| *r -= alpha * a);
        now = dot(&r, &r).sqrt();
        if now <= STABLE_TOL * bnorm {
            return Some(Ok(()));
        }
        z = t(&r);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        p.par_iter_mut().with_min_len(GRAIN).zip(&z).for_each(|(p, z)| *p = z + beta * *p);
    }
    None
}

/// The most a buckling factor's low end is checked to (`buckle`): a hundred times the loads, past
/// any margin a design asks of buckling.
const STABLE_TO: f64 = 100.0;

/// `stable_at`'s residual, over its right-hand side's: far under any mode's share of a random vector
/// (one over the square root of the unknowns: 7e-4 of two million).
const STABLE_TOL: f64 = 1e-5;
/// The steps `stable_at` goes on for without its residual halving before it stops, as the static
/// solve's own guard (a thin fin held too stiff takes it 1100 steps, and the check as many).
const STABLE_ITER: usize = 500;

/// Vectors the buckling eigensolver keeps in its block (the first its
/// answer, the second speeding it where modes crowd: a square column's
/// two equal ones, a shell's many close ones).
const BUCKLING_BLOCK: usize = 2;
/// The buckling eigensolver stops where the first pair's residual, in
/// the norm of K^-1 (the preconditioner's), is under this share of its
/// value nu: the distance from nu to an eigenvalue is at most that norm,
/// so the factor is resolved to this share of itself, whatever the load's
/// scale.
const BUCKLING_TOL: f64 = 1e-2;
/// Resolved enough where the residual bounds the least factor at 1 / this
/// or above (1000: loads a thousandth of what buckles the part). Where the
/// loads hardly compress it (a bar pulled, a beam whose compression is
/// local) its values crowd at zero, where a share of themselves is not
/// reached in the iterations; that bound is what its claims need. A
/// pulled bar's search stops at 20 steps, three times its static solve; a
/// floor of 1e-5 finds its clamp's corner mode (about 16,300) and fails to
/// resolve it in 100.
const BUCKLING_NU_FLOOR: f64 = 1e-3;
/// Steps of the buckling eigensolve before that bound may end it.
const BUCKLING_SEARCH: usize = 20;
/// Most iterations of the buckling eigensolve (the tests' tube takes about 50).
const BUCKLING_ITER: usize = 100;

/// A free turn's geometric force is none where it is under this share of
/// what the prestress does to an arbitrary motion of the turn's size
/// (`buckling`'s free turns): a stress that is none is known no closer
/// than the solve and the cells' place leave it, 1.7e-6 the most measured
/// (a bar pulled between mirrors, across its pull, on cells of 0.6 mm),
/// where a turn the prestress works on reads 0.03 to 0.3.
const NO_PRESTRESS: f64 = 1e-5;

/// Per axis pair (a, b) and nodes (i, j), the coefficients of
/// dN_i/dxi_a dN_j/dxi_b over the 27 monomials (`poly_basis`): a cell's
/// geometric stiffness from its material's moments.
fn gradient_basis() -> &'static Vec<Vec<f64>> {
    static B: std::sync::OnceLock<Vec<Vec<f64>>> = std::sync::OnceLock::new();
    B.get_or_init(|| {
        poly_basis(|xi| {
            let g: Vec<V3> = (0..8).map(|i| shape_dn(i, xi)).collect();
            (0..9 * 64).map(|k| g[(k % 64) / 8][k / 192] * g[k % 8][(k / 64) % 3]).collect()
        })
    })
}

/// A cell's geometric stiffness per displacement component (8 x 8; the
/// element's is it on each of x, y, z): the integral over its material
/// of sigma_ab dN_i/dx_a dN_j/dx_b, for a stress constant over it
/// (Voigt), from its gradients' integrals `g` (`gradient_basis` combined
/// with its moments). By columns (entry i, j at j * 8 + i), as the
/// stiffness's are (`System::new`).
fn geometric(g: &[f64], sigma: &[f64; 6], h: f64) -> [f64; 64] {
    let s = [[sigma[0], sigma[3], sigma[5]], [sigma[3], sigma[1], sigma[4]], [sigma[5], sigma[4], sigma[2]]];
    std::array::from_fn(|ji| {
        let ij = (ji % 8) * 8 + ji / 8;
        0.5 * h * (0..9).map(|ab| s[ab / 3][ab % 3] * g[ab * 64 + ij]).sum::<f64>()
    })
}

/// A cell's held material, for its geometric stiffness: the moments (of
/// unit coefficient) of its material outside every region, and on each axis
/// a region holds, outside those that hold it; and the axes held.
pub(crate) type HeldGeometry = ([f64; 27], [[f64; 27]; 3], [bool; 3]);

/// Linear buckling of a solved static elastic study: the lowest load
/// factors lambda of K phi = -lambda K_G phi, as the most negative nu of
/// K_G phi = nu K phi (lambda = -1 / nu: K positive definite, held by its
/// supports, so LOBPCG with K as the inner product and multigrid, an
/// approximate K^-1, as preconditioner; K^-1 K_G's buckling end is well
/// separated, as an inverse iteration's). K_G is each cell's stress at
/// its material's centroid (`cell_field`'s at a full cell's centre) times
/// its material's gradient integrals. With its mode, 1 mm at its peak.
/// A cell holding held material (`held`: its free material's moments and
/// the axes held) is its free material's on those axes: the held
/// material's motion is the hold's. (On every axis, a part held along one
/// all through, a beam's guide, had nothing left to buckle.)
///
/// It runs on one of the solver's threads: from any other, each of the
/// hundreds of parallel steps of a multigrid cycle waits for a thread to
/// wake and take it up (4% of benchmark 05's search and check).
#[allow(clippy::too_many_arguments)]
pub(crate) fn buckle<'a>(
    sys: &System,
    mg: &Multigrid,
    u: &[f64],
    cell_field: &[[f64; 6]],
    held: &[Option<Box<HeldGeometry>>],
    h: f64,
    pos: &[V3],
    rotations: impl Iterator<Item = &'a Vec<f64>>,
) -> (Buckle, Vec<Vec<f64>>) {
    let rotations: Vec<&Vec<f64>> = rotations.collect();
    rayon::scope(|_| buckle_here(sys, mg, u, cell_field, held, h, pos, &rotations))
}

#[allow(clippy::too_many_arguments)]
fn buckle_here(
    sys: &System,
    mg: &Multigrid,
    u: &[f64],
    cell_field: &[[f64; 6]],
    held: &[Option<Box<HeldGeometry>>],
    h: f64,
    pos: &[V3],
    rotations: &[&Vec<f64>],
) -> (Buckle, Vec<Vec<f64>>) {
    let t0 = Instant::now();
    let basis = element_basis(Physics::Elastic { nu: 0.3 });
    // Each cell's stress in its geometric stiffness: a cut cell's at its
    // material's centroid, a full cell's at its centre (`cell_field`); a
    // cell holding held material, its free material's at its centroid, and
    // none where all of it is held on every axis.
    let stress: Vec<[f64; 6]> = sys
        .cells
        .par_iter()
        .zip(cell_field)
        .zip(held)
        .map(|((c, &centre), held)| {
            let m = match held.as_deref() {
                Some((f, _, _)) if f[0] > 0.0 => f,
                Some((_, _, axes)) if axes.iter().all(|a| *a) => return [0.0; 6],
                _ => match c.material.as_deref() {
                    Some(m) => m,
                    None => return centre,
                },
            };
            let size = h / f64::from(1u32 << c.level);
            let xi: V3 = [m[1], m[3], m[9]].map(|v| (v / m[0].max(1e-300)).clamp(-1.0, 1.0));
            let ue: [f64; 24] = std::array::from_fn(|j| u[c.nodes[j / 3] as usize * 3 + j % 3]);
            if let Some(lam) = &c.lam {
                return laminate_strains(lam, m, &ue, &[xi], size).map_or(centre, |e| laminate_mean(lam, true, &e[0]).map(|v| v * c.mult));
            }
            match basis.strains(c.nu, m, &ue, &[xi], size) {
                Some(e) => {
                    let d = d_matrix(c.nu);
                    std::array::from_fn(|r| c.coef * c.mult * (0..6).map(|q| d[r][q] * e[0][q]).sum::<f64>())
                }
                None => centre,
            }
        })
        .collect();
    let no_mode = |low| (Buckle { factor: f64::INFINITY, low, ms: t0.elapsed().as_secs_f64() * 1e3, ..Default::default() }, Vec::new());
    let has_rotations = !rotations.is_empty();
    // With no free turns to check, an uncompressed study needs only its
    // stress. Leave the geometric stiffness unassembled.
    if !has_rotations && let Some(low) = uncompressed(sys, &stress, h, pos) {
        return no_mode(low);
    }

    // The prestress stays fixed through the search and stability check.
    // Assemble full cells' geometric stiffness once too, reusing each
    // matrix on every application.
    let grads = gradient_basis();
    let full = combine(grads, &full_moments(), 1.0);
    let matrices: Vec<[f64; 64]> = sys
        .cells
        .par_iter()
        .zip(&stress)
        .map(|(c, sigma)| {
            let size = h / f64::from(1u32 << c.level);
            match c.material.as_deref() {
                Some(m) => geometric(&combine(grads, m, 1.0), sigma, size),
                None => geometric(&full, sigma, size),
            }
        })
        .collect();
    // Only held cells carry separate matrices for their held axes.
    let free: Vec<_> = sys
        .cells
        .par_iter()
        .zip(&stress)
        .zip(held)
        .map(|((c, sigma), held)| {
            held.as_deref().map(|&(_, ref per, axes)| {
                let size = h / f64::from(1u32 << c.level);
                let g: [[f64; 64]; 3] = std::array::from_fn(|d| geometric(&combine(grads, &per[d], 1.0), sigma, size));
                (Box::new(g), axes)
            })
        })
        .collect();
    // A cell's geometric stiffness applied to `x`, added to `y`: a split
    // cell's at its own size, its hanging corners through what they hang
    // on (`System::weighted_dofs`), as the stiffness is. Each row's sum is
    // in order from -0.0, as its own `sum` would be, the rows of an axis
    // side by side (`System::element`).
    let cell_kg = |ci: usize, x: &[f64], y: &Ptr| {
        let c = &sys.cells[ci];
        let dofs = sys.hangs(c).then(|| sys.weighted_dofs(c));
        let ue: [f64; 24] = match &dofs {
            None => std::array::from_fn(|j| x[c.nodes[j / 3] as usize * 3 + j % 3]),
            Some(dofs) => std::array::from_fn(|j| dofs[j].iter().map(|&(n, w)| w * x[n as usize]).sum()),
        };
        let mut v = [[-0.0; 8]; 3];
        for (d, v) in v.iter_mut().enumerate() {
            let g = match &free[ci] {
                Some((g, axes)) if axes[d] => &g[d],
                _ => &matrices[ci],
            };
            for (column, j) in g.as_chunks::<8>().0.iter().zip(0..8) {
                let u = ue[j * 3 + d];
                for (v, g) in v.iter_mut().zip(column) {
                    *v += g * u;
                }
            }
        }
        for i in 0..8 {
            for (d, v) in v.iter().enumerate() {
                match &dofs {
                    // SAFETY: as below.
                    None => unsafe { *y.0.add(c.nodes[i] as usize * 3 + d) += v[i] },
                    Some(dofs) => {
                        for &(n, w) in &dofs[i * 3 + d] {
                            // SAFETY: cells of one colour (a grid cell's split
                            // cells, one run, in order) share no nodes, nor any
                            // their corners hang on, so no two threads write the
                            // same dof in this phase.
                            unsafe { *y.0.add(n as usize) += w * v[i] };
                        }
                    }
                }
            }
        }
    };
    let geometric_force = |x: &[f64]| -> Vec<f64> {
        let mut y = vec![0.0; x.len()];
        let yp = Ptr(y.as_mut_ptr());
        for color in &sys.colors {
            color.par_iter().with_min_len(CELL_GRAIN).for_each(|&ci| cell_kg(ci as usize, x, &yp));
        }
        for runs in &sys.runs {
            runs.par_iter().for_each(|&(a, b)| (a as usize..b as usize).for_each(|ci| cell_kg(ci, x, &yp)));
        }
        y
    };
    // K has a null space of free rigid motions. Translations have zero
    // gradient. Rotations can carry prestress energy: choosing their
    // amplitude by the static solver's rest pin would change the factor.
    // Minimise geometric energy over them first (Schur complement).
    // G-orthonormal turns r give G_eff = G - sum (Gr)(Gr)^T. Thus the
    // arbitrary rest pins in K only pick representatives of elastic modes.
    let mut turns: Vec<(Vec<f64>, Vec<f64>)> = Vec::new();
    // What the prestress does to an arbitrary motion, for each unit of it:
    // what a turn's force is measured against.
    let prestress = std::cell::OnceCell::new();
    let arbitrary = || {
        let mut seed = 0x9E37_79B9_7F4A_7C15u64;
        let z: Vec<f64> = (0..sys.ndofs)
            .map(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
            })
            .collect();
        let g = geometric_force(&z);
        (dot(&g, &g) / dot(&z, &z).max(1e-300)).sqrt()
    };
    for &rigid in rotations {
        let mut r = rigid.clone();
        for (basis, force) in &turns {
            let coefficient = dot(&r, force);
            for (x, b) in r.iter_mut().zip(basis) {
                *x -= coefficient * b;
            }
        }
        let mut force = geometric_force(&r);
        // A null vector of both operators is just a free coordinate:
        // neither its amplitude nor its coupling enters the eigenproblem.
        // Zero rotational energy alone would not establish this; the
        // prestress doing nothing to it does, and nothing is what the
        // solve leaves of a stress that is none (`NO_PRESTRESS`). A block
        // on a roller floor between two mirrors under its own weight has
        // no stress across its weight, and a turn about it reads 1e-15 to
        // 1e-11 of what an arbitrary motion does, of either sign: held to
        // exactly zero, the same study's buckling was resolved by one solve
        // and "not resolved" by another solved closer, and a claim on it
        // settled or not by which of the two it was given.
        if dot(&force, &force).sqrt() <= NO_PRESTRESS * *prestress.get_or_init(arbitrary) * dot(&r, &r).sqrt() {
            continue;
        }
        let energy = dot(&r, &force);
        let scale = dot(&r, &r).sqrt() * dot(&force, &force).sqrt();
        if !(energy > 1e-10 * scale && energy.is_finite()) {
            return (
                Buckle {
                    failure: Some(
                        "prestress gives a free rigid rotation no resolved positive restoring energy. Model its physical restraint or check the whole part's stability".into(),
                    ),
                    ..Default::default()
                },
                Vec::new(),
            );
        }
        let norm = energy.sqrt();
        r.iter_mut().chain(&mut force).for_each(|v| *v /= norm);
        turns.push((r, force));
    }
    // No cell compressed, nothing buckles, however short a search would
    // stop. A cell's K_G is the integral over its material of
    // sum_d grad(phi_d)^T sigma grad(phi_d), so it is positive semidefinite
    // where sigma is (no negative principal stress); a hanging corner's
    // weights are a congruence, which keeps it so; and the turns taken
    // out above are G-orthonormal, so G - sum (G r)(G r)^T is G^(1/2) times
    // a projector times G^(1/2), semidefinite too. With K semidefinite,
    // K phi = -lambda K_G phi then has no lambda > 0 but for a free
    // coordinate of both (a motion neither resists: not a buckle).
    if has_rotations && let Some(low) = uncompressed(sys, &stress, h, pos) {
        return no_mode(low);
    }
    let kg = |x: &[f64]| {
        let mut y = geometric_force(x);
        for (_, force) in &turns {
            let coefficient = dot(force, x);
            for (y, f) in y.iter_mut().zip(force) {
                *y -= coefficient * f;
            }
        }
        y
    };
    let stiffness = |v: &[f64]| {
        let mut o = vec![0.0; v.len()];
        sys.apply(v, &mut o);
        o
    };
    let precond = |r: &[f64]| mg.vcycle(0, sys, r);
    // The residual's K^-1 norm (the preconditioner's) for a mode of unit
    // K norm: its value's distance from an eigenvalue, at most.
    let distance = |r: &[f64], w: &[f64]| dot(r, w).max(0.0).sqrt();
    // (The bound at 1000 only after BUCKLING_SEARCH steps: a few steps in,
    // the block may not yet have found the lowest mode it bounds.)
    let steps = std::cell::Cell::new(0usize);
    let error = |nu: f64, r: &[f64], w: &[f64], _: &[f64]| {
        steps.set(steps.get() + 1);
        let d = distance(r, w);
        let resolved = d / (BUCKLING_TOL * nu.abs());
        if steps.get() >= BUCKLING_SEARCH { resolved.min((d - nu) / BUCKLING_NU_FLOOR) } else { resolved }
    };
    let (nu, vecs, iterations, converged) = lobpcg(sys.ndofs, &kg, &stiffness, &|r| mg.vcycle_each(sys, r), (BUCKLING_BLOCK, 1), &error, BUCKLING_ITER, &[]);
    if vecs.is_empty() {
        let failure = Some("the buckling eigensolver broke down: no direction left".to_string());
        return (Buckle { antisymmetric: [false; 3], factor: f64::NAN, low: f64::NAN, at: [0.0; 3], iterations, ms: t0.elapsed().as_secs_f64() * 1e3, failure }, Vec::new());
    }
    let delta = {
        let (x, kx) = (&vecs[0], stiffness(&vecs[0]));
        let r: Vec<f64> = kg(x).iter().zip(&kx).map(|(a, k)| a - nu[0] * k).collect();
        distance(&r, &precond(&r))
    };
    let shape: Vec<f64> = {
        let mut x = vecs[0].clone();
        for (r, force) in &turns {
            let coefficient = dot(force, &vecs[0]);
            for (v, r) in x.iter_mut().zip(r) {
                *v -= coefficient * r;
            }
        }
        sys.expand(&mut x, 3);
        let peak = (0..x.len() / 3).map(|n| (0..3).map(|d| x[n * 3 + d].powi(2)).sum::<f64>().sqrt()).fold(0.0, f64::max).max(1e-300);
        x.iter().map(|v| v / peak).collect()
    };
    // Converged, the lowest nu is within delta of the one found (taking
    // it as the lowest), so the least factor is at least 1 / (delta - nu).
    // Where no direction it found has a negative nu (a bar pulled), it
    // found no buckling mode: the factor is infinite, and the bound what
    // its search excluded, no more (a pulled bar's clamp compresses its
    // corners, which buckle at a factor in the ten thousands, a mode a
    // longer search finds).
    let failure = if !(nu[0].is_finite() && delta.is_finite() && shape.iter().all(|v| v.is_finite())) {
        Some("the buckling eigensolver broke down: a value not finite".to_string())
    } else if !converged {
        Some(format!(
            "the buckling eigensolver did not resolve the lowest buckling factor to {}% in {iterations} iterations: the loads' \
             lowest modes too close for it to part",
            BUCKLING_TOL * 100.0
        ))
    } else {
        None
    };
    let none = nu[0] >= 0.0;
    let norm2 = |n: usize| (0..3).map(|d| shape[n * 3 + d].powi(2)).sum::<f64>();
    let at = (0..shape.len() / 3).max_by(|&a, &b| norm2(a).total_cmp(&norm2(b))).map_or([0.0; 3], |n| pos[n]);
    // The Rayleigh quotient bounds the lowest nu from above, so its factor
    // the lowest factor from above; the residual bounds it below.
    let low = if delta > nu[0] { 1.0 / (delta - nu[0]) } else { f64::INFINITY };
    let factor = if none { f64::INFINITY } else { -1.0 / nu[0] };
    // That low end takes the mode found as the lowest, by a residual read through the
    // preconditioner. It is checked by what assumes neither (`stable_at`): the loads times it
    // must leave the part stable, and where they do not, the largest of its halves that does is
    // the low end (none in ten: nothing excluded). So the low end is what was checked.
    let (mut low, mut factor, mut checks, mut unchecked) = (low, factor, 0, false);
    if failure.is_none() && low.is_finite() && low > 0.0 {
        let stable = |x: f64| stable_at(sys.ndofs, &|v: &[f64]| stiffness(v).iter().zip(&kg(v)).map(|(k, g)| k + x * g).collect(), &precond, &|b: &mut Vec<f64>| sys.fold(b, 3));
        // Checked no further than `STABLE_TO` times the loads: past it the low end is that, and the
        // search's own is not said (at its own factor the operator is near singular along the mode
        // found, and the check then costs as much again as the solve).
        let mut x = low.min(STABLE_TO);
        low = 0.0;
        for _ in 0..10 {
            checks += 1;
            match stable(x) {
                Some(Ok(())) => {
                    low = x;
                    break;
                }
                // A direction the loads times x buckle: the lowest factor is under x.
                Some(Err(())) => factor = factor.min(x),
                // Not converged: nothing is known of its low end.
                None => {
                    unchecked = true;
                    break;
                }
            }
            x *= 0.5;
        }
    }
    let iterations = iterations + checks;
    let failure = failure.or_else(|| unchecked.then(|| "the check of the buckling factor's low end did not converge".to_string()));
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    match failure {
        Some(f) => (Buckle { antisymmetric: [false; 3], factor: f64::NAN, low: f64::NAN, at, iterations, ms, failure: Some(f) }, Vec::new()),
        None => (Buckle { antisymmetric: [false; 3], factor, low, at, iterations, ms, failure: None }, if none { Vec::new() } else { vec![shape] }),
    }
}

/// Where no cell's stress (`buckle`'s) is compressive, a low end of the
/// least buckling factor: infinite where none has any compression at all,
/// no direction the loads compress existing to buckle; finite where some
/// has a compression under roundoff's allowance. None where a cell is
/// compressed beyond it.
///
/// The solved stress of a part pulled without lateral restraint is exact
/// tension plus roundoff: its zero principal values read about 1e-15 of
/// the peak (a uniaxial bar on 1 mm cells: 1.7e-15). A compressive value
/// is ignored only under both of: 1e-12 of the peak principal stress
/// (hundreds of times that roundoff), and a thousandth of the least
/// stress that buckles any member the grid holds. That stress is bounded
/// by its slenderest: thickness t at least a hundredth of the finest cell,
/// length or width at most the grid's diagonal D, in the softest cell's
/// E. A strut's Euler stress pi^2 E / 12 (t / L)^2 at twice its length (a
/// cantilever), a plate's at k 0.425 (one edge free) and a thin section's
/// torsion (shear modulus E / 3) are all above 0.05 E (t / D)^2, so a
/// compression a thousandth of that, everywhere it acts, buckles nothing
/// below a thousand times the loads. A real compression, however local,
/// beyond it leaves the search to decide (and, unresolved, unknown).
fn uncompressed(sys: &System, stress: &[[f64; 6]], h: f64, pos: &[V3]) -> Option<f64> {
    let principal = |s: &[f64; 6]| dense_sym_eigen(&[s[0], s[3], s[5], s[3], s[1], s[4], s[5], s[4], s[2]], 3).0;
    let values: Vec<Vec<f64>> = stress.par_iter().map(principal).collect();
    let peak = values.iter().flatten().fold(0.0, |m: f64, v| m.max(v.abs()));
    if !peak.is_finite() {
        return None;
    }
    let soft = sys.cells.iter().map(|c| c.coef * c.mult).fold(f64::INFINITY, f64::min);
    let finest = h / f64::from(1u32 << sys.cells.iter().map(|c| c.level).max().unwrap_or(0));
    let (lo, hi) =
        pos.iter().fold(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]), |(lo, hi), p| (std::array::from_fn(|d| lo[d].min(p[d])), std::array::from_fn(|d| hi[d].max(p[d]))));
    let diagonal = (0..3).map(|d| (hi[d] - lo[d]).powi(2)).sum::<f64>().sqrt();
    let buckles = 0.05 * soft * (0.01 * finest / diagonal).powi(2);
    let allowance = (1e-12 * peak).min(1e-3 * buckles);
    // A NaN sorts last and `max` skips it: every value is asked finite.
    if !(soft > 0.0 && diagonal > 0.0 && allowance.is_finite() && values.iter().flatten().all(|v| v.is_finite())) {
        return None;
    }
    // The most compression ignored. None at all: nothing buckles at any
    // load. Some, under the allowance: only what the argument above gives,
    // nothing below (the least stress that buckles) / (that compression),
    // a thousand or more.
    let ignored = values.iter().map(|v| -v[0]).fold(0.0, f64::max);
    (ignored <= allowance).then(|| if ignored > 0.0 { buckles / ignored } else { f64::INFINITY })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `stable_at` takes no mode as the lowest: one buckling direction among 40,000, with no part in a
    /// smooth start, is met (Err), and with none it converges (Ok). The eigen-search's low end rests
    /// on the mode it found being the lowest; this is what checks it.
    #[test]
    fn a_single_buckling_direction_among_many_is_met() {
        let n = 40_000;
        // A stiffness of springs 1 to 100 on each unknown, and a load that softens one of them past
        // nothing (its own factor 0.5) and every other by a thousandth.
        let k = |i: usize| 1.0 + 99.0 * (i as f64 / n as f64);
        let g = |i: usize| if i == 31_337 { -2.0 * k(i) } else { -1e-3 * k(i) };
        let at = |x: f64| move |v: &[f64]| -> Vec<f64> { v.iter().enumerate().map(|(i, v)| (k(i) + x * g(i)) * v).collect() };
        // A preconditioner off by a factor that varies across the unknowns, as a multigrid cycle is.
        let t = |r: &[f64]| -> Vec<f64> { r.iter().enumerate().map(|(i, r)| r / (k(i) * (1.0 + 0.3 * ((i % 7) as f64 - 3.0) / 3.0))).collect() };
        assert_eq!(stable_at(n, &at(1.0), &t, &|_| ()), Some(Err(())), "the loads buckle one direction");
        assert_eq!(stable_at(n, &at(0.6), &t, &|_| ()), Some(Err(())), "and at 0.6 of them: its factor is 0.5");
        assert_eq!(stable_at(n, &at(0.45), &t, &|_| ()), Some(Ok(())), "and not at 0.45");
    }
}
