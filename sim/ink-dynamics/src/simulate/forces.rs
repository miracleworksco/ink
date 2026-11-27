use super::linalg::{Lu, solve};
use super::types::{Fail, Layout};
use crate::model::Model;
use crate::model::springs::At;

/// The accelerations at a state, and what the step reads of it: the
/// simulated coordinates' individual generalized inertias (a loop
/// folding at a dead point makes it grow without bound), the inertia along
/// the stepped components, the springs' stiffness along the stepped
/// deviations (those without inertia condensed out), each spring line's
/// tension, and each released joint's commanded value, rate and
/// acceleration.
#[derive(Clone)]
pub(super) struct Eval {
    pub(super) acc: Vec<f64>,
    pub(super) inertia: Vec<f64>,
    pub(super) mss: Vec<Vec<f64>>,
    pub(super) keff: Vec<Vec<f64>>,
    pub(super) two_sided: bool,
    pub(super) tensions: Vec<f64>,
    pub(super) drives: Vec<[f64; 3]>,
}

/// The simulated coordinates' accelerations at moment t, positions `q` and
/// speeds `v` (every released joint's deviation stepped after them), from
/// `guess` (the last ones). Err: why not (Fail::said).
pub fn accelerations(model: &Model, t: f64, q: &[f64], v: &[f64], guess: &[f64]) -> Result<Vec<f64>, String> {
    forces(model, &Layout::of(model), t, q, v, guess).map(|x| x.acc).map_err(|f| said(model, t, v, f))
}

/// What a failure says, at moment t and speeds v: at rest, no acceleration
/// agreeing with its losses is static friction holding it; moving, a
/// self-locking coupling driven back has a negative effective inertia.
pub(super) fn said(model: &Model, t: f64, v: &[f64], f: Fail) -> String {
    let names = model.simulated.iter().map(|&(k, _)| model.coords[k].name.as_str()).collect::<Vec<_>>().join(", ");
    match f {
        Fail::Refused(why) => why,
        Fail::Unbalanced { .. } if v.iter().all(|x| *x == 0.0) => format!(
            "{names} at rest at t = {t} is held by friction: no acceleration agrees with the direction its couplings' losses turn with (what drives it \
             lies within what their friction holds still, or a coupling that self-locks is driven back), and sticking is not modelled"
        ),
        Fail::Unbalanced { locked: true, .. } => format!(
            "a coupling along {names} that self-locks (its efficiency driven back, EB, at most 0: a lead screw or worm of efficiency at most 0.5) is \
             driven back at speed at t = {t}: its effective inertia there is negative (I + EB m r^2 below 0), so a rigid model's motion is not unique: \
             write the motion and read its effort (ink effort), or command it past what drives it back"
        ),
        Fail::Unbalanced { why, .. } => format!("no accelerations of {names} balance the forces at t = {t}: {why}"),
    }
}

/// The released joints carrying no inertia (`lay.y`) where the springs
/// leave them, the others' deviations in `d`: the springs' pull along each
/// balancing what it carries needs (`fr`: loads on it), by Newton's method
/// on their pull, piecewise linear (a line slack or taut), until the lines
/// taut stay so. Err: one nothing holds.
fn held(model: &Model, st: &At, fr: &[f64], lay: &Layout, d: &mut [f64]) -> Result<(), String> {
    if lay.y.is_empty() {
        return Ok(());
    }
    let taut = |d: &[f64]| -> Vec<bool> { (0..st.kl.len()).map(|l| st.t0[l] + st.kl[l] * st.g[l].iter().zip(d).map(|(g, x)| g * x).sum::<f64>() >= 0.0).collect() };
    for &y in &lay.y {
        d[y] = 0.0;
    }
    for _ in 0..16 {
        let before = taut(d);
        let fs = st.force(d);
        let k = st.stiffness(d);
        let kyy: Vec<Vec<f64>> = lay.y.iter().map(|&a| lay.y.iter().map(|&b| k[a][b]).collect()).collect();
        let ry: Vec<f64> = lay.y.iter().map(|&y| fs[y] - fr[y]).collect();
        let Some(step) = Lu::of(kyy).map(|lu| lu.solve(&ry)) else {
            let names: Vec<&str> = lay.y.iter().map(|&y| model.released()[y].0).collect();
            return Err(format!("{} carry no inertia along their own motion and no spring holds them there: their motion is not determined", names.join(", ")));
        };
        for (i, &y) in lay.y.iter().enumerate() {
            d[y] += step[i];
        }
        if taut(d) == before {
            return Ok(());
        }
    }
    Err("the springs' lines on the released joints carrying no inertia go slack and taut in turn: where they hold them is not found".into())
}

/// The springs' stiffness along the stepped deviations, those carrying no
/// inertia following (Guyan): K_ww - K_wy K_yy^-1 K_yw.
fn condensed(k: &[Vec<f64>], lay: &Layout) -> Vec<Vec<f64>> {
    let kww: Vec<Vec<f64>> = lay.w.iter().map(|&a| lay.w.iter().map(|&b| k[a][b]).collect()).collect();
    if lay.y.is_empty() {
        return kww;
    }
    let kyy: Vec<Vec<f64>> = lay.y.iter().map(|&a| lay.y.iter().map(|&b| k[a][b]).collect()).collect();
    let lu = Lu::of(kyy);
    let cols: Vec<Vec<f64>> =
        lay.w.iter().map(|&b| lu.as_ref().map_or_else(|| vec![0.0; lay.y.len()], |lu| lu.solve(&lay.y.iter().map(|&a| k[a][b]).collect::<Vec<_>>()))).collect();
    (0..lay.w.len()).map(|i| (0..lay.w.len()).map(|j| kww[i][j] - lay.y.iter().enumerate().map(|(c, &y)| k[lay.w[i]][y] * cols[j][c]).sum::<f64>()).collect()).collect()
}

/// The accelerations at a state (accelerations), and what the step reads of
/// it (Eval).
pub(super) fn forces(model: &Model, lay: &Layout, t: f64, q: &[f64], v: &[f64], guess: &[f64]) -> Result<Eval, Fail> {
    let (ns, n) = (lay.ns, lay.n());
    if q.iter().chain(v).any(|x| !x.is_finite()) {
        return Err(Fail::Refused(format!("its state is not finite at t = {t} (positions {q:?}, speeds {v:?})")));
    }
    let nq = model.coords.len();
    let sims = &model.simulated;
    // The rows the state's components read in the model's: a simulated
    // coordinate's own, a released joint's after the coordinates.
    let rows: Vec<usize> = sims.iter().map(|s| s.0).chain(lay.w.iter().map(|&j| nq + j)).collect();
    let mut d = vec![0.0; lay.nr];
    for (j, &w) in lay.w.iter().enumerate() {
        d[w] = q[ns + j];
    }
    // The residual at accelerations acc, each one's rounding, M, whether a
    // coupling that self-locks is driven back there, and the springs'.
    type Springy = (Vec<Vec<f64>>, bool, Vec<f64>, Vec<[f64; 3]>);
    let mut eval = |acc: &[f64]| -> Result<(Vec<f64>, Vec<f64>, Vec<Vec<f64>>, bool, Springy), Fail> {
        // At rest a coupling's loss turns with the direction the motion
        // starts in: read at the speed an instant later.
        let sim: Vec<[f64; 3]> = (0..ns).map(|i| [q[i], if v[i] == 0.0 { 1e-12 * acc[i] } else { v[i] }, acc[i]]).collect();
        let m = model.state_near(t, &model.beside(t, &sim).map_err(Fail::Refused)?);
        if let Some(why) = &m.unsolved {
            return Err(Fail::Refused(why.clone()));
        }
        if let Some(why) = model.forward_refusal(&m) {
            return Err(Fail::Refused(why));
        }
        // A site whose back efficiency is at most 0 (a screw that
        // self-locks) passing power back from what it drives.
        let locked = m.sites.iter().any(|s| s.eff <= 0.0 && s.total * s.rate < 0.0);
        // Each actuator's command, clamped to its limits at this state.
        let mut tau = vec![0.0; model.actuators.len()];
        for (a, act) in model.actuators.iter().enumerate() {
            if !sims.iter().any(|&(k, _)| m.bm[a][k] != 0.0) {
                continue;
            }
            // One moving a written coordinate too would have to supply what
            // that motion needs besides its command.
            if let Some(k) = (0..m.q.len()).find(|&k| m.bm[a][k] != 0.0 && model.coords[k].moves() && !sims.iter().any(|s| s.0 == k)) {
                return Err(Fail::Refused(format!(
                    "{} moves the simulated coordinates and {}, whose motion is written: simulate it too, or drive it by another actuator",
                    act.path, model.coords[k].name
                )));
            }
            let Some(c) = m.command[a] else { continue };
            if m.lo[a] > m.hi[a] {
                return Err(Fail::Refused(format!("{} gives no effort at its speed, {:.6}: its limits cross there (past its curve's last speed)", act.path, m.rate[a])));
            }
            if c.is_nan() {
                return Err(Fail::Refused(format!("the command on {} is not a number at t = {t}", act.path)));
            }
            tau[a] = c.clamp(m.lo[a], m.hi[a]);
            if !tau[a].is_finite() {
                return Err(Fail::Refused(format!(
                    "the command on {} is {c} at t = {t}, and its limits there, [{}, {}], do not bound it: no finite effort",
                    act.path, m.lo[a], m.hi[a]
                )));
            }
        }
        let (mut r, mut tol) = (vec![0.0; n], vec![0.0; n]);
        for (i, &(k, _)) in sims.iter().enumerate() {
            let q: f64 = (0..tau.len()).map(|a| m.bm[a][k] * tau[a]).sum();
            r[i] = q - m.f[k];
            tol[i] = 16.0 * m.rounding(k) + 1e-14 * q.abs();
        }
        // The springs' pull along each released joint carrying inertia,
        // those carrying none held where the springs leave them.
        let mut springy: Springy = (Vec::new(), true, Vec::new(), Vec::new());
        if lay.nr > 0 {
            let st = m.springs.as_ref().expect("the springs at a state of a model with springs");
            held(model, st, &m.fr, lay, &mut d).map_err(Fail::Refused)?;
            let fs = st.force(&d);
            for (j, &w) in lay.w.iter().enumerate() {
                r[ns + j] = fs[w] - m.fr[w];
                tol[ns + j] = 16.0 * m.noise_r[w] + 1e-14 * fs[w].abs();
            }
            let (tangent, two_sided) = st.tangent(&d);
            springy = (condensed(&tangent, lay), two_sided, st.tensions(&d), st.drive.clone());
        }
        // The deviations' accelerations move what they carry: no state of
        // the model reads them, so their inertia enters by M.
        for a in 0..n {
            let terms: Vec<f64> = lay.w.iter().enumerate().map(|(j, &w)| m.gm[rows[a]][nq + w] * acc[ns + j]).collect();
            r[a] -= 1e-6 * terms.iter().sum::<f64>();
            tol[a] += 1e-15 * terms.iter().map(|x| x.abs()).sum::<f64>();
        }
        if r.iter().any(|x| !x.is_finite()) {
            return Err(Fail::Refused(format!(
                "the forces at t = {t} are not finite (the efforts needed: {:?})",
                rows.iter().map(|&k| if k < nq { m.f[k] } else { m.fr[k - nq] }).collect::<Vec<_>>()
            )));
        }
        let mss = rows.iter().map(|&k| rows.iter().map(|&l| m.gm[k][l]).collect()).collect();
        Ok((r, tol, mss, locked, springy))
    };
    let (mut r, tol0, mss, mut locked, (keff, two_sided, tensions, drives)) = eval(guess)?;
    let inertia: Vec<_> = (0..ns).map(|i| mss[i][i]).collect();
    let done = |acc: Vec<f64>| Eval { acc, inertia: inertia.clone(), mss: mss.clone(), keff: keff.clone(), two_sided, tensions: tensions.clone(), drives: drives.clone() };
    let mut acc = guess.to_vec();
    // M's step, taken from any guess (one within the rounding kept as it is
    // would make the accelerations stick to the last state's, which the
    // integrator's Newton iterations and differences read as no change):
    // exact where the effort is affine in the accelerations.
    // Nothing with mass along them and nothing pushing them: their
    // accelerations are not determined, and they keep their speed.
    let pushed = !r.iter().zip(&tol0).all(|(x, e)| x.abs() <= *e);
    let step = solve(&mss, &r.iter().map(|x| 1e6 * x).collect::<Vec<_>>()).or_else(|| (!pushed).then(|| vec![0.0; n])).ok_or_else(|| {
        let names: Vec<&str> = sims.iter().map(|&(k, _)| model.coords[k].name.as_str()).chain(lay.w.iter().map(|&w| model.released()[w].0)).collect();
        let springs = if lay.nr > 0 {
            "; a simulated coordinate whose every body rides joints released for springs moves nothing of its own: write its motion, which commands them"
        } else {
            ""
        };
        Fail::Refused(format!("nothing with mass moves along {} at t = {t} (no inertia{springs}), so their accelerations are not determined", names.join(", ")))
    })?;
    for i in 0..n {
        acc[i] += step[i];
    }
    // With no simulated coordinate, no state of the model reads the
    // accelerations (the deviations' enter by M alone): the step is exact.
    if ns == 0 {
        return Ok(done(acc));
    }
    let mut last = step;
    // Whether a self-locking coupling was driven back on the way: where no
    // root is found, that is why.
    let mut ever = locked;
    for _ in 0..30 {
        let (r1, tol, _, l, _) = eval(&acc)?;
        (r, locked) = (r1, l);
        ever |= l;
        if r.iter().zip(&tol).all(|(x, e)| x.abs() <= *e) {
            return if locked { Err(Fail::Unbalanced { locked, why: String::new() }) } else { Ok(done(acc)) };
        }
        // r's own Jacobian, by differences (exact on a branch of a loss
        // that turns with the power's direction), and Newton's step.
        let mut jac = vec![vec![0.0; n]; n];
        for j in 0..n {
            let d = 1e-4 * acc[j].abs().max(last[j].abs()).max(1.0);
            let mut moved = acc.clone();
            moved[j] += d;
            let (rj, _, _, _, _) = eval(&moved)?;
            for i in 0..n {
                jac[i][j] = (rj[i] - r[i]) / d;
            }
        }
        let Some(step) = Lu::of(jac).map(|lu| lu.solve(&r.iter().map(|x| -x).collect::<Vec<_>>())) else {
            return Err(Fail::Unbalanced {
                locked: ever,
                why: "the efforts do not change with them there (a coupling's loss making them independent of the accelerations)".into(),
            });
        };
        for i in 0..n {
            acc[i] += step[i];
        }
        last = step;
    }
    Err(Fail::Unbalanced {
        locked: ever,
        why: format!("Newton's method on them, the couplings' losses turning with the power's direction, did not converge in 30 steps (what is left: {r:?})"),
    })
}
