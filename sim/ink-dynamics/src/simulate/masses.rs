//! What the bodies' measured masses move a simulated motion by.
//!
//! Each body's mass, centre and inertia come from its ten measured moments,
//! each known only to its proven bounds (Body::measure): a source; so is
//! each carried body's mass per mm (Strip::err). The
//! motion is not linear in them, but near the measured masses it moves by
//! its first-order sensitivity to each, S_j, the change of the stepped
//! state y (positions, then speeds) were source j alone at the far side of
//! its bounds. S_j obeys the variational equations along the motion,
//! S_j' = J S_j + (0, g_j), J the motion's Jacobian there and g_j the
//! change source j makes to the accelerations at the state itself:
//! 1e6 M^-1 times its change of the equations' residual (Model::
//! mass_sources), widened by the couplings' worst loss as an effort's
//! band is. Each starts at 0 but a released joint's deviation, whose rest
//! moves by -K^-1 times what the source changes there (K the springs'
//! condensed stiffness). Each step carries S_j from knot to knot by the
//! exponential midpoint rule, S_j1 = e^X S_j0 + h phi1(X) gbar_j, with
//! X = h (J0 + J1) / 2 and gbar_j the two knots' forcing's mean: exact where the
//! Jacobian and the forcing hold over the step, second order otherwise,
//! and an undamped mode's sensitivity neither damped nor grown by a step
//! longer than its period. Between knots each S_j is read as the motion
//! is, by the quintic Hermite through its ends' changes of position, speed
//! and acceleration (J S_j + g_j), and the band of each stepped component
//! and of its first four derivatives at a moment is the sum of their sizes
//! there over the sources: what the moments anywhere in their bounds move
//! the motion by, independently, as an effort's band takes them
//! (Moment::mass), to first order. At a turning point a source that only
//! shifts the motion in time moves it by nothing.
//!
//! Sources whose forcing and change stay in proportion (the nuts and
//! screws riding one carriage, each moving it by its mass alone) are
//! carried as one group, its members' changes their shares of its own
//! (S_j = c_j S_g, so their sizes sum to S_g's times the sum of the
//! shares' sizes), until a member's forcing parts from the group's at a
//! knot; it is carried on its own from there (exact: the shares held over
//! every step before). A gantry's 1380 sources, 179 of them moving it,
//! are 14 to 16 groups. A knot where the sources or the Jacobian are not found
//! loses the band from there on (infinite). Note: the rule is exact where
//! J and g hold continuously over a step; sudden acceleration jumps (such as
//! a lossy coupling reversing at a turning point) are not jump-corrected.

use super::forces::Eval;
use super::linalg::{Lu, exp_phi};
use super::stepping::jacobian_at;
use super::trajectory::quintic;
use super::types::Layout;
use crate::model::Model;

/// Sources carried as one: its members (Model::mass_sources' order) and
/// each one's share of its change, and at the last knot its change of the
/// stepped state (positions, then speeds), of the accelerations, and its
/// forcing, the change it makes to the accelerations at the state alone.
struct Group {
    members: Vec<(usize, f64)>,
    s: Vec<f64>,
    a: Vec<f64>,
    g: Vec<f64>,
}

/// The motion's sensitivity to the masses at the last knot: its groups of
/// sources, and the motion's Jacobian there.
pub(super) struct Sensitivity {
    groups: Vec<Group>,
    jac: Vec<Vec<f64>>,
}

/// A step's band (Trajectory::mass_band): each group's quintic over it
/// per stepped component, and the sum of its members' shares' sizes.
pub(super) type Pieces = Vec<(f64, Vec<[f64; 6]>)>;

/// Each source's forcing at state y, moment t, read there as `ev`: 1e6
/// M^-1 times its change of the residual along the stepped components,
/// widened by Eval::lossy. Err: not found there, why.
fn forcing(model: &Model, lay: &Layout, t: f64, y: &[f64], ev: &Eval) -> Result<Vec<Vec<f64>>, String> {
    let (ns, n, nq) = (lay.ns, lay.n(), model.coords.len());
    let sim: Vec<[f64; 3]> = (0..ns).map(|i| [y[i], y[n + i], ev.acc[i]]).collect();
    let mut released = vec![0.0; lay.nr];
    for (j, &w) in lay.w.iter().enumerate() {
        released[w] = ev.acc[ns + j];
    }
    let rows: Vec<usize> = model.simulated.iter().map(|s| s.0).chain(lay.w.iter().map(|&w| nq + w)).collect();
    let sources = model.mass_sources(t, &model.beside(t, &sim)?, &released)?;
    // (A coordinate nothing with mass moves keeps its speed, and no source
    // changes what it needs: its inertia, singular, is not needed.)
    let mut lu = None;
    sources
        .iter()
        .map(|d| {
            let r: Vec<f64> = rows.iter().map(|&k| -1e6 * ev.lossy * d[k]).collect();
            if r.iter().any(|x| !x.is_finite()) {
                return Err("a moment of a body the simulation moves is not measured".to_string());
            }
            if r.iter().all(|x| *x == 0.0) {
                return Ok(r);
            }
            if lu.is_none() {
                lu = Some(Lu::of(ev.mss.clone()).ok_or("the stepped components' generalized inertia is singular there")?);
            }
            Ok(lu.as_ref().expect("factored").solve(&r))
        })
        .collect()
}

fn times(m: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    m.iter().map(|row| row.iter().zip(x).map(|(a, b)| a * b).sum()).collect()
}

/// `key`'s share of `lead`'s, where it is that share of it in every
/// component, to a billionth of each, or of a trillionth of the largest
/// in its block of `n` (one quantity: positions, speeds, accelerations),
/// rounding apart; a zero lead's every zero key at share 1.
fn share(key: &[f64], lead: &[f64], n: usize) -> Option<f64> {
    let (kk, kl) = (key.iter().zip(lead).map(|(a, b)| a * b).sum::<f64>(), lead.iter().map(|b| b * b).sum::<f64>());
    let c = if kl > 0.0 { kk / kl } else { 1.0 };
    key.chunks(n)
        .zip(lead.chunks(n))
        .all(|(ka, la)| {
            let most = ka.iter().zip(la).map(|(a, b)| a.abs().max((c * b).abs())).fold(0.0, f64::max);
            ka.iter().zip(la).all(|(a, b)| (a - c * b).abs() <= 1e-9 * (a.abs() + (c * b).abs()) + 1e-12 * most)
        })
        .then_some(c)
}

/// Sources gathered into groups by their keys (each one's change and
/// forcing, in blocks of n), each joining the first group whose lead's key
/// its own is a share of: per group, its members with their shares of the
/// lead's.
fn gather(keys: &[(usize, Vec<f64>)], n: usize) -> Vec<Vec<(usize, f64)>> {
    let mut groups: Vec<(usize, Vec<(usize, f64)>)> = Vec::new();
    for (k, (j, key)) in keys.iter().enumerate() {
        match groups.iter_mut().find_map(|(lead, members)| share(key, &keys[*lead].1, n).map(|c| (members, c))) {
            Some((members, c)) => members.push((*j, c)),
            None => groups.push((k, vec![(*j, 1.0)])),
        }
    }
    groups.into_iter().map(|(_, m)| m).collect()
}

/// Groups whose changes and forcing at the knot are in proportion, one:
/// each member's share of the lead's group's (exact from there while the
/// forcing stays in proportion, which each step checks).
fn merged(groups: Vec<Group>, n: usize) -> Vec<Group> {
    let keys: Vec<(usize, Vec<f64>)> = groups.iter().enumerate().map(|(k, gr)| (k, gr.s.iter().chain(&gr.a).chain(&gr.g).copied().collect())).collect();
    let mut slots: Vec<Option<Group>> = groups.into_iter().map(Some).collect();
    gather(&keys, n)
        .into_iter()
        .map(|same| {
            let mut lead = slots[same[0].0].take().expect("a group");
            for &(k, c) in &same[1..] {
                lead.members.extend(slots[k].take().expect("a group").members.into_iter().map(|(j, x)| (j, c * x)));
            }
            lead
        })
        .collect()
}

impl Sensitivity {
    /// At the start, state y read as `ev` at moment t, the deviations at
    /// rest where the springs hold the mechanism still (`still`, the
    /// coordinates' state there): each source's change of that rest, and
    /// none of the rest of the state. None: no source (every mass exact).
    pub(super) fn start(model: &Model, lay: &Layout, t: f64, y: &[f64], ev: &Eval, still: &[[f64; 3]]) -> Result<Option<Sensitivity>, String> {
        let (ns, n, nq) = (lay.ns, lay.n(), model.coords.len());
        let at_rest = model.mass_sources(t, still, &vec![0.0; lay.nr])?;
        if at_rest.is_empty() {
            return Ok(None);
        }
        let lu = (!lay.w.is_empty()).then(|| Lu::of(ev.keff.clone()).ok_or("the springs' stiffness is singular at the start")).transpose()?;
        let mut s = Vec::with_capacity(at_rest.len());
        for d in &at_rest {
            let mut sj = vec![0.0; 2 * n];
            if let Some(lu) = &lu {
                let moved = lu.solve(&lay.w.iter().map(|&w| -d[nq + w]).collect::<Vec<_>>());
                if moved.iter().any(|x| !x.is_finite()) {
                    return Err("a moment of a body the springs move is not measured".into());
                }
                sj[ns..n].copy_from_slice(&moved);
            }
            s.push(sj);
        }
        let jac = jacobian_at(model, lay, t, y, ev, &mut 0).map_err(|e| e.2)?;
        let g = forcing(model, lay, t, y, ev)?;
        let keys: Vec<(usize, Vec<f64>)> = s.iter().zip(&g).enumerate().map(|(j, (s, g))| (j, s.iter().chain(g).copied().collect())).collect();
        let groups = gather(&keys, n)
            .into_iter()
            .map(|members| {
                let lead = members[0].0;
                let a = times(&jac[n..], &s[lead]).iter().zip(&g[lead]).map(|(x, y)| x + y).collect();
                Group { members, s: s[lead].clone(), a, g: g[lead].clone() }
            })
            .collect();
        Ok(Some(Sensitivity { groups: merged(groups, n), jac }))
    }

    /// Carried over a step h long to the next knot, state y read there as
    /// `ev` at moment t: the step's band (Pieces). Err: not found there
    /// (the band lost), why.
    pub(super) fn step(&mut self, model: &Model, lay: &Layout, t: f64, h: f64, y: &[f64], ev: &Eval) -> Result<Pieces, String> {
        let n = lay.n();
        let jac = jacobian_at(model, lay, t, y, ev, &mut 0).map_err(|e| e.2)?;
        let g = forcing(model, lay, t, y, ev)?;
        // A member whose forcing parts from its share of its lead's is
        // carried on its own from the group's state at the last knot, with
        // those parting alike.
        let mut groups = Vec::with_capacity(self.groups.len());
        for gr in std::mem::take(&mut self.groups) {
            let lead = &g[gr.members[0].0];
            let (stay, part): (Vec<(usize, f64)>, Vec<(usize, f64)>) =
                gr.members.iter().partition(|&&(j, c)| share(&g[j], &lead.iter().map(|x| c * x).collect::<Vec<_>>(), n).is_some_and(|x| (x - 1.0).abs() <= 1e-9));
            let keys: Vec<(usize, Vec<f64>)> = part.iter().map(|&(j, c)| (j, gr.s.iter().map(|x| c * x).chain(g[j].iter().copied()).collect())).collect();
            for members in gather(&keys, n) {
                let c = part.iter().find(|m| m.0 == members[0].0).expect("a member parting").1;
                let at = |v: &[f64]| v.iter().map(|x| c * x).collect::<Vec<f64>>();
                groups.push(Group { members, s: at(&gr.s), a: at(&gr.a), g: at(&gr.g) });
            }
            if !stay.is_empty() {
                groups.push(Group { members: stay, ..gr });
            }
        }
        let x: Vec<Vec<f64>> = jac.iter().zip(&self.jac).map(|(r1, r0)| r1.iter().zip(r0).map(|(a, b)| 0.5 * h * (a + b)).collect()).collect();
        let (e, phi) = exp_phi(&x);
        // h phi1(X) along the speeds: what the mean forcing adds.
        let push: Vec<Vec<f64>> = phi.iter().map(|row| row[n..].iter().map(|v| h * v).collect()).collect();
        let mut pieces = Pieces::with_capacity(groups.len());
        for gr in &mut groups {
            let g1 = &g[gr.members[0].0];
            let mean: Vec<f64> = g1.iter().zip(&gr.g).map(|(a, b)| 0.5 * (a + b)).collect();
            let s1: Vec<f64> = times(&e, &gr.s).iter().zip(times(&push, &mean)).map(|(a, b)| a + b).collect();
            let a1: Vec<f64> = times(&jac[n..], &s1).iter().zip(g1).map(|(a, b)| a + b).collect();
            if s1.iter().chain(&a1).any(|x| !x.is_finite()) {
                return Err(format!("the motion's sensitivity to the masses is not finite at t = {t}"));
            }
            let each = (0..n).map(|i| quintic([gr.s[i], gr.s[n + i], gr.a[i]], [s1[i], s1[n + i], a1[i]], s1[i] - gr.s[i], h)).collect();
            pieces.push((gr.members.iter().map(|m| m.1.abs()).sum(), each));
            (gr.s, gr.a, gr.g) = (s1, a1, g1.clone());
        }
        (self.groups, self.jac) = (merged(groups, n), jac);
        Ok(pieces)
    }
}
