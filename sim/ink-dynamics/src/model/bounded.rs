//! An actuator's effort bounded over a step of the motion, not read at its
//! moments: inverse dynamics in outward-rounded interval arithmetic over
//! every moment of the step at once (ink_core::interval).
//!
//! Each body's share is the virtual work of its mass moving with its
//! placement x = R u + p (u a point of it, in its frame): along coordinate
//! k, the integral of (dR_k u + dp_k) . (R'' u + p'' - g) over its mass,
//! which is linear in its ten moments, its mass m, first moment S = m c and
//! second moment J = (integral of u u^T), and reads only its placement's
//! second rates along t and its partials along the coordinates:
//!
//!   f_k = 1e-6 sum over world axes l of
//!         dR_k[l] . (J R''[l]) + (dR_k[l] . S) a_l + dp_k[l] (R''[l] . S) + m dp_k[l] a_l,
//!   a_l = p''[l] - g[l].
//!
//! (For a rotation this is Kane's m (a - g) . J_v + (I alpha + omega x I
//! omega) . J_omega, which Model::compute sums at a moment.) A load's share
//! is its force and torque through its point's partial velocity, as at a
//! moment. Every number of the model is read over the step as its bound
//! there, its bound at the step's middle and its rate's bound over the step
//! (Over), so each sum is the meet of its plain bound and its mean value
//! form, which shrinks with the step's square. A body's measured moments,
//! each anywhere in its proven bounds (Body::measure), move the sums by
//! their changes' bounds over the step, as Moment::mass sums them at a
//! moment. The split between actuators is the drives' constant ratios
//! inverted in intervals (a block with as many drives as coordinates, whose
//! split rests on no limit), and the limits are bounded over the step, so a
//! limit that follows the drive's speed (a motor's torque-speed line) is
//! read over every speed the step passes.
//!
//! Only what this covers is bounded (Model::bounded): a mechanism moving
//! through no solve (no loop closed at each moment, nothing simulated or
//! sprung), with no coupling and no carried material, its drives' ratios
//! constant, every body's measures known. Anything else is read at its
//! moments as before.

use super::{Model, Num, ROUNDING, Sets};
use crate::steps::MAX_STEPS;
use ink_core::eval::{Ctx, Inputs, eval_interval, eval_interval_rates, span_mid};
use ink_core::interval::{I, RAD, add_o, add_up, div_o, hull, mag, meet, mul_o, mul_up, sub_o};
use ink_lang::Steps;
use std::sync::Arc;

const ZERO: I = I { lo: 0.0, hi: 0.0 };

/// A number over a step of the motion: its bound over the step, its bound
/// at the step's middle, and its rate's bound over the step (ENTIRE where
/// it may jump there). Sums and products carry each (the rate by the
/// product rule over the step's bounds), so a result's rate bounds its
/// change from its middle by the mean value theorem.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Over {
    v: I,
    m: I,
    d: I,
}

impl Over {
    const ZERO: Over = Over { v: ZERO, m: ZERO, d: ZERO };

    /// A number that does not move over the step.
    fn exact(x: I) -> Over {
        Over { v: x, m: x, d: ZERO }
    }

    fn add(self, o: Over) -> Over {
        Over { v: add_o(self.v, o.v), m: add_o(self.m, o.m), d: add_o(self.d, o.d) }
    }

    fn sub(self, o: Over) -> Over {
        Over { v: sub_o(self.v, o.v), m: sub_o(self.m, o.m), d: sub_o(self.d, o.d) }
    }

    fn mul(self, o: Over) -> Over {
        Over { v: mul_o(self.v, o.v), m: mul_o(self.m, o.m), d: add_o(mul_o(self.d, o.v), mul_o(self.v, o.d)) }
    }

    /// Times a number that does not move over the step.
    fn times(self, k: I) -> Over {
        Over { v: mul_o(self.v, k), m: mul_o(self.m, k), d: mul_o(self.d, k) }
    }

    /// Its bound over a step reaching `r` either side of its middle: the
    /// plain bound met with the middle's plus the rate's over the step times
    /// how far from the middle (both hold every value; so does their meet).
    fn bound(self, r: f64) -> I {
        meet(self.v, add_o(self.m, mul_o(self.d, I::new(-r, r))))
    }
}

fn dot(a: [Over; 3], b: [Over; 3]) -> Over {
    a[0].mul(b[0]).add(a[1].mul(b[1])).add(a[2].mul(b[2]))
}

fn dot_i(a: [I; 3], b: [I; 3]) -> I {
    add_o(add_o(mul_o(a[0], b[0]), mul_o(a[1], b[1])), mul_o(a[2], b[2]))
}

/// A model's number over the step `over` spans, its middle `at` (Over).
/// None through a solve: what it sets has no bound over a step here.
fn read(n: &Num, over: &Ctx, at: &Ctx) -> Option<Over> {
    match n {
        Num::C(c) => Some(Over::exact(I::point(*c))),
        Num::T(t) => {
            let origin = [I::point(0.0); 3];
            let (v, d) = eval_interval_rates(t, over, origin, &Inputs { t: true, params: Vec::new() });
            let m = eval_interval(t, at, origin, None);
            Some(Over { v, m, d: d.first().copied().unwrap_or(I::ENTIRE) })
        }
        Num::S(..) => None,
    }
}

/// Whether a number is continuous in time over `step` (no jump in position
/// or first rate). `eval_interval_rates_jump` reports 0 jump AND finite
/// rates on smooth expressions (including trig/rotations) and non-straddling
/// selects (note a non-monotone straddled select isn't split and gives jump
/// 0 with unbounded rates, so finiteness is required too); exact rational
/// `continuous_in_time` covers piecewise-continuous joins like a trapezoid's
/// velocity corner. A velocity jump inside a placement (`abs`, `max`,
/// `select`) fails both, so symbolic `dt2` cannot drop an impulse.
fn smooth(n: &Num, over: &Ctx, base: &Ctx, step: (f64, f64)) -> bool {
    match n {
        Num::C(_) => true,
        Num::T(t) => {
            let origin = [I::point(0.0); 3];
            let (_, d, jump) = ink_core::eval::eval_interval_rates_jump(t, over, origin, &Inputs { t: true, params: Vec::new() });
            (jump == 0.0 && d.iter().all(|x| x.lo.is_finite() && x.hi.is_finite())) || ink_core::eval::continuous_in_time(t, base, origin, step)
        }
        Num::S(..) => false,
    }
}

fn read_all<const N: usize>(ns: &[Num], over: &Ctx, at: &Ctx) -> Option<[Over; N]> {
    let v: Vec<Over> = ns.iter().map(|n| read(n, over, at)).collect::<Option<_>>()?;
    v.try_into().ok()
}

/// 10^-6 exactly: g mm/s^2 and g mm^2/s^2 are that of N and N mm.
fn micro() -> I {
    I::new(1e-6f64.next_down(), 1e-6f64.next_up())
}

/// A body's ten moments from its mass (g), centre (mm) and inertia about
/// that centre (g mm^2), each in its frame: m, S = m c, and the second
/// moment about its frame's origin J = C + m c c^T, C = (tr I / 2) 1 - I
/// (each product of inertia read either way it is written).
type Moments = (I, [I; 3], [[I; 3]; 3]);

fn moments(m: I, c: [I; 3], inertia: [[I; 3]; 3]) -> Moments {
    let s = c.map(|x| mul_o(m, x));
    let half = mul_o(I::point(0.5), add_o(add_o(inertia[0][0], inertia[1][1]), inertia[2][2]));
    let j = std::array::from_fn(|i| {
        std::array::from_fn(|k| {
            let own = if i == k { half } else { ZERO };
            add_o(sub_o(own, hull(inertia[i][k], inertia[k][i])), mul_o(m, mul_o(c[i], c[k])))
        })
    });
    (m, s, j)
}

fn point3(v: [f64; 3]) -> [I; 3] {
    v.map(I::point)
}

/// What the bodies and loads sum into over a step along each coordinate
/// (N mm per unit of it): the effort, as Over; how far its bodies'
/// measured moments can move it, each anywhere in its proven bounds,
/// summed (as Moment::mass sums them at a moment); and its rounding's
/// scale, ROUNDING of the magnitudes summed (Moment::noise_f's), within
/// which a moment's float reading may lie off its real value.
pub(crate) struct Swept {
    f: Vec<Over>,
    mass: Vec<f64>,
    noise: Vec<f64>,
}

/// What actuator a's effort is bounded with (Model::bounded): its share of
/// each coordinate's effort (the drives' constant ratios inverted, in
/// intervals), and the actuators whose limits crossing would refuse it (a
/// first, then those sharing a coordinate with it).
pub(crate) struct Plan {
    row: Vec<(usize, I)>,
    limits: Vec<usize>,
}

/// A step judged (Model::judge).
#[derive(Clone, Copy, Debug, PartialEq)]
enum Judged {
    /// The effort's bound over the whole step lies inside its limits there.
    Proven,
    /// The claim fails at the step's middle, as its claim reads it there.
    Fails,
    /// Not even at its middle is the bound inside: what the masses are
    /// known to (or the rounding) straddles a limit there, so no halving
    /// settles it.
    Short,
    /// The bound over the step crosses a limit; halved, it may not.
    Open,
}

/// An effort's steps bounded over a span (Model::bound_steps): the steps,
/// each a step the bound settled halved as it was; the steps of those given
/// it did not settle, read at their moments as before; how many steps it
/// proved, and the halvings; and whether the claim fails at a moment read.
pub(crate) struct Bounded {
    pub steps: Steps,
    pub unsettled: Vec<(f64, f64)>,
    pub proven: usize,
    pub halved: usize,
    pub failed: bool,
}

impl Model {
    /// Whether actuator a's effort is bounded over a step here, and with
    /// what (Plan); Err: why not, what lies outside what the bound covers.
    pub(crate) fn bounded(&self, a: usize) -> Result<Plan, String> {
        let act = &self.actuators[a];
        if self.springs.is_some() || !self.simulated.is_empty() || self.states.is_some() {
            return Err("its motion is simulated, not written".into());
        }
        if !self.solves.is_empty() || self.reads_moment {
            return Err("its mechanism moves through a solve made again at each moment (a loop's closure), which has no bound over a step here".into());
        }
        if !self.couplings.is_empty() || !self.balances.is_empty() {
            return Err("couplings' tensions are a solve whose coefficients move, not bounded over a step here".into());
        }
        if let Some(s) = self.strips.first() {
            return Err(format!("{}'s material moves round its path, not bounded over a step here", s.label));
        }
        if act.coupled.is_empty() || self.commanded(a) || !act.leaves_out.is_empty() {
            return Err("its effort is not the motion's inverse dynamics".into());
        }
        if matches!((&act.lo, &act.hi), (Num::C(lo), Num::C(hi)) if lo == hi) {
            return Err("its rating is one number, held only to its rounding".into());
        }
        if let Some(b) = self.bodies.iter().find(|b| b.unknown != [0.0; 2]) {
            return Err(format!("{}'s departure from round about its spin is not known within its measure", b.label));
        }
        if let Some(b) = self.bodies.iter().find(|b| b.measure.iter().any(|(dm, dc, di)| !dm.is_finite() || !dc.iter().chain(di.iter().flatten()).all(|x| x.is_finite()))) {
            return Err(format!("what {}'s mass is known to is not", b.label));
        }
        let nq = self.coords.len();
        // Each actuator's drive along each coordinate, radians for a turn: constant, or the split moves.
        let mut bm: Vec<Vec<I>> = Vec::with_capacity(self.actuators.len());
        for x in &self.actuators {
            let scale = if x.turn { RAD } else { I::point(1.0) };
            let row = x.along.iter().map(|n| match n {
                Num::C(c) => Some(mul_o(scale, I::point(*c))),
                _ => None,
            });
            bm.push(row.collect::<Option<_>>().ok_or("a drive's ratio to its coordinates moves along the motion")?);
        }
        let moves = |d: usize, k: usize| bm[d][k] != ZERO;
        let mut sets = Sets((0..nq).collect());
        for d in 0..bm.len() {
            let axes: Vec<usize> = (0..nq).filter(|&k| moves(d, k)).collect();
            sets.join(&axes);
        }
        let root: Vec<usize> = (0..nq).map(|k| sets.find(k)).collect();
        // A block: coordinates the drives join, and those drives. Each a
        // coordinate the effort rests on is in must split its effort by the
        // drives' ratios alone: as many drives as coordinates, independent.
        let block = |r: usize| -> (Vec<usize>, Vec<usize>) {
            let axes: Vec<usize> = (0..nq).filter(|&k| root[k] == r).collect();
            let drives: Vec<usize> = (0..bm.len()).filter(|&d| axes.iter().any(|&k| moves(d, k))).collect();
            (axes, drives)
        };
        let inverse = |r: usize| -> Result<(Vec<usize>, Vec<usize>, Vec<Vec<I>>), String> {
            let (axes, drives) = block(r);
            if drives.len() != axes.len() {
                return Err(format!("{} drives share {} coordinates, so their split rests on what each can give, not on their ratios alone", drives.len(), axes.len()));
            }
            // Coordinate k's balance: sum over drives d of B[d][k] tau_d = f_k.
            let m: Vec<Vec<I>> = axes.iter().map(|&k| drives.iter().map(|&d| bm[d][k]).collect()).collect();
            let inv = invert(m).ok_or("its drives' directions are too close to separate")?;
            Ok((axes, drives, inv))
        };
        let own = |b: usize| (0..nq).filter(|&k| moves(b, k)).collect::<Vec<usize>>();
        let mine = own(a);
        let first = *mine.first().ok_or("its drive moves no coordinate")?;
        for &k in act.coupled.iter().filter(|&&k| self.needs_drive(k)) {
            if block(root[k]).1.is_empty() {
                return Err(format!("{} moves and no actuator drives it", self.coords[k].name));
            }
            inverse(root[k])?;
        }
        let (axes, drives, inv) = inverse(root[first])?;
        let i = drives.iter().position(|&d| d == a).expect("its own block");
        let row = axes.iter().zip(&inv[i]).map(|(&k, &p)| (k, p)).collect();
        let limits = std::iter::once(a).chain((0..self.actuators.len()).filter(|&b| b != a && own(b).iter().any(|k| mine.contains(k)))).collect();
        Ok(Plan { row, limits })
    }

    /// What the bodies and loads sum into over `step` (Swept), once per
    /// step. None: something it reads has no bound over a step.
    fn swept(&self, step: (f64, f64)) -> Option<Arc<Swept>> {
        let key = (step.0.to_bits(), step.1.to_bits());
        if let Some(s) = self.swept.lock().unwrap().get(&key) {
            return s.clone();
        }
        let s = self.sweep(step).map(Arc::new);
        self.swept.lock().unwrap().insert(key, s.clone());
        s
    }

    fn sweep(&self, step: (f64, f64)) -> Option<Swept> {
        let mid = span_mid(step);
        let over = Ctx { t: mid, params: self.params.clone(), t_span: Some(step), p_spans: Default::default() };
        let at = Ctx { t_span: Some((mid, mid)), ..over.clone() };
        let nq = self.coords.len();
        let (mut f, mut mass, mut scale) = (vec![Over::ZERO; nq], vec![0.0; nq], vec![0.0; nq]);
        let micro = micro();
        let base_ctx = Ctx::new(0.0, self.params.clone());
        let gv: [I; 3] = self.gravity.map(|g| mul_o(I::point(1000.0), I::point(g)));
        for b in self.bodies.iter().filter(|b| !b.at.dq.is_empty()) {
            // A placement or its first rate jumping inside the step implies an
            // impulse that symbolic dt2 drops; leave the step unsettled.
            if !b.at.r.iter().chain(&b.at.p).chain(&b.at.dr).chain(&b.at.dp).all(|n| smooth(n, &over, &base_ctx, step)) {
                return None;
            }
            let dt2: [Over; 12] = read_all(&b.at.dt2, &over, &at)?;
            let rdd = |l: usize| [dt2[3 * l], dt2[3 * l + 1], dt2[3 * l + 2]];
            let acc: [Over; 3] = std::array::from_fn(|l| dt2[9 + l].sub(Over::exact(gv[l])));
            let (m, s, j) = moments(I::point(b.mass), point3(b.center), b.inertia.map(point3));
            // Per world axis l: J R''[l] and R''[l] . S.
            let jr: [[Over; 3]; 3] = std::array::from_fn(|l| std::array::from_fn(|i| dot(rdd(l), j[i].map(Over::exact))));
            let rs: [Over; 3] = std::array::from_fn(|l| dot(rdd(l), s.map(Over::exact)));
            // Each measured moment's change at the far side of its bounds, as its ten moments move.
            let base = (m, s, j);
            let changes: Vec<Moments> = b
                .measure
                .iter()
                .map(|(dm, dc, di)| {
                    let moved = moments(
                        add_o(I::point(b.mass), I::point(*dm)),
                        std::array::from_fn(|k| add_o(I::point(b.center[k]), I::point(dc[k]))),
                        std::array::from_fn(|r| std::array::from_fn(|c| add_o(I::point(b.inertia[r][c]), I::point(di[r][c])))),
                    );
                    (I::point(*dm), std::array::from_fn(|k| sub_o(moved.1[k], base.1[k])), std::array::from_fn(|r| std::array::from_fn(|c| sub_o(moved.2[r][c], base.2[r][c]))))
                })
                .collect();
            for (k, d) in &b.at.dq {
                let d: [Over; 12] = read_all(d, &over, &at)?;
                let (mut sum, mut size) = (Over::ZERO, 0.0);
                for l in 0..3 {
                    let (dr, dp) = ([d[3 * l], d[3 * l + 1], d[3 * l + 2]], d[9 + l]);
                    let terms = [dot(dr, jr[l]), dot(dr, s.map(Over::exact)).mul(acc[l]), dp.mul(rs[l]), dp.times(m).mul(acc[l])];
                    for x in terms {
                        sum = sum.add(x);
                        size = add_up(size, mag(x.v));
                    }
                }
                f[*k] = f[*k].add(sum.times(micro));
                scale[*k] = add_up(scale[*k], mul_up(size, micro.hi));
                // Linear in the moments: each change's share over the step, its bound's size.
                for (dm, ds, dj) in &changes {
                    let mut x = ZERO;
                    for l in 0..3 {
                        let (dr, dp, r) = ([d[3 * l].v, d[3 * l + 1].v, d[3 * l + 2].v], d[9 + l].v, rdd(l).map(|o| o.v));
                        let djr: [I; 3] = std::array::from_fn(|i| dot_i(r, dj[i]));
                        let a = acc[l].v;
                        for y in [dot_i(dr, djr), mul_o(dot_i(dr, *ds), a), mul_o(dp, dot_i(r, *ds)), mul_o(mul_o(*dm, dp), a)] {
                            x = add_o(x, y);
                        }
                    }
                    mass[*k] = add_up(mass[*k], mul_up(mag(x), micro.hi));
                }
            }
        }
        // A load: its force and torque in the world's axes at its point x,
        // which moves with what it acts on (placed R p + T): through the
        // point's partial velocity dR R^T (x - T) + dT and the angular one,
        // vee(dR R^T), an effort it supplies is one the actuators need not.
        for l in self.loads.iter().filter(|l| !l.dq.is_empty()) {
            if !l.at.iter().chain(&l.on).all(|n| smooth(n, &over, &base_ctx, step)) {
                return None;
            }
            let (x, force, torque): ([Over; 3], [Over; 3], [Over; 3]) = (read_all(&l.at, &over, &at)?, read_all(&l.force, &over, &at)?, read_all(&l.torque, &over, &at)?);
            let on: [Over; 12] = read_all(&l.on, &over, &at)?;
            let arm: [Over; 3] = std::array::from_fn(|i| x[i].sub(on[9 + i]));
            for (k, d) in &l.dq {
                let d: [Over; 12] = read_all(d, &over, &at)?;
                // W = dR R^T, row by row.
                let w: [[Over; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| dot([d[3 * i], d[3 * i + 1], d[3 * i + 2]], [on[3 * j], on[3 * j + 1], on[3 * j + 2]])));
                let jv: [Over; 3] = std::array::from_fn(|i| dot(w[i], arm).add(d[9 + i]));
                let half = I::point(0.5);
                let jw = [w[2][1].sub(w[1][2]).times(half), w[0][2].sub(w[2][0]).times(half), w[1][0].sub(w[0][1]).times(half)];
                let mut size = 0.0;
                for i in 0..3 {
                    for x in [jv[i].mul(force[i]), jw[i].mul(torque[i])] {
                        f[*k] = f[*k].sub(x);
                        size = add_up(size, mag(x.v));
                    }
                }
                scale[*k] = add_up(scale[*k], size);
            }
        }
        let noise = scale.iter().map(|s| mul_up(ROUNDING, *s)).collect();
        Some(Swept { f, mass, noise })
    }

    /// Actuator a's effort over `step` as Plan splits it: its bound over the
    /// step (the motion's alone, its masses as measured, and the rounding a
    /// moment's reading may carry), its bound at the step's middle, and how
    /// far its bodies' measured moments can move it (Moment::mass over the
    /// step). None: it rests on what has no bound over a step.
    pub(crate) fn effort_over(&self, plan: &Plan, step: (f64, f64)) -> Option<(I, I, f64)> {
        let sw = self.swept(step)?;
        let tau = plan.row.iter().fold(Over::ZERO, |x, &(k, p)| x.add(sw.f[k].times(p)));
        let noise = plan.row.iter().fold(mul_up(ROUNDING, mag(tau.v)), |x, &(k, p)| add_up(x, mul_up(mag(p), sw.noise[k])));
        let mass = plan.row.iter().fold(0.0, |x, &(k, p)| add_up(x, mul_up(mag(p), sw.mass[k])));
        let mid = span_mid(step);
        let r = (mid - step.0).max(step.1 - mid).next_up();
        let wide = |x: I| add_o(x, I::new(-noise, noise));
        Some((wide(tau.bound(r)), wide(tau.m), mass))
    }

    /// Whether every actuator of `plan.limits` has its limits apart over the
    /// context's span (none crossing: no effort lies within crossed ones),
    /// and actuator a's lie either side of `e` there.
    fn inside(&self, a: usize, plan: &Plan, e: I, c: &Ctx) -> bool {
        let apart = plan.limits.iter().all(|&b| self.actuators[b].lo.over(c).hi <= self.actuators[b].hi.over(c).lo);
        let act = &self.actuators[a];
        apart && e.lo > act.lo.over(c).hi && e.hi < act.hi.over(c).lo
    }

    /// Step `step` of actuator a's claim judged by its effort's bound over
    /// it (Judged), and at its middle as the claim reads it there.
    fn judge(&self, a: usize, plan: &Plan, step: (f64, f64)) -> Judged {
        let mid = span_mid(step);
        // Its claim at the middle, as it reads it (Effort::provide, its band
        // Moment::mass): past a limit by more than its masses' band and the
        // bodies read round's bound, a witness.
        let fails = || {
            let m = self.read(mid);
            m.unsolved.is_none() && self.refusal(a, &m).is_none() && m.margin(a) < -(m.mass[a] + m.bound[a])
        };
        let Some((over, at, mass)) = self.effort_over(plan, step) else {
            return if fails() { Judged::Fails } else { Judged::Short };
        };
        let band = |x: I| add_o(x, I::new(-mass, mass));
        let ctx = Ctx { t: mid, params: self.params.clone(), t_span: Some(step), p_spans: Default::default() };
        if self.inside(a, plan, band(over), &ctx) {
            return Judged::Proven;
        }
        if fails() {
            return Judged::Fails;
        }
        if !self.inside(a, plan, band(at), &Ctx { t_span: Some((mid, mid)), ..ctx }) {
            return Judged::Short;
        }
        Judged::Open
    }

    /// `steps` (Model::steps over a span) with actuator a's effort bounded
    /// over each step (Plan), each step whose bound crosses a limit halved
    /// while that may settle it, at most MAX_STEPS halvings in all. A step
    /// is settled where every piece it was halved into is proven (Judged::
    /// Proven), or one fails at its middle: its pieces are steps then, each
    /// read. One not settled (a piece crossing at the budget, or at its
    /// middle within what the masses are known to of a limit) is kept whole,
    /// as given, its halvings undone, for its moments to decide as before.
    /// Where one fails, the halving stops: the claim fails, and nothing is
    /// loose. Unchecked steps stay so.
    pub(crate) fn bound_steps(&self, a: usize, plan: &Plan, steps: &Steps) -> Bounded {
        use rayon::prelude::*;
        let roots: Vec<(f64, f64)> = steps.ends.windows(2).map(|w| (w[0], w[1])).filter(|s| !steps.unchecked.contains(s)).collect();
        let mut leaves: Vec<Vec<((f64, f64), Judged)>> = vec![Vec::new(); roots.len()];
        let mut todo: Vec<((f64, f64), usize)> = roots.iter().copied().zip(0..).collect();
        let (mut room, mut failed, mut halved) = (MAX_STEPS, false, 0);
        // Do not stop halving on a single limit crossing: a claim with
        // `(accuracy x)` may still hold at that crossing, so other open steps
        // must keep halving down to any real bump beyond accuracy (`MAX_STEPS`
        // caps total work).
        while !todo.is_empty() && !ink_kernel::stop::reached() {
            let judged: Vec<Judged> = todo.par_iter().map(|&(s, _)| self.judge(a, plan, s)).collect();
            let mut next = Vec::new();
            for ((s, i), j) in todo.drain(..).zip(judged) {
                let m = span_mid(s);
                if j == Judged::Open && room > 0 && s.0 < m && m < s.1 {
                    next.extend([((s.0, m), i), ((m, s.1), i)]);
                    (room, halved) = (room - 1, halved + 1);
                } else {
                    failed |= j == Judged::Fails;
                    leaves[i].push((s, j));
                }
            }
            todo = next;
        }
        // Pieces left when it failed or stopped: not settled.
        for (s, i) in todo {
            leaves[i].push((s, Judged::Open));
        }
        let mut ends = steps.ends.clone();
        let (mut unsettled, mut proven) = (Vec::new(), 0);
        for (&root, pieces) in roots.iter().zip(&leaves) {
            if pieces.iter().all(|p| p.1 == Judged::Proven) {
                ends.extend(pieces.iter().flat_map(|p| [p.0.0, p.0.1]));
                proven += pieces.len();
            } else if pieces.iter().any(|p| p.1 == Judged::Fails) {
                ends.extend(pieces.iter().flat_map(|p| [p.0.0, p.0.1]));
                for &(s, j) in pieces {
                    if j == Judged::Proven {
                        proven += 1;
                    } else if j != Judged::Fails {
                        unsettled.push(s);
                    }
                }
            } else {
                unsettled.push(root);
            }
        }
        ends.sort_by(f64::total_cmp);
        ends.dedup();
        // Emit only unsettled pieces lying inside a loose root so `loose`
        // matches actual windows of `ends` and never swallows Proven/Fails pieces.
        let loose = unsettled.iter().copied().filter(|u| steps.loose.iter().any(|s| u.0 >= s.0 && u.1 <= s.1)).collect();
        Bounded { steps: Steps { ends, unchecked: steps.unchecked.clone(), loose, basis: steps.basis.clone() }, unsettled, proven, halved, failed }
    }

    /// Actuator a's effort bounded over `step` where Model::bounded covers
    /// it: its bound over the step with its masses as measured and the
    /// rounding a moment's reading carries, and what its masses' measures
    /// can move it by. For tests reading the bound against dense moments.
    #[cfg(test)]
    pub(crate) fn effort_bound(&self, a: usize, step: (f64, f64)) -> Result<(I, f64), String> {
        let plan = self.bounded(a)?;
        let (over, _, mass) = self.effort_over(&plan, step).ok_or("no bound over the step")?;
        Ok((over, mass))
    }
}

/// The inverse of a square matrix of intervals, every matrix within them
/// inverted (Gauss-Jordan, outward rounded, the pivot of least magnitude
/// largest in its column). None where a pivot may be 0, or is under 1e-9 of
/// the matrix's largest entry (directions too close to separate, which the
/// split at a moment takes as shared: drives::effort_split).
fn invert(m: Vec<Vec<I>>) -> Option<Vec<Vec<I>>> {
    let n = m.len();
    let least = |x: I| {
        if x.lo > 0.0 {
            x.lo
        } else if x.hi < 0.0 {
            -x.hi
        } else {
            0.0
        }
    };
    let largest = m.iter().flatten().map(|x| mag(*x)).fold(0.0, f64::max);
    let mut a: Vec<Vec<I>> = m.into_iter().enumerate().map(|(i, row)| row.into_iter().chain((0..n).map(|j| I::point(if i == j { 1.0 } else { 0.0 }))).collect()).collect();
    for col in 0..n {
        let p = (col..n).max_by(|&x, &y| least(a[x][col]).total_cmp(&least(a[y][col])))?;
        if !(least(a[p][col]) > 1e-9 * largest) {
            return None;
        }
        a.swap(col, p);
        let pivot = a[col][col];
        for j in 0..2 * n {
            a[col][j] = if j == col { I::point(1.0) } else { div_o(a[col][j], pivot) };
        }
        for r in (0..n).filter(|&r| r != col) {
            let k = a[r][col];
            if k == ZERO {
                continue;
            }
            for j in 0..2 * n {
                a[r][j] = if j == col { ZERO } else { sub_o(a[r][j], mul_o(k, a[col][j])) };
            }
        }
    }
    Some(a.into_iter().map(|row| row[n..].to_vec()).collect())
}
