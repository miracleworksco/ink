use super::{Model, Moment, Num, PathStep, Placement};
use ink_core::eval::{Ctx, eval_point};
use ink_core::graph::{Graph, NodeId, Var};
use ink_core::interval::I;
use ink_lang::Doc;

/// What an effort reads of a placement's 12 numbers (Placement), with
/// their derivatives along t and along the coordinates `along` (by index;
/// eps their params) that move them (a partial zero in every entry
/// dropped).
pub(super) fn placement(d: &mut Doc, entries: &[NodeId], along: &[usize], eps: &[u32]) -> Placement {
    let dr: Vec<NodeId> = entries[..9].iter().map(|&e| d.derivative(e, Var::T)).collect();
    let dp: Vec<NodeId> = entries[9..].iter().map(|&e| d.derivative(e, Var::T)).collect();
    let p = nums(&d.g, &entries[9..]);
    let dt2: Vec<NodeId> = dr.iter().chain(&dp).map(|&de| d.derivative(de, Var::T)).collect();
    let dq = partials(d, entries, along, eps);
    let g = &d.g;
    Placement { r: nums(g, &entries[..9]), dr: nums(g, &dr), p, dp: nums(g, &dp), dt2: nums(g, &dt2), dq }
}

/// Numbers of the model, compiled.
pub(super) fn nums<const N: usize>(g: &Graph, ns: &[NodeId]) -> [Num; N] {
    std::array::from_fn(|k| Num::of(g, ns[k]))
}

/// A placement's 12 numbers' derivatives along each coordinate of `along`
/// (by index; eps their params), a partial zero in every entry dropped.
pub(super) fn partials(d: &mut Doc, entries: &[NodeId], along: &[usize], eps: &[u32]) -> Vec<(usize, [Num; 12])> {
    let mut dq = Vec::new();
    for &k in along {
        let dk: Vec<NodeId> = entries.iter().map(|&e| d.derivative(e, Var::Param(eps[k]))).collect();
        if dk.iter().any(|&x| d.g.constant_value(x) != Some(0.0)) {
            dq.push((k, nums(&d.g, &dk)));
        }
    }
    dq
}

/// A frame's 12 numbers: R row major, then T.
pub(super) fn xf_entries(x: &ink_lang::Xf) -> [NodeId; 12] {
    std::array::from_fn(|k| if k < 9 { x.r[k / 3][k % 3] } else { x.t[k - 9] })
}

impl Model {
    /// Whether actuator a moves simulated coordinates alone: in a
    /// simulation its effort is its command, clamped to its limits.
    pub fn commanded(&self, a: usize) -> bool {
        let moved: Vec<usize> = (0..self.coords.len()).filter(|&k| !matches!(self.actuators[a].along[k], Num::C(x) if x == 0.0)).collect();
        !moved.is_empty() && moved.iter().all(|k| self.simulated.iter().any(|s| s.0 == *k))
    }

    /// The mechanism at a state, in a model of a load with states
    /// (Doc::with_states): at moment t, what reads t otherwise than through
    /// the coordinates read there, coordinate k at value q[k][0], speed
    /// q[k][1] and acceleration q[k][2] (one held at a pose stays at it),
    /// every body as measured.
    pub fn state(&self, t: f64, q: &[[f64; 3]]) -> Moment {
        let mut m = self.compute_in(t, self.ctx_with(t, self.state_params(t, q)), &vec![false; self.groups.len()]);
        m.q = q.to_vec();
        m
    }

    /// The mechanism at a state of a simulation (Model::state): a loop's
    /// solve on t found at the state by its own solver from its root at
    /// the state before (Model::roots), not followed along t from t = 0
    /// (which would follow each state's own motion back there: coordinate
    /// partitioning solves the dependent coordinates at each state); its
    /// rates at that root, certified regular there, as at any moment.
    pub fn state_near(&self, t: f64, q: &[[f64; 3]]) -> Moment {
        let mut m = self.compute_in(t, self.near_ctx(t, q), &vec![false; self.groups.len()]);
        m.q = q.to_vec();
        m
    }

    /// The context of a state of a simulation (Model::state_near), its
    /// loops solved there from their roots at the state before.
    pub(super) fn near_ctx(&self, t: f64, q: &[[f64; 3]]) -> Result<Ctx, String> {
        let mut params = self.state_params(t, q);
        match &self.follow {
            None => Ok(Ctx::new(t, params)),
            Some(f) => (|| {
                let mut roots = self.roots.lock().unwrap_or_else(|e| e.into_inner());
                for &i in &self.solves {
                    let (followed, owns, label) = f.followed(i);
                    if !followed {
                        f.provide(i, &mut params, t)?;
                        continue;
                    }
                    let near = roots.get(&i);
                    if let Some(r) = near {
                        for (&k, &v) in owns.iter().zip(r) {
                            params[k] = v;
                        }
                    }
                    // A loop that does not close near its root at the state
                    // before has reached a dead point (or where it cannot
                    // be assembled).
                    f.solve_near(i, &mut params, t).map_err(String::from).map_err(|why| match near {
                        // The solver's own hint (check the equations) points
                        // away from what happened: the loop folds here.
                        Some(_) => format!(
                            "a dead point of {label} at t = {t}: its loop does not close near its root at the state before ({})",
                            why.split("; check the").next().unwrap_or(&why)
                        ),
                        None => why,
                    })?;
                    roots.insert(i, owns.iter().map(|&k| params[k]).collect());
                }
                Ok(Ctx::new(t, params))
            })(),
        }
    }

    /// Over step `s` of the written motion, each coordinate's bound on its
    /// rate (None: its value may jump in the step: a select or floor on t
    /// switches there), bound on its rate's rate (None: its rate may jump),
    /// and its rate's bounds over the step. Err: a coordinate set through a
    /// solve on t, whose bounds over a step are not found.
    pub fn path_step(&self, s: (f64, f64)) -> Result<Vec<PathStep>, String> {
        let ctx = Ctx { t: ink_core::eval::span_mid(s), params: self.params.clone(), t_span: Some(s), p_spans: Default::default() };
        self.coords
            .iter()
            .map(|c| match (&c.q, &c.dq, &c.ddq) {
                (Num::S(..), _, _) | (_, Num::S(..), _) | (_, _, Num::S(..)) => {
                    Err(format!("{} as passed reads what a solve on t sets, so where its path turns or stops is not found", c.name))
                }
                (q, dq, ddq) => Ok(PathStep { speed: self.rate(q, s), accel: self.rate(dq, s), jerk: self.rate(ddq, s), rate: dq.over(&ctx), acc: ddq.over(&ctx) }),
            })
            .collect()
    }

    /// Each coordinate's written rate and its next three derivatives at t,
    /// as enclosures (outward rounded).
    pub fn rates_at(&self, t: f64) -> Result<Vec<[I; 4]>, String> {
        let ctx = Ctx { t, params: self.params.clone(), t_span: Some((t, t)), p_spans: Default::default() };
        self.coords
            .iter()
            .map(|c| match [&c.dq, &c.ddq, &c.dddq, &c.d4q] {
                n if n.iter().any(|n| matches!(n, Num::S(..))) => Err(format!("{} as passed reads what a solve on t sets, so where its path turns or stops is not found", c.name)),
                n => Ok(n.map(|n| n.over(&ctx))),
            })
            .collect()
    }

    /// Each coordinate's value, rate and acceleration as written at t.
    pub fn written(&self, t: f64) -> Result<Vec<[f64; 3]>, String> {
        let ctx = self.ctx(t)?;
        Ok(self.coords.iter().map(|c| [c.q.at(&ctx), c.dq.at(&ctx), c.ddq.at(&ctx)]).collect())
    }

    /// Each coordinate's state at moment t in a simulation (simulate.rs):
    /// a simulated one's `sim` (in Model::simulated's order: position,
    /// speed, acceleration), each other's as written. Err: one written
    /// through a solve on t, which the state's own solves would follow.
    pub(crate) fn beside(&self, t: f64, sim: &[[f64; 3]]) -> Result<Vec<[f64; 3]>, String> {
        let mut params = self.params.clone();
        for (&(_, [x, v, a, j, sn, t0]), s) in self.simulated.iter().zip(sim) {
            (params[x], params[v], params[a], params[j], params[sn], params[t0]) = (s[0], s[1], s[2], 0.0, 0.0, t);
        }
        let ctx = Ctx::new(t, params);
        let mut out = Vec::with_capacity(self.coords.len());
        for (k, c) in self.coords.iter().enumerate() {
            if let Some(i) = self.simulated.iter().position(|s| s.0 == k) {
                out.push(sim[i]);
                continue;
            }
            if [&c.q, &c.dq, &c.ddq].iter().any(|n| matches!(n, Num::S(..))) {
                return Err(format!("{} is written through a solve on t, beside a simulated coordinate: write its motion as an expression of t, or simulate it too", c.name));
            }
            out.push([c.q.at(&ctx), c.dq.at(&ctx), c.ddq.at(&ctx)]);
        }
        Ok(out)
    }

    /// The params of a state (Model::state).
    pub(super) fn state_params(&self, t: f64, q: &[[f64; 3]]) -> Vec<f64> {
        let mut params = self.params.clone();
        // None where no coordinate reads t (a mechanism held at a pose, its
        // springs simulated): nothing to set but the simulated ones.
        if let Some(st) = &self.states {
            for (&(k, jet), s) in st.of.iter().zip(q) {
                if let Some((r, a)) = jet {
                    (params[k], params[r], params[a]) = (s[0], s[1], s[2]);
                }
            }
            params[st.t0] = t;
        }
        // A simulated coordinate's value as a def reads it (a solve, a
        // command) is the state too: its coordinate param (set above) the
        // position, its Taylor terms the speed and acceleration.
        for &(k, [x, v, a, j, s, t0]) in &self.simulated {
            (params[x], params[v], params[a], params[j], params[s], params[t0]) = (0.0, q[k][1], q[k][2], 0.0, 0.0, t);
        }
        params
    }

    /// Each actuator's limits (lo, hi), and each joint's speed limits with
    /// its rate (rate, lo, hi), at a state (Model::state): what a torque-
    /// speed line or a published curve gives at that speed.
    pub fn limits_at(&self, t: f64, q: &[[f64; 3]]) -> Result<(Vec<(f64, f64)>, Vec<(f64, f64, f64)>), String> {
        let ctx = self.ctx_with(t, self.state_params(t, q))?;
        let efforts = self.actuators.iter().map(|a| (a.lo.at(&ctx), a.hi.at(&ctx))).collect();
        let speeds = self.speeds.iter().map(|s| (s.rate.at(&ctx), s.lo.at(&ctx), s.hi.at(&ctx))).collect();
        Ok((efforts, speeds))
    }

    /// Each actuator's limits (lo, hi) bounded over a box of states at
    /// moment t: each coordinate's position, rate and acceleration anywhere
    /// in its interval (Model::state's, as intervals), each limit's least
    /// and most there by the interval evaluation an effort's bound over a
    /// step reads them with (bounded.rs), rounded outward. So a limit that
    /// follows the drive's speed (a torque-speed line, a pull-out curve, a
    /// band of speeds it gives nothing in) is bounded over every speed the
    /// box holds, not read at one. Unbounded through a solve on t.
    pub fn limits_over(&self, t: f64, q: &[[I; 3]]) -> Result<Vec<(I, I)>, String> {
        let mid: Vec<[f64; 3]> = q.iter().map(|c| c.map(|x| x.mid())).collect();
        let mut ctx = self.ctx_with(t, self.state_params(t, &mid))?;
        let mut spans: Vec<(u32, I)> = Vec::new();
        if let Some(st) = &self.states {
            for (&(k, jet), s) in st.of.iter().zip(q) {
                if let Some((r, a)) = jet {
                    spans.extend([(k, s[0]), (r, s[1]), (a, s[2])].map(|(p, x)| (p as u32, x)));
                }
            }
        }
        for &(k, [_, v, a, ..]) in &self.simulated {
            spans.extend([(v as u32, q[k][1]), (a as u32, q[k][2])]);
        }
        if let Some((_, x)) = spans.iter().find(|(_, x)| !(x.lo <= x.hi && x.lo.is_finite() && x.hi.is_finite())) {
            return Err(format!("a state's bound [{}, {}] is not a finite interval", x.lo, x.hi));
        }
        ctx.p_spans = ink_core::eval::ParamSpans::new(spans.into_iter().filter(|(_, x)| x.lo < x.hi).map(|(p, x)| (p, (x.lo, x.hi))));
        Ok(self.actuators.iter().map(|a| (a.lo.over(&ctx), a.hi.over(&ctx))).collect())
    }

    /// Whether each actuator's limits read the state (a coordinate's
    /// position, rate or acceleration as Model::state passes it): where
    /// not, limits_over is limits_at's at any state.
    pub fn limits_read_state(&self) -> Vec<bool> {
        let mut state: Vec<u32> = Vec::new();
        if let Some(st) = &self.states {
            for &(k, jet) in &st.of {
                if let Some((r, a)) = jet {
                    state.extend([k, r, a].map(|p| p as u32));
                }
            }
        }
        for &(_, [_, v, a, ..]) in &self.simulated {
            state.extend([v, a].map(|p| p as u32));
        }
        let reads = |n: &Num| match n {
            Num::C(_) => false,
            Num::T(t) => t.params_read().iter().any(|p| state.contains(p)),
            Num::S(..) => true,
        };
        self.actuators.iter().map(|a| reads(&a.lo) || reads(&a.hi)).collect()
    }

    /// Its context at moment t: the solves on t it reads solved there.
    /// Err: one that does not solve there, or whose rates are not defined
    /// there, why.
    pub(super) fn ctx(&self, t: f64) -> Result<Ctx, String> {
        self.ctx_with(t, self.params.clone())
    }

    /// Its context at moment t with these params (a state's).
    pub(super) fn ctx_with(&self, t: f64, mut params: Vec<f64>) -> Result<Ctx, String> {
        if let Some(f) = &self.follow {
            for &i in &self.solves {
                f.provide(i, &mut params, t)?;
            }
        }
        Ok(Ctx::new(t, params))
    }

    /// The largest |d/dt| of n over step `s` (None where it has no bound:
    /// a jump). Through a solve, the largest of its rate at the step's
    /// ends and middle (exact there), as no bound over the step is known;
    /// None where it is not defined at one of them.
    pub(super) fn rate(&self, n: &Num, s: (f64, f64)) -> Option<f64> {
        match n {
            Num::C(_) => Some(0.0),
            Num::T(t) => crate::steps::rate(t, &Ctx::new(0.0, self.params.clone()), [I::point(0.0); 3], s),
            Num::S(_, r) => [s.0, ink_core::eval::span_mid(s), s.1].into_iter().try_fold(0.0f64, |m, t| Some(m.max(eval_point(r, &self.ctx(t).ok()?, [0.0; 3]).abs()))),
        }
    }

    /// A rate bound already proves continuity. When it is unavailable,
    /// compare the exact values where the coordinate's formula switches.
    pub(super) fn may_jump(&self, n: &Num, s: (f64, f64)) -> bool {
        self.rate(n, s).is_none()
            && match n {
                Num::T(tape) => !ink_core::eval::continuous_in_time(tape, &Ctx::new(0.0, self.params.clone()), [I::point(0.0); 3], s),
                _ => true,
            }
    }

    /// Whether, over the moments `span` of a simulation from state `sim` at
    /// its start, what drives it may jump in t: an actuator's command or
    /// limits, a load, or a written coordinate's position, speed or
    /// acceleration (a step in a commanded angle), whose rate over the span
    /// has no bound (a select or floor on t switching there) and whose
    /// bounds over it are wider than 1e-6 (and its rounding, 1e-12 of its
    /// size): a jump keeps its height however short the span, whatever its
    /// origin (a step of 1.8 deg after thousands of turns); a slope with no
    /// bound, a square root at 0, does not.
    pub(crate) fn switches(&self, span: (f64, f64), sim: &[[f64; 3]]) -> bool {
        let Ok(q) = self.beside(span.0, sim) else { return false };
        let ctx = Ctx::new(span.0, self.state_params(span.0, &q));
        let over = Ctx { t_span: Some(span), ..ctx.clone() };
        let jumps = |n: &Num| match n {
            Num::T(tape) => {
                crate::steps::rate(tape, &ctx, [I::point(0.0); 3], span).is_none() && {
                    let b = ink_core::eval::eval_interval(tape, &over, [I::point(0.0); 3], None);
                    !(b.hi - b.lo <= 1e-6 + 1e-12 * n.at(&ctx).abs())
                }
            }
            _ => false,
        };
        self.actuators.iter().any(|a| a.command.iter().chain([&a.lo, &a.hi]).any(jumps))
            || self.loads.iter().any(|l| l.force.iter().chain(&l.torque).any(jumps))
            || self.coords.iter().enumerate().any(|(k, c)| c.moves && !self.simulated.iter().any(|s| s.0 == k) && [&c.q, &c.dq, &c.ddq].into_iter().any(jumps))
    }
}
