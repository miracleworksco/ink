use super::{Coordinate, Model, Moment, Num, ROUNDING, V3};
use ink_core::eval::Ctx;
use ink_core::graph::{NodeId, Var};
use ink_core::interval::I;
use ink_core::tape::Tape;
use ink_core::v3::norm;
use ink_lang::{Doc, Steps};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// A quantity of a placement an effort rests on, whose steps keep it
/// within 1% inside each: its components, their rates in t (whose rates
/// over a step bound its second derivative), and its scale where it has
/// its own (a rotation's unit entries; else its peak at the moments
/// measured).
pub(super) struct Tracked {
    /// f″, f‴, ... for a partial velocity whose components cancel.
    /// A centred enclosure can resolve its own peak on a step where
    /// bounding f″ directly loses the cancellation.
    curvature: Vec<Tape>,
    values: Vec<Num>,
    rates: Vec<Num>,
    unit: Option<f64>,
    /// Tracked only where this group's spin is resolved (Model::choose).
    group: Option<usize>,
    /// A body only reactions carry: tracked only in the steps of a
    /// claimed one of these (by index); none, in every choice.
    serves: Vec<usize>,
}

impl Tracked {
    /// Taylor's bound for f″ about the centre of this interval. For each
    /// order, enclose its coefficient at the centre and its next derivative
    /// everywhere on the interval. A finite rate at every order used is
    /// essential: differentiating either side of a jump can give zero.
    fn curvature_fits(&self, step: (f64, f64), params: &[f64], fits: impl Fn(f64) -> bool) -> bool {
        use ink_core::interval::{add_up, mag, mul_up};
        if self.curvature.is_empty() {
            return false;
        }
        let t = ink_core::eval::span_mid(step);
        let centre = Ctx { t_span: Some((t, t)), ..Ctx::new(t, params.to_vec()) };
        let radius = (t - step.0).max(step.1 - t).next_up();
        let origin = [I::point(0.0); 3];
        let (mut sum, mut factor) = (0.0, 1.0); // factor = radius^k / k!
        for (k, derivative) in self.curvature.iter().enumerate() {
            let coefficient = mag(ink_core::eval::eval_interval(derivative, &centre, origin, None));
            sum = add_up(sum, mul_up(factor, coefficient));
            if !fits(sum) {
                return false;
            }
            let Some(next) = crate::steps::rate(derivative, &centre, origin, step) else { return false };
            factor = mul_up(factor, (radius / (k + 1) as f64).next_up());
            if fits(add_up(sum, mul_up(factor, next))) {
                return true;
            }
        }
        false
    }
}

/// A claim on a reaction its own steps serve (Model::steps_for): the
/// reaction's index, the claim's (its provider's) key, and the frame the
/// claim reads it in at each moment (None: its path's).
pub struct Serve<'a> {
    pub reaction: usize,
    pub claim: usize,
    pub frame: Option<&'a (dyn Fn(f64) -> [f64; 12] + Sync)>,
}

impl Serve<'_> {
    fn wrench(&self, m: &Moment) -> (V3, V3) {
        let w = &m.reactions[self.reaction];
        match self.frame {
            Some(frame) => w.in_frame(&frame(m.t)),
            None => w.local(),
        }
    }
}

/// Steps chosen over a span, the groups read with their spin at rest
/// there, and over each step what those groups' departure from round can
/// move the effort along each coordinate, summed (Group::moves).
pub(super) struct Chosen {
    steps: Steps,
    rest: Vec<bool>,
    bounds: Vec<Vec<f64>>,
}

/// How steps over a span are chosen (Model::steps); what is read with its
/// spin at rest for its departure from round follows it.
pub const STEPS: &str = "measured at the span's two ends and at the middle of each step, each moment an actuator's limit may jump (a select \
    or floor on t) isolated by bisection to the resolution of t in a step of its own, so jumps apart in t are found apart; a step a limit may \
    jump in whose range over it, by intervals, can cross the effort is loose; steps split until, inside each, no coordinate's speed or \
    acceleration, no moving body's rotation, translation's acceleration or partial velocity, no carried body's shift's speed or acceleration, and no \
    load's force or torque, can stray from its value at the middle by more than 1% (of its peak; of a turn, of a radian's sine): max|f''| h^2 / 8 from \
    bounds of their rates over the step, the bodies spinning with one link read with the spin at rest where they are round about its axis, or where \
    their departure from round (their mass off the axis, their inertia not round about it), bounded over each step with the spin's and its frame's \
    rates, can move no effort by more than 1% of its peak, or, their spin past the steps' budget, no actuator's effort past its claim's margin: read \
    round, each effort then within that bound of the real one at any angle; then until each effort lies within 4% of its peak of the line through \
    its neighbours; and a step of its own where an actuator's effort, or its margin to a limit, turns back near its extreme between the moments \
    measured (by bisection on the sign of its rate); a directly ranged rating that can still cross the measured effort inside its step leaves it loose: an estimate between steps";

/// What a claim on a reaction's own steps add to the efforts' (STEPS).
pub const REACTION_STEPS: &str = "for this claim on a reaction, besides: what only the reaction carries (a part t moves) tracked so too, and each of its \
    components in the claim's frame within 4% of its peak of the line through its neighbours, and what is not known of a spinning body it \
    carries moving it by no more than that over a step, with a step at each located component extreme";

/// The share of its peak by which a peak between steps may be missed.
const SHARE: f64 = 0.01;

/// The share of its peak by which an effort at a moment may lie off the
/// line through its neighbours' before its steps are halved (STEPS): what
/// the steps may read a peak between their moments low by. A rule the
/// steps are chosen by, not a bound.
pub const OFF_LINE: f64 = 4.0 * SHARE;

/// What an effort rests on that moves, each once (Tracked), as compiled.
#[derive(Default)]
pub(super) struct Tracker {
    seen: HashMap<(Vec<NodeId>, Option<usize>, Vec<usize>), usize>,
    out: Vec<Tracked>,
    /// Each partial velocity as a vector, including constant components for
    /// its peak. An index names a tracked component eligible for Taylor bounds.
    partials: Vec<Vec<(NodeId, Option<usize>)>>,
}

impl Tracker {
    /// values (against `unit`, else their peak), where they read t, their
    /// rates their t-derivatives; tracked only where `group`'s spin is
    /// resolved, if it names one, and only for the reactions it `serves`,
    /// if any.
    pub(super) fn add(&mut self, d: &mut Doc, values: &[NodeId], unit: Option<f64>, group: Option<usize>, serves: &[usize]) {
        if unit.is_none() && values.len() > 1 {
            for n in values {
                self.add(d, std::slice::from_ref(n), None, group, serves);
            }
            return;
        }
        self.track(d, values, unit, group, serves);
    }

    fn track(&mut self, d: &mut Doc, values: &[NodeId], unit: Option<f64>, group: Option<usize>, serves: &[usize]) -> Option<usize> {
        if !values.iter().any(|&n| changes_in_t(d, n)) {
            return None;
        }
        let index = self.seen.entry((values.to_vec(), group, serves.to_vec())).or_insert_with(|| {
            let rates: Vec<NodeId> = values.iter().map(|&n| d.derivative(n, Var::T)).collect();
            let values = values.iter().map(|&n| Num::of(&d.g, n)).collect();
            let rates = rates.iter().map(|&n| Num::timed(d, n)).collect();
            self.out.push(Tracked { curvature: Vec::new(), values, rates, unit, group, serves: serves.to_vec() });
            self.out.len() - 1
        });
        Some(*index)
    }

    pub(super) fn partial(&mut self, d: &mut Doc, values: &[NodeId], group: Option<usize>, serves: &[usize]) {
        let parts: Vec<_> = values
            .iter()
            .map(|&n| {
                let index = self.track(d, &[n], None, group, serves).filter(|_| !d.solves_read(&[n]).iter().any(|&i| d.solves[i].on_t));
                (n, index)
            })
            .collect();
        if parts.iter().any(|(_, i)| i.is_some()) {
            self.partials.push(parts);
        }
    }

    pub(super) fn finish(mut self, d: &mut Doc, model: &Model, t: f64) -> Vec<Tracked> {
        // Placement derivatives through a solve can append rate parameters.
        // Sample only after the model has all of them and can solve them at t.
        debug_assert_eq!(model.params.len(), d.params.len(), "tracking requires the compiled document's full parameter table");
        if self.partials.is_empty() {
            return self.out;
        }
        let Ok(ctx) = model.ctx(t) else { return self.out }; // no sample at a dead point
        for parts in self.partials {
            // The vector peak only selects where to try a tighter enclosure.
            // Every component still has to meet its own 1% test.
            let samples: Vec<_> = parts.iter().map(|&(n, _)| d.scalar(n, &ctx).abs()).collect();
            let peak = samples.iter().copied().fold(0.0, f64::max);
            for ((n, index), sample) in parts.into_iter().zip(samples) {
                let Some(i) = index else { continue };
                let q = &mut self.out[i];
                if peak == 0.0 || sample > ROUNDING * peak || !q.curvature.is_empty() {
                    continue;
                }
                // Limit the cost of the optional enclosure. If it cannot
                // resolve the step, normal subdivision and refusal apply.
                let mut derivative = d.g.derivative(n, Var::T);
                for _ in 0..5 {
                    if d.g.reachable(&[derivative]).len() > 8192 {
                        break;
                    }
                    derivative = d.g.derivative(derivative, Var::T);
                    q.curvature.push(Tape::new(&d.g, derivative));
                }
            }
        }
        self.out
    }
}

/// Whether `n` changes with t: it reads t, or what a solve on t sets.
fn changes_in_t(d: &Doc, n: NodeId) -> bool {
    d.g.depends_on(n, |x| x == Var::T) || d.solves_read(&[n]).iter().any(|&i| d.solves[i].on_t)
}

/// The larger of two bounds, NaN (no bound) where either is.
fn worst(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() { f64::NAN } else { a.max(b) }
}

#[cfg(test)]
thread_local! {
    /// How many times this thread has chosen a motion's steps (Model::choose).
    pub(crate) static CHOSEN: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

impl Model {
    /// The model at moment t as the steps over `span` read it (Model::
    /// steps): the groups read there with their spin at rest read round,
    /// with how far each effort may be from the real one at any angle of
    /// their spins (`bound`, `bound_f`).
    pub fn over(&self, span: (f64, f64), t: f64) -> Arc<Moment> {
        let c = self.chosen(span, None);
        let m = self.moment(t, &c.rest);
        let ends = &c.steps.ends;
        let j = ends.partition_point(|&e| e <= t).clamp(1, ends.len() - 1) - 1;
        if c.bounds[j].iter().all(|b| *b == 0.0) {
            return m;
        }
        let mut m = (*m).clone();
        m.bound_f = c.bounds[j].clone();
        m.bound = m.pm.iter().map(|p| p.iter().zip(&m.bound_f).map(|(p, d)| p.abs() * d).sum()).collect();
        m.bound_u = m.uncovered.iter().map(|q| q.iter().zip(&m.bound_f).map(|(q, d)| q * d).sum()).collect();
        Arc::new(m)
    }

    /// Actuator a's RMS effort over `span`, sqrt of the mean of its square
    /// over the steps (Model::steps), each step's middle standing for it
    /// (the midpoint rule), the moments as the steps read them
    /// (Model::over); with the most the groups read round there, and what
    /// the bodies' measures are known to (Moment::mass), can move
    /// it (an RMS moves by no more than its samples do). A step left loose
    /// (coarser than the steps' rule, at the budget) is read at its middle
    /// as the rest. Err: a moment it is refused at, or a step left unchecked
    /// (a jump in a coordinate's position or speed: an impulse), why.
    pub fn rms(&self, a: usize, span: (f64, f64)) -> Result<(f64, f64), String> {
        let steps = self.steps(span);
        if let Some(&(x, y)) = steps.unchecked.first() {
            return Err(format!("its effort over [{x}, {y}] is not measured (a coordinate's position or speed may jump there: an impulse), so its RMS over the span is not found"));
        }
        let (mut sum, mut time, mut band) = (0.0, 0.0, 0.0f64);
        for w in steps.ends.windows(2) {
            let m = self.over(span, ink_core::eval::span_mid((w[0], w[1])));
            if let Some(why) = self.refusal(a, &m) {
                return Err(why);
            }
            sum += (w[1] - w[0]) * m.effort[a] * m.effort[a];
            time += w[1] - w[0];
            band = band.max(m.bound[a] + m.mass[a]);
        }
        Ok(((sum / time).sqrt(), band))
    }

    /// The model at moment t as its effort claims read it: inside the
    /// file's span as the steps over it read it (Model::over), elsewhere at
    /// one moment (Model::at). A function of t alone, so a moment's value
    /// is the same whichever claim asks first.
    pub fn read(&self, t: f64) -> Arc<Moment> {
        match self.span {
            Some(s) if s.0 <= t && t <= s.1 => self.over(s, t),
            _ => self.at(t),
        }
    }

    /// The steps over `span` an effort is measured at (STEPS), with the
    /// span's two ends: from 16, halved where a coordinate's speed or
    /// acceleration, or a tracked quantity of a moving body (Tracked: a
    /// round body's with its spin at rest, Spin), could stray inside the
    /// step from its value at the middle by more than 1% (of its peak at
    /// the steps' middles; of a turn's unit entries): within max|f''| h^2 /
    /// 8, f'' bounded over the step (a coordinate's jerk and the jerk's
    /// rate, a rotation's R''); then where an effort
    /// (each actuator's, else each coordinate's) is off the line through
    /// its neighbours by more than 4% of its peak and its rounding. At most
    /// MAX_STEPS: those left wanting are `loose`; those where a coordinate's
    /// position or speed may jump (their rates have no bound however short),
    /// `unchecked`.
    pub fn steps(&self, span: (f64, f64)) -> Steps {
        self.chosen(span, None).steps.clone()
    }

    /// The steps over `span` a claim on a reaction is measured at: the
    /// efforts' rules, and besides what only the reaction carries, each of
    /// its components in the claim's frame in the line test, and the groups
    /// it carries read at rest where their bound on it is small or their
    /// spin past the budget (Model::shaken). The claim's own choice: no
    /// other claim's steps follow it.
    pub fn steps_for(&self, span: (f64, f64), serve: &Serve) -> Steps {
        self.chosen(span, Some(serve)).steps.clone()
    }

    /// Reaction r at moment t as the steps over `span` for `serve` read it
    /// (the groups read at rest there read round), and, for each such group
    /// r carries, the most its departure from round can move r's force (N)
    /// and its moment about the group's axis point (N mm) at any angle of
    /// its spin over the step holding t, with that point at t.
    pub fn shaken(&self, span: (f64, f64), serve: Option<&Serve>, r: usize, t: f64) -> (Arc<Moment>, Vec<(f64, f64, V3)>) {
        let c = self.chosen(span, serve);
        let m = self.moment(t, &c.rest);
        let ends = &c.steps.ends;
        let j = ends.partition_point(|&e| e <= t).clamp(1, ends.len() - 1) - 1;
        let step = (ends[j], ends[j + 1]);
        let ctx = Ctx::new(t, self.params.clone());
        let over = Ctx { t: ink_core::eval::span_mid(step), t_span: Some(step), ..ctx.clone() };
        let shakes = (0..self.groups.len())
            .filter(|&g| c.rest[g] && self.carries(r, g))
            .map(|g| {
                let gr = &self.groups[g];
                let (f, mo) = gr.shakes(self.departure(gr), &over, self.gravity);
                (f, mo, std::array::from_fn(|k| gr.frame[9 + k].at(&ctx)))
            })
            // And what is not known of every group it carries, read at rest or not.
            .chain(self.unknown_shakes(r, &over, &ctx))
            .collect();
        (m, shakes)
    }

    fn chosen(&self, span: (f64, f64), serve: Option<&Serve>) -> Arc<Chosen> {
        let key = (span.0.to_bits(), span.1.to_bits(), serve.map(|s| s.claim));
        if let Some(c) = self.chosen.lock().unwrap().get(&key) {
            return c.clone();
        }
        let c = Arc::new(self.choose(span, serve));
        self.chosen.lock().unwrap().insert(key, c.clone());
        c
    }

    fn choose(&self, (a, b): (f64, f64), serving: Option<&Serve>) -> Chosen {
        #[cfg(test)]
        CHOSEN.set(CHOSEN.get() + 1);
        let serve = serving.map(|s| s.reaction);
        use crate::steps::{MAX_STEPS, START_STEPS};
        use rayon::prelude::*;
        let mid = |s: (f64, f64)| ink_core::eval::span_mid(s);
        // From even steps, and a step of its own at each moment an
        // actuator's limit may jump (a drive disabled for a moment): each
        // side of it is measured, however short.
        let mut ends: Vec<f64> = (0..=START_STEPS).map(|k| if k == START_STEPS { b } else { a + (b - a) * (k as f64 / START_STEPS as f64) }).collect();
        let limits: Vec<&Num> = self.actuators.iter().flat_map(|x| [&x.lo, &x.hi]).collect();
        ends.extend(self.jumps(&limits, (a, b)).into_iter().flat_map(|s| [s.0, s.1]));
        ends.sort_by(f64::total_cmp);
        ends.dedup();
        let mut steps: Vec<(f64, f64)> = ends.windows(2).map(|w| (w[0], w[1])).collect();
        let shortest = (b - a) / MAX_STEPS as f64;
        let (mut unchecked, mut loose) = (Vec::new(), Vec::new());
        // The motion: each coordinate's speed and acceleration, each turning
        // body's rotation. A step where a coordinate's position or speed
        // has no rate bound is split to the finest and left unchecked if
        // it may still jump there: an effort needs an impulse. Where only
        // acceleration or jerk jumps (a trapezoidal profile's corner),
        // the step is split to the finest and measured either side.
        let miss = |f2: f64, s: (f64, f64)| f2 * (s.1 - s.0) * (s.1 - s.0) / 8.0;
        let finest = |s: (f64, f64)| s.1 - s.0 <= shortest;
        // Each group's spin is read at rest until the steps are chosen
        // without it; then those whose departure from round can move an
        // effort by more than the steps keep are resolved too.
        let (mut resolved, mut decided, mut at_rest) = (vec![false; self.groups.len()], self.groups.is_empty(), Vec::new());
        let mut rest = vec![false; self.groups.len()];
        loop {
            let tracked: Vec<&Tracked> =
                self.tracked.iter().filter(|q| q.group.is_none_or(|g| resolved[g]) && (q.serves.is_empty() || serve.is_some_and(|r| q.serves.contains(&r)))).collect();
            // The moments the steps are read at (where the solves solve).
            let middles: Vec<Option<Ctx>> = steps.par_iter().map(|&s| self.ctx(mid(s)).ok()).collect();
            let at = |c: &Coordinate, x: &Ctx| (c.dq.at(x).abs(), c.ddq.at(x).abs());
            let peaks: Vec<(f64, f64)> =
                self.coords.iter().map(|c| middles.iter().flatten().map(|x| at(c, x)).fold((0.0f64, 0.0f64), |p, v| (p.0.max(v.0), p.1.max(v.1)))).collect();
            // Each tracked quantity's scale: its own, else its peak at the
            // span's ends and the steps' middles.
            let moments: Vec<Ctx> = [a, b].into_iter().filter_map(|t| self.ctx(t).ok()).chain(middles.iter().flatten().cloned()).collect();
            let scales: Vec<f64> =
                tracked.par_iter().map(|q| q.unit.unwrap_or_else(|| moments.iter().flat_map(|m| q.values.iter().map(move |n| n.at(m).abs())).fold(0.0, f64::max))).collect();
            // (a jump; a peak inside may be missed by more than 1%; a rate
            // with no bound, measured either side at the finest)
            let state = |(s, middle): ((f64, f64), &Option<Ctx>)| -> (bool, bool, bool) {
                // A middle where the solves do not solve is not read: halved,
                // to find where they do. Not where its ends do not solve either
                // (past where a simulation stopped, no moment does): it is
                // unchecked whole (below), not halved to the finest.
                if middle.is_none() {
                    return (false, self.ctx(s.0).is_ok() || self.ctx(s.1).is_ok(), false);
                }
                let (mut jump, mut want, mut unbounded) = (false, false, false);
                // Whether the step is too long for a second derivative
                // bounded by `bound` against `scale`; None: not bounded.
                let short = |bound: Option<f64>, scale: f64| bound.map(|b| miss(b, s) > SHARE * scale);
                let mut over = |r: Option<bool>| match r {
                    None => unbounded = true,
                    Some(w) => want |= w,
                };
                for (c, &(v, acc)) in self.coords.iter().zip(&peaks) {
                    jump |= self.rate(&c.q, s).is_none() || self.rate(&c.dq, s).is_none();
                    over(short(self.rate(&c.ddq, s), v));
                    over(short(self.rate(&c.dddq, s), acc));
                }
                for (q, &scale) in tracked.iter().zip(&scales) {
                    let bound = q.rates.iter().try_fold(0.0f64, |m, n| Some(m.max(self.rate(n, s)?)));
                    match short(bound, scale) {
                        None => unbounded = true,
                        Some(true) => {
                            // Bound the same step's curvature more tightly before
                            // splitting it. At the last step its two halves can
                            // bound the curvature separately, still over the full width.
                            use ink_core::interval::mul_up;
                            let width = (s.1 - s.0).next_up();
                            let factor = mul_up(mul_up(width, 0.125), width);
                            let fits = |piece| q.curvature_fits(piece, &self.params, |bound| mul_up(bound, factor) <= SHARE * scale);
                            let half = mid(s);
                            if !(fits(s) || finest(s) && fits((s.0, half)) && fits((half, s.1))) {
                                return (jump, true, unbounded);
                            }
                        }
                        Some(false) => {}
                    }
                }
                (jump, want, unbounded)
            };
            let states: Vec<(bool, bool, bool)> = steps.par_iter().copied().zip(&middles).map(state).collect();
            let wanted: HashSet<(u64, u64)> = steps.iter().zip(&states).filter(|(_, (j, w, u))| *j || *w || *u).map(|(s, _)| (s.0.to_bits(), s.1.to_bits())).collect();
            let (next, split) = self.split(&steps, shortest, |s| wanted.contains(&(s.0.to_bits(), s.1.to_bits())));
            if !split {
                if !decided {
                    decided = true;
                    (resolved, at_rest) = self.round_enough(&steps, (a, b), serve);
                    rest = resolved.iter().map(|r| !r).collect();
                    if resolved.contains(&true) {
                        continue;
                    }
                }
                for (&s, &(jump, want, _)) in steps.iter().zip(&states) {
                    if jump && self.coords.iter().any(|c| self.may_jump(&c.q, s) || self.may_jump(&c.dq, s)) {
                        unchecked.push(s);
                    } else if s.0.next_up() < s.1 && (want || jump && !finest(s)) {
                        loose.push(s);
                    }
                }
                break;
            }
            steps = next;
        }
        // The efforts, at the span's ends and the steps' middles: where one
        // bends between its neighbours, the steps there are halved.
        loop {
            let live: Vec<usize> = (0..steps.len()).filter(|&j| !unchecked.contains(&steps[j])).collect();
            // Each moment measured, with the step it stands for.
            let moments: Vec<(f64, usize)> = std::iter::once((a, 0)).chain(live.iter().map(|&j| (mid(steps[j]), j))).chain([(b, steps.len() - 1)]).collect();
            let at: Vec<Arc<Moment>> = moments.par_iter().map(|&(t, _)| self.moment(t, &rest)).collect();
            // Those where the solves do not solve are not measured.
            let (moments, at): (Vec<(f64, usize)>, Vec<Arc<Moment>>) = moments.into_iter().zip(at).filter(|(_, m)| m.unsolved.is_none()).unzip();
            let found: Vec<usize> = (0..self.actuators.len()).filter(|&k| !self.actuators[k].coupled.is_empty()).collect();
            let mut series: Vec<(Vec<f64>, Vec<f64>)> = match found.is_empty() {
                false => found.iter().map(|&k| at.iter().map(|m| (m.effort[k], m.noise[k])).unzip()).collect(),
                true => (0..self.coords.len()).map(|k| at.iter().map(|m| (m.f[k], m.noise_f[k])).unzip()).collect(),
            };
            // The reaction they serve: its force and moment in the claim's
            // frame.
            if let Some(s) = serving {
                let read: Vec<(V3, V3, (f64, f64))> = at
                    .iter()
                    .map(|m| {
                        let w = &m.reactions[s.reaction];
                        let (f, mo) = s.wrench(m);
                        (f, mo, w.noise)
                    })
                    .collect();
                for c in 0..6 {
                    series.push(read.iter().map(|x| if c < 3 { (x.0[c], x.2.0) } else { (x.1[c - 3], x.2.1) }).unzip());
                }
            }
            let mut bent: HashSet<usize> = HashSet::new();
            for (v, noise) in &series {
                let peak = v.iter().fold(0.0f64, |p, x| p.max(x.abs()));
                for i in 1..v.len().saturating_sub(1) {
                    let (t0, t1, t2) = (moments[i - 1].0, moments[i].0, moments[i + 1].0);
                    let line = v[i - 1] + (v[i + 1] - v[i - 1]) * (t1 - t0) / (t2 - t0);
                    let floor = (OFF_LINE * peak).max(noise[i - 1]).max(noise[i]).max(noise[i + 1]);
                    if !((v[i] - line).abs() <= floor) {
                        bent.extend([moments[i - 1].1, moments[i].1, moments[i + 1].1]);
                    }
                }
            }
            for (&(t, j), m) in moments.iter().zip(&at) {
                let s = steps[j];
                if t == mid(s) && found.iter().any(|&k| m.margin(k) >= 0.0 && self.rating_dip(k, s, m.effort[k])) {
                    bent.insert(j);
                }
            }
            // What is not known of the spinning bodies the reaction carries moves it, over a step, by
            // what that unknown gives at the spin the step bounds: a step where that is more than the
            // line test's share of the reaction's peak force is halved, as a bent one is. Read round
            // or not, a smaller unknown then leaves no wider band.
            if let Some(serve) = serving {
                let peak = at.iter().map(|m| ink_core::v3::norm(serve.wrench(m).0)).fold(0.0f64, f64::max);
                let floor = at.iter().map(|m| m.reactions[serve.reaction].noise.0).fold(OFF_LINE * peak, f64::max);
                let wide: Vec<usize> = live
                    .par_iter()
                    .copied()
                    .filter(|&j| {
                        let s = steps[j];
                        let ctx = Ctx::new(mid(s), self.params.clone());
                        let over = Ctx { t_span: Some(s), ..ctx.clone() };
                        let band: f64 = self.unknown_shakes(serve.reaction, &over, &ctx).iter().map(|x| x.0).sum();
                        !(band <= floor)
                    })
                    .collect();
                bent.extend(wide);
            }
            let bent: HashSet<(u64, u64)> = bent.iter().map(|&j| (steps[j].0.to_bits(), steps[j].1.to_bits())).collect();
            let (next, split) = self.split(&steps, shortest, |s| bent.contains(&(s.0.to_bits(), s.1.to_bits())));
            if !split {
                let more: Vec<(f64, f64)> = steps.iter().copied().filter(|s| bent.contains(&(s.0.to_bits(), s.1.to_bits())) && !loose.contains(s) && !finest(*s)).collect();
                loose.extend(more);
                break;
            }
            steps = next;
        }
        // A step an actuator's limit may jump in across its effort is loose:
        // its middle decides nothing there (Model::limit_undecided).
        for &s in &steps {
            // (A step with no float inside it is its two ends, each read.)
            if loose.contains(&s) || unchecked.contains(&s) || s.0.next_up() >= s.1 {
                continue;
            }
            let m = self.moment(mid(s), &rest);
            if m.unsolved.is_none() && (0..self.actuators.len()).any(|k| self.limit_undecided(k, s, m.effort[k])) {
                loose.push(s);
            }
        }
        // Where an actuator's effort, or its margin to either limit, turns
        // back between the moments measured (a peak the steps read up to
        // their 4% low), found as a simulated coordinate's extremes are
        // (steps::turn_back), its rate by differences: a step of its own
        // there, so its claims and peaks read it.
        steps = self.turns(steps, (a, b), &rest, &mut loose, &unchecked, serving);
        // A step whose middle, or an end of the span, where the solves do
        // not solve (or their rates are not defined: a dead point) is
        // unchecked: nothing is found there.
        for s in steps.iter().map(|&s| (s, mid(s))).chain([((a, a), a), ((b, b), b)]) {
            if !unchecked.contains(&s.0) && self.at(s.1).unsolved.is_some() {
                loose.retain(|&l| l != s.0);
                unchecked.push(s.0);
            }
        }
        // An end no moment can be read at (a time under ink's least, 1e-100, but 0: where halving toward a jump at
        // 0 ends) is no end: its steps are one, unchecked or loose if any part was.
        let unread = |t: f64| t != 0.0 && t.abs() < 1e-100;
        let mut merged: Vec<(f64, f64)> = Vec::with_capacity(steps.len());
        for &s in &steps {
            match merged.last_mut() {
                Some(last) if unread(last.1) => {
                    let whole = (last.0, s.1);
                    for list in [&mut unchecked, &mut loose] {
                        if list.contains(last) || list.contains(&s) {
                            list.retain(|x| *x != *last && *x != s);
                            list.push(whole);
                        }
                    }
                    *last = whole;
                }
                _ => merged.push(s),
            }
        }
        let steps = merged;
        let mut ends = vec![a];
        ends.extend(steps.iter().map(|s| s.1));
        let mut basis = match at_rest.is_empty() {
            true => STEPS.to_string(),
            false => format!("{STEPS}; {}", at_rest.join("; ")),
        };
        if serving.is_some() {
            basis = format!("{basis}; {REACTION_STEPS}");
        }
        if !self.made_round.is_empty() {
            let each: Vec<String> = self.made_round.iter().map(|m| format!("{} ({:.3e} g mm, {:.3e} g mm2)", m.body, m.first, m.second)).collect();
            basis = format!(
                "{basis}; spinning bodies whose departure from round lies within what their measure is known to (twice what its last refinement moved \
                 its first moments, its second): read round there, the measure's error being no imbalance, and an imbalance under that is unknown, so a \
                 reaction that carries one is known only to what such an imbalance would add at its spin: {}",
                each.join(", ")
            );
        }
        // Over each step, what the groups read at rest can move each
        // effort, summed.
        let nq = self.coords.len();
        let bounds: Vec<Vec<f64>> = match rest.contains(&true) {
            false => vec![vec![0.0; nq]; steps.len()],
            true => self.bounds(&steps, &rest).into_iter().map(|per| per.iter().fold(vec![0.0; nq], |s, d| s.iter().zip(d).map(|(x, y)| x + y).collect())).collect(),
        };
        // And what is not known of every group's departure, at rest or not, over each step.
        let ctx = Ctx::new(0.0, self.params.clone());
        let bounds = bounds
            .into_iter()
            .zip(&steps)
            .map(|(b, &st)| {
                let over = Ctx { t: ink_core::eval::span_mid(st), t_span: Some(st), ..ctx.clone() };
                b.iter().zip(self.unknown_moves(&over)).map(|(x, y)| x + y).collect()
            })
            .collect();
        Chosen { steps: Steps { ends, unchecked, loose, basis }, rest, bounds }
    }

    /// `steps` with a step of its own (its middle there) where an actuator's
    /// effort, or its margin to either limit, turns back within 8% of its
    /// peak of its extreme among the span's ends and the steps' middles.
    /// A step whose extreme cannot fit the budget is loose.
    fn turns(
        &self,
        mut steps: Vec<(f64, f64)>,
        (a, b): (f64, f64),
        rest: &[bool],
        loose: &mut Vec<(f64, f64)>,
        unchecked: &[(f64, f64)],
        serving: Option<&Serve>,
    ) -> Vec<(f64, f64)> {
        use rayon::prelude::*;
        let mid = |s: (f64, f64)| ink_core::eval::span_mid(s);
        let ts: Vec<f64> = std::iter::once(a).chain(steps.iter().filter(|s| !unchecked.contains(s)).map(|&s| mid(s))).chain([b]).collect();
        let at: Vec<Arc<Moment>> = ts.par_iter().map(|&t| self.moment(t, rest)).collect();
        if at.iter().any(|m| m.unsolved.is_some()) {
            return steps;
        }
        #[derive(Clone, Copy)]
        enum Series {
            Effort(usize),
            Lower(usize),
            Upper(usize),
            Reaction(usize),
        }
        let series: Vec<_> = match serving {
            Some(_) => (0..6).map(Series::Reaction).collect(),
            None => (0..self.actuators.len())
                .filter(|&k| !self.actuators[k].coupled.is_empty() && !self.commanded(k))
                .flat_map(|k| [Series::Effort(k), Series::Lower(k), Series::Upper(k)])
                .collect(),
        };
        let read = |series, m: &Moment| match series {
            Series::Effort(k) => m.effort[k],
            Series::Lower(k) => m.effort[k] - m.lo[k],
            Series::Upper(k) => m.hi[k] - m.effort[k],
            Series::Reaction(k) => {
                let (f, mo) = serving.unwrap().wrench(m);
                if k < 3 { f[k] } else { mo[k - 3] }
            }
        };
        let mut brackets = Vec::new();
        for field in series {
            let v: Vec<f64> = at.iter().map(|m| read(field, m)).collect();
            if v.iter().any(|x| !x.is_finite()) {
                continue;
            }
            let (hi, lo) = (v.iter().copied().fold(f64::NEG_INFINITY, f64::max), v.iter().copied().fold(f64::INFINITY, f64::min));
            let margin = matches!(field, Series::Lower(_) | Series::Upper(_));
            let peak = match field {
                Series::Lower(k) | Series::Upper(k) => at.iter().map(|m| m.effort[k].abs()).fold(0.0, f64::max),
                _ => hi.abs().max(lo.abs()),
            };
            let near = 8.0 * SHARE * peak;
            for i in 1..v.len() - 1 {
                let up = !margin && v[i] > v[i - 1] && v[i] >= v[i + 1] && v[i] >= hi - near;
                let down = v[i] < v[i - 1] && v[i] <= v[i + 1] && v[i] <= lo + near;
                if up || down {
                    brackets.push((field, ts[i - 1], ts[i + 1], down));
                }
            }
        }
        let found: Vec<f64> = brackets
            .par_iter()
            .filter_map(|&(field, t0, t2, falling)| {
                let h = 1e-6 * (t2 - t0);
                let value = |t: f64| {
                    let m = self.moment(t, rest);
                    m.unsolved.is_none().then(|| read(field, &m))
                };
                let rate = |t: f64| match (value((t + h).max(t.next_up()).min(b)), value((t - h).min(t.next_down()).max(a))) {
                    (Some(p), Some(q)) => p - q,
                    _ => f64::NAN,
                };
                let t = crate::steps::turn_back(t0, t2, falling, rate);
                (t > t0 && t < t2 && value(t).is_some()).then_some(t)
            })
            .collect();
        for t in found {
            // Straddle an existing boundary too. Only its adjacent floats
            // are needed to make the measured middle the extreme itself.
            let (lo, hi) = (t.next_down().max(a), t.next_up().min(b));
            let first = steps.partition_point(|s| s.1 <= lo);
            let last = steps.partition_point(|s| s.0 < hi);
            let across = &steps[first..last];
            if across.is_empty() || across.iter().any(|s| mid(*s) == t || loose.contains(s) || unchecked.contains(s)) {
                continue;
            }
            // Reaction refinement may already use the whole adaptive grid.
            // These are isolated witnesses from its six components, at most
            // one per component per bracket; none triggers further splitting.
            if serving.is_none() && steps.len() + 3 - across.len() > crate::steps::MAX_STEPS {
                loose.extend_from_slice(across);
                continue;
            }
            let (left, right) = (across[0].0, across[across.len() - 1].1);
            steps.splice(first..last, [(left, lo), (lo, hi), (hi, right)].into_iter().filter(|s| s.0 < s.1));
        }
        steps
    }

    /// The moments over `span` where one of `fields` may jump (its rate over
    /// a step unbounded: a select or floor on t switching), each isolated by
    /// bisection to the resolution of t (a step between adjacent floats, or
    /// as near as its rate's bound narrows): two switches any distance apart
    /// in t are found apart, the moments between them a step of their own.
    /// At most MAX_STEPS / 8 of them; the steps holding the rest are judged
    /// as any step a limit may jump in is (Model::limit_undecided).
    fn jumps(&self, fields: &[&Num], span: (f64, f64)) -> Vec<(f64, f64)> {
        let most = crate::steps::MAX_STEPS / 8;
        let mut out = Vec::new();
        for n in fields {
            let mut open = vec![span];
            while let Some(s) = open.pop() {
                if out.len() >= most {
                    return out;
                }
                if self.rate(n, s).is_some() {
                    continue;
                }
                let m = ink_core::eval::span_mid(s);
                if !(s.0 < m && m < s.1) {
                    out.push(s);
                    continue;
                }
                open.extend([(m, s.1), (s.0, m)]);
            }
        }
        out
    }

    /// Smooth rating crossings left unresolved for this actuator alone.
    pub(crate) fn rating_gaps(&self, a: usize, span: (f64, f64), steps: &Steps) -> Vec<(f64, f64)> {
        if !(matches!(self.actuators[a].lo, Num::T(_)) || matches!(self.actuators[a].hi, Num::T(_))) {
            return Vec::new();
        }
        steps
            .ends
            .windows(2)
            .map(|w| (w[0], w[1]))
            // (A step with no float inside it is its two ends, each read: its rating cannot cross between them.)
            .filter(|s| !steps.unchecked.contains(s) && !steps.loose.contains(s) && s.0.next_up() < s.1)
            .filter(|&s| {
                let m = self.over(span, ink_core::eval::span_mid(s));
                m.margin(a) >= 0.0 && self.rating_dip(a, s, m.effort[a])
            })
            .collect()
    }

    /// A rating given directly as a field of time can cross the sampled
    /// effort even when both its sampled margins and its rate look flat.
    pub(crate) fn rating_dip(&self, k: usize, step: (f64, f64), effort: f64) -> bool {
        let act = &self.actuators[k];
        if !(matches!(act.lo, Num::T(_)) || matches!(act.hi, Num::T(_))) || !effort.is_finite() {
            return false;
        }
        let ctx = Ctx { t_span: Some(step), ..Ctx::new(ink_core::eval::span_mid(step), self.params.clone()) };
        (matches!(act.lo, Num::T(_)) && act.lo.over(&ctx).hi > effort) || (matches!(act.hi, Num::T(_)) && act.hi.over(&ctx).lo < effort)
    }

    /// Whether actuator k's limits may jump inside step s (a limit's rate
    /// over it unbounded) and their range over all of it, bounded by
    /// intervals, can cross `effort` (its effort, read at the step's
    /// middle): the middle's reading then decides nothing about them, and
    /// the step is loose. A limit that cannot cross there decides the step
    /// as its middle does; one that surely crosses fails at the middle too.
    pub(crate) fn limit_undecided(&self, k: usize, s: (f64, f64), effort: f64) -> bool {
        let act = &self.actuators[k];
        let mid = ink_core::eval::span_mid(s);
        let ctx = Ctx::new(mid, self.params.clone());
        let over = Ctx { t_span: Some(s), ..ctx };
        let crosses = |n: &Num, below: bool| {
            self.rate(n, s).is_none() && {
                let r = n.over(&over);
                if below { !(r.hi <= effort) } else { !(r.lo >= effort) }
            }
        };
        crosses(&act.lo, true) || crosses(&act.hi, false)
    }

    /// Which groups' spins must be resolved over `steps` (by index), and
    /// what is said of each read at rest. The efforts are measured at the
    /// span's ends and the steps' middles, every group read round: the
    /// actuators' where any is found, else the coordinates'. Each group's
    /// departure from round moves the effort along each coordinate by at
    /// most its bound over the step (Group::moves), an actuator's through
    /// the split. Groups are read at rest, least first, while for each
    /// effort the bounds of those read so, summed, stay within 1% of its
    /// peak (or its rounding), the share the steps keep the rest within. A
    /// group whose spin the budget cannot resolve is read at rest too where
    /// the efforts are all actuators' and the bounds, summed, stay under
    /// each claim's margin at every moment, so that carried as its margin
    /// they leave it decided (a fan held at speed, its weight off the axis
    /// all its effort, far from its limits); else it is resolved. A group
    /// no effort reads (spun by t, only reactions carry it) is read at rest.
    /// For a reaction a claim reads (`serve`), the groups it carries are
    /// read at rest where, least first, the bounds on its force and on its
    /// moment about its frame's origin (Group::shakes) of those so read,
    /// summed, stay within 1% of its peaks, or where their spin is past the
    /// budget; its claims decide with that band (Model::shaken).
    fn round_enough(&self, steps: &[(f64, f64)], (a, b): (f64, f64), serve: Option<usize>) -> (Vec<bool>, Vec<String>) {
        use rayon::prelude::*;
        let n = self.groups.len();
        let per = self.bounds(steps, &vec![true; n]);
        // Each moment, and the step it is in.
        let moments: Vec<(f64, usize)> = [(a, 0), (b, steps.len() - 1)].into_iter().chain(steps.iter().enumerate().map(|(j, &s)| (ink_core::eval::span_mid(s), j))).collect();
        let at: Vec<Arc<Moment>> = moments.par_iter().map(|&(t, _)| self.moment(t, &vec![true; n])).collect();
        // Each effort: its name, its limits' claim where it is an
        // actuator's, and at each moment its value, rounding, claim margin
        // and each group's bound on it.
        let found: Vec<usize> = (0..self.actuators.len()).filter(|&k| !self.actuators[k].coupled.is_empty()).collect();
        type At = (f64, f64, f64, Vec<f64>);
        let series: Vec<(String, bool, Vec<At>)> = match found.is_empty() {
            true => (0..self.coords.len())
                .map(|k| {
                    (self.coords[k].name.clone(), false, at.iter().zip(&moments).map(|(m, &(_, j))| (m.f[k], m.noise_f[k], 0.0, per[j].iter().map(|d| d[k]).collect())).collect())
                })
                .collect(),
            false => found
                .iter()
                .map(|&a| {
                    let on = |m: &Moment, d: &[f64]| m.pm[a].iter().zip(d).map(|(p, d)| p.abs() * d).sum();
                    (
                        self.actuators[a].path.clone(),
                        true,
                        at.iter().zip(&moments).map(|(m, &(_, j))| (m.effort[a], m.noise[a], m.margin(a), per[j].iter().map(|d| on(m, d)).collect())).collect(),
                    )
                })
                .collect(),
        };
        let allowed: Vec<f64> = series.iter().map(|(_, _, v)| v.iter().fold(0.0f64, |p, x| p.max(SHARE * x.0.abs()).max(x.1))).collect();
        // A bound that is none (NaN: an infinite one times a zero weight)
        // stays none, never 0.
        let most = |g: usize, s: usize| series[s].2.iter().map(|x| x.3[g]).fold(0.0, worst);
        let share = |g: usize| (0..series.len()).map(|s| most(g, s) / allowed[s]).fold(0.0, worst);
        // Whether resolving its spin is past the budget: steps keeping its
        // rotation within 1% (|R''| <= omega^2 + |alpha|) are at most
        // sqrt(8 1% / (omega^2 + |alpha|)) long, and halving lands each
        // between that and half of it.
        let past = |g: usize| {
            let need: f64 = steps
                .iter()
                .map(|&s| {
                    let over = Ctx { t: ink_core::eval::span_mid(s), t_span: Some(s), ..Ctx::new(0.0, self.params.clone()) };
                    let (w, al) = (ink_core::interval::mag(self.groups[g].rate.over(&over)), ink_core::interval::mag(self.groups[g].acc.over(&over)));
                    2.0 * (s.1 - s.0) * ((w * w + al) / (8.0 * SHARE)).sqrt()
                })
                .sum();
            !(need <= crate::steps::MAX_STEPS as f64)
        };
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&x, &y| share(x).total_cmp(&share(y)));
        // What those read at rest move each effort by: at most over the
        // moments, and at each.
        let mut used: Vec<f64> = vec![0.0; series.len()];
        let mut sum: Vec<Vec<f64>> = series.iter().map(|(_, _, v)| vec![0.0; v.len()]).collect();
        // One behind a coupling (Model::behind) is resolved too.
        let mut resolved: Vec<bool> = self.groups.iter().zip(&self.behind).map(|(g, b)| g.effort || *b).collect();
        let mut said = Vec::new();
        // A group through a solve has no bound over a step: resolved.
        for g in order.into_iter().filter(|&g| self.groups[g].effort && !self.groups[g].solved() && !self.behind[g]) {
            let small = (0..series.len()).all(|s| used[s] + most(g, s) <= allowed[s]);
            let decides = || series.iter().zip(&sum).all(|((_, claimed, v), sum)| *claimed && v.iter().zip(sum).all(|(x, u)| u + x.3[g] < x.2.abs()));
            if !(small || past(g) && decides()) {
                continue;
            }
            resolved[g] = false;
            for s in 0..series.len() {
                used[s] += most(g, s);
                for (u, x) in sum[s].iter_mut().zip(&series[s].2) {
                    *u += x.3[g];
                }
            }
            let gr = &self.groups[g];
            let (h, d2, d1) = self.departure(gr);
            let others: Vec<&str> = gr.members.iter().map(|m| m.0.as_str()).filter(|m| *m != gr.link).collect();
            let with = if others.is_empty() { String::new() } else { format!(" with {}", others.join(", ")) };
            let by: Vec<String> = (0..series.len()).map(|s| most(g, s)).map(|m| if m == 0.0 { "0".into() } else { format!("{m:.2e}") }).collect();
            said.push(format!("{}{with} ({h:.3}, {:.3}): {}", gr.link, d2.hypot(d1), by.join(", ")));
        }
        if !said.is_empty() {
            let names: Vec<&str> = series.iter().map(|(n, _, _)| n.as_str()).collect();
            said[0] = format!(
                "read with the spin at rest and their mass round about its axis, the bodies spinning with each link, by their departure from round (g mm \
                 of mass off its axis, g mm2 of inertia not round about it) and the most it moves each effort ({}, in its units) over the span: {}",
                names.join(", "),
                said[0]
            );
        }
        if let Some(r) = serve {
            // Its peaks at the moments, every group read round, and each
            // carried group's most on its force and on its moment about
            // its frame's origin over the steps.
            let ctx = Ctx::new(0.0, self.params.clone());
            let allowed = at.iter().fold((0.0f64, 0.0f64), |p, m| {
                let w = &m.reactions[r];
                let (f, mo) = w.local();
                (p.0.max(SHARE * norm(f)).max(w.noise.0), p.1.max(SHARE * norm(mo)).max(w.noise.1))
            });
            let most = |g: usize| {
                let gr = &self.groups[g];
                let dep = self.departure(gr);
                steps.iter().fold((0.0f64, 0.0f64), |p, &s| {
                    let over = Ctx { t: ink_core::eval::span_mid(s), t_span: Some(s), ..ctx.clone() };
                    let (f, mo) = gr.shakes(dep, &over, self.gravity);
                    let at = Ctx { t: ink_core::eval::span_mid(s), ..ctx.clone() };
                    let arm = norm(std::array::from_fn(|k| gr.frame[9 + k].at(&at) - self.reactions[r].frame[9 + k].at(&at)));
                    (worst(p.0, f), worst(p.1, mo + arm * f))
                })
            };
            // A group through a solve has no bound over a step: resolved.
            let mut carried: Vec<(usize, (f64, f64))> = (0..n).filter(|&g| self.carries(r, g) && !self.groups[g].solved() && !self.behind[g]).map(|g| (g, most(g))).collect();
            carried.sort_by(|x, y| (x.1.0 / allowed.0).max(x.1.1 / allowed.1).total_cmp(&(y.1.0 / allowed.0).max(y.1.1 / allowed.1)));
            let mut used = (0.0, 0.0);
            let mut held = Vec::new();
            for (g, (f, mo)) in carried {
                let small = used.0 + f <= allowed.0 && used.1 + mo <= allowed.1;
                resolved[g] = !(small || past(g));
                if !resolved[g] {
                    used = (used.0 + f, used.1 + mo);
                    held.push(format!("{} ({f:.4} N, {mo:.4} N mm)", self.groups[g].link));
                }
            }
            if !held.is_empty() {
                said.push(format!(
                    "for this reaction, the groups it carries read at rest, round, each with the most its departure from round moves the reaction's force and its \
                     moment about the reaction's frame's origin over the span: {}; its claims hold where they hold for every value within that band and fail \
                     at a moment's real value",
                    held.join(", ")
                ));
            }
        }
        (resolved, said)
    }

    /// `steps` with each that `wants` halved, while the budget lasts and it
    /// is longer than `shortest`; and whether any was.
    fn split(&self, steps: &[(f64, f64)], shortest: f64, wants: impl Fn((f64, f64)) -> bool) -> (Vec<(f64, f64)>, bool) {
        let mut out = Vec::with_capacity(steps.len());
        let mut room = crate::steps::MAX_STEPS.saturating_sub(steps.len());
        let mut split = false;
        for &s in steps {
            let m = ink_core::eval::span_mid(s);
            if room > 0 && s.1 - s.0 > shortest && s.0 < m && m < s.1 && wants(s) {
                out.extend([(s.0, m), (m, s.1)]);
                room -= 1;
                split = true;
            } else {
                out.push(s);
            }
        }
        (out, split)
    }
}
