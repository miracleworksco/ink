//! Providers for actuator effort and RMS effort claims.

use crate::providers::{Shared, Sprung, spring_effort, spring_steps};
use crate::{model, simulate};
use ink_lang::{Doc, Form, ParamProvider, Refused, Steps};
use std::sync::Arc;

pub(super) const RMS_USAGE: &str = "(effort-rms PATH): the RMS over the file's (t-span A B) of the effort the actuator PATH supplies (a joint with effort-limits, by \
    its instance path, root/arm/s/output), a torque a turn or a force a slide, as ink effort's rms: compare it with a record's continuous rating, \
    its unit said, e.g. (value (effort-rms root/arm/s/output) < (* 1N*m (get (get m.spec \"gearing\") \"rated-torque-Nm\")))";

/// The form keeping one actuator's effort provider (its index among the
/// shared model's actuators' paths, and the param it sets).
#[derive(Clone)]
pub(super) struct Effort {
    pub(super) shared: Arc<Shared>,
    pub(super) actuator: usize,
    pub(super) owns: [usize; 1],
}

impl Form for Effort {
    fn name(&self) -> &str {
        "effort"
    }

    fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(self.clone()))
    }
}

impl ParamProvider for Effort {
    fn name(&self) -> &str {
        "the actuator's effort along the motion (ink effort; its effort-limits claim)"
    }

    fn owns(&self) -> &[usize] {
        &self.owns
    }

    fn uses_t(&self) -> bool {
        true
    }

    /// Computed only while a claim reading an effort is checked.
    fn on_demand(&self) -> bool {
        true
    }

    fn warm_starts(&self) -> bool {
        false
    }

    /// What `ink check` asks of it (ParamProvider::refusal): why it is
    /// refused at every moment: the model not built, its load running
    /// through a simulation's springs (a sprung effort's, as `provide`), or
    /// bodies left out (Model::leaves_out). Not why it is refused at one
    /// moment of the motion: that needs the whole motion stepped
    /// (Model::read), which only its claim reads, and its claim's answer
    /// says it.
    fn refusal(&self, params: &[f64], t: f64) -> Option<Refused> {
        match self.shared.sprung(self.actuator) {
            Err(why) => Some(why),
            Ok(Sprung::Spring { .. }) => self.provide(&mut params.to_vec(), t).err(),
            Ok(Sprung::No) => match self.shared.model() {
                Err(why) => Some(why),
                Ok((model, index)) => model.missing_mass(index[self.actuator]).or_else(|| model.left_out(index[self.actuator])),
            },
        }
    }

    /// Its effort at t (Model::read): inside the file's span as the steps
    /// over it read it, the groups read there with their spin at rest read
    /// round, the effort within `bound` of the real one at any angle of
    /// their spins. None where it is not this motion's (Model::refusal), or
    /// where that bound, taken as its claim's margin, leaves the claim
    /// undecided: the claim holds at a moment where effort +- bound lies
    /// within the limits, and fails where the effort is past one by more
    /// than it, or where the effort at that moment's own angles is (a
    /// witness, the moment's real effort: that one is given).
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
        // Its joint released for a simulation's springs: what it supplies
        // is its spring's force, against its deviation.
        if let Sprung::Spring { i, k } = self.shared.sprung(self.actuator)? {
            let sim = self.shared.springs().expect("a simulation of springs")?;
            params[self.owns[0]] = spring_effort(&sim, i, k, t)?;
            return Ok(None);
        }
        let (model, index) = self.shared.model()?;
        let a = index[self.actuator];
        let m = model.read(t);
        if let Some(why) = model.refused(a, &m) {
            return Err(why);
        }
        if m.fixed_margin_unresolved(a) {
            return Err(Refused::from(format!("the required effort {} differs from the fixed rating {} only within rounding of its load", m.effort[a], m.lo[a])));
        }
        let (x, b) = (m.margin(a), m.bound[a]);
        let mut effort = m.effort[a];
        // (An effort that is its command, clamped to its limits, lies within them whatever the bodies' departure.)
        if b > 0.0 && -b <= x && x < b && m.command[a].is_none() {
            let real = model.at(t);
            if !(real.margin(a) < 0.0) {
                return Err(Refused::from(format!(
                    "its effort, {:.6} with the bodies read with their spin at rest read round, is {x:.6} inside its limits [{:.6}, {:.6}] where the \
                     departure from round of those bodies can move it by {b:.6} at any angle of their spins (motion.basis): undecided here; a \
                     margin beyond that decides",
                    m.effort[a], m.lo[a], m.hi[a]
                )));
            }
            effort = real.effort[a];
        }
        params[self.owns[0]] = effort;
        Ok((!model.unsettled.is_empty()).then(|| format!("masses it rests on are not measured to 2%: {}", model.unsettled.join("; "))))
    }

    /// Its spring's force is an estimate known either way to within its
    /// stiffness times its following error's band at the span's end (the
    /// widest: the band grows with the largest deviation so far and the
    /// periods stepped) and what the masses are known to move it by at t.
    fn error(&self, _params: &[f64], t: f64) -> Option<Vec<f64>> {
        let Ok(Sprung::Spring { i, k }) = self.shared.sprung(self.actuator) else {
            // What its bodies' measures are known to (Moment::mass).
            let (model, index) = self.shared.model().ok()?;
            let mass = model.read(t).mass[index[self.actuator]];
            return (mass > 0.0).then(|| vec![mass]);
        };
        let sim = self.shared.springs()?.ok()?;
        let end = self.shared.model().ok()?.0.span?.1;
        Some(vec![k * (sim.band(i, end) + sim.mass_band(i, t)[0])])
    }

    /// What the model reads, as the document numbers it; none where it
    /// cannot be built (its refusal is then the answer).
    fn reads(&self) -> Option<Vec<usize>> {
        let (model, _) = self.shared.model().ok()?;
        let mut names: Vec<&String> = model.reads.iter().collect();
        // Its spring's force reads what the simulation does.
        let sim = self.shared.springs().and_then(Result::ok);
        if let Some(sim) = &sim {
            names.extend(&sim.model().reads);
        }
        let mut out: Vec<usize> = names.into_iter().filter_map(|n| self.shared.names.iter().position(|m| m == n)).collect();
        out.sort_unstable();
        out.dedup();
        Some(out)
    }

    /// Where its couplings' loss rests on tensions that leave bodies out:
    /// the claim is decided only where its margin passes what that can be
    /// off by (Model::refusal).
    fn note(&self) -> Option<String> {
        if let Ok(Sprung::Spring { .. }) = self.shared.sprung(self.actuator) {
            return Some(
                "its joint is released for the springs of the file's simulation: what it supplies is its spring's force, its stiffness times its \
                 following error (a stepper's holding torque against its load angle), an estimate known either way to within its stiffness times \
                 the following error's band at the span's end (measured, no bound: the simulation's note says how) and what the masses move it \
                 by, its limits read at the commanded speed"
                    .into(),
            );
        }
        let (model, index) = self.shared.model.get()?.as_ref().ok()?;
        let loss = model.loss_leaves_out(index[self.actuator]).then(|| {
            format!(
                "the loss its couplings charge to this effort rests on tensions that leave out what {} exert (moving with the joints whose balances \
                 give them, carried by no actuated or balancing joint): it can be off by up to its worst efficiency's share of what they need along \
                 its coordinates, and the claim is unknown where its margin is within that",
                model.unseparated.join(", ")
            )
        });
        let notes: Vec<String> = loss.into_iter().chain(model.moved_by_t_note()).collect();
        (!notes.is_empty()).then(|| notes.join("; "))
    }

    /// The model's steps over `span`, and, where groups are read there with
    /// their spin at rest, the most they move this actuator's effort: its
    /// claim's margin.
    fn steps(&self, _params: &[f64], span: (f64, f64)) -> Option<Result<Steps, Refused>> {
        // Its spring's force: at the efforts' steps and where it turns back.
        if let Ok(Sprung::Spring { i, k }) = self.shared.sprung(self.actuator) {
            let read = |sim: Arc<simulate::Simulation>| self.shared.model().map(|(m, index)| spring_steps(m, index[self.actuator], &sim, i, k, span));
            return Some(self.shared.springs().expect("a simulation of springs").map_err(Refused::from).and_then(read));
        }
        Some(self.shared.model().and_then(|(m, index)| {
            let a = index[self.actuator];
            if let Some(why) = m.missing_mass(a) {
                return Err(why);
            }
            let mut s = match m.bounded(a) {
                Ok(plan) => bounded(m, a, span, &plan, &m.steps(span)),
                Err(why) => {
                    let mut s = m.steps(span);
                    s.basis = format!("{}; this effort is not bounded between its moments ({why}): an estimate there", s.basis);
                    estimated(m, a, span, &mut s, None);
                    s
                }
            };
            let middles = s.ends.windows(2).map(|w| ink_core::eval::span_mid((w[0], w[1])));
            let most = [span.0, span.1].into_iter().chain(middles).map(|t| m.read(t).bound[a]).fold(0.0, f64::max);
            if most > 0.0 {
                s.basis = format!(
                    "{}; this effort is read with them round, within {most:.6} of the real one at any angle of their spins (at most, over the steps): \
                     its claim holds at a moment where effort +- that bound lies within its limits, fails where the effort is past one by more, and \
                     is undecided between",
                    s.basis
                );
            }
            Ok(s)
        }))
    }
}

/// Actuator a's steps over a span with its effort bounded over each
/// (Model::bound_steps): a claim holding at them holds between them too,
/// where every step was settled; those not, read at their moments as before
/// (`estimated`).
fn bounded(m: &model::Model, a: usize, span: (f64, f64), plan: &model::Plan, steps: &Steps) -> Steps {
    let b = m.bound_steps(a, plan, steps);
    let mut s = b.steps;
    let how = "this effort is bounded over each step, not only read at its moments: interval arithmetic rounded outward over every moment of the \
               step at once (each moving body's acceleration and partial velocities, its mass, centre and inertia anywhere in their proven bounds, \
               the loads, the drives' constant split, its limits over the step's speeds), a step whose bound crosses a limit halved while that may \
               settle it";
    let tail = match (b.failed, b.unsettled.len(), s.unchecked.is_empty()) {
        (true, 0, _) => format!(" ({} halvings): it crosses a limit at a moment read, past what its masses are known to", b.halved),
        (true, n, _) => {
            format!(" ({} halvings): a moment read crossed a limit past what its masses are known to; {n} unsettled steps read at their moments alone: an estimate there", b.halved)
        }
        (false, 0, true) => format!(" ({} halvings, {} steps): within its limits at every moment of the span, not only at those read", b.halved, b.proven),
        (false, 0, false) => format!(" ({} halvings): within its limits at every moment of its {} steps, not only at those read, the unchecked ones aside", b.halved, b.proven),
        (false, n, _) => format!(
            " ({} halvings): within its limits throughout {} steps; {n} it did not settle (its bound still crossing a limit at the budget, as a \
             turn too fast for a step's bound makes it, or at the step's middle within what the masses or the rounding are known to of one), \
             read at their moments alone: an estimate there",
            b.halved, b.proven
        ),
    };
    s.basis = format!("{}; {how}{tail}", s.basis);
    if !b.unsettled.is_empty() {
        estimated(m, a, span, &mut s, Some(&b.unsettled));
    }
    s
}

/// Actuator a's steps over a span read at their moments alone (its effort
/// not bounded over them): an estimate between them, a step where it lies
/// inside its limits by less than a peak between moments may add decided
/// by nothing, a rating crossing it inside a step loose; of those `only`
/// holds, where it is given.
fn estimated(m: &model::Model, a: usize, span: (f64, f64), s: &mut Steps, only: Option<&[(f64, f64)]>) {
    let mine = |step: &(f64, f64)| only.is_none_or(|o| o.contains(step));
    let gaps: Vec<(f64, f64)> = m.rating_gaps(a, span, s).into_iter().filter(mine).collect();
    s.loose.extend(gaps);
    // Between the moments read nothing is bounded: the steps are
    // chosen so that a peak is read at most OFF_LINE of itself low.
    // A step whose effort lies inside its limits by less than that
    // decides nothing (a 10.00001 N peak between moments held a
    // rating of 10.000005 the moments read lay 0.000005 inside).
    let reads: Vec<((f64, f64), Vec<f64>)> = s
        .ends
        .windows(2)
        .enumerate()
        .map(|(i, w)| {
            let at = [ink_core::eval::span_mid((w[0], w[1]))].into_iter().chain((i == 0).then_some(span.0)).chain((i + 2 == s.ends.len()).then_some(span.1));
            ((w[0], w[1]), at.map(|t| m.read(t)).flat_map(|r| [r.effort[a].abs(), r.margin(a)]).collect())
        })
        .collect();
    let between = model::OFF_LINE * reads.iter().flat_map(|r| r.1.iter().step_by(2)).fold(0.0f64, |p, &x| p.max(x));
    let tight: Vec<(f64, f64)> = reads
        .iter()
        .filter(|(step, r)| mine(step) && r.iter().skip(1).step_by(2).any(|&x| x >= 0.0 && x < between) && !s.unchecked.contains(step) && !s.loose.contains(step))
        .map(|r| r.0)
        .collect();
    if !tight.is_empty() {
        s.basis = format!(
            "{}; this effort lies inside its limits by less than {between:.6}, the {}% of its peak a peak between the moments read may lie over \
             them, in {} of its steps read at their moments alone: those decide nothing",
            s.basis,
            100.0 * model::OFF_LINE,
            tight.len()
        );
        s.loose.extend(tight);
    }
}

/// The form keeping one RMS effort's provider: its actuator's index among
/// the shared model's actuators' paths, and the param it sets.
#[derive(Clone)]
pub(super) struct Rmsing {
    pub(super) shared: Arc<Shared>,
    pub(super) actuator: usize,
    pub(super) owns: [usize; 1],
}

impl Form for Rmsing {
    fn name(&self) -> &str {
        "effort-rms"
    }

    fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(self.clone()))
    }
}

impl ParamProvider for Rmsing {
    fn name(&self) -> &str {
        "the RMS of an actuator's effort over the file's span (ink effort's rms; (effort-rms PATH))"
    }

    fn owns(&self) -> &[usize] {
        &self.owns
    }

    fn uses_t(&self) -> bool {
        false
    }

    fn on_demand(&self) -> bool {
        true
    }

    fn warm_starts(&self) -> bool {
        false
    }

    /// What `ink check` asks of it (ParamProvider::refusal): what `provide`
    /// refuses before stepping the motion (Effort::refusal).
    fn refusal(&self, params: &[f64], t: f64) -> Option<Refused> {
        let (model, index) = match self.shared.model() {
            Ok(m) => m,
            Err(why) => return Some(why),
        };
        match (model.span, self.shared.sprung(self.actuator)) {
            (None, _) | (_, Ok(Sprung::Spring { .. })) => self.provide(&mut params.to_vec(), t).err(),
            (_, Err(why)) => Some(why),
            (Some(_), Ok(Sprung::No)) => model.missing_mass(index[self.actuator]).or_else(|| model.left_out(index[self.actuator])),
        }
    }

    /// Model::rms over the file's span, once per model.
    fn provide(&self, params: &mut [f64], _t: f64) -> Result<Option<String>, Refused> {
        let (model, index) = self.shared.model()?;
        let a = index[self.actuator];
        if let Some(why) = model.missing_mass(a) {
            return Err(why);
        }
        let span = model.span.ok_or_else(|| format!("an RMS effort is over the file's span: it declares no (t-span A B) ({RMS_USAGE})"))?;
        // Its joint released for a simulation's springs: its spring's force's.
        if let Sprung::Spring { i, k } = self.shared.sprung(self.actuator)? {
            let sim = self.shared.springs().expect("a simulation of springs")?;
            params[self.owns[0]] = k * sim.rms(i, span)?;
            return Ok(None);
        }
        params[self.owns[0]] = model.rms(a, span)?.0;
        Ok((!model.unsettled.is_empty()).then(|| format!("masses it rests on are not measured to 2%: {}", model.unsettled.join("; "))))
    }

    fn reads(&self) -> Option<Vec<usize>> {
        let (model, _) = self.shared.model().ok()?;
        Some(model.reads.iter().filter_map(|n| self.shared.names.iter().position(|m| m == n)).collect())
    }

    /// Where its steps read groups round: the most that moves the effort at
    /// a step, which moves an RMS by no more; a spring's force's, its
    /// stiffness times its following error's band at the span's end and
    /// the most the masses move it over the span (an RMS moves no more than
    /// the most its samples do).
    fn error(&self, _: &[f64], _t: f64) -> Option<Vec<f64>> {
        if let Ok(Sprung::Spring { i, k }) = self.shared.sprung(self.actuator) {
            let sim = self.shared.springs()?.ok()?;
            return Some(vec![k * sim.widest(i, self.shared.model().ok()?.0.span?.1)]);
        }
        let (model, index) = self.shared.model().ok()?;
        let band = model.rms(index[self.actuator], model.span?).ok()?.1;
        (band > 0.0).then(|| vec![band])
    }

    fn note(&self) -> Option<String> {
        if let Ok(Sprung::Spring { .. }) = self.shared.sprung(self.actuator) {
            return Some(
                "its joint is released for the springs of the file's simulation: its RMS effort is its spring's force's, its stiffness times its \
                 following error's RMS (Gauss-Legendre over each step of the simulation), an estimate known either way to within its stiffness \
                 times the following error's band at the span's end (measured, no bound: the simulation's note says how)"
                    .into(),
            );
        }
        Some(
            "an RMS effort is an estimate: the effort at the middles of the efforts' steps over the file's span (the midpoint rule), each step \
              keeping it within 4% of its peak of the line through its neighbours"
                .into(),
        )
    }
}
