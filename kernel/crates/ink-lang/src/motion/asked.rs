use super::*;

/// How a claim over a motion was checked.
#[derive(Clone, Debug)]
pub struct Stepped {
    pub motion: Motion,
    /// Its answer is a proof (t ranged over the span, or asked once);
    /// else it is measured at the steps' middles.
    pub bound: bool,
    /// t was a ranged param of the one search (`ranged`).
    pub ranged: bool,
    /// The body only moves rigidly, and what is measured was asked once in
    /// its own frame (carried::over_span): a thickness at one moment, an
    /// integral's or a mean's own frame's integrals carried by its
    /// placement with t ranged.
    pub rigid: bool,
    /// The uniform step a measured claim was checked at (its configured
    /// size, however many steps were checked before a disproof).
    pub step: Option<f64>,
    /// Measured at steps a provider of what it reads chose, and how
    /// (ParamProvider::steps), not at uniform ones.
    pub basis: Option<String>,
    /// Measured, why t ranged did not decide it (`ranged`), and where it
    /// searched and was left undecided, its own unknown answer (its hint,
    /// undecided cell, the params it ranged).
    pub left: Option<String>,
    pub left_answer: Option<Box<Check>>,
    /// Not settled: it holds at each moment read, and what lies between them is not covered (`between`),
    /// in the words every answer from moments says it in (`not_followed`).
    pub between: Option<String>,
}

impl Stepped {
    /// Where it is not settled for what lies between its moments (`between`): that, and what would
    /// settle it: for a contact or overlap, which no range of t shows, a claim that can hold at every
    /// moment; else why t ranged left the steps (`left`).
    pub fn unsettled(&self) -> Option<String> {
        // Read at a provider's steps, some of which it does not follow.
        let m = &self.motion;
        if self.between.is_none() && self.basis.is_some() && m.status == ink_kernel::ranged::Status::Unknown {
            let mut gaps: Vec<(f64, f64)> = m.unchecked.iter().chain(&m.loose).copied().collect();
            gaps.sort_by(|a, b| a.0.total_cmp(&b.0));
            return not_followed(MODEL_GAP, &gaps);
        }
        let between = self.between.as_ref()?;
        Some(match (self.left_answer.as_ref().is_some_and(|a| a.note == Some(claims::AT_ONE_VALUE)), &self.left) {
            (true, _) => format!(
                "{between}. A contact found at one moment says nothing of the next, so no step between them is shown: claim what can hold at every moment (a clearance, an overlap, disjoint, a value), or ask it at rest with why the motion can't change it"
            ),
            (false, Some(left)) => format!("{between}. There, {left}"),
            (false, None) => between.clone(),
        })
    }
    pub(super) fn of(motion: Motion, bound: bool, step: Option<f64>) -> Stepped {
        Stepped { motion, bound, ranged: false, rigid: false, step, basis: None, left: None, left_answer: None, between: None }
    }
    /// Asked of a body that only moves rigidly (`rigid`).
    pub(crate) fn rigid(motion: Motion) -> Stepped {
        Stepped { rigid: true, ..Stepped::of(motion, true, None) }
    }
}

/// A kernel search a check waits on, built: it reads its tapes and its
/// context alone, no view, so the searches of many checks run side by
/// side (`Asked`): an emptiness with t ranged (motion::ranged), as the
/// tape, its context, the region and the proof's options.
pub struct Search(pub(super) std::sync::Arc<Tape>, pub(super) Ctx, pub(super) (V3, V3), pub(super) ink_kernel::ranged::ProofOptions);

/// What may replace an answer given already (a moving pair's sampled one:
/// `over_span`), begun once every answer is given, with the time left.
pub struct Upgrade {
    pub start: Box<dyn FnOnce(&mut ModelView) -> Result<Later, InkError>>,
}

/// What an Upgrade's start gives: the answer replacing the one given (None:
/// it stands), or the search that decides it and what it then gives.
pub enum Later {
    Now(Option<Check>),
    Wait(Search, Box<dyn FnOnce(&mut ModelView, Found) -> Result<Option<Check>, InkError>>),
}

/// What a search found.
pub type Found = ink_kernel::ranged::Proof;

impl Search {
    /// Its question's key in the answer cache (None: not cached): the one
    /// ranged::prove stores its proof under, so a search another in flight
    /// asks (every copy of a sub-assembly) waits for that one's proof.
    pub fn key(&self) -> Option<u128> {
        let Search(tape, ctx, (lo, hi), opts) = self;
        ink_kernel::ranged::prove_key(tape, None, ctx, *lo, *hi, opts)
    }

    pub fn run(&self) -> Found {
        let Search(tape, ctx, (lo, hi), opts) = self;
        ink_kernel::ranged::prove(tape, ctx, *lo, *hi, opts)
    }
}

/// What a check that may wait on a search gives: its answer, or the
/// search and what is done with what it finds (with the view again: a
/// clash's moment, an exact question, another search). Its question is
/// built in the view (which adds to the graph), one check after another;
/// its search needs no view, so a caller with many (ink interfere's moving
/// pairs) runs their searches side by side while it builds the next.
pub enum Asked<T> {
    Done(T),
    Wait(Search, Box<dyn FnOnce(&mut ModelView, Found) -> Result<Asked<T>, InkError>>),
}

impl<T: 'static> Asked<T> {
    /// Its answer, its searches run here, one after another.
    pub fn wait(mut self, v: &mut ModelView) -> Result<T, InkError> {
        loop {
            match self {
                Asked::Done(x) => return Ok(x),
                Asked::Wait(s, then) => self = then(v, s.run())?,
            }
        }
    }

    /// `f` of its answer, once it has one.
    pub fn then<U: 'static>(self, v: &mut ModelView, f: impl FnOnce(&mut ModelView, T) -> Result<Asked<U>, InkError> + 'static) -> Result<Asked<U>, InkError> {
        match self {
            Asked::Done(x) => f(v, x),
            Asked::Wait(s, then) => Ok(Asked::Wait(s, Box::new(move |v, found| then(v, found)?.then(v, f)))),
        }
    }
}
