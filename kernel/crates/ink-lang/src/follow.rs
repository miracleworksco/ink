//! What the forms that provide params (ParamProvider: a simulated motion,
//! a measurement) set at a moment: each sets its params itself. None gives
//! equations, so nothing is a root followed along t, and what one sets is
//! measured at each moment, never bounded over a span of t. The solver
//! that was (`(solve ...)`), with this file as it followed and certified
//! a root, is at commit b0bfe6a306.

use crate::Refused;
use crate::ext::ParamProvider;
use std::sync::Arc;

/// A document's providers in their order (Doc::solves).
pub struct Follow {
    solves: Vec<Solve>,
}

struct Solve {
    p: Arc<dyn ParamProvider>,
    /// What it is, for messages ("effort 'rotor'").
    label: String,
    /// It is on t (ext::Solve::on_t).
    on_t: bool,
    /// The params it rests on that no provider on t sets.
    fixed: Vec<usize>,
}

impl Follow {
    /// The providers in their order: each, the params it reads
    /// (Solve::reads), whether it is on t (Solve::on_t), and what it is.
    pub fn new(solves: Vec<(Arc<dyn ParamProvider>, Vec<usize>, bool, String)>) -> Follow {
        let on_t: Vec<bool> = solves.iter().map(|s| s.2).collect();
        let owner = |k: usize| solves.iter().position(|s| s.0.owns().contains(&k));
        let mut out: Vec<Solve> = Vec::with_capacity(solves.len());
        for (i, (p, reads, t, label)) in solves.iter().enumerate() {
            let mut upstream: Vec<usize> = reads.iter().filter_map(|&k| owner(k)).filter(|&j| j < i && on_t[j]).collect();
            upstream.sort_unstable();
            upstream.dedup();
            let mut fixed: Vec<usize> = reads.iter().copied().filter(|&k| owner(k).is_none_or(|j| !on_t[j])).collect();
            for &j in &upstream {
                fixed.extend(&out[j].fixed);
            }
            fixed.sort_unstable();
            fixed.dedup();
            out.push(Solve { p: p.clone(), label: label.clone(), on_t: *t, fixed });
        }
        Follow { solves: out }
    }

    /// Whether provider `i` is a root followed along t (none is: no form
    /// gives equations), the params it sets, and what it is.
    pub fn followed(&self, i: usize) -> (bool, &[usize], &str) {
        (false, self.solves[i].p.owns(), &self.solves[i].label)
    }

    /// Set provider `i`'s params at moment t from their values in `params`
    /// (a simulation's state from the last one's).
    pub fn solve_near(&self, i: usize, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
        self.solves[i].p.provide(params, t)
    }

    /// Set provider `i`'s params at moment t, as it sets them from params'
    /// values. What falls short of the accuracy it states, why; Err: why it
    /// sets nothing.
    pub fn provide(&self, i: usize, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
        self.provide_with_bounds(i, params, t, &[], &mut Vec::new())
    }

    pub(crate) fn provide_with_bounds(
        &self,
        i: usize,
        params: &mut [f64],
        t: f64,
        inputs: &[crate::ProvidedParam],
        bounds: &mut Vec<crate::ProvidedParam>,
    ) -> Result<Option<String>, Refused> {
        let s = &self.solves[i];
        if s.on_t && !t.is_finite() {
            return Err(Refused::from(format!("{} reads t, and t = {t} is not a moment", s.label)));
        }
        s.p.provide_with_bounds(params, t, inputs, bounds)
    }

    /// The params provider `i` rests on that no provider on t sets.
    pub fn rests_on(&self, i: usize) -> &[usize] {
        &self.solves[i].fixed
    }
}
