//! Provides the tensions of the lines in a coupling.

use crate::providers::{Shared, unseparated};
use ink_lang::{Doc, Form, ParamProvider, Refused, Steps};
use std::sync::Arc;

/// The form keeping one tension's provider: the coupling's ends, and the
/// param it sets.
#[derive(Clone)]
pub(super) struct Tensing {
    pub(super) shared: Arc<Shared>,
    pub(super) a: String,
    pub(super) b: String,
    pub(super) owns: Vec<usize>,
}

impl Form for Tensing {
    fn name(&self) -> &str {
        "tension"
    }

    fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(self.clone()))
    }
}

impl ParamProvider for Tensing {
    fn name(&self) -> &str {
        "a coupling's lines' tensions along the motion (the balance of the passive joints it reaches; (tension A B))"
    }

    fn owns(&self) -> &[usize] {
        &self.owns
    }

    fn uses_t(&self) -> bool {
        true
    }

    fn on_demand(&self) -> bool {
        true
    }

    fn warm_starts(&self) -> bool {
        false
    }

    /// Its lines' tensions at t, every body as measured (a group whose
    /// departure from round could move them is never read at rest:
    /// Model::behind); none where they, or a balance they rest on, are not
    /// found (Model::tension_refusal).
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
        let model = &self.shared.model()?.0;
        if let Some(sim) = self.shared.springs()
            && sim?.model().sprung(&self.a, &self.b)
        {
            return Err(Refused::from(format!(
                "{}~{} is a coupling of springs in the file's simulation: its lines' tensions are their springs' (ink effort, simulated.springs.spans), \
                 not the rigid model's balance",
                self.a, self.b
            )));
        }
        let m = model.at(t);
        let k = model.coupling_index(&self.a, &self.b).ok_or_else(|| format!("no coupling relates {} and {}", self.a, self.b))?;
        if let Some(why) = model.tension_refusal_strict(&m, k) {
            return Err(Refused::from(why));
        }
        let lines = model.couplings()[k].1.clone();
        for (k, v) in self.owns.iter().zip(&m.tensions[lines]) {
            params[*k] = *v;
        }
        let unsettled: Vec<&str> = model.unsettled.iter().chain(&model.unsettled_carried).map(String::as_str).collect();
        Ok((!unsettled.is_empty()).then(|| format!("masses it rests on are not measured to 2%: {}", unsettled.join("; "))))
    }

    fn reads(&self) -> Option<Vec<usize>> {
        let (model, _) = self.shared.model().ok()?;
        Some(model.reads.iter().chain(&model.reads_carried).filter_map(|n| self.shared.names.iter().position(|m| m == n)).collect())
    }

    /// The efforts' steps, which resolve every group a tension rests on.
    fn steps(&self, _params: &[f64], span: (f64, f64)) -> Option<Result<Steps, Refused>> {
        Some(self.shared.model().map(|(m, _)| m.steps(span)))
    }

    /// Where a fixed-effort drive's required effort differs from its rating
    /// within the rounding of its load, its tensions are known to the
    /// shortfall's share of their rounding either way
    /// (Model::tension_shortfall).
    fn error(&self, _params: &[f64], t: f64) -> Option<Vec<f64>> {
        let (model, _) = self.shared.model().ok()?;
        let k = model.coupling_index(&self.a, &self.b)?;
        model.tension_shortfall(&model.at(t), k)
    }

    /// What the balances giving the tensions leave unseparated.
    fn note(&self) -> Option<String> {
        let (model, _) = self.shared.model.get()?.as_ref().ok()?;
        unseparated(model)
    }
}
