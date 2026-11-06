//! What makes an answer measured rather than proven: the fields it reads
//! read sampled data (an imported mesh, a baked field, a package's result)
//! or a solved value (a param a solve sets, not certified as a root over
//! exact inputs). Both show in the graph, so one walk finds them
//! (`Doc::inexact`). Nothing hides them at load time: a number made at
//! load time from either is an error (`Doc::load_number`).

use super::*;
use std::sync::Arc;

pub(super) const USE_IN_FIELDS: &str = "use it in fields (placements, shapes, claims), where answers that read it are measured";

/// The inexact inputs a value reads.
#[derive(Clone, Debug, Default)]
pub struct Inexact {
    /// The sampled fields it reads, each once: what each was sampled from,
    /// and at what spacing.
    pub samples: Vec<Arc<crate::grid::Source>>,
    /// The solved params it reads that are not certified roots (by index).
    pub solved: Vec<u32>,
    /// Those of them whose provider said, where it last set them, that
    /// they fall short of the accuracy it states, and why.
    pub unsettled: Vec<(u32, String)>,
    /// What their providers say of them beyond their names
    /// (ParamProvider::note), each once.
    pub notes: Vec<String>,
}

impl Inexact {
    /// Whether an answer reading it is measured.
    pub fn any(&self) -> bool {
        !self.samples.is_empty() || !self.solved.is_empty()
    }
}

impl Doc {
    /// The inexact inputs of the fields at `roots`: the sampled fields they
    /// read, and the solved params they read, through the inputs of the
    /// solves setting them. `certified`: solved params a view derives as
    /// certified roots (ModelView::certified), exact if their solve's
    /// inputs are.
    pub fn inexact(&self, roots: &[NodeId], certified: &[u32]) -> Inexact {
        let mut all = roots.to_vec();
        let mut todo = self.params_read(roots);
        let mut seen = HashSet::new();
        let (mut solved, mut unsettled, mut notes) = (Vec::new(), Vec::new(), Vec::new());
        while let Some(k) = todo.pop() {
            if !seen.insert(k) {
                continue;
            }
            for s in self.solves.iter().filter(|s| s.owns.contains(&k)) {
                if !certified.contains(&(k as u32)) {
                    solved.push(k as u32);
                    if let Some(why) = self.provided.borrow().get(&s.ext).and_then(|p| p.unsettled.clone()) {
                        unsettled.push((k as u32, why));
                    }
                    if let Some(n) = self.ext[s.ext].form.provider(self).and_then(|p| p.note()).filter(|n| !notes.contains(n)) {
                        notes.push(n);
                    }
                }
                let inputs = self.form_inputs(&self.ext[s.ext]);
                todo.extend(self.params_read(&inputs));
                all.extend(inputs);
            }
        }
        solved.sort_unstable();
        unsettled.sort_unstable();
        Inexact { samples: self.g.sources(&all), solved, unsettled, notes }
    }

    /// A field as a number at load time (`load_value`), where the model's
    /// structure needs one (a `for` list, an `if`, an index, a count): an
    /// error if it reads a solved or sampled value, which would hide in the
    /// number and leave what is built from it looking exact.
    pub(super) fn load_number(&self, e: &Sexp, f: NodeId) -> R<Option<f64>> {
        self.exact_at_load(e, f)?;
        Ok(self.load_value(f))
    }

    /// An error if `f` reads a solved or sampled value (see `load_number`).
    pub(super) fn exact_at_load(&self, e: &Sexp, f: NodeId) -> R<()> {
        let x = self.inexact(&[f], &[]);
        if !x.any() {
            return Ok(());
        }
        let what = if !x.solved.is_empty() { "a solved value" } else { "sampled data" };
        Err(self.err_hint(e, format!("this reads {what}, which cannot shape the model at load time (a list, a condition, a count, a box)"), Some(USE_IN_FIELDS.into())))
    }

    /// Record that these expressions were evaluated at load time, their
    /// values kept as numbers or samples: their params are materialized
    /// there (`materialized`: they cannot be ranged).
    pub(super) fn bake(&self, nodes: &[NodeId]) {
        let ks = self.params_read(nodes);
        {
            let mut m = self.materialized.borrow_mut();
            for &k in &ks {
                m.entry(k).or_insert_with(|| self.form_at.clone());
            }
        }
        self.baked.borrow_mut().extend(ks);
    }

    /// The fields a check of the claim evaluates: its own (Claim::fields);
    /// a package metric's, what its form read (it is solved from those).
    pub fn claim_roots(&self, claim: &Claim) -> Vec<NodeId> {
        match claim {
            Claim::Metric { metric, .. } => self.metric_roots(metric),
            c => c.fields(),
        }
    }

    /// The fields a package's metric is solved from: what its form read,
    /// or its shape.
    pub fn metric_roots(&self, metric: &FormMetric) -> Vec<NodeId> {
        if let Some(f) = metric.shape {
            return vec![f];
        }
        self.ext_form(Some(metric.package), &metric.form).map(|e| self.form_inputs(e)).unwrap_or_default()
    }
}
