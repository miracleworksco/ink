//! What makes an answer measured rather than proven: the fields it reads
//! read sampled data (an imported mesh, a baked field, a package's result)
//! or a solved value (a param a solve sets, not certified as a root over
//! exact inputs). Both show in the graph, so one walk finds them
//! (`Doc::inexact`). Nothing hides them at load time: a number made at
//! load time from either is an error (`Doc::load_number`).

use super::*;
use std::sync::Arc;

pub(super) const USE_IN_FIELDS: &str = "use it in fields (placements, shapes, claims), where answers that read it are measured";

/// What sampled data is, for a message: each source's kind, what it was
/// sampled from, and where.
pub(super) fn sampled_from(samples: &[Arc<crate::grid::Source>]) -> String {
    samples.iter().map(|s| format!("{} of {} at {}", s.kind, s.what, s.at)).collect::<Vec<_>>().join("; ")
}

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
        // What a sample found from samples was found from, after it, and
        // what that was found from (a bake of a bake), to the first.
        let mut samples = self.g.sources(&all);
        let mut grew = true;
        while grew {
            grew = false;
            for (found, from) in &self.found_from {
                if samples.iter().any(|x| Arc::ptr_eq(x, found)) {
                    for f in from {
                        if !samples.iter().any(|x| Arc::ptr_eq(x, f)) {
                            samples.push(f.clone());
                            grew = true;
                        }
                    }
                }
            }
        }
        Inexact { samples, solved, unsettled, notes }
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
        let what = if !x.solved.is_empty() { "a solved value".to_string() } else { format!("sampled data ({})", sampled_from(&x.samples)) };
        Err(self.err_hint(e, format!("this reads {what}, which cannot shape the model at load time (a list, a condition, a count, a form's setting)"), Some(USE_IN_FIELDS.into())))
    }

    /// The box `(bounds f)` finds at load time of a field `what` that reads
    /// sampled data (`from`): a bound of the samples, so measured, not
    /// proven, of what they were sampled from. Its numbers stay samples: read
    /// from a grid of their own whose source says what they were found from,
    /// so what is built from them is measured as what reads `from` is, and
    /// an `if` or a count made from them is refused (`exact_at_load`). Held
    /// as f32, as samples are, rounded outward.
    pub(super) fn measured_box(&mut self, e: &Sexp, what: &Sexp, lo: V3, hi: V3, from: &[Arc<crate::grid::Source>]) -> Value {
        let out = |v: f64, up: bool| {
            let s = v as f32;
            match (s as f64).partial_cmp(&v) {
                Some(std::cmp::Ordering::Less) if up => s.next_up(),
                Some(std::cmp::Ordering::Greater) if !up => s.next_down(),
                _ => s,
            }
        };
        let six: Vec<f32> = (0..6).map(|i| if i < 3 { out(lo[i], false) } else { out(hi[i - 3], true) }).collect();
        // Six points along x (two along y and z, as a grid has at least two
        // per axis), number i at x = i.
        let values = (0..24).map(|l| six[l % 6]).collect();
        let mut text = what.to_string();
        if text.chars().count() > 80 {
            text = text.chars().take(77).chain("...".chars()).collect();
        }
        let source = crate::grid::Source {
            what: text,
            at: format!("{}:{}", self.file, e.line),
            kind: "bounds",
            res: from.iter().map(|s| s.res).fold(0.0, f64::max),
            band: None,
            measured: None,
            note: Some(format!("a box found from sampled data ({}): no error bound to what that was sampled from", sampled_from(from))),
            unit: "",
            span: None,
            component: None,
            undecided: undecided_from(from),
        };
        let data = Arc::new(crate::grid::GridData::new([0.0; 3], 1.0, [6, 2, 2], values).holding().with_source(source));
        if let Some(s) = data.source() {
            self.found_from.push((s.clone(), from.to_vec()));
        }
        let zero = self.g.constant(0.0);
        let mut n: Vec<Value> = (0..6)
            .map(|i| {
                let x = self.g.constant(i as f64);
                Value::Field(self.g.grid(data.clone(), [x, zero, zero]), None)
            })
            .collect();
        let hi = n.split_off(3);
        Value::list(vec![Value::list(n), Value::list(hi)])
    }

    /// Record that these expressions were evaluated at load time, their
    /// values kept as numbers or samples: their params are materialized
    /// there (`materialized`: they cannot be ranged).
    pub(super) fn bake(&self, nodes: &[NodeId]) {
        self.pin(&self.params_read(nodes));
    }

    /// Record that these params' values became numbers at load time
    /// (`bake`).
    pub(super) fn pin(&self, ks: &[usize]) {
        {
            let mut m = self.materialized.borrow_mut();
            for &k in ks {
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

/// Why no claim is decided on what is found from samples `from`: the first
/// of theirs (grid::Source::undecided), so a grid baked, redistanced or
/// bounded from a result its form gave no band keeps that refusal as its
/// own, whatever it is called.
pub(crate) fn undecided_from(from: &[Arc<crate::grid::Source>]) -> Option<String> {
    from.iter().find_map(|s| s.undecided.clone())
}
