//! An empirical expression remains a field. A derived zero carries the
//! condition under which it applies and the estimate's standing through
//! arithmetic, without freezing the expression at one parameter value.
use super::{Doc, Env, Form, Param, ParamProvider, R, Reads, Refused, Sexp, Value};
use crate::eval::{Ctx, eval_interval, eval_point};
use crate::graph::{NodeId, Var};
use crate::tape::Tape;
use ink_core::interval::I;
use std::collections::HashSet;
use std::sync::Arc;

#[derive(Clone)]
pub(super) struct Estimate {
    owns: [usize; 1],
    root: NodeId,
    condition: Tape,
    on_t: bool,
    reason: String,
}

pub(super) fn read(doc: &mut Doc, env: &mut Env, e: &Sexp, args: &[Sexp]) -> R<Value> {
    let [value, condition, reason] = args else { return Err(doc.err(e, "expected (estimate-when value condition reason)")) };
    let c = doc.field(condition, env)?;
    let Value::Str(why) = doc.eval(reason, env)? else { return Err(doc.err(reason, "an estimate needs a reason for refusing")) };
    let Value::Field(f, unit) = doc.eval(value, env)? else { return Err(doc.err(value, "an estimate needs a number expression")) };
    let guard = marker(doc, e, c, why)?;
    Ok(Value::Field(doc.g.add(f, guard), unit))
}

fn marker(doc: &mut Doc, site: &Sexp, root: NodeId, reason: String) -> R<NodeId> {
    if doc.g.depends_on(root, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
        return Err(doc.err(site, "an estimate's condition must be independent of position"));
    }
    let k = doc.add_param(Param { name: "conditional estimate".into(), default: f64::NAN, tol: None });
    let estimate = Estimate { owns: [k], root, condition: Tape::new(&doc.g, root), on_t: doc.g.depends_on(root, |v| v == Var::T), reason };
    doc.solved_params.insert(k);
    // Only the condition is provided. The value stays a field, including over ranges.
    doc.keep_form("std", "estimate-when", site, Arc::new(estimate), Reads { fields: vec![root], ..Reads::default() })?;
    // A marker: the estimate less itself, or 0 times it, still reads it (Graph::marker).
    Ok(doc.g.marker(k as u32))
}

/// A time substitution must reach the condition hidden behind each marker.
/// Visit outer guards first to find nested ones, then rebuild in dependency order.
pub(super) fn retime(doc: &mut Doc, site: &Sexp, value: &Value, map: &mut Vec<(Var, NodeId)>) -> R<()> {
    let mut roots = Vec::new();
    value.fields_into(&mut roots, false);
    let mut needed: HashSet<_> = doc.params_read(&roots).into_iter().collect();
    let mut guards = Vec::new();
    for form in doc.ext.iter().rev() {
        if let Some(guard) = Estimate::of(form)
            && needed.contains(&guard.owns[0])
        {
            needed.extend(doc.params_read(&[guard.root]));
            guards.push((guard.owns[0], guard.root, guard.reason.clone()));
        }
    }
    for (k, root, reason) in guards.into_iter().rev() {
        let at_time = doc.g.subst(root, map);
        if at_time != root {
            map.push((Var::Param(k as u32), marker(doc, site, at_time, reason)?));
        }
    }
    Ok(())
}

impl Estimate {
    pub(super) fn of(form: &super::ExtForm) -> Option<&Self> {
        (form.form.as_ref() as &dyn std::any::Any).downcast_ref()
    }

    /// The marker is identically zero wherever this whole scope is applicable.
    /// Proving applicability never certifies the empirical field it marks.
    pub(super) fn provide_over(&self, ctx: &mut Ctx) -> Result<(), String> {
        let ok = eval_interval(&self.condition, ctx, [I::point(0.0); 3], None);
        if ok.lo.is_finite() && ok.hi.is_finite() && (ok.lo > 0.0 || ok.hi < 0.0) {
            ctx.params[self.owns[0]] = 0.0;
            Ok(())
        } else {
            Err(self.reason.clone())
        }
    }
}

/// What is said where an estimate's condition (`why`: its reason) is not shown at every moment of
/// `span`, the times what reads it reads it at.
pub(super) fn not_over(why: &str, span: (f64, f64)) -> String {
    let sig = crate::claims::sig;
    format!("an estimate it reads applies only under a condition, which is read at every moment of t {} to {} and not shown at each: {why}", sig(span.0), sig(span.1))
}

impl Doc {
    /// Of solve `i`, where it is an estimate's guard (estimate-when): whether its condition keeps one sign
    /// over all of `over` (its params' ranges, a band's among them), so the guard is 0 at every value there
    /// and what it reads moves nothing through it: a formula, not a solve. Err: its reason, where the
    /// condition may change sign there (the estimate does not apply over all of it). None: not a guard.
    pub(crate) fn guard_over(&self, i: usize, over: &Ctx) -> Option<Result<(), String>> {
        let guard = Estimate::of(&self.ext[self.solves[i].ext])?;
        Some(guard.provide_over(&mut over.clone()))
    }

    /// Of solve `i`, where it is an estimate's guard: its condition shown at `ctx`'s moment over its
    /// params' ranges (interval bounds), and where readers read it on clocks of their own, over those
    /// times too (`provide_one`). The one check every way of reaching a moment over ranges makes (want,
    /// set_moment, scope_spans), so they agree. Ok: its marker 0 there. None: not a guard.
    pub(crate) fn show_guard(&self, i: usize, ctx: &Ctx) -> Option<Result<(), String>> {
        let s = &self.solves[i];
        let guard = Estimate::of(&self.ext[s.ext])?;
        let at = Ctx { t_span: None, ..ctx.clone() };
        Some(guard.provide_over(&mut at.clone()).and_then(|()| match self.guard_clock(i).filter(|_| s.on_t) {
            None => Ok(()),
            Some(Err(r)) => Err(r.why),
            Some(Ok(span)) => {
                let mut over = Ctx { t: ink_core::eval::span_mid(span), t_span: Some(span), ..at };
                guard.provide_over(&mut over).map_err(|why| not_over(&why, span))
            }
        }))
    }

    /// Whether every estimate the fields `roots` read under a condition on t (estimate-when) applies at
    /// every moment of `span`, the params at `params`: what reads them at other times than the one its
    /// params were filled at (a model stepping the whole motion) asks this first, since a condition
    /// filled at one moment says nothing of another. Err: the first that is not shown to, with its reason.
    pub fn estimates_apply_over(&self, roots: &[NodeId], params: &[f64], span: (f64, f64)) -> Result<(), String> {
        for i in self.solves_read(roots) {
            if let Some(guard) = Estimate::of(&self.ext[self.solves[i].ext]).filter(|g| g.on_t) {
                let mut over = Ctx { t_span: Some(span), ..Ctx::new(ink_core::eval::span_mid(span), params.to_vec()) };
                guard.provide_over(&mut over).map_err(|why| not_over(&why, span))?;
            }
        }
        Ok(())
    }
}

impl Form for Estimate {
    fn name(&self) -> &str {
        "conditional estimate"
    }
    fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(self.clone()))
    }
}

impl ParamProvider for Estimate {
    fn name(&self) -> &str {
        Form::name(self)
    }
    fn owns(&self) -> &[usize] {
        &self.owns
    }
    fn uses_t(&self) -> bool {
        self.on_t
    }
    fn warm_starts(&self) -> bool {
        false
    }
    fn on_demand(&self) -> bool {
        true
    }
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
        let ok = eval_point(&self.condition, &Ctx::new(t, params.to_vec()), [0.; 3]);
        if !ok.is_finite() || ok == 0.0 {
            return Err(Refused::from(self.reason.clone()));
        }
        params[self.owns[0]] = 0.0;
        Ok(None)
    }
}

/// Differentiating a sum drops its constant term. The estimate marker is
/// numerically zero, but its applicability and standing also govern the rate.
pub(super) fn keep_conditions(doc: &mut Doc, root: NodeId, mut rate: NodeId) -> NodeId {
    let mut estimates = doc.forms::<Estimate>().peekable();
    if estimates.peek().is_none() {
        return rate;
    }
    let inputs: HashSet<_> = doc.params_read(&[root]).into_iter().collect();
    let kept: HashSet<_> = doc.params_read(&[rate]).into_iter().collect();
    let missing: Vec<_> = estimates.map(|e| e.owns[0]).filter(|k| inputs.contains(k) && !kept.contains(k)).collect();
    for k in missing {
        let marker = doc.g.var(Var::Param(k as u32));
        rate = doc.g.add(rate, marker);
    }
    rate
}
