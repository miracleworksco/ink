//! Loads and compiles mechanisms for a moment, a motion, or a simulation.

use crate::supports::supports_near;
use crate::{Ctx, Model, R};
use ink_core::graph::{NodeId, Var};
use ink_lang::view::ModelView;
use ink_lang::{Doc, InkError, Reload};
use std::sync::Arc;

/// The paths a whole answer reports reactions for: the target, then each
/// joint under it, by instance path.
pub fn joints(doc: &mut Doc, target: &str) -> Result<Vec<String>, InkError> {
    let items = doc.items(target).map_err(InkError::msg)?;
    let mut out = vec![target.to_string()];
    for l in items.iter().flat_map(|it| &it.links).filter(|l| l.joint.is_some()) {
        if !out.contains(&l.path) {
            out.push(l.path.clone());
        }
    }
    Ok(out)
}

/// The model of `target` in the file `reload` loads again with its
/// coordinates, its masses measured at moment t (a deforming body's is
/// that moment's), with the reactions of `reactions` (`claimed`: claims
/// read them: Model::compile).
pub fn build(reload: &Reload, target: &str, gravity: [f64; 3], t: f64, reactions: &[String], claimed: bool) -> R<Model> {
    build_as(reload, target, gravity, t, reactions, claimed, false, &[])
}

/// A joint's drive's rate (deg/s on a turn, mm/s on a slide) as a field of params and t: the written motion's
/// own, which needs no mass and no gravity. None where it reads a value worked out at a run (a simulated
/// coordinate, a solve): then it is the run's to give.
pub fn drive_rate(doc: &mut Doc, joint: &ink_lang::Joint) -> Option<NodeId> {
    let rate = doc.derivative(joint.drive, Var::T);
    doc.solves_read(&[rate]).is_empty().then_some(rate)
}

/// The bounds of a drive's rate (`drive_rate`) over every t of `span`, by the kernel's interval reading of
/// the formula: (0, 0) proves the joint at rest throughout, whatever happens between any moments sampled.
pub fn rate_bounds(doc: &Doc, rate: NodeId, ctx: &Ctx, span: (f64, f64)) -> (f64, f64) {
    let scope = Ctx { t: 0.5 * (span.0 + span.1), t_span: Some(span), ..ctx.clone() };
    let i = ink_core::eval::eval_interval_scope(&ink_core::tape::Tape::new(&doc.g, rate), &scope, [ink_core::interval::I::point(0.0); 3], None);
    (i.lo, i.hi)
}

/// `build`, leaving out the moving bodies with no mass and naming them (`Model::compile_kinematic`): for one
/// who reads the motion whatever the masses, and each effort only where no such body would change it.
pub fn build_kinematic(reload: &Reload, target: &str, gravity: [f64; 3], t: f64, reactions: &[String]) -> R<Model> {
    build_as(reload, target, gravity, t, reactions, false, true, &[])
}

pub(super) fn build_as(reload: &Reload, target: &str, gravity: [f64; 3], t: f64, reactions: &[String], claimed: bool, kinematic: bool, held: &[(String, Vec<String>)]) -> R<Model> {
    let doc = reload.with_coordinates(target)?;
    let span = doc.t_span;
    if doc.assembly(target).is_none() {
        return Err(InkError::msg(format!("'{target}' is not a top-level assembly of the file: an effort is of the mechanism an assembly's instances pass arguments to")));
    }
    let (ctx, unsolved) = doc.ctx_at(&reload.overrides, t).map_err(InkError::msg)?;
    let mut v = ModelView::new(doc, ctx);
    v.overrides = reload.overrides.clone();
    v.unsolved = unsolved;
    supports_near(&mut v, target, held)?;
    estimates_apply(&mut v, target, span)?;
    let mut model = if kinematic { Model::compile_partial(&mut v, target, gravity, reactions, claimed)? } else { Model::compile(&mut v, target, gravity, reactions, claimed)? };
    model.span = span;
    model.twin = Some((reload.clone(), target.to_string(), Arc::default()));
    Ok(model)
}

/// A model steps the whole motion with its params as they are at one moment, and an estimate's
/// condition on t (estimate-when) is filled at that moment alone: so every one that `target`'s motion,
/// its joints' limits or a form of this package reads is shown over all of `span` first, or there is
/// no model.
fn estimates_apply(v: &mut ModelView, target: &str, span: Option<(f64, f64)>) -> R<()> {
    let Some(span) = span else { return Ok(()) };
    let mut roots: Vec<NodeId> = Vec::new();
    for it in v.doc.items(target)? {
        for l in &it.links {
            roots.extend(l.at.r.iter().flatten().chain(&l.at.t));
            if let Some(j) = &l.joint {
                roots.push(j.drive);
                roots.extend(j.limits.iter().chain(&j.speed_limits).flat_map(|&(a, b)| [a, b]));
            }
        }
    }
    for e in v.doc.ext.iter().filter(|e| e.package == "dynamics") {
        roots.extend(v.doc.form_inputs(e));
    }
    roots.extend(v.doc.coords.iter().map(|c| c.value));
    v.doc.estimates_apply_over(&roots, &v.ctx.params, span).map_err(InkError::msg)
}

/// The model of `target` in the file `reload` loads again with its
/// coordinates built as states (Doc::with_states): the mechanism at any
/// position, speed and acceleration along them (Model::state).
pub fn build_states(reload: &Reload, target: &str, gravity: [f64; 3]) -> R<Model> {
    states_model(reload.with_states(target)?, reload, target, gravity)
}

/// The model of `target` in a load with states (`with_states`) that a
/// simulation of its springs steps: its springs released
/// (Model::compile_springs).
pub(super) fn springs_model(doc: Doc, reload: &Reload, target: &str, gravity: [f64; 3]) -> R<Model> {
    if doc.assembly(target).is_none() {
        return Err(InkError::msg(format!("'{target}' is not a top-level assembly of the file: a simulation is of the mechanism an assembly's instances pass arguments to")));
    }
    let span = doc.t_span;
    let (ctx, unsolved) = doc.ctx_at(&reload.overrides, 0.0).map_err(InkError::msg)?;
    let mut v = ModelView::new(doc, ctx);
    v.overrides = reload.overrides.clone();
    v.unsolved = unsolved;
    estimates_apply(&mut v, target, span)?;
    let mut model = Model::compile_springs(&mut v, target, gravity)?;
    model.span = span;
    Ok(model)
}

/// `build_states` of a load with states already made.
pub(super) fn states_model(doc: Doc, reload: &Reload, target: &str, gravity: [f64; 3]) -> R<Model> {
    if doc.assembly(target).is_none() {
        return Err(InkError::msg(format!("'{target}' is not a top-level assembly of the file: an effort is of the mechanism an assembly's instances pass arguments to")));
    }
    let span = doc.t_span;
    let (ctx, unsolved) = doc.ctx_at(&reload.overrides, 0.0).map_err(InkError::msg)?;
    let mut v = ModelView::new(doc, ctx);
    v.overrides = reload.overrides.clone();
    v.unsolved = unsolved;
    estimates_apply(&mut v, target, span)?;
    let mut model = Model::compile(&mut v, target, gravity, &[], false)?;
    model.span = span;
    Ok(model)
}
