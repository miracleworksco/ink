//! Applied loads and diagnostics for loads that the mechanism does not carry.

use crate::{R, clause, model};
use ink_core::graph::{NodeId, Var};
use ink_lang::sexp::Sexp;
use ink_lang::{Doc, Form, Quantity, Reader, Value, Xf};
use std::sync::Arc;

/// `(load NAME PATH (force fx fy fz) [(torque tx ty tz)] [(at FRAME)])`: a
/// force (N) and torque (N mm), in the world's axes, on what the instance
/// path PATH places, acting at FRAME's origin (a frame in the world;
/// default PATH's placement's origin), moving with it; fields of params and
/// t, which may read frames. The file run's alone count.
pub struct Load {
    pub name: String,
    pub path: String,
    pub force: [NodeId; 3],
    pub torque: [NodeId; 3],
    pub at: Option<Xf>,
    pub run: bool,
    /// Of an imported file's: the bodies its PATH places in that file (their definitions' names), by which a
    /// design that places one of them too is told the load is not on it (`loads_not_applied`).
    places: Vec<String>,
    pub(super) form: Sexp,
}

/// A load an imported file states on its own root, on a part this design places too. It is that file's
/// machine's, not the part's: the design does not carry it, and an answer about the design must say so.
pub struct NotApplied {
    pub name: String,
    /// The path as its file writes it, and its line there.
    pub path: String,
    pub line: u32,
    /// Where this design places what it loads there.
    pub here: Vec<String>,
}

impl NotApplied {
    pub fn said(&self) -> String {
        format!(
            "the load '{}' (line {} of an imported file, on its own {}) is not applied: a file's top-level load is on that file's root, which this design does not place, \
             though it places what the load is on ({}). A load that goes wherever its assembly is placed is written inside it: (assembly NAME ... (load {} PART (force ...)))",
            self.name,
            self.line,
            self.path,
            self.here.join(", "),
            self.name
        )
    }
}

/// Every load the design applies under `target`, as it is written: the instance path it is on and its
/// force's components in the world's axes (fields of params and t). For one who adds them up apart from
/// the model.
pub fn loads_written(doc: &mut Doc, target: &str) -> Vec<(String, [NodeId; 3])> {
    let mut on: Vec<(String, [NodeId; 3])> = doc.forms::<Load>().filter(|l| l.run).map(|l| (l.path.clone(), l.force)).collect();
    on.extend(clause::loads_under(doc, target).into_iter().map(|l| (l.path, l.force)));
    on
}

/// The loads under `target` on a part that is on no joint, each said: nothing that moves or is rated carries
/// such a load. It goes into the ground through a part the model takes as fixed, so a requirement of
/// ratings checks nothing against it.
pub fn loads_to_ground(doc: &mut Doc, target: &str) -> Vec<String> {
    let mut on: Vec<(String, String)> = doc.forms::<Load>().filter(|l| l.run).map(|l| (l.name.clone(), l.path.clone())).collect();
    on.extend(clause::loads_under(doc, target).into_iter().map(|l| (l.name, l.path)));
    let under = |path: &str| path == target || path.strip_prefix(target).is_some_and(|rest| rest.starts_with('/'));
    on.into_iter()
        .filter(|(_, path)| under(path) && model::riding(doc, target, path).is_none())
        .map(|(name, path)| {
            format!(
                "the load '{name}' is on {path}, which is on no joint: it goes into the ground through a part taken as fixed, and no part here carries it. Put {path} \
                 on what holds it ((frame-of ITS-PART FRAME) in its WHERE, a drum on a motor's shaft), or a study of {path} answers for it"
            )
        })
        .collect()
}

/// The imported files' top-level loads on parts `target` places too (`NotApplied`).
pub fn loads_not_applied(doc: &mut Doc, target: &str) -> Vec<NotApplied> {
    let left: Vec<(String, String, u32, Vec<String>)> =
        doc.forms::<Load>().filter(|l| !l.run && !l.places.is_empty()).map(|l| (l.name.clone(), l.path.clone(), l.form.line, l.places.clone())).collect();
    if left.is_empty() {
        return Vec::new();
    }
    let Ok(items) = doc.items(target) else { return Vec::new() };
    let mut out = Vec::new();
    for (name, path, line, places) in left {
        let here: Vec<String> = items.iter().filter(|it| it.instance.as_ref().is_some_and(|(body, _)| places.contains(body))).map(|it| it.label.clone()).collect();
        if !here.is_empty() {
            out.push(NotApplied { name, path, line, here });
        }
    }
    out
}

impl Form for Load {
    fn name(&self) -> &str {
        &self.name
    }
}

pub(super) const LOAD_USAGE: &str = "(load NAME PATH (force fx fy fz) [(torque tx ty tz)] [(at FRAME)]): a force in N and a torque in N mm, in the world's axes, on \
    what the instance path PATH places (root/arm/fore), acting at FRAME's origin (a frame in the world: (at (frame-of root/arm/fore tip)) moves with the \
    part, (at (move 0 0 150)) stays where it is; default PATH's own origin); components are fields of params and t, and may read frames (a spring: along the line between two parts)";

/// `(load NAME PATH (force fx fy fz) [(torque tx ty tz)] [(at FRAME)])`.
pub(super) fn read_load(r: &mut Reader, form: &Sexp, args: &[Sexp]) -> R<Arc<dyn Form>> {
    let [name, path, rest @ ..] = args else { return Err(r.err_hint(form, "bad load", Some(LOAD_USAGE.into()))) };
    let (name, path) = (r.name(name)?, r.name(path)?);
    // A load's name is its own: the same form read twice would be summed, a force doubled with nothing said.
    let run = r.in_file_run();
    if let Some(line) = r.forms::<Load>().into_iter().find(|l| l.name == name && l.run == run).map(|l| l.form.line) {
        return Err(r.err_hint(
            form,
            format!("the load '{name}' is already stated, at line {line}: a second one would be added to it"),
            Some("one load, one name: delete one, or give a second force on the same part a name of its own".into()),
        ));
    }
    let zero = r.graph().constant(0.0);
    let (mut force, mut torque, mut at) = (None, None, None);
    for c in rest {
        let l = c.list().unwrap_or(&[]);
        match (l.first().and_then(|h| h.sym()), l.get(1..).unwrap_or(&[])) {
            (Some(k @ ("force" | "torque")), [a, b, cc]) => {
                let (q, what) = if k == "force" { (Quantity::Force, "force") } else { (Quantity::Torque, "torque") };
                let mut v = [zero; 3];
                for (i, e) in [a, b, cc].into_iter().enumerate() {
                    v[i] = r.field_of(e, q, what)?;
                    if r.graph().depends_on(v[i], |x| matches!(x, Var::X | Var::Y | Var::Z)) {
                        return Err(r.err_hint(e, "a load's components are fields of params and t, not of x y z", Some(LOAD_USAGE.into())));
                    }
                }
                let slot = if k == "force" { &mut force } else { &mut torque };
                if slot.replace(v).is_some() {
                    return Err(r.err(c, format!("one ({k} ...) per load")));
                }
            }
            // One frame (`(at (frame-of root/arm tip))`), or the frame the clause itself writes, as an
            // instance's WHERE does: `(at (move 0 0 150))`.
            (Some("at"), [f, ..]) => match (l.len() == 2).then(|| r.value(f)).and_then(Result::ok) {
                Some(Value::Frame(x)) => at = Some(x),
                _ => match r.value(c)? {
                    Value::Frame(x) => at = Some(x),
                    other => return Err(r.err_hint(f, format!("a load's (at FRAME) takes a frame, not a {}", other.kind()), Some(LOAD_USAGE.into()))),
                },
            },
            _ => return Err(r.err_hint(c, "bad load clause", Some(LOAD_USAGE.into()))),
        }
    }
    if force.is_none() && torque.is_none() {
        return Err(r.err_hint(form, "a load needs (force fx fy fz) or (torque tx ty tz)", Some(LOAD_USAGE.into())));
    }
    // An imported file's: what it loads there, by the bodies' definitions.
    let places = match run {
        true => Vec::new(),
        false => r.items(&path).map(|items| items.iter().filter_map(|it| it.instance.as_ref().map(|(body, _)| body.clone())).collect()).unwrap_or_default(),
    };
    Ok(Arc::new(Load { name, path, force: force.unwrap_or([zero; 3]), torque: torque.unwrap_or([zero; 3]), at, run, places, form: form.clone() }))
}
