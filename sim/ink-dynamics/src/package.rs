//! Reads dynamics forms and registers the providers and claims they request.

use crate::effort::{Effort, RMS_USAGE, Rmsing};
use crate::loads::{LOAD_USAGE, Load, loads_not_applied, read_load};
use crate::providers::Shared;
use crate::reactions::{Reacting, compiled, divided};
use crate::simulation::{COMMAND_USAGE, Command, FOLLOWING_USAGE, Following, SIMULATE_USAGE, SPRINGS_USAGE, Simulated, Springs, run_of, springs_said, taylor};
use crate::supports::{Support, carries, load_frame};
use crate::tension::Tensing;
use crate::{Dynamics, R, clause, command, coupling, joint, split};
use ink_core::graph::{NodeId, Var};
use ink_lang::sexp::Sexp;
use ink_lang::view::ModelView;
use ink_lang::{AssertDef, Claim, Cmp, Doc, Form, InkError, Package, ParamProvider, Quantity, Reader, Refused, Value, Xf};
use std::sync::{Arc, OnceLock};

/// `(gravity gx gy gz)`: the world's gravity, m/s^2; the file run's alone.
pub struct Gravity {
    pub g: [f64; 3],
    /// Written in the file run, not an import.
    pub run: bool,
    at: Sexp,
}

impl Form for Gravity {
    fn name(&self) -> &str {
        "gravity"
    }
}

const USAGE: &str = "(gravity gx gy gz): the world's gravity in m/s2, once, in the file run, e.g. (gravity 0 0 -9.81m/s2) for z up: an \
    effort, a reaction and every load the model works out need it, and a joint's (effort-limits LO HI) is claimed only with it";

/// The run file's gravity (m/s^2), if it declares one.
pub fn gravity(doc: &Doc) -> Option<[f64; 3]> {
    doc.forms::<Gravity>().find(|g| g.run).map(|g| g.g)
}

/// The error for a file without gravity, for an answer that needs it.
pub fn no_gravity() -> InkError {
    InkError {
        hint: Some(USAGE.into()),
        ..InkError::msg(
            "the file declares no (gravity gx gy gz): an effort needs the world's gravity, and none is assumed, since nothing in a file says which way is up. One line: \
             (gravity 0 0 -9.81m/s2) where z is up, (gravity 0 -9.81m/s2 0) where y is",
        )
    }
}

const SUPPORTS_USAGE: &str = "(supports PATH S ...): the parts S (instance paths: bearings, carriages, a thrust bearing) hold what the joint at instance \
    path PATH carries, so each one's load is its share of that joint's reaction, by statics, and the joint's own part carries none of it. For a body \
    held by parts that are not its joint (a shaft in two bearings): a part on a joint (a carriage on its slide) needs no line, and is never its own \
    support; a support touches what it holds. At the top \
    level, paths from root ((supports root/h/motor/rotor root/h/ba root/h/bb)); in an assembly, its own instances' ((supports motor/rotor ba bb))";

/// `(supports PATH S ...)` at the top level: named by its path.
struct Supports {
    path: String,
}

impl Form for Supports {
    fn name(&self) -> &str {
        &self.path
    }
}

/// `(supports PATH S ...)` in an assembly: its paths as that assembly's instances'.
struct AssemblySupports {
    path: String,
    by: Vec<String>,
}

/// `(head item ...)`, written where `at` is.
fn listed(at: &Sexp, head: &str, items: &[Sexp]) -> Sexp {
    let sym = Sexp { k: ink_lang::sexp::SK::Sym(head.into()), line: at.line, col: at.col };
    Sexp { k: ink_lang::sexp::SK::List(std::iter::once(sym).chain(items.iter().cloned()).collect::<Vec<_>>().into()), line: at.line, col: at.col }
}

/// `(reaction PATH [FRAME])` in an expression: its six numbers, provided
/// at each moment, and the frame they are in (None: PATH's own); or
/// `(reaction PATH (by S ...))`, six for each support S, each in the frame
/// its share is in (None: its placement).
struct Want {
    path: String,
    frame: Option<Xf>,
    by: Vec<(String, Option<Xf>)>,
    owns: Vec<usize>,
    run: bool,
    at: Sexp,
    /// Said by `(supports ...)`, read by no expression: the division is named, its shares facts for ratings.
    named: bool,
}

impl Form for Want {
    fn name(&self) -> &str {
        &self.path
    }
}

const REACTION_USAGE: &str = "(reaction PATH [FRAME]): what the rest of the model exerts on the bodies the instance path PATH carries (those under it \
    and those placed on them, (frame-of ...) in their WHERE), (vec fx fy fz mx my mz), force in N and moment in N mm about FRAME's origin in its axes \
    (a frame in the world: (at) the world's); default PATH's joint's frame (its WHERE moved by its motion), else its placement; root's, the world's. \
    (reaction PATH (by S ...)): its share each support S carries (a placed part whose record's load-ratings say what it carries, a bearing, a \
    carriage, at the frame its model declares, load), (vec W ...), each W in S's load frame about its origin, or in F where S is written (S F)";

/// `(tension A B)` in an expression: the params it is provided as, the
/// tension of each line of the coupling relating A and B at each moment.
struct Pull {
    a: String,
    b: String,
    owns: Vec<usize>,
    run: bool,
    at: Sexp,
}

impl Form for Pull {
    fn name(&self) -> &str {
        &self.a
    }
}

const TENSION_USAGE: &str = "(tension A B): the tension of each line of the coupling relating the instance paths A and B, N, at each moment, \
    (vec T ...) in its lines' order (a belt's two spans, a gear's two flanks, a screw's push and pull): their relation's (couple ...), as std's \
    belt-drive, chain-drive, belt-span, gear-mesh, bevel-mesh, rack-mesh and lead-screw give one; its root assembly placed above";

/// `(effort-rms PATH)` in an expression: the param it is provided as, the
/// RMS over the file's span of the effort of the actuator PATH.
struct Rms {
    path: String,
    owns: usize,
    run: bool,
    at: Sexp,
}

impl Form for Rms {
    fn name(&self) -> &str {
        &self.path
    }
}

/// Whether `path` names an instance `items` place.
fn placed(items: &[ink_lang::Item], path: &str) -> bool {
    items.iter().any(|it| it.links.iter().any(|l| l.path == path))
}

/// The actuators placed under `root`, by instance path, with their joints.
fn actuators(items: &[ink_lang::Item]) -> Vec<(String, ink_lang::Joint)> {
    let mut out: Vec<(String, ink_lang::Joint)> = Vec::new();
    for l in items.iter().flat_map(|it| &it.links) {
        if let Some(j) = l.joint.as_ref().filter(|j| joint::effort(j).is_some())
            && !out.iter().any(|(p, _)| *p == l.path)
        {
            out.push((l.path.clone(), j.clone()));
        }
    }
    out
}

/// The effort claims a file without gravity makes none of: each actuator
/// under its root assembly, `effort-limits:PATH`.
pub fn not_checked(doc: &mut Doc) -> Vec<String> {
    Dynamics.unchecked_claims(doc).unwrap_or_default().into_iter().map(|c| c.name).collect()
}

const EFFORT_CLAIM: &str =
    "lo <= effort <= hi, the effort the actuator supplies along its motion (inverse dynamics of the bodies it moves), as value = min(effort - lo, hi - effort) >= 0";

/// Why `not_checked` claims are not made.
pub const NOT_CHECKED: &str = "the file declares no (gravity gx gy gz), so its actuators' effort-limits claims are not made: add one (z up: (gravity 0 0 -9.81m/s2))";

impl Package for Dynamics {
    fn name(&self) -> &'static str {
        "dynamics"
    }

    fn forms(&self) -> &'static [&'static str] {
        &["gravity", "load", "command", "supports", "springs"]
    }

    fn commands(&self) -> &'static [ink_lang::command::Command] {
        &command::COMMANDS
    }

    fn clauses(&self) -> &'static [&'static str] {
        &["load", "supports"]
    }

    fn joint_clauses(&self) -> &'static [&'static str] {
        &["effort-limits", "actuator", "stiffness"]
    }

    fn relation_clauses(&self) -> &'static [&'static str] {
        &["couple"]
    }

    /// `(load NAME PATH ...)` in an assembly (clause.rs), `(effort-limits
    /// LO HI)` and `(stiffness K)` on a joint (joint.rs), `(couple ...)` in
    /// a relation (coupling.rs).
    fn clause(&self, r: &mut Reader, c: &Sexp, head: &str, args: &[Sexp]) -> R<Arc<dyn std::any::Any + Send + Sync>> {
        match head {
            "couple" => coupling::read(r, c, args),
            "effort-limits" | "actuator" | "stiffness" => joint::read(r, c, head, args),
            "supports" => {
                // In an assembly: its instance paths, made whole where the assembly is placed (finish).
                let [path, by @ ..] = args else { return Err(r.err_hint(c, "bad supports", Some(SUPPORTS_USAGE.into()))) };
                if by.is_empty() {
                    return Err(r.err_hint(c, "supports names the parts that hold it", Some(SUPPORTS_USAGE.into())));
                }
                let path = r.instance_path(path).map_err(|e| InkError { hint: Some(SUPPORTS_USAGE.into()), ..e })?;
                let by = by.iter().map(|s| r.instance_path(s)).collect::<R<Vec<String>>>().map_err(|e| InkError { hint: Some(SUPPORTS_USAGE.into()), ..e })?;
                Ok(Arc::new(AssemblySupports { path, by }))
            }
            _ => clause::read(r, c, args),
        }
    }

    /// What the design applies and how it says it is held, for `ink check`: its gravity, each load (the file's
    /// and its assemblies'), in the world's axes at the view's moment, and each reaction it divides among supports.
    fn check(&self, s: &mut ModelView) -> serde_json::Map<String, serde_json::Value> {
        use serde_json::json;
        let mut out = serde_json::Map::new();
        let at = |s: &ModelView, v: &[NodeId; 3]| json!(v.map(|n| s.doc.scalar(n, &s.ctx)));
        // The file's own: an imported file's top-level load is on that file's root, and is no load of this design.
        let mut loads: Vec<serde_json::Value> = (s.doc.forms::<Load>().filter(|l| l.run))
            .map(|l| json!({"name": l.name, "on": l.path, "force_N": at(s, &l.force), "torque_Nmm": at(s, &l.torque), "line": l.form.line}))
            .collect();
        let left = if s.doc.assembly("root").is_some() { loads_not_applied(&mut s.doc, "root") } else { Vec::new() };
        let placed = if s.doc.assembly("root").is_some() { clause::loads_under(&mut s.doc, "root") } else { Vec::new() };
        loads.extend(placed.iter().map(|l| json!({"name": l.name, "on": l.path, "force_N": at(s, &l.force), "torque_Nmm": at(s, &l.torque), "file": l.file, "line": l.line})));
        let supports: Vec<serde_json::Value> = divided(&s.doc).iter().map(|d| json!({"reaction": d.path, "by": d.supports})).collect();
        if loads.is_empty() && supports.is_empty() && left.is_empty() {
            return out;
        }
        if !left.is_empty() {
            out.insert(
                "loads_not_applied".into(),
                json!(left.iter().map(|l| json!({"name": l.name, "on": l.path, "line": l.line, "here": l.here, "why": l.said()})).collect::<Vec<_>>()),
            );
        }
        out.insert("gravity".into(), json!(gravity(&s.doc)));
        out.insert("loads".into(), json!(loads));
        out.insert("supports".into(), json!(supports));
        out
    }

    fn unchecked_claims(&self, doc: &mut Doc) -> R<Vec<ink_lang::UncheckedClaim>> {
        if gravity(doc).is_some() || doc.assembly("root").is_none() {
            return Ok(Vec::new());
        }
        let items = doc.items("root").map_err(InkError::msg)?;
        Ok(actuators(&items)
            .into_iter()
            .filter_map(|(path, joint)| Some((path, joint::effort(&joint).filter(|limit| limit.claimed)?)))
            .map(|(path, limit)| ink_lang::UncheckedClaim {
                name: format!("effort-limits:{path}"),
                file: limit.file,
                line: limit.line,
                reason: NOT_CHECKED,
                claim: EFFORT_CLAIM,
                declared_at: Some(limit.declared_at),
            })
            .collect())
    }

    fn values(&self) -> &'static [&'static str] {
        &["reaction", "tension", "effort-rms", "simulate", "following-error"]
    }

    /// `(reaction PATH [FRAME])`: six params this package provides at each
    /// moment, their frame kept for the provider made once the file is read.
    fn value(&self, r: &mut Reader, e: &Sexp, head: &str, args: &[Sexp]) -> R<Value> {
        if head == "simulate" {
            let (x0, v0) = match args {
                [x] => (x, None),
                [x, v] => (x, Some(v)),
                _ => return Err(r.err_hint(e, "bad simulate", Some(SIMULATE_USAGE.into()))),
            };
            let x0 = r.field(x0)?;
            let v0 = match v0 {
                Some(v) => r.field(v)?,
                None => r.graph().constant(0.0),
            };
            for (f, at) in [(x0, &args[0]), (v0, args.get(1).unwrap_or(e))] {
                if r.graph().depends_on(f, |x| matches!(x, Var::X | Var::Y | Var::Z | Var::T)) {
                    return Err(r.err_hint(at, "a simulation starts from a position and speed of params, not of x y z or t", Some(SIMULATE_USAGE.into())));
                }
            }
            let index = r.forms::<Simulated>().len();
            let run = run_of(r);
            let (owns, field) = taylor(r, &format!("simulated #{}", index + 1));
            // Plus its coordinate's own param (0 here): passed as an
            // argument, the def is that coordinate, so a solve or placement
            // reading it moves along it (its bodies weigh on it).
            let field = r.as_coordinate(field, &format!("coordinate of simulation #{}", index + 1));
            let speed = r.graph().diff(field, Var::T);
            r.keep("dynamics", "simulate", e, Arc::new(Simulated { field, speed, owns, line: e.line, x0, v0, run, index }))?;
            return Ok(Value::Field(field, None));
        }
        if head == "following-error" {
            let [p] = args else { return Err(r.err_hint(e, "bad following-error", Some(FOLLOWING_USAGE.into()))) };
            // A deviation from the written motion is a spring's: without (springs) there is none to read.
            if !springs_said(r) && r.in_file_run() {
                return Err(r.err_hint(
                    e,
                    "a following error is a joint's deviation from its written motion on its spring, and the file says no (springs): say (springs) first, \
                     so the joints its stiffnesses hold are simulated as springs for every answer, not only this one",
                    Some(SPRINGS_USAGE.into()),
                ));
            }
            let path = r.name(p)?;
            let n = r.forms::<Following>().into_iter().filter(|f| f.path == path).count();
            let again = if n == 0 { String::new() } else { format!("#{}", n + 1) };
            let run = run_of(r);
            let (owns, field) = taylor(r, &format!("following-error:{path}{again}"));
            r.keep("dynamics", "following-error", e, Arc::new(Following { path, owns, run }))?;
            return Ok(Value::Field(field, None));
        }
        if head == "effort-rms" {
            let [p] = args else { return Err(r.err_hint(e, "bad effort-rms", Some(RMS_USAGE.into()))) };
            let path = r.name(p)?;
            // A turn's effort is a torque, a slide's a force: the actuator's joint says which.
            let items = r.items("root").map_err(|m| r.err(e, m))?;
            let Some(turn) = actuators(&items).into_iter().find(|(a, _)| *a == path).map(|(_, j)| j.turn) else {
                return Err(r.err_hint(e, format!("{path} is no actuator under root (a joint with effort-limits)"), Some(RMS_USAGE.into())));
            };
            let n = r.forms::<Rms>().into_iter().filter(|x| x.path == path).count();
            let again = if n == 0 { String::new() } else { format!("#{}", n + 1) };
            let owns = r.derived(&format!("effort-rms:{path}{again}"));
            r.solved(owns);
            let run = r.in_file_run();
            r.keep("dynamics", "effort-rms", e, Arc::new(Rms { path, owns, run, at: e.clone() }))?;
            let f = r.graph().var(Var::Param(owns as u32));
            return Ok(r.worked(f, if turn { Quantity::Torque } else { Quantity::Force }, e));
        }
        if head == "tension" {
            let [a, b] = args else { return Err(r.err_hint(e, "bad tension", Some(TENSION_USAGE.into()))) };
            let (a, b) = (r.name(a)?, r.name(b)?);
            // Its lines, one param each.
            let related = r.relations_under("root").map_err(|m| r.err(e, m))?;
            let Some(lines) = related
                .iter()
                .filter_map(|x| Some((x, coupling::Couple::of(x.kind.as_deref()?)?)))
                .find(|(x, _)| (x.a == a && x.b == b) || (x.a == b && x.b == a))
                .map(|(_, c)| c.lines.len())
            else {
                return Err(InkError { hint: Some(TENSION_USAGE.into()), ..r.err(e, format!("no coupling relates {a} and {b}")) });
            };
            let n = r.forms::<Pull>().into_iter().filter(|p| p.a == a && p.b == b).count();
            let again = if n == 0 { String::new() } else { format!("#{}", n + 1) };
            let owns: Vec<usize> = (1..=lines).map(|i| r.derived(&format!("tension:{a}~{b}{again}.{i}"))).collect();
            owns.iter().for_each(|&k| r.solved(k));
            let run = r.in_file_run();
            r.keep("dynamics", "tension", e, Arc::new(Pull { a, b, owns: owns.clone(), run, at: e.clone() }))?;
            let lines: Vec<Value> = owns
                .iter()
                .map(|&k| {
                    let f = r.graph().var(Var::Param(k as u32));
                    r.worked(f, Quantity::Force, e)
                })
                .collect();
            return Ok(Value::list(lines));
        }
        let (p, frame) = match args {
            [p] => (p, None),
            [p, f] => (p, Some(f)),
            _ => return Err(r.err_hint(e, "bad reaction", Some(REACTION_USAGE.into()))),
        };
        let path = r.name(p)?;
        let mut by: Vec<(String, Option<Xf>)> = Vec::new();
        let frame = match frame {
            Some(f) if f.head_args().is_some_and(|(h, _)| h == "by") => {
                let (_, ss) = f.head_args().unwrap();
                for s in ss {
                    let (name, at) = match s.list() {
                        Some([n, x]) => match r.value(x)? {
                            Value::Frame(x) => (r.name(n)?, Some(x)),
                            other => return Err(r.err_hint(x, format!("a support's frame must be a frame, not a {}", other.kind()), Some(REACTION_USAGE.into()))),
                        },
                        _ => (r.name(s)?, None),
                    };
                    if by.iter().any(|b| b.0 == name) {
                        return Err(r.err(s, format!("{name} is named twice among the supports")));
                    }
                    by.push((name, at));
                }
                if by.is_empty() {
                    return Err(r.err_hint(f, "(by S ...) names the supports the reaction is split among", Some(REACTION_USAGE.into())));
                }
                None
            }
            Some(f) => match r.value(f)? {
                Value::Frame(x) => Some(x),
                other => return Err(r.err_hint(e, format!("a reaction's FRAME must be a frame, not a {}", other.kind()), Some(REACTION_USAGE.into()))),
            },
            None => None,
        };
        let same = |w: &&Want| w.path == path && w.by.len() == by.len() && w.by.iter().zip(&by).all(|(a, b)| a.0 == b.0);
        let n = r.forms::<Want>().into_iter().filter(|w| same(w)).count();
        let again = if n == 0 { String::new() } else { format!("#{}", n + 1) };
        let names: Vec<String> = if by.is_empty() { vec![format!("reaction:{path}{again}")] } else { by.iter().map(|(s, _)| format!("reaction:{path}{again}@{s}")).collect() };
        let mut owns = Vec::new();
        for name in &names {
            for c in ["fx", "fy", "fz", "mx", "my", "mz"] {
                let k = r.derived(&format!("{name}.{c}"));
                r.solved(k);
                owns.push(k);
            }
        }
        let run = r.in_file_run();
        let split = !by.is_empty();
        r.keep("dynamics", "reaction", e, Arc::new(Want { path, frame, by, owns: owns.clone(), run, at: e.clone(), named: head == "supports" }))?;
        // Three forces (N), then three moments (N mm).
        let six: Vec<Value> = owns
            .chunks(6)
            .map(|w| {
                let each: Vec<Value> = w
                    .iter()
                    .enumerate()
                    .map(|(c, &k)| {
                        let f = r.graph().var(Var::Param(k as u32));
                        r.worked(f, if c < 3 { Quantity::Force } else { Quantity::Torque }, e)
                    })
                    .collect();
                Value::list(each)
            })
            .collect();
        Ok(if split { Value::list(six) } else { six.into_iter().next().unwrap() })
    }

    fn read(&self, r: &mut Reader, form: &Sexp, head: &str, args: &[Sexp]) -> R<Arc<dyn Form>> {
        if head == "load" {
            return read_load(r, form, args);
        }
        if head == "supports" {
            // What holds a joint's load, said where the design is written: the same division of its reaction a
            // claim reads with (reaction PATH (by S ...)), wanted here by being said.
            let [path, by @ ..] = args else { return Err(r.err_hint(form, "bad supports", Some(SUPPORTS_USAGE.into()))) };
            if by.is_empty() {
                return Err(r.err_hint(form, "supports names the parts that hold it", Some(SUPPORTS_USAGE.into())));
            }
            let name = r.name(path)?;
            self.value(r, form, "supports", &[path.clone(), listed(form, "by", by)])?;
            return Ok(Arc::new(Supports { path: name }));
        }
        if head == "springs" {
            // The file's stated stiffnesses simulated as springs: what a joint they hold supplies, carries and
            // lags is its spring's, whatever the file reads of it.
            if !args.is_empty() {
                return Err(r.err_hint(form, "springs takes nothing", Some(SPRINGS_USAGE.into())));
            }
            let in_run = r.in_file_run();
            if in_run && r.forms::<Springs>().iter().any(|s| s.run) {
                return Err(r.err(form, "the file says (springs) twice"));
            }
            return Ok(Arc::new(Springs { run: in_run, sim: run_of(r) }));
        }
        if head == "command" {
            let [path, effort] = args else { return Err(r.err_hint(form, "bad command", Some(COMMAND_USAGE.into()))) };
            let path = r.name(path)?;
            let effort = r.field(effort)?;
            if r.graph().depends_on(effort, |x| matches!(x, Var::X | Var::Y | Var::Z)) {
                return Err(r.err_hint(form, "a command is a field of params, t and the simulated coordinates, not of x y z", Some(COMMAND_USAGE.into())));
            }
            return Ok(Arc::new(Command { path, effort, run: r.in_file_run(), form: form.clone() }));
        }
        let [x, y, z] = args else { return Err(r.err_hint(form, "bad gravity", Some(USAGE.into()))) };
        let mut n = |e| r.number_of(e, Quantity::Acceleration, "gravity");
        let g = [n(x)?, n(y)?, n(z)?];
        if !g.iter().all(|v| v.is_finite()) {
            return Err(r.err_hint(form, "gravity must be finite", Some(USAGE.into())));
        }
        let run = r.in_file_run();
        if run && r.forms::<Gravity>().iter().any(|g| g.run) {
            return Err(r.err(form, "the file sets (gravity ...) twice"));
        }
        Ok(Arc::new(Gravity { g, run, at: form.clone() }))
    }

    /// Each actuator under the root assembly, in a file with gravity: its
    /// effort, provided at each moment, claimed within its limits; each
    /// reaction asked for, provided at each moment. Loads and reactions
    /// name paths placed under the root.
    fn finish(&self, r: &mut Reader) -> R<()> {
        // A reaction or tension an import asks for is answered in its own
        // file, with its gravity, loads and span: here its claims are
        // refused, saying so.
        let mut imported: Vec<(Vec<usize>, Sexp)> = r.forms::<Want>().into_iter().filter(|w| !w.run).map(|w| (w.owns.clone(), w.at.clone())).collect();
        imported.extend(r.forms::<Pull>().into_iter().filter(|p| !p.run).map(|p| (p.owns.clone(), p.at.clone())));
        imported.extend(r.forms::<Rms>().into_iter().filter(|x| !x.run).map(|x| (vec![x.owns], x.at.clone())));
        for (owns, at) in imported {
            r.keep("dynamics", "reaction", &at, Arc::new(Elsewhere { owns }))?;
        }
        // The supports an assembly names, where the root places it: each a division of its joint's reaction.
        {
            let placed = r.clauses_under("root", "dynamics").unwrap_or_default();
            for p in placed {
                let Some(named) = p.clause.kept.downcast_ref::<AssemblySupports>() else { continue };
                let at = Sexp { k: ink_lang::sexp::SK::List(Vec::new().into()), line: p.clause.line, col: 1 };
                let sym = |path: &str| Sexp { k: ink_lang::sexp::SK::Sym(format!("{}/{path}", p.path)), line: at.line, col: 1 };
                let by: Vec<Sexp> = named.by.iter().map(|s| sym(s)).collect();
                self.value(r, &at, "supports", &[sym(&named.path), listed(&at, "by", &by)])?;
            }
        }
        let gravity = r.forms::<Gravity>().into_iter().find(|g| g.run).map(|g| (g.g, g.at.clone()));
        type Asked = (String, Option<Xf>, Vec<usize>, Sexp, Vec<(String, Option<Xf>)>, bool);
        let wants: Vec<Asked> = r.forms::<Want>().into_iter().filter(|w| w.run).map(|w| (w.path.clone(), w.frame, w.owns.clone(), w.at.clone(), w.by.clone(), w.named)).collect();
        let loads: Vec<(String, Sexp)> = r.forms::<Load>().into_iter().filter(|l| l.run).map(|l| (l.path.clone(), l.form.clone())).collect();
        let pulls: Vec<(String, String, Vec<usize>, Sexp)> =
            r.forms::<Pull>().into_iter().filter(|p| p.run).map(|p| (p.a.clone(), p.b.clone(), p.owns.clone(), p.at.clone())).collect();
        let rms: Vec<(String, usize, Sexp)> = r.forms::<Rms>().into_iter().filter(|x| x.run).map(|x| (x.path.clone(), x.owns, x.at.clone())).collect();
        // A simulation's provider reads what its model reads, by name.
        let names = r.param_names();
        if !r.forms::<Simulated>().is_empty() || !r.forms::<Following>().is_empty() || springs_said(r) {
            let _ = run_of(r).names.set(names);
        }
        // A command is an actuator's, moving a simulated coordinate.
        let commands: Vec<(String, NodeId, Sexp)> = r.forms::<Command>().into_iter().filter(|c| c.run).map(|c| (c.path.clone(), c.effort, c.form.clone())).collect();
        if !r.forms::<Simulated>().is_empty() || !r.forms::<Following>().is_empty() || springs_said(r) {
            let run = run_of(r);
            let _ = r.items("root");
            let _ = r.relations_under("root");
            let _ = clause::loads_read(r, "root");
            let mut inputs = r.motion_inputs();
            // A simulated coordinate uses rates of the prescribed motion.
            // A springs-only response also accepts a position step: the
            // integrator transfers the jump to the following error while
            // preserving the joint's physical position and speed. Explicit
            // derivatives in those inputs still belong to motion_inputs.
            if !r.forms::<Simulated>().is_empty() {
                inputs.extend(r.coordinate_rates("root"));
            }
            inputs.sort_unstable();
            inputs.dedup();
            let _ = run.inputs.set(inputs);
        }
        if !commands.is_empty() {
            // What a dynamics answer provides (a reaction, a tension, an RMS
            // effort) is no state a simulation steps through.
            let mut answers: Vec<u32> = r.forms::<Want>().into_iter().flat_map(|w| w.owns.clone()).map(|k| k as u32).collect();
            answers.extend(r.forms::<Pull>().into_iter().flat_map(|p| p.owns.clone()).map(|k| k as u32));
            answers.extend(r.forms::<Rms>().into_iter().map(|x| x.owns as u32));
            answers.extend(r.forms::<Following>().into_iter().flat_map(|f| f.owns).map(|k| k as u32));
            for (path, effort, at) in &commands {
                if r.graph().depends_on(*effort, |v| matches!(v, Var::Param(k) if answers.contains(&k))) {
                    return Err(InkError {
                        hint: Some(COMMAND_USAGE.into()),
                        ..r.err(
                            at,
                            format!("the command on {path} reads a reaction, a tension, an RMS effort or a following error, which a simulation does not provide at its states"),
                        )
                    });
                }
            }
            let simulated: Vec<u32> = r.forms::<Simulated>().iter().flat_map(|s| s.owns).map(|k| k as u32).collect();
            let items = r.items("root").unwrap_or_default();
            let acts = actuators(&items);
            for (path, _, at) in &commands {
                let Some((_, j)) = acts.iter().find(|(p, _)| p == path) else {
                    let hint = ink_lang::closest(path, acts.iter().map(|(p, _)| p.as_str())).unwrap_or_else(|| COMMAND_USAGE.into());
                    return Err(InkError { hint: Some(hint), ..r.err(at, format!("{path} is no actuator placed under root (a joint with effort-limits)")) });
                };
                if !r.graph().depends_on(j.drive, |x| matches!(x, Var::Param(k) if simulated.contains(&k))) {
                    return Err(InkError {
                        hint: Some(COMMAND_USAGE.into()),
                        ..r.err(at, format!("the command on {path}: its drive moves no simulated coordinate, so no simulation reads it"))
                    });
                }
            }
            if commands.iter().enumerate().any(|(i, c)| commands[..i].iter().any(|d| d.0 == c.0)) {
                return Err(r.err(&commands[0].2, "an actuator is commanded twice"));
            }
        }
        if gravity.is_none() && wants.is_empty() && loads.is_empty() && pulls.is_empty() && rms.is_empty() {
            return Ok(());
        }
        // What the model reads cannot read a reaction, which sums it: a
        // load (friction as mu times a normal force), a placement, a drive,
        // a limit, a density.
        let mut provided: Vec<u32> = r.forms::<Want>().into_iter().flat_map(|w| w.owns.iter().copied()).map(|k| k as u32).collect();
        provided.extend(r.forms::<Pull>().into_iter().flat_map(|p| p.owns.iter().copied()).map(|k| k as u32));
        if !provided.is_empty() {
            let reads = |r: &mut Reader, ns: &[NodeId]| ns.iter().any(|&n| r.graph().depends_on(n, |v| matches!(v, Var::Param(k) if provided.contains(&k))));
            let why = "a reaction sums the loads and bodies it carries, so none of them can read one: whether a support's friction holds is std's \
                       (slip R) against its coefficient, (value (slip (reaction root (at))) < 0.8)";
            let loaded: Vec<(String, Vec<NodeId>, Sexp)> = r
                .forms::<Load>()
                .into_iter()
                .filter(|l| l.run)
                .map(|l| {
                    (
                        l.name.clone(),
                        l.force.iter().chain(&l.torque).copied().chain(l.at.iter().flat_map(|x| x.r.iter().flatten().chain(&x.t).copied().collect::<Vec<_>>())).collect(),
                        l.form.clone(),
                    )
                })
                .collect();
            for (name, ns, at) in loaded {
                if reads(r, &ns) {
                    return Err(InkError { hint: Some(why.into()), ..r.err(&at, format!("the load {name} reads a reaction")) });
                }
            }
            let at = r
                .forms::<Want>()
                .into_iter()
                .next()
                .map(|w| w.at.clone())
                .or_else(|| r.forms::<Pull>().into_iter().next().map(|p| p.at.clone()))
                .expect("a reaction or tension asked for");
            for l in clause::loads_read(r, "root").map_err(|m| r.err(&at, m))? {
                let ns: Vec<NodeId> =
                    l.force.iter().chain(&l.torque).copied().chain(l.at.iter().flat_map(|x| x.r.iter().flatten().chain(&x.t).copied().collect::<Vec<_>>())).collect();
                if reads(r, &ns) {
                    return Err(InkError { hint: Some(why.into()), ..r.err(&at, format!("the load {} ({}:{}) reads a reaction", l.name, l.file, l.line)) });
                }
            }
            if let Ok(items) = r.items("root") {
                for it in &items {
                    let mut ns: Vec<NodeId> = vec![it.shape];
                    ns.extend(it.instance.iter().flat_map(|(_, x)| x.r.iter().flatten().chain(&x.t).copied().collect::<Vec<_>>()));
                    ns.extend(it.attrs.iter().filter_map(|(_, v)| if let Value::Field(f, _) = v { Some(*f) } else { None }));
                    for j in it.links.iter().filter_map(|l| l.joint.as_ref()) {
                        ns.push(j.drive);
                        ns.extend([j.limits, j.speed_limits, joint::effort(j).map(|e| (e.lo, e.hi))].into_iter().flatten().flat_map(|(a, b)| [a, b]));
                    }
                    if reads(r, &ns) {
                        let at = r
                            .forms::<Want>()
                            .into_iter()
                            .next()
                            .map(|w| w.at.clone())
                            .or_else(|| r.forms::<Pull>().into_iter().next().map(|p| p.at.clone()))
                            .expect("a reaction or tension asked for");
                        return Err(InkError {
                            hint: Some(why.into()),
                            ..r.err(&at, format!("{} reads a reaction (its shape, placement, motion, limits or attributes)", it.label))
                        });
                    }
                }
            }
        }
        let Ok(items) = r.items("root") else {
            return match wants.first().map(|w| &w.3).or(loads.first().map(|l| &l.1)).or(pulls.first().map(|p| &p.3)) {
                Some(at) => Err(r.err(at, "loads and reactions are of what the root assembly places, and the file has none")),
                None => Ok(()),
            };
        };
        // A reaction may be root's (the support's); a load and a support
        // are of an instance.
        let reacting = wants.iter().map(|w| (&w.0, &w.3, false, REACTION_USAGE));
        let supports = wants.iter().flat_map(|w| w.4.iter().map(move |s| (&s.0, &w.3, true, REACTION_USAGE)));
        let ends = pulls.iter().flat_map(|p| [(&p.0, &p.3, true, TENSION_USAGE), (&p.1, &p.3, true, TENSION_USAGE)]);
        for (path, at, instance, usage) in reacting.chain(loads.iter().map(|l| (&l.0, &l.1, true, LOAD_USAGE))).chain(supports).chain(ends) {
            if (instance || path != "root") && !placed(&items, path) {
                let names = items.iter().flat_map(|it| &it.links).map(|l| l.path.as_str());
                let hint = ink_lang::closest(path, names).unwrap_or_else(|| usage.into());
                return Err(InkError { hint: Some(hint), ..r.err(at, format!("no instance {path} is placed under root")) });
            }
        }
        let Some((g, at)) = gravity else {
            let at = &wants.first().map(|w| &w.3).or(loads.first().map(|l| &l.1)).or(pulls.first().map(|p| &p.3)).unwrap_or_else(|| &rms[0].2);
            return Err(InkError {
                hint: Some(USAGE.into()),
                ..r.err(at, "an effort, a reaction or a tension needs the world's gravity: the file declares no (gravity gx gy gz)")
            });
        };
        let acts = actuators(&items);
        if acts.is_empty() && wants.is_empty() && pulls.is_empty() && rms.is_empty() {
            return Ok(());
        }
        let reload = r.reloader().ok_or_else(|| r.err(&at, "effort claims need the file's source"))?;
        let mut owns = Vec::new();
        for (path, j) in &acts {
            let k = r.derived(&format!("effort:{path}"));
            r.solved(k);
            owns.push(k);
            let lim = joint::effort(j).unwrap();
            // An actuator that states no limits claims none: its effort is worked out where something reads it.
            if !lim.claimed {
                continue;
            }
            let rate = r.derivative(j.drive, Var::T);
            let gr = r.graph();
            let e = gr.var(Var::Param(k as u32));
            let (lo, hi) = (gr.sub(e, lim.lo), gr.sub(lim.hi, e));
            let claim = Claim::Value(gr.min(lo, hi), Cmp::Ge, 0.0);
            let quantity = Some(if j.turn { Quantity::Torque } else { Quantity::Force });
            let terms = vec![("effort".into(), e), ("rate".into(), rate), ("lo".into(), lim.lo), ("hi".into(), lim.hi)];
            r.assert_from(
                &at,
                AssertDef {
                    name: format!("effort-limits:{path}"),
                    claim,
                    text: EFFORT_CLAIM.into(),
                    bound: None,
                    res: None,
                    accuracy: lim.accuracy,
                    quantity,
                    file: lim.file,
                    // The run supplies gravity, loads and motion for every placed actuator.
                    scope: r.scope(),
                    line: lim.line,
                    terms,
                    fixed: Vec::new(),
                },
                Some(lim.declared_at),
            )?;
        }
        // One provider per actuator and per reaction, sharing the model:
        // each claim is decided, or refused, on its own.
        let names = r.param_names();
        let mut paths: Vec<String> = Vec::new();
        for (p, ..) in &wants {
            if !paths.contains(p) {
                paths.push(p.clone());
            }
        }
        // A file simulating its springs (it says so) answers what they touch from them; another
        // simulation's refusal touches none of these answers.
        let run = springs_said(r).then(|| run_of(r));
        let held = wants.iter().filter(|w| !w.4.is_empty()).map(|w| (w.0.clone(), w.4.iter().map(|(s, _)| s.clone()).collect())).collect();
        let shared = Arc::new(Shared { reload, gravity: g, paths: acts.into_iter().map(|(p, _)| p).collect(), reactions: paths.clone(), held, names, model: OnceLock::new(), run });
        for (a, k) in owns.into_iter().enumerate() {
            r.keep("dynamics", "effort-limits", &at, Arc::new(Effort { shared: shared.clone(), actuator: a, owns: [k] }))?;
        }
        for (path, k, at) in rms {
            let Some(actuator) = shared.paths.iter().position(|p| *p == path) else {
                let hint = ink_lang::closest(&path, shared.paths.iter().map(String::as_str)).unwrap_or_else(|| RMS_USAGE.into());
                return Err(InkError { hint: Some(hint), ..r.err(&at, format!("{path} is no actuator placed under root (a joint with effort-limits)")) });
            };
            r.keep("dynamics", "effort-rms", &at, Arc::new(Rmsing { shared: shared.clone(), actuator, owns: [k] }))?;
        }
        for (claim, (path, frame, owns, at, by, named)) in wants.into_iter().enumerate() {
            // Its frame, if it names one, and its supports': compiled here,
            // with the params they read.
            let n = shared.names.len();
            let frame = frame.map(|x| compiled(r, &x, n));
            let mut supports = Vec::new();
            // A part on a joint is held by that joint: its load is the joint's reaction, and it is no support
            // of itself. (`(reaction P (by P))` alone, a claim's, reads that reaction in P's load frame.)
            let own =
                |s: &str| s == path || s.strip_prefix(path.as_str()).is_some_and(|rest| rest.starts_with('/')) || path.strip_prefix(s).is_some_and(|rest| rest.starts_with('/'));
            if let Some((s, _)) = by.iter().find(|(s, _)| own(s)).filter(|_| named && by.iter().all(|(s, _)| own(s))) {
                return Err(InkError {
                    hint: Some(SUPPORTS_USAGE.into()),
                    ..r.err(
                        &at,
                        format!(
                            "{s} cannot support {path}: it is the part on that joint (or part of it), and a part on a joint is held by the joint. Its load is the \
                             joint's reaction, which a requirement of its ratings reads with no supports line: delete this one. supports names other parts that hold \
                             what turns or slides there (the bearings a shaft turns in)"
                        ),
                    )
                });
            }
            for (s, given) in by {
                let (carries, pulls, rating) = carries(r, &items, &s, &at)?;
                let load = load_frame(r, &items, &s, &at)?;
                supports.push(Support { path: s, carries, pulls, rating, load: compiled(r, &load, n), given: given.map(|x| compiled(r, &x, n)) });
            }
            let drive = (!supports.is_empty())
                .then(|| items.iter().flat_map(|it| &it.links).find(|l| l.path == path).and_then(|l| l.joint.as_ref()))
                .flatten()
                .map(|j| split::Drive { turn: j.turn, axis: j.axis });
            let index = paths.iter().position(|p| *p == path).unwrap();
            r.keep("dynamics", "reaction", &at, Arc::new(Reacting { shared: shared.clone(), index, claim, frame, owns, supports, drive, named }))?;
        }
        for (a, b, owns, at) in pulls {
            r.keep("dynamics", "tension", &at, Arc::new(Tensing { shared: shared.clone(), a, b, owns }))?;
        }
        Ok(())
    }
}

/// A reaction or tension an imported file asks for: its
/// params, never set here.
#[derive(Clone)]
struct Elsewhere {
    owns: Vec<usize>,
}

impl Form for Elsewhere {
    fn name(&self) -> &str {
        "reaction"
    }

    fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(self.clone()))
    }
}

impl ParamProvider for Elsewhere {
    fn name(&self) -> &str {
        "a reaction or tension an imported file asks for"
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

    fn provide(&self, _: &mut [f64], _: f64) -> Result<Option<String>, Refused> {
        Err("a reaction or a tension is answered in the file that asks for it, with that file's (gravity ...), loads and span: check the imported file itself".into())
    }
}
