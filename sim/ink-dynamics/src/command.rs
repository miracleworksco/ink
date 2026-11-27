//! `ink effort`: the effort each actuator supplies along the written motion,
//! at a moment or over the span.

use crate::{self as dynamics, Model, Moment};
use ink_lang::InkError;
use ink_lang::command::{Args, Command, Opt, bounds, n, t_span};
use ink_lang::view::ModelView;
use serde_json::{Value as J, json};

/// The package's commands.
pub(crate) static COMMANDS: [Command; 1] = [Command {
    name: "effort",
    options: &[Opt::reads("target"), Opt::reads("t-span"), Opt::reads("least-time"), Opt::reads("modes")],
    help: "The effort each actuator (a joint with effort-limits) supplies along the written motion, by inverse
                dynamics, at --t or over --t-span (default the file's): effort, rate, limits, margin and reflected
                inertia (over a span: peak, RMS, least margin); the effort along each coordinate (an argument the
                target's instances pass that a joint's drive reads); what no actuator supplies (unpowered).
                Needs (gravity gx gy gz). An estimate. [--target ASSEMBLY] (default root)
                [--modes] (natural frequencies at the pose at --t, from stated stiffness; needs no gravity)
                example: ink effort examples/arm.ink",
    tool: r#"The effort each actuator (a joint with effort-limits) supplies along the written motion, by inverse dynamics, at t or over a span: effort, rate, limits, margin, reflected inertia (over a span: peak, RMS, least margin), effort along each coordinate, what no actuator supplies (unpowered). A rate needs nothing more; an effort needs (gravity gx gy gz) and its bodies' masses, else `lacks` says why. modes: the natural frequencies at t's pose, lowest first, each with the joints moving in it, from stated stiffness (a belt's EA, a stepper's holding stiffness). An estimate. {"target": "root", "t-span": "0:2"} or {"modes": true, "t": 0}"#,
    run: ink_lang::command::Run::Model(cmd_effort),
}];

fn cmd_effort(s: &mut ModelView, a: &Args) -> Result<J, InkError> {
    // An option this release holds back, read so it is refused by name, not as one ink effort lacks.
    if a.get("least-time").is_some() {
        return Err(InkError::msg(
            "--least-time is not in this release of ink: the least time to run a path within its actuators' limits is held back until it is ready to rely on",
        ));
    }
    let target = a.get("target").unwrap_or("root");
    let reload = s.doc.reloader().ok_or_else(|| InkError::msg("ink effort needs the file's source"))?;
    if a.flag("modes")? {
        return modes(s, &reload, target);
    }
    let gravity = dynamics::gravity(&s.doc);
    let span = t_span(a)?.or(s.doc.t_span.filter(|_| a.get("t").is_none()));
    // Each number from the least it needs, as the claims' plug gives them (ink-effort): the motion (a rate,
    // a coordinate's value) whatever the masses and the gravity; an effort, a reaction or a tension where
    // the moving bodies it rests on have a mass and a gravity is declared, and otherwise null, saying why.
    let Read { model, steps, at, gaps, method, weightless, .. } = read_kinematic(s, target, span)?;
    let mut out = json!({"target": target, "gravity": gravity.map(|g| g.map(n))});
    let moments = at.clone();
    let (mut actuators, mut coordinates, mut unpowered, mut reactions, mut couplings) = match (span, steps) {
        (Some(span), Some(mut steps)) => {
            // An actuator's rating gaps are loose steps of the answer's.
            steps.loose.extend(gaps.into_iter().flatten());
            steps.loose.sort_by(|a, b| a.0.total_cmp(&b.0));
            let spans = |v: &[(f64, f64)]| J::Array(joined(v).into_iter().map(|(a, b)| bounds(a, b)).collect());
            // Where a step is unchecked the peaks, margins and RMS below are
            // the rest's, and where one is loose they read it coarser than
            // the rule: either way partial.
            let middles = at.iter().filter(|(width, _)| *width > 0.0).count();
            out["motion"] = json!({"span": bounds(span.0, span.1), "steps": middles, "basis": steps.basis, "unchecked": spans(&steps.unchecked),
                "loose": spans(&steps.loose), "partial": !(steps.unchecked.is_empty() && steps.loose.is_empty())});
            over_span(&model, span, &at)
        }
        _ => {
            let m = &at[0].1;
            out["t"] = n(m.t);
            at_moment(&model, m)
        }
    };
    withhold(&model, weightless, &moments, [&mut actuators, &mut coordinates, &mut unpowered, &mut reactions, &mut couplings]);
    if !model.massless.is_empty() {
        out["without_mass"] = json!(model.massless.iter().map(|(label, ..)| label).collect::<Vec<_>>());
    }
    out["actuators"] = J::Array(actuators);
    out["coordinates"] = J::Array(coordinates);
    out["reactions"] = J::Array(reactions);
    if !couplings.is_empty() {
        out["couplings"] = J::Array(couplings);
    }
    if !unpowered.is_empty() {
        out["unpowered"] = J::Array(unpowered);
    }
    out["driven_by_nothing"] = json!(model.driven_by_nothing);
    out["carried"] = J::Array(model.carried.iter().map(|(t, m)| json!({"target": t, "mass": n(*m)})).collect());
    out["deforming"] = J::Array(model.deforming.iter().map(|(t, m)| json!({"target": t, "mass": m.map(n)})).collect());
    out["method"] = json!(method);
    out["accuracy"] = J::Null;
    // A motion not written but simulated says so: how far it was stepped
    // and why it stopped short, if it did, and how each joint its springs
    // hold followed what the written motion commands.
    let simulation = dynamics::simulation(&s.doc);
    let sprung = dynamics::with_springs(simulation.clone()).and_then(Result::ok);
    supplied(&mut out["actuators"], &model, sprung.as_deref(), span.ok_or(s.ctx.t));
    if !model.simulated.is_empty() || simulation.is_some() {
        let names: Vec<&str> = model.simulated.iter().map(|&(k, _)| model.coords[k].name.as_str()).collect();
        out["simulated"] = json!({"coordinates": names});
        // The efforts above stand however the simulation goes: where it
        // does not start, it says why here.
        match simulation {
            Some(Ok(sim)) => {
                let r = sim.stepped();
                let stopped = r.stopped.as_ref().map_or(J::Null, |(t, why)| json!({"t": n(*t), "why": why}));
                for (key, v) in [("end", n(r.end)), ("stopped", stopped), ("steps", json!(r.steps)), ("rejected", json!(r.rejected)), ("evaluations", json!(r.evaluations))] {
                    out["simulated"][key] = v;
                }
                out["simulated"]["method"] = json!(sim.note());
                out["simulated"]["accuracy"] = J::Null;
                let r = sim.springs();
                if !r.joints.is_empty() || !r.held.is_empty() {
                    out["simulated"]["springs"] = springs(&r);
                }
            }
            Some(Err(why)) => out["simulated"]["refused"] = json!(why),
            None => {}
        }
    }
    let res = if model.measured > 0 { json!([n(model.res.0), n(model.res.1)]) } else { J::Null };
    out["masses"] = json!({"bodies": model.measured, "published": model.published, "spin_inertia": model.spin_inertia, "res": res, "unsettled": model.unsettled.iter().chain(&model.unsettled_carried).collect::<Vec<_>>()});
    Ok(out)
}

/// Takes out of `ink effort`'s rows (actuators, coordinates, unpowered, reactions, couplings, as `at_moment`
/// or `over_span` made them) every number the model does not have: with no gravity declared every effort,
/// reaction and tension; with moving bodies left out for having no mass (`Model::massless`), those of the
/// actuators, coordinates, reactions and couplings one of them would change. Each row says why (`lacks`),
/// and keeps its motion: an actuator's rate, at the moment or its fastest over the span.
fn withhold(model: &Model, weightless: bool, at: &[(f64, std::sync::Arc<Moment>)], rows: [&mut Vec<J>; 5]) {
    let why = |without: Vec<&str>| match (weightless, without.is_empty()) {
        (true, _) => Some(dynamics::no_gravity().message),
        (false, false) => Some(format!("these moving bodies it rests on have no mass: {}", model.massless_said(&without))),
        (false, true) => None,
    };
    let [actuators, coordinates, unpowered, reactions, couplings] = rows;
    let strip = |j: &mut J, keep: &[&str], why: String| {
        let kept: serde_json::Map<String, J> = j.as_object().into_iter().flatten().filter(|(k, _)| keep.contains(&k.as_str())).map(|(k, v)| (k.clone(), v.clone())).collect();
        *j = J::Object(kept);
        j["lacks"] = json!(why);
    };
    for (k, j) in actuators.iter_mut().enumerate() {
        let Some(why) = why(model.blocked(k)) else { continue };
        strip(j, &["path", "rate"], why);
        if at.len() > 1
            && let Some((_, m)) = peak(at.iter(), |x| x.1.rate[k].abs())
        {
            j["peak"] = json!({"rate": n(m.rate[k]), "t": n(m.t)});
        }
    }
    let blocked: Vec<Option<String>> = (0..model.coords.len()).map(|k| why(model.blocks_coord(k))).collect();
    for (k, j) in coordinates.iter_mut().enumerate() {
        let Some(why) = blocked[k].clone() else { continue };
        strip(j, &["name", "value", "rate", "acceleration"], why);
    }
    unpowered.retain(|j| !model.coords.iter().zip(&blocked).any(|(c, why)| why.is_some() && j["coordinate"] == json!(c.name)));
    for (r, j) in reactions.iter_mut().enumerate() {
        let Some(why) = why(model.blocks_reaction(r)) else { continue };
        strip(j, &["path"], why);
    }
    for (c, j) in couplings.iter_mut().enumerate() {
        let Some(why) = why(model.blocks_coupling(c)) else { continue };
        strip(j, &["relation"], why);
    }
}

/// What `ink effort FILE --target TARGET` reads before it prints: the rigid
/// model, and its moments, each with its step's width. Over a span (the
/// file's, or the one asked): the span's two ends (no share of its time)
/// and each step's middle, every step but an unchecked one (nothing is read
/// there; a loose one is read at its middle, coarser than the rule), with
/// the steps as chosen and each actuator's own loose steps besides (where
/// its rating may cross its effort inside one). With no span: the one
/// moment the view is at.
pub struct Read {
    pub model: Model,
    pub span: Option<(f64, f64)>,
    pub steps: Option<ink_lang::Steps>,
    pub at: Vec<(f64, std::sync::Arc<Moment>)>,
    pub gaps: Vec<Vec<(f64, f64)>>,
    /// Each actuator: whether the file's simulation releases its joint for
    /// its springs (its effort is then its spring's, not the rigid model's:
    /// `ink effort`'s basis `spring`), or why it is refused there.
    pub sprung: Vec<Result<bool, String>>,
    /// Read with no gravity declared (`read_kinematic`): its motion is the motion's, and none of its
    /// efforts or reactions is one.
    pub weightless: bool,
    /// How the efforts are found, in words (`ink effort`'s `method`).
    pub method: &'static str,
}

/// `ink effort`'s reading of `target` (an assembly or instance path) of a
/// loaded design, at the view's moment or over `span`: the command's own
/// entry, for a caller that wants the moments and not its words. Err: no
/// gravity, a model that cannot be built, or with no span a moment the
/// solves do not solve at.
pub fn read(s: &mut ModelView, target: &str, span: Option<(f64, f64)>) -> Result<Read, InkError> {
    read_as(s, target, span, false)
}

/// `read`, of the motion whatever the masses and the gravity (`build_kinematic`): the moments and each
/// coordinate's value and rate at them are the motion's own; with no gravity declared none is taken (the
/// model's is zero, and `Read::weightless` says so: no effort or reaction of it is one), and the moving bodies
/// with no mass are left out and named (`Model::massless`: what each would change is refused).
pub fn read_kinematic(s: &mut ModelView, target: &str, span: Option<(f64, f64)>) -> Result<Read, InkError> {
    read_as(s, target, span, true)
}

fn read_as(s: &mut ModelView, target: &str, span: Option<(f64, f64)>, kinematic: bool) -> Result<Read, InkError> {
    let reload = s.doc.reloader().ok_or_else(|| InkError::msg("ink effort needs the file's source"))?;
    let declared = dynamics::gravity(&s.doc);
    let gravity = match (declared, kinematic) {
        (Some(g), _) => g,
        (None, true) => [0.0; 3],
        (None, false) => return Err(dynamics::no_gravity()),
    };
    let held: Vec<(String, Vec<String>)> = dynamics::divided(&s.doc).into_iter().map(|d| (d.path, d.supports)).collect();
    dynamics::supports_near(s, target, &held)?;
    let paths = dynamics::joints(&mut s.doc, target)?;
    let model = match kinematic {
        true => dynamics::build_kinematic(&reload, target, gravity, span.map_or(s.ctx.t, |s| s.0), &paths)?,
        false => dynamics::build(&reload, target, gravity, span.map_or(s.ctx.t, |s| s.0), &paths, false)?,
    };
    let (steps, at, gaps) = match span {
        None => {
            let m = model.at(s.ctx.t);
            if let Some(why) = &m.unsolved {
                return Err(InkError::msg(why.clone()));
            }
            (None, vec![(0.0, m)], Vec::new())
        }
        Some(span) => {
            let steps = model.steps(span);
            let gaps: Vec<Vec<(f64, f64)>> = (0..model.actuators.len()).map(|a| model.rating_gaps(a, span, &steps)).collect();
            let skipped = |w: &[f64]| steps.unchecked.contains(&(w[0], w[1]));
            let middles = steps.ends.windows(2).filter(|w| !skipped(w)).map(|w| (w[1] - w[0], ink_core::eval::span_mid((w[0], w[1]))));
            // An end where the solves do not solve is unchecked, as a step.
            let ends = |t: f64| (!steps.unchecked.contains(&(t, t))).then_some((0.0, t));
            let at = ends(span.0).into_iter().chain(middles).chain(ends(span.1)).map(|(w, t)| (w, model.over(span, t))).collect();
            (Some(steps), at, gaps)
        }
    };
    let sprung = dynamics::with_springs(dynamics::simulation(&s.doc)).and_then(Result::ok);
    let sprung = (0..model.actuators.len())
        .map(|a| match &sprung {
            None => Ok(false),
            Some(sim) => dynamics::Sprung::of(sim, model.commanded(a), &model.actuators[a].path).map(|s| matches!(s, dynamics::Sprung::Spring { .. })),
        })
        .collect();
    Ok(Read { model, span, steps, at, gaps, sprung, weightless: declared.is_none(), method: if span.is_some() { EFFORT_SPAN } else { EFFORT_MOMENT } })
}

/// How each answer is found, in its JSON (`method`): every one an estimate,
/// `accuracy` null.
const EFFORT_MOMENT: &str = "inverse dynamics of the rigid mechanism at moment t, from the written motion's position, speed and acceleration there \
    and the masses measured (masses): an estimate";
const EFFORT_SPAN: &str = "inverse dynamics of the rigid mechanism at the span's ends and the middle of each step (motion.basis): each peak, margin \
    and RMS is of those moments, an estimate between them; a `bound` is how far bodies read round may move it, at any angle of their spins";
const MODES: &str = "the mechanism linearised about its pose at t: its generalized inertia there on the stiffnesses stated (springs), each \
    joint they hold given a coordinate of its own (released), the rest rigid; gravity's and tensions' stiffness, damping and belts' mass are \
    left out; the bodies' masses, centres and inertias read as measured, exactly (no band for what their measures are known to): an \
    estimate";

const SPRING: &str = "its spring's force against its following error in the file's simulation (simulated.method), as its claims read it, its \
    limits at the commanded speed; over a span at the efforts' steps and where the spring's force turns back (`loose` where a limit may jump \
    across it inside a step), the RMS by Gauss-Legendre over the simulation's steps; `band` its stiffness times the following error's band at \
    the span's end, measured (band_kind), not a bound: an estimate";

/// The `springs` of `ink effort`'s `simulated`: each released joint's
/// following error, the joints held, the lines' tensions.
fn springs(r: &dynamics::simulate::Springs) -> J {
    let joints: Vec<J> = r
        .joints
        .iter()
        .map(|f| {
            let ringing = f.ringing.as_ref().map_or(J::Null, |g| json!({"from": n(g.from), "amplitude": n(g.amplitude), "frequency": g.frequency.map(n), "share": g.share.map(n)}));
            json!({"joint": f.path, "peak": n(f.peak), "at": n(f.at), "end": n(f.end), "ringing": ringing, "effort": f.effort.map(n)})
        })
        .collect();
    let spans: Vec<J> = r
        .spans
        .iter()
        .map(|sp| json!({"line": sp.name, "site": sp.site, "least": n(sp.least.1), "least_at": n(sp.least.0), "most": n(sp.most.1), "most_at": n(sp.most.0)}))
        .collect();
    json!({"joints": joints, "held": r.held, "spans": spans, "still": r.still.map(n)})
}

/// Each actuator's basis: `rigid` (the rigid model's inverse dynamics),
/// `command` (it moves simulated coordinates alone: its command), or, where
/// the file's simulation releases its joint for its springs, `spring`: what
/// it supplies is its spring's force against its deviation, as its claims
/// read it (Sprung), at the moment t (Err) or over the span (Ok), its limits
/// at the commanded speed; the rigid model's demand then under `rigid`. One
/// whose coordinates move what the springs release is refused.
fn supplied(actuators: &mut J, model: &Model, sim: Option<&dynamics::simulate::Simulation>, at: Result<(f64, f64), f64>) {
    let Some(list) = actuators.as_array_mut() else { return };
    for (a, entry) in list.iter_mut().enumerate() {
        let commanded = model.commanded(a);
        let Some(sim) = sim else {
            entry["basis"] = json!(if commanded { "command" } else { "rigid" });
            continue;
        };
        let (i, k) = match dynamics::Sprung::of(sim, commanded, &model.actuators[a].path) {
            Err(why) => {
                entry["basis"] = json!("rigid");
                entry["refused"] = json!(why);
                continue;
            }
            Ok(dynamics::Sprung::No) => {
                entry["basis"] = json!(if commanded { "command" } else { "rigid" });
                continue;
            }
            Ok(dynamics::Sprung::Spring { i, k }) => (i, k),
        };
        let mut rigid = entry.take();
        rigid.as_object_mut().map(|o| o.remove("path"));
        let mut j = json!({"path": model.actuators[a].path, "basis": "spring"});
        // Its effort, rate, limits and margin at moment t, as the spring reads it.
        let read = |t: f64| -> Result<(f64, std::sync::Arc<Moment>), String> {
            let m = model.at(t);
            if let Some(why) = &m.unsolved {
                return Err(why.clone());
            }
            Ok((dynamics::spring_effort(sim, i, k, t)?, m))
        };
        let margin = |e: f64, m: &Moment| (e - m.lo[a]).min(m.hi[a] - e);
        let row = |e: f64, m: &Moment| json!({"effort": n(e), "t": n(m.t), "rate": n(m.rate[a]), "lo": n(m.lo[a]), "hi": n(m.hi[a]), "margin": n(margin(e, m))});
        match at {
            Err(t) => match read(t) {
                Ok((e, m)) => {
                    for (key, v) in row(e, &m).as_object().unwrap() {
                        j[key] = v.clone();
                    }
                }
                Err(why) => j["refused"] = json!(why),
            },
            Ok(span) => {
                let steps = dynamics::spring_steps(model, a, sim, i, k, span);
                let mids = steps.ends.windows(2).map(|w| ink_core::eval::span_mid((w[0], w[1])));
                let moments: Result<Vec<(f64, std::sync::Arc<Moment>)>, (f64, String)> =
                    [span.0].into_iter().chain(mids).chain([span.1]).map(|t| read(t).map_err(|why| (t, why))).collect();
                match moments {
                    Err((t, why)) => j["refused"] = json!({"t": n(t), "why": why}),
                    Ok(v) => {
                        if let Some((e, m)) = peak(v.iter(), |x| x.0.abs()) {
                            j["peak"] = row(*e, m);
                        }
                        if let Some((e, m)) = v
                            .iter()
                            .filter(|(_, m)| m.lo[a] <= m.hi[a] && (m.lo[a].is_finite() || m.hi[a].is_finite()))
                            .min_by(|x, y| margin(x.0, &x.1).total_cmp(&margin(y.0, &y.1)))
                        {
                            j["margin"] = row(*e, m);
                        }
                        j["rms"] = sim.rms(i, span).map_or(J::Null, |r| n(k * r));
                        // Known either way to within its stiffness times its
                        // following error's band at the span's end (measured).
                        j["band"] = n(k * sim.band(i, span.1));
                        j["band_kind"] = json!("measured");
                        j["loose"] = J::Array(joined(&steps.loose).into_iter().map(|(a, b)| bounds(a, b)).collect());
                        j["partial"] = json!(!steps.loose.is_empty());
                    }
                }
            }
        }
        j["method"] = json!(SPRING);
        j["rigid"] = rigid;
        *entry = j;
    }
}

/// Intervals sorted, those meeting or overlapping joined: a run of
/// adjacent steps is one interval.
fn joined(v: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut v = v.to_vec();
    v.sort_by(|x, y| x.0.total_cmp(&y.0));
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(v.len());
    for (a, b) in v {
        match out.last_mut() {
            Some(last) if a <= last.1 => last.1 = last.1.max(b),
            _ => out.push((a, b)),
        }
    }
    out
}

/// `ink effort --modes [--t T]`: the mechanism's natural frequencies at the
/// pose at moment t (dynamics::modes), lowest first, each with what moves.
fn modes(s: &mut ModelView, reload: &ink_lang::Reload, target: &str) -> Result<J, InkError> {
    let m = dynamics::modes::modes(reload, target, s.ctx.t)?;
    // Each entry's own unit: a turn's motion in degrees, a slide's in mm; a
    // turn's stiffness a torque per radian, a slide's N/mm.
    let modes: Vec<J> = m
        .modes
        .iter()
        .map(|x| {
            let moves: Vec<J> = x
                .moves
                .iter()
                .map(|v| match v.motion {
                    Some(mo) => json!({"joint": v.joint, "motion": n(mo), "energy": n(v.energy)}),
                    None => json!({"body": v.joint, "energy": n(v.energy)}),
                })
                .collect();
            json!({"frequency": n(x.frequency), "moves": moves})
        })
        .collect();
    let springs: Vec<J> = m
        .springs
        .iter()
        .map(|sp| {
            let per = if sp.turn { json!("radian") } else { J::Null };
            json!({"what": sp.what, "stiffness": n(sp.stiffness), "per": per, "site": sp.site})
        })
        .collect();
    let model = &m.model;
    let res = if model.measured > 0 { json!([n(model.res.0), n(model.res.1)]) } else { J::Null };
    Ok(json!({
        "target": target,
        "t": n(m.t),
        "pose": m.pose.iter().map(|(name, v)| json!({"name": name, "value": n(*v)})).collect::<Vec<_>>(),
        "modes": modes,
        "springs": springs,
        "released": m.released,
        "rigid": m.rigid,
        "carried": model.carried.iter().map(|(t, mass)| json!({"target": t, "mass": n(*mass)})).collect::<Vec<_>>(),
        "deforming": model.deforming.iter().map(|(t, mass)| json!({"target": t, "mass": mass.map(n)})).collect::<Vec<_>>(),
        "method": MODES,
        "accuracy": J::Null,
        "masses": {"bodies": model.measured, "published": model.published, "spin_inertia": model.spin_inertia, "res": res, "unsettled": model.unsettled},
    }))
}

/// Over a span where groups are read with their spin at rest, round
/// (motion.basis): the most their departure from round moves this effort,
/// at any angle of their spins.
fn bound(j: &mut J, over: impl Iterator<Item = f64>) {
    let b = over.fold(0.0, f64::max);
    if b > 0.0 {
        j["bound"] = n(b);
    }
}

/// What an actuator's effort leaves out (deforming bodies moving with its
/// coordinates, whose inertia is not counted): its claim is refused.
fn leaves_out(j: &mut J, a: &dynamics::Actuator) {
    if !a.leaves_out.is_empty() {
        j["leaves_out"] = json!(a.leaves_out);
    }
}

/// The couplings passing power for an actuator that give no efficiency:
/// its effort leaves their loss out, and its claim is refused.
fn loss_unknown(j: &mut J, couplings: Vec<&str>) {
    if !couplings.is_empty() {
        j["loss_unknown"] = json!(couplings);
    }
}

/// A reaction's path, what it is not found for or leaves out, what moves
/// with it unseparated.
fn reaction(model: &Model, r: usize) -> J {
    let re = &model.reactions[r];
    let mut j = json!({"path": re.path});
    for (key, list) in [("without_mass", &re.without_mass), ("unbounded", &re.unbounded), ("leaves_out", &re.leaves_out), ("moves_with", &re.moves_with)] {
        if !list.is_empty() {
            j[key] = json!(list);
        }
    }
    j
}

/// A coupling's lines' tensions at moment m, or why they are not found;
/// whether its loss is unknown there.
fn coupling(model: &Model, m: &Moment, k: usize, lines: &std::ops::Range<usize>) -> J {
    let mut j = match model.tension_refusal(m, k) {
        Some(why) => json!({"refused": why}),
        None => json!({"tensions": m.tensions[lines.clone()].iter().copied().map(n).collect::<Vec<_>>()}),
    };
    if m.lossless.contains(&k) {
        j["loss"] = json!(if model.gives_efficiency(model.couplings()[k].0) {
            "unknown: its parts ride a carrier that moves, and its loss in that frame is not found"
        } else {
            "unknown: no efficiency given"
        });
    }
    j
}

/// What the couplings' tensions leave out: bodies moving with the joints
/// whose balances give them that no actuated or balancing joint carries.
fn unseparated(j: &mut J, model: &Model) {
    if !model.unseparated.is_empty() {
        j["unseparated"] = json!(model.unseparated);
    }
}

/// Per actuator, per coordinate, the coordinates with unpowered effort,
/// each reaction and each coupling, at one moment.
fn at_moment(model: &Model, m: &Moment) -> (Vec<J>, Vec<J>, Vec<J>, Vec<J>, Vec<J>) {
    let acts = model
        .actuators
        .iter()
        .enumerate()
        .map(|(k, a)| {
            let mut j = json!({"path": a.path, "effort": n(m.effort[k]), "rate": n(m.rate[k]), "lo": n(m.lo[k]), "hi": n(m.hi[k]), "margin": n(m.margin(k)),
                   "reflected_inertia": n(m.reflected[k])});
            leaves_out(&mut j, a);
            loss_unknown(&mut j, model.unknown_losses(k, m));
            if m.loss_bound[k] > 0.0 {
                j["loss_bound"] = n(m.loss_bound[k]);
            }
            if let Some(why) = model.refusal(k, m) {
                j["margin"] = J::Null;
                j["refused"] = json!(why);
            }
            j
        })
        .collect();
    let coords = model
        .coords
        .iter()
        .enumerate()
        .map(|(k, c)| json!({"name": c.name, "value": n(m.q[k][0]), "rate": n(m.q[k][1]), "acceleration": n(m.q[k][2]), "effort": n(m.f[k])}))
        .collect();
    let unpowered = model.coords.iter().zip(&m.unpowered).filter_map(|(c, u)| Some(json!({"coordinate": c.name, "effort": n((*u)?)}))).collect();
    let reactions = (0..model.reactions.len())
        .map(|r| {
            let mut j = reaction(model, r);
            match model.reaction_refusal(r, m) {
                Some(why) => (j["force"], j["moment"], j["refused"]) = (J::Null, J::Null, json!(why)),
                None => {
                    let (f, mo) = m.reactions[r].local();
                    (j["force"], j["moment"]) = (json!(f.map(n)), json!(mo.map(n)));
                }
            }
            j
        })
        .collect();
    let couplings = model
        .couplings()
        .into_iter()
        .enumerate()
        .map(|(k, (name, lines))| {
            let mut j = coupling(model, m, k, &lines);
            j["relation"] = json!(name);
            unseparated(&mut j, model);
            j
        })
        .collect();
    (acts, coords, unpowered, reactions, couplings)
}

/// The moment of moments in time order where `size` is largest: the
/// earliest of those equal within their rounding, so which of two equal
/// peaks (a swing's, mirrored) is shown does not turn on their last bits.
fn peak<T>(moments: impl Iterator<Item = T>, size: impl Fn(&T) -> f64) -> Option<T> {
    moments
        .fold(None, |best: Option<(f64, T)>, x| {
            let s = size(&x);
            match best {
                Some((b, y)) if !(s > b + crate::model::ROUNDING * b.abs()) => Some((b, y)),
                _ => Some((s, x)),
            }
        })
        .map(|(_, x)| x)
}

/// The same over a span's steps (each its width and the moment at its
/// middle): each actuator's peak effort (with its moment and rate), its
/// RMS, its least margin; each coordinate's peak effort; the peak of
/// each unpowered effort; each reaction where its force, and where its
/// moment, is largest, as the steps read it (Model::shaken), with the most
/// the groups read at rest there can move it, or the first moment it is
/// refused at.
fn over_span(model: &Model, span: (f64, f64), at: &[(f64, std::sync::Arc<Moment>)]) -> (Vec<J>, Vec<J>, Vec<J>, Vec<J>, Vec<J>) {
    let peak_of = |v: &dyn Fn(&Moment) -> Option<f64>| peak(at.iter().filter_map(|(_, m)| Some((v(m)?, m.clone()))), |x| x.0.abs());
    let acts = model
        .actuators
        .iter()
        .enumerate()
        .map(|(k, a)| {
            let mut j = json!({"path": a.path});
            leaves_out(&mut j, a);
            let mut unknown: Vec<&str> = Vec::new();
            for (_, m) in at {
                unknown.extend(model.unknown_losses(k, m).into_iter().filter(|c| !unknown.contains(c)).collect::<Vec<_>>());
            }
            loss_unknown(&mut j, unknown);
            let most = at.iter().map(|(_, m)| m.loss_bound[k]).fold(0.0, f64::max);
            if most > 0.0 {
                j["loss_bound"] = n(most);
            }
            let refusal = at.iter().find_map(|(_, m)| Some((m.t, model.refusal(k, m)?)));
            let refused = refusal.is_some();
            if let Some((t, why)) = refusal {
                j["refused"] = json!({"t": n(t), "why": why});
            }
            if let Some((e, m)) = peak_of(&|m: &Moment| Some(m.effort[k])) {
                j["peak"] = json!({"effort": n(e), "t": n(m.t), "rate": n(m.rate[k]), "reflected_inertia": n(m.reflected[k])});
            }
            // The RMS of the same moments, each weighing its step's width: a
            // mean needs the whole span, so none where a step is unchecked (a
            // peak over the rest is still at least that much).
            let (sum, time) = at.iter().fold((0.0, 0.0), |(s, w), (dt, m)| (s + dt * m.effort[k] * m.effort[k], w + dt));
            let whole = time >= (span.1 - span.0) * (1.0 - 1e-12);
            j["rms"] = if refused || !whole { J::Null } else { n((sum / time).sqrt()) };
            // Where its limits cross it has no margin (refused, saying so), nor where it states none (`(actuator)`).
            if !refused
                && let Some((_, m)) =
                    at.iter().filter(|(_, m)| m.lo[k] <= m.hi[k] && (m.lo[k].is_finite() || m.hi[k].is_finite())).map(|(_, m)| (m.margin(k), m)).min_by(|x, y| x.0.total_cmp(&y.0))
            {
                j["margin"] = json!({"margin": n(m.margin(k)), "t": n(m.t), "effort": n(m.effort[k]), "rate": n(m.rate[k]), "lo": n(m.lo[k]), "hi": n(m.hi[k])});
            }
            bound(&mut j, at.iter().map(|(_, m)| m.bound[k]));
            j
        })
        .collect();
    let coords = model
        .coords
        .iter()
        .enumerate()
        .map(|(k, c)| {
            let mut j = json!({"name": c.name});
            if let Some((e, m)) = peak_of(&|m: &Moment| Some(m.f[k])) {
                j["peak"] = json!({"effort": n(e), "t": n(m.t), "value": n(m.q[k][0]), "rate": n(m.q[k][1]), "acceleration": n(m.q[k][2])});
            }
            bound(&mut j, at.iter().map(|(_, m)| m.bound_f[k]));
            j
        })
        .collect();
    let unpowered = model
        .coords
        .iter()
        .enumerate()
        .filter_map(|(k, c)| {
            let (e, m) = peak_of(&|m: &Moment| m.unpowered[k])?;
            let mut j = json!({"coordinate": c.name, "peak": {"effort": n(e), "t": n(m.t)}});
            bound(&mut j, at.iter().map(|(_, m)| m.bound_u[k]));
            Some(j)
        })
        .collect();
    let reactions = (0..model.reactions.len())
        .map(|r| {
            let mut j = reaction(model, r);
            if let Some((t, why)) = at.iter().find_map(|(_, m)| Some((m.t, model.reaction_refusal(r, &model.at(m.t))?))) {
                j["refused"] = json!({"t": n(t), "why": why});
                return j;
            }
            // Each moment's wrench in its path's frame, and its bound there.
            let read: Vec<(f64, [f64; 3], [f64; 3], f64, f64)> = at
                .iter()
                .map(|(_, m)| {
                    let (w, shakes) = model.shaken(span, None, r, m.t);
                    let w = &w.reactions[r];
                    let (f, mo) = w.local();
                    let o = [w.frame[9], w.frame[10], w.frame[11]];
                    let bf: f64 = shakes.iter().map(|s| s.0).sum();
                    let bm: f64 = shakes.iter().map(|&(f, mo, p)| mo + f * ink_core::v3::norm([p[0] - o[0], p[1] - o[1], p[2] - o[2]])).sum();
                    (m.t, f, mo, bf, bm)
                })
                .collect();
            for (key, part) in [("peak_force", 0), ("peak_moment", 1)] {
                let size = |x: &(f64, [f64; 3], [f64; 3], f64, f64)| ink_core::v3::norm([x.1, x.2][part]);
                if let Some(x) = peak(read.iter(), |x| size(x)) {
                    j[key] = json!({"t": n(x.0), "force": x.1.map(n), "moment": x.2.map(n)});
                }
            }
            let (bf, bm) = read.iter().fold((0.0f64, 0.0f64), |b, x| (b.0.max(x.3), b.1.max(x.4)));
            if bf > 0.0 || bm > 0.0 {
                j["bound"] = json!({"force": n(bf), "moment": n(bm)});
            }
            j
        })
        .collect();
    // Each coupling where its tightest line is tightest, or where it is
    // first refused.
    let couplings = model
        .couplings()
        .into_iter()
        .enumerate()
        .map(|(k, (name, lines))| {
            let mut j = json!({"relation": name});
            if let Some((t, why)) = at.iter().find_map(|(_, m)| Some((m.t, model.tension_refusal(m, k)?))) {
                j["refused"] = json!({"t": n(t), "why": why});
                return j;
            }
            let tightest = |m: &Moment| m.tensions[lines.clone()].iter().fold(0.0f64, |p, x| p.max(*x));
            if let Some((_, m)) = peak(at.iter(), |x| tightest(&x.1)) {
                j["peak"] = json!({"t": n(m.t), "tensions": m.tensions[lines.clone()].iter().copied().map(n).collect::<Vec<_>>()});
            }
            if at.iter().any(|(_, m)| m.lossless.contains(&k)) {
                j["loss"] = json!(if model.gives_efficiency(model.couplings()[k].0) {
                    "unknown: its parts ride a carrier that moves, and its loss in that frame is not found"
                } else {
                    "unknown: no efficiency given"
                });
            }
            unseparated(&mut j, model);
            j
        })
        .collect();
    (acts, coords, unpowered, reactions, couplings)
}
