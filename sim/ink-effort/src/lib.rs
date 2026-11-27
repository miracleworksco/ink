//! The effort plug: what an actuator (a joint with `effort-limits`) needs
//! through the written motion, as facts for the claims core.
//!
//! The solver is ink-dynamics, entered as `ink effort FILE --target TARGET`
//! enters it (`ink_dynamics::read`): a run is the model and the moments
//! that command reads, over the file's span or at its one moment. This
//! adapter starts it the first time a fact it offers is asked for, keeps
//! it, and reads every fact off it. On an actuator at instance path P that turns it offers
//! `peak-torque on P` (N*mm) and `speed-at-peak-torque on P` (deg/s); on one
//! that slides, `peak-force on P` (N) and `speed-at-peak-force on P` (mm/s):
//! the language's units, which are the solver's own, so a part's rating
//! (ink-ratings) can be held against them. Each is a size, without its sign.
//!
//! A peak is an estimate, and its band is what the solver itself says of
//! it. Either way, how far the effort it read may be from the real one
//! (`bound`: bodies read with their spin at rest, round; `loss_bound`: a
//! coupling's loss resting on tensions that leave bodies out): within that
//! of a limit the old claim on it is unknown, and so is a rating here.
//! Upward only, over a motion, what a peak between the moments read may add:
//! the solver bounds nothing between them, and chooses them so that an
//! effort lies within 4% of its peak of the line through its neighbours
//! (`ink_dynamics::OFF_LINE`; where it turns back near its extreme, a
//! moment of its own), reading a peak at most that low. So a rating within
//! 4% above the peak read is not settled, where `ink require`, which looks
//! at the moments alone, says it holds. An effort the solver's own claim
//! refuses (bodies or a loss it leaves out) is open here too.
//!
//! Beside the peaks, each moment read is a tool of its own (`Effort::at`):
//! on an actuator's path P, and on a frame its part has (`P shaft`: the
//! joint's effort is taken as what the part it turns or slides passes,
//! wherever on it a rating names), `torque` and `speed` (or `force` and
//! `speed`) at that moment, with their sign. A moment's speed is a band:
//! from its own to its neighbours', so a limit read at it is read over
//! every speed the joint passes between the moments, and a speed stands
//! where the effort there is refused.
//!
//! A moment's tool also offers what the model works out at a part a maker
//! rates by load: on the part a joint moves, the joint's reaction
//! (`radial-load`, `axial-load`, and `moment-x`, `moment-y`, `moment-z`
//! about that instance's own axes: a carriage's); on a named support its
//! share; on a coupling's part `tension`, the tightest of its couplings'
//! lines (a screw's thrust on its nut). A rating is then data alone.
//!
//! Each fact is offered from the least it needs. A speed is the written
//! motion's: it needs no mass and no gravity, and is read from the drives
//! alone at moments of the motion's own (the span's ends and the middles of
//! 64 even steps, `Kin`) as long as no effort has been asked; nothing is
//! solved for it. An effort or a load needs the masses it rests on and a
//! gravity. The run is made whatever is missing (`ink_dynamics::
//! read_kinematic`), and what is missing takes only its own: a moving body
//! with no mass leaves open the efforts of the actuators it moves with (its
//! own joints' and those coupled to them) and the reactions that carry it,
//! each naming the body; an effort no such body is on is the one a design
//! without that body would get, to the digit. With no gravity declared no
//! effort is given. `ink effort`, the command, gives the same: a row with
//! no effort says why (`lacks`) and keeps its rate.
//!
//! A run is over the span the design was loaded to be asked over (a
//! command's `--t-span`, else the file's), and at one value of each param.
//! Where the model reads a ranged param no fact is given: one value says
//! nothing of the rest of the range, as the old claim on an effort says. A rating whose limit moves with speed (a
//! motor's curve, read at the shaft's speed) is then answered at every
//! moment read, by the core's `over`, and holds only where it holds at
//! each: the peak and the speed at the peak say nothing of a motor nearest
//! its curve elsewhere in the motion.

use ink_claims::{About, Claims, Covers, Fact, Strength, Tool, Unavailable, Unsettled, Verdict};
use ink_dynamics::{Moment, OFF_LINE, Read};
use ink_lang::view::ModelView;
use std::cell::{Cell, OnceCell, RefCell};
use std::collections::HashMap;

pub struct Effort<'a> {
    view: &'a RefCell<ModelView>,
    target: String,
    /// Each actuator under the target, by instance path, whether it turns,
    /// and the frames of the part it moves, where it moves a body.
    actuators: Vec<(String, bool, Vec<String>)>,
    /// Each joint under the target, by instance path, with the frames of the part it moves: where a rating
    /// reads the load its joint holds.
    joints: Vec<(String, Vec<String>, usize)>,
    /// Each actuator's part's frames placed in its joint's frame (R row major, then T), by actuator path and
    /// frame name: where a frame's effort is read. A frame that is no constant placement is not here.
    placed: HashMap<(String, String), [f64; 12]>,
    /// The ends of the couplings under the target (a screw and its nut, a belt's pulleys, meshing gears), by
    /// instance path: where a `tension` is offered.
    ends: Vec<String>,
    /// Each of `joints`: its drive's rate as a field of t (`ink_dynamics::drive_rate`), and whether it turns.
    moved: Vec<(Option<ink_dynamics::Rate>, bool)>,
    /// The reactions the design divides among supports (`(reaction PATH (by S ...))`), but one "divided"
    /// among its own part alone (`(by PATH)`: a claim's way to read a joint's reaction in its part's load
    /// frame), which is the joint's.
    divided: Vec<ink_dynamics::Divided>,
    /// What the model works out at a moment, kept: a joint's reaction by its path, a division's shares by its place.
    reactions: RefCell<HashMap<(String, usize), Result<([f64; 3], [f64; 3]), String>>>,
    shares: RefCell<HashMap<(usize, usize), Result<(Vec<[f64; 6]>, Option<f64>), String>>>,
    run: OnceCell<Result<Read, String>>,
    runs: Cell<usize>,
    /// How many joints at how many moments `books` added up the balance of apart from the solve (the rest
    /// it could not count: a carried belt, a jump in acceleration, a refused reaction).
    balanced: Cell<usize>,
    /// And of those, how many had their moments balanced too.
    balanced_moments: Cell<usize>,
    /// Each actuator's drive's rate as a field of t (`ink_dynamics::drive_rate`), and the motion read from
    /// them alone (`Kin`), made when a speed is first asked and no run has started.
    /// What the run leaves out of each actuator's facts, once it is asked (`begun`).
    left: RefCell<HashMap<usize, Vec<String>>>,
    rates: Vec<Option<ink_dynamics::Rate>>,
    kin: OnceCell<Option<Kin>>,
}

/// The written motion at moments of its own: their times, and each actuator's rate at each. It needs no
/// model, so no mass and no gravity: a speed is the motion's whatever else the design lacks.
struct Kin {
    t: Vec<f64>,
    rate: Vec<Vec<f64>>,
    span: Option<(f64, f64)>,
}

/// The even steps the motion is read at for its speeds alone: their middles, and the span's two ends.
const KIN_STEPS: usize = 64;

/// One moment of an effort's run, as a tool (`Effort::at`).
pub struct AtMoment<'a> {
    effort: &'a Effort<'a>,
    moment: usize,
    /// It was made before the run: its moment is the written motion's own (`Kin`), not the run's.
    early: bool,
}

impl<'a> Effort<'a> {
    /// The plug for the assembly or instance path `target` of a loaded
    /// design, which it shares with the other tools asked beside it. It
    /// reads which joints are actuators, and solves nothing.
    pub fn new(view: &'a RefCell<ModelView>, target: &str) -> Result<Effort<'a>, String> {
        let mut actuators: Vec<(String, bool, Vec<String>)> = vec![];
        let mut rates = vec![];
        let mut joints: Vec<(String, Vec<String>, usize)> = vec![];
        let mut placed = HashMap::new();
        let mut moved = vec![];
        let ctx = view.borrow().ctx.clone();
        let doc = &mut view.borrow_mut().doc;
        let items = doc.items(target)?;
        for link in items.iter().flat_map(|item| &item.links) {
            let Some(joint) = link.joint.as_ref().filter(|_| !joints.iter().any(|(path, ..)| *path == link.path)) else { continue };
            // The body placed at the joint's own path, if one is.
            let body = items.iter().find(|item| item.label == link.path).and_then(|item| item.instance.as_ref()).and_then(|(name, _)| doc.body(name));
            let frames: Vec<String> = body.map_or(vec![], |b| b.frames.iter().map(|f| f.0.clone()).collect());
            if ink_dynamics::joint::effort(joint).is_some() {
                // Each frame of the part, in the joint's frame: the part's placement below the motion, then the frame,
                // at the design's params.
                let num = |x: &ink_lang::Xf| -> Option<[f64; 12]> { Some(ink_dynamics::placement_at(doc, x, &ctx)).filter(|v| v.iter().all(|e| e.is_finite())) };
                if let (Some(b), Some(below)) = (body, num(&joint.below)) {
                    for (name, x) in &b.frames {
                        if let Some(f) = num(x) {
                            let at = |e: &[f64; 12], i: usize, j: usize| e[3 * i + j];
                            let mut out = [0.0; 12];
                            for i in 0..3 {
                                for j in 0..3 {
                                    out[3 * i + j] = (0..3).map(|k| at(&below, i, k) * at(&f, k, j)).sum();
                                }
                                out[9 + i] = (0..3).map(|k| at(&below, i, k) * f[9 + k]).sum::<f64>() + below[9 + i];
                            }
                            placed.insert((link.path.clone(), name.clone()), out);
                        }
                    }
                }
                actuators.push((link.path.clone(), joint.turn, frames.clone()));
                rates.push(ink_dynamics::drive_rate(doc, joint));
            }
            joints.push((link.path.clone(), frames, joint.axis));
            moved.push((ink_dynamics::drive_rate(doc, joint), joint.turn));
        }
        let divided = ink_dynamics::divided(doc).into_iter().filter(|d| d.supports != [d.path.clone()]).collect();
        let ends = ink_dynamics::coupling::related(doc, target).iter().filter(|r| r.couple.is_some()).flat_map(|r| [r.a.clone(), r.b.clone()]).collect();
        Ok(Effort {
            view,
            target: target.into(),
            actuators,
            joints,
            placed,
            moved,
            ends,
            divided,
            reactions: RefCell::default(),
            shares: RefCell::default(),
            run: OnceCell::new(),
            runs: Cell::new(0),
            balanced: Cell::new(0),
            balanced_moments: Cell::new(0),
            left: RefCell::default(),
            rates,
            kin: OnceCell::new(),
        })
    }

    /// The time each moment of its run weighs in a mean over the motion: the step it is the middle of, and
    /// none for the span's two ends (the midpoint rule, as its RMS). Empty where the design has no span, or
    /// the run does not start: at one moment nothing is a mean. It starts the run.
    pub fn steps(&self) -> Vec<f64> {
        match self.run() {
            Ok(run) if run.span.is_some() => run.at.iter().map(|(step, _)| *step).collect(),
            _ => vec![],
        }
    }

    /// The parts its model puts a demand on that a maker rates, each with why: what drives the motion (an
    /// actuator's joint), what a coupling passes force through (a screw and its nut, a belt's pulleys, meshing
    /// gears), and what the design names as holding something (`(supports ...)`). Read off the design: nothing
    /// is solved. A printed bracket, a plain shaft or a carriage is none of these: the model works out no
    /// load a datasheet could be held against.
    pub fn loaded(&self) -> Vec<(String, String)> {
        let mut out: Vec<(String, String)> = self.actuators.iter().map(|(path, ..)| (path.clone(), "it drives the motion".to_string())).collect();
        let doc = &mut self.view.borrow_mut().doc;
        for related in ink_dynamics::coupling::related(doc, &self.target).iter().filter(|r| r.couple.is_some()) {
            for (part, other) in [(&related.a, &related.b), (&related.b, &related.a)] {
                out.push((part.clone(), format!("a coupling passes force between it and {other}")));
            }
        }
        for division in &self.divided {
            out.extend(division.supports.iter().map(|s| (s.clone(), format!("it holds what {} carries", division.path))));
        }
        let mut seen: Vec<(String, String)> = vec![];
        for (path, why) in out {
            if !seen.iter().any(|(p, _)| *p == path) {
                seen.push((path, why));
            }
        }
        seen
    }

    /// The joints that hold a coupling's pull, each with why: the one each end of a coupling rides (a motor's
    /// rotor under the pulley on its shaft, a carriage a belt is clamped to), but one whose load the design
    /// divides among supports, which then hold it. The joint's reaction has the pull: a rating of its part
    /// that reads a load across or along it (`radial-load`, `axial-load`, a moment) checks it; a part
    /// rated for its torque alone is not checked for it.
    pub fn holds(&self) -> Vec<(String, String)> {
        let doc = &mut self.view.borrow_mut().doc;
        let mut out: Vec<(String, String)> = vec![];
        for related in ink_dynamics::coupling::related(doc, &self.target).iter().filter(|r| r.couple.is_some()) {
            for end in [&related.a, &related.b] {
                let Some(joint) = ink_dynamics::riding(doc, &self.target, end) else { continue };
                if self.divided.iter().any(|d| d.path == joint) || out.iter().any(|(j, _)| *j == joint) {
                    continue;
                }
                let on = if joint == *end { "it".to_string() } else { format!("{end}, which it carries") };
                out.push((joint, format!("a coupling pulls on {on}, across or along its joint")));
            }
        }
        let _ = doc;
        // And every actuator's own joint: what it carries bears on its output across its axis, a force and a
        // bending moment the model works out. Said where the bending is over a quarter of the most torque
        // it turns with (or that is nothing), with the numbers: a hobby servo's spline is not a bearing.
        if let Some(Ok(run)) = self.run.get() {
            for (a, act) in run.model.actuators.iter().enumerate() {
                if self.divided.iter().any(|d| d.path == act.path) || out.iter().any(|(j, _)| *j == act.path) {
                    continue;
                }
                let (Some(r), Some((_, _, axis))) = (run.model.reactions.iter().position(|r| r.path == act.path), self.joints.iter().find(|(p, ..)| *p == act.path)) else {
                    continue;
                };
                let mut most: Option<(f64, f64, f64)> = None;
                for (_, m) in run.at.iter().filter(|(_, m)| m.unsolved.is_none()) {
                    if run.model.reaction_refused_lossless(r, m).is_some() || !run.model.blocks_reaction(r).is_empty() {
                        most = None;
                        break;
                    }
                    let (force, moment) = m.reactions[r].local();
                    let across = |v: [f64; 3]| (0..3).filter(|k| k != axis).map(|k| v[k]).fold(0.0, f64::hypot);
                    if most.is_none_or(|b| across(moment) > b.1) {
                        most = Some((across(force), across(moment), m.t));
                    }
                }
                let turned = run.at.iter().map(|(_, m)| m.effort[a].abs()).fold(0.0, f64::max);
                if let Some((force, bending, t)) = most.filter(|b| b.1 > 0.25 * turned && b.1 > 1.0) {
                    let round = |x: f64| (x * 1e3).round() / 1e3;
                    out.push((
                        act.path.clone(),
                        format!(
                            "what it carries bends its output with {} N*mm across its axis and pushes {} N across it at t = {}, beside the {} N*mm at most it turns with",
                            round(bending),
                            round(force),
                            round(t),
                            round(turned)
                        ),
                    ));
                }
            }
        }
        out
    }

    /// What the motion needs that nothing in the design supplies or holds, each said: a coordinate along
    /// which the motion as written needs an effort no actuator gives (a weight on a slide that nothing
    /// drives and nothing joins to what does: it would fall, or the motor that should lift it does no
    /// work for it), with the most it needs over the moments read; a body moved by a formula of t through
    /// no joint; and a load on a part that is on no joint. With any of these the efforts and loads the
    /// model gives are of a design less than the one meant. It starts the run where a gravity is declared.
    pub fn unheld(&self) -> Vec<String> {
        let mut out = vec![];
        if ink_dynamics::gravity(&self.view.borrow().doc).is_some()
            && let Ok(run) = self.run()
        {
            for (k, coord) in run.model.coords.iter().enumerate() {
                let most = run.at.iter().filter_map(|(_, m)| Some((m.unpowered[k]?, m.t))).max_by(|a, b| a.0.abs().total_cmp(&b.0.abs()));
                if let Some((effort, t)) = most {
                    out.push(format!(
                        "the motion as written needs {} (N on a slide, N*mm per degree turned) along {} at t = {} that no actuator supplies: nothing drives what moves \
                         along it, and nothing joins it to what does, so the efforts found leave it out. Join it: place it on the frame of the part that moves it, write \
                         its motion as a formula of that part's ((slide y (* r (/ angle 1rad))) for what a drum of radius r winds) and relate the two by what is \
                         between them ((cord 10mm)), or give its joint an (actuator)",
                        ((effort * 1e3).round() / 1e3).abs(),
                        coord.name,
                        (t * 1e3).round() / 1e3
                    ));
                }
            }
            // A joint that is no actuator's and does not move is free along its axis still: what it carries
            // must need nothing there, or it would slide or swing away (a weight on a slide at rest, an arm on
            // a hinge posed at an angle). Read at the first moment.
            let doc = &self.view.borrow().doc;
            for (j, (path, _, axis)) in self.joints.iter().enumerate() {
                let (rate, turn) = self.moved[j];
                let still = rate.is_some_and(|rate| doc.g.as_const(rate) == Some(0.0));
                if !still || self.actuators.iter().any(|(p, ..)| p == path) {
                    continue;
                }
                let Some(r) = run.model.reactions.iter().position(|r| r.path == *path) else { continue };
                let Some((_, m)) = run.at.first() else { continue };
                if run.model.reaction_refused_lossless(r, m).is_some() || !run.model.blocks_reaction(r).is_empty() {
                    continue;
                }
                let (force, moment) = m.reactions[r].local();
                let (along, all) =
                    if turn { (moment[*axis], moment.iter().map(|x| x * x).sum::<f64>().sqrt()) } else { (force[*axis], force.iter().map(|x| x * x).sum::<f64>().sqrt()) };
                // (Past what a mass measured to a part in a thousand may leave: a balanced roller is not unheld.)
                if along.abs() > 1e-3 * all + 0.01 {
                    out.push(format!(
                        "{path} is free to {} and does not move, and what it carries needs {} {} there that nothing supplies: nothing holds it, and it would not stay. Say \
                         what holds it: an (actuator) or (effort-limits LO HI) on its joint (a motor, a brake), a coupling to what does, or place it without the joint",
                        if turn { "turn" } else { "slide" },
                        ((along * 1e3).round() / 1e3).abs(),
                        if turn { "N*mm about its axis" } else { "N along it" }
                    ));
                }
            }
            // An actuator that turns or slides through the whole motion with no effort at all, while something
            // in the design has mass: it drives nothing. What it was meant to move is on another part (its
            // case, not its output).
            // (Nothing: under a hundred-thousandth of the largest force any joint holds, at a millimetre.)
            let held = run.at.iter().flat_map(|(_, m)| m.reactions.iter()).map(|w| w.force.iter().map(|x| x * x).sum::<f64>().sqrt()).fold(0.0, f64::max);
            let massive = held > 1e-9;
            for (a, act) in run.model.actuators.iter().enumerate() {
                let (effort, rate) = run.at.iter().fold((0.0f64, 0.0f64), |(e, r), (_, m)| (e.max(m.effort[a].abs()), r.max(m.rate[a].abs())));
                if massive && rate > 0.0 && effort <= 1e-5 * held && run.at.iter().all(|(_, m)| run.model.refusal_lossless(a, m).is_none()) {
                    out.push(format!(
                        "{} moves through the whole motion with no effort to speak of ({effort:.1e} at most): this motor drives nothing. What it was meant to move is not on the part it turns or \
                         slides: place that on the motor's moving frame ((frame-of MOTOR output), (frame-of MOTOR shaft)), not on its case",
                        act.path
                    ));
                }
            }
            if !run.model.driven_by_nothing.is_empty() {
                out.push(format!(
                    "these bodies are moved by a formula of t through no joint, so no actuator is charged for moving them: {}. Put each on a joint ((turn ...), (slide ...))",
                    run.model.driven_by_nothing.join(", ")
                ));
            }
        }
        // A joint moved by a motor's coordinate with nothing between them: the motion is right, the motor is
        // charged for what the joint moves, but whatever carries the force (a belt, a cord, a mesh) is not
        // in the model. Said wherever the run was made, with or without a gravity.
        if let Some(Ok(run)) = self.run.get() {
            for (joint, by) in &run.model.unlinked {
                out.push(format!(
                    "{joint} moves by {}'s coordinate, and no part between them is in the model: nothing relates them (a belt, a cord, a mesh) and it does not ride \
                     {}'s moving frame. What such a part adds (its mass, its tension, its stiffness, its losses and its own rating) is not checked. Relate \
                     {joint} to what moves it by what is between them: (relate DRUM {joint} (cord 10mm)), (belt-span ...), (gear-mesh ...), (rack-mesh ...); \
                     or, if it is fixed to the motor's output, place it there ((frame-of MOTOR shaft)) with no joint of its own",
                    by.join(" and "),
                    by.join(" or ")
                ));
            }
        }
        out.extend(ink_dynamics::loads_to_ground(&mut self.view.borrow_mut().doc, &self.target));
        // A support in mid-air: it is given a share of a load and touches nothing that could take it on.
        let held: Vec<(String, Vec<String>)> = self.divided.iter().map(|d| (d.path.clone(), d.supports.clone())).collect();
        for (support, joint) in ink_dynamics::supports_floating(&mut self.view.borrow_mut(), &self.target, &held) {
            out.push(format!(
                "{support} is named as holding {joint} and touches nothing else: its box is clear of every other part by more than 1 mm, so nothing passes its load on. \
                 Place it in what holds it (a bearing in its housing), or model that part"
            ));
        }
        out
    }

    /// The answer's own books, added up from the numbers already worked out, each imbalance said: which
    /// part or joint, how much is unaccounted, and in what direction. Three sums, none made by the split
    /// or the solve they check:
    ///
    /// - Of each joint whose load the design divides among supports, at each moment read: the shares handed
    ///   out, turned into the world and added, against the joint's reaction as the model has it, force and
    ///   moment (the thrust the locating bearings carry among them is in no share, and is allowed for).
    /// - Of each support: the load its facts offer (`radial-load`, asked as a rating asks it) against the
    ///   sum of every share any division gave it. A part named in two divisions and offered one of them
    ///   shows here.
    /// - Of the drive: the power the actuators supply (effort by rate) against the rate of work the motion
    ///   takes along its coordinates. With no efficiency stated they are equal; with one, the supply is no
    ///   less.
    ///
    /// To the accuracy the facts carry: a part in a thousand of the sums, and the masses' 2% where a
    /// rounding of theirs enters. Empty where the run has not started. What it does not add up: each
    /// body's own balance of weight, inertia and loads against its joint's reaction, which is how the
    /// solver finds that reaction and would be the same sum twice.
    pub fn books(&self) -> Vec<String> {
        let Some(Ok(run)) = self.run.get() else { return vec![] };
        if run.weightless {
            return vec![];
        }
        let params = self.view.borrow().ctx.params.clone();
        let size = |v: [f64; 3]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
        let sub = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] { std::array::from_fn(|i| a[i] - b[i]) };
        let round = |x: f64| (x * 1e3).round() / 1e3;
        let dir = |v: [f64; 3]| {
            let n = size(v).max(1e-300);
            format!("({}, {}, {})", round(v[0] / n), round(v[1] / n), round(v[2] / n))
        };
        // The worst of each kind over the moments read: (how much, at t, said).
        let mut worst: Vec<(String, f64, String)> = vec![];
        let mut note = |key: String, much: f64, said: String| match worst.iter_mut().find(|w| w.0 == key) {
            Some(w) if w.1 >= much => {}
            Some(w) => (w.1, w.2) = (much, said),
            None => worst.push((key, much, said)),
        };
        let design = ink_claims::Design::default();
        for (moment, (_, m)) in run.at.iter().enumerate() {
            if m.unsolved.is_some() {
                continue;
            }
            for division in &self.divided {
                let Some(Ok(([given, needed, given_m, needed_m], axis))) = division.books(run, m.t, &params) else { continue };
                // Less what the joint's own drive carries along its axis: a turn's torque, a slide's force.
                let across = |v: [f64; 3], turns: bool| match axis {
                    Some((turn, u)) if turn == turns => {
                        let along = (0..3).map(|i| v[i] * u[i]).sum::<f64>();
                        std::array::from_fn(|i| v[i] - along * u[i])
                    }
                    _ => v,
                };
                let thrust = division.shares(run, m.t, &params).and_then(Result::ok).and_then(|(_, thrust)| thrust).unwrap_or(0.0);
                // (Less what its masses' proven bounds can move the thrust by.)
                let moved: f64 = division.share_moves(run, m.t, &params).iter().map(|(_, t)| t.abs()).sum();
                let (left, scale) = (across(sub(needed, given), false), size(needed).max(size(given)));
                if (size(left) - thrust).abs() > 1e-3 * scale + moved + 1e-6 {
                    note(
                        format!("f {}", division.path),
                        size(left),
                        format!(
                            "what {} needs held and what its supports ({}) were given differ by {} N along {} at t = {} (needed {} N, given {} N)",
                            division.path,
                            division.supports.join(", "),
                            round(size(left)),
                            dir(left),
                            round(m.t),
                            round(size(needed)),
                            round(size(given))
                        ),
                    );
                }
                // (A moment: where no thrust is left to the locating bearings, whose line is not the split's to say.)
                let (left, scale) = (across(sub(needed_m, given_m), true), size(needed_m).max(size(given_m)));
                if thrust == 0.0 && size(left) > 1e-3 * scale + 1e-3 * size(needed) + 1e-6 {
                    note(
                        format!("m {}", division.path),
                        size(left),
                        format!(
                            "the moment {} needs held and the moment its supports ({}) were given differ by {} N*mm about {} at t = {}",
                            division.path,
                            division.supports.join(", "),
                            round(size(left)),
                            dir(left),
                            round(m.t)
                        ),
                    );
                }
            }
            // Each support: what its fact offers against every share it was given.
            let mut seen: Vec<&String> = vec![];
            for (division, k) in self.divided.iter().flat_map(|d| (0..d.supports.len()).map(move |k| (d, k))) {
                let support = &division.supports[k];
                if seen.contains(&support) {
                    continue;
                }
                seen.push(support);
                let mut given = [0.0; 6];
                let mut whole = true;
                for (d2, k2) in self.divided.iter().filter_map(|d2| Some((d2, d2.supports.iter().position(|s| s == support)?))) {
                    match d2.given(run, m.t, &params, k2) {
                        Some(share) => (0..6).for_each(|i| given[i] += share[i]),
                        None => whole = false,
                    }
                }
                if !whole {
                    continue;
                }
                let assigned = if division.guides[k] { given[1].hypot(given[2]) } else { given[0].hypot(given[1]) };
                let claims = ink_claims::Claims::new(&design, &[]);
                let Some(offered) = self.load(&ink_claims::about("radial-load", support), moment, false, &claims).and_then(|f| f.value) else { continue };
                if (offered - assigned).abs() > 1e-3 * assigned.max(offered) + 1e-6 {
                    note(
                        format!("s {support}"),
                        (offered - assigned).abs(),
                        format!(
                            "{support} is offered {} N across it and was given {} N in all by the joints it is named as holding: {} N {} at t = {}",
                            round(offered),
                            round(assigned),
                            round((offered - assigned).abs()),
                            if offered < assigned { "is not counted" } else { "is counted twice" },
                            round(m.t)
                        ),
                    );
                }
            }
            // The drive: power in against the rate of work the motion takes.
            let supplied: f64 = (0..run.model.actuators.len())
                .map(|a| {
                    m.effort[a] * m.rate[a] * if self.actuators.iter().any(|(p, turn, _)| *p == run.model.actuators[a].path && *turn) { std::f64::consts::PI / 180.0 } else { 1.0 }
                })
                .sum();
            let taken: f64 = (0..m.f.len()).map(|k| m.f[k] * m.q[k][1]).sum();
            let scale: f64 = (0..m.f.len()).map(|k| (m.f[k] * m.q[k][1]).abs()).sum::<f64>().max(supplied.abs());
            let slack = 1e-3 * scale + (0..run.model.actuators.len()).map(|a| ((m.bound[a] + m.loss_bound[a]) * m.rate[a]).abs()).sum::<f64>() + 1e-6;
            let powered = m.unpowered.iter().all(Option::is_none)
                && supplied.is_finite()
                && taken.is_finite()
                && (0..run.model.actuators.len()).all(|a| run.model.refusal_lossless(a, m).is_none());
            let off = if run.model.lossy() { (taken - supplied).max(0.0) } else { (taken - supplied).abs() };
            if powered && off > slack {
                note(
                    "p".into(),
                    off,
                    format!(
                        "the actuators supply {} N*mm/s and the motion takes {} N*mm/s along its coordinates at t = {}: {} N*mm/s {}",
                        round(supplied),
                        round(taken),
                        round(m.t),
                        round(off),
                        if supplied < taken { "comes from nowhere" } else { "goes nowhere, with no coupling stating a loss" }
                    ),
                );
            }
        }
        // Each joint, apart from how its reaction was found: the weight of every body it carries, the loads
        // written on them, what couplings pull on them with from outside, and their masses times the
        // accelerations of their centres (by differences of where the centres are, a little before and after),
        // against the reaction. Read at no more than 33 of the moments. What it cannot count is left to the
        // joints it says nothing of: a carried belt, a body whose shape moves, a moment where the motion's
        // acceleration jumps.
        let gravity = ink_dynamics::gravity(&self.view.borrow().doc).unwrap_or([0.0; 3]);
        let loads = ink_dynamics::loads_written(&mut self.view.borrow_mut().doc, &self.target);
        let view = &*self.view.borrow();
        let every = (run.at.len() / 33).max(1);
        let h = run.span.map_or(1e-4, |(a, b)| 1e-4 * (b - a));
        let add = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] { std::array::from_fn(|i| a[i] + b[i]) };
        for (_, m) in run.at.iter().step_by(every).filter(|(_, m)| m.unsolved.is_none()) {
            for (r, reaction) in run.model.reactions.iter().enumerate() {
                if run.model.reaction_refused_lossless(r, m).is_some() || !run.model.blocks_reaction(r).is_empty() || reaction.path == self.target {
                    continue;
                }
                let Some((bodies, other)) = run.model.carried_masses(r, m.t) else { continue };
                if !other.is_empty() || bodies.is_empty() {
                    continue;
                }
                let Some(pull) = run.model.coupling_pull(r, m) else { continue };
                // m a of what it carries: second differences of the centres, at h and at 2 h, which agree
                // where the motion is smooth there.
                let at = |t: f64| run.model.carried_masses(r, t).map(|(bodies, _)| bodies);
                let (Some(b2), Some(b1), Some(a1), Some(a2)) = (at(m.t - 2.0 * h), at(m.t - h), at(m.t + h), at(m.t + 2.0 * h)) else { continue };
                // (centred, from before alone, from after alone: where an acceleration jumps at this moment
                // they differ, and none of them is the moment's.)
                let second = |x: &[(f64, [f64; 3])], y: &[(f64, [f64; 3])], z: &[(f64, [f64; 3])]| -> [f64; 3] {
                    std::array::from_fn(|i| (0..bodies.len()).map(|k| bodies[k].0 * 1e-6 * (x[k].1[i] - 2.0 * y[k].1[i] + z[k].1[i]) / (h * h)).sum())
                };
                let (ma, before, after) = (second(&b1, &bodies, &a1), second(&b2, &b1, &bodies), second(&bodies, &a1, &a2));
                let weight: [f64; 3] = std::array::from_fn(|i| bodies.iter().map(|(mass, _)| mass * 1e-3 * gravity[i]).sum());
                let on = |path: &str| reaction.carries.iter().any(|c| c == path || c.strip_prefix(path).is_some_and(|rest| rest.starts_with('/')));
                let ctx = ink_dynamics::Ctx { t: m.t, ..view.ctx.clone() };
                let applied = loads.iter().filter(|(path, _)| on(path)).fold([0.0; 3], |sum, (_, force)| add(sum, force.map(|n| view.doc.scalar(n, &ctx))));
                let held = m.reactions[r].force;
                let scale = size(held) + size(weight) + size(applied) + size(pull) + size(ma);
                let jump = 2e-2 * size(ma).max(size(before)).max(size(after)) + 1e-6 * scale + 1e-9;
                if !applied.iter().all(|x| x.is_finite()) || size(sub(ma, before)) > jump || size(sub(ma, after)) > jump {
                    continue;
                }
                // The joint's reaction is what the rest exerts on what it carries: with the weight, the loads
                // and the pulls it is their mass times acceleration.
                self.balanced.set(self.balanced.get() + 1);
                let left = sub(add(add(add(held, weight), applied), pull), ma);
                // (To what the solver says its read of this reaction may be off by: what spins, read round.)
                let shaken = ink_dynamics::reaction_bound(run, &reaction.path, m.t).0;
                if size(left) > 1e-3 * scale + shaken + 1e-6 {
                    note(
                        format!("b {}", reaction.path),
                        size(left),
                        format!(
                            "what {} carries weighs {} N, is loaded with {} N, is pulled by couplings with {} N and takes {} N to speed up, and its joint's reaction is {} N: {} N \
                             along {} is unaccounted at t = {}",
                            reaction.path,
                            round(size(weight)),
                            round(size(applied)),
                            round(size(pull)),
                            round(size(ma)),
                            round(size(held)),
                            round(size(left)),
                            dir(left),
                            round(m.t)
                        ),
                    );
                }
                // The same of moments, about the world's origin: the reaction's moment, the weights at their
                // centres, the loads at their points and the couplings' lines, against the rate the carried
                // bodies' angular momentum changes at (their inertias turning with them: what makes a
                // gyroscope push). By differences again, centred and from each side, which agree where it is smooth.
                let at = |t: f64| run.model.carried_momentum(r, t, h);
                let (Some(l0), Some(l1), Some(l2), Some(acting)) = (at(m.t - 2.0 * h), at(m.t), at(m.t + 2.0 * h), run.model.carried_moment(r, m, gravity)) else { continue };
                let diff = |x: [f64; 3], y: [f64; 3], over: f64| -> [f64; 3] { std::array::from_fn(|i| (y[i] - x[i]) / over) };
                let (dl, before, after) = (diff(l0, l2, 4.0 * h), diff(l0, l1, 2.0 * h), diff(l1, l2, 2.0 * h));
                let held_m = m.reactions[r].moment;
                let scale_m = size(held_m) + size(acting) + size(dl);
                let jump = 2e-2 * size(dl).max(size(before)).max(size(after)) + 1e-6 * scale_m + 1e-9;
                if size(sub(dl, before)) > jump || size(sub(dl, after)) > jump {
                    continue;
                }
                self.balanced_moments.set(self.balanced_moments.get() + 1);
                let left = sub(add(held_m, acting), dl);
                let shaken_m = ink_dynamics::reaction_bound(run, &reaction.path, m.t).1;
                if size(left) > 2e-3 * scale_m + shaken_m + 1e-6 {
                    note(
                        format!("bm {}", reaction.path),
                        size(left),
                        format!(
                            "about the world's origin, what acts on what {} carries (weights, loads, couplings) turns it with {} N*mm, its motion takes {} N*mm, and its \
                             joint's reaction gives {} N*mm: {} N*mm about {} is unaccounted at t = {}",
                            reaction.path,
                            round(size(acting)),
                            round(size(dl)),
                            round(size(held_m)),
                            round(size(left)),
                            dir(left),
                            round(m.t)
                        ),
                    );
                }
            }
        }
        worst.into_iter().map(|(_, _, said)| said).collect()
    }

    /// How many joint-moments `books` balanced apart from the solve, once it has run.
    pub fn balanced(&self) -> (usize, usize) {
        (self.balanced.get(), self.balanced_moments.get())
    }

    /// What its answers assume because the design does not say, each said once: a coupling that passes
    /// power and states no efficiency is taken as lossless (an effort through it is then the least the
    /// motion can need). Read off the run, where it has started.
    pub fn assumed(&self) -> Vec<String> {
        let Some(Ok(run)) = self.run.get() else { return vec![] };
        let mut names: Vec<&str> = vec![];
        for (_, m) in &run.at {
            for k in 0..run.model.actuators.len() {
                for name in run.model.unknown_losses(k, m) {
                    if !names.contains(&name) {
                        names.push(name);
                    }
                }
            }
        }
        // Those that give an efficiency whose loss is still not found: a mesh on a carrier that moves.
        let (carried, names): (Vec<&str>, Vec<&str>) =
            names.into_iter().partition(|name| run.model.gives_efficiency(name.rsplit_once(" (").map_or(name, |(n, _)| n)) || run.model.gives_efficiency(name));
        let on_carrier = (!carried.is_empty()).then(|| {
            format!(
                "no loss in {}: its parts ride a joint that moves (a planet on its carrier, a pulley on a carriage), and its loss, on their turning in that carrier's frame, is \
                 not found, so the efforts through it are the least the motion can need",
                carried.join(", ")
            )
        });
        if names.is_empty() {
            return on_carrier.into_iter().collect();
        }
        // Each kind's own form, with where the efficiency goes: last.
        let mut kinds: Vec<&str> = names.iter().filter_map(|name| name.rsplit_once(" (").map(|(_, kind)| kind.trim_end_matches(')'))).collect();
        kinds.dedup();
        let forms: Vec<String> = kinds.iter().map(|kind| format!("({kind} ... 0.95)")).collect();
        let unstated = format!(
            "no loss in {}: {} no efficiency, so the efforts through {} are the least the motion can need. To count its loss, write the share of power it passes last \
             in its kind, after what is there: {}, with the maker's number in place of 0.95",
            names.join(", "),
            if names.len() == 1 { "it states" } else { "they state" },
            if names.len() == 1 { "it" } else { "them" },
            forms.join(", ")
        );
        std::iter::once(unstated).chain(on_carrier).collect()
    }

    /// The imported files' loads this design does not carry though it places what they are on, each said.
    pub fn loads_left(&self) -> Vec<String> {
        ink_dynamics::loads_not_applied(&mut self.view.borrow_mut().doc, &self.target).iter().map(|l| l.said()).collect()
    }

    /// The span its run was over, if it ran over one; and, where the written motion is not over at the
    /// span's end, why: some coordinate is not at rest there, or moves again past it (as far again as the span
    /// is long, read at 8 moments), and its speed past the end does not repeat what it was over the span (a
    /// cycle; a machine running steadily repeats). An answer over the span then says nothing of what comes
    /// after: a lift read while it cruises has not yet stopped. Nothing where the run never started, or the
    /// motion past the end cannot be read.
    pub fn span(&self) -> (Option<(f64, f64)>, Option<String>) {
        // Read for its speeds alone: the same judgement, of the actuators' rates.
        if !matches!(self.run.get(), Some(Ok(_))) {
            let Some((kin, (a, b))) = self.kin.get().and_then(|k| k.as_ref()).and_then(|k| Some((k, k.span?))) else { return (None, None) };
            let view = &*self.view.borrow();
            let rates = |t: f64| -> Vec<f64> {
                let ctx = ink_dynamics::Ctx { t, ..view.ctx.clone() };
                self.rates.iter().flatten().map(|&rate| view.doc.scalar(rate, &ctx)).collect()
            };
            let end = rates(b);
            let near = |k: usize, x: f64, y: f64| (x - y).abs() <= 1e-6 * kin.rate[k].iter().fold(end[k].abs(), |m, r| m.max(r.abs())) + 1e-9;
            let moving = (0..end.len()).find(|&k| !near(k, end[k], 0.0));
            let (mut holds, mut repeats, mut first) = (moving.is_none(), true, None);
            for i in 1..=8 {
                let s = (b - a) * i as f64 / 8.0;
                let (after, again) = (rates(b + s), rates(a + s));
                for k in 0..end.len() {
                    if !near(k, after[k], end[k]) {
                        holds = false;
                        first = first.or(Some((k, b + s, after[k])));
                    }
                    if moving == Some(k) && i == 8 {
                        first = first.or(Some((k, b + s, after[k])));
                    }
                    repeats &= near(k, after[k], again[k]);
                }
            }
            let why = first.filter(|_| !holds && !repeats).map(|(k, t, rate)| {
                format!(
                    "the written motion goes on after the span's end: {} moves at {} a second at t = {b} and at {} at t = {}, neither at rest nor repeating the span. \
                     The answer is of t {a} to {b} alone: give a span that covers the motion ((t-span A B), or --t-span A:B)",
                    self.actuators[k].0,
                    (end[k] * 1e3).round() / 1e3,
                    (rate * 1e3).round() / 1e3,
                    (t * 1e3).round() / 1e3
                )
            });
            return (Some((a, b)), why.or_else(|| if holds { self.beyond(view, a, b) } else { None }));
        }
        let Some(Ok(run)) = self.run.get() else { return (None, None) };
        let Some((a, b)) = run.span else { return (None, None) };
        let rates = |t: f64| {
            let m = run.model.at(t);
            m.unsolved.is_none().then(|| m.q.iter().map(|q| q[1]).collect::<Vec<f64>>())
        };
        let Some(end) = rates(b) else { return (run.span, None) };
        let most = |k: usize| run.at.iter().map(|(_, m)| m.q[k][1].abs()).fold(end[k].abs(), f64::max);
        let near = |k: usize, x: f64, y: f64| (x - y).abs() <= 1e-6 * most(k) + 1e-9;
        // At rest at the end: then what moves it again after is the first thing said, else that it still moves.
        let moving = (0..end.len()).find(|&k| !near(k, end[k], 0.0));
        let (mut holds, mut repeats, mut first) = (moving.is_none(), true, None);
        for i in 1..=8 {
            let s = (b - a) * i as f64 / 8.0;
            let Some(after) = rates(b + s) else { return (run.span, None) };
            let again = rates(a + s);
            for k in 0..end.len() {
                if !near(k, after[k], end[k]) {
                    holds = false;
                    first = first.or(Some((k, b + s, after[k])));
                }
                if moving == Some(k) && i == 8 {
                    first = first.or(Some((k, b + s, after[k])));
                }
                repeats &= again.as_ref().is_some_and(|again| near(k, after[k], again[k]));
            }
        }
        let why = first.filter(|_| !holds && !repeats).map(|(k, t, rate)| {
            format!(
                "the written motion goes on after the span's end: {} moves at {} a second at t = {b} and at {} at t = {}, neither at rest nor repeating the span. \
                 The answer is of t {a} to {b} alone: give a span that covers the motion ((t-span A B), or --t-span A:B)",
                run.model.coords[k].name,
                (end[k] * 1e3).round() / 1e3,
                (rate * 1e3).round() / 1e3,
                (t * 1e3).round() / 1e3
            )
        });
        (run.span, why.or_else(|| if holds { self.beyond(&self.view.borrow(), a, b) } else { None }))
    }

    /// A motion at rest at the span's end and through as long again: whether it is over, by its formulas.
    /// None: every joint's drive is proven at rest at every t after b (the kernel's bounds of its rate over
    /// all of it: a profile with a duration, a clamp or a select on t, has run out). Else where it moves
    /// again, found by halving what is not proven at rest and reading the rate there; or, where nothing
    /// is found moving and nothing proves it at rest, that it is not known: either way the answer is of the
    /// span alone.
    fn beyond(&self, view: &ModelView, a: f64, b: f64) -> Option<String> {
        if !view.ctx.p_spans.is_empty() {
            return None;
        }
        let ctx = &view.ctx;
        let rest = |k: usize, from: f64, to: f64| match self.moved[k].0 {
            Some(rate) => {
                let (low, high) = ink_dynamics::rate_bounds(&view.doc, rate, ctx, (from, to));
                low.abs() <= 1e-9 && high.abs() <= 1e-9
            }
            None => false,
        };
        // For as long as anything is written to last: a thousand million seconds.
        // (From just past the end: a profile that ends at the span's end is read there as either piece.)
        let (b0, far) = (b + 1e-6 * (b - a), b + 1e9);
        let open: Vec<usize> = (0..self.moved.len()).filter(|&k| !rest(k, b0, far)).collect();
        let k = *open.first()?;
        let alone = format!("The answer is of t {a} to {b} alone: give a span that covers the motion ((t-span A B), or --t-span A:B)");
        let Some(rate) = self.moved[k].0 else {
            return Some(format!(
                "it is not known whether the written motion goes on after the span's end: {} is moved by a value solved or simulated, which is not read past t = {b}. {alone}",
                self.joints[k].0
            ));
        };
        // The earliest stretch not proven at rest, halved down to a thousandth of the span.
        let (mut from, mut to) = (b0, far);
        while to - from > 1e-3 * (b - a) {
            let mid = 0.5 * (from + to);
            if rest(k, from, mid) {
                from = mid;
            } else {
                to = mid;
            }
        }
        let read = |t: f64| view.doc.scalar(rate, &ink_dynamics::Ctx { t, ..ctx.clone() });
        let found = (0..=64).map(|i| from + (to - from + 0.25 * (b - a)) * i as f64 / 64.0).map(|t| (t, read(t))).find(|(_, r)| r.abs() > 1e-9);
        Some(match found {
            Some((t, r)) => format!(
                "the written motion goes on after the span's end: {} is at rest at t = {b} and moves again, at {} a second at t = {}. {alone}",
                self.joints[k].0,
                (r * 1e3).round() / 1e3,
                (t * 1e3).round() / 1e3
            ),
            None => format!(
                "it is not known whether the written motion goes on after the span's end: nothing proves {} at rest for good after t = {b} (its drive's formula is not \
                 one whose end the kernel can read: a profile with a duration, a clamp or a select on t is), and it was not found moving. {alone}",
                self.joints[k].0
            ),
        })
    }

    /// How many times the solver has been started.
    pub fn runs(&self) -> usize {
        self.runs.get()
    }

    /// The tool for the run's moment `moment` (0, 1, ...), None past its
    /// last. Asking for one starts nothing: before the run is made the
    /// moments are the written motion's own (`Kin`), at which a speed is
    /// read; an effort asked at one starts the run, and the moments after
    /// are the run's (the first of either is the span's start).
    pub fn at(&'a self, moment: usize) -> Option<AtMoment<'a>> {
        // (A run that did not start has no moments: the motion's own stand, and its speeds with them.)
        let known = match self.run.get() {
            Some(Ok(run)) => run.at.len(),
            _ => self.kin().map_or(1, |kin| kin.t.len()),
        };
        (moment < known).then_some(AtMoment { effort: self, moment, early: !matches!(self.run.get(), Some(Ok(_))) })
    }

    /// The written motion read for its speeds (`Kin`): None where an actuator's rate is the run's to give
    /// (it reads a simulated coordinate or a solve), where a param is ranged, or where a rate is no number.
    fn kin(&self) -> Option<&Kin> {
        self.kin
            .get_or_init(|| {
                let view = &*self.view.borrow();
                let rates: Vec<ink_dynamics::Rate> = self.rates.iter().copied().collect::<Option<_>>()?;
                if !view.ctx.p_spans.is_empty() {
                    return None;
                }
                let span = view.motion.or(view.doc.t_span);
                let t: Vec<f64> = match span {
                    Some((a, b)) => {
                        let h = (b - a) / KIN_STEPS as f64;
                        std::iter::once(a).chain((0..KIN_STEPS).map(|i| a + h * (i as f64 + 0.5))).chain([b]).collect()
                    }
                    None => vec![view.ctx.t],
                };
                let at = |t: f64| -> Vec<f64> {
                    let ctx = ink_dynamics::Ctx { t, ..view.ctx.clone() };
                    rates.iter().map(|&rate| view.doc.scalar(rate, &ctx)).collect()
                };
                let by_moment: Vec<Vec<f64>> = t.iter().map(|&t| at(t)).collect();
                let rate: Vec<Vec<f64>> = (0..rates.len()).map(|k| by_moment.iter().map(|m| m[k]).collect()).collect();
                rate.iter().flatten().all(|x| x.is_finite()).then_some(Kin { t, rate, span })
            })
            .as_ref()
    }

    /// Why a tool made before the run, for the written motion's moment `moment`, has no effort or load of the
    /// run's moment of that number (None: it has): the two are one moment only at the span's start. So no
    /// speed of one moment is ever read beside an effort of another, and no stretch between the two lists
    /// goes unread; whoever asks moment by moment from the first (`over`) starts the run there and never
    /// meets this.
    fn other_moment(&self, run: &Read, moment: usize) -> Option<String> {
        let kin = self.kin.get()?.as_ref()?;
        let same = moment == 0 && run.at.first().is_some_and(|(_, m)| m.t == kin.t[0]);
        (!same).then(|| {
            format!(
                "this moment (t = {}) was read from the written motion alone, before any effort was asked; the efforts' run reads moments of its own. \
                 Ask the motion's moments in order from the first",
                kin.t.get(moment).copied().unwrap_or(f64::NAN)
            )
        })
    }

    /// Actuator `k`'s speed at the motion's own moment `moment`, where no run has started: its rate there,
    /// with over a motion every speed from that moment's to its neighbours'.
    fn kin_speed(&self, about: &About, k: usize, moment: usize, unit: &str) -> Option<Fact> {
        if matches!(self.run.get(), Some(Ok(_))) {
            return None;
        }
        let kin = self.kin()?;
        let (rate, t) = (*kin.rate[k].get(moment)?, kin.t[moment]);
        let mut fact = Fact { about: about.clone(), unit: unit.into(), from: format!("the written motion of {}", self.target), ..Fact::default() };
        (fact.value, fact.strength, fact.covers) = (Some(rate), Some(Strength::Estimated), vec![("t".into(), Covers::One)]);
        fact.notes.push(format!("at t = {t}"));
        fact.notes.push(match kin.span {
            Some(_) => format!("the rate of its drive as written, read at the span's ends and the middles of {KIN_STEPS} even steps: an estimate between them"),
            None => "the rate of its drive as written, at that moment".into(),
        });
        let near = [moment.checked_sub(1), Some(moment + 1)].into_iter().flatten().filter_map(|i| kin.rate[k].get(i));
        let (low, high) = near.fold((rate, rate), |(low, high), &r| (low.min(r), high.max(r)));
        if low < high {
            fact.band = Some((low, high));
            fact.notes.push("its band reaches to its speed at the moments before and after".into());
        }
        Some(fact)
    }

    fn run(&self) -> &Result<Read, String> {
        self.run.get_or_init(|| {
            self.runs.set(self.runs.get() + 1);
            let view = &mut *self.view.borrow_mut();
            // Over the span the design is asked over (a command's --t-span), else the file's; with neither, at
            // the view's moment.
            let span = view.motion.or(view.doc.t_span);
            // The motion whatever the masses and the gravity: a speed needs neither, and an effort or a
            // reaction is given only where both are there for it (`lacks`).
            ink_dynamics::read_kinematic(view, &self.target, span).map_err(|e| e.message)
        })
    }

    /// What is asked of the actuator at `on` (its path, or with `frames` a
    /// frame its part has): which it is, and as what unit an effort and a
    /// rate are.
    fn asked(&self, on: &str, frames: bool) -> Option<(&str, &'static str, &'static str, &'static str)> {
        let frame = |path: &str, has: &[String]| frames && on.strip_prefix(path).and_then(|rest| rest.strip_prefix(' ')).is_some_and(|frame| has.iter().any(|f| f == frame));
        let (path, turn, _) = self.actuators.iter().find(|(path, _, has)| path == on || frame(path, has))?;
        // The solver's own units, and the language's: N mm and deg/s on a turn, N and mm/s on a slide.
        Some(if *turn { (path, "torque", "N*mm", "deg/s") } else { (path, "force", "N", "mm/s") })
    }

    /// The fact about `about`, in `unit`, before its number: the run, the
    /// actuator's place in it and what every fact of it says, or the fact
    /// open with why there is none.
    /// Why the run, read for its motion alone, has no effort of actuator `k` (None: it has): no gravity is
    /// declared, or moving bodies its effort rests on have no mass, each named.
    fn lacks(run: &Read, k: usize) -> Option<String> {
        if run.weightless {
            return Some(ink_dynamics::no_gravity().message);
        }
        let without = run.model.blocked(k);
        (!without.is_empty()).then(|| format!("these moving bodies its effort rests on have no mass: {}", run.model.massless_said(&without)))
    }

    /// `effort`: the fact is an effort (a torque, a force, a peak), not the motion's speed.
    fn begun(&self, about: &About, path: &str, unit: &str, effort: bool) -> Result<(&Read, usize, Fact), Fact> {
        let mut fact = Fact { about: about.clone(), unit: unit.into(), from: format!("ink effort on {}", self.target), ..Fact::default() };
        let open = |mut fact: Fact, why, note: String| {
            (fact.held, fact.notes) = (Some(why), vec![note]);
            fact
        };
        let (outside, undecided) = (Verdict::NotAvailable(Unavailable::Outside), Verdict::NotSettled(Unsettled::Undecided));
        // An effort with no gravity declared is not found: said before anything is built for it.
        if effort && ink_dynamics::gravity(&self.view.borrow().doc).is_none() {
            return Err(open(fact, Verdict::NotSettled(Unsettled::Lacks), ink_dynamics::no_gravity().message));
        }
        let run = match self.run() {
            Ok(run) => run,
            Err(why) if ink_dynamics::gravity(&self.view.borrow().doc).is_none() => return Err(open(fact, Verdict::NotSettled(Unsettled::Lacks), why.clone())),
            // Why, the solver says only in words (a moving body with no mass): a source that could not decide.
            Err(why) => return Err(open(fact, undecided, why.clone())),
        };
        let Some(k) = run.model.actuators.iter().position(|a| a.path == path) else {
            return Err(open(fact, outside, "the run names no such actuator".into()));
        };
        if let Some(why) = self.ranged(run) {
            return Err(open(fact, outside, why));
        }
        if let Some(why) = Self::lacks(run, k).filter(|_| effort) {
            return Err(open(fact, Verdict::NotSettled(Unsettled::Lacks), why));
        }
        match &run.sprung[k] {
            Err(why) => return Err(open(fact, undecided, format!("refused: {why}"))),
            Ok(true) => return Err(open(fact, outside, "its effort is its spring's in a simulation, whose moments this plug does not read".into())),
            Ok(false) => {}
        }
        fact.strength = Some(Strength::Estimated);
        // What the run leaves out of this actuator's facts: the same at every moment, so written once.
        let mut left = self.left.borrow_mut();
        let left = left.entry(k).or_insert_with(|| {
            let mut left = vec![];
            if let Some(steps) = &run.steps {
                // The steps the model did not read or read loosely, and this actuator's own: where its rating
                // may cross its effort inside one. Said as every answer from moments says it.
                let mut gaps: Vec<(f64, f64)> = steps.unchecked.iter().chain(&steps.loose).chain(&run.gaps[k]).copied().collect();
                gaps.sort_by(|a, b| a.0.total_cmp(&b.0));
                left.extend(ink_lang::motion::not_followed(ink_lang::motion::MODEL_GAP, &gaps));
            }
            left.extend(run.model.unsettled.iter().chain(&run.model.unsettled_carried).map(|m| format!("a mass that did not settle: {m}")));
            left
        });
        fact.left_out.extend(left.iter().cloned());
        fact.notes.push(run.method.into());
        Ok((run, k, fact))
    }
}

impl Effort<'_> {
    /// Why the run says nothing of the design as asked, if it does not: it is at one value of each param, and
    /// of a ranged one its model reads that value says nothing of the rest.
    fn ranged(&self, run: &Read) -> Option<String> {
        let view = self.view.borrow();
        let ranged: Vec<String> = (run.model.reads.iter().chain(&run.model.reads_carried))
            .filter_map(|name| Some((name, view.ctx.p_spans.get(view.doc.param_index(name)? as u32)?)))
            .map(|(name, (low, high))| format!("'{name}' ({low} to {high})"))
            .collect();
        (!ranged.is_empty()).then(|| {
            format!(
                "it is found at one value of each param, and its model reads ranged {}: that value says nothing of the rest of the range. Ask at one value (--param NAME=V)",
                ranged.join(", ")
            )
        })
    }

    /// The load the model works out at a place at the run's moment `moment`, where it does: on a support of a
    /// reaction the design divides (`(reaction PATH (by S ...))`), `radial-load` and `axial-load` (N: its
    /// share's force across and along its load frame's z) and, where PATH is an actuator that turns, `speed`
    /// (deg/s, PATH's own); on the part a joint moves (its path, or a frame the part has), `radial-load` and
    /// `axial-load`: the joint's reaction across and along its axis, and nothing where the design divides that
    /// reaction among supports, which then carry all of it. Each a size, as it is at that moment.
    /// And on the part a joint moves, `moment-x`, `moment-y` and `moment-z` (N*mm): the size of the joint's
    /// reaction moment about each axis of that instance's frame (a carriage's roll, pitch and yaw, by its
    /// part's frame).
    fn load(&self, about: &About, moment: usize, early: bool, claims: &Claims) -> Option<Fact> {
        let what = ["radial-load", "axial-load", "speed", "moment-x", "moment-y", "moment-z"].iter().position(|w| *w == about.what)?;
        let within = |path: &str| about.on == path || about.on.strip_prefix(path).is_some_and(|rest| rest.starts_with(['/', ' ']));
        let support = self.divided.iter().enumerate().find_map(|(d, division)| Some((d, division.supports.iter().position(|s| within(s))?)));
        let joint = self.joints.iter().find(|(path, frames, _)| {
            about.on == *path || about.on.strip_prefix(path).and_then(|rest| rest.strip_prefix(' ')).is_some_and(|frame| frames.iter().any(|f| f == frame))
        });
        // A speed: of the joint its part is on, or of the joint it supports. An actuator's own is its tool's
        // (`AtMoment::ask`); any other joint's is its drive's rate as written, which the couplings and the
        // design's formulas fix (an idler a belt turns, a carriage it pulls).
        if what == 2 {
            let path = match (support, joint) {
                (Some((d, _)), _) => &self.divided[d].path,
                (None, Some((path, ..))) if !self.actuators.iter().any(|(p, ..)| p == path) => path,
                _ => return None,
            };
            let j = self.joints.iter().position(|(p, ..)| p == path)?;
            let said = support.map(|(d, _)| format!("the speed of {path}, which it carries: (supports {path} {})", self.divided[d].supports.join(" ")));
            return self.joint_speed(about, j, moment, early, said);
        }
        match (support, joint, what) {
            // A moment about an axis: of a joint's own instance, or of a guide's carriage among the supports.
            (Some((d, k)), _, 3..) if !self.divided[d].guides[k] => return None,
            (None, None, 3..) => return None,
            // A load on a part that is no joint's and that nothing names as holding one: said, with what to write.
            (None, None, _) if about.on.starts_with(&format!("{}/", self.target)) => {
                let place = about.on.split(' ').next().unwrap_or(&about.on);
                // Where it sits on a turning joint's axis, the line itself, with the other bearings there.
                let line = ink_dynamics::sits_on(&mut self.view.borrow_mut(), &self.target, place).map(|(joint, others)| {
                    let all: Vec<&str> = std::iter::once(place).chain(others.iter().map(String::as_str)).collect();
                    format!(
                        " It sits on the axis of {joint}: if it holds what turns there, the line is (supports {joint} {}), at the top level, or inside their assembly with its \
                         own instance names",
                        all.join(" ")
                    )
                });
                return Some(Fact {
                    about: about.clone(),
                    held: Some(Verdict::NotAvailable(Unavailable::NoSource)),
                    notes: vec![format!(
                        "nothing states it, and the model works out no load on {place}: it is on no joint (a part on a joint is held by that joint, and its load is the \
                         joint's reaction, with nothing more to write), and nothing says it holds what is on one. If it holds a body that turns or slides on another \
                         instance's joint (a bearing a shaft turns in), say so beside that joint: (supports JOINT {place} ...), JOINT the instance that turns or slides, \
                         never {place} itself; or state the load.{}",
                        line.unwrap_or_default()
                    )],
                    ..Fact::default()
                });
            }
            (None, None, _) => return None,
            _ => {}
        }
        let unit = ["N", "N", "deg/s", "N*mm", "N*mm", "N*mm"][what];
        let mut fact = Fact { about: about.clone(), unit: unit.into(), from: format!("ink effort on {}", self.target), ..Fact::default() };
        let open = |mut fact: Fact, why, note: String| {
            (fact.held, fact.notes) = (Some(why), vec![note]);
            Some(fact)
        };
        let (outside, undecided) = (Verdict::NotAvailable(Unavailable::Outside), Verdict::NotSettled(Unsettled::Undecided));
        // A load with no gravity declared is not found: said before anything is built for it.
        if what != 2 && ink_dynamics::gravity(&self.view.borrow().doc).is_none() {
            return open(fact, Verdict::NotSettled(Unsettled::Lacks), format!("nothing states it, and the model does not work it out: {}", ink_dynamics::no_gravity().message));
        }
        let run = match self.run() {
            Ok(run) => run,
            Err(why) if ink_dynamics::gravity(&self.view.borrow().doc).is_none() => {
                return open(fact, Verdict::NotSettled(Unsettled::Lacks), format!("nothing states it, and the model does not work it out: {why}"));
            }
            Err(why) => return open(fact, undecided, why.clone()),
        };
        if let Some(why) = self.ranged(run) {
            return open(fact, outside, why);
        }
        // A load is of the masses and the gravity: with none declared, or with a moving body that has no
        // mass among what its reaction rests on, it is not found. (A support's speed is the motion's.)
        if what != 2 {
            let lacks = Verdict::NotSettled(Unsettled::Lacks);
            if run.weightless {
                return open(fact, lacks, format!("nothing states it, and the model does not work it out: {}", ink_dynamics::no_gravity().message));
            }
            let path = support.map(|(d, _)| &self.divided[d].path).or(joint.map(|(path, ..)| path));
            let without = path.and_then(|p| run.model.reactions.iter().position(|r| r.path == *p)).map(|r| run.model.blocks_reaction(r)).unwrap_or_default();
            if !without.is_empty() {
                return open(fact, lacks, format!("these moving bodies its reaction rests on have no mass: {}", run.model.massless_said(&without)));
            }
        }
        if let Some(why) = self.other_moment(run, moment).filter(|_| early) {
            return open(fact, undecided, why);
        }
        let m = &run.at.get(moment)?.1;
        if let Some(why) = &m.unsolved {
            return open(fact, undecided, format!("refused at t {}: {why}", m.t));
        }
        fact.notes.push(format!("at t = {}", m.t));
        fact.left_out.extend(run.model.unsettled.iter().chain(&run.model.unsettled_carried).map(|m| format!("a mass that did not settle: {m}")));
        if run.span.is_some() {
            fact.left_out.push("a reaction's peak between the moments read: they are chosen for the efforts".into());
        }
        (fact.strength, fact.covers) = (Some(Strength::Estimated), vec![("t".into(), Covers::One)]);
        // A share's size: across the support's z and along it; a guide's carriage across its x, its travel
        // (it carries none along it), and its moment about each of its axes.
        let guide = support.is_some_and(|(d, k)| self.divided[d].guides[k]);
        let size = |share: &[f64]| match (guide, what) {
            (_, 3..) => share[what].abs(),
            (true, 0) => share[1].hypot(share[2]),
            (true, _) => 0.0,
            (false, 0) => share[0].hypot(share[1]),
            (false, _) => share[2].abs(),
        };
        // A load read off a reaction is as good as the masses under it: each body's measured moments
        // anywhere in their proven bounds move it (ink_dynamics::reaction_moves, Divided::share_moves), and
        // its band says how far, either way.
        // And as good as the solver's own read of it: over a span, bodies that spin are read round, and
        // their departure from round may move the reaction by its `bound` at some angle of the spin, either
        // way. A share of it may be off by as much.
        let held = support.map(|(d, _)| &self.divided[d].path).or(joint.map(|(path, ..)| path));
        let (shaken_f, shaken_m) = held.map_or((0.0, 0.0), |path| ink_dynamics::reaction_bound(run, path, m.t));
        let shaken = if what >= 3 { shaken_m } else { shaken_f };
        let within = |mut fact: Fact, load: f64, masses: f64| {
            fact.value = Some(load);
            if masses > 0.0 || shaken > 0.0 {
                fact.band = Some(((load - masses - shaken).max(0.0), load + masses + shaken));
            }
            if masses > 0.0 {
                fact.notes.push(format!(
                    "within {} {} either way: what the masses it rests on can move it, each within its proven bounds",
                    (masses * 1e3).ceil() / 1e3,
                    if what >= 3 { "N*mm" } else { "N" }
                ));
            }
            if shaken > 0.0 {
                fact.notes.push(format!(
                    "and within {} {} more either way: over the motion what spins is read round, at rest, and that is the most its departure from round moves this at any angle of its spin",
                    (shaken * 1e3).ceil() / 1e3,
                    if what >= 3 { "N*mm" } else { "N" }
                ));
            }
            fact
        };
        if let Some((d, k)) = support {
            let division = &self.divided[d];
            let said = format!("(supports {} {})", division.path, division.supports.join(" "));
            let shares = self
                .shares
                .borrow_mut()
                .entry((d, moment))
                .or_insert_with(|| division.shares(run, m.t, &self.view.borrow().ctx.params).unwrap_or_else(|| Err(format!("the run reads no reaction of {}", division.path))))
                .clone();
            // A part named as holding several joints (two shafts in one pair of bearings) carries its share
            // of each: they are summed, each in the part's own load frame.
            let mut shares = shares;
            let mut said = said;
            let mut several = false;
            for (d2, division2) in self.divided.iter().enumerate().filter(|(d2, _)| *d2 != d) {
                let under = |path: &str| about.on == path || about.on.strip_prefix(path).is_some_and(|rest| rest.starts_with(['/', ' ']));
                let Some(k2) = division2.supports.iter().position(|s| under(s)) else { continue };
                let also = self
                    .shares
                    .borrow_mut()
                    .entry((d2, moment))
                    .or_insert_with(|| {
                        division2.shares(run, m.t, &self.view.borrow().ctx.params).unwrap_or_else(|| Err(format!("the run reads no reaction of {}", division2.path)))
                    })
                    .clone();
                shares = match (shares, also) {
                    (Ok((mut shares, thrust)), Ok((also, thrust2))) => {
                        (0..6).for_each(|i| shares[k][i] += also[k2][i]);
                        let thrust = match (thrust, thrust2.filter(|_| division2.locates[k2])) {
                            (a, None) => a,
                            (a, Some(b)) => Some(a.unwrap_or(0.0) + b),
                        };
                        Ok((shares, thrust))
                    }
                    (Err(why), _) | (_, Err(why)) => Err(why),
                };
                said = format!("{said} and (supports {} {})", division2.path, division2.supports.join(" "));
                several = true;
            }
            // What the masses can move this support's share by: each change of the reaction(s) split as the share is.
            let params = self.view.borrow().ctx.params.clone();
            let mut moves: Vec<(Vec<[f64; 6]>, f64)> = division.share_moves(run, m.t, &params);
            for (_, division2) in self.divided.iter().enumerate().filter(|(d2, _)| *d2 != d) {
                let under = |path: &str| about.on == path || about.on.strip_prefix(path).is_some_and(|rest| rest.starts_with(['/', ' ']));
                let Some(k2) = division2.supports.iter().position(|s| under(s)) else { continue };
                moves.extend(division2.share_moves(run, m.t, &params).into_iter().map(|(d, t)| {
                    let mut only = vec![[0.0; 6]; d.len().max(k + 1)];
                    only[k] = d[k2];
                    (only, if division2.locates[k2] { t } else { 0.0 })
                }));
            }
            return match shares {
                Ok((shares, thrust)) => {
                    let masses: f64 =
                        moves.iter().filter_map(|(d, _)| d.get(k)).map(|d| (size(&std::array::from_fn::<f64, 6, _>(|i| shares[k][i] + d[i])) - size(&shares[k])).abs()).sum();
                    let thrust_masses: f64 = moves.iter().map(|(_, t)| t.abs()).sum();
                    fact.notes.push(match several {
                        true => format!("the sum of its shares of what it is named as holding, by statics: {said}"),
                        false => format!("its share of what {} carries, by statics: {said}", division.path),
                    });
                    match thrust.filter(|total| what == 1 && *total > 1e-9 * (1.0 + size(&shares[k]))) {
                        None => Some(within(fact, size(&shares[k]), masses)),
                        // The thrust the locating bearings carry among them: statics does not say which of them.
                        Some(_) if !division.locates[k] => Some(within(fact, 0.0, 0.0)),
                        Some(total) => {
                            // What the design states of the others' thrust leaves this one: the total less the most
                            // they are said to carry, to the total less the least.
                            let (mut least, mut most) = (total, total);
                            for (other, _) in division.supports.iter().enumerate().filter(|(j, _)| *j != k && division.locates[*j]) {
                                let on = about.on.replacen(division.supports[k].as_str(), division.supports[other].as_str(), 1);
                                let said: Vec<f64> =
                                    claims.stated(&About { what: about.what.clone(), on }).iter().filter(|f| f.unit == "N").filter_map(|f| f.value).map(f64::abs).collect();
                                match said.as_slice() {
                                    [] => least = 0.0,
                                    said => {
                                        least -= said.iter().copied().fold(0.0, f64::max);
                                        most -= said.iter().copied().fold(f64::INFINITY, f64::min);
                                    }
                                }
                            }
                            let (least, most) = (least.max(0.0), most.max(0.0));
                            // Held beside what the design states of this one, statics says only the least it must
                            // carry (what the others are not said to): which bearing locates the shaft is the
                            // design's to say, and saying this one carries none is no claim against the model.
                            let most = if claims.weighing() { least } else { most };
                            fact.value = Some(most);
                            fact.band = (least < most || thrust_masses > 0.0).then_some(((least - thrust_masses).max(0.0), most + thrust_masses));
                            let locating: Vec<&str> = division.supports.iter().zip(&division.locates).filter(|(_, l)| **l).map(|(s, _)| s.as_str()).collect();
                            fact.notes.push(format!(
                                "the bearings that can locate the shaft ({}) carry {} N along it between them, and statics does not say which: the design does, by \
                                 stating each one's ((fact axial-load ...): the whole of it on the one that locates the shaft, 0N on one free to float)",
                                locating.join(", "),
                                (total * 1e3).round() / 1e3
                            ));
                            Some(fact)
                        }
                    }
                }
                Err(why) => open(Fact { strength: None, covers: vec![], left_out: vec![], ..fact }, undecided, format!("refused at t {}: {why}", m.t)),
            };
        }
        let (path, _, axis) = joint?;
        let axis = *axis;
        if let Some(division) = self.divided.iter().find(|d| d.path == *path) {
            fact.value = Some(0.0);
            fact.notes.push(format!("what {path} carries is held by its supports, which the design names: (supports {path} {})", division.supports.join(" ")));
            return Some(fact);
        }
        let index = run.model.reactions.iter().position(|r| r.path == *path)?;
        let closed = Fact { strength: None, covers: vec![], left_out: vec![], ..fact.clone() };
        // What a coupling pulls with on what the joint carries (a belt on its pulley, clamped to its
        // carriage) is the joint's to hold with the rest: the model's reaction has it.
        if run.model.coupled(index) {
            fact.notes.push(format!("with what the couplings on what {path} carries pull with (a belt, a mesh, a screw): the joint holds that too"));
        }
        let reaction =
            self.reactions.borrow_mut().entry((path.clone(), moment)).or_insert_with(|| ink_dynamics::reaction(run, path, m.t).expect("its reaction is in the run")).clone();
        // What the masses can move the reaction by, each change in its path's frame (as `reaction` reads it).
        let moves = ink_dynamics::reaction_moves(run, path, m.t, None);
        match reaction {
            Ok((_, moment)) if what >= 3 => {
                let axis = ["x", "y", "z"][what - 3];
                fact.notes.push(format!("the model's reaction at {path}, its moment about {path}'s own {axis}: (reaction {path})"));
                let masses = moves.iter().map(|(_, dm)| dm[what - 3].abs()).sum();
                Some(within(fact, moment[what - 3].abs(), masses))
            }
            Ok((force, moment)) => {
                let across_of = |f: [f64; 3]| (0..3).filter(|k| *k != axis).map(|k| f[k]).fold(0.0, f64::hypot);
                let (along, across) = (force[axis].abs(), across_of(force));
                fact.notes.push(format!("the model's reaction at {path}, {} its axis: (reaction {path})", if what == 0 { "across" } else { "along" }));
                let masses: f64 = moves
                    .iter()
                    .map(|(df, _)| match what {
                        0 => (across_of(std::array::from_fn(|k| force[k] + df[k])) - across).abs(),
                        _ => df[axis].abs(),
                    })
                    .sum();
                let mut fact = within(fact, if what == 0 { across } else { along }, masses);
                let tilt = (0..3).filter(|k| *k != axis).map(|k| moment[k]).fold(0.0, f64::hypot);
                if what == 0 && tilt > 1e-6 * (1.0 + across) {
                    // Said once whatever the moment: its size there is a note of that moment's.
                    fact.left_out.push("the moment across its axis: where along the shaft the load acts is not held against the rating's point".into());
                    fact.notes.push(format!("the moment across its axis there is {} N*mm about the joint's frame", (tilt * 1e3).round() / 1e3));
                }
                Some(fact)
            }
            Err(why) => open(closed, undecided, format!("refused at t {}: {why}", m.t)),
        }
    }
}

impl Effort<'_> {
    /// `tension` (N) on a coupling's part, an end of one or a part within an end, at the run's moment
    /// `moment`: the largest of its couplings' lines' tensions there (a screw's thrust on it and its nut, a
    /// belt's tight span on its pulleys, a tooth's force on meshing gears), as the model's balances give them.
    fn tension(&self, about: &About, moment: usize, early: bool) -> Option<Fact> {
        let within = |end: &str| about.on == end || about.on.strip_prefix(end).is_some_and(|rest| rest.starts_with(['/', ' ']));
        if about.what != "tension" || !self.ends.iter().any(|end| within(end)) {
            return None;
        }
        let fact = Fact { about: about.clone(), unit: "N".into(), from: format!("ink effort on {}", self.target), ..Fact::default() };
        let open = |mut fact: Fact, why, note: String| {
            (fact.held, fact.notes) = (Some(why), vec![note]);
            Some(fact)
        };
        let (lacks, outside, undecided) = (Verdict::NotSettled(Unsettled::Lacks), Verdict::NotAvailable(Unavailable::Outside), Verdict::NotSettled(Unsettled::Undecided));
        if ink_dynamics::gravity(&self.view.borrow().doc).is_none() {
            return open(fact, lacks, format!("nothing states it, and the model does not work it out: {}", ink_dynamics::no_gravity().message));
        }
        let run = match self.run() {
            Ok(run) => run,
            Err(why) => return open(fact, undecided, why.clone()),
        };
        if let Some(why) = self.ranged(run) {
            return open(fact, outside, why);
        }
        if let Some(why) = self.other_moment(run, moment).filter(|_| early) {
            return open(fact, undecided, why);
        }
        let m = &run.at.get(moment)?.1;
        if let Some(why) = &m.unsolved {
            return open(fact, undecided, format!("refused at t {}: {why}", m.t));
        }
        let mut fact = fact;
        match run.model.pull(&about.on, m)? {
            Err(why) if why.contains("have no mass") => open(fact, lacks, why),
            Err(why) => open(fact, undecided, format!("refused at t {}: {why}", m.t)),
            Ok((tension, coupling)) => {
                (fact.value, fact.strength, fact.covers) = (Some(tension), Some(Strength::Estimated), vec![("t".into(), Covers::One)]);
                fact.notes.push(format!("at t = {}", m.t));
                fact.notes.push(format!("the tightest line of {coupling}, as the model's balances give it: (tension ...)"));
                let masses = run.model.pull_mass(&about.on, m);
                if masses > 0.0 {
                    fact.band = Some(((tension - masses).max(0.0), tension + masses));
                    fact.notes.push(format!("within {} N either way: what the masses it rests on can move it, each within its proven bounds", (masses * 1e3).ceil() / 1e3));
                }
                fact.left_out.extend(run.model.unsettled.iter().chain(&run.model.unsettled_carried).map(|m| format!("a mass that did not settle: {m}")));
                if !run.model.unseparated.is_empty() {
                    fact.left_out.push(format!("what these exert, moving with its joints through nothing declared: {}", run.model.unseparated.join(", ")));
                }
                if run.span.is_some() {
                    fact.left_out.push("a tension's peak between the moments read: they are chosen for the efforts".into());
                }
                Some(fact)
            }
        }
    }
}

impl Effort<'_> {
    /// `speed` at moment `moment` of joint `j`, one that is no actuator: the rate of its drive as written
    /// (deg/s on a turn, mm/s on a slide), with as its band every speed from there to its neighbours'. The
    /// written motion's: it needs no mass and no gravity. `said`: whose speed it is, where the fact is of a
    /// part that supports the joint's.
    fn joint_speed(&self, about: &About, j: usize, moment: usize, early: bool, said: Option<String>) -> Option<Fact> {
        let (rate, turn) = self.moved[j];
        let path = &self.joints[j].0;
        let mut fact = Fact { about: about.clone(), unit: if turn { "deg/s" } else { "mm/s" }.into(), from: format!("the written motion of {}", self.target), ..Fact::default() };
        let open = |mut fact: Fact, note: String| {
            (fact.held, fact.notes) = (Some(Verdict::NotSettled(Unsettled::Undecided)), vec![note]);
            Some(fact)
        };
        let Some(rate) = rate else {
            return open(fact, format!("the speed of {path} is not a formula of t alone: its drive reads a value solved or simulated at a run"));
        };
        // The moments' times: the run's where it has started, else the written motion's own.
        let times: Vec<f64> = match self.run.get() {
            Some(Ok(run)) => {
                if let Some(why) = self.other_moment(run, moment).filter(|_| early) {
                    return open(fact, why);
                }
                run.at.iter().map(|(_, m)| m.t).collect()
            }
            _ => match self.kin() {
                Some(kin) => kin.t.clone(),
                None => self.run().as_ref().ok()?.at.iter().map(|(_, m)| m.t).collect(),
            },
        };
        let view = &*self.view.borrow();
        if !view.ctx.p_spans.is_empty() {
            return None;
        }
        let at = |i: usize| times.get(i).map(|&t| view.doc.scalar(rate, &ink_dynamics::Ctx { t, ..view.ctx.clone() }));
        let speed = at(moment).filter(|x| x.is_finite())?;
        (fact.value, fact.strength, fact.covers) = (Some(speed), Some(Strength::Estimated), vec![("t".into(), Covers::One)]);
        fact.notes.push(format!("at t = {}", times[moment]));
        fact.notes.push(said.unwrap_or_else(|| format!("the rate of {path}'s drive as written")));
        let near = [moment.checked_sub(1), Some(moment + 1)].into_iter().flatten().filter_map(at);
        let (low, high) = near.fold((speed, speed), |(low, high), r| (low.min(r), high.max(r)));
        if low < high {
            fact.band = Some((low, high));
            fact.notes.push("its band reaches to its speed at the moments before and after".into());
        }
        Some(fact)
    }
}

/// Why actuator `k`'s effort at moment `m` is not the motion's (the solver's
/// own claim is refused there), if it is not.
fn refused(run: &Read, k: usize, m: &Moment) -> Option<String> {
    // (A coupling that states no efficiency is taken as lossless, which the answer says once: `Effort::assumed`.)
    run.model.refusal_lossless(k, m)
}

impl Tool for Effort<'_> {
    fn measures(&self) -> Option<ink_claims::Worse> {
        Some(ink_claims::Worse::More)
    }

    fn ask(&self, about: &About, _: &mut Claims) -> Option<Fact> {
        let (path, effort, unit, rate_unit) = self.asked(&about.on, false)?;
        let speed = about.what == format!("speed-at-peak-{effort}");
        if !speed && about.what != format!("peak-{effort}") {
            return None;
        }
        // (A peak's speed is of the peak: with no effort there is none.)
        let (run, k, mut fact) = match self.begun(about, path, if speed { rate_unit } else { unit }, true) {
            Ok(begun) => begun,
            Err(open) => return Some(open),
        };
        if let Some((t, why)) = run.at.iter().find_map(|(_, m)| Some((m.t, refused(run, k, m)?))) {
            (fact.held, fact.strength, fact.notes) = (Some(Verdict::NotSettled(Unsettled::Undecided)), None, vec![format!("refused at t {t}: {why}")]);
            return Some(fact);
        }
        // The largest effort among the moments read, either way round.
        let m = &run.at.iter().max_by(|x, y| x.1.effort[k].abs().total_cmp(&y.1.effort[k].abs())).expect("a run reads a moment").1;
        let (peak, rate) = (m.effort[k], m.rate[k]);
        // The span's two ends and each step's middle, those read.
        fact.covers = vec![("t".into(), if run.span.is_some() { Covers::These(run.at.len()) } else { Covers::One })];
        fact.notes.push(format!("the peak is {peak} at t {}, moving at {rate}, in the solver's units", m.t));
        if speed {
            fact.value = Some(rate.abs());
        } else {
            fact.value = Some(peak.abs());
            // How far the solver says the real effort may be from the one it read, either way (bodies read
            // round, a coupling's loss, and what its bodies' proven mass, centre and inertia can move it:
            // Moment::mass); and upward, over a motion, how low its steps may read a peak between their moments.
            let most = |of: &dyn Fn(&Moment) -> f64| run.at.iter().map(|(_, m)| of(m)).fold(0.0, f64::max);
            let off = most(&|m| m.bound[k]) + most(&|m| m.loss_bound[k]) + most(&|m| m.mass[k]);
            let between = if run.span.is_some() { OFF_LINE * peak.abs() } else { 0.0 };
            fact.band = (off > 0.0 || between > 0.0).then(|| ((peak.abs() - off).max(0.0), peak.abs() + off + between));
            fact.notes.extend((between > 0.0).then(between_note));
        }
        Some(fact)
    }
}

fn between_note() -> String {
    format!("between the moments read nothing is bounded: they are chosen so that a peak is read at most {}% low, which its band's upper end allows", 100.0 * OFF_LINE)
}

impl Tool for AtMoment<'_> {
    fn measures(&self) -> Option<ink_claims::Worse> {
        Some(ink_claims::Worse::More)
    }

    fn ask(&self, about: &About, claims: &mut Claims) -> Option<Fact> {
        if let Some(load) = self.effort.load(about, self.moment, self.early, claims) {
            return Some(load);
        }
        if let Some(tension) = self.effort.tension(about, self.moment, self.early) {
            return Some(tension);
        }
        let (path, effort, unit, rate_unit) = self.effort.asked(&about.on, true)?;
        let speed = about.what == "speed";
        if !speed && about.what != effort {
            return None;
        }
        // A speed is the written motion's: read from it alone where no run has started.
        if speed
            && let Some(k) = self.effort.actuators.iter().position(|(p, ..)| p == path)
            && let Some(fact) = self.effort.kin_speed(about, k, self.moment, rate_unit)
        {
            return Some(fact);
        }
        // On a frame of the part not along the joint's axis (a frame across it, or turned against it), the
        // effort it passes is not the joint's: the moment (a force, on a slide) about that frame's z, of what
        // the joint's reaction is about the frame's origin.
        let across = (!speed).then(|| about.on.strip_prefix(path).and_then(|rest| rest.strip_prefix(' '))).flatten().and_then(|name| {
            let axis = self.effort.joints.iter().find(|(p, ..)| p == path)?.2;
            let local = self.effort.placed.get(&(path.to_string(), name.to_string()));
            let along = local.is_some_and(|l| l[3 * axis + 2] >= 1.0 - 1e-12 && (0..3).all(|k| k == axis || l[9 + k].abs() <= 1e-12 * (1.0 + l[9 + axis].abs())));
            (!along).then(|| (name.to_string(), local.copied()))
        });
        let (run, k, mut fact) = match self.effort.begun(about, path, if speed { rate_unit } else { unit }, !speed) {
            Ok(begun) => begun,
            Err(open) => return Some(open),
        };
        if let Some((name, local)) = across {
            let turn = effort == "torque";
            let Some(local) = local else {
                (fact.held, fact.strength, fact.notes) = (
                    Some(Verdict::NotAvailable(Unavailable::Outside)),
                    None,
                    vec![format!("its frame {name} is not placed on its part by numbers alone, so the moment about it is not read")],
                );
                return Some(fact);
            };
            if let Some(why) = self.effort.other_moment(run, self.moment).filter(|_| self.early) {
                (fact.held, fact.strength, fact.notes) = (Some(Verdict::NotSettled(Unsettled::Undecided)), None, vec![why]);
                return Some(fact);
            }
            let m = &run.at.get(self.moment)?.1;
            match ink_dynamics::reaction_in(run, path, m.t, &local)? {
                Err(why) => {
                    (fact.held, fact.strength, fact.notes) = (Some(Verdict::NotSettled(Unsettled::Undecided)), None, vec![format!("refused at t {}: {why}", m.t)]);
                }
                Ok((force, moment)) => {
                    let x = if turn { moment[2] } else { force[2] };
                    let masses: f64 = ink_dynamics::reaction_moves(run, path, m.t, Some(&local)).iter().map(|(df, dm)| if turn { dm[2] } else { df[2] }.abs()).sum();
                    fact.value = Some(x);
                    fact.band = Some((x - masses, x + masses));
                    fact.covers = vec![("t".into(), Covers::One)];
                    fact.notes.insert(0, format!("at t = {}", m.t));
                    fact.notes.push(format!(
                        "the {} about its frame {name}'s z of the model's reaction at {path}, there: (reaction {path} (frame-of {path} {name})), within {} {} either way, what the masses it rests on can move it, each within its proven bounds",
                        if turn { "moment" } else { "force along" },
                        (masses * 1e3).ceil() / 1e3,
                        if turn { "N*mm" } else { "N" }
                    ));
                    if run.span.is_some() {
                        fact.left_out.push("a reaction's peak between the moments read: they are chosen for the efforts".into());
                    }
                }
            }
            return Some(fact);
        }
        if let Some(why) = self.effort.other_moment(run, self.moment).filter(|_| self.early) {
            (fact.held, fact.strength, fact.notes) = (Some(Verdict::NotSettled(Unsettled::Undecided)), None, vec![why]);
            return Some(fact);
        }
        let m = &run.at.get(self.moment)?.1;
        fact.notes.insert(0, format!("at t = {}", m.t));
        // A speed is the motion's, whatever the effort there: refused only where the motion itself is not found.
        if let Some(why) = if speed { m.unsolved.clone() } else { refused(run, k, m) } {
            (fact.held, fact.strength, fact.notes) = (Some(Verdict::NotSettled(Unsettled::Undecided)), None, vec![format!("refused at t {}: {why}", m.t)]);
            return Some(fact);
        }
        fact.covers = vec![("t".into(), Covers::One)];
        if speed {
            fact.value = Some(m.rate[k]);
            // Over a motion, every speed from this moment's to its neighbours': what the joint passes between them.
            let near = [self.moment.checked_sub(1), Some(self.moment + 1)].into_iter().flatten().filter_map(|i| run.at.get(i)).filter(|(_, m)| m.unsolved.is_none());
            let (low, high) = near.map(|(_, m)| m.rate[k]).fold((m.rate[k], m.rate[k]), |(low, high), rate| (low.min(rate), high.max(rate)));
            if low < high {
                fact.band = Some((low, high));
                fact.notes.push("its band reaches to its speed at the moments before and after".into());
            }
            return Some(fact);
        }
        let x = m.effort[k];
        fact.value = Some(x);
        // Either way, what the solver says of this moment's read; away from nothing, over a motion, what
        // an effort between this moment and the next may add: the share of the motion's peak a peak is read low by.
        let off = m.bound[k] + m.loss_bound[k] + m.mass[k];
        let peak = run.at.iter().map(|(_, m)| m.effort[k].abs()).fold(0.0, f64::max);
        let between = if run.span.is_some() { OFF_LINE * peak } else { 0.0 };
        if off > 0.0 || between > 0.0 {
            fact.band = Some((x - off - if x < 0.0 { between } else { 0.0 }, x + off + if x >= 0.0 { between } else { 0.0 }));
            fact.notes.extend((between > 0.0).then(between_note));
        }
        Some(fact)
    }
}
