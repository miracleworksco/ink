use super::bodies::{path_numbers, strip};
use super::couplings::belt_pieces;
use super::kinematics::{nums, partials, placement, xf_entries};
use super::sampling::Tracker;
use super::spin::{about_axis, known_round, own_round, round_part, spin, spun, turning};
use super::{
    Actuator, Applied, Balancing, Body, Coordinate, Coupled, End, Group, Input, M3, MadeRound, Model, NO_MASS, Num, Reaction, Rides, Sets, Speed, SpinCentre, States, V3, mv,
};
use ink_carried::Carried;
use ink_core::graph::{Graph, Node, NodeId, Var};
use ink_core::v3::{dot, norm};
use ink_lang::view::ModelView;
use ink_lang::{Doc, InkError, Value};
use std::collections::{HashMap, HashSet};
use std::sync::Mutex;

/// The placement in the world of what instance path `path` places, from
/// an item under it: its links down to it composed.
fn path_placement(v: &mut ModelView, items: &[ink_lang::Item], path: &str) -> Option<ink_lang::Xf> {
    let (it, i) = items.iter().find_map(|it| Some((it, it.links.iter().position(|l| l.path == path)?)))?;
    let mut x = v.doc.identity();
    for l in &it.links[..=i] {
        x = v.doc.compose(&l.at, &x);
    }
    Some(x)
}

/// A reaction's frame in the world (Reaction): the target's, the world's;
/// a path's joint's, its WHERE moved by its motion (the links above it,
/// then its placement without the child's own below it: Joint::below);
/// else its placement.
pub(super) fn path_frame(v: &mut ModelView, items: &[ink_lang::Item], target: &str, path: &str) -> Option<[NodeId; 12]> {
    if path == target {
        let id = v.doc.identity();
        return Some(xf_entries(&id));
    }
    let (it, i) = items.iter().find_map(|it| Some((it, it.links.iter().position(|l| l.path == path)?)))?;
    let mut above = v.doc.identity();
    for l in &it.links[..i] {
        above = v.doc.compose(&l.at, &above);
    }
    let l = &it.links[i];
    let own = match &l.joint {
        Some(j) => {
            let back = v.doc.inverse(&j.below);
            v.doc.compose(&back, &l.at)
        }
        None => l.at,
    };
    let x = v.doc.compose(&own, &above);
    Some(xf_entries(&x))
}

/// The joint at instance path `path`: whether it turns, and its axis in the world at the view's moment, a
/// point on it and its direction.
pub(crate) fn joint_axis(v: &mut ModelView, items: &[ink_lang::Item], target: &str, path: &str) -> Option<(bool, V3, V3)> {
    let joint = items.iter().flat_map(|it| &it.links).find(|l| l.path == path)?.joint.as_ref()?;
    let (turn, axis) = (joint.turn, joint.axis);
    let frame = path_frame(v, items, target, path)?.map(|n| v.doc.scalar(n, &v.ctx));
    Some((turn, [frame[9], frame[10], frame[11]], [frame[axis], frame[3 + axis], frame[6 + axis]]))
}

/// The joint the instance at `path` under `target` rides (its own, or the one it is placed on), by path.
pub fn riding(doc: &mut ink_lang::Doc, target: &str, path: &str) -> Option<String> {
    let items = doc.items(target).ok()?;
    rides(&items, path, target).0.map(|(joint, _)| joint)
}

/// The joint instance path `path` rides, its path and drive, the nearest
/// found going out from it (its own, then what it is placed on, then what
/// it is under), and what supports it: the paths so passed, up to the
/// target.
pub(super) fn rides(items: &[ink_lang::Item], path: &str, target: &str) -> (Option<(String, NodeId)>, HashSet<String>) {
    let links: HashMap<&str, &ink_lang::Link> = items.iter().flat_map(|it| &it.links).map(|l| (l.path.as_str(), l)).collect();
    let (mut drive, mut supports, mut todo) = (None, HashSet::new(), vec![path.to_string()]);
    while let Some(q) = todo.pop() {
        if q == target || !supports.insert(q.clone()) {
            continue;
        }
        let Some(l) = links.get(q.as_str()) else { continue };
        if drive.is_none() {
            drive = l.joint.as_ref().map(|j| (q.clone(), j.drive));
        }
        let parent = q.rsplit_once('/').map_or(String::new(), |x| x.0.to_string());
        todo.push(parent.clone());
        todo.extend(l.on.iter().map(|o| format!("{parent}/{o}")));
    }
    supports.remove(path);
    (drive, supports)
}

/// The instance paths each of `paths` carries: those under it, and those
/// placed on what it carries (each instance their WHERE reads,
/// Instance::on, carried), with the instances under those; and the bodies
/// placed partly on what it carries (some of what their WHERE reads
/// carried, some not).
pub(crate) fn carried(items: &[ink_lang::Item], paths: &[String]) -> (Vec<HashSet<String>>, Vec<Vec<String>>) {
    let mut links: Vec<(&str, &[String])> = Vec::new();
    let mut seen = HashSet::new();
    for l in items.iter().flat_map(|it| &it.links) {
        if seen.insert(l.path.as_str()) {
            links.push((&l.path, &l.on));
        }
    }
    paths
        .iter()
        .map(|p| {
            let below = format!("{p}/");
            let mut set: HashSet<String> =
                items.iter().map(|it| it.label.as_str()).chain(links.iter().map(|l| l.0)).filter(|q| *q == p || q.starts_with(&below)).map(String::from).collect();
            loop {
                let before = set.len();
                for &(q, on) in &links {
                    let parent = q.rsplit_once('/').map_or("", |x| x.0);
                    if !set.contains(q) && (set.contains(parent) || (!on.is_empty() && on.iter().all(|o| set.contains(&format!("{parent}/{o}"))))) {
                        set.insert(q.to_string());
                    }
                }
                if set.len() == before {
                    break;
                }
            }
            // The bodies under what it carries.
            for it in items {
                if !set.contains(&it.label) && it.links.iter().any(|l| set.contains(&l.path)) {
                    set.insert(it.label.clone());
                }
            }
            let partly = items
                .iter()
                .filter(|it| {
                    !set.contains(&it.label)
                        && it.links.iter().any(|l| {
                            let parent = l.path.rsplit_once('/').map_or("", |x| x.0);
                            l.on.iter().any(|o| set.contains(&format!("{parent}/{o}")))
                        })
                })
                .map(|it| it.label.clone())
                .collect();
            (set, partly)
        })
        .unzip()
}

/// Whether `n` reads t other than through the nodes `through` (the
/// coordinates as passed): a value a solve sets reads what its equations
/// do.
fn reads_t_besides(d: &Doc, n: NodeId, through: &HashSet<NodeId>) -> bool {
    let mut seen = HashSet::new();
    let mut todo = vec![n];
    while let Some(x) = todo.pop() {
        if through.contains(&x) || !seen.insert(x) {
            continue;
        }
        match d.g.nodes[x as usize] {
            Node::Var(Var::T) => return true,
            Node::Var(Var::Param(k)) => {
                if let Some(s) = d.solves.iter().find(|s| s.owns.contains(&(k as usize))) {
                    todo.extend(d.form_motion_inputs(&d.ext[s.ext]));
                }
            }
            node => {
                let (ch, k) = node.children();
                todo.extend(&ch[..k]);
            }
        }
    }
    false
}

impl Model {
    /// The mechanism under `target` in a view of a load with coordinates
    /// (Doc::with_coordinates): its moving bodies' masses measured (cached
    /// per body), their placements' derivatives compiled (through solves,
    /// the implicit function's: Doc::derivative), its loads (the file
    /// run's on paths under the target), and the reactions of the paths
    /// `wanted` (`claimed`: asked for by claims), every body they carry
    /// measured. Err: a body that moves along a coordinate without mass, an
    /// actuator whose drive reads no coordinate, a path placed nowhere
    /// under the target.
    pub fn compile(v: &mut ModelView, target: &str, gravity: [f64; 3], wanted: &[String], claimed: bool) -> Result<Model, InkError> {
        Model::compile_in(v, target, gravity, wanted, claimed, false, false, |_, _, _| ()).map(|(model, ())| model)
    }

    /// `compile`, for one who reads the motion whatever the masses: a moving body with no mass does not
    /// refuse the model, it is left out and named (`massless`), and what it would change is refused.
    pub fn compile_kinematic(v: &mut ModelView, target: &str, gravity: [f64; 3], wanted: &[String]) -> Result<Model, InkError> {
        Self::compile_partial(v, target, gravity, wanted, false)
    }

    /// `compile_kinematic`, preserving `claimed` for claimed reactions.
    pub fn compile_partial(v: &mut ModelView, target: &str, gravity: [f64; 3], wanted: &[String], claimed: bool) -> Result<Model, InkError> {
        Model::compile_in(v, target, gravity, wanted, claimed, false, true, |_, _, _| ()).map(|(model, ())| model)
    }

    /// Model::compile; `sprung`: the couplings stating a stiffness are left
    /// out, springs (Model::compile_springs) and no relations.
    /// Compile the caller's joint derivatives before binding the model's
    /// parameter slots and solve order; these derivatives can add rates of their own.
    pub(super) fn compile_in<T>(
        v: &mut ModelView,
        target: &str,
        gravity: [f64; 3],
        wanted: &[String],
        claimed: bool,
        sprung: bool,
        kinematic: bool,
        joint_derivatives: impl FnOnce(&mut ModelView, &[u32], Vec<(String, ink_lang::Joint)>) -> T,
    ) -> Result<(Model, T), InkError> {
        let items = v.doc.items(target).map_err(InkError::msg)?;
        let coords = v.doc.coords.clone();
        let mut joints: Vec<(String, ink_lang::Joint)> = Vec::new();
        for l in items.iter().flat_map(|it| &it.links) {
            if let Some(j) = &l.joint
                && !joints.iter().any(|(p, _)| *p == l.path)
            {
                joints.push((l.path.clone(), j.clone()));
            }
        }
        // Params a solve on t sets: what reads one moves.
        let timed: Vec<u32> = (0..v.doc.solves.len()).filter(|&i| v.doc.solves[i].on_t).flat_map(|i| v.doc.solves[i].owns.iter().map(|&k| k as u32)).collect();
        let g = &v.doc.g;
        let moves = |c: &ink_lang::Coord| g.depends_on(c.value, |x| x == Var::T || matches!(x, Var::Param(k) if timed.contains(&k)));
        // A coordinate is an argument some joint's drive reads. One passed
        // as a constant that a body's shape reads is a size (a length, a
        // lead): it does not move, so it is no coordinate.
        let used: Vec<(&ink_lang::Coord, bool)> = coords
            .iter()
            .filter(|c| {
                let e = Var::Param(c.param as u32);
                joints.iter().any(|(_, j)| g.depends_on(j.drive, |x| x == e)) && (moves(c) || !items.iter().any(|it| g.depends_on(it.shape, |x| x == e)))
            })
            .map(|c| (c, moves(c)))
            .collect();
        let eps: Vec<u32> = used.iter().map(|(c, _)| c.param as u32).collect();
        // The coordinates n moves with: those it reads, directly or
        // through the solves it rests on.
        let along_of = |d: &Doc, n: NodeId| -> Vec<usize> {
            let solves = d.solves_read(&[n]);
            (0..eps.len()).filter(|&k| d.g.depends_on(n, |x| x == Var::Param(eps[k])) || solves.iter().any(|&i| d.solves[i].reads.contains(&(eps[k] as usize)))).collect()
        };
        // The solves the mechanism moves through: on t, or reading a
        // coordinate (a loop closed at a pose held still moves with it),
        // directly or through what they read; and the params they set.
        let through: Vec<usize> = (0..v.doc.solves.len())
            .filter(|&i| v.doc.solves[i].on_t || v.doc.solves_under(&v.doc.solves[i].owns).iter().any(|&j| v.doc.solves[j].reads.iter().any(|k| eps.contains(&(*k as u32)))))
            .collect();
        let solved: Vec<u32> = through.iter().flat_map(|&i| v.doc.solves[i].owns.iter().map(|&k| k as u32)).collect();
        // What changes over the motion reads t or a coordinate: directly
        // (a constant one a shape reads is a size), or through such a
        // solve.
        let changes = |g: &Graph, n: NodeId| g.depends_on(n, |x| x == Var::T || matches!(x, Var::Param(k) if solved.contains(&k) || eps.contains(&k)));
        // The commands an actuator supplies in a simulation, and the
        // coordinates simulated: each `(simulate ...)` is the argument a
        // coordinate is passed as, itself.
        let commands: Vec<(String, NodeId)> = v.doc.forms::<crate::Command>().filter(|c| c.run).map(|c| (c.path.clone(), c.effort)).collect();
        // A simulation is of root's mechanism: another target's model reads
        // its own coordinates as written.
        let mut simulated = Vec::new();
        for s in v.doc.forms::<crate::Simulated>().filter(|_| target == "root") {
            let Some(k) = used.iter().position(|(c, _)| c.value == s.field) else {
                return Err(InkError {
                    hint: Some("pass it as an argument of an instance the target places, itself, read by a joint's drive: (instance arm (arm (simulate 30)))".into()),
                    ..InkError::msg(format!("(simulate ...) at line {} is no coordinate's motion: it is the motion of the argument it is passed as", s.line))
                });
            };
            simulated.push((k, s.owns));
        }
        let d = &mut v.doc;
        let passed: HashSet<NodeId> = used.iter().map(|(c, _)| c.passed).collect();
        let mut coupling = Sets((0..used.len()).collect());
        let mut actuators = Vec::new();
        // Actuators held at a pose a body's shape follows (an argument both
        // the drive and a shape read, so a size), each with those arguments'
        // params: their efforts leave those bodies out.
        let mut shaped: Vec<(usize, Vec<u32>)> = Vec::new();
        for (path, j, lim) in joints.iter().filter_map(|(p, j)| Some((p, j, crate::joint::effort(j)?))) {
            let along: Vec<NodeId> = eps.iter().map(|&k| d.derivative(j.drive, Var::Param(k))).collect();
            // Those its drive reads (a rate not folded to 0 of one it does
            // not read moves nothing: another machine's coordinate).
            let reads = along_of(d, j.drive);
            let moved: Vec<usize> = (0..along.len()).filter(|&k| d.g.constant_value(along[k]) != Some(0.0) && reads.contains(&k)).collect();
            let sizes: Vec<u32> = coords.iter().map(|c| c.param as u32).filter(|k| !eps.contains(k) && d.g.depends_on(j.drive, |x| x == Var::Param(*k))).collect();
            if moved.is_empty() && !sizes.is_empty() {
                shaped.push((actuators.len(), sizes));
            } else if moved.is_empty() {
                return Err(InkError {
                    hint: Some(format!(
                        "pass the motion as an argument of an assembly the target places, read by the joint's drive: (assembly (arm (angle 0)) ... (turn z angle) (effort-limits ...)) and (instance a (arm {})) in {target}",
                        j.text
                    )),
                    ..InkError::msg(format!("actuator {path}: its drive, {}, reads no coordinate (an argument {target}'s instances pass), so no effort is along it", j.text))
                });
            }
            coupling.join(&moved);
            let rate = d.derivative(j.drive, Var::T);
            let num = |n| Num::of(&d.g, n);
            let source = format!("{}:{}", lim.file, lim.line);
            let command = commands.iter().find(|c| c.0 == *path).map(|c| num(c.1));
            actuators.push(Actuator {
                path: path.clone(),
                turn: j.turn,
                along: along.iter().map(|&a| num(a)).collect(),
                rate: num(rate),
                lo: num(lim.lo),
                hi: num(lim.hi),
                coupled: moved,
                leaves_out: Vec::new(),
                source,
                command,
            });
        }
        // The relations' couplings under the target, and the passive joints
        // (no effort-limits) that carry one end of one, whether this motion
        // moves them or holds them (a CoreXY's gantry on a move along x):
        // each passes nothing along its motion but through its couplings,
        // so its balance gives their tensions. Their reactions are found
        // with those asked for.
        let relations: Vec<crate::coupling::Related> =
            crate::coupling::related(d, target).into_iter().filter(|r| r.couple.as_ref().is_some_and(|c| !(sprung && c.lines.iter().any(|l| l.stiffness.is_some())))).collect();
        let passive: Vec<String> = joints.iter().filter(|(_, j)| crate::joint::effort(j).is_none()).map(|(p, _)| p.clone()).collect();
        let (passive_sets, _) = carried(&items, &passive);
        let reaches = |set: &HashSet<String>, r: &crate::coupling::Related| set.contains(&r.a) != set.contains(&r.b);
        let mut all: Vec<String> = wanted.to_vec();
        let mut balancing: Vec<(String, bool, usize)> = Vec::new();
        for (p, set) in passive.iter().zip(&passive_sets) {
            if relations.iter().any(|r| reaches(set, r)) {
                let j = &joints.iter().find(|(q, _)| q == p).unwrap().1;
                balancing.push((p.clone(), j.turn, j.axis));
                if !all.contains(p) {
                    all.push(p.clone());
                }
            }
        }
        // The passive joints moved by an actuator's coordinate, read in their drive directly (not through a
        // loop's closure, which a coupler between them closes), that no relation reaches: no part carrying
        // the force between them (a belt, a cord, a mesh) is in the model. Each with the actuators moving it.
        let linked: Vec<crate::coupling::Related> = crate::coupling::related(d, target).into_iter().filter(|r| r.couple.is_some()).collect();
        let mut unlinked: Vec<(String, Vec<String>)> = Vec::new();
        for (p, set) in passive.iter().zip(&passive_sets) {
            let j = &joints.iter().find(|(q, _)| q == p).unwrap().1;
            let direct: Vec<usize> = (0..eps.len())
                .filter(|&k| {
                    let along = d.derivative(j.drive, Var::Param(eps[k]));
                    d.g.depends_on(j.drive, |x| x == Var::Param(eps[k])) && d.g.constant_value(along) != Some(0.0)
                })
                .collect();
            let by: Vec<String> = actuators.iter().filter(|a| a.coupled.iter().any(|k| direct.contains(k))).map(|a| a.path.clone()).collect();
            if !by.is_empty() && !linked.iter().any(|r| reaches(set, r)) {
                unlinked.push((p.clone(), by));
            }
        }
        let wanted: &[String] = &all;
        // What the effort reads: the coordinates as passed, the actuators'
        // drives and limits, and (below) each moving body's shape, density,
        // placement and published spin.
        let mut read: Vec<NodeId> =
            used.iter().map(|(c, _)| c.value).chain(joints.iter().filter_map(|(_, j)| crate::joint::effort(j).map(|e| [j.drive, e.lo, e.hi])).flatten()).collect();
        read.extend(commands.iter().map(|c| c.1));
        // What only the reactions read (the bodies they carry that no
        // effort does: still, or moved by t alone).
        let mut read_carried: Vec<NodeId> = Vec::new();
        let mut model = Model {
            gravity,
            coords: used
                .iter()
                .map(|&(c, moves)| {
                    let dq = d.derivative(c.value, Var::T);
                    let ddq = d.derivative(dq, Var::T);
                    let dddq = d.derivative(ddq, Var::T);
                    let d4q = d.derivative(dddq, Var::T);
                    Coordinate {
                        name: c.name.clone(),
                        q: Num::timed(d, c.value),
                        dq: Num::timed(d, dq),
                        ddq: Num::timed(d, ddq),
                        dddq: Num::timed(d, dddq),
                        d4q: Num::timed(d, d4q),
                        moves,
                    }
                })
                .collect(),
            actuators,
            bodies: Vec::new(),
            strips: Vec::new(),
            loads: Vec::new(),
            reactions: Vec::new(),
            couplings: Vec::new(),
            balances: Vec::new(),
            unseparated: Vec::new(),
            trains: Vec::new(),
            behind: Vec::new(),
            tracked: Vec::new(),
            groups: Vec::new(),
            reads: Vec::new(),
            reads_carried: Vec::new(),
            span: None,
            params: Vec::new(),
            driven_by_nothing: Vec::new(),
            held_by_nothing: Vec::new(),
            unlinked,
            made_round: Vec::new(),
            spin_centres: Vec::new(),
            deforming: Vec::new(),
            carried: Vec::new(),
            measured: 0,
            published: 0,
            spin_inertia: Vec::new(),
            res: (f64::INFINITY, 0.0),
            unsettled: Vec::new(),
            unsettled_carried: Vec::new(),
            massless: Vec::new(),
            sets: Vec::new(),
            solves: Vec::new(),
            follow: None,
            states: d.state_t.map(|t0| States { of: used.iter().map(|(c, _)| (c.param, c.jet)).collect(), t0 }),
            simulated,
            roots: Mutex::new(HashMap::new()),
            leaves: Vec::new(),
            speeds: Vec::new(),
            springs: None,
            moments: Mutex::new(HashMap::new()),
            chosen: Mutex::new(HashMap::new()),
            twin: None,
        };
        // At a state, every joint's speed limits bound the path's speed.
        if model.states.is_some() {
            for (path, j) in &joints {
                if let Some((lo, hi)) = j.speed_limits {
                    let rate = d.derivative(j.drive, Var::T);
                    model.speeds.push(Speed { path: path.clone(), rate: Num::of(&d.g, rate), lo: Num::of(&d.g, lo), hi: Num::of(&d.g, hi) });
                }
            }
        }
        // The reactions asked for: each path's frame, and which bodies it
        // carries.
        for p in wanted {
            let frame = path_frame(v, &items, target, p).ok_or_else(|| InkError::msg(format!("a reaction of {p}: no instance of that path is placed under {target}")))?;
            let frame = nums(&v.doc.g, &frame);
            model.reactions.push(Reaction {
                path: p.clone(),
                frame,
                carries: Vec::new(),
                without_mass: Vec::new(),
                leaves_out: Vec::new(),
                coords: Vec::new(),
                moves_with: Vec::new(),
                moved_by_t: Vec::new(),
                belts: Vec::new(),
                unbounded: Vec::new(),
            });
        }
        let (carried_paths, partly) = carried(&items, wanted);
        for (r, set) in model.reactions.iter_mut().zip(&carried_paths) {
            r.carries = items.iter().filter(|it| set.contains(&it.label)).map(|it| it.label.clone()).collect();
        }
        // Each body moving along a coordinate, with them (Reaction::moves_with).
        let mut along_coords: Vec<(String, Vec<usize>)> = Vec::new();
        let carriers: Vec<Vec<usize>> = items.iter().map(|it| (0..wanted.len()).filter(|&r| carried_paths[r].contains(&it.label)).collect()).collect();
        // What moves, each with its density, and what a reaction carries.
        // A body moving along a coordinate without one refuses the answer
        // before anything is measured; one only t moves adds nothing to an
        // effort. A body whose `carried` says its material moves (its path
        // or shift moves, its profile still) is counted so, whatever its
        // shape does; another whose shape moves is left out (deforming). A
        // reaction counts every body it carries, still or moved by t alone
        // too; one without mass leaves that reaction unfound.
        let mut moving = Vec::new();
        let mut without: Vec<(String, Vec<String>)> = Vec::new();
        for (i, it) in items.iter().enumerate() {
            let of = carriers[i].clone();
            let Some((name, xf)) = it.instance.clone() else { continue };
            let entries: Vec<NodeId> = xf.r.iter().flatten().chain(&xf.t).copied().collect();
            let attrs = v.doc.body(&name).map(|b| b.attrs.clone()).unwrap_or_default();
            let carried = Carried::in_attrs(&v.doc.g, &attrs)
                .map_err(|e| InkError::msg(format!("{}: {e}", it.label)))?
                .filter(|c| !changes(&v.doc.g, c.profile) && path_numbers(c).iter().any(|&n| changes(&v.doc.g, n)));
            let deforms = carried.is_none() && changes(&v.doc.g, it.shape);
            // What moves it: its placement, and a carried body's path.
            let mut moved_by = entries.clone();
            moved_by.extend(carried.iter().flat_map(path_numbers));
            let moves = deforms || moved_by.iter().any(|&e| changes(&v.doc.g, e));
            if !moves && of.is_empty() {
                continue;
            }
            let mut along: Vec<usize> = moved_by.iter().flat_map(|&e| along_of(&v.doc, e)).collect();
            along.sort_unstable();
            along.dedup();
            if moved_by.iter().any(|&e| reads_t_besides(&v.doc, e, &passed)) {
                model.driven_by_nothing.push(it.label.clone());
                // Held by nothing declared: placed by a formula of t, riding
                // no joint (one riding a joint, whatever moves it or the
                // joint's frame, is held by that joint).
                if rides(&items, &it.label, target).0.is_none() {
                    model.held_by_nothing.push(it.label.clone());
                    for &r in &of {
                        model.reactions[r].moved_by_t.push(it.label.clone());
                    }
                }
            }
            if along.is_empty() && !deforms && of.is_empty() {
                continue;
            }
            let who = v.doc.named(it);
            let rho = ink_lang::materials::mass_density(&mut v.doc, &who, &attrs, None)?;
            if rho.is_none() && !deforms {
                if along.is_empty() {
                    for &r in &of {
                        model.reactions[r].without_mass.push(it.label.clone());
                    }
                    continue;
                }
                if kinematic {
                    // Left out, and named: its coordinates rest on each other as with any body, so every
                    // effort it would change is found by them (`blocked`).
                    for &r in &of {
                        model.reactions[r].without_mass.push(it.label.clone());
                        model.reactions[r].coords.extend(&along);
                    }
                    coupling.join(&along);
                    along_coords.push((it.label.clone(), along.clone()));
                    model.massless.push((it.label.clone(), name.clone(), along));
                    continue;
                }
                match without.iter_mut().find(|(b, _)| *b == name) {
                    Some((_, ls)) => ls.push(it.label.clone()),
                    None => without.push((name.clone(), vec![it.label.clone()])),
                }
            }
            for &r in &of {
                model.reactions[r].coords.extend(&along);
            }
            if !along.is_empty() {
                along_coords.push((it.label.clone(), along.clone()));
            }
            moving.push((it, xf, who, attrs, entries, deforms, carried, along, rho, of));
        }
        if !without.is_empty() {
            let list: Vec<String> = without.iter().map(|(b, ls)| format!("{b} ({})", ls.join(", "))).collect();
            return Err(InkError {
                hint: Some(NO_MASS.into()),
                open: Some(ink_lang::Open::Lacks),
                ..InkError::msg(format!("these moving bodies have no mass: {}", list.join("; ")))
            });
        }
        // Read for its motion (`kinematic`): a body whose coordinates rest on a massless body's is not
        // measured. Every effort and reaction it is in is refused for that body (`blocked`,
        // `blocks_reaction`: its coordinates are joined to theirs here, as measuring it would join them), so
        // its inertia is read by nothing, and weighing a robot's every part to refuse them all is the cost.
        if kinematic && !model.massless.is_empty() {
            // The sets as they will be at the least: the massless bodies' and the actuators' joins made, and
            // each body's own.
            let mut sets = Sets(coupling.0.clone());
            moving.iter().for_each(|m: &(_, _, _, _, _, _, _, Vec<usize>, _, _)| sets.join(&m.7));
            let without: Vec<usize> = model.massless.iter().flat_map(|(_, _, along)| along.clone()).map(|k| sets.find(k)).collect();
            let (kept, dropped): (Vec<_>, Vec<_>) = moving.into_iter().partition(|m| !m.7.iter().any(|&k| without.contains(&sets.find(k))));
            for m in &dropped {
                coupling.join(&m.7);
            }
            moving = kept;
        }
        let mut tracker = Tracker::default();
        let (mut body_labels, mut strip_labels, mut load_paths): (Vec<String>, Vec<String>, Vec<String>) = (Vec::new(), Vec::new(), Vec::new());
        let mut strip_ends: Vec<Vec<[[NodeId; 3]; 2]>> = Vec::new();
        let turning = turning(v, &items, &changes);
        // Deforming bodies, each with the coordinates it moves with and its
        // mass.
        let mut deforming: Vec<(String, Vec<usize>, Option<f64>)> = Vec::new();
        for (it, xf, who, attrs, entries, deforms, carried, along, rho, of) in moving {
            let rho_node = rho.as_ref().map(|r| r.0);
            if let Some((rho, _)) = &rho {
                // The state model measures the body's inertia once, then uses
                // it at every simulated pose. A density that changes with
                // the motion cannot be validated or represented by that mass.
                if v.doc.states() && changes(&v.doc.g, *rho) {
                    return Err(InkError::msg(format!("{who}: density depends on time or a moving coordinate; this simulation uses one mass and inertia throughout the motion")));
                }
                v.want(&[*rho])?;
            }
            if let Some(c) = carried {
                let (rho, what) = rho.expect("a density");
                coupling.join(&along);
                read.extend(path_numbers(&c));
                read.extend([c.profile, c.length, rho]);
                read.extend(&entries);
                let Some((strip, ends, m)) = strip(v, &it.label, &c, &xf, rho, &along, &eps, of, &mut tracker) else {
                    return Err(InkError::msg(format!("{who}: its carried profile is not proven bounded, so its mass is not measured")));
                };
                model.measured += 1;
                model.published += usize::from(what.starts_with("published mass"));
                model.res = (model.res.0.min(m.res), model.res.1.max(m.res));
                if let Some(why) = m.unsettled() {
                    model.unsettled.push(format!("{}: {why}", it.label));
                }
                model.carried.push((it.label.clone(), m.sums[0]));
                model.strips.push(strip);
                strip_ends.push(ends);
                strip_labels.push(it.label.clone());
                continue;
            }
            let found = rho.map(|(rho, what)| (v.inertia(it.shape, Some(rho), None), what));
            if deforms {
                let mut reads = along_of(&v.doc, it.shape);
                reads.extend(along);
                coupling.join(&reads);
                let mass = found.and_then(|(m, _)| m).map(|m| m.sums[0]);
                for &r in of.iter().filter(|_| mass != Some(0.0)) {
                    model.reactions[r].leaves_out.push(match mass {
                        Some(m) => format!("{} ({m:.4} g)", it.label),
                        None => format!("{} (no mass)", it.label),
                    });
                }
                model.deforming.push((it.label.clone(), mass));
                deforming.push((it.label.clone(), reads, mass));
                continue;
            }
            coupling.join(&along);
            // One moving along no coordinate is a reaction's alone: no
            // effort reads it.
            let (read, unsettled) = if along.is_empty() { (&mut read_carried, &mut model.unsettled_carried) } else { (&mut read, &mut model.unsettled) };
            read.push(it.shape);
            read.extend(rho_node);
            read.extend(&entries);
            let j = attrs.iter().find_map(|(k, v)| match v {
                Value::Field(f, _) if k == "spin-inertia" => Some(*f),
                _ => None,
            });
            read.extend(j);
            let (Some(m), what) = found.expect("a density") else { return Err(InkError::msg(format!("{who}: not proven bounded, so its mass is not measured"))) };
            // One reaching past every box (a ground half-space) is measured
            // over the window alone: no mass to sum.
            if !v.whole(it.shape) {
                if !along.is_empty() {
                    return Err(InkError::msg(format!("{who} ({}): not proven bounded, so its mass is measured over the window alone", it.label)));
                }
                for &r in &of {
                    model.reactions[r].unbounded.push(it.label.clone());
                }
                continue;
            }
            model.measured += 1;
            model.published += usize::from(what.starts_with("published mass"));
            model.res = (model.res.0.min(m.res), model.res.1.max(m.res));
            if let Some(why) = m.unsettled() {
                unsettled.push(format!("{}: {why}", it.label));
            }
            let Some(p) = ink_kernel::rigid::mass_props(&m.moments()) else { continue };
            // Each of its ten moments anywhere in its proven bounds (Measured::error, the farther side). A
            // published mass is its own, exact, and its other moments are the measured ones scaled to it by
            // the measured volume: each one's own reach, times the most that scale can be, and its inertia
            // scaled by as much as the volume is known to (its mass and centre kept).
            let published = what.starts_with("published mass");
            let scale = if published { m.error(0) / p.mass.abs() } else { 0.0 };
            let mut measure: Vec<(f64, V3, M3)> = (usize::from(published)..10)
                .filter_map(|j| {
                    let e = m.error(j) * (1.0 + scale);
                    if !e.is_finite() {
                        return Some((f64::INFINITY, [0.0; 3], [[0.0; 3]; 3]));
                    }
                    let mut s = m.moments();
                    s.m[j] += e;
                    let q = ink_kernel::rigid::mass_props(&s).filter(|_| e > 0.0)?;
                    Some((q.mass - p.mass, std::array::from_fn(|k| q.center[k] - p.center[k]), std::array::from_fn(|r| std::array::from_fn(|c| q.inertia[r][c] - p.inertia[r][c]))))
                })
                .collect();
            if published && scale > 0.0 {
                measure.push((0.0, [0.0; 3], p.inertia.map(|row| row.map(|x| scale * x))));
            }
            // What its first and second moments are known to, by their proven bounds (Measured::error), the
            // second carried from the point they were measured about to its centre.
            let moved = |k: std::ops::Range<usize>| k.map(|j| m.error(j).powi(2)).sum::<f64>().sqrt();
            let off = norm(std::array::from_fn(|k| p.center[k] - m.about[k]));
            let known = [moved(1..4), moved(4..10) + 2.0 * off * moved(1..4)];
            // A body round about its spin has its centre on the axis and
            // the same inertia about every direction across it, exactly:
            // its measures are made so.
            let spin = spin(v, it, rho_node, j.is_some(), &changes, &turning, &eps, &solved);
            // A body only reactions carry is tracked only for the claimed
            // ones carrying it, and never where none is claimed (ink
            // effort's steps are the efforts').
            let serves: Vec<usize> = if along.is_empty() { of.clone() } else { Vec::new() };
            let track = along.is_empty().then_some(claimed);
            let moves_coords = !along.is_empty();
            let mut unknown = [0.0; 2];
            // Whether a body shown repeating about its line repeats every turn (else every half turn).
            let mut every = true;
            let (center, mut inertia, round) = match &spin {
                Some(s) if s.round => {
                    let (c, i) = about_axis(p.mass, p.center, &p.inertia, s.axis);
                    (c, i, [true; 2])
                }
                // It repeats about a line of its own along its spin (a wheel moved a hair off its
                // shaft, a pulley but for its set screw's hole): that much of it is exact, and what is
                // not known of it is what is not known of its small parts.
                Some(s) if let Some((c, i, first, second, all)) = own_round(v, it.shape, rho_node, &m, s.axis).filter(|o| o.2 <= known[0]) => {
                    unknown = [first, second];
                    if unknown != [0.0; 2] {
                        model.made_round.push(MadeRound { body: it.label.clone(), first, second });
                    }
                    every = all;
                    if first > 0.0 && p.mass > 0.0 {
                        model.spin_centres.push(SpinCentre { body: it.label.clone(), within: first / p.mass, of: of.clone() });
                    }
                    (c, i, [true; 2])
                }
                Some(s) => {
                    let (c, i, made) = known_round(p.mass, p.center, &p.inertia, s.axis, known);
                    if known[0] > 0.0 && p.mass > 0.0 {
                        model.spin_centres.push(SpinCentre { body: it.label.clone(), within: known[0] / p.mass, of: of.clone() });
                    }
                    // Only what was made round is unknown: a part read as measured is its measure.
                    unknown = [if made[0] { known[0] } else { 0.0 }, if made[1] { known[1] } else { 0.0 }];
                    if made != [false; 2] {
                        model.made_round.push(MadeRound { body: it.label.clone(), first: unknown[0], second: unknown[1] });
                    }
                    (c, i, made)
                }
                None => (p.center, p.inertia, [false; 2]),
            };
            // A body round about its spin, or repeating about a line of its own along it, keeps its
            // measures' changes round too: what departs from round is what `unknown` (own_round) carries.
            // One made round within what its measure is known to (known_round) keeps round the parts
            // made so, their departure being `unknown`'s (counted once), and the rest as measured.
            if let Some(s) = spin.as_ref().filter(|_| round != [false; 2]) {
                // Linear in each change, about the line's direction alone: its centre's change along it,
                // its second moment's change with no product with it, and the same every way across it
                // where it repeats every turn (round_part).
                let u = s.axis.1;
                for d in measure.iter_mut().filter(|d| d.0.is_finite()) {
                    if round[0] {
                        d.1 = u.map(|x| x * dot(d.1, u));
                    }
                    if round[1] {
                        d.2 = match every {
                            true => round_part(&d.2, u),
                            false => {
                                let du: V3 = mv(&d.2, &u);
                                let uu = dot(u, du);
                                std::array::from_fn(|a| std::array::from_fn(|b| d.2[a][b] - u[a] * du[b] - du[a] * u[b] + 2.0 * uu * u[a] * u[b]))
                            }
                        };
                    }
                }
            }
            if let Some(j) = j {
                spun(p.mass, &center, &mut inertia, v.doc.scalar(j, &v.ctx));
                model.spin_inertia.push(it.label.clone());
            }
            let d = &mut v.doc;
            // A spinning body not shown round joins its link's group, made
            // with the first: the frame the link turns in, and the spin.
            let group = spin.as_ref().filter(|s| !s.round).map(|s| match model.groups.iter().position(|gr| gr.link == s.link) {
                Some(k) => k,
                None => {
                    let frame: Vec<NodeId> = s.frame.r.iter().flatten().chain(&s.frame.t).copied().collect();
                    let dt: Vec<NodeId> = frame.iter().map(|&e| d.derivative(e, Var::T)).collect();
                    let dt2: Vec<NodeId> = dt.iter().map(|&e| d.derivative(e, Var::T)).collect();
                    // Their bounds over a step: none through a solve.
                    let timed = |d: &mut Doc, ns: &[NodeId]| -> [Num; 12] {
                        let v: Vec<Num> = ns.iter().map(|&n| Num::timed(d, n)).collect();
                        v.try_into().unwrap_or_else(|_| unreachable!("a placement's 12 numbers"))
                    };
                    let dq = (0..eps.len())
                        .filter_map(|k| {
                            let dk: Vec<NodeId> = frame.iter().map(|&e| d.derivative(e, Var::Param(eps[k]))).collect();
                            dk.iter().any(|&x| d.g.constant_value(x) != Some(0.0)).then(|| (k, timed(d, &dk)))
                        })
                        .collect();
                    let rad = d.g.constant(s.radians);
                    let angle = d.g.mul(rad, s.angle);
                    let rate = d.derivative(angle, Var::T);
                    let acc = d.derivative(rate, Var::T);
                    let along = (0..eps.len())
                        .filter_map(|k| {
                            let dk = d.derivative(angle, Var::Param(eps[k]));
                            (d.g.constant_value(dk) != Some(0.0)).then(|| (k, Num::timed(d, dk)))
                        })
                        .collect();
                    model.groups.push(Group {
                        link: s.link.clone(),
                        members: Vec::new(),
                        axis: s.turn_axis,
                        frame: timed(d, &frame),
                        dt: timed(d, &dt),
                        dt2: timed(d, &dt2),
                        dq,
                        rate: Num::timed(d, rate),
                        acc: Num::timed(d, acc),
                        along,
                        effort: moves_coords,
                    });
                    model.groups.len() - 1
                }
            });
            if let (Some(k), Some(s)) = (group, &spin) {
                model.groups[k].members.push((it.label.clone(), model.bodies.len(), s.within.0, s.within.1));
            }
            // Its effort follows its placement; a spinning body's, its
            // placement with the spin at rest and the spin's angle, and,
            // where it is not shown round, its placement too, tracked where
            // its group's spin is resolved.
            let mut views = vec![match &spin {
                Some(s) => (s.rest.to_vec(), Some((s.angle, s.radians)), None),
                None => (entries.clone(), None, None),
            }];
            if group.is_some() {
                views.push((entries.clone(), None, group));
            }
            for (placed, angle, group) in views.iter().filter(|_| track != Some(false)) {
                let acc: Vec<NodeId> = placed[9..]
                    .iter()
                    .map(|&e| {
                        let de = d.derivative(e, Var::T);
                        d.derivative(de, Var::T)
                    })
                    .collect();
                tracker.add(d, &placed[..9], Some(1.0), *group, &serves);
                tracker.add(d, &acc, None, *group, &serves);
                // The spin's rate and acceleration, where they are not the
                // coordinates' (tracked as theirs): a spin t drives.
                if let Some((a, _)) = angle.filter(|&(a, _)| reads_t_besides(d, a, &passed)) {
                    let w = d.derivative(a, Var::T);
                    let wd = d.derivative(w, Var::T);
                    tracker.add(d, &[w], None, *group, &serves);
                    tracker.add(d, &[wd], None, *group, &serves);
                }
                // Its partial velocity: a spinning body's, its rest's and
                // its spin's (radians, as a rotation's entries are).
                for &k in &along {
                    let along_k = |d: &mut Doc, ns: &[NodeId]| -> Vec<NodeId> { ns.iter().map(|&e| d.derivative(e, Var::Param(eps[k]))).collect() };
                    let mut dk = along_k(d, placed);
                    if let Some((a, radians)) = *angle {
                        let rad = d.g.constant(radians);
                        dk.extend(along_k(d, &[a]).into_iter().map(|da| d.g.mul(rad, da)));
                    }
                    tracker.partial(d, &dk, *group, &serves);
                }
            }
            model.bodies.push(Body {
                label: it.label.clone(),
                mass: p.mass,
                center,
                inertia,
                at: placement(d, &entries, &along, &eps),
                spin: group.zip(spin.as_ref()).map(|(k, s)| (k, s.axis)),
                unknown,
                measure,
                of,
                apart: false,
                released: Vec::new(),
            });
            body_labels.push(it.label.clone());
        }
        // The file run's loads on what the target places, and the loads its
        // assemblies declare, placed with them (clause::loads_under), each
        // acting at its frame's origin, the point moving with what its path
        // places.
        let mut loads: Vec<(String, [NodeId; 3], [NodeId; 3], Option<ink_lang::Xf>)> =
            v.doc.forms::<crate::Load>().filter(|l| l.run).map(|l| (l.path.clone(), l.force, l.torque, l.at)).collect();
        loads.extend(crate::clause::loads_under(&mut v.doc, target).into_iter().map(|l| (l.path, l.force, l.torque, l.at)));
        for (path, force, torque, at) in loads {
            if !path.starts_with(&format!("{target}/")) {
                continue;
            }
            let on = path_placement(v, &items, &path).ok_or_else(|| InkError::msg(format!("the load on {path}: no instance of that path is placed under {target}")))?;
            // It acts on what moves as PATH's placement does: an assembly
            // whose parts move inside it (a joint, a placement reading t or
            // a coordinate) has no one point it rides.
            let mut inside: Vec<&str> = Vec::new();
            for it in items.iter().filter(|it| it.label.starts_with(&format!("{path}/"))) {
                for l in it.links.iter().filter(|l| l.path.len() > path.len() && l.path.starts_with(&path)) {
                    if l.at.r.iter().flatten().chain(&l.at.t).any(|&e| changes(&v.doc.g, e)) && !inside.contains(&l.path.as_str()) {
                        inside.push(&l.path);
                    }
                }
            }
            if !inside.is_empty() {
                return Err(InkError {
                    hint: Some(
                        "put it on the part its point rides, the instance FRAME is a frame of: (load NAME root/arm/fore (force ...) (at (frame-of root/arm/fore tip)))".into(),
                    ),
                    ..InkError::msg(format!(
                        "the load on {path}: {path} moves inside ({}), so what it acts on does not move as one piece: name the part it acts on",
                        inside.join(", ")
                    ))
                });
            }
            let (on, frame) = (xf_entries(&on), xf_entries(&at.unwrap_or(on)));
            let mut along: Vec<usize> = on.iter().flat_map(|&e| along_of(&v.doc, e)).collect();
            along.sort_unstable();
            along.dedup();
            let of: Vec<usize> = (0..wanted.len()).filter(|&r| carried_paths[r].contains(&path)).collect();
            for &r in &of {
                model.reactions[r].coords.extend(&along);
            }
            read.extend(force.iter().chain(&torque).chain(&frame[9..]).chain(&on));
            let d = &mut v.doc;
            tracker.add(d, &force, None, None, &[]);
            tracker.add(d, &torque, None, None, &[]);
            let dq = partials(d, &on, &along, &eps);
            let g = &d.g;
            model.loads.push(Applied {
                force: nums(g, &force),
                torque: nums(g, &torque),
                at: nums(g, &frame[9..]),
                on: nums(g, &on),
                dq,
                of,
                apart: false,
                path: path.clone(),
                released: Vec::new(),
            });
            load_paths.push(path);
        }
        // Each coupling: its lines in the world (their geometry follows the
        // parts' placements, which the steps track), the side each reaction
        // has it on, and its ends: where each is placed, and the joint it
        // rides.
        let mut span_ends: Vec<(usize, usize, [V3; 2])> = Vec::new();
        for r in &relations {
            let c = r.couple.as_ref().unwrap();
            // A line with a length runs from a point on B to one on A (a belt's span, tangent to tangent). On
            // an end that turns (a pulley, a sprocket) that point is on the wheel: within its width along its
            // axis. One drawn past the wheel's faces by more than a millimetre is a belt beside its pulley.
            for (l, length) in c.lines.iter().filter_map(|l| Some((l, l.length?))) {
                let num = |v: &ModelView, n: NodeId| v.doc.scalar(n, &v.ctx);
                let (q, u, len) = (l.at.map(|n| num(v, n)), l.dir.map(|n| num(v, n)), num(v, length));
                let size = norm(u);
                if !(size > 0.0 && len.is_finite()) {
                    continue;
                }
                for (path, point) in [(&r.b, q), (&r.a, std::array::from_fn(|k| q[k] + u[k] / size * len))] {
                    let Some((jp, _)) = rides(&items, path, target).0 else { continue };
                    let Some(j) = joints.iter().find(|(p, _)| *p == jp).map(|(_, j)| j).filter(|j| j.turn) else { continue };
                    let Some(frame) = path_frame(v, &items, target, &jp) else { continue };
                    let frame = frame.map(|n| num(v, n));
                    let (axis, origin): (V3, V3) = ([frame[j.axis], frame[3 + j.axis], frame[6 + j.axis]], [frame[9], frame[10], frame[11]]);
                    let below = format!("{path}/");
                    let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
                    for it in items.iter().filter(|it| it.label == *path || it.label.starts_with(&below)) {
                        let Some((lo, hi, _)) = v.item_bounds(it) else { continue };
                        for corner in 0..8 {
                            let p: V3 = std::array::from_fn(|k| if corner >> k & 1 == 0 { lo[k] } else { hi[k] });
                            let h = (0..3).map(|k| (p[k] - origin[k]) * axis[k]).sum::<f64>();
                            (low, high) = (low.min(h), high.max(h));
                        }
                    }
                    let h = (0..3).map(|k| (point[k] - origin[k]) * axis[k]).sum::<f64>();
                    if low <= high && (h < low - 1.0 || h > high + 1.0) {
                        let name = r.kind.as_ref().map_or("", |k| k.name.as_str());
                        return Err(InkError::msg(format!(
                            "the line of {}~{} ({name}) meets {path} {:.3} mm along its axis from its joint's origin, and {path} reaches from {:.3} to {:.3} there: the line is \
                             beside the wheel, not on it. Put the line in the wheel's plane (a belt's span at its pulley's belt centre), or the wheel where the line is",
                            r.a, r.b, h, low, high
                        )));
                    }
                }
            }
            let first: usize = model.couplings.iter().map(|c| c.lines.len()).sum();
            // Each span's two ends here, A's then B's: where a carried belt's
            // pieces are matched to it (Model::belts).
            for (i, l) in c.lines.iter().enumerate() {
                let Some(length) = l.length else { continue };
                let num = |v: &ModelView, n: NodeId| v.doc.scalar(n, &v.ctx);
                let (q, u, len) = (l.at.map(|n| num(v, n)), l.dir.map(|n| num(v, n)), num(v, length));
                let size = norm(u);
                if size > 0.0 && len.is_finite() {
                    span_ends.push((model.couplings.len(), i, [std::array::from_fn(|k| q[k] + u[k] / size * len), q]));
                }
            }
            let d = &mut v.doc;
            let mut lines = Vec::new();
            for l in &c.lines {
                read.extend(l.at.iter().chain(&l.dir).chain([&l.pitch]).chain(&l.length));
                let length = l.length.map_or(Num::C(1.0), |n| Num::of(&d.g, n));
                lines.push((nums(&d.g, &l.at), nums(&d.g, &l.dir), Num::of(&d.g, l.pitch), length));
            }
            read.push(c.pretension);
            read.extend(c.efficiency.iter().flat_map(|e| [e.0, e.1]));
            let side = carried_paths.iter().map(|s| i8::from(s.contains(&r.b)) - i8::from(s.contains(&r.a))).collect();
            let carries = carried_paths.iter().map(|s| [s.contains(&r.a), s.contains(&r.b)]).collect();
            let mut ends = Vec::new();
            let joint_input = |v: &mut ModelView, jp: &str, drive: NodeId| -> Input {
                let j = &joints.iter().find(|(p, _)| *p == jp).unwrap().1;
                let (axis, turn) = (j.axis, j.turn);
                let frame = path_frame(v, &items, target, jp).expect("a joint placed under the target");
                let d = &mut v.doc;
                let rate = d.derivative(drive, Var::T);
                let along = (0..eps.len())
                    .filter_map(|k| {
                        let dk = d.derivative(drive, Var::Param(eps[k]));
                        (d.g.constant_value(dk) != Some(0.0)).then(|| (k, Num::of(&d.g, dk)))
                    })
                    .collect();
                Input { joint: jp.to_string(), frame: nums(&d.g, &frame), axis, turn, rate: Num::of(&d.g, rate), along }
            };
            for path in [&r.a, &r.b] {
                let on = xf_entries(&path_placement(v, &items, path).expect("a relation's end placed under the target"));
                let (own, outward) = rides(&items, path, target);
                let input = own.as_ref().map(|(jp, drive)| joint_input(v, jp, *drive));
                let mut outward: Vec<String> = outward.into_iter().filter(|q| own.as_ref().is_none_or(|o| o.0 != *q)).collect();
                outward.sort();
                let moving: Vec<(String, NodeId)> = outward
                    .iter()
                    .filter_map(|q| joints.iter().find(|(p, _)| p == q).map(|(p, j)| (p.clone(), j.drive)))
                    .filter(|(_, drive)| !along_of(&v.doc, *drive).is_empty())
                    .collect();
                let carriers: Vec<Input> = moving.into_iter().map(|(p, drive)| joint_input(v, &p, drive)).collect();
                let d = &mut v.doc;
                let mut along: Vec<usize> = on.iter().flat_map(|&e| along_of(d, e)).collect();
                along.sort_unstable();
                along.dedup();
                let dq = partials(d, &on, &along, &eps);
                ends.push(End { path: path.clone(), on: nums(&d.g, &on), dq, input, carriers });
            }
            let ends: [End; 2] = ends.try_into().unwrap_or_else(|_| unreachable!("a relation's two ends"));
            let g = &v.doc.g;
            let name = format!("{}~{} ({})", r.a, r.b, r.kind.as_ref().map_or("", |k| k.name.as_str()));
            let efficiency = c.efficiency.map(|(e, b)| (Num::of(g, e), Num::of(g, b)));
            model.couplings.push(Coupled { name, lines, first, pretension: Num::of(g, c.pretension), held: c.held, efficiency, side, carries, ends });
        }
        // A carried belt whose every piece lies on a span or wraps a wheel a
        // span meets is separated: its material's inertia is its parts'.
        for (k, ends) in strip_ends.iter().enumerate() {
            let at: Vec<[V3; 2]> = ends.iter().map(|e| e.map(|p| p.map(|n| v.doc.scalar(n, &v.ctx)))).collect();
            let arcs: Vec<bool> = model.strips[k].sweeps.iter().map(Option::is_some).collect();
            if let Some(on) = belt_pieces(&at, &arcs, &span_ends) {
                model.strips[k].on = on;
            }
        }
        model.balances = balancing.iter().map(|(p, turn, axis)| Balancing { reaction: wanted.iter().position(|q| q == p).unwrap(), axis: *axis, turn: *turn }).collect();
        let nc = model.couplings.len();
        let mut trains = Sets((0..nc).collect());
        for b in &model.balances {
            trains.join(&(0..nc).filter(|&c| model.couplings[c].side[b.reaction] != 0).collect::<Vec<_>>());
        }
        for c in 0..nc {
            for d in 0..c {
                if model.couplings[c].ends.iter().any(|e| model.couplings[d].ends.iter().any(|f| f.path == e.path)) {
                    trains.join(&[c, d]);
                }
            }
        }
        model.trains = (0..nc).map(|c| trains.find(c)).collect();
        // Coordinates coupled through a coupling's ends or the balancing reactions of its train
        // load each other: missing mass on one blocks the others.
        for c in &model.couplings {
            let mut coords: Vec<usize> = c.ends.iter().flat_map(|e| e.dq.iter().map(|(k, _)| *k).chain(e.input.iter().flat_map(|i| i.along.iter().map(|(k, _)| *k)))).collect();
            for b in model.balances.iter().filter(|b| c.side[b.reaction] != 0) {
                coords.extend(&model.reactions[b.reaction].coords);
            }
            coupling.join(&coords);
        }
        // The coordinates the drive of the joint each reaction's path rides
        // reads (riding none, a part placed by a formula, those its own
        // placement reads: a linkage's rocker placed by the crank's
        // angle), and what supports the path.
        let rode: Vec<(Vec<usize>, HashSet<String>)> = wanted
            .iter()
            .map(|p| {
                let (drive, supports) = rides(&items, p, target);
                let own = match drive {
                    Some((_, d)) => along_of(&v.doc, d),
                    None => path_placement(v, &items, p).map_or_else(Vec::new, |x| xf_entries(&x).iter().flat_map(|&e| along_of(&v.doc, e)).collect()),
                };
                (own, supports)
            })
            .collect();
        // Where couplings are, what an actuated joint or a balancing one
        // carries is accounted for: its pull on the rest goes through that
        // joint or its couplings. What else moves (a linkage's coupler on a
        // driven crank, a belt's own mass) pushes on what it touches
        // unseparated.
        let mut accounted: HashSet<String> = HashSet::new();
        if !model.couplings.is_empty() {
            let driven: Vec<String> = model.actuators.iter().map(|a| a.path.clone()).collect();
            accounted.extend(carried(&items, &driven).0.into_iter().flatten());
            accounted.extend(model.balances.iter().flat_map(|b| carried_paths[b.reaction].iter().cloned()));
        }
        let separated: Vec<String> = model.strips.iter().filter(|s| !s.on.is_empty()).map(|s| s.label.clone()).collect();
        let mut belts: Vec<Vec<String>> = vec![Vec::new(); model.reactions.len()];
        for s in model.strips.iter().filter(|s| !s.on.is_empty()) {
            for rides in &s.on {
                let (Rides::Span(c, _) | Rides::Wrap(c, _, _)) = *rides;
                for (r, b) in belts.iter_mut().enumerate() {
                    if model.couplings[c].carries[r].iter().any(|&x| x) && !s.of.contains(&r) && !b.contains(&s.label) {
                        b.push(s.label.clone());
                    }
                }
            }
        }
        for (k, r) in model.reactions.iter_mut().enumerate() {
            r.belts = std::mem::take(&mut belts[k]);
            r.coords.sort_unstable();
            r.coords.dedup();
            let (own, supports) = &rode[k];
            // A coupling reaching what it carries separates what moves with
            // it through accounted joints, and so does its riding one whose
            // couplings reach none of it (a servo bolted to a belted carriage:
            // the belt pulls the carriage, not the servo): what else moves
            // with it, and what is placed partly on it, stays unseparated.
            let coupled = model.couplings.iter().any(|c| c.side[k] != 0) || accounted.contains(&wanted[k]);
            // A separated belt's push is on the parts it lies on (belts),
            // so on what carries them, and on nothing else.
            let apart = |l: &String| (coupled && accounted.contains(l)) || carried_paths[k].contains(l) || supports.contains(l) || separated.contains(l);
            let mut with: Vec<String> = along_coords.iter().filter(|(l, a)| !apart(l) && a.iter().any(|c| own.contains(c))).map(|(l, _)| l.clone()).collect();
            for l in &partly[k] {
                if !with.contains(l) {
                    with.push(l.clone());
                }
            }
            r.moves_with = with;
        }
        // What the balances giving the tensions leave unseparated, its
        // effort summed apart too: a coupling's loss rests on tensions that
        // leave it out.
        for b in &model.balances {
            for l in &model.reactions[b.reaction].moves_with {
                if !model.unseparated.contains(l) {
                    model.unseparated.push(l.clone());
                }
            }
        }
        for (b, l) in model.bodies.iter_mut().zip(&body_labels) {
            b.apart = model.unseparated.contains(l);
        }
        for (s, l) in model.strips.iter_mut().zip(&strip_labels) {
            s.apart = model.unseparated.contains(l);
        }
        for (l, p) in model.loads.iter_mut().zip(&load_paths) {
            l.apart = model.unseparated.iter().any(|u| u == p || u.starts_with(&format!("{p}/")));
        }
        model.behind =
            (0..model.groups.len()).map(|g| model.departure(&model.groups[g]) != (0.0, 0.0, 0.0) && model.balances.iter().any(|b| model.carries(b.reaction, g))).collect();
        model.leaves = deforming.iter().filter(|(_, _, mass)| *mass != Some(0.0)).map(|(l, reads, _)| (l.clone(), reads.clone())).collect();
        // Each actuator rests on the coordinates coupled with those it
        // moves, and leaves out the deforming bodies moving with them.
        model.sets = (0..eps.len()).map(|k| coupling.find(k)).collect();
        for a in &mut model.actuators {
            let sets: Vec<usize> = a.coupled.iter().map(|&k| coupling.find(k)).collect();
            a.coupled = (0..eps.len()).filter(|&k| sets.contains(&coupling.find(k))).collect();
            a.leaves_out = deforming
                .iter()
                .filter(|(_, reads, mass)| reads.iter().any(|k| a.coupled.contains(k)) && *mass != Some(0.0))
                .map(|(l, _, mass)| match mass {
                    Some(m) => format!("{l} ({m:.4} g)"),
                    None => format!("{l} (no mass)"),
                })
                .collect();
        }
        // An actuator held at a pose a body's shape follows leaves that
        // body out: its effort is not found.
        for (a, sizes) in shaped {
            let follows: Vec<&ink_lang::Item> = items.iter().filter(|it| v.doc.g.depends_on(it.shape, |x| matches!(x, Var::Param(k) if sizes.contains(&k)))).collect();
            for it in follows {
                let attrs = v.doc.body(&it.instance.as_ref().map_or(String::new(), |i| i.0.clone())).map(|b| b.attrs.clone()).unwrap_or_default();
                let who = v.doc.named(it);
                let mass = match ink_lang::materials::mass_density(&mut v.doc, &who, &attrs, None)? {
                    Some((rho, _)) => {
                        v.want(&[rho])?;
                        v.moments(it.shape, Some(rho), None).map(|m| m.sums[0])
                    }
                    None => None,
                };
                model.actuators[a].leaves_out.push(match mass {
                    Some(m) => format!("{} ({m:.4} g)", it.label),
                    None => format!("{} (no mass)", it.label),
                });
            }
        }
        // The params it reads, by name: the document it answers for
        // numbers its params otherwise (no coordinates).
        let own: Vec<usize> = coords.iter().map(|c| c.param).collect();
        model.reads = v.doc.params_read(&read).into_iter().filter(|k| !own.contains(k)).map(|k| v.doc.params[k].name.clone()).collect();
        model.reads_carried = v.doc.params_read(&read_carried).into_iter().filter(|k| !own.contains(k)).map(|k| v.doc.params[k].name.clone()).collect();
        let derivatives = joint_derivatives(v, &eps, joints);
        // The rates made since the view's context have their values at its
        // moment there (solved again at each moment, as the solves on t).
        model.params = v.ctx.params.clone();
        model.params.extend(v.doc.params[model.params.len()..].iter().map(|p| p.default));
        // Solved again at each moment: those the mechanism moves through,
        // and their rates (made since `through`), on t or reading a
        // coordinate as their roots do.
        let eps_read = |i: usize| v.doc.solves_under(&v.doc.solves[i].owns).iter().any(|&j| v.doc.solves[j].reads.iter().any(|k| eps.contains(&(*k as u32))));
        model.solves = (0..v.doc.solves.len()).filter(|&i| (v.doc.solves[i].on_t || eps_read(i)) && !v.doc.solves[i].on_demand).collect();
        model.follow = Some(v.doc.follow());
        model.tracked = tracker.finish(&mut v.doc, &model, v.ctx.t);
        Ok((model, derivatives))
    }
}
