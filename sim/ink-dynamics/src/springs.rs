//! A mechanism's springs in a simulation (simulate.rs): the model its modes
//! linearise (modes.rs), stepped forward. Each joint a stiff coupling line's
//! end rests on (the one it rides and those supporting it), and each joint
//! stating `(stiffness K)`, is released: it gets a coordinate of its own,
//! its deviation from what the coordinates command it to (a small turn,
//! radians, or slide, mm, beside what they move it by), but an actuator
//! commanded (its command moves its coordinate: the coordinate is its
//! motion) or stating no stiffness (the written motion holds it, its loop
//! taken as rigid, as an effort takes it). The bodies it carries are
//! placed where the coordinates put them, and its deviation moves them by
//! its unit motion there: their
//! inertia along it is the modes' (their partial velocities), and the
//! effort it needs is what those bodies' motion and the loads on them need
//! along it, as along a coordinate (Model::compute_in). A stiff line
//! carries its pretension plus its stiffness times its stretch, its ends'
//! points drawing apart along it as the deviations move them (the modes'
//! stretch per unit of each), and pulls nothing below 0 (a line that would
//! push goes slack, as a coupling's does: couple.rs); a joint's stiffness
//! holds its deviation at 0, at what its drive commands (a stepper's rotor
//! at its commanded angle). Its coupling is no rigid relation there: the
//! springs carry it.

use super::*;

/// A joint released for the springs: its path, kind, frame in the world
/// (its axis the frame's column `axis`), its drive's value, rate and
/// acceleration as the coordinates command it (deg a turn, mm a slide), its
/// drive's rate along each coordinate (radians a turn), and its stiffness
/// with where it is stated.
pub(crate) struct Sprung {
    pub path: String,
    pub turn: bool,
    axis: usize,
    frame: [Num; 12],
    drive: [Num; 3],
    along: Vec<(usize, Num)>,
    stiffness: Option<(Num, String)>,
}

/// A coupling line stating a stiffness: its name, point, direction, pitch,
/// stiffness and its coupling's pretension, where it is stated, and each
/// end's placement in the world with the released joints carrying it.
pub(crate) struct Line {
    pub name: String,
    pair: (String, String),
    p: [Num; 3],
    u: [Num; 3],
    h: Num,
    k: Num,
    t0: Num,
    pub site: String,
    /// The released joints carrying each end.
    ends: [Vec<usize>; 2],
}

/// The springs of a model of states (Model::compile_springs).
pub(crate) struct Springs {
    pub joints: Vec<Sprung>,
    pub lines: Vec<Line>,
    /// What the actuators that keep their coordinates carry (commanded, or
    /// stating no stiffness: the gauge the deviations are measured from).
    kept: Vec<String>,
}

/// The springs at a state: each line's stretch per unit of each released
/// joint's deviation (N mm of work per N of tension), stiffness and
/// pretension; each joint's stiffness (0 where none is stated), its drive's
/// value, rate and acceleration (radians a turn), its heft (what it
/// carries: g mm^2 of their inertia about its axis's origin a turn, their
/// mass g a slide), against which its inertia along itself is none, and
/// the mass it carries (g).
pub struct At {
    pub g: Vec<Vec<f64>>,
    pub kl: Vec<f64>,
    pub t0: Vec<f64>,
    pub kj: Vec<f64>,
    pub drive: Vec<[f64; 3]>,
    pub heft: Vec<f64>,
    pub mass: Vec<f64>,
}

/// What a spring's stretch along a motion reads as none below: this share
/// of what its stretch's terms would sum to, its stretch per unit of each
/// coordinate times the motion's part along it, each in units scaled alike
/// (1e-12 of a stiffness, squared). Each spring is judged on its own, so a
/// stiff one elsewhere hides no soft one's hold.
pub(crate) const HELD: f64 = 1e-6;

/// Whether spring row h (its stretch per unit of each coordinate) stretches
/// along motion z, both in the same scaled coordinates: past HELD of |h| |z|,
/// which a motion it does not reach (its rounding, of z and of h) is not.
pub(crate) fn stretches(h: &[f64], z: &[f64]) -> bool {
    let size = h.iter().map(|x| x * x).sum::<f64>().sqrt() * z.iter().map(|x| x * x).sum::<f64>().sqrt();
    dot_n(h, z).abs() > HELD * size
}

/// A line pulls on extension only. At the threshold, equal lines with
/// opposite stretch have a single two-sided tangent: one line's stiffness.
/// For an unmatched line keep its tensile-side tangent for the solver's
/// iteration, and name it so a modal answer can decline that linearisation.
pub(super) fn line_tangents(g: &[Vec<f64>], lines: impl Iterator<Item = (f64, f64)>) -> (Vec<f64>, Option<usize>) {
    // Compare each stretch with its counterpart in the same joint units.
    // Roundoff in a long moment arm cannot hide a different slide component.
    let same = |a: f64, b: f64| (a - b).abs() <= 64.0 * f64::EPSILON * a.abs().max(b.abs());
    let mut tangent = Vec::new();
    let mut unmatched: Vec<usize> = Vec::new();
    for (i, (k, tension)) in lines.enumerate() {
        let mut active = tension >= 0.0;
        if tension == 0.0 && k > 0.0 && g[i].iter().any(|&x| x != 0.0) {
            if let Some(j) = unmatched.iter().position(|&j| same(tangent[j], k) && g[j].iter().zip(&g[i]).all(|(a, b)| same(*a, -*b))) {
                unmatched.swap_remove(j);
                active = false;
            } else {
                unmatched.push(i);
            }
        }
        tangent.push(if active { k } else { 0.0 });
    }
    (tangent, unmatched.first().copied())
}

impl At {
    /// Each line's tension at deviations `d` (all released joints'): its
    /// pretension plus its stiffness times its stretch, at least 0.
    pub fn tensions(&self, d: &[f64]) -> Vec<f64> {
        (0..self.kl.len()).map(|l| (self.t0[l] + self.kl[l] * dot_n(&self.g[l], d)).max(0.0)).collect()
    }

    /// What the springs exert along each released joint at deviations `d`
    /// (N a slide, N mm per radian a turn): each line's tension beyond its
    /// pretension (the pretensions balanced at rest), each joint's
    /// stiffness times its deviation, both against it.
    pub fn force(&self, d: &[f64]) -> Vec<f64> {
        let mut out: Vec<f64> = (0..d.len()).map(|j| -self.kj[j] * d[j]).collect();
        for (l, t) in self.tensions(d).into_iter().enumerate() {
            let pull = t - self.t0[l];
            for (o, g) in out.iter_mut().zip(&self.g[l]) {
                *o -= pull * g;
            }
        }
        out
    }

    /// Stiffness for a solver iteration, using a one-sided tangent where
    /// the force has no two-sided derivative.
    pub fn stiffness(&self, d: &[f64]) -> Vec<Vec<f64>> {
        self.tangent(d).0
    }

    /// The stiffness matrix and whether it is a two-sided linearisation.
    pub(crate) fn tangent(&self, d: &[f64]) -> (Vec<Vec<f64>>, bool) {
        let n = d.len();
        let (springs, two_sided) = self.springs(d);
        let mut k = vec![vec![0.0; n]; n];
        for (g, s) in &springs {
            (0..n).for_each(|i| (0..n).for_each(|j| k[i][j] += s * g[i] * g[j]));
        }
        (k, two_sided)
    }

    /// Each spring pulling at deviations `d`, its stretch per unit of each
    /// released joint's deviation and its tangent stiffness (a joint's
    /// stiffness, its own deviation; a line's, where it pulls), and whether
    /// that is a two-sided linearisation.
    pub(crate) fn springs(&self, d: &[f64]) -> (Vec<(Vec<f64>, f64)>, bool) {
        let n = d.len();
        let mut out: Vec<(Vec<f64>, f64)> = (0..n).filter(|&j| self.kj[j] > 0.0).map(|j| ((0..n).map(|i| f64::from(u8::from(i == j))).collect(), self.kj[j])).collect();
        let (lines, unmatched) = line_tangents(&self.g, self.kl.iter().enumerate().map(|(l, &s)| (s, self.t0[l] + s * dot_n(&self.g[l], d))));
        out.extend(self.g.iter().zip(lines).filter(|(_, s)| *s > 0.0).map(|(g, s)| (g.clone(), s)));
        (out, unmatched.is_none())
    }
}

/// The joints a part's motion rests on: the one it rides and every one
/// supporting it, up to the target (a pulley's turn and the carriage
/// carrying it), those with joints (`jointed`).
pub(crate) fn resting(items: &[ink_lang::Item], jointed: &HashSet<&str>, path: &str, target: &str) -> Vec<String> {
    let (own, supports) = rides(items, path, target);
    let mut out: Vec<String> = own.map(|(p, _)| p).into_iter().collect();
    let mut more: Vec<String> = supports.into_iter().filter(|p| jointed.contains(p.as_str()) && !out.contains(p)).collect();
    more.sort();
    out.extend(more);
    out
}

/// Whether a coupling's line states a stiffness.
fn stiff(c: &crate::coupling::Couple) -> bool {
    c.lines.iter().any(|l| l.stiffness.is_some())
}

/// The joints a simulation releases for its springs under `target`: each a
/// stiff line's end rests on, then each stating a stiffness, but those of
/// `commanded` actuators, and actuators stating no stiffness: the written
/// motion holds those (a servo's or a motor's loop taken as rigid, as an
/// effort takes it). Err: a coupling with stiff and rigid lines, or a rigid
/// line reaching a released joint (a constraint between deviations a
/// simulation does not hold yet).
pub(crate) fn released(items: &[ink_lang::Item], relations: &[crate::coupling::Related], target: &str, commanded: &[String]) -> Result<Vec<String>, InkError> {
    let jointed: HashSet<&str> = items.iter().flat_map(|it| &it.links).filter(|l| l.joint.is_some()).map(|l| l.path.as_str()).collect();
    let held: Vec<String> = items
        .iter()
        .flat_map(|it| &it.links)
        .filter(|l| l.joint.as_ref().is_some_and(|j| crate::joint::effort(j).is_some() && crate::joint::stiffness(j).is_none()))
        .map(|l| l.path.clone())
        .collect();
    let mut out: Vec<String> = Vec::new();
    let mut add = |p: String| {
        if !out.contains(&p) && !commanded.contains(&p) && !held.contains(&p) {
            out.push(p);
        }
    };
    for r in relations {
        let Some(c) = &r.couple else { continue };
        if !stiff(c) {
            continue;
        }
        if c.lines.iter().any(|l| l.stiffness.is_none()) {
            return Err(InkError::msg(format!(
                "the coupling {}~{} states a stiffness on some of its lines and not others: in a simulation each of its lines is a spring or none is",
                r.a, r.b
            )));
        }
        for e in resting(items, &jointed, &r.a, target).into_iter().chain(resting(items, &jointed, &r.b, target)) {
            add(e);
        }
    }
    for l in items.iter().flat_map(|it| &it.links) {
        if l.joint.as_ref().is_some_and(|j| crate::joint::stiffness(j).is_some()) {
            add(l.path.clone());
        }
    }
    let (sets, _) = carried(items, &out);
    for r in relations {
        let Some(c) = &r.couple else { continue };
        if stiff(c) {
            continue;
        }
        if let Some(i) = (0..out.len()).find(|&i| sets[i].contains(&r.a) != sets[i].contains(&r.b)) {
            let kind = r.kind.as_ref().map_or("couple", |k| k.name.as_str());
            return Err(InkError {
                hint: Some("state the coupling's stiffness ((tensile-stiffness EA) on a belt; (vec P U H L K) on a line), or simulate without the springs".into()),
                ..InkError::msg(format!(
                    "{}~{} ({kind}) ties {} rigidly to what the springs move apart: a simulation does not yet hold a rigid line between joints released \
                     for springs (the modes do)",
                    r.a, r.b, out[i]
                ))
            });
        }
    }
    Ok(out)
}

impl Model {
    /// The mechanism under `target` in a view of a load with states, its
    /// springs released (springs::released): the stiff couplings are no
    /// relations of it, and each released joint carries its bodies and
    /// loads along its deviation (Model::compute_in). Err: as
    /// Model::compile, or as `released`.
    pub(crate) fn compile_springs(v: &mut ModelView, target: &str, gravity: [f64; 3]) -> Result<Model, InkError> {
        let items = v.doc.items(target).map_err(InkError::msg)?;
        let relations: Vec<crate::coupling::Related> = crate::coupling::related(&mut v.doc, target).into_iter().filter(|r| r.couple.is_some()).collect();
        let commanded: Vec<String> = v.doc.forms::<crate::Command>().filter(|c| c.run).map(|c| c.path.clone()).collect();
        let released = released(&items, &relations, target, &commanded)?;
        let (mut model, (joints, sprung)) = Model::compile_in(v, target, gravity, &released, false, true, false, |v, eps, joints| {
            let mut sprung = Vec::new();
            for p in &released {
                let j = &joints.iter().find(|(q, _)| q == p).expect("a released joint").1;
                let frame = path_frame(v, &items, target, p).expect("a joint placed under the target");
                let d = &mut v.doc;
                let rate = d.derivative(j.drive, Var::T);
                let acc = d.derivative(rate, Var::T);
                let unit = if j.turn { std::f64::consts::PI / 180.0 } else { 1.0 };
                let along = (0..eps.len())
                    .filter_map(|k| {
                        let dk = d.derivative(j.drive, Var::Param(eps[k]));
                        let u = d.g.constant(unit);
                        let dk = d.g.mul(u, dk);
                        (d.g.constant_value(dk) != Some(0.0)).then(|| (k, Num::of(&d.g, dk)))
                    })
                    .collect();
                let stiffness = crate::joint::stiffness(j).map(|s| (Num::of(&d.g, s.k), s.site.clone()));
                sprung.push(Sprung {
                    path: p.clone(),
                    turn: j.turn,
                    axis: j.axis,
                    frame: nums(&d.g, &frame),
                    drive: [j.drive, rate, acc].map(|n| Num::of(&d.g, n)),
                    along,
                    stiffness,
                });
            }
            (joints, sprung)
        })?;
        if released.is_empty() {
            return Ok(model);
        }
        let (sets, _) = carried(&items, &released);
        let carriers = |path: &str| -> Vec<usize> { (0..released.len()).filter(|&i| sets[i].contains(path)).collect() };
        for b in &mut model.bodies {
            b.released = carriers(&b.label);
        }
        for l in &mut model.loads {
            l.released = carriers(&l.path);
        }
        // What the springs read besides the model's own: the stiffnesses,
        // and the stiff couplings' lines and pretensions.
        let mut read: Vec<NodeId> = joints.iter().filter(|(p, _)| released.contains(p)).filter_map(|(_, j)| crate::joint::stiffness(j).map(|s| s.k)).collect();
        let mut lines = Vec::new();
        for r in relations.iter().filter(|r| stiff(r.couple.as_ref().unwrap())) {
            read.extend(r.couple.as_ref().unwrap().fields());
            let c = r.couple.as_ref().unwrap();
            let name = format!("{}~{} ({})", r.a, r.b, r.kind.as_ref().map_or("", |k| k.name.as_str()));
            let ends = [&r.a, &r.b].map(|path| carriers(path));
            for (i, l) in c.lines.iter().enumerate() {
                let g = &v.doc.g;
                lines.push(Line {
                    name: format!("line {} of {name}", i + 1),
                    pair: (r.a.clone(), r.b.clone()),
                    p: nums(g, &l.at),
                    u: nums(g, &l.dir),
                    h: Num::of(g, l.pitch),
                    k: Num::of(g, l.stiffness.expect("a stiff line")),
                    t0: Num::of(g, c.pretension),
                    site: c.site.clone(),
                    ends: ends.clone(),
                });
            }
        }
        let own: Vec<usize> = v.doc.coords.iter().map(|c| c.param).collect();
        for k in v.doc.params_read(&read) {
            let name = &v.doc.params[k].name;
            if !own.contains(&k) && !model.reads.contains(name) {
                model.reads.push(name.clone());
            }
        }
        let kept: Vec<String> = joints.iter().filter(|(p, j)| crate::joint::effort(j).is_some() && !released.contains(p)).map(|(p, _)| p.clone()).collect();
        let kept = carried(&items, &kept).0.into_iter().flatten().collect();
        // A stiffness is at least 0 and finite, as a mode reads it (modes): a
        // negative spring pushes its stretch further, and none is a number.
        let ctx = Ctx::new(0.0, model.params.clone());
        let wrong = |k: f64| !(k >= 0.0 && k.is_finite());
        let bad =
            lines.iter().map(|l| (l.name.clone(), l.k.at(&ctx))).chain(sprung.iter().filter_map(|j| j.stiffness.as_ref().map(|(k, _)| (format!("joint {}", j.path), k.at(&ctx)))));
        if let Some((name, k)) = bad.into_iter().find(|&(_, k)| wrong(k)) {
            return Err(InkError::msg(format!("the stiffness of {name} is {k}: a stiffness is at least 0 and finite")));
        }
        model.springs = Some(Springs { joints: sprung, lines, kept });
        Ok(model)
    }

    /// Joint i's unit motion at a placement (R row major, then T): its
    /// rate dR, dT per radian of a turn about its axis, per mm of a slide
    /// along it; `frames` each released joint's axis and origin.
    pub(crate) fn unit_motion(&self, frames: &[(V3, V3, bool)], i: usize, r: &[f64], p: &[f64]) -> [f64; 12] {
        let (a, o, turn) = frames[i];
        let mut out = [0.0; 12];
        if turn {
            let w: M3 = [[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]];
            let dr = mul(&w, &rot(r));
            for x in 0..3 {
                for y in 0..3 {
                    out[3 * x + y] = dr[x][y];
                }
            }
            let arm: V3 = std::array::from_fn(|k| p[k] - o[k]);
            out[9..].copy_from_slice(&mv(&w, &arm));
        } else {
            out[9..].copy_from_slice(&a);
        }
        out
    }

    /// Each released joint's axis, origin and kind at `ctx`.
    pub(crate) fn spring_frames(&self, ctx: &Ctx) -> Vec<(V3, V3, bool)> {
        let Some(s) = &self.springs else { return Vec::new() };
        s.joints
            .iter()
            .map(|j| {
                let f: Vec<f64> = j.frame.iter().map(|n| n.at(ctx)).collect();
                ([f[j.axis], f[3 + j.axis], f[6 + j.axis]], [f[9], f[10], f[11]], j.turn)
            })
            .collect()
    }

    /// The springs at `ctx` (At), given each joint's heft and mass.
    pub(crate) fn springs_at(&self, ctx: &Ctx, frames: &[(V3, V3, bool)], heft: Vec<f64>, mass: Vec<f64>) -> At {
        let s = self.springs.as_ref().expect("a model with springs");
        let nr = s.joints.len();
        let mut at = At { g: Vec::new(), kl: Vec::new(), t0: Vec::new(), kj: vec![0.0; nr], drive: Vec::new(), heft, mass };
        for l in &s.lines {
            let p: V3 = l.p.each_ref().map(|e| e.at(ctx));
            let u: V3 = l.u.each_ref().map(|e| e.at(ctx));
            let len = norm(u);
            let u = u.map(|x| x / len);
            let h = l.h.at(ctx);
            let mut g = vec![0.0; nr];
            // Each end's point on the line moving along it, and turning
            // about it times its pitch, per unit of each joint carrying it.
            for (carriers, sign) in l.ends.iter().zip([1.0, -1.0]) {
                for &i in carriers {
                    let (a, o, turn) = frames[i];
                    let (vel, w) = if turn { (cross(a, std::array::from_fn(|k| p[k] - o[k])), a) } else { (a, [0.0; 3]) };
                    g[i] += sign * (dot(u, vel) + h * dot(u, w));
                }
            }
            at.g.push(g);
            at.kl.push(l.k.at(ctx));
            at.t0.push(l.t0.at(ctx));
        }
        for (i, j) in s.joints.iter().enumerate() {
            if let Some((k, _)) = &j.stiffness {
                at.kj[i] = k.at(ctx);
            }
            let unit = if j.turn { std::f64::consts::PI / 180.0 } else { 1.0 };
            at.drive.push(j.drive.each_ref().map(|n| unit * n.at(ctx)));
        }
        at
    }

    /// The bodies a written motion moves otherwise than through the
    /// released joints carrying them, where it commands one of those too:
    /// tied to them only by the drives' expressions, which releasing them
    /// cuts (the modes refuse it alike), each with that coordinate's name.
    /// A coordinate held at a pose moves nothing, and a simulated one is the
    /// motion of what it moves; an actuator keeping its coordinate (commanded,
    /// or stating no stiffness) carries what it carries on it.
    pub(crate) fn springs_cut(&self, t: f64, q: &[[f64; 3]]) -> Result<Vec<(String, String)>, String> {
        let Some(s) = &self.springs else { return Ok(Vec::new()) };
        let ctx = &self.ctx_with(t, self.state_params(t, q))?;
        let frames = self.spring_frames(ctx);
        let vals = |n: &[Num]| -> Vec<f64> { n.iter().map(|e| e.at(ctx)).collect() };
        let mut cut = Vec::new();
        for b in self.bodies.iter().filter(|b| !s.kept.contains(&b.label)) {
            let (r, p) = (vals(&b.at.r), vals(&b.at.p));
            for (k, d) in &b.at.dq {
                if !self.coords[*k].moves() || self.simulated.iter().any(|x| x.0 == *k) || !s.joints.iter().any(|j| j.along.iter().any(|(l, _)| l == k)) {
                    continue;
                }
                let mut e = vals(d);
                let mut size: f64 = e.iter().map(|x| x.abs()).sum();
                for &i in &b.released {
                    let Some((_, a)) = s.joints[i].along.iter().find(|(l, _)| l == k) else { continue };
                    let a = a.at(ctx);
                    let u = self.unit_motion(&frames, i, &r, &p);
                    for (x, y) in e.iter_mut().zip(u) {
                        *x -= a * y;
                    }
                    size += a.abs() * u.iter().map(|x| x.abs()).sum::<f64>();
                }
                if e.iter().map(|x| x.abs()).sum::<f64>() > 1e-9 * size {
                    cut.push((b.label.clone(), self.coords[*k].name.clone()));
                }
            }
        }
        Ok(cut)
    }

    /// The released joints' paths, kinds and stiffness's sites, in order
    /// (none without springs).
    pub fn released(&self) -> Vec<(&str, bool)> {
        self.springs.as_ref().map_or_else(Vec::new, |s| s.joints.iter().map(|j| (j.path.as_str(), j.turn)).collect())
    }

    /// The spring lines' names, and where each is stated.
    pub fn spring_lines(&self) -> Vec<(&str, &str)> {
        self.springs.as_ref().map_or_else(Vec::new, |s| s.lines.iter().map(|l| (l.name.as_str(), l.site.as_str())).collect())
    }

    /// Released joint j's stiffness (0 where it states none: N/mm a slide,
    /// N mm per radian a turn; of params alone).
    pub fn released_stiffness(&self, j: usize) -> f64 {
        let s = self.springs.as_ref().expect("a model with springs");
        s.joints[j].stiffness.as_ref().map_or(0.0, |(k, _)| k.at(&Ctx::new(0.0, self.params.clone())))
    }

    /// The joints released for springs whose parts move along the
    /// coordinates actuator `path`'s effort rests on: its load runs through
    /// them, so its effort is not the rigid model's.
    pub fn springs_behind(&self, path: &str) -> Vec<String> {
        let (Some(s), Some(a)) = (&self.springs, self.actuators.iter().find(|a| a.path == path)) else { return Vec::new() };
        let mut out: Vec<String> = Vec::new();
        for b in self.bodies.iter().filter(|b| b.at.dq.iter().any(|(k, _)| a.coupled.contains(k))) {
            for &j in &b.released {
                if !out.contains(&s.joints[j].path) {
                    out.push(s.joints[j].path.clone());
                }
            }
        }
        out
    }

    /// Whether a coupling of springs relates a and b.
    pub fn sprung(&self, a: &str, b: &str) -> bool {
        self.springs.as_ref().is_some_and(|s| s.lines.iter().any(|l| (l.pair.0 == a && l.pair.1 == b) || (l.pair.0 == b && l.pair.1 == a)))
    }

    /// The bodies the released joints carry.
    pub fn sprung_bodies(&self) -> Vec<String> {
        self.bodies.iter().filter(|b| !b.released.is_empty()).map(|b| b.label.clone()).collect()
    }

    /// Which of `labels` reaction r carries.
    pub fn carried_by(&self, r: usize, labels: &[String]) -> Vec<String> {
        self.bodies.iter().filter(|b| b.of.contains(&r) && labels.contains(&b.label)).map(|b| b.label.clone()).collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::simulate::Simulation;
    use ink_lang::Doc;
    use std::path::Path;
    use std::sync::Arc;

    static REGISTRY: ink_lang::Registry = ink_lang::Registry::new(&[&crate::Dynamics]);

    fn doc(src: &str) -> Doc {
        Doc::load_str_in(&REGISTRY, src, "springs.ink", Path::new("."), &[]).unwrap_or_else(|e| panic!("{} at {}", e.message, e.line))
    }

    /// The file's simulation (it reads a following error).
    fn sim(d: &Doc) -> Arc<Simulation> {
        crate::simulation(d).expect("a simulation").unwrap_or_else(|e| panic!("{e}"))
    }

    fn hz(w2: f64) -> f64 {
        w2.sqrt() / (2.0 * std::f64::consts::PI)
    }

    /// A 1 kg block on a slide its (stiffness 100N/mm) holds at its command,
    /// the command stepping from 5 to 0 mm at 0.01 s: released from a
    /// stretch of 5 mm, it rings at sqrt(k / m), 50.33 Hz, the frequency
    /// its modes give, 5 cos(omega (t - 0.01)) after the step, its energy
    /// (k x^2 + m v^2) / 2 held.
    #[test]
    fn a_mass_on_a_spring_rings_at_its_mode_and_keeps_its_energy() {
        let d = doc(r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.1s)
            (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly (axis (c 0)) (instance b b (at) (slide x c) (stiffness 100N/mm)))
            (assembly root (instance a (axis (select (- (/ t 1s) 0.01) 5 0))))
            (springs) (def lag (following-error root/a/b))"#);
        let s = sim(&d);
        let (i, turn) = s.released("root/a/b").unwrap();
        assert!(!turn);
        let modes = crate::model::modes::modes(&d.reloader().unwrap(), "root", 0.02).unwrap();
        let mass = modes.model.bodies.iter().find(|b| b.label == "root/a/b").unwrap().mass;
        let w = (1e6 * 100.0 / mass).sqrt();
        assert!(modes.modes.len() == 1 && (modes.modes[0].frequency / (w / (2.0 * std::f64::consts::PI)) - 1.0).abs() < 1e-9, "{}", modes.modes[0].frequency);
        assert!((modes.modes[0].frequency - 50.329).abs() < 1e-3);
        // Still before the step, then the cosine.
        for t in [0.0, 0.005, 0.0099] {
            let x = s.at(i, t).unwrap();
            assert!(x[0].abs() < 1e-12 && x[1].abs() < 1e-9, "at {t}: {x:?}");
        }
        let energy = |x: [f64; 5]| 0.5 * 100.0 * x[0] * x[0] + 0.5 * 1e-6 * mass * x[1] * x[1];
        let mut worst = 0.0f64;
        for k in 1..=90 {
            let t = 0.01 + k as f64 * 0.001;
            let x = s.at(i, t).unwrap();
            let want = 5.0 * (w * (t - 0.01)).cos();
            worst = worst.max((x[0] - want).abs());
            assert!((energy(x) / 1250.0 - 1.0).abs() < 3e-5, "its energy at {t}: {}", energy(x));
        }
        assert!(worst < 1e-4, "{worst} mm from 5 cos(omega t)");
        // Its steps follow the ringing (some 25 a period at its tolerance,
        // 1e-6 of its 5 mm swing).
        assert!(s.reach(0.1, |tr| tr.t.len()) < 500);
        // Its half period between the following error's turning points.
        let turns = s.moments(i, (0.011, 0.1));
        let half = (turns[turns.len() - 1] - turns[0]) / (turns.len() - 1) as f64;
        assert!((0.5 / half / modes.modes[0].frequency - 1.0).abs() < 1e-6, "{} Hz against {}", 0.5 / half, modes.modes[0].frequency);
    }

    /// modes.rs's printer axis (a stepper's rotor and pulley, 5478.7 g mm2,
    /// on its holding stiffness 22500 N mm/rad; an idler without inertia;
    /// the carriage between spans of EA 16903 N), its command stepping one
    /// full step, 1.8 deg (0.2 mm of belt), at 2 ms. By hand, in carriage
    /// mm (the rotor's rim x_r = R theta): M = diag(J / R^2, m), K = [[k_h /
    /// R^2 + k_b, -k_b], [-k_b, k_b]], k_b = EA / 295.2 + EA / (600 +
    /// 294.8) at the new pose (the idler's spans in series); each joint
    /// still at the old angle, 0.2 mm behind: x(t) = sum phi_i (phi_i^T M
    /// x0) cos(omega_i t). The carriage's and the rotor's following errors
    /// match it, and its modes there are the two roots.
    #[test]
    fn a_stepped_belt_axis_follows_its_two_mass_response() {
        let r = 20.0 / std::f64::consts::PI;
        let src = r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.04s)
            (def R (/ 20 pi))
            (def EA 16903N)
            (body cart (shape (box 40 30 20)) (attr density 1.27g/cm3))
            (body rotor (shape (cylinder-z 6 10)) (attr density 7.85g/cm3) (attr spin-inertia 5478.7g*mm2))
            (body idler (shape (cylinder-z 6 10)) (attr density 2.7g/cm3) (attr spin-inertia 0))
            (assembly (axis (s 300))
              (instance cart cart (at) (slide x s))
              (instance rotor rotor (at) (turn z (rad->deg (/ s R))) (effort-limits -450N*mm 450N*mm) (stiffness 22500N*mm))
              (instance idler idler (at (move 600 0 0)) (turn z (rad->deg (/ s R))))
              (relate rotor idler (belt-span (vec (vec 0 R 0) (vec 600 R 0)) 50N 0.95 (tensile-stiffness EA)))
              (relate idler cart (belt-span (vec (vec 600 (neg R) 0) (vec (+ s 5) (neg R) 0)) 50N 0.95 (tensile-stiffness EA)))
              (relate cart rotor (belt-span (vec (vec (- s 5) (neg R) 0) (vec 0 (neg R) 0)) 50N 0.95 (tensile-stiffness EA))))
            (assembly root (instance a (axis (select (- (/ t 1s) 0.002) 300 300.2))))
            (springs) (def cart-lag (following-error root/a/cart))
            (def rotor-lag (following-error root/a/rotor))"#;
        let d = doc(src);
        let s = sim(&d);
        let (ci, _) = s.released("root/a/cart").unwrap();
        let (ri, turn) = s.released("root/a/rotor").unwrap();
        assert!(turn && s.released("root/a/idler").unwrap_err().contains("carries no inertia"));
        assert_eq!(s.model().released().len(), 3);
        let mc = crate::model::modes::modes(&d.reloader().unwrap(), "root", 0.0).unwrap().model.bodies.iter().find(|b| b.label == "root/a/cart").unwrap().mass;
        let ea = 16903.0;
        let kb = ea / 295.2 + ea / (600.0 + 294.8);
        let (m1, m2, k11, k12, k22) = (5478.7 / (r * r), mc, 22500.0 / (r * r) + kb, -kb, kb);
        // The two modes, omega^2 (in 1e6 / s^2: N/mm over g), and shapes.
        let (a, b, c) = (m1 * m2, -(k11 * m2 + k22 * m1), k11 * k22 - k12 * k12);
        let sq = (b * b - 4.0 * a * c).sqrt();
        let w2 = [(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)];
        let phi: Vec<[f64; 2]> = w2.iter().map(|&w| [-k12, k11 - w * m1]).collect();
        let x0 = [-0.2, -0.2];
        let amp: Vec<f64> = phi.iter().map(|p| (p[0] * m1 * x0[0] + p[1] * m2 * x0[1]) / (p[0] * m1 * p[0] + p[1] * m2 * p[1])).collect();
        let hand = |t: f64| -> [f64; 2] {
            let mut x = [0.0; 2];
            for i in 0..2 {
                for k in 0..2 {
                    x[k] += phi[i][k] * amp[i] * ((1e6 * w2[i]).sqrt() * t).cos();
                }
            }
            x
        };
        let modes = crate::model::modes::modes(&d.reloader().unwrap(), "root", 0.01).unwrap();
        let f: Vec<f64> = modes.modes.iter().map(|x| x.frequency).collect();
        assert!(
            f.len() == 2 && (f[0] / hz(1e6 * w2[0]) - 1.0).abs() < 1e-9 && (f[1] / hz(1e6 * w2[1]) - 1.0).abs() < 1e-9,
            "{f:?} against {} {}",
            hz(1e6 * w2[0]),
            hz(1e6 * w2[1])
        );
        let mut worst = [0.0f64; 2];
        for k in 0..380 {
            let t = 0.002 + 1e-9 + k as f64 * 1e-4;
            let want = hand(t - 0.002);
            let (xc, xr) = (s.at(ci, t).unwrap()[0], s.at(ri, t).unwrap()[0] * r);
            worst = [worst[0].max((xc - want[1]).abs()), worst[1].max((xr - want[0]).abs())];
        }
        assert!(worst[0] < 2e-5 && worst[1] < 2e-5, "{worst:?} mm from the two-mass response");
        // Its ringing, by mode (before the fix: the turning points of two
        // modes ringing together gave a frequency of neither): each joint's
        // the mode carrying most of it, with its share, as the closed form
        // splits it.
        let r = s.springs();
        for (k, path) in [(1, "root/a/cart"), (0, "root/a/rotor")] {
            let parts: Vec<f64> = (0..2).map(|i| (amp[i] * phi[i][k]).abs()).collect();
            let top = if parts[0] > parts[1] { 0 } else { 1 };
            let ring = r.joints.iter().find(|f| f.path == path).unwrap().ringing.as_ref().unwrap();
            let (f, share) = (ring.frequency.unwrap(), ring.share.unwrap());
            assert!(
                (f / hz(1e6 * w2[top]) - 1.0).abs() < 1e-6 && (share - parts[top] / (parts[0] + parts[1])).abs() < 2e-3,
                "{path}: {f} Hz, {share} against {} Hz, {}",
                hz(1e6 * w2[top]),
                parts[top] / (parts[0] + parts[1])
            );
        }
        // Still before the step.
        assert!(s.at(ci, 0.001).unwrap()[0].abs() < 1e-12);
    }

    /// The same axis with its rotor (no stiffness) commanded 10 N mm and its
    /// coordinate simulated: the rotor keeps the coordinate (its command
    /// moves it), the carriage is released on its spans. By hand, from rest,
    /// the drive's force F = 10 / R on the rotor's rim and the carriage
    /// (J / R^2 + m) accelerate at a = F / (J / R^2 + m), and the carriage
    /// lags by delta = -(m a / k_b) (1 - cos omega t), omega^2 = k_b (1 / m +
    /// R^2 / J): 3.8 um at 278 Hz, matched over its first 5 ms (as it moves
    /// 0.12 mm, its spans' stiffness changes by 4e-4).
    #[test]
    fn a_commanded_motor_drags_its_sprung_carriage_behind() {
        let r = 20.0 / std::f64::consts::PI;
        let d = doc(r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.005s)
            (def R (/ 20 pi))
            (def EA 16903N)
            (body cart (shape (box 40 30 20)) (attr density 1.27g/cm3))
            (body rotor (shape (cylinder-z 6 10)) (attr density 7.85g/cm3) (attr spin-inertia 5478.7g*mm2))
            (body idler (shape (cylinder-z 6 10)) (attr density 2.7g/cm3) (attr spin-inertia 0))
            (assembly (axis (s 300))
              (instance cart cart (at) (slide x s))
              (instance rotor rotor (at) (turn z (rad->deg (/ s R))) (effort-limits -450N*mm 450N*mm))
              (instance idler idler (at (move 600 0 0)) (turn z (rad->deg (/ s R))))
              (relate rotor idler (belt-span (vec (vec 0 R 0) (vec 600 R 0)) 50N 0.95 (tensile-stiffness EA)))
              (relate idler cart (belt-span (vec (vec 600 (neg R) 0) (vec (+ s 5) (neg R) 0)) 50N 0.95 (tensile-stiffness EA)))
              (relate cart rotor (belt-span (vec (vec (- s 5) (neg R) 0) (vec 0 (neg R) 0)) 50N 0.95 (tensile-stiffness EA))))
            (assembly root (instance a (axis (simulate 300))))
            (command root/a/rotor 10)
            (springs) (def lag (following-error root/a/cart))"#);
        let s = sim(&d);
        assert_eq!(s.model().released().iter().map(|r| r.0).collect::<Vec<_>>(), ["root/a/idler", "root/a/cart"]);
        let (ci, _) = s.released("root/a/cart").unwrap();
        let mc = s.model().bodies.iter().find(|b| b.label == "root/a/cart").unwrap().mass;
        let kb = 16903.0 * (1.0 / 295.0 + 1.0 / 895.0);
        let (j, f) = (5478.7, 10.0 / r);
        let a = 1e6 * f / (j / (r * r) + mc);
        let w = (1e6 * kb * (1.0 / mc + r * r / j)).sqrt();
        let mut worst = 0.0f64;
        for k in 1..=50 {
            let t = k as f64 * 1e-4;
            let want = -(1e-6 * mc * a / kb) * (1.0 - (w * t).cos());
            worst = worst.max((s.at(ci, t).unwrap()[0] - want).abs());
        }
        assert!(
            worst < 2e-6 && (1e-6 * mc * a / kb - 3.8e-3).abs() < 1e-4 && (w / (2.0 * std::f64::consts::PI) - 278.0).abs() < 1.0,
            "{worst} mm off; lag {}, {} Hz",
            1e-6 * mc * a / kb,
            w / (2.0 * std::f64::consts::PI)
        );
        // Its coordinate is the rotor's rim: half a t^2, less the carriage's
        // share of the lag.
        let x = s.at(0, 0.005).unwrap()[0];
        let want = 300.0 + 0.5 * a * 0.005 * 0.005 - mc / (mc + j / (r * r)) * (-(1e-6 * mc * a / kb) * (1.0 - (w * 0.005).cos()));
        assert!((x - want).abs() < 2e-6, "{x} against {want}");
    }

    /// A rigid belt between joints the springs release (modes.rs's rotor
    /// geared rigidly to a shaft that a stiff belt joins to a flywheel) is
    /// refused in a simulation, naming it (the modes hold it); a following error of a joint no spring
    /// releases names those that are.
    #[test]
    fn a_simulation_refuses_a_rigid_line_between_released_joints_and_names_its_joints() {
        let d = doc(r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.1s)
            (body rotor (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 1000g*mm2))
            (body shaft (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 3000g*mm2))
            (body fly (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 8000g*mm2))
            (assembly (m (s 0))
              (instance rotor rotor (at) (turn z s) (stiffness 20000N*mm))
              (instance shaft shaft (at (move 100 0 0)) (turn z (/ s 2)))
              (instance fly fly (at (move 100 200 0)) (turn z (/ s 2)))
              (relate rotor shaft (belt-drive 10 20 50N 0.95))
              (relate shaft fly (belt-drive 10 10 50N 0.95 (tensile-stiffness 10000N))))
            (assembly root (instance a (m 0)))
            (springs) (def lag (following-error root/a/fly))"#);
        let e = crate::simulation(&d).unwrap().err().expect("refused");
        assert!(e.contains("root/a/rotor~root/a/shaft (belt-drive) ties") && e.contains("does not yet hold a rigid line"), "{e}");
        let d = doc(r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.1s)
            (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly (axis (c 0) (e 0)) (instance b b (at) (slide x c) (stiffness 100N/mm)) (instance f b (at (move 0 300 0)) (slide x e)))
            (assembly root (instance a (axis (* 10 (/ t 1s)) (* 5 (/ t 1s)))))
            (springs) (def lag (following-error root/a/f))"#);
        let e = sim(&d).released("root/a/f").unwrap_err();
        assert!(e.contains("root/a/f is no joint released for springs (those are: root/a/b)"), "{e}");
    }

    /// A simulation refuses a negative or infinite tensile stiffness, as a
    /// mode does (it stepped a spring pushing its stretch further): a
    /// param's where it starts, one known at load as it is read.
    #[test]
    fn a_simulation_refuses_a_stiffness_below_0_or_not_finite() {
        let src = |ea: &str, param: &str| {
            format!(
                r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.1s) {param}
                (body shaft (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 3000g*mm2))
                (body fly (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 8000g*mm2))
                (assembly (m (s 0))
                  (instance shaft shaft (at) (turn z s) (stiffness 20000N*mm))
                  (instance fly fly (at (move 0 200 0)) (turn z s))
                  (relate shaft fly (belt-drive 10 10 50N 0.95 (tensile-stiffness {ea}))))
                (assembly root (instance a (m 0)))
                (springs) (def lag (following-error root/a/fly))"#
            )
        };
        assert!(crate::simulation(&doc(&src("10000N", ""))).unwrap().is_ok());
        let e = Doc::load_str_in(&REGISTRY, &src("-1000N", ""), "springs.ink", Path::new("."), &[]).err().expect("loaded: -1000N");
        assert!(e.message.contains("a stiffness is at least 0 and finite"), "{}", e.message);
        for (ea, param) in [("ea", "(param ea -1000N)"), ("(* ea (exp 1000))", "(param ea 1N)")] {
            let e = crate::simulation(&doc(&src(ea, param))).unwrap().err().unwrap_or_else(|| panic!("{ea} {param} simulated"));
            assert!(e.contains("a stiffness is at least 0 and finite"), "{ea} {param}: {e}");
        }
    }

    /// A 1 kg slide along gravity on (stiffness 100N/mm), its command held:
    /// it starts where its spring holds it, m g / k = 0.0981 mm down, and
    /// stays there (from its command's place it would ring by as much).
    #[test]
    fn a_spring_starts_where_it_holds_its_load() {
        let d = doc(r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.05s)
            (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly (axis (c 0)) (instance b b (at) (slide z c) (stiffness 100N/mm)))
            (assembly root (instance a (axis 0)))
            (springs) (def lag (following-error root/a/b))"#);
        let s = sim(&d);
        let (i, _) = s.released("root/a/b").unwrap();
        let sag = -1e-6 * s.model().bodies.iter().find(|b| b.label == "root/a/b").unwrap().mass * 9810.0 / 100.0;
        for t in [0.0, 0.013, 0.05] {
            let x = s.at(i, t).unwrap();
            assert!((x[0] - sag).abs() < 1e-9 && x[1].abs() < 1e-6, "at {t}: {x:?} against {sag}");
        }
    }

    /// A body a written motion moves beside a joint its springs release,
    /// tied to it by their drives alone (a second slide on the same
    /// command), is refused in a simulation, named, as the modes refuse it.
    #[test]
    fn a_body_tied_by_drives_alone_to_a_released_joint_is_refused() {
        let d = doc(r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.1s)
            (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly (axis (c 0)) (instance b b (at) (slide x c) (stiffness 100N/mm)) (instance f b (at (move 0 300 0)) (slide x c)))
            (assembly root (instance a (axis (* 10 (/ t 1s)))))
            (springs) (def lag (following-error root/a/b))"#);
        let e = crate::simulation(&d).unwrap().err().expect("refused");
        assert!(e.contains("root/a/f (with a.c) move with a coordinate that also commands joints released for springs"), "{e}");
    }

    /// A spring line pulls with its pretension plus its stiffness times its
    /// stretch, never below 0: slack, it holds nothing and its pretension's
    /// pull on the joint goes with it.
    #[test]
    fn a_slack_line_pulls_nothing() {
        let at = super::At { g: vec![vec![1.0, -1.0]], kl: vec![10.0], t0: vec![5.0], kj: vec![0.0, 2.0], drive: vec![[0.0; 3]; 2], heft: vec![1.0; 2], mass: vec![1.0; 2] };
        // Taut: 5 + 10 (0.2 - 0) = 7 N, pulling joint 0 back by its 2 N of stretch.
        assert_eq!(at.tensions(&[0.2, 0.0]), [7.0]);
        assert_eq!(at.force(&[0.2, 0.0]), [-2.0, 2.0]);
        assert_eq!(at.stiffness(&[0.2, 0.0]), [[10.0, -10.0], [-10.0, 12.0]]);
        // Slack past -0.5: none, and its 5 N of pretension no longer balanced.
        assert_eq!(at.tensions(&[-1.0, 0.0]), [0.0]);
        assert_eq!(at.force(&[-1.0, 0.0]), [5.0, -5.0]);
        assert_eq!(at.stiffness(&[-1.0, 0.0]), [[0.0, 0.0], [0.0, 2.0]]);
    }

    /// A rotor of 1e6 g mm2 on 1e6 N mm/rad stepped one full step, 1.8
    /// deg, at 0.01 s rings as -1.8 deg cos(1000 (t - 0.01)) whatever its
    /// command's origin: after 5556 turns (2e6 deg) alike (before the fix
    /// the step went missing there, its jump judged against the angle).
    #[test]
    fn a_step_rings_alike_whatever_its_commands_origin() {
        let run = |from: f64| {
            let d = doc(&format!(
                r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.02s)
                (body r (shape (cylinder-z 10 10)) (attr density 1g/cm3) (attr spin-inertia 1000000g*mm2))
                (assembly (m (a 0)) (instance r r (at) (turn z a) (stiffness 1000000N*mm)))
                (assembly root (instance m (m (select (- (/ t 1s) 0.01) {from} (+ {from} 1.8)))))
                (springs) (def lag (following-error root/m/r))"#
            ));
            let s = sim(&d);
            let (i, _) = s.released("root/m/r").unwrap();
            [0.012, 0.015, 0.02].map(|t| s.at(i, t).unwrap()[0])
        };
        let (near, far) = (run(0.0), run(2e6));
        for (k, t) in [0.012f64, 0.015, 0.02].into_iter().enumerate() {
            let want = -1.8f64.to_radians() * (1000.0f64 * (t - 0.01)).cos();
            assert!((near[k] - want).abs() < 1e-6 && (far[k] - want).abs() < 1e-6, "at {t}: {} and {} against {want}", near[k], far[k]);
        }
    }

    /// A 1000 mm bar of 1000 g (333333333 g mm2 about its pivot) on as many
    /// N mm/rad, stepped 1e-5 rad at 0.01 s, rings as 1e-5 cos(1000 (t -
    /// 0.01)) rad, 10 um at its tip, to the span's end at 1 s (before the
    /// fix, a microradian's tolerance let the steps smooth it away, -3e-15 rad
    /// at 1 s): its angle is held to a nanometre at the bar's radius of
    /// gyration: over its last 50 ms it swings by 1e-5 rad within 5% (9.78e-6
    /// after 160 periods at its tolerance) at 1000 rad/s within 1e-3.
    #[test]
    fn a_long_arms_small_ringing_is_kept() {
        let d = doc(r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 1s)
            (body bar (shape (move (box 1000 10 10) 500 0 0)) (attr mass 1000g) (attr spin-inertia 333333333.3333333g*mm2))
            (assembly (m (a 0)) (instance b bar (at) (turn z a) (stiffness 333333333.3333333N*mm)))
            (assembly root (instance m (m (select (- (/ t 1s) 0.01) (rad->deg 1e-5) 0))))
            (springs) (def lag (following-error root/m/b))"#);
        let s = sim(&d);
        let (i, _) = s.released("root/m/b").unwrap();
        // Over its last 50 ms: its swing, and its frequency from its
        // turning points (the phase drifting over 160 periods).
        let turns = s.moments(i, (0.95, 1.0));
        let swing = turns.iter().map(|&t| s.at(i, t).unwrap()[0].abs()).fold(0.0, f64::max);
        let f = (turns.len() - 1) as f64 / (2.0 * (turns[turns.len() - 1] - turns[0]));
        assert!((swing / 1e-5 - 1.0).abs() < 0.05 && (f * 2.0 * std::f64::consts::PI / 1000.0 - 1.0).abs() < 1e-3, "{swing} rad at {f} Hz");
    }

    /// A 1 kg slide along gravity hung on a line of 100 N/mm with 5 N of
    /// pretension, pulling it the way gravity does: slack, the line takes
    /// 5 N off at most, and 9.81 N is needed, so no rest exists and the
    /// simulation says so (before the fix it started from the taut line's
    /// sag, 0.0981 mm, accelerating).
    #[test]
    fn a_start_with_no_rest_is_refused() {
        let d = doc(r#"(import "std") (gravity 0 0 9.81m/s2) (t-span 0 0.1s)
            (body anchor (shape (box 1 1 1)))
            (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (fn elastic (fa fb sa sb) (relation elastic (couple (vec (vec (vec 0 0 0) (vec 0 0 1) 0 100 100)) 5N)))
            (assembly root (instance anchor anchor) (instance b b (at) (slide z 0)) (relate anchor b elastic))
            (springs) (def lag (following-error root/b))"#);
        let e = crate::simulation(&d).unwrap().err().expect("refused");
        assert!(e.contains("the springs do not hold root/b at rest") && e.contains("9.81"), "{e}");
    }

    /// Two 1 kg bodies stacked on z slides, on 100 and 1000 N/mm, their
    /// command stepping 0.1 mm: the top's ringing is split alike with
    /// gravity along the slides (the cart sagging 0.196 mm) as without,
    /// about the springs' static rest (before the fix the sag, not removed,
    /// weighed on the modes: 35.14 Hz for 227.94). By hand the split is even
    /// (M = [[2, 1], [1, 1]] kg, K = diag(100, 1000) N/mm, x0 = (-0.1 mm, 0):
    /// 4.9938e-3 mm in each mode, far above the band; at a 1 um step, 50 nm
    /// in each, half the band, the 228 Hz part is smoothed by 5%).
    #[test]
    fn a_ringings_mode_ignores_its_sag() {
        let mode = |g: &str| {
            let d = doc(&format!(
                r#"(import "std") (gravity 0 0 {g}) (t-span 0 0.1s)
                (body b (shape (box 100 100 100)) (attr density 1g/cm3))
                (assembly (axis (c 0))
                  (instance cart b (at) (slide z c) (stiffness 100N/mm))
                  (instance top b (at (move 0 0 200) (frame-of cart)) (slide z 0) (stiffness 1000N/mm)))
                (assembly root (instance a (axis (select (- (/ t 1s) 0.01) 0 0.1))))
                (springs) (def lag (following-error root/a/cart))"#
            ));
            let r = sim(&d).springs();
            let top = r.joints.iter().find(|f| f.path == "root/a/top").unwrap().ringing.as_ref().unwrap();
            (top.frequency.unwrap(), top.share.unwrap())
        };
        let (still, sagging) = (mode("0"), mode("-9.81m/s2"));
        assert!((still.0 / sagging.0 - 1.0).abs() < 1e-9 && (still.1 - sagging.1).abs() < 1e-4 && (still.1 - 0.5).abs() < 1e-3, "{still:?} {sagging:?}");
    }

    /// 1 kg on a slide its (stiffness K) holds, its command stepping from
    /// A to 0 mm at 0.01 s: it rings as A cos omega (t - 0.01), omega^2 = K
    /// / m.
    fn ringing(k: f64, a: f64, span: f64) -> Doc {
        doc(&format!(
            r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 {span}s)
            (body b (shape (box 100 100 100)) (attr mass 1000g))
            (assembly (axis (c 0)) (instance b b (at) (slide x c) (stiffness {k}N/mm)))
            (assembly root (instance a (axis (select (- (/ t 1s) 0.01) {a} 0))))
            (springs) (def lag (following-error root/a/b))"#
        ))
    }

    /// Radau IIA damps an undamped ringing a little each period, near the
    /// tolerance about 0.2 tolerances a period (0.15 um at
    /// 335 Hz, 1 nm a tolerance, gone by 2.4 s, 150 tolerances from the
    /// closed form, past the band's 100 then): the band widens by 0.3 a
    /// period of the springs' fastest mode stepped, so the closed form
    /// stays within it over 1000 periods at 335 Hz and 750 at 50 Hz.
    #[test]
    fn a_long_ringing_stays_within_its_widening_band() {
        let a = 1.5e-4;
        for (k, span) in [(4430.0, 3.0), (100.0, 15.0)] {
            let s = sim(&ringing(k, a, span));
            let (i, _) = s.released("root/a/b").unwrap();
            let w = (1e6 * k / 1000.0).sqrt();
            let (mut worst, mut lost) = (0.0f64, 0.0f64);
            let n = 40000;
            for m in 0..=n {
                let t = 0.0101 + (span - 0.0102) * m as f64 / n as f64;
                let e = (s.at(i, t).unwrap()[0] - a * (w * (t - 0.01)).cos()).abs();
                (worst, lost) = (worst.max(e / s.band(i, t)), lost.max(e));
            }
            assert!(worst < 1.0 && lost > 1.2e-4, "{k} N/mm over {span} s: {lost} mm from the closed form, {worst} of its band");
        }
    }

    /// A simulation of springs is stepped by Radau IIA throughout:
    /// Dormand-Prince's stability function passes 1 on the
    /// imaginary axis, so its steps grow an undamped ringing near the
    /// tolerance (0.1 um at 335 Hz swung 2.27 times as wide by 0.3 s).
    /// Radau's damp it, by less than its band.
    #[test]
    fn a_springs_ringing_is_stepped_by_radau_throughout() {
        let s = sim(&ringing(4430.0, 1e-4, 0.3));
        let (i, _) = s.released("root/a/b").unwrap();
        let swing = (0..=2000).map(|m| s.at(i, 0.25 + 0.05 * m as f64 / 2000.0).unwrap()[0].abs()).fold(0.0, f64::max);
        assert!(swing < 1e-4 && 1e-4 - swing < s.band(i, 0.25), "{swing} mm against 1e-4, its band {}", s.band(i, 0.25));
    }

    /// The stepped belt axis with its rotor a plain written turn (no effort
    /// limits, no stiffness): released (the spans' ends rest on it), it and
    /// the cart make a motion no spring resists, so their following error
    /// would drift with the whole command (before the fix the cart read
    /// -1.24 mm, the command's own travel): refused, naming them.
    #[test]
    fn a_motion_nothing_holds_is_refused_naming_its_joints() {
        let d = doc(r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.008s)
            (def R (/ 20 pi))
            (def EA 16903N)
            (body cart (shape (box 40 30 20)) (attr density 1.27g/cm3))
            (body rotor (shape (cylinder-z 6 10)) (attr density 7.85g/cm3) (attr spin-inertia 5478.7g*mm2))
            (body idler (shape (cylinder-z 6 10)) (attr density 2.7g/cm3) (attr spin-inertia 0))
            (assembly (axis (s 300))
              (instance cart cart (at) (slide x s))
              (instance rotor rotor (at) (turn z (rad->deg (/ s R))))
              (instance idler idler (at (move 600 0 0)) (turn z (rad->deg (/ s R))))
              (relate rotor idler (belt-span (vec (vec 0 R 0) (vec 600 R 0)) 50N 0.95 (tensile-stiffness EA)))
              (relate idler cart (belt-span (vec (vec 600 (neg R) 0) (vec (+ s 5) (neg R) 0)) 50N 0.95 (tensile-stiffness EA)))
              (relate cart rotor (belt-span (vec (vec (- s 5) (neg R) 0) (vec 0 (neg R) 0)) 50N 0.95 (tensile-stiffness EA))))
            (assembly root (instance a (axis (+ 300 (* 10 (- 1 (cos (* 20 180deg/s t))))))))
            (springs) (def lag (following-error root/a/cart))"#);
        let e = crate::simulation(&d).unwrap().err().expect("refused");
        assert!(e.contains("nothing holds root/a/rotor, root/a/cart along a motion they make together"), "{e}");
    }

    /// Two 1 kg sliders on their own springs, 1 and 1e13 N/mm: the soft one
    /// is held by its own (its mode sqrt(1e6 x 1 / 1000) / 2 pi, 5.03 Hz),
    /// however stiff the other, and the simulation starts (it was refused as
    /// held by nothing, its mode under 1e-12 of the stiff one's).
    #[test]
    fn a_soft_spring_beside_a_stiff_one_holds_its_joint() {
        let d = doc(r#"(import "std") (gravity 0 0 -9.81m/s2) (t-span 0 0.001s)
            (body b (shape (box 50 50 50)) (attr mass 1000g))
            (assembly root (instance soft b (at) (slide x 0) (stiffness 1N/mm))
              (instance stiff b (at (move 0 100 0)) (slide x 0) (stiffness 1e13N/mm)))
            (springs) (def lag (following-error root/soft))"#);
        let s = sim(&d);
        assert!(s.released("root/soft").is_ok() && s.released("root/stiff").is_ok());
    }
}
